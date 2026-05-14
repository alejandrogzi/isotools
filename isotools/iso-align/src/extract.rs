// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Stage-3 extraction: stream the user's FASTQ/FASTA, emit FASTA, FASTQ, or
//! a sorted name list for the candidate set.
//!
//! This stage is single-threaded by design: needletail's decompression-and-
//! parse plus a hashset lookup will saturate disk on most setups, so feeding
//! a thread pool here is more complexity than it's worth.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{anyhow, Context, Result};
use log::{info, warn};
use needletail::parse_fastx_file;
use rustc_hash::FxHashSet;

use crate::logging::elapsed;
use crate::types::OutputFormat;

/// Output buffer size — large enough that small per-record writes don't churn
/// the kernel, small enough not to balloon RSS.
const OUTPUT_BUFFER_BYTES: usize = 4 * 1024 * 1024;

/// Set of candidate read names (raw bytes, not UTF-8-validated).
pub type CandidateSet = FxHashSet<Vec<u8>>;

/// Result of stage 3 — useful for the caller's end-of-run summary.
#[derive(Clone, Copy, Debug, Default)]
pub struct ExtractStats {
    /// Number of records emitted to the output.
    pub written: u64,
    /// Number of candidate names that were *not* found in the FASTQ.
    /// A non-zero value usually means the BAM and FASTQ don't correspond.
    pub missing: u64,
    /// Number of FASTQ/FASTA records scanned.
    pub scanned: u64,
}

/// Writes the candidate names — sorted, deduplicated, one per line — directly
/// to `output`. Does not consult the reads file at all; the candidate set
/// already contains the canonical names from the BAM.
pub fn write_names(output: &Path, candidates: &CandidateSet) -> Result<ExtractStats> {
    let mut sorted: Vec<&[u8]> = candidates.iter().map(Vec::as_slice).collect();
    sorted.sort_unstable();

    let file = File::create(output)
        .with_context(|| format!("failed to create {}", output.display()))?;
    let mut writer = BufWriter::with_capacity(OUTPUT_BUFFER_BYTES, file);

    let mut written = 0_u64;
    for name in sorted {
        writer.write_all(name)?;
        writer.write_all(b"\n")?;
        written += 1;
    }
    writer.flush()?;

    info!("[{}] wrote {written} names to {}", elapsed(), output.display());
    Ok(ExtractStats {
        written,
        missing: 0,
        scanned: 0,
    })
}

/// Streams `reads` and writes a FASTA or FASTQ subset containing every record
/// whose id (truncated at the first whitespace, BAM-style) is in `candidates`.
///
/// `format` must be either [`OutputFormat::Fasta`] or [`OutputFormat::Fastq`].
/// FASTQ output requires the input file to actually carry quality scores;
/// FASTA inputs rendered as FASTQ would have no qualities to write.
pub fn write_sequences(
    reads: &Path,
    output: &Path,
    candidates: &CandidateSet,
    format: OutputFormat,
) -> Result<ExtractStats> {
    if matches!(format, OutputFormat::Names) {
        return Err(anyhow!(
            "write_sequences called with OutputFormat::Names; use write_names instead"
        ));
    }

    let mut reader = parse_fastx_file(reads)
        .with_context(|| format!("failed to open reads {}", reads.display()))?;

    let file = File::create(output)
        .with_context(|| format!("failed to create {}", output.display()))?;
    let mut writer = BufWriter::with_capacity(OUTPUT_BUFFER_BYTES, file);

    let mut seen: FxHashSet<Vec<u8>> = FxHashSet::default();
    let mut written = 0_u64;
    let mut scanned = 0_u64;

    while let Some(rec) = reader.next() {
        let rec = rec.with_context(|| format!("failed to parse a record in {}", reads.display()))?;
        scanned += 1;

        let id = rec.id();
        // FASTQ/FASTA headers may include a comment after the id; BAM QNAMEs
        // never do. Truncate at the first whitespace so the lookup matches.
        let qname = id_prefix(id);

        if !candidates.contains(qname) {
            continue;
        }
        // Track only matched names so we can report unfound candidates.
        seen.insert(qname.to_vec());

        match format {
            OutputFormat::Fasta => write_fasta(&mut writer, qname, &rec.seq())?,
            OutputFormat::Fastq => {
                let qual = rec.qual().ok_or_else(|| {
                    anyhow!(
                        "FASTQ output requested but {} has no quality scores",
                        reads.display()
                    )
                })?;
                write_fastq(&mut writer, qname, &rec.seq(), qual)?;
            }
            OutputFormat::Names => unreachable!(),
        }
        written += 1;

        if scanned % 1_000_000 == 0 {
            info!(
                "[{}] extract progress scanned={scanned} written={written}",
                elapsed()
            );
        }
    }

    writer.flush()?;

    let missing = candidates.len() as u64 - seen.len() as u64;
    if missing > 0 {
        warn!(
            "[{}] {missing} candidate name(s) not found in {} — usually means BAM and reads file don't correspond",
            elapsed(),
            reads.display()
        );
    }

    info!(
        "[{}] extract done scanned={scanned} written={written} missing={missing}",
        elapsed()
    );

    Ok(ExtractStats {
        written,
        missing,
        scanned,
    })
}

/// Returns the substring of `id` up to the first ASCII whitespace, matching
/// the QNAME convention from BAM. Operates on raw bytes — no UTF-8 cost.
pub fn id_prefix(id: &[u8]) -> &[u8] {
    match id.iter().position(|b| matches!(b, b' ' | b'\t')) {
        Some(p) => &id[..p],
        None => id,
    }
}

fn write_fasta<W: Write>(writer: &mut W, id: &[u8], seq: &[u8]) -> Result<()> {
    writer.write_all(b">")?;
    writer.write_all(id)?;
    writer.write_all(b"\n")?;
    writer.write_all(seq)?;
    writer.write_all(b"\n")?;
    Ok(())
}

fn write_fastq<W: Write>(writer: &mut W, id: &[u8], seq: &[u8], qual: &[u8]) -> Result<()> {
    writer.write_all(b"@")?;
    writer.write_all(id)?;
    writer.write_all(b"\n")?;
    writer.write_all(seq)?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(qual)?;
    writer.write_all(b"\n")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Read, Write};
    use tempfile::tempdir;

    fn write_file(dir: &Path, name: &str, body: &[u8]) -> std::path::PathBuf {
        let path = dir.join(name);
        std::fs::File::create(&path)
            .unwrap()
            .write_all(body)
            .unwrap();
        path
    }

    fn read_file(path: &Path) -> Vec<u8> {
        let mut buf = Vec::new();
        std::fs::File::open(path).unwrap().read_to_end(&mut buf).unwrap();
        buf
    }

    fn candset(items: &[&[u8]]) -> CandidateSet {
        items.iter().map(|s| s.to_vec()).collect()
    }

    #[test]
    fn id_prefix_truncates_at_whitespace() {
        assert_eq!(id_prefix(b"read/1 comment"), b"read/1");
        assert_eq!(id_prefix(b"read\textra"), b"read");
        assert_eq!(id_prefix(b"plain"), b"plain");
    }

    #[test]
    fn names_output_is_sorted_and_deduped() {
        let dir = tempdir().unwrap();
        let out = dir.path().join("names.txt");
        let cands = candset(&[b"r3", b"r1", b"r2"]);
        write_names(&out, &cands).unwrap();
        let body = read_file(&out);
        assert_eq!(body, b"r1\nr2\nr3\n");
    }

    #[test]
    fn fasta_extraction_filters_and_drops_comments() {
        let dir = tempdir().unwrap();
        let reads = write_file(
            dir.path(),
            "reads.fastq",
            b"@r1 some comment\nACGT\n+\n!!!!\n@r2\nTTTT\n+\n####\n@r3\nGGGG\n+\n$$$$\n",
        );
        let out = dir.path().join("out.fa");
        let cands = candset(&[b"r1", b"r3"]);
        let stats = write_sequences(&reads, &out, &cands, OutputFormat::Fasta).unwrap();
        assert_eq!(stats.written, 2);
        assert_eq!(stats.missing, 0);
        assert_eq!(stats.scanned, 3);
        let body = read_file(&out);
        assert_eq!(body, b">r1\nACGT\n>r3\nGGGG\n");
    }

    #[test]
    fn fastq_extraction_preserves_quality() {
        let dir = tempdir().unwrap();
        let reads = write_file(
            dir.path(),
            "reads.fastq",
            b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\n####\n",
        );
        let out = dir.path().join("out.fq");
        let cands = candset(&[b"r2"]);
        let stats = write_sequences(&reads, &out, &cands, OutputFormat::Fastq).unwrap();
        assert_eq!(stats.written, 1);
        let body = read_file(&out);
        assert_eq!(body, b"@r2\nTTTT\n+\n####\n");
    }

    #[test]
    fn missing_candidates_are_counted_not_fatal() {
        let dir = tempdir().unwrap();
        let reads = write_file(dir.path(), "reads.fastq", b"@r1\nACGT\n+\nIIII\n");
        let out = dir.path().join("out.fa");
        let cands = candset(&[b"r1", b"ghost"]);
        let stats = write_sequences(&reads, &out, &cands, OutputFormat::Fasta).unwrap();
        assert_eq!(stats.written, 1);
        assert_eq!(stats.missing, 1);
    }

    #[test]
    fn fastq_output_from_fasta_input_errors() {
        let dir = tempdir().unwrap();
        let reads = write_file(dir.path(), "reads.fa", b">r1\nACGT\n");
        let out = dir.path().join("out.fq");
        let cands = candset(&[b"r1"]);
        let err = write_sequences(&reads, &out, &cands, OutputFormat::Fastq).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("no quality"), "unexpected error: {msg}");
    }

    #[test]
    fn names_format_rejected_by_write_sequences() {
        let dir = tempdir().unwrap();
        let reads = write_file(dir.path(), "reads.fa", b">r1\nACGT\n");
        let out = dir.path().join("names.txt");
        let cands = candset(&[b"r1"]);
        let err = write_sequences(&reads, &out, &cands, OutputFormat::Names).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("write_names"), "unexpected error: {msg}");
    }
}
