// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Stage-3 extraction: emit candidate names, or re-scan the BAM and emit
//! candidate read sequences as FASTA.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::num::NonZeroUsize;
use std::path::Path;

use anyhow::{Context, Result};
use log::{info, warn};
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use rustc_hash::FxHashSet;

use crate::logging::elapsed;

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
    /// Number of candidate names that did not have a usable BAM sequence.
    pub missing: u64,
    /// Number of BAM records scanned while extracting FASTA.
    pub scanned: u64,
}

/// Writes the candidate names — sorted, deduplicated, one per line — directly
/// to `output`.
pub fn write_names(output: &Path, candidates: &CandidateSet) -> Result<ExtractStats> {
    let mut sorted: Vec<&[u8]> = candidates.iter().map(Vec::as_slice).collect();
    sorted.sort_unstable();

    let file =
        File::create(output).with_context(|| format!("failed to create {}", output.display()))?;
    let mut writer = BufWriter::with_capacity(OUTPUT_BUFFER_BYTES, file);

    let mut written = 0_u64;
    for name in sorted {
        writer.write_all(name)?;
        writer.write_all(b"\n")?;
        written += 1;
    }
    writer.flush()?;

    info!(
        "[{}] wrote {written} names to {}",
        elapsed(),
        output.display()
    );
    Ok(ExtractStats {
        written,
        missing: 0,
        scanned: 0,
    })
}

/// Re-scans `bam_path` and writes FASTA records for candidate QNAMEs.
///
/// Each candidate is written at most once. Secondary records are ignored;
/// supplementary records are allowed because they normally carry the same
/// full read sequence and can appear before the primary record. Records with
/// an empty BAM SEQ field are skipped, allowing a later record with the same
/// QNAME to provide the sequence.
pub fn write_fasta_from_bam(
    bam_path: &Path,
    output: &Path,
    candidates: &CandidateSet,
    bgzf_workers: NonZeroUsize,
) -> Result<ExtractStats> {
    let file = File::open(bam_path)
        .with_context(|| format!("failed to open BAM {}", bam_path.display()))?;
    let bgzf_reader = bgzf::MultithreadedReader::with_worker_count(bgzf_workers, file);
    let mut reader = bam::io::Reader::from(bgzf_reader);
    let _header = reader
        .read_header()
        .with_context(|| format!("failed to read BAM header from {}", bam_path.display()))?;

    let file =
        File::create(output).with_context(|| format!("failed to create {}", output.display()))?;
    let mut writer = BufWriter::with_capacity(OUTPUT_BUFFER_BYTES, file);

    info!(
        "[{}] fasta-extract begin bam={} output={} bgzf_workers={}",
        elapsed(),
        bam_path.display(),
        output.display(),
        bgzf_workers.get()
    );

    let mut seen: CandidateSet = FxHashSet::default();
    let mut record = bam::Record::default();
    let mut scanned = 0_u64;
    let mut written = 0_u64;

    loop {
        let bytes = reader.read_record(&mut record).with_context(|| {
            format!("failed to read BAM record #{scanned} during FASTA extraction")
        })?;
        if bytes == 0 {
            break;
        }
        scanned += 1;

        let flags = record.flags();
        if flags.is_secondary() {
            continue;
        }

        let Some(name) = record.name() else {
            continue;
        };
        let qname: &[u8] = name.as_ref();

        if !candidates.contains(qname) || seen.contains(qname) {
            continue;
        }

        let sequence = record.sequence();
        if sequence.is_empty() {
            continue;
        }

        let seq = original_read_sequence(&record);
        write_fasta(&mut writer, qname, &seq)?;
        seen.insert(qname.to_vec());
        written += 1;

        if scanned % 1_000_000 == 0 {
            info!(
                "[{}] fasta-extract progress scanned={scanned} written={written}",
                elapsed()
            );
        }
    }

    writer.flush()?;

    let missing = candidates.len() as u64 - seen.len() as u64;
    if missing > 0 {
        warn!(
            "[{}] {missing} candidate name(s) had no usable sequence in {}",
            elapsed(),
            bam_path.display()
        );
    }

    info!(
        "[{}] fasta-extract done scanned={scanned} written={written} missing={missing}",
        elapsed()
    );

    Ok(ExtractStats {
        written,
        missing,
        scanned,
    })
}

fn original_read_sequence(record: &bam::Record) -> Vec<u8> {
    let bases = record.sequence().iter();
    if record.flags().is_reverse_complemented() {
        bases.rev().map(complement_base).collect()
    } else {
        bases.collect()
    }
}

fn complement_base(b: u8) -> u8 {
    match b {
        b'=' => b'=',
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'U' => b'A',
        b'W' => b'W',
        b'S' => b'S',
        b'M' => b'K',
        b'K' => b'M',
        b'R' => b'Y',
        b'Y' => b'R',
        b'B' => b'V',
        b'D' => b'H',
        b'H' => b'D',
        b'V' => b'B',
        b'N' => b'N',
        _ => b'N',
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    use noodles_sam::{
        self as sam,
        alignment::{
            io::Write as _,
            record::Flags,
            record_buf::{QualityScores, Sequence},
            RecordBuf,
        },
    };
    use tempfile::tempdir;

    fn read_file(path: &Path) -> Vec<u8> {
        let mut buf = Vec::new();
        std::fs::File::open(path)
            .unwrap()
            .read_to_end(&mut buf)
            .unwrap();
        buf
    }

    fn candset(items: &[&[u8]]) -> CandidateSet {
        items.iter().map(|s| s.to_vec()).collect()
    }

    fn record(name: &[u8], seq: &[u8], flags: Flags) -> RecordBuf {
        RecordBuf::builder()
            .set_name(name.to_vec())
            .set_flags(flags)
            .set_sequence(Sequence::from(seq))
            .set_quality_scores(QualityScores::from(vec![30; seq.len()]))
            .build()
    }

    fn write_bam(dir: &Path, name: &str, records: &[RecordBuf]) -> std::path::PathBuf {
        let path = dir.join(name);
        let file = std::fs::File::create(&path).unwrap();
        let mut writer = bam::io::Writer::new(file);
        let header = sam::Header::default();
        writer.write_header(&header).unwrap();
        for record in records {
            writer.write_alignment_record(&header, record).unwrap();
        }
        writer.try_finish().unwrap();
        path
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
    fn fasta_from_bam_filters_candidates() {
        let dir = tempdir().unwrap();
        let bam = write_bam(
            dir.path(),
            "reads.bam",
            &[
                record(b"r1", b"ACGT", Flags::UNMAPPED),
                record(b"r2", b"TTTT", Flags::UNMAPPED),
                record(b"r3", b"GGGG", Flags::UNMAPPED),
            ],
        );
        let out = dir.path().join("out.fa");
        let cands = candset(&[b"r1", b"r3"]);

        let stats =
            write_fasta_from_bam(&bam, &out, &cands, NonZeroUsize::new(1).unwrap()).unwrap();

        assert_eq!(stats.written, 2);
        assert_eq!(stats.missing, 0);
        assert_eq!(stats.scanned, 3);
        assert_eq!(read_file(&out), b">r1\nACGT\n>r3\nGGGG\n");
    }

    #[test]
    fn supplementary_duplicate_writes_once() {
        let dir = tempdir().unwrap();
        let bam = write_bam(
            dir.path(),
            "reads.bam",
            &[
                record(b"r1", b"ACGT", Flags::SUPPLEMENTARY),
                record(b"r1", b"GGGG", Flags::UNMAPPED),
            ],
        );
        let out = dir.path().join("out.fa");
        let cands = candset(&[b"r1"]);

        let stats =
            write_fasta_from_bam(&bam, &out, &cands, NonZeroUsize::new(1).unwrap()).unwrap();

        assert_eq!(stats.written, 1);
        assert_eq!(stats.missing, 0);
        assert_eq!(read_file(&out), b">r1\nACGT\n");
    }

    #[test]
    fn secondary_records_are_ignored() {
        let dir = tempdir().unwrap();
        let bam = write_bam(
            dir.path(),
            "reads.bam",
            &[record(b"r1", b"ACGT", Flags::SECONDARY)],
        );
        let out = dir.path().join("out.fa");
        let cands = candset(&[b"r1"]);

        let stats =
            write_fasta_from_bam(&bam, &out, &cands, NonZeroUsize::new(1).unwrap()).unwrap();

        assert_eq!(stats.written, 0);
        assert_eq!(stats.missing, 1);
        assert!(read_file(&out).is_empty());
    }

    #[test]
    fn reverse_complemented_records_emit_original_orientation() {
        let dir = tempdir().unwrap();
        let bam = write_bam(
            dir.path(),
            "reads.bam",
            &[record(
                b"r1",
                b"ACGA",
                Flags::UNMAPPED | Flags::REVERSE_COMPLEMENTED,
            )],
        );
        let out = dir.path().join("out.fa");
        let cands = candset(&[b"r1"]);

        write_fasta_from_bam(&bam, &out, &cands, NonZeroUsize::new(1).unwrap()).unwrap();

        assert_eq!(read_file(&out), b">r1\nTCGT\n");
    }

    #[test]
    fn empty_sequence_is_missing_unless_later_record_has_sequence() {
        let dir = tempdir().unwrap();
        let bam = write_bam(
            dir.path(),
            "reads.bam",
            &[
                record(b"r1", b"", Flags::UNMAPPED),
                record(b"r1", b"ACGT", Flags::UNMAPPED),
            ],
        );
        let out = dir.path().join("out.fa");
        let cands = candset(&[b"r1"]);

        let stats =
            write_fasta_from_bam(&bam, &out, &cands, NonZeroUsize::new(1).unwrap()).unwrap();

        assert_eq!(stats.written, 1);
        assert_eq!(stats.missing, 0);
        assert_eq!(read_file(&out), b">r1\nACGT\n");
    }

    #[test]
    fn empty_sequence_without_later_sequence_counts_missing() {
        let dir = tempdir().unwrap();
        let bam = write_bam(
            dir.path(),
            "reads.bam",
            &[record(b"r1", b"", Flags::UNMAPPED)],
        );
        let out = dir.path().join("out.fa");
        let cands = candset(&[b"r1"]);

        let stats =
            write_fasta_from_bam(&bam, &out, &cands, NonZeroUsize::new(1).unwrap()).unwrap();

        assert_eq!(stats.written, 0);
        assert_eq!(stats.missing, 1);
        assert!(read_file(&out).is_empty());
    }
}
