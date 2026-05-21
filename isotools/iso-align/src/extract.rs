// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Stage-3 extraction: emit candidate names, or re-scan the BAM and emit
//! candidate read sequences as FASTA.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use log::{info, warn};
use needletail::parse_fastx_file;
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::alignment::record::cigar::op::Kind;
use rustc_hash::{FxHashMap, FxHashSet};

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
    /// Number of candidate names that did not have a usable sequence.
    pub missing: u64,
    /// Number of input records scanned while extracting FASTA.
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
/// Each candidate is written at most once. Secondary records and records with
/// an empty BAM SEQ field are ignored. For each QNAME, all usable records are
/// considered and the best available sequence is emitted: records without hard
/// clipping outrank hard-clipped records, non-supplementary records outrank
/// supplementary records, and longer sequences outrank shorter sequences.
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

    let mut best: FxHashMap<Vec<u8>, BestBamSequence> = FxHashMap::default();
    let mut order: Vec<Vec<u8>> = Vec::new();
    let mut record = bam::Record::default();
    let mut scanned = 0_u64;

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

        if !candidates.contains(qname) {
            continue;
        }

        let seq_len = record.sequence().len();
        if seq_len == 0 {
            continue;
        }

        let score = BamSeqScore {
            no_hard_clip: !has_hard_clip(&record)
                .with_context(|| format!("invalid CIGAR at record #{scanned}"))?,
            is_primary: !flags.is_supplementary(),
            seq_len,
        };

        let should_update = match best.get(qname) {
            Some(current) => score > current.score,
            None => true,
        };

        if should_update {
            let seq = original_read_sequence(&record);
            if !best.contains_key(qname) {
                order.push(qname.to_vec());
            }
            best.insert(qname.to_vec(), BestBamSequence { seq, score });
        }

        if scanned % 1_000_000 == 0 {
            info!(
                "[{}] fasta-extract progress scanned={scanned} selected={}",
                elapsed(),
                best.len()
            );
        }
    }

    let mut written = 0_u64;
    for qname in &order {
        if let Some(selected) = best.get(qname) {
            write_fasta(&mut writer, qname, &selected.seq)?;
            written += 1;
        }
    }

    writer.flush()?;

    let missing = candidates.len() as u64 - best.len() as u64;
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

/// Reads original FASTA/FASTQ records and writes FASTA records for candidate
/// QNAMEs. This is exact when the input reads files are the source used to
/// build the pass-1 BAM.
pub fn write_fasta_from_reads(
    reads_paths: &[PathBuf],
    output: &Path,
    candidates: &CandidateSet,
) -> Result<ExtractStats> {
    let file =
        File::create(output).with_context(|| format!("failed to create {}", output.display()))?;
    let mut writer = BufWriter::with_capacity(OUTPUT_BUFFER_BYTES, file);

    info!(
        "[{}] fasta-extract begin reads={} output={}",
        elapsed(),
        display_paths(reads_paths),
        output.display()
    );

    let mut seen: CandidateSet = FxHashSet::default();
    let mut scanned = 0_u64;
    let mut written = 0_u64;

    for reads_path in reads_paths {
        let mut reader = parse_fastx_file(reads_path)
            .with_context(|| format!("failed to open FASTA/FASTQ {}", reads_path.display()))?;
        let mut file_scanned = 0_u64;

        while let Some(record) = reader.next() {
            scanned += 1;
            file_scanned += 1;
            let seqrec = record.with_context(|| {
                format!(
                    "failed to parse FASTA/FASTQ record #{file_scanned} from {}",
                    reads_path.display()
                )
            })?;

            let qname = read_name_token(seqrec.id());
            if qname.is_empty() || !candidates.contains(qname) || seen.contains(qname) {
                continue;
            }

            let seq = seqrec.seq();
            write_fasta(&mut writer, qname, seq.as_ref())?;
            seen.insert(qname.to_vec());
            written += 1;

            if seen.len() == candidates.len() {
                break;
            }

            if scanned % 1_000_000 == 0 {
                info!(
                    "[{}] fasta-extract progress scanned={scanned} written={written}",
                    elapsed()
                );
            }
        }

        if seen.len() == candidates.len() {
            break;
        }
    }

    writer.flush()?;

    let missing = candidates.len() as u64 - seen.len() as u64;
    if missing > 0 {
        warn!(
            "[{}] {missing} candidate name(s) were not found in {}",
            elapsed(),
            display_paths(reads_paths)
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

fn display_paths(paths: &[PathBuf]) -> String {
    paths
        .iter()
        .map(|path| path.display().to_string())
        .collect::<Vec<_>>()
        .join(",")
}

#[derive(Debug)]
struct BestBamSequence {
    seq: Vec<u8>,
    score: BamSeqScore,
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct BamSeqScore {
    no_hard_clip: bool,
    is_primary: bool,
    seq_len: usize,
}

fn has_hard_clip(record: &bam::Record) -> Result<bool> {
    for op in record.cigar().iter() {
        if op?.kind() == Kind::HardClip {
            return Ok(true);
        }
    }
    Ok(false)
}

fn read_name_token(header: &[u8]) -> &[u8] {
    let start = header
        .iter()
        .position(|b| !b.is_ascii_whitespace())
        .unwrap_or(header.len());
    let rest = &header[start..];
    let end = rest
        .iter()
        .position(|b| b.is_ascii_whitespace())
        .unwrap_or(rest.len());
    &rest[..end]
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

    use flate2::{write::GzEncoder, Compression};
    use noodles_sam::{
        self as sam,
        alignment::{
            io::Write as _,
            record::{
                cigar::{op::Kind as CigarKind, Op},
                Flags,
            },
            record_buf::{Cigar, QualityScores, Sequence},
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

    fn record_with_cigar(name: &[u8], seq: &[u8], flags: Flags, ops: &[Op]) -> RecordBuf {
        let cigar: Cigar = ops.iter().copied().collect();
        RecordBuf::builder()
            .set_name(name.to_vec())
            .set_flags(flags)
            .set_cigar(cigar)
            .set_sequence(Sequence::from(seq))
            .set_quality_scores(QualityScores::from(vec![30; seq.len()]))
            .build()
    }

    fn write_text(dir: &Path, name: &str, body: &[u8]) -> std::path::PathBuf {
        let path = dir.join(name);
        std::fs::write(&path, body).unwrap();
        path
    }

    fn write_gzip(dir: &Path, name: &str, body: &[u8]) -> std::path::PathBuf {
        let path = dir.join(name);
        let file = std::fs::File::create(&path).unwrap();
        let mut encoder = GzEncoder::new(file, Compression::default());
        std::io::Write::write_all(&mut encoder, body).unwrap();
        encoder.finish().unwrap();
        path
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
    fn bam_best_effort_prefers_full_primary_over_hard_clipped_supplementary() {
        let dir = tempdir().unwrap();
        let bam = write_bam(
            dir.path(),
            "reads.bam",
            &[
                record_with_cigar(
                    b"r1",
                    b"ACGT",
                    Flags::SUPPLEMENTARY,
                    &[
                        Op::new(CigarKind::HardClip, 2),
                        Op::new(CigarKind::Match, 4),
                    ],
                ),
                record_with_cigar(
                    b"r1",
                    b"TTACGTAA",
                    Flags::UNMAPPED,
                    &[
                        Op::new(CigarKind::SoftClip, 2),
                        Op::new(CigarKind::Match, 4),
                        Op::new(CigarKind::SoftClip, 2),
                    ],
                ),
            ],
        );
        let out = dir.path().join("out.fa");
        let cands = candset(&[b"r1"]);

        let stats =
            write_fasta_from_bam(&bam, &out, &cands, NonZeroUsize::new(1).unwrap()).unwrap();

        assert_eq!(stats.written, 1);
        assert_eq!(stats.missing, 0);
        assert_eq!(read_file(&out), b">r1\nTTACGTAA\n");
    }

    #[test]
    fn bam_best_effort_prefers_longer_equally_ranked_sequence() {
        let dir = tempdir().unwrap();
        let bam = write_bam(
            dir.path(),
            "reads.bam",
            &[
                record(b"r1", b"AC", Flags::UNMAPPED),
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
    fn reverse_complemented_best_record_emits_original_orientation() {
        let dir = tempdir().unwrap();
        let bam = write_bam(
            dir.path(),
            "reads.bam",
            &[
                record(b"r1", b"AA", Flags::UNMAPPED),
                record(
                    b"r1",
                    b"ACGA",
                    Flags::UNMAPPED | Flags::REVERSE_COMPLEMENTED,
                ),
            ],
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

    #[test]
    fn fasta_from_reads_extracts_candidates_by_first_header_token() {
        let dir = tempdir().unwrap();
        let reads = write_text(
            dir.path(),
            "reads.fa",
            b">r2 ignored\nTTTT\n>r1 full header text\nACGT\n",
        );
        let out = dir.path().join("out.fa");
        let cands = candset(&[b"r1", b"r3"]);

        let stats = write_fasta_from_reads(&[reads], &out, &cands).unwrap();

        assert_eq!(stats.written, 1);
        assert_eq!(stats.missing, 1);
        assert_eq!(stats.scanned, 2);
        assert_eq!(read_file(&out), b">r1\nACGT\n");
    }

    #[test]
    fn fastq_from_reads_extracts_exact_sequence() {
        let dir = tempdir().unwrap();
        let reads = write_text(
            dir.path(),
            "reads.fq",
            b"@r1 full header\nACGT\n+\n!!!!\n@r2\nTT\n+\n!!\n",
        );
        let out = dir.path().join("out.fa");
        let cands = candset(&[b"r1"]);

        let stats = write_fasta_from_reads(&[reads], &out, &cands).unwrap();

        assert_eq!(stats.written, 1);
        assert_eq!(stats.missing, 0);
        assert_eq!(stats.scanned, 1);
        assert_eq!(read_file(&out), b">r1\nACGT\n");
    }

    #[test]
    fn duplicate_read_id_emits_once_from_reads() {
        let dir = tempdir().unwrap();
        let reads = write_text(dir.path(), "reads.fa", b">r1\nAAAA\n>r1\nTTTT\n");
        let out = dir.path().join("out.fa");
        let cands = candset(&[b"r1", b"r2"]);

        let stats = write_fasta_from_reads(&[reads], &out, &cands).unwrap();

        assert_eq!(stats.written, 1);
        assert_eq!(stats.missing, 1);
        assert_eq!(stats.scanned, 2);
        assert_eq!(read_file(&out), b">r1\nAAAA\n");
    }

    #[test]
    fn reads_extraction_scans_multiple_files_in_order() {
        let dir = tempdir().unwrap();
        let reads1 = write_text(dir.path(), "reads1.fa", b">r1\nAAAA\n");
        let reads2 = write_text(dir.path(), "reads2.fa", b">r2\nCCCC\n");
        let out = dir.path().join("out.fa");
        let cands = candset(&[b"r1", b"r2"]);

        let stats = write_fasta_from_reads(&[reads1, reads2], &out, &cands).unwrap();

        assert_eq!(stats.written, 2);
        assert_eq!(stats.missing, 0);
        assert_eq!(stats.scanned, 2);
        assert_eq!(read_file(&out), b">r1\nAAAA\n>r2\nCCCC\n");
    }

    #[test]
    fn duplicate_read_id_across_files_emits_first_occurrence() {
        let dir = tempdir().unwrap();
        let reads1 = write_text(dir.path(), "reads1.fa", b">r1\nAAAA\n");
        let reads2 = write_text(dir.path(), "reads2.fa", b">r1\nTTTT\n>r2\nCCCC\n");
        let out = dir.path().join("out.fa");
        let cands = candset(&[b"r1", b"r2", b"r3"]);

        let stats = write_fasta_from_reads(&[reads1, reads2], &out, &cands).unwrap();

        assert_eq!(stats.written, 2);
        assert_eq!(stats.missing, 1);
        assert_eq!(stats.scanned, 3);
        assert_eq!(read_file(&out), b">r1\nAAAA\n>r2\nCCCC\n");
    }

    #[test]
    fn gzipped_fasta_from_reads_is_supported() {
        let dir = tempdir().unwrap();
        let reads = write_gzip(dir.path(), "reads.fa.gz", b">r1\nACGT\n");
        let out = dir.path().join("out.fa");
        let cands = candset(&[b"r1"]);

        let stats = write_fasta_from_reads(&[reads], &out, &cands).unwrap();

        assert_eq!(stats.written, 1);
        assert_eq!(stats.missing, 0);
        assert_eq!(stats.scanned, 1);
        assert_eq!(read_file(&out), b">r1\nACGT\n");
    }

    #[test]
    fn gzipped_fastq_from_reads_is_supported() {
        let dir = tempdir().unwrap();
        let reads = write_gzip(dir.path(), "reads.fq.gz", b"@r1\nACGT\n+\n!!!!\n");
        let out = dir.path().join("out.fa");
        let cands = candset(&[b"r1"]);

        let stats = write_fasta_from_reads(&[reads], &out, &cands).unwrap();

        assert_eq!(stats.written, 1);
        assert_eq!(stats.missing, 0);
        assert_eq!(stats.scanned, 1);
        assert_eq!(read_file(&out), b">r1\nACGT\n");
    }
}
