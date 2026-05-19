// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Stage-1 scan of the pass-1 BAM into a per-QNAME map of [`MiniAln`].
//!
//! Bgzf decompression is the dominant cost when scanning a long-read BAM at
//! this rate, so we wrap the file in [`bgzf::MultithreadedReader`] and let
//! the caller pick the worker count. Record consumption stays single-threaded
//! because parallelizing the per-QNAME hashmap insert would cost more in
//! contention than it gains. We use the lazy `bam::Record` so that SEQ / QUAL
//! / aux are never materialized — only flags, position and CIGAR are touched.

use std::fs::File;
use std::io;
use std::num::NonZeroUsize;
use std::path::Path;

use anyhow::{Context, Result};
use log::info;
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::alignment::record::cigar::op::Kind;
use rustc_hash::FxHashMap;
use smallvec::SmallVec;

use crate::logging::elapsed;
use crate::types::MiniAln;

/// Per-QNAME alignment list. Two alignments is the common case for the
/// pass-1-split signal we're after, so a small inline buffer of 2 fits most
/// reads without a heap allocation.
pub type AlnGroups = FxHashMap<Vec<u8>, SmallVec<[MiniAln; 2]>>;

/// Scans the BAM at `path` and returns per-QNAME alignment groups.
///
/// `bgzf_workers` controls the size of the bgzf decompression worker pool.
/// `1` keeps decompression on the calling thread; values >1 spin up that many
/// background workers. Unmapped, secondary and supplementary records are
/// dropped. QNAMEs are stored as raw bytes (never UTF-8-validated) — this
/// is both faster and strictly correct, since SAM permits any printable ASCII.
pub fn scan(path: &Path, bgzf_workers: NonZeroUsize) -> Result<AlnGroups> {
    let file =
        File::open(path).with_context(|| format!("failed to open BAM {}", path.display()))?;
    let bgzf_reader = bgzf::MultithreadedReader::with_worker_count(bgzf_workers, file);
    let mut reader = bam::io::Reader::from(bgzf_reader);
    let _header = reader
        .read_header()
        .with_context(|| format!("failed to read BAM header from {}", path.display()))?;

    info!(
        "[{}] bam-scan begin path={} bgzf_workers={}",
        elapsed(),
        path.display(),
        bgzf_workers.get()
    );

    let mut groups: AlnGroups = FxHashMap::default();
    let mut record = bam::Record::default();
    let mut total: u64 = 0;
    let mut kept: u64 = 0;

    loop {
        let bytes = reader
            .read_record(&mut record)
            .with_context(|| format!("failed to read BAM record #{total}"))?;
        if bytes == 0 {
            break;
        }
        total += 1;

        let flags = record.flags();
        if flags.is_unmapped() || flags.is_secondary() {
            continue;
        }

        let Some(rid_res) = record.reference_sequence_id() else {
            continue;
        };
        let rid =
            rid_res.with_context(|| format!("invalid reference id at record #{total}"))? as u32;

        let Some(start_res) = record.alignment_start() else {
            continue;
        };
        let start = usize::from(
            start_res.with_context(|| format!("invalid alignment_start at record #{total}"))?,
        ) as i64;

        let (left_clip, right_clip, ref_consumed) = cigar_summary(record.cigar().iter())
            .with_context(|| format!("invalid CIGAR at record #{total}"))?;
        let ref_end = start + ref_consumed;

        let Some(name) = record.name() else {
            // Unmapped reads can lack QNAME, but we already filtered them.
            // A mapped record without a name is malformed; skip defensively.
            continue;
        };
        let qname: &[u8] = name.as_ref();

        let aln = MiniAln {
            rname: rid,
            ref_start: start,
            ref_end,
            is_reverse: flags.is_reverse_complemented(),
            left_clip,
            right_clip,
        };

        groups.entry(qname.to_vec()).or_default().push(aln);
        kept += 1;

        if total % 1_000_000 == 0 {
            info!(
                "[{}] bam-scan progress total={total} kept={kept} groups={}",
                elapsed(),
                groups.len()
            );
        }
    }

    info!(
        "[{}] bam-scan done total={total} kept={kept} distinct_qnames={}",
        elapsed(),
        groups.len()
    );
    Ok(groups)
}

/// Walks a CIGAR op stream once and returns
/// `(leading_clip, trailing_clip, ref_consumed)`.
///
/// Leading / trailing clips include both `S` and `H`. Reference-consuming
/// ops are `M`, `D`, `N`, `=`, `X` per the SAM spec.
pub fn cigar_summary<I>(ops: I) -> io::Result<(u32, u32, i64)>
where
    I: IntoIterator<Item = io::Result<noodles_sam::alignment::record::cigar::Op>>,
{
    let mut left_clip: u32 = 0;
    let mut right_clip: u32 = 0;
    let mut ref_consumed: i64 = 0;
    let mut seen_body = false;

    for op in ops {
        let op = op?;
        let len = op.len() as u32;
        match op.kind() {
            Kind::SoftClip | Kind::HardClip => {
                if !seen_body {
                    left_clip = left_clip.saturating_add(len);
                } else {
                    right_clip = right_clip.saturating_add(len);
                }
            }
            other => {
                seen_body = true;
                // Reset the trailing-clip accumulator: a clip that turned out
                // to be followed by a body op is *not* a trailing clip.
                right_clip = 0;
                if matches!(
                    other,
                    Kind::Match
                        | Kind::Deletion
                        | Kind::Skip
                        | Kind::SequenceMatch
                        | Kind::SequenceMismatch
                ) {
                    ref_consumed += len as i64;
                }
            }
        }
    }

    Ok((left_clip, right_clip, ref_consumed))
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam::alignment::record::cigar::Op;

    fn ops(items: &[(Kind, usize)]) -> Vec<io::Result<Op>> {
        items.iter().map(|(k, l)| Ok(Op::new(*k, *l))).collect()
    }

    #[test]
    fn no_clips_only_match() {
        let summary = cigar_summary(ops(&[(Kind::Match, 100)])).unwrap();
        assert_eq!(summary, (0, 0, 100));
    }

    #[test]
    fn soft_clip_both_ends() {
        let summary = cigar_summary(ops(&[
            (Kind::SoftClip, 5),
            (Kind::Match, 90),
            (Kind::SoftClip, 7),
        ]))
        .unwrap();
        assert_eq!(summary, (5, 7, 90));
    }

    #[test]
    fn hard_then_soft_clip_at_start() {
        let summary = cigar_summary(ops(&[
            (Kind::HardClip, 3),
            (Kind::SoftClip, 5),
            (Kind::Match, 90),
            (Kind::SoftClip, 2),
            (Kind::HardClip, 1),
        ]))
        .unwrap();
        assert_eq!(summary, (3 + 5, 2 + 1, 90));
    }

    #[test]
    fn skip_and_deletion_consume_reference() {
        // 5S 30M 100N 20M 5S → ref_consumed = 30 + 100 + 20 = 150, clips 5/5.
        let summary = cigar_summary(ops(&[
            (Kind::SoftClip, 5),
            (Kind::Match, 30),
            (Kind::Skip, 100),
            (Kind::Match, 20),
            (Kind::SoftClip, 5),
        ]))
        .unwrap();
        assert_eq!(summary, (5, 5, 150));
    }

    #[test]
    fn insertion_does_not_consume_reference() {
        let summary = cigar_summary(ops(&[
            (Kind::Match, 50),
            (Kind::Insertion, 10),
            (Kind::Match, 50),
        ]))
        .unwrap();
        assert_eq!(summary, (0, 0, 100));
    }

    #[test]
    fn equal_and_diff_count_as_reference_consuming() {
        let summary = cigar_summary(ops(&[
            (Kind::SequenceMatch, 60),
            (Kind::SequenceMismatch, 5),
            (Kind::SequenceMatch, 35),
        ]))
        .unwrap();
        assert_eq!(summary, (0, 0, 100));
    }

    #[test]
    fn empty_cigar() {
        let summary =
            cigar_summary::<std::iter::Empty<io::Result<Op>>>(std::iter::empty()).unwrap();
        assert_eq!(summary, (0, 0, 0));
    }
}
