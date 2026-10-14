// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Stage-1 scan of the pass-1 BAM into a per-QNAME map of [`MiniAln`].
//!
//! Bgzf decompression is the dominant cost when scanning a long-read BAM at
//! this rate, so we wrap the file in [`bgzf::MultithreadedReader`] and let
//! the caller pick the worker count. Record consumption stays single-threaded
//! because parallelizing the per-QNAME hashmap insert would cost more in
//! contention than it gains. We use the lazy `bam::Record` so that SEQ / QUAL
//! are never materialized — only flags, position, CIGAR and the presence of
//! the `SA` tag are touched.

use std::fs::File;
use std::io;
use std::num::NonZeroUsize;
use std::path::Path;

use anyhow::{Context, Result};
use log::info;
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::data::field::Tag;
use rustc_hash::FxHashMap;
use smallvec::SmallVec;

use crate::logging::elapsed;
use crate::types::{MiniAln, ScanConfig};

/// Per-QNAME alignment list. Two alignments is the common case for the
/// pass-1-split signal we're after, so a small inline buffer of 2 fits most
/// reads without a heap allocation.
pub type AlnGroups = FxHashMap<Vec<u8>, SmallVec<[MiniAln; 2]>>;

/// Scans the BAM at `path` and returns per-QNAME alignment groups.
///
/// `bgzf_workers` controls the size of the bgzf decompression worker pool.
/// `1` keeps decompression on the calling thread; values >1 spin up that many
/// background workers. Unmapped and secondary records are dropped;
/// supplementary records are kept (the split rule needs them). QNAMEs are
/// stored as raw bytes (never UTF-8-validated) — this is both faster and
/// strictly correct, since SAM permits any printable ASCII.
pub fn scan(path: &Path, bgzf_workers: NonZeroUsize, cfg: &ScanConfig) -> Result<AlnGroups> {
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

        let feats = cigar_features(record.cigar().iter(), cfg)
            .with_context(|| format!("invalid CIGAR at record #{total}"))?;
        let ref_end = start + feats.ref_consumed;

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
            left_clip: feats.left_clip,
            right_clip: feats.right_clip,
            is_supplementary: flags.is_supplementary(),
            has_sa: record.data().get(&Tag::OTHER_ALIGNMENTS).is_some(),
            max_skip: feats.max_skip,
            ins_near_skip: feats.ins_near_skip,
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

/// Single-pass CIGAR summary kept per record.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct CigarFeatures {
    pub left_clip: u32,
    pub right_clip: u32,
    pub ref_consumed: i64,
    pub max_skip: u32,
    pub ins_near_skip: bool,
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
    let f = cigar_features(ops, &ScanConfig::disabled())?;
    Ok((f.left_clip, f.right_clip, f.ref_consumed))
}

/// Walks a CIGAR op stream once and returns clips, reference span, the longest
/// `N`, and whether an insertion sits next to a long `N`.
///
/// Insertion test: each insertion of at least `cfg.ins_min_len` is paired with
/// its nearest `N` (distance = min over both `N` boundaries of
/// |insertion ref position - boundary|; the first `N` wins ties). The test
/// passes when that distance is at most `cfg.ins_max_dist` and that `N` is at
/// least `cfg.ins_min_skip` long. A record without `N` never passes.
pub fn cigar_features<I>(ops: I, cfg: &ScanConfig) -> io::Result<CigarFeatures>
where
    I: IntoIterator<Item = io::Result<noodles_sam::alignment::record::cigar::Op>>,
{
    let mut left_clip: u32 = 0;
    let mut right_clip: u32 = 0;
    let mut ref_consumed: i64 = 0;
    let mut max_skip: u32 = 0;
    let mut seen_body = false;
    // INFO: reference offsets of long-enough insertions and of every N.
    let mut insertions: SmallVec<[i64; 4]> = SmallVec::new();
    let mut skips: SmallVec<[(i64, i64); 16]> = SmallVec::new();

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
                match other {
                    Kind::Skip => {
                        skips.push((ref_consumed, ref_consumed + len as i64));
                        max_skip = max_skip.max(len);
                        ref_consumed += len as i64;
                    }
                    Kind::Insertion => {
                        if len >= cfg.ins_min_len {
                            insertions.push(ref_consumed);
                        }
                    }
                    Kind::Match | Kind::Deletion | Kind::SequenceMatch | Kind::SequenceMismatch => {
                        ref_consumed += len as i64;
                    }
                    _ => {}
                }
            }
        }
    }

    let ins_near_skip = insertions.iter().any(|&p| {
        let mut best: Option<(i64, i64)> = None; // (distance, N length)
        for &(s, e) in &skips {
            let d = (p - s).abs().min((p - e).abs());
            if best.map_or(true, |(bd, _)| d < bd) {
                best = Some((d, e - s));
            }
        }
        matches!(best, Some((d, n)) if d <= cfg.ins_max_dist as i64 && n >= cfg.ins_min_skip as i64)
    });

    Ok(CigarFeatures {
        left_clip,
        right_clip,
        ref_consumed,
        max_skip,
        ins_near_skip,
    })
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

    fn ins_cfg() -> ScanConfig {
        ScanConfig {
            ins_min_len: 20,
            ins_max_dist: 10,
            ins_min_skip: 20_000,
        }
    }

    #[test]
    fn bridge_insertion_next_to_long_skip() {
        // Kcnh7-like: 264,397N + 7-nt pseudo-exon + 42,573N, with the skipped
        // exon left as a 557-nt insertion right after the second N.
        let f = cigar_features(
            ops(&[
                (Kind::Match, 300),
                (Kind::Skip, 264_397),
                (Kind::Match, 7),
                (Kind::Skip, 42_573),
                (Kind::Insertion, 557),
                (Kind::Match, 900),
            ]),
            &ins_cfg(),
        )
        .unwrap();
        assert!(f.ins_near_skip);
        assert_eq!(f.max_skip, 264_397);
        assert_eq!(f.ref_consumed, 300 + 264_397 + 7 + 42_573 + 900);
    }

    #[test]
    fn insertion_too_short_or_too_far_or_short_skip() {
        let short = ops(&[
            (Kind::Match, 100),
            (Kind::Skip, 50_000),
            (Kind::Insertion, 19),
            (Kind::Match, 100),
        ]);
        assert!(!cigar_features(short, &ins_cfg()).unwrap().ins_near_skip);

        let far = ops(&[
            (Kind::Match, 100),
            (Kind::Skip, 50_000),
            (Kind::Match, 11),
            (Kind::Insertion, 40),
            (Kind::Match, 100),
        ]);
        assert!(!cigar_features(far, &ins_cfg()).unwrap().ins_near_skip);

        let at_ten = ops(&[
            (Kind::Match, 100),
            (Kind::Skip, 50_000),
            (Kind::Match, 10),
            (Kind::Insertion, 40),
            (Kind::Match, 100),
        ]);
        assert!(cigar_features(at_ten, &ins_cfg()).unwrap().ins_near_skip);

        let small_skip = ops(&[
            (Kind::Match, 100),
            (Kind::Skip, 19_999),
            (Kind::Insertion, 40),
            (Kind::Match, 100),
        ]);
        assert!(
            !cigar_features(small_skip, &ins_cfg())
                .unwrap()
                .ins_near_skip
        );

        let no_skip = ops(&[
            (Kind::Match, 100),
            (Kind::Insertion, 40),
            (Kind::Match, 100),
        ]);
        assert!(!cigar_features(no_skip, &ins_cfg()).unwrap().ins_near_skip);
    }

    #[test]
    fn insertion_judged_by_its_nearest_skip_only() {
        // The nearest N (3 nt away) is short; a long N 200 nt away does not count.
        let f = cigar_features(
            ops(&[
                (Kind::Match, 100),
                (Kind::Skip, 100_000),
                (Kind::Match, 200),
                (Kind::Insertion, 30),
                (Kind::Match, 3),
                (Kind::Skip, 500),
                (Kind::Match, 100),
            ]),
            &ins_cfg(),
        )
        .unwrap();
        assert!(!f.ins_near_skip);
    }

    #[test]
    fn insertion_tie_goes_to_first_skip() {
        // Insertion exactly between two Ns, 5 nt from each; the first N is long.
        let f = cigar_features(
            ops(&[
                (Kind::Match, 100),
                (Kind::Skip, 30_000),
                (Kind::Match, 5),
                (Kind::Insertion, 25),
                (Kind::Match, 5),
                (Kind::Skip, 400),
                (Kind::Match, 100),
            ]),
            &ins_cfg(),
        )
        .unwrap();
        assert!(f.ins_near_skip);
    }

    #[test]
    fn empty_cigar() {
        let summary =
            cigar_summary::<std::iter::Empty<io::Result<Op>>>(std::iter::empty()).unwrap();
        assert_eq!(summary, (0, 0, 0));
    }
}
