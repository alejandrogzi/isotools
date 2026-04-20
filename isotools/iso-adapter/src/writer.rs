// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Trimming of a soft-clip range from a BAM record.
//!
//! This module is deliberately *policy-free*: the caller computes which
//! absolute range of clip bases must be spliced out, and this module
//! performs the mechanical edit of `SEQ`, `QUAL` and `CIGAR`. The aligned
//! core is never altered, so `POS` and aligned CIGAR ops are preserved.
//!
//! The caller's policy (adapter-outward vs. polyA-keep-inward) lives in
//! `engine::plan_trim`.

use noodles_sam::alignment::record::cigar::{op::Kind, Op};
use noodles_sam::alignment::record_buf::{Cigar as CigarBuf, QualityScores as QualityBuf, Sequence};
use noodles_sam::alignment::RecordBuf;

use crate::detector::ClipEnd;

/// Outcome of a single `trim_clip` call.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TrimOutcome {
    /// Bases were removed and the record was updated in place.
    Trimmed { bases_removed: usize },
    /// The cut range was degenerate (empty / out-of-bounds) and the record
    /// was left untouched.
    Unchanged,
}

/// Splice `[cut_start, cut_end)` out of the given end's soft-clip.
///
/// * `end` — which clip the cut is relative to.
/// * `clip_len` — length of that clip before any changes.
/// * `cut_start`, `cut_end` — half-open range *within the clip* to remove.
///
/// On success, the surrounding hard-clips are preserved and the soft-clip op
/// is shrunk (or removed entirely when it reaches zero length).
pub fn trim_clip(
    record: &mut RecordBuf,
    end: ClipEnd,
    clip_len: usize,
    cut_start: usize,
    cut_end: usize,
) -> TrimOutcome {
    if clip_len == 0 || cut_start >= cut_end || cut_end > clip_len {
        return TrimOutcome::Unchanged;
    }

    let bases_removed = cut_end - cut_start;
    let seq_len = record.sequence().as_ref().len();
    if clip_len > seq_len {
        return TrimOutcome::Unchanged;
    }

    let (removed_start, removed_end) = match end {
        ClipEnd::FivePrime => (cut_start, cut_end),
        ClipEnd::ThreePrime => {
            let offset = seq_len - clip_len;
            (offset + cut_start, offset + cut_end)
        }
    };

    let new_seq = splice_out(record.sequence().as_ref(), removed_start, removed_end);
    let quals = record.quality_scores().as_ref();
    let new_qual = if quals.is_empty() {
        Vec::new()
    } else if quals.len() == seq_len {
        splice_out(quals, removed_start, removed_end)
    } else {
        quals.to_vec()
    };

    let old_ops: Vec<Op> = record.cigar().as_ref().to_vec();
    let new_ops = match adjust_softclip(&old_ops, end, bases_removed) {
        Some(ops) => ops,
        None => return TrimOutcome::Unchanged,
    };

    *record.sequence_mut() = Sequence::from(new_seq);
    *record.quality_scores_mut() = QualityBuf::from(new_qual);
    *record.cigar_mut() = CigarBuf::from(new_ops);

    TrimOutcome::Trimmed { bases_removed }
}

/// Returns a new `Vec<u8>` with bytes in `[start, end)` removed.
fn splice_out(input: &[u8], start: usize, end: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(input.len() - (end - start));
    out.extend_from_slice(&input[..start]);
    out.extend_from_slice(&input[end..]);
    out
}

/// Shrinks the leading or trailing soft-clip op by `bases_removed`.
///
/// Surrounding hard-clip ops are preserved. If the soft-clip's new length is
/// zero, the op is dropped. Returns `None` when no soft-clip is present on
/// the requested end or when the removal would exceed the clip's length.
fn adjust_softclip(ops: &[Op], end: ClipEnd, bases_removed: usize) -> Option<Vec<Op>> {
    if ops.is_empty() {
        return None;
    }

    let mut out = ops.to_vec();
    let range: Box<dyn Iterator<Item = usize>> = match end {
        ClipEnd::FivePrime => Box::new(0..out.len()),
        ClipEnd::ThreePrime => Box::new((0..out.len()).rev()),
    };

    for idx in range {
        match out[idx].kind() {
            Kind::HardClip => continue,
            Kind::SoftClip => {
                let len = out[idx].len();
                if bases_removed > len {
                    return None;
                }
                let new_len = len - bases_removed;
                if new_len == 0 {
                    out.remove(idx);
                } else {
                    out[idx] = Op::new(Kind::SoftClip, new_len);
                }
                return Some(out);
            }
            _ => return None,
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(cigar: Vec<Op>, seq: &[u8], qual: &[u8]) -> RecordBuf {
        let mut rec = RecordBuf::default();
        *rec.cigar_mut() = CigarBuf::from(cigar);
        *rec.sequence_mut() = Sequence::from(seq.to_vec());
        *rec.quality_scores_mut() = QualityBuf::from(qual.to_vec());
        rec
    }

    fn cigar_of(rec: &RecordBuf) -> Vec<(Kind, usize)> {
        rec.cigar()
            .as_ref()
            .iter()
            .map(|op| (op.kind(), op.len()))
            .collect()
    }

    #[test]
    fn trim_5p_outward_range_removes_leading_bytes() {
        let mut rec = record(
            vec![Op::new(Kind::SoftClip, 10), Op::new(Kind::Match, 5)],
            b"ADAPTERCLPMATCH",
            b"!!!!!!!!!!#####",
        );
        // Cut [0 .. 7) within the 5' clip (adapter + any outward overhang).
        let outcome = trim_clip(&mut rec, ClipEnd::FivePrime, 10, 0, 7);
        assert!(matches!(outcome, TrimOutcome::Trimmed { bases_removed: 7 }));
        assert_eq!(rec.sequence().as_ref(), b"CLPMATCH");
        assert_eq!(rec.quality_scores().as_ref(), b"!!!#####");
        assert_eq!(
            cigar_of(&rec),
            vec![(Kind::SoftClip, 3), (Kind::Match, 5)]
        );
    }

    #[test]
    fn trim_5p_entire_softclip_drops_op() {
        let mut rec = record(
            vec![Op::new(Kind::SoftClip, 5), Op::new(Kind::Match, 3)],
            b"AAAAAMMM",
            b"!!!!!!!!",
        );
        let outcome = trim_clip(&mut rec, ClipEnd::FivePrime, 5, 0, 5);
        assert!(matches!(outcome, TrimOutcome::Trimmed { bases_removed: 5 }));
        assert_eq!(rec.sequence().as_ref(), b"MMM");
        assert_eq!(cigar_of(&rec), vec![(Kind::Match, 3)]);
    }

    #[test]
    fn trim_3p_outward_range_preserves_hard_clip() {
        let mut rec = record(
            vec![
                Op::new(Kind::Match, 4),
                Op::new(Kind::SoftClip, 8),
                Op::new(Kind::HardClip, 2),
            ],
            b"MMMMCCAAAAAA",
            b"!!!!????????",
        );
        // 3' adapter trim semantics: cut [2 .. 8) inside the 3' clip
        // (adapter and everything outward).
        let outcome = trim_clip(&mut rec, ClipEnd::ThreePrime, 8, 2, 8);
        assert!(matches!(outcome, TrimOutcome::Trimmed { bases_removed: 6 }));
        assert_eq!(rec.sequence().as_ref(), b"MMMMCC");
        assert_eq!(rec.quality_scores().as_ref(), b"!!!!??");
        assert_eq!(
            cigar_of(&rec),
            vec![
                (Kind::Match, 4),
                (Kind::SoftClip, 2),
                (Kind::HardClip, 2)
            ]
        );
    }

    #[test]
    fn trim_3p_polya_keeps_run_trims_outward() {
        // Read:  MMMM | AAAAAAAATTTT  (3' clip = AAAAAAAATTTT, polyA at [0..8))
        // polyA keep-inward policy: cut [8..12) — outward of the run.
        let mut rec = record(
            vec![Op::new(Kind::Match, 4), Op::new(Kind::SoftClip, 12)],
            b"MMMMAAAAAAAATTTT",
            b"!!!!############",
        );
        let outcome = trim_clip(&mut rec, ClipEnd::ThreePrime, 12, 8, 12);
        assert!(matches!(outcome, TrimOutcome::Trimmed { bases_removed: 4 }));
        assert_eq!(rec.sequence().as_ref(), b"MMMMAAAAAAAA");
        assert_eq!(cigar_of(&rec), vec![(Kind::Match, 4), (Kind::SoftClip, 8)]);
    }

    #[test]
    fn trim_ignores_zero_length_range() {
        let mut rec = record(
            vec![Op::new(Kind::SoftClip, 4), Op::new(Kind::Match, 2)],
            b"ACGTGG",
            b"!!!!!!",
        );
        let outcome = trim_clip(&mut rec, ClipEnd::FivePrime, 4, 2, 2);
        assert_eq!(outcome, TrimOutcome::Unchanged);
        assert_eq!(rec.sequence().as_ref(), b"ACGTGG");
    }

    #[test]
    fn trim_rejects_out_of_bounds_cut() {
        let mut rec = record(
            vec![Op::new(Kind::SoftClip, 3), Op::new(Kind::Match, 2)],
            b"AAAGG",
            b"!!!!!",
        );
        let outcome = trim_clip(&mut rec, ClipEnd::FivePrime, 3, 0, 5);
        assert_eq!(outcome, TrimOutcome::Unchanged);
        assert_eq!(rec.sequence().as_ref(), b"AAAGG");
    }
}
