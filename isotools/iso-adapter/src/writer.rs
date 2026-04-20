// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Adapter trimming of BAM records.
//!
//! Trimming only touches the soft-clip portion of a read. The aligned core
//! is never altered, so `pos`, aligned CIGAR ops, and all BAM tags other
//! than SEQ/QUAL/CIGAR remain untouched.
//!
//! The function `trim_record` takes a record, the end it should trim on,
//! the length of that end's soft-clip, and the `AdapterMatch` (whose
//! `clip_range` is relative to the clip). It removes only the adapter bases;
//! any non-adapter overhang in the clip survives as a shorter soft-clip.

use noodles_sam::alignment::record::cigar::{op::Kind, Op};
use noodles_sam::alignment::record_buf::{Cigar as CigarBuf, QualityScores as QualityBuf, Sequence};
use noodles_sam::alignment::RecordBuf;

use crate::detector::{AdapterMatch, ClipEnd};

/// Outcome of a single `trim_record` call.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TrimOutcome {
    /// The adapter was removed and the record was updated in place.
    Trimmed { bases_removed: usize },
    /// The match was degenerate (zero-length, or inconsistent with the clip)
    /// and the record was left untouched.
    Unchanged,
}

/// Trims the adapter region out of a record's soft-clip without touching
/// the aligned core.
pub fn trim_record(
    record: &mut RecordBuf,
    end: ClipEnd,
    clip_len: usize,
    m: &AdapterMatch,
) -> TrimOutcome {
    if clip_len == 0
        || m.clip_range.start >= m.clip_range.end
        || m.clip_range.end > clip_len
    {
        return TrimOutcome::Unchanged;
    }

    let bases_removed = m.clip_range.end - m.clip_range.start;
    if bases_removed == 0 {
        return TrimOutcome::Unchanged;
    }

    let seq_len = record.sequence().as_ref().len();
    if clip_len > seq_len {
        return TrimOutcome::Unchanged;
    }

    // Derive absolute [removed_start, removed_end) within the read sequence.
    let (removed_start, removed_end) = match end {
        ClipEnd::FivePrime => (m.clip_range.start, m.clip_range.end),
        ClipEnd::ThreePrime => {
            let offset = seq_len - clip_len;
            (offset + m.clip_range.start, offset + m.clip_range.end)
        }
    };

    // SEQ / QUAL.
    let new_seq = splice_out(record.sequence().as_ref(), removed_start, removed_end);
    let quals = record.quality_scores().as_ref();
    let new_qual = if quals.is_empty() {
        Vec::new()
    } else if quals.len() == seq_len {
        splice_out(quals, removed_start, removed_end)
    } else {
        quals.to_vec()
    };

    // CIGAR: shorten the affected soft-clip by `bases_removed`, dropping it
    // entirely when its new length reaches zero.
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

/// Shrinks the leading / trailing soft-clip op by `bases_removed`.
///
/// Hard-clip operations surrounding the soft-clip are preserved. If the
/// soft-clip's new length is zero, the op is removed entirely. Returns
/// `None` when the CIGAR has no soft-clip on the requested end or the
/// removal would exceed the clip's length.
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

    fn record(
        cigar: Vec<Op>,
        seq: &[u8],
        qual: &[u8],
    ) -> RecordBuf {
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
    fn trim_5p_adapter_at_start() {
        let mut rec = record(
            vec![
                Op::new(Kind::SoftClip, 10),
                Op::new(Kind::Match, 5),
            ],
            b"ADAPTERCLPMATCH",
            b"!!!!!!!!!!#####",
        );
        let m = AdapterMatch {
            label: "x",
            clip_range: 0..7,
            edit_distance: 0,
            on_reverse_strand: false,
        };
        let outcome = trim_record(&mut rec, ClipEnd::FivePrime, 10, &m);
        assert!(matches!(outcome, TrimOutcome::Trimmed { bases_removed: 7 }));
        assert_eq!(rec.sequence().as_ref(), b"CLPMATCH");
        assert_eq!(rec.quality_scores().as_ref(), b"!!!#####");
        assert_eq!(
            cigar_of(&rec),
            vec![(Kind::SoftClip, 3), (Kind::Match, 5)]
        );
    }

    #[test]
    fn trim_5p_removes_entire_softclip() {
        let mut rec = record(
            vec![
                Op::new(Kind::SoftClip, 5),
                Op::new(Kind::Match, 3),
            ],
            b"AAAAAMMM",
            b"!!!!!!!!",
        );
        let m = AdapterMatch {
            label: "x",
            clip_range: 0..5,
            edit_distance: 0,
            on_reverse_strand: false,
        };
        let outcome = trim_record(&mut rec, ClipEnd::FivePrime, 5, &m);
        assert!(matches!(outcome, TrimOutcome::Trimmed { bases_removed: 5 }));
        assert_eq!(rec.sequence().as_ref(), b"MMM");
        assert_eq!(cigar_of(&rec), vec![(Kind::Match, 3)]);
    }

    #[test]
    fn trim_3p_adapter_at_end_preserves_hard_clip() {
        let mut rec = record(
            vec![
                Op::new(Kind::Match, 4),
                Op::new(Kind::SoftClip, 8),
                Op::new(Kind::HardClip, 2),
            ],
            b"MMMMCCAAAAAA",
            b"!!!!????????",
        );
        let m = AdapterMatch {
            label: "x",
            clip_range: 2..8,
            edit_distance: 0,
            on_reverse_strand: false,
        };
        let outcome = trim_record(&mut rec, ClipEnd::ThreePrime, 8, &m);
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
    fn trim_ignores_zero_length_range() {
        let mut rec = record(
            vec![Op::new(Kind::SoftClip, 4), Op::new(Kind::Match, 2)],
            b"ACGTGG",
            b"!!!!!!",
        );
        let m = AdapterMatch {
            label: "x",
            clip_range: 2..2,
            edit_distance: 0,
            on_reverse_strand: false,
        };
        let outcome = trim_record(&mut rec, ClipEnd::FivePrime, 4, &m);
        assert_eq!(outcome, TrimOutcome::Unchanged);
        assert_eq!(rec.sequence().as_ref(), b"ACGTGG");
    }

    #[test]
    fn trim_rejects_out_of_bounds() {
        let mut rec = record(
            vec![Op::new(Kind::SoftClip, 3), Op::new(Kind::Match, 2)],
            b"AAAGG",
            b"!!!!!",
        );
        let m = AdapterMatch {
            label: "x",
            clip_range: 0..5,
            edit_distance: 0,
            on_reverse_strand: false,
        };
        let outcome = trim_record(&mut rec, ClipEnd::FivePrime, 3, &m);
        assert_eq!(outcome, TrimOutcome::Unchanged);
        assert_eq!(rec.sequence().as_ref(), b"AAAGG");
    }
}
