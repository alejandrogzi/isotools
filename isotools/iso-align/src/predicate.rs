// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Stage-2 predicate: decide whether a read's alignment set looks like a
//! pass-1 split that would benefit from re-alignment with a larger intron cap.
//!
//! A read is a candidate when *all* of:
//!
//! 1. it has at least two alignments,
//! 2. all alignments share the same reference and same strand,
//! 3. the largest reference-coord gap between adjacent segments (sorted by
//!    `ref_start`) is at least `min_gap` and no greater than `max_gap`,
//! 4. at least one flank clip on the side facing the gap is at least
//!    `min_flank_clip` bases (or anywhere on the read, with `FlankSide::Any`).
//!
//! The function is deliberately pure / branch-explicit so it can be unit-tested
//! against hand-built `MiniAln` vectors without touching any I/O.

use smallvec::SmallVec;

use crate::types::{FlankSide, MiniAln, PredicateConfig};

/// Returns `true` iff the alignment set is a re-align candidate.
pub fn is_candidate(alns: &[MiniAln], cfg: &PredicateConfig) -> bool {
    if alns.len() < 2 {
        return false;
    }

    let rname = alns[0].rname;
    if !alns.iter().all(|a| a.rname == rname) {
        return false;
    }

    let rev = alns[0].is_reverse;
    if !alns.iter().all(|a| a.is_reverse == rev) {
        return false;
    }

    let mut sorted: SmallVec<[&MiniAln; 4]> = alns.iter().collect();
    sorted.sort_unstable_by_key(|a| a.ref_start);

    let max_gap = sorted
        .windows(2)
        .map(|w| w[1].ref_start - w[0].ref_end)
        .max()
        .unwrap_or(0);

    if max_gap < cfg.min_gap as i64 {
        return false;
    }
    if let Some(mx) = cfg.max_gap {
        if max_gap > mx as i64 {
            return false;
        }
    }

    let min_clip = cfg.min_flank_clip;
    let clip_ok = match cfg.flank_side {
        FlankSide::Any => alns
            .iter()
            .any(|a| a.left_clip >= min_clip || a.right_clip >= min_clip),
        // For both strands the clip lengths are reported in reference
        // orientation by the CIGAR, so checking the upstream segment's right
        // clip and the downstream segment's left clip is correct regardless
        // of the read's strand.
        FlankSide::Inner => sorted
            .windows(2)
            .any(|w| w[0].right_clip >= min_clip || w[1].left_clip >= min_clip),
    };
    if !clip_ok {
        return false;
    }

    true
}

#[cfg(test)]
mod tests {
    use super::*;

    fn aln(rname: u32, start: i64, end: i64, rev: bool, lc: u32, rc: u32) -> MiniAln {
        MiniAln {
            rname,
            ref_start: start,
            ref_end: end,
            is_reverse: rev,
            left_clip: lc,
            right_clip: rc,
        }
    }

    fn cfg(min_gap: u64, max_gap: Option<u64>, min_clip: u32, side: FlankSide) -> PredicateConfig {
        PredicateConfig {
            min_gap,
            max_gap,
            min_flank_clip: min_clip,
            flank_side: side,
        }
    }

    #[test]
    fn single_alignment_rejected() {
        let alns = [aln(0, 100, 200, false, 50, 0)];
        assert!(!is_candidate(
            &alns,
            &cfg(300_000, None, 20, FlankSide::Inner)
        ));
    }

    #[test]
    fn cross_chromosome_rejected() {
        let alns = [
            aln(0, 100, 200, false, 0, 50),
            aln(1, 600_000, 600_100, false, 50, 0),
        ];
        assert!(!is_candidate(
            &alns,
            &cfg(300_000, None, 20, FlankSide::Inner)
        ));
    }

    #[test]
    fn opposite_strand_rejected() {
        let alns = [
            aln(0, 100, 200, false, 0, 50),
            aln(0, 600_000, 600_100, true, 50, 0),
        ];
        assert!(!is_candidate(
            &alns,
            &cfg(300_000, None, 20, FlankSide::Inner)
        ));
    }

    #[test]
    fn gap_below_min_rejected() {
        let alns = [
            aln(0, 100, 200, false, 0, 50),
            aln(0, 1_000, 1_100, false, 50, 0),
        ];
        assert!(!is_candidate(
            &alns,
            &cfg(300_000, None, 20, FlankSide::Inner)
        ));
    }

    #[test]
    fn gap_above_max_rejected() {
        let alns = [
            aln(0, 100, 200, false, 0, 50),
            aln(0, 2_000_200, 2_000_300, false, 50, 0),
        ];
        assert!(!is_candidate(
            &alns,
            &cfg(300_000, Some(1_600_000), 20, FlankSide::Inner)
        ));
    }

    #[test]
    fn clip_below_threshold_rejected_inner() {
        // Inner mode: only the upstream-right and downstream-left clips matter.
        let alns = [
            aln(0, 100, 200, false, 50, 0),
            aln(0, 600_000, 600_100, false, 0, 50),
        ];
        assert!(!is_candidate(
            &alns,
            &cfg(300_000, None, 20, FlankSide::Inner)
        ));
    }

    #[test]
    fn clip_below_threshold_inner_passes_under_any() {
        // Same record set as above — Any mode picks up the outer clips.
        let alns = [
            aln(0, 100, 200, false, 50, 0),
            aln(0, 600_000, 600_100, false, 0, 50),
        ];
        assert!(is_candidate(&alns, &cfg(300_000, None, 20, FlankSide::Any)));
    }

    #[test]
    fn valid_candidate_inner() {
        let alns = [
            aln(0, 100, 200, false, 0, 50),
            aln(0, 600_000, 600_100, false, 50, 0),
        ];
        assert!(is_candidate(
            &alns,
            &cfg(300_000, None, 20, FlankSide::Inner)
        ));
    }

    #[test]
    fn unsorted_input_handled() {
        // Same as `valid_candidate_inner` but with the alignments reversed in
        // the input; the predicate must sort by ref_start internally.
        let alns = [
            aln(0, 600_000, 600_100, false, 50, 0),
            aln(0, 100, 200, false, 0, 50),
        ];
        assert!(is_candidate(
            &alns,
            &cfg(300_000, None, 20, FlankSide::Inner)
        ));
    }

    #[test]
    fn three_segments_picks_max_gap() {
        // Two small gaps, one large — predicate looks at the *max* gap.
        let alns = [
            aln(0, 100, 200, false, 0, 50),
            aln(0, 1_000, 1_100, false, 5, 30),
            aln(0, 600_000, 600_100, false, 50, 0),
        ];
        assert!(is_candidate(
            &alns,
            &cfg(300_000, None, 20, FlankSide::Inner)
        ));
    }
}
