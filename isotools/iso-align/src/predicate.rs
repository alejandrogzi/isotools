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
//! That is the `split` rule. [`candidate_rules`] adds three rules that look at
//! the primary record alone (pass-1 records that were not split):
//!
//! - `insertion`: an insertion next to a long `N` (bridge pseudo-exon: the
//!   true exon is left as an insertion between two long introns),
//! - `long-intron-clip`: a clip on either end and a very long `N`,
//! - `orphan-clip`: a long clip and no `SA` tag.
//!
//! The functions are deliberately pure / branch-explicit so they can be
//! unit-tested against hand-built `MiniAln` vectors without touching any I/O.

use smallvec::SmallVec;

use crate::types::{FlankSide, MiniAln, PredicateConfig, Rule, RuleSet};

/// Returns the enabled rules that fire for this read's alignment set; empty
/// when the read is not a candidate.
pub fn candidate_rules(alns: &[MiniAln], cfg: &PredicateConfig) -> RuleSet {
    let mut fired = RuleSet::default();
    if cfg.rules.contains(Rule::Split) && is_candidate(alns, cfg) {
        fired.insert(Rule::Split);
    }

    let Some(primary) = alns.iter().find(|a| !a.is_supplementary) else {
        return fired;
    };
    let max_clip = primary.left_clip.max(primary.right_clip);

    if cfg.rules.contains(Rule::Insertion) && primary.ins_near_skip {
        fired.insert(Rule::Insertion);
    }
    if cfg.rules.contains(Rule::LongIntronClip)
        && primary.max_skip >= cfg.long_intron_min_skip
        && max_clip >= cfg.long_intron_min_clip
    {
        fired.insert(Rule::LongIntronClip);
    }
    if cfg.rules.contains(Rule::OrphanClip) && !primary.has_sa && max_clip >= cfg.orphan_min_clip {
        fired.insert(Rule::OrphanClip);
    }
    fired
}

/// Returns `true` iff the alignment set passes the `split` rule.
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
            is_supplementary: false,
            has_sa: false,
            max_skip: 0,
            ins_near_skip: false,
        }
    }

    fn cfg(min_gap: u64, max_gap: Option<u64>, min_clip: u32, side: FlankSide) -> PredicateConfig {
        PredicateConfig {
            min_gap,
            max_gap,
            min_flank_clip: min_clip,
            flank_side: side,
            rules: RuleSet::from_rules(&[Rule::Split]),
            long_intron_min_skip: 150_000,
            long_intron_min_clip: 21,
            orphan_min_clip: 50,
        }
    }

    fn all_rules() -> PredicateConfig {
        PredicateConfig {
            rules: RuleSet::from_rules(&Rule::ALL),
            ..cfg(150_000, Some(1_600_000), 20, FlankSide::Inner)
        }
    }

    #[test]
    fn split_only_config_ignores_single_record_rules() {
        let mut a = aln(0, 100, 400_000, false, 80, 0);
        a.ins_near_skip = true;
        a.max_skip = 300_000;
        let c = cfg(150_000, Some(1_600_000), 20, FlankSide::Inner);
        assert!(candidate_rules(&[a], &c).is_empty());
    }

    #[test]
    fn insertion_rule_fires_on_primary() {
        let mut a = aln(0, 100, 400_000, false, 0, 0);
        a.ins_near_skip = true;
        a.max_skip = 264_397;
        let fired = candidate_rules(&[a], &all_rules());
        assert_eq!(fired.labels(), "insertion");
    }

    #[test]
    fn long_intron_clip_needs_both_clip_and_skip() {
        let mut a = aln(0, 100, 400_000, false, 21, 0);
        a.max_skip = 150_000;
        a.has_sa = true;
        assert_eq!(
            candidate_rules(&[a], &all_rules()).labels(),
            "long-intron-clip"
        );

        a.left_clip = 20;
        assert!(candidate_rules(&[a], &all_rules()).is_empty());

        a.left_clip = 21;
        a.max_skip = 149_999;
        assert!(candidate_rules(&[a], &all_rules()).is_empty());
    }

    #[test]
    fn orphan_clip_needs_no_sa() {
        let mut a = aln(0, 100, 5_000, true, 0, 50);
        assert_eq!(candidate_rules(&[a], &all_rules()).labels(), "orphan-clip");
        a.has_sa = true;
        assert!(candidate_rules(&[a], &all_rules()).is_empty());
        a.has_sa = false;
        a.right_clip = 49;
        assert!(candidate_rules(&[a], &all_rules()).is_empty());
    }

    #[test]
    fn single_record_rules_read_the_primary_not_the_supplementary() {
        let primary = aln(0, 100, 5_000, false, 0, 0);
        let mut supp = aln(0, 900_000, 901_000, false, 70, 0);
        supp.is_supplementary = true;
        supp.ins_near_skip = true;
        supp.max_skip = 200_000;
        // The pair is a split (895 kb gap, 70-nt gap-facing clip); only the
        // split rule may fire.
        assert_eq!(candidate_rules(&[supp, primary], &all_rules()).labels(), "split");
    }

    #[test]
    fn rules_combine_and_split_is_reported() {
        let mut p = aln(0, 100, 200, false, 0, 50);
        p.has_sa = true;
        let mut s = aln(0, 600_000, 600_100, false, 50, 0);
        s.is_supplementary = true;
        s.has_sa = true;
        assert_eq!(candidate_rules(&[p, s], &all_rules()).labels(), "split");
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
