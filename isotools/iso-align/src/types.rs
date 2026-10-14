// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Shared, lightweight types used across the iso-align stages.

use clap::ValueEnum;

/// Compact projection of a BAM record retained during the candidate scan.
///
/// We keep only the fields required by the predicate so that ~5 M reads with
/// up to a handful of split alignments each fit in a few GB of RAM. Notably
/// QUAL / MD / aux are never touched during the predicate scan; FASTA output
/// comes from a second BAM pass over candidate QNAMEs.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct MiniAln {
    /// Reference sequence index (BAM tid), not the textual chromosome.
    pub rname: u32,
    /// 1-based inclusive start on the reference.
    pub ref_start: i64,
    /// 1-based exclusive end on the reference (start + ref-consuming op lengths).
    pub ref_end: i64,
    /// True iff the alignment is on the reverse strand (flag 0x10).
    pub is_reverse: bool,
    /// Sum of leading S+H lengths (CIGAR is in reference orientation).
    pub left_clip: u32,
    /// Sum of trailing S+H lengths.
    pub right_clip: u32,
    /// True iff the record is supplementary (flag 0x800).
    pub is_supplementary: bool,
    /// True iff the record carries an `SA` tag.
    pub has_sa: bool,
    /// Longest `N` operation in the CIGAR (bp).
    pub max_skip: u32,
    /// True iff an insertion of at least `ins_min_len` lies within
    /// `ins_max_dist` of its nearest `N`, and that `N` is at least
    /// `ins_min_skip` long (see [`ScanConfig`]).
    pub ins_near_skip: bool,
}

/// CIGAR thresholds that have to be resolved while the record is in hand,
/// because [`MiniAln`] keeps a summary of the CIGAR, not the CIGAR itself.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ScanConfig {
    pub ins_min_len: u32,
    pub ins_max_dist: u32,
    pub ins_min_skip: u32,
}

impl ScanConfig {
    /// Thresholds that never flag an insertion.
    pub const fn disabled() -> Self {
        ScanConfig {
            ins_min_len: u32::MAX,
            ins_max_dist: 0,
            ins_min_skip: u32::MAX,
        }
    }
}

/// One selection rule. A read is a candidate when any enabled rule fires.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash, ValueEnum)]
pub enum Rule {
    /// Split alignment across a gap in `[min_gap, max_gap]` with a gap-facing
    /// clip (the original iso-align predicate, on all records of the read).
    Split,
    /// Primary record: an insertion next to a long `N`, the shape of an exon
    /// left unaligned between two long introns (bridge pseudo-exon).
    Insertion,
    /// Primary record: a clip on either end plus one very long `N` (an exon
    /// beyond the intron cap was clipped off instead of split off).
    LongIntronClip,
    /// Primary record: a long clip and no supplementary alignment (`SA`),
    /// i.e. a terminal exon minimap2 could not place at all.
    OrphanClip,
}

impl Rule {
    /// Bit used in a [`RuleSet`].
    pub const fn bit(self) -> u8 {
        match self {
            Rule::Split => 1,
            Rule::Insertion => 2,
            Rule::LongIntronClip => 4,
            Rule::OrphanClip => 8,
        }
    }

    /// Label written to the report.
    pub const fn label(self) -> &'static str {
        match self {
            Rule::Split => "split",
            Rule::Insertion => "insertion",
            Rule::LongIntronClip => "long-intron-clip",
            Rule::OrphanClip => "orphan-clip",
        }
    }

    pub const ALL: [Rule; 4] = [
        Rule::Split,
        Rule::Insertion,
        Rule::LongIntronClip,
        Rule::OrphanClip,
    ];
}

/// Bit set of [`Rule`]s.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash)]
pub struct RuleSet(pub u8);

impl RuleSet {
    pub fn from_rules(rules: &[Rule]) -> Self {
        RuleSet(rules.iter().fold(0, |acc, r| acc | r.bit()))
    }

    pub fn contains(self, rule: Rule) -> bool {
        self.0 & rule.bit() != 0
    }

    pub fn insert(&mut self, rule: Rule) {
        self.0 |= rule.bit();
    }

    pub fn is_empty(self) -> bool {
        self.0 == 0
    }

    /// Comma-joined labels, in [`Rule::ALL`] order.
    pub fn labels(self) -> String {
        Rule::ALL
            .iter()
            .filter(|r| self.contains(**r))
            .map(|r| r.label())
            .collect::<Vec<_>>()
            .join(",")
    }
}

/// Where the flank-clip threshold is checked, relative to the inter-segment gap.
#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
pub enum FlankSide {
    /// Only clips on the side of each segment that *faces the gap* count: the
    /// upstream segment's right-flank clip or the downstream segment's
    /// left-flank clip. This is the semantically correct setting for
    /// "alignment was split because the intron-cap was hit".
    Inner,
    /// Any flank clip on any segment is allowed to satisfy the threshold —
    /// includes end-of-read clipping unrelated to the split, which is mostly
    /// noise for this use case.
    Any,
}

/// Output mode for stage 3.
#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
pub enum OutputFormat {
    /// FASTA — the recommended default for feeding minimap2 `-x splice:hq`
    /// on the second pass.
    Fasta,
    /// Sorted, deduplicated list of candidate read names, one per line.
    /// Composes with `seqkit grep -f`, `samtools view -N`, etc.
    Names,
}

/// Predicate parameters resolved from the CLI, in a form that's easy to test
/// independently of clap.
#[derive(Clone, Copy, Debug)]
pub struct PredicateConfig {
    pub min_gap: u64,
    pub max_gap: Option<u64>,
    pub min_flank_clip: u32,
    pub flank_side: FlankSide,
    /// Enabled rules; `Split` alone reproduces iso-align <= 0.0.4.
    pub rules: RuleSet,
    /// `LongIntronClip`: minimum longest `N` on the primary record (bp).
    pub long_intron_min_skip: u32,
    /// `LongIntronClip`: minimum clip on either end of the primary record (bp).
    pub long_intron_min_clip: u32,
    /// `OrphanClip`: minimum clip on either end of the primary record (bp).
    pub orphan_min_clip: u32,
}

/// Acceptance thresholds for a realigned record (the iso-segment hq test,
/// without polyA trimming of the 3' clip).
#[derive(Clone, Copy, Debug)]
pub struct AcceptConfig {
    pub min_identity: f32,
    pub max_clip: u32,
}
