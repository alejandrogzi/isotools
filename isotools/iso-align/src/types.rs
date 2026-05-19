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
}
