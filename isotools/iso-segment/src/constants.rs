// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Core module for segmenting reads based on their polyA features
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main functions for grouping reads
//! and processing components based on polyA features in parallel.
//!
//! In short, this modules provides iso-segment. The segment module
//! filters reads based on alignment quality and predicts the polyA
//! tail using a two-state HMM model.

/// Minimum percentage of identity required for valid alignment
pub const MIN_PER_ID: usize = 98;

/// Maximum allowed 5' soft or hard clip length
pub const MAX_CLIP5: usize = 20;

/// Maximum allowed 3' soft or hard clip length
pub const MAX_CLIP3: usize = 20;

/// Number of bases at the end of the read to consider as tail suffix
pub const POLYA_SUFFIX: usize = 30;

/// Step size in bp to move backwards when searching for tail
pub const SUFFIX_STEP_SIZE: usize = 50;

/// Minimum identity threshold for accepting reads (percent)
pub const IDENTITY_THRESHOLD: f32 = 98.0;

/// Minimum identity to discard extremely low quality reads (percent)
pub const MINIMUM_IDENTITY: f32 = 60.0;

/// Transition probability of staying in polyA state
// [HMM parameters]
// INFO: transition prob for polyA tail
pub const P2P: f64 = 0.9;
/// Emission probability of 'A' in polyA state
pub const EMIT_A: f64 = 0.99;

// file name prefixes
/// Prefix for high-quality (accepted) output files
pub const PASS_PREFIX: &str = "hq";
/// Prefix for low-quality (rejected) output files
pub const FAIL_PREFIX: &str = "lq";

// color schema
/// RGB color for accepted reads in BED format
pub const RGB_ACCEPT: &str = "43,118,219";
/// RGB color for rejected reads in BED format
pub const RGB_REJECT: &str = "213,67,67";

// config derived constants
/// Separator between read number and chromosome in tagged read names
pub const BIG_SEP: &str = "__";
/// Separator between fields in tagged read names
pub const SEP: &str = "#";

// collections
/// Nucleotide complement lookup table
pub const COMPLEMENT: [u8; 128] = {
    let mut nt = [0; 128];
    nt[b'A' as usize] = b'T';
    nt[b'T' as usize] = b'A';
    nt[b'C' as usize] = b'G';
    nt[b'G' as usize] = b'C';
    nt[b'a' as usize] = b't';
    nt[b't' as usize] = b'a';
    nt[b'c' as usize] = b'g';
    nt[b'g' as usize] = b'c';
    nt[b'N' as usize] = b'N';
    nt[b'n' as usize] = b'n';
    nt
};
