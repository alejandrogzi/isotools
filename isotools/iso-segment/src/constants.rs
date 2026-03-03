//! Core module for segmenting reads based on their polyA features
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main functions for grouping reads
//! and processing components based on polyA features in parallel.
//!
//! In short, this modules provides iso-segment. The segment module
//! filters reads based on alignment quality and predicts the polyA
//! tail using a two-state HMM model.

// [segment parameters]
pub const MIN_PER_ID: usize = 98;
pub const MAX_CLIP5: usize = 20;
pub const MAX_CLIP3: usize = 20;
pub const POLYA_SUFFIX: usize = 30;
pub const SUFFIX_STEP_SIZE: usize = 50;
pub const IDENTITY_THRESHOLD: f32 = 98.0;
pub const MINIMUM_IDENTITY: f32 = 60.0;

// [HMM parameters]
// INFO: transition prob for polyA tail
pub const P2P: f64 = 0.9;
// INFO: emission prob for A in polyA state
pub const EMIT_A: f64 = 0.99;

// file name prefixes
pub const PASS_PREFIX: &str = "hq";
pub const FAIL_PREFIX: &str = "lq";

// color schema
pub const RGB_ACCEPT: &str = "43,118,219";
pub const RGB_REJECT: &str = "213,67,67";

// config derived constants
pub const BIG_SEP: &str = "__";
pub const SEP: &str = "#";

// collections
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
