//! Core module for intron classification
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module provides a comprehensive sub-pipeline for classifying
//! introns within genomic sequences. It different data sources to
//! categorize introns based on their predicted splicing potential
//! and structural characteristics.
//!
//! In essence, this module identifies and characterizes introns from
//! input long-read sequencing data. It performs data integration,
//! collecting splice site prediction scores (from tools like SpliceAI
//! and MaxEntScan), analyzing genomic sequence context, and detecting
//! specific sequence patterns such as RT repeats and NAG motifs. Through
//! a parallel processing approach, each intron is evaluated to determine
//! its "support type", indicating whether it is likely to be a genuine
//! spliced intron, an RT-driven event, or an unclear case requiring
//! further investigation. The final output is a detailed, classified list
//! of introns, enabling deeper insights into alternative splicing and
//! RNA processing.

pub(crate) mod cli;
pub(crate) mod core;
pub(crate) mod utils;

pub use cli::{Args, SubArgs};
pub use core::classify_introns;
