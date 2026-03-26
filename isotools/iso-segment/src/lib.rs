//! Core module for segmenting reads based on their polyA features
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main functions for grouping reads
//! and processing components based on polyA features in parallel.
//!
//! In short, this modules provides iso-segment. The segment module
//! filters reads based on alignment quality and predicts the polyA
//! tail using a two-state HMM model.

//! Command-line argument parser
pub mod cli;
/// Constant values for configuration
pub mod constants;
/// Core segmentation logic
pub mod core;

/// Command-line arguments parsed by clap
pub use cli::Args;

/// Re-export of constant values
pub use constants::*;

/// Main segmentation function for processing BAM files
pub use core::segment;
