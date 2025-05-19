//! Core module for scoring, segmenting and clustering reads
//! based on their polyA features
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main functions for grouping reads
//! and processing components based on polyA features in parallel.
//!
//! In short, this modules provides three subtools, namely: aparent,
//! caller and segment. Each one with a specific goal. The first one
//! runs APARENT, a machine-learning model, to score each read's end.
//! The segment module filters reads based on alignment quality and
//! predicts the polyA tail using a two-state HMM model. Finally, the
//! caller module groups all the previous information and tries to
//! determine the intraprimming potential for each read,

pub mod cli;
pub mod core;
pub mod utils;

use config::ModuleMap;
use dashmap::DashMap;
use std::sync::Arc;

pub fn lib_iso_polya(args: Arc<Vec<String>>) -> DashMap<String, Box<dyn ModuleMap>> {
    let args = cli::CallerArgs::from(args);
    let descriptor = crate::core::pas::pas_caller(args).expect("ERROR: Failed to run PAS caller!");

    return descriptor;
}

pub fn lib_iso_segment(args: Vec<String>) {
    let args = cli::SegmentArgs::from(args);
    let _ = crate::core::segment::segment(args).expect("ERROR: Failed to segment reads!");

    log::info!("SUCCESS: iso-polya segment ran succesfully!");
}
