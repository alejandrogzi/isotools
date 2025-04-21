//! Core module for detecting intron retentions in a query set of reads
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main function for detecting intron retentions
//! and processing the components of reads and introns in parallel.
//!
//! In short, each read is checked for the presence of intron retentions
//! or RT introns. If a read has a true intron retention, it is discarded.
//! If a read has an RT intron, it is also discarded. The veracity of an
//! 'intron' is determined by 'iso-classify', using machine-learning models,
//! ab initio gene prediction, and other heuristics. The process is heavily
//! parallelized to offer fast performance on large datasets.

use config::ModuleMap;
use dashmap::DashMap;
use std::sync::Arc;

pub mod cli;
pub mod core;
pub mod utils;

pub fn lib_iso_intron(args: Arc<Vec<String>>) -> DashMap<String, Box<dyn ModuleMap>> {
    let args = cli::Args::from(args);
    let descriptor = crate::core::detect_intron_retentions(args)
        .expect("ERROR: Failed to detect intron retentions");

    return descriptor;
}
