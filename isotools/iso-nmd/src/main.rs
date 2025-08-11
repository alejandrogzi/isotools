//! Core module for detecting non-mediated decays in a query set of reads
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main function for detecting non-mediated decays (NMDs)
//! and processing the components of reads and introns in parallel.
//!
//! In short, identifies and categorizes transcripts based on nonsense-mediated decay (NMD)
//! rules for each read. It processes reads in parallel, filtering out blacklisted entries.
//! For each transcript, it calculates key metrics like the length of the coding sequence
//! and the 3' UTR. Using these metrics and predefined thresholds, it assigns a tag: strong NMD,
//! weak NMD, or no NMD. The final output separates NMD-free reads from NMD reads,
//! with the latter group tagged and color-coded for easy visualization.

use clap::Parser;
use iso_nmd::{cli::Args, core::classify_nmd};
use log::{error, info, Level};
use simple_logger::init_with_level;

fn main() {
    let start = std::time::Instant::now();
    init_with_level(Level::Info).unwrap();

    let args: Args = Args::parse();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .unwrap();

    classify_nmd(args).unwrap_or_else(|e| {
        error!("{}", e);
        std::process::exit(1);
    });

    let elapsed = start.elapsed();
    info!("Elapsed time: {:.3?}", elapsed);
}
