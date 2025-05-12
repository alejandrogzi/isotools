//! Core module for detecting fusions in a query set of reads
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main function for detecting fusions
//! and processing the components of reads and introns in parallel.
//!
//! In short, each read is checked for the presence of fusions using
//! a set of reference coordinates. If a read has a true fusion, it
//! is discarded. If a read has an RT intron, it is also discarded.
//! The veracity of a 'fusion' is determined by the exact match of at
//! least 1 splice site in the reference. If a read has a fake fusion, it is
//! marked as such. The results are written to a set of files.

use clap::{self, Parser};
use config::ArgCheck;
use log::{error, info, Level};
use simple_logger::init_with_level;

use iso_fusion::{
    cli::Args,
    core::{detect_fusions, detect_fusions_with_mapping},
};

fn main() {
    let start = std::time::Instant::now();
    init_with_level(Level::Info).unwrap();

    let args: Args = Args::parse();
    args.check().unwrap_or_else(|e| {
        error!("{}", e);
        std::process::exit(1);
    });

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .unwrap();

    if !args.map {
        info!("Detecting fusions in default mode...");
        std::fs::create_dir_all(&args.prefix).expect(&format!(
            "ERROR: Failed to create output directory -> {}",
            args.prefix.display()
        ));

        detect_fusions(args).unwrap_or_else(|e| {
            error!("{}", e);
            std::process::exit(1);
        });
    } else {
        info!("Detecting fusions in isoform mapping mode...");

        detect_fusions_with_mapping(args).unwrap_or_else(|e| {
            error!("{}", e);
            std::process::exit(1);
        });
    }

    let elapsed = start.elapsed();
    info!("Elapsed time: {:?}", elapsed);
}
