//! Core module for detecting real 3' ends on isoseq reads
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main function for detecting real 3' ends
//! on isoseq reads and processing the components of reads in parallel.
//!
//! In short, uses APARENT and HMM polyA scoring systems to define real ends.  
//! Finally, the caller module groups all the previous information and tries to
//! determine the intraprimming potential for each read,

use clap::{self, Parser};
use log::{error, info};
use simple_logger::init_with_level;

use iso_pas::{cli::Args, core::pas_caller};

fn main() {
    let start = std::time::Instant::now();

    let args: Args = Args::parse();
    init_with_level(args.level).unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .unwrap();

    pas_caller(args).unwrap_or_else(|e| {
        error!("{}", e);
        std::process::exit(1);
    });

    let elapsed = start.elapsed();
    info!("Elapsed time: {:.3?}", elapsed);
}
