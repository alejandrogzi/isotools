//! Core module for detecting intron retentions in a query set of reads
//! Alejandro Gonzales-Irribarren, 2026
//!
//! This module contains the main algorithm for detecting truncations
//! in a query set of reads.
//!
//! In short it takes a set of query reads and a set of reference reads and
//! detects truncations in the query reads. It does this by checking if the
//! query reads overlap with any middle exon from the reference set of reads.
//! A recovery step can be performed by evaluating the support of the middle
//! exon in the reference set of reads.

use clap::{self, Parser};
use log::{error, info};
use simple_logger::init_with_level;

use iso_utr::{cli::Args, core::detect_truncations};

fn main() {
    let start = std::time::Instant::now();

    let args: Args = Args::parse();
    init_with_level(args.level).unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .unwrap();

    detect_truncations(args).unwrap_or_else(|e| {
        error!("{}", e);
        std::process::exit(1);
    });

    let elapsed = start.elapsed();
    info!("Elapsed time: {:?}", elapsed);
}
