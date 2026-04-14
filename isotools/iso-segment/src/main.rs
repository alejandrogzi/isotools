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

use clap::{self, Parser};
use log::{error, info, Level};
use simple_logger::init_with_level;

use iso_segment::*;

fn main() {
    let start = std::time::Instant::now();
    init_with_level(Level::Info).unwrap();

    let args: Args = Args::parse();

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .unwrap();

    pool.install(|| {
        segment(args).unwrap_or_else(|e| {
            error!("{}", e);
            std::process::exit(1);
        });
    });

    let elapsed = start.elapsed();
    info!("Elapsed time: {:.3?}", elapsed);
}
