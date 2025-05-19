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

use clap::{self, Parser};
use config::ArgCheck;
use log::{error, info, Level};
use simple_logger::init_with_level;

use iso_polya::{
    cli::{Args, SubArgs},
    core::{
        apa::{calculate_polya, simulate_polya_reads},
        pas::pas_caller,
        segment::segment,
    },
};

#[allow(unused_variables)]
fn main() {
    let start = std::time::Instant::now();
    init_with_level(Level::Info).unwrap();

    let args: Args = Args::parse();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .unwrap();

    match args.command {
        SubArgs::Aparent { args } => {
            if !args.simulate {
                args.check().unwrap_or_else(|e| {
                    error!("{}", e);
                    std::process::exit(1);
                });

                calculate_polya(args).unwrap_or_else(|e| {
                    error!("{}", e);
                    std::process::exit(1);
                });
            } else {
                info!("INFO: Running iso-polya aparent in simulation mode...");
                info!("INFO: args: {:?}", args);

                simulate_polya_reads(args).unwrap_or_else(|e| {
                    error!("{}", e);
                    std::process::exit(1);
                });
            }
        }
        SubArgs::Segment { args } => {
            segment(args).unwrap_or_else(|e| {
                error!("{}", e);
                std::process::exit(1);
            });
        }
        SubArgs::Caller { args } => {
            pas_caller(args).unwrap_or_else(|e| {
                error!("{}", e);
                std::process::exit(1);
            });
        }
    }

    let elapsed = start.elapsed();
    info!("Elapsed time: {:.3?}", elapsed);
}
