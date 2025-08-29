//! Core module for splitting a .fa/.fq file into chunks
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main function for splitting .fa/.fq files
//! based on custom requirements in parallel.
//!
//! In short, the module accepts any type of .fa or .fq file
//! and process the reads or sequences inside them in parallel
//! when is possible. Compressed files are also accepted. The
//! user has the ability to specify is the splitting process should
//! be done based on specific chunk sizes or number of files, and
//! the amount of parallelization that should be used in the process.

use anyhow::{Ok, Result};
use clap::{self, Parser};
use iso_split::cli::Args;
use log::{info, Level};
use simple_logger::init_with_level;

use iso_split::*;

fn main() -> Result<()> {
    let start = std::time::Instant::now();
    init_with_level(Level::Info).unwrap();

    let args: Args = Args::parse();

    let _ = dispatch!(&args.file, {
        "fa.gz" => split_fa_gz(&args)?,
        "fasta.gz" =>  split_fa_gz(&args)?,
        "fq.gz" =>  split_fq(&args)?,
        "fastq.gz" =>  split_fq(&args)?,
        "fa" =>  split_fa(&args)?,
        "fasta" =>  split_fa(&args)?,
    });

    let elapsed = start.elapsed();
    info!("Elapsed time: {:?}", elapsed);

    Ok(())
}
