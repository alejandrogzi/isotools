//! Core module for detecting orphans in a query set of reads
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the core functions for detecting orphans in a query set of reads
//! and processing the components in a parallel fashion.
//!
//! In short, each component of reads is subjected to any of the two modes: guided or
//! self-guided. Guided mode is the default mode and is used when the user provides a TOGA file
//! as input. Self-guided mode is used when the user does not provide a TOGA file as input.
//! Both, guided and self-guided, cover an extensive amount of curated oprhan cases under
//! the assumption that they do not represent a valid source of evidence for transcription.
//! The process is heavily parallellized to offer fast performance on large datasets.
//!
//! # Usage
//!
//! The program can be run using the following command:
//!
//! ```bash
//! iso-orphan <ARGUMENTS>
//! ```
//!
//! ## Arguments
//!
//! - `-b`, `--bed`: Paths to BED12 files delimited by comma
//!     - Required: true
//!     - Value name: PATHS
//!     - Value delimiter: comma
//!     - Number of values: 1..
//!     - Help: Paths to BED12 files delimited by comma
//! - `-t`, `--toga`: Paths to TOGA BED12 files delimited by comma
//!     - Required: false
//!     - Value name: PATHS
//!     - Value delimiter: comma
//!     - Number of values: 1..
//!     - Help: Paths to TOGA BED12 files delimited by comma
//! - `-A`, `--all`: Complete mode, accounting for overlapping and non-overlapping reads (default: true)
//!     - Required: false
//!     - Value name: FLAG
//!     - Help: Complete mode, accounting for overlapping and non-overlapping reads (default: true)
//! - `-O`, `--overlapping`: Overlapping mode, accounting for reads that overlap a TOGA reference (requires TOGA arg)
//!     - Required: false
//!     - Value name: FLAG
//!     - Help: Overlapping mode, accounting for reads that overlap a TOGA reference (requires TOGA arg)
//! - `-N`, `--non-overlapping`: Non-overlapping mode, triggers self-guided algorithm to discard orphans
//!     - Required: false
//!     - Value name: FLAG
//!     - Help: Non-overlapping mode, triggers self-guided algorithm to discard orphans
//! - `-K`, `--keep`: Writes orphan reads to a separate file (default: false)
//!     - Required: false
//!     - Value name: FLAG
//!     - Help: Writes orphan reads to a separate file (default: false)
//! - `-m`, `--min-read-num`: Threshold for minimum number of reads to be considered a component when self-guided (default: 5)
//!     - Required: false
//!     - Value name: NUM
//!     - Default value: 5
//!     - Help: Threshold for minimum number of reads to be considered a component when self-guided (default: 5)
//! - `-o`, `--outdir`: Path to output directory (/orphans)
//!     - Required: false
//!     - Default value: .
//!     - Help: Path to output directory (/orphans)
//! - `-n`, `--name`: Name of the output file without any extension
//!     - Required: false
//!     - Default value: output
//!     - Help: Name of the output file without any extension
//! - `-t`, `--threads`: Number of threads
//!     - Required: false
//!     - Default value: num_cpus::get()
//!     - Help: Number of threads
//! - `-L`, `--level`: Logging level
//!     - Required: false
//!     - Default value: info
//!     - Help: Logging level
//!
//! ## Examples
//!
//! ```bash
//! iso-orphan --bed tests/data/mm39.bed --toga tests/data/mm39.toga.bed --outdir tests/data --name test --threads 1
//! ```
//!
//! ```bash
//! iso-orphan --bed tests/data/mm39.bed --all --overlapping --outdir tests/data --name test --threads 1
//! ```
//!
//! ```bash
//!  iso-orphan -b tests/data/mm39.bed -t tests/data/mm39.toga.bed -A -O -N test -T 1
//! ```

use clap::{self, Parser};
use log::{info, Level};
use simple_logger::init_with_level;

use iso_orphan::{cli::Args, core::__detect_orphans};

fn main() {
    let start = std::time::Instant::now();
    init_with_level(Level::Info).unwrap();

    let args: Args = Args::parse();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .unwrap();

    __detect_orphans(args);

    let elapsed = start.elapsed();
    info!("Elapsed time: {:.3?}", elapsed);
}
