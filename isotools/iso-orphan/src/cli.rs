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

use clap::{ArgAction, ArgGroup, Parser};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[clap(
    name="iso-orphan",
    author = env!("CARGO_PKG_AUTHORS"),
    version = env!("CARGO_PKG_VERSION"),
    about = "detect oprhan reads and background transcription",
    long_about = None)]
#[command(
    group(
        ArgGroup::new("mode")
        .required(true)
        .multiple(false)
        .args(&["all", "overlapping", "non_overlapping"])
    )
)]
#[command(
    group(
        ArgGroup::new("guided")
        .required(true)
        .multiple(true)
        .args(&["refs", "overlapping"])
    )
)]
pub struct Args {
    #[arg(
        short = 'q',
        long = "query",
        required = true,
        value_name = "PATH",
        help = "Path to BED12 file"
    )]
    pub query: PathBuf,

    #[arg(
        short = 'r',
        long = "ref",
        required = false,
        value_name = "PATHS",
        value_delimiter = ',',
        num_args = 1..,
        help = "Paths to reference BED12 files delimited by comma"
    )]
    pub refs: Option<Vec<PathBuf>>,

    #[arg(
        short = 'A',
        long = "all",
        required = false,
        value_name = "FLAG",
        help = "Complete mode, accounting for overlapping and non-overlapping reads (default: true)",
        action = ArgAction::SetTrue,
    )]
    pub all: bool,

    #[arg(
        short = 'O',
        long = "overlapping",
        required = false,
        value_name = "FLAG",
        help = "Overlapping mode, accounting for reads that overlap a TOGA reference (requires TOGA arg)",
        action = ArgAction::SetTrue,
    )]
    pub overlapping: bool,

    #[arg(
        short = 'N',
        long = "non-overlapping",
        required = false,
        value_name = "FLAG",
        help = "Non-overlapping mode, triggers self-guided algorithm to discard orphans",
        action = ArgAction::SetTrue,
        conflicts_with = "toga"
    )]
    pub non_overlapping: bool,

    #[arg(
        short = 'm',
        long = "min-read-num-denovo",
        help = "Threshold for minimum number of reads to be considered a component when self-guided (default: 5)",
        value_name = "NUM",
        default_value_t = 5
    )]
    pub min_read_num_denovo: usize,

    #[arg(
        short = 'P',
        long = "min-discard-percent",
        help = "Threshold for minimum percentage of discards to apply splice match rescueing",
        value_name = "NUM",
        default_value_t = 0.5
    )]
    pub min_discard_percentage: f32,

    #[arg(
        short = 'o',
        long = "outdir",
        help = "Path to output directory (/orphans)",
        value_name = "PATH",
        required = false,
        default_value = "."
    )]
    pub outdir: PathBuf,

    #[arg(
        short = 'p',
        long = "prefix",
        help = "Name of the output file without any extension",
        value_name = "NAME",
        required = false,
        default_value = "output"
    )]
    pub prefix: String,

    #[arg(
        short = 't',
        long = "threads",
        help = "Number of threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get()
    )]
    pub threads: usize,

    #[arg(
        short = 'L',
        long = "level",
        help = "Logging level",
        value_name = "LEVEL",
        default_value_t = log::Level::Info,
    )]
    pub level: log::Level,
}

pub enum Mode {
    Guided,
    DeNovo,
}

impl Mode {
    pub fn from(args: &Args) -> Self {
        if args.all || args.overlapping {
            Mode::Guided
        } else if args.non_overlapping {
            Mode::DeNovo
        } else {
            unreachable!()
        }
    }
}
