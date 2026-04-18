// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

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
        help = "Overlapping mode, accounting for reads that overlap a reference (requires --ref arg)",
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

    #[arg(
        short = 'M',
        long = "overlap-mode",
        help = "Pack overlapping mode (default: CDS) (possible values: CDS, exon, bounds)",
        value_name = "MODE",
        default_value_t = String::from("CDS")
    )]
    pub overlap_mode: String,

    #[arg(
        short = 'J',
        long = "junction-tolerance",
        help = "Tolerance in bp for splice junction matching (default: 0 = exact match)",
        value_name = "NUM",
        default_value_t = 0
    )]
    pub junction_tolerance: u64,

    #[arg(
        short = 'j',
        long = "min-junction-frac",
        help = "Minimum fraction of query junctions that must match a reference transcript (default: 0.5)",
        value_name = "NUM",
        default_value_t = 0.5
    )]
    pub min_junction_frac: f64,

    #[arg(
        short = 'f',
        long = "min-overlap-frac",
        help = "Minimum reciprocal overlap for single-exon read scoring (default: 0.5)",
        value_name = "NUM",
        default_value_t = 0.5
    )]
    pub min_overlap_frac: f64,

    #[arg(
        short = 'e',
        long = "end-tolerance",
        help = "Tolerance in bp for terminal exon boundary matching (default: 50)",
        value_name = "NUM",
        default_value_t = 0
    )]
    pub end_tolerance: u64,

    #[arg(
        short = 'S',
        long = "splicing-scores",
        help = "Path to directory with 4 SpliceAI BigWig files (donor/acceptor × plus/minus)",
        value_name = "DIR",
        required = false
    )]
    pub splicing_scores: Option<PathBuf>,

    #[arg(
        long = "min-intron-support-frac",
        help = "Minimum fraction of a read's introns that must be individually supported in de novo mode (default: 0.5)",
        value_name = "NUM",
        default_value_t = 0.5
    )]
    pub min_intron_support_frac: f64,

    #[arg(
        long = "intron-support-threshold",
        help = "Minimum fraction of component reads containing an intron for it to count as supported (default: 0.5)",
        value_name = "NUM",
        default_value_t = 0.5
    )]
    pub intron_support_threshold: f64,

    #[arg(
        long = "min-splice-score",
        help = "Minimum median splice-site score from BigWig to keep a de novo read (default: 0.5)",
        value_name = "NUM",
        default_value_t = 0.5
    )]
    pub min_splice_score: f32,
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
