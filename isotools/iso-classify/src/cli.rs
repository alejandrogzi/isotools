// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Core module for intron classification
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module provides a comprehensive sub-pipeline for classifying
//! introns within genomic sequences. It different data sources to
//! categorize introns based on their predicted splicing potential
//! and structural characteristics.
//!
//! In essence, this module identifies and characterizes introns from
//! input long-read sequencing data. It performs data integration,
//! collecting splice site prediction scores (from tools like SpliceAI
//! and MaxEntScan), analyzing genomic sequence context, and detecting
//! specific sequence patterns such as RT repeats and NAG motifs. Through
//! a parallel processing approach, each intron is evaluated to determine
//! its "support type", indicating whether it is likely to be a genuine
//! spliced intron, an RT-driven event, or an unclear case requiring
//! further investigation. The final output is a detailed, classified list
//! of introns, enabling deeper insights into alternative splicing and
//! RNA processing.

use clap::{ArgAction, Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[clap(
    name = "iso-classify",
    version = env!("CARGO_PKG_VERSION"),
    author = "Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>",
    about = "Classify introns based on their splice site prediction and structural characteristics",
)]
pub struct Args {
    #[command(subcommand)]
    pub command: SubArgs,

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

impl Args {}

#[derive(Debug, Subcommand)]
pub enum SubArgs {
    #[command(name = "intron")]
    Intron {
        #[command(flatten)]
        args: Box<IntronArgs>,
    },
    #[command(name = "exon")]
    Exon {
        #[command(flatten)]
        args: ExonArgs,
    },
}

#[derive(Debug, Parser)]
pub struct IntronArgs {
    #[arg(
        short = 'i',
        long = "isoseq",
        required = true,
        value_name = "PATH",
        help = "Paths to IsoSeq BED12 file"
    )]
    pub isoseq: PathBuf,

    #[arg(
        short = 'w',
        long = "bigwig",
        required = false,
        value_name = "PATH",
        num_args = 1,
        help = "Path to spliceAI directory [will asume 2 files per strand: acceptor and donor .bw]"
    )]
    pub spliceai: Option<PathBuf>,

    #[arg(
        short = 's',
        long = "sequence",
        required = false,
        value_name = "PATH",
        help = "Path to genome 2bit/fa/fa.gz file"
    )]
    pub sequence: Option<PathBuf>,

    #[arg(
        short = 't',
        long = "toga",
        required = false,
        value_name = "PATH",
        value_delimiter = ',',
        num_args = 1..,
        help = "Path to TOGA annotation .bed file"
    )]
    pub toga: Option<Vec<PathBuf>>,

    #[arg(
        long = "scan",
        required = false,
        value_name = "FLAG",
        help = "Use MaxEntScan for splice site prediction",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
        requires("sequence")
    )]
    pub scan: bool,

    #[arg(
        long = "nag",
        required = false,
        value_name = "FLAG",
        help = "Use TOGA-nag for splice site prediction",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
        requires("sequence"),
        requires("toga")
    )]
    pub nag: bool,

    #[arg(
        short = 'I',
        long = "iic",
        required = false,
        value_name = "PATH",
        help = "Path to intron intronIC .fa for intron type classification [U2, U12]"
    )]
    pub iic: Option<PathBuf>,

    #[arg(
        short = 'R',
        long = "repeats",
        required = false,
        value_name = "PATH",
        help = "Path to BED3 genomic repeats file"
    )]
    pub repeats: Option<PathBuf>,

    #[arg(
        short = 'F',
        long = "rt_freq_threshold",
        required = false,
        value_name = "FLOAT",
        num_args = 1,
        help = "RT intron frequency threshold",
        default_value_t = 0.5
    )]
    pub rt_freq_threshold: f32,

    #[arg(
        short = 'S',
        long = "spliceai_min_ss_signal",
        required = false,
        value_name = "FLOAT",
        num_args = 1,
        help = "SpliceAI minimum splice signal threshold",
        default_value_t = 0.02
    )]
    pub spliceai_min_ss_signal: f32,

    #[arg(
        short = 'f',
        long = "intron_freq_threshold",
        required = false,
        value_name = "FLOAT",
        num_args = 1,
        help = "Intron frequency threshold",
        default_value_t = 0.5
    )]
    pub intron_freq_threshold: f32,

    #[arg(
        short = 'M',
        long = "maxent_min_ss_signal",
        required = false,
        value_name = "FLOAT",
        num_args = 1,
        help = "MaxEnt minimum splice signal threshold",
        default_value_t = 1.5
    )]
    pub maxent_min_ss_signal: f32,

    #[arg(
        short = 'p',
        long = "prefix",
        required = false,
        value_name = "PATH",
        num_args = 1,
        help = "Prefix for output files"
    )]
    pub prefix: Option<String>,

    #[arg(
        long = "intron-track",
        required = false,
        value_name = "FLAG",
        help = "Flag to output artifacts/unclear/RT introns as a BED4 track",
        action = ArgAction::SetTrue,
    )]
    pub do_intron_track: bool,

    #[arg(
        long = "outdir",
        required = false,
        value_name = "PATH",
        num_args = 1,
        help = "Path to output directory",
        default_value = "."
    )]
    pub outdir: PathBuf,
}

#[derive(Debug, Parser)]
pub struct ExonArgs {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_intron_outdir_defaults_to_current_directory() {
        let args = Args::parse_from(["iso-classify", "intron", "--isoseq", "input.bed"]);

        let SubArgs::Intron { args } = args.command else {
            panic!("ERROR: Failed to parse intron subcommand arguments");
        };

        assert_eq!(args.outdir, PathBuf::from("."));
    }

    #[test]
    fn test_scan_requires_sequence_argument() {
        let err =
            Args::try_parse_from(["iso-classify", "intron", "--isoseq", "input.bed", "--scan"])
                .expect_err("ERROR: --scan should require --sequence");

        assert!(err.to_string().contains("--sequence"));
    }
}
