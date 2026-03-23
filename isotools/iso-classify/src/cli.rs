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
#[command(version, about, long_about = None)]
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
        args: IntronArgs,
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
        requires("twobit")
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
        requires("twobit"),
        requires("toga")
    )]
    pub nag: bool,

    #[arg(
        short = 'I',
        long = "iic",
        required = false,
        value_name = "FLAG",
        help = "Path to intron intronIC .fa for intron type classification [U2, U12]"
    )]
    pub iic: Option<PathBuf>,

    #[arg(
        long = "rt_freq_threshold",
        required = false,
        value_name = "FLOAT",
        num_args = 1,
        help = "RT intron frequency threshold",
        default_value_t = 0.5,
        default_missing_value("0.5")
    )]
    pub rt_freq_threshold: f64,

    #[arg(
        long = "outdir",
        required = false,
        value_name = "PATH",
        num_args = 1,
        help = "Path to output directory",
        default_value = env!("CARGO_MANIFEST_DIR"),
    )]
    pub outdir: PathBuf,
}

#[derive(Debug, Parser)]
pub struct ExonArgs {}
