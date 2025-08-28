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
use config::ArgCheck;
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
        long = "iso",
        required = true,
        value_name = "PATHS",
        value_delimiter = ',',
        num_args = 1..,
        help = "Paths to IsoSeq's BED12 file(s) delimited by comma"
    )]
    pub iso: Vec<PathBuf>,

    #[arg(
        short = 'b',
        long = "blacklist",
        required = false,
        value_name = "PATH",
        value_delimiter = ',',
        num_args = 1..,
        help = "Path to BED4 file with blacklisted introns"
    )]
    pub blacklist: Vec<PathBuf>,

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
        long = "twobit",
        required = false,
        value_name = "PATH",
        num_args = 1,
        help = "Path to genome 2bit file"
    )]
    pub twobit: Option<PathBuf>,

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
        long = "outdir",
        required = false,
        value_name = "PATH",
        num_args = 1,
        help = "Path to output directory",
        default_value = env!("CARGO_MANIFEST_DIR"),
    )]
    pub outdir: PathBuf,
}

impl ArgCheck for IntronArgs {
    fn get_blacklist(&self) -> &Vec<PathBuf> {
        &self.blacklist
    }

    fn get_ref(&self) -> &Vec<PathBuf> {
        &self.iso
    }

    // filled without validation
    fn get_query(&self) -> &Vec<PathBuf> {
        &self.iso
    }
}

impl IntronArgs {
    pub fn from(args: Vec<String>) -> Self {
        let mut local_args = Vec::new();
        let mut iter = args.iter().peekable();

        while let Some(arg) = iter.next() {
            // INFO: skipping --aparent + value
            if arg == "--aparent" {
                iter.next();
                continue;
            }

            if arg == "--query" {
                local_args.push("--iso".to_string());
            } else {
                local_args.push(arg.clone());
            }
        }

        let mut full_args = vec![env!("CARGO_PKG_NAME").to_string()];
        full_args.extend(local_args);
        full_args.push("--scan".to_string());
        full_args.push("--nag".to_string());

        IntronArgs::parse_from(full_args)
    }
}

#[derive(Debug, Parser)]
pub struct ExonArgs {}
