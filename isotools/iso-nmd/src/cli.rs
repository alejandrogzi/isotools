//! Core module for detecting non-mediated decays in a query set of reads
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main function for detecting non-mediated decays (NMDs)
//! and processing the components of reads and introns in parallel.
//!
//! In short, identifies and categorizes transcripts based on nonsense-mediated decay (NMD)
//! rules for each read. It processes reads in parallel, filtering out blacklisted entries.
//! For each transcript, it calculates key metrics like the length of the coding sequence
//! and the 3' UTR. Using these metrics and predefined thresholds, it assigns a tag: strong NMD,
//! weak NMD, or no NMD. The final output separates NMD-free reads from NMD reads,
//! with the latter group tagged and color-coded for easy visualization.

use clap::{ArgAction, Parser};
use std::path::PathBuf;

pub const NMD_DISTANCE: u64 = 55; // 55 bp
pub const WEAK_NMD_DISTANCE: i64 = 80; // 85 bp
pub const ATG_DISTANCE: u64 = 100; // 100 bp
pub const BIG_EXON_DIST_TO_EJ: u64 = 400; // 400 bp

#[derive(Debug, Parser)]
pub struct Args {
    #[arg(
        short = 'r',
        long = "ref",
        required = true,

        value_name = "PATHS",
        value_delimiter = ',',
        num_args = 1..,
        help = "Path to .bed files"
    )]
    pub refs: Vec<PathBuf>,

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
        short = 'o',
        long = "outdir",
        required = true,
        value_name = "PATH",
        help = "Path to BED4 file with blacklisted introns"
    )]
    pub outdir: PathBuf,

    #[arg(
        short = 't',
        long = "threads",
        help = "Number of threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get(),
    )]
    pub threads: usize,

    #[arg(
        short = 'n',
        long = "nmd-distance",
        help = "Distance to consider NMD",
        value_name = "DISTANCE",
        default_value_t = NMD_DISTANCE,
        action = ArgAction::Set
    )]
    pub nmd_distance: u64,

    #[arg(
        short = 'w',
        long = "weak-nmd-distance",
        help = "Distance to consider weak NMD",
        value_name = "DISTANCE",
        default_value_t = WEAK_NMD_DISTANCE,
        action = ArgAction::Set
    )]
    pub weak_nmd_distance: i64,

    #[arg(
        short = 'a',
        long = "atg-distance",
        help = "Distance to consider ATG",
        value_name = "DISTANCE",
        default_value_t = ATG_DISTANCE,
        action = ArgAction::Set
    )]
    pub atg_distance: u64,

    #[arg(
        short = 'e',
        long = "big-exon-dist-to-ej",
        help = "Distance to consider big exon to exon junction",
        value_name = "DISTANCE",
        default_value_t = BIG_EXON_DIST_TO_EJ,
        action = ArgAction::Set
    )]
    pub big_exon_dist_to_ej: u64,

    #[arg(
        short = 'p',
        long = "prefix",
        required = false,
        value_name = "VALUE",
        help = "Prefix to append at the beginning of the read name"
    )]
    pub prefix: Option<String>,
}
