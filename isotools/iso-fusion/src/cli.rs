// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use clap::{ArgAction, Parser};
use packbed::OverlapType;
use std::path::PathBuf;

pub const FUSION_RATIO_THRESHOLD: f32 = 0.5;

#[derive(Debug, Parser)]
#[command(
    author = env!("CARGO_PKG_AUTHORS"),
    version = env!("CARGO_PKG_VERSION"),
    about = env!("CARGO_PKG_DESCRIPTION"),
)]
pub struct Args {
    #[arg(
        short = 'r',
        long = "ref",
        required = true,
        value_name = "PATH",
        help = "Path to BED12/Isoform [--map] file with rule transcripts"
    )]
    pub refs: Vec<PathBuf>,

    #[arg(
        short = 'q',
        long = "query",
        required = true,
        value_name = "PATHS",
        value_delimiter = ',',
        num_args = 1..,
        help = "Path to BED12 file to classify"
    )]
    pub query: Vec<PathBuf>,

    #[arg(
        short = 't',
        long = "threads",
        help = "Number of threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get()
    )]
    pub threads: usize,

    #[arg(
        short = 'T',
        long = "threshold",
        required = false,
        value_name = "VALUE",
        help = "Threshold for fusion component ratio",
        default_value_t = FUSION_RATIO_THRESHOLD
    )]
    pub threshold: f32,

    #[arg(
        long = "recover",
        help = "Flag to recover from disputed fusions",
        value_name = "FLAG",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
    )]
    pub recover: bool,

    #[arg(
        long = "intron-match",
        help = "Flag to intron-specific match instead splicing match",
        value_name = "FLAG",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
    )]
    pub intron_match: bool,

    #[arg(
        short = 'm',
        long = "map",
        help = "Flag to read an isoforms file [tx->gene] instead of BED12",
        value_name = "FLAG",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
    )]
    pub map: bool,

    #[arg(
        short = 'o',
        long = "overlap-type",
        help = "Type of overlap to consider",
        value_name = "OVERLAP TYPE",
        required = false,
        default_value("exon")
    )]
    pub overlap_type: OverlapType,

    #[arg(
        short = 'b',
        long = "blacklist",
        required = false,
        value_name = "PATH",
        value_delimiter = ',',
        num_args = 1..,
        help = "Path to BED12 file with blacklisted reads"
    )]
    pub blacklist: Vec<PathBuf>,

    #[arg(
        short = 'p',
        long = "prefix",
        required = false,
        value_name = "PATH",
        help = "Prefix for output files",
        default_value("fusion_results")
    )]
    pub prefix: PathBuf,

    #[arg(
        long = "tag",
        required = false,
        value_name = "FLAG",
        help = "Flag to tag fake fusions with :FK",
        action = ArgAction::SetTrue,
    )]
    pub tag: bool,

    #[arg(
        short = 'D',
        long = "descriptor",
        help = "Flag to output descriptor as a .tsv file",
        value_name = "FLAG",
        action = ArgAction::SetTrue,
    )]
    pub descriptor: bool,

    #[arg(
        long = "fuzzy-exclusion",
        help = "Mask same-locus reference parents sharing exons/introns/gene-tags as non-fusible",
        value_name = "FLAG",
        action = ArgAction::SetTrue,
    )]
    pub fuzzy_exclusion: bool,
}
