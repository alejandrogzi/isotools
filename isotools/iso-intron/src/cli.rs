// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use clap::{ArgAction, Parser, ValueEnum};
use std::path::PathBuf;

pub const RETENTION_RATIO_THRESHOLD: f32 = 0.5;

/// Which introns of the `--introns` table can be retained
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum RetentionSource {
    /// Every intron in the table
    All,
    /// Only introns spliced by at least one record of the --query BED
    Query,
}

#[derive(Debug, Parser)]
#[clap(
    name = "iso-intron",
    version = env!("CARGO_PKG_VERSION"),
    author = "Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>",
    about = "Detects intron retentions in a query set of reads"
)]
pub struct Args {
    #[arg(
        short = 'i',
        long = "introns",
        required = true,
        value_name = "PATH",
        help = "Paths to reference_introns TSV file produced by iso-classify"
    )]
    pub introns: PathBuf,

    #[arg(
        short = 'q',
        long = "query",
        required = true,
        value_name = "PATH",
        help = "Path to BED12 file to classify"
    )]
    pub query: PathBuf,

    #[arg(
        short = 'b',
        long = "blacklist",
        required = false,
        value_name = "PATH",
        help = "Path to BED4 file with blacklisted introns"
    )]
    pub blacklist: Option<PathBuf>,

    #[arg(
        short = 't',
        long = "threads",
        help = "Number of threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get()
    )]
    pub threads: usize,

    #[arg(
        long = "recover",
        help = "Flag to recover from disputed retentions",
        value_name = "FLAG",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
    )]
    pub recover: bool,

    #[arg(
        short = 'r',
        long = "ratio_threshold",
        help = "Ratio threshold for intron retentions",
        value_name = "RATIO",
        default_value_t = RETENTION_RATIO_THRESHOLD
    )]
    pub ratio_threshold: f32,

    #[arg(
        short = 'p',
        long = "prefix",
        required = false,
        value_name = "PATH",
        help = "Prefix for output files",
        default_value = "isotools"
    )]
    pub prefix: String,

    #[arg(
        long = "outdir",
        short = 'o',
        required = false,
        value_name = "PATH",
        num_args = 1,
        help = "Path to output directory",
        default_value = "."
    )]
    pub outdir: PathBuf,

    #[arg(
        short = 'u',
        long = "ignore-utr",
        help = "Flag to ignore intron retentions in UTR regions",
        value_name = "FLAG",
        action = ArgAction::SetTrue
    )]
    pub ignore_utr: bool,

    #[arg(
        long = "allow-missing",
        help = "Allow read/retention introns missing from reference (warn instead of panic)",
        value_name = "FLAG",
        action = ArgAction::SetTrue
    )]
    pub allow_missing: bool,

    #[arg(
        short = 's',
        long = "retention-source",
        help = "Introns a record can retain: 'all' table introns, or only those spliced by at least one --query record (the table can hold introns of reads that never reach the query set)",
        value_name = "SOURCE",
        value_enum,
        default_value_t = RetentionSource::All
    )]
    pub retention_source: RetentionSource,
}
