// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Command-line interface for iso-align.

use std::num::NonZeroUsize;
use std::path::PathBuf;

use clap::{Parser, ValueEnum};
use log::LevelFilter;

use crate::types::{FlankSide, OutputFormat, PredicateConfig};

/// Default minimum gap between adjacent segments to flag a candidate (bp).
pub const DEFAULT_MIN_GAP: u64 = 300_000;
/// Default minimum flank-clip length to count toward the clip criterion (bp).
pub const DEFAULT_MIN_FLANK_CLIP: u32 = 20;

#[derive(Debug, Clone, Parser)]
#[command(
    author = env!("CARGO_PKG_AUTHORS"),
    version = env!("CARGO_PKG_VERSION"),
    about = "Identify long reads whose pass-1 split alignment suggests a re-align with a larger intron cap, and emit them as FASTA / names"
)]
pub struct Cli {
    /// Pass-1 BAM (alignment selection signal).
    #[arg(short = 'b', long = "bam", value_name = "PATH")]
    pub bam: PathBuf,

    /// Output path. Format is determined by `--output-format`.
    #[arg(short = 'o', long = "output", value_name = "PATH")]
    pub output: PathBuf,

    /// Minimum reference-coord gap between adjacent split segments (bp).
    #[arg(long = "min-gap", default_value_t = DEFAULT_MIN_GAP)]
    pub min_gap: u64,

    /// Maximum reference-coord gap between adjacent split segments (bp).
    /// Defaults to unbounded; supply e.g. 1_600_000 to cap at the next-pass
    /// minimap2 ceiling.
    #[arg(long = "max-gap")]
    pub max_gap: Option<u64>,

    /// Minimum flank-clip length to count toward the clip criterion (bp).
    #[arg(long = "min-flank-clip", default_value_t = DEFAULT_MIN_FLANK_CLIP)]
    pub min_flank_clip: u32,

    /// Whether the clip threshold is checked on the gap-facing flank only
    /// (`inner`, recommended) or anywhere on the read (`any`).
    #[arg(long = "flank-side", value_enum, default_value_t = FlankSide::Inner)]
    pub flank_side: FlankSide,

    /// Output format. `fasta` is the recommended default; `names` emits a
    /// sorted, deduplicated list of read names.
    #[arg(long = "output-format", value_enum, default_value_t = OutputFormat::Fasta)]
    pub output_format: OutputFormat,

    /// Optional sidecar TSV with per-candidate diagnostics
    /// (`name rname max_gap left_clip right_clip n_alns`). Useful for
    /// spot-checking in IGV. Disabled by default.
    #[arg(long = "report", value_name = "PATH")]
    pub report: Option<PathBuf>,

    /// Number of worker threads. Drives the bgzf decompression worker pool
    /// for BAM scanning and BAM-backed FASTA extraction. `1` keeps
    /// decompression on the calling thread.
    #[arg(short = 'T', long = "threads", default_value_t = num_cpus::get())]
    pub threads: usize,

    /// Log level.
    #[arg(short = 'L', long = "level", value_enum, default_value_t = LogLevel::Info)]
    pub level: LogLevel,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
pub enum LogLevel {
    Trace,
    Debug,
    Info,
    Warn,
    Error,
    Off,
}

impl From<LogLevel> for LevelFilter {
    fn from(value: LogLevel) -> Self {
        match value {
            LogLevel::Trace => LevelFilter::Trace,
            LogLevel::Debug => LevelFilter::Debug,
            LogLevel::Info => LevelFilter::Info,
            LogLevel::Warn => LevelFilter::Warn,
            LogLevel::Error => LevelFilter::Error,
            LogLevel::Off => LevelFilter::Off,
        }
    }
}

impl Cli {
    /// Returns the predicate parameters as the testable struct used by
    /// `predicate::is_candidate`.
    pub fn predicate(&self) -> PredicateConfig {
        PredicateConfig {
            min_gap: self.min_gap,
            max_gap: self.max_gap,
            min_flank_clip: self.min_flank_clip,
            flank_side: self.flank_side,
        }
    }

    /// Number of bgzf decompression workers. Always at least 1.
    pub fn bgzf_workers(&self) -> NonZeroUsize {
        NonZeroUsize::new(self.threads.max(1)).expect("max(1) is non-zero")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_sane() {
        let cli =
            Cli::try_parse_from(["iso-align", "--bam", "in.bam", "--output", "out.fa"]).unwrap();
        assert_eq!(cli.min_gap, DEFAULT_MIN_GAP);
        assert!(cli.max_gap.is_none());
        assert_eq!(cli.min_flank_clip, DEFAULT_MIN_FLANK_CLIP);
        assert_eq!(cli.flank_side, FlankSide::Inner);
        assert_eq!(cli.output_format, OutputFormat::Fasta);
    }

    #[test]
    fn explicit_fasta_output_format() {
        let cli = Cli::try_parse_from([
            "iso-align",
            "--bam",
            "in.bam",
            "--output",
            "out.fa",
            "--output-format",
            "fasta",
        ])
        .unwrap();
        assert_eq!(cli.output_format, OutputFormat::Fasta);
    }

    #[test]
    fn names_format_is_supported() {
        let cli = Cli::try_parse_from([
            "iso-align",
            "--bam",
            "in.bam",
            "--output",
            "names.txt",
            "--output-format",
            "names",
        ])
        .unwrap();
        assert_eq!(cli.output_format, OutputFormat::Names);
    }

    #[test]
    fn reads_arg_is_rejected() {
        let err = Cli::try_parse_from([
            "iso-align",
            "--bam",
            "in.bam",
            "--reads",
            "in.fq",
            "--output",
            "out.fa",
        ])
        .unwrap_err();
        assert_eq!(err.kind(), clap::error::ErrorKind::UnknownArgument);
    }

    #[test]
    fn fastq_output_format_is_rejected() {
        let err = Cli::try_parse_from([
            "iso-align",
            "--bam",
            "in.bam",
            "--output",
            "out.fq",
            "--output-format",
            "fastq",
        ])
        .unwrap_err();
        assert_eq!(err.kind(), clap::error::ErrorKind::InvalidValue);
    }

    #[test]
    fn predicate_propagates_max_gap() {
        let cli = Cli::try_parse_from([
            "iso-align",
            "--bam",
            "in.bam",
            "--output",
            "out.fa",
            "--max-gap",
            "1600000",
            "--flank-side",
            "any",
        ])
        .unwrap();
        let p = cli.predicate();
        assert_eq!(p.max_gap, Some(1_600_000));
        assert_eq!(p.flank_side, FlankSide::Any);
    }
}
