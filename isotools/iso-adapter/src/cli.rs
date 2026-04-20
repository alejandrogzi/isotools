// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Command-line interface for iso-adapter.

use std::path::PathBuf;

use clap::{ArgAction, Parser, ValueEnum};
use log::LevelFilter;

/// Default minimum clip length considered for fuzzy matching.
pub const DEFAULT_MIN_CLIP_LEN: usize = 10;
/// Default maximum edit distance for fuzzy fallback.
pub const DEFAULT_MAX_EDIT_DIST: u32 = 3;
/// Default frequency threshold used to flag novel clips.
pub const DEFAULT_FREQ_THRESHOLD: u32 = 5;

#[derive(Debug, Clone, Parser)]
#[command(
    author = env!("CARGO_PKG_AUTHORS"),
    version = env!("CARGO_PKG_VERSION"),
    about = "Detect and optionally remove adapter sequences from soft-clipped long-read BAM alignments"
)]
pub struct Cli {
    #[arg(value_name = "INPUT_BAM", help = "Path to aligned BAM file")]
    pub input: PathBuf,

    #[arg(
        short = 'R',
        long = "remove-adapters",
        action = ArgAction::SetTrue,
        help = "Remove detected adapter sequences (excluding homopolymers) and write a trimmed BAM"
    )]
    pub remove_adapters: bool,

    #[arg(
        short = 'P',
        long = "trim-polya",
        action = ArgAction::SetTrue,
        help = "Trim polyA (3' end) / polyT (5' end) homopolymer runs and write a trimmed BAM. \
                3' polyA: keep the run, trim everything outward of it. \
                5' polyT: remove the run and everything outward."
    )]
    pub trim_polya: bool,

    #[arg(
        short = 'o',
        long = "out-bam",
        value_name = "PATH",
        help = "Output BAM path or directory (only used with -R). If omitted, \
                <input_stem>.without_adapters.bam is written next to the input."
    )]
    pub out_bam: Option<PathBuf>,

    #[arg(
        short = 'T',
        long = "threads",
        default_value_t = num_cpus::get(),
        help = "Number of worker threads"
    )]
    pub threads: usize,

    #[arg(
        short = 'L',
        long = "level",
        value_enum,
        default_value_t = LogLevel::Info,
        help = "Log level"
    )]
    pub level: LogLevel,

    #[arg(
        long = "max-edit-dist",
        default_value_t = DEFAULT_MAX_EDIT_DIST,
        help = "Maximum edit distance for fuzzy adapter matching"
    )]
    pub max_edit_dist: u32,

    #[arg(
        long = "min-clip-len",
        default_value_t = DEFAULT_MIN_CLIP_LEN,
        help = "Minimum soft-clip length to consider for fuzzy matching"
    )]
    pub min_clip_len: usize,

    #[arg(
        long = "freq-threshold",
        default_value_t = DEFAULT_FREQ_THRESHOLD,
        help = "Unknown-clip occurrence count at which a clip is flagged as a candidate novel adapter"
    )]
    pub freq_threshold: u32,
}

/// Log level enum for CLI.
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
    /// Returns the resolved output BAM path, honoring the naming rules in INSTRUCTIONS.md.
    pub fn resolved_output_bam(&self) -> PathBuf {
        let stem = self
            .input
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "output".to_string());
        let derived = format!("{stem}.without_adapters.bam");

        match &self.out_bam {
            None => {
                let mut path = self
                    .input
                    .parent()
                    .map(PathBuf::from)
                    .unwrap_or_else(|| PathBuf::from("."));
                path.push(derived);
                path
            }
            Some(explicit) => {
                let is_dir_like = explicit.is_dir()
                    || explicit
                        .to_string_lossy()
                        .ends_with(std::path::MAIN_SEPARATOR);
                if is_dir_like {
                    explicit.join(derived)
                } else {
                    explicit.clone()
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cli_defaults_sane() {
        let cli = Cli::try_parse_from(["iso-adapter", "input.bam"]).unwrap();
        assert!(!cli.remove_adapters);
        assert_eq!(cli.max_edit_dist, DEFAULT_MAX_EDIT_DIST);
        assert_eq!(cli.min_clip_len, DEFAULT_MIN_CLIP_LEN);
        assert_eq!(cli.freq_threshold, DEFAULT_FREQ_THRESHOLD);
        assert_eq!(cli.level, LogLevel::Info);
    }

    #[test]
    fn cli_remove_flag_parses() {
        let cli = Cli::try_parse_from(["iso-adapter", "-R", "input.bam"]).unwrap();
        assert!(cli.remove_adapters);
    }

    #[test]
    fn cli_polya_flag_parses() {
        let cli = Cli::try_parse_from(["iso-adapter", "-P", "input.bam"]).unwrap();
        assert!(cli.trim_polya);
        assert!(!cli.remove_adapters);
    }

    #[test]
    fn cli_both_trim_flags_independent() {
        let cli = Cli::try_parse_from(["iso-adapter", "-R", "-P", "input.bam"]).unwrap();
        assert!(cli.remove_adapters);
        assert!(cli.trim_polya);
    }

    #[test]
    fn cli_default_output_derives_from_stem() {
        let cli = Cli::try_parse_from(["iso-adapter", "/tmp/sample.bam"]).unwrap();
        let out = cli.resolved_output_bam();
        assert_eq!(out, PathBuf::from("/tmp/sample.without_adapters.bam"));
    }

    #[test]
    fn cli_explicit_file_kept_verbatim() {
        let cli =
            Cli::try_parse_from(["iso-adapter", "-o", "/tmp/out.bam", "/tmp/sample.bam"]).unwrap();
        assert_eq!(cli.resolved_output_bam(), PathBuf::from("/tmp/out.bam"));
    }
}
