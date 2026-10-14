// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Command-line interface for iso-align.

use std::num::NonZeroUsize;
use std::path::PathBuf;

use clap::{Parser, ValueEnum};
use log::LevelFilter;

use crate::types::{
    AcceptConfig, FlankSide, OutputFormat, PredicateConfig, Rule, RuleSet, ScanConfig,
};

/// Default minimum gap between adjacent segments to flag a candidate (bp).
pub const DEFAULT_MIN_GAP: u64 = 300_000;
/// Default minimum flank-clip length to count toward the clip criterion (bp).
pub const DEFAULT_MIN_FLANK_CLIP: u32 = 20;
/// `insertion` rule: minimum insertion length (bp).
pub const DEFAULT_INS_MIN_LEN: u32 = 20;
/// `insertion` rule: maximum distance from the insertion to its nearest `N` (bp).
pub const DEFAULT_INS_MAX_DIST: u32 = 10;
/// `insertion` rule: minimum length of that `N` (bp).
pub const DEFAULT_INS_MIN_SKIP: u32 = 20_000;
/// `long-intron-clip` rule: minimum longest `N` (bp).
pub const DEFAULT_LONG_INTRON_MIN_SKIP: u32 = 150_000;
/// `long-intron-clip` rule: minimum clip on either end (bp).
pub const DEFAULT_LONG_INTRON_MIN_CLIP: u32 = 21;
/// `orphan-clip` rule: minimum clip on either end (bp).
pub const DEFAULT_ORPHAN_MIN_CLIP: u32 = 50;
/// Accept mode: minimum identity of the realigned record (iso-segment hq).
pub const DEFAULT_MIN_IDENTITY: f32 = 98.0;
/// Accept mode: maximum clip on either end of the realigned record (bp).
pub const DEFAULT_MAX_CLIP: u32 = 20;

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

    /// Original read FASTA/FASTQ file(s), comma-delimited. When set with
    /// FASTA output, sequences are emitted from these files instead of
    /// best-effort reconstruction from BAM.
    #[arg(long = "reads", value_name = "PATH", value_delimiter = ',')]
    pub reads: Option<Vec<PathBuf>>,

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

    /// Selection rules, comma-separated; a read is selected when any fires.
    /// `split` alone is the original predicate. `insertion`,
    /// `long-intron-clip` and `orphan-clip` read the primary record only.
    #[arg(
        long = "rules",
        value_enum,
        value_delimiter = ',',
        default_value = "split"
    )]
    pub rules: Vec<Rule>,

    /// `insertion` rule: minimum insertion length (bp).
    #[arg(long = "ins-min-len", default_value_t = DEFAULT_INS_MIN_LEN)]
    pub ins_min_len: u32,

    /// `insertion` rule: maximum distance from the insertion to its nearest
    /// `N` (bp).
    #[arg(long = "ins-max-dist", default_value_t = DEFAULT_INS_MAX_DIST)]
    pub ins_max_dist: u32,

    /// `insertion` rule: minimum length of the nearest `N` (bp).
    #[arg(long = "ins-min-skip", default_value_t = DEFAULT_INS_MIN_SKIP)]
    pub ins_min_skip: u32,

    /// `long-intron-clip` rule: minimum longest `N` on the primary (bp).
    #[arg(long = "long-intron-min-skip", default_value_t = DEFAULT_LONG_INTRON_MIN_SKIP)]
    pub long_intron_min_skip: u32,

    /// `long-intron-clip` rule: minimum clip on either end of the primary (bp).
    #[arg(long = "long-intron-min-clip", default_value_t = DEFAULT_LONG_INTRON_MIN_CLIP)]
    pub long_intron_min_clip: u32,

    /// `orphan-clip` rule: minimum clip on either end of a primary without an
    /// `SA` tag (bp).
    #[arg(long = "orphan-min-clip", default_value_t = DEFAULT_ORPHAN_MIN_CLIP)]
    pub orphan_min_clip: u32,

    /// Accept mode. BAM of the realigned reads; `--bam` must then be the
    /// pass-1 BAM (it holds the pass-1 `AS`). Writes the names of reads whose
    /// realigned primary has a higher `AS` and passes `--min-identity` /
    /// `--max-clip` to `--output`; `--reads`, `--output-format` and the
    /// selection options are ignored.
    #[arg(long = "realigned-bam", value_name = "PATH")]
    pub realigned_bam: Option<PathBuf>,

    /// Accept mode: minimum identity (%) of the realigned primary, computed
    /// as iso-segment does.
    #[arg(long = "min-identity", default_value_t = DEFAULT_MIN_IDENTITY)]
    pub min_identity: f32,

    /// Accept mode: maximum clip on either end of the realigned primary (bp).
    #[arg(long = "max-clip", default_value_t = DEFAULT_MAX_CLIP)]
    pub max_clip: u32,

    /// Optional sidecar TSV with per-read diagnostics. Selection mode:
    /// `name rname max_gap left_clip right_clip n_alns rules`. Accept mode:
    /// pass-1 / realigned `AS`, identity, clips and the decision. Useful for
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
            rules: RuleSet::from_rules(&self.rules),
            long_intron_min_skip: self.long_intron_min_skip,
            long_intron_min_clip: self.long_intron_min_clip,
            orphan_min_clip: self.orphan_min_clip,
        }
    }

    /// CIGAR thresholds resolved during the BAM scan.
    pub fn scan_config(&self) -> ScanConfig {
        ScanConfig {
            ins_min_len: self.ins_min_len,
            ins_max_dist: self.ins_max_dist,
            ins_min_skip: self.ins_min_skip,
        }
    }

    /// Accept-mode thresholds.
    pub fn accept_config(&self) -> AcceptConfig {
        AcceptConfig {
            min_identity: self.min_identity,
            max_clip: self.max_clip,
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
        assert!(cli.reads.is_none());
        assert_eq!(cli.rules, vec![Rule::Split]);
        assert!(cli.realigned_bam.is_none());
        let p = cli.predicate();
        assert_eq!(p.rules, RuleSet::from_rules(&[Rule::Split]));
    }

    #[test]
    fn rules_parse_comma_delimited() {
        let cli = Cli::try_parse_from([
            "iso-align",
            "--bam",
            "in.bam",
            "--output",
            "out.fa",
            "--rules",
            "split,insertion,long-intron-clip,orphan-clip",
        ])
        .unwrap();
        assert_eq!(cli.predicate().rules, RuleSet::from_rules(&Rule::ALL));
        assert_eq!(
            cli.scan_config(),
            ScanConfig {
                ins_min_len: 20,
                ins_max_dist: 10,
                ins_min_skip: 20_000
            }
        );
    }

    #[test]
    fn unknown_rule_is_rejected() {
        let err = Cli::try_parse_from([
            "iso-align",
            "--bam",
            "in.bam",
            "--output",
            "o",
            "--rules",
            "split,bogus",
        ])
        .unwrap_err();
        assert_eq!(err.kind(), clap::error::ErrorKind::InvalidValue);
    }

    #[test]
    fn accept_mode_args() {
        let cli = Cli::try_parse_from([
            "iso-align",
            "--bam",
            "pass1.bam",
            "--realigned-bam",
            "frag.bam",
            "--output",
            "accepted.txt",
            "--min-identity",
            "97.5",
        ])
        .unwrap();
        assert_eq!(cli.realigned_bam, Some(PathBuf::from("frag.bam")));
        let a = cli.accept_config();
        assert_eq!(a.min_identity, 97.5);
        assert_eq!(a.max_clip, DEFAULT_MAX_CLIP);
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
    fn reads_arg_is_supported() {
        let cli = Cli::try_parse_from([
            "iso-align",
            "--bam",
            "in.bam",
            "--reads",
            "in.fq",
            "--output",
            "out.fa",
        ])
        .unwrap();
        assert_eq!(cli.reads, Some(vec![PathBuf::from("in.fq")]));
    }

    #[test]
    fn reads_arg_accepts_comma_delimited_paths() {
        let cli = Cli::try_parse_from([
            "iso-align",
            "--bam",
            "in.bam",
            "--reads",
            "reads1.fa,reads2.fasta,reads3.fq.gz",
            "--output",
            "out.fa",
        ])
        .unwrap();
        assert_eq!(
            cli.reads,
            Some(vec![
                PathBuf::from("reads1.fa"),
                PathBuf::from("reads2.fasta"),
                PathBuf::from("reads3.fq.gz"),
            ])
        );
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
