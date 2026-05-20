// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Core module for segmenting reads based on their polyA features
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main functions for grouping reads
//! and processing components based on polyA features in parallel.
//!
//! In short, this modules provides iso-segment. The segment module
//! filters reads based on alignment quality and predicts the polyA
//! tail using a two-state HMM model.

use clap::{ArgAction, Parser};
use std::path::PathBuf;

use crate::constants::*;

#[derive(Debug, Parser, Clone)]
#[clap(
    name = "iso-segment",
    version = env!("CARGO_PKG_VERSION"),
    author = "Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>",
    about = "segment isoseq polyA tails",
)]
pub struct Args {
    #[arg(
        short = 'b',
        long = "bam",
        required = true,
        value_name = "PATH",
        value_delimiter = ',',
        num_args = 1,
        help = "Path to .bam file"
    )]
    pub bam: PathBuf,

    #[arg(
        short = 'I',
        long = "identity",
        help = "Min %id (computed without the 5' and 3' clip). Must be [0-100] (percent).",
        value_name = "VALUE",
        default_value_t = IDENTITY_THRESHOLD,
    )]
    pub identity: f32,

    #[arg(
        short = 'i',
        long = "min-identity",
        help = "Mininum % identity, used to discard extremely low reads
        (computed without the 5' and 3' clip). Must be [0-100] (percent).",
        value_name = "VALUE",
        default_value_t = MINIMUM_IDENTITY,
    )]
    pub min_identity: f32,

    #[arg(
        short = 's',
        long = "tail-suffix",
        help = "Suffix of the tail. Determines how many bases to consider at the end of the read",
        value_name = "VALUE",
        default_value_t = POLYA_SUFFIX,
    )]
    pub tail_suffix: usize,

    #[arg(
        short = 'S',
        long = "step-size",
        help = "Determines how many bp to move backwards in the read to find the tail",
        value_name = "VALUE",
        default_value_t = SUFFIX_STEP_SIZE,
    )]
    pub suffix_step_size: usize,

    #[arg(
        long = "clip5",
        help = "Max 5' soft or hard clip",
        value_name = "VALUE",
        default_value_t = MAX_CLIP5,
    )]
    pub max_clip_five: usize,

    #[arg(
        long = "clip3",
        help = "Max 3' soft or hard clip",
        value_name = "VALUE",
        default_value_t = MAX_CLIP3,
    )]
    pub max_clip_three: usize,

    #[arg(
        short = 'P',
        long = "p2p",
        help = "Transition probability of looping in the polyA state ",
        value_name = "VALUE",
        default_value_t = P2P,
    )]
    pub p2p: f64,

    #[arg(
        short = 'E',
        long = "emit-a",
        help = "Probability of emitting A in the polyA state",
        value_name = "VALUE",
        default_value_t = EMIT_A,
    )]
    pub emit_a: f64,

    #[arg(
        short = 'T',
        long = "tag",
        help = "Flag to tag read names with polyA information",
        value_name = "FLAG",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
    )]
    pub tag: bool,

    #[arg(
        short = 'B',
        long = "bed",
        help = "Flag to convert bam to bed",
        value_name = "FLAG",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
    )]
    pub bed: bool,

    #[arg(
        short = 'X',
        long = "split",
        help = "Split output files by chromosome",
        value_name = "FLAG",
        action = ArgAction::SetTrue,
    )]
    pub split: bool,

    #[arg(
        short = 'd',
        long = "delimiter",
        help = "Delimiter inserted between chromosome and prefix in split output names (only used with --split)",
        value_name = "DELIMITER",
        default_value = "@"
    )]
    pub delimiter: String,

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
        long = "prefix",
        required = false,
        value_name = "FILE_PREFIX",
        num_args = 1,
        help = "File name prefix",
        default_value = "file"
    )]
    pub prefix: PathBuf,

    #[arg(
        short = 'b',
        long = "batch",
        required = false,
        value_name = "BATCH",
        help = "Index to use as a prefix for each batch when run in parallel [NUMBER]",
        default_value("")
    )]
    pub batch: String,

    #[arg(
        short = 'S',
        long = "singleton",
        help = "Flag to add singleton (SG) tag to read",
        value_name = "FLAG",
        action = ArgAction::SetTrue,
    )]
    pub singleton: bool,

    #[arg(
        global = true,
        short = 't',
        long = "threads",
        help = "Number of threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get()
    )]
    pub threads: usize,

    #[arg(
        short = 'F',
        long = "fragmented",
        help = "Flag to add fragmenteed (FG) tag to read name",
        value_name = "FLAG",
        action = ArgAction::SetTrue,
    )]
    pub fragmented: bool,
}
