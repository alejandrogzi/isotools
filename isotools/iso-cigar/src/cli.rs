// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use std::path::PathBuf;

use clap::{ArgAction, Parser, ValueEnum};
use log::LevelFilter;

#[derive(Debug, Clone, Parser)]
#[command(
    author = env!("CARGO_PKG_AUTHORS"),
    version = env!("CARGO_PKG_VERSION"),
    about = "Rescue missed 3' splice junctions by rewriting transcript CIGARs"
)]
pub struct Cli {
    #[arg(
        short = 'b',
        long,
        value_name = "PATH",
        required = true,
        help = "Path to BAM file"
    )]
    pub bam: PathBuf,

    #[arg(
        long,
        short = 'a',
        value_name = "PATH",
        required = true,
        help = "Path to annotation file"
    )]
    pub annotation: PathBuf,

    #[arg(
        long,
        short = 's',
        value_name = "PATH",
        required = true,
        help = "Path to sequence file"
    )]
    pub sequence: PathBuf,

    #[arg(
        long,
        short = 'o',
        value_name = "PATH",
        required = false,
        default_value = ".",
        help = "Path to output directory"
    )]
    pub outdir: PathBuf,

    #[arg(
        long, 
        short = 'S', 
        action = ArgAction::SetTrue, 
        help = "Split output BAM file by category: extended, aligned"
    )]
    pub split_bam: bool,

    #[arg(
        long, 
        short = 't', 
        default_value_t = num_cpus::get(), 
        help = "Number of threads", 
        required = false
    )]
    pub threads: usize,

    #[arg(
        long, 
        short = 'L', 
        value_enum, 
        default_value_t = LogLevel::Info
    )]
    pub level: LogLevel,

    #[arg(
        long,
        short = 'w',
        default_value_t = 10,
        required = false,
        help = "Wiggle room for splice junction correction in bases"
    )]
    pub wiggle: u32,

    #[arg(
        long = "clip-cutoff",
        short = 'c',
        default_value_t = 5,
        required = false,
        help = "Clipping cutoff in bases"
    )]
    pub clip_cutoff: u32,

    #[arg(
        long,
        short = 'K',
        action = ArgAction::SetTrue,
        help = "Keep additional corrections"
    )]
    pub keep_additional_corrections: bool,
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
    /// Converts LogLevel to log::LevelFilter.
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
