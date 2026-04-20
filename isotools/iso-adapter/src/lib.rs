// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! iso-adapter: detect and optionally remove adapter sequences from
//! soft-clipped regions of long-read BAM alignments.
//!
//! Long-read cDNA / direct-RNA sequencing pipelines (PacBio Iso-Seq, ONT)
//! occasionally leave residual adapter, primer or polyA/polyT tails inside
//! the soft-clipped flanks of minimap2 alignments. This tool scans those
//! flanks, reports which known adapters are present, tracks unknown
//! high-frequency clips as potential novel contamination, and — when
//! requested with `-R / --remove-adapters` — rewrites the BAM so that
//! adapter bases are stripped but genuine overhang sequence is preserved.

pub mod adapters;
pub mod cli;
pub mod detector;
pub mod engine;
pub mod error;
pub mod logging;
pub mod stats;
pub mod writer;

use anyhow::Result;
use clap::Parser;

/// Runs iso-adapter.
pub fn run() -> Result<()> {
    let cli = cli::Cli::parse();
    logging::init(cli.level.into())?;
    engine::run(cli)
}
