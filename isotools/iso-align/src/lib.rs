// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! iso-align: select long reads whose pass-1 split alignment suggests a
//! re-align with a larger minimap2 intron cap, and emit them as
//! FASTA / names.
//!
//! The tool runs in three stages:
//!
//! 1. **BAM scan** ([`bam_scan`]) — single-pass over the pass-1 BAM, builds a
//!    per-QNAME map of compact [`types::MiniAln`] records (no SEQ / QUAL).
//! 2. **Predicate** ([`predicate`]) — flags QNAMEs whose alignment set is
//!    a same-chrom / same-strand split with a sufficiently large gap and
//!    a flank clip on the gap-facing side.
//! 3. **Extract** ([`extract`]) — re-scans the BAM, writing matching records
//!    out as FASTA (default), or writes just a name list.

pub mod bam_scan;
pub mod cli;
pub mod engine;
pub mod error;
pub mod extract;
pub mod logging;
pub mod predicate;
pub mod types;

use anyhow::Result;
use clap::Parser;

/// Runs iso-align from CLI arguments.
pub fn run() -> Result<()> {
    let cli = cli::Cli::parse();
    logging::init(cli.level.into())?;
    engine::run(cli)
}
