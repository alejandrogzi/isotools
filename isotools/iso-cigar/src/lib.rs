// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

mod annotation;
mod cli;
mod engine;
mod genome;
mod logging;

use anyhow::Result;
use clap::Parser;

/// Runs iso-cigar
pub fn run() -> Result<()> {
    let cli = cli::Cli::parse();
    logging::init(cli.level.into())?;
    engine::run(cli)
}
