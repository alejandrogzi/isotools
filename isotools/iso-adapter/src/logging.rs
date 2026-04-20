// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Logging initialization and elapsed-time helper.

use std::sync::OnceLock;
use std::time::Instant;

use anyhow::Result;
use log::LevelFilter;
use simple_logger::SimpleLogger;

static START: OnceLock<Instant> = OnceLock::new();

/// Initializes the logging system with the specified level.
pub fn init(level: LevelFilter) -> Result<()> {
    START.get_or_init(Instant::now);
    SimpleLogger::new().with_level(level).init()?;
    Ok(())
}

/// Returns elapsed time since initialization as a formatted string.
pub fn elapsed() -> String {
    let elapsed = START
        .get()
        .map(Instant::elapsed)
        .unwrap_or_default()
        .as_secs_f64();
    format!("{elapsed:8.3}s")
}
