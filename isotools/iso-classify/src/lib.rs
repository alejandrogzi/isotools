use anyhow::Result;
use core::classify_introns;
use std::path::PathBuf;

pub mod cli;
pub mod core;
pub mod utils;

pub fn lib_iso_classify(args: Vec<String>) -> Result<PathBuf> {
    let args = cli::IntronArgs::from(args);
    let introns = classify_introns(args);

    return introns;
}
