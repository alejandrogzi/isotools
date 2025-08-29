//! Core module for splitting a .fa/.fq file into chunks
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main function for splitting .fa/.fq files
//! based on custom requirements in parallel.
//!
//! In short, the module accepts any type of .fa or .fq file
//! and process the reads or sequences inside them in parallel
//! when is possible. Compressed files are also accepted. The
//! user has the ability to specify is the splitting process should
//! be done based on specific chunk sizes or number of files, and
//! the amount of parallelization that should be used in the process.

use clap::Parser;
use config::SplitMode;
use std::path::PathBuf;

#[derive(Debug, Parser)]
pub struct Args {
    #[arg(
        short = 'f',
        long = "file",
        required = true,
        value_name = "PATH",
        help = ".fa/.fq file to split"
    )]
    pub file: PathBuf,

    #[arg(
        short = 'c',
        long = "chunks",
        required = false,
        value_name = "CHUNKS",
        conflicts_with("files"),
        help = "Number of chunks [amount of records in each splitted file]"
    )]
    pub chunks: Option<usize>,

    #[arg(
        short = 'F',
        long = "files",
        required = false,
        value_name = "FILES",
        conflicts_with("chunks"),
        help = "Number of files to split the input in"
    )]
    pub files: Option<usize>,

    #[arg(
        short = 't',
        long = "threads",
        help = "Number of threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get()
    )]
    pub threads: usize,

    #[arg(
        short = 'o',
        long = "outdir",
        required = false,
        value_name = "PATH",
        help = "Output directory path",
        default_value("chunks")
    )]
    pub outdir: PathBuf,

    #[arg(
        short = 's',
        long = "suffix",
        required = false,
        value_name = "VALUE",
        help = "Suffix to append at the end of the chunk file"
    )]
    pub suffix: Option<String>,
}

impl Args {
    pub fn from(args: Vec<String>) -> Self {
        let mut full_args = vec![env!("CARGO_PKG_NAME").to_string()];
        full_args.extend(args);

        Args::parse_from(full_args)
    }

    pub fn mode(&self) -> anyhow::Result<SplitMode> {
        match (self.chunks, self.files) {
            (Some(n), None) => Ok(SplitMode::ChunkSize(n)),
            (None, Some(n)) => Ok(SplitMode::NumFiles(n)),
            _ => Err(anyhow::anyhow!(
                "You must provide either --chunks or --files"
            )),
        }
    }
}
