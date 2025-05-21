use anyhow::{Ok, Result};
use clap::{self, Parser};
use iso_split::cli::Args;
use log::{info, Level};
use simple_logger::init_with_level;

use iso_split::*;

fn main() -> Result<()> {
    let start = std::time::Instant::now();
    init_with_level(Level::Info).unwrap();

    let args: Args = Args::parse();

    let _ = dispatch!(&args.file, {
        "fa.gz" => split_fa_gz(&args)?,
        "fasta.gz" =>  split_fa_gz(&args)?,
        "fq.gz" =>  split_fq(&args)?,
        "fastq.gz" =>  split_fq(&args)?,
        "fa" =>  split_fa(&args)?,
        "fasta" =>  split_fa(&args)?,
    });

    let elapsed = start.elapsed();
    info!("Elapsed time: {:?}", elapsed);

    Ok(())
}
