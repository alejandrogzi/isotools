// 1. read --raw file.bed --calls file.bed, calls needs to have special parsing (Bed8 -> strip name)
// 2. pack them
// 3. iterate over the pack and a) merge information, b) keep track of non-paired reads
// 4. send output to orf_reads.bed and unpaired_reads.bed

use clap::{self, Parser};
use config::ArgCheck;
use log::{error, info, Level};
use simple_logger::init_with_level;

use iso_orf::{cli::Args, core::call_orfs};

fn main() {
    let start = std::time::Instant::now();
    init_with_level(Level::Info).unwrap();

    let args: Args = Args::parse();
    args.check().unwrap_or_else(|e| {
        error!("{}", e);
        std::process::exit(1);
    });

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .unwrap();

    call_orfs(args).unwrap_or_else(|e| {
        error!("{}", e);
        std::process::exit(1);
    });

    let elapsed = start.elapsed();
    info!("Elapsed time: {:.3?}", elapsed);
}
