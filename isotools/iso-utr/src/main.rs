use clap::{self, Parser};
use config::ArgCheck;
use log::{error, info, Level};
use simple_logger::init_with_level;

use iso_utr::{cli::Args, core::detect_truncations};

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

    detect_truncations(args).unwrap_or_else(|e| {
        error!("{}", e);
        std::process::exit(1);
    });

    let elapsed = start.elapsed();
    info!("Elapsed time: {:?}", elapsed);
}
