use clap::{self, Parser};
use config::ArgCheck;
use log::{error, info, Level};
use simple_logger::init_with_level;

use iso_fusion::{
    cli::Args,
    core::{detect_fusions, detect_fusions_with_mapping},
};

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

    if !args.map {
        info!("Detecting fusions in default mode...");

        detect_fusions(args).unwrap_or_else(|e| {
            error!("{}", e);
            std::process::exit(1);
        });
    } else {
        info!("Detecting fusions in isoform mapping mode...");

        detect_fusions_with_mapping(args).unwrap_or_else(|e| {
            error!("{}", e);
            std::process::exit(1);
        });
    }

    let elapsed = start.elapsed();
    info!("Elapsed time: {:?}", elapsed);
}
