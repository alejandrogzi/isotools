use clap::{self, Parser};
use config::ArgCheck;
use log::{error, info, Level};
use simple_logger::init_with_level;

use iso_polya::{
    cli::{Args, SubArgs},
    core::{apa::calculate_polya, filter::filter_minimap, pas::pas_caller},
};

#[allow(unused_variables)]
fn main() {
    let start = std::time::Instant::now();
    init_with_level(Level::Info).unwrap();

    let args: Args = Args::parse();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .unwrap();

    match args.command {
        SubArgs::Aparent { args } => {
            args.check().unwrap_or_else(|e| {
                error!("{}", e);
                std::process::exit(1);
            });

            calculate_polya(args).unwrap_or_else(|e| {
                error!("{}", e);
                std::process::exit(1);
            });
        }
        SubArgs::Filter { args } => {
            args.check().unwrap_or_else(|e| {
                error!("{}", e);
                std::process::exit(1);
            });

            filter_minimap(args).unwrap_or_else(|e| {
                error!("{}", e);
                std::process::exit(1);
            });
        }
        SubArgs::Caller { args } => {
            pas_caller(args).unwrap_or_else(|e| {
                error!("{}", e);
                std::process::exit(1);
            });
        }
    }

    let elapsed = start.elapsed();
    info!("Elapsed time: {:.3?}", elapsed);
}
