use clap::{self, Parser};
use config::ArgCheck;
use log::{error, info, Level};
use simple_logger::init_with_level;

use iso_classify::cli::{Args, SubArgs};

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
        SubArgs::Intron { args } => {
            use iso_classify::core::classify_introns;

            args.check().unwrap_or_else(|e| {
                error!("{}", e);
                std::process::exit(1);
            });

            classify_introns(args).unwrap_or_else(|e| {
                error!("{}", e);
                std::process::exit(1);
            });
        }
        SubArgs::Exon { args } => {
            todo!()
        }
    }

    let elapsed = start.elapsed();
    info!("Elapsed time: {:.3?}", elapsed);
}
