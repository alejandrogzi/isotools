/// isotools: tools for long-read sequencing data analysis
///
/// This is the entry point for the isotools CLI.
/// It is responsible for parsing the CLI arguments
/// and executing the appropriate subcommand [iso-tool].
///
/// This wrapper offers 7 different subcommands:
/// - iso-fusion
/// - iso-polya
/// - iso-intron
/// - iso-utr
/// - iso-classify
///
/// Each subcommand/submodule offers different functionalities,
/// such as detecting fusions, polyadenylation sites, introns,
/// UTRs, coverage, classifying introns/exons, and ORFs. There
/// is an internal dependence between some of the subcommands.
/// In addition to the latter, isotools also includes to main
/// hidden submodules: 'iso-pack' and 'config'. The former is
/// a modified version of 'packbed' [see https://github.com/alejandrogzi/packbed]
/// and the latter is a configuration file generator with universal
/// constants for the isotools pipeline.
///
/// To get help on the subcommands, you can run:
///
/// ```shell
/// isotools iso-fusion -- --help
/// ```
///
use clap::{Args, Parser, Subcommand};
use isotools::lib;
use log::{error, info, Level};
use simple_logger::init_with_level;

use std::process::Command;

const ENTRY: &str = env!("CARGO_MANIFEST_DIR");
const RELEASES: &str = "target/release";

const HELP: &str = r#"
Usage: isotools run --query <PATHS> --toga <PATH> --aparent <PATH> --bigwig <DIR_PATH> --twobit <PATH>

 Options:
  --query <PATHS>...          Path to Iso-Seqs BED12 file to classify
  --toga <PATH>               Path to TOGA BED12 file to classify
  --aparent <PATH>            Path to APARENT BED12 file to classify [iso-polya aparent output]
  --bigwig <DIR_PATH>         Path to BigWig file to classify
  --twobit <PATH>             Path to TwoBit file to classify
  --outdir <DIR>              Output directory for the classification results
  -h, --help                  Print help
"#;

#[derive(Parser)]
#[command(name = "isotools")]
#[command(about = "isotools: tools for long-read sequencing data analysis")]
#[command(version = env!("CARGO_PKG_VERSION"))]
#[command(author = "Alejandro Gonzales-Irribarren, 2025")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    #[command(name = "iso-fusion")]
    Fusion(IsoArgs),
    #[command(name = "iso-polya")]
    Polya(IsoArgs),
    #[command(name = "iso-intron")]
    Intron(IsoArgs),
    #[command(name = "iso-utr")]
    Utr(IsoArgs),
    #[command(name = "iso-classify")]
    Classify(IsoArgs),
    #[command(name = "run")]
    Run(IsoArgs),
}

#[derive(Args)]
struct IsoArgs {
    #[arg(trailing_var_arg = true, allow_hyphen_values = true, help = HELP)]
    args: Vec<String>,
}

fn main() {
    init_with_level(Level::Info).unwrap();
    let cli = Cli::parse();

    init();

    let (cmd, args) = match cli.command {
        Commands::Fusion(args) => ("iso-fusion", args.args),
        Commands::Polya(args) => ("iso-polya", args.args),
        Commands::Intron(args) => ("iso-intron", args.args),
        Commands::Utr(args) => ("iso-utr", args.args),
        Commands::Classify(args) => ("iso-classify", args.args),
        Commands::Run(args) => ("run", args.args),
    };

    match cmd {
        "run" => lib(args),
        _ => {
            let package = std::path::Path::new(ENTRY)
                .parent()
                .expect("ERROR: Could not get parent dir")
                .join(RELEASES)
                .join(cmd);

            if args.contains(&"--help".to_string()) || args.contains(&"-h".to_string()) {
                let output = Command::new(package)
                    .arg("--help")
                    .output()
                    .expect("ERROR: Failed to execute process");

                check_output(output);
            } else {
                let output = Command::new(package)
                    .args(args)
                    .output()
                    .expect("ERROR: Failed to execute process");

                check_output(output);
            }
        }
    }
}

fn check_output(output: std::process::Output) {
    if output.status.success() {
        info!("{}", String::from_utf8_lossy(&output.stdout));
    } else {
        error!("{}", String::from_utf8_lossy(&output.stderr));
        std::process::exit(1);
    }
}

fn init() {
    let message = format!(
        r#"

        isotools: tools for long-read sequencing data analysis

        this is the entry point for the isotools CLI
        and it is responsible for parsing the CLI arguments
        for each iso-tool:

        - iso-fusion
        - iso-polya
        - iso-intron
        - iso-utr
        - iso-classify

        > version: {}
        > author: alejandro gonzales-irribarren, 2025
        > repository: github.com/alejandrogzi/isotools

        for any bug, please open an issue on the repository.

        * to get help on the subcommands, run:
            isotools <SUBCOMMAND> -- --help

        "#,
        env!("CARGO_PKG_VERSION")
    );

    println!("{}", message);
}
