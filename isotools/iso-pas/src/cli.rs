use clap::{ArgAction, Parser};
use std::path::PathBuf;

pub const POLYA_LENGTH_THRESHOLD: usize = 50;
pub const GENOMIC_POLYA_THRESHOLD: usize = 5; // INFO: 5 A's in genome
pub const APARENT_THRESHOLD: f32 = 0.01;
pub const WIGGLE: usize = 2;
pub const INTRAPRIMING_RATIO_THRESHOLD: f32 = 0.5;

#[derive(Debug, Parser)]
#[clap(
    name = "iso-pas",
    version = env!("CARGO_PKG_VERSION"),
    author = "Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>",
    about = "polyA tail detection and localization in isoseq reads"
)]

pub struct Args {
    #[arg(
        short = 'F',
        long = "aparent-forward",
        required = true,
        value_name = "PATH",
        help = "Path to output APARENT forward bigWig file"
    )]
    pub aparent_forward: PathBuf,

    #[arg(
        short = 'R',
        long = "aparent-reverse",
        required = true,
        value_name = "PATH",
        help = "Path to output APARENT reverse bigWig file"
    )]
    pub aparent_reverse: PathBuf,

    #[arg(
        short = 'q',
        long = "query",
        required = true,
        value_name = "PATH",
        help = "Path to .bed file [note that specific naming convention is assumed]"
    )]
    pub query: PathBuf,

    #[arg(
        short = 'r',
        long = "refs",
        required = true,
        value_name = "PATH",
        help = "Path to output TOGA projections to localize polyA tails in the genome [CDS, UTR, CDS/UTR]"
    )]
    pub refs: Option<PathBuf>,

    #[arg(
        short = 'w',
        long = "wiggle",
        required = false,
        help = "Wiggle room for polyA tail length",
        value_name = "VALUE",
        default_value_t = WIGGLE,
    )]
    pub wiggle: usize,

    #[arg(
        long = "recover",
        help = "Flag to recover from disputed components where discard ratio is bigger than threshold",
        value_name = "FLAG",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
    )]
    pub recover: bool,

    #[arg(
        long = "max-gpa-length",
        required = false,
        help = "Genomic polyA tail length threshold [max length allowed]",
        value_name = "VALUE",
        default_value_t = GENOMIC_POLYA_THRESHOLD,
    )]
    pub max_gpa_length: usize,

    #[arg(
        long = "min-polya-length",
        required = false,
        help = "PolyA tail length threshold [min length allowed]",
        value_name = "VALUE",
        default_value_t = POLYA_LENGTH_THRESHOLD,
    )]
    pub min_polya_length: usize,

    #[arg(
        long = "aparent-threshold",
        required = false,
        help = "APARENT threshold [min score allowed]",
        value_name = "VALUE",
        default_value_t = APARENT_THRESHOLD,
    )]
    pub aparent_threshold: f32,

    #[arg(
        long = "outdir",
        short = 'o',
        required = false,
        value_name = "PATH",
        num_args = 1,
        help = "Path to output directory",
        default_value = "."
    )]
    pub outdir: PathBuf,

    #[arg(
        short = 'p',
        long = "prefix",
        required = false,
        value_name = "PATH",
        help = "Prefix for output files",
        default_value = "isotools"
    )]
    pub prefix: PathBuf,

    #[arg(
        global = true,
        short = 't',
        long = "threads",
        help = "Number of threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get()
    )]
    pub threads: usize,

    #[arg(
        short = 'L',
        long = "level",
        help = "Logging level",
        value_name = "LEVEL",
        default_value_t = log::Level::Info,
    )]
    pub level: log::Level,
}
