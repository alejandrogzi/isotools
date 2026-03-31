use clap::{ArgAction, Parser};
use std::path::PathBuf;

pub const RETENTION_RATIO_THRESHOLD: f32 = 0.5;

#[derive(Debug, Parser)]
pub struct Args {
    #[arg(
        short = 'r',
        long = "ref",
        required = true,
        value_name = "PATH",
        help = "Paths to reference_introns TSV file produced by iso-classify"
    )]
    pub refs: PathBuf,

    #[arg(
        short = 'q',
        long = "query",
        required = true,
        value_name = "PATH",
        help = "Path to BED12 file to classify"
    )]
    pub query: PathBuf,

    #[arg(
        short = 'b',
        long = "blacklist",
        required = false,
        value_name = "PATH",
        help = "Path to BED4 file with blacklisted introns"
    )]
    pub blacklist: Option<PathBuf>,

    #[arg(
        short = 't',
        long = "threads",
        help = "Number of threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get()
    )]
    pub threads: usize,

    #[arg(
        long = "recover",
        help = "Flag to recover from disputed retentions",
        value_name = "FLAG",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
    )]
    pub recover: bool,

    #[arg(
        long = "im",
        long = "in-memory",
        help = "Flag to avoid writing output files",
        value_name = "FLAG",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
    )]
    pub in_memory: bool,

    #[arg(
        short = 'p',
        long = "prefix",
        required = false,
        value_name = "PATH",
        help = "Prefix for output files"
    )]
    pub prefix: Option<PathBuf>,
}
