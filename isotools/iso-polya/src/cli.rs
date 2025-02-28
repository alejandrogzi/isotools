use clap::{ArgAction, Parser};
use config::ArgCheck;
use std::path::PathBuf;

#[derive(Debug, Parser)]
pub struct Args {
    #[arg(
        short = 'b',
        long = "bed",
        required = true,
        value_name = "PATHS",
        value_delimiter = ',',
        num_args = 1..,
        help = "Path to raw Iso-Seq's BED12 files"
    )]
    pub bed: Vec<PathBuf>,

    #[arg(
        long = "twobit",
        required = true,
        value_name = "PATH",
        num_args = 1,
        help = "Path to genome 2bit file"
    )]
    pub twobit: PathBuf,

    #[arg(
        short = 'b',
        long = "blacklist",
        required = false,
        value_name = "PATH",
        value_delimiter = ',',
        num_args = 1..,
        help = "Path to BED4 file with blacklisted reads"
    )]
    pub blacklist: Vec<PathBuf>,

    #[arg(
        short = 't',
        long = "threads",
        help = "Number of threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get()
    )]
    pub threads: usize,

    #[arg(
        long = "threshold",
        help = "Number of threads",
        value_name = "THREADS",
        default_value_t = 0.01
    )]
    pub threshold: f32,

    #[arg(
        long = "use-pf",
        required = false,
        value_name = "FLAG",
        help = "Use APARENT peak finder from the scipy package",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
    )]
    pub pf: bool,

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
        long = "para",
        required = false,
        value_name = "FLAG",
        help = "Send jobs to Hillerlab cluster",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
    )]
    pub para: bool,

    #[arg(
        long = "use-max-peak",
        required = false,
        value_name = "FLAG",
        help = "Use maximum peak instead of the APARENT peak",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
    )]
    pub use_max_peak: bool,
}

impl ArgCheck for Args {
    fn get_blacklist(&self) -> &Vec<PathBuf> {
        &self.blacklist
    }

    fn get_ref(&self) -> &Vec<PathBuf> {
        &self.bed
    }

    // WARN: placeholder
    fn get_query(&self) -> &Vec<PathBuf> {
        &self.bed
    }
}
