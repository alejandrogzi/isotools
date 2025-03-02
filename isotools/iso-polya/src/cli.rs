use clap::{ArgAction, Parser, Subcommand};
use config::ArgCheck;
use std::path::PathBuf;

pub const MIN_PER_ID: usize = 98;
pub const MAX_CLIP5: usize = 20;
pub const MAX_CLIP3: usize = 20;

// HMM parameters
pub const P2P: f32 = 0.9; // INFO: transition prob for polyA tail
pub const EMIT_A: f32 = 0.99; // INFO: emission prob for A in polyA state

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    #[command(subcommand)]
    pub command: SubArgs,

    #[arg(
        short = 't',
        long = "threads",
        help = "Number of threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get()
    )]
    pub threads: usize,
}

impl Args {}

#[derive(Debug, Subcommand)]
pub enum SubArgs {
    #[command(name = "aparent")]
    Aparent {
        #[command(flatten)]
        args: AparentArgs,
    },
    #[command(name = "filter")]
    Filter {
        #[command(flatten)]
        args: FilterArgs,
    },

    #[command(name = "caller")]
    Caller {
        #[command(flatten)]
        args: CallerArgs,
    },
}

#[derive(Debug, Parser)]
pub struct AparentArgs {
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

impl ArgCheck for AparentArgs {
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

#[derive(Debug, Parser)]
pub struct FilterArgs {
    #[arg(
        short = 's',
        long = "sam",
        required = true,
        value_name = "PATH",
        value_delimiter = ',',
        num_args = 1,
        help = "Path to .sam file"
    )]
    pub sam: Vec<PathBuf>,

    #[arg(
        long = "per_id",
        help = "Min %id (computed without the 5' and 3' clip). Must be [0-100] (percent).",
        value_name = "VALUE",
        default_value_t = MIN_PER_ID,
    )]
    pub per_id: usize,

    #[arg(
        long = "clip5",
        help = "Max 5' soft or hard clip",
        value_name = "VALUE",
        default_value_t = MAX_CLIP5,
    )]
    pub clip5: usize,

    #[arg(
        long = "clip3",
        help = "Max 3' soft or hard clip",
        value_name = "VALUE",
        default_value_t = MAX_CLIP3,
    )]
    pub clip3: usize,

    #[arg(
        long = "p2p",
        help = "Transition probability of looping in the polyA state ",
        value_name = "VALUE",
        default_value_t = P2P,
    )]
    pub p2p: f32,

    #[arg(
        long = "emit-a",
        help = "Probability of emitting A in the polyA state",
        value_name = "VALUE",
        default_value_t = EMIT_A,
    )]
    pub emit_a: f32,
}

impl ArgCheck for FilterArgs {
    // WARN: placeholder
    fn get_blacklist(&self) -> &Vec<PathBuf> {
        &self.sam
    }

    // WARN: placeholder
    fn get_ref(&self) -> &Vec<PathBuf> {
        &self.sam
    }

    // WARN: placeholder
    fn get_query(&self) -> &Vec<PathBuf> {
        &self.sam
    }
}

#[derive(Debug, Parser)]
pub struct CallerArgs {}
