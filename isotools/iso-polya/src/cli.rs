use clap::{ArgAction, Parser, Subcommand};
use config::ArgCheck;
use std::path::PathBuf;
use std::sync::Arc;

// Aparent parameters
pub const CHUNK_SIZE: usize = 500;

// Read parameters
pub const MIN_PER_ID: usize = 98;
pub const MAX_CLIP5: usize = 20;
pub const MAX_CLIP3: usize = 20;

// HMM parameters
pub const P2P: f32 = 0.9; // INFO: transition prob for polyA tail
pub const EMIT_A: f32 = 0.99; // INFO: emission prob for A in polyA state

// PASCaller parameters
pub const POLYA_LENGTH_THRESHOLD: usize = 50;
pub const GENOMIC_POLYA_THRESHOLD: usize = 5; // INFO: 5 A's in genome
pub const APARENT_THRESHOLD: f32 = 0.01;
pub const WIGGLE: usize = 2;

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

#[derive(Debug, Parser, Clone)]
pub struct AparentArgs {
    #[arg(
        short = 'b',
        long = "bed",
        required_if_eq("simulate","false"),
        value_name = "PATHS",
        value_delimiter = ',',
        num_args = 1..,
        help = "Path to raw Iso-Seq's BED12 files"
    )]
    pub bed: Vec<PathBuf>,

    #[arg(
        long = "twobit",
        conflicts_with = "simulate",
        required = false,
        value_name = "PATH",
        num_args = 1,
        help = "Path to genome 2bit file"
    )]
    pub twobit: Option<PathBuf>,

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

    #[arg(
        short = 'c',
        long = "chunk-size",
        help = "Chunk size for parallel processing",
        value_name = "CHUNK_SIZE",
        default_value_t = CHUNK_SIZE
    )]
    pub chunk_size: usize,

    #[arg(
        long = "simulate",
        required = false,
        value_name = "FLAG",
        help = "Simulate reads and run APARENT on them",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
        conflicts_with("bed"),
        conflicts_with("twobit"),
    )]
    pub simulate: bool,

    #[arg(
        long = "number-of-reads",
        value_name = "VALUE",
        default_value_t = 1000,
        help = "Number of reads to simulate",
        conflicts_with("bed"),
        conflicts_with("twobit")
    )]
    pub number_of_reads: usize,

    #[arg(
        long = "polya-range",
        value_name = "VALUE",
        help = "PolyA range for simulation",
        default_value_t = 10,
        conflicts_with("bed"),
        conflicts_with("twobit")
    )]
    pub polya_range: usize,

    #[arg(
        long = "read-length",
        value_name = "VALUE",
        help = "Read length for simulation",
        default_value_t = 300,
        conflicts_with("bed"),
        conflicts_with("twobit")
    )]
    pub read_length: usize,

    #[arg(
        long = "stranded",
        value_name = "FLAG",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        help = "Flag to simulate stranded reads",
        require_equals(true),
        action = ArgAction::Set,
        conflicts_with("bed"),
        conflicts_with("twobit")
    )]
    pub stranded: bool,
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

#[derive(Debug, Parser, Clone)]
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

    #[arg(
        long = "stat",
        required = false,
        value_name = "FLAG",
        help = "If set output a {input}.tsv file that contains statistics of all reads readID {tab} perID {tab} 5'clip {tab} 3'clip (non-PolyA part){tab} polyA tail length of 3'clip [optionally]{tab} polyA tail length of read",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
    )]
    pub stat: bool,

    #[arg(
        long = "keep",
        required = false,
        value_name = "FLAG",
        help = "If set, produce an {input}.TESBad5Prime.bed file which contains reads that align well but only have a 5' clip above our threshold. Can be used for PAScaller",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
    )]
    pub keep: bool,

    #[arg(
        long = "suffix",
        required = false,
        value_name = "FLAG",
        help = "If >0, compute and output the length of the polyA tail in the polyAReadSuffix + 3'clip_len suffix of the read (default don't do that)"
    )]
    pub suffix: Option<usize>,

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
        short = 'q',
        long = "queue",
        required = false,
        value_name = "QUEUE",
        help = "Queue to send jobs to",
        requires("para")
    )]
    pub queue: Option<String>,

    #[arg(
        short = 'm',
        long = "mem",
        required = false,
        value_name = "QUEUE",
        help = "Memory to send jobs to",
        requires("para")
    )]
    pub mem: Option<usize>,

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
pub struct CallerArgs {
    #[arg(
        short = 'a',
        long = "aparent",
        required = true,
        value_name = "PATH",
        help = "Path to output APARENT .bed file"
    )]
    pub aparent: PathBuf,

    #[arg(
        short = 'b',
        long = "bed",
        required = true,
        value_name = "PATH",
        help = "Path to output filterMinimapQuality.perl output .bed file"
    )]
    pub bed: PathBuf,

    #[arg(
        long = "toga",
        required = true,
        value_name = "PATH",
        help = "Path to output TOGA projections to localize polyA tails in the genome [CDS, UTR, CDS/UTR]"
    )]
    pub toga: Option<Vec<PathBuf>>,

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
        long = "filter",
        help = "Flag to filter out reads above/below a certain polyA tail length + APARENT score threshold",
        value_name = "FLAG",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
        conflicts_with("recover"),
    )]
    pub filter: bool,

    #[arg(
        long = "filter-type",
        help = "Filter out reads above/below a certain polyA tail length + APARENT score threshold",
        value_name = "FLAG",
        conflicts_with("recover"),
        default_value("above")
    )]
    pub filter_side: config::FilterSide,

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
}

impl CallerArgs {
    pub fn from(args: Arc<Vec<String>>) -> Self {
        let drop = vec!["--introns", "--bigwig", "--twobit"];

        let mut local_args = Vec::new();
        let mut iter = args.iter().peekable();

        while let Some(arg) = iter.next() {
            // INFO: skipping useless args + value
            if drop.contains(&arg.as_str()) {
                iter.next();
                continue;
            }

            if arg == "--query" {
                local_args.push("--bed".to_string());
            } else {
                local_args.push(arg.clone());
            }
        }

        let mut full_args = vec![env!("CARGO_PKG_NAME").to_string()];
        full_args.extend(local_args);
        full_args.push("--recover".to_string());

        CallerArgs::parse_from(full_args)
    }
}
