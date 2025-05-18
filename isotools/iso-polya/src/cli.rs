use clap::{ArgAction, Parser, Subcommand};
use config::*;
use std::path::PathBuf;
use std::sync::Arc;

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
    #[command(name = "segment")]
    Segment {
        #[command(flatten)]
        args: SegmentArgs,
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
        value_name = "THRESHOLD",
        default_value_t = 0.005
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

impl AparentArgs {
    pub fn from(args: Vec<String>) -> Self {
        AparentArgs::parse_from(args)
    }
}

#[derive(Debug, Parser, Clone)]
pub struct SegmentArgs {
    #[arg(
        short = 'b',
        long = "bam",
        required = true,
        value_name = "PATH",
        value_delimiter = ',',
        num_args = 1,
        help = "Path to .bam file"
    )]
    pub bam: PathBuf,

    #[arg(
        short = 'I',
        long = "identity",
        help = "Min %id (computed without the 5' and 3' clip). Must be [0-100] (percent).",
        value_name = "VALUE",
        default_value_t = IDENTITY_THRESHOLD,
    )]
    pub identity: f32,

    #[arg(
        short = 'i',
        long = "min-identity",
        help = "Mininum % identity, used to discard extremely low reads
        (computed without the 5' and 3' clip). Must be [0-100] (percent).",
        value_name = "VALUE",
        default_value_t = MINIMUM_IDENTITY,
    )]
    pub min_identity: f32,

    #[arg(
        short = 's',
        long = "tail-suffix",
        help = "Suffix of the tail. Determines how many bases to consider at the end of the read",
        value_name = "VALUE",
        default_value_t = POLYA_SUFFIX,
    )]
    pub tail_suffix: usize,

    #[arg(
        short = 'S',
        long = "step-size",
        help = "Determines how many bp to move backwards in the read to find the tail",
        value_name = "VALUE",
        default_value_t = SUFFIX_STEP_SIZE,
    )]
    pub suffix_step_size: usize,

    #[arg(
        short = 'f',
        long = "clip5",
        help = "Max 5' soft or hard clip",
        value_name = "VALUE",
        default_value_t = MAX_CLIP5,
    )]
    pub max_clip_five: usize,

    #[arg(
        short = 't',
        long = "clip3",
        help = "Max 3' soft or hard clip",
        value_name = "VALUE",
        default_value_t = MAX_CLIP3,
    )]
    pub max_clip_three: usize,

    #[arg(
        short = 'P',
        long = "p2p",
        help = "Transition probability of looping in the polyA state ",
        value_name = "VALUE",
        default_value_t = P2P,
    )]
    pub p2p: f64,

    #[arg(
        short = 'E',
        long = "emit-a",
        help = "Probability of emitting A in the polyA state",
        value_name = "VALUE",
        default_value_t = EMIT_A,
    )]
    pub emit_a: f64,

    #[arg(
        short = 'T',
        long = "tag",
        help = "Flag to tag read names with polyA information",
        value_name = "FLAG",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
    )]
    pub tag: bool,

    #[arg(
        short = 'B',
        long = "bed",
        help = "Flag to convert bam to bed",
        value_name = "FLAG",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
    )]
    pub bed: bool,

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
        long = "prefix",
        required = false,
        value_name = "FILE_PREFIX",
        num_args = 1,
        help = "File name prefix",
        default_value = "file"
    )]
    pub prefix: PathBuf,
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

    #[arg(
        short = 'p',
        long = "prefix",
        required = false,
        value_name = "PATH",
        help = "Prefix for output files"
    )]
    pub prefix: Option<PathBuf>,
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
