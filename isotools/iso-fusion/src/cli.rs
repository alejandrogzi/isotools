use clap::{ArgAction, Parser};
use config::{validate, ArgCheck, CliError, OverlapType};
use std::path::PathBuf;
use std::sync::Arc;

#[derive(Debug, Parser)]
pub struct Args {
    #[arg(
        short = 'r',
        long = "ref",
        required = true,
        value_name = "PATH",
        help = "Path to BED12/Isoform [--map] file with rule transcripts"
    )]
    pub refs: Vec<PathBuf>,

    #[arg(
        short = 'q',
        long = "query",
        required = true,
        value_name = "PATHS",
        value_delimiter = ',',
        num_args = 1..,
        help = "Path to BED12 file to classify"
    )]
    pub query: Vec<PathBuf>,

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
        help = "Flag to recover from disputed fusions",
        value_name = "FLAG",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
    )]
    pub recover: bool,

    #[arg(
        long = "intron-match",
        help = "Flag to intron-specific match instead splicing match",
        value_name = "FLAG",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
    )]
    pub intron_match: bool,

    #[arg(
        short = 'm',
        long = "map",
        help = "Flag to read an isoforms file [tx->gene] instead of BED12",
        value_name = "FLAG",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
    )]
    pub map: bool,

    #[arg(
        short = 'o',
        long = "overlap-type",
        help = "Type of overlap to consider",
        value_name = "OVERLAP TYPE",
        required = false,
        requires_if("map", "true"),
        default_value("exon")
    )]
    pub overlap_type: OverlapType,

    #[arg(
        short = 'b',
        long = "blacklist",
        required = false,
        value_name = "PATH",
        value_delimiter = ',',
        num_args = 1..,
        help = "Path to BED12 file with blacklisted reads"
    )]
    pub blacklist: Vec<PathBuf>,

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
        help = "Prefix for output files",
        default_value("")
    )]
    pub prefix: PathBuf,
}

impl ArgCheck for Args {
    fn check_dbs(&self) -> Result<(), CliError> {
        if self.get_ref().is_empty() {
            log::warn!("No reference hints provided. Running in frequency mode...");
        } else {
            for db in self.get_ref() {
                validate(db)?;
            }
        }

        if self.get_query().is_empty() {
            let err = "No query file provided".to_string();
            return Err(CliError::InvalidInput(err));
        }
        for query in self.get_query() {
            validate(query)?;
        }

        Ok(())
    }

    fn validate_args(&self) -> Result<(), CliError> {
        self.check_dbs()?;

        if !self.get_blacklist().is_empty() {
            self.check_blacklist()?;
        } else {
            log::warn!("No blacklist provided. Skipping...");
        };

        Ok(())
    }

    fn get_blacklist(&self) -> &Vec<PathBuf> {
        &self.blacklist
    }

    fn get_ref(&self) -> &Vec<PathBuf> {
        &self.refs
    }

    fn get_query(&self) -> &Vec<PathBuf> {
        &self.query
    }
}

impl Args {
    pub fn from(args: Arc<Vec<String>>) -> Self {
        let drop = vec!["--introns", "--aparent", "--bigwig", "--twobit"];

        let mut local_args = Vec::new();
        let mut iter = args.iter().peekable();

        while let Some(arg) = iter.next() {
            // INFO: skipping useless args + value
            if drop.contains(&arg.as_str()) {
                iter.next();
                continue;
            }

            if arg == "--toga" {
                local_args.push("--ref".to_string());
            } else {
                local_args.push(arg.clone());
            }
        }

        let mut full_args = vec![env!("CARGO_PKG_NAME").to_string()];
        full_args.extend(local_args);
        full_args.push("--recover".to_string());

        Args::parse_from(full_args)
    }
}
