use clap::{ArgAction, Parser};
use config::ArgCheck;
use std::path::PathBuf;
use std::sync::Arc;

#[derive(Debug, Parser)]
pub struct Args {
    #[arg(
        short = 'r',
        long = "ref",
        required = true,
        value_name = "PATHS",
        value_delimiter = ',',
        num_args = 1..,
        help = "Paths to BED12 files delimited by comma"
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
        long = "ignore-exon",
        help = "Flag to ignore upstream 5' end",
        value_name = "FLAG",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
    )]
    pub skip_exon: bool,

    #[arg(
        short = 'b',
        long = "blacklist",
        required = false,
        value_name = "PATH",
        value_delimiter = ',',
        num_args = 1..,
        help = "Path to BED12 file with blacklisted transcripts"
    )]
    pub blacklist: Vec<PathBuf>,

    #[arg(
        long = "recover",
        help = "Flag to recover from disputed truncations",
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

impl ArgCheck for Args {
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
        let drop = vec!["--introns", "--toga", "--aparent", "--bigwig", "--twobit"];

        let mut local_args = Vec::new();
        let mut iter = args.iter().peekable();

        while let Some(arg) = iter.next() {
            // INFO: skipping useless args + value
            if drop.contains(&arg.as_str()) {
                iter.next();
                continue;
            }

            if arg == "--query" {
                local_args.push("--query".to_string());

                let val = iter
                    .next()
                    .expect("ERROR: Missing value for --query")
                    .to_owned();

                local_args.push(val.clone());
                local_args.push("--ref".to_string());
                local_args.push(val);
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
