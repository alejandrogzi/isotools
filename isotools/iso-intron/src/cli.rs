use clap::{ArgAction, Parser};
use config::ArgCheck;
use std::path::PathBuf;

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
        short = 'b',
        long = "blacklist",
        required = false,
        value_name = "PATH",
        value_delimiter = ',',
        num_args = 1..,
        help = "Path to BED4 file with blacklisted introns"
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
        short = 'p',
        long = "plot",
        help = "Flag to output retentions in a BED4 file",
        value_name = "FLAG",
        default_value = "false"
    )]
    pub plot: bool,

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
        long = "spliceai",
        required = false,
        value_name = "PATH",
        num_args = 1,
        help = "Path to spliceAI directory [will asume 2 files per strand: acceptor and donor .bw]"
    )]
    pub splice_scores: Option<PathBuf>,

    #[arg(
        long = "toga",
        required = false,
        value_name = "PATH",
        num_args = 1..,
        help = "Path to BED12 TOGA annotation file"
    )]
    pub toga: Option<PathBuf>,
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
