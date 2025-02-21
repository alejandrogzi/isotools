use clap::Parser;
use config::ArgCheck;
use std::path::PathBuf;

#[derive(Debug, Parser)]
pub struct Args {
    #[arg(
        short = 'r',
        long = "raw",
        required = true,
        value_name = "PATHS",
        value_delimiter = ',',
        num_args = 1..,
        help = "Paths to raw Iso-Seq's BED12 files"
    )]
    pub raw: Vec<PathBuf>,

    #[arg(
        short = 'c',
        long = "calls",
        required = true,
        value_name = "PATHS",
        value_delimiter = ',',
        num_args = 1..,
        help = "Path to BED8 files with ORF calls"
    )]
    pub calls: Vec<PathBuf>,

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
}

impl ArgCheck for Args {
    fn get_blacklist(&self) -> &Vec<PathBuf> {
        &self.blacklist
    }

    fn get_ref(&self) -> &Vec<PathBuf> {
        &self.raw
    }

    fn get_query(&self) -> &Vec<PathBuf> {
        &self.calls
    }
}
