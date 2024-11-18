use clap::{ArgAction, Parser};
use config::{validate, ArgCheck, CliError};
use std::path::PathBuf;

#[derive(Debug, Parser)]
pub struct Args {
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
        short = 'h',
        long = "hint",
        required = false,
        value_name = "PATH",
        value_delimiter = ',',
        num_args = 1..,
        help = "Path to BED12 file with rule transcripts"
    )]
    pub hint: Vec<PathBuf>,

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
        long = "cds",
        help = "Flag to skip UTRs in chimeric regions when --hint",
        value_name = "FLAG",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
        requires("hint"),
    )]
    pub cds: bool,

    #[arg(
        short = 'w',
        long = "write",
        help = "Flag to write intergenic.bed with intergenic regions",
        value_name = "FLAG",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
    )]
    pub write: bool,
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
        &self.hint
    }

    fn get_query(&self) -> &Vec<PathBuf> {
        &self.query
    }
}
