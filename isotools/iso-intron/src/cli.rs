use clap::Parser;
use std::path::PathBuf;
use thiserror::Error;

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
}

impl Args {
    pub fn check(&self) -> Result<(), CliError> {
        self.validate_args()
    }

    fn validate_args(&self) -> Result<(), CliError> {
        self.check_dbs()?;
        self.check_query()?;

        if !self.blacklist.is_empty() {
            self.check_blacklist()?;
        } else {
            log::warn!("No blacklist provided. Skipping...");
        }
        Ok(())
    }

    fn check_dbs(&self) -> Result<(), CliError> {
        if self.refs.is_empty() {
            let err = "No reference files provided".to_string();
            return Err(CliError::InvalidInput(err));
        }
        for db in &self.refs {
            validate(db)?;
        }

        if self.query.is_empty() {
            let err = "No query file provided".to_string();
            return Err(CliError::InvalidInput(err));
        }
        for query in &self.query {
            validate(query)?;
        }

        Ok(())
    }

    fn check_query(&self) -> Result<(), CliError> {
        validate(&self.query[0])
    }

    fn check_blacklist(&self) -> Result<(), CliError> {
        for bl in &self.blacklist {
            validate(bl)?;
        }
        Ok(())
    }
}

#[derive(Debug, Error)]
pub enum CliError {
    #[error("Invalid input: {0}")]
    InvalidInput(String),
    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),
}

pub fn validate(arg: &PathBuf) -> Result<(), CliError> {
    if !arg.exists() {
        return Err(CliError::InvalidInput(format!("{:?} does not exist", arg)));
    }

    if !arg.is_file() {
        return Err(CliError::InvalidInput(format!("{:?} is not a file", arg)));
    }

    match arg.extension() {
        Some(ext) if ext == "bed" => (),
        _ => {
            return Err(CliError::InvalidInput(format!(
                "file {:?} is not a BED file",
                arg
            )))
        }
    }

    match std::fs::metadata(arg) {
        Ok(metadata) if metadata.len() == 0 => {
            Err(CliError::InvalidInput(format!("file {:?} is empty", arg)))
        }
        Ok(_) => Ok(()),
        Err(e) => Err(CliError::IoError(e)),
    }
}
