use dashmap::DashSet;
use indicatif::{ProgressBar, ProgressStyle};
use num_traits::{Num, NumCast};
use thiserror::Error;

use std::borrow::Borrow;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::time::Duration;

// os
#[cfg(not(windows))]
const TICK_SETTINGS: (&str, u64) = ("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏ ", 80);
#[cfg(windows)]
const TICK_SETTINGS: (&str, u64) = (r"+-x| ", 200);

/// return a pre-configured progress bar
pub fn get_progress_bar(length: u64, msg: &str) -> ProgressBar {
    let progressbar_style = ProgressStyle::default_spinner()
        .tick_chars(TICK_SETTINGS.0)
        .template(" {spinner} {msg:<30} {wide_bar} ETA {eta_precise} ")
        .expect("no template error");

    let progress_bar = ProgressBar::new(length);

    progress_bar.set_style(progressbar_style);
    progress_bar.enable_steady_tick(Duration::from_millis(TICK_SETTINGS.1));
    progress_bar.set_message(msg.to_owned());

    progress_bar
}

/// write a DashSet to a file
pub fn write_objs<T>(data: &DashSet<T>, fname: &str)
where
    T: AsRef<str> + Sync + Send + Eq + std::hash::Hash,
{
    log::info!("Reads in {}: {:?}. Writing...", fname, data.len());
    let f = match File::create(fname) {
        Ok(f) => f,
        Err(e) => panic!("Error creating file: {}", e),
    };
    let mut writer = BufWriter::new(f);

    for line in data.iter() {
        writeln!(writer, "{}", line.as_ref()).unwrap_or_else(|e| {
            panic!("Error writing to file: {}", e);
        });
    }
}

/// write any collection to a file
pub fn write_collection(data: &Vec<String>, fname: &str) {
    log::info!("Reads in {}: {:?}. Writing...", fname, data.len());
    let f = match File::create(fname) {
        Ok(f) => f,
        Err(e) => panic!("Error creating file: {}", e),
    };
    let mut writer = BufWriter::new(f);

    for line in data.iter() {
        writeln!(writer, "{}", line).unwrap_or_else(|e| {
            panic!("Error writing to file: {}", e);
        });
    }
}

/// argument checker for all subcommands
pub trait ArgCheck {
    fn check(&self) -> Result<(), CliError> {
        self.validate_args()
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

    fn check_dbs(&self) -> Result<(), CliError> {
        if self.get_ref().is_empty() {
            let err = "No reference files provided".to_string();
            return Err(CliError::InvalidInput(err));
        }
        for db in self.get_ref() {
            validate(db)?;
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

    fn check_blacklist(&self) -> Result<(), CliError> {
        for bl in self.get_blacklist() {
            validate(bl)?;
        }
        Ok(())
    }

    fn get_blacklist(&self) -> &Vec<PathBuf>;
    fn get_ref(&self) -> &Vec<PathBuf>;
    fn get_query(&self) -> &Vec<PathBuf>;
}

/// error handling for CLI
#[derive(Debug, Error)]
pub enum CliError {
    #[error("Invalid input: {0}")]
    InvalidInput(String),
    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),
}

/// argument validation
pub fn validate(arg: &PathBuf) -> Result<(), CliError> {
    if !arg.exists() {
        return Err(CliError::InvalidInput(format!(
            "ERROR: {:?} does not exist",
            arg
        )));
    }

    if !arg.is_file() {
        return Err(CliError::InvalidInput(format!(
            "ERROR: {:?} is not a file",
            arg
        )));
    }

    match arg.extension() {
        Some(ext) if ext == "bed" => (),
        _ => {
            return Err(CliError::InvalidInput(format!(
                "ERROR: file {:?} is not a BED file",
                arg
            )))
        }
    }

    match std::fs::metadata(arg) {
        Ok(metadata) if metadata.len() == 0 => Err(CliError::InvalidInput(format!(
            "ERROR: file {:?} is empty",
            arg
        ))),
        Ok(_) => Ok(()),
        Err(e) => Err(CliError::IoError(e)),
    }
}

// quality of life improvement fns
#[inline(always)]
pub fn exonic_overlap<N, I>(exons_a: &I, exons_b: &I) -> bool
where
    N: Num + NumCast + Copy + PartialOrd,
    I: IntoIterator,
    I::Item: Borrow<(N, N)>,
    for<'a> &'a I: IntoIterator<Item = &'a I::Item>,
{
    let mut iter_a = exons_a.into_iter();
    let mut iter_b = exons_b.into_iter();

    let mut exon_a = iter_a.next();
    let mut exon_b = iter_b.next();

    loop {
        match (exon_a, exon_b) {
            (Some(start_end_a), Some(start_end_b)) => {
                let (start_a, end_a) = start_end_a.borrow();
                let (start_b, end_b) = start_end_b.borrow();

                if *start_a < *end_b && *start_b < *end_a {
                    return true;
                }

                if *end_a < *end_b {
                    exon_a = iter_a.next();
                } else {
                    exon_b = iter_b.next();
                }
            }
            _ => break,
        }
    }

    false
}
