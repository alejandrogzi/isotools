use anyhow::{bail, Result};
use dashmap::{DashMap, DashSet};
use hashbrown::{HashMap, HashSet};
use indicatif::{ProgressBar, ProgressStyle};
use log::{info, warn};
use num_traits::{Num, NumCast};
use rayon::prelude::*;
use thiserror::Error;

use std::borrow::Borrow;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::str::FromStr;
use std::sync::Arc;
use std::time::Duration;

use crate::{BedColumn, BedColumnValue, BedParser, CoordType, MatchType, OverlapType, TsvParser};

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
        Some(ext) if ext == "bed" || ext == "tsv" => (),
        _ => {
            return Err(CliError::InvalidInput(format!(
                "ERROR: file {:?} is not a BED or TSV file",
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
/// Exonic overlap between two sets of exons.
///
/// This function is used to determine if two
/// sets of exons overlap by comparing the start
/// and end positions of each exon. The sets of
/// exons could be any type of collection.
///
/// # Arguments
/// * `exons_a` - a collection of exons
/// * `exons_b` - a collection of exons
///
/// # Returns
/// * `bool` - true if there is an overlap, false otherwise
///
/// # Example
/// ```rust
/// use isotools::config::fns::exonic_overlap;
///
/// let exons_a = vec![(1, 10), (20, 30)];
/// let exons_b = vec![(5, 15), (25, 35)];
/// assert_eq!(exonic_overlap(&exons_a, &exons_b), true);
/// ```
#[inline(always)]
pub fn exonic_overlap<N, I1, I2>(exons_a: &I1, exons_b: &I2) -> bool
where
    N: Num + NumCast + Copy + PartialOrd,
    I1: IntoIterator,
    I2: IntoIterator,
    I1::Item: Borrow<(N, N)>,
    I2::Item: Borrow<(N, N)>,
    for<'a> &'a I1: IntoIterator<Item = &'a I1::Item>,
    for<'a> &'a I2: IntoIterator<Item = &'a I2::Item>,
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

/// Splice site overlap between two sets of introns.
///
/// This function is used to determine if two sets of
/// introns overlap by comparing the start and end
/// positions of each intron. The latter will depend
/// on the match type.
///
/// # Arguments
/// * `introns_a` - a collection of introns
/// * `introns_b` - a collection of introns
/// * `match_type` - the match type [Intron, SpliceSite]
///
/// # Returns
/// * `bool` - true if there is an overlap, false otherwise
///
/// # Example
/// ```rust
/// use isotools::config::fns::splice_site_overlap;
/// use hashbrown::HashSet;
///
/// let introns_a = vec![(1, 10), (20, 30)];
/// let introns_b = [(5, 15), (25, 35)].iter().cloned().collect::<HashSet<_>>();
/// assert_eq!(splice_site_overlap(&introns_a, &introns_b), true);
/// ```
#[inline(always)]
pub fn splice_site_overlap<N>(
    introns_a: &Vec<(N, N)>,
    introns_b: &HashSet<(N, N)>,
    match_type: MatchType,
) -> bool
where
    N: Num + NumCast + Copy + PartialOrd + Eq + std::hash::Hash,
{
    match match_type {
        MatchType::Intron => introns_a.iter().any(|a| introns_b.contains(a)),
        MatchType::SpliceSite => {
            let splice_sites: HashSet<N> = introns_b.iter().flat_map(|&(x, y)| [x, y]).collect();
            introns_a
                .iter()
                .any(|&(x, y)| splice_sites.contains(&x) || splice_sites.contains(&y))
        }
    }
}

/// Convert a BED file into a HashMap of chrom-dependent intervals
///
/// This function is used to parse a BED file and convert it into a
/// HashMap of chromosome-dependent intervals. The intervals are
/// stored as a HashSet of tuples containing the start and end
/// positions of each interval.
///
/// # Arguments
///
/// * `contents` - the contents of the BED file
/// * `hint` - the coordinate type [Bounds, Intronic, Exonic]
///
/// # Returns
///
/// * `HashMap<String, HashSet<(u64, u64)>>` - a HashMap of chromosome-dependent intervals
///
/// # Example
///
/// ```rust, no_run
/// use isotools::config::fns::bed_to_map;
/// use isotools::config::BedParser;
/// use isotools::config::CoordType;
/// use isotools::config::OverlapType;
/// use std::sync::Arc;
/// use hashbrown::{HashSet, HashMap};
///
/// let contents = Arc::new("chr1\t10\t20\nchr2\t30\t40".to_string());
/// let hint = CoordType::Bounds;
/// let result = bed_to_map::<BedParser>(contents, hint).unwrap();
///
/// let mut expected = HashMap::new();
/// expected.insert("chr1".to_string(), [(10, 20)].iter().cloned().collect::<HashSet<_>>());
///
/// assert_eq!(result, expected);
/// ```
pub fn bed_to_map<T>(
    contents: Arc<String>,
    hint: CoordType,
) -> Result<HashMap<String, HashSet<(u64, u64)>>, anyhow::Error>
where
    T: BedParser,
{
    let pb = get_progress_bar(contents.lines().count() as u64, "Parsing BED files...");
    let tracks = contents
        .par_lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| {
            T::parse(line, OverlapType::Exon, false) // INFO: placeholders
                .map_err(|e| warn!("Error parsing {}: {}", line, e))
                .ok()
        })
        .fold(
            || HashMap::new(),
            |mut acc: HashMap<String, HashSet<(u64, u64)>>, record| {
                let entry = acc.entry(record.chrom().to_owned()).or_default();
                match hint {
                    CoordType::Bounds => {
                        entry.insert(record.coord());
                        ()
                    }
                    CoordType::Intronic => entry.extend(record.intronic_coords()),
                    CoordType::Exonic => entry.extend(record.exonic_coords()),
                };

                pb.inc(1);
                acc
            },
        )
        .reduce(
            || HashMap::new(),
            |mut acc, map| {
                for (k, v) in map {
                    acc.entry(k).or_default().extend(v);
                }
                acc
            },
        );

    pb.finish_and_clear();

    if tracks.is_empty() {
        anyhow::bail!("No tracks found in the provided file!");
    }
    info!(
        "Parsed {} blacklisted intervals.",
        tracks.values().flatten().count()
    );
    Ok(tracks)
}

/// Convert a BED file into a row-wise nested structure
///
/// This function is used to parse a BED file and convert it into a
/// DashMap of chromosome-dependent HashMaps. Each inner HashMap stores
/// the 'attribute' of each row in the BED file as a `BedColumnValue`.
/// The user specifies which attribute should be stored.
///
/// # Arguments
///
/// * `contents` - the contents of the BED file
/// * `attribute` - the attribute to be stored in the inner HashMap
///
/// # Returns
///
/// * `DashMap<String, HashMap<String, BedColumnValue>>` - a DashMap where the key is the chromosome
///   and the value is a HashMap mapping unique region names to the chosen `BedColumnValue`.
pub fn bed_to_nested_map<T>(
    contents: Arc<String>,
    attribute: BedColumn,
) -> Result<DashMap<String, HashMap<String, BedColumnValue>>, anyhow::Error>
where
    T: BedParser,
{
    let pb = get_progress_bar(contents.lines().count() as u64, "Parsing BED files...");
    let tracks = DashMap::new();

    contents
        .par_lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| {
            T::parse(line, OverlapType::Exon, false) // WARN: placeholders
                .map_err(|e| warn!("Error parsing {}: {}", line, e))
                .ok()
        })
        .for_each(|record| {
            let chrom = record.chrom().to_owned();
            let name = record.name().to_owned(); // INFO: name is unique (or should be)!

            let value = match attribute {
                BedColumn::Chrom => BedColumnValue::Chrom(record.chrom().to_string()),
                BedColumn::Start => BedColumnValue::Start(record.start()),
                BedColumn::End => BedColumnValue::End(record.end()),
                BedColumn::Name => BedColumnValue::Name(record.name().to_string()),
                BedColumn::Score => BedColumnValue::Score(vec![record.score()]),
                BedColumn::Strand => BedColumnValue::Strand(record.strand()),
                BedColumn::ThickStart => BedColumnValue::ThickStart(record.cds_start()),
                BedColumn::ThickEnd => BedColumnValue::ThickEnd(record.cds_end()),
                BedColumn::ItemRgb => BedColumnValue::ItemRgb(record.rgb().to_string()),
                BedColumn::BlockCount => BedColumnValue::BlockCount(record.block_count()),
                BedColumn::BlockSizes => BedColumnValue::BlockSizes(record.block_sizes().clone()),
                BedColumn::BlockStarts => {
                    BedColumnValue::BlockStarts(record.block_starts().clone())
                }
            };

            let mut entry = tracks.entry(chrom).or_insert_with(HashMap::new);
            entry
                .entry(name)
                .and_modify(|existing| match (existing, &value) {
                    // INFO: trick to handle duplicated rows with different scores!
                    (BedColumnValue::Score(scores), BedColumnValue::Score(new_scores)) => {
                        scores.extend(new_scores);
                    }
                    _ => {} // INFO: do nothing for other attributes
                })
                .or_insert(value);

            pb.inc(1);
        });

    pb.finish_and_clear();

    if tracks.is_empty() {
        anyhow::bail!("ERROR: No tracks found in the provided file!");
    }

    info!(
        "INFO: Parsed {} rows from the BED file!",
        tracks
            .iter()
            .map(|entry| entry.value().len())
            .sum::<usize>()
    );

    Ok(tracks)
}

/// Convert a TSV file into a HashMap of key-dependent values
///
/// This function is used to parse a TSV file and convert it into a
/// HashMap of key-dependent values. The user specifies which columns
/// should be used as the key and value.
///
/// # Arguments
///
/// * `contents` - the contents of the TSV file
/// * `key` - the column index to be used as the key
/// * `value` - the column index to be used as the value
///
/// # Returns
///
/// * `HashMap<String, Vec<V>>` - a HashMap of key-dependent values
///
/// # Example
///
/// ```rust, no_run
/// use isotools::config::fns::tsv_to_map;
/// use isotools::config::TsvParser;
/// use std::sync::Arc;
/// use hashbrown::HashMap;
///
/// struct MyParser;
///
/// impl TsvParser for MyParser {
///    fn parse(line: &str) -> Result<Self, anyhow::Error> {
///       let fields: Vec<&str> = line.split('\t').collect();
///      Ok(MyParser)
///   }
///
///  fn key(&self, idx: usize) -> &str {
///     "key"
/// }
///
/// fn value<V>(&self, idx: usize) -> Result<V, anyhow::Error>
/// where
///    V: FromStr,
///   <V as FromStr>::Err: std::fmt::Debug,
/// {
///    Ok("value".parse().unwrap())
/// }
///
/// }
///
/// let contents = Arc::new("key\tvalue".to_string());
/// let result = tsv_to_map::<MyParser, String>(contents, 0, 1).unwrap();
///
/// let mut expected = HashMap::new();
/// expected.insert("key".to_string(), vec!["value".to_string()]);
///
/// assert_eq!(result, expected);
/// ```
pub fn tsv_to_map<T, V>(
    contents: Arc<String>,
    key: usize,
    value: usize,
) -> Result<HashMap<String, Vec<V>>, anyhow::Error>
where
    T: TsvParser + Send + Sync,
    V: FromStr + Send + Sync,
    <V as FromStr>::Err: std::fmt::Debug,
{
    let pb = get_progress_bar(contents.lines().count() as u64, "Parsing TSV file...");

    let tracks = contents
        .par_lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|row| match T::parse(row) {
            Ok(record) => match record.value::<V>(value) {
                Ok(val) => Some((record.key(key).to_string(), val)),
                Err(e) => {
                    warn!("ERROR: Could not parse value from '{}': {:?}", row, e);
                    None
                }
            },
            Err(e) => {
                warn!("ERROR: Could not parse '{}': {}", row, e);
                None
            }
        })
        .fold(
            || HashMap::new(),
            |mut acc: HashMap<String, Vec<V>>, (k, v)| {
                acc.entry(k).or_default().push(v);
                pb.inc(1);
                acc
            },
        )
        .reduce(
            || HashMap::new(),
            |mut acc, map| {
                for (k, v) in map {
                    acc.entry(k).or_default().extend(v);
                }
                acc
            },
        );

    pb.finish_and_clear();

    if tracks.is_empty() {
        bail!("No tracks found in the provided file!");
    }

    info!(
        "Parsed {} records from file!",
        tracks.values().map(Vec::len).sum::<usize>()
    );

    Ok(tracks)
}
