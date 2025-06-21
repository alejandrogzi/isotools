use anyhow::{bail, Result};
use dashmap::{DashMap, DashSet};
use hashbrown::{HashMap, HashSet};
use indicatif::{ProgressBar, ProgressStyle};
use log::{info, warn};
use num_traits::{Num, NumCast};
use rayon::prelude::*;
use serde_json::{Map, Value};
use thiserror::Error;

use std::borrow::Borrow;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::str::FromStr;
use std::sync::Arc;
use std::time::Duration;

use crate::{
    BedColumn, BedColumnValue, BedParser, CoordType, MatchType, ModuleMap, OverlapType,
    ParallelCollector, TsvParser,
};

// os
#[cfg(not(windows))]
const TICK_SETTINGS: (&str, u64) = ("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏ ", 80);
#[cfg(windows)]
const TICK_SETTINGS: (&str, u64) = (r"+-x| ", 200);

/// Return a pre-configured progress bar
///
/// # Arguments
///
/// * `length` - The length of the progress bar
/// * `msg` - The message to display
///
/// # Returns
///
/// * `ProgressBar` - A pre-configured progress bar
///
/// # Example
///
/// ```rust, no_run
/// use isotools::config::fns::get_progress_bar;
///
/// let length = 100;
/// let msg = "Processing...";
///
/// let progress_bar = get_progress_bar(length, msg);
///
/// assert_eq!(progress_bar.length(), length);
/// assert_eq!(progress_bar.message(), msg);
/// ```
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

/// Write a DashSet to a file
///
/// # Arguments
///
/// * `data` - The DashSet to write
/// * `fname` - The name of the file to write to
///
/// # Example
///
/// ```rust, no_run
/// use isotools::config::fns::write_objs;
///
/// let data = DashSet::new();
/// data.insert("line1".to_string());
/// data.insert("line2".to_string());
///
/// let fname = "output.txt";
///
/// write_objs(&data, fname);
/// ```
pub fn write_objs<T>(data: &DashSet<T>, fname: &str)
where
    T: AsRef<[u8]> + Sync + Send + Eq + std::hash::Hash,
{
    log::info!("Reads in {}: {:?}. Writing...", fname, data.len());
    let f = match File::create(fname) {
        Ok(f) => f,
        Err(e) => panic!("Error creating file: {}", e),
    };
    let mut writer = BufWriter::new(f);

    for line in data.iter() {
        let bytes = line.as_ref();
        writer
            .write_all(bytes)
            .unwrap_or_else(|e| panic!("ERROR: Error writing to file -> {e}"));
        writer
            .write_all(b"\n")
            .unwrap_or_else(|e| panic!("ERROR: Error newline to file -> {e}"));
    }
}

/// Write any collection to a file
///
/// # Arguments
///
/// * `data` - The collection to write
/// * `fname` - The name of the file to write to
///
/// # Example
///
/// ```rust, no_run
/// use isotools::config::fns::write_collection;
///
/// let data = vec!["line1".to_string(), "line2".to_string()];
/// let fname = "output.txt";
///
/// write_collection(&data, fname);
/// ```
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

/// Argument checker for all subcommands
pub trait ArgCheck {
    /// Check the arguments
    ///
    /// # Returns
    ///
    /// * `Result<(), CliError>` - Ok if the arguments are valid, Err if not
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isotools::config::fns::ArgCheck;
    ///
    /// let args = MyArgs::new();
    /// let result = args.check();
    ///
    /// assert!(result.is_ok());
    /// ```
    fn check(&self) -> Result<(), CliError> {
        self.validate_args()
    }

    /// Inner function to validate the arguments
    ///
    /// # Returns
    ///
    /// * `Result<(), CliError>` - Ok if the arguments are valid, Err if not
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isotools::config::fns::ArgCheck;
    ///
    /// let args = MyArgs::new();
    /// let result = args.validate_args();
    ///
    /// assert!(result.is_ok());
    /// ```
    fn validate_args(&self) -> Result<(), CliError> {
        self.check_dbs()?;

        if !self.get_blacklist().is_empty() {
            self.check_blacklist()?;
        } else {
            log::warn!("No blacklist provided. Skipping...");
        };

        Ok(())
    }

    /// Check input files for validity
    ///
    /// # Returns
    ///
    /// * `Result<(), CliError>` - Ok if the files are valid, Err if not
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isotools::config::fns::ArgCheck;
    ///
    /// let args = MyArgs::new();
    /// let result = args.check_dbs();
    ///
    /// assert!(result.is_ok());
    /// ```
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

    /// Check the blacklist files for validity
    ///
    /// # Returns
    ///
    /// * `Result<(), CliError>` - Ok if the files are valid, Err if not
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isotools::config::fns::ArgCheck;
    ///
    /// let args = MyArgs::new();
    /// let result = args.check_blacklist();
    ///
    /// assert!(result.is_ok());
    /// ```
    fn check_blacklist(&self) -> Result<(), CliError> {
        for bl in self.get_blacklist() {
            validate(bl)?;
        }
        Ok(())
    }

    /// Get the blacklist files
    ///
    /// # Returns
    ///
    /// * `&Vec<PathBuf>` - A reference to the vector of PathBufs representing the blacklist files
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isotools::config::fns::ArgCheck;
    ///
    /// let args = MyArgs::new();
    /// let blacklist = args.get_blacklist();
    ///
    /// assert!(!blacklist.is_empty());
    /// ```
    fn get_blacklist(&self) -> &Vec<PathBuf>;

    /// Get the reference files
    ///
    /// # Returns
    ///
    /// * `&Vec<PathBuf>` - A reference to the vector of PathBufs representing the reference files
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isotools::config::fns::ArgCheck;
    ///
    /// let args = MyArgs::new();
    /// let reference = args.get_ref();
    ///
    /// assert!(!reference.is_empty());
    /// ```
    fn get_ref(&self) -> &Vec<PathBuf>;

    /// Get the query files
    ///
    /// # Returns
    ///
    /// * `&Vec<PathBuf>` - A reference to the vector of PathBufs representing the query files
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use isotools::config::fns::ArgCheck;
    ///
    /// let args = MyArgs::new();
    /// let query = args.get_query();
    ///
    /// assert!(!query.is_empty());
    /// ```
    fn get_query(&self) -> &Vec<PathBuf>;
}

/// Rrror handling for CLI
///
/// This enum is used to handle
/// errors that may occur during
/// CLI operations.
#[derive(Debug, Error)]
pub enum CliError {
    #[error("Invalid input: {0}")]
    InvalidInput(String),
    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),
}

/// Argument validation
///
/// # Arguments
///
/// * `arg` - Path to the file to be validated
///
/// # Returns
///
/// * `Result<(), CliError>` - Ok if the file is valid, Err if not
///
/// # Example
///
/// ```rust, no_run
/// use isotools::config::fns::validate;
///
/// let path = PathBuf::from("path/to/file.bed");
/// let result = validate(&path);
///
/// if result.is_ok() {
///    println!("File is valid");
/// } else {
///   println!("File is invalid: {:?}", result.err());
/// }
/// ```
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
/// * `key` - the column to be used as the key for the outer HashMap
/// * `attribute` - the attribute to be stored in the inner HashMap
///
/// # Returns
///
/// * `DashMap<String, HashMap<String, BedColumnValue>>` - a DashMap where the key is the chromosome
/// * and the value is a HashMap of key-dependent values. User can specify which column should be used
/// * as the key and which attribute should be stored.
pub fn bed_to_nested_map<T>(
    contents: Arc<String>,
    key: BedColumn,
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

            // INFO: we need to use the key column to create a unique key for each entry
            // INFO: in case of duplicates, only BedColumn::Score stores multiple values!
            let inner_key = match key {
                BedColumn::Chrom => record.chrom().to_string(),
                BedColumn::Start => record.start().to_string(),
                BedColumn::End => record.end().to_string(),
                BedColumn::Name => record.name().to_string(),
                BedColumn::Score => record.score().to_string(),
                BedColumn::Strand => record.strand().to_string(),
                BedColumn::ThickStart => record.cds_start().to_string(),
                BedColumn::ThickEnd => record.cds_end().to_string(),
                _ => {
                    warn!(
                        "ERROR: Invalid key column for BED file. Is not possible
                        to use as a key either because is redundant or is not meaningful!"
                    );
                    return;
                }
            };

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
                .entry(inner_key)
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

/// Convert a BED file into a row-wise nested collection structure
///
/// This function is used to parse a BED file and convert it into a
/// DashMap of chromosome-dependent HashMaps. Each inner HashMap stores
/// the 'attributes' of each row in the BED file as a `Vec<BedColumnValue>`.
/// The user specifies which attributes should be stored.
///
/// # Arguments
///
/// * `contents` - the contents of the BED file
/// * `key` - the column to be used as the key for the outer HashMap
/// * `attributes` - the attributes to be stored in the inner HashMap
///
/// # Returns
///
/// * `DashMap<String, HashMap<String, Vec<BedColumnValue>>>` - a DashMap where the key is the chromosome
/// * and the value is a HashMap of key-dependent values. User can specify which column should be used
/// * as the key and which attribute should be stored.
pub fn bed_to_nested_collection<T>(
    contents: Arc<String>,
    key: BedColumn,
    attributes: Vec<BedColumn>,
) -> Result<DashMap<String, HashMap<String, Vec<BedColumnValue>>>, anyhow::Error>
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

            // INFO: we need to use the key column to create a unique key for each entry
            // INFO: in case of duplicates, only BedColumn::Score stores multiple values!
            let inner_key = match key {
                BedColumn::Chrom => record.chrom().to_string(),
                BedColumn::Start => record.start().to_string(),
                BedColumn::End => record.end().to_string(),
                BedColumn::Name => record.name().to_string(),
                BedColumn::Score => record.score().to_string(),
                BedColumn::Strand => record.strand().to_string(),
                BedColumn::ThickStart => record.cds_start().to_string(),
                BedColumn::ThickEnd => record.cds_end().to_string(),
                _ => {
                    warn!(
                        "ERROR: Invalid key column for BED file. Is not possible
                        to use as a key either because is redundant or is not meaningful!"
                    );
                    return;
                }
            };

            let mut values = Vec::new();
            for attribute in &attributes {
                if !attribute.is_valid() {
                    warn!(
                        "ERROR: Invalid attribute column for BED file: {:?}. Skipping...",
                        attribute
                    );
                    continue;
                }

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
                    BedColumn::BlockSizes => {
                        BedColumnValue::BlockSizes(record.block_sizes().clone())
                    }
                    BedColumn::BlockStarts => {
                        BedColumnValue::BlockStarts(record.block_starts().clone())
                    }
                };

                values.push(value);
            }

            let mut entry = tracks.entry(chrom).or_insert_with(HashMap::new);
            entry
                .entry(inner_key)
                .and_modify(|existing: &mut Vec<BedColumnValue>| existing.extend(values.clone()))
                .or_insert(values);

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

/// Writes the descriptor to a JSON file
///
/// # Arguments
///
/// * `descriptor` - Descriptor to write
/// * `path` - Path to the JSON file
///
/// # Example
///
/// ```rust, no_run
/// let descriptor = DashMap::new();
/// write_descriptor(&descriptor, "path/to/descriptor.json");
/// ```
pub fn write_descriptor_as_json(descriptor: &DashMap<String, Box<dyn ModuleMap>>, path: &str) {
    let mut json_map = Map::with_capacity(descriptor.len());

    descriptor.iter().for_each(|entry| {
        let (key, value) = entry.pair();
        json_map.insert(key.clone(), value.to_json());
    });

    // let json_map: Map<String, Value> = descriptor
    //     .iter()
    //     .par_bridge()
    //     .map(|entry| {
    //         let (key, value) = entry.pair();
    //         (key.clone(), value.to_json())
    //     })
    //     .collect::<HashMap<_, _>>() // INFO: intermediate parallel-friendly structure
    //     .into_iter()
    //     .collect(); // INFO: convert to serde_json::Map

    let file = File::create(path).expect("Unable to create file");
    let writer = BufWriter::new(file);

    serde_json::to_writer(writer, &Value::Object(json_map)).expect("Failed to write JSON");
}

/// Write results from a ParallelAccumulator to files in parallel
///
/// # Arguments
///
/// * `accumulator` - The accumulator containing the results
/// * `filenames` - A vector of PathBufs representing the filenames
/// * `outdir` - An optional PathBuf representing the output directory
///
/// # Example
///
/// ```rust, no_run
/// use isotools::config::fns::par_write_results;
///
/// let accumulator = ParallelAccumulator::default();
/// let filenames = vec![PathBuf::from("file1.txt"), PathBuf::from("file2.txt")];
/// let outdir = Some(PathBuf::from("output"));
///
/// par_write_results(accumulator, filenames, outdir);
/// ```
pub fn par_write_results<K: ParallelCollector>(
    accumulator: &K,
    filenames: Vec<PathBuf>,
    outdir: Option<PathBuf>,
) {
    if accumulator.len() != filenames.len() {
        log::error!("ERROR: Number of filenames does not match the number of results!");
        std::process::exit(1);
    }

    let collections = accumulator
        .get_collections()
        .expect("ERROR: Could not get collections!");

    if let Some(outdir) = outdir {
        let files = filenames
            .iter()
            .map(|filename| outdir.join(filename))
            .collect::<Vec<_>>();

        write_pairs(collections, files);
    } else {
        write_pairs(collections, filenames);
    }
}

/// Inner function to write pairs of collections to files
///
/// # Arguments
///
/// * `collections` - A vector of DashSet<String> representing the collections
/// * `filenames` - A vector of PathBufs representing the filenames
///
/// # Example
///
/// ```rust, no_run
/// use isotools::config::fns::write_pairs;
///
/// let collections = vec![DashSet::new(), DashSet::new()];
/// let filenames = vec![PathBuf::from("file1.txt"), PathBuf::from("file2.txt")];
///
/// write_pairs(collections, filenames);
/// ```
fn write_pairs(collections: Vec<&DashSet<String>>, filenames: Vec<PathBuf>) {
    collections
        .par_iter()
        .zip(filenames.par_iter())
        .for_each(|(collection, path)| {
            if collection.is_empty() {
                log::warn!("WARN: A collection from the accumulator is empty! Skipping...");
                return;
            }

            write_objs(
                collection,
                path.to_str().expect("ERROR: Invalid path to write!"),
            );
        });
}

/// Write a descriptor to a .tsv file
///
/// # Arguments
///
/// * `descriptor` - Descriptor to write
/// * `path` - Path to the .tsv file
///
/// # Example
///
/// ```rust, no_run
/// let descriptor = DashMap::new();
/// write_descriptor(&descriptor, "path/to/descriptor.tsv");
/// ```
pub fn write_descriptor(descriptor: &DashMap<String, Box<dyn ModuleMap>>, path: &str) {
    // Get all keys and column names (assume all descriptors share same schema)
    let mut rows = Vec::with_capacity(descriptor.len());
    let mut all_columns: Vec<String> = Vec::new();

    for entry in descriptor.iter() {
        let (id, value) = entry.pair();
        let json = value.to_json();

        if let Value::Object(obj) = json {
            // INFO: cache the first object's keys as columns
            if all_columns.is_empty() {
                all_columns = obj.keys().cloned().collect();
            }

            let row = all_columns
                .iter()
                .map(|k| json_get_field(&obj, k))
                .collect::<Vec<String>>();

            rows.push((id.clone(), row));
        } else {
            panic!("Expected object in descriptor JSON");
        }
    }

    let file = File::create(path).expect("Unable to create file");
    let mut writer = BufWriter::new(file);

    writeln!(writer, "id\t{}", all_columns.join("\t")).expect("Write failed");

    for (id, values) in rows {
        writeln!(writer, "{}\t{}", id, values.join("\t")).expect("Write failed");
    }
}

/// Convert a JSON value to string, escaping tabs and newlines
///
/// # Arguments
///
/// * `obj` - The JSON object
/// * `key` - The key to retrieve the value
///
/// # Returns
///
/// * `String` - The string representation of the value
///
/// # Example
///
/// ```rust, no_run
/// use isotools::config::fns::json_get_field;
///
/// let obj = serde_json::json!({
///    "key1": "value1",
///    "key2": 42,
///    "key3": true,
///    "key4": null,
///    "key5": ["array", "of", "values"],
///    "key6": {"nested": "object"}
/// });
///
/// json_get_field(&obj.as_object().unwrap(), "key1");
/// ```
fn json_get_field(obj: &Map<String, Value>, key: &str) -> String {
    let value = obj.get(key).unwrap_or(&Value::Null);
    let s = match value {
        Value::Null => "NULL".to_string(),
        Value::Bool(b) => b.to_string(),
        Value::Number(n) => n.to_string(),
        Value::String(s) => s.clone(),
        _ => serde_json::to_string(value).unwrap_or_default(), // INFO: for arrays or objects
    };
    s.replace('\t', "   ").replace('\n', " ")
}
