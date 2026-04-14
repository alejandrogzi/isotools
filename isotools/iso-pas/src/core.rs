// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Core module for detecting real 3' ends on isoseq reads
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main function for detecting real 3' ends
//! on isoseq reads and processing the components of reads in parallel.
//!
//! In short, uses APARENT and HMM polyA scoring systems to define real ends.  
//! Finally, the caller module groups all the previous information and tries to
//! determine the intraprimming potential for each read,

use bigtools::{utils::reopen::Reopen, BigWigRead, Value};
use dashmap::{DashMap, DashSet};
use genepred::{GenePred, Strand};
use hashbrown::{HashMap, HashSet};
use packbed::{OverlapType, Role};
use rayon::prelude::*;

use std::{
    fs::File,
    io::{BufWriter, Write},
    str::from_utf8_unchecked,
    sync::{
        atomic::{AtomicU32, Ordering},
        Mutex,
    },
};

use crate::cli::*;

// tag separators
pub const SEP: &str = "#"; // WARN: should change to ':' -> breaks ORF caller [translationai hdf5]
pub const BIG_SEP: &str = "__"; // WARN: should change to '::' -> breaks ORF caller [translationai hdf5]

pub fn pas_caller(args: Args) -> Result<(), Box<dyn std::error::Error>> {
    log::info!("INFO: Starting PAS caller...");

    let tracks = if let Some(toga) = &args.refs {
        packbed::pack(
            vec![toga, &args.query],
            vec![Role::Reference, Role::Query],
            OverlapType::Exon,
        )
    } else {
        packbed::pack(vec![&args.query], vec![Role::Query], OverlapType::Exon)
    }
    .unwrap_or_else(|e| {
        log::error!("{}", e);
        std::process::exit(1);
    });

    let mut chrs = Vec::with_capacity(tracks.len());
    tracks.iter().for_each(|bucket| {
        log::debug!(
            "DEBUG: Bucket key: {}, value size: {}",
            bucket.key(),
            bucket.value().len()
        );

        if !bucket.value().is_empty() {
            let bind = bucket.key().clone();
            let chr = bind
                .split(':')
                .next()
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: Could not get chromosome from key -> {:?}",
                        bucket.key()
                    )
                })
                .as_bytes()
                .to_vec();

            chrs.push(chr);
        }
    });

    log::info!(
        "INFO: Collecting data for {} chromosome keys...",
        chrs.len()
    );

    let aparent_scores = load_aparent(&args.aparent_forward, &args.aparent_reverse, chrs);
    let counter = ParallelCounter::default();
    let config = Config::new(&args);
    let accumulator = ParallelAccumulator::default();
    let empty_scores = HashMap::new();

    tracks.into_par_iter().for_each(|(chr, components)| {
        counter.inc_comp(components.len() as u32);

        let local_scores = aparent_scores.get(chr.as_bytes()).unwrap_or_else(|| {
            log::debug!(
                "DEBUG: No APARENT scores found for chromosome {:?}; defaulting to empty map",
                chr
            );
            &empty_scores
        });

        process_components(components, local_scores, &counter, &config, &accumulator);
    });

    log::info!(
        "Detected potential intraprimings: {}",
        counter.load_intrapriming()
    );

    if args.recover {
        log::warn!(
            "Number of dirty components in query reads: {:?} ({:.3}%)",
            counter.num_of_reviews,
            counter.load_ratio()
        );
    }

    let file = args
        .outdir
        .join(format!("{}.pas.tsv", args.prefix.display()));
    let mut writer = BufWriter::new(File::create(file).unwrap());

    accumulator.lines.into_iter().for_each(|line| {
        writer
            .write_all(&line)
            .unwrap_or_else(|e| panic!("ERROR: Could not write schema -> {e}!"));
    });

    Ok(())
}

/// Config struct for the parallel processing
///
/// # Fields
///
/// * `aparent_threshold` - APARENT threshold [min score allowed]
/// * `filter` - Flag to filter components based on polyA tail length and APARENT score
/// * `filter_side` - Side to filter components based on polyA tail length and APARENT score
/// * `max_gpa_length` - Genomic polyA tail length threshold [max length allowed]
/// * `min_polya_length` - PolyA tail length threshold [min length allowed]
/// * `recover` - Flag to recover from disputed components where discard ratio is bigger than threshold
/// * `wiggle` - Wiggle room for polyA tail length
///
/// # Example
///
/// ```ignore
/// let config = Config::new(args);
///
/// assert_eq!(config.aparent_threshold, 0.01);
/// ```
#[derive(Debug, Clone)]
struct Config {
    aparent_threshold: f32,
    max_gpa_length: usize,
    min_polya_length: usize,
    recover: bool,
    wiggle: usize,
}

/// Constructor for the Config struct
///
/// # Arguments
///
/// * `args` - Args struct from the CLI
///
/// # Returns
///
/// * `Config` - Config struct with the given arguments
///
/// # Example
///
/// ```ignore
/// let config = Config::new(args);
/// ```
impl Config {
    fn new(args: &Args) -> Self {
        Self {
            aparent_threshold: args.aparent_threshold,
            max_gpa_length: args.max_gpa_length,
            min_polya_length: args.min_polya_length,
            recover: args.recover,
            wiggle: args.wiggle,
        }
    }
}

/// Distributes the components to the processing function
///
/// # Arguments
///
/// * `counter` - Counter for the parallel processing
/// * `components` - Vector of components to process
/// * `scores` - Nested map with chromosome as key and a map of read name and aparent score as value
/// * `accumulator` - Accumulator for the parallel processing
/// * `recover` - Flag to recover from disputed components where discard ratio is bigger than threshold
/// * `wiggle` - Wiggle room for polyA tail length
/// * `max_gpa_length` - Genomic polyA tail length threshold [max length allowed]
/// * `min_polya_length` - PolyA tail length threshold [min length allowed]
/// * `aparent_threshold` - APARENT threshold [min score allowed]
/// * `filter` - Flag to filter components based on polyA tail length and APARENT score
/// * `filter_side` - Side to filter components based on polyA tail length and APARENT score
///
/// # Example
///
/// ```ignore
/// distribute(args);
/// ```
#[inline(always)]
fn process_components(
    components: Vec<Vec<GenePred>>,
    scores: &HashMap<usize, f32>,
    counter: &ParallelCounter,
    config: &Config,
    accumulator: &ParallelAccumulator,
) {
    components.into_par_iter().for_each(|comp| {
        let (refs, queries) = comp.into_iter().partition(|record| {
            let role = record
                .get_extra(b"role")
                .unwrap_or_else(|| panic!("ERROR: Could not get role from record!"))
                .clone()
                .into_scalar();

            role == Some(b"reference".to_vec())
        });

        let result = process_component(refs, queries, scores, config, counter);

        result.into_iter().for_each(|descriptor| {
            if !descriptor.is_empty() {
                accumulator.lines.insert(descriptor);
            }
        });
    });
}

/// Parallel counter for the processing function
///
/// # Fields
///
/// * `num_of_comps` - AtomicU32 to store the number of components
/// * `num_of_dirty` - AtomicU32 to store the number of dirty components
///
/// # Example
///
/// ```ignore
/// let counter = ParallelCounter::default();
///
/// assert_eq!(counter.num_of_comps.load(Ordering::Relaxed), 0);
/// ```
struct ParallelCounter {
    num_of_comps: AtomicU32,
    num_of_reviews: AtomicU32,
    num_of_intrapriming: AtomicU32,
}

impl Default for ParallelCounter {
    fn default() -> Self {
        Self {
            num_of_comps: AtomicU32::new(0),
            num_of_reviews: AtomicU32::new(0),
            num_of_intrapriming: AtomicU32::new(0),
        }
    }
}

impl ParallelCounter {
    /// Increment the component counter
    ///
    /// # Arguments
    ///
    /// * `value` - Value to increment the counter
    ///
    /// # Example
    ///
    /// ```ignore
    /// let counter = ParallelCounter::default();
    /// counter.inc_comp(1);
    ///
    /// assert_eq!(counter.num_of_comps.load(Ordering::Relaxed), 1);
    /// ```
    fn inc_comp(&self, value: u32) {
        self.num_of_comps.fetch_add(value, Ordering::Relaxed);
    }

    /// Increment the intrapriming counter
    ///
    /// # Arguments
    ///
    /// * `value` - Value to increment the counter
    ///
    /// # Example
    ///
    /// ```ignore
    /// let counter = ParallelCounter::default();
    /// counter.inc_intrapriming();
    ///
    /// assert_eq!(counter.num_of_intrapriming.load(Ordering::Relaxed), 1);
    /// ```
    fn inc_intrapriming(&self) {
        self.num_of_intrapriming.fetch_add(1, Ordering::Relaxed);
    }

    /// Increment the review counter
    ///
    /// # Arguments
    ///
    /// * `value` - Value to increment the counter
    ///
    /// # Example
    ///
    /// ```ignore
    /// let counter = ParallelCounter::default();
    /// counter.inc_review(1);
    ///
    /// assert_eq!(counter.num_of_reviews.load(Ordering::Relaxed), 1);
    /// ```
    fn inc_review(&self, value: u32) {
        self.num_of_reviews.fetch_add(value, Ordering::Relaxed);
    }

    /// Load the ratio of reviews to components
    ///
    /// # Example
    ///
    /// ```ignore
    /// let counter = ParallelCounter::default();
    /// counter.inc_comp(10);
    /// counter.inc_review(2);
    ///
    /// assert_eq!(counter.load_ratio(), 0.2);
    /// ```
    fn load_ratio(&self) -> f64 {
        self.num_of_reviews.load(Ordering::Relaxed) as f64
            / self.num_of_comps.load(Ordering::Relaxed) as f64
    }

    /// Get the number of intraprimming events
    fn load_intrapriming(&self) -> u32 {
        self.num_of_intrapriming.load(Ordering::Relaxed)
    }
}

/// Macro to check if any of the conditions are true
///
/// # Example
///
/// ```ignore
/// let x = any_true!(false, false, true);
///
/// assert_eq!(x, true);
/// ```
macro_rules! any_true {
    ($($x:expr),*) => {
        false $(|| $x)*
    };
}

/// Process Vec<PolyApred> to detect intrapriming events
///
/// # Arguments
///
/// * `component` - Vector of PolyAPred to process
/// * `scores` - Nested map with chromosome as key and a map of read name and aparent score as value
/// * 'recover' - Flag to recover from disputed components where discard ratio is bigger than threshold
/// * `wiggle` - Wiggle room for polyA tail length
/// * `max_gpa_length` - Genomic polyA tail length threshold [max length allowed]
/// * `min_polya_length` - PolyA tail length threshold [min length allowed]
/// * `aparent_threshold` - APARENT threshold [min score allowed]
///
/// # Returns
///
/// * `Vec<&String>` - Vector of pass reads
/// * `Vec<&String>` - Vector of intrapriming reads
/// * `Option<Vec<&String>>` - Vector of review reads
fn process_component(
    refs: Vec<GenePred>,
    queries: Vec<GenePred>,
    scores: &HashMap<usize, f32>,
    config: &Config,
    counter: &ParallelCounter,
) -> Vec<Vec<u8>> {
    if queries.is_empty() {
        // log::debug!("DEBUG: No queries found, skipping component");
        return Vec::new();
    }

    let mut descriptor = HashMap::new();
    let mut accumulator = Vec::new();

    let mut seen_positions = HashSet::new();
    let mut seen_reads = HashSet::new();

    let (mut count, totals) = (0_f32, queries.len() as f32);

    queries.iter().for_each(|read| {
        let end = match read.strand() {
            Some(Strand::Forward) => read.end(),
            Some(Strand::Reverse) => read.start(),
            Some(Strand::Unknown) | None => {
                log::error!("ERROR: Could not get strand from read {:?}", read.name());
                std::process::exit(1);
            }
        };

        let mut schema = Schema::default();
        let score = scores.get(&(end as usize)).copied().unwrap_or_else(|| {
            // INFO: we cannot panic because peak_threshold in APARENT side drops
            // very low scores
            log::warn!(
                "WARN: Could not get aparent score for read {:?} in scores",
                read,
            );
            0.0
        });

        let (clip, clipped_a, poly_a, gpa) = get_polya_stats(
            read.name()
                .unwrap_or_else(|| panic!("ERROR: Could not get read name {:?}", read.name())),
        );

        let location = get_location(read.end, &refs);
        schema.position = location;

        if any_true!(
            gpa < config.max_gpa_length as u32,
            poly_a >= config.min_polya_length as u32,
            seen_positions.contains(&end),
            seen_positions.contains(&(end - config.wiggle as u64)),
            seen_positions.contains(&(end + config.wiggle as u64)),
            score > config.aparent_threshold,
            location == Position::UTR
        ) {
            // passes.push(read.line.clone());
            schema.status = b"PASS";
            schema.code = b"A";
            schema.support = b"POLYA_SUPPORTED";
            schema.forced = b"NOT_FORCED";

            seen_positions.insert(end); // INFO: only considering 'good' reads
            seen_positions.insert(end - config.wiggle as u64); // INFO: inserting wiggle room!
            seen_positions.insert(end + config.wiggle as u64);

            seen_reads.insert(&read.name); // INFO: storing 'good' reads for later escape
        }

        schema.id = read.name().unwrap_or_else(|| {
            panic!("ERROR: Read has no name! -> {:?}", read.name());
        });
        schema.score = score;
        schema.position = location;
        schema.clip = clip;
        schema.clipped_a = clipped_a;
        schema.poly_a = poly_a;
        schema.gpa = gpa;
        schema.size = queries.len();

        descriptor.insert(read.name(), schema);
    });

    // INFO: necessary second iteration to check for intrapriming events
    // INFO: since the only filter to doble check is the wiggle + genomic pos
    // INFO: we can iterate over the component again and avoid the rest of steps
    queries.iter().for_each(|read| {
        let end = match read.strand() {
            Some(Strand::Forward) => read.end(),
            Some(Strand::Reverse) => read.start(),
            Some(Strand::Unknown) | None => {
                panic!("ERROR: Could not get strand from read {:?}", read.name())
            }
        };

        let schema = descriptor.get_mut(&read.name()).unwrap_or_else(|| {
            panic!(
                "ERROR: Could not get handle for read {:?}",
                read.name()
                    .unwrap_or_else(|| panic!("ERROR: Could not get read name {:?}", read.name()))
            )
        });

        if seen_reads.contains(&&read.name) {
            log::debug!("DEBUG: Read {:?} classified as PASS in first iteration, skipping intrapriming check -> {:?}", read.name(), schema);
            return;
        }

        if any_true!(
            seen_positions.contains(&end),
            seen_positions.contains(&(end - config.wiggle as u64)),
            seen_positions.contains(&(end + config.wiggle as u64))
        ) {
            seen_positions.insert(end); // INFO: only considering 'good' reads
            seen_positions.insert(end - config.wiggle as u64); // INFO: inserting wiggle room!
            seen_positions.insert(end + config.wiggle as u64);

            schema.status = b"PASS";
            schema.code = b"AF";

            // INFO: set to false because the evidence is based on other reads!
            schema.support = b"POLYA_NOT_SUPPORTED";
            schema.forced = b"FORCED";
        } else {
            counter.inc_intrapriming();
            count += 1.0;

            schema.status = b"INTRAPRIMING";
            schema.support = b"POLYA_NOT_SUPPORTED";
            schema.forced = b"NOT_FORCED";
            schema.code = b"I";
        }
    });

    let ratio = count / totals;

    // INFO: third pass to fill out ratios and evaluate recovery
    queries.iter().for_each(|read| {
        let schema = descriptor.get_mut(&read.name()).unwrap_or_else(|| {
            panic!(
                "ERROR: Could not get handle for read {:?}",
                read.name()
                    .unwrap_or_else(|| panic!("ERROR: Could not get read name {:?}", read.name()))
            )
        });

        schema.ratio = ratio;
        schema.events = count as usize;

        if config.recover && ratio > INTRAPRIMING_RATIO_THRESHOLD {
            counter.inc_review(1);
            schema.status = b"REVIEW";
        }

        accumulator.push(schema.to_line());
    });

    accumulator
}

/// Struct to hold the schema for the PASCaller
#[derive(Debug, Default, Clone)]
struct Schema<'a> {
    id: &'a [u8],
    status: &'a [u8],
    code: &'a [u8],
    support: &'a [u8],
    forced: &'a [u8],
    clip: u32,
    clipped_a: u32,
    poly_a: u32,
    gpa: u32,
    position: Position,
    score: f32,
    ratio: f32,
    size: usize,
    events: usize,
}

impl Schema<'_> {
    pub fn to_line(&self) -> Vec<u8> {
        let mut body = String::new();

        body.push_str("<h2>Intraprimming</h2><br>");
        body.push_str(&format!(
            "- status: {}<br>",
            std::str::from_utf8(&self.status).unwrap_or("NULL")
        ));
        body.push_str(&format!(
            "- code: {} [I: intrapriming, A: no events, F: forced]<br>",
            std::str::from_utf8(&self.code).unwrap_or("NULL")
        ));
        body.push_str(&format!(
            "- support: {}<br>",
            std::str::from_utf8(&self.support).unwrap_or("NULL")
        ));
        body.push_str(&format!(
            "- forced: {}<br>",
            std::str::from_utf8(&self.forced).unwrap_or("NULL")
        ));

        body.push_str(&format!("- component size: {}<br>", self.size));
        body.push_str(&format!("- component events: {}<br>", self.events));
        body.push_str(&format!("- ratio: {}<br>", self.ratio));
        body.push_str(&format!("- APARENT score: {}<br>", self.score));
        body.push_str(&format!("- GPA: {}<br>", self.gpa));
        body.push_str(&format!("- polyA: {}<br>", self.poly_a));
        body.push_str(&format!("- clipped A: {}<br>", self.clipped_a));
        body.push_str(&format!("- clip: {}<br>", self.clip));
        body.push_str(&format!("- position: {}<br>", self.position));

        // fmt: id\tstatus\tcode\thtml  — all on one line
        let mut out = Vec::new();
        out.extend_from_slice(&self.id);
        out.push(b'\t');
        out.extend_from_slice(&self.status);
        out.push(b'\t');
        out.extend_from_slice(&self.code);
        out.push(b'\t');
        out.extend_from_slice(body.as_bytes());
        out.push(b'\n'); // single trailing newline for the TSV row

        out
    }
}

/// Enum to represent the location of the polyA tail
#[derive(Debug, PartialEq, Clone, Copy, Default)]
enum Position {
    CDS,
    UTR,
    #[default]
    Unknown,
}

impl std::fmt::Display for Position {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Position::CDS => write!(f, "CDS"),
            Position::UTR => write!(f, "UTR"),
            Position::Unknown => write!(f, "UNKNOWN"),
        }
    }
}

/// Get the read location using TOGA
///
/// # Arguments
///
/// * `read_end` - End position of the read
/// * `toga` - Vector of GenePred to process
/// * `handle` - ModuleMap to set the value
///
/// # Example
///
/// ```ignore
/// let read_end = 100;
/// let mut toga = vec![GenePred::default()];
/// let mut handle = ModuleMap::default();
///
/// get_location(read_end, &mut toga, &mut handle);
/// ```
fn get_location(read_end: u64, refs: &Vec<GenePred>) -> Position {
    let location = Position::Unknown;

    for projection in refs.iter() {
        if projection.start < read_end && read_end < projection.end {
            let utr = projection.three_prime_utr();
            if utr
                .iter()
                .any(|(utr_start, utr_end)| *utr_start < read_end && read_end < *utr_end)
            {
                return Position::UTR;
            }

            return Position::CDS;
        }
    }

    location
}

/// Parallel accumulator
struct ParallelAccumulator {
    lines: DashSet<Vec<u8>>,
}

/// Default implementation for ParallelAccumulator
impl Default for ParallelAccumulator {
    fn default() -> Self {
        Self {
            lines: DashSet::new(),
        }
    }
}

/// Loads aparent scores from a bigWig fieles.
///
/// # Arguments
///
/// * `forward`: Path to the forward bigWig file.
/// * `reverse`: Path to the reverse bigWig file.
/// * `chrs`: A `Vec<String>` of chromosome names to process.
///
/// # Example
///
/// ```rust,ignore
/// let aparent_scores = load_aparent(PathBuf::from("aparent.bw"));
/// ```
pub fn load_aparent<T: AsRef<std::path::Path> + std::fmt::Debug + Sync>(
    forward: T,
    reverse: T,
    chrs: Vec<Vec<u8>>,
) -> HashMap<Vec<u8>, HashMap<usize, f32>> {
    log::info!("INFO: Parsing BigWigs...");
    log::info!("INFO: BigWig files: {forward:?}, {reverse:?}");
    let scores = bigwig_to_map(vec![forward, reverse], &chrs);
    scores
}

/// Converts a vector of BigWig files into a vector of thread-safe maps.
///
/// This function is designed to be run in parallel for plus and minus strands. It iterates through
/// the BigWig files, reads the chromosome data, and populates `DashMap`s with scores that
/// meet a certain significance threshold.
///
/// # Arguments
///
/// * `bigwigs`: A `Vec` of paths to the BigWig files (e.g., donor and acceptor).
/// * `chrs`: A slice of `String` representing the chromosomes to be processed.
///
/// # Returns
///
/// * A `Vec<DashMap<String, DashMap<usize, f32>>>` where the outer vector corresponds to donor/acceptor
///   sites, the middle map keys are chromosome names, and the inner map keys are genomic positions with their scores.
///
/// # Example
///
/// ```rust,ignore
/// let splice_maps = bigwig_to_map(bigwigs, &chrs);
/// ```
fn bigwig_to_map<T: AsRef<std::path::Path> + std::fmt::Debug + Sized + Sync>(
    bigwigs: Vec<T>,
    chrs: &Vec<Vec<u8>>,
) -> HashMap<Vec<u8>, HashMap<usize, f32>> {
    let total_count = AtomicU32::new(0);
    let accumulator = Mutex::new(HashMap::new());

    log::debug!(
        "Will extract scores from the following chromosomes: {:?}",
        chrs.iter()
            .map(|chr| std::str::from_utf8(chr).unwrap())
            .collect::<Vec<_>>()
    );

    // [donor, acceptor]
    bigwigs
        .into_par_iter()
        .zip(vec![genepred::Strand::Forward, genepred::Strand::Reverse])
        .for_each(|(bigwig, strand)| {
            log::debug!(
                "Will process bigwig file: {:?} in strand {:?}",
                bigwig,
                strand
            );
            let acc = DashMap::new();

            let bwread = BigWigRead::open_file(bigwig).expect("ERROR: Cannot open BigWig file");
            let chroms: Vec<_> = bwread.chroms().to_vec();

            chroms.into_par_iter().for_each(|chr| {
                log::debug!("Will process chromosome {:?}", chr.name);
                let mut bwread =
                    BigWigRead::reopen(&bwread).expect("ERROR: Cannot re-open BigWig file");

                if !chrs.contains(&chr.name.as_bytes().to_vec()) {
                    log::debug!("Skipping chromosome {:?}", chr.name);
                    return; // INFO: skip chromosomes not in records
                }

                let name = format!("{}:{}", chr.name, strand);
                let length = chr.length;
                let values = bwread
                    .get_interval(&chr.name, 0, length)
                    .expect("ERROR: Cannot read intervals from BigWig!");

                let (positions, local_count) = collect_bigwig_values(
                    values.map(|value| value.expect("ERROR: Cannot read values from BigWig!")),
                );

                if positions.is_empty() {
                    log::warn!("WARN: No bigwig values found for {:?}", chr.name);
                    return;
                }

                acc.insert(name.as_bytes().to_vec(), positions);
                total_count.fetch_add(local_count, Ordering::Relaxed);
            });

            // INFO: lock main hash and diffuse local dash
            let mut guard = accumulator.lock().expect("ERROR: Cannot lock mutex");
            acc.into_iter().for_each(|(k, v)| {
                guard.entry(k).or_insert(v);
            });
        });

    log::info!(
        "INFO: Parsed and combined {} significant splicing scores from BigWigs!",
        total_count.load(Ordering::Relaxed)
    );

    accumulator
        .into_inner()
        .unwrap_or_else(|e| panic!("ERROR: Could not unwrap accumulator -> {e}!"))
}

fn collect_bigwig_values<I>(values: I) -> (HashMap<usize, f32>, u32)
where
    I: IntoIterator<Item = Value>,
{
    let mut positions = HashMap::new();
    let mut count = 0;

    for value in values {
        for pos in value.start as usize..value.end as usize {
            positions.entry(pos).or_insert(value.value);
            count += 1;
        }
    }

    (positions, count)
}

/// Extracts the polyA stats from the read name
///
/// # Arguments
///
/// * `read` - A string slice that holds the read name
///
/// # Returns
///
/// A tuple with the following values:
///
/// * `clip` - The number of clipped bases
/// * `clipped_a` - The number of clipped A's
/// * `read_a` - The number of A's in the read
/// * `gpa` - The number of genomic A's
///
/// # Example
///
/// ```ignore
/// let read = "m54164U_210309_085211/74646562/ccs_PerID0.995_5Clip0_3Clip0_PolyA49_PolyARead50";
/// let stats = get_polya_stats(read);
///
/// assert_eq!(stats, (0, 49, 50, 1));
/// ```
fn get_polya_stats(read: &[u8]) -> (u32, u32, u32, u32) {
    let tags = get_tags(read);

    let clip3 = *tags.get(b"TC").unwrap() as u32;
    let clipped_a = *tags.get(b"PA").unwrap() as u32;
    let read_a = *tags.get(b"PR").unwrap() as u32;

    // INFO: the whole polyA is clipped!
    if clipped_a == read_a {
        return (clip3, clipped_a, read_a, 0);
    } else if clipped_a > read_a {
        log::error!("ERROR: clipped A's is greater than read A's");
        std::process::exit(1);
    }

    let gpa = read_a - (clip3 + clipped_a);

    (clip3, clipped_a, read_a, gpa)
}

/// Extracts two-letter key, usize value tags from a read's name field.
///
/// This function is designed to parse a specific tag format, where tags are separated by a
/// `BIG_SEP` and then a `SEP` character, with each tag being two characters followed by a
/// numeric value.
///
/// # Arguments
///
/// * `read` - A string slice containing the full read name, including tags.
///
/// # Returns
///
/// A `HashMap<String, usize>` containing the parsed tags.
///
/// # Panics
///
/// This function assumes the presence of `BIG_SEP` and `SEP` constants and that tag keys are
/// exactly two characters long. The capacity is enforced at 5.
///
fn get_tags(read: &[u8]) -> HashMap<[u8; 2], usize> {
    let mut map = HashMap::with_capacity(10); // WARN: enforcing 5 tags -> fusion tags come after this step!
    let read = unsafe { from_utf8_unchecked(read) };

    if let Some((_, tags)) = read.split_once(BIG_SEP) {
        // WARN: enforcing two letter tag!
        for tag in tags.split(SEP) {
            if tag.len() >= 3 {
                // Safe slicing because ASCII: two-letter key + numeric value
                let key: [u8; 2] = [tag.as_bytes()[0], tag.as_bytes()[1]];
                let val = &tag[2..];
                if let Ok(parsed_val) = val.parse::<usize>() {
                    map.insert(key, parsed_val);
                }
            }
        }
    }

    map
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn collect_bigwig_values_only_materializes_covered_bases() {
        let (positions, count) = collect_bigwig_values(vec![
            Value {
                start: 2,
                end: 4,
                value: 0.5,
            },
            Value {
                start: 7,
                end: 8,
                value: 0.8,
            },
        ]);

        assert_eq!(count, 3);
        assert_eq!(positions.len(), 3);
        assert_eq!(positions.get(&2), Some(&0.5));
        assert_eq!(positions.get(&3), Some(&0.5));
        assert_eq!(positions.get(&7), Some(&0.8));
        assert!(!positions.contains_key(&0));
        assert!(!positions.contains_key(&6));
    }

    #[test]
    fn get_location_detects_three_prime_utr_without_panic() {
        let refs = vec![GenePred {
            chrom: b"chr1".to_vec(),
            start: 100,
            end: 200,
            name: Some(b"tx1".to_vec()),
            strand: Some(Strand::Forward),
            thick_start: Some(120),
            thick_end: Some(150),
            block_count: Some(1),
            block_starts: Some(vec![100]),
            block_ends: Some(vec![200]),
            extras: std::collections::HashMap::new(),
        }];

        assert_eq!(get_location(175, &refs), Position::UTR);
        assert_eq!(get_location(130, &refs), Position::CDS);
    }
}
