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

    if args.genome.is_none() {
        log::warn!(
            "WARN: No --genome given; the genomic A-run (GPA) falls back to the read-name tags [PR - TC - PA]"
        );
    }

    tracks.into_par_iter().for_each(|(chr, components)| {
        counter.inc_comp(components.len() as u32);

        let local_scores = aparent_scores.get(chr.as_bytes()).unwrap_or_else(|| {
            log::debug!(
                "DEBUG: No APARENT scores found for chromosome {:?}; defaulting to empty map",
                chr
            );
            &empty_scores
        });

        // INFO: one chromosome sequence per task; keys are '<chr>:<strand>'
        let sequence = args
            .genome
            .as_ref()
            .and_then(|genome| load_chromosome(genome, chr.split(':').next().unwrap_or(&chr)));

        process_components(
            components,
            local_scores,
            sequence.as_deref(),
            &counter,
            &config,
            &accumulator,
        );
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
    aparent_window: usize,
    gpa_window: usize,
    max_gpa_length: usize,
    min_polya_length: usize,
    recover: bool,
    tes_window: usize,
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
            aparent_window: args.aparent_window,
            gpa_window: args.gpa_window,
            max_gpa_length: args.max_gpa_length,
            min_polya_length: args.min_polya_length,
            recover: args.recover,
            tes_window: args.tes_window,
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
    sequence: Option<&[u8]>,
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

        let result = process_component(refs, queries, scores, sequence, config, counter);

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
    sequence: Option<&[u8]>,
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

    // INFO: a `#CN<n>` record stands for n identical reads
    let size = queries.iter().map(cn_weight).sum::<usize>();
    let (mut count, totals) = (0_f32, size as f32);

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
        // INFO: max over a window around the 3' end, not a point lookup: a cleavage
        // site a few nt away from the model end otherwise reads as 0. Missing positions
        // are expected (peak_threshold in APARENT drops very low scores)
        let score = get_aparent_score(scores, end, config.aparent_window);

        let (clip, clipped_a, poly_a, tag_gpa) = get_polya_stats(
            read.name()
                .unwrap_or_else(|| panic!("ERROR: Could not get read name {:?}", read.name())),
        );

        // INFO: with a genome, the genomic A-run is read downstream of the model's
        // 3' end, not taken from the tags of the one read that named the model
        let gpa = match sequence {
            Some(seq) => get_genomic_a_run(seq, end, read.strand, config.gpa_window),
            None => tag_gpa,
        };

        // WARN: strand-aware 3' end; `read.end` is the 5' end on the minus strand
        let location = get_location(end, read.strand, &refs, config.tes_window);
        schema.position = location;

        if any_true!(
            gpa < config.max_gpa_length as u32,
            poly_a >= config.min_polya_length as u32,
            seen_positions.contains(&end),
            seen_positions.contains(&(end - config.wiggle as u64)),
            seen_positions.contains(&(end + config.wiggle as u64)),
            score > config.aparent_threshold,
            matches!(location, Position::UTR | Position::TES)
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
        schema.gpa_source = if sequence.is_some() { "genome" } else { "tags" };
        schema.size = size;

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
            count += cn_weight(read) as f32;

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
    gpa_source: &'static str,
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
        body.push_str(&format!("- GPA: {} [{}]<br>", self.gpa, self.gpa_source));
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
    TES, // INFO: within --tes-window of an annotated transcript end
    #[default]
    Unknown,
}

impl std::fmt::Display for Position {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Position::CDS => write!(f, "CDS"),
            Position::UTR => write!(f, "UTR"),
            Position::TES => write!(f, "TES"),
            Position::Unknown => write!(f, "UNKNOWN"),
        }
    }
}

/// Get the location of a strand-aware 3' end against the references
///
/// Priority: inside an annotated 3' UTR (UTR), else within `tes_window` nt of an
/// annotated transcript end (TES), else inside a reference body (CDS), else UNKNOWN.
/// The TES test is what rescues annotated ends when references carry no UTRs
/// (TOGA2 projections: thickEnd == txEnd), where the UTR test cannot fire.
///
/// # Arguments
///
/// * `three_prime` - 3' end of the model [end on '+', start on '-']
/// * `strand` - Strand of the model
/// * `refs` - Vector of reference GenePred in the component
/// * `tes_window` - Max distance to an annotated transcript end
///
/// # Example
///
/// ```ignore
/// let location = get_location(175, Some(Strand::Forward), &refs, 50);
/// ```
fn get_location(
    three_prime: u64,
    strand: Option<Strand>,
    refs: &[GenePred],
    tes_window: usize,
) -> Position {
    let mut location = Position::Unknown;

    // WARN: intervals are half-open [start, end); `three_prime` is the exclusive end on
    // '+' and the inclusive start on '-'. Test the last transcribed base, so that a model
    // ending exactly at an annotated txEnd/UTR end is inside it [the old strict `<` on
    // both sides left such ends UNKNOWN with --tes-window 0]
    let base = match strand {
        Some(Strand::Reverse) => three_prime,
        _ => three_prime.saturating_sub(1),
    };

    for projection in refs.iter() {
        if projection.start <= base && base < projection.end {
            let utr = projection.three_prime_utr();
            if utr
                .iter()
                .any(|(utr_start, utr_end)| *utr_start <= base && base < *utr_end)
            {
                return Position::UTR;
            }

            location = Position::CDS;
        }
    }

    let near_tes = refs.iter().any(|projection| {
        let tes = match projection.strand {
            Some(Strand::Reverse) => projection.start,
            _ => projection.end,
        };
        tes.abs_diff(three_prime) <= tes_window as u64
    });

    if near_tes {
        return Position::TES;
    }

    location
}

/// Max APARENT score over `[end - window, end + window]`; 0.0 when no peak is there
///
/// # Example
///
/// ```ignore
/// let score = get_aparent_score(&scores, 1000, 25);
/// ```
fn get_aparent_score(scores: &HashMap<usize, f32>, end: u64, window: usize) -> f32 {
    let end = end as usize;
    (end.saturating_sub(window)..=end.saturating_add(window))
        .filter_map(|pos| scores.get(&pos).copied())
        .fold(0.0, f32::max)
}

/// Longest genomic A-run downstream of a strand-aware 3' end
///
/// On '+' scans `seq[end..end + window]` for A; on '-' scans
/// `seq[start - window..start]` for T (A on the transcript strand).
///
/// # Example
///
/// ```ignore
/// let gpa = get_genomic_a_run(b"CCCCAAAAAAG", 4, Some(Strand::Forward), 20);
/// assert_eq!(gpa, 6);
/// ```
fn get_genomic_a_run(seq: &[u8], end: u64, strand: Option<Strand>, window: usize) -> u32 {
    let end = (end as usize).min(seq.len());
    let (slice, base) = match strand {
        Some(Strand::Reverse) => (&seq[end.saturating_sub(window)..end], b'T'),
        _ => (&seq[end..(end + window).min(seq.len())], b'A'),
    };

    let (mut best, mut run) = (0, 0);
    for nt in slice {
        if nt.to_ascii_uppercase() == base {
            run += 1;
            best = best.max(run);
        } else {
            run = 0;
        }
    }

    best
}

/// Loads one chromosome sequence from a .2bit file; None when it is absent
fn load_chromosome(genome: &std::path::Path, chr: &str) -> Option<Vec<u8>> {
    let mut reader = twobit::TwoBitFile::open(genome)
        .unwrap_or_else(|e| panic!("ERROR: Cannot open 2bit file {:?} -> {e}", genome));

    match reader.read_sequence(chr, ..) {
        Ok(seq) => Some(seq.into_bytes()),
        Err(e) => {
            log::warn!(
                "WARN: No sequence for {chr} in {:?} -> {e}; using read tags",
                genome
            );
            None
        }
    }
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

/// Support of a read: the `#CN<n>` tag of its name, 1 when absent.
fn cn_weight(read: &GenePred) -> usize {
    let tags = get_tags(read.name().unwrap_or_default());
    tags.get(b"CN").copied().unwrap_or(1)
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

        let fwd = Some(Strand::Forward);
        assert_eq!(get_location(175, fwd, &refs, 0), Position::UTR);
        assert_eq!(get_location(130, fwd, &refs, 0), Position::CDS);
        assert_eq!(get_location(130, fwd, &refs, 50), Position::CDS);
    }

    #[test]
    fn end_exactly_at_the_annotated_tx_end_is_inside_the_utr() {
        // half-open coordinates: a '+' model ending at txEnd (200) has its last base (199)
        // in the 3' UTR [150, 200); a '-' model starting at txStart (100) has its first
        // base in the 3' UTR [100, 120). Both read UNKNOWN with the old strict `<` tests
        let utr_ref = |strand| GenePred {
            chrom: b"chr1".to_vec(),
            start: 100,
            end: 200,
            name: Some(b"tx1".to_vec()),
            strand: Some(strand),
            thick_start: Some(120),
            thick_end: Some(150),
            block_count: Some(1),
            block_starts: Some(vec![100]),
            block_ends: Some(vec![200]),
            extras: std::collections::HashMap::new(),
        };

        let fwd = vec![utr_ref(Strand::Forward)];
        assert_eq!(
            get_location(200, Some(Strand::Forward), &fwd, 0),
            Position::UTR
        );
        // one past txEnd is outside the record
        assert_eq!(
            get_location(201, Some(Strand::Forward), &fwd, 0),
            Position::Unknown
        );
        // first base of the UTR [end = 151 -> base 150] is UTR, the last CDS base is not
        assert_eq!(
            get_location(151, Some(Strand::Forward), &fwd, 0),
            Position::UTR
        );
        assert_eq!(
            get_location(150, Some(Strand::Forward), &fwd, 0),
            Position::CDS
        );

        let rev = vec![utr_ref(Strand::Reverse)];
        assert_eq!(
            get_location(100, Some(Strand::Reverse), &rev, 0),
            Position::UTR
        );
        assert_eq!(
            get_location(119, Some(Strand::Reverse), &rev, 0),
            Position::UTR
        );
        assert_eq!(
            get_location(120, Some(Strand::Reverse), &rev, 0),
            Position::CDS
        );
        assert_eq!(
            get_location(99, Some(Strand::Reverse), &rev, 0),
            Position::Unknown
        );
    }

    fn test_config() -> Config {
        Config {
            aparent_threshold: 0.01,
            aparent_window: APARENT_WINDOW,
            gpa_window: GPA_WINDOW,
            max_gpa_length: 5,
            min_polya_length: 50,
            recover: false,
            tes_window: TES_WINDOW,
            wiggle: 2,
        }
    }

    fn model(name: &str, start: u64, end: u64, strand: Strand) -> GenePred {
        let mut read = GenePred::from_coords(b"chr1".to_vec(), start, end, Default::default());
        read.set_name(Some(name.as_bytes().to_vec()));
        read.set_strand(Some(strand));
        read
    }

    fn status(
        refs: Vec<GenePred>,
        query: GenePred,
        scores: &HashMap<usize, f32>,
        sequence: Option<&[u8]>,
        config: &Config,
    ) -> String {
        let counter = ParallelCounter::default();
        let line = process_component(refs, vec![query], scores, sequence, config, &counter);
        String::from_utf8(line[0].clone()).unwrap()
    }

    // Q4 3a: AP2D-like model, GPA 21 from the tags, APARENT ~0, UNKNOWN position
    const IP_NAME: &str = "R16@53__TC0#PA0#PR21#CN33";

    #[test]
    fn minus_strand_end_inside_three_prime_utr_is_rescued() {
        // '-' reference: 3' UTR is [1000, 1200), CDS [1200, 1800)
        let mut reference = model("ref", 1000, 2000, Strand::Reverse);
        reference.thick_start = Some(1200);
        reference.thick_end = Some(1800);
        // the model's 3' end (start = 1100) is in the UTR; its 5' end (2500) is not
        let query = model(IP_NAME, 1100, 2500, Strand::Reverse);
        let config = Config {
            tes_window: 0,
            ..test_config()
        };

        let line = status(vec![reference], query, &HashMap::new(), None, &config);
        assert!(line.contains("\tPASS\t"), "{line}");
        assert!(line.contains("position: UTR"), "{line}");
    }

    #[test]
    fn reference_without_utr_rescues_ends_at_the_annotated_tes() {
        // TOGA2-like '+' projection: thickEnd == txEnd, no 3' UTR to test
        let mut reference = model("ref", 1000, 2000, Strand::Forward);
        reference.thick_start = Some(1000);
        reference.thick_end = Some(2000);
        let near = model(IP_NAME, 1000, 2030, Strand::Forward);
        let far = model(IP_NAME, 1000, 2300, Strand::Forward);

        let line = status(
            vec![reference.clone()],
            near,
            &HashMap::new(),
            None,
            &test_config(),
        );
        assert!(
            line.contains("\tPASS\t") && line.contains("position: TES"),
            "{line}"
        );
        let line = status(vec![reference], far, &HashMap::new(), None, &test_config());
        assert!(line.contains("\tINTRAPRIMING\t"), "{line}");
    }

    #[test]
    fn tes_window_is_inclusive_and_strand_aware() {
        // '+' TES = txEnd (2000, exclusive), '-' TES = txStart (1000)
        let fwd = vec![model("ref", 1000, 2000, Strand::Forward)];
        let rev = vec![model("ref", 1000, 2000, Strand::Reverse)];
        let f = Some(Strand::Forward);
        let r = Some(Strand::Reverse);

        assert_eq!(get_location(2050, f, &fwd, 50), Position::TES);
        assert_eq!(get_location(2051, f, &fwd, 50), Position::Unknown);
        assert_eq!(get_location(950, r, &rev, 50), Position::TES);
        assert_eq!(get_location(949, r, &rev, 50), Position::Unknown);
        // a '-' end near txEnd (its 5' end) is not a TES
        assert_eq!(get_location(1990, r, &rev, 50), Position::CDS);
        // window 0: only the exact TES
        assert_eq!(get_location(2001, f, &fwd, 0), Position::Unknown);
        assert_eq!(get_location(1000, r, &rev, 0), Position::TES);
    }

    #[test]
    fn aparent_peak_a_few_nt_from_the_end_is_found() {
        let query = || model(IP_NAME, 1000, 2000, Strand::Forward);
        let scores: HashMap<usize, f32> = [(2003, 0.4)].into_iter().collect();

        let line = status(vec![], query(), &scores, None, &test_config());
        assert!(
            line.contains("\tPASS\t") && line.contains("APARENT score: 0.4"),
            "{line}"
        );
        // the old point lookup misses it
        let config = Config {
            aparent_window: 0,
            ..test_config()
        };
        let line = status(vec![], query(), &scores, None, &config);
        assert!(line.contains("\tINTRAPRIMING\t"), "{line}");
    }

    #[test]
    fn load_chromosome_reads_a_two_bit_sequence() {
        // minimal .2bit: one sequence "chr1" = CCCCAAAAAAGT, no N/mask blocks
        let seq = b"CCCCAAAAAAGT";
        let code = |nt: u8| match nt {
            b'T' => 0u8,
            b'C' => 1,
            b'A' => 2,
            _ => 3,
        };
        let mut bytes = Vec::new();
        for word in [0x1A412743_u32, 0, 1, 0] {
            bytes.extend(word.to_le_bytes());
        }
        bytes.push(4);
        bytes.extend(b"chr1");
        bytes.extend(((bytes.len() + 4) as u32).to_le_bytes());
        for word in [seq.len() as u32, 0, 0, 0] {
            bytes.extend(word.to_le_bytes());
        }
        bytes.extend(seq.chunks(4).map(|chunk| {
            chunk
                .iter()
                .enumerate()
                .fold(0u8, |acc, (i, nt)| acc | (code(*nt) << (6 - 2 * i)))
        }));

        let path = std::env::temp_dir().join(format!("iso_pas_{}.2bit", std::process::id()));
        std::fs::write(&path, bytes).unwrap();
        let chr1 = load_chromosome(&path, "chr1");
        let missing = load_chromosome(&path, "chr2");
        std::fs::remove_file(&path).unwrap();

        assert_eq!(chr1.as_deref(), Some(&seq[..]));
        assert_eq!(missing, None);
        assert_eq!(
            get_genomic_a_run(seq, 4, Some(Strand::Forward), GPA_WINDOW),
            6
        );
    }

    #[test]
    fn genomic_a_run_comes_from_the_genome_downstream_of_the_end() {
        // tags claim 21 genomic A; the genome downstream of the '+' end has none
        let mut genome = vec![b'C'; 3000];
        let query = || model(IP_NAME, 1000, 2000, Strand::Forward);
        let line = status(
            vec![],
            query(),
            &HashMap::new(),
            Some(&genome),
            &test_config(),
        );
        assert!(
            line.contains("\tPASS\t") && line.contains("GPA: 0 [genome]"),
            "{line}"
        );

        // an A-run of 8 starting 3 nt downstream is intrapriming
        genome[2003..2011].copy_from_slice(b"AAAAAAAA");
        let line = status(
            vec![],
            query(),
            &HashMap::new(),
            Some(&genome),
            &test_config(),
        );
        assert!(
            line.contains("\tINTRAPRIMING\t") && line.contains("GPA: 8 [genome]"),
            "{line}"
        );

        // '-' strand: the run is upstream in genome coordinates and made of T
        genome[990..998].copy_from_slice(b"TTTTTTTT");
        let query = model(IP_NAME, 1000, 2000, Strand::Reverse);
        let line = status(
            vec![],
            query,
            &HashMap::new(),
            Some(&genome),
            &test_config(),
        );
        assert!(
            line.contains("\tINTRAPRIMING\t") && line.contains("GPA: 8 [genome]"),
            "{line}"
        );
    }

    #[test]
    fn cn_tag_weighs_a_read_like_that_many_identical_reads() {
        let config = Config {
            recover: true,
            ..test_config()
        };
        let read = |name: &str, end: u64| {
            let mut read = GenePred::from_coords(b"chr1".to_vec(), 100, end, Default::default());
            read.set_name(Some(name.as_bytes().to_vec()));
            read.set_strand(Some(Strand::Forward));
            read
        };
        let run = |queries: Vec<GenePred>| -> Vec<String> {
            let counter = ParallelCounter::default();
            process_component(
                Vec::new(),
                queries,
                &HashMap::new(),
                None,
                &config,
                &counter,
            )
            .into_iter()
            .map(|line| String::from_utf8(line).unwrap())
            .collect()
        };

        // the whole polyA of `pass` is clipped; `ip` has 20 genomic A and no support
        let pass = || read("pass__TC0#PA20#PR20", 1000);
        let weighted = run(vec![pass(), read("ip__TC0#PA0#PR20#CN5", 2000)]);
        let mut identical = vec![pass()];
        identical.extend((0..5).map(|i| read(&format!("ip{i}__TC0#PA0#PR20"), 2000)));
        let identical = run(identical);

        // 5 of 6 reads are intrapriming; counting records, 1 of 2 would not pass `> 0.5`
        assert!(weighted[0].contains("REVIEW"), "{}", weighted[0]);
        assert_eq!(weighted[0], identical[0]);
        let tail = |line: &str| line.split_once('\t').unwrap().1.to_string();
        assert_eq!(tail(&weighted[1]), tail(&identical[1]));
    }
}
