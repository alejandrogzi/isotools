// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Core module for detecting intron retentions in a query set of reads
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main function for detecting intron retentions
//! and processing the components of reads and introns in parallel.
//!
//! In short, each read is checked for the presence of intron retentions
//! or RT introns. If a read has a true intron retention, it is discarded.
//! If a read has an RT intron, it is also discarded. The veracity of an
//! 'intron' is determined by 'iso-classify', using machine-learning models,
//! ab initio gene prediction, and other heuristics. The process is heavily
//! parallelized to offer fast performance on large datasets.

use genepred::Bed4;
use hashbrown::{HashMap, HashSet};
use rust_lapper::{Interval, Lapper};

use std::{
    fmt::Debug,
    fs::File,
    io::Read,
    path::{Path, PathBuf},
    sync::atomic::{AtomicU32, Ordering},
};

pub type Iv = Interval<u64, ()>;
pub type IntronIndex = HashMap<Vec<u8>, Lapper<u64, ()>>;

/// Parallel counter for the processing function
///
/// # Fields
///
/// - `dirties`: Atomic counter for dirty items.
/// - `components`: Atomic counter for components.
/// - `retentions`: Atomic counter for retentions.
///
/// # Example
///
/// ```rust, no_run
/// let counter = ParallelCounter::default();
///
/// assert_eq!(counter.num_of_comps.load(Ordering::Relaxed), 0);
/// ```
pub struct ParallelCounter {
    pub dirties: AtomicU32,
    pub components: AtomicU32,
    pub retentions: AtomicU32,
}

impl ParallelCounter {
    /// Constructor for ParallelCounter
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    ///
    /// assert_eq!(counter.dirties.load(Ordering::Relaxed), 0);
    /// assert_eq!(counter.components.load(Ordering::Relaxed), 0);
    /// assert_eq!(counter.retentions.load(Ordering::Relaxed), 0);
    /// ```
    fn new() -> Self {
        Self {
            dirties: AtomicU32::new(0),
            components: AtomicU32::new(0),
            retentions: AtomicU32::new(0),
        }
    }

    /// Increment the number of components
    ///
    /// # Parameters
    ///
    /// - `count`: The number of components to add.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_components(5);
    ///
    /// assert_eq!(counter.components.load(Ordering::Relaxed), 5);
    /// ```
    pub fn inc_components(&self, count: u32) {
        self.components.fetch_add(count, Ordering::Relaxed);
    }

    /// Increment the number of dirty items
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_dirty();
    ///
    /// assert_eq!(counter.dirties.load(Ordering::Relaxed), 1);
    /// ```
    pub fn inc_dirty(&self) {
        self.dirties.fetch_add(1, Ordering::Relaxed);
    }

    /// Get the number of dirty items
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_dirty();
    ///
    /// assert_eq!(counter.get_dirty(), 1);
    /// ```
    pub fn get_counters(&self) -> (f64, f64) {
        (
            self.dirties.load(Ordering::Relaxed) as f64,
            self.components.load(Ordering::Relaxed) as f64,
        )
    }

    /// Get the number of dirty items
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_dirty();
    ///
    /// assert_eq!(counter.get_dirty(), 1);
    /// ```
    pub fn get_stat(&self) -> (f64, f64) {
        let (dirties, components) = self.get_counters();
        (dirties, (dirties / components) * 100.0)
    }

    /// Increment the number of retentions
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_retentions();
    ///
    /// assert_eq!(counter.get_retentions(), 1);
    /// ```
    pub fn inc_retentions(&self) {
        self.retentions.fetch_add(1, Ordering::Relaxed);
    }

    /// Get the number of retentions
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let counter = ParallelCounter::new();
    /// counter.inc_retentions();
    ///
    /// assert_eq!(counter.get_retentions(), 1);
    /// ```
    pub fn num_retentions(&self) -> u32 {
        self.retentions.load(Ordering::Relaxed)
    }
}

/// Default implementation for ParallelCounter
///
/// # Example
///
/// ```rust, no_run
/// let counter = ParallelCounter::default();
///
/// assert_eq!(counter.dirties.load(Ordering::Relaxed), 0);
/// assert_eq!(counter.components.load(Ordering::Relaxed), 0);
/// assert_eq!(counter.retentions.load(Ordering::Relaxed), 0);
/// ```
impl Default for ParallelCounter {
    fn default() -> Self {
        Self::new()
    }
}

/// Unpack blacklist from a vector of .bed paths
///
/// # Parameters
///
/// - `paths`: A vector of PathBuf representing the paths to the .bed files.
///
/// # Returns
///
/// - An Option containing a HashMap of strings to HashSets of tuples (u64, u64) if the paths are not empty.
///
/// # Example
///
/// ```rust, no_run
/// let paths = vec![PathBuf::from("path/to/bed1.bed"), PathBuf::from("path/to/bed2.bed")];
/// let result = unpack_blacklist(paths);
///
/// assert!(result.is_some());
/// ```
pub fn unpack_blacklist(path: Option<PathBuf>) -> Option<HashMap<Vec<u8>, HashSet<(u64, u64)>>> {
    if path.is_none() {
        return None;
    }

    let path = path.unwrap();

    let mut blacklist: HashMap<Vec<u8>, HashSet<(u64, u64)>> = HashMap::new();
    let mut tracks = genepred::Reader::<Bed4>::from_mmap(&path).unwrap_or_else(|e| {
        panic!("ERROR: Could not read blacklist from {:?} -> {e}!", path);
    });

    tracks.records().for_each(|record| {
        let record = record.unwrap_or_else(|e| panic!("ERROR: Could not parse record -> {e}!"));

        let chrom = record.chrom();
        let start = record.start();
        let end = record.end();

        blacklist
            .entry(chrom.to_vec())
            .or_default()
            .insert((start, end));
    });

    Some(blacklist)
}

/// Splice site type (U2, U12, Unknown)
///
/// This enum is used to store the type of splice site.
///
/// # Example
///
/// ```rust, no_run
/// use iso::SpliceSite;
///
/// let donor = SpliceSite::Donor;
/// let acceptor = SpliceSite::Acceptor;
/// ```
#[derive(Debug, Eq, PartialEq, Clone)]
pub enum USpliceType {
    U2,
    U12,
    Unknown,
}

/// Implements std::fmt::Display for USpliceType
impl std::fmt::Display for USpliceType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            USpliceType::U2 => write!(f, "U2"),
            USpliceType::U12 => write!(f, "U12"),
            USpliceType::Unknown => write!(f, "UNKNOWN"),
        }
    }
}

/// Implements std::str::FromStr for USpliceType
impl std::str::FromStr for USpliceType {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "U2" => Ok(USpliceType::U2),
            "U12" => Ok(USpliceType::U12),
            "UNKNOWN" => Ok(USpliceType::Unknown),
            _ => Err(format!("ERROR: Unknown splice type -> {s}")),
        }
    }
}

/// Position of the intron in the reference
#[derive(Debug, Eq, PartialEq, Clone)]
#[allow(clippy::upper_case_acronyms)]
pub enum Position {
    CDS,
    UTR,
    Unknown,
}

/// Implements std::fmt::Display for Position
impl std::fmt::Display for Position {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Position::CDS => write!(f, "CDS"),
            Position::UTR => write!(f, "UTR"),
            Position::Unknown => write!(f, "NEW"),
        }
    }
}

/// Implements std::str::FromStr for Position
impl std::str::FromStr for Position {
    type Err = Box<dyn std::error::Error>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "CDS" => Ok(Position::CDS),
            "UTR" => Ok(Position::UTR),
            "NEW" => Ok(Position::Unknown),
            _ => Err("ERROR: Cannot parse position!".into()),
        }
    }
}

/// Support type
///
/// This enum is used to store the type of support.
///
/// # Example
///
/// ```rust, no_run
/// use iso::SupportType;
///
/// let spliced = SupportType::Splicing;
/// let rt = SupportType::RT;
/// let unclear = SupportType::Unclear;
/// ```
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum SupportType {
    Splicing,
    StrongRT,
    WeakRT,
    Unclear,
    Artifact,
}

impl std::fmt::Display for SupportType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SupportType::Splicing => write!(f, "SPLICED"),
            SupportType::StrongRT => write!(f, "STRONG_RT"),
            SupportType::WeakRT => write!(f, "WEAK_RT"),
            SupportType::Unclear => write!(f, "UNCLEAR"),
            SupportType::Artifact => write!(f, "ARTIFACT"),
        }
    }
}

impl std::str::FromStr for SupportType {
    type Err = Box<dyn std::error::Error>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "SPLICED" => Ok(SupportType::Splicing),
            "STRONG_RT" => Ok(SupportType::StrongRT),
            "WEAK_RT" => Ok(SupportType::WeakRT),
            "UNCLEAR" => Ok(SupportType::Unclear),
            "ARTIFACT" => Ok(SupportType::Artifact),
            _ => Err("ERROR: Cannot parse support type!".into()),
        }
    }
}

/// Repeat span type
///
/// This enum is used to store the type of support.
///
/// # Example
///
/// ```rust, no_run
/// use iso::SpanRepeat;
///
/// let within = SpanRepeat::Within;
/// let spans = SpanRepeat::Spans;
/// let null = SpanRepeat::Null;
/// ```
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum SpanRepeat {
    Within,
    Spans,
    Null,
}

/// Implements std::fmt::Display for SpanRepeat
impl std::fmt::Display for SpanRepeat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SpanRepeat::Within => write!(f, "WITHIN_REPEAT"),
            SpanRepeat::Spans => write!(f, "SPANS_REPEAT"),
            SpanRepeat::Null => write!(f, "NO_REPEAT"),
        }
    }
}

/// Implements std::str::FromStr for SpanRepeat
impl std::str::FromStr for SpanRepeat {
    type Err = Box<dyn std::error::Error>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "WITHIN_REPEAT" => Ok(SpanRepeat::Within),
            "SPANS_REPEAT" => Ok(SpanRepeat::Spans),
            "NO_REPEAT" => Ok(SpanRepeat::Null),
            _ => Err("ERROR: Cannot parse span repeat!".into()),
        }
    }
}

/// Holds a variety of statistical and contextual data for a predicted intron.
///
/// This struct contains metrics from different prediction tools and classification
/// systems, used to evaluate the quality and nature of a predicted intron.
#[derive(Debug, PartialEq, Clone)]
pub struct Intron {
    /// Chromosome name
    pub chrom: Vec<u8>,
    /// Start position of the intron (1-based)
    pub start: u64,
    /// End position of the intron (1-based)
    pub end: u64,
    /// Strand of the intron
    pub strand: genepred::Strand,
    /// The frequency of how many reads contain this intron.
    pub seen: usize,
    /// The frequency of how many reads span this intron.
    pub spanned: usize,
    /// SpliceAI score for the donor site.
    pub splice_ai_donor: f32,
    /// SpliceAI score for the acceptor site.
    pub splice_ai_acceptor: f32,
    /// MaxEntScan score for the donor site.
    pub max_ent_donor: f32,
    /// MaxEntScan score for the acceptor site.
    pub max_ent_acceptor: f32,
    /// The sequence around the donor site.
    pub donor_sequence: Vec<u8>,
    /// The sequence around the acceptor site.
    pub acceptor_sequence: Vec<u8>,
    /// The MaxEntScan 9-mer donor context sequence.
    pub donor_context: Vec<u8>,
    /// The MaxEntScan 23-mer acceptor context sequence.
    pub acceptor_context: Vec<u8>,
    /// The classification of the intron's position according to TOGA.
    pub intron_position: Position,
    /// A boolean indicating if the intron is supported by TOGA.
    pub is_toga_supported: Vec<u8>,
    /// A boolean indicating if the intron maintains the reading frame.
    pub is_in_frame: Vec<u8>,
    /// The RT-switch context sequence for the donor site.
    pub donor_rt_context: Vec<u8>,
    /// The RT-switch context sequence for the acceptor site.
    pub acceptor_rt_context: Vec<u8>,
    /// A boolean indicating if the intron is an RT-switch intron.
    pub is_rt_intron: Vec<u8>,
    /// A boolean indicating if the intron is a TOGA-nag intron.
    pub is_nag_intron: Vec<u8>,
    /// A classification of the intron's splice type.
    pub splice_u_type: USpliceType,
    /// A boolean indicating if the intron is within a repeat.
    pub within_repeat: SpanRepeat,
    /// Intron length
    pub intron_len: u64,
    /// A classification of the intron's support type.
    pub support: SupportType,
}

impl Intron {
    pub fn from(record: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let mut fields = record.split('\t');

        let chrom = fields
            .next()
            .unwrap_or_else(|| panic!("ERROR: Could not get chromosome from {}!", record))
            .into();
        let start = fields
            .next()
            .unwrap_or_else(|| panic!("ERROR: Could not get start from {}!", record))
            .parse::<u64>()
            .unwrap_or_else(|_| panic!("ERROR: Could not parse start from {}!", record));
        let end = fields
            .next()
            .unwrap_or_else(|| panic!("ERROR: Could not get end from {}!", record))
            .parse::<u64>()
            .unwrap_or_else(|_| panic!("ERROR: Could not parse end from {}!", record));
        let strand = fields
            .next()
            .unwrap_or_else(|| panic!("ERROR: Could not get strand from {}!", record));
        let seen = fields
            .next()
            .unwrap_or_else(|| panic!("ERROR: Could not get seen from {}!", record))
            .parse::<usize>()
            .unwrap_or_else(|_| panic!("ERROR: Could not parse seen from {}!", record));
        let spanned = fields
            .next()
            .unwrap_or_else(|| panic!("ERROR: Could not get spanned from {}!", record))
            .parse::<usize>()
            .unwrap_or_else(|_| panic!("ERROR: Could not parse spanned from {}!", record));
        let splice_ai_donor = fields
            .next()
            .unwrap_or_else(|| panic!("ERROR: Could not get splice_ai_donor from {}!", record))
            .parse::<f32>()
            .unwrap_or_else(|_| panic!("ERROR: Could not parse splice_ai_donor from {}!", record));
        let splice_ai_acceptor = fields
            .next()
            .unwrap_or_else(|| panic!("ERROR: Could not get splice_ai_acceptor from {}!", record))
            .parse::<f32>()
            .unwrap_or_else(|_| {
                panic!("ERROR: Could not parse splice_ai_acceptor from {}!", record)
            });
        let max_ent_donor = fields
            .next()
            .unwrap_or_else(|| panic!("ERROR: Could not get max_ent_donor from {}!", record))
            .parse::<f32>()
            .unwrap_or_else(|_| panic!("ERROR: Could not parse max_ent_donor from {}!", record));
        let max_ent_acceptor = fields
            .next()
            .unwrap_or_else(|| panic!("ERROR: Could not get max_ent_acceptor from {}!", record))
            .parse::<f32>()
            .unwrap_or_else(|_| panic!("ERROR: Could not parse max_ent_acceptor from {}!", record));
        let donor_sequence = fields
            .next()
            .unwrap_or_else(|| panic!("ERROR: Could not get donor_sequence from {}!", record))
            .into();
        let acceptor_sequence = fields
            .next()
            .unwrap_or_else(|| panic!("ERROR: Could not get acceptor_sequence from {}!", record))
            .into();
        let donor_context = fields
            .next()
            .unwrap_or_else(|| panic!("ERROR: Could not get donor_context from {}!", record))
            .into();
        let acceptor_context = fields
            .next()
            .unwrap_or_else(|| panic!("ERROR: Could not get acceptor_context from {}!", record))
            .into();
        let intron_position = fields
            .next()
            .unwrap_or_else(|| panic!("ERROR: Could not get intron_position from {}!", record))
            .parse::<Position>()
            .unwrap_or_else(|_| panic!("ERROR: Could not parse intron_position from {}!", record));
        // INFO: match TOGA_SUPPORT or NOT_TOGA_SUPPORT
        let is_toga_supported = fields
            .next()
            .unwrap_or_else(|| panic!("ERROR: Could not get is_toga_supported from {}!", record))
            .into();
        let is_in_frame = fields
            .next()
            .unwrap_or_else(|| panic!("ERROR: Could not get is_in_frame from {}!", record))
            .into();
        let donor_rt_context = fields
            .next()
            .unwrap_or_else(|| panic!("ERROR: Could not get donor_rt_context from {}!", record))
            .into();
        let acceptor_rt_context = fields
            .next()
            .unwrap_or_else(|| panic!("ERROR: Could not get acceptor_rt_context from {}!", record))
            .into();
        let is_rt_intron = fields
            .next()
            .unwrap_or_else(|| panic!("ERROR: Could not get is_rt_intron from {}!", record))
            .into();
        let is_nag_intron = fields
            .next()
            .unwrap_or_else(|| panic!("ERROR: Could not get is_nag_intron from {}!", record))
            .into();
        let splice_u_type = fields
            .next()
            .unwrap_or_else(|| panic!("ERROR: Could not get splice_u_type from {}!", record))
            .parse::<USpliceType>()
            .unwrap_or_else(|_| panic!("ERROR: Could not parse splice_u_type from {}!", record));
        let within_repeat = fields
            .next()
            .unwrap_or_else(|| panic!("ERROR: Could not get within_repeat from {}!", record))
            .parse::<SpanRepeat>()
            .unwrap_or_else(|_| panic!("ERROR: Could not parse within_repeat from {}!", record));
        let intron_len = fields
            .next()
            .unwrap_or_else(|| panic!("ERROR: Could not get intron_len from {}!", record))
            .parse::<u64>()
            .unwrap_or_else(|_| panic!("ERROR: Could not parse intron_len from {}!", record));
        let support = fields
            .next()
            .unwrap_or_else(|| panic!("ERROR: Could not get support from {}!", record))
            .parse::<SupportType>()
            .unwrap_or_else(|_| panic!("ERROR: Could not parse support from {}!", record));

        let strand = match strand {
            "+" => genepred::Strand::Forward,
            "-" => genepred::Strand::Reverse,
            _ => panic!("ERROR: Unknown strand -> {strand} from {}", record),
        };

        Ok(Self {
            chrom,
            start,
            end,
            strand,
            seen,
            spanned,
            splice_ai_donor,
            splice_ai_acceptor,
            max_ent_donor,
            max_ent_acceptor,
            donor_sequence,
            acceptor_sequence,
            donor_context,
            acceptor_context,
            intron_position,
            is_toga_supported,
            is_in_frame,
            donor_rt_context,
            acceptor_rt_context,
            is_rt_intron,
            is_nag_intron,
            splice_u_type,
            within_repeat,
            intron_len,
            support,
        })
    }
}

impl std::fmt::Display for Intron {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            std::str::from_utf8(&self.chrom).unwrap_or("NULL"),
            self.start,
            self.end,
            self.strand,
            self.seen,
            self.spanned,
            self.splice_ai_donor,
            self.splice_ai_acceptor,
            self.max_ent_donor,
            self.max_ent_acceptor,
            std::str::from_utf8(&self.donor_sequence).unwrap_or("NULL"),
            std::str::from_utf8(&self.acceptor_sequence).unwrap_or("NULL"),
            std::str::from_utf8(&self.donor_context).unwrap_or("NULL"),
            std::str::from_utf8(&self.acceptor_context).unwrap_or("NULL"),
            self.intron_position,
            std::str::from_utf8(&self.is_toga_supported).unwrap_or("NULL"),
            std::str::from_utf8(&self.is_in_frame).unwrap_or("NULL"),
            std::str::from_utf8(&self.donor_rt_context).unwrap_or("NULL"),
            std::str::from_utf8(&self.acceptor_rt_context).unwrap_or("NULL"),
            std::str::from_utf8(&self.is_rt_intron).unwrap_or("NULL"),
            std::str::from_utf8(&self.is_nag_intron).unwrap_or("NULL"),
            self.splice_u_type,
            self.within_repeat,
            self.intron_len,
            self.support
        )
    }
}

/// Loads intron DB
///
/// This function reads a BED3 file containing genomic repeats and builds an index
/// of sorted intervals per chromosome. The index is returned as a `HashMap` where
/// the keys are chromosome names and the values are `Lapper` instances.
///
/// # Arguments
///
/// * `path`: Path to the BED3 file containing repeats.
///
/// # Returns
///
/// * An `Option<RepeatIndex>` containing the repeat index.
///
/// # Example
///
/// ```rust,ignore
/// let repeats = load_repeats(PathBuf::from("repeats.bed3"));
/// ```
pub fn load_introns<P: AsRef<Path> + Debug + Copy>(
    path: P,
) -> Option<(IntronIndex, HashMap<Vec<u8>, Intron>)> {
    let mut counter = 0_usize;

    // INFO: build two collections: 1) lookup Intron hash, 2) Lapper index
    let mut iv_collector: HashMap<Vec<u8>, Vec<Iv>> = HashMap::new();
    let mut introns: HashMap<Vec<u8>, Intron> = HashMap::new();

    let lines = reader(path)
        .unwrap_or_else(|e| panic!("ERROR: Could not read introns file {path:?} -> {e}!"));
    lines.lines().for_each(|line| {
        // INFO: convert line to Intron
        let record = Intron::from(line)
            .unwrap_or_else(|e| panic!("ERROR: Could not parse line {line:?} -> {e}!"));

        // INFO: add to lookup, key fmt: chrom:start-end(strand) as bytes
        let mut key = Vec::new();
        key.extend_from_slice(&record.chrom);
        key.extend_from_slice(b":");
        key.extend_from_slice(&record.start.to_string().as_bytes());
        key.extend_from_slice(b"-");
        key.extend_from_slice(&record.end.to_string().as_bytes());
        key.extend_from_slice(b"(");
        key.extend_from_slice(&record.strand.to_string().as_bytes());
        key.extend_from_slice(b")");

        // INFO: add to interval collector for Lapper index
        let iv = Iv {
            start: record.start,
            stop: record.end,
            val: (),
        };

        let mut chr_stranded = Vec::new();
        chr_stranded.extend_from_slice(&record.chrom);
        chr_stranded.push(b':');
        chr_stranded.extend_from_slice(&record.strand.to_string().as_bytes());
        iv_collector.entry(chr_stranded).or_default().push(iv);

        introns.insert(key, record);

        counter += 1;
    });

    log::info!("INFO: Read {} introns from {:?}", counter, path);

    // Second pass: build a sorted Lapper per chrom (Lapper::new sorts internally)
    let index = iv_collector
        .into_iter()
        .map(|(chrom, ivs)| (chrom, Lapper::new(ivs)))
        .collect();

    log::info!("INFO: Built intron index from {:?}", path);
    Some((index, introns))
}

/// Returns true if the intron falls ENTIRELY within a single repeat interval.
///
/// # Arguments
///
/// * `index`: A `RepeatIndex` containing sorted intervals per chromosome.
/// * `chrom`: A slice of the chromosome name.
/// * `intron`: A tuple of the intron coordinates (start, end).
///
/// # Returns
///
/// * A boolean indicating whether the intron is within a repeat interval.
///
/// # Example
///
/// ```rust,ignore
/// let index = load_repeats(PathBuf::from("repeats.bed3")).unwrap();
/// let chrom = b"chr1";
/// let intron = (100, 200);
/// let within_repeat = intron_within_repeat(&index, &chrom, intron);
/// ```
#[inline]
pub fn is_intron_retention(index: &Lapper<u64, ()>, exon: (u64, u64)) -> bool {
    let (start, end) = exon;
    index
        .find(start, end)
        .any(|iv| start <= iv.start && iv.stop <= end)
}

/// Returns an iterator over the intron retentions within an exon
///
/// # Arguments
///
/// * `index` - The index to search for intron retentions.
/// * `exon` - The exon to search for intron retentions.
///
/// # Returns
///
/// * An iterator over the intron retentions within the exon.
///
/// # Example
///
/// ```rust, no_run
/// let index = load_repeats(PathBuf::from("repeats.bed3")).unwrap();
/// let chrom = b"chr1";
/// let exon = (100, 200);
/// let introns = intron_retentions(&index, &exon);
///
/// for intron in introns {
///     println!("Intron retention: {:?}", intron);
/// }
/// ```
#[inline]
pub fn intron_retentions<'a>(
    index: &'a Lapper<u64, ()>,
    exon: (u64, u64),
) -> impl Iterator<Item = (u64, u64)> + 'a {
    let (start, end) = exon;
    index
        .find(start, end)
        .filter(move |iv| start <= iv.start && iv.stop <= end)
        .map(|iv| (iv.start, iv.stop))
}

/// Reads the entire content of a file into a `String`.
///
/// This function provides a basic utility for synchronously reading a file's
/// contents. It's generic over any type that can be converted to a `Path` and
/// is `Debug` printable.
///
/// # Arguments
///
/// * `file` - The path to the file to read.
///
/// # Returns
///
/// A `Result<String, Box<dyn std::error::Error>>` containing the file's
/// contents on success, or an error if the file cannot be opened or read.
///
pub fn reader<P: AsRef<Path> + Debug>(file: P) -> Result<String, Box<dyn std::error::Error>> {
    let mut file = File::open(file)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    Ok(contents)
}
