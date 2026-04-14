// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Core module for intron classification
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module provides a comprehensive sub-pipeline for classifying
//! introns within genomic sequences. It different data sources to
//! categorize introns based on their predicted splicing potential
//! and structural characteristics.
//!
//! In essence, this module identifies and characterizes introns from
//! input long-read sequencing data. It performs data integration,
//! collecting splice site prediction scores (from tools like SpliceAI
//! and MaxEntScan), analyzing genomic sequence context, and detecting
//! specific sequence patterns such as RT repeats and NAG motifs. Through
//! a parallel processing approach, each intron is evaluated to determine
//! its "support type", indicating whether it is likely to be a genuine
//! spliced intron, an RT-driven event, or an unclear case requiring
//! further investigation. The final output is a detailed, classified list
//! of introns, enabling deeper insights into alternative splicing and
//! RNA processing.

use anyhow::Result;
use bigtools::{utils::reopen::Reopen, BigWigRead};
use dashmap::DashMap;
use flate2::read::MultiGzDecoder;
use genepred::Bed3;
use hashbrown::HashMap;
use log::info;
use rayon::prelude::*;
use rust_lapper::{Interval, Lapper};
use twobit::TwoBitFile;

use std::borrow::Borrow;
use std::fmt::Debug;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::Mutex;
use std::{
    fs::File,
    io::{Read, Seek},
};

// spliceai-related names/vars
pub const ACCEPTOR_MINUS: &str = "spliceAiAcceptorMinus.bw";
pub const ACCEPTOR_PLUS: &str = "spliceAiAcceptorPlus.bw";
pub const DONOR_MINUS: &str = "spliceAiDonorMinus.bw";
pub const DONOR_PLUS: &str = "spliceAiDonorPlus.bw";
pub const SPLICE_AI_SCORE_RECOVERY_THRESHOLD: f32 = 0.001;

// maxentscan-related names
pub const MINIMUM_ACCEPTOR_LENGTH: usize = 23;
pub const MAXENTSCAN_ACCEPTOR_DB: &str = "db.tsv";
pub const MAXENTSCAN_DONOR_DB: &str = "donor.tsv";
pub const MAXENTSCAN_ACCEPTOR_DB_CONTENTS: &str = include_str!("../assets/db.tsv");
pub const MAXENTSCAN_DONOR_DB_CONTENTS: &str = include_str!("../assets/donor.tsv");

// types
pub type StrandSpliceMap = DashMap<String, DashMap<usize, f32>>;
pub type SharedSpliceMap = (Option<DashMap<usize, f32>>, Option<DashMap<usize, f32>>);
pub type SpliceScores = (Vec<StrandSpliceMap>, Vec<StrandSpliceMap>);
pub type SpliceScoreMap = HashMap<Sequence, Vec<f64>>;
pub type Iv = Interval<u64, ()>;
pub type RepeatIndex = HashMap<Vec<u8>, Lapper<u64, ()>>;

// filanames
pub const INTRON_CLASSIFICATION: &str = "reference_introns.tsv";

// collections
pub const COMPLEMENT: [u8; 128] = {
    let mut nt = [0; 128];
    nt[b'A' as usize] = b'T';
    nt[b'T' as usize] = b'A';
    nt[b'C' as usize] = b'G';
    nt[b'G' as usize] = b'C';
    nt[b'a' as usize] = b't';
    nt[b't' as usize] = b'a';
    nt[b'c' as usize] = b'g';
    nt[b'g' as usize] = b'c';
    nt[b'N' as usize] = b'N';
    nt[b'n' as usize] = b'n';
    nt
};

pub const BGD: [f64; 128] = {
    let mut bgd = [0.0; 128];
    bgd[b'A' as usize] = 0.27;
    bgd[b'T' as usize] = 0.27;
    bgd[b'C' as usize] = 0.23;
    bgd[b'G' as usize] = 0.23;
    bgd
};

pub const CONS1: [f64; 128] = {
    let mut bgd = [0.0; 128];
    bgd[b'A' as usize] = 0.9903;
    bgd[b'C' as usize] = 0.0032;
    bgd[b'G' as usize] = 0.0034;
    bgd[b'T' as usize] = 0.0030;
    bgd
};

pub const CONS2: [f64; 128] = {
    let mut bgd = [0.0; 128];
    bgd[b'A' as usize] = 0.0027;
    bgd[b'C' as usize] = 0.0037;
    bgd[b'G' as usize] = 0.9905;
    bgd[b'T' as usize] = 0.0030;
    bgd
};

/// Splice site type
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
pub enum SpliceSite {
    Donor,
    Acceptor,
}

/// Creates `StrandSpliceMap`s for both plus and minus strands by parsing BigWig files.
///
/// This function takes a directory containing BigWig files for donor and acceptor splice scores for both
/// strands. It then uses `rayon` to parallelize the parsing of these files into `DashMap`s,
/// which are a thread-safe hash map, and returns the results.
///
/// # Arguments
///
/// * `dir`: The path to the directory containing the BigWig files.
/// * `chrs`: A `Vec<String>` of chromosome names to process.
///
/// # Returns
///
/// * A tuple of two `Vec<StrandSpliceMap>`, where the first vector is for the plus strand
///   (donor and acceptor) and the second is for the minus strand.
///
/// # Example
///
/// ```rust,ignore
/// let (plus_scores, minus_scores) = make_splice_map(dir, chrs);
/// ```
pub fn make_splice_map<T: AsRef<std::path::Path> + std::fmt::Debug>(
    dir: T,
    chrs: Vec<String>,
) -> (Vec<StrandSpliceMap>, Vec<StrandSpliceMap>) {
    let plus = vec![
        dir.as_ref().join(DONOR_PLUS),
        dir.as_ref().join(ACCEPTOR_PLUS),
    ];
    let minus = vec![
        dir.as_ref().join(DONOR_MINUS),
        dir.as_ref().join(ACCEPTOR_MINUS),
    ];

    info!("INFO: Parsing BigWigs...");
    info!("INFO: BigWig files: {plus:?}, {minus:?}");
    let (plus, minus) = rayon::join(
        || bigwig_to_map(plus, &chrs),
        || bigwig_to_map(minus, &chrs),
    );

    (plus, minus)
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
    chrs: &Vec<String>,
) -> Vec<DashMap<String, DashMap<usize, f32>>> {
    let total_count = AtomicU32::new(0);
    let rs = Mutex::new(vec![DashMap::new(), DashMap::new()]);

    log::debug!("Will extract scores from the following chromosomes: {chrs:?}");

    // [donor, acceptor]
    bigwigs
        .into_par_iter()
        .zip(vec![SpliceSite::Donor, SpliceSite::Acceptor])
        .for_each(|(bigwig, site)| {
            let acc = DashMap::new();

            let bwread = BigWigRead::open_file(bigwig).expect("ERROR: Cannot open BigWig file");
            let chroms: Vec<_> = bwread.chroms().to_vec();

            chroms.into_par_iter().for_each(|chr| {
                let mut bwread =
                    BigWigRead::reopen(&bwread).expect("ERROR: Cannot re-open BigWig file");

                if !chrs.contains(&chr.name) {
                    log::debug!("Skipping chromosome {chr:?}");
                    return; // INFO: skip chromosomes not in records
                }

                let name = chr.name.clone();
                let length = chr.length;
                let values = bwread
                    .values(&name, 0, length)
                    .expect("ERROR: Cannot read values from BigWig!");

                let mapper = DashMap::new();
                let local_count = AtomicU32::new(0);

                values.into_iter().enumerate().for_each(|(i, v)| {
                    if v >= SPLICE_AI_SCORE_RECOVERY_THRESHOLD {
                        let pos = i;
                        mapper.entry(pos).or_insert(v);
                        local_count.fetch_add(1, Ordering::Relaxed);
                    }
                });

                acc.insert(name, mapper);
                total_count.fetch_add(local_count.load(Ordering::Relaxed), Ordering::Relaxed);
            });

            let mut guard = rs.lock().expect("ERROR: Cannot lock mutex");
            match site {
                SpliceSite::Donor => guard[0] = acc,
                SpliceSite::Acceptor => guard[1] = acc,
            }
        });

    info!(
        "INFO: Parsed and combined {} significant splicing scores from BigWigs!",
        total_count.load(Ordering::Relaxed)
    );

    rs.into_inner()
        .expect("ERROR: Cannot unwrap collection of SpliceAI scores!")
}

/// Fetches and processes splice scores from BigWig files.
///
/// This is a convenience function that wraps `make_splice_map`, handling the case where
/// no splice score files are provided. If files are provided, it calls `make_splice_map` to
/// load the data; otherwise, it returns empty maps.
///
/// # Arguments
///
/// * `splice_scores`: An `Option<T>` with the path to the directory containing splice score BigWigs.
/// * `chrs`: A `Vec<String>` of chromosome names to process.
///
/// # Returns
///
/// * A `SpliceScores` type, which is a tuple of two vectors of `StrandSpliceMap`.
///
/// # Example
///
/// ```rust,ignore
/// let scores = get_splice_scores(splice_scores, chrs);
/// ```
pub fn get_splice_scores<T: AsRef<std::path::Path> + std::fmt::Debug>(
    splice_scores: Option<T>,
    chrs: Vec<String>,
) -> SpliceScores {
    if let Some(splice_scores) = splice_scores {
        make_splice_map(splice_scores, chrs)
    } else {
        log::warn!("No splice scores provided, skipping splice score processing...");
        (vec![DashMap::new()], vec![DashMap::new()])
    }
}

/// Creates a `SharedSpliceMap` for a specific chromosome.
///
/// This function extracts the splice score data for a given chromosome from the full `StrandSpliceMap`s.
/// It is used to get the relevant scores for a single chromosome being processed in parallel.
///
/// # Arguments
///
/// * `chr`: The chromosome name as a string slice.
/// * `splice_plus`: A slice of `StrandSpliceMap` for the plus strand.
/// * `splice_minus`: A slice of `StrandSpliceMap` for the minus strand.
///
/// # Returns
///
/// * A tuple of `SharedSpliceMap` for both plus and minus strands for the specified chromosome.
///
/// # Example
///
/// ```rust,ignore
/// let (plus_chr_scores, minus_chr_scores) = create_splice_map(chr, &splice_plus, &splice_minus);
/// ```
pub fn create_splice_map(
    chr: &str,
    splice_plus: &[StrandSpliceMap],
    splice_minus: &[StrandSpliceMap],
) -> (SharedSpliceMap, SharedSpliceMap) {
    let get_splice_values = |splices: &[DashMap<String, DashMap<usize, f32>>]| {
        (
            splices
                .first()
                .and_then(|s| s.get(chr).map(|v| v.value().clone())),
            splices
                .get(1)
                .and_then(|s| s.get(chr).map(|v| v.value().clone())),
        )
    };
    (
        get_splice_values(splice_plus),
        get_splice_values(splice_minus),
    )
}

/// A trait for parsing MaxEnt scan scores from a line of text.
///
/// This trait defines the behavior for parsing donor and acceptor scores from a TSV-formatted file.
/// `DonorScores` and `AcceptorScores` implement this trait to handle their specific data formats.
///
/// # Methods
///
/// * `parse(line: &str)`: Parses a single line of text into a struct implementing the trait.
/// * `get_sequence()`: Returns the genomic sequence associated with the score.
/// * `get_scores()`: Returns a `Vec<f64>` of the parsed scores.
pub trait SpliceEntropy: Sized + Send + Sync {
    fn parse(line: &str) -> Result<Self, anyhow::Error>
    where
        Self: Sized;
    fn get_sequence(&self) -> Sequence;
    fn get_scores(self) -> Vec<f64>;
}

/// Represents the scores for an acceptor splice site MaxEnt scan.
///
/// This struct holds the sequence and a vector of MaxEnt scores for an acceptor splice site, as parsed
/// from a MaxEntScan database file.
pub struct AcceptorScores {
    seq: Sequence,
    scores: Vec<f64>,
}

/// Implements `SpliceEntropy` for `AcceptorScores`.
///
/// This implementation provides the logic to parse a line from the acceptor score database file.
impl SpliceEntropy for AcceptorScores {
    fn parse(line: &str) -> Result<Self, anyhow::Error> {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 10 {
            return Err(anyhow::anyhow!("ERROR: {} has less than 10 fields!", line));
        }

        let sequence = Sequence::new(parts[0].as_bytes());
        let scores: Vec<f64> = parts[1..]
            .iter()
            .map(|s| s.parse::<f64>().unwrap_or(0.0))
            .collect();

        Ok(AcceptorScores {
            seq: sequence,
            scores,
        })
    }

    fn get_sequence(&self) -> Sequence {
        self.seq.clone()
    }

    fn get_scores(self) -> Vec<f64> {
        self.scores
    }
}

/// Represents the score for a donor splice site MaxEnt scan.
///
/// This struct holds the sequence and a single MaxEnt score for a donor splice site.
pub struct DonorScores {
    seq: Sequence,
    score: f64,
}

/// Implements `SpliceEntropy` for `DonorScores`.
///
/// This implementation provides the logic to parse a line from the donor score database file.
impl SpliceEntropy for DonorScores {
    fn parse(line: &str) -> Result<Self, anyhow::Error> {
        let mut parts = line.split('\t');

        let (sequence, score) = (
            parts
                .next()
                .expect("ERROR: Cannot parse sequence!")
                .as_bytes(),
            parts
                .next()
                .expect("ERROR: Cannot parse donor score!")
                .parse::<f64>()
                .unwrap_or(0.0),
        );

        Ok(DonorScores {
            seq: Sequence::new(sequence),
            score,
        })
    }

    fn get_sequence(&self) -> Sequence {
        self.seq.clone()
    }

    fn get_scores(self) -> Vec<f64> {
        vec![self.score]
    }
}

/// Parses a TSV file into a `SpliceScoreMap`.
///
/// This function reads the content of a TSV file, processes it in parallel, and aggregates the
/// results into a `SpliceScoreMap`. It uses the `SpliceEntropy` trait to handle different
/// formats for donor and acceptor scores.
///
/// # Arguments
///
/// * `contents`: A `String` containing the full content of the TSV file.
///
/// # Type Parameters
///
/// * `K`: A type that implements the `SpliceEntropy` trait.
///
/// # Returns
///
/// * A `Result<SpliceScoreMap, anyhow::Error>` containing the parsed scores.
///
/// # Example
///
/// ```rust,ignore
/// let donor_scores = parse_tsv::<DonorScores>(contents).expect("Failed to parse donor scores");
/// ```
pub fn parse_tsv<K>(contents: String) -> Result<SpliceScoreMap, anyhow::Error>
where
    K: SpliceEntropy,
{
    let tracks = contents
        .par_lines()
        .filter(|row| !row.starts_with("#"))
        .filter_map(|row| K::parse(row).ok())
        .fold(HashMap::new, |mut acc: SpliceScoreMap, splice_site| {
            let entry = acc.entry(splice_site.get_sequence()).or_default();
            entry.extend(splice_site.get_scores());

            acc
        })
        .reduce(HashMap::new, |mut acc, map| {
            for (k, v) in map {
                let acc_v = acc.entry(k).or_insert(Vec::new());
                acc_v.extend(v);
            }
            acc
        });

    info!(
        "INFO: Records parsed: {}",
        tracks.values().flatten().count()
    );
    Ok(tracks)
}

/// Loads MaxEnt scan scores from pre-computed database assets embedded in the binary.
///
/// This function parses the embedded MaxEntScan donor and acceptor tables into
/// `SpliceScoreMap`s using `parse_tsv`.
///
/// # Returns
///
/// * An `Option<(SpliceScoreMap, SpliceScoreMap)>` containing the donor and acceptor score maps.
///
/// # Example
///
/// ```rust,ignore
/// let maxent_scores = load_scan_scores().expect("Failed to load MaxEnt scores");
/// ```
pub fn load_scan_scores() -> Option<(SpliceScoreMap, SpliceScoreMap)> {
    let acceptor_scores = parse_tsv::<AcceptorScores>(MAXENTSCAN_ACCEPTOR_DB_CONTENTS.to_owned())
        .unwrap_or_else(|e| {
            panic!(
                "ERROR: Could not parse embedded acceptor scores from {MAXENTSCAN_ACCEPTOR_DB}! -> {e}"
            )
        });

    let donor_scores = parse_tsv::<DonorScores>(MAXENTSCAN_DONOR_DB_CONTENTS.to_owned())
        .unwrap_or_else(|e| {
            panic!(
                "ERROR: Could not parse embedded donor scores from {MAXENTSCAN_DONOR_DB}! -> {e}"
            )
        });

    Some((donor_scores, acceptor_scores))
}

/// Calculates the combined MaxEnt and consensus score for an acceptor splice site sequence.
///
/// This function computes a score for a 23-base acceptor sequence by combining a MaxEnt score and a
/// consensus sequence score. It returns the natural logarithm of their product.
///
/// # Arguments
///
/// * `seq`: The `Sequence` of the acceptor splice site.
/// * `tables`: A `SpliceScoreMap` containing the MaxEnt score tables.
///
/// # Returns
///
/// * A `f64` representing the combined score.
///
/// # Example
///
/// ```rust,ignore
/// let score = calculate_acceptor_score(&seq, &tables);
/// ```
pub fn calculate_acceptor_score(seq: &Sequence, tables: &SpliceScoreMap) -> f64 {
    if seq.len() != MINIMUM_ACCEPTOR_LENGTH {
        let msg = format!(
            "ERROR: Sequence must be a 23-mer for acceptor score calculation, yours is {}!",
            seq.len()
        );
        log::error!("{}", msg);
        std::process::exit(1);
    }

    let me_score = score_max_ent(seq, tables);
    let c_score = score_consensus_seq(seq);

    if me_score == 0.0 {
        return 0.0;
    }

    (c_score * me_score).log2()
}

/// Calculates the MaxEnt score for an acceptor splice site.
///
/// This function uses the MaxEnt score tables to calculate a score based on specific subsequences
/// of the acceptor splice site sequence.
///
/// # Arguments
///
/// * `seq`: The `Sequence` of the acceptor splice site.
/// * `tables`: A `SpliceScoreMap` containing the MaxEnt score tables.
///
/// # Returns
///
/// * A `f64` representing the MaxEnt score.
///
/// # Example
///
/// ```rust,ignore
/// let me_score = score_max_ent(&seq, &tables);
/// ```
pub fn score_max_ent(seq: &Sequence, tables: &SpliceScoreMap) -> f64 {
    let seq = seq.skip(18, 20);

    let binding = vec![0.0];
    let scores = [
        tables.get(&seq.slice(0, 7)).unwrap_or(&binding).first(),
        tables.get(&seq.slice(7, 14)).unwrap_or(&binding).get(1),
        tables.get(&seq.slice(14, 21)).unwrap_or(&binding).get(2),
        tables.get(&seq.slice(4, 11)).unwrap_or(&binding).get(3),
        tables.get(&seq.slice(11, 18)).unwrap_or(&binding).get(4),
        tables
            .get(&seq.slice_as_seq(4, 7).fill(4))
            .unwrap_or(&binding)
            .get(5),
        tables
            .get(&seq.slice_as_seq(7, 11).fill(3))
            .unwrap_or(&binding)
            .get(6),
        tables
            .get(&seq.slice_as_seq(11, 14).fill(4))
            .unwrap_or(&binding)
            .get(7),
        tables
            .get(&seq.slice_as_seq(14, 18).fill(3))
            .unwrap_or(&binding)
            .get(8),
    ];

    let num: f64 = scores[..5].iter().map(|s| s.unwrap_or(&0.0)).product();
    let den: f64 = scores[5..].iter().map(|s| s.unwrap_or(&0.0)).product();

    num / den
}

/// Calculates the consensus sequence score for an acceptor splice site.
///
/// This function computes a score based on the consensus sequence of the splice site, using pre-defined
/// background probabilities.
///
/// # Arguments
///
/// * `seq`: The `Sequence` of the acceptor splice site.
///
/// # Returns
///
/// * A `f64` representing the consensus score.
///
/// # Example
///
/// ```rust,ignore
/// let c_score = score_consensus_seq(&seq);
/// ```
pub fn score_consensus_seq(seq: &Sequence) -> f64 {
    let nt1 = seq.at_as_bytes(18);
    let nt2 = seq.at_as_bytes(19);

    CONS1[nt1] * CONS2[nt2] / (BGD[nt1] * BGD[nt2])
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

/// Reads multiple files in parallel and concatenates their contents into a single `String`.
///
/// This function leverages Rayon for parallel processing, making it efficient for
/// reading many files concurrently. It panics if any individual file fails to read.
///
/// # Arguments
///
/// * `files` - A vector of paths to the files to read.
///
/// # Returns
///
/// A `Result<String, anyhow::Error>` containing the concatenated contents of all files
/// on success, or an error if parallel reading fails.
///
#[allow(dead_code)]
pub fn par_reader<P: AsRef<Path> + Debug + Sync + Send>(
    files: Vec<P>,
) -> Result<String, anyhow::Error> {
    let contents: Vec<String> = files
        .par_iter()
        .map(|path| {
            reader(path)
                .unwrap_or_else(|e| panic!("ERROR: Could not read file: {:?} -> {:?}", e, path))
        })
        .collect();

    Ok(contents.concat())
}

/// Sequence struct
///
/// This struct is used to store a sequence of nucleotides.
///
/// # Example
/// ```rust, no_run
/// use iso::Sequence;
///
/// let seq = Sequence::new(b"ATCG");
/// assert_eq!(seq.len(), 4);
/// assert_eq!(seq.is_empty(), false);
/// assert_eq!(seq.as_bytes(), b"ATCG");
/// assert_eq!(seq.as_str(), "ATCG");
/// assert_eq!(seq.to_string(), "ATCG");
/// assert_eq!(seq.to_uppercase(), "ATCG");
/// assert_eq!(seq.to_lowercase(), "atcg");
/// assert_eq!(seq.reverse_complement().to_string(), "CGAT");
/// assert_eq!(seq.slice(0, 2), "AT");
/// assert_eq!(seq.slice_as_seq(0, 2).to_string(), "AT");
/// assert_eq!(seq.slice_as_bytes(0, 2), b"AT");
/// assert_eq!(seq.at_as_bytes(0), 65);
/// assert_eq!(seq.fill(2), "AATCG");
/// assert_eq!(seq.skip(1, 3).to_string(), "ACG");
/// ```
#[derive(Debug, PartialEq, Clone, Eq, Hash)]
pub struct Sequence {
    pub seq: Vec<u8>,
}

impl Sequence {
    /// Create a new sequence
    ///
    /// # Example
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.len(), 4);
    /// ```
    pub fn new(seq: &[u8]) -> Self {
        Self { seq: seq.to_vec() }
    }

    /// Decode a sequence from bytes
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.decode(b"ACGT"), "ACGT");
    /// ```
    #[allow(dead_code)]
    pub fn decode(seq: &[u8]) -> Self {
        let base_count = seq.len() * 2;
        let mut capacity = Vec::with_capacity(base_count);

        for i in 0..base_count {
            let byte = seq[i / 2];
            let base_code = if i % 2 == 0 { byte >> 4 } else { byte & 0x0F };

            let base = Sequence::__decode_base(base_code);
            if base != b'=' {
                capacity.push(base);
            }
        }

        Self { seq: capacity }
    }

    /// Get the length of the sequence
    ///
    /// # Example
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.len(), 4);
    /// ```
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    /// Check if the sequence is empty
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.is_empty(), false);
    /// ```
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    /// Get the sequence as a string
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.to_string(), String::from("ATCG"));
    /// ```
    #[allow(dead_code)]
    #[allow(clippy::inherent_to_string_shadow_display)]
    pub fn to_string(&self) -> String {
        String::from_utf8_lossy(&self.seq).to_string()
    }

    /// Get the sequence as uppercase
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"atcg");
    /// assert_eq!(seq.to_uppercase(), "ATCG");
    /// ```
    #[allow(dead_code)]
    pub fn to_uppercase(&self) -> Vec<u8> {
        self.seq.to_ascii_uppercase()
    }

    /// Get the sequence as lowercase
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.to_lowercase(), "atcg");
    /// ```
    #[allow(dead_code)]
    pub fn to_lowercase(&self) -> Vec<u8> {
        self.seq.to_ascii_lowercase()
    }

    /// Get the complement of the sequence
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.complement().to_string(), "GCTA");
    /// ```
    #[allow(dead_code)]
    pub fn complement(&self) -> Self {
        let mut comp = self.seq.to_vec();
        comp.iter_mut().for_each(|c| *c = COMPLEMENT[*c as usize]);

        Self { seq: comp }
    }

    /// Get the reverse complement of the sequence
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.reverse_complement().to_string(), "CGAT");
    /// ```
    pub fn reverse_complement(&self) -> Self {
        let mut rev = self.seq.to_vec();
        rev.reverse();
        rev.make_ascii_uppercase();

        rev.iter_mut().for_each(|c| *c = COMPLEMENT[*c as usize]);

        Self { seq: rev }
    }

    /// Get a slice of the sequence
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.slice(0, 2), "AT");
    /// ```
    pub fn slice(&self, start: usize, end: usize) -> Vec<u8> {
        self.seq[start..end].to_vec()
    }

    /// Get a slice of the sequence as a Sequence struct
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.slice_as_seq(0, 2).to_string(), "AT");
    /// ```
    pub fn slice_as_seq(&self, start: usize, end: usize) -> Self {
        Self {
            seq: self.seq[start..end].to_vec(),
        }
    }

    /// Get a slice of the sequence as bytes
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.slice_as_bytes(0, 2), b"AT");
    /// ```
    #[allow(dead_code)]
    pub fn slice_as_bytes(&self, start: usize, end: usize) -> &[u8] {
        &self.seq[start..end]
    }

    /// Get the ASCII value of a nucleotide at a given index
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.at_as_bytes(0), 65);
    /// ```
    pub fn at_as_bytes(&self, idx: usize) -> usize {
        self.seq[idx] as usize
    }

    /// Skip a given range of the sequence
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.skip(1, 3).to_string(), "ACG");
    /// ```
    pub fn skip(&self, from: usize, to: usize) -> Sequence {
        let mut seq = self.seq[..from].to_vec();
        seq.extend(self.seq[to..].iter());

        Sequence { seq }
    }

    /// Decode a base to a u8
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let base = Sequence::__decode_base(1);
    /// assert_eq!(base, b'A');
    /// ```
    pub fn __decode_base(nt: u8) -> u8 {
        match nt & 0x0f {
            0 => b'=',
            1 => b'A',
            2 => b'C',
            3 => b'M',
            4 => b'G',
            5 => b'R',
            6 => b'S',
            7 => b'V',
            8 => b'T',
            9 => b'W',
            10 => b'Y',
            11 => b'H',
            12 => b'K',
            13 => b'D',
            14 => b'B',
            15 => b'N',
            _ => panic!(
                "{}",
                format!("ERROR: invalid character in sequence: {}", nt)
            ),
        }
    }

    /// Encode a base to a u8
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let base = Sequence::__encode_base(b'A');
    /// assert_eq!(base, 1);
    /// ```
    pub fn __encode_base(nt: u8) -> u8 {
        match nt {
            b'=' => 0,
            b'A' => 1,
            b'C' => 2,
            b'M' => 3,
            b'G' => 4,
            b'R' => 5,
            b'S' => 6,
            b'V' => 7,
            b'T' => 8,
            b'W' => 9,
            b'Y' => 10,
            b'H' => 11,
            b'K' => 12,
            b'D' => 13,
            b'B' => 14,
            _ => panic!(
                "{}",
                format!("ERROR: invalid character in sequence: {}", nt)
            ),
        }
    }

    /// Encode a cannonical base to a u8
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let base = Sequence::__encode_base_2(b'A');
    /// assert_eq!(base, 0);
    /// ```
    pub fn __encode_base_2(nt: u8) -> Option<u8> {
        match nt {
            b'A' => Some(0),
            b'C' => Some(1),
            b'T' => Some(2),
            b'G' => Some(3),
            _ => None,
        }
    }

    /// Encode the reverse of a sequence to a usize
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.reverse_encode(0, 4), vec![3, 2, 1, 0]);
    /// ```
    #[allow(dead_code)]
    pub fn reverse_encode(&self, start: usize, end: usize) -> Vec<usize> {
        self.slice_as_bytes(start, end)
            .iter()
            .rev()
            .filter_map(|b| Self::__encode_base_2(*b))
            .map(|nt| nt as usize)
            .collect::<Vec<usize>>()
    }

    /// Encode the reverse of a sequence to a u8
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"ATCG");
    /// assert_eq!(seq.reverse_encode_u8(0, 4), vec![3, 2, 1, 0]);
    /// ```
    #[allow(dead_code)]
    pub fn reverse_encode_u8(&self, start: usize, end: usize) -> Vec<u8> {
        self.slice_as_bytes(start, end)
            .iter()
            .rev()
            .filter_map(|b| Self::__encode_base_2(*b))
            .collect::<Vec<u8>>()
    }

    /// Fill the sequence with a given kmer
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// use iso::Sequence;
    ///
    /// let seq = Sequence::new(b"TCG");
    /// assert_eq!(seq.fill(2), "AATCG");
    /// ```
    pub fn fill(&self, kmer: usize) -> Vec<u8> {
        let mut seq = b"A".repeat(kmer);
        seq.extend(self.seq.iter());

        seq
    }
}

impl std::fmt::Display for Sequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", String::from_utf8_lossy(&self.seq))
    }
}

impl Borrow<Vec<u8>> for Sequence {
    fn borrow(&self) -> &Vec<u8> {
        &self.seq
    }
}

/// Loads genome sequences from a file (2bit or FASTA format).
///
/// # Arguments
///
/// - `sequence`: Path to the genome file (.fa, .fa.gz, or .2bit)
///
/// # Example
///
/// ```rust,ignore
/// use std::path::PathBuf;
///
/// let genome = get_sequences(PathBuf::from("genome.2bit"));
/// let genome = get_sequences(PathBuf::from("genome.fa"));
/// let genome = get_sequences(PathBuf::from("genome.fa.gz"));
/// ```
pub fn get_sequences(sequence: PathBuf) -> HashMap<Vec<u8>, Vec<u8>> {
    info!("INFO: Reading sequences from file {}", sequence.display());
    match sequence.extension() {
        Some(ext) => match ext.to_str() {
            Some("2bit") => from_2bit(sequence),
            Some("fa") | Some("fasta") | Some("fna") | Some("gz") => from_fa(sequence),
            _ => panic!("ERROR: Unsupported file format"),
        },
        None => panic!("ERROR: No file extension"),
    }
}

/// Loads genome sequences from a 2bit compressed format file.
///
/// # Arguments
///
/// - `twobit`: Path to the 2bit file
///
/// # Example
///
/// ```rust,ignore
/// use std::path::PathBuf;
///
/// let sequences = from_2bit(PathBuf::from("genome.2bit"));
/// let chr1 = sequences.get(b"chr1");
/// ```
fn from_2bit(twobit: PathBuf) -> HashMap<Vec<u8>, Vec<u8>> {
    let genome = TwoBitFile::open_and_read(&twobit).expect("ERROR: Cannot open 2bit file");
    let source = format!("file {}", twobit.display());
    collect_2bit_sequences(genome, &source)
}

/// Loads genome sequences from a 2bit compressed format file.
///
/// # Arguments
///
/// - `genome`: TwoBitFile struct
/// - `source`: Path to the 2bit file
///
/// # Example
///
/// ```rust,ignore
/// use std::path::PathBuf;
///
/// let genome = TwoBitFile::open_and_read(PathBuf::from("genome.2bit")).unwrap();
/// let sequences = collect_2bit_sequences(genome, PathBuf::from("genome.2bit"));
/// let chr1 = sequences.get(b"chr1");
/// ```
fn collect_2bit_sequences<R: Read + Seek>(
    mut genome: TwoBitFile<R>,
    source: &str,
) -> HashMap<Vec<u8>, Vec<u8>> {
    let mut sequences = HashMap::new();
    genome.chrom_names().iter().for_each(|chr| {
        let seq = genome
            .read_sequence(chr, ..)
            .unwrap_or_else(|e| panic!("ERROR: {}", e))
            .as_bytes()
            .to_vec();

        sequences.insert(chr.as_bytes().to_vec(), seq);
    });

    info!("INFO: Read {} sequences from {}", sequences.len(), source);

    sequences
}

/// Loads genome sequences from a FASTA format file (optionally gzipped).
///
/// # Arguments
///
/// - `f`: Path to the FASTA file (.fa or .fa.gz)
///
/// # Example
///
/// ```rust,ignore
/// use std::path::PathBuf;
///
/// let sequences = from_fa(PathBuf::from("genome.fa"));
/// let sequences = from_fa(PathBuf::from("genome.fa.gz"));
/// let chr1 = sequences.get(b"chr1");
/// ```
pub fn from_fa<P: AsRef<Path>>(f: P) -> HashMap<Vec<u8>, Vec<u8>> {
    let path = f.as_ref();
    let file = File::open(path)
        .unwrap_or_else(|e| panic!("ERROR: cannot open FASTA {}: {}", path.display(), e));

    let reader: Box<dyn BufRead> = match path.extension().and_then(|ext| ext.to_str()) {
        Some("gz") => Box::new(BufReader::new(MultiGzDecoder::new(file))),
        _ => Box::new(BufReader::new(file)),
    };

    let source = format!("file {}", path.display());
    parse_fasta_reader(reader, &source)
}

/// Parses a FASTA file into a HashMap of sequences.
///
/// # Arguments
///
/// - `reader`: A BufRead object representing the FASTA file
/// - `source`: Path to the FASTA file
///
/// # Example
///
/// ```rust,ignore
/// use std::path::PathBuf;
///
/// let reader = BufReader::new(File::open(PathBuf::from("genome.fa")).unwrap());
/// let sequences = parse_fasta_reader(reader, PathBuf::from("genome.fa"));
/// let chr1 = sequences.get(b"chr1");
/// ```
fn parse_fasta_reader<R: BufRead>(mut reader: R, source: &str) -> HashMap<Vec<u8>, Vec<u8>> {
    let mut acc = HashMap::new();
    let mut line = Vec::new();
    let mut header: Option<Vec<u8>> = None;
    let mut seq = Vec::new();

    loop {
        line.clear();
        let bytes_read = reader
            .read_until(b'\n', &mut line)
            .unwrap_or_else(|e| panic!("ERROR: cannot read FASTA {}: {}", source, e));

        if bytes_read == 0 {
            break;
        }

        if line.ends_with(b"\n") {
            line.pop();
        }

        if line.ends_with(b"\r") {
            line.pop();
        }

        if line.is_empty() {
            continue;
        }

        if line[0] == b'>' {
            if let Some(prev_header) = header.replace(line[1..].to_vec()) {
                acc.insert(prev_header, std::mem::take(&mut seq));
            }
        } else {
            seq.extend_from_slice(&line);
        }
    }

    if let Some(last_header) = header {
        acc.insert(last_header, seq);
    }

    info!("INFO: Read {} sequences from {}", acc.len(), source);

    acc
}

/// Loads genomic repeats from a BED3 file.
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
pub fn load_repeats(path: PathBuf) -> Option<RepeatIndex> {
    // First pass: collect raw intervals per chrom
    let mut raw: HashMap<Vec<u8>, Vec<Iv>> = HashMap::new();
    let mut counter = 0_usize;

    genepred::Reader::<Bed3>::from_mmap(&path)
        .unwrap_or_else(|e| panic!("ERROR: Could not read repeats file -> {e}!"))
        .records()
        .for_each(|record| {
            let record = record.unwrap_or_else(|e| panic!("ERROR: Could not read record -> {e}!"));
            raw.entry(record.chrom().to_vec()).or_default().push(Iv {
                start: record.start(),
                stop: record.end(),
                val: (),
            });

            counter += 1;
        });
    info!("INFO: Read {} repeats from {}", counter, path.display());

    // Second pass: build a sorted Lapper per chrom (Lapper::new sorts internally)
    let index = raw
        .into_iter()
        .map(|(chrom, ivs)| (chrom, Lapper::new(ivs)))
        .collect();

    info!("INFO: Built repeat index from {}", path.display());
    Some(index)
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
pub fn is_intron_within_repeat(index: &Lapper<u64, ()>, intron: (u64, u64)) -> bool {
    let (start, end) = intron;
    index
        .find(start, end)
        .any(|iv| iv.start <= start && iv.stop >= end)
}

/// Returns true if ≥ `threshold` fraction (e.g. 0.5) of the intron
/// is covered by repeats (handles overlapping/adjacent repeat intervals).
///
/// # Arguments
///
/// * `index`: A `RepeatIndex` containing sorted intervals per chromosome.
/// * `chrom`: A slice of the chromosome name.
/// * `intron`: A tuple of the intron coordinates (start, end).
/// * `threshold`: A float between 0 and 1 representing the minimum fraction of the intron covered by repeats.
///
/// # Returns
///
/// * A boolean indicating whether the intron is covered by repeats.
///
/// # Example
///
/// ```rust,ignore
/// let index = load_repeats(PathBuf::from("repeats.bed3")).unwrap();
/// let chrom = b"chr1";
/// let intron = (100, 200);
/// let covered = intron_covered_by_repeats(&index, &chrom, intron, 0.5);
/// ```
#[inline]
pub fn is_intron_covered_by_repeat(
    index: &Lapper<u64, ()>,
    intron: (u64, u64),
    threshold: f64,
) -> bool {
    let (start, end) = intron;
    let intron_len = end - start;
    if intron_len == 0 {
        return false;
    }

    // INFO: merge overlapping hits on-the-fly to avoid double-counting
    let mut covered: u64 = 0;
    let mut current_end: u64 = 0; // tracks merged interval right edge

    // INFO: find() returns only intervals overlapping [start, end)
    // Ww sort by start to merge correctly — Lapper guarantees sorted order
    for iv in index.find(start, end) {
        let iv_start = iv.start.max(start); // clamp to intron bounds
        let iv_stop = iv.stop.min(end);

        if iv_start >= current_end {
            // INFO: new non-overlapping segment
            covered += iv_stop - iv_start;
            current_end = iv_stop;
        } else if iv_stop > current_end {
            // INFO: extends the current merged segment
            covered += iv_stop - current_end;
            current_end = iv_stop;
        }
        // else: fully contained in already-counted region, skip
    }

    covered as f64 / intron_len as f64 >= threshold
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_load_scan_scores_returns_embedded_tables() {
        let (donor_scores, acceptor_scores) =
            load_scan_scores().expect("ERROR: Could not load scan scores!");

        assert!(!donor_scores.is_empty());
        assert!(!acceptor_scores.is_empty());
    }

    #[test]
    fn test_calculate_consensus_score() {
        let seq = Sequence::new(b"AAAAAAAAAAAAAAAAAAGTAAA");
        let c_score = score_consensus_seq(&seq);

        assert_eq!(c_score, 0.00016425120772946855);

        let seq = Sequence::new(b"TTTAAAAAAAAAAAAAAAGTAAT");
        let c_score = score_consensus_seq(&seq);

        assert_eq!(c_score, 0.00016425120772946854);
    }

    #[test]
    fn test_calculate_max_ent_score() {
        let seq = Sequence::new(b"AAAAAAAAAAAAAAAAAAGTAAA");
        let tables = load_scan_scores().expect("ERROR: Could not load scan scores!");
        let me_score = score_max_ent(&seq, &tables.1);

        assert_eq!(me_score, 0.0003461868180847604);

        let seq = Sequence::new(b"TTTAAAAAAAAAAAAAAAGTAAT");
        let tables = load_scan_scores().expect("ERROR: Could not load scan scores!");
        let me_score = score_max_ent(&seq, &tables.1);

        assert_eq!(me_score, 0.001167787963137549);
    }

    #[test]
    fn test_calculate_acceptor_score() {
        let seq = Sequence::new(b"AAAAAAAAAAAAAAAAAAGTAAA");
        let tables = load_scan_scores().expect("ERROR: Could not load scan scores!");
        let score = calculate_acceptor_score(&seq, &tables.1);

        assert_eq!(score, -24.067969988875006);

        let seq = Sequence::new(b"TTTAAAAAAAAAAAAAAAGTAAT");
        let score = calculate_acceptor_score(&seq, &tables.1);

        assert_eq!(score, -22.313814339694286);
    }

    #[test]
    fn test_calculate_donor_score() {
        let seq = Sequence::new(b"AAGGAAAAA");
        let tables = load_scan_scores().expect("ERROR: Could not load scan scores!");

        let score = tables
            .0
            .get(&seq)
            .expect("ERROR: Could not get donor scores!")
            .first()
            .expect("ERROR: Could not get donor scores!");

        assert_eq!(score, &0.192);
    }
}
