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
use hashbrown::{HashMap, HashSet};
use log::info;
use packbed::{par_reader, reader, record::Bed4};
use rayon::prelude::*;
use twobit::TwoBitFile;

use std::path::PathBuf;
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::Arc;
use std::sync::Mutex;

use config::{
    bed_to_map, get_progress_bar, CoordType, Sequence, SharedSpliceMap, SpliceScores, SpliceSite,
    StrandSpliceMap, ACCEPTOR_MINUS, ACCEPTOR_PLUS, BGD, CLASSIFY_ASSETS, CONS1, CONS2,
    DONOR_MINUS, DONOR_PLUS, MAXENTSCAN_ACCEPTOR_DB, MAXENTSCAN_DONOR_DB,
};

pub const MINIMUM_ACCEPTOR_LENGTH: usize = 23;

pub type SpliceScoreMap = HashMap<Sequence, Vec<f64>>;

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

    info!("Parsing BigWigs...");
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

    // [donor, acceptor]
    bigwigs
        .into_par_iter()
        .zip(vec![SpliceSite::Donor, SpliceSite::Acceptor])
        .for_each(|(bigwig, site)| {
            let acc = DashMap::new();

            let bwread = BigWigRead::open_file(&bigwig).expect("ERROR: Cannot open BigWig file");
            let chroms: Vec<_> = bwread.chroms().to_vec();

            chroms.into_par_iter().for_each(|chr| {
                let mut bwread =
                    BigWigRead::reopen(&bwread).expect("ERROR: Cannot re-open BigWig file");

                if !chrs.contains(&chr.name) {
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
                    if v >= config::SPLICE_AI_SCORE_RECOVERY_THRESHOLD {
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
        "Parsed and combined {} significant splicing scores from BigWigs!",
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
                .get(0)
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

/// Unpacks a list of blacklist files into a single `HashMap`.
///
/// This function reads one or more blacklist files (in a format like BED4) and parses their contents
/// into a `HashMap` where keys are chromosome names and values are sets of blacklisted intron coordinates.
///
/// # Arguments
///
/// * `paths`: A `Vec<PathBuf>` of paths to the blacklist files.
///
/// # Returns
///
/// * An `Option<HashMap<String, HashSet<(u64, u64)>>>` containing the blacklisted intervals, or `None` if
///   no paths are provided.
///
/// # Example
///
/// ```rust,ignore
/// let blacklist = unpack_blacklist(paths).unwrap_or_default();
/// ```
pub fn unpack_blacklist<'a>(paths: Vec<PathBuf>) -> Option<HashMap<String, HashSet<(u64, u64)>>> {
    if paths.is_empty() {
        return None;
    }

    let contents = Arc::new(par_reader(paths).unwrap());
    let tracks = bed_to_map::<Bed4>(contents, CoordType::Bounds).unwrap();

    Some(tracks)
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
    let pb = get_progress_bar(contents.lines().count() as u64, "Parsing BED12 files");
    let tracks = contents
        .par_lines()
        .filter(|row| !row.starts_with("#"))
        .filter_map(|row| K::parse(row).ok())
        .fold(
            || HashMap::new(),
            |mut acc: SpliceScoreMap, splice_site| {
                let entry = acc.entry(splice_site.get_sequence()).or_default();
                entry.extend(splice_site.get_scores());

                pb.inc(1);
                acc
            },
        )
        .reduce(
            || HashMap::new(),
            |mut acc, map| {
                for (k, v) in map {
                    let acc_v = acc.entry(k).or_insert(Vec::new());
                    acc_v.extend(v);
                }
                acc
            },
        );

    pb.finish_and_clear();
    info!("Records parsed: {}", tracks.values().flatten().count());

    Ok(tracks)
}

/// Loads MaxEnt scan scores from pre-computed database files.
///
/// This function locates the MaxEntScan database files for donors and acceptors, reads them,
/// and parses their contents into `SpliceScoreMap`s using `parse_tsv`.
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
    let assets = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(CLASSIFY_ASSETS);

    let acceptor_scores = parse_tsv::<AcceptorScores>(
        reader(assets.join(MAXENTSCAN_ACCEPTOR_DB)).expect(
            format!(
                "ERROR: Cannot read acceptor scores from {:?}!",
                assets.join(MAXENTSCAN_ACCEPTOR_DB)
            )
            .as_str(),
        ),
    )
    .expect("ERROR: Could not parse acceptor scores!");

    let donor_scores = parse_tsv::<DonorScores>(
        reader(assets.join(MAXENTSCAN_DONOR_DB)).expect("ERROR: Cannot read donor scores!"),
    )
    .expect("ERROR: Could not parse donor scores!");

    Some((donor_scores, acceptor_scores))
}

/// Reads sequences from a 2bit file into a thread-safe map.
///
/// This function opens a 2bit file, reads all chromosome sequences, and stores them in a `DashMap`
/// for efficient, thread-safe access.
///
/// # Arguments
///
/// * `twobit`: The `PathBuf` to the 2bit file.
///
/// # Returns
///
/// * An `Option<DashMap<String, Vec<u8>>>` where keys are chromosome names and values are the sequences.
///
/// # Example
///
/// ```rust,ignore
/// let genome_sequences = get_sequences(twobit_path).expect("Failed to read genome sequences");
/// ```
pub fn get_sequences(twobit: PathBuf) -> Option<DashMap<String, Vec<u8>>> {
    let mut genome = TwoBitFile::open_and_read(twobit).expect("ERROR: Cannot open 2bit file");

    let sequences = DashMap::new();
    genome.chrom_names().iter().for_each(|chr| {
        let seq = genome
            .read_sequence(chr, ..)
            .expect("ERROR: Could not read sequence from .2bit!")
            .as_bytes()
            .to_vec();
        sequences.insert(chr.to_string(), seq);
    });

    Some(sequences)
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
    let scores = vec![
        tables.get(&seq.slice(0, 7)).unwrap_or(&binding).get(0),
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

    let me_score = num / den;

    me_score
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

#[cfg(test)]
mod tests {
    use super::*;

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
            .get(0)
            .expect("ERROR: Could not get donor scores!");

        assert_eq!(score, &0.192);
    }
}
