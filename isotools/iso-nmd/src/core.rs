//! Core module for detecting non-mediated decays in a query set of reads
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main function for detecting non-mediated decays (NMDs)
//! and processing the components of reads and introns in parallel.
//!
//! In short, identifies and categorizes transcripts based on nonsense-mediated decay (NMD)
//! rules for each read. It processes reads in parallel, filtering out blacklisted entries.
//! For each transcript, it calculates key metrics like the length of the coding sequence
//! and the 3' UTR. Using these metrics and predefined thresholds, it assigns a tag: strong NMD,
//! weak NMD, or no NMD. The final output separates NMD-free reads from NMD reads,
//! with the latter group tagged and color-coded for easy visualization.

use dashmap::DashSet;
use genepred::{Bed12, Bed4, GenePred};
use hashbrown::{HashMap, HashSet};
use log::info;
use rayon::prelude::*;

use std::io::BufWriter;
use std::path::PathBuf;
use std::{fs::File, io::Write};

use crate::cli::Args;

pub const NO_NMD: &str = "NN";
pub const STRONG_NMD: &str = "SN";
pub const WEAK_NMD: &str = "WN";
pub const SN_COLOR: &str = "125,40,40"; // dark-red
pub const WN_COLOR: &str = "213,67,67"; // red
pub const SCALE: u64 = 100_000_000_000; // 100Gb

pub const NMD_FREE_READS: &str = "reads.bed";
pub const NMD_READS: &str = "nmd.bed";

// tag separators
pub const SEP: &str = "#"; // WARN: should change to ':' -> breaks ORF caller [translationai hdf5]
pub const BIG_SEP: &str = "__"; // WARN: should change to '::' -> breaks ORF caller [translationai hdf5]

/// Classifies transcript reads as having strong, weak, or no nonsense-mediated decay (NMD) potential.
///
/// This function orchestrates the entire NMD classification process. It first organizes the input
/// genomic features (e.g., transcripts) into "tracks" or buckets, typically by chromosome,
/// to facilitate parallel processing. It also loads a list of blacklisted reads to be excluded
/// from analysis.
///
/// For each genomic track, it processes the individual transcript components in parallel.
/// The classification is based on several criteria, including the distance of the
/// terminal exon-exon junction from the stop codon, the length of the coding sequence, and
/// the total length of the 3' untranslated region (3' UTR).
///
/// Finally, the classified reads are written to output files, separating those with NMD potential
/// from those without.
///
/// # Arguments
///
/// * `args` - An `Args` struct containing all necessary command-line arguments, such as
///   input file paths, output directory, and various distance thresholds for NMD classification.
///
/// # Returns
///
/// * `Result<(), String>` - Returns `Ok(())` on successful execution.
///   Returns an `Err` with a descriptive string if any step, such as file packing or
///   parallel processing, fails.
pub fn classify_nmd(args: Args) -> Result<(), Box<dyn std::error::Error>> {
    info!("INFO: Classifying NMD classes...");

    let accumulator = ParallelAccumulator::default();
    let blacklist = if let Some(path) = args.blacklist {
        unpack_blacklist(path).unwrap_or_default()
    } else {
        HashMap::new()
    };

    genepred::Reader::<Bed12>::from_mmap(args.bed)?
        .par_records()?
        .for_each(|record| {
            let record = record.expect("Error reading record");

            let chrom =
                std::str::from_utf8(&record.chrom).expect("ERROR: could not convert chrom to str");

            let binding = HashSet::new();
            let banned = blacklist.get(chrom).unwrap_or(&binding);

            process_read(
                record,
                banned,
                args.nmd_distance,
                args.weak_nmd_distance,
                args.atg_distance,
                args.big_exon_dist_to_ej,
                &accumulator,
            );
        });

    info!("Reads categorized as NMDs: {}", accumulator.num_nmds());

    let (no_nmd, nmd) = if let Some(prefix) = args.prefix {
        info!("INFO: Prefix for output files: {prefix}");
        (
            format!("{prefix}_{NMD_FREE_READS}").into(),
            format!("{prefix}_{NMD_READS}").into(),
        )
    } else {
        info!("INFO: No prefix provided, using default output names.");
        (NMD_FREE_READS.into(), NMD_READS.into())
    };

    std::fs::create_dir_all(&args.outdir).map_err(|e| {
        format!(
            "ERROR: Could not create output directory {:?} -> {e}",
            &args.outdir
        )
    })?;
    par_write_results(&accumulator, vec![no_nmd, nmd], Some(args.outdir));

    Ok(())
}

/// Classifies a single genomic component (a set of transcripts) for NMD potential.
///
/// This function performs the core logic of NMD classification. It iterates through each
/// `GenePred` transcript within a given component. It first skips blacklisted and non-coding
/// transcripts.
///
/// For each valid transcript, it calculates several key metrics:
///   - `utr_len`: The length of the 3' UTR.
///   - `cds_len`: The length of the coding sequence (CDS).
///   - `bp_utr_to_last_ex_ex_jct`: The distance from the stop codon to the last exon-exon junction.
///
/// Based on these metrics and the provided distance thresholds, it assigns a classification tag:
/// `NO_NMD`, `STRONG_NMD`, or `WEAK_NMD`. Transcripts with NMD potential are tagged and their
/// names are modified to reflect the classification.
///
/// # Arguments
///
/// * `component` - A `Box` containing a tuple `(Vec<GenePred>, Vec<GenePred>)`, where the first
///   vector contains the transcripts to be classified.
/// * `banned` - A `HashSet` of coordinates for blacklisted reads.
/// * `nmd_distance` - The minimum 3' UTR length for a strong NMD classification.
/// * `weak_nmd_distance` - The maximum distance from the stop codon to the last exon-exon junction
///   for a weak NMD classification.
/// * `atg_distance` - The maximum CDS length for a weak NMD classification.
/// * `big_exon_dist_to_ej` - The maximum distance from the stop codon to the next splice junction
///   for a big exon test, used in weak NMD classification.
///
/// # Returns
///
/// * `(Vec<String>, Vec<String>)` - A tuple containing two vectors of strings. The first vector
///   holds the lines of reads classified as having no NMD potential, and the second
///   contains the lines of reads with NMD potential (both strong and weak).
fn process_read(
    mut read: GenePred,
    banned: &HashSet<(u64, u64)>,
    nmd_distance: u64,
    weak_nmd_distance: i64,
    atg_distance: u64,
    big_exon_dist_to_ej: u64,
    accumulator: &ParallelAccumulator,
) {
    let mut nmd = Vec::new();
    let mut no_nmd = Vec::new();

    // INFO: skip blacklisted reads
    if banned.contains(&(read.start, read.end)) {
        return;
    }

    let mut cds_start = read
        .thick_start()
        .unwrap_or_else(|| panic!("ERROR: thick_start is None for read {:?}", read));

    let mut cds_end = read
        .thick_end()
        .unwrap_or_else(|| panic!("ERROR: thick_end is None for read {:?}", read));

    // INFO: noncoding transcripts
    if cds_start == cds_end {
        // INFO: label = "noNMD", ex_ex_junction_utr = 0, bpUTRtoLastEEJ = 0
        return;
    }

    let mut nmd_count: i64 = -1;
    let mut _ex_ex_junction_utr: i64 = -1;
    let mut dist_stop_to_next_sj = 0; // INFO: for big exon test
    let mut in_utr = false;
    let mut utr_len = 0;
    let mut cds_len = 0;
    let mut bp_utr_to_last_ex_ex_jct = 0;

    let mut exons = read.exons(); // INFO: already in ascending order

    exons = match read.strand() {
        Some(genepred::Strand::Forward) => exons,
        Some(genepred::Strand::Reverse) => {
            let tmp = SCALE - cds_start;
            cds_start = SCALE - cds_end;
            cds_end = tmp;

            // INFO: we scale the coordinates and then sort them
            // INFO: this is equivalent to reversing the exons and reversing the coordinates
            let mut exons = exons
                .iter()
                .map(|(start, end)| (SCALE - end, SCALE - start))
                .collect::<Vec<(u64, u64)>>();

            exons.sort_by(|a, b| a.0.cmp(&b.0));
            exons
        }
        Some(genepred::Strand::Unknown) | None => {
            panic!("ERROR: unexpected strand value: {:?}", read.strand())
        }
    };

    println!("exons: {:?}", exons);
    println!("cds_start: {:?}, cds_end: {:?}", cds_start, cds_end);

    for (i, exon) in exons.iter().enumerate() {
        let exon_start = exon.0;
        let exon_end = exon.1;

        // INFO: Count EEJs in 3'UTR
        if exon_end >= cds_end {
            _ex_ex_junction_utr += 1;

            // INFO: first exon containing stop codon
            if dist_stop_to_next_sj == 0 {
                dist_stop_to_next_sj = exon_end - cds_end;
            }

            if !in_utr {
                utr_len += exon_end - cds_end;
                in_utr = true;
            } else {
                utr_len += exon_end - exon_start;
            }

            if utr_len >= nmd_distance {
                nmd_count += 1;
            }

            // If last exon, compute bpUTRtoLastEEJ
            if i == exons.len() - 1 {
                bp_utr_to_last_ex_ex_jct = utr_len as i64 - (exon_end as i64 - exon_start as i64);
            }
        }

        // INFO: CDS length accumulation
        if exon_end < cds_start || exon_start > cds_end {
            continue; // INFO: skip pure UTR exons
        }

        // INFO: first coding exon
        if exon_end >= cds_start && cds_start >= exon_start {
            if exon_end >= cds_end {
                cds_len += cds_end - cds_start;
            } else {
                cds_len += exon_end - cds_start;
            }
        }
        // INFO: internal coding exon
        else if exon_start > cds_start && exon_end < cds_end {
            cds_len += exon_end - exon_start;
        }
        // INFO: last coding exon
        else if exon_start > cds_start && exon_end >= cds_end {
            cds_len += cds_end - exon_start;
        }
    }

    // INFO: final classification -> tag [NN: no_nmd, SN: strong_nmd, WN: weak_nmd]
    let tag = if nmd_count == 0 || nmd_count == -1 {
        NO_NMD.to_string()
    } else {
        let mut lbl = format!("{STRONG_NMD}{}", nmd_count);
        if bp_utr_to_last_ex_ex_jct <= weak_nmd_distance
            || dist_stop_to_next_sj >= big_exon_dist_to_ej
            || cds_len <= atg_distance
        {
            lbl = format!("{WEAK_NMD}{}", nmd_count);
        }
        lbl
    };

    // INFO: append tags to read name
    match &tag[..2] {
        NO_NMD => {
            // INFO: send to accumulator no_nmd
            let query_line = read.to_bed::<Bed12>();
            no_nmd.push(query_line);
        }
        STRONG_NMD => {
            let query_name = unsafe { std::str::from_utf8_unchecked(read.name().unwrap()) };
            let name = format!("{}{SEP}{tag}", query_name).as_bytes().to_vec();
            read.set_name(Some(name));
            // read.set_item_rgb(SN_COLOR);

            nmd.push(read.to_bed::<Bed12>());
        }
        WEAK_NMD => {
            let query_name = unsafe { std::str::from_utf8_unchecked(read.name().unwrap()) };
            let name = format!("{}{SEP}{tag}", query_name).as_bytes().to_vec();
            read.set_name(Some(name));
            // read.set_item_rgb(WN_COLOR);

            nmd.push(read.to_bed::<Bed12>());
        }
        _ => {
            log::error!("ERROR: unrecognized tag detected {tag} -> this is a bug");
            std::process::exit(1);
        }
    }

    accumulator.add(no_nmd, nmd);
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
pub fn unpack_blacklist<'a>(path: PathBuf) -> Option<HashMap<String, HashSet<(u64, u64)>>> {
    let mut blacklist = HashMap::new();

    genepred::Reader::<Bed4>::from_mmap(path)
        .unwrap_or_else(|e| panic!("ERROR: could not read blacklist -> {e}"))
        .records()
        .for_each(|record| {
            let record = record.unwrap_or_else(|e| panic!("ERROR: could not read record -> {e}"));
            let chrom = std::str::from_utf8(&record.chrom)
                .unwrap_or_else(|e| panic!("ERROR: could not convert chrom to str -> {e}"));
            let start = record.start();
            let end = record.end();

            blacklist
                .entry(chrom.to_string())
                .or_insert_with(HashSet::new)
                .insert((start, end));
        });

    Some(blacklist)
}

/// Parallel accumulator for the processing function
///
/// # Fields
///
/// - `retentions`: A set of strings representing the retentions.
/// - `non_retentions`: A set of strings representing the non-retentions.
/// - `miscellaneous`: A set of strings representing miscellaneous items.
/// - `descriptor`: A map of strings to boxed `ModuleMap` trait objects.
///
/// # Example
///
/// ```rust, no_run
/// let accumulator = ParallelAccumulator::default();
///
/// assert_eq!(accumulator.pass.len(), 0);
/// ```
pub struct ParallelAccumulator {
    pub no_nmd: DashSet<Vec<u8>>,
    pub nmd: DashSet<Vec<u8>>,
}

/// ParallelAccumulator constructor
///
/// # Example
///
/// ```rust, no_run
/// let accumulator = ParallelAccumulator::default();
///
/// assert_eq!(accumulator.no_nmd.len(), 0);
/// assert_eq!(accumulator.nmd.len(), 0);
/// ```
impl Default for ParallelAccumulator {
    fn default() -> Self {
        Self {
            no_nmd: DashSet::new(),
            nmd: DashSet::new(),
        }
    }
}

impl ParallelAccumulator {
    /// Number of fields in the accumulator of type DashSet<String>
    pub const NUM_FIELDS: usize = 2;

    /// Get the number of retentions
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let accumulator = ParallelAccumulator::default();
    /// assert_eq!(accumulator.num_nmds(), 0);
    /// ```
    pub fn num_nmds(&self) -> usize {
        self.nmd.len()
    }

    /// Add items to the accumulator
    ///
    /// # Parameters
    ///
    /// - `keep`: A vector of strings to be retained.
    /// - `discard`: A vector of strings to be discarded.
    /// - `descriptor`: A HashMap of strings to boxed `ModuleMap` trait objects.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let mut accumulator = ParallelAccumulator::default();
    /// accumulator.add(vec!["item1".to_string()], vec!["item2".to_string()], HashMap::new());
    ///
    /// assert_eq!(accumulator.num_retentions(), 1);
    /// ```
    pub fn add(&self, keep: Vec<Vec<u8>>, discard: Vec<Vec<u8>>) {
        for item in keep {
            self.no_nmd.insert(item);
        }
        for item in discard {
            self.nmd.insert(item);
        }
    }
}

pub trait ParallelCollector {
    fn len(&self) -> usize;
    fn get_collections(&self) -> Result<Vec<&DashSet<Vec<u8>>>, Box<dyn std::error::Error>>;
}

/// ParallelCollector trait for ParallelAccumulator
impl<'a> ParallelCollector for ParallelAccumulator {
    /// Get the number of fields in the accumulator
    fn len(&self) -> usize {
        ParallelAccumulator::NUM_FIELDS
    }

    /// Get the a collection of items from the accumulator
    fn get_collections(&self) -> Result<Vec<&DashSet<Vec<u8>>>, Box<dyn std::error::Error>> {
        let mut collections = Vec::with_capacity(ParallelAccumulator::NUM_FIELDS);

        collections.push(&self.no_nmd);
        collections.push(&self.nmd);

        Ok(collections)
    }
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
fn write_pairs(collections: Vec<&DashSet<Vec<u8>>>, filenames: Vec<PathBuf>) {
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
