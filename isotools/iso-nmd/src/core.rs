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
use hashbrown::{HashMap, HashSet};
use log::info;
use packbed::{packbed, par_reader, record::Bed4, BedPackage, GenePred};
use rayon::prelude::*;

use std::path::PathBuf;
use std::sync::Arc;

use config::{
    bed_to_map, get_progress_bar, par_write_results, BedColumn, CoordType, OverlapType,
    ParallelCollector, NMD_FREE_READS, NMD_READS, SEP,
};

use crate::cli::Args;

pub const NO_NMD: &str = "NN";
pub const STRONG_NMD: &str = "SN";
pub const WEAK_NMD: &str = "WN";
pub const SN_COLOR: &str = "125,40,40"; // dark-red
pub const WN_COLOR: &str = "213,67,67"; // red

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
pub fn classify_nmd(args: Args) -> Result<(), String> {
    info!("INFO: Classifying NMD classes...");

    let tracks = packbed(
        args.refs.clone(),
        None,
        OverlapType::Exon,
        packbed::PackMode::Paired,
    )
    .unwrap_or_else(|e| panic!("ERROR: failed to pack records in {:?} -> {e}", &args.refs));
    let blacklist = unpack_blacklist(args.blacklist).unwrap_or_default();

    let pb = get_progress_bar(tracks.len() as u64, "Processing...");
    let accumulator = ParallelAccumulator::default();

    tracks.into_par_iter().for_each(|bucket| {
        let chr = bucket.0;
        let components = bucket.1;

        let binding = HashSet::new();
        let banned = blacklist.get(&chr).unwrap_or(&binding);

        process_components(
            components,
            banned,
            args.nmd_distance,
            args.weak_nmd_distance,
            args.atg_distance,
            args.big_exon_dist_to_ej,
            &accumulator,
        );

        pb.inc(1);
    });

    pb.finish_and_clear();
    info!("Reads categorized as NMDs: {}", accumulator.num_nmds());

    par_write_results(
        &accumulator,
        vec![NMD_FREE_READS.into(), NMD_READS.into()],
        Some(args.outdir),
    );

    Ok(())
}

/// Processes a collection of `BedPackage` components in parallel to classify them for NMD.
///
/// This function is the parallel workhorse of the NMD classification. It takes a vector of
/// `BedPackage` objects (representing genomic regions and their associated transcripts) and
/// processes each one concurrently.
///
/// It downcasts each `BedPackage` to a specific type, `(Vec<GenePred>, Vec<GenePred>)`, which
/// is expected to contain the transcript and intron predictions for a given genomic locus.
/// It then calls `process_component` for each individual component to perform the actual
/// NMD classification. The results from each parallel task are accumulated in a thread-safe
/// `ParallelAccumulator`.
///
/// # Arguments
///
/// * `components` - A `Vec` of `Box<dyn BedPackage>`, where each box contains a set of
///   genomic features (e.g., transcripts) for a specific locus.
/// * `banned` - A `HashSet` of genomic coordinates representing blacklisted regions or reads
///   to be excluded from the analysis.
/// * `nmd_distance` - A `u64` representing the distance threshold for strong NMD classification.
/// * `weak_nmd_distance` - An `i64` for the distance threshold for weak NMD classification.
/// * `atg_distance` - A `u64` for the ATG distance threshold used in weak NMD classification.
/// * `big_exon_dist_to_ej` - A `u64` for the big exon distance threshold to the exon-exon junction,
///   used in weak NMD classification.
/// * `accumulator` - A reference to a `ParallelAccumulator` to store the classified reads
///   from all parallel tasks.
fn process_components(
    components: Vec<Box<dyn BedPackage>>,
    banned: &HashSet<(u64, u64)>,
    nmd_distance: u64,
    weak_nmd_distance: i64,
    atg_distance: u64,
    big_exon_dist_to_ej: u64,
    accumulator: &ParallelAccumulator,
) {
    components.into_par_iter().for_each(|comp| {
        let (no_nmd, nmd) = {
            let comp = comp
                .as_any_owned()
                .downcast::<(Vec<GenePred>, Vec<GenePred>)>()
                .expect("ERROR: Could not downcast to IntronPred and GenePred!");

            process_component(
                comp,
                banned,
                nmd_distance,
                weak_nmd_distance,
                atg_distance,
                big_exon_dist_to_ej,
            )
        };

        accumulator.add(no_nmd, nmd);
    });
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
fn process_component(
    component: Box<(Vec<GenePred>, Vec<GenePred>)>,
    banned: &HashSet<(u64, u64)>,
    nmd_distance: u64,
    weak_nmd_distance: i64,
    atg_distance: u64,
    big_exon_dist_to_ej: u64,
) -> (Vec<String>, Vec<String>) {
    let reads = component.0;
    let mut nmd = Vec::new();
    let mut no_nmd = Vec::new();

    for mut read in reads {
        // INFO: skip blacklisted reads
        if banned.contains(&(read.start, read.end)) {
            continue;
        }

        // INFO: noncoding transcripts
        if read.cds_start == read.cds_end {
            // INFO: label = "noNMD", ex_ex_junction_utr = 0, bpUTRtoLastEEJ = 0
            continue;
        }

        let cds_start = read.cds_start;
        let cds_end = read.cds_end;
        let mut nmd_count: i64 = -1;
        let mut _ex_ex_junction_utr: i64 = -1;
        let mut dist_stop_to_next_sj = 0; // INFO: for big exon test
        let mut in_utr = false;
        let mut utr_len = 0;
        let mut cds_len = 0;
        let mut bp_utr_to_last_ex_ex_jct = 0;

        let exons = read.get_exons(); // INFO: already in ascending order

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
                    bp_utr_to_last_ex_ex_jct =
                        utr_len as i64 - (exon_end as i64 - exon_start as i64);
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
                no_nmd.push(read.line);
            }
            STRONG_NMD => {
                let name = format!("{}{SEP}{tag}", read.name);

                read.modify_field(BedColumn::Name.into(), &name);
                read.modify_field(BedColumn::ItemRgb.into(), SN_COLOR);

                nmd.push(read.line);
            }
            WEAK_NMD => {
                let name = format!("{}{SEP}{tag}", read.name);

                read.modify_field(BedColumn::Name.into(), &name);
                read.modify_field(BedColumn::ItemRgb.into(), WN_COLOR);

                nmd.push(read.line);
            }
            _ => {
                log::error!("ERROR: unrecognized tag detected {tag} -> this is a bug");
                std::process::exit(1);
            }
        }
    }

    (no_nmd, nmd)
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
pub fn unpack_blacklist<'a>(paths: Vec<PathBuf>) -> Option<HashMap<String, HashSet<(u64, u64)>>> {
    if paths.is_empty() {
        return None;
    }

    let contents = Arc::new(par_reader(paths).unwrap());
    let tracks = bed_to_map::<Bed4>(contents, CoordType::Bounds).unwrap();

    Some(tracks)
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
    pub no_nmd: DashSet<String>,
    pub nmd: DashSet<String>,
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
    pub fn add(&self, keep: Vec<String>, discard: Vec<String>) {
        for item in keep {
            self.no_nmd.insert(item);
        }
        for item in discard {
            self.nmd.insert(item);
        }
    }
}

/// ParallelCollector trait for ParallelAccumulator
impl<'a> ParallelCollector for ParallelAccumulator {
    /// Get the number of fields in the accumulator
    fn len(&self) -> usize {
        ParallelAccumulator::NUM_FIELDS
    }

    /// Get the a collection of items from the accumulator
    fn get_collections(&self) -> Result<Vec<&DashSet<String>>, Box<dyn std::error::Error>> {
        let mut collections = Vec::with_capacity(ParallelAccumulator::NUM_FIELDS);

        collections.push(&self.no_nmd);
        collections.push(&self.nmd);

        Ok(collections)
    }
}
