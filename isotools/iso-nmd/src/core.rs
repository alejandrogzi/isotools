use hashbrown::{HashMap, HashSet};
use log::info;
use packbed::{packbed, par_reader, record::Bed4, BedPackage, GenePred};
use rayon::prelude::*;

use std::path::PathBuf;
use std::sync::Arc;

use config::{bed_to_map, get_progress_bar, par_write_results, CoordType, OverlapType};

use crate::cli::Args;

pub const NO_NMD: &str = "NN";
pub const STRONG_NMD: &str = "SN";
pub const WEAK_NMD: &str = "WN";
pub const SN_COLOR: &str = ""; // dark-red
pub const WN_COLOR: &str = ""; // red

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

    // let accumulator = ParallelAccumulator::default();
    // let counter = ParallelCounter::default();

    tracks.into_par_iter().for_each(|bucket| {
        let chr = bucket.0;
        let components = bucket.1;

        counter.inc_components(components.len() as u32);

        let binding = HashSet::new();
        let banned = blacklist.get(&chr).unwrap_or(&binding);

        process_components(
            components,
            banned,
            args.nmd_distance,
            args.weak_nmd_distance,
            args.atg_distance,
            args.big_exon_dist_to_ej,
        );

        pb.inc(1);
    });

    pb.finish_and_clear();
    // info!(
    //     "Reads with retained introns: {}",
    //     accumulator.num_retentions()
    // );

    Ok(())
}

fn process_components(
    components: Vec<Box<dyn BedPackage>>,
    banned: &HashSet<String>,
    nmd_distance: u64,
    weak_nmd_distance: i64,
    atg_distance: u64,
    big_exon_dist_to_ej: u64,
) {
    components.into_par_iter().for_each(|mut comp| {
        let comp = comp
            .as_any_mut()
            .downcast_mut::<(Vec<GenePred>, Vec<GenePred>)>()
            .expect("ERROR: Could not downcast to IntronPred and GenePred!");

        // let (keep, discard, review, descriptor) = process_component(comp, banned, counter, recover);
        // accumulator.add(keep, discard, review, descriptor);

        let _ = process_component(
            comp,
            banned,
            nmd_distance,
            weak_nmd_distance,
            atg_distance,
            big_exon_dist_to_ej,
        );
    });
}

fn process_component(
    component: &mut (Vec<GenePred>, Vec<GenePred>),
    banned: &HashSet<String>,
    nmd_distance: u64,
    weak_nmd_distance: i64,
    atg_distance: u64,
    big_exon_dist_to_ej: u64,
) {
    let reads = &mut component.0;

    for read in reads {
        if banned.contains(&read.name) {
            // INFO: skip blacklisted reads
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
        let mut ex_ex_junction_utr: i64 = -1;
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
                ex_ex_junction_utr += 1;

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

        // let line = format!(
        //     "{}\t{}\t{}\t{}",
        //     read.name, tag, ex_ex_junction_utr, bp_utr_to_last_ex_ex_jct,
        // );
        // println!("{}", line);

        // INFO: append tags to read name
        match &tag[..2] {
            NO_NMD => {
                // send to accumulator no_nmd
            }
            STRONG_NMD => {
                let name = format!("{}{SEP}{tag}", read.name);

                query.modify_field(BedColumn::Name.into(), &name);
                query.modify_field(BedColumn::ItemRgb.into(), SN_COLOR);
            }
            WEAK_NMD => {
                let name = format!("{}{SEP}{tag}", read.name);

                query.modify_field(BedColumn::Name.into(), &name);
                query.modify_field(BedColumn::ItemRgb.into(), WN_COLOR);
            }
        }

        // At this point: label, ex_ex_junction_utr, bp_utr_to_last_ex_ex_jct
        // match Perl output exactly (without verbose)
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
pub fn unpack_blacklist<'a>(paths: Vec<PathBuf>) -> Option<HashMap<String, HashSet<(u64, u64)>>> {
    if paths.is_empty() {
        return None;
    }

    let contents = Arc::new(par_reader(paths).unwrap());
    let tracks = bed_to_map::<Bed4>(contents, CoordType::Bounds).unwrap();

    Some(tracks)
}
