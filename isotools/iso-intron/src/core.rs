//! Core module for detecting intron retentions in a query set of reads
//!
//!
//! expected output type: HashMap/DashMap<ReadName, RetentionDescriptor>
//! expected output structure: ReadName -> <T> where T: ModuleDescriptor
//!     - ReadName is the name of the read
//!     - <T> is a struct implementing ModuleDescriptor trait
//!
//! expected output example:
//!     {
//!        "read1": {
//!           "intron_retention": true,
//!           "retains_rt_intron": true,
//!           "is_retention_supported": false,
//!           "is_retention_supported_map": [false, false],
//!           "intron_support_ratio": [0.23, 0.4],
//!           "exon_support_ratio": [0.65, 0.43],
//!           "number_of_retentions": 2,
//!           "number_of_unsupported_retentions": 2,
//!           "number_of_support_retentions": 0,
//!           "location_of_retention": ["chr1:1000-2000", "chr1:3000-4000"],
//!           "retention_in_cds": [false, true],
//!           "retention_in_utr": [true, false],
//!         }
//!     }

use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use anyhow::Result;
use hashbrown::{HashMap, HashSet};
use log::{info, warn};
use packbed::record::IntronPosition;
use packbed::{packbed, BedPackage, GenePred, IntronPred, RefGenePred};
use rayon::prelude::*;
use serde_json::Value;

use crate::cli::Args;
use crate::utils::{unpack_blacklist, write_results, ParallelAccumulator, ParallelCounter};

use config::{
    get_progress_bar, IntronRetentionValue, ModuleDescriptor, ModuleMap, ModuleType, OverlapType,
    INTRON_RETENTION_RECOVERY_THRESHOLD, RETENTION_RATIO_THRESHOLD, SCALE,
    SPLICE_AI_SCORE_RECOVERY_THRESHOLD,
};

/// Detects intron retentions in an query set of reads
/// based on a provided reference set.
///
/// # Example
/// ```ignore
/// use iso_intron::cli::Args;
/// use iso_intron::core::detect_intron_retentions;
///
/// let args = Args {
///   refs: vec![],
///   query: vec![],
///   blacklist: vec![],
///   plot: true,
///   threads: 4,
/// };
///
/// detect_intron_retentions(args);
/// ```
///
/// # Description
///
/// Entry point for detecting intron retentions in a query set of reads.
/// The function packs the reference and query transcripts into a collection
/// and processes each bucket of overlapping transcripts to identify
/// retained introns. Here we interpret intron retentions as the complete
/// exon overlap with a reference intron:
///
/// ```text
///     5'                          3'
///     XXXXX-----XXXXXX--------XXXXXX
///                     ^^^^^^^^
///     XXXXX-----XXXXXXXXXXXXXXXXXXXX
/// ```
pub fn detect_intron_retentions(args: Args) -> Result<()> {
    info!("Detecting intron retentions...");

    let tracks = packbed(
        args.refs,
        Some(args.query),
        OverlapType::Exon,
        packbed::PackMode::Query,
    )?;
    let blacklist = unpack_blacklist(args.blacklist).unwrap_or_default();

    let pb = get_progress_bar(tracks.len() as u64, "Processing...");

    let accumulator = ParallelAccumulator::default();
    let counter = ParallelCounter::default();

    tracks.into_par_iter().for_each(|bucket| {
        let chr = bucket.0;
        let components = bucket.1;

        counter.inc_components(components.len() as u32);

        let binding = HashSet::new();
        let banned = blacklist.get(&chr).unwrap_or(&binding);

        process_components(components, banned, &accumulator, &counter);

        pb.inc(1);
    });

    pb.finish_and_clear();
    info!(
        "Reads with retained introns: {}",
        accumulator.num_retentions()
    );

    write_results(&accumulator);

    Ok(())
}

#[inline(always)]
fn process_components(
    components: Vec<Box<dyn BedPackage>>,
    banned: &HashSet<(u64, u64)>,
    accumulator: &ParallelAccumulator,
    counter: &ParallelCounter,
) {
    components.into_par_iter().for_each(|mut comp| {
        let comp = comp
            .as_any_mut()
            .downcast_mut::<(Vec<IntronPred>, Vec<GenePred>)>()
            .expect("ERROR: Could not downcast to IntronPred!");

        let (keep, discard) = process_component(comp, banned, counter);

        discard.into_iter().for_each(|r| {
            accumulator.retentions.insert(r);
        });

        keep.into_iter().for_each(|nr| {
            accumulator.non_retentions.insert(nr);
        });
    });
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum IntronModuleReadCategory {
    RT,      // RT-intron
    Unclear, // Artifact
    Clean,   // None
    IR,      // Retention
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum IntronModuleReadAction {
    // Review,  // Unclear
    Discard, // RT-intron
    Keep,
}

#[inline(always)]
pub fn process_component(
    comp: &mut (Vec<IntronPred>, Vec<GenePred>),
    ban: &HashSet<(u64, u64)>,
    counter: &ParallelCounter,
) -> (
    Vec<String>,
    Vec<String>,
    // HashMap<String, Box<dyn ModuleMap>>,
) {
    let mut keep = Vec::new();
    let mut discard = Vec::new();

    let introns = &comp.0;
    let reads = &comp.1;

    // INFO: convert Vec<IntronPred> into HashMap<(u64, u64), IntronPred>
    let intron_map = introns
        .iter()
        .map(|intron| ((intron.start, intron.end), intron))
        .collect::<HashMap<_, _>>();

    // INFO: for every read -> see which introns are retained and which the read has!
    for read in reads {
        // INFO: Determine if read has RT introns
        let read_introns = read.get_introns();
        let mut intronic_status = IntronModuleReadAction::Keep;
        for intron in read_introns {
            let hit = intron_map
                .get(&intron)
                .expect("ERROR: Intron not found, this is likely a bug!");

            match hit.stats.support {
                config::SupportType::RT => {
                    // INFO: if read is RT, discard
                    intronic_status = IntronModuleReadAction::Discard;
                }
                config::SupportType::Splicing | config::SupportType::Unclear => {} // INFO: Do nothing
            };

            if intronic_status == IntronModuleReadAction::Discard {
                break;
            }
        }

        let read_exons = read.get_exons();
        let mut exonic_status = IntronModuleReadAction::Keep;
        for exon in read_exons {
            let exon_start = exon.0;
            let exon_end = exon.1;

            for (intron, stats) in intron_map.iter() {
                let intron_start = intron.0;
                let intron_end = intron.1;

                // INFO: early exit to avoid unnecessary checks
                if intron_end < exon_start {
                    break;
                }

                if exon_start < intron_start && intron_end < exon_end {
                    match stats.stats.support {
                        config::SupportType::RT => {} // INFO: retaining a false intron is not an IR, do nothing
                        config::SupportType::Unclear | config::SupportType::Splicing => {
                            match stats.stats.intron_position {
                                IntronPosition::CDS => {
                                    // INFO: if intron in frame keep read --> splice variant producing a longer protein
                                    // INFO: else do not annotate read (exon structure is flawed)
                                    if !stats.stats.is_in_frame {
                                        exonic_status = IntronModuleReadAction::Discard;
                                    }
                                }
                                // INFO: keep read --> variant would not affect CDS
                                IntronPosition::UTR
                                | IntronPosition::Unknown
                                | IntronPosition::Mixed => {}
                            }
                        }
                    }
                }
            }
        }

        match (intronic_status, exonic_status) {
            (IntronModuleReadAction::Discard, _) | (_, IntronModuleReadAction::Discard) => {
                counter.inc_retentions();
                discard.push(read.line().to_owned());
            }
            _ => {
                keep.push(read.line().to_owned());
            }
        };
    }

    (keep, discard)
}

// #[inline(always)]
// pub fn process_component(
//     comp: &mut Vec<GenePred>,
//     ban: &HashSet<(u64, u64)>,
//     plot: bool,
//     recover: bool,
//     counter: &ParallelCounter,
// ) -> (
//     Vec<String>,
//     Vec<String>,
//     Option<String>,
//     HashMap<String, Box<dyn ModuleMap>>,
//     bool,
// ) {
//     let mut hits = Vec::new();
//     let mut pass = Vec::new();
//     let mut blocks = if plot { Some(String::new()) } else { None };

//     let mut tmp_dirt = Vec::new();
//     let mut exon_owners = BTreeSet::new();
//     let mut introns_owned = BTreeMap::new();

//     let mut descriptor = HashMap::new();

//     let refs = comp.0;
//     let queries = comp.1;

//     let ref_size = refs.reads.len();
//     let query_size = queries.len();
//     let comp_size = ref_size + query_size;
//     let (mut count, totals) = (0_f32, queries.len() as f32);

//     for query in queries.iter() {
//         descriptor.insert(
//             query.name.clone(),
//             ModuleDescriptor::with_schema(ModuleType::IntronRetention),
//         );

//         let query = Arc::new(query);
//         let mut hit: bool = false;

//         let five_utr = query.get_five_utr();
//         let three_utr = query.get_three_utr();

//         let mut number_of_retentions = 0_u32;
//         let mut number_of_true_retentions = 0_u32;
//         let mut number_of_partial_retentions = 0_u32;
//         let mut number_of_false_retentions = 0_u32;

//         let mut location_of_retentions = Vec::new();
//         let mut retention_in_cds = Vec::new();
//         let mut retention_in_utr = Vec::new();
//         let mut retention_in_frame = Vec::new();
//         let mut donor_scores = Vec::new();
//         let mut acceptor_scores = Vec::new();

//         for exon in &query.exons {
//             for intron in &refs.introns {
//                 if exon.1 < intron.0 {
//                     break;
//                 }

//                 if ban.contains(intron) || ban.contains(&(SCALE - intron.1, SCALE - intron.0)) {
//                     continue;
//                 }

//                 if exon.0 < intron.0 && exon.1 > intron.1 {
//                     if !hit {
//                         count += 1.0;
//                         tmp_dirt.push(query.clone());
//                     }

//                     hit = true;
//                     exon_owners.insert((exon.0, exon.1));
//                     introns_owned
//                         .entry((intron.0, intron.1))
//                         .or_insert((0.0, 0.0));

//                     number_of_retentions += 1;

//                     if (five_utr.0 <= intron.0 && intron.0 <= five_utr.1)
//                         || (three_utr.0 <= intron.1 && intron.1 <= three_utr.1)
//                     {
//                         retention_in_cds.push(false);
//                         retention_in_utr.push(true);
//                     } else {
//                         retention_in_cds.push(true);
//                         retention_in_utr.push(false);
//                     }

//                     if (intron.1 - intron.0) % 3 == 0 {
//                         retention_in_frame.push(true);
//                     } else {
//                         retention_in_frame.push(false);
//                     }

//                     match query.strand {
//                         '+' => {
//                             location_of_retentions
//                                 .push(format!("{}:{}-{}", query.chrom, intron.0, intron.1,));

//                             if plot {
//                                 Bed4::from(
//                                     query.chrom.clone(),
//                                     intron.0 - 1,
//                                     intron.1 + 1,
//                                     query.name.clone(),
//                                 )
//                                 .send(&mut blocks.as_mut().unwrap())
//                             }
//                         }
//                         '-' => {
//                             location_of_retentions.push(format!(
//                                 "{}:{}-{}",
//                                 query.chrom,
//                                 SCALE - intron.1,
//                                 SCALE - intron.0,
//                             ));

//                             if plot {
//                                 Bed4::from(
//                                     query.chrom.clone(),
//                                     SCALE - intron.1 - 1,
//                                     SCALE - intron.0 + 1,
//                                     query.name.clone(),
//                                 )
//                                 .send(&mut blocks.as_mut().unwrap())
//                             }
//                         }
//                         _ => panic!("Invalid strand"),
//                     }

//                     if let Some(scores) = splice_scores {
//                         if let Some(donor_score_map) = scores.0.as_ref() {
//                             let acceptor_score_map = scores.1.as_ref().unwrap();

//                             let (intron_donor, intron_acceptor) = match query.strand {
//                                 // donor [-1 to match bigtools coords]
//                                 '+' => (intron.0 as usize - 1, intron.1 as usize),
//                                 // acceptor [-1 to match bigtools coords]
//                                 '-' => {
//                                     ((SCALE - intron.0) as usize, (SCALE - intron.1) as usize - 1)
//                                 }
//                                 _ => panic!("ERROR: Invalid strand in query component!"),
//                             };

//                             let (donor_score, acceptor_score) = (
//                                 donor_score_map
//                                     .get(&intron_donor)
//                                     .map(|r| *r)
//                                     .unwrap_or(0.0),
//                                 acceptor_score_map
//                                     .get(&intron_acceptor)
//                                     .map(|r| *r)
//                                     .unwrap_or(0.0),
//                             );

//                             donor_scores.push(donor_score);
//                             acceptor_scores.push(acceptor_score);

//                             if donor_score >= 0.02 && acceptor_score >= 0.02 {
//                                 number_of_true_retentions += 1;
//                             } else if donor_score >= 0.01 && donor_score >= 0.01 {
//                                 number_of_partial_retentions += 1;
//                             } else {
//                                 number_of_false_retentions += 1;
//                             }

//                             introns_owned.entry((intron.0, intron.1)).and_modify(|e| {
//                                 *e = (donor_score, acceptor_score);
//                             });
//                         }
//                     }

//                     continue;
//                 } else {
//                     continue;
//                 }
//             }
//         }

//         let handle = descriptor.get_mut(&query.name).unwrap();

//         handle
//             .set_value(
//                 Box::new(IntronRetentionValue::ComponentSize),
//                 Value::Number(comp_size.into()),
//             )
//             .ok();
//         handle
//             .set_value(
//                 Box::new(IntronRetentionValue::RefComponentSize),
//                 Value::Number(ref_size.into()),
//             )
//             .ok();
//         handle
//             .set_value(
//                 Box::new(IntronRetentionValue::QueryComponentSize),
//                 Value::Number(query_size.into()),
//             )
//             .ok();

//         if hit {
//             let line = query.line().to_owned();
//             hits.push(line);

//             handle
//                 .set_value(
//                     Box::new(IntronRetentionValue::IsIntronRetention),
//                     Value::Bool(true),
//                 )
//                 .ok();
//             handle
//                 .set_value(
//                     Box::new(IntronRetentionValue::NumberOfRetentions),
//                     Value::Number(number_of_retentions.into()),
//                 )
//                 .ok();
//             handle
//                 .set_value(
//                     Box::new(IntronRetentionValue::RetentionLocation),
//                     Value::Array(
//                         location_of_retentions
//                             .into_iter()
//                             .map(Value::String)
//                             .collect(),
//                     ),
//                 )
//                 .ok();
//             handle
//                 .set_value(
//                     Box::new(IntronRetentionValue::IsRetentionInCds),
//                     Value::Array(retention_in_cds.into_iter().map(Value::Bool).collect()),
//                 )
//                 .ok();
//             handle
//                 .set_value(
//                     Box::new(IntronRetentionValue::IsRetentionInUtr),
//                     Value::Array(retention_in_utr.into_iter().map(Value::Bool).collect()),
//                 )
//                 .ok();
//             handle
//                 .set_value(
//                     Box::new(IntronRetentionValue::IsIntronRetainedInFrame),
//                     Value::Array(retention_in_frame.into_iter().map(Value::Bool).collect()),
//                 )
//                 .ok();
//             handle
//                 .set_value(
//                     Box::new(IntronRetentionValue::RetentionDonorScore),
//                     serde_json::json!(donor_scores),
//                 )
//                 .ok();
//             handle
//                 .set_value(
//                     Box::new(IntronRetentionValue::RetentionAcceptorScore),
//                     serde_json::json!(acceptor_scores),
//                 )
//                 .ok();
//             handle
//                 .set_value(
//                     Box::new(IntronRetentionValue::NumberOfTrueRetentions),
//                     Value::Number(number_of_true_retentions.into()),
//                 )
//                 .ok();
//             handle
//                 .set_value(
//                     Box::new(IntronRetentionValue::NumberOfPartialRetentions),
//                     Value::Number(number_of_partial_retentions.into()),
//                 )
//                 .ok();
//             handle
//                 .set_value(
//                     Box::new(IntronRetentionValue::NumberOfFalseRetentions),
//                     Value::Number(number_of_false_retentions.into()),
//                 )
//                 .ok();
//         } else {
//             let line = query.line().to_owned();
//             pass.push(line);
//         }
//     }

//     // after classying reads, we check bucket frequencies
//     // of the total reads in the bucket, we consider the bucket
//     // to be dirty and return the reads for recovery if args.recover
//     let ratio = count / totals;
//     let mut is_dirty = false;
//     if recover {
//         if ratio >= RETENTION_RATIO_THRESHOLD {
//             // log::warn!("Bucket {:?} is dirty -> {}", refs.reads, ratio);
//             let dirt = tmp_dirt;
//             is_dirty = true;

//             let new_passes = recover_from_dirt(
//                 dirt,
//                 exon_owners,
//                 &introns_owned,
//                 &refs,
//                 &mut descriptor,
//                 ratio,
//                 toga_introns,
//                 counter,
//             );
//             new_passes.iter().for_each(|p| {
//                 hits.retain(|x| x != p);
//                 pass.push(p.to_owned());
//             });
//         } else {
//             for query in queries.iter() {
//                 let handle = descriptor.get_mut(&query.name).unwrap();

//                 handle
//                     .set_value(
//                         Box::new(IntronRetentionValue::ComponentRetentionRatio),
//                         serde_json::json!(ratio),
//                     )
//                     .ok();
//             }
//         }
//     }

//     // dbg!(&descriptor);

//     (hits, pass, blocks, descriptor, is_dirty)
// }

// pub fn recover_from_dirt(
//     dirt: Vec<Arc<&GenePred>>,
//     exon_owners: BTreeSet<(u64, u64)>,
//     introns_owned: &BTreeMap<(u64, u64), (f32, f32)>,
//     refs: &RefGenePred,
//     descriptor: &mut HashMap<String, Box<dyn ModuleMap>>,
//     ratio: f32,
//     toga_introns: &Option<&HashSet<(u64, u64)>>,
//     counter: &ParallelCounter,
// ) -> Vec<String> {
//     let mut local_passes = vec![];
//     let mut owner_ratios = HashMap::new();
//     let mut owned_ratios = HashMap::new();

//     let mut is_toga_contained = false;

//     // exons: 1) build ratios
//     for owner in exon_owners.iter() {
//         let mut background = 0_f32;
//         let mut count = 0_f32;
//         for read in refs.reads.iter() {
//             let ref_exons = &read.exons;
//             let (owner_start, owner_end) = owner;

//             if ref_exons.contains(&(*owner_start, *owner_end)) {
//                 count += 1.0;
//             }

//             // any read that spans the owner is considered background
//             if (read.start < *owner_start) && (*owner_end < read.end) {
//                 background += 1.0;
//             }
//         }

//         let ratio = count / background;
//         owner_ratios.insert(owner, ratio);
//     }

//     // introns: 1) build ratios and 2) remove introns with low ratios
//     let mut supported_introns = vec![];
//     for (owned, score) in introns_owned.iter() {
//         let mut background = 0.0;
//         let mut count = 0.0;
//         for read in refs.reads.iter() {
//             let ref_introns = read.get_introns();
//             let (owned_start, owned_end) = owned;

//             if ref_introns.contains(&(*owned_start, *owned_end)) {
//                 count += 1.0;
//             }

//             // any read that spans the owned intron is considered background
//             if (read.start < *owned_start) && (*owned_end < read.end) {
//                 background += 1.0;
//             }
//         }

//         let ratio = count / background;
//         owned_ratios.insert(owned, ratio);

//         let donor_score = score.0;
//         let acceptor_score = score.1;

//         is_toga_contained = toga_introns.map(|t| t.contains(owned)).unwrap_or(false);

//         // if the intron is supported by the majority of reads (>=50%) or
//         // the splice scores are high enough, we consider the intron to be
//         // supported and consider it a true intron
//         if (ratio >= INTRON_RETENTION_RECOVERY_THRESHOLD)
//             || (donor_score >= SPLICE_AI_SCORE_RECOVERY_THRESHOLD
//                 && acceptor_score >= SPLICE_AI_SCORE_RECOVERY_THRESHOLD)
//             || is_toga_contained
//         {
//             if donor_score >= 0.02 && acceptor_score >= 0.02 {
//                 counter.inc_true_retentions();
//             } else if donor_score >= 0.01 && donor_score >= 0.01 {
//                 counter.inc_partial_retentions();
//             } else {
//                 // if refs.strand == '+' {
//                 //     dbg!(
//                 //         &dirt[0].chrom,
//                 //         owned.0,
//                 //         owned.1,
//                 //         donor_score,
//                 //         acceptor_score,
//                 //         ratio
//                 //     )
//                 // } else {
//                 //     dbg!(
//                 //         &dirt[0].chrom,
//                 //         SCALE - owned.1,
//                 //         SCALE - owned.0,
//                 //         donor_score,
//                 //         acceptor_score,
//                 //         ratio
//                 //     )
//                 // };
//                 counter.inc_false_retentions();
//             }

//             supported_introns.push(owned);
//         }
//     }

//     for read in dirt.clone().iter_mut() {
//         let handle = descriptor.get_mut(&read.name).unwrap();

//         handle
//             .set_value(
//                 Box::new(IntronRetentionValue::IsDirtyComponent),
//                 Value::Bool(true),
//             )
//             .ok();
//         handle
//             .set_value(
//                 Box::new(IntronRetentionValue::ComponentRetentionRatio),
//                 serde_json::json!(ratio),
//             )
//             .ok();

//         let mut number_of_unsupported_retentions = 0;
//         let mut number_of_supported_retentions = 0;
//         let mut number_of_read_false_supported_retentions = 0;
//         let mut number_of_read_partial_supported_retentions = 0;
//         let mut number_of_read_true_supported_retentions = 0;

//         let mut intron_support_ratio = vec![];
//         let mut exon_support_ratio = vec![];
//         let mut is_retention_supported_map = vec![];

//         for exon in &read.exons {
//             for (intron, score) in introns_owned.iter() {
//                 if exon.1 < intron.0 {
//                     break;
//                 }

//                 if exon.0 < intron.0 && intron.1 < exon.1 {
//                     intron_support_ratio.push(
//                         owned_ratios
//                             .get(intron)
//                             .expect("ERROR: Intron not found in intron_owned_ratios"),
//                     );
//                     exon_support_ratio.push(
//                         owner_ratios
//                             .get(exon)
//                             .expect("ERROR: Exon not found in exon_owner_ratios"),
//                     );

//                     if supported_introns.contains(&intron) {
//                         // exon eats a 'likely' real intron, fill up the following descriptors:
//                         // unrecover,
//                         // is_retention_supported_map
//                         number_of_supported_retentions += 1;
//                         is_retention_supported_map.push(true);

//                         let donor_score = score.0;
//                         let acceptor_score = score.1;

//                         if donor_score >= 0.02 && acceptor_score >= 0.02 {
//                             number_of_read_true_supported_retentions += 1;
//                         } else if donor_score >= 0.01 && donor_score >= 0.01 {
//                             number_of_read_partial_supported_retentions += 1;
//                         } else {
//                             number_of_read_false_supported_retentions += 1;
//                         }
//                     } else {
//                         // exon eats a 'likely' false intron, fill up the following descriptors:
//                         // is_retention_supported_map
//                         number_of_unsupported_retentions += 1;
//                         is_retention_supported_map.push(false);
//                     }
//                 }
//             }
//         }

//         if number_of_supported_retentions < 1
//             && number_of_unsupported_retentions + number_of_supported_retentions > 0
//         {
//             local_passes.push(read.line().to_owned());

//             handle
//                 .set_value(
//                     Box::new(IntronRetentionValue::IsRetentionSupported),
//                     serde_json::json!(false),
//                 )
//                 .ok();
//         } else {
//             handle
//                 .set_value(
//                     Box::new(IntronRetentionValue::IsRetentionSupported),
//                     serde_json::json!(true),
//                 )
//                 .ok();
//         }

//         handle
//             .set_value(
//                 Box::new(IntronRetentionValue::NumberOfUnrecovers),
//                 serde_json::json!(number_of_supported_retentions),
//             )
//             .ok();
//         handle
//             .set_value(
//                 Box::new(IntronRetentionValue::NumberOfRecovers),
//                 serde_json::json!(number_of_unsupported_retentions),
//             )
//             .ok();
//         handle
//             .set_value(
//                 Box::new(IntronRetentionValue::IntronSupportRatio),
//                 serde_json::json!(intron_support_ratio),
//             )
//             .ok();
//         handle
//             .set_value(
//                 Box::new(IntronRetentionValue::ExonSupportRatio),
//                 serde_json::json!(exon_support_ratio),
//             )
//             .ok();
//         handle
//             .set_value(
//                 Box::new(IntronRetentionValue::IsRetentionSupportedMap),
//                 serde_json::json!(is_retention_supported_map),
//             )
//             .ok();
//         handle
//             .set_value(
//                 Box::new(IntronRetentionValue::NumberOfFalseRententionsSupported),
//                 serde_json::json!(number_of_read_false_supported_retentions),
//             )
//             .ok();
//         handle
//             .set_value(
//                 Box::new(IntronRetentionValue::NumberOfTrueRententionsSupported),
//                 serde_json::json!(number_of_read_true_supported_retentions),
//             )
//             .ok();
//         handle
//             .set_value(
//                 Box::new(IntronRetentionValue::NumberOfPartialRententionsSupported),
//                 serde_json::json!(number_of_read_partial_supported_retentions),
//             )
//             .ok();
//         handle
//             .set_value(
//                 Box::new(IntronRetentionValue::IsTogaIntron),
//                 serde_json::json!(is_toga_contained),
//             )
//             .ok();
//     }

//     local_passes
// }

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use std::io::Write;
//     use tempfile;

//     #[test]
//     fn test_detection_fn_with_tempfile() {
//         let mut file = tempfile::NamedTempFile::new().unwrap();
//         let path = file.path().to_path_buf();
//         write!(
//             file,
//             "s1/t778845/t804002/tm54164U_210309_085211/65275776/ccs_PerID0.996_5Clip0_3Clip0_PolyA274_PolyARead275/t60/t-/t778845/t804002/t255,0,0/t14/t1268,142,88,200,253,203,254,167,120,142,218,197,132,302/t0,14396,14736,14943,16461,16846,17598,18073,18511,18890,20330,21290,22627,24855"
//         ).unwrap();

//         write!(
//             file,
//             "s1/t778870/t803968/tm54164U_210309_085211/92276372/ccs_PerID1.000_5Clip0_3Clip0_PolyA71_PolyARead72/t60/t-/t778870/t803968/t255,0,0/t11/t1243,142,88,200,253,1006,167,218,197,132,268/t0,14371,14711,14918,16436,16821,18048,20305,21265,22602,24830"
//         ).unwrap();

//         let args = Args {
//             refs: [path.clone()].to_vec(),
//             query: [path].to_vec(),
//             threads: 1,
//             blacklist: Vec::new(),
//             plot: false,
//             recover: false,
//             splice_scores: Some(std::path::PathBuf::new()),
//             toga: None,
//         };

//         assert!(detect_intron_retentions(args).is_ok());
//     }
// }
