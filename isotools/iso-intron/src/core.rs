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

use std::collections::BTreeMap;
use std::path::PathBuf;

use anyhow::Result;
use hashbrown::{HashMap, HashSet};
use log::info;
use packbed::record::IntronPosition;
use packbed::{packbed, BedPackage, GenePred, IntronPred};
use rayon::prelude::*;
use serde_json::Value;

use crate::cli::Args;
use crate::utils::{unpack_blacklist, ParallelAccumulator, ParallelCounter};

use config::{
    get_progress_bar, par_write_results, write_descriptor, IntronRetentionValue, ModuleDescriptor,
    ModuleMap, ModuleType, OverlapType, INTRON_RETENTIONS, INTRON_RETENTION_DESCRIPTOR,
    INTRON_RETENTION_FREE, INTRON_RETENTION_REVIEW, RETENTION_RATIO_THRESHOLD, SCALE,
};

/// Detects intron retentions in a query set of reads
///
/// # Arguments
///
/// * `args` - The command line arguments
///
/// # Returns
///
/// * `Result<()>` - The result of the operation
///
/// # Example
///
/// ```rust, no_run
/// let args = Args::new();
/// detect_intron_retentions(args).unwrap();
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

        process_components(components, banned, &accumulator, &counter, args.recover);

        pb.inc(1);
    });

    pb.finish_and_clear();
    info!(
        "Reads with retained introns: {}",
        accumulator.num_retentions()
    );

    write_descriptor(&accumulator.descriptor, INTRON_RETENTION_DESCRIPTOR);
    par_write_results(
        accumulator,
        vec![
            PathBuf::from(INTRON_RETENTIONS),
            PathBuf::from(INTRON_RETENTION_FREE),
            PathBuf::from(INTRON_RETENTION_REVIEW),
        ],
        None,
    );

    Ok(())
}

/// Processes the components of reads and introns in parallel
///
/// # Arguments
///
/// * `components` - The components to process per chromosome
/// * `banned` - the set of banned introns for the chromosome
/// * `accumulator` - the accumulator to use
/// * `counter` - The counter to use
/// * `recover` - Indicates if the component should be recovered
///
/// # Returns
///
/// * None
///
/// # Example
///
/// ```rust, no_run
/// let components = vec![];
/// let banned = HashSet::new();
/// let accumulator = ParallelAccumulator::default();
/// let counter = ParallelCounter::default();
///
/// process_components(components, &banned, &accumulator, &counter, false);
///
/// assert_eq!(accumulator.num_retentions(), 0);
/// assert_eq!(counter.num_components(), 0);
/// ```
#[inline(always)]
fn process_components(
    components: Vec<Box<dyn BedPackage>>,
    banned: &HashSet<(u64, u64)>,
    accumulator: &ParallelAccumulator,
    counter: &ParallelCounter,
    recover: bool,
) {
    components.into_par_iter().for_each(|mut comp| {
        let comp = comp
            .as_any_mut()
            .downcast_mut::<(Vec<IntronPred>, Vec<GenePred>)>()
            .expect("ERROR: Could not downcast to IntronPred and GenePred!");

        // if comp is len 1 OR comp is len <=5 and no TOGA, continue
        // if comp.1.len() <= 5 {
        //     counter.inc_skipped();
        //     return;
        // }

        let (keep, discard, review, descriptor) = process_component(comp, banned, counter, recover);
        accumulator.add(keep, discard, review, descriptor);
    });
}

/// IntronModuleReadAction enum
///
/// This enum is used to determine the action to take for a read
///
/// * Discard: The read has an RT intron and should be discarded
/// * Keep: The read does not have an RT intron or retains anything and should be kept
/// * Unclear: The read has an unclear status and should be reviewed
///
/// # Example
///
/// ```rust, no_run
/// let action = IntronModuleReadAction::Discard;
///
/// match action {
///     IntronModuleReadAction::Discard => println!("Discarding read"),
///     IntronModuleReadAction::Keep => println!("Keeping read"),
///     IntronModuleReadAction::Unclear => println!("Unclear read"),
/// }
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum IntronModuleReadAction {
    // Review,  // Unclear
    Discard, // RT-intron
    Keep,
    Unclear, // Artifact
}

/// Display trait for IntronModuleReadAction
///
/// This trait is used to display the action as a string
///
/// # Example
///
/// ```rust, no_run
/// let action = IntronModuleReadAction::Discard;
///
/// println!("{}", action);
/// ```
impl std::fmt::Display for IntronModuleReadAction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            IntronModuleReadAction::Discard => write!(f, "DISCARD"),
            IntronModuleReadAction::Keep => write!(f, "KEEP"),
            IntronModuleReadAction::Unclear => write!(f, "UNCLEAR"),
        }
    }
}

/// Processes a component of reads and introns
///
/// # Arguments
///
/// * `comp` - The component to process
/// * `ban` - The set of banned introns
/// * `counter` - The counter to use
/// * `recover` - Indicates if the component should be recovered
///
/// # Returns
///
/// * `Vec<String>` - The vector of reads to keep
/// * `Vec<String>` - The vector of reads to discard
/// * `Option<Vec<String>>` - The vector of reads to review
/// * `HashMap<String, Box<dyn ModuleMap>>` - The descriptor to fill up
///
/// # Example
///
/// ```rust, no_run
/// let mut comp = (vec![], vec![]);
/// let ban = HashSet::new();
/// let counter = ParallelCounter::default();
/// let recover = false;
///
/// let (keep, discard, review, descriptor) = process_component(&mut comp, &ban, &counter, recover);
///
/// assert_eq!(keep.len(), 0);
/// assert_eq!(discard.len(), 0);
/// assert_eq!(review, None);
/// assert_eq!(descriptor.len(), 0);
/// ```
#[inline(always)]
pub fn process_component(
    comp: &mut (Vec<IntronPred>, Vec<GenePred>),
    ban: &HashSet<(u64, u64)>,
    counter: &ParallelCounter,
    recover: bool,
) -> (
    Vec<String>,
    Vec<String>,
    Option<Vec<String>>,
    HashMap<String, Box<dyn ModuleMap>>,
) {
    let mut descriptor = HashMap::new();

    let mut keep = Vec::new();
    let mut discard = Vec::new();

    let introns = &comp.0;
    let reads = &comp.1;

    let (mut count, totals) = (0_f32, reads.len() as f32);

    // INFO: convert Vec<IntronPred> into BTreeMap<(u64, u64), IntronPred>
    // WARN: conserving order is important to avoid early false breaks!
    let ref_introns = introns
        .iter()
        .map(|intron| ((intron.start, intron.end), intron))
        .collect::<BTreeMap<_, _>>();

    // INFO: for every read -> see which introns are retained and which the read has!
    for read in reads {
        let mut schema = RetentionSchema::default();

        schema.ref_introns_component_size = Value::Number(introns.len().into());
        schema.query_component_size = Value::Number(reads.len().into());

        detect_rt_intron(read, ban, &ref_introns, &mut schema);
        detect_retention(read, &ref_introns, &mut schema);

        if schema.intron_retention {
            count += 1.0;
        }

        match (schema.intronic_status, schema.exonic_status) {
            (IntronModuleReadAction::Discard, _) | (_, IntronModuleReadAction::Discard) => {
                counter.inc_retentions();
                discard.push(read.line().to_owned());
            }
            _ => {
                keep.push(read.line().to_owned());
            }
        };

        schema.difuse(&mut descriptor, read);
    }

    let ratio = count / totals;

    // INFO: second pass to fill up the descriptor with comp ratio
    descriptor.iter_mut().for_each(|(_, handle)| {
        handle
            .set_value(
                Box::new(IntronRetentionValue::ComponentRetentionRatio),
                serde_json::json!(ratio),
            )
            .ok();
    });

    if recover {
        if ratio > RETENTION_RATIO_THRESHOLD {
            let review = recover_component(&mut keep, &mut discard, reads, &mut descriptor);
            return (keep, discard, Some(review), descriptor);
        }
    }

    (keep, discard, None, descriptor)
}

/// Recovers the component by moving reads from discard to review
/// if the retention ratio is above the threshold
///
/// # Arguments
///
/// * `keep` - The vector of reads to keep
/// * `discard` - The vector of reads to discard
/// * `reads` - The vector of reads to process
/// * `descriptor` - The descriptor to fill up
///
/// # Returns
///
/// * `Vec<String>` - The vector of reads to review
///
/// # Example
///
/// ```rust, no_run
/// let mut keep = vec![];
/// let mut discard = vec![];
/// let reads = vec![];
/// let mut descriptor = HashMap::new();
///
/// let review = recover_component(&mut keep, &mut discard, &reads, &mut descriptor);
///
/// assert_eq!(review.len(), 0);
/// ```
fn recover_component(
    keep: &mut Vec<String>,
    discard: &mut Vec<String>,
    reads: &Vec<GenePred>,
    descriptor: &mut HashMap<String, Box<dyn ModuleMap>>,
) -> Vec<String> {
    let mut review = vec![];

    let k = keep.drain(..);
    let d = discard.drain(..);

    review.extend(k);
    review.extend(d);

    reads.iter().for_each(|read| {
        let handle = descriptor.get_mut(&read.name).unwrap();

        handle
            .set_value(
                Box::new(IntronRetentionValue::IsDirtyComponent),
                serde_json::json!(true),
            )
            .ok();
    });

    return review;
}

/// RetentionSchema struct
///
/// This struct is used to store the retention schema for a read.
/// Mimics IntronRetentionValue enum and fields of IntronModuleDescriptor
///
/// # Fields
///
/// * `intron_retention`: bool - Indicates if the read has an intron retention
/// * `retention_support_type`: Vec<Value> - The support type of the retention
/// * `number_of_retentions`: Value - The number of retentions
/// * `coords_of_retention`: Vec<Value> - The coordinates of the retention
/// * `location_of_retention`: Vec<Value> - The location of the retention
/// * `is_intron_retained_in_frame`: Vec<Value> - Indicates if the retention is in frame
/// * `retains_rt_intron`: Value - Indicates if the retention retains an RT intron
/// * `retains_rt_map`: Vec<Value> - The map of the retained RT intron
/// * `has_rt_intron`: Value - Indicates if the read has an RT intron
/// * `has_rt_intron_map`: Vec<Value> - The map of the RT intron
/// * `retention_acceptor_score`: Vec<Value> - The acceptor score of the retention
/// * `retention_donor_score`: Vec<Value> - The donor score of the retention
/// * `ref_introns_component_size`: Value - The size of the reference introns component
/// * `query_component_size`: Value - The size of the query component
/// * `component_retention_ratio`: Value - The retention ratio of the component
/// * `is_dirty_component`: Value - Indicates if the component is dirty
/// * `exonic_status`: IntronModuleReadAction - The exonic status of the read
/// * `intronic_status`: IntronModuleReadAction - The intronic status of the read
struct RetentionSchema {
    pub intron_retention: bool,
    pub retention_support_type: Vec<Value>,
    pub number_of_retentions: Value,
    pub coords_of_retention: Vec<Value>,
    pub location_of_retention: Vec<Value>,
    pub is_intron_retained_in_frame: Vec<Value>,
    pub retains_rt_intron: Value,
    pub retains_rt_map: Vec<Value>,
    pub has_rt_intron: Value,
    pub has_rt_intron_map: Vec<Value>,
    pub retention_acceptor_score: Vec<Value>,
    pub retention_donor_score: Vec<Value>,
    pub ref_introns_component_size: Value,
    pub query_component_size: Value,
    pub component_retention_ratio: Value,
    pub is_dirty_component: Value,
    pub exonic_status: IntronModuleReadAction,
    pub intronic_status: IntronModuleReadAction,
}

impl RetentionSchema {
    /// Fills up the descriptor with the values from the schema
    ///
    /// # Arguments
    ///
    /// * `descriptor` - The descriptor to fill up
    /// * `read` - The read to fill up the descriptor with
    ///
    /// # Returns
    ///
    /// * None
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let mut descriptor = HashMap::new();
    /// let read = GenePred::new();
    ///
    /// schema.difuse(&mut descriptor, &read);
    /// ```
    fn difuse(&self, descriptor: &mut HashMap<String, Box<dyn ModuleMap>>, read: &GenePred) {
        descriptor.insert(
            read.name.clone(),
            ModuleDescriptor::with_schema(ModuleType::IntronRetention),
        );
        let handle = descriptor.get_mut(&read.name).unwrap();

        for (field, value) in self.as_vals() {
            handle.set_value(Box::new(field), value).ok();
        }
    }

    /// Converts the schema into an iterable collection
    ///
    /// # Returns
    ///
    /// * Vec<(IntronRetentionValue, Value)> - The vector of tuples
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let schema = RetentionSchema::default();
    ///
    /// let vals = schema.as_vals();
    ///
    /// for (field, value) in vals {
    ///     println!("{}: {}", field, value);
    /// }
    /// ```
    pub fn as_vals(&self) -> Vec<(IntronRetentionValue, Value)> {
        vec![
            (
                IntronRetentionValue::IsIntronRetention,
                serde_json::json!(self.intron_retention),
            ),
            (
                IntronRetentionValue::RetentionSupportType,
                Value::Array(self.retention_support_type.clone()),
            ),
            (
                IntronRetentionValue::NumberOfRetentions,
                self.number_of_retentions.clone(),
            ),
            (
                IntronRetentionValue::RetentionCoords,
                Value::Array(self.coords_of_retention.clone()),
            ),
            (
                IntronRetentionValue::RetentionLocation,
                Value::Array(self.location_of_retention.clone()),
            ),
            (
                IntronRetentionValue::IsIntronRetainedInFrame,
                Value::Array(self.is_intron_retained_in_frame.clone()),
            ),
            (
                IntronRetentionValue::RetainsRtIntron,
                self.retains_rt_intron.clone(),
            ),
            (
                IntronRetentionValue::RetainsRtIntronMap,
                Value::Array(self.retains_rt_map.clone()),
            ),
            (
                IntronRetentionValue::HasRTIntron,
                self.has_rt_intron.clone(),
            ),
            (
                IntronRetentionValue::HasRTIntronMap,
                Value::Array(self.has_rt_intron_map.clone()),
            ),
            (
                IntronRetentionValue::RetentionAcceptorScore,
                Value::Array(self.retention_acceptor_score.clone()),
            ),
            (
                IntronRetentionValue::RetentionDonorScore,
                Value::Array(self.retention_donor_score.clone()),
            ),
            (
                IntronRetentionValue::RefIntronsComponentSize,
                self.ref_introns_component_size.clone(),
            ),
            (
                IntronRetentionValue::QueryComponentSize,
                self.query_component_size.clone(),
            ),
            (
                IntronRetentionValue::ComponentRetentionRatio,
                self.component_retention_ratio.clone(),
            ),
            (
                IntronRetentionValue::IsDirtyComponent,
                self.is_dirty_component.clone(),
            ),
            (
                IntronRetentionValue::ExonicStatus,
                Value::String(self.exonic_status.to_string()),
            ),
            (
                IntronRetentionValue::IntronicStatus,
                Value::String(self.intronic_status.to_string()),
            ),
        ]
    }
}

/// Implements the Default trait for RetentionSchema
///
/// This trait is used to create a default instance of the RetentionSchema struct
///
/// # Example
///
/// ```rust, no_run
/// let schema = RetentionSchema::default();
///
/// assert_eq!(schema.intron_retention, false);
/// assert_eq!(schema.retention_support_type.len(), 0);
/// ```
impl Default for RetentionSchema {
    fn default() -> Self {
        RetentionSchema {
            intron_retention: false,
            retention_support_type: vec![],
            number_of_retentions: Value::Null,
            coords_of_retention: vec![],
            location_of_retention: vec![],
            is_intron_retained_in_frame: vec![],
            retains_rt_intron: Value::Null,
            retains_rt_map: vec![],
            has_rt_intron: Value::Null,
            has_rt_intron_map: vec![],
            retention_acceptor_score: vec![],
            retention_donor_score: vec![],
            ref_introns_component_size: Value::Null,
            query_component_size: Value::Null,
            component_retention_ratio: Value::Null,
            is_dirty_component: Value::Bool(false),
            exonic_status: IntronModuleReadAction::Keep,
            intronic_status: IntronModuleReadAction::Keep,
        }
    }
}

/// Detects if a read has any RT introns
///
/// # Arguments
///
/// * `read` - The read to check for RT introns
/// * `handle` - The handle to set the values in
/// * `ban` - The set of banned introns
/// * `ref_introns` - The reference introns
///
/// # Returns
///
/// * `IntronModuleReadAction` - The action to take for the read
///
/// # Example
///
/// ```rust, no_run
/// let status = detect_rt_intron(&read, &mut handle, &ban, &ref_introns);
///
/// match status {
///   IntronModuleReadAction::Discard => println!("Discarding read"),
///   IntronModuleReadAction::Keep => println!("Keeping read"),
///   IntronModuleReadAction::Unclear => println!("Unclear read"),
/// }
/// ```
fn detect_rt_intron(
    read: &GenePred,
    ban: &HashSet<(u64, u64)>,
    ref_introns: &BTreeMap<(u64, u64), &IntronPred>,
    schema: &mut RetentionSchema,
) {
    // INFO: Determine if read has RT introns -> limit to only first hit
    // INFO: will not count ALL RT introns but just the first one!
    let read_introns = read.get_introns();
    let mut action = IntronModuleReadAction::Keep;

    for intron in read_introns {
        if ban.contains(&(intron.0, intron.1)) {
            continue;
        }

        let hit = ref_introns
            .get(&intron)
            .expect("ERROR: Intron not found, this is likely a bug!");

        match hit.stats.support {
            config::SupportType::RT => {
                // INFO: if read is RT, discard
                action = IntronModuleReadAction::Discard;
                schema.has_rt_intron = Value::Bool(true);

                match read.strand {
                    config::Strand::Forward => {
                        let coord = format!("{}:{}-{}", read.chrom, intron.0, intron.1);
                        schema.has_rt_intron_map.push(Value::String(coord));
                    }
                    config::Strand::Reverse => {
                        let coord =
                            format!("{}:{}-{}", read.chrom, SCALE - intron.1, SCALE - intron.0);
                        schema.has_rt_intron_map.push(Value::String(coord));
                    }
                }
            }
            config::SupportType::Splicing | config::SupportType::Unclear => {} // INFO: Do nothing
        };
    }

    if action == IntronModuleReadAction::Keep {
        schema.has_rt_intron = Value::Bool(false);
    }

    schema.intronic_status = action;
}

/// Determine if any read exons retain any intron
///
/// # Arguments
///
/// * `read` - The read to check for intron retention
/// * `ref_introns` - The reference introns
/// * `schema` - The schema to set the values in
///
/// # Returns
///
/// * None -> the schema is modified in place
///
/// # Example
///
/// ```rust, no_run
/// let mut schema = RetentionSchema::default();
/// detect_retention(&read, &ref_introns, &mut schema);
///
/// match schema.intron_retention {
///     true => println!("Intron retention detected"),
///     false => println!("No intron retention detected"),
/// }
/// ```
fn detect_retention(
    read: &GenePred,
    ref_introns: &BTreeMap<(u64, u64), &IntronPred>,
    schema: &mut RetentionSchema,
) {
    let mut action = IntronModuleReadAction::Keep;
    let read_exons = read.get_exons();

    for exon in read_exons {
        let exon_start = exon.0;
        let exon_end = exon.1;

        for (intron, stats) in ref_introns.iter() {
            let intron_start = intron.0;
            let intron_end = intron.1;

            // INFO: early exit to avoid unnecessary checks -> enforces BTreeMap instead of HashMap!
            if intron_end < exon_start {
                continue;
            }

            if exon_start < intron_start && intron_end < exon_end {
                match stats.stats.support {
                    config::SupportType::RT => {
                        // INFO: retaining a false intron is not an IR, do nothing
                        match read.strand {
                            config::Strand::Forward => {
                                let coord =
                                    format!("{}:{}-{}", read.chrom, intron_start, intron_end);
                                schema.retains_rt_map.push(Value::String(coord));
                            }
                            config::Strand::Reverse => {
                                let coord = format!(
                                    "{}:{}-{}",
                                    read.chrom,
                                    SCALE - intron_end,
                                    SCALE - intron_start
                                );
                                schema.retains_rt_map.push(Value::String(coord));
                            }
                        }

                        schema.retains_rt_intron = Value::Bool(true);
                    }
                    config::SupportType::Unclear | config::SupportType::Splicing => {
                        // INFO: any other variant is considered an intron retention!
                        if !schema.intron_retention {
                            schema.intron_retention = true;
                        }

                        schema
                            .retention_support_type
                            .push(Value::String(stats.stats.support.to_string()));
                        schema
                            .is_intron_retained_in_frame
                            .push(Value::Bool(stats.stats.is_in_frame));
                        schema
                            .retention_donor_score
                            .push(serde_json::json!(stats.stats.splice_ai_donor));
                        schema
                            .retention_acceptor_score
                            .push(serde_json::json!(stats.stats.splice_ai_acceptor));

                        match read.strand {
                            config::Strand::Forward => {
                                let coord =
                                    format!("{}:{}-{}", read.chrom, intron_start, intron_end);
                                schema.coords_of_retention.push(Value::String(coord));
                            }
                            config::Strand::Reverse => {
                                let coord = format!(
                                    "{}:{}-{}",
                                    read.chrom,
                                    SCALE - intron_end,
                                    SCALE - intron_start
                                );
                                schema.coords_of_retention.push(Value::String(coord));
                            }
                        }

                        match stats.stats.intron_position {
                            IntronPosition::CDS => {
                                // INFO: if intron in frame keep read --> splice variant producing a longer protein
                                // INFO: else do not annotate read (exon structure is flawed)
                                if !stats.stats.is_in_frame {
                                    if stats.stats.support == config::SupportType::Unclear {
                                        // INFO: if support is unclear, keep read
                                        // INFO: edge case where IR is in CDS + not in-frame
                                        // INFO: but we do not know if intron is real or not
                                        action = IntronModuleReadAction::Unclear;
                                    } else {
                                        // INFO: clear point where IR is in CDS + not in-frame + splicing
                                        // INFO: if support is splicing, discard read
                                        action = IntronModuleReadAction::Discard;
                                    }
                                }

                                schema
                                    .location_of_retention
                                    .push(Value::String(stats.stats.intron_position.to_string()));
                            }
                            // INFO: keep read --> variant would not affect CDS
                            IntronPosition::UTR | IntronPosition::Mixed => {
                                schema
                                    .location_of_retention
                                    .push(Value::String(stats.stats.intron_position.to_string()));
                            }
                            IntronPosition::Unknown => {
                                schema
                                    .location_of_retention
                                    .push(Value::String(stats.stats.intron_position.to_string()));

                                // INFO: logic changed here -> if intron is not in frame, discard read
                                // INFO: if intron is in frame, keep read -> flag it as UNCLEAR
                                if !stats.stats.is_in_frame {
                                    action = IntronModuleReadAction::Discard;
                                } else {
                                    action = IntronModuleReadAction::Unclear;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    schema.number_of_retentions = Value::Number(schema.coords_of_retention.len().into());
    schema.exonic_status = action;
}
