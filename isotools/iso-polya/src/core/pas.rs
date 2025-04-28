use config::{
    bed_to_nested_map, get_progress_bar, write_descriptor, write_objs, BedColumn, BedColumnValue,
    ModuleDescriptor, ModuleMap, ModuleType, PolyAPredictionValue, INTRAPRIMING_RATIO_THRESHOLD,
    INTRAPRIMING_REVIEW, POLYA_DESCRIPTOR, POLYA_INTRAPRIMING, POLYA_PASS,
};
use dashmap::{DashMap, DashSet};
use hashbrown::{HashMap, HashSet};
use packbed::{
    packbed, reader,
    record::{Bed6, PolyAPred},
    BedPackage, GenePred, PackMode,
};
use rayon::prelude::*;

lazy_static::lazy_static! {
    static ref _HASH: HashMap<String, BedColumnValue> = HashMap::new();
}

use std::{
    path::PathBuf,
    sync::{
        atomic::{AtomicU32, Ordering},
        Arc,
    },
};

use crate::cli::CallerArgs;

pub fn pas_caller(
    args: CallerArgs,
) -> Result<DashMap<String, Box<dyn ModuleMap>>, Box<dyn std::error::Error>> {
    log::info!("INFO: Starting PAS caller...");

    let isoseqs = packbed(
        vec![args.bed.clone()],
        args.toga,
        config::OverlapType::Exon,
        PackMode::PolyA,
    )?;

    // WARN: duplicated aparent predictions are taken as max
    let aparent_scores = get_aparent_scores(args.aparent);

    let pb = get_progress_bar(isoseqs.len() as u64, "Processing reads...");
    let accumulator = ParallelAccumulator::default();
    let counter = ParallelCounter::default();

    isoseqs.into_par_iter().for_each(|(chr, components)| {
        counter.inc_comp(components.len() as u32);

        let aparent_ref = aparent_scores.get(&chr);
        let scores = match &aparent_ref {
            Some(r) => &**r,
            None => &_HASH,
        };

        distribute(
            components,
            scores,
            &accumulator,
            &counter,
            args.recover,
            args.wiggle,
            args.max_gpa_length,
            args.min_polya_length,
            args.aparent_threshold,
            args.filter,
            args.filter_side,
        );

        pb.inc(1);
    });

    pb.finish_and_clear();
    log::info!(
        "Detected potential intraprimings: {}",
        accumulator.intrapriming.len()
    );
    log::info!("Intrapriming-free reads: {}", accumulator.pass.len());

    if args.recover {
        log::warn!(
            "Number of dirty components in query reads: {:?} ({:.3}%)",
            counter.num_of_reviews,
            counter.load_ratio()
        );

        if !args.in_memory {
            if accumulator.review.len() > 0 {
                write_objs(&accumulator.review, INTRAPRIMING_REVIEW);
            }
        }
    }

    if !args.in_memory {
        log::info!("INFO: Writing results to disk...");

        [&accumulator.pass, &accumulator.intrapriming]
            .par_iter()
            .zip(
                [
                    args.outdir
                        .join(POLYA_PASS)
                        .to_str()
                        .expect("ERROR: Could not convert path!"),
                    args.outdir
                        .join(POLYA_INTRAPRIMING)
                        .to_str()
                        .expect("ERROR: Could not convert path!"),
                ]
                .par_iter(),
            )
            .for_each(|(rx, path)| write_objs(&rx, path));

        write_descriptor(&accumulator.descriptor, POLYA_DESCRIPTOR);
    }

    Ok(accumulator.descriptor)
}

/// Reads the APARENT file and returns a nested map
///
/// # Arguments
///
/// * `file` - Path to the APARENT file
///
/// # Returns
///
/// * `DashMap<String, HashMap<String, BedColumnValue>>`
///    Nested map with chromosome as key and a map of read
///    name and aparent score as value
///
/// # Example
///
/// ```rust, no_run
/// let aparent_scores = get_aparent_scores(PathBuf::from("path/to/aparent.bed"));
///
/// assert_eq!(aparent_scores.len(), 25);
/// ```
fn get_aparent_scores(file: PathBuf) -> DashMap<String, HashMap<String, BedColumnValue>> {
    let content = reader(file).expect("ERROR: Could not read aparent file!");
    let aparent_scores =
        bed_to_nested_map::<Bed6>(Arc::new(content), BedColumn::End, BedColumn::Score)
            .expect("ERROR: Could not build mapper from aparent file!");

    aparent_scores
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
/// ```rust, no_run
/// distribute(args);
/// ```
#[inline(always)]
fn distribute(
    components: Vec<Box<dyn BedPackage>>,
    scores: &HashMap<String, BedColumnValue>,
    accumulator: &ParallelAccumulator,
    counter: &ParallelCounter,
    recover: bool,
    wiggle: usize,
    max_gpa_length: usize,
    min_polya_length: usize,
    aparent_threshold: f32,
    filter: bool,
    filter_side: config::FilterSide,
) {
    components.into_par_iter().for_each(|mut comp| {
        let comp = comp
            .as_any_mut()
            .downcast_mut::<(Vec<PolyAPred>, Vec<GenePred>)>()
            .expect("ERROR: Could not downcast to PolyAPred!");

        if !filter {
            let result = process_component(
                comp,
                scores,
                recover,
                wiggle,
                max_gpa_length,
                min_polya_length,
                aparent_threshold,
            );

            if result.2.is_some() {
                counter.inc_review(1);
            }

            accumulator.add(result);
        } else {
            let rs = filter_component(
                comp,
                scores,
                min_polya_length,
                aparent_threshold,
                filter_side,
            );

            accumulator.add_filter(rs);
        }
    });
}

/// Parallel accumulator for the processing function
///
/// # Fields
///
/// * `pass` - DashSet to store the pass reads
/// * `intrapriming` - DashSet to store the intrapriming reads
/// * `review` - DashSet to store the review reads
///
/// # Example
///
/// ```rust, no_run
/// let accumulator = ParallelAccumulator::default();
///
/// assert_eq!(accumulator.pass.len(), 0);
/// ```
struct ParallelAccumulator {
    pass: DashSet<String>,
    intrapriming: DashSet<String>,
    review: DashSet<String>,
    descriptor: DashMap<String, Box<dyn ModuleMap>>,
}

impl Default for ParallelAccumulator {
    fn default() -> Self {
        Self {
            pass: DashSet::new(),
            intrapriming: DashSet::new(),
            review: DashSet::new(),
            descriptor: DashMap::new(),
        }
    }
}

impl ParallelAccumulator {
    /// Adds a result to the accumulator
    ///
    /// # Arguments
    ///
    /// * `result` - Tuple with the pass, intrapriming and review reads
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let mut accumulator = ParallelAccumulator::default();
    /// accumulator.add((vec!["read1", "read2"], vec!["read3"], Some(vec!["read4"])));
    ///
    /// assert_eq!(accumulator.pass.len(), 2);
    /// ```
    fn add(
        &self,
        result: (
            Vec<String>,
            Vec<String>,
            Option<Vec<String>>,
            HashMap<String, Box<dyn ModuleMap>>,
        ),
    ) {
        let (passes, intrapriming, review, descriptor) = result;

        intrapriming.into_iter().for_each(|ipm| {
            self.intrapriming.insert(ipm);
        });
        passes.into_iter().for_each(|pass| {
            self.pass.insert(pass);
        });
        if let Some(r) = review {
            r.into_iter().for_each(|r| {
                self.review.insert(r);
            });
        }

        descriptor.into_iter().for_each(|(k, v)| {
            self.descriptor.insert(k, v);
        });
    }

    fn add_filter(&self, result: Vec<String>) {
        result.into_iter().for_each(|r| {
            self.pass.insert(r);
        });
    }
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
/// ```rust, no_run
/// let counter = ParallelCounter::default();
///
/// assert_eq!(counter.num_of_comps.load(Ordering::Relaxed), 0);
/// ```
struct ParallelCounter {
    num_of_comps: AtomicU32,
    num_of_reviews: AtomicU32,
}

impl Default for ParallelCounter {
    fn default() -> Self {
        Self {
            num_of_comps: AtomicU32::new(0),
            num_of_reviews: AtomicU32::new(0),
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
    /// ```rust, no_run
    /// let counter = ParallelCounter::default();
    /// counter.inc_comp(1);
    ///
    /// assert_eq!(counter.num_of_comps.load(Ordering::Relaxed), 1);
    /// ```
    fn inc_comp(&self, value: u32) {
        self.num_of_comps.fetch_add(value, Ordering::Relaxed);
    }

    /// Increment the review counter
    ///
    /// # Arguments
    ///
    /// * `value` - Value to increment the counter
    ///
    /// # Example
    ///
    /// ```rust, no_run
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
    /// ```rust, no_run
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
}

/// Macro to check if any of the conditions are true
///
/// # Example
///
/// ```rust, no_run
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
    component: &mut (Vec<PolyAPred>, Vec<GenePred>),
    scores: &HashMap<String, BedColumnValue>,
    recover: bool,
    wiggle: usize,
    max_gpa_length: usize,
    min_polya_length: usize,
    aparent_threshold: f32,
) -> (
    Vec<String>,
    Vec<String>,
    Option<Vec<String>>,
    HashMap<String, Box<dyn ModuleMap>>,
) {
    let mut descriptor = HashMap::new();

    let reads = &mut component.0;
    let toga = &mut component.1;

    let mut seen_positions = HashSet::new();
    let mut seen_reads = HashSet::new();

    let mut passes = Vec::new();
    let mut intrapraming = Vec::new();
    let mut review = Vec::new();

    let (mut intrapriming_count, totals) = (0_f32, reads.len() as f32);

    // dbg!(&component);

    reads.iter().for_each(|read| {
        let end_of_read = match read.strand {
            config::Strand::Forward => (read.end).to_string(),
            config::Strand::Reverse => (read.end + 1).to_string(),
        };

        let aparent_score = if let Some(score) = scores.get(&end_of_read) {
            let score = score
                .max_score()
                .expect("ERROR: Could not get aparent score as float!");

            score
        } else {
            0.0 // WARN: default to 0.0 if no score is found
        };

        // INFO: filling up descriptor for each read
        descriptor.insert(
            read.name.clone(),
            ModuleDescriptor::with_schema(ModuleType::PolyAPrediction),
        );
        let handle = descriptor.get_mut(&read.name).unwrap();

        handle
            .set_value(
                Box::new(PolyAPredictionValue::PolyAScore),
                serde_json::json!(aparent_score),
            )
            .ok();
        handle
            .set_value(
                Box::new(PolyAPredictionValue::GenomicPolyA),
                serde_json::json!(read.gpa),
            )
            .ok();
        handle
            .set_value(
                Box::new(PolyAPredictionValue::WholePolyALength),
                serde_json::json!(read.poly_a),
            )
            .ok();

        let location = get_read_toga_location(read.end, toga, handle);

        if any_true!(
            read.gpa < max_gpa_length as u32,
            read.poly_a >= min_polya_length as u32,
            seen_positions.contains(&read.end),
            seen_positions.contains(&(read.end - wiggle as u64)),
            seen_positions.contains(&(read.end + wiggle as u64)),
            aparent_score > aparent_threshold,
            location == PolyATogaLocation::UTR
        ) {
            passes.push(read.line.clone());

            seen_positions.insert(read.end); // INFO: only considering 'good' reads
            seen_positions.insert(read.end - wiggle as u64); // INFO: inserting wiggle room!
            seen_positions.insert(read.end + wiggle as u64);

            seen_reads.insert(&read.name); // INFO: storing 'good' reads for later escape
        }
    });

    // INFO: necessary second iteration to check for intrapriming events
    // INFO: since the only filter to doble check is the wiggle + genomic pos
    // INFO: we can iterate over the component again and avoid the rest of steps
    reads.iter().for_each(|read| {
        let handle = descriptor.get_mut(&read.name).unwrap();

        if seen_reads.contains(&&read.name) {
            handle
                .set_value(
                    Box::new(PolyAPredictionValue::IsIntrapriming),
                    serde_json::json!(false),
                )
                .ok();
            handle
                .set_value(
                    Box::new(PolyAPredictionValue::IsPolyASupported),
                    serde_json::json!(true),
                )
                .ok();
            handle
                .set_value(
                    Box::new(PolyAPredictionValue::ForcedPolyAPass),
                    serde_json::json!(false),
                )
                .ok();

            return;
        }

        if any_true!(
            seen_positions.contains(&read.end),
            seen_positions.contains(&(read.end - wiggle as u64)),
            seen_positions.contains(&(read.end + wiggle as u64))
        ) {
            passes.push(read.line.clone());

            seen_positions.insert(read.end); // INFO: only considering 'good' reads
            seen_positions.insert(read.end - wiggle as u64); // INFO: inserting wiggle room!
            seen_positions.insert(read.end + wiggle as u64);

            handle
                .set_value(
                    Box::new(PolyAPredictionValue::IsIntrapriming),
                    serde_json::json!(false),
                )
                .ok();

            // INFO: set to false because the evidence is based on other reads!
            handle
                .set_value(
                    Box::new(PolyAPredictionValue::IsPolyASupported),
                    serde_json::json!(false),
                )
                .ok();
            handle
                .set_value(
                    Box::new(PolyAPredictionValue::ForcedPolyAPass),
                    serde_json::json!(true),
                )
                .ok();
        } else {
            intrapriming_count += 1.0;
            intrapraming.push(read.line.clone());

            handle
                .set_value(
                    Box::new(PolyAPredictionValue::IsIntrapriming),
                    serde_json::json!(true),
                )
                .ok();
            handle
                .set_value(
                    Box::new(PolyAPredictionValue::IsPolyASupported),
                    serde_json::json!(false),
                )
                .ok();
        }
    });

    let ratio = intrapriming_count / totals;

    if recover {
        if ratio > INTRAPRIMING_RATIO_THRESHOLD {
            let p = passes.drain(..);
            let i = intrapraming.drain(..);

            review.extend(p);
            review.extend(i);

            reads.iter().for_each(|read| {
                let handle = descriptor.get_mut(&read.name).unwrap();

                handle
                    .set_value(
                        Box::new(PolyAPredictionValue::IsDirtyPolyAComponent),
                        serde_json::json!(true),
                    )
                    .ok();
                handle
                    .set_value(
                        Box::new(PolyAPredictionValue::IntraprimingComponentRatio),
                        serde_json::json!(ratio),
                    )
                    .ok();
            });

            return (passes, intrapraming, Some(review), descriptor);
        }
    }

    return (passes, intrapraming, None, descriptor);
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
/// ```rust, no_run
/// let read_end = 100;
/// let mut toga = vec![GenePred::default()];
/// let mut handle = ModuleMap::default();
///
/// get_read_toga_location(read_end, &mut toga, &mut handle);
/// ```
fn get_read_toga_location(
    read_end: u64,
    toga: &mut Vec<GenePred>,
    handle: &mut Box<dyn ModuleMap>,
) -> PolyATogaLocation {
    let location = PolyATogaLocation::CDS;

    for projection in toga.iter() {
        if projection.start < read_end && read_end < projection.end {
            handle
                .set_value(
                    Box::new(PolyAPredictionValue::PolyALocation),
                    serde_json::json!("CDS"),
                )
                .ok();

            return location;
        }
    }

    // INFO: only if we didn’t return early
    handle
        .set_value(
            Box::new(PolyAPredictionValue::PolyALocation),
            serde_json::json!("UTR"),
        )
        .ok();

    PolyATogaLocation::UTR
}

/// Enum to represent the location of the polyA tail
///
/// # Example
///
/// ```rust, no_run
/// let location = PolyATogaLocation::CDS;
///
/// assert_eq!(location, PolyATogaLocation::CDS);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
enum PolyATogaLocation {
    CDS,
    UTR,
}

/// Filter components based on polyA tail length and APARENT score
///
/// # Arguments
///
/// * `component` - Vector of PolyAPred to process
/// * `scores` - Nested map with chromosome as key and a map of read name and aparent score as value
/// * `min_polya_length` - PolyA tail length threshold [min length allowed]
/// * `aparent_threshold` - APARENT threshold [min score allowed]
///
/// # Returns
///
/// * `Vec<String>` - Vector of reads that passed the filter
///
/// # Example
///
/// ```rust, no_run
/// let mut component = vec![PolyAPred::default()];
/// let scores = HashMap::new();
/// let min_polya_length = 10;
/// let aparent_threshold = 0.5;
/// let side = FilterSide::Above;
///
/// let filtered = filter_component(&mut component, &scores, min_polya_length, aparent_threshold, side);
///
/// assert_eq!(filtered.len(), 0);
/// ```
fn filter_component(
    component: &mut (Vec<PolyAPred>, Vec<GenePred>),
    scores: &HashMap<String, BedColumnValue>,
    min_polya_length: usize,
    aparent_threshold: f32,
    side: config::FilterSide,
) -> Vec<String> {
    let reads = &mut component.0;
    let _ = &component.1;

    reads.sort_by(|a, b| a.end.cmp(&b.end));

    let mut above = Vec::new();
    let mut below = Vec::new();

    reads.iter().for_each(|read| {
        let aparent_score = if let Some(score) = scores.get(&read.name) {
            let score = score
                .max_score()
                .expect("ERROR: Could not get aparent score as float!");

            score
        } else {
            0.0 // WARN: default to 0.0 if no score is found
        };

        if any_true!(
            read.poly_a >= min_polya_length as u32,
            aparent_score > aparent_threshold
        ) {
            above.push(read.line.clone());
        } else {
            below.push(read.line.clone());
        }
    });

    match side {
        config::FilterSide::Above => return above,
        config::FilterSide::Below => return below,
    }
}
