use config::{
    bed_to_nested_map, get_progress_bar, write_objs, BedColumn, BedColumnValue,
    INTRAPRIMING_RATIO_THRESHOLD, INTRAPRIMING_REVIEW, POLYA_INTRAPRIMING, POLYA_PASS,
};
use dashmap::{DashMap, DashSet};
use hashbrown::{HashMap, HashSet};
use packbed::{
    packbed, reader,
    record::{Bed6, PolyAPred},
    BedPackage, PackMode,
};
use rayon::prelude::*;

use std::{
    path::PathBuf,
    sync::{
        atomic::{AtomicU32, Ordering},
        Arc,
    },
};

use crate::cli::CallerArgs;

pub fn pas_caller(args: CallerArgs) -> Result<(), Box<dyn std::error::Error>> {
    let isoseqs = packbed(
        vec![args.bed],
        None,
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

        if let Some(aparent) = aparent_scores.get(&chr) {
            distribute(
                components,
                &*aparent,
                &accumulator,
                &counter,
                args.recover,
                args.wiggle,
                args.max_gpa_length,
                args.min_polya_length,
                args.aparent_threshold,
            );
        } else {
            log::warn!("No APARENT scores found for chromosome: {}", chr);

            distribute(
                components,
                &HashMap::new(),
                &accumulator,
                &counter,
                args.recover,
                args.wiggle,
                args.max_gpa_length,
                args.min_polya_length,
                args.aparent_threshold,
            );
        }

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

        if accumulator.review.len() > 0 {
            write_objs(&accumulator.review, INTRAPRIMING_REVIEW);
        }
    }

    [&accumulator.pass, &accumulator.intrapriming]
        .par_iter()
        .zip([POLYA_PASS, POLYA_INTRAPRIMING].par_iter())
        .for_each(|(rx, path)| write_objs(&rx, path));

    Ok(())
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
    let aparent_scores = bed_to_nested_map::<Bed6>(Arc::new(content), BedColumn::Score)
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
) {
    components.into_par_iter().for_each(|mut comp| {
        let comp = comp
            .as_any_mut()
            .downcast_mut::<Vec<PolyAPred>>()
            .expect("ERROR: Could not downcast to PolyAPred!");

        let rs = process_component(
            comp,
            scores,
            recover,
            wiggle,
            max_gpa_length,
            min_polya_length,
            aparent_threshold,
        );

        if rs.2.is_some() {
            counter.inc_review(1);
        }

        accumulator.add(rs);
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
}

impl Default for ParallelAccumulator {
    fn default() -> Self {
        Self {
            pass: DashSet::new(),
            intrapriming: DashSet::new(),
            review: DashSet::new(),
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
    fn add(&self, result: (Vec<String>, Vec<String>, Option<Vec<String>>)) {
        let (passes, intrapriming, review) = result;

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
    component: &mut Vec<PolyAPred>,
    scores: &HashMap<String, BedColumnValue>,
    recover: bool,
    wiggle: usize,
    max_gpa_length: usize,
    min_polya_length: usize,
    aparent_threshold: f32,
) -> (Vec<String>, Vec<String>, Option<Vec<String>>) {
    let mut seen_positions = HashSet::new();
    component.sort_by(|a, b| a.end.cmp(&b.end));

    let mut passes = Vec::new();
    let mut intrapraming = Vec::new();
    let mut review = Vec::new();

    let (mut intrapriming_count, totals) = (0_f32, component.len() as f32);

    // dbg!(&component);

    component.iter().for_each(|read| {
        let aparent_score = if let Some(score) = scores.get(&read.name) {
            let score = score
                .max_score()
                .expect("ERROR: Could not get aparent score as float!");

            score
        } else {
            0.0 // WARN: default to 0.0 if no score is found
        };

        if any_true!(
            read.gpa < max_gpa_length as u32,
            read.poly_a > min_polya_length as u32,
            seen_positions.contains(&read.end),
            seen_positions.contains(&(read.end - wiggle as u64)),
            seen_positions.contains(&(read.end + wiggle as u64)),
            aparent_score > aparent_threshold
        ) {
            passes.push(read.line.clone());
            seen_positions.insert(read.end); // INFO: only considering 'good' reads
        } else {
            intrapriming_count += 1.0;
            intrapraming.push(read.line.clone());
        }
    });

    let ratio = intrapriming_count / totals;

    if recover {
        if ratio > INTRAPRIMING_RATIO_THRESHOLD {
            let p = passes.drain(..);
            let i = intrapraming.drain(..);

            review.extend(p);
            review.extend(i);

            return (passes, intrapraming, Some(review));
        }
    }

    return (passes, intrapraming, None);
}
