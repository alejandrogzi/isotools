use dashmap::DashSet;
use hashbrown::{HashMap, HashSet};
use iso_classify::utils::{bed_to_map, Bed4};
use packbed::par_reader;
use rayon::prelude::*;

use std::path::PathBuf;
use std::sync::Arc;

use std::sync::atomic::{AtomicU32, Ordering};

use config::{write_objs, CoordType, BED3, INTRON_RETENTIONS, INTRON_RETENTION_FREE};

pub struct ParallelCounter {
    pub dirties: AtomicU32,
    pub components: AtomicU32,
    pub true_retentions: AtomicU32,
    pub false_retentions: AtomicU32,
    pub partial_retentions: AtomicU32,
}

impl ParallelCounter {
    fn new() -> Self {
        Self {
            dirties: AtomicU32::new(0),
            components: AtomicU32::new(0),
            true_retentions: AtomicU32::new(0),
            false_retentions: AtomicU32::new(0),
            partial_retentions: AtomicU32::new(0),
        }
    }

    pub fn inc_components(&self, count: u32) {
        self.components.fetch_add(count, Ordering::Relaxed);
    }

    pub fn inc_dirty(&self) {
        self.dirties.fetch_add(1, Ordering::Relaxed);
    }

    pub fn get_counters(&self) -> (f64, f64) {
        (
            self.dirties.load(Ordering::Relaxed) as f64,
            self.components.load(Ordering::Relaxed) as f64,
        )
    }

    pub fn get_stat(&self) -> (f64, f64) {
        let (dirties, components) = self.get_counters();
        (dirties, (dirties / components) * 100.0)
    }

    pub fn inc_true_retentions(&self) {
        self.true_retentions.fetch_add(1, Ordering::Relaxed);
    }

    pub fn inc_false_retentions(&self) {
        self.false_retentions.fetch_add(1, Ordering::Relaxed);
    }

    pub fn inc_partial_retentions(&self) {
        self.partial_retentions.fetch_add(1, Ordering::Relaxed);
    }

    fn load_retentions(&self) -> (u32, u32, u32) {
        (
            self.true_retentions.load(Ordering::Relaxed),
            self.false_retentions.load(Ordering::Relaxed),
            self.partial_retentions.load(Ordering::Relaxed),
        )
    }
}

impl Default for ParallelCounter {
    fn default() -> Self {
        Self::new()
    }
}

pub struct ParallelAccumulator {
    pub retentions: DashSet<String>,
    pub non_retentions: DashSet<String>,
    pub miscellaneous: DashSet<String>,
}

impl Default for ParallelAccumulator {
    fn default() -> Self {
        Self {
            retentions: DashSet::new(),
            non_retentions: DashSet::new(),
            miscellaneous: DashSet::new(),
        }
    }
}

impl ParallelAccumulator {
    pub fn num_retentions(&self) -> usize {
        self.retentions.len()
    }
}

pub fn unpack_blacklist<'a>(paths: Vec<PathBuf>) -> Option<HashMap<String, HashSet<(u64, u64)>>> {
    if paths.is_empty() {
        return None;
    }

    let contents = Arc::new(par_reader(paths).unwrap());
    let tracks = bed_to_map::<Bed4>(contents, CoordType::Bounds).unwrap();

    Some(tracks)
}

pub fn write_results(accumulator: &ParallelAccumulator, plot: bool) {
    [&accumulator.retentions, &accumulator.non_retentions]
        .par_iter()
        .zip([INTRON_RETENTIONS, INTRON_RETENTION_FREE].par_iter())
        .for_each(|(acc, path)| write_objs(acc, path));

    if plot {
        write_objs(&accumulator.miscellaneous, BED3);
    }
}
