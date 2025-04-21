use std::path::PathBuf;
use std::sync::atomic::{AtomicU32, Ordering};

use dashmap::{DashMap, DashSet};
use hashbrown::HashSet;
use packbed::par_reader;
use rayon::iter::ParallelIterator;
use rayon::prelude::*;
use rayon::str::ParallelString;

use config::{write_objs, ModuleMap, TRUNCATIONS, TRUNCATION_FREE};

pub fn unpack_blacklist<'a>(paths: Vec<PathBuf>) -> Option<HashSet<String>> {
    if paths.is_empty() {
        return None;
    }

    let contents = par_reader(paths).unwrap();
    let tracks = contents
        .par_lines()
        .filter_map(|line| {
            if line.is_empty() {
                return None;
            };

            Some(line.to_string())
        })
        .collect::<HashSet<String>>();

    Some(tracks)
}

pub struct ParallelCounter {
    pub dirties: AtomicU32,
    pub components: AtomicU32,
}

impl ParallelCounter {
    fn new() -> Self {
        Self {
            dirties: AtomicU32::new(0),
            components: AtomicU32::new(0),
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
}

impl Default for ParallelCounter {
    fn default() -> Self {
        Self::new()
    }
}

pub struct ParallelAccumulator {
    pub truncations: DashSet<String>,
    pub no_truncations: DashSet<String>,
    pub miscellaneous: DashSet<String>,
    pub descriptor: DashMap<String, Box<dyn ModuleMap>>,
}

impl Default for ParallelAccumulator {
    fn default() -> Self {
        Self {
            truncations: DashSet::new(),
            no_truncations: DashSet::new(),
            miscellaneous: DashSet::new(),
            descriptor: DashMap::new(),
        }
    }
}

impl ParallelAccumulator {
    pub fn num_truncations(&self) -> usize {
        self.truncations.len()
    }
}

pub fn write_results(accumulator: &ParallelAccumulator) {
    [&accumulator.truncations, &accumulator.no_truncations]
        .par_iter()
        .zip([TRUNCATIONS, TRUNCATION_FREE].par_iter())
        .for_each(|(acc, path)| write_objs(acc, path));
}
