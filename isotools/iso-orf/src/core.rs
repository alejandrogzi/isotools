use anyhow::Result;
use config::{get_progress_bar, write_objs, OverlapType, ORF_ASSIGNED_READS};
use dashmap::DashSet;
use hashbrown::{HashMap, HashSet};
use iso_utr::utils::unpack_blacklist;
use packbed::{buckerize, combine, record::Bed6, unpack, GenePred, PackMode};
use rayon::prelude::*;

use std::sync::atomic::{AtomicU32, Ordering};

use crate::cli::Args;

pub fn call_orfs(args: Args) -> Result<()> {
    let raw = unpack::<GenePred, _>(args.raw, OverlapType::Exon, true)
        .expect("ERROR: Failed to unpack raw tracks");
    let calls = unpack::<Bed6, _>(args.calls, OverlapType::Exon, false)
        .expect("ERROR: Failed to unpack calls");

    let (tracks, n) = combine(calls, raw);
    let buckets = buckerize(tracks, OverlapType::Exon, n, PackMode::Paired); // WARN: overlap does not matter
    let blacklist = unpack_blacklist(args.blacklist).unwrap_or_default();

    let accumulator = ParallelAccumulator::default();
    let counter = ParallelCounter::default();
    let pb = get_progress_bar(buckets.len() as u64, "Processing...");

    buckets.into_par_iter().for_each(|bucket| {
        let components = bucket.1;
        counter.inc_components(components.len() as u32);

        components.into_par_iter().for_each(|mut comp| {
            let comp = comp
                .as_any_mut()
                .downcast_mut::<(Vec<GenePred>, Vec<GenePred>)>()
                .expect("ERROR: Failed to downcast to GenePred");

            let (assigned, unassigned) = process_component(comp, &blacklist);

            assigned.into_iter().for_each(|read| {
                accumulator.orfs.insert(read);
            });

            unassigned.into_iter().for_each(|read| {
                accumulator.unassigned.insert(read);
                counter.inc_unassgined();
            });
        });

        pb.inc(1);
    });

    pb.finish_and_clear();
    write_objs(&accumulator.orfs, ORF_ASSIGNED_READS);

    Ok(())
}

pub struct ParallelCounter {
    pub unassigned: AtomicU32,
    pub components: AtomicU32,
}

impl ParallelCounter {
    fn new() -> Self {
        Self {
            unassigned: AtomicU32::new(0),
            components: AtomicU32::new(0),
        }
    }

    pub fn inc_components(&self, count: u32) {
        self.components.fetch_add(count, Ordering::Relaxed);
    }

    pub fn inc_unassgined(&self) {
        self.unassigned.fetch_add(1, Ordering::Relaxed);
    }

    pub fn get_counters(&self) -> (f64, f64) {
        (
            self.unassigned.load(Ordering::Relaxed) as f64,
            self.components.load(Ordering::Relaxed) as f64,
        )
    }

    pub fn get_stat(&self) -> (f64, f64) {
        let (unassigned, components) = self.get_counters();
        (unassigned, (unassigned / components) * 100.0)
    }
}

impl Default for ParallelCounter {
    fn default() -> Self {
        Self::new()
    }
}

pub struct ParallelAccumulator {
    pub orfs: DashSet<String>,
    pub unassigned: DashSet<String>,
}

impl Default for ParallelAccumulator {
    fn default() -> Self {
        Self {
            orfs: DashSet::new(),
            unassigned: DashSet::new(),
        }
    }
}

fn process_component(
    component: &mut (Vec<GenePred>, Vec<GenePred>),
    banned: &HashSet<String>,
) -> (Vec<String>, Vec<String>) {
    let (mut assigned, mut unassigned) = (Vec::new(), Vec::new());
    let (raw_reads, orf_calls) = component;

    let orfs = orf_calls
        .iter()
        .map(|orf| {
            let name = orf.name.split('_').collect::<Vec<_>>();
            let name = name[..name.len() - 2].join("_");

            (name, (orf.start, orf.end))
        })
        .collect::<HashMap<_, _>>();

    raw_reads.iter_mut().for_each(|read| {
        if banned.contains(&read.name) {
            return;
        }

        if let Some((start, end)) = orfs.get(&read.name) {
            // INFO: mutate read information
            assigned.push(read.insert_utr(*start, *end).construct_bed_line());
        } else {
            unassigned.push(read.line.clone());
        }
    });

    (assigned, unassigned)
}
