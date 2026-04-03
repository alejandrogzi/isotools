//! Core module for detecting intron retentions in a query set of reads
//! Alejandro Gonzales-Irribarren, 2026
//!
//! This module contains the main algorithm for detecting truncations
//! in a query set of reads.
//!
//! In short it takes a set of query reads and a set of reference reads and
//! detects truncations in the query reads. It does this by checking if the
//! query reads overlap with any middle exon from the reference set of reads.
//! A recovery step can be performed by evaluating the support of the middle
//! exon in the reference set of reads.

use dashmap::DashSet;
use std::sync::atomic::{AtomicU32, Ordering};

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

    pub fn load_components(&self) -> u32 {
        self.components.load(Ordering::Relaxed)
    }
}

impl Default for ParallelCounter {
    fn default() -> Self {
        Self::new()
    }
}

pub struct ParallelAccumulator {
    pub lines: DashSet<Vec<u8>>,
}

impl Default for ParallelAccumulator {
    fn default() -> Self {
        Self {
            lines: DashSet::new(),
        }
    }
}
