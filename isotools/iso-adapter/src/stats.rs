// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Thread-safe run-time statistics for iso-adapter.
//!
//! `ClipStats` lives in an `Arc` and is updated from every worker without
//! locking the whole structure — only the individual per-adapter and
//! per-bin counters use atomics. A short `report()` emits an INFO-level
//! table at the end of the run.

use std::collections::HashMap;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Mutex;

use log::info;

use crate::detector::{AdapterMatch, ClipEnd};
use crate::logging::elapsed;

/// Number of log2-spaced bins: covers clip lengths up to 2^(NUM_BINS-1).
pub const NUM_LENGTH_BINS: usize = 16;

/// Aggregate statistics for an iso-adapter run.
pub struct ClipStats {
    /// Total BAM records inspected.
    records_total: AtomicU64,
    /// Records that were unmapped / secondary / supplementary and therefore
    /// skipped for clip analysis.
    records_skipped: AtomicU64,
    /// Total 5' soft-clips seen (including zero-length: they are filtered out
    /// upstream, so this only counts non-zero).
    clips_5p: AtomicU64,
    /// Total 3' soft-clips seen.
    clips_3p: AtomicU64,
    /// Records that had at least one adapter detected anywhere.
    records_with_adapter: AtomicU64,
    /// Records that had a clip but no adapter could be assigned.
    records_with_unknown_clip: AtomicU64,

    /// log2-binned histogram of clip lengths. `length_hist[i]` counts clips
    /// whose length falls in `[2^i, 2^(i+1))`.
    length_hist: [AtomicU64; NUM_LENGTH_BINS],

    /// Per-adapter-label hit counters. A Mutex here is fine because inserts
    /// only happen the first time a label is seen; after that the hot path
    /// is a single atomic fetch_add.
    per_adapter: Mutex<HashMap<&'static str, AtomicU64>>,

    /// Fuzzy-match count (edit distance ≥ 1).
    fuzzy_hits: AtomicU64,
    /// Exact-match count (edit distance == 0).
    exact_hits: AtomicU64,

    /// Records successfully trimmed (only incremented when -R is used).
    records_trimmed: AtomicU64,
}

impl Default for ClipStats {
    fn default() -> Self {
        Self {
            records_total: AtomicU64::new(0),
            records_skipped: AtomicU64::new(0),
            clips_5p: AtomicU64::new(0),
            clips_3p: AtomicU64::new(0),
            records_with_adapter: AtomicU64::new(0),
            records_with_unknown_clip: AtomicU64::new(0),
            length_hist: std::array::from_fn(|_| AtomicU64::new(0)),
            per_adapter: Mutex::new(HashMap::new()),
            fuzzy_hits: AtomicU64::new(0),
            exact_hits: AtomicU64::new(0),
            records_trimmed: AtomicU64::new(0),
        }
    }
}

impl ClipStats {
    /// Records one BAM record having been seen.
    pub fn record_seen(&self) {
        self.records_total.fetch_add(1, Ordering::Relaxed);
    }

    /// Records a skipped BAM record (unmapped / secondary / etc.).
    pub fn record_skipped(&self) {
        self.records_skipped.fetch_add(1, Ordering::Relaxed);
    }

    /// Records one observed clip of a given length and end.
    pub fn observe_clip(&self, end: ClipEnd, len: usize) {
        if len == 0 {
            return;
        }
        match end {
            ClipEnd::FivePrime => &self.clips_5p,
            ClipEnd::ThreePrime => &self.clips_3p,
        }
        .fetch_add(1, Ordering::Relaxed);

        let bin = log2_bin(len);
        self.length_hist[bin].fetch_add(1, Ordering::Relaxed);
    }

    /// Records a successful adapter hit.
    pub fn observe_match(&self, m: &AdapterMatch) {
        if m.edit_distance == 0 {
            self.exact_hits.fetch_add(1, Ordering::Relaxed);
        } else {
            self.fuzzy_hits.fetch_add(1, Ordering::Relaxed);
        }

        let mut map = self
            .per_adapter
            .lock()
            .expect("per_adapter mutex poisoned");
        map.entry(m.label)
            .or_insert_with(|| AtomicU64::new(0))
            .fetch_add(1, Ordering::Relaxed);
    }

    /// Records that a record had any adapter detected.
    pub fn record_with_adapter(&self) {
        self.records_with_adapter.fetch_add(1, Ordering::Relaxed);
    }

    /// Records that a record had an unknown clip (no DB hit).
    pub fn record_with_unknown_clip(&self) {
        self.records_with_unknown_clip
            .fetch_add(1, Ordering::Relaxed);
    }

    /// Records that a record was trimmed in the output BAM.
    pub fn record_trimmed(&self) {
        self.records_trimmed.fetch_add(1, Ordering::Relaxed);
    }

    /// Total records seen.
    pub fn total(&self) -> u64 {
        self.records_total.load(Ordering::Relaxed)
    }

    /// Total records with at least one adapter hit.
    pub fn matched(&self) -> u64 {
        self.records_with_adapter.load(Ordering::Relaxed)
    }

    /// Emits a structured INFO log table at end of run.
    pub fn report(&self, novel_candidates: &[(Vec<u8>, u32)]) {
        let total = self.records_total.load(Ordering::Relaxed);
        let skipped = self.records_skipped.load(Ordering::Relaxed);
        let with_adapter = self.records_with_adapter.load(Ordering::Relaxed);
        let with_unknown = self.records_with_unknown_clip.load(Ordering::Relaxed);
        let trimmed = self.records_trimmed.load(Ordering::Relaxed);
        let exact = self.exact_hits.load(Ordering::Relaxed);
        let fuzzy = self.fuzzy_hits.load(Ordering::Relaxed);
        let c5 = self.clips_5p.load(Ordering::Relaxed);
        let c3 = self.clips_3p.load(Ordering::Relaxed);

        info!(
            "[{}] summary records={total} skipped={skipped} with_adapter={with_adapter} \
             with_unknown_clip={with_unknown} trimmed={trimmed}",
            elapsed()
        );
        info!(
            "[{}] clips 5p={c5} 3p={c3} exact_hits={exact} fuzzy_hits={fuzzy}",
            elapsed(),
        );

        info!("[{}] clip length distribution (log2 bins):", elapsed());
        for (i, bin) in self.length_hist.iter().enumerate() {
            let count = bin.load(Ordering::Relaxed);
            if count == 0 {
                continue;
            }
            let lo = 1u64 << i;
            let hi = 1u64 << (i + 1);
            info!(
                "[{}]   len[{lo:>5}..{hi:>5}) = {count}",
                elapsed()
            );
        }

        let map = self
            .per_adapter
            .lock()
            .expect("per_adapter mutex poisoned");
        let mut ordered: Vec<(&&'static str, u64)> = map
            .iter()
            .map(|(k, v)| (k, v.load(Ordering::Relaxed)))
            .collect();
        ordered.sort_unstable_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(b.0)));

        info!("[{}] per-adapter hits:", elapsed());
        for (label, count) in &ordered {
            info!("[{}]   {label:<40} {count}", elapsed());
        }

        if !novel_candidates.is_empty() {
            info!(
                "[{}] novel-adapter candidates (≥ threshold, top 20):",
                elapsed()
            );
            for (clip, count) in novel_candidates.iter().take(20) {
                let display = String::from_utf8_lossy(clip);
                info!("[{}]   {count:>6}× {display}", elapsed());
            }
        }
    }
}

/// Returns the log2 bin index for a clip length, clamped to the last bin.
fn log2_bin(len: usize) -> usize {
    if len == 0 {
        return 0;
    }
    let floor = (usize::BITS - 1 - len.leading_zeros()) as usize;
    floor.min(NUM_LENGTH_BINS - 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log2_bin_monotone() {
        assert_eq!(log2_bin(1), 0);
        assert_eq!(log2_bin(2), 1);
        assert_eq!(log2_bin(3), 1);
        assert_eq!(log2_bin(4), 2);
        assert_eq!(log2_bin(255), 7);
        assert_eq!(log2_bin(usize::MAX), NUM_LENGTH_BINS - 1);
    }

    #[test]
    fn stats_counts() {
        let stats = ClipStats::default();
        stats.record_seen();
        stats.record_seen();
        stats.observe_clip(ClipEnd::FivePrime, 12);
        stats.observe_clip(ClipEnd::ThreePrime, 200);
        stats.record_with_adapter();
        assert_eq!(stats.total(), 2);
        assert_eq!(stats.matched(), 1);
        assert_eq!(stats.clips_5p.load(Ordering::Relaxed), 1);
        assert_eq!(stats.clips_3p.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn stats_records_per_adapter() {
        use crate::detector::AdapterKind;
        let stats = ClipStats::default();
        let m = AdapterMatch {
            label: "pacbio:isoseq:primer_5p_tso",
            clip_range: 0..25,
            edit_distance: 0,
            on_reverse_strand: false,
            kind: AdapterKind::Adapter,
        };
        stats.observe_match(&m);
        stats.observe_match(&m);
        let map = stats.per_adapter.lock().unwrap();
        assert_eq!(
            map.get("pacbio:isoseq:primer_5p_tso")
                .unwrap()
                .load(Ordering::Relaxed),
            2
        );
    }
}
