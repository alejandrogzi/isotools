// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Soft-clip extraction and adapter matching.
//!
//! The detector is built once per run from the static `ADAPTER_DB` and then
//! shared across worker threads. Each worker extracts 5' / 3' soft-clip
//! byte slices straight out of `noodles_sam::alignment::RecordBuf` and
//! looks them up. Matching is two-stage:
//!
//! 1. **Exact**: Aho–Corasick over all adapters and their reverse complements.
//!    O(n) in the clip length, returns the leftmost longest match.
//! 2. **Fuzzy fallback**: `triple_accel::levenshtein_search` with a bounded
//!    edit distance. Only invoked when the clip is at least `min_clip_len`
//!    bases long and the exact pass came up empty.
//!
//! Unknown clips (those that fail both stages) flow into `UnknownClipStore`,
//! a sharded concurrent counter that surfaces candidate novel adapters at
//! the end of the run.

use std::sync::Arc;

use aho_corasick::{AhoCorasick, AhoCorasickKind, MatchKind};
use dashmap::DashMap;
use noodles_sam::alignment::record::cigar::{op::Kind, Op};
use noodles_sam::alignment::RecordBuf;
use triple_accel::levenshtein::levenshtein_search_simd_with_opts;
use triple_accel::{Match as EditMatch, SearchType};

use crate::adapters::{ADAPTER_DB, MIN_ADAPTER_LEN};
use crate::error::AdapterError;

/// Which end of the read a clip came from.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum ClipEnd {
    FivePrime,
    ThreePrime,
}

/// A successful adapter hit inside a clip.
#[derive(Clone, Debug)]
pub struct AdapterMatch {
    /// Human-readable label from `ADAPTER_DB`.
    pub label: &'static str,
    /// Offset range within the clip (half-open, [start, end)).
    pub clip_range: std::ops::Range<usize>,
    /// Edit distance of the hit; 0 for exact matches.
    pub edit_distance: u32,
    /// Whether the match was found as the stored sequence or its reverse
    /// complement.
    pub on_reverse_strand: bool,
}

/// Immutable adapter matcher shared across threads.
pub struct AdapterDb {
    automaton: AhoCorasick,
    /// Parallel to `automaton` pattern IDs: `(label, is_reverse_complement,
    /// pattern_bytes)`. Stored sequences are kept so the fuzzy fallback can
    /// reuse them.
    entries: Vec<DbEntry>,
    /// Minimum clip length for fuzzy matching to run at all.
    min_clip_len: usize,
    /// Maximum edit distance permitted for a fuzzy hit.
    max_edit_dist: u32,
}

#[derive(Clone, Debug)]
struct DbEntry {
    label: &'static str,
    sequence: Vec<u8>,
    reverse_complement: bool,
}

impl AdapterDb {
    /// Builds the matcher from the static `ADAPTER_DB`.
    pub fn new(min_clip_len: usize, max_edit_dist: u32) -> Result<Self, AdapterError> {
        let mut entries: Vec<DbEntry> = Vec::with_capacity(ADAPTER_DB.len() * 2);
        let mut patterns: Vec<Vec<u8>> = Vec::with_capacity(ADAPTER_DB.len() * 2);

        for (seq, label) in ADAPTER_DB {
            if seq.len() < MIN_ADAPTER_LEN {
                continue;
            }
            let forward = seq.to_vec();
            patterns.push(forward.clone());
            entries.push(DbEntry {
                label,
                sequence: forward,
                reverse_complement: false,
            });

            let rc = reverse_complement(seq);
            if rc != *seq {
                patterns.push(rc.clone());
                entries.push(DbEntry {
                    label,
                    sequence: rc,
                    reverse_complement: true,
                });
            }
        }

        let automaton = AhoCorasick::builder()
            .match_kind(MatchKind::LeftmostLongest)
            .kind(Some(AhoCorasickKind::DFA))
            .ascii_case_insensitive(false)
            .build(&patterns)
            .map_err(|e| AdapterError::Matcher(e.to_string()))?;

        Ok(Self {
            automaton,
            entries,
            min_clip_len,
            max_edit_dist,
        })
    }

    /// Returns the number of enrolled patterns (including RCs).
    pub fn pattern_count(&self) -> usize {
        self.entries.len()
    }

    /// Attempts to locate an adapter in `clip`.
    ///
    /// Exact matches are preferred; when none exist the fuzzy fallback runs
    /// only if the clip is long enough and `max_edit_dist > 0`.
    pub fn match_adapter(&self, clip: &[u8]) -> Option<AdapterMatch> {
        if clip.is_empty() {
            return None;
        }

        if let Some(m) = self.automaton.find(clip) {
            let entry = &self.entries[m.pattern().as_usize()];
            return Some(AdapterMatch {
                label: entry.label,
                clip_range: m.start()..m.end(),
                edit_distance: 0,
                on_reverse_strand: entry.reverse_complement,
            });
        }

        if self.max_edit_dist == 0 || clip.len() < self.min_clip_len {
            return None;
        }

        let mut best: Option<(EditMatch, &DbEntry)> = None;
        for entry in &self.entries {
            if entry.sequence.len() > clip.len() + self.max_edit_dist as usize {
                continue;
            }
            let hits = levenshtein_search_simd_with_opts(
                &entry.sequence,
                clip,
                self.max_edit_dist,
                SearchType::Best,
                triple_accel::levenshtein::LEVENSHTEIN_COSTS,
                false,
            );
            for hit in hits {
                match &best {
                    None => best = Some((hit, entry)),
                    Some((current, _)) if hit.k < current.k => best = Some((hit, entry)),
                    _ => {}
                }
            }
        }

        best.map(|(m, entry)| AdapterMatch {
            label: entry.label,
            clip_range: m.start..m.end,
            edit_distance: m.k,
            on_reverse_strand: entry.reverse_complement,
        })
    }
}

/// Extracts the 5' and 3' soft-clip byte slices from a record.
///
/// Hard-clip (`H`) operations are ignored because they carry no sequence.
/// Returns `(five_prime, three_prime)` where each is `Some(&[u8])` only when
/// a soft-clip operation of non-zero length exists at that end.
pub fn extract_clips(record: &RecordBuf) -> (Option<&[u8]>, Option<&[u8]>) {
    let ops: &[Op] = record.cigar().as_ref();
    if ops.is_empty() {
        return (None, None);
    }

    let seq: &[u8] = record.sequence().as_ref();
    if seq.is_empty() {
        return (None, None);
    }

    let five_prime = leading_softclip(ops).and_then(|len| {
        if len == 0 || len > seq.len() {
            None
        } else {
            Some(&seq[..len])
        }
    });

    let three_prime = trailing_softclip(ops).and_then(|len| {
        if len == 0 || len > seq.len() {
            None
        } else {
            Some(&seq[seq.len() - len..])
        }
    });

    (five_prime, three_prime)
}

/// Returns the length of the 5' soft-clip, skipping leading hard-clips.
fn leading_softclip(ops: &[Op]) -> Option<usize> {
    for op in ops {
        match op.kind() {
            Kind::HardClip => continue,
            Kind::SoftClip => return Some(op.len()),
            _ => return None,
        }
    }
    None
}

/// Returns the length of the 3' soft-clip, skipping trailing hard-clips.
fn trailing_softclip(ops: &[Op]) -> Option<usize> {
    for op in ops.iter().rev() {
        match op.kind() {
            Kind::HardClip => continue,
            Kind::SoftClip => return Some(op.len()),
            _ => return None,
        }
    }
    None
}

/// Reverse-complement of a DNA sequence (A/C/G/T only; other bases pass
/// through unchanged to preserve length).
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            other => other,
        })
        .collect()
}

/// Thread-safe store of unknown clips and their occurrence counts.
///
/// Clips are trimmed to at most `MAX_TRACKED_CLIP_LEN` bytes before being
/// stored, both to bound memory and to avoid treating genomic tails as
/// distinct keys just because of their length.
pub struct UnknownClipStore {
    counts: DashMap<Vec<u8>, u32>,
    /// Clip length threshold above which we ignore the clip. Genuine adapters
    /// are ≤ 50 nt; treating very long clips as adapter candidates is noise.
    max_tracked_len: usize,
}

pub const MAX_TRACKED_CLIP_LEN: usize = 64;

impl Default for UnknownClipStore {
    fn default() -> Self {
        Self {
            counts: DashMap::new(),
            max_tracked_len: MAX_TRACKED_CLIP_LEN,
        }
    }
}

impl UnknownClipStore {
    /// Records an unknown clip occurrence.
    pub fn record(&self, clip: &[u8]) {
        if clip.is_empty() || clip.len() > self.max_tracked_len {
            return;
        }
        *self.counts.entry(clip.to_vec()).or_insert(0) += 1;
    }

    /// Returns `(clip, count)` entries whose count is ≥ `threshold`,
    /// sorted by descending count.
    pub fn novel_candidates(&self, threshold: u32) -> Vec<(Vec<u8>, u32)> {
        let mut out: Vec<(Vec<u8>, u32)> = self
            .counts
            .iter()
            .filter(|e| *e.value() >= threshold)
            .map(|e| (e.key().clone(), *e.value()))
            .collect();
        out.sort_unstable_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
        out
    }

    /// Total number of distinct unknown clips tracked.
    pub fn distinct(&self) -> usize {
        self.counts.len()
    }
}

/// Shared handle type for the detector components.
pub type SharedAdapterDb = Arc<AdapterDb>;
pub type SharedUnknownClips = Arc<UnknownClipStore>;

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam::alignment::record_buf::{Cigar as CigarBuf, Sequence as SequenceBuf};

    fn record_with(cigar: Vec<Op>, seq: &[u8]) -> RecordBuf {
        let mut record = RecordBuf::default();
        *record.cigar_mut() = CigarBuf::from(cigar);
        *record.sequence_mut() = SequenceBuf::from(seq.to_vec());
        record
    }

    #[test]
    fn reverse_complement_basics() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AATT"), b"AATT");
        assert_eq!(reverse_complement(b"AAGCAG"), b"CTGCTT");
    }

    #[test]
    fn extract_clips_handles_no_clips() {
        let record = record_with(vec![Op::new(Kind::Match, 10)], b"ACACACACAC");
        let (five, three) = extract_clips(&record);
        assert!(five.is_none());
        assert!(three.is_none());
    }

    #[test]
    fn extract_clips_both_ends() {
        let record = record_with(
            vec![
                Op::new(Kind::SoftClip, 3),
                Op::new(Kind::Match, 4),
                Op::new(Kind::SoftClip, 2),
            ],
            b"ACGTTTTGG",
        );
        let (five, three) = extract_clips(&record);
        assert_eq!(five, Some(&b"ACG"[..]));
        assert_eq!(three, Some(&b"GG"[..]));
    }

    #[test]
    fn extract_clips_skips_hard_clips() {
        let record = record_with(
            vec![
                Op::new(Kind::HardClip, 5),
                Op::new(Kind::SoftClip, 2),
                Op::new(Kind::Match, 4),
                Op::new(Kind::SoftClip, 2),
                Op::new(Kind::HardClip, 3),
            ],
            b"ACGGGGAA",
        );
        let (five, three) = extract_clips(&record);
        assert_eq!(five, Some(&b"AC"[..]));
        assert_eq!(three, Some(&b"AA"[..]));
    }

    #[test]
    fn extract_clips_zero_length_seq() {
        let record = record_with(vec![Op::new(Kind::SoftClip, 0)], b"");
        let (five, three) = extract_clips(&record);
        assert!(five.is_none());
        assert!(three.is_none());
    }

    #[test]
    fn db_exact_match_finds_primer() {
        let db = AdapterDb::new(10, 3).unwrap();
        // Known Iso-Seq 3' primer stem.
        let clip = b"AAGCAGTGGTATCAACGCAGAGTAC";
        let m = db.match_adapter(clip).expect("exact hit");
        assert_eq!(m.edit_distance, 0);
        assert!(
            m.label.contains("primer_3p")
                || m.label.contains("nebnext")
                || m.label.contains("primer_5p")
        );
    }

    #[test]
    fn db_reverse_complement_hit_is_flagged() {
        let db = AdapterDb::new(10, 3).unwrap();
        let forward = b"AAGCAGTGGTATCAACGCAGAGTAC";
        let rc = reverse_complement(forward);
        let m = db.match_adapter(&rc).expect("rc hit");
        assert_eq!(m.edit_distance, 0);
        assert!(m.on_reverse_strand);
    }

    #[test]
    fn db_fuzzy_tolerates_small_edit() {
        let db = AdapterDb::new(10, 3).unwrap();
        // Introduce a single mismatch into the polyT tail.
        let clip = b"TTTTTATTTTTTTTTTTT";
        let m = db.match_adapter(clip);
        assert!(m.is_some());
        assert!(m.unwrap().edit_distance <= 3);
    }

    #[test]
    fn db_returns_none_for_random_clip() {
        let db = AdapterDb::new(10, 3).unwrap();
        let clip = b"GCGCGCGCGCGCGCGCGC";
        let m = db.match_adapter(clip);
        assert!(m.is_none(), "unexpected match: {m:?}");
    }

    #[test]
    fn unknown_clip_store_counts_and_sorts() {
        let store = UnknownClipStore::default();
        store.record(b"AAAATTTT");
        store.record(b"AAAATTTT");
        store.record(b"CCCCGGGG");
        store.record(b"AAAATTTT");
        let out = store.novel_candidates(2);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].0, b"AAAATTTT".to_vec());
        assert_eq!(out[0].1, 3);
    }

    #[test]
    fn unknown_clip_store_drops_long_clips() {
        let store = UnknownClipStore::default();
        let long = vec![b'A'; MAX_TRACKED_CLIP_LEN + 1];
        store.record(&long);
        assert_eq!(store.distinct(), 0);
    }
}
