// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Scoring module for evaluating transcript-level evidence.
//!
//! Provides reusable functions for:
//! - Junction-based scoring of multi-exon reads against reference transcripts
//! - Reciprocal overlap scoring for single-exon reads
//! - Exon boundary matching as a fallback
//! - Intron chain extraction and comparison for clustering
//! - Reference-bounds checking for guided mode partitioning
//! - Per-intron support counting for de novo evaluation

use genepred::GenePred;
use std::collections::{HashMap, HashSet};

/// Parameters controlling scoring and classification thresholds.
#[derive(Debug, Clone)]
pub struct ScoringParams {
    /// Tolerance in bp for matching splice junctions (default: 0 = exact)
    pub junction_tolerance: u64,
    /// Minimum fraction of query junctions that must match a reference (default: 0.5)
    pub min_junction_frac: f64,
    /// Minimum reciprocal overlap for single-exon reads (default: 0.5)
    pub min_overlap_frac: f64,
    /// Minimum reads per de novo cluster or component to keep (default: 5)
    pub min_cluster_support: usize,
    /// Tolerance for terminal exon boundary matching in bp (default: 50)
    pub end_tolerance: u64,
    /// Minimum fraction of a read's introns that must be individually supported (default: 0.5)
    pub min_intron_support_frac: f64,
    /// Minimum fraction of component reads containing an intron for it to count as "supported" (default: 0.5)
    pub intron_support_threshold: f64,
    /// Minimum median splice-site score from bigWig to keep a read (default: 0.5)
    pub min_splice_score: f32,
}

impl Default for ScoringParams {
    fn default() -> Self {
        Self {
            junction_tolerance: 0,
            min_junction_frac: 0.5,
            min_overlap_frac: 0.5,
            min_cluster_support: 5,
            end_tolerance: 50,
            min_intron_support_frac: 0.5,
            intron_support_threshold: 0.5,
            min_splice_score: 0.5,
        }
    }
}

/// Result of scoring a multi-exon query against a reference transcript by junction comparison.
#[derive(Debug, Clone)]
pub struct JunctionScore {
    pub junction_matches: usize,
    pub query_junction_count: usize,
    pub ref_junction_count: usize,
    pub query_junction_frac: f64,
    pub ref_junction_frac: f64,
}

impl JunctionScore {
    /// Whether this score meets the minimum junction fraction threshold.
    pub fn passes(&self, min_frac: f64) -> bool {
        self.query_junction_frac >= min_frac
    }
}

/// Score a multi-exon query against a single reference transcript by intron comparison.
///
/// For each query intron (donor-acceptor pair), checks if a matching reference intron
/// exists within `tolerance` bp at both sites. Each reference intron is matched at most once.
///
/// # Examples
///
/// ```text
/// Query:     [100-200]---[300-400]---[500-600]  → introns: (200,300), (400,500)
/// Reference: [100-200]---[300-400]---[500-600]  → introns: (200,300), (400,500)
/// → junction_matches=2, query_junction_frac=1.0, ref_junction_frac=1.0
/// ```
pub fn score_junctions(query: &GenePred, reference: &GenePred, tolerance: u64) -> JunctionScore {
    let q_introns = query.introns();
    let r_introns = reference.introns();

    if q_introns.is_empty() || r_introns.is_empty() {
        return JunctionScore {
            junction_matches: 0,
            query_junction_count: q_introns.len(),
            ref_junction_count: r_introns.len(),
            query_junction_frac: 0.0,
            ref_junction_frac: 0.0,
        };
    }

    let mut used = vec![false; r_introns.len()];
    let mut matches = 0usize;

    for (qs, qe) in &q_introns {
        for (j, (rs, re)) in r_introns.iter().enumerate() {
            if !used[j] && qs.abs_diff(*rs) <= tolerance && qe.abs_diff(*re) <= tolerance {
                matches += 1;
                used[j] = true;
                break;
            }
        }
    }

    let q_frac = matches as f64 / q_introns.len() as f64;
    let r_frac = matches as f64 / r_introns.len() as f64;

    JunctionScore {
        junction_matches: matches,
        query_junction_count: q_introns.len(),
        ref_junction_count: r_introns.len(),
        query_junction_frac: q_frac,
        ref_junction_frac: r_frac,
    }
}

/// Score a multi-exon query against multiple reference transcripts, returning the best score.
///
/// Only considers references with introns (multi-exon). Returns None if no multi-exon
/// references exist.
pub fn best_junction_score(
    query: &GenePred,
    references: &[GenePred],
    tolerance: u64,
) -> Option<JunctionScore> {
    references
        .iter()
        .filter(|r| !r.introns().is_empty())
        .map(|r| score_junctions(query, r, tolerance))
        .max_by(|a, b| {
            a.query_junction_frac
                .partial_cmp(&b.query_junction_frac)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
}

/// Check if any exon boundary of the query matches any reference exon boundary.
///
/// Used as a fallback for multi-exon queries that fail junction scoring but may still
/// share splice sites with the reference, indicating partial structural support.
pub fn has_boundary_match(query: &GenePred, references: &[GenePred], tolerance: u64) -> bool {
    let ref_boundaries: HashSet<u64> = references
        .iter()
        .flat_map(|r| r.exons())
        .flat_map(|(s, e)| [s, e])
        .collect();

    if tolerance == 0 {
        query
            .exons()
            .iter()
            .any(|(s, e)| ref_boundaries.contains(s) || ref_boundaries.contains(e))
    } else {
        query.exons().iter().any(|(s, e)| {
            ref_boundaries
                .iter()
                .any(|&rb| s.abs_diff(rb) <= tolerance || e.abs_diff(rb) <= tolerance)
        })
    }
}

/// Compute reciprocal overlap between two intervals [a_start, a_end) and [b_start, b_end).
///
/// Returns min(overlap / a_length, overlap / b_length).
/// Returns 0.0 if intervals don't overlap or have zero length.
pub fn reciprocal_overlap(a_start: u64, a_end: u64, b_start: u64, b_end: u64) -> f64 {
    let overlap_start = a_start.max(b_start);
    let overlap_end = a_end.min(b_end);

    if overlap_start >= overlap_end {
        return 0.0;
    }

    let overlap = (overlap_end - overlap_start) as f64;
    let a_len = (a_end - a_start) as f64;
    let b_len = (b_end - b_start) as f64;

    if a_len == 0.0 || b_len == 0.0 {
        return 0.0;
    }

    (overlap / a_len).min(overlap / b_len)
}

/// Score a single-exon query against all reference exons by reciprocal overlap.
///
/// Returns the best reciprocal overlap of the query's single exon with any reference exon.
pub fn best_single_exon_overlap(query: &GenePred, references: &[GenePred]) -> f64 {
    let q_exons = query.exons();
    if q_exons.is_empty() {
        return 0.0;
    }
    let (qs, qe) = q_exons[0];

    references
        .iter()
        .flat_map(|r| r.exons())
        .map(|(rs, re)| reciprocal_overlap(qs, qe, rs, re))
        .fold(0.0f64, f64::max)
}

/// Extract the intron chain from a record (list of intron intervals).
pub fn intron_chain(record: &GenePred) -> Vec<(u64, u64)> {
    record.introns()
}

/// Check if two intron chains match within a given tolerance at each junction.
pub fn intron_chains_match(a: &[(u64, u64)], b: &[(u64, u64)], tolerance: u64) -> bool {
    if a.len() != b.len() {
        return false;
    }
    a.iter().zip(b.iter()).all(|((as_, ae), (bs, be))| {
        as_.abs_diff(*bs) <= tolerance && ae.abs_diff(*be) <= tolerance
    })
}

// ---------------------------------------------------------------------------
// Reference-overlap helpers
// ---------------------------------------------------------------------------

/// Check if any exon of the query overlaps any exon of any reference transcript.
///
/// This is the guided-vs-de-novo partition test: a read whose exons touch at least
/// one reference exon is scored in guided mode; a read with no exon overlap at all
/// is redirected to de novo evaluation.
pub fn overlaps_any_reference_exon(read: &GenePred, references: &[GenePred]) -> bool {
    let read_exons = read.exons();
    if read_exons.is_empty() {
        return false;
    }
    for (rs, re) in &read_exons {
        for reference in references {
            for (es, ee) in reference.exons() {
                // Two intervals overlap iff start_a < end_b AND start_b < end_a
                if *rs < ee && es < *re {
                    return true;
                }
            }
        }
    }
    false
}

// ---------------------------------------------------------------------------
// Per-intron support helpers
// ---------------------------------------------------------------------------

/// Count how many reads contain each intron (exact coordinates).
///
/// Returns a map from intron (start, end) to the number of reads containing it.
pub fn compute_intron_support(reads: &[&GenePred]) -> HashMap<(u64, u64), usize> {
    let mut support: HashMap<(u64, u64), usize> = HashMap::new();
    for read in reads {
        for intron in read.introns() {
            *support.entry(intron).or_insert(0) += 1;
        }
    }
    support
}

/// Compute the fraction of a read's introns that are "supported".
///
/// An intron is supported if it appears in >= `threshold` fraction of the `total_reads`.
/// Returns 0.0 for reads with no introns.
pub fn intron_support_fraction(
    read: &GenePred,
    support: &HashMap<(u64, u64), usize>,
    total_reads: usize,
    threshold: f64,
) -> f64 {
    let introns = read.introns();
    if introns.is_empty() || total_reads == 0 {
        return 0.0;
    }

    let supported = introns
        .iter()
        .filter(|intron| {
            let count = support.get(intron).copied().unwrap_or(0);
            (count as f64 / total_reads as f64) >= threshold
        })
        .count();

    supported as f64 / introns.len() as f64
}

/// Compute the median of a mutable f32 slice. Returns None for empty slices.
pub fn median_f32(values: &mut [f32]) -> Option<f32> {
    if values.is_empty() {
        return None;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mid = values.len() / 2;
    if values.len() % 2 == 0 {
        Some((values[mid - 1] + values[mid]) / 2.0)
    } else {
        Some(values[mid])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_record(
        name: &[u8],
        start: u64,
        end: u64,
        block_starts: Vec<u64>,
        block_ends: Vec<u64>,
    ) -> GenePred {
        GenePred {
            chrom: b"chr1".to_vec(),
            start,
            end,
            name: Some(name.to_vec()),
            strand: None,
            thick_start: Some(start),
            thick_end: Some(end),
            block_count: Some(block_starts.len() as u32),
            block_starts: Some(block_starts),
            block_ends: Some(block_ends),
            extras: HashMap::new(),
        }
    }

    // -- reciprocal overlap tests --

    #[test]
    fn test_reciprocal_overlap_full() {
        assert_eq!(reciprocal_overlap(100, 200, 100, 200), 1.0);
    }

    #[test]
    fn test_reciprocal_overlap_half() {
        assert_eq!(reciprocal_overlap(100, 200, 150, 250), 0.5);
    }

    #[test]
    fn test_reciprocal_overlap_none() {
        assert_eq!(reciprocal_overlap(100, 200, 300, 400), 0.0);
    }

    #[test]
    fn test_reciprocal_overlap_contained() {
        let ov = reciprocal_overlap(100, 400, 150, 200);
        assert!((ov - 50.0 / 300.0).abs() < 1e-10);
    }

    // -- junction scoring tests --

    #[test]
    fn test_score_junctions_exact_match() {
        let query = make_record(b"q1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let reference = make_record(b"r1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);

        let score = score_junctions(&query, &reference, 0);
        assert_eq!(score.junction_matches, 2);
        assert_eq!(score.query_junction_frac, 1.0);
        assert_eq!(score.ref_junction_frac, 1.0);
        assert!(score.passes(0.5));
    }

    #[test]
    fn test_score_junctions_partial_match() {
        let query = make_record(b"q1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let reference = make_record(
            b"r1",
            100,
            800,
            vec![100, 300, 500, 700],
            vec![200, 400, 600, 800],
        );

        let score = score_junctions(&query, &reference, 0);
        assert_eq!(score.junction_matches, 2);
        assert_eq!(score.query_junction_frac, 1.0);
        assert!((score.ref_junction_frac - 2.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_score_junctions_no_match() {
        let query = make_record(b"q1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let reference = make_record(b"r1", 100, 600, vec![100, 350, 550], vec![250, 450, 600]);

        let score = score_junctions(&query, &reference, 0);
        assert_eq!(score.junction_matches, 0);
        assert!(!score.passes(0.5));
    }

    #[test]
    fn test_score_junctions_with_tolerance() {
        let query = make_record(b"q1", 100, 600, vec![100, 302, 500], vec![200, 398, 600]);
        let reference = make_record(b"r1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);

        let score = score_junctions(&query, &reference, 3);
        assert_eq!(score.junction_matches, 2);
        assert!(score.passes(0.5));
    }

    #[test]
    fn test_best_junction_score_prefers_best_match() {
        let query = make_record(b"q1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let ref1 = make_record(b"r1", 100, 600, vec![100, 350, 550], vec![250, 450, 600]);
        let ref2 = make_record(b"r2", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);

        let best = best_junction_score(&query, &[ref1, ref2], 0).unwrap();
        assert_eq!(best.junction_matches, 2);
        assert_eq!(best.query_junction_frac, 1.0);
    }

    #[test]
    fn test_best_junction_score_no_multi_exon_refs() {
        let query = make_record(b"q1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let se_ref = make_record(b"r1", 100, 600, vec![100], vec![600]);

        assert!(best_junction_score(&query, &[se_ref], 0).is_none());
    }

    #[test]
    fn test_best_single_exon_overlap_exact() {
        let query = make_record(b"q1", 100, 200, vec![100], vec![200]);
        let ref1 = make_record(b"r1", 100, 200, vec![100], vec![200]);
        let ref2 = make_record(b"r2", 300, 400, vec![300], vec![400]);

        assert_eq!(best_single_exon_overlap(&query, &[ref1, ref2]), 1.0);
    }

    #[test]
    fn test_has_boundary_match_positive() {
        let query = make_record(b"q1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let reference = make_record(b"r1", 100, 400, vec![100, 300], vec![200, 400]);
        assert!(has_boundary_match(&query, &[reference], 0));
    }

    #[test]
    fn test_has_boundary_match_negative() {
        let query = make_record(b"q1", 100, 300, vec![100, 250], vec![200, 300]);
        let reference = make_record(b"r1", 500, 800, vec![500, 700], vec![600, 800]);

        assert!(!has_boundary_match(&query, &[reference], 0));
    }

    #[test]
    fn test_intron_chains_match_exact() {
        assert!(intron_chains_match(
            &[(200, 300), (400, 500)],
            &[(200, 300), (400, 500)],
            0
        ));
    }

    #[test]
    fn test_intron_chains_match_different_length() {
        assert!(!intron_chains_match(
            &[(200, 300)],
            &[(200, 300), (400, 500)],
            0
        ));
    }

    #[test]
    fn test_intron_chains_match_with_tolerance() {
        let a = vec![(200, 302), (398, 500)];
        let b = vec![(200, 300), (400, 500)];
        assert!(intron_chains_match(&a, &b, 3));
        assert!(!intron_chains_match(&a, &b, 1));
    }

    #[test]
    fn test_junction_score_beats_weak_overlap() {
        let reference = make_record(b"r1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let good_query = make_record(b"q1", 110, 590, vec![110, 300, 500], vec![200, 400, 590]);
        let bad_query = make_record(b"q2", 100, 600, vec![100, 350, 550], vec![250, 450, 600]);

        let good_score = score_junctions(&good_query, &reference, 0);
        let bad_score = score_junctions(&bad_query, &reference, 0);

        assert!(good_score.query_junction_frac > bad_score.query_junction_frac);
        assert!(good_score.passes(0.5));
        assert!(!bad_score.passes(0.5));
    }

    #[test]
    fn test_multi_isoform_reference() {
        let ref1 = make_record(b"r1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let ref2 = make_record(b"r2", 100, 800, vec![100, 350, 700], vec![300, 600, 800]);
        let query = make_record(b"q1", 100, 800, vec![100, 350, 700], vec![300, 600, 800]);

        let best = best_junction_score(&query, &[ref1, ref2], 0).unwrap();
        assert_eq!(best.junction_matches, 2);
        assert_eq!(best.query_junction_frac, 1.0);
        assert!(best.passes(0.5));
    }

    #[test]
    fn test_single_exon_vs_multi_exon_ref() {
        let query = make_record(b"q1", 300, 400, vec![300], vec![400]);
        let reference = make_record(b"r1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);

        let overlap = best_single_exon_overlap(&query, &[reference]);
        assert_eq!(overlap, 1.0);
    }

    // -- reference overlap tests --

    #[test]
    fn test_overlaps_ref_exon_exact() {
        let read = make_record(b"r", 300, 400, vec![300], vec![400]);
        let refs = vec![make_record(b"ref", 100, 600, vec![100, 300, 500], vec![200, 400, 600])];
        // Read exon [300,400) overlaps ref exon [300,400)
        assert!(overlaps_any_reference_exon(&read, &refs));
    }

    #[test]
    fn test_overlaps_ref_exon_partial() {
        let read = make_record(b"r", 150, 250, vec![150], vec![250]);
        let refs = vec![make_record(b"ref", 100, 600, vec![100, 300, 500], vec![200, 400, 600])];
        // Read exon [150,250) overlaps ref exon [100,200)
        assert!(overlaps_any_reference_exon(&read, &refs));
    }

    #[test]
    fn test_no_overlap_intronic() {
        // Read sits entirely within a reference intron
        let read = make_record(b"r", 220, 280, vec![220], vec![280]);
        let refs = vec![make_record(b"ref", 100, 600, vec![100, 300, 500], vec![200, 400, 600])];
        // Ref exons: [100,200), [300,400), [500,600). Read [220,280) is in intron [200,300).
        assert!(!overlaps_any_reference_exon(&read, &refs));
    }

    #[test]
    fn test_no_overlap_intergenic() {
        let read = make_record(b"r", 700, 800, vec![700], vec![800]);
        let refs = vec![make_record(b"ref", 100, 600, vec![100, 300, 500], vec![200, 400, 600])];
        assert!(!overlaps_any_reference_exon(&read, &refs));
    }

    #[test]
    fn test_overlaps_multi_exon_read_partial() {
        // Multi-exon read: one exon overlaps reference, the other does not
        let read = make_record(b"r", 150, 700, vec![150, 650], vec![250, 700]);
        let refs = vec![make_record(b"ref", 100, 600, vec![100, 300, 500], vec![200, 400, 600])];
        // Read exon [150,250) overlaps ref exon [100,200) → true
        assert!(overlaps_any_reference_exon(&read, &refs));
    }

    // -- intron support tests --

    #[test]
    fn test_compute_intron_support() {
        // 3 reads sharing intron (200,300); 2 sharing (400,500); 1 unique (600,700)
        let r1 = make_record(
            b"r1",
            100,
            800,
            vec![100, 300, 500, 700],
            vec![200, 400, 600, 800],
        );
        let r2 = make_record(b"r2", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let r3 = make_record(b"r3", 100, 400, vec![100, 300], vec![200, 400]);

        let reads: Vec<&GenePred> = vec![&r1, &r2, &r3];
        let support = compute_intron_support(&reads);

        assert_eq!(support[&(200, 300)], 3);
        assert_eq!(support[&(400, 500)], 2);
        assert_eq!(support[&(600, 700)], 1);
    }

    #[test]
    fn test_intron_support_fraction_all_supported() {
        let r1 = make_record(b"r1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        // Support: both introns appear in 5/5 reads
        let mut support = HashMap::new();
        support.insert((200, 300), 5);
        support.insert((400, 500), 5);

        let frac = intron_support_fraction(&r1, &support, 5, 0.5);
        assert_eq!(frac, 1.0);
    }

    #[test]
    fn test_intron_support_fraction_nine_of_ten() {
        // Read with 10 introns, 9 are well-supported, 1 is not
        let mut starts = vec![100u64];
        let mut ends = vec![];
        for i in 0..11 {
            let s = 100 + i * 200;
            let e = s + 100;
            if i > 0 {
                starts.push(s);
            }
            ends.push(e);
        }
        let read = make_record(b"r1", starts[0], *ends.last().unwrap(), starts, ends);
        let introns = read.introns();
        assert_eq!(introns.len(), 10);

        let mut support = HashMap::new();
        for (i, intron) in introns.iter().enumerate() {
            if i < 9 {
                support.insert(*intron, 10); // 10/10 reads
            } else {
                support.insert(*intron, 1); // 1/10 reads — not supported at threshold 0.5
            }
        }

        let frac = intron_support_fraction(&read, &support, 10, 0.5);
        assert!((frac - 0.9).abs() < 1e-10);
        assert!(frac >= 0.5); // passes min_intron_support_frac
    }

    #[test]
    fn test_intron_support_fraction_weak() {
        // Read with 4 introns, only 1 is supported
        let read = make_record(
            b"r1",
            100,
            1000,
            vec![100, 300, 500, 700, 900],
            vec![200, 400, 600, 800, 1000],
        );
        let introns = read.introns();
        assert_eq!(introns.len(), 4);

        let mut support = HashMap::new();
        support.insert(introns[0], 8); // supported
        support.insert(introns[1], 1); // not supported
        support.insert(introns[2], 1); // not supported
        support.insert(introns[3], 1); // not supported

        let frac = intron_support_fraction(&read, &support, 10, 0.5);
        assert!((frac - 0.25).abs() < 1e-10);
        assert!(frac < 0.5); // fails min_intron_support_frac
    }

    // -- median tests --

    #[test]
    fn test_median_f32_odd() {
        assert_eq!(median_f32(&mut [3.0, 1.0, 2.0]), Some(2.0));
    }

    #[test]
    fn test_median_f32_even() {
        assert_eq!(median_f32(&mut [4.0, 1.0, 3.0, 2.0]), Some(2.5));
    }

    #[test]
    fn test_median_f32_single() {
        assert_eq!(median_f32(&mut [5.0]), Some(5.0));
    }

    #[test]
    fn test_median_f32_empty() {
        assert_eq!(median_f32(&mut []), None);
    }
}
