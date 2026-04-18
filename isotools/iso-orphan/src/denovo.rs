// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! De novo clustering module for grouping reads without reference annotation.
//!
//! Provides intron-chain-based clustering for multi-exon reads and
//! reciprocal-overlap-based clustering for single-exon reads.
//!
//! Clustering is performed within packbed components — reads are already
//! pre-grouped by genomic overlap. This module sub-clusters by structural
//! similarity (intron chain for multi-exon, interval overlap for single-exon).

use std::collections::HashMap;

use genepred::GenePred;

use crate::scoring::{intron_chain, intron_chains_match, reciprocal_overlap};

/// Cluster multi-exon reads by intron chain similarity.
///
/// Reads with matching intron chains (within `tolerance` bp at each junction)
/// are grouped together. Uses greedy first-match assignment against cluster
/// representatives (first member of each cluster).
///
/// Returns a map from cluster ID to member indices (into the input slice).
pub fn cluster_multi_exon(reads: &[&GenePred], tolerance: u64) -> HashMap<usize, Vec<usize>> {
    if reads.is_empty() {
        return HashMap::new();
    }

    let chains: Vec<Vec<(u64, u64)>> = reads.iter().map(|r| intron_chain(r)).collect();

    let mut clusters: Vec<(Vec<(u64, u64)>, Vec<usize>)> = Vec::new();

    log::debug!("Clustering {} reads", chains.len());

    for (idx, chain) in chains.iter().enumerate() {
        let mut found = false;
        for (rep_chain, members) in clusters.iter_mut() {
            if intron_chains_match(chain, rep_chain, tolerance) {
                members.push(idx);
                found = true;
                break;
            }
        }
        if !found {
            clusters.push((chain.clone(), vec![idx]));
        }
    }

    log::debug!(
        "DEBUG: Clustered {} reads into {} clusters",
        chains.len(),
        clusters.len()
    );

    clusters
        .into_iter()
        .enumerate()
        .map(|(id, (_, members))| (id, members))
        .collect()
}

/// Cluster single-exon reads by reciprocal overlap of their exonic span.
///
/// Uses greedy clustering against each cluster's representative (first member).
/// This avoids the "snowball" effect that would occur if the cluster interval
/// grew with each new member.
///
/// Returns a map from cluster ID to member indices (into the input slice).
pub fn cluster_single_exon(reads: &[&GenePred], min_overlap: f64) -> HashMap<usize, Vec<usize>> {
    if reads.is_empty() {
        return HashMap::new();
    }

    let spans: Vec<(u64, u64)> = reads
        .iter()
        .map(|r| {
            let exons = r.exons();
            if exons.is_empty() {
                (r.start(), r.end())
            } else {
                (exons[0].0, exons[0].1)
            }
        })
        .collect();

    // Sort by start for better locality
    let mut order: Vec<usize> = (0..reads.len()).collect();
    order.sort_by_key(|&i| spans[i].0);

    let mut clusters: Vec<((u64, u64), Vec<usize>)> = Vec::new();

    for &idx in &order {
        let (rs, re) = spans[idx];
        let mut best_cluster = None;
        let mut best_ov = 0.0f64;

        for (cid, (rep_span, _)) in clusters.iter().enumerate() {
            let ov = reciprocal_overlap(rs, re, rep_span.0, rep_span.1);
            if ov >= min_overlap && ov > best_ov {
                best_ov = ov;
                best_cluster = Some(cid);
            }
        }

        if let Some(cid) = best_cluster {
            clusters[cid].1.push(idx);
        } else {
            clusters.push(((rs, re), vec![idx]));
        }
    }

    clusters
        .into_iter()
        .enumerate()
        .map(|(id, (_, members))| (id, members))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap as StdHashMap;

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
            extras: StdHashMap::new(),
        }
    }

    #[test]
    fn test_cluster_multi_exon_same_chain() {
        // 3 reads with identical intron chain (200,300), (400,500) → 1 cluster
        let r1 = make_record(b"r1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let r2 = make_record(b"r2", 110, 590, vec![110, 300, 500], vec![200, 400, 590]);
        let r3 = make_record(b"r3", 105, 595, vec![105, 300, 500], vec![200, 400, 595]);

        let reads: Vec<&GenePred> = vec![&r1, &r2, &r3];
        let clusters = cluster_multi_exon(&reads, 0);

        assert_eq!(clusters.len(), 1);
        assert_eq!(clusters.values().next().unwrap().len(), 3);
    }

    #[test]
    fn test_cluster_multi_exon_different_chains() {
        let r1 = make_record(b"r1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        let r2 = make_record(b"r2", 100, 600, vec![100, 350, 500], vec![250, 450, 600]);

        let reads: Vec<&GenePred> = vec![&r1, &r2];
        let clusters = cluster_multi_exon(&reads, 0);

        assert_eq!(clusters.len(), 2);
    }

    #[test]
    fn test_cluster_multi_exon_tolerance_merges() {
        let r1 = make_record(b"r1", 100, 600, vec![100, 300, 500], vec![200, 400, 600]);
        // Introns differ by 2bp → within tolerance 3
        let r2 = make_record(b"r2", 100, 600, vec![100, 298, 500], vec![202, 402, 600]);

        let reads: Vec<&GenePred> = vec![&r1, &r2];

        assert_eq!(cluster_multi_exon(&reads, 0).len(), 2);
        assert_eq!(cluster_multi_exon(&reads, 3).len(), 1);
    }

    #[test]
    fn test_cluster_single_exon_overlapping() {
        let r1 = make_record(b"r1", 100, 300, vec![100], vec![300]);
        let r2 = make_record(b"r2", 120, 320, vec![120], vec![320]);
        let r3 = make_record(b"r3", 110, 310, vec![110], vec![310]);

        let reads: Vec<&GenePred> = vec![&r1, &r2, &r3];
        let clusters = cluster_single_exon(&reads, 0.5);

        assert_eq!(clusters.len(), 1);
    }

    #[test]
    fn test_cluster_single_exon_separate() {
        let r1 = make_record(b"r1", 100, 200, vec![100], vec![200]);
        let r2 = make_record(b"r2", 500, 600, vec![500], vec![600]);

        let reads: Vec<&GenePred> = vec![&r1, &r2];
        let clusters = cluster_single_exon(&reads, 0.5);

        assert_eq!(clusters.len(), 2);
    }

    #[test]
    fn test_cluster_preserves_all_reads() {
        let reads_data: Vec<GenePred> = (0..10)
            .map(|i| {
                let offset = i * 5;
                make_record(
                    format!("r{}", i).as_bytes(),
                    100 + offset,
                    300 + offset,
                    vec![100 + offset],
                    vec![300 + offset],
                )
            })
            .collect();

        let reads: Vec<&GenePred> = reads_data.iter().collect();
        let clusters = cluster_single_exon(&reads, 0.5);

        let total: usize = clusters.values().map(|v| v.len()).sum();
        assert_eq!(total, 10, "All reads must be assigned to a cluster");

        // Verify disjointness: no read in multiple clusters
        let mut seen = std::collections::HashSet::new();
        for members in clusters.values() {
            for &idx in members {
                assert!(seen.insert(idx), "Read {} in multiple clusters", idx);
            }
        }
    }

    #[test]
    fn test_cluster_multi_exon_preserves_all_reads() {
        let reads_data: Vec<GenePred> = (0..6)
            .map(|i| {
                // Alternate between two intron chains
                if i % 2 == 0 {
                    make_record(
                        format!("r{}", i).as_bytes(),
                        100,
                        600,
                        vec![100, 300, 500],
                        vec![200, 400, 600],
                    )
                } else {
                    make_record(
                        format!("r{}", i).as_bytes(),
                        100,
                        600,
                        vec![100, 350, 500],
                        vec![250, 450, 600],
                    )
                }
            })
            .collect();

        let reads: Vec<&GenePred> = reads_data.iter().collect();
        let clusters = cluster_multi_exon(&reads, 0);

        let total: usize = clusters.values().map(|v| v.len()).sum();
        assert_eq!(total, 6, "All reads must be assigned to a cluster");
        assert_eq!(
            clusters.len(),
            2,
            "Two distinct intron chains → two clusters"
        );
    }
}
