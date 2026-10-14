// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Pipeline orchestration: BAM scan → predicate → extract.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use log::info;
use rustc_hash::{FxHashMap, FxHashSet};
use smallvec::SmallVec;

use crate::accept::run_accept;
use crate::bam_scan::{scan, AlnGroups};
use crate::cli::Cli;
use crate::extract::{write_fasta_from_bam, write_fasta_from_reads, write_names, CandidateSet};
use crate::logging::elapsed;
use crate::predicate::candidate_rules;
use crate::types::{MiniAln, OutputFormat, PredicateConfig, RuleSet};

/// Selected read names with the rules that fired for each.
pub type CandidateRules = FxHashMap<Vec<u8>, RuleSet>;

/// Runs the iso-align pipeline end-to-end.
pub fn run(cli: Cli) -> Result<()> {
    if let Some(realigned) = cli.realigned_bam.as_deref() {
        info!(
            "[{}] start accept pass1={} realigned={} output={} min_identity={} max_clip={}",
            elapsed(),
            cli.bam.display(),
            realigned.display(),
            cli.output.display(),
            cli.min_identity,
            cli.max_clip
        );
        run_accept(
            &cli.bam,
            realigned,
            &cli.output,
            cli.report.as_deref(),
            &cli.accept_config(),
            cli.bgzf_workers(),
        )?;
        info!("[{}] done!", elapsed());
        return Ok(());
    }

    info!(
        "[{}] start bam={} output={} reads={:?} format={:?} threads={} min_gap={} max_gap={:?} min_flank_clip={} flank_side={:?} rules={} scan={:?} long_intron_min_skip={} long_intron_min_clip={} orphan_min_clip={}",
        elapsed(),
        cli.bam.display(),
        cli.output.display(),
        cli.reads.as_deref(),
        cli.output_format,
        cli.threads,
        cli.min_gap,
        cli.max_gap,
        cli.min_flank_clip,
        cli.flank_side,
        RuleSet::from_rules(&cli.rules).labels(),
        cli.scan_config(),
        cli.long_intron_min_skip,
        cli.long_intron_min_clip,
        cli.orphan_min_clip
    );

    let groups = scan(&cli.bam, cli.bgzf_workers(), &cli.scan_config())?;

    let predicate_cfg = cli.predicate();
    let selected = select_candidate_rules(&groups, &predicate_cfg);
    let mut per_rule: Vec<String> = Vec::new();
    for rule in crate::types::Rule::ALL {
        if predicate_cfg.rules.contains(rule) {
            let n = selected.values().filter(|r| r.contains(rule)).count();
            per_rule.push(format!("{}={n}", rule.label()));
        }
    }
    info!(
        "[{}] predicate done candidates={} of {} qnames ({})",
        elapsed(),
        selected.len(),
        groups.len(),
        per_rule.join(" ")
    );

    if let Some(report) = cli.report.as_deref() {
        write_report(report, &groups, &selected, &predicate_cfg)?;
    }
    let candidates: CandidateSet = selected.keys().cloned().collect();

    match cli.output_format {
        OutputFormat::Names => {
            write_names(&cli.output, &candidates)?;
        }
        OutputFormat::Fasta => match cli.reads.as_deref() {
            Some(reads) if !reads.is_empty() => {
                write_fasta_from_reads(reads, &cli.output, &candidates)?;
            }
            _ => {
                write_fasta_from_bam(&cli.bam, &cli.output, &candidates, cli.bgzf_workers())?;
            }
        },
    }

    info!("[{}] done!", elapsed());
    Ok(())
}

/// Applies the predicate to every QNAME group and returns the set of
/// candidate names. Pulled out of `run` so it's straightforward to test.
pub fn select_candidates(groups: &AlnGroups, cfg: &PredicateConfig) -> CandidateSet {
    let mut out: CandidateSet = FxHashSet::default();
    out.extend(select_candidate_rules(groups, cfg).into_keys());
    out
}

/// Like [`select_candidates`], keeping the rules that fired for each read.
pub fn select_candidate_rules(groups: &AlnGroups, cfg: &PredicateConfig) -> CandidateRules {
    let mut out: CandidateRules = FxHashMap::default();
    for (name, alns) in groups {
        let fired = candidate_rules(alns, cfg);
        if !fired.is_empty() {
            out.insert(name.clone(), fired);
        }
    }
    out
}

/// Writes a per-candidate diagnostic TSV that's useful for spot-checking the
/// flagged reads in IGV.
fn write_report(
    path: &Path,
    groups: &AlnGroups,
    candidates: &CandidateRules,
    _cfg: &PredicateConfig,
) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "name\trname\tmax_gap\tleft_clip_outer\tright_clip_outer\tn_alns\trules"
    )?;

    let mut rows: Vec<&Vec<u8>> = candidates.keys().collect();
    rows.sort_unstable();

    for name in rows {
        let alns = match groups.get(name) {
            Some(a) => a,
            None => continue,
        };
        let mut sorted: SmallVec<[&MiniAln; 4]> = alns.iter().collect();
        sorted.sort_unstable_by_key(|a| a.ref_start);

        let max_gap = sorted
            .windows(2)
            .map(|w| w[1].ref_start - w[0].ref_end)
            .max()
            .unwrap_or(0);
        let left_outer = sorted.first().map(|a| a.left_clip).unwrap_or(0);
        let right_outer = sorted.last().map(|a| a.right_clip).unwrap_or(0);
        let rname = sorted.first().map(|a| a.rname).unwrap_or(0);

        writer.write_all(name)?;
        writeln!(
            writer,
            "\t{rname}\t{max_gap}\t{left_outer}\t{right_outer}\t{}\t{}",
            alns.len(),
            candidates[name].labels()
        )?;
    }
    writer.flush()?;
    info!("[{}] wrote report {}", elapsed(), path.display());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{FlankSide, MiniAln};
    use smallvec::smallvec;

    fn aln(rname: u32, start: i64, end: i64, lc: u32, rc: u32) -> MiniAln {
        MiniAln {
            rname,
            ref_start: start,
            ref_end: end,
            is_reverse: false,
            left_clip: lc,
            right_clip: rc,
            is_supplementary: false,
            has_sa: false,
            max_skip: 0,
            ins_near_skip: false,
        }
    }

    #[test]
    fn select_candidates_keeps_only_passing() {
        let mut groups: AlnGroups = Default::default();
        groups.insert(
            b"good".to_vec(),
            smallvec![aln(0, 100, 200, 0, 50), aln(0, 600_000, 600_100, 50, 0)],
        );
        groups.insert(b"single".to_vec(), smallvec![aln(0, 100, 200, 0, 0)]);
        groups.insert(
            b"small_gap".to_vec(),
            smallvec![aln(0, 100, 200, 0, 50), aln(0, 1_000, 1_100, 50, 0)],
        );

        let cfg = PredicateConfig {
            min_gap: 300_000,
            max_gap: None,
            min_flank_clip: 20,
            flank_side: FlankSide::Inner,
            rules: RuleSet::from_rules(&[crate::types::Rule::Split]),
            long_intron_min_skip: 150_000,
            long_intron_min_clip: 21,
            orphan_min_clip: 50,
        };
        let cands = select_candidates(&groups, &cfg);
        assert_eq!(cands.len(), 1);
        assert!(cands.contains(b"good".as_slice()));
    }
}
