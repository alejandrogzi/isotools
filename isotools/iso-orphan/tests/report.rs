// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! End-to-end tests for `<outdir>/orphans/<prefix>.report.tsv`.
//!
//! Each test runs the real binary on a synthetic BED12 fixture and checks the report
//! against the BED files produced by the same run, so a divergence between the two
//! outputs fails the suite.

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use tempfile::TempDir;

/// Column order the report must always use.
const COLUMNS: [&str; 26] = [
    "read_id",
    "chrom",
    "start",
    "end",
    "strand",
    "exon_count",
    "read_type",
    "mode",
    "evaluation_path",
    "group_size",
    "reference_count",
    "cluster_size",
    "overlaps_reference_exon",
    "best_reference_id",
    "best_reference_overlap",
    "junction_matches",
    "junction_count",
    "junction_fraction",
    "boundary_match",
    "boundary_evidence",
    "intron_support_fraction",
    "median_splice_score",
    "weakest_splice_score",
    "applied_thresholds",
    "reasons",
    "decision",
];

/// Build a BED12 line spanning `exons`.
fn bed12(name: &str, exons: &[(u64, u64)]) -> String {
    let start = exons[0].0;
    let end = exons[exons.len() - 1].1;

    let sizes: Vec<String> = exons.iter().map(|(s, e)| (e - s).to_string()).collect();
    let starts: Vec<String> = exons.iter().map(|(s, _)| (s - start).to_string()).collect();

    format!(
        "chr1\t{start}\t{end}\t{name}\t0\t+\t{start}\t{end}\t0,0,0\t{}\t{},\t{},",
        exons.len(),
        sizes.join(","),
        starts.join(",")
    )
}

/// Reference transcript with exons [1000,1200), [1400,1600), [1800,2000).
fn reference_bed() -> String {
    format!(
        "{}\n",
        bed12("refA", &[(1000, 1200), (1400, 1600), (1800, 2000)])
    )
}

/// Query set covering one read per guided and de novo terminal branch.
///
/// The two `q_dup` records are byte-identical, so the BED outputs deduplicate them
/// while the report must still hold one row each.
fn query_bed() -> String {
    let mut records = vec![
        // exact intron chain of the reference
        bed12(
            "q_junction_pass",
            &[(1000, 1200), (1400, 1600), (1800, 2000)],
        ),
        // byte-identical duplicates, both matching the first reference intron
        bed12("q_dup", &[(1000, 1200), (1400, 1600)]),
        bed12("q_dup", &[(1000, 1200), (1400, 1600)]),
        // single exon matching a reference exon exactly
        bed12("q_se_pass", &[(1000, 1200)]),
        // shares the complete junction (1200,1400) with refA but only 1 of its 3
        // junctions, so it needs the coherent rescue
        bed12(
            "q_junction_rescue",
            &[(1100, 1200), (1400, 1500), (1700, 1800), (1900, 1950)],
        ),
        // its only agreement with refA is the single coordinate 1000: not a rescue
        bed12("q_single_coordinate", &[(1000, 1100), (1250, 1350)]),
        // no matching junction and no shared boundary
        bed12("q_low_junction", &[(1050, 1150), (1250, 1350)]),
        // single exon overlapping a reference exon by 0.200 reciprocally
        bed12("q_se_scrap", &[(1150, 1400)]),
    ];

    // five identical chains inside the reference intron [1200,1400): no exon overlap,
    // so they are redirected to de novo and form a dominant cluster
    for i in 0..5 {
        records.push(bed12(&format!("q_oob{}", i), &[(1220, 1250), (1300, 1330)]));
    }

    // a lone read in a component of its own
    records.push(bed12("q_lone", &[(50000, 50200), (50400, 50600)]));

    format!("{}\n", records.join("\n"))
}

/// Number of query records in the fixture.
fn query_record_count() -> usize {
    query_bed().lines().filter(|l| !l.trim().is_empty()).count()
}

struct Run {
    _dir: TempDir,
    outdir: PathBuf,
}

impl Run {
    /// Run the binary on the fixture, with `RAYON_NUM_THREADS` pinned to `threads`.
    fn new(mode: &str, threads: usize) -> Self {
        let dir = TempDir::new().expect("temp dir");
        let query = dir.path().join("query.bed");
        let reference = dir.path().join("ref.bed");
        let outdir = dir.path().join("out");

        fs::write(&query, query_bed()).expect("write query");
        fs::write(&reference, reference_bed()).expect("write reference");

        let output = Command::new(env!("CARGO_BIN_EXE_iso-orphan"))
            .arg("--query")
            .arg(&query)
            // INFO: --ref is required by the CLI argument groups even in de novo mode,
            // where the references are ignored
            .arg("--ref")
            .arg(&reference)
            .arg(mode)
            .arg("--outdir")
            .arg(&outdir)
            .arg("--prefix")
            .arg("test")
            .arg("--threads")
            .arg(threads.to_string())
            .env("RAYON_NUM_THREADS", threads.to_string())
            .output()
            .expect("run iso-orphan");

        assert!(
            output.status.success(),
            "iso-orphan failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );

        Self { _dir: dir, outdir }
    }

    fn path(&self, name: &str) -> PathBuf {
        self.outdir.join("orphans").join(name)
    }

    fn report(&self) -> Report {
        Report::read(&self.path("test.report.tsv"))
    }

    fn bed_lines(&self, name: &str) -> Vec<String> {
        read_lines(&self.path(name))
    }
}

struct Report {
    header: Vec<String>,
    rows: Vec<Vec<String>>,
    raw: String,
}

impl Report {
    fn read(path: &Path) -> Self {
        let raw = fs::read_to_string(path)
            .unwrap_or_else(|e| panic!("cannot read {:?}: {e}", path.display()));
        let mut lines = raw.lines();

        let header: Vec<String> = lines
            .next()
            .expect("report has a header")
            .split('\t')
            .map(str::to_string)
            .collect();

        let rows: Vec<Vec<String>> = lines
            .map(|line| line.split('\t').map(str::to_string).collect())
            .collect();

        Self { header, rows, raw }
    }

    fn index(&self, column: &str) -> usize {
        self.header
            .iter()
            .position(|name| name == column)
            .unwrap_or_else(|| panic!("missing column {column}"))
    }

    fn field(&self, row: &[String], column: &str) -> String {
        row[self.index(column)].clone()
    }

    /// Reason / decision pairs keyed by read id, allowing repeated ids.
    fn outcomes(&self) -> HashMap<String, Vec<(String, String)>> {
        let mut out: HashMap<String, Vec<(String, String)>> = HashMap::new();

        for row in &self.rows {
            out.entry(self.field(row, "read_id"))
                .or_default()
                .push((self.field(row, "reasons"), self.field(row, "decision")));
        }

        out
    }
}

fn read_lines(path: &Path) -> Vec<String> {
    fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("cannot read {:?}: {e}", path.display()))
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(str::to_string)
        .collect()
}

/// `chrom:start-end:name` identity of a BED line.
fn bed_identity(line: &str) -> String {
    let fields: Vec<&str> = line.split('\t').collect();
    format!("{}:{}-{}:{}", fields[0], fields[1], fields[2], fields[3])
}

/// Same identity, taken from a report row.
fn row_identity(report: &Report, row: &[String]) -> String {
    format!(
        "{}:{}-{}:{}",
        report.field(row, "chrom"),
        report.field(row, "start"),
        report.field(row, "end"),
        report.field(row, "read_id")
    )
}

// ---------------------------------------------------------------------------
// Shape of the file
// ---------------------------------------------------------------------------

#[test]
fn report_is_written_next_to_the_bed_outputs() {
    let run = Run::new("--all", 4);

    assert!(run.path("test.hq.bed").is_file());
    assert!(run.path("test.scraps.bed").is_file());
    assert!(
        run.path("test.report.tsv").is_file(),
        "the report is generated without a dedicated flag"
    );
}

#[test]
fn report_header_matches_the_required_schema() {
    let run = Run::new("--all", 4);
    let report = run.report();

    assert_eq!(report.header, COLUMNS.to_vec());
    assert_eq!(report.header[report.header.len() - 2], "reasons");
    assert_eq!(report.header[report.header.len() - 1], "decision");
    assert!(
        !report.raw.starts_with('#'),
        "no comment or metadata line precedes the header"
    );
}

#[test]
fn report_has_one_row_per_query_record() {
    for mode in ["--all", "--non-overlapping"] {
        let run = Run::new(mode, 4);
        let report = run.report();

        assert_eq!(
            report.rows.len(),
            query_record_count(),
            "{mode}: one row per valid query BED12 record"
        );

        for row in &report.rows {
            assert_eq!(row.len(), COLUMNS.len(), "{mode}: every row is complete");
        }
    }
}

#[test]
fn report_never_contains_reference_records() {
    let run = Run::new("--all", 4);
    let report = run.report();

    for row in &report.rows {
        assert_ne!(
            report.field(row, "read_id"),
            "refA",
            "reference records must not appear in the report"
        );
    }
}

// ---------------------------------------------------------------------------
// Values
// ---------------------------------------------------------------------------

#[test]
fn report_decisions_are_only_pass_or_scrap() {
    for mode in ["--all", "--non-overlapping"] {
        let run = Run::new(mode, 4);
        let report = run.report();

        for row in &report.rows {
            let decision = report.field(row, "decision");
            assert!(
                decision == "PASS" || decision == "SCRAP",
                "{mode}: unexpected decision {decision}"
            );
        }
    }
}

#[test]
fn report_values_are_formatted_as_specified() {
    for mode in ["--all", "--non-overlapping"] {
        let run = Run::new(mode, 4);
        let report = run.report();

        let float_columns = [
            "best_reference_overlap",
            "junction_fraction",
            "intron_support_fraction",
            "median_splice_score",
            "weakest_splice_score",
        ];
        let bool_columns = ["overlaps_reference_exon", "boundary_match"];

        for row in &report.rows {
            for (idx, field) in row.iter().enumerate() {
                assert!(
                    !field.is_empty(),
                    "{mode}: empty field in column {}",
                    COLUMNS[idx]
                );
                assert!(
                    !field.contains("NaN") && !field.contains("inf"),
                    "{mode}: non-finite value in column {}",
                    COLUMNS[idx]
                );
                assert!(
                    !field.contains("Some(") && !field.contains('{'),
                    "{mode}: debug representation in column {}",
                    COLUMNS[idx]
                );
            }

            for column in float_columns {
                let value = report.field(row, column);
                if value != "." {
                    let decimals = value.split('.').nth(1).unwrap_or("");
                    assert_eq!(
                        decimals.len(),
                        3,
                        "{mode}: {column}={value} must use three decimals"
                    );
                }
            }

            for column in bool_columns {
                let value = report.field(row, column);
                assert!(
                    ["true", "false", "."].contains(&value.as_str()),
                    "{mode}: {column}={value} must be lowercase or '.'"
                );
            }

            let reasons = report.field(row, "reasons");
            assert!(
                !reasons.contains(' '),
                "{mode}: reasons must not hold spaces"
            );
            assert!(reasons != ".", "{mode}: every row explains its decision");

            let thresholds = report.field(row, "applied_thresholds");
            assert!(
                !thresholds.contains(' '),
                "{mode}: applied_thresholds must not hold spaces"
            );

            let read_type = report.field(row, "read_type");
            assert!(
                read_type == "SINGLE_EXON" || read_type == "MULTI_EXON",
                "{mode}: unexpected read_type {read_type}"
            );
        }
    }
}

#[test]
fn report_uses_dot_for_unavailable_values() {
    let run = Run::new("--all", 4);
    let report = run.report();

    // a guided single-exon read is scored by overlap, never by junctions
    let se_row = report
        .rows
        .iter()
        .find(|row| report.field(row, "read_id") == "q_se_pass")
        .expect("q_se_pass row");

    assert_eq!(report.field(se_row, "junction_matches"), ".");
    assert_eq!(report.field(se_row, "junction_count"), ".");
    assert_eq!(report.field(se_row, "junction_fraction"), ".");
    assert_eq!(report.field(se_row, "cluster_size"), ".");
    assert_eq!(report.field(se_row, "median_splice_score"), ".");
    assert_eq!(report.field(se_row, "intron_support_fraction"), ".");

    // a de novo read never has a best reference
    let oob_row = report
        .rows
        .iter()
        .find(|row| report.field(row, "read_id") == "q_oob0")
        .expect("q_oob0 row");

    assert_eq!(report.field(oob_row, "best_reference_id"), ".");
    assert_eq!(report.field(oob_row, "best_reference_overlap"), ".");
    assert_eq!(report.field(oob_row, "boundary_match"), ".");
}

// ---------------------------------------------------------------------------
// Branch coverage
// ---------------------------------------------------------------------------

#[test]
fn report_maps_every_guided_branch() {
    let run = Run::new("--all", 4);
    let report = run.report();
    let outcomes = report.outcomes();

    let expected: [(&str, &str, &str); 6] = [
        ("q_junction_pass", "JUNCTION_SUPPORT_PASS", "PASS"),
        ("q_se_pass", "REFERENCE_OVERLAP_PASS", "PASS"),
        // a complete junction shared with one reference rescues a low junction fraction
        (
            "q_junction_rescue",
            "LOW_JUNCTION_SUPPORT/BOUNDARY_MATCH_RESCUE",
            "PASS",
        ),
        // reference evidence is insufficient, and the other reads do not corroborate it
        (
            "q_se_scrap",
            "LOW_REFERENCE_OVERLAP/LOW_SINGLE_EXON_CLUSTER_SUPPORT",
            "SCRAP",
        ),
        (
            "q_single_coordinate",
            "LOW_JUNCTION_SUPPORT/NO_BOUNDARY_MATCH/LOW_INTRON_CHAIN_CLUSTER_SUPPORT/LOW_INTRON_SUPPORT/SPLICE_SCORE_UNAVAILABLE",
            "SCRAP",
        ),
        (
            "q_low_junction",
            "LOW_JUNCTION_SUPPORT/NO_BOUNDARY_MATCH/LOW_INTRON_CHAIN_CLUSTER_SUPPORT/LOW_INTRON_SUPPORT/SPLICE_SCORE_UNAVAILABLE",
            "SCRAP",
        ),
    ];

    for (read, reasons, decision) in expected {
        let rows = outcomes
            .get(read)
            .unwrap_or_else(|| panic!("no row for {read}"));
        assert_eq!(
            rows,
            &vec![(reasons.to_string(), decision.to_string())],
            "unexpected outcome for {read}"
        );
    }

    // redirected reads keep the contextual reason first
    for i in 0..5 {
        let read = format!("q_oob{}", i);
        let rows = outcomes.get(&read).expect("redirected row");
        assert_eq!(
            rows,
            &vec![(
                "NO_REFERENCE_EXON_OVERLAP/DOMINANT_INTRON_CHAIN_CLUSTER".to_string(),
                "PASS".to_string()
            )],
            "unexpected outcome for {read}"
        );
    }

    let redirected = report
        .rows
        .iter()
        .find(|row| report.field(row, "read_id") == "q_oob0")
        .expect("q_oob0 row");
    assert_eq!(report.field(redirected, "mode"), "GUIDED");
    assert_eq!(
        report.field(redirected, "evaluation_path"),
        "GUIDED_TO_DE_NOVO"
    );
    assert_eq!(report.field(redirected, "overlaps_reference_exon"), "false");
    assert_eq!(report.field(redirected, "reference_count"), "1");
    assert_eq!(
        report.field(redirected, "group_size"),
        "13",
        "the de novo evidence base is the whole component, not the redirected subset"
    );
    assert_eq!(report.field(redirected, "cluster_size"), "5");

    // the rescue names the one reference that carried the whole feature
    let rescued = report
        .rows
        .iter()
        .find(|row| report.field(row, "read_id") == "q_junction_rescue")
        .expect("q_junction_rescue row");
    assert_eq!(report.field(rescued, "boundary_match"), "true");
    assert_eq!(
        report.field(rescued, "boundary_evidence"),
        "SPLICE_JUNCTION"
    );
    assert_eq!(report.field(rescued, "best_reference_id"), "refA");
    assert_eq!(report.field(rescued, "evaluation_path"), "GUIDED_REFERENCE");

    // a single shared coordinate is not a rescue, and the read continues to de novo
    let coordinate = report
        .rows
        .iter()
        .find(|row| report.field(row, "read_id") == "q_single_coordinate")
        .expect("q_single_coordinate row");
    assert_eq!(report.field(coordinate, "boundary_match"), "false");
    assert_eq!(report.field(coordinate, "boundary_evidence"), ".");
    assert_eq!(
        report.field(coordinate, "evaluation_path"),
        "GUIDED_TO_DE_NOVO"
    );
    assert_eq!(
        report.field(coordinate, "overlaps_reference_exon"),
        "true",
        "it overlaps a reference exon but its reference evidence was insufficient"
    );

    let guided = report
        .rows
        .iter()
        .find(|row| report.field(row, "read_id") == "q_junction_pass")
        .expect("q_junction_pass row");
    assert_eq!(report.field(guided, "mode"), "GUIDED");
    assert_eq!(report.field(guided, "evaluation_path"), "GUIDED_REFERENCE");
    assert_eq!(report.field(guided, "overlaps_reference_exon"), "true");
    assert_eq!(report.field(guided, "best_reference_id"), "refA");
    assert_eq!(report.field(guided, "junction_matches"), "2");
    assert_eq!(report.field(guided, "junction_count"), "2");
    assert_eq!(report.field(guided, "junction_fraction"), "1.000");
    assert_eq!(
        report.field(guided, "applied_thresholds"),
        "min_junction_frac=0.500"
    );
    assert_eq!(report.field(guided, "exon_count"), "3");
    assert_eq!(report.field(guided, "read_type"), "MULTI_EXON");
    assert_eq!(report.field(guided, "strand"), "+");
}

#[test]
fn report_maps_de_novo_branches() {
    let run = Run::new("--non-overlapping", 4);
    let report = run.report();
    let outcomes = report.outcomes();

    for row in &report.rows {
        assert_eq!(report.field(row, "mode"), "DE_NOVO");
        assert_eq!(report.field(row, "evaluation_path"), "DE_NOVO");
        assert_eq!(
            report.field(row, "overlaps_reference_exon"),
            ".",
            "reference-exon overlap is not applicable in de novo mode"
        );
        assert_eq!(report.field(row, "reference_count"), "0");
    }

    // the lone read is still rejected, but only after every criterion was evaluated:
    // component size alone no longer decides it
    assert_eq!(
        outcomes.get("q_lone").expect("q_lone row"),
        &vec![(
            "SINGLE_READ_COMPONENT/LOW_INTRON_CHAIN_CLUSTER_SUPPORT/INTRON_SUPPORT_UNAVAILABLE/SPLICE_SCORE_UNAVAILABLE"
                .to_string(),
            "SCRAP".to_string()
        )]
    );

    // single-exon reads are clustered by reciprocal overlap
    for read in ["q_se_pass", "q_se_scrap"] {
        let rows = outcomes.get(read).expect("single-exon row");
        assert_eq!(rows.len(), 1);
        assert!(
            rows[0].0 == "SUPPORTED_SINGLE_EXON_CLUSTER"
                || rows[0].0 == "LOW_SINGLE_EXON_CLUSTER_SUPPORT",
            "{read}: unexpected reasons {}",
            rows[0].0
        );
    }

    // the five identical chains form a dominant intron-chain cluster
    for i in 0..5 {
        let read = format!("q_oob{}", i);
        assert_eq!(
            outcomes.get(&read).expect("cluster row"),
            &vec![(
                "DOMINANT_INTRON_CHAIN_CLUSTER".to_string(),
                "PASS".to_string()
            )],
            "unexpected outcome for {read}"
        );
    }
}

// ---------------------------------------------------------------------------
// Agreement with the BED outputs
// ---------------------------------------------------------------------------

#[test]
fn report_decisions_agree_with_the_bed_outputs() {
    for mode in ["--all", "--non-overlapping"] {
        let run = Run::new(mode, 4);
        let report = run.report();

        let hq: Vec<String> = run
            .bed_lines("test.hq.bed")
            .iter()
            .map(|line| bed_identity(line))
            .collect();
        let scraps: Vec<String> = run
            .bed_lines("test.scraps.bed")
            .iter()
            .map(|line| bed_identity(line))
            .collect();

        let mut pass = 0usize;
        let mut scrap = 0usize;

        for row in &report.rows {
            let identity = row_identity(&report, row);
            match report.field(row, "decision").as_str() {
                "PASS" => {
                    pass += 1;
                    assert!(
                        hq.contains(&identity),
                        "{mode}: {identity} is PASS but missing from the HQ BED"
                    );
                    assert!(
                        !scraps.contains(&identity),
                        "{mode}: {identity} is PASS but present in the scraps BED"
                    );
                }
                "SCRAP" => {
                    scrap += 1;
                    assert!(
                        scraps.contains(&identity),
                        "{mode}: {identity} is SCRAP but missing from the scraps BED"
                    );
                    assert!(
                        !hq.contains(&identity),
                        "{mode}: {identity} is SCRAP but present in the HQ BED"
                    );
                }
                other => panic!("{mode}: unexpected decision {other}"),
            }
        }

        // and the other way around: no BED line is missing from the report
        let pass_ids: Vec<String> = report
            .rows
            .iter()
            .filter(|row| report.field(row, "decision") == "PASS")
            .map(|row| row_identity(&report, row))
            .collect();
        let scrap_ids: Vec<String> = report
            .rows
            .iter()
            .filter(|row| report.field(row, "decision") == "SCRAP")
            .map(|row| row_identity(&report, row))
            .collect();

        for identity in &hq {
            assert!(
                pass_ids.contains(identity),
                "{mode}: {identity} is in the HQ BED without a PASS row"
            );
        }
        for identity in &scraps {
            assert!(
                scrap_ids.contains(identity),
                "{mode}: {identity} is in the scraps BED without a SCRAP row"
            );
        }

        // the BED sets deduplicate identical lines, the report does not
        assert_eq!(
            hq.len(),
            unique(&hq).len(),
            "{mode}: the HQ BED holds no duplicate line"
        );
        assert_eq!(
            scraps.len(),
            unique(&scraps).len(),
            "{mode}: the scraps BED holds no duplicate line"
        );
        assert!(
            pass >= hq.len(),
            "{mode}: PASS rows ({pass}) cover every HQ line ({})",
            hq.len()
        );
        assert!(
            scrap >= scraps.len(),
            "{mode}: SCRAP rows ({scrap}) cover every scraps line ({})",
            scraps.len()
        );
        assert_eq!(pass + scrap, report.rows.len());
    }
}

fn unique(values: &[String]) -> Vec<String> {
    let mut seen: Vec<String> = Vec::new();
    for value in values {
        if !seen.contains(value) {
            seen.push(value.clone());
        }
    }
    seen
}

#[test]
fn report_keeps_identical_duplicate_records_as_separate_rows() {
    let run = Run::new("--all", 4);
    let report = run.report();
    let outcomes = report.outcomes();

    let duplicates = outcomes.get("q_dup").expect("q_dup rows");
    assert_eq!(
        duplicates.len(),
        2,
        "byte-identical query records keep one report row each"
    );
    assert_eq!(duplicates[0], duplicates[1]);

    let hq = run.bed_lines("test.hq.bed");
    let dup_lines = hq.iter().filter(|line| line.contains("q_dup")).count();
    assert_eq!(
        dup_lines, 1,
        "the HQ BED keeps its existing deduplication behavior"
    );
}

// ---------------------------------------------------------------------------
// Determinism
// ---------------------------------------------------------------------------

#[test]
fn report_is_identical_with_one_and_many_threads() {
    for mode in ["--all", "--non-overlapping"] {
        let single = Run::new(mode, 1);
        let many = Run::new(mode, 8);

        let single_report = fs::read_to_string(single.path("test.report.tsv")).unwrap();
        let many_report = fs::read_to_string(many.path("test.report.tsv")).unwrap();

        assert_eq!(
            single_report, many_report,
            "{mode}: thread scheduling must not change the report"
        );
    }
}

#[test]
fn report_rows_are_sorted_by_the_documented_key() {
    let run = Run::new("--all", 4);
    let report = run.report();

    let keys: Vec<(String, u64, u64, String, String, String)> = report
        .rows
        .iter()
        .map(|row| {
            (
                report.field(row, "chrom"),
                report.field(row, "start").parse().unwrap(),
                report.field(row, "end").parse().unwrap(),
                report.field(row, "read_id"),
                report.field(row, "decision"),
                report.field(row, "reasons"),
            )
        })
        .collect();

    let mut sorted = keys.clone();
    sorted.sort();

    assert_eq!(keys, sorted, "rows follow chrom, start, end, read_id, ...");
}
