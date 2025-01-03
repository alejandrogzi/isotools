use anyhow::Result;
use dashmap::DashSet;
use log::{info, warn};
use packbed::{packbed, GenePred, RefGenePred};
use rayon::prelude::*;

use std::sync::atomic::AtomicUsize;

use crate::cli::Args;
use config::{exonic_overlap, get_progress_bar, write_objs, OVERLAP_CDS, OVERLAP_EXON};

type Component = (RefGenePred, Vec<GenePred>);
type Components = Vec<Component>;

const COVERAGE: &str = "coverage.tsv";

pub fn calculate_coverage(args: Args) -> Result<()> {
    info!("Calculating coverage of your Iso-Seqs...");

    let tracks = packbed(args.refs, Some(args.query), OVERLAP_CDS, OVERLAP_EXON)?;
    let pb = get_progress_bar(tracks.len() as u64, "Processing...");

    let accumulator = ParallelAccumulator::default();

    tracks.par_iter().for_each(|bucket| {
        let components = bucket.value().to_owned();

        process_components(components, &accumulator);

        pb.inc(1);
    });

    pb.finish_and_clear();
    info!(
        "Total number of unassigned queries: {}",
        accumulator
            .unassigned
            .load(std::sync::atomic::Ordering::Relaxed)
    );

    write_results(&accumulator);

    warn!("Coverage file written to: {}", COVERAGE);
    info!("Done!");

    Ok(())
}

fn process_components(components: Components, accumulator: &ParallelAccumulator) {
    components.into_par_iter().for_each(|comp| {
        let (rows, unassigned) = get_component_coverage(comp);

        rows.into_iter().for_each(|r| {
            accumulator.coverage.insert(r);
        });

        accumulator
            .unassigned
            .fetch_add(unassigned, std::sync::atomic::Ordering::Relaxed);
    });
}

// should return the coverage of each ref gene in the component:
// coverage = (number of queries that overlap with the ref gene)
fn get_component_coverage(component: Component) -> (DashSet<String>, usize) {
    let acc = DashSet::new();
    let mut unassigned = 0;

    let refs = component.0;
    let queries = component.1;

    let genes = refs.get_names_split();

    if genes.len() > 1 {
        // iterate over all queries and check if
        // they overlap with the each one of the refs
        let chr = refs.reads[0].chrom.as_str();
        let ref_groups = refs.smash_exons_by_name_split();
        for (ref_gene, ref_gene_exons) in ref_groups.iter() {
            let mut hits = 0;
            for query in queries.iter() {
                let query_exons = &query.exons;

                if exonic_overlap(ref_gene_exons, query_exons) {
                    // assign query to ref_gene
                    hits += 1;
                }
            }

            let row = format!("{}\t{}\t{}", chr, ref_gene, hits);
            acc.insert(row);
        }
    } else if genes.is_empty() {
        // species-specific | intergenic | background transcription
        unassigned = queries.len();
    } else {
        // single gene -> assign all queries to this gene
        let chr = refs.reads[0].chrom.as_str();
        let gene = genes
            .iter()
            .next()
            .expect("ERROR: No gene found. This is a bug!");
        // let name: Vec<&str> = gene.splitn(3, '.').collect();
        let hits = queries.len();
        let row = format!("{}\t{}\t{}", chr, gene, hits);
        acc.insert(row);
    }

    (acc, unassigned)
}

struct ParallelAccumulator {
    coverage: DashSet<String>,
    unassigned: AtomicUsize,
}

impl Default for ParallelAccumulator {
    fn default() -> Self {
        Self {
            coverage: DashSet::new(),
            unassigned: AtomicUsize::new(0),
        }
    }
}

fn write_results(accumulator: &ParallelAccumulator) {
    [&accumulator.coverage]
        .par_iter()
        .zip([COVERAGE].par_iter())
        .for_each(|(acc, path)| write_objs(acc, path));
}
