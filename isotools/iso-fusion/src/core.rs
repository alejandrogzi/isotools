use anyhow::Result;
use config::{
    get_progress_bar, write_objs, ModuleMap, FUSIONS, FUSION_FREE, OVERLAP_CDS, OVERLAP_EXON,
};
use dashmap::DashSet;
use hashbrown::{HashMap, HashSet};
use log::info;
use packbed::{buckerize, combine, parse_tracks, unpack, GenePred, RefGenePred};
use rayon::prelude::*;

use std::sync::atomic::{AtomicU32, Ordering};

use crate::cli::Args;
use crate::utils::{exonic_overlap, prepare_refs, unpack_blacklist};

pub fn detect_fusions(args: Args) -> Result<()> {
    info!("Preparing files for fusion detection...");
    let ref_str = prepare_refs(args.refs)?;
    let refs = parse_tracks(ref_str.as_str(), OVERLAP_CDS, true)
        .expect("Failed to parse reference transcripts");
    let query = unpack(args.query, OVERLAP_CDS, false).expect("Failed to unpack query tracks");
    let (tracks, n) = combine(refs, query);
    let buckets = buckerize(tracks, OVERLAP_CDS, OVERLAP_EXON, n);

    let blacklist = if !args.blacklist.is_empty() {
        unpack_blacklist(args.blacklist, true, false).ok()
    } else {
        None
    };
    let pb = get_progress_bar(buckets.len() as u64, "Detecting fusions...");

    let n_comps = AtomicU32::new(0);
    let dirty_count = AtomicU32::new(0);

    let hits: DashSet<String> = DashSet::new();
    let passes: DashSet<String> = DashSet::new();

    buckets.par_iter().for_each(|bucket| {
        let chr = bucket.key();
        let components = bucket.value().to_owned();
        n_comps.fetch_add(components.len() as u32, Ordering::Relaxed);

        let binding = HashSet::new();
        let banned = if let Some(bl) = blacklist.as_ref() {
            bl.get(chr).unwrap_or(&binding)
        } else {
            &binding
        };

        components.into_par_iter().for_each(|comp| {
            let (fusions, no_fusions, descriptor, is_dirty) =
                process_component(comp, banned, args.recover);

            fusions.into_iter().for_each(|hit| {
                hits.insert(hit);
            });
            no_fusions.into_iter().for_each(|p| {
                passes.insert(p);
            });

            if is_dirty {
                dirty_count.fetch_add(1, Ordering::Relaxed);
            }
        });

        pb.inc(1);
    });

    pb.finish_and_clear();
    info!("Detected fusions: {}", hits.len());
    info!("Fusion-free reads: {}", passes.len());

    [&hits, &passes]
        .par_iter()
        .zip([FUSIONS, FUSION_FREE].par_iter())
        .for_each(|(rx, path)| write_objs(&rx, path));

    Ok(())
}

fn process_component(
    component: (RefGenePred, Vec<GenePred>),
    banned: &HashSet<String>,
    recover: bool,
) -> (
    Vec<String>,
    Vec<String>,
    HashMap<String, Box<dyn ModuleMap>>,
    bool,
) {
    let mut descriptor = HashMap::new();

    let mut fusions = Vec::new();
    let mut no_fusions = Vec::new();

    let mut is_dirty = false;

    let refs = component.0;
    let queries = component.1;

    let genes = refs.get_names();
    if genes.len() > 1 {
        // fusion
        // log::warn!("Multiple ref genes in component: {:?}", genes);
        let refs = refs.smash_exons_by_name();

        queries.into_iter().for_each(|query| {
            if banned.contains(&query.name) {
                no_fusions.push(query.line);
                return;
            }

            let mut count = 0;

            // check if query read overlaps any of the genes
            for ref_exons in refs.iter() {
                if exonic_overlap(&query.exons, &ref_exons) {
                    count += 1;
                }
            }

            if count > 1 {
                fusions.push(query.line);
            } else {
                no_fusions.push(query.line);
            }
        });
    } else if genes.is_empty() {
        // species-specific gene | intergenic region | missing gene in refs
        // log::warn!("No ref genes in component, check loci: {:?}", queries);
        queries.into_iter().for_each(|query| {
            no_fusions.push(query.line);
        });
    } else {
        // all good, no fusions here
        queries.into_iter().for_each(|query| {
            no_fusions.push(query.line);
        });
    }

    (fusions, no_fusions, descriptor, is_dirty)
}
