use config::{bed_to_nested_map, get_progress_bar, BedColumn, BedColumnValue};
use dashmap::{DashMap, DashSet};
use hashbrown::HashMap;
use packbed::{
    packbed, reader,
    record::{Bed6, PolyAPred},
    BedPackage, PackMode,
};
use rayon::prelude::*;

use std::{
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
    sync::{
        atomic::{AtomicUsize, Ordering},
        Arc,
    },
};

use crate::cli::CallerArgs;

// 3. pack filter result -> new struct that includes polyA tail length

pub fn pas_caller(args: CallerArgs) -> Result<(), Box<dyn std::error::Error>> {
    let isoseqs = packbed(
        vec![args.bed],
        None,
        config::OverlapType::Exon,
        PackMode::PolyA,
    )?;

    // WARN: need to fix duplicated aparent predictions here!
    let aparent_scores = get_aparent_scores(args.aparent);

    let pb = get_progress_bar(isoseqs.len() as u64, "Processing reads...");
    let accumulator = ParallelAccumulator::default();

    isoseqs.into_par_iter().for_each(|(chr, reads)| {
        if let Some(aparent) = aparent_scores.get(&chr) {
            distribute(reads, &*aparent, &accumulator); // Pass reference
        } else {
            panic!("ERROR: Could not get aparent scores for chromosome!");
        }

        pb.inc(1);
    });

    pb.finish_and_clear();

    Ok(())
}

fn get_aparent_scores(file: PathBuf) -> DashMap<String, HashMap<String, BedColumnValue>> {
    let content = reader(file).expect("ERROR: Could not read aparent file!");
    let aparent_scores = bed_to_nested_map::<Bed6>(Arc::new(content), BedColumn::Score)
        .expect("ERROR: Could not build mapper from aparent file!");

    aparent_scores
}

#[inline(always)]
fn distribute(
    components: Vec<Box<dyn BedPackage>>,
    scores: &HashMap<String, BedColumnValue>,
    accumulator: &ParallelAccumulator,
) {
    components.into_par_iter().for_each(|comp| {
        let comp = comp
            .as_any()
            .downcast_ref::<Vec<PolyAPred>>()
            .expect("ERROR: Could not downcast to PolyAPred!");

        // let info = process_component(comp, banned, splice_map, scan_scores, genome, nag);

        // info.into_iter().for_each(|(_, intron_descriptor)| {
        //     if !intron_descriptor.is_empty() {
        //         accumulator.introns.insert(intron_descriptor);
        //     }
        // });
    });
}

struct ParallelAccumulator {
    lines: DashSet<String>,
    paths: DashSet<String>,
}

impl Default for ParallelAccumulator {
    fn default() -> Self {
        Self {
            lines: DashSet::new(),
            paths: DashSet::new(),
        }
    }
}
