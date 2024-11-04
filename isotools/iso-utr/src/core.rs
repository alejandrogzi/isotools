use anyhow::Result;
use config::{get_progress_bar, write_objs, COLORIZE, OVERLAP};
use dashmap::DashSet;
use hashbrown::HashSet;
use log::info;
use packbed::{packbed, GenePred};
use rayon::prelude::*;
use std::sync::Arc;

use crate::cli::Args;
use crate::utils::unpack_blacklist;

pub fn detect_truncations(args: Args) -> Result<()> {
    let tracks = packbed(args.refs, args.query, OVERLAP, COLORIZE)?;
    let blacklist = unpack_blacklist(args.blacklist).unwrap_or_default();

    let hit_acc: DashSet<&String> = DashSet::new();
    let pass_acc: DashSet<&String> = DashSet::new();
    let dirt_acc: DashSet<String> = DashSet::new();

    info!("Detecting intron retentions...");
    let pb = get_progress_bar(tracks.len() as u64, "Processing...");

    tracks.par_iter().for_each(|(chr, buckets)| {
        let binding = HashSet::new();

        buckets.par_iter().for_each(|bucket| {
            // let (hits, pass, dirt) = process_bucket(bucket, &blacklist);

            // hits.iter().for_each(|hit| {
            //     hit_acc.insert(hit);
            // });
            // pass.iter().for_each(|p| {
            //     pass_acc.insert(p);
            // });

            // if args.recover {
            //     // let (new_hits, new_pass) = recover_from_dirt(dirt, &blacklist);

            //     new_hits.iter().for_each(|hit| {
            //         hit_acc.insert(hit);
            //     });

            //     new_pass.iter().for_each(|p| {
            //         pass_acc.insert(p);
            //     });
            // };
        });

        pb.inc(1);
    });

    pb.finish_and_clear();
    info!("Reads with retained introns: {}", hit_acc.len());

    [hit_acc, pass_acc]
        .par_iter()
        .zip([HIT, PASS, DIRTY].par_iter())
        .for_each(|(rx, path)| write_objs(&rx, path));

    Ok(())
}

pub fn process_bucket<'a>(bucket: &'a Vec<Arc<GenePred>>, ban: &HashSet<String>) {
    todo!();
}

pub fn recover_from_dirt(dirt: HashSet<String>, ban: &HashSet<String>) {
    todo!();
}
