use std::fmt::Debug;
use std::path::Path;
use std::str::FromStr;

use anyhow::{Ok, Result};
use hashbrown::{HashMap, HashSet};
use log::info;
use packbed::{packbed, par_reader, Bed12, GenePred, RefGenePred};
use rayon::prelude::*;

use config::{get_progress_bar, OverlapType, TsvParser};

pub fn unpack_blacklist<P: AsRef<Path> + Debug + Sync + Send>(
    files: Vec<P>,
    cds_overlap: OverlapType,
    is_ref: bool,
) -> Result<HashMap<String, HashSet<String>>, anyhow::Error> {
    let contents = par_reader(files)?;
    let tracks = parse_tracks(&contents, cds_overlap, is_ref)?;

    Ok(tracks)
}

pub fn prepare_refs<P: AsRef<Path> + Debug + Sync + Send>(
    refs: Vec<P>,
) -> Result<String, anyhow::Error> {
    let tracks = packbed(
        refs,
        None,
        config::OverlapType::Exon,
        packbed::PackMode::Default,
    )?;

    let results: Vec<String> = tracks
        .into_par_iter()
        .map(|bucket| {
            let components = bucket.1;
            let mut acc = String::new();

            for comp in components {
                let refs = comp
                    .as_any()
                    .downcast_ref::<(RefGenePred, Vec<GenePred>)>()
                    .expect("ERROR: Failed to downcast to GenePred")
                    .0
                    .clone();

                let loci = refs.merge_names();
                refs.reads.into_iter().for_each(|mut record| {
                    let mut line = record.line.clone();

                    // if loci has more than 1 gene name -> fusion
                    if loci.split('.').count() > 1 {
                        line = record.mut_name_from_line(&loci);
                    }

                    acc += format!("{}\n", line).as_str();
                });
            }

            acc
        })
        .collect();

    info!("Reference transcripts fixed for fusions!");
    Ok(results.concat())
}

fn parse_tracks<'a>(
    contents: &'a str,
    cds_overlap: OverlapType,
    is_ref: bool,
) -> Result<HashMap<String, HashSet<String>>, anyhow::Error> {
    let pb = get_progress_bar(
        contents.lines().count() as u64,
        "Parsing BED12 blacklist...",
    );
    let tracks = contents
        .par_lines()
        .filter(|row| !row.starts_with("#"))
        .filter_map(|row| Bed12::read(row, cds_overlap, is_ref).ok())
        .fold(
            || HashMap::new(),
            |mut acc: HashMap<String, HashSet<String>>, record| {
                acc.entry(record.chrom.clone())
                    .or_default()
                    .insert(record.name);
                pb.inc(1);
                acc
            },
        )
        .reduce(
            || HashMap::new(),
            |mut acc, map| {
                for (k, v) in map {
                    let acc_v = acc.entry(k).or_insert(HashSet::new());
                    acc_v.extend(v);
                }
                acc
            },
        );

    pb.finish_and_clear();
    info!("Records parsed: {}", tracks.values().flatten().count());

    Ok(tracks)
}

#[derive(Debug)]
pub struct IsoformParser {
    fields: Vec<String>,
}

impl TsvParser for IsoformParser {
    fn parse(line: &str) -> Result<Self, anyhow::Error> {
        Ok(Self {
            fields: line.split('\t').map(|s| s.to_string()).collect(),
        })
    }

    fn key(&self, index: usize) -> &str {
        &self.fields[index]
    }

    fn value<V: FromStr>(&self, index: usize) -> Result<V, V::Err> {
        self.fields[index].parse::<V>()
    }
}
