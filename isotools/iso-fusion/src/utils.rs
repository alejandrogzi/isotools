use std::collections::BTreeSet;
use std::fmt::Debug;
use std::path::Path;

use anyhow::{Ok, Result};
use dashmap::DashMap;
use hashbrown::{HashMap, HashSet};
use log::info;
use packbed::{packbed, par_reader, Bed12, GenePred};
use rayon::prelude::*;

use config::get_progress_bar;

pub fn unpack_blacklist<P: AsRef<Path> + Debug + Sync + Send>(
    files: Vec<P>,
    cds_overlap: bool,
    is_ref: bool,
) -> Result<HashMap<String, HashSet<String>>, anyhow::Error> {
    let contents = par_reader(files)?;
    let tracks = parse_tracks(&contents, cds_overlap, is_ref)?;

    Ok(tracks)
}

pub fn prepare_refs<P: AsRef<Path> + Debug + Sync + Send>(
    refs: Vec<P>,
) -> Result<String, anyhow::Error> {
    let tracks = packbed(refs, None, true, false)?;

    let results: Vec<String> = tracks
        .into_par_iter()
        .map(|bucket| {
            let components = bucket.1;
            let mut acc = String::new();

            for comp in components {
                let refs = comp.0;
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
    cds_overlap: bool,
    is_ref: bool,
) -> Result<HashMap<String, HashSet<String>>, anyhow::Error> {
    let pb = get_progress_bar(
        contents.lines().count() as u64,
        "Parsing BED12 blacklist...",
    );
    let tracks = contents
        .par_lines()
        .filter(|x| !x.starts_with("#"))
        .filter_map(|x| Bed12::parse(x, cds_overlap, is_ref).ok())
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

#[allow(dead_code)]
fn get_intergenic_regions(
    refs: HashMap<String, Vec<GenePred>>,
    cds: bool,
) -> Result<DashMap<String, BTreeSet<(u64, u64, char)>>> {
    let pb = get_progress_bar(
        refs.values().flatten().count() as u64,
        "Parsing intergenic regions...",
    );

    let regions = refs
        .into_par_iter()
        .fold(
            || DashMap::new(),
            |acc, record| {
                let chr = record.0;
                let bucket = record.1;
                let mut ranges = BTreeSet::new();

                let (plus, minus): (Vec<_>, Vec<_>) = bucket.iter().partition(|t| t.strand == '+');

                for transcripts in [plus, minus] {
                    if transcripts.is_empty() {
                        continue;
                    }

                    let mut regions: BTreeSet<(u64, u64, char)> = BTreeSet::new();
                    let mut current = 0;
                    let mut next = current + 1;
                    let mut last_region_start = 0;

                    while current < transcripts.len() - 1 {
                        // 1. establish region boundaries
                        //    - left edge starts being the end of the first transcript

                        let mut region_start = if cds {
                            transcripts[current].cds_end
                        } else {
                            transcripts[current].end
                        };
                        let mut region_end = if cds {
                            transcripts[next].cds_start
                        } else {
                            transcripts[next].start
                        };
                        if current < 1 {
                            last_region_start = region_start
                        };

                        // 2. if transcript is inside next or prevoius transcript
                        //   skip to the next transcript

                        if next < transcripts.len() - 1 {
                            if current > 0 {
                                if (transcripts[next].start < transcripts[current].end
                                    && transcripts[current].end < transcripts[next].end)
                                    || (transcripts[current - 1].start < transcripts[current].end
                                        && transcripts[current].end < transcripts[current - 1].end)
                                {
                                    current += 1;
                                    next = current + 1;

                                    continue;
                                }
                            } else {
                                if transcripts[next].start < transcripts[current].end
                                    && transcripts[current].end < transcripts[next].end
                                {
                                    current += 1;
                                    next = current + 1;

                                    continue;
                                }
                            }
                        }

                        // 3. check if region boundaries are valid, if they are not
                        //    increase next until finding a valid region
                        //    - region start should be less than region end
                        //    - region start should be greater than left edge

                        while region_start > region_end && next < transcripts.len() - 1 {
                            next += 1;

                            region_end = if cds {
                                transcripts[next].cds_start
                            } else {
                                transcripts[next].start
                            };

                            if region_start < last_region_start {
                                current += 1;
                                next = current + 1;
                                if current >= transcripts.len() - 1 {
                                    break;
                                }

                                region_start = if cds {
                                    transcripts[current].cds_end
                                } else {
                                    transcripts[current].end
                                };
                                region_end = if cds {
                                    transcripts[next].cds_start
                                } else {
                                    transcripts[next].start
                                };
                                last_region_start = region_start;
                            }
                        }

                        // 3. check if region boundaries are valid

                        if region_start < last_region_start || region_start > region_end {
                            if current + 1 >= transcripts.len() - 1 {
                                break;
                            }

                            current += 1;
                            next = current + 1;
                            continue;
                        }

                        // 4. if next < transcripts.len() - 1, iterate over that space and check
                        //    if intergenic region is inside any of the transcripts

                        let mut is_retained = false;
                        if next + 1 < transcripts.len() - 1 {
                            for i in current + 1..next + 1 {
                                if transcripts[i].start < region_start
                                    && region_end < transcripts[i].end
                                {
                                    is_retained = true;
                                    break;
                                }
                            }
                        }

                        last_region_start = region_start;

                        if !is_retained {
                            regions.insert((region_start, region_end, transcripts[current].strand));
                        }

                        current += 1;
                        next = current + 1;
                        pb.inc(1);
                    }

                    ranges.extend(regions);
                }

                acc.insert(chr, ranges);
                acc
            },
        )
        .reduce(
            || DashMap::new(),
            |map1, map2| {
                for entry in map2.into_iter() {
                    map1.entry(entry.0.clone())
                        .or_insert_with(BTreeSet::new)
                        .extend(entry.1.iter().cloned());
                }
                map1
            },
        );

    pb.finish_and_clear();
    info!(
        "Intergenic regions parsed: {}",
        regions.iter().map(|chr| chr.value().len()).sum::<usize>()
    );

    Ok(regions)
}
