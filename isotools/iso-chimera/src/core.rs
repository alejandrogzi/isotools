use anyhow::Result;
use config::{get_progress_bar, write_objs, CHIMERAS, CHIMERIC_FREE, INTERGENIC_REGIONS};
use dashmap::{DashMap, DashSet};
use hashbrown::HashMap;
use log::info;
use packbed::{unpack, GenePred};
use rayon::prelude::*;

use std::collections::BTreeSet;
use std::io::Write;

use crate::cli::Args;

const CHIMERA_CDS_OVERLAP: bool = true;
const CHIMERA_REF: bool = false;

pub fn detect_chimeras(args: Args) -> Result<()> {
    let reads = unpack(args.query, CHIMERA_CDS_OVERLAP, CHIMERA_REF)?;
    let regions = if !args.hint.is_empty() {
        let mut refs = unpack(args.hint, CHIMERA_CDS_OVERLAP, CHIMERA_REF)?;
        refs.par_iter_mut().for_each(|(_, v)| {
            v.par_sort_unstable_by(|a, b| a.start.cmp(&b.start));
        });

        get_intergenic_regions(refs, args.cds).ok()
    } else {
        None
    };

    if args.write {
        let mut count = 0;
        if let Some(regions) = regions.as_ref() {
            let mut bed = std::fs::File::create(INTERGENIC_REGIONS)?;

            for bucket in regions.iter() {
                let chr = bucket.key();
                let regions = bucket.value();

                for (mut start, mut end, strand) in regions.into_iter() {
                    if *strand == '-' {
                        let tmp = start.clone();
                        start = config::SCALE - end;
                        end = config::SCALE - tmp;
                    }

                    writeln!(
                        bed,
                        "{}\t{}\t{}\t{}\t{}\t{}",
                        chr, start, end, count, "0", strand
                    )?;
                    count += 1;
                }
            }
        }
    }

    if let Some(regions) = regions.as_ref() {
        info!("Detecting chimeras...");
        let chimeras = DashSet::new();
        let pb = get_progress_bar(
            reads.values().flatten().count() as u64,
            "Detecting chimeras...",
        );

        reads.par_iter().for_each(|(chr, bucket)| {
            for read in bucket {
                let read_start = read.start;
                let read_end = read.end;
                let read_strand = read.strand;

                if let Some(intergenic) = regions.get(chr) {
                    for (region_start, region_end, region_strand) in intergenic.iter() {
                        if read_strand != *region_strand {
                            continue;
                        }

                        // TODO: check how make this work -> optimize
                        // if region_start > &read_end || region_end < &read_start {
                        //     break;
                        // }

                        if &read_start < region_start && region_end < &read_end {
                            chimeras.insert(read.line().to_owned());

                            break;
                        }
                    }
                }
            }

            pb.inc(1);
        });

        pb.finish_and_clear();
        info!("Detected {} chimeras", chimeras.len());

        let clean: HashMap<_, Vec<_>> = reads
            .into_iter()
            .map(|(chr, bucket)| {
                let filtered_bucket: Vec<_> = bucket
                    .into_iter()
                    .filter(|read| !chimeras.contains(read.line()))
                    .collect();
                (chr, filtered_bucket)
            })
            .collect();

        info!(
            "Chimeric-free dataset size: {}",
            clean.values().flatten().count()
        );

        let no_chimeras = DashSet::new();
        for (_, bucket) in clean.iter() {
            for read in bucket {
                no_chimeras.insert(read.line().to_owned());
            }
        }

        [&chimeras, &no_chimeras]
            .par_iter()
            .zip([CHIMERAS, CHIMERIC_FREE].par_iter())
            .for_each(|(rx, path)| write_objs(&rx, path));
    } else {
        info!("No hints provided, detecting chimeras on frequency mode...");
        let pb = get_progress_bar(
            reads.values().flatten().count() as u64,
            "Detecting chimeras...",
        );

        todo!()
    };

    Ok(())
}

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
