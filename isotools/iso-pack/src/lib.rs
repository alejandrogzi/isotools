use std::fmt::Debug;
use std::fs::File;
use std::io::Read;
use std::path::Path;
use std::sync::{Arc, Mutex};

use config::get_progress_bar;
use dashmap::DashMap;
use hashbrown::HashMap;
use log::info;
use rand::Rng;
use rayon::prelude::*;

pub mod record;
pub use record::{Bed12, GenePred};

pub type GenePredMap = HashMap<String, Vec<GenePred>>;

pub const RGB: [&str; 10] = [
    "255,0,0",    // red
    "0,255,0",    // green
    "0,0,255",    // blue
    "58,134,47",  // dark-green
    "255,0,255",  // magenta
    "0,255,255",  // cyan
    "255,128,0",  // orange
    "51,153,255", // sky-blue
    "118,115,15", // dark-yellow
    "172,126,0",  // brown
];

fn reader<P: AsRef<Path> + Debug>(file: P) -> Result<String, Box<dyn std::error::Error>> {
    let mut file = File::open(file)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    Ok(contents)
}

pub fn par_reader<P: AsRef<Path> + Debug + Sync + Send>(
    files: Vec<P>,
) -> Result<String, anyhow::Error> {
    let contents: Vec<String> = files
        .par_iter()
        .map(|path| reader(path).unwrap_or_else(|e| panic!("Error reading file: {:?}", e)))
        .collect();

    Ok(contents.concat())
}

fn unpack<P: AsRef<Path> + Debug + Sync + Send>(
    files: Vec<P>,
    cds_overlap: bool,
    is_ref: bool,
) -> Result<GenePredMap, anyhow::Error> {
    let contents = par_reader(files)?;
    let tracks = parse_tracks(&contents, cds_overlap, is_ref)?;

    Ok(tracks)
}

fn parse_tracks<'a>(
    contents: &'a str,
    cds_overlap: bool,
    is_ref: bool,
) -> Result<GenePredMap, anyhow::Error> {
    let pb = get_progress_bar(contents.lines().count() as u64, "Parsing BED12 files");
    let mut tracks = contents
        .par_lines()
        .filter(|x| !x.starts_with("#"))
        .filter_map(|x| Bed12::parse(x, cds_overlap, is_ref).ok())
        .fold(
            || HashMap::new(),
            |mut acc: GenePredMap, record| {
                acc.entry(record.chrom.clone()).or_default().push(record);
                pb.inc(1);
                acc
            },
        )
        .reduce(
            || HashMap::new(),
            |mut acc, map| {
                for (k, v) in map {
                    let acc_v = acc.entry(k).or_insert(Vec::new());
                    acc_v.extend(v);
                }
                acc
            },
        );

    // sort by start/end in descending order
    tracks.par_iter_mut().for_each(|(_, v)| {
        v.par_sort_unstable_by(|a, b| a.start.cmp(&b.start).then(b.end.cmp(&a.end)));
    });

    pb.finish_and_clear();
    info!("Records parsed: {}", tracks.values().flatten().count());

    Ok(tracks)
}

#[inline(always)]
fn exonic_overlap(exons_a: &Vec<(u64, u64)>, exons_b: &Vec<(u64, u64)>) -> bool {
    let mut i = 0;
    let mut j = 0;

    while i < exons_a.len() && j < exons_b.len() {
        let (start_a, end_a) = exons_a[i];
        let (start_b, end_b) = exons_b[j];

        if start_a < end_b && start_b < end_a {
            return true;
        }

        if end_a < end_b {
            i += 1;
        } else {
            j += 1;
        }
    }

    false
}

fn buckerize(
    tracks: GenePredMap,
    overlap_cds: bool,
    colorize: bool,
    amount: usize,
) -> HashMap<String, Vec<Vec<Arc<GenePred>>>> {
    info!("Packing transcripts...");
    let pb = get_progress_bar(amount as u64, "Buckerizing transcripts");
    let cmap = DashMap::new();

    tracks.into_par_iter().for_each(|(chr, records)| {
        let mut acc: Vec<(u64, u64, Vec<Arc<GenePred>>, &str)> = Vec::new();

        for tx in records {
            pb.inc(1);
            let tx_start = if overlap_cds { tx.cds_start } else { tx.start };
            let tx_end = if overlap_cds { tx.cds_end } else { tx.end };

            let mut added = false;

            for (ref mut group_start, ref mut group_end, txs, group_color) in &mut acc {
                if tx_start < *group_end && tx_end > *group_start {
                    if overlap_cds {
                        let exon_overlap = txs
                            .iter()
                            .any(|group_tx| exonic_overlap(&group_tx.exons, &tx.exons));

                        if exon_overlap {
                            *group_start = (*group_start).min(tx_start);
                            *group_end = (*group_end).max(tx_end);
                            let tx_arc = Arc::new(tx.clone());

                            if colorize {
                                txs.push(tx_arc.colorline(*group_color));
                            } else {
                                txs.push(tx_arc);
                            }

                            added = true;
                            break;
                        }
                    } else {
                        *group_start = (*group_start).min(tx_start);
                        *group_end = (*group_end).max(tx_end);

                        let tx_arc = Arc::new(tx.clone());

                        if colorize {
                            txs.push(tx_arc.colorline(*group_color));
                        } else {
                            txs.push(tx_arc);
                        }

                        added = true;
                        break;
                    }
                }
            }

            if !added {
                let color = choose_color();
                let tx_arc = Arc::new(tx);

                if colorize {
                    acc.push((tx_start, tx_end, vec![tx_arc.colorline(color)], color));
                } else {
                    acc.push((tx_start, tx_end, vec![tx_arc], color));
                }
            }
        }

        let acc_map: Vec<Vec<Arc<GenePred>>> = acc.into_iter().map(|(_, _, txs, _)| txs).collect();
        cmap.insert(chr.to_string(), acc_map);
    });

    pb.finish_and_clear();
    info!("Transcripts packed!");
    cmap.into_iter().collect()
}

fn choose_color<'a>() -> &'a str {
    let mut rng = rand::thread_rng();
    let idx = rng.gen_range(0..RGB.len());
    RGB[idx]
}

pub fn packbed<T: AsRef<Path> + Debug + Send + Sync>(
    refs: Vec<T>,
    queries: Vec<T>,
    overlap_cds: bool,
    colorize: bool,
) -> Result<HashMap<String, Vec<Vec<Arc<GenePred>>>>, anyhow::Error> {
    let refs = unpack(refs, overlap_cds, true).unwrap();
    let query = unpack(queries, overlap_cds, false).unwrap();

    let (tracks, n) = combine(refs, query);
    let buckets = buckerize(tracks, overlap_cds, colorize, n);

    Ok(buckets)
}

fn combine(refs: GenePredMap, queries: GenePredMap) -> (GenePredMap, usize) {
    info!("Combining reference and query tracks...");
    let pb = get_progress_bar(queries.values().len() as u64, "Combining tracks");
    let count = queries.values().flatten().count() + refs.values().flatten().count();
    let tracks = Arc::new(Mutex::new(refs));

    queries.into_par_iter().for_each(|(chr, records)| {
        let mut tracks = tracks.lock().expect("Mutex lock failed");
        let acc = tracks.entry(chr).or_default();
        pb.inc(1);

        acc.extend(records);
    });

    let tracks = Arc::try_unwrap(tracks)
        .expect("Arc has more than one reference")
        .into_inner()
        .unwrap();

    pb.finish_and_clear();
    info!("Number of transcripts combined: {}", count);

    (tracks, count)
}
