use anyhow::Result;
use bigtools::{utils::reopen::Reopen, BigWigRead};
use dashmap::DashMap;
use hashbrown::{HashMap, HashSet};
use log::{info, warn};
use packbed::par_reader;
use rayon::prelude::*;
use std::path::PathBuf;

use std::sync::atomic::{AtomicU32, Ordering};

use config::{get_progress_bar, SpliceMap, ACCEPTOR_MINUS, ACCEPTOR_PLUS, DONOR_MINUS, DONOR_PLUS};

#[derive(Debug, PartialEq, Clone)]
pub struct Bed4<'a> {
    pub chrom: &'a str,
    pub coord: (u64, u64),
    pub id: &'a str,
}

impl<'a> Bed4<'a> {
    pub fn new(line: &str) -> Result<Bed4, &'static str> {
        if line.is_empty() {
            return Err("Empty line");
        }

        let mut fields = line.split('\t');
        let get = |field: &str| field.parse::<u64>().map_err(|_| "Cannot parse field");

        let (chrom, start, end, id) = (
            fields.next().ok_or("Cannot parse chrom")?,
            get(fields.next().ok_or("Cannot parse start")?)?,
            get(fields.next().ok_or("Cannot parse end")?)?,
            fields.next().ok_or("Cannot parse id")?,
        );

        Ok(Bed4 {
            chrom,
            coord: (start + 1, end - 1),
            id,
        })
    }

    pub fn from(chrom: &'a str, start: u64, end: u64, id: &'a str) -> Bed4<'a> {
        Bed4 {
            chrom,
            coord: (start, end),
            id,
        }
    }

    pub fn send(&self, acc: &mut String) {
        acc.push_str(&format!(
            "{}\t{}\t{}\n",
            self.chrom, self.coord.0, self.coord.1
        ));
    }
}

pub fn unpack_blacklist<'a>(paths: Vec<PathBuf>) -> Option<HashMap<String, HashSet<(u64, u64)>>> {
    if paths.is_empty() {
        return None;
    }

    let contents = par_reader(paths).unwrap();
    let tracks = parse_bed4(&contents).unwrap();

    Some(tracks)
}

pub fn parse_bed4<'a>(
    contents: &'a str,
) -> Result<HashMap<String, HashSet<(u64, u64)>>, anyhow::Error> {
    let pb = get_progress_bar(contents.lines().count() as u64, "Parsing BED4 files...");
    let tracks = contents
        .par_lines()
        .filter(|x| !x.starts_with("#"))
        .filter_map(|x| Bed4::new(x).map_err(|e| warn!("{} from: {}. ", x, e)).ok())
        .fold(
            || HashMap::new(),
            |mut acc: HashMap<String, HashSet<(u64, u64)>>, record| {
                acc.entry(record.chrom.to_string())
                    .or_default()
                    .insert(record.coord);
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
    match tracks.is_empty() {
        true => {
            anyhow::bail!("Blacklist file provided but no tracks found!")
        }
        false => {
            info!(
                "Parsed {} blacklisted introns!",
                tracks.values().flatten().count()
            );

            Ok(tracks)
        }
    }
}

pub fn get_splice_scores<T: AsRef<std::path::Path> + std::fmt::Debug>(dir: T) -> SpliceMap {
    let plus = vec![
        dir.as_ref().join(DONOR_PLUS),
        dir.as_ref().join(ACCEPTOR_PLUS),
    ];
    let minus = vec![
        dir.as_ref().join(DONOR_MINUS),
        dir.as_ref().join(ACCEPTOR_MINUS),
    ];

    info!("Parsing BigWigs...");
    let (plus, minus) = rayon::join(|| bigwig_to_map(plus), || bigwig_to_map(minus));

    (plus, minus)
}

fn bigwig_to_map<T: AsRef<std::path::Path> + std::fmt::Debug + Sized + Sync>(
    bigwigs: Vec<T>,
) -> DashMap<String, DashMap<usize, f32>> {
    let acc = DashMap::new();
    let total_count = AtomicU32::new(0);

    bigwigs
        .into_par_iter()
        .enumerate()
        .for_each(|(hint, bigwig)| {
            let bwread = BigWigRead::open_file(&bigwig).expect("ERROR: Cannot open BigWig file");
            let chroms: Vec<_> = bwread.chroms().to_vec();

            chroms.into_par_iter().for_each(|chr| {
                let mut bwread =
                    BigWigRead::reopen(&bwread).expect("ERROR: Cannot re-open BigWig file");

                let name = chr.name.clone();
                let length = chr.length;
                let values = bwread
                    .values(&name, 0, length)
                    .expect("ERROR: Cannot read values from BigWig!");

                let local_mapper = DashMap::new();
                let local_count = AtomicU32::new(0);

                values.into_iter().enumerate().for_each(|(i, v)| {
                    if v > 0.0 {
                        // donor [+1 to match isotools/bigtools coords]
                        let key = if hint == 0 { i + 1 } else { i };
                        local_mapper.entry(key).or_insert(v);
                        local_count.fetch_add(1, Ordering::Relaxed);
                    }
                });

                // merge local_mapper into the global acc
                let global_entry = acc.entry(name).or_insert_with(DashMap::new);
                local_mapper.into_iter().for_each(|(k, v)| {
                    global_entry.entry(k).or_insert(v);
                });

                total_count.fetch_add(local_count.load(Ordering::Relaxed), Ordering::Relaxed);
            });
        });

    info!(
        "Parsed and combined {} scores from BigWigs!",
        total_count.load(Ordering::Relaxed)
    );

    acc
}
