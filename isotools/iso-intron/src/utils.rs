use anyhow::Result;
use config::get_progress_bar;
use hashbrown::{HashMap, HashSet};
use log::{info, warn};
use packbed::par_reader;
use rayon::prelude::*;
use std::path::PathBuf;

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
            coord: (start, end),
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
