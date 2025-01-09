use anyhow::Result;
use bigtools::{utils::reopen::Reopen, BigWigRead};
use dashmap::DashMap;
use hashbrown::{HashMap, HashSet};
use log::{info, warn};
use packbed::par_reader;
use rayon::prelude::*;

use std::path::PathBuf;
use std::sync::{Arc, Mutex};

use std::sync::atomic::{AtomicU32, Ordering};

use config::{
    get_progress_bar, CoordType, SpliceSite, Strand, StrandSpliceMap, ACCEPTOR_MINUS,
    ACCEPTOR_PLUS, DONOR_MINUS, DONOR_PLUS, SCALE,
};

pub trait BedRecord: Send + Sync {
    fn parse(line: String) -> Result<Self, Box<dyn std::error::Error>>
    where
        Self: Sized;
    fn chrom(&self) -> &str;
    fn coord(&self) -> (u64, u64);
    fn intronic_coords(&self) -> HashSet<(u64, u64)>;
    fn exonic_coords(&self) -> HashSet<(u64, u64)>;
}

#[derive(Debug, PartialEq, Clone)]
pub struct Bed4 {
    pub chrom: String,
    pub coord: (u64, u64),
    pub id: String,
}

impl Bed4 {
    pub fn new(line: String) -> Result<Bed4, Box<dyn std::error::Error>> {
        if line.is_empty() {
            return Err("Empty line".into());
        }

        let mut fields = line.split('\t');
        let get = |field: &str| field.parse::<u64>().map_err(|_| "Cannot parse field");

        let (chrom, start, end, id) = (
            fields.next().ok_or("Cannot parse chrom")?.to_string(),
            get(fields.next().ok_or("Cannot parse start")?)?,
            get(fields.next().ok_or("Cannot parse end")?)?,
            fields.next().ok_or("Cannot parse id")?.to_string(),
        );

        Ok(Bed4 {
            chrom,
            coord: (start + 1, end - 1), // 0-based to 1-based
            id,
        })
    }

    pub fn from(chrom: String, start: u64, end: u64, id: String) -> Bed4 {
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

impl BedRecord for Bed4 {
    fn parse(line: String) -> Result<Self, Box<dyn std::error::Error>> {
        Bed4::new(line)
    }

    fn chrom(&self) -> &str {
        &self.chrom
    }

    fn coord(&self) -> (u64, u64) {
        self.coord
    }
}

// WARN: will cover any high-order bed file [6,8,12]
#[derive(Debug, PartialEq, Clone)]
pub struct Bed6 {
    pub chrom: String,
    pub coord: (u64, u64),
    pub id: String,
    pub strand: Strand,
}

impl Bed6 {
    pub fn new(line: String) -> Result<Bed6, Box<dyn std::error::Error>> {
        if line.is_empty() {
            return Err("ERROR: Empty line in .bed!".into());
        }

        let mut fields = line.split('\t');
        let get = |field: &str| field.parse::<u64>().map_err(|_| "Cannot parse field");

        let (chrom, start, end, id, _, strand) = (
            fields.next().ok_or("Cannot parse chrom")?.to_string(),
            get(fields.next().ok_or("Cannot parse start")?)?,
            get(fields.next().ok_or("Cannot parse end")?)?,
            fields.next().ok_or("Cannot parse id")?.to_string(),
            fields.next().ok_or("Cannot parse score")?,
            fields
                .next()
                .ok_or("ERROR: Cannot parse strand!")?
                .parse::<Strand>()?,
        );

        let (start, end) = match strand {
            Strand::Forward => {
                if start > end {
                    return Err("ERROR: Start is greater than end!".into());
                }
                (start, end)
            }
            Strand::Reverse => {
                if start < end {
                    return Err("ERROR: Start is less than end!".into());
                }
                (SCALE - end, SCALE - start)
            }
        };

        Ok(Bed6 {
            chrom,
            coord: (start, end),
            id,
            strand,
        })
    }

    pub fn from(chrom: String, start: u64, end: u64, id: String, strand: Strand) -> Bed6 {
        Bed6 {
            chrom,
            coord: (start, end),
            id,
            strand,
        }
    }

    pub fn send(&self, acc: &mut String) {
        acc.push_str(&format!(
            "{}\t{}\t{}\n",
            self.chrom, self.coord.0, self.coord.1
        ));
    }
}

impl BedRecord for Bed6 {
    fn parse(line: String) -> Result<Self, Box<dyn std::error::Error>> {
        Bed6::new(line)
    }

    fn chrom(&self) -> &str {
        &self.chrom
    }

    fn coord(&self) -> (u64, u64) {
        self.coord
    }
}

pub fn bed_to_map<T>(
    contents: Arc<String>,
    hint: CoordType,
) -> Result<HashMap<String, HashSet<(u64, u64)>>, anyhow::Error>
where
    T: BedRecord,
{
    let pb = get_progress_bar(contents.lines().count() as u64, "Parsing BED files...");
    let tracks = contents
        .par_lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| {
            T::parse(line.to_string())
                .map_err(|e| warn!("Error parsing {}: {}", line, e))
                .ok()
        })
        .fold(
            || HashMap::new(),
            |mut acc: HashMap<String, HashSet<(u64, u64)>>, record| {
                acc.entry(record.chrom().to_owned())
                    .or_default()
                    .insert(match hint {
                        CoordType::Bounds => record.coord(),
                        CoordType::Intronic => record.intronic_coords(),
                        CoordType::Exonic => record.exonic_coords(),
                    });
                pb.inc(1);
                acc
            },
        )
        .reduce(
            || HashMap::new(),
            |mut acc, map| {
                for (k, v) in map {
                    acc.entry(k).or_default().extend(v);
                }
                acc
            },
        );

    pb.finish_and_clear();
    if tracks.is_empty() {
        anyhow::bail!("No tracks found in the provided file!");
    }
    info!(
        "Parsed {} blacklisted intervals.",
        tracks.values().flatten().count()
    );
    Ok(tracks)
}

pub fn unpack_blacklist<'a>(paths: Vec<PathBuf>) -> Option<HashMap<String, HashSet<(u64, u64)>>> {
    if paths.is_empty() {
        return None;
    }

    let contents = Arc::new(par_reader(paths).unwrap());
    let tracks = bed_to_map::<Bed4>(contents).unwrap();

    Some(tracks)
}

pub fn get_toga_coords(toga: Option<PathBuf>) -> Option<HashMap<String, HashSet<(u64, u64)>>> {
    if let Some(toga) = toga {
        let contents =
            Arc::new(par_reader(vec![toga]).expect("ERROR: Cannot read TOGA annotation file!"));
        let tracks = bed_to_map::<Bed6>(contents).expect("ERROR: Cannot parse TOGA annotation!");

        Some(tracks)
    } else {
        return None;
    }
}

pub fn get_splice_scores<T: AsRef<std::path::Path> + std::fmt::Debug>(
    dir: T,
) -> (Vec<StrandSpliceMap>, Vec<StrandSpliceMap>) {
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
) -> Vec<DashMap<String, DashMap<usize, f32>>> {
    let total_count = AtomicU32::new(0);
    let rs = Mutex::new(vec![DashMap::new(), DashMap::new()]);

    // [donor, acceptor]
    bigwigs
        .into_par_iter()
        .zip(vec![SpliceSite::Donor, SpliceSite::Acceptor])
        .for_each(|(bigwig, site)| {
            let acc = DashMap::new();

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

                let mapper = DashMap::new();
                let local_count = AtomicU32::new(0);

                values.into_iter().enumerate().for_each(|(i, v)| {
                    if v >= config::SPLICE_AI_SCORE_RECOVERY_THRESHOLD {
                        let pos = i;
                        mapper.entry(pos).or_insert(v);
                        local_count.fetch_add(1, Ordering::Relaxed);
                    }
                });

                acc.insert(name, mapper);
                total_count.fetch_add(local_count.load(Ordering::Relaxed), Ordering::Relaxed);
            });

            let mut guard = rs.lock().expect("ERROR: Cannot lock mutex");
            match site {
                SpliceSite::Donor => guard[0] = acc,
                SpliceSite::Acceptor => guard[1] = acc,
            }
        });

    info!(
        "Parsed and combined {} significant splicing scores from BigWigs!",
        total_count.load(Ordering::Relaxed)
    );

    rs.into_inner()
        .expect("ERROR: Cannot unwrap collection of SpliceAI scores!")
}
