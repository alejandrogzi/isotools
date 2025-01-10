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
pub use record::{Bed12, GenePred, IntronPred, RefGenePred};

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

pub trait BedPackage: Send + Sync + Debug {}

pub type DefaultBucket = (RefGenePred, Vec<GenePred>);
impl BedPackage for DefaultBucket {}
impl BedPackage for IntronPred {}

#[derive(Debug, Clone)]
pub enum PackMode {
    Default,
    Intron,
    Exon,
}

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
        .map(|path| reader(path).unwrap_or_else(|e| panic!("ERROR: Could not read file: {:?}", e)))
        .collect();

    Ok(contents.concat())
}

pub fn unpack<P: AsRef<Path> + Debug + Sync + Send>(
    files: Vec<P>,
    cds_overlap: bool,
    is_ref: bool,
) -> Result<GenePredMap, anyhow::Error> {
    let contents = par_reader(files)?;
    let tracks = parse_tracks(&contents, cds_overlap, is_ref)?;

    Ok(tracks)
}

pub fn parse_tracks<'a>(
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

#[derive(Debug, Clone)]
struct UnionFind {
    parent: Vec<usize>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
        }
    }

    #[inline(always)]
    fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            self.parent[x] = self.find(self.parent[x]);
        }
        self.parent[x]
    }

    #[inline(always)]
    fn union(&mut self, x: usize, y: usize) {
        let root_x = self.find(x);
        let root_y = self.find(y);
        if root_x != root_y {
            self.parent[root_y] = root_x;
        }
    }
}

// Vec<(RefGenePred, Vec<GenePred>)>

pub fn buckerize(
    tracks: GenePredMap,
    overlap_cds: bool,
    overlap_exon: bool,
    amount: usize,
    mode: PackMode,
) -> DashMap<String, Vec<Box<dyn BedPackage>>> {
    let cmap = DashMap::new();
    let pb = get_progress_bar(amount as u64, "Buckerizing transcripts");

    tracks.into_par_iter().for_each(|(chr, transcripts)| {
        let mut exons = Vec::new();
        let mut id_map = HashMap::new();
        let mut uf = UnionFind::new(transcripts.len());

        // if base mode, tx boundaries will behave as exons ranges
        for (i, transcript) in transcripts.iter().enumerate() {
            id_map.insert(i, transcript);

            if !overlap_exon && !overlap_cds {
                exons.push((transcript.start, transcript.end, i));
            } else {
                for &(start, end) in &transcript.exons {
                    exons.push((start, end, i));
                }
            }
        }

        exons.par_sort_unstable_by(|a, b| a.0.cmp(&b.0));

        let mut prev_end = exons[0].1;
        let mut prev_idx = exons[0].2;
        for &(start, end, idx) in &exons[1..] {
            if start < prev_end {
                uf.union(prev_idx, idx);
                prev_end = prev_end.max(end);
            } else {
                // no overlap, update prev_end and prev_idx
                prev_end = end;
                prev_idx = idx;
            }
        }

        let mut groups = HashMap::new();
        for i in 0..transcripts.len() {
            pb.inc(1);
            let root = uf.find(i);
            groups
                .entry(root)
                .or_insert_with(Vec::new)
                .push(id_map[&i].clone());
        }

        let comps = groups
            .into_iter()
            .map(|(_, group)| {
                // separate refs from queries, use refs to build RefGenePred
                let (refs, queries): (Vec<_>, Vec<_>) = group.into_iter().partition(|x| x.is_ref);

                match mode {
                    PackMode::Default => {
                        let refs = RefGenePred::from(refs);
                        return Box::new((refs, queries)) as Box<dyn BedPackage>;
                    }
                    PackMode::Intron => {
                        // queries are TOGA introns
                        let refs = IntronPred::from(refs, queries);
                        return Box::new(refs) as Box<dyn BedPackage>;
                    }
                    PackMode::Exon => {
                        // let refs = ExonPred::from(refs);
                        // return refs;
                        todo!()
                    }
                }
            })
            .collect::<Vec<_>>();

        cmap.insert(chr, comps);
    });

    pb.finish_and_clear();
    info!("Transcripts packed!");
    cmap
}

#[allow(dead_code)]
fn choose_color<'a>() -> &'a str {
    let mut rng = rand::thread_rng();
    let idx = rng.gen_range(0..RGB.len());
    RGB[idx]
}

pub fn packbed<T: AsRef<Path> + Debug + Send + Sync>(
    refs: Vec<T>,
    queries: Option<Vec<T>>,
    overlap_cds: bool,
    overlap_exon: bool,
    mode: PackMode,
) -> Result<DashMap<String, Vec<Box<dyn BedPackage>>>, anyhow::Error> {
    let refs = unpack(refs, overlap_cds, true).expect("ERROR: Failed to unpack reference tracks");

    let (tracks, n) = if let Some(query) = queries {
        let query =
            unpack(query, overlap_cds, false).expect("ERROR: Failed to unpack query tracks");
        combine(refs, query)
    } else {
        let n = refs.values().flatten().count();
        (refs, n)
    };

    let buckets = buckerize(tracks, overlap_cds, overlap_exon, n, mode);
    Ok(buckets)
}

pub fn combine(refs: GenePredMap, queries: GenePredMap) -> (GenePredMap, usize) {
    info!("Combining reference and query tracks...");
    let pb = get_progress_bar(queries.values().len() as u64, "Combining tracks");
    let count = queries.values().flatten().count() + refs.values().flatten().count();
    let tracks = Arc::new(Mutex::new(refs));

    queries.into_par_iter().for_each(|(chr, records)| {
        let mut tracks = tracks.lock().expect("ERROR: Mutex lock failed");
        let acc = tracks.entry(chr).or_default();
        pb.inc(1);

        acc.extend(records);
    });

    let tracks = Arc::try_unwrap(tracks)
        .expect("ERROR: Arc has more than one reference")
        .into_inner()
        .unwrap();

    pb.finish_and_clear();
    info!("Number of transcripts combined: {}", count);

    (tracks, count)
}
