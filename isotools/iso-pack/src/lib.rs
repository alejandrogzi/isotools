use std::any::Any;
use std::fmt::Debug;
use std::fs::File;
use std::io::Read;
use std::path::Path;
use std::sync::{Arc, Mutex};

use config::{get_progress_bar, BedParser, OverlapType};
use dashmap::DashMap;
use hashbrown::HashMap;
use log::info;
use rand::Rng;
use rayon::prelude::*;

pub mod record;
pub use record::{Bed12, GenePred, IntronBucket, IntronPred, PolyAPred, RefGenePred};

pub type GenePredMap = HashMap<String, Box<dyn BedPackage>>;

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

pub trait BedPackage: Send + Sync + Debug + Any {
    fn as_any(&self) -> &dyn Any;
    fn as_any_mut(&mut self) -> &mut dyn Any;
    fn as_any_owned(self: Box<Self>) -> Box<dyn Any>;
}

pub type DefaultBucket = (RefGenePred, Vec<GenePred>);

impl BedPackage for DefaultBucket {
    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }

    fn as_any_owned(self: Box<Self>) -> Box<dyn Any> {
        self
    }
}

impl BedPackage for IntronBucket {
    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }

    fn as_any_owned(self: Box<Self>) -> Box<dyn Any> {
        self
    }
}

impl BedPackage for Vec<GenePred> {
    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }

    fn as_any_owned(self: Box<Self>) -> Box<dyn Any> {
        self
    }
}

impl BedPackage for Vec<IntronPred> {
    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }

    fn as_any_owned(self: Box<Self>) -> Box<dyn Any> {
        self
    }
}

impl BedPackage for (Vec<IntronPred>, Vec<GenePred>) {
    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }

    fn as_any_owned(self: Box<Self>) -> Box<dyn Any> {
        self
    }
}

impl BedPackage for (Vec<GenePred>, Vec<GenePred>) {
    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }

    fn as_any_owned(self: Box<Self>) -> Box<dyn Any> {
        self
    }
}

impl BedPackage for Vec<PolyAPred> {
    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }

    fn as_any_owned(self: Box<Self>) -> Box<dyn Any> {
        self
    }
}

impl BedPackage for (Vec<PolyAPred>, Vec<GenePred>) {
    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }

    fn as_any_owned(self: Box<Self>) -> Box<dyn Any> {
        self
    }
}

#[derive(Debug, Clone)]
pub enum PackMode {
    Default, // RefGenePred + Vec<GenePred>
    Intron,  // IntronPred
    Exon,    // ExonPred
    Query,   // Vec<GenePred> + Vec<IntronPred>
    Paired,  // Vec<GenePred> + Vec<GenePred>
    PolyA,   // Vec<PolyAPred>
}

pub fn reader<P: AsRef<Path> + Debug>(file: P) -> Result<String, Box<dyn std::error::Error>> {
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
        .map(|path| {
            reader(path)
                .unwrap_or_else(|e| panic!("ERROR: Could not read file: {:?} -> {:?}", e, path))
        })
        .collect();

    Ok(contents.concat())
}

pub fn unpack<K, P>(
    files: Vec<P>,
    overlap: OverlapType,
    is_ref: bool,
) -> Result<HashMap<String, Vec<K>>, anyhow::Error>
where
    K: BedParser + Debug + Send + Sync,
    P: AsRef<Path> + Debug + Sync + Send,
{
    let contents = par_reader(files)?;
    let tracks = parse_tracks::<K>(&contents, overlap, is_ref)?;

    Ok(tracks)
}

pub fn parse_tracks<'a, K>(
    contents: &'a str,
    overlap: OverlapType,
    is_ref: bool,
) -> Result<HashMap<String, Vec<K>>, anyhow::Error>
where
    K: BedParser + Debug + Send + Sync,
{
    let pb = get_progress_bar(contents.lines().count() as u64, "Parsing BED12 files");
    let mut tracks = contents
        .par_lines()
        .filter(|row| !row.starts_with("#"))
        .filter_map(|row| K::parse(row, overlap, is_ref).ok())
        .fold(
            || HashMap::new(),
            |mut acc: HashMap<String, Vec<K>>, record| {
                acc.entry(record.chrom().to_string())
                    .or_insert_with(Vec::new)
                    .push(record);

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

    // INFO: sort by start/end in descending order
    tracks.par_iter_mut().for_each(|(_, v)| {
        // v.par_sort_unstable_by(|a, b| a.start.cmp(&b.start).then(b.end.cmp(&a.end)));
        v.par_sort_unstable_by(|a, b| {
            let (a_start, a_end) = a.coord();
            let (b_start, b_end) = b.coord();

            a_start.cmp(&b_start).then(b_end.cmp(&a_end))
        });
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

pub fn buckerize(
    tracks: HashMap<String, Vec<GenePred>>,
    overlap: OverlapType,
    amount: usize,
    mode: PackMode,
) -> DashMap<String, Vec<Box<dyn BedPackage>>> {
    let cmap = DashMap::new();
    let pb = get_progress_bar(amount as u64, "Buckerizing transcripts");

    tracks.into_par_iter().for_each(|(chr, transcripts)| {
        let mut exons = Vec::new();
        let mut id_map = HashMap::new();
        let mut uf = UnionFind::new(transcripts.len());

        for (i, transcript) in transcripts.iter().enumerate() {
            id_map.insert(i, transcript);

            match overlap {
                OverlapType::Boundary => {
                    // INFO: if no overlap choice, use tx boundaries as exons
                    exons.push((transcript.start, transcript.end, i));
                }
                OverlapType::CDSBound => {
                    for &(start, end) in &transcript.exons {
                        // INFO: UTRs are not considered in the overlap
                        if end < transcript.cds_start || start > transcript.cds_end {
                            continue;
                        }

                        // INFO: taking care of 5' UTR-CDS nested exons
                        if start < transcript.cds_start && end < transcript.cds_end {
                            exons.push((transcript.cds_start, end, i));
                            continue;
                        }

                        // INFO: taking care of CDS-3' UTR nested exons
                        if start > transcript.cds_start && end > transcript.cds_end {
                            exons.push((start, transcript.cds_end, i));
                            continue;
                        }

                        exons.push((start, end, i));
                    }
                }
                OverlapType::Exon | OverlapType::CDS => {
                    for &(start, end) in &transcript.exons {
                        exons.push((start, end, i));
                    }
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
            .filter_map(|(_, group)| {
                // INFO: separate refs from queries and match mode
                let (refs, queries): (Vec<_>, Vec<_>) = group.into_iter().partition(|x| x.is_ref);

                match mode {
                    PackMode::Default => {
                        // INFO: separate refs from queries, use refs to build RefGenePred
                        // INFO: used in iso-utr
                        let refs = RefGenePred::from(refs);
                        return Some(Box::new((refs, queries)) as Box<dyn BedPackage>);
                    }
                    PackMode::Intron => {
                        // WARN: queries are TOGA introns -> avoid empty refs!
                        // INFO: used in iso-classify
                        if refs.is_empty() {
                            return None;
                        }

                        let refs = IntronBucket::from(refs, queries);
                        return Some(Box::new(refs) as Box<dyn BedPackage>);
                    }
                    PackMode::Query => {
                        // INFO: introns are refs, reads are queries [both Vec<GenePred>]
                        // INFO: used in iso-intron
                        // WARN: avoid empty queries or refs -> should be at least one of each
                        // because we are using one to build the other!
                        if queries.is_empty() | refs.is_empty() {
                            return None;
                        }

                        // INFO: Vec<GenePred> -> Vec<IntronPred>
                        let introns = refs
                            .into_iter()
                            .map(|read| IntronPred::from(read))
                            .collect::<Vec<_>>();

                        return Some(Box::new((introns, queries)) as Box<dyn BedPackage>);
                    }
                    PackMode::Exon => {
                        // let refs = ExonPred::from(refs);
                        // return refs;
                        todo!()
                    }
                    PackMode::Paired => {
                        // INFO: both refs and queries are Vec<GenePred>
                        // INFO: used in iso-orf
                        return Some(Box::new((refs, queries)) as Box<dyn BedPackage>);
                    }
                    PackMode::PolyA => {
                        // INFO: queries are an optional TOGA projection file to determine polyA pos
                        // INFO: used in pas-caller [iso-polya caller]
                        if refs.is_empty() {
                            return None;
                        }

                        // INFO: Vec<GenePred> -> Vec<PolyAPred>
                        let reads = refs
                            .into_iter()
                            .map(|read| PolyAPred::from(read))
                            .collect::<Vec<_>>();

                        return Some(Box::new((reads, queries)) as Box<dyn BedPackage>);
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
    overlap: OverlapType,
    mode: PackMode,
) -> Result<DashMap<String, Vec<Box<dyn BedPackage>>>, anyhow::Error> {
    let (tracks, n) = match mode {
        PackMode::Default
        | PackMode::Intron
        | PackMode::Exon
        | PackMode::Paired
        | PackMode::PolyA => {
            let refs = unpack::<GenePred, _>(refs, overlap, true)
                .expect("ERROR: Failed to unpack reference tracks");

            if let Some(query) = queries {
                let query =
                    unpack(query, overlap, false).expect("ERROR: Failed to unpack query tracks");
                combine::<GenePred, GenePred>(refs, query)
            } else {
                let n = refs.values().flatten().count();
                (refs, n)
            }
        }
        PackMode::Query => {
            let refs = unpack::<IntronPred, _>(refs, overlap, true)
                .expect("ERROR: Failed to unpack reference tracks");
            let query = unpack::<GenePred, _>(
                queries.expect("ERROR: queries were not provided!"),
                overlap,
                false,
            )
            .expect("ERROR: Failed to unpack query tracks");

            combine::<IntronPred, GenePred>(refs, query)
        }
    };

    let buckets = buckerize(tracks, overlap, n, mode);
    Ok(buckets)
}

pub fn combine<P, K>(
    refs: HashMap<String, Vec<P>>,
    queries: HashMap<String, Vec<K>>,
) -> (HashMap<String, Vec<K>>, usize)
where
    K: BedParser + Debug + Send + Sync + 'static,
    P: BedParser + Debug + Send + Sync + Into<K> + 'static,
{
    info!("Combining reference and query tracks...");
    let pb = get_progress_bar(queries.values().len() as u64, "Combining tracks");
    let count = queries.values().flatten().count() + refs.values().flatten().count();
    let tracks = Arc::new(Mutex::new(queries));

    refs.into_par_iter().for_each(|(chr, records)| {
        let mut tracks = tracks.lock().expect("ERROR: Mutex lock failed");
        let acc = tracks.entry(chr).or_default();
        pb.inc(1);

        // INFO: checking if K and P are the same type
        if std::any::TypeId::of::<K>() == std::any::TypeId::of::<P>() {
            // INFO: K and P are the same type, so this is safe
            let records: Vec<K> = unsafe { std::mem::transmute(records) };
            acc.extend(records);
        } else {
            acc.extend(records.into_iter().map(|p| p.into()));
        }
    });

    let tracks = Arc::try_unwrap(tracks)
        .expect("ERROR: Arc has more than one reference")
        .into_inner()
        .unwrap();

    pb.finish_and_clear();
    info!("Number of transcripts combined: {}", count);

    (tracks, count)
}
