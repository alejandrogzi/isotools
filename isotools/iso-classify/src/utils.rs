use anyhow::Result;
use bigtools::{utils::reopen::Reopen, BigWigRead};
use dashmap::DashMap;
use hashbrown::{HashMap, HashSet};
use log::{info, warn};
use packbed::{par_reader, reader};
use rayon::prelude::*;
use twobit::TwoBitFile;

use std::path::PathBuf;
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::Arc;
use std::sync::Mutex;

use config::{
    get_progress_bar, BedParser, CoordType, Sequence, SharedSpliceMap, SpliceScores, SpliceSite,
    Strand, StrandSpliceMap, ACCEPTOR_MINUS, ACCEPTOR_PLUS, BGD, CLASSIFY_ASSETS, CONS1, CONS2,
    DONOR_MINUS, DONOR_PLUS, MAXENTSCAN_ACCEPTOR_DB, MAXENTSCAN_DONOR_DB, OVERLAP_CDS, SCALE,
};

pub const MINIMUM_ACCEPTOR_LENGTH: usize = 23;

pub type SpliceScoreMap = HashMap<Sequence, Vec<f64>>;

pub fn make_splice_map<T: AsRef<std::path::Path> + std::fmt::Debug>(
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

pub fn get_splice_scores<T: AsRef<std::path::Path> + std::fmt::Debug>(
    splice_scores: Option<T>,
) -> SpliceScores {
    if let Some(splice_scores) = splice_scores {
        make_splice_map(splice_scores)
    } else {
        log::warn!("No splice scores provided, skipping splice score processing...");
        (vec![DashMap::new()], vec![DashMap::new()])
    }
}

pub fn create_splice_map(
    chr: &str,
    splice_plus: &[StrandSpliceMap],
    splice_minus: &[StrandSpliceMap],
) -> (SharedSpliceMap, SharedSpliceMap) {
    let get_splice_values = |splices: &[DashMap<String, DashMap<usize, f32>>]| {
        (
            splices
                .get(0)
                .and_then(|s| s.get(chr).map(|v| v.value().clone())),
            splices
                .get(1)
                .and_then(|s| s.get(chr).map(|v| v.value().clone())),
        )
    };
    (
        get_splice_values(splice_plus),
        get_splice_values(splice_minus),
    )
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

impl BedParser for Bed4 {
    fn parse(
        line: &str,
        _cds_overlap: bool,
        _is_ref: bool,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        Bed4::new(line.to_string())
    }

    fn chrom(&self) -> &str {
        &self.chrom
    }

    fn coord(&self) -> (u64, u64) {
        self.coord
    }

    // WARN: placeholder for Bed4
    fn intronic_coords(&self) -> HashSet<(u64, u64)> {
        let mut introns = HashSet::new();
        introns.insert(self.coord);
        introns
    }

    // WARN: placeholder for Bed4
    fn exonic_coords(&self) -> HashSet<(u64, u64)> {
        let mut exons = HashSet::new();
        exons.insert(self.coord);
        exons
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

impl BedParser for Bed6 {
    fn parse(
        line: &str,
        _cds_overlap: bool,
        _is_ref: bool,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        Bed6::new(line.to_string())
    }

    fn chrom(&self) -> &str {
        &self.chrom
    }

    fn coord(&self) -> (u64, u64) {
        self.coord
    }

    // WARN: placeholder for Bed6
    fn intronic_coords(&self) -> HashSet<(u64, u64)> {
        let mut introns = HashSet::new();
        introns.insert(self.coord);
        introns
    }

    // WARN: placeholder for Bed6
    fn exonic_coords(&self) -> HashSet<(u64, u64)> {
        let mut exons = HashSet::new();
        exons.insert(self.coord);
        exons
    }
}

pub fn bed_to_map<T>(
    contents: Arc<String>,
    hint: CoordType,
) -> Result<HashMap<String, HashSet<(u64, u64)>>, anyhow::Error>
where
    T: BedParser,
{
    let pb = get_progress_bar(contents.lines().count() as u64, "Parsing BED files...");
    let tracks = contents
        .par_lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| {
            T::parse(line, OVERLAP_CDS, false) // INFO: placeholders
                .map_err(|e| warn!("Error parsing {}: {}", line, e))
                .ok()
        })
        .fold(
            || HashMap::new(),
            |mut acc: HashMap<String, HashSet<(u64, u64)>>, record| {
                let entry = acc.entry(record.chrom().to_owned()).or_default();
                match hint {
                    CoordType::Bounds => {
                        entry.insert(record.coord());
                        ()
                    }
                    CoordType::Intronic => entry.extend(record.intronic_coords()),
                    CoordType::Exonic => entry.extend(record.exonic_coords()),
                };

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
    let tracks = bed_to_map::<Bed4>(contents, CoordType::Bounds).unwrap();

    Some(tracks)
}

pub trait SpliceEntropy: Sized + Send + Sync {
    fn parse(line: &str) -> Result<Self, anyhow::Error>
    where
        Self: Sized;
    fn get_sequence(&self) -> Sequence;
    fn get_scores(self) -> Vec<f64>;
}

pub struct AcceptorScores {
    seq: Sequence,
    scores: Vec<f64>,
}

impl SpliceEntropy for AcceptorScores {
    fn parse(line: &str) -> Result<Self, anyhow::Error> {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 10 {
            return Err(anyhow::anyhow!("ERROR: {} has less than 10 fields!", line));
        }

        let sequence = Sequence::new(parts[0].as_bytes());
        let scores: Vec<f64> = parts[1..]
            .iter()
            .map(|s| s.parse::<f64>().unwrap_or(0.0))
            .collect();

        Ok(AcceptorScores {
            seq: sequence,
            scores,
        })
    }

    fn get_sequence(&self) -> Sequence {
        self.seq.clone()
    }

    fn get_scores(self) -> Vec<f64> {
        self.scores
    }
}

pub struct DonorScores {
    seq: Sequence,
    score: f64,
}

impl SpliceEntropy for DonorScores {
    fn parse(line: &str) -> Result<Self, anyhow::Error> {
        let mut parts = line.split('\t');

        let (sequence, score) = (
            parts
                .next()
                .expect("ERROR: Cannot parse sequence!")
                .as_bytes(),
            parts
                .next()
                .expect("ERROR: Cannot parse donor score!")
                .parse::<f64>()
                .unwrap_or(0.0),
        );

        Ok(DonorScores {
            seq: Sequence::new(sequence),
            score,
        })
    }

    fn get_sequence(&self) -> Sequence {
        self.seq.clone()
    }

    fn get_scores(self) -> Vec<f64> {
        vec![self.score]
    }
}

pub fn parse_tsv<K>(contents: String) -> Result<SpliceScoreMap, anyhow::Error>
where
    K: SpliceEntropy,
{
    let pb = get_progress_bar(contents.lines().count() as u64, "Parsing BED12 files");
    let tracks = contents
        .par_lines()
        .filter(|row| !row.starts_with("#"))
        .filter_map(|row| K::parse(row).ok())
        .fold(
            || HashMap::new(),
            |mut acc: SpliceScoreMap, splice_site| {
                let entry = acc.entry(splice_site.get_sequence()).or_default();
                entry.extend(splice_site.get_scores());

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

    pb.finish_and_clear();
    info!("Records parsed: {}", tracks.values().flatten().count());

    Ok(tracks)
}

pub fn load_scan_scores() -> Option<(SpliceScoreMap, SpliceScoreMap)> {
    let mut assets = std::env::current_dir().expect("Failed to get executable path");

    if !assets.ends_with("iso-classify") {
        let rest = PathBuf::from("iso-classify").join(CLASSIFY_ASSETS);
        assets.push(rest);
    } else {
        assets = assets.join(CLASSIFY_ASSETS);
    }

    let acceptor_scores = parse_tsv::<AcceptorScores>(
        reader(assets.join(MAXENTSCAN_ACCEPTOR_DB)).expect(
            format!(
                "ERROR: Cannot read acceptor scores from {:?}!",
                assets.join(MAXENTSCAN_ACCEPTOR_DB)
            )
            .as_str(),
        ),
    )
    .expect("ERROR: Could not parse acceptor scores!");

    let donor_scores = parse_tsv::<DonorScores>(
        reader(assets.join(MAXENTSCAN_DONOR_DB)).expect("ERROR: Cannot read donor scores!"),
    )
    .expect("ERROR: Could not parse donor scores!");

    Some((donor_scores, acceptor_scores))
}

pub fn get_sequences(twobit: PathBuf) -> Option<DashMap<String, Vec<u8>>> {
    let mut genome = TwoBitFile::open_and_read(twobit).expect("ERROR: Cannot open 2bit file");

    let sequences = DashMap::new();
    genome.chrom_names().iter().for_each(|chr| {
        let seq = genome
            .read_sequence(chr, ..)
            .expect("ERROR: Could not read sequence from .2bit!")
            .as_bytes()
            .to_vec();
        sequences.insert(chr.to_string(), seq);
    });

    Some(sequences)
}

pub fn calculate_acceptor_score(seq: &Sequence, tables: &SpliceScoreMap) -> f64 {
    if seq.len() != MINIMUM_ACCEPTOR_LENGTH {
        let msg = format!(
            "ERROR: Sequence must be a 23-mer for acceptor score calculation, yours is {}!",
            seq.len()
        );
        log::error!("{}", msg);
        std::process::exit(1);
    }

    let me_score = score_max_ent(seq, tables);
    let c_score = score_consensus_seq(seq);

    if me_score == 0.0 {
        return 0.0;
    }

    (c_score * me_score).log2()
}

pub fn score_max_ent(seq: &Sequence, tables: &SpliceScoreMap) -> f64 {
    let seq = seq.skip(18, 20);

    let scores = vec![
        tables.get(&seq.slice(0, 7)).unwrap_or(&vec![0.0])[0],
        tables.get(&seq.slice(7, 14)).unwrap_or(&vec![0.0])[1],
        tables.get(&seq.slice(14, 21)).unwrap_or(&vec![0.0])[2],
        tables.get(&seq.slice(4, 11)).unwrap_or(&vec![0.0])[3],
        tables.get(&seq.slice(11, 18)).unwrap_or(&vec![0.0])[4],
        tables
            .get(&seq.slice_as_seq(4, 7).fill(4))
            .unwrap_or(&vec![0.0])[5],
        tables
            .get(&seq.slice_as_seq(7, 11).fill(3))
            .unwrap_or(&vec![0.0])[6],
        tables
            .get(&seq.slice_as_seq(11, 14).fill(4))
            .unwrap_or(&vec![0.0])[7],
        tables
            .get(&seq.slice_as_seq(14, 18).fill(3))
            .unwrap_or(&vec![0.0])[8],
    ];

    let num: f64 = scores[0..5].iter().product();
    let den: f64 = scores[5..].iter().product();

    let me_score = num / den;

    me_score
}

pub fn score_consensus_seq(seq: &Sequence) -> f64 {
    let nt1 = seq.at_as_bytes(18);
    let nt2 = seq.at_as_bytes(19);

    CONS1[nt1] * CONS2[nt2] / (BGD[nt1] * BGD[nt2])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calculate_consensus_score() {
        let seq = Sequence::new(b"AAAAAAAAAAAAAAAAAAGTAAA");
        let c_score = score_consensus_seq(&seq);

        assert_eq!(c_score, 0.00016425120772946855);

        let seq = Sequence::new(b"TTTAAAAAAAAAAAAAAAGTAAT");
        let c_score = score_consensus_seq(&seq);

        assert_eq!(c_score, 0.00016425120772946854);
    }

    #[test]
    fn test_calculate_max_ent_score() {
        let seq = Sequence::new(b"AAAAAAAAAAAAAAAAAAGTAAA");
        let tables = load_scan_scores().expect("ERROR: Could not load scan scores!");
        let me_score = score_max_ent(&seq, &tables.1);

        assert_eq!(me_score, 0.0003461868180847604);

        let seq = Sequence::new(b"TTTAAAAAAAAAAAAAAAGTAAT");
        let tables = load_scan_scores().expect("ERROR: Could not load scan scores!");
        let me_score = score_max_ent(&seq, &tables.1);

        assert_eq!(me_score, 0.001167787963137549);
    }

    #[test]
    fn test_calculate_acceptor_score() {
        let seq = Sequence::new(b"AAAAAAAAAAAAAAAAAAGTAAA");
        let tables = load_scan_scores().expect("ERROR: Could not load scan scores!");
        let score = calculate_acceptor_score(&seq, &tables.1);

        assert_eq!(score, -24.067969988875006);

        let seq = Sequence::new(b"TTTAAAAAAAAAAAAAAAGTAAT");
        let score = calculate_acceptor_score(&seq, &tables.1);

        assert_eq!(score, -22.313814339694286);
    }

    #[test]
    fn test_calculate_donor_score() {
        let seq = Sequence::new(b"AAGGAAAAA");
        let tables = load_scan_scores().expect("ERROR: Could not load scan scores!");

        let score = tables
            .0
            .get(&seq)
            .expect("ERROR: Could not get donor scores!")
            .get(0)
            .expect("ERROR: Could not get donor scores!");

        assert_eq!(score, &0.192);
    }
}
