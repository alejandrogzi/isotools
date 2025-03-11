use config::{BedParser, OverlapType, Strand};
use dashmap::DashMap;
use hashbrown::HashSet;
use packbed::{record::abs_pos, unpack};
use serde::{Deserialize, Serialize};
use twobit::TwoBitFile;

use std::collections::HashMap;
use std::error::Error;
use std::path::PathBuf;
use std::str::FromStr;

pub const EXPANSION_SIZE: u64 = 150; // 150bp
pub const ISO_POLYA: &str = "iso-polya";
pub const ASSETS: &str = "assets";

#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct MiniPolyAPred {
    pub name: String,
    pub chrom: String,
    pub strand: Strand,
    pub start: u64,
    pub end: u64,
}

impl MiniPolyAPred {
    #[inline(always)]
    pub fn read(line: &str, _: OverlapType, _: bool) -> Result<Self, &'static str> {
        if line.is_empty() {
            return Err("Empty line");
        }

        let mut fields = line.split('\t');
        let (chrom, tx_start, tx_end, name, _, strand) = (
            fields.next().ok_or("Cannot parse chrom")?,
            fields.next().ok_or("Cannot parse tx_start")?,
            fields.next().ok_or("Cannot parse tx_end")?,
            fields.next().ok_or("Cannot parse name")?,
            fields.next().ok_or("Cannot parse score")?,
            fields
                .next()
                .ok_or("Cannot parse strand")?
                .chars()
                .next()
                .ok_or("Cannot parse strand as char")?,
        );

        let get = |field: &str| {
            field
                .parse::<u64>()
                .map_err(|_| "ERROR: Cannot parse field")
        };

        let strand = match strand {
            '+' => Strand::Forward,
            '-' => Strand::Reverse,
            _ => return Err("ERROR: Strand is not + or -"),
        };

        // WARN: not doing any coord conversion!
        let (start, end) = match strand {
            Strand::Forward => (get(tx_end)? - EXPANSION_SIZE, get(tx_end)? + EXPANSION_SIZE),
            Strand::Reverse => (
                get(tx_start)? - EXPANSION_SIZE,
                get(tx_start)? + EXPANSION_SIZE,
            ),
        };

        Ok(MiniPolyAPred {
            name: name.into(),
            chrom: chrom.into(),
            strand,
            start,
            end,
        })
    }
}

impl BedParser for MiniPolyAPred {
    fn parse(
        line: &str,
        overlap: OverlapType,
        is_ref: bool,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let data = MiniPolyAPred::read(line, overlap, is_ref).expect("ERROR: Cannot parse line");
        Ok(data)
    }

    fn chrom(&self) -> &str {
        &self.chrom.as_str()
    }

    fn coord(&self) -> (u64, u64) {
        (self.start, self.end)
    }

    fn name(&self) -> &str {
        &self.name.as_str()
    }

    fn strand(&self) -> Strand {
        self.strand.clone()
    }

    fn start(&self) -> u64 {
        self.start
    }

    fn end(&self) -> u64 {
        self.end
    }

    fn cds_start(&self) -> u64 {
        self.start
    }

    fn cds_end(&self) -> u64 {
        self.end
    }

    // WARN: placeholder for trait
    fn intronic_coords(&self) -> HashSet<(u64, u64)> {
        HashSet::new()
    }
    fn exonic_coords(&self) -> HashSet<(u64, u64)> {
        HashSet::new()
    }
    fn score(&self) -> f32 {
        0.0
    }
    fn block_sizes(&self) -> Vec<u64> {
        vec![]
    }
    fn block_starts(&self) -> Vec<u64> {
        vec![]
    }
    fn block_count(&self) -> u64 {
        0
    }
    fn rgb(&self) -> &str {
        ""
    }
}

pub fn get_sequences<'a>(
    twobit: PathBuf,
) -> Option<(DashMap<String, Vec<u8>>, HashMap<String, u32>)> {
    let mut genome = TwoBitFile::open_and_read(twobit).expect("ERROR: Cannot open 2bit file");

    let chrom_sizes = genome
        .chrom_names()
        .into_iter()
        .zip(genome.chrom_sizes().into_iter())
        .map(|(chr, size)| (chr, size as u32))
        .collect();

    let sequences = DashMap::new();
    genome.chrom_names().iter().for_each(|chr| {
        let seq = genome
            .read_sequence(chr, ..)
            .expect("ERROR: Could not read sequence from .2bit!")
            .as_bytes()
            .to_vec();
        sequences.insert(chr.to_string(), seq);
    });

    Some((sequences, chrom_sizes))
}

#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct BedGraph {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub score: f64,
    pub strand: Strand,
}

impl BedGraph {
    pub fn read(line: &str) -> Result<BedGraph, &'static str> {
        if line.is_empty() {
            return Err("Empty line");
        }

        let mut fields = line.split('\t');
        let (chrom, start, end, score, strand) = (
            fields.next().ok_or("ERROR: Cannot parse chrom")?,
            fields.next().ok_or("ERROR: Cannot parse start")?,
            fields.next().ok_or("ERROR: Cannot parse end")?,
            fields.next().ok_or("ERROR: Cannot parse score")?,
            fields.next().ok_or("ERROR: Cannot parse strand")?,
        );

        let get = |field: &str| {
            field
                .parse::<u64>()
                .map_err(|_| "ERROR: Cannot parse field")
        };

        let get_score = |field: &str| {
            field
                .parse::<f64>()
                .map_err(|_| "ERROR: Cannot parse score")
        };

        Ok(BedGraph {
            chrom: chrom.into(),
            start: get(start)?,
            end: get(end)?,
            score: get_score(score)?,
            strand: Strand::from_str(strand).expect("ERROR: Cannot parse strand"),
        })
    }
}

impl BedParser for BedGraph {
    fn parse(line: &str, _: OverlapType, _: bool) -> Result<Self, Box<dyn std::error::Error>> {
        let data = BedGraph::read(line).expect("ERROR: Cannot parse line");
        Ok(data)
    }

    fn chrom(&self) -> &str {
        &self.chrom.as_str()
    }

    fn coord(&self) -> (u64, u64) {
        (self.start, self.end)
    }

    fn strand(&self) -> Strand {
        self.strand.clone()
    }

    fn start(&self) -> u64 {
        self.start
    }

    fn end(&self) -> u64 {
        self.end
    }

    fn cds_start(&self) -> u64 {
        self.start
    }

    fn cds_end(&self) -> u64 {
        self.end
    }

    fn score(&self) -> f32 {
        self.score as f32
    }

    // WARN: placeholder for trait
    fn intronic_coords(&self) -> HashSet<(u64, u64)> {
        HashSet::new()
    }
    fn exonic_coords(&self) -> HashSet<(u64, u64)> {
        HashSet::new()
    }
    fn name(&self) -> &str {
        ""
    }
    fn block_sizes(&self) -> Vec<u64> {
        vec![]
    }
    fn block_starts(&self) -> Vec<u64> {
        vec![]
    }
    fn block_count(&self) -> u64 {
        0
    }
    fn rgb(&self) -> &str {
        ""
    }
}

pub fn bg_par_reader(bgs: Vec<PathBuf>) -> Result<(String, String), Box<dyn Error>> {
    let mut mapper: hashbrown::HashMap<String, Vec<BedGraph>> =
        unpack::<BedGraph, _>(bgs, config::OverlapType::Exon, true)
            .expect("ERROR: Could not unpack bed file!");

    let mut chroms: Vec<_> = mapper.keys().cloned().collect();
    chroms.sort_unstable();

    let mut plus = String::new();
    let mut minus = String::new();

    log::info!("INFO: Processing bedGraph fragments...");

    for chrom in chroms {
        if let Some(entries) = mapper.get_mut(&chrom) {
            let mut max_scores_plus: HashMap<(u64, u64), f64> = HashMap::new();
            let mut max_scores_minus: HashMap<(u64, u64), f64> = HashMap::new();

            for bed in entries.iter() {
                let key = (bed.start, bed.end);

                match bed.strand {
                    Strand::Forward => max_scores_plus
                        .entry(key)
                        .and_modify(|s| *s = f64::max(*s, bed.score))
                        .or_insert(bed.score),
                    Strand::Reverse => max_scores_minus
                        .entry(key)
                        .and_modify(|s| *s = f64::max(*s, bed.score))
                        .or_insert(bed.score),
                };
            }

            let mut merged_plus_entries: Vec<_> = max_scores_plus.into_iter().collect();
            merged_plus_entries.sort_unstable_by_key(|&(start, _)| start);

            let mut merged_minus_entries: Vec<_> = max_scores_minus.into_iter().collect();
            merged_minus_entries.sort_unstable_by_key(|&(start, _)| start);

            for ((start, end), score) in merged_plus_entries {
                plus.push_str(&format!("{}\t{}\t{}\t{}\n", chrom, start, end, score));
            }

            for ((start, end), score) in merged_minus_entries {
                minus.push_str(&format!("{}\t{}\t{}\t{}\n", chrom, start, end, score));
            }
        }
    }

    Ok((plus, minus))
}

pub fn get_assets_dir() -> PathBuf {
    let mut assets = std::env::current_dir().expect("Failed to get executable path");

    if !assets.ends_with(ISO_POLYA) {
        let rest = PathBuf::from(ISO_POLYA).join(ASSETS);
        assets.push(rest);

        return assets;
    } else {
        return assets.join(ASSETS);
    }
}
