use config::{BedParser, OverlapType, Strand};
use hashbrown::HashSet;
use serde::{Deserialize, Serialize};

pub const EXPANSION_SIZE: u64 = 150; // 150bp

#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct PolyAPred {
    pub name: String,
    pub chrom: String,
    pub strand: Strand,
    pub start: u64,
    pub end: u64,
}

impl PolyAPred {
    #[inline(always)]
    pub fn read(line: &str, _: OverlapType, _: bool) -> Result<PolyAPred, &'static str> {
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

        Ok(PolyAPred {
            name: name.into(),
            chrom: chrom.into(),
            strand,
            start,
            end,
        })
    }
}

impl BedParser for PolyAPred {
    fn parse(
        line: &str,
        overlap: OverlapType,
        is_ref: bool,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let data = PolyAPred::read(line, overlap, is_ref).expect("ERROR: Cannot parse line");
        Ok(data)
    }

    fn chrom(&self) -> &str {
        &self.chrom.as_str()
    }

    fn coord(&self) -> (u64, u64) {
        (self.start, self.end)
    }

    // WARN: placeholder for trait
    fn intronic_coords(&self) -> HashSet<(u64, u64)> {
        HashSet::new()
    }

    // WARN: placeholder for trait
    fn exonic_coords(&self) -> HashSet<(u64, u64)> {
        HashSet::new()
    }
}
