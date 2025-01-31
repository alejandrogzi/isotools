use config::{BedParser, OverlapType, Sequence, Strand, SupportType, SCALE};
use hashbrown::{HashMap, HashSet};
use serde::{Deserialize, Serialize};

use std::collections::BTreeSet;

#[derive(Debug, PartialEq, Clone)]
pub struct Bed12 {
    data: GenePred,
}

#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct GenePred {
    pub name: String,
    pub chrom: String,
    pub strand: Strand,
    pub start: u64,
    pub end: u64,
    pub cds_start: u64,
    pub cds_end: u64,
    pub exons: Vec<(u64, u64)>,
    pub introns: Vec<(u64, u64)>,
    pub exon_count: usize,
    pub line: String,
    pub is_ref: bool,
}

#[derive(Debug, PartialEq, Clone)]
pub struct RefGenePred {
    pub reads: Vec<GenePred>,
    pub starts: BTreeSet<(u64, u64)>,
    pub middles: BTreeSet<(u64, u64)>,
    pub introns: BTreeSet<(u64, u64)>,
    pub bounds: (u64, u64),
    pub strand: Strand,
}

#[derive(Debug, PartialEq, Clone)]
pub struct IntronBucket {
    pub chrom: String,
    pub strand: Strand,
    pub introns: HashMap<(u64, u64), IntronPredStats>,
}

#[derive(Debug, PartialEq, Clone)]
pub struct IntronPred {
    pub chrom: String,
    pub strand: Strand,
    pub start: u64,
    pub end: u64,
    pub stats: IntronPredStats,
}

#[derive(Debug, PartialEq, Clone)]
pub struct IntronPredStats {
    pub seen: usize,                     // Frequency
    pub spanned: usize,                  // Frequency
    pub splice_ai_donor: f32,            // SpliceAi
    pub splice_ai_acceptor: f32,         // SpliceAi
    pub max_ent_donor: f32,              // MaxEntScan
    pub max_ent_acceptor: f32,           // MaxEntScan
    pub donor_sequence: String,          // MaxEntScan/TOGA-nag
    pub acceptor_sequence: String,       // MaxEntScan/TOGA-nag
    pub donor_context: Sequence,         // MaxEntScan 9-mer
    pub acceptor_context: Sequence,      // MaxEntScan 23-mer
    pub intron_position: IntronPosition, // TOGA-dependent
    pub is_toga_supported: bool,         // TOGA-dependent
    pub is_in_frame: bool,               // Miscellaneous
    pub donor_rt_context: String,        // RT-switch
    pub acceptor_rt_context: String,     // RT-switch
    pub is_rt_intron: bool,              // RT-switch
    pub is_nag_intron: bool,             // TOGA-nag
    pub support: SupportType,            // Classification
}

impl BedParser for IntronPred {
    fn parse(
        line: &str,
        _overlap: OverlapType,
        _is_ref: bool,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let data = IntronPred::read(line).expect("ERROR: Cannot parse line!");
        Ok(data)
    }

    fn chrom(&self) -> &str {
        &self.chrom.as_str()
    }

    fn coord(&self) -> (u64, u64) {
        (self.start, self.end)
    }

    // WARN: placeholder producing dummy data
    fn intronic_coords(&self) -> HashSet<(u64, u64)> {
        vec![(self.start, self.end)].iter().cloned().collect()
    }

    // WARN: placeholder producing dummy data
    fn exonic_coords(&self) -> HashSet<(u64, u64)> {
        vec![(self.start, self.end)].iter().cloned().collect()
    }
}

impl IntronPred {
    pub fn read(line: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let mut data = line.split('\t');

        let (
            chrom,
            start,
            end,
            strand,
            seen,
            spanned,
            splice_ai_donor,
            splice_ai_acceptor,
            max_ent_donor,
            max_ent_acceptor,
            donor_sequence,
            acceptor_sequence,
            donor_context,
            acceptor_context,
            intron_position,
            is_toga_supported,
            is_in_frame,
            donor_rt_context,
            acceptor_rt_context,
            is_rt_intron,
            is_nag_intron,
            support,
        ) = (
            data.next().expect("ERROR: Cannot parse chrom"),
            data.next().expect("ERROR: Cannot parse start"),
            data.next().expect("ERROR: Cannot parse end"),
            data.next().expect("ERROR: Cannot parse strand"),
            data.next().expect("ERROR: Cannot parse seen"),
            data.next().expect("ERROR: Cannot parse spanned"),
            data.next().expect("ERROR: Cannot parse splice_ai_donor"),
            data.next().expect("ERROR: Cannot parse splice_ai_acceptor"),
            data.next().expect("ERROR: Cannot parse max_ent_donor"),
            data.next().expect("ERROR: Cannot parse max_ent_acceptor"),
            data.next().expect("ERROR: Cannot parse donor_sequence"),
            data.next().expect("ERROR: Cannot parse acceptor_sequence"),
            data.next().expect("ERROR: Cannot parse donor_context"),
            data.next().expect("ERROR: Cannot parse acceptor_context"),
            data.next().expect("ERROR: Cannot parse intron_position"),
            data.next().expect("ERROR: Cannot parse is_toga_supported"),
            data.next().expect("ERROR: Cannot parse is_in_frame"),
            data.next().expect("ERROR: Cannot parse donor_rt_context"),
            data.next()
                .expect("ERROR: Cannot parse acceptor_rt_context"),
            data.next().expect("ERROR: Cannot parse is_rt_intron"),
            data.next().expect("ERROR: Cannot parse is_nag_intron"),
            data.next().expect("ERROR: Cannot parse support"),
        );

        let strand = match strand {
            "+" => Strand::Forward,
            "-" => Strand::Reverse,
            _ => panic!("ERROR: Strand is not + or -"),
        };

        let (start, end) = match strand {
            Strand::Forward => {
                let start = start.parse::<u64>().expect("ERROR: Cannot parse start");
                let end = end.parse::<u64>().expect("ERROR: Cannot parse end");

                (start, end)
            }
            Strand::Reverse => {
                let start = start.parse::<u64>().expect("ERROR: Cannot parse start");
                let end = end.parse::<u64>().expect("ERROR: Cannot parse end");

                (SCALE - end, SCALE - start)
            }
        };

        let stats = IntronPredStats::from(vec![
            seen,
            spanned,
            splice_ai_donor,
            splice_ai_acceptor,
            max_ent_donor,
            max_ent_acceptor,
            donor_sequence,
            acceptor_sequence,
            donor_context,
            acceptor_context,
            intron_position,
            is_toga_supported,
            is_in_frame,
            donor_rt_context,
            acceptor_rt_context,
            is_rt_intron,
            is_nag_intron,
            support,
        ]);

        Ok(Self {
            chrom: chrom.into(),
            strand,
            start,
            end,
            stats,
        })
    }
}

impl From<GenePred> for IntronPred {
    fn from(read: GenePred) -> Self {
        let stats = read.line.split('\t').collect::<Vec<_>>()[4..].to_vec();

        IntronPred {
            chrom: read.chrom.clone(),
            strand: read.strand.clone(),
            start: read.start,
            end: read.end,
            stats: IntronPredStats::from(stats),
        }
    }
}

impl IntronPredStats {
    pub fn new() -> Self {
        Self {
            seen: 0,
            spanned: 0,
            splice_ai_donor: 0.0,
            splice_ai_acceptor: 0.0,
            max_ent_donor: 0.0,
            max_ent_acceptor: 0.0,
            donor_sequence: String::new(),
            acceptor_sequence: String::new(),
            donor_context: Sequence::new(&[]),
            acceptor_context: Sequence::new(&[]),
            intron_position: IntronPosition::Unknown,
            is_toga_supported: false,
            is_in_frame: false,
            donor_rt_context: String::new(),
            acceptor_rt_context: String::new(),
            is_rt_intron: false,
            is_nag_intron: false,
            support: SupportType::Unclear,
        }
    }

    pub fn from(data: Vec<&str>) -> Self {
        let (
            seen,
            spanned,
            splice_ai_donor,
            splice_ai_acceptor,
            max_ent_donor,
            max_ent_acceptor,
            donor_sequence,
            acceptor_sequence,
            donor_context,
            acceptor_context,
            intron_position,
            is_toga_supported,
            is_in_frame,
            donor_rt_context,
            acceptor_rt_context,
            is_rt_intron,
            is_nag_intron,
            support,
        ) = (
            data[0].parse::<usize>().expect("ERROR: Cannot parse seen"),
            data[1]
                .parse::<usize>()
                .expect("ERROR: Cannot parse spanned"),
            data[2]
                .parse::<f32>()
                .expect("ERROR: Cannot parse splice_ai_donor"),
            data[3]
                .parse::<f32>()
                .expect("ERROR: Cannot parse splice_ai_acceptor"),
            data[4]
                .parse::<f32>()
                .expect("ERROR: Cannot parse max_ent_donor"),
            data[5]
                .parse::<f32>()
                .expect("ERROR: Cannot parse max_ent_acceptor"),
            data[6].into(),
            data[7].into(),
            Sequence::new(data[8].as_bytes()),
            Sequence::new(data[9].as_bytes()),
            match data[10] {
                "UTR" => IntronPosition::UTR,
                "CDS" => IntronPosition::CDS,
                "Mixed" => IntronPosition::Mixed,
                "Unknown" => IntronPosition::Unknown,
                _ => panic!("ERROR: Cannot parse intron_position"),
            },
            data[11]
                .parse::<bool>()
                .expect("ERROR: Cannot parse is_toga_supported"),
            data[12]
                .parse::<bool>()
                .expect("ERROR: Cannot parse is_in_frame"),
            data[13].into(),
            data[14].into(),
            data[15]
                .parse::<bool>()
                .expect("ERROR: Cannot parse is_rt_intron"),
            data[16]
                .parse::<bool>()
                .expect("ERROR: Cannot parse is_nag_intron"),
            data[17]
                .parse::<SupportType>()
                .expect("ERROR: Cannot parse support"),
        );

        Self {
            seen,
            spanned,
            splice_ai_donor,
            splice_ai_acceptor,
            max_ent_donor,
            max_ent_acceptor,
            donor_sequence,
            acceptor_sequence,
            donor_context,
            acceptor_context,
            intron_position,
            is_toga_supported,
            is_in_frame,
            donor_rt_context,
            acceptor_rt_context,
            is_rt_intron,
            is_nag_intron,
            support,
        }
    }

    pub fn fmt(&self, chr: &String, strand: &Strand, start: u64, end: u64) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            chr,
            start,
            end,
            strand,
            self.seen,
            self.spanned,
            self.splice_ai_donor,
            self.splice_ai_acceptor,
            self.max_ent_donor,
            self.max_ent_acceptor,
            self.donor_sequence,
            self.acceptor_sequence,
            self.donor_context,
            self.acceptor_context,
            self.intron_position,
            self.is_toga_supported,
            self.is_in_frame,
            self.donor_rt_context,
            self.acceptor_rt_context,
            self.is_rt_intron,
            self.is_nag_intron,
            self.support
        )
    }
}

#[derive(Debug, PartialEq, Clone, Eq, Hash)]
pub enum IntronPosition {
    UTR,
    CDS,
    Mixed,
    Unknown,
}

impl std::fmt::Display for IntronPosition {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            IntronPosition::UTR => write!(f, "UTR"),
            IntronPosition::CDS => write!(f, "CDS"),
            IntronPosition::Mixed => write!(f, "Mixed"),
            IntronPosition::Unknown => write!(f, "Unknown"),
        }
    }
}

impl GenePred {
    pub fn line(&self) -> &String {
        &self.line
    }

    pub fn name(&self) -> &String {
        &self.name
    }

    pub fn is_ref(&self) -> bool {
        self.is_ref
    }

    pub fn line_mut(&mut self) -> &mut String {
        &mut self.line
    }

    #[inline(always)]
    pub fn get_first_exon(&self) -> (u64, u64) {
        self.exons[0].clone()
    }

    #[inline(always)]
    pub fn get_middle_exons(&self) -> BTreeSet<(u64, u64)> {
        self.exons[1..].iter().cloned().collect()
    }

    #[inline(always)]
    pub fn get_exons(&self) -> BTreeSet<(u64, u64)> {
        self.exons.iter().cloned().collect()
    }

    #[inline(always)]
    pub fn get_introns(&self) -> BTreeSet<(u64, u64)> {
        self.introns.iter().cloned().collect()
    }

    pub fn get_five_utr(&self) -> (u64, u64) {
        (self.start, self.cds_start)
    }

    pub fn get_three_utr(&self) -> (u64, u64) {
        (self.cds_end, self.end)
    }

    pub fn get_split_name(&self) -> &str {
        let id: Vec<&str> = self.name.splitn(3, ".").collect();
        id[1]
    }

    #[inline(always)]
    pub fn colorline(self: Self, color: &str) -> Self {
        let nline = self.line.clone();
        let mut fields = nline.split('\t').collect::<Vec<_>>();
        fields[8] = color;
        let new_line = fields.join("\t");

        GenePred {
            line: new_line,
            name: self.name.clone(),
            chrom: self.chrom.clone(),
            strand: self.strand,
            start: self.start,
            end: self.end,
            cds_start: self.cds_start,
            cds_end: self.cds_end,
            exons: self.exons.clone(),
            introns: self.introns.clone(),
            exon_count: self.exon_count,
            is_ref: self.is_ref,
        }
    }

    pub fn mut_name_from_line(&mut self, name: &str) -> String {
        let line = self.line.clone();
        let mut fields = line.split('\t').collect::<Vec<_>>();
        fields[3] = name;

        fields.join("\t")
    }
}

impl BedParser for GenePred {
    fn parse(
        line: &str,
        overlap: OverlapType,
        is_ref: bool,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let data = Bed12::read(line, overlap, is_ref).expect("ERROR: Cannot parse line");
        Ok(data)
    }

    fn chrom(&self) -> &str {
        &self.chrom.as_str()
    }

    fn coord(&self) -> (u64, u64) {
        (self.start, self.end)
    }

    fn intronic_coords(&self) -> HashSet<(u64, u64)> {
        self.introns.iter().cloned().collect()
    }

    fn exonic_coords(&self) -> HashSet<(u64, u64)> {
        self.exons.iter().cloned().collect()
    }
}

impl From<IntronPred> for GenePred {
    fn from(intron: IntronPred) -> Self {
        GenePred {
            name: "".to_string(),
            chrom: intron.chrom.clone(),
            strand: intron.strand.clone(),
            start: intron.start,
            end: intron.end,
            cds_start: 0,
            cds_end: 0,
            exons: vec![(intron.start - 2, intron.end + 2)], // WARN: 2 +/- offset to mimic exon bounds
            introns: vec![],
            exon_count: 0,
            line: intron
                .stats
                .fmt(&intron.chrom, &intron.strand, intron.start, intron.end), // INFO: holds IntronPredStats
            is_ref: true, // INFO: always true since it is a reference file being read
        }
    }
}

impl Bed12 {
    #[inline(always)]
    pub fn read(line: &str, overlap: OverlapType, is_ref: bool) -> Result<GenePred, &'static str> {
        if line.is_empty() {
            return Err("Empty line");
        }

        let mut fields = line.split('\t');
        let (
            chrom,
            tx_start,
            tx_end,
            name,
            _,
            strand,
            cds_start,
            cds_end,
            _,
            _,
            exon_sizes,
            exon_starts,
        ) = (
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
            fields.next().ok_or("Cannot parse cds_start")?,
            fields.next().ok_or("Cannot parse cds_end")?,
            fields.next().ok_or("Cannot parse rgb")?,
            fields.next().ok_or("Cannot parse block_count")?,
            fields.next().ok_or("Cannot parse exon_sizes")?,
            fields.next().ok_or("Cannot parse exon_starts")?,
        );

        let get = |field: &str| field.parse::<u64>().map_err(|_| "Cannot parse field");
        let (tx_start, tx_end, cds_start, cds_end) =
            abs_pos(tx_start, tx_end, cds_start, cds_end, strand, get)?;

        let (exons, introns) = get_coords(
            exon_starts,
            exon_sizes,
            tx_start,
            tx_end,
            cds_start,
            cds_end,
            strand,
            overlap,
        )?;

        let mut exons = exons.iter().cloned().collect::<Vec<_>>();
        exons.sort_unstable();

        let mut introns = introns.iter().cloned().collect::<Vec<_>>();
        introns.sort_unstable();

        let exon_count = exons.len();

        let strand = match strand {
            '+' => Strand::Forward,
            '-' => Strand::Reverse,
            _ => return Err("ERROR: Strand is not + or -"),
        };

        Ok(GenePred {
            name: name.into(),
            chrom: chrom.into(),
            strand,
            start: tx_start,
            end: tx_end,
            cds_start,
            cds_end,
            exons,
            introns,
            exon_count,
            line: line.to_string(),
            is_ref,
        })
    }
}

impl BedParser for Bed12 {
    fn parse(
        line: &str,
        overlap: OverlapType,
        is_ref: bool,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let data = Bed12::read(line, overlap, is_ref).expect("ERROR: Cannot parse line");
        Ok(Bed12 { data })
    }

    fn chrom(&self) -> &str {
        &self.data.chrom.as_str()
    }

    fn coord(&self) -> (u64, u64) {
        (self.data.start, self.data.end)
    }

    fn intronic_coords(&self) -> HashSet<(u64, u64)> {
        self.data.introns.iter().cloned().collect()
    }

    fn exonic_coords(&self) -> HashSet<(u64, u64)> {
        self.data.exons.iter().cloned().collect()
    }
}

#[inline(always)]
fn get_coords(
    starts: &str,
    sizes: &str,
    tx_start: u64,
    tx_end: u64,
    cds_start: u64,
    cds_end: u64,
    strand: char,
    overlap: OverlapType,
) -> Result<(HashSet<(u64, u64)>, HashSet<(u64, u64)>), &'static str> {
    let group = |field: &str| -> Result<Vec<u64>, &'static str> {
        field
            .split(',')
            .filter_map(|num| {
                if !num.is_empty() {
                    Some(num.parse::<u64>().expect("Cannot parse number"))
                } else {
                    None
                }
            })
            .map(|num| Ok(num))
            .collect()
    };

    let ss = group(starts)?;
    let sz = group(sizes)?;

    if ss.len() != sz.len() {
        return Err("Exon start and end vectors have different lengths");
    }

    let offset = match strand {
        '+' => tx_start,
        '-' => tx_end,
        _ => return Err("Strand is not + or -"),
    };

    let exons = ss
        .iter()
        .zip(&sz)
        .map(|(&s, &z)| match strand {
            '+' => match overlap {
                OverlapType::Exon | OverlapType::Boundary => Ok((s + offset, s + z + offset)),
                OverlapType::CDS => {
                    if s + z + offset < cds_start || s + offset > cds_end {
                        return Err("ERROR: UTRs are not allowed in CDS exons!");
                    } else if s + offset < cds_start {
                        if s + z + offset > cds_end {
                            return Ok((cds_start, cds_end));
                        }
                        return Ok((cds_start, s + z + offset));
                    } else if s + z + offset > cds_end {
                        return Ok((s + offset, cds_end));
                    } else {
                        Ok((s + offset, s + z + offset))
                    }
                }
            },
            '-' => match overlap {
                OverlapType::Exon | OverlapType::Boundary => Ok((offset - s - z, offset - s)),
                OverlapType::CDS => {
                    if offset - s < cds_start || offset - s - z > cds_end {
                        return Err("UTRs are not allowed in CDS exons");
                    } else if offset - s - z < cds_start {
                        if offset - s > cds_end {
                            return Ok((cds_start, cds_end));
                        }
                        return Ok((cds_start, offset - s));
                    } else if offset - s > cds_end {
                        return Ok((offset - s - z, cds_end));
                    } else {
                        Ok((offset - s - z, offset - s))
                    }
                }
            },
            _ => return Err("Strand is not + or -"),
        })
        .filter_map(Result::ok)
        .collect::<HashSet<_>>();

    let introns = gapper(&exons);

    Ok((exons, introns))
}

#[inline(always)]
fn abs_pos(
    tx_start: &str,
    tx_end: &str,
    cds_start: &str,
    cds_end: &str,
    strand: char,
    get: impl Fn(&str) -> Result<u64, &'static str>,
) -> Result<(u64, u64, u64, u64), &'static str> {
    match strand {
        '+' => {
            let tx_start = get(tx_start)?;
            let tx_end = get(tx_end)?;
            let cds_start = get(cds_start)?;
            let cds_end = get(cds_end)?;

            Ok((tx_start, tx_end, cds_start, cds_end))
        }
        '-' => {
            let tx_start = get(tx_start)?;
            let tx_end = get(tx_end)?;
            let cds_start = get(cds_start)?;
            let cds_end = get(cds_end)?;

            Ok((
                SCALE - tx_end,
                SCALE - tx_start,
                SCALE - cds_end,
                SCALE - cds_start,
            ))
        }
        _ => Err("Strand is not + or -"),
    }
}

fn gapper(intervals: &HashSet<(u64, u64)>) -> HashSet<(u64, u64)> {
    let mut vintervals: Vec<(u64, u64)> = intervals.iter().copied().collect();
    vintervals.sort_by(|a, b| a.0.cmp(&b.0));

    let mut gaps = HashSet::with_capacity(vintervals.len());
    for window in vintervals.windows(2) {
        if let [prev, next] = window {
            let gap_start = prev.1 + 1;
            let gap_end = next.0 - 1;

            if gap_start < gap_end {
                gaps.insert((gap_start, gap_end));
            }
        }
    }

    gaps
}

impl RefGenePred {
    pub fn new(
        reads: Vec<GenePred>,
        starts: BTreeSet<(u64, u64)>,
        middles: BTreeSet<(u64, u64)>,
        introns: BTreeSet<(u64, u64)>,
        bounds: (u64, u64),
        strand: Strand,
    ) -> Self {
        Self {
            reads,
            starts,
            middles,
            introns,
            bounds,
            strand,
        }
    }

    pub fn get_names(&self) -> HashSet<String> {
        let mut names = HashSet::new();
        self.reads.iter().for_each(|read| {
            names.insert(read.name.clone());
        });

        names
    }

    pub fn get_names_split(&self) -> HashSet<&str> {
        let mut names = HashSet::new();
        self.reads.iter().for_each(|read| {
            let id: Vec<&str> = read.name.splitn(3, ".").collect();
            names.insert(id[1]);
        });

        names
    }

    pub fn merge_names(&self) -> String {
        let names = self.get_names_split().into_iter().collect::<Vec<&str>>();
        names.join(".")
    }

    pub fn smash_exons_by_name(&self) -> Vec<Vec<(u64, u64)>> {
        let names = self.get_names_split();
        let mut smashed = Vec::new();
        for name in names {
            let exons = self
                .reads
                .iter()
                .filter(|read| read.get_split_name() == name)
                .map(|read| read.exons.clone())
                .flatten()
                .collect::<BTreeSet<_>>();

            smashed.push(exons.into_iter().collect());
        }

        smashed
    }

    pub fn smash_exons_by_name_split(&self) -> HashMap<&str, Vec<(u64, u64)>> {
        let names = self.get_names_split();
        let mut smashed = HashMap::new();
        for name in names {
            let exons = self
                .reads
                .iter()
                .filter(|read| read.get_split_name() == name)
                .map(|read| read.exons.clone())
                .flatten()
                .collect::<BTreeSet<_>>();

            smashed.insert(name, exons.into_iter().collect());
        }

        smashed
    }

    #[inline(always)]
    pub fn from(reads: Vec<GenePred>) -> Self {
        let mut starts = BTreeSet::new();
        let mut middles = BTreeSet::new();
        let mut introns = BTreeSet::new();
        let mut bounds = (u64::MAX, 0);

        let strand = if !reads.is_empty() {
            reads[0].strand.clone()
        } else {
            // WARN: case where trying to create a RefGenePred from an empty Vec<GenePred>
            Strand::Forward
        };

        for read in &reads {
            bounds.0 = bounds.0.min(read.start);
            bounds.1 = bounds.1.max(read.end);

            let read_start = read.get_first_exon();
            introns.extend(read.get_introns());

            if !middles.is_empty() {
                let overlap = middles.iter().any(|middle| {
                    let (mid_exon_start, mid_exon_end) = middle;
                    (read_start.0 >= *mid_exon_start && read_start.0 < *mid_exon_end)
                        || (read_start.1 > *mid_exon_start && read_start.1 <= *mid_exon_end)
                        || (read_start.0 < *mid_exon_start && read_start.1 > *mid_exon_end)
                });

                if !overlap {
                    starts.insert(read_start);
                } else {
                    starts.clone().iter().for_each(|start| {
                        let (five_end_start, five_end_end) = start;
                        if read_start.0 >= *five_end_start && read_start.0 < *five_end_end {
                            starts.insert(read_start);
                        }
                    });
                }
            } else {
                starts.insert(read_start);
            }

            if read.exon_count > 1 {
                middles.extend(read.get_middle_exons());
            }
        }

        Self::new(reads, starts, middles, introns, bounds, strand)
    }
}

impl IntronBucket {
    // Vec<GenePred> -> IntronBucket
    #[inline(always)]
    #[allow(unused_assignments)]
    pub fn from(reads: Vec<GenePred>, toga: Vec<GenePred>) -> Self {
        let mut introns = HashMap::new();
        let mut toga_introns = HashSet::new();

        let mut chr = String::new();
        let mut strand = Strand::Forward;

        let mut cds_start = 0;
        let mut cds_end = 0;

        if !toga.is_empty() {
            chr = toga[0].chrom.clone();
            strand = toga[0].strand.clone();

            for projection in toga {
                let introns = projection.introns;

                cds_start = cds_start.min(projection.cds_start);
                cds_end = cds_end.max(projection.cds_end);

                for intron in introns {
                    // INFO: fmt -> intron coord -> (CDS start, CDS end)
                    toga_introns.insert(intron);
                }
            }
        } else {
            chr = reads[0].chrom.clone();
            strand = reads[0].strand.clone();
        }

        for read in &reads {
            for ref_intron in &read.introns {
                if introns.contains_key(ref_intron) {
                    let stats: &mut IntronPredStats = introns
                        .get_mut(ref_intron)
                        .expect("ERROR: Cannot get intron handle. This is a bug!");
                    stats.seen += 1;
                } else {
                    let is_toga_supported = toga_introns.contains(ref_intron);
                    if is_toga_supported {
                        toga_introns.remove(ref_intron);
                    }

                    let position = if is_toga_supported {
                        IntronPosition::CDS
                    } else if (ref_intron.0 > cds_end || ref_intron.1 < cds_start) && cds_end > 0 {
                        IntronPosition::UTR
                    } else if (ref_intron.0 < cds_start && ref_intron.1 < cds_end)
                        || (ref_intron.0 < cds_end && ref_intron.1 > cds_end)
                    {
                        IntronPosition::Mixed
                    } else {
                        IntronPosition::Unknown
                    };

                    let in_frame = (ref_intron.1 - ref_intron.0) % 3 == 0;

                    let stats = IntronPredStats {
                        seen: 1,
                        spanned: 0,
                        splice_ai_donor: 0.0,
                        splice_ai_acceptor: 0.0,
                        max_ent_donor: 0.0,
                        max_ent_acceptor: 0.0,
                        donor_sequence: String::new(),
                        acceptor_sequence: String::new(),
                        donor_context: Sequence::new(&[]),
                        acceptor_context: Sequence::new(&[]),
                        intron_position: position,
                        is_toga_supported,
                        is_in_frame: in_frame,
                        donor_rt_context: String::new(),
                        acceptor_rt_context: String::new(),
                        is_rt_intron: false,
                        is_nag_intron: false,
                        support: SupportType::Unclear,
                    };

                    introns.insert(ref_intron.clone(), stats);
                }
            }

            for (intron, stats) in introns.iter_mut() {
                if read.start <= intron.0 && read.end >= intron.1 {
                    stats.spanned += 1;
                }
            }
        }

        // INFO: add remaining TOGA introns if any and include flag to indentify them
        if !toga_introns.is_empty() {
            for toga_intron in toga_introns {
                let stats = IntronPredStats {
                    seen: 0,
                    spanned: 0,
                    splice_ai_donor: 0.0,
                    splice_ai_acceptor: 0.0,
                    max_ent_donor: 0.0,
                    max_ent_acceptor: 0.0,
                    donor_sequence: String::new(),
                    acceptor_sequence: String::new(),
                    donor_context: Sequence::new(&[]),
                    acceptor_context: Sequence::new(&[]),
                    intron_position: IntronPosition::Unknown,
                    is_toga_supported: true, // INFO: flag to identify TOGA introns
                    is_in_frame: false,
                    donor_rt_context: String::new(),
                    acceptor_rt_context: String::new(),
                    is_rt_intron: false,
                    is_nag_intron: false,
                    support: SupportType::Unclear, // INFO: flag to identify TOGA introns
                };

                introns.insert(toga_intron, stats);
            }
        }

        Self {
            chrom: chr,
            strand,
            introns,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bed12_abs_pos_plus() {
        let tx_start = "10";
        let tx_end = "20";
        let cds_start = "10";
        let cds_end = "20";
        let strand = '+';

        let (tx_start, tx_end, cds_start, cds_end) =
            abs_pos(tx_start, tx_end, cds_start, cds_end, strand, |x| {
                Ok(x.parse().unwrap())
            })
            .unwrap();

        assert_eq!(tx_start, 10);
        assert_eq!(tx_end, 20);
        assert_eq!(cds_start, 10);
        assert_eq!(cds_end, 20);
    }

    #[test]
    fn test_bed12_abs_pos_minus() {
        let tx_start = "10";
        let tx_end = "20";
        let cds_start = "10";
        let cds_end = "20";
        let strand = '-';

        let (tx_start, tx_end, cds_start, cds_end) =
            abs_pos(tx_start, tx_end, cds_start, cds_end, strand, |x| {
                Ok(x.parse().unwrap())
            })
            .unwrap();

        assert_eq!(tx_start, SCALE - 20);
        assert_eq!(tx_end, SCALE - 10);
        assert_eq!(cds_start, SCALE - 20);
        assert_eq!(cds_end, SCALE - 10);
    }

    #[test]
    fn test_bed12_get_coords_and_gapper_plus() {
        let start = "0,30";
        let size = "10,10";
        let tx_start = 10;
        let tx_end = 50;
        let cds_start = 15;
        let cds_end = 45;
        let strand = '+';
        let overlap = OverlapType::CDS;

        let (exons, introns) = get_coords(
            start, size, tx_start, tx_end, cds_start, cds_end, strand, overlap,
        )
        .unwrap();

        let mut exons = exons.iter().cloned().collect::<Vec<_>>();
        exons.sort_unstable_by(|a, b| a.cmp(b));

        let mut introns = introns.iter().cloned().collect::<Vec<_>>();
        introns.sort_unstable_by(|a, b| a.cmp(b));

        assert_eq!(
            exons,
            [(15, 20), (40, 45)].iter().cloned().collect::<Vec<_>>()
        );
        assert_eq!(introns, [(21, 39)].iter().cloned().collect::<Vec<_>>());
    }

    #[test]
    fn test_bed12_get_coords_and_gapper_minus() {
        let start = "0,20,40,60,80";
        let size = "10,10,10,10,10";
        let tx_start = "10";
        let tx_end = "100";
        let cds_start = "30";
        let cds_end = "80";
        let strand = '-';
        let overlap = OverlapType::CDS;

        let get = |field: &str| field.parse::<u64>().map_err(|_| "Cannot parse field");
        let (tx_start, tx_end, cds_start, cds_end) =
            abs_pos(tx_start, tx_end, cds_start, cds_end, strand, get).unwrap();

        let (exons, introns) = get_coords(
            start, size, tx_start, tx_end, cds_start, cds_end, strand, overlap,
        )
        .unwrap();

        let mut exons = exons.iter().cloned().collect::<Vec<_>>();
        exons.sort_unstable_by(|a, b| a.cmp(b));

        let mut introns = introns.iter().cloned().collect::<Vec<_>>();
        introns.sort_unstable_by(|a, b| a.cmp(b));

        assert_eq!(
            exons,
            [
                (99999999920, 99999999930),
                (99999999940, 99999999950),
                (99999999960, 99999999970)
            ]
            .iter()
            .cloned()
            .collect::<Vec<_>>()
        );
        assert_eq!(
            introns,
            [(99999999931, 99999999939), (99999999951, 99999999959)]
                .iter()
                .cloned()
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_bed12_get_coords_and_gapper_plus_nested_utr() {
        let start = "0,20,40,60,80";
        let size = "10,10,10,10,10";
        let tx_start = 10;
        let tx_end = 100;
        let cds_start = 15;
        let cds_end = 95;
        let strand = '+';
        let overlap = OverlapType::CDS;

        let (exons, introns) = get_coords(
            start, size, tx_start, tx_end, cds_start, cds_end, strand, overlap,
        )
        .unwrap();

        let mut exons = exons.iter().cloned().collect::<Vec<_>>();
        exons.sort_unstable_by(|a, b| a.cmp(b));

        let mut introns = introns.iter().cloned().collect::<Vec<_>>();
        introns.sort_unstable_by(|a, b| a.cmp(b));

        assert_eq!(
            exons,
            [(15, 20), (30, 40), (50, 60), (70, 80), (90, 95)]
                .iter()
                .cloned()
                .collect::<Vec<_>>()
        );
        assert_eq!(
            introns,
            [(21, 29), (41, 49), (61, 69), (81, 89)]
                .iter()
                .cloned()
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_bed12_get_coords_and_gapper_minus_nested_utr() {
        let start = "0,20,40,60";
        let size = "10,10,10,10";
        let tx_start = "10";
        let tx_end = "80";
        let cds_start = "15";
        let cds_end = "75";
        let strand = '-';
        let overlap = OverlapType::CDS;

        let get = |field: &str| field.parse::<u64>().map_err(|_| "Cannot parse field");
        let (tx_start, tx_end, cds_start, cds_end) =
            abs_pos(tx_start, tx_end, cds_start, cds_end, strand, get).unwrap();

        let (exons, introns) = get_coords(
            start, size, tx_start, tx_end, cds_start, cds_end, strand, overlap,
        )
        .unwrap();

        let mut exons = exons.iter().cloned().collect::<Vec<_>>();
        exons.sort_unstable_by(|a, b| a.cmp(b));

        let mut introns = introns.iter().cloned().collect::<Vec<_>>();
        introns.sort_unstable_by(|a, b| a.cmp(b));

        assert_eq!(
            exons,
            [
                (99999999925, 99999999930),
                (99999999940, 99999999950),
                (99999999960, 99999999970),
                (99999999980, 99999999985)
            ]
            .iter()
            .cloned()
            .collect::<Vec<_>>()
        );
        assert_eq!(
            introns,
            [
                (99999999931, 99999999939),
                (99999999951, 99999999959),
                (99999999971, 99999999979)
            ]
            .iter()
            .cloned()
            .collect::<Vec<_>>()
        );
    }
}
