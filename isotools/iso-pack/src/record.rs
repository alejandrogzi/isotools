use config::{BedParser, OverlapType, Sequence, Strand, SupportType, BIG_SEP, SCALE, SEP};
use hashbrown::{HashMap, HashSet};
use serde::{Deserialize, Serialize};

use std::collections::BTreeSet;

#[derive(Debug, PartialEq, Clone)]
pub struct Bed12 {
    pub data: GenePred,
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

impl GenePred {
    /// Maps ORF transcript coordinates to absolute genomic coordinates (strand-aware).
    ///
    /// This function takes ORF coordinates defined in **linear transcript space**
    /// (as produced by a CDS/ORF predictor) and returns the corresponding genomic
    /// start and end coordinates of the coding region, accounting for exon structure
    /// and strand orientation. For transcripts on the reverse strand, it assumes
    /// coordinates have been **strand-normalized** using a fixed scale (`SCALE`) and
    /// reverses them appropriately.
    ///
    /// # Arguments
    ///
    /// * `orf_start` - The start of the ORF in transcript (spliced) coordinates.
    /// * `orf_end` - The end of the ORF in transcript (spliced) coordinates. This is exclusive.
    ///
    /// # Returns
    ///
    /// * `Option<(u64, u64)>` - A tuple of genomic start and end coordinates corresponding
    ///   to the CDS region. Returns `Some((start, end))` if mapping succeeds; otherwise,
    ///   may return `None` if the coordinates do not map within the exon structure.
    ///
    /// # Strand Behavior
    ///
    /// * On the `Strand::Forward`, transcript coordinates map left-to-right across exons.
    /// * On the `Strand::Reverse`, the exons are reversed and mapped using:
    ///   `SCALE - (exon_end - offset)`, where `SCALE` is a user-defined upper genomic bound
    ///   used for coordinate normalization.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let mut gene = GenePred { /* initialized */ };
    /// let orf_start = 15;
    /// let orf_end = 45;
    /// if let Some((cds_start, cds_end)) = gene.map_absolute_cds(orf_start, orf_end) {
    ///     println!("Genomic CDS coordinates: {} - {}", cds_start, cds_end);
    /// }
    /// ```
    pub fn map_absolute_cds(&mut self, orf_start: u64, orf_end: u64) -> Option<(u64, u64)> {
        let exons = match self.strand {
            Strand::Forward => &self.exons,
            Strand::Reverse => &self.exons.iter().rev().map(|&(s, e)| (s, e)).collect(),
        };

        let mut transcript_offset = 0;
        let mut cds_start = 0;
        let mut cds_end = 0;

        for (exon_start, exon_end) in exons {
            let exon_len = exon_end - exon_start;
            let exon_tx_start = transcript_offset;
            let exon_tx_end = transcript_offset + exon_len;

            // INFO: ORF start falls in this exon
            if orf_start >= exon_tx_start && orf_start < exon_tx_end {
                let offset = orf_start - exon_tx_start;
                let g_start = match self.strand {
                    Strand::Forward => exon_start + offset,
                    Strand::Reverse => SCALE - (exon_end - offset), // INFO: reverse strand -> walk from right to left
                };
                cds_start = g_start;
            }

            // INFO: ORF end falls in this exon
            if orf_end > exon_tx_start && orf_end <= exon_tx_end {
                let offset = orf_end - exon_tx_start;
                let g_end = match self.strand {
                    Strand::Forward => exon_start + offset,
                    Strand::Reverse => SCALE - (exon_end - offset),
                };
                cds_end = g_end;
                break; // INFO: early exit, done mapping
            }

            transcript_offset += exon_len;
        }

        // WARN: the only case where both coords should be discarded is when
        // WARN: the ORF start is out of bounds
        if cds_start == 0 {
            // INFO: ORF start or end is not within the exon structure
            log::warn!(
                "WARN: CDS start ({}) is 0 for ORF coordinates: {}-{} in {} strand. Exons: {:?}. Will skip this ORF!",
                cds_start,
                orf_start,
                orf_end,
                self.strand,
                self.exons
            );
            return None;
        } else if cds_end == 0 {
            // INFO: ORF end is out of bounds but we extend the CDS end to the transcript end
            log::warn!("ERROR: CDS start ({}) is greater than or equal to CDS end ({}) for ORF coordinates: {}-{} in {} strand. Exons: {:?}",
                    cds_start, cds_end, orf_start, orf_end, self.strand, self.exons);
            return Some((cds_start, self.end));
        }

        Some((cds_start, cds_end))
    }
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

    fn start(&self) -> u64 {
        self.start
    }

    fn end(&self) -> u64 {
        self.end
    }

    fn strand(&self) -> Strand {
        self.strand.clone()
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
    fn name(&self) -> &str {
        ""
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
    pub fn colorline(self: &Self, color: &str) -> Self {
        let nline = self.line.clone();
        let mut fields = nline.split('\t').collect::<Vec<_>>();
        fields[8] = color;
        let new_line = fields.join("\t");

        GenePred {
            line: new_line,
            name: self.name.clone(),
            chrom: self.chrom.clone(),
            strand: self.strand.clone(),
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

    /// Modifies the field at the given index with the provided value.
    ///
    /// # Arguments
    ///
    /// * `field` - The index of the field to modify (0-based) -> BedColumn.
    /// * `value` - The new value to set for the field.
    ///
    /// # Example
    ///
    /// ```rust
    /// let mut gene_pred = GenePred::new();
    /// gene_pred.modify_field(3, "new_value");
    ///
    /// assert_eq!(gene_pred.line(), "new_value");
    /// ```
    #[inline(always)]
    pub fn modify_field(&mut self, field: usize, value: &str) {
        let mut tab_indices = Vec::with_capacity(field + 2);
        for (i, b) in self.line.bytes().enumerate() {
            if b == b'\t' {
                tab_indices.push(i);
                if tab_indices.len() > field + 1 {
                    break;
                }
            }
        }

        let start = if field == 0 {
            0
        } else {
            tab_indices
                .get(field - 1)
                .map_or(self.line.len(), |&i| i + 1)
        };

        let end = tab_indices
            .get(field)
            .copied()
            .unwrap_or_else(|| self.line.len());

        let mut line = String::with_capacity(self.line.len() - (end - start) + value.len());
        line.push_str(&self.line[..start]);
        line.push_str(value);
        line.push_str(&self.line[end..]);

        self.line = line.clone();
    }

    pub fn mut_name_from_line(&mut self, name: &str) -> String {
        let line = self.line.clone();
        let mut fields = line.split('\t').collect::<Vec<_>>();
        fields[3] = name;

        fields.join("\t")
    }

    // WARN: modifies CDStart and CDEnd
    pub fn insert_utr(&mut self, start: u64, end: u64) -> &mut Self {
        self.cds_start = start;
        self.cds_end = end;

        self
    }

    // WARN: this needs to be refactored -> current impl is just a solution for iso-orf!
    pub fn construct_bed_line(&mut self) -> String {
        let line = self.line.clone();
        let mut fields = line.split('\t').collect::<Vec<_>>();

        match self.strand {
            Strand::Forward => {
                let mut s_binding = self.cds_start.to_string();
                let mut e_binding = self.cds_end.to_string();

                if self.cds_start < self.start {
                    s_binding = self.start.to_string();
                }

                if self.cds_end > self.end {
                    e_binding = self.end.to_string();
                }

                fields[6] = s_binding.as_str();
                fields[7] = e_binding.as_str();

                return fields.join("\t");
            }
            Strand::Reverse => {
                let mut s_binding = (SCALE - self.cds_end).to_string();
                let mut e_binding = (SCALE - self.cds_start).to_string();

                if self.cds_start < self.start {
                    s_binding = self.start.to_string();
                }

                if self.cds_end > self.end {
                    e_binding = self.end.to_string();
                }

                fields[6] = s_binding.as_str();
                fields[7] = e_binding.as_str();

                return fields.join("\t");
            }
        }
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

    fn start(&self) -> u64 {
        self.start
    }

    fn end(&self) -> u64 {
        self.end
    }

    fn cds_start(&self) -> u64 {
        self.cds_start
    }

    fn cds_end(&self) -> u64 {
        self.cds_end
    }

    fn strand(&self) -> Strand {
        self.strand.clone()
    }

    fn name(&self) -> &str {
        &self.name
    }

    fn intronic_coords(&self) -> HashSet<(u64, u64)> {
        self.introns.iter().cloned().collect()
    }

    fn exonic_coords(&self) -> HashSet<(u64, u64)> {
        self.exons.iter().cloned().collect()
    }

    fn block_sizes(&self) -> Vec<u64> {
        self.exons.iter().map(|(s, e)| e - s).collect()
    }

    fn block_starts(&self) -> Vec<u64> {
        self.exons.iter().map(|(s, _)| *s).collect()
    }

    fn block_count(&self) -> u64 {
        self.exon_count as u64
    }

    // WARN: placeholder for trait
    fn score(&self) -> f32 {
        0.0
    }
    fn rgb(&self) -> &str {
        ""
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

impl From<Bed6> for GenePred {
    fn from(record: Bed6) -> Self {
        GenePred {
            name: record.id,
            chrom: record.chrom,
            strand: record.strand,
            start: record.coord.0,
            end: record.coord.1,
            cds_start: record.coord.0,
            cds_end: record.coord.1,
            exons: vec![record.coord],
            introns: vec![],
            exon_count: 0,
            line: "".to_string(),
            is_ref: false,
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

    fn start(&self) -> u64 {
        self.data.start
    }

    fn end(&self) -> u64 {
        self.data.end
    }

    fn cds_start(&self) -> u64 {
        self.data.cds_start
    }

    fn cds_end(&self) -> u64 {
        self.data.cds_end
    }

    fn strand(&self) -> Strand {
        self.data.strand.clone()
    }

    fn name(&self) -> &str {
        &self.data.name
    }

    fn block_sizes(&self) -> Vec<u64> {
        self.data.exons.iter().map(|(s, e)| e - s).collect()
    }

    fn block_starts(&self) -> Vec<u64> {
        self.data.exons.iter().map(|(s, _)| *s).collect()
    }

    fn block_count(&self) -> u64 {
        self.data.exon_count as u64
    }

    // WARN: placeholder for trait
    fn score(&self) -> f32 {
        0.0
    }
    fn rgb(&self) -> &str {
        ""
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
                    Some(num.parse::<u64>().unwrap_or_else(|e| {
                        panic!(
                            "ERROR: Cannot parse number -> {e}. {starts}
                                {sizes}
                                {tx_start}
                                {tx_end}
                                {strand}"
                        )
                    }))
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
                OverlapType::Exon | OverlapType::Boundary | OverlapType::CDSBound => {
                    Ok((s + offset, s + z + offset))
                }
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
                OverlapType::Exon | OverlapType::Boundary | OverlapType::CDSBound => {
                    Ok((offset - s - z, offset - s))
                }
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

/// Convert relative positions to absolute positions
///
/// # Arguments
///
/// * `tx_start` - Start of transcript
/// * `tx_end` - End of transcript
/// * `cds_start` - Start of CDS
/// * `cds_end` - End of CDS
/// * `strand` - Strand of transcript
/// * `get` - Function to parse string to u64
///
/// # Returns
///
/// * `Ok` - Tuple of absolute positions
/// * `Err` - Error message
///
/// # Example
///
/// ```rust, no_run
/// let tx_start = "0";
/// let tx_end = "100";
/// let cds_start = "10";
/// let cds_end = "90";
/// let strand = '+';
/// let get = |field: &str| field.parse::<u64>().map_err(|_| "Cannot parse field");
///
/// let (tx_start, tx_end, cds_start, cds_end) = abs_pos(tx_start, tx_end, cds_start, cds_end, strand, get).unwrap();
///
/// assert_eq!(tx_start, 0);
/// assert_eq!(tx_end, 100);
/// assert_eq!(cds_start, 10);
/// assert_eq!(cds_end, 90);
/// ```
#[inline(always)]
pub fn abs_pos(
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

/// Calculate intron coordinates based on exon boundaries
///
/// # Arguments
///
/// * `intervals` - HashSet of exon boundaries
///
/// # Returns
///
/// * `HashSet` of intron boundaries
///
/// # Example
///
/// ```rust, no_run
/// let intervals = HashSet::new();
/// intervals.insert((0, 10));
/// intervals.insert((20, 30));
/// intervals.insert((40, 50));
///
/// let introns = gapper(&intervals);
///
/// assert_eq!(introns.len(), 2);
/// ```
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

    // WARN: returns HashSet instead of Vec in inner collection!
    pub fn smash_introns_by_name(&self) -> Vec<HashSet<(u64, u64)>> {
        let names = self.get_names_split();
        let mut smashed = Vec::new();
        for name in names {
            let introns = self
                .reads
                .iter()
                .filter(|read| read.get_split_name() == name)
                .map(|read| read.introns.clone())
                .flatten()
                .collect::<HashSet<_>>();

            // smashed.push(introns.into_iter().collect());
            smashed.push(introns);
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

        let mut cds_start = u64::MAX;
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
                        // INFO: removing TOGA intron from the set
                        // INFO: to allow later inserting of TOGA unique introns
                        toga_introns.remove(ref_intron);
                    }

                    // WARN: an edge case for TOGA vs Iso comprise an unseen intron
                    // WARN: that falls inside the CDS region but does not match any
                    // WARN: TOGA splice sites
                    let position = if is_toga_supported
                        || (cds_start < ref_intron.0 && ref_intron.1 < cds_end)
                    {
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
        }

        // INFO: counting spanned introns requires an additional loop,
        // INFO: problem arises when an intron appears after k-reads
        // INFO: without it have been seen already, leading to low and wrong
        // INFO: counts
        for read in &reads {
            for (intron, stats) in introns.iter_mut() {
                if read.start <= intron.0 && intron.1 <= read.end {
                    stats.spanned += 1;
                }
            }
        }

        // INFO: add remaining TOGA introns if any and include flag to indentify them
        if !toga_introns.is_empty() {
            for toga_intron in toga_introns {
                let stats = IntronPredStats {
                    seen: 0,    // INFO: is TOGA unique
                    spanned: 0, // INFO: is TOGA unique
                    splice_ai_donor: 0.0,
                    splice_ai_acceptor: 0.0,
                    max_ent_donor: 0.0,
                    max_ent_acceptor: 0.0,
                    donor_sequence: String::new(),
                    acceptor_sequence: String::new(),
                    donor_context: Sequence::new(&[]),
                    acceptor_context: Sequence::new(&[]),
                    intron_position: IntronPosition::CDS, // INFO: is always CDS -> flag
                    is_toga_supported: true,              // INFO: flag to identify TOGA introns
                    is_in_frame: false,
                    donor_rt_context: String::new(),
                    acceptor_rt_context: String::new(),
                    is_rt_intron: false,
                    is_nag_intron: false,
                    // WARN: inside iso-intron this will be interpreted as SPLICED because of TOGA
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

/// Bed4
///
/// A struct to hold the Bed4 data
///
/// # Attributes
///
/// * `chrom` - The chromosome of the read
/// * `coord` - The coordinates of the read
/// * `id` - The name of the read
///
/// # Example
///
/// ```rust, no_run
/// let line = "chr1\t1000\t2000\tread1";
/// let bed4 = Bed4::new(line).unwrap();
///
/// assert_eq!(bed4.chrom, "chr1");
/// assert_eq!(bed4.coord, (1000, 2000));
/// assert_eq!(bed4.id, "read1");
/// ```
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
        _overlap: OverlapType,
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

    fn start(&self) -> u64 {
        self.coord.0
    }

    fn end(&self) -> u64 {
        self.coord.1
    }

    fn cds_start(&self) -> u64 {
        self.coord.0
    }

    fn cds_end(&self) -> u64 {
        self.coord.1
    }

    fn name(&self) -> &str {
        &self.id
    }

    // WARN: placeholder for Bed4
    fn intronic_coords(&self) -> HashSet<(u64, u64)> {
        let mut introns = HashSet::new();
        introns.insert(self.coord);
        introns
    }
    fn exonic_coords(&self) -> HashSet<(u64, u64)> {
        let mut exons = HashSet::new();
        exons.insert(self.coord);
        exons
    }
    fn strand(&self) -> Strand {
        Strand::Forward
    }
    fn block_sizes(&self) -> Vec<u64> {
        vec![self.coord.1 - self.coord.0]
    }
    fn block_starts(&self) -> Vec<u64> {
        vec![self.coord.0]
    }
    fn block_count(&self) -> u64 {
        1
    }
    fn score(&self) -> f32 {
        0.0
    }
    fn rgb(&self) -> &str {
        ""
    }
}

/// Bed6
///
/// A struct to hold the Bed6 data
///
/// # Attributes
///
/// * `chrom` - The chromosome of the read
/// * `coord` - The coordinates of the read
/// * `id` - The name of the read
/// * `strand` - The strand of the read
/// * `score` - The score of the read
///
/// # Example
///
/// ```rust, no_run
/// let line = "chr1\t1000\t2000\tread1\t0.0\t+";
/// let bed6 = Bed6::new(line).unwrap();
///
/// assert_eq!(bed6.chrom, "chr1");
/// assert_eq!(bed6.coord, (1000, 2000));
/// assert_eq!(bed6.id, "read1");
/// assert_eq!(bed6.strand, Strand::Forward);
/// assert_eq!(bed6.score, 0.0);
/// ```
#[derive(Debug, PartialEq, Clone)]
// WARN: will cover any high-order bed file [6,8,12]
pub struct Bed6 {
    pub chrom: String,
    pub coord: (u64, u64),
    pub id: String,
    pub strand: Strand,
    pub score: f32,
}

impl Bed6 {
    pub fn new(line: String) -> Result<Bed6, Box<dyn std::error::Error>> {
        if line.is_empty() {
            return Err("ERROR: Empty line in .bed!".into());
        }

        let mut fields = line.split('\t');
        let get = |field: &str| field.parse::<u64>().map_err(|_| "Cannot parse field");

        let (chrom, start, end, id, score, strand) = (
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
                let (start, end) = (SCALE - end, SCALE - start);
                if start > end {
                    return Err("ERROR: Start is less than end!".into());
                }

                (start, end)
            }
        };

        let score = score
            .parse::<f32>()
            .map_err(|_| "ERROR: Cannot parse score")?;

        Ok(Bed6 {
            chrom,
            coord: (start, end),
            id,
            strand,
            score,
        })
    }

    pub fn from(chrom: String, start: u64, end: u64, id: String, strand: Strand) -> Bed6 {
        Bed6 {
            chrom,
            coord: (start, end),
            id,
            strand,
            score: 0.0,
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
        _overlap: OverlapType,
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

    fn start(&self) -> u64 {
        self.coord.0
    }

    fn end(&self) -> u64 {
        self.coord.1
    }

    fn cds_start(&self) -> u64 {
        self.coord.0
    }

    fn cds_end(&self) -> u64 {
        self.coord.1
    }

    fn name(&self) -> &str {
        &self.id
    }

    fn strand(&self) -> Strand {
        self.strand.clone()
    }

    fn score(&self) -> f32 {
        self.score
    }

    // WARN: placeholder for Bed4
    fn intronic_coords(&self) -> HashSet<(u64, u64)> {
        let mut introns = HashSet::new();
        introns.insert(self.coord);
        introns
    }
    fn exonic_coords(&self) -> HashSet<(u64, u64)> {
        let mut exons = HashSet::new();
        exons.insert(self.coord);
        exons
    }
    fn block_sizes(&self) -> Vec<u64> {
        vec![self.coord.1 - self.coord.0]
    }
    fn block_starts(&self) -> Vec<u64> {
        vec![self.coord.0]
    }
    fn block_count(&self) -> u64 {
        1
    }
    fn rgb(&self) -> &str {
        ""
    }
}

/// PolyAPred
///
/// A struct to hold the PolyAPred data
///
/// # Attributes
///
/// * `name` - The name of the read
/// * `chrom` - The chromosome of the read
/// * `strand` - The strand of the read
/// * `start` - The start position of the read
/// * `end` - The end position of the read
/// * `cds_start` - The start position of the CDS
/// * `cds_end` - The end position of the CDS
/// * `clip` - The number of clipped bases
/// * `clipped_a` - The number of clipped A's
/// * `poly_a` - The number of polyA bases
/// * `gpa` - The number of genomic A's
/// * `line` - The original line of the read
///
/// # Example
///
/// ```rust, no_run
/// let line = "chr1\t1000\t2000\tread1_3Clip0_PolyA50_PolyARead52\t0\t+\t1000\t2000\t0\t0\t0\t0";
/// let data = PolyAPred::read(line, OverlapType::Exon, false).unwrap();
///
/// assert_eq!(data, PolyAPred {
///   name: "read1".to_string(),
///   chrom: "chr1".to_string(),
///   strand: Strand::Forward,
///   start: 1000,
///   end: 2000,
///   cds_start: 1000,
///   cds_end: 2000,
///   clip: 0,
///   clipped_a: 50,
///   poly_a: 52,
///   gpa: 2,
///   line: line.to_string(),
/// });
/// ```
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct PolyAPred {
    pub name: String,
    pub chrom: String,
    pub strand: Strand,
    pub start: u64,
    pub end: u64,
    pub cds_start: u64,
    pub cds_end: u64,
    pub clip: u32,
    pub clipped_a: u32,
    pub poly_a: u32,
    pub gpa: u32,
    pub line: String,
}

impl PolyAPred {
    #[inline(always)]
    pub fn read(line: &str, _: OverlapType, _: bool) -> Result<PolyAPred, &'static str> {
        if line.is_empty() {
            return Err("Empty line");
        }

        let get = |field: &str| field.parse::<u64>().map_err(|_| "Cannot parse field");

        let mut fields = line.split('\t');
        let (chrom, tx_start, tx_end, name, _, strand, cds_start, cds_end) = (
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
        );

        let (tx_start, tx_end, cds_start, cds_end) =
            abs_pos(tx_start, tx_end, cds_start, cds_end, strand, get)?;

        let strand = match strand {
            '+' => Strand::Forward,
            '-' => Strand::Reverse,
            _ => return Err("ERROR: Strand is not + or -"),
        };

        let (clip, clipped_a, poly_a, gpa) = get_polya_stats(name);

        Ok(PolyAPred {
            name: name.into(),
            chrom: chrom.into(),
            strand,
            start: tx_start,
            end: tx_end,
            cds_start,
            cds_end,
            clip,
            clipped_a,
            poly_a,
            gpa,
            line: line.into(),
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
        self.cds_start
    }

    fn cds_end(&self) -> u64 {
        self.cds_end
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

/// Extracts the polyA stats from the read name
///
/// # Arguments
///
/// * `read` - A string slice that holds the read name
///
/// # Returns
///
/// A tuple with the following values:
///
/// * `clip` - The number of clipped bases
/// * `clipped_a` - The number of clipped A's
/// * `read_a` - The number of A's in the read
/// * `gpa` - The number of genomic A's
///
/// # Example
///
/// ```
/// let read = "m54164U_210309_085211/74646562/ccs_PerID0.995_5Clip0_3Clip0_PolyA49_PolyARead50";
/// let stats = get_polya_stats(read);
///
/// assert_eq!(stats, (0, 49, 50, 1));
/// ```
fn get_polya_stats(read: &str) -> (u32, u32, u32, u32) {
    let tags = get_tags_from_read(read);

    let clip3 = *tags.get("TC").unwrap() as u32;
    let clipped_a = *tags.get("PA").unwrap() as u32;
    let read_a = *tags.get("PR").unwrap() as u32;

    // INFO: the whole polyA is clipped!
    if clipped_a == read_a {
        return (clip3, clipped_a, read_a, 0);
    } else if clipped_a > read_a {
        log::error!("ERROR: clipped A's is greater than read A's");
        std::process::exit(1);
    }

    let gpa = read_a - (clip3 + clipped_a);

    (clip3, clipped_a, read_a, gpa)
}

fn get_tags_from_read(read: &str) -> HashMap<String, usize> {
    let mut map = HashMap::with_capacity(5); // WARN: enforcing 5 tags -> fusion tags come after this step!

    if let Some((_, tags_part)) = read.split_once(BIG_SEP) {
        // WARN: enforcing two letter tag!
        for tag in tags_part.split(SEP) {
            if tag.len() >= 3 {
                // Safe slicing because ASCII: two-letter key + numeric value
                let key = &tag[..2];
                let val = &tag[2..];
                if let Ok(parsed_val) = val.parse::<usize>() {
                    map.insert(key.to_string(), parsed_val);
                }
            }
        }
    }

    map
}

impl From<GenePred> for PolyAPred {
    fn from(read: GenePred) -> Self {
        let (clip, clipped_a, poly_a, gpa) = get_polya_stats(&read.name);

        PolyAPred {
            name: read.name,
            chrom: read.chrom,
            strand: read.strand,
            start: read.start,
            end: read.end,
            cds_start: read.cds_start,
            cds_end: read.cds_end,
            clip,
            clipped_a,
            poly_a,
            gpa,
            line: read.line,
        }
    }
}

#[cfg(test)]
mod tests {
    use config::BedColumn;

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

    #[test]
    fn test_modify_line_color() {
        let mut gp = GenePred {
            name: "read1".to_string(),
            chrom: "chr1".to_string(),
            strand: Strand::Forward,
            start: 1000,
            end: 2000,
            cds_start: 1000,
            cds_end: 2000,
            exon_count: 1,
            exons: vec![(1000, 2000)],
            introns: vec![],
            line: "chr1\t1000\t2000\tread1\t0\t+\t1000\t2000\t0\t0\t0\t0".to_string(),
            is_ref: false,
        };

        gp.modify_field(BedColumn::ItemRgb.into(), "0,0,255");

        assert_eq!(
            gp.line,
            "chr1\t1000\t2000\tread1\t0\t+\t1000\t2000\t0,0,255\t0\t0\t0"
        );
    }

    #[test]
    fn test_modify_line_name() {
        let mut gp = GenePred {
            name: "read1".to_string(),
            chrom: "chr1".to_string(),
            strand: Strand::Forward,
            start: 1000,
            end: 2000,
            cds_start: 1000,
            cds_end: 2000,
            exon_count: 1,
            exons: vec![(1000, 2000)],
            introns: vec![],
            line: "chr1\t1000\t2000\tread1\t0\t+\t1000\t2000\t0\t0\t0\t0".to_string(),
            is_ref: false,
        };

        let name = format!("{}_{}", gp.name, "SNG");
        gp.modify_field(BedColumn::Name.into(), &name);

        assert_eq!(
            gp.line,
            "chr1\t1000\t2000\tread1_SNG\t0\t+\t1000\t2000\t0\t0\t0\t0"
        );
    }

    #[test]
    fn test_absolute_cds_mapping_forward() {
        let line = "chr6\t8259278\t8593709\tR441_chr6__FC23#TC31#PA0#PR0#IY1000\t60\t+\t8259278\t8593709\t43,118,219\t10\t215,109,62,215,152,87,117,153,211,2165\t0,11194,114627,167692,278538,299221,313903,320308,323301,332266";
        let mut gp = Bed12::read(line, OverlapType::Exon, false)
            .unwrap_or_else(|e| panic!("ERROR: could not parse line into GenePred: {e}"));

        let orf_start = 237;
        let orf_end = 420;

        let (predicted_cds_start, predicted_cds_end) =
            gp.map_absolute_cds(orf_start, orf_end).unwrap();

        assert_eq!(predicted_cds_start, 8270494);
        assert_eq!(predicted_cds_end, 8427004);
    }

    #[test]
    fn test_absolute_cds_mapping_reverse() {
        let line = "chr6\t8259278\t8593709\tR441_chr6__FC23#TC31#PA0#PR0#IY1000\t60\t-\t8259278\t8593709\t43,118,219\t10\t215,109,62,215,152,87,117,153,211,2165\t0,11194,114627,167692,278538,299221,313903,320308,323301,332266";
        let mut gp = Bed12::read(line, OverlapType::Exon, false)
            .unwrap_or_else(|e| panic!("ERROR: could not parse line into GenePred: {e}"));

        let orf_start = 237;
        let orf_end = 420;

        let (predicted_cds_start, predicted_cds_end) =
            gp.map_absolute_cds(orf_start, orf_end).unwrap();

        assert_eq!(predicted_cds_start, 8270494);
        assert_eq!(predicted_cds_end, 8427004);
    }
}
