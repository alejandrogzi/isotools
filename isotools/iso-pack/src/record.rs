use config::{BedParser, OverlapType, Sequence, Strand, SupportType, BIG_SEP, SCALE, SEP};
use hashbrown::{HashMap, HashSet};
use serde::{Deserialize, Serialize};

use std::collections::BTreeSet;

/// Represents a single genomic feature in the BED12 format, which contains
/// detailed information about a gene's structure.
///
/// This struct acts as a wrapper for the `GenePred` struct, simplifying its
/// representation to align with the BED12 standard.
#[derive(Debug, PartialEq, Clone)]
pub struct Bed12 {
    pub data: GenePred,
}

/// Represents a gene prediction record, containing comprehensive information
/// about a gene's genomic location, structure, and sequence.
///
/// This struct is a core data structure for storing gene-centric information
/// and is designed to be serialized and deserialized for data exchange.
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct GenePred {
    /// The name of the gene or transcript.
    pub name: String,
    /// The chromosome on which the gene is located.
    pub chrom: String,
    /// The strand of the chromosome (e.g., forward or reverse).
    pub strand: Strand,
    /// The start position of the gene on the chromosome (0-indexed).
    pub start: u64,
    /// The end position of the gene on the chromosome (1-indexed).
    pub end: u64,
    /// The start position of the coding sequence (CDS) on the chromosome.
    pub cds_start: u64,
    /// The end position of the coding sequence (CDS) on the chromosome.
    pub cds_end: u64,
    /// A vector of tuples, where each tuple represents the start and end
    /// coordinates of an exon.
    pub exons: Vec<(u64, u64)>,
    /// A vector of tuples, where each tuple represents the start and end
    /// coordinates of an intron.
    pub introns: Vec<(u64, u64)>,
    /// The total length of all exons combined.
    pub exon_len: u64,
    /// The number of exons in the gene.
    pub exon_count: usize,
    /// The original line from which this record was parsed.
    pub line: String,
    /// A boolean flag indicating if this is a reference gene prediction.
    pub is_ref: bool,
}

impl GenePred {
    /// Given a position within the concatenated exonic sequence, this function
    /// maps it back to its corresponding genomic coordinate.
    ///
    /// The function iterates through the gene's exons, calculating the length of
    /// each exon. It keeps a running total of the length of the exons processed so
    /// far (`current_pos`). If the target position (`pos`) falls within the current
    /// exon, it calculates the offset from the exon's start and returns the absolute
    /// genomic coordinate. If the target position is beyond all exons, it returns `None`.
    ///
    /// # Arguments
    ///
    /// * `pos` - A `u64` representing the position within the concatenated exonic sequence (0-indexed).
    ///
    /// # Returns
    ///
    /// An `Option<u64>` which is:
    /// * `Some(genomic_position)` if the position is found within the exons.
    /// * `None` if the position is outside the total length of all exons.
    ///
    /// # Panics
    ///
    /// This function does not panic under normal conditions.
    ///
    fn get_pos_in_exons(&self, pos: u64) -> Option<u64> {
        let mut current_pos = 0;

        for exon in self.exons.iter() {
            let block = exon.1 - exon.0;
            if current_pos + block < pos {
                current_pos += block
            } else {
                let remainder = pos - current_pos;
                return Some(exon.0 + remainder);
            }
        }

        None
    }

    /// Maps the start and end positions of an Open Reading Frame (ORF) from the
    /// concatenated exonic sequence to absolute genomic coordinates.
    ///
    /// This function is crucial for converting ORF coordinates, which are relative to the
    /// combined sequence of all exons, into the real genomic coordinates (chromosome, start,
    /// end). It calls `get_pos_in_exons` for both the ORF's start and end positions.
    ///
    /// The function also accounts for the gene's strand. For the forward strand, the
    /// converted coordinates are returned directly. For the reverse strand, the coordinates
    /// are flipped and adjusted based on a `SCALE` value to correctly reflect their position
    /// on the reverse complement.
    ///
    /// # Arguments
    ///
    /// * `orf_start` - A `u64` representing the starting position of the ORF within the
    ///   concatenated exonic sequence.
    /// * `orf_end` - A `u64` representing the ending position of the ORF within the
    ///   concatenated exonic sequence.
    ///
    /// # Returns
    ///
    /// A tuple `(u64, u64)` containing the absolute genomic start and end coordinates of the ORF.
    ///
    /// # Panics
    ///
    /// This function will panic if `get_pos_in_exons` returns `None` for either `orf_start`
    /// or `orf_end`, indicating that the ORF coordinates are outside the exonic sequence.
    ///
    pub fn get_cds_from_pos(&self, orf_start: u64, orf_end: u64) -> (u64, u64) {
        let cds_start = self.get_pos_in_exons(orf_start).unwrap_or_else(|| {
            panic!("ERROR: could not map ORF_START to exonic coordinates: {orf_start} -> {self:?}")
        });

        let cds_end = self.get_pos_in_exons(orf_end).unwrap_or_else(|| {
            panic!("ERROR: could not map ORF_END to exonic coordinates: {orf_end} -> {self:?}")
        });

        match self.strand {
            Strand::Forward => (cds_start, cds_end),
            Strand::Reverse => (SCALE - cds_end, SCALE - cds_start),
        }
    }

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
    /// * `(u64, u64)` - A tuple of genomic start and end coordinates corresponding
    ///   to the CDS region. Returns `(start, end)` if mapping succeeds; otherwise,
    ///   may return `Error` if the coordinates do not map within the exon structure.
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
    pub fn map_absolute_cds(&mut self, orf_start: u64, orf_end: u64) -> (u64, u64) {
        let (cds_start, cds_end) = self.get_cds_from_pos(orf_start, orf_end);

        assert!(
            cds_start > 0,
            "ERROR: {orf_start} is not bigger than 0! -> {self:?}"
        );
        assert!(
            cds_end > 0,
            "ERROR: {orf_end} is not bigget than 0! -> {self:?}"
        );

        (cds_start, cds_end)
    }
}

/// Represents a collection of gene prediction records associated with a single
/// reference gene.
///
/// This struct is used to group multiple `GenePred` instances that are derived from
/// a common reference gene, along with their structural features. It is
/// particularly useful for analyzing and comparing the predicted gene structures.
#[derive(Debug, PartialEq, Clone)]
pub struct RefGenePred {
    /// A vector of `GenePred` instances that are associated with the reference.
    pub reads: Vec<GenePred>,
    /// A sorted set of exon start and end coordinates, providing a quick way to check
    /// for overlaps and positions.
    pub starts: BTreeSet<(u64, u64)>,
    /// A sorted set of intron middle coordinates, used for efficient searching.
    pub middles: BTreeSet<(u64, u64)>,
    /// A sorted set of intron start and end coordinates.
    pub introns: BTreeSet<(u64, u64)>,
    /// The overall genomic bounds of the reference gene (start and end).
    pub bounds: (u64, u64),
    /// The strand of the reference gene.
    pub strand: Strand,
}

/// A container for intron-related data, used to group introns by chromosome and strand.
///
/// This struct acts as a bucket for all introns found on a specific chromosome and strand,
/// storing detailed statistics and metadata for each intron in a hash map.
#[derive(Debug, PartialEq, Clone)]
pub struct IntronBucket {
    /// The chromosome identifier.
    pub chrom: String,
    /// The strand (forward or reverse) of the introns in this bucket.
    pub strand: Strand,
    /// A `HashMap` where the key is a tuple representing the intron's
    /// start and end coordinates, and the value is `IntronPredStats`.
    pub introns: HashMap<(u64, u64), IntronPredStats>,
}

/// Represents a single predicted intron with its genomic location and associated statistics.
///
/// This struct combines the genomic coordinates of an intron with the
/// `IntronPredStats` struct, providing a complete record for a single intron.
#[derive(Debug, PartialEq, Clone)]
pub struct IntronPred {
    /// The chromosome identifier.
    pub chrom: String,
    /// The strand of the intron.
    pub strand: Strand,
    /// The start position of the intron.
    pub start: u64,
    /// The end position of the intron.
    pub end: u64,
    /// A struct containing various statistical and contextual data about the intron.
    pub stats: IntronPredStats,
}

/// Holds a variety of statistical and contextual data for a predicted intron.
///
/// This struct contains metrics from different prediction tools and classification
/// systems, used to evaluate the quality and nature of a predicted intron.
#[derive(Debug, PartialEq, Clone)]
pub struct IntronPredStats {
    /// The frequency of how many reads contain this intron.
    pub seen: usize,
    /// The frequency of how many reads span this intron.
    pub spanned: usize,
    /// SpliceAI score for the donor site.
    pub splice_ai_donor: f32,
    /// SpliceAI score for the acceptor site.
    pub splice_ai_acceptor: f32,
    /// MaxEntScan score for the donor site.
    pub max_ent_donor: f32,
    /// MaxEntScan score for the acceptor site.
    pub max_ent_acceptor: f32,
    /// The sequence around the donor site.
    pub donor_sequence: String,
    /// The sequence around the acceptor site.
    pub acceptor_sequence: String,
    /// The MaxEntScan 9-mer donor context sequence.
    pub donor_context: Sequence,
    /// The MaxEntScan 23-mer acceptor context sequence.
    pub acceptor_context: Sequence,
    /// The classification of the intron's position according to TOGA.
    pub intron_position: IntronPosition,
    /// A boolean indicating if the intron is supported by TOGA.
    pub is_toga_supported: bool,
    /// A boolean indicating if the intron maintains the reading frame.
    pub is_in_frame: bool,
    /// The RT-switch context sequence for the donor site.
    pub donor_rt_context: String,
    /// The RT-switch context sequence for the acceptor site.
    pub acceptor_rt_context: String,
    /// A boolean indicating if the intron is an RT-switch intron.
    pub is_rt_intron: bool,
    /// A boolean indicating if the intron is a TOGA-nag intron.
    pub is_nag_intron: bool,
    /// A classification of the intron's support type.
    pub support: SupportType,
}

/// BedParser implementation for IntronPred
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

/// Implements the `IntronPred` struct, providing methods for creating instances from
/// a string or a `GenePred` object.
impl IntronPred {
    /// Parses a tab-delimited line of text into a new `IntronPred` instance.
    ///
    /// This function is the primary method for constructing an `IntronPred` from a
    /// data file line. It splits the input line by tabs and attempts to parse each
    /// field into the corresponding type. It also handles strand-specific coordinate
    /// conversion for the reverse strand.
    ///
    /// # Arguments
    ///
    /// * `line` - A string slice (`&str`) representing a single line from an intron prediction
    ///   data file.
    ///
    /// # Returns
    ///
    /// A `Result<Self, Box<dyn std::error::Error>>` which is:
    /// * `Ok(IntronPred)` if the parsing is successful.
    /// * `Err(Box<dyn std::error::Error>)` if there is an error during parsing, such as
    ///   a malformed line or an invalid numeric value.
    ///
    /// # Panics
    ///
    /// The function will panic if it fails to parse any of the fields or if a field is
    /// missing. It also panics if the strand is not `+` or `-`.
    ///
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

/// Implements a conversion from `GenePred` to `IntronPred`.
impl From<GenePred> for IntronPred {
    /// Creates a new `IntronPred` instance from a `GenePred` object.
    ///
    /// This conversion method is useful when you have a `GenePred` record and need to
    /// create a corresponding `IntronPred` with its statistics. It extracts the
    /// relevant statistical data from the `GenePred`'s line and uses it to construct
    /// the `IntronPredStats` component.
    ///
    /// # Arguments
    ///
    /// * `read` - A `GenePred` instance.
    ///
    /// # Returns
    ///
    /// A new `IntronPred` instance.
    ///
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

/// Implements the `IntronPredStats` struct, providing methods for creation and formatting.
impl IntronPredStats {
    /// Creates a new `IntronPredStats` instance with all fields initialized to default values.
    ///
    /// This is a convenience constructor for creating a blank statistics object.
    ///
    /// # Returns
    ///
    /// A new `IntronPredStats` instance with default values.
    ///
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

    /// Creates a new `IntronPredStats` instance from a vector of string slices.
    ///
    /// This method is designed to parse a slice of `&str` fields from a data line and
    /// convert them into the appropriate types for the `IntronPredStats` struct.
    ///
    /// # Arguments
    ///
    /// * `data` - A `Vec<&str>` containing the string representations of the statistical fields.
    ///
    /// # Returns
    ///
    /// A new `IntronPredStats` instance.
    ///
    /// # Panics
    ///
    /// The function will panic if it fails to parse any of the fields into the
    /// correct type (e.g., parsing a non-numeric string as `f32`).
    ///
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

    /// Formats the `IntronPredStats` into a tab-separated string suitable for writing to a file.
    ///
    /// This method combines the statistical data with the intron's genomic coordinates
    /// and strand to produce a complete, formatted line.
    ///
    /// # Arguments
    ///
    /// * `chr` - The chromosome identifier.
    /// * `strand` - The strand of the intron.
    /// * `start` - The start position of the intron.
    /// * `end` - The end position of the intron.
    ///
    /// # Returns
    ///
    /// A `String` containing the tab-separated line of data.
    ///
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

/// Represents the positional context of an intron within a gene.
///
/// This enum is used to classify whether an intron is located in a UTR, a CDS,
/// spans both, or its position is unknown. It implements `Display` to provide
/// a user-friendly string representation.
#[derive(Debug, PartialEq, Clone, Eq, Hash)]
pub enum IntronPosition {
    /// Intron is located in an untranslated region.
    UTR,
    /// Intron is located within the coding sequence.
    CDS,
    /// Intron spans both a UTR and a CDS.
    Mixed,
    /// The position of the intron is unknown.
    Unknown,
}

/// Implements `std::fmt::Display` for `IntronPosition`, allowing it to be
/// easily formatted as a string.
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

/// Implements various utility methods for accessing and manipulating data within the `GenePred` struct.
impl GenePred {
    /// Returns a shared reference to the original line of text from which the `GenePred` was parsed.
    ///
    /// # Arguments
    ///
    /// * `self` - A reference to the `GenePred` instance.
    ///
    /// # Returns
    ///
    /// A `&String` containing the original line of data.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let gene_pred: GenePred = // ... create a GenePred instance ...;
    /// let line = gene_pred.line();
    /// println!("Original line: {}", line);
    /// ```
    pub fn line(&self) -> &String {
        &self.line
    }

    /// Returns a shared reference to the name of the gene.
    ///
    /// # Arguments
    ///
    /// * `self` - A reference to the `GenePred` instance.
    ///
    /// # Returns
    ///
    /// A `&String` containing the gene's name.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let gene_pred: GenePred = // ... create a GenePred instance ...;
    /// let name = gene_pred.name();
    /// println!("Gene name: {}", name);
    /// ```
    pub fn name(&self) -> &String {
        &self.name
    }

    /// Returns a boolean indicating whether the `GenePred` instance represents a reference gene.
    ///
    /// # Arguments
    ///
    /// * `self` - A reference to the `GenePred` instance.
    ///
    /// # Returns
    ///
    /// A `bool` which is `true` if the gene is a reference, `false` otherwise.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let gene_pred: GenePred = // ... create a GenePred instance ...;
    /// if gene_pred.is_ref() {
    ///     println!("This is a reference gene.");
    /// }
    /// ```
    pub fn is_ref(&self) -> bool {
        self.is_ref
    }

    /// Returns a mutable reference to the original line of text.
    ///
    /// This method allows for in-place modification of the raw data line.
    ///
    /// # Arguments
    ///
    /// * `self` - A mutable reference to the `GenePred` instance.
    ///
    /// # Returns
    ///
    /// A `&mut String` that can be used to modify the original line.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let mut gene_pred: GenePred = // ... create a GenePred instance ...;
    /// *gene_pred.line_mut() = "a new line of data".to_string();
    /// ```
    pub fn line_mut(&mut self) -> &mut String {
        &mut self.line
    }

    /// Returns a clone of the coordinates of the first exon.
    ///
    /// This method provides quick access to the start and end coordinates of the initial exon.
    ///
    /// # Arguments
    ///
    /// * `self` - A reference to the `GenePred` instance.
    ///
    /// # Returns
    ///
    /// A tuple `(u64, u64)` containing the start and end coordinates of the first exon.
    ///
    /// # Panics
    ///
    /// This function will panic if the `exons` vector is empty.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let gene_pred: GenePred = // ... create a GenePred instance with exons ...;
    /// let first_exon = gene_pred.get_first_exon();
    /// println!("First exon: {:?}", first_exon);
    /// ```
    #[inline(always)]
    pub fn get_first_exon(&self) -> (u64, u64) {
        self.exons[0].clone()
    }

    /// Returns a sorted set of all exon coordinates except for the first one.
    ///
    /// This is useful for analyzing the internal structure of multi-exon genes.
    ///
    /// # Arguments
    ///
    /// * `self` - A reference to the `GenePred` instance.
    ///
    /// # Returns
    ///
    /// A `BTreeSet<(u64, u64)>` containing the coordinates of the middle exons.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let gene_pred: GenePred = // ... create a GenePred instance with multiple exons ...;
    /// let middle_exons = gene_pred.get_middle_exons();
    /// println!("Middle exons: {:?}", middle_exons);
    /// ```
    #[inline(always)]
    pub fn get_middle_exons(&self) -> BTreeSet<(u64, u64)> {
        self.exons[1..].iter().cloned().collect()
    }

    /// Returns a sorted set of all exon coordinates.
    ///
    /// # Arguments
    ///
    /// * `self` - A reference to the `GenePred` instance.
    ///
    /// # Returns
    ///
    /// A `BTreeSet<(u64, u64)>` containing the coordinates of all exons.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let gene_pred: GenePred = // ... create a GenePred instance ...;
    /// let all_exons = gene_pred.get_exons();
    /// println!("All exons: {:?}", all_exons);
    /// ```
    #[inline(always)]
    pub fn get_exons(&self) -> BTreeSet<(u64, u64)> {
        self.exons.iter().cloned().collect()
    }

    /// Returns a sorted set of all intron coordinates.
    ///
    /// # Arguments
    ///
    /// * `self` - A reference to the `GenePred` instance.
    ///
    /// # Returns
    ///
    /// A `BTreeSet<(u64, u64)>` containing the coordinates of all introns.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let gene_pred: GenePred = // ... create a GenePred instance ...;
    /// let all_introns = gene_pred.get_introns();
    /// println!("All introns: {:?}", all_introns);
    /// ```
    #[inline(always)]
    pub fn get_introns(&self) -> BTreeSet<(u64, u64)> {
        self.introns.iter().cloned().collect()
    }

    /// Returns the genomic coordinates of the 5' untranslated region (UTR).
    ///
    /// # Arguments
    ///
    /// * `self` - A reference to the `GenePred` instance.
    ///
    /// # Returns
    ///
    /// A tuple `(u64, u64)` representing the start and end coordinates of the 5' UTR.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let gene_pred: GenePred = // ... create a GenePred instance ...;
    /// let five_utr = gene_pred.get_five_utr();
    /// println!("5' UTR: {:?}", five_utr);
    /// ```
    pub fn get_five_utr(&self) -> (u64, u64) {
        (self.start, self.cds_start)
    }

    /// Returns the genomic coordinates of the 3' untranslated region (UTR).
    ///
    /// # Arguments
    ///
    /// * `self` - A reference to the `GenePred` instance.
    ///
    /// # Returns
    ///
    /// A tuple `(u64, u64)` representing the start and end coordinates of the 3' UTR.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let gene_pred: GenePred = // ... create a GenePred instance ...;
    /// let three_utr = gene_pred.get_three_utr();
    /// println!("3' UTR: {:?}", three_utr);
    /// ```
    pub fn get_three_utr(&self) -> (u64, u64) {
        (self.cds_end, self.end)
    }

    /// Splits the gene's name and returns the second part.
    ///
    /// This method assumes the name has a specific format with at least two
    /// parts separated by a ".".
    ///
    /// # Arguments
    ///
    /// * `self` - A reference to the `GenePred` instance.
    ///
    /// # Returns
    ///
    /// A `&str` containing the second part of the split name.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let gene_pred_name = "ENST00000456328.2";
    /// // let gene_pred = ...;
    /// // let split_name = gene_pred.get_split_name();
    /// // assert_eq!(split_name, "2");
    /// ```
    pub fn get_split_name(&self) -> &str {
        let id: Vec<&str> = self.name.splitn(3, ".").collect();
        id[1]
    }

    /// Creates a new `GenePred` instance with a modified color field in the original line.
    ///
    /// This method is useful for updating the visual representation of a gene in a
    /// format like BED12 without modifying the other genomic data fields. It clones
    /// all fields and then updates the `line` string.
    ///
    /// # Arguments
    ///
    /// * `self` - A reference to the `GenePred` instance.
    /// * `color` - A string slice (`&str`) representing the new color to be used.
    ///
    /// # Returns
    ///
    /// A new `GenePred` instance with the updated color line.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let gene_pred: GenePred = // ... create a GenePred instance ...;
    /// let colored_gene = gene_pred.colorline("255,0,0"); // Red color
    /// println!("New colored line: {}", colored_gene.line());
    /// ```
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
            exon_len: self.exon_len,
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

    /// Modifies the name field (the fourth tab-separated field, at index 3) of the
    /// gene prediction's line and returns the new line.
    ///
    /// This method is useful for generating a new output line with an updated name
    /// without modifying the `name` field within the struct itself. It clones the original
    /// line, updates the name, and returns the modified string.
    ///
    /// # Arguments
    ///
    /// * `self` - A mutable reference to the `GenePred` instance.
    /// * `name` - A string slice (`&str`) containing the new name.
    ///
    /// # Returns
    ///
    /// A `String` containing the newly formatted line with the updated name.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let mut gene_pred = GenePred {
    ///     // ... fields with a name "old_name"
    ///     line: "chr1\t100\t200\told_name\t...".to_string(),
    ///     // ... other fields
    /// };
    /// let new_line = gene_pred.mut_name_from_line("new_name");
    /// // The `gene_pred.line` remains unchanged, but `new_line` is now "chr1\t100\t200\tnew_name\t...".
    /// ```
    pub fn mut_name_from_line(&mut self, name: &str) -> String {
        let line = self.line.clone();
        let mut fields = line.split('\t').collect::<Vec<_>>();
        fields[3] = name;

        fields.join("\t")
    }

    /// Inserts new start and end coordinates for the Coding Sequence (CDS).
    ///
    /// This method is used to define the boundaries of the CDS, effectively "inserting"
    /// UTRs (Untranslated Regions) by setting the `cds_start` and `cds_end` fields.
    ///
    /// **Warning:** This method directly modifies the `cds_start` and `cds_end` fields of
    /// the `GenePred` struct, overwriting any previous values.
    ///
    /// # Arguments
    ///
    /// * `self` - A mutable reference to the `GenePred` instance.
    /// * `start` - A `u64` representing the new CDS start position.
    /// * `end` - A `u64` representing the new CDS end position.
    ///
    /// # Returns
    ///
    /// A mutable reference to the `GenePred` instance (`&mut Self`), allowing for method chaining.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let mut gene_pred: GenePred = // ... create a GenePred instance ...;
    /// // Define the new CDS start and end
    /// gene_pred.insert_utr(100, 500);
    ///
    /// assert_eq!(gene_pred.cds_start, 100);
    /// assert_eq!(gene_pred.cds_end, 500);
    /// ```
    pub fn insert_utr(&mut self, start: u64, end: u64) -> &mut Self {
        // WARN: modifies CDStart and CDEnd
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

/// BedParser implementation for GenePred
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

/// Implements a conversion from `IntronPred` to `GenePred`.
impl From<IntronPred> for GenePred {
    /// Creates a new `GenePred` instance from an `IntronPred` object.
    ///
    /// This conversion is used to represent an intron as if it were a single-exon gene,
    /// a common practice for tools that primarily work with gene structures.
    ///
    /// **Warning:** The `exons` field is initialized with a simulated exon that is 2 base pairs
    /// larger than the intron's bounds on both sides (`intron.start - 2` to `intron.end + 2`).
    /// The `exon_len` and `exon_count` are also not correctly calculated. This is a deliberate
    /// simplification for a specific use case and may not be suitable for general-purpose
    /// gene prediction.
    ///
    /// # Arguments
    ///
    /// * `intron` - The `IntronPred` object to convert.
    ///
    /// # Returns
    ///
    /// A new `GenePred` instance.
    ///
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
            exon_len: 0, // WARN: unimplemented! -> should be the diff exon_s and exon_e
            exon_count: 0,
            line: intron
                .stats
                .fmt(&intron.chrom, &intron.strand, intron.start, intron.end), // INFO: holds IntronPredStats
            is_ref: true, // INFO: always true since it is a reference file being read
        }
    }
}

/// Implements a conversion from `Bed6` to `GenePred`.
impl From<Bed6> for GenePred {
    /// Creates a new `GenePred` instance from a `Bed6` record.
    ///
    /// This conversion is used to parse a simple BED6 record, which lacks detailed exon
    /// information, into a `GenePred` struct. It treats the entire BED6 record as a
    /// single exon.
    ///
    /// # Arguments
    ///
    /// * `record` - The `Bed6` object to convert.
    ///
    /// # Returns
    ///
    /// A new `GenePred` instance.
    ///
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
            exon_len: record.coord.1 - record.coord.0,
            exon_count: 0,
            line: "".to_string(),
            is_ref: false,
        }
    }
}

/// Implements methods for the `Bed12` struct, primarily for parsing.
impl Bed12 {
    /// Parses a line from a BED12-formatted file into a `GenePred` struct.
    ///
    /// This is a complex parsing function that handles the various fields of a BED12
    /// record, including coordinate parsing, exon/intron reconstruction, and strand
    /// handling. It relies on helper functions `abs_pos` and `get_coords` to
    /// perform the coordinate calculations.
    ///
    /// # Arguments
    ///
    /// * `line` - A string slice (`&str`) representing a single line from a BED12 file.
    /// * `overlap` - An `OverlapType` enum indicating the type of overlap to consider.
    /// * `is_ref` - A boolean flag indicating if the record is a reference gene.
    ///
    /// # Returns
    ///
    /// A `Result<GenePred, &'static str>` which is:
    /// * `Ok(GenePred)` if the parsing is successful.
    /// * `Err(&'static str)` if the line is empty or a parsing error occurs.
    ///
    /// # Panics
    ///
    /// This function does not panic. Instead, it returns an error string.
    ///
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

        let (exons, introns, total_exon_len) = get_coords(
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
            exon_len: total_exon_len,
            exon_count,
            line: line.to_string(),
            is_ref,
        })
    }
}

/// BedParser implemention for Bed12
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

/// Parses the exon and intron coordinates from a BED12 record's `exon_starts` and `exon_sizes` fields.
///
/// This function calculates the absolute genomic coordinates of exons and introns, taking into
/// account the gene's strand and a specified overlap type. It is a core component of the
/// `Bed12::read` method.
///
/// # Arguments
///
/// * `starts` - A string slice containing a comma-separated list of relative exon start positions.
/// * `sizes` - A string slice containing a comma-separated list of exon sizes.
/// * `tx_start` - The transcription start position of the gene.
/// * `tx_end` - The transcription end position of the gene.
/// * `cds_start` - The coding sequence start position.
/// * `cds_end` - The coding sequence end position.
/// * `strand` - A character representing the gene's strand ('+' or '-').
/// * `overlap` - An `OverlapType` enum value that defines how to handle overlapping regions.
///
/// # Returns
///
/// A `Result<(HashSet<(u64, u64)>, HashSet<(u64, u64)>, u64), &'static str>`.
/// If successful, it returns a tuple containing:
/// * A `HashSet` of exon coordinates.
/// * A `HashSet` of intron coordinates.
/// * A `u64` representing the total length of all exons.
///
/// Returns an `Err` if parsing fails, if the number of starts and sizes don't match, or if the strand is invalid.
///
/// # Panics
///
/// This function will panic if any of the numeric fields cannot be parsed into a `u64`.
///
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
) -> Result<(HashSet<(u64, u64)>, HashSet<(u64, u64)>, u64), &'static str> {
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

    let total_exonic_len: u64 = sz.iter().sum();

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

    Ok((exons, introns, total_exonic_len))
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

/// Implements methods for the `RefGenePred` struct.
impl RefGenePred {
    /// Creates a new `RefGenePred` instance.
    ///
    /// This constructor is used to manually build a `RefGenePred` from its constituent parts.
    ///
    /// # Arguments
    ///
    /// * `reads` - A `Vec<GenePred>` containing the gene predictions.
    /// * `starts` - A `BTreeSet<(u64, u64)>` of start exon coordinates.
    /// * `middles` - A `BTreeSet<(u64, u64)>` of middle exon coordinates.
    /// * `introns` - A `BTreeSet<(u64, u64)>` of intron coordinates.
    /// * `bounds` - A tuple `(u64, u64)` for the genomic start and end of the gene.
    /// * `strand` - The `Strand` of the gene.
    ///
    /// # Returns
    ///
    /// A new `RefGenePred` instance.
    ///
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

    /// Gathers and returns a unique set of full gene names.
    ///
    /// # Arguments
    ///
    /// * `self` - A reference to the `RefGenePred` instance.
    ///
    /// # Returns
    ///
    /// A `HashSet<String>` containing the unique names of the genes.
    ///
    pub fn get_names(&self) -> HashSet<String> {
        let mut names = HashSet::new();
        self.reads.iter().for_each(|read| {
            names.insert(read.name.clone());
        });

        names
    }

    /// Gathers and returns a unique set of gene names after splitting them by ".".
    ///
    /// This method is useful for extracting a common identifier or transcript version from
    /// gene names that follow a specific `name.version.sub-version` pattern.
    ///
    /// # Arguments
    ///
    /// * `self` - A reference to the `RefGenePred` instance.
    ///
    /// # Returns
    ///
    /// A `HashSet<&str>` containing the unique split names.
    ///
    pub fn get_names_split(&self) -> HashSet<&str> {
        let mut names = HashSet::new();
        self.reads.iter().for_each(|read| {
            let id: Vec<&str> = read.name.splitn(3, ".").collect();
            names.insert(id[1]);
        });

        names
    }

    /// Merges the unique split names into a single, period-separated string.
    ///
    /// # Arguments
    ///
    /// * `self` - A reference to the `RefGenePred` instance.
    ///
    /// # Returns
    ///
    /// A `String` with the unique names joined by ".".
    ///
    pub fn merge_names(&self) -> String {
        let names = self.get_names_split().into_iter().collect::<Vec<&str>>();
        names.join(".")
    }

    /// Groups and merges all exons by their split name.
    ///
    /// This method aggregates all exons for each unique gene name (e.g., all exons from
    /// "gene1.1" and "gene1.2" are combined under "gene1").
    ///
    /// # Arguments
    ///
    /// * `self` - A reference to the `RefGenePred` instance.
    ///
    /// # Returns
    ///
    /// A `Vec<Vec<(u64, u64)>>` where each inner `Vec` contains the sorted, unique exon coordinates
    /// for a single gene (identified by its split name).
    ///
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

    /// Groups and merges all introns by their split name.
    ///
    /// **Warning:** This method returns a `HashSet` as the inner collection, which does not
    /// guarantee order.
    ///
    /// # Arguments
    ///
    /// * `self` - A reference to the `RefGenePred` instance.
    ///
    /// # Returns
    ///
    /// A `Vec<HashSet<(u64, u64)>>` where each `HashSet` contains the unique intron coordinates
    /// for a single gene (identified by its split name).
    pub fn smash_introns_by_name(&self) -> Vec<HashSet<(u64, u64)>> {
        // WARN: returns HashSet instead of Vec in inner collection!
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

    /// Groups and merges all exons by their split name, returning a `HashMap`.
    ///
    /// This is an alternative to `smash_exons_by_name` that provides a direct mapping from
    /// the gene's split name to its aggregated exons.
    ///
    /// # Arguments
    ///
    /// * `self` - A reference to the `RefGenePred` instance.
    ///
    /// # Returns
    ///
    /// A `HashMap<&str, Vec<(u64, u64)>>` where keys are the split names and values are
    /// the sorted exon coordinates.
    ///
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

    /// Creates a new `RefGenePred` instance from a `Vec<GenePred>`.
    ///
    /// This is the primary constructor, which automatically processes a collection of
    /// `GenePred` reads to build the `RefGenePred` fields, including bounds, exons,
    /// and introns. It also handles edge cases for empty input and strand determination.
    ///
    /// # Arguments
    ///
    /// * `reads` - A `Vec<GenePred>` to be used for building the `RefGenePred`.
    ///
    /// # Returns
    ///
    /// A new `RefGenePred` instance.
    ///
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

/// Implements methods for the `IntronBucket` struct.
impl IntronBucket {
    /// Creates a new `IntronBucket` from a collection of `GenePred` reads and a TOGA reference.
    ///
    /// This is a complex function that populates the `IntronBucket` with intron data from a set
    /// of reads, and then cross-references this data with a provided set of TOGA gene predictions.
    /// It determines various statistics and properties for each intron, such as whether it
    /// is supported by the TOGA reference and its position relative to the CDS.
    ///
    /// # Arguments
    ///
    /// * `reads` - A `Vec<GenePred>` representing the sequencing reads.
    /// * `toga` - A `Vec<GenePred>` representing the TOGA reference genes.
    ///
    /// # Returns
    ///
    /// A new `IntronBucket` instance.
    ///
    /// # Side Effects
    ///
    /// This method modifies the `toga_introns` set by removing introns as they are matched
    /// to introns from the `reads` vector.
    ///
    #[inline(always)]
    #[allow(unused_assignments)]
    pub fn from(reads: Vec<GenePred>, toga: Vec<GenePred>) -> Self {
        // INFO: Vec<GenePred> -> IntronBucket

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
    /// Creates a new `Bed4` instance by parsing a BED4-formatted string.
    ///
    /// This constructor parses a tab-separated line containing chromosome, start coordinate,
    /// end coordinate, and an ID. It converts the coordinates from the 0-based, half-open
    /// standard of the BED format to a 1-based, inclusive system by adjusting the start
    /// and end positions.
    ///
    /// # Arguments
    ///
    /// * `line` - A `String` containing the BED4 record.
    ///
    /// # Returns
    ///
    /// A `Result<Bed4, Box<dyn std::error::Error>>` which is `Ok(Bed4)` on success, or
    /// `Err` if the line is empty or fields cannot be parsed.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let line = "chr1\t100\t200\tgeneA".to_string();
    /// let bed_record = Bed4::new(line).unwrap();
    /// assert_eq!(bed_record.chrom, "chr1");
    /// assert_eq!(bed_record.coord, (101, 199));
    /// ```
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

    /// Creates a new `Bed4` instance from its component parts.
    ///
    /// This is a simple constructor that directly builds a `Bed4` struct without any
    /// coordinate conversion or parsing from a string. It is useful when the data is
    /// already in the correct format.
    ///
    /// # Arguments
    ///
    /// * `chrom` - A `String` for the chromosome name.
    /// * `start` - A `u64` for the start coordinate.
    /// * `end` - A `u64` for the end coordinate.
    /// * `id` - A `String` for the record ID.
    ///
    /// # Returns
    ///
    /// A new `Bed4` instance.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let bed_record = Bed4::from("chrX".to_string(), 500, 600, "geneB".to_string());
    /// assert_eq!(bed_record.coord, (500, 600));
    /// ```
    pub fn from(chrom: String, start: u64, end: u64, id: String) -> Bed4 {
        Bed4 {
            chrom,
            coord: (start, end),
            id,
        }
    }

    /// Appends the `Bed4` record to a mutable string in a simplified BED format.
    ///
    /// This method formats the `Bed4` record as `chrom\tstart\tend\n` and appends it to a
    /// provided string accumulator. The record's ID is not included in the output.
    ///
    /// # Arguments
    ///
    /// * `self` - A reference to the `Bed4` instance.
    /// * `acc` - A mutable reference to a `String` to which the output will be appended.
    ///
    /// # Example
    ///
    /// ```rust, no_run
    /// let bed_record = Bed4::from("chr3".to_string(), 10, 20, "region_1".to_string());
    /// let mut output = String::new();
    /// bed_record.send(&mut output);
    /// assert_eq!(output, "chr3\t10\t20\n");
    /// ```
    pub fn send(&self, acc: &mut String) {
        acc.push_str(&format!(
            "{}\t{}\t{}\n",
            self.chrom, self.coord.0, self.coord.1
        ));
    }
}

/// BedParser implementation for Bed4
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

/// Implements methods for the `Bed6` struct.
impl Bed6 {
    /// Creates a new `Bed6` instance by parsing a BED6-formatted string.
    ///
    /// This constructor parses a tab-separated line containing chromosome, start, end, ID, score,
    /// and strand. It handles coordinate conversion for reverse strand records, assuming a `SCALE`
    /// constant is available.
    ///
    /// # Arguments
    ///
    /// * `line` - A `String` containing the BED6 record.
    ///
    /// # Returns
    ///
    /// A `Result<Bed6, Box<dyn std::error::Error>>` which is `Ok(Bed6)` on success, or
    /// `Err` if the line is empty or fields cannot be parsed.
    ///
    /// # Panics
    ///
    /// The `SCALE` constant is not defined here; it must be a globally accessible `u64`
    /// representing the size of the genome or chromosome.
    ///
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

    /// Creates a new `Bed6` instance from its component parts.
    ///
    /// This is a simple constructor that builds a `Bed6` struct from its individual fields,
    /// setting the score to `0.0` by default.
    ///
    /// # Arguments
    ///
    /// * `chrom` - A `String` for the chromosome name.
    /// * `start` - A `u64` for the start coordinate.
    /// * `end` - A `u64` for the end coordinate.
    /// * `id` - A `String` for the record ID.
    /// * `strand` - A `Strand` enum value.
    ///
    /// # Returns
    ///
    /// A new `Bed6` instance.
    ///
    pub fn from(chrom: String, start: u64, end: u64, id: String, strand: Strand) -> Bed6 {
        Bed6 {
            chrom,
            coord: (start, end),
            id,
            strand,
            score: 0.0,
        }
    }

    /// Appends the `Bed6` record to a mutable string in a simplified BED format.
    ///
    /// This method formats the `Bed6` record as `chrom\tstart\tend\n` and appends it to a
    /// provided string accumulator. The record's ID, score, and strand are not included.
    ///
    /// # Arguments
    ///
    /// * `self` - A reference to the `Bed6` instance.
    /// * `acc` - A mutable reference to a `String` to which the output will be appended.
    ///
    pub fn send(&self, acc: &mut String) {
        acc.push_str(&format!(
            "{}\t{}\t{}\n",
            self.chrom, self.coord.0, self.coord.1
        ));
    }
}

/// BedParser implementation for Bed6
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

/// Implements methods for the `PolyAPred` struct.
impl PolyAPred {
    /// Parses a line of a file into a new `PolyAPred` instance.
    ///
    /// This function reads a tab-separated line, extracts fields like chromosome,
    /// coordinates, name, and strand, and then converts these into a `PolyAPred` struct.
    /// It handles coordinate and strand parsing, and extracts additional poly(A)
    /// statistics from the record's name field.
    ///
    /// # Arguments
    ///
    /// * `line` - The line to be parsed.
    /// * `_` - An `OverlapType` (unused).
    /// * `_` - A boolean (unused).
    ///
    /// # Returns
    ///
    /// A `Result<PolyAPred, &'static str>`, returning `Ok(PolyAPred)` on success
    /// or `Err` if parsing fails.
    ///
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

/// BedParser implementation for PolyAPred
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

/// Extracts two-letter key, usize value tags from a read's name field.
///
/// This function is designed to parse a specific tag format, where tags are separated by a
/// `BIG_SEP` and then a `SEP` character, with each tag being two characters followed by a
/// numeric value.
///
/// # Arguments
///
/// * `read` - A string slice containing the full read name, including tags.
///
/// # Returns
///
/// A `HashMap<String, usize>` containing the parsed tags.
///
/// # Panics
///
/// This function assumes the presence of `BIG_SEP` and `SEP` constants and that tag keys are
/// exactly two characters long. The capacity is enforced at 5.
///
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
    /// Converts a `GenePred` struct into a `PolyAPred` struct.
    ///
    /// This `From` implementation takes a `GenePred` instance and re-uses its fields to
    /// construct a `PolyAPred` instance. It also calls `get_polya_stats` to populate the
    /// additional poly(A) specific fields.
    ///
    /// # Arguments
    ///
    /// * `read` - The `GenePred` instance to be converted.
    ///
    /// # Returns
    ///
    /// A new `PolyAPred` instance.
    ///
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

        let (exons, introns, _) = get_coords(
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

        let (exons, introns, _) = get_coords(
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

        let (exons, introns, _) = get_coords(
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

        let (exons, introns, _) = get_coords(
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
            exon_len: 0,
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
            exon_len: 0,
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

        let (predicted_cds_start, predicted_cds_end) = gp.map_absolute_cds(orf_start, orf_end);

        assert_eq!(predicted_cds_start, 8270494);
        assert_eq!(predicted_cds_end, 8427004);
    }

    #[test]
    fn test_absolute_cds_mapping_reverse() {
        let line = "chr7\t45482351\t45520392\tR206671_chr7__FC0#TC23#PA0#PR0#IY998\t60\t-\t45482351\t45520392\t43,118,219\t13\t1162,233,188,161,232,126,154,171,210,117,620,482,1256\t0,1321,1711,20275,21300,21646,23812,24550,24947,25522,29018,33186,36785";
        let mut gp = Bed12::read(line, OverlapType::Exon, false)
            .unwrap_or_else(|e| panic!("ERROR: could not parse line into GenePred: {e}"));

        let orf_start = 1284;
        let orf_end = 3321;

        let (predicted_cds_start, predicted_cds_end) = gp.map_absolute_cds(orf_start, orf_end);
        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start, 45503698);
        assert_eq!(predicted_cds_end, 45515991);
    }

    #[test]
    fn test_get_cds_from_pos_forward() {
        let line = "chr12\t102521030\t102531285\tR909465_chr12__FC30#TC0#PA164#PR164#IY896\t60\t+\t102521030\t102531285\t43,118,219\t8\t399,47,94,69,159,417,476,417\t0,1113,3694,4517,6573,6895,7826,9838";
        let gp = Bed12::read(line, OverlapType::Exon, false)
            .unwrap_or_else(|e| panic!("ERROR: could not parse line into GenePred: {e}"));

        let orf_start = 1530;
        let orf_end = 1674;

        let (predicted_cds_start, predicted_cds_end) = gp.get_cds_from_pos(orf_start, orf_end);

        assert_eq!(predicted_cds_start, 102529201);
        assert_eq!(predicted_cds_end, 102530881);
    }

    #[test]
    fn test_get_pos_in_exons_forward() {
        let line = "chr12\t102521030\t102531285\tR909465_chr12__FC30#TC0#PA164#PR164#IY896\t60\t+\t102521030\t102531285\t43,118,219\t8\t399,47,94,69,159,417,476,417\t0,1113,3694,4517,6573,6895,7826,9838";
        let gp = Bed12::read(line, OverlapType::Exon, false)
            .unwrap_or_else(|e| panic!("ERROR: could not parse line into GenePred: {e}"));

        let orf_start = 1530;
        let orf_end = 1674;

        let predicted_cds_start = gp.get_pos_in_exons(orf_start);
        let predicted_cds_end = gp.get_pos_in_exons(orf_end);

        assert_eq!(predicted_cds_start, Some(102529201));
        assert_eq!(predicted_cds_end, Some(102530881));
    }

    #[test]
    fn test_get_cds_from_pos_reverse() {
        let line = "chr7\t45482351\t45520392\tR206671_chr7__FC0#TC23#PA0#PR0#IY998\t60\t-\t45482351\t45520392\t43,118,219\t13\t1162,233,188,161,232,126,154,171,210,117,620,482,1256\t0,1321,1711,20275,21300,21646,23812,24550,24947,25522,29018,33186,36785";
        let gp = Bed12::read(line, OverlapType::Exon, false)
            .unwrap_or_else(|e| panic!("ERROR: could not parse line into GenePred: {e}"));

        let orf_start = 358;
        let orf_end = 553;

        let (predicted_cds_start, predicted_cds_end) = gp.get_cds_from_pos(orf_start, orf_end);

        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start, 45519839);
        assert_eq!(predicted_cds_end, 45520034);
    }

    #[test]
    fn test_get_pos_in_exons_reverse() {
        let line = "chr7\t45482351\t45520392\tR206671_chr7__FC0#TC23#PA0#PR0#IY998\t60\t-\t45482351\t45520392\t43,118,219\t13\t1162,233,188,161,232,126,154,171,210,117,620,482,1256\t0,1321,1711,20275,21300,21646,23812,24550,24947,25522,29018,33186,36785";
        let gp = Bed12::read(line, OverlapType::Exon, false)
            .unwrap_or_else(|e| panic!("ERROR: could not parse line into GenePred: {e}"));

        let orf_start = 1284;
        let orf_end = 3321;

        let predicted_cds_end = gp.get_pos_in_exons(orf_start);
        let predicted_cds_start = gp.get_pos_in_exons(orf_end);

        dbg!(
            SCALE - predicted_cds_start.unwrap(),
            SCALE - predicted_cds_end.unwrap()
        );

        assert_eq!(SCALE - predicted_cds_start.unwrap(), 45503698);
        assert_eq!(SCALE - predicted_cds_end.unwrap(), 45515991);
    }
}
