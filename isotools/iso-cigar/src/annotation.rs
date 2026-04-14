// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Core modul for rescuing missed 3' splice junctions by CIGAR matching
//! Alejandro Gonzales-Irribarren, 2026
//!
//! This module contains the main function for rescuing missed 3' splice
//! junctions by CIGAR matching. It is a wrapper around the `engine` module
//! that handles the command-line interface and logging. In short, it reads
//! a BAM file and a GTF file, and writes a new BAM file with rescued junctions.
//!
//! A minimap2-aligned transcript can end exactly at an exon boundary
//! but carry a small 3' soft-clipped sequence that was not placed anywhere.
//! The hypothesis is: that clipped sequence is the start of the next downstream
//! exon, but the aligner missed the intron.
//!
//! In short, iso-cigar identifies candidate reads: those ending ±wiggle bp from a
//! known internal exon boundary, with soft-clip ≥ cutoff; retrieves the sequence of
//! the immediately downstream exon from the reference genome; checks whether the
//! soft-clipped bases match the beginning of that downstream exon; if yes: rewrites
//! the CIGAR to add an N intron gap and converts the soft-clip into a = match

use std::collections::HashMap;
use std::path::Path;

use anyhow::{anyhow, bail, Context, Result};
use genepred::{Bed12, GenePred, Gff, Gtf, Reader, Strand as GenePredStrand};

/// Strand orientation.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum Strand {
    Forward,
    Reverse,
}

/// An exon coordinate range.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct Exon {
    pub start: u64,
    pub end: u64,
}

/// A junction between two adjacent exons.
#[derive(Clone, Debug)]
pub struct Junction {
    pub id: usize,
    pub contig_id: u32,
    pub strand: Strand,
    pub upstream: Exon,
    pub downstream: Exon,
    pub boundary: u64,
    pub intron_len: u64,
}

/// Key for looking up junctions by genomic boundary.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct BoundaryKey {
    pub contig_id: u32,
    pub strand: Strand,
    pub boundary: u64,
}

/// Internal geometry key for junction deduplication.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
struct JunctionGeometry {
    contig_id: u32,
    strand: Strand,
    upstream: Exon,
    downstream: Exon,
    boundary: u64,
}

/// An index of genomic annotations for junction lookup.
#[derive(Clone, Debug)]
pub struct AnnotationIndex {
    contigs: Vec<Vec<u8>>,
    contig_ids: HashMap<Vec<u8>, u32>,
    junctions: Vec<Junction>,
    boundaries: HashMap<BoundaryKey, Vec<usize>>,
}

impl AnnotationIndex {
    /// Loads an annotation index from BED12, GTF, or GFF file.
    ///
    /// # Arguments
    /// * `path` - Path to annotation file
    ///
    /// # Example
    /// ```rust,ignore
    /// let index = AnnotationIndex::load(Path::new("annotation.gtf")).unwrap();
    /// ```
    pub fn load(path: &Path) -> Result<Self> {
        let mut builder = IndexBuilder::default();

        match detect_annotation_format(path)? {
            AnnotationFormat::Bed12 => {
                let mut reader = Reader::<Bed12>::from_mmap(path).with_context(|| {
                    format!("failed to open BED12 annotation {}", path.display())
                })?;
                for record in reader.records() {
                    let record = record.with_context(|| {
                        format!("failed to parse BED12 record from {}", path.display())
                    })?;
                    builder.push_record(record)?;
                }
            }
            AnnotationFormat::Gtf => {
                let mut reader = Reader::<Gtf>::from_mmap(path)
                    .with_context(|| format!("failed to open GTF annotation {}", path.display()))?;
                for record in reader.records() {
                    let record = record.with_context(|| {
                        format!("failed to parse GTF record from {}", path.display())
                    })?;
                    builder.push_record(record)?;
                }
            }
            AnnotationFormat::Gff => {
                let mut reader = Reader::<Gff>::from_mmap(path)
                    .with_context(|| format!("failed to open GFF annotation {}", path.display()))?;
                for record in reader.records() {
                    let record = record.with_context(|| {
                        format!("failed to parse GFF record from {}", path.display())
                    })?;
                    builder.push_record(record)?;
                }
            }
        }

        builder.finish()
    }

    /// Looks up a contig ID by name.
    pub fn contig_id(&self, name: &[u8]) -> Option<u32> {
        self.contig_ids.get(name).copied()
    }

    /// Gets a contig name by ID.
    pub fn contig_name(&self, contig_id: u32) -> &[u8] {
        &self.contigs[contig_id as usize]
    }

    /// Finds junctions near a genomic boundary.
    ///
    /// # Arguments
    /// * `contig_id` - Contig ID
    /// * `strand` - Strand orientation
    /// * `boundary` - Genomic boundary position
    /// * `wiggle` - Search radius
    pub fn junctions_near(
        &self,
        contig_id: u32,
        strand: Strand,
        boundary: u64,
        wiggle: u32,
    ) -> Vec<&Junction> {
        let start = boundary.saturating_sub(u64::from(wiggle));
        let end = boundary.saturating_add(u64::from(wiggle));
        let mut hits = Vec::new();

        for coord in start..=end {
            let key = BoundaryKey {
                contig_id,
                strand,
                boundary: coord,
            };
            if let Some(ids) = self.boundaries.get(&key) {
                hits.extend(ids.iter().map(|id| &self.junctions[*id]));
            }
        }

        hits
    }

    /// Returns the number of contigs in the index.
    pub fn contig_count(&self) -> usize {
        self.contigs.len()
    }

    /// Returns the number of junctions in the index.
    pub fn junction_count(&self) -> usize {
        self.junctions.len()
    }
}

/// Builder for constructing an AnnotationIndex.
#[derive(Default)]
struct IndexBuilder {
    contigs: Vec<Vec<u8>>,
    contig_ids: HashMap<Vec<u8>, u32>,
    junctions: Vec<Junction>,
    boundaries: HashMap<BoundaryKey, Vec<usize>>,
    seen: HashMap<JunctionGeometry, usize>,
}

impl IndexBuilder {
    /// Pushes a GenePred record to build the index.
    fn push_record(&mut self, record: GenePred) -> Result<()> {
        let strand = match record.strand() {
            Some(GenePredStrand::Forward) => Strand::Forward,
            Some(GenePredStrand::Reverse) => Strand::Reverse,
            _ => return Ok(()),
        };

        let exons = record.exons();
        if exons.len() < 2 {
            return Ok(());
        }

        let contig_id = self.intern_contig(record.chrom.as_slice());
        let mut transcript_exons: Vec<Exon> = exons
            .into_iter()
            .map(|(start, end)| Exon { start, end })
            .collect();

        if strand == Strand::Reverse {
            transcript_exons.reverse();
        }

        for pair in transcript_exons.windows(2) {
            let upstream = pair[0];
            let downstream = pair[1];
            let boundary = match strand {
                Strand::Forward => upstream.end,
                Strand::Reverse => upstream.start,
            };
            let intron_len = match strand {
                Strand::Forward => downstream.start.saturating_sub(upstream.end),
                Strand::Reverse => upstream.start.saturating_sub(downstream.end),
            };

            if intron_len == 0 {
                continue;
            }

            let geometry = JunctionGeometry {
                contig_id,
                strand,
                upstream,
                downstream,
                boundary,
            };

            let junction_id = if let Some(id) = self.seen.get(&geometry).copied() {
                id
            } else {
                let id = self.junctions.len();
                self.seen.insert(geometry, id);
                self.junctions.push(Junction {
                    id,
                    contig_id,
                    strand,
                    upstream,
                    downstream,
                    boundary,
                    intron_len,
                });
                id
            };

            let key = BoundaryKey {
                contig_id,
                strand,
                boundary,
            };
            self.boundaries.entry(key).or_default().push(junction_id);
        }

        Ok(())
    }

    /// Finishes building the index.
    fn finish(self) -> Result<AnnotationIndex> {
        if self.junctions.is_empty() {
            bail!("annotation did not yield any exon-to-next-exon junctions");
        }

        Ok(AnnotationIndex {
            contigs: self.contigs,
            contig_ids: self.contig_ids,
            junctions: self.junctions,
            boundaries: self.boundaries,
        })
    }

    /// Interns a contig name, returning its ID.
    fn intern_contig(&mut self, contig: &[u8]) -> u32 {
        if let Some(id) = self.contig_ids.get(contig) {
            *id
        } else {
            let id = self.contigs.len() as u32;
            let owned = contig.to_vec();
            self.contigs.push(owned.clone());
            self.contig_ids.insert(owned, id);
            id
        }
    }
}

/// Annotation file format.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum AnnotationFormat {
    Bed12,
    Gtf,
    Gff,
}

/// Detects annotation format from file extension.
fn detect_annotation_format(path: &Path) -> Result<AnnotationFormat> {
    let ext = normalized_extension(path)?;

    match ext.as_str() {
        "bed" => Ok(AnnotationFormat::Bed12),
        "gtf" => Ok(AnnotationFormat::Gtf),
        "gff" | "gff3" => Ok(AnnotationFormat::Gff),
        other => Err(anyhow!(
            "unsupported annotation extension .{other} for {}",
            path.display()
        )),
    }
}

/// Returns normalized file extension (handles compression).
fn normalized_extension(path: &Path) -> Result<String> {
    let ext = path
        .extension()
        .and_then(|ext| ext.to_str())
        .ok_or_else(|| anyhow!("cannot determine file extension for {}", path.display()))?;

    if matches!(
        ext.to_ascii_lowercase().as_str(),
        "gz" | "bgz" | "bz2" | "zst"
    ) {
        let inner = Path::new(path.file_stem().and_then(|stem| stem.to_str()).ok_or_else(
            || {
                anyhow!(
                    "cannot determine compressed inner extension for {}",
                    path.display()
                )
            },
        )?);
        return inner
            .extension()
            .and_then(|inner| inner.to_str())
            .map(|value| value.to_ascii_lowercase())
            .ok_or_else(|| {
                anyhow!(
                    "cannot determine compressed inner extension for {}",
                    path.display()
                )
            });
    }

    Ok(ext.to_ascii_lowercase())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalized_extension_handles_gzip_suffix() {
        let ext = normalized_extension(Path::new("annotations.gtf.gz")).unwrap();
        assert_eq!(ext, "gtf");
    }
}
