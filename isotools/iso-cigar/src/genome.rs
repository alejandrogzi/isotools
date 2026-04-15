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
//! exon (or any other downstream), but the aligner missed the intron.
//!
//! In short, iso-cigar identifies candidate reads: those ending ±wiggle bp from a
//! known internal exon boundary, with soft-clip ≥ cutoff; retrieves the sequence of
//! the immediately downstream exon from the reference genome; checks whether the
//! soft-clipped bases match the beginning of that downstream exon; if yes: rewrites
//! the CIGAR to add an N intron gap and converts the soft-clip into a = match

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{anyhow, bail, Context, Result};
use flate2::read::MultiGzDecoder;
use twobit::TwoBitFile;

#[derive(Clone)]
/// A genome source that can be opened from 2bit or FASTA files.
///
/// # Example
/// ```rust,ignore
/// let source = GenomeSource::open(Path::new("genome.2bit")).unwrap();
/// ```
pub struct GenomeSource {
    inner: Arc<GenomeInner>,
}

/// Internal genome source variants.
enum GenomeInner {
    TwoBit {
        path: PathBuf,
    },
    Loaded {
        path: PathBuf,
        sequences: Arc<HashMap<Vec<u8>, Vec<u8>>>,
    },
}

/// A session for reading sequences from a genome source.
pub struct GenomeSession {
    inner: SessionInner,
}

/// Internal genome session variants.
enum SessionInner {
    TwoBit {
        path: PathBuf,
        reader: twobit::TwoBitPhysicalFile,
    },
    Loaded {
        path: PathBuf,
        sequences: Arc<HashMap<Vec<u8>, Vec<u8>>>,
    },
}

impl GenomeSource {
    /// Opens a genome source from a 2bit or FASTA file.
    ///
    /// # Arguments
    /// * `path` - Path to the genome file (.2bit, .fa, .fasta, .fna, or gzipped FASTA)
    ///
    /// # Example
    /// ```rust,ignore
    /// let source = GenomeSource::open(Path::new("genome.fa")).unwrap();
    /// ```
    pub fn open(path: &Path) -> Result<Self> {
        match detect_sequence_format(path)? {
            SequenceFormat::TwoBit => Ok(Self {
                inner: Arc::new(GenomeInner::TwoBit {
                    path: path.to_path_buf(),
                }),
            }),
            SequenceFormat::Fasta => Ok(Self {
                inner: Arc::new(GenomeInner::Loaded {
                    path: path.to_path_buf(),
                    sequences: Arc::new(load_fasta(path)?),
                }),
            }),
        }
    }

    /// Creates a new genome session for reading sequences.
    pub fn session(&self) -> Result<GenomeSession> {
        match self.inner.as_ref() {
            GenomeInner::TwoBit { path } => {
                let reader = TwoBitFile::open(path)
                    .with_context(|| format!("failed to open 2bit genome {}", path.display()))?
                    .enable_softmask(true);

                Ok(GenomeSession {
                    inner: SessionInner::TwoBit {
                        path: path.clone(),
                        reader,
                    },
                })
            }
            GenomeInner::Loaded { path, sequences } => Ok(GenomeSession {
                inner: SessionInner::Loaded {
                    path: path.clone(),
                    sequences: Arc::clone(sequences),
                },
            }),
        }
    }

    /// Returns a static string describing the genome format ("2bit" or "fasta").
    pub fn describe(&self) -> &'static str {
        match self.inner.as_ref() {
            GenomeInner::TwoBit { .. } => "2bit",
            GenomeInner::Loaded { .. } => "fasta",
        }
    }
}

impl GenomeSession {
    /// Fetches a sequence interval from a contig.
    ///
    /// # Arguments
    /// * `contig` - Contig name as bytes
    /// * `start` - Start position (0-based)
    /// * `end` - End position (exclusive)
    ///
    /// # Example
    /// ```rust,ignore
    /// let mut session = source.session().unwrap();
    /// let seq = session.fetch(b"chr1", 0, 100).unwrap();
    /// ```
    pub fn fetch(&mut self, contig: &[u8], start: u64, end: u64) -> Result<Vec<u8>> {
        if start > end {
            bail!("invalid sequence interval: {start} > {end}");
        }

        match &mut self.inner {
            SessionInner::TwoBit { path, reader } => {
                let contig = std::str::from_utf8(contig).with_context(|| {
                    format!("2bit contig names must be valid UTF-8: {:?}", contig)
                })?;
                let seq = reader
                    .read_sequence(contig, start as usize..end as usize)
                    .with_context(|| {
                        format!(
                            "failed to fetch {}:{}-{} from {}",
                            contig,
                            start,
                            end,
                            path.display()
                        )
                    })?;
                let mut seq = seq.into_bytes();
                seq.make_ascii_uppercase();
                Ok(seq)
            }
            SessionInner::Loaded { path, sequences } => {
                let sequence = sequences.get(contig).ok_or_else(|| {
                    anyhow!(
                        "contig {} not present in genome {}",
                        String::from_utf8_lossy(contig),
                        path.display()
                    )
                })?;

                let start = start as usize;
                let end = end as usize;
                if end > sequence.len() {
                    bail!(
                        "requested {}:{}-{} beyond genome length {} in {}",
                        String::from_utf8_lossy(contig),
                        start,
                        end,
                        sequence.len(),
                        path.display()
                    );
                }

                Ok(sequence[start..end].to_vec())
            }
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum SequenceFormat {
    TwoBit,
    Fasta,
}

/// Detects genome file format from extension.
fn detect_sequence_format(path: &Path) -> Result<SequenceFormat> {
    let ext = path
        .extension()
        .and_then(|ext| ext.to_str())
        .ok_or_else(|| anyhow!("cannot determine genome extension for {}", path.display()))?;

    if ext.eq_ignore_ascii_case("2bit") {
        return Ok(SequenceFormat::TwoBit);
    }

    if is_fasta_extension(ext) {
        return Ok(SequenceFormat::Fasta);
    }

    if ext.eq_ignore_ascii_case("gz") {
        let inner = Path::new(path.file_stem().and_then(|stem| stem.to_str()).ok_or_else(
            || {
                anyhow!(
                    "cannot determine compressed genome extension for {}",
                    path.display()
                )
            },
        )?);
        let inner_ext = inner
            .extension()
            .and_then(|value| value.to_str())
            .ok_or_else(|| {
                anyhow!(
                    "cannot determine compressed genome extension for {}",
                    path.display()
                )
            })?;

        if is_fasta_extension(inner_ext) {
            return Ok(SequenceFormat::Fasta);
        }
    }

    bail!(
        "unsupported genome format for {} (expected .2bit, .fa, .fasta, .fna, or gzipped FASTA)",
        path.display()
    )
}

/// Loads a FASTA genome file into a hash map.
fn load_fasta(path: &Path) -> Result<HashMap<Vec<u8>, Vec<u8>>> {
    let file = File::open(path)
        .with_context(|| format!("failed to open FASTA genome {}", path.display()))?;
    let gzip = path
        .extension()
        .and_then(|ext| ext.to_str())
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"));

    let mut reader: Box<dyn BufRead> = if gzip {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut sequences = HashMap::new();
    let mut line = Vec::new();
    let mut current_name: Option<Vec<u8>> = None;
    let mut current_sequence = Vec::new();

    loop {
        line.clear();
        let bytes = reader.read_until(b'\n', &mut line)?;
        if bytes == 0 {
            break;
        }

        trim_newline(&mut line);
        if line.is_empty() {
            continue;
        }

        if line[0] == b'>' {
            if let Some(name) = current_name.replace(parse_record_name(&line[1..]).to_vec()) {
                current_sequence.make_ascii_uppercase();
                sequences.insert(name, std::mem::take(&mut current_sequence));
            }
        } else {
            if current_name.is_none() {
                bail!(
                    "invalid FASTA {}: sequence before first header",
                    path.display()
                );
            }
            current_sequence.extend_from_slice(&line);
        }
    }

    if let Some(name) = current_name {
        current_sequence.make_ascii_uppercase();
        sequences.insert(name, current_sequence);
    }

    if sequences.is_empty() {
        bail!("no FASTA records found in {}", path.display());
    }

    Ok(sequences)
}

/// Trims newline characters from a line.
fn trim_newline(line: &mut Vec<u8>) {
    if line.ends_with(b"\n") {
        line.pop();
    }
    if line.ends_with(b"\r") {
        line.pop();
    }
}

/// Parses a FASTA record name (first non-whitespace sequence).
fn parse_record_name(line: &[u8]) -> &[u8] {
    let start = line
        .iter()
        .position(|byte| !byte.is_ascii_whitespace())
        .unwrap_or(line.len());
    let rest = &line[start..];
    let end = rest
        .iter()
        .position(|byte| byte.is_ascii_whitespace())
        .unwrap_or(rest.len());
    &rest[..end]
}

/// Returns true if extension is a FASTA variant.
fn is_fasta_extension(ext: &str) -> bool {
    ext.eq_ignore_ascii_case("fa")
        || ext.eq_ignore_ascii_case("fasta")
        || ext.eq_ignore_ascii_case("fna")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_record_name_stops_at_whitespace() {
        assert_eq!(parse_record_name(b"chr1 some desc"), b"chr1");
    }
}
