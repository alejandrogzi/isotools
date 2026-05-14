// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Error types for iso-align.

use std::path::PathBuf;

use thiserror::Error;

#[derive(Debug, Error)]
pub enum AlignError {
    #[error("I/O error on {path}: {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("output format {0:?} requires a different code path; this is a bug")]
    OutputFormatMismatch(&'static str),

    #[error("FASTQ output requested but {path} carries no quality scores (FASTA input?)")]
    FastqWithoutQuality { path: PathBuf },
}
