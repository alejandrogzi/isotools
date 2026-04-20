// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Error types for iso-adapter.

use std::path::PathBuf;

use thiserror::Error;

/// Errors produced by the iso-adapter pipeline.
#[derive(Debug, Error)]
pub enum AdapterError {
    #[error("I/O error on {path}: {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("failed to read BAM header from {path}")]
    BamHeader { path: PathBuf },

    #[error("failed to read BAM record at index {index}: {source}")]
    BamRecord {
        index: u64,
        #[source]
        source: std::io::Error,
    },

    #[error("failed to construct adapter matcher: {0}")]
    Matcher(String),

    #[error("reader/worker/writer channel closed unexpectedly")]
    ChannelClosed,

    #[error("output path {0} resolved to an existing directory without trailing separator")]
    OutputPath(PathBuf),
}

impl AdapterError {
    /// Wraps a std::io::Error with a path for context.
    pub fn io<P: Into<PathBuf>>(path: P, source: std::io::Error) -> Self {
        Self::Io {
            path: path.into(),
            source,
        }
    }
}
