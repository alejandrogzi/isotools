// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Core module for detecting real 3' ends on isoseq reads
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main function for detecting real 3' ends
//! on isoseq reads and processing the components of reads in parallel.
//!
//! In short, uses APARENT and HMM polyA scoring systems to define real ends.  
//! Finally, the caller module groups all the previous information and tries to
//! determine the intraprimming potential for each read,

pub mod cli;
pub mod core;
