// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! iso-adapter: detect and optionally remove adapter sequences from
//! soft-clipped regions of long-read BAM alignments.

fn main() {
    if let Err(err) = iso_adapter::run() {
        eprintln!("error: {err:#}");
        std::process::exit(1);
    }
}
