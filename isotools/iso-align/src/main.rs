// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! iso-align: select long reads whose pass-1 split alignment suggests a
//! re-align with a larger minimap2 intron cap.

fn main() {
    if let Err(err) = iso_align::run() {
        eprintln!("error: {err:#}");
        std::process::exit(1);
    }
}
