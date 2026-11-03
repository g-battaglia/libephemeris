// SPDX-License-Identifier: AGPL-3.0
//! Native executable adapter for the shared command-line application.
fn main() {
    std::process::exit(libephemeris_core::cli::run(std::env::args().skip(1)));
}
