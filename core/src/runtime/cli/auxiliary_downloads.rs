// SPDX-License-Identifier: AGPL-3.0
//! Explicit acquisition commands for ASSIST force data and IERS observations.
//!
//! Option parsing is pure. Execution resolves local destinations, applies the
//! explicit download policy and reports native transport results. The parent
//! dispatches commands and retains the shared data-directory resolver.

use super::{Command, EXIT_FAILED, EXIT_OK, EXIT_USAGE, resolve_data_dir};
use crate::context::EffectiveNetworkPolicy;
use crate::download;
use crate::net::TransportConfig;
use std::path::{Path, PathBuf};

/// Parse ASSIST file selections and destination/output options without I/O.
///
/// Both standard datasets are selected unless `--no-planets` or
/// `--no-asteroids` removes one. `--target-dir` takes precedence over the
/// data root. `--no-progress` suppresses per-file lines; `--quiet` suppresses
/// all normal output. The optional deep-time planetary file is not selected.
pub(super) fn parse_download_assist_args<'a>(
    args: &mut std::iter::Peekable<impl Iterator<Item = &'a str>>,
) -> Result<Command, String> {
    let mut target_dir: Option<String> = None;
    let mut data_dir: Option<String> = None;
    let mut planets = true;
    let mut asteroids = true;
    let mut force = false;
    let mut quiet = false;
    let mut no_progress = false;
    while let Some(arg) = args.next() {
        if arg == "--force" || arg == "-f" {
            force = true;
        } else if arg == "--quiet" || arg == "-q" {
            quiet = true;
        } else if arg == "--no-progress" {
            no_progress = true;
        } else if arg == "--no-planets" {
            planets = false;
        } else if arg == "--no-asteroids" {
            asteroids = false;
        } else if arg == "--target-dir" {
            match args.next() {
                Some(dir) => target_dir = Some(dir.to_string()),
                None => return Err("--target-dir needs a value".to_string()),
            }
        } else if arg == "--data-dir" {
            match args.next() {
                Some(dir) => data_dir = Some(dir.to_string()),
                None => return Err("--data-dir needs a value".to_string()),
            }
        } else {
            return Err(format!("unexpected argument '{arg}'"));
        }
    }
    Ok(Command::DownloadAssist {
        target_dir,
        data_dir,
        planets,
        asteroids,
        force,
        quiet,
        no_progress,
    })
}

/// Parse IERS refresh options without I/O.
///
/// `--force` refreshes even a fresh cache, `--quiet` suppresses normal output,
/// and `--data-dir` overrides the configured data root.
pub(super) fn parse_download_iers_args<'a>(
    args: &mut std::iter::Peekable<impl Iterator<Item = &'a str>>,
) -> Result<Command, String> {
    let mut data_dir: Option<String> = None;
    let mut force = false;
    let mut quiet = false;
    while let Some(arg) = args.next() {
        if arg == "--force" || arg == "-f" {
            force = true;
        } else if arg == "--quiet" || arg == "-q" {
            quiet = true;
        } else if arg == "--data-dir" {
            match args.next() {
                Some(dir) => data_dir = Some(dir.to_string()),
                None => return Err("--data-dir needs a value".to_string()),
            }
        } else {
            return Err(format!("unexpected argument '{arg}'"));
        }
    }
    Ok(Command::DownloadIers {
        data_dir,
        force,
        quiet,
    })
}

/// Resolve an explicit ASSIST destination or the data root's `assist` directory.
fn resolve_assist_dir(target_dir: Option<&str>, data_dir: &Path) -> PathBuf {
    match target_dir {
        Some(dir) => PathBuf::from(dir),
        None => data_dir.join("assist"),
    }
}

pub(super) fn run_download_assist(
    target_dir: Option<String>,
    data_dir: Option<String>,
    planets: bool,
    asteroids: bool,
    force: bool,
    quiet: bool,
    no_progress: bool,
) -> i32 {
    let data_dir = match resolve_data_dir(
        data_dir.as_deref(),
        crate::config::ConfigurationSources::discover(),
    ) {
        Ok(dir) => dir,
        Err(message) => {
            eprintln!("libephemeris: {message}");
            return EXIT_USAGE;
        }
    };
    let dest = resolve_assist_dir(target_dir.as_deref(), &data_dir);
    let config = match TransportConfig::new() {
        Ok(config) => config,
        Err(error) => {
            eprintln!("libephemeris: cannot configure transport: {error}");
            return EXIT_FAILED;
        }
    };
    let show_lines = !quiet && !no_progress;
    if !quiet {
        println!("libephemeris download assist → {}", dest.display());
    }
    // The DE441 deep-time file is not requested by the release CLI
    // spelling either, so it stays deselected here by construction.
    match download::download_assist_data(
        &dest,
        planets,
        asteroids,
        false,
        force,
        EffectiveNetworkPolicy::Allow,
        &config,
        None,
    ) {
        Ok(report) => {
            if show_lines {
                if report.downloaded == 0 && report.skipped == 0 {
                    println!("  [--] nothing requested (--no-planets --no-asteroids)");
                } else {
                    println!(
                        "  [OK] assist: {} downloaded, {} kept",
                        report.downloaded, report.skipped
                    );
                }
            }
            if !quiet {
                println!("download assist complete");
            }
            EXIT_OK
        }
        Err(error) => {
            eprintln!("libephemeris: download failed: {error}");
            EXIT_FAILED
        }
    }
}

/// Refresh the three standard IERS files under the data root's `iers_cache`.
///
/// Per-file lines identify both reused and refreshed files. `--quiet`
/// suppresses normal output while errors remain visible.
pub(super) fn run_download_iers(data_dir: Option<String>, force: bool, quiet: bool) -> i32 {
    let data_dir = match resolve_data_dir(
        data_dir.as_deref(),
        crate::config::ConfigurationSources::discover(),
    ) {
        Ok(dir) => dir,
        Err(message) => {
            eprintln!("libephemeris: {message}");
            return EXIT_USAGE;
        }
    };
    let cache_dir = download::iers_cache_dir(&data_dir);
    let config = match TransportConfig::new() {
        Ok(config) => config,
        Err(error) => {
            eprintln!("libephemeris: cannot configure transport: {error}");
            return EXIT_FAILED;
        }
    };
    if !quiet {
        println!("libephemeris download iers → {}", cache_dir.display());
    }
    match download::refresh_iers_cache(&cache_dir, force, EffectiveNetworkPolicy::Allow, &config) {
        Ok(report) => {
            if !quiet {
                for file in &report.files {
                    println!("  [OK] {}: {}", file.name, file.path.display());
                }
                println!("download iers complete");
            }
            EXIT_OK
        }
        Err(error) => {
            eprintln!("libephemeris: download failed: {error}");
            EXIT_FAILED
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn assist_dir_prefers_target_dir_over_the_data_root() {
        assert_eq!(
            resolve_assist_dir(Some("/explicit"), &PathBuf::from("/data")),
            PathBuf::from("/explicit")
        );
        assert_eq!(
            resolve_assist_dir(None, &PathBuf::from("/data")),
            PathBuf::from("/data/assist")
        );
    }
}
