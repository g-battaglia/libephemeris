// SPDX-License-Identifier: AGPL-3.0
//! Read-only SPK availability reports for a requested minor-body TDB grid.
//!
//! This command reports the largest matching file and its coverage envelope;
//! it does not fit elements or resolve internal descriptor gaps.

use super::element_sources::{
    DEFAULT_END_YEAR, DEFAULT_SPACING, DEFAULT_START_YEAR, approx_year_from_jd,
    collect_spk_candidates, file_stem_for_id, generate_epoch_grid, is_element_source_descriptor,
    resolve_asteroid_target, select_bodies,
};
use super::{Command, EXIT_FAILED, EXIT_OK, resolve_data_dir};
use crate::spk;
use std::path::{Path, PathBuf};

/// Largest candidate, with ties resolved by the lexicographically smallest
/// path. This reporting policy does not choose the element generator's source.
fn best_candidate(candidates: &[(PathBuf, u64)]) -> Option<&Path> {
    candidates
        .iter()
        .max_by(|a, b| a.1.cmp(&b.1).then_with(|| b.0.cmp(&a.0)))
        .map(|(path, _)| path.as_path())
}

/// Outer coverage and descriptor compatibility of one resolved target.
#[derive(Debug, Clone, PartialEq)]
struct BodyKernelRange {
    /// Minimum Julian Day (TDB) over the target's merged coverage.
    jd_min: f64,
    /// Maximum Julian Day (TDB) over the target's merged coverage.
    jd_max: f64,
    /// Whether any target segment has a native type-2/3/21 PV layout with
    /// center 10 (Sun) and frame 1 (J2000 equatorial).
    has_element_source: bool,
}

/// Resolve the configured legacy asteroid ID against either declared alias,
/// then union that target's decoded coverage into TDB Julian days. Invalid
/// images, absent targets or non-finite conversions return `None`. The union
/// includes all supported segment types; the separate flag identifies whether
/// any descriptor has the generator's supported Sun/J2000 layout.
fn kernel_body_range(bytes: &[u8], legacy_naif_id: i32) -> Option<BodyKernelRange> {
    let index = spk::SpkIndex::parse(bytes).ok()?;
    let naif_id = resolve_asteroid_target(&index, legacy_naif_id)?;
    let spans = index.coverage(naif_id);
    if spans.is_empty() {
        return None;
    }
    let mut jd_min = f64::INFINITY;
    let mut jd_max = f64::NEG_INFINITY;
    for (start_et, end_et) in spans {
        let low = spk::et_to_jd_tdb(start_et).ok()?;
        let high = spk::et_to_jd_tdb(end_et).ok()?;
        jd_min = jd_min.min(low);
        jd_max = jd_max.max(high);
    }
    let has_element_source = index
        .segments_for(naif_id)
        .any(|segment| is_element_source_descriptor(&segment.descriptor));
    Some(BodyKernelRange {
        jd_min,
        jd_max,
        has_element_source,
    })
}

/// Grid epochs strictly inside the kernel range with a
/// one-day guard on each side. This outer-envelope count does not resolve
/// internal descriptor gaps.
fn count_covered_epochs(grid: &[f64], jd_min: f64, jd_max: f64) -> usize {
    grid.iter()
        .filter(|jd| jd_min + 1.0 < **jd && **jd < jd_max - 1.0)
        .count()
}

/// One availability row rendered by [`render_dry_run_text`].
struct DryRunRow {
    name: String,
    const_name: String,
    /// Discovery path of the reported kernel, if any file matched.
    path: Option<String>,
    /// Decoded target range, if the kernel parsed and covers the body.
    range: Option<BodyKernelRange>,
    covered: usize,
}

/// Render the availability report over computed rows: header
/// with the grid parameters, the epoch grid in JD TDB, then one line
/// per body. This function performs no filesystem reads or state evaluations.
fn render_dry_run_text(
    start: i32,
    end: i32,
    spacing: i32,
    grid: &[f64],
    search_dirs: &[PathBuf],
    rows: &[DryRunRow],
) -> String {
    let mut out = String::from("Multi-epoch element generation (native dry-run)\n");
    out.push_str(&format!("  Epoch range: {start}–{end}\n"));
    out.push_str(&format!("  Spacing: {spacing} years\n"));
    out.push_str(&format!("  Epochs: {}\n", grid.len()));
    out.push_str(&format!("  Bodies: {}\n", rows.len()));
    out.push('\n');
    out.push_str("Epoch grid (JD TDB):\n");
    for jd in grid {
        out.push_str(&format!("  {:.0}: JD {jd:.1}\n", approx_year_from_jd(*jd)));
    }
    out.push('\n');
    let search = search_dirs
        .iter()
        .map(|dir| dir.display().to_string())
        .collect::<Vec<_>>()
        .join(", ");
    out.push_str(&format!("Body SPK availability (search: {search}):\n"));
    for row in rows {
        match (&row.path, &row.range) {
            (Some(_), Some(range)) => {
                let acceptance = if range.has_element_source {
                    "; type-2/3/21 heliocentric J2000"
                } else {
                    "; no type-2/3/21 heliocentric-J2000 segment"
                };
                out.push_str(&format!(
                    "  {:14} {:20} SPK: {:.0}–{:.0} ({}/{}{acceptance})\n",
                    row.name,
                    row.const_name,
                    approx_year_from_jd(range.jd_min),
                    approx_year_from_jd(range.jd_max),
                    row.covered,
                    grid.len(),
                ));
            }
            (Some(_), None) => {
                out.push_str(&format!(
                    "  {:14} {:20} SPK: (range unknown)\n",
                    row.name, row.const_name
                ));
            }
            (None, _) => {
                out.push_str(&format!("  {:14} {:20} NO SPK\n", row.name, row.const_name));
            }
        }
    }
    out.push_str("\nCoverage per body is the NAIF-target union over the natively decoded segment types (2/3/21); type-2/3/21 heliocentric-J2000 availability is reported separately. This command performs no fitting.\n");
    out.push_str("Generate elements with 'libephemeris generate keplerian-elements' and the same source/grid options.\n");
    out
}

/// Parse the words after `generate keplerian-dry-run`.
///
/// `--spk-dir` repeats (default: `<data-dir>/spk`, else `./spk`);
/// `--start`/`--end`/`--spacing` override the uniform TDB grid labels;
/// `--body` repeats a case-insensitive exact body-name filter.
pub(super) fn parse_dry_run_args<'a>(
    args: &mut std::iter::Peekable<impl Iterator<Item = &'a str>>,
) -> Result<Command, String> {
    let mut spk_dirs: Vec<String> = Vec::new();
    let mut data_dir: Option<String> = None;
    let mut start = DEFAULT_START_YEAR;
    let mut end = DEFAULT_END_YEAR;
    let mut spacing = DEFAULT_SPACING;
    let mut bodies: Vec<String> = Vec::new();
    while let Some(arg) = args.next() {
        if arg == "--spk-dir" {
            match args.next() {
                Some(dir) => spk_dirs.push(dir.to_string()),
                None => return Err("--spk-dir needs a value".to_string()),
            }
        } else if arg == "--data-dir" {
            match args.next() {
                Some(dir) => data_dir = Some(dir.to_string()),
                None => return Err("--data-dir needs a value".to_string()),
            }
        } else if arg == "--start" {
            match args.next() {
                Some(value) => match value.parse::<i32>() {
                    Ok(year) => start = year,
                    Err(_) => {
                        return Err(format!("--start needs an integer year, got '{value}'"));
                    }
                },
                None => return Err("--start needs a value".to_string()),
            }
        } else if arg == "--end" {
            match args.next() {
                Some(value) => match value.parse::<i32>() {
                    Ok(year) => end = year,
                    Err(_) => return Err(format!("--end needs an integer year, got '{value}'")),
                },
                None => return Err("--end needs a value".to_string()),
            }
        } else if arg == "--spacing" {
            match args.next() {
                Some(value) => match value.parse::<i32>() {
                    Ok(step) => spacing = step,
                    Err(_) => {
                        return Err(format!("--spacing needs an integer, got '{value}'"));
                    }
                },
                None => return Err("--spacing needs a value".to_string()),
            }
        } else if arg == "--body" {
            match args.next() {
                Some(name) => bodies.push(name.to_string()),
                None => return Err("--body needs a value".to_string()),
            }
        } else {
            return Err(format!("unexpected argument '{arg}'"));
        }
    }
    generate_epoch_grid(start, end, spacing)?;
    select_bodies(&bodies)?;
    Ok(Command::GenerateKeplerianDryRun {
        spk_dirs,
        data_dir,
        start,
        end,
        spacing,
        bodies,
    })
}

/// Run `generate keplerian-dry-run`: decode local kernels with the
/// native SPK reader and print the requested TDB epoch grid plus per-body
/// availability. Read-only over the search dirs; never creates,
/// downloads or modifies anything. Exit is 0 even when no kernel is
/// found (a dry-run reports, it does not provision).
pub(super) fn run_generate_keplerian_dry_run(
    spk_dirs: &[String],
    data_dir: Option<String>,
    start: i32,
    end: i32,
    spacing: i32,
    bodies: &[String],
) -> i32 {
    let grid = match generate_epoch_grid(start, end, spacing) {
        Ok(grid) => grid,
        Err(message) => {
            eprintln!("libephemeris: {message}");
            return EXIT_FAILED;
        }
    };
    let selected = match select_bodies(bodies) {
        Ok(selected) => selected,
        Err(message) => {
            eprintln!("libephemeris: {message}");
            return EXIT_FAILED;
        }
    };
    let mut search: Vec<PathBuf> = spk_dirs.iter().map(PathBuf::from).collect();
    if search.is_empty() {
        match resolve_data_dir(
            data_dir.as_deref(),
            crate::config::ConfigurationSources::discover(),
        ) {
            Ok(root) => search.push(root.join("spk")),
            Err(_) => search.push(PathBuf::from("spk")),
        }
    }
    let mut rows: Vec<DryRunRow> = Vec::with_capacity(selected.len());
    for body in &selected {
        let stems = vec![
            file_stem_for_id(body.horizons_id),
            body.horizons_id.to_string(),
            body.name.to_lowercase(),
        ];
        let candidates = collect_spk_candidates(&search, &stems);
        let found = best_candidate(&candidates).map(|path| path.to_path_buf());
        let (path, range, covered) = match &found {
            None => (None, None, 0),
            Some(path) => {
                let display = path.display().to_string();
                match std::fs::read(path) {
                    Ok(bytes) => match kernel_body_range(&bytes, body.legacy_naif_id) {
                        Some(body_range) => {
                            let covered =
                                count_covered_epochs(&grid, body_range.jd_min, body_range.jd_max);
                            (Some(display), Some(body_range), covered)
                        }
                        None => (Some(display), None, 0),
                    },
                    Err(_) => (Some(display), None, 0),
                }
            }
        };
        rows.push(DryRunRow {
            name: body.name.clone(),
            const_name: body.const_name.clone(),
            path,
            range,
            covered,
        });
    }
    let report = render_dry_run_text(start, end, spacing, &grid, &search, &rows);
    print!("{report}");
    EXIT_OK
}

#[cfg(test)]
mod tests {
    use super::super::element_sources::hand_kernel;
    use super::*;

    fn scratch_dir(tag: &str) -> PathBuf {
        let dir =
            std::env::temp_dir().join(format!("leph-generate-test-{}-{tag}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn spk_candidate_discovery_prefers_the_largest_file() {
        let dir = scratch_dir("candidates");
        std::fs::write(dir.join("ceres_small.bsp"), vec![0u8; 10]).unwrap();
        std::fs::write(dir.join("ceres_big.bsp"), vec![0u8; 20]).unwrap();
        std::fs::write(dir.join("ceres_notes.txt"), vec![0u8; 30]).unwrap();
        std::fs::write(dir.join("vesta_big.bsp"), vec![0u8; 40]).unwrap();
        let found = collect_spk_candidates(std::slice::from_ref(&dir), &["ceres".to_string()]);
        assert_eq!(found.len(), 2);
        assert_eq!(
            best_candidate(&found),
            Some(dir.join("ceres_big.bsp").as_path())
        );
        assert!(collect_spk_candidates(&[dir.join("absent")], &["ceres".to_string()]).is_empty());
        assert!(best_candidate(&[]).is_none());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn covered_epoch_count_applies_the_one_day_guard() {
        let grid = vec![100.0, 101.0, 102.0, 199.0, 200.0];
        assert_eq!(count_covered_epochs(&grid, 99.0, 201.0), 3);
        assert_eq!(count_covered_epochs(&grid, 100.0, 200.0), 1);
        assert_eq!(count_covered_epochs(&grid, 150.0, 160.0), 0);
    }

    #[test]
    fn kernel_body_range_checks_native_layouts_and_rejects_other_centers_and_frames() {
        let img = hand_kernel(&[(-1_000_000.0, 1_000_000.0, 2000001, 10, 1, 21)]);
        let range = kernel_body_range(&img, 2000001).unwrap();
        assert!(range.has_element_source);
        let expected_min = spk::et_to_jd_tdb(-1_000_000.0).unwrap();
        let expected_max = spk::et_to_jd_tdb(1_000_000.0).unwrap();
        assert_eq!(range.jd_min, expected_min);
        assert_eq!(range.jd_max, expected_max);
        for (kind, center, frame, accepted) in [
            (2, 10, 1, true),
            (3, 10, 1, true),
            (21, 0, 1, false),
            (2, 0, 1, false),
            (3, 10, 2, false),
        ] {
            let image = hand_kernel(&[(-1_000_000.0, 1_000_000.0, 2000001, center, frame, kind)]);
            let range = kernel_body_range(&image, 2000001).unwrap();
            assert_eq!(range.has_element_source, accepted);
            assert_eq!(range.jd_min, expected_min);
            assert_eq!(range.jd_max, expected_max);
        }
        assert!(kernel_body_range(&img, 999).is_none());
        assert!(kernel_body_range(b"not a kernel", 2000001).is_none());
        assert!(kernel_body_range(&[], 2000001).is_none());
        let truncated = hand_kernel(&[(-1_000_000.0, 1_000_000.0, 2000001, 10, 1, 21)]);
        assert!(kernel_body_range(&truncated[..1024], 2000001).is_none());
    }

    #[test]
    fn dry_run_resolves_declared_aliases_without_merging_other_targets() {
        let modern = (-1_000_000.0, 1_000_000.0, 20_000_001, 10, 1, 21);
        let legacy = (-2_000_000.0, 2_000_000.0, 2_000_001, 10, 1, 21);
        let expected_min = spk::et_to_jd_tdb(modern.0).unwrap();
        let expected_max = spk::et_to_jd_tdb(modern.1).unwrap();
        for descriptors in [vec![modern], vec![legacy, modern]] {
            let image = hand_kernel(&descriptors);
            let range = kernel_body_range(&image, 2_000_001).unwrap();
            assert!(range.has_element_source);
            assert_eq!(range.jd_min, expected_min);
            assert_eq!(range.jd_max, expected_max);
        }
        let different_object = hand_kernel(&[(-1_000_000.0, 1_000_000.0, 20_000_002, 10, 1, 21)]);
        assert!(kernel_body_range(&different_object, 2_000_001).is_none());
    }

    #[test]
    fn dry_run_report_marks_found_and_missing_kernels() {
        let grid = generate_epoch_grid(2000, 2010, 10).unwrap();
        let rows = vec![
            DryRunRow {
                name: "Ceres".to_string(),
                const_name: "CERES".to_string(),
                path: Some("/kernels/ceres_x.bsp".to_string()),
                range: Some(BodyKernelRange {
                    jd_min: grid[0] - 100.0,
                    jd_max: grid[1] + 100.0,
                    has_element_source: true,
                }),
                covered: 2,
            },
            DryRunRow {
                name: "Vesta".to_string(),
                const_name: "VESTA".to_string(),
                path: None,
                range: None,
                covered: 0,
            },
        ];
        let text = render_dry_run_text(2000, 2010, 10, &grid, &[PathBuf::from("/kernels")], &rows);
        assert!(text.contains("Epochs: 2"));
        assert!(text.contains("Bodies: 2"));
        assert!(text.contains("type-2/3/21 heliocentric J2000"));
        assert!(text.contains("NO SPK"));
    }
}
