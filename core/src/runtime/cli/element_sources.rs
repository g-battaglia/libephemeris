// SPDX-License-Identifier: AGPL-3.0
//! Shared minor-body catalog, uniform TDB epochs and local SPK source keys.
//!
//! Discovery only enumerates candidates. Each command retains its own policy
//! for choosing a kernel or refusing ambiguity.

use crate::constants::SPK_BODY_NAME_MAP;
use crate::{minor_bodies, spk};
use std::path::PathBuf;

/// Defaults of the multi-epoch generation grid:
/// Julian-year labels 1650-2450 at ten-year spacing, 81 epochs.
pub(super) const DEFAULT_START_YEAR: i32 = 1650;
/// Grid end label (see [`DEFAULT_START_YEAR`]).
pub(super) const DEFAULT_END_YEAR: i32 = 2450;
/// Grid spacing in Julian years (see [`DEFAULT_START_YEAR`]).
pub(super) const DEFAULT_SPACING: i32 = 10;
/// J2000.0 reference of the uniform generation grid.
const GRID_J2000_JD: f64 = 2451545.0;
/// Julian-year length of the uniform generation grid.
const GRID_DAYS_PER_JULIAN_YEAR: f64 = 365.25;

/// TDB Julian day for a uniform Julian-year label:
/// `2451545.0 + (year - 2000) * 365.25`. Labels are not civil-calendar instants.
fn grid_year_to_jd(year: i32) -> f64 {
    GRID_J2000_JD + (f64::from(year) - 2000.0) * GRID_DAYS_PER_JULIAN_YEAR
}

/// Uniform Julian-year label for a Julian day, used only to render reports.
pub(super) fn approx_year_from_jd(jd: f64) -> f64 {
    2000.0 + (jd - GRID_J2000_JD) / GRID_DAYS_PER_JULIAN_YEAR
}

/// Uniform TDB epochs between inclusive Julian-year labels
/// spaced by `spacing` label years; the end is emitted only when it lies on
/// the grid. Nonpositive spacing or reversed bounds return a usage error.
pub(super) fn generate_epoch_grid(start: i32, end: i32, spacing: i32) -> Result<Vec<f64>, String> {
    if spacing <= 0 {
        return Err("spacing must be a positive number of years".to_string());
    }
    if end < start {
        return Err(format!("end year {end} precedes start year {start}"));
    }
    let mut epochs = Vec::new();
    let mut year = start;
    loop {
        epochs.push(grid_year_to_jd(year));
        match year.checked_add(spacing) {
            Some(next) if next <= end => year = next,
            _ => break,
        }
    }
    Ok(epochs)
}

/// Discovery stem formed by lowercasing a Horizons ID, replacing characters
/// other than alphanumerics, `-` and `_`, and stripping trailing underscores.
pub(super) fn file_stem_for_id(horizons_id: &str) -> String {
    horizons_id
        .to_lowercase()
        .chars()
        .map(|c| {
            if c.is_alphanumeric() || c == '-' || c == '_' {
                c
            } else {
                '_'
            }
        })
        .collect::<String>()
        .trim_end_matches('_')
        .to_string()
}

/// Exported element constant for a display name, avoiding planet-name collisions.
fn element_constant_name(body_name: &str) -> String {
    match body_name {
        "Europa" => "EUROPA_AST".to_string(),
        "Pandora" => "PANDORA_AST".to_string(),
        "Lilith" => "LILITH_AST".to_string(),
        _ => body_name.to_uppercase(),
    }
}

/// One configured minor body and its local-kernel discovery keys.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct ElementBody {
    pub(super) body_id: i32,
    pub(super) name: String,
    pub(super) const_name: String,
    pub(super) horizons_id: &'static str,
    pub(super) legacy_naif_id: i32,
}

/// Select configured minor bodies by case-insensitive exact display name.
/// The empty filter selects all bodies; results are sorted by name. Unknown
/// names return a usage error that lists the available choices.
pub(super) fn select_bodies(filter: &[String]) -> Result<Vec<ElementBody>, String> {
    let mut bodies: Vec<ElementBody> = SPK_BODY_NAME_MAP
        .iter()
        .map(|(body_id, horizons_id, legacy_naif_id)| {
            let name = minor_bodies::minor_body_elements(*body_id)
                .map(|elements| elements.name().to_string())
                .unwrap_or_else(|| horizons_id.to_string());
            ElementBody {
                body_id: *body_id,
                const_name: element_constant_name(&name),
                name,
                horizons_id,
                legacy_naif_id: *legacy_naif_id,
            }
        })
        .collect();
    if !filter.is_empty() {
        let wanted: Vec<String> = filter.iter().map(|name| name.to_lowercase()).collect();
        for name in &wanted {
            if !bodies.iter().any(|body| body.name.to_lowercase() == *name) {
                let mut available: Vec<&str> =
                    bodies.iter().map(|body| body.name.as_str()).collect();
                available.sort_unstable();
                return Err(format!(
                    "unknown body '{name}' (available: {})",
                    available.join(", ")
                ));
            }
        }
        bodies.retain(|body| wanted.contains(&body.name.to_lowercase()));
    }
    bodies.sort_by(|a, b| a.name.cmp(&b.name));
    Ok(bodies)
}

/// Candidate kernels on disk: `(path, byte size)` pairs whose filename
/// matches one of the discovery stems (`{safe_id}_*.bsp`,
/// `{horizons_id}_*.bsp`, `{name}_*.bsp`) with a `.bsp` suffix.
/// Missing/unreadable directories and unreadable entries are skipped; failed
/// metadata reads contribute a size of zero. Paths are sorted and deduplicated,
/// keeping their first size. This helper does not select a preferred candidate.
pub(super) fn collect_spk_candidates(dirs: &[PathBuf], stems: &[String]) -> Vec<(PathBuf, u64)> {
    let prefixes: Vec<String> = stems.iter().map(|stem| format!("{stem}_")).collect();
    let mut out: Vec<(PathBuf, u64)> = Vec::new();
    for dir in dirs {
        let entries = match std::fs::read_dir(dir) {
            Ok(entries) => entries,
            Err(_) => continue,
        };
        for entry in entries.flatten() {
            let name = entry.file_name().to_string_lossy().into_owned();
            if !name.ends_with(".bsp") {
                continue;
            }
            if !prefixes.iter().any(|prefix| name.starts_with(prefix)) {
                continue;
            }
            let size = entry.metadata().map(|meta| meta.len()).unwrap_or(0);
            out.push((entry.path(), size));
        }
    }
    out.sort_by(|a, b| a.0.cmp(&b.0));
    out.dedup_by(|a, b| a.0 == b.0);
    out
}

/// Resolve only the two declared codes of one numbered asteroid.
/// The shared SPK policy prefers the Horizons code when both aliases exist.
pub(super) fn resolve_asteroid_target(index: &spk::SpkIndex, legacy_naif_id: i32) -> Option<i32> {
    let number = i64::from(legacy_naif_id) - spk::NAIF_ASTEROID_OFFSET_LEGACY;
    let targets: Vec<_> = index.all_segments().map(|s| s.descriptor.target).collect();
    spk::find_naif_id_for_asteroid(&targets, number)
}

/// Native PV layouts suitable for osculation: type 2/3/21, Sun-relative,
/// J2000 equatorial. Target identity, selected-segment priority and TDB
/// coverage must still be checked separately by the caller.
pub(super) fn is_element_source_descriptor(descriptor: &spk::Descriptor) -> bool {
    matches!(descriptor.seg_type, 2 | 3 | 21) && descriptor.center == 10 && descriptor.frame == 1
}

/// Minimal little-endian DAF image: file + summary + name + element
/// records with the given `(start_et, end_et, target, center, frame,
/// type)` segments. Element words are never read by the dry-run
/// path, so the addresses only need to land inside the image.
#[cfg(test)]
pub(super) fn hand_kernel(segments: &[(f64, f64, i32, i32, i32, i32)]) -> Vec<u8> {
    let mut img = vec![0u8; 4 * 1024];
    img[0..8].copy_from_slice(b"DAF/SPK ");
    img[8..12].copy_from_slice(&2i32.to_le_bytes());
    img[12..16].copy_from_slice(&6i32.to_le_bytes());
    img[76..80].copy_from_slice(&2i32.to_le_bytes());
    img[80..84].copy_from_slice(&2i32.to_le_bytes());
    img[88..96].copy_from_slice(b"LTL-IEEE");
    let base = 1024usize;
    img[base..base + 8].copy_from_slice(&0f64.to_le_bytes());
    img[base + 8..base + 16].copy_from_slice(&0f64.to_le_bytes());
    img[base + 16..base + 24].copy_from_slice(&(segments.len() as f64).to_le_bytes());
    for (i, (start, end, target, center, frame, seg_type)) in segments.iter().enumerate() {
        let off = base + (3 + i * 5) * 8;
        img[off..off + 8].copy_from_slice(&start.to_le_bytes());
        img[off + 8..off + 16].copy_from_slice(&end.to_le_bytes());
        for (k, v) in [*target, *center, *frame, *seg_type, 385i32, 394i32]
            .iter()
            .enumerate()
        {
            let w = off + 16 + (k / 2) * 8 + (k % 2) * 4;
            img[w..w + 4].copy_from_slice(&v.to_le_bytes());
        }
    }
    img
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn generation_grid_matches_the_reviewed_labels() {
        assert_eq!(grid_year_to_jd(2000), 2451545.0);
        assert_eq!(grid_year_to_jd(1650), 2323707.5);
        assert_eq!(grid_year_to_jd(2450), 2615907.5);
        assert_eq!(approx_year_from_jd(2451545.0), 2000.0);
        let grid = generate_epoch_grid(1650, 2450, 10).unwrap();
        assert_eq!(grid.len(), 81);
        assert_eq!(grid[0], 2323707.5);
        assert_eq!(grid[80], 2615907.5);
        assert!(generate_epoch_grid(2000, 1990, 10).is_err());
        assert!(generate_epoch_grid(2000, 2010, 0).is_err());
    }

    #[test]
    fn file_stems_follow_the_release_discovery_shape() {
        assert_eq!(file_stem_for_id("Ceres;"), "ceres");
        assert_eq!(file_stem_for_id("2060"), "2060");
        assert_eq!(file_stem_for_id("Bennu;"), "bennu");
        assert_eq!(file_stem_for_id("10199"), "10199");
    }

    #[test]
    fn element_constant_names_avoid_planet_name_collisions() {
        assert_eq!(element_constant_name("Europa"), "EUROPA_AST");
        assert_eq!(element_constant_name("Pandora"), "PANDORA_AST");
        assert_eq!(element_constant_name("Lilith"), "LILITH_AST");
        assert_eq!(element_constant_name("Ceres"), "CERES");
    }

    #[test]
    fn dry_run_bodies_cover_37_sorted_names() {
        let bodies = select_bodies(&[]).unwrap();
        assert_eq!(bodies.len(), 37);
        let names: Vec<&str> = bodies.iter().map(|body| body.name.as_str()).collect();
        let mut sorted = names.clone();
        sorted.sort_unstable();
        assert_eq!(names, sorted);
        assert!(names.contains(&"Ceres"));
        let filtered = select_bodies(&["cErEs".to_string()]).unwrap();
        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered[0].name, "Ceres");
        let unknown = select_bodies(&["Nope".to_string()]).unwrap_err();
        assert!(unknown.contains("available:"));
        assert!(unknown.contains("Ceres"));
    }
}
