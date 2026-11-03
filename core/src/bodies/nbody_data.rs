// SPDX-License-Identifier: AGPL-3.0
//! ASSIST dataset locations and structural installation checks.

use std::path::Path;

/// JPL DE440 planetary integration dataset.
pub const ASSIST_PLANETS_URL: &str =
    "https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440/linux_p1550p2650.440";
/// JPL DE441 planetary dataset for deep-time integrations.
pub const ASSIST_PLANETS_DE441_URL: &str =
    "https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de441/linux_m13000p17000.441";
/// JPL sixteen-body asteroid perturber ephemeris.
pub const ASSIST_ASTEROIDS_URL: &str =
    "https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/sb441-n16.bsp";
/// Expected bytes for the DE440, DE441 and asteroid files, respectively.
pub const ASSIST_EXPECTED_SIZES: [u64; 3] = [102_272_352, 2_788_676_624, 645_727_232];
/// Installation size floor; coefficient evaluation additionally validates records.
pub const ASSIST_MIN_SIZE_RATIO: f64 = 0.90;

/// Inspect a nonempty local dataset without changing it or downloading data.
/// Known names must reach ninety percent of the published file size; SPK files
/// must also pass the native structural reader. These checks detect truncation,
/// not cryptographic identity or the scientific accuracy of all coefficients.
pub fn verify_assist_file(path: &Path, logical_name: &str) -> bool {
    let Ok(metadata) = std::fs::metadata(path) else {
        return false;
    };
    let size = metadata.len();
    if size == 0 {
        return false;
    }
    let expected = match logical_name {
        "linux_p1550p2650.440" => ASSIST_EXPECTED_SIZES[0],
        "linux_m13000p17000.441" => ASSIST_EXPECTED_SIZES[1],
        "sb441-n16.bsp" => ASSIST_EXPECTED_SIZES[2],
        _ => 0,
    };
    if size < (expected as f64 * ASSIST_MIN_SIZE_RATIO) as u64 {
        return false;
    }
    !logical_name.ends_with(".bsp") || crate::download::is_valid_bsp_file(path)
}
