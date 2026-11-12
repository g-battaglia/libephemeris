//! Configured symbolic lunar models and local fictitious-orbit datasets.

use super::{Runtime, Session};
use crate::errors::{Error, Result};
use crate::hypothetical::{self, HypoState, OrbitalElements};
use std::path::{Path, PathBuf};

impl Runtime {
    /// Uniform Waldemath position at TT with the configured Delta-T anchor.
    ///
    /// Longitude/latitude are degrees, radius AU and rates per day. The anchor
    /// UT1 epoch is converted with this runtime's selected model/override; its
    /// mean-ecliptic origin removes IAU 2006/2000A nutation at that TT epoch.
    /// Non-finite epochs, clock and nutation failures propagate.
    pub fn calc_waldemath(&self, jd_tt: f64) -> Result<HypoState> {
        let session = Session::from_context(self.effective_context()?);
        let anchor_tt = hypothetical::WALDEMATH_ANCHOR_JD_UT
            + session.deltat(hypothetical::WALDEMATH_ANCHOR_JD_UT)?;
        let (dpsi, _) = crate::erfa::nut06a(anchor_tt)?;
        hypothetical::waldemath_state_with_anchor(
            jd_tt,
            anchor_tt,
            hypothetical::WALDEMATH_ANCHOR_APPARENT_LON_DEG - dpsi.to_degrees(),
        )
    }

    /// Dispatch a hypothetical model at TT, in degrees/AU and daily rates.
    /// Waldemath uses the configured clock; other models use their published
    /// epoch/equinox conventions. Unsupported ids and model errors propagate.
    pub fn calc_hypothetical_position(&self, body_id: i32, jd_tt: f64) -> Result<HypoState> {
        if !hypothetical::is_hypothetical_body(body_id) {
            return Err(Error::InvalidBody {
                body_id,
                operation: "calc_hypothetical_position",
            });
        }
        if body_id == crate::constants::WALDEMATH {
            self.calc_waldemath(jd_tt)
        } else {
            hypothetical::calc_hypothetical_position(body_id, jd_tt)
        }
    }

    /// Locate the installed project-owned fictitious-orbit CSV dataset.
    /// The caller supplies its package directory. Missing/unreadable files
    /// return the actual filesystem error; no replacement is downloaded.
    pub fn bundled_fictitious_orbits_path(package_dir: &Path) -> Result<PathBuf> {
        let path = package_dir.join("data/fictitious_orbits.csv");
        std::fs::metadata(&path).map_err(|error| Error::EphemerisLoaderIo {
            path: path.display().to_string(),
            kind: error.kind(),
            errno: error.raw_os_error(),
            message: error.to_string(),
        })?;
        Ok(path)
    }

    /// Parse a user-selected text orbital-elements file into numeric records.
    /// UTF-8, syntax, epoch and physical-parameter errors propagate; a missing
    /// or inaccessible file retains its ordinary filesystem failure.
    pub fn parse_orbital_elements_file(path: &Path) -> Result<Vec<OrbitalElements>> {
        let text = std::fs::read_to_string(path).map_err(|error| Error::EphemerisLoaderIo {
            path: path.display().to_string(),
            kind: error.kind(),
            errno: error.raw_os_error(),
            message: error.to_string(),
        })?;
        hypothetical::parse_orbital_elements_text(&text)
    }

    /// Parse the installed CSV including its actual line-number metadata.
    /// The package directory is application supplied, with the same file and
    /// UTF-8 errors as [`Self::bundled_fictitious_orbits_path`].
    pub fn load_bundled_fictitious_orbits(package_dir: &Path) -> Result<Vec<OrbitalElements>> {
        let path = Self::bundled_fictitious_orbits_path(package_dir)?;
        let text = std::fs::read_to_string(&path).map_err(|error| Error::EphemerisLoaderIo {
            path: path.display().to_string(),
            kind: error.kind(),
            errno: error.raw_os_error(),
            message: error.to_string(),
        })?;
        hypothetical::parse_fictitious_csv_text(&text)
    }
}
