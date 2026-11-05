// SPDX-License-Identifier: AGPL-3.0
//! Explicit physical pole conventions for harmonic force models.
//!
//! Solar spin: IAU/WGCCRE pole RA286.13°, declination63.87° (Archinal et al.
//! 2018; Park et al. 2021, section3.4). The runtime Earth zonal model uses
//! the fixed ICRF north axis, preserving the exposed prescribed-force model;
//! it is not a claim that the instantaneous physical Earth pole is fixed.
//! A distinct dated terrestrial diagnostic follows Park et al.
//! 2021 equations20–26, the 18.6-year nutation and Vondrak mean precession.
//! The inertial realization here uses the canonical IERS frame bias, **not**
//! JPL's separately fitted DE pole-offset parameters. It is consequently an
//! explicitly named physical model, not a claim of exact JPL-program identity
//! or of reproducing an unspecified legacy default. Final qualification must
//! assess any resulting model differences with independent references.
//!
//! <https://ssd.jpl.nasa.gov/doc/Park.2021.AJ.DE440.pdf>

use super::catalogue::HarmonicPoles;
#[cfg(test)]
use super::state::Sum;
use crate::errors::{Error, Result};
use crate::types::TwoPartEpoch;

#[cfg(test)]
const J2000: f64 = 2451545.0;
#[cfg(test)]
const ARCSECOND: f64 = std::f64::consts::PI / (180.0 * 3600.0);

/// Explicit model conventions: one runtime prescription and a distinct dated
/// diagnostic. No automatic model switching or force-output fitting occurs.
#[derive(Clone, Copy, Debug)]
pub(super) enum PoleModel {
    /// ICRF-north terrestrial figure axis and IAU/WGCCRE solar figure axis.
    IcrfTerrestrialIauSolar,
    /// Separate dated terrestrial model with canonical IERS bias. It remains
    /// a diagnostic, not a silently substituted runtime force convention.
    #[cfg(test)]
    PublishedDateIers,
}

impl PoleModel {
    /// Unit pole directions in ICRF at two-part TDB JD. For the slow terrestrial
    /// orientation, TDB is used as the TT approximation. The published nutation
    /// is truncated, not a full Earth-orientation/IERS tidal realization. The
    /// precession series has its existing native long-term domain/uncertainty;
    /// nonfinite/overflowing angle or matrix evaluation is an explicit error.
    pub(super) fn at(self, epoch: [f64; 2]) -> Result<HarmonicPoles> {
        for value in epoch {
            if !value.is_finite() {
                return Err(Error::NonFiniteArgument { value });
            }
        }
        let normalized = TwoPartEpoch {
            hi: epoch[0],
            lo: epoch[1],
        }
        .add_offset(0.0)?;
        if !normalized.value().is_finite() {
            return Err(Error::MathDomain {
                value: normalized.value(),
            });
        }
        let sun = solar_pole();
        let earth = match self {
            Self::IcrfTerrestrialIauSolar => [0.0, 0.0, 1.0],
            #[cfg(test)]
            Self::PublishedDateIers => {
                let interval = normalized.add_offset(-J2000)?;
                let centuries = interval.value() / 36525.0;
                let datum = terrestrial_datum_pole(centuries)?;
                let matrix = crate::positions::vondrak_precession_matrix(normalized.value(), true)?;
                // Native P maps ICRF column vectors to mean date. Its transpose
                // brings the published nutated pole-of-date back to ICRF.
                std::array::from_fn(|axis| {
                    let mut value = Sum::default();
                    for row in 0..3 {
                        value.add(matrix[row][axis] * datum[row]);
                    }
                    value.value()
                })
            }
        };
        if !sun.into_iter().chain(earth).all(f64::is_finite) {
            return Err(Error::MathDomain { value: f64::NAN });
        }
        Ok(HarmonicPoles { sun, earth })
    }
}

fn solar_pole() -> [f64; 3] {
    let (sin_ra, cos_ra) = 286.13_f64.to_radians().sin_cos();
    let (sin_dec, cos_dec) = 63.87_f64.to_radians().sin_cos();
    [cos_dec * cos_ra, cos_dec * sin_ra, sin_dec]
}

/// Published nutated terrestrial pole expressed in mean-equatorial date axes.
/// The angles are the explicit arcsecond polynomials in Park2021, not a table
/// of coefficients extracted from any software implementation.
#[cfg(test)]
fn terrestrial_datum_pole(t: f64) -> Result<[f64; 3]> {
    let node = (((0.008 * t + 7.455) * t - 6962890.539) * t + 450160.280) * ARCSECOND;
    let obliquity = (((0.001813 * t - 0.00059) * t - 46.815) * t + 84381.448) * ARCSECOND;
    for value in [node, obliquity] {
        if !value.is_finite() {
            return Err(Error::MathDomain { value });
        }
    }
    let (sin_node, cos_node) = node.sin_cos();
    let longitude = -17.1996 * ARCSECOND * sin_node;
    let delta_obliquity = 9.2025 * ARCSECOND * cos_node;
    let (sin_mean, cos_mean) = obliquity.sin_cos();
    let (sin_true, cos_true) = (obliquity + delta_obliquity).sin_cos();
    let (sin_longitude, cos_longitude) = longitude.sin_cos();
    Ok([
        sin_longitude * sin_true,
        cos_mean * cos_longitude * sin_true - sin_mean * cos_true,
        sin_mean * cos_longitude * sin_true + cos_mean * cos_true,
    ])
}

#[cfg(test)]
#[path = "poles_tests.rs"]
mod tests;
