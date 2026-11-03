//! Explicit modern IAU 2006 sidereal observables with independent clocks.

use crate::errors::{Error, Result};
use crate::positions::{compensated_dot3, earth_rotation_angle_rad};

/// Greenwich IAU 2006 mean or apparent sidereal time, in hours `[0,24)`.
///
/// `jd_ut1` drives the compensated IAU 2000 Earth Rotation Angle; `jd_tt`
/// independently drives precession/nutation. With `apparent=false`, add the
/// Capitaine–Wallace–Chapront precession-in-RA polynomial to ERA. Otherwise
/// subtract the equation of origins from the IAU 2006/2000A true frame and
/// its IERS `s06` CIO locator. No Delta-T model, ephemeris, network or session
/// state is acquired. This is a standard-model utility, not the VCW origin.
///
/// Its declared numerical comparison domain is TT JD
/// `2396757.5..=2469808.5`, the qualified modern interval and one-day clock
/// envelope, and UT1 within 200,000 Julian years of J2000. These limits apply
/// only to this explicit utility, not existing public sidereal entries.
/// Finite inputs outside that domain return [`Error::MathDomain`]; nonfinite
/// inputs return [`Error::NonFiniteArgument`]. Model uncertainty is distinct
/// from numerical accuracy; no remote P03 extrapolation is performed.
pub fn sidereal_time_iau2006_hours(jd_ut1: f64, jd_tt: f64, apparent: bool) -> Result<f64> {
    for value in [jd_ut1, jd_tt] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    if !(2_396_757.5..=2_469_808.5).contains(&jd_tt) {
        return Err(Error::MathDomain { value: jd_tt });
    }
    if !(2_451_545.0 - 73_050_000.0..=2_451_545.0 + 73_050_000.0).contains(&jd_ut1) {
        return Err(Error::MathDomain { value: jd_ut1 });
    }
    let era = earth_rotation_angle_rad(jd_ut1)?;
    let t = (jd_tt - 2_451_545.0) / 36_525.0;
    let angle = if apparent {
        let frame = crate::erfa::pnm06a(jd_tt)?;
        let [x, y, _] = frame[2];
        let locator = crate::erfa::cio_locator_s(t, x, y);
        era - equation_of_origins_radians(frame, locator)?
    } else {
        era + precession_in_ra_arcsec(t) * (std::f64::consts::PI / 648_000.0)
    };
    Ok(reduced_hours(angle))
}

fn reduced_hours(angle: f64) -> f64 {
    let hours = angle.rem_euclid(std::f64::consts::TAU) * (12.0 / std::f64::consts::PI);
    // The rounded endpoint is the same phase and must retain the half-open
    // public interval, including a tiny negative phase reduced to one turn.
    if hours >= 24.0 || hours == 0.0 {
        0.0
    } else {
        hours
    }
}

/// IAU 2006 precession in right ascension, arcseconds, at TT centuries.
pub(super) fn precession_in_ra_arcsec(t: f64) -> f64 {
    0.014506
        + t * (4612.156534
            + t * (1.3915817 + t * (-0.00000044 + t * (-0.000029956 + t * -0.0000000368))))
}

/// IERS equation of origins `s−atan2(y·a,x·a)`, in radians.
///
/// The rows are a proper equinox frame; its pole defines the nonsingular
/// north-chart meridian `a`. This shared geometry is independent of the
/// particular precession model. Undefined/nonfinite geometry returns an error.
pub(super) fn equation_of_origins_radians(frame: [[f64; 3]; 3], locator: f64) -> Result<f64> {
    let [x, y, z] = frame[2];
    let denominator = 1.0 + z;
    if !denominator.is_finite() || denominator <= 0.0 {
        return Err(Error::MathDomain { value: denominator });
    }
    let meridian = [1.0 - x * x / denominator, -x * y / denominator, -x];
    let p = compensated_dot3(frame[0], meridian);
    let q = compensated_dot3(frame[1], meridian);
    if p == 0.0 && q == 0.0 {
        return Err(Error::SingularMatrix);
    }
    let value = locator - q.atan2(p);
    if !value.is_finite() {
        return Err(Error::MathDomain { value });
    }
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hour_reduction_preserves_the_half_open_interval_at_rounded_endpoints() {
        for angle in [0.0, -0.0, -f64::from_bits(1), std::f64::consts::TAU] {
            assert_eq!(reduced_hours(angle).to_bits(), 0.0f64.to_bits());
        }
        for angle in [-1e-12, 1e-12, -1.0, 1.0, -100.0, 100.0] {
            assert!((0.0..24.0).contains(&reduced_hours(angle)));
        }
    }

    #[test]
    fn explicit_standard_matches_independent_erfa_mean_and_apparent_controls() {
        let cases = [
            (
                2_451_545.0,
                2_451_545.0,
                18.69737482696563,
                18.697138155539015,
            ),
            (
                2_396_758.5,
                2_396_758.5 + 64.0 / 86400.0,
                6.686144931370595,
                6.685992047171471,
            ),
            (
                2_469_807.5,
                2_469_807.5 + 64.0 / 86400.0,
                6.723047170931536,
                6.723304912637069,
            ),
            (
                2_453_736.5,
                2_453_736.5,
                6.700454827709062,
                6.700421083569526,
            ),
        ];
        for (ut1, tt, mean, apparent) in cases {
            assert!((sidereal_time_iau2006_hours(ut1, tt, false).unwrap() - mean).abs() < 1e-11);
            assert!((sidereal_time_iau2006_hours(ut1, tt, true).unwrap() - apparent).abs() < 1e-11);
        }
    }

    #[test]
    fn held_tt_consumes_only_independent_ut1_rotation() {
        for ut1 in [2_396_757.5, 2_451_545.0, 2_469_808.5] {
            for apparent in [false, true] {
                let a = sidereal_time_iau2006_hours(ut1, 2_451_545.0, apparent).unwrap();
                let b = sidereal_time_iau2006_hours(ut1 + 1.0, 2_451_545.0, apparent).unwrap();
                let advance = (b - a).rem_euclid(24.0);
                assert!((advance - 24.0 * 0.00273781191135448).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn explicit_domain_rejects_extrapolation_without_changing_existing_entries() {
        for tt in [2_396_757.5, 2_469_808.5] {
            for apparent in [false, true] {
                let value = sidereal_time_iau2006_hours(2_451_545.0, tt, apparent).unwrap();
                assert!((0.0..24.0).contains(&value));
            }
        }
        assert!(matches!(
            sidereal_time_iau2006_hours(2_451_545.0, 2_396_757.49, true),
            Err(Error::MathDomain { .. })
        ));
        assert!(matches!(
            sidereal_time_iau2006_hours(2_451_545.0, 2_469_808.51, false),
            Err(Error::MathDomain { .. })
        ));
        assert!(matches!(
            sidereal_time_iau2006_hours(2_451_545.0 + 73_050_001.0, 2_451_545.0, false),
            Err(Error::MathDomain { .. })
        ));
        assert!(matches!(
            sidereal_time_iau2006_hours(f64::NAN, 2_451_545.0, true),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            sidereal_time_iau2006_hours(2_451_545.0, f64::INFINITY, false),
            Err(Error::NonFiniteArgument { .. })
        ));
    }
}
