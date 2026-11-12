//! Scalar public sidereal acquisition with independently consumed model clocks.
//!
//! Inputs are UT1 Julian days; the borrowed callback supplies TT-minus-UT1.
//! Exactly the requested clock is consumed. ERA uses UT1 alone; the continuous
//! VCW celestial origin uses TT alone, with no automatic boundary calibration.
//!
//! This public acquisition preserves IEEE NaN propagation. The checked stored
//! and house kernels keep their separate finite-input contracts. Infinite
//! arguments to trigonometric series and zero precession norms remain errors.
use super::{SIDEREAL_J2000_JD, ecliptic_pole_vondrak, equator_pole_vondrak};
use crate::errors::{Error, Result};

/// Floating remainder reduced to a nonnegative angle for a positive period.
/// Infinite or NaN dividends retain the IEEE remainder result; no finite
/// normalization is used to invent a number. Divisors here are 360 or 24.
fn positive_mod(value: f64, divisor: f64) -> f64 {
    let remainder = value % divisor;
    if remainder == 0.0 {
        0.0
    } else if remainder < 0.0 {
        remainder + divisor
    } else {
        remainder
    }
}

/// Report a consumed trigonometric domain failure as an ordinary core error.
fn trigonometric_error(value: f64) -> Error {
    Error::MathDomain { value }
}

/// Refuse an infinite Vondrák periodic-series argument at its consumed stage.
/// NaN is permitted because sine, cosine and inverse cosine propagate it.
fn pole_domain(tt: f64) -> Result<()> {
    let angle = 2.0 * std::f64::consts::PI * ((tt - SIDEREAL_J2000_JD) / 36525.0);
    if angle.is_infinite() {
        Err(trigonometric_error(angle))
    } else {
        Ok(())
    }
}

/// Mean equator/ecliptic pole angle in radians at scalar TT Julian Day.
/// Explicit comparison clamps preserve NaN, while rounded finite dot products
/// beyond the cosine interval are clamped to the corresponding endpoint.
fn obliquity(tt: f64) -> Result<f64> {
    pole_domain(tt)?;
    let e = ecliptic_pole_vondrak(tt);
    let q = equator_pole_vondrak(tt);
    let d = (e[0] * q[0] + e[1] * q[1] + e[2] * q[2]).clamp(-1.0, 1.0);
    Ok(d.acos())
}

/// Canonical mean degrees using exactly the requested model clock.
/// UT1 drives ERA alone and TT the continuous mean VCW origin. No boundary
/// clock, calibration or stored offset is acquired. Provider errors propagate
/// at the requested epoch; IEEE NaN and consumed trigonometric-domain behavior
/// retain the public scalar contract. Borrowed callbacks never escape.
pub(crate) fn mean_deg(ut: f64, delta: &dyn Fn(f64) -> Result<f64>) -> Result<f64> {
    let tt = ut + delta(ut)?;
    pole_domain(tt)?;
    if ut.is_nan() {
        return Ok(f64::NAN);
    }
    super::canonical_sidereal::mean_degrees(ut, tt)
}

/// Public Greenwich apparent sidereal hours with caller-supplied angles.
/// `eps` and `nut` are degrees. The radians/degrees round-trip precedes the
/// equation of equinoxes. Infinite obliquity fails at cosine after model-clock
/// acquisition; NaN angles and model offsets propagate through modulo 24.
pub(crate) fn hours(
    ut: f64,
    eps: f64,
    nut: f64,
    delta: &dyn Fn(f64) -> Result<f64>,
) -> Result<f64> {
    hours_with_longitude(ut, 0.0, eps, nut, delta)
}

/// Local apparent sidereal hours with east-positive longitude in degrees.
/// Longitude/15 is added to the unreduced Greenwich apparent hours before
/// modulo 24. Clock acquisition, equation of equinoxes and raw IEEE domains
/// otherwise follow the Greenwich helper; no extra boundary read occurs.
pub(crate) fn hours_with_longitude(
    ut: f64,
    longitude: f64,
    eps: f64,
    nut: f64,
    delta: &dyn Fn(f64) -> Result<f64>,
) -> Result<f64> {
    let gmst = mean_deg(ut, delta)?.to_radians().to_degrees() / 15.0;
    let eps_rad = eps.to_radians();
    if eps_rad.is_infinite() {
        return Err(trigonometric_error(eps_rad));
    }
    let eoe = nut * eps_rad.cos() / 15.0;
    let gast = gmst + eoe;
    let last = gast + longitude / 15.0;
    let mut last = positive_mod(last, 24.0);
    if last < 0.0 {
        last += 24.0;
    }
    Ok(last)
}

/// Configured scalar House ARMC and true obliquity at original UT1 JD.
/// The central model TT feeds raw nutation then Vondrak mean obliquity before
/// the separately consumed Greenwich mean-sidereal model clocks. Return
/// `(ARMC0, true_obliquity)` in degrees. No observer, catalog or image is read.
/// Canonical mean time consumes only its requested central clock at every
/// epoch, without boundary calibration. Raw NaN remains IEEE;
/// infinite pole/cosine arguments retain their owned builtin primitive role.
/// This private acquisition leaves existing checked House APIs unchanged.
pub(crate) fn house_frame_at(ut: f64, delta: &dyn Fn(f64) -> Result<f64>) -> Result<(f64, f64)> {
    let tt = ut + delta(ut)?;
    let (dpsi, deps) = crate::erfa::source_nut06a(tt);
    let eps = obliquity(tt)?.to_degrees() + deps.to_degrees();
    let mean = mean_deg(ut, delta)?;
    let eps_rad = eps.to_radians();
    if eps_rad.is_infinite() {
        return Err(trigonometric_error(eps_rad));
    }
    let armc = positive_mod(mean + dpsi.to_degrees() * eps_rad.cos() + 0.0, 360.0);
    Ok((armc, eps))
}

#[cfg(test)]
mod tests {
    use super::super::{SIDEREAL_MODERN_END_JD, SIDEREAL_MODERN_START_JD};
    use super::*;
    use std::cell::RefCell;
    #[test]
    fn modern_consumes_only_the_requested_clock() {
        let calls = RefCell::new(Vec::new());
        let result = mean_deg(2451545.0, &|jd| {
            calls.borrow_mut().push(jd);
            if jd != 2451545.0 {
                return Err(Error::DeltaTModelUnavailable {
                    model: "unused-boundary",
                });
            }
            Ok(0.001)
        })
        .unwrap();
        assert!(result.is_finite());
        assert_eq!(*calls.borrow(), vec![2451545.0]);
    }
    #[test]
    fn outside_and_exact_boundaries_consume_only_the_requested_clock() {
        for (jd, _boundary) in [
            (2390000.0, SIDEREAL_MODERN_START_JD),
            (SIDEREAL_MODERN_START_JD, SIDEREAL_MODERN_START_JD),
            (SIDEREAL_MODERN_END_JD, SIDEREAL_MODERN_END_JD),
            (2470000.0, SIDEREAL_MODERN_END_JD),
        ] {
            let calls = RefCell::new(Vec::new());
            let value = mean_deg(jd, &|x| {
                calls.borrow_mut().push(x);
                Ok(0.001)
            })
            .unwrap();
            assert!(value.is_finite());
            assert_eq!(*calls.borrow(), vec![jd]);
        }
    }
    #[test]
    fn requested_clock_failure_precedes_any_boundary() {
        let calls = RefCell::new(Vec::new());
        let error = mean_deg(2390000.0, &|x| {
            calls.borrow_mut().push(x);
            Err(Error::DeltaTModelUnavailable { model: "request" })
        })
        .unwrap_err();
        assert!(matches!(&error, Error::DeltaTModelUnavailable { .. }));
        assert_eq!(*calls.borrow(), vec![2390000.0]);
    }
    #[test]
    fn unused_boundary_failure_cannot_enter_the_canonical_clock() {
        let requested = 2390000.0;
        let calls = RefCell::new(Vec::new());
        let value = mean_deg(requested, &|x| {
            calls.borrow_mut().push(x);
            if x != requested {
                return Err(Error::DeltaTModelUnavailable {
                    model: "unused-boundary",
                });
            }
            Ok(0.001)
        })
        .unwrap();
        assert!(value.is_finite());
        assert_eq!(*calls.borrow(), vec![requested]);
    }
    #[test]
    fn nan_model_propagates_in_each_selected_branch() {
        for jd in [2390000.0, 2451545.0, 2470000.0] {
            assert!(mean_deg(jd, &|_| Ok(f64::NAN)).unwrap().is_nan());
            assert!(hours(jd, 23.4, 0.001, &|_| Ok(f64::NAN)).unwrap().is_nan());
        }
    }
    #[test]
    fn infinite_obliquity_refuses_after_the_consumed_clock() {
        let calls = RefCell::new(Vec::new());
        assert!(matches!(
            hours(2451545.0, f64::INFINITY, 0.0, &|x| {
                calls.borrow_mut().push(x);
                Ok(0.001)
            }),
            Err(error) if matches!(error.primary_error(), Error::MathDomain { .. })
        ));
        assert_eq!(*calls.borrow(), vec![2451545.0]);
        assert!(
            hours(2451545.0, f64::NAN, 0.0, &|_| Ok(0.001))
                .unwrap()
                .is_nan()
        );
    }
    #[test]
    fn finite_shared_realization_preserves_checked_kernel_words() {
        let (lo, hi) = super::super::sidereal_continuity_offsets(0.001, 0.001).unwrap();
        let dt = super::super::SiderealDeltaT::try_new(0.001, lo, hi).unwrap();
        for jd in [2390000.0, 2451545.0, 2470000.0] {
            assert_eq!(
                mean_deg(jd, &|_| Ok(0.001)).unwrap().to_bits(),
                super::super::mean_sidereal_time_deg(jd, &dt)
                    .unwrap()
                    .to_bits()
            );
        }
    }
}

#[cfg(test)]
mod domain_tests {
    use super::*;
    use std::cell::Cell;

    #[test]
    fn infinite_obliquity_reports_the_consumed_math_domain() {
        for input in [f64::INFINITY, f64::NEG_INFINITY] {
            let calls = Cell::new(0);
            let error = hours(2451545.0, input, 0.0, &|_| {
                calls.set(calls.get() + 1);
                Ok(0.001)
            })
            .unwrap_err();
            assert_eq!(calls.get(), 1);
            assert!(matches!(error, Error::MathDomain { value } if value == input));
        }
    }

    #[test]
    fn infinite_periodic_argument_is_a_native_error() {
        assert!(matches!(
            pole_domain(f64::INFINITY),
            Err(Error::MathDomain { .. })
        ));
    }

    #[test]
    fn model_failure_propagates_without_a_presentation_layer() {
        let error = hours(2451545.0, f64::INFINITY, 0.0, &|_| {
            Err(Error::DeltaTModelUnavailable {
                model: "unavailable",
            })
        })
        .unwrap_err();
        assert!(matches!(error, Error::DeltaTModelUnavailable { .. }));
    }
}
