//! A continuous mean sidereal origin from the Vondrák mean poles.
//!
//! The celestial origin uses the mean equator and ecliptic of Vondrák,
//! Capitaine & Wallace (2011, A&A 534, A22; corrected coefficient tables),
//! and the nonrotating-origin connection of Capitaine et al. (2003,
//! A&A 412, 567). The J2000 datum and linear correction are the published
//! IAU 2006 precession-in-right-ascension convention, P03 equation 42.
//! No nutation or Delta-T model is evaluated here.

use crate::errors::Result;
use crate::positions::{
    ARCSEC_TO_RAD, VONDRAK_XYPER_PERIODS, VONDRAK_XYPER_X_COS, VONDRAK_XYPER_X_SIN,
    VONDRAK_XYPER_Y_COS, VONDRAK_XYPER_Y_SIN, VONDRAK_XYPOL, earth_rotation_angle_rad,
    ecliptic_pole_vondrak, equator_pole_vondrak,
};
use std::f64::consts::{PI, TAU};

use super::canonical_sidereal_table::{CIO_S, GAUSS16, NODE_START_CENTURIES, NODE_STEP_CENTURIES};

#[cfg(test)]
use super::canonical_sidereal_table::GAUSS32;

const J2000: f64 = 2_451_545.0;
const DAYS_PER_CENTURY: f64 = 36_525.0;
const J2000_DATUM_ARCSEC: f64 = 0.014_506;
const LINEAR_CORRECTION_ARCSEC: f64 = -0.003_872;

/// Greenwich mean sidereal angle, in degrees in `[0,360)`.
///
/// `jd_ut1` and `jd_tt` are independent Julian Days. UT1 drives only the
/// existing compensated Earth Rotation Angle. TT drives only the Vondrák
/// mean poles and their celestial origin. For TT centuries `T` from J2000,
///
/// ```text
/// s(T) = -integral_0^T (X Y' - Y X')/(1+Z) dT,
/// relative_EO(T) = s(T) - (phase(T) - phase(0)),
/// GMST = ERA(UT1) - relative_EO(T)
///        + (0.014506 - 0.003872 T) arcsec.
/// ```
///
/// Here `(X,Y,Z)` is the unit mean-equator pole, and `phase` is the angle
/// of its north-pole tangent origin in the mean equinox frame. The stored
/// locator at the nearest ten-century node plus a fixed 16-point Gauss
/// correction evaluates `s`. The equinox phase can wrap by a full turn;
/// the final sidereal angle is circular and has no such discontinuity.
///
/// The model comparison domain is `|T| <= 2000` (200,000 Julian years).
/// Outside it the same rule selects the nearest available endpoint and
/// applies the same fixed correction; neither model accuracy nor quadrature
/// convergence is qualified there. No date admission limit is introduced.
/// A nonfinite UT1 follows [`earth_rotation_angle_rad`]'s error contract;
/// undefined celestial arithmetic propagates IEEE NaN. This function uses
/// no heap allocation, ephemeris data, network, or mutable configuration.
pub(super) fn mean_degrees(jd_ut1: f64, jd_tt: f64) -> Result<f64> {
    mean_with_rule(jd_ut1, jd_tt, &GAUSS16)
}

/// The same mean angle using the fixed 32-point correction rule.
///
/// This independent-order comparison is a numerical control, not a
/// runtime adaptive choice or an alternative sidereal convention.
#[cfg(test)]
pub(super) fn mean_degrees_order32(jd_ut1: f64, jd_tt: f64) -> Result<f64> {
    mean_with_rule(jd_ut1, jd_tt, &GAUSS32)
}

fn mean_with_rule(jd_ut1: f64, jd_tt: f64, rule: &[(f64, f64)]) -> Result<f64> {
    let era = earth_rotation_angle_rad(jd_ut1)?;
    let t = (jd_tt - J2000) / DAYS_PER_CENTURY;
    let relative_eo = locator(t, rule) - (equinox_phase(jd_tt) - equinox_phase(J2000));
    let datum = (J2000_DATUM_ARCSEC + LINEAR_CORRECTION_ARCSEC * t) * ARCSEC_TO_RAD;
    let degrees = ((era - relative_eo + datum) * (180.0 / PI)).rem_euclid(360.0);
    // Rust's remainder can retain -0; callers receive the usual canonical zero.
    Ok(if degrees == 0.0 { 0.0 } else { degrees })
}

fn nearest_node(t: f64) -> usize {
    let position = (t - NODE_START_CENTURIES) / NODE_STEP_CENTURIES;
    // At an exact tie choose the lower index. Saturating casts also keep an
    // undefined TT from becoming an indexing failure; its arithmetic is NaN.
    let index = (position - 0.5).ceil().max(0.0) as usize;
    index.min(CIO_S.len() - 1)
}

fn locator(t: f64, rule: &[(f64, f64)]) -> f64 {
    let index = nearest_node(t);
    let node = NODE_START_CENTURIES + index as f64 * NODE_STEP_CENTURIES;
    let half = (t - node) * 0.5;
    let centre = node + half;
    let mut correction = 0.0;
    for &(abscissa, weight) in rule {
        correction += weight * locator_rate(centre + half * abscissa);
    }
    CIO_S[index] + half * correction
}

/// Values and analytic first derivatives in radians per TT century.
///
/// Coefficients and summation order are those of the engine's published
/// Vondrák equator series. Differentiating each sine/cosine and polynomial
/// term avoids both a finite difference and a round trip through Julian Days.
fn equator_xy_and_rate(t: f64) -> (f64, f64, f64, f64) {
    let (mut x, mut y, mut dx, mut dy) = (0.0, 0.0, 0.0, 0.0);
    for i in 0..VONDRAK_XYPER_PERIODS.len() {
        let angle = TAU * t / VONDRAK_XYPER_PERIODS[i];
        let (s, c) = angle.sin_cos();
        let frequency = TAU / VONDRAK_XYPER_PERIODS[i];
        x += c * VONDRAK_XYPER_X_COS[i] + s * VONDRAK_XYPER_X_SIN[i];
        y += c * VONDRAK_XYPER_Y_COS[i] + s * VONDRAK_XYPER_Y_SIN[i];
        dx += frequency * (-s * VONDRAK_XYPER_X_COS[i] + c * VONDRAK_XYPER_X_SIN[i]);
        dy += frequency * (-s * VONDRAK_XYPER_Y_COS[i] + c * VONDRAK_XYPER_Y_SIN[i]);
    }
    let mut power = 1.0;
    for &(px, py) in &VONDRAK_XYPOL {
        x += px * power;
        y += py * power;
        power *= t;
    }
    let mut power = 1.0;
    for (degree, &(px, py)) in VONDRAK_XYPOL.iter().enumerate().skip(1) {
        dx += degree as f64 * px * power;
        dy += degree as f64 * py * power;
        power *= t;
    }
    (
        x * ARCSEC_TO_RAD,
        y * ARCSEC_TO_RAD,
        dx * ARCSEC_TO_RAD,
        dy * ARCSEC_TO_RAD,
    )
}

fn locator_rate(t: f64) -> f64 {
    let (x, y, dx, dy) = equator_xy_and_rate(t);
    let z2 = 1.0 - x * x - y * y;
    let z = if z2 > 0.0 { z2.sqrt() } else { 0.0 };
    // For the published domain n=1 analytically. This equivalent normalized
    // expression also follows the engine's clamped pole outside that domain;
    // no derivative of its normalization appears in the antisymmetric term.
    let n = (x * x + y * y + z * z).sqrt();
    -(x * dy - y * dx) / (n * (n + z))
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn unit(v: [f64; 3]) -> [f64; 3] {
    let norm = dot(v, v).sqrt();
    [v[0] / norm, v[1] / norm, v[2] / norm]
}

fn equinox_phase(jd_tt: f64) -> f64 {
    let q = unit(equator_pole_vondrak(jd_tt));
    let e = unit(ecliptic_pole_vondrak(jd_tt));
    let x = unit(cross(q, e));
    let y = cross(q, x);
    let denominator = 1.0 + q[2];
    let a = [
        1.0 - q[0] * q[0] / denominator,
        -q[0] * q[1] / denominator,
        -q[0],
    ];
    dot(y, a).atan2(dot(x, a))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cyclic_degrees(a: f64, b: f64) -> f64 {
        ((a - b + 180.0).rem_euclid(360.0) - 180.0).abs()
    }

    #[test]
    fn analytic_equator_values_keep_the_engine_clock() {
        for t in [-2000.0, -200.0, -1.0, 0.0, 0.25, 200.0, 2000.0] {
            let tt = J2000 + t * DAYS_PER_CENTURY;
            let (x, y, _, _) = equator_xy_and_rate((tt - J2000) / DAYS_PER_CENTURY);
            let q = equator_pole_vondrak(tt);
            assert_eq!(x.to_bits(), q[0].to_bits());
            assert_eq!(y.to_bits(), q[1].to_bits());
        }
    }

    #[test]
    fn analytic_rates_match_symmetric_differences() {
        let step = 1.0 / 1024.0;
        for t in [-1000.0, -2.0, 0.0, 3.0, 1000.0] {
            let (_, _, dx, dy) = equator_xy_and_rate(t);
            let (xp, yp, _, _) = equator_xy_and_rate(t + step);
            let (xm, ym, _, _) = equator_xy_and_rate(t - step);
            assert!((dx - (xp - xm) / (2.0 * step)).abs() < 2e-12);
            assert!((dy - (yp - ym) / (2.0 * step)).abs() < 2e-12);
        }
    }

    #[test]
    fn j2000_is_the_published_datum() {
        let actual = mean_degrees(J2000, J2000).unwrap();
        let expected = (earth_rotation_angle_rad(J2000).unwrap() * (180.0 / PI)
            + J2000_DATUM_ARCSEC / 3600.0)
            .rem_euclid(360.0);
        assert!(cyclic_degrees(actual, expected) < 1e-12);
        assert_eq!(locator(0.0, &GAUSS16), 0.0);
    }

    #[test]
    fn independent_ut1_day_is_only_earth_rotation() {
        // Fifteen microarcseconds, expressed in degrees.
        let tolerance = 15e-6 / 3600.0;
        for t in [-2000.0, -10.0, 0.0, 10.0, 2000.0] {
            let tt = J2000 + t * DAYS_PER_CENTURY;
            let ut1 = tt + 0.125;
            let first = mean_degrees(ut1, tt).unwrap();
            let next = mean_degrees(ut1 + 1.0, tt).unwrap();
            let era_step = (earth_rotation_angle_rad(ut1 + 1.0).unwrap()
                - earth_rotation_angle_rad(ut1).unwrap())
                * (180.0 / PI);
            assert!(cyclic_degrees(next - first, era_step) <= tolerance);
        }
    }

    #[test]
    fn fixed_orders_agree_across_the_model_domain() {
        let tolerance = 1e-6 / 3600.0;
        for t in [
            -2000.0, -1997.5, -1005.0, -5.0, 0.0, 5.0, 1005.0, 1997.5, 2000.0,
        ] {
            let tt = J2000 + t * DAYS_PER_CENTURY;
            let a = mean_degrees(tt, tt).unwrap();
            let b = mean_degrees_order32(tt, tt).unwrap();
            assert!(cyclic_degrees(a, b) <= tolerance, "T={t}: {a} versus {b}");
        }
    }

    #[test]
    fn nearest_node_ties_and_joins_are_continuous() {
        for index in [0, 1, 199, 200, 399] {
            let t = NODE_START_CENTURIES + (index as f64 + 0.5) * NODE_STEP_CENTURIES;
            assert_eq!(nearest_node(t), index);
            let epsilon = 1e-9;
            // A change of node must agree after allowing the actual small
            // elapsed celestial connection, independently of Earth rotation.
            let left = locator(t - epsilon, &GAUSS16);
            let right = locator(t + epsilon, &GAUSS16);
            let expected_step = 2.0 * epsilon * locator_rate(t);
            assert!((right - left - expected_step).abs() < 5e-14);
        }
    }

    #[test]
    fn finite_dates_outside_model_domain_keep_the_same_rule() {
        assert_eq!(nearest_node(-1e8), 0);
        assert_eq!(nearest_node(1e8), CIO_S.len() - 1);
        for t in [-10_000.0, -2001.0, 2001.0, 10_000.0] {
            let tt = J2000 + t * DAYS_PER_CENTURY;
            assert!(mean_degrees(tt, tt).is_ok());
        }
        assert!(mean_degrees(J2000, f64::NAN).unwrap().is_nan());
    }

    /// Private, opt-in numerical probe bridge; never called by the hot path.
    ///
    /// Input is `{"clocks":[[ut1_jd,tt_jd],...]}`. Output preserves those
    /// clocks and reports both fixed quadrature orders in degrees. Running
    /// this test requires both named paths and does not read any prior result.
    #[test]
    #[ignore = "explicit GMST_CANONICAL_INPUT/OUTPUT numerical qualification"]
    fn qualification_dump() {
        let input = std::env::var_os("GMST_CANONICAL_INPUT").expect("GMST_CANONICAL_INPUT");
        let output = std::env::var_os("GMST_CANONICAL_OUTPUT").expect("GMST_CANONICAL_OUTPUT");
        let document: serde_json::Value =
            serde_json::from_reader(std::fs::File::open(input).unwrap()).unwrap();
        let clocks = document["clocks"].as_array().expect("clocks array");
        let mut rows = Vec::with_capacity(clocks.len());
        for clock in clocks {
            let pair = clock.as_array().expect("two-element clock pair");
            assert_eq!(pair.len(), 2);
            let ut1 = pair[0].as_f64().expect("UT1 number");
            let tt = pair[1].as_f64().expect("TT number");
            let t = (tt - J2000) / DAYS_PER_CENTURY;
            let delta = tt - ut1;
            assert_eq!(
                (ut1 + delta).to_bits(),
                tt.to_bits(),
                "public TT reconstruction"
            );
            let public_hours =
                super::super::public_sidereal::hours(ut1, 23.4, 0.0, &|_| Ok(delta)).unwrap();
            rows.push(serde_json::json!({
                "ut1": ut1,
                "tt": tt,
                "public_hours": public_hours,
                "centuries": t,
                "mean16_degrees": mean_degrees(ut1, tt).unwrap(),
                "mean32_degrees": mean_degrees_order32(ut1, tt).unwrap(),
                "locator16_radians": locator(t, &GAUSS16),
                "locator32_radians": locator(t, &GAUSS32),
                "relative_phase_radians": equinox_phase(tt) - equinox_phase(J2000),
            }));
        }
        serde_json::to_writer(
            std::io::BufWriter::new(std::fs::File::create(output).unwrap()),
            &serde_json::json!({"rows": rows}),
        )
        .unwrap();
    }
}
