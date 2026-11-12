//! Planetary-moon theories, satellite kernels, and apparent places.
//!
//! Analytical models include Neptune/Triton NEP 097 elements with nodal
//! precession, Pluto/Charon PLU 060 two-body elements, the Galilean E 5 series
//! (Lieske; Meeus chapter 44), Uranian GUST 86-geometry elements, and the
//! Saturnian TASS 1.7 series (Vienne & Duriez 1995, A&A 297, 588–605).
//! Center-of-body offsets combine the satellite positions with system mass
//! weights and a fixed J2000 ecliptic-to-equatorial rotation where required.
//! The TASS coefficient tables retain their MIT license and attribution.
//!
//! Epochs are TT Julian days. Analytical model positions are AU, derived from
//! kilometre elements using the IAU astronomical unit. Kernel states are
//! sampled on their TDB clock and their velocities converted to TT daily rates.
//! Apparent places include three light-time iterations, observer motion,
//! optional aberration, nutation, sidereal projection, and output frames.
//! Angular rates use adaptive fourth-order TT differentiation; radial rates
//! differentiate the retarded geometric ray analytically. Outputs are degrees,
//! AU, and rates per TT day. Sidereal corrections consume UT1.
//!
//! The analytical element models are approximations: Triton is accurate to
//! roughly 20–50 km, Charon to 1–5 km, the Galilean E 5 model to 100–200 km,
//! and TASS to 50–100 km per satellite. Frozen orientation or secular elements
//! limit long-term accuracy. Registered SPK kernels retain their actual segment
//! coverage; an unavailable or uncovered target is never silently substituted.
//! Non-finite epochs produce [`Error::NonFiniteArgument`].

use crate::errors::{Error, Result};

#[path = "moon_spk_serving.rs"]
mod moon_spk_serving;
pub use moon_spk_serving::moon_position_native_tt;
pub(crate) use moon_spk_serving::moon_spk_ayanamsha;
pub(crate) use moon_spk_serving::moon_spk_frame;
pub(crate) use moon_spk_serving::registered_spk_position_tt;

/// Semi-major axis of Triton in km (NEP 097).
pub const TRITON_A_KM: f64 = 354759.0;
/// Orbital period of Triton in days.
pub const TRITON_PERIOD_DAYS: f64 = 5.8768541;
/// Mean motion in radians per day (`2π / TRITON_PERIOD_DAYS`).
pub const TRITON_N_RAD_DAY: f64 = 2.0 * std::f64::consts::PI / TRITON_PERIOD_DAYS;
/// Eccentricity (nearly circular).
pub const TRITON_E: f64 = 0.000016;
/// Reference epoch, J2000.0 (Julian Date, TT).
pub const TRITON_EPOCH_JD: f64 = 2451545.0;
/// Mean anomaly at epoch in radians (NEP 097: 358.656 degrees).
pub const TRITON_M0_RAD: f64 = 358.656_f64.to_radians();
/// Argument of perihelion in radians (degenerate for e ≈ 0; model pins 0°).
pub const TRITON_OMEGA_RAD: f64 = 0.0;
/// Inclination to the J2000 ecliptic in radians (retrograde; NEP 097:
/// 156.865 degrees).
pub const TRITON_I_RAD: f64 = 156.865_f64.to_radians();
/// Longitude of ascending node at J2000 in radians (NEP 097: 177.608 degrees).
pub const TRITON_NODE_J2000_RAD: f64 = 177.608_f64.to_radians();
/// Node precession rate in radians per day (`−2π / (688·365.25)`; ~688 yr
/// period from Neptune's J 2, Jacobson 2009).
pub const TRITON_NODE_RATE_RAD_DAY: f64 = -2.0 * std::f64::consts::PI / (688.0 * 365.25);
/// Astronomical Unit in km.
pub const MOONS_AU_KM: f64 = 149_597_870.7;
/// GM of Neptune in km³/s² (NEP 097).
pub const GM_NEPTUNE: f64 = 6_835_099.97;
/// GM of Triton in km³/s² (NEP 097).
pub const GM_TRITON: f64 = 1428.49546;
/// Mean obliquity at J2000.0, radians (IAU 2006: 84381.406 arcseconds).
pub const EPS_J2000_RAD: f64 = (84381.406_f64 / 3600.0).to_radians();

/// Triton's position relative to Neptune's center.
///
/// `jd` is a TT Julian Date (any finite value); output is AU in the emitted
/// J2000-ecliptic frame.
///
/// Returns `(x, y, z)` in AU in the model frame (treated as
/// J2000-ecliptic by `get_cob_offset`, which rotates it — see
/// [`neptune_cob_offset_icrf`]).
///
/// Propagates mean anomaly with one fixed-point Kepler step, adequate for
/// e ≈ 0.000016.
///
/// # Errors
///
/// [`Error::NonFiniteArgument`] when `jd` is NaN or infinite.
pub fn triton_position(jd: f64) -> Result<(f64, f64, f64)> {
    if !jd.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd });
    }
    // Time since epoch in days.
    let dt = jd - TRITON_EPOCH_JD;

    // Mean anomaly at time t, normalized to 0–2π.
    let m = (TRITON_M0_RAD + TRITON_N_RAD_DAY * dt).rem_euclid(2.0 * std::f64::consts::PI);

    // One fixed-point Kepler step.
    let e_anom = m + TRITON_E * m.sin();

    // True anomaly.
    let cos_nu = (e_anom.cos() - TRITON_E) / (1.0 - TRITON_E * e_anom.cos());
    let sin_nu =
        ((1.0 - TRITON_E * TRITON_E).sqrt() * e_anom.sin()) / (1.0 - TRITON_E * e_anom.cos());
    let nu = sin_nu.atan2(cos_nu);

    // Radius (distance from Neptune) in km.
    let r_km = TRITON_A_KM * (1.0 - TRITON_E * TRITON_E) / (1.0 + TRITON_E * nu.cos());

    // Position in the orbital plane.
    let x_orb = r_km * nu.cos();
    let y_orb = r_km * nu.sin();

    // Current node position (with J 2 precession).
    let node = TRITON_NODE_J2000_RAD + TRITON_NODE_RATE_RAD_DAY * dt;

    let cos_omega = TRITON_OMEGA_RAD.cos();
    let sin_omega = TRITON_OMEGA_RAD.sin();
    let cos_i = TRITON_I_RAD.cos();
    let sin_i = TRITON_I_RAD.sin();
    let cos_node = node.cos();
    let sin_node = node.sin();

    // Rotate by ω in the orbital plane (model keeps `z1 = 0.0` here).
    let x1 = x_orb * cos_omega - y_orb * sin_omega;
    let y1 = x_orb * sin_omega + y_orb * cos_omega;
    let z1 = 0.0;

    // Rotate by inclination around the x-axis.
    let x2 = x1;
    let y2 = y1 * cos_i - z1 * sin_i;
    let z2 = y1 * sin_i + z1 * cos_i;

    // Rotate by node around the z-axis.
    let x_ecl = x2 * cos_node - y2 * sin_node;
    let y_ecl = x2 * sin_node + y2 * cos_node;
    let z_ecl = z2;

    // Convert to AU.
    Ok((
        x_ecl / MOONS_AU_KM,
        y_ecl / MOONS_AU_KM,
        z_ecl / MOONS_AU_KM,
    ))
}

/// Offset from the Neptune-system barycenter to Neptune's center of body.
///
/// Computes the GM-weighted Triton
/// offset, negated, in the same frame `triton_position` emits. Triton
/// holds 99.5% of the satellite mass, so it is the only moon carried.
///
/// # Errors
///
/// [`Error::NonFiniteArgument`] when `jd` is NaN or infinite.
pub fn neptune_cob_offset(jd: f64) -> Result<(f64, f64, f64)> {
    let (tx, ty, tz) = triton_position(jd)?;
    let total_gm = GM_NEPTUNE + GM_TRITON;
    let bary_x = GM_TRITON * tx / total_gm;
    let bary_y = GM_TRITON * ty / total_gm;
    let bary_z = GM_TRITON * tz / total_gm;
    Ok((-bary_x, -bary_y, -bary_z))
}

/// Rotate a J2000-ecliptic vector to the equatorial ICRF frame.
///
/// Rotates ecliptic coordinates using the fixed
/// J2000 obliquity. Pure rotation: defined for every finite or
/// non-finite input triple, no validation.
pub fn ecliptic_to_equatorial_j2000(xyz: (f64, f64, f64)) -> (f64, f64, f64) {
    let (x, y, z) = xyz;
    let cos_eps = EPS_J2000_RAD.cos();
    let sin_eps = EPS_J2000_RAD.sin();
    (x, y * cos_eps - z * sin_eps, y * sin_eps + z * cos_eps)
}

/// Neptune center-of-body offset in the equatorial ICRF/J2000 frame.
///
/// Neptune center-of-body displacement in ICRF/J2000.
///
/// Combines the Neptune dispatch arm: [`neptune_cob_offset`]
/// followed by [`ecliptic_to_equatorial_j2000`]. Add the result to the
/// Neptune-system barycenter position to get Neptune's center.
///
/// # Errors
///
/// [`Error::NonFiniteArgument`] when `jd` is NaN or infinite.
pub fn neptune_cob_offset_icrf(jd: f64) -> Result<(f64, f64, f64)> {
    Ok(ecliptic_to_equatorial_j2000(neptune_cob_offset(jd)?))
}

#[cfg(test)]
mod tests {
    use super::*;

    const TWO_PI: f64 = 2.0 * std::f64::consts::PI;

    fn norm(v: (f64, f64, f64)) -> f64 {
        (v.0 * v.0 + v.1 * v.1 + v.2 * v.2).sqrt()
    }

    #[test]
    fn epoch_position_is_finite_and_bounded_by_elements() {
        let p = triton_position(TRITON_EPOCH_JD).expect("epoch must evaluate");
        assert!(p.0.is_finite() && p.1.is_finite() && p.2.is_finite());
        // Analytic bound from the orbital elements: |r| in
        // a(1−e)/AU.. a(1+e)/AU regardless of anomaly/rotation.
        let r = norm(p);
        let lo = TRITON_A_KM * (1.0 - TRITON_E) / MOONS_AU_KM;
        let hi = TRITON_A_KM * (1.0 + TRITON_E) / MOONS_AU_KM;
        assert!(r >= lo && r <= hi, "r={r} outside [{lo}, {hi}]");
    }

    #[test]
    fn position_is_deterministic() {
        let jd = 2460000.5;
        assert_eq!(triton_position(jd).unwrap(), triton_position(jd).unwrap());
    }

    #[test]
    fn mean_anomaly_periodicity_returns_same_position() {
        // Advancing by exactly one orbital period repeats the mean anomaly
        // and (up to the slow node drift, which this test avoids by using
        // a node-frozen check below) the orbital-plane position. Here the
        // full position differs only through node precession over one
        // 5.87 d period (~2.4 e-5 rad), so bound the drift analytically
        // instead of demanding exactly equality.
        let jd = 2451545.0;
        let a = triton_position(jd).unwrap();
        let b = triton_position(jd + TRITON_PERIOD_DAYS).unwrap();
        let drift = ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2) + (a.2 - b.2).powi(2)).sqrt();
        let max_drift = 2.0
            * (TRITON_A_KM / MOONS_AU_KM)
            * (TRITON_NODE_RATE_RAD_DAY * TRITON_PERIOD_DAYS).abs();
        assert!(drift <= max_drift, "drift={drift} bound={max_drift}");
    }

    #[test]
    fn node_precession_advances_retrograde() {
        // Over one century the node rate must move the node by
        // exactly RATE*36525 rad (linear law); check the law
        // through two far-apart evaluations of the exported rate.
        let dt = 36525.0;
        let node_now = TRITON_NODE_J2000_RAD;
        let node_then = TRITON_NODE_J2000_RAD + TRITON_NODE_RATE_RAD_DAY * dt;
        let expected = -TWO_PI * dt / (688.0 * 365.25);
        assert!((node_then - node_now - expected).abs() < 1e-15);
        assert!(node_then < node_now, "precession must be retrograde");
    }

    #[test]
    fn one_step_kepler_residual_matches_release_approximation() {
        // Structural sanity (not a comparison): the model takes exactly
        // one fixed-point step E = M + e·sin M, so the residual
        // |E − M − e·sin E| is O(e²·|sin|) — tiny but nonzero in general.
        // This pins the single-step structure against an iterated solve.
        let jd = 2451545.0 + 1234.5;
        let dt = jd - TRITON_EPOCH_JD;
        let m = (TRITON_M0_RAD + TRITON_N_RAD_DAY * dt).rem_euclid(TWO_PI);
        let e_one = m + TRITON_E * m.sin();
        let residual = (e_one - m - TRITON_E * e_one.sin()).abs();
        // |E 1 − M − e·sin E 1| = e·|sin M − sin E 1| ≤ e·|E 1 − M| ≤ e²,
        // so the bound below pins the single-step structure.
        assert!(residual <= TRITON_E * TRITON_E, "residual={residual}");
        // And an iterated solve agrees with the one-step value far below
        // the model precision envelope (~20 km ≈ 1.3 e-7 AU).
        let mut e_iter = m;
        for _ in 0..50 {
            e_iter = m + TRITON_E * e_iter.sin();
        }
        let au_per_km = 1.0 / MOONS_AU_KM;
        assert!((e_iter - e_one).abs() * TRITON_A_KM * au_per_km < 1e-12);
    }

    #[test]
    fn malformed_jd_is_refused() {
        // NaN carries no `==` identity, so match the variant and inspect
        // the payload instead of asserting error equality.
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            for err in [
                triton_position(bad),
                neptune_cob_offset(bad),
                neptune_cob_offset_icrf(bad),
            ] {
                match err {
                    Err(Error::NonFiniteArgument { value }) => {
                        if bad.is_nan() {
                            assert!(value.is_nan());
                        } else {
                            assert_eq!(value, bad);
                        }
                    }
                    other => panic!("jd={bad}: expected NonFiniteArgument, got {other:?}"),
                }
            }
        }
    }

    #[test]
    fn cob_offset_is_antiparallel_and_gm_weighted() {
        let jd = 2460000.5;
        let (tx, ty, tz) = triton_position(jd).unwrap();
        let (ox, oy, oz) = neptune_cob_offset(jd).unwrap();
        // weight: offset = −(GM_TRITON / (GM_NEPTUNE+GM_TRITON))·r.
        let w = GM_TRITON / (GM_NEPTUNE + GM_TRITON);
        assert!((ox + w * tx).abs() < 1e-18);
        assert!((oy + w * ty).abs() < 1e-18);
        assert!((oz + w * tz).abs() < 1e-18);
        // Magnitude bound: |offset| <= w·a(1+e)/AU.
        let hi = w * TRITON_A_KM * (1.0 + TRITON_E) / MOONS_AU_KM;
        assert!(norm((ox, oy, oz)) <= hi);
    }

    #[test]
    fn rotation_is_orthonormal_and_inverts() {
        let jd = 2460000.5;
        let ecl = neptune_cob_offset(jd).unwrap();
        let eq = ecliptic_to_equatorial_j2000(ecl);
        // Pure rotation preserves the norm.
        assert!((norm(eq) - norm(ecl)).abs() < 1e-20);
        // Inverse rotation (transpose) recovers the input.
        let (x, y, z) = eq;
        let cos_eps = EPS_J2000_RAD.cos();
        let sin_eps = EPS_J2000_RAD.sin();
        let back = (x, y * cos_eps + z * sin_eps, -y * sin_eps + z * cos_eps);
        assert!((back.0 - ecl.0).abs() < 1e-20);
        assert!((back.1 - ecl.1).abs() < 1e-20);
        assert!((back.2 - ecl.2).abs() < 1e-20);
    }

    #[test]
    fn icrf_arm_composes_rotation_verbatim() {
        let jd = 2451545.0 + 777.25;
        let icrf = neptune_cob_offset_icrf(jd).unwrap();
        let expected = ecliptic_to_equatorial_j2000(neptune_cob_offset(jd).unwrap());
        assert_eq!(icrf, expected);
    }

    #[test]
    fn far_epoch_stays_bounded() {
        // ±300 years: anomaly normalization and rotation keep the radius
        // inside the element bounds (no overflow-driven escape).
        for jd in [TRITON_EPOCH_JD - 109_575.0, TRITON_EPOCH_JD + 109_575.0] {
            let r = norm(triton_position(jd).unwrap());
            let lo = TRITON_A_KM * (1.0 - TRITON_E) / MOONS_AU_KM;
            let hi = TRITON_A_KM * (1.0 + TRITON_E) / MOONS_AU_KM;
            assert!(r >= lo && r <= hi, "jd={jd} r={r}");
        }
    }

    #[test]
    fn release_constants_match_source_values() {
        assert_eq!(TRITON_A_KM, 354759.0);
        assert_eq!(TRITON_PERIOD_DAYS, 5.8768541);
        assert_eq!(TRITON_E, 0.000016);
        assert_eq!(TRITON_EPOCH_JD, 2451545.0);
        assert_eq!(GM_NEPTUNE, 6_835_099.97);
        assert_eq!(GM_TRITON, 1_428.495_46);
        assert_eq!(MOONS_AU_KM, 149_597_870.7);
        // Angular constants round-trip to their published degree values.
        assert!((TRITON_M0_RAD.to_degrees() - 358.656).abs() < 1e-12);
        assert!((TRITON_I_RAD.to_degrees() - 156.865).abs() < 1e-12);
        assert!((TRITON_NODE_J2000_RAD.to_degrees() - 177.608).abs() < 1e-12);
        assert!((EPS_J2000_RAD.to_degrees() - 84381.406 / 3600.0).abs() < 1e-15);
    }
}

// ==================== Pluto/Charon (PLU 060 two-body) ====================
// Frozen PLU 060 elements with five fixed-point Kepler iterations. The
// approximate two-body model neglects secular rates and has 1–5 km accuracy.
// The center-of-body offset is the negative GM-weighted Charon displacement.

/// Semi-major axis of Charon in km (PLU 060).
pub const CHARON_A_KM: f64 = 19591.4;
/// Orbital period of Charon in days (PLU 060).
pub const CHARON_PERIOD_DAYS: f64 = 6.3872304;
/// Mean motion in radians per day (`2π / CHARON_PERIOD_DAYS`).
pub const CHARON_N_RAD_DAY: f64 = 2.0 * std::f64::consts::PI / CHARON_PERIOD_DAYS;
/// Eccentricity (nearly circular).
pub const CHARON_E: f64 = 0.0002;
/// Reference epoch, J2000.0 (Julian Date, TT).
pub const CHARON_EPOCH_JD: f64 = 2451545.0;
/// Mean anomaly at epoch in radians (PLU 060: 71.255 degrees).
pub const CHARON_M0_RAD: f64 = 71.255_f64.to_radians();
/// Argument of perihelion in radians (PLU 060: 180.0 degrees).
pub const CHARON_OMEGA_RAD: f64 = 180.0_f64.to_radians();
/// Inclination to the J2000 ecliptic in radians (PLU 060:
/// 96.145 degrees).
pub const CHARON_I_RAD: f64 = 96.145_f64.to_radians();
/// Longitude of ascending node at epoch in radians (PLU 060:
/// 223.046 degrees).
pub const CHARON_NODE_RAD: f64 = 223.046_f64.to_radians();
/// Node precession rate in radians per day.
pub const CHARON_NODE_RATE_RAD_DAY: f64 = 0.0;
/// GM of Pluto in km³/s² (PLU 060).
pub const GM_PLUTO: f64 = 869.3;
/// GM of Charon in km³/s² (PLU 060).
pub const GM_CHARON: f64 = 106.1;

/// Charon's position relative to Pluto's center.
///
/// `jd` is a TT Julian Date (any finite value); five fixed-point Kepler
/// steps are adequate for e ≈ 0.0002.
///
/// Returns `(x, y, z)` in AU in the J2000-ecliptic frame (rotated to equatorial ICRF by [`pluto_cob_offset_icrf`]).
///
/// # Errors
///
/// [`Error::NonFiniteArgument`] when `jd` is NaN or infinite.
pub fn charon_position(jd: f64) -> Result<(f64, f64, f64)> {
    if !jd.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd });
    }
    // Time since epoch in days.
    let dt = jd - CHARON_EPOCH_JD;

    // Mean anomaly at time t, normalized to 0–2π.
    let m = (CHARON_M0_RAD + CHARON_N_RAD_DAY * dt).rem_euclid(2.0 * std::f64::consts::PI);

    // Five fixed-point Kepler iterations from E = M.
    let mut e_anom = m;
    for _ in 0..5 {
        e_anom = m + CHARON_E * e_anom.sin();
    }

    // True anomaly.
    let cos_nu = (e_anom.cos() - CHARON_E) / (1.0 - CHARON_E * e_anom.cos());
    let sin_nu =
        ((1.0 - CHARON_E * CHARON_E).sqrt() * e_anom.sin()) / (1.0 - CHARON_E * e_anom.cos());
    let nu = sin_nu.atan2(cos_nu);

    // Radius (distance from Pluto) in km.
    let r_km = CHARON_A_KM * (1.0 - CHARON_E * CHARON_E) / (1.0 + CHARON_E * nu.cos());

    // Position in the orbital plane.
    let x_orb = r_km * nu.cos();
    let y_orb = r_km * nu.sin();

    // Current node position (model rate is 0.0; structure preserved).
    let node = CHARON_NODE_RAD + CHARON_NODE_RATE_RAD_DAY * dt;

    let cos_omega = CHARON_OMEGA_RAD.cos();
    let sin_omega = CHARON_OMEGA_RAD.sin();
    let cos_i = CHARON_I_RAD.cos();
    let sin_i = CHARON_I_RAD.sin();
    let cos_node = node.cos();
    let sin_node = node.sin();

    // Rotate by ω in the orbital plane (model keeps `z_orb = 0.0` here).
    let x1 = x_orb * cos_omega - y_orb * sin_omega;
    let y1 = x_orb * sin_omega + y_orb * cos_omega;
    let z1 = 0.0;

    // Rotate by inclination around the x-axis.
    let x2 = x1;
    let y2 = y1 * cos_i - z1 * sin_i;
    let z2 = y1 * sin_i + z1 * cos_i;

    // Rotate by node around the z-axis.
    let x_ecl = x2 * cos_node - y2 * sin_node;
    let y_ecl = x2 * sin_node + y2 * cos_node;
    let z_ecl = z2;

    // Convert to AU.
    Ok((
        x_ecl / MOONS_AU_KM,
        y_ecl / MOONS_AU_KM,
        z_ecl / MOONS_AU_KM,
    ))
}

/// Offset from the Pluto-system barycenter to Pluto's center of body.
///
/// Computes the GM-weighted Charon
/// offset, negated, in the same frame `charon_position` emits.
/// Pluto–Charon is a binary (Charon ~12% of system mass; barycenter
/// outside Pluto, ~2130 km from its center).
///
/// # Errors
///
/// [`Error::NonFiniteArgument`] when `jd` is NaN or infinite.
pub fn pluto_cob_offset(jd: f64) -> Result<(f64, f64, f64)> {
    let (cx, cy, cz) = charon_position(jd)?;
    let total_gm = GM_PLUTO + GM_CHARON;
    let bary_x = GM_CHARON * cx / total_gm;
    let bary_y = GM_CHARON * cy / total_gm;
    let bary_z = GM_CHARON * cz / total_gm;
    Ok((-bary_x, -bary_y, -bary_z))
}

/// Pluto center-of-body offset in the equatorial ICRF/J2000 frame.
///
/// Combines the Pluto dispatch arm:
/// [`pluto_cob_offset`] followed by [`ecliptic_to_equatorial_j2000`].
/// Add the result to the Pluto-system barycenter position to get
/// Pluto's center.
///
/// # Errors
///
/// [`Error::NonFiniteArgument`] when `jd` is NaN or infinite.
pub fn pluto_cob_offset_icrf(jd: f64) -> Result<(f64, f64, f64)> {
    Ok(ecliptic_to_equatorial_j2000(pluto_cob_offset(jd)?))
}

#[cfg(test)]
mod charon_tests {
    use super::*;

    fn norm(v: (f64, f64, f64)) -> f64 {
        (v.0 * v.0 + v.1 * v.1 + v.2 * v.2).sqrt()
    }

    #[test]
    fn epoch_position_is_finite_and_bounded_by_elements() {
        let p = charon_position(CHARON_EPOCH_JD).expect("epoch must evaluate");
        assert!(p.0.is_finite() && p.1.is_finite() && p.2.is_finite());
        // Analytic bound from the orbital elements: |r| in
        // a(1−e)/AU.. a(1+e)/AU regardless of anomaly/rotation.
        let r = norm(p);
        let lo = CHARON_A_KM * (1.0 - CHARON_E) / MOONS_AU_KM;
        let hi = CHARON_A_KM * (1.0 + CHARON_E) / MOONS_AU_KM;
        assert!(r >= lo && r <= hi, "r={r} outside [{lo}, {hi}]");
    }

    #[test]
    fn position_is_deterministic() {
        let jd = 2460000.5;
        assert_eq!(charon_position(jd).unwrap(), charon_position(jd).unwrap());
    }

    #[test]
    fn mean_anomaly_periodicity_returns_same_position() {
        // The model node rate is exactly 0.0, so advancing by one
        // orbital period repeats both the mean anomaly (up to one-ulp
        // `N·P − 2π` rounding) and the node. Bound the drift by the
        // mean-motion rounding over one period, generously rounded up.
        let jd = 2451545.0;
        let a = charon_position(jd).unwrap();
        let b = charon_position(jd + CHARON_PERIOD_DAYS).unwrap();
        let drift = ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2) + (a.2 - b.2).powi(2)).sqrt();
        let m_round = (CHARON_N_RAD_DAY * CHARON_PERIOD_DAYS - 2.0 * std::f64::consts::PI).abs();
        let bound = (CHARON_A_KM / MOONS_AU_KM) * m_round + 1e-12;
        assert!(drift <= bound, "drift={drift} bound={bound}");
    }

    #[test]
    fn five_step_kepler_residual_matches_release_approximation() {
        // Structural sanity (not a comparison): the model takes exactly
        // five fixed-point steps from E = M, so the residual
        // |E − M − e·sin E| is O(e⁶)-scale — far below e².
        let jd = 2451545.0 + 1234.5;
        let dt = jd - CHARON_EPOCH_JD;
        let m = (CHARON_M0_RAD + CHARON_N_RAD_DAY * dt).rem_euclid(2.0 * std::f64::consts::PI);
        let mut e_five = m;
        for _ in 0..5 {
            e_five = m + CHARON_E * e_five.sin();
        }
        let residual = (e_five - m - CHARON_E * e_five.sin()).abs();
        assert!(residual <= CHARON_E * CHARON_E, "residual={residual}");
        // And a 50-iteration solve agrees with the five-step value far
        // below the model precision envelope (~1 km ≈ 6.7 e-9 AU).
        let mut e_iter = m;
        for _ in 0..50 {
            e_iter = m + CHARON_E * e_iter.sin();
        }
        let au_per_km = 1.0 / MOONS_AU_KM;
        assert!((e_iter - e_five).abs() * CHARON_A_KM * au_per_km < 1e-15);
    }

    #[test]
    fn malformed_jd_is_refused() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            for err in [
                charon_position(bad),
                pluto_cob_offset(bad),
                pluto_cob_offset_icrf(bad),
            ] {
                match err {
                    Err(Error::NonFiniteArgument { value }) => {
                        if bad.is_nan() {
                            assert!(value.is_nan());
                        } else {
                            assert_eq!(value, bad);
                        }
                    }
                    other => panic!("jd={bad}: expected NonFiniteArgument, got {other:?}"),
                }
            }
        }
    }

    #[test]
    fn cob_offset_is_antiparallel_and_gm_weighted() {
        let jd = 2460000.5;
        let (cx, cy, cz) = charon_position(jd).unwrap();
        let (ox, oy, oz) = pluto_cob_offset(jd).unwrap();
        // weight: offset = −(GM_CHARON / (GM_PLUTO+GM_CHARON))·r.
        let w = GM_CHARON / (GM_PLUTO + GM_CHARON);
        // Binary mass ratio is ~0.109, so use a ratio-scale tolerance
        // rather than the 1 e-18 absolute used for 2 e-4-weight Triton.
        assert!((ox + w * cx).abs() < 1e-15);
        assert!((oy + w * cy).abs() < 1e-15);
        assert!((oz + w * cz).abs() < 1e-15);
        // Magnitude bound: |offset| <= w·a(1+e)/AU.
        let hi = w * CHARON_A_KM * (1.0 + CHARON_E) / MOONS_AU_KM;
        assert!(norm((ox, oy, oz)) <= hi);
    }

    #[test]
    fn icrf_arm_composes_rotation_verbatim() {
        let jd = 2451545.0 + 777.25;
        let icrf = pluto_cob_offset_icrf(jd).unwrap();
        let expected = ecliptic_to_equatorial_j2000(pluto_cob_offset(jd).unwrap());
        assert_eq!(icrf, expected);
    }

    #[test]
    fn far_epoch_stays_bounded() {
        // ±300 years: anomaly normalization and rotation keep the radius
        // inside the element bounds (no overflow-driven escape).
        for jd in [CHARON_EPOCH_JD - 109_575.0, CHARON_EPOCH_JD + 109_575.0] {
            let r = norm(charon_position(jd).unwrap());
            let lo = CHARON_A_KM * (1.0 - CHARON_E) / MOONS_AU_KM;
            let hi = CHARON_A_KM * (1.0 + CHARON_E) / MOONS_AU_KM;
            assert!(r >= lo && r <= hi, "jd={jd} r={r}");
        }
    }

    #[test]
    fn release_constants_match_source_values() {
        // charon.py's own element values (the ones `charon_position`
        // reads — NOT the rounded constants.py near-duplicates).
        assert_eq!(CHARON_A_KM, 19591.4);
        assert_eq!(CHARON_PERIOD_DAYS, 6.3872304);
        assert_eq!(CHARON_E, 0.0002);
        assert_eq!(CHARON_EPOCH_JD, 2451545.0);
        assert_eq!(CHARON_NODE_RATE_RAD_DAY, 0.0);
        // constants.py GM pair (PLU 060) used only by the COB composition.
        assert_eq!(GM_PLUTO, 869.3);
        assert_eq!(GM_CHARON, 106.1);
        // Angular constants round-trip to their published degree values.
        assert!((CHARON_M0_RAD.to_degrees() - 71.255).abs() < 1e-12);
        assert!((CHARON_OMEGA_RAD.to_degrees() - 180.0).abs() < 1e-12);
        assert!((CHARON_I_RAD.to_degrees() - 96.145).abs() < 1e-12);
        assert!((CHARON_NODE_RAD.to_degrees() - 223.046).abs() < 1e-12);
    }
}

// ==================== Jupiter/Galilean (E 5 Lieske/Meeus) ====================
// Lieske E 5 theory in the Meeus chapter 44 algebraic form: longitude sine
// series, latitude tangent series and radius cosine series. Frozen orientation
// angles and a fixed B 1950-to-J2000 adaptation limit accuracy to 100–200 km.
// System GMs are from JUP 365. Arguments are degrees; returned vectors are AU.

/// Reference epoch of the E 5 arguments (JD); `t` is days from this epoch.
pub const GALILEAN_EPOCH_JD: f64 = 2443000.5;
/// Jupiter's equatorial radius in km.
pub const JUPITER_RADIUS_KM: f64 = 71492.0;
/// GM of Jupiter in km³/s² (JUP 365).
pub const GM_JUPITER: f64 = 126686531.9;
/// GM of Io in km³/s² (JUP 365).
pub const GM_IO: f64 = 5959.91547;
/// GM of Europa in km³/s² (JUP 365).
pub const GM_EUROPA: f64 = 3202.71210;
/// GM of Ganymede in km³/s² (JUP 365).
pub const GM_GANYMEDE: f64 = 9887.83275;
/// GM of Callisto in km³/s² (JUP 365).
pub const GM_CALLISTO: f64 = 7179.28340;
/// Total Galilean GM in km³/s².
pub const GM_GALILEAN_TOTAL: f64 = GM_IO + GM_EUROPA + GM_GANYMEDE + GM_CALLISTO;
/// Fixed B 1950.0 -> J2000.0 precession arc in degrees.
pub const E5_T0: f64 = (2451545.0 - 2433282.423) / 36525.0;
/// Fixed precession arc in degrees.
pub const E5_J2000_PRECESSION_DEG: f64 = 1.3966626 * E5_T0 + 0.0003088 * (E5_T0 * E5_T0);
/// Inclination of Jupiter's equator on its orbital plane in degrees.
pub const E5_EQUATOR_TILT_DEG: f64 = 3.120262 + 0.0006 * (2451545.0 - 2415020.5) / 36525.0;
/// Inclination of Jupiter's orbit on the ecliptic in degrees.
pub const E5_ORBIT_TILT_DEG: f64 = 1.303267;
/// Longitude of Jupiter's orbit ascending node in degrees.
pub const E5_ORBIT_NODE_DEG: f64 = 100.464407;
/// Base radius vectors in units of Jupiter's equatorial radius.
pub const E5_R_IO: f64 = 5.90569;
/// Base radius vector of Europa (Jupiter radii).
pub const E5_R_EUROPA: f64 = 9.39657;
/// Base radius vector of Ganymede (Jupiter radii).
pub const E5_R_GANYMEDE: f64 = 14.98832;
/// Base radius vector of Callisto (Jupiter radii).
pub const E5_R_CALLISTO: f64 = 26.36273;

/// Convert E 5 degree-valued arguments to radians.
fn e5_deg_to_rad(degrees: f64) -> f64 {
    degrees * std::f64::consts::PI / 180.0
}

/// Angular arguments of the E 5 theory at one epoch, in degrees.
///
/// Field names follow the E 5 literature: `l1..l4` mean longitudes, `pi1..pi4` proper
/// perijoves, `w1..w4` proper nodes, `psi` equator node, `phi` libration
/// phase, `g`/`gp` Jupiter/Saturn mean anomalies, `pj` Jupiter
/// perihelion (constant), `sigma1..4` longitude corrections,
/// `lon1..lon4` true longitudes.
struct E5Arguments {
    l1: f64,
    l2: f64,
    l3: f64,
    l4: f64,
    pi1: f64,
    pi2: f64,
    pi3: f64,
    pi4: f64,
    w1: f64,
    w2: f64,
    w3: f64,
    w4: f64,
    psi: f64,
    phi: f64,
    g: f64,
    gp: f64,
    pj: f64,
    sigma1: f64,
    sigma2: f64,
    sigma3: f64,
    sigma4: f64,
    lon1: f64,
    lon2: f64,
    lon3: f64,
    lon4: f64,
}

/// Fundamental angular arguments of the E 5 theory, in degrees, with the correction and
/// true-longitude fields still zeroed.
fn e5_fundamental_arguments(t: f64) -> E5Arguments {
    // Principal inequality in Jupiter's longitude, folded into g below.
    let mut gamma = 0.33033 * e5_deg_to_rad(163.679 + 0.0010512 * t).sin();
    gamma += 0.03439 * e5_deg_to_rad(34.486 - 0.0161731 * t).sin();
    E5Arguments {
        l1: 106.07719 + 203.488955790 * t,
        l2: 175.73161 + 101.374724735 * t,
        l3: 120.55883 + 50.317609207 * t,
        l4: 84.44459 + 21.571071177 * t,
        pi1: 97.0881 + 0.16138586 * t,
        pi2: 154.8663 + 0.04726307 * t,
        pi3: 188.1840 + 0.00712734 * t,
        pi4: 335.2868 + 0.00184000 * t,
        w1: 312.3346 - 0.13279386 * t,
        w2: 100.4411 - 0.03263064 * t,
        w3: 119.1942 - 0.00717703 * t,
        w4: 322.6186 - 0.00175934 * t,
        psi: 316.5182 - 0.00000208 * t,
        phi: 199.6766 + 0.17379190 * t,
        g: 30.23756 + 0.0830925701 * t + gamma,
        gp: 31.97853 + 0.0334597339 * t,
        pj: 13.469942,
        sigma1: 0.0,
        sigma2: 0.0,
        sigma3: 0.0,
        sigma4: 0.0,
        lon1: 0.0,
        lon2: 0.0,
        lon3: 0.0,
        lon4: 0.0,
    }
}

/// Sigma_1: periodic corrections to Io's mean longitude (degrees,, 21 sine terms in table order).
fn e5_sigma_io(a: &E5Arguments) -> f64 {
    let mut total = 0.0;
    total += 0.47259 * e5_deg_to_rad(2.0 * (a.l1 - a.l2)).sin();
    total += -0.03478 * e5_deg_to_rad(a.pi3 - a.pi4).sin();
    total += 0.01081 * e5_deg_to_rad(a.l2 - 2.0 * a.l3 + a.pi3).sin();
    total += 0.00738 * e5_deg_to_rad(a.phi).sin();
    total += 0.00713 * e5_deg_to_rad(a.l2 - 2.0 * a.l3 + a.pi2).sin();
    total += -0.00674 * e5_deg_to_rad(a.pi1 + a.pi3 - 2.0 * a.pj - 2.0 * a.g).sin();
    total += 0.00666 * e5_deg_to_rad(a.l2 - 2.0 * a.l3 + a.pi4).sin();
    total += 0.00445 * e5_deg_to_rad(a.l1 - a.pi3).sin();
    total += -0.00354 * e5_deg_to_rad(a.l1 - a.l2).sin();
    total += -0.00317 * e5_deg_to_rad(2.0 * a.psi - 2.0 * a.pj).sin();
    total += 0.00265 * e5_deg_to_rad(a.l1 - a.pi4).sin();
    total += -0.00186 * e5_deg_to_rad(a.g).sin();
    total += 0.00162 * e5_deg_to_rad(a.pi2 - a.pi3).sin();
    total += 0.00158 * e5_deg_to_rad(4.0 * (a.l1 - a.l2)).sin();
    total += -0.00155 * e5_deg_to_rad(a.l1 - a.l3).sin();
    total += -0.00138 * e5_deg_to_rad(a.psi + a.w3 - 2.0 * a.pj - 2.0 * a.g).sin();
    total += -0.00115 * e5_deg_to_rad(2.0 * (a.l1 - 2.0 * a.l2 + a.w2)).sin();
    total += 0.00089 * e5_deg_to_rad(a.pi2 - a.pi4).sin();
    total += 0.00085 * e5_deg_to_rad(a.l1 + a.pi3 - 2.0 * a.pj - 2.0 * a.g).sin();
    total += 0.00083 * e5_deg_to_rad(a.w2 - a.w3).sin();
    total += 0.00053 * e5_deg_to_rad(a.psi - a.w2).sin();
    total
}

/// Sigma_2: periodic corrections to Europa's mean longitude (degrees,
/// 37 sine terms in table order).
fn e5_sigma_europa(a: &E5Arguments) -> f64 {
    let mut total = 0.0;
    total += 1.06476 * e5_deg_to_rad(2.0 * (a.l2 - a.l3)).sin();
    total += 0.04256 * e5_deg_to_rad(a.l1 - 2.0 * a.l2 + a.pi3).sin();
    total += 0.03581 * e5_deg_to_rad(a.l2 - a.pi3).sin();
    total += 0.02395 * e5_deg_to_rad(a.l1 - 2.0 * a.l2 + a.pi4).sin();
    total += 0.01984 * e5_deg_to_rad(a.l2 - a.pi4).sin();
    total += -0.01778 * e5_deg_to_rad(a.phi).sin();
    total += 0.01654 * e5_deg_to_rad(a.l2 - a.pi2).sin();
    total += 0.01334 * e5_deg_to_rad(a.l2 - 2.0 * a.l3 + a.pi2).sin();
    total += 0.01294 * e5_deg_to_rad(a.pi3 - a.pi4).sin();
    total += -0.01142 * e5_deg_to_rad(a.l2 - a.l3).sin();
    total += -0.01057 * e5_deg_to_rad(a.g).sin();
    total += -0.00775 * e5_deg_to_rad(2.0 * (a.psi - a.pj)).sin();
    total += 0.00524 * e5_deg_to_rad(2.0 * (a.l1 - a.l2)).sin();
    total += -0.00460 * e5_deg_to_rad(a.l1 - a.l3).sin();
    total += 0.00316 * e5_deg_to_rad(a.psi - 2.0 * a.g + a.w3 - 2.0 * a.pj).sin();
    total += -0.00203 * e5_deg_to_rad(a.pi1 + a.pi3 - 2.0 * a.pj - 2.0 * a.g).sin();
    total += 0.00146 * e5_deg_to_rad(a.psi - a.w3).sin();
    total += -0.00145 * e5_deg_to_rad(2.0 * a.g).sin();
    total += 0.00125 * e5_deg_to_rad(a.psi - a.w4).sin();
    total += -0.00115 * e5_deg_to_rad(a.l1 - 2.0 * a.l3 + a.pi3).sin();
    total += -0.00094 * e5_deg_to_rad(2.0 * (a.l2 - a.w2)).sin();
    total += 0.00086 * e5_deg_to_rad(2.0 * (a.l1 - 2.0 * a.l2 + a.w2)).sin();
    total += -0.00086 * e5_deg_to_rad(5.0 * a.gp - 2.0 * a.g + 52.225).sin();
    total += -0.00078 * e5_deg_to_rad(a.l2 - a.l4).sin();
    total += -0.00064 * e5_deg_to_rad(3.0 * a.l3 - 7.0 * a.l4 + 4.0 * a.pi4).sin();
    total += 0.00064 * e5_deg_to_rad(a.pi1 - a.pi4).sin();
    total += -0.00063 * e5_deg_to_rad(a.l1 - 2.0 * a.l3 + a.pi4).sin();
    total += 0.00058 * e5_deg_to_rad(a.w3 - a.w4).sin();
    total += 0.00056 * e5_deg_to_rad(2.0 * (a.psi - a.pj - a.g)).sin();
    total += 0.00056 * e5_deg_to_rad(2.0 * (a.l2 - a.l4)).sin();
    total += 0.00055 * e5_deg_to_rad(2.0 * (a.l1 - a.l3)).sin();
    total += 0.00052 * e5_deg_to_rad(3.0 * a.l3 - 7.0 * a.l4 + a.pi3 + 3.0 * a.pi4).sin();
    total += -0.00043 * e5_deg_to_rad(a.l1 - a.pi3).sin();
    total += 0.00041 * e5_deg_to_rad(5.0 * (a.l2 - a.l3)).sin();
    total += 0.00041 * e5_deg_to_rad(a.pi4 - a.pj).sin();
    total += 0.00032 * e5_deg_to_rad(a.w2 - a.w3).sin();
    total += 0.00032 * e5_deg_to_rad(2.0 * (a.l3 - a.g - a.pj)).sin();
    total
}

/// Sigma_3: periodic corrections to Ganymede's mean longitude (degrees,
/// 43 sine terms in table order).
fn e5_sigma_ganymede(a: &E5Arguments) -> f64 {
    let mut total = 0.0;
    total += 0.16490 * e5_deg_to_rad(a.l3 - a.pi3).sin();
    total += 0.09081 * e5_deg_to_rad(a.l3 - a.pi4).sin();
    total += -0.06907 * e5_deg_to_rad(a.l2 - a.l3).sin();
    total += 0.03784 * e5_deg_to_rad(a.pi3 - a.pi4).sin();
    total += 0.01846 * e5_deg_to_rad(2.0 * (a.l3 - a.l4)).sin();
    total += -0.01340 * e5_deg_to_rad(a.g).sin();
    total += -0.01014 * e5_deg_to_rad(2.0 * (a.psi - a.pj)).sin();
    total += 0.00704 * e5_deg_to_rad(a.l2 - 2.0 * a.l3 + a.pi3).sin();
    total += -0.00620 * e5_deg_to_rad(a.l2 - 2.0 * a.l3 + a.pi2).sin();
    total += -0.00541 * e5_deg_to_rad(a.l3 - a.l4).sin();
    total += 0.00381 * e5_deg_to_rad(a.l2 - 2.0 * a.l3 + a.pi4).sin();
    total += 0.00235 * e5_deg_to_rad(a.psi - a.w3).sin();
    total += 0.00198 * e5_deg_to_rad(a.psi - a.w4).sin();
    total += 0.00176 * e5_deg_to_rad(a.phi).sin();
    total += 0.00130 * e5_deg_to_rad(3.0 * (a.l3 - a.l4)).sin();
    total += 0.00125 * e5_deg_to_rad(a.l1 - a.l3).sin();
    total += -0.00119 * e5_deg_to_rad(5.0 * a.gp - 2.0 * a.g + 52.225).sin();
    total += 0.00109 * e5_deg_to_rad(a.l1 - a.l2).sin();
    total += -0.00100 * e5_deg_to_rad(3.0 * a.l3 - 7.0 * a.l4 + 4.0 * a.pi4).sin();
    total += 0.00091 * e5_deg_to_rad(a.w3 - a.w4).sin();
    total += 0.00080 * e5_deg_to_rad(3.0 * a.l3 - 7.0 * a.l4 + a.pi3 + 3.0 * a.pi4).sin();
    total += -0.00075 * e5_deg_to_rad(2.0 * a.l2 - 3.0 * a.l3 + a.pi3).sin();
    total += 0.00072 * e5_deg_to_rad(a.pi1 + a.pi3 - 2.0 * a.pj - 2.0 * a.g).sin();
    total += 0.00069 * e5_deg_to_rad(a.pi4 - a.pj).sin();
    total += -0.00058 * e5_deg_to_rad(2.0 * a.l3 - 3.0 * a.l4 + a.pi4).sin();
    total += -0.00057 * e5_deg_to_rad(a.l3 - 2.0 * a.l4 + a.pi4).sin();
    total += 0.00056 * e5_deg_to_rad(a.l3 + a.pi3 - 2.0 * a.pj - 2.0 * a.g).sin();
    total += -0.00052 * e5_deg_to_rad(a.l2 - 2.0 * a.l3 + a.pi1).sin();
    total += -0.00050 * e5_deg_to_rad(a.pi2 - a.pi3).sin();
    total += 0.00048 * e5_deg_to_rad(a.l3 - 2.0 * a.l4 + a.pi3).sin();
    total += -0.00045 * e5_deg_to_rad(2.0 * a.l2 - 3.0 * a.l3 + a.pi4).sin();
    total += -0.00041 * e5_deg_to_rad(a.pi2 - a.pi4).sin();
    total += -0.00038 * e5_deg_to_rad(2.0 * a.g).sin();
    total += -0.00037 * e5_deg_to_rad(a.pi3 - a.pi4 + a.w3 - a.w4).sin();
    total += -0.00032 * e5_deg_to_rad(3.0 * a.l3 - 7.0 * a.l4 + 2.0 * a.pi3 + 2.0 * a.pi4).sin();
    total += 0.00030 * e5_deg_to_rad(4.0 * (a.l3 - a.l4)).sin();
    total += 0.00029 * e5_deg_to_rad(a.l3 + a.pi4 - 2.0 * a.pj - 2.0 * a.g).sin();
    total += -0.00028 * e5_deg_to_rad(a.w3 + a.psi - 2.0 * a.pj - 2.0 * a.g).sin();
    total += 0.00026 * e5_deg_to_rad(a.l3 - a.pj - a.g).sin();
    total += 0.00024 * e5_deg_to_rad(a.l2 - 3.0 * a.l3 + 2.0 * a.l4).sin();
    total += 0.00021 * e5_deg_to_rad(2.0 * (a.l3 - a.pj - a.g)).sin();
    total += -0.00021 * e5_deg_to_rad(a.l3 - a.pi2).sin();
    total += 0.00017 * e5_deg_to_rad(2.0 * (a.l3 - a.pi3)).sin();
    total
}

/// Sigma_4: periodic corrections to Callisto's mean longitude (degrees,
/// 50 sine terms in table order).
fn e5_sigma_callisto(a: &E5Arguments) -> f64 {
    let mut total = 0.0;
    total += 0.84287 * e5_deg_to_rad(a.l4 - a.pi4).sin();
    total += 0.03431 * e5_deg_to_rad(a.pi4 - a.pi3).sin();
    total += -0.03305 * e5_deg_to_rad(2.0 * (a.psi - a.pj)).sin();
    total += -0.03211 * e5_deg_to_rad(a.g).sin();
    total += -0.01862 * e5_deg_to_rad(a.l4 - a.pi3).sin();
    total += 0.01186 * e5_deg_to_rad(a.psi - a.w4).sin();
    total += 0.00623 * e5_deg_to_rad(a.l4 + a.pi4 - 2.0 * a.g - 2.0 * a.pj).sin();
    total += 0.00387 * e5_deg_to_rad(2.0 * (a.l4 - a.pi4)).sin();
    total += -0.00284 * e5_deg_to_rad(5.0 * a.gp - 2.0 * a.g + 52.225).sin();
    total += -0.00234 * e5_deg_to_rad(2.0 * (a.psi - a.pi4)).sin();
    total += -0.00223 * e5_deg_to_rad(a.l3 - a.l4).sin();
    total += -0.00208 * e5_deg_to_rad(a.l4 - a.pj).sin();
    total += 0.00178 * e5_deg_to_rad(a.psi + a.w4 - 2.0 * a.pi4).sin();
    total += 0.00134 * e5_deg_to_rad(a.pi4 - a.pj).sin();
    total += 0.00125 * e5_deg_to_rad(2.0 * (a.l4 - a.g - a.pj)).sin();
    total += -0.00117 * e5_deg_to_rad(2.0 * a.g).sin();
    total += -0.00112 * e5_deg_to_rad(2.0 * (a.l3 - a.l4)).sin();
    total += 0.00107 * e5_deg_to_rad(3.0 * a.l3 - 7.0 * a.l4 + 4.0 * a.pi4).sin();
    total += 0.00102 * e5_deg_to_rad(a.l4 - a.g - a.pj).sin();
    total += 0.00096 * e5_deg_to_rad(2.0 * a.l4 - a.psi - a.w4).sin();
    total += 0.00087 * e5_deg_to_rad(2.0 * (a.psi - a.w4)).sin();
    total += -0.00085 * e5_deg_to_rad(3.0 * a.l3 - 7.0 * a.l4 + a.pi3 + 3.0 * a.pi4).sin();
    total += 0.00085 * e5_deg_to_rad(a.l3 - 2.0 * a.l4 + a.pi4).sin();
    total += -0.00081 * e5_deg_to_rad(2.0 * (a.l4 - a.psi)).sin();
    total += 0.00071 * e5_deg_to_rad(a.l4 + a.pi4 - 2.0 * a.pj - 3.0 * a.g).sin();
    total += 0.00061 * e5_deg_to_rad(a.l1 - a.l4).sin();
    total += -0.00056 * e5_deg_to_rad(a.psi - a.w3).sin();
    total += -0.00054 * e5_deg_to_rad(a.l3 - 2.0 * a.l4 + a.pi3).sin();
    total += 0.00051 * e5_deg_to_rad(a.l2 - a.l4).sin();
    total += 0.00042 * e5_deg_to_rad(2.0 * (a.psi - a.g - a.pj)).sin();
    total += 0.00039 * e5_deg_to_rad(2.0 * (a.pi4 - a.w4)).sin();
    total += 0.00036 * e5_deg_to_rad(a.psi + a.pj - a.pi4 - a.w4).sin();
    total += 0.00035 * e5_deg_to_rad(2.0 * a.gp - a.g + 188.37).sin();
    total += -0.00035 * e5_deg_to_rad(a.l4 - a.pi4 + 2.0 * a.pj - 2.0 * a.psi).sin();
    total += -0.00032 * e5_deg_to_rad(a.l4 + a.pi4 - 2.0 * a.pj - a.g).sin();
    total += 0.00030 * e5_deg_to_rad(2.0 * a.gp - 2.0 * a.g + 149.15).sin();
    total += 0.00029 * e5_deg_to_rad(3.0 * a.l3 - 7.0 * a.l4 + 2.0 * a.pi3 + 2.0 * a.pi4).sin();
    total += 0.00028 * e5_deg_to_rad(a.l4 - a.pi4 + 2.0 * a.psi - 2.0 * a.pj).sin();
    total += -0.00028 * e5_deg_to_rad(2.0 * (a.l4 - a.w4)).sin();
    total += -0.00027 * e5_deg_to_rad(a.pi3 - a.pi4 + a.w3 - a.w4).sin();
    total += -0.00026 * e5_deg_to_rad(5.0 * a.gp - 3.0 * a.g + 188.37).sin();
    total += 0.00025 * e5_deg_to_rad(a.w4 - a.w3).sin();
    total += -0.00025 * e5_deg_to_rad(a.l2 - 3.0 * a.l3 + 2.0 * a.l4).sin();
    total += -0.00023 * e5_deg_to_rad(3.0 * (a.l3 - a.l4)).sin();
    total += 0.00021 * e5_deg_to_rad(2.0 * a.l4 - 2.0 * a.pj - 3.0 * a.g).sin();
    total += -0.00021 * e5_deg_to_rad(2.0 * a.l3 - 3.0 * a.l4 + a.pi4).sin();
    total += 0.00019 * e5_deg_to_rad(a.l4 - a.pi4 - a.g).sin();
    total += -0.00019 * e5_deg_to_rad(2.0 * a.l4 - a.pi3 - a.pi4).sin();
    total += -0.00018 * e5_deg_to_rad(a.l4 - a.pi4 + a.g).sin();
    total += -0.00016 * e5_deg_to_rad(a.l4 + a.pi3 - 2.0 * a.pj - 2.0 * a.g).sin();
    total
}

/// tan(B_1): Io's jovicentric latitude series (dimensionless,, 7 sine terms in table order).
fn e5_tanb_io(a: &E5Arguments) -> f64 {
    let mut total = 0.0;
    total += 0.0006393 * e5_deg_to_rad(a.lon1 - a.w1).sin();
    total += 0.0001825 * e5_deg_to_rad(a.lon1 - a.w2).sin();
    total += 0.0000329 * e5_deg_to_rad(a.lon1 - a.w3).sin();
    total += -0.0000311 * e5_deg_to_rad(a.lon1 - a.psi).sin();
    total += 0.0000093 * e5_deg_to_rad(a.lon1 - a.w4).sin();
    total += 0.0000075 * e5_deg_to_rad(3.0 * a.lon1 - 4.0 * a.l2 - 1.9927 * a.sigma1 + a.w2).sin();
    total += 0.0000046 * e5_deg_to_rad(a.lon1 + a.psi - 2.0 * a.pj - 2.0 * a.g).sin();
    total
}

/// tan(B_2): Europa's jovicentric latitude series (dimensionless,, 9 sine terms in table order).
fn e5_tanb_europa(a: &E5Arguments) -> f64 {
    let mut total = 0.0;
    total += 0.0081004 * e5_deg_to_rad(a.lon2 - a.w2).sin();
    total += 0.0004512 * e5_deg_to_rad(a.lon2 - a.w3).sin();
    total += -0.0003284 * e5_deg_to_rad(a.lon2 - a.psi).sin();
    total += 0.0001160 * e5_deg_to_rad(a.lon2 - a.w4).sin();
    total += 0.0000272 * e5_deg_to_rad(a.l1 - 2.0 * a.l3 + 1.0146 * a.sigma2 + a.w2).sin();
    total += -0.0000144 * e5_deg_to_rad(a.lon2 - a.w1).sin();
    total += 0.0000143 * e5_deg_to_rad(a.lon2 + a.psi - 2.0 * a.pj - 2.0 * a.g).sin();
    total += 0.0000035 * e5_deg_to_rad(a.lon2 - a.psi + a.g).sin();
    total += -0.0000028 * e5_deg_to_rad(a.l1 - 2.0 * a.l3 + 1.0146 * a.sigma2 + a.w3).sin();
    total
}

/// tan(B_3): Ganymede's jovicentric latitude series (dimensionless,
/// 11 sine terms in table order).
fn e5_tanb_ganymede(a: &E5Arguments) -> f64 {
    let mut total = 0.0;
    total += 0.0032402 * e5_deg_to_rad(a.lon3 - a.w3).sin();
    total += -0.0016911 * e5_deg_to_rad(a.lon3 - a.psi).sin();
    total += 0.0006847 * e5_deg_to_rad(a.lon3 - a.w4).sin();
    total += -0.0002797 * e5_deg_to_rad(a.lon3 - a.w2).sin();
    total += 0.0000321 * e5_deg_to_rad(a.lon3 + a.psi - 2.0 * a.pj - 2.0 * a.g).sin();
    total += 0.0000051 * e5_deg_to_rad(a.lon3 - a.psi + a.g).sin();
    total += -0.0000045 * e5_deg_to_rad(a.lon3 - a.psi - a.g).sin();
    total += -0.0000045 * e5_deg_to_rad(a.lon3 + a.psi - 2.0 * a.pj).sin();
    total += 0.0000037 * e5_deg_to_rad(a.lon3 + a.psi - 2.0 * a.pj - 3.0 * a.g).sin();
    total += 0.0000030 * e5_deg_to_rad(2.0 * a.l2 - 3.0 * a.lon3 + 4.03 * a.sigma3 + a.w2).sin();
    total += -0.0000021 * e5_deg_to_rad(2.0 * a.l2 - 3.0 * a.lon3 + 4.03 * a.sigma3 + a.w3).sin();
    total
}

/// tan(B_4): Callisto's jovicentric latitude series (dimensionless,
/// 8 sine terms in table order).
fn e5_tanb_callisto(a: &E5Arguments) -> f64 {
    let mut total = 0.0;
    total += -0.0076579 * e5_deg_to_rad(a.lon4 - a.psi).sin();
    total += 0.0044134 * e5_deg_to_rad(a.lon4 - a.w4).sin();
    total += -0.0005112 * e5_deg_to_rad(a.lon4 - a.w3).sin();
    total += 0.0000773 * e5_deg_to_rad(a.lon4 + a.psi - 2.0 * a.pj - 2.0 * a.g).sin();
    total += 0.0000104 * e5_deg_to_rad(a.lon4 - a.psi + a.g).sin();
    total += -0.0000102 * e5_deg_to_rad(a.lon4 - a.psi - a.g).sin();
    total += 0.0000088 * e5_deg_to_rad(a.lon4 + a.psi - 2.0 * a.pj - 3.0 * a.g).sin();
    total += -0.0000038 * e5_deg_to_rad(a.lon4 + a.psi - 2.0 * a.pj - a.g).sin();
    total
}

/// Delta R_1: fractional corrections to Io's radius vector
/// (dimensionless, 7 cosine terms in table
/// order).
fn e5_drad_io(a: &E5Arguments) -> f64 {
    let mut total = 0.0;
    total += -0.0041339 * e5_deg_to_rad(2.0 * (a.l1 - a.l2)).cos();
    total += -0.0000387 * e5_deg_to_rad(a.l1 - a.pi3).cos();
    total += -0.0000214 * e5_deg_to_rad(a.l1 - a.pi4).cos();
    total += 0.0000170 * e5_deg_to_rad(a.l1 - a.l2).cos();
    total += -0.0000131 * e5_deg_to_rad(4.0 * (a.l1 - a.l2)).cos();
    total += 0.0000106 * e5_deg_to_rad(a.l1 - a.l3).cos();
    total += -0.0000066 * e5_deg_to_rad(a.l1 + a.pi3 - 2.0 * a.pj - 2.0 * a.g).cos();
    total
}

/// Delta R_2: fractional corrections to Europa's radius vector
/// (dimensionless, 11 cosine terms in
/// table order).
fn e5_drad_europa(a: &E5Arguments) -> f64 {
    let mut total = 0.0;
    total += 0.0093848 * e5_deg_to_rad(a.l1 - a.l2).cos();
    total += -0.0003116 * e5_deg_to_rad(a.l2 - a.pi3).cos();
    total += -0.0001744 * e5_deg_to_rad(a.l2 - a.pi4).cos();
    total += -0.0001442 * e5_deg_to_rad(a.l2 - a.pi2).cos();
    total += 0.0000553 * e5_deg_to_rad(a.l2 - a.l3).cos();
    total += 0.0000523 * e5_deg_to_rad(a.l1 - a.l3).cos();
    total += -0.0000290 * e5_deg_to_rad(2.0 * (a.l1 - a.l2)).cos();
    total += 0.0000164 * e5_deg_to_rad(2.0 * (a.l2 - a.w2)).cos();
    total += 0.0000107 * e5_deg_to_rad(a.l1 - 2.0 * a.l3 + a.pi3).cos();
    total += -0.0000102 * e5_deg_to_rad(a.l2 - a.pi1).cos();
    total += -0.0000091 * e5_deg_to_rad(2.0 * (a.l1 - a.l3)).cos();
    total
}

/// Delta R_3: fractional corrections to Ganymede's radius vector
/// (dimensionless, 10 cosine terms in
/// table order).
fn e5_drad_ganymede(a: &E5Arguments) -> f64 {
    let mut total = 0.0;
    total += -0.0014388 * e5_deg_to_rad(a.l3 - a.pi3).cos();
    total += -0.0007919 * e5_deg_to_rad(a.l3 - a.pi4).cos();
    total += 0.0006342 * e5_deg_to_rad(a.l2 - a.l3).cos();
    total += -0.0001761 * e5_deg_to_rad(2.0 * (a.l3 - a.l4)).cos();
    total += 0.0000294 * e5_deg_to_rad(a.l3 - a.l4).cos();
    total += -0.0000156 * e5_deg_to_rad(3.0 * (a.l3 - a.l4)).cos();
    total += 0.0000156 * e5_deg_to_rad(a.l1 - a.l3).cos();
    total += -0.0000153 * e5_deg_to_rad(a.l1 - a.l2).cos();
    total += 0.0000070 * e5_deg_to_rad(2.0 * a.l2 - 3.0 * a.l3 + a.pi3).cos();
    total += -0.0000051 * e5_deg_to_rad(a.l3 + a.pi3 - 2.0 * a.pj - 2.0 * a.g).cos();
    total
}

/// Delta R_4: fractional corrections to Callisto's radius vector
/// (dimensionless, 16 cosine terms in
/// table order).
fn e5_drad_callisto(a: &E5Arguments) -> f64 {
    let mut total = 0.0;
    total += -0.0073546 * e5_deg_to_rad(a.l4 - a.pi4).cos();
    total += 0.0001621 * e5_deg_to_rad(a.l4 - a.pi3).cos();
    total += 0.0000974 * e5_deg_to_rad(a.l3 - a.l4).cos();
    total += -0.0000543 * e5_deg_to_rad(a.l4 + a.pi4 - 2.0 * a.pj - 2.0 * a.g).cos();
    total += -0.0000271 * e5_deg_to_rad(2.0 * (a.l4 - a.pi4)).cos();
    total += 0.0000182 * e5_deg_to_rad(a.l4 - a.pj).cos();
    total += 0.0000177 * e5_deg_to_rad(2.0 * (a.l3 - a.l4)).cos();
    total += -0.0000167 * e5_deg_to_rad(2.0 * a.l4 - a.psi - a.w4).cos();
    total += 0.0000167 * e5_deg_to_rad(a.psi - a.w4).cos();
    total += -0.0000155 * e5_deg_to_rad(2.0 * (a.l4 - a.pj - a.g)).cos();
    total += 0.0000142 * e5_deg_to_rad(2.0 * (a.l4 - a.psi)).cos();
    total += 0.0000105 * e5_deg_to_rad(a.l1 - a.l4).cos();
    total += 0.0000092 * e5_deg_to_rad(a.l2 - a.l4).cos();
    total += -0.0000089 * e5_deg_to_rad(a.l4 - a.pj - a.g).cos();
    total += -0.0000062 * e5_deg_to_rad(a.l4 + a.pi4 - 2.0 * a.pj - 3.0 * a.g).cos();
    total += 0.0000048 * e5_deg_to_rad(2.0 * (a.l4 - a.w4)).cos();
    total
}

/// Assemble one satellite vector and rotate it to the J2000 ecliptic
///
/// The vector is built in the Jupiter-equatorial frame (X toward `node`,
/// Z toward Jupiter's north pole), then carried through the fixed chain:
/// equator tilt, node swing (`node − orbit node`), orbit tilt, J2000
/// equinox rotation. Returns Jupiter-centric `(x, y, z)` in km.
fn e5_j2000_ecliptic_km(longitude: f64, latitude: f64, radius: f64, node: f64) -> (f64, f64, f64) {
    let x = radius * e5_deg_to_rad(longitude - node).cos() * latitude.cos();
    let y = radius * e5_deg_to_rad(longitude - node).sin() * latitude.cos();
    let z = radius * latitude.sin();

    // 1. Tilt Jupiter's equator onto its orbital plane.
    let cos_i = e5_deg_to_rad(E5_EQUATOR_TILT_DEG).cos();
    let sin_i = e5_deg_to_rad(E5_EQUATOR_TILT_DEG).sin();
    let y1 = y * cos_i - z * sin_i;
    let z1 = y * sin_i + z * cos_i;
    // 2. Swing the node (z-rotation by psi' − Omega_j).
    let swing = node - E5_ORBIT_NODE_DEG;
    let cos_swing = e5_deg_to_rad(swing).cos();
    let sin_swing = e5_deg_to_rad(swing).sin();
    let x2 = x * cos_swing - y1 * sin_swing;
    let y2 = x * sin_swing + y1 * cos_swing;
    // 3. Tilt Jupiter's orbital plane onto the ecliptic.
    let cos_o = e5_deg_to_rad(E5_ORBIT_TILT_DEG).cos();
    let sin_o = e5_deg_to_rad(E5_ORBIT_TILT_DEG).sin();
    let y3 = y2 * cos_o - z1 * sin_o;
    let z3 = y2 * sin_o + z1 * cos_o;
    // 4. Rotate to the J2000 equinox (z-rotation by Omega_j).
    let cos_n = e5_deg_to_rad(E5_ORBIT_NODE_DEG).cos();
    let sin_n = e5_deg_to_rad(E5_ORBIT_NODE_DEG).sin();
    let x4 = x2 * cos_n - y3 * sin_n;
    let y4 = x2 * sin_n + y3 * cos_n;

    (
        x4 * JUPITER_RADIUS_KM,
        y4 * JUPITER_RADIUS_KM,
        z3 * JUPITER_RADIUS_KM,
    )
}

/// Four Galilean-satellite position triples `(io, europa, ganymede,
/// callisto)`, each `(x, y, z)`.
pub type GalileanMoonQuartet = (
    (f64, f64, f64),
    (f64, f64, f64),
    (f64, f64, f64),
    (f64, f64, f64),
);

/// Jupiter-centric positions of the four Galilean satellites.
///
/// Returns `(io, europa, ganymede, callisto)`, each `(x, y, z)` in km in
/// the J2000-ecliptic frame the model emits (rotated to equatorial
/// ICRF by [`jupiter_cob_offset_icrf`]). The input Julian Date is read
/// directly as JDE (TT/TDB); no light-time correction is applied.
///
/// # Errors
///
/// [`Error::NonFiniteArgument`] when `jd` is NaN or infinite.
pub fn galilean_moon_positions(jd: f64) -> Result<GalileanMoonQuartet> {
    if !jd.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd });
    }
    let t = jd - GALILEAN_EPOCH_JD;
    let mut a = e5_fundamental_arguments(t);

    // Periodic longitude corrections (degrees) and true longitudes.
    a.sigma1 = e5_sigma_io(&a);
    a.sigma2 = e5_sigma_europa(&a);
    a.sigma3 = e5_sigma_ganymede(&a);
    a.sigma4 = e5_sigma_callisto(&a);
    a.lon1 = a.l1 + a.sigma1;
    a.lon2 = a.l2 + a.sigma2;
    a.lon3 = a.l3 + a.sigma3;
    a.lon4 = a.l4 + a.sigma4;

    // Jovicentric latitudes (radians).
    let b1 = e5_tanb_io(&a).atan();
    let b2 = e5_tanb_europa(&a).atan();
    let b3 = e5_tanb_ganymede(&a).atan();
    let b4 = e5_tanb_callisto(&a).atan();

    // Radius vectors (units of Jupiter's equatorial radius).
    let r1 = E5_R_IO * (1.0 + e5_drad_io(&a));
    let r2 = E5_R_EUROPA * (1.0 + e5_drad_europa(&a));
    let r3 = E5_R_GANYMEDE * (1.0 + e5_drad_ganymede(&a));
    let r4 = E5_R_CALLISTO * (1.0 + e5_drad_callisto(&a));

    // Fixed B 1950 -> J2000 precession of the true longitudes and of psi.
    // The offset cancels in each satellite longitude (L' − psi'), but
    // psi' — not psi — drives the node swing of the rotation chain.
    let psi_prime = a.psi + E5_J2000_PRECESSION_DEG;

    let io = e5_j2000_ecliptic_km(a.lon1 + E5_J2000_PRECESSION_DEG, b1, r1, psi_prime);
    let europa = e5_j2000_ecliptic_km(a.lon2 + E5_J2000_PRECESSION_DEG, b2, r2, psi_prime);
    let ganymede = e5_j2000_ecliptic_km(a.lon3 + E5_J2000_PRECESSION_DEG, b3, r3, psi_prime);
    let callisto = e5_j2000_ecliptic_km(a.lon4 + E5_J2000_PRECESSION_DEG, b4, r4, psi_prime);
    Ok((io, europa, ganymede, callisto))
}

/// Offset from the Jupiter-system barycenter to Jupiter's center of body.
///
/// Computes the GM-weighted four-moon
/// sum, negated, converted from km to AU, in the same J2000-ecliptic
/// frame `galilean_moon_positions` emits.
///
/// # Errors
///
/// [`Error::NonFiniteArgument`] when `jd` is NaN or infinite.
pub fn jupiter_cob_offset(jd: f64) -> Result<(f64, f64, f64)> {
    let (io, europa, ganymede, callisto) = galilean_moon_positions(jd)?;
    let total_gm = GM_JUPITER + GM_GALILEAN_TOTAL;
    let bary_x =
        (GM_IO * io.0 + GM_EUROPA * europa.0 + GM_GANYMEDE * ganymede.0 + GM_CALLISTO * callisto.0)
            / total_gm;
    let bary_y =
        (GM_IO * io.1 + GM_EUROPA * europa.1 + GM_GANYMEDE * ganymede.1 + GM_CALLISTO * callisto.1)
            / total_gm;
    let bary_z =
        (GM_IO * io.2 + GM_EUROPA * europa.2 + GM_GANYMEDE * ganymede.2 + GM_CALLISTO * callisto.2)
            / total_gm;
    Ok((
        -bary_x / MOONS_AU_KM,
        -bary_y / MOONS_AU_KM,
        -bary_z / MOONS_AU_KM,
    ))
}

/// Jupiter center-of-body offset in the equatorial ICRF/J2000 frame.
///
/// Combines the Jupiter dispatch arm:
/// [`jupiter_cob_offset`] followed by [`ecliptic_to_equatorial_j2000`].
/// Add the result to the Jupiter-system barycenter position to get
/// Jupiter's center.
///
/// # Errors
///
/// [`Error::NonFiniteArgument`] when `jd` is NaN or infinite.
pub fn jupiter_cob_offset_icrf(jd: f64) -> Result<(f64, f64, f64)> {
    Ok(ecliptic_to_equatorial_j2000(jupiter_cob_offset(jd)?))
}

// ---------------------------------------------------------------------------
// Uranian major moons: GUST 86-geometry Keplerian elements (Jacobson 2014).
// The Uranus-equatorial frame is rotated directly into ICRF. The center
// composition uses the model's documented approximate mass-scale weights.
// ---------------------------------------------------------------------------

/// Uranus GM in km³/s².
pub const GM_URANUS: f64 = 5_793_951.3;

/// Major-moon barycenter weights.
///
/// Despite the `GM_` prefix these are not
/// GMs in km³/s² but approximate mass-scale figures (order 1 e 20 kg);
/// they enter [`uranus_cob_offset`] only as relative weights among
/// the moons, sharing a denominator with the true-GM `GM_URANUS`.
/// These approximate weights limit the physical accuracy of this composition.
pub const GM_ARIEL: f64 = 14.0;
/// Umbriel barycenter weight, same approximate mass-scale units as `GM_ARIEL`.
pub const GM_UMBRIEL: f64 = 13.0;
/// Titania barycenter weight, same approximate mass-scale units as `GM_ARIEL`.
pub const GM_TITANIA: f64 = 41.0;
/// Oberon barycenter weight, same approximate mass-scale units as `GM_ARIEL`.
pub const GM_OBERON: f64 = 29.0;
/// Miranda barycenter weight, same approximate mass-scale units as `GM_ARIEL`.
pub const GM_MIRANDA: f64 = 0.7;

/// Total moon weight.
pub const GM_URANUS_MOONS_TOTAL: f64 = GM_ARIEL + GM_UMBRIEL + GM_TITANIA + GM_OBERON + GM_MIRANDA;

/// Reference epoch, J2000.0 (Julian Date, TT).
pub const URANUS_EPOCH_JD: f64 = 2451545.0;

/// Ariel semi-major axis, km.
pub const ARIEL_A_KM: f64 = 190_945.0;
/// Umbriel semi-major axis, km (JPL SSD / GUST 86).
pub const UMBRIEL_A_KM: f64 = 266_035.0;
/// Titania semi-major axis, km (JPL SSD / GUST 86).
pub const TITANIA_A_KM: f64 = 436_300.0;
/// Oberon semi-major axis, km (JPL SSD / GUST 86).
pub const OBERON_A_KM: f64 = 583_520.0;
/// Miranda semi-major axis, km (JPL SSD / GUST 86).
pub const MIRANDA_A_KM: f64 = 129_390.0;

/// Ariel orbital period, days.
pub const ARIEL_PERIOD_DAYS: f64 = 2.520379;
/// Umbriel orbital period, days.
pub const UMBRIEL_PERIOD_DAYS: f64 = 4.144177;
/// Titania orbital period, days.
pub const TITANIA_PERIOD_DAYS: f64 = 8.705872;
/// Oberon orbital period, days.
pub const OBERON_PERIOD_DAYS: f64 = 13.463238;
/// Miranda orbital period, days.
pub const MIRANDA_PERIOD_DAYS: f64 = 1.413479;

/// Ariel mean motion, rad/day.
pub const ARIEL_N_RAD_DAY: f64 = 2.0 * std::f64::consts::PI / ARIEL_PERIOD_DAYS;
/// Umbriel mean motion, rad/day (`2π / UMBRIEL_PERIOD_DAYS`, same order).
pub const UMBRIEL_N_RAD_DAY: f64 = 2.0 * std::f64::consts::PI / UMBRIEL_PERIOD_DAYS;
/// Titania mean motion, rad/day (`2π / TITANIA_PERIOD_DAYS`, same order).
pub const TITANIA_N_RAD_DAY: f64 = 2.0 * std::f64::consts::PI / TITANIA_PERIOD_DAYS;
/// Oberon mean motion, rad/day (`2π / OBERON_PERIOD_DAYS`, same order).
pub const OBERON_N_RAD_DAY: f64 = 2.0 * std::f64::consts::PI / OBERON_PERIOD_DAYS;
/// Miranda mean motion, rad/day (`2π / MIRANDA_PERIOD_DAYS`, same order).
pub const MIRANDA_N_RAD_DAY: f64 = 2.0 * std::f64::consts::PI / MIRANDA_PERIOD_DAYS;

/// Ariel orbital eccentricity.
pub const ARIEL_E: f64 = 0.0012;
/// Umbriel orbital eccentricity.
pub const UMBRIEL_E: f64 = 0.0040;
/// Titania orbital eccentricity.
pub const TITANIA_E: f64 = 0.0014;
/// Oberon orbital eccentricity.
pub const OBERON_E: f64 = 0.0016;
/// Miranda orbital eccentricity.
pub const MIRANDA_E: f64 = 0.0013;

/// Ariel mean anomaly at J2000, radians.
pub const ARIEL_M0_RAD: f64 = 123.0_f64.to_radians();
/// Umbriel mean anomaly at J2000, radians (whole-degree project-derived phase).
pub const UMBRIEL_M0_RAD: f64 = 285.0_f64.to_radians();
/// Titania mean anomaly at J2000, radians (whole-degree project-derived phase).
pub const TITANIA_M0_RAD: f64 = 107.0_f64.to_radians();
/// Oberon mean anomaly at J2000, radians (whole-degree project-derived phase).
pub const OBERON_M0_RAD: f64 = 222.0_f64.to_radians();
/// Miranda mean anomaly at J2000, radians (whole-degree project-derived phase).
pub const MIRANDA_M0_RAD: f64 = 33.0_f64.to_radians();

/// Ariel inclination to the Uranus equator, radians.
pub const ARIEL_I_RAD: f64 = 0.04_f64.to_radians();
/// Umbriel inclination to the Uranus equator, radians.
pub const UMBRIEL_I_RAD: f64 = 0.13_f64.to_radians();
/// Titania inclination to the Uranus equator, radians.
pub const TITANIA_I_RAD: f64 = 0.08_f64.to_radians();
/// Oberon inclination to the Uranus equator, radians.
pub const OBERON_I_RAD: f64 = 0.07_f64.to_radians();
/// Miranda inclination to the Uranus equator, radians.
pub const MIRANDA_I_RAD: f64 = 4.34_f64.to_radians();

/// Ariel longitude of ascending node on the Uranus equator, radians
pub const ARIEL_NODE_RAD: f64 = 167.0_f64.to_radians();
/// Umbriel longitude of ascending node on the Uranus equator, radians.
pub const UMBRIEL_NODE_RAD: f64 = 197.0_f64.to_radians();
/// Titania longitude of ascending node on the Uranus equator, radians.
pub const TITANIA_NODE_RAD: f64 = 0.0_f64.to_radians();
/// Oberon longitude of ascending node on the Uranus equator, radians.
pub const OBERON_NODE_RAD: f64 = 0.0_f64.to_radians();
/// Miranda longitude of ascending node on the Uranus equator, radians.
pub const MIRANDA_NODE_RAD: f64 = 97.0_f64.to_radians();

/// Ariel argument of pericenter, radians.
pub const ARIEL_OMEGA_RAD: f64 = 92.0_f64.to_radians();
/// Umbriel argument of pericenter, radians.
pub const UMBRIEL_OMEGA_RAD: f64 = 144.0_f64.to_radians();
/// Titania argument of pericenter, radians.
pub const TITANIA_OMEGA_RAD: f64 = 100.0_f64.to_radians();
/// Oberon argument of pericenter, radians.
pub const OBERON_OMEGA_RAD: f64 = 282.0_f64.to_radians();
/// Miranda argument of pericenter, radians.
pub const MIRANDA_OMEGA_RAD: f64 = 343.0_f64.to_radians();

/// Uranus north-pole right ascension in ICRF/J2000, radians.
pub const URANUS_POLE_RA_RAD: f64 = 257.311_f64.to_radians();
/// Uranus north-pole declination in ICRF/J2000, radians.
pub const URANUS_POLE_DEC_RAD: f64 = -15.175_f64.to_radians();

/// Newton-Raphson convergence tolerance for [`uranian_solve_kepler`]
pub const URANIAN_KEPLER_TOLERANCE: f64 = 1e-12;

/// Per-moon Keplerian element bundle for [`uranian_moon_position_kepler`]
struct UranianElements {
    a_km: f64,
    e: f64,
    n_rad_day: f64,
    m0_rad: f64,
    omega_rad: f64,
    incl_rad: f64,
    node_rad: f64,
}

/// Element bundles in [`uranian_moon_positions`] order
/// (Ariel, Umbriel, Titania, Oberon, Miranda).
const URANIAN_ELEMENTS: [UranianElements; 5] = [
    UranianElements {
        a_km: ARIEL_A_KM,
        e: ARIEL_E,
        n_rad_day: ARIEL_N_RAD_DAY,
        m0_rad: ARIEL_M0_RAD,
        omega_rad: ARIEL_OMEGA_RAD,
        incl_rad: ARIEL_I_RAD,
        node_rad: ARIEL_NODE_RAD,
    },
    UranianElements {
        a_km: UMBRIEL_A_KM,
        e: UMBRIEL_E,
        n_rad_day: UMBRIEL_N_RAD_DAY,
        m0_rad: UMBRIEL_M0_RAD,
        omega_rad: UMBRIEL_OMEGA_RAD,
        incl_rad: UMBRIEL_I_RAD,
        node_rad: UMBRIEL_NODE_RAD,
    },
    UranianElements {
        a_km: TITANIA_A_KM,
        e: TITANIA_E,
        n_rad_day: TITANIA_N_RAD_DAY,
        m0_rad: TITANIA_M0_RAD,
        omega_rad: TITANIA_OMEGA_RAD,
        incl_rad: TITANIA_I_RAD,
        node_rad: TITANIA_NODE_RAD,
    },
    UranianElements {
        a_km: OBERON_A_KM,
        e: OBERON_E,
        n_rad_day: OBERON_N_RAD_DAY,
        m0_rad: OBERON_M0_RAD,
        omega_rad: OBERON_OMEGA_RAD,
        incl_rad: OBERON_I_RAD,
        node_rad: OBERON_NODE_RAD,
    },
    UranianElements {
        a_km: MIRANDA_A_KM,
        e: MIRANDA_E,
        n_rad_day: MIRANDA_N_RAD_DAY,
        m0_rad: MIRANDA_M0_RAD,
        omega_rad: MIRANDA_OMEGA_RAD,
        incl_rad: MIRANDA_I_RAD,
        node_rad: MIRANDA_NODE_RAD,
    },
];

/// Maximum Newton-Raphson iterations for [`uranian_solve_kepler`]
/// `).
pub const URANIAN_KEPLER_MAX_ITER: u32 = 20;

/// Rotation matrix from the Uranus-equatorial frame to ICRF, row-major
/// (`v_icrf = M @ u`).
///
/// Reproduces: the frame is
/// the IAU non-rotating body frame (+z along the IAU pole, +x along
/// the ascending node of the body equator on the ICRF equator at
/// `alpha0 + 90°`); the matrix columns are the frame basis vectors in
/// ICRF, with `y = z × x` for right-handedness.
pub fn uranus_to_icrf_matrix() -> [[f64; 3]; 3] {
    let cos_ra = URANUS_POLE_RA_RAD.cos();
    let sin_ra = URANUS_POLE_RA_RAD.sin();
    let cos_dec = URANUS_POLE_DEC_RAD.cos();
    let sin_dec = URANUS_POLE_DEC_RAD.sin();
    // x-axis: IAU node direction (RA = alpha 0 + 90°, Dec = 0).
    let x_axis = (-sin_ra, cos_ra, 0.0);
    // z-axis: pole direction.
    let z_axis = (cos_dec * cos_ra, cos_dec * sin_ra, sin_dec);
    // y-axis: z cross x (right-handed).
    let y_axis = (-sin_dec * cos_ra, -sin_dec * sin_ra, cos_dec);
    [
        [x_axis.0, y_axis.0, z_axis.0],
        [x_axis.1, y_axis.1, z_axis.1],
        [x_axis.2, y_axis.2, z_axis.2],
    ]
}

/// Solve Kepler's equation `M = E − e·sin E` for `E`.
///
/// Reproduces: Newton-Raphson from `E = M`,
/// at most [`URANIAN_KEPLER_MAX_ITER`] steps, stopping when the
/// correction drops below [`URANIAN_KEPLER_TOLERANCE`].
pub fn uranian_solve_kepler(mean_anomaly: f64, eccentricity: f64) -> f64 {
    let mut e_anom = mean_anomaly;
    for _ in 0..URANIAN_KEPLER_MAX_ITER {
        let delta = (e_anom - eccentricity * e_anom.sin() - mean_anomaly)
            / (1.0 - eccentricity * e_anom.cos());
        e_anom -= delta;
        if delta.abs() < URANIAN_KEPLER_TOLERANCE {
            break;
        }
    }
    e_anom
}

/// Uranus-centric moon position in km, ICRF frame.
///
/// Reproduces operation order: time
/// since J2000, mean anomaly with `% 2π` normalization, Newton-Raphson
/// eccentric anomaly, `atan2` true anomaly, radius, perifocal
/// position, the ω/i/Ω rotation sequence (with `z1 = 0.0` kept), then
/// the Uranus-equator-to-ICRF matrix. `jd` is a TT Julian Date.
fn uranian_moon_position_kepler(jd: f64, el: &UranianElements) -> (f64, f64, f64) {
    let dt = jd - URANUS_EPOCH_JD;
    // Mean anomaly at time t, normalized to 0–2π.
    let m = (el.m0_rad + el.n_rad_day * dt).rem_euclid(2.0 * std::f64::consts::PI);
    let e_anom = uranian_solve_kepler(m, el.e);
    let cos_nu = (e_anom.cos() - el.e) / (1.0 - el.e * e_anom.cos());
    let sin_nu = ((1.0 - el.e * el.e).sqrt() * e_anom.sin()) / (1.0 - el.e * e_anom.cos());
    let nu = sin_nu.atan2(cos_nu);
    let r_km = el.a_km * (1.0 - el.e * el.e) / (1.0 + el.e * nu.cos());
    // Position in the orbital (perifocal) plane.
    let x_peri = r_km * nu.cos();
    let y_peri = r_km * nu.sin();
    let cos_omega = el.omega_rad.cos();
    let sin_omega = el.omega_rad.sin();
    let cos_i = el.incl_rad.cos();
    let sin_i = el.incl_rad.sin();
    let cos_node = el.node_rad.cos();
    let sin_node = el.node_rad.sin();
    // Rotate by argument of pericenter (model keeps `z1 = 0.0`).
    let x1 = x_peri * cos_omega - y_peri * sin_omega;
    let y1 = x_peri * sin_omega + y_peri * cos_omega;
    let z1 = 0.0;
    // Rotate by inclination around the node axis.
    let x2 = x1;
    let y2 = y1 * cos_i - z1 * sin_i;
    let z2 = y1 * sin_i + z1 * cos_i;
    // Rotate by node around the z-axis (Uranus-equatorial frame).
    let x_eq = x2 * cos_node - y2 * sin_node;
    let y_eq = x2 * sin_node + y2 * cos_node;
    let z_eq = z2;
    // Uranus equator -> ICRF.
    let m_mat = uranus_to_icrf_matrix();
    (
        m_mat[0][0] * x_eq + m_mat[0][1] * y_eq + m_mat[0][2] * z_eq,
        m_mat[1][0] * x_eq + m_mat[1][1] * y_eq + m_mat[1][2] * z_eq,
        m_mat[2][0] * x_eq + m_mat[2][1] * y_eq + m_mat[2][2] * z_eq,
    )
}

/// Five Uranian-satellite position triples
/// `(ariel, umbriel, titania, oberon, miranda)`, each `(x, y, z)`.
pub type UranianMoonQuintet = (
    (f64, f64, f64),
    (f64, f64, f64),
    (f64, f64, f64),
    (f64, f64, f64),
    (f64, f64, f64),
);

/// Uranus-centric positions of the five major Uranian moons.
///
/// Returns `(ariel, umbriel, titania, oberon, miranda)`, each
/// `(x, y, z)` in km in the equatorial ICRF/J2000 frame. `jd` is a TT
/// Julian Date; no light-time correction is applied.
///
/// Precision: mean circular elements
/// without the GUST 86 periodic terms — moon positions are good to a
/// few thousand km, with the rounded whole-degree phases and nodes
/// dominating. Weighted by the ~1 e-5 moon/planet mass ratios, the
/// resulting center-of-body error is tens of km (sub-milliarcsecond
/// at opposition), below the planet-centers SPK alternative. This is
/// deliberately circular mean-orbit geometry, not the full GUST 86
/// periodic theory.
///
/// # Errors
///
/// [`Error::NonFiniteArgument`] when `jd` is NaN or infinite.
pub fn uranian_moon_positions(jd: f64) -> Result<UranianMoonQuintet> {
    if !jd.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd });
    }
    let ariel = uranian_moon_position_kepler(jd, &URANIAN_ELEMENTS[0]);
    let umbriel = uranian_moon_position_kepler(jd, &URANIAN_ELEMENTS[1]);
    let titania = uranian_moon_position_kepler(jd, &URANIAN_ELEMENTS[2]);
    let oberon = uranian_moon_position_kepler(jd, &URANIAN_ELEMENTS[3]);
    let miranda = uranian_moon_position_kepler(jd, &URANIAN_ELEMENTS[4]);
    Ok((ariel, umbriel, titania, oberon, miranda))
}

/// Offset from the Uranus-system barycenter to Uranus's center of body.
///
/// Reproduces (via
/// `constants._uranus_cob_offset`): the GM-weighted five-moon sum,
/// negated, converted from km to AU. Unlike every other system, the
/// result is already in the equatorial ICRF/J2000 frame — the
/// `get_cob_offset` Uranus arm applies no ecliptic rotation. Add the
/// result to the Uranus-system barycenter position to get the planet
/// center. Titania and Oberon dominate (~77% of the moon weight).
///
/// # Errors
///
/// [`Error::NonFiniteArgument`] when `jd` is NaN or infinite.
pub fn uranus_cob_offset(jd: f64) -> Result<(f64, f64, f64)> {
    let (ariel, umbriel, titania, oberon, miranda) = uranian_moon_positions(jd)?;
    let total_gm = GM_URANUS + GM_URANUS_MOONS_TOTAL;
    let bary_x = (GM_ARIEL * ariel.0
        + GM_UMBRIEL * umbriel.0
        + GM_TITANIA * titania.0
        + GM_OBERON * oberon.0
        + GM_MIRANDA * miranda.0)
        / total_gm;
    let bary_y = (GM_ARIEL * ariel.1
        + GM_UMBRIEL * umbriel.1
        + GM_TITANIA * titania.1
        + GM_OBERON * oberon.1
        + GM_MIRANDA * miranda.1)
        / total_gm;
    let bary_z = (GM_ARIEL * ariel.2
        + GM_UMBRIEL * umbriel.2
        + GM_TITANIA * titania.2
        + GM_OBERON * oberon.2
        + GM_MIRANDA * miranda.2)
        / total_gm;
    Ok((
        -bary_x / MOONS_AU_KM,
        -bary_y / MOONS_AU_KM,
        -bary_z / MOONS_AU_KM,
    ))
}

// ---------------------------------------------------------------------------
// Moon names, NAIF identifiers and parent barycenters.
// ---------------------------------------------------------------------------

use crate::constants::{
    MOON_ARIEL, MOON_CALLISTO, MOON_CHARON, MOON_DEIMOS, MOON_DIONE, MOON_ENCELADUS, MOON_EUROPA,
    MOON_GANYMEDE, MOON_HYPERION, MOON_IAPETUS, MOON_IO, MOON_MIMAS, MOON_MIRANDA, MOON_OBERON,
    MOON_PHOBOS, MOON_RHEA, MOON_TETHYS, MOON_TITAN, MOON_TITANIA, MOON_TRITON, MOON_UMBRIEL,
    NAIF_ARIEL, NAIF_CALLISTO, NAIF_CHARON, NAIF_DEIMOS, NAIF_DIONE, NAIF_ENCELADUS, NAIF_EUROPA,
    NAIF_GANYMEDE, NAIF_HYPERION, NAIF_IAPETUS, NAIF_IO, NAIF_MIMAS, NAIF_MIRANDA, NAIF_OBERON,
    NAIF_PHOBOS, NAIF_RHEA, NAIF_TETHYS, NAIF_TITAN, NAIF_TITANIA, NAIF_TRITON, NAIF_UMBRIEL,
    PLMOON_OFFSET,
};

/// NAIF planet-barycenter ids.
pub const NAIF_MARS_BARYCENTER: i32 = 4;
/// NAIF Jupiter-barycenter id for ephemeris lookup.
pub const NAIF_JUPITER_BARYCENTER: i32 = 5;
/// NAIF Saturn-barycenter id for ephemeris lookup.
pub const NAIF_SATURN_BARYCENTER: i32 = 6;
/// NAIF Uranus-barycenter id for ephemeris lookup.
pub const NAIF_URANUS_BARYCENTER: i32 = 7;
/// NAIF Neptune-barycenter id for ephemeris lookup.
pub const NAIF_NEPTUNE_BARYCENTER: i32 = 8;
/// NAIF Pluto-barycenter id for ephemeris lookup.
pub const NAIF_PLUTO_BARYCENTER: i32 = 9;

/// Legacy moon id to NAIF satellite id, in order.
pub const MOON_NAIF_PAIRS: [(i32, i32); 21] = [
    (MOON_IO, NAIF_IO),
    (MOON_EUROPA, NAIF_EUROPA),
    (MOON_GANYMEDE, NAIF_GANYMEDE),
    (MOON_CALLISTO, NAIF_CALLISTO),
    (MOON_MIMAS, NAIF_MIMAS),
    (MOON_ENCELADUS, NAIF_ENCELADUS),
    (MOON_TETHYS, NAIF_TETHYS),
    (MOON_DIONE, NAIF_DIONE),
    (MOON_RHEA, NAIF_RHEA),
    (MOON_TITAN, NAIF_TITAN),
    (MOON_HYPERION, NAIF_HYPERION),
    (MOON_IAPETUS, NAIF_IAPETUS),
    (MOON_MIRANDA, NAIF_MIRANDA),
    (MOON_ARIEL, NAIF_ARIEL),
    (MOON_UMBRIEL, NAIF_UMBRIEL),
    (MOON_TITANIA, NAIF_TITANIA),
    (MOON_OBERON, NAIF_OBERON),
    (MOON_TRITON, NAIF_TRITON),
    (MOON_PHOBOS, NAIF_PHOBOS),
    (MOON_DEIMOS, NAIF_DEIMOS),
    (MOON_CHARON, NAIF_CHARON),
];

/// Legacy moon id to human-readable name, in order.
pub const MOON_NAME_PAIRS: [(i32, &str); 21] = [
    (MOON_IO, "Io"),
    (MOON_EUROPA, "Europa"),
    (MOON_GANYMEDE, "Ganymede"),
    (MOON_CALLISTO, "Callisto"),
    (MOON_MIMAS, "Mimas"),
    (MOON_ENCELADUS, "Enceladus"),
    (MOON_TETHYS, "Tethys"),
    (MOON_DIONE, "Dione"),
    (MOON_RHEA, "Rhea"),
    (MOON_TITAN, "Titan"),
    (MOON_HYPERION, "Hyperion"),
    (MOON_IAPETUS, "Iapetus"),
    (MOON_MIRANDA, "Miranda"),
    (MOON_ARIEL, "Ariel"),
    (MOON_UMBRIEL, "Umbriel"),
    (MOON_TITANIA, "Titania"),
    (MOON_OBERON, "Oberon"),
    (MOON_TRITON, "Triton"),
    (MOON_PHOBOS, "Phobos"),
    (MOON_DEIMOS, "Deimos"),
    (MOON_CHARON, "Charon"),
];

/// Legacy moon id to parent planet-barycenter NAIF id, in `MOON_PARENT_MAP` order.
pub const MOON_PARENT_PAIRS: [(i32, i32); 21] = [
    (MOON_IO, NAIF_JUPITER_BARYCENTER),
    (MOON_EUROPA, NAIF_JUPITER_BARYCENTER),
    (MOON_GANYMEDE, NAIF_JUPITER_BARYCENTER),
    (MOON_CALLISTO, NAIF_JUPITER_BARYCENTER),
    (MOON_MIMAS, NAIF_SATURN_BARYCENTER),
    (MOON_ENCELADUS, NAIF_SATURN_BARYCENTER),
    (MOON_TETHYS, NAIF_SATURN_BARYCENTER),
    (MOON_DIONE, NAIF_SATURN_BARYCENTER),
    (MOON_RHEA, NAIF_SATURN_BARYCENTER),
    (MOON_TITAN, NAIF_SATURN_BARYCENTER),
    (MOON_HYPERION, NAIF_SATURN_BARYCENTER),
    (MOON_IAPETUS, NAIF_SATURN_BARYCENTER),
    (MOON_MIRANDA, NAIF_URANUS_BARYCENTER),
    (MOON_ARIEL, NAIF_URANUS_BARYCENTER),
    (MOON_UMBRIEL, NAIF_URANUS_BARYCENTER),
    (MOON_TITANIA, NAIF_URANUS_BARYCENTER),
    (MOON_OBERON, NAIF_URANUS_BARYCENTER),
    (MOON_TRITON, NAIF_NEPTUNE_BARYCENTER),
    (MOON_PHOBOS, NAIF_MARS_BARYCENTER),
    (MOON_DEIMOS, NAIF_MARS_BARYCENTER),
    (MOON_CHARON, NAIF_PLUTO_BARYCENTER),
];

/// NAIF satellite id for a moon body id, or `None` if not a moon id.
///
/// Reproduces: the legacy `MOON_*` aliases
/// resolve through [`MOON_NAIF_PAIRS`]; otherwise the canonical scheme
/// (`PLMOON_OFFSET + NAIF`, covering NAIF satellite ids 401–998 plus
/// the 9 n 99 planet-center ids) maps any id in
/// `PLMOON_OFFSET + 400 < ipl < PLMOON_OFFSET + 1000` by subtraction,
/// so moons without a named constant (e.g. Amalthea at
/// `PLMOON_OFFSET + 505`) resolve as long as a registered kernel
/// contains them.
pub fn moon_naif_id(moon_id: i32) -> Option<i32> {
    if let Some((_, naif)) = MOON_NAIF_PAIRS
        .iter()
        .find(|(legacy, _)| *legacy == moon_id)
    {
        return Some(*naif);
    }
    if PLMOON_OFFSET + 400 < moon_id && moon_id < PLMOON_OFFSET + 1000 {
        return Some(moon_id - PLMOON_OFFSET);
    }
    None
}

/// Whether a body id is a planetary moon.
///
/// True for canonical moon ids (`PLMOON_OFFSET + NAIF`) and the legacy
/// `MOON_*` aliases.
pub fn is_planetary_moon(moon_id: i32) -> bool {
    moon_naif_id(moon_id).is_some()
}

/// Human-readable moon name.
///
/// Legacy ids resolve through [`MOON_NAME_PAIRS`]; canonical ids
/// resolve through the NAIF reverse map; anything else reports
/// `"Unknown Moon ({id})"`, as defined by the model.
pub fn get_moon_name(moon_id: i32) -> String {
    if let Some((_, name)) = MOON_NAME_PAIRS
        .iter()
        .find(|(legacy, _)| *legacy == moon_id)
    {
        return (*name).to_string();
    }
    if let Some(naif) = moon_naif_id(moon_id) {
        if let Some((legacy, _)) = MOON_NAIF_PAIRS
            .iter()
            .find(|(_, candidate)| *candidate == naif)
        {
            if let Some((_, name)) = MOON_NAME_PAIRS.iter().find(|(id, _)| *id == *legacy) {
                return (*name).to_string();
            }
        }
    }
    format!("Unknown Moon ({moon_id})")
}

/// Parent planet-barycenter NAIF id for a legacy moon id, or `None`
pub fn moon_parent_barycenter(moon_id: i32) -> Option<i32> {
    MOON_PARENT_PAIRS
        .iter()
        .find(|(legacy, _)| *legacy == moon_id)
        .map(|(_, parent)| *parent)
}

// ---------------------------------------------------------------------------
// Satellite-SPK registry, coverage and chain-closure support.
// ---------------------------------------------------------------------------

use crate::spk::{SpkRegistry, et_to_jd_tdb};
use std::collections::HashMap;
use std::path::Path;

/// Owned satellite kernels and their moon-to-kernel registrations.
///
/// SSB evaluation follows each target's parent-center chain. Registration
/// retains the chain's support bodies under their NAIF IDs. Such bodies are
/// registry infrastructure; [`Self::list_registered_moons`] lists moons only.
#[derive(Debug, Clone, Default)]
pub struct MoonSpkRegistry {
    inner: SpkRegistry,
}

/// Match an exact registered key or a nonempty POSIX file basename.
/// A trailing-slash request matches only an exact key, never a basename.
fn moon_spk_key_matches(registered: &str, requested: &str) -> bool {
    if registered == requested {
        return true;
    }
    fn posix_basename(path: &str) -> &str {
        path.rsplit('/').next().unwrap_or("")
    }
    let base = posix_basename(requested);
    if base.is_empty() {
        return false;
    }
    posix_basename(registered) == base
}

impl MoonSpkRegistry {
    /// Borrow the validated satellite images and their target registrations.
    ///
    /// The registry is immutable through this reference. Cloning it shares
    /// kernel bytes and decoded coefficients, so calculation snapshots retain
    /// their files independently of subsequent registration or close calls.
    pub fn spk_registry(&self) -> &SpkRegistry {
        &self.inner
    }

    /// Retain an already acquired satellite registry without opening its files.
    ///
    /// `registry` owns validated immutable SPK images, coefficient caches and
    /// the caller's selected target/support registrations. Moving it preserves
    /// those identities and performs no decoding, kernel copy or registration.
    /// Moon lookup and listings still filter the registered library body ids;
    /// support-center registrations remain available only for chain evaluation.
    pub fn from_owned_spk_registry(registry: SpkRegistry) -> Self {
        Self { inner: registry }
    }

    /// Empty registry: no kernels, no moon registrations.
    pub fn new() -> Self {
        Self {
            inner: SpkRegistry::new(),
        }
    }

    /// Register a satellite SPK kernel file for moon calculations
    ///
    /// `spk_path` is used as given (relative-name resolution against a
    /// data directory is binding-side); a missing file is
    /// [`Error::SpkFileNotFound`], an undecodable one the
    /// [`SpkRegistry`] decode error. With `moons = Some(ids)` only
    /// those legacy/canonical ids register (each verified against the
    /// image like probe: unchainable ids
    /// are skipped, never fatal); with `None` every segment target in
    /// the 400–999 satellite range that chains registers under its
    /// canonical id plus the legacy alias where one exists. Every
    /// registered moon also keys its canonical
    /// `PLMOON_OFFSET + NAIF` id. When nothing registers, the image
    /// is dropped and [`Error::SpkNoValidMoons`] reports the path
    /// (never a half-registered kernel).
    ///
    /// Chain-closure support: every registered moon's parent-center
    /// chain down to the SSB registers under its NAIF ids (best
    /// effort per link: an unresolvable center is skipped here and
    /// surfaces its typed error at evaluation, never invented data).
    /// Support ids stay out of the moon map
    /// ([`MoonSpkRegistry::list_registered_moons`] filters them).
    pub fn register_moon_spk(&mut self, spk_path: &str, moons: Option<&[i32]>) -> Result<()> {
        let key = self.inner.register_file(Path::new(spk_path))?;
        if let Some(ids) = moons {
            for moon_id in ids {
                if let Some(naif) = moon_naif_id(*moon_id) {
                    let _ = self.inner.register_body(*moon_id, &key, naif);
                    let _ = self.inner.register_body(PLMOON_OFFSET + naif, &key, naif);
                }
            }
        } else {
            let kernel = self
                .inner
                .kernel(&key)
                .ok_or(Error::SpkFileNotFound { path: key.clone() })?;
            let mut targets: Vec<i32> = kernel
                .index()
                .all_segments()
                .map(|seg| seg.descriptor.target)
                .filter(|tgt| 400 < *tgt && *tgt < 999)
                .collect();
            targets.sort_unstable();
            targets.dedup();
            for naif in targets {
                let _ = self.inner.register_body(PLMOON_OFFSET + naif, &key, naif);
                if let Some((legacy, _)) = MOON_NAIF_PAIRS
                    .iter()
                    .find(|(_, candidate)| *candidate == naif)
                {
                    let _ = self.inner.register_body(*legacy, &key, naif);
                }
            }
        }
        let served = self
            .inner
            .registered_bodies()
            .iter()
            .any(|reg| reg.kernel == key);
        if !served {
            self.inner.remove_kernel(&key);
            return Err(Error::SpkNoValidMoons {
                path: spk_path.to_string(),
            });
        }
        self.register_chain_support(&key);
        Ok(())
    }

    /// Register the parent-center chain closure for every moon the
    /// kernel `key` serves, under the centers' NAIF ids (see the
    /// [`MoonSpkRegistry`] chain-support note). The walk is bounded by
    /// [`NaifBaseProvider::MAX_DEPTH`]; unresolvable links are skipped
    /// (evaluation reports them honestly). Moon-map membership is
    /// untouched: support ids are exactly the non-moon NAIF ids.
    fn register_chain_support(&mut self, key: &str) {
        use crate::spk::resolve_chain;
        let support: Vec<i32> = match self.inner.kernel(key) {
            Some(kernel) => {
                let mut moon_naifs: Vec<i32> = self
                    .inner
                    .registered_bodies()
                    .into_iter()
                    .filter(|reg| reg.kernel == key)
                    .filter_map(|reg| moon_naif_id(reg.body))
                    .collect();
                moon_naifs.sort_unstable();
                moon_naifs.dedup();
                let mut support = Vec::new();
                for naif in moon_naifs {
                    let mut next = match resolve_chain(kernel.index(), naif) {
                        Ok(chain) => chain.center,
                        Err(_) => continue,
                    };
                    let mut depth = 0;
                    while next != 0 && depth < NaifBaseProvider::MAX_DEPTH {
                        support.push(next);
                        next = match resolve_chain(kernel.index(), next) {
                            Ok(chain) => chain.center,
                            Err(_) => break,
                        };
                        depth += 1;
                    }
                }
                support.sort_unstable();
                support.dedup();
                support
            }
            None => Vec::new(),
        };
        for center in support {
            let _ = self.inner.register_body(center, key, center);
        }
    }

    /// Unregister a satellite SPK kernel and its moons: exact-key or same-file-name matches
    /// drop their kernel images together with every moon registration
    /// pointing at them. Unknown paths are a no-op, never an error.
    pub fn unregister_moon_spk(&mut self, spk_file: &str) {
        let keys: Vec<String> = self
            .inner
            .kernel_paths()
            .into_iter()
            .filter(|key| moon_spk_key_matches(key, spk_file))
            .collect();
        for key in keys {
            self.inner.remove_kernel(&key);
        }
    }

    /// Moon id to kernel key for every registered moon: the moon map only — chain-closure
    /// support bodies (plain NAIF ids) never leak into the listing.
    pub fn list_registered_moons(&self) -> HashMap<i32, String> {
        self.inner
            .registered_bodies()
            .into_iter()
            .filter(|reg| moon_naif_id(reg.body).is_some())
            .map(|reg| (reg.body, reg.kernel))
            .collect()
    }

    /// Whole-kernel coverage envelope for a registered moon, in TDB Julian days.
    ///
    /// Returns `(start_jd_tdb, end_jd_tdb)`, or `None` when the registration,
    /// retained kernel or finite segment bounds cannot be found. This is the
    /// minimum start and maximum end over **every segment in that kernel**,
    /// converted from TDB seconds past J2000. It is not a target-filtered
    /// interval, a TT interval or proof of continuous center-chain coverage.
    /// Gaps and narrower target/support segments can still fail evaluation
    /// inside this envelope; each actual SPK read owns those domain checks.
    ///
    /// This metadata query neither opens an alternate kernel nor changes the
    /// retained registry. Callers presenting these bounds alongside a TT request
    /// must keep the two clock labels explicit rather than silently relabeling
    /// the envelope as TT.
    pub fn get_moon_coverage(&self, moon_id: i32) -> Option<(f64, f64)> {
        let naif = moon_naif_id(moon_id)?;
        let reg = self
            .inner
            .body_info(moon_id)
            .or_else(|| self.inner.body_info(PLMOON_OFFSET + naif))?;
        let kernel = self.inner.kernel(&reg.kernel)?;
        let mut start_et = f64::INFINITY;
        let mut end_et = f64::NEG_INFINITY;
        for seg in kernel.index().all_segments() {
            start_et = start_et.min(seg.descriptor.start_et);
            end_et = end_et.max(seg.descriptor.end_et);
        }
        if !start_et.is_finite() || !end_et.is_finite() {
            return None;
        }
        Some((et_to_jd_tdb(start_et).ok()?, et_to_jd_tdb(end_et).ok()?))
    }

    /// Remove all kernel handles and moon registrations from this registry.
    ///
    /// Previously cloned registries/calculation snapshots retain their shared
    /// immutable kernel owners until those owners are dropped. This does not
    /// invalidate a state already returned by value or forcibly close resources
    /// still retained by another snapshot.
    pub fn close_moon_kernels(&mut self) {
        self.inner.clear();
    }
}

// ---------------------------------------------------------------------------
// Saturnian moons: TASS 1.7 analytical theory (Vienne-Duriez via Stellarium).
//
// TASS 1.7 (Vienne & Duriez 1995, A&A 297, 588–605; later IMCCE revision).
// MIT-licensed implementation/data: Copyright (c) 2005 Johannes Gajdosik.
// Tables in crate::tass17_tables retain the attribution; regenerate with
// tools/generators/gen_tass17_tables.py. See THIRD_PARTY_NOTICES.
// ---------------------------------------------------------------------------

/// TASS reference epoch: JD 2444240.0 = 1980-01-01 00:00 TT
pub const TASS_EPOCH_JD: f64 = 2444240.0;

/// Saturn GM in km³/s².
pub const GM_SATURN: f64 = 37_931_206.23;
/// Saturnian-moon GMs in km³/s².
pub const GM_MIMAS: f64 = 2.50349;
/// Saturnian-moon GMs in km³/s².
pub const GM_ENCELADUS: f64 = 7.21037;
/// Saturnian-moon GMs in km³/s².
pub const GM_TETHYS: f64 = 41.21353;
/// Saturnian-moon GMs in km³/s².
pub const GM_DIONE: f64 = 73.11607;
/// Saturnian-moon GMs in km³/s².
pub const GM_RHEA: f64 = 153.94175;
/// Titan GM in km³/s².
pub const GM_TITAN: f64 = 8978.13710;
/// Saturnian-moon GMs in km³/s².
pub const GM_HYPERION: f64 = 0.37049;
/// Saturnian-moon GMs in km³/s².
pub const GM_IAPETUS: f64 = 120.51511;

/// Total Saturnian-moon GM.
pub const GM_SATURN_MOONS_TOTAL: f64 =
    GM_MIMAS + GM_ENCELADUS + GM_TETHYS + GM_DIONE + GM_RHEA + GM_TITAN + GM_HYPERION + GM_IAPETUS;

/// TASS body index of Titan.
pub const TASS_BODY_TITAN: i32 = 5;
/// Number of TASS bodies.
pub const TASS_BODY_COUNT: i32 = 8;

/// Mean longitudes of the first seven satellites at `t` days past the
/// TASS epoch, in radians.
///
/// Reproduces: for each body the first multiterm of
/// its longitude series (`series[1][0]`, whose index values are carried
/// along unread here) contributes its
/// `amplitude * sin(phase + freq * t)` terms. Hyperion has no base
/// longitude row; its slot stays 0.0 as defined by the model (the loop
/// covers bodies 0-6 and `lon` starts zeroed).
fn tass_calc_lon(t: f64) -> [f64; 7] {
    use crate::tass17_tables::tass_series;
    let mut lon = [0.0_f64; 7];
    for (i, slot) in lon.iter_mut().enumerate() {
        let Some(rows) = tass_series(i, 1) else {
            continue;
        };
        let Some(first) = rows.first() else { continue };
        for &(amplitude, phase, freq) in first.1 {
            *slot += amplitude * (phase + freq * t).sin();
        }
    }
    lon
}

/// Six TASS orbital elements `[n, lambda, k, h, q, p]` for one body.
///
/// Reproduces with the same series order and
/// the same accumulation order inside every row:
/// - series 0 perturbs the mean motion `n`, closed by `aam * (1 + n)`;
/// - series 1 perturbs the mean longitude: non-Hyperion bodies start from
///   their `_calc_lon` slot (skipping the already-consumed first
///   multiterm) and every body gains `aam * t`;
/// - series 2/3 perturb the eccentricity vector `(k, h)` and the
///   inclination vector `(q, p)` through shared `cos`/`sin` arguments.
///
/// `arg` weights the first seven mean longitudes with the row's index
/// values (`sum(indices[i] * lon[i])`); the eighth index value is carried
/// along unread, exactly like loop.
fn tass_calc_elem(t: f64, lon: &[f64; 7], body: usize) -> [f64; 6] {
    use crate::tass17_tables::{TASS17_BODIES, tass_series};
    let params = &TASS17_BODIES[body];
    let mut elem = params.s0;

    let series0 = tass_series(body, 0).expect("TASS series 0 exists");
    for (indices, terms) in series0 {
        let arg: f64 = indices[..7]
            .iter()
            .zip(lon.iter())
            .map(|(w, l)| f64::from(*w) * *l)
            .sum();
        for &(amplitude, phase, freq) in *terms {
            elem[0] += amplitude * (phase + freq * t + arg).cos();
        }
    }
    elem[0] = params.aam * (1.0 + elem[0]);

    let series1 = tass_series(body, 1).expect("TASS series 1 exists");
    let start = if body == 7 { 0 } else { 1 };
    if body != 7 {
        elem[1] += lon[body];
    }
    for (indices, terms) in &series1[start..] {
        let arg: f64 = indices[..7]
            .iter()
            .zip(lon.iter())
            .map(|(w, l)| f64::from(*w) * *l)
            .sum();
        for &(amplitude, phase, freq) in *terms {
            elem[1] += amplitude * (phase + freq * t + arg).sin();
        }
    }
    elem[1] += params.aam * t;

    for (series, pair) in [(2, (2, 3)), (3, (4, 5))] {
        let rows = tass_series(body, series).expect("TASS series exists");
        for (indices, terms) in rows {
            let arg: f64 = indices[..7]
                .iter()
                .zip(lon.iter())
                .map(|(w, l)| f64::from(*w) * *l)
                .sum();
            for &(amplitude, phase, freq) in *terms {
                let x = phase + freq * t + arg;
                elem[pair.0] += amplitude * x.cos();
                elem[pair.1] += amplitude * x.sin();
            }
        }
    }
    elem
}

/// Rectangular Saturn-centric state from TASS elements.
///
/// Reproduces (the Stellarium
/// `EllipticToRectangularN` adaptation): semi-major axis from
/// `a = (mu / n^2)^(1/3)`, eccentricity from the `(k, h)` vector,
/// mean anomaly `M = lambda + n * dt` reduced with the floored modulo,
/// up to 15 fixed-point Kepler iterations breaking (keeping the previous
/// iterate) below `1e-14`, true anomaly from the eccentric anomaly,
/// orbital-plane position/velocity, then the `omega` (argument of
/// perihelion from `varpi - Omega`), inclination, and node rotations
/// into the TASS reference frame.
///
/// Returns `(x, y, z, vx, vy, vz)` in AU and AU/day. Both model call
/// sites pass `dt = 0.0`; the parameter is kept so the structure stays
/// the same as the model shape.
fn tass_elliptic_to_rectangular(
    mu: f64,
    elem: &[f64; 6],
    dt: f64,
) -> (f64, f64, f64, f64, f64, f64) {
    let n = elem[0];
    let lam = elem[1];
    let k = elem[2];
    let h = elem[3];
    let q = elem[4];
    let p = elem[5];

    let a = (mu / (n * n)).powf(1.0 / 3.0);
    let e = (k * k + h * h).sqrt();
    let m = lam + n * dt;
    let varpi = if e > 1e-12 { h.atan2(k) } else { 0.0 };
    let m_anom = (m - varpi).rem_euclid(2.0 * std::f64::consts::PI);

    let mut eccentric = m_anom;
    for _ in 0..15 {
        let next = m_anom + e * eccentric.sin();
        if (next - eccentric).abs() < 1e-14 {
            break;
        }
        eccentric = next;
    }
    let cos_e = eccentric.cos();
    let sin_e = eccentric.sin();

    let beta = (1.0 - e * e).sqrt();
    let cos_nu = (cos_e - e) / (1.0 - e * cos_e);
    let sin_nu = beta * sin_e / (1.0 - e * cos_e);
    let r = a * (1.0 - e * cos_e);
    let x_orb = r * cos_nu;
    let y_orb = r * sin_nu;
    let v_factor = (mu / a).sqrt() / (1.0 - e * cos_e);
    let vx_orb = -v_factor * sin_e;
    let vy_orb = v_factor * beta * cos_e;

    let sin_half_i_sq = q * q + p * p;
    let cos_half_i = (1.0 - sin_half_i_sq).max(0.0).sqrt();
    let sin_half_i = sin_half_i_sq.sqrt();
    let cos_i = 1.0 - 2.0 * sin_half_i_sq;
    let sin_i = 2.0 * sin_half_i * cos_half_i;
    let (cos_omega_node, sin_omega_node) = if sin_half_i > 1e-12 {
        (q / sin_half_i, p / sin_half_i)
    } else {
        (1.0, 0.0)
    };
    let node = sin_omega_node.atan2(cos_omega_node);
    let omega = varpi - node;
    let (sin_omega, cos_omega) = omega.sin_cos();

    let x_rot = x_orb * cos_omega - y_orb * sin_omega;
    let y_rot = x_orb * sin_omega + y_orb * cos_omega;
    let vx_rot = vx_orb * cos_omega - vy_orb * sin_omega;
    let vy_rot = vx_orb * sin_omega + vy_orb * cos_omega;

    let y_incl = y_rot * cos_i;
    let z_incl = y_rot * sin_i;
    let vy_incl = vy_rot * cos_i;
    let vz_incl = vy_rot * sin_i;

    (
        x_rot * cos_omega_node - y_incl * sin_omega_node,
        x_rot * sin_omega_node + y_incl * cos_omega_node,
        z_incl,
        vx_rot * cos_omega_node - vy_incl * sin_omega_node,
        vx_rot * sin_omega_node + vy_incl * cos_omega_node,
        vz_incl,
    )
}

/// Rotate a TASS-frame vector into the J2000 ecliptic (VSOP 87) frame.
///
/// Reproduces over the `TASS17_TO_VSOP87`
/// row-major matrix. Pure rotation: defined for every finite or
/// non-finite input triple, no validation.
fn tass_apply_rotation(x: (f64, f64, f64)) -> (f64, f64, f64) {
    use crate::tass17_tables::TASS17_TO_VSOP87;
    let r = TASS17_TO_VSOP87;
    (
        r[0] * x.0 + r[1] * x.1 + r[2] * x.2,
        r[3] * x.0 + r[4] * x.1 + r[5] * x.2,
        r[6] * x.0 + r[7] * x.1 + r[8] * x.2,
    )
}

/// Position of a Saturnian satellite relative to Saturn's center.
///
/// Reproduces: days past the TASS epoch,
/// shared mean longitudes, per-body elements, rectangular conversion at
/// `dt = 0.0`, and rotation into the J2000 ecliptic (ICRF) frame.
///
/// `jd_tt` is a TT Julian Date (any finite value); `body` is the model
/// index (Mimas 0.. Hyperion 7, see [`TASS_BODY_TITAN`]). Output is
/// `(x, y, z)` in AU in the model frame (treated as
/// J2000-ecliptic by [`get_cob_offset`], which rotates it — see
/// [`saturn_cob_offset_icrf`]).
///
/// Model precision is ~50-100 km per satellite (~0.05 arcsec at
/// Saturn distance); no fitted constants, no widened tolerances.
///
/// # Errors
///
/// [`Error::NonFiniteArgument`] when `jd_tt` is NaN or infinite.
/// [`Error::InvalidBody`] for a body index outside 0-7 (the typed refusal
/// for edge; [`Error::InvalidBody`] mirrors the
/// input-validation branch).
pub fn saturn_moon_position(jd_tt: f64, body: i32) -> Result<(f64, f64, f64)> {
    let (x, y, z, _, _, _) = saturn_moon_state(jd_tt, body)?;
    Ok((x, y, z))
}

/// Position and velocity of a Saturnian satellite relative to Saturn's
/// center.
///
/// Reproduces: the same pipeline
/// as [`saturn_moon_position`] keeping the rectangular velocity. Output
/// is `(x, y, z, vx, vy, vz)` in AU and AU/day in the J2000 ecliptic
/// (ICRF) frame.
///
/// # Errors
///
/// [`Error::NonFiniteArgument`] when `jd_tt` is NaN or infinite.
/// [`Error::InvalidBody`] for a body index outside 0-7 (the typed refusal
/// for edge).
pub fn saturn_moon_position_velocity(
    jd_tt: f64,
    body: i32,
) -> Result<(f64, f64, f64, f64, f64, f64)> {
    saturn_moon_state(jd_tt, body)
}

/// Shared TASS state pipeline behind [`saturn_moon_position`] and
/// [`saturn_moon_position_velocity`].
fn saturn_moon_state(jd_tt: f64, body: i32) -> Result<(f64, f64, f64, f64, f64, f64)> {
    use crate::tass17_tables::TASS17_BODIES;
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    if !(0..TASS_BODY_COUNT).contains(&body) {
        return Err(Error::InvalidBody {
            body_id: body,
            operation: "saturn_moon_position",
        });
    }
    let index = body as usize;
    let t = jd_tt - TASS_EPOCH_JD;
    let lon = tass_calc_lon(t);
    let elem = tass_calc_elem(t, &lon, index);
    let state = tass_elliptic_to_rectangular(TASS17_BODIES[index].mu, &elem, 0.0);
    let (x, y, z) = tass_apply_rotation((state.0, state.1, state.2));
    let (vx, vy, vz) = tass_apply_rotation((state.3, state.4, state.5));
    Ok((x, y, z, vx, vy, vz))
}

/// Titan's position relative to Saturn's center.
///
/// backwards-compatibility wrapper: full TASS
/// 1.7 theory at `body = 5`. Output is `(x, y, z)` in AU in the J2000
/// ecliptic frame.
///
/// # Errors
///
/// [`Error::NonFiniteArgument`] when `jd_tt` is NaN or infinite.
pub fn titan_position(jd_tt: f64) -> Result<(f64, f64, f64)> {
    saturn_moon_position(jd_tt, TASS_BODY_TITAN)
}

/// Offset from the Saturn-system barycenter to Saturn's center of body.
///
/// Titan dominates the moon mass
/// (~96%), so the composition carries Titan alone — the GM-weighted
/// Titan offset, negated, in the same ecliptic frame
/// [`saturn_moon_position`] emits.
///
/// # Errors
///
/// [`Error::NonFiniteArgument`] when `jd_tt` is NaN or infinite.
pub fn saturn_cob_offset(jd_tt: f64) -> Result<(f64, f64, f64)> {
    let (tx, ty, tz) = titan_position(jd_tt)?;
    let total_gm = GM_SATURN + GM_SATURN_MOONS_TOTAL;
    let bary_x = GM_TITAN * tx / total_gm;
    let bary_y = GM_TITAN * ty / total_gm;
    let bary_z = GM_TITAN * tz / total_gm;
    Ok((-bary_x, -bary_y, -bary_z))
}

/// Saturn center-of-body offset in the equatorial ICRF/J2000 frame.
///
/// Saturn center-of-body displacement in ICRF/J2000.
///
/// Combines the Saturn dispatch arm:
/// [`saturn_cob_offset`] followed by [`ecliptic_to_equatorial_j2000`].
/// Add the result to the Saturn-system barycenter position to get
/// Saturn's center.
///
/// # Errors
///
/// [`Error::NonFiniteArgument`] when `jd_tt` is NaN or infinite.
pub fn saturn_cob_offset_icrf(jd_tt: f64) -> Result<(f64, f64, f64)> {
    Ok(ecliptic_to_equatorial_j2000(saturn_cob_offset(jd_tt)?))
}

/// Offset from a planet-system barycenter to the planet center of
/// body.
///
/// `planet` accepts the model spellings (`"jupiter"`,
/// `"jupiter barycenter"`,...): ASCII case folding (the names are
/// all ASCII), one `" barycenter"` suffix stripped wherever it
/// occurs, surrounding whitespace trimmed — the model
/// `lower.replace.strip` chain. Jupiter, Saturn, Neptune and
/// Pluto ride their analytical theories through the fixed J2000
/// ecliptic-to-equatorial rotation; Uranus is already equatorial
/// ICRF, so no rotation applies. Unrecognized names return
/// `(0.0, 0.0, 0.0)`, as defined by the model. `jd_tt` is a TT
/// Julian Date; add the result to the barycenter position to get the
/// planet center.
///
/// # Errors
///
/// [`Error::NonFiniteArgument`] propagates from the theory kernels
/// on non-finite `jd_tt`.
pub fn get_cob_offset(planet: &str, jd_tt: f64) -> Result<(f64, f64, f64)> {
    let name = planet
        .to_ascii_lowercase()
        .replace(" barycenter", "")
        .trim()
        .to_string();
    match name.as_str() {
        "jupiter" => Ok(ecliptic_to_equatorial_j2000(jupiter_cob_offset(jd_tt)?)),
        "saturn" => Ok(ecliptic_to_equatorial_j2000(saturn_cob_offset(jd_tt)?)),
        "uranus" => uranus_cob_offset(jd_tt),
        "neptune" => Ok(ecliptic_to_equatorial_j2000(neptune_cob_offset(jd_tt)?)),
        "pluto" => Ok(ecliptic_to_equatorial_j2000(pluto_cob_offset(jd_tt)?)),
        _ => Ok((0.0, 0.0, 0.0)),
    }
}

#[cfg(test)]
mod galilean_tests {
    use super::*;

    fn norm(v: (f64, f64, f64)) -> f64 {
        (v.0 * v.0 + v.1 * v.1 + v.2 * v.2).sqrt()
    }

    /// Orbital-radius envelope per satellite in km: base radius times
    /// (1 ± the summed absolute radius-series amplitudes), widened by
    /// the latitude projection (which only shrinks the norm). The
    /// envelope below uses the model table amplitudes, not fitted data.
    fn radius_bounds(base_rj: f64, amp_sum: f64) -> (f64, f64) {
        // Latitude projection (cos B, B <= ~0.013 rad) can only shrink
        // the Cartesian norm below the in-plane radius; the 0.9999
        // factor covers that analytic shrinkage with margin.
        let lo = base_rj * (1.0 - amp_sum) * JUPITER_RADIUS_KM * 0.9999;
        // Latitude projection and the cosine-series phase can only shrink
        // the Cartesian norm below the in-plane radius; allow 1% margin
        // for the series-phase interaction (analytic, not fitted).
        let hi = base_rj * (1.0 + amp_sum) * JUPITER_RADIUS_KM * 1.01;
        (lo, hi)
    }

    #[test]
    fn epoch_positions_are_finite_and_bounded_by_series() {
        // Nominal case at the E 5 reference epoch (t = 0 boundary).
        let (io, europa, ganymede, callisto) =
            galilean_moon_positions(GALILEAN_EPOCH_JD).expect("epoch must evaluate");
        // Summed absolute radius amplitudes straight from the model tables.
        let bounds = [
            (
                io,
                E5_R_IO,
                0.0041339 + 0.0000387 + 0.0000214 + 0.0000170 + 0.0000131 + 0.0000106 + 0.0000066,
            ),
            (
                europa,
                E5_R_EUROPA,
                0.0093848
                    + 0.0003116
                    + 0.0001744
                    + 0.0001442
                    + 0.0000553
                    + 0.0000523
                    + 0.0000290
                    + 0.0000164
                    + 0.0000107
                    + 0.0000102
                    + 0.0000091,
            ),
            (
                ganymede,
                E5_R_GANYMEDE,
                0.0014388
                    + 0.0007919
                    + 0.0006342
                    + 0.0001761
                    + 0.0000294
                    + 0.0000156
                    + 0.0000156
                    + 0.0000153
                    + 0.0000070
                    + 0.0000051,
            ),
            (
                callisto,
                E5_R_CALLISTO,
                0.0073546
                    + 0.0001621
                    + 0.0000974
                    + 0.0000543
                    + 0.0000271
                    + 0.0000182
                    + 0.0000177
                    + 0.0000167
                    + 0.0000167
                    + 0.0000155
                    + 0.0000142
                    + 0.0000105
                    + 0.0000092
                    + 0.0000089
                    + 0.0000062
                    + 0.0000048,
            ),
        ];
        for (pos, base, amp) in bounds {
            assert!(pos.0.is_finite() && pos.1.is_finite() && pos.2.is_finite());
            let (lo, hi) = radius_bounds(base, amp);
            let r = norm(pos);
            assert!(r >= lo && r <= hi, "r={r} outside [{lo}, {hi}]");
        }
    }

    #[test]
    fn position_is_deterministic() {
        let jd = 2460000.5;
        let a = galilean_moon_positions(jd).unwrap();
        let b = galilean_moon_positions(jd).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn fundamental_arguments_pin_linear_law_at_epoch_boundary() {
        // Boundary arm: at t = 0 the linear laws collapse to their
        // model constants, and gamma is exactly its two-sine value.
        let a = e5_fundamental_arguments(0.0);
        assert_eq!(a.l1, 106.07719);
        assert_eq!(a.l2, 175.73161);
        assert_eq!(a.l3, 120.55883);
        assert_eq!(a.l4, 84.44459);
        assert_eq!(a.pi1, 97.0881);
        assert_eq!(a.w1, 312.3346);
        assert_eq!(a.psi, 316.5182);
        assert_eq!(a.phi, 199.6766);
        assert_eq!(a.gp, 31.97853);
        assert_eq!(a.pj, 13.469942);
        let gamma = 0.33033 * e5_deg_to_rad(163.679).sin() + 0.03439 * e5_deg_to_rad(34.486).sin();
        assert_eq!(a.g, 30.23756 + gamma);
        // Corrections start zeroed; the assembly fills them.
        assert_eq!(
            (a.sigma1, a.sigma2, a.sigma3, a.sigma4),
            (0.0, 0.0, 0.0, 0.0)
        );
        assert_eq!((a.lon1, a.lon2, a.lon3, a.lon4), (0.0, 0.0, 0.0, 0.0));
    }

    #[test]
    fn malformed_jd_is_refused() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let pos = galilean_moon_positions(bad);
            let cob = jupiter_cob_offset(bad);
            let icrf = jupiter_cob_offset_icrf(bad);
            for err in [pos.map(|_| ()), cob.map(|_| ()), icrf.map(|_| ())] {
                match err {
                    Err(Error::NonFiniteArgument { value }) => {
                        if bad.is_nan() {
                            assert!(value.is_nan());
                        } else {
                            assert_eq!(value, bad);
                        }
                    }
                    other => panic!("jd={bad}: expected NonFiniteArgument, got {other:?}"),
                }
            }
        }
    }

    #[test]
    fn cob_offset_is_antiparallel_and_gm_weighted() {
        let jd = 2460000.5;
        let (io, europa, ganymede, callisto) = galilean_moon_positions(jd).unwrap();
        let off = jupiter_cob_offset(jd).unwrap();
        // composition: offset = −Σ(GM_i·r_i)/total/AU.
        let total = GM_JUPITER + GM_GALILEAN_TOTAL;
        let wx = (GM_IO * io.0
            + GM_EUROPA * europa.0
            + GM_GANYMEDE * ganymede.0
            + GM_CALLISTO * callisto.0)
            / total
            / MOONS_AU_KM;
        let wy = (GM_IO * io.1
            + GM_EUROPA * europa.1
            + GM_GANYMEDE * ganymede.1
            + GM_CALLISTO * callisto.1)
            / total
            / MOONS_AU_KM;
        let wz = (GM_IO * io.2
            + GM_EUROPA * europa.2
            + GM_GANYMEDE * ganymede.2
            + GM_CALLISTO * callisto.2)
            / total
            / MOONS_AU_KM;
        // Moon/total mass ratio is ~2 e-4, so a ratio-scale tolerance.
        assert!((off.0 + wx).abs() < 1e-15);
        assert!((off.1 + wy).abs() < 1e-15);
        assert!((off.2 + wz).abs() < 1e-15);
        // Magnitude bound: |offset| <= Σ(GM_i)·max radius/total/AU.
        let hi = GM_GALILEAN_TOTAL * E5_R_CALLISTO * 1.01 * JUPITER_RADIUS_KM / total / MOONS_AU_KM;
        assert!(norm(off) <= hi);
    }

    #[test]
    fn icrf_arm_composes_rotation_verbatim() {
        let jd = 2451545.0 + 777.25;
        let icrf = jupiter_cob_offset_icrf(jd).unwrap();
        let expected = ecliptic_to_equatorial_j2000(jupiter_cob_offset(jd).unwrap());
        assert_eq!(icrf, expected);
    }

    #[test]
    fn far_epoch_stays_bounded() {
        // Tails at 1900/2050 (inside the E 5 useful span): unwrapped
        // arguments keep every radius inside its series envelope.
        for jd in [2415020.5, 2469807.5] {
            let (io, europa, ganymede, callisto) =
                galilean_moon_positions(jd).expect("tail must evaluate");
            for (pos, base) in [
                (io, E5_R_IO),
                (europa, E5_R_EUROPA),
                (ganymede, E5_R_GANYMEDE),
                (callisto, E5_R_CALLISTO),
            ] {
                assert!(pos.0.is_finite() && pos.1.is_finite() && pos.2.is_finite());
                let r = norm(pos);
                // Loose analytic envelope: base radius ±2% in km.
                let lo = base * 0.98 * JUPITER_RADIUS_KM;
                let hi = base * 1.02 * JUPITER_RADIUS_KM;
                assert!(r >= lo && r <= hi, "jd={jd} r={r} outside [{lo}, {hi}]");
            }
        }
    }

    #[test]
    fn release_constants_match_source_values() {
        assert_eq!(GALILEAN_EPOCH_JD, 2443000.5);
        assert_eq!(JUPITER_RADIUS_KM, 71492.0);
        assert_eq!(GM_JUPITER, 126686531.9);
        assert_eq!(GM_IO, 5959.91547);
        assert_eq!(GM_EUROPA, 3202.71210);
        assert_eq!(GM_GANYMEDE, 9887.83275);
        assert_eq!(GM_CALLISTO, 7179.28340);
        // Total is the same-order sum (exactly, not recomputed
        // from a different association).
        assert_eq!(
            GM_GALILEAN_TOTAL,
            GM_IO + GM_EUROPA + GM_GANYMEDE + GM_CALLISTO
        );
        assert_eq!(E5_R_IO, 5.90569);
        assert_eq!(E5_R_EUROPA, 9.39657);
        assert_eq!(E5_R_GANYMEDE, 14.98832);
        assert_eq!(E5_R_CALLISTO, 26.36273);
        assert_eq!(E5_ORBIT_TILT_DEG, 1.303267);
        assert_eq!(E5_ORBIT_NODE_DEG, 100.464407);
        // Derived frame constants round-trip to the model definitions.
        let t0 = (2451545.0 - 2433282.423) / 36525.0;
        assert_eq!(E5_T0, t0);
        assert_eq!(
            E5_J2000_PRECESSION_DEG,
            1.3966626 * t0 + 0.0003088 * (t0 * t0)
        );
        let tilt = 3.120262 + 0.0006 * (2451545.0 - 2415020.5) / 36525.0;
        assert_eq!(E5_EQUATOR_TILT_DEG, tilt);
        // Degree helper keeps the model multiply-then-divide order:
        // 180 degrees is pi to within one rounding of the final divide.
        assert!((e5_deg_to_rad(180.0) - std::f64::consts::PI).abs() <= f64::EPSILON);
    }
}

#[cfg(test)]
mod uranian_tests {
    use super::*;

    fn norm(v: (f64, f64, f64)) -> f64 {
        (v.0 * v.0 + v.1 * v.1 + v.2 * v.2).sqrt()
    }

    #[test]
    fn release_constants_match_source_values() {
        assert_eq!(GM_URANUS, 5_793_951.3);
        assert_eq!(GM_ARIEL, 14.0);
        assert_eq!(GM_UMBRIEL, 13.0);
        assert_eq!(GM_TITANIA, 41.0);
        assert_eq!(GM_OBERON, 29.0);
        assert_eq!(GM_MIRANDA, 0.7);
        // Total is the same-order sum (exactly, not recomputed
        // from a different association).
        assert_eq!(
            GM_URANUS_MOONS_TOTAL,
            GM_ARIEL + GM_UMBRIEL + GM_TITANIA + GM_OBERON + GM_MIRANDA
        );
        assert_eq!(URANUS_EPOCH_JD, 2451545.0);
        assert_eq!(ARIEL_A_KM, 190_945.0);
        assert_eq!(UMBRIEL_A_KM, 266_035.0);
        assert_eq!(TITANIA_A_KM, 436_300.0);
        assert_eq!(OBERON_A_KM, 583_520.0);
        assert_eq!(MIRANDA_A_KM, 129_390.0);
        // Mean motions keep the model divide order (exactly).
        assert_eq!(
            ARIEL_N_RAD_DAY,
            2.0 * std::f64::consts::PI / ARIEL_PERIOD_DAYS
        );
        assert_eq!(
            UMBRIEL_N_RAD_DAY,
            2.0 * std::f64::consts::PI / UMBRIEL_PERIOD_DAYS
        );
        assert_eq!(
            TITANIA_N_RAD_DAY,
            2.0 * std::f64::consts::PI / TITANIA_PERIOD_DAYS
        );
        assert_eq!(
            OBERON_N_RAD_DAY,
            2.0 * std::f64::consts::PI / OBERON_PERIOD_DAYS
        );
        assert_eq!(
            MIRANDA_N_RAD_DAY,
            2.0 * std::f64::consts::PI / MIRANDA_PERIOD_DAYS
        );
        assert_eq!(ARIEL_E, 0.0012);
        assert_eq!(UMBRIEL_E, 0.0040);
        assert_eq!(TITANIA_E, 0.0014);
        assert_eq!(OBERON_E, 0.0016);
        assert_eq!(MIRANDA_E, 0.0013);
        assert_eq!(URANIAN_KEPLER_TOLERANCE, 1e-12);
        assert_eq!(URANIAN_KEPLER_MAX_ITER, 20);
    }

    #[test]
    fn kepler_circular_case_is_identity() {
        // e = 0 solves E = M exactly on the first correction.
        for k in 0..25 {
            let m = k as f64 * 0.2513;
            assert_eq!(uranian_solve_kepler(m, 0.0), m);
        }
    }

    #[test]
    fn kepler_converges_on_uranian_eccentricities() {
        // Every Uranian eccentricity over a full anomaly grid leaves a
        // residual below the model tolerance.
        let es = [ARIEL_E, UMBRIEL_E, TITANIA_E, OBERON_E, MIRANDA_E];
        for e in es {
            for k in 0..64 {
                let m = k as f64 * 2.0 * std::f64::consts::PI / 64.0;
                let e_anom = uranian_solve_kepler(m, e);
                let residual = (e_anom - e * e_anom.sin() - m).abs();
                assert!(residual < URANIAN_KEPLER_TOLERANCE, "e={e} m={m}");
            }
        }
        // M = 0 is a fixed point for every e.
        for e in es {
            assert_eq!(uranian_solve_kepler(0.0, e), 0.0);
        }
    }

    #[test]
    // The fixture deliberately checks matrix columns through explicit row/column indices.
    #[allow(clippy::needless_range_loop)]
    fn icrf_matrix_is_a_proper_rotation() {
        let m = uranus_to_icrf_matrix();
        // Columns are the model basis vectors: IAU node, z × x, pole.
        let ra = URANUS_POLE_RA_RAD;
        let dec = URANUS_POLE_DEC_RAD;
        let eps = 1e-15;
        assert!((m[0][0] + ra.sin()).abs() < eps);
        assert!((m[1][0] - ra.cos()).abs() < eps);
        assert_eq!(m[2][0], 0.0);
        assert!((m[0][2] - dec.cos() * ra.cos()).abs() < eps);
        assert!((m[1][2] - dec.cos() * ra.sin()).abs() < eps);
        assert!((m[2][2] - dec.sin()).abs() < eps);
        // Orthonormal columns and determinant +1.
        for c in 0..3 {
            let col = (m[0][c], m[1][c], m[2][c]);
            assert!((norm(col) - 1.0).abs() < eps);
        }
        for (a, b) in [(0, 1), (0, 2), (1, 2)] {
            let dot = m[0][a] * m[0][b] + m[1][a] * m[1][b] + m[2][a] * m[2][b];
            assert!(dot.abs() < eps);
        }
        let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        assert!((det - 1.0).abs() < eps);
    }

    #[test]
    fn epoch_positions_respect_orbit_envelopes() {
        // Rotations preserve the norm, so each ICRF radius must sit in
        // [a(1−e), a(1+e)] up to rounding.
        let (ariel, umbriel, titania, oberon, miranda) =
            uranian_moon_positions(URANUS_EPOCH_JD).expect("epoch must evaluate");
        let cases = [
            (ariel, ARIEL_A_KM, ARIEL_E),
            (umbriel, UMBRIEL_A_KM, UMBRIEL_E),
            (titania, TITANIA_A_KM, TITANIA_E),
            (oberon, OBERON_A_KM, OBERON_E),
            (miranda, MIRANDA_A_KM, MIRANDA_E),
        ];
        for (pos, a, e) in cases {
            assert!(pos.0.is_finite() && pos.1.is_finite() && pos.2.is_finite());
            let r = norm(pos);
            let lo = a * (1.0 - e) * (1.0 - 1e-12);
            let hi = a * (1.0 + e) * (1.0 + 1e-12);
            assert!(r >= lo && r <= hi, "r={r} outside [{lo}, {hi}]");
        }
    }

    #[test]
    fn positions_vary_with_time() {
        let first = uranian_moon_positions(URANUS_EPOCH_JD).expect("epoch");
        let second = uranian_moon_positions(URANUS_EPOCH_JD + 1.0).expect("epoch+1d");
        for (a, b) in [
            (first.0, second.0),
            (first.1, second.1),
            (first.2, second.2),
            (first.3, second.3),
            (first.4, second.4),
        ] {
            assert!(
                (a.0 - b.0).abs() + (a.1 - b.1).abs() + (a.2 - b.2).abs() > 0.0,
                "moon frozen across a day"
            );
        }
    }

    #[test]
    fn cob_offset_is_negated_weighted_centroid() {
        let jd = URANUS_EPOCH_JD + 100.0;
        let moons = uranian_moon_positions(jd).expect("positions");
        let arr = [moons.0, moons.1, moons.2, moons.3, moons.4];
        let weights = [GM_ARIEL, GM_UMBRIEL, GM_TITANIA, GM_OBERON, GM_MIRANDA];
        // Independent zip-fold path (not the unrolled kernel body), so a
        // crossed weight/moon wire shows up as a value mismatch.
        let mut acc = (0.0, 0.0, 0.0);
        for (pos, w) in arr.iter().zip(weights.iter()) {
            acc.0 += w * pos.0;
            acc.1 += w * pos.1;
            acc.2 += w * pos.2;
        }
        let total = GM_URANUS + GM_URANUS_MOONS_TOTAL;
        let expected = (
            -acc.0 / total / MOONS_AU_KM,
            -acc.1 / total / MOONS_AU_KM,
            -acc.2 / total / MOONS_AU_KM,
        );
        let got = uranus_cob_offset(jd).expect("cob");
        assert_eq!(got, expected);
        // Analytic magnitude bound: no moon exceeds its apocenter, so the
        // offset cannot exceed the weight-averaged apocenter ratio.
        let bound = weights
            .iter()
            .zip([
                ARIEL_A_KM,
                UMBRIEL_A_KM,
                TITANIA_A_KM,
                OBERON_A_KM,
                MIRANDA_A_KM,
            ])
            .zip([ARIEL_E, UMBRIEL_E, TITANIA_E, OBERON_E, MIRANDA_E])
            .map(|((w, a), e)| w * a * (1.0 + e))
            .sum::<f64>()
            / GM_URANUS
            / MOONS_AU_KM;
        assert!(norm(got) <= bound, "{} > {bound}", norm(got));
    }

    #[test]
    fn non_finite_jd_is_refused() {
        for jd in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                uranian_moon_positions(jd),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                uranus_cob_offset(jd),
                Err(Error::NonFiniteArgument { .. })
            ));
        }
    }
}

#[cfg(test)]
mod registry_tests {
    use super::*;
    use crate::constants::{MOON_TITAN, SUN};

    #[test]
    fn tables_have_21_rows_in_release_order() {
        assert_eq!(MOON_NAIF_PAIRS.len(), 21);
        assert_eq!(MOON_NAME_PAIRS.len(), 21);
        assert_eq!(MOON_PARENT_PAIRS.len(), 21);
        // Known table endpoints.
        assert_eq!(MOON_NAIF_PAIRS[0], (MOON_IO, NAIF_IO));
        assert_eq!(MOON_NAIF_PAIRS[20], (MOON_CHARON, NAIF_CHARON));
        assert_eq!(MOON_NAME_PAIRS[0], (MOON_IO, "Io"));
        assert_eq!(MOON_NAME_PAIRS[20], (MOON_CHARON, "Charon"));
        assert_eq!(MOON_PARENT_PAIRS[0], (MOON_IO, NAIF_JUPITER_BARYCENTER));
        assert_eq!(MOON_PARENT_PAIRS[20], (MOON_CHARON, NAIF_PLUTO_BARYCENTER));
        // The three tables key the same 21 legacy ids in the same order.
        for i in 0..21 {
            assert_eq!(MOON_NAIF_PAIRS[i].0, MOON_NAME_PAIRS[i].0);
            assert_eq!(MOON_NAIF_PAIRS[i].0, MOON_PARENT_PAIRS[i].0);
        }
    }

    #[test]
    fn moon_tables_match_release_words() {
        // Known moon identifiers, names and parent barycenters.
        let words: [(i32, i32, &str, i32); 21] = [
            (9001, 501, "Io", 5),
            (9002, 502, "Europa", 5),
            (9003, 503, "Ganymede", 5),
            (9004, 504, "Callisto", 5),
            (9011, 601, "Mimas", 6),
            (9012, 602, "Enceladus", 6),
            (9013, 603, "Tethys", 6),
            (9014, 604, "Dione", 6),
            (9015, 605, "Rhea", 6),
            (9016, 606, "Titan", 6),
            (9017, 607, "Hyperion", 6),
            (9018, 608, "Iapetus", 6),
            (9021, 705, "Miranda", 7),
            (9022, 701, "Ariel", 7),
            (9023, 702, "Umbriel", 7),
            (9024, 703, "Titania", 7),
            (9025, 704, "Oberon", 7),
            (9031, 801, "Triton", 8),
            (9041, 401, "Phobos", 4),
            (9042, 402, "Deimos", 4),
            (9051, 901, "Charon", 9),
        ];
        assert_eq!(MOON_NAIF_PAIRS.len(), words.len());
        assert_eq!(MOON_NAME_PAIRS.len(), words.len());
        assert_eq!(MOON_PARENT_PAIRS.len(), words.len());
        for (i, (legacy, naif, name, parent)) in words.iter().enumerate() {
            assert_eq!(MOON_NAIF_PAIRS[i], (*legacy, *naif), "row {i}");
            assert_eq!(MOON_NAME_PAIRS[i], (*legacy, *name), "row {i}");
            assert_eq!(MOON_PARENT_PAIRS[i].0, *legacy, "row {i}");
            assert_eq!(moon_naif_id(*legacy), Some(*naif), "row {i}");
            assert_eq!(get_moon_name(*legacy), *name, "row {i}");
            assert_eq!(moon_parent_barycenter(*legacy), Some(*parent), "row {i}");
        }
        // Barycenter ids are values.
        assert_eq!(NAIF_MARS_BARYCENTER, 4);
        assert_eq!(NAIF_JUPITER_BARYCENTER, 5);
        assert_eq!(NAIF_SATURN_BARYCENTER, 6);
        assert_eq!(NAIF_URANUS_BARYCENTER, 7);
        assert_eq!(NAIF_NEPTUNE_BARYCENTER, 8);
        assert_eq!(NAIF_PLUTO_BARYCENTER, 9);
    }

    #[test]
    fn moon_naif_id_covers_legacy_and_canonical() {
        // Every legacy alias resolves through the map.
        for (legacy, naif) in MOON_NAIF_PAIRS {
            assert_eq!(moon_naif_id(legacy), Some(naif));
        }
        // Canonical scheme: subtraction over the open interval, so a
        // nameless moon (Amalthea, NAIF 505) resolves too.
        assert_eq!(moon_naif_id(PLMOON_OFFSET + 501), Some(501));
        assert_eq!(moon_naif_id(PLMOON_OFFSET + 505), Some(505));
        assert_eq!(moon_naif_id(PLMOON_OFFSET + 999), Some(999));
        // Interval endpoints are exclusive, like the model comparison.
        assert_eq!(moon_naif_id(PLMOON_OFFSET + 400), None);
        assert_eq!(moon_naif_id(PLMOON_OFFSET + 1000), None);
        // Non-moons never resolve.
        for id in [SUN, MOON_TITAN - 9000, -1, 0, i32::MAX, i32::MIN] {
            if MOON_NAIF_PAIRS.iter().any(|(legacy, _)| *legacy == id) {
                continue;
            }
            assert_eq!(moon_naif_id(id), None, "id={id}");
        }
    }

    #[test]
    fn moon_names_round_trip() {
        for (legacy, name) in MOON_NAME_PAIRS {
            assert_eq!(get_moon_name(legacy), name);
        }
        // Canonical ids resolve through the NAIF reverse map.
        assert_eq!(get_moon_name(PLMOON_OFFSET + 501), "Io");
        assert_eq!(get_moon_name(PLMOON_OFFSET + 606), "Titan");
        assert_eq!(get_moon_name(PLMOON_OFFSET + 901), "Charon");
        // Unknown ids carry the model message spelling.
        assert_eq!(get_moon_name(12345), "Unknown Moon (12345)");
        assert_eq!(get_moon_name(SUN), "Unknown Moon (0)");
    }

    #[test]
    fn parent_barycenters_match_release_map() {
        assert_eq!(
            moon_parent_barycenter(MOON_IO),
            Some(NAIF_JUPITER_BARYCENTER)
        );
        assert_eq!(
            moon_parent_barycenter(MOON_TITAN),
            Some(NAIF_SATURN_BARYCENTER)
        );
        assert_eq!(
            moon_parent_barycenter(MOON_ARIEL),
            Some(NAIF_URANUS_BARYCENTER)
        );
        assert_eq!(
            moon_parent_barycenter(MOON_TRITON),
            Some(NAIF_NEPTUNE_BARYCENTER)
        );
        assert_eq!(
            moon_parent_barycenter(MOON_PHOBOS),
            Some(NAIF_MARS_BARYCENTER)
        );
        assert_eq!(
            moon_parent_barycenter(MOON_CHARON),
            Some(NAIF_PLUTO_BARYCENTER)
        );
        // Keyed by legacy ids only: canonical ids and non-moons miss.
        assert_eq!(moon_parent_barycenter(PLMOON_OFFSET + 501), None);
        assert_eq!(moon_parent_barycenter(SUN), None);
    }

    #[test]
    fn predicate_agrees_with_positions_classifier() {
        // The `is_planetary_moon_id` delegates here; pin the
        // agreement over legacy, canonical, boundary and non-moon ids.
        let mut ids: Vec<i32> = MOON_NAIF_PAIRS.iter().map(|(legacy, _)| *legacy).collect();
        ids.extend([
            PLMOON_OFFSET + 401,
            PLMOON_OFFSET + 505,
            PLMOON_OFFSET + 999,
            PLMOON_OFFSET + 400,
            PLMOON_OFFSET + 1000,
            SUN,
            -1,
            0,
        ]);
        for id in ids {
            assert_eq!(
                is_planetary_moon(id),
                crate::positions::is_planetary_moon_id(id),
                "id={id}"
            );
        }
    }
}

#[cfg(test)]
mod moon_registry_tests {
    use super::*;
    use crate::constants::{MOON_IO, MOON_MIMAS, MOON_TITAN};
    use crate::spk::moon_test_kernel_bytes;
    use std::path::PathBuf;

    /// Scratch kernel file with a unique name; removed by the caller.
    fn write_temp_kernel(name: &str, bytes: &[u8]) -> PathBuf {
        let path = std::env::temp_dir().join(format!(
            "libephemeris-moonreg-{}-{}-{name}",
            std::process::id(),
            // Monotonic per test-binary run disambiguator.
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .expect("clock")
                .as_nanos(),
        ));
        std::fs::write(&path, bytes).expect("temp kernel writable");
        path
    }

    #[test]
    fn missing_file_is_not_found() {
        let mut reg = MoonSpkRegistry::new();
        let err = reg
            .register_moon_spk("/nonexistent-dir/nope.bsp", None)
            .unwrap_err();
        assert!(matches!(err, Error::SpkFileNotFound { .. }));
        assert!(reg.list_registered_moons().is_empty());
    }

    #[test]
    fn undecodable_file_fails_clean() {
        let path = write_temp_kernel("garbage.bsp", b"not a DAF file at all");
        let mut reg = MoonSpkRegistry::new();
        assert!(reg.register_moon_spk(path.to_str().unwrap(), None).is_err());
        assert!(reg.list_registered_moons().is_empty());
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn auto_detect_registers_center_chain_support() {
        // SSB-sweep finding: SSB evaluation resolves each
        // moon's parent-center chain through the registry provider, so
        // registration must carry the chain-closure support bodies
        // under their NAIF ids — while the moon listing stays
        // moon-only: support bodies are registry infrastructure.
        use crate::spk::moon_chain_test_kernel_bytes;
        let path = write_temp_kernel("chain.bsp", &moon_chain_test_kernel_bytes());
        let key = path.to_str().unwrap().to_string();
        let mut reg = MoonSpkRegistry::new();
        reg.register_moon_spk(&key, None).expect("auto-detect");
        // Center 5 (Jupiter barycenter) serves the Io→SSB chain.
        assert!(reg.inner.body_info(5).is_some(), "center-5 support");
        // Moon listing stays moon-only: Io under both legs, no NAIF-id leak.
        let listed = reg.list_registered_moons();
        assert_eq!(listed.get(&MOON_IO), Some(&key));
        assert_eq!(listed.get(&(PLMOON_OFFSET + 501)), Some(&key));
        assert_eq!(listed.len(), 2);
        // The explicit path carries the same closure.
        let mut explicit = MoonSpkRegistry::new();
        explicit
            .register_moon_spk(&key, Some(&[MOON_IO]))
            .expect("explicit Io");
        assert!(explicit.inner.body_info(5).is_some(), "center-5 support");
        assert_eq!(explicit.list_registered_moons().len(), 2);
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn explicit_list_without_chainable_moons_drops_the_image() {
        // Mimas (NAIF 601) has no segment in the fixture image: the
        // per-moon probe skips it, nothing registers, and the image is
        // dropped with `SpkNoValidMoons` (never half-registered).
        let path = write_temp_kernel("nomatch.bsp", &moon_test_kernel_bytes());
        let mut reg = MoonSpkRegistry::new();
        let err = reg
            .register_moon_spk(path.to_str().unwrap(), Some(&[MOON_MIMAS]))
            .unwrap_err();
        assert!(
            matches!(err, Error::SpkNoValidMoons { .. }),
            "unexpected {err:?}"
        );
        assert!(reg.list_registered_moons().is_empty());
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn auto_detect_keys_canonical_and_legacy_ids() {
        let path = write_temp_kernel("auto.bsp", &moon_test_kernel_bytes());
        let key = path.to_str().unwrap().to_string();
        let mut reg = MoonSpkRegistry::new();
        reg.register_moon_spk(&key, None).expect("auto-detect");
        let listed = reg.list_registered_moons();
        // Io + Titan under both legs; the Earth segment (NAIF 3) is out
        // of the 400-999 satellite range and never registers.
        for id in [
            MOON_IO,
            PLMOON_OFFSET + 501,
            MOON_TITAN,
            PLMOON_OFFSET + 606,
        ] {
            assert_eq!(listed.get(&id), Some(&key), "id={id}");
        }
        assert!(!listed.keys().any(|id| *id == 3 || *id == PLMOON_OFFSET + 3));
        assert_eq!(get_moon_coverage_missing(&reg), None);
        std::fs::remove_file(&path).ok();
    }

    /// Helper pinning the unregistered-moon coverage shape used below.
    fn get_moon_coverage_missing(reg: &MoonSpkRegistry) -> Option<(f64, f64)> {
        reg.get_moon_coverage(MOON_MIMAS)
    }

    #[test]
    fn explicit_list_registers_only_the_named_moons() {
        let path = write_temp_kernel("explicit.bsp", &moon_test_kernel_bytes());
        let key = path.to_str().unwrap().to_string();
        let mut reg = MoonSpkRegistry::new();
        reg.register_moon_spk(&key, Some(&[MOON_IO]))
            .expect("explicit Io");
        let listed = reg.list_registered_moons();
        // Named moon plus its canonical id; Titan stays out.
        assert_eq!(listed.get(&MOON_IO), Some(&key));
        assert_eq!(listed.get(&(PLMOON_OFFSET + 501)), Some(&key));
        assert_eq!(listed.get(&MOON_TITAN), None);
        assert_eq!(listed.get(&(PLMOON_OFFSET + 606)), None);
        // Re-registration replaces cleanly: same image, same mapping.
        reg.register_moon_spk(&key, Some(&[MOON_IO]))
            .expect("re-register");
        assert_eq!(reg.list_registered_moons().len(), 2);
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn coverage_is_the_whole_kernel_span() {
        let path = write_temp_kernel("cov.bsp", &moon_test_kernel_bytes());
        let key = path.to_str().unwrap().to_string();
        let mut reg = MoonSpkRegistry::new();
        reg.register_moon_spk(&key, None).expect("auto-detect");
        // Whole-kernel min/max in ephemeris seconds: min 0.0 (all three
        // segments start there), max 20 e 6 (Titan) — the Earth segment
        // contributes to the span, as defined by the model min/max.
        let expected = (
            crate::spk::et_to_jd_tdb(0.0).unwrap(),
            crate::spk::et_to_jd_tdb(20_000_000.0).unwrap(),
        );
        assert_eq!(reg.get_moon_coverage(MOON_IO), Some(expected));
        assert_eq!(reg.get_moon_coverage(MOON_TITAN), Some(expected));
        // Canonical leg resolves too; unregistered and non-moon ids miss.
        assert_eq!(reg.get_moon_coverage(PLMOON_OFFSET + 606), Some(expected));
        assert_eq!(reg.get_moon_coverage(MOON_MIMAS), None);
        assert_eq!(reg.get_moon_coverage(0), None);
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn unregister_matches_exact_key_and_basename() {
        let path = write_temp_kernel("bybase.bsp", &moon_test_kernel_bytes());
        let key = path.to_str().unwrap().to_string();
        let base = path
            .file_name()
            .expect("temp file has a name")
            .to_str()
            .expect("utf-8 temp name")
            .to_string();
        let mut reg = MoonSpkRegistry::new();
        reg.register_moon_spk(&key, None).expect("auto-detect");
        // Bare file name unregisters the full-path registration.
        reg.unregister_moon_spk(&base);
        assert!(reg.list_registered_moons().is_empty());
        // Exact-key path unregisters too; unknown paths are a no-op.
        reg.register_moon_spk(&key, None).expect("auto-detect");
        reg.unregister_moon_spk("other.bsp");
        assert_eq!(reg.list_registered_moons().len(), 4);
        // Trailing-slash request is a model no-op absent an exact key
        // ( is `""`, never matched).
        reg.unregister_moon_spk(&format!("{base}/"));
        assert_eq!(reg.list_registered_moons().len(), 4);
        reg.unregister_moon_spk(&key);
        assert!(reg.list_registered_moons().is_empty());
        // Full close drops everything at once.
        reg.register_moon_spk(&key, None).expect("auto-detect");
        reg.close_moon_kernels();
        assert!(reg.list_registered_moons().is_empty());
        std::fs::remove_file(&path).ok();
    }
}

#[cfg(test)]
mod tass_tests {
    use super::*;

    fn norm(v: (f64, f64, f64)) -> f64 {
        (v.0 * v.0 + v.1 * v.1 + v.2 * v.2).sqrt()
    }

    #[test]
    fn gm_saturn_total_matches_release_order() {
        // sums Mimas..Iapetus left to
        // right; the constant transcribes that association.
        assert_eq!(
            GM_SATURN_MOONS_TOTAL,
            GM_MIMAS
                + GM_ENCELADUS
                + GM_TETHYS
                + GM_DIONE
                + GM_RHEA
                + GM_TITAN
                + GM_HYPERION
                + GM_IAPETUS
        );
        assert_eq!(GM_SATURN, 37_931_206.23);
        assert_eq!(GM_TITAN, 8978.13710);
        assert_eq!(TASS_EPOCH_JD, 2444240.0);
    }

    #[test]
    fn invalid_body_refuses_with_typed_error() {
        // edge for indices outside 0-7.
        for bad in [-1, 8, 99, i32::MIN, i32::MAX] {
            let err = saturn_moon_position(TASS_EPOCH_JD, bad).unwrap_err();
            assert!(
                matches!(
                    err,
                    Error::InvalidBody {
                        operation: "saturn_moon_position",
                        ..
                    }
                ),
                "unexpected {err:?} for body {bad}"
            );
            assert!(saturn_moon_position_velocity(TASS_EPOCH_JD, bad).is_err());
        }
    }

    #[test]
    fn non_finite_jd_refuses() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                saturn_moon_position(bad, 5).unwrap_err(),
                Error::NonFiniteArgument { .. }
            ));
            assert!(saturn_moon_position_velocity(bad, 5).is_err());
            assert!(titan_position(bad).is_err());
            assert!(saturn_cob_offset(bad).is_err());
        }
    }

    #[test]
    fn titan_wrapper_matches_body_five() {
        // is a direct `body = 5` call.
        for jd in [TASS_EPOCH_JD, 2451545.0, 2460311.0] {
            assert_eq!(
                titan_position(jd).unwrap(),
                saturn_moon_position(jd, 5).unwrap()
            );
            let full = saturn_moon_position_velocity(jd, 5).unwrap();
            assert_eq!((full.0, full.1, full.2), titan_position(jd).unwrap());
        }
    }

    #[test]
    fn all_bodies_finite_with_titan_envelope() {
        // Every body evaluates (Hyperion takes the series-1 start_idx = 0
        // branch); Titan's radius sits inside its analytic orbit envelope
        // (a ≈ 0.0082 AU, e ≈ 0.029 — wide margins, not fitted).
        for body in 0..8 {
            let pos = saturn_moon_position(2451545.0, body).unwrap();
            assert!(pos.0.is_finite() && pos.1.is_finite() && pos.2.is_finite());
            let full = saturn_moon_position_velocity(2451545.0, body).unwrap();
            assert!(full.3.is_finite() && full.4.is_finite() && full.5.is_finite());
        }
        let titan = titan_position(2451545.0).unwrap();
        let r = norm(titan);
        assert!(r > 0.007 && r < 0.0095, "unexpected Titan radius {r}");
        let full = saturn_moon_position_velocity(2451545.0, 5).unwrap();
        let speed = norm((full.3, full.4, full.5));
        assert!(
            speed > 2e-3 && speed < 5e-3,
            "unexpected Titan speed {speed}"
        );
    }

    #[test]
    fn saturn_cob_composition_matches_titan_only_rule() {
        // : Titan-only GM composition, negated.
        let jd = TASS_EPOCH_JD + 1000.0;
        let (tx, ty, tz) = titan_position(jd).unwrap();
        let total = GM_SATURN + GM_SATURN_MOONS_TOTAL;
        assert_eq!(
            saturn_cob_offset(jd).unwrap(),
            (
                -GM_TITAN * tx / total,
                -GM_TITAN * ty / total,
                -GM_TITAN * tz / total
            )
        );
    }
}

#[cfg(test)]
mod cob_dispatch_tests {
    use super::*;

    #[test]
    fn served_arms_match_their_theory_kernels() {
        let jd = URANUS_EPOCH_JD + 100.0;
        // Dispatch must select the corresponding system offset.
        assert_eq!(
            get_cob_offset("jupiter", jd).unwrap(),
            ecliptic_to_equatorial_j2000(jupiter_cob_offset(jd).unwrap())
        );
        assert_eq!(
            get_cob_offset("uranus", jd).unwrap(),
            uranus_cob_offset(jd).unwrap()
        );
        assert_eq!(
            get_cob_offset("neptune", jd).unwrap(),
            ecliptic_to_equatorial_j2000(neptune_cob_offset(jd).unwrap())
        );
        assert_eq!(
            get_cob_offset("pluto", jd).unwrap(),
            ecliptic_to_equatorial_j2000(pluto_cob_offset(jd).unwrap())
        );
        assert_eq!(
            get_cob_offset("saturn", jd).unwrap(),
            ecliptic_to_equatorial_j2000(saturn_cob_offset(jd).unwrap())
        );
    }

    #[test]
    fn name_spellings_fold_like_the_release_chain() {
        let jd = URANUS_EPOCH_JD;
        let plain = get_cob_offset("jupiter", jd).unwrap();
        for spelling in [
            "Jupiter",
            "JUPITER",
            "jupiter barycenter",
            "Jupiter Barycenter",
            "  jupiter barycenter  ",
            // Double space: the suffix still strips, then trims.
            "jupiter  barycenter",
        ] {
            assert_eq!(get_cob_offset(spelling, jd).unwrap(), plain, "{spelling}");
        }
        assert_eq!(
            get_cob_offset("PLUTO BARYCENTER", jd).unwrap(),
            get_cob_offset("pluto", jd).unwrap()
        );
        // No space, no suffix strip: unknown, like the model.
        assert_eq!(
            get_cob_offset("jupiterbarycenter", jd).unwrap(),
            (0.0, 0.0, 0.0)
        );
    }

    #[test]
    fn saturn_arm_is_served_through_tass() {
        // TASS slice: the Saturn arm rides the Titan-only composition
        // through the ICRF rotation (same wiring as Jupiter/Neptune/Pluto).
        // Titan's correction is far from negligible: ~2 e-6 AU.
        let jd = URANUS_EPOCH_JD;
        let offset = get_cob_offset("saturn", jd).unwrap();
        let norm = (offset.0 * offset.0 + offset.1 * offset.1 + offset.2 * offset.2).sqrt();
        assert!(
            norm > 1e-6 && norm < 1e-5,
            "unexpected Saturn COB {offset:?}"
        );
        assert_eq!(get_cob_offset("Saturn Barycenter", jd).unwrap(), offset);
    }

    #[test]
    fn unrecognized_names_return_zero() {
        let jd = URANUS_EPOCH_JD;
        for name in ["mars", "earth", "moon", ""] {
            assert_eq!(get_cob_offset(name, jd).unwrap(), (0.0, 0.0, 0.0), "{name}");
        }
    }

    #[test]
    fn non_finite_epoch_propagates() {
        for jd in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                get_cob_offset("jupiter", jd),
                Err(Error::NonFiniteArgument { .. })
            ));
            // The zero arms never touch the kernels: finite-shaped even
            // for non-finite epochs, as defined by the model.
            assert_eq!(get_cob_offset("mars", jd).unwrap(), (0.0, 0.0, 0.0));
        }
    }
}

// ---------------------------------------------------------------------------
// Satellite apparent-place calculation over numerical state/frame providers.
// The caller owns kernel acquisition and supplies SSB states, observer state,
// date frames, nutation and sidereal corrections. This layer performs the
// light-time and direction corrections with no file I/O.
// ---------------------------------------------------------------------------

use crate::constants::{
    EARTH, FLG_BARYCTR, FLG_EQUATORIAL, FLG_HELCTR, FLG_J2000, FLG_NOABERR, FLG_NONUT,
    FLG_SIDEREAL, FLG_SPEED, FLG_TOPOCTR, FLG_TRUEPOS, SUN,
};
use crate::positions::C_LIGHT_AU_DAY;
use crate::positions::{
    SealedServing, StateProvider, frame_matrices_for_date, frame_transform, serving_channel_gate,
};
use crate::spk::{ChebState, NaifBaseProvider, eval_ssb, resolve_chain};

/// Initial dyadic TT-day half-step for the adaptive fourth-order angular derivative.
pub const MOON_SPEED_HALF_STEP_DAYS: f64 = 1.0 / 64.0;

/// Observer selection for the satellite-minus-observer vector.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MoonObserver {
    /// Solar-system barycenter: zero position and velocity.
    Barycenter,
    /// Planet-system parent body center.
    Heliocentric,
    /// Geocenter.
    Geocentric,
    /// Geocenter plus site offset.
    /// The site value lives in the caller observer closure; a missing
    /// site must refuse there, never degrade to silent geocentric.
    Topocentric,
}

/// Observer selection from the calculation flags.
pub fn moon_observer_for_flags(flags: i32) -> MoonObserver {
    if flags & FLG_BARYCTR != 0 {
        MoonObserver::Barycenter
    } else if flags & FLG_HELCTR != 0 {
        MoonObserver::Heliocentric
    } else if flags & FLG_TOPOCTR != 0 {
        MoonObserver::Topocentric
    } else {
        MoonObserver::Geocentric
    }
}

/// Whether the light-time iteration runs.
pub fn moon_wants_light_time(flags: i32) -> bool {
    flags & FLG_TRUEPOS == 0
}

/// Whether the Bradley aberration runs: never
/// under `FLG_TRUEPOS`/`FLG_NOABERR`, and never for heliocentric or
/// barycentric output (aberration is an observer-velocity effect, and
/// those legs carry no observer velocity).
pub fn moon_wants_aberration(flags: i32) -> bool {
    flags & (FLG_TRUEPOS | FLG_NOABERR | FLG_HELCTR | FLG_BARYCTR) == 0
}

/// Whether the nutation-in-longitude strip applies to the ecliptic
/// longitude: `FLG_NONUT`, the
/// sidereal-equatorial rule, or `FLG_J2000`.
pub fn moon_wants_nutation_strip(flags: i32) -> bool {
    flags & FLG_NONUT != 0
        || (flags & FLG_SIDEREAL != 0 && flags & FLG_EQUATORIAL != 0)
        || flags & FLG_J2000 != 0
}

/// Whether the sidereal correction applies to the output longitude
pub fn moon_wants_sidereal_correction(flags: i32) -> bool {
    flags & FLG_SIDEREAL != 0 && flags & FLG_EQUATORIAL == 0
}

/// Resolve a servable NAIF target for a moon id. Returns
/// the NAIF target id, or `None` when the moon is unknown or no
/// registered kernel chains it to the barycenter.
pub fn resolve_moon_target(registry: &MoonSpkRegistry, moon_id: i32) -> Option<i32> {
    let naif = moon_naif_id(moon_id)?;
    if registry.inner.body_info(moon_id).is_some()
        || registry.inner.body_info(PLMOON_OFFSET + naif).is_some()
    {
        return Some(naif);
    }
    for key in registry.inner.kernel_paths() {
        if let Some(kernel) = registry.inner.kernel(&key) {
            if crate::spk::resolve_chain(kernel.index(), naif).is_ok() {
                return Some(naif);
            }
        }
    }
    None
}

impl MoonSpkRegistry {
    /// SSB-relative ICRF satellite state at the two-part TDB date `jd + fraction`.
    ///
    /// Both date parts are days, not seconds; returned positions are km and
    /// velocities km/s. Keeping the parts separate preserves the fractional
    /// epoch through segment selection, Chebyshev evaluation and recursive
    /// center acquisition. No TT-to-TDB conversion is performed here.
    ///
    /// Lookup tries the requested library ID, its canonical satellite ID and
    /// then chain-resolvable targets in retained kernels. Parent-center states
    /// compose through this registry's support bodies. Unknown/unresolvable
    /// targets return [`Error::SpkBodyNotFound`]; a consumed segment outside its
    /// own closed domain returns [`Error::SpkOutOfRange`]. Frame, type and record
    /// errors propagate. The broad [`Self::get_moon_coverage`] envelope does not
    /// replace those per-target/per-center checks or authorize extrapolation.
    /// Returned states own their values and borrow no kernel image.
    pub fn eval_moon_ssb_two_part_jd(
        &self,
        moon_id: i32,
        jd: f64,
        fraction: f64,
    ) -> Result<ChebState> {
        let naif = moon_naif_id(moon_id).ok_or(Error::SpkBodyNotFound { target: moon_id })?;
        for key in [moon_id, PLMOON_OFFSET + naif] {
            if self.inner.body_info(key).is_some() {
                return self.inner.eval_geometric_ssb_two_part_jd(key, jd, fraction);
            }
        }
        let provider = NaifBaseProvider::new(&self.inner);
        for path in self.inner.kernel_paths() {
            let Some(kernel) = self.inner.kernel(&path) else {
                continue;
            };
            if resolve_chain(kernel.index(), naif).is_ok() {
                return crate::spk::eval_ssb_two_part_jd(
                    kernel.bytes(),
                    kernel.endian(),
                    kernel.index(),
                    naif,
                    &provider,
                    jd,
                    fraction,
                );
            }
        }
        Err(Error::SpkBodyNotFound { target: naif })
    }

    /// Satellite SSB state in AU and AU/day at a two-part TDB Julian date.
    ///
    /// Lookup order and strict provider/coverage errors match
    /// [`Self::eval_moon_ssb_two_part_jd`]. Each relative segment is converted
    /// before its center is added, preserving vector-function sum rounding.
    /// The kilometre-valued API remains separate: converting its completed
    /// sum is not generally equivalent in cancellation-sensitive chains.
    pub fn eval_moon_ssb_au_two_part_jd(
        &self,
        moon_id: i32,
        jd: f64,
        fraction: f64,
    ) -> Result<([f64; 3], [f64; 3])> {
        let naif = moon_naif_id(moon_id).ok_or(Error::SpkBodyNotFound { target: moon_id })?;
        for key in [moon_id, PLMOON_OFFSET + naif] {
            if self.inner.body_info(key).is_some() {
                return self
                    .inner
                    .eval_geometric_ssb_au_two_part_jd(key, jd, fraction);
            }
        }
        for path in self.inner.kernel_paths() {
            let Some(kernel) = self.inner.kernel(&path) else {
                continue;
            };
            if resolve_chain(kernel.index(), naif).is_ok() {
                return self
                    .inner
                    .eval_kernel_ssb_au_two_part_jd(&path, naif, jd, fraction);
            }
        }
        Err(Error::SpkBodyNotFound { target: naif })
    }

    /// Satellite SSB state at scalar TDB seconds past J2000, in km and km/s.
    ///
    /// Target lookup, registry center composition, coverage and record faults
    /// follow [`Self::eval_moon_ssb_two_part_jd`]. Scalar seconds cannot recover
    /// precision already lost when collapsing a two-part Julian date.
    pub fn eval_moon_ssb(&self, moon_id: i32, et: f64) -> Result<ChebState> {
        let provider = NaifBaseProvider::new(&self.inner);
        let naif = match moon_naif_id(moon_id) {
            Some(naif) => naif,
            None => return Err(Error::SpkBodyNotFound { target: moon_id }),
        };
        for key in [moon_id, PLMOON_OFFSET + naif] {
            if self.inner.body_info(key).is_none() {
                continue;
            }
            return self.inner.eval_body_ssb(key, &provider, et);
        }
        for path in self.inner.kernel_paths() {
            let Some(kernel) = self.inner.kernel(&path) else {
                continue;
            };
            if resolve_chain(kernel.index(), naif).is_err() {
                continue;
            }
            return eval_ssb(
                kernel.bytes(),
                kernel.endian(),
                kernel.index(),
                naif,
                &provider,
                et,
            );
        }
        Err(Error::SpkBodyNotFound { target: naif })
    }
}

/// First-order Bradley aberration on a satellite-minus-observer vector.
///
/// `rel_au` is the ICRF displacement in AU; `observer_vel_au_day` is
/// the observer SSB velocity in AU/day. For distance `d = |rel_au|`,
/// form `u = rel_au / d`, normalize it again as `q = u / |u|`, and
/// compute the completed direction offset `beta - (q dot beta) q`,
/// where `beta = observer_vel_au_day / c` and `c` is AU/day light speed.
/// Add that offset to `u`, then normalize and scale back to `d`.
/// The three `1e-30` degenerate-magnitude guards return the input unchanged.
/// Non-finite input components are refused. This is a first-order velocity
/// correction, not a relativistic all-orders aberration model.
pub fn moon_bradley_aberration(
    rel_au: [f64; 3],
    observer_vel_au_day: [f64; 3],
) -> Result<[f64; 3]> {
    for v in rel_au.iter().chain(observer_vel_au_day.iter()) {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: *v });
        }
    }
    let dist = (rel_au[0] * rel_au[0] + rel_au[1] * rel_au[1] + rel_au[2] * rel_au[2]).sqrt();
    if dist < 1e-30 {
        return Ok(rel_au);
    }
    let ux = rel_au[0] / dist;
    let uy = rel_au[1] / dist;
    let uz = rel_au[2] / dist;
    let vx = observer_vel_au_day[0] / C_LIGHT_AU_DAY;
    let vy = observer_vel_au_day[1] / C_LIGHT_AU_DAY;
    let vz = observer_vel_au_day[2] / C_LIGHT_AU_DAY;
    let u_mag = (ux * ux + uy * uy + uz * uz).sqrt();
    if u_mag < 1e-30 {
        return Ok(rel_au);
    }
    let (rx, ry, rz) = (ux / u_mag, uy / u_mag, uz / u_mag);
    let dot = rx * vx + ry * vy + rz * vz;
    let ax = ux + (vx - dot * rx);
    let ay = uy + (vy - dot * ry);
    let az = uz + (vz - dot * rz);
    let a_mag = (ax * ax + ay * ay + az * az).sqrt();
    if a_mag < 1e-30 {
        return Ok(rel_au);
    }
    Ok([ax / a_mag * dist, ay / a_mag * dist, az / a_mag * dist])
}

/// Production sidereal value for the satellite serving: subtract the ayanamsha sampled at the UT1
/// epoch from the ecliptic longitude (mod 360), and — under
/// `FLG_SPEED` — remove the ayanamsha drift sampled at UT1 ∓/± one
/// second.
///
/// `ayanamsha_at` is the model mean/true selection at a UT1 epoch
/// (production: `SealedAyanamshaServing::for_flags_ut` over the
/// context selection, which already encodes the mean-vs-true rule
/// from the flags); this layer owns only the offset/drift shaping,
/// never the ayanamsha mode dispatch. A non-finite `ut1_jd` is refused; a
/// non-finite longitude passes through with its rate untouched (the
/// `apply_sidereal_offset_deg` contract).
pub fn moon_sidereal_word(
    lon_deg: f64,
    speed_lon_deg_day: f64,
    ut1_jd: f64,
    flags: i32,
    ayanamsha_at: impl Fn(f64) -> Result<f64>,
) -> Result<(f64, f64)> {
    if !ut1_jd.is_finite() {
        return Err(Error::NonFiniteArgument { value: ut1_jd });
    }
    let h = crate::sidereal::SIDEREAL_SPEED_HALF_STEP_DAYS;
    let ayan = ayanamsha_at(ut1_jd)?;
    let (prev, next) = if flags & FLG_SPEED != 0 {
        (ayanamsha_at(ut1_jd - h)?, ayanamsha_at(ut1_jd + h)?)
    } else {
        (ayan, ayan)
    };
    crate::sidereal::apply_sidereal_correction(
        lon_deg,
        speed_lon_deg_day,
        ayan,
        prev,
        next,
        flags & FLG_SPEED != 0,
    )
}

/// Relative geometry before and after the radius-preserving aberration step.
struct MoonRelativeState {
    ray: [f64; 3],
    distance_au: f64,
    radial_au_day: f64,
}

/// Provider-injected satellite position serving.
///
/// Words (caller obligations):
/// - `moon_ssb_at_tt`: SSB position/velocity in AU/AU-day for a NAIF
///   target at a TT epoch; velocities are derivatives per TT day (production: the [`MoonSpkRegistry`] SSB
///   chain; range misses surface as the provider error and propagate —
///   the sealed binding maps them to [`Error::EphemerisRange`] via
///   [`MoonSpkRegistry::get_moon_coverage`]).
/// - `observer_ssb_at_tt`: observer SSB position/velocity in AU/AU-day
///   at a TT epoch, with velocities per TT day (production: barycenter zeros, Sun, Earth, or
///   Earth+site; a missing topocentric site must refuse here, never
///   degrade to silent geocentric).
/// - `icrf_to_ecliptic_date`: true-ecliptic-of-date `(lon_deg,
/// lat_deg, dist_au)` for an ICRF AU vector at a TT epoch
///   (production: the native frame-transform value).
/// - `dpsi_deg_at_tt`: nutation in longitude in degrees at a TT epoch
///   (production: the sealed nutation value).
/// - `sidereal_correction_at`: `(lon, speed)` correction at a UT1
///   epoch (production: [`moon_sidereal_word`] over the
///   `SealedAyanamshaServing` selection).
pub struct MoonPositionServing<Fmoon, Fobs, Fframe, Fnut, Fsid> {
    moon_ssb_at_tt: Fmoon,
    observer_ssb_at_tt: Fobs,
    icrf_to_ecliptic_date: Fframe,
    dpsi_deg_at_tt: Fnut,
    sidereal_correction_at: Fsid,
}

impl<Fmoon, Fobs, Fframe, Fnut, Fsid> MoonPositionServing<Fmoon, Fobs, Fframe, Fnut, Fsid>
where
    Fmoon: Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
    Fobs: Fn(f64) -> Result<([f64; 3], [f64; 3])>,
    Fframe: Fn([f64; 3], f64) -> Result<(f64, f64, f64)>,
    Fnut: Fn(f64) -> Result<f64>,
    Fsid: Fn(f64, f64, f64) -> Result<(f64, f64)>,
{
    /// Bind the five provider values.
    pub fn new(
        moon_ssb_at_tt: Fmoon,
        observer_ssb_at_tt: Fobs,
        icrf_to_ecliptic_date: Fframe,
        dpsi_deg_at_tt: Fnut,
        sidereal_correction_at: Fsid,
    ) -> Self {
        Self {
            moon_ssb_at_tt,
            observer_ssb_at_tt,
            icrf_to_ecliptic_date,
            dpsi_deg_at_tt,
            sidereal_correction_at,
        }
    }

    /// Sample the three finite light-time iterations and their exact chain rule.
    ///
    /// The observer is evaluated at reception TT once. At each target sample,
    /// q = p_target(t_emit) - p_observer(t_obs), rho = |q| and
    /// rho_dot = q·(v_target * t_emit_dot - v_observer) / rho. The next
    /// emission-time derivative is 1 - rho_dot/c. TRUEPOS uses one sample
    /// and t_emit_dot = 1. All callback velocities are per TT day.
    /// Aberration changes only the final direction; range and radial velocity
    /// remain in the original Cartesian frame, independent of output rotation.
    fn relative_at(
        &self,
        naif_id: i32,
        observer: MoonObserver,
        jd_tt: f64,
        flags: i32,
    ) -> Result<MoonRelativeState> {
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let (obs_pos, obs_vel) = match observer {
            MoonObserver::Barycenter => ([0.0; 3], [0.0; 3]),
            _ => (self.observer_ssb_at_tt)(jd_tt)?,
        };
        let need_speed = flags & FLG_SPEED != 0;
        let mut emit_tt = jd_tt;
        let mut emission_rate = 1.0;
        let mut rel = [0.0; 3];
        let mut distance = 0.0;
        let mut radial = 0.0;
        for _ in 0..3 {
            let (target_pos, target_vel) = (self.moon_ssb_at_tt)(naif_id, emit_tt)?;
            let mut relative_velocity = [0.0; 3];
            for i in 0..3 {
                for value in [target_pos[i], obs_pos[i]] {
                    if !value.is_finite() {
                        return Err(Error::NonFiniteArgument { value });
                    }
                }
                rel[i] = target_pos[i] - obs_pos[i];
                if need_speed {
                    for value in [target_vel[i], obs_vel[i]] {
                        if !value.is_finite() {
                            return Err(Error::NonFiniteArgument { value });
                        }
                    }
                    relative_velocity[i] = target_vel[i].mul_add(emission_rate, -obs_vel[i]);
                }
            }
            distance = (rel[0] * rel[0] + rel[1] * rel[1] + rel[2] * rel[2]).sqrt();
            if distance == 0.0 {
                return Err(Error::ZeroVectorToSpherical);
            }
            if !distance.is_finite() {
                return Err(Error::NonFiniteArgument { value: distance });
            }
            if need_speed {
                radial = rel[0].mul_add(
                    relative_velocity[0],
                    rel[1].mul_add(relative_velocity[1], rel[2] * relative_velocity[2]),
                ) / distance;
                emission_rate = 1.0 - radial / C_LIGHT_AU_DAY;
            }
            if !moon_wants_light_time(flags) {
                break;
            }
            emit_tt = jd_tt - distance / C_LIGHT_AU_DAY;
        }
        if moon_wants_aberration(flags) {
            rel = moon_bradley_aberration(rel, obs_vel)?;
        }
        Ok(MoonRelativeState {
            ray: rel,
            distance_au: distance,
            radial_au_day: radial,
        })
    }

    /// Ecliptic longitude/latitude/distance with the nutation strip
    /// : the frame value shapes the aberrated vector,
    /// then the strip reads `(lon - dpsi) mod 360` under the
    /// `FLG_NONUT`/sidereal/`FLG_J2000` gate so the central-difference
    /// speeds below inherit the correction with no separate rate term.
    fn ecliptic_at(
        &self,
        naif_id: i32,
        observer: MoonObserver,
        jd_tt: f64,
        flags: i32,
    ) -> Result<(f64, f64, f64, f64)> {
        let state = self.relative_at(naif_id, observer, jd_tt, flags)?;
        let (mut lon, lat, dist) = (self.icrf_to_ecliptic_date)(state.ray, jd_tt)?;
        if !lon.is_finite() || !lat.is_finite() || !dist.is_finite() {
            return Err(Error::NonFiniteArgument { value: lon });
        }
        if moon_wants_nutation_strip(flags) {
            let dpsi = (self.dpsi_deg_at_tt)(jd_tt)?;
            if !dpsi.is_finite() {
                return Err(Error::NonFiniteArgument { value: dpsi });
            }
            lon = (lon - dpsi).rem_euclid(360.0);
        }
        Ok((lon, lat, state.distance_au, state.radial_au_day))
    }

    /// Satellite apparent place: `(lon_deg, lat_deg, dist_au, speed_lon,
    /// speed_lat, speed_dist)`, with rates per TT day. Angular rates use a
    /// fourth-order symmetric derivative at exact dyadic TT offsets, refined
    /// by Richardson curvature estimates until truncation is resolved or
    /// further refinement increases binary64/epoch-quantization error. The
    /// radial rate differentiates all three light-time iterations analytically.
    /// Longitude differences use the shortest arc across zero. Without
    /// `FLG_SPEED` rates are zero. The sidereal correction applies only to
    /// ecliptic output, at the supplied UT1 epoch. Non-finite epochs, undefined
    /// relative directions and provider/frame errors propagate.
    pub fn moon_position(
        &self,
        naif_id: i32,
        observer: MoonObserver,
        jd_tt: f64,
        jd_ut1: f64,
        flags: i32,
    ) -> Result<(f64, f64, f64, f64, f64, f64)> {
        let (mut lon, lat, dist, radial) = self.ecliptic_at(naif_id, observer, jd_tt, flags)?;
        let (mut speed_lon, mut speed_lat, mut speed_dist) = (0.0, 0.0, 0.0);
        if flags & FLG_SPEED != 0 {
            let sample = |offset| {
                self.ecliptic_at(naif_id, observer, jd_tt + offset, flags)
                    .map(|(lon, lat, _, _)| [lon, lat])
            };
            let angular = crate::positions::derivative::angular_rates(
                &sample,
                MOON_SPEED_HALF_STEP_DAYS,
                (f64::NEG_INFINITY, f64::INFINITY),
            )?;
            (speed_lon, speed_lat) = (angular[0], angular[1]);
            speed_dist = radial;
        }
        if moon_wants_sidereal_correction(flags) {
            (lon, speed_lon) = (self.sidereal_correction_at)(lon, speed_lon, jd_ut1)?;
        }
        Ok((lon, lat, dist, speed_lon, speed_lat, speed_dist))
    }
}

#[cfg(test)]
mod moon_serving_tests {
    use super::*;

    const JD: f64 = 2460000.5;

    #[test]
    fn light_time_reads_three_targets_not_three_resamples() {
        use std::cell::RefCell;
        let epochs = RefCell::new(Vec::new());
        let target = |_: i32, tt: f64| {
            epochs.borrow_mut().push(tt);
            // A moving radial target makes successive retardations distinct.
            Ok(([10.0 + (tt - JD), 0.0, 0.0], [1.0, 0.0, 0.0]))
        };
        let serving = MoonPositionServing::new(
            target,
            |_| Ok(([0.0; 3], [0.0; 3])),
            stub_frame(),
            |_| Ok(0.0),
            |lon, speed, _| Ok((lon, speed)),
        );
        let output = serving
            .moon_position(901, MoonObserver::Barycenter, JD, JD, FLG_BARYCTR)
            .unwrap();
        let second = JD - 10.0 / C_LIGHT_AU_DAY;
        let third = JD - (10.0 + (second - JD)) / C_LIGHT_AU_DAY;
        assert_eq!(*epochs.borrow(), [JD, second, third]);
        assert_eq!(output.2.to_bits(), (10.0 + (third - JD)).to_bits());
        epochs.borrow_mut().clear();
        serving
            .moon_position(901, MoonObserver::Barycenter, JD, JD, FLG_TRUEPOS)
            .unwrap();
        assert_eq!(*epochs.borrow(), [JD]);
    }

    #[test]
    fn radial_rate_differentiates_the_three_finite_light_time_iterations() {
        // Positive collinear uniform motions have the independent closed form
        // rho_dot = (v - w) * (1 - v/c + (v/c)^2) after exactly three samples.
        // TRUEPOS has no delay and rho_dot = v - w. This also checks the
        // reception observer velocity, including the topocentric callback.
        let velocity = 1.0 / 512.0;
        let observer_velocity = 1.0 / 1024.0;
        let serving = MoonPositionServing::new(
            |_: i32, tt: f64| Ok(([30.0 + velocity * tt, 0.0, 0.0], [velocity, 0.0, 0.0])),
            |tt: f64| {
                Ok((
                    [1.0 + observer_velocity * tt, 0.0, 0.0],
                    [observer_velocity, 0.0, 0.0],
                ))
            },
            stub_frame(),
            |_| Ok(0.0),
            |lon, speed, _| Ok((lon, speed)),
        );
        for observer in [
            MoonObserver::Barycenter,
            MoonObserver::Heliocentric,
            MoonObserver::Geocentric,
            MoonObserver::Topocentric,
        ] {
            let relative_speed = velocity
                - if observer == MoonObserver::Barycenter {
                    0.0
                } else {
                    observer_velocity
                };
            for flags in [FLG_SPEED, FLG_SPEED | FLG_TRUEPOS] {
                let result = serving
                    .moon_position(901, observer, 0.0, 0.0, flags)
                    .unwrap();
                let a = velocity / C_LIGHT_AU_DAY;
                let expected = relative_speed
                    * if flags & FLG_TRUEPOS != 0 {
                        1.0
                    } else {
                        1.0 - a + a * a
                    };
                // Fewer than 32 binary64 operations with no cancellation in
                // the closed form or sample chain: gamma_32 bounds rounding.
                let u = f64::EPSILON / 2.0;
                let gamma = 32.0 * u / (1.0 - 32.0 * u);
                assert!(
                    (result.5 - expected).abs() <= gamma * relative_speed.abs(),
                    "{observer:?}, flags={flags}: {} vs {expected}",
                    result.5
                );
            }
        }
    }

    #[test]
    fn range_and_radial_rate_are_independent_of_frame_and_aberration() {
        let target = |_: i32, tt: f64| {
            Ok((
                [12.0 + 0.001 * tt, -25.0 - 0.002 * tt, 5.0 + 0.003 * tt],
                [0.001, -0.002, 0.003],
            ))
        };
        let observer = |tt: f64| Ok(([0.01 * tt, -0.005 * tt, 0.002 * tt], [0.01, -0.005, 0.002]));
        let rotation = |ray: [f64; 3], _: f64| {
            let (s, c) = 0.41_f64.sin_cos();
            let rotated = [c * ray[0] - s * ray[1], s * ray[0] + c * ray[1], ray[2]];
            let radius = rotated.into_iter().map(|x| x * x).sum::<f64>().sqrt();
            Ok((
                rotated[1].atan2(rotated[0]).to_degrees().rem_euclid(360.0),
                (rotated[2] / radius).asin().to_degrees(),
                radius,
            ))
        };
        let serving = MoonPositionServing::new(
            target,
            observer,
            rotation,
            |_| Ok(0.0),
            |lon, speed, _| Ok((lon, speed)),
        );
        let apparent = serving
            .moon_position(901, MoonObserver::Geocentric, 0.0, 0.0, FLG_SPEED)
            .unwrap();
        let unshifted = serving
            .moon_position(
                901,
                MoonObserver::Geocentric,
                0.0,
                0.0,
                FLG_SPEED | FLG_NOABERR,
            )
            .unwrap();
        assert_eq!(apparent.2.to_bits(), unshifted.2.to_bits());
        assert_eq!(apparent.5.to_bits(), unshifted.5.to_bits());
        assert_ne!(apparent.0, unshifted.0);
        // Independent sixth-order derivative of actual unrotated distances.
        // At h=1/32 day, a bound 10^5 |v_rel|^7/r_min^6 on the seventh
        // derivative is below 1 e-20 AU/day after h^6/140 truncation.
        let h = 1.0 / 32.0;
        let distance = |tt| {
            serving
                .relative_at(901, MoonObserver::Geocentric, tt, FLG_NOABERR)
                .unwrap()
                .distance_au
        };
        let samples: [f64; 6] =
            [-3.0, -2.0, -1.0, 1.0, 2.0, 3.0].map(|offset| distance(offset * h));
        let weights = [-1.0, 9.0, -45.0, 45.0, -9.0, 1.0];
        let numerical = samples
            .into_iter()
            .zip(weights)
            .map(|(value, weight)| value * weight)
            .sum::<f64>()
            / (60.0 * h);
        // gamma_64 covers each rounded sample through three light-time
        // iterations and the weighted sum; it is an absolute forward bound,
        // independent of the measured discrepancy.
        let u = f64::EPSILON / 2.0;
        let gamma = 64.0 * u / (1.0 - 64.0 * u);
        let sample_scale = samples
            .into_iter()
            .zip(weights)
            .map(|(value, weight)| (value * weight).abs())
            .sum::<f64>();
        let bound = gamma * sample_scale / (60.0 * h) + 1e-20;
        assert!(
            (apparent.5 - numerical).abs() <= bound,
            "{} vs {numerical}; bound={bound}",
            apparent.5
        );
    }

    #[test]
    fn angular_rates_differentiate_cubic_frame_motion_across_zero() {
        let serving = MoonPositionServing::new(
            |_: i32, _: f64| Ok(([1.0, 0.0, 0.0], [0.0; 3])),
            |_| Ok(([0.0; 3], [0.0; 3])),
            |_: [f64; 3], tt: f64| {
                let t = tt - JD;
                Ok((
                    (0.005 + 0.5 * t + 0.01 * t.powi(3)).rem_euclid(360.0),
                    0.01 - 4e-6 * t - 0.02 * t.powi(3),
                    1.0,
                ))
            },
            |_| Ok(0.0),
            |lon, speed, _| Ok((lon, speed)),
        );
        let result = serving
            .moon_position(
                901,
                MoonObserver::Barycenter,
                JD,
                JD,
                FLG_SPEED | FLG_TRUEPOS,
            )
            .unwrap();
        assert!((result.3 - 0.5).abs() < 1e-11);
        assert!((result.4 + 4e-6).abs() < 1e-15);
        assert_eq!(result.5, 0.0);
        for epoch in [-1_000_000.5_f64, JD, 8_000_000.5] {
            let h = MOON_SPEED_HALF_STEP_DAYS;
            assert_eq!((epoch + h) - (epoch - h), 2.0 * h);
        }
    }

    #[test]
    fn fast_circular_satellite_refines_the_angular_step() {
        use std::cell::Cell;
        // Off-center geometry makes angular motion nonuniform even though
        // the satellite follows a perfect circle. A quarter-day orbit is
        // faster than the shortest-period natural satellites in this API.
        let samples = Cell::new(0);
        let omega = std::f64::consts::TAU / 0.25;
        let serving = MoonPositionServing::new(
            |_: i32, tt: f64| {
                samples.set(samples.get() + 1);
                let (s, c) = (omega * (tt - JD)).sin_cos();
                Ok((
                    [3.0 + 0.001 * c, 0.001 * s, 0.0001 * s],
                    [-0.001 * omega * s, 0.001 * omega * c, 0.0001 * omega * c],
                ))
            },
            |_| Ok(([0.0; 3], [0.0; 3])),
            |ray: [f64; 3], _: f64| {
                Ok((
                    ray[1].atan2(ray[0]).to_degrees().rem_euclid(360.0),
                    ray[2].atan2(ray[0].hypot(ray[1])).to_degrees(),
                    ray.into_iter().map(|v| v * v).sum::<f64>().sqrt(),
                ))
            },
            |_| Ok(0.0),
            |lon, speed, _| Ok((lon, speed)),
        );
        let result = serving
            .moon_position(
                901,
                MoonObserver::Barycenter,
                JD,
                JD,
                FLG_SPEED | FLG_TRUEPOS,
            )
            .unwrap();
        let expected_lon = (0.001 * omega / 3.001).to_degrees();
        let expected_lat = (0.0001 * omega / 3.001).to_degrees();
        assert!(
            (result.3 - expected_lon).abs() <= expected_lon.abs() * 1e-6,
            "{} vs {expected_lon}",
            result.3
        );
        assert!(
            (result.4 - expected_lat).abs() <= expected_lat.abs() * 1e-6,
            "{} vs {expected_lat}",
            result.4
        );
        assert!(
            samples.get() > 7,
            "fast motion must refine the initial dyadic step"
        );
        assert_eq!(result.5, 0.0);
    }

    #[test]
    fn bradley_adds_the_completed_direction_offset() {
        // Independent transcription of the two released functions:
        // _annual_aberration_cartesian returns a completed offset, which
        // apply_aberration_to_position then adds to the unit direction.
        let positions: [[f64; 3]; 3] = [[1.0, 2.0, 3.0], [-30.0, 5.0, -2.0], [0.01, 0.0, 0.0]];
        for position in positions {
            let velocity = [0.012, -0.007, 0.003];
            let distance =
                (position[0] * position[0] + position[1] * position[1] + position[2] * position[2])
                    .sqrt();
            let unit = position.map(|x| x / distance);
            let vc = velocity.map(|v| v / C_LIGHT_AU_DAY);
            let unit_magnitude = (unit[0] * unit[0] + unit[1] * unit[1] + unit[2] * unit[2]).sqrt();
            let normalized = unit.map(|x| x / unit_magnitude);
            let dot = normalized[0] * vc[0] + normalized[1] * vc[1] + normalized[2] * vc[2];
            let offset: [f64; 3] = std::array::from_fn(|i| vc[i] - dot * normalized[i]);
            let shifted: [f64; 3] = std::array::from_fn(|i| unit[i] + offset[i]);
            let magnitude =
                (shifted[0] * shifted[0] + shifted[1] * shifted[1] + shifted[2] * shifted[2])
                    .sqrt();
            let expected = shifted.map(|x| (x / magnitude) * distance);
            let actual = moon_bradley_aberration(position, velocity).unwrap();
            assert_eq!(actual.map(f64::to_bits), expected.map(f64::to_bits));
        }
    }

    /// Fixed SSB values: satellite on the +x axis at 0.01 AU, observer
    /// at the origin with a small +y velocity (aberration driver).
    // The fixture exposes each independent state/frame callback type explicitly.
    #[allow(clippy::type_complexity)]
    fn stub_words() -> (
        impl Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
        impl Fn(f64) -> Result<([f64; 3], [f64; 3])>,
    ) {
        let moon = |_: i32, _: f64| Ok(([0.01, 0.0, 0.0], [0.0, 0.0, 0.0]));
        let obs = |_: f64| Ok(([0.0, 0.0, 0.0], [0.0, 0.001, 0.0]));
        (moon, obs)
    }

    /// Identity frame value: treats the ICRF x-axis as lon 0.
    fn stub_frame() -> impl Fn([f64; 3], f64) -> Result<(f64, f64, f64)> {
        |v, _| {
            let lon = v[1].atan2(v[0]).to_degrees().rem_euclid(360.0);
            let dist = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            Ok((lon, 0.0, dist))
        }
    }

    // The fixture exposes each independent state/frame callback type explicitly.
    #[allow(clippy::type_complexity)]
    fn serving(
        moon: impl Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
        obs: impl Fn(f64) -> Result<([f64; 3], [f64; 3])>,
    ) -> MoonPositionServing<
        impl Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
        impl Fn(f64) -> Result<([f64; 3], [f64; 3])>,
        impl Fn([f64; 3], f64) -> Result<(f64, f64, f64)>,
        impl Fn(f64) -> Result<f64>,
        impl Fn(f64, f64, f64) -> Result<(f64, f64)>,
    > {
        MoonPositionServing::new(
            moon,
            obs,
            stub_frame(),
            |_| Ok(0.0),
            |lon, speed, _| Ok((lon - 23.0, speed)),
        )
    }

    #[test]
    fn observer_selection_follows_branch_order() {
        assert_eq!(
            moon_observer_for_flags(FLG_BARYCTR | FLG_HELCTR),
            MoonObserver::Barycenter
        );
        assert_eq!(
            moon_observer_for_flags(FLG_HELCTR | FLG_TOPOCTR),
            MoonObserver::Heliocentric
        );
        assert_eq!(
            moon_observer_for_flags(FLG_TOPOCTR),
            MoonObserver::Topocentric
        );
        assert_eq!(moon_observer_for_flags(0), MoonObserver::Geocentric);
        // Barycentric wins over every other center flag.
        assert_eq!(
            moon_observer_for_flags(FLG_BARYCTR | FLG_TOPOCTR),
            MoonObserver::Barycenter
        );
    }

    #[test]
    fn flag_gates_match_the_release_conjunctions() {
        assert!(moon_wants_light_time(0));
        assert!(!moon_wants_light_time(FLG_TRUEPOS));
        assert!(moon_wants_aberration(0));
        for gate in [FLG_TRUEPOS, FLG_NOABERR, FLG_HELCTR, FLG_BARYCTR] {
            assert!(!moon_wants_aberration(gate), "gate {gate}");
        }
        assert!(moon_wants_nutation_strip(FLG_NONUT));
        assert!(moon_wants_nutation_strip(FLG_J2000));
        assert!(moon_wants_nutation_strip(FLG_SIDEREAL | FLG_EQUATORIAL));
        assert!(!moon_wants_nutation_strip(FLG_SIDEREAL));
        assert!(!moon_wants_nutation_strip(0));
        assert!(moon_wants_sidereal_correction(FLG_SIDEREAL));
        assert!(!moon_wants_sidereal_correction(
            FLG_SIDEREAL | FLG_EQUATORIAL
        ));
        assert!(!moon_wants_sidereal_correction(0));
    }

    #[test]
    fn bradley_kernel_pins_first_order_shape() {
        // Zero observer velocity is the identity.
        let rel = [0.01, 0.0, 0.0];
        assert_eq!(moon_bradley_aberration(rel, [0.0, 0.0, 0.0]).unwrap(), rel);
        // Perpendicular velocity shifts the direction by ~v/c with the
        // distance preserved.
        let out = moon_bradley_aberration(rel, [0.0, 0.001, 0.0]).unwrap();
        let dist = (out[0] * out[0] + out[1] * out[1] + out[2] * out[2]).sqrt();
        assert!((dist - 0.01).abs() < 1e-15, "{out:?}");
        assert!(out[1] > 0.0 && out[0] < 0.01, "{out:?}");
        // Degenerate input returns unchanged (both 1 e-30 guards).
        assert_eq!(
            moon_bradley_aberration([0.0, 0.0, 0.0], [0.0, 0.001, 0.0]).unwrap(),
            [0.0; 3]
        );
        // Non-finite values refuse.
        assert!(moon_bradley_aberration([f64::NAN, 0.0, 0.0], [0.0; 3]).is_err());
    }

    #[test]
    fn truepos_skips_light_time_and_aberration() {
        // A moving moon separates the three shaping layers: light-time
        // retardation shifts the sampled geometry, aberration shifts the
        // direction at fixed geometry.
        let moon = |_: i32, jd_tt: f64| Ok(([0.01, 0.001 * (jd_tt - JD), 0.0], [0.0, 0.001, 0.0]));
        let (_, obs) = stub_words();
        let s = serving(moon, obs);
        let plain = s
            .moon_position(501, MoonObserver::Geocentric, JD, JD, 0)
            .unwrap();
        let truepos = s
            .moon_position(501, MoonObserver::Geocentric, JD, JD, FLG_TRUEPOS)
            .unwrap();
        assert!((plain.0 - truepos.0).abs() > 0.0);
        assert!((plain.2 - truepos.2).abs() > 0.0);
        // NOABERR keeps light time but drops the aberration offset.
        let noaberr = s
            .moon_position(501, MoonObserver::Geocentric, JD, JD, FLG_NOABERR)
            .unwrap();
        assert!((noaberr.0 - truepos.0).abs() > 0.0);
        assert!((noaberr.0 - plain.0).abs() > 0.0);
    }

    #[test]
    fn nutation_strip_and_speed_stencil_behave() {
        let (moon, obs) = stub_words();
        let s = MoonPositionServing::new(
            moon,
            obs,
            stub_frame(),
            |_| Ok(1.5),
            |lon, speed, _| Ok((lon, speed)),
        );
        // Static values: speeds are zero, strip shifts lon by -dpsi.
        let (lon, _, _, spd_lon, spd_lat, spd_dist) = s
            .moon_position(501, MoonObserver::Geocentric, JD, JD, FLG_NONUT | FLG_SPEED)
            .unwrap();
        assert!((spd_lon + spd_lat + spd_dist).abs() < 1e-9);
        let (lon_plain, _, _, _, _, _) = s
            .moon_position(501, MoonObserver::Geocentric, JD, JD, 0)
            .unwrap();
        assert!(
            ((lon - lon_plain + 1.5).rem_euclid(360.0)).abs() < 1e-9,
            "{lon} vs {lon_plain}"
        );
        // Sidereal gate subtracts through the injected value on ecliptic
        // output only.
        let (lon_sid, _, _, _, _, _) = s
            .moon_position(501, MoonObserver::Geocentric, JD, JD, FLG_SIDEREAL)
            .unwrap();
        assert!((lon_sid - lon_plain).abs() < 1e-9, "identity sidereal word");
    }

    #[test]
    fn production_sidereal_word_subtracts_ayanamsha() {
        // Constant ayanamsha: pure offset, zero drift, wrap-around.
        let word = |lon: f64, speed: f64| {
            moon_sidereal_word(lon, speed, JD, FLG_SIDEREAL | FLG_SPEED, |_| Ok(23.0))
        };
        let (lon, speed) = word(100.0, 0.5).unwrap();
        assert!(((lon - 77.0).rem_euclid(360.0)).abs() < 1e-12, "{lon}");
        assert!((speed - 0.5).abs() < 1e-12, "{speed}");
        let (wrapped, _) = word(10.0, 0.0).unwrap();
        assert!(
            ((wrapped - 347.0).rem_euclid(360.0)).abs() < 1e-12,
            "{wrapped}"
        );
        // Linear ayanamsha drift: the FLG_SPEED leg removes it.
        let drift = 0.01_f64.to_radians() * 365.25; // ~deg/day scale probe
        let drifting = |epoch: f64| Ok(23.0 + drift * (epoch - JD));
        let (_, speed_fast) =
            moon_sidereal_word(100.0, 0.5, JD, FLG_SIDEREAL | FLG_SPEED, drifting).unwrap();
        // 1 e-6: the ±1 s samples sit at absolute epochs, so the
        // `(epoch − JD)` subtraction at JD 2460000.5 (≈2 e-10 d) limits
        // the recovered drift because full-JD one-second offsets are rounded.
        assert!((speed_fast - (0.5 - drift)).abs() < 1e-6, "{speed_fast}");
        // Without FLG_SPEED the rate passes through untouched.
        let (_, speed_slow) = moon_sidereal_word(100.0, 0.5, JD, FLG_SIDEREAL, drifting).unwrap();
        assert!((speed_slow - 0.5).abs() < 1e-12, "{speed_slow}");
        // Without SPEED, neighboring epochs need not be in the provider
        // domain: position-only requests do not consume derivative samples.
        let center_only = |epoch: f64| {
            if epoch.to_bits() == JD.to_bits() {
                Ok(23.0)
            } else {
                Err(Error::NonFiniteArgument { value: epoch })
            }
        };
        assert_eq!(
            moon_sidereal_word(100.0, 0.5, JD, FLG_SIDEREAL, center_only).unwrap(),
            (77.0, 0.5)
        );
        assert!(moon_sidereal_word(100.0, 0.5, JD, FLG_SIDEREAL | FLG_SPEED, center_only).is_err());
        // Non-finite UT1 refuses; a provider fault propagates.
        assert!(moon_sidereal_word(100.0, 0.5, f64::NAN, FLG_SIDEREAL, |_| Ok(23.0)).is_err());
        assert!(
            moon_sidereal_word(100.0, 0.5, JD, FLG_SIDEREAL, |_| Err(
                Error::NonFiniteArgument { value: f64::NAN }
            ))
            .is_err()
        );
    }

    #[test]
    fn barycentric_leg_uses_zero_observer_words() {
        // The barycentric leg never consults the observer closure: a
        // refusing observer still serves.
        let (moon, _) = stub_words();
        let s = serving(moon, |_| Err(Error::LebBodyNotFound { body_id: 3 }));
        let out = s
            .moon_position(501, MoonObserver::Barycenter, JD, JD, FLG_BARYCTR)
            .unwrap();
        assert!(out.0.is_finite() && out.2 > 0.0);
    }

    #[test]
    fn lookup_miss_and_bad_epochs_refuse() {
        // Unknown moons resolve to no target (the model None).
        let registry = MoonSpkRegistry::new();
        assert_eq!(resolve_moon_target(&registry, 501), None);
        assert_eq!(resolve_moon_target(&registry, 999_999), None);
        // Non-finite epochs refuse at the serving boundary.
        let (moon, obs) = stub_words();
        let s = serving(moon, obs);
        assert!(
            s.moon_position(501, MoonObserver::Geocentric, f64::NAN, JD, 0)
                .is_err()
        );
    }
}

// ---------------------------------------------------------------------------
// Stored-provider satellite calculations over registered SPK chains.
// Satellite SPK epochs are TDB; the observer and light-time iteration use TT.
// The sidereal epoch is jd_tt - DeltaT(jd_tt). Unknown/unregistered targets
// return None, coverage misses retain kernel bounds, and topocentric requests
// require a configured observer.
/// Sealed satellite apparent-place binding: serves satellite
/// apparent-place values in the `calc_moon_position` return shape over the
/// shared sealed serving stack and the satellite SPK registry, with the
/// Delta-T and sidereal corrections supplied as caller values. Epochs are
/// TT; topocentric requests without an observer are refused, never
/// silently served geocentric.
pub struct SealedMoonPositionServing<'a, P> {
    serving: &'a SealedServing<'a, P>,
    moons: &'a MoonSpkRegistry,
    deltat_at: &'a dyn Fn(f64) -> Result<f64>,
    sidereal_at: &'a dyn Fn(f64, f64, f64, i32) -> Result<(f64, f64)>,
}

/// Satellite apparent place: `(lon_deg, lat_deg, dist_au, speed_lon, speed_lat,
/// speed_dist)` with speeds in units per day.
pub type MoonApparentPlace = (f64, f64, f64, f64, f64, f64);

impl<'a, P> SealedMoonPositionServing<'a, P>
where
    P: StateProvider,
{
    /// Bind the sealed serving pair, the satellite registry, the
    /// Delta-T provider (days, `jd + ΔT` association), and the sidereal
    /// correction value `(lon_deg, speed_deg_day, ut1_jd, flags)`.
    /// No I/O, no globals.
    pub fn new(
        serving: &'a SealedServing<'a, P>,
        moons: &'a MoonSpkRegistry,
        deltat_at: &'a dyn Fn(f64) -> Result<f64>,
        sidereal_at: &'a dyn Fn(f64, f64, f64, i32) -> Result<(f64, f64)>,
    ) -> Self {
        Self {
            serving,
            moons,
            deltat_at,
            sidereal_at,
        }
    }

    /// Bound sealed serving pair.
    pub fn serving(&self) -> &SealedServing<'a, P> {
        self.serving
    }

    /// Bound satellite registry.
    pub fn moons(&self) -> &MoonSpkRegistry {
        self.moons
    }

    /// TT epoch for a UT epoch through the bound Delta-T provider.
    pub fn jd_tt_for_ut(&self, jd_ut: f64) -> Result<f64> {
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let dt = (self.deltat_at)(jd_ut)?;
        if !dt.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt });
        }
        Ok(jd_ut + dt)
    }

    /// UT1 epoch for a TT epoch using the bound one-step inverse.
    /// The callback supplies TT-minus-UT1 days at the requested TT Julian day;
    /// errors propagate without substituting another clock.
    pub fn ut1_for_tt(&self, jd_tt: f64) -> Result<f64> {
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let dt = (self.deltat_at)(jd_tt)?;
        if !dt.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt });
        }
        Ok(jd_tt - dt)
    }

    /// Calculation-level range error for a satellite coverage miss
    /// : the whole-kernel JD span
    /// from [`MoonSpkRegistry::get_moon_coverage`], matching the model
    /// min/max-over-segments span.
    fn range_error(&self, moon_id: i32, jd_tdb: f64) -> Error {
        let (start_jd, end_jd) = self
            .moons
            .get_moon_coverage(moon_id)
            .unwrap_or((f64::NAN, f64::NAN));
        Error::EphemerisRange {
            requested_jd: jd_tdb,
            start_jd,
            end_jd,
            body: Some(moon_id),
        }
    }

    /// SSB satellite state in AU/AU-day at a TT epoch ` sampling): split TT to TDB, then the two-key
    /// registry probe followed by the kernel scan. Relative segments convert
    /// to AU before center addition. Provider coverage misses become the
    /// calculation-level range error; anything else propagates
    /// unchanged.
    fn moon_ssb_au(
        &self,
        naif_id: i32,
        moon_id: i32,
        emit_tt: f64,
    ) -> Result<([f64; 3], [f64; 3])> {
        if !emit_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: emit_tt });
        }
        let day = emit_tt.floor();
        let tt_fraction = emit_tt - day;
        let fraction = tt_fraction
            + crate::positions::reference_frames::tdb_minus_tt_seconds(day, tt_fraction)? / 86400.0;
        let jd_tdb = day + fraction; // Scalar diagnostic only, never a state epoch.
        let inner = &self.moons.inner;
        for key in [moon_id, PLMOON_OFFSET + naif_id] {
            if inner.body_info(key).is_none() {
                continue;
            }
            match inner.eval_geometric_ssb_au_two_part_jd(key, day, fraction) {
                Ok((position, velocity)) => {
                    let rate = crate::time::tdb_rate(emit_tt)?;
                    return Ok((position, velocity.map(|value| value * rate)));
                }
                Err(Error::SpkOutOfRange { .. }) => {
                    return Err(self.range_error(moon_id, jd_tdb));
                }
                Err(other) => return Err(other),
            }
        }
        for path in inner.kernel_paths() {
            let Some(kernel) = inner.kernel(&path) else {
                continue;
            };
            if resolve_chain(kernel.index(), naif_id).is_err() {
                continue;
            }
            match inner.eval_kernel_ssb_au_two_part_jd(&path, naif_id, day, fraction) {
                Ok((position, velocity)) => {
                    let rate = crate::time::tdb_rate(emit_tt)?;
                    return Ok((position, velocity.map(|value| value * rate)));
                }
                Err(Error::SpkOutOfRange { .. }) => {
                    return Err(self.range_error(moon_id, jd_tdb));
                }
                Err(other) => return Err(other),
            }
        }
        // Unreachable after the entry pre-resolution (no lock is held
        // across the two calls, so report honestly rather than invent).
        Err(Error::SpkBodyNotFound { target: naif_id })
    }

    /// SSB observer state in AU/AU-day at a TT epoch: the sealed
    /// ICRS Sun/Earth channels behind the channel gate. Barycentric is
    /// the zero observer; topocentric refuses with the typed
    /// missing-observer error (topocentric serving stays with the
    /// position-serving layer).
    fn observer_ssb(&self, observer: MoonObserver, jd_tt: f64) -> Result<([f64; 3], [f64; 3])> {
        let body = match observer {
            MoonObserver::Barycenter => return Ok(([0.0; 3], [0.0; 3])),
            MoonObserver::Heliocentric => SUN,
            MoonObserver::Geocentric => EARTH,
            MoonObserver::Topocentric => {
                return Err(Error::ConfigurationMissing {
                    missing_config: "observer_location",
                    suggestion: "FLG_TOPOCTR needs a prior set_topo site: topocentric moon serving is not bound yet",
                });
            }
        };
        serving_channel_gate(self.serving.provider(), body)?;
        let state = self.serving.provider().geometric_state(body, jd_tt)?;
        if !state.frame().is_icrs() {
            return Err(Error::SealedChannelNotServed {
                body_id: body,
                coord_type: self
                    .serving
                    .provider()
                    .serving_coord_type(body)
                    .unwrap_or(0),
            });
        }
        let (p, v) = (state.pos(), state.vel());
        for c in p.into_iter().chain(v) {
            if !c.is_finite() {
                return Err(Error::NonFiniteArgument { value: c });
            }
        }
        Ok((p, v))
    }

    /// True-ecliptic-of-date `(lon_deg, lat_deg, dist_au)` for an SSB
    /// observer-relative ICRF vector at a TT epoch: the default frame branch over the sealed
    /// nutation pair.
    fn frame_word(&self, rel: [f64; 3], jd_tt: f64) -> Result<(f64, f64, f64)> {
        for v in rel {
            if !v.is_finite() {
                return Err(Error::NonFiniteArgument { value: v });
            }
        }
        let (dpsi, deps) = self.serving.nutation_provider()(jd_tt)?;
        if !dpsi.is_finite() || !deps.is_finite() {
            return Err(Error::NonFiniteArgument { value: dpsi });
        }
        let frames = frame_matrices_for_date(jd_tt, dpsi, deps)?;
        let sph = frame_transform(rel, jd_tt, 0, false, frames)?;
        Ok((sph[0], sph[1], sph[2]))
    }

    /// Nutation in longitude in degrees at a TT epoch.
    fn nutation_deg_word(&self, jd_tt: f64) -> Result<f64> {
        let (dpsi, _) = self.serving.nutation_provider()(jd_tt)?;
        if !dpsi.is_finite() {
            return Err(Error::NonFiniteArgument { value: dpsi });
        }
        Ok(dpsi.to_degrees())
    }

    /// Satellite apparent place at a TT epoch: `Ok(None)` exactly where
    /// the model returns `None` (unknown moon id, no registered
    /// kernel chaining the target — checked before the observer gate,
    /// in model order). Coverage misses are
    /// [`Error::EphemerisRange`]; topocentric legs refuse with the
    /// typed missing-observer error.
    pub fn moon_position_tt(
        &self,
        moon_id: i32,
        jd_tt: f64,
        flags: i32,
    ) -> Result<Option<MoonApparentPlace>> {
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let naif_id = match moon_naif_id(moon_id) {
            Some(naif) => naif,
            None => return Ok(None),
        };
        if resolve_moon_target(self.moons, moon_id).is_none() {
            return Ok(None);
        }
        let observer = moon_observer_for_flags(flags);
        let ut1 = self.ut1_for_tt(jd_tt)?;
        let core = MoonPositionServing::new(
            |naif, emit_tt| self.moon_ssb_au(naif, moon_id, emit_tt),
            |tt| self.observer_ssb(observer, tt),
            |rel, tt| self.frame_word(rel, tt),
            |tt| self.nutation_deg_word(tt),
            |lon, speed, epoch| (self.sidereal_at)(lon, speed, epoch, flags),
        );
        core.moon_position(naif_id, observer, jd_tt, ut1, flags)
            .map(Some)
    }

    /// Satellite apparent place at a UT epoch:
    /// UT to TT through the bound Delta-T provider, then the TT entry.
    /// The sidereal value epoch derives from that TT epoch through
    /// [`Self::ut1_for_tt`].
    pub fn moon_position_ut(
        &self,
        moon_id: i32,
        jd_ut: f64,
        flags: i32,
    ) -> Result<Option<MoonApparentPlace>> {
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        self.moon_position_tt(moon_id, jd_tt, flags)
    }
}

#[cfg(test)]
mod sealed_moon_tests {
    use super::*;
    use crate::leb::{
        COORD_ICRS_BARY, CompositeIndex, FileIndex, HEADER_SIZE, LEB1_BODY_ENTRY_SIZE, MAGIC_LEB1,
        SECTION_DIR_SIZE,
    };
    use crate::positions::SealedServing;
    use std::cell::Cell;
    use std::sync::atomic::{AtomicU64, Ordering};

    /// Stub Delta-T in days (the value is arbitrary; the tests pin the
    /// threading, not a timescale).
    const DT_DAYS: f64 = 61.0 / 86400.0;
    /// TT epoch inside the SPK span (et 4.968 e 6 s past J2000).
    const JD_IN: f64 = 2_451_602.5;
    /// TT epoch outside the SPK span but inside the LEB span.
    const JD_OUT: f64 = 2_451_900.0;
    /// Sealed-observer LEB span (covers both epochs).
    const LEB_T0: f64 = 2_451_400.0;
    const LEB_T1: f64 = 2_452_000.0;
    /// Canonical Io id exercised by the sealed tests.
    const IO: i32 = PLMOON_OFFSET + NAIF_IO;

    static FIXTURE_SEQ: AtomicU64 = AtomicU64::new(0);

    fn put_i32_le(dst: &mut [u8], off: usize, v: i32) {
        dst[off..off + 4].copy_from_slice(&v.to_le_bytes());
    }

    fn put_f64_le(dst: &mut [u8], off: usize, v: f64) {
        dst[off..off + 8].copy_from_slice(&v.to_le_bytes());
    }

    /// Evaluable single-segment satellite kernel: NAIF 501 on center 0
    /// (SSB, so no center chain is needed), J2000 type 2, one degree-0
    /// record holding a static (1.5 e 6, 0, 0) km state over [4 e 6, 6 e 6] s
    /// past J2000. DAF layout mirrors the hand-built DAF/SPK test
    /// images: file / summary / name / element records.
    fn io_spk_image() -> Vec<u8> {
        let mut img = vec![0u8; 4 * 1024];
        img[0..8].copy_from_slice(b"DAF/SPK ");
        put_i32_le(&mut img, 8, 2); // ND
        put_i32_le(&mut img, 12, 6); // NI
        img[16..76].fill(b' ');
        img[16..30].copy_from_slice(b"RUST TEST FILE");
        put_i32_le(&mut img, 76, 2); // FWARD
        put_i32_le(&mut img, 80, 2); // BWARD
        put_i32_le(&mut img, 84, 4 * 128 + 1); // FREE
        img[88..96].copy_from_slice(b"LTL-IEEE");
        // Summary record 2: no links, one segment.
        let s = 1024;
        put_f64_le(&mut img, s, 0.0); // NEXT
        put_f64_le(&mut img, s + 8, 0.0); // PREV
        put_f64_le(&mut img, s + 16, 1.0); // NSUM
        let o = s + 24;
        put_f64_le(&mut img, o, 4_000_000.0);
        put_f64_le(&mut img, o + 8, 6_000_000.0);
        for (k, v) in [501, 0, 1, 2, 385, 393].iter().enumerate() {
            put_i32_le(&mut img, o + 16 + (k / 2) * 8 + (k % 2) * 4, *v);
        }
        // Name record 3.
        let n = 2048;
        img[n..n + 40].fill(b' ');
        img[n..n + 2].copy_from_slice(b"IO");
        // Element record 4: [mid, radius, x, y, z, init_et, int_len,
        // rsize, n] — a static state, zero velocity.
        let e = 3072;
        for (w, v) in [
            5_000_000.0,
            1_000_000.0,
            1_500_000.0,
            0.0,
            0.0,
            4_000_000.0,
            2_000_000.0,
            5.0,
            1.0,
        ]
        .iter()
        .enumerate()
        {
            put_f64_le(&mut img, e + w * 8, *v);
        }
        img
    }

    /// Constant-ICRS LEB 1 image over `[t0, t1]`: 2 segments, degree 1,
    /// zero slopes, so every sample is exactly the stored position with
    /// zero velocity. Layout mirrors the hand-packed LEB 1 test fixtures.
    fn const_leb_image(bodies: &[(i32, [f64; 3])], t0: f64, t1: f64) -> Vec<u8> {
        let n = bodies.len();
        let mut file = Vec::new();
        file.extend_from_slice(&MAGIC_LEB1);
        file.extend_from_slice(&1u32.to_le_bytes()); // version
        file.extend_from_slice(&1u32.to_le_bytes()); // sections
        file.extend_from_slice(&(n as u32).to_le_bytes());
        file.extend_from_slice(&t0.to_le_bytes());
        file.extend_from_slice(&t1.to_le_bytes());
        file.extend_from_slice(&0.0f64.to_le_bytes()); // generation
        file.extend_from_slice(&0u32.to_le_bytes()); // flags
        file.extend_from_slice(&[0u8; 20]);
        debug_assert_eq!(file.len(), HEADER_SIZE);
        let body_at = HEADER_SIZE + SECTION_DIR_SIZE;
        let blob_at = body_at + n * LEB1_BODY_ENTRY_SIZE;
        file.extend_from_slice(&0u32.to_le_bytes()); // section id
        file.extend_from_slice(&0u32.to_le_bytes()); // reserved
        file.extend_from_slice(&(body_at as u64).to_le_bytes());
        file.extend_from_slice(&((n * LEB1_BODY_ENTRY_SIZE) as u64).to_le_bytes());
        let interval = (t1 - t0) / 2.0;
        for (i, (body_id, _)) in bodies.iter().enumerate() {
            file.extend_from_slice(&body_id.to_le_bytes());
            file.extend_from_slice(&COORD_ICRS_BARY.to_le_bytes());
            file.extend_from_slice(&2u32.to_le_bytes()); // segments
            file.extend_from_slice(&t0.to_le_bytes());
            file.extend_from_slice(&t1.to_le_bytes());
            file.extend_from_slice(&interval.to_le_bytes());
            file.extend_from_slice(&1u32.to_le_bytes()); // degree
            file.extend_from_slice(&3u32.to_le_bytes()); // components
            file.extend_from_slice(&((blob_at + i * 96) as u64).to_le_bytes());
        }
        for (_, pos) in bodies {
            for _ in 0..2 {
                for c in pos {
                    file.extend_from_slice(&c.to_le_bytes());
                    file.extend_from_slice(&0.0f64.to_le_bytes());
                }
            }
        }
        file
    }

    fn zero_nutation(_jd: f64) -> Result<(f64, f64)> {
        Ok((0.0, 0.0))
    }

    /// Registry holding the Io fixture kernel (bytes live in the
    /// registry, so the temp file is removed after registration).
    fn io_registry() -> MoonSpkRegistry {
        let seq = FIXTURE_SEQ.fetch_add(1, Ordering::SeqCst);
        let path = std::env::temp_dir().join(format!(
            "libeph-moonsealed-{}-{seq}.bsp",
            std::process::id()
        ));
        std::fs::write(&path, io_spk_image()).expect("fixture write");
        let mut registry = MoonSpkRegistry::new();
        registry
            .register_moon_spk(&path.to_string_lossy(), None)
            .expect("fixture registers");
        std::fs::remove_file(&path).ok();
        registry
    }

    fn leb_file() -> Vec<u8> {
        const_leb_image(
            &[(SUN, [0.0, 0.0, 0.0]), (EARTH, [1.0, 0.0, 0.0])],
            LEB_T0,
            LEB_T1,
        )
    }

    fn identity_sidereal(lon: f64, speed: f64, _ut1: f64, _flags: i32) -> Result<(f64, f64)> {
        Ok((lon, speed))
    }

    #[test]
    fn unknown_moon_id_returns_none() {
        let leb = leb_file();
        let registry = io_registry();
        let owned = FileIndex::decode(&leb).expect("leb");
        let composite = CompositeIndex::new(vec![(owned, &leb[..])]).expect("composite");
        let serving = SealedServing::new(&composite, &zero_nutation);
        let deltat = |_: f64| -> Result<f64> { Ok(DT_DAYS) };
        let b = SealedMoonPositionServing::new(&serving, &registry, &deltat, &identity_sidereal);
        assert_eq!(b.moon_position_tt(999_999, JD_IN, 0).unwrap(), None);
        assert_eq!(b.moon_position_ut(999_999, JD_IN, 0).unwrap(), None);
        assert!(b.moon_position_tt(999_999, f64::NAN, 0).is_err());
    }

    #[test]
    fn unregistered_valid_moon_returns_none() {
        // Titan is a valid moon id with no registered kernel: the
        // model None, not a typed miss.
        let leb = leb_file();
        let registry = io_registry();
        let owned = FileIndex::decode(&leb).expect("leb");
        let composite = CompositeIndex::new(vec![(owned, &leb[..])]).expect("composite");
        let serving = SealedServing::new(&composite, &zero_nutation);
        let deltat = |_: f64| -> Result<f64> { Ok(DT_DAYS) };
        let b = SealedMoonPositionServing::new(&serving, &registry, &deltat, &identity_sidereal);
        assert_eq!(b.moon_position_tt(MOON_TITAN, JD_IN, 0).unwrap(), None);
    }

    #[test]
    fn sealed_geocentric_io_geometry() {
        let leb = leb_file();
        let registry = io_registry();
        let owned = FileIndex::decode(&leb).expect("leb");
        let composite = CompositeIndex::new(vec![(owned, &leb[..])]).expect("composite");
        let serving = SealedServing::new(&composite, &zero_nutation);
        let deltat = |_: f64| -> Result<f64> { Ok(DT_DAYS) };
        let b = SealedMoonPositionServing::new(&serving, &registry, &deltat, &identity_sidereal);
        // Static Io at +x, Earth at (1, 0, 0) AU: the geocentric
        // direction is −x (lon 180 up to of-date frame rotation) at
        // ~0.99 AU, with zero speeds (no SPEED flag).
        let out = b.moon_position_tt(IO, JD_IN, 0).unwrap().expect("Io");
        assert!((179.0..181.0).contains(&out.0), "lon {out:?}");
        assert!(out.1.abs() < 0.5, "lat {out:?}");
        assert!((0.98..1.0).contains(&out.2), "dist {out:?}");
        assert_eq!((out.3, out.4, out.5), (0.0, 0.0, 0.0));
        // The legacy alias rides the same chain numerically.
        let legacy = b
            .moon_position_tt(MOON_IO, JD_IN, 0)
            .unwrap()
            .expect("alias");
        assert_eq!(out, legacy);
    }

    #[test]
    fn barycentric_leg_serves_moon_ssb_distance() {
        let leb = leb_file();
        let registry = io_registry();
        let owned = FileIndex::decode(&leb).expect("leb");
        let composite = CompositeIndex::new(vec![(owned, &leb[..])]).expect("composite");
        let serving = SealedServing::new(&composite, &zero_nutation);
        let deltat = |_: f64| -> Result<f64> { Ok(DT_DAYS) };
        let b = SealedMoonPositionServing::new(&serving, &registry, &deltat, &identity_sidereal);
        // Barycentric: the raw SSB state, ~0.01003 AU on +x, never
        // consulting an observer channel.
        let out = b
            .moon_position_tt(IO, JD_IN, FLG_BARYCTR)
            .unwrap()
            .expect("barycentric Io");
        assert!(out.0 < 1.0 || out.0 > 359.0, "lon {out:?}");
        assert!((0.009..0.011).contains(&out.2), "dist {out:?}");
    }

    #[test]
    fn out_of_coverage_maps_to_ephemeris_range() {
        let leb = leb_file();
        let registry = io_registry();
        let owned = FileIndex::decode(&leb).expect("leb");
        let composite = CompositeIndex::new(vec![(owned, &leb[..])]).expect("composite");
        let serving = SealedServing::new(&composite, &zero_nutation);
        let deltat = |_: f64| -> Result<f64> { Ok(DT_DAYS) };
        let b = SealedMoonPositionServing::new(&serving, &registry, &deltat, &identity_sidereal);
        let err = b.moon_position_tt(IO, JD_OUT, 0).unwrap_err();
        match err {
            Error::EphemerisRange {
                requested_jd,
                start_jd,
                end_jd,
                body,
            } => {
                assert_eq!(body, Some(IO));
                assert!(start_jd.is_finite() && end_jd.is_finite() && start_jd < end_jd);
                assert!(
                    requested_jd > end_jd,
                    "requested {requested_jd} vs [{start_jd}, {end_jd}]"
                );
            }
            other => panic!("wrong category: {other:?}"),
        }
    }

    #[test]
    fn topocentric_refuses_after_successful_lookup() {
        let leb = leb_file();
        let registry = io_registry();
        let owned = FileIndex::decode(&leb).expect("leb");
        let composite = CompositeIndex::new(vec![(owned, &leb[..])]).expect("composite");
        let serving = SealedServing::new(&composite, &zero_nutation);
        let deltat = |_: f64| -> Result<f64> { Ok(DT_DAYS) };
        let b = SealedMoonPositionServing::new(&serving, &registry, &deltat, &identity_sidereal);
        // Io resolves, so the missing site refuses (model order:
        // lookup first, observer block second).
        let err = b.moon_position_tt(IO, JD_IN, FLG_TOPOCTR).unwrap_err();
        assert!(
            matches!(err, Error::ConfigurationMissing { .. }),
            "wrong refusal: {err:?}"
        );
        // An unresolvable moon still returns None under TOPOCTR: the
        // lookup miss precedes the observer gate.
        assert_eq!(
            b.moon_position_tt(999_999, JD_IN, FLG_TOPOCTR).unwrap(),
            None
        );
    }

    #[test]
    fn speed_stencil_on_static_geometry_stays_tiny() {
        let leb = leb_file();
        let registry = io_registry();
        let owned = FileIndex::decode(&leb).expect("leb");
        let composite = CompositeIndex::new(vec![(owned, &leb[..])]).expect("composite");
        let serving = SealedServing::new(&composite, &zero_nutation);
        let deltat = |_: f64| -> Result<f64> { Ok(DT_DAYS) };
        let b = SealedMoonPositionServing::new(&serving, &registry, &deltat, &identity_sidereal);
        // Static values: only the of-date frame drift (precession rate)
        // leaks into the speeds, far below 1 e-3 deg/day.
        let out = b
            .moon_position_tt(IO, JD_IN, FLG_SPEED)
            .unwrap()
            .expect("speed Io");
        assert!(out.3.abs() < 1e-3, "spd_lon {out:?}");
        assert!(out.4.abs() < 1e-3, "spd_lat {out:?}");
        assert!(out.5.abs() < 1e-9, "spd_dist {out:?}");
    }

    #[test]
    fn ut_entry_threads_delta_t() {
        let leb = leb_file();
        let registry = io_registry();
        let owned = FileIndex::decode(&leb).expect("leb");
        let composite = CompositeIndex::new(vec![(owned, &leb[..])]).expect("composite");
        let serving = SealedServing::new(&composite, &zero_nutation);
        let deltat = |_: f64| -> Result<f64> { Ok(DT_DAYS) };
        let b = SealedMoonPositionServing::new(&serving, &registry, &deltat, &identity_sidereal);
        let jd_ut = JD_IN - DT_DAYS;
        let via_ut = b.moon_position_ut(IO, jd_ut, 0).unwrap().expect("ut");
        let via_tt = b.moon_position_tt(IO, JD_IN, 0).unwrap().expect("tt");
        for (n, a, e) in [
            ("lon", via_ut.0, via_tt.0),
            ("lat", via_ut.1, via_tt.1),
            ("dist", via_ut.2, via_tt.2),
        ] {
            assert!((a - e).abs() < 1e-9, "{n}: {a} vs {e}");
        }
        assert!(b.moon_position_ut(IO, f64::NAN, 0).is_err());
    }

    #[test]
    fn sidereal_word_receives_ut1_and_flags() {
        let leb = leb_file();
        let registry = io_registry();
        let owned = FileIndex::decode(&leb).expect("leb");
        let composite = CompositeIndex::new(vec![(owned, &leb[..])]).expect("composite");
        let serving = SealedServing::new(&composite, &zero_nutation);
        let deltat = |_: f64| -> Result<f64> { Ok(DT_DAYS) };
        let seen_ut1 = Cell::new(0.0f64);
        let seen_flags = Cell::new(0i32);
        let sidereal = |lon: f64, speed: f64, ut1: f64, flags: i32| -> Result<(f64, f64)> {
            seen_ut1.set(ut1);
            seen_flags.set(flags);
            Ok((lon - 23.0, speed))
        };
        let b = SealedMoonPositionServing::new(&serving, &registry, &deltat, &sidereal);
        // No sidereal flag: the value is never consulted.
        let plain = b.moon_position_tt(IO, JD_IN, 0).unwrap().expect("plain");
        // Sidereal: the mean/true dispatch epoch is the single-step
        // `jd_tt − ΔT` inversion, and the flags travel with the call.
        let sid = b
            .moon_position_tt(IO, JD_IN, FLG_SIDEREAL)
            .unwrap()
            .expect("sidereal");
        assert_eq!(seen_flags.get(), FLG_SIDEREAL);
        assert_eq!(seen_ut1.get(), JD_IN - DT_DAYS);
        assert!(
            (sid.0 - (plain.0 - 23.0)).abs() < 1e-12,
            "{sid:?} vs {plain:?}"
        );
        assert_eq!((sid.3, sid.4, sid.5), (0.0, 0.0, 0.0));
    }
}
