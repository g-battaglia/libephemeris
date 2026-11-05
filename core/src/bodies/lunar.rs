//! Mean, osculating and interpolated lunar nodes and apsides.
//!
//! Mean-point geometry uses the Simon et al. (1994) Delaunay polynomials adopted
//! by IERS Conventions (2003), with fused Horner evaluation. ERFA numeric models
//! are SOFA-derived, Copyright NumFOCUS, under three-clause BSD-style ERFA
//! conditions; see `THIRD_PARTY_NOTICES.md` and `docs/sources.md`.
//! Conventional orbit parameters use Meeus/ELP values and the IAU astronomical
//! unit. No measured position is fitted at runtime.
//!
//! True points reduce caller-owned same-epoch Moon/Earth position and velocity
//! to angular momentum and eccentricity geometry. Interpolated apsides use the
//! bundled DE440-fitted periodic and residual tables. Model-construction helpers
//! include basis evaluation, cubic splines, candidate scans and passage refinement;
//! they do not acquire a kernel implicitly.
//!
//! Epochs are TT Julian days. Public angles are degrees, ranges AU and rates
//! corresponding units per TT day. Mean points begin on the mean ecliptic of date;
//! flag-aware composition applies the requested output frame. Explicit
//! [`LunarCoverage`] selects whether a provider window is enforced. Nonfinite
//! inputs and required coverage failures return typed errors. Analytic validity
//! bands and interpolated-table windows are model limits, not accuracy guarantees.

use crate::constants::{
    FLG_BARYCTR, FLG_EQUATORIAL, FLG_HELCTR, FLG_J2000, FLG_NONUT, FLG_SIDEREAL, FLG_SPEED,
    FLG_SPEED3, FLG_TOPOCTR, INTP_APOG, INTP_PERG, MEAN_APOG, MEAN_NODE, OSCU_APOG, TRUE_NODE,
};
use crate::erfa::{
    FA_D03 as POLY_D, FA_F03 as POLY_F, FA_L03 as POLY_L, FA_LP03 as POLY_LP,
    FA_OM03 as POLY_OMEGA, fa_poly_rad as iers_argument_rad,
};
use crate::errors::Error;
#[path = "lunar_context.rs"]
mod context;
pub use context::{
    ecliptic_moon_state_of, mean_lunar_state_with_providers, sealed_true_state_at,
    true_lunar_state_with_providers,
};

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// TT Julian Day of J2000.0.
const J2000_JD: f64 = 2451545.0;
/// Days per Julian century.
const CENTURY_DAYS: f64 = 36525.0;
/// Half-span of the centered argument-rate stencil, in days.
const RATE_HALF_SPAN_DAYS: f64 = 0.5;
/// Exact IAU 2012 astronomical unit in km.
const AU_KM: f64 = 149_597_870.7;
/// Conventional mean lunar distance in km.
const MEAN_LUNAR_DISTANCE_KM: f64 = 384_400.0;
/// Conventional mean lunar eccentricity.
const MEAN_LUNAR_ECCENTRICITY: f64 = 0.0549;

/// Mean inclination of the lunar orbit to the ecliptic, in degrees
pub const INCLINATION_DEG: f64 = 5.145396;
/// Mean lunar distance in AU.
pub const MEAN_NODE_DISTANCE_AU: f64 = MEAN_LUNAR_DISTANCE_KM / AU_KM;
/// Mean lunar-apogee distance in AU.
pub const MEAN_APOGEE_DISTANCE_AU: f64 =
    MEAN_LUNAR_DISTANCE_KM * (1.0 + MEAN_LUNAR_ECCENTRICITY) / AU_KM;

/// Mean Earth-orbit eccentricity at J2000.
pub const EARTH_ECCENTRICITY_J2000: f64 = 0.016708634;
/// Linear Earth-eccentricity coefficient per century.
pub const EARTH_ECCENTRICITY_T: f64 = -0.000042037;
/// Quadratic Earth-eccentricity coefficient per century squared.
pub const EARTH_ECCENTRICITY_T2: f64 = -0.0000001267;

/// Optimal Meeus-polynomial window, in centuries from J2000.
pub const MEEUS_OPTIMAL_CENTURIES: f64 = 2.0;
/// Usable Meeus-polynomial window, in centuries from J2000.
pub const MEEUS_VALID_CENTURIES: f64 = 10.0;
/// Outer Meeus-polynomial window, in centuries from J2000.
pub const MEEUS_MAX_CENTURIES: f64 = 20.0;

// IAU fundamental arguments and frame models.

// ---------------------------------------------------------------------------
// Range policy
// ---------------------------------------------------------------------------

/// Explicit coverage policy for the mean-point range check in the
/// mean-kernel layer.
///
/// the model consults mutable global session state
/// (`_active_ephemeris_range`); this enum carries the resolved answer:
/// [`LunarCoverage::Unenforced`] mirrors the Horizons/LEB/auto-sealed
/// branches (the analytic model is not sampled from a kernel, so no
/// window narrows its validity), while [`LunarCoverage::Window`]
/// mirrors the loaded-JPL-kernel branch. A `(0, 0)` window is
/// unenforced, using the same physical model falsy `start_jd and end_jd` branch.
/// The caller (session/state layer) resolves which variant applies; the
/// geometry here never loads a kernel to answer a range question.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LunarCoverage {
    /// No kernel window narrows the analytic validity.
    Unenforced,
    /// Enforce a loaded-kernel window, in TT Julian days.
    Window {
        /// Inclusive lower bound in TT Julian days; two zero bounds disable the check.
        start_jd: f64,
        /// Inclusive upper bound in TT Julian days; two zero bounds disable the check.
        end_jd: f64,
    },
}

/// (non-finite epoch refusal plus the
/// window-miss escalation naming the requested mean body).
fn check_lunar_range(jd_tt: f64, body_id: i32, coverage: LunarCoverage) -> Result<(), Error> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    if let LunarCoverage::Window { start_jd, end_jd } = coverage {
        if start_jd != 0.0 || end_jd != 0.0 {
            if !start_jd.is_finite() {
                return Err(Error::NonFiniteArgument { value: start_jd });
            }
            if !end_jd.is_finite() {
                return Err(Error::NonFiniteArgument { value: end_jd });
            }
            if !(start_jd <= jd_tt && jd_tt <= end_jd) {
                return Err(Error::EphemerisRange {
                    requested_jd: jd_tt,
                    start_jd,
                    end_jd,
                    body: Some(body_id),
                });
            }
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// IERS fundamental arguments
// ---------------------------------------------------------------------------

/// TT Julian Day to Julian centuries since J2000.0.
pub fn centuries_since_j2000(jd_tt: f64) -> f64 {
    (jd_tt - J2000_JD) / CENTURY_DAYS
}

/// `(l, F, Omega)` in radians.
fn fundamental_l_f_node(jd_tt: f64) -> (f64, f64, f64) {
    let t = centuries_since_j2000(jd_tt);
    (
        iers_argument_rad(t, POLY_L),
        iers_argument_rad(t, POLY_F),
        iers_argument_rad(t, POLY_OMEGA),
    )
}

/// IERS 2003 lunar arguments `(D, M, M', F)` in radians.
/// No range check on this path, using the same physical model.
/// The explicit (D, M, M', F) order matches the trigonometric basis of the passage-fitted apsis model below.
pub fn lunar_delaunay_arguments(jd_tt: f64) -> (f64, f64, f64, f64) {
    let t = centuries_since_j2000(jd_tt);
    (
        iers_argument_rad(t, POLY_D),
        iers_argument_rad(t, POLY_LP),
        iers_argument_rad(t, POLY_L),
        iers_argument_rad(t, POLY_F),
    )
}

/// Signed principal angular difference in radians.
fn angle_delta(after: f64, before: f64) -> f64 {
    (after - before).sin().atan2((after - before).cos())
}

/// Smooth IERS argument rates `(l, F, Omega)` in radians per day
/// : centered one-day stencil,
/// far shorter than the 18.6-year nodal period and long enough to avoid
/// cancellation in binary64.
fn fundamental_argument_rates(jd_tt: f64) -> (f64, f64, f64) {
    let before = fundamental_l_f_node(jd_tt - RATE_HALF_SPAN_DAYS);
    let after = fundamental_l_f_node(jd_tt + RATE_HALF_SPAN_DAYS);
    let span = 2.0 * RATE_HALF_SPAN_DAYS;
    (
        angle_delta(after.0, before.0) / span,
        angle_delta(after.1, before.1) / span,
        angle_delta(after.2, before.2) / span,
    )
}

// ---------------------------------------------------------------------------
// Mean apse plane / node longitudes
// ---------------------------------------------------------------------------

/// Mean-apogee plane longitude and mean-node longitude in degrees
/// without a range check.
///
/// The apogee plane follows from the Delaunay identities `F = L - Omega`
/// and `l = L - varpi`: `varpi_apogee = Omega + F - l + pi`. Rust
/// `rem_euclid` is the nonnegative Euclidean remainder,
/// so the wrap lands in `[0, 360)` on both sides.
fn mean_apse_plane_node_unchecked(jd_tt: f64) -> (f64, f64) {
    let (lunar_anomaly, argument_latitude, node) = fundamental_l_f_node(jd_tt);
    let apogee_plane = node + argument_latitude - lunar_anomaly + std::f64::consts::PI;
    (
        apogee_plane.to_degrees().rem_euclid(360.0),
        node.to_degrees().rem_euclid(360.0),
    )
}

/// Mean-apogee plane longitude and mean-node longitude in degrees
/// , range-checked for the
/// `MEAN_APOG` body.
pub fn mean_lunar_apse_plane_node(
    jd_tt: f64,
    coverage: LunarCoverage,
) -> Result<(f64, f64), Error> {
    check_lunar_range(jd_tt, MEAN_APOG, coverage)?;
    Ok(mean_apse_plane_node_unchecked(jd_tt))
}

// ---------------------------------------------------------------------------
// Mean-apogee Cartesian state and polar shaping
// ---------------------------------------------------------------------------

/// Conventional mean-apogee Cartesian state of date: ecliptic AU and
/// AU/day.
///
/// The unit direction tilts the apogee-plane longitude by
/// [`INCLINATION_DEG`] about the node line; the velocity applies the
/// chain rule over the node rate and the argument rate
/// (`F_rate - l_rate`), all in the calculation order.
fn mean_apogee_cartesian_state(jd_tt: f64) -> (f64, f64, f64, f64, f64, f64) {
    let (lunar_anomaly, argument_latitude, node) = fundamental_l_f_node(jd_tt);
    let (lunar_anomaly_rate, argument_latitude_rate, node_rate) = fundamental_argument_rates(jd_tt);
    let argument = argument_latitude - lunar_anomaly + std::f64::consts::PI;
    let argument_rate = argument_latitude_rate - lunar_anomaly_rate;

    let inclination = INCLINATION_DEG.to_radians();
    let cos_i = inclination.cos();
    let sin_i = inclination.sin();
    let cos_node = node.cos();
    let sin_node = node.sin();
    let cos_argument = argument.cos();
    let sin_argument = argument.sin();

    let x_unit = cos_node * cos_argument - sin_node * sin_argument * cos_i;
    let y_unit = sin_node * cos_argument + cos_node * sin_argument * cos_i;
    let z_unit = sin_argument * sin_i;

    let dx_dnode = -y_unit;
    let dy_dnode = x_unit;
    let dx_dargument = -cos_node * sin_argument - sin_node * cos_argument * cos_i;
    let dy_dargument = -sin_node * sin_argument + cos_node * cos_argument * cos_i;
    let dz_dargument = cos_argument * sin_i;

    let distance = MEAN_APOGEE_DISTANCE_AU;
    (
        distance * x_unit,
        distance * y_unit,
        distance * z_unit,
        distance * (dx_dnode * node_rate + dx_dargument * argument_rate),
        distance * (dy_dnode * node_rate + dy_dargument * argument_rate),
        distance * dz_dargument * argument_rate,
    )
}

/// Cartesian position/velocity to spherical state and rates.
///
/// Inputs are ecliptic Cartesian AU and AU/day; outputs are longitude,
/// latitude, distance, degree/day angular rates and AU/day radial rate.
/// The latitude uses two-argument `atan2(z, rho)` exactly like the
/// the model (not the one-argument `atan` of the quotient: the two spell
/// different roundings at isolated inputs). The radial slot is `0.0`: the model computes the radial
/// speed only as an intermediate for the latitude rate and reports zero.
/// A vanishing cylindrical radius would divide by zero on both sides;
/// it is unreachable for this model (the tilted unit circle never meets
/// the ecliptic pole: `x = y = 0` would need `cos(arg) = 0` and
/// `sin(arg) = 0` together).
pub(crate) fn cartesian_to_polar_state(
    state: (f64, f64, f64, f64, f64, f64),
) -> (f64, f64, f64, f64, f64, f64) {
    let (x, y, z, vx, vy, vz) = state;
    let rho_squared = x * x + y * y;
    let rho = rho_squared.sqrt();
    let distance = (rho_squared + z * z).sqrt();
    let radial_speed = (x * vx + y * vy + z * vz) / distance;
    let longitude_speed = (x * vy - y * vx) / rho_squared;
    let latitude_speed = (vz * distance - z * radial_speed) / (distance * rho);
    (
        y.atan2(x).to_degrees().rem_euclid(360.0),
        z.atan2(rho).to_degrees(),
        distance,
        longitude_speed.to_degrees(),
        latitude_speed.to_degrees(),
        0.0,
    )
}

// ---------------------------------------------------------------------------
// Checked mean positions and states
// ---------------------------------------------------------------------------

/// Mean lunar-apogee longitude, latitude and distance without a range
/// check.
pub fn mean_lunar_apogee_position_unchecked(jd_tt: f64) -> (f64, f64, f64) {
    let (lon, lat, dist, _, _, _) = cartesian_to_polar_state(mean_apogee_cartesian_state(jd_tt));
    (lon, lat, dist)
}

/// Mean lunar-apogee longitude, latitude and distance.
pub fn mean_lunar_apogee_position(
    jd_tt: f64,
    coverage: LunarCoverage,
) -> Result<(f64, f64, f64), Error> {
    check_lunar_range(jd_tt, MEAN_APOG, coverage)?;
    Ok(mean_lunar_apogee_position_unchecked(jd_tt))
}

/// Standards-derived mean lunar-apogee polar state
/// `(lon, lat, dist, dlon, dlat, ddist)`.
///
/// `speed3` is accepted for the shared backend call signature and
/// ignored: this model returns its smooth standards-defined
/// instantaneous derivative, while public SPEED 3 callers rebuild the
/// derivative from fully transformed positions.
pub fn mean_lunar_apogee_state(
    jd_tt: f64,
    coverage: LunarCoverage,
    speed3: bool,
) -> Result<(f64, f64, f64, f64, f64, f64), Error> {
    let _ = speed3;
    check_lunar_range(jd_tt, MEAN_APOG, coverage)?;
    Ok(cartesian_to_polar_state(mean_apogee_cartesian_state(jd_tt)))
}

/// Mean lunar-node longitude, zero latitude and mean distance without a
/// range check.
pub fn mean_lunar_node_position_unchecked(jd_tt: f64) -> (f64, f64, f64) {
    let (_, _, node) = fundamental_l_f_node(jd_tt);
    (
        node.to_degrees().rem_euclid(360.0),
        0.0,
        MEAN_NODE_DISTANCE_AU,
    )
}

/// Mean lunar-node longitude, zero latitude and mean distance.
pub fn mean_lunar_node_position(
    jd_tt: f64,
    coverage: LunarCoverage,
) -> Result<(f64, f64, f64), Error> {
    check_lunar_range(jd_tt, MEAN_NODE, coverage)?;
    Ok(mean_lunar_node_position_unchecked(jd_tt))
}

/// ERFA/IERS mean lunar-node polar state. `speed3` is accepted and ignored, as for
/// the apogee state.
pub fn mean_lunar_node_state(
    jd_tt: f64,
    coverage: LunarCoverage,
    speed3: bool,
) -> Result<(f64, f64, f64, f64, f64, f64), Error> {
    let _ = speed3;
    check_lunar_range(jd_tt, MEAN_NODE, coverage)?;
    let (_, _, node) = fundamental_l_f_node(jd_tt);
    let (_, _, node_rate) = fundamental_argument_rates(jd_tt);
    Ok((
        node.to_degrees().rem_euclid(360.0),
        0.0,
        MEAN_NODE_DISTANCE_AU,
        node_rate.to_degrees(),
        0.0,
        0.0,
    ))
}

// ---------------------------------------------------------------------------
// Thin lunar.py mean wrappers
// ---------------------------------------------------------------------------

/// Independent ERFA/IERS mean lunar-node longitude.
pub fn calc_mean_lunar_node(jd_tt: f64, coverage: LunarCoverage) -> Result<f64, Error> {
    Ok(mean_lunar_node_position(jd_tt, coverage)?.0)
}

/// Independent ERFA/IERS mean lunar-apogee longitude.
pub fn calc_mean_lilith(jd_tt: f64, coverage: LunarCoverage) -> Result<f64, Error> {
    Ok(mean_lunar_apogee_position(jd_tt, coverage)?.0)
}

/// Independent mean lunar-apogee longitude and latitude.
pub fn calc_mean_lilith_with_latitude(
    jd_tt: f64,
    coverage: LunarCoverage,
) -> Result<(f64, f64), Error> {
    let (lon, lat, _) = mean_lunar_apogee_position(jd_tt, coverage)?;
    Ok((lon, lat))
}

/// Independent mean-node polar state on the mean ecliptic.
pub fn calc_mean_lunar_node_state(
    jd_tt: f64,
    coverage: LunarCoverage,
    speed3: bool,
) -> Result<(f64, f64, f64, f64, f64, f64), Error> {
    mean_lunar_node_state(jd_tt, coverage, speed3)
}

/// Independent mean-apogee polar state on the mean ecliptic.
pub fn calc_mean_lilith_state(
    jd_tt: f64,
    coverage: LunarCoverage,
    speed3: bool,
) -> Result<(f64, f64, f64, f64, f64, f64), Error> {
    mean_lunar_apogee_state(jd_tt, coverage, speed3)
}

// ---------------------------------------------------------------------------
// Meeus validity bands
// ---------------------------------------------------------------------------

/// Project validity band for a Julian-century distance from J2000.0
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MeeusBand {
    /// Within ±10 centuries: optimal polynomial range.
    Optimal,
    /// Beyond ±10 centuries up to ±20: degraded precision, still usable.
    Degraded,
    /// Beyond ±20 centuries: error may exceed one degree.
    OutsideValid,
}

/// Classify a Julian-century distance from J2000.0 with the exact
/// the model thresholds (boundary values stay in the inner band:
/// `> 10.0`, then `> 20.0`).
pub fn meeus_band(centuries_from_j2000: f64) -> MeeusBand {
    let abs_t = centuries_from_j2000.abs();
    if abs_t > MEEUS_MAX_CENTURIES {
        MeeusBand::OutsideValid
    } else if abs_t > MEEUS_VALID_CENTURIES {
        MeeusBand::Degraded
    } else {
        MeeusBand::Optimal
    }
}

// Osculating node and apsis geometry from a geocentric lunar state.
//
// Input position/velocity is AU and AU/day in the true ecliptic of date.
// The native state provider owns acquisition and frame selection. The conic
// reduction uses the Earth-Moon gravitational parameter below. Nonfinite
// states return a typed error; a degenerate zero state produces the documented
// undefined apsis values or node fallback rather than a physical orbit.

/// Earth-Moon system gravitational parameter in AU^3/day^2.
///
/// The two inputs are the JPL DE430 ephemeris constants (Folkner et
/// al. 2014, IPN Progress Report 42-196): `GM_Earth = 398600.435436`
/// km^3/s^2 and Earth/Moon mass ratio `81.3005691`. DE440 revised these
/// marginally (`GM_Earth = 398600.435507`, EMRAT `= 81.3005682...`), a
/// ~2 e-10 relative change the model calls negligible for this
/// osculating-element math, which is its single consumer. The spelling
/// keeps the calculation order (`**3`/`**2` evaluate to the
/// same values as the written products on this toolchain; the value is
/// pinned by test).
pub const GM_EARTH_MOON_AU3_DAY2: f64 = (398_600.435436 * (1.0 + 1.0 / 81.3005691))
    / (149_597_870.7 * 149_597_870.7 * 149_597_870.7)
    * (86_400.0 * 86_400.0);

/// Eccentricity floor for a defined node direction: at or below it the orbit is treated as circular and the
/// node distance falls back to the conventional mean distance.
pub const CIRCULAR_ECCENTRICITY_FLOOR: f64 = 1e-15;

/// Ordinary `TRUE_NODE` speed half-step in days.
pub const TRUE_NODE_HALF_STEP_DAYS: f64 = 0.05;
/// Ordinary `OSCU_APOG` speed half-step in days.
pub const OSCULATING_APSE_HALF_STEP_DAYS: f64 = 0.002;
/// Fine speed half-step in days for `SPEED3`/topocentric requests
pub const FINE_POINT_HALF_STEP_DAYS: f64 = 0.0001;

/// Select the true-point speed half-step: `fine` mirrors the
/// `FLG_SPEED3 | FLG_TOPOCTR` branch, `ordinary_half_step_days` is the
/// per-point ordinary step above. Non-finite steps are refused; the
/// caller owns positivity.
pub fn select_true_point_half_step(fine: bool, ordinary_half_step_days: f64) -> Result<f64, Error> {
    if !ordinary_half_step_days.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: ordinary_half_step_days,
        });
    }
    Ok(if fine {
        FINE_POINT_HALF_STEP_DAYS
    } else {
        ordinary_half_step_days
    })
}

/// Refuse a non-finite true-ecliptic Moon state value, naming the first
/// bad value.
fn check_ecliptic_state(r: (f64, f64, f64), v: (f64, f64, f64)) -> Result<(), Error> {
    for value in [r.0, r.1, r.2, v.0, v.1, v.2] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    Ok(())
}

/// Angular momentum `h = r x v` in the model component order
/// (shared by the node, apogee and perigee paths).
fn angular_momentum(r: (f64, f64, f64), v: (f64, f64, f64)) -> (f64, f64, f64) {
    (
        r.1 * v.2 - r.2 * v.1,
        r.2 * v.0 - r.0 * v.2,
        r.0 * v.1 - r.1 * v.0,
    )
}

/// True (osculating) lunar ascending-node longitude, zero latitude and
/// node distance from a true-ecliptic-of-date Moon state.
///
/// The node direction is `n = k x h = (-h_y, h_x, 0)`, so the longitude
/// is `atan2(h_x, -h_y)` wrapped to `[0, 360)`. The distance comes from the
/// osculating elements (`p = |h|^2/mu`, eccentricity vector `e = (v x
/// h)/mu - r/|r|`, argument of perigee from the node direction,
/// `r = p / (1 + e cos(2pi - omega))`); a vanishing node magnitude or
/// a circular orbit (`ecc <= 1e-15`) falls back to the conventional
/// mean distance, exactly like the model.
pub fn true_lunar_node_from_state(
    r: (f64, f64, f64),
    v: (f64, f64, f64),
) -> Result<(f64, f64, f64), Error> {
    check_ecliptic_state(r, v)?;
    let (h_x, h_y, h_z) = angular_momentum(r, v);
    let node_lon = h_x.atan2(-h_y).to_degrees().rem_euclid(360.0);

    let (r_x, r_y, r_z) = r;
    let (v_x, v_y, v_z) = v;
    let r_mag = (r_x * r_x + r_y * r_y + r_z * r_z).sqrt();
    let h_mag = (h_x * h_x + h_y * h_y + h_z * h_z).sqrt();
    let p_orbit = h_mag * h_mag / GM_EARTH_MOON_AU3_DAY2;

    let vxh_x = v_y * h_z - v_z * h_y;
    let vxh_y = v_z * h_x - v_x * h_z;
    let vxh_z = v_x * h_y - v_y * h_x;
    let e_x = vxh_x / GM_EARTH_MOON_AU3_DAY2 - r_x / r_mag;
    let e_y = vxh_y / GM_EARTH_MOON_AU3_DAY2 - r_y / r_mag;
    let e_z = vxh_z / GM_EARTH_MOON_AU3_DAY2 - r_z / r_mag;
    let ecc = (e_x * e_x + e_y * e_y + e_z * e_z).sqrt();

    let n_mag = (h_y * h_y + h_x * h_x).sqrt();
    let dist = if n_mag > 0.0 && ecc > CIRCULAR_ECCENTRICITY_FLOOR {
        let mut cos_omega = (-h_y * e_x + h_x * e_y) / (n_mag * ecc);
        cos_omega = cos_omega.clamp(-1.0, 1.0);
        let mut omega = cos_omega.acos();
        if e_z < 0.0 {
            omega = 2.0 * std::f64::consts::PI - omega;
        }
        let f_asc = 2.0 * std::f64::consts::PI - omega;
        p_orbit / (1.0 + ecc * f_asc.cos())
    } else {
        MEAN_NODE_DISTANCE_AU
    };
    Ok((node_lon, 0.0, dist))
}

/// True (osculating) lunar-apogee longitude, latitude and distance
/// from a true-ecliptic-of-date Moon state.
///
/// The eccentricity vector points at perigee, so the apogee direction
/// is `-e`; the longitude/latitude are its spherical angles of date
/// and the distance is `p / (1 - e)` with `p = |h|^2/mu`.
pub fn true_lilith_from_state(
    r: (f64, f64, f64),
    v: (f64, f64, f64),
) -> Result<(f64, f64, f64), Error> {
    check_ecliptic_state(r, v)?;
    let (e_x, e_y, e_z, e_mag, p) = apse_elements(r, v);
    // The apogee direction is exactly `-e` (negation is exact in
    // binary64); the model converts that negated vector to spherical
    // angles and divides by the un-negated eccentricity magnitude.
    let (apogee_x, apogee_y, apogee_z) = (-e_x, -e_y, -e_z);
    let apogee_mag = (apogee_x * apogee_x + apogee_y * apogee_y + apogee_z * apogee_z).sqrt();
    let longitude = apogee_y.atan2(apogee_x).to_degrees().rem_euclid(360.0);
    let latitude = (apogee_z / apogee_mag).asin().to_degrees();
    let r_apogee = p / (1.0 - e_mag);
    Ok((longitude, latitude, r_apogee))
}

/// Shared apse-element core in the model component order: angular
/// momentum, eccentricity vector, its magnitude and the semi-latus
/// rectum. Returns `(e_x, e_y, e_z, e_mag, p)`.
fn apse_elements(r: (f64, f64, f64), v: (f64, f64, f64)) -> (f64, f64, f64, f64, f64) {
    let r_mag = (r.0 * r.0 + r.1 * r.1 + r.2 * r.2).sqrt();
    let (h_x, h_y, h_z) = angular_momentum(r, v);
    let mu = GM_EARTH_MOON_AU3_DAY2;
    let vxh_x = v.1 * h_z - v.2 * h_y;
    let vxh_y = v.2 * h_x - v.0 * h_z;
    let vxh_z = v.0 * h_y - v.1 * h_x;
    let e_x = vxh_x / mu - r.0 / r_mag;
    let e_y = vxh_y / mu - r.1 / r_mag;
    let e_z = vxh_z / mu - r.2 / r_mag;
    let e_mag = (e_x * e_x + e_y * e_y + e_z * e_z).sqrt();
    let h_mag = (h_x * h_x + h_y * h_y + h_z * h_z).sqrt();
    let p = h_mag * h_mag / mu;
    (e_x, e_y, e_z, e_mag, p)
}

/// True lunar apogee via the classical orbital-elements spelling
pub fn true_lilith_orbital_elements_from_state(
    r: (f64, f64, f64),
    v: (f64, f64, f64),
) -> Result<(f64, f64, f64), Error> {
    true_lilith_from_state(r, v)
}

/// Osculating lunar-perigee longitude, latitude and distance from a
/// true-ecliptic-of-date Moon state.
///
/// Same eccentricity-vector core as the apogee, but the perigee
/// direction is `+e` itself and the distance is `p / (1 + e)`.
pub fn osculating_perigee_from_state(
    r: (f64, f64, f64),
    v: (f64, f64, f64),
) -> Result<(f64, f64, f64), Error> {
    check_ecliptic_state(r, v)?;
    let (e_x, e_y, e_z, e_mag, p) = apse_elements(r, v);
    let longitude = e_y.atan2(e_x).to_degrees().rem_euclid(360.0);
    let latitude = (e_z / e_mag).asin().to_degrees();
    let r_perigee = p / (1.0 + e_mag);
    Ok((longitude, latitude, r_perigee))
}

/// Wrap-aware central-difference polar rates over two side samples
///
/// The longitude difference folds through the 0/360 seam with the
/// same strict `> 180` / `< -180` (`elif`) branches; latitude and
/// distance difference plainly. The span is `2 * half_step_days`.
/// The float flow is unshaped: no range or finiteness shaping beyond the
/// arithmetic itself.
pub fn central_difference_polar_rates(
    before: (f64, f64, f64),
    after: (f64, f64, f64),
    half_step_days: f64,
) -> (f64, f64, f64) {
    let mut lon_diff = after.0 - before.0;
    if lon_diff > 180.0 {
        lon_diff -= 360.0;
    } else if lon_diff < -180.0 {
        lon_diff += 360.0;
    }
    let span = 2.0 * half_step_days;
    (
        lon_diff / span,
        (after.1 - before.1) / span,
        (after.2 - before.2) / span,
    )
}

/// Central-difference true-point speeds over a position provider
///
/// Samples the provider at `jd_tt - half_step_days` and `jd_tt +
/// half_step_days` and folds them through
/// [`central_difference_polar_rates`]. A provider failure on either
/// side yields `(0.0, 0.0, 0.0)`; provider misses
/// propagate no further. A non-finite epoch is a typed refusal: the
/// the model would fail inside its timescale construction instead.
pub fn stencil_true_point_speed(
    provider: &dyn Fn(f64) -> Result<(f64, f64, f64), Error>,
    jd_tt: f64,
    half_step_days: f64,
) -> Result<(f64, f64, f64), Error> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    if !half_step_days.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: half_step_days,
        });
    }
    let before = provider(jd_tt - half_step_days);
    let after = provider(jd_tt + half_step_days);
    match (before, after) {
        (Ok(before), Ok(after)) => Ok(central_difference_polar_rates(
            before,
            after,
            half_step_days,
        )),
        _ => Ok((0.0, 0.0, 0.0)),
    }
}

// Osculating lunar geometry from a shared geocentric state.
//
// The Moon-minus-Earth position and velocity are rotated from ICRS to the
// requested ecliptic frame before the conic reduction. TT epochs, AU distances
// and AU/day velocities are supplied by the bound native state provider.
// Acquisition and frame errors propagate without changing the requested data.

use crate::lunar_apse_tables::{
    APOGEE_CORRECTIONS, APOGEE_CT_COUNT, APOGEE_CT_JD_START, APOGEE_CT_STEP_DAYS,
    APOGEE_LAT_COEFFS, APOGEE_MEAN_DIST_AU, APOGEE_TRIG_TERMS, PERIGEE_CORRECTIONS,
    PERIGEE_CT_COUNT, PERIGEE_CT_JD_START, PERIGEE_CT_STEP_DAYS, PERIGEE_LAT_COEFFS,
    PERIGEE_MEAN_DIST_AU, PERIGEE_TRIG_TERMS,
};

/// Linear taper window applied to correction-table edge values outside
/// the table range, in days.
/// The taper keeps positions continuous at the table limits and decays the edge correction to zero within a year of the edge.
pub const APSE_EDGE_TAPER_DAYS: f64 = 365.25;

/// Centered half-step of the interpolated-apsis speed stencil, in days
pub const INTERP_APSE_SPEED_HALF_STEP_DAYS: f64 = 0.002;

/// Evaluate one DE440 passage-fitted trig series at `jd_tt`.
///
/// The basis is `coeff * E**e_power * sin(cD*D + cM*M + cMp*Mp + cF*F)`
/// over the IERS Delaunay arguments `(D, M, M', F)` in radians, with
/// `E` the normalized mean Earth-orbit eccentricity; the summation and
/// the per-term `e_power` branches are the calculation order.
fn de440_passage_terms(jd_tt: f64, terms: &[(f64, i32, i32, i32, i32, i32)]) -> f64 {
    let t = centuries_since_j2000(jd_tt);
    de440_passage_terms_acquired(t, jd_tt, terms)
}

/// Evaluate the fitted series after independent scalar epoch acquisitions.
fn de440_passage_terms_acquired(
    t: f64,
    jd_tt: f64,
    terms: &[(f64, i32, i32, i32, i32, i32)],
) -> f64 {
    de440_passage_terms_from_arguments(t, t * t, lunar_delaunay_arguments(jd_tt), terms)
}

/// Fitted-series arithmetic after argument and square primitives complete.
fn de440_passage_terms_from_arguments(
    t: f64,
    t_squared: f64,
    arguments: (f64, f64, f64, f64),
    terms: &[(f64, i32, i32, i32, i32, i32)],
) -> f64 {
    let (delaunay_d, m_sun, m_moon, arg_f) = arguments;
    let eccentricity =
        EARTH_ECCENTRICITY_J2000 + EARTH_ECCENTRICITY_T * t + EARTH_ECCENTRICITY_T2 * t_squared;
    let e_norm = eccentricity / EARTH_ECCENTRICITY_J2000;
    let mut perturbation = 0.0;
    for (coeff, c_d, c_m, c_mp, c_f, e_pow) in terms {
        let angle = (*c_d as f64) * delaunay_d
            + (*c_m as f64) * m_sun
            + (*c_mp as f64) * m_moon
            + (*c_f as f64) * arg_f;
        let mut term = *coeff * angle.sin();
        if *e_pow == 1 {
            term *= e_norm;
        } else if *e_pow == 2 {
            term *= e_norm * e_norm;
        }
        perturbation += term;
    }
    perturbation
}

/// Apogee passage-fitted trig perturbation in degrees.
/// Trig-only residuals against the DE440 apogee passages have stated RMS error; the residual table supplies the remaining low-frequency correction.
pub fn de440_apogee_passage_terms(jd_tt: f64) -> f64 {
    de440_passage_terms(jd_tt, APOGEE_TRIG_TERMS)
}

/// Perigee passage-fitted trig perturbation in degrees.
/// The perigee basis is generated independently of the apogee basis; its trig-only residuals have their own stated RMS error, likewise completed by the residual table.
pub fn de440_perigee_passage_terms(jd_tt: f64) -> f64 {
    de440_passage_terms(jd_tt, PERIGEE_TRIG_TERMS)
}

/// Fitted apsis terms in degrees after acquiring their separate input clocks.
///
/// `centuries` is `(raw_jd - J2000) / 36525`, while `delaunay_jd_tt` is
/// the separately acquired scalar TT used by the IERS arguments. Apogee and
/// perigee select independent fixed coefficient tables. No epoch conversion
/// or finite check is repeated; IEEE non-finite arithmetic propagates.
pub fn interpolated_apse_passage_terms(centuries: f64, delaunay_jd_tt: f64, apogee: bool) -> f64 {
    let terms = if apogee {
        APOGEE_TRIG_TERMS
    } else {
        PERIGEE_TRIG_TERMS
    };
    de440_passage_terms_acquired(centuries, delaunay_jd_tt, terms)
}

/// Passage series in degrees from separately acquired centuries and arguments.
///
/// `centuries_squared` is the result of the caller's actual power primitive,
/// acquired after the IERS Delaunay arguments `(D, M, M', F)` in radians.
/// This permits a language boundary to preserve power overflow and operator
/// order while the fixed harmonic series and eccentricity arithmetic remain
/// in this kernel. Non-finite values propagate without a finite guard.
pub fn interpolated_apse_passage_from_arguments(
    centuries: f64,
    centuries_squared: f64,
    arguments: (f64, f64, f64, f64),
    apogee: bool,
) -> f64 {
    let terms = if apogee {
        APOGEE_TRIG_TERMS
    } else {
        PERIGEE_TRIG_TERMS
    };
    de440_passage_terms_from_arguments(centuries, centuries_squared, arguments, terms)
}

/// Fixed residual-table origin and step, in TT Julian days and days.
pub fn interpolated_apse_table_grid(apogee: bool) -> (f64, f64) {
    if apogee {
        (APOGEE_CT_JD_START, APOGEE_CT_STEP_DAYS)
    } else {
        (PERIGEE_CT_JD_START, PERIGEE_CT_STEP_DAYS)
    }
}

/// Arcsecond residual after the caller acquires the raw table-index quotient.
///
/// Infinite indices reach the tapered edge. NaN refuses only at the interior
/// integer-index stage. Tables, counts and positive steps are fixed owned
/// generated data; callers cannot introduce out-of-bounds table metadata.
pub fn interpolated_apse_correction_from_index(index: f64, apogee: bool) -> Result<f64, Error> {
    let (_, step) = interpolated_apse_table_grid(apogee);
    let (corrections, count) = if apogee {
        (APOGEE_CORRECTIONS, APOGEE_CT_COUNT)
    } else {
        (PERIGEE_CORRECTIONS, PERIGEE_CT_COUNT)
    };
    interpolate_apse_correction_index(index, step, corrections, count)
}

/// Longitude in degrees from the mean, passage degrees and arcsecond residual.
/// Addition and arcsecond conversion precede normalization to `[0, 360)`.
pub fn interpolated_apse_longitude_from_parts(mean: f64, passage: f64, residual: f64) -> f64 {
    ((mean + passage) + residual / 3600.0).rem_euclid(360.0)
}

/// Latitude in degrees and conventional distance in AU from acquired angles.
/// `node` is the separately acquired mean-node longitude in degrees. The
/// fixed two-harmonic latitude fit and mean distance are selected per apsis.
pub fn interpolated_apse_latitude_distance(lon: f64, node: f64, apogee: bool) -> (f64, f64) {
    let (coeffs, distance) = if apogee {
        (APOGEE_LAT_COEFFS, APOGEE_MEAN_DIST_AU)
    } else {
        (PERIGEE_LAT_COEFFS, PERIGEE_MEAN_DIST_AU)
    };
    (interpolated_apse_latitude(lon, node, coeffs), distance)
}

/// Interpolate an apse residual correction in arcseconds.
///
/// Linear interpolation between JD-indexed table entries. Outside the
/// table range the edge value tapers linearly to 0.0 over
/// [`APSE_EDGE_TAPER_DAYS`] (a hard step here used to jump positions by
/// the edge correction, up to a few hundred arcsec at the table
/// limits). `n` mirrors the model table-length argument; a NaN epoch
/// is refused.
pub fn interpolate_apse_correction(
    jd_tt: f64,
    jd_start: f64,
    step_days: f64,
    corrections: &[f64],
    n: usize,
) -> Result<f64, Error> {
    if corrections.is_empty() || n == 0 {
        return Ok(0.0);
    }
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    interpolate_apse_correction_source(jd_tt, jd_start, step_days, corrections, n)
}

/// Residual-table acquisition without a public finite-epoch guard.
///
/// Infinite epochs reach the appropriate tapered edge and return zero. NaN
/// reaches the interior integer-index stage after both comparisons and is
/// refused there, allowing a binding to acquire its native `int(NaN)` error.
/// The table must contain at least `n` entries and the step must be positive;
/// the owning generated apsis tables establish those private invariants.
fn interpolate_apse_correction_source(
    jd_tt: f64,
    jd_start: f64,
    step_days: f64,
    corrections: &[f64],
    n: usize,
) -> Result<f64, Error> {
    if corrections.is_empty() || n == 0 {
        return Ok(0.0);
    }
    let idx_float = (jd_tt - jd_start) / step_days;
    interpolate_apse_correction_index(idx_float, step_days, corrections, n)
}

/// Fixed-table correction after raw subtraction/division has completed.
fn interpolate_apse_correction_index(
    idx_float: f64,
    step_days: f64,
    corrections: &[f64],
    n: usize,
) -> Result<f64, Error> {
    if idx_float < 0.0 {
        let dist_days = -idx_float * step_days;
        let weight = 1.0 - dist_days / APSE_EDGE_TAPER_DAYS;
        if weight <= 0.0 {
            return Ok(0.0);
        }
        return Ok(corrections[0] * weight);
    }
    if idx_float >= n as f64 - 1.0 {
        let dist_days = (idx_float - (n as f64 - 1.0)) * step_days;
        let weight = 1.0 - dist_days / APSE_EDGE_TAPER_DAYS;
        if weight <= 0.0 {
            return Ok(0.0);
        }
        return Ok(corrections[n - 1] * weight);
    }
    if idx_float.is_nan() {
        return Err(Error::NonFiniteArgument { value: idx_float });
    }
    let idx_low = idx_float as usize;
    let frac = idx_float - idx_low as f64;
    Ok(corrections[idx_low] + frac * (corrections[idx_low + 1] - corrections[idx_low]))
}

/// Shared longitude core: mean baseline + trig perturbation + residual
/// table.
fn interpolated_apse_longitude(
    jd_tt: f64,
    mean_baseline_lon: f64,
    passage_terms: f64,
    jd_start: f64,
    step_days: f64,
    corrections: &[f64],
    n: usize,
) -> Result<f64, Error> {
    let mut interp_lon = mean_baseline_lon + passage_terms;
    interp_lon +=
        interpolate_apse_correction_source(jd_tt, jd_start, step_days, corrections, n)? / 3600.0;
    Ok(interp_lon.rem_euclid(360.0))
}

/// Shared latitude core: two-harmonic inclined-plane fit
/// `a1*sin(w) + a3*sin(3w)` with `w = lon - node` in radians.
fn interpolated_apse_latitude(interp_lon: f64, node_lon: f64, lat_coeffs: (f64, f64)) -> f64 {
    let omega = (interp_lon - node_lon).to_radians();
    lat_coeffs.0 * omega.sin() + lat_coeffs.1 * (3.0 * omega).sin()
}

/// Independently generated interpolated lunar-apogee position
/// `(lon, lat, dist)` in degrees and AU on the true ecliptic of date
///
/// Beyond the fitted table interval the residual tapers out and the trig-only track carries stated RMS error; this abstract coordinate is not a lunar-distance event finder.
pub fn calc_interpolated_apogee(jd_tt: f64) -> Result<(f64, f64, f64), Error> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    calc_interpolated_apogee_source(jd_tt)
}

/// Interpolated apogee from an acquired scalar TT Julian day.
///
/// Degrees/AU, baseline, passage harmonics and tapered correction table follow
/// [`calc_interpolated_apogee`]. There is no entry finite guard: infinities
/// propagate to NaN angles with the conventional finite distance; NaN fails
/// at the table's integer-index primitive after analytical terms are consumed.
pub fn calc_interpolated_apogee_source(jd_tt: f64) -> Result<(f64, f64, f64), Error> {
    let mean_apogee = mean_lunar_apogee_position_unchecked(jd_tt).0;
    let interp_lon = interpolated_apse_longitude(
        jd_tt,
        mean_apogee,
        de440_apogee_passage_terms(jd_tt),
        APOGEE_CT_JD_START,
        APOGEE_CT_STEP_DAYS,
        APOGEE_CORRECTIONS,
        APOGEE_CT_COUNT,
    )?;
    let node_lon = mean_lunar_node_position_unchecked(jd_tt).0;
    Ok((
        interp_lon,
        interpolated_apse_latitude(interp_lon, node_lon, APOGEE_LAT_COEFFS),
        APOGEE_MEAN_DIST_AU,
    ))
}

/// Independently generated interpolated lunar-perigee position
/// `(lon, lat, dist)` in degrees and AU on the true ecliptic of date
///
/// Beyond the fitted table interval the residual tapers out and the trig-only track carries stated RMS error; this abstract coordinate is not a lunar-distance or tidal event finder.
pub fn calc_interpolated_perigee(jd_tt: f64) -> Result<(f64, f64, f64), Error> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    calc_interpolated_perigee_source(jd_tt)
}

/// Interpolated perigee from an acquired scalar TT Julian day.
///
/// Units and independent fit follow [`calc_interpolated_perigee`]. As for
/// [`calc_interpolated_apogee_source`], infinities propagate through the
/// analytical fit and NaN is refused only at the correction-table index stage.
pub fn calc_interpolated_perigee_source(jd_tt: f64) -> Result<(f64, f64, f64), Error> {
    let mean_perigee = (mean_lunar_apogee_position_unchecked(jd_tt).0 + 180.0).rem_euclid(360.0);
    let interp_lon = interpolated_apse_longitude(
        jd_tt,
        mean_perigee,
        de440_perigee_passage_terms(jd_tt),
        PERIGEE_CT_JD_START,
        PERIGEE_CT_STEP_DAYS,
        PERIGEE_CORRECTIONS,
        PERIGEE_CT_COUNT,
    )?;
    let node_lon = mean_lunar_node_position_unchecked(jd_tt).0;
    Ok((
        interp_lon,
        interpolated_apse_latitude(interp_lon, node_lon, PERIGEE_LAT_COEFFS),
        PERIGEE_MEAN_DIST_AU,
    ))
}

/// Centered polar-state stencil over a position function.
///
/// Samples `position` at `jd_tt - half_step_days` and `jd_tt +
/// half_step_days`; the longitude difference folds through the 0/360
/// seam with the `(delta + 180) % 360 - 180` rule, latitude
/// and distance difference plainly, over the `2 * half_step_days`
/// span.
/// The 0.002-day-scale half-step matches the osculating-apogee stencil: short-period structure needs it, while the derivative is stable across a wide plateau around it.
pub fn finite_difference_polar_state(
    position: &dyn Fn(f64) -> Result<(f64, f64, f64), Error>,
    jd_tt: f64,
    half_step_days: f64,
) -> Result<(f64, f64, f64, f64, f64, f64), Error> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    if !half_step_days.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: half_step_days,
        });
    }
    let (longitude, latitude, distance) = position(jd_tt)?;
    let before = position(jd_tt - half_step_days)?;
    let after = position(jd_tt + half_step_days)?;
    let longitude_delta = (after.0 - before.0 + 180.0).rem_euclid(360.0) - 180.0;
    let span = 2.0 * half_step_days;
    Ok((
        longitude,
        latitude,
        distance,
        longitude_delta / span,
        (after.1 - before.1) / span,
        (after.2 - before.2) / span,
    ))
}

/// Independently generated interpolated lunar-apsis polar state
///
/// `body_id` is [`INTP_APOG`] or
/// [`INTP_PERG`]; anything else is
/// [`Error::UnknownBody`]. `speed3` is accepted and ignored: the
/// stencil derivative is the smooth instantaneous rate the model
/// serves for both SPEED and SPEED 3.
pub fn calc_interpolated_apse_state(
    jd_tt: f64,
    body_id: i32,
    speed3: bool,
) -> Result<(f64, f64, f64, f64, f64, f64), Error> {
    let _ = speed3;
    if body_id == INTP_APOG {
        finite_difference_polar_state(
            &calc_interpolated_apogee,
            jd_tt,
            INTERP_APSE_SPEED_HALF_STEP_DAYS,
        )
    } else if body_id == INTP_PERG {
        finite_difference_polar_state(
            &calc_interpolated_perigee,
            jd_tt,
            INTERP_APSE_SPEED_HALF_STEP_DAYS,
        )
    } else {
        Err(Error::UnknownBody { body_id })
    }
}

// Osculating lunar geometry from a shared geocentric state.
//
// The Moon-minus-Earth position and velocity are rotated from ICRS to the
// requested ecliptic frame before the conic reduction. TT epochs, AU distances
// and AU/day velocities are supplied by the bound native state provider.
// Acquisition and frame errors propagate without changing the requested data.

/// One Delaunay sine-term descriptor
/// `(cD, cM, cMp, cF, e_power)` for
/// `E**e_power * sin(cD*D + cM*M + cMp*Mp + cF*F)`.
pub type ApseBasisTerm = (i32, i32, i32, i32, i32);

/// Construct the project-selected Delaunay sine-term inventory
/// : the `1..=max_k` principal evection ladder, the
/// three first-order solar terms, the `sideband_ks` first-order solar
/// sideband pairs, the six second-order terms, and the `2(F - Mp)`
/// inclination term, in the calculation order.
///
/// The ladder and its truncation are project model-design choices
/// validated against physical DE440 apsis passages.
pub fn apse_basis_terms(max_k: u32, sideband_ks: &[u32]) -> Vec<ApseBasisTerm> {
    let mut terms = Vec::new();
    for k in 1..=max_k {
        let k = k as i32;
        terms.push((k, 0, -k, 0, 0));
    }
    terms.push((0, 1, 0, 0, 1));
    terms.push((1, 1, -1, 0, 1));
    terms.push((1, -1, -1, 0, 1));
    for &k in sideband_ks {
        let k = k as i32;
        terms.push((k, -1, -k, 0, 1));
        terms.push((k, 1, -k, 0, 1));
    }
    terms.push((0, 2, 0, 0, 2));
    terms.push((2, -2, -2, 0, 2));
    terms.push((2, 2, -2, 0, 2));
    terms.push((4, -2, -4, 0, 2));
    terms.push((4, 2, -4, 0, 2));
    terms.push((6, -2, -6, 0, 2));
    terms.push((0, 0, -2, 2, 0));
    terms
}

/// `APOGEE_TERMS` inventory (`_basis(14, [2, 4, 6, 8, 10])`):
/// 34 descriptors.
pub fn apogee_basis_terms() -> Vec<ApseBasisTerm> {
    apse_basis_terms(14, &[2, 4, 6, 8, 10])
}

/// `PERIGEE_TERMS` inventory (`_basis(30, 2..15)`):
/// 66 descriptors.
pub fn perigee_basis_terms() -> Vec<ApseBasisTerm> {
    apse_basis_terms(30, &[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
}

/// One generator design-matrix row over explicit Delaunay arguments
/// : `(E**e) * sin(arg)` with
/// `arg = cD*D + cM*M + cMp*Mp + cF*F` in radians and `e_norm` the
/// Simon mean Earth-orbit eccentricity normalized to J2000. Columns
/// follow `terms` exactly. Non-finite inputs are typed
/// [`Error::NonFiniteArgument`] refusals.
///
/// This is a linear-model design row only; it carries no observations
/// from any external ephemeris product.
pub fn apse_design_row(
    delaunay_d: f64,
    m_sun: f64,
    m_moon: f64,
    arg_f: f64,
    e_norm: f64,
    terms: &[ApseBasisTerm],
) -> Result<Vec<f64>, Error> {
    for value in [delaunay_d, m_sun, m_moon, arg_f, e_norm] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    let mut row = Vec::with_capacity(terms.len());
    for (c_d, c_m, c_mp, c_f, e_pow) in terms {
        let angle = (*c_d as f64) * delaunay_d
            + (*c_m as f64) * m_sun
            + (*c_mp as f64) * m_moon
            + (*c_f as f64) * arg_f;
        row.push(e_norm.powi(*e_pow) * angle.sin());
    }
    Ok(row)
}

/// Natural cubic spline with the model construction: zero endpoint second derivatives, Thomas
/// solve for the interior coefficients, per-interval `(a, b, c, d)`
/// storage. Fields stay private with accessors.
pub struct NaturalCubicSpline {
    knots: Vec<f64>,
    a: Vec<f64>,
    b: Vec<f64>,
    c: Vec<f64>,
    d: Vec<f64>,
}

impl NaturalCubicSpline {
    /// Build the spline over strictly increasing finite `knots` with
    /// finite `values`.
    /// Shape/order violations are [`Error::InvalidGeneratorInput`]
    ///
    /// non-finite knots or values are [`Error::NonFiniteArgument`].
    pub fn try_new(knots: &[f64], values: &[f64]) -> Result<Self, Error> {
        if knots.len() != values.len() || knots.len() < 2 {
            return Err(Error::InvalidGeneratorInput {
                entry: "NaturalCubicSpline::try_new",
                detail: "knot and value slices must share a length >= 2",
            });
        }
        for &value in knots.iter().chain(values.iter()) {
            if !value.is_finite() {
                return Err(Error::NonFiniteArgument { value });
            }
        }
        for window in knots.windows(2) {
            if window[1] <= window[0] {
                return Err(Error::InvalidGeneratorInput {
                    entry: "NaturalCubicSpline::try_new",
                    detail: "knots must be strictly increasing",
                });
            }
        }
        let n = knots.len() - 1;
        let mut step = Vec::with_capacity(n);
        for window in knots.windows(2) {
            step.push(window[1] - window[0]);
        }
        // `alpha[0]` and `alpha[n]` stay 0.0 and are never read.
        let mut alpha = vec![0.0; n + 1];
        for i in 1..n {
            alpha[i] = 3.0
                * ((values[i + 1] - values[i]) / step[i]
                    - (values[i] - values[i - 1]) / step[i - 1]);
        }
        let mut line = vec![1.0; n + 1];
        let mut mu = vec![0.0; n + 1];
        let mut z = vec![0.0; n + 1];
        for i in 1..n {
            line[i] = 2.0 * (knots[i + 1] - knots[i - 1]) - step[i - 1] * mu[i - 1];
            mu[i] = step[i] / line[i];
            z[i] = (alpha[i] - step[i - 1] * z[i - 1]) / line[i];
        }
        let mut c_full = vec![0.0; n + 1];
        let mut b = vec![0.0; n];
        let mut d = vec![0.0; n];
        for j in (0..n).rev() {
            c_full[j] = z[j] - mu[j] * c_full[j + 1];
            b[j] = (values[j + 1] - values[j]) / step[j]
                - step[j] * (c_full[j + 1] + 2.0 * c_full[j]) / 3.0;
            d[j] = (c_full[j + 1] - c_full[j]) / (3.0 * step[j]);
        }
        Ok(Self {
            knots: knots.to_vec(),
            a: values[..n].to_vec(),
            b,
            c: c_full[..n].to_vec(),
            d,
        })
    }

    /// Evaluate at one query point:
    /// left-sided knot search (`np.searchsorted` side `'left'`) minus
    /// one, clamped to the interval range, then the cubic Horner form.
    /// Queries outside the knots extrapolate the edge interval on both
    /// sides. A non-finite query is [`Error::NonFiniteArgument`].
    pub fn eval(&self, xq: f64) -> Result<f64, Error> {
        if !xq.is_finite() {
            return Err(Error::NonFiniteArgument { value: xq });
        }
        let found = self.knots.partition_point(|&knot| knot < xq);
        let idx = found.saturating_sub(1).min(self.b.len() - 1);
        let dx = xq - self.knots[idx];
        Ok(self.a[idx] + dx * (self.b[idx] + dx * (self.c[idx] + dx * self.d[idx])))
    }

    /// Evaluate over several query points sharing one construction.
    pub fn eval_many(&self, queries: &[f64]) -> Result<Vec<f64>, Error> {
        queries.iter().map(|&xq| self.eval(xq)).collect()
    }

    /// Knot abscissae (strictly increasing).
    pub fn knots(&self) -> &[f64] {
        &self.knots
    }

    /// Interval constant terms (`values[..n]`).
    pub fn a(&self) -> &[f64] {
        &self.a
    }

    /// Interval linear coefficients.
    pub fn b(&self) -> &[f64] {
        &self.b
    }

    /// Interval quadratic coefficients.
    pub fn c(&self) -> &[f64] {
        &self.c
    }

    /// Interval cubic coefficients.
    pub fn d(&self) -> &[f64] {
        &self.d
    }
}

/// Assembled residual correction grid with the model table-error
/// diagnostics (post-fit half of `build_side`).
pub struct ApseGridAssembly {
    corrections_arcsec: Vec<f64>,
    table_error_deg: Vec<f64>,
    table_rms_arcsec: f64,
    table_max_arcsec: f64,
}

impl ApseGridAssembly {
    /// Grid corrections in arcseconds (length `ct_count`).
    pub fn corrections_arcsec(&self) -> &[f64] {
        &self.corrections_arcsec
    }

    /// `residual - interpolated` at each passage, in degrees.
    pub fn table_error_deg(&self) -> &[f64] {
        &self.table_error_deg
    }

    /// RMS table error in arcseconds (diagnostic).
    pub fn table_rms_arcsec(&self) -> f64 {
        self.table_rms_arcsec
    }

    /// Maximum absolute table error in arcseconds (diagnostic).
    pub fn table_max_arcsec(&self) -> f64 {
        self.table_max_arcsec
    }
}

/// Assemble a uniform lunar-apse residual grid from passage samples.
///
/// Inside the passage span, spline residuals are converted to arcseconds.
/// Outside it, the terminal residual tapers linearly to zero across the stated
/// edge interval. Grid epochs are ct_start+ct_step_days*i. RMS and maximum
/// residuals describe this construction, not a complete orbit-error bound.
pub fn assemble_correction_grid(
    passages: &[f64],
    residuals_deg: &[f64],
    ct_start: f64,
    ct_step_days: f64,
    ct_count: usize,
    edge_taper_days: f64,
) -> Result<ApseGridAssembly, Error> {
    if passages.len() != residuals_deg.len() || passages.len() < 2 {
        return Err(Error::InvalidGeneratorInput {
            entry: "assemble_correction_grid",
            detail: "passage and residual slices must share a length >= 2",
        });
    }
    if ct_count == 0 {
        return Err(Error::InvalidGeneratorInput {
            entry: "assemble_correction_grid",
            detail: "correction grid must hold >= 1 entry",
        });
    }
    for value in [ct_start, ct_step_days, edge_taper_days] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    if ct_step_days <= 0.0 || edge_taper_days <= 0.0 {
        return Err(Error::InvalidGeneratorInput {
            entry: "assemble_correction_grid",
            detail: "grid step and edge taper must be finite and positive",
        });
    }
    // Non-finite passages/residuals are refused here with
    // `NonFiniteArgument`.
    let spline = NaturalCubicSpline::try_new(passages, residuals_deg)?;
    let lo = passages[0];
    let hi = passages[passages.len() - 1];
    let mut corrections = Vec::with_capacity(ct_count);
    for i in 0..ct_count {
        let grid_jd = ct_start + ct_step_days * (i as f64);
        if grid_jd >= lo && grid_jd <= hi {
            corrections.push(spline.eval(grid_jd)? * 3600.0);
        } else if grid_jd < lo {
            let weight = (1.0 - (lo - grid_jd) / edge_taper_days).clamp(0.0, 1.0);
            corrections.push(residuals_deg[0] * weight * 3600.0);
        } else {
            let weight = (1.0 - (grid_jd - hi) / edge_taper_days).clamp(0.0, 1.0);
            corrections.push(residuals_deg[residuals_deg.len() - 1] * weight * 3600.0);
        }
    }
    let mut table_error = Vec::with_capacity(passages.len());
    for (k, &jd) in passages.iter().enumerate() {
        let interpolated = interp_uniform_grid(ct_start, ct_step_days, &corrections, jd);
        table_error.push(residuals_deg[k] - interpolated / 3600.0);
    }
    let mut sum_sq: f64 = 0.0;
    let mut max_abs: f64 = 0.0;
    for &error in &table_error {
        sum_sq += error * error;
        max_abs = max_abs.max(error.abs());
    }
    Ok(ApseGridAssembly {
        corrections_arcsec: corrections,
        table_error_deg: table_error,
        table_rms_arcsec: (sum_sq / passages.len() as f64).sqrt() * 3600.0,
        table_max_arcsec: max_abs * 3600.0,
    })
}

/// Piecewise-linear interpolation over a uniform grid with endpoint
/// fill:
/// queries at or outside the grid ends return the edge values, interior
/// queries interpolate the bracketing segment. At an exact interior
/// knot both bracketing segments agree up to the final rounding, so the
/// the final segment is selected by the actual knot comparison.
fn interp_uniform_grid(ct_start: f64, ct_step_days: f64, values: &[f64], xq: f64) -> f64 {
    let n = values.len();
    let last = ct_start + ct_step_days * ((n - 1) as f64);
    if xq <= ct_start {
        return values[0];
    }
    if xq >= last {
        return values[n - 1];
    }
    let pos = (xq - ct_start) / ct_step_days;
    let i = (pos.floor() as usize).min(n - 2);
    let x0 = ct_start + ct_step_days * (i as f64);
    values[i] + (values[i + 1] - values[i]) * ((xq - x0) / ct_step_days)
}

/// Split one distance sample series into apogee/perigee candidate index
/// sets: maxima satisfy `r[i] >= r[i-1] && r[i] > r[i+1]`, minima
/// the mirror. Series shorter than three samples hold no inner point
/// and yield two empty sets. NaN samples never satisfy a predicate on
/// either side, so they are skipped exactly like the model
/// comparisons skip them.
pub fn scan_passage_candidates(distances: &[f64]) -> (Vec<usize>, Vec<usize>) {
    let mut maxima = Vec::new();
    let mut minima = Vec::new();
    if distances.len() < 3 {
        return (maxima, minima);
    }
    for i in 1..distances.len() - 1 {
        let (before, at, after) = (distances[i - 1], distances[i], distances[i + 1]);
        if at >= before && at > after {
            maxima.push(i);
        }
        if at <= before && at < after {
            minima.push(i);
        }
    }
    (maxima, minima)
}

/// Parabolic passage refinement: six rounds of symmetric `±h` sampling around `t0` with the
/// `0.5 * h * (rm - rp) / denom` vertex shift clipped to `[-h, h]` and
/// `h /= 4` per round. `distance_at` is the injected provider; provider failures
/// propagate fail-closed. A zero denominator — or a non-finite one,
/// whose `abs > 0.0` test is false on both sides — leaves `t0` in
/// place for that round, exactly like the model `np.where` arm.
pub fn refine_passage(
    t0: f64,
    step_days: f64,
    distance_at: &dyn Fn(f64) -> Result<f64, Error>,
) -> Result<f64, Error> {
    if !t0.is_finite() {
        return Err(Error::NonFiniteArgument { value: t0 });
    }
    if !step_days.is_finite() || step_days <= 0.0 {
        return Err(Error::InvalidGeneratorInput {
            entry: "refine_passage",
            detail: "step_days must be finite and positive",
        });
    }
    let mut t = t0;
    let mut h = step_days;
    for _ in 0..6 {
        let rm = distance_at(t - h)?;
        let r0 = distance_at(t)?;
        let rp = distance_at(t + h)?;
        let denom = rm - 2.0 * r0 + rp;
        let shift = if denom.abs() > 0.0 {
            0.5 * h * (rm - rp) / denom
        } else {
            0.0
        };
        t += shift.clamp(-h, h);
        h /= 4.0;
    }
    Ok(t)
}

// Lunar-point family dispatch and output frames.
//
// Mean, osculating and interpolated points share native geometry and output
// flag handling. The state provider supplies geocentric lunar states; other
// providers supply nutation, ayanamsha and coordinate transforms. Central-state
// failures propagate. Side-sample state failures produce zero rates under the
// existing boundary policy. Meeus validity ranges are advisory and may emit
// MeeusPolynomialWarning. Undefined conics retain NaN values; heliocentric and
// barycentric lunar points use the documented zero sentinel.

/// Shaped lunar output state: ecliptic `(lon, lat, dist, dlon, dlat,
/// ddist)` in degrees, AU and degrees/day (AU/day radial). Plain alias
/// for the 6-tuples used across the mean/true/interpolated point kernels;
/// it only keeps the dispatch-layer signatures readable.
pub type LunarState6 = (f64, f64, f64, f64, f64, f64);
/// Raw lunar position triple: ecliptic `(lon, lat, dist)` in degrees
/// and AU.
pub type LunarPosition3 = (f64, f64, f64);
/// True-ecliptic-of-date Moon state: geocentric position in AU and
/// velocity in AU/day.
pub type EclipticMoonState = (LunarPosition3, LunarPosition3);

/// True when `body_id` is one of the six served lunar points
/// (`MEAN_NODE`, `TRUE_NODE`, `MEAN_APOG`, `OSCU_APOG`, `INTP_APOG`,
/// `INTP_PERG`).
pub fn is_lunar_point(body_id: i32) -> bool {
    matches!(
        body_id,
        MEAN_NODE | TRUE_NODE | MEAN_APOG | OSCU_APOG | INTP_APOG | INTP_PERG
    )
}

/// True for the analytic mean points (`MEAN_NODE`, `MEAN_APOG`).
pub fn is_mean_lunar_point(body_id: i32) -> bool {
    matches!(body_id, MEAN_NODE | MEAN_APOG)
}

/// True for the JPL-state osculating points (`TRUE_NODE`, `OSCU_APOG`).
pub fn is_true_lunar_point(body_id: i32) -> bool {
    matches!(body_id, TRUE_NODE | OSCU_APOG)
}

/// True for the generated-model interpolated points (`INTP_APOG`,
/// `INTP_PERG`).
pub fn is_interp_lunar_point(body_id: i32) -> bool {
    matches!(body_id, INTP_APOG | INTP_PERG)
}

/// True for the ascending-node points (`MEAN_NODE`, `TRUE_NODE`).
pub fn is_lunar_node(body_id: i32) -> bool {
    matches!(body_id, MEAN_NODE | TRUE_NODE)
}

/// True for the descending ("south") node aliases (`-MEAN_NODE`,
/// `-TRUE_NODE`).
///
/// the model public entries intercept these before `_calc_body` and
/// derive the antipode via `planets._south_node_from_north` (the
/// center-composition layer's `south_node_from_north`); the dispatch below only ever receives the
/// positive (north) id.
pub fn is_south_lunar_node(body_id: i32) -> bool {
    body_id == -MEAN_NODE || body_id == -TRUE_NODE
}

/// Map a south-node alias to its north source id, or `None` for any
/// other body.
pub fn south_node_north_source(body_id: i32) -> Option<i32> {
    if body_id == -MEAN_NODE {
        Some(MEAN_NODE)
    } else if body_id == -TRUE_NODE {
        Some(TRUE_NODE)
    } else {
        None
    }
}

/// Moon-derived abstract points have no heliocentric or barycentric
/// place: under
/// `FLG_HELCTR`/`FLG_BARYCTR` the position is an exact zero 6-tuple
/// with the normal echoed flags.
pub fn wants_lunar_helio_bary_zero(body_id: i32, flags: i32) -> bool {
    is_lunar_point(body_id) && (flags & (FLG_HELCTR | FLG_BARYCTR)) != 0
}

/// Same `_no_nut` rule the model applies to all six lunar points:
/// the output stays on the mean ecliptic (no `dpsi` shaping) under
/// `FLG_NONUT`/`FLG_J2000`, or under `SIDEREAL+EQUATORIAL` (which the
/// the model converts with the mean obliquity, so the longitude must
/// stay mean too).
pub fn no_nutation_frame(flags: i32) -> bool {
    (flags & (FLG_NONUT | FLG_J2000)) != 0
        || ((flags & FLG_SIDEREAL) != 0 && (flags & FLG_EQUATORIAL) != 0)
}

/// Same `_sid_j2000_ecl` rule as the model: `SIDEREAL+J2000` without
/// `EQUATORIAL` composes by expressing the tropical state on the J2000
/// ecliptic first and only then subtracting the mean ayanamsha (the
/// rotations do not commute; this ordering is specific to the lunar
/// points).
pub fn is_sidereal_j2000_ecliptic(flags: i32) -> bool {
    (flags & FLG_SIDEREAL) != 0 && (flags & FLG_J2000) != 0 && (flags & FLG_EQUATORIAL) == 0
}

/// Output-frame application order for a shaped lunar 6-tuple.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LunarFrameOrder {
    /// `SIDEREAL+J2000` ecliptic: equatorial-convert first (J2000
    /// precession), then apply the sidereal correction in that frame.
    EquatorialThenSidereal,
    /// Plain sidereal ecliptic: sidereal-correct first, then convert.
    SiderealThenEquatorial,
    /// Tropical (or sidereal-equatorial): convert only.
    EquatorialOnly,
}

/// Select the model frame-application order for `flags`.
pub fn lunar_frame_order(flags: i32) -> LunarFrameOrder {
    if is_sidereal_j2000_ecliptic(flags) {
        LunarFrameOrder::EquatorialThenSidereal
    } else if (flags & FLG_SIDEREAL) != 0 && (flags & FLG_EQUATORIAL) == 0 {
        LunarFrameOrder::SiderealThenEquatorial
    } else {
        LunarFrameOrder::EquatorialOnly
    }
}

/// Apply the model frame-application order over injected providers
///
/// `equatorial_convert` maps `(state, jd_tt, flags)` through the
/// J2000-precession/equatorial tail; `sidereal_correct` maps
/// `(lon, dlon)` through the ayanamsha subtraction. Both run even when
/// their flag is absent; the order is the only
/// thing this function owns. The always-call shape is a Rust-side
/// normalization: the model always runs the equatorial conversion
/// but gates the sidereal correction on its flag, so a sidereal
/// provider must be the identity off its flag.
pub fn apply_lunar_frame_order(
    state: LunarState6,
    jd_tt: f64,
    flags: i32,
    equatorial_convert: &dyn Fn(LunarState6, f64, i32) -> LunarState6,
    sidereal_correct: &dyn Fn(f64, f64) -> (f64, f64),
) -> LunarState6 {
    match lunar_frame_order(flags) {
        LunarFrameOrder::EquatorialThenSidereal => {
            let converted = equatorial_convert(state, jd_tt, flags);
            let (lon, dlon) = sidereal_correct(converted.0, converted.3);
            (
                lon,
                converted.1,
                converted.2,
                dlon,
                converted.4,
                converted.5,
            )
        }
        LunarFrameOrder::SiderealThenEquatorial => {
            let (lon, dlon) = sidereal_correct(state.0, state.3);
            let corrected = (lon, state.1, state.2, dlon, state.4, state.5);
            equatorial_convert(corrected, jd_tt, flags)
        }
        LunarFrameOrder::EquatorialOnly => equatorial_convert(state, jd_tt, flags),
    }
}

/// Shape a mean-point state onto the requested ecliptic.
///
/// `state` is the analytic mean-ecliptic-of-date state from the mean-point
/// kernel. Unless
/// [`no_nutation_frame`] holds, the nutation in longitude is added to
/// `lon` (`% 360` floored, hence `rem_euclid`) and — under `FLG_SPEED`
/// only — the nutation rate over the `0.001`-day stencil half
/// step is added to `dlon`. Without `FLG_SPEED` all three speeds are
/// zeroed. A nutation-provider failure propagates.
pub fn shape_mean_point_state(
    state: LunarState6,
    jd_tt: f64,
    flags: i32,
    dpsi_rad_at: &dyn Fn(f64) -> Result<f64, Error>,
    nut_rate_at: &dyn Fn(f64, f64) -> Result<f64, Error>,
) -> Result<LunarState6, Error> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let no_nut = no_nutation_frame(flags);
    let (mut lon, lat, dist, mut dlon, mut dlat, mut ddist) = state;
    if !no_nut {
        lon = (lon + dpsi_rad_at(jd_tt)?.to_degrees()).rem_euclid(360.0);
    }
    if (flags & FLG_SPEED) != 0 {
        if !no_nut {
            dlon += nut_rate_at(jd_tt, 0.001)?;
        }
    } else {
        dlon = 0.0;
        dlat = 0.0;
        ddist = 0.0;
    }
    Ok((lon, lat, dist, dlon, dlat, ddist))
}

/// Full mean-point dispatch: analytic state plus flag shaping
///
/// `speed3` is accepted for the shared backend call signature and
/// ignored by the mean-point dispatch layer; the model mean dispatch never sets it.
/// Anything but `MEAN_NODE`/`MEAN_APOG` is [`Error::UnknownBody`];
/// range enforcement is the caller's explicit [`LunarCoverage`].
pub fn mean_lunar_point_state(
    body_id: i32,
    jd_tt: f64,
    flags: i32,
    coverage: LunarCoverage,
    dpsi_rad_at: &dyn Fn(f64) -> Result<f64, Error>,
    nut_rate_at: &dyn Fn(f64, f64) -> Result<f64, Error>,
) -> Result<LunarState6, Error> {
    let state = if body_id == MEAN_NODE {
        mean_lunar_node_state(jd_tt, coverage, false)?
    } else if body_id == MEAN_APOG {
        mean_lunar_apogee_state(jd_tt, coverage, false)?
    } else {
        return Err(Error::UnknownBody { body_id });
    };
    shape_mean_point_state(state, jd_tt, flags, dpsi_rad_at, nut_rate_at)
}

/// Shape a true-point central position plus stencil speeds onto the
/// requested ecliptic.
///
/// `central` is the raw true-ecliptic-of-date position at `jd_tt`
/// (nutation included); `position_at` serves raw side-sample positions
/// for the speed stencil. Under [`no_nutation_frame`] the nutation in
/// longitude is subtracted from the central longitude. Under `FLG_SPEED` the
/// speeds come from the [`stencil_true_point_speed`] boundary-fallback
/// stencil at `half_step_days`, minus the nutation rate at that same
/// half step when the output is mean-ecliptic; otherwise speeds are
/// zeroed.
pub fn shape_true_point_state(
    central: LunarPosition3,
    position_at: &dyn Fn(f64) -> Result<LunarPosition3, Error>,
    jd_tt: f64,
    flags: i32,
    half_step_days: f64,
    dpsi_rad_at: &dyn Fn(f64) -> Result<f64, Error>,
    nut_rate_at: &dyn Fn(f64, f64) -> Result<f64, Error>,
) -> Result<LunarState6, Error> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let no_nut = no_nutation_frame(flags);
    let (mut lon, lat, dist) = central;
    if no_nut {
        lon = (lon - dpsi_rad_at(jd_tt)?.to_degrees()).rem_euclid(360.0);
    }
    let (mut dlon, mut dlat, mut ddist) = (0.0, 0.0, 0.0);
    if (flags & FLG_SPEED) != 0 {
        (dlon, dlat, ddist) = stencil_true_point_speed(position_at, jd_tt, half_step_days)?;
        if no_nut {
            dlon -= nut_rate_at(jd_tt, half_step_days)?;
        }
    }
    Ok((lon, lat, dist, dlon, dlat, ddist))
}

/// Full true-point dispatch over a true-ecliptic-of-date Moon-state
/// provider.
///
/// The ordinary speed half step is per-point (`0.05` day for the
/// smooth node curve, `0.002` day for the fast-structured apogee);
/// `FLG_SPEED3`/`FLG_TOPOCTR` requests use the finer `0.0001`-day
/// public sampling convention ([`select_true_point_half_step`]).
/// Side samples re-derive raw positions through the same true-point
/// geometry, so a provider miss on either side yields zero speeds
/// while a central miss propagates fail-closed. Anything but
/// `TRUE_NODE`/`OSCU_APOG` is [`Error::UnknownBody`].
pub fn true_lunar_point_state(
    body_id: i32,
    jd_tt: f64,
    flags: i32,
    true_state_at: &dyn Fn(f64) -> Result<EclipticMoonState, Error>,
    dpsi_rad_at: &dyn Fn(f64) -> Result<f64, Error>,
    nut_rate_at: &dyn Fn(f64, f64) -> Result<f64, Error>,
) -> Result<LunarState6, Error> {
    let ordinary_half_step_days = if body_id == TRUE_NODE {
        TRUE_NODE_HALF_STEP_DAYS
    } else if body_id == OSCU_APOG {
        OSCULATING_APSE_HALF_STEP_DAYS
    } else {
        return Err(Error::UnknownBody { body_id });
    };
    let fine = (flags & (FLG_SPEED3 | FLG_TOPOCTR)) != 0;
    let half_step_days = select_true_point_half_step(fine, ordinary_half_step_days)?;
    let raw_at = |jd: f64| -> Result<LunarPosition3, Error> {
        let (r, v) = true_state_at(jd)?;
        if body_id == TRUE_NODE {
            true_lunar_node_from_state(r, v)
        } else {
            true_lilith_from_state(r, v)
        }
    };
    let central = raw_at(jd_tt)?;
    shape_true_point_state(
        central,
        &raw_at,
        jd_tt,
        flags,
        half_step_days,
        dpsi_rad_at,
        nut_rate_at,
    )
}

/// Same interpolated-apsis speed-mode selection as the model
/// (`_calc_body` `INTP_APOG`/`INTP_PERG` branch):
/// `(use_speed3_state, use_speed3_velocity)`. `use_speed3_state`
/// mirrors the `speed3=` argument; `use_speed3_velocity` requests the
/// finer public sampling for the speeds even when the central state
/// spellings stays ordinary (SPEED + TOPOCTR without SPEED 3).
pub fn select_interp_speed_mode(flags: i32) -> (bool, bool) {
    let use_speed3_state = (flags & FLG_SPEED3) != 0;
    let use_speed3_velocity =
        (flags & FLG_SPEED) != 0 && (use_speed3_state || (flags & FLG_TOPOCTR) != 0);
    (use_speed3_state, use_speed3_velocity)
}

/// Shape an interpolated-apsis state onto the requested ecliptic
///
/// `state` is the stencil state for the requested speed mode and
/// `speed3_speeds` carries the finer-sampling speeds when
/// [`select_interp_speed_mode`] requests them without the SPEED 3
/// state spelling. The model is defined on the true ecliptic of date,
/// so under [`no_nutation_frame`] the nutation in longitude is
/// subtracted (not added, unlike the mean points). Under `FLG_SPEED`
/// the mean-ecliptic output additionally loses the nutation rate over
/// the `0.001`-day stencil half step; otherwise speeds are
/// zeroed.
pub fn shape_interp_point_state(
    state: LunarState6,
    speed3_speeds: Option<LunarPosition3>,
    jd_tt: f64,
    flags: i32,
    dpsi_rad_at: &dyn Fn(f64) -> Result<f64, Error>,
    nut_rate_at: &dyn Fn(f64, f64) -> Result<f64, Error>,
) -> Result<LunarState6, Error> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let no_nut = no_nutation_frame(flags);
    let (mut lon, lat, dist, mut dlon, mut dlat, mut ddist) = state;
    if let Some((s3_dlon, s3_dlat, s3_ddist)) = speed3_speeds {
        dlon = s3_dlon;
        dlat = s3_dlat;
        ddist = s3_ddist;
    }
    if no_nut {
        lon = (lon - dpsi_rad_at(jd_tt)?.to_degrees()).rem_euclid(360.0);
    }
    if (flags & FLG_SPEED) != 0 {
        if no_nut {
            dlon -= nut_rate_at(jd_tt, 0.001)?;
        }
    } else {
        dlon = 0.0;
        dlat = 0.0;
        ddist = 0.0;
    }
    Ok((lon, lat, dist, dlon, dlat, ddist))
}

/// Full interpolated-apsis dispatch: generated-model state plus flag
/// shaping.
///
/// The `speed3` stencil spelling is the smooth instantaneous rate the
/// the model serves for both SPEED and SPEED 3 on the interpolated-apsis path; the
/// state/velocity split of [`select_interp_speed_mode`] is preserved
/// as written even though both spellings evaluate the same stencil.
/// Anything but `INTP_APOG`/`INTP_PERG` is [`Error::UnknownBody`].
/// There is no range gate on this path; the residual-grid edge
/// taper owns the out-of-model behavior.
pub fn interp_lunar_point_state(
    body_id: i32,
    jd_tt: f64,
    flags: i32,
    dpsi_rad_at: &dyn Fn(f64) -> Result<f64, Error>,
    nut_rate_at: &dyn Fn(f64, f64) -> Result<f64, Error>,
) -> Result<LunarState6, Error> {
    if !is_interp_lunar_point(body_id) {
        return Err(Error::UnknownBody { body_id });
    }
    let (use_speed3_state, use_speed3_velocity) = select_interp_speed_mode(flags);
    let state = calc_interpolated_apse_state(jd_tt, body_id, use_speed3_state)?;
    let speed3_speeds = if use_speed3_velocity && !use_speed3_state {
        let speed3_state = calc_interpolated_apse_state(jd_tt, body_id, true)?;
        Some((speed3_state.3, speed3_state.4, speed3_state.5))
    } else {
        None
    };
    shape_interp_point_state(state, speed3_speeds, jd_tt, flags, dpsi_rad_at, nut_rate_at)
}

/// Top-level lunar-point dispatch: the HELCTR/BARYCTR zero sentinel,
/// then the per-family state plus flag shaping.
///
/// `coverage` is the caller's resolved mean-point window (mean-kernel
/// range-check layer);
/// `true_state_at` serves true-ecliptic-of-date Moon states for the
/// osculating points (JPL-kernel or sealed-LEB acquisition stays
/// provider-side). Frame application (sidereal correction,
/// equatorial/J2000 conversion) stays with
/// [`apply_lunar_frame_order`]; south-node aliases never reach here
/// (see [`is_south_lunar_node`]). A non-lunar `body_id` is
/// [`Error::UnknownBody`].
pub fn lunar_point_state(
    body_id: i32,
    jd_tt: f64,
    flags: i32,
    coverage: LunarCoverage,
    true_state_at: &dyn Fn(f64) -> Result<EclipticMoonState, Error>,
    dpsi_rad_at: &dyn Fn(f64) -> Result<f64, Error>,
    nut_rate_at: &dyn Fn(f64, f64) -> Result<f64, Error>,
) -> Result<LunarState6, Error> {
    if wants_lunar_helio_bary_zero(body_id, flags) {
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        return Ok((0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
    }
    if is_mean_lunar_point(body_id) {
        mean_lunar_point_state(body_id, jd_tt, flags, coverage, dpsi_rad_at, nut_rate_at)
    } else if is_true_lunar_point(body_id) {
        true_lunar_point_state(
            body_id,
            jd_tt,
            flags,
            true_state_at,
            dpsi_rad_at,
            nut_rate_at,
        )
    } else if is_interp_lunar_point(body_id) {
        interp_lunar_point_state(body_id, jd_tt, flags, dpsi_rad_at, nut_rate_at)
    } else {
        Err(Error::UnknownBody { body_id })
    }
}

/// Advisory fitted-model window for an interpolated apsidal body
/// : the exact
/// inclusive TT window from the generated residual-grid metadata
/// (`start + step * (count - 1)`), or `None` for another body. This
/// describes the model, not a binary file, and never refuses: past
/// the last sample the correction tapers rather than steps. The
/// grids are linked unconditionally in Rust, so the model
/// None-when-unavailable arm cannot trigger here.
pub fn interpolated_apse_model_window(body_id: i32) -> Option<(f64, f64)> {
    if body_id == INTP_APOG {
        let end = APOGEE_CT_JD_START + APOGEE_CT_STEP_DAYS * (APOGEE_CT_COUNT as f64 - 1.0);
        Some((APOGEE_CT_JD_START, end))
    } else if body_id == INTP_PERG {
        let end = PERIGEE_CT_JD_START + PERIGEE_CT_STEP_DAYS * (PERIGEE_CT_COUNT as f64 - 1.0);
        Some((PERIGEE_CT_JD_START, end))
    } else {
        None
    }
}

/// Residual-table grid and taper used by an interpolated apsidal source.
///
/// Returns `(first_tt_jd, last_tt_jd, step_days, taper_days)` for interpolated
/// apogee/perigee, or `None` for another body. The residual is piecewise linear
/// between these knots, and tapers linearly to zero outside the table. Offline
/// polynomial fitting must respect those derivative breaks, including the two
/// outer taper endpoints; this metadata is not a physical accuracy guarantee.
/// No file, provider, epoch query or error-producing operation is consumed.
pub fn interpolated_apse_source_grid(body_id: i32) -> Option<(f64, f64, f64, f64)> {
    let (start, end) = interpolated_apse_model_window(body_id)?;
    let (_, step) = interpolated_apse_table_grid(body_id == INTP_APOG);
    Some((start, end, step, APSE_EDGE_TAPER_DAYS))
}

/// Diagnostic comparison of the two True Lilith spellings over one
/// provider state.
///
/// Both spellings are the same maintained alias (the true-point layer's
/// `true_lilith_orbital_elements_from_state`), so the differences are
/// exactly zero by construction; the struct pins the model dict
/// shape (`eccentricity_vector`, `orbital_elements`, `lon_diff`,
/// `lat_diff`, `dist_diff`) for the serving binding — the third slot
/// carries apogee distance in AU, not eccentricity. The longitude fold
/// keeps the two-`if` (non-`elif`) the model branches.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TrueLilithComparison {
    /// `(lon, lat, dist)` from the eccentricity-vector spelling.
    pub eccentricity_vector: LunarPosition3,
    /// `(lon, lat, dist)` from the orbital-elements alias spelling.
    pub orbital_elements: LunarPosition3,
    /// Wrapped longitude difference in degrees.
    pub lon_diff: f64,
    /// Latitude difference in degrees.
    pub lat_diff: f64,
    /// Third-slot (distance) difference in AU.
    pub dist_diff: f64,
}

/// Run the [`TrueLilithComparison`] diagnostic over a provider state.
pub fn compare_true_lilith_methods_via(
    true_state_at: &dyn Fn(f64) -> Result<EclipticMoonState, Error>,
    jd_tt: f64,
) -> Result<TrueLilithComparison, Error> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let (r, v) = true_state_at(jd_tt)?;
    let eccentricity_vector = true_lilith_from_state(r, v)?;
    let orbital_elements = true_lilith_orbital_elements_from_state(r, v)?;
    let mut lon_diff = eccentricity_vector.0 - orbital_elements.0;
    if lon_diff > 180.0 {
        lon_diff -= 360.0;
    }
    if lon_diff < -180.0 {
        lon_diff += 360.0;
    }
    Ok(TrueLilithComparison {
        eccentricity_vector,
        orbital_elements,
        lon_diff,
        lat_diff: eccentricity_vector.1 - orbital_elements.1,
        dist_diff: eccentricity_vector.2 - orbital_elements.2,
    })
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::erfa::{DAS2R, TURN_ARCSEC};

    const J2000: f64 = 2451545.0;

    #[test]
    fn fundamental_polynomial_pins_zero_epoch_words() {
        // At t = 0 every Horner level collapses to c 0, so each argument
        // is exactly fmod(c 0, TURN) * DAS 2 R with no rounding beyond the
        // two final operations.
        let (l, f, node) = fundamental_l_f_node(J2000);
        assert_eq!(l, (485_868.249036 % TURN_ARCSEC) * DAS2R);
        assert_eq!(f, (335_779.526232 % TURN_ARCSEC) * DAS2R);
        assert_eq!(node, (450_160.398036 % TURN_ARCSEC) * DAS2R);
        let (d, m, mp, ff) = lunar_delaunay_arguments(J2000);
        assert_eq!(d, (1_072_260.703692 % TURN_ARCSEC) * DAS2R);
        assert_eq!(m, (1_287_104.793048 % TURN_ARCSEC) * DAS2R);
        assert_eq!(mp, l);
        assert_eq!(ff, f);
    }

    #[test]
    fn apse_plane_identity_holds_over_a_nodal_cycle() {
        // varpi_apogee = Omega + F - l + pi (mod 360), structurally.
        let mut jd = J2000 - 4000.0;
        while jd <= J2000 + 4000.0 {
            let (l, f, node) = fundamental_l_f_node(jd);
            let (plane, node_deg) = mean_apse_plane_node_unchecked(jd);
            let expect = (node + f - l + std::f64::consts::PI)
                .to_degrees()
                .rem_euclid(360.0);
            assert_eq!(plane, expect, "plane identity at {jd}");
            assert_eq!(node_deg, node.to_degrees().rem_euclid(360.0));
            assert!((0.0..360.0).contains(&plane));
            jd += 137.0;
        }
    }

    #[test]
    fn node_position_contract_at_j2000() {
        let (lon, lat, dist) =
            mean_lunar_node_position(J2000, LunarCoverage::Unenforced).expect("finite");
        assert!((0.0..360.0).contains(&lon));
        assert_eq!(lat, 0.0);
        assert_eq!(dist, MEAN_NODE_DISTANCE_AU);
        assert_eq!(dist, 384_400.0 / 149_597_870.7);
    }

    #[test]
    fn apogee_position_contract_at_j2000() {
        let (lon, lat, dist) =
            mean_lunar_apogee_position(J2000, LunarCoverage::Unenforced).expect("finite");
        assert!((0.0..360.0).contains(&lon));
        assert!(lat.abs() <= INCLINATION_DEG + 1e-9, "tilt bound: {lat}");
        assert_eq!(dist, MEAN_APOGEE_DISTANCE_AU);
        // States agree with positions on the shared slots.
        let st = mean_lunar_apogee_state(J2000, LunarCoverage::Unenforced, false).expect("state");
        assert_eq!((st.0, st.1, st.2), (lon, lat, dist));
        let nst = mean_lunar_node_state(J2000, LunarCoverage::Unenforced, false).expect("state");
        let npos = mean_lunar_node_position(J2000, LunarCoverage::Unenforced).expect("pos");
        assert_eq!((nst.0, nst.1, nst.2), npos);
    }

    #[test]
    fn nodal_and_apsidal_rates_have_release_sign_and_scale() {
        // The node regresses (~18.6-yr period) while the apogee plane
        // advances (~8.85-yr perigee cycle).
        let node = mean_lunar_node_state(J2000, LunarCoverage::Unenforced, false).expect("state");
        assert!(
            (-0.06..-0.045).contains(&node.3),
            "node rate deg/day: {}",
            node.3
        );
        assert_eq!(node.4, 0.0);
        assert_eq!(node.5, 0.0);
        let apogee =
            mean_lunar_apogee_state(J2000, LunarCoverage::Unenforced, false).expect("state");
        assert!(
            (0.08..0.14).contains(&apogee.3),
            "apogee rate deg/day: {}",
            apogee.3
        );
        assert_eq!(apogee.5, 0.0);
    }

    #[test]
    fn full_cycle_returns_mean_longitudes() {
        // independent theory-consistency check of the swept entries. The periods below are derived algebraically from the
        // published Simon et al. (1994) linear coefficients (arcsec/cy),
        // not fitted: P = 360*3600*36525/|rate|. After one full nodal
        // (~18.61 yr) or apsidal (~8.85 yr) cycle the mean node longitude
        // (`calc_mean_lunar_node`) and the mean apogee *plane* longitude
        // (`mean_lunar_apse_plane_node`) must return to their starting
        // values up to the higher-order polynomial content over that span
        // (node 0.26 arcsec ~ 7.2 e-5 deg, plane 0.29 arcsec ~ 8.1 e-5 deg;
        // binary64 rounding over ~1 e 4-day spans is ~1 e-12 deg). The 5 e-4
        // deg bound therefore keeps >2 x analytic margin even off J2000.
        //
        // The projected apogee ecliptic longitude/latitude are
        // deliberately not asserted here: they also depend on the node,
        // which advances ~half a cycle over one apsidal period, so they
        // do not return after a single apsidal cycle by construction.
        const NODE_PERIOD_DAYS: f64 = 360.0 * 3600.0 * 36525.0 / 6_962_890.543_1;
        const APOGEE_PERIOD_DAYS: f64 = 360.0 * 3600.0 * 36525.0
            / (1_739_527_262.847_8 - 1_717_915_923.217_8 - 6_962_890.543_1);
        const RETURN_BOUND_DEG: f64 = 5e-4;
        fn wrapped_delta(after: f64, before: f64) -> f64 {
            let d = (after - before).rem_euclid(360.0);
            if d > 180.0 { d - 360.0 } else { d }
        }
        let cov = LunarCoverage::Unenforced;
        for jd in [J2000, J2000 + 1000.0, J2000 - 1000.0, J2000 + 2000.0] {
            let node0 = calc_mean_lunar_node(jd, cov).expect("node");
            let node1 = calc_mean_lunar_node(jd + NODE_PERIOD_DAYS, cov).expect("node");
            let dnode = wrapped_delta(node1, node0);
            assert!(
                dnode.abs() <= RETURN_BOUND_DEG,
                "node return at {jd}: {dnode}"
            );
            let (plane0, _) = mean_lunar_apse_plane_node(jd, cov).expect("plane");
            let (plane1, _) =
                mean_lunar_apse_plane_node(jd + APOGEE_PERIOD_DAYS, cov).expect("plane");
            let dplane = wrapped_delta(plane1, plane0);
            assert!(
                dplane.abs() <= RETURN_BOUND_DEG,
                "apogee-plane return at {jd}: {dplane}"
            );
        }
    }

    #[test]
    fn speed3_flag_is_accepted_and_ignored() {
        for jd in [J2000, J2000 + 100.0, J2000 - 1000.0] {
            assert_eq!(
                mean_lunar_node_state(jd, LunarCoverage::Unenforced, false),
                mean_lunar_node_state(jd, LunarCoverage::Unenforced, true),
            );
            assert_eq!(
                mean_lunar_apogee_state(jd, LunarCoverage::Unenforced, false),
                mean_lunar_apogee_state(jd, LunarCoverage::Unenforced, true),
            );
        }
    }

    #[test]
    fn thin_wrappers_agree_with_mean_functions() {
        for jd in [J2000, J2000 + 45.25, J2000 - 800.5] {
            let cov = LunarCoverage::Unenforced;
            assert_eq!(
                calc_mean_lunar_node(jd, cov).expect("node"),
                mean_lunar_node_position(jd, cov).expect("pos").0
            );
            assert_eq!(
                calc_mean_lilith(jd, cov).expect("lilith"),
                mean_lunar_apogee_position(jd, cov).expect("pos").0
            );
            let (lon, lat) = calc_mean_lilith_with_latitude(jd, cov).expect("pair");
            let pos = mean_lunar_apogee_position(jd, cov).expect("pos");
            assert_eq!((lon, lat), (pos.0, pos.1));
            assert_eq!(
                calc_mean_lunar_node_state(jd, cov, false).expect("st"),
                mean_lunar_node_state(jd, cov, false).expect("st"),
            );
            assert_eq!(
                calc_mean_lilith_state(jd, cov, true).expect("st"),
                mean_lunar_apogee_state(jd, cov, true).expect("st"),
            );
        }
    }

    #[test]
    fn window_miss_names_the_mean_body() {
        let window = LunarCoverage::Window {
            start_jd: 2_451_545.0,
            end_jd: 2_451_910.0,
        };
        assert_eq!(
            mean_lunar_node_position(2_451_911.0, window),
            Err(Error::EphemerisRange {
                requested_jd: 2_451_911.0,
                start_jd: 2_451_545.0,
                end_jd: 2_451_910.0,
                body: Some(MEAN_NODE),
            })
        );
        assert_eq!(
            mean_lunar_apogee_position(2_451_544.0, window),
            Err(Error::EphemerisRange {
                requested_jd: 2_451_544.0,
                start_jd: 2_451_545.0,
                end_jd: 2_451_910.0,
                body: Some(MEAN_APOG),
            })
        );
        // Closed endpoints serve; the (0, 0) window is unenforced like
        // the model falsy branch.
        assert!(mean_lunar_node_position(2_451_545.0, window).is_ok());
        assert!(mean_lunar_node_position(2_451_910.0, window).is_ok());
        let zero = LunarCoverage::Window {
            start_jd: 0.0,
            end_jd: 0.0,
        };
        assert!(mean_lunar_node_position(2_000_000.0, zero).is_ok());
    }

    #[test]
    fn non_finite_epochs_are_typed_refusals() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            // NaN never equals itself, so match the variant structurally.
            let expect_nan = bad.is_nan();
            for err in [
                mean_lunar_node_position(bad, LunarCoverage::Unenforced).map(|_| ()),
                mean_lunar_apogee_position(bad, LunarCoverage::Unenforced).map(|_| ()),
                mean_lunar_apse_plane_node(bad, LunarCoverage::Unenforced).map(|_| ()),
            ] {
                match err {
                    Err(Error::NonFiniteArgument { value }) => {
                        assert_eq!(value.is_nan(), expect_nan);
                        if !expect_nan {
                            assert_eq!(value, bad);
                        }
                    }
                    _ => panic!("expected NonFiniteArgument for {bad:?}, got {err:?}"),
                }
            }
            match mean_lunar_apogee_state(bad, LunarCoverage::Unenforced, false) {
                Err(Error::NonFiniteArgument { value }) => {
                    assert_eq!(value.is_nan(), expect_nan);
                }
                other => panic!("expected NonFiniteArgument, got {other:?}"),
            }
        }
        // The unchecked kernels propagate non-finite inputs without
        // panicking, using the same physical model float flow.
        let (lon, _, _) = mean_lunar_node_position_unchecked(f64::NAN);
        assert!(lon.is_nan());
    }

    #[test]
    fn meeus_bands_follow_release_thresholds() {
        assert_eq!(meeus_band(0.0), MeeusBand::Optimal);
        assert_eq!(meeus_band(2.0), MeeusBand::Optimal);
        assert_eq!(meeus_band(-10.0), MeeusBand::Optimal);
        assert_eq!(meeus_band(10.0001), MeeusBand::Degraded);
        assert_eq!(meeus_band(-15.0), MeeusBand::Degraded);
        assert_eq!(meeus_band(20.0), MeeusBand::Degraded);
        assert_eq!(meeus_band(20.0001), MeeusBand::OutsideValid);
        assert_eq!(meeus_band(-25.0), MeeusBand::OutsideValid);
        // The sweep domain (1990-2035) sits deep inside Optimal.
        assert_eq!(
            meeus_band(centuries_since_j2000(2_451_545.0)),
            MeeusBand::Optimal
        );
    }

    #[test]
    fn apogee_latitude_respects_the_tilt_over_time() {
        // The tilted-circle latitude can never exceed the inclination.
        let mut jd = J2000 - 5000.0;
        while jd <= J2000 + 5000.0 {
            let (_, lat, _) = mean_lunar_apogee_position_unchecked(jd);
            assert!(
                lat.abs() <= INCLINATION_DEG + 1e-9,
                "tilt bound at {jd}: {lat}"
            );
            jd += 61.0;
        }
    }

    // -- true-point geometry ----------------------------------

    #[test]
    fn gm_constant_matches_the_release_word() {
        // The Rust spelling of the model expression must evaluate to
        // the model value (checked ephemerally against the pinned
        // the model: 0 x 3 e 0 ee 9 dadb 7211 c 1).
        assert_eq!(
            GM_EARTH_MOON_AU3_DAY2.to_bits(),
            0x3e0e_e9da_db72_11c1,
            "GM spelling diverged from the release evaluation"
        );
        // The node fallback distance is the single-division mean
        // distance, identical to the constant.
        assert_eq!(MEAN_NODE_DISTANCE_AU.to_bits(), 0x3f65_0cbf_7deb_2640);
    }

    #[test]
    fn equatorial_circular_orbit_hits_the_node_fallback() {
        // Equatorial prograde circle: no node direction exists, so the
        // the model branch falls back to the mean distance while the
        // longitude lands on atan 2(0, -0) = pi exactly.
        let d = MEAN_NODE_DISTANCE_AU;
        let n = 0.23; // any positive mean motion in rad/day
        let (lon, lat, dist) =
            true_lunar_node_from_state((d, 0.0, 0.0), (0.0, n * d, 0.0)).expect("finite");
        assert_eq!(lon, 180.0);
        assert_eq!(lat, 0.0);
        assert_eq!(dist, MEAN_NODE_DISTANCE_AU);
    }

    /// Rotate a perifocal vector to the ecliptic through
    /// `R3(node) R1(incl) R3(argp)` (independent construction path for
    /// the synthetic-element tests: plain rotation matrices, never the
    /// eccentricity-vector method under test).
    fn perifocal_to_ecliptic(
        v: (f64, f64, f64),
        node: f64,
        incl: f64,
        argp: f64,
    ) -> (f64, f64, f64) {
        let apply_r3 = |v: (f64, f64, f64), a: f64| {
            (
                v.0 * a.cos() - v.1 * a.sin(),
                v.0 * a.sin() + v.1 * a.cos(),
                v.2,
            )
        };
        let apply_r1 = |v: (f64, f64, f64), a: f64| {
            (
                v.0,
                v.1 * a.cos() - v.2 * a.sin(),
                v.1 * a.sin() + v.2 * a.cos(),
            )
        };
        apply_r3(apply_r1(apply_r3(v, argp), incl), node)
    }

    /// Geocentric Moon state on a Kepler ellipse with the given
    /// elements (angles in radians, `a` in AU), via the matrix path.
    fn kepler_state(
        a: f64,
        e: f64,
        incl: f64,
        node: f64,
        argp: f64,
        anomaly: f64,
    ) -> ((f64, f64, f64), (f64, f64, f64)) {
        let p = a * (1.0 - e * e);
        let radius = p / (1.0 + e * anomaly.cos());
        let r_pf = (radius * anomaly.cos(), radius * anomaly.sin(), 0.0);
        let rate = (GM_EARTH_MOON_AU3_DAY2 / p).sqrt();
        let v_pf = (rate * -anomaly.sin(), rate * (e + anomaly.cos()), 0.0);
        (
            perifocal_to_ecliptic(r_pf, node, incl, argp),
            perifocal_to_ecliptic(v_pf, node, incl, argp),
        )
    }

    #[test]
    fn synthetic_elements_recover_node_and_apses() {
        let deg = std::f64::consts::PI / 180.0;
        let (a, e) = (0.00257, 0.0549);
        let (incl, node, argp) = (5.145396 * deg, 125.08 * deg, 318.06 * deg);
        let (r, v) = kepler_state(a, e, incl, node, argp, 1.0);

        // Ascending node: the model longitude must recover Omega.
        let (lon, lat, dist) = true_lunar_node_from_state(r, v).expect("finite");
        assert!((lon - 125.08).abs() < 1e-9, "node lon: {lon}");
        assert_eq!(lat, 0.0);
        let p = a * (1.0 - e * e);
        let expect_dist = p / (1.0 + e * (2.0 * std::f64::consts::PI - argp).cos());
        assert!(
            (dist - expect_dist).abs() / expect_dist < 1e-12,
            "node dist: {dist} vs {expect_dist}"
        );

        // Apogee: perifocal (-1, 0, 0) rotated to the ecliptic, at a(1+e).
        let apogee_dir = perifocal_to_ecliptic((-1.0, 0.0, 0.0), node, incl, argp);
        let expect_lon = apogee_dir
            .1
            .atan2(apogee_dir.0)
            .to_degrees()
            .rem_euclid(360.0);
        let expect_lat = apogee_dir.2.asin().to_degrees();
        let (lon, lat, dist) = true_lilith_from_state(r, v).expect("finite");
        assert!(
            (lon - expect_lon).abs() < 1e-9,
            "apogee lon: {lon} vs {expect_lon}"
        );
        assert!(
            (lat - expect_lat).abs() < 1e-9,
            "apogee lat: {lat} vs {expect_lat}"
        );
        assert!(
            (dist - a * (1.0 + e)).abs() / (a * (1.0 + e)) < 1e-12,
            "apogee dist: {dist}"
        );

        // The orbital-elements alias is exactly the primary path.
        assert_eq!(
            true_lilith_orbital_elements_from_state(r, v).expect("finite"),
            (lon, lat, dist)
        );

        // Perigee: opposite direction at a(1-e).
        let perigee_dir = perifocal_to_ecliptic((1.0, 0.0, 0.0), node, incl, argp);
        let expect_lon = perigee_dir
            .1
            .atan2(perigee_dir.0)
            .to_degrees()
            .rem_euclid(360.0);
        let expect_lat = perigee_dir.2.asin().to_degrees();
        let (lon, lat, dist) = osculating_perigee_from_state(r, v).expect("finite");
        assert!(
            (lon - expect_lon).abs() < 1e-9,
            "perigee lon: {lon} vs {expect_lon}"
        );
        assert!(
            (lat - expect_lat).abs() < 1e-9,
            "perigee lat: {lat} vs {expect_lat}"
        );
        assert!(
            (dist - a * (1.0 - e)).abs() / (a * (1.0 - e)) < 1e-12,
            "perigee dist: {dist}"
        );
    }

    #[test]
    fn non_finite_true_point_states_are_typed_refusals() {
        let good = (0.0025, 0.0001, 0.0002);
        let good_v = (-0.0001, 0.0006, 0.00001);
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            for r in [
                (bad, good.1, good.2),
                (good.0, bad, good.2),
                (good.0, good.1, bad),
            ] {
                for err in [
                    true_lunar_node_from_state(r, good_v).map(|_| ()),
                    true_lilith_from_state(r, good_v).map(|_| ()),
                    osculating_perigee_from_state(r, good_v).map(|_| ()),
                ] {
                    assert!(
                        matches!(err, Err(Error::NonFiniteArgument { value }) if value == bad || (value.is_nan() && bad.is_nan())),
                        "expected NonFiniteArgument, got {err:?}"
                    );
                }
            }
            for v in [
                (bad, good_v.1, good_v.2),
                (good_v.0, bad, good_v.2),
                (good_v.0, good_v.1, bad),
            ] {
                assert!(true_lunar_node_from_state(good, v).is_err());
                assert!(true_lilith_from_state(good, v).is_err());
                assert!(osculating_perigee_from_state(good, v).is_err());
            }
        }
    }

    #[test]
    fn zero_state_degeneracy_is_loud_never_a_panic() {
        // Exactly-zero state: the node path takes the model fallback
        // branch (`n_mag == 0`, the NaN eccentricity never selects).
        let zero = (0.0, 0.0, 0.0);
        assert_eq!(
            true_lunar_node_from_state(zero, zero).expect("node flows"),
            (180.0, 0.0, MEAN_NODE_DISTANCE_AU)
        );
        // The apse paths divide by the vanishing eccentricity where the
        // the model raises ZeroDivisionError; Rust IEEE division propagates
        // NaN instead (documented divergence, binding-side mapping).
        let (lon, lat, dist) = true_lilith_from_state(zero, zero).expect("apse flows");
        assert!(lon.is_nan() && lat.is_nan() && dist.is_nan());
        let (lon, lat, dist) = osculating_perigee_from_state(zero, zero).expect("apse flows");
        assert!(lon.is_nan() && lat.is_nan() && dist.is_nan());
    }

    #[test]
    fn stencil_folds_the_seam_and_scales_plain_slots() {
        // Forward seam crossing: +0.2 deg over a 0.1-day span.
        let (dlon, dlat, ddist) =
            central_difference_polar_rates((359.9, 1.0, 0.0025), (0.1, 1.2, 0.0027), 0.05);
        assert!((dlon - 2.0).abs() < 1e-12, "folded dlon: {dlon}");
        assert!((dlat - 2.0).abs() < 1e-12, "dlat: {dlat}");
        assert!((ddist - 0.002).abs() < 1e-15, "ddist: {ddist}");
        // Backward seam crossing folds the other way.
        let (dlon, _, _) = central_difference_polar_rates((0.1, 0.0, 0.0), (359.9, 0.0, 0.0), 0.05);
        assert!((dlon + 2.0).abs() < 1e-12, "backward dlon: {dlon}");
        // The strict branches leave exactly-180 differences unfolded,
        // like the model `if/elif`.
        let (dlon, _, _) = central_difference_polar_rates((0.0, 0.0, 0.0), (180.0, 0.0, 0.0), 0.05);
        assert!((dlon - 1800.0).abs() < 1e-9, "boundary dlon: {dlon}");
    }

    #[test]
    fn half_step_selector_mirrors_the_flag_branch() {
        assert_eq!(
            select_true_point_half_step(false, TRUE_NODE_HALF_STEP_DAYS).expect("finite"),
            0.05
        );
        assert_eq!(
            select_true_point_half_step(false, OSCULATING_APSE_HALF_STEP_DAYS).expect("finite"),
            0.002
        );
        assert_eq!(
            select_true_point_half_step(true, TRUE_NODE_HALF_STEP_DAYS).expect("fine"),
            FINE_POINT_HALF_STEP_DAYS
        );
        assert_eq!(
            select_true_point_half_step(true, OSCULATING_APSE_HALF_STEP_DAYS).expect("fine"),
            0.0001
        );
        assert!(select_true_point_half_step(false, f64::NAN).is_err());
    }

    #[test]
    fn stencil_provider_failure_falls_back_to_zero() {
        let provider = |jd: f64| -> Result<(f64, f64, f64), Error> {
            Ok((100.0 + 0.5 * (jd - J2000), 2.0, 0.0026))
        };
        // Healthy provider matches the pure stencil core.
        let rates = stencil_true_point_speed(&provider, J2000, 0.05).expect("finite");
        assert_eq!(
            rates,
            central_difference_polar_rates(
                provider(J2000 - 0.05).expect("sample"),
                provider(J2000 + 0.05).expect("sample"),
                0.05
            )
        );
        // Healthy provider matches the pure stencil core (tolerance
        // covers the (jd - J2000) large-epoch cancellation: ulp(2.4 e 6)
        // ~= 4.7 e-10 day per side, amplified by the 0.1-day span).
        assert!((rates.0 - 0.5).abs() < 1e-8);
        // A miss on either side yields
        // zeros rather than propagating.
        let failing = |jd: f64| -> Result<(f64, f64, f64), Error> {
            if jd < J2000 {
                return Err(Error::EphemerisRange {
                    requested_jd: jd,
                    start_jd: J2000,
                    end_jd: J2000 + 10.0,
                    body: None,
                });
            }
            provider(jd)
        };
        assert_eq!(
            stencil_true_point_speed(&failing, J2000, 0.05).expect("fallback"),
            (0.0, 0.0, 0.0)
        );
        // Non-finite epochs are refused, not sampled.
        assert!(stencil_true_point_speed(&provider, f64::NAN, 0.05).is_err());
        assert!(stencil_true_point_speed(&provider, J2000, f64::INFINITY).is_err());
    }

    // -- interpolated apsides -----------------------------------

    #[test]
    fn apse_table_metadata_matches_the_fitted_windows() {
        use crate::leb::interpolated_apse_model_window;
        use crate::lunar_apse_tables as tables;

        assert_eq!(tables::APOGEE_TRIG_TERMS.len(), 34);
        assert_eq!(tables::PERIGEE_TRIG_TERMS.len(), 66);
        assert_eq!(tables::APOGEE_CORRECTIONS.len(), 40250);
        assert_eq!(tables::PERIGEE_CORRECTIONS.len(), 201249);
        assert_eq!(tables::APOGEE_CT_COUNT, 40250);
        assert_eq!(tables::PERIGEE_CT_COUNT, 201249);
        assert_eq!(
            tables::APOGEE_LAT_COEFFS,
            (5.148886522214669, 0.0036009601970691863)
        );
        assert_eq!(
            tables::PERIGEE_LAT_COEFFS,
            (5.090870865743963, 0.0035299982729599375)
        );
        assert_eq!(tables::APOGEE_MEAN_DIST_AU, 0.0027099776307836668);
        assert_eq!(tables::PERIGEE_MEAN_DIST_AU, 0.0024236498902102843);
        // Grid ends agree with the LEB-09 fitted-window derivation:
        // apogee 2286820.5 + 10 * 40249 = 2689310.5,
        // perigee 2286820.5 + 2 * 201248 = 2689316.5.
        assert_eq!(
            tables::APOGEE_CT_JD_START + tables::APOGEE_CT_STEP_DAYS * 40249.0,
            2689310.5
        );
        assert_eq!(
            tables::PERIGEE_CT_JD_START + tables::PERIGEE_CT_STEP_DAYS * 201248.0,
            2689316.5
        );
        assert_eq!(
            interpolated_apse_model_window(INTP_APOG),
            Some((2286820.5, 2689310.5))
        );
        assert_eq!(
            interpolated_apse_model_window(INTP_PERG),
            Some((2286820.5, 2689316.5))
        );
    }

    #[test]
    fn table_correction_pins_first_and_last_words() {
        use crate::lunar_apse_tables as tables;

        // the model edge values (first/last table entries).
        assert_eq!(tables::APOGEE_CORRECTIONS[0], -0.0);
        assert_eq!(tables::APOGEE_CORRECTIONS[40249], 2.553);
        assert_eq!(tables::PERIGEE_CORRECTIONS[0], 0.0);
        assert_eq!(tables::PERIGEE_CORRECTIONS[201248], 0.0);
        // Exact grid nodes reproduce the table value exactly: the
        // fractional part is 0.0, so `c + 0.0 * _` is exact.
        assert_eq!(
            interpolate_apse_correction(
                2286820.5,
                2286820.5,
                10.0,
                tables::APOGEE_CORRECTIONS,
                40250
            )
            .expect("node"),
            tables::APOGEE_CORRECTIONS[0]
        );
        assert_eq!(
            interpolate_apse_correction(
                2286820.5 + 10.0,
                2286820.5,
                10.0,
                tables::APOGEE_CORRECTIONS,
                40250
            )
            .expect("node"),
            tables::APOGEE_CORRECTIONS[1]
        );
        // Interior midpoints are the exact average (binary midpoint of
        // two finite values halved once is exact up to the final add).
        let mid = interpolate_apse_correction(
            2286820.5 + 5.0,
            2286820.5,
            10.0,
            tables::APOGEE_CORRECTIONS,
            40250,
        )
        .expect("mid");
        let expect = tables::APOGEE_CORRECTIONS[0]
            + 0.5 * (tables::APOGEE_CORRECTIONS[1] - tables::APOGEE_CORRECTIONS[0]);
        assert_eq!(mid, expect);
    }

    #[test]
    fn table_edge_taper_follows_the_release_branches() {
        use crate::lunar_apse_tables as tables;

        // One day before the apogee grid: the edge value scaled by
        // 1 - 1/365.25 (association order).
        let before = interpolate_apse_correction(
            2286820.5 - 1.0,
            2286820.5,
            10.0,
            tables::APOGEE_CORRECTIONS,
            40250,
        )
        .expect("taper");
        assert_eq!(before, tables::APOGEE_CORRECTIONS[0] * (1.0 - 1.0 / 365.25));
        // A full taper year out the correction is exactly 0.0, and stays
        // 0.0 beyond (weight <= 0 branch).
        assert_eq!(
            interpolate_apse_correction(
                2286820.5 - 365.25,
                2286820.5,
                10.0,
                tables::APOGEE_CORRECTIONS,
                40250
            )
            .expect("edge"),
            0.0
        );
        assert_eq!(
            interpolate_apse_correction(
                2286820.5 - 1000.0,
                2286820.5,
                10.0,
                tables::APOGEE_CORRECTIONS,
                40250
            )
            .expect("far"),
            0.0
        );
        // Past the last apogee sample the live edge value (2.553") tapers.
        // The reference evaluates the recurrence in the calculation order
        // order: `(idx - (n-1)) * step` is not exactly 2.0 in binary64.
        let after = interpolate_apse_correction(
            2689310.5 + 2.0,
            2286820.5,
            10.0,
            tables::APOGEE_CORRECTIONS,
            40250,
        )
        .expect("taper");
        let idx = (2689310.5 + 2.0 - 2286820.5) / 10.0;
        let dist = (idx - (40250.0 - 1.0)) * 10.0;
        assert_eq!(after, 2.553 * (1.0 - dist / 365.25));
        assert_eq!(
            interpolate_apse_correction(
                2689310.5 + 400.0,
                2286820.5,
                10.0,
                tables::APOGEE_CORRECTIONS,
                40250
            )
            .expect("far"),
            0.0
        );
        // Empty tables yield 0.0.
        assert_eq!(
            interpolate_apse_correction(2451545.0, 2286820.5, 10.0, &[], 0).expect("empty"),
            0.0
        );
        // NaN epochs are refused.
        assert!(matches!(
            interpolate_apse_correction(
                f64::NAN,
                2286820.5,
                10.0,
                tables::APOGEE_CORRECTIONS,
                40250
            ),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[test]
    fn interpolated_positions_share_state_slots_and_mean_distances() {
        for jd in [2286820.5, 2451545.0, 2500000.0, 2689316.5] {
            let apogee = calc_interpolated_apogee(jd).expect("finite");
            let perigee = calc_interpolated_perigee(jd).expect("finite");
            assert!((0.0..360.0).contains(&apogee.0), "apogee lon at {jd}");
            assert!((0.0..360.0).contains(&perigee.0), "perigee lon at {jd}");
            assert!(apogee.1.abs() < 10.0 && perigee.1.abs() < 10.0);
            assert_eq!(apogee.2, APOGEE_MEAN_DIST_AU);
            assert_eq!(perigee.2, PERIGEE_MEAN_DIST_AU);
            assert!(apogee.2 > perigee.2);
            // States agree with positions on the shared slots, and
            // speed 3 is accepted and ignored like the model.
            let apogee_state = calc_interpolated_apse_state(jd, INTP_APOG, false).expect("state");
            assert_eq!((apogee_state.0, apogee_state.1, apogee_state.2), apogee);
            assert_eq!(
                calc_interpolated_apse_state(jd, INTP_APOG, true).expect("speed3"),
                apogee_state
            );
            let perigee_state = calc_interpolated_apse_state(jd, INTP_PERG, false).expect("state");
            assert_eq!((perigee_state.0, perigee_state.1, perigee_state.2), perigee);
            // The two curves are fitted independently: not antipodal
            let mut separation = (apogee.0 - perigee.0).rem_euclid(360.0);
            if separation > 180.0 {
                separation = 360.0 - separation;
            }
            assert!(
                (separation - 180.0).abs() > 0.5,
                "unexpected antipode at {jd}: {separation}"
            );
        }
    }

    #[test]
    fn interpolated_trig_terms_stay_in_release_plausibility() {
        // Trig-only RMS is 0.045 deg (apogee) / 0.112 deg (perigee)
        // against the passage track; the instantaneous perturbation
        // stays at the few-degree scale across the whole window.
        let mut jd = 2286820.5;
        while jd <= 2689316.5 {
            assert!(
                de440_apogee_passage_terms(jd).abs() < 10.0,
                "apogee at {jd}"
            );
            assert!(
                de440_perigee_passage_terms(jd).abs() < 30.0,
                "perigee at {jd}"
            );
            jd += 1000.0;
        }
    }

    #[test]
    fn interpolated_stencil_folds_the_seam_and_refuses_bad_bodies() {
        // A synthetic provider straddling 0 deg: the wrap fold keeps the
        // rate on the short way (10 deg over the 0.004-day span).
        let provider = |jd: f64| -> Result<(f64, f64, f64), Error> {
            Ok((if jd < J2000 { 355.0 } else { 5.0 }, 1.0, 0.0026))
        };
        let state = finite_difference_polar_state(&provider, J2000, 0.002).expect("fold");
        assert_eq!(state.0, 5.0);
        assert!((state.3 - 2500.0).abs() < 1e-9, "folded rate: {}", state.3);
        // Unsupported bodies name the caller's key; non-finite epochs
        // and steps are typed refusals.
        assert_eq!(
            calc_interpolated_apse_state(J2000, MEAN_NODE, false),
            Err(Error::UnknownBody { body_id: MEAN_NODE })
        );
        assert!(matches!(
            calc_interpolated_apse_state(f64::NAN, INTP_APOG, false),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            finite_difference_polar_state(&provider, J2000, f64::INFINITY),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(calc_interpolated_apogee(f64::NAN).is_err());
        assert!(calc_interpolated_perigee(f64::INFINITY).is_err());
    }

    // ------------------------------------------------------------------
    // dispatch tests (analytic stub providers, no kernels)
    // ------------------------------------------------------------------

    /// Fixed stub nutation: dpsi = +0.01 deg in radians, rate = +2 e-6
    /// deg/day regardless of epoch (scale-correct, keeps every shaping
    /// branch arithmetically exact).
    fn stub_dpsi_rad(_jd: f64) -> Result<f64, Error> {
        Ok(0.01_f64.to_radians())
    }

    fn stub_nut_rate(_jd: f64, _half_step: f64) -> Result<f64, Error> {
        Ok(2e-6)
    }

    /// Circular equatorial-orbit stub state: node at lon 0, apogee at
    /// lon 90 with e = 0.05. Radius 0.0025 AU, circular speed
    /// `sqrt(mu/r)`.
    fn stub_true_state(_jd: f64) -> Result<EclipticMoonState, Error> {
        let r = (0.0025, 0.0, 0.0);
        let speed = (GM_EARTH_MOON_AU3_DAY2 / 0.0025).sqrt();
        Ok((r, (0.0, speed, 0.0)))
    }

    #[test]
    fn dispatch_predicates_cover_six_points_and_south_aliases() {
        use crate::constants::{INTP_APOG, INTP_PERG, MEAN_APOG, MEAN_NODE, OSCU_APOG, TRUE_NODE};
        for body in [
            MEAN_NODE, TRUE_NODE, MEAN_APOG, OSCU_APOG, INTP_APOG, INTP_PERG,
        ] {
            assert!(is_lunar_point(body), "{body}");
        }
        for body in [0, 1, 10 + 100, MEAN_NODE + 100] {
            assert!(!is_lunar_point(body), "{body}");
        }
        assert!(is_mean_lunar_point(MEAN_NODE) && is_mean_lunar_point(MEAN_APOG));
        assert!(!is_mean_lunar_point(TRUE_NODE));
        assert!(is_true_lunar_point(TRUE_NODE) && is_true_lunar_point(OSCU_APOG));
        assert!(!is_true_lunar_point(INTP_APOG));
        assert!(is_interp_lunar_point(INTP_APOG) && is_interp_lunar_point(INTP_PERG));
        assert!(!is_interp_lunar_point(MEAN_APOG));
        assert!(is_lunar_node(MEAN_NODE) && is_lunar_node(TRUE_NODE));
        assert!(!is_lunar_node(MEAN_APOG));
        // South aliases route to their north source; nothing else does.
        assert!(is_south_lunar_node(-MEAN_NODE) && is_south_lunar_node(-TRUE_NODE));
        assert!(!is_south_lunar_node(MEAN_NODE));
        assert_eq!(south_node_north_source(-MEAN_NODE), Some(MEAN_NODE));
        assert_eq!(south_node_north_source(-TRUE_NODE), Some(TRUE_NODE));
        assert_eq!(south_node_north_source(MEAN_NODE), None);
        // Model window is advisory and exact from the grid metadata.
        let (start, end) = interpolated_apse_model_window(INTP_APOG).expect("window");
        assert_eq!(start, APOGEE_CT_JD_START);
        assert_eq!(
            end,
            APOGEE_CT_JD_START + APOGEE_CT_STEP_DAYS * (APOGEE_CT_COUNT as f64 - 1.0)
        );
        assert!(interpolated_apse_model_window(MEAN_NODE).is_none());
    }

    #[test]
    fn helio_bary_sentinel_returns_exact_zeros() {
        use crate::constants::{FLG_BARYCTR, FLG_HELCTR};
        for flags in [
            FLG_HELCTR,
            FLG_BARYCTR,
            FLG_HELCTR | FLG_BARYCTR | FLG_SPEED,
        ] {
            assert!(wants_lunar_helio_bary_zero(MEAN_NODE, flags));
            let state = lunar_point_state(
                MEAN_NODE,
                J2000,
                flags,
                LunarCoverage::Unenforced,
                &stub_true_state,
                &stub_dpsi_rad,
                &stub_nut_rate,
            )
            .expect("sentinel");
            assert_eq!(state, (0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
        }
        assert!(!wants_lunar_helio_bary_zero(MEAN_NODE, FLG_SPEED));
        assert!(!wants_lunar_helio_bary_zero(1, FLG_HELCTR));
        // Non-lunar bodies are not served here at all.
        assert_eq!(
            lunar_point_state(
                1,
                J2000,
                0,
                LunarCoverage::Unenforced,
                &stub_true_state,
                &stub_dpsi_rad,
                &stub_nut_rate
            ),
            Err(Error::UnknownBody { body_id: 1 })
        );
    }

    #[test]
    fn no_nutation_frame_rule_matches_release_branches() {
        use crate::constants::{FLG_EQUATORIAL, FLG_J2000, FLG_NONUT, FLG_SIDEREAL, FLG_SPEED};
        assert!(!no_nutation_frame(FLG_SPEED));
        assert!(no_nutation_frame(FLG_NONUT));
        assert!(no_nutation_frame(FLG_J2000));
        assert!(no_nutation_frame(FLG_SIDEREAL | FLG_EQUATORIAL));
        assert!(!no_nutation_frame(FLG_SIDEREAL));
        assert!(!no_nutation_frame(FLG_EQUATORIAL));
        // Frame order: SIDEREAL+J2000 ecliptic converts first.
        assert_eq!(
            lunar_frame_order(FLG_SIDEREAL | FLG_J2000),
            LunarFrameOrder::EquatorialThenSidereal
        );
        assert_eq!(
            lunar_frame_order(FLG_SIDEREAL),
            LunarFrameOrder::SiderealThenEquatorial
        );
        assert_eq!(lunar_frame_order(0), LunarFrameOrder::EquatorialOnly);
        assert_eq!(
            lunar_frame_order(FLG_SIDEREAL | FLG_EQUATORIAL),
            LunarFrameOrder::EquatorialOnly
        );
        // The order is the only thing owned: identity providers expose
        // the input in every branch, swapped providers show the order.
        let state = (10.0, 1.0, 0.0026, 0.5, 0.01, 0.0);
        let ident = |s: (f64, f64, f64, f64, f64, f64), _: f64, _: i32| s;
        let sid = |lon: f64, dlon: f64| (lon + 1.0, dlon + 0.1);
        for flags in [0, FLG_SIDEREAL, FLG_SIDEREAL | FLG_J2000] {
            let out = apply_lunar_frame_order(state, J2000, flags, &ident, &sid);
            // The tropical branch converts only (identity here); both
            // sidereal branches correct (order differs, values agree on
            // identity conversion).
            let expect = if flags == 0 { (10.0, 0.5) } else { (11.0, 0.6) };
            assert_eq!((out.0, out.3), expect, "flags {flags}");
            assert_eq!((out.1, out.2, out.4, out.5), (1.0, 0.0026, 0.01, 0.0));
        }
        // The order is structural: with non-commuting providers the two
        // sidereal branches disagree exactly as the model branch order
        // dictates (convert-then-correct vs correct-then-convert).
        let times_two = |s: (f64, f64, f64, f64, f64, f64), _: f64, _: i32| {
            (s.0 * 2.0, s.1, s.2, s.3, s.4, s.5)
        };
        let out = apply_lunar_frame_order(state, J2000, FLG_SIDEREAL | FLG_J2000, &times_two, &sid);
        assert_eq!((out.0, out.3), (21.0, 0.6));
        let out = apply_lunar_frame_order(state, J2000, FLG_SIDEREAL, &times_two, &sid);
        assert_eq!((out.0, out.3), (22.0, 0.6));
    }

    #[test]
    fn mean_dispatch_adds_nutation_and_gates_speed() {
        use crate::constants::FLG_SPEED;
        // SPEED: raw state + dpsi on lon, + rate on dlon.
        let out = mean_lunar_point_state(
            MEAN_NODE,
            J2000,
            FLG_SPEED,
            LunarCoverage::Unenforced,
            &stub_dpsi_rad,
            &stub_nut_rate,
        )
        .expect("mean speed");
        let raw = mean_lunar_node_state(J2000, LunarCoverage::Unenforced, false).expect("raw");
        assert_eq!(out.0, (raw.0 + 0.01).rem_euclid(360.0));
        assert_eq!((out.1, out.2), (raw.1, raw.2));
        assert_eq!(out.3, raw.3 + 2e-6);
        assert_eq!((out.4, out.5), (raw.4, raw.5));
        // No SPEED: speeds zeroed, position still nutated.
        let out = mean_lunar_point_state(
            MEAN_APOG,
            J2000,
            0,
            LunarCoverage::Unenforced,
            &stub_dpsi_rad,
            &stub_nut_rate,
        )
        .expect("mean pos");
        let raw = mean_lunar_apogee_state(J2000, LunarCoverage::Unenforced, false).expect("raw");
        assert_eq!(out.0, (raw.0 + 0.01).rem_euclid(360.0));
        assert_eq!((out.3, out.4, out.5), (0.0, 0.0, 0.0));
        // NONUT: mean ecliptic passes through untouched (SPEED included).
        let out = mean_lunar_point_state(
            MEAN_NODE,
            J2000,
            FLG_SPEED | FLG_NONUT,
            LunarCoverage::Unenforced,
            &stub_dpsi_rad,
            &stub_nut_rate,
        )
        .expect("nonut");
        assert_eq!(out, raw_for_node());
        // Family refusal + range escalation + non-finite refusal.
        assert_eq!(
            mean_lunar_point_state(
                TRUE_NODE,
                J2000,
                0,
                LunarCoverage::Unenforced,
                &stub_dpsi_rad,
                &stub_nut_rate
            ),
            Err(Error::UnknownBody { body_id: TRUE_NODE })
        );
        assert!(matches!(
            mean_lunar_point_state(
                MEAN_NODE,
                J2000,
                0,
                LunarCoverage::Window {
                    start_jd: J2000 + 1.0,
                    end_jd: J2000 + 2.0
                },
                &stub_dpsi_rad,
                &stub_nut_rate
            ),
            Err(Error::EphemerisRange { .. })
        ));
        assert!(matches!(
            mean_lunar_point_state(
                MEAN_NODE,
                f64::NAN,
                0,
                LunarCoverage::Unenforced,
                &stub_dpsi_rad,
                &stub_nut_rate
            ),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    fn raw_for_node() -> (f64, f64, f64, f64, f64, f64) {
        mean_lunar_node_state(J2000, LunarCoverage::Unenforced, false).expect("raw")
    }

    #[test]
    fn true_dispatch_strips_nutation_with_per_point_steps() {
        use crate::constants::{FLG_NONUT, FLG_SPEED, FLG_SPEED3, FLG_TOPOCTR};
        // Default SPEED: 0.05-day stencil for the node.
        let out = true_lunar_point_state(
            TRUE_NODE,
            J2000,
            FLG_SPEED,
            &stub_true_state,
            &stub_dpsi_rad,
            &stub_nut_rate,
        )
        .expect("node speed");
        let (r, v) = stub_true_state(J2000).expect("stub");
        let central = true_lunar_node_from_state(r, v).expect("central");
        assert_eq!((out.0, out.1, out.2), central);
        // The stub state is epoch-constant, so side samples agree and
        // speeds are exactly zero (stencil identity, not a fallback).
        assert_eq!((out.3, out.4, out.5), (0.0, 0.0, 0.0));
        // NONUT strips dpsi from the central longitude and the rate
        // from dlon (zero here, so the strip is directly visible).
        let out = true_lunar_point_state(
            TRUE_NODE,
            J2000,
            FLG_SPEED | FLG_NONUT,
            &stub_true_state,
            &stub_dpsi_rad,
            &stub_nut_rate,
        )
        .expect("node nonut");
        assert_eq!(out.0, (central.0 - 0.01).rem_euclid(360.0));
        assert_eq!(out.3, -2e-6);
        // OSCU_APOG uses the 0.002-day ordinary step: a provider that
        // fails off-epoch exercises the zero-speed boundary fallback.
        let edge = |jd: f64| -> Result<EclipticMoonState, Error> {
            if jd == J2000 {
                stub_true_state(jd)
            } else {
                Err(Error::EphemerisRange {
                    requested_jd: jd,
                    start_jd: J2000,
                    end_jd: J2000,
                    body: Some(OSCU_APOG),
                })
            }
        };
        let out = true_lunar_point_state(
            OSCU_APOG,
            J2000,
            FLG_SPEED,
            &edge,
            &stub_dpsi_rad,
            &stub_nut_rate,
        )
        .expect("apse edge");
        assert_eq!((out.3, out.4, out.5), (0.0, 0.0, 0.0));
        // But a central miss propagates fail-closed, never zeros.
        assert!(matches!(
            true_lunar_point_state(
                OSCU_APOG,
                J2000 + 5.0,
                FLG_SPEED,
                &edge,
                &stub_dpsi_rad,
                &stub_nut_rate
            ),
            Err(Error::EphemerisRange { .. })
        ));
        // SPEED 3 selects the fine step; verify via the selector wiring
        // (the stub is epoch-constant, so all stencils read zero).
        for flags in [FLG_SPEED | FLG_SPEED3, FLG_SPEED | FLG_TOPOCTR] {
            let out = true_lunar_point_state(
                TRUE_NODE,
                J2000,
                flags,
                &stub_true_state,
                &stub_dpsi_rad,
                &stub_nut_rate,
            )
            .expect("fine");
            assert_eq!((out.3, out.4, out.5), (0.0, 0.0, 0.0));
        }
        assert_eq!(
            select_true_point_half_step(true, TRUE_NODE_HALF_STEP_DAYS).expect("fine"),
            FINE_POINT_HALF_STEP_DAYS
        );
        // No SPEED: raw position, zeroed speeds.
        let out = true_lunar_point_state(
            TRUE_NODE,
            J2000,
            0,
            &stub_true_state,
            &stub_dpsi_rad,
            &stub_nut_rate,
        )
        .expect("pos");
        assert_eq!((out.0, out.1, out.2), central);
        assert_eq!((out.3, out.4, out.5), (0.0, 0.0, 0.0));
        // Family refusal.
        assert_eq!(
            true_lunar_point_state(
                MEAN_NODE,
                J2000,
                0,
                &stub_true_state,
                &stub_dpsi_rad,
                &stub_nut_rate
            ),
            Err(Error::UnknownBody { body_id: MEAN_NODE })
        );
    }

    #[test]
    fn interp_dispatch_keeps_verbatim_speed3_split() {
        use crate::constants::{FLG_NONUT, FLG_SPEED, FLG_SPEED3, FLG_TOPOCTR};
        assert_eq!(select_interp_speed_mode(0), (false, false));
        assert_eq!(select_interp_speed_mode(FLG_SPEED), (false, false));
        assert_eq!(
            select_interp_speed_mode(FLG_SPEED | FLG_TOPOCTR),
            (false, true)
        );
        assert_eq!(select_interp_speed_mode(FLG_SPEED3), (true, false));
        assert_eq!(
            select_interp_speed_mode(FLG_SPEED | FLG_SPEED3),
            (true, true)
        );
        // Default (no SPEED): central position, zeroed speeds.
        let out = interp_lunar_point_state(INTP_APOG, J2000, 0, &stub_dpsi_rad, &stub_nut_rate)
            .expect("intp pos");
        let central = calc_interpolated_apogee(J2000).expect("central");
        assert_eq!((out.0, out.1, out.2), central);
        assert_eq!((out.3, out.4, out.5), (0.0, 0.0, 0.0));
        // SPEED: stencil speeds carried through on the true ecliptic.
        let out =
            interp_lunar_point_state(INTP_APOG, J2000, FLG_SPEED, &stub_dpsi_rad, &stub_nut_rate)
                .expect("intp speed");
        let full = calc_interpolated_apse_state(J2000, INTP_APOG, false).expect("full");
        assert_eq!(out, full);
        // NONUT strips dpsi (subtracted, unlike the mean add) and the
        // 0.001-day rate from dlon.
        let out = interp_lunar_point_state(
            INTP_APOG,
            J2000,
            FLG_SPEED | FLG_NONUT,
            &stub_dpsi_rad,
            &stub_nut_rate,
        )
        .expect("intp nonut");
        assert_eq!(out.0, (full.0 - 0.01).rem_euclid(360.0));
        assert_eq!(out.3, full.3 - 2e-6);
        // SPEED+TOPOCTR takes the speed 3-velocity arm (identical values:
        // the stencil spelling is already the smooth rate).
        let out = interp_lunar_point_state(
            INTP_PERG,
            J2000,
            FLG_SPEED | FLG_TOPOCTR,
            &stub_dpsi_rad,
            &stub_nut_rate,
        )
        .expect("topo");
        let full = calc_interpolated_apse_state(J2000, INTP_PERG, false).expect("full");
        assert_eq!((out.0, out.1, out.2), (full.0, full.1, full.2));
        assert_eq!((out.3, out.4, out.5), (full.3, full.4, full.5));
        // Family refusal + non-finite refusal.
        assert_eq!(
            interp_lunar_point_state(MEAN_NODE, J2000, 0, &stub_dpsi_rad, &stub_nut_rate),
            Err(Error::UnknownBody { body_id: MEAN_NODE })
        );
        assert!(matches!(
            interp_lunar_point_state(INTP_APOG, f64::INFINITY, 0, &stub_dpsi_rad, &stub_nut_rate),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[test]
    fn lilith_diagnostic_compares_identical_alias_spellings() {
        // Both spellings are the maintained alias, so every difference
        // is exactly zero by construction; the struct pins the model
        // dict shape for the serving binding.
        let cmp = compare_true_lilith_methods_via(&stub_true_state, J2000).expect("diag");
        assert_eq!(cmp.eccentricity_vector, cmp.orbital_elements);
        assert_eq!((cmp.lon_diff, cmp.lat_diff, cmp.dist_diff), (0.0, 0.0, 0.0));
        let (r, v) = stub_true_state(J2000).expect("stub");
        assert_eq!(
            cmp.eccentricity_vector,
            true_lilith_from_state(r, v).expect("geo")
        );
        assert!(matches!(
            compare_true_lilith_methods_via(&stub_true_state, f64::NAN),
            Err(Error::NonFiniteArgument { .. })
        ));
        let failing = |_: f64| -> Result<EclipticMoonState, Error> {
            Err(Error::UnknownBody { body_id: 99 })
        };
        assert_eq!(
            compare_true_lilith_methods_via(&failing, J2000),
            Err(Error::UnknownBody { body_id: 99 })
        );
    }

    // -- model-regeneration tooling ------------------------------

    #[test]
    fn basis_inventories_match_the_release_selection() {
        // Counts construction: apogee
        // 14 + 3 + 2 * 5 + 6 + 1 = 34, perigee 30 + 3 + 2 * 13 + 6 + 1
        // = 66.
        let apo = apogee_basis_terms();
        let per = perigee_basis_terms();
        assert_eq!(apo.len(), 34);
        assert_eq!(per.len(), 66);
        // Ladder heads/tails, solar triple, first sideband pair,
        // inclination tail.
        assert_eq!(
            &apo[..3],
            &[(1, 0, -1, 0, 0), (2, 0, -2, 0, 0), (3, 0, -3, 0, 0)]
        );
        assert_eq!(apo[13], (14, 0, -14, 0, 0));
        assert_eq!(
            &apo[14..17],
            &[(0, 1, 0, 0, 1), (1, 1, -1, 0, 1), (1, -1, -1, 0, 1)]
        );
        assert_eq!(&apo[17..19], &[(2, -1, -2, 0, 1), (2, 1, -2, 0, 1)]);
        assert_eq!(apo[33], (0, 0, -2, 2, 0));
        assert_eq!(per[29], (30, 0, -30, 0, 0));
        assert_eq!(
            &per[30..33],
            &[(0, 1, 0, 0, 1), (1, 1, -1, 0, 1), (1, -1, -1, 0, 1)]
        );
        assert_eq!(per[65], (0, 0, -2, 2, 0));
        // Full-inventory cross-check against the checked-in generated
        // tables: stripping the fitted coefficients must reproduce the
        // generator basis exactly (no hand edit on either side).
        use crate::lunar_apse_tables as tables;

        let stripped: Vec<ApseBasisTerm> = tables::APOGEE_TRIG_TERMS
            .iter()
            .map(|&(_, c_d, c_m, c_mp, c_f, e)| (c_d, c_m, c_mp, c_f, e))
            .collect();
        assert_eq!(stripped, apo);
        let stripped: Vec<ApseBasisTerm> = tables::PERIGEE_TRIG_TERMS
            .iter()
            .map(|&(_, c_d, c_m, c_mp, c_f, e)| (c_d, c_m, c_mp, c_f, e))
            .collect();
        assert_eq!(stripped, per);
    }

    #[test]
    fn design_row_follows_the_generator_association() {
        let terms = apogee_basis_terms();
        // Zero angles: every sine is 0, so the row is all zeros
        // regardless of the eccentricity powers.
        let row = apse_design_row(0.0, 0.0, 0.0, 0.0, 1.03, &terms).expect("finite");
        assert_eq!(row.len(), 34);
        assert!(row.iter().all(|&v| v == 0.0));
        // Powers of two scale exactly: with E = 2 the e=1 solar term
        // (0,1,0,0,1) is exactly 2*sin(M) and the e=2 term
        // (0,2,0,0,2) exactly 4*sin(2 M); the e=0 ladder term is the
        // bare sine.
        let (d, m, mp, f) = (0.7, -1.3, 2.1, 0.2);
        let row = apse_design_row(d, m, mp, f, 2.0, &terms).expect("finite");
        assert_eq!(row[0], (d - mp).sin());
        let solar = row[14];
        assert_eq!(solar, 2.0 * m.sin());
        let second = row[27];
        assert_eq!(second, 4.0 * (2.0 * m).sin());
        // Empty term slices give empty rows; non-finite inputs are
        // typed refusals.
        assert_eq!(
            apse_design_row(d, m, mp, f, 1.0, &[]).expect("empty"),
            Vec::<f64>::new()
        );
        assert!(matches!(
            apse_design_row(f64::NAN, m, mp, f, 1.0, &terms),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            apse_design_row(d, m, mp, f, f64::INFINITY, &terms),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[test]
    fn spline_reproduces_linear_data_and_knot_values() {
        // Linear data: the natural spline is the line itself, exactly
        // (all second differences are exact zeros on small integers).
        let spline = NaturalCubicSpline::try_new(&[0.0, 1.0, 2.0, 4.0], &[1.0, 3.0, 5.0, 9.0])
            .expect("linear");
        for q in [0.0, 0.5, 1.0, 2.5, 4.0, -1.0e6] {
            assert_eq!(spline.eval(q).expect("q"), 2.0 * q + 1.0, "q = {q}");
        }
        // eval_many shares one construction; knots evaluate through
        // the left interval, so knot identity holds up to
        // the interval Horner rounding, not exactly.
        let xs: [f64; 5] = [0.0, 0.5, 1.5, 3.0, 3.5];
        let ys: Vec<f64> = xs.iter().map(|x| (*x).sin()).collect();
        let spline = NaturalCubicSpline::try_new(&xs, &ys).expect("sine");
        let back = spline.eval_many(&xs).expect("knots");
        for (got, want) in back.iter().zip(ys.iter()) {
            assert!(
                (got - want).abs() <= 8.0 * f64::EPSILON * want.abs().max(1.0),
                "knot eval {got} vs {want}"
            );
        }
        assert_eq!(spline.knots(), xs);
        assert_eq!(spline.a().len(), 4);
        // Out-of-knot queries extrapolate the edge interval instead of
        // failing.
        assert!(spline.eval(-2.0).expect("left").is_finite());
        assert!(spline.eval(9.0).expect("right").is_finite());
    }

    #[test]
    fn spline_construction_refuses_malformed_knots() {
        assert!(matches!(
            NaturalCubicSpline::try_new(&[0.0], &[1.0]),
            Err(Error::InvalidGeneratorInput { .. })
        ));
        assert!(matches!(
            NaturalCubicSpline::try_new(&[], &[]),
            Err(Error::InvalidGeneratorInput { .. })
        ));
        assert!(matches!(
            NaturalCubicSpline::try_new(&[0.0, 1.0], &[1.0]),
            Err(Error::InvalidGeneratorInput { .. })
        ));
        assert!(matches!(
            NaturalCubicSpline::try_new(&[0.0, 1.0, 1.0], &[0.0, 0.0, 0.0]),
            Err(Error::InvalidGeneratorInput { .. })
        ));
        assert!(matches!(
            NaturalCubicSpline::try_new(&[1.0, 0.0], &[0.0, 0.0]),
            Err(Error::InvalidGeneratorInput { .. })
        ));
        assert!(matches!(
            NaturalCubicSpline::try_new(&[0.0, f64::NAN], &[0.0, 0.0]),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            NaturalCubicSpline::try_new(&[0.0, 1.0], &[0.0, f64::INFINITY]),
            Err(Error::NonFiniteArgument { .. })
        ));
        let spline = NaturalCubicSpline::try_new(&[0.0, 1.0], &[0.0, 1.0]).expect("ok");
        assert!(matches!(
            spline.eval(f64::NAN),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            spline.eval_many(&[0.5, f64::NAN]),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[test]
    fn grid_assembly_tapers_edges_and_hits_passages() {
        // Constant unit residual: the spline is exactly 1.0 everywhere,
        // so in-span corrections are exactly 3600 arcsec and the table
        // error is exactly zero.
        let passages = [10.0, 20.0, 30.0];
        let residuals = [1.0, 1.0, 1.0];
        let assembly =
            assemble_correction_grid(&passages, &residuals, 0.0, 5.0, 9, 365.25).expect("grid");
        assert_eq!(assembly.corrections_arcsec().len(), 9);
        assert_eq!(assembly.corrections_arcsec()[2], 3600.0);
        assert_eq!(assembly.corrections_arcsec()[4], 3600.0);
        // Taper arms use the clipped linear weight.
        let weight = 1.0 - 10.0 / 365.25;
        assert_eq!(assembly.corrections_arcsec()[0], 1.0 * weight * 3600.0);
        assert_eq!(assembly.corrections_arcsec()[8], 1.0 * weight * 3600.0);
        assert_eq!(assembly.table_error_deg(), &[0.0, 0.0, 0.0]);
        assert_eq!(assembly.table_rms_arcsec(), 0.0);
        assert_eq!(assembly.table_max_arcsec(), 0.0);
        // Beyond the taper the correction is exactly 0.0.
        let wide = assemble_correction_grid(&passages, &residuals, -1000.0, 5.0, 801, 365.25)
            .expect("wide");
        assert_eq!(wide.corrections_arcsec()[0], 0.0);
        assert_eq!(wide.corrections_arcsec()[800], 0.0);
        // Grid knots coinciding with passages reproduce them: the
        // table error stays exactly zero on non-constant data too.
        let passages = [0.0, 10.0, 20.0];
        let residuals = [0.0, 1.0, 0.0];
        let assembly =
            assemble_correction_grid(&passages, &residuals, 0.0, 10.0, 3, 365.25).expect("aligned");
        assert_eq!(assembly.table_error_deg(), &[0.0, 0.0, 0.0]);
        assert_eq!(assembly.table_rms_arcsec(), 0.0);
    }

    #[test]
    fn grid_assembly_refuses_malformed_inputs() {
        let ok = ([10.0, 20.0, 30.0], [1.0, 1.0, 1.0]);
        assert!(matches!(
            assemble_correction_grid(&[10.0], &[1.0], 0.0, 5.0, 9, 365.25),
            Err(Error::InvalidGeneratorInput { .. })
        ));
        assert!(matches!(
            assemble_correction_grid(&ok.0, &[1.0], 0.0, 5.0, 9, 365.25),
            Err(Error::InvalidGeneratorInput { .. })
        ));
        assert!(matches!(
            assemble_correction_grid(&ok.0, &ok.1, 0.0, 5.0, 0, 365.25),
            Err(Error::InvalidGeneratorInput { .. })
        ));
        assert!(matches!(
            assemble_correction_grid(&ok.0, &ok.1, 0.0, 0.0, 9, 365.25),
            Err(Error::InvalidGeneratorInput { .. })
        ));
        assert!(matches!(
            assemble_correction_grid(&ok.0, &ok.1, 0.0, -5.0, 9, 365.25),
            Err(Error::InvalidGeneratorInput { .. })
        ));
        assert!(matches!(
            assemble_correction_grid(&ok.0, &ok.1, 0.0, 5.0, 9, 0.0),
            Err(Error::InvalidGeneratorInput { .. })
        ));
        assert!(matches!(
            assemble_correction_grid(&[10.0, f64::NAN, 30.0], &ok.1, 0.0, 5.0, 9, 365.25),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[test]
    fn passage_scan_pins_release_extremum_predicates() {
        // Cosine samples at 0.5 spacing over [0, 10]: interior minimum
        // near pi (i = 6), interior maximum near 2 pi (i = 13),
        // interior minimum near 3 pi (i = 19). The global-maximum
        // endpoint i = 0 is never reported.
        let distances: Vec<f64> = (0..21).map(|i| (i as f64 * 0.5).cos()).collect();
        let (maxima, minima) = scan_passage_candidates(&distances);
        assert_eq!(maxima, vec![13]);
        assert_eq!(minima, vec![6, 19]);
        // `>=` / `>` asymmetry: a flat top reports its right
        // edge as the maximum, never the left.
        let (maxima, minima) = scan_passage_candidates(&[1.0, 2.0, 2.0, 1.0]);
        assert_eq!(maxima, vec![2]);
        assert!(minima.is_empty());
        // Flat bottom mirrors the rule.
        let (maxima, minima) = scan_passage_candidates(&[2.0, 1.0, 1.0, 2.0]);
        assert!(maxima.is_empty());
        assert_eq!(minima, vec![2]);
        // Short series hold no inner point; NaN samples never win a
        // predicate on either side.
        assert_eq!(
            scan_passage_candidates(&[1.0, 2.0]),
            (Vec::new(), Vec::new())
        );
        assert_eq!(
            scan_passage_candidates(&[1.0, f64::NAN, 1.0]),
            (Vec::new(), Vec::new())
        );
    }

    #[test]
    fn passage_refinement_finds_parabola_vertices() {
        // A parabola refines to its vertex in one round; the assert
        // allows rounding slack over the remaining five.
        let bowl = |t: f64| -> Result<f64, Error> { Ok((t - 5.0) * (t - 5.0)) };
        let found = refine_passage(5.3, 0.5, &bowl).expect("min");
        assert!((found - 5.0).abs() < 1e-12, "found = {found}");
        let cap = |t: f64| -> Result<f64, Error> { Ok(-(t - 5.0) * (t - 5.0)) };
        let found = refine_passage(4.6, 0.5, &cap).expect("max");
        assert!((found - 5.0).abs() < 1e-12, "found = {found}");
        // Zero curvature leaves the seed in place; provider failures propagate fail-closed.
        assert_eq!(refine_passage(7.0, 0.5, &|_| Ok(3.0)).expect("flat"), 7.0);
        assert_eq!(
            refine_passage(7.0, 0.5, &|_| Err(Error::UnknownBody { body_id: 1 })),
            Err(Error::UnknownBody { body_id: 1 })
        );
        assert!(matches!(
            refine_passage(f64::NAN, 0.5, &bowl),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            refine_passage(5.0, 0.0, &bowl),
            Err(Error::InvalidGeneratorInput { .. })
        ));
        assert!(matches!(
            refine_passage(5.0, -0.5, &bowl),
            Err(Error::InvalidGeneratorInput { .. })
        ));
    }

    #[test]
    fn checked_in_corrections_fit_milliarcsecond_storage() {
        use crate::lunar_apse_tables as tables;

        // Every checked-in correction value round-trips through the
        // the model `{v:.3f}` milliarcsecond storage resolution: parsing
        // the 3-decimal rendering reproduces the value. (Exact halfway
        // ties between adjacent 3-decimal strings are unrepresentable
        // in binary64, so the two correctly-rounded formatters cannot
        // disagree.)
        for words in [tables::APOGEE_CORRECTIONS, tables::PERIGEE_CORRECTIONS] {
            for &word in words.iter() {
                let text = format!("{word:.3}");
                assert_eq!(text.parse::<f64>().expect("3dp"), word, "word {word}");
            }
        }
    }
}
