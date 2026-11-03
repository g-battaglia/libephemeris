//! Serving implementation for the parent module.
//! Public entry points and types are exposed through the parent API.

use super::*;

/// Sine of a degree argument.
///
/// Formulas in this section follow the literature in degrees
/// (semi-arcs, meridian distances, poles); these six helpers keep
/// degree arguments and returns so call sites read like the published
/// constructions, with the inverse pair clamped to the principal domain
/// against few-ULP overshoot. Scalar on both sides.
#[inline]
pub(super) fn sin_deg(x_deg: f64) -> f64 {
    x_deg.to_radians().sin()
}

/// Cosine of a degree argument.
#[inline]
pub(super) fn cos_deg(x_deg: f64) -> f64 {
    x_deg.to_radians().cos()
}

/// Tangent of a degree argument.
#[inline]
pub(super) fn tan_deg(x_deg: f64) -> f64 {
    x_deg.to_radians().tan()
}

/// Principal arctangent in degrees.
#[inline]
pub(super) fn atan_deg(x: f64) -> f64 {
    x.atan().to_degrees()
}

/// Arcsine in degrees after the roundoff-only unit-domain clamp
#[inline]
pub(super) fn asin_deg(x: f64) -> f64 {
    clamp_unit(x).asin().to_degrees()
}

/// Arccosine in degrees after the roundoff-only unit-domain clamp
#[inline]
pub(super) fn acos_deg(x: f64) -> f64 {
    clamp_unit(x).acos().to_degrees()
}

/// Tilt a spherical position into a frame rotated about the x-axis
/// (Meeus Ch. 13): rotation by `-tilt`, `atan2` longitude, clamped-`asin`
/// latitude, floored `% 360` normalization, radius dropped.
///
/// The rotation is spelled here in the same operations and association
/// rather than routing through [`crate::positions::cotrans`], whose
/// `r`-division and re-associated latitude sum can round the last ulp
/// differently, and whose finiteness gate would refuse inputs the
/// release propagates. `lon_deg`/`lat_deg`/`tilt_deg` are degrees in;
/// returns `(lon_out, lat_out)` in degrees, longitude in `[0, 360)`.
/// Scalar on both sides.
pub(super) fn rotate_frame(lon_deg: f64, lat_deg: f64, tilt_deg: f64) -> (f64, f64) {
    let lon_rad = lon_deg.to_radians();
    let lat_rad = lat_deg.to_radians();
    let eps_rad = (-tilt_deg).to_radians();
    let (sin_eps, cos_eps) = sin_cos(eps_rad);
    let (sin_lat, cos_lat) = sin_cos(lat_rad);
    let (sin_lon, cos_lon) = sin_cos(lon_rad);
    let sin_new_lat = clamp_unit(sin_lat * cos_eps + cos_lat * sin_eps * sin_lon);
    let new_lat = sin_new_lat.asin().to_degrees();
    let y = cos_lat * sin_lon * cos_eps - sin_lat * sin_eps;
    let x = cos_lat * cos_lon;
    let new_lon = y.atan2(x).to_degrees().rem_euclid(360.0);
    (new_lon, new_lat)
}

/// Ascendant longitude anchored to the rising intersection (Smart/Meeus horizon construction).
///
/// The closed-form ascendant silently picks the horizon-ecliptic
/// intersection nearer the meridian, which inside the polar circles
/// can be the setting point — a six-house error for Asc-anchored
/// systems; the MC-declination altitude test
/// (`δ_MC = atan(sin ARMC·tan ε)`, northern altitude `90 − φ + δ_MC`,
/// southern `−(−90 − φ + δ_MC)`, negative altitude → `+180°` flip)
/// restores the eastern (rising) point. Strict `< 0.0` / `> 0.0` edges
/// on both sides. Arguments are degrees in (closed-form Ascendant,
/// ARMC, obliquity, latitude), eastern-horizon Ascendant out in
/// `[0, 360)`; scalar on both sides.
pub fn ascendant_on_eastern_horizon(
    asc_deg: f64,
    armc_deg: f64,
    eps_deg: f64,
    geolat_deg: f64,
) -> f64 {
    let delta_mc = atan_deg(sin_deg(armc_deg) * tan_deg(eps_deg));
    let mut asc = asc_deg;
    if geolat_deg >= 0.0 && 90.0 - geolat_deg + delta_mc < 0.0 {
        asc = (asc + 180.0).rem_euclid(360.0);
    }
    if geolat_deg < 0.0 && -90.0 - geolat_deg + delta_mc > 0.0 {
        asc = (asc + 180.0).rem_euclid(360.0);
    }
    asc
}

/// Signed above-horizon indicator: `tan(dec) * tan(lat) + cos(md)` from the altitude
/// formula divided through by `cos(lat) * cos(dec)`. Non-negative
/// means on or above the horizon, valid where the divided-through
/// divisor is positive (`|lat|`, `|dec|` below 90°).
pub(super) fn above_horizon_test(dec_deg: f64, geolat_deg: f64, merid_dist_deg: f64) -> f64 {
    tan_deg(dec_deg) * tan_deg(geolat_deg) + cos_deg(merid_dist_deg)
}

/// Ascensional difference saturated at the circumpolar limit (Meeus Ch. 15): `AD = asin(tan φ·tan δ)`, the
/// rising-time correction; unit-excess saturates to ±90° so the
/// semi-arc becomes 180/0°. Arguments are degrees in (declination,
/// latitude), saturated AD degrees out; scalar on both sides.
pub(super) fn asc_diff_saturated(dec_deg: f64, geolat_deg: f64) -> f64 {
    let s = tan_deg(dec_deg) * tan_deg(geolat_deg);
    if s >= 1.0 {
        90.0
    } else if s <= -1.0 {
        -90.0
    } else {
        asin_deg(s)
    }
}

// JD acquisition, sidereal post-processing, ex2 walks, fallbacks, context.

/// Validate a JD-entry coordinate pair: the JD-entry lat/lon gate —
/// latitude in `[-90, 90]`, longitude in `[-180, 360]` accepting the
/// signed and east-positive conventions, both finite with NaN explicitly
/// refused (a bare range test lets NaN through). Mirrors the release
/// order (latitude first) and the release `CoordinateError` with
/// name/value/bounds + `func_name` prefix → [`Error::InvalidObserver`]
/// contract; the Python message prose stays binding-side. Scalar gate
/// on both sides.
pub fn validate_houses_coordinates(lat_deg: f64, lon_deg: f64) -> Result<()> {
    validate_houses_armc_latitude(lat_deg)?;
    validate_ecliptic_longitude(lon_deg)
}

/// Longitude validation for the JD-entry gate: finite longitude in `[-180, 360]`
/// accepting both the signed convention (`[-180, 180]`, west negative)
/// and the east-positive convention (`[0, 360]`), NaN explicitly
/// refused — mapped to [`Error::InvalidObserver`] with
/// `field == "longitude"` (the `CoordinateError` → `InvalidObserver`
/// contract in [`crate::errors`], message prose binding-side). Scalar
/// gate on both sides.
pub fn validate_ecliptic_longitude(lon_deg: f64) -> Result<()> {
    if !lon_deg.is_finite() || lon_deg < -180.0 || lon_deg > 360.0 {
        return Err(Error::InvalidObserver {
            field: "longitude",
            value: lon_deg,
        });
    }
    Ok(())
}

/// Classify a latitude and obliquity in degrees for house calculation.
///
/// The extreme flag uses the project's `|latitude| >= 80°` warning convention.
/// The polar flag uses [`is_polar_circle`], preserving the exact input side
/// even when the displayed threshold rounds onto that latitude. The remaining
/// fields provide N/S hemisphere, rounded thresholds and the affected,
/// numerically unstable and all-latitude-stable house-system lists.
pub fn extreme_latitude_info(lat_deg: f64, obliquity_deg: f64) -> ExtremeLatitudeInfo {
    let polar_threshold = polar_latitude_threshold_deg(obliquity_deg);
    let is_polar = is_polar_circle(lat_deg, obliquity_deg);
    let is_extreme = is_extreme_latitude(lat_deg);
    ExtremeLatitudeInfo {
        latitude: lat_deg,
        is_extreme,
        is_polar_circle: is_polar,
        polar_threshold,
        extreme_threshold: EXTREME_LATITUDE_THRESHOLD_DEG,
        hemisphere: if lat_deg >= 0.0 { 'N' } else { 'S' },
        affected_systems: if is_polar {
            POLAR_AFFECTED_HSYS.to_vec()
        } else {
            Vec::new()
        },
        unstable_systems: if is_extreme {
            EXTREME_UNSTABLE_HSYS.to_vec()
        } else {
            Vec::new()
        },
        stable_systems: STABLE_ALL_LATITUDE_HSYS.to_vec(),
    }
}

/// Display name for fallback warning messages.
pub fn fallback_system_name(code: char) -> String {
    match code {
        'P' => "Placidus".to_string(),
        'K' => "Koch".to_string(),
        'G' => "Gauquelin".to_string(),
        'O' => "Porphyry".to_string(),
        'E' => "Equal".to_string(),
        'W' => "Whole Sign".to_string(),
        'R' => "Regiomontanus".to_string(),
        'C' => "Campanus".to_string(),
        'T' => "Topocentric".to_string(),
        'B' => "Alcabitius".to_string(),
        'M' => "Morinus".to_string(),
        'X' => "Meridian".to_string(),
        'H' => "Horizontal".to_string(),
        'V' => "Vehlow".to_string(),
        'U' => "Krusinski".to_string(),
        'F' => "Carter".to_string(),
        'Y' => "APC".to_string(),
        'N' => "Natural Gradient".to_string(),
        other => other.to_string(),
    }
}

/// Attempt one house system, with an explicitly selected fallback.
///
/// `sample(folded_selector)` returns owned cusp degrees and eight auxiliary
/// angles in the order on [`HouseFallbackOutcome`]. Latitude and informational
/// obliquity are degrees. Selectors are case-folded before sampling.
///
/// A successful primary result is checked by [`validate_cusps`] when `validate`
/// is true, at any latitude: a defect immediately invokes the fallback and
/// reports `used_fallback=true`. Otherwise an unstable extreme-latitude system
/// can retain its primary result with an advisory warning and `used_fallback=false`.
/// Independently, a primary error whose cause is [`Error::PolarCircle`] invokes
/// the fallback. Every other primary error propagates.
///
/// The fallback result is not recursively retried or validated here; its errors
/// propagate. The optional warning reports which branch was taken. No callback
/// reference or internal sample cache survives the call; returned cusp/angle
/// values are owned by the result. No ephemeris/clock is acquired except through
/// the supplied callback. Fallback geometry is a caller
/// choice, not a silent claim that the requested system was defined.
pub fn houses_fallback_via(
    sample: &dyn Fn(char) -> Result<(Vec<f64>, [f64; N_ASCMC])>,
    lat_deg: f64,
    hsys: char,
    fallback_hsys: char,
    validate: bool,
    info_obliquity_deg: f64,
) -> Result<HouseFallbackOutcome> {
    let code = fold_hsys_case(hsys);
    let fallback_code = fold_hsys_case(fallback_hsys);
    let primary_name = fallback_system_name(code);
    let fallback_name = fallback_system_name(fallback_code);
    let extreme = is_extreme_latitude(lat_deg);
    match sample(code) {
        Ok((cusps, ascmc)) => {
            if validate {
                if let Err(defect) = validate_cusps(&cusps) {
                    let (fallback_cusps, fallback_ascmc) = sample(fallback_code)?;
                    return Ok(HouseFallbackOutcome {
                        cusps: fallback_cusps,
                        ascmc: fallback_ascmc,
                        used_fallback: true,
                        warning: Some(format!(
                            "{primary_name} house system produced invalid cusps at latitude {:.2}° ({defect}). Using {fallback_name} as fallback.",
                            lat_deg.abs()
                        )),
                    });
                }
            }
            if extreme {
                let info = extreme_latitude_info(lat_deg, info_obliquity_deg);
                if info.unstable_systems.contains(&code) {
                    return Ok(HouseFallbackOutcome {
                        cusps,
                        ascmc,
                        used_fallback: false,
                        warning: Some(format!(
                            "{primary_name} house system may have reduced accuracy at extreme latitude {:.2}°{}. Consider using a stable system like Porphyry, Equal, or Whole Sign.",
                            lat_deg.abs(),
                            info.hemisphere
                        )),
                    });
                }
            }
            Ok(HouseFallbackOutcome {
                cusps,
                ascmc,
                used_fallback: false,
                warning: None,
            })
        }
        Err(error) if matches!(error.primary_error(), Error::PolarCircle { .. }) => {
            let Error::PolarCircle { threshold, .. } = error.primary_error() else {
                unreachable!("polar fallback guard selected the primary")
            };
            let (fallback_cusps, fallback_ascmc) = sample(fallback_code)?;
            Ok(HouseFallbackOutcome {
                cusps: fallback_cusps,
                ascmc: fallback_ascmc,
                used_fallback: true,
                warning: Some(format!(
                    "{primary_name} house system unavailable at latitude {:.2}° (polar circle threshold: {threshold:.2}°). Using {fallback_name} as fallback.",
                    lat_deg.abs()
                )),
            })
        }
        Err(other) => Err(other),
    }
}

/// Calculate ARMC-based house cusps with an explicit polar-safe fallback.
///
/// ARMC, geographic latitude, true obliquity and the Sunshine-only Sun
/// declination are degrees. The primary system is attempted first. Polar
/// failures or invalid cusps select `fallback_hsys`; other numerical failures
/// propagate. Output contains 12 cusps (36 for uppercase Gauquelin), eight
/// angles, a fallback flag and an optional informational message. Both systems
/// use the same native geometry, and no ephemeris or time-scale data are read.
#[allow(clippy::too_many_arguments)]
pub fn houses_armc_with_fallback(
    armc: f64,
    lat_deg: f64,
    eps_deg: f64,
    hsys: char,
    fallback_hsys: char,
    validate: bool,
    ascmc9_deg: f64,
) -> Result<HouseFallbackOutcome> {
    let sample = |system| super::houses_armc(armc, lat_deg, eps_deg, system, ascmc9_deg);
    houses_fallback_via(&sample, lat_deg, hsys, fallback_hsys, validate, eps_deg)
}

/// Whether a failed Sunshine Sun fetch degrades to `0.0`:
/// value-domain failures of the sealed Sun computation. Range misses
/// take the analytic fallback instead ([`sealed_sun_declination_deg`]);
/// data/config failures propagate, like the uncaught release arms.
pub(super) fn is_sun_fetch_degraded(error: &Error) -> bool {
    matches!(
        error,
        Error::NonFiniteArgument { .. }
            | Error::UnknownBody { .. }
            | Error::InvalidBody { .. }
            | Error::InvalidObserver { .. }
            | Error::ConvergenceFailed { .. }
    )
}

/// Sunshine Sun declination with the release fallback policy
/// (`houses` Sunshine fetch; same policy in
/// `_houses_fixed_epoch_sidereal` and `_houses_sidbit_projection`): the
/// fetched declination when the ephemeris covers the date, the analytic
/// Meeus value on a typed range miss, `0.0` on a degraded value-domain
/// failure, propagation otherwise.
///
/// Python→Rust error-set mapping: `EphemerisRangeError` →
/// `Error::EphemerisRange` (analytic fallback);
/// `(IndexError, TypeError, ValueError, CalculationError)` →
/// `NonFiniteArgument | UnknownBody | InvalidBody | InvalidObserver |
/// ConvergenceFailed` (the degraded-set predicate, `0.0` degraded arm);
/// every other arm uncaught and propagated. `tjdut` is JD UT and
/// `ephemeris_flags` the ephemeris-bit-masked flags in, declination
/// degrees out; the fetch and ΔT legs arrive injected. Scalar on both
/// sides.
pub fn sealed_sun_declination_deg(
    tjdut: f64,
    ephemeris_flags: i32,
    sun_fetch: &dyn Fn(f64, i32) -> Result<f64>,
    deltat_at: &dyn Fn(f64) -> Result<f64>,
) -> Result<f64> {
    match sun_fetch(tjdut, ephemeris_flags) {
        Ok(declination) => Ok(declination),
        Err(Error::EphemerisRange { .. }) => Ok(sun_declination_analytic(tjdut, deltat_at(tjdut)?)),
        Err(other) if is_sun_fetch_degraded(&other) => Ok(0.0),
        Err(other) => Err(other),
    }
}

/// JD-based house cusps and ASCMC: JD-based cusps
/// + ASCMC with reference-API compatibility — coordinate validation,
///   ARMC = `armc0 + lon`, the Sunshine Sun-declination fetch gated on the
///   folded `I`/`i` selector with the release ephemeris-bit mask, then
///   the [`houses_armc`] entry (whose prologue, MC-flip rule, APC/Sunshine
///   re-syncs, refusals and `G`/`g` shapes the JD entry shares) and the
///   `degnorm` output shaping.
///
/// `tjdut` is JD UT1, `lat_deg`/`lon_deg` signed geographic degrees,
/// `hsys` the folded selector, `iflag` the propagated ephemeris flags;
/// acquisition arrives through the injected providers. The entry owns
/// validation, the `armc0 + lon` composition, the masked Sun gate and
/// the shaping (shared output-shaping bank); the kernels take the composed angles as
/// given. Non-finite epochs refuse with [`Error::MathDomain`] before
/// any provider runs.
///
/// Returns Rust `Vec` + `[f64; 8]`: 12 cusps (36 sectors for `G`, first-12 prefix for
/// `g`) + 8 ASCMC.
pub fn houses_jd_via(
    tjdut: f64,
    lat_deg: f64,
    lon_deg: f64,
    hsys: char,
    iflag: i32,
    providers: &HousesJdProviders<'_>,
) -> Result<(Vec<f64>, [f64; N_ASCMC])> {
    houses_jd_with_epoch_policy(tjdut, lat_deg, lon_deg, hsys, iflag, providers, true)
        .map(|(cusps, angles)| (cusps.into_vec(), angles))
}

/// Consume the real model frame after coordinate admission. The frame owner
/// controls calendar/operator failures for arbitrary epoch words; no separate
/// finite-epoch refusal precedes that genuine acquisition.
pub(crate) fn houses_jd_via_model_clock(
    tjdut: f64,
    lat_deg: f64,
    lon_deg: f64,
    hsys: char,
    iflag: i32,
    providers: &HousesJdProviders<'_>,
) -> Result<(Vec<f64>, [f64; N_ASCMC])> {
    houses_jd_with_epoch_policy(tjdut, lat_deg, lon_deg, hsys, iflag, providers, false)
        .map(|(cusps, angles)| (cusps.into_vec(), angles))
}

/// Shared geometry after the selected epoch/clock admission policy.
fn houses_jd_with_epoch_policy(
    tjdut: f64,
    lat_deg: f64,
    lon_deg: f64,
    hsys: char,
    iflag: i32,
    providers: &HousesJdProviders<'_>,
    reject_nonfinite_epoch: bool,
) -> Result<(CuspBuffer, [f64; N_ASCMC])> {
    use crate::constants::{FLG_JPLEPH, FLG_SWIEPH};
    validate_houses_coordinates(lat_deg, lon_deg)?;
    let code = fold_hsys_case(hsys);
    if reject_nonfinite_epoch && !tjdut.is_finite() {
        return Err(Error::MathDomain { value: tjdut });
    }
    let (armc0_deg, eps_true_deg) = (providers.frame_at)(tjdut)?;
    let armc_deg = (armc0_deg + lon_deg).rem_euclid(360.0);
    let sun_dec_deg = if code == 'I' || code == 'i' {
        sealed_sun_declination_deg(
            tjdut,
            iflag & (FLG_JPLEPH | FLG_SWIEPH),
            providers.sun_declination_at,
            providers.deltat_days_at,
        )?
    } else {
        0.0
    };
    houses_armc_buffer(armc_deg, lat_deg, eps_true_deg, code, sun_dec_deg)
}

/// Passive x-axis rotation matching the release `_rx`/`_rot_x`
/// association (`y`/`z` by `+angle` in the `(1,1)`/`(2,2)` slots,
/// antisymmetric sign on `(1,2)`). Same-argument pairs go through
/// [`sin_cos`].
pub(super) fn rot_x_local(angle_rad: f64) -> [[f64; 3]; 3] {
    let (sin_a, cos_a) = sin_cos(angle_rad);
    [[1.0, 0.0, 0.0], [0.0, cos_a, sin_a], [0.0, -sin_a, cos_a]]
}

/// `houses_ex` output shaping ( tail): `degnorm` on every
/// cusp and angle, then the `FLG_RADIANS` position-only conversion
/// (`math.radians`, i.e. `x * (π / 180)`, speeds untouched — the ex2
/// stencil stays in the degree layer).
pub(super) fn shape_houses_ex(
    cusps_deg: CuspBuffer,
    ascmc_deg: [f64; N_ASCMC],
    flags: i32,
) -> (CuspBuffer, [f64; N_ASCMC]) {
    use crate::constants::FLG_RADIANS;
    let cusps: CuspBuffer = cusps_deg.into_iter().map(degnorm_deg).collect();
    let ascmc: [f64; N_ASCMC] = ascmc_deg.map(degnorm_deg);
    if (flags & FLG_RADIANS) != 0 {
        let to_radians = std::f64::consts::PI / 180.0;
        (
            cusps.into_iter().map(|c| c * to_radians).collect(),
            ascmc.map(|a| a * to_radians),
        )
    } else {
        (cusps, ascmc)
    }
}

/// Construct fixed-epoch sidereal houses on their defining mean ecliptic.
///
/// `tjdut` is a UT1 Julian day; `lat_deg` is geographic latitude in degrees,
/// `code` is the folded house-system selector, and `trop_ascmc` contains eight
/// tropical angles in degrees. The selected fixed-mode table supplies `t0` as
/// a TT Julian day and an ecliptic-longitude zero-point offset in degrees.
/// USER and Valens defining-epoch callbacks are not consumed by this path.
///
/// Observation TT is `tjdut + deltat_days_at(tjdut)`. The nutation callback
/// supplies radians at that TT. With `PN` the true-equator-of-date rotation,
/// `V(t0)` the precession rotation and `eps_t0` the mean obliquity in radians,
/// `M = PN @ V(t0)^T @ Rx(-eps_t0)` maps mean-ecliptic-of-t0 components into
/// the true equator of date. The plane node rebases ARMC; its inclination is
/// the construction obliquity. The anchor arc is oriented by observation
/// `TT - t0`, while reported ARMC is not reflected.
///
/// Only Sunshine `I`/`i` consumes Sun declination at the original `tjdut`,
/// with the ephemeris mask passed to [`sealed_sun_declination_deg`]. Its
/// range/analytic and degraded-error fallback rules remain in that helper.
/// `i` retains its distinct Makransky geometry; `N` has zero-Aries cusps.
/// Returns unshaped degree cusps and eight degree angles, with slot 2 equal
/// to node-rebased ARMC; the caller applies [`shape_houses_ex`]. Missing mode,
/// model/frame/callback and house-geometry errors propagate. A degenerate
/// plane node retains IEEE nonfinite arithmetic rather than a fitted plane.
#[allow(clippy::too_many_arguments)]
pub(super) fn fixed_epoch_houses(
    tjdut: f64,
    lat_deg: f64,
    code: char,
    trop_ascmc: &[f64; N_ASCMC],
    ephemeris_flags: i32,
    base_mode: i32,
    jd: &HousesJdProviders<'_>,
    sidereal: &HousesSiderealCallbacks<'_>,
) -> Result<(CuspBuffer, [f64; N_ASCMC])> {
    use crate::positions::{mat3_mul, mat3_vec3, vondrak_pn_matrix, vondrak_precession_matrix};
    use crate::sidereal::{fixed_epoch_lon_offset_deg, fixed_epoch_t0_jd};
    let t0 =
        fixed_epoch_t0_jd(base_mode).ok_or(Error::AyanamshaAnchorUnserved { mode: base_mode })?;
    let lon_offset = fixed_epoch_lon_offset_deg(base_mode);
    let jd_tt = tjdut + (jd.deltat_days_at)(tjdut)?;
    let (dpsi_rad, deps_rad) = (sidereal.nutation_at)(jd_tt)?;
    let (pn, _) = vondrak_pn_matrix(jd_tt, dpsi_rad, deps_rad, true)?;
    let v_t0 = vondrak_precession_matrix(t0, true)?;
    let eps_t0 = crate::erfa::obl06(t0)?.to_radians();
    let m = mat3_mul(mat3_mul(pn, transpose3_local(v_t0))?, rot_x_local(-eps_t0))?;
    // Third row of M: the t0-ecliptic pole in equator-of-date components
    // (`m.T @ [0, 0, 1]` picks the third row, not the third column).
    let z = [m[2][0], m[2][1], m[2][2]];
    // NaN passes through (no panic: the bounds are finite literals), so
    // a pathological provider degrades to NaN words instead of trapping.
    let z2 = z[2].clamp(-1.0, 1.0);
    let eps_p = z2.acos().to_degrees();
    let mut node = [z[1], -z[0], 0.0];
    let norm = (node[0] * node[0] + node[1] * node[1] + 0.0 * 0.0).sqrt();
    node[0] /= norm;
    node[1] /= norm;
    let lon_node = node[1].atan2(node[0]).to_degrees().rem_euclid(360.0);
    let n_eq = mat3_vec3(m, node)?;
    let ang_n = n_eq[1].atan2(n_eq[0]).to_degrees().rem_euclid(360.0);
    let armc_p = (trop_ascmc[ASCMC_ARMC] - ang_n).rem_euclid(360.0);
    // Compatibility contract: the anchor arc is oriented by the sign of
    // `(t − t0)`; the reported ARMC above is not reflected.
    let d_node = if lon_node <= 180.0 {
        lon_node
    } else {
        lon_node - 360.0
    };
    let lon_node_oriented = if jd_tt >= t0 {
        -d_node.abs()
    } else {
        d_node.abs()
    };
    let engine_hsys = code;
    let ascmc9 = if engine_hsys == 'I' || engine_hsys == 'i' {
        sealed_sun_declination_deg(
            tjdut,
            ephemeris_flags,
            jd.sun_declination_at,
            jd.deltat_days_at,
        )?
    } else {
        0.0
    };
    let (eng_cusps, eng_ascmc) = houses_armc_buffer(armc_p, lat_deg, eps_p, engine_hsys, ascmc9)?;
    let cusps = if code == 'N' {
        (0..N_CUSPS).map(|i| i as f64 * 30.0).collect()
    } else {
        eng_cusps
            .into_iter()
            .map(|c| (lon_node_oriented - lon_offset + c).rem_euclid(360.0))
            .collect()
    };
    let mut ascmc = *trop_ascmc;
    for i in [
        ASCMC_ASC,
        ASCMC_MC,
        ASCMC_VERTEX,
        ASCMC_EQUASC,
        ASCMC_COASC_KOCH,
        ASCMC_COASC_MUNKASEY,
        ASCMC_POLAR_ASC,
    ] {
        ascmc[i] = (lon_node_oriented - lon_offset + eng_ascmc[i]).rem_euclid(360.0);
    }
    ascmc[ASCMC_ARMC] = armc_p;
    Ok((cusps, ascmc))
}

/// Construct sidereal houses on a defining-ecliptic or invariable plane.
///
/// `tjdut` is a UT1 Julian day, `lat_deg` is geographic latitude in degrees,
/// `code` is the folded selector, and `trop_ascmc` contains eight tropical
/// angles in degrees. Observation TT is `tjdut + deltat_days_at(tjdut)`;
/// nutation is acquired in radians there. `SIDBIT_ECL_T0` takes precedence
/// over `SIDBIT_SSY_PLANE`. The defining `t0` used for the ECL plane is a TT
/// Julian day: USER or Valens callbacks are consumed only when selected by
/// that branch; other modes use their table epoch. SSY consumes neither
/// defining-epoch callback.
///
/// The zero-point mean-ayanamsha callback receives the scalar defining `t0`
/// through its UT entry, or the scalar J2000 epoch for SSY. This distinct
/// callback clock is retained independently of the TT plane definition and
/// configured observation-TT clock. With `base` mapping J2000 mean-ecliptic
/// components into the selected plane, `PN` the true-equator-of-date rotation
/// and `eps_J2000` the mean obliquity in radians,
/// `M = PN @ V(J2000)^T @ Rx(-eps_J2000) @ base^T` maps plane components into
/// the true equator of date. The node rebases ARMC and the plane inclination
/// supplies the construction obliquity.
///
/// Only Sunshine `I`/`i` consumes the original-UT Sun callback via
/// [`sealed_sun_declination_deg`], including its range/analytic and degraded
/// fallback rules. `i` keeps its distinct Makransky geometry; `N` has
/// zero-Aries cusps. Returns unshaped degree cusps and eight degree angles;
/// slot 2 is node-rebased ARMC. The caller applies [`shape_houses_ex`]. Missing
/// projection bits and selected frame, epoch, zero-point or geometry failures
/// propagate in acquisition order. A degenerate node retains IEEE nonfinite
/// arithmetic. On SSY, a NaN derived by the selected mean-zero-point model
/// propagates into cusps and ecliptic angles while the independent finite
/// construction still supplies ARMC. No constant plane or zero point
/// substitutes for a failed feed; callback errors retain acquisition order.
#[allow(clippy::too_many_arguments)]
pub(super) fn sidbit_houses(
    tjdut: f64,
    lat_deg: f64,
    code: char,
    trop_ascmc: &[f64; N_ASCMC],
    ephemeris_flags: i32,
    base_mode: i32,
    sid_bits: i32,
    jd: &HousesJdProviders<'_>,
    sidereal: &HousesSiderealCallbacks<'_>,
) -> Result<(CuspBuffer, [f64; N_ASCMC])> {
    use crate::constants::{SIDBIT_ECL_T0, SIDBIT_SSY_PLANE};
    use crate::positions::{mat3_mul, mat3_vec3, vondrak_pn_matrix, vondrak_precession_matrix};
    use crate::sidereal::{
        J2000_TT, ecliptic_of_t0_matrix, invariable_plane_matrix, ssy_plane_zero_point_deg_ieee,
    };
    let jd_tt = tjdut + (jd.deltat_days_at)(tjdut)?;
    let (dpsi_rad, deps_rad) = (sidereal.nutation_at)(jd_tt)?;
    let (pn, _) = vondrak_pn_matrix(jd_tt, dpsi_rad, deps_rad, true)?;
    let v_j2k = vondrak_precession_matrix(J2000_TT, true)?;
    let eps_j2000 = crate::erfa::obl06(J2000_TT)?.to_radians();
    let (base_mat, zero_point) = if (sid_bits & SIDBIT_ECL_T0) != 0 {
        let t0_jd = sidereal.defining_epoch_tt(base_mode)?;
        (
            ecliptic_of_t0_matrix(t0_jd)?,
            (sidereal.mean_ut_plain_at)(t0_jd)?,
        )
    } else if (sid_bits & SIDBIT_SSY_PLANE) != 0 {
        (
            invariable_plane_matrix(),
            ssy_plane_zero_point_deg_ieee((sidereal.mean_ut_plain_at)(J2000_TT)?)?,
        )
    } else {
        return Err(Error::AyanamshaAnchorUnserved { mode: base_mode });
    };
    let m = mat3_mul(
        mat3_mul(
            mat3_mul(pn, transpose3_local(v_j2k))?,
            rot_x_local(-eps_j2000),
        )?,
        transpose3_local(base_mat),
    )?;
    // Third row of M (`m.T @ [0, 0, 1]`), like the fixed-epoch arm.
    let z = [m[2][0], m[2][1], m[2][2]];
    // NaN passes through (no panic: the bounds are finite literals), so
    // a pathological provider degrades to NaN words instead of trapping.
    let z2 = z[2].clamp(-1.0, 1.0);
    let eps_p = z2.acos().to_degrees();
    let mut node = [z[1], -z[0], 0.0];
    let norm = (node[0] * node[0] + node[1] * node[1] + 0.0 * 0.0).sqrt();
    node[0] /= norm;
    node[1] /= norm;
    let lon_node = node[1].atan2(node[0]).to_degrees().rem_euclid(360.0);
    let n_eq = mat3_vec3(m, node)?;
    let ang_n = n_eq[1].atan2(n_eq[0]).to_degrees().rem_euclid(360.0);
    let armc_p = (trop_ascmc[ASCMC_ARMC] - ang_n).rem_euclid(360.0);
    let engine_hsys = code;
    let ascmc9 = if engine_hsys == 'I' || engine_hsys == 'i' {
        sealed_sun_declination_deg(
            tjdut,
            ephemeris_flags,
            jd.sun_declination_at,
            jd.deltat_days_at,
        )?
    } else {
        0.0
    };
    let (eng_cusps, eng_ascmc) = houses_armc_buffer(armc_p, lat_deg, eps_p, engine_hsys, ascmc9)?;
    let off = lon_node - zero_point;
    let cusps = if code == 'N' {
        (0..N_CUSPS).map(|i| i as f64 * 30.0).collect()
    } else {
        eng_cusps
            .into_iter()
            .map(|c| (off + c).rem_euclid(360.0))
            .collect()
    };
    let mut ascmc = *trop_ascmc;
    for i in [
        ASCMC_ASC,
        ASCMC_MC,
        ASCMC_VERTEX,
        ASCMC_EQUASC,
        ASCMC_COASC_KOCH,
        ASCMC_COASC_MUNKASEY,
        ASCMC_POLAR_ASC,
    ] {
        ascmc[i] = (off + eng_ascmc[i]).rem_euclid(360.0);
    }
    ascmc[ASCMC_ARMC] = armc_p;
    Ok((cusps, ascmc))
}

/// Extended house calculation with sidereal support: the tropical pair (seeded from `I` for Sunshine `i`
/// in an ayanamsha zodiac, kept on the Makransky path for the
/// fixed-epoch modes), then the fixed-epoch frame request, the
/// `SIDBIT_ECL_T0`/`SIDBIT_SSY_PLANE` plane reconstruction (inert for
/// the suppress-listed live-direction modes, which take the ayanamsha
/// path), else the mean-equinox ayanamsha shift with the
/// `TRUEPOS`/`NOABERR`/`NOGDEFL` anchor toggles and the
/// Ascendant-based `W`/`E`/`A`/`V`/`N` arms (`W` plain floor, no 1e-9
/// nudge — sidereal path, unlike the tropical kernel).
/// `selection` is the context-threaded sidereal state (no globals).
///
/// `tjdut` is JD UT1, `lat_deg`/`lon_deg` signed geographic degrees,
/// `hsys` the folded selector, `flags` the release flag word;
/// `selection` the explicitly threaded sidereal state; JD/sidereal legs
/// injected. Returns owned cusp and eight-angle arrays, in the order on
/// [`HouseFallbackOutcome`], with 12 cusps or 36 for uppercase Gauquelin `G`.
/// Lowercase `g` retains the first 12 sectors. Angles are normalized degrees,
/// or radians under `FLG_RADIANS`; this entry returns no speeds. Coordinate,
/// construction, polar-domain and consumed callback errors propagate, subject
/// to the documented Sunshine Sun-declination policy.
#[allow(clippy::too_many_arguments)] // positional arity mirrors the release entries
pub fn houses_ex_via(
    tjdut: f64,
    lat_deg: f64,
    lon_deg: f64,
    hsys: char,
    flags: i32,
    selection: crate::sidereal::SiderealSelection,
    jd: &HousesJdProviders<'_>,
    sidereal: &HousesSiderealProviders<'_>,
) -> Result<(Vec<f64>, [f64; N_ASCMC])> {
    houses_ex_buffer_via(
        tjdut, lat_deg, lon_deg, hsys, flags, selection, jd, sidereal,
    )
    .map(|(cusps, angles)| (cusps.into_vec(), angles))
}

/// Same legacy provider-to-callback acquisition, with stack cusp storage.
#[allow(clippy::too_many_arguments)]
pub(crate) fn houses_ex_buffer_via(
    tjdut: f64,
    lat_deg: f64,
    lon_deg: f64,
    hsys: char,
    flags: i32,
    selection: crate::sidereal::SiderealSelection,
    jd: &HousesJdProviders<'_>,
    sidereal: &HousesSiderealProviders<'_>,
) -> Result<(CuspBuffer, [f64; N_ASCMC])> {
    let user_t0_tt_at = || Ok(sidereal.user_t0_tt);
    let valens_epoch_tt_at = || Ok(sidereal.valens_epoch_tt);
    let callbacks = HousesSiderealCallbacks {
        ayanamsha_ex_at: sidereal.ayanamsha_ex_at,
        mean_ut_plain_at: sidereal.mean_ut_plain_at,
        nutation_at: sidereal.nutation_at,
        user_t0_tt_at: &user_t0_tt_at,
        valens_epoch_tt_at: &valens_epoch_tt_at,
    };
    houses_ex_buffer_via_callbacks(
        tjdut, lat_deg, lon_deg, hsys, flags, selection, jd, &callbacks,
    )
}

/// Extended house calculation with lazily acquired defining epochs.
///
/// Geometry, flags, units, layouts and shaping follow [`houses_ex_via`].
/// Tropical seed acquisition runs first. Fixed-epoch and plain ayanamsha
/// routes consume only their selected feeds; the SSY plane never reads
/// USER/Valens plane epochs. ECL_T0 requests the selected defining epoch
/// after observation-frame acquisition and before its plane and zero point.
/// No image or body is required unless a consumed callback requires it.
/// Provider errors propagate before a result is returned.
#[allow(clippy::too_many_arguments)]
pub fn houses_ex_via_callbacks(
    tjdut: f64,
    lat_deg: f64,
    lon_deg: f64,
    hsys: char,
    flags: i32,
    selection: crate::sidereal::SiderealSelection,
    jd: &HousesJdProviders<'_>,
    sidereal: &HousesSiderealCallbacks<'_>,
) -> Result<(Vec<f64>, [f64; N_ASCMC])> {
    houses_ex_with_epoch_policy(
        tjdut, lat_deg, lon_deg, hsys, flags, selection, jd, sidereal, true,
    )
    .map(|(cusps, angles)| (cusps.into_vec(), angles))
}

/// Extended geometry whose tropical seed consumes its genuine model frame
/// before any independent epoch refusal. All subsequent sidereal branches,
/// callbacks, flag masks and result shaping are shared with the public entry.
#[allow(clippy::too_many_arguments)]
pub(crate) fn houses_ex_via_callbacks_model_clock(
    tjdut: f64,
    lat_deg: f64,
    lon_deg: f64,
    hsys: char,
    flags: i32,
    selection: crate::sidereal::SiderealSelection,
    jd: &HousesJdProviders<'_>,
    sidereal: &HousesSiderealCallbacks<'_>,
) -> Result<(Vec<f64>, [f64; N_ASCMC])> {
    houses_ex_with_epoch_policy(
        tjdut, lat_deg, lon_deg, hsys, flags, selection, jd, sidereal, false,
    )
    .map(|(cusps, angles)| (cusps.into_vec(), angles))
}

/// Same selected acquisition and projection as the public Vec entry, retaining
/// stack-owned cusps for a fixed-layout SDK result.
#[allow(clippy::too_many_arguments)]
pub(crate) fn houses_ex_buffer_via_callbacks(
    tjdut: f64,
    lat_deg: f64,
    lon_deg: f64,
    hsys: char,
    flags: i32,
    selection: crate::sidereal::SiderealSelection,
    jd: &HousesJdProviders<'_>,
    sidereal: &HousesSiderealCallbacks<'_>,
) -> Result<(CuspBuffer, [f64; N_ASCMC])> {
    houses_ex_with_epoch_policy(
        tjdut, lat_deg, lon_deg, hsys, flags, selection, jd, sidereal, true,
    )
}

/// The existing five-sample speed stencil, with identical evaluation and
/// arithmetic order and fixed cusp storage instead of intermediate Vecs.
pub(crate) fn houses_ex2_buffer_via<F>(
    sample: &F,
    tjdut: f64,
    folded_hsys: char,
) -> Result<(CuspBuffer, [f64; N_ASCMC], CuspBuffer, [f64; N_ASCMC])>
where
    F: Fn(f64) -> Result<(CuspBuffer, [f64; N_ASCMC])>,
{
    super::specialized_systems::houses_ex2_storage(sample, tjdut, folded_hsys)
}

/// Shared projection orchestration with only tropical epoch admission selected.
#[allow(clippy::too_many_arguments)]
fn houses_ex_with_epoch_policy(
    tjdut: f64,
    lat_deg: f64,
    lon_deg: f64,
    hsys: char,
    flags: i32,
    selection: crate::sidereal::SiderealSelection,
    jd: &HousesJdProviders<'_>,
    sidereal: &HousesSiderealCallbacks<'_>,
    reject_nonfinite_epoch: bool,
) -> Result<(CuspBuffer, [f64; N_ASCMC])> {
    use crate::constants::{
        FLG_NOABERR, FLG_NOGDEFL, FLG_SIDEREAL, FLG_TRUEPOS, SIDBIT_ECL_T0, SIDBIT_SSY_PLANE,
    };
    use crate::sidereal::{
        SIDBIT_PROJECTION_SUPPRESS_MODES, is_fixed_epoch_mode, strip_to_base_mode,
    };
    let code = fold_hsys_case(hsys);
    let base_mode = strip_to_base_mode(selection.mode());
    let sid_bits = selection.bits();
    let tropical_hsys =
        if (flags & FLG_SIDEREAL) != 0 && code == 'i' && !is_fixed_epoch_mode(base_mode) {
            'I'
        } else {
            hsys
        };
    let (trop_cusps, trop_ascmc) = houses_jd_with_epoch_policy(
        tjdut,
        lat_deg,
        lon_deg,
        tropical_hsys,
        flags,
        jd,
        reject_nonfinite_epoch,
    )?;
    if (flags & FLG_SIDEREAL) == 0 {
        return Ok(shape_houses_ex(trop_cusps, trop_ascmc, flags));
    }
    if is_fixed_epoch_mode(base_mode) {
        let (cusps, ascmc) = fixed_epoch_houses(
            tjdut,
            lat_deg,
            code,
            &trop_ascmc,
            flags & (crate::constants::FLG_JPLEPH | crate::constants::FLG_SWIEPH),
            base_mode,
            jd,
            sidereal,
        )?;
        return Ok(shape_houses_ex(cusps, ascmc, flags));
    }
    if (sid_bits & (SIDBIT_ECL_T0 | SIDBIT_SSY_PLANE)) != 0
        && !SIDBIT_PROJECTION_SUPPRESS_MODES.contains(&base_mode)
    {
        let (cusps, ascmc) = sidbit_houses(
            tjdut,
            lat_deg,
            code,
            &trop_ascmc,
            flags & (crate::constants::FLG_JPLEPH | crate::constants::FLG_SWIEPH),
            base_mode,
            sid_bits,
            jd,
            sidereal,
        )?;
        return Ok(shape_houses_ex(cusps, ascmc, flags));
    }
    let ayanamsa =
        (sidereal.ayanamsha_ex_at)(tjdut, flags & (FLG_TRUEPOS | FLG_NOABERR | FLG_NOGDEFL))?.1;
    let mut ascmc = trop_ascmc;
    for i in [
        ASCMC_ASC,
        ASCMC_MC,
        ASCMC_VERTEX,
        ASCMC_EQUASC,
        ASCMC_COASC_KOCH,
        ASCMC_COASC_MUNKASEY,
        ASCMC_POLAR_ASC,
    ] {
        ascmc[i] = (ascmc[i] - ayanamsa).rem_euclid(360.0);
    }
    let sid_asc = ascmc[ASCMC_ASC];
    let cusps = match code {
        'W' => {
            // Defensive: the tropical layer refuses non-finite epochs
            // before this arm runs, but a NaN Ascendant would reach the
            // release `math.floor` refusal; fail the same way rather
            // than emitting a NaN wheel.
            if !sid_asc.is_finite() {
                return Err(Error::MathDomain { value: sid_asc });
            }
            let start = (sid_asc / 30.0).floor() * 30.0;
            (0..N_CUSPS)
                .map(|i| (start + i as f64 * 30.0).rem_euclid(360.0))
                .collect()
        }
        'A' | 'E' => (0..N_CUSPS)
            .map(|i| (sid_asc + i as f64 * 30.0).rem_euclid(360.0))
            .collect(),
        'V' => {
            let start = (sid_asc - 15.0).rem_euclid(360.0);
            (0..N_CUSPS)
                .map(|i| (start + i as f64 * 30.0).rem_euclid(360.0))
                .collect()
        }
        'N' => (0..N_CUSPS).map(|i| i as f64 * 30.0).collect(),
        _ => trop_cusps
            .into_iter()
            .map(|c| (c - ayanamsa).rem_euclid(360.0))
            .collect(),
    };
    Ok(shape_houses_ex(cusps, ascmc, flags))
}

/// Read an owned sidereal selection for context-threaded calculations.
///
/// Copies the stored base, defining Julian Day, angular offset in degrees and
/// independent projection/timescale flags without recombining their integers.
/// Negative bases and raw epoch/offset IEEE words remain unchanged. A tropical
/// context supplies the default unset selection. This acquires no data or
/// clock and performs no validation or warning delivery: the diagnostic vector
/// is empty because setters already own diagnostics. Tropical entry points
/// ignore the selection; sidereal entry points consume it without shared state.
pub fn sidereal_selection_for_context(
    context: &crate::context::EphemerisContext,
) -> (
    crate::sidereal::SiderealSelection,
    Vec<crate::errors::Warning>,
) {
    let selection = match context.sidereal() {
        crate::context::SiderealConfig::Sidereal(selection) => selection,
        crate::context::SiderealConfig::Tropical => crate::sidereal::SiderealSelection::unset(),
    };
    (selection, Vec::new())
}

/// Extended houses using an owned snapshot of the context's sidereal selection.
///
/// `tjdut` is JD UT1; geographic coordinates are signed degrees with east-positive
/// longitude. Other arguments, output shape/units and errors follow
/// [`houses_ex_via`]. The context is not mutated and no global sidereal state
/// is consulted. Tropical requests ignore the selection; sidereal requests use
/// it with the supplied providers. Setter diagnostics remain with their owner.
#[allow(clippy::too_many_arguments)] // same arity as houses_ex_via plus the context
pub fn houses_ex_with_context(
    context: &crate::context::EphemerisContext,
    tjdut: f64,
    lat_deg: f64,
    lon_deg: f64,
    hsys: char,
    flags: i32,
    jd: &HousesJdProviders<'_>,
    sidereal: &HousesSiderealProviders<'_>,
) -> Result<(Vec<f64>, [f64; N_ASCMC])> {
    let (selection, _) = sidereal_selection_for_context(context);
    houses_ex_via(
        tjdut, lat_deg, lon_deg, hsys, flags, selection, jd, sidereal,
    )
}

#[cfg(test)]
mod armc_fallback_tests {
    use super::*;

    #[test]
    fn polar_primary_uses_the_requested_safe_geometry() {
        let result = houses_armc_with_fallback(123.45, 70.0, 23.44, 'P', 'O', true, 0.0).unwrap();
        let expected = super::super::houses_armc(123.45, 70.0, 23.44, 'O', 0.0).unwrap();
        assert_eq!((result.cusps, result.ascmc), expected);
        assert!(result.used_fallback);
        assert!(result.warning.unwrap().contains("Porphyry"));
    }

    #[test]
    fn valid_primary_does_not_activate_the_fallback() {
        let result = houses_armc_with_fallback(123.45, 41.9, 23.44, 'E', 'O', true, 0.0).unwrap();
        let expected = super::super::houses_armc(123.45, 41.9, 23.44, 'E', 0.0).unwrap();
        assert_eq!((result.cusps, result.ascmc), expected);
        assert!(!result.used_fallback);
        assert!(result.warning.is_none());
    }
}
