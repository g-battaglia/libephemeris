//! Context-bound lunar acquisition and six-component result projection.
//! No alternate model or ephemeris fallback is introduced.
#[cfg(test)]
#[path = "lunar_context_tests.rs"]
mod tests;

use crate::constants::{EARTH, MOON};
use crate::context::ProviderMode;
use crate::errors::{Error as LebError, Result as LebResult};
use crate::lunar::{
    EclipticMoonState, LunarCoverage, mean_lunar_point_state, true_lunar_point_state,
};
use crate::positions::{
    StateProvider, escalate_sealed_provider_miss, frame_matrices_for_date, mat3_vec3,
    rotate_equatorial_to_ecliptic,
};
/// Shaped mean-point state over injected nutation providers.
///
/// Pure wiring around the public
/// [`crate::lunar::mean_lunar_point_state`] dispatch: no
/// provider contact except through `dpsi_rad_at`/`nut_rate_at`, so
/// unit tests pin the flag shaping with stub providers and no bound
/// image. `flags` are the computation flags with the sidereal and
/// equatorial bits already masked out (both tails are caller-side);
/// `FLG_SPEED` gates the speeds, like the model branch.
pub fn mean_lunar_state_with_providers(
    body: i32,
    jd_tt: f64,
    flags: i32,
    dpsi_rad_at: &dyn Fn(f64) -> LebResult<f64>,
    nut_rate_at: &dyn Fn(f64, f64) -> LebResult<f64>,
) -> LebResult<[f64; 6]> {
    let (lon, lat, dist, dlon, dlat, ddist) = mean_lunar_point_state(
        body,
        jd_tt,
        flags,
        LunarCoverage::Unenforced,
        dpsi_rad_at,
        nut_rate_at,
    )?;
    Ok([lon, lat, dist, dlon, dlat, ddist])
}

/// Shaped true-point state over an injected Moon-state provider.
///
/// Pure wiring around the public
/// [`crate::lunar::true_lunar_point_state`] dispatch: the only
/// provider contact is `true_state_at` (TT epoch in, Cartesian
/// true-ecliptic-of-date Moon state out) plus the nutation pair, so
/// unit tests pin the flag shaping with stub states and no bound
/// image. `flags` are the computation flags with the sidereal and
/// equatorial bits already masked out (both tails are caller-side);
/// `FLG_SPEED` gates the boundary-fallback stencil speeds (per-point
/// half steps, finer under `FLG_SPEED3`/`FLG_TOPOCTR`), like the
/// the model branch. Anything but the two true ids is the dispatch's
/// `UnknownBody` refusal.
pub fn true_lunar_state_with_providers(
    body: i32,
    jd_tt: f64,
    flags: i32,
    true_state_at: &dyn Fn(f64) -> LebResult<EclipticMoonState>,
    dpsi_rad_at: &dyn Fn(f64) -> LebResult<f64>,
    nut_rate_at: &dyn Fn(f64, f64) -> LebResult<f64>,
) -> LebResult<[f64; 6]> {
    let (lon, lat, dist, dlon, dlat, ddist) =
        true_lunar_point_state(body, jd_tt, flags, true_state_at, dpsi_rad_at, nut_rate_at)?;
    Ok([lon, lat, dist, dlon, dlat, ddist])
}

/// Cartesian true-ecliptic-of-date Moon state over injected vectors.
///
/// Pure composition mirroring the producer's own sealed acquisition
/// : difference the geocentric ICRS vectors,
/// rotate by the biased Vondrak precession-nutation matrix, then to
/// the ecliptic at the true obliquity of date. Position and velocity
/// take the same rotation, exactly like the producer's sealed leg.
/// No kernel is reimplemented: every matrix comes from the public
/// [`frame_matrices_for_date`] assembly over the caller's nutation
/// pair. Only a non-finite epoch refuses here; vector finiteness is
/// the frame assembly's contract, and state finiteness gates at the
/// dispatch (`check_ecliptic_state`).
pub fn ecliptic_moon_state_of(
    moon_pos: [f64; 3],
    moon_vel: [f64; 3],
    earth_pos: [f64; 3],
    earth_vel: [f64; 3],
    jd_tt: f64,
    dpsi_rad: f64,
    deps_rad: f64,
) -> LebResult<EclipticMoonState> {
    if !jd_tt.is_finite() {
        return Err(LebError::NonFiniteArgument { value: jd_tt });
    }
    let frames = frame_matrices_for_date(jd_tt, dpsi_rad, deps_rad)?;
    let to_ecl = |v: [f64; 3]| -> LebResult<[f64; 3]> {
        let eq = mat3_vec3(frames.pn(), v)?;
        rotate_equatorial_to_ecliptic(eq, frames.eps_true_rad())
    };
    let r = to_ecl([
        moon_pos[0] - earth_pos[0],
        moon_pos[1] - earth_pos[1],
        moon_pos[2] - earth_pos[2],
    ])?;
    let v = to_ecl([
        moon_vel[0] - earth_vel[0],
        moon_vel[1] - earth_vel[1],
        moon_vel[2] - earth_vel[2],
    ])?;
    Ok(((r[0], r[1], r[2]), (v[0], v[1], v[2])))
}

/// Sealed Moon-state acquisition over the bound provider.
///
/// Acquire geometric LEB Moon/Earth ICRS states at TT `jd_tt` in
/// AU and AU/day, with no alternate source:
/// (support-read misses escalate through the public sealed-miss
/// mapping, so a support miss keeps its identity instead of naming
/// the requested point), differenced and rotated by
/// [`ecliptic_moon_state_of`] over the bound stored nutation pair.
/// `nutation_at` is the serving stack's stored-nutation feed (never an
/// invented value). Shared by the sealed position arms and the
/// Moon-osculating nodes/apsides composition.
pub fn sealed_true_state_at<P: StateProvider>(
    provider: &P,
    nutation_at: &dyn Fn(f64) -> LebResult<(f64, f64)>,
    jd_tt: f64,
) -> LebResult<EclipticMoonState> {
    if !jd_tt.is_finite() {
        return Err(LebError::NonFiniteArgument { value: jd_tt });
    }
    let moon = provider.geometric_state(MOON, jd_tt).map_err(|e| {
        escalate_sealed_provider_miss(ProviderMode::SealedLeb, provider, MOON, jd_tt, e)
    })?;
    let earth = provider.geometric_state(EARTH, jd_tt).map_err(|e| {
        escalate_sealed_provider_miss(ProviderMode::SealedLeb, provider, EARTH, jd_tt, e)
    })?;
    let (dpsi_rad, deps_rad) = nutation_at(jd_tt)?;
    ecliptic_moon_state_of(
        moon.pos(),
        moon.vel(),
        earth.pos(),
        earth.vel(),
        jd_tt,
        dpsi_rad,
        deps_rad,
    )
}
