//! Observer-relative stellar samples, physical corrections and sample rates.
//!
//! Stored scalar-TT and topocentric projected-TDB rays have distinct clock and
//! frame laws. They share the reported-sample contract, not an interchangeable
//! propagation route. Row ownership and request serving stay in the facade.

use super::{
    Error, STAR_GEOCENTRIC_VELOCITY_STEP_DAYS, StarData, StarFrame, StarObserver, StarSample,
    fixed_star_data, star_geocentric_icrs, star_icrs_initial_state_catalog_angles, star_icrs_state,
    star_icrs_state_catalog_angles,
};

/// Apparent ecliptic place `(lon_deg, lat_deg, dist_au)` for a catalog
/// row.
///
/// `observer` is the already-selected center state (Earth, Sun, or the
/// SSB origin); `deflect` serves each deflection-layer deflector position (only
/// consulted when `nogdefl` is false); `frame` supplies the of-date
/// matrix pair (only consulted when `j2000_frame` is false).
/// `jd_tt` is a TT Julian day. The deflector callback receives a body ID and
/// TT Julian day and returns barycentric ICRS position in AU; its errors
/// propagate at the consumed deflection stage. The result is spherical
/// degrees/degrees/AU in the selected J2000 or true-of-date ecliptic frame.
/// Deflection is skipped only when `nogdefl` and aberration only when
/// `noaberr`, independently; errors come from the deflect closure, the
/// missing-frame configuration, and the rotation/spherical calls.
#[allow(clippy::too_many_arguments)]
pub fn star_apparent_ecliptic_via(
    data: &StarData,
    jd_tt: f64,
    observer: &StarObserver,
    noaberr: bool,
    nogdefl: bool,
    j2000_frame: bool,
    deflect: impl Fn(i32, f64) -> Result<[f64; 3], Error>,
    frame: Option<StarFrame>,
) -> Result<[f64; 3], Error> {
    star_apparent_sample_via(
        data,
        jd_tt,
        observer,
        noaberr,
        nogdefl,
        j2000_frame,
        deflect,
        frame,
    )
    .map(|sample| sample.position)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn star_apparent_sample_via(
    data: &StarData,
    jd_tt: f64,
    observer: &StarObserver,
    noaberr: bool,
    nogdefl: bool,
    j2000_frame: bool,
    deflect: impl Fn(i32, f64) -> Result<[f64; 3], Error>,
    frame: Option<StarFrame>,
) -> Result<StarSample, Error> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    // Steps 4-5 (): rigorous ICRS state, geocentric composition
    // against the selected observer, one-way light time.
    let (geo, _dist, lt) = star_geocentric_icrs(data, jd_tt, observer.pos_au());
    let (_, velocity) = star_icrs_state(data, jd_tt);
    star_apparent_ecliptic_ray_via(
        geo,
        std::array::from_fn(|i| velocity[i] - observer.vel_au_day()[i]),
        lt,
        jd_tt,
        observer,
        noaberr,
        nogdefl,
        j2000_frame,
        deflect,
        || Ok(frame),
    )
}

/// Reduce an acquired stored constructor through the ordinary apparent tail.
/// Initial position/velocity are ICRS AU/AU-day; scalar TT propagation and
/// observer subtraction preserve the existing association. `frame` is lazy:
/// it is acquired only after deflection/aberration and only for the of-date arm.
/// No row, observer or alternate provider is selected by this helper.
#[allow(clippy::too_many_arguments)]
pub(crate) fn star_apparent_ecliptic_initial_via(
    initial: ([f64; 3], [f64; 3]),
    jd_tt: f64,
    observer: &StarObserver,
    noaberr: bool,
    nogdefl: bool,
    j2000_frame: bool,
    deflect: impl Fn(i32, f64) -> Result<[f64; 3], Error>,
    frame: impl FnOnce() -> Result<Option<StarFrame>, Error>,
) -> Result<StarSample, Error> {
    let (pos, vel) = initial;
    let span = jd_tt - crate::constants::J2000;
    let obs = observer.pos_au();
    let geo = std::array::from_fn(|i| pos[i] + vel[i] * span - obs[i]);
    let distance = crate::positions::vec3_dist(geo);
    let lt = if distance > 0.0 {
        distance / crate::positions::C_LIGHT_AU_DAY
    } else {
        0.0
    };
    star_apparent_ecliptic_ray_via(
        geo,
        std::array::from_fn(|i| vel[i] - observer.vel_au_day()[i]),
        lt,
        jd_tt,
        observer,
        noaberr,
        nogdefl,
        j2000_frame,
        deflect,
        frame,
    )
}

/// Common correction/frame tail over an acquired observer-relative AU ray.
/// Light time is days; the frame callback runs at its consumed of-date stage.
#[allow(clippy::too_many_arguments)]
fn star_apparent_ecliptic_ray_via(
    mut geo: [f64; 3],
    relative_velocity: [f64; 3],
    lt: f64,
    jd_tt: f64,
    observer: &StarObserver,
    noaberr: bool,
    nogdefl: bool,
    j2000_frame: bool,
    deflect: impl Fn(i32, f64) -> Result<[f64; 3], Error>,
    frame: impl FnOnce() -> Result<Option<StarFrame>, Error>,
) -> Result<StarSample, Error> {
    use crate::positions::{
        apply_aberration, apply_gravitational_deflection, cartesian_to_spherical_deg, mat3_vec3,
        rotate_equatorial_to_ecliptic, rotate_icrs_to_ecliptic_j2000,
    };
    let (distance, radial_rate) = star_geometric_range_state(geo, relative_velocity);
    // Step 6: PPN deflection, skipped only when `nogdefl`.
    if !nogdefl {
        geo = apply_gravitational_deflection(geo, observer.pos_au(), jd_tt, lt, deflect)?;
    }
    // Step 7: aberration, skipped only when `noaberr`.
    if !noaberr {
        geo = apply_aberration(geo, observer.vel_au_day(), lt)?;
    }
    // Step 8: frame rotation.
    let mut position = if j2000_frame {
        let ecl = rotate_icrs_to_ecliptic_j2000(geo)?;
        cartesian_to_spherical_deg(ecl)
    } else {
        let frame = frame()?.ok_or(Error::ConfigurationMissing {
            missing_config: "star_frame",
            suggestion: "of-date star places need the pn/eps_true frame pair",
        })?;
        let geo_eq = mat3_vec3(frame.pn(), geo)?;
        let ecl = rotate_equatorial_to_ecliptic(geo_eq, frame.eps_true_rad())?;
        cartesian_to_spherical_deg(ecl)
    }?;
    // Range is invariant under the orthogonal coordinate-frame rotation.
    position[2] = distance;
    Ok(StarSample {
        position,
        radial_rate,
    })
}

/// Propagation epoch for an observer-dependent stellar ray.
///
/// The epoch is TDB plus the projection of BCRS observer position onto the
/// J2000 catalog direction, divided by light speed. The geocentric catalog
/// path uses TT directly; this topocentric construction retains the propagation
/// epoch correction. Nonfinite epochs/observer states and a degenerate catalog
/// position return NonFiniteArgument.
pub fn topo_star_epoch_days(
    data: &StarData,
    jd_tt: f64,
    observer_pos_au: [f64; 3],
) -> Result<f64, Error> {
    use crate::constants::J2000;
    use crate::positions::FRAME_LIGHT_SPEED_AU_DAY;
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    for v in observer_pos_au {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    let t_tdb = crate::time::tt_to_tdb_jd(jd_tt)?;
    let (pos0, _vel) = star_icrs_state_catalog_angles(data, J2000);
    let n0 = (pos0[0] * pos0[0] + pos0[1] * pos0[1] + pos0[2] * pos0[2]).sqrt();
    if !n0.is_finite() || n0 == 0.0 {
        return Err(Error::NonFiniteArgument { value: n0 });
    }
    let dt_extra = (pos0[0] / n0 * observer_pos_au[0]
        + pos0[1] / n0 * observer_pos_au[1]
        + pos0[2] / n0 * observer_pos_au[2])
        / FRAME_LIGHT_SPEED_AU_DAY;
    Ok(t_tdb + dt_extra)
}

/// Apparent ecliptic place for a topocentric observer, with projected-TDB
/// propagation and equinox-based frames rather than the scalar-TT propagation
/// and stored frame pair of [`star_apparent_ecliptic_via`].
///
/// Sequence: the catalog star propagates to [`topo_star_epoch_days`]
/// (TDB plus the observer-projection light-time difference, not
/// `jd_tt`); the observer-relative vector takes its light time over
/// [`crate::positions::FRAME_LIGHT_SPEED_AU_DAY`]; Sun/Jupiter/Saturn PPN deflection unless
/// `nogdefl`; the Earth's own deflection for every non-geometric arm
///
/// evaluated on the Earth-centered site vector carried on the observer
///
/// full special-relativistic aberration over the defined light-speed constant
/// unless `noaberr`; then the equinox-based frame rotation — the fixed
/// J2000 rotation for `j2000_frame`, else the true-ecliptic-of-date
/// rotation, never the
/// LEB-direct `pn`/`eps_true` pair. Deflection is skipped only when
/// `nogdefl` and aberration only when `noaberr`, independently; errors
/// come from the epoch helper, the deflect closure, the bundled
/// timescale/TDB acquisition, and the rotation/spherical calls.
#[allow(clippy::too_many_arguments)]
pub fn star_apparent_ecliptic_topo_via(
    data: &StarData,
    jd_tt: f64,
    observer: &StarObserver,
    noaberr: bool,
    nogdefl: bool,
    j2000_frame: bool,
    deflect: impl Fn(i32, f64) -> Result<[f64; 3], Error>,
) -> Result<[f64; 3], Error> {
    star_apparent_topo_sample_via(
        data,
        jd_tt,
        observer,
        noaberr,
        nogdefl,
        j2000_frame,
        deflect,
    )
    .map(|sample| sample.position)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn star_apparent_topo_sample_via(
    data: &StarData,
    jd_tt: f64,
    observer: &StarObserver,
    noaberr: bool,
    nogdefl: bool,
    j2000_frame: bool,
    deflect: impl Fn(i32, f64) -> Result<[f64; 3], Error>,
) -> Result<StarSample, Error> {
    use crate::positions::reference_frames::{
        FRAME_ASEC2RAD, FRAME_TENTH_USEC_2_RAD, ecliptic_date_matrix, equinox_matrix,
        frame_latlon_deg, iau2000a_tenth_microarcsec, mean_obliquity_arcsec, nutation_matrix,
        precession_matrix,
    };
    use crate::positions::{
        FRAME_LIGHT_SPEED_AU_DAY, apply_aberration_with_c, apply_earth_deflection,
        apply_gravitational_deflection, cartesian_to_spherical_deg, rotate_icrs_to_ecliptic_j2000,
    };
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let obs_pos = observer.pos_au();
    let star_epoch = topo_star_epoch_days(data, jd_tt, obs_pos)?;
    let (star_pos, star_velocity) = star_icrs_state_catalog_angles(data, star_epoch);
    let (initial, _) = star_icrs_initial_state_catalog_angles(data);
    let initial_range = crate::positions::vec3_dist(initial);
    let epoch_rate = crate::time::tdb_rate(jd_tt)?
        + (0..3)
            .map(|i| initial[i] / initial_range * observer.vel_au_day()[i])
            .sum::<f64>()
            / FRAME_LIGHT_SPEED_AU_DAY;
    let relative_velocity =
        std::array::from_fn(|i| star_velocity[i] * epoch_rate - observer.vel_au_day()[i]);
    let mut geo = [
        star_pos[0] - obs_pos[0],
        star_pos[1] - obs_pos[1],
        star_pos[2] - obs_pos[2],
    ];
    let dist = (geo[0] * geo[0] + geo[1] * geo[1] + geo[2] * geo[2]).sqrt();
    let lt = if dist > 0.0 {
        dist / FRAME_LIGHT_SPEED_AU_DAY
    } else {
        0.0
    };
    let (distance, radial_rate) = star_geometric_range_state(geo, relative_velocity);
    // Step 6: PPN deflection, skipped only when `nogdefl`.
    if !nogdefl {
        geo = apply_gravitational_deflection(geo, obs_pos, jd_tt, lt, &deflect)?;
    }
    // Earth term: every non-geometric arm, on the Earth-centered site
    // vector.
    if !(noaberr && nogdefl) {
        let site_gcrs = observer.site_gcrs_au().ok_or(Error::ConfigurationMissing {
            missing_config: "observer_gcrs",
            suggestion: "topocentric star places need the GCRS site vector",
        })?;
        geo = apply_earth_deflection(geo, site_gcrs)?;
    }
    // Step 7: aberration, skipped only when `noaberr`.
    if !noaberr {
        geo = apply_aberration_with_c(geo, observer.vel_au_day(), lt, FRAME_LIGHT_SPEED_AU_DAY)?;
    }
    // Step 8: frame rotation.
    let mut position = if j2000_frame {
        let ecl = rotate_icrs_to_ecliptic_j2000(geo)?;
        cartesian_to_spherical_deg(ecl)
    } else {
        let jd_tdb = crate::time::tt_to_tdb_jd(jd_tt)?;
        let (dpsi_u, deps_u) = iau2000a_tenth_microarcsec(jd_tt)?;
        let dpsi = dpsi_u * FRAME_TENTH_USEC_2_RAD;
        let deps = deps_u * FRAME_TENTH_USEC_2_RAD;
        let mean_obl = mean_obliquity_arcsec(jd_tdb)? * FRAME_ASEC2RAD;
        let nut = nutation_matrix(mean_obl, mean_obl + deps, dpsi)?;
        let prec = precession_matrix(jd_tdb)?;
        let equinox = equinox_matrix(nut, prec)?;
        let ecl_matrix = ecliptic_date_matrix(equinox, mean_obl, deps)?;
        frame_latlon_deg(ecl_matrix, geo)
    }?;
    position[2] = distance;
    Ok(StarSample {
        position,
        radial_rate,
    })
}

/// Astrometric range and its derivative, invariant under direction corrections.
/// The ICRS ray is AU and its derivative is AU per TT day. For positive range
/// R, dR/dTT = (r/R)·dr/dTT. Computing this dot product avoids subtraction of
/// separately rounded million-AU positions or norms. At coincident points the
/// direction is undefined and the existing zero-range convention returns zero.
pub(crate) fn star_geometric_range_state(ray: [f64; 3], velocity: [f64; 3]) -> (f64, f64) {
    let range = crate::positions::vec3_dist(ray);
    let rate = if range == 0.0 {
        0.0
    } else {
        crate::positions::compensated_dot3(ray.map(|value| value / range), velocity)
    };
    (range, rate)
}

/// Angular finite differences with the central sample's analytic range rate.
/// Epochs and half-step are TT Julian days; no radial difference is formed.
pub(crate) fn star_velocity_samples_via(
    sample_at: impl Fn(f64) -> Result<StarSample, Error>,
    jd_tt: f64,
    h_days: f64,
) -> Result<[f64; 6], Error> {
    let central_rate = std::cell::Cell::new(None);
    let mut result = star_velocity_via(
        |epoch| {
            let sample = sample_at(epoch)?;
            if central_rate.get().is_none() {
                central_rate.set(Some(sample.radial_rate));
            }
            Ok(sample.position)
        },
        jd_tt,
        h_days,
    )?;
    result[5] = central_rate.get().expect("central stellar sample");
    Ok(result)
}

/// Central-difference velocity for the star pipeline: positions at
/// `jd_tt ± h_days` through `position_at`, longitude wrapped to the
/// shortest arc before dividing by `2h`. Returns `(lon, lat, dist,
/// speed_lon, speed_lat, speed_dist)` with only longitude wrapped; the
/// distance speed is the centered derivative of the observer-relative
/// distance.
///
/// `h` rule (`0.005` for topocentric requests, `0.5`
/// otherwise) lives at the entry layer; this core takes `h_days` as
/// given without selecting a different stencil.
pub fn star_velocity_via(
    position_at: impl Fn(f64) -> Result<[f64; 3], Error>,
    jd_tt: f64,
    h_days: f64,
) -> Result<[f64; 6], Error> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    if !(h_days.is_finite() && h_days > 0.0) {
        return Err(Error::NonFiniteArgument { value: h_days });
    }
    let center = position_at(jd_tt)?;
    let prev = position_at(jd_tt - h_days)?;
    let next = position_at(jd_tt + h_days)?;
    let mut speed_lon = next[0] - prev[0];
    if speed_lon > 180.0 {
        speed_lon -= 360.0;
    } else if speed_lon < -180.0 {
        speed_lon += 360.0;
    }
    Ok([
        center[0],
        center[1],
        center[2],
        speed_lon / (2.0 * h_days),
        (next[1] - prev[1]) / (2.0 * h_days),
        (next[2] - prev[2]) / (2.0 * h_days),
    ])
}

/// Sealed apparent star place over a bound LEB reader.
///
/// Observer states come from the serving provider (`EARTH` normally;
/// `SUN` for `FLG_HELCTR`; the SSB origin for `FLG_BARYCTR`), the
/// of-date frame from the serving nutation plus the native Vondrák
/// bundle, deflectors from the same provider. `center_flag` is `0`,
/// `FLG_HELCTR` or `FLG_BARYCTR` (see [`super::fixstar_observation_center`]);
/// `noaberr`/`nogdefl` travel as given (entry-level helio/bary
/// forcing is honored by the heliocentric/barycentric entry handling).
/// Errors: unknown-id refusal,
/// observer/deflector/frame provider passthrough, and the downstream
/// correction and shaping calls.
#[allow(clippy::too_many_arguments)]
pub fn sealed_star_position<P: crate::positions::StateProvider>(
    provider: &P,
    nutation_at: &dyn Fn(f64) -> Result<(f64, f64), Error>,
    star_id: i32,
    jd_tt: f64,
    noaberr: bool,
    nogdefl: bool,
    j2000_frame: bool,
    center_flag: i32,
) -> Result<[f64; 3], Error> {
    sealed_star_sample(
        provider,
        nutation_at,
        star_id,
        jd_tt,
        noaberr,
        nogdefl,
        j2000_frame,
        center_flag,
    )
    .map(|sample| sample.position)
}

#[allow(clippy::too_many_arguments)]
fn sealed_star_sample<P: crate::positions::StateProvider>(
    provider: &P,
    nutation_at: &dyn Fn(f64) -> Result<(f64, f64), Error>,
    star_id: i32,
    jd_tt: f64,
    noaberr: bool,
    nogdefl: bool,
    j2000_frame: bool,
    center_flag: i32,
) -> Result<StarSample, Error> {
    use crate::constants::{EARTH, FLG_BARYCTR, FLG_HELCTR, SUN};
    let data = fixed_star_data(star_id).ok_or_else(|| Error::StarNotFound {
        star_id: star_id.to_string(),
        search_type: Some("sequential"),
    })?;
    let observer = if center_flag == FLG_HELCTR {
        let state = provider.geometric_state(SUN, jd_tt)?;
        StarObserver::new(state.pos(), state.vel())?
    } else if center_flag == FLG_BARYCTR {
        StarObserver::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0])?
    } else {
        let state = provider.geometric_state(EARTH, jd_tt)?;
        StarObserver::new(state.pos(), state.vel())?
    };
    let frame = if j2000_frame {
        None
    } else {
        let (dpsi, deps) = nutation_at(jd_tt)?;
        let frames = crate::positions::frame_matrices_for_date(jd_tt, dpsi, deps)?;
        Some(StarFrame::new(frames.pn(), frames.eps_true_rad())?)
    };
    star_apparent_sample_via(
        &data,
        jd_tt,
        &observer,
        noaberr,
        nogdefl,
        j2000_frame,
        |id, jd| Ok(provider.geometric_state(id, jd)?.pos()),
        frame,
    )
}

/// Sealed six-slot star place with the geocentric `h = 0.5`
/// central-difference velocity (sample epochs re-evaluate the observer,
/// deflectors and frame through [`sealed_star_position`]).
#[allow(clippy::too_many_arguments)]
pub fn sealed_star_velocity<P: crate::positions::StateProvider>(
    provider: &P,
    nutation_at: &dyn Fn(f64) -> Result<(f64, f64), Error>,
    star_id: i32,
    jd_tt: f64,
    noaberr: bool,
    nogdefl: bool,
    j2000_frame: bool,
    center_flag: i32,
) -> Result<[f64; 6], Error> {
    star_velocity_samples_via(
        |jd| {
            sealed_star_sample(
                provider,
                nutation_at,
                star_id,
                jd,
                noaberr,
                nogdefl,
                j2000_frame,
                center_flag,
            )
        },
        jd_tt,
        STAR_GEOCENTRIC_VELOCITY_STEP_DAYS,
    )
}
