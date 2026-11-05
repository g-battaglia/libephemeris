//! Catalog constructors, linear space motion and geometric ICRS projection.
//!
//! This leaf owns the angle conventions and propagation arithmetic, not clock
//! acquisition or apparent-place corrections. Native keyword angles and catalog
//! degree angles intentionally keep their separate conversion expressions.

use super::{
    AcquiredStarRow, Error, NativeStarInputs, STAR_ASEC2RAD, STAR_AU_KM, STAR_C_M_PER_S,
    STAR_PARALLAX_FLOOR_MAS, STAR_POLE_COS_DEC_FLOOR, STAR_SECONDS_PER_DAY, StarData,
    star_entry_by_id,
};

/// Native sine and cosine for an angle in radians.
#[inline]
pub(super) fn sin_cos(x: f64) -> (f64, f64) {
    x.sin_cos()
}

/// Right-ascension increment in degrees from an acquired proper-motion prefix.
///
/// `dt_years` is the Julian-year interval, `pm_ra_deg_per_year` is the
/// cosine-weighted angular rate, and `cos_dec` is the reference declination's
/// cosine. The division precedes multiplication. At an undefined pole
/// (`abs(cos_dec) <= 1e-10`) or a NaN cosine, return `None`, leaving the
/// caller's original right ascension untouched. Other non-finite inputs
/// propagate through the selected IEEE arithmetic; no finite check is added.
pub fn proper_motion_ra_increment(
    dt_years: f64,
    pm_ra_deg_per_year: f64,
    cos_dec: f64,
) -> Option<f64> {
    (cos_dec.abs() > STAR_POLE_COS_DEC_FLOOR).then(|| (pm_ra_deg_per_year / cos_dec) * dt_years)
}

/// Propagate a star position by linear proper motion between two epochs
///
/// `pm_ra_cosdec_arcsec_yr` already includes the `cos(dec)` factor
/// (Hipparcos `mu_alpha*` convention); the RA change divides it back out
/// at the reference declination. Time runs in Julian years of exactly
/// `DAYS_PER_JULIAN_YEAR` days. Declination clamps to `[-90, 90]`; right
/// ascension normalizes with `rem_euclid(360.0)`, followed by a
/// negative-angle guard. At the poles (`|cos(dec)| <= 1e-10`) RA is undefined and
/// kept unchanged.
///
/// Linear first-order model in RA/Dec; over long spans high-proper-motion
/// stars need the rigorous 3D space-motion path used by
/// [`star_icrs_state`].
pub fn propagate_proper_motion(
    ra_epoch_deg: f64,
    dec_epoch_deg: f64,
    pm_ra_cosdec_arcsec_yr: f64,
    pm_dec_arcsec_yr: f64,
    from_jd: f64,
    to_jd: f64,
) -> (f64, f64) {
    use crate::constants::DAYS_PER_JULIAN_YEAR;
    let dt_years = (to_jd - from_jd) / DAYS_PER_JULIAN_YEAR;
    let pm_ra_deg_per_year = pm_ra_cosdec_arcsec_yr / 3600.0;
    let pm_dec_deg_per_year = pm_dec_arcsec_yr / 3600.0;
    let dec_target = (dec_epoch_deg + pm_dec_deg_per_year * dt_years).clamp(-90.0, 90.0);
    let cos_dec = dec_epoch_deg.to_radians().cos();
    let mut ra_target = proper_motion_ra_increment(dt_years, pm_ra_deg_per_year, cos_dec)
        .map_or(ra_epoch_deg, |increment| ra_epoch_deg + increment);
    ra_target = ra_target.rem_euclid(360.0);
    if ra_target < 0.0 {
        ra_target += 360.0;
    }
    (ra_target, dec_target)
}

/// Barycentric ICRS state of a catalog star at `jd_tt`.
///
/// Returns `(position_au, velocity_au_per_day)`: the J2000 position from
/// `(ra, dec, parallax)` and the orthogonal space-velocity components
/// from `(pm_ra, pm_dec, radial velocity)` with the Doppler `k` factor,
/// propagated linearly by `jd_tt - J2000`. The `pm_ra` catalog value
/// already includes `cos(dec)`, matching the calculation order
/// `(-pmr * sra - pmd * sdc * cra + rvl * cdc * cra, ...)`.
///
/// Velocity scales by the Doppler factor `k = 1/(1 - RV/c)`: proper-motion
/// terms via parallax and days-per-Julian-year, radial term via 86400 s/day
/// and the AU in km; non-positive parallaxes resolve at
/// [`STAR_PARALLAX_FLOOR_MAS`].
pub fn star_icrs_state(data: &StarData, jd_tt: f64) -> ([f64; 3], [f64; 3]) {
    use crate::constants::{DAYS_PER_JULIAN_YEAR, J2000};
    let ra_rad = data.ra_j2000_deg().to_radians();
    let dec_rad = data.dec_j2000_deg().to_radians();
    star_icrs_state_from_radians(data, ra_rad, dec_rad, jd_tt, J2000, DAYS_PER_JULIAN_YEAR)
}

/// Raw J2000 ICRS AU/AU-day vectors from acquired native keyword units.
/// Native angle values are hours/24*tau and degrees/360*tau. Proper motion
/// remains mas/year; parallax and Doppler conventions follow the ordinary
/// constructor. This arithmetic retains derived IEEE values and has no new
/// finite-input guard. The adapter owns native constructor input admission.
pub fn star_initial_state_from_native(inputs: NativeStarInputs) -> ([f64; 3], [f64; 3]) {
    use crate::positions::reference_frames::FRAME_TAU;
    let px = if inputs.parallax_mas > 0.0 {
        inputs.parallax_mas
    } else {
        STAR_PARALLAX_FLOOR_MAS
    };
    let ra = inputs.ra_hours / 24.0 * FRAME_TAU;
    let dec = inputs.dec_degrees / 360.0 * FRAME_TAU;
    let cra = ra.cos();
    let sra = ra.sin();
    let cdc = dec.cos();
    let sdc = dec.sin();
    let distance = 1.0 / (px * 1.0e-3 * STAR_ASEC2RAD).sin();
    let pos = [distance * cdc * cra, distance * cdc * sra, distance * sdc];
    let k = 1.0 / (1.0 - inputs.radial_km_per_s / STAR_C_M_PER_S * 1000.0);
    let pmr = inputs.ra_mas_per_year / (px * crate::constants::DAYS_PER_JULIAN_YEAR) * k;
    let pmd = inputs.dec_mas_per_year / (px * crate::constants::DAYS_PER_JULIAN_YEAR) * k;
    let rvl = (inputs.radial_km_per_s * STAR_SECONDS_PER_DAY) / STAR_AU_KM * k;
    (pos, space_velocity(sra, cra, sdc, cdc, pmr, pmd, rvl))
}

/// Acquire one native catalog tuple and compute its J2000 ICRS state.
///
/// Results are position AU and velocity AU/day. Catalog acquisition errors
/// propagate before the common native space-motion calculation; no arithmetic
/// callbacks or mutable per-operand reads are involved.
pub fn star_initial_state_from_stored(
    row: &dyn AcquiredStarRow,
) -> Result<([f64; 3], [f64; 3]), Error> {
    Ok(star_initial_state_from_native(row.native_inputs()?))
}

/// Resolve orthogonal RA, declination and radial amplitudes in the ICRS basis.
/// Inputs/outputs are AU/day; association preserves the scalar constructor values.
fn space_velocity(
    sra: f64,
    cra: f64,
    sdc: f64,
    cdc: f64,
    pmr: f64,
    pmd: f64,
    rvl: f64,
) -> [f64; 3] {
    [
        -pmr * sra - pmd * sdc * cra + rvl * cdc * cra,
        pmr * cra - pmd * sdc * sra + rvl * cdc * sra,
        pmd * cdc + rvl * sdc,
    ]
}

/// ICRS catalog space motion with right ascension in hours and declination
/// in degrees converted to turns, then radians.
///
/// The finite TT Julian day and catalog motion supply the same propagation
/// model as [`star_icrs_state`]. Results are position in AU and velocity in
/// AU/day. The caller validates catalog data and owns its applicability span;
/// this calculation includes no observer or apparent-place corrections.
pub fn star_icrs_state_catalog_angles(data: &StarData, jd_tt: f64) -> ([f64; 3], [f64; 3]) {
    use crate::constants::{DAYS_PER_JULIAN_YEAR, J2000};
    use crate::positions::reference_frames::FRAME_TAU;
    let ra_rad = data.ra_j2000_deg() / 15.0 / 24.0 * FRAME_TAU;
    let dec_rad = data.dec_j2000_deg() / 360.0 * FRAME_TAU;
    star_icrs_state_from_radians(data, ra_rad, dec_rad, jd_tt, J2000, DAYS_PER_JULIAN_YEAR)
}

/// Raw J2000 ICRS constructor state using stellar hour/degree angles.
///
/// Returns position AU and velocity AU/day before any epoch propagation.
/// RA converts through degrees/15/24 turns and declination through degrees/
/// 360 turns. Parallax is milliarcseconds with the 1 e-6 mas floor for zero
/// or negative values; proper motion is arcsec/year and radial velocity km/s.
/// The Doppler factor and orthogonal components retain the space-motion
/// constructor arithmetic. No velocity-times-zero addition is made here,
/// so the sign of a raw zero coordinate remains available to propagation.
/// Callers validate consumed inputs and propagate exactly once using their
/// appropriate scalar TT or TDB-plus-projected-light-time epoch.
pub(crate) fn star_icrs_initial_state_catalog_angles(data: &StarData) -> ([f64; 3], [f64; 3]) {
    use crate::positions::reference_frames::FRAME_TAU;
    star_icrs_initial_state_from_radians(
        data,
        data.ra_j2000_deg() / 15.0 / 24.0 * FRAME_TAU,
        data.dec_j2000_deg() / 360.0 * FRAME_TAU,
        crate::constants::DAYS_PER_JULIAN_YEAR,
    )
}

/// Constructor arithmetic shared with the catalog propagation path.
///
/// `ra_rad` and `dec_rad` are the selected angle convention in radians;
/// `days_per_year` is its positive day count. Catalog proper motions are
/// arcsec/year, parallax milliarcseconds and radial velocity km/s. Nonpositive
/// parallax uses [`STAR_PARALLAX_FLOOR_MAS`]. For selected parallax `p` in
/// mas, position radius is `1/sin(p*1e-3*STAR_ASEC2RAD)` AU. With
/// `k = 1/(1-RV*1000/STAR_C_M_PER_S)`, the orthogonal motion amplitudes
/// are `pm*1000/(p*days_per_year)*k` AU/day; the radial amplitude is
/// `RV*STAR_SECONDS_PER_DAY/STAR_AU_KM*k` AU/day. The proper-motion RA
/// value already includes cosine of declination. Resolve these amplitudes
/// along the RA, declination and radial basis as on [`star_icrs_state`].
/// Returns J2000 ICRS position AU and velocity AU/day before any epoch offset.
/// This arithmetic helper acquires no clock, frame or provider and returns no
/// typed errors; IEEE singularities/overflow are retained. Callers own input
/// admission and subsequent propagation. It claims no distance accuracy for
/// the conventional zero-parallax floor.
fn star_icrs_initial_state_from_radians(
    data: &StarData,
    ra_rad: f64,
    dec_rad: f64,
    days_per_year: f64,
) -> ([f64; 3], [f64; 3]) {
    let px_mas = if data.parallax_mas() > 0.0 {
        data.parallax_mas()
    } else {
        STAR_PARALLAX_FLOOR_MAS
    };
    let (sra, cra) = sin_cos(ra_rad);
    let (sdc, cdc) = sin_cos(dec_rad);
    let dist_internal = 1.0 / (px_mas * 1.0e-3 * STAR_ASEC2RAD).sin();
    let pos0 = [
        dist_internal * cdc * cra,
        dist_internal * cdc * sra,
        dist_internal * sdc,
    ];
    let k = 1.0 / (1.0 - data.radial_velocity_km_s() / STAR_C_M_PER_S * 1000.0);
    let pmr = (data.pm_ra_arcsec_yr() * 1000.0) / (px_mas * days_per_year) * k;
    let pmd = (data.pm_dec_arcsec_yr() * 1000.0) / (px_mas * days_per_year) * k;
    let rvl = data.radial_velocity_km_s() * STAR_SECONDS_PER_DAY / STAR_AU_KM * k;
    let vel = space_velocity(sra, cra, sdc, cdc, pmr, pmd, rvl);
    (pos0, vel)
}

/// Propagate a constructor state over the caller-selected scalar epoch.
/// [`star_icrs_initial_state_from_radians`] supplies J2000 ICRS AU/AU-day;
/// `jd_tt` and `j2000` are Julian days in the selected catalog clock. Add
/// `velocity*(jd_tt-j2000)` componentwise in the existing order and return
/// position AU plus unchanged velocity AU/day. No new finite-input check,
/// velocity-times-zero shortcut or epoch normalization is introduced.
fn star_icrs_state_from_radians(
    data: &StarData,
    ra_rad: f64,
    dec_rad: f64,
    jd_tt: f64,
    j2000: f64,
    days_per_year: f64,
) -> ([f64; 3], [f64; 3]) {
    let (pos0, vel) = star_icrs_initial_state_from_radians(data, ra_rad, dec_rad, days_per_year);
    let dt_days = jd_tt - j2000;
    let pos = [
        pos0[0] + vel[0] * dt_days,
        pos0[1] + vel[1] * dt_days,
        pos0[2] + vel[2] * dt_days,
    ];
    (pos, vel)
}

/// Geocentric ICRS vector with light time: geo = star - earth_pos_au,
/// dist = |geo|, lt = dist / C_LIGHT_AU_DAY (0.0 at exact-zero
/// coincidence); proper motion is evaluated at jd_tt, not the retarded
/// epoch.
///
/// `earth_pos_au` is the barycentric Earth (or Sun/SSB for the
/// helio/bary centers) ICRS position in AU. Returns `(geo_au, dist_au,
/// light_time_days)`.
pub fn star_geocentric_icrs(
    data: &StarData,
    jd_tt: f64,
    earth_pos_au: [f64; 3],
) -> ([f64; 3], f64, f64) {
    use crate::positions::C_LIGHT_AU_DAY;
    let (pos, _vel) = star_icrs_state(data, jd_tt);
    let geo = [
        pos[0] - earth_pos_au[0],
        pos[1] - earth_pos_au[1],
        pos[2] - earth_pos_au[2],
    ];
    let dist = (geo[0] * geo[0] + geo[1] * geo[1] + geo[2] * geo[2]).sqrt();
    let lt = if dist > 0.0 {
        dist / C_LIGHT_AU_DAY
    } else {
        0.0
    };
    (geo, dist, lt)
}

/// Geometric J2000 ecliptic place `(lon_deg, lat_deg, dist_au)` for a
/// catalog row:
/// [`star_geocentric_icrs`] composed with the frame-projection layer's
/// J2000 rotation and spherical shaping; errors come from the rotation
/// and spherical-shaping calls. Aberration, deflection, of-date frames
/// and the noaberr/nogdefl branch rules stay with the apparent-star pipeline.
pub fn star_geometric_j2000_ecliptic(
    data: &StarData,
    jd_tt: f64,
    earth_pos_au: [f64; 3],
) -> Result<[f64; 3], Error> {
    use crate::positions::{cartesian_to_spherical_deg, rotate_icrs_to_ecliptic_j2000};
    let (geo, _dist, _lt) = star_geocentric_icrs(data, jd_tt, earth_pos_au);
    let ecl = rotate_icrs_to_ecliptic_j2000(geo)?;
    let mut position = cartesian_to_spherical_deg(ecl)?;
    position[2] = crate::positions::vec3_dist(geo);
    Ok(position)
}

/// Identifier-resolving form of [`star_geometric_j2000_ecliptic`]: an
/// unknown numeric id returns the typed [`Error::StarNotFound`] refusal,
/// never a silent zero.
pub fn star_geometric_j2000_by_id(
    star_id: i32,
    jd_tt: f64,
    earth_pos_au: [f64; 3],
) -> Result<[f64; 3], Error> {
    match star_entry_by_id(star_id) {
        Some(entry) => star_geometric_j2000_ecliptic(&entry.data(), jd_tt, earth_pos_au),
        None => Err(Error::StarNotFound {
            star_id: star_id.to_string(),
            search_type: Some("sequential"),
        }),
    }
}
