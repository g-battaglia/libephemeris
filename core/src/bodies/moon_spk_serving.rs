//! Local satellite-SPK apparent places with independent model and frame clocks.

use super::{
    MoonObserver, MoonPositionServing, MoonSpkRegistry, moon_observer_for_flags,
    resolve_moon_target,
};
use crate::context::EphemerisContext;
use crate::errors::{Error, Result};
use crate::sidereal::{SealedAyanamshaServing, SiderealSelection};
use crate::spk::SpkRegistry;

/// Satellite apparent place at scalar TT using the context's registered DE states.
///
/// `jd_tt` is a finite TT Julian Date, not UT1. `moon_id` is a legacy or
/// canonical satellite id. The separate satellite registry supplies the SSB
/// target chain; the context registry supplies Earth, Sun and live sidereal
/// anchors. Earth/Sun registrations are selected by NAIF identity, regardless
/// of their library body keys. Center chains resolve from registered NAIF
/// identities or segments in the selected image, never unrelated images.
/// Explicit sealed LEB mode refuses JPL source access, even for a
/// barycentric observation. Other modes use these local states and never
/// substitute a bound LEB image or perform network acquisition.
///
/// Returns longitude/latitude in degrees, distance in AU and their daily
/// speeds. Unknown or unregistered targets return `None` before source access.
/// Observer precedence is barycentric, heliocentric, topocentric, geocentric.
/// A topocentric request requires a configured site. Flags control light time,
/// aberration, longitude nutation stripping, sidereal reduction and SPEED.
/// This internal ecliptic-of-date entry does not perform public equatorial or
/// J2000 output conversion. SPEED uses a fourth-order dyadic TT stencil for
/// angles and an analytic derivative of the light-time chain for range.
///
/// `delta_cfg` must be the caller-resolved context Delta-T configuration. It
/// affects model-defined sidereal epochs only, not the TT input, SPK TDB
/// acquisition, frame epoch or bundled TT/UT1 clock. Non-finite epochs,
/// missing source/site configuration and SPK coverage/provider errors propagate.
// The native site model uses the IERS 2010 ellipsoid and equinox-based Earth rotation.
pub fn moon_position_native_tt(
    ctx: &EphemerisContext,
    moons: &MoonSpkRegistry,
    moon_id: i32,
    jd_tt: f64,
    flags: i32,
    delta_cfg: &crate::time::DeltaTConfig<'_>,
) -> Result<Option<[f64; 6]>> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    if resolve_moon_target(moons, moon_id).is_none() {
        return Ok(None);
    }
    crate::positions::forbid_jpl_in_sealed_mode(ctx.mode())?;
    if ctx.spk_registry().is_empty() {
        return Err(Error::ConfigurationMissing {
            missing_config: "jpl_ephemeris",
            suggestion: "register the selected local DE kernel in the context before calculating satellite positions",
        });
    }
    let topo = ctx
        .observer()
        .map(|site| (site.lon_deg(), site.lat_deg(), site.height_m()));
    let (selection, _) = crate::houses::sidereal_selection_for_context(ctx);
    let deltat = |ut| crate::time::deltat_with_config(ut, delta_cfg);
    moon_position_spk_tt(
        ctx.spk_registry(),
        moons,
        moon_id,
        jd_tt,
        flags,
        topo,
        selection,
        &deltat,
    )
}

/// Borrowed local-SPK satellite kernel shared by production and comparisons.
///
/// All state inputs are AU in ICRF; SPK velocities per TDB day are scaled
/// by d(TDB)/d(TT) before adding a site velocity per TT day. Scalar TT is split before its TDB
/// correction; satellite and observer evaluations retain those parts. Model
/// Delta-T is consumed only by sidereal defining epochs. The true frame,
/// nutation strip and site rotation retain their independent acquisition routes.
/// Returns six ecliptic-of-date values or `None`; errors are never substituted.
#[allow(clippy::too_many_arguments)]
pub(crate) fn moon_position_spk_tt(
    de: &SpkRegistry,
    moons: &MoonSpkRegistry,
    moon_id: i32,
    jd_tt: f64,
    flags: i32,
    topo: Option<(f64, f64, f64)>,
    selection: SiderealSelection,
    deltat_days_at: &dyn Fn(f64) -> Result<f64>,
) -> Result<Option<[f64; 6]>> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let Some(naif_id) = resolve_moon_target(moons, moon_id) else {
        return Ok(None);
    };
    let satellite = |_: i32, tt| {
        let (day, fraction) = tt_to_tdb_parts(tt)?;
        moons
            .eval_moon_ssb_au_two_part_jd(moon_id, day, fraction)
            .map_err(|error| {
                if matches!(error, Error::SpkOutOfRange { .. }) {
                    let (start_jd, end_jd) = moons
                        .get_moon_coverage(moon_id)
                        .unwrap_or((f64::NAN, f64::NAN));
                    Error::EphemerisRange {
                        requested_jd: tt,
                        start_jd,
                        end_jd,
                        body: Some(moon_id),
                    }
                } else {
                    error
                }
            })
            .and_then(|(position, velocity)| {
                let rate = crate::time::tdb_rate(tt)?;
                Ok((position, velocity.map(|value| value * rate)))
            })
    };
    position_spk_tt(
        de,
        naif_id,
        jd_tt,
        flags,
        topo,
        selection,
        deltat_days_at,
        &satellite,
        None,
    )
}

type TargetStateAt<'a> = dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])> + 'a;
type EclipticAt<'a> = dyn Fn([f64; 3], f64) -> Result<(f64, f64, f64)> + 'a;

fn tt_to_tdb_parts(tt: f64) -> Result<(f64, f64)> {
    let day = tt.floor();
    let fraction = tt - day;
    let correction =
        crate::positions::reference_frames::tdb_minus_tt_seconds(day, fraction)? / 86400.0;
    Ok((day, fraction + correction))
}

#[allow(clippy::too_many_arguments)]
fn position_spk_tt(
    de: &SpkRegistry,
    naif_id: i32,
    jd_tt: f64,
    flags: i32,
    topo: Option<(f64, f64, f64)>,
    selection: SiderealSelection,
    deltat_days_at: &dyn Fn(f64) -> Result<f64>,
    target_at_tt: &TargetStateAt<'_>,
    frame_override: Option<&EclipticAt<'_>>,
) -> Result<Option<[f64; 6]>> {
    use crate::positions::reference_frames::frame_latlon_deg;
    let anchor_states =
        crate::positions::LocalSpkStateProvider::new(de, crate::context::ProviderMode::NativeJpl)?;
    let anchors = crate::stars::SealedAnchorProvider::new(&anchor_states);
    let ayanamsha = SealedAyanamshaServing::new(selection, deltat_days_at, Some(&anchors));
    let observer = moon_observer_for_flags(flags);
    let site = if observer == MoonObserver::Topocentric {
        topo.map(|(lon, lat, height)| crate::stars::StarTopoSite::new(lon, lat, height))
            .transpose()?
    } else {
        None
    };
    let observer_at = |tt| {
        let target = match observer {
            MoonObserver::Barycenter => return Ok(([0.0; 3], [0.0; 3])),
            MoonObserver::Heliocentric => 10,
            MoonObserver::Geocentric | MoonObserver::Topocentric => 399,
        };
        if observer == MoonObserver::Topocentric && site.is_none() {
            return Err(Error::ConfigurationMissing {
                missing_config: "observer_location",
                suggestion: "FLG_TOPOCTR needs a prior set_topo site",
            });
        }
        let body = de
            .registered_bodies()
            .into_iter()
            .find(|entry| entry.naif_id == target)
            .map(|entry| entry.body)
            .ok_or(Error::SpkUnregisteredBody { body: target })?;
        let (day, fraction) = tt_to_tdb_parts(tt)?;
        let (pos, vel) = de.eval_geometric_ssb_au_two_part_jd(body, day, fraction)?;
        let rate = crate::time::tdb_rate(tt)?;
        let vel = vel.map(|value| value * rate);
        if let Some(site) = site {
            moon_spk_topo_observer(pos, vel, site, tt)
        } else {
            Ok((pos, vel))
        }
    };
    let frame_at = |rel, tt| {
        if let Some(frame) = frame_override {
            return frame(rel, tt);
        }
        let [lon, lat, distance] = frame_latlon_deg(moon_spk_frame(tt)?, rel)?;
        Ok((lon, lat, distance))
    };
    let strip = |tt| crate::erfa::nut06a(tt).map(|(dpsi, _)| dpsi.to_degrees());
    let sidereal = |lon, speed, ut| {
        super::moon_sidereal_word(lon, speed, ut, flags, |epoch| {
            moon_spk_ayanamsha(epoch, flags, &ayanamsha)
        })
    };
    let serving = MoonPositionServing::new(target_at_tt, observer_at, frame_at, strip, sidereal);
    let ut1 = if super::moon_wants_sidereal_correction(flags) {
        crate::time::default_timescale_ut1_from_tt(jd_tt)?
    } else {
        f64::NAN
    };
    let (lon, lat, dist, dlon, dlat, ddist) =
        serving.moon_position(naif_id, observer, jd_tt, ut1, flags)?;
    Ok(Some([lon, lat, dist, dlon, dlat, ddist]))
}

/// Direct explicitly registered SPK position, independent of coefficient mode.
///
/// TT epochs, AU states and TDB/TT velocity conversion share the satellite
/// observation pipeline. Body/frame selectors follow the direct SPK interface:
/// type 2/3 outputs are ecliptic-of-date; type 21 also supports J2000. Cartesian,
/// radians and barycentric output flags belong to the outer public calc API.
/// Unknown registrations return None; actual coverage/observer errors propagate.
#[allow(clippy::too_many_arguments)]
pub(crate) fn registered_spk_position_tt(
    de: &SpkRegistry,
    targets: &SpkRegistry,
    body: i32,
    jd_tt: f64,
    flags: i32,
    topo: Option<(f64, f64, f64)>,
    selection: SiderealSelection,
    deltat_days_at: &dyn Fn(f64) -> Result<f64>,
) -> Result<Option<[f64; 6]>> {
    use crate::constants::*;
    let Some(registration) = targets.body_info(body) else {
        return Ok(None);
    };
    if !crate::positions::registered_spk_position_probe(targets, jd_tt, body, flags)? {
        return Ok(None);
    }
    let type21 = targets.kernel(&registration.kernel).is_some_and(|kernel| {
        kernel
            .index()
            .all_segments()
            .any(|segment| segment.descriptor.seg_type == crate::spk::TYPE_EMDA)
    });
    let mut registry = de.clone();
    registry.merge_shared_from(targets);
    let target = |_: i32, tt| {
        let (whole, fraction) = tt_to_tdb_parts(tt)?;
        // The selected image owns every center link that it contains. Only a
        // center absent from that image is acquired from the main registry.
        // A DE system barycenter must not replace a different ephemeris carried
        // by a self-contained target image.
        let (position, velocity) =
            registry.eval_geometric_ssb_au_two_part_jd(body, whole, fraction)?;
        let rate = crate::time::tdb_rate(tt)?;
        Ok((position, velocity.map(|value| value * rate)))
    };
    let strip = flags & (FLG_NONUT | FLG_J2000) != 0
        || flags & (FLG_SIDEREAL | FLG_EQUATORIAL) == (FLG_SIDEREAL | FLG_EQUATORIAL);
    let frame = |ray: [f64; 3], tt: f64| -> Result<(f64, f64, f64)> {
        if !type21 {
            let [mut lon, lat, radius] =
                crate::positions::reference_frames::frame_latlon_deg(moon_spk_frame(tt)?, ray)?;
            if strip {
                lon = (lon - crate::erfa::nut06a(tt)?.0.to_degrees()).rem_euclid(360.0);
            }
            return Ok((lon, lat, radius));
        }
        let (sin, cos) = crate::spk::OBLIQUITY_J2000_DEG.to_radians().sin_cos();
        let [x, y, z] = [
            ray[0],
            ray[1] * cos + ray[2] * sin,
            -ray[1] * sin + ray[2] * cos,
        ];
        let radius = (x * x + y * y + z * z).sqrt();
        let longitude = y.atan2(x).to_degrees().rem_euclid(360.0);
        let latitude = (z / radius).asin().to_degrees();
        let (mut longitude, latitude) = if flags & FLG_J2000 != 0 {
            (longitude, latitude)
        } else {
            crate::hypothetical::precess_ecliptic_source_angles(longitude, latitude, 2451545.0, tt)?
        };
        if !strip {
            longitude = (longitude + crate::erfa::nut00b(tt)?.0.to_degrees()).rem_euclid(360.0);
        }
        Ok((longitude, latitude, radius))
    };
    let mut ray_flags = flags
        & !(FLG_NONUT
            | FLG_J2000
            | FLG_SIDEREAL
            | FLG_EQUATORIAL
            | FLG_BARYCTR
            | FLG_XYZ
            | FLG_RADIANS);
    if type21 {
        ray_flags &= !FLG_TOPOCTR;
    }
    let mut result = position_spk_tt(
        de,
        registration.naif_id,
        jd_tt,
        ray_flags,
        topo,
        selection,
        deltat_days_at,
        &target,
        Some(&frame),
    )?;
    if super::moon_wants_sidereal_correction(flags) {
        let states = crate::positions::LocalSpkStateProvider::new(
            de,
            crate::context::ProviderMode::NativeJpl,
        )?;
        let anchors = crate::stars::SealedAnchorProvider::new(&states);
        let aya = SealedAyanamshaServing::new(selection, deltat_days_at, Some(&anchors));
        let ut = crate::time::default_timescale_ut1_from_tt(jd_tt)?;
        if let Some(value) = &mut result {
            (value[0], value[3]) =
                super::moon_sidereal_word(value[0], value[3], ut, flags, |epoch| {
                    moon_spk_ayanamsha(epoch, flags, &aya)
                })?;
        }
    }
    Ok(result)
}

/// Add an IERS 2010 geodetic site's GCRS position/velocity to an Earth
/// barycentric state in AU and AU/day at scalar TT.
///
/// Semi-major axis is 6,378,136.6 m and inverse flattening 298.25642.
/// Geodetic latitude/longitude are radians after conversion from degrees;
/// height is metres. Earth rotation uses the bundled TT-minus-UT1 curve
/// and IAU equinox/GAST primitives, with no polar-motion table.
/// Non-finite coordinates and frame/provider errors propagate.
// The native site model uses the IERS 2010 ellipsoid and equinox-based Earth rotation.
fn moon_spk_topo_observer(
    earth_pos: [f64; 3],
    earth_vel: [f64; 3],
    site: crate::stars::StarTopoSite,
    tt: f64,
) -> Result<([f64; 3], [f64; 3])> {
    use crate::positions::reference_frames::{
        FRAME_ASEC2RAD, FRAME_AU_M, FRAME_TAU, FRAME_TENTH_USEC_2_RAD, equinox_matrix, gast_hours,
        iau2000a_tenth_microarcsec, itrs_rotation_matrix, itrs_to_gcrs, mean_obliquity_arcsec,
        nutation_matrix, precession_matrix, scalar_cos, scalar_sin,
    };
    let lat = site.lat_deg() / 360.0 * FRAME_TAU;
    let lon = site.lon_deg() / 360.0 * FRAME_TAU;
    let radius = 6378136.6 / FRAME_AU_M;
    let height = site.height_m() / FRAME_AU_M;
    let sinphi = scalar_sin(lat);
    let cosphi = scalar_cos(lat);
    let omf = (298.25642 - 1.0) / 298.25642;
    let omf2 = omf * omf;
    let c = 1.0 / (cosphi * cosphi + sinphi * sinphi * omf2).sqrt();
    let s = omf2 * c;
    let xy = (radius * c + height) * cosphi;
    let itrs = [
        xy * scalar_cos(lon),
        xy * scalar_sin(lon),
        (radius * s + height) * sinphi,
    ];
    let tdb = crate::time::tt_to_tdb_jd(tt)?;
    let (psi, deps) = iau2000a_tenth_microarcsec(tt)?;
    let psi = psi * FRAME_TENTH_USEC_2_RAD;
    let deps = deps * FRAME_TENTH_USEC_2_RAD;
    let eps = mean_obliquity_arcsec(tdb)? * FRAME_ASEC2RAD;
    let nut = nutation_matrix(eps, eps + deps, psi)?;
    let equinox = equinox_matrix(nut, precession_matrix(tdb)?)?;
    let dt_seconds = crate::smh2016::bundled_deltat_seconds(tt)?;
    let gast = gast_hours(tt, dt_seconds, tdb)?;
    let rotation = itrs_rotation_matrix(equinox, gast)?;
    let (off_pos, off_vel) = itrs_to_gcrs(rotation, itrs)?;
    Ok((
        std::array::from_fn(|i| earth_pos[i] + off_pos[i]),
        std::array::from_fn(|i| earth_vel[i] + off_vel[i]),
    ))
}

/// Flag-aware scalar ayanamsha in degrees at UT1 for SPK moon composition.
///
/// Ordinary evaluation TT comes from the bundled timescale. Caller-model
/// TT remains distinct for user-UT and Valens defining epochs. The
/// NONUT/J2000 flags select mean ayanamsha; otherwise native IAU 2006/2000 A
/// nutation is added without renormalizing the sum. The bound serving
/// owns selection, model-defined epochs and the DE-backed live anchors;
/// ordinary evaluation TT never uses its caller-model Delta-T callback.
pub(crate) fn moon_spk_ayanamsha(
    ut: f64,
    flags: i32,
    aya: &SealedAyanamshaServing<'_>,
) -> Result<f64> {
    let tt = crate::time::default_timescale_tt_from_ut1(ut)?;
    let (noaberr, nogdefl, use_true) = crate::sidereal::resolve_ayanamsha_flags(flags);
    let (mean, _) = aya.mean_tt_with(tt, noaberr, nogdefl)?;
    if use_true {
        Ok(mean + crate::erfa::nut06a(tt)?.0.to_degrees())
    } else {
        Ok(mean)
    }
}

/// True-ecliptic-of-date matrix for a scalar TT observation epoch.
///
/// Nutation is IAU 2000A at TT; P 03 precession and mean obliquity use
/// the corresponding scalar TDB epoch. The equinox matrix includes frame
/// bias and nutation, followed by the true-obliquity ecliptic rotation.
/// Returns a row-major ICRF-to-ecliptic rotation; non-finite epochs and
/// underlying frame errors propagate. This is the registered-SPK route,
/// not the separate ERFA longitude-strip or scalar-JD sidereal route.
pub(crate) fn moon_spk_frame(tt: f64) -> Result<[[f64; 3]; 3]> {
    use crate::positions::reference_frames::{
        FRAME_ASEC2RAD, FRAME_TENTH_USEC_2_RAD, ecliptic_date_matrix, equinox_matrix,
        iau2000a_tenth_microarcsec, mean_obliquity_arcsec, nutation_matrix, precession_matrix,
    };
    let tdb = crate::time::tt_to_tdb_jd(tt)?;
    let (dpsi, deps) = iau2000a_tenth_microarcsec(tt)?;
    let dpsi = dpsi * FRAME_TENTH_USEC_2_RAD;
    let deps = deps * FRAME_TENTH_USEC_2_RAD;
    let eps = mean_obliquity_arcsec(tdb)? * FRAME_ASEC2RAD;
    let nut = nutation_matrix(eps, eps + deps, dpsi)?;
    let equinox = equinox_matrix(nut, precession_matrix(tdb)?)?;
    ecliptic_date_matrix(equinox, eps, deps)
}
