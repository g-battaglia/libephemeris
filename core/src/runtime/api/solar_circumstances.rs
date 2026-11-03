//! Fixed-instant solar circumstances over owned public acquisition channels.
//! Caller UT1, local observer acquisition and site-free shadow geometry remain
//! separate. The returned arrays use the public twenty/ten-slot layouts.
use super::*;
use crate::events::{self as solar, SolHowCircumstances, WhereGeometry};

/// Synchronous configuration-only acquisition at a consumed native main-image
/// boundary. No body words, time values or provider borrow escape the callback.
pub(super) type PositionResolver<'a> =
    dyn Fn(&EphemerisContext, i32) -> Result<EphemerisContext, Error> + 'a;
/// Call-local main image retention. No callback runs under a cache borrow or
/// mutex, and the cell is dropped only after all local/global reads complete.
pub(super) type SelectedMain = std::cell::OnceCell<Session>;

impl Session {
    /// Solar eclipse circumstances at original UT1 Julian Day and a longitude
    /// east/latitude north degree, height-metre site. Returns character bits and
    /// twenty attributes: magnitude, diameter ratio, obscuration, global core
    /// shadow diameter km, south-westward azimuth, true/apparent altitude,
    /// separation degrees, NASA magnitude, Saros series/member, reserved zeros.
    /// Apparent Sun altitude <=0 clears character bits; a zero character clears
    /// slots0,1,2,3,8,9,10 and retains4..7. No event search is performed.
    ///
    /// Native local positions use the shared direct WGS84 Sun-then-Moon apparent
    /// pair at bundled UT1 Time, independently of caller output/correction flags.
    /// Local projection uses independent frame TT then configured sidereal
    /// clocks. Site-free shadow geometry follows: Moon then Sun public calc_ut
    /// at original UT1 with only ephemeris bits plus EQUATORIAL|XYZ, followed by
    /// public sidtime at the geometry's consumed rotation stage. One-try is
    /// stripped. Explicit LEB/Auto retain the unchanged sealed serving/guards.
    /// Body/image/observer/vector/clock errors propagate without category folds;
    /// nonfinite caller sites/epochs follow those actual acquisition guards.
    pub fn sol_eclipse_how(
        &self,
        jd_ut: f64,
        site: [f64; 3],
        flags: i32,
    ) -> Result<(i32, [f64; 20]), Error> {
        self.sol_eclipse_how_inner(jd_ut, site, flags, None)
    }

    /// Where the shadow axis intersects or comes closest to Earth at UT1 JD.
    /// Return global character, ten geoposition words (east longitude/north
    /// latitude degrees then eight zeros), and twenty local attributes at that
    /// actual point with height0. Shadow geometry precedes local Sun/Moon and
    /// horizon acquisition. Character0 is valid closest-approach data; there is
    /// no local horizon gate or zeroing. Units/clocks/errors follow how above.
    pub fn sol_eclipse_where(
        &self,
        jd_ut: f64,
        flags: i32,
    ) -> Result<(i32, [f64; 10], [f64; 20]), Error> {
        self.sol_eclipse_where_inner(jd_ut, flags, None)
    }

    /// How with configuration resolved only at consumed Native main-image
    /// boundaries. Supplied registries retain priority; no context observer or
    /// global setting is mutated, and all callbacks remain call-scoped.
    pub fn sol_eclipse_how_with_position_context(
        &self,
        jd_ut: f64,
        site: [f64; 3],
        flags: i32,
        resolve: &PositionResolver<'_>,
    ) -> Result<(i32, [f64; 20]), Error> {
        self.sol_eclipse_how_inner(jd_ut, site, flags, Some(resolve))
    }

    /// Configuration-only counterpart of where. Global public position and
    /// independent site-vector acquisition remain their distinct channels.
    pub fn sol_eclipse_where_with_position_context(
        &self,
        jd_ut: f64,
        flags: i32,
        resolve: &PositionResolver<'_>,
    ) -> Result<(i32, [f64; 10], [f64; 20]), Error> {
        self.sol_eclipse_where_inner(jd_ut, flags, Some(resolve))
    }

    /// Retain sealed acquisition and native stage order without an eager clock
    /// or observer preflight. Native callbacks resolve only what they consume.
    fn sol_eclipse_how_inner(
        &self,
        jd: f64,
        site: [f64; 3],
        flags: i32,
        resolve: Option<&PositionResolver<'_>>,
    ) -> Result<(i32, [f64; 20]), Error> {
        if self.ctx.mode() != ProviderMode::NativeJpl {
            return self.with_serving_stack(|stack| {
                let dt = |ut| crate::time::deltat_with_config(ut, stack.delta_cfg);
                let bound = solar::SealedSolarLocalServing::try_new(
                    stack.serving,
                    &dt,
                    solar::strip_one_try_bit(flags),
                    site[0],
                    site[1],
                    site[2],
                )?;
                solar::sol_how_sealed(&bound, jd)
            });
        }
        let selected = SelectedMain::new();
        solar_how_via(
            jd,
            site,
            flags,
            &mut |ut, gp, f| self.native_solar_local(ut, gp, f, resolve, &selected),
            &mut |ut, f| self.native_solar_global(ut, f, resolve, &selected),
        )
    }

    /// The stored where route remains the original global/center local serving;
    /// native geometry and actual-center local sampling are consumed in order.
    fn sol_eclipse_where_inner(
        &self,
        jd: f64,
        flags: i32,
        resolve: Option<&PositionResolver<'_>>,
    ) -> Result<(i32, [f64; 10], [f64; 20]), Error> {
        if self.ctx.mode() != ProviderMode::NativeJpl {
            return self.with_serving_stack(|stack| {
                let dt = |ut| crate::time::deltat_with_config(ut, stack.delta_cfg);
                let bound = solar::SealedSolarGlobalServing::try_new(
                    stack.serving,
                    &dt,
                    solar::strip_one_try_bit(flags),
                )?;
                solar::sol_eclipse_where_sealed(&bound, jd)
            });
        }
        let selected = SelectedMain::new();
        solar_where_via(
            jd,
            flags,
            &mut |ut, f| self.native_solar_global(ut, f, resolve, &selected),
            &mut |ut, gp, f| self.native_solar_local(ut, gp, f, resolve, &selected),
        )
    }

    /// Public Moon then Sun Cartesian true-equatorial positions; caller UT1
    /// remains exact. SPEED and unrelated flags are not inserted. Sidereal time
    /// is acquired lazily after the two-pass cone/surface calculation.
    pub(super) fn native_solar_global(
        &self,
        jd: f64,
        flags: i32,
        resolve: Option<&PositionResolver<'_>>,
        selected: &SelectedMain,
    ) -> Result<WhereGeometry, Error> {
        let base = (flags & (FLG_JPLEPH | FLG_SWIEPH | FLG_MOSEPH)) | FLG_EQUATORIAL | FLG_XYZ;
        let moon = solar_position(self, jd, MOON, base, resolve, selected)?;
        let sun = solar_position(self, jd, SUN, base, resolve, selected)?;
        solar::eclipse_where_core_lazy(
            (moon[0], moon[1], moon[2]),
            (sun[0], sun[1], sun[2]),
            &mut || self.sidtime(jd).map(|h| h * 15.0),
            solar::ECL_RSUN_AU,
        )
    }

    /// Direct local vector acquisition ignores caller flags. Disc geometry
    /// precedes the horizontal clock/refraction stage, then Saros acquisition.
    pub(super) fn native_solar_local(
        &self,
        jd: f64,
        site: [f64; 3],
        _flags: i32,
        resolve: Option<&PositionResolver<'_>>,
        selected: &SelectedMain,
    ) -> Result<SolHowCircumstances, Error> {
        let context = solar_selected(self, SUN, resolve, selected)?;
        let (sun, moon) = crate::positions::context_serving::native_lunar_eclipse_topocentric_pair(
            &context.ctx,
            jd,
            site,
        )?;
        solar::sol_how_core_lazy(
            (sun[0], sun[1], sun[2]),
            (moon[0], moon[1], moon[2]),
            solar::ECL_RSUN_AU,
            true,
            false,
            &mut || solar_horizontal(self, jd, site, sun),
            site[2],
            &mut || Ok(solar::saros_info(jd, solar::EclipseKind::Solar)),
        )
    }
}

/// Preserve a supplied registry. A missing binding configuration is acquired
/// only for an empty Native registry; returned images/settings stay owned.
fn solar_context(
    context: &EphemerisContext,
    body: i32,
    resolve: Option<&PositionResolver<'_>>,
) -> Result<EphemerisContext, Error> {
    if context.spk_registry().is_empty() {
        if let Some(resolve) = resolve {
            let selected = resolve(context, body)?;
            if selected.mode() != ProviderMode::NativeJpl {
                return Err(Error::ConfigurationMissing {
                    missing_config: "native_solar_provider",
                    suggestion: "preserve NativeJpl in the position configuration resolver",
                });
            }
            return Ok(selected);
        }
    }
    Ok(context.clone())
}

/// Validate actual public flags before resolving otherwise-unused paths; the
/// source's centre/ephemeris priorities live in the shared validator. Public
/// calc_ut retains its genuine nonfinite/scalar range constructor and clocks.
pub(super) fn solar_position(
    session: &Session,
    jd: f64,
    body: i32,
    flags: i32,
    resolve: Option<&PositionResolver<'_>>,
    selected: &SelectedMain,
) -> Result<[f64; 6], Error> {
    if !jd.is_finite() || crate::positions::validate_barycentric_moseph(flags, "calc_ut").is_err() {
        return session.calc_ut(jd, body, flags).map(|p| p.values);
    }
    solar_selected(session, body, resolve, selected)?
        .calc_ut(jd, body, flags)
        .map(|p| p.values)
}

/// Acquire the actual main once at its first consumed boundary, retaining every
/// registered mapping and shared immutable image through the complete query.
/// A caller-provided registry is used directly without a cache clone or loader.
pub(super) fn solar_selected<'a>(
    session: &'a Session,
    body: i32,
    resolve: Option<&PositionResolver<'_>>,
    selected: &'a SelectedMain,
) -> Result<&'a Session, Error> {
    if !session.ctx.spk_registry().is_empty() {
        return Ok(session);
    }
    if selected.get().is_none() {
        let context = solar_context(&session.ctx, body, resolve)?;
        let context = crate::positions::context_serving::native_eclipse_main_context(&context)?;
        let _ = selected.set(Session::from_context(context));
    }
    Ok(selected.get().expect("consumed main acquisition retained"))
}

/// Ecliptic-to-horizontal projection with frame TT before configured sidereal
/// acquisition. Longitude is added before modulo24; pressure is estimated and
/// temperature10C. Existing cotrans/refraction kernels retain their units and
/// arithmetic. This private finite acquisition does not use checked Houses
/// continuity-offset preflight or caller model TT as the frame epoch.
pub(super) fn solar_horizontal(
    session: &Session,
    jd: f64,
    site: [f64; 3],
    place: [f64; 3],
) -> Result<(f64, f64, f64), Error> {
    let frame_tt = crate::time::default_timescale_tt_from_ut1(jd)?;
    let eps0 = crate::positions::vondrak_mean_obliquity_deg(frame_tt)?;
    let (dpsi, deps) = crate::erfa::nut06a(frame_tt)?;
    let eps = eps0 + deps.to_degrees();
    let (ra, dec, _) = solar::cotrans(place[0], place[1], place[2], -eps);
    let hours = crate::time::public_sidereal::hours_with_longitude(
        jd,
        site[0],
        eps,
        dpsi.to_degrees(),
        &|ut| session.deltat(ut),
    )?;
    let ha = (hours * 15.0 - ra).rem_euclid(360.0).to_radians();
    let dec = dec.to_radians();
    let lat = site[1].to_radians();
    let sin_alt = lat.sin() * dec.sin() + lat.cos() * dec.cos() * ha.cos();
    // Source max(-1,min(1,x)) retains its first operand when x is NaN.
    let limited = if sin_alt < 1.0 { sin_alt } else { 1.0 };
    let limited = if limited > -1.0 { limited } else { -1.0 };
    let true_alt = limited.asin().to_degrees();
    let y = ha.sin() * dec.cos();
    let x = ha.cos() * dec.cos() * lat.sin() - dec.sin() * lat.cos();
    let az = y.atan2(x).to_degrees().rem_euclid(360.0);
    let pressure = crate::refraction::estimate_pressure_mbar(site[2]);
    let lapse = session
        .ctx
        .lapse_rate_override()
        .unwrap_or(solar::RISE_SET_DIP_LAPSE_RATE_K_M);
    let (app, _) =
        crate::refraction::refrac_extended_true_to_app(true_alt, site[2], pressure, 10.0, lapse);
    Ok((az, true_alt, app))
}

/// Compose the local core before the independent global geometry at original
/// UT1. Horizon/no-event zeroing belongs to the existing composition kernel.
fn solar_how_via(
    jd_ut: f64,
    site: [f64; 3],
    flags: i32,
    local: &mut dyn FnMut(f64, [f64; 3], i32) -> Result<SolHowCircumstances, Error>,
    global: &mut dyn FnMut(f64, i32) -> Result<WhereGeometry, Error>,
) -> Result<(i32, [f64; 20]), Error> {
    let flags = solar::strip_one_try_bit(flags);
    let how = local(jd_ut, site, flags)?;
    let geometry = global(jd_ut, flags)?;
    Ok(solar::sol_eclipse_how_parts(&how, &geometry))
}

/// Acquire site-free shadow geometry before local circumstances at the actual
/// shadow center, with zero height. The local character is discarded; the
/// global result stays valid when the local apparent Sun is below the horizon.
fn solar_where_via(
    jd_ut: f64,
    flags: i32,
    global: &mut dyn FnMut(f64, i32) -> Result<WhereGeometry, Error>,
    local: &mut dyn FnMut(f64, [f64; 3], i32) -> Result<SolHowCircumstances, Error>,
) -> Result<(i32, [f64; 10], [f64; 20]), Error> {
    let flags = solar::strip_one_try_bit(flags);
    let geometry = global(jd_ut, flags)?;
    let how = local(jd_ut, [geometry.lon_deg(), geometry.lat_deg(), 0.0], flags)?;
    Ok(solar::sol_eclipse_where_parts(&geometry, &how))
}

#[cfg(test)]
#[path = "solar_circumstance_tests.rs"]
mod tests;
