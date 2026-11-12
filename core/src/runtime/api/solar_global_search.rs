//! Global solar contact search over the public scalar, direct vector and frame
//! acquisition channels. No observer is installed in the caller's context.
use super::solar_circumstances::{PositionResolver, SelectedMain, solar_position, solar_selected};
use super::*;
use crate::events::{self as solar, EclipseSearchDirection, SolarGlobalEphemeris, SunMoonState};

impl Session {
    /// Find the next global solar eclipse, or the previous one when backwards
    /// is true, starting at UT1 Julian Day. Flags are the ephemeris bitmask;
    /// mask0 accepts any type. Return character bits and ten UT1 JDs:
    /// maximum, central noon, penumbral begin/end, umbral begin/end, central-line
    /// begin/end and two zero hybrid-transition slots. Missing contacts are0.
    /// Impossible type masks refuse before position acquisition. A maximum at
    /// the start is skipped using the existing epoch margin; the existing
    /// lunation horizon, golden brackets and contact bisections are retained.
    /// Body/configuration/clock errors propagate; exhausted coverage has the
    /// global search's no-event carrier. Supplied registries remain owned.
    pub fn sol_eclipse_when_glob(
        &self,
        jd_ut: f64,
        flags: i32,
        mask: i32,
        backwards: bool,
    ) -> Result<(i32, [f64; 10]), Error> {
        let direction = if backwards {
            EclipseSearchDirection::Backward
        } else {
            EclipseSearchDirection::Forward
        };
        self.sol_eclipse_when_glob_direction_inner(jd_ut, flags, mask, direction, None)
    }

    /// Configuration-only native main resolver. It is invoked at the first
    /// actual position/vector acquisition, without injecting time or position
    /// words. A supplied registry has priority and stays retained until return.
    pub fn sol_eclipse_when_glob_with_position_context(
        &self,
        jd_ut: f64,
        flags: i32,
        mask: i32,
        backwards: bool,
        resolve: &PositionResolver<'_>,
    ) -> Result<(i32, [f64; 10]), Error> {
        let direction = if backwards {
            EclipseSearchDirection::Backward
        } else {
            EclipseSearchDirection::Forward
        };
        self.sol_eclipse_when_glob_direction_inner(jd_ut, flags, mask, direction, Some(resolve))
    }

    /// Rust search-direction interface, including the bidirectional15-day
    /// probe. Boolean-direction methods select only forward/backward.
    /// Numerical limits, contact slots, clocks and errors match the scalar route.
    pub fn sol_eclipse_when_glob_direction(
        &self,
        jd_ut: f64,
        flags: i32,
        mask: i32,
        direction: EclipseSearchDirection,
    ) -> Result<(i32, [f64; 10]), Error> {
        self.sol_eclipse_when_glob_direction_inner(jd_ut, flags, mask, direction, None)
    }

    /// Retain the existing sealed serving and lazily acquire native words.
    /// Type-mask validation lives in the walk, before unused native resources.
    pub(super) fn sol_eclipse_when_glob_direction_inner(
        &self,
        jd_ut: f64,
        flags: i32,
        mask: i32,
        direction: EclipseSearchDirection,
        resolve: Option<&PositionResolver<'_>>,
    ) -> Result<(i32, [f64; 10]), Error> {
        let flags = solar::strip_one_try_bit(flags);
        if self.ctx.mode() != ProviderMode::NativeJpl {
            return self.with_serving_stack(|stack| {
                let dt = |ut| crate::time::deltat_with_config(ut, stack.delta_cfg);
                let bound = solar::SealedSolarGlobalServing::try_new(stack.serving, &dt, flags)?;
                solar::sol_eclipse_when_glob_sealed(&bound, jd_ut, mask, direction)
            });
        }
        let native = NativeSolarGlobal {
            session: self,
            flags,
            resolve,
            selected: SelectedMain::new(),
            besselian: std::cell::RefCell::new(std::collections::HashMap::new()),
        };
        solar::sol_eclipse_when_glob_source_via(&native, jd_ut, mask, direction)
    }
}

/// One search retains its selected immutable main. Default-speed scalar calls
/// and caller-flags-only Moon probes remain distinct; gamma uses raw ICRS
/// apparent rays at independent bundled UT1 Time, not public J2000 output.
struct NativeSolarGlobal<'a> {
    session: &'a Session,
    flags: i32,
    resolve: Option<&'a PositionResolver<'a>>,
    selected: SelectedMain,
    besselian: std::cell::RefCell<std::collections::HashMap<u64, f64>>,
}
impl NativeSolarGlobal<'_> {
    /// Geocentric Sun and Moon position/speed at the caller UT1 epoch.
    fn fixed_state(&self, jd: f64) -> Result<SunMoonState, Error> {
        let sun = solar_position(
            self.session,
            jd,
            SUN,
            FLG_SPEED,
            self.resolve,
            &self.selected,
        )?;
        let moon = solar_position(
            self.session,
            jd,
            MOON,
            FLG_SPEED,
            self.resolve,
            &self.selected,
        )?;
        SunMoonState::try_new(
            sun[0], sun[1], sun[2], sun[3], moon[0], moon[1], moon[2], moon[3],
        )
    }
    /// Global true-equatorial Moon then Sun geometry followed by public GAST.
    fn geometry(&self, jd: f64) -> Result<solar::WhereGeometry, Error> {
        self.session
            .native_solar_global(jd, FLG_SWIEPH, self.resolve, &self.selected)
    }
}
impl SolarGlobalEphemeris for NativeSolarGlobal<'_> {
    fn sun_moon_state(&self, jd: f64) -> Result<SunMoonState, Error> {
        self.fixed_state(jd)
    }
    fn candidate_moon_longitude_deg(&self, jd: f64) -> Result<f64, Error> {
        Ok(solar_position(
            self.session,
            jd,
            MOON,
            self.flags | FLG_SPEED,
            self.resolve,
            &self.selected,
        )?[0])
    }
    fn mean_node_lon_deg(&self, jd: f64) -> Result<f64, Error> {
        let tt = crate::time::default_timescale_tt_from_ut1(jd)?;
        crate::lunar::calc_mean_lunar_node(tt, crate::lunar::LunarCoverage::Unenforced)
    }
    fn gamma_er(&self, jd: f64) -> Result<f64, Error> {
        let selected = solar_selected(self.session, SUN, self.resolve, &self.selected)?;
        let (sun, moon) =
            crate::positions::context_serving::native_solar_gamma_pair(&selected.ctx, jd)?;
        solar::gamma_from_vectors((sun[0], sun[1], sun[2]), (moon[0], moon[1], moon[2]))
    }
    fn shadow_limits_er(&self, jd: f64) -> Result<(f64, f64), Error> {
        // Penumbral and umbral radii use independently acquired states.
        let pen = self.fixed_state(jd)?;
        let l1 = solar::penumbra_limit_er(pen.sun_dist_au(), pen.moon_dist_au())?;
        let umb = self.fixed_state(jd)?;
        let l2 = solar::umbra_limit_er(umb.sun_dist_au(), umb.moon_dist_au())?;
        Ok((l1, l2))
    }
    fn where_gate(&self, jd: f64) -> Result<(i32, f64), Error> {
        let geometry = self.geometry(jd)?;
        let site = [geometry.lon_deg(), geometry.lat_deg(), 0.];
        let local =
            self.session
                .native_solar_local(jd, site, FLG_SWIEPH, self.resolve, &self.selected)?;
        Ok((geometry.retc(), local.attr()[0]))
    }
    fn where_retflag(&self, jd: f64) -> Result<i32, Error> {
        Ok(self.geometry(jd)?.retc())
    }
    fn where_magnitude(&self, jd: f64) -> Result<f64, Error> {
        self.where_gate(jd).map(|pair| pair.1)
    }
    fn contact_residuals_km(&self, jd: f64) -> Result<(f64, f64, f64), Error> {
        let geom = self.geometry(jd)?;
        let r = solar::besselian_contact_residuals(geom.dcore())?;
        Ok((r.pen_km, r.umb_km, r.cen_km))
    }
    fn besselian_x_er(&self, jd: f64) -> Result<f64, Error> {
        let key = jd.to_bits();
        let cached = { self.besselian.borrow().get(&key).copied() };
        if let Some(value) = cached {
            return Ok(value);
        }
        let value = native_besselian_x(self.session, jd, self.resolve, &self.selected)?;
        let mut cache = self.besselian.borrow_mut();
        if cache.len() > 512 {
            cache.clear();
        }
        cache.insert(key, value);
        Ok(value)
    }
    fn ground_core_width_km(&self, jd: f64) -> Result<f64, Error> {
        Ok(self.geometry(jd)?.dcore()[0])
    }
}

/// Besselian Sun then Moon spherical true-equatorial default-ephemeris words;
/// scalar-to-Cartesian conversion occurs before public sidtime is consumed.
fn native_besselian_x(
    session: &Session,
    jd: f64,
    resolve: Option<&PositionResolver<'_>>,
    retained: &SelectedMain,
) -> Result<f64, Error> {
    let flags = FLG_SWIEPH | FLG_EQUATORIAL;
    let sun = solar_position(session, jd, SUN, flags, resolve, retained)?;
    let moon = solar_position(session, jd, MOON, flags, resolve, retained)?;
    let vector = |p: [f64; 6]| {
        let ra = p[0].to_radians();
        let dec = p[1].to_radians();
        let r = p[2];
        (
            r * dec.cos() * ra.cos(),
            r * dec.cos() * ra.sin(),
            r * dec.sin(),
        )
    };
    let sun = vector(sun);
    let moon = vector(moon);
    solar::besselian_core_lazy(sun, moon, &mut || session.sidtime(jd).map(|h| h * 15.0))
        .map(|b| b.x())
}

#[cfg(test)]
mod tests;
