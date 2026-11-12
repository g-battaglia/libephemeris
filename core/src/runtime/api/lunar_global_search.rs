//! Global lunar eclipses over distinct public UT1 position acquisitions.
use super::solar_circumstances::{PositionResolver, SelectedMain, solar_position};
use super::*;
use crate::events::{self as lunar, LunarGlobalEphemeris, SunMoonState};

impl Session {
    /// Next or previous global lunar eclipse at an original UT1 Julian day.
    /// Flags select ephemerides; mask0 accepts any lunar type. Returns character
    /// bits and ten UT1 JDs: maximum, reserved0, partial begin/end, total
    /// begin/end, penumbral begin/end and two reserved zeros. Missing phases
    /// are0. Type normalization precedes native image acquisition. The default
    /// SPEED Sun/Moon full-moon walk, caller-flags Moon-only node prefilter,
    /// independently clocked mean node, and Cartesian depth/shadow acquisition
    /// remain separate. Maximum uses a0.3-day half-window, sixty golden steps
    /// and1e-7-day convergence; contacts use four-hour brackets and sixty
    /// bisections with1e-8-day convergence. Epoch-margin and forward preprobe
    /// ordering are preserved. Boundary exhaustion has its own base Error;
    /// other numerical and provider failures propagate as native errors.
    /// The first default position acquisition owns epoch admission; there is
    /// no preceding replacement of nonfinite epochs. Supplied registries and lazily acquired images
    /// remain owned until return; caller observer/settings are never changed.
    /// LEB/Auto retain their checked sealed serving and flag/domain gates.
    /// This search does not certify independent physical accuracy.
    pub fn lun_eclipse_when(
        &self,
        jd_ut: f64,
        flags: i32,
        mask: i32,
        backwards: bool,
    ) -> Result<(i32, [f64; 10]), Error> {
        self.lun_eclipse_when_inner(jd_ut, flags, mask, backwards, None)
    }

    /// Global lunar search with a configuration-only native main resolver.
    /// The resolver is synchronous and consumed at the first actual public
    /// body acquisition when the supplied registry is empty. It returns owned
    /// context configuration, never a replacement position or clock. Arguments,
    /// outputs, errors and domains follow [`Self::lun_eclipse_when`]. Other
    /// provider modes never consume the resolver.
    pub fn lun_eclipse_when_with_position_context(
        &self,
        jd_ut: f64,
        flags: i32,
        mask: i32,
        backwards: bool,
        resolve: &PositionResolver<'_>,
    ) -> Result<(i32, [f64; 10]), Error> {
        self.lun_eclipse_when_inner(jd_ut, flags, mask, backwards, Some(resolve))
    }

    /// Own one native acquisition cache while retaining the sealed dispatcher.
    fn lun_eclipse_when_inner(
        &self,
        jd_ut: f64,
        flags: i32,
        mask: i32,
        backwards: bool,
        resolve: Option<&PositionResolver<'_>>,
    ) -> Result<(i32, [f64; 10]), Error> {
        let flags = lunar::strip_one_try_bit(flags);
        if self.ctx.mode() != ProviderMode::NativeJpl {
            return self.with_serving_stack(|stack| {
                let dt = |ut| crate::time::deltat_with_config(ut, stack.delta_cfg);
                let bound = lunar::SealedLunarGlobalServing::try_new(stack.serving, &dt, flags)?;
                lunar::lun_eclipse_when_glob_sealed(&bound, jd_ut, mask, backwards)
            });
        }
        let native = NativeLunarGlobal {
            session: self,
            flags,
            resolve,
            selected: SelectedMain::new(),
        };
        lunar::lun_eclipse_when_source_via(&native, jd_ut, mask, backwards)
    }
}

/// One synchronous search retains an immutable native image snapshot. Scalar
/// clocks and flags are those of each consumed public calc_ut call; the mean
/// node uses its separate bundled global Time TT at the same original UT1.
struct NativeLunarGlobal<'a> {
    session: &'a Session,
    flags: i32,
    resolve: Option<&'a PositionResolver<'a>>,
    selected: SelectedMain,
}
impl NativeLunarGlobal<'_> {
    /// Moon then Sun equatorial Cartesian AU, with SPEED only for the shadow
    /// model. The depth-minimum route deliberately omits that flag.
    fn xyz(&self, jd: f64, speed: bool) -> Result<lunar::MoonSunXyzAu, Error> {
        native_xyz_via(jd, self.flags, speed, &|ut, body, flags| {
            solar_position(self.session, ut, body, flags, self.resolve, &self.selected)
        })
    }
}
impl LunarGlobalEphemeris for NativeLunarGlobal<'_> {
    fn sun_moon_state(&self, jd: f64) -> Result<SunMoonState, Error> {
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
    fn geocentric_xyz_au(&self, jd: f64) -> Result<lunar::MoonSunXyzAu, Error> {
        self.xyz(jd, false)
    }
    fn lun_shadow(&self, jd: f64) -> Result<(i32, lunar::LunarShadowDcore), Error> {
        let (moon, sun) = self.xyz(jd, true)?;
        let result =
            lunar::lun_how_core_via(moon, sun, lunar::saros_info(jd, lunar::EclipseKind::Lunar))?;
        Ok((
            result.retc(),
            lunar::LunarShadowDcore::from_lun_how_dcore(result.dcore())?,
        ))
    }
}

/// Discard non-ephemeris caller flags for Cartesian shadow acquisition. Return
/// the first three AU words of each actual Moon then Sun position, preserving
/// original UT1 and the distinction between depth and canonical shadow flags.
fn native_xyz_via(
    jd: f64,
    flags: i32,
    speed: bool,
    position: &dyn Fn(f64, i32, i32) -> Result<[f64; 6], Error>,
) -> Result<lunar::MoonSunXyzAu, Error> {
    let base =
        lunar::ecl_eph_flags(flags) | FLG_EQUATORIAL | FLG_XYZ | if speed { FLG_SPEED } else { 0 };
    let moon = position(jd, MOON, base)?;
    let sun = position(jd, SUN, base)?;
    Ok(((moon[0], moon[1], moon[2]), (sun[0], sun[1], sun[2])))
}

#[cfg(test)]
mod tests;
