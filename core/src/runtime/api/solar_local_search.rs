//! Local solar contacts over an owned observer and consumed public channels.
use super::solar_circumstances::{
    PositionResolver, SelectedMain, solar_horizontal, solar_selected,
};
use super::*;
use crate::events::{self as solar, SolarLocalEphemeris};

/// Local solar character, ten UT1 contact/maximum/horizon times and twenty
/// circumstance words. Absent contacts and reserved times are positive zero.
pub type SolarLocalResult = (i32, [f64; 10], [f64; 20]);

impl Session {
    /// Next visible local solar eclipse from UT1 Julian Day at east longitude,
    /// north latitude degrees and WGS84 height metres. `backwards` selects the
    /// backward collection arm; there is no local type mask or bidirectional
    /// public argument. One-try is stripped; caller flags select only global
    /// scalar geometry/horizon ephemeris bits. Direct local Topos rays ignore
    /// output/correction flags, using independent bundled UT1/TT/ITRS clocks.
    ///
    /// Returns character/visibility bits, ten UT1 JDs (reported maximum, four
    /// geometric contacts, sunrise/sunset within contacts, three zeros), and
    /// twenty attributes. Apparent altitude >0 sets visibility. When the
    /// geometric maximum is hidden, a horizon crossing can become the reported
    /// maximum and its character is recomputed. Attributes use that reported
    /// instant; ground shadow width and NONCENTRAL come from one global cone.
    /// Observer settings are never published into the caller's context.
    ///
    /// The global candidate uses the bidirectional search default. A candidate
    /// epoch has a 0.35-day selection window; maximum refinement uses +/-3 hours,
    /// outer contacts +/-2 hours then +/-4 hours, and inner contacts +/-0.5 hour.
    /// The maximum resolves the representable UT1 clock; contact bisections
    /// retain their sixty-iteration limit. Exhaustion is a no-event error after
    /// 120 outer iterations. Backward collection can inspect more than one
    /// global candidate inside an outer iteration. Global
    /// Configuration and coverage failures are caught only at that stage;
    /// local geometry and clock failures propagate. Explicit LEB/Auto retains
    /// its configured serving.
    pub fn sol_eclipse_when_loc(
        &self,
        jd_ut: f64,
        site: [f64; 3],
        flags: i32,
        backwards: bool,
    ) -> Result<SolarLocalResult, Error> {
        self.sol_eclipse_when_loc_inner(jd_ut, site, flags, backwards, None)
    }

    /// Configure a missing native main only at a consumed acquisition. The
    /// resolver returns owned settings/registry, not time or astronomical words.
    /// Its selected snapshot survives all global/local/contact/horizon reads.
    pub fn sol_eclipse_when_loc_with_position_context(
        &self,
        jd_ut: f64,
        site: [f64; 3],
        flags: i32,
        backwards: bool,
        resolve: &PositionResolver<'_>,
    ) -> Result<SolarLocalResult, Error> {
        self.sol_eclipse_when_loc_inner(jd_ut, site, flags, backwards, Some(resolve))
    }

    /// Choose only the bound provider family; no native clock, observer or
    /// image is acquired ahead of the actual global/local callback.
    fn sol_eclipse_when_loc_inner(
        &self,
        jd_ut: f64,
        site: [f64; 3],
        flags: i32,
        backwards: bool,
        resolve: Option<&PositionResolver<'_>>,
    ) -> Result<SolarLocalResult, Error> {
        let flags = solar::strip_one_try_bit(flags);
        if self.ctx.mode() != ProviderMode::NativeJpl {
            return self.with_serving_stack(|stack| {
                let dt = |ut| crate::time::deltat_with_config(ut, stack.delta_cfg);
                let serving = solar::SealedSolarLocalServing::try_new(
                    stack.serving,
                    &dt,
                    flags,
                    site[0],
                    site[1],
                    site[2],
                )?;
                solar::sol_eclipse_when_loc_sealed(&serving, jd_ut, backwards)
            });
        }
        let native = NativeSolarLocal {
            session: self,
            site,
            flags,
            resolve,
            selected: SelectedMain::new(),
        };
        solar::sol_eclipse_when_loc_source_via(
            &native,
            jd_ut,
            backwards,
            &caught_global,
            &|active| local_exhausted(jd_ut, site, active),
        )
    }
}

/// Read-only per-call observer acquisition, with no mutable borrow across any
/// provider or configuration callback. Local rays retain the real WGS84
/// surface, independently of the site-free gamma channel in the global search.
struct NativeSolarLocal<'a> {
    session: &'a Session,
    site: [f64; 3],
    flags: i32,
    resolve: Option<&'a PositionResolver<'a>>,
    selected: SelectedMain,
}
impl NativeSolarLocal<'_> {
    /// Direct Sun then Moon apparent true-ecliptic-of-date degrees/AU. The
    /// existing Topos helper retains Some(actual_surface) deflection semantics.
    fn pair(&self, jd: f64) -> Result<([f64; 3], [f64; 3]), Error> {
        let owned = solar_selected(self.session, SUN, self.resolve, &self.selected)?;
        crate::positions::context_serving::native_lunar_eclipse_topocentric_pair(
            &owned.ctx, jd, self.site,
        )
    }
}
impl SolarLocalEphemeris for NativeSolarLocal<'_> {
    fn next_global_forward(&self, jd: f64) -> Result<(i32, [f64; 10]), Error> {
        // This trait's legacy method name does not select direction. The local
        // public call's private candidate keeps its bidirectional default.
        let retain = |_: &EphemerisContext, body| {
            solar_selected(self.session, body, self.resolve, &self.selected)
                .map(|selected| selected.ctx.clone())
        };
        self.session.sol_eclipse_when_glob_direction_inner(
            jd,
            self.flags,
            0,
            solar::EclipseSearchDirection::Bidirectional,
            Some(&retain),
        )
    }
    fn separation_deg(&self, jd: f64) -> Result<f64, Error> {
        let (sun, moon) = self.pair(jd)?;
        solar::vincenty_separation_deg(sun[0], sun[1], moon[0], moon[1])
    }
    fn radii_sep_deg(&self, jd: f64, inner: bool) -> Result<(f64, f64, f64), Error> {
        let (sun, moon) = self.pair(jd)?;
        let sun_radius = radius_asin(solar::ECL_RSUN_AU / sun[2])?;
        let radius = if inner {
            solar::ECL_RMOON_INNER_AU
        } else {
            solar::ECL_RMOON_AU
        };
        let moon_radius = radius_asin(radius / moon[2])?;
        let separation = solar::vincenty_separation_deg(sun[0], sun[1], moon[0], moon[1])?;
        Ok((sun_radius, moon_radius, separation))
    }
    fn sun_altaz_deg(&self, jd: f64) -> Result<(f64, f64, f64), Error> {
        let (sun, _) = self.pair(jd)?;
        let (az, true_alt, apparent_alt) = solar_horizontal(self.session, jd, self.site, sun)?;
        Ok((true_alt, az, apparent_alt))
    }
    fn topo_diameters_deg(&self, jd: f64) -> Result<(f64, f64), Error> {
        let (sun, moon) = self.pair(jd)?;
        Ok((
            2.0 * radius_asin(solar::ECL_RMOON_AU / moon[2])?,
            2.0 * radius_asin(solar::ECL_RSUN_AU / sun[2])?,
        ))
    }
    fn sun_horizon_crossings(&self, jd: f64) -> Result<solar::BodyHorizonCrossings, Error> {
        let owned = solar_selected(self.session, SUN, self.resolve, &self.selected)?;
        let flags = self.flags & (FLG_JPLEPH | FLG_SWIEPH | FLG_MOSEPH);
        let (rise_rc, rise) = owned.rise_trans(
            jd,
            SUN,
            crate::constants::CALC_RISE | crate::constants::BIT_DISC_BOTTOM,
            self.site,
            0.0,
            0.0,
            flags,
        )?;
        let (set_rc, set) = owned.rise_trans(
            jd,
            SUN,
            crate::constants::CALC_SET | crate::constants::BIT_DISC_BOTTOM,
            self.site,
            0.0,
            0.0,
            flags,
        )?;
        Ok(solar::BodyHorizonCrossings {
            rise_rc,
            rise_jd: rise[0],
            set_rc,
            set_jd: set[0],
        })
    }
    fn how_core(&self, jd: f64) -> Result<(i32, [f64; 20]), Error> {
        let how = self.session.native_solar_local(
            jd,
            self.site,
            self.flags,
            self.resolve,
            &self.selected,
        )?;
        Ok((how.retc(), how.attr()))
    }
    fn where_core(&self, jd: f64) -> Result<(i32, [f64; 10]), Error> {
        let geometry =
            self.session
                .native_solar_global(jd, self.flags, self.resolve, &self.selected)?;
        Ok((geometry.retc(), geometry.dcore()))
    }
}

/// Apparent angular radius in degrees from a finite radius/distance ratio.
/// Ratios outside [-1, 1] return the native mathematical domain error.
fn radius_asin(ratio: f64) -> Result<f64, Error> {
    if !ratio.is_finite() || !(-1.0..=1.0).contains(&ratio) {
        return Err(Error::MathDomain { value: ratio });
    }
    Ok(ratio.asin().to_degrees())
}

/// Native configuration, body and ephemeris failures end a global candidate
/// walk. Numerical-domain and arithmetic failures remain visible to the caller.
fn caught_global(error: &Error) -> bool {
    matches!(
        error.primary_error(),
        Error::ConfigurationMissing { .. }
            | Error::NetworkSealed { .. }
            | Error::EphemerisRange { .. }
            | Error::UnknownBody { .. }
            | Error::UnsupportedFlagCombination { .. }
            | Error::InvalidEclipseSearch { .. }
            | Error::EclipseSearchExhausted { .. }
    )
}

/// Visibility-search exhaustion with an optional directly owned native cause.
fn local_exhausted(jd: f64, site: [f64; 3], active: Option<Error>) -> Error {
    let error = Error::EclipseSearchExhausted {
        entry: "sol_eclipse_when_loc",
        detail: format!(
            "No solar eclipse visible from longitude {}, latitude {} within 50 years of JD {jd}",
            site[0], site[1]
        ),
    };
    match active {
        Some(cause) => error.with_cause(cause),
        None => error,
    }
}

#[cfg(test)]
mod tests;
