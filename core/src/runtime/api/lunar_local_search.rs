//! Locally visible lunar eclipses over owned global, site and horizon channels.
use super::solar_circumstances::{PositionResolver, SelectedMain, solar_position, solar_selected};
use super::*;
use crate::events::{self as lunar, LunarLocalEphemeris};

/// Lunar type/visibility bits, ten UT1 phase/horizon times and twenty attributes.
pub type LunarLocalResult = lunar::LunarLocalSearchResult;

impl Session {
    /// Find the next or previous lunar eclipse visible at east longitude and
    /// north latitude in degrees, WGS84 height in metres, from original UT1 JD.
    /// One-try is removed; scalar shadow/global positions retain the caller's
    /// ephemeris flags. Local apparent Sun/Moon rays use the actual observer
    /// surface and independent bundled UT1/TT rotation clocks. Moon altitude
    /// uses configured sidereal terms and the independent ecliptic frame.
    ///
    /// Returns type and per-phase visibility bits; ten UT1 Julian days hold
    /// maximum, reserved zero, partial begin/end, total begin/end, penumbral
    /// begin/end, and Moon rise/set inside the eclipse. Hidden contacts become
    /// zero; a horizon crossing can replace maximum. Twenty attributes use the
    /// reported maximum: magnitudes0/1, azimuth4, true/apparent altitude5/6,
    /// opposition distance7, repeated umbral magnitude8, lunar Saros9/10;
    /// remaining slots are zero. Local phase refraction estimates pressure and
    /// uses 10 Celsius; rise/set requests use zero pressure/temperature.
    ///
    /// At most 250 candidates are considered, each with a 0.35-day offset against
    /// direction. Apparent altitude greater than zero sets visibility. The
    /// local window and geocentric maximum apply the epoch margin before final
    /// how. Every consumed global/phase/rise/set/final failure propagates with
    /// its owned carrier. Only exhausted visibility search creates its own base
    /// error. No early clock/image/site check or caller-observer publication is
    /// introduced. Explicit LEB/Auto retains the checked sealed serving.
    /// The supplied geometry requires finite, nondegenerate vectors. This
    /// method accepts fixed Rust scalar/site values and returns owned arrays.
    pub fn lun_eclipse_when_loc(
        &self,
        jd_ut: f64,
        site: [f64; 3],
        flags: i32,
        backwards: bool,
    ) -> Result<LunarLocalResult, Error> {
        self.lun_eclipse_when_loc_inner(jd_ut, site, flags, backwards, None)
    }

    /// Resolve missing native-main configuration at its consumed public body
    /// boundary. The callback borrows a snapshot synchronously and returns an
    /// owned context, never computed words. Images and registry owners survive
    /// through each call. Sealed/Auto never consumes this callback; settings
    /// and observers are not published back into the caller's context.
    pub fn lun_eclipse_when_loc_with_position_context(
        &self,
        jd_ut: f64,
        site: [f64; 3],
        flags: i32,
        backwards: bool,
        resolve: &PositionResolver<'_>,
    ) -> Result<LunarLocalResult, Error> {
        self.lun_eclipse_when_loc_inner(jd_ut, site, flags, backwards, Some(resolve))
    }

    /// Select the bound provider without acquiring a native image, clock or
    /// observer before the actual global-candidate callback.
    fn lun_eclipse_when_loc_inner(
        &self,
        jd_ut: f64,
        site: [f64; 3],
        flags: i32,
        backwards: bool,
        resolve: Option<&PositionResolver<'_>>,
    ) -> Result<LunarLocalResult, Error> {
        let flags = lunar::strip_one_try_bit(flags);
        if self.ctx.mode() != ProviderMode::NativeJpl {
            return self.with_serving_stack(|stack| {
                let dt = |ut| crate::time::deltat_with_config(ut, stack.delta_cfg);
                let serving = lunar::SealedLunarLocalServing::try_new(
                    stack.serving,
                    &dt,
                    flags,
                    site[0],
                    site[1],
                    site[2],
                )?;
                lunar::lun_eclipse_when_loc_sealed(&serving, jd_ut, backwards)
            });
        }
        let native = NativeLunarLocal {
            session: self,
            site,
            flags,
            resolve,
            selected: SelectedMain::new(),
        };
        lunar::lun_eclipse_when_loc_source_via(&native, jd_ut, backwards, &|| {
            local_exhausted(jd_ut, site)
        })
    }
}

/// One local walk borrows immutable caller settings and owns its acquired local
/// main snapshot. Global search retains its own accepted acquisition channels;
/// Phase and horizon requests borrow immutable native configuration.
struct NativeLunarLocal<'a> {
    session: &'a Session,
    site: [f64; 3],
    flags: i32,
    resolve: Option<&'a PositionResolver<'a>>,
    selected: SelectedMain,
}
impl LunarLocalEphemeris for NativeLunarLocal<'_> {
    fn next_global(&self, jd: f64, backwards: bool) -> Result<(i32, [f64; 10]), Error> {
        match self.resolve {
            Some(resolve) => self
                .session
                .lun_eclipse_when_with_position_context(jd, self.flags, 0, backwards, resolve),
            None => self.session.lun_eclipse_when(jd, self.flags, 0, backwards),
        }
    }

    fn how_attr(&self, jd: f64) -> Result<(i32, [f64; 20]), Error> {
        let flags = lunar::ecl_eph_flags(self.flags) | FLG_EQUATORIAL | FLG_XYZ | FLG_SPEED;
        let moon = solar_position(self.session, jd, MOON, flags, self.resolve, &self.selected)?;
        let sun = solar_position(self.session, jd, SUN, flags, self.resolve, &self.selected)?;
        let shadow = lunar::lun_how_core_via(
            (moon[0], moon[1], moon[2]),
            (sun[0], sun[1], sun[2]),
            lunar::saros_info(jd, lunar::EclipseKind::Lunar),
        )?;
        let owned = solar_selected(self.session, SUN, self.resolve, &self.selected)?;
        let (_sun, moon) =
            crate::positions::context_serving::native_lunar_eclipse_topocentric_pair(
                owned.context(),
                jd,
                self.site,
            )?;
        let horizontal = self.session.azalt(
            jd,
            lunar::ECL2HOR,
            self.site[0],
            self.site[1],
            self.site[2],
            0.0,
            10.0,
            moon[0],
            moon[1],
        )?;
        Ok(lunar::lun_eclipse_how_parts(
            &shadow,
            (
                horizontal.azimuth,
                horizontal.true_altitude,
                horizontal.apparent_altitude,
            ),
        ))
    }

    fn moon_horizon_crossings(&self, jd: f64) -> Result<lunar::BodyHorizonCrossings, Error> {
        let owned = solar_selected(self.session, MOON, self.resolve, &self.selected)?;
        moon_crossings_via(
            jd,
            self.site,
            lunar::ecl_eph_flags(self.flags),
            &|jd, rsmi, site, flags| owned.rise_trans(jd, MOON, rsmi, site, 0.0, 0.0, flags),
        )
    }
}

/// Lower-limb Moon rise before set, both from the same original UT1 and site.
/// A failure prevents later calls; negative circumpolar statuses remain values.
type MoonRiseCall<'a> = dyn Fn(f64, i32, [f64; 3], i32) -> Result<(i32, [f64; 10]), Error> + 'a;

/// Acquire the two lower-limb events through the same synchronous provider.
fn moon_crossings_via(
    jd: f64,
    site: [f64; 3],
    flags: i32,
    rise: &MoonRiseCall<'_>,
) -> Result<lunar::BodyHorizonCrossings, Error> {
    let (rise_rc, rises) = rise(jd, CALC_RISE | BIT_DISC_BOTTOM, site, flags)?;
    let (set_rc, sets) = rise(jd, CALC_SET | BIT_DISC_BOTTOM, site, flags)?;
    Ok(lunar::BodyHorizonCrossings {
        rise_rc,
        rise_jd: rises[0],
        set_rc,
        set_jd: sets[0],
    })
}

/// Native exhausted visibility-search error containing the observer coordinates.
fn local_exhausted(jd: f64, site: [f64; 3]) -> Error {
    Error::EclipseSearchExhausted {
        entry: "lun_eclipse_when_loc",
        detail: format!(
            "No lunar eclipse visible from longitude {}, latitude {} within the search range of JD {jd}",
            site[0], site[1]
        ),
    }
}

#[cfg(test)]
mod tests;
