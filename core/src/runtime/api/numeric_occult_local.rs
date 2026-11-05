//! Numeric local lunar occultations over one immutable main and observer.
use super::numeric_occultation::{NativeNumericOccultWhere, numeric_occultation_target};
use super::solar_circumstances::PositionResolver;
use super::*;
use crate::events::{self as event, OccultLocalEphemeris, WhereGeometry};

/// Character, ten UT1-JD contact/horizon words and twenty local attributes.
pub type NumericOccultLocalResult = (i32, [f64; 10], [f64; 20]);

impl Session {
    /// Search for a numeric lunar occultation visible at east longitude,
    /// north latitude degrees and WGS84 height metres. Negative bodies select
    /// Sun; Moon-self refuses first. Flags default at bindings to2; bit0 in
    /// backwards selects earlier events and ECL_ONE_TRY examines one conjunction.
    /// There is no local type mask or atmospheric argument. Body/Sun horizon
    /// calls use zero pressure/temperature and retain their normal refraction.
    ///
    /// Outputs are local phase/visibility/daylight/global NONCENTRAL bits,
    /// ten UT1 dates (maximum, first..fourth contact, body rise/set, three
    /// zeros) and twenty circumstances. One-try misses return character0,
    /// continuation JD in slot0, and zero attributes. Topocentric target
    /// precedes a fresh Sun/Moon pair at every radius sample. Inner contacts
    /// use the smaller lunar radius; each visible-contact altitude is a fresh
    /// local circumstance acquisition. Final local attributes precede ground
    /// geometry. Model/scalar clocks and independent frame/Topos clocks stay
    /// distinct; the observer and selected images survive the full call.
    ///
    /// Scalar acquisition owns epoch admission. Conjunction/maximum/contact
    /// walks keep300/60/60 limits and their fixed width tolerances. First-target
    /// coverage exhaustion becomes an owned base Error; all other acquired
    /// failures preserve their channels. Sealed/LEB modes enforce their checked epoch/flag/site domain.
    /// No caller observer is published. String/catalog targets use the separate
    /// stellar entry point.
    pub fn lun_occult_when_loc(
        &self,
        jd: f64,
        body: i32,
        site: [f64; 3],
        flags: i32,
        backwards: i32,
    ) -> Result<NumericOccultLocalResult, Error> {
        self.numeric_occult_local(jd, body, site, flags, backwards, None)
    }

    /// Configuration-only native main acquisition at its first consumed
    /// boundary. Registered images retain priority; no position words are
    /// supplied by the callback and no borrow survives the completed call.
    #[allow(clippy::too_many_arguments)]
    pub fn lun_occult_when_loc_with_position_context(
        &self,
        jd: f64,
        body: i32,
        site: [f64; 3],
        flags: i32,
        backwards: i32,
        resolve: &PositionResolver<'_>,
    ) -> Result<NumericOccultLocalResult, Error> {
        self.numeric_occult_local(jd, body, site, flags, backwards, Some(resolve))
    }

    /// Select provider family without consuming an unused native epoch/site.
    fn numeric_occult_local(
        &self,
        jd: f64,
        body: i32,
        site: [f64; 3],
        flags: i32,
        backwards: i32,
        resolve: Option<&PositionResolver<'_>>,
    ) -> Result<NumericOccultLocalResult, Error> {
        let body = numeric_occultation_target(body)?;
        let backward = backwards & 1 != 0;
        let one_try = backwards & ECL_ONE_TRY != 0;
        let flags = event::strip_one_try_bit(flags);
        if self.ctx.mode() != ProviderMode::NativeJpl {
            self.with_serving_stack(|stack| {
                let dt = |ut| crate::time::deltat_with_config(ut, stack.delta_cfg);
                let (selection, _) = crate::houses::sidereal_selection_for_context(&self.ctx);
                let stars = crate::stars::SealedStarServing::new(
                    stack.serving.provider(),
                    stack.nutation_at,
                    &dt,
                    selection,
                )
                .with_default_timescale();
                let bound = event::SealedOccultLocalServing::try_new_planet(
                    stack.serving,
                    &stars,
                    &dt,
                    stack.nutation_at,
                    flags,
                    body,
                    site[0],
                    site[1],
                    site[2],
                )?;
                event::lun_occult_when_loc_sealed(&bound, jd, backward, one_try)
            })
        } else {
            event::lun_occult_when_loc_source_via(
                &NativeLocal {
                    body,
                    site,
                    flags,
                    geometry: NativeNumericOccultWhere::new(self, resolve),
                },
                jd,
                body,
                backward,
                one_try,
            )
        }
    }
}

/// Read-only local acquisition owner. No position result is cached or shared
/// between scalar and independent Topos samples.
struct NativeLocal<'a> {
    body: i32,
    site: [f64; 3],
    flags: i32,
    geometry: NativeNumericOccultWhere<'a>,
}
impl NativeLocal<'_> {
    /// Fresh rise followed by set, over the retained main and temporary owned
    /// rise observer. A nonnegative status contributes its first date.
    fn horizon(
        &self,
        jd: f64,
        body: i32,
        bottom: bool,
    ) -> Result<(Option<f64>, Option<f64>), Error> {
        let read = |flag| {
            self.geometry
                .rise_at(
                    jd,
                    body,
                    flag | if bottom { BIT_DISC_BOTTOM } else { 0 },
                    self.site,
                    event::ecl_eph_flags(self.flags),
                )
                .map(|(rc, t)| if rc >= 0 { Some(t[0]) } else { None })
        };
        Ok((read(CALC_RISE)?, read(CALC_SET)?))
    }
}
impl OccultLocalEphemeris for NativeLocal<'_> {
    fn body_lonlat_deg(&self, jd: f64) -> Result<(f64, f64), Error> {
        let p = self
            .geometry
            .position_at(jd, self.body, event::ecl_eph_flags(self.flags))?;
        Ok((p[0], p[1]))
    }
    fn moon_lonlat_deg(&self, jd: f64) -> Result<(f64, f64), Error> {
        let p = self
            .geometry
            .position_at(jd, MOON, event::ecl_eph_flags(self.flags))?;
        Ok((p[0], p[1]))
    }
    fn topo_radii_sep_deg(&self, jd: f64, inner: bool) -> Result<(f64, f64, f64), Error> {
        let target = self
            .geometry
            .target_at(jd, self.body, self.flags, self.site)?;
        let (_, moon) = self.geometry.pair_at(jd, self.site)?;
        let radius = event::occ_body_radius_au(Some(self.body));
        let rb = if radius != 0.0 {
            checked_asin(radius / target[2])?
        } else {
            0.0
        };
        let rm = checked_asin(
            (if inner {
                event::ECL_RMOON_INNER_AU
            } else {
                event::ECL_RMOON_AU
            }) / moon[2],
        )?;
        let separation =
            event::numeric_occult_local_separation(target[0], target[1], moon[0], moon[1])?;
        Ok((rb, rm, separation))
    }
    fn body_alt_deg(&self, jd: f64) -> Result<f64, Error> {
        Ok(self
            .geometry
            .local_at(jd, self.body, self.flags, self.site)?
            .attr()[6])
    }
    fn local_how_attr(&self, jd: f64) -> Result<[f64; 20], Error> {
        Ok(self
            .geometry
            .local_at(jd, self.body, self.flags, self.site)?
            .attr())
    }
    fn where_core(&self, jd: f64) -> Result<WhereGeometry, Error> {
        self.geometry.global_at(jd, self.body, self.flags)
    }
    fn body_rise_set(&self, jd: f64) -> Result<(Option<f64>, Option<f64>), Error> {
        self.horizon(jd, self.body, true)
    }
    fn sun_rise_set(&self, jd: f64) -> Result<(Option<f64>, Option<f64>), Error> {
        self.horizon(jd, SUN, false)
    }
}

/// Consumed radius inverse sine in degrees, preserving its actual operand.
fn checked_asin(value: f64) -> Result<f64, Error> {
    if value > 1.0 || value < -1.0 {
        return Err(Error::MathDomain { value });
    }
    Ok(value.asin().to_degrees())
}

#[cfg(test)]
#[path = "numeric_occult_local/tests.rs"]
mod tests;
