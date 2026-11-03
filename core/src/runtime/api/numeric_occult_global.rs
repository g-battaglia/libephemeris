//! Global numeric lunar-occultation search over call-owned scalar acquisitions.
use super::solar_circumstances::{PositionResolver, SelectedMain, solar_position};
use super::*;
use crate::events::{self as event, OccultGlobalEphemeris, WhereGeometry};

impl Session {
    /// Search globally for a lunar occultation of a numeric body from UT1 JD.
    /// Negative ids fold to Sun; Moon-self refuses first. `backwards` bit0
    /// selects the earlier event and ECL_ONE_TRY examines one conjunction.
    /// Ephemeris flag bits alone select scalar target/Moon positions; the type
    /// mask is normalized before images or epochs are acquired. The result is
    /// `(character, times10)`: maximum, RA noon-transit, outer begin/end, total
    /// begin/end, center-line begin/end, and two reserved zeros. Missing phases
    /// are zero. A one-try miss retains its continuation JD in slot zero.
    ///
    /// Conjunctions use 300 mean-rate steps; maximum and contact solves use 60
    /// steps with 1e-7/1e-8-day widths. Provider failures propagate except a
    /// first-target range boundary, which ends with the owned no-event Error.
    /// Model/scalar and independent GAST clocks remain separately consumed.
    /// The immutable selected main lives through all reads; observer state is
    /// unchanged. Sealed/LEB modes retain checked serving and flag/domain gates.
    /// Numeric positive-distance discs form the numerical domain. Integer
    /// catalog and string-star targets are excluded. Undefined numeric
    /// geometry returns an error instead of a no-event result.
    pub fn lun_occult_when_glob(
        &self,
        jd_ut1: f64,
        body: i32,
        flags: i32,
        mask: i32,
        backwards: i32,
    ) -> Result<(i32, [f64; 10]), Error> {
        self.numeric_occult_global(jd_ut1, body, flags, mask, backwards, None)
    }

    /// Same global search with a synchronous configuration-only main resolver.
    /// Supplied registry entries keep priority. Resolution occurs at the first
    /// position acquisition and returns an owned context, never position words.
    /// Non-Native modes do not invoke it; output/time/error contracts are those
    /// of [`Self::lun_occult_when_glob`].
    #[allow(clippy::too_many_arguments)]
    pub fn lun_occult_when_glob_with_position_context(
        &self,
        jd_ut1: f64,
        body: i32,
        flags: i32,
        mask: i32,
        backwards: i32,
        resolve: &PositionResolver<'_>,
    ) -> Result<(i32, [f64; 10]), Error> {
        self.numeric_occult_global(jd_ut1, body, flags, mask, backwards, Some(resolve))
    }

    /// Own one main selection while preserving target/filter acquisition order.
    #[allow(clippy::too_many_arguments)]
    fn numeric_occult_global(
        &self,
        jd: f64,
        body: i32,
        flags: i32,
        mask: i32,
        backwards: i32,
        resolve: Option<&PositionResolver<'_>>,
    ) -> Result<(i32, [f64; 10]), Error> {
        let body = super::numeric_occultation::numeric_occultation_target(body)?;
        event::normalize_occultation_filter(mask, body == SUN, &body.to_string())?;
        let backward = backwards & 1 != 0;
        let one_try = backwards & ECL_ONE_TRY != 0;
        let flags = event::ecl_eph_flags(flags);
        if self.ctx.mode() != ProviderMode::NativeJpl {
            return self.with_serving_stack(|stack| {
                let dt = |ut| crate::time::deltat_with_config(ut, stack.delta_cfg);
                let (selection, _) = crate::houses::sidereal_selection_for_context(&self.ctx);
                let stars = crate::stars::SealedStarServing::new(
                    stack.serving.provider(),
                    stack.nutation_at,
                    &dt,
                    selection,
                )
                .with_default_timescale();
                let bound = event::SealedOccultGlobalServing::try_new_planet(
                    stack.serving,
                    &stars,
                    &dt,
                    flags,
                    body,
                )?;
                event::lun_occult_when_glob_sealed(&bound, jd, mask, backward, one_try)
            });
        }
        event::lun_occult_when_glob_source_via(
            &NativeNumericOccultGlobal {
                session: self,
                body,
                flags,
                resolve,
                selected: SelectedMain::new(),
            },
            jd,
            body,
            mask,
            backward,
            one_try,
        )
    }
}

/// Retained Session/main pair for all target, Moon, geometry and RA samples.
struct NativeNumericOccultGlobal<'a> {
    session: &'a Session,
    body: i32,
    flags: i32,
    resolve: Option<&'a PositionResolver<'a>>,
    selected: SelectedMain,
}
impl NativeNumericOccultGlobal<'_> {
    /// Public calc_ut scalar words at caller UT1 with exact acquisition flags.
    fn position(&self, jd: f64, body: i32, shape: i32) -> Result<[f64; 6], Error> {
        solar_position(
            self.session,
            jd,
            body,
            self.flags | shape,
            self.resolve,
            &self.selected,
        )
    }
}
impl OccultGlobalEphemeris for NativeNumericOccultGlobal<'_> {
    fn body_lonlat_deg(&self, jd: f64) -> Result<(f64, f64), Error> {
        let p = self.position(jd, self.body, 0)?;
        Ok((p[0], p[1]))
    }
    fn moon_lonlat_deg(&self, jd: f64) -> Result<(f64, f64), Error> {
        let p = self.position(jd, MOON, 0)?;
        Ok((p[0], p[1]))
    }
    fn sep_minus_radii_deg(&self, jd: f64) -> Result<f64, Error> {
        let body = self.position(jd, self.body, FLG_EQUATORIAL | FLG_XYZ)?;
        let moon = self.position(jd, MOON, FLG_EQUATORIAL | FLG_XYZ)?;
        event::numeric_occult_separation_source(
            [body[0], body[1], body[2]],
            [moon[0], moon[1], moon[2]],
            event::occ_body_radius_au(Some(self.body)),
        )
    }
    fn where_core(&self, jd: f64) -> Result<WhereGeometry, Error> {
        let moon = self.position(jd, MOON, FLG_EQUATORIAL | FLG_XYZ)?;
        let body = self.position(jd, self.body, FLG_EQUATORIAL | FLG_XYZ)?;
        event::eclipse_where_core_lazy(
            (moon[0], moon[1], moon[2]),
            (body[0], body[1], body[2]),
            &mut || self.session.sidtime(jd).map(|h| h * 15.0),
            event::occ_body_radius_au(Some(self.body)),
        )
    }
    fn ra_gap_deg(&self, jd: f64) -> Result<f64, Error> {
        let body = self.position(jd, self.body, FLG_EQUATORIAL)?;
        let moon = self.position(jd, MOON, FLG_EQUATORIAL)?;
        Ok(event::wrap_ra_gap_deg(body[0] - moon[0]))
    }
}

#[cfg(test)]
#[path = "numeric_occult_global/tests.rs"]
mod tests;
