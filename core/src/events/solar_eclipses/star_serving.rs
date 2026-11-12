//! Resolved-star state and sidereal adapters for event searches.
//!
//! The bound catalogue serving retains star resolution and place acquisition;
//! this adapter supplies the existing UT/TT, horizon and transit inputs.

use super::{Error, Result, SealedStarEventServing, TransitSolution, calculate_transit};

impl<'a, P> SealedStarEventServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Checked constructor: resolves the star name through the v1
    /// (`fixstar_ut`) grammar and the sidereal-continuity offsets from
    /// the Delta-T provider's own boundary words. No I/O, no globals.
    pub fn try_new(
        star_serving: &'a crate::stars::SealedStarServing<'a, P>,
        star_name: &str,
        deltat_at: &'a dyn Fn(f64) -> Result<f64>,
        nutation_at: &'a dyn Fn(f64) -> Result<(f64, f64)>,
    ) -> Result<Self> {
        use crate::time::{
            SIDEREAL_MODERN_END_JD, SIDEREAL_MODERN_START_JD, sidereal_continuity_offsets,
        };
        let star = match crate::stars::resolve_star_ref(star_name) {
            Ok(matched) => matched,
            Err(crate::stars::StarResolveError::Miss(failure)) => {
                return Err(failure.to_error(star_name));
            }
            Err(crate::stars::StarResolveError::Refused(err)) => return Err(err),
        };
        let dt_lo = deltat_at(SIDEREAL_MODERN_START_JD)?;
        if !dt_lo.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_lo });
        }
        let dt_hi = deltat_at(SIDEREAL_MODERN_END_JD)?;
        if !dt_hi.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_hi });
        }
        let (offset_lo_deg, offset_hi_deg) = sidereal_continuity_offsets(dt_lo, dt_hi)?;
        Ok(Self {
            star_serving,
            star,
            deltat_at,
            nutation_at,
            offset_lo_deg,
            offset_hi_deg,
        })
    }

    /// Bound resolved star (id plus canonical `"Name,nomenclature"`).
    pub fn star(&self) -> &crate::stars::StarMatch {
        &self.star
    }

    /// Bound sealed nutation pair (radians at a TT epoch).
    pub fn nutation_provider(&self) -> &dyn Fn(f64) -> Result<(f64, f64)> {
        self.nutation_at
    }

    /// TT epoch for a UT epoch through the bound Delta-T provider.
    pub fn jd_tt_for_ut(&self, jd_ut: f64) -> Result<f64> {
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let dt = (self.deltat_at)(jd_ut)?;
        if !dt.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt });
        }
        Ok(jd_ut + dt)
    }

    /// Apparent geocentric ecliptic triple `(lon deg, lat deg, dist AU)`
    /// at a UT epoch.
    pub fn star_ecliptic_at_ut(&self, jd_ut: f64) -> Result<[f64; 3]> {
        use crate::constants::FLG_SPEED;
        let (place, _, _) = self
            .star_serving
            .fixstar_ut(&self.star, jd_ut, FLG_SPEED, None)?;
        Ok([place[0], place[1], place[2]])
    }

    /// Apparent `(RA hours, Dec deg)` at a UT epoch.
    pub fn star_ra_dec_at_ut(&self, jd_ut: f64) -> Result<(f64, f64)> {
        use crate::constants::{FLG_EQUATORIAL, FLG_SPEED};
        let (place, _, _) =
            self.star_serving
                .fixstar_ut(&self.star, jd_ut, FLG_EQUATORIAL | FLG_SPEED, None)?;
        Ok((place[0] / 15.0, place[1]))
    }

    /// Greenwich apparent sidereal time in hours `[0, 24)` at a UT epoch
    /// the
    /// native sidereal-time composition over the sealed nutation pair at the TT
    /// epoch and the Vondrak mean obliquity the frame pipeline uses.
    pub fn gast_hours_at_ut(&self, jd_ut: f64) -> Result<f64> {
        use crate::time::{SiderealDeltaT, greenwich_apparent_sidereal_time_hours};
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let dt_days = (self.deltat_at)(jd_ut)?;
        if !dt_days.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_days });
        }
        let jd_tt = jd_ut + dt_days;
        let (dpsi_rad, deps_rad) = (self.nutation_at)(jd_tt)?;
        let dpsi_deg = dpsi_rad.to_degrees();
        let eps_mean_deg = crate::positions::vondrak_mean_obliquity_deg(jd_tt)?;
        let eps_true_deg = eps_mean_deg + deps_rad.to_degrees();
        let sdt = SiderealDeltaT::try_new(dt_days, self.offset_lo_deg, self.offset_hi_deg)?;
        greenwich_apparent_sidereal_time_hours(jd_ut, eps_true_deg, dpsi_deg, &sdt)
    }

    /// Local apparent sidereal time in hours `[0, 24)` at a UT epoch for
    /// an observer longitude (same `_sidtime_internal` composition as the
    /// release): the `azalt` LST acquisition the sealed
    /// rise/set altitude engine uses.
    pub fn last_hours_at_ut(&self, jd_ut: f64, lon_deg_east: f64) -> Result<f64> {
        use crate::time::{SiderealDeltaT, local_apparent_sidereal_time_hours};
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        if !lon_deg_east.is_finite() {
            return Err(Error::InvalidObserver {
                field: "lon_deg",
                value: lon_deg_east,
            });
        }
        let dt_days = (self.deltat_at)(jd_ut)?;
        if !dt_days.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_days });
        }
        let jd_tt = jd_ut + dt_days;
        let (dpsi_rad, deps_rad) = (self.nutation_at)(jd_tt)?;
        let dpsi_deg = dpsi_rad.to_degrees();
        let eps_mean_deg = crate::positions::vondrak_mean_obliquity_deg(jd_tt)?;
        let eps_true_deg = eps_mean_deg + deps_rad.to_degrees();
        let sdt = SiderealDeltaT::try_new(dt_days, self.offset_lo_deg, self.offset_hi_deg)?;
        local_apparent_sidereal_time_hours(jd_ut, lon_deg_east, eps_true_deg, dpsi_deg, &sdt)
    }

    /// Meridian transit over the sealed star serving (same
    /// `_calculate_transit` acquisition as the release; the
    /// [`calculate_transit`] transit core). `lower_transit` selects the
    /// anti-culmination (`CALC_ITRANSIT` 180-degree shift). The Hindu bit
    /// is ignored for transits, like the release.
    pub fn star_transit_at_ut(
        &self,
        jd_start: f64,
        lon_deg_east: f64,
        lower_transit: bool,
    ) -> Result<TransitSolution> {
        if !jd_start.is_finite() {
            return Err(Error::NonFiniteSearchValue { x: jd_start });
        }
        let mut get_ra_dec = |jd: f64| self.star_ra_dec_at_ut(jd);
        let mut gast_hours = |jd: f64| self.gast_hours_at_ut(jd);
        calculate_transit(
            jd_start,
            lon_deg_east,
            lower_transit,
            &mut get_ra_dec,
            &mut gast_hours,
        )
    }
}
