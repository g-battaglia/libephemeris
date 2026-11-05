//! Minor-body event serving and the retained channel-less Keplerian route.
//!
//! Channel selection, remapped identifiers, tropical-frame admission and the
//! family-specific apparent-place policy stay with this provider adapter.

use super::{
    Error, Result, SealedMinorEventServing, SealedTransitServing, TransitSolution,
    calculate_transit, cotrans,
};

/// True for the out-of-map minor-body ids served on this path (after
/// the release `_remap_ast_offset` rule): the 37 release
/// `MINOR_BODY_ELEMENTS` keys — dedicated ids 15-20 plus the 31
/// `AST_OFFSET + N` numbers.
pub fn is_minor_event_body(body_id: i32) -> bool {
    let remapped = crate::positions::remap_ast_offset(body_id);
    crate::minor_bodies::minor_body_elements(remapped).is_some()
}

impl<'a, P> SealedMinorEventServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Checked constructor: refuses the degenerate `EARTH`/`ECL_NUT`
    /// ids (`InvalidBody`, the release `undefined` arm), non-table ids
    /// (`UnknownBody`, the release `illegal planet number` probe arm),
    /// and flag selections outside the verified tropical basis. No I/O,
    /// no globals.
    pub fn try_new(
        serving: &'a crate::positions::SealedServing<'a, P>,
        deltat_at: &'a dyn Fn(f64) -> Result<f64>,
        body_id: i32,
        calc_flags: i32,
    ) -> Result<Self> {
        use crate::constants::{
            EARTH, ECL_NUT, FLG_BARYCTR, FLG_HELCTR, FLG_J2000, FLG_NONUT, FLG_RADIANS,
            FLG_SIDEREAL, FLG_XYZ,
        };
        if body_id == EARTH || body_id == ECL_NUT {
            return Err(Error::InvalidBody {
                body_id,
                operation: "minor_rise_set_ut",
            });
        }
        let minor_id = crate::positions::remap_ast_offset(body_id);
        if crate::minor_bodies::minor_body_elements(minor_id).is_none() {
            return Err(Error::UnknownBody { body_id });
        }
        if (calc_flags & (FLG_SIDEREAL | FLG_J2000 | FLG_NONUT | FLG_XYZ | FLG_RADIANS)) != 0
            || (calc_flags & (FLG_HELCTR | FLG_BARYCTR)) != 0
        {
            return Err(Error::UnsupportedFlagCombination {
                entry: "minor_event_serving",
                flags: calc_flags,
            });
        }
        let dt_lo = deltat_at(crate::time::SIDEREAL_MODERN_START_JD)?;
        if !dt_lo.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_lo });
        }
        let dt_hi = deltat_at(crate::time::SIDEREAL_MODERN_END_JD)?;
        if !dt_hi.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_hi });
        }
        let (offset_lo_deg, offset_hi_deg) =
            crate::time::sidereal_continuity_offsets(dt_lo, dt_hi)?;
        let leb_channeled = serving.provider().has_body(minor_id);
        Ok(Self {
            serving,
            deltat_at,
            body_id,
            minor_id,
            calc_flags,
            leb_channeled,
            offset_lo_deg,
            offset_hi_deg,
        })
    }

    /// Requested body id (an `AST_OFFSET + N` alias keeps its requested
    /// form; the remapped id drives every read).
    pub fn body_id(&self) -> i32 {
        self.body_id
    }

    /// Remapped table id behind an `AST_OFFSET + N` alias (identity for
    /// dedicated ids and direct numbers).
    pub fn minor_id(&self) -> i32 {
        self.minor_id
    }

    /// Whether the bound composite provider carries an LEB channel for
    /// the remapped id (the 36 channeled table ids take the shared
    /// serving; the channel-less remainder takes the
    /// Keplerian-fallback composition below).
    pub fn is_leb_channeled(&self) -> bool {
        self.leb_channeled
    }

    /// Bound sealed serving pair.
    pub fn serving(&self) -> &crate::positions::SealedServing<'a, P> {
        self.serving
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

    /// Map a sealed provider miss onto the public range contract
    /// a missing channel inside asset
    /// coverage stays the typed miss, past coverage it names the
    /// requested body with the asset window.
    fn escalate(&self, jd_tt: f64, body_id: i32, err: Error) -> Error {
        match err {
            Error::LebBodyNotFound { .. } | Error::LebOutOfRange { .. } => {
                crate::positions::escalate_sealed_provider_miss(
                    crate::context::ProviderMode::SealedLeb,
                    self.serving.provider(),
                    body_id,
                    jd_tt,
                    err,
                )
            }
            other => other,
        }
    }

    /// Shared transit serving over the bound pair for the
    /// LEB-channeled path.
    pub(in crate::events) fn inner_transit(&self) -> Result<SealedTransitServing<'a, P>> {
        SealedTransitServing::try_new(self.serving, self.deltat_at, self.minor_id, self.calc_flags)
    }

    /// Geocentric ecliptic `(lon, lat, dist)` of the bound minor body
    /// at a TT epoch through the same `_keplerian_position_at`
    /// non-heliocentric arm as the release: heliocentric true-ecliptic-of-date
    /// words differenced with the
    /// sealed LEB Earth heliocentric place in the same frame (no
    /// aberration, no topocentric parallax — the fallback ignores
    /// `FLG_TOPOCTR`). The Earth vector is retarded by its own
    /// light-time. Only table ids without an LEB
    /// channel reach this composition; everything else runs the shared
    /// Pipeline-A serving.
    pub fn kepler_ecliptic_at_tt(&self, jd_tt: f64) -> Result<(f64, f64, f64)> {
        use crate::constants::{EARTH, SUN};
        use crate::positions::C_LIGHT_AU_DAY;
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let (dpsi_rad, deps_rad) = (self.serving.nutation_provider())(jd_tt)?;
        if !dpsi_rad.is_finite() || !deps_rad.is_finite() {
            return Err(Error::NonFiniteArgument { value: dpsi_rad });
        }
        // Heliocentric true-ecliptic-of-date words.
        let (lon_h, lat_h, r_h) = crate::minor_bodies::calc_minor_body_heliocentric(
            self.minor_id,
            jd_tt,
            false,
            None,
            dpsi_rad,
        )
        .map_err(|e| self.escalate(jd_tt, self.minor_id, e))?;
        // Sealed LEB Earth heliocentric place in the true ecliptic of
        // date: ICRS SSB states differenced, then the same
        //  true-of-date rotation the Pipeline-A shaping uses.
        let provider = self.serving.provider();
        let sun = provider
            .geometric_state(SUN, jd_tt)
            .map_err(|e| self.escalate(jd_tt, SUN, e))?;
        let sp = sun.pos();
        let mut earth_pos = provider
            .geometric_state(EARTH, jd_tt)
            .map_err(|e| self.escalate(jd_tt, EARTH, e))?
            .pos();
        for _ in 0..3 {
            let dx = earth_pos[0] - sp[0];
            let dy = earth_pos[1] - sp[1];
            let dz = earth_pos[2] - sp[2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            if dist == 0.0 {
                break;
            }
            let lt_days = dist / C_LIGHT_AU_DAY;
            earth_pos = provider
                .geometric_state(EARTH, jd_tt - lt_days)
                .map_err(|e| self.escalate(jd_tt, EARTH, e))?
                .pos();
        }
        let ep = earth_pos;
        let helio_icrs = [ep[0] - sp[0], ep[1] - sp[1], ep[2] - sp[2]];
        let frames = crate::positions::frame_matrices_for_date(jd_tt, dpsi_rad, deps_rad)?;
        let to_ecl = |v: [f64; 3]| -> Result<[f64; 3]> {
            let eq = crate::positions::mat3_vec3(frames.pn(), v)?;
            crate::positions::rotate_equatorial_to_ecliptic(eq, frames.eps_true_rad())
        };
        let (ex, ey, ez) = {
            let e = to_ecl(helio_icrs)?;
            (e[0], e[1], e[2])
        };
        // Heliocentric spherical to Cartesian ( association),
        // geocentric differencing, back to spherical.
        let lon_r = lon_h.to_radians();
        let lat_r = lat_h.to_radians();
        let cl = lat_r.cos();
        let xh = r_h * cl * lon_r.cos();
        let yh = r_h * cl * lon_r.sin();
        let zh = r_h * lat_r.sin();
        let xg = xh - ex;
        let yg = yh - ey;
        let zg = zh - ez;
        let r_geo = (xg * xg + yg * yg + zg * zg).sqrt();
        let lon = yg.atan2(xg).to_degrees().rem_euclid(360.0);
        let lat = if r_geo > 0.0 {
            (zg / r_geo).asin().to_degrees()
        } else {
            0.0
        };
        Ok((lon, lat, r_geo))
    }

    /// Apparent `(RA hours, Dec deg)` at a TT epoch over the fallback
    /// words.
    pub fn kepler_ra_dec_at_tt(&self, jd_tt: f64) -> Result<(f64, f64)> {
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let (lon, lat, _) = self.kepler_ecliptic_at_tt(jd_tt)?;
        let (_, deps_rad) = (self.serving.nutation_provider())(jd_tt)?;
        let eps_true_deg =
            crate::positions::vondrak_mean_obliquity_deg(jd_tt)? + deps_rad.to_degrees();
        let (ra_deg, dec_deg, _) = cotrans(lon, lat, 1.0, -eps_true_deg);
        Ok((ra_deg / 15.0, dec_deg))
    }

    /// Apparent `(RA hours, Dec deg)` at a UT epoch over the fallback
    /// words.
    pub fn kepler_ra_dec_at_ut(&self, jd_ut: f64) -> Result<(f64, f64)> {
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        self.kepler_ra_dec_at_tt(jd_tt)
    }

    /// Greenwich apparent sidereal time in hours `[0, 24)` at a UT epoch
    /// (same acquisition as [`SealedTransitServing::gast_hours_at_ut`]).
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
        let (dpsi_rad, deps_rad) = (self.serving.nutation_provider())(jd_tt)?;
        let dpsi_deg = dpsi_rad.to_degrees();
        let eps_mean_deg = crate::positions::vondrak_mean_obliquity_deg(jd_tt)?;
        let eps_true_deg = eps_mean_deg + deps_rad.to_degrees();
        let sdt = SiderealDeltaT::try_new(dt_days, self.offset_lo_deg, self.offset_hi_deg)?;
        greenwich_apparent_sidereal_time_hours(jd_ut, eps_true_deg, dpsi_deg, &sdt)
    }

    /// Local apparent sidereal time in hours `[0, 24)` at a UT epoch for
    /// an observer longitude (same acquisition as
    /// [`SealedTransitServing::last_hours_at_ut`]).
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
        let (dpsi_rad, deps_rad) = (self.serving.nutation_provider())(jd_tt)?;
        let dpsi_deg = dpsi_rad.to_degrees();
        let eps_mean_deg = crate::positions::vondrak_mean_obliquity_deg(jd_tt)?;
        let eps_true_deg = eps_mean_deg + deps_rad.to_degrees();
        let sdt = SiderealDeltaT::try_new(dt_days, self.offset_lo_deg, self.offset_hi_deg)?;
        local_apparent_sidereal_time_hours(jd_ut, lon_deg_east, eps_true_deg, dpsi_deg, &sdt)
    }

    /// Meridian transit over the fallback words (same
    /// `_calculate_transit` acquisition as the release; the
    /// [`calculate_transit`] transit core). `lower_transit` selects the
    /// anti-culmination (`CALC_ITRANSIT` 180-degree shift).
    pub fn kepler_transit_at_ut(
        &self,
        jd_start: f64,
        lon_deg_east: f64,
        lower_transit: bool,
    ) -> Result<TransitSolution> {
        if !jd_start.is_finite() {
            return Err(Error::NonFiniteSearchValue { x: jd_start });
        }
        let mut get_ra_dec = |jd: f64| self.kepler_ra_dec_at_ut(jd);
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
