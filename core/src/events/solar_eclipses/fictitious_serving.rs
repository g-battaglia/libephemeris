//! Analytical fictitious-body positions and sidereal inputs for event searches.
//!
//! The serving object borrows its state, nutation and Delta-T providers. This
//! module retains the family-specific reductions and propagates their errors;
//! it does not acquire a provider or change the caller's configuration.

use super::{
    add_of_date_nutation_lon, bradley_aberration_xyz, cotrans, ecl_xyz_to_spherical,
    fictitious_family, first_order_aberration_xyz, rotate_vector_to_of_date,
};
use crate::errors::{Error, Result};
use crate::events::{FictFamily, SealedFictitiousEventServing, TransitSolution, calculate_transit};

impl<'a, P> SealedFictitiousEventServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Bind an analytical fictitious-body model to borrowed state and nutation providers.
    ///
    /// `deltat_at` takes a UT Julian date and returns TT minus UT in days. It is
    /// called at both modern sidereal-time boundaries during construction; its
    /// errors and nonfinite results are propagated. The providers remain
    /// borrowed for the lifetime of this serving object.
    ///
    /// Returns `InvalidBody` for Earth or the ecliptic/nutation pseudo-body,
    /// `UnknownBody` for ids without a supported hypothetical model, and
    /// `UnsupportedFlagCombination` for sidereal, J2000, non-nutated, Cartesian,
    /// radian, heliocentric or barycentric selections. The supported analytical
    /// models are hypothetical; this interface does not establish physical
    /// accuracy for their trajectories.
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
                operation: "fictitious_rise_set_ut",
            });
        }
        if !crate::hypothetical::is_hypothetical_body(body_id)
            || fictitious_family(body_id).is_none()
        {
            return Err(Error::UnknownBody { body_id });
        }
        if (calc_flags & (FLG_SIDEREAL | FLG_J2000 | FLG_NONUT | FLG_XYZ | FLG_RADIANS)) != 0
            || (calc_flags & (FLG_HELCTR | FLG_BARYCTR)) != 0
        {
            return Err(Error::UnsupportedFlagCombination {
                entry: "fictitious_event_serving",
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
        Ok(Self {
            serving,
            deltat_at,
            body_id,
            calc_flags,
            offset_lo_deg,
            offset_hi_deg,
        })
    }

    /// Requested fictitious-body identifier.
    pub fn body_id(&self) -> i32 {
        self.body_id
    }

    /// Caller flags retained at construction after the tropical-basis checks.
    /// The family-specific apparent-place reductions are applied independently
    /// of additional flags that the constructor accepts.
    pub fn calc_flags(&self) -> i32 {
        self.calc_flags
    }

    /// Borrowed sealed state and nutation providers used by this object.
    pub fn serving(&self) -> &crate::positions::SealedServing<'a, P> {
        self.serving
    }

    /// Convert a UT Julian date to TT using the bound Delta-T callback in days.
    ///
    /// Returns an error for a nonfinite input or callback result and propagates
    /// callback errors. The result is the binary64 sum of the input and Delta-T.
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

    /// Preserve a missing-channel error inside the provider's coverage and
    /// identify the requested body and coverage window for out-of-range errors.
    /// Other provider errors pass through unchanged.
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

    /// Retarded heliocentric Earth position in the J2000 ecliptic frame, in AU.
    ///
    /// The Sun stays at the requested TT Julian date. Earth is retarded by its
    /// own light-time for at most ten iterations, stopping when successive
    /// delays differ by less than `1e-14` day. The ICRS difference is then rotated
    /// to the J2000 ecliptic frame. Provider and nonfinite-input errors propagate.
    fn earth_helio_ecl_j2000(&self, jd_tt: f64) -> Result<[f64; 3]> {
        self.earth_helio_ecl_j2000_with_origin(jd_tt)
            .map(|observation| observation.0)
    }

    /// Retain the reception Sun in ICRS AU without another provider read.
    fn earth_helio_ecl_j2000_with_origin(&self, jd_tt: f64) -> Result<([f64; 3], [f64; 3])> {
        use crate::constants::{EARTH, SUN};
        use crate::positions::C_LIGHT_AU_DAY;
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let provider = self.serving.provider();
        let sun = provider
            .geometric_state(SUN, jd_tt)
            .map_err(|e| self.escalate(jd_tt, SUN, e))?
            .pos();
        let mut earth = provider
            .geometric_state(EARTH, jd_tt)
            .map_err(|e| self.escalate(jd_tt, EARTH, e))?
            .pos();
        let mut dx = earth[0] - sun[0];
        let mut dy = earth[1] - sun[1];
        let mut dz = earth[2] - sun[2];
        let mut lt_days = (dx * dx + dy * dy + dz * dz).sqrt() / C_LIGHT_AU_DAY;
        for _ in 0..10 {
            earth = provider
                .geometric_state(EARTH, jd_tt - lt_days)
                .map_err(|e| self.escalate(jd_tt, EARTH, e))?
                .pos();
            dx = earth[0] - sun[0];
            dy = earth[1] - sun[1];
            dz = earth[2] - sun[2];
            let lt_next = (dx * dx + dy * dy + dz * dz).sqrt() / C_LIGHT_AU_DAY;
            let converged = (lt_next - lt_days).abs() < 1e-14;
            lt_days = lt_next;
            if converged {
                break;
            }
        }
        let rot = crate::positions::icrs_to_ecl_j2000_matrix();
        Ok((crate::positions::mat3_vec3(rot, [dx, dy, dz])?, sun))
    }

    /// Earth's SSB velocity at the requested TT Julian date, rotated to the
    /// J2000 ecliptic frame, in AU/day. This supplies the annual-aberration term.
    fn earth_bary_vel_ecl_j2000(&self, jd_tt: f64) -> Result<[f64; 3]> {
        use crate::constants::EARTH;
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let vel = self
            .serving
            .provider()
            .geometric_state(EARTH, jd_tt)
            .map_err(|e| self.escalate(jd_tt, EARTH, e))?
            .vel();
        let rot = crate::positions::icrs_to_ecl_j2000_matrix();
        crate::positions::mat3_vec3(rot, vel)
    }

    /// Analytical heliocentric position in the J2000 ecliptic frame, in AU.
    ///
    /// Uranian, Transpluto and predicted-planet models supply J2000 coordinates.
    /// Vulcan and Proserpina supply ecliptic-of-date coordinates, which are
    /// precessed to J2000. Model and nonfinite-input errors propagate.
    fn helio_j2000_xyz(&self, jd_tt: f64) -> Result<[f64; 3]> {
        use crate::constants::{PROSERPINA, VULCAN};
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let state = crate::hypothetical::calc_hypothetical_position(self.body_id, jd_tt)?;
        let (lon_j, lat_j) = if self.body_id == VULCAN || self.body_id == PROSERPINA {
            crate::hypothetical::precess_ecliptic(
                state.lon_deg,
                state.lat_deg,
                jd_tt,
                crate::hypothetical::J2000_JD,
            )?
        } else {
            (state.lon_deg, state.lat_deg)
        };
        let lon_r = lon_j.to_radians();
        let lat_r = lat_j.to_radians();
        let cl = lat_r.cos();
        Ok([
            state.dist_au * cl * lon_r.cos(),
            state.dist_au * cl * lon_r.sin(),
            state.dist_au * lat_r.sin(),
        ])
    }

    /// Apparent geocentric `(longitude, latitude, distance)` at a TT Julian date.
    ///
    /// Angles are in degrees and distance is in AU. The ecliptic-of-date result
    /// includes each family's prescribed light-time, aberration, precession
    /// and longitude nutation; the predicted-heliocentric family also applies
    /// gravitational deflection. The symbolic lunar models have no precession
    /// step. This method does not apply topocentric parallax.
    ///
    /// Nonfinite epochs, unavailable provider data, unsupported models and
    /// errors from the family-specific reductions are returned to the caller.
    pub fn fict_geo_ecliptic_at_tt(&self, jd_tt: f64) -> Result<(f64, f64, f64)> {
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        match fictitious_family(self.body_id) {
            Some(FictFamily::Uranian) => self.uranian_geo_of_date(jd_tt),
            Some(FictFamily::Transpluto) => self.transpluto_geo_of_date(jd_tt),
            Some(FictFamily::PredictedHelio) => self.predicted_geo_of_date(jd_tt),
            Some(FictFamily::WhiteMoon) => self.white_moon_geo_of_date(jd_tt),
            Some(FictFamily::Waldemath) => self.waldemath_geo_of_date(jd_tt),
            None => Err(Error::UnknownBody {
                body_id: self.body_id,
            }),
        }
    }

    /// Uranian reduction: heliocentric differencing, Bradley aberration with
    /// Earth's barycentric velocity, precession from J2000 to the mean ecliptic
    /// of date, then longitude nutation. No target light-time is applied.
    fn uranian_geo_of_date(&self, jd_tt: f64) -> Result<(f64, f64, f64)> {
        let helio = self.helio_j2000_xyz(jd_tt)?;
        let earth = self.earth_helio_ecl_j2000(jd_tt)?;
        let vel = self.earth_bary_vel_ecl_j2000(jd_tt)?;
        let geo = [
            helio[0] - earth[0],
            helio[1] - earth[1],
            helio[2] - earth[2],
        ];
        let ab = bradley_aberration_xyz(geo, vel);
        let (lon_j, lat_j, dist) = ecl_xyz_to_spherical(ab);
        let (lon, lat) = crate::hypothetical::precess_ecliptic(
            lon_j,
            lat_j,
            crate::hypothetical::J2000_JD,
            jd_tt,
        )?;
        Ok((
            add_of_date_nutation_lon(self.serving, lon, jd_tt)?,
            lat,
            dist,
        ))
    }

    /// Transpluto reduction with the retarded Earth position held fixed and
    /// three target light-time iterations. First-order aberration, J2000
    /// spherical conversion, precession and longitude nutation follow.
    fn transpluto_geo_of_date(&self, jd_tt: f64) -> Result<(f64, f64, f64)> {
        use crate::positions::C_LIGHT_AU_DAY;
        let earth = self.earth_helio_ecl_j2000(jd_tt)?;
        let vel = self.earth_bary_vel_ecl_j2000(jd_tt)?;
        let (mut xg, mut yg, mut zg) = (0.0, 0.0, 0.0);
        let mut lt_days = 0.0;
        for _ in 0..3 {
            let h = self.helio_j2000_xyz(jd_tt - lt_days)?;
            xg = h[0] - earth[0];
            yg = h[1] - earth[1];
            zg = h[2] - earth[2];
            lt_days = (xg * xg + yg * yg + zg * zg).sqrt() / C_LIGHT_AU_DAY;
        }
        let mut rg = (xg * xg + yg * yg + zg * zg).sqrt();
        if rg > 0.0 {
            let ab = first_order_aberration_xyz([xg, yg, zg], vel, rg);
            xg = ab[0];
            yg = ab[1];
            zg = ab[2];
            rg = (xg * xg + yg * yg + zg * zg).sqrt();
        }
        let (lon_j, lat_j, _) = ecl_xyz_to_spherical([xg, yg, zg]);
        let (lon, lat) = crate::hypothetical::precess_ecliptic(
            lon_j,
            lat_j,
            crate::hypothetical::J2000_JD,
            jd_tt,
        )?;
        Ok((add_of_date_nutation_lon(self.serving, lon, jd_tt)?, lat, rg))
    }

    /// Predicted-heliocentric reduction: three target light-time iterations
    /// with a fixed Earth position, gravitational deflection, first-order
    /// aberration, J2000 spherical conversion, precession and longitude nutation.
    fn predicted_geo_of_date(&self, jd_tt: f64) -> Result<(f64, f64, f64)> {
        use crate::positions::C_LIGHT_AU_DAY;
        let (earth, sun_reception) = self.earth_helio_ecl_j2000_with_origin(jd_tt)?;
        let vel = self.earth_bary_vel_ecl_j2000(jd_tt)?;
        let (mut xg, mut yg, mut zg) = (0.0, 0.0, 0.0);
        let mut lt_days = 0.0;
        for _ in 0..3 {
            let h = self.helio_j2000_xyz(jd_tt - lt_days)?;
            xg = h[0] - earth[0];
            yg = h[1] - earth[1];
            zg = h[2] - earth[2];
            lt_days = (xg * xg + yg * yg + zg * zg).sqrt() / C_LIGHT_AU_DAY;
        }
        let mut rg = (xg * xg + yg * yg + zg * zg).sqrt();
        if rg > 0.0 {
            // Match the event ray's IAU 2006 ecliptic axes (including bias)
            // and its fixed reception-Sun origin at both PPN query epochs.
            let lookup = |id: i32, jd: f64| -> Result<[f64; 3]> {
                let position = self.serving.provider().geometric_state(id, jd)?.pos();
                crate::positions::deflector_in_frame(position, sun_reception, |relative| {
                    crate::positions::mat3_vec3(
                        crate::positions::icrs_to_ecl_j2000_matrix(),
                        relative,
                    )
                })
            };
            let deflected = crate::positions::apply_gravitational_deflection(
                [xg, yg, zg],
                earth,
                jd_tt,
                rg / C_LIGHT_AU_DAY,
                lookup,
            )?;
            xg = deflected[0];
            yg = deflected[1];
            zg = deflected[2];
            rg = (xg * xg + yg * yg + zg * zg).sqrt();
        }
        if rg > 0.0 {
            let ab = first_order_aberration_xyz([xg, yg, zg], vel, rg);
            xg = ab[0];
            yg = ab[1];
            zg = ab[2];
            rg = (xg * xg + yg * yg + zg * zg).sqrt();
        }
        let (lon_j, lat_j, _) = ecl_xyz_to_spherical([xg, yg, zg]);
        let (lon, lat) = crate::hypothetical::precess_ecliptic(
            lon_j,
            lat_j,
            crate::hypothetical::J2000_JD,
            jd_tt,
        )?;
        Ok((add_of_date_nutation_lon(self.serving, lon, jd_tt)?, lat, rg))
    }

    /// White Moon's symbolic geocentric ecliptic-of-date model, with
    /// first-order aberration and longitude nutation but no precession step.
    fn white_moon_geo_of_date(&self, jd_tt: f64) -> Result<(f64, f64, f64)> {
        let model = crate::hypothetical::calc_white_moon_position(jd_tt, false)?;
        self.symbolic_geo_of_date(model.lon_deg, model.lat_deg, model.dist_au, jd_tt)
    }

    /// Waldemath's symbolic model anchored at 1898 February 2, 00:00 UT.
    ///
    /// The bound Delta-T callback supplies the TT anchor. Subtracting the bound
    /// nutation in longitude from the stored apparent anchor gives its mean
    /// longitude. Uniform tropical motion is followed by first-order
    /// aberration and longitude nutation, without a precession step.
    fn waldemath_geo_of_date(&self, jd_tt: f64) -> Result<(f64, f64, f64)> {
        let anchor_ut = crate::hypothetical::WALDEMATH_ANCHOR_JD_UT;
        let dt_anchor = (self.deltat_at)(anchor_ut)?;
        if !dt_anchor.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_anchor });
        }
        let anchor_tt = anchor_ut + dt_anchor;
        let (dpsi_anchor_rad, _) = (self.serving.nutation_provider())(anchor_tt)?;
        if !dpsi_anchor_rad.is_finite() {
            return Err(Error::NonFiniteArgument {
                value: dpsi_anchor_rad,
            });
        }
        let anchor_mean_lon =
            crate::hypothetical::WALDEMATH_ANCHOR_APPARENT_LON_DEG - dpsi_anchor_rad.to_degrees();
        let state =
            crate::hypothetical::waldemath_state_with_anchor(jd_tt, anchor_tt, anchor_mean_lon)?;
        self.symbolic_geo_of_date(state.lon_deg, state.lat_deg, state.dist_au, jd_tt)
    }

    /// Convert symbolic ecliptic-of-date longitude/latitude in degrees and
    /// distance in AU through first-order annual aberration and longitude
    /// nutation. Nonfinite model coordinates and provider errors are returned.
    fn symbolic_geo_of_date(
        &self,
        lon_deg: f64,
        lat_deg: f64,
        dist_au: f64,
        jd_tt: f64,
    ) -> Result<(f64, f64, f64)> {
        if !lon_deg.is_finite() || !lat_deg.is_finite() || !dist_au.is_finite() {
            return Err(Error::NonFiniteArgument {
                value: if !lon_deg.is_finite() {
                    lon_deg
                } else if !lat_deg.is_finite() {
                    lat_deg
                } else {
                    dist_au
                },
            });
        }
        let lo_r = lon_deg.to_radians();
        let la_r = lat_deg.to_radians();
        let cb = la_r.cos();
        let geo = [
            dist_au * cb * lo_r.cos(),
            dist_au * cb * lo_r.sin(),
            dist_au * la_r.sin(),
        ];
        let vel = self.earth_bary_vel_ecl_j2000(jd_tt)?;
        let rg0 = (geo[0] * geo[0] + geo[1] * geo[1] + geo[2] * geo[2]).sqrt();
        let ab = first_order_aberration_xyz(geo, vel, rg0);
        let (lon, lat, rg) = ecl_xyz_to_spherical(ab);
        Ok((add_of_date_nutation_lon(self.serving, lon, jd_tt)?, lat, rg))
    }

    /// Observer position offset in the mean ecliptic-of-date frame, in AU.
    ///
    /// `jd_tt` and `jd_ut1` are Julian dates; longitude is east-positive degrees,
    /// latitude is degrees and height is metres. The WGS84 site uses `c2t06a`
    /// with zero polar motion, then J2000 and of-date rotations. Only position
    /// is consumed by the event search. Coordinate and rotation errors propagate.
    fn hypo_topo_offset_of_date(
        &self,
        jd_tt: f64,
        jd_ut1: f64,
        lon_deg_east: f64,
        lat_deg: f64,
        height_m: f64,
    ) -> Result<[f64; 3]> {
        let c2t = crate::erfa::c2t06a(jd_tt, jd_ut1, 0.0, 0.0)?;
        // This event route uses a WGS84 observer position; the offset's
        // velocity is not consumed by the altitude or transit searches.
        let off = crate::positions::topocentric_offset(
            lon_deg_east,
            lat_deg,
            height_m,
            c2t,
            crate::positions::SiteEllipsoid::Wgs84,
        )?;
        let (op, _ov) = off.pair();
        let rot = crate::positions::icrs_to_ecl_j2000_matrix();
        let oj = crate::positions::mat3_vec3(rot, op)?;
        rotate_vector_to_of_date(oj, jd_tt)
    }

    /// Apparent topocentric `(longitude, latitude, distance)` at a UT Julian date.
    ///
    /// The site uses east-positive longitude and latitude in degrees, height in
    /// metres, WGS84 and zero polar motion. Output angles are degrees and distance
    /// is AU. The of-date observer offset is subtracted from the geocentric
    /// apparent vector; a zero resulting distance returns the geocentric tuple.
    ///
    /// Input, Delta-T, model, provider and observer-transformation errors
    /// propagate from the corresponding calculations.
    pub fn fict_topo_ecliptic_at_ut(
        &self,
        jd_ut: f64,
        lon_deg_east: f64,
        lat_deg: f64,
        height_m: f64,
    ) -> Result<(f64, f64, f64)> {
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let (lon, lat, dist) = self.fict_geo_ecliptic_at_tt(jd_tt)?;
        let off = self.hypo_topo_offset_of_date(jd_tt, jd_ut, lon_deg_east, lat_deg, height_m)?;
        let lon_r = lon.to_radians();
        let lat_r = lat.to_radians();
        let cl = lat_r.cos();
        let bx = dist * cl * lon_r.cos();
        let by = dist * cl * lon_r.sin();
        let bz = dist * lat_r.sin();
        let tx = bx - off[0];
        let ty = by - off[1];
        let tz = bz - off[2];
        let r = (tx * tx + ty * ty + tz * tz).sqrt();
        if r == 0.0 {
            return Ok((lon, lat, dist));
        }
        let lon_t = ty.atan2(tx).to_degrees().rem_euclid(360.0);
        let lat_t = (tz / r).clamp(-1.0, 1.0).asin().to_degrees();
        Ok((lon_t, lat_t, r))
    }

    /// Apparent geocentric `(right ascension in hours, declination in degrees)`
    /// at a TT Julian date, using the bound nutation and true obliquity of date.
    /// Nonfinite-input, model, provider and obliquity errors propagate.
    pub fn fict_ra_dec_at_tt(&self, jd_tt: f64) -> Result<(f64, f64)> {
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let (lon, lat, _) = self.fict_geo_ecliptic_at_tt(jd_tt)?;
        let (_, deps_rad) = (self.serving.nutation_provider())(jd_tt)?;
        let eps_true_deg =
            crate::positions::vondrak_mean_obliquity_deg(jd_tt)? + deps_rad.to_degrees();
        let (ra_deg, dec_deg, _) = cotrans(lon, lat, 1.0, -eps_true_deg);
        Ok((ra_deg / 15.0, dec_deg))
    }

    /// Apparent geocentric `(right ascension in hours, declination in degrees)`
    /// at a UT Julian date. Uses the bound Delta-T callback and propagates its
    /// errors and those of the TT calculation.
    pub fn fict_ra_dec_at_ut(&self, jd_ut: f64) -> Result<(f64, f64)> {
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        self.fict_ra_dec_at_tt(jd_tt)
    }

    /// Greenwich apparent sidereal time in hours `[0, 24)` at a UT Julian date.
    ///
    /// Uses the bound Delta-T and nutation providers, true obliquity of date and
    /// the continuity offsets fixed at construction. Nonfinite-input, provider
    /// and sidereal-time calculation errors propagate.
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

    /// Local apparent sidereal time in hours `[0, 24)` at a UT Julian date.
    ///
    /// Longitude is east-positive degrees. Uses the bound Delta-T and nutation
    /// providers and the continuity offsets fixed at construction. A nonfinite
    /// longitude returns `InvalidObserver`; other input and calculation errors
    /// propagate.
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

    /// Find a geocentric apparent meridian transit from a UT Julian date.
    ///
    /// Longitude is east-positive degrees. `lower_transit` selects the opposite
    /// meridian, a 180-degree shift. The search delegates to [`calculate_transit`]
    /// using this object's right ascension and Greenwich apparent sidereal time;
    /// no topocentric parallax is applied. A nonfinite start returns
    /// `NonFiniteSearchValue`; callback and transit-search errors propagate.
    pub fn fict_transit_at_ut(
        &self,
        jd_start: f64,
        lon_deg_east: f64,
        lower_transit: bool,
    ) -> Result<TransitSolution> {
        if !jd_start.is_finite() {
            return Err(Error::NonFiniteSearchValue { x: jd_start });
        }
        let mut get_ra_dec = |jd: f64| self.fict_ra_dec_at_ut(jd);
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
