//! Site-bound solar-eclipse provider adapters and public compositions.
//!
//! This route retains the caller's Delta-T and bound nutation provider for its
//! local frame. It is intentionally independent of the global ground-geometry
//! adapter and of the existing path/search source-acquisition modules.

use super::{
    BodyHorizonCrossings, ECL_RMOON_AU, ECL_RMOON_INNER_AU, ECL_RSUN_AU, EclipseHowDetails,
    EclipseHowDetailsEphemeris, EclipseKind, EclipseSearchDirection, Error,
    RISE_SET_DIP_LAPSE_RATE_K_M, Result, RiseSetRequest, SealedSolarLocalServing,
    SealedTransitServing, SolHowCircumstances, SolarGlobalEphemeris, SolarLocalEphemeris,
    SunMoonState, SunMoonXyz, WhereGeometry, besselian_contact_residuals, besselian_core_via,
    calc_local_eclipse_max_time_via, cotrans, eclipse_disc_angular_radii_deg,
    eclipse_how_details_via, eclipse_where_core_via, gamma_from_vectors, penumbra_limit_er,
    rise_set_ut_serving, saros_info, sol_eclipse_how_parts, sol_eclipse_magnitude_at_loc_via,
    sol_eclipse_when_glob_via, sol_eclipse_when_loc_via, sol_how_core_via,
    solar_gamma_vectors_at_tt, umbra_limit_er, validate_geopos, vincenty_separation_deg,
};

impl<'a, P> SealedSolarLocalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Borrow a solar serving, model clock and explicit surface observer.
    ///
    /// The callback maps UT1 Julian day to TT-minus-UT1 in days. Longitude is
    /// east-positive and latitude north-positive, in degrees; height is metres.
    /// Observer admission checks finiteness only, not geographic ranges.
    /// SIDEREAL, J2000, NONUT, XYZ, RADIANS, HELCTR and BARYCTR are refused
    /// before observer checks. Construction then eagerly reads the two modern
    /// sidereal boundary offsets from the callback. Consumed callback errors
    /// and nonfinite offsets propagate. No body state is acquired; the serving
    /// and callback remain borrowed, with no ownership transfer.
    pub fn try_new(
        serving: &'a crate::positions::SealedServing<'a, P>,
        deltat_at: &'a dyn Fn(f64) -> Result<f64>,
        calc_flags: i32,
        lon_deg_east: f64,
        lat_deg: f64,
        alt_m: f64,
    ) -> Result<Self> {
        use crate::constants::{
            FLG_BARYCTR, FLG_HELCTR, FLG_J2000, FLG_NONUT, FLG_RADIANS, FLG_SIDEREAL, FLG_XYZ,
        };
        if (calc_flags & (FLG_SIDEREAL | FLG_J2000 | FLG_NONUT | FLG_XYZ | FLG_RADIANS)) != 0
            || (calc_flags & (FLG_HELCTR | FLG_BARYCTR)) != 0
        {
            return Err(Error::UnsupportedFlagCombination {
                entry: "sol_loc_serving",
                flags: calc_flags,
            });
        }
        validate_geopos(lon_deg_east, lat_deg, alt_m)?;
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
            calc_flags,
            offset_lo_deg,
            offset_hi_deg,
            lon_deg_east,
            lat_deg,
            alt_m,
        })
    }

    /// Bound sealed serving pair.
    pub fn serving(&self) -> &crate::positions::SealedServing<'a, P> {
        self.serving
    }

    /// Bound base calculation flags (without the per-acquisition shaping
    /// bits the release call sites add).
    pub fn calc_flags(&self) -> i32 {
        self.calc_flags
    }

    /// Bound observer triple `(longitude_deg_east, latitude_deg,
    /// altitude_m)` in release `geopos` order.
    pub fn observer_lon_lat_alt(&self) -> (f64, f64, f64) {
        (self.lon_deg_east, self.lat_deg, self.alt_m)
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
    /// a missing Sun/Moon/Earth channel
    /// inside asset coverage stays the typed miss, past coverage it
    /// names the requested body with the asset window.
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

    /// Apparent ecliptic words `(lon, lat, dist, dlon, dlat, ddist)` for
    /// one body at a TT epoch.
    pub(in crate::events) fn ecliptic_words_at_tt(
        &self,
        body_id: i32,
        jd_tt: f64,
    ) -> Result<[f64; 6]> {
        use crate::constants::FLG_SPEED;
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let (pos, _) = crate::positions::serving_compute_tt(
            self.serving,
            jd_tt,
            body_id,
            self.calc_flags | FLG_SPEED,
        )
        .map_err(|e| self.escalate(jd_tt, body_id, e))?;
        Ok(pos)
    }

    /// Geocentric apparent Sun/Moon Cartesian vectors in AU at TT.
    /// Direct biased J2000 mean-equator vectors serve gamma when
    /// `j2000` is true; otherwise true-of-date vectors reconstructed
    /// from spherical words serve ground-shadow geometry.
    fn equatorial_xyz_at_tt(&self, jd_tt: f64, j2000: bool) -> Result<SunMoonXyz> {
        if j2000 {
            return solar_gamma_vectors_at_tt(self.serving, jd_tt);
        }
        use crate::constants::{FLG_EQUATORIAL, FLG_XYZ, MOON, SUN};
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let shape = self.calc_flags | FLG_EQUATORIAL;
        let vector = |body_id: i32| -> Result<(f64, f64, f64)> {
            let (pos, _) =
                crate::positions::serving_compute_tt(self.serving, jd_tt, body_id, shape)
                    .map_err(|e| self.escalate(jd_tt, body_id, e))?;
            let xyz = crate::positions::apply_output_flags(pos, FLG_XYZ);
            Ok((xyz[0], xyz[1], xyz[2]))
        };
        Ok((vector(SUN)?, vector(MOON)?))
    }

    /// Greenwich apparent sidereal time in degrees at a UT epoch
    /// (same `sidtime(jd) * 15.0` acquisition as the release over the
    /// bound Delta-T provider and nutation pair).
    fn gast_deg_at_ut(&self, jd_ut: f64) -> Result<f64> {
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
        Ok(greenwich_apparent_sidereal_time_hours(jd_ut, eps_true_deg, dpsi_deg, &sdt)? * 15.0)
    }

    /// Shadow-axis ground geometry at a UT epoch (same `_eclipse_where_core`
    /// acquisition as the release over the sealed true-of-date vectors
    /// and GAST, `ECL_RSUN_AU` occulted radius).
    fn where_geometry_at_ut(&self, jd_ut: f64) -> Result<WhereGeometry> {
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let (sun_xyz, moon_xyz) = self.equatorial_xyz_at_tt(jd_tt, false)?;
        let gast_deg = self.gast_deg_at_ut(jd_ut)?;
        eclipse_where_core_via(moon_xyz, sun_xyz, gast_deg, ECL_RSUN_AU)
    }

    /// Topocentric apparent ecliptic `(lon, lat, dist)` for one body at
    /// a UT epoch over the bound observer triple.
    pub(in crate::events) fn topo_ecliptic_at_ut(
        &self,
        body_id: i32,
        jd_ut: f64,
    ) -> Result<(f64, f64, f64)> {
        use crate::constants::{FLG_SPEED, FLG_TOPOCTR};
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let site = (self.lon_deg_east, self.lat_deg, self.alt_m);
        let c2t_at = |tt: f64, ut1: f64| crate::erfa::c2t06a(tt, ut1, 0.0, 0.0);
        let (pos, _) = crate::positions::serving_compute_tt_topocentric(
            self.serving,
            jd_tt,
            jd_ut,
            body_id,
            self.calc_flags | FLG_TOPOCTR | FLG_SPEED,
            site,
            &c2t_at,
            false,
            crate::positions::SiteEllipsoid::Wgs84,
        )
        .map_err(|e| self.escalate(jd_tt, body_id, e))?;
        Ok((pos[0], pos[1], pos[2]))
    }

    /// Topocentric apparent ecliptic `(lon, lat, dist)` for one body at
    /// a UT epoch over an explicit site triple.
    fn topo_ecliptic_at_site_ut(
        &self,
        body_id: i32,
        jd_ut: f64,
        site: (f64, f64, f64),
    ) -> Result<(f64, f64, f64)> {
        use crate::constants::{FLG_SPEED, FLG_TOPOCTR};
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let c2t_at = |tt: f64, ut1: f64| crate::erfa::c2t06a(tt, ut1, 0.0, 0.0);
        let (pos, _) = crate::positions::serving_compute_tt_topocentric(
            self.serving,
            jd_tt,
            jd_ut,
            body_id,
            self.calc_flags | FLG_TOPOCTR | FLG_SPEED,
            site,
            &c2t_at,
            false,
            crate::positions::SiteEllipsoid::Wgs84,
        )
        .map_err(|e| self.escalate(jd_tt, body_id, e))?;
        Ok((pos[0], pos[1], pos[2]))
    }

    /// Topocentric apparent distances `(sun_dist_au, moon_dist_au)` at
    /// a UT epoch over the bound observer triple.
    pub(in crate::events) fn topo_distances_au_at_ut(&self, jd_ut: f64) -> Result<(f64, f64)> {
        use crate::constants::{MOON, SUN};
        let (_, _, sun_dist) = self.topo_ecliptic_at_ut(SUN, jd_ut)?;
        let (_, _, moon_dist) = self.topo_ecliptic_at_ut(MOON, jd_ut)?;
        Ok((sun_dist, moon_dist))
    }

    /// Where-magnitude at a UT epoch (same `_where_eclipse_magnitude`
    /// composition as the release: `_sol_how_core` `attr[0]` at the
    /// shadow-center ground point — the solar/lunar disc radii and the
    /// Vincenty separation of the sealed topocentric places).
    fn where_magnitude_at_ut(&self, jd_ut: f64) -> Result<f64> {
        use crate::constants::{MOON, SUN};
        let geom = self.where_geometry_at_ut(jd_ut)?;
        let site = (geom.lon_deg(), geom.lat_deg(), 0.0);
        let (sun_lon, sun_lat, sun_dist) = self.topo_ecliptic_at_site_ut(SUN, jd_ut, site)?;
        let (moon_lon, moon_lat, moon_dist) = self.topo_ecliptic_at_site_ut(MOON, jd_ut, site)?;
        let rsun = (ECL_RSUN_AU / sun_dist).asin().to_degrees();
        let rmoon = (ECL_RMOON_AU / moon_dist).asin().to_degrees();
        let dctr = vincenty_separation_deg(sun_lon, sun_lat, moon_lon, moon_lat)?;
        Ok((rsun + rmoon - dctr) / (2.0 * rsun))
    }

    /// Solar horizontal triple `(true_altitude_deg, azimuth_deg,
    /// apparent_altitude_deg)` at a UT epoch.
    fn sun_altaz_at_ut(&self, jd_ut: f64) -> Result<(f64, f64, f64)> {
        use crate::constants::SUN;
        use crate::time::{SiderealDeltaT, local_apparent_sidereal_time_hours};
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let (ecl_lon, ecl_lat, dist) = self.topo_ecliptic_at_ut(SUN, jd_ut)?;
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let (dpsi_rad, deps_rad) = (self.serving.nutation_provider())(jd_tt)?;
        let eps_true_deg =
            crate::positions::vondrak_mean_obliquity_deg(jd_tt)? + deps_rad.to_degrees();
        let (ra_deg, dec_deg, _) = cotrans(ecl_lon, ecl_lat, dist, -eps_true_deg);
        let dt_days = (self.deltat_at)(jd_ut)?;
        if !dt_days.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_days });
        }
        let sdt = SiderealDeltaT::try_new(dt_days, self.offset_lo_deg, self.offset_hi_deg)?;
        let last_hours = local_apparent_sidereal_time_hours(
            jd_ut,
            self.lon_deg_east,
            eps_true_deg,
            dpsi_rad.to_degrees(),
            &sdt,
        )?;
        let ha_deg = (last_hours * 15.0 - ra_deg).rem_euclid(360.0);
        let ha_rad = ha_deg.to_radians();
        let dec_rad = dec_deg.to_radians();
        let lat_rad = self.lat_deg.to_radians();
        let sin_alt = lat_rad.sin() * dec_rad.sin() + lat_rad.cos() * dec_rad.cos() * ha_rad.cos();
        let true_alt_deg = sin_alt.clamp(-1.0, 1.0).asin().to_degrees();
        let y = ha_rad.sin() * dec_rad.cos();
        let x = ha_rad.cos() * dec_rad.cos() * lat_rad.sin() - dec_rad.sin() * lat_rad.cos();
        let azimuth_deg = y.atan2(x).to_degrees().rem_euclid(360.0);
        let pressure_mbar = crate::refraction::estimate_pressure_mbar(self.alt_m);
        let (apparent_alt_deg, _) = crate::refraction::refrac_extended_true_to_app(
            true_alt_deg,
            self.alt_m,
            pressure_mbar,
            10.0,
            RISE_SET_DIP_LAPSE_RATE_K_M,
        );
        if !apparent_alt_deg.is_finite() {
            return Err(Error::NonFiniteArgument {
                value: apparent_alt_deg,
            });
        }
        Ok((true_alt_deg, azimuth_deg, apparent_alt_deg))
    }
}

impl<'a, P> SolarGlobalEphemeris for SealedSolarLocalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Apparent Sun/Moon ecliptic state with longitude speeds: the
    /// candidate-lunation acquisition behind
    /// [`SolarLocalEphemeris::next_global_forward`].
    fn sun_moon_state(&self, jd_ut: f64) -> Result<SunMoonState> {
        use crate::constants::{MOON, SUN};
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let sun = self.ecliptic_words_at_tt(SUN, jd_tt)?;
        let moon = self.ecliptic_words_at_tt(MOON, jd_tt)?;
        SunMoonState::try_new(
            sun[0], sun[1], sun[2], sun[3], moon[0], moon[1], moon[2], moon[3],
        )
    }

    /// Analytic LUN-01 mean lunar-node longitude.
    fn mean_node_lon_deg(&self, jd_ut: f64) -> Result<f64> {
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        crate::lunar::calc_mean_lunar_node(jd_tt, crate::lunar::LunarCoverage::Unenforced)
    }

    /// Gamma in Earth radii.
    fn gamma_er(&self, jd_ut: f64) -> Result<f64> {
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let (sun_xyz, moon_xyz) = self.equatorial_xyz_at_tt(jd_tt, true)?;
        gamma_from_vectors(sun_xyz, moon_xyz)
    }

    /// Fundamental-plane shadow radii `(l1, l2)` in Earth radii.
    fn shadow_limits_er(&self, jd_ut: f64) -> Result<(f64, f64)> {
        let state = self.sun_moon_state(jd_ut)?;
        Ok((
            penumbra_limit_er(state.sun_dist_au(), state.moon_dist_au())?,
            umbra_limit_er(state.sun_dist_au(), state.moon_dist_au())?,
        ))
    }

    /// Where-core retc plus where-magnitude (the graze-boundary gate
    /// words behind [`super::classify_solar_eclipse`]).
    fn where_gate(&self, jd_ut: f64) -> Result<(i32, f64)> {
        let geom = self.where_geometry_at_ut(jd_ut)?;
        let magnitude = self.where_magnitude_at_ut(jd_ut)?;
        Ok((geom.retc(), magnitude))
    }

    /// `(pen, umb, cen)` contact residuals in km.
    fn contact_residuals_km(&self, jd_ut: f64) -> Result<(f64, f64, f64)> {
        let geom = self.where_geometry_at_ut(jd_ut)?;
        let residuals = besselian_contact_residuals(geom.dcore())?;
        Ok((residuals.pen_km, residuals.umb_km, residuals.cen_km))
    }

    /// Besselian-x in Earth radii.
    fn besselian_x_er(&self, jd_ut: f64) -> Result<f64> {
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let (sun_xyz, moon_xyz) = self.equatorial_xyz_at_tt(jd_tt, false)?;
        let gast_deg = self.gast_deg_at_ut(jd_ut)?;
        Ok(besselian_core_via(sun_xyz, moon_xyz, gast_deg)?.x())
    }

    /// Ground-point core-shadow width (`dcore[0]`) in km for the hybrid
    /// probe.
    fn ground_core_width_km(&self, jd_ut: f64) -> Result<f64> {
        Ok(self.where_geometry_at_ut(jd_ut)?.dcore()[0])
    }
}

impl<'a, P> SolarLocalEphemeris for SealedSolarLocalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Next global solar maximum at or after `jd`: a
    /// coverage miss surfaces through the walk's exhaustion conversion.
    fn next_global_forward(&self, jd: f64) -> Result<(i32, [f64; 10])> {
        sol_eclipse_when_glob_via(self, jd, 0, EclipseSearchDirection::Forward)
    }

    /// Topocentric Sun-Moon center separation, degrees.
    fn separation_deg(&self, jd: f64) -> Result<f64> {
        use crate::constants::{MOON, SUN};
        let (sun_lon, sun_lat, _) = self.topo_ecliptic_at_ut(SUN, jd)?;
        let (moon_lon, moon_lat, _) = self.topo_ecliptic_at_ut(MOON, jd)?;
        vincenty_separation_deg(sun_lon, sun_lat, moon_lon, moon_lat)
    }

    /// Apparent radii plus center separation `(rsun, rmoon, sep)`,
    /// degrees.
    fn radii_sep_deg(&self, jd: f64, inner_moon: bool) -> Result<(f64, f64, f64)> {
        use crate::constants::{MOON, SUN};
        let (sun_lon, sun_lat, sun_dist) = self.topo_ecliptic_at_ut(SUN, jd)?;
        let (moon_lon, moon_lat, moon_dist) = self.topo_ecliptic_at_ut(MOON, jd)?;
        let rmoon_au = if inner_moon {
            ECL_RMOON_INNER_AU
        } else {
            ECL_RMOON_AU
        };
        let rsun = (ECL_RSUN_AU / sun_dist).asin().to_degrees();
        let rmoon = (rmoon_au / moon_dist).asin().to_degrees();
        let sep = vincenty_separation_deg(sun_lon, sun_lat, moon_lon, moon_lat)?;
        Ok((rsun, rmoon, sep))
    }

    /// Solar `(true_altitude_deg, azimuth_deg, apparent_altitude_deg)`
    fn sun_altaz_deg(&self, jd: f64) -> Result<(f64, f64, f64)> {
        self.sun_altaz_at_ut(jd)
    }

    /// Topocentric apparent diameters `(moon_diam_deg, sun_diam_deg)`
    /// (full-disc local eclipse layout
    /// path; the `_sol_eclipse_when_loc_impl` walk ported by
    /// [`sol_eclipse_when_loc_via`] does not consume this word — the
    /// diameters travel inside [`SolarLocalEphemeris::how_core`]).
    fn topo_diameters_deg(&self, jd: f64) -> Result<(f64, f64)> {
        use crate::constants::{MOON, SUN};
        let (_, _, sun_dist) = self.topo_ecliptic_at_ut(SUN, jd)?;
        let (_, _, moon_dist) = self.topo_ecliptic_at_ut(MOON, jd)?;
        let sun_diam = 2.0 * (ECL_RSUN_AU / sun_dist).asin().to_degrees();
        let moon_diam = 2.0 * (ECL_RMOON_AU / moon_dist).asin().to_degrees();
        Ok((moon_diam, sun_diam))
    }

    /// Lower-limb solar horizon crossings from `jd_from`.
    fn sun_horizon_crossings(&self, jd_from: f64) -> Result<BodyHorizonCrossings> {
        use crate::constants::{BIT_DISC_BOTTOM, CALC_RISE, CALC_SET, SUN};
        let sun =
            SealedTransitServing::try_new(self.serving, self.deltat_at, SUN, self.calc_flags)?;
        let request = |rsmi: i32| RiseSetRequest {
            jd_start_ut: jd_from,
            lon_deg_east: self.lon_deg_east,
            lat_deg: self.lat_deg,
            height_m: self.alt_m,
            atpress_mbar: 0.0,
            attemp_c: 0.0,
            horhgt_deg: 0.0,
            rsmi,
        };
        let (rise_rc, rise_tret) =
            rise_set_ut_serving(&sun, &request(CALC_RISE | BIT_DISC_BOTTOM))?;
        let (set_rc, set_tret) = rise_set_ut_serving(&sun, &request(CALC_SET | BIT_DISC_BOTTOM))?;
        Ok(BodyHorizonCrossings {
            rise_rc,
            rise_jd: rise_tret[0],
            set_rc,
            set_jd: set_tret[0],
        })
    }

    /// Solar local-circumstance core at `jd`: `(retc, attr20)` in the
    /// `_sol_how_core` layout.
    fn how_core(&self, jd: f64) -> Result<(i32, [f64; 20])> {
        use crate::constants::{MOON, SUN};
        let sun_topo = self.topo_ecliptic_at_ut(SUN, jd)?;
        let moon_topo = self.topo_ecliptic_at_ut(MOON, jd)?;
        let (true_alt, azimuth, apparent_alt) = self.sun_altaz_at_ut(jd)?;
        let saros = saros_info(jd, EclipseKind::Solar);
        let how = sol_how_core_via(
            sun_topo,
            moon_topo,
            ECL_RSUN_AU,
            true,
            false,
            (azimuth, true_alt, apparent_alt),
            self.alt_m,
            saros,
        )?;
        Ok((how.retc(), how.attr()))
    }

    /// Shadow-axis core at `jd`: `(retc, dcore10)`.
    fn where_core(&self, jd: f64) -> Result<(i32, [f64; 10])> {
        let geom = self.where_geometry_at_ut(jd)?;
        Ok((geom.retc(), geom.dcore()))
    }
}

/// Next solar eclipse visible from the bound observer triple, from
/// `jd_start_ut` over the sealed serving: the release
/// `sol_eclipse_when_loc` `(ecl_type, tret, attr)` composition —
/// `tret` is the 10-slot local phase layout (all Julian Days, UT):
/// `[0]` local maximum (or the horizon crossing when the geometric
/// maximum is down), `[1]`/`[4]` outer contacts, `[2]`/`[3]` inner
/// contacts (`0.0` when none), `[5]`/`[6]` sunrise/sunset inside the
/// contacts (`0.0` when none), `[7]`–`[9]` reserved (`0.0`); `attr`
/// is the 20-slot `_sol_how_core` array with `attr[3]` carrying the
/// ground core-shadow width in km. `backwards` selects the forward
/// (`false`) or backward (`true`) walk. A walk that finds nothing returns
/// [`Error::EclipseSearchExhausted`]. Any forward global-candidate failure is
/// converted to that error; any backward candidate failure ends its current
/// collection sweep, retaining candidates already found. Coverage failure in
/// local-maximum refinement also becomes exhaustion. Other failures consumed
/// by local refinement, contacts, visibility and circumstance assembly propagate.
pub fn sol_eclipse_when_loc_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarLocalServing<'_, P>,
    jd_start_ut: f64,
    backwards: bool,
) -> Result<(i32, [f64; 10], [f64; 20])> {
    sol_eclipse_when_loc_via(serving, jd_start_ut, backwards)
}

/// Sealed local-circumstance serving: the release `sol_eclipse_how`
/// public entry (`_sol_eclipse_how_impl` tail) over the sealed local
/// solar serving. Returns `(retflag, attr[20])`: the local eclipse
/// character (`ECL_TOTAL` / `ECL_ANNULAR` / `ECL_PARTIAL`, plus
/// `ECL_VISIBLE` when part of the eclipsed Sun stands above the local
/// horizon, plus the global `ECL_CENTRAL` / `ECL_NONCENTRAL` character
/// at this instant) with the 20-slot circumstance array — magnitude
/// as solar-diameter fraction in `attr[0]`, the lunar/solar diameter
/// ratio in `attr[1]`, obscuration in `attr[2]`, the
/// observer-independent core-shadow diameter at the global
/// shadow-center point in `attr[3]`, the Sun azalt triple in
/// `attr[4]`/`[5]`/`[6]`, the Moon-Sun separation in `attr[7]`, the
/// NASA magnitude in `attr[8]` and the saros pair in `attr[9]`/`[10]`.
/// Below the apparent horizon (`attr[6] <= 0.0`) nothing of the
/// eclipse is observable and `retflag` is 0; when `retflag` is 0 the
/// slots 0, 1, 2, 3, 8, 9, 10 are zeroed while the
/// azimuth/altitude/separation values in slots 4..7 are kept.
/// All timestamps are UT1 Julian days. Angles, including separation, are
/// degrees; azimuth is south-zero westward, true altitude is unrefracted and
/// apparent altitude includes atmospheric refraction. `attr[2]` is uncapped
/// under the solar totality convention and can exceed one; it is not uniformly
/// a covered-area fraction bounded by one. Slots 11–19 are reserved zeros.
/// Nonfinite instants and consumed provider/frame/clock failures return errors.
pub fn sol_how_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarLocalServing<'_, P>,
    jd_ut: f64,
) -> Result<(i32, [f64; 20])> {
    if !jd_ut.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_ut });
    }
    let (retc, attr) = serving.how_core(jd_ut)?;
    let geom = serving.where_geometry_at_ut(jd_ut)?;
    Ok(sol_eclipse_how_parts(
        &SolHowCircumstances { retc, attr },
        &geom,
    ))
}

// Sealed solar fixed-instant local servings

/// Sealed fixed-instant solar-eclipse magnitude at the bound site:
/// the release `sol_eclipse_magnitude_at_loc` public entry over the
/// sealed local solar serving. Returns the fraction of the Sun's
/// diameter covered by the Moon (`0.0` when the Sun stands below the
/// `-1.0` deg gate or the discs do not overlap, `>= 1.0` at and past
/// totality). Non-finite instants refuse.
pub fn sol_eclipse_magnitude_at_loc_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarLocalServing<'_, P>,
    jd_ut: f64,
) -> Result<f64> {
    if !jd_ut.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_ut });
    }
    let (true_alt_deg, _, _) = serving.sun_altaz_deg(jd_ut)?;
    let separation_deg = serving.separation_deg(jd_ut)?;
    let (sun_dist_au, moon_dist_au) = serving.topo_distances_au_at_ut(jd_ut)?;
    sol_eclipse_magnitude_at_loc_via(true_alt_deg, separation_deg, sun_dist_au, moon_dist_au)
}

/// Sealed local solar maximum time near `jd_approx` for the bound site:
/// the release `sol_eclipse_max_time` local leg (golden-section minimum
/// of the topocentric Sun-Moon separation over `jd_approx ±
/// search_range_days`, 60 iterations, `1e-8` day width stop). Returns
/// `(jd_maximum, min_separation_degrees)`, with the date in UT1 Julian days.
/// A finite approximate epoch and finite positive half-window are required;
/// invalid windows, nonfinite sampled separations and provider failures return
/// typed errors. Minimizing separation does not itself prove an eclipse exists.
pub fn sol_local_max_time_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarLocalServing<'_, P>,
    jd_approx: f64,
    search_range_days: f64,
) -> Result<(f64, f64)> {
    use crate::constants::{MOON, SUN};
    let mut separation_at = |jd: f64| {
        let (sun_lon, sun_lat, _) = serving.topo_ecliptic_at_ut(SUN, jd)?;
        let (moon_lon, moon_lat, _) = serving.topo_ecliptic_at_ut(MOON, jd)?;
        vincenty_separation_deg(sun_lon, sun_lat, moon_lon, moon_lat)
    };
    calc_local_eclipse_max_time_via(&mut separation_at, jd_approx, search_range_days)
}

// Sealed solar how-details serving

impl<'a, P> SealedSolarLocalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Topocentric apparent equatorial `(ra_deg, dec_deg)` for one body
    /// at a UT epoch: the release `_calc_position_angle` acquisition
    /// order (`calc_ut(jd, body, FLG_EQUATORIAL | FLG_J2000 |
    /// FLG_TOPOCTR | FLG_SPEED)` words over the bound site), not the
    /// of-date `cotrans` projection — the contact position angles are
    /// measured in the J2000 frame on both sides.
    fn topo_equatorial_at_ut(&self, body_id: i32, jd_ut: f64) -> Result<(f64, f64)> {
        use crate::constants::{FLG_EQUATORIAL, FLG_J2000, FLG_SPEED, FLG_TOPOCTR};
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let site = (self.lon_deg_east, self.lat_deg, self.alt_m);
        let c2t_at = |tt: f64, ut1: f64| crate::erfa::c2t06a(tt, ut1, 0.0, 0.0);
        let (pos, _) = crate::positions::serving_compute_tt_topocentric(
            self.serving,
            jd_tt,
            jd_ut,
            body_id,
            self.calc_flags | FLG_EQUATORIAL | FLG_J2000 | FLG_TOPOCTR | FLG_SPEED,
            site,
            &c2t_at,
            false,
            crate::positions::SiteEllipsoid::Wgs84,
        )
        .map_err(|e| self.escalate(jd_tt, body_id, e))?;
        Ok((pos[0], pos[1]))
    }
}

impl<'a, P> EclipseHowDetailsEphemeris for SealedSolarLocalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// `(eclipse_type, attr20)` at the requested instant: the release
    /// `_sol_eclipse_how_impl` tail over the raw core words.
    /// The details composition consumes the public-entry pair, not the
    /// `_sol_how_core` layout — the global central/noncentral
    /// character ORs into a nonzero type, `attr[3]` is the where-core
    /// shadow diameter (signed: negative when the umbra reaches the
    /// surface), and below the apparent horizon the type zeroes with
    /// the magnitude/ratio/obscuration/width/saros slots. (`how_core`
    /// itself stays raw: the when_loc walk composes these same words
    /// at its own call site.)
    fn details_initial(&self, jd_ut: f64) -> Result<(i32, [f64; 20])> {
        use crate::constants::{ECL_CENTRAL, ECL_NONCENTRAL};
        let (retc, mut attr) = self.how_core(jd_ut)?;
        let (rc_where, dcore) = self.where_core(jd_ut)?;
        let mut retflag = retc;
        if retflag != 0 {
            retflag |= rc_where & (ECL_CENTRAL | ECL_NONCENTRAL);
        }
        attr[3] = dcore[0];
        // Below the apparent horizon nothing of the eclipse is
        // observable (`NaN <= 0.0` is false, so a non-finite altitude
        // keeps the type exactly like the release comparison).
        if attr[6] <= 0.0 {
            retflag = 0;
        }
        if retflag == 0 {
            for i in [0, 1, 2, 3, 8, 9, 10] {
                attr[i] = 0.0;
            }
        }
        Ok((retflag, attr))
    }

    /// Topocentric Sun-Moon separation in degrees.
    fn details_separation_deg(&self, jd_ut: f64) -> Result<f64> {
        self.separation_deg(jd_ut)
    }

    /// True Sun `(altitude, azimuth)` in degrees.
    fn details_sun_alt_az_deg(&self, jd_ut: f64) -> Result<(f64, f64)> {
        let (true_alt, azimuth, _) = self.sun_altaz_at_ut(jd_ut)?;
        Ok((true_alt, azimuth))
    }

    /// Apparent `(sun_radius, moon_radius)` in degrees.
    fn details_disc_radii_deg(&self, jd_ut: f64) -> Result<(f64, f64)> {
        use crate::constants::{MOON, SUN};
        let (_, _, sun_dist) = self.topo_ecliptic_at_ut(SUN, jd_ut)?;
        let (_, _, moon_dist) = self.topo_ecliptic_at_ut(MOON, jd_ut)?;
        eclipse_disc_angular_radii_deg(sun_dist, moon_dist)
    }

    /// Equatorial `(sun_ra, sun_dec, moon_ra, moon_dec)` in degrees for
    /// the contact position angles.
    fn details_equatorial_deg(&self, jd_ut: f64) -> Result<(f64, f64, f64, f64)> {
        use crate::constants::{MOON, SUN};
        let (sun_ra, sun_dec) = self.topo_equatorial_at_ut(SUN, jd_ut)?;
        let (moon_ra, moon_dec) = self.topo_equatorial_at_ut(MOON, jd_ut)?;
        Ok((sun_ra, sun_dec, moon_ra, moon_dec))
    }
}

/// Comprehensive local solar-eclipse circumstances at `jd_ut` over the
/// site-bound sealed serving. The site travels inside the bound serving.
pub fn sol_eclipse_how_details_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarLocalServing<'_, P>,
    jd_ut: f64,
) -> Result<EclipseHowDetails> {
    let (lon, lat, alt) = serving.observer_lon_lat_alt();
    eclipse_how_details_via(serving, jd_ut, lon, lat, alt)
}
