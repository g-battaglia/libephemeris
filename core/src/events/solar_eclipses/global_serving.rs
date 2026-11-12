//! Sealed global solar-eclipse provider adapters and public compositions.
//!
//! Global ground geometry intentionally acquires its frame nutation through
//! `ecl_nut_result_ut`, independently of the bound ephemeris Delta-T model.
//! Local/path-search adapters elsewhere have different acquisition contracts;
//! this module does not consolidate those routes or their Sun/Moon reads.

use super::{
    BesselianCore, ECL_RMOON_AU, ECL_RSUN_AU, EclipseLimitEphemeris, EclipsePathEphemeris,
    EclipsePathWidthEphemeris, EclipseSearchDirection, Error, PathWidthSample, Result,
    SealedSolarGlobalServing, SealedSolarLocalServing, SolHowCircumstances, SolarDurationEphemeris,
    SolarGlobalEphemeris, SolarLocalEphemeris, SunMoonState, SunMoonXyz, WhereGeometry,
    altaz_of_topo_ecliptic_impl, besselian_contact_residuals, besselian_core_via,
    calculate_local_eclipse_phases_via, eclipse_central_line_via, eclipse_northern_limit_via,
    eclipse_path_width_via, eclipse_southern_limit_via, eclipse_where_core_via, gamma_from_vectors,
    penumbra_limit_er, sol_eclipse_when_glob_via, sol_eclipse_where_parts, solar_contact_c1_via,
    solar_contact_c2_via, solar_contact_c3_via, solar_contact_c4_via,
    solar_eclipse_duration_min_via, solar_gamma_vectors_at_tt, umbra_limit_er,
    vincenty_separation_deg,
};

impl<'a, P> SealedSolarGlobalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Borrow a solar serving and a UT1-JD → TT-minus-UT1-days callback.
    ///
    /// SIDEREAL, J2000, NONUT, XYZ, RADIANS, HELCTR and BARYCTR return
    /// `UnsupportedFlagCombination`. Other calculation flags remain attached;
    /// this admission does not independently force every apparent correction.
    /// After flag admission, construction eagerly reads Delta-T at both modern
    /// sidereal-model boundary epochs to retain continuity offsets. Callback
    /// errors and nonfinite offsets propagate. No body state is acquired here;
    /// any callback I/O belongs to its supplied implementation.
    ///
    /// Geometry uses model TT for positions but an independent bundled frame
    /// clock for ground rotation/nutation, not a single interchangeable epoch.
    pub fn try_new(
        serving: &'a crate::positions::SealedServing<'a, P>,
        deltat_at: &'a dyn Fn(f64) -> Result<f64>,
        calc_flags: i32,
    ) -> Result<Self> {
        use crate::constants::{
            FLG_BARYCTR, FLG_HELCTR, FLG_J2000, FLG_NONUT, FLG_RADIANS, FLG_SIDEREAL, FLG_XYZ,
        };
        if (calc_flags & (FLG_SIDEREAL | FLG_J2000 | FLG_NONUT | FLG_XYZ | FLG_RADIANS)) != 0
            || (calc_flags & (FLG_HELCTR | FLG_BARYCTR)) != 0
        {
            return Err(Error::UnsupportedFlagCombination {
                entry: "sol_glob_serving",
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
            calc_flags,
            offset_lo_deg,
            offset_hi_deg,
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
    fn ecliptic_words_at_tt(&self, body_id: i32, jd_tt: f64) -> Result<[f64; 6]> {
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
    /// `j2000` selects direct mean-equator J2000 vectors for gamma,
    /// including frame bias and avoiding a spherical round trip.
    /// Otherwise returns true-of-date vectors reconstructed from the
    /// spherical words used by the ground-shadow and Besselian geometry.
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

    /// Greenwich apparent sidereal time in degrees at a UT1 epoch.
    /// Nutation and true obliquity use IAU 2006/2000A and Vondrák's
    /// long-term epsilon-A at the bundled-timescale TT epoch, not the
    /// caller's ephemeris model TT or a LEB nutation channel. The
    /// caller-selected Delta-T still drives the sidereal model and
    /// continuity offsets. This deliberate separation preserves the
    /// independent acquisition contracts of the two time/frame terms.
    fn gast_deg_at_ut(&self, jd_ut: f64) -> Result<f64> {
        use crate::time::{SiderealDeltaT, greenwich_apparent_sidereal_time_hours};
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let dt_days = (self.deltat_at)(jd_ut)?;
        if !dt_days.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_days });
        }
        let (nutation, _) = crate::positions::ecl_nut_result_ut(jd_ut, 0)?;
        let eps_true_deg = nutation[0];
        let dpsi_deg = nutation[2];
        let sdt = SiderealDeltaT::try_new(dt_days, self.offset_lo_deg, self.offset_hi_deg)?;
        Ok(greenwich_apparent_sidereal_time_hours(jd_ut, eps_true_deg, dpsi_deg, &sdt)? * 15.0)
    }

    /// Shadow-axis ground geometry at a UT epoch (site-independent
    /// global acquisition: release `_eclipse_where_core` over the
    /// sealed true-of-date vectors and GAST, `ECL_RSUN_AU` occulted
    /// radius). Public so the sealed `sol_eclipse_where` binding can
    /// name the shadow-center point.
    pub fn shadow_center_geometry_at_ut(&self, jd_ut: f64) -> Result<WhereGeometry> {
        self.where_geometry_at_ut(jd_ut)
    }

    /// Shadow-axis ground geometry at a UT epoch.
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
    /// a UT epoch over an explicit site triple.
    fn topo_ecliptic_at_ut(
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

    /// Where-magnitude at a UT epoch (same `_where_eclipse_magnitude`
    /// composition as the release: `_sol_how_core` `attr[0]` at the
    /// shadow-center ground point). Only the `attr[0]` words are composed
    /// here — the solar disc radius, the lunar disc radius and the
    /// Vincenty separation of the sealed topocentric places — because
    /// the `azalt` triple and the saros words feed the other `attr`
    /// slots, never the magnitude.
    fn where_magnitude_at_ut(&self, jd_ut: f64) -> Result<f64> {
        use crate::constants::{MOON, SUN};
        let geom = self.where_geometry_at_ut(jd_ut)?;
        let site = (geom.lon_deg(), geom.lat_deg(), 0.0);
        let (sun_lon, sun_lat, sun_dist) = self.topo_ecliptic_at_ut(SUN, jd_ut, site)?;
        let (moon_lon, moon_lat, moon_dist) = self.topo_ecliptic_at_ut(MOON, jd_ut, site)?;
        let rsun = (ECL_RSUN_AU / sun_dist).asin().to_degrees();
        let rmoon = (ECL_RMOON_AU / moon_dist).asin().to_degrees();
        let dctr = vincenty_separation_deg(sun_lon, sun_lat, moon_lon, moon_lat)?;
        Ok((rsun + rmoon - dctr) / (2.0 * rsun))
    }

    /// Fundamental-plane shadow geometry at a UT epoch.
    pub(in crate::events) fn besselian_core_at_ut(&self, jd_ut: f64) -> Result<BesselianCore> {
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let (sun_xyz, moon_xyz) = self.equatorial_xyz_at_tt(jd_tt, false)?;
        let gast_deg = self.gast_deg_at_ut(jd_ut)?;
        besselian_core_via(sun_xyz, moon_xyz, gast_deg)
    }
}

impl<'a, P> SolarGlobalEphemeris for SealedSolarGlobalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Apparent Sun/Moon ecliptic state with longitude speeds.
    fn sun_moon_state(&self, jd_ut: f64) -> Result<SunMoonState> {
        use crate::constants::{MOON, SUN};
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let sun = self.ecliptic_words_at_tt(SUN, jd_tt)?;
        let moon = self.ecliptic_words_at_tt(MOON, jd_tt)?;
        SunMoonState::try_new(
            sun[0], sun[1], sun[2], sun[3], moon[0], moon[1], moon[2], moon[3],
        )
    }

    /// Analytic LUN-01 mean lunar-node longitude: the analytic model is not sampled
    /// from a kernel, so the sealed LEB branch takes
    /// [`crate::lunar::LunarCoverage::Unenforced`].
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

    /// Fundamental-plane shadow radii `(l1, l2)` in Earth radii
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
        Ok(self.besselian_core_at_ut(jd_ut)?.x())
    }

    /// Ground-point core-shadow width (`dcore[0]`) in km for the hybrid
    /// probe.
    fn ground_core_width_km(&self, jd_ut: f64) -> Result<f64> {
        Ok(self.where_geometry_at_ut(jd_ut)?.dcore()[0])
    }
}

/// Next global solar eclipse from `jd_start_ut` over the sealed serving:
/// the release `sol_eclipse_when_glob` `(retflag, tret)` composition —
/// `retflag` is the classified type bitmask (`ECL_TOTAL` / `ECL_ANNULAR`
/// / `ECL_PARTIAL` / `ECL_ANNULAR_TOTAL` with `ECL_CENTRAL` /
/// `ECL_NONCENTRAL`; `0` accepts any type via `mask = 0`) and `tret` is
/// the 10-slot global phase layout (all Julian Days, UT):
/// `[0]` maximum eclipse, `[1]` central at local apparent noon,
/// `[2]`/`[3]` begin/end (penumbra touches/leaves Earth anywhere, `0`
/// when not found), `[4]`/`[5]` totality/annularity begin/end anywhere
/// (`0` when none), `[6]`/`[7]` central-line begin/end (`0` when none),
/// `[8]`/`[9]` hybrid total-to-annular transitions (`0`; not
/// implemented in the release). `direction` selects the forward,
/// backward or bidirectional walk. A walk that steps off the ephemeris
/// is [`Error::EclipseSearchExhausted`] (never a fabricated event); any
/// other provider failure propagates.
pub fn sol_eclipse_when_glob_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarGlobalServing<'_, P>,
    jd_start_ut: f64,
    mask: i32,
    direction: EclipseSearchDirection,
) -> Result<(i32, [f64; 10])> {
    sol_eclipse_when_glob_via(serving, jd_start_ut, mask, direction)
}

// Sealed solar-eclipse path servings

impl<'a, P> EclipsePathEphemeris for SealedSolarGlobalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Shadow-axis ground point.
    fn path_where(&self, jd_ut: f64) -> Result<(i32, f64, f64)> {
        let geom = self.where_geometry_at_ut(jd_ut)?;
        Ok((geom.retc(), geom.lon_deg(), geom.lat_deg()))
    }

    /// Besselian `(x, y)` in Earth radii.
    fn path_besselian_xy(&self, jd_ut: f64) -> Result<(f64, f64)> {
        let core = self.besselian_core_at_ut(jd_ut)?;
        Ok((core.x(), core.y()))
    }
}

impl<'a, P> EclipseLimitEphemeris for SealedSolarGlobalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// `(x_er, y_er, d_deg, mu_deg, l2_er)` at `jd_ut`.
    fn limit_besselian(&self, jd_ut: f64) -> Result<(f64, f64, f64, f64, f64)> {
        let core = self.besselian_core_at_ut(jd_ut)?;
        Ok((core.x(), core.y(), core.d_deg(), core.mu_deg(), core.l2()))
    }
}

/// Geographic points along the central line of a solar eclipse over the
/// sealed serving.
/// Returns `(times_jd_ut, latitudes_deg, longitudes_deg)`; provider
/// failures propagate.
pub fn sol_eclipse_central_line_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarGlobalServing<'_, P>,
    jd_start_ut: f64,
    jd_end_ut: f64,
    step_minutes: f64,
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>)> {
    eclipse_central_line_via(serving, jd_start_ut, jd_end_ut, step_minutes)
}

/// Points along the northern limit of a central eclipse path over the
/// sealed serving.
pub fn sol_eclipse_northern_limit_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarGlobalServing<'_, P>,
    jd_start_ut: f64,
    jd_end_ut: f64,
    step_minutes: f64,
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>)> {
    eclipse_northern_limit_via(serving, jd_start_ut, jd_end_ut, step_minutes)
}

/// Points along the southern limit of a central eclipse path over the
/// sealed serving.
pub fn sol_eclipse_southern_limit_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarGlobalServing<'_, P>,
    jd_start_ut: f64,
    jd_end_ut: f64,
    step_minutes: f64,
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>)> {
    eclipse_southern_limit_via(serving, jd_start_ut, jd_end_ut, step_minutes)
}

// Sealed solar path-width / duration servings

impl<'a, P> EclipsePathWidthEphemeris for SealedSolarGlobalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Shadow-axis ground point.
    fn path_width_where(&self, jd_ut: f64) -> Result<(i32, f64, f64)> {
        let geom = self.where_geometry_at_ut(jd_ut)?;
        Ok((geom.retc(), geom.lon_deg(), geom.lat_deg()))
    }

    /// Topocentric Sun/Moon sample at `(lon_deg, lat_deg)`, altitude 0
    fn path_width_sample(&self, jd_ut: f64, lon_deg: f64, lat_deg: f64) -> Result<PathWidthSample> {
        use crate::constants::{MOON, SUN};
        let site = (lon_deg, lat_deg, 0.0);
        let (sun_lon, sun_lat, sun_dist) = self.topo_ecliptic_at_ut(SUN, jd_ut, site)?;
        let (moon_lon, moon_lat, moon_dist) = self.topo_ecliptic_at_ut(MOON, jd_ut, site)?;
        let (_, sun_alt_true_deg, _) = altaz_of_topo_ecliptic_impl(
            (sun_lon, sun_lat, sun_dist),
            jd_ut,
            lon_deg,
            lat_deg,
            0.0,
            self.deltat_at,
            self.offset_lo_deg,
            self.offset_hi_deg,
        )?;
        let separation_deg = vincenty_separation_deg(sun_lon, sun_lat, moon_lon, moon_lat)?;
        Ok(PathWidthSample {
            sun_dist_au: sun_dist,
            moon_dist_au: moon_dist,
            sun_alt_true_deg,
            separation_deg,
        })
    }
}

impl<'a, P> SolarDurationEphemeris for SealedSolarGlobalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Shadow-axis ground point.
    fn duration_where(&self, jd_ut: f64) -> Result<(i32, f64, f64)> {
        let geom = self.where_geometry_at_ut(jd_ut)?;
        Ok((geom.retc(), geom.lon_deg(), geom.lat_deg()))
    }

    /// Local `(second_contact_jd, third_contact_jd)` at the site
    fn duration_local_contacts(
        &self,
        jd_ut: f64,
        lat_deg: f64,
        lon_deg: f64,
    ) -> Result<(f64, f64)> {
        use crate::constants::{MOON, SUN};
        let site = (lon_deg, lat_deg, 0.0);
        let mut separation = |jd: f64| -> Result<f64> {
            let (sun_lon, sun_lat, _) = self.topo_ecliptic_at_ut(SUN, jd, site)?;
            let (moon_lon, moon_lat, _) = self.topo_ecliptic_at_ut(MOON, jd, site)?;
            vincenty_separation_deg(sun_lon, sun_lat, moon_lon, moon_lat)
        };
        let mut sun_altaz = |jd: f64| -> Result<(f64, f64)> {
            let sun_topo = self.topo_ecliptic_at_ut(SUN, jd, site)?;
            let (_, true_alt, _) = altaz_of_topo_ecliptic_impl(
                sun_topo,
                jd,
                lon_deg,
                lat_deg,
                0.0,
                self.deltat_at,
                self.offset_lo_deg,
                self.offset_hi_deg,
            )?;
            // `calculate_local_eclipse_phases_via` consumes
            // `(true_altitude, azimuth)`; the azimuth word is unused on
            // this path, so only the altitude travels.
            Ok((true_alt, 0.0))
        };
        let mut distances = |jd: f64| -> Result<(f64, f64)> {
            let (_, _, sun_dist) = self.topo_ecliptic_at_ut(SUN, jd, site)?;
            let (_, _, moon_dist) = self.topo_ecliptic_at_ut(MOON, jd, site)?;
            Ok((sun_dist, moon_dist))
        };
        let phases = calculate_local_eclipse_phases_via(
            jd_ut,
            &mut separation,
            &mut sun_altaz,
            &mut distances,
        )?;
        Ok((phases.slots[2], phases.slots[3]))
    }
}

/// Approximate totality/annularity path width in km at a UT1 Julian day.
///
/// Both site coordinates are degrees, latitude north-positive and longitude
/// east-positive. A missing coordinate selects the central-line location.
/// The simplified projection, cap, zero-width sentinels and error policy of
/// [`eclipse_path_width_via`] apply unchanged: central-line acquisition errors
/// propagate, while selected sample failures degrade to zero. In particular,
/// the sample classifier examines direct range/corruption variants, not every
/// wrapped cause. This is not a certified perpendicular ground-track width.
pub fn sol_eclipse_path_width_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarGlobalServing<'_, P>,
    jd_ut: f64,
    lat_deg: Option<f64>,
    lon_deg: Option<f64>,
) -> Result<f64> {
    eclipse_path_width_via(serving, jd_ut, lat_deg, lon_deg)
}

/// Central-line maximum duration of totality/annularity in minutes at
/// `jd_max` over the sealed serving.
pub fn sol_eclipse_duration_min_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarGlobalServing<'_, P>,
    jd_max: f64,
) -> Result<f64> {
    solar_eclipse_duration_min_via(serving, jd_max)
}

/// Sealed global first-contact serving: the release
/// `calc_eclipse_first_contact_c1` public entry (the slot-0
/// ellipsoidal solve before maximum, matching `sol_eclipse_when_glob`
/// `tret[2]`) over the sealed solar contact-residual acquisition.
/// `jd_max` and the returned C1 are UT1 Julian days. This is first contact
/// anywhere on Earth, not first contact at an observer's site. Zero also marks
/// an absent bracket or a selected degraded acquisition/refinement, not only
/// a physically absent eclipse. Nonfinite `jd_max` and uncaught provider errors
/// propagate; residual units, bounded search and catches follow
/// [`solar_contact_c1_via`].
pub fn sol_contact_c1_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarGlobalServing<'_, P>,
    jd_max: f64,
) -> Result<f64> {
    if !jd_max.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_max });
    }
    solar_contact_c1_via(jd_max, &mut |jd| serving.contact_residuals_km(jd))
}

/// Sealed global second-contact serving: the release
/// `calc_eclipse_second_contact_c2` public entry (the slot-1
/// umbral/antumbral solve before maximum, matching
/// `sol_eclipse_when_glob` `tret[4]`) over the sealed solar
/// contact-residual acquisition. Input/output are UT1 Julian days. Zero
/// denotes absent or unresolved inner contact, including penumbra-only partials
/// and selected degraded failures; it does not describe an observer's local C2.
/// Error and bounded-search behavior follow [`solar_contact_c2_via`].
pub fn sol_contact_c2_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarGlobalServing<'_, P>,
    jd_max: f64,
) -> Result<f64> {
    if !jd_max.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_max });
    }
    solar_contact_c2_via(jd_max, &mut |jd| serving.contact_residuals_km(jd))
}

/// Sealed global third-contact serving: the release
/// `calc_eclipse_third_contact_c3` public entry (the slot-1
/// umbral/antumbral solve after maximum, matching
/// `sol_eclipse_when_glob` `tret[5]`) over the sealed solar
/// contact-residual acquisition. Input/output are UT1 Julian days. Zero
/// denotes absent or unresolved inner contact, including penumbra-only partials
/// and selected degraded failures; it does not describe an observer's local C3.
/// Error and bounded-search behavior follow [`solar_contact_c3_via`].
pub fn sol_contact_c3_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarGlobalServing<'_, P>,
    jd_max: f64,
) -> Result<f64> {
    if !jd_max.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_max });
    }
    solar_contact_c3_via(jd_max, &mut |jd| serving.contact_residuals_km(jd))
}

/// Sealed global fourth-contact serving: the release
/// `calc_eclipse_fourth_contact_c4` public entry (the slot-0
/// ellipsoidal solve after maximum, matching `sol_eclipse_when_glob`
/// `tret[3]`) over the sealed solar contact-residual acquisition.
/// Input/output are UT1 Julian days. This is last contact anywhere on Earth,
/// not a local observer contact. Zero marks an absent or unresolved contact;
/// selected degraded failures also use that sentinel. Nonfinite input and
/// uncaught provider errors propagate as on [`solar_contact_c4_via`].
pub fn sol_contact_c4_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarGlobalServing<'_, P>,
    jd_max: f64,
) -> Result<f64> {
    if !jd_max.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_max });
    }
    solar_contact_c4_via(jd_max, &mut |jd| serving.contact_residuals_km(jd))
}

/// Global solar-shadow location and local circumstances at a UT1 Julian day.
///
/// Returns `(retflag, geopos, attr)`. `retflag` is the global character from
/// [`WhereGeometry`]; `geopos[0..2]` are east-positive longitude and geodetic
/// latitude in degrees, followed by eight zeros. The twenty attributes follow
/// [`SolHowCircumstances`] at that ground point and zero height, with signed
/// core-shadow diameter in kilometres in `attr[3]`. Solar azimuth is south-zero
/// westward, not the heliacal north-zero convention. Raw no-event attribute
/// sentinels are retained; this is not the horizon-zeroing local-how wrapper.
/// Nonfinite input and consumed position/frame/clock errors propagate.
pub fn sol_eclipse_where_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarGlobalServing<'_, P>,
    jd_ut: f64,
) -> Result<(i32, [f64; 10], [f64; 20])> {
    if !jd_ut.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_ut });
    }
    let geom = serving.where_geometry_at_ut(jd_ut)?;
    let center = SealedSolarLocalServing::try_new(
        serving.serving,
        serving.deltat_at,
        serving.calc_flags,
        geom.lon_deg(),
        geom.lat_deg(),
        0.0,
    )?;
    let (retc, attr) = center.how_core(jd_ut)?;
    Ok(sol_eclipse_where_parts(
        &geom,
        &SolHowCircumstances { retc, attr },
    ))
}
