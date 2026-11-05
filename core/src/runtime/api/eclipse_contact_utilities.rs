//! Fixed-maximum lunar contacts and central-line solar duration.
use super::solar_circumstances::{SelectedMain, solar_position, solar_selected};
use super::*;
use crate::events::{self, LunarShadowDcore};

/// Provider stage at which a solar-duration calculation failed.
/// Native error values retain the actual file, domain and numerical cause.
#[derive(Debug)]
pub enum SolarDurationError {
    /// The genuine public where calculation failed before local phases.
    Where(Error),
    /// The fixed-local reader branch failed, after actual probe acquisition.
    Reader {
        /// Actual primary data or numerical error.
        error: Error,
        /// True only for the forced stored-reader mode of this acquisition.
        sealed: bool,
        /// Actual bound provider kind, absent if reader construction failed.
        tiered: Option<bool>,
    },
    /// The direct native local observer/phase calculation failed.
    Native(Error),
}

impl Session {
    /// Lunar contact at fixed maximum UT1 JD. Slots6/7 are P1/P4, slots2/3
    /// U1/U4 and slots4/5 U2/U3. Each call classifies the canonical geocentric
    /// shadow, then computes the complete original ten-slot contact array.
    /// Missing events/phases return positive zero. Moon then Sun public
    /// Cartesian EQUATORIAL|XYZ|SPEED reads precede Danjon geometry and Saros.
    /// Original ephemeris bits, scalar-model clocks and acquired errors remain
    /// unchanged. One selected native image survives every bisection sample.
    pub fn lunar_eclipse_contact_at_ut(
        &self,
        jd_max: f64,
        flags: i32,
        slot: usize,
    ) -> Result<f64, Error> {
        let selected = SelectedMain::new();
        events::lun_contact_via_core(jd_max, slot, &mut |jd| {
            self.utility_lunar_shadow(jd, flags, &selected)
        })
    }

    /// Totality U3−U2 in minutes at the supplied maximum UT1 JD. U2 is
    /// obtained by a complete phase calculation; zero stops before U3.
    /// Nonzero U2 then reacquires a complete U3 calculation over the same
    /// owned main image. The quotient's final day-to-minute conversion is
    /// `(U3-U2)*24*60`, preserving its two multiplication stages.
    pub fn lunar_eclipse_total_duration_at_ut(
        &self,
        jd_max: f64,
        flags: i32,
    ) -> Result<f64, Error> {
        let selected = SelectedMain::new();
        events::lunar_total_duration_min_via(jd_max, &mut |jd| {
            self.utility_lunar_shadow(jd, flags, &selected)
        })
    }

    /// Full umbral U4−U1 duration in minutes. U1 then U4 each use the full
    /// original phase solve, with the same zero short circuits and actual
    /// error/clock acquisition as totality. No event is searched or shifted.
    pub fn lunar_eclipse_umbral_duration_at_ut(
        &self,
        jd_max: f64,
        flags: i32,
    ) -> Result<f64, Error> {
        let selected = SelectedMain::new();
        events::lunar_umbral_duration_min_via(jd_max, &mut |jd| {
            self.utility_lunar_shadow(jd, flags, &selected)
        })
    }

    /// Full geocentric shadow at original UT1, including the unused speed
    /// acquisition and Saros metadata of the underlying public how operation.
    /// Neither observer settings nor sidereal continuity clocks are consumed.
    fn utility_lunar_shadow(
        &self,
        jd: f64,
        flags: i32,
        selected: &SelectedMain,
    ) -> Result<(i32, LunarShadowDcore), Error> {
        let base = (flags & 7) | FLG_EQUATORIAL | FLG_XYZ | FLG_SPEED;
        let position = |body| {
            if self.ctx.mode() == ProviderMode::NativeJpl {
                solar_position(self, jd, body, base, None, selected)
            } else {
                self.calc_ut(jd, body, base).map(|p| p.values)
            }
        };
        lunar_shadow_from_positions(jd, &position)
    }

    /// Central-line local C3−C2 duration in minutes at fixed maximum UT1 JD.
    /// Complete public where precedes the central-bit test. A central event
    /// uses its actual shadow-center longitude/latitude at height0, then the
    /// accepted fixed-local phase kernel; there is no next-event search.
    /// Caller flags belong only to where. Local scalar LEB uses one common
    /// configured TT per paired callback and no SPEED. Direct native callbacks
    /// use fresh independent UT1 Time and WGS84 rays; altitude acquires only
    /// Sun and the SDK ITRS horizon. Reader failures retain their local role
    /// for attached mode-aware public exception construction.
    pub fn solar_eclipse_duration_at_ut(
        &self,
        jd_max: f64,
        flags: i32,
    ) -> Result<f64, SolarDurationError> {
        let selected = SelectedMain::new();
        let (retflag, geopos, _) = if self.ctx.mode() == ProviderMode::NativeJpl {
            self.sol_eclipse_where_with_position_context(jd_max, flags, &|_, body| {
                solar_selected(self, body, None, &selected).map(|s| s.ctx.clone())
            })
        } else {
            self.sol_eclipse_where(jd_max, flags)
        }
        .map_err(SolarDurationError::Where)?;
        if retflag & ECL_CENTRAL == 0 {
            return Ok(0.0);
        }
        let owner = selected.get().unwrap_or(self);
        let phases = owner.utility_fixed_solar_phases(jd_max, [geopos[0], geopos[1], 0.0])?;
        Ok(central_duration_minutes(phases[2], phases[3]))
    }

    /// Auto retries the local solve only after a genuine stored range/body
    /// miss. Corruption and unrelated failures
    /// never become sentinels or fallback, and forced LEB never selects JPL.
    fn utility_fixed_solar_phases(
        &self,
        jd: f64,
        site: [f64; 3],
    ) -> Result<[f64; 10], SolarDurationError> {
        let reader = self
            .utility_reader_owner()
            .map_err(|error| SolarDurationError::Reader {
                error,
                sealed: self.ctx.mode() == ProviderMode::SealedLeb,
                tiered: None,
            })?;
        if let Some(reader) = reader {
            let attempt = reader.utility_reader_phases(jd, site);
            match attempt {
                Ok(phases) => return Ok(phases),
                Err(error)
                    if self.ctx.mode() != ProviderMode::SealedLeb && fast_reader_miss(&error) => {}
                Err(error) => {
                    return Err(SolarDurationError::Reader {
                        error,
                        sealed: self.ctx.mode() == ProviderMode::SealedLeb,
                        tiered: reader
                            .with_serving_stack(|stack| {
                                Ok(matches!(
                                    stack.serving.provider(),
                                    crate::positions::context_serving::BoundProvider::Tiered(_)
                                ))
                            })
                            .ok(),
                    });
                }
            }
        }
        let context = self.ctx.clone().with_mode(ProviderMode::NativeJpl);
        let context = crate::positions::context_serving::native_eclipse_main_context(&context)
            .map_err(SolarDurationError::Native)?;
        let observer =
            crate::positions::context_serving::NativeFixedSolarObserver::new(&context, site)
                .map_err(SolarDurationError::Native)?;
        events::calculate_local_eclipse_phases_via(
            jd,
            &mut |ut| observer.separation(ut),
            &mut |ut| observer.sun_altaz(ut),
            &mut |ut| observer.distances(ut),
        )
        .map(|phases| phases.slots)
        .map_err(SolarDurationError::Native)
    }

    /// Retain existing reader mappings, or consume the resolved explicit file
    /// once. Reader directory/tier discovery is a separate provisioner role;
    /// absence is not presented as a fictitious reader. Forced missing readers
    /// retain the existing typed missing-data refusal before any native loader.
    fn utility_reader_owner(&self) -> Result<Option<Self>, Error> {
        if !matches!(
            self.ctx.mode(),
            ProviderMode::Auto | ProviderMode::SealedLeb
        ) {
            return Ok(None);
        }
        if !self.ctx.leb_files().is_empty() {
            return Ok(Some(self.clone()));
        }
        if let Some(path) = self.ctx.leb_file() {
            let file = self.ctx.open_leb_file(path)?;
            let mut context = self.ctx.clone();
            context.attach_leb_file(file);
            return Ok(Some(Self::from_context(context)));
        }
        if self.ctx.mode() == ProviderMode::SealedLeb {
            crate::positions::require_leb(false)?;
        }
        Ok(None)
    }

    /// Low-level fast topocentric ecliptic place at already acquired scalar
    /// model TT and the original UT1. No velocity stencil is requested.
    pub(super) fn utility_reader_place(
        &self,
        ut: f64,
        tt: f64,
        site: [f64; 3],
        body: i32,
    ) -> Result<[f64; 6], Error> {
        self.with_serving_stack(|stack| {
            crate::positions::serving_compute_tt_topocentric(
                stack.serving,
                tt,
                ut,
                body,
                FLG_TOPOCTR,
                (site[0], site[1], site[2]),
                &|tt, ut| crate::erfa::c2t06a(tt, ut, 0.0, 0.0),
                false,
                crate::positions::SiteEllipsoid::Wgs84,
            )
            .map(|(place, _)| place)
        })
    }

    /// Reader local closures each acquire their consumed model TT once. Both
    /// Sun and Moon share that scalar TT in separation/distance callbacks.
    /// The altitude callback reads only Sun, then independent frame TT and
    /// configured sidereal branches in the public horizontal primitive order.
    fn utility_reader_phases(&self, jd: f64, site: [f64; 3]) -> Result<[f64; 10], Error> {
        let pair = |ut| {
            let tt = ut + self.deltat(ut)?;
            let sun = self.utility_reader_place(ut, tt, site, SUN)?;
            let moon = self.utility_reader_place(ut, tt, site, MOON)?;
            Ok::<_, Error>((sun, moon))
        };
        events::calculate_local_eclipse_phases_via(
            jd,
            &mut |ut| {
                let (sun, moon) = pair(ut)?;
                events::vincenty_separation_deg(sun[0], sun[1], moon[0], moon[1])
            },
            &mut |ut| {
                let tt = ut + self.deltat(ut)?;
                let sun = self.utility_reader_place(ut, tt, site, SUN)?;
                fixed_reader_sun_altaz(self, ut, site, sun)
            },
            &mut |ut| {
                let (sun, moon) = pair(ut)?;
                Ok((sun[2], moon[2]))
            },
        )
        .map(|phases| phases.slots)
    }
}

/// Actual stored-file range/body misses eligible for an automatic provider retry.
fn fast_reader_miss(error: &Error) -> bool {
    matches!(
        error.primary_error(),
        Error::LebOutOfRange { .. } | Error::LebBodyNotFound { .. }
    )
}

/// Moon/Sun acquisition before any shadow reduction, including its metadata.
fn lunar_shadow_from_positions(
    jd: f64,
    position: &dyn Fn(i32) -> Result<[f64; 6], Error>,
) -> Result<(i32, LunarShadowDcore), Error> {
    let moon = position(MOON)?;
    let sun = position(SUN)?;
    let how = events::lun_how_core_via(
        (moon[0], moon[1], moon[2]),
        (sun[0], sun[1], sun[2]),
        events::saros_info(jd, events::EclipseKind::Lunar),
    )?;
    Ok((
        how.retc(),
        LunarShadowDcore::from_lun_how_dcore(how.dcore())?,
    ))
}

/// Local zero gate and source day→hour→minute rounding, without a cached event.
fn central_duration_minutes(c2: f64, c3: f64) -> f64 {
    if c2 <= 0.0 || c3 <= 0.0 {
        0.0
    } else {
        (c3 - c2) * 24.0 * 60.0
    }
}

/// Public ecliptic horizon tail with independent frame TT before the
/// configured apparent sidereal clock. Temperature0 and estimated pressure
/// are consumed even though only true altitude/south-west azimuth are returned.
fn fixed_reader_sun_altaz(
    session: &Session,
    jd: f64,
    site: [f64; 3],
    sun: [f64; 6],
) -> Result<(f64, f64), Error> {
    let frame_tt = crate::time::default_timescale_tt_from_ut1(jd)?;
    let eps0 = crate::positions::vondrak_mean_obliquity_deg(frame_tt)?;
    let (dpsi, deps) = crate::erfa::nut06a(frame_tt)?;
    let eps = eps0 + deps.to_degrees();
    let (ra, dec, _) = events::cotrans(sun[0], sun[1], sun[2], -eps);
    let hours = crate::time::public_sidereal::hours_with_longitude(
        jd,
        site[0],
        eps,
        dpsi.to_degrees(),
        &|ut| session.deltat(ut),
    )?;
    let ha = (hours * 15.0 - ra).rem_euclid(360.0).to_radians();
    let dec = dec.to_radians();
    let lat = site[1].to_radians();
    let sin_alt = lat.sin() * dec.sin() + lat.cos() * dec.cos() * ha.cos();
    let limited = if sin_alt < 1.0 { sin_alt } else { 1.0 };
    let limited = if limited > -1.0 { limited } else { -1.0 };
    let altitude = limited.asin().to_degrees();
    let y = ha.sin() * dec.cos();
    let x = ha.cos() * dec.cos() * lat.sin() - dec.sin() * lat.cos();
    let azimuth = y.atan2(x).to_degrees().rem_euclid(360.0);
    let pressure = crate::refraction::estimate_pressure_mbar(site[2]);
    let lapse = session
        .ctx
        .lapse_rate_override()
        .unwrap_or(events::RISE_SET_DIP_LAPSE_RATE_K_M);
    let _ = crate::refraction::refrac_extended_true_to_app(altitude, site[2], pressure, 0.0, lapse);
    Ok((altitude, azimuth))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::RefCell;

    #[test]
    fn lunar_later_sun_error_precedes_checked_moon_geometry() {
        let trace = RefCell::new(Vec::new());
        let error = lunar_shadow_from_positions(2451545.0, &|body| {
            trace.borrow_mut().push(body);
            if body == MOON {
                Ok([f64::NAN; 6])
            } else {
                Err(Error::UnknownBody { body_id: 719 })
            }
        })
        .unwrap_err();
        assert!(matches!(error, Error::UnknownBody { body_id: 719 }));
        assert_eq!(*trace.borrow(), [MOON, SUN]);
    }

    #[test]
    fn duration_uses_both_original_multiplications_and_nonpositive_gate() {
        assert_eq!(central_duration_minutes(0.0, 3.0), 0.0);
        assert_eq!(central_duration_minutes(3.0, -0.0), 0.0);
        let (c2, c3) = (2460409.76234, 2460409.76351);
        assert_eq!(
            central_duration_minutes(c2, c3).to_bits(),
            ((c3 - c2) * 24.0 * 60.0).to_bits()
        );
    }

    #[test]
    fn only_actual_reader_misses_are_retryable() {
        assert!(fast_reader_miss(&Error::LebBodyNotFound { body_id: 1 }));
        assert!(!fast_reader_miss(&Error::LebCorrupt {
            what: "directory".into(),
            offset: 1
        }));
        assert!(!fast_reader_miss(&Error::UnknownBody { body_id: 1 }));
    }
}
