//! Fixed-time solar details, raw-bracket maxima and planetary occultations.
//! Stored scalar/model acquisition and native observed-vector acquisition are
//! distinct. Each operation retains its immutable image and explicit site.
use super::*;
use crate::events::{self as events, EclipseHowDetailsEphemeris};
use crate::positions::context_serving::NativeFixedSolarObserver;

impl Session {
    /// Model TT Julian day from an original UT1 sample and this context's Delta-T.
    pub(super) fn event_model_tt(&self, ut: f64) -> Result<f64, Error> {
        self.deltat(ut).map(|delta| ut + delta)
    }

    /// Topocentric ecliptic place over retained strict files at explicit clocks.
    pub(super) fn event_place_at(
        &self,
        ut: f64,
        tt: f64,
        site: [f64; 3],
        body: i32,
    ) -> Result<[f64; 6], Error> {
        self.utility_reader_place(ut, tt, site, body)
    }

    /// Public native position, optionally with an owned temporary observer.
    pub(super) fn event_public_position(
        &self,
        ut: f64,
        body: i32,
        flags: i32,
        site: Option<[f64; 3]>,
    ) -> Result<[f64; 6], Error> {
        match site {
            None => self.calc_ut(ut, body, flags).map(|result| result.values),
            Some(site) => {
                let context = self
                    .ctx
                    .clone()
                    .with_observer(Observer::new(site[0], site[1], site[2])?);
                Session::from_context(context)
                    .calc_ut(ut, body, flags)
                    .map(|result| result.values)
            }
        }
    }

    /// Global planetary occultation search at UT1 Julian day.
    /// Negative direction searches backward; otherwise forward. A nonempty
    /// star name selects the immutable stellar catalog when no target id is given.
    /// Returns character flags and ten contact Julian-day slots.
    pub fn planet_occult_when_glob(
        &self,
        jd_ut: f64,
        occulting: i32,
        occulted: Option<i32>,
        starname: &str,
        direction: i32,
    ) -> Result<(i32, [f64; 10]), Error> {
        events::validate_planet_occult_request(occulting, occulted, starname)?;
        let provider = PlanetaryPlaces {
            session: self,
            occulting,
            occulted,
            starname,
            site: None,
            native: None,
        };
        events::planet_occult_when_glob_via(
            &provider, jd_ut, occulting, occulted, starname, direction,
        )
    }

    /// Find a planetary occultation visible from an explicit WGS84 site.
    /// Global candidates retain their geocentric scalar channel; local places
    /// retain stored model TT or native bundled Time. All eleven window samples
    /// and maximum are consumed, then occulting altitude, target altitude and
    /// local geometry. Both altitudes must exceed -0.5°. Return character bits,
    /// ten locally ordered UT1 contact slots and twenty attributes, including
    /// original reserved zeros. Only public Error/RuntimeError ends a candidate;
    /// arithmetic/lookup/value failures propagate. No time limit means no event.
    pub fn planet_occult_when_loc(
        &self,
        jd_ut: f64,
        occulting: i32,
        occulted: Option<i32>,
        starname: &str,
        site: [f64; 3],
    ) -> Result<(i32, [f64; 10], [f64; 20]), Error> {
        let result =
            self.planet_occult_when_loc_result(jd_ut, occulting, occulted, starname, site)?;
        Ok((result.flags, result.times, result.attributes))
    }

    /// Local occultation as an owned record of character, contacts and geometry.
    /// Images and the WGS84 observer survive the entire native solve; resource
    /// or numerical failures propagate without fabricating an event.
    pub fn planet_occult_when_loc_result(
        &self,
        jd_ut: f64,
        occulting: i32,
        occulted: Option<i32>,
        starname: &str,
        site: [f64; 3],
    ) -> Result<events::PlanetOccultLocalResult, Error> {
        events::validate_planet_occult_request(occulting, occulted, starname)?;
        let native = if self.ctx.mode() == ProviderMode::NativeJpl {
            let context =
                crate::positions::context_serving::native_eclipse_main_context(&self.ctx)?;
            Some(crate::positions::context_serving::NativePlanetObserver::new(&context, site)?)
        } else {
            None
        };
        let provider = PlanetaryPlaces {
            session: self,
            occulting,
            occulted,
            starname,
            site: Some(site),
            native,
        };
        events::planet_occult_when_loc_result_via(
            &provider, jd_ut, occulting, occulted, starname, site[1], site[0],
        )
    }

    /// Minimize an acquired global gamma objective over a UT1-JD bracket.
    ///
    /// `gamma` returns the shadow-axis distance in Earth radii at each consumed
    /// UT1 sample. It owns acquisition and fallback policy, including any live
    /// configuration. This search acquires no provider or model clock. It uses
    /// the original 60-step golden loop, stops after an update when the bracket
    /// width is below 1e-8 day, then independently rereads the returned midpoint.
    /// Zero, reversed and IEEE nonfinite bounds retain their actual samples.
    /// The first objective error terminates the search without a later read.
    pub fn sol_eclipse_max_time_with_gamma(
        low: f64,
        high: f64,
        mut gamma: &mut dyn FnMut(f64) -> Result<f64, Error>,
    ) -> Result<(f64, f64), Error> {
        minimum_from_bracket(&mut gamma, low, high)
    }

    /// Read apparent geocentric Sun then Moon J2000 reader vectors in AU.
    ///
    /// `tt` is the already acquired scalar model TT Julian day. This method
    /// neither evaluates Delta-T nor applies an outer fallback policy. The
    /// strict immutable files remain owned throughout the coefficient reads.
    /// Coverage, body and data failures propagate as native typed errors.
    pub fn solar_gamma_reader_vectors_at_tt(&self, tt: f64) -> Result<([f64; 3], [f64; 3]), Error> {
        self.with_serving_stack(|stack| {
            events::solar_gamma_reader_vectors_at_tt(stack.serving, tt)
                .map(|(sun, moon)| ([sun.0, sun.1, sun.2], [moon.0, moon.1, moon.2]))
        })
    }

    /// Acquire apparent geocentric Sun then Moon ICRS rays in AU at UT1 JD.
    ///
    /// The NativeJpl context owns Main selection, independent bundled split
    /// Time and Earth state before the two light-time observations. Configured
    /// model Delta-T, sidereal/output flags and surface vectors are not consumed.
    /// Provider directory, coverage, acquisition and observation failures
    /// propagate at their actual stage. There is no reader retry or gamma
    /// geometry inside this method.
    pub fn solar_gamma_native_vectors_at_ut(&self, ut: f64) -> Result<([f64; 3], [f64; 3]), Error> {
        crate::positions::context_serving::native_solar_gamma_pair(&self.ctx, ut)
    }

    /// Minimize geocentric gamma or local Sun-Moon separation over a UT1 bracket.
    /// The 60-step golden solve uses a 1e-8-day width target. A site selects
    /// angular separation in degrees; otherwise the objective is Earth radii.
    pub fn sol_eclipse_max_time_from_bracket(
        &self,
        low: f64,
        high: f64,
        site: Option<[f64; 3]>,
    ) -> Result<(f64, f64), Error> {
        if let Some(site) = site {
            if self.ctx.mode() != ProviderMode::NativeJpl {
                let stored = StoredDetails {
                    session: self,
                    site,
                    flags: 0,
                };
                return minimum_from_bracket(
                    &mut |ut| stored.details_separation_deg(ut),
                    low,
                    high,
                );
            }
            let context =
                crate::positions::context_serving::native_eclipse_main_context(&self.ctx)?;
            let observer = NativeFixedSolarObserver::new(&context, site)?;
            return minimum_from_bracket(&mut |ut| observer.separation(ut), low, high);
        }
        if self.ctx.mode() != ProviderMode::NativeJpl {
            return self.with_serving_stack(|stack| {
                minimum_from_bracket(
                    &mut |ut| {
                        let tt = ut + crate::time::deltat_with_config(ut, stack.delta_cfg)?;
                        let (sun, moon) =
                            events::solar_gamma_reader_vectors_at_tt(stack.serving, tt)?;
                        events::gamma_from_vectors(sun, moon)
                    },
                    low,
                    high,
                )
            });
        }
        let context = crate::positions::context_serving::native_eclipse_main_context(&self.ctx)?;
        minimum_from_bracket(
            &mut |ut| {
                let (sun, moon) =
                    crate::positions::context_serving::native_solar_gamma_pair(&context, ut)?;
                events::gamma_from_vectors((sun[0], sun[1], sun[2]), (moon[0], moon[1], moon[2]))
            },
            low,
            high,
        )
    }

    /// Solar maximum near an approximate UT1 JD. Site and result units match
    /// [`Self::sol_eclipse_max_time_from_bracket`]. Subtraction then addition
    /// construct the bracket; no extra finite or positive-range check is used.
    pub fn sol_eclipse_max_time(
        &self,
        jd_approx: f64,
        site: Option<[f64; 3]>,
        search_range: f64,
    ) -> Result<(f64, f64), Error> {
        let low = jd_approx - search_range;
        let high = jd_approx + search_range;
        self.sol_eclipse_max_time_from_bracket(low, high, site)
    }

    /// Comprehensive local solar circumstances at UT1 JD: type/predicates,
    /// instant and maximum magnitude/obscuration, radii and separation degrees,
    /// contact UT1 JDs, true Sun altitudes/south-westward azimuths, ICRS position
    /// angles and phase durations in minutes. Missing contacts are zero.
    /// Caller flags affect only the initial public how; subsequent local
    /// observations use their independent acquisition channels. The existing
    /// contact isolation and search tolerances are retained. Data/time/frame
    /// errors propagate except at the explicitly protected contact blocks.
    pub fn sol_eclipse_how_details(
        &self,
        jd_ut: f64,
        site: [f64; 3],
        flags: i32,
    ) -> Result<events::EclipseHowDetails, Error> {
        if self.ctx.mode() != ProviderMode::NativeJpl {
            return events::eclipse_how_details_via(
                &StoredDetails {
                    session: self,
                    site,
                    flags,
                },
                jd_ut,
                site[0],
                site[1],
                site[2],
            );
        }
        let context = crate::positions::context_serving::native_eclipse_main_context(&self.ctx)?;
        let observer = NativeFixedSolarObserver::new(&context, site)?;
        events::eclipse_how_details_via(
            &NativeDetails {
                session: self,
                site,
                flags,
                observer,
            },
            jd_ut,
            site[0],
            site[1],
            site[2],
        )
    }

    /// Instant visible solar magnitude, as the covered fraction of the Sun's
    /// diameter. Site units match maximum above. True Sun altitude below -1°
    /// and nonoverlapping discs return zero; a total eclipse may exceed one.
    /// Caller flags are not consumed. Stored places share one configured model
    /// TT; native rays use independent bundled UT1 Time and a real WGS84 site.
    pub fn sol_eclipse_magnitude_at_loc(&self, jd_ut: f64, site: [f64; 3]) -> Result<f64, Error> {
        self.solar_instant_fraction(jd_ut, site, false)
    }

    /// Instant visible fraction of the Sun's disc area covered by the Moon.
    /// Visibility, observations, units and errors match magnitude above. Full
    /// coverage is one; annular and partial overlap use the original circle
    /// intersection arithmetic with its unchanged ordered clamps.
    pub fn sol_eclipse_obscuration_at_loc(&self, jd_ut: f64, site: [f64; 3]) -> Result<f64, Error> {
        self.solar_instant_fraction(jd_ut, site, true)
    }

    /// Protect only the observed-place acquisition. Horizon and overlap
    /// arithmetic follow afterward and retain their errors; range misses are
    /// never turned into a zero observation.
    fn solar_instant_fraction(
        &self,
        jd: f64,
        site: [f64; 3],
        obscuration: bool,
    ) -> Result<f64, Error> {
        let sample = if self.ctx.mode() != ProviderMode::NativeJpl {
            let pair = (|| {
                let tt = self.event_model_tt(jd)?;
                let sun = self.event_place_at(jd, tt, site, SUN)?;
                let moon = self.event_place_at(jd, tt, site, MOON)?;
                Ok::<_, Error>((sun, moon))
            })();
            let (sun, moon) = match pair {
                Ok(pair) => pair,
                Err(error) => {
                    let caught = instant_observation_caught(&error);
                    if caught {
                        return Ok(0.0);
                    }
                    return Err(error);
                }
            };
            let (alt, _) = stored_altaz(self, jd, site, sun)?;
            if alt < -1.0 {
                return Ok(0.0);
            }
            (
                alt,
                events::vincenty_separation_deg(sun[0], sun[1], moon[0], moon[1])?,
                sun[2],
                moon[2],
            )
        } else {
            let context =
                crate::positions::context_serving::native_eclipse_main_context(&self.ctx)?;
            let observer = NativeFixedSolarObserver::new(&context, site)?;
            match observer.instant_geometry(jd, &instant_observation_caught)? {
                Some(sample) => sample,
                None => return Ok(0.0),
            }
        };
        if obscuration {
            events::sol_eclipse_obscuration_at_loc_via(sample.0, sample.1, sample.2, sample.3)
        } else {
            events::sol_eclipse_magnitude_at_loc_via(sample.0, sample.1, sample.2, sample.3)
        }
    }
}

/// Original 60-step golden loop and final independent objective read. Checked
/// general minimizers are separate contracts and retain their admission guards.
fn minimum_from_bracket(
    objective: &mut impl FnMut(f64) -> Result<f64, Error>,
    low: f64,
    high: f64,
) -> Result<(f64, f64), Error> {
    let jd = events::golden_section_min_source(objective, low, high, 1e-8, 60, &mut Vec::new())?;
    Ok((jd, objective(jd)?))
}

/// Numerical observation-domain failures handled as no visible geometry.
/// File/body/coverage failures propagate for provider selection by the caller.
fn instant_observation_caught(error: &Error) -> bool {
    if matches!(
        error.primary_error(),
        Error::LebOutOfRange { .. } | Error::EphemerisRange { .. }
    ) {
        return false;
    }

    matches!(
        error.primary_error(),
        Error::MathDomain { .. }
            | Error::NonFiniteArgument { .. }
            | Error::IntegerOverflow
            | Error::RationalDivisionByZero
            | Error::NumericalOverflow { .. }
    )
}

/// Stored reader closures retain their original one-model-TT pair order and
/// public horizontal frame clock, without a continuity-offset preflight.
struct StoredDetails<'a> {
    session: &'a Session,
    site: [f64; 3],
    flags: i32,
}
impl StoredDetails<'_> {
    /// Acquire Sun then Moon at the same already selected model TT.
    fn pair(&self, jd: f64) -> Result<([f64; 6], [f64; 6]), Error> {
        let tt = self.session.event_model_tt(jd)?;
        let sun = self.session.event_place_at(jd, tt, self.site, SUN)?;
        let moon = self.session.event_place_at(jd, tt, self.site, MOON)?;
        Ok((sun, moon))
    }
}
impl EclipseHowDetailsEphemeris for StoredDetails<'_> {
    fn details_initial(&self, jd: f64) -> Result<(i32, [f64; 20]), Error> {
        self.session.sol_eclipse_how(jd, self.site, self.flags)
    }
    fn details_separation_deg(&self, jd: f64) -> Result<f64, Error> {
        let (s, m) = self.pair(jd)?;
        events::vincenty_separation_deg(s[0], s[1], m[0], m[1])
    }
    fn details_sun_alt_az_deg(&self, jd: f64) -> Result<(f64, f64), Error> {
        let tt = self.session.event_model_tt(jd)?;
        let sun = self.session.event_place_at(jd, tt, self.site, SUN)?;
        stored_altaz(self.session, jd, self.site, sun)
    }
    fn details_disc_radii_deg(&self, jd: f64) -> Result<(f64, f64), Error> {
        let (s, m) = self.pair(jd)?;
        events::eclipse_disc_angular_radii_deg(s[2], m[2])
    }
    fn details_equatorial_deg(&self, jd: f64) -> Result<(f64, f64, f64, f64), Error> {
        let mut ctx = self.session.ctx.clone();
        ctx.set_observer(Observer::new(self.site[0], self.site[1], self.site[2])?);
        let owner = Session::from_context(ctx);
        let flags = FLG_EQUATORIAL | FLG_J2000 | FLG_TOPOCTR | FLG_SPEED;
        let s = owner.calc_ut(jd, SUN, flags)?.values;
        let m = owner.calc_ut(jd, MOON, flags)?.values;
        Ok((s[0], s[1], m[0], m[1]))
    }
}

/// Native details retain one real observer; every closure constructs a fresh
/// independent Time and observes the bodies it actually consumes.
struct NativeDetails<'a> {
    session: &'a Session,
    site: [f64; 3],
    flags: i32,
    observer: NativeFixedSolarObserver,
}
impl EclipseHowDetailsEphemeris for NativeDetails<'_> {
    fn details_initial(&self, jd: f64) -> Result<(i32, [f64; 20]), Error> {
        self.session.sol_eclipse_how(jd, self.site, self.flags)
    }
    fn details_separation_deg(&self, jd: f64) -> Result<f64, Error> {
        self.observer.separation(jd)
    }
    fn details_sun_alt_az_deg(&self, jd: f64) -> Result<(f64, f64), Error> {
        self.observer.pair_sun_altaz(jd)
    }
    fn details_disc_radii_deg(&self, jd: f64) -> Result<(f64, f64), Error> {
        let (s, m) = self.observer.distances(jd)?;
        events::eclipse_disc_angular_radii_deg(s, m)
    }
    fn details_equatorial_deg(&self, jd: f64) -> Result<(f64, f64, f64, f64), Error> {
        self.observer.equatorial(jd)
    }
    fn details_position_angle_deg(&self, jd: f64) -> Result<f64, Error> {
        self.observer.position_angle(jd)
    }
}

/// Genuine ECL2HOR projection after the stored place. Frame TT remains
/// independent from model TT; the result is true altitude then azimuth.
fn stored_altaz(
    session: &Session,
    jd: f64,
    site: [f64; 3],
    place: [f64; 6],
) -> Result<(f64, f64), Error> {
    let (az, alt, _) = super::solar_circumstances::solar_horizontal(
        session,
        jd,
        site,
        [place[0], place[1], place[2]],
    )?;
    Ok((alt, az))
}

/// Separate global public positions and local observed-place channels. An
/// explicit native observer is absent for a global query and retained once
/// for the local walk; no configuration is borrowed.
struct PlanetaryPlaces<'a> {
    session: &'a Session,
    occulting: i32,
    occulted: Option<i32>,
    starname: &'a str,
    site: Option<[f64; 3]>,
    native: Option<crate::positions::context_serving::NativePlanetObserver>,
}
impl PlanetaryPlaces<'_> {
    /// Resolve a star at its actual target-read boundary, preserving the
    /// documented ValueError for a name miss and its original search string.
    fn star_data(&self) -> Result<crate::stars::StarData, Error> {
        use crate::stars::StarResolveError;
        let found =
            crate::stars::resolve_star_ref(self.starname).map_err(|failure| match failure {
                StarResolveError::Miss(failure) => Error::InvalidEclipseSearch {
                    entry: "planet_occult",
                    detail: failure.message,
                },
                StarResolveError::Refused(error) => error,
            })?;
        crate::stars::fixed_star_data(found.star_id).ok_or(Error::UnknownBody {
            body_id: found.star_id,
        })
    }
    /// Linear catalog proper motion in degrees at the requested UT1 epoch.
    /// The observer applies apparent-place corrections to this stellar direction.
    fn local_star(
        &self,
        jd: f64,
    ) -> Result<crate::positions::context_serving::PlanetObservedTarget, Error> {
        let data = self.star_data()?;
        let years = (jd - 2451545.0) / 365.25;
        Ok(
            crate::positions::context_serving::PlanetObservedTarget::Star {
                ra: data.ra_j2000_deg() + (data.pm_ra_arcsec_yr() * years) / 3600.0,
                dec: data.dec_j2000_deg() + (data.pm_dec_arcsec_yr() * years) / 3600.0,
            },
        )
    }
    /// Fast topocentric SPEED position at one configured scalar model TT.
    /// Caller output flags are not mixed into this independent local request.
    fn stored_body(&self, jd: f64, body: i32) -> Result<[f64; 6], Error> {
        let tt = self.session.event_model_tt(jd)?;
        self.stored_body_at_tt(jd, tt, body)
    }
    /// A held scalar model TT shared by the maximum's two reader positions.
    /// This entry does not resample Delta-T between occulting and target.
    fn stored_body_at_tt(&self, jd: f64, tt: f64, body: i32) -> Result<[f64; 6], Error> {
        let site = self.site.expect("local provider owns site");

        self.session.with_serving_stack(|stack| {
            crate::positions::serving_compute_tt_topocentric(
                stack.serving,
                tt,
                jd,
                body,
                FLG_TOPOCTR | FLG_SPEED,
                (site[0], site[1], site[2]),
                &|tt, ut| crate::erfa::c2t06a(tt, ut, 0., 0.),
                false,
                crate::positions::SiteEllipsoid::Wgs84,
            )
            .map(|value| value.0)
        })
    }
}
impl events::PlanetOccultEphemeris for PlanetaryPlaces<'_> {
    fn occulting_position(&self, jd: f64) -> Result<(f64, f64, f64), Error> {
        let p = self.session.event_public_position(
            jd,
            self.occulting,
            FLG_EQUATORIAL | FLG_SPEED,
            None,
        )?;
        Ok((
            p[0],
            p[1],
            events::calc_planet_angular_radius_deg(self.occulting, p[2]),
        ))
    }
    fn target_position(&self, jd: f64) -> Result<(f64, f64, f64), Error> {
        if let Some(body) = self.occulted {
            let p =
                self.session
                    .event_public_position(jd, body, FLG_EQUATORIAL | FLG_SPEED, None)?;
            Ok((
                p[0],
                p[1],
                events::calc_planet_angular_radius_deg(body, p[2]),
            ))
        } else {
            self.star_data()?;
            let p = self
                .session
                .fixstar_ut(self.starname, jd, FLG_EQUATORIAL)?
                .values;
            Ok((p[0], p[1], events::OCC_POINT_SOURCE_RADIUS_DEG))
        }
    }
    fn occulting_ra_dec_deg(&self, jd: f64) -> Result<(f64, f64), Error> {
        let (a, b, _) = self.occulting_position(jd)?;
        Ok((a, b))
    }
    fn target_ra_dec_deg(&self, jd: f64) -> Result<(f64, f64), Error> {
        let (a, b, _) = self.target_position(jd)?;
        Ok((a, b))
    }
    fn occulting_radius_deg(&self, jd: f64) -> Result<f64, Error> {
        Ok(self.occulting_position(jd)?.2)
    }
    fn target_radius_deg(&self, jd: f64) -> Result<f64, Error> {
        Ok(self.target_position(jd)?.2)
    }
}
impl events::PlanetOccultLocalEphemeris for PlanetaryPlaces<'_> {
    fn occulting_alt_az_deg(&self, jd: f64) -> Result<(f64, f64), Error> {
        if let Some(observer) = &self.native {
            observer
                .place(
                    jd,
                    crate::positions::context_serving::PlanetObservedTarget::Body(self.occulting),
                )?
                .altaz()
        } else {
            stored_altaz(
                self.session,
                jd,
                self.site.unwrap(),
                self.stored_body(jd, self.occulting)?,
            )
        }
    }
    fn target_alt_az_deg(&self, jd: f64) -> Result<(f64, f64), Error> {
        if let Some(observer) = &self.native {
            let (epochs, tt) =
                crate::positions::reference_frames::VectorRotationEpochs::from_ut1(jd)?;
            let target = match self.occulted {
                Some(body) => crate::positions::context_serving::PlanetObservedTarget::Body(body),
                None => self.local_star(jd)?,
            };
            observer.place_after_time(epochs, tt, target)?.altaz()
        } else {
            let p = match self.occulted {
                Some(body) => self.stored_body(jd, body)?,
                None => {
                    let _unused_tt = self.session.event_model_tt(jd)?;
                    self.session
                        .fixstar_ut(self.starname, jd, FLG_SPEED)?
                        .values
                }
            };
            stored_altaz(self.session, jd, self.site.unwrap(), p)
        }
    }
    fn local_radii_sep_deg(&self, jd: f64) -> Result<(f64, f64, f64), Error> {
        if let Some(observer) = &self.native {
            use crate::positions::context_serving::{
                PlanetObservedTarget, fixed_solar_vector_separation,
            };
            let (epochs, tt) =
                crate::positions::reference_frames::VectorRotationEpochs::from_ut1(jd)?;
            let occ = observer.place_after_time(
                epochs,
                tt,
                PlanetObservedTarget::Body(self.occulting),
            )?;
            let target = match self.occulted {
                Some(body) => PlanetObservedTarget::Body(body),
                None => self.local_star(jd)?,
            };
            let tgt = observer.place_after_time(epochs, tt, target)?;
            let separation = fixed_solar_vector_separation(occ.ray, tgt.ray);
            let occ_radius = events::calc_planet_angular_radius_deg(self.occulting, occ.distance);
            let tgt_radius = self
                .occulted
                .map_or(events::OCC_POINT_SOURCE_RADIUS_DEG, |body| {
                    events::calc_planet_angular_radius_deg(body, tgt.distance)
                });
            Ok((occ_radius, tgt_radius, separation))
        } else {
            let tt = self.session.event_model_tt(jd)?;
            let occ = self.stored_body_at_tt(jd, tt, self.occulting)?;
            let (tgt, tgt_radius) = match self.occulted {
                Some(body) => {
                    let p = self.stored_body_at_tt(jd, tt, body)?;
                    (p, events::calc_planet_angular_radius_deg(body, p[2]))
                }
                None => (
                    self.session
                        .fixstar_ut(self.starname, jd, FLG_SPEED)?
                        .values,
                    events::OCC_POINT_SOURCE_RADIUS_DEG,
                ),
            };
            let separation = events::vincenty_separation_deg(occ[0], occ[1], tgt[0], tgt[1])?;
            let occ_radius = events::calc_planet_angular_radius_deg(self.occulting, occ[2]);
            Ok((occ_radius, tgt_radius, separation))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn raw_brackets_acquire_before_any_window_refusal() {
        for (low, high) in [(1., 1.), (2., 1.), (f64::NAN, 1.)] {
            let mut reads = 0;
            let error = minimum_from_bracket(
                &mut |_| {
                    reads += 1;
                    Err(Error::UnknownBody { body_id: 71 })
                },
                low,
                high,
            )
            .unwrap_err();
            assert!(matches!(error, Error::UnknownBody { body_id: 71 }));
            assert_eq!(reads, 1);
        }
    }
    #[test]
    fn zero_window_still_updates_then_rereads() {
        let mut reads = Vec::new();
        let result = minimum_from_bracket(
            &mut |jd| {
                reads.push(jd);
                Ok(0.)
            },
            2451545.,
            2451545.,
        )
        .unwrap();
        assert_eq!(result, (2451545., 0.));
        assert_eq!(reads, vec![2451545.; 4]);
    }
    #[test]
    fn injected_gamma_rereads_the_final_sample_with_current_owner_state() {
        let mut reads = 0;
        let result = Session::sol_eclipse_max_time_with_gamma(2451545., 2451545., &mut |_| {
            reads += 1;
            Ok(reads as f64)
        })
        .unwrap();
        assert_eq!(reads, 4);
        assert_eq!(result, (2451545., 4.));
    }
    #[test]
    fn injected_gamma_retains_the_first_error_and_stops_acquisition() {
        let mut reads = 0;
        let error = Session::sol_eclipse_max_time_with_gamma(2., 1., &mut |_| {
            reads += 1;
            Err(Error::UnknownBody { body_id: 71 })
        })
        .unwrap_err();
        assert_eq!(reads, 1);
        assert!(matches!(error, Error::UnknownBody { body_id: 71 }));
    }
}
