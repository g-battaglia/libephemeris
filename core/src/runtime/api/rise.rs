//! Context-owned planetary rise, set and meridian-transit searches.
use super::*;
use crate::events::{self, RiseSetEvent, RiseSetOutcome, RiseSetRequest};

mod search;
use search::native_rise_set_search_with_policy;

/// Synchronous configuration-only resolver borrowed for one event search.
/// The context is shared and the returned snapshot owns its configuration and
/// registry; the callback retains no input reference or astronomical values.
type PositionContextResolver<'a> =
    dyn Fn(&EphemerisContext, i32) -> Result<EphemerisContext, Error> + 'a;

/// Boundary policy for latitude, which meridian transit never consumes.
/// Checked callers retain all-scalar finite admission; operand-compatible
/// callers retain the converted latitude until an actual rise/set consumes it.
#[derive(Clone, Copy)]
enum EventInputPolicy {
    Checked,
    ConsumedLatitude,
}

impl Session {
    /// Next planetary or hypothetical-model rise, set or transit at a UT1
    /// Julian day, using a zero horizon. `geopos` is east longitude/latitude
    /// in degrees and altitude in metres; pressure is mbar (zero estimates
    /// it), temperature is Celsius.
    /// `rsmi` selects SET > RISE > ITRANSIT > MTRANSIT, or rise if absent, plus
    /// disc, refraction, twilight and Hindu projection modifiers.
    /// Returns `(0, tret)` for an event, or `(-2, [0;10])` for no crossing;
    /// `tret[0]` is UT1 JD and the remaining nine words are positive zero.
    /// See [`Self::rise_trans_true_hor`] for clocks, errors and search limits.
    #[allow(clippy::too_many_arguments)]
    pub fn rise_trans(
        &self,
        jd_ut: f64,
        body: i32,
        rsmi: i32,
        geopos: [f64; 3],
        pressure: f64,
        temperature: f64,
        flags: i32,
    ) -> Result<(i32, [f64; 10]), Error> {
        self.rise_trans_true_hor(jd_ut, body, rsmi, geopos, pressure, temperature, 0.0, flags)
    }

    /// Planetary horizon event at UT1 with `horhgt` in degrees. Negative
    /// horizons are admitted; `-100` means the refracted sea-horizon dip plus
    /// 1e-4 degree. Sun twilight selects -6/-12/-18 degrees and forces its
    /// geometric centre. Hindu rise/set uses the geocentric ecliptic place
    /// with latitude zeroed; meridian transit ignores that modifier.
    ///
    /// NativeJpl samples the genuine public [`Self::calc_ut`] route: apparent
    /// spherical TOPOCTR|SPEED for altitude, EQUATORIAL|SPEED for transit and
    /// ordinary circumpolar declination, SPEED for Hindu projection. These
    /// fixed callback flags do not consume caller `flags`; the one-try bit is
    /// inert. Scalar public UT1, split provider samples, derivative scalar-TT
    /// sides and independent Topos clocks remain owned by that route.
    /// [`Self::azalt`] supplies the shared spherical ECL2HOR projection:
    /// configured Delta-T supplies its sidereal terms, bundled global TT its
    /// frame. Meridian GAST instead uses fresh bundled UT1 Time parts and
    /// an explicit original UT1 fraction, independently of the configured model.
    /// No native vector-altaz substitute is used. Explicit LEB/Auto retains
    /// its sealed planetary or analytical-model event serving and its
    /// calculation-flag gate.
    ///
    /// The ordinary rise/set search owns a temporary observer context; Hindu
    /// and transit never validate or replace an unused observer. The caller's
    /// context/settings/registry stay unchanged. Out-of-map integers undergo
    /// the preliminary public SPEED probe: its typed range error propagates,
    /// other positionability failures carry the requested UnknownBody key.
    /// Earth has no apparent direction; ECL_NUT first undergoes that probe.
    /// Neither becomes an invented event. Star names require a separate API.
    ///
    /// Rise/set scans two days with grazing-aware 10/20/60-minute steps and
    /// at most 60 bisections, the existing 0.03-second slope-aware exit and
    /// conditional Newton refinement. Transit uses at most 30 hour-angle
    /// iterations and returns its final estimate when not converged. Data,
    /// provider, support, original sample range and observer errors propagate;
    /// public nonfinite scalar inputs refuse. These bounded searches establish
    /// no certified physical accuracy or wider integer-family acceptance.
    #[allow(clippy::too_many_arguments)]
    pub fn rise_trans_true_hor(
        &self,
        jd_ut: f64,
        body: i32,
        rsmi: i32,
        geopos: [f64; 3],
        pressure: f64,
        temperature: f64,
        horhgt: f64,
        flags: i32,
    ) -> Result<(i32, [f64; 10]), Error> {
        self.rise_trans_true_hor_inner(
            jd_ut,
            body,
            rsmi,
            geopos,
            pressure,
            temperature,
            horhgt,
            flags,
            None,
            EventInputPolicy::Checked,
        )
    }

    /// Horizon event with latitude admission at its consuming event branch.
    ///
    /// Units, flags, search bounds, output layout, ownership and acquisition
    /// follow [`Self::rise_trans_true_hor`]. Rise/set retains its finite
    /// latitude check. Meridian transit preserves an unused latitude,
    /// including NaN or infinity, without constructing an observer or replacing
    /// that value. Every other scalar retains its finite check. Pressure and
    /// sea-horizon setup still run in their ordinary order before the transit
    /// search; no altitude, refraction or latitude callback is consumed there.
    /// This boundary serves callers that have already converted geographic
    /// operands; the checked Rust and C entry points remain unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn rise_trans_true_hor_consumed(
        &self,
        jd_ut: f64,
        body: i32,
        rsmi: i32,
        geopos: [f64; 3],
        pressure: f64,
        temperature: f64,
        horhgt: f64,
        flags: i32,
    ) -> Result<(i32, [f64; 10]), Error> {
        self.rise_trans_true_hor_inner(
            jd_ut,
            body,
            rsmi,
            geopos,
            pressure,
            temperature,
            horhgt,
            flags,
            None,
            EventInputPolicy::ConsumedLatitude,
        )
    }

    /// Horizon event with a caller-owned NativeJpl configuration resolver.
    ///
    /// Arguments, output layout, units, flags and search limits are identical
    /// to [`Self::rise_trans_true_hor`]. This configuration seam supplies the
    /// selected main-image path at a binding boundary that owns process
    /// settings; the core itself reads no environment or home directory.
    /// The resolver is consumed only immediately before a genuine public
    /// [`Self::calc_ut`] acquisition. It receives that call's context (including
    /// the temporary observer when consumed) and actual requested body key,
    /// and returns an owned configuration snapshot preserving the provider,
    /// registry, clocks and settings except for its resolved main-image inputs.
    /// An already supplied registry must remain authoritative.
    ///
    /// Nonfinite arguments, undefined Earth direction and invalid consumed
    /// observer coordinates fail before this resolver. An out-of-map body
    /// consumes it only at its preliminary positionability probe. Non-Native
    /// modes retain their existing serving and do not consume the resolver.
    /// Resolver errors propagate at that acquisition boundary before position
    /// clocks/data; source position, frame and independent Topos clocks still
    /// belong to `calc_ut`. The callback is synchronous, retains no borrowed
    /// reference, and neither changes the caller's session nor transfers any
    /// registry/image ownership. It does not provide astronomical values.
    #[allow(clippy::too_many_arguments)]
    pub fn rise_trans_true_hor_with_position_context(
        &self,
        jd_ut: f64,
        body: i32,
        rsmi: i32,
        geopos: [f64; 3],
        pressure: f64,
        temperature: f64,
        horhgt: f64,
        flags: i32,
        resolve: &dyn Fn(&EphemerisContext, i32) -> Result<EphemerisContext, Error>,
    ) -> Result<(i32, [f64; 10]), Error> {
        self.rise_trans_true_hor_inner(
            jd_ut,
            body,
            rsmi,
            geopos,
            pressure,
            temperature,
            horhgt,
            flags,
            Some(resolve),
            EventInputPolicy::Checked,
        )
    }

    /// Validate the public event request and compose its owned callback search.
    /// The optional resolver changes configuration acquisition only; the actual
    /// position call remains public `calc_ut` at each original UT1 sample.
    /// Default Session calls retain their direct borrowed-context route.
    /// Ordinary observer construction and direction-probe catch order precede
    /// the same callbacks as without a resolver. Units, finite gates, search
    /// bounds, event/no-event layout and errors follow the public method.
    /// All callback borrows end before this method returns.
    #[allow(clippy::too_many_arguments)]
    fn rise_trans_true_hor_inner(
        &self,
        jd_ut: f64,
        body: i32,
        rsmi: i32,
        geopos: [f64; 3],
        pressure: f64,
        temperature: f64,
        horhgt: f64,
        flags: i32,
        resolve: Option<&PositionContextResolver<'_>>,
        admission: EventInputPolicy,
    ) -> Result<(i32, [f64; 10]), Error> {
        let req = RiseSetRequest {
            jd_start_ut: jd_ut,
            lon_deg_east: geopos[0],
            lat_deg: geopos[1],
            height_m: geopos[2],
            atpress_mbar: pressure,
            attemp_c: temperature,
            horhgt_deg: horhgt,
            rsmi,
        };
        match admission {
            EventInputPolicy::Checked => validate_event_inputs(&req)?,
            EventInputPolicy::ConsumedLatitude => validate_consumed_event_inputs(&req)?,
        }
        let flags = events::strip_one_try_bit(flags);
        if self.ctx.mode() != ProviderMode::NativeJpl {
            return self.with_serving_stack(|stack| {
                let dt = |ut| crate::time::deltat_with_config(ut, stack.delta_cfg);
                if events::is_fictitious_event_body(body) {
                    // These analytical models consume Sun/Earth support, not
                    // a stored channel for the requested hypothetical body.
                    let serving = events::SealedFictitiousEventServing::try_new(
                        stack.serving,
                        &dt,
                        body,
                        flags,
                    )?;
                    return match events::resolve_event_kind(rsmi) {
                        RiseSetEvent::Rise | RiseSetEvent::Set => {
                            events::fictitious_rise_set_ut_serving(&serving, &req)
                        }
                        RiseSetEvent::UpperTransit | RiseSetEvent::LowerTransit => {
                            events::fictitious_transit_ut_serving(
                                &serving,
                                jd_ut,
                                geopos[0],
                                events::resolve_event_kind(rsmi) == RiseSetEvent::LowerTransit,
                            )
                            .map(|jd| (0, events::make_tret(jd)))
                        }
                    };
                }
                let serving =
                    events::SealedTransitServing::try_new(stack.serving, &dt, body, flags)?;
                match events::resolve_event_kind(rsmi) {
                    RiseSetEvent::Rise | RiseSetEvent::Set => {
                        events::rise_set_ut_serving(&serving, &req)
                    }
                    RiseSetEvent::UpperTransit | RiseSetEvent::LowerTransit => {
                        events::transit_ut_serving(
                            &serving,
                            jd_ut,
                            geopos[0],
                            events::resolve_event_kind(rsmi) == RiseSetEvent::LowerTransit,
                        )
                        .map(|jd| (0, events::make_tret(jd)))
                    }
                }
            });
        }
        validate_direction(body, jd_ut, &|ut, id| {
            position_with_context(self, ut, id, FLG_SPEED, resolve).map(|_| ())
        })?;
        let kind = events::resolve_event_kind(rsmi);
        // Horizon/atmosphere setup precedes the consumed observer scope.
        let parameters = event_parameters(body, &req);
        let mut owned = self.clone();
        if matches!(kind, RiseSetEvent::Rise | RiseSetEvent::Set)
            && rsmi & BIT_GEOCTR_NO_ECL_LAT == 0
        {
            // Geographic validation consumes latitude before longitude.
            if !(-90.0..=90.0).contains(&geopos[1]) {
                return Err(Error::InvalidObserver {
                    field: "lat_deg",
                    value: geopos[1],
                });
            }
            owned.ctx = owned
                .ctx
                .clone()
                .with_observer(Observer::new(geopos[0], geopos[1], geopos[2])?);
        }
        let calc = |ut, callback_flags| {
            position_with_context(&owned, ut, body, callback_flags, resolve).map(|p| p.values)
        };
        let altitude = |ut, pos: [f64; 3]| {
            self.azalt(
                ut,
                crate::events::ECL2HOR,
                geopos[0],
                geopos[1],
                geopos[2],
                0.0,
                0.0,
                pos[0],
                pos[1],
            )
            .map(|horizontal| horizontal.true_altitude)
        };
        native_event_callbacks(
            body,
            &req,
            &parameters,
            &calc,
            &altitude,
            &projected_ra_dec,
            &default_gast_hours,
        )
    }
}

/// Acquire one public position after an optional caller-owned configuration
/// snapshot. `ut` is the original scalar UT1 JD and `flags` the actual callback
/// word; output is the public six components/retflags in its requested frame.
/// A resolver error precedes any consuming position clock/provider call.
/// With no resolver the original Session/context is borrowed directly. Both
/// paths retain calc_ut's error, split-state and independent Topos contracts;
/// no numeric value is injected and no callback borrow or snapshot escapes.
fn position_with_context(
    session: &Session,
    ut: f64,
    body: i32,
    flags: i32,
    resolve: Option<&PositionContextResolver<'_>>,
) -> Result<Position, Error> {
    match resolve {
        Some(resolve) => {
            Session::from_context(resolve(&session.ctx, body)?).calc_ut(ut, body, flags)
        }
        None => session.calc_ut(ut, body, flags),
    }
}

/// The classical planet map includes Sun through Pluto and Earth; other IDs
/// require the public positionability probe before direction validation.
fn classical_planet(body: i32) -> bool {
    (SUN..=PLUTO).contains(&body) || body == EARTH
}

/// Validate that the body has a defined direction.
/// Bodies outside the classical map are probed through the supplied UT1
/// position calculation. Coverage/corruption errors retain their native
/// identity; unknown objects retain the original cause. Earth and ECL_NUT
/// have no direction for this operation. The synchronous probe is not retained.
fn validate_direction(
    body: i32,
    ut: f64,
    probe: &dyn Fn(f64, i32) -> Result<(), Error>,
) -> Result<(), Error> {
    if !classical_planet(body) {
        probe(ut, body).map_err(|error| {
            // The public preliminary catch preserves fatal coverage/corruption.
            // Other caught failures belong to the source IllegalBody raise-from;
            // retain both owned channels for its binding constructor.
            if matches!(
                error.primary_error(),
                Error::EphemerisRange { .. }
                    | Error::LebOutOfRange { .. }
                    | Error::LebCorrupt { .. }
            ) {
                error
            } else {
                Error::UnknownBody { body_id: body }.with_cause(error)
            }
        })?;
    }
    if matches!(body, EARTH | ECL_NUT) {
        return Err(Error::InvalidBody {
            body_id: body,
            operation: "rise_trans",
        });
    }
    Ok(())
}

/// Checked public value arguments, including atmosphere/horizon inputs ignored
/// by transit. Geographic range checks occur only on a consumed Topos observer.
pub(super) fn validate_event_inputs(req: &RiseSetRequest) -> Result<(), Error> {
    validate_event_inputs_with_latitude(req, true)
}

/// Retain unused transit latitude while checking every other scalar in the
/// original order. Event precedence decides whether latitude is consumed.
fn validate_consumed_event_inputs(req: &RiseSetRequest) -> Result<(), Error> {
    let consumes_latitude = matches!(
        events::resolve_event_kind(req.rsmi),
        RiseSetEvent::Rise | RiseSetEvent::Set
    );
    validate_event_inputs_with_latitude(req, consumes_latitude)
}

/// Shared ordered finite admission; skipping latitude neither changes it nor
/// skips the height, atmosphere or horizon checks that follow it.
fn validate_event_inputs_with_latitude(
    req: &RiseSetRequest,
    check_latitude: bool,
) -> Result<(), Error> {
    for (index, value) in [
        req.jd_start_ut,
        req.lon_deg_east,
        req.lat_deg,
        req.height_m,
        req.atpress_mbar,
        req.attemp_c,
        req.horhgt_deg,
    ]
    .into_iter()
    .enumerate()
    {
        if (index != 2 || check_latitude) && !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    Ok(())
}

/// GAST in [0,24) hours from a fresh bundled UT1 Time. The original UT1
/// scalar remains `whole`; two bundled Delta-T evaluations supply the TT
/// fraction, while the UT1 fraction is explicitly `0*ut` (positive zero at
/// positive UT1). There is
/// no third Delta-T read or inverse reconstruction of the UT1 fraction.
/// Nutation/complementary terms use collapsed TT; mean obliquity and GMST
/// precession use TDB derived from the retained TT split. The configured
/// model, observer and caller position clock are independent. Nonfinite
/// clock/frame values and acquisition errors propagate.
pub(super) fn default_gast_hours(ut: f64) -> Result<f64, Error> {
    default_gast_with(
        ut,
        &crate::smh2016::timescale_parts_from_ut1,
        &gast_at_parts,
    )
}

/// Pass a freshly acquired UT1 Time split to the sidereal tail without an
/// inverse clock operation. `parts` returns `(original UT1 whole, TT fraction)`
/// in days after two table reads; `gast` receives the same words and an
/// explicit `0*ut` UT1 fraction, returning hours. Both callbacks are
/// synchronous borrowed acquisitions, retain no foreign state and propagate
/// their original errors; a failed clock acquisition cannot consume GAST.
fn default_gast_with(
    ut: f64,
    parts: &dyn Fn(f64) -> Result<(f64, f64), Error>,
    gast: &dyn Fn(f64, f64, f64) -> Result<f64, Error>,
) -> Result<f64, Error> {
    let (whole, tt_fraction) = parts(ut)?;
    gast(whole, tt_fraction, 0.0 * ut)
}

/// Apparent sidereal tail for explicit whole, TT-fraction and UT1-fraction
/// Julian-day words. GMST retains both fractions; seven-term TDB supplies
/// P03 mean obliquity. IAU2000A dpsi and complementary terms use scalar TT.
/// With angles in radians, GAST=mod24(GMST+(dpsi*cos(eps)+C)*24/tau).
/// Returns hours in [0,24); finite inputs/collapsed TT/TDB are required.
/// There is no provider, polar motion, local-site or Delta-T acquisition.
fn gast_at_parts(whole: f64, tt_fraction: f64, ut1_fraction: f64) -> Result<f64, Error> {
    use crate::positions::reference_frames as sf;
    let jd_tt = whole + tt_fraction;
    let jd_tdb =
        whole + (tt_fraction + crate::time::tdb_minus_tt_seconds(whole, tt_fraction)? / 86400.0);
    for value in [whole, tt_fraction, ut1_fraction, jd_tt, jd_tdb] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    let gmst = sf::gmst_hours_two_part(whole, tt_fraction, ut1_fraction)?;
    let (dpsi_u, _) = sf::iau2000a_tenth_microarcsec(jd_tt)?;
    let dpsi = dpsi_u * sf::FRAME_TENTH_USEC_2_RAD;
    let mean_obl = sf::mean_obliquity_arcsec(jd_tdb)? * sf::FRAME_ASEC2RAD;
    let c_terms = sf::equinoxes_complimentary_terms_radians(jd_tt)?;
    let eq_eq = dpsi * mean_obl.cos() + c_terms;
    Ok(sf::floored_mod(gmst + eq_eq / sf::FRAME_TAU * 24.0, 24.0))
}

/// Project an apparent geocentric ecliptic place onto latitude zero, then to
/// true equator-of-date at bundled scalar TT acquired from original UT1.
/// Returns RA hours/Dec degrees; input longitude is degrees and distance AU.
/// Vondrak mean obliquity plus IAU2006 nutation in obliquity supplies `-eps`
/// to cotrans. This is independent of the body's configured model clock.
pub(super) fn projected_ra_dec(ut: f64, pos: [f64; 3]) -> Result<(f64, f64), Error> {
    let tt = crate::time::default_timescale_tt_from_ut1(ut)?;
    let (_, deps) = crate::erfa::nut06a(tt)?;
    let eps = crate::positions::vondrak_mean_obliquity_deg(tt)? + deps.to_degrees();
    let (ra, dec, _) = crate::events::cotrans(pos[0], 0.0, pos[2], -eps);
    Ok((ra / 15.0, dec))
}

/// Atmosphere/disc/horizon choices resolved before an ordinary observer is
/// constructed. Pressure is mbar, horizon degrees and disc radius AU; effective
/// modifiers preserve centre/twilight precedence. This owns no provider/time
/// acquisition. Only horizon -100 consumes dip here; per-height refraction
/// acquires its dip after the actual body altitude. Floating intermediate words
/// are subsequently checked by the event kernels; no physical bound is proved.
#[derive(Clone, Copy)]
pub(super) struct EventParameters {
    pub(super) pressure: f64,
    pub(super) horizon: f64,
    pub(super) effective: i32,
    pub(super) disc: events::DiscPoint,
    pub(super) radius: f64,
    pub(super) lift: f64,
}

/// Resolve pressure0=1013.25*max(1-0.0065*h/288,0)^5.255 and the literal
/// horizon (-100 becomes dip(max(h,0)) +1e-4 degree). Sun twilight replaces
/// that horizon and forces an unrefracted centre. Fixed-disc choice changes
/// distance only when evaluating an actual height. Transit receives the same
/// preliminary setup but consumes no disc/refraction callback.
fn event_parameters(body: i32, req: &RiseSetRequest) -> EventParameters {
    let pressure = if req.atpress_mbar == 0.0 {
        events::pressure_from_altitude_m(req.height_m)
    } else {
        req.atpress_mbar
    };
    let mut horizon = req.horhgt_deg;
    if horizon == -100.0 {
        let dip = crate::refraction::calc_dip(
            req.height_m.max(0.0),
            events::RISE_SET_DIP_LAPSE_RATE_K_M,
            pressure,
            req.attemp_c,
        );
        horizon = events::resolve_horizon_altitude_deg(horizon, dip);
    }
    let mut effective = req.rsmi;
    if let Some(twilight) = events::twilight_override(req.rsmi, body == SUN) {
        effective |= BIT_NO_REFRACTION | BIT_DISC_CENTER;
        horizon = twilight.horizon_alt_deg;
    }
    let disc = events::select_disc_point(effective, false);
    let radius = events::body_disc_radius_au(body);
    let mut lift = 0.0;
    if effective & BIT_NO_REFRACTION == 0 {
        lift += 0.75;
    }
    if disc != events::DiscPoint::Center {
        lift += 0.30;
    }
    EventParameters {
        pressure,
        horizon,
        effective,
        disc,
        radius,
        lift,
    }
}

/// Shared Native public callback orchestration. `calc` returns the six actual
/// public spherical components (degrees/AU and their daily rates) at each
/// original UT1 sample and the fixed callback flag word. `altitude` consumes
/// lon/lat/dist in degrees/degrees/AU and returns true altitude degrees via
/// spherical azalt; `projected` returns Hindu RA hours/Dec degrees. `gast`
/// returns independent default GAST hours. All borrows are synchronous and
/// retained only through this search; acquisition errors propagate unchanged.
///
/// Disc lift is signed asin(min(1,R/dist)) in degrees; fixed disc uses Sun1 AU
/// or Moon0.00257 AU. Sinclair refraction solves app=true+R(app) through the
/// existing ten-step Newton inversion and dip-floor branch. Pressure0 uses
/// 1013.25*max(1-0.0065*h/288,0)^5.255. The refraction/transit primitives
/// remain shared; the private Native scan owns its literal sample count.
pub(super) fn native_event_callbacks(
    body: i32,
    req: &RiseSetRequest,
    parameters: &EventParameters,
    calc: &dyn Fn(f64, i32) -> Result<[f64; 6], Error>,
    altitude: &dyn Fn(f64, [f64; 3]) -> Result<f64, Error>,
    projected: &dyn Fn(f64, [f64; 3]) -> Result<(f64, f64), Error>,
    gast: &dyn Fn(f64) -> Result<f64, Error>,
) -> Result<(i32, [f64; 10]), Error> {
    native_event_callbacks_with_policy(
        body, req, parameters, calc, altitude, projected, gast, false,
    )
}

/// Stellar source-stage policy: setup is already resolved, then catalog/place
/// callbacks consume their epochs before any solver-derived checks. Transit
/// never consumes latitude/height/atmosphere; rise/set begins with declination.
/// All callbacks are synchronous; checked planetary callers retain their guards.
#[allow(clippy::too_many_arguments)]
pub(super) fn native_star_event_callbacks(
    body: i32,
    req: &RiseSetRequest,
    parameters: &EventParameters,
    calc: &dyn Fn(f64, i32) -> Result<[f64; 6], Error>,
    altitude: &dyn Fn(f64, [f64; 3]) -> Result<f64, Error>,
    projected: &dyn Fn(f64, [f64; 3]) -> Result<(f64, f64), Error>,
    gast: &dyn Fn(f64) -> Result<f64, Error>,
) -> Result<(i32, [f64; 10]), Error> {
    native_event_callbacks_with_policy(body, req, parameters, calc, altitude, projected, gast, true)
}

/// Shared callback arithmetic and output layout. `source_star` selects only
/// consumed-stage admission; false retains the checked planetary API exactly.
/// Public places/altitude/frames own their errors; no value or epoch is replaced.
#[allow(clippy::too_many_arguments)]
fn native_event_callbacks_with_policy(
    body: i32,
    req: &RiseSetRequest,
    parameters: &EventParameters,
    calc: &dyn Fn(f64, i32) -> Result<[f64; 6], Error>,
    altitude: &dyn Fn(f64, [f64; 3]) -> Result<f64, Error>,
    projected: &dyn Fn(f64, [f64; 3]) -> Result<(f64, f64), Error>,
    gast: &dyn Fn(f64) -> Result<f64, Error>,
    source_star: bool,
) -> Result<(i32, [f64; 10]), Error> {
    let kind = events::resolve_event_kind(req.rsmi);
    let ra_dec = |ut| calc(ut, FLG_EQUATORIAL | FLG_SPEED).map(|p| (p[0] / 15.0, p[1]));
    if matches!(
        kind,
        RiseSetEvent::UpperTransit | RiseSetEvent::LowerTransit
    ) {
        if source_star {
            return super::star_rise::star_calculate_transit(
                req.jd_start_ut,
                req.lon_deg_east,
                kind == RiseSetEvent::LowerTransit,
                &ra_dec,
                gast,
            );
        }
        return events::calculate_transit(
            req.jd_start_ut,
            req.lon_deg_east,
            kind == RiseSetEvent::LowerTransit,
            &mut |ut| ra_dec(ut),
            &mut |ut| gast(ut),
        )
        .map(|solution| (0, events::make_tret(solution.jd_ut)));
    }
    let EventParameters {
        pressure,
        horizon,
        effective,
        disc,
        radius,
        lift,
    } = *parameters;
    let refract = effective & BIT_NO_REFRACTION == 0;
    let hindu = effective & BIT_GEOCTR_NO_ECL_LAT != 0;
    let mut height = |ut| {
        let flags = if hindu {
            FLG_SPEED
        } else {
            FLG_TOPOCTR | FLG_SPEED
        };
        let p = calc(ut, flags)?;
        let place = [p[0], if hindu { 0.0 } else { p[1] }, p[2]];
        let mut distance = p[2];
        let mut alt = altitude(ut, place)?;
        if effective & BIT_FIXED_DISC_SIZE != 0 {
            if body == SUN {
                distance = events::FIXED_DISC_SUN_DIST_AU;
            } else if body == MOON {
                distance = events::FIXED_DISC_MOON_DIST_AU;
            }
        }
        alt += events::disc_altitude_term_deg(disc, radius, distance)?;
        if refract {
            let dip = crate::refraction::calc_dip(
                if source_star {
                    super::star_rise::star_nonnegative(req.height_m)
                } else {
                    req.height_m.max(0.0)
                },
                events::RISE_SET_DIP_LAPSE_RATE_K_M,
                pressure,
                req.attemp_c,
            );
            alt = events::true_to_apparent_altitude(
                alt,
                dip,
                &mut |app| {
                    Ok(crate::refraction::sinclair_refraction_deg(
                        app,
                        pressure,
                        req.attemp_c,
                    ))
                },
                ut,
            )?;
        }
        Ok(alt - horizon)
    };
    let mut precheck = |ut| {
        if hindu {
            let p = calc(ut, FLG_SPEED)?;
            projected(ut, [p[0], 0.0, p[2]])
        } else {
            ra_dec(ut)
        }
    };
    match native_rise_set_search_with_policy(
        req.jd_start_ut,
        req.lat_deg,
        kind == RiseSetEvent::Rise,
        &mut height,
        &mut precheck,
        horizon,
        lift,
        !source_star,
    )? {
        RiseSetOutcome::Found { jd_ut } => Ok((0, events::make_tret(jd_ut))),
        RiseSetOutcome::Circumpolar => Ok((-2, events::make_tret(0.0))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::RefCell;

    fn stored_model_session(mode: ProviderMode) -> Session {
        let path = crate::test_data::core_leb_path();
        let image = crate::leb::LebFile::open(&path).expect("pinned core fixture");
        let mut context = EphemerisContext::new()
            .with_mode(mode)
            .with_network_policy(Some(NetworkPolicy::Sealed))
            .with_user_delta_t(Some(64.0 / 86400.0))
            .with_iers_delta_t(Some(false));
        context.attach_leb_file(image);
        Session::from_context(context)
    }

    #[test]
    fn stored_harrington_transit_uses_analytical_model() {
        let session = stored_model_session(ProviderMode::SealedLeb);
        let (status, times) = session
            .rise_trans(
                2451545.0,
                HARRINGTON,
                CALC_MTRANSIT,
                [12.5, 41.9, 150.0],
                0.0,
                0.0,
                FLG_SWIEPH,
            )
            .unwrap();
        assert_eq!(status, 0);
        assert!(times[0].is_finite() && (2451545.0..2451547.0).contains(&times[0]));
        assert!(times[1..].iter().all(|value| value.to_bits() == 0));
    }

    fn assert_event_layout(start: f64, result: (i32, [f64; 10])) {
        let (status, times) = result;
        assert!(times[1..].iter().all(|value| value.to_bits() == 0));
        match status {
            0 => assert!(times[0].is_finite() && (start..start + 2.0).contains(&times[0])),
            -2 => assert_eq!(times[0].to_bits(), 0),
            _ => panic!("unexpected event status {status}"),
        }
    }

    #[test]
    fn stored_supported_models_reach_each_event_kind_in_leb_and_auto() {
        let bodies = (CUPIDO..=WALDEMATH)
            .filter(|body| events::is_fictitious_event_body(*body))
            .collect::<Vec<_>>();
        assert_eq!(bodies.len(), 18);
        let kinds = [CALC_RISE, CALC_SET, CALC_MTRANSIT, CALC_ITRANSIT];
        for mode in [ProviderMode::SealedLeb, ProviderMode::Auto] {
            let session = stored_model_session(mode);
            for (index, body) in bodies.iter().copied().enumerate() {
                let result = session
                    .rise_trans(
                        2451545.0,
                        body,
                        kinds[index % kinds.len()],
                        [12.5, 41.9, 150.0],
                        0.0,
                        10.0,
                        FLG_SWIEPH,
                    )
                    .unwrap_or_else(|error| panic!("{mode:?}, body {body}: {error:?}"));
                assert_event_layout(2451545.0, result);
            }
        }
    }

    #[test]
    fn stored_model_horizon_and_atmosphere_retain_crossing_direction() {
        let session = stored_model_session(ProviderMode::SealedLeb);
        // Begin at culmination so different setting horizons select the same
        // daily crossing instead of straddling the start of a search window.
        let start = session
            .rise_trans(
                2451545.0,
                HARRINGTON,
                CALC_MTRANSIT,
                [12.5, 41.9, 1500.0],
                0.0,
                10.0,
                FLG_SWIEPH,
            )
            .unwrap()
            .1[0];
        let event = |kind, pressure, temperature, horizon| {
            let result = session
                .rise_trans_true_hor(
                    start,
                    HARRINGTON,
                    kind | BIT_DISC_CENTER,
                    [12.5, 41.9, 1500.0],
                    pressure,
                    temperature,
                    horizon,
                    FLG_SWIEPH,
                )
                .unwrap();
            assert_eq!(result.0, 0);
            assert_event_layout(start, result);
            result.1[0]
        };
        // A higher geometric horizon delays rise and advances set. These
        // inequalities test request forwarding without a stored event answer.
        let rise = CALC_RISE | BIT_NO_REFRACTION;
        let set = CALC_SET | BIT_NO_REFRACTION;
        assert!(event(rise, 900.0, 10.0, -2.0) < event(rise, 900.0, 10.0, 2.0));
        assert!(event(set, 900.0, 10.0, 2.0) < event(set, 900.0, 10.0, -2.0));
        assert!(event(rise, 900.0, 10.0, -100.0) < event(rise, 900.0, 10.0, 0.0));
        // More pressure or colder air increases refraction and advances rise.
        assert!(event(CALC_RISE, 1000.0, 10.0, 0.0) < event(CALC_RISE, 400.0, 10.0, 0.0));
        assert!(event(CALC_RISE, 900.0, -20.0, 0.0) < event(CALC_RISE, 900.0, 30.0, 0.0));
        assert!(event(CALC_RISE, 0.0, 10.0, 0.0) < event(rise, 0.0, 10.0, 0.0));
        assert_event_layout(
            2451545.0,
            session
                .rise_trans(
                    2451545.0,
                    HARRINGTON,
                    CALC_RISE | BIT_HINDU_RISING,
                    [12.5, 41.9, 150.0],
                    0.0,
                    10.0,
                    FLG_SWIEPH,
                )
                .unwrap(),
        );
    }

    #[test]
    fn stored_model_transit_keeps_kind_precedence_and_unused_latitude() {
        let session = stored_model_session(ProviderMode::SealedLeb);
        let event = |kind, latitude, flags| {
            session.rise_trans_true_hor_consumed(
                2451545.0,
                HARRINGTON,
                kind,
                [12.5, latitude, 150.0],
                900.0,
                10.0,
                -100.0,
                flags,
            )
        };
        let upper = event(CALC_MTRANSIT, 41.9, FLG_SWIEPH).unwrap();
        let lower = event(CALC_ITRANSIT, 41.9, FLG_SWIEPH).unwrap();
        assert_ne!(upper.1[0], lower.1[0]);
        assert_eq!(
            event(CALC_MTRANSIT, f64::NAN, FLG_SWIEPH | ECL_ONE_TRY).unwrap(),
            upper
        );
        assert_eq!(
            event(CALC_MTRANSIT | CALC_ITRANSIT, 41.9, FLG_SWIEPH).unwrap(),
            lower
        );
        let set = event(CALC_SET, 41.9, FLG_SWIEPH).unwrap();
        assert_eq!(
            event(CALC_SET | CALC_RISE | CALC_MTRANSIT, 41.9, FLG_SWIEPH).unwrap(),
            set
        );
        assert!(matches!(
            event(CALC_RISE, f64::NAN, FLG_SWIEPH),
            Err(Error::NonFiniteArgument { value }) if value.is_nan()
        ));
    }

    #[test]
    fn stored_model_event_flags_and_provider_failures_remain_errors() {
        let session = stored_model_session(ProviderMode::SealedLeb);
        let event = |session: &Session, jd, body, flags| {
            session.rise_trans(
                jd,
                body,
                CALC_MTRANSIT,
                [12.5, 41.9, 150.0],
                0.0,
                10.0,
                flags,
            )
        };
        for flag in [
            FLG_SIDEREAL,
            FLG_J2000,
            FLG_NONUT,
            FLG_XYZ,
            FLG_RADIANS,
            FLG_HELCTR,
            FLG_BARYCTR,
        ] {
            assert!(matches!(
                event(&session, 2451545.0, HARRINGTON, FLG_SWIEPH | flag),
                Err(Error::UnsupportedFlagCombination { .. })
            ));
        }
        for body in [NIBIRU, 59, 69] {
            let error = event(&session, 2451545.0, body, FLG_SWIEPH).unwrap_err();
            let expected = if body == NIBIRU {
                Error::UnknownBody { body_id: body }
            } else {
                Error::LebBodyNotFound { body_id: body }
            };
            assert_eq!(error.primary_error(), &expected);
        }
        let range = event(&session, 1_000_000_000.0, HARRINGTON, FLG_SWIEPH).unwrap_err();
        assert!(matches!(
            range.primary_error(),
            Error::EphemerisRange { .. } | Error::LebOutOfRange { .. }
        ));
        let mut context = session.context().clone();
        context.release_leb_images();
        let empty = Session::from_context(context);
        assert!(matches!(
            event(&empty, 2451545.0, HARRINGTON, FLG_SWIEPH),
            Err(Error::ConfigurationMissing {
                missing_config: "leb_reader",
                ..
            })
        ));
        assert!(matches!(
            event(&empty, f64::NAN, HARRINGTON, FLG_SWIEPH),
            Err(Error::NonFiniteArgument { value }) if value.is_nan()
        ));
        // Closed input admission still precedes the selected model's flags.
        let path = crate::test_data::core_leb_path();
        let mut image = crate::leb::LebFile::open(&path).unwrap();
        image.close();
        let mut context = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        context.attach_leb_file(image);
        let closed = Session::from_context(context);
        assert!(matches!(
            event(&closed, 2451545.0, HARRINGTON, FLG_SWIEPH | FLG_XYZ),
            Err(Error::LebReaderClosed)
        ));
    }

    fn event(
        body: i32,
        req: &RiseSetRequest,
        calc: &dyn Fn(f64, i32) -> Result<[f64; 6], Error>,
        altitude: &dyn Fn(f64, [f64; 3]) -> Result<f64, Error>,
        projected: &dyn Fn(f64, [f64; 3]) -> Result<(f64, f64), Error>,
        gast: &dyn Fn(f64) -> Result<f64, Error>,
    ) -> Result<(i32, [f64; 10]), Error> {
        native_event_callbacks(
            body,
            req,
            &event_parameters(body, req),
            calc,
            altitude,
            projected,
            gast,
        )
    }
    fn req(rsmi: i32) -> RiseSetRequest {
        RiseSetRequest {
            jd_start_ut: 2451545.125,
            lon_deg_east: 0.0,
            lat_deg: 0.0,
            height_m: 0.0,
            atpress_mbar: 1013.25,
            attemp_c: 0.0,
            horhgt_deg: 0.0,
            rsmi,
        }
    }
    #[test]
    fn earth_refuses_without_probe_but_ecl_nut_keeps_prior_range() {
        assert!(matches!(
            validate_direction(EARTH, 2451545.125, &|_, _| panic!("Earth acquired")),
            Err(Error::InvalidBody { body_id: EARTH, .. })
        ));
        let error = Error::EphemerisRange {
            requested_jd: 2451545.125,
            start_jd: 2451546.0,
            end_jd: 2451600.0,
            body: Some(ECL_NUT),
        };
        assert_eq!(
            validate_direction(ECL_NUT, 2451545.125, &|ut, id| {
                assert_eq!(ut, 2451545.125);
                assert_eq!(id, ECL_NUT);
                Err(error.clone())
            })
            .unwrap_err(),
            error
        );
    }
    #[test]
    fn out_of_map_nonrange_probe_failure_keeps_requested_key() {
        assert_eq!(
            validate_direction(999999, 2451545.125, &|_, _| Err(
                Error::SpkUnregisteredBody { body: SUN }
            ))
            .unwrap_err(),
            Error::UnknownBody { body_id: 999999 }
                .with_cause(Error::SpkUnregisteredBody { body: SUN })
        );
    }

    #[test]
    fn preliminary_wrapped_range_keeps_typed_words_and_context() {
        let error = Error::EphemerisRange {
            requested_jd: -0.0,
            start_jd: 0.0,
            end_jd: 1.0,
            body: Some(999999),
        }
        .with_cause(Error::UnknownBody { body_id: 999998 });
        let actual = validate_direction(999999, 2451545.0, &|_, _| Err(error.clone())).unwrap_err();
        assert_eq!(actual, error);
        assert!(
            matches!(actual.primary_error(), Error::EphemerisRange { requested_jd, .. }
            if requested_jd.to_bits() == (-0.0_f64).to_bits())
        );
    }

    #[test]
    fn preliminary_leb_corruption_is_not_illegal_body() {
        let error = Error::LebCorrupt {
            what: "controlled broken directory".into(),
            offset: 24,
        };
        let actual = validate_direction(999999, 2451545.0, &|_, _| Err(error.clone())).unwrap_err();
        assert_eq!(actual, error);
        assert!(actual.cause().is_none());
    }
    #[test]
    fn ordinary_callback_flags_and_no_event_zero_words() {
        let calls = RefCell::new(Vec::new());
        let calc = |_: f64, flags| {
            calls.borrow_mut().push(flags);
            Ok([37.0, 12.0, 1.0, 0.0, 0.0, 0.0])
        };
        let result = event(
            MARS,
            &req(CALC_RISE | BIT_DISC_CENTER | BIT_NO_REFRACTION),
            &calc,
            &|_, _| Ok(-90.0),
            &|_, _| panic!("ordinary path projected"),
            &|_| panic!("rise path acquired GAST"),
        )
        .unwrap();
        assert_eq!(result.0, -2);
        assert!(result.1.iter().all(|v| v.to_bits() == 0));
        assert_eq!(calls.borrow()[0], FLG_EQUATORIAL | FLG_SPEED);
        assert!(calls.borrow().contains(&(FLG_TOPOCTR | FLG_SPEED)));
        assert!(
            calls
                .borrow()
                .iter()
                .all(|f| *f == FLG_EQUATORIAL | FLG_SPEED || *f == FLG_TOPOCTR | FLG_SPEED)
        );
    }
    #[test]
    fn hindu_projection_drops_latitude_and_topocentric_flag() {
        let calls = RefCell::new(Vec::new());
        let calc = |_: f64, flags| {
            calls.borrow_mut().push(flags);
            Ok([37.0, 12.0, 1.0, 0.0, 0.0, 0.0])
        };
        let result = event(
            MOON,
            &req(CALC_RISE | BIT_HINDU_RISING),
            &calc,
            &|_, p| {
                assert_eq!(p[1].to_bits(), 0);
                Ok(-90.0)
            },
            &|_, p| {
                assert_eq!(p[1].to_bits(), 0);
                Ok((1.0, 0.0))
            },
            &|_| panic!("Hindu acquired GAST"),
        )
        .unwrap();
        assert_eq!(result.0, -2);
        assert!(calls.borrow().iter().all(|flags| *flags == FLG_SPEED));
    }
    #[test]
    fn transit_uses_geocentric_equatorial_before_gast_and_ignores_hindu() {
        let order = RefCell::new(Vec::new());
        let calc = |_: f64, flags| {
            assert_eq!(flags, FLG_EQUATORIAL | FLG_SPEED);
            order.borrow_mut().push("calc");
            Ok([15.0, 0.0, 1.0, 0.0, 0.0, 0.0])
        };
        let result = event(
            MARS,
            &req(CALC_MTRANSIT | BIT_HINDU_RISING),
            &calc,
            &|_, _| panic!("transit altitude"),
            &|_, _| panic!("transit projection"),
            &|_| {
                order.borrow_mut().push("gast");
                Ok(1.0)
            },
        )
        .unwrap();
        assert_eq!(result.0, 0);
        assert_eq!(result.1[0], 2451545.125);
        assert_eq!(&order.borrow()[..2], &["calc", "gast"]);
    }
    #[test]
    fn acquired_support_error_is_not_no_event() {
        let error = Error::SpkUnregisteredBody { body: JUPITER };
        let result = event(
            MARS,
            &req(CALC_RISE),
            &|_, _| Err(error.clone()),
            &|_, _| panic!("after error"),
            &|_, _| panic!("after error"),
            &|_| panic!("after error"),
        );
        assert_eq!(result.unwrap_err(), error);
    }
    #[test]
    fn nonfinite_horizon_is_checked_even_for_transit() {
        let mut r = req(CALC_MTRANSIT);
        r.horhgt_deg = f64::NAN;
        assert!(matches!(
            validate_event_inputs(&r),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[test]
    fn consumed_admission_retains_only_unused_transit_latitude() {
        for kind in [CALC_MTRANSIT, CALC_ITRANSIT, CALC_MTRANSIT | CALC_ITRANSIT] {
            for latitude in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                let mut r = req(kind);
                r.lat_deg = latitude;
                validate_consumed_event_inputs(&r).unwrap();
                assert!(matches!(
                    validate_event_inputs(&r),
                    Err(Error::NonFiniteArgument { .. })
                ));
                let outcome = event(
                    SUN,
                    &r,
                    &|_, flags| {
                        assert_eq!(flags, FLG_EQUATORIAL | FLG_SPEED);
                        Ok([15.0, 0.0, 1.0, 0.0, 0.0, 0.0])
                    },
                    &|_, _| panic!("transit consumed altitude"),
                    &|_, _| panic!("transit consumed projection"),
                    &|_| Ok(1.0),
                )
                .unwrap();
                assert_eq!(outcome.0, 0);
                assert_eq!(r.lat_deg.to_bits(), latitude.to_bits());
            }
        }
        for kind in [0, CALC_RISE, CALC_SET, CALC_RISE | CALC_MTRANSIT] {
            let mut r = req(kind);
            r.lat_deg = f64::NAN;
            assert!(matches!(
                validate_consumed_event_inputs(&r),
                Err(Error::NonFiniteArgument { .. })
            ));
        }
    }

    #[test]
    fn consumed_transit_admission_keeps_other_finite_guards_in_order() {
        for index in [0, 1, 3, 4, 5, 6] {
            let mut r = req(CALC_MTRANSIT);
            r.lat_deg = f64::NAN;
            let refused = if index < 2 {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            };
            match index {
                0 => r.jd_start_ut = refused,
                1 => r.lon_deg_east = refused,
                3 => r.height_m = refused,
                4 => r.atpress_mbar = refused,
                5 => r.attemp_c = refused,
                6 => r.horhgt_deg = refused,
                _ => unreachable!(),
            }
            let error = validate_consumed_event_inputs(&r).unwrap_err();
            assert!(
                matches!(error, Error::NonFiniteArgument { value } if value.to_bits() == refused.to_bits())
            );
        }
    }

    #[test]
    fn checked_and_consumed_session_admission_remain_distinct() {
        let session =
            Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
        for kind in [CALC_MTRANSIT, CALC_ITRANSIT] {
            let checked = session.rise_trans_true_hor(
                2451545.0,
                EARTH,
                kind,
                [0.0, f64::NAN, 0.0],
                0.0,
                0.0,
                0.0,
                FLG_SWIEPH,
            );
            assert!(matches!(checked, Err(Error::NonFiniteArgument { .. })));
            let consumed = session.rise_trans_true_hor_consumed(
                2451545.0,
                EARTH,
                kind,
                [0.0, f64::NAN, 0.0],
                0.0,
                0.0,
                0.0,
                FLG_SWIEPH,
            );
            assert!(matches!(
                consumed,
                Err(Error::InvalidBody { body_id: EARTH, .. })
            ));
        }
    }

    #[test]
    fn fresh_gast_passes_original_ut_whole_and_positive_zero_fraction() {
        let order = RefCell::new(Vec::new());
        let got = default_gast_with(
            1000.25,
            &|u| {
                order.borrow_mut().push("parts");
                assert_eq!(u, 1000.25);
                Ok((u, 12.0 / 86400.0))
            },
            &|whole, tt_fraction, ut1_fraction| {
                order.borrow_mut().push("gast");
                assert_eq!(whole, 1000.25);
                assert_eq!(tt_fraction, 12.0 / 86400.0);
                assert_eq!(ut1_fraction.to_bits(), 0);
                Ok(3.0)
            },
        )
        .unwrap();
        assert_eq!(got, 3.0);
        assert_eq!(&order.borrow()[..], &["parts", "gast"]);
    }
    #[test]
    fn failed_fresh_time_cannot_acquire_gast() {
        let error = Error::SpkUnregisteredBody { body: SUN };
        assert_eq!(
            default_gast_with(1000.25, &|_| Err(error.clone()), &|_, _, _| panic!(
                "GAST acquired after clock error"
            ))
            .unwrap_err(),
            error
        );
    }

    #[test]
    fn native_context_resolver_is_lazy_for_early_refusals() {
        let session =
            Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
        let calls = std::cell::Cell::new(0);
        let resolve = |_: &EphemerisContext, _: i32| {
            calls.set(calls.get() + 1);
            Err(Error::SpkUnregisteredBody { body: 987654 })
        };
        for (body, lat, ut) in [
            (EARTH, 41.9, 2451545.125),
            (MARS, 91.0, 2451545.125),
            (MARS, 41.9, f64::NAN),
        ] {
            assert!(
                session
                    .rise_trans_true_hor_with_position_context(
                        ut,
                        body,
                        1,
                        [12.5, lat, 0.0],
                        0.0,
                        0.0,
                        0.0,
                        FLG_SWIEPH,
                        &resolve,
                    )
                    .is_err()
            );
        }
        assert_eq!(calls.get(), 0);
    }

    #[test]
    fn native_context_resolver_receives_consumed_temporary_observer() {
        let session =
            Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
        let calls = std::cell::Cell::new(0);
        let resolve = |ctx: &EphemerisContext, body| {
            calls.set(calls.get() + 1);
            assert_eq!(body, MARS);
            assert!(ctx.observer().is_some());
            Err(Error::SpkUnregisteredBody { body: 987654 })
        };
        let error = session
            .rise_trans_true_hor_with_position_context(
                2451545.125,
                MARS,
                1,
                [12.5, 41.9, 0.0],
                0.0,
                0.0,
                0.0,
                FLG_SWIEPH,
                &resolve,
            )
            .unwrap_err();
        assert!(matches!(error, Error::SpkUnregisteredBody { body: 987654 }));
        assert_eq!(calls.get(), 1);
        assert!(session.context().observer().is_none());
    }

    #[test]
    fn native_context_resolver_transit_does_not_install_unused_observer() {
        let session =
            Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
        let resolve = |ctx: &EphemerisContext, body| {
            assert_eq!(body, MARS);
            assert!(ctx.observer().is_none());
            Err(Error::SpkUnregisteredBody { body: 987654 })
        };
        let error = session
            .rise_trans_true_hor_with_position_context(
                2451545.125,
                MARS,
                4,
                [12.5, 91.0, 0.0],
                0.0,
                0.0,
                0.0,
                FLG_SWIEPH,
                &resolve,
            )
            .unwrap_err();
        assert!(matches!(error, Error::SpkUnregisteredBody { body: 987654 }));
    }

    #[test]
    fn native_context_resolver_probe_preserves_source_range_catch() {
        let session =
            Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
        let resolve = |_: &EphemerisContext, body| {
            assert_eq!(body, 999999);
            Err(Error::EphemerisRange {
                requested_jd: 2451545.125,
                start_jd: 2451546.0,
                end_jd: 2451547.0,
                body: Some(body),
            })
        };
        let error = session
            .rise_trans_true_hor_with_position_context(
                2451545.125,
                999999,
                4,
                [12.5, 41.9, 0.0],
                0.0,
                0.0,
                0.0,
                FLG_SWIEPH,
                &resolve,
            )
            .unwrap_err();
        assert!(matches!(
            error,
            Error::EphemerisRange {
                body: Some(999999),
                ..
            }
        ));
    }
}
