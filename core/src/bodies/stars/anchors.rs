//! Dynamic sidereal anchor serving over borrowed providers and live row leases.
//!
//! This boundary owns provider priority, per-node cache publication, projected
//! propagation clocks and outer-request range errors as one operation. Static
//! catalog stars and live defining rows remain distinct acquisition channels.

use super::{
    AcquiredAnchorRayAt, AnchorRangeErrorAt, AnchorRayAt, AnchorRayRole, AnchorRowAcquisition,
    AnchorSample, AnchorSampleOptions, AnchorStateAt, Error, SealedAnchorProvider, StarData,
    StarFrame, StarObserver, star_icrs_state_catalog_angles, star_initial_state_from_native,
    topo_star_epoch_days,
};

/// Mean-ecliptic-of-date frame for the live ayanamsha anchors: the biased
/// Vondrák 2011 ICRS-to-mean-equator matrix, including frame bias,
/// with the Vondrák mean obliquity. Nutation never enters: star-based mean
/// ayanamshas use this longitude and the true arms add nutation in
/// longitude exactly once downstream.
fn anchor_mean_frame(jd_tt: f64) -> Result<StarFrame, Error> {
    let pn = crate::positions::vondrak_precession_matrix(jd_tt, true)?;
    let eps = crate::positions::vondrak_mean_obliquity_rad(jd_tt)?;
    StarFrame::new(pn, eps)
}

/// Sample a periodic anchor longitude in degrees at a finite TT Julian Day.
///
/// With active deflection, evaluate exactly at the requested epoch: solar
/// conjunction can produce a sharp feature that must not be interpolated.
/// Otherwise sample the three 8-day nodes enclosing the 16-day grid cell,
/// unwrap the outer values around the central longitude, and evaluate the
/// Lagrange quadratic in its written operation order. The result is reduced
/// to `[0, 360)`. Every node is queried even at a grid endpoint; provider
/// errors propagate without an exact-date fallback or extrapolation.
fn anchor_longitude_sampled_via(
    jd_tt: f64,
    nogdefl: bool,
    at: impl Fn(f64) -> Result<f64, Error>,
) -> Result<f64, Error> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    if !nogdefl {
        return at(jd_tt);
    }
    let h = 8.0;
    let base = (jd_tt / (2.0 * h)).floor() * (2.0 * h);
    let mut y0 = at(base)?;
    let y1 = at(base + h)?;
    let mut y2 = at(base + 2.0 * h)?;
    for value in [y0, y1, y2] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    if y0 - y1 > 180.0 {
        y0 -= 360.0;
    } else if y0 - y1 < -180.0 {
        y0 += 360.0;
    }
    if y2 - y1 > 180.0 {
        y2 -= 360.0;
    } else if y2 - y1 < -180.0 {
        y2 += 360.0;
    }
    let x = (jd_tt - base) / h;
    let lon = y0 * (x - 1.0) * (x - 2.0) / 2.0 - y1 * x * (x - 2.0) + y2 * x * (x - 1.0) / 2.0;
    Ok(lon.rem_euclid(360.0))
}

/// Reduce an astrometric ICRS vector in AU with independent anchor flags.
///
/// Geometric requests return the input; deflection-free requests aberrate
/// only. Otherwise deflect then aberrate. NOABERR alone subtracts the
/// aberration displacement evaluated on the original astrometric vector
/// from the full apparent vector, not from the deflected vector. Closure
/// errors propagate and unused reductions are never evaluated.
fn anchor_reduced_vector_via(
    astrometric: [f64; 3],
    noaberr: bool,
    nogdefl: bool,
    deflect: impl Fn([f64; 3]) -> Result<[f64; 3], Error>,
    aberrate: impl Fn([f64; 3]) -> Result<[f64; 3], Error>,
) -> Result<[f64; 3], Error> {
    if noaberr && nogdefl {
        return Ok(astrometric);
    }
    if nogdefl {
        return aberrate(astrometric);
    }
    let apparent = aberrate(deflect(astrometric)?)?;
    if noaberr {
        let aberrated = aberrate(astrometric)?;
        Ok(std::array::from_fn(|i| {
            apparent[i] - (aberrated[i] - astrometric[i])
        }))
    } else {
        Ok(apparent)
    }
}

impl<'a, P: crate::positions::StateProvider> SealedAnchorProvider<'a, P> {
    /// Borrow a sealed geometric-state provider for anchor evaluation.
    pub fn new(provider: &'a P) -> Self {
        Self {
            provider,
            native_state_at: None,
            native_ray_at: None,
            range_error_at: None,
            rows: None,
            acquired_ray_at: None,
            dynamic_override: None,
        }
    }

    /// Bind strict native acquisition and original-request range promotion.
    /// Both callbacks are borrowed for this provider's lifetime. Earth and
    /// deflector reads use the same selected native snapshot at their actual
    /// sample TT, including every interpolation node and closest approach.
    /// A consumed failure stops acquisition; raw range failures are wrapped
    /// once at the outer request TT. Geometric or deflection-free reductions
    /// do not read unused deflectors. With no callbacks, the existing sealed
    /// provider and optional LEB deflector policy are retained unchanged.
    pub fn with_native_acquisition(
        mut self,
        state_at: &'a AnchorStateAt<'a>,
        range_error_at: &'a AnchorRangeErrorAt<'a>,
    ) -> Self {
        self.native_state_at = Some(state_at);
        self.native_ray_at = None;
        self.range_error_at = Some(range_error_at);
        self
    }

    /// Bind a complete native ICRS anchor ray and outer range promotion.
    /// Both callbacks borrow the selected immutable provider for this value's
    /// lifetime. The ray callback owns fresh scalar TT/source-Time acquisition,
    /// catalog propagation and strict apparent reductions; this provider keeps
    /// its mean Vondrák frame, exact-date or three-node interpolation, and
    /// original-request error wrapper. This replaces a previously bound state
    /// callback; without either callback the stored optional-deflector policy applies.
    pub fn with_native_ray(
        mut self,
        ray_at: &'a AnchorRayAt<'a>,
        range_error_at: &'a AnchorRangeErrorAt<'a>,
    ) -> Self {
        self.native_ray_at = Some(ray_at);
        self.native_state_at = None;
        self.range_error_at = Some(range_error_at);
        self
    }

    /// Use actual defining STARS rows and their cache/constructor callbacks.
    /// Selection remains lazy until a dynamic anchor is consumed. No row or
    /// field is read here, and no callback failure falls back to baked data.
    pub fn with_live_rows(mut self, rows: &'a dyn AnchorRowAcquisition) -> Self {
        self.rows = Some(rows);
        self
    }

    /// Delegate dynamic anchors to an independently owned computation provider.
    ///
    /// No row, state or clock is acquired here. The supplied provider owns its
    /// live rows, cache, Time/main ordering and range errors. Its refusal is
    /// terminal; this provider's static rows and states are never substituted.
    pub fn with_dynamic_override(
        mut self,
        provider: &'a dyn crate::sidereal::DynamicAnchorProvider,
    ) -> Self {
        self.dynamic_override = Some(provider);
        self
    }

    /// Bind complete source-Time rays for already constructed live anchors.
    /// The callback and range mapper borrow the same immutable resource owner.
    /// The helper role preserves ephemeris-before-Time versus Time-before-
    /// ephemeris acquisition. Static callback paths remain separately owned.
    pub fn with_acquired_ray(
        mut self,
        ray_at: &'a AcquiredAnchorRayAt<'a>,
        range_error_at: &'a AnchorRangeErrorAt<'a>,
    ) -> Self {
        self.acquired_ray_at = Some(ray_at);
        self.range_error_at = Some(range_error_at);
        self
    }

    /// Acquire the geocentric observer at a finite scalar TT Julian day.
    /// A selected native state callback reads Earth from its bound immutable
    /// snapshot; otherwise the stored geometric provider supplies Earth.
    /// Both return barycentric ICRS position in AU and velocity in AU/day,
    /// admitted by StarObserver's finite-state checks. Consumed acquisition
    /// and observer validation errors propagate without a substitute state.
    fn earth_observer(&self, jd_tt: f64) -> Result<StarObserver, Error> {
        if let Some(state_at) = self.native_state_at {
            let (position, velocity) = state_at(crate::constants::EARTH, jd_tt)?;
            return StarObserver::new(position, velocity);
        }
        let state = self
            .provider
            .geometric_state(crate::constants::EARTH, jd_tt)?;
        StarObserver::new(state.pos(), state.vel())
    }

    /// Acquire one stored-provider deflector position at scalar TT, in ICRS AU.
    /// This closure propagates the geometric-state error unchanged. Its enclosing
    /// stored apparent-reduction lookup owns the optional-deflector policy;
    /// strict native ray/state callbacks bypass that policy. No state is cached
    /// or replaced by this helper.
    fn deflect_at(&self, id: i32, jd: f64) -> Result<[f64; 3], Error> {
        Ok(self.provider.geometric_state(id, jd)?.pos())
    }

    /// Resolve a built-in stellar or galactic anchor by its catalog key.
    /// The owned row carries J2000 direction, proper motion, parallax and radial
    /// velocity in the units of [`StarData`]. An absent dispatched key is a typed
    /// configuration fault, never a zero direction or an invented anchor.
    fn anchor_row(&self, key: &str) -> Result<StarData, Error> {
        crate::sidereal::anchor_star(key).ok_or(Error::ConfigurationMissing {
            missing_config: "anchor_star",
            suggestion: "baked ANCHOR_STARS row missing for a dispatched anchor",
        })
    }

    /// Mean-ecliptic longitude in degrees of a catalog anchor at TT.
    ///
    /// Deflection-free samples use the 8-day quadratic; active deflection
    /// uses exact-date evaluation. Callers restrict toggles to responsive
    /// modes and force both reductions off for the geometric pole. Each
    /// sampled epoch obtains its own observer, propagation and mean frame.
    fn anchor_ecliptic_lon(
        &self,
        key: &str,
        jd_tt: f64,
        noaberr: bool,
        nogdefl: bool,
    ) -> Result<f64, Error> {
        if let Some(rows) = self.rows {
            return self.live_anchor_ecliptic_lon(key, jd_tt, noaberr, nogdefl, false, rows);
        }
        let data = self.anchor_row(key)?;
        anchor_longitude_sampled_via(jd_tt, nogdefl, |tt| {
            let geo = self.anchor_icrs_at(&data, tt, noaberr, nogdefl)?;
            let frame = anchor_mean_frame(tt)?;
            let eq = crate::positions::mat3_vec3(frame.pn(), geo)?;
            let eps = frame.eps_true_rad();
            let ce = crate::positions::reference_frames::scalar_cos(eps);
            let se = crate::positions::reference_frames::scalar_sin(eps);
            let ye = eq[1] * ce + eq[2] * se;
            Ok(ye.atan2(eq[0]).to_degrees().rem_euclid(360.0))
        })
    }

    /// Hold one defining row while each node rereads its actual cache operands.
    /// Hits return before construction/Time/provider/frame work; misses retain
    /// original keyword units, calculate, then publish only a complete result.
    fn live_anchor_ecliptic_lon(
        &self,
        key: &str,
        jd_tt: f64,
        noaberr: bool,
        nogdefl: bool,
        geometric: bool,
        rows: &dyn AnchorRowAcquisition,
    ) -> Result<f64, Error> {
        let row = rows.acquire_row(key)?;
        let options = AnchorSampleOptions {
            eps_true: 0.0,
            nonut: true,
            geometric,
            noaberr,
            nogdefl,
        };
        anchor_longitude_sampled_via(jd_tt, geometric || nogdefl, |tt| {
            let sample = match row.sample(tt, options)? {
                AnchorSample::Cached(longitude) => return Ok(longitude),
                AnchorSample::Pending(sample) => sample,
            };
            let inputs = sample.native_inputs()?;
            let initial = star_initial_state_from_native(inputs);
            let geo = self.acquired_anchor_icrs_at(
                initial,
                tt,
                geometric || noaberr,
                geometric || nogdefl,
                AnchorRayRole::EclipticLongitude,
            )?;
            let frame = anchor_mean_frame(tt)?;
            let eq = crate::positions::mat3_vec3(frame.pn(), geo)?;
            let eps = frame.eps_true_rad();
            let ce = crate::positions::reference_frames::scalar_cos(eps);
            let se = crate::positions::reference_frames::scalar_sin(eps);
            let ye = eq[1] * ce + eq[2] * se;
            let longitude = ye.atan2(eq[0]).to_degrees().rem_euclid(360.0);
            sample.publish(longitude)?;
            Ok(longitude)
        })
    }

    /// Observe an owned native constructor without a degree/hour roundtrip.
    /// Complete ray callbacks retain their actual source-Time acquisition.
    /// The independent stored default uses the same geometric states and
    /// established optional-deflector policy as its static provider.
    fn acquired_anchor_icrs_at(
        &self,
        (initial, velocity): ([f64; 3], [f64; 3]),
        jd_tt: f64,
        noaberr: bool,
        nogdefl: bool,
        role: AnchorRayRole,
    ) -> Result<[f64; 3], Error> {
        if let Some(ray_at) = self.acquired_ray_at {
            return ray_at((initial, velocity), jd_tt, noaberr, nogdefl, role);
        }
        let observer = self.earth_observer(jd_tt)?;
        let norm = crate::positions::vec3_dist(initial);
        let unit: [f64; 3] = std::array::from_fn(|i| initial[i] / (norm + f64::MIN_POSITIVE));
        let projection = (unit[0] * observer.pos_au()[0]
            + unit[1] * observer.pos_au()[1]
            + unit[2] * observer.pos_au()[2])
            / crate::positions::FRAME_LIGHT_SPEED_AU_DAY;
        let hi = jd_tt.floor();
        let time = crate::types::TwoPartEpoch { hi, lo: jd_tt - hi };
        let tdb = time.hi
            + (time.lo + crate::time::source_tdb_minus_tt_seconds(time.hi, time.lo) / 86400.0);
        let span = tdb + projection - crate::constants::J2000;
        let astro = std::array::from_fn(|i| initial[i] + velocity[i] * span - observer.pos_au()[i]);
        let light_time =
            crate::positions::vec3_dist(astro) / crate::positions::FRAME_LIGHT_SPEED_AU_DAY;
        self.reduce_anchor_ray(&observer, astro, jd_tt, light_time, noaberr, nogdefl)
    }

    /// Catalog anchor's reduced observer-relative ICRS vector in AU at TT.
    ///
    /// The catalog state propagates to TDB plus the observer's light-time
    /// projection on the J2000 direction. The same Earth velocity and
    /// astrometric light time drive all aberration reductions, including
    /// the NOABERR post-correction. Provider and reduction faults propagate.
    fn anchor_icrs_at(
        &self,
        data: &StarData,
        jd_tt: f64,
        noaberr: bool,
        nogdefl: bool,
    ) -> Result<[f64; 3], Error> {
        if let Some(ray_at) = self.native_ray_at {
            return ray_at(data, jd_tt, noaberr, nogdefl);
        }
        use crate::positions::FRAME_LIGHT_SPEED_AU_DAY;
        let observer = self.earth_observer(jd_tt)?;
        let epoch = topo_star_epoch_days(data, jd_tt, observer.pos_au())?;
        let (pos, _) = star_icrs_state_catalog_angles(data, epoch);
        let astro: [f64; 3] = std::array::from_fn(|i| pos[i] - observer.pos_au()[i]);
        let distance = (astro[0] * astro[0] + astro[1] * astro[1] + astro[2] * astro[2]).sqrt();
        let lt = distance / FRAME_LIGHT_SPEED_AU_DAY;
        self.reduce_anchor_ray(&observer, astro, jd_tt, lt, noaberr, nogdefl)
    }

    /// Reduce an already constructed observer-relative anchor ray in ICRS AU.
    ///
    /// `jd_tt` selects deflector states and `light_time` is the astrometric
    /// travel time in days. Native state callbacks require every consumed
    /// deflector; stored reads retain their optional-deflector policy. Both use
    /// the acquired observer's velocity and the same astrometric light time for
    /// all aberration evaluations, including the NOABERR displacement. This
    /// helper acquires no observer or clock and does not wrap provider errors.
    fn reduce_anchor_ray(
        &self,
        observer: &StarObserver,
        astrometric: [f64; 3],
        jd_tt: f64,
        light_time: f64,
        noaberr: bool,
        nogdefl: bool,
    ) -> Result<[f64; 3], Error> {
        use crate::positions::{
            FRAME_LIGHT_SPEED_AU_DAY, apply_aberration_with_c, apply_gravitational_deflection,
        };
        anchor_reduced_vector_via(
            astrometric,
            noaberr,
            nogdefl,
            |vector| {
                if let Some(state_at) = self.native_state_at {
                    return crate::positions::apply_gravitational_deflection_with_lookup(
                        vector,
                        observer.pos_au(),
                        jd_tt,
                        light_time,
                        |id, sample_tt| state_at(id, sample_tt).map(|(position, _)| Some(position)),
                    );
                }
                apply_gravitational_deflection(
                    vector,
                    observer.pos_au(),
                    jd_tt,
                    light_time,
                    |id, jd| self.deflect_at(id, jd),
                )
            },
            |vector| {
                apply_aberration_with_c(
                    vector,
                    observer.vel_au_day(),
                    light_time,
                    FRAME_LIGHT_SPEED_AU_DAY,
                )
            },
        )
    }

    /// Preserve the consumed raw node/deflector identity until public wrapping.
    /// Only SpkOutOfRange invokes the caller's main-image range wrapper;
    /// `request_tt` is the outer scalar anchor epoch, never the sampled node
    /// or retarded deflector epoch. Other errors retain their typed identity.
    fn public_anchor_error(&self, request_tt: f64, error: Error) -> Error {
        if matches!(&error, Error::SpkOutOfRange { .. }) {
            if let Some(wrap) = self.range_error_at {
                return wrap(request_tt, error);
            }
        }
        error
    }

    /// Apparent Sgr A* right ascension on the mean equator of date, in radians.
    /// Apply the propagation and selected aberration/deflection reductions of
    /// [`Self::anchor_icrs_at`], then biased Vondrak precession and atan 2(y,x).
    /// With both reductions disabled the ray is geometric. Consumed catalog,
    /// provider and frame errors propagate to the enclosing range wrapper.
    fn galcent_ra_of_date(&self, jd_tt: f64, noaberr: bool, nogdefl: bool) -> Result<f64, Error> {
        use crate::positions::{mat3_vec3, vondrak_precession_matrix};
        let row = self
            .rows
            .map(|rows| rows.acquire_row("GAL_CENTER"))
            .transpose()?;
        let geo = if let Some(row) = &row {
            let inputs = row.galactic_inputs()?;
            self.acquired_anchor_icrs_at(
                star_initial_state_from_native(inputs),
                jd_tt,
                noaberr,
                nogdefl,
                AnchorRayRole::GalacticRightAscension,
            )?
        } else {
            let data = self.anchor_row("GAL_CENTER")?;
            self.anchor_icrs_at(&data, jd_tt, noaberr, nogdefl)?
        };
        let pn = vondrak_precession_matrix(jd_tt, true)?;
        let eq = mat3_vec3(pn, geo)?;
        Ok(eq[1].atan2(eq[0]))
    }

    /// Sgr A* hour-circle projection onto the mean ecliptic of date, in degrees.
    /// For right ascension alpha and mean obliquity eps (both radians), evaluate
    /// atan 2(sin(alpha), cos(alpha)*cos(eps)), convert to degrees and reduce to
    /// [0,360). The observer, apparent corrections and mean frame are those of
    /// [`Self::galcent_ra_of_date`]; no substitute angle is supplied on failure.
    fn mula_wilhelm_lon(&self, jd_tt: f64, noaberr: bool, nogdefl: bool) -> Result<f64, Error> {
        let alpha = self.galcent_ra_of_date(jd_tt, noaberr, nogdefl)?;
        let eps = crate::positions::vondrak_mean_obliquity_rad(jd_tt)?;
        Ok((alpha.sin().atan2(alpha.cos() * eps.cos()))
            .to_degrees()
            .rem_euclid(360.0))
    }
}

impl<P: crate::positions::StateProvider> crate::sidereal::DynamicAnchorProvider
    for SealedAnchorProvider<'_, P>
{
    fn has_analytic_anchor_rates(&self) -> bool {
        if let Some(provider) = self.dynamic_override {
            return provider.has_analytic_anchor_rates();
        }
        self.rows.is_none()
            && self.native_ray_at.is_none()
            && self.native_state_at.is_none()
            && self.provider.has_analytic_kinematics()
    }
    fn anchor_rate_deg_day(
        &self,
        anchor: crate::sidereal::DynamicAnchor,
        jd_tt: f64,
        noaberr: bool,
        nogdefl: bool,
    ) -> Result<Option<f64>, Error> {
        use crate::sidereal::DynamicAnchor as A;
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        if let Some(provider) = self.dynamic_override {
            return provider.anchor_rate_deg_day(anchor, jd_tt, noaberr, nogdefl);
        }
        if self.rows.is_some()
            || self.native_ray_at.is_some()
            || self.native_state_at.is_some()
            || !self.provider.has_analytic_kinematics()
        {
            return Ok(None);
        }
        let (key, noaberr, nogdefl, wilhelm) = match anchor {
            A::Spica => ("SPICA", noaberr, nogdefl, false),
            A::Revati => ("REVATI", noaberr, nogdefl, false),
            A::Pushya => ("PUSHYA", noaberr, nogdefl, false),
            A::Mula => ("MULA", noaberr, nogdefl, false),
            A::Aldebaran => ("ALDEBARAN", noaberr, nogdefl, false),
            A::GalCenter => ("GAL_CENTER", noaberr, nogdefl, false),
            A::GalNorthPole => ("GAL_NORTH_POLE", true, true, false),
            A::MulaWilhelm => ("GAL_CENTER", noaberr, nogdefl, true),
        };
        let data = self.anchor_row(key)?;
        let at = |tt| {
            crate::positions::static_anchor_motion(
                self.provider,
                &data,
                tt,
                noaberr,
                nogdefl,
                wilhelm,
            )
        };
        let result = if nogdefl && !wilhelm {
            // Differentiate the actual three fixed eight-day node polynomial.
            // This is an interpolation law, not a widened derivative stencil.
            let base = (jd_tt / 16.0).floor() * 16.0;
            let a = at(base)?.0;
            let b = at(base + 8.0)?.0;
            let c = at(base + 16.0)?.0;
            let unwrap = |v: f64| {
                if v - b > 180.0 {
                    v - 360.0
                } else if v - b < -180.0 {
                    v + 360.0
                } else {
                    v
                }
            };
            let x = (jd_tt - base) / 8.0;
            Ok((unwrap(a) * (x - 1.5) - b * (2.0 * x - 2.0) + unwrap(c) * (x - 0.5)) / 8.0)
        } else {
            at(jd_tt).map(|(_, rate)| rate)
        };
        result
            .map(Some)
            .map_err(|error| self.public_anchor_error(jd_tt, error))
    }

    /// Mean-ecliptic anchor longitude at scalar TT, in degrees.
    /// The sidereal dispatcher owns target-angle subtraction and normalization.
    /// Star/center anchors retain the requested correction toggles; the north
    /// galactic pole is geometric, and Wilhelm projects Sgr A*'s hour circle.
    /// Nonfinite epochs refuse. Consumed raw SPK ranges are wrapped only at the
    /// original anchor request TT; other provider/frame errors retain identity.
    fn anchor_longitude_deg(
        &self,
        anchor: crate::sidereal::DynamicAnchor,
        jd_tt: f64,
        noaberr: bool,
        nogdefl: bool,
    ) -> Result<f64, Error> {
        use crate::sidereal::DynamicAnchor as A;
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        if let Some(provider) = self.dynamic_override {
            return provider.anchor_longitude_deg(anchor, jd_tt, noaberr, nogdefl);
        }
        let result = match anchor {
            A::Spica => self.anchor_ecliptic_lon("SPICA", jd_tt, noaberr, nogdefl),
            A::Revati => self.anchor_ecliptic_lon("REVATI", jd_tt, noaberr, nogdefl),
            A::Pushya => self.anchor_ecliptic_lon("PUSHYA", jd_tt, noaberr, nogdefl),
            A::Mula => self.anchor_ecliptic_lon("MULA", jd_tt, noaberr, nogdefl),
            A::Aldebaran => self.anchor_ecliptic_lon("ALDEBARAN", jd_tt, noaberr, nogdefl),
            A::GalCenter => self.anchor_ecliptic_lon("GAL_CENTER", jd_tt, noaberr, nogdefl),
            A::GalNorthPole => {
                if let Some(rows) = self.rows {
                    self.live_anchor_ecliptic_lon("GAL_NORTH_POLE", jd_tt, false, false, true, rows)
                } else {
                    self.anchor_ecliptic_lon("GAL_NORTH_POLE", jd_tt, true, true)
                }
            }
            A::MulaWilhelm => self.mula_wilhelm_lon(jd_tt, noaberr, nogdefl),
        };
        result.map_err(|error| self.public_anchor_error(jd_tt, error))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn anchor_quadratic_unwraps_periodic_nodes_in_source_order() {
        use std::cell::RefCell;
        let seen = RefCell::new(Vec::new());
        // A periodic synthetic quadratic crosses zero between the nodes.
        let longitude = anchor_longitude_sampled_via(4.0, true, |tt| {
            seen.borrow_mut().push(tt);
            let x = tt / 8.0;
            Ok((358.0 + 3.0 * x + x * x).rem_euclid(360.0))
        })
        .unwrap();
        assert_eq!(seen.into_inner(), [0.0, 8.0, 16.0]);
        assert_eq!(longitude, 359.75);
        for tt in [0.0, 8.0, 16.0] {
            let at = |epoch: f64| Ok((358.0 + 0.5 * epoch).rem_euclid(360.0));
            assert_eq!(
                anchor_longitude_sampled_via(tt, true, at).unwrap(),
                at(tt).unwrap()
            );
        }
    }

    #[test]
    // Keep the fixture Lagrange-weight arithmetic and its operation order unchanged.
    #[allow(clippy::neg_multiply)]
    fn anchor_deflection_free_sinusoid_uses_quadratic_not_exact_date() {
        let at = |tt: f64| Ok((359.0 + 20.0 * (tt / 20.0).sin()).rem_euclid(360.0));
        let value = anchor_longitude_sampled_via(4.0, true, at).unwrap();
        // The Lagrange weights at x=1/2 are 3/8, 3/4, -1/8.
        let expected = (-1.0 * 3.0 / 8.0 + (20.0 * (8.0_f64 / 20.0).sin() - 1.0) * 3.0 / 4.0
            - (20.0 * (16.0_f64 / 20.0).sin() - 1.0) / 8.0)
            .rem_euclid(360.0);
        assert!((value - expected).abs() < 1e-12);
        assert!((value - at(4.0).unwrap()).abs() > 0.01);
    }

    #[test]
    fn anchor_active_deflection_is_exact_and_nodes_propagate_coverage() {
        use std::cell::RefCell;
        let seen = RefCell::new(Vec::new());
        let sampled = |tt| {
            seen.borrow_mut().push(tt);
            if tt > 8.0 {
                Err(Error::EphemerisRange {
                    requested_jd: tt,
                    start_jd: 0.0,
                    end_jd: 8.0,
                    body: None,
                })
            } else {
                Ok(42.0)
            }
        };
        assert_eq!(
            anchor_longitude_sampled_via(4.0, false, sampled).unwrap(),
            42.0
        );
        assert_eq!(*seen.borrow(), [4.0]);
        seen.borrow_mut().clear();
        // Even an endpoint reads the complete source three-node stencil.
        assert!(matches!(anchor_longitude_sampled_via(0.0, true, sampled),
            Err(Error::EphemerisRange { requested_jd, .. }) if requested_jd == 16.0));
        assert_eq!(*seen.borrow(), [0.0, 8.0, 16.0]);
        assert!(matches!(
            anchor_longitude_sampled_via(f64::NAN, true, |_| {
                panic!("non-finite epochs must refuse before sampling")
            }),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[test]
    fn anchor_flags_use_source_apparent_minus_astrometric_aberration() {
        use crate::constants::{FLG_BARYCTR, FLG_HELCTR, FLG_NOABERR, FLG_NOGDEFL, FLG_TRUEPOS};
        use std::cell::RefCell;
        let astro = [1.0, 2.0, 3.0];
        for (flags, expected, operations) in [
            (0, [14.0, 28.0, 42.0], vec!["deflect", "aberrate"]),
            (
                FLG_NOABERR,
                [3.0, 6.0, 9.0],
                vec!["deflect", "aberrate", "aberrate"],
            ),
            (FLG_NOGDEFL, [12.0, 24.0, 36.0], vec!["aberrate"]),
            (FLG_TRUEPOS, astro, vec![]),
            (FLG_NOABERR | FLG_NOGDEFL, astro, vec![]),
            (FLG_HELCTR, astro, vec![]),
            (FLG_BARYCTR, astro, vec![]),
        ] {
            let seen = RefCell::new(Vec::new());
            let (noaberr, nogdefl, _) = crate::sidereal::resolve_ayanamsha_flags(flags);
            let reduced = anchor_reduced_vector_via(
                astro,
                noaberr,
                nogdefl,
                |v| {
                    seen.borrow_mut().push("deflect");
                    Ok([v[0] + 1.0, v[1] + 2.0, v[2] + 3.0])
                },
                |v| {
                    seen.borrow_mut().push("aberrate");
                    Ok([2.0 * v[0] + 10.0, 2.0 * v[1] + 20.0, 2.0 * v[2] + 30.0])
                },
            )
            .unwrap();
            assert_eq!(reduced, expected, "flags={flags}");
            assert_eq!(seen.into_inner(), operations, "flags={flags}");
        }
    }
}
