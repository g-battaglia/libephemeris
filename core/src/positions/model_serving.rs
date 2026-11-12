//! Provider-backed analytical-model acquisition.
//!
//! Model elements remain independent of the support ephemeris. The provider
//! supplies the Sun/Earth and deflector states, while the analytical branches
//! retain their own light-time stencils, frames and output reduction.

use super::*;

/// Calculate a curated minor body or hypothetical model over stored support.
///
/// `request_ut1` and `global_tt` are independently retained UT1/TT observation
/// JDs. Configured `model_deltat` supplies days at UT1, including Waldemath's
/// defining anchor; `ayanamsha_at` supplies degrees for the selected sidereal
/// model. Sun/Earth/deflector states come only from the supplied provider,
/// in barycentric ICRS AU and AU/day.
///
/// Curated bodies use epoch-blended Keplerian elements. Their fallback echoes
/// TOPOCTR without adding site parallax; hypothetical model reductions can
/// consume a site in east-positive degrees and metres. Output follows the
/// requested spherical/Cartesian and degree/radian representation, with rates
/// per day. A caller applying representation shaping afterward must omit XYZ
/// and RADIANS here. Missing support, unavailable observers, model-domain,
/// coverage, corrupt-data and non-finite errors propagate. Element propagation
/// is an approximation whose physical error grows away from defining epochs.
#[allow(clippy::too_many_arguments)]
pub fn serving_model_position<P: StateProvider>(
    provider: &P,
    request_ut1: f64,
    global_tt: f64,
    body: i32,
    flags: i32,
    site: Option<(f64, f64, f64)>,
    model_deltat: &dyn Fn(f64) -> Result<f64>,
    ayanamsha_at: &dyn Fn(f64, i32) -> Result<f64>,
    _sidereal_mode: i32,
) -> Result<(PositionResult, i32)> {
    for value in [request_ut1, global_tt] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    if crate::minor_bodies::minor_body_elements(body).is_some() {
        let mut state = curated_keplerian_state(provider, global_tt, body, flags)?;
        if flags & FLG_SIDEREAL != 0 && flags & FLG_EQUATORIAL == 0 {
            state = horizons_output::local_sidereal(state, request_ut1, flags, ayanamsha_at)?;
        }
        state = horizons_output::local_convert(state, global_tt, flags, false)?;
        return Ok((apply_output_flags(state, flags), flags));
    }
    if !(CUPIDO..=WALDEMATH).contains(&body) || body == NIBIRU {
        return Err(Error::UnknownBody { body_id: body });
    }
    let state_at = |id, tt| icrs_state(provider, id, tt);
    let earth_helio_at = |tt| observed_earth_minus_sun(&state_at, tt);
    let uranian_earth_at = |tt| leb_earth_minus_sun(&state_at, tt);
    let already_requested_frame = (CUPIDO..=ISIS).contains(&body);
    let mut rate_origin = ModelRateOrigin::LegacyStencil;
    let state = if (CUPIDO..=POSEIDON).contains(&body) {
        uranian_center_state(&state_at, &uranian_earth_at, global_tt, body, flags, false)?
    } else if body == ISIS {
        isis_center_state(&state_at, &earth_helio_at, global_tt, flags)?
    } else if matches!(body, WHITE_MOON | WALDEMATH) {
        let symbolic = if body == WALDEMATH {
            let anchor_ut = crate::hypothetical::WALDEMATH_ANCHOR_JD_UT;
            let anchor_tt = anchor_ut + model_deltat(anchor_ut)?;
            let mean_lon = crate::hypothetical::WALDEMATH_ANCHOR_APPARENT_LON_DEG
                - crate::erfa::nut06a(anchor_tt)?.0.to_degrees();
            crate::hypothetical::waldemath_state_with_anchor(global_tt, anchor_tt, mean_lon)?
        } else {
            crate::hypothetical::calc_white_moon_position(global_tt, false)?
        };
        symbolic_center_state(&state_at, &earth_helio_at, global_tt, flags, symbolic)?
    } else if native_predicted_rates_available(provider, body, flags) {
        let (state, origin) = native_provider_predicted_state(provider, global_tt, body, flags)?;
        rate_origin = origin;
        state
    } else {
        predicted_center_state(
            &state_at,
            &|tt| observed_earth_with_origin(&state_at, tt),
            true,
            global_tt,
            body,
            flags,
        )?
    };
    let state = finish_model_state_with_rate_origin(
        &state_at,
        state,
        &|_| Ok(request_ut1),
        global_tt,
        flags,
        already_requested_frame,
        site,
        ayanamsha_at,
        rate_origin,
    )?;
    Ok((apply_output_flags(state, flags), flags))
}

/// Calculate a Uranian or Isis model using explicit support callbacks.
///
/// `tt` is the observation TT JD. A lazy `ut_at` supplies UT1 only for ecliptic
/// sidereal output. `state_at` returns finite SSB ICRS AU/AU-day states;
/// `earth_helio_at` returns Earth's retarded ICRS AU ray from a Sun fixed at
/// reception. Models use J2000 ecliptic elements, their own retardation and
/// center-specific rate stencil, then the shared nutation/site/sidereal tail.
///
/// Site coordinates are east longitude/latitude in degrees and height metres.
/// Output is spherical degrees/AU with per-day rates; external entries apply
/// representation shaping. Unsupported identifiers, invalid clocks, missing
/// support or site, and frame/model/coverage errors propagate without provider
/// substitution. Physical accuracy remains limited by the analytical model.
#[allow(clippy::too_many_arguments)]
#[cfg(test)]
pub(super) fn classic_model_state(
    state_at: &dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
    earth_helio_at: &dyn Fn(f64) -> Result<[f64; 3]>,
    ut_at: &dyn Fn(f64) -> Result<f64>,
    tt: f64,
    body: i32,
    flags: i32,
    site: Option<(f64, f64, f64)>,
    aya: &dyn Fn(f64, i32) -> Result<f64>,
) -> Result<PositionResult> {
    classic_model_state_with_models(
        None,
        state_at,
        earth_helio_at,
        ut_at,
        tt,
        body,
        flags,
        site,
        aya,
    )
}

/// Evaluate a classical analytical model through an already validated
/// native observation entry. The callbacks supply Sun/Earth support and
/// preserve their typed coverage errors. Output uses degrees, AU and per-day
/// rates with the same center and frame reduction as [`classic_model_state`].
#[allow(clippy::too_many_arguments)]
#[cfg(test)]
pub(super) fn classic_model_state_source_time(
    state_at: &dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
    earth_helio_at: &dyn Fn(f64) -> Result<[f64; 3]>,
    ut_at: &dyn Fn(f64) -> Result<f64>,
    tt: f64,
    body: i32,
    flags: i32,
    site: Option<(f64, f64, f64)>,
    aya: &dyn Fn(f64, i32) -> Result<f64>,
) -> Result<PositionResult> {
    classic_model_state_source_time_with_models(
        None,
        state_at,
        earth_helio_at,
        ut_at,
        tt,
        body,
        flags,
        site,
        aya,
    )
}

/// Curated Keplerian fallback in ecliptic degrees/AU and per-day rates.
/// The model's heliocentric place is independent of the provider; geocentric
/// and BARYCTR requests resolve Earth then Sun aliases before observing the
/// retarded Earth through scalar-TT ICRS support. Frame TT/TDB belong to the
/// independent scalar epoch object, not configured fast TT or stored nutation.
/// Missing aliases and support ranges retain their public body identity;
/// corrupt central/retarded states propagate without another fallback.
fn curated_keplerian_state<P: StateProvider>(
    provider: &P,
    tt: f64,
    body: i32,
    flags: i32,
) -> Result<PositionResult> {
    let earth_at = |jd: f64| {
        // Alias resolution precedes Sun-first vector evaluation. In
        // particular, a missing Earth wins over a missing/outside Sun.
        for endpoint in [EARTH, SUN] {
            if !provider.has_body(endpoint) {
                return Err(Error::UnknownBody { body_id: endpoint });
            }
        }
        let state_at = |endpoint, epoch| {
            icrs_state(provider, endpoint, epoch).map_err(|error| {
                escalate_sealed_provider_miss(
                    ProviderMode::SealedLeb,
                    provider,
                    endpoint,
                    epoch,
                    error,
                )
            })
        };
        let relative = observed_earth_minus_sun(&state_at, jd)?;
        let whole = jd.floor();
        let fraction = jd - whole;
        let tdb =
            whole + (fraction + reference_frames::tdb_minus_tt_seconds(whole, fraction)? / 86400.0);
        true_ecliptic_vector_at(jd, tdb, relative)
    };
    keplerian_fallback_with_earth(tt, body, flags, &earth_at)
}

/// Acquire one scalar-TT barycentric ICRS state, retaining the provider's
/// owned entry selection, center composition and scalar range checks.
fn icrs_state<P: StateProvider>(provider: &P, body: i32, tt: f64) -> Result<([f64; 3], [f64; 3])> {
    serving_channel_gate(provider, body)?;
    let state = provider.geometric_state(body, tt)?;
    Ok((state.pos(), state.vel()))
}

/// Euclidean norm in AU, with the model branches' left-associated sum.
pub(in crate::positions) fn model_norm(p: [f64; 3]) -> f64 {
    (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt()
}

/// Model spherical degrees/AU to ecliptic Cartesian AU.
pub(in crate::positions) fn model_xyz(lon: f64, lat: f64, radius: f64) -> [f64; 3] {
    let (lr, br) = (lon.to_radians(), lat.to_radians());
    let cb = br.cos();
    [
        radius * cb * lr.cos(),
        radius * cb * lr.sin(),
        radius * br.sin(),
    ]
}

/// Ecliptic Cartesian AU to longitude [0,360), latitude degrees and AU.
/// A zero radius has zero latitude; the longitude retains atan2 semantics.
pub(in crate::positions) fn model_spherical(p: [f64; 3]) -> [f64; 3] {
    let radius = model_norm(p);
    [
        p[1].atan2(p[0]).to_degrees().rem_euclid(360.0),
        if radius > 0.0 {
            (p[2] / radius).clamp(-1.0, 1.0).asin().to_degrees()
        } else {
            0.0
        },
        radius,
    ]
}

/// Centered spherical state. `step` is days; `wrapped` selects the fast-body
/// modulo rule rather than the slow-model post-division 0/360 rescue.
pub(in crate::positions) fn model_stencil(
    center: [f64; 3],
    prev: [f64; 3],
    next: [f64; 3],
    step: f64,
    wrapped: bool,
) -> PositionResult {
    let denominator = 2.0 * step;
    let mut dlon = if wrapped {
        ((next[0] - prev[0] + 180.0).rem_euclid(360.0) - 180.0) / denominator
    } else {
        (next[0] - prev[0]) / denominator
    };
    if !wrapped {
        if dlon > 180.0 / denominator {
            dlon -= 360.0 / denominator;
        } else if dlon < -180.0 / denominator {
            dlon += 360.0 / denominator;
        }
    }
    [
        center[0],
        center[1],
        center[2],
        dlon,
        (next[1] - prev[1]) / denominator,
        (next[2] - prev[2]) / denominator,
    ]
}

/// Angular derivative on the already acquired native support domain.
/// `None` rejects a side sample before it reads a different support record.
fn native_prediction_angular_rates(
    sample: &dyn Fn(f64) -> Result<Option<[f64; 2]>>,
    support: (f64, f64),
) -> Result<[f64; 2]> {
    super::super::derivative::angular_rates_on_owned_domain(
        sample,
        crate::orbital::NODAPS_SPEED_STEP_DAYS,
        support,
    )
}

/// Spherical velocity of the same apparent place, in degrees/day and AU/day.
/// Radius uses the Cartesian derivative's relative resolution, without the
/// angular near-zero floor. All components retain the same owned support.
pub(super) fn native_prediction_rates(
    sample: &dyn Fn(f64) -> Result<Option<[f64; 3]>>,
    support: (f64, f64),
) -> Result<[f64; 3]> {
    let angular = native_prediction_angular_rates(
        &|offset| sample(offset).map(|value| value.map(|value| [value[0], value[1]])),
        support,
    )?;
    let radial = super::super::derivative::cartesian_rates_on_owned_domain(
        &|offset| sample(offset).map(|value| value.map(|value| [value[2], 0.0, 0.0])),
        crate::orbital::NODAPS_SPEED_STEP_DAYS,
        support,
    )?;
    Ok([angular[0], angular[1], radial[0]])
}

/// Optional spherical-rate refinement after the original three places succeed.
pub(super) type PredictedRateAt<'a> =
    dyn Fn(&dyn Fn(f64) -> Result<[f64; 3]>, f64, f64) -> Result<Option<[f64; 3]>> + 'a;

/// Velocity method actually consumed by a model state, independent of its id.
/// Only a successful owned-place derivative admits an instantaneous frame
/// rate. Unrefined states retain their existing combined stencil convention.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(in crate::positions) enum ModelRateOrigin {
    LegacyStencil,
    OwnedDerivative,
}

/// Record successful refinement within this synchronous calculation only.
/// Eligibility, metadata availability and a declined `None` do not qualify a
/// state. The borrowed callback and every error keep their original order.
pub(super) fn with_model_rate_origin(
    rates_at: &PredictedRateAt<'_>,
    calculate: impl FnOnce(&PredictedRateAt<'_>) -> Result<PositionResult>,
) -> Result<(PositionResult, ModelRateOrigin)> {
    let origin = std::cell::Cell::new(ModelRateOrigin::LegacyStencil);
    let observed = |place: &dyn Fn(f64) -> Result<[f64; 3]>, epoch, width| {
        let rates = rates_at(place, epoch, width)?;
        origin.set(if rates.is_some() {
            ModelRateOrigin::OwnedDerivative
        } else {
            ModelRateOrigin::LegacyStencil
        });
        Ok(rates)
    };
    let state = calculate(&observed)?;
    Ok((state, origin.get()))
}

type PredictedStateAt<'a> = dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])> + 'a;
type PredictedDeflectorAt<'a> = dyn Fn(i32, f64) -> Result<[f64; 3]> + 'a;

const PREDICTED_SUPPORT_BODIES: [i32; 4] = [
    SUN,
    EARTH,
    crate::constants::JUPITER,
    crate::constants::SATURN,
];

/// A bounded ordered sequence of one built-in place's support roles only.
/// Overflow declines refinement; it never truncates or broadens admission.
pub(super) const PREDICTED_MAX_READ_ROLES: usize = 32;

/// Ordered central epochs only: no stored states, coefficient copies or cache.
/// The original three places complete before any optional rate work begins.
struct PredictedSupportReads {
    roles: [Option<(i32, f64)>; PREDICTED_MAX_READ_ROLES],
    count: usize,
    recording: bool,
    complete: bool,
}

impl PredictedSupportReads {
    fn new() -> Self {
        Self {
            roles: [None; PREDICTED_MAX_READ_ROLES],
            count: 0,
            recording: true,
            complete: true,
        }
    }

    fn note(&mut self, body: i32, epoch: f64, success: bool) {
        if !self.recording {
            return;
        }
        if !success
            || !epoch.is_finite()
            || !PREDICTED_SUPPORT_BODIES.contains(&body)
            || self.count == self.roles.len()
        {
            self.complete = false;
            return;
        }
        self.roles[self.count] = Some((body, epoch));
        self.count += 1;
    }

    /// Each actual central role owns its own record, even when another role
    /// reads the same body in a different record. Consult metadata only after
    /// all three original places, retaining their physical error ordering.
    fn owned_domains<P: StateProvider>(&self, provider: &P) -> Option<PredictedSupportDomains> {
        if !self.complete {
            return None;
        }
        let mut roles = [None; PREDICTED_MAX_READ_ROLES];
        for (index, role) in self.roles.iter().take(self.count).enumerate() {
            let (body, epoch) = (*role)?;
            let interval = provider.smooth_interval_tt(body, epoch).ok().flatten()?;
            if !(interval.0 < epoch && epoch < interval.1) {
                return None;
            }
            roles[index] = Some(PredictedSupportDomain {
                body,
                epoch,
                interval,
            });
        }
        Some(PredictedSupportDomains {
            roles,
            count: self.count,
        })
    }
}

#[derive(Clone, Copy)]
struct PredictedSupportDomain {
    body: i32,
    epoch: f64,
    interval: (f64, f64),
}

#[derive(Clone, Copy)]
struct PredictedSupportDomains {
    roles: [Option<PredictedSupportDomain>; PREDICTED_MAX_READ_ROLES],
    count: usize,
}

impl PredictedSupportDomains {
    fn contains(self, role: usize, body: i32, epoch: f64) -> bool {
        self.roles
            .get(role)
            .copied()
            .flatten()
            .is_some_and(|domain| {
                domain.body == body && domain.interval.0 < epoch && epoch < domain.interval.1
            })
    }

    fn direction_hint(self) -> (f64, f64) {
        self.roles
            .iter()
            .flatten()
            .fold((f64::INFINITY, f64::INFINITY), |(left, right), domain| {
                (
                    left.min(domain.epoch - domain.interval.0),
                    right.min(domain.interval.1 - domain.epoch),
                )
            })
    }
}

/// Built-in predicted models sharing the owned-domain apparent-place method.
/// The moving date frame of Proserpina is part of the sampled place itself.
pub(in crate::positions) fn is_builtin_predicted_rate_body(body: i32) -> bool {
    (HARRINGTON..=VULCAN).contains(&body) || body == PROSERPINA
}

fn native_predicted_rates_available<P: StateProvider>(provider: &P, body: i32, flags: i32) -> bool {
    flags & FLG_SPEED != 0
        && flags & (FLG_HELCTR | FLG_BARYCTR) == 0
        && is_builtin_predicted_rate_body(body)
        && provider.has_analytic_kinematics()
}

/// Refine native predicted-model spherical rates without changing their place.
/// Opaque providers and callback models use the original reduction. The old
/// center/previous/next reads and errors precede supplementary differentiation.
/// Every new physical read must remain inside an already consumed record;
/// unsupported or joined domains retain the original derivative method.
fn native_provider_predicted_state<P: StateProvider>(
    provider: &P,
    tt: f64,
    body: i32,
    flags: i32,
) -> Result<(PositionResult, ModelRateOrigin)> {
    let errors = super::PredictedSpkReads::new(true);
    with_owned_predicted_support(provider, &errors, |state_at, rates_at, central_complete| {
        with_model_rate_origin(rates_at, |rates_at| {
            predicted_center_state_source_time_with_models_and_rates(
                None,
                state_at,
                &|epoch| observed_earth_with_origin(state_at, epoch),
                true,
                tt,
                body,
                flags,
                Some(rates_at),
                Some(central_complete),
            )
        })
    })
}

/// Differentiate the same scalar serving recipe over actual acquired LEB and
/// optional SPK-offset domains, without claiming analytic provider velocities.
#[allow(clippy::too_many_arguments)]
pub(in crate::positions) fn owned_provider_predicted_position<P: StateProvider>(
    provider: &P,
    original_deflector_provider: &P,
    position_provider: &P,
    offsets: &super::PredictedSpkReads,
    request_ut1: f64,
    tt: f64,
    body: i32,
    flags: i32,
    site: Option<(f64, f64, f64)>,
    ayanamsha_at: &dyn Fn(f64, i32) -> Result<f64>,
) -> Result<(PositionResult, i32)> {
    with_owned_predicted_support_and_positions(
        provider,
        Some((original_deflector_provider, position_provider)),
        offsets,
        |state_at, rates_at, central_complete, deflector_at| {
            owned_predicted_model_position(
                state_at,
                deflector_at,
                rates_at,
                central_complete,
                request_ut1,
                tt,
                body,
                flags,
                site,
                ayanamsha_at,
            )
        },
    )
}

/// Share source admission and first-error retention for the two owned routes.
/// The ledger retains no SPK query when only LEB records are involved.
fn with_owned_predicted_support<P: StateProvider, T>(
    provider: &P,
    offsets: &super::PredictedSpkReads,
    consume: impl FnOnce(&PredictedStateAt<'_>, &PredictedRateAt<'_>, &dyn Fn()) -> Result<T>,
) -> Result<T> {
    with_owned_predicted_support_and_positions(
        provider,
        None,
        offsets,
        |state, rates, complete, _| consume(state, rates, complete),
    )
}

/// The optional second view supplies only the position consumed by PPN. All
/// original states still use the full provider; Earth physical velocity is
/// always acquired through that original route, including aberration samples.
fn with_owned_predicted_support_and_positions<P: StateProvider, T>(
    provider: &P,
    deflector_providers: Option<(&P, &P)>,
    offsets: &super::PredictedSpkReads,
    consume: impl FnOnce(
        &PredictedStateAt<'_>,
        &PredictedRateAt<'_>,
        &dyn Fn(),
        &PredictedDeflectorAt<'_>,
    ) -> Result<T>,
) -> Result<T> {
    use std::cell::{Cell, RefCell};
    let reads = RefCell::new(PredictedSupportReads::new());
    let active = Cell::new(None::<PredictedSupportDomains>);
    let cursor = Cell::new(0);
    let acquire = |id, epoch, position_only| {
        offsets.before_read(epoch)?;
        if let Some(domains) = active.get() {
            if !domains.contains(cursor.get(), id, epoch) {
                return offsets.refuse(epoch);
            }
        }
        if active.get().is_some() {
            cursor.set(cursor.get() + 1);
        }
        let source = match (position_only, deflector_providers) {
            (true, Some((original, position))) => {
                if active.get().is_some() {
                    position
                } else {
                    original
                }
            }
            _ => provider,
        };
        let result = icrs_state(source, id, epoch);
        if active.get().is_some() {
            offsets.physical_result(&result);
        } else {
            reads.borrow_mut().note(id, epoch, result.is_ok());
        }
        result
    };
    let state_at = |id, epoch| acquire(id, epoch, false);
    let deflector_at = |id, epoch| acquire(id, epoch, true).map(|state| state.0);
    let rates_at = |place: &dyn Fn(f64) -> Result<[f64; 3]>, epoch, outer_half_width| {
        let Some(domains) = reads.borrow().owned_domains(provider) else {
            return Ok(None);
        };
        active.set(Some(domains));
        let sampled = |at| {
            cursor.set(0);
            let value = place(at);
            offsets.before_read(at)?;
            if cursor.get() != domains.count {
                return offsets.refuse(at);
            }
            value
        };
        let result = offsets.rates_with_direction_hint(
            &sampled,
            epoch,
            outer_half_width,
            true,
            domains.direction_hint(),
        );
        active.set(None);
        result
    };
    let central_complete = || {
        reads.borrow_mut().recording = false;
        offsets.finish_central();
    };
    // The hook is installed by the callers below without reacquiring the
    // central place. Keep source callbacks outside every bookkeeping borrow.
    consume(&state_at, &rates_at, &central_complete, &deflector_at)
}

/// Owned-source analogue of the generic predicted entry. The caller admits
/// supplementary reads from actual LEB/SPK records; no provider capability is
/// inferred. Original center/previous/next acquisitions occur in the shared
/// model reduction before it invokes the optional rate callback.
#[allow(clippy::too_many_arguments)]
fn owned_predicted_model_position(
    state_at: &dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
    deflector_at: &PredictedDeflectorAt<'_>,
    rates_at: &PredictedRateAt<'_>,
    central_complete: &dyn Fn(),
    request_ut1: f64,
    tt: f64,
    body: i32,
    flags: i32,
    site: Option<(f64, f64, f64)>,
    ayanamsha_at: &dyn Fn(f64, i32) -> Result<f64>,
) -> Result<(PositionResult, i32)> {
    for value in [request_ut1, tt] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    let (state, rate_origin) = with_model_rate_origin(rates_at, |rates_at| {
        predicted_center_state_source_time_with_models_rates_and_positions(
            None,
            state_at,
            &|epoch| observed_earth_with_origin(state_at, epoch),
            true,
            tt,
            body,
            flags,
            Some(rates_at),
            Some(central_complete),
            Some(deflector_at),
        )
    })?;
    let state = finish_model_state_with_rate_origin(
        state_at,
        state,
        &|_| Ok(request_ut1),
        tt,
        flags,
        false,
        site,
        ayanamsha_at,
        rate_origin,
    )?;
    Ok((apply_output_flags(state, flags), flags))
}

#[cfg(test)]
mod native_prediction_rate_tests {
    use super::*;
    use std::cell::{Cell, RefCell};

    #[test]
    fn rate_origin_requires_a_consumed_some_and_is_local_to_each_state() {
        let state = [120.0, 2.0, 1.0, 0.25, 0.125, 0.0625];
        for supplied in [Some([0.1, 0.2, 0.3]), None, Some([0.4, 0.5, 0.6])] {
            let calls = Cell::new(0);
            let rates = |place: &dyn Fn(f64) -> Result<[f64; 3]>, tt, width| {
                calls.set(calls.get() + 1);
                assert_eq!(tt, 2451545.0);
                assert_eq!(width, 0.1);
                assert_eq!(place(tt)?, [120.0, 2.0, 1.0]);
                Ok(supplied)
            };
            // Merely supplying an eligible hook does not mean it was used.
            let untouched = with_model_rate_origin(&rates, |_| Ok(state)).unwrap();
            assert_eq!(untouched, (state, ModelRateOrigin::LegacyStencil));
            assert_eq!(calls.get(), 0);
            let (actual, origin) = with_model_rate_origin(&rates, |rates| {
                let mut result = state;
                if let Some(rates) = rates(&|_| Ok([120.0, 2.0, 1.0]), 2451545.0, 0.1)? {
                    result[3..].copy_from_slice(&rates);
                }
                Ok(result)
            })
            .unwrap();
            assert_eq!(actual[..3], state[..3]);
            assert_eq!(
                actual[3..],
                supplied.unwrap_or([state[3], state[4], state[5]])
            );
            assert_eq!(calls.get(), 1);
            assert_eq!(
                origin,
                if supplied.is_some() {
                    ModelRateOrigin::OwnedDerivative
                } else {
                    ModelRateOrigin::LegacyStencil
                }
            );
        }
    }

    #[test]
    fn rate_origin_preserves_first_error_and_rate_callback_order() {
        let calls = RefCell::new(Vec::new());
        let rates = |_: &dyn Fn(f64) -> Result<[f64; 3]>, _, _| {
            calls.borrow_mut().push("rates");
            Err(Error::MathDomain { value: 123.0 })
        };
        let first = with_model_rate_origin(&rates, |_| {
            calls.borrow_mut().push("model");
            Err(Error::MathDomain { value: 456.0 })
        });
        assert!(matches!(first, Err(Error::MathDomain { value: 456.0 })));
        assert_eq!(calls.take(), ["model"]);
        let rate_failure = with_model_rate_origin(&rates, |rates| {
            calls.borrow_mut().push("model");
            rates(
                &|_| panic!("rate callback failed before sampling"),
                2451545.0,
                0.1,
            )?;
            panic!("rate error must stop the calculation")
        });
        assert!(matches!(
            rate_failure,
            Err(Error::MathDomain { value: 123.0 })
        ));
        assert_eq!(calls.take(), ["model", "rates"]);
    }

    #[test]
    fn unrefined_model_tail_retains_its_half_day_nutation_secant() {
        // An unrefined input still carries its declared finite-difference
        // velocity. Changing only its frame derivative is not a qualified
        // derivative of the complete place.
        let tt = 2451545.0;
        let state = [120.0, 2.0, 1.0, 0.25, 0.125, 0.0625];
        let actual = finish_model_state(
            &|_, _| panic!("unused support"),
            state,
            &|_| panic!("unused UT"),
            tt,
            FLG_SPEED,
            true,
            None,
            &|_, _| panic!("unused ayanamsha"),
        )
        .unwrap();
        let secant = (crate::erfa::nut06a(tt + 0.5).unwrap().0
            - crate::erfa::nut06a(tt - 0.5).unwrap().0)
            .to_degrees();
        assert_eq!(actual[3].to_bits(), (state[3] + secant).to_bits());
        assert_eq!(actual[4..], state[4..]);
    }

    #[test]
    fn stationary_mean_ecliptic_ray_has_the_instantaneous_nutation_rotation_rate() {
        // A ray fixed in mean ecliptic axes moves only through longitude
        // nutation. Its Cartesian tangent is z-hat cross the nutated ray.
        // The reference is the independent PyERFA three-step derivative at
        // J2000; the former half-day secant violates the AU/day floor here.
        let independent_rate = 3.5641518015717425e-8;
        let flags = FLG_SPEED | FLG_XYZ;
        let spherical = finish_model_state_with_rate_origin(
            &|_, _| panic!("no support for a supplied geocentric ray"),
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            &|_| panic!("unused UT"),
            2451545.0,
            flags,
            false,
            None,
            &|_, _| panic!("unused ayanamsha"),
            ModelRateOrigin::OwnedDerivative,
        )
        .unwrap();
        let state = apply_output_flags(spherical, flags);
        let expected = [
            -state[1] * independent_rate,
            state[0] * independent_rate,
            0.0,
        ];
        for (value, oracle) in state[3..].iter().zip(expected) {
            assert!((value - oracle).abs() <= 1e-12);
        }
    }

    #[test]
    fn closed_latitude_oscillation_resolves_the_two_point_truncation() {
        // A spherical trajectory has a known instantaneous angular velocity.
        // The coarse model stencil attenuates it by sin(omega*h)/(omega*h).
        let amplitude = 0.01;
        let omega = 1.25;
        let sample = |offset: f64| Ok(Some([120.0, amplitude * (omega * offset).sin()]));
        let received = native_prediction_angular_rates(&sample, (-0.1, 0.1)).unwrap();
        let expected = amplitude * omega;
        assert!(
            (received[1] - expected).abs() < expected * 1e-6,
            "angular derivative: received={}, expected={expected}",
            received[1],
        );
    }

    #[test]
    fn native_jpl_legacy_source_secant_exposes_closed_latitude_truncation() {
        // The NativeJpl entry passes these two scalar-TT callbacks with
        // public_of_date=false and no owned-domain rate callback. Construct
        // an observer whose subtraction leaves a closed spherical motion:
        // latitude=A*sin(omega*(TT-J2000)), so beta'(J2000)=A*omega.
        // This preserves the demonstrated legacy 0.1-day source secant.
        // An opaque callback's domain is not certified and gains no reads.
        let epoch = 2451545.0;
        let amplitude = 0.01;
        let omega = 1.25;
        let body = crate::constants::NEPTUNE_LEVERRIER;
        let observations = RefCell::new(Vec::new());
        let velocities = RefCell::new(Vec::new());
        let (sine, cosine) = (84381.448_f64 / 3600.0).to_radians().sin_cos();
        let observation = |tt: f64| {
            observations.borrow_mut().push(tt.to_bits());
            let model = crate::hypothetical::calc_fictitious_position(body, tt)?;
            let target = model_xyz(model.lon_deg, model.lat_deg, model.dist_au);
            let latitude = amplitude * (omega * (tt - epoch)).sin();
            let desired_ray = model_xyz(120.0, latitude, 40.0);
            let relative: [f64; 3] = std::array::from_fn(|axis| target[axis] - desired_ray[axis]);
            Ok(RetardedEarthObservation {
                relative: [
                    relative[0],
                    cosine * relative[1] - sine * relative[2],
                    sine * relative[1] + cosine * relative[2],
                ],
                sun_reception: [0.0; 3],
            })
        };
        let state = |id, tt: f64| {
            assert_eq!(id, EARTH);
            velocities.borrow_mut().push(tt.to_bits());
            // TRUEPOS retains this acquisition but consumes neither channel.
            Ok(([0.0; 3], [0.0; 3]))
        };
        let actual = predicted_center_state_source_time_with_models(
            None,
            &state,
            &observation,
            false,
            epoch,
            body,
            FLG_J2000 | FLG_SPEED | FLG_TRUEPOS,
        )
        .unwrap();
        let expected_epochs = [epoch, epoch - 0.1, epoch + 0.1].map(f64::to_bits);
        assert_eq!(*observations.borrow(), expected_epochs);
        assert_eq!(*velocities.borrow(), expected_epochs);
        let expected = amplitude * omega;
        assert!(
            (actual[4] - expected).abs() > expected.abs() * 1e-6,
            "legacy source secant unexpectedly resolved the closed derivative",
        );
        let coarse = amplitude
            * ((omega * ((epoch + 0.1) - epoch)).sin() - (omega * ((epoch - 0.1) - epoch)).sin())
            / 0.2;
        assert!((actual[4] - coarse).abs() < expected.abs() * 1e-9);
    }

    #[test]
    fn closed_radial_oscillation_resolves_the_two_point_truncation() {
        // This positive range has an exact instantaneous AU/day derivative.
        // The old two-point radial stencil has the same sinc attenuation as
        // an angular sinusoid, which can cancel an independent angular error
        // when a spherical velocity is projected to Cartesian components.
        let amplitude = 0.01;
        let omega = 1.25;
        let sample = |offset: f64| [120.0, 1.0, 40.0 + amplitude * (omega * offset).sin()];
        let coarse = model_stencil(sample(0.0), sample(-0.1), sample(0.1), 0.1, false)[5];
        let received =
            native_prediction_rates(&|offset| Ok(Some(sample(offset))), (-0.1, 0.1)).unwrap()[2];
        let expected = amplitude * omega;
        assert!((coarse - expected).abs() > expected * 1e-3);
        assert!(
            (received - expected).abs() < expected * 1e-6,
            "radial derivative: received={received}, expected={expected}",
        );
    }

    #[test]
    fn fourth_order_rejects_unowned_offsets_before_acquisition_and_wraps_longitude() {
        let acquired = RefCell::new(Vec::new());
        let rejected = Cell::new(0);
        let sample = |offset: f64| {
            if offset.abs() >= 0.02 {
                rejected.set(rejected.get() + 1);
                return Ok(None);
            }
            acquired.borrow_mut().push(offset);
            Ok(Some([
                (359.96875 + 4.0 * offset).rem_euclid(360.0),
                0.01 * (1.25 * offset).sin(),
                40.0 + 0.01 * (1.25 * offset).sin(),
            ]))
        };
        let rates = native_prediction_rates(&sample, (-0.1, 0.1)).unwrap();
        assert_eq!(rates[0], 4.0);
        assert!((rates[1] - 0.0125).abs() < 0.0125 * 1e-6);
        assert!((rates[2] - 0.0125).abs() < 0.0125 * 1e-6);
        assert!(rejected.get() > 0);
        assert!(acquired.borrow().iter().all(|offset| offset.abs() < 0.02));
    }

    #[derive(Clone, Copy)]
    enum Domain {
        Smooth,
        Joined,
        Separate,
        Unavailable,
        Invalid,
    }

    struct IntrinsicSupport {
        native: bool,
        domain: Domain,
        fail_earth: bool,
        fail_earth_at: Option<f64>,
        fail_optional_after_metadata: bool,
        stationary_earth: bool,
        earth_velocity_scale: f64,
        reads: RefCell<Vec<(i32, u64)>>,
        metadata: Cell<usize>,
    }

    impl IntrinsicSupport {
        fn new(native: bool, domain: Domain) -> Self {
            Self {
                native,
                domain,
                fail_earth: false,
                fail_earth_at: None,
                fail_optional_after_metadata: false,
                stationary_earth: false,
                earth_velocity_scale: 1.0,
                reads: RefCell::new(Vec::new()),
                metadata: Cell::new(0),
            }
        }
    }

    impl StateProvider for IntrinsicSupport {
        fn has_analytic_kinematics(&self) -> bool {
            self.native
        }
        fn provider_name(&self) -> &'static str {
            "closed-intrinsic-model-support"
        }
        fn has_body(&self, body: i32) -> bool {
            PREDICTED_SUPPORT_BODIES.contains(&body)
        }
        fn body_coverage(&self, _: i32) -> Option<(f64, f64)> {
            Some((2451541.0, 2451549.0))
        }
        fn smooth_interval_tt(&self, _: i32, epoch: f64) -> Result<Option<(f64, f64)>> {
            self.metadata.set(self.metadata.get() + 1);
            match self.domain {
                Domain::Smooth => Ok(Some((2451541.0, 2451549.0))),
                Domain::Joined => {
                    let start = (epoch * 16.0).floor() / 16.0;
                    Ok(Some((start, start + 1.0 / 16.0)))
                }
                Domain::Separate => {
                    let start = ((epoch + 0.001953125) * 16.0).floor() / 16.0 - 0.001953125;
                    Ok(Some((start, start + 0.0625)))
                }
                Domain::Unavailable => Ok(None),
                Domain::Invalid => Err(Error::MathDomain { value: epoch }),
            }
        }
        fn serving_coord_type(&self, _: i32) -> Option<u32> {
            Some(COORD_ICRS_BARY)
        }
        fn eval_raw(&self, body: i32, epoch: f64) -> Result<(Vec<f64>, Vec<f64>)> {
            self.reads.borrow_mut().push((body, epoch.to_bits()));
            if body == crate::constants::JUPITER
                && self.fail_optional_after_metadata
                && self.metadata.get() != 0
            {
                return Err(Error::UnknownBody { body_id: body });
            }
            if body == EARTH && (self.fail_earth || self.fail_earth_at == Some(epoch)) {
                return Err(Error::MathDomain { value: epoch });
            }
            let (position, velocity) = match body {
                SUN => ([0.0; 3], [0.0; 3]),
                EARTH if self.stationary_earth => ([0.0; 3], [0.0; 3]),
                EARTH => {
                    let phase = 1.25 * (epoch - 2451545.0);
                    (
                        [1.0, 0.0, 0.001 * phase.sin()],
                        [0.0, 0.0, 0.00125 * phase.cos()],
                    )
                }
                crate::constants::JUPITER => ([5.0, 1.0, 0.0], [0.0; 3]),
                crate::constants::SATURN => ([9.0, -1.0, 0.0], [0.0; 3]),
                _ => return Err(Error::UnknownBody { body_id: body }),
            };
            let velocity = if body == EARTH {
                velocity.map(|value| value * self.earth_velocity_scale)
            } else {
                velocity
            };
            Ok((position.to_vec(), velocity.to_vec()))
        }
        fn eval_raw_two_part(&self, _: i32, _: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)> {
            panic!("predicted scalar support must retain its acquisition clock")
        }
        fn serving_entry_snapshot(&self, _: i32) -> Result<EntrySnapshot> {
            Ok(EntrySnapshot {
                coord_type: COORD_ICRS_BARY,
            })
        }
    }

    fn run(provider: &IntrinsicSupport, flags: i32, models_entry: bool) -> Result<PositionResult> {
        let delta = |_: f64| Ok(0.0);
        let aya = |_: f64, _: i32| Ok(0.0);
        let result = if models_entry {
            serving_model_position_with_models(
                None,
                provider,
                2451545.0,
                2451545.0,
                crate::constants::NEPTUNE_ADAMS,
                flags,
                None,
                &delta,
                &aya,
                0,
            )
        } else {
            serving_model_position(
                provider,
                2451545.0,
                2451545.0,
                crate::constants::NEPTUNE_ADAMS,
                flags,
                None,
                &delta,
                &aya,
                0,
            )
        };
        result.map(|(state, _)| state)
    }

    #[test]
    fn native_refinement_preserves_place_and_original_acquisition_prefix() {
        let flags = FLG_J2000 | FLG_SPEED | FLG_TRUEPOS;
        for models_entry in [false, true] {
            let opaque = IntrinsicSupport::new(false, Domain::Smooth);
            let native = IntrinsicSupport::new(true, Domain::Smooth);
            let before = run(&opaque, flags, models_entry).unwrap();
            let after = run(&native, flags, models_entry).unwrap();
            for component in [0, 1, 2] {
                assert_eq!(after[component].to_bits(), before[component].to_bits());
            }
            assert_ne!(after[4].to_bits(), before[4].to_bits());
            assert_ne!(after[5].to_bits(), before[5].to_bits());
            assert!(native.reads.borrow().starts_with(&opaque.reads.borrow()));
            assert!(native.reads.borrow().len() > opaque.reads.borrow().len());
            assert_eq!(opaque.metadata.get(), 0);
            assert!(native.metadata.get() > 0);
        }
    }

    #[test]
    fn central_leb_roles_refine_when_legacy_sides_consume_other_records() {
        let flags = FLG_J2000 | FLG_SPEED | FLG_TRUEPOS;
        for models_entry in [false, true] {
            let opaque = IntrinsicSupport::new(false, Domain::Separate);
            let before = run(&opaque, flags, models_entry).unwrap();
            let native = IntrinsicSupport::new(true, Domain::Separate);
            let after = run(&native, flags, models_entry).unwrap();
            let control = run(
                &IntrinsicSupport::new(true, Domain::Smooth),
                flags,
                models_entry,
            )
            .unwrap();
            assert_eq!(after[..3], before[..3]);
            assert!(native.reads.borrow().starts_with(&opaque.reads.borrow()));
            assert!(native.reads.borrow().len() > opaque.reads.borrow().len());
            assert_ne!(after[4].to_bits(), before[4].to_bits());
            assert_ne!(after[5].to_bits(), before[5].to_bits());
            for (actual, expected) in after[3..].iter().zip(&control[3..]) {
                assert!((actual - expected).abs() <= expected.abs() * 1e-6);
            }
        }
    }

    #[test]
    fn declined_refinement_keeps_the_legacy_tail_for_generic_and_owned_entries() {
        let flags = FLG_SPEED | FLG_TRUEPOS;
        for domain in [Domain::Joined, Domain::Unavailable, Domain::Invalid] {
            for models_entry in [false, true] {
                let opaque = IntrinsicSupport::new(false, domain);
                let before = run(&opaque, flags, models_entry).unwrap();
                let native = IntrinsicSupport::new(true, domain);
                let after = run(&native, flags, models_entry).unwrap();
                assert_eq!(after.map(f64::to_bits), before.map(f64::to_bits));
                assert_eq!(*native.reads.borrow(), *opaque.reads.borrow());
                assert!(native.metadata.get() > 0);
            }
            let opaque = IntrinsicSupport::new(false, domain);
            let before = run(&opaque, flags, true).unwrap();
            let owned = IntrinsicSupport::new(false, domain);
            let reads = super::super::PredictedSpkReads::new(true);
            let (after, returned_flags) = owned_provider_predicted_position(
                &owned,
                &owned,
                &owned,
                &reads,
                2451545.0,
                2451545.0,
                crate::constants::NEPTUNE_ADAMS,
                flags,
                None,
                &|_, _| panic!("unused ayanamsha"),
            )
            .unwrap();
            assert_eq!(returned_flags, flags);
            assert_eq!(after.map(f64::to_bits), before.map(f64::to_bits));
            assert_eq!(*owned.reads.borrow(), *opaque.reads.borrow());
        }
    }

    #[test]
    fn supplied_models_keep_legacy_tail_and_callbacks_with_intrinsic_support() {
        let mut results = Vec::new();
        for capability in [false, true] {
            let provider = IntrinsicSupport::new(capability, Domain::Smooth);
            let models = RefCell::new(Vec::new());
            let model = |family, tt: f64, source_time| {
                models
                    .borrow_mut()
                    .push((family, tt.to_bits(), source_time));
                Ok(Some(crate::hypothetical::HypoState {
                    lon_deg: 45.0 + 0.01 * (tt - 2451545.0),
                    lat_deg: 3.0,
                    dist_au: 2.0,
                    dlon: 0.01,
                    dlat: 0.0,
                    ddist: 0.0,
                }))
            };
            let output = serving_model_position_with_models(
                Some(&model),
                &provider,
                2451545.0,
                2451545.0,
                VULCAN,
                FLG_SPEED | FLG_TRUEPOS,
                None,
                &|_| panic!("unused Delta-T"),
                &|_, _| panic!("unused ayanamsha"),
                0,
            )
            .unwrap();
            assert!(!models.borrow().is_empty());
            assert_eq!(provider.metadata.get(), 0);
            results.push((
                output.0.map(f64::to_bits),
                output.1,
                models.into_inner(),
                provider.reads.into_inner(),
            ));
        }
        assert_eq!(results[0], results[1]);
    }

    #[test]
    fn joined_or_unsupported_domains_preserve_original_words_and_reads() {
        let flags = FLG_J2000 | FLG_SPEED | FLG_TRUEPOS;
        for domain in [Domain::Joined, Domain::Unavailable, Domain::Invalid] {
            let opaque = IntrinsicSupport::new(false, domain);
            let native = IntrinsicSupport::new(true, domain);
            assert_eq!(
                run(&native, flags, true).unwrap().map(f64::to_bits),
                run(&opaque, flags, true).unwrap().map(f64::to_bits),
            );
            assert_eq!(*native.reads.borrow(), *opaque.reads.borrow());
        }
    }

    #[test]
    fn explicit_opaque_model_preserves_words_callbacks_and_support_on_native_provider() {
        let evaluate = |native, body| {
            let provider = IntrinsicSupport::new(native, Domain::Smooth);
            let model_calls = RefCell::new(Vec::new());
            let model = |family, epoch: f64, source_time| {
                model_calls
                    .borrow_mut()
                    .push((family, epoch.to_bits(), source_time));
                Ok(Some(crate::hypothetical::HypoState {
                    lon_deg: 123.0 + 0.01 * (epoch - 2451545.0),
                    lat_deg: 1.0,
                    dist_au: 2.0,
                    dlon: 0.01,
                    dlat: 0.0,
                    ddist: 0.0,
                }))
            };
            let (state, _) = serving_model_position_with_models(
                Some(&model),
                &provider,
                2451545.0,
                2451545.0,
                body,
                FLG_J2000 | FLG_SPEED | FLG_TRUEPOS,
                None,
                &|_| Ok(0.0),
                &|_, _| Ok(0.0),
                0,
            )
            .unwrap();
            (
                state.map(f64::to_bits),
                provider.reads.into_inner(),
                model_calls.into_inner(),
                provider.metadata.get(),
            )
        };
        for body in [crate::constants::NEPTUNE_ADAMS, VULCAN, PROSERPINA] {
            let opaque = evaluate(false, body);
            let intrinsic = evaluate(true, body);
            assert_eq!(intrinsic, opaque);
            assert_eq!(intrinsic.3, 0);
            if body == VULCAN {
                assert!(!intrinsic.2.is_empty());
            } else if body == PROSERPINA {
                // This built-in formula ignores model callbacks, but an explicit
                // callback still does not authorize supplementary support reads.
                assert!(intrinsic.2.is_empty());
            }
        }
    }

    #[test]
    fn predicted_builtin_family_keeps_capability_and_flag_admission() {
        let native = IntrinsicSupport::new(true, Domain::Smooth);
        let opaque = IntrinsicSupport::new(false, Domain::Smooth);
        for body in HARRINGTON..=PROSERPINA {
            let eligible = (HARRINGTON..=VULCAN).contains(&body) || body == PROSERPINA;
            assert_eq!(
                native_predicted_rates_available(&native, body, FLG_SPEED),
                eligible
            );
            assert!(!native_predicted_rates_available(&opaque, body, FLG_SPEED));
            for flags in [0, FLG_SPEED | FLG_HELCTR, FLG_SPEED | FLG_BARYCTR] {
                assert!(!native_predicted_rates_available(&native, body, flags));
            }
        }
        assert_eq!(native.metadata.get(), 0);
        assert!(native.reads.borrow().is_empty());
    }

    #[test]
    fn nonconsuming_flags_and_original_errors_do_not_acquire_refinement_support() {
        for flags in [
            FLG_J2000,
            FLG_J2000 | FLG_SPEED | FLG_HELCTR,
            FLG_J2000 | FLG_SPEED | FLG_BARYCTR,
        ] {
            let opaque = IntrinsicSupport::new(false, Domain::Smooth);
            let native = IntrinsicSupport::new(true, Domain::Smooth);
            assert_eq!(
                run(&native, flags, true).unwrap().map(f64::to_bits),
                run(&opaque, flags, true).unwrap().map(f64::to_bits),
            );
            assert_eq!(*native.reads.borrow(), *opaque.reads.borrow());
            assert_eq!(native.metadata.get(), 0);
        }
        for native in [false, true] {
            let mut provider = IntrinsicSupport::new(native, Domain::Smooth);
            provider.fail_earth = true;
            assert!(matches!(
                run(&provider, FLG_J2000 | FLG_SPEED, true),
                Err(Error::MathDomain { value: 2451545.0 })
            ));
            assert_eq!(provider.reads.borrow().len(), 2);
            assert_eq!(provider.metadata.get(), 0);
        }
    }

    fn closed_vulcan_rates(epoch: f64) -> [f64; 3] {
        let elements = crate::hypothetical::vulcan_elements();
        let centuries = (epoch - elements.epoch) / 36525.0;
        let mean = (elements.m0 + elements.n_century * centuries)
            .rem_euclid(360.0)
            .to_radians();
        let mut eccentric = mean;
        for _ in 0..8 {
            eccentric -= (eccentric - elements.e * eccentric.sin() - mean)
                / (1.0 - elements.e * eccentric.cos());
        }
        let (sine, cosine) = eccentric.sin_cos();
        let denominator = 1.0 - elements.e * cosine;
        let eccentric_rate = (elements.n_century / 36525.0).to_radians() / denominator;
        let root = (1.0 - elements.e * elements.e).sqrt();
        let anomaly = (root * sine).atan2(cosine - elements.e);
        let argument = anomaly + (elements.omega0 + elements.omega_rate * centuries).to_radians();
        let argument_rate =
            root * eccentric_rate / denominator + (elements.omega_rate / 36525.0).to_radians();
        let (sin_i, cos_i) = elements.i.to_radians().sin_cos();
        let sin_latitude = sin_i * argument.sin();
        let horizontal_squared = 1.0 - sin_latitude * sin_latitude;
        [
            elements.node_rate / 36525.0
                + (cos_i * argument_rate / horizontal_squared).to_degrees(),
            (sin_i * argument.cos() * argument_rate / horizontal_squared.sqrt()).to_degrees(),
            elements.a * elements.e * sine * eccentric_rate,
        ]
    }

    #[test]
    fn native_vulcan_velocity_differentiates_the_declared_moving_ellipse() {
        // With Earth at the Sun and TRUEPOS, the reported mean-of-date orbit
        // is the declared ellipse itself. Kepler's implicit derivative and
        // the prescribed apsidal/node rates give all three exact rates;
        // neither a position stencil nor another engine defines this check.
        let epoch = 2451545.0;
        let expected = closed_vulcan_rates(epoch);
        let mut provider = IntrinsicSupport::new(true, Domain::Smooth);
        provider.stationary_earth = true;
        let (state, _) = serving_model_position(
            &provider,
            epoch,
            epoch,
            VULCAN,
            FLG_SPEED | FLG_TRUEPOS | FLG_NONUT,
            None,
            &|_| Ok(0.0),
            &|_, _| Ok(0.0),
            0,
        )
        .unwrap();
        for component in 0..3 {
            assert!(
                (state[component + 3] - expected[component]).abs()
                    <= expected[component].abs() * 1e-6,
                "component {component}: actual={}, exact={}",
                state[component + 3],
                expected[component],
            );
        }
    }

    #[test]
    fn auto_style_sampled_vulcan_velocity_requires_owned_domain_refinement() {
        // Auto WithCenters has this same false capability. The private owned
        // entry differentiates places, rather than trusting stored velocities.
        // This synthetic provider is not a claim about an actual LEB/SPK file.
        // With Earth at the Sun and TRUEPOS, the reported mean-of-date orbit
        // is the declared ellipse itself. Kepler's implicit derivative and
        // the prescribed apsidal/node rates give all three exact rates;
        // neither a position stencil nor another engine defines this check.
        let epoch = 2451545.0;
        let expected = closed_vulcan_rates(epoch);
        let mut provider = IntrinsicSupport::new(false, Domain::Smooth);
        provider.stationary_earth = true;
        assert!(!provider.has_analytic_kinematics());
        let reads = super::super::PredictedSpkReads::new(true);
        let (state, _) = owned_provider_predicted_position(
            &provider,
            &provider,
            &provider,
            &reads,
            epoch,
            epoch,
            VULCAN,
            FLG_SPEED | FLG_TRUEPOS | FLG_NONUT,
            None,
            &|_, _| Ok(0.0),
        )
        .unwrap();
        for component in 0..3 {
            assert!(
                (state[component + 3] - expected[component]).abs()
                    <= expected[component].abs() * 1e-6,
                "component {component}: actual={}, exact={}",
                state[component + 3],
                expected[component],
            );
        }
    }

    #[test]
    fn owned_auto_entry_retains_opaque_reads_errors_and_metadata_refusals() {
        let flags = FLG_J2000 | FLG_SPEED | FLG_TRUEPOS;
        for domain in [
            Domain::Smooth,
            Domain::Joined,
            Domain::Unavailable,
            Domain::Invalid,
        ] {
            let opaque = IntrinsicSupport::new(false, domain);
            let before = run(&opaque, flags, true).unwrap();
            let owned = IntrinsicSupport::new(false, domain);
            let reads = super::super::PredictedSpkReads::new(true);
            let (after, _) = owned_provider_predicted_position(
                &owned,
                &owned,
                &owned,
                &reads,
                2451545.0,
                2451545.0,
                crate::constants::NEPTUNE_ADAMS,
                flags,
                None,
                &|_, _| Ok(0.0),
            )
            .unwrap();
            assert!(!owned.has_analytic_kinematics());
            assert!(owned.reads.borrow().starts_with(&opaque.reads.borrow()));
            assert_eq!(opaque.metadata.get(), 0);
            if matches!(domain, Domain::Smooth) {
                assert_eq!(after[..3], before[..3]);
                assert_ne!(after[4].to_bits(), before[4].to_bits());
                assert!(owned.reads.borrow().len() > opaque.reads.borrow().len());
            } else {
                assert_eq!(after.map(f64::to_bits), before.map(f64::to_bits));
                assert_eq!(*owned.reads.borrow(), *opaque.reads.borrow());
            }
        }
        let mut owned = IntrinsicSupport::new(false, Domain::Smooth);
        owned.fail_earth = true;
        let reads = super::super::PredictedSpkReads::new(true);
        assert!(matches!(
            owned_provider_predicted_position(
                &owned,
                &owned,
                &owned,
                &reads,
                2451545.0,
                2451545.0,
                crate::constants::NEPTUNE_ADAMS,
                flags,
                None,
                &|_, _| Ok(0.0),
            ),
            Err(Error::MathDomain { value: 2451545.0 })
        ));
        assert_eq!(owned.reads.borrow().len(), 2);
        assert_eq!(owned.metadata.get(), 0);
    }

    #[test]
    fn owned_deflector_view_never_replaces_earth_aberration_velocity() {
        let calculate = |earth_scale, decoy_earth| {
            let mut full = IntrinsicSupport::new(false, Domain::Smooth);
            full.earth_velocity_scale = earth_scale;
            let mut positions = IntrinsicSupport::new(false, Domain::Smooth);
            positions.stationary_earth = decoy_earth;
            let reads = super::super::PredictedSpkReads::new(true);
            let (state, _) = owned_provider_predicted_position(
                &full,
                &positions,
                &positions,
                &reads,
                2451545.0,
                2451545.0,
                crate::constants::NEPTUNE_ADAMS,
                FLG_SPEED | FLG_J2000,
                None,
                &|_, _| panic!("unused ayanamsha"),
            )
            .unwrap();
            let actual = positions.reads.borrow();
            assert!(
                !actual.is_empty(),
                "the position-only view must actually run"
            );
            assert!(actual.iter().all(|&(body, _)| body != EARTH));
            state
        };
        let ordinary = calculate(1.0, false);
        let decoy = calculate(1.0, true);
        assert_eq!(ordinary.map(f64::to_bits), decoy.map(f64::to_bits));
        // Earth's physical velocity is a real input to annual aberration,
        // including its time dependence; it must not use the decoy view.
        let changed = calculate(2.0, true);
        assert_ne!(ordinary[1].to_bits(), changed[1].to_bits());
        assert_ne!(ordinary[4].to_bits(), changed[4].to_bits());
    }

    #[test]
    fn failed_optional_reads_do_not_admit_supplementary_states() {
        let provider = IntrinsicSupport::new(true, Domain::Smooth);
        let mut reads = PredictedSupportReads::new();
        reads.note(SUN, 2451545.0, true);
        reads.note(EARTH, 2451545.0, true);
        reads.note(crate::constants::JUPITER, 2451545.0, false);
        assert!(reads.owned_domains(&provider).is_none());
        assert_eq!(provider.metadata.get(), 0);
    }

    #[test]
    fn forced_leb_keeps_optional_read_error_before_later_domain_refusal() {
        for physical_first in [true, false] {
            let mut provider = IntrinsicSupport::new(true, Domain::Smooth);
            provider.fail_optional_after_metadata = true;
            let errors = super::super::PredictedSpkReads::new(true);
            let epoch = 2451545.0;
            let result = with_owned_predicted_support(
                &provider,
                &errors,
                |state_at, rates_at, central_complete| {
                    for body in [crate::constants::JUPITER, EARTH] {
                        state_at(body, epoch)?;
                    }
                    central_complete();
                    rates_at(
                        &|tt| {
                            let optional = || {
                                crate::positions::frames::deflector_position(
                                    &|body, epoch| state_at(body, epoch).map(|state| state.0),
                                    crate::constants::JUPITER,
                                    tt,
                                )
                            };
                            if physical_first {
                                // The real PPN policy consumes UnknownBody as None.
                                assert_eq!(optional()?, None);
                                let _ = state_at(EARTH, tt + 100.0);
                            } else {
                                let _ = state_at(EARTH, tt + 100.0);
                                let _ = optional();
                            }
                            Ok([120.0, 0.0, 40.0])
                        },
                        epoch,
                        0.1,
                    )
                },
            );
            let acquired = provider.reads.borrow();
            let original: Vec<_> = [crate::constants::JUPITER, EARTH]
                .into_iter()
                .map(|body| (body, epoch.to_bits()))
                .collect();
            assert!(acquired.starts_with(&original));
            if physical_first {
                assert_eq!(
                    result,
                    Err(Error::UnknownBody {
                        body_id: crate::constants::JUPITER
                    })
                );
                assert_eq!(acquired.len(), original.len() + 1);
                assert_eq!(acquired.last().unwrap().0, crate::constants::JUPITER);
            } else {
                assert_eq!(result.unwrap(), None);
                assert_eq!(*acquired, original);
            }
        }
    }

    #[test]
    fn supplementary_physical_error_propagates_after_original_acquisition() {
        let flags = FLG_J2000 | FLG_SPEED | FLG_TRUEPOS;
        let opaque = IntrinsicSupport::new(false, Domain::Smooth);
        run(&opaque, flags, true).unwrap();
        let failure_epoch = 2451545.0 + crate::orbital::NODAPS_SPEED_STEP_DAYS;
        assert!(
            !opaque
                .reads
                .borrow()
                .contains(&(EARTH, failure_epoch.to_bits()))
        );
        let mut native = IntrinsicSupport::new(true, Domain::Smooth);
        native.fail_earth_at = Some(failure_epoch);
        assert!(matches!(
            run(&native, flags, true),
            Err(Error::MathDomain { value }) if value == failure_epoch
        ));
        assert!(native.reads.borrow().starts_with(&opaque.reads.borrow()));
        assert!(native.metadata.get() > 0);
    }
}

/// Precess model angular positions and, only when requested, moving-frame
/// angular rates. Radius/radial rate retain their units and values.
fn precess_model_state(
    mut state: PositionResult,
    tt: f64,
    flags: i32,
    to_j2000: bool,
) -> Result<PositionResult> {
    if flags & FLG_SPEED != 0 {
        return horizons_output::local_precess_state(state, tt, to_j2000);
    }
    let (from, to) = if to_j2000 {
        (tt, crate::hypothetical::J2000_JD)
    } else {
        (crate::hypothetical::J2000_JD, tt)
    };
    (state[0], state[1]) = crate::hypothetical::precess_ecliptic(state[0], state[1], from, to)?;
    state[3..].fill(0.0);
    Ok(state)
}

/// Scalar LEB Earth-from-Sun astrometric acquisition used by the Uranian
/// reader branch. `state_at` reads SSB ICRS AU/AU-day states at scalar TT JD;
/// the central Earth read precedes Sun, so a missing Earth prevents Sun access.
/// Sun is then frozen while up to ten Earth updates subtract the geometric
/// light time `|Earth-Sun| / C_LIGHT_AU_DAY` from that same scalar TT. The
/// 1e-14-day successive-light-time test can stop early; exhaustion returns the
/// last relative ICRS vector in AU. Callback errors propagate at their actual
/// support epoch. Native vector-observe TDB acquisition is a separate route.
fn leb_earth_minus_sun(
    state_at: &dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
    tt: f64,
) -> Result<[f64; 3]> {
    let (earth, _) = state_at(EARTH, tt)?;
    let (sun, _) = state_at(SUN, tt)?;
    let mut relative = std::array::from_fn(|i| earth[i] - sun[i]);
    let mut lt = model_norm(relative) / C_LIGHT_AU_DAY;
    for _ in 0..10 {
        let (earth, _) = state_at(EARTH, tt - lt)?;
        relative = std::array::from_fn(|i| earth[i] - sun[i]);
        let next = model_norm(relative) / C_LIGHT_AU_DAY;
        let converged = (next - lt).abs() < 1e-14;
        lt = next;
        if converged {
            break;
        }
    }
    Ok(relative)
}

/// Retarded Earth and its fixed reception origin, both in ICRS AU.
/// The origin is retained from the same acquisition, never reread for PPN.
pub(super) struct RetardedEarthObservation {
    pub(super) relative: [f64; 3],
    pub(super) sun_reception: [f64; 3],
}

/// Observe Earth from a Sun position fixed at reception through a scalar-TT
/// provider. Convert the TDB emission epoch back to TT for each state lookup.
/// The fixed-point solve uses at most ten iterations and a 1e-12-day
/// light-time update threshold. Positions are ICRS AU; provider errors propagate.
fn observed_earth_minus_sun(
    state_at: &dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
    tt: f64,
) -> Result<[f64; 3]> {
    observed_earth_with_origin(state_at, tt).map(|observation| observation.relative)
}

/// Retain the reception Sun from the same scalar-TT Earth observation.
fn observed_earth_with_origin(
    state_at: &dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
    tt: f64,
) -> Result<RetardedEarthObservation> {
    let (sun, _) = state_at(SUN, tt)?;
    let (mut earth, _) = state_at(EARTH, tt)?;
    let whole = tt.floor();
    let fraction = tt - whole;
    let tdb_fraction =
        fraction + reference_frames::tdb_minus_tt_seconds(whole, fraction)? / 86400.0;
    let mut previous_lt = 0.0;
    for _ in 0..10 {
        let relative = std::array::from_fn(|i| earth[i] - sun[i]);
        let lt = model_norm(relative) / reference_frames::FRAME_LIGHT_SPEED_AU_DAY_LOCAL;
        if (lt - previous_lt).abs() < 1e-12 {
            break;
        }
        let retarded_tdb = tdb_fraction - lt;
        let stored_tt =
            retarded_tdb - reference_frames::tdb_minus_tt_seconds(whole, retarded_tdb)? / 86400.0;
        earth = state_at(EARTH, whole + stored_tt)?.0;
        previous_lt = lt;
    }
    Ok(RetardedEarthObservation {
        relative: std::array::from_fn(|i| earth[i] - sun[i]),
        sun_reception: sun,
    })
}

/// Normalized first-order Bradley aberration on a direction. Both 1e-30
/// magnitude guards and the unit-direction re-normalization are retained.
pub(in crate::positions) fn bradley_aberration(position: [f64; 3], velocity: [f64; 3]) -> [f64; 3] {
    let distance = model_norm(position);
    if distance < 1e-30 {
        return position;
    }
    let unit = position.map(|x| x / distance);
    let unit_norm = model_norm(unit);
    if unit_norm < 1e-30 {
        return position;
    }
    let direction = unit.map(|x| x / unit_norm);
    let v = velocity.map(|x| x / C_LIGHT_AU_DAY);
    let dot = direction[0] * v[0] + direction[1] * v[1] + direction[2] * v[2];
    let apparent = std::array::from_fn(|i| unit[i] + (v[i] - dot * direction[i]));
    let magnitude = model_norm(apparent);
    if magnitude < 1e-30 {
        return position;
    }
    apparent.map(|x| x / magnitude * distance)
}

/// Uranian center composition from geometric J2000 model coordinates.
///
/// `tt` and each centered one-day model sample are scalar TT JD. `state_at`
/// supplies SSB ICRS Sun/Earth position and velocity in AU/AU-day; barycentric
/// composition uses Sun, and apparent geocentric samples use Earth's velocity.
/// `earth_helio_at` independently returns an astrometric ICRS Earth-from-Sun
/// position in AU at each sample: sealed callers retain their direct-reader
/// scalar TT iteration, while native callers retain observation-TDB retardation.
/// Heliocentric requests acquire neither support callback. Geocentric apparent
/// samples apply normalized first-order Bradley aberration unless TRUEPOS or
/// NOABERR suppresses it. SPEED selects centered angular/distance differences
/// with longitude unwrapping on geocentric output; barycentric composition
/// always samples its Sun/model stencil before its caller zeros absent-SPEED
/// rates. Returned rates without SPEED are positive zeros.
///
/// Caller-normalized BARYCTR precedes HELCTR. The returned spherical state is
/// degrees/AU and per-day rates in the requested J2000 or mean-of-date ecliptic
/// frame, before nutation, site, sidereal and equatorial conversion. Invalid
/// model domains, missing support, non-finite states, ranges and frame errors
/// propagate. These conventional analytical models are not physical bodies.
fn uranian_center_state(
    state_at: &dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
    earth_helio_at: &dyn Fn(f64) -> Result<[f64; 3]>,
    tt: f64,
    body: i32,
    flags: i32,
    source_time: bool,
) -> Result<PositionResult> {
    if flags & FLG_BARYCTR != 0 {
        let mut state = uranian_bary_time_state(state_at, tt, body, source_time)?;
        if flags & FLG_SPEED == 0 {
            state[3..].fill(0.0);
        }
        if flags & FLG_J2000 == 0 {
            state = precess_model_state(state, tt, flags, false)?;
        }
        return Ok(state);
    }
    let sample = |jd| -> Result<[f64; 3]> {
        let h = if source_time {
            crate::hypothetical::calc_uranian_planet_source_time(body, jd)?
        } else {
            crate::hypothetical::calc_uranian_planet(body, jd)?
        };
        let helio = model_xyz(h.lon_deg, h.lat_deg, h.dist_au);
        if flags & FLG_HELCTR != 0 {
            return Ok([h.lon_deg, h.lat_deg, h.dist_au]);
        }
        let earth = rotate_icrs_to_iau1976_ecliptic_j2000(earth_helio_at(jd)?);
        let mut geo = std::array::from_fn(|i| helio[i] - earth[i]);
        if flags & (FLG_TRUEPOS | FLG_NOABERR) == 0 {
            let velocity = rotate_icrs_to_iau1976_ecliptic_j2000(state_at(EARTH, jd)?.1);
            geo = bradley_aberration(geo, velocity);
        }
        Ok(model_spherical(geo))
    };
    let center = sample(tt)?;
    // Centered branches always compute their model state; geocentric speed
    // has an explicit acquisition gate in the public contract.
    let mut state = if flags & FLG_SPEED != 0 || flags & (FLG_HELCTR | FLG_BARYCTR) != 0 {
        model_stencil(center, sample(tt - 1.0)?, sample(tt + 1.0)?, 1.0, false)
    } else {
        [center[0], center[1], center[2], 0.0, 0.0, 0.0]
    };
    if flags & FLG_SPEED == 0 {
        state[3..].fill(0.0);
    }
    if flags & FLG_J2000 == 0 {
        state = precess_model_state(state, tt, flags, false)?;
    }
    Ok(state)
}

/// Shared Uranian barycentric stencil with explicit model Time admission.
///
/// At scalar TT JD `tt` and its two one-day neighbors, first acquire the
/// heliocentric model, then `state_at(SUN, jd)` in ICRS AU/AU-day. Rotate its
/// position to J2000 ecliptic and add it to the model's AU vector. Returns
/// barycentric J2000-ecliptic longitude/latitude in degrees, distance AU and
/// centered rates in degrees/day or AU/day; longitude wraps through
/// [`model_stencil`]. `source_time` selects internal IEEE model admission,
/// while checked callers keep their direct model guard. Unknown model,
/// consumed support and defining-frame errors propagate. This helper neither
/// shapes public flags nor certifies physical accuracy of the model.
fn uranian_bary_time_state(
    state_at: &dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
    tt: f64,
    body: i32,
    source_time: bool,
) -> Result<PositionResult> {
    let sample = |jd| -> Result<[f64; 3]> {
        let h = if source_time {
            crate::hypothetical::calc_uranian_planet_source_time(body, jd)?
        } else {
            crate::hypothetical::calc_uranian_planet(body, jd)?
        };
        let helio = model_xyz(h.lon_deg, h.lat_deg, h.dist_au);
        let sun = rotate_icrs_to_iau1976_ecliptic_j2000(state_at(SUN, jd)?.0);
        Ok(model_spherical(std::array::from_fn(|i| helio[i] + sun[i])))
    };
    Ok(model_stencil(
        sample(tt)?,
        sample(tt - 1.0)?,
        sample(tt + 1.0)?,
        1.0,
        false,
    ))
}

/// Finish an acquired Uranian barycentric J2000 spherical state.
///
/// `state` contains degrees/AU and per-day rates from the caller's unchanged
/// Sun/model one-day stencil. `tt` is scalar TT JD for frame/nutation reduction;
/// caller-normalized flags select BARYCTR without TOPOCTR. Rates are zeroed
/// without SPEED before the same requested J2000-to-date precession. The common
/// model tail then adds nutation, ecliptic sidereal correction and equatorial
/// conversion, retaining their order and arithmetic. No support is reacquired.
/// `ut_at` lazily supplies the independent UT1 JD only at the consuming
/// ecliptic sidereal stage; `aya` supplies flag-aware ayanamsha in degrees.
/// Frame, non-finite and consumed clock/ayanamsha errors propagate. The public
/// position entry separately owns XYZ/RADIANS shaping and returned flags.
pub(super) fn finish_uranian_bary_state(
    mut state: PositionResult,
    ut_at: &dyn Fn(f64) -> Result<f64>,
    tt: f64,
    flags: i32,
    aya: &dyn Fn(f64, i32) -> Result<f64>,
) -> Result<PositionResult> {
    if flags & FLG_SPEED == 0 {
        state[3..].fill(0.0);
    }
    if flags & FLG_J2000 == 0 {
        state = precess_model_state(state, tt, flags, false)?;
    }
    let no_support = |body, _: f64| Err(Error::SpkUnregisteredBody { body });
    finish_model_state(&no_support, state, ut_at, tt, flags, true, None, aya)
}

/// Isis center composition in J2000 or mean-of-date ecliptic coordinates.
///
/// `tt` and speed samples are scalar TT JD. `state_at` supplies SSB ICRS
/// Sun/Earth AU/AU-day states; `earth_helio_at` separately supplies the
/// astrometric ICRS Earth-from-Sun AU vector with the caller's observation
/// clock and provider channel. Heliocentric output consumes no support.
/// Non-TRUEPOS target acquisition performs exactly three geometric light-time
/// updates using `radius / C_LIGHT_AU_DAY`, with the observer fixed. Geocentric
/// apparent output adds first-order observer velocity unless TRUEPOS/NOABERR;
/// barycentric output adds the observed-epoch Sun after target retardation.
///
/// SPEED uses a centered one-day geocentric stencil, a centered 0.01-day
/// barycentric stencil, or the heliocentric model's reported rates. Without
/// SPEED, rates are positive zeros. Outputs are degrees/AU and per-day rates
/// before shared nutation/surface observer/sidereal/equatorial reduction. Caller-normalized
/// center flags and finite model times are required; model, callback, coverage,
/// registration, non-finite and frame errors propagate without substitution.
fn isis_center_state(
    state_at: &dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
    earth_helio_at: &dyn Fn(f64) -> Result<[f64; 3]>,
    tt: f64,
    flags: i32,
) -> Result<PositionResult> {
    let helio_lt = |jd| -> Result<f64> {
        if flags & FLG_TRUEPOS != 0 {
            return Ok(0.0);
        }
        let mut lt = 0.0;
        for _ in 0..3 {
            lt = crate::hypothetical::calc_transpluto(jd - lt)?.dist_au / C_LIGHT_AU_DAY;
        }
        Ok(lt)
    };
    if flags & FLG_HELCTR != 0 {
        let h = crate::hypothetical::calc_transpluto(tt - helio_lt(tt)?)?;
        let mut state = [h.lon_deg, h.lat_deg, h.dist_au, h.dlon, h.dlat, h.ddist];
        if flags & FLG_SPEED == 0 {
            state[3..].fill(0.0);
        }
        if flags & FLG_J2000 == 0 {
            state = precess_model_state(state, tt, flags, false)?;
        }
        return Ok(state);
    }
    let bary = flags & FLG_BARYCTR != 0;
    let sample = |jd| -> Result<[f64; 3]> {
        if bary {
            let h = crate::hypothetical::calc_transpluto(jd - helio_lt(jd)?)?;
            let helio = model_xyz(h.lon_deg, h.lat_deg, h.dist_au);
            let sun = rotate_icrs_to_iau1976_ecliptic_j2000(state_at(SUN, jd)?.0);
            let mut result = model_spherical(std::array::from_fn(|i| helio[i] + sun[i]));
            if flags & FLG_J2000 == 0 {
                (result[0], result[1]) = crate::hypothetical::precess_ecliptic(
                    result[0],
                    result[1],
                    crate::hypothetical::J2000_JD,
                    jd,
                )?;
            }
            return Ok(result);
        }
        let earth = rotate_icrs_to_iau1976_ecliptic_j2000(earth_helio_at(jd)?);
        let mut lt = 0.0;
        let mut geo = [0.0; 3];
        for _ in 0..3 {
            let h = crate::hypothetical::calc_transpluto(jd - lt)?;
            let target = model_xyz(h.lon_deg, h.lat_deg, h.dist_au);
            geo = std::array::from_fn(|i| target[i] - earth[i]);
            if flags & FLG_TRUEPOS != 0 {
                break;
            }
            lt = model_norm(geo) / C_LIGHT_AU_DAY;
        }
        let radius = model_norm(geo);
        if flags & (FLG_TRUEPOS | FLG_NOABERR) == 0 && radius > 0.0 {
            let velocity = rotate_icrs_to_iau1976_ecliptic_j2000(state_at(EARTH, jd)?.1);
            geo = std::array::from_fn(|i| geo[i] + velocity[i] * radius / C_LIGHT_AU_DAY);
        }
        Ok(model_spherical(geo))
    };
    let center = sample(tt)?;
    let step = if bary { 0.01 } else { 1.0 };
    let mut state = if flags & FLG_SPEED != 0 {
        model_stencil(center, sample(tt - step)?, sample(tt + step)?, step, bary)
    } else {
        [center[0], center[1], center[2], 0.0, 0.0, 0.0]
    };
    if !bary && flags & FLG_J2000 == 0 {
        state = precess_model_state(state, tt, flags, false)?;
    }
    Ok(state)
}

/// Recenter a symbolic mean-of-date state using the public model branch's
/// J2000 support vectors; rates retain the symbolic model's own convention.
pub(super) fn symbolic_center_state(
    state_at: &dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
    earth_helio_at: &dyn Fn(f64) -> Result<[f64; 3]>,
    tt: f64,
    flags: i32,
    symbolic: crate::hypothetical::HypoState,
) -> Result<PositionResult> {
    // A geometric/NOABERR geocentric request keeps the symbolic spherical
    // state directly. Reconstructing its Cartesian vector would introduce
    // rounding into the fixed radius and change near-apsis element branches.
    if flags & (FLG_HELCTR | FLG_BARYCTR) == 0 && flags & (FLG_TRUEPOS | FLG_NOABERR) != 0 {
        return Ok([
            symbolic.lon_deg,
            symbolic.lat_deg,
            symbolic.dist_au,
            symbolic.dlon,
            symbolic.dlat,
            symbolic.ddist,
        ]);
    }
    let mut xyz = model_xyz(symbolic.lon_deg, symbolic.lat_deg, symbolic.dist_au);
    if flags & (FLG_HELCTR | FLG_BARYCTR) != 0 {
        let origin = if flags & FLG_HELCTR != 0 {
            earth_helio_at(tt)?
        } else {
            state_at(EARTH, tt)?.0
        };
        let origin = rotate_icrs_to_iau1976_ecliptic_j2000(origin);
        xyz = std::array::from_fn(|i| xyz[i] + origin[i]);
    } else if flags & (FLG_TRUEPOS | FLG_NOABERR) == 0 {
        let velocity = rotate_icrs_to_iau1976_ecliptic_j2000(state_at(EARTH, tt)?.1);
        let radius = model_norm(xyz);
        xyz = std::array::from_fn(|i| xyz[i] + velocity[i] * radius / C_LIGHT_AU_DAY);
    }
    let s = model_spherical(xyz);
    Ok([
        s[0],
        s[1],
        s[2],
        symbolic.dlon,
        symbolic.dlat,
        symbolic.ddist,
    ])
}

/// Apply the shared public model tail at its actual intermediate frame.
/// Complete a spherical model state in degrees/AU and per-day rates.
///
/// `tt` is scalar TT JD for nutation, model parallax and output frames.
/// Uranian/Isis states have already selected J2000 or mean-of-date ecliptic
/// axes; predictions/symbolic states remain of-date until final conversion.
/// The tail first zeros rates without SPEED, adds of-date nutation where
/// applicable (with its analytic TT nutation rate), then applies surface observer,
/// ecliptic sidereal correction and equatorial conversion in that order.
///
/// `state_at` supplies SSB ICRS Earth AU/AU-day states at scalar TT for the
/// separate observer and geocenter reads. `site` is east-positive longitude
/// and latitude in degrees with height in metres; source surface observer keeps its own
/// Earth-rotation clock. `ut_at` lazily supplies independent UT1 JD only when
/// SIDEREAL is set without EQUATORIAL. `aya` then returns flag-aware ayanamsha
/// in degrees at that UT1. Unused clocks and support are not acquired.
/// Missing observer, callback/support, non-finite, range and frame failures
/// propagate at the first consumed operation. XYZ/RADIANS shaping is external.
#[allow(clippy::too_many_arguments)]
fn finish_model_state(
    state_at: &dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
    state: PositionResult,
    ut_at: &dyn Fn(f64) -> Result<f64>,
    tt: f64,
    flags: i32,
    already_requested_frame: bool,
    site: Option<(f64, f64, f64)>,
    aya: &dyn Fn(f64, i32) -> Result<f64>,
) -> Result<PositionResult> {
    finish_model_state_with_rate_origin(
        state_at,
        state,
        ut_at,
        tt,
        flags,
        already_requested_frame,
        site,
        aya,
        ModelRateOrigin::LegacyStencil,
    )
}

/// Apply an instantaneous nutation rate only to a genuinely refined place.
/// This does not qualify an unrefined model or change its supplied velocity.
#[allow(clippy::too_many_arguments)]
fn finish_model_state_with_rate_origin(
    state_at: &dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
    mut state: PositionResult,
    ut_at: &dyn Fn(f64) -> Result<f64>,
    tt: f64,
    flags: i32,
    already_requested_frame: bool,
    site: Option<(f64, f64, f64)>,
    aya: &dyn Fn(f64, i32) -> Result<f64>,
    rate_origin: ModelRateOrigin,
) -> Result<PositionResult> {
    if flags & FLG_SPEED == 0 {
        state[3..].fill(0.0);
    }
    let nut_free = flags & (FLG_NONUT | FLG_J2000) != 0
        || flags & (FLG_SIDEREAL | FLG_EQUATORIAL) == (FLG_SIDEREAL | FLG_EQUATORIAL);
    if !nut_free {
        state[0] = (state[0] + crate::erfa::nut06a(tt)?.0.to_degrees()).rem_euclid(360.0);
        if flags & FLG_SPEED != 0 {
            state[3] += match rate_origin {
                ModelRateOrigin::OwnedDerivative => {
                    crate::erfa::nut06a_longitude_rate(tt)?.to_degrees()
                }
                ModelRateOrigin::LegacyStencil => (crate::erfa::nut06a(tt + 0.5)?.0
                    - crate::erfa::nut06a(tt - 0.5)?.0)
                    .to_degrees(),
            };
        }
    }
    if let Some(site) = hypothetical_topo_site(flags, site)? {
        let offset = source_model_topo_offset(state_at, tt, site)?;
        state = model_parallax(
            state,
            tt,
            already_requested_frame && flags & FLG_J2000 != 0,
            offset,
        )?;
    }
    if flags & FLG_SIDEREAL != 0 && flags & FLG_EQUATORIAL == 0 {
        let ut = ut_at(tt)?;
        if !ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: ut });
        }
        state = horizons_output::local_sidereal(state, ut, flags, aya)?;
    }
    horizons_output::local_convert(state, tt, flags, already_requested_frame)
}

/// Public surface observer subtraction in AU/AU-day retains the observer's Earth read,
/// site addition, then the independent geocenter read and subtraction before
/// rotating through the model helper's fixed 23.4392911-degree axes. Scalar
/// TT is retained by the support callback; the site has degrees/metre units.
/// Missing Earth support precedes surface observer geometry; either read can fail.
fn source_model_topo_offset(
    state_at: &dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
    tt: f64,
    site: (f64, f64, f64),
) -> Result<([f64; 3], [f64; 3])> {
    let (earth, earth_v) = state_at(EARTH, tt)?;
    let day = tt.floor();
    let (site_p, site_v) = local_jpl_topo(day, tt - day, site)?;
    let observer: [f64; 3] = std::array::from_fn(|i| earth[i] + site_p[i]);
    let observer_v: [f64; 3] = std::array::from_fn(|i| earth_v[i] + site_v[i]);
    let (geocenter, geocenter_v) = state_at(EARTH, tt)?;
    let position = std::array::from_fn(|i| observer[i] - geocenter[i]);
    let velocity = std::array::from_fn(|i| observer_v[i] - geocenter_v[i]);
    let eps = 23.4392911_f64.to_radians();
    let rotate = |v: [f64; 3]| {
        [
            v[0],
            v[1] * eps.cos() + v[2] * eps.sin(),
            -v[1] * eps.sin() + v[2] * eps.cos(),
        ]
    };
    Ok((rotate(position), rotate(velocity)))
}

/// Subtract the observer state in the result's ecliptic frame. The J2000
/// branch skips the direction precession; of-date parallax uses the same
/// magnitude-preserving precession as the public hypothetical helper.
pub(in crate::positions) fn model_parallax(
    state: PositionResult,
    tt: f64,
    j2000: bool,
    mut offset: ([f64; 3], [f64; 3]),
) -> Result<PositionResult> {
    if !j2000 {
        let rotate = |v: [f64; 3]| -> Result<[f64; 3]> {
            let radius = model_norm(v);
            if radius == 0.0 {
                return Ok([0.0; 3]);
            }
            let lon = v[1].atan2(v[0]).to_degrees();
            let lat = (v[2] / radius).clamp(-1.0, 1.0).asin().to_degrees();
            let (lon, lat) =
                crate::hypothetical::precess_ecliptic(lon, lat, crate::hypothetical::J2000_JD, tt)?;
            Ok(model_xyz(lon, lat, radius))
        };
        offset = (rotate(offset.0)?, rotate(offset.1)?);
    }
    let [lon, lat, dist, dlon, dlat, ddist] = state;
    let (lr, br, dlr, dbr) = (
        lon.to_radians(),
        lat.to_radians(),
        dlon.to_radians(),
        dlat.to_radians(),
    );
    let (cl, sl, cb, sb) = (lr.cos(), lr.sin(), br.cos(), br.sin());
    let body = [dist * cb * cl, dist * cb * sl, dist * sb];
    let body_v = [
        ddist * cb * cl - dist * sb * dbr * cl - dist * cb * sl * dlr,
        ddist * cb * sl - dist * sb * dbr * sl + dist * cb * cl * dlr,
        ddist * sb + dist * cb * dbr,
    ];
    let p: [f64; 3] = std::array::from_fn(|i| body[i] - offset.0[i]);
    let v: [f64; 3] = std::array::from_fn(|i| body_v[i] - offset.1[i]);
    let rxy2 = p[0] * p[0] + p[1] * p[1];
    let r = (rxy2 + p[2] * p[2]).sqrt();
    if r == 0.0 {
        return Ok(state);
    }
    let rxy = rxy2.sqrt();
    let dr = (p[0] * v[0] + p[1] * v[1] + p[2] * v[2]) / r;
    Ok([
        p[1].atan2(p[0]).to_degrees().rem_euclid(360.0),
        (p[2] / r).clamp(-1.0, 1.0).asin().to_degrees(),
        r,
        if rxy2 > 0.0 {
            ((p[0] * v[1] - p[1] * v[0]) / rxy2).to_degrees()
        } else {
            0.0
        },
        if rxy > 0.0 {
            ((v[2] - (p[2] / r) * dr) / rxy).to_degrees()
        } else {
            0.0
        },
        dr,
    ])
}

/// Predicted-model acquisition through independent ICRS and Earth-observe callbacks.
///
/// `public_of_date` selects the public model branch's mean-of-date intermediate
/// frame; the native model-frame caller retains its established tail convention.
/// Scalar epochs are TT JD; states are degrees/AU with per-day rates. The
/// Earth observation retains both the retarded ICRS ray and its already acquired
/// reception Sun; PPN deflectors use that fixed origin and the model's IAU 1976
/// ecliptic axes. Acquisition, model-domain and precession failures propagate.
pub(super) fn predicted_center_state(
    state_at: &dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
    earth_observation_at: &dyn Fn(f64) -> Result<RetardedEarthObservation>,
    public_of_date: bool,
    tt: f64,
    body: i32,
    flags: i32,
) -> Result<PositionResult> {
    if !tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: tt });
    }
    predicted_center_state_source_time(
        state_at,
        earth_observation_at,
        public_of_date,
        tt,
        body,
        flags,
    )
}

/// Predicted center reduction at the scalar TT of a guarded public Time.
/// The geocentric observation is acquired before model lookup/propagation,
/// including when internal clock arithmetic is nonfinite. Direct model
/// domains and the ordinary checked center entry retain their refusals.
/// Heliocentric/barycentric branches continue their model-first acquisition.
/// Callback/frame/units, flags, stencils and consumed errors otherwise follow
/// [`predicted_center_state`]; this variant changes acquired-Time dispatch.
pub(super) fn predicted_center_state_source_time(
    state_at: &dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
    earth_observation_at: &dyn Fn(f64) -> Result<RetardedEarthObservation>,
    public_of_date: bool,
    tt: f64,
    body: i32,
    flags: i32,
) -> Result<PositionResult> {
    let of_date = body == VULCAN || body == PROSERPINA;
    let model_at = |jd: f64| -> Result<crate::hypothetical::HypoState> {
        if body == VULCAN {
            crate::hypothetical::calc_vulcan(jd)
        } else if body == PROSERPINA {
            crate::hypothetical::calc_proserpina(jd)
        } else {
            crate::hypothetical::calc_fictitious_position(body, jd)
        }
    };
    // Heliocentric J2000-ecliptic Cartesian model place at `jd`: Vulcan and
    // Proserpina use the equinox of date, so their angles rotate to J2000
    // before the Earth subtraction is frame-consistent.
    let helio_xyz_j2000 = |jd: f64| -> Result<[f64; 3]> {
        let state = model_at(jd)?;
        let (lon, lat) = if of_date {
            crate::hypothetical::precess_ecliptic(
                state.lon_deg,
                state.lat_deg,
                jd,
                crate::hypothetical::J2000_JD,
            )?
        } else {
            (state.lon_deg, state.lat_deg)
        };
        let (lon_r, lat_r) = (lon.to_radians(), lat.to_radians());
        let cl = lat_r.cos();
        Ok([
            state.dist_au * cl * lon_r.cos(),
            state.dist_au * cl * lon_r.sin(),
            state.dist_au * lat_r.sin(),
        ])
    };
    let cartesian_norm =
        |xyz: [f64; 3]| -> f64 { (xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]).sqrt() };
    let to_spherical = |xyz: [f64; 3]| -> [f64; 3] {
        let r = cartesian_norm(xyz);
        [
            xyz[1].atan2(xyz[0]).to_degrees().rem_euclid(360.0),
            if r > 0.0 {
                (xyz[2] / r).clamp(-1.0, 1.0).asin().to_degrees()
            } else {
                0.0
            },
            r,
        ]
    };
    // Centered longitude rate with the angular difference reduced to a half turn.
    let wrap_rate = |next_lon: f64, prev_lon: f64, h: f64| -> f64 {
        ((next_lon - prev_lon + 180.0).rem_euclid(360.0) - 180.0) / (2.0 * h)
    };
    let truepos = flags & FLG_TRUEPOS != 0;
    // Sun-to-body light time in days at `jd` (0 under TRUEPOS), iterated
    // like the geocentric case: the heliocentric distance is itself
    // evaluated at the retarded epoch.
    let helio_light_time = |jd: f64| -> Result<f64> {
        if truepos {
            return Ok(0.0);
        }
        let mut lt = 0.0;
        for _ in 0..3 {
            lt = cartesian_norm(helio_xyz_j2000(jd - lt)?) / C_LIGHT_AU_DAY;
        }
        Ok(lt)
    };
    // Rotate a J2000-ecliptic place to the mean ecliptic at TT Julian Day `jd`.
    let to_of_date = |lon: f64, lat: f64, jd: f64| -> Result<(f64, f64)> {
        crate::hypothetical::precess_ecliptic(lon, lat, crate::hypothetical::J2000_JD, jd)
    };
    if flags & FLG_BARYCTR != 0 && flags & FLG_HELCTR == 0 {
        // Barycentric place: the model is heliocentric, so the Sun's
        // barycentric vector is added in J2000-ecliptic axes.
        let place_j2000 = |jd: f64| -> Result<[f64; 3]> {
            let lt = helio_light_time(jd)?;
            let helio = helio_xyz_j2000(jd - lt)?;
            let (sun, _) = state_at(SUN, jd)?;
            let sun = rotate_icrs_to_iau1976_ecliptic_j2000(sun);
            Ok(to_spherical([
                helio[0] + sun[0],
                helio[1] + sun[1],
                helio[2] + sun[2],
            ]))
        };
        let center = place_j2000(tt)?;
        let (lon, lat) = if of_date || public_of_date {
            to_of_date(center[0], center[1], tt)?
        } else {
            (center[0], center[1])
        };
        if public_of_date && flags & FLG_SPEED == 0 {
            return Ok([lon, lat, center[2], 0.0, 0.0, 0.0]);
        }
        let h = 0.01;
        let prev = place_j2000(tt - h)?;
        let next = place_j2000(tt + h)?;
        let (p_lon, p_lat) = if of_date || public_of_date {
            to_of_date(prev[0], prev[1], tt - h)?
        } else {
            (prev[0], prev[1])
        };
        let (n_lon, n_lat) = if of_date || public_of_date {
            to_of_date(next[0], next[1], tt + h)?
        } else {
            (next[0], next[1])
        };
        return Ok([
            lon,
            lat,
            center[2],
            wrap_rate(n_lon, p_lon, h),
            (n_lat - p_lat) / (2.0 * h),
            (next[2] - prev[2]) / (2.0 * h),
        ]);
    }
    if flags & FLG_HELCTR != 0 {
        // Apparent heliocentric place: the model at the Sun-to-body
        // retarded epoch. TRUEPOS zeroes that light time. Speeds
        // central-difference the reported quantity so they stay the
        // derivative of the reported position.
        let h = 0.01;
        let reported = |jd: f64| -> Result<[f64; 3]> {
            let state = model_at(jd - helio_light_time(jd)?)?;
            // Vulcan/Proserpina models already use the equinox of date; the
            // predicted-planet propagation already returns J2000. Both match
            // the tail convention without further rotation here.
            Ok([state.lon_deg, state.lat_deg, state.dist_au])
        };
        let center = reported(tt)?;
        // Vulcan/Proserpina models already use the equinox of date; the
        // predicted-planet propagation already returns J2000. Both match
        // the tail convention without further rotation here.
        if public_of_date && flags & FLG_SPEED == 0 {
            let mut state = [center[0], center[1], center[2], 0.0, 0.0, 0.0];
            if !of_date {
                state = precess_model_state(state, tt, flags, false)?;
            }
            return Ok(state);
        }
        let prev = reported(tt - h)?;
        let next = reported(tt + h)?;
        let mut state = [
            center[0],
            center[1],
            center[2],
            wrap_rate(next[0], prev[0], h),
            (next[1] - prev[1]) / (2.0 * h),
            (next[2] - prev[2]) / (2.0 * h),
        ];
        if public_of_date && !of_date {
            state = precess_model_state(state, tt, flags, false)?;
        }
        return Ok(state);
    }
    // Construct a geocentric apparent place in the mean ecliptic of date.
    // Retard target and Earth, then apply PPN deflection and annual aberration
    // unless geometric flags suppress them. Reevaluate the observer at samples.
    let noaberr = flags & (FLG_TRUEPOS | FLG_NOABERR) != 0;
    let nogdefl = flags & (FLG_TRUEPOS | FLG_NOGDEFL) != 0;
    let geo_of_date = |jd: f64| -> Result<[f64; 3]> {
        let observation = earth_observation_at(jd)?;
        let earth_j2000 = rotate_icrs_to_iau1976_ecliptic_j2000(observation.relative);
        let (_, earth_vel_icrs) = state_at(EARTH, jd)?;
        let earth_vel = rotate_icrs_to_iau1976_ecliptic_j2000(earth_vel_icrs);
        let mut lt = 0.0;
        let mut geo = [0.0; 3];
        for _ in 0..3 {
            let target = helio_xyz_j2000(jd - lt)?;
            geo = [
                target[0] - earth_j2000[0],
                target[1] - earth_j2000[1],
                target[2] - earth_j2000[2],
            ];
            if truepos {
                break;
            }
            lt = cartesian_norm(geo) / C_LIGHT_AU_DAY;
        }
        let mut r = cartesian_norm(geo);
        if !nogdefl && r > 0.0 {
            // All PPN vectors use the model's bias-free IAU 1976 ecliptic
            // axes and the same Sun origin fixed at this reception epoch.
            geo = apply_gravitational_deflection(
                geo,
                earth_j2000,
                jd,
                if public_of_date {
                    r / C_LIGHT_AU_DAY
                } else {
                    lt
                },
                |id, epoch| {
                    let position = state_at(id, epoch)?.0;
                    deflector_in_frame(position, observation.sun_reception, |relative| {
                        Ok(rotate_icrs_to_iau1976_ecliptic_j2000(relative))
                    })
                },
            )?;
            r = cartesian_norm(geo);
        }
        if !noaberr && r > 0.0 {
            geo = [
                geo[0] + earth_vel[0] * r / C_LIGHT_AU_DAY,
                geo[1] + earth_vel[1] * r / C_LIGHT_AU_DAY,
                geo[2] + earth_vel[2] * r / C_LIGHT_AU_DAY,
            ];
            r = cartesian_norm(geo);
        }
        let spherical = to_spherical(geo);
        if of_date || public_of_date {
            let (lon, lat) = to_of_date(spherical[0], spherical[1], jd)?;
            Ok([lon, lat, r])
        } else {
            Ok(spherical)
        }
    };
    // Vulcan's short orbit moves degrees per day geocentrically, so its
    // central difference uses a tighter step than the slow outer bodies.
    let dt_v = if body == VULCAN { 0.01 } else { 0.1 };
    let center = geo_of_date(tt)?;
    let prev = geo_of_date(tt - dt_v)?;
    let next = geo_of_date(tt + dt_v)?;
    if public_of_date {
        Ok(model_stencil(center, prev, next, dt_v, false))
    } else {
        Ok([
            center[0],
            center[1],
            center[2],
            wrap_rate(next[0], prev[0], dt_v),
            (next[1] - prev[1]) / (2.0 * dt_v),
            (next[2] - prev[2]) / (2.0 * dt_v),
        ])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::RefCell;

    #[test]
    fn classic_geocentric_models_use_observation_channel_and_gate_speed_samples() {
        let tt = 2451545.125;
        let no_support = |_: i32, _: f64| -> Result<([f64; 3], [f64; 3])> {
            panic!("TRUEPOS without Topos needs no barycentric velocity")
        };
        let no_aya = |_: f64, _: i32| -> Result<f64> { panic!("tropical output") };
        let no_ut = |_: f64| -> Result<f64> { panic!("tropical output needs no inverse UT") };
        for body in CUPIDO..=ISIS {
            for speed in [false, true] {
                let samples = RefCell::new(Vec::new());
                let observation = |epoch| {
                    samples.borrow_mut().push(epoch);
                    Ok([1.0, 0.1, 0.0])
                };
                let flags = FLG_J2000 | FLG_TRUEPOS | if speed { FLG_SPEED } else { 0 };
                let state = classic_model_state(
                    &no_support,
                    &observation,
                    &no_ut,
                    tt,
                    body,
                    flags,
                    None,
                    &no_aya,
                )
                .unwrap();
                assert!(state.into_iter().all(f64::is_finite));
                assert_eq!(
                    *samples.borrow(),
                    if speed {
                        vec![tt, tt - 1.0, tt + 1.0]
                    } else {
                        vec![tt]
                    },
                    "body {body}"
                );
                if !speed {
                    assert!(state[3..].iter().all(|word| word.to_bits() == 0));
                }
            }
        }
    }

    #[test]
    fn classic_heliocentric_models_bypass_support_and_keep_sidereal_ut() {
        let no_support = |_: i32, _: f64| -> Result<([f64; 3], [f64; 3])> {
            panic!("heliocentric model must not acquire DE support")
        };
        let no_observation = |_: f64| -> Result<[f64; 3]> { panic!("no Earth observation") };
        let tt: f64 = 2451545.125;
        let ut = 2451544.875;
        for body in CUPIDO..=ISIS {
            let calls = RefCell::new(Vec::new());
            let aya = |clock, flags| {
                calls.borrow_mut().push((clock, flags));
                Ok(24.0)
            };
            let flags = FLG_HELCTR | FLG_J2000 | FLG_SIDEREAL;
            let state = classic_model_state(
                &no_support,
                &no_observation,
                &|sample_tt| {
                    assert_eq!(sample_tt.to_bits(), tt.to_bits());
                    Ok(ut)
                },
                tt,
                body,
                flags,
                None,
                &aya,
            )
            .unwrap();
            assert!(state.into_iter().all(f64::is_finite));
            assert!(state[3..].iter().all(|word| word.to_bits() == 0));
            assert_eq!(*calls.borrow(), [(ut, flags)]);
        }
    }

    #[test]
    fn classic_topos_offset_retains_independent_geocenter_read() {
        let calls = RefCell::new(Vec::new());
        let support = |body, tt| {
            let mut calls = calls.borrow_mut();
            calls.push((body, tt));
            let x = if calls.len() == 1 { 1e12 } else { 1e12 + 1.0 };
            Ok(([x, 0.0, 0.0], [0.0; 3]))
        };
        let tt = 2451545.125;
        let offset = source_model_topo_offset(&support, tt, (13.405, 52.52, 34.0)).unwrap();
        assert_eq!(*calls.borrow(), [(EARTH, tt), (EARTH, tt)]);
        // The tiny site x addition rounds to the first synthetic barycentric
        // x. Subtracting the second geocenter retains their exact one-AU
        // difference; cancelling a single reused Earth state would lose it.
        assert_eq!(offset.0[0], -1.0);
    }

    #[test]
    fn classic_topos_missing_earth_precedes_site_geometry() {
        let no_earth = |body, _: f64| -> Result<([f64; 3], [f64; 3])> {
            Err(Error::SpkUnregisteredBody { body })
        };
        assert_eq!(
            source_model_topo_offset(&no_earth, 2451545.125, (f64::NAN, 0.0, 0.0)),
            Err(Error::SpkUnregisteredBody { body: EARTH })
        );
    }

    #[test]
    fn classic_equatorial_sidereal_does_not_consume_ut_or_ayanamsha() {
        let no_support = |_: i32, _: f64| -> Result<([f64; 3], [f64; 3])> {
            panic!("heliocentric output needs no support")
        };
        let no_observation = |_: f64| -> Result<[f64; 3]> { panic!("heliocentric output") };
        let no_ut = |_: f64| -> Result<f64> { panic!("equatorial output needs no UT") };
        let no_aya = |_: f64, _: i32| -> Result<f64> { panic!("equatorial output needs no aya") };
        for body in CUPIDO..=ISIS {
            let state = classic_model_state(
                &no_support,
                &no_observation,
                &no_ut,
                2451545.125,
                body,
                FLG_HELCTR | FLG_J2000 | FLG_SIDEREAL | FLG_EQUATORIAL,
                None,
                &no_aya,
            )
            .unwrap();
            assert!(state.into_iter().all(f64::is_finite));
        }
    }

    #[test]
    fn uranian_leb_reads_earth_before_sun_and_retains_scalar_retardation() {
        let calls = RefCell::new(Vec::new());
        let tt = 2451545.125;
        let support = |body, epoch| {
            calls.borrow_mut().push((body, epoch));
            Ok((
                if body == EARTH {
                    [1.0, 0.0, 0.0]
                } else {
                    [0.0; 3]
                },
                [0.0; 3],
            ))
        };
        assert_eq!(leb_earth_minus_sun(&support, tt).unwrap(), [1.0, 0.0, 0.0]);
        let calls = calls.borrow();
        assert_eq!(calls[0], (EARTH, tt));
        assert_eq!(calls[1], (SUN, tt));
        assert_eq!(calls.len(), 3);
        assert_eq!(calls[2], (EARTH, tt - 1.0 / C_LIGHT_AU_DAY));
    }

    #[test]
    fn uranian_leb_missing_earth_prevents_sun_acquisition() {
        let calls = RefCell::new(Vec::new());
        let support = |body, _: f64| -> Result<([f64; 3], [f64; 3])> {
            calls.borrow_mut().push(body);
            Err(Error::SpkUnregisteredBody { body })
        };
        assert_eq!(
            leb_earth_minus_sun(&support, 2451545.125),
            Err(Error::SpkUnregisteredBody { body: EARTH })
        );
        assert_eq!(*calls.borrow(), [EARTH]);
    }

    #[test]
    fn uranian_bary_sidereal_tail_consumes_original_ut_at_ayanamsha() {
        let tt: f64 = 2451545.125;
        let original_ut: f64 = 2451544.875;
        let clocks = RefCell::new(Vec::new());
        let ayanamsha_clocks = RefCell::new(Vec::new());
        let ut_at = |epoch: f64| {
            clocks.borrow_mut().push(epoch);
            Ok(original_ut)
        };
        let aya = |ut: f64, _: i32| {
            ayanamsha_clocks.borrow_mut().push(ut);
            Ok(24.0)
        };
        let state = finish_uranian_bary_state(
            [120.0, 2.0, 50.0, 0.01, 0.002, 0.003],
            &ut_at,
            tt,
            FLG_BARYCTR | FLG_J2000 | FLG_SIDEREAL | FLG_SPEED,
            &aya,
        )
        .unwrap();
        assert!(state.into_iter().all(f64::is_finite));
        assert_eq!(*clocks.borrow(), [tt]);
        assert_eq!(
            ayanamsha_clocks.borrow()[0].to_bits(),
            original_ut.to_bits()
        );
    }

    #[test]
    fn uranian_bary_nonconsuming_tail_preserves_legacy_frame_words() {
        let no_ut = |_: f64| -> Result<f64> { panic!("nonconsuming tail acquired UT") };
        let no_aya = |_: f64, _: i32| -> Result<f64> { panic!("nonconsuming tail acquired aya") };
        let delta = |_: f64| Ok(0.0);
        let legacy = horizons_output::HorizonsNumerics::new(
            &delta,
            &no_aya,
            0,
            crate::lunar::LunarCoverage::Unenforced,
        );
        let state = [120.0, 2.0, 50.0, 0.01, 0.002, 0.003];
        for tt in [2451545.125, 2460400.25] {
            for extra in [
                0,
                FLG_J2000,
                FLG_NONUT,
                FLG_SPEED,
                FLG_EQUATORIAL | FLG_SPEED,
                FLG_EQUATORIAL | FLG_J2000 | FLG_SPEED,
                FLG_EQUATORIAL | FLG_SIDEREAL | FLG_SPEED,
            ] {
                let flags = FLG_BARYCTR | extra;
                let previous = legacy
                    .hypothetical_state(tt, tt, CUPIDO, flags, state)
                    .unwrap()
                    .0;
                let current = finish_uranian_bary_state(state, &no_ut, tt, flags, &no_aya).unwrap();
                assert_eq!(current.map(f64::to_bits), previous.map(f64::to_bits));
            }
        }
    }

    struct Support {
        queries: RefCell<Vec<(i32, f64)>>,
        available: bool,
    }

    impl StateProvider for Support {
        fn provider_name(&self) -> &'static str {
            "model-support-test"
        }
        fn has_body(&self, body: i32) -> bool {
            self.available && matches!(body, SUN | EARTH)
        }
        fn body_coverage(&self, _: i32) -> Option<(f64, f64)> {
            Some((2400000.0, 2500000.0))
        }
        fn serving_coord_type(&self, _: i32) -> Option<u32> {
            Some(COORD_ICRS_BARY)
        }
        fn eval_raw(&self, body: i32, tt: f64) -> Result<(Vec<f64>, Vec<f64>)> {
            self.queries.borrow_mut().push((body, tt));
            if !self.has_body(body) {
                return Err(Error::LebBodyNotFound { body_id: body });
            }
            let position = if body == EARTH {
                vec![1.0, 0.2, 0.1]
            } else {
                vec![0.01, 0.0, 0.0]
            };
            Ok((position, vec![0.0; 3]))
        }
        fn eval_raw_two_part(&self, _: i32, _: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)> {
            panic!("public model support queries are scalar TT")
        }
        fn serving_entry_snapshot(&self, _: i32) -> Result<EntrySnapshot> {
            Ok(EntrySnapshot {
                coord_type: COORD_ICRS_BARY,
            })
        }
    }

    #[test]
    fn curated_model_uses_scalar_retarded_support_and_ignores_parallax() {
        let provider = Support {
            queries: RefCell::new(vec![]),
            available: true,
        };
        let no_clock = |_: f64| -> Result<f64> {
            panic!("curated fallback never queries configured model TT")
        };
        let no_aya =
            |_: f64, _: i32| -> Result<f64> { panic!("tropical output never queries ayanamsha") };
        let tt = 2451545.125;
        let calculate = |flags, site| {
            serving_model_position(
                &provider, 2451545.0, tt, PHOLUS, flags, site, &no_clock, &no_aya, 0,
            )
            .unwrap()
        };
        let base = calculate(0, None).0;
        let queries = provider.queries.borrow().clone();
        assert_eq!(queries[0], (SUN, tt));
        assert_eq!(queries[1], (EARTH, tt));
        assert!(
            queries[2..]
                .iter()
                .all(|&(body, jd)| body == EARTH && jd < tt)
        );
        assert!(queries.len() > 2);
        assert!(base[3..].iter().all(|x| x.to_bits() == 0));
        for site in [Some((13.405, 52.52, 34.0)), Some((151.21, -33.86, 43.0))] {
            let flags = FLG_TOPOCTR | FLG_TRUEPOS | FLG_NOABERR;
            let (topo, ret) = calculate(flags, site);
            assert_eq!(topo.map(f64::to_bits), base.map(f64::to_bits));
            assert_eq!(ret, flags);
        }
        assert_eq!(calculate(FLG_BARYCTR, None).0, base);
    }

    #[test]
    fn curated_heliocentric_bypasses_support_and_sidereal_uses_original_ut() {
        let provider = Support {
            queries: RefCell::new(vec![]),
            available: false,
        };
        let no_clock = |_: f64| -> Result<f64> {
            panic!("curated fallback never queries configured model TT")
        };
        let aya_calls = RefCell::new(vec![]);
        let aya = |ut: f64, flags| -> Result<f64> {
            aya_calls.borrow_mut().push((ut, flags));
            Ok(24.0)
        };
        let flags = FLG_HELCTR | FLG_SIDEREAL | FLG_SPEED;
        let ut = 2451545.0;
        let (state, ret) = serving_model_position(
            &provider,
            ut,
            2451545.125,
            PHOLUS,
            flags,
            None,
            &no_clock,
            &aya,
            0,
        )
        .unwrap();
        assert!(state.into_iter().all(f64::is_finite));
        assert_eq!(ret, flags);
        assert_eq!(
            &*aya_calls.borrow(),
            &[
                (ut, flags),
                (ut - 1.0 / 86400.0, flags),
                (ut + 1.0 / 86400.0, flags)
            ]
        );
        assert!(provider.queries.borrow().is_empty());
    }

    #[test]
    fn curated_missing_earth_alias_precedes_sun_and_state_acquisition() {
        let provider = Support {
            queries: RefCell::new(vec![]),
            available: false,
        };
        let unused_clock = |_: f64| -> Result<f64> { panic!("no model clock") };
        let unused_aya = |_: f64, _: i32| -> Result<f64> { panic!("no ayanamsha") };
        assert_eq!(
            serving_model_position(
                &provider,
                2451545.0,
                2451545.125,
                PHOLUS,
                0,
                None,
                &unused_clock,
                &unused_aya,
                0,
            ),
            Err(Error::UnknownBody { body_id: EARTH })
        );
        assert!(provider.queries.borrow().is_empty());
    }

    #[test]
    fn unknown_model_and_heliocentric_models_need_no_support_acquisition() {
        let provider = Support {
            queries: RefCell::new(vec![]),
            available: false,
        };
        let no_clock =
            |_: f64| -> Result<f64> { panic!("model branch must not query fast-path Delta-T") };
        let no_aya =
            |_: f64, _: i32| -> Result<f64> { panic!("tropical output must not query ayanamsha") };
        assert!(matches!(
            serving_model_position(
                &provider,
                2451545.0,
                2451545.125,
                NIBIRU,
                FLG_HELCTR,
                None,
                &no_clock,
                &no_aya,
                0
            ),
            Err(Error::UnknownBody { body_id: NIBIRU })
        ));
        for body in CUPIDO..=PROSERPINA {
            if matches!(body, NIBIRU | WHITE_MOON) {
                continue;
            }
            let (state, _) = serving_model_position(
                &provider,
                2451545.0,
                2451545.125,
                body,
                FLG_HELCTR | FLG_J2000,
                None,
                &no_clock,
                &no_aya,
                0,
            )
            .unwrap();
            assert!(state.into_iter().all(f64::is_finite), "body {body}");
            assert_eq!(&state[3..], &[0.0; 3], "body {body}");
        }
        assert!(provider.queries.borrow().is_empty());
    }

    #[test]
    fn waldemath_anchor_keeps_configured_delta_t_separate_from_request_clock() {
        let provider = Support {
            queries: RefCell::new(vec![]),
            available: false,
        };
        let anchors = RefCell::new(vec![]);
        let delta_t = |ut| -> Result<f64> {
            anchors.borrow_mut().push(ut);
            Ok(0.002)
        };
        let no_aya =
            |_: f64, _: i32| -> Result<f64> { panic!("tropical output must not query ayanamsha") };
        let first = serving_model_position(
            &provider,
            2451544.0,
            2451545.125,
            WALDEMATH,
            FLG_TRUEPOS | FLG_NONUT,
            None,
            &delta_t,
            &no_aya,
            0,
        )
        .unwrap()
        .0;
        let second = serving_model_position(
            &provider,
            2451543.0,
            2451545.125,
            WALDEMATH,
            FLG_TRUEPOS | FLG_NONUT,
            None,
            &delta_t,
            &no_aya,
            0,
        )
        .unwrap()
        .0;
        assert_eq!(first, second);
        assert_eq!(
            &*anchors.borrow(),
            &[crate::hypothetical::WALDEMATH_ANCHOR_JD_UT; 2]
        );
        assert!(provider.queries.borrow().is_empty());
    }

    #[test]
    fn topocentric_symbolic_point_requires_site_and_uses_independent_clock() {
        let provider = Support {
            queries: RefCell::new(vec![]),
            available: true,
        };
        let no_clock =
            |_: f64| -> Result<f64> { panic!("White Moon never needs configured Delta-T") };
        let no_aya =
            |_: f64, _: i32| -> Result<f64> { panic!("tropical output must not query ayanamsha") };
        let flags = FLG_TRUEPOS | FLG_TOPOCTR | FLG_NONUT;
        assert!(matches!(
            serving_model_position(
                &provider,
                2451545.0,
                2451545.125,
                WHITE_MOON,
                flags,
                None,
                &no_clock,
                &no_aya,
                0
            ),
            Err(Error::ConfigurationMissing {
                missing_config: "observer_location",
                ..
            })
        ));
        let first = serving_model_position(
            &provider,
            2451545.0,
            2451545.125,
            WHITE_MOON,
            flags,
            Some((13.405, 52.52, 34.0)),
            &no_clock,
            &no_aya,
            0,
        )
        .unwrap()
        .0;
        let second = serving_model_position(
            &provider,
            2451535.0,
            2451545.125,
            WHITE_MOON,
            flags,
            Some((13.405, 52.52, 34.0)),
            &no_clock,
            &no_aya,
            0,
        )
        .unwrap()
        .0;
        assert_eq!(first, second);
        assert!(
            provider
                .queries
                .borrow()
                .iter()
                .all(|&(body, epoch)| body == EARTH && epoch == 2451545.125)
        );
    }
}

#[cfg(test)]
mod source_time_domain_tests {
    use super::*;
    use std::cell::Cell;

    #[test]
    fn internal_classic_time_reaches_consumed_earth_observation() {
        let observations = Cell::new(0);
        let earth = |tt: f64| {
            assert!(tt.is_nan());
            observations.set(observations.get() + 1);
            Err(Error::SpkOutOfRange {
                target: 399,
                et: f64::NAN,
                start: -1.0,
                end: 1.0,
            })
        };
        let unused_state =
            |_: i32, _: f64| -> Result<([f64; 3], [f64; 3])> { panic!("unused velocity") };
        let unused_ut = |_: f64| -> Result<f64> { panic!("unused inverse") };
        let unused_aya = |_: f64, _: i32| -> Result<f64> { panic!("unused sidereal") };
        for body in [CUPIDO, ISIS] {
            let result = classic_model_state_source_time(
                &unused_state,
                &earth,
                &unused_ut,
                f64::NAN,
                body,
                FLG_SWIEPH,
                None,
                &unused_aya,
            );
            assert!(matches!(
                result,
                Err(Error::SpkOutOfRange { target: 399, .. })
            ));
        }
        assert_eq!(observations.get(), 2);
        assert!(matches!(
            classic_model_state(
                &unused_state,
                &earth,
                &unused_ut,
                f64::NAN,
                CUPIDO,
                FLG_SWIEPH,
                None,
                &unused_aya
            ),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert_eq!(observations.get(), 2);
    }

    #[test]
    fn predicted_internal_time_observes_before_missing_model() {
        let earth = |tt: f64| -> Result<RetardedEarthObservation> {
            assert!(tt.is_nan());
            Err(Error::SpkOutOfRange {
                target: 399,
                et: f64::NAN,
                start: -1.0,
                end: 1.0,
            })
        };
        let unused = |_: i32, _: f64| -> Result<([f64; 3], [f64; 3])> {
            panic!("observation failure must win")
        };
        assert!(matches!(
            predicted_center_state_source_time(
                &unused,
                &earth,
                false,
                f64::NAN,
                NIBIRU,
                FLG_SWIEPH
            ),
            Err(Error::SpkOutOfRange { target: 399, .. })
        ));
    }
}

/// Per-call owned model variant of the corresponding static reduction.
/// Model callbacks run only where its geometric state is consumed. They return
/// heliocentric degrees/AU and per-day rates; all support clocks, frames, flags,
/// stencils and callback errors keep the same ordering. None uses static rows.
#[allow(clippy::too_many_arguments)]
pub fn serving_model_position_with_models<P: StateProvider>(
    models: Option<&AnalyticalModelCallback<'_>>,
    provider: &P,
    request_ut1: f64,
    global_tt: f64,
    body: i32,
    flags: i32,
    site: Option<(f64, f64, f64)>,
    model_deltat: &dyn Fn(f64) -> Result<f64>,
    ayanamsha_at: &dyn Fn(f64, i32) -> Result<f64>,
    _sidereal_mode: i32,
) -> Result<(PositionResult, i32)> {
    for value in [request_ut1, global_tt] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    if crate::minor_bodies::minor_body_elements(body).is_some() {
        let mut state = curated_keplerian_state(provider, global_tt, body, flags)?;
        if flags & FLG_SIDEREAL != 0 && flags & FLG_EQUATORIAL == 0 {
            state = horizons_output::local_sidereal(state, request_ut1, flags, ayanamsha_at)?;
        }
        state = horizons_output::local_convert(state, global_tt, flags, false)?;
        return Ok((apply_output_flags(state, flags), flags));
    }
    if !(CUPIDO..=WALDEMATH).contains(&body) || body == NIBIRU {
        return Err(Error::UnknownBody { body_id: body });
    }
    let state_at = |id, tt| icrs_state(provider, id, tt);
    let earth_helio_at = |tt| observed_earth_minus_sun(&state_at, tt);
    let uranian_earth_at = |tt| leb_earth_minus_sun(&state_at, tt);
    let already_requested_frame = (CUPIDO..=ISIS).contains(&body);
    let mut rate_origin = ModelRateOrigin::LegacyStencil;
    let state = if (CUPIDO..=POSEIDON).contains(&body) {
        uranian_center_state_with_models(
            models,
            &state_at,
            &uranian_earth_at,
            global_tt,
            body,
            flags,
            false,
        )?
    } else if body == ISIS {
        isis_center_state_with_models(models, &state_at, &earth_helio_at, global_tt, flags)?
    } else if matches!(body, WHITE_MOON | WALDEMATH) {
        let symbolic = if body == WALDEMATH {
            let anchor_ut = crate::hypothetical::WALDEMATH_ANCHOR_JD_UT;
            let anchor_tt = anchor_ut + model_deltat(anchor_ut)?;
            let mean_lon = crate::hypothetical::WALDEMATH_ANCHOR_APPARENT_LON_DEG
                - crate::erfa::nut06a(anchor_tt)?.0.to_degrees();
            crate::hypothetical::waldemath_state_with_anchor(global_tt, anchor_tt, mean_lon)?
        } else {
            crate::hypothetical::calc_white_moon_position(global_tt, false)?
        };
        symbolic_center_state(&state_at, &earth_helio_at, global_tt, flags, symbolic)?
    } else if models.is_none() && native_predicted_rates_available(provider, body, flags) {
        let (state, origin) = native_provider_predicted_state(provider, global_tt, body, flags)?;
        rate_origin = origin;
        state
    } else {
        predicted_center_state_with_models(
            models,
            &state_at,
            &|tt| observed_earth_with_origin(&state_at, tt),
            true,
            global_tt,
            body,
            flags,
        )?
    };
    let state = finish_model_state_with_rate_origin(
        &state_at,
        state,
        &|_| Ok(request_ut1),
        global_tt,
        flags,
        already_requested_frame,
        site,
        ayanamsha_at,
        rate_origin,
    )?;
    Ok((apply_output_flags(state, flags), flags))
}

/// Per-call owned model variant of the corresponding static reduction.
/// Model callbacks run only where its geometric state is consumed. They return
/// heliocentric degrees/AU and per-day rates; all support clocks, frames, flags,
/// stencils and callback errors keep the same ordering. None uses static rows.
#[allow(clippy::too_many_arguments)]
pub(super) fn classic_model_state_with_models(
    models: Option<&AnalyticalModelCallback<'_>>,
    state_at: &dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
    earth_helio_at: &dyn Fn(f64) -> Result<[f64; 3]>,
    ut_at: &dyn Fn(f64) -> Result<f64>,
    tt: f64,
    body: i32,
    flags: i32,
    site: Option<(f64, f64, f64)>,
    aya: &dyn Fn(f64, i32) -> Result<f64>,
) -> Result<PositionResult> {
    if !tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: tt });
    }
    classic_model_state_core_with_models(
        models,
        state_at,
        earth_helio_at,
        ut_at,
        tt,
        body,
        flags,
        site,
        aya,
        false,
    )
}

/// Per-call owned model variant of the corresponding static reduction.
/// Model callbacks run only where its geometric state is consumed. They return
/// heliocentric degrees/AU and per-day rates; all support clocks, frames, flags,
/// stencils and callback errors keep the same ordering. None uses static rows.
#[allow(clippy::too_many_arguments)]
pub(super) fn classic_model_state_source_time_with_models(
    models: Option<&AnalyticalModelCallback<'_>>,
    state_at: &dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
    earth_helio_at: &dyn Fn(f64) -> Result<[f64; 3]>,
    ut_at: &dyn Fn(f64) -> Result<f64>,
    tt: f64,
    body: i32,
    flags: i32,
    site: Option<(f64, f64, f64)>,
    aya: &dyn Fn(f64, i32) -> Result<f64>,
) -> Result<PositionResult> {
    if tt.is_finite() {
        return classic_model_state_with_models(
            models,
            state_at,
            earth_helio_at,
            ut_at,
            tt,
            body,
            flags,
            site,
            aya,
        );
    }
    classic_model_state_core_with_models(
        models,
        state_at,
        earth_helio_at,
        ut_at,
        tt,
        body,
        flags,
        site,
        aya,
        true,
    )
}

/// Per-call owned model variant of the corresponding static reduction.
/// Model callbacks run only where its geometric state is consumed. They return
/// heliocentric degrees/AU and per-day rates; all support clocks, frames, flags,
/// stencils and callback errors keep the same ordering. None uses static rows.
#[allow(clippy::too_many_arguments)]
fn classic_model_state_core_with_models(
    models: Option<&AnalyticalModelCallback<'_>>,
    state_at: &dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
    earth_helio_at: &dyn Fn(f64) -> Result<[f64; 3]>,
    ut_at: &dyn Fn(f64) -> Result<f64>,
    tt: f64,
    body: i32,
    flags: i32,
    site: Option<(f64, f64, f64)>,
    aya: &dyn Fn(f64, i32) -> Result<f64>,
    source_time: bool,
) -> Result<PositionResult> {
    let state = if (CUPIDO..=POSEIDON).contains(&body) {
        uranian_center_state_with_models(
            models,
            state_at,
            earth_helio_at,
            tt,
            body,
            flags,
            source_time,
        )?
    } else if body == ISIS {
        isis_center_state_with_models(models, state_at, earth_helio_at, tt, flags)?
    } else {
        return Err(Error::UnknownBody { body_id: body });
    };
    finish_model_state(state_at, state, ut_at, tt, flags, true, site, aya)
}

/// Per-call owned model variant of the corresponding static reduction.
/// Model callbacks run only where its geometric state is consumed. They return
/// heliocentric degrees/AU and per-day rates; all support clocks, frames, flags,
/// stencils and callback errors keep the same ordering. None uses static rows.
#[allow(clippy::too_many_arguments)]
fn uranian_center_state_with_models(
    models: Option<&AnalyticalModelCallback<'_>>,
    state_at: &dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
    earth_helio_at: &dyn Fn(f64) -> Result<[f64; 3]>,
    tt: f64,
    body: i32,
    flags: i32,
    source_time: bool,
) -> Result<PositionResult> {
    if flags & FLG_BARYCTR != 0 {
        let mut state =
            uranian_bary_time_state_with_models(models, state_at, tt, body, source_time)?;
        if flags & FLG_SPEED == 0 {
            state[3..].fill(0.0);
        }
        if flags & FLG_J2000 == 0 {
            state = precess_model_state(state, tt, flags, false)?;
        }
        return Ok(state);
    }
    let sample = |jd| -> Result<[f64; 3]> {
        let h = crate::positions::model_acquisition::selected_state(
            models,
            AnalyticalModel::Uranian(body),
            jd,
            source_time,
            || {
                if source_time {
                    crate::hypothetical::calc_uranian_planet_source_time(body, jd)
                } else {
                    crate::hypothetical::calc_uranian_planet(body, jd)
                }
            },
        )?;
        let helio = model_xyz(h.lon_deg, h.lat_deg, h.dist_au);
        if flags & FLG_HELCTR != 0 {
            return Ok([h.lon_deg, h.lat_deg, h.dist_au]);
        }
        let earth = rotate_icrs_to_iau1976_ecliptic_j2000(earth_helio_at(jd)?);
        let mut geo = std::array::from_fn(|i| helio[i] - earth[i]);
        if flags & (FLG_TRUEPOS | FLG_NOABERR) == 0 {
            let velocity = rotate_icrs_to_iau1976_ecliptic_j2000(state_at(EARTH, jd)?.1);
            geo = bradley_aberration(geo, velocity);
        }
        Ok(model_spherical(geo))
    };
    let center = sample(tt)?;
    // Centered branches always compute their model state; geocentric speed
    // has an explicit acquisition gate in the public contract.
    let mut state = if flags & FLG_SPEED != 0 || flags & (FLG_HELCTR | FLG_BARYCTR) != 0 {
        model_stencil(center, sample(tt - 1.0)?, sample(tt + 1.0)?, 1.0, false)
    } else {
        [center[0], center[1], center[2], 0.0, 0.0, 0.0]
    };
    if flags & FLG_SPEED == 0 {
        state[3..].fill(0.0);
    }
    if flags & FLG_J2000 == 0 {
        state = precess_model_state(state, tt, flags, false)?;
    }
    Ok(state)
}

/// Per-call owned model variant of the corresponding static reduction.
/// Model callbacks run only where its geometric state is consumed. They return
/// heliocentric degrees/AU and per-day rates; all support clocks, frames, flags,
/// stencils and callback errors keep the same ordering. None uses static rows.
#[allow(clippy::too_many_arguments)]
pub(super) fn uranian_bary_source_time_state_with_models(
    models: Option<&AnalyticalModelCallback<'_>>,
    state_at: &dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
    tt: f64,
    body: i32,
) -> Result<PositionResult> {
    uranian_bary_time_state_with_models(models, state_at, tt, body, true)
}

/// Per-call owned model variant of the corresponding static reduction.
/// Model callbacks run only where its geometric state is consumed. They return
/// heliocentric degrees/AU and per-day rates; all support clocks, frames, flags,
/// stencils and callback errors keep the same ordering. None uses static rows.
#[allow(clippy::too_many_arguments)]
fn uranian_bary_time_state_with_models(
    models: Option<&AnalyticalModelCallback<'_>>,
    state_at: &dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
    tt: f64,
    body: i32,
    source_time: bool,
) -> Result<PositionResult> {
    let sample = |jd| -> Result<[f64; 3]> {
        let h = crate::positions::model_acquisition::selected_state(
            models,
            AnalyticalModel::Uranian(body),
            jd,
            source_time,
            || {
                if source_time {
                    crate::hypothetical::calc_uranian_planet_source_time(body, jd)
                } else {
                    crate::hypothetical::calc_uranian_planet(body, jd)
                }
            },
        )?;
        let helio = model_xyz(h.lon_deg, h.lat_deg, h.dist_au);
        let sun = rotate_icrs_to_iau1976_ecliptic_j2000(state_at(SUN, jd)?.0);
        Ok(model_spherical(std::array::from_fn(|i| helio[i] + sun[i])))
    };
    Ok(model_stencil(
        sample(tt)?,
        sample(tt - 1.0)?,
        sample(tt + 1.0)?,
        1.0,
        false,
    ))
}

/// Per-call owned model variant of the corresponding static reduction.
/// Model callbacks run only where its geometric state is consumed. They return
/// heliocentric degrees/AU and per-day rates; all support clocks, frames, flags,
/// stencils and callback errors keep the same ordering. None uses static rows.
#[allow(clippy::too_many_arguments)]
fn isis_center_state_with_models(
    models: Option<&AnalyticalModelCallback<'_>>,
    state_at: &dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
    earth_helio_at: &dyn Fn(f64) -> Result<[f64; 3]>,
    tt: f64,
    flags: i32,
) -> Result<PositionResult> {
    let helio_lt = |jd| -> Result<f64> {
        if flags & FLG_TRUEPOS != 0 {
            return Ok(0.0);
        }
        let mut lt = 0.0;
        for _ in 0..3 {
            lt = crate::positions::model_acquisition::selected_state(
                models,
                AnalyticalModel::Transpluto,
                jd - lt,
                false,
                || crate::hypothetical::calc_transpluto(jd - lt),
            )?
            .dist_au
                / C_LIGHT_AU_DAY;
        }
        Ok(lt)
    };
    if flags & FLG_HELCTR != 0 {
        let model_tt = tt - helio_lt(tt)?;
        let h = crate::positions::model_acquisition::selected_state(
            models,
            AnalyticalModel::Transpluto,
            model_tt,
            false,
            || crate::hypothetical::calc_transpluto(model_tt),
        )?;
        let mut state = [h.lon_deg, h.lat_deg, h.dist_au, h.dlon, h.dlat, h.ddist];
        if flags & FLG_SPEED == 0 {
            state[3..].fill(0.0);
        }
        if flags & FLG_J2000 == 0 {
            state = precess_model_state(state, tt, flags, false)?;
        }
        return Ok(state);
    }
    let bary = flags & FLG_BARYCTR != 0;
    let sample = |jd| -> Result<[f64; 3]> {
        if bary {
            let model_tt = jd - helio_lt(jd)?;
            let h = crate::positions::model_acquisition::selected_state(
                models,
                AnalyticalModel::Transpluto,
                model_tt,
                false,
                || crate::hypothetical::calc_transpluto(model_tt),
            )?;
            let helio = model_xyz(h.lon_deg, h.lat_deg, h.dist_au);
            let sun = rotate_icrs_to_iau1976_ecliptic_j2000(state_at(SUN, jd)?.0);
            let mut result = model_spherical(std::array::from_fn(|i| helio[i] + sun[i]));
            if flags & FLG_J2000 == 0 {
                (result[0], result[1]) = crate::hypothetical::precess_ecliptic(
                    result[0],
                    result[1],
                    crate::hypothetical::J2000_JD,
                    jd,
                )?;
            }
            return Ok(result);
        }
        let earth = rotate_icrs_to_iau1976_ecliptic_j2000(earth_helio_at(jd)?);
        let mut lt = 0.0;
        let mut geo = [0.0; 3];
        for _ in 0..3 {
            let h = crate::positions::model_acquisition::selected_state(
                models,
                AnalyticalModel::Transpluto,
                jd - lt,
                false,
                || crate::hypothetical::calc_transpluto(jd - lt),
            )?;
            let target = model_xyz(h.lon_deg, h.lat_deg, h.dist_au);
            geo = std::array::from_fn(|i| target[i] - earth[i]);
            if flags & FLG_TRUEPOS != 0 {
                break;
            }
            lt = model_norm(geo) / C_LIGHT_AU_DAY;
        }
        let radius = model_norm(geo);
        if flags & (FLG_TRUEPOS | FLG_NOABERR) == 0 && radius > 0.0 {
            let velocity = rotate_icrs_to_iau1976_ecliptic_j2000(state_at(EARTH, jd)?.1);
            geo = std::array::from_fn(|i| geo[i] + velocity[i] * radius / C_LIGHT_AU_DAY);
        }
        Ok(model_spherical(geo))
    };
    let center = sample(tt)?;
    let step = if bary { 0.01 } else { 1.0 };
    let mut state = if flags & FLG_SPEED != 0 {
        model_stencil(center, sample(tt - step)?, sample(tt + step)?, step, bary)
    } else {
        [center[0], center[1], center[2], 0.0, 0.0, 0.0]
    };
    if !bary && flags & FLG_J2000 == 0 {
        state = precess_model_state(state, tt, flags, false)?;
    }
    Ok(state)
}

/// Per-call owned model variant of the corresponding static reduction.
/// Model callbacks run only where its geometric state is consumed. They return
/// heliocentric degrees/AU and per-day rates; all support clocks, frames, flags,
/// stencils and callback errors keep the same ordering. None uses static rows.
#[allow(clippy::too_many_arguments)]
pub(super) fn predicted_center_state_with_models(
    models: Option<&AnalyticalModelCallback<'_>>,
    state_at: &dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
    earth_observation_at: &dyn Fn(f64) -> Result<RetardedEarthObservation>,
    public_of_date: bool,
    tt: f64,
    body: i32,
    flags: i32,
) -> Result<PositionResult> {
    if !tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: tt });
    }
    predicted_center_state_source_time_with_models(
        models,
        state_at,
        earth_observation_at,
        public_of_date,
        tt,
        body,
        flags,
    )
}

/// Per-call owned model variant of the corresponding static reduction.
/// Model callbacks run only where its geometric state is consumed. They return
/// heliocentric degrees/AU and per-day rates; all support clocks, frames, flags,
/// stencils and callback errors keep the same ordering. None uses static rows.
#[allow(clippy::too_many_arguments)]
pub(super) fn predicted_center_state_source_time_with_models(
    models: Option<&AnalyticalModelCallback<'_>>,
    state_at: &dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
    earth_observation_at: &dyn Fn(f64) -> Result<RetardedEarthObservation>,
    public_of_date: bool,
    tt: f64,
    body: i32,
    flags: i32,
) -> Result<PositionResult> {
    predicted_center_state_source_time_with_models_and_rates(
        models,
        state_at,
        earth_observation_at,
        public_of_date,
        tt,
        body,
        flags,
        None,
        None,
    )
}

#[allow(clippy::too_many_arguments)]
/// Shared model reduction with an optional caller-owned support-domain sampler.
/// `None` retains the original callbacks, positions, errors and rate stencil.
pub(super) fn predicted_center_state_source_time_with_models_and_rates(
    models: Option<&AnalyticalModelCallback<'_>>,
    state_at: &dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
    earth_observation_at: &dyn Fn(f64) -> Result<RetardedEarthObservation>,
    public_of_date: bool,
    tt: f64,
    body: i32,
    flags: i32,
    spherical_rates: Option<&PredictedRateAt<'_>>,
    central_complete: Option<&dyn Fn()>,
) -> Result<PositionResult> {
    predicted_center_state_source_time_with_models_rates_and_positions(
        models,
        state_at,
        earth_observation_at,
        public_of_date,
        tt,
        body,
        flags,
        spherical_rates,
        central_complete,
        None,
    )
}

#[allow(clippy::too_many_arguments)]
fn predicted_center_state_source_time_with_models_rates_and_positions(
    models: Option<&AnalyticalModelCallback<'_>>,
    state_at: &PredictedStateAt<'_>,
    earth_observation_at: &dyn Fn(f64) -> Result<RetardedEarthObservation>,
    public_of_date: bool,
    tt: f64,
    body: i32,
    flags: i32,
    spherical_rates: Option<&PredictedRateAt<'_>>,
    central_complete: Option<&dyn Fn()>,
    deflector_at: Option<&PredictedDeflectorAt<'_>>,
) -> Result<PositionResult> {
    let of_date = body == VULCAN || body == PROSERPINA;
    let model_at = |jd: f64| -> Result<crate::hypothetical::HypoState> {
        if body == VULCAN {
            crate::positions::model_acquisition::selected_state(
                models,
                AnalyticalModel::Vulcan,
                jd,
                false,
                || crate::hypothetical::calc_vulcan(jd),
            )
        } else if body == PROSERPINA {
            crate::hypothetical::calc_proserpina(jd)
        } else {
            crate::hypothetical::calc_fictitious_position(body, jd)
        }
    };
    // Heliocentric J2000-ecliptic Cartesian model place at `jd`: Vulcan and
    // Proserpina use the equinox of date, so their angles rotate to J2000
    // before the Earth subtraction is frame-consistent.
    let helio_xyz_j2000 = |jd: f64| -> Result<[f64; 3]> {
        let state = model_at(jd)?;
        let (lon, lat) = if of_date {
            crate::hypothetical::precess_ecliptic(
                state.lon_deg,
                state.lat_deg,
                jd,
                crate::hypothetical::J2000_JD,
            )?
        } else {
            (state.lon_deg, state.lat_deg)
        };
        let (lon_r, lat_r) = (lon.to_radians(), lat.to_radians());
        let cl = lat_r.cos();
        Ok([
            state.dist_au * cl * lon_r.cos(),
            state.dist_au * cl * lon_r.sin(),
            state.dist_au * lat_r.sin(),
        ])
    };
    let cartesian_norm =
        |xyz: [f64; 3]| -> f64 { (xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]).sqrt() };
    let to_spherical = |xyz: [f64; 3]| -> [f64; 3] {
        let r = cartesian_norm(xyz);
        [
            xyz[1].atan2(xyz[0]).to_degrees().rem_euclid(360.0),
            if r > 0.0 {
                (xyz[2] / r).clamp(-1.0, 1.0).asin().to_degrees()
            } else {
                0.0
            },
            r,
        ]
    };
    // Centered longitude rate with the angular difference reduced to a half turn.
    let wrap_rate = |next_lon: f64, prev_lon: f64, h: f64| -> f64 {
        ((next_lon - prev_lon + 180.0).rem_euclid(360.0) - 180.0) / (2.0 * h)
    };
    let truepos = flags & FLG_TRUEPOS != 0;
    // Sun-to-body light time in days at `jd` (0 under TRUEPOS), iterated
    // like the geocentric case: the heliocentric distance is itself
    // evaluated at the retarded epoch.
    let helio_light_time = |jd: f64| -> Result<f64> {
        if truepos {
            return Ok(0.0);
        }
        let mut lt = 0.0;
        for _ in 0..3 {
            lt = cartesian_norm(helio_xyz_j2000(jd - lt)?) / C_LIGHT_AU_DAY;
        }
        Ok(lt)
    };
    // Rotate a J2000-ecliptic place to the mean ecliptic at TT Julian Day `jd`.
    let to_of_date = |lon: f64, lat: f64, jd: f64| -> Result<(f64, f64)> {
        crate::hypothetical::precess_ecliptic(lon, lat, crate::hypothetical::J2000_JD, jd)
    };
    if flags & FLG_BARYCTR != 0 && flags & FLG_HELCTR == 0 {
        // Barycentric place: the model is heliocentric, so the Sun's
        // barycentric vector is added in J2000-ecliptic axes.
        let place_j2000 = |jd: f64| -> Result<[f64; 3]> {
            let lt = helio_light_time(jd)?;
            let helio = helio_xyz_j2000(jd - lt)?;
            let (sun, _) = state_at(SUN, jd)?;
            let sun = rotate_icrs_to_iau1976_ecliptic_j2000(sun);
            Ok(to_spherical([
                helio[0] + sun[0],
                helio[1] + sun[1],
                helio[2] + sun[2],
            ]))
        };
        let center = place_j2000(tt)?;
        let (lon, lat) = if of_date || public_of_date {
            to_of_date(center[0], center[1], tt)?
        } else {
            (center[0], center[1])
        };
        if public_of_date && flags & FLG_SPEED == 0 {
            return Ok([lon, lat, center[2], 0.0, 0.0, 0.0]);
        }
        let h = 0.01;
        let prev = place_j2000(tt - h)?;
        let next = place_j2000(tt + h)?;
        let (p_lon, p_lat) = if of_date || public_of_date {
            to_of_date(prev[0], prev[1], tt - h)?
        } else {
            (prev[0], prev[1])
        };
        let (n_lon, n_lat) = if of_date || public_of_date {
            to_of_date(next[0], next[1], tt + h)?
        } else {
            (next[0], next[1])
        };
        return Ok([
            lon,
            lat,
            center[2],
            wrap_rate(n_lon, p_lon, h),
            (n_lat - p_lat) / (2.0 * h),
            (next[2] - prev[2]) / (2.0 * h),
        ]);
    }
    if flags & FLG_HELCTR != 0 {
        // Apparent heliocentric place: the model at the Sun-to-body
        // retarded epoch. TRUEPOS zeroes that light time. Speeds
        // central-difference the reported quantity so they stay the
        // derivative of the reported position.
        let h = 0.01;
        let reported = |jd: f64| -> Result<[f64; 3]> {
            let state = model_at(jd - helio_light_time(jd)?)?;
            // Vulcan/Proserpina models already use the equinox of date; the
            // predicted-planet propagation already returns J2000. Both match
            // the tail convention without further rotation here.
            Ok([state.lon_deg, state.lat_deg, state.dist_au])
        };
        let center = reported(tt)?;
        // Vulcan/Proserpina models already use the equinox of date; the
        // predicted-planet propagation already returns J2000. Both match
        // the tail convention without further rotation here.
        if public_of_date && flags & FLG_SPEED == 0 {
            let mut state = [center[0], center[1], center[2], 0.0, 0.0, 0.0];
            if !of_date {
                state = precess_model_state(state, tt, flags, false)?;
            }
            return Ok(state);
        }
        let prev = reported(tt - h)?;
        let next = reported(tt + h)?;
        let mut state = [
            center[0],
            center[1],
            center[2],
            wrap_rate(next[0], prev[0], h),
            (next[1] - prev[1]) / (2.0 * h),
            (next[2] - prev[2]) / (2.0 * h),
        ];
        if public_of_date && !of_date {
            state = precess_model_state(state, tt, flags, false)?;
        }
        return Ok(state);
    }
    // Construct a geocentric apparent place in the mean ecliptic of date.
    // Retard target and Earth, then apply PPN deflection and annual aberration
    // unless geometric flags suppress them. Reevaluate the observer at samples.
    let noaberr = flags & (FLG_TRUEPOS | FLG_NOABERR) != 0;
    let nogdefl = flags & (FLG_TRUEPOS | FLG_NOGDEFL) != 0;
    let geo_of_date = |jd: f64| -> Result<[f64; 3]> {
        let observation = earth_observation_at(jd)?;
        let earth_j2000 = rotate_icrs_to_iau1976_ecliptic_j2000(observation.relative);
        let (_, earth_vel_icrs) = state_at(EARTH, jd)?;
        let earth_vel = rotate_icrs_to_iau1976_ecliptic_j2000(earth_vel_icrs);
        let mut lt = 0.0;
        let mut geo = [0.0; 3];
        for _ in 0..3 {
            let target = helio_xyz_j2000(jd - lt)?;
            geo = [
                target[0] - earth_j2000[0],
                target[1] - earth_j2000[1],
                target[2] - earth_j2000[2],
            ];
            if truepos {
                break;
            }
            lt = cartesian_norm(geo) / C_LIGHT_AU_DAY;
        }
        let mut r = cartesian_norm(geo);
        if !nogdefl && r > 0.0 {
            // All PPN vectors use the model's bias-free IAU 1976 ecliptic
            // axes and the same Sun origin fixed at this reception epoch.
            geo = apply_gravitational_deflection(
                geo,
                earth_j2000,
                jd,
                if public_of_date {
                    r / C_LIGHT_AU_DAY
                } else {
                    lt
                },
                |id, epoch| {
                    let position = match deflector_at {
                        Some(position_at) => position_at(id, epoch)?,
                        None => state_at(id, epoch)?.0,
                    };
                    deflector_in_frame(position, observation.sun_reception, |relative| {
                        Ok(rotate_icrs_to_iau1976_ecliptic_j2000(relative))
                    })
                },
            )?;
            r = cartesian_norm(geo);
        }
        if !noaberr && r > 0.0 {
            geo = [
                geo[0] + earth_vel[0] * r / C_LIGHT_AU_DAY,
                geo[1] + earth_vel[1] * r / C_LIGHT_AU_DAY,
                geo[2] + earth_vel[2] * r / C_LIGHT_AU_DAY,
            ];
            r = cartesian_norm(geo);
        }
        let spherical = to_spherical(geo);
        if of_date || public_of_date {
            let (lon, lat) = to_of_date(spherical[0], spherical[1], jd)?;
            Ok([lon, lat, r])
        } else {
            Ok(spherical)
        }
    };
    // Vulcan's short orbit moves degrees per day geocentrically, so its
    // central difference uses a tighter step than the slow outer bodies.
    let dt_v = if body == VULCAN { 0.01 } else { 0.1 };
    let center = geo_of_date(tt)?;
    if let Some(complete) = central_complete {
        complete();
    }
    let prev = geo_of_date(tt - dt_v)?;
    let next = geo_of_date(tt + dt_v)?;
    let mut result = if public_of_date {
        model_stencil(center, prev, next, dt_v, false)
    } else {
        [
            center[0],
            center[1],
            center[2],
            wrap_rate(next[0], prev[0], dt_v),
            (next[1] - prev[1]) / (2.0 * dt_v),
            (next[2] - prev[2]) / (2.0 * dt_v),
        ]
    };
    if let Some(rates_at) = spherical_rates {
        if let Some(rates) = rates_at(&geo_of_date, tt, dt_v)? {
            result[3..].copy_from_slice(&rates);
        }
    }
    Ok(result)
}

#[cfg(test)]
mod mutable_model_tests {
    use super::*;
    use std::cell::Cell;

    #[test]
    fn public_uranian_stencil_consumes_only_three_raw_callbacks() {
        let calls = std::cell::RefCell::new(Vec::new());
        let model = |family, tt: f64, source_time: bool| {
            assert_eq!(family, AnalyticalModel::Uranian(CUPIDO));
            assert!(!source_time);
            let index = calls.borrow().len();
            calls.borrow_mut().push(tt);
            Ok(Some(crate::hypothetical::HypoState {
                lon_deg: if index == 0 { 1.0 } else { 2.0 },
                lat_deg: 0.0,
                dist_au: 1.0,
                dlon: 0.0,
                dlat: 0.0,
                ddist: 0.0,
            }))
        };
        let support = |_, _| panic!("heliocentric raw samples need no support");
        let earth = |_| panic!("heliocentric raw samples need no Earth");
        let state = uranian_center_state_with_models(
            Some(&model),
            &support,
            &earth,
            2451545.0,
            CUPIDO,
            FLG_HELCTR | FLG_J2000 | FLG_SPEED,
            false,
        )
        .unwrap();
        assert_eq!(*calls.borrow(), vec![2451545.0, 2451544.0, 2451546.0]);
        assert_eq!(state[0], 1.0);
        assert_eq!(state[3], 0.0);
    }

    #[test]
    fn public_uranian_raw_reselection_failure_stops_before_next_read() {
        let reads = Cell::new(0);
        let model = |_, _, _| {
            reads.set(reads.get() + 1);
            if reads.get() == 2 {
                return Err(Error::UnknownBody { body_id: 777 });
            }
            Ok(Some(crate::hypothetical::HypoState {
                lon_deg: 1.0,
                lat_deg: 0.0,
                dist_au: 1.0,
                dlon: 0.0,
                dlat: 0.0,
                ddist: 0.0,
            }))
        };
        let support = |_, _| panic!("heliocentric raw samples need no support");
        let earth = |_| panic!("heliocentric raw samples need no Earth");
        assert_eq!(
            uranian_center_state_with_models(
                Some(&model),
                &support,
                &earth,
                2451545.0,
                CUPIDO,
                FLG_HELCTR | FLG_J2000 | FLG_SPEED,
                false
            ),
            Err(Error::UnknownBody { body_id: 777 })
        );
        assert_eq!(reads.get(), 2);
    }

    #[test]
    fn uranian_model_failure_precedes_earth_support() {
        let support = Cell::new(0);
        let state_at = |_, _| {
            support.set(support.get() + 1);
            Ok(([0.0; 3], [0.0; 3]))
        };
        let earth = |_| {
            support.set(support.get() + 1);
            Ok([0.0; 3])
        };
        let model = |_, _, _| Err(Error::UnknownBody { body_id: 777 });
        assert_eq!(
            uranian_center_state_with_models(
                Some(&model),
                &state_at,
                &earth,
                2451545.0,
                CUPIDO,
                FLG_J2000,
                false
            ),
            Err(Error::UnknownBody { body_id: 777 })
        );
        assert_eq!(support.get(), 0);
    }

    #[test]
    fn transpluto_geocentric_support_precedes_mutable_model() {
        let models = Cell::new(0);
        let model = |_, _, _| {
            models.set(models.get() + 1);
            Ok(None)
        };
        let state_at = |_, _| unreachable!("failed Earth observation is terminal");
        let earth = |_| Err(Error::UnknownBody { body_id: EARTH });
        assert_eq!(
            isis_center_state_with_models(Some(&model), &state_at, &earth, 2451545.0, FLG_J2000),
            Err(Error::UnknownBody { body_id: EARTH })
        );
        assert_eq!(models.get(), 0);
    }

    #[test]
    fn vulcan_geocentric_support_precedes_mutable_model() {
        let models = Cell::new(0);
        let model = |_, _, _| {
            models.set(models.get() + 1);
            Ok(None)
        };
        let state_at = |_, _| unreachable!("failed Earth observation is terminal");
        let earth = |_| Err(Error::UnknownBody { body_id: EARTH });
        assert_eq!(
            predicted_center_state_source_time_with_models(
                Some(&model),
                &state_at,
                &earth,
                false,
                2451545.0,
                VULCAN,
                0
            ),
            Err(Error::UnknownBody { body_id: EARTH })
        );
        assert_eq!(models.get(), 0);
    }
}

#[cfg(test)]
mod deflection_origin_tests {
    use super::*;
    use crate::constants::{JUPITER, NEPTUNE_ADAMS, SATURN};
    use std::cell::RefCell;

    #[test]
    fn predicted_deflection_is_invariant_under_common_support_translation() {
        // Synthetic stationary support. All positions and translations are
        // dyadic, so subtracting the common origin is exact in binary64.
        // Model coordinates are heliocentric and must not follow an arbitrary
        // translation of every barycentric support position.
        let run = |translation: [f64; 3], with_models: bool, public_of_date: bool, flags| {
            let reads = RefCell::new(Vec::new());
            let state_at = |body, epoch: f64| -> Result<([f64; 3], [f64; 3])> {
                reads.borrow_mut().push((body, epoch.to_bits()));
                let position = match body {
                    SUN => [0.125, -0.25, 0.5],
                    EARTH => [1.125, -0.75, 0.75],
                    JUPITER => [4.125, 2.25, -1.5],
                    SATURN => [-3.875, 7.75, 2.5],
                    _ => return Err(Error::UnknownBody { body_id: body }),
                };
                Ok((
                    std::array::from_fn(|axis| position[axis] + translation[axis]),
                    [0.0; 3],
                ))
            };
            // The same Sun-relative Earth observation in either coordinate
            // origin. Zero support velocities make retardation immaterial.
            let earth_at = |_: f64| {
                Ok(RetardedEarthObservation {
                    relative: [1.0, -0.5, 0.25],
                    sun_reception: std::array::from_fn(|axis| {
                        [0.125, -0.25, 0.5][axis] + translation[axis]
                    }),
                })
            };
            let result = if with_models {
                predicted_center_state_source_time_with_models(
                    None,
                    &state_at,
                    &earth_at,
                    public_of_date,
                    2451545.0,
                    NEPTUNE_ADAMS,
                    flags,
                )
            } else {
                predicted_center_state_source_time(
                    &state_at,
                    &earth_at,
                    public_of_date,
                    2451545.0,
                    NEPTUNE_ADAMS,
                    flags,
                )
            }
            .unwrap();
            assert!(result.into_iter().all(f64::is_finite));
            (result.map(f64::to_bits), reads.into_inner())
        };
        for with_models in [false, true] {
            for public_of_date in [false, true] {
                let flags = FLG_J2000 | FLG_SPEED | FLG_NOABERR;
                // The no-deflection control locates any failure specifically
                // in the PPN path, rather than the model or spherical stencil.
                let omitted = flags | FLG_NOGDEFL;
                assert_eq!(
                    run([0.0; 3], with_models, public_of_date, omitted),
                    run([16.0, -8.0, 4.0], with_models, public_of_date, omitted),
                );
                let original = run([0.0; 3], with_models, public_of_date, flags);
                let translated = run([16.0, -8.0, 4.0], with_models, public_of_date, flags);
                assert_eq!(
                    original.0, translated.0,
                    "common-origin invariant: models={with_models}, public={public_of_date}"
                );
                assert_eq!(
                    original.1, translated.1,
                    "closest-approach epochs must also be origin independent"
                );
                assert_eq!(original.1.len(), 21);
                for sample in original.1.chunks_exact(7) {
                    assert_eq!(
                        sample.iter().map(|read| read.0).collect::<Vec<_>>(),
                        [EARTH, SUN, SUN, JUPITER, JUPITER, SATURN, SATURN]
                    );
                }
            }
        }
    }

    #[test]
    fn scalar_observation_retains_reception_sun_and_retarded_earth() {
        let tt = 2451545.0;
        let reads = RefCell::new(Vec::new());
        let state = |body, epoch| {
            reads.borrow_mut().push((body, epoch));
            let position = match body {
                SUN => {
                    assert_eq!(epoch, tt);
                    [8.0, -4.0, 2.0]
                }
                EARTH if epoch == tt => [9.0, -4.0, 2.0],
                EARTH => [10.0, -4.0, 2.0],
                _ => panic!("unexpected acquisition"),
            };
            Ok((position, [0.0; 3]))
        };
        let observation = observed_earth_with_origin(&state, tt).unwrap();
        assert_eq!(observation.sun_reception, [8.0, -4.0, 2.0]);
        assert_eq!(observation.relative, [2.0, 0.0, 0.0]);
        let reads = reads.into_inner();
        assert_eq!(
            reads.iter().map(|read| read.0).collect::<Vec<_>>(),
            [SUN, EARTH, EARTH, EARTH]
        );
        assert_eq!(reads[0].1, tt);
        assert_eq!(reads[1].1, tt);
        assert!(reads[2].1 < tt && reads[3].1 < reads[2].1);
    }

    struct EventSupport {
        translation: [f64; 3],
        reads: RefCell<Vec<(i32, f64)>>,
    }
    impl StateProvider for EventSupport {
        fn provider_name(&self) -> &'static str {
            "synthetic-deflection-support"
        }
        fn has_body(&self, body: i32) -> bool {
            matches!(body, SUN | EARTH | JUPITER | SATURN)
        }
        fn body_coverage(&self, _: i32) -> Option<(f64, f64)> {
            Some((2400000.0, 2500000.0))
        }
        fn serving_coord_type(&self, _: i32) -> Option<u32> {
            Some(COORD_ICRS_BARY)
        }
        fn eval_raw(&self, body: i32, tt: f64) -> Result<(Vec<f64>, Vec<f64>)> {
            self.reads.borrow_mut().push((body, tt));
            let position = match body {
                SUN => [0.125, -0.25, 0.5],
                EARTH => [1.125, -0.75, 0.75],
                JUPITER => [4.125, 2.25, -1.5],
                SATURN => [-3.875, 7.75, 2.5],
                _ => return Err(Error::LebBodyNotFound { body_id: body }),
            };
            Ok((
                (0..3)
                    .map(|axis| position[axis] + self.translation[axis])
                    .collect(),
                vec![0.0; 3],
            ))
        }
        fn eval_raw_two_part(&self, _: i32, _: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)> {
            panic!("event model uses scalar TT")
        }
        fn serving_entry_snapshot(&self, _: i32) -> Result<EntrySnapshot> {
            Ok(EntrySnapshot {
                coord_type: COORD_ICRS_BARY,
            })
        }
    }

    #[test]
    fn event_predicted_deflection_preserves_common_origin_and_acquisition_order() {
        let run = |translation| {
            let provider = EventSupport {
                translation,
                reads: RefCell::new(Vec::new()),
            };
            let nutation = |_| Ok((0.0, 0.0));
            let serving = SealedServing::new(&provider, &nutation);
            let delta_t = |_| Ok(0.0);
            let event = crate::events::SealedFictitiousEventServing::try_new(
                &serving,
                &delta_t,
                NEPTUNE_ADAMS,
                FLG_SWIEPH,
            )
            .unwrap();
            let result = event.fict_geo_ecliptic_at_tt(2451545.0).unwrap();
            let words = [result.0, result.1, result.2].map(f64::to_bits);
            (words, provider.reads.into_inner())
        };
        let plain = run([0.0; 3]);
        let translated = run([16.0, -8.0, 4.0]);
        assert_eq!(plain, translated);
        assert_eq!(
            plain.1.iter().map(|read| read.0).collect::<Vec<_>>(),
            [
                SUN, EARTH, EARTH, EARTH, SUN, SUN, JUPITER, JUPITER, SATURN, SATURN
            ]
        );
        assert_eq!(plain.1[0].1, 2451545.0);
        assert_eq!(plain.1[1].1, 2451545.0);
        assert!(plain.1[2].1 < 2451545.0);
        assert_eq!(plain.1[3].1, 2451545.0); // Velocity at reception, not retarded Earth.
    }

    #[test]
    fn supplied_model_keeps_heliocentric_words_and_acquisition_order() {
        let run = |translation: [f64; 3], public_of_date| {
            let calls = RefCell::new(Vec::new());
            let state = |body, tt: f64| {
                calls.borrow_mut().push((false, body, tt.to_bits()));
                let position = match body {
                    SUN => [0.125, -0.25, 0.5],
                    EARTH => [1.125, -0.75, 0.75],
                    JUPITER => [4.125, 2.25, -1.5],
                    SATURN => [-3.875, 7.75, 2.5],
                    _ => return Err(Error::UnknownBody { body_id: body }),
                };
                Ok((
                    std::array::from_fn(|axis| position[axis] + translation[axis]),
                    [0.0; 3],
                ))
            };
            let model = |family, tt: f64, source_time: bool| {
                assert_eq!(family, AnalyticalModel::Vulcan);
                assert!(!source_time);
                calls.borrow_mut().push((true, VULCAN, tt.to_bits()));
                Ok(Some(crate::hypothetical::HypoState {
                    lon_deg: 45.0,
                    lat_deg: 3.0,
                    dist_au: 2.0,
                    dlon: 0.0,
                    dlat: 0.0,
                    ddist: 0.0,
                }))
            };
            let result = predicted_center_state_source_time_with_models(
                Some(&model),
                &state,
                &|tt| observed_earth_with_origin(&state, tt),
                public_of_date,
                2451545.0,
                VULCAN,
                FLG_J2000 | FLG_SPEED | FLG_NOABERR,
            )
            .unwrap();
            assert!(result.into_iter().all(f64::is_finite));
            (result.map(f64::to_bits), calls.into_inner())
        };
        for public in [false, true] {
            let original = run([0.0; 3], public);
            assert_eq!(original, run([16.0, -8.0, 4.0], public));
            assert_eq!(original.1.len(), 39);
            for sample in original.1.chunks_exact(13) {
                assert_eq!(
                    sample
                        .iter()
                        .map(|read| (read.0, read.1))
                        .collect::<Vec<_>>(),
                    [
                        (false, SUN),
                        (false, EARTH),
                        (false, EARTH),
                        (false, EARTH),
                        (true, VULCAN),
                        (true, VULCAN),
                        (true, VULCAN),
                        (false, SUN),
                        (false, SUN),
                        (false, JUPITER),
                        (false, JUPITER),
                        (false, SATURN),
                        (false, SATURN),
                    ]
                );
            }
        }
    }
}
