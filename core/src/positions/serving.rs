//! Serving implementation for the parent module.
//! Public entry points and types are exposed through the parent API.

use super::*;

use super::routing;

mod rates;

mod flags;
pub use flags::{
    calc_tt_epheflag_echo, consume_non_echoed_flags, degrade_jpl_echo_on_fallback,
    echo_request_bits, exclusive_ephemeris_bit, implied_retflag_bits, normalize_calc_flags,
    strip_output_flags, validate_barycentric_moseph,
};

mod planet_centered;
pub(crate) use planet_centered::pctr_missing_vector_alias;
#[cfg(test)]
use planet_centered::{
    PctrRatePolicy, pctr_tt_parts, pctr_vector_result, pctr_vector_result_acquired,
};
pub use planet_centered::{
    calc_pctr, calc_pctr_context_ut1_parts_in_context, calc_pctr_in_context,
    calc_pctr_in_context_with_anchors, is_pctr_supported_body, pctr_apparent_vector,
    pctr_central_speed, pctr_position_no_speed, pctr_speed_half_step_days,
};

/// Rebuild an apparent ICRS AU ray at a local stencil epoch.
///
/// Linearly advance the geometric ray, replace sampled site curvature, and
/// recompute light time before enabled deflection and aberration. The Sun's
/// central deflection increment stays fixed; other deflectors are reacquired
/// through `deflect` at `jd_sample` (TT JD). Observer velocities are AU/day.
/// Non-finite inputs and correction/provider errors propagate.
pub fn sample_apparent_vector(
    center: &StencilCenter,
    sample_offset_days: f64,
    earth_vel_sample: [f64; 3],
    jd_sample: f64,
    deflect: impl Fn([f64; 3], [f64; 3], f64, f64) -> Result<[f64; 3]>,
) -> Result<[f64; 3]> {
    if !sample_offset_days.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: sample_offset_days,
        });
    }
    if !jd_sample.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_sample });
    }
    check_finite_vec3(earth_vel_sample)?;
    let mut g = extrapolate_geometric_vector(center.geo_geom, center.geo_vel, sample_offset_days);
    if let Some(topo) = center.topo {
        if let Some((sample_pos, _)) = topo.sample {
            g = remove_topo_curvature(
                g,
                topo.central_pos,
                topo.central_vel,
                sample_pos,
                sample_offset_days,
            );
        }
    }
    let lt_s = sample_light_time_days(g);
    if center.do_defl {
        if center.target_body_id == SUN {
            // The Sun is its own dominant deflector, so the PPN geometry
            // is degenerate on the extrapolated vector (a huge spurious
            // rate); the term is physically constant over ±_VEL_H, so the
            // position-path delta is frozen instead of re-evaluated.
            g = [
                g[0] + center.defl_delta[0],
                g[1] + center.defl_delta[1],
                g[2] + center.defl_delta[2],
            ];
        } else {
            let obs_s = extrapolate_geometric_vector(
                center.observer,
                center.observer_vel,
                sample_offset_days,
            );
            g = deflect(g, obs_s, jd_sample, lt_s)?;
        }
    }
    if center.do_aber {
        // Annual and diurnal aberration use the complete sampled observer velocity.
        // Without a supplied site sample, retain its central nominal velocity.
        let ev_s = match center.topo {
            Some(topo) => {
                let tv = topo.sample.map(|(_, vel)| vel).unwrap_or(topo.central_vel);
                [
                    earth_vel_sample[0] + tv[0],
                    earth_vel_sample[1] + tv[1],
                    earth_vel_sample[2] + tv[2],
                ]
            }
            None => earth_vel_sample,
        };
        g = apply_aberration(g, ev_s, lt_s)?;
    }
    Ok(g)
}

/// Differentiate two output triples at the fixed half-step [`VEL_H_DAYS`].
///
/// Spherical triples are degrees/degrees/AU; longitude differences use the
/// short signed arc across a branch cut. Cartesian triples differentiate each
/// axis in AU. Output rates are per TT day. Non-finite components are errors.
pub fn difference_stencil_rate(
    w_minus: (f64, f64, f64),
    w_plus: (f64, f64, f64),
    want_xyz: bool,
) -> Result<(f64, f64, f64)> {
    for v in [
        w_minus.0, w_minus.1, w_minus.2, w_plus.0, w_plus.1, w_plus.2,
    ] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    let inv2h = 1.0 / (2.0 * VEL_H_DAYS);
    if want_xyz {
        return Ok((
            (w_plus.0 - w_minus.0) * inv2h,
            (w_plus.1 - w_minus.1) * inv2h,
            (w_plus.2 - w_minus.2) * inv2h,
        ));
    }
    let mut dlon = (w_plus.0 - w_minus.0).rem_euclid(360.0);
    if dlon > 180.0 {
        dlon -= 360.0;
    }
    Ok((
        dlon * inv2h,
        (w_plus.1 - w_minus.1) * inv2h,
        (w_plus.2 - w_minus.2) * inv2h,
    ))
}

/// Central-difference an apparent place and its output frame at
/// `jd_tt ± VEL_H_DAYS`.
///
/// `frame_at` projects each sampled ICRS ray and `deflect` supplies its enabled
/// PPN correction. Observer velocities are AU/day. Returns angular rates in
/// degrees/TT-day and radial rate in AU/TT-day, or Cartesian AU/TT-day when
/// `want_xyz`. Epoch, provider, frame and correction errors propagate.
pub fn stencil_velocity(
    center: &StencilCenter,
    jd_tt: f64,
    earth_vel_minus: [f64; 3],
    earth_vel_plus: [f64; 3],
    want_xyz: bool,
    deflect: impl Fn([f64; 3], [f64; 3], f64, f64) -> Result<[f64; 3]>,
    frame_at: impl Fn([f64; 3], f64) -> Result<(f64, f64, f64)>,
) -> Result<(f64, f64, f64)> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let jd_m = jd_tt - VEL_H_DAYS;
    let jd_p = jd_tt + VEL_H_DAYS;
    let g_m = sample_apparent_vector(center, -VEL_H_DAYS, earth_vel_minus, jd_m, &deflect)?;
    // Borrow the same deflector callback for both stencil samples.
    let g_p = sample_apparent_vector(center, VEL_H_DAYS, earth_vel_plus, jd_p, &deflect)?;
    let w_m = frame_at(g_m, jd_m)?;
    let w_p = frame_at(g_p, jd_p)?;
    difference_stencil_rate(w_m, w_p, want_xyz)
}

/// Identity of one requested native motion within the same provider/feed borrow.
#[derive(Clone, Copy, PartialEq, Eq)]
struct OutputMotionQuery {
    tt: u64,
    split: Option<[u64; 2]>,
    body: i32,
    flags: i32,
    xyz: bool,
}
impl OutputMotionQuery {
    fn new(tt: f64, split: Option<(f64, f64)>, body: i32, flags: i32, xyz: bool) -> Self {
        Self {
            tt: tt.to_bits(),
            split: split.map(|(a, b)| [a.to_bits(), b.to_bits()]),
            body,
            flags,
            xyz,
        }
    }
}
#[derive(Clone, Copy)]
struct OutputMotionRequest {
    query: OutputMotionQuery,
    motion: Option<super::derivative::OutputWordMotion>,
}
/// Local synchronous carrier; no motion survives the requesting invocation.
#[derive(Default)]
pub(super) struct OutputMotionCapture {
    request: std::cell::Cell<Option<OutputMotionRequest>>,
}
struct RestoreOutputMotion<'a> {
    slot: &'a std::cell::Cell<Option<OutputMotionRequest>>,
    previous: Option<OutputMotionRequest>,
}
impl Drop for RestoreOutputMotion<'_> {
    fn drop(&mut self) {
        self.slot.set(self.previous);
    }
}
impl OutputMotionCapture {
    fn with_motion<R>(
        &self,
        query: OutputMotionQuery,
        calculate: impl FnOnce() -> Result<R>,
    ) -> Result<(R, Option<super::derivative::OutputWordMotion>)> {
        let previous = self.request.replace(Some(OutputMotionRequest {
            query,
            motion: None,
        }));
        let _restore = RestoreOutputMotion {
            slot: &self.request,
            previous,
        };
        let result = calculate()?;
        let motion = self
            .request
            .get()
            .filter(|r| r.query == query)
            .and_then(|r| r.motion);
        Ok((result, motion))
    }
    fn publish(&self, query: OutputMotionQuery, motion: super::derivative::OutputWordMotion) {
        if let Some(mut request) = self.request.get() {
            if request.query == query {
                request.motion = Some(motion);
                self.request.set(Some(request));
            }
        }
    }
}

/// A rotating surface observer's position and velocity in AU and AU/TT-day.
type SiteStateAt<'a> = dyn Fn(f64) -> Result<([f64; 3], [f64; 3])> + 'a;
type SiteClock<'a> = (&'a dyn Fn(f64) -> Result<f64>, (f64, f64));

/// Actual output place and its retained small physical correction, in AU.
#[derive(Clone, Copy, Default)]
struct AngularModelSample {
    output: [f64; 3],
    apparent: [f64; 3],
    correction: [f64; 3],
    reads: [Option<CorrectionRead>; DEFLECTORS.len() * 2],
}

/// One successful PPN lookup at its actual closest-approach or observing epoch.
#[derive(Clone, Copy)]
struct CorrectionRead {
    body: i32,
    epoch: f64,
    interval: Option<(f64, f64)>,
}

/// One successful scalar acquisition, with its actual input words retained.
#[derive(Clone, Copy)]
struct PpnScalarSample {
    body: i32,
    tt: u64,
    state: GeometricState,
}

/// Fixed scratch for one position's scalar PPN reads, borrowed by its SPEED tail.
/// No source owner, failed lookup, neighboring place or result survives the call.
#[derive(Default)]
pub(super) struct PpnScalarSamples {
    slots: [std::cell::Cell<Option<PpnScalarSample>>; DEFLECTORS.len() * 2],
}

impl PpnScalarSamples {
    fn record(&self, index: usize, body: i32, tt: f64, state: GeometricState) {
        // The existing CorrectionRead admission already checked this same slot.
        self.slots[index].set(Some(PpnScalarSample {
            body,
            tt: tt.to_bits(),
            state,
        }));
    }

    pub(super) fn state_at(&self, body: i32, tt: f64) -> Option<GeometricState> {
        self.slots
            .iter()
            .filter_map(std::cell::Cell::get)
            .find_map(|sample| {
                (sample.body == body
                    && sample.tt == tt.to_bits()
                    && sample.state.body_id() == body
                    && sample.state.jd_tt().to_bits() == tt.to_bits())
                .then_some(sample.state)
            })
    }
}

impl AngularModelSample {
    fn correction_support(self, mut support: (f64, f64)) -> (f64, f64) {
        for read in self.reads.into_iter().flatten() {
            if let Some((start, end)) = read.interval {
                support.0 = support.0.max(start - read.epoch);
                support.1 = support.1.min(end - read.epoch);
            }
        }
        support
    }

    fn owns_correction_sample(self, sample: Self) -> bool {
        self.reads
            .iter()
            .zip(sample.reads)
            .all(|(central, sample)| match (central, sample) {
                (None, None) => true,
                (Some(central), Some(sample)) => {
                    central.body == sample.body
                        && central.interval == sample.interval
                        && central
                            .interval
                            .is_none_or(|(start, end)| start <= sample.epoch && sample.epoch <= end)
                }
                _ => false,
            })
    }
}

/// Complete query identity for a central sample inside one borrowed serving.
/// The provider and stored feed are fixed by that borrow, never by a name or
/// callback address. Keep original TT parts and every non-speed flag exact.
#[derive(Clone, Copy, PartialEq, Eq)]
struct CentralAngularQuery {
    tt: u64,
    parts: Option<[u64; 2]>,
    body: i32,
    flags: i32,
    order: CentralStateOrder,
}

impl CentralAngularQuery {
    fn new(
        tt: f64,
        parts: Option<(f64, f64)>,
        body: i32,
        flags: i32,
        order: CentralStateOrder,
    ) -> Self {
        Self {
            tt: tt.to_bits(),
            parts: parts.map(|(hi, lo)| [hi.to_bits(), lo.to_bits()]),
            body,
            flags: flags & !(FLG_SPEED | FLG_SPEED3),
            order,
        }
    }
}

/// Successful central geometry, borrowed only until its SPEED tail finishes.
/// Raw correction reads defer metadata validation to the original sample point.
struct CentralAngularPlace<'a> {
    query: CentralAngularQuery,
    apparent: [f64; 3],
    correction: [f64; 3],
    frames: &'a FrameMatrices,
    reads: [Option<CorrectionRead>; DEFLECTORS.len() * 2],
}

/// Reuse only identical inputs from a producer-owned immutable serving.
/// Projection, correction-domain failures and tracing retain their place in
/// the sampler. Arbitrary providers/feeds always execute the original route.
fn central_angular_sample<P: StateProvider>(
    serving: &SealedServing<'_, P>,
    central: &CentralAngularPlace<'_>,
    query: CentralAngularQuery,
) -> Result<Option<AngularModelSample>> {
    if !serving.reuse_immutable_central || central.query != query {
        return Ok(None);
    }
    let output = frame_transform(
        central.apparent,
        f64::from_bits(query.tt),
        query.flags,
        true,
        *central.frames,
    )?;
    let reads = correction_read_domains(serving.provider, central.reads)?;
    super::routing::record_stored_trace_source(query.body)?;
    Ok(Some(AngularModelSample {
        output,
        apparent: central.apparent,
        correction: central.correction,
        reads,
    }))
}

/// Record successful correction lookups without changing optional-missing policy.
fn capture_correction_read<P: StateProvider>(
    provider: &P,
    body: i32,
    epoch: f64,
    capture: Option<&std::cell::Cell<[Option<CorrectionRead>; DEFLECTORS.len() * 2]>>,
    scalar_capture: Option<&PpnScalarSamples>,
) -> Result<[f64; 3]> {
    let state = provider.geometric_state(body, epoch)?;
    if let Some(capture) = capture {
        let mut reads = capture.get();
        let next = reads
            .iter()
            .position(Option::is_none)
            .ok_or(Error::MathDomain { value: epoch })?;
        reads[next] = Some(CorrectionRead {
            body,
            epoch,
            interval: None,
        });
        capture.set(reads);
        if let Some(scalar_capture) = scalar_capture {
            scalar_capture.record(next, body, epoch, state);
        }
    }
    Ok(state.pos())
}

/// Resolve metadata after correction policy has accepted the successful reads.
/// Corrupt domain metadata is not reclassified as an optional missing body.
fn correction_read_domains<P: StateProvider>(
    provider: &P,
    mut reads: [Option<CorrectionRead>; DEFLECTORS.len() * 2],
) -> Result<[Option<CorrectionRead>; DEFLECTORS.len() * 2]> {
    for read in reads.iter_mut().flatten() {
        read.interval = provider.smooth_interval_tt(read.body, read.epoch)?;
    }
    Ok(reads)
}

/// Differentiate the actual Cartesian place in its central angular tangent plane.
///
/// The fixed orthonormal longitude/latitude directions give the exact spherical
/// Jacobian at the central ray. Projecting Cartesian displacements before the
/// derivative avoids subtracting rounded angles on short owned intervals.
/// Samples remain full physical places, in the same output frame and AU units;
/// no trajectory or polynomial is extrapolated. Rates are degrees per TT day.
#[cfg(test)]
fn cartesian_angular_velocity(
    place_at: &dyn Fn(f64) -> Result<[f64; 3]>,
    support: (f64, f64),
) -> Result<[f64; 2]> {
    super::derivative::cartesian_angular_rates_on_owned_domain(
        &|offset| place_at(offset).map(Some),
        1.0 / 64.0,
        support,
    )
}

/// Differentiate the same owned local geometry in every serving branch.
/// Angular rates use exact dyadic fourth-order samples. The established
/// radial/Cartesian stencil remains independent of this angular improvement.
#[allow(clippy::too_many_arguments)]
fn local_model_velocity<P: StateProvider>(
    serving: &SealedServing<'_, P>,
    center: &StencilCenter,
    geometric_velocity: [f64; 3],
    jd_tt: f64,
    flags: i32,
    want_xyz: bool,
    required_lookup: bool,
    topo_at: Option<&SiteStateAt<'_>>,
    angular_at: &dyn Fn(f64) -> Result<AngularModelSample>,
    geometric_support: (f64, f64),
) -> Result<(f64, f64, f64)> {
    let sample = |offset: f64| -> Result<[f64; 3]> {
        let epoch = jd_tt + offset;
        let mut local = *center;
        if let (Some(topo), Some(at)) = (&mut local.topo, topo_at) {
            topo.sample = Some(at(epoch)?);
        }
        let earth_velocity = if center.do_aber {
            serving.provider.geometric_state(EARTH, epoch)?.vel()
        } else {
            [0.0; 3]
        };
        let vector = sample_apparent_vector(
            &local,
            offset,
            earth_velocity,
            epoch,
            |g, observer, time, light_time| {
                if required_lookup {
                    apply_gravitational_deflection_with_lookup(
                        g,
                        observer,
                        time,
                        light_time,
                        |id, jd| {
                            serving
                                .provider
                                .geometric_state(id, jd)
                                .map(|state| Some(state.pos()))
                        },
                    )
                } else {
                    apply_gravitational_deflection(g, observer, time, light_time, |id, jd| {
                        serving
                            .provider
                            .geometric_state(id, jd)
                            .map(|state| state.pos())
                    })
                }
            },
        )?;
        let (dpsi, deps) = (serving.nutation_at)(epoch)?;
        frame_transform(
            vector,
            epoch,
            flags,
            want_xyz,
            frame_matrices_for_date(epoch, dpsi, deps)?,
        )
    };
    let mut support = geometric_support;
    let mut include = |span: Option<(f64, f64)>| {
        if let Some((start, end)) = span {
            support.0 = support.0.max(start - jd_tt);
            support.1 = support.1.min(end - jd_tt);
        }
    };
    include(serving.provider.nutation_smooth_interval_tt(jd_tt)?);
    if center.do_aber {
        include(serving.provider.body_coverage(EARTH));
        include(serving.provider.smooth_interval_tt(EARTH, jd_tt)?);
    }
    if center.do_defl && center.target_body_id != SUN {
        for deflector in [SUN, JUPITER, SATURN] {
            if serving.provider.has_body(deflector) {
                include(serving.provider.body_coverage(deflector));
                include(serving.provider.smooth_interval_tt(deflector, jd_tt)?);
            }
        }
    }
    if want_xyz {
        let velocity = super::derivative::cartesian_rates(
            &|offset| angular_at(offset).map(|sample| sample.output),
            1.0 / 64.0,
            support,
        )?;
        return Ok((velocity[0], velocity[1], velocity[2]));
    }
    let central = angular_at(0.0)?;
    let correction_support = central.correction_support(support);
    let correction_rate = super::derivative::cartesian_rates_on_owned_domain(
        &|offset| {
            let sample = angular_at(offset)?;
            Ok(central
                .owns_correction_sample(sample)
                .then_some(sample.correction))
        },
        1.0 / 64.0,
        correction_support,
    )?;
    let apparent_velocity = std::array::from_fn(|i| geometric_velocity[i] + correction_rate[i]);
    let angular = apparent_output_rates(
        serving,
        center.target_body_id,
        None,
        jd_tt,
        flags,
        central.apparent,
        central.output,
        apparent_velocity,
        false,
        false,
        None,
    )?;
    let minus = sample(-VEL_H_DAYS)?;
    let plus = sample(VEL_H_DAYS)?;
    Ok((
        angular.0,
        angular.1,
        (plus[2] - minus[2]) / (2.0 * VEL_H_DAYS),
    ))
}

/// Apply the frame derivative and the Cartesian/spherical output Jacobian.
/// Position arithmetic stays with the original scalar place. Rates are per TT
/// day, in AU or degrees; the range rate differentiates the returned norm.
#[allow(clippy::too_many_arguments)]
fn apparent_output_rates<P: StateProvider>(
    serving: &SealedServing<'_, P>,
    body_id: i32,
    obs_split: Option<(f64, f64)>,
    jd_tt: f64,
    flags: i32,
    apparent: [f64; 3],
    output: [f64; 3],
    apparent_velocity: [f64; 3],
    want_xyz: bool,
    analytic_frame: bool,
    expanded: Option<(
        [super::derivative::AngularWord; 3],
        [super::derivative::AngularWord; 3],
    )>,
) -> Result<(f64, f64, f64)> {
    let nutation = (serving.nutation_at)(jd_tt)?;
    let stored_rate = if analytic_frame && serving.stored_nutation_feed && flags & FLG_J2000 == 0 {
        serving.provider.nutation_rate_tt(jd_tt)?
    } else {
        None
    };
    let nutation_rate = if flags & FLG_J2000 != 0 {
        [0.0; 3]
    } else if let Some((psi, deps)) = stored_rate {
        check_finite_vec3([psi, deps, 0.0])?;
        [psi, deps, 0.0]
    } else {
        let nutation_support = serving
            .provider
            .nutation_smooth_interval_tt(jd_tt)?
            .map_or((f64::NEG_INFINITY, f64::INFINITY), |(start, end)| {
                (start - jd_tt, end - jd_tt)
            });
        super::derivative::cartesian_rates(
            &|offset| (serving.nutation_at)(jd_tt + offset).map(|(psi, deps)| [psi, deps, 0.0]),
            1.0 / 64.0,
            nutation_support,
        )?
    };
    if let Some((position, velocity)) = expanded {
        let (matrix, matrix_rate) =
            super::frame_velocity::analytic_output_frame_words_with_nutation(
                jd_tt,
                flags,
                nutation,
                (nutation_rate[0], nutation_rate[1]),
            )?;
        let motion =
            super::derivative::expanded_output_motion(matrix, matrix_rate, position, velocity)?;
        let rates = super::derivative::expanded_motion_rates(motion, want_xyz)?;
        // A failed enclosing call never consumes this local value: its scope
        // restores the previous slot on Result errors, unwind and nested calls.
        serving.output_motion.publish(
            OutputMotionQuery::new(jd_tt, obs_split, body_id, flags, want_xyz),
            motion,
        );
        return Ok(rates);
    }
    let frame = if analytic_frame {
        super::frame_velocity::analytic_output_frame_matrix_rate_with_nutation
    } else {
        super::frame_velocity::output_frame_matrix_rate_with_nutation
    };
    let (matrix, matrix_rate) =
        frame(jd_tt, flags, nutation, (nutation_rate[0], nutation_rate[1]))?;
    let output_velocity = {
        let motion = mat3_vec3(matrix, apparent_velocity)?;
        let rotation = mat3_vec3(matrix_rate, apparent)?;
        std::array::from_fn(|i| motion[i] + rotation[i])
    };
    if want_xyz {
        return Ok((output_velocity[0], output_velocity[1], output_velocity[2]));
    }
    let angular = if analytic_frame {
        super::derivative::spherical_angular_velocity_direct(output, output_velocity)?
    } else {
        super::derivative::spherical_angular_velocity(output, output_velocity)?
    };
    Ok((
        angular[0],
        angular[1],
        super::planet_range::RangeState::from_cartesian(output, output_velocity).rate,
    ))
}

#[cfg(test)]
mod physical_angular_model_tests {
    use super::*;

    struct InertialModel {
        origin: [f64; 3],
        position: [f64; 3],
        velocity: [f64; 3],
    }

    impl StateProvider for InertialModel {
        fn provider_name(&self) -> &'static str {
            "inertial-test-model"
        }
        fn has_body(&self, body: i32) -> bool {
            matches!(body, SUN | EARTH | MERCURY)
        }
        fn body_coverage(&self, body: i32) -> Option<(f64, f64)> {
            self.has_body(body)
                .then_some((FRAME_J2000 - 10.0, FRAME_J2000 + 10.0))
        }
        fn serving_coord_type(&self, body: i32) -> Option<u32> {
            self.has_body(body).then_some(COORD_ICRS_BARY)
        }
        fn serving_entry_snapshot(&self, body: i32) -> Result<EntrySnapshot> {
            if self.has_body(body) {
                Ok(EntrySnapshot {
                    coord_type: COORD_ICRS_BARY,
                })
            } else {
                Err(Error::LebBodyNotFound { body_id: body })
            }
        }
        fn eval_raw(&self, body: i32, epoch: f64) -> Result<(Vec<f64>, Vec<f64>)> {
            match body {
                SUN | EARTH => Ok((self.origin.to_vec(), vec![0.0; 3])),
                MERCURY => Ok((
                    self.position
                        .iter()
                        .zip(self.velocity)
                        .enumerate()
                        .map(|(i, (p, v))| self.origin[i] + p + v * (epoch - FRAME_J2000))
                        .collect(),
                    self.velocity.to_vec(),
                )),
                _ => Err(Error::LebBodyNotFound { body_id: body }),
            }
        }
        fn eval_relative_position_two_part(
            &self,
            target: i32,
            epoch: TwoPartEpoch,
            center: i32,
            _center_epoch: TwoPartEpoch,
        ) -> Result<Option<[f64; 3]>> {
            if target == MERCURY && matches!(center, SUN | EARTH) {
                Ok(Some(std::array::from_fn(|i| {
                    self.position[i] + self.velocity[i] * ((epoch.hi - FRAME_J2000) + epoch.lo)
                })))
            } else {
                Ok(None)
            }
        }
        fn eval_raw_two_part(
            &self,
            body: i32,
            epoch: TwoPartEpoch,
        ) -> Result<(Vec<f64>, Vec<f64>)> {
            match body {
                SUN | EARTH => Ok((self.origin.to_vec(), vec![0.0; 3])),
                MERCURY => Ok((
                    self.position
                        .iter()
                        .zip(self.velocity)
                        .enumerate()
                        .map(|(i, (p, v))| {
                            self.origin[i] + p + v * ((epoch.hi - FRAME_J2000) + epoch.lo)
                        })
                        .collect(),
                    self.velocity.to_vec(),
                )),
                _ => Err(Error::LebBodyNotFound { body_id: body }),
            }
        }
    }

    fn angular_jacobian(ray: [f64; 3], rate: [f64; 3], flags: i32) -> [f64; 2] {
        let frames = frame_matrices_for_date(FRAME_J2000, 0.0, 0.0).unwrap();
        let p = frame_transform(ray, FRAME_J2000, flags, true, frames).unwrap();
        let v = frame_transform(rate, FRAME_J2000, flags, true, frames).unwrap();
        let horizontal = p[0].hypot(p[1]);
        let radius2 = horizontal * horizontal + p[2] * p[2];
        let horizontal_rate = (p[0] * v[0] + p[1] * v[1]) / horizontal;
        [
            ((p[0] * v[1] - p[1] * v[0]) / (horizontal * horizontal)).to_degrees(),
            ((horizontal * v[2] - p[2] * horizontal_rate) / radius2).to_degrees(),
        ]
    }

    /// Counters observe calls without changing the affine physical model.
    struct CountedModel {
        model: InertialModel,
        central_reads: std::cell::Cell<usize>,
        domains: std::cell::Cell<usize>,
        fail_domains: bool,
        reuse_states: bool,
        split_reads: std::cell::Cell<usize>,
        proof_calls: std::cell::Cell<usize>,
        fail_split_read: Option<usize>,
        split_queries: std::cell::RefCell<Vec<[u64; 3]>>,
        relative_queries: std::cell::RefCell<Vec<[u64; 6]>>,
    }

    impl CountedModel {
        fn new() -> Self {
            Self {
                model: InertialModel {
                    origin: [1.0, -0.8, 0.2],
                    position: [2.0, 0.2, 0.3],
                    velocity: [0.005, 0.002, -0.001],
                },
                central_reads: std::cell::Cell::new(0),
                domains: std::cell::Cell::new(0),
                fail_domains: false,
                reuse_states: false,
                split_reads: std::cell::Cell::new(0),
                proof_calls: std::cell::Cell::new(0),
                fail_split_read: None,
                split_queries: std::cell::RefCell::new(Vec::new()),
                relative_queries: std::cell::RefCell::new(Vec::new()),
            }
        }

        fn reset_queries(&self) {
            self.central_reads.set(0);
            self.split_reads.set(0);
            self.proof_calls.set(0);
            self.split_queries.borrow_mut().clear();
            self.relative_queries.borrow_mut().clear();
        }
    }

    impl StateProvider for CountedModel {
        fn provider_name(&self) -> &'static str {
            "counted-affine-test-model"
        }
        fn has_body(&self, body: i32) -> bool {
            self.model.has_body(body)
        }
        fn body_coverage(&self, body: i32) -> Option<(f64, f64)> {
            self.model.body_coverage(body)
        }
        fn serving_coord_type(&self, body: i32) -> Option<u32> {
            self.model.serving_coord_type(body)
        }
        fn serving_entry_snapshot(&self, body: i32) -> Result<EntrySnapshot> {
            self.model.serving_entry_snapshot(body)
        }
        fn smooth_interval_tt(&self, _body: i32, epoch: f64) -> Result<Option<(f64, f64)>> {
            self.domains.set(self.domains.get() + 1);
            if self.fail_domains {
                return Err(Error::MathDomain { value: epoch });
            }
            Ok(None)
        }
        fn eval_raw(&self, body: i32, epoch: f64) -> Result<(Vec<f64>, Vec<f64>)> {
            if body == MERCURY && epoch.to_bits() == FRAME_J2000.to_bits() {
                self.central_reads.set(self.central_reads.get() + 1);
            }
            self.model.eval_raw(body, epoch)
        }
        fn can_reuse_scalar_state(&self, _body: i32, jd: f64, epoch: TwoPartEpoch) -> bool {
            self.proof_calls.set(self.proof_calls.get() + 1);
            let Ok(canonical) = TwoPartEpoch::split(jd) else {
                return false;
            };
            // This affine test provider has its own exact input argument;
            // it does not inherit any LEB or center-composition capability.
            self.reuse_states
                && epoch.hi.to_bits() == canonical.hi.to_bits()
                && epoch.lo.to_bits() == canonical.lo.to_bits()
                && epoch.value().to_bits() == jd.to_bits()
                && (jd - FRAME_J2000).to_bits() == ((epoch.hi - FRAME_J2000) + epoch.lo).to_bits()
        }

        fn eval_raw_two_part(
            &self,
            body: i32,
            epoch: TwoPartEpoch,
        ) -> Result<(Vec<f64>, Vec<f64>)> {
            self.split_queries.borrow_mut().push([
                body as u64,
                epoch.hi.to_bits(),
                epoch.lo.to_bits(),
            ]);
            if body == MERCURY {
                let read = self.split_reads.get() + 1;
                self.split_reads.set(read);
                if self.fail_split_read == Some(read) {
                    return Err(Error::MathDomain {
                        value: epoch.value(),
                    });
                }
            }
            self.model.eval_raw_two_part(body, epoch)
        }
        fn eval_relative_position_two_part(
            &self,
            target: i32,
            epoch: TwoPartEpoch,
            center: i32,
            center_epoch: TwoPartEpoch,
        ) -> Result<Option<[f64; 3]>> {
            self.relative_queries.borrow_mut().push([
                target as u64,
                epoch.hi.to_bits(),
                epoch.lo.to_bits(),
                center as u64,
                center_epoch.hi.to_bits(),
                center_epoch.lo.to_bits(),
            ]);
            self.model
                .eval_relative_position_two_part(target, epoch, center, center_epoch)
        }
    }

    fn light_time_words(solution: LightTimeSolution) -> Vec<u64> {
        [
            solution.geo_pos,
            solution.geo_vel,
            solution.retarded_pos,
            solution.retarded_vel,
        ]
        .into_iter()
        .flatten()
        .chain([solution.light_time_days, solution.retarded_jd_tt])
        .map(f64::to_bits)
        .collect()
    }

    fn test_observer() -> PropagationObserver {
        PropagationObserver {
            position: [1.0, -0.8, 0.2],
            velocity: [0.0; 3],
            body: Some(EARTH),
            site_offset: [0.0; 3],
        }
    }

    #[test]
    fn observation_state_reuse_removes_only_the_first_split_read() {
        let mut model = CountedModel::new();
        model.reuse_states = true;
        let jd = FRAME_J2000 + 0.25;
        let canonical = TwoPartEpoch::split(jd).unwrap();
        let initial = model.geometric_state(MERCURY, jd).unwrap();
        for (parts, reusable) in [
            (None, true),
            (Some((canonical.hi, canonical.lo)), true),
            (Some((jd, 1e-12)), false),
        ] {
            for apply in [false, true] {
                model.reset_queries();
                let expected =
                    provider_light_time(&model, MERCURY, jd, parts, test_observer(), apply)
                        .unwrap();
                let expected_queries = model.split_queries.borrow().clone();
                let expected_relative = model.relative_queries.borrow().clone();
                assert_eq!(model.split_reads.get(), if apply { 4 } else { 1 });
                model.reset_queries();
                let actual = provider_light_time_with_initial(
                    &model,
                    MERCURY,
                    jd,
                    parts,
                    test_observer(),
                    apply,
                    Some(initial),
                )
                .unwrap();
                assert_eq!(light_time_words(actual), light_time_words(expected));
                let remaining = if reusable {
                    &expected_queries[1..]
                } else {
                    &expected_queries[..]
                };
                assert_eq!(&*model.split_queries.borrow(), remaining);
                assert_eq!(*model.relative_queries.borrow(), expected_relative);
                assert_eq!(model.proof_calls.get(), 1);
            }
        }
    }

    #[test]
    fn observation_state_reuse_keeps_observer_and_epoch_validation_before_proof() {
        let mut model = CountedModel::new();
        model.reuse_states = true;
        let initial = model.geometric_state(MERCURY, FRAME_J2000).unwrap();
        let observer = test_observer();
        let bad = f64::from_bits(0x7ff8_1234_5678_0001);
        let nonfinite_word = |error| match error {
            Error::NonFiniteArgument { value } => value.to_bits(),
            other => panic!("expected nonfinite argument, got {other:?}"),
        };
        for (parts, observer) in [
            (Some((bad, 0.0)), observer),
            (Some((FRAME_J2000, bad)), observer),
            (
                None,
                PropagationObserver {
                    position: [bad, 0.0, 0.0],
                    ..observer
                },
            ),
            (
                None,
                PropagationObserver {
                    velocity: [0.0, bad, 0.0],
                    ..observer
                },
            ),
            (
                None,
                PropagationObserver {
                    site_offset: [0.0, 0.0, bad],
                    ..observer
                },
            ),
        ] {
            model.reset_queries();
            let expected =
                provider_light_time(&model, MERCURY, FRAME_J2000, parts, observer, false)
                    .unwrap_err();
            let actual = provider_light_time_with_initial(
                &model,
                MERCURY,
                FRAME_J2000,
                parts,
                observer,
                false,
                Some(initial),
            )
            .unwrap_err();
            assert_eq!(nonfinite_word(actual), nonfinite_word(expected));
            assert_eq!(model.proof_calls.get(), 0);
            assert_eq!(model.split_reads.get(), 0);
        }
    }

    #[test]
    fn observation_state_reuse_rejects_seed_body_epoch_and_signed_zero_changes() {
        let mut model = CountedModel::new();
        model.reuse_states = true;
        let initial = model.geometric_state(MERCURY, FRAME_J2000).unwrap();
        for initial in [
            GeometricState {
                body_id: VENUS,
                ..initial
            },
            GeometricState {
                jd_tt: f64::from_bits(FRAME_J2000.to_bits() ^ 1),
                ..initial
            },
        ] {
            model.reset_queries();
            let expected =
                provider_light_time(&model, MERCURY, FRAME_J2000, None, test_observer(), false)
                    .unwrap();
            model.reset_queries();
            let actual = provider_light_time_with_initial(
                &model,
                MERCURY,
                FRAME_J2000,
                None,
                test_observer(),
                false,
                Some(initial),
            )
            .unwrap();
            assert_eq!(light_time_words(actual), light_time_words(expected));
            assert_eq!(model.split_reads.get(), 1);
            assert_eq!(model.proof_calls.get(), 0);
        }
        let negative_zero = model.geometric_state(MERCURY, -0.0).unwrap();
        model.reset_queries();
        provider_light_time_with_initial(
            &model,
            MERCURY,
            0.0,
            None,
            test_observer(),
            false,
            Some(negative_zero),
        )
        .unwrap();
        assert_eq!(model.split_reads.get(), 1);
        assert_eq!(model.proof_calls.get(), 0);
    }

    #[test]
    fn observation_state_reuse_keeps_the_opaque_second_target_read_failure() {
        let mut model = CountedModel::new();
        model.reuse_states = true;
        model.fail_split_read = Some(1);
        let nutation = |_| Ok((0.0, 0.0));
        // Even an opted-in provider cannot bypass its ordinary calls through
        // the public callback-bound constructor.
        let serving = SealedServing::new(&model, &nutation);
        let (outcome, events) = crate::trace::capture_events(|| {
            serving_compute_tt(&serving, FRAME_J2000, MERCURY, FLG_TRUEPOS | FLG_J2000)
        });
        assert!(matches!(
            outcome,
            Err(Error::MathDomain { value: FRAME_J2000 })
        ));
        assert_eq!(model.central_reads.get(), 1);
        assert_eq!(model.split_reads.get(), 1);
        assert_eq!(model.proof_calls.get(), 0);
        assert!(events.is_empty());
    }

    #[test]
    fn observation_state_reuse_preserves_core_words_flags_and_trace_sequence() {
        let original = CountedModel::new();
        let mut optimized = CountedModel::new();
        optimized.reuse_states = true;
        let nutation = |_| Ok((0.0, 0.0));
        let baseline = SealedServing::new_stored(&original, &nutation);
        let candidate = SealedServing::new_stored(&optimized, &nutation);
        // The central-angular capability is identical on both sides; only
        // observation-state reuse differs in this comparison.
        for epoch in [FRAME_J2000 - 0.25, FRAME_J2000, FRAME_J2000 + 0.25] {
            let canonical = TwoPartEpoch::split(epoch).unwrap();
            for flags in [
                FLG_TRUEPOS | FLG_J2000 | FLG_ICRS,
                FLG_SPEED | FLG_TRUEPOS | FLG_J2000,
                FLG_SPEED | FLG_NOABERR | FLG_NOGDEFL,
                FLG_SPEED | FLG_NOABERR | FLG_NOGDEFL | FLG_EQUATORIAL,
                FLG_SPEED | FLG_TRUEPOS | FLG_EQUATORIAL | FLG_NONUT,
                FLG_SPEED | FLG_SPEED3 | FLG_TRUEPOS | FLG_EQUATORIAL,
            ] {
                for order in [
                    CentralStateOrder::TargetFirst,
                    CentralStateOrder::ObserverFirst,
                ] {
                    for xyz in [false, true] {
                        for parts in [
                            None,
                            Some((canonical.hi, canonical.lo)),
                            Some((epoch, 1e-12)),
                        ] {
                            let run = |serving| {
                                crate::trace::capture_events(|| {
                                    serving_position_core_ordered(
                                        serving, epoch, MERCURY, flags, xyz, parts, order,
                                    )
                                    .unwrap()
                                })
                            };
                            let (expected, expected_events) = run(&baseline);
                            let (actual, actual_events) = run(&candidate);
                            assert_eq!(actual.0.map(f64::to_bits), expected.0.map(f64::to_bits));
                            assert_eq!(actual.1, expected.1);
                            assert_eq!(actual_events, expected_events);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn central_angular_reuse_requires_explicit_capability_and_identical_query_words() {
        let model = CountedModel::new();
        let nutation = |_| Ok((0.0, 0.0));
        let opaque = SealedServing::new(&model, &nutation);
        let immutable = SealedServing::new_stored(&model, &nutation);
        let flags = FLG_SPEED | FLG_TRUEPOS | FLG_J2000;
        let order = CentralStateOrder::TargetFirst;
        let query = CentralAngularQuery::new(FRAME_J2000, None, MERCURY, flags, order);
        let frames = frame_matrices_for_date(FRAME_J2000, 0.0, 0.0).unwrap();
        let central = CentralAngularPlace {
            query,
            apparent: [2.0, 0.2, 0.3],
            correction: [0.0; 3],
            frames: &frames,
            reads: [None; DEFLECTORS.len() * 2],
        };
        assert!(
            central_angular_sample(&opaque, &central, query)
                .unwrap()
                .is_none()
        );
        assert!(
            central_angular_sample(&immutable, &central, query)
                .unwrap()
                .is_some()
        );
        for changed in [
            CentralAngularQuery {
                tt: query.tt ^ 1,
                ..query
            },
            CentralAngularQuery {
                parts: Some([query.tt, 0]),
                ..query
            },
            CentralAngularQuery {
                body: VENUS,
                ..query
            },
            CentralAngularQuery {
                flags: query.flags | FLG_EQUATORIAL,
                ..query
            },
            CentralAngularQuery {
                order: CentralStateOrder::ObserverFirst,
                ..query
            },
        ] {
            assert!(
                central_angular_sample(&immutable, &central, changed)
                    .unwrap()
                    .is_none()
            );
        }
        let mut zero = query;
        zero.tt = 0.0_f64.to_bits();
        let mut minus_zero = zero;
        minus_zero.tt = (-0.0_f64).to_bits();
        assert!(zero != minus_zero);
        zero.parts = Some([query.tt, 0.0_f64.to_bits()]);
        minus_zero = zero;
        minus_zero.parts = Some([query.tt, (-0.0_f64).to_bits()]);
        assert!(zero != minus_zero);
        let stripped =
            CentralAngularQuery::new(FRAME_J2000, None, MERCURY, flags | FLG_SPEED3, order);
        assert!(stripped == query);
    }

    #[test]
    fn central_angular_reuse_removes_only_the_identical_central_scalar_read() {
        let model = CountedModel::new();
        let nutation = |_| Ok((0.0, 0.0));
        let flags = FLG_SPEED | FLG_TRUEPOS | FLG_J2000;
        for parts in [None, Some((FRAME_J2000, 0.0)), Some((FRAME_J2000, 1e-12))] {
            let opaque = SealedServing::new(&model, &nutation);
            let immutable = SealedServing::new_stored(&model, &nutation);
            let run = |serving| {
                model.central_reads.set(0);
                let value = serving_position_core_ordered(
                    serving,
                    FRAME_J2000,
                    MERCURY,
                    flags,
                    false,
                    parts,
                    CentralStateOrder::TargetFirst,
                )
                .unwrap();
                (value, model.central_reads.get())
            };
            let (expected, original_reads) = run(&opaque);
            let (actual, reused_reads) = run(&immutable);
            assert_eq!(actual.0.map(f64::to_bits), expected.0.map(f64::to_bits));
            assert_eq!(actual.1, expected.1);
            assert_eq!(original_reads, 2);
            assert_eq!(
                reused_reads,
                if parts == Some((FRAME_J2000, 1e-12)) {
                    2
                } else {
                    1
                }
            );
        }
        // Non-geocentric centers retain the original route even when the
        // immutable capability is present.
        for extra in [FLG_HELCTR, FLG_BARYCTR] {
            model.central_reads.set(0);
            serving_compute_tt(
                &SealedServing::new_stored(&model, &nutation),
                FRAME_J2000,
                MERCURY,
                flags | extra,
            )
            .unwrap();
            assert_eq!(model.central_reads.get(), 2);
        }
    }

    #[test]
    fn central_angular_reuse_preserves_all_output_words_and_trace_events() {
        let model = CountedModel::new();
        let nutation = |_| Ok((0.0, 0.0));
        let opaque = SealedServing::new(&model, &nutation);
        let immutable = SealedServing::new_stored(&model, &nutation);
        for epoch in [FRAME_J2000 - 0.5, FRAME_J2000, FRAME_J2000 + 0.5] {
            for flags in [
                FLG_SPEED | FLG_TRUEPOS,
                FLG_SPEED | FLG_TRUEPOS | FLG_J2000 | FLG_ICRS,
                FLG_SPEED | FLG_NOABERR | FLG_NOGDEFL,
                FLG_SPEED | FLG_NOABERR | FLG_NOGDEFL | FLG_EQUATORIAL,
                FLG_SPEED | FLG_TRUEPOS | FLG_EQUATORIAL | FLG_NONUT,
            ] {
                for order in [
                    CentralStateOrder::TargetFirst,
                    CentralStateOrder::ObserverFirst,
                ] {
                    for xyz in [false, true] {
                        for parts in [None, Some((epoch.floor(), epoch - epoch.floor()))] {
                            let run = |serving| {
                                crate::trace::capture_events(|| {
                                    serving_position_core_ordered(
                                        serving, epoch, MERCURY, flags, xyz, parts, order,
                                    )
                                    .unwrap()
                                })
                            };
                            let (expected, expected_events) = run(&opaque);
                            let (actual, actual_events) = run(&immutable);
                            assert_eq!(actual.0.map(f64::to_bits), expected.0.map(f64::to_bits));
                            assert_eq!(actual.1, expected.1);
                            assert_eq!(actual_events, expected_events);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn central_angular_reuse_keeps_opaque_second_callback_failure() {
        let model = CountedModel::new();
        let central_calls = std::cell::Cell::new(0);
        let nutation = |epoch: f64| {
            if epoch.to_bits() == FRAME_J2000.to_bits() {
                central_calls.set(central_calls.get() + 1);
                if central_calls.get() == 2 {
                    return Err(Error::MathDomain { value: epoch });
                }
            }
            Ok((0.0, 0.0))
        };
        let serving = SealedServing::new(&model, &nutation);
        assert!(matches!(
            serving_compute_tt(
                &serving,
                FRAME_J2000,
                MERCURY,
                FLG_SPEED | FLG_TRUEPOS | FLG_J2000
            ),
            Err(Error::MathDomain { .. })
        ));
        assert_eq!(central_calls.get(), 2);
        assert_eq!(model.central_reads.get(), 2);
    }

    #[test]
    fn central_angular_reuse_does_not_skip_an_opaque_topocentric_clock() {
        let mut model = CountedModel::new();
        model.reuse_states = true;
        let nutation = |_| Ok((0.0, 0.0));
        for serving in [
            SealedServing::new(&model, &nutation),
            SealedServing::new_stored(&model, &nutation),
        ] {
            let central_calls = std::cell::Cell::new(0);
            let c2t = |tt: f64, _ut: f64| {
                if tt.to_bits() == FRAME_J2000.to_bits() {
                    central_calls.set(central_calls.get() + 1);
                    if central_calls.get() == 2 {
                        return Err(Error::MathDomain { value: tt });
                    }
                }
                Ok([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
            };
            let outcome = serving_compute_tt_topocentric(
                &serving,
                FRAME_J2000,
                FRAME_J2000,
                MERCURY,
                FLG_SPEED | FLG_TRUEPOS | FLG_J2000 | FLG_TOPOCTR,
                (10.0, 45.0, 100.0),
                &c2t,
                false,
                SiteEllipsoid::Wgs84,
            );
            assert!(matches!(outcome, Err(Error::MathDomain { .. })));
            assert_eq!(central_calls.get(), 2);
            assert_eq!(model.proof_calls.get(), 0);
        }
    }

    #[test]
    fn central_angular_reuse_defers_domain_errors_and_propagates_trace_refusal() {
        let mut model = CountedModel::new();
        model.fail_domains = true;
        let nutation = |_| Ok((0.0, 0.0));
        let serving = SealedServing::new_stored(&model, &nutation);
        let query = CentralAngularQuery::new(
            FRAME_J2000,
            None,
            MERCURY,
            FLG_SPEED,
            CentralStateOrder::TargetFirst,
        );
        let frames = frame_matrices_for_date(FRAME_J2000, 0.0, 0.0).unwrap();
        let mut central = CentralAngularPlace {
            query,
            apparent: [2.0, 0.2, 0.3],
            correction: [0.0; 3],
            frames: &frames,
            reads: [None; DEFLECTORS.len() * 2],
        };
        central.reads[0] = Some(CorrectionRead {
            body: SUN,
            epoch: FRAME_J2000,
            interval: None,
        });
        assert_eq!(model.domains.get(), 0);
        assert!(
            central_angular_sample(
                &serving,
                &central,
                CentralAngularQuery {
                    tt: query.tt ^ 1,
                    ..query
                }
            )
            .unwrap()
            .is_none()
        );
        assert_eq!(model.domains.get(), 0);
        let (outcome, events) =
            crate::trace::capture_events(|| central_angular_sample(&serving, &central, query));
        assert!(matches!(outcome, Err(Error::MathDomain { .. })));
        assert_eq!(model.domains.get(), 1);
        assert!(events.is_empty());
        central.reads = [None; DEFLECTORS.len() * 2];
        let sink: crate::trace::FallibleTraceEventSink =
            std::sync::Arc::new(|_, _| Err(Error::MathDomain { value: 17.0 }));
        let (outcome, events) = crate::trace::capture_events_with_fallible_sink(
            || central_angular_sample(&serving, &central, query),
            sink,
        );
        assert!(matches!(outcome, Err(Error::MathDomain { value: 17.0 })));
        assert_eq!(
            events,
            vec![crate::trace::TraceEvent::Record {
                body: MERCURY,
                source: "LEB".to_owned()
            }]
        );
    }

    #[test]
    fn repeated_speed_keeps_the_actual_nutation_callback_and_its_failures() {
        let model = InertialModel {
            origin: [0.0; 3],
            position: [2.0, 0.2, 0.3],
            velocity: [0.005, 0.002, -0.001],
        };
        let calls = std::cell::Cell::new(0);
        let fail = std::cell::Cell::new(false);
        let nutation = |epoch| {
            calls.set(calls.get() + 1);
            if fail.get() {
                return Err(Error::MathDomain { value: epoch });
            }
            Ok((0.0, 0.0))
        };
        let serving = SealedServing::new(&model, &nutation);
        let flags = FLG_SPEED | FLG_TRUEPOS;
        let first = serving_compute_tt(&serving, FRAME_J2000, MERCURY, flags).unwrap();
        let first_calls = calls.get();
        let second = serving_compute_tt(&serving, FRAME_J2000, MERCURY, flags).unwrap();
        assert!(calls.get() > first_calls);
        assert_eq!(first.0.map(f64::to_bits), second.0.map(f64::to_bits));
        assert_eq!(first.1, second.1);
        fail.set(true);
        assert!(matches!(
            serving_compute_tt(&serving, FRAME_J2000, MERCURY, flags),
            Err(Error::MathDomain { .. })
        ));
    }

    #[test]
    fn apparent_rates_differentiate_each_retarded_target_evaluation() {
        // An inertial subluminal trajectory has an exact finite-iteration
        // light-time chain rule, independent of any angular difference.
        let model = InertialModel {
            origin: [0.0; 3],
            position: [2.0, 0.2, 0.3],
            velocity: [0.5, 0.2, -0.1],
        };
        let serving = SealedServing::new(&model, &|_| Ok((0.0, 0.0)));
        let flags = FLG_SPEED | FLG_NOABERR | FLG_NOGDEFL | FLG_J2000 | FLG_ICRS;
        let actual = serving_compute_tt(&serving, FRAME_J2000, MERCURY, flags)
            .unwrap()
            .0;
        let mut delay = 0.0;
        let mut emission_rate = 1.0;
        for _ in 0..LIGHT_TIME_ITERATIONS {
            let p: [f64; 3] =
                std::array::from_fn(|i| model.position[i] - model.velocity[i] * delay);
            let radius = p.iter().map(|x| x * x).sum::<f64>().sqrt();
            let delay_rate = p
                .iter()
                .zip(model.velocity)
                .map(|(p, v)| p * v * emission_rate)
                .sum::<f64>()
                / (radius * C_LIGHT_AU_DAY);
            delay = radius / C_LIGHT_AU_DAY;
            emission_rate = 1.0 - delay_rate;
        }
        let ray = std::array::from_fn(|i| model.position[i] - model.velocity[i] * delay);
        let expected = angular_jacobian(ray, model.velocity.map(|v| v * emission_rate), flags);
        for component in 0..2 {
            assert!(
                (actual[component + 3] - expected[component]).abs()
                    <= expected[component].abs() * 1e-6,
                "component {component}: actual {}, exact {}",
                actual[component + 3],
                expected[component]
            );
        }
    }

    #[test]
    fn retarded_epoch_evaluates_the_fraction_below_a_full_julian_day_ulp() {
        let model = InertialModel {
            origin: [0.0; 3],
            position: [0.0027, 0.0001, -0.0003],
            velocity: [0.017, -0.005, 0.001],
        };
        let actual = provider_light_time(
            &model,
            MERCURY,
            FRAME_J2000,
            None,
            PropagationObserver {
                position: [0.0; 3],
                velocity: [0.0; 3],
                body: Some(EARTH),
                site_offset: [0.0; 3],
            },
            true,
        )
        .unwrap();
        // This independent affine trajectory uses only the small time
        // displacement, so none of its propagation epochs can lose a day ULP.
        let mut delay = 0.0;
        let mut exact = model.position;
        for _ in 0..LIGHT_TIME_ITERATIONS {
            delay = exact.iter().map(|x| x * x).sum::<f64>().sqrt() / C_LIGHT_AU_DAY;
            exact = std::array::from_fn(|i| model.position[i] - model.velocity[i] * delay);
        }
        for (actual, expected) in actual.geo_pos().into_iter().zip(exact) {
            assert!((actual - expected).abs() <= 32.0 * f64::EPSILON * expected.abs());
        }
        let rounded_delay = FRAME_J2000 - (FRAME_J2000 - delay);
        assert!((rounded_delay - delay).abs() > 1e-12);
        assert!(
            model
                .velocity
                .iter()
                .any(|v| (v * (rounded_delay - delay)).abs() > 1e-15)
        );
    }

    #[test]
    fn paired_propagation_keeps_small_geometry_and_subtracts_the_actual_site_once() {
        let model = InertialModel {
            origin: [1.0, -0.8, 0.2],
            position: [0.0027, 0.0001, -0.0003],
            velocity: [0.017, -0.005, 0.001],
        };
        let site_offset = [0.00003, -0.00002, 0.00001];
        let observer = PropagationObserver {
            position: std::array::from_fn(|i| model.origin[i] + site_offset[i]),
            velocity: [0.0; 3],
            body: Some(EARTH),
            site_offset,
        };
        let actual =
            provider_light_time(&model, MERCURY, FRAME_J2000, None, observer, true).unwrap();
        let mut exact: [f64; 3] = std::array::from_fn(|i| model.position[i] - site_offset[i]);
        for _ in 0..LIGHT_TIME_ITERATIONS {
            let delay = exact.iter().map(|x| x * x).sum::<f64>().sqrt() / C_LIGHT_AU_DAY;
            exact = std::array::from_fn(|i| {
                model.position[i] - model.velocity[i] * delay - site_offset[i]
            });
        }
        for (actual, expected) in actual.geo_pos().into_iter().zip(exact) {
            assert!((actual - expected).abs() <= 32.0 * f64::EPSILON * expected.abs());
        }
    }

    #[test]
    fn surface_rates_include_the_supplied_ut1_clock_jacobian() {
        let model = InertialModel {
            origin: [0.0; 3],
            position: [1.0, 0.2, 0.3],
            velocity: [0.0; 3],
        };
        let serving = SealedServing::new(&model, &|_| Ok((0.0, 0.0)));
        let flags = FLG_SPEED | FLG_TRUEPOS | FLG_J2000 | FLG_ICRS | FLG_TOPOCTR;
        let ut0 = FRAME_J2000 - 64.0 / 86400.0;
        let clock_rate = 1.001;
        let ut1_at = |tt| Ok(ut0 + (tt - FRAME_J2000) * clock_rate);
        let omega = std::f64::consts::TAU;
        let c2t_at = |_, ut: f64| {
            let (sin, cos) = (omega * (ut - ut0)).sin_cos();
            Ok([[cos, -sin, 0.0], [sin, cos, 0.0], [0.0, 0.0, 1.0]])
        };
        let actual = serving_compute_topocentric_with_clock(
            &serving,
            FRAME_J2000,
            ut0,
            MERCURY,
            flags,
            (0.0, 0.0, 0.0),
            &c2t_at,
            false,
            &ut1_at,
            (f64::NEG_INFINITY, f64::INFINITY),
        )
        .unwrap()
        .0;
        let radius = 6378137.0 / 149597870700.0;
        let ray = [
            model.position[0] - radius,
            model.position[1],
            model.position[2],
        ];
        let expected = angular_jacobian(ray, [0.0, omega * radius * clock_rate, 0.0], flags);
        for component in 0..2 {
            assert!(
                (actual[component + 3] - expected[component]).abs()
                    <= expected[component].abs() * 1e-6,
                "component {component}: actual {}, exact {}",
                actual[component + 3],
                expected[component]
            );
        }
    }
}

#[cfg(test)]
#[path = "fast_angular_derivative_tests.rs"]
mod fast_angular_derivative_tests;

#[cfg(test)]
#[path = "ppn_scalar_replay_tests.rs"]
mod ppn_scalar_replay_tests;

/// Immutable observing center and its actual surface displacement.
#[derive(Clone, Copy)]
struct PropagationObserver {
    position: [f64; 3],
    velocity: [f64; 3],
    body: Option<i32>,
    site_offset: [f64; 3],
}

/// Solve propagation without collapsing a fractional emission epoch.
///
/// All bodies use the same three fixed-point iterations. The observation's
/// exact binary64 value (or supplied two-part value) is retained, and each
/// light time is subtracted through error-free epoch addition before the
/// provider evaluates its coefficients. This avoids full-Julian-day rounding
/// of subsecond propagation changes. State vectors are AU and AU/TT-day;
/// propagation time is in days. Every epoch/provider/domain error propagates.
fn provider_light_time<P: StateProvider>(
    provider: &P,
    body: i32,
    jd_tt: f64,
    observation_parts: Option<(f64, f64)>,
    observer: PropagationObserver,
    apply: bool,
) -> Result<LightTimeSolution> {
    provider_light_time_with_initial(
        provider,
        body,
        jd_tt,
        observation_parts,
        observer,
        apply,
        None,
    )
}

/// Reuse only the observation target; every relative and retarded read keeps
/// the established route. The seed belongs to this same immutable serving.
fn provider_light_time_with_initial<P: StateProvider>(
    provider: &P,
    body: i32,
    jd_tt: f64,
    observation_parts: Option<(f64, f64)>,
    observer: PropagationObserver,
    apply: bool,
    initial: Option<GeometricState>,
) -> Result<LightTimeSolution> {
    let observation = match observation_parts {
        Some((hi, lo)) => {
            if !hi.is_finite() || !lo.is_finite() {
                return Err(Error::NonFiniteArgument { value: hi + lo });
            }
            TwoPartEpoch { hi, lo }
        }
        None => TwoPartEpoch::split(jd_tt)?,
    };
    super::frames::check_finite_vec3(observer.position)?;
    super::frames::check_finite_vec3(observer.velocity)?;
    super::frames::check_finite_vec3(observer.site_offset)?;
    let relative = |state: GeometricState, epoch| -> Result<[f64; 3]> {
        if let Some(center) = observer.body {
            if let Some(ray) =
                provider.eval_relative_position_two_part(body, epoch, center, observation)?
            {
                return Ok(std::array::from_fn(|i| ray[i] - observer.site_offset[i]));
            }
        }
        Ok(std::array::from_fn(|i| {
            state.pos()[i] - observer.position[i]
        }))
    };
    // Keep all observation/observer validation before this optional reuse.
    // The scalar read already succeeded through the same immutable provider;
    // the concrete route must additionally prove normalized-input identity.
    let mut state = match initial {
        Some(state)
            if state.body_id() == body
                && state.jd_tt().to_bits() == jd_tt.to_bits()
                && provider.can_reuse_scalar_state(body, jd_tt, observation) =>
        {
            state
        }
        _ => provider.geometric_state_two_part(body, observation)?,
    };
    let mut geo = relative(state, observation)?;
    let mut velocity = std::array::from_fn(|i| state.vel()[i] - observer.velocity[i]);
    let mut delay = 0.0;
    let mut delay_rate = 0.0;
    let mut emission = observation;
    if apply {
        for _ in 0..LIGHT_TIME_ITERATIONS {
            let distance = super::frames::vec3_dist(geo);
            if distance == 0.0 {
                break;
            }
            delay_rate = super::planet_range::RangeState::from_cartesian(geo, velocity).rate
                / C_LIGHT_AU_DAY;
            let emission_rate = 1.0 - delay_rate;
            delay = distance / C_LIGHT_AU_DAY;
            emission = observation.add_offset(-delay)?;
            state = provider.geometric_state_two_part(body, emission)?;
            geo = relative(state, emission)?;
            velocity =
                std::array::from_fn(|i| state.vel()[i] * emission_rate - observer.velocity[i]);
        }
    }
    Ok(LightTimeSolution {
        geo_pos: geo,
        geo_vel: velocity,
        retarded_pos: state.pos(),
        retarded_vel: state.vel(),
        light_time_days: delay,
        light_time_rate: delay_rate,
        retarded_jd_tt: emission.value(),
    })
}

/// Continuous observation-epoch offsets for the actual acquired geometry.
/// A target is read both at observation and emission; their selected polynomial
/// intervals must remain owned by the same local model. Observer and deflector
/// intervals are in observation TT. Nutation support is included by the caller.
fn angular_geometry_support<P: StateProvider>(
    provider: &P,
    body: i32,
    observer: Option<i32>,
    jd_tt: f64,
    emission_tt: f64,
) -> Result<(f64, f64)> {
    let mut support = (f64::NEG_INFINITY, f64::INFINITY);
    let mut include = |body, query| -> Result<()> {
        for (start, end) in [
            provider.body_coverage(body),
            provider.smooth_interval_tt(body, query)?,
        ]
        .into_iter()
        .flatten()
        {
            support.0 = support.0.max(start - query);
            support.1 = support.1.min(end - query);
        }
        Ok(())
    };
    include(body, jd_tt)?;
    include(body, emission_tt)?;
    if let Some(observer) = observer {
        include(observer, jd_tt)?;
    }
    Ok(support)
}

/// Apply the spherical-to-Cartesian position and velocity Jacobian.
///
/// Inputs are longitude/latitude in degrees, radius AU, angular rates in
/// degrees/day and radial rate AU/day. Output positions are AU and velocities
/// AU/day; all-zero input rates yield exactly zero Cartesian velocity.
pub fn spherical_to_cartesian_with_velocity(
    lon_deg: f64,
    lat_deg: f64,
    dist_au: f64,
    dlon_deg_day: f64,
    dlat_deg_day: f64,
    ddist_au_day: f64,
) -> PositionResult {
    let lon_r = lon_deg.to_radians();
    let lat_r = lat_deg.to_radians();
    let cos_lat = lat_r.cos();
    let sin_lat = lat_r.sin();
    let cos_lon = lon_r.cos();
    let sin_lon = lon_r.sin();
    let x = dist_au * cos_lat * cos_lon;
    let y = dist_au * cos_lat * sin_lon;
    let z = dist_au * sin_lat;
    let (vx, vy, vz) = if dlon_deg_day != 0.0 || dlat_deg_day != 0.0 || ddist_au_day != 0.0 {
        let dlon_r = dlon_deg_day.to_radians();
        let dlat_r = dlat_deg_day.to_radians();
        (
            -dist_au * cos_lat * sin_lon * dlon_r - dist_au * sin_lat * cos_lon * dlat_r
                + cos_lat * cos_lon * ddist_au_day,
            dist_au * cos_lat * cos_lon * dlon_r - dist_au * sin_lat * sin_lon * dlat_r
                + cos_lat * sin_lon * ddist_au_day,
            dist_au * cos_lat * dlat_r + sin_lat * ddist_au_day,
        )
    } else {
        (0.0, 0.0, 0.0)
    };
    [x, y, z, vx, vy, vz]
}

/// Apply the requested representation to a spherical result.
///
/// XYZ takes precedence and returns Cartesian AU/AU-day. Otherwise RADIANS
/// converts only angular slots and rates. Nutation pseudo-body results bypass
/// this transformation.
pub fn apply_output_flags(result: PositionResult, flags: i32) -> PositionResult {
    if (flags & FLG_XYZ) != 0 {
        return spherical_to_cartesian_with_velocity(
            result[0], result[1], result[2], result[3], result[4], result[5],
        );
    }
    if (flags & FLG_RADIANS) != 0 {
        return [
            result[0].to_radians(),
            result[1].to_radians(),
            result[2],
            result[3].to_radians(),
            result[4].to_radians(),
            result[5],
        ];
    }
    result
}

/// Apply representation shaping and include its implied public return flags.
pub fn finalize_output_flags(
    pos: PositionResult,
    retflag: i32,
    flags: i32,
) -> (PositionResult, i32) {
    (
        apply_output_flags(pos, flags),
        retflag | (flags & (FLG_XYZ | FLG_RADIANS)) | implied_retflag_bits(flags),
    )
}

/// Construct the descending-node antipode from the formatted ascending node.
///
/// Cartesian output negates position and velocity. Spherical output adds a
/// half turn to longitude and negates latitude and its rate, leaving radius
/// and radial rate unchanged. `flags` specifies the output representation.
pub fn south_node_from_north(north: PositionResult, flags: i32) -> PositionResult {
    if (flags & FLG_XYZ) != 0 {
        return [
            -north[0], -north[1], -north[2], -north[3], -north[4], -north[5],
        ];
    }
    let half_turn = if (flags & FLG_RADIANS) != 0 {
        std::f64::consts::PI
    } else {
        180.0
    };
    let full_turn = 2.0 * half_turn;
    [
        (north[0] + half_turn).rem_euclid(full_turn),
        -north[1],
        north[2],
        north[3],
        -north[4],
        north[5],
    ]
}

/// Return the six-zero sentinel for a topocentric Earth request.
/// This public origin convention does not expose the geocenter-to-site offset.
pub fn degenerate_origin_result() -> PositionResult {
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
}

/// Rebuild rates from positions at equal TT offsets using a centered
/// difference. Spherical longitude is unwrapped about the central position;
/// Cartesian output differentiates each axis independently.
pub fn apply_central_speed3_stencil(
    current: PositionResult,
    previous: PositionResult,
    following: PositionResult,
    flags: i32,
) -> PositionResult {
    let (v0, v1, v2) = if (flags & FLG_XYZ) != 0 {
        (
            (following[0] - previous[0]) / (2.0 * SPEED3_STEP_DAYS),
            (following[1] - previous[1]) / (2.0 * SPEED3_STEP_DAYS),
            (following[2] - previous[2]) / (2.0 * SPEED3_STEP_DAYS),
        )
    } else {
        let half_cycle = if (flags & FLG_RADIANS) != 0 {
            std::f64::consts::PI
        } else {
            180.0
        };
        let full_cycle = 2.0 * half_cycle;
        let previous_delta =
            (previous[0] - current[0] + half_cycle).rem_euclid(full_cycle) - half_cycle;
        let following_delta =
            (following[0] - current[0] + half_cycle).rem_euclid(full_cycle) - half_cycle;
        (
            (following_delta - previous_delta) / (2.0 * SPEED3_STEP_DAYS),
            (following[1] - previous[1]) / (2.0 * SPEED3_STEP_DAYS),
            (following[2] - previous[2]) / (2.0 * SPEED3_STEP_DAYS),
        )
    };
    [current[0], current[1], current[2], v0, v1, v2]
}

/// Whether a body uses the SPEED3 three-position derivative convention
/// (`planets._rebuild_mean_point_speed3` gate): the mean lunar node and
/// mean apogee only. All other bodies keep the core-computed speeds.
pub fn uses_mean_point_speed3(body_id: i32) -> bool {
    body_id == MEAN_NODE || body_id == MEAN_APOG
}

/// Differentiate a mean lunar point after complete output-frame reduction.
///
/// `sample` maps a signed day offset to the final no-speed result, keeping
/// nutation, sidereal and representation shaping inside the same stencil.
/// The callback must disable SPEED and SPEED3. Non-mean points pass through;
/// provider and frame errors propagate.
pub fn rebuild_mean_point_speed3(
    current: PositionResult,
    body_id: i32,
    flags: i32,
    sample: impl Fn(f64) -> Result<PositionResult>,
) -> Result<PositionResult> {
    if (flags & FLG_SPEED3) == 0 || !uses_mean_point_speed3(body_id) {
        return Ok(current);
    }
    let previous = sample(-SPEED3_STEP_DAYS)?;
    let following = sample(SPEED3_STEP_DAYS)?;
    Ok(apply_central_speed3_stencil(
        current, previous, following, flags,
    ))
}

/// Shared calculation tail after the observation has been expressed in TT.
/// Run the origin sentinel or provider, shape output, rebuild mean-point
/// SPEED3 when requested and return the original request's flag labels.
pub(super) fn calc_tt_body(
    engine: &CalcEngine<'_>,
    jd_tt: f64,
    body_id: i32,
    raw_flags: i32,
    norm_flags: i32,
) -> Result<(PositionResult, i32)> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    calc_tt_body_source_time(engine, jd_tt, body_id, raw_flags, norm_flags)
}

/// Reduce a model after the public observation epoch has been validated.
///
/// An internally derived TT is passed to the owning model, provider and frame
/// checks; the public finite-input guard has already run. The shared tail
/// applies output representation, SPEED3 and returned flags. General checked
/// engine callers use [`calc_tt_body`].
pub(super) fn calc_tt_body_source_time(
    engine: &CalcEngine<'_>,
    jd_tt: f64,
    body_id: i32,
    raw_flags: i32,
    norm_flags: i32,
) -> Result<(PositionResult, i32)> {
    calc_tt_body_source_time_with_echo(engine, jd_tt, body_id, raw_flags, norm_flags, false)
}

/// Same numerical tail while an attached caller owns the final raw echo.
pub(super) fn calc_tt_body_source_time_with_echo(
    engine: &CalcEngine<'_>,
    jd_tt: f64,
    body_id: i32,
    raw_flags: i32,
    norm_flags: i32,
    defer_echo: bool,
) -> Result<(PositionResult, i32)> {
    let calc_iflag = strip_output_flags(norm_flags);
    // Degenerate topocentric Earth (compatibility contract): the exact
    // zero vector. Routed past the LEB fast path, which would otherwise
    // return the geocenter-from-topocenter offset.
    let (pos, retflag) = if body_id == EARTH && (norm_flags & FLG_TOPOCTR) != 0 {
        super::routing::mark_native_trace_source("Analytical");
        (degenerate_origin_result(), calc_iflag)
    } else {
        (engine.compute_tt)(jd_tt, body_id, calc_iflag)?
    };
    let (pos, retflag) = finalize_output_flags(pos, retflag, norm_flags);
    // SPEED3 differentiates the complete no-speed reduction at symmetric TT offsets.
    let sample_flags = norm_flags & !(FLG_SPEED | FLG_SPEED3);
    let sample_calc_flags = strip_output_flags(sample_flags);
    let pos = rebuild_mean_point_speed3(pos, body_id, norm_flags, |offset| {
        let (sample_pos, sample_retflag) =
            (engine.compute_tt)(jd_tt + offset, body_id, sample_calc_flags)?;
        Ok(finalize_output_flags(sample_pos, sample_retflag, sample_flags).0)
    })?;
    Ok((
        pos,
        if defer_echo {
            retflag
        } else {
            echo_request_bits(retflag, raw_flags)
        },
    ))
}

/// Typed flag words already consumed for one numerical entry.
/// `raw` follows CENTER_BODY/JPLHOR consumption; `normalized` retains the
/// single calculation normalization. A caller may defer restoring request
/// bits until it receives the numerical/backend flags. Construction performs
/// no additional validation, normalization or calculation.
#[derive(Clone, Copy, Debug)]
pub struct PreparedCalcFlags {
    /// Post-consumption request word used when restoring public flag bits.
    pub raw: i32,
    /// Single-normalized computation/output word, including implied SPEED.
    pub normalized: i32,
    /// Return numerical/backend flags without restoring request bits.
    /// False restores them using `raw`; true leaves that operation to the caller.
    pub defer_echo: bool,
}

/// Calculate a body at a UT1 Julian Day using the bound engine.
///
/// The body identifier includes nutation pseudo-body and descending-node
/// aliases. The engine supplies Delta-T, provider selection and sidereal
/// models. The result is the selected six-component layout and public flags.
/// Configuration, epoch, model and provider errors propagate.
pub fn calc_ut(
    engine: &CalcEngine<'_>,
    tjdut: f64,
    planet: i32,
    flags: i32,
) -> Result<(PositionResult, i32)> {
    calc_ut_flag_entry(engine, tjdut, planet, flags, None)
}

/// UT1 calculation with request flags already validated and normalized.
/// The time, body, output and SPEED3 reduction is the same as [`calc_ut`];
/// this entry avoids repeating the flag prelude.
pub fn calc_ut_prepared(
    engine: &CalcEngine<'_>,
    tjdut: f64,
    planet: i32,
    flags: PreparedCalcFlags,
) -> Result<(PositionResult, i32)> {
    calc_ut_flag_entry(engine, tjdut, planet, flags.raw, Some(flags))
}

fn calc_ut_flag_entry(
    engine: &CalcEngine<'_>,
    tjdut: f64,
    planet: i32,
    flags: i32,
    prepared: Option<PreparedCalcFlags>,
) -> Result<(PositionResult, i32)> {
    if !tjdut.is_finite() {
        return Err(Error::NonFiniteArgument { value: tjdut });
    }
    if prepared.is_none() {
        validate_barycentric_moseph(flags, "calc_ut")?;
    }
    // ECL_NUT runs on raw flags before normalization: one exclusive
    // ephemeris bit, no SPEED3 remap beyond the SPEED-priority strip, no
    // output shaping — the nutation tuple is returned as computed, with
    // only the implied bits added to the echo.
    if planet == ECL_NUT {
        let mut res_flags = resolve_center_flags(flags);
        if (res_flags & FLG_SPEED) != 0 {
            res_flags &= !FLG_SPEED3;
        }
        let nut_flags = exclusive_ephemeris_bit(res_flags, engine.jpl_fallback_active);
        let (pos_nut, rf_nut) = (engine.nutation_obliquity_ut)(tjdut, nut_flags)?;
        crate::trace::try_record(planet, "ERFA")?;
        return Ok((pos_nut, rf_nut | implied_retflag_bits(res_flags)));
    }
    let requested_planet = planet;
    let raw_flags = prepared.map_or_else(
        || consume_non_echoed_flags(flags, planet),
        |flags| flags.raw,
    );
    let norm_flags = prepared.map_or_else(
        || normalize_calc_flags(raw_flags, engine.jpl_fallback_active),
        |flags| flags.normalized,
    );
    // Resolve built-in asteroid aliases before selecting their physical provider.
    let planet = remap_ast_offset(planet);
    // South nodes derive from the north node through the same dispatch
    // path (backend consistency), then the antipode map in the
    // normalized output representation.
    if planet == -MEAN_NODE || planet == -TRUE_NODE {
        // Move the recursive ascending-node source record to the requested
        // ascending/descending alias and restore the previous ascending record.
        // Inactive tracing performs no publication.
        let north = planet.abs();
        let previous = crate::trace::try_snapshot_record(north)?;
        let (north_result, retflag) = match calc_ut(engine, tjdut, north, norm_flags) {
            Ok(result) => result,
            Err(error) => {
                crate::trace::with_caught_error(&error, || {
                    crate::trace::try_restore_record(north, previous)
                })?;
                return Err(error);
            }
        };
        crate::trace::try_record_alias(planet, north, previous)?;
        return Ok((
            south_node_from_north(north_result, norm_flags),
            if prepared.is_some_and(|flags| flags.defer_echo) {
                retflag
            } else {
                echo_request_bits(retflag, raw_flags)
            },
        ));
    }
    super::routing::with_native_trace_body(requested_planet, planet, || {
        let jd_tt = (engine.ut1_to_tt_jd)(tjdut)?;
        let result = if prepared.is_some_and(|flags| flags.defer_echo) {
            calc_tt_body_source_time_with_echo(engine, jd_tt, planet, raw_flags, norm_flags, true)?
        } else {
            calc_tt_body(engine, jd_tt, planet, raw_flags, norm_flags)?
        };
        super::routing::publish_native_trace_source(requested_planet)?;
        Ok(result)
    })
}

/// Calculate a body at a TT Julian Day using the bound engine.
///
/// The six-component layout and numerical error policy match [`calc_ut`].
/// TT output returns only the caller's ephemeris-selection labels.
pub fn calc(
    engine: &CalcEngine<'_>,
    tjdet: f64,
    planet: i32,
    flags: i32,
) -> Result<(PositionResult, i32)> {
    calc_tt_flag_entry(engine, tjdet, planet, flags, None)
}

/// TT counterpart of calc_ut_prepared. The already normalized words preserve
/// lone SPEED3 instead of treating its implied SPEED as an explicit request.
pub fn calc_prepared(
    engine: &CalcEngine<'_>,
    tjdet: f64,
    planet: i32,
    flags: PreparedCalcFlags,
) -> Result<(PositionResult, i32)> {
    calc_tt_flag_entry(engine, tjdet, planet, flags.raw, Some(flags))
}

fn calc_tt_flag_entry(
    engine: &CalcEngine<'_>,
    tjdet: f64,
    planet: i32,
    flags: i32,
    prepared: Option<PreparedCalcFlags>,
) -> Result<(PositionResult, i32)> {
    if !tjdet.is_finite() {
        return Err(Error::NonFiniteArgument { value: tjdet });
    }
    if prepared.is_none() {
        validate_barycentric_moseph(flags, "calc")?;
    }
    if planet == ECL_NUT {
        let mut res_flags = resolve_center_flags(flags);
        if (res_flags & FLG_SPEED) != 0 {
            res_flags &= !FLG_SPEED3;
        }
        let nut_flags = exclusive_ephemeris_bit(res_flags, engine.jpl_fallback_active);
        let (pos_nut, rf_nut) = (engine.nutation_obliquity_tt)(tjdet, nut_flags)?;
        crate::trace::try_record(planet, "ERFA")?;
        let retflag = calc_tt_epheflag_echo(rf_nut | implied_retflag_bits(res_flags), flags);
        return Ok((pos_nut, retflag));
    }
    let requested_planet = planet;
    let raw_flags = prepared.map_or_else(
        || consume_non_echoed_flags(flags, planet),
        |flags| flags.raw,
    );
    let norm_flags = prepared.map_or_else(
        || normalize_calc_flags(raw_flags, engine.jpl_fallback_active),
        |flags| flags.normalized,
    );
    let planet = remap_ast_offset(planet);
    if planet == -MEAN_NODE || planet == -TRUE_NODE {
        // Publish the descending-node trace while restoring its ascending-node record.
        let north = planet.abs();
        let previous = crate::trace::try_snapshot_record(north)?;
        let (north_result, retflag) = match calc(engine, tjdet, north, norm_flags) {
            Ok(result) => result,
            Err(error) => {
                crate::trace::with_caught_error(&error, || {
                    crate::trace::try_restore_record(north, previous)
                })?;
                return Err(error);
            }
        };
        crate::trace::try_record_alias(planet, north, previous)?;
        return Ok((
            south_node_from_north(north_result, norm_flags),
            if prepared.is_some_and(|flags| flags.defer_echo) {
                retflag
            } else {
                calc_tt_epheflag_echo(echo_request_bits(retflag, raw_flags), raw_flags)
            },
        ));
    }
    super::routing::with_native_trace_body(requested_planet, planet, || {
        let defer_echo = prepared.is_some_and(|flags| flags.defer_echo);
        let (pos, retflag) = if defer_echo {
            calc_tt_body_source_time_with_echo(engine, tjdet, planet, raw_flags, norm_flags, true)?
        } else {
            calc_tt_body(engine, tjdet, planet, raw_flags, norm_flags)?
        };
        let result = (
            pos,
            if defer_echo {
                retflag
            } else {
                calc_tt_epheflag_echo(retflag, raw_flags)
            },
        );
        super::routing::publish_native_trace_source(requested_planet)?;
        Ok(result)
    })
}

/// TT position entry with an optional direct Cartesian computation.
///
/// Normalization, ephemeris echo, ECL_NUT, lunar-node recursion and mean-point
/// SPEED3 keep the [`calc`] contract. For non-sidereal XYZ requests the
/// additional provider returns the six Cartesian AU/AU-day components
/// without a spherical roundtrip. `None` selects the existing spherical
/// provider and its output conversion, as required by non-Cartesian channels.
/// XYZ takes precedence over RADIANS, while both requested bits are echoed.
/// Errors and model/data requirements propagate unchanged.
pub fn calc_with_cartesian(
    engine: &CalcEngine<'_>,
    cartesian: &CartesianCompute<'_>,
    jd_tt: f64,
    body_id: i32,
    flags: i32,
) -> Result<(PositionResult, i32)> {
    calc_with_cartesian_entry(engine, cartesian, jd_tt, body_id, flags, false, None)
}

/// UT1 position entry with direct Cartesian output where available.
///
/// The time provider and flag/error order are those of [`calc_ut`]; output
/// routing follows [`calc_with_cartesian`]. The direct provider receives the
/// acquired TT epoch, never UT1. Existing spherical callbacks are unchanged.
pub fn calc_ut_with_cartesian(
    engine: &CalcEngine<'_>,
    cartesian: &CartesianCompute<'_>,
    jd_ut1: f64,
    body_id: i32,
    flags: i32,
) -> Result<(PositionResult, i32)> {
    calc_with_cartesian_entry(engine, cartesian, jd_ut1, body_id, flags, true, None)
}

/// Adapt the legacy spherical entry without changing its callback contract.
///
/// At TT epochs, substitute direct Cartesian components (or the spherical
/// channel's inverse-Jacobian output) as the entry's unformatted six slots.
/// Remove only format flags during that entry, then restore their echo.
/// Node antipodes, mean-point SPEED3, ECL_NUT and sidereal projections retain
/// their original spherical/special-body dispatch and shaping instead.
/// Direct-Cartesian numerical adapter with one already completed flag prelude.
/// UT1 versus TT selects only the original clock; all provider/format/SPEED3
/// decisions retain the same adapter body as the legacy entry.
pub fn calc_with_cartesian_prepared(
    engine: &CalcEngine<'_>,
    cartesian: &CartesianCompute<'_>,
    jd: f64,
    body: i32,
    flags: PreparedCalcFlags,
    is_ut: bool,
) -> Result<(PositionResult, i32)> {
    calc_with_cartesian_entry(
        engine,
        cartesian,
        jd,
        body,
        flags.normalized,
        is_ut,
        Some(flags),
    )
}

fn calc_with_cartesian_entry(
    engine: &CalcEngine<'_>,
    cartesian: &CartesianCompute<'_>,
    jd: f64,
    body_id: i32,
    flags: i32,
    is_ut1: bool,
    prepared: Option<PreparedCalcFlags>,
) -> Result<(PositionResult, i32)> {
    let entry = |provider: &CalcEngine<'_>, request: i32| {
        if let Some(mut flags) = prepared {
            // The Cartesian adapter masks only its own output format. Keep the
            // single-normalized SPEED3 marker and post-consumption raw echo.
            let removed = flags.normalized & !request;
            flags.raw &= !removed;
            flags.normalized = request;
            if is_ut1 {
                calc_ut_prepared(provider, jd, body_id, flags)
            } else {
                calc_prepared(provider, jd, body_id, flags)
            }
        } else if is_ut1 {
            calc_ut(provider, jd, body_id, request)
        } else {
            calc(provider, jd, body_id, request)
        }
    };
    let body = remap_ast_offset(body_id);
    if flags & FLG_XYZ == 0
        || flags & FLG_SIDEREAL != 0
        || body == ECL_NUT
        || matches!(body, MEAN_NODE | TRUE_NODE | MEAN_APOG)
        || body == -MEAN_NODE
        || body == -TRUE_NODE
    {
        return entry(engine, flags);
    }
    let compute_tt = |tt, id, calc_flags| {
        if let Some(value) = cartesian(tt, id, calc_flags)? {
            Ok(value)
        } else {
            let (position, retflag) = (engine.compute_tt)(tt, id, calc_flags)?;
            Ok((apply_output_flags(position, FLG_XYZ), retflag))
        }
    };
    let adapted = CalcEngine {
        compute_tt: &compute_tt,
        ut1_to_tt_jd: engine.ut1_to_tt_jd,
        nutation_obliquity_ut: engine.nutation_obliquity_ut,
        nutation_obliquity_tt: engine.nutation_obliquity_tt,
        jpl_fallback_active: engine.jpl_fallback_active,
    };
    let format = flags & (FLG_XYZ | FLG_RADIANS);
    let (position, retflag) = entry(&adapted, flags & !format)?;
    Ok((position, retflag | format))
}

/// Look up the sidereal offset or refuse (entry-layer contract): sidereal
/// planet-centric output needs the sidereal provider; inventing a zero
/// offset would silently return tropical positions.
pub(super) fn pctr_ayanamsa_deg(engine: &PctrEngine<'_>, jd_tt: f64) -> Result<f64> {
    match engine.ayanamsha_deg {
        Some(provider) => provider(jd_tt),
        None => Err(Error::ConfigurationMissing {
            missing_config: "sidereal_ayanamsha",
            suggestion: "Provide an ayanamsha source (SID-03) for sidereal planet-centric output.",
        }),
    }
}

#[cfg(test)]
mod context_pctr_parts_tests {
    use super::*;
    use std::cell::RefCell;

    #[test]
    fn central_parts_and_scalar_tt_side_parts_follow_original_alias_order() {
        let reads = RefCell::new(Vec::new());
        let aliases = || {
            reads.borrow_mut().push(("alias", MOON, 0, 0));
            reads.borrow_mut().push(("alias", EARTH, 0, 0));
            Ok(())
        };
        let state = |body, epoch: TwoPartEpoch| {
            reads
                .borrow_mut()
                .push(("state", body, epoch.hi.to_bits(), epoch.lo.to_bits()));
            Ok((
                if body == MOON {
                    [2.0, 0.4, 0.1]
                } else {
                    [1.0, 0.0, 0.0]
                },
                [0.0; 3],
            )
                .into())
        };
        let epoch = TwoPartEpoch {
            hi: 2451545.125,
            lo: 1.0e-10,
        };
        let flags = FLG_TRUEPOS | FLG_SPEED | FLG_J2000 | FLG_NONUT;
        let result = pctr_vector_result_acquired(
            &state,
            None,
            &|_, _| panic!("nonsidereal"),
            epoch,
            MOON,
            EARTH,
            flags,
            Some(&aliases),
            PctrRatePolicy::Sampled,
        )
        .unwrap();
        assert!(result.values.iter().all(|word| word.is_finite()));
        let half = MOON_SPEED_HALF_STEP_DAYS;
        let expected = [
            epoch,
            pctr_tt_parts(epoch.value(), -half),
            pctr_tt_parts(epoch.value(), half),
        ];
        let reads = reads.borrow();
        assert_eq!(reads.len(), 12);
        for (sample, epoch) in reads.chunks_exact(4).zip(expected) {
            assert_eq!(sample[0], ("alias", MOON, 0, 0));
            assert_eq!(sample[1], ("alias", EARTH, 0, 0));
            assert_eq!(
                sample[2],
                ("state", MOON, epoch.hi.to_bits(), epoch.lo.to_bits())
            );
            assert_eq!(
                sample[3],
                ("state", EARTH, epoch.hi.to_bits(), epoch.lo.to_bits())
            );
        }
        assert_ne!(epoch.lo.to_bits(), ((epoch.value()) - epoch.hi).to_bits());
    }

    #[test]
    fn source_rereads_observer_before_three_retarded_targets_without_center_mapping() {
        for center_flag in [FLG_HELCTR, FLG_BARYCTR] {
            let reads = RefCell::new(Vec::new());
            let aliases = || {
                reads.borrow_mut().push(("aliases", MARS));
                Ok(())
            };
            let state = |body, _: TwoPartEpoch| {
                reads.borrow_mut().push(("tt", body));
                Ok((
                    if body == MARS {
                        [2.0, 0.4, 0.1]
                    } else {
                        [1.0, 0.0, 0.0]
                    },
                    [0.0; 3],
                )
                    .into())
            };
            let retarded = |body, epoch: TwoPartEpoch| {
                assert_eq!(epoch.hi.to_bits(), 2451545.125f64.to_bits());
                reads.borrow_mut().push(("tdb", body));
                Ok(([2.0, 0.4, 0.1], [0.0; 3]).into())
            };
            let flags = center_flag | FLG_NOGDEFL | FLG_NOABERR | FLG_J2000 | FLG_NONUT;
            pctr_vector_result_acquired(
                &state,
                Some(&retarded),
                &|_, _| panic!("nonsidereal"),
                TwoPartEpoch {
                    hi: 2451545.125,
                    lo: 1.0e-10,
                },
                MARS,
                EARTH,
                flags,
                Some(&aliases),
                PctrRatePolicy::Sampled,
            )
            .unwrap();
            assert_eq!(
                *reads.borrow(),
                [
                    ("aliases", MARS),
                    ("tt", MARS),
                    ("tt", EARTH),
                    ("tt", EARTH),
                    ("tdb", MARS),
                    ("tdb", MARS),
                    ("tdb", MARS),
                ]
            );
        }
    }

    #[test]
    fn alias_or_target_failure_stops_before_time_and_observer_work() {
        let marker = Error::UnknownBody { body_id: 719 };
        let alias = || Err(marker.clone());
        let result = pctr_vector_result_acquired(
            &|_, _| panic!("alias failed before state"),
            None,
            &|_, _| panic!("unused offset"),
            TwoPartEpoch {
                hi: f64::NAN,
                lo: f64::NAN,
            },
            MARS,
            EARTH,
            0,
            Some(&alias),
            PctrRatePolicy::Sampled,
        );
        assert_eq!(result, Err(marker.clone()));
        let reads = RefCell::new(Vec::new());
        let state = |body, _: TwoPartEpoch| {
            reads.borrow_mut().push(body);
            Err(marker.clone())
        };
        let result = pctr_vector_result_acquired(
            &state,
            None,
            &|_, _| panic!("unused offset"),
            TwoPartEpoch {
                hi: 2451545.0,
                lo: 0.001,
            },
            MARS,
            EARTH,
            0,
            Some(&|| Ok(())),
            PctrRatePolicy::Sampled,
        );
        assert_eq!(result, Err(marker));
        assert_eq!(*reads.borrow(), [MARS]);
    }

    #[test]
    fn module_tt_delegate_keeps_single_observer_and_original_parts() {
        let reads = RefCell::new(Vec::new());
        let state = |body, epoch: TwoPartEpoch| {
            reads
                .borrow_mut()
                .push((body, epoch.hi.to_bits(), epoch.lo.to_bits()));
            Ok((
                if body == MARS {
                    [2.0, 0.4, 0.1]
                } else {
                    [1.0, 0.0, 0.0]
                },
                [0.0; 3],
            )
                .into())
        };
        let tt = 2451545.125;
        pctr_vector_result(
            &state,
            Some(&state),
            &|_, _| panic!("unused offset"),
            tt,
            MARS,
            EARTH,
            FLG_NOGDEFL | FLG_NOABERR | FLG_J2000 | FLG_NONUT,
            PctrRatePolicy::Sampled,
        )
        .unwrap();
        let reads = reads.borrow();
        assert_eq!(
            &reads[..2],
            [
                (MARS, 2451545.0f64.to_bits(), 0.125f64.to_bits()),
                (EARTH, 2451545.0f64.to_bits(), 0.125f64.to_bits())
            ]
        );
        assert_eq!(reads.iter().filter(|read| read.0 == EARTH).count(), 1);
        assert_eq!(reads.len(), 5);
    }
}

/// Batch `calc_ut` over one body/flag set and many UT1 epochs (batch-entry layer).
///
/// One scalar [`calc_ut`] call per epoch, in input order, failing fast
/// with the first element's error. No segment, stencil, or flag path of
/// its own — element `i` is bit-identical to a repeated scalar call
/// with `(tjds[i], body_id, flags)`.
pub fn calc_ut_many(
    engine: &CalcEngine<'_>,
    tjds_ut: &[f64],
    body_id: i32,
    flags: i32,
) -> Result<Vec<(PositionResult, i32)>> {
    let mut out = Vec::with_capacity(tjds_ut.len());
    for &jd in tjds_ut {
        out.push(calc_ut(engine, jd, body_id, flags)?);
    }
    Ok(out)
}

/// Batch `calc` over one body/flag set and many TT epochs (batch-entry layer).
///
/// Same fail-fast/input-order contract as [`calc_ut_many`]; element `i`
/// equals a repeated scalar [`calc`] call with
/// `(tjds_tt[i], body_id, flags)`.
pub fn calc_many(
    engine: &CalcEngine<'_>,
    tjds_tt: &[f64],
    body_id: i32,
    flags: i32,
) -> Result<Vec<(PositionResult, i32)>> {
    let mut out = Vec::with_capacity(tjds_tt.len());
    for &jd in tjds_tt {
        out.push(calc(engine, jd, body_id, flags)?);
    }
    Ok(out)
}

/// Heterogeneous batch `calc_ut` over per-request epochs, bodies, and
/// flags (batch-entry layer).
///
/// Each element runs the full scalar [`calc_ut`] entry — per-request
/// flag normalization/echo, asteroid remap, south-node recursion, and
/// the `ECL_NUT` branch — so a mixed batch agrees element-wise with
/// repeated scalar calls. Fails fast with the first element's error;
/// an empty batch returns an empty vector.
pub fn calc_ut_requests(
    engine: &CalcEngine<'_>,
    requests: &[CalcRequest],
) -> Result<Vec<(PositionResult, i32)>> {
    let mut out = Vec::with_capacity(requests.len());
    for request in requests {
        out.push(calc_ut(engine, request.jd, request.body_id, request.flags)?);
    }
    Ok(out)
}

/// Heterogeneous batch `calc` over per-request TT epochs, bodies, and
/// flags (batch-entry layer).
///
/// TT analogue of [`calc_ut_requests`]: element `i` equals a repeated
/// scalar [`calc`] call with the request's `(jd, body_id, flags)`,
/// including the TT ephemeris-bit echo convention.
pub fn calc_requests(
    engine: &CalcEngine<'_>,
    requests: &[CalcRequest],
) -> Result<Vec<(PositionResult, i32)>> {
    let mut out = Vec::with_capacity(requests.len());
    for request in requests {
        out.push(calc(engine, request.jd, request.body_id, request.flags)?);
    }
    Ok(out)
}

/// Batch [`calc_pctr`] over one target/center/flag set and many TT
/// epochs (batch-entry layer).
///
/// One scalar call per epoch, in input order, failing fast with the
/// first element's error.
pub fn calc_pctr_many(
    engine: &PctrEngine<'_>,
    tjds_tt: &[f64],
    target: i32,
    center: i32,
    flags: i32,
) -> Result<Vec<(PositionResult, i32)>> {
    let mut out = Vec::with_capacity(tjds_tt.len());
    for &jd in tjds_tt {
        out.push(calc_pctr(engine, jd, target, center, flags)?);
    }
    Ok(out)
}

/// Heterogeneous batch [`calc_pctr`] over per-request epochs,
/// endpoints, and flags (batch-entry layer).
///
/// Element `i` equals a repeated scalar [`calc_pctr`] call with the
/// request's `(jd_tt, target, center, flags)`; fails fast with the
/// first element's error.
pub fn calc_pctr_requests(
    engine: &PctrEngine<'_>,
    requests: &[PctrRequest],
) -> Result<Vec<(PositionResult, i32)>> {
    let mut out = Vec::with_capacity(requests.len());
    for request in requests {
        out.push(calc_pctr(
            engine,
            request.jd_tt,
            request.target,
            request.center,
            request.flags,
        )?);
    }
    Ok(out)
}

/// Resolve numbered aliases for Ceres, Pallas, Juno, Vesta, Chiron and Pholus
/// to their dedicated body identifiers; all other identifiers pass unchanged.
pub fn remap_ast_offset(body_id: i32) -> i32 {
    if body_id > AST_OFFSET {
        match body_id - AST_OFFSET {
            1 => CERES,
            2 => PALLAS,
            3 => JUNO,
            4 => VESTA,
            2060 => CHIRON,
            5145 => PHOLUS,
            _ => body_id,
        }
    } else {
        body_id
    }
}

/// Whether the body is a standard Sun/Moon/planet/Earth ephemeris target.
/// Analytical lunar, hypothetical and element-based minor models are classified
/// separately. Actual epoch coverage remains the selected provider's decision.
pub fn body_uses_jpl_ephemeris(body_id: i32) -> bool {
    matches!(
        body_id,
        SUN | MOON | MERCURY | VENUS | MARS | JUPITER | SATURN | URANUS | NEPTUNE | PLUTO | EARTH
    )
}

/// Whether an identifier names a legacy satellite alias or a recognized
/// NAIF satellite/planet-center ID. The authoritative mapping is
/// [`crate::planetary_moons::moon_naif_id`]; this predicate only classifies.
pub fn is_planetary_moon_id(body_id: i32) -> bool {
    crate::planetary_moons::moon_naif_id(body_id).is_some()
}

/// Standard body name, including the six built-in numbered asteroid aliases.
/// Unrecognized identifiers return an empty string. Registry and satellite
/// name tables are resolved by their respective public family APIs.
pub fn planet_name(body_id: i32) -> &'static str {
    match body_id {
        SUN => "Sun",
        MOON => "Moon",
        MERCURY => "Mercury",
        VENUS => "Venus",
        MARS => "Mars",
        JUPITER => "Jupiter",
        SATURN => "Saturn",
        URANUS => "Uranus",
        NEPTUNE => "Neptune",
        PLUTO => "Pluto",
        MEAN_NODE => "mean Node",
        TRUE_NODE => "true Node",
        MEAN_APOG => "mean Apogee",
        OSCU_APOG => "osc. Apogee",
        INTP_APOG => "intp. Apogee",
        INTP_PERG => "intp. Perigee",
        EARTH => "Earth",
        ISIS => "Isis-Transpluto",
        NIBIRU => "Nibiru",
        HARRINGTON => "Harrington",
        NEPTUNE_LEVERRIER => "Leverrier (Neptune)",
        NEPTUNE_ADAMS => "Adams (Neptune)",
        PLUTO_LOWELL => "Lowell (Pluto)",
        PLUTO_PICKERING => "Pickering (Pluto)",
        CUPIDO => "Cupido",
        HADES => "Hades",
        ZEUS => "Zeus",
        KRONOS => "Kronos",
        APOLLON => "Apollon",
        ADMETOS => "Admetos",
        VULKANUS => "Vulcanus",
        POSEIDON => "Poseidon",
        VULCAN => "Vulcan",
        WHITE_MOON => "Selena/White Moon",
        PROSERPINA => "Proserpina",
        WALDEMATH => "Waldemath",
        CHIRON => "Chiron",
        PHOLUS => "Pholus",
        CERES => "Ceres",
        PALLAS => "Pallas",
        JUNO => "Juno",
        VESTA => "Vesta",
        _ => {
            if body_id > AST_OFFSET {
                match body_id - AST_OFFSET {
                    1 => "Ceres",
                    2 => "Pallas",
                    3 => "Juno",
                    4 => "Vesta",
                    2060 => "Chiron",
                    5145 => "Pholus",
                    _ => "",
                }
            } else {
                ""
            }
        }
    }
}

/// Owning task of a serving family, for gap tracking.
pub fn body_kind_owner(kind: BodyKind) -> &'static str {
    match kind {
        BodyKind::Planet => "POS-02..09 provider path (served)",
        BodyKind::LunarNode | BodyKind::LunarApse => "LUN-02",
        BodyKind::InterpolatedApse => "LUN-03",
        BodyKind::Uranian | BodyKind::Transpluto | BodyKind::FictHelio => "BODY-03",
        BodyKind::PlanetaryMoon => "BODY-05",
        BodyKind::MinorBody => "BODY-04",
        BodyKind::FixedStar => "STAR-02..05",
        BodyKind::Angle => "HOU-07",
        BodyKind::ArabicPart => "BODY-06",
        BodyKind::Unknown => "none (unknown body)",
    }
}

/// Classify a body after numbered-alias and descending-node resolution.
/// Satellite aliases precede lunar points, hypothetical models and minor
/// bodies; standard planets and unsupported identifiers complete the selection.
pub fn classify_body(body_id: i32) -> BodyKind {
    if is_planetary_moon_id(body_id) {
        return BodyKind::PlanetaryMoon;
    }
    if body_id == MEAN_NODE || body_id == TRUE_NODE {
        return BodyKind::LunarNode;
    }
    if body_id == MEAN_APOG || body_id == OSCU_APOG {
        return BodyKind::LunarApse;
    }
    if body_id == INTP_APOG || body_id == INTP_PERG {
        return BodyKind::InterpolatedApse;
    }
    // The eight Uranian model bodies occupy one contiguous identifier interval.
    if (CUPIDO..=POSEIDON).contains(&body_id) {
        return BodyKind::Uranian;
    }
    if body_id == ISIS {
        return BodyKind::Transpluto;
    }
    if body_id == WHITE_MOON
        || body_id == VULCAN
        || body_id == PROSERPINA
        || body_id == WALDEMATH
        || body_id == NIBIRU
        || body_id == HARRINGTON
        || body_id == NEPTUNE_LEVERRIER
        || body_id == NEPTUNE_ADAMS
        || body_id == PLUTO_LOWELL
        || body_id == PLUTO_PICKERING
    {
        return BodyKind::FictHelio;
    }
    if matches!(body_id, CHIRON | PHOLUS | CERES | PALLAS | JUNO | VESTA)
        || (AST_OFFSET < body_id && body_id < FIXSTAR_OFFSET)
    {
        return BodyKind::MinorBody;
    }
    if FIXSTAR_OFFSET < body_id && body_id < ANGLE_OFFSET {
        return BodyKind::FixedStar;
    }
    if (ANGLE_OFFSET..ARABIC_OFFSET).contains(&body_id) {
        return BodyKind::Angle;
    }
    if (ARABIC_OFFSET..ARABIC_OFFSET + 100).contains(&body_id) {
        return BodyKind::ArabicPart;
    }
    if body_uses_jpl_ephemeris(body_id) {
        return BodyKind::Planet;
    }
    BodyKind::Unknown
}

/// Validate a body channel for stored ICRS apparent-place reduction.
///
/// Analytical hypothetical models are not persisted-channel targets. Absent
/// bodies produce a typed missing-channel error; non-ICRS channels produce
/// [`Error::SealedChannelNotServed`]. The returned [`StateFrame`] determines
/// whether physical center-minus-barycenter composition is required.
pub fn serving_channel_gate<P: StateProvider>(provider: &P, body_id: i32) -> Result<StateFrame> {
    if (FICT_OFFSET..=WALDEMATH).contains(&body_id) {
        return Err(Error::UnknownBody { body_id });
    }
    if !provider.has_body(body_id) {
        return Err(Error::LebBodyNotFound { body_id });
    }
    let coord_type = provider
        .serving_coord_type(body_id)
        .ok_or(Error::LebBodyNotFound { body_id })?;
    let frame = StateFrame::from_coord_type(coord_type)?;
    match frame {
        StateFrame::IcrsBarycentric | StateFrame::IcrsSystemBarycenter => Ok(frame),
        StateFrame::EclipticDirect
        | StateFrame::GeocentricEcliptic
        | StateFrame::HeliocentricEclipticRetired => Err(Error::SealedChannelNotServed {
            body_id,
            coord_type,
        }),
    }
}

impl<'a, P> SealedServing<'a, P> {
    /// Borrow motion from one trusted native scalar call, with rollback on every exit.
    pub(crate) fn with_output_motion<R>(
        &self,
        tt: f64,
        body: i32,
        flags: i32,
        calculate: impl FnOnce() -> Result<R>,
    ) -> Result<(R, Option<super::derivative::OutputWordMotion>)> {
        self.output_motion.with_motion(
            OutputMotionQuery::new(tt, None, body, flags, false),
            calculate,
        )
    }

    /// Borrow the identical source/feed with its complete previous sampled law.
    /// An opaque external anchor has no derivative capability; changing only
    /// the tropical portion would introduce unqualified cancellation changes.
    pub(crate) fn with_sampled_derivatives(&self) -> Self {
        Self {
            provider: self.provider,
            nutation_at: self.nutation_at,
            stored_nutation_feed: self.stored_nutation_feed,
            reuse_immutable_central: self.reuse_immutable_central,
            force_sampled_derivatives: true,
            owned_apparent_id: None,
            output_motion: OutputMotionCapture::default(),
        }
    }

    /// Checked constructor: a provider/nutation pair serving one
    /// Pipeline-A computation. No I/O, no globals.
    pub fn new(provider: &'a P, nutation_at: &'a dyn Fn(f64) -> Result<(f64, f64)>) -> Self {
        Self {
            provider,
            nutation_at,
            reuse_immutable_central: false,
            stored_nutation_feed: false,
            output_motion: OutputMotionCapture::default(),
            force_sampled_derivatives: false,
            owned_apparent_id: None,
        }
    }

    /// Bind a producer-owned immutable reader and its actual stored nutation.
    /// Callers must own the reader snapshot for the complete borrowed serving;
    /// custom providers or replacement feeds use the public constructor instead.
    pub(super) fn new_stored(
        provider: &'a P,
        nutation_at: &'a dyn Fn(f64) -> Result<(f64, f64)>,
    ) -> Self {
        Self {
            provider,
            nutation_at,
            reuse_immutable_central: true,
            stored_nutation_feed: true,
            output_motion: OutputMotionCapture::default(),
            force_sampled_derivatives: false,
            owned_apparent_id: None,
        }
    }

    /// Lend a producer token only at the native immutable Owned Session entry.
    /// This fresh invocation owns no inherited output-motion capture scope.
    pub(crate) fn with_owned_apparent_token(
        &self,
        token: Option<std::num::NonZeroU64>,
    ) -> SealedServing<'_, P> {
        SealedServing {
            provider: self.provider,
            nutation_at: self.nutation_at,
            reuse_immutable_central: self.reuse_immutable_central,
            stored_nutation_feed: self.stored_nutation_feed,
            output_motion: OutputMotionCapture::default(),
            force_sampled_derivatives: self.force_sampled_derivatives,
            owned_apparent_id: token,
        }
    }

    /// Serving geometric-state provider (shared borrow; the
    /// pre-existing by-value form was unusable through a shared
    /// reference, which is how every serving binding holds it).
    pub fn provider(&self) -> &P
    where
        P: StateProvider,
    {
        self.provider
    }

    /// Nutation-angle provider (radians at a TT epoch).
    pub fn nutation_at(self) -> &'a dyn Fn(f64) -> Result<(f64, f64)> {
        self.nutation_at
    }

    /// Shared-borrow view of the nutation-angle provider (the by-value
    /// `nutation_at` form is unusable through the shared references the
    /// serving bindings hold).
    pub fn nutation_provider(&self) -> &dyn Fn(f64) -> Result<(f64, f64)> {
        self.nutation_at
    }
}

impl<P: StateProvider> SealedServing<'_, P> {
    /// Physical vector geometry at an original observation TT, in ecliptic
    /// degrees and AU with degree/day and AU/day rates. The observing center
    /// is evaluated before the target, while scalar stored-position entries
    /// retain target-first acquisition. `obs_split` keeps the observation TT
    /// parts through target retardation; provider and frame failures propagate.
    pub(crate) fn physical_vector_at_tt(
        &self,
        jd_tt: f64,
        body_id: i32,
        calc_flags: i32,
        obs_split: (f64, f64),
    ) -> Result<(PositionResult, i32)> {
        serving_position_core_ordered(
            self,
            jd_tt,
            body_id,
            calc_flags,
            false,
            Some(obs_split),
            CentralStateOrder::ObserverFirst,
        )
    }

    /// Physical surface-observer geometry with independent TT/UT1 clocks.
    /// Earth is evaluated before the surface rotation and target. The retained
    /// IERS2010 site, Earth-mass deflection and split observation TT use the
    /// existing vector numerical reduction. Output units match
    /// [`Self::physical_vector_at_tt`]; site, provider and frame errors retain
    /// their consumed order without a scalar-position retry.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn physical_vector_topocentric_at_tt(
        &self,
        jd_tt: f64,
        jd_ut1: f64,
        body_id: i32,
        calc_flags: i32,
        geopos: (f64, f64, f64),
        c2t_at: &dyn Fn(f64, f64) -> Result<[[f64; 3]; 3]>,
        obs_split: (f64, f64),
    ) -> Result<(PositionResult, i32)> {
        serving_topocentric_core_ordered(
            self,
            jd_tt,
            jd_ut1,
            body_id,
            calc_flags,
            geopos,
            c2t_at,
            false,
            true,
            SiteEllipsoid::Iers2010,
            Some(obs_split),
            CentralStateOrder::ObserverFirst,
            None,
        )
    }
}

/// Central numerical acquisition order, independent of epoch representation.
#[derive(Clone, Copy, PartialEq, Eq)]
enum CentralStateOrder {
    TargetFirst,
    ObserverFirst,
}

/// Compute an apparent place from stored geometric states at a TT Julian Day.
///
/// The provider supplies actual body-center geometry and nutation. The shared
/// pipeline composes the observer, retards the target, applies enabled PPN
/// and Lorentz corrections and rotates into the requested frame. Rates use
/// the physical light-path, correction and frame derivatives. The result is
/// an unshaped spherical degrees/degrees/AU core with per-TT-day rates, plus
/// calculation flags. The public entry applies representation shaping.
/// Missing channels, coverage, invalid states and correction errors propagate.
pub fn serving_compute_tt<P: StateProvider>(
    serving: &SealedServing<'_, P>,
    jd_tt: f64,
    body_id: i32,
    calc_flags: i32,
) -> Result<(PositionResult, i32)> {
    serving_position_with_vector_fallback(serving, jd_tt, body_id, calc_flags, false, None)
}

/// Compute an apparent place from stored geometric states at a TT Julian Day.
///
/// The provider supplies actual body-center geometry and nutation. The shared
/// pipeline composes the observer, retards the target, applies enabled PPN
/// and Lorentz corrections and rotates into the requested frame. Rates use
/// the physical light-path, correction and frame derivatives. The result is
/// an unshaped spherical degrees/degrees/AU core with per-TT-day rates, plus
/// calculation flags. The public entry applies representation shaping.
/// Missing channels, coverage, invalid states and correction errors propagate.
pub fn serving_compute_tt_split<P: StateProvider>(
    serving: &SealedServing<'_, P>,
    jd_tt: f64,
    body_id: i32,
    calc_flags: i32,
    obs_split: Option<(f64, f64)>,
) -> Result<(PositionResult, i32)> {
    serving_position_with_vector_fallback(serving, jd_tt, body_id, calc_flags, false, obs_split)
}

/// Spherical stored-state computation without a derived-range retry.
///
/// Inputs and degree/AU/day outputs match [`serving_compute_tt_split`].
/// `obs_split` retains optional observation TT parts for light time. Every
/// support, retarded and speed-stencil read propagates its typed failure.
/// This lets an enclosing multi-leg calculation own its complete fallback
/// and clock selection instead of retrying an individual position leg.
pub fn serving_compute_tt_direct<P: StateProvider>(
    serving: &SealedServing<'_, P>,
    jd_tt: f64,
    body_id: i32,
    calc_flags: i32,
    obs_split: Option<(f64, f64)>,
) -> Result<(PositionResult, i32)> {
    serving_position_core(serving, jd_tt, body_id, calc_flags, false, obs_split)
}

/// Apparent Cartesian position and velocity at TT.
///
/// Returns `(x, y, z, vx, vy, vz)` in AU and AU/day in the frame selected
/// by `calc_flags`, plus the computation retflag. With SPEED cleared the
/// last three components are zero. Eligible native kinematics differentiate
/// the solved apparent state and output frame analytically; other requests
/// differentiate the actual sampled places on their source domains. Cartesian
/// projection avoids spherical conversion, angular wrapping and inverse-Jacobian
/// reconstruction in either route. Geometric acquisition, corrections, model
/// validity and errors follow [`serving_compute_tt`].
pub fn serving_compute_xyz_tt<P: StateProvider>(
    serving: &SealedServing<'_, P>,
    jd_tt: f64,
    body_id: i32,
    calc_flags: i32,
) -> Result<(PositionResult, i32)> {
    serving_position_with_vector_fallback(serving, jd_tt, body_id, calc_flags, true, None)
}

/// Apparent Cartesian position in AU at a TT Julian date, without speeds.
///
/// Uses the same geometric states, light-time iteration, deflection,
/// aberration, and output-frame selection as [`serving_compute_tt`], but
/// returns the rotated vector directly. No spherical conversion or inverse
/// conversion is performed. Equatorial plus J2000 flags select the biased
/// J2000 mean equator; adding ICRS removes that bias. The speed flag is ignored.
///
/// This geocentric/heliocentric/barycentric path requires ICRS state channels.
/// Non-finite epochs, unavailable channels, coverage misses and topocentric
/// requests without an observer return the same typed errors as the spherical
/// path. Distances retain the provider's physical model and accuracy limits.
/// Nutation data are still required by the shared frame assembly.
pub fn serving_cartesian_tt<P: StateProvider>(
    serving: &SealedServing<'_, P>,
    jd_tt: f64,
    body_id: i32,
    calc_flags: i32,
) -> Result<[f64; 3]> {
    let (position, _) =
        serving_position_core(serving, jd_tt, body_id, calc_flags & !FLG_SPEED, true, None)?;
    Ok([position[0], position[1], position[2]])
}

/// Try the stored-state fast pipeline, retaining the full vector fallback.
///
/// A derived coverage miss can leave the coarse fast SPEED stencil outside
/// coverage even when a derivative within the central owned polynomial remains
/// available. Only a typed stored-range miss triggers this
/// alternative; central coverage and every fallback read remain checked. No
/// polynomial extrapolation, new tolerance or corruption fallback is admitted.
/// Sample clocks retain split TT, but each LEB-vector state query collapses
/// `whole + fraction` at its reader boundary. Vector frames and apparent
/// corrections follow the complete vector observation. Output units/flags
/// match the fast core, including Cartesian public shaping when requested.
fn serving_position_with_vector_fallback<P: StateProvider>(
    serving: &SealedServing<'_, P>,
    jd_tt: f64,
    body_id: i32,
    calc_flags: i32,
    want_xyz: bool,
    obs_split: Option<(f64, f64)>,
) -> Result<(PositionResult, i32)> {
    match serving_position_core(serving, jd_tt, body_id, calc_flags, want_xyz, obs_split) {
        Err(Error::LebOutOfRange { .. })
            if vector_fallback_central_covered(serving.provider, body_id, jd_tt) =>
        {
            serving_vector_state_fallback(serving, jd_tt, body_id, calc_flags, None, want_xyz)
        }
        outcome => outcome,
    }
}

/// Whether stored metadata admits a central epoch for the full vector route.
///
/// This guard excludes central-outside and unknown-coverage requests, retaining
/// the fast core's original typed error and acquisition order for those calls.
/// Bounds are closed; no epsilon or extrapolation is applied.
fn vector_fallback_central_covered<P: StateProvider>(provider: &P, body: i32, jd: f64) -> bool {
    provider
        .body_coverage(body)
        .is_some_and(|(start, end)| start <= jd && jd <= end)
}

/// Full LEB-vector position after a fast-path derived-range miss.
///
/// Central and every retarded/stencil body state remain range checked. Angular
/// samples remain on the actual central geometric polynomials. Each LEB read
/// consumes scalar TT. Optional sites use the IERS2010 geoid and independent
/// bundled UT1/frame clock; the fast c2t/model partner is not reused here.
/// Native series supply full-vector nutation; no stored polynomial is extended.
/// Returns spherical components/rates or public Cartesian shaping, plus flags.
fn serving_vector_state_fallback<P: StateProvider>(
    serving: &SealedServing<'_, P>,
    jd_tt: f64,
    body_id: i32,
    calc_flags: i32,
    site: Option<(f64, f64, f64)>,
    want_xyz: bool,
) -> Result<(PositionResult, i32)> {
    // The public sealed guard retains central stored-body coverage, including
    // self-observation; a zero result does not permit outside central reads.
    let _ = serving.provider.geometric_state(body_id, jd_tt)?;
    let mut position = super::stored_vector_retry::stored_vector_position(
        serving.provider,
        &|error| error,
        jd_tt,
        None,
        body_id,
        calc_flags,
        site,
        false,
    )?;
    if want_xyz {
        position = apply_output_flags(position, FLG_XYZ);
    }
    super::routing::record_stored_trace_source(body_id)?;
    Ok((position, calc_flags))
}

/// Shared [`serving_compute_tt`] apparent-position pipeline. `want_xyz`
/// selects direct AU/AU-day output instead of spherical degrees/degrees/AU
/// and their rates. The same selection applies to the central value and
/// the analytic or sampled rate route. Sampled spherical rates retain their
/// angular-wrap convention; Cartesian projection does not wrap angles.
pub(super) fn serving_position_core<P: StateProvider>(
    serving: &SealedServing<'_, P>,
    jd_tt: f64,
    body_id: i32,
    calc_flags: i32,
    want_xyz: bool,
    obs_split: Option<(f64, f64)>,
) -> Result<(PositionResult, i32)> {
    serving_position_core_ordered(
        serving,
        jd_tt,
        body_id,
        calc_flags,
        want_xyz,
        obs_split,
        CentralStateOrder::TargetFirst,
    )
}

/// Shared numerical reduction after choosing the scalar or vector central
/// acquisition order. Changing the order does not change support selection,
/// target retardation, frame arithmetic or velocity formulas.
fn serving_position_core_ordered<P: StateProvider>(
    serving: &SealedServing<'_, P>,
    jd_tt: f64,
    body_id: i32,
    calc_flags: i32,
    want_xyz: bool,
    obs_split: Option<(f64, f64)>,
    order: CentralStateOrder,
) -> Result<(PositionResult, i32)> {
    serving_position_core_recorded(
        serving, jd_tt, body_id, calc_flags, want_xyz, obs_split, order, None,
    )
}

/// Original physical-query words. EQUATORIAL selects only the output frame
/// in this geocentric stored core; no other flag, epoch part or layout is erased.
fn owned_apparent_query(
    producer: std::num::NonZeroU64,
    tt: f64,
    parts: Option<(f64, f64)>,
    body: i32,
    flags: i32,
    order: CentralStateOrder,
    xyz: bool,
) -> [u64; 9] {
    let (present, hi, lo) = parts.map_or((0, 0, 0), |(hi, lo)| (1, hi.to_bits(), lo.to_bits()));
    [
        producer.get(),
        tt.to_bits(),
        present,
        hi,
        lo,
        body as u32 as u64,
        (flags & !FLG_EQUATORIAL) as u32 as u64,
        u64::from(order == CentralStateOrder::ObserverFirst),
        u64::from(xyz),
    ]
}

type OwnedApparentMotion = (
    [super::derivative::AngularWord; 3],
    [super::derivative::AngularWord; 3],
);

fn owned_apparent_words(position: [f64; 3], motion: OwnedApparentMotion) -> [f64; 15] {
    std::array::from_fn(|index| match index {
        0..=2 => position[index],
        3..=8 => {
            let word = motion.0[(index - 3) / 2];
            if index % 2 == 1 { word.hi } else { word.lo }
        }
        _ => {
            let word = motion.1[(index - 9) / 2];
            if index % 2 == 1 { word.hi } else { word.lo }
        }
    })
}

fn owned_apparent_from_words(words: [f64; 15]) -> ([f64; 3], OwnedApparentMotion) {
    use super::derivative::AngularWord as Word;
    let values = |start| {
        std::array::from_fn(|axis| Word {
            hi: words[start + axis * 2],
            lo: words[start + axis * 2 + 1],
        })
    };
    ([words[0], words[1], words[2]], (values(3), values(9)))
}

/// Output-specific projection remains outside the physical state memo. Both
/// paths retain exactly the old velocity rounding and frame/Jacobian order.
#[allow(clippy::too_many_arguments)]
fn projected_apparent_rates<P: StateProvider>(
    serving: &SealedServing<'_, P>,
    body: i32,
    split: Option<(f64, f64)>,
    tt: f64,
    flags: i32,
    position: [f64; 3],
    shaped: [f64; 3],
    frames: FrameMatrices,
    xyz: bool,
    motion: OwnedApparentMotion,
) -> Result<(f64, f64, f64)> {
    use super::derivative::AngularWord as Word;
    let velocity = motion.1.map(Word::rounded);
    let output = if xyz {
        shaped
    } else {
        frame_transform(position, tt, flags, true, frames)?
    };
    apparent_output_rates(
        serving,
        body,
        split,
        tt,
        flags,
        position,
        output,
        velocity,
        xyz,
        true,
        Some(motion),
    )
}

/// Retain the actual correction increment while evaluating the same place.
#[allow(clippy::too_many_arguments)]
fn serving_position_core_recorded<P: StateProvider>(
    serving: &SealedServing<'_, P>,
    jd_tt: f64,
    body_id: i32,
    calc_flags: i32,
    want_xyz: bool,
    obs_split: Option<(f64, f64)>,
    order: CentralStateOrder,
    capture: Option<&std::cell::Cell<AngularModelSample>>,
) -> Result<(PositionResult, i32)> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let _channel = serving_channel_gate(serving.provider, body_id)?;
    // This geocentric entry rejects requests requiring a site. Geocentrically
    // defined lunar points have no parallax and retain the geocentric path.
    match resolve_topo_request(calc_flags, body_id, None, None)? {
        TopoPlan::Geocentric | TopoPlan::TopoNoop => {}
        TopoPlan::Topocentric { .. } => {
            return Err(Error::ConfigurationMissing {
                missing_config: "observer_location",
                suggestion: "topocentric serving is not bound yet (POS-13)",
            });
        }
    }
    let center = observation_center(resolve_center_flags(calc_flags));
    let owned_key = serving
        .owned_apparent_id
        .filter(|_| {
            capture.is_none()
                && !serving.force_sampled_derivatives
                && calc_flags & FLG_SPEED != 0
                && calc_flags & FLG_TOPOCTR == 0
                && matches!(center, ObservationCenter::Geocentric)
        })
        .map(|producer| {
            owned_apparent_query(
                producer, jd_tt, obs_split, body_id, calc_flags, order, want_xyz,
            )
        });
    if let Some(words) = owned_key.and_then(crate::pure_cache::owned_apparent) {
        let (position, motion) = owned_apparent_from_words(words);
        let (dpsi, deps) = (serving.nutation_at)(jd_tt)?;
        let frames = frame_matrices_for_date(jd_tt, dpsi, deps)?;
        let shaped = frame_transform(position, jd_tt, calc_flags, want_xyz, frames)?;
        let (dlon, dlat, ddist) = projected_apparent_rates(
            serving, body_id, obs_split, jd_tt, calc_flags, position, shaped, frames, want_xyz,
            motion,
        )?;
        // Direct and public-scoped trace policy is unchanged, including failure.
        super::routing::record_stored_trace_source(body_id)?;
        return Ok((
            [shaped[0], shaped[1], shaped[2], dlon, dlat, ddist],
            calc_flags,
        ));
    }
    let mut pending_owned = None;
    let observer_body = match center {
        ObservationCenter::Heliocentric => Some(SUN),
        ObservationCenter::Geocentric => Some(EARTH),
        _ => None,
    };
    let held_observer = if order == CentralStateOrder::ObserverFirst {
        observer_body
            .map(|body| serving.provider.geometric_state(body, jd_tt))
            .transpose()?
    } else {
        None
    };
    // Scalar positions consume the target first. A physical vector place
    // instead holds the observing center before it acquires the target.
    let target = serving.provider.geometric_state(body_id, jd_tt)?;
    // Read only support bodies consumed by the selected center. An unused
    // observer cannot introduce an unrelated coverage failure.
    let (sun_pos, sun_vel) = match center {
        ObservationCenter::Heliocentric => {
            let state = match held_observer {
                Some(state) => state,
                None => serving.provider.geometric_state(SUN, jd_tt)?,
            };
            (state.pos(), state.vel())
        }
        _ => ([0.0; 3], [0.0; 3]),
    };
    let (earth_pos, earth_vel) = match center {
        ObservationCenter::Geocentric => {
            let state = match held_observer {
                Some(state) => state,
                None => serving.provider.geometric_state(EARTH, jd_tt)?,
            };
            (state.pos(), state.vel())
        }
        _ => ([0.0; 3], [0.0; 3]),
    };
    // Self-observation is the exact zero sentinel, with no defined angular
    // direction or rate. Check the central stored state first, but do not
    // acquire frame or SPEED-neighbor data for an identically zero result.
    // In particular this preserves valid central epochs at coverage edges.
    if matches!(
        (center, body_id),
        (ObservationCenter::Geocentric, EARTH) | (ObservationCenter::Heliocentric, SUN)
    ) {
        super::routing::record_stored_trace_source(body_id)?;
        return Ok(([0.0; 6], calc_flags));
    }
    let (obs_pos, obs_vel) =
        observer_state_for_center(center, sun_pos, sun_vel, earth_pos, earth_vel, None);
    let rel = compose_relative_state(
        target.pos(),
        target.vel(),
        sun_pos,
        sun_vel,
        earth_pos,
        earth_vel,
        center,
        None,
    )?;
    let lt = provider_light_time_with_initial(
        serving.provider,
        body_id,
        jd_tt,
        obs_split,
        PropagationObserver {
            position: obs_pos,
            velocity: obs_vel,
            body: observer_body,
            site_offset: [0.0; 3],
        },
        wants_light_time(calc_flags),
        (serving.reuse_immutable_central
            && calc_flags & FLG_TOPOCTR == 0
            && matches!(center, ObservationCenter::Geocentric))
        .then_some(target),
    )?;
    let reuse_central = serving.reuse_immutable_central
        && calc_flags & FLG_SPEED != 0
        && calc_flags & FLG_TOPOCTR == 0
        && matches!(center, ObservationCenter::Geocentric);
    let correction_reads = std::cell::Cell::new([None; DEFLECTORS.len() * 2]);
    // Capture only producer-owned scalar samples. Lend them only after the
    // existing complete analytic source/acceleration branch has been selected.
    let ppn_scalar_samples =
        (reuse_central && !serving.force_sampled_derivatives).then(PpnScalarSamples::default);
    let corrected = super::frames::apply_deflection_and_aberration_with_vector_route(
        lt.geo_pos(),
        rel.aberration_vel(),
        obs_pos,
        jd_tt,
        lt.light_time_days(),
        body_id,
        calc_flags,
        |id, jd| {
            capture_correction_read(
                serving.provider,
                id,
                jd,
                (capture.is_some() || reuse_central).then_some(&correction_reads),
                ppn_scalar_samples.as_ref(),
            )
        },
        None,
        obs_split.is_some(),
    )?;
    let (dpsi, deps) = (serving.nutation_at)(jd_tt)?;
    let frames = frame_matrices_for_date(jd_tt, dpsi, deps)?;
    let shaped = frame_transform(corrected.pos(), jd_tt, calc_flags, want_xyz, frames)?;
    if let Some(capture) = capture {
        capture.set(AngularModelSample {
            output: shaped,
            apparent: corrected.pos(),
            correction: corrected.correction_increment(),
            reads: correction_read_domains(serving.provider, correction_reads.get())?,
        });
    }
    let (dlon, dlat, ddist) = if (calc_flags & FLG_SPEED) != 0 {
        let do_defl = wants_deflection(calc_flags, body_id, lt.light_time_days())
            || (obs_split.is_some()
                && body_id == MOON
                && calc_flags & (FLG_NOGDEFL | FLG_HELCTR | FLG_BARYCTR | FLG_TRUEPOS) == 0
                && lt.light_time_days() > 0.0);
        let do_aber = wants_aberration(calc_flags);
        let stencil_center = StencilCenter::new(
            lt.geo_pos(),
            lt.geo_vel(),
            obs_pos,
            obs_vel,
            corrected.deflection_delta(),
            body_id,
            do_defl,
            do_aber,
            None,
        )?;
        let acceleration = if do_aber
            && !serving.force_sampled_derivatives
            && serving.provider.has_analytic_kinematics()
        {
            serving.provider.geometric_acceleration_tt(EARTH, jd_tt)?
        } else {
            None
        };
        // Borrow only the already solved observation. Eligibility and the
        // fallible acceleration read above stay outside both rate routes.
        let rate_observation = rates::RateObservation {
            jd_tt,
            obs_split,
            observer_body,
            light_time: &lt,
            stencil: &stencil_center,
        };
        if let Some(acceleration) = acceleration {
            // This borrowed read view is created only after the original
            // completeness/acceleration gate. It cannot select another route.
            let source = super::state_reads::ProviderReads::new(serving.provider);
            let motion = rates::analytic_apparent_motion(
                &source,
                &rate_observation,
                acceleration,
                ppn_scalar_samples.as_ref(),
            )?;
            if owned_key.is_some() {
                pending_owned = Some(owned_apparent_words(corrected.pos(), motion));
            }
            projected_apparent_rates(
                serving,
                body_id,
                obs_split,
                jd_tt,
                calc_flags,
                corrected.pos(),
                shaped,
                frames,
                want_xyz,
                motion,
            )?
        } else {
            let central = reuse_central.then(|| CentralAngularPlace {
                query: CentralAngularQuery::new(jd_tt, obs_split, body_id, calc_flags, order),
                apparent: corrected.pos(),
                correction: corrected.correction_increment(),
                frames: &frames,
                reads: correction_reads.get(),
            });
            rates::sampled_apparent_rates(
                serving,
                &rate_observation,
                calc_flags,
                want_xyz,
                order,
                central.as_ref(),
            )?
        }
    } else {
        (0.0, 0.0, 0.0)
    };
    let out = (
        [shaped[0], shaped[1], shaped[2], dlon, dlat, ddist],
        calc_flags,
    );
    // Direct serving retains its legacy LEB hook. An explicit public branch
    // defers the inner record until its output/speed/frame tail completes;
    // support reads must not leak a remapped body or replace a fallback label.
    super::routing::record_stored_trace_source(body_id)?;
    if let (Some(key), Some(words)) = (owned_key, pending_owned) {
        if words.iter().all(|value| value.is_finite()) {
            crate::pure_cache::retain_owned_apparent(key, words);
        }
    }
    Ok(out)
}

/// Compute a stored apparent place for a geodetic observer.
///
/// `jd_tt` and `jd_ut1` are paired Julian Days. `geopos` contains east longitude,
/// latitude in degrees and height in metres; `c2t_at` supplies their
/// celestial-to-terrestrial orientation. `site_ellipsoid` selects published
/// WGS84 or IERS 2010 geodetic parameters. `earth_term` optionally enables
/// Earth's surface-observer deflection for the vector-phenomena model.
///
/// The full site displacement enters observer and apparent-place geometry.
/// Without a real parallax request, delegate to [`serving_compute_tt`].
/// The default sampler holds central TT-minus-UT1 fixed; context-aware entry
/// points can supply the configured variable-clock callback. Results use the
/// shared degrees/AU layout and per-TT-day rates. Provider, clock, orientation,
/// site and correction errors propagate.
#[allow(clippy::too_many_arguments)]
pub fn serving_compute_tt_topocentric<P: StateProvider>(
    serving: &SealedServing<'_, P>,
    jd_tt: f64,
    jd_ut1: f64,
    body_id: i32,
    calc_flags: i32,
    geopos: (f64, f64, f64),
    c2t_at: &dyn Fn(f64, f64) -> Result<[[f64; 3]; 3]>,
    earth_term: bool,
    site_ellipsoid: SiteEllipsoid,
) -> Result<(PositionResult, i32)> {
    match serving_topocentric_core(
        serving,
        jd_tt,
        jd_ut1,
        body_id,
        calc_flags,
        geopos,
        c2t_at,
        false,
        earth_term,
        site_ellipsoid,
        None,
        None,
    ) {
        Err(Error::LebOutOfRange { .. })
            if vector_fallback_central_covered(serving.provider, body_id, jd_tt) =>
        {
            serving_vector_state_fallback(serving, jd_tt, body_id, calc_flags, Some(geopos), false)
        }
        outcome => outcome,
    }
}

/// Direct topocentric Cartesian position and velocity, in AU and AU/day.
///
/// The site, TT/UT1 partners, observer correction, c2t provider, frozen
/// sample-timescale difference, site ellipsoid, domains and errors follow
/// [`serving_compute_tt_topocentric`]. Topocentric rates retain their sampled
/// differentiation, with central and support vectors kept Cartesian.
/// Non-topocentric requests delegate to [`serving_compute_xyz_tt`], including
/// that path's analytic eligibility rather than forcing a sampled convention.
#[allow(clippy::too_many_arguments)]
pub fn serving_compute_xyz_tt_topocentric<P: StateProvider>(
    serving: &SealedServing<'_, P>,
    jd_tt: f64,
    jd_ut1: f64,
    body_id: i32,
    calc_flags: i32,
    geopos: (f64, f64, f64),
    c2t_at: &dyn Fn(f64, f64) -> Result<[[f64; 3]; 3]>,
    earth_term: bool,
    site_ellipsoid: SiteEllipsoid,
) -> Result<(PositionResult, i32)> {
    match serving_topocentric_core(
        serving,
        jd_tt,
        jd_ut1,
        body_id,
        calc_flags,
        geopos,
        c2t_at,
        true,
        earth_term,
        site_ellipsoid,
        None,
        None,
    ) {
        Err(Error::LebOutOfRange { .. })
            if vector_fallback_central_covered(serving.provider, body_id, jd_tt) =>
        {
            serving_vector_state_fallback(serving, jd_tt, body_id, calc_flags, Some(geopos), true)
        }
        outcome => outcome,
    }
}

/// Compute a topocentric stored place with an optional split TT observation.
///
/// Observer, ellipsoid, flags, output and error policy are those of
/// [`serving_compute_tt_topocentric`]. The optional TT whole/fraction parts
/// retain short light-time offsets without full-Julian-Day cancellation.
#[allow(clippy::too_many_arguments)]
pub fn serving_compute_tt_topocentric_split<P: StateProvider>(
    serving: &SealedServing<'_, P>,
    jd_tt: f64,
    jd_ut1: f64,
    body_id: i32,
    calc_flags: i32,
    geopos: (f64, f64, f64),
    c2t_at: &dyn Fn(f64, f64) -> Result<[[f64; 3]; 3]>,
    earth_term: bool,
    site_ellipsoid: SiteEllipsoid,
    obs_split: Option<(f64, f64)>,
) -> Result<(PositionResult, i32)> {
    match serving_topocentric_core(
        serving,
        jd_tt,
        jd_ut1,
        body_id,
        calc_flags,
        geopos,
        c2t_at,
        false,
        earth_term,
        site_ellipsoid,
        obs_split,
        None,
    ) {
        Err(Error::LebOutOfRange { .. })
            if vector_fallback_central_covered(serving.provider, body_id, jd_tt) =>
        {
            serving_vector_state_fallback(serving, jd_tt, body_id, calc_flags, Some(geopos), false)
        }
        outcome => outcome,
    }
}

/// Topocentric spherical stored-state computation without a range retry.
///
/// TT/UT1 clocks, east-positive geodetic degrees/metres, rotation callback,
/// ellipsoid, optional Earth deflection and observation split retain the
/// contracts of [`serving_compute_tt_topocentric_split`]. Returned units are
/// degrees/degrees/AU and daily rates. The central and all derived support
/// failures propagate so a multi-leg caller can retry its entire calculation
/// using its own fallback clock. No individual leg selects a vector fallback.
#[allow(clippy::too_many_arguments)]
pub fn serving_compute_tt_topocentric_direct<P: StateProvider>(
    serving: &SealedServing<'_, P>,
    jd_tt: f64,
    jd_ut1: f64,
    body_id: i32,
    calc_flags: i32,
    geopos: (f64, f64, f64),
    c2t_at: &dyn Fn(f64, f64) -> Result<[[f64; 3]; 3]>,
    earth_term: bool,
    site_ellipsoid: SiteEllipsoid,
    obs_split: Option<(f64, f64)>,
) -> Result<(PositionResult, i32)> {
    serving_topocentric_core(
        serving,
        jd_tt,
        jd_ut1,
        body_id,
        calc_flags,
        geopos,
        c2t_at,
        false,
        earth_term,
        site_ellipsoid,
        obs_split,
        None,
    )
}

/// Stored topocentric place with the caller's actual TT-to-UT1 model.
/// The clock is evaluated at every angular sample, while the central UT1 word
/// remains supplied explicitly. Geometry is the same native place calculation
/// with SPEED disabled; it never substitutes a linearized orbit or fixed delay.
#[allow(clippy::too_many_arguments)]
pub(super) fn serving_compute_topocentric_with_clock<P: StateProvider>(
    serving: &SealedServing<'_, P>,
    jd_tt: f64,
    jd_ut1: f64,
    body_id: i32,
    calc_flags: i32,
    geopos: (f64, f64, f64),
    c2t_at: &dyn Fn(f64, f64) -> Result<[[f64; 3]; 3]>,
    want_xyz: bool,
    ut1_at: &dyn Fn(f64) -> Result<f64>,
    clock_support: (f64, f64),
) -> Result<(PositionResult, i32)> {
    match serving_topocentric_core(
        serving,
        jd_tt,
        jd_ut1,
        body_id,
        calc_flags,
        geopos,
        c2t_at,
        want_xyz,
        false,
        SiteEllipsoid::Wgs84,
        None,
        Some((ut1_at, clock_support)),
    ) {
        Err(Error::LebOutOfRange { .. })
            if vector_fallback_central_covered(serving.provider, body_id, jd_tt) =>
        {
            serving_vector_state_fallback(
                serving,
                jd_tt,
                body_id,
                calc_flags,
                Some(geopos),
                want_xyz,
            )
        }
        outcome => outcome,
    }
}

/// Shared stored observer and apparent-place reduction.
///
/// `want_xyz` selects Cartesian AU/AU-day rather than spherical degrees/AU.
/// `earth_term` enables surface-Earth deflection; `site_ellipsoid` selects
/// geodetic parameters. Optional `obs_split` retains the reception TT parts
/// for light-time iteration. Central and sampled observers share the same
/// configuration and clock policy.
#[allow(clippy::too_many_arguments)]
fn serving_topocentric_core<P: StateProvider>(
    serving: &SealedServing<'_, P>,
    jd_tt: f64,
    jd_ut1: f64,
    body_id: i32,
    calc_flags: i32,
    geopos: (f64, f64, f64),
    c2t_at: &dyn Fn(f64, f64) -> Result<[[f64; 3]; 3]>,
    want_xyz: bool,
    earth_term: bool,
    site_ellipsoid: SiteEllipsoid,
    obs_split: Option<(f64, f64)>,
    ut1_at: Option<SiteClock<'_>>,
) -> Result<(PositionResult, i32)> {
    serving_topocentric_core_ordered(
        serving,
        jd_tt,
        jd_ut1,
        body_id,
        calc_flags,
        geopos,
        c2t_at,
        want_xyz,
        earth_term,
        site_ellipsoid,
        obs_split,
        CentralStateOrder::TargetFirst,
        ut1_at,
    )
}

/// Topocentric numerical reduction with an explicit central acquisition order.
/// Vector observation holds the center before the surface frame; scalar
/// positions keep their original surface-frame, target, then observer stages.
#[allow(clippy::too_many_arguments)]
fn serving_topocentric_core_ordered<P: StateProvider>(
    serving: &SealedServing<'_, P>,
    jd_tt: f64,
    jd_ut1: f64,
    body_id: i32,
    calc_flags: i32,
    geopos: (f64, f64, f64),
    c2t_at: &dyn Fn(f64, f64) -> Result<[[f64; 3]; 3]>,
    want_xyz: bool,
    earth_term: bool,
    site_ellipsoid: SiteEllipsoid,
    obs_split: Option<(f64, f64)>,
    order: CentralStateOrder,
    ut1_at: Option<SiteClock<'_>>,
) -> Result<(PositionResult, i32)> {
    serving_topocentric_core_recorded(
        serving,
        jd_tt,
        jd_ut1,
        body_id,
        calc_flags,
        geopos,
        c2t_at,
        want_xyz,
        earth_term,
        site_ellipsoid,
        obs_split,
        order,
        ut1_at,
        None,
    )
}

/// Surface-observer place with a retained algebraic correction increment.
#[allow(clippy::too_many_arguments)]
fn serving_topocentric_core_recorded<P: StateProvider>(
    serving: &SealedServing<'_, P>,
    jd_tt: f64,
    jd_ut1: f64,
    body_id: i32,
    calc_flags: i32,
    geopos: (f64, f64, f64),
    c2t_at: &dyn Fn(f64, f64) -> Result<[[f64; 3]; 3]>,
    want_xyz: bool,
    earth_term: bool,
    site_ellipsoid: SiteEllipsoid,
    obs_split: Option<(f64, f64)>,
    order: CentralStateOrder,
    ut1_at: Option<SiteClock<'_>>,
    capture: Option<&std::cell::Cell<AngularModelSample>>,
) -> Result<(PositionResult, i32)> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    if !jd_ut1.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_ut1 });
    }
    let _channel = serving_channel_gate(serving.provider, body_id)?;
    let (lon_deg, lat_deg, height_m) =
        match resolve_topo_request(calc_flags, body_id, Some(geopos), None)? {
            TopoPlan::Geocentric | TopoPlan::TopoNoop => {
                let split = match order {
                    CentralStateOrder::TargetFirst => None,
                    CentralStateOrder::ObserverFirst => obs_split,
                };
                return serving_position_core_recorded(
                    serving, jd_tt, body_id, calc_flags, want_xyz, split, order, capture,
                );
            }
            TopoPlan::Topocentric {
                lon_deg,
                lat_deg,
                height_m,
            } => (lon_deg, lat_deg, height_m),
        };
    let observer_body = if calc_flags & FLG_HELCTR != 0 {
        Some(SUN)
    } else if calc_flags & FLG_BARYCTR != 0 {
        None
    } else {
        Some(EARTH)
    };
    let held_observer = if order == CentralStateOrder::ObserverFirst {
        observer_body
            .map(|body| serving.provider.geometric_state(body, jd_tt))
            .transpose()?
    } else {
        None
    };
    let c2t = c2t_at(jd_tt, jd_ut1)?;
    let central = topocentric_offset(lon_deg, lat_deg, height_m, c2t, site_ellipsoid)?;
    let (topo_pos, topo_vel) = central.pair();
    let target = serving.provider.geometric_state(body_id, jd_tt)?;
    // Only a geocentric topocentric center consumes the site offset.
    let (obs_pos, obs_vel, earth_vel, compose_center, compose_topo) =
        if (calc_flags & FLG_HELCTR) != 0 {
            let state = match held_observer {
                Some(state) => state,
                None => serving.provider.geometric_state(SUN, jd_tt)?,
            };
            (
                state.pos(),
                state.vel(),
                [0.0; 3],
                ObservationCenter::Heliocentric,
                None,
            )
        } else if (calc_flags & FLG_BARYCTR) != 0 {
            (
                [0.0; 3],
                [0.0; 3],
                [0.0; 3],
                ObservationCenter::Barycentric,
                None,
            )
        } else {
            let state = match held_observer {
                Some(state) => state,
                None => serving.provider.geometric_state(EARTH, jd_tt)?,
            };
            let (ep, ev) = (state.pos(), state.vel());
            let obs = [
                ep[0] + topo_pos[0],
                ep[1] + topo_pos[1],
                ep[2] + topo_pos[2],
            ];
            let obs_vel = [
                ev[0] + topo_vel[0],
                ev[1] + topo_vel[1],
                ev[2] + topo_vel[2],
            ];
            (
                obs,
                obs_vel,
                obs_vel,
                ObservationCenter::Topocentric,
                Some((topo_pos, topo_vel)),
            )
        };
    // Heliocentric work reuses the Sun observer. Barycentric work needs neither
    // Sun nor Earth for center subtraction; other centers acquire both supports.
    let (sun_pos, sun_vel) = match compose_center {
        ObservationCenter::Heliocentric => (obs_pos, obs_vel),
        ObservationCenter::Barycentric => ([0.0; 3], [0.0; 3]),
        _ => {
            let sun = serving.provider.geometric_state(SUN, jd_tt)?;
            (sun.pos(), sun.vel())
        }
    };
    let (earth_pos, earth_vel_state) = match compose_center {
        ObservationCenter::Barycentric | ObservationCenter::Heliocentric => ([0.0; 3], [0.0; 3]),
        _ => {
            let earth = serving.provider.geometric_state(EARTH, jd_tt)?;
            (earth.pos(), earth.vel())
        }
    };
    let rel = compose_relative_state(
        target.pos(),
        target.vel(),
        sun_pos,
        sun_vel,
        earth_pos,
        earth_vel_state,
        compose_center,
        compose_topo,
    )?;
    // Observer selection and relative-state composition must agree.
    debug_assert_eq!(obs_pos, {
        let (op, _) = observer_state_for_center(
            compose_center,
            sun_pos,
            sun_vel,
            earth_pos,
            earth_vel_state,
            compose_topo,
        );
        op
    });
    let lt = provider_light_time(
        serving.provider,
        body_id,
        jd_tt,
        obs_split,
        PropagationObserver {
            position: obs_pos,
            velocity: obs_vel,
            body: observer_body,
            site_offset: compose_topo
                .map(|(position, _)| position)
                .unwrap_or([0.0; 3]),
        },
        wants_light_time(calc_flags),
    )?;
    let correction_reads = std::cell::Cell::new([None; DEFLECTORS.len() * 2]);
    let corrected = apply_deflection_and_aberration(
        lt.geo_pos(),
        rel.aberration_vel(),
        obs_pos,
        jd_tt,
        lt.light_time_days(),
        body_id,
        calc_flags,
        |id, jd| {
            capture_correction_read(
                serving.provider,
                id,
                jd,
                capture.map(|_| &correction_reads),
                None,
            )
        },
        // Enable Earth's surface-deflection contribution only for the requested
        // apparent vector model with a real site displacement.
        if earth_term && (calc_flags & FLG_TRUEPOS) == 0 {
            Some(topo_pos)
        } else {
            None
        },
    )?;
    // `earth_vel` drives aberration: the observer velocity on the
    // topocentric branch (diurnal aberration), zeros elsewhere.
    debug_assert_eq!(rel.aberration_vel(), earth_vel);
    let (dpsi, deps) = (serving.nutation_at)(jd_tt)?;
    let frames = frame_matrices_for_date(jd_tt, dpsi, deps)?;
    let shaped = frame_transform(corrected.pos(), jd_tt, calc_flags, want_xyz, frames)?;
    if let Some(capture) = capture {
        capture.set(AngularModelSample {
            output: shaped,
            apparent: corrected.pos(),
            correction: corrected.correction_increment(),
            reads: correction_read_domains(serving.provider, correction_reads.get())?,
        });
    }
    let (dlon, dlat, ddist) = if (calc_flags & FLG_SPEED) != 0 {
        let do_defl = wants_deflection(calc_flags, body_id, lt.light_time_days())
            // Sampled lunar correction gates must match the central surface-vector model.
            || (earth_term
                && body_id == MOON
                && (calc_flags & (FLG_NOGDEFL | FLG_HELCTR | FLG_BARYCTR | FLG_TRUEPOS)) == 0
                && lt.light_time_days() > 0.0);
        let do_aber = wants_aberration(calc_flags);
        // Construct the site sampler with the central TT-minus-UT1 difference.
        let sampler = TopoSampler::new(
            lon_deg,
            lat_deg,
            height_m,
            frozen_delta_t_days(jd_tt, jd_ut1)?,
            site_ellipsoid,
        )?;
        let topo_at = |jd_sample: f64| -> Result<([f64; 3], [f64; 3])> {
            let jd_ut1_sample = sampler.ut1_at(jd_sample)?;
            let sample_c2t = c2t_at(jd_sample, jd_ut1_sample)?;
            Ok(sampler.offset_at(jd_sample, sample_c2t)?.pair())
        };
        let stencil_center = StencilCenter::new(
            lt.geo_pos(),
            lt.geo_vel(),
            obs_pos,
            obs_vel,
            corrected.deflection_delta(),
            body_id,
            do_defl,
            do_aber,
            Some(StencilTopo::new(topo_pos, topo_vel, None)?),
        )?;
        let angular_at = |offset: f64| {
            let epoch = jd_tt + offset;
            let sample_ut1 = match (offset == 0.0, ut1_at) {
                (true, _) => jd_ut1,
                (false, Some((at, _))) => at(epoch)?,
                (false, None) => sampler.ut1_at(epoch)?,
            };
            let split = obs_split.map(|_| {
                let whole = epoch.floor();
                (whole, epoch - whole)
            });
            let capture = std::cell::Cell::new(AngularModelSample::default());
            serving_topocentric_core_recorded(
                serving,
                epoch,
                sample_ut1,
                body_id,
                calc_flags & !(FLG_SPEED | FLG_SPEED3),
                geopos,
                c2t_at,
                true,
                earth_term,
                site_ellipsoid,
                split,
                order,
                ut1_at,
                Some(&capture),
            )
            .map(|_| capture.get())
        };
        let mut support = angular_geometry_support(
            serving.provider,
            body_id,
            observer_body,
            jd_tt,
            lt.retarded_jd_tt(),
        )?;
        if let Some((_, clock_support)) = ut1_at {
            support.0 = support.0.max(clock_support.0);
            support.1 = support.1.min(clock_support.1);
        }
        // The derivative of a rotating site follows the actual supplied
        // observation clock and rotation model. Its smooth support is
        // independent of the target's much shorter retarded record.
        let geometric_velocity = if compose_topo.is_some() {
            let site_support =
                ut1_at.map_or((f64::NEG_INFINITY, f64::INFINITY), |(_, support)| support);
            let site_velocity = super::derivative::cartesian_rates(
                &|offset| {
                    let epoch = jd_tt + offset;
                    let ut = match (offset == 0.0, ut1_at) {
                        (true, _) => jd_ut1,
                        (false, Some((at, _))) => at(epoch)?,
                        (false, None) => sampler.ut1_at(epoch)?,
                    };
                    Ok(topocentric_offset(
                        lon_deg,
                        lat_deg,
                        height_m,
                        c2t_at(epoch, ut)?,
                        site_ellipsoid,
                    )?
                    .pair()
                    .0)
                },
                1.0 / 64.0,
                site_support,
            )?;
            provider_light_time(
                serving.provider,
                body_id,
                jd_tt,
                obs_split,
                PropagationObserver {
                    position: obs_pos,
                    velocity: std::array::from_fn(|i| earth_vel_state[i] + site_velocity[i]),
                    body: observer_body,
                    site_offset: compose_topo.map_or([0.0; 3], |(position, _)| position),
                },
                wants_light_time(calc_flags),
            )?
            .geo_vel()
        } else {
            lt.geo_vel()
        };
        local_model_velocity(
            serving,
            &stencil_center,
            geometric_velocity,
            jd_tt,
            calc_flags,
            want_xyz,
            earth_term,
            Some(&topo_at),
            &angular_at,
            support,
        )?
    } else {
        (0.0, 0.0, 0.0)
    };
    let out = (
        [shaped[0], shaped[1], shaped[2], dlon, dlat, ddist],
        calc_flags,
    );
    // Trace hook, topocentric-arm counterpart of the `serving_compute_tt`
    // hook above: this sealed LEB leg served `body_id`.
    super::routing::record_stored_trace_source(body_id)?;
    Ok(out)
}

/// Nutation and obliquity at a UT1 Julian date using the bundled timescale.
///
/// Converts UT1 to TT with [`crate::time::default_timescale_tt_from_ut1`],
/// then evaluates IAU 2006/2000A nutation and Vondrák long-term mean
/// obliquity. Returns the six-component degree tuple
/// `(true_obliquity, mean_obliquity, nutation_longitude, nutation_obliquity,
/// 0, 0)` and echoes `nut_flags`. This is a time/frame quantity, not a body
/// coordinate: it does not consume a LEB nutation channel or the user's
/// Delta-T override/model. Non-finite and unrepresentable epochs propagate
/// typed errors. Series evaluation is not a physical-accuracy claim outside
/// the underlying model's astronomical validity.
pub fn ecl_nut_result_ut(jd_ut: f64, nut_flags: i32) -> Result<(PositionResult, i32)> {
    let jd_tt = crate::time::default_timescale_tt_from_ut1(jd_ut)?;
    ecl_nut_result(jd_tt, nut_flags, &crate::erfa::nut06a)
}

/// Nutation and obliquity at a TT Julian date over an explicit nutation source.
///
/// Returns `(true_obliquity, mean_obliquity, nutation_longitude,
/// nutation_obliquity, 0, 0)` in degrees and echoes `nut_flags`. Mean
/// obliquity is the directly evaluated Vondrák et al. (2011) epsilon-A
/// series, not a short-range polynomial extrapolation; true obliquity is
/// `mean_obliquity + nutation_obliquity`. The supplied `nutation_at` returns
/// longitude/obliquity nutation in radians at the TT argument.
///
/// Native calculation entries use [`crate::erfa::nut06a`] (IAU 2006/2000A).
/// Explicit callers may supply another pair; its numerical/physical model
/// limits remain their responsibility. No ephemeris-body output shaping or
/// velocity computation applies. Non-finite epochs and pair components
/// return typed errors. See [`vondrak_mean_obliquity_deg`] for the long-term
/// mean-obliquity construction and domain.
pub fn ecl_nut_result(
    jd_tt: f64,
    nut_flags: i32,
    nutation_at: &dyn Fn(f64) -> Result<(f64, f64)>,
) -> Result<(PositionResult, i32)> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let mean_obliquity = vondrak_mean_obliquity_deg(jd_tt)?;
    let (dpsi_rad, deps_rad) = nutation_at(jd_tt)?;
    for v in [dpsi_rad, deps_rad] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    let dpsi_deg = dpsi_rad.to_degrees();
    let deps_deg = deps_rad.to_degrees();
    Ok((
        [
            mean_obliquity + deps_deg,
            mean_obliquity,
            dpsi_deg,
            deps_deg,
            0.0,
            0.0,
        ],
        nut_flags,
    ))
}

#[cfg(test)]
mod source_time_domain_tests {
    use super::*;
    use std::cell::Cell;

    #[test]
    fn source_time_tail_dispatches_after_internal_ieee_clock() {
        let calls = Cell::new(0);
        let compute = |tt: f64, body: i32, flags: i32| {
            assert!(tt.is_nan());
            assert_eq!(body, WHITE_MOON);
            assert_eq!(flags, FLG_SWIEPH);
            calls.set(calls.get() + 1);
            Ok(([f64::NAN, 0.0, f64::NAN, 0.0, 0.0, 0.0], flags))
        };
        let engine = CalcEngine {
            ut1_to_tt_jd: &|_| panic!("already acquired"),
            compute_tt: &compute,
            nutation_obliquity_ut: &|_, _| panic!("not ECL_NUT"),
            nutation_obliquity_tt: &|_, _| panic!("not ECL_NUT"),
            jpl_fallback_active: false,
        };
        let (state, flags) =
            calc_tt_body_source_time(&engine, f64::NAN, WHITE_MOON, FLG_SWIEPH, FLG_SWIEPH)
                .unwrap();
        assert_eq!(calls.get(), 1);
        assert_eq!(flags, FLG_SWIEPH);
        assert_eq!(
            state.map(f64::to_bits),
            [0x7ff8_0000_0000_0000, 0, 0x7ff8_0000_0000_0000, 0, 0, 0]
        );
        assert!(matches!(
            calc_tt_body(&engine, f64::NAN, WHITE_MOON, FLG_SWIEPH, FLG_SWIEPH),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert_eq!(calls.get(), 1);
    }
}

#[cfg(test)]
mod prepared_flag_tests {
    use super::*;
    use std::cell::Cell;

    #[test]
    fn prepared_speed3_keeps_single_normalization_and_legacy_numerical_tail() {
        for is_ut in [false, true] {
            let calls = Cell::new(0);
            let compute = |tt: f64, _body, flags| {
                calls.set(calls.get() + 1);
                Ok(([tt - FRAME_J2000, 0.0, 1.0, 0.0, 0.0, 0.0], flags))
            };
            let engine = CalcEngine {
                ut1_to_tt_jd: &|ut| Ok(ut),
                compute_tt: &compute,
                nutation_obliquity_ut: &|_, _| panic!("unused nutation"),
                nutation_obliquity_tt: &|_, _| panic!("unused nutation"),
                jpl_fallback_active: false,
            };
            let raw = FLG_SPEED3;
            let normalized = normalize_calc_flags(raw, false);
            let old = if is_ut {
                calc_ut(&engine, FRAME_J2000, MEAN_NODE, raw)
            } else {
                calc(&engine, FRAME_J2000, MEAN_NODE, raw)
            }
            .unwrap();
            let old_calls = calls.replace(0);
            let flags = PreparedCalcFlags {
                raw,
                normalized,
                defer_echo: false,
            };
            let prepared = if is_ut {
                calc_ut_prepared(&engine, FRAME_J2000, MEAN_NODE, flags)
            } else {
                calc_prepared(&engine, FRAME_J2000, MEAN_NODE, flags)
            }
            .unwrap();
            assert_eq!(old.0.map(f64::to_bits), prepared.0.map(f64::to_bits));
            assert_eq!(old.1, prepared.1);
            assert_eq!(calls.get(), old_calls);
            assert!(
                old_calls > 1,
                "SPEED3 samples were dropped by a second normalization"
            );
            assert_ne!(prepared.1 & FLG_SPEED3, 0);
            assert_eq!(prepared.1 & FLG_SPEED, 0);
        }
    }

    #[test]
    fn prepared_deferred_echo_keeps_backend_bits_until_language_success_tail() {
        let compute = |_: f64, _: i32, word: i32| Ok(([12.0, 0.0, 1.0, 2.0, 0.0, 0.0], word));
        let engine = CalcEngine {
            ut1_to_tt_jd: &|ut| Ok(ut),
            compute_tt: &compute,
            nutation_obliquity_ut: &|_, _| panic!("unused nutation"),
            nutation_obliquity_tt: &|_, _| panic!("unused nutation"),
            jpl_fallback_active: false,
        };
        // A raw SPEED3|MOSEPH request computes SPEED|SWIEPH, and the
        // language owner performs the actual raw echo only after success.
        let raw = FLG_SPEED3 | FLG_MOSEPH;
        let normalized = normalize_calc_flags(raw, false);
        let prepared = PreparedCalcFlags {
            raw: normalized,
            normalized,
            defer_echo: true,
        };
        let result = calc_prepared(&engine, FRAME_J2000, MARS, prepared).unwrap();
        assert_eq!(result.1, normalized);
        let old = calc(&engine, FRAME_J2000, MARS, raw).unwrap();
        assert_eq!(result.0, old.0);
        assert_eq!(echo_request_bits(result.1, raw), old.1);
        let zero = PreparedCalcFlags {
            raw: FLG_SWIEPH,
            normalized: FLG_SWIEPH,
            defer_echo: true,
        };
        assert_eq!(
            calc_prepared(&engine, FRAME_J2000, MARS, zero).unwrap().1,
            FLG_SWIEPH
        );
        assert_eq!(calc(&engine, FRAME_J2000, MARS, 0).unwrap().1, 0);
    }

    #[test]
    fn prepared_cartesian_uses_existing_adapter_and_restores_format_echo() {
        let flags = FLG_XYZ | FLG_RADIANS | FLG_SPEED;
        let compute = |_, _, _| panic!("unused spherical provider");
        let engine = CalcEngine {
            ut1_to_tt_jd: &|ut| Ok(ut),
            compute_tt: &compute,
            nutation_obliquity_ut: &|_, _| panic!("unused nutation"),
            nutation_obliquity_tt: &|_, _| panic!("unused nutation"),
            jpl_fallback_active: false,
        };
        let cartesian = |_, _, word| Ok(Some(([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], word)));
        let old = calc_with_cartesian(&engine, &cartesian, FRAME_J2000, MARS, flags).unwrap();
        let prepared = calc_with_cartesian_prepared(
            &engine,
            &cartesian,
            FRAME_J2000,
            MARS,
            PreparedCalcFlags {
                raw: flags,
                normalized: normalize_calc_flags(flags, false),
                defer_echo: false,
            },
            false,
        )
        .unwrap();
        assert_eq!(old, prepared);
        assert_eq!(prepared.0, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        assert_eq!(prepared.1 & (FLG_XYZ | FLG_RADIANS), FLG_XYZ | FLG_RADIANS);
    }
}

#[cfg(test)]
mod output_motion_capture_tests {
    use super::super::derivative::{AngularWord as Word, OutputWordMotion};
    use super::*;
    use std::cell::{Cell, RefCell};

    fn motion(value: f64) -> OutputWordMotion {
        OutputWordMotion {
            position: [Word::from(value); 3],
            velocity: [Word::from(value + 1.0); 3],
        }
    }
    #[test]
    fn output_motion_identity_keeps_exact_tt_parts_body_flags_and_layout() {
        let capture = OutputMotionCapture::default();
        let query = OutputMotionQuery::new(0.0, None, MARS, FLG_SPEED, false);
        let wrong = [
            OutputMotionQuery::new(-0.0, None, MARS, FLG_SPEED, false),
            OutputMotionQuery::new(0.0, Some((0.0, 0.0)), MARS, FLG_SPEED, false),
            OutputMotionQuery::new(0.0, Some((0.0, -0.0)), MARS, FLG_SPEED, false),
            OutputMotionQuery::new(0.0, None, JUPITER, FLG_SPEED, false),
            OutputMotionQuery::new(0.0, None, MARS, FLG_SPEED | FLG_NONUT, false),
            OutputMotionQuery::new(0.0, None, MARS, FLG_SPEED, true),
        ];
        let (_, captured) = capture
            .with_motion(query, || {
                for key in wrong {
                    capture.publish(key, motion(7.0));
                }
                Ok(())
            })
            .unwrap();
        assert!(captured.is_none());
        assert!(capture.request.get().is_none());
    }
    #[test]
    fn output_motion_nesting_errors_and_unwind_restore_the_previous_scope() {
        let capture = OutputMotionCapture::default();
        let outer = OutputMotionQuery::new(2451545.0, None, MARS, FLG_SPEED, false);
        let inner = OutputMotionQuery::new(2451546.0, None, JUPITER, FLG_SPEED, false);
        let (_, got) = capture
            .with_motion(outer, || {
                capture.publish(outer, motion(1.0));
                let (_, got) = capture.with_motion(inner, || {
                    capture.publish(inner, motion(3.0));
                    Ok(())
                })?;
                assert_eq!(got.unwrap().position[0].rounded(), 3.0);
                let error = capture.with_motion(inner, || {
                    capture.publish(inner, motion(4.0));
                    Err::<(), _>(Error::LebBodyNotFound { body_id: JUPITER })
                });
                assert!(error.is_err());
                let panic = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    let _ = capture.with_motion::<()>(inner, || {
                        capture.publish(inner, motion(5.0));
                        panic!("intentional unwind")
                    });
                }));
                assert!(panic.is_err());
                assert_eq!(
                    capture.request.get().unwrap().motion.unwrap().position[0].rounded(),
                    1.0
                );
                Ok(())
            })
            .unwrap();
        assert_eq!(got.unwrap().position[0].rounded(), 1.0);
        assert!(capture.request.get().is_none());
    }

    struct NativeModel {
        reads: RefCell<Vec<(i32, u64)>>,
    }
    impl StateProvider for NativeModel {
        fn body_coverage(&self, body: i32) -> Option<(f64, f64)> {
            self.has_body(body).then_some((2451500.0, 2451600.0))
        }
        fn eval_raw_two_part(
            &self,
            body: i32,
            epoch: TwoPartEpoch,
        ) -> Result<(Vec<f64>, Vec<f64>)> {
            self.eval_raw(body, epoch.value())
        }
        fn serving_entry_snapshot(&self, body: i32) -> Result<EntrySnapshot> {
            if self.has_body(body) {
                Ok(EntrySnapshot {
                    coord_type: COORD_ICRS_BARY,
                })
            } else {
                Err(Error::LebBodyNotFound { body_id: body })
            }
        }
        fn provider_name(&self) -> &'static str {
            "native-word-motion-test"
        }
        fn has_body(&self, body: i32) -> bool {
            matches!(body, EARTH | MARS | CERES)
        }
        fn serving_coord_type(&self, body: i32) -> Option<u32> {
            self.has_body(body).then_some(COORD_ICRS_BARY)
        }
        fn has_analytic_kinematics(&self) -> bool {
            true
        }
        fn geometric_acceleration_tt(&self, body: i32, _: f64) -> Result<Option<[f64; 3]>> {
            Ok(self.has_body(body).then_some([0.0; 3]))
        }
        fn eval_raw(&self, body: i32, tt: f64) -> Result<(Vec<f64>, Vec<f64>)> {
            self.reads.borrow_mut().push((body, tt.to_bits()));
            let (p, v) = match body {
                EARTH => ([0.1, -0.5, 0.2], [0.001, 0.0001, 0.0]),
                MARS | CERES => ([1.0, 2.0, 0.5], [0.001, 0.003, -0.002]),
                _ => return Err(Error::LebBodyNotFound { body_id: body }),
            };
            Ok((
                std::array::from_fn::<_, 3, _>(|i| p[i] + v[i] * (tt - 2451545.0)).to_vec(),
                v.to_vec(),
            ))
        }
    }
    #[test]
    fn output_motion_public_entry_canonicalizes_aliases_and_default_ephemeris_before_keying() {
        let model = NativeModel {
            reads: RefCell::new(Vec::new()),
        };
        let nutation = |_: f64| Ok((0.0, 0.0));
        let serving = SealedServing::new(&model, &nutation);
        let hits = Cell::new(0);
        let compute = |tt: f64, body: i32, flags: i32| {
            assert_eq!(body, CERES);
            assert_eq!(flags & (FLG_SWIEPH | FLG_JPLEPH | FLG_MOSEPH), FLG_SWIEPH);
            let (result, motion) = serving.with_output_motion(tt, body, flags, || {
                serving_compute_tt(&serving, tt, body, flags)
            })?;
            assert!(
                motion.is_some(),
                "the actual normalized query must fill the carrier"
            );
            hits.set(hits.get() + 1);
            Ok(result)
        };
        let engine = CalcEngine {
            ut1_to_tt_jd: &|_| panic!("TT entry"),
            compute_tt: &compute,
            nutation_obliquity_ut: &|_, _| panic!("not ECL_NUT"),
            nutation_obliquity_tt: &|_, _| panic!("not ECL_NUT"),
            jpl_fallback_active: false,
        };
        for ephemeris in [0, FLG_SWIEPH] {
            let flags = ephemeris | FLG_SPEED | FLG_J2000 | FLG_NONUT | FLG_NOGDEFL;
            let canonical = calc(&engine, 2451545.0, CERES, flags).unwrap();
            let alias = calc(&engine, 2451545.0, AST_OFFSET + 1, flags).unwrap();
            assert_eq!(alias.0.map(f64::to_bits), canonical.0.map(f64::to_bits));
            assert_eq!(alias.1, canonical.1);
        }
        assert_eq!(hits.get(), 4);
    }

    #[test]
    fn output_motion_reuses_the_same_successful_native_call_and_read_sequence() {
        let model = NativeModel {
            reads: RefCell::new(Vec::new()),
        };
        let nutation_calls = Cell::new(0);
        let nutation = |_: f64| {
            nutation_calls.set(nutation_calls.get() + 1);
            Ok((0.0, 0.0))
        };
        let serving = SealedServing::new(&model, &nutation);
        let flags = FLG_SPEED | FLG_J2000 | FLG_NONUT | FLG_NOGDEFL;
        let expected = serving_compute_tt(&serving, 2451545.0, MARS, flags).unwrap();
        let reads = model.reads.take();
        let nuts = nutation_calls.replace(0);
        let (actual, motion) = serving
            .with_output_motion(2451545.0, MARS, flags, || {
                serving_compute_tt(&serving, 2451545.0, MARS, flags)
            })
            .unwrap();
        assert_eq!(actual.0.map(f64::to_bits), expected.0.map(f64::to_bits));
        assert_eq!(actual.1, expected.1);
        assert_eq!(model.reads.take(), reads);
        assert_eq!(nutation_calls.get(), nuts);
        assert!(motion.is_some());
        assert!(serving.output_motion.request.get().is_none());
    }
}

#[cfg(test)]
mod owned_apparent_operand_tests {
    use super::*;

    #[test]
    fn physical_query_erases_only_equatorial_and_preserves_split_and_sidereal() {
        let producer = std::num::NonZeroU64::new(9).unwrap();
        let query = |tt, parts, flags, order, xyz| {
            owned_apparent_query(producer, tt, parts, MARS, flags, order, xyz)
        };
        let base = query(
            2451545.0,
            None,
            FLG_SPEED,
            CentralStateOrder::TargetFirst,
            false,
        );
        assert_eq!(
            base,
            query(
                2451545.0,
                None,
                FLG_SPEED | FLG_EQUATORIAL,
                CentralStateOrder::TargetFirst,
                false
            )
        );
        assert_ne!(
            base,
            query(
                2451545.0,
                None,
                FLG_SPEED | FLG_SIDEREAL,
                CentralStateOrder::TargetFirst,
                false
            )
        );
        assert_ne!(
            base,
            query(
                2451545.0,
                Some((2451545.0, 0.0)),
                FLG_SPEED,
                CentralStateOrder::TargetFirst,
                false
            )
        );
        let parts = query(
            2451545.0,
            Some((2451545.0, 0.0)),
            FLG_SPEED,
            CentralStateOrder::TargetFirst,
            false,
        );
        assert_ne!(
            parts,
            query(
                2451545.0,
                Some((2451545.0, -0.0)),
                FLG_SPEED,
                CentralStateOrder::TargetFirst,
                false
            )
        );
        assert_ne!(
            base,
            query(
                2451545.0,
                None,
                FLG_SPEED,
                CentralStateOrder::ObserverFirst,
                false
            )
        );
        assert_ne!(
            base,
            query(
                2451545.0,
                None,
                FLG_SPEED,
                CentralStateOrder::TargetFirst,
                true
            )
        );
        for flag in [
            FLG_NONUT,
            FLG_J2000,
            FLG_TRUEPOS,
            FLG_NOGDEFL,
            FLG_NOABERR,
            FLG_SPEED3,
            FLG_RADIANS,
            FLG_XYZ,
            FLG_HELCTR,
            FLG_BARYCTR,
        ] {
            assert_ne!(
                base,
                query(
                    2451545.0,
                    None,
                    FLG_SPEED | flag,
                    CentralStateOrder::TargetFirst,
                    false
                )
            );
        }
    }

    #[test]
    fn physical_payload_keeps_scalar_and_expanded_values_distinct_and_bit_exact() {
        let bits: [u64; 15] = std::array::from_fn(|i| (i as u64) | (1_u64 << 63));
        let values = bits.map(f64::from_bits);
        let (position, motion) = owned_apparent_from_words(values);
        assert_eq!(
            owned_apparent_words(position, motion).map(f64::to_bits),
            bits
        );
        assert_ne!(position[0].to_bits(), motion.0[0].hi.to_bits());
    }
}

#[cfg(test)]
mod owned_apparent_capability_tests {
    use super::*;
    use std::cell::Cell;

    #[test]
    fn public_opaque_feed_declines_reuse_and_keeps_its_next_refusal() {
        let bytes = crate::test_data::core_bytes();
        let image = LebFile::from_bytes("base_core.leb2".into(), bytes.to_vec()).unwrap();
        let provider = SingleLeb::from_file(&image).unwrap();
        let calls = Cell::new(0);
        let refuse = Cell::new(false);
        let feed = |_| {
            calls.set(calls.get() + 1);
            if refuse.get() {
                Err(Error::MathDomain { value: 9.0 })
            } else {
                Ok((0.0, 0.0))
            }
        };
        let opaque = SealedServing::new(&provider, &feed);
        assert!(opaque.owned_apparent_id.is_none());
        let flags = FLG_SWIEPH | FLG_SPEED;
        serving_compute_tt(&opaque, 2451545.0, MARS, flags).unwrap();
        let first_calls = calls.get();
        refuse.set(true);
        assert!(matches!(
            serving_compute_tt(&opaque, 2451545.0, MARS, flags | FLG_EQUATORIAL),
            Err(Error::MathDomain { value: 9.0 })
        ));
        assert!(calls.get() > first_calls);
        let native = SealedServing::new_stored(&provider, &feed);
        let trusted =
            native.with_owned_apparent_token(crate::pure_cache::new_owned_apparent_identity());
        assert!(trusted.owned_apparent_id.is_some());
        assert!(
            trusted
                .with_sampled_derivatives()
                .owned_apparent_id
                .is_none()
        );
    }
}
