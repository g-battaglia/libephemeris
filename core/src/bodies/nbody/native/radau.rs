// SPDX-License-Identifier: AGPL-3.0
//! Original eight-node left-Radau collocation for x'' = a(t, x, v).
//!
//! The acceleration is interpolated in the Newton basis at the zeros of
//! `P_8 + P_7` (mapped to `[0,1]`). Exact polynomial antiderivatives yield the
//! trial position and velocity. Implicit stages are solved by fixed-point
//! iteration; a failure is explicit and can only trigger a smaller trial step.
//! Equations: Rein & Spiegel (2015), section 2.1, equations 1–7.
//!
//! This is a numerical foundation, not a claim to reproduce IAS15's full
//! predictor, convergence policy, timestep controller or long-arc accuracy.
//! Timestep selection uses the acceleration/jerk/snap/crackle Aarseth timescale
//! discussed in Pham, Rein & Spiegel (2024), equations6–7:
//! <https://arxiv.org/html/2401.02849v1>.
//! The unchanged dimensionless budget is calibrated by7!=5040 for a harmonic
//! orbit, not by comparison output. Lower derivatives avoid mistaking noisy
//! seventh differences for physical curvature near offset close approaches.
//! No published vendor table or implementation is an input.

use super::radau_tables::{NODES, POSITION_INTEGRALS, VELOCITY_INTEGRALS};
use super::state::{AccumulatedState, NumericalError, State, Sum};

const STAGES: usize = 8;
const ENDPOINT: usize = STAGES;

/// Distinguish solver failures from the caller's original force/data error.
/// Only a numerical collocation convergence failure permits a smaller trial.
/// Even a callback error with an identical name must propagate without retry.
#[derive(Debug)]
pub(super) enum IntegrationError<E> {
    Numerical(NumericalError),
    Force(E),
}

impl<E> From<NumericalError> for IntegrationError<E> {
    fn from(error: NumericalError) -> Self {
        Self::Numerical(error)
    }
}

impl IntegrationError<NumericalError> {
    fn flatten(self) -> NumericalError {
        match self {
            Self::Numerical(error) | Self::Force(error) => error,
        }
    }
}

/// Explicit limits and dimensionless acceleration-coefficient error control.
#[derive(Clone, Copy, Debug)]
pub(super) struct Settings {
    pub(super) coefficient_tolerance: f64,
    pub(super) initial_step: f64,
    pub(super) maximum_step: f64,
    pub(super) maximum_attempts: usize,
    pub(super) maximum_iterations: usize,
}

impl Settings {
    fn checked(self) -> Result<Self, NumericalError> {
        if !self.coefficient_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.coefficient_tolerance)
            || self.coefficient_tolerance == 0.0
            || !self.initial_step.is_finite()
            || self.initial_step <= 0.0
            || !self.maximum_step.is_finite()
            || self.maximum_step <= 0.0
            || self.maximum_attempts == 0
            || self.maximum_iterations == 0
        {
            return Err(NumericalError::InvalidSettings);
        }
        Ok(self)
    }
}

/// Counts real attempts and evaluations, including rejected implicit trials.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct Statistics {
    pub(super) accepted_steps: usize,
    pub(super) rejected_steps: usize,
    pub(super) force_evaluations: usize,
}

/// Acceleration and an optional independently known local sampling timescale.
/// A positive infinite limit leaves step selection to the numerical controller.
/// The limit is evaluated at committed and converged trial states in either
/// time direction, including the endpoint omitted by left-Radau force nodes;
/// it is not a global error bound or an alternate force/propagation route.
pub(super) trait ForceProvider {
    type Error;

    fn acceleration(&mut self, time: f64, state: &State) -> Result<[f64; 3], Self::Error>;

    fn maximum_step(&mut self, _time: f64, _state: &State) -> Result<f64, Self::Error> {
        Ok(f64::INFINITY)
    }
}

impl<E, F: FnMut(f64, &State) -> Result<[f64; 3], E>> ForceProvider for F {
    type Error = E;

    fn acceleration(&mut self, time: f64, state: &State) -> Result<[f64; 3], E> {
        self(time, state)
    }
}

/// One converged step; it has not yet been accepted by the adaptive controller.
#[derive(Clone, Copy, Debug)]
pub(super) struct Trial {
    pub(super) increment: State,
    /// Harmonic-equivalent coefficient ratio `(h/tau_A85)^7/7!`, with all
    /// derivatives evaluated consistently at the greatest-acceleration stage.
    pub(super) coefficient_ratio: f64,
    /// Smallest physical sampling interval at the converged nodes and endpoint.
    maximum_step: f64,
    /// Diagnostic-only polynomial; absent from every runtime build.
    #[cfg(test)]
    acceleration_coefficients: [[f64; 3]; STAGES],
}

/// Accepted numerical polynomial data for an independent IVP error reference.
///
/// Time/step are relative days; state is SSB AU/AU-day for astronomical probes.
/// Coefficients multiply the Newton basis at the stored normalized Radau nodes
/// and have acceleration units. A reference must integrate that basis itself,
/// retain knot/time/roundoff differences and prove the complete force defect.
/// Capturing these words alone does not certify a trajectory.
#[cfg(test)]
#[derive(Clone, Copy, Debug)]
pub(super) struct AcceptedStep {
    pub(super) start_time: f64,
    pub(super) end_time: f64,
    pub(super) step: f64,
    pub(super) initial: State,
    pub(super) endpoint: State,
    pub(super) acceleration_coefficients: [[f64; 3]; STAGES],
}

#[cfg(test)]
const MAX_CAPTURED_STEPS: usize = 8192;

/// Newton divided differences for the acceleration at all eight nodes.
///
/// Subtract the common initial acceleration once before the recurrence; the
/// nonconstant polynomial then interpolates small acceleration changes without
/// repeatedly subtracting a large common component at each divided difference.
fn coefficients(samples: &[[f64; 3]; STAGES]) -> [[f64; 3]; STAGES] {
    let mut result = std::array::from_fn(|stage| {
        std::array::from_fn(|axis| samples[stage][axis] - samples[0][axis])
    });
    for width in 1..STAGES {
        for stage in (width..STAGES).rev() {
            let previous = result[stage - 1];
            let spacing = NODES[stage] - NODES[stage - width];
            for (value, preceding) in result[stage].iter_mut().zip(previous) {
                *value = (*value - preceding) / spacing;
            }
        }
    }
    result[0] = samples[0];
    result
}

/// Twice/once integrate the interpolating acceleration to one stage or endpoint.
fn increment(
    initial: State,
    coefficients: &[[f64; 3]; STAGES],
    step: f64,
    stage: usize,
) -> Result<State, NumericalError> {
    let fraction = if stage == ENDPOINT { 1.0 } else { NODES[stage] };
    let mut position = [0.0; 3];
    let mut velocity = [0.0; 3];
    for axis in 0..3 {
        let mut distance = Sum::new(initial.velocity[axis] * (step * fraction));
        let mut speed = Sum::default();
        for (order, coefficient) in coefficients.iter().enumerate() {
            // Scale each integral before multiplying by h a second time, rather
            // than squaring h (which can overflow even for a finite increment).
            distance.add(((coefficient[axis] * POSITION_INTEGRALS[stage][order]) * step) * step);
            speed.add((coefficient[axis] * VELOCITY_INTEGRALS[stage][order]) * step);
        }
        position[axis] = distance.value();
        velocity[axis] = speed.value();
    }
    State { position, velocity }.checked()
}

/// Solve one implicit collocation step without changing its initial state.
///
/// Time and step are relative to the caller's fixed epoch, not a large JD.
/// The callback sees actual finite stage states and may depend on velocity.
/// Callback errors propagate unchanged: they are not evidence that a smaller
/// step or a different force model can serve the request.
#[cfg(test)]
pub(super) fn trial(
    initial: AccumulatedState,
    time: f64,
    step: f64,
    maximum_iterations: usize,
    force: &mut impl FnMut(f64, &State) -> Result<[f64; 3], NumericalError>,
    statistics: &mut Statistics,
) -> Result<Trial, NumericalError> {
    trial_with_error(initial, time, step, maximum_iterations, force, statistics)
        .map_err(IntegrationError::flatten)
}

/// The same collocation solve with the caller's typed force failure retained.
#[cfg(test)]
pub(super) fn trial_with_error<E>(
    initial: AccumulatedState,
    time: f64,
    step: f64,
    maximum_iterations: usize,
    force: &mut impl FnMut(f64, &State) -> Result<[f64; 3], E>,
    statistics: &mut Statistics,
) -> Result<Trial, IntegrationError<E>> {
    trial_with_provider(initial, time, step, maximum_iterations, force, statistics)
}

fn trial_with_provider<P: ForceProvider>(
    initial: AccumulatedState,
    time: f64,
    step: f64,
    maximum_iterations: usize,
    provider: &mut P,
    statistics: &mut Statistics,
) -> Result<Trial, IntegrationError<P::Error>> {
    if !time.is_finite() || !step.is_finite() || !(time + step).is_finite() {
        return Err(NumericalError::NonFinite.into());
    }
    if step == 0.0 {
        return Err(NumericalError::InvalidStep.into());
    }
    if time + step == time {
        return Err(NumericalError::StepUnderflow.into());
    }
    let state = initial.state().checked()?;
    let first = checked_force(time, &state, provider, statistics)?;
    let mut samples = [first; STAGES];
    for _ in 0..maximum_iterations {
        let polynomial = coefficients(&samples);
        let mut next = samples;
        let mut change = 0.0_f64;
        let mut scale = 0.0_f64;
        for stage in 1..STAGES {
            let prediction = initial
                .advanced(increment(state, &polynomial, step, stage)?)?
                .state();
            let stage_time = time + step * NODES[stage];
            next[stage] = checked_force(stage_time, &prediction, provider, statistics)?;
        }
        for stage in 0..STAGES {
            for axis in 0..3 {
                change = change.max((next[stage][axis] - samples[stage][axis]).abs());
                scale = scale
                    .max(next[stage][axis].abs())
                    .max(samples[stage][axis].abs());
            }
        }
        samples = next;
        // A roundoff-level stage residual, not a fitted trajectory tolerance.
        // A separate coefficient budget controls the accepted step length.
        if change <= 32.0 * f64::EPSILON * scale {
            let polynomial = coefficients(&samples);
            let ratio = timestep_ratio(&polynomial);
            if !ratio.is_finite() {
                return Err(NumericalError::NonFinite.into());
            }
            // Left-Radau has no endpoint force node. Check its endpoint too:
            // a weak narrow encounter can otherwise occur wholly between the
            // last force sample and the committed state. These are converged
            // polynomial predictions, never unconverged iteration guesses.
            let mut maximum_step = f64::INFINITY;
            for (stage, fraction) in NODES
                .into_iter()
                .enumerate()
                .chain(std::iter::once((ENDPOINT, 1.0)))
            {
                let prediction = initial
                    .advanced(increment(state, &polynomial, step, stage)?)?
                    .state();
                maximum_step = maximum_step.min(checked_limit(
                    provider
                        .maximum_step(time + step * fraction, &prediction)
                        .map_err(IntegrationError::Force)?,
                )?);
            }
            return Ok(Trial {
                increment: increment(state, &polynomial, step, ENDPOINT)?,
                coefficient_ratio: ratio,
                maximum_step,
                #[cfg(test)]
                acceleration_coefficients: polynomial,
            });
        }
    }
    Err(NumericalError::IterationLimit.into())
}

/// Evaluate the Newton polynomial and its first three normalized-time
/// derivatives. These all have acceleration units; dividing by h, h², h³
/// would give physical jerk, snap and crackle. Nested product differentiation
/// avoids a high-order finite difference or a rounded epoch.
fn scaled_derivatives(polynomial: &[[f64; 3]; STAGES], fraction: f64) -> [[f64; 3]; 4] {
    let mut value = polynomial[STAGES - 1];
    let mut first = [0.0_f64; 3];
    let mut second = [0.0_f64; 3];
    let mut third = [0.0_f64; 3];
    for order in (0..STAGES - 1).rev() {
        let offset = fraction - NODES[order];
        for axis in 0..3 {
            third[axis] = third[axis].mul_add(offset, 3.0 * second[axis]);
            second[axis] = second[axis].mul_add(offset, 2.0 * first[axis]);
            first[axis] = first[axis].mul_add(offset, value[axis]);
            value[axis] = value[axis].mul_add(offset, polynomial[order][axis]);
        }
    }
    [value, first, second, third]
}

/// Choose one consistent physical anchor for all derivative norms. Greatest
/// acceleration over collocation points/endpoint is rotation-invariant and
/// resolves pericentre; no separate time is chosen for each derivative.
fn derivative_norms(polynomial: &[[f64; 3]; STAGES]) -> [f64; 4] {
    fn magnitude(vector: [f64; 3]) -> f64 {
        vector[0].hypot(vector[1]).hypot(vector[2])
    }
    let mut selected = [[0.0; 3]; 4];
    let mut acceleration = 0.0_f64;
    for fraction in NODES.into_iter().chain(std::iter::once(1.0)) {
        let derivatives = scaled_derivatives(polynomial, fraction);
        let norm = magnitude(derivatives[0]);
        if norm > acceleration {
            acceleration = norm;
            selected = derivatives;
        }
    }
    selected.map(magnitude)
}

/// A85: tau²=(|a||s|+|j|²)/(|j||c|+|s|²), using norms, not dot
/// products. The extra derivative retains a finite timescale at a genuine
/// acceleration zero with nonzero jerk (including velocity-dependent forces).
/// It remains a two-body timescale criterion, not a rigorous arbitrary-force
/// global error bound; independent trajectory qualification is still required.
fn timestep_ratio(polynomial: &[[f64; 3]; STAGES]) -> f64 {
    a85_ratio(derivative_norms(polynomial))
}

fn a85_ratio(norms: [f64; 4]) -> f64 {
    // A common scale cancels algebraically and keeps squared norms finite.
    let scale = norms.into_iter().fold(0.0_f64, f64::max);
    if scale == 0.0 {
        return 0.0;
    }
    let [acceleration, jerk, snap, crackle] = norms.map(|value| value / scale);
    let numerator = jerk * crackle + snap * snap;
    let denominator = acceleration * snap + jerk * jerk;
    if numerator == 0.0 {
        // Constant acceleration/jerk is integrated exactly by this collocation
        // polynomial. This is a removable dynamical-timescale limit.
        return 0.0;
    }
    let phase_squared = numerator / denominator;
    phase_squared.sqrt().powi(7) / 5040.0
}

fn checked_force<P: ForceProvider>(
    time: f64,
    state: &State,
    provider: &mut P,
    statistics: &mut Statistics,
) -> Result<[f64; 3], IntegrationError<P::Error>> {
    if !time.is_finite() {
        return Err(NumericalError::NonFinite.into());
    }
    statistics.force_evaluations += 1;
    let result = provider
        .acceleration(time, state)
        .map_err(IntegrationError::Force)?;
    if result.into_iter().all(f64::is_finite) {
        Ok(result)
    } else {
        Err(NumericalError::NonFinite.into())
    }
}

fn checked_limit(limit: f64) -> Result<f64, NumericalError> {
    if limit.is_nan() || limit < 0.0 {
        Err(NumericalError::InvalidSettings)
    } else {
        Ok(limit)
    }
}

/// Retained per-trajectory state; no globals, ephemeris cache or callback owner.
#[derive(Clone)]
pub(super) struct Radau {
    state: AccumulatedState,
    time: f64,
    next_step: f64,
    settings: Settings,
    #[cfg(test)]
    captured_steps: Option<Vec<AcceptedStep>>,
    #[cfg(test)]
    capture_overflow: bool,
}

impl Radau {
    pub(super) fn new(state: State, time: f64, settings: Settings) -> Result<Self, NumericalError> {
        if !time.is_finite() {
            return Err(NumericalError::NonFinite);
        }
        let settings = settings.checked()?;
        Ok(Self {
            state: AccumulatedState::new(state)?,
            time,
            next_step: settings.initial_step.min(settings.maximum_step),
            settings,
            #[cfg(test)]
            captured_steps: None,
            #[cfg(test)]
            capture_overflow: false,
        })
    }

    pub(super) fn state(&self) -> State {
        self.state.state()
    }

    pub(super) fn time(&self) -> f64 {
        self.time
    }

    /// Enable a bounded diagnostic capture without changing force evaluations.
    /// Only test builds contain this storage; disabled captures allocate nothing.
    #[cfg(test)]
    pub(super) fn enable_step_capture(&mut self) {
        self.captured_steps = Some(Vec::new());
        self.capture_overflow = false;
    }

    /// Borrow accepted records and whether the bounded capture stayed complete.
    /// An overflow truncates diagnostics only; numerical integration is unchanged.
    #[cfg(test)]
    pub(super) fn step_capture(&self) -> Option<(&[AcceptedStep], bool)> {
        self.captured_steps
            .as_deref()
            .map(|steps| (steps, !self.capture_overflow))
    }

    /// Advance to the exact represented endpoint, retaining only accepted steps.
    ///
    /// The harmonic-equivalent timescale ratio scales as h^7; its budget and
    /// seventh-root control are dimensionless. Fourfold growth bounds surprises.
    /// A nonconvergent trial may be retried at half length. No force/data error
    /// is retried as success, and no analytic output replaces a failed solve.
    pub(super) fn advance(
        &mut self,
        target: f64,
        force: &mut impl FnMut(f64, &State) -> Result<[f64; 3], NumericalError>,
    ) -> Result<Statistics, NumericalError> {
        self.advance_with_error(target, force)
            .map_err(IntegrationError::flatten)
    }

    /// Advance through exactly the same adaptive numerical route, retaining
    /// arbitrary typed force/data errors without converting them to convergence
    /// failures, replacing sources, or publishing a fabricated endpoint state.
    pub(super) fn advance_with_error<E>(
        &mut self,
        target: f64,
        force: &mut impl FnMut(f64, &State) -> Result<[f64; 3], E>,
    ) -> Result<Statistics, IntegrationError<E>> {
        self.advance_with_provider(target, force)
    }

    /// The same adaptive solve, capped by the force provider's physical local
    /// sampling interval. Provider failures propagate once, without retry.
    pub(super) fn advance_with_provider<P: ForceProvider>(
        &mut self,
        target: f64,
        provider: &mut P,
    ) -> Result<Statistics, IntegrationError<P::Error>> {
        if !target.is_finite() || !(target - self.time).is_finite() {
            return Err(NumericalError::NonFinite.into());
        }
        let mut statistics = Statistics::default();
        for _ in 0..self.settings.maximum_attempts {
            if self.time == target {
                return Ok(statistics);
            }
            let remaining = target - self.time;
            let local_limit = checked_limit(
                provider
                    .maximum_step(self.time, &self.state.state())
                    .map_err(IntegrationError::Force)?,
            )?;
            let length = self
                .next_step
                .min(remaining.abs())
                .min(self.settings.maximum_step)
                .min(local_limit);
            let step = length.copysign(remaining);
            if step == 0.0 || self.time + step == self.time {
                return Err(NumericalError::StepUnderflow.into());
            }
            let candidate = match trial_with_provider(
                self.state,
                self.time,
                step,
                self.settings.maximum_iterations,
                provider,
                &mut statistics,
            ) {
                Ok(candidate) => candidate,
                Err(IntegrationError::Numerical(NumericalError::IterationLimit)) => {
                    statistics.rejected_steps += 1;
                    self.next_step = length * 0.5;
                    continue;
                }
                Err(error) => return Err(error),
            };
            if length > candidate.maximum_step {
                statistics.rejected_steps += 1;
                self.next_step = candidate.maximum_step;
                continue;
            }
            let factor = if candidate.coefficient_ratio == 0.0 {
                4.0
            } else {
                (0.9 * (self.settings.coefficient_tolerance / candidate.coefficient_ratio)
                    .powf(1.0 / 7.0))
                .clamp(0.25, 4.0)
            };
            self.next_step = (length * factor).min(self.settings.maximum_step);
            if candidate.coefficient_ratio > self.settings.coefficient_tolerance {
                statistics.rejected_steps += 1;
                continue;
            }
            #[cfg(test)]
            let capture_start = (self.time, self.state.state());
            self.state = self.state.advanced(candidate.increment)?;
            self.time = if length == remaining.abs() {
                target
            } else {
                self.time + step
            };
            #[cfg(test)]
            if let Some(steps) = &mut self.captured_steps {
                if steps.len() < MAX_CAPTURED_STEPS {
                    steps.push(AcceptedStep {
                        start_time: capture_start.0,
                        end_time: self.time,
                        step,
                        initial: capture_start.1,
                        endpoint: self.state.state(),
                        acceleration_coefficients: candidate.acceleration_coefficients,
                    });
                } else {
                    self.capture_overflow = true;
                }
            }
            statistics.accepted_steps += 1;
        }
        if self.time == target {
            Ok(statistics)
        } else {
            Err(NumericalError::StepLimit.into())
        }
    }
}

#[cfg(test)]
#[path = "radau_error_tests.rs"]
mod error_tests;

#[cfg(test)]
#[path = "radau_timescale_tests.rs"]
mod timescale_tests;

#[cfg(test)]
#[path = "radau_capture_tests.rs"]
mod capture_tests;
