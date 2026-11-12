// SPDX-License-Identifier: AGPL-3.0
//! Retained raw-SSB test-particle propagation over selected force resources.
//!
//! The fixed origin remains two-part TDB JD; solver times are relative days.
//! Data/epoch/pole failures retain their original category. Each requested
//! endpoint is transactional: an unsuccessful request may update immutable
//! coefficient caches, but never publishes a partially advanced solver state.
//! The orientation policy is an explicit independently defined function, not
//! an inferred legacy default or a per-input fallback.

use super::catalogue::HarmonicPoles;
use super::ephemeris::Ephemeris;
use super::evaluator::{Evaluator, ForceError};
use super::radau::{ForceProvider, IntegrationError, Radau, Settings};
use super::state::{NumericalError, State};
use crate::errors::Error;
use crate::types::TwoPartEpoch;
use std::sync::Arc;

/// Pure physical orientation policy evaluated at the actual two-part TDB JD.
pub(super) type PoleProvider = fn([f64; 2]) -> Result<HarmonicPoles, Error>;
pub(super) type PropagationError = IntegrationError<ForceError>;

/// One target's owned integration state and reusable force workspace.
pub(super) struct ForcedTrajectory {
    origin: TwoPartEpoch,
    integrator: Radau,
    evaluator: Evaluator,
    poles: PoleProvider,
    requested_step: f64,
}

/// One borrow supplies both the physical sampling limit and complete force.
/// Epoch splitting and data failures follow the same typed boundary in each.
struct TrajectoryForce<'a> {
    origin: TwoPartEpoch,
    evaluator: &'a mut Evaluator,
    poles: PoleProvider,
}

impl ForceProvider for TrajectoryForce<'_> {
    type Error = ForceError;

    fn acceleration(&mut self, relative: f64, state: &State) -> Result<[f64; 3], ForceError> {
        let epoch = self
            .origin
            .add_offset(relative)
            .map_err(ForceError::Frame)?;
        let epoch = [epoch.hi, epoch.lo];
        let orientation = (self.poles)(epoch).map_err(ForceError::Frame)?;
        self.evaluator.acceleration(epoch, *state, orientation)
    }

    fn maximum_step(&mut self, relative: f64, state: &State) -> Result<f64, ForceError> {
        let epoch = self
            .origin
            .add_offset(relative)
            .map_err(ForceError::Frame)?;
        self.evaluator.maximum_step([epoch.hi, epoch.lo], *state)
    }
}

impl ForcedTrajectory {
    /// Construct from raw SSB ICRF AU/AU-day state. `origin` is TDB JD;
    /// the optional step is finite days, negative accepted as a magnitude.
    /// A zero step is valid only while no nonzero interval is requested.
    /// No mutable trajectory state is shared with another target/worker.
    pub(super) fn new(
        ephemeris: Arc<Ephemeris>,
        initial: State,
        origin: [f64; 2],
        step: Option<f64>,
        nongravitational: Option<[f64; 3]>,
        exclude: Option<i32>,
        poles: PoleProvider,
    ) -> Result<Self, PropagationError> {
        initial.checked()?;
        if !origin.into_iter().all(f64::is_finite) || step.is_some_and(|value| !value.is_finite()) {
            return Err(NumericalError::NonFinite.into());
        }
        let origin = TwoPartEpoch {
            hi: origin[0],
            lo: origin[1],
        }
        .add_offset(0.0)
        .map_err(frame_error)?;
        let initial_poles = poles([origin.hi, origin.lo]).map_err(frame_error)?;
        let evaluator = Evaluator::new(ephemeris, initial_poles, nongravitational, exclude)
            .map_err(IntegrationError::Force)?;
        let requested_step = step.unwrap_or(0.001).abs();
        let initial_step = if requested_step == 0.0 {
            0.001
        } else {
            requested_step
        };
        let integrator = Radau::new(
            initial,
            0.0,
            Settings {
                coefficient_tolerance: 1e-9,
                initial_step,
                maximum_step: 32.0,
                maximum_attempts: 20_000_000,
                maximum_iterations: 64,
            },
        )?;
        Ok(Self {
            origin,
            integrator,
            evaluator,
            poles,
            requested_step,
        })
    }

    pub(super) fn state(&self) -> State {
        self.integrator.state()
    }
    pub(super) fn time(&self) -> f64 {
        self.integrator.time()
    }

    /// Evaluate the retained complete force at the committed endpoint.
    ///
    /// Returns SSB ICRF AU/day² without advancing time or changing the state.
    /// Selected data/pole/domain errors propagate exactly as during integration;
    /// only reusable force scratch and immutable coefficient caches may change.
    pub(super) fn acceleration(&mut self) -> Result<[f64; 3], PropagationError> {
        let epoch = self.origin.add_offset(self.time()).map_err(frame_error)?;
        let epoch = [epoch.hi, epoch.lo];
        let poles = (self.poles)(epoch).map_err(frame_error)?;
        let state = self.state();
        self.evaluator
            .acceleration(epoch, state, poles)
            .map_err(IntegrationError::Force)
    }

    /// Reach a finite relative TDB-day endpoint, in either direction. Output is
    /// raw SSB ICRF AU/AU-day, without TT conversion, heliocentric translation,
    /// apparent-place correction, sorting, duplication removal or extrapolation.
    /// A failed request retains the previous committed state/time exactly.
    pub(super) fn advance(&mut self, target: f64) -> Result<State, PropagationError> {
        if !target.is_finite() || !(target - self.time()).is_finite() {
            return Err(NumericalError::NonFinite.into());
        }
        if target == self.time() {
            return Ok(self.state());
        }
        if self.requested_step == 0.0 {
            return Err(NumericalError::InvalidStep.into());
        }
        let mut candidate = self.integrator.clone();
        candidate.advance_with_provider(
            target,
            &mut TrajectoryForce {
                origin: self.origin,
                poles: self.poles,
                evaluator: &mut self.evaluator,
            },
        )?;
        self.integrator = candidate;
        Ok(self.state())
    }
}

fn frame_error(error: Error) -> PropagationError {
    IntegrationError::Force(ForceError::Frame(error))
}

#[cfg(test)]
#[path = "trajectory_tests.rs"]
mod tests;

#[cfg(all(test, feature = "nbody"))]
#[path = "trajectory_data_tests.rs"]
mod data_tests;
