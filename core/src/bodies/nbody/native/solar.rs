// SPDX-License-Identifier: AGPL-3.0
//! Sun + massless target evolution with independently authored numerical maps.
//!
//! The supported public problem has one stationary massive primary and one
//! massless target, not a mutable general N-body system. In this limit the
//! Wisdom-Holman interaction, central jump and center-of-mass maps vanish.
//! Mercurius has no planet-pair encounter map; moving a zero jump/interaction
//! term between TRACE's sub-Hamiltonians leaves the same Kepler flow, including
//! its central-pericentre regime. Their applicable maps therefore share the
//! exact conic flow, not an unrelated numerical method disguised by a name.
//! Leapfrog retains its distinct symmetric drift-kick-drift splitting, and
//! the adaptive selector retains its original Gauss-Radau numerical solve.
//!
//! References: Wisdom & Holman (1991); Rein et al. (2019), hybrid map;
//! TRACE (2024), section3 and appendixA. These are mathematical reductions
//! of the exposed problem, not claims to contain those upstream programs.

use super::conic::{Elements, initial_state};
use super::newtonian::{PointMass, acceleration};
use super::radau::{Radau, Settings};
use super::state::{AccumulatedState, NumericalError, State};

/// Public method labels translated once into native mathematical capabilities.
#[derive(Clone, Copy, Debug)]
pub(super) enum Method {
    Adaptive,
    WisdomHolman,
    Leapfrog,
    Mercurius,
    Trace,
}

/// Limits are explicit and independent of target identity or comparison output.
const MAXIMUM_STEPS: usize = 20_000_000;

/// One owned Solar-system test-particle trajectory in its supplied inertial plane.
pub(super) struct SolarTrajectory {
    initial: State,
    gm: f64,
    method: Method,
    fixed_step: f64,
    state: AccumulatedState,
    time: f64,
    adaptive: Option<Radau>,
}

impl SolarTrajectory {
    /// Construct from independently validated elements; steps use the caller's
    /// time unit. The optional step is an initial adaptive step or a fixed map
    /// step, not an accuracy promise. Zero is admissible only for no evolution.
    pub(super) fn new(
        elements: Elements,
        gm: f64,
        method: Method,
        step: Option<f64>,
    ) -> Result<Self, NumericalError> {
        if step.is_some_and(|value| !value.is_finite()) {
            return Err(NumericalError::NonFinite);
        }
        let initial_state = initial_state(elements, gm)?;
        let selected_step = step.unwrap_or(0.001).abs();
        let adaptive = if matches!(method, Method::Adaptive) {
            // A zero requested step is kept in fixed_step for the advance
            // guard; constructor-only/zero-duration calls still return PV.
            Some(Radau::new(
                initial_state,
                0.0,
                Settings {
                    coefficient_tolerance: 1e-9,
                    initial_step: if selected_step == 0.0 {
                        0.001
                    } else {
                        selected_step
                    },
                    maximum_step: 32.0,
                    maximum_attempts: MAXIMUM_STEPS,
                    maximum_iterations: 64,
                },
            )?)
        } else {
            None
        };
        Ok(Self {
            initial: initial_state,
            gm,
            method,
            fixed_step: selected_step,
            state: AccumulatedState::new(initial_state)?,
            time: 0.0,
            adaptive,
        })
    }

    pub(super) fn state(&self) -> State {
        self.state.state()
    }

    #[cfg(feature = "nbody")]
    pub(super) fn time(&self) -> f64 {
        self.time
    }

    /// Reach an exact finite relative endpoint; failure never invents an
    /// endpoint state or exchanges the requested force/method for another.
    pub(super) fn advance(&mut self, target: f64) -> Result<State, NumericalError> {
        if !target.is_finite() || !(target - self.time).is_finite() {
            return Err(NumericalError::NonFinite);
        }
        if target == self.time {
            return Ok(self.state());
        }
        if self.fixed_step == 0.0 {
            return Err(NumericalError::InvalidStep);
        }
        let state = match self.method {
            Method::Adaptive => {
                let source = PointMass {
                    body: 10,
                    position: [0.0; 3],
                    gravitational_parameter: self.gm,
                };
                let mut adaptive = self
                    .adaptive
                    .as_ref()
                    .expect("adaptive state follows method")
                    .clone();
                adaptive.advance(target, &mut |_, state| {
                    acceleration(state.position, &[source], None)
                })?;
                let state = adaptive.state();
                self.adaptive = Some(adaptive);
                state
            }
            Method::Leapfrog => {
                let state = self.leapfrog(target)?;
                self.state = state;
                self.time = target;
                return Ok(state.state());
            }
            // These are proven reductions of three separately named published
            // maps in this fixed, two-body public problem. No many-body or
            // encounter capability is advertised by this private solver.
            Method::WisdomHolman => super::kepler::advance(self.initial, self.gm, target)?,
            Method::Mercurius => super::kepler::advance(self.initial, self.gm, target)?,
            Method::Trace => super::kepler::advance(self.initial, self.gm, target)?,
        };
        self.state = AccumulatedState::new(state)?;
        self.time = target;
        Ok(state)
    }

    fn leapfrog(&self, target: f64) -> Result<AccumulatedState, NumericalError> {
        let source = PointMass {
            body: 10,
            position: [0.0; 3],
            gravitational_parameter: self.gm,
        };
        let mut state = self.state;
        let mut time = self.time;
        for _ in 0..MAXIMUM_STEPS {
            if time == target {
                return Ok(state);
            }
            let remaining = target - time;
            let step = self.fixed_step.min(remaining.abs()).copysign(remaining);
            if time + step == time {
                return Err(NumericalError::StepUnderflow);
            }
            let current = state.state();
            let first_drift = State {
                position: current.velocity.map(|value| value * (step / 2.0)),
                velocity: [0.0; 3],
            };
            let midpoint = state.advanced(first_drift)?;
            let force = acceleration(midpoint.state().position, &[source], None)?;
            let kicked = midpoint.advanced(State {
                position: [0.0; 3],
                velocity: force.map(|value| value * step),
            })?;
            let final_drift = State {
                position: kicked.state().velocity.map(|value| value * (step / 2.0)),
                velocity: [0.0; 3],
            };
            state = kicked.advanced(final_drift)?;
            time = if step.abs() == remaining.abs() {
                target
            } else {
                time + step
            };
        }
        if time == target {
            Ok(state)
        } else {
            Err(NumericalError::StepLimit)
        }
    }
}

#[cfg(test)]
#[path = "solar_tests.rs"]
mod tests;
