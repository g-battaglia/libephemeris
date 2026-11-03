// SPDX-License-Identifier: AGPL-3.0
//! Private finite Cartesian states and compensated state accumulation.
//!
//! Positions and velocities have caller-selected consistent units; numerical
//! integration adds no frame or timescale conversion. Astronomical callers will
//! select SSB ICRF AU, AU/day and relative TDB days at their own boundary.

/// Failures of the independently derived numerical kernels.
///
/// These private errors are not yet mapped into the public compatibility layer.
/// A rejected numerical step never changes its initial state.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum NumericalError {
    /// A supplied or derived scalar is not finite.
    NonFinite,
    /// A requested integration step is zero.
    InvalidStep,
    /// Numerical settings do not define positive finite bounds.
    InvalidSettings,
    /// A finite conic description has no admitted physical orbit.
    InvalidOrbit(&'static str),
    /// A nonexcluded massive source coincides with the target.
    SingularForce { body: i32 },
    /// A source gravitational parameter is negative.
    InvalidMass { body: i32 },
    /// The selected source list contains a repeated physical identifier.
    DuplicateSource { body: i32 },
    /// A configured physical force requires an unavailable source state.
    MissingSource { body: i32 },
    /// The implicit collocation equations did not converge.
    IterationLimit,
    /// The configured bound on attempted steps was reached.
    StepLimit,
    /// A smaller step can no longer advance the represented time.
    StepUnderflow,
}

/// Finite Cartesian position and velocity in one explicitly selected frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct State {
    pub(super) position: [f64; 3],
    pub(super) velocity: [f64; 3],
}

impl State {
    /// Refuse nonfinite states instead of allowing a solver to report success.
    pub(super) fn checked(self) -> Result<Self, NumericalError> {
        if self
            .position
            .into_iter()
            .chain(self.velocity)
            .all(f64::is_finite)
        {
            Ok(self)
        } else {
            Err(NumericalError::NonFinite)
        }
    }
}

/// Neumaier summation, retaining low terms when a new term exceeds the sum.
///
/// This is the elementary compensated-addition identity, not a vendor-specific
/// summation schedule. Both the sum and correction belong to this value.
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct Sum {
    high: f64,
    low: f64,
}

impl Sum {
    pub(super) fn new(value: f64) -> Self {
        Self {
            high: value,
            low: 0.0,
        }
    }

    pub(super) fn add(&mut self, value: f64) {
        let next = self.high + value;
        self.low += if self.high.abs() >= value.abs() {
            (self.high - next) + value
        } else {
            (value - next) + self.high
        };
        self.high = next;
    }

    pub(super) fn value(self) -> f64 {
        self.high + self.low
    }

    /// Normalize the retained sum and its error-free residual without dropping
    /// the correction before a cancellation-sensitive physical invariant uses it.
    pub(super) fn parts(self) -> [f64; 2] {
        let high = self.high + self.low;
        let recovered_low = high - self.high;
        let low = (self.high - (high - recovered_low)) + (self.low - recovered_low);
        [high, low]
    }
}

/// Running state with persistent summation compensation between accepted steps.
#[derive(Clone, Copy, Debug)]
pub(super) struct AccumulatedState {
    position: [Sum; 3],
    velocity: [Sum; 3],
}

impl AccumulatedState {
    pub(super) fn new(state: State) -> Result<Self, NumericalError> {
        let state = state.checked()?;
        Ok(Self {
            position: state.position.map(Sum::new),
            velocity: state.velocity.map(Sum::new),
        })
    }

    pub(super) fn state(self) -> State {
        State {
            position: self.position.map(Sum::value),
            velocity: self.velocity.map(Sum::value),
        }
    }

    /// Apply a complete trial increment to a copy; callers commit only success.
    pub(super) fn advanced(mut self, increment: State) -> Result<Self, NumericalError> {
        increment.checked()?;
        for axis in 0..3 {
            self.position[axis].add(increment.position[axis]);
            self.velocity[axis].add(increment.velocity[axis]);
        }
        self.state().checked()?;
        Ok(self)
    }
}
