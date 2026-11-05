//! Astrometric range and its TT derivative before direction corrections.
//!
//! A physical ray `r = p_target(t_emit) - p_observer(t_obs)` has range
//! `R = |r|` and rate `R' = r·r'/R`. Velocities supplied here are derivatives
//! with respect to the outer TT day. Aberration, gravitational deflection and
//! rigid output frames change the direction, not this physical range.

use super::compensated_dot3;

/// Range in AU and its analytic derivative in AU per TT day.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct RangeState {
    /// Euclidean length of the uncorrected observer-relative ICRS ray, AU.
    pub(super) distance: f64,
    /// `ray·velocity/|ray|`, AU per TT day; zero at a degenerate origin.
    pub(super) rate: f64,
}

impl RangeState {
    /// Differentiate one finite Cartesian ray without differencing radii.
    ///
    /// `ray` is in AU and `velocity` is its full derivative per outer TT day,
    /// including observer motion, the provider clock and emission-time
    /// Jacobian. Compensated products limit cancellation in a nearly tangential
    /// state. A zero ray uses the public zero-origin convention. The caller
    /// validates finite inputs; no provider or force model is selected here.
    pub(super) fn from_cartesian(ray: [f64; 3], velocity: [f64; 3]) -> Self {
        let distance = compensated_dot3(ray, ray).sqrt();
        let rate = if distance == 0.0 {
            0.0
        } else {
            compensated_dot3(ray, velocity) / distance
        };
        Self { distance, rate }
    }

    /// Differentiate the next retarded target coordinate epoch.
    ///
    /// For `t_emit = TDB(TT) - R/c`, the emission derivative is
    /// `dTDB/dTT - R'/c`. `observation_rate` is the provider coordinate clock
    /// per TT day, and `light_speed` is positive AU per coordinate day.
    /// Applying this at every actual light-time sample differentiates both
    /// finite iteration maps and converged light paths with the same formula.
    pub(super) fn next_emission_rate(self, observation_rate: f64, light_speed: f64) -> f64 {
        observation_rate - self.rate / light_speed
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn radial_derivative_preserves_a_small_tangential_residual() {
        let state = RangeState::from_cartesian([1.0; 3], [1e16, 1.0, -1e16]);
        assert_eq!(state.distance, 3.0_f64.sqrt());
        assert_eq!(state.rate, 1.0 / 3.0_f64.sqrt());
        assert_eq!(RangeState::from_cartesian([0.0; 3], [1.0; 3]).rate, 0.0);
    }

    #[test]
    fn range_rate_is_invariant_under_an_orthogonal_coordinate_rotation() {
        let ray = [3.0, -4.0, 12.0];
        let velocity = [0.003, 0.009, -0.001];
        let state = RangeState::from_cartesian(ray, velocity);
        let rotate = |v: [f64; 3]| [-v[1], v[0], v[2]];
        assert_eq!(
            state,
            RangeState::from_cartesian(rotate(ray), rotate(velocity))
        );
    }

    #[test]
    fn converged_light_time_rate_matches_the_closed_linear_solution() {
        let speed = 173.1446326846693;
        let target_position = 7.0;
        let observer_position = 2.0;
        let target_velocity = 0.01;
        let observer_velocity = 0.003;
        let clock = 1.0000000003;
        let mut emission_rate = clock;
        let mut delay = 0.0;
        let mut state = RangeState::from_cartesian([0.0; 3], [0.0; 3]);
        for _ in 0..10 {
            state = RangeState::from_cartesian(
                [
                    target_position - target_velocity * delay - observer_position,
                    0.0,
                    0.0,
                ],
                [
                    target_velocity * emission_rate - observer_velocity * clock,
                    0.0,
                    0.0,
                ],
            );
            delay = state.distance / speed;
            emission_rate = state.next_emission_rate(clock, speed);
        }
        let exact = clock * (target_velocity - observer_velocity) / (1.0 + target_velocity / speed);
        assert!((state.rate - exact).abs() <= 2.0 * f64::EPSILON * exact.abs());
    }
}
