// SPDX-License-Identifier: AGPL-3.0
//! Prepared Solar-map checks; execute with the final integration qualification.

use super::*;

fn elements() -> Elements {
    Elements {
        semimajor_axis: 1.0,
        eccentricity: 0.0,
        inclination: 0.0,
        ascending_node: 0.0,
        periapsis_argument: 0.0,
        mean_anomaly: 0.0,
    }
}

#[test]
fn every_exposed_method_preserves_its_zero_duration_state() {
    for method in [
        Method::Adaptive,
        Method::WisdomHolman,
        Method::Leapfrog,
        Method::Mercurius,
        Method::Trace,
    ] {
        let mut trajectory = SolarTrajectory::new(elements(), 1.0, method, Some(0.0)).unwrap();
        let initial = trajectory.state();
        assert_eq!(trajectory.advance(0.0).unwrap(), initial);
        assert_eq!(trajectory.advance(1.0), Err(NumericalError::InvalidStep));
        assert_eq!(trajectory.state(), initial);
    }
}

#[test]
fn reduced_hamiltonian_maps_agree_with_the_exact_circular_flow() {
    for method in [Method::WisdomHolman, Method::Mercurius, Method::Trace] {
        let mut trajectory = SolarTrajectory::new(elements(), 1.0, method, Some(0.3)).unwrap();
        for target in [0.5_f64, -0.7, 0.0, 10.0] {
            let actual = trajectory.advance(target).unwrap();
            let expected = [
                target.cos(),
                target.sin(),
                0.0,
                -target.sin(),
                target.cos(),
                0.0,
            ];
            for (actual, expected) in actual
                .position
                .into_iter()
                .chain(actual.velocity)
                .zip(expected)
            {
                assert!((actual - expected).abs() <= 64.0 * f64::EPSILON);
            }
        }
    }
}

#[test]
fn small_motion_is_not_swallowed_by_a_large_original_mean_anomaly() {
    let elements = Elements {
        mean_anomaly: 1e22,
        ..elements()
    };
    for method in [Method::WisdomHolman, Method::Mercurius, Method::Trace] {
        let mut trajectory = SolarTrajectory::new(elements, 1.0, method, Some(0.001)).unwrap();
        let initial = trajectory.state();
        let interval = 0.125_f64;
        let (sine, cosine) = interval.sin_cos();
        let rotate = |value: [f64; 3]| {
            [
                cosine * value[0] - sine * value[1],
                sine * value[0] + cosine * value[1],
                value[2],
            ]
        };
        let actual = trajectory.advance(interval).unwrap();
        assert_ne!(actual, initial);
        for (actual, expected) in actual.position.into_iter().chain(actual.velocity).zip(
            rotate(initial.position)
                .into_iter()
                .chain(rotate(initial.velocity)),
        ) {
            assert!((actual - expected).abs() <= 512.0 * f64::EPSILON);
        }
    }
}

#[test]
fn numerical_maps_retain_a_circular_orbit_without_a_method_fallback() {
    for method in [Method::Adaptive, Method::Leapfrog] {
        let mut trajectory = SolarTrajectory::new(elements(), 1.0, method, Some(0.001)).unwrap();
        let actual = trajectory.advance(1.0).unwrap();
        let radius = actual.position[0]
            .hypot(actual.position[1])
            .hypot(actual.position[2]);
        let momentum =
            actual.position[0] * actual.velocity[1] - actual.position[1] * actual.velocity[0];
        assert!((radius - 1.0).abs() < 1e-6);
        assert!((momentum - 1.0).abs() < 1e-6);
        assert!(
            actual
                .position
                .into_iter()
                .chain(actual.velocity)
                .all(f64::is_finite)
        );
    }
}

#[test]
fn nonfinite_epoch_and_step_do_not_advance_a_solar_trajectory() {
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            SolarTrajectory::new(elements(), 1.0, Method::Leapfrog, Some(bad)),
            Err(NumericalError::NonFinite)
        ));
        let mut trajectory = SolarTrajectory::new(elements(), 1.0, Method::Adaptive, None).unwrap();
        let initial = trajectory.state();
        assert_eq!(trajectory.advance(bad), Err(NumericalError::NonFinite));
        assert_eq!(trajectory.state(), initial);
    }
}
