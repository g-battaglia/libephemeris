//! Analytic tangent-plane derivatives on short owned observation intervals.

use super::*;
use crate::positions::derivative::spherical_angular_velocity;

#[test]
fn affine_transverse_motion_keeps_small_rates_on_short_owned_domains() {
    // Dyadic transverse coordinates have exact affine derivatives. Their
    // angular Jacobian at the x axis is velocity/radius, independently of
    // the derivative stencil and the absolute distance scale.
    for radius in [1.0 / 1024.0, 1.0, 1048576.0] {
        for sign in [-1.0, 1.0] {
            let sample = |offset: f64| {
                Ok([
                    sign * radius,
                    radius * offset / 64.0,
                    radius * offset / 65536.0,
                ])
            };
            for support in [(-1.0, 1.0), (0.0, 1.0 / 65536.0), (-1.0 / 65536.0, 0.0)] {
                let actual = cartesian_angular_velocity(&sample, support).unwrap();
                let exact = [
                    (sign / 64.0_f64).to_degrees(),
                    (1.0 / 65536.0_f64).to_degrees(),
                ];
                for (value, expected) in actual.into_iter().zip(exact) {
                    assert!((value - expected).abs() <= 16.0 * f64::EPSILON * expected.abs());
                }
            }
        }
    }
}

#[test]
fn stationary_pole_and_origin_keep_the_spherical_zero_rate_convention() {
    for ray in [[0.0; 3], [0.0, 0.0, 1.0], [0.0, 0.0, -1.0]] {
        assert_eq!(
            cartesian_angular_velocity(&|_| Ok(ray), (-1.0, 1.0)).unwrap(),
            [0.0; 2]
        );
    }
}

#[test]
fn physical_sample_errors_propagate_from_tangent_differentiation() {
    let sample = |offset| {
        if offset == 0.0 {
            Ok([1.0, 0.0, 0.0])
        } else {
            Err(Error::LebBodyNotFound { body_id: MERCURY })
        }
    };
    assert!(matches!(
        cartesian_angular_velocity(&sample, (-1.0, 1.0)),
        Err(Error::LebBodyNotFound { body_id: MERCURY })
    ));
}

#[test]
fn spherical_chain_rule_recovers_independent_rotation_and_range_motion() {
    let longitude_rate = 1.0 / 64.0;
    let latitude_rate = 1.0 / 128.0;
    let range_rate = 1.0 / 1024.0;
    let radial = [3.0 / 13.0, 4.0 / 13.0, 12.0 / 13.0];
    let east = [-4.0 / 5.0, 3.0 / 5.0, 0.0];
    let north = [-36.0 / 65.0, -48.0 / 65.0, 5.0 / 13.0];
    let velocity = std::array::from_fn(|i| {
        range_rate * radial[i] + 5.0 * longitude_rate * east[i] + 13.0 * latitude_rate * north[i]
    });
    let actual = spherical_angular_velocity([3.0, 4.0, 12.0], velocity).unwrap();
    for (value, expected) in actual.into_iter().zip([longitude_rate, latitude_rate]) {
        assert!((value - expected.to_degrees()).abs() <= 32.0 * f64::EPSILON);
    }
}

#[test]
fn correction_derivative_owns_the_actual_closest_approach_record() {
    let mut central = AngularModelSample::default();
    central.reads[0] = Some(CorrectionRead {
        body: SUN,
        epoch: 0.01,
        interval: Some((0.0, 1.0)),
    });
    central.reads[1] = Some(CorrectionRead {
        body: SUN,
        epoch: -0.1,
        interval: Some((-1.0, 0.0)),
    });
    // The observing record permits a long forward step, but the distinct
    // closest-approach record ends only 0.1 day after its actual query.
    assert_eq!(central.correction_support((-2.0, 2.0)), (-0.01, 0.1));
    let mut inside = central;
    inside.reads[0].as_mut().unwrap().epoch += 0.05;
    inside.reads[1].as_mut().unwrap().epoch += 0.05;
    assert!(central.owns_correction_sample(inside));
    let mut crossed = inside;
    crossed.reads[1] = Some(CorrectionRead {
        body: SUN,
        epoch: 0.05,
        interval: Some((0.0, 1.0)),
    });
    assert!(!central.owns_correction_sample(crossed));
    assert!(AngularModelSample::default().owns_correction_sample(AngularModelSample::default()));
    assert!(!AngularModelSample::default().owns_correction_sample(inside));
}
