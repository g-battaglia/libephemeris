// SPDX-License-Identifier: AGPL-3.0
//! Prepared independently defined orientation checks, not legacy-model claims.
use super::*;
fn norm(value: [f64; 3]) -> f64 {
    value[0].hypot(value[1]).hypot(value[2])
}

#[test]
fn iau_solar_pole_and_icrf_terrestrial_policy_have_their_explicit_definitions() {
    let first = PoleModel::IcrfTerrestrialIauSolar.at([J2000, 0.0]).unwrap();
    let other = PoleModel::IcrfTerrestrialIauSolar
        .at([J2000, 100000.0])
        .unwrap();
    assert_eq!(first.earth, [0.0, 0.0, 1.0]);
    assert_eq!(first.sun, other.sun);
    let ra = first.sun[1]
        .atan2(first.sun[0])
        .to_degrees()
        .rem_euclid(360.0);
    let dec = first.sun[2]
        .atan2(first.sun[0].hypot(first.sun[1]))
        .to_degrees();
    assert!((ra - 286.13).abs() <= 64.0 * f64::EPSILON * 286.13);
    assert!((dec - 63.87).abs() <= 64.0 * f64::EPSILON * 63.87);
    assert!((norm(first.sun) - 1.0).abs() <= 64.0 * f64::EPSILON);
}

#[test]
fn dated_pole_retains_unit_geometry_and_is_not_a_frozen_j2000_axis() {
    let mut previous = None;
    for offset in [-36525.0_f64, 0.0, 36525.0] {
        let poles = PoleModel::PublishedDateIers.at([J2000, offset]).unwrap();
        assert!((norm(poles.earth) - 1.0).abs() <= 256.0 * f64::EPSILON);
        assert!(poles.earth[0] != 0.0 || poles.earth[1] != 0.0);
        if let Some(value) = previous {
            assert_ne!(poles.earth, value);
        }
        previous = Some(poles.earth);
    }
}

#[test]
fn precessing_the_inertial_pole_back_to_date_recovers_the_published_datum() {
    for t in [-4.0_f64, 0.0, 6.0] {
        let date = J2000 + t * 36525.0;
        let poles = PoleModel::PublishedDateIers
            .at([J2000, t * 36525.0])
            .unwrap();
        let matrix = crate::positions::vondrak_precession_matrix(date, true).unwrap();
        let expected = terrestrial_datum_pole(t).unwrap();
        for (coefficients, reference) in matrix.iter().zip(expected) {
            let mut actual = Sum::default();
            for (coefficient, component) in coefficients.iter().zip(poles.earth) {
                actual.add(coefficient * component);
            }
            // The inherited IERS realization is first order, not an exact
            // orthogonal matrix; its quadratic bias remainder stays visible.
            assert!((actual.value() - reference).abs() <= 512.0 * f64::EPSILON);
        }
    }
}

#[test]
fn nonfinite_epoch_parts_fail_without_selecting_a_different_pole_policy() {
    for model in [
        PoleModel::IcrfTerrestrialIauSolar,
        PoleModel::PublishedDateIers,
    ] {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                model.at([J2000, bad]),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                model.at([bad, 0.0]),
                Err(Error::NonFiniteArgument { .. })
            ));
        }
    }
}
