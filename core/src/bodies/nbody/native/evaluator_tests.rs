// SPDX-License-Identifier: AGPL-3.0
//! Prepared retained-force wiring checks, not astronomical qualification.
use super::super::classic_de::{
    DeError,
    tests::{Fixture, Temporary},
};
use super::super::radau::{IntegrationError, Radau, Settings};
use super::*;

fn data() -> (Arc<Ephemeris>, Temporary) {
    let fixture = Fixture::force_ready();
    let (_reader, temporary) = fixture.open();
    let ephemeris = Arc::new(Ephemeris::open(temporary.path(), None).unwrap());
    (ephemeris, temporary)
}
fn poles() -> HarmonicPoles {
    HarmonicPoles {
        sun: [0.0, 0.0, 1.0],
        earth: [0.0, 0.0, 1.0],
    }
}
fn state() -> State {
    State {
        position: [2.0, 3.0, 4.0],
        velocity: [0.01, 0.02, 0.03],
    }
}

#[test]
fn retained_evaluation_composes_the_full_selected_frame_and_every_force_family() {
    let (data, _temporary) = data();
    let parameters = Some([1e-8, -2e-8, 3e-8]);
    let mut evaluator = Evaluator::new(Arc::clone(&data), poles(), parameters, None).unwrap();
    let mut frame = Vec::new();
    data.frame([1008.0, 0.0], &mut frame).unwrap();
    let model = data.catalogue.force_model(poles(), parameters).unwrap();
    let expected = model.acceleration(state(), &frame, None).unwrap();
    assert_eq!(
        evaluator
            .acceleration([1008.0, 0.0], state(), poles())
            .unwrap(),
        expected
    );
    assert_eq!(evaluator.sources.len(), data.particle_count());
    assert!(evaluator.model.relativity.is_some());
    assert_eq!(evaluator.model.zonals.len(), 2);
    assert!(evaluator.model.nongravitational.is_some());
}

#[test]
fn new_epoch_poles_change_only_orientation_not_coefficients_or_source_selection() {
    let (data, _temporary) = data();
    let mut evaluator = Evaluator::new(Arc::clone(&data), poles(), None, None).unwrap();
    let original = evaluator
        .acceleration([1008.0, 0.0], state(), poles())
        .unwrap();
    let changed = HarmonicPoles {
        sun: [1.0, 0.0, 0.0],
        earth: [0.0, 1.0, 0.0],
    };
    let actual = evaluator
        .acceleration([1008.0, 0.0], state(), changed)
        .unwrap();
    let mut frame = Vec::new();
    data.frame([1008.0, 0.0], &mut frame).unwrap();
    let expected = data
        .catalogue
        .force_model(changed, None)
        .unwrap()
        .acceleration(state(), &frame, None)
        .unwrap();
    assert_eq!(actual, expected);
    assert_ne!(actual, original);
    assert_eq!(
        evaluator
            .acceleration([1008.0, 0.0], state(), poles())
            .unwrap(),
        original
    );
}

#[test]
fn data_and_domain_failures_keep_their_original_typed_categories() {
    let (data, _temporary) = data();
    let mut evaluator = Evaluator::new(data, poles(), None, None).unwrap();
    assert!(matches!(
        evaluator.acceleration([999.0, 0.0], state(), poles()),
        Err(ForceError::Ephemeris(EphemerisError::De(
            DeError::OutsideCoverage
        )))
    ));
    let invalid = State {
        position: [f64::NAN, 0.0, 0.0],
        ..state()
    };
    assert!(matches!(
        evaluator.acceleration([1008.0, 0.0], invalid, poles()),
        Err(ForceError::Numerical(NumericalError::NonFinite))
    ));
    let invalid_pole = HarmonicPoles {
        earth: [0.0; 3],
        ..poles()
    };
    assert!(matches!(
        evaluator.acceleration([1008.0, 0.0], state(), invalid_pole),
        Err(ForceError::Numerical(NumericalError::InvalidSettings))
    ));
}

#[test]
fn ephemeris_error_reaches_radau_once_without_force_replacement_or_a_fake_endpoint() {
    let (data, _temporary) = data();
    let mut evaluator = Evaluator::new(data, poles(), None, None).unwrap();
    let mut integrator = Radau::new(
        state(),
        0.0,
        Settings {
            coefficient_tolerance: 1e-9,
            initial_step: 0.1,
            maximum_step: 1.0,
            maximum_attempts: 100,
            maximum_iterations: 32,
        },
    )
    .unwrap();
    let mut calls = 0;
    let result = integrator.advance_with_error(1.0, &mut |time, state| {
        calls += 1;
        evaluator.acceleration([999.0, time], *state, poles())
    });
    assert!(matches!(
        result,
        Err(IntegrationError::Force(ForceError::Ephemeris(
            EphemerisError::De(DeError::OutsideCoverage)
        )))
    ));
    assert_eq!(calls, 1);
    assert_eq!(integrator.time(), 0.0);
    assert_eq!(integrator.state(), state());
}

#[test]
fn explicit_self_identifier_does_not_disable_other_sources_or_force_families() {
    let (data, _temporary) = data();
    let primary = data.particle(0, [1008.0, 0.0]).unwrap();
    let state = State {
        position: primary.point.position,
        velocity: primary.velocity,
    };
    let mut singular = Evaluator::new(Arc::clone(&data), poles(), None, None).unwrap();
    assert!(matches!(
        singular.acceleration([1008.0, 0.0], state, poles()),
        Err(ForceError::Numerical(NumericalError::SingularForce {
            body: 10
        }))
    ));
    let mut selected = Evaluator::new(data, poles(), None, Some(10)).unwrap();
    let actual = selected
        .acceleration([1008.0, 0.0], state, poles())
        .unwrap();
    assert!(actual.into_iter().any(|value| value != 0.0));
    assert_eq!(selected.sources.len(), 11);
    assert!(selected.model.relativity.is_some());
    assert_eq!(selected.model.zonals.len(), 2);
}
