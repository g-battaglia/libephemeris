// SPDX-License-Identifier: AGPL-3.0
//! Prepared safe-adapter shape/ownership checks, not final API qualification.
use super::super::classic_de::tests::Fixture;
use super::*;

fn orbit() -> Orbit {
    Orbit {
        semi_major_axis: 1.0,
        eccentricity: 0.0,
        inclination: 0.0,
        perihelion_argument: 0.0,
        ascending_node: 0.0,
        mean_anomaly: 0.0,
    }
}

#[test]
fn all_public_selector_labels_construct_the_native_solar_problem_honestly() {
    for method in [
        ReboundIntegrator::Ias15,
        ReboundIntegrator::Whfast,
        ReboundIntegrator::Leapfrog,
        ReboundIntegrator::Mercurius,
        ReboundIntegrator::Trace,
    ] {
        let mut simulation = Simulation::two_body(orbit(), method, Some(0.001)).unwrap();
        assert_eq!(simulation.state(0).unwrap().position, [0.0; 3]);
        let initial = simulation.state(1).unwrap();
        assert_eq!(initial.position, [1.0, 0.0, 0.0]);
        assert!(initial.velocity[1] > 0.0);
        assert!(simulation.integrate(0.0).unwrap().is_empty());
        assert_eq!(simulation.state(1).unwrap(), initial);
        assert!(simulation.state(2).is_err());
    }
    assert!(versions().0.starts_with("libephemeris-native-"));
    assert!(versions().1.starts_with("libephemeris-native-"));
    assert!(!versions().0.starts_with("4.6.0"));
    assert!(!versions().1.starts_with("1.2.3"));
}

#[test]
fn prescribed_staging_retains_one_raw_target_and_reports_real_dataset_values() {
    let fixture = Fixture::force_ready();
    let (_reader, temporary) = fixture.open();
    let ephemeris = Ephemeris::open(temporary.path(), None).unwrap();
    assert_eq!(ephemeris.coverage(), [(1000.0, 1064.0)]);
    assert_eq!(
        ephemeris.curated_initial_state(2_000_001, 1008.0).unwrap(),
        None
    );
    assert!(
        ephemeris
            .curated_initial_state(2_000_001, f64::NAN)
            .is_err()
    );
    let time = 1008.0 - ephemeris.reference_jd();
    let source = ephemeris.state(0, time).unwrap();
    assert_eq!(source.position, [48.5, 49.0, 49.5]);
    assert_eq!(ephemeris.position(0, time).unwrap(), source.position);
    assert_eq!(ephemeris.gravitational_parameter(0, time).unwrap(), 1.0);
    assert!(
        ephemeris
            .sun_acceleration(time)
            .unwrap()
            .into_iter()
            .all(f64::is_finite)
    );
    let mut simulation = Simulation::assist(&ephemeris).unwrap();
    let initial = CartesianState {
        position: [2.0, 3.0, 4.0],
        velocity: [0.01, 0.02, 0.03],
    };
    simulation.add(initial);
    simulation.set_time(time);
    simulation.non_gravitational_parameters([0.0; 3]);
    simulation.exclude_asteroid(11); // no asteroid file: no physical source removed.
    assert_eq!(simulation.state(0).unwrap(), initial);
    assert_eq!(simulation.asteroid_acceleration().unwrap(), [0.0; 3]);
    assert!(simulation.integrate(time).unwrap().is_empty());
    assert_eq!(simulation.state(0).unwrap(), initial);
    assert!(simulation.state(1).is_err());
    assert!(
        simulation
            .acceleration()
            .unwrap()
            .into_iter()
            .all(f64::is_finite)
    );
    assert_eq!(simulation.state(0).unwrap(), initial);
}

#[test]
fn invalid_native_requests_do_not_select_legacy_computation() {
    let mut simulation =
        Simulation::two_body(orbit(), ReboundIntegrator::Ias15, Some(0.0)).unwrap();
    assert!(simulation.integrate(1.0).is_err());
    assert!(simulation.integrate(f64::NAN).is_err());
    let invalid = Orbit {
        eccentricity: 1.0,
        ..orbit()
    };
    assert!(matches!(
        Simulation::two_body(invalid, ReboundIntegrator::Whfast, Some(0.1)),
        Err(Error::InvalidOrbitElements { .. })
    ));
    assert!(matches!(
        Ephemeris::open(Path::new("invalid\0force-path"), None),
        Err(Error::InvalidOrbitElements { .. })
    ));
}

fn assert_terminal_propagation_error(error: Error, expected_diagnostic: &str) {
    match &error {
        Error::NBodyNative {
            stage,
            code,
            diagnostic,
        } => {
            assert_eq!((*stage, *code), (4, 1));
            assert_eq!(diagnostic, expected_diagnostic);
        }
        other => panic!("propagation failure lost its native phase: {other:?}"),
    }
    assert_eq!(error.category(), crate::errors::ErrorCategory::Calculation);
    assert!(
        !crate::nbody::assist_fallback(&error),
        "a consumed solver/force failure must not select a weaker model"
    );
}

#[test]
fn every_numerical_propagation_domain_stays_terminal_with_its_real_diagnostic() {
    for numerical in [
        NumericalError::NonFinite,
        NumericalError::InvalidStep,
        NumericalError::InvalidSettings,
        NumericalError::InvalidOrbit("undefined propagated orbital plane"),
        NumericalError::SingularForce { body: 2_000_001 },
        NumericalError::InvalidMass { body: 2_000_001 },
        NumericalError::DuplicateSource { body: 2_000_001 },
        NumericalError::MissingSource { body: 10 },
        NumericalError::IterationLimit,
        NumericalError::StepLimit,
        NumericalError::StepUnderflow,
    ] {
        let expected = format!("native numerical failure: {numerical:?}");
        assert_terminal_propagation_error(
            propagation_error(IntegrationError::Numerical(numerical)),
            &expected,
        );
    }

    // Exercise the actual Solar propagation boundary as well as the direct
    // typed mapping. A valid constructor cannot make its later zero-step
    // failure eligible for another propagation model.
    let mut simulation =
        Simulation::two_body(orbit(), ReboundIntegrator::Ias15, Some(0.0)).unwrap();
    assert_terminal_propagation_error(
        simulation.integrate(1.0).unwrap_err(),
        "native numerical failure: InvalidStep",
    );
}

#[test]
fn selected_force_errors_preserve_their_diagnostics_and_cannot_become_input_errors() {
    use super::super::classic_de::DeError;
    use super::super::ephemeris::EphemerisError;
    use super::super::evaluator::ForceError;

    for force in [
        ForceError::Numerical(NumericalError::InvalidMass { body: 2_000_001 }),
        ForceError::Frame(Error::InvalidOrbitElements {
            entry: "selected-force",
            detail: "invalid consumed orientation".into(),
        }),
        ForceError::Ephemeris(EphemerisError::De(DeError::Malformed(
            "selected coefficient record",
        ))),
        ForceError::Ephemeris(EphemerisError::De(DeError::Io(std::io::Error::new(
            std::io::ErrorKind::UnexpectedEof,
            "selected record read failed",
        )))),
        ForceError::Ephemeris(EphemerisError::De(DeError::Io(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            "consumed source became unavailable",
        )))),
        ForceError::Ephemeris(EphemerisError::Spk(Error::SpkFileNotFound {
            path: "consumed-source.bsp".into(),
        })),
    ] {
        let expected = force.to_string();
        assert_terminal_propagation_error(
            propagation_error(IntegrationError::Force(force)),
            &expected,
        );
    }
}

#[test]
fn invalid_initial_orbit_domains_retain_their_input_classification() {
    for (numerical, expected_detail) in [
        (
            NumericalError::InvalidOrbit("inconsistent initial elements"),
            "inconsistent initial elements",
        ),
        (
            NumericalError::InvalidMass { body: 10 },
            "invalid GM for body 10",
        ),
    ] {
        let error = initial_orbit_error(numerical);
        assert!(matches!(
            &error,
            Error::InvalidOrbitElements { entry: "nbody", detail } if detail == expected_detail
        ));
        assert_eq!(
            error.category(),
            crate::errors::ErrorCategory::InputValidation
        );
        assert!(crate::nbody::assist_fallback(&error));
    }
}

#[test]
fn acquisition_admits_only_structured_missing_files_without_hiding_initialization_failures() {
    let planets = Path::new("selected-planets.440");
    for (source, expected_path) in [
        (
            EphemerisError::De(DeError::Io(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                "missing selected planet file",
            ))),
            "selected-planets.440",
        ),
        (
            EphemerisError::Spk(Error::SpkFileNotFound {
                path: "selected-asteroids.bsp".into(),
            }),
            "selected-asteroids.bsp",
        ),
    ] {
        let expected_message = source.to_string();
        let error = acquisition_error(source, planets);
        assert!(matches!(
            &error,
            Error::EphemerisLoaderIo { path, kind: std::io::ErrorKind::NotFound, errno: None, message }
                if path == expected_path && message == &expected_message
        ));
        assert!(crate::nbody::assist_fallback(&error));
    }

    for source in [
        EphemerisError::De(DeError::Malformed("selected DE header")),
        EphemerisError::De(DeError::Io(std::io::Error::new(
            std::io::ErrorKind::PermissionDenied,
            "selected file access denied",
        ))),
        EphemerisError::De(DeError::Io(std::io::Error::new(
            std::io::ErrorKind::UnexpectedEof,
            "selected header is truncated",
        ))),
    ] {
        let expected_diagnostic = source.to_string();
        let error = acquisition_error(source, planets);
        assert!(matches!(
            &error,
            Error::NBodyNative { stage: 2, code: 1, diagnostic } if diagnostic == &expected_diagnostic
        ));
        assert!(!crate::nbody::assist_fallback(&error));
    }
}

#[test]
fn nonfinite_computed_output_stays_terminal_without_reclassifying_initial_arguments() {
    for state in [
        State {
            position: [f64::NAN, 0.0, 0.0],
            velocity: [0.0; 3],
        },
        State {
            position: [0.0; 3],
            velocity: [f64::INFINITY, 0.0, 0.0],
        },
    ] {
        let initial = CartesianState {
            position: state.position,
            velocity: state.velocity,
        };
        let argument_error = initial.check().unwrap_err();
        let expected_diagnostic = argument_error.to_string();
        assert!(matches!(argument_error, Error::NonFiniteArgument { .. }));
        assert!(crate::nbody::assist_fallback(&argument_error));

        let output = external(state).unwrap_err();
        assert!(matches!(
            &output,
            Error::NBodyNative { stage: 5, code: 1, diagnostic }
                if diagnostic == &expected_diagnostic
        ));
        assert_eq!(output.category(), crate::errors::ErrorCategory::Calculation);
        assert!(!crate::nbody::assist_fallback(&output));
    }
}

#[cfg(feature = "nbody")]
#[test]
fn coincident_asteroid_identity_survives_reordering_and_more_than_sixteen_sources() {
    use super::super::super::{
        AssistEphemConfig, AssistEphemeris, assist_barycentric_initial, rotated, shifted,
    };
    use super::super::classic_de::tests::Temporary;
    use crate::errors::Warning;

    let mut fixture = Fixture::force_ready();
    for number in (1..=20).filter(|number| ![1, 2, 4].contains(number)) {
        fixture = fixture.with_constant(&format!("MA{number:04}"), number as f64 / 1024.0);
    }
    let (_reader, planets) = fixture.open();
    // Put Ceres first, then after the former sixteen-source search boundary.
    for selected in [0, 19] {
        let mut targets: Vec<i32> = (2_000_002..=2_000_020).collect();
        targets.insert(selected, 2_000_001);
        let mut bytes = crate::spk::force_test_image(
            crate::spk::DafEndian::Little,
            crate::spk::TYPE_CHEBYSHEV_POS,
            targets.len(),
            (1000.0 - 2451545.0) * 86400.0,
            32.0 * 86400.0,
        );
        for (index, target) in targets.into_iter().enumerate() {
            let offset = 1024 + (3 + index * 5) * 8 + 16;
            bytes[offset..offset + 4].copy_from_slice(&target.to_le_bytes());
        }
        let asteroids = Temporary::new(&bytes);
        let ephemeris = Ephemeris::open(planets.path(), Some(asteroids.path())).unwrap();
        let time = 1008.0 - ephemeris.reference_jd();
        let initial = ephemeris.state(selected as i32 + 11, time).unwrap();
        let (mut simulation, warnings) =
            assist_barycentric_initial(initial, time, &ephemeris, None, Some(2_000_001)).unwrap();
        assert!(
            matches!(
                warnings.as_slice(),
                [Warning::ReboundSelfInteraction { perturber }] if perturber == "Ceres"
            ),
            "selected ordinal={selected}; warnings={warnings:?}"
        );
        let acceleration = simulation.asteroid_acceleration().unwrap();
        assert!(acceleration.into_iter().all(f64::is_finite));
        assert!(acceleration.into_iter().any(|value| value != 0.0));
        let mut unexcluded = Simulation::assist(&ephemeris).unwrap();
        unexcluded.add(initial);
        unexcluded.set_time(time);
        assert!(
            unexcluded.asteroid_acceleration().is_err(),
            "self-force control must be singular"
        );

        let config = AssistEphemConfig {
            planets_file: Some(planets.path().to_owned()),
            asteroids_file: Some(asteroids.path().to_owned()),
            data_dir: None,
        };
        let public = AssistEphemeris::open(Some(&config)).unwrap();
        // This fixture deliberately declares a nonphysical AU of 2 km.
        // Form the public ecliptic seed in its selected force units rather
        // than using a physical-SPK convenience conversion with the real AU.
        let seed = rotated(
            shifted(initial, public.inner.state(0, time).unwrap(), false),
            false,
        );
        let outcome = public
            .propagate_state_for_target(seed, 1008.0, 1008.001, None, Some(2_000_001))
            .unwrap();
        assert!(matches!(
            outcome.warnings.as_slice(),
            [Warning::ReboundSelfInteraction { perturber }] if perturber == "Ceres"
        ));
        assert_eq!(outcome.value.jd_tt, 1008.001);

        // A distinct particle inside the former 1e-4 AU exclusion radius
        // still feels Ceres. Anonymous public propagation carries no self ID;
        // inspect its actual staged force without advancing a stiff encounter.
        let nearby = CartesianState {
            position: [
                initial.position[0] + 1e-5,
                initial.position[1],
                initial.position[2],
            ],
            velocity: [
                initial.velocity[0] + 0.125,
                initial.velocity[1],
                initial.velocity[2],
            ],
        };
        let anonymous_seed = rotated(
            shifted(nearby, public.inner.state(0, time).unwrap(), false),
            false,
        );
        let anonymous = public
            .propagate_state(anonymous_seed, 1008.0, 1008.0, None)
            .unwrap();
        assert!(anonymous.warnings.is_empty());
        let observed = public
            .propagated_state_with(
                anonymous_seed,
                1008.0,
                1008.0,
                None,
                None,
                |simulation, ephemeris, time| {
                    let target = simulation.state(0)?;
                    let actual = simulation.asteroid_acceleration()?;
                    let (mut excluded, _) =
                        assist_barycentric_initial(target, time, ephemeris, None, Some(2_000_001))?;
                    let without_self = excluded.asteroid_acceleration()?;
                    let source = ephemeris.state(selected as i32 + 11, time)?;
                    let offset: [f64; 3] =
                        std::array::from_fn(|axis| source.position[axis] - target.position[axis]);
                    let radius = offset[0].hypot(offset[1]).hypot(offset[2]);
                    assert!(radius > 0.0 && radius < 1e-4);
                    let gm = ephemeris.gravitational_parameter(selected as i32 + 11, time)?;
                    for axis in 0..3 {
                        let expected = (gm / radius / radius) * (offset[axis] / radius);
                        let contribution = actual[axis] - without_self[axis];
                        let rounding = 16.0 * f64::EPSILON * actual[axis].abs().max(expected.abs());
                        assert!((contribution - expected).abs() <= rounding);
                    }
                    Ok(())
                },
            )
            .unwrap();
        assert!(observed.warnings.is_empty());
    }
}
