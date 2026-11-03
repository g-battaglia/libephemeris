// SPDX-License-Identifier: AGPL-3.0
//! Prepared raw-trajectory ownership/endpoint checks using synthetic data.
use super::super::classic_de::{
    DeError,
    tests::{Fixture, Temporary},
};
use super::super::ephemeris::EphemerisError;
use super::*;

fn data() -> (Arc<Ephemeris>, Temporary) {
    let fixture = Fixture::force_ready();
    let (_reader, temporary) = fixture.open();
    let data = Arc::new(Ephemeris::open(temporary.path(), None).unwrap());
    (data, temporary)
}
fn poles(_epoch: [f64; 2]) -> Result<HarmonicPoles, Error> {
    Ok(HarmonicPoles {
        sun: [0.0, 0.0, 1.0],
        earth: [0.0, 0.0, 1.0],
    })
}
fn initial() -> State {
    State {
        position: [2.0, 3.0, 4.0],
        velocity: [0.01, 0.02, 0.03],
    }
}
fn norm(value: [f64; 3]) -> f64 {
    value[0].hypot(value[1]).hypot(value[2])
}

#[test]
fn zero_duration_is_exact_and_zero_step_cannot_serve_a_nonzero_interval() {
    let (data, _temporary) = data();
    let mut trajectory =
        ForcedTrajectory::new(data, initial(), [1008.0, 0.0], Some(0.0), None, None, poles)
            .unwrap();
    assert_eq!(trajectory.advance(0.0).unwrap(), initial());
    assert!(matches!(
        trajectory.advance(1.0),
        Err(IntegrationError::Numerical(NumericalError::InvalidStep))
    ));
    assert_eq!(trajectory.time(), 0.0);
    assert_eq!(trajectory.state(), initial());
}

#[test]
fn raw_endpoints_retain_units_and_allow_both_time_directions() {
    let (data, _temporary) = data();
    let mut trajectory = ForcedTrajectory::new(
        data,
        initial(),
        [1008.0, 0.0],
        Some(0.001),
        None,
        None,
        poles,
    )
    .unwrap();
    let endpoint = trajectory.advance(0.01).unwrap();
    assert_eq!(trajectory.time(), 0.01);
    assert_ne!(endpoint, initial());
    let returned = trajectory.advance(0.0).unwrap();
    assert_eq!(trajectory.time(), 0.0);
    for (actual, expected) in [
        (returned.position, initial().position),
        (returned.velocity, initial().velocity),
    ] {
        let scale = norm(expected);
        for axis in 0..3 {
            assert!((actual[axis] - expected[axis]).abs() <= 512.0 * f64::EPSILON * scale);
        }
    }
}

#[test]
fn request_failure_discards_even_previously_accepted_local_steps() {
    let (data, _temporary) = data();
    let mut trajectory = ForcedTrajectory::new(
        Arc::clone(&data),
        initial(),
        [1063.5, 0.0],
        Some(0.001),
        None,
        None,
        poles,
    )
    .unwrap();
    let result = trajectory.advance(10.0);
    assert!(matches!(
        result,
        Err(IntegrationError::Force(ForceError::Ephemeris(
            EphemerisError::De(DeError::OutsideCoverage)
        )))
    ));
    assert_eq!(trajectory.time(), 0.0);
    assert_eq!(trajectory.state(), initial());
    assert_eq!(trajectory.advance(0.0).unwrap(), initial());
    // Recovery must perform real nonzero propagation, not just return the
    // already committed zero-duration state. Reused scratch/cache contents
    // after the failed attempt must not alter the clean trajectory's result.
    let mut fresh = ForcedTrajectory::new(
        data,
        initial(),
        [1063.5, 0.0],
        Some(0.001),
        None,
        None,
        poles,
    )
    .unwrap();
    assert_eq!(
        trajectory.advance(0.125).unwrap(),
        fresh.advance(0.125).unwrap()
    );
    assert_eq!(trajectory.time(), 0.125);
}

#[test]
fn different_immutable_force_owners_remain_distinct_across_threads_and_failed_calls() {
    use std::sync::{Arc, Barrier};

    let first = Fixture::force_ready().with_solar_gm(1.0);
    let second = Fixture::force_ready().with_solar_gm(2.0);
    let (_first_reader, first_file) = first.open();
    let (_second_reader, second_file) = second.open();
    let owners = [
        Arc::new(Ephemeris::open(first_file.path(), None).unwrap()),
        Arc::new(Ephemeris::open(second_file.path(), None).unwrap()),
    ];
    assert_eq!(owners[0].catalogue.major_masses[0], 1.0);
    assert_eq!(owners[1].catalogue.major_masses[0], 2.0);

    let make = |index| {
        ForcedTrajectory::new(
            Arc::clone(&owners[index]),
            initial(),
            [1008.0, 0.0],
            Some(0.001),
            None,
            None,
            poles,
        )
        .unwrap()
    };
    let forces = [0, 1].map(|index| make(index).acceleration().unwrap());
    assert_ne!(forces[0], forces[1], "changed GM must alter actual force");
    let expected = [0, 1].map(|index| make(index).advance(0.01).unwrap());
    assert_ne!(
        expected[0], expected[1],
        "force change must reach the endpoint"
    );

    let requests = [1, 0, 1, 0, 0, 1, 0, 1];
    let rendezvous = Arc::new(Barrier::new(requests.len()));
    std::thread::scope(|scope| {
        let handles: Vec<_> = requests
            .into_iter()
            .map(|index| {
                let rendezvous = Arc::clone(&rendezvous);
                let mut trajectory = make(index);
                scope.spawn(move || {
                    rendezvous.wait();
                    let result = trajectory.advance(0.01).unwrap();
                    // The failed nonfinite call must neither publish a partial
                    // state nor change the owner's subsequent answer.
                    assert!(trajectory.advance(f64::NAN).is_err());
                    assert_eq!(trajectory.state(), result);
                    assert_eq!(trajectory.time(), 0.01);
                    result
                })
            })
            .collect();
        for (index, handle) in requests.into_iter().zip(handles) {
            assert_eq!(handle.join().unwrap(), expected[index]);
        }
    });
    for (index, owner) in owners.into_iter().enumerate() {
        let mut independent = ForcedTrajectory::new(
            owner,
            initial(),
            [1008.0, 0.0],
            Some(0.001),
            None,
            None,
            poles,
        )
        .unwrap();
        assert_eq!(independent.advance(0.01).unwrap(), expected[index]);
    }
}

#[cfg(feature = "nbody")]
#[test]
fn public_sessions_keep_explicit_synthetic_datasets_and_reject_wrong_owner_results() {
    use crate::nbody::{AssistEphemConfig, AssistEphemeris, AssistSession};

    let (_reader_a, file_a) = Fixture::force_ready().with_solar_gm(1.0).open();
    let (_reader_b, file_b) = Fixture::force_ready().with_solar_gm(2.0).open();
    let asteroids = Temporary::new(&crate::spk::force_test_image(
        crate::spk::DafEndian::Little,
        crate::spk::TYPE_CHEBYSHEV_POS,
        1,
        (1000.0 - 2451545.0) * 86400.0,
        32.0 * 86400.0,
    ));
    // Both paths are explicit. In particular, an omitted asteroid path would
    // permit HOME/ASSIST_DIR discovery and spoil this synthetic control.
    let configs = [file_a.path(), file_b.path()].map(|path| AssistEphemConfig {
        planets_file: Some(path.to_owned()),
        asteroids_file: Some(asteroids.path().to_owned()),
        data_dir: None,
    });
    let orbit = crate::minor_bodies::OrbitalElements::from_nbody_source_record(
        "synthetic owner".into(),
        1008.0,
        2.7,
        0.12,
        7.0,
        40.0,
        65.0,
        30.0,
        0.222,
    );
    let endpoint = 1008.01;
    let words = |value: crate::nbody::PropagationResult| {
        [
            value.x,
            value.y,
            value.z,
            value.vx,
            value.vy,
            value.vz,
            value.jd_tt,
        ]
        .map(f64::to_bits)
    };
    let expected = configs.each_ref().map(|config| {
        let direct = AssistEphemeris::open(Some(config)).unwrap();
        words(
            direct
                .propagate(&orbit, 1008.0, endpoint, None)
                .unwrap()
                .value,
        )
    });
    assert_ne!(expected[0], expected[1]);
    let matches_owner = |index: usize, value| expected[index] == value;
    // Positive discrimination: substituting the other owner's valid result
    // must fail, even though both requests have identical orbit/time fields.
    assert!(!matches_owner(0, expected[1]));
    assert!(!matches_owner(1, expected[0]));

    let sessions = configs.map(AssistSession::new);
    let requests = [1, 0, 0, 1, 1, 0, 1, 0];
    std::thread::scope(|scope| {
        let handles: Vec<_> = requests
            .into_iter()
            .map(|index| {
                let session = &sessions[index];
                let orbit = &orbit;
                scope.spawn(move || {
                    words(
                        session
                            .propagate(orbit, 1008.0, endpoint, None)
                            .unwrap()
                            .value,
                    )
                })
            })
            .collect();
        for (index, handle) in requests.into_iter().zip(handles) {
            assert!(matches_owner(index, handle.join().unwrap()));
        }
    });
    for (index, session) in sessions.iter().enumerate() {
        assert_eq!(
            session.config().planets_file.as_deref(),
            Some([file_a.path(), file_b.path()][index])
        );
        assert!(session.propagate(&orbit, 1008.0, f64::NAN, None).is_err());
        assert!(matches_owner(
            index,
            words(
                session
                    .propagate(&orbit, 1008.0, endpoint, None)
                    .unwrap()
                    .value
            )
        ));
    }
}

#[test]
fn concurrent_integrations_enter_force_callbacks_together_without_shared_state() {
    use std::sync::{Condvar, Mutex};
    use std::time::Duration;

    static ARRIVALS: (Mutex<usize>, Condvar) = (Mutex::new(0), Condvar::new());
    fn overlapping_poles(epoch: [f64; 2]) -> Result<HarmonicPoles, Error> {
        let (lock, condition) = &ARRIVALS;
        let mut count = lock.lock().unwrap();
        if *count < 2 {
            *count += 1;
            condition.notify_all();
            let (_count, timeout) = condition
                .wait_timeout_while(count, Duration::from_secs(30), |value| *value < 2)
                .unwrap();
            if timeout.timed_out() {
                return Err(Error::SingularMatrix);
            }
        }
        poles(epoch)
    }

    let (data, _temporary) = data();
    let construct = || {
        ForcedTrajectory::new(
            Arc::clone(&data),
            initial(),
            [1008.0, 0.0],
            Some(0.001),
            None,
            None,
            poles,
        )
        .unwrap()
    };
    let expected = construct().advance(0.01).unwrap();
    let mut workers = [construct(), construct()];
    for worker in &mut workers {
        // Arm only after construction, so the rendezvous proves two advance
        // calls overlap inside their actual force callbacks. The pole values
        // remain the same pure orientation. No runtime hook or timing ratio
        // is introduced; a serialized pair fails the bounded rendezvous.
        worker.poles = overlapping_poles;
    }
    assert_eq!(*ARRIVALS.0.lock().unwrap(), 0);
    std::thread::scope(|scope| {
        let handles: Vec<_> = workers
            .into_iter()
            .map(|mut worker| scope.spawn(move || worker.advance(0.01).unwrap()))
            .collect();
        for handle in handles {
            assert_eq!(handle.join().unwrap(), expected);
        }
    });
    assert_eq!(*ARRIVALS.0.lock().unwrap(), 2);
}

#[test]
fn nonfinite_inputs_and_orientation_failures_keep_distinct_categories() {
    let (data, _temporary) = data();
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            ForcedTrajectory::new(
                Arc::clone(&data),
                initial(),
                [1008.0, 0.0],
                Some(bad),
                None,
                None,
                poles
            ),
            Err(IntegrationError::Numerical(NumericalError::NonFinite))
        ));
        assert!(matches!(
            ForcedTrajectory::new(
                Arc::clone(&data),
                initial(),
                [1008.0, bad],
                None,
                None,
                None,
                poles
            ),
            Err(IntegrationError::Numerical(NumericalError::NonFinite))
        ));
    }
    fn unavailable(_epoch: [f64; 2]) -> Result<HarmonicPoles, Error> {
        Err(Error::SingularMatrix)
    }
    assert!(matches!(
        ForcedTrajectory::new(
            data,
            initial(),
            [1008.0, 0.0],
            None,
            None,
            None,
            unavailable
        ),
        Err(IntegrationError::Force(ForceError::Frame(
            Error::SingularMatrix
        )))
    ));
}
