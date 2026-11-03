//! Provider-free acquisition, branch, residual and error-ownership controls.
use super::*;
use crate::constants::{FLG_SWIEPH, MARS, MOON};

fn position(longitude: f64, speed: f64) -> [f64; 6] {
    [longitude, 0.0, 1.0, speed, 0.0, 0.0]
}

fn range(jd: f64) -> Error {
    Error::EphemerisRange {
        requested_jd: jd,
        start_jd: 900.0,
        end_jd: 1100.0,
        body: Some(MARS),
    }
}

fn io_failure() -> Error {
    Error::EphemerisLoaderIo {
        path: "/owned/fake-root".into(),
        kind: std::io::ErrorKind::AlreadyExists,
        errno: Some(17),
        message: "fake filesystem refusal".into(),
    }
}

#[test]
fn raw_scan_reacquires_start_then_sample_before_exact_root_return() {
    let mut calls = Vec::new();
    let mut get = |jd, flags| {
        calls.push((jd, flags));
        Ok(position(0.0, -1.0))
    };
    assert_eq!(
        planet(&mut get, MARS, 0.0, 1000.0, FLG_SWIEPH).unwrap(),
        1000.0
    );
    assert_eq!(
        calls,
        vec![
            (1000.0, FLG_SWIEPH | FLG_SPEED),
            (1000.0, FLG_SWIEPH),
            (1000.5, FLG_SWIEPH)
        ]
    );
}

#[test]
fn scan_bisection_keeps_raw_flags_and_has_no_rate_guard_tail() {
    let mut calls = Vec::new();
    let mut get = |jd, flags| {
        calls.push((jd, flags));
        Ok(position(jd - 1000.25, -1.0))
    };
    assert_eq!(
        planet(&mut get, MARS, 0.0, 1000.0, FLG_SWIEPH).unwrap(),
        1000.25
    );
    assert_eq!(calls.len(), 64);
    assert_eq!(calls[0], (1000.0, FLG_SWIEPH | FLG_SPEED));
    assert_eq!(calls[3], (1000.0, FLG_SWIEPH));
    assert!(calls[1..].iter().all(|call| call.1 == FLG_SWIEPH));
}

#[test]
fn initial_range_failure_owns_stage_and_original_cause() {
    let original = range(1000.0);
    let mut get = |_, _| Err(original.clone());
    let error = luminary(&mut get, SUN, 10.0, 1000.0, 0, SearchDirection::Forward).unwrap_err();
    match &error {
        Error::PublicCrossing {
            message,
            cause: Some(cause),
            ..
        } => {
            assert_eq!(
                message,
                &format!("Failed to calculate Sun position: {original}")
            );
            assert_eq!(cause.as_ref(), &original);
        }
        other => panic!("wrong public stage {other:?}"),
    }
    assert_eq!(
        std::error::Error::source(&error).unwrap().to_string(),
        original.to_string()
    );
}

#[test]
fn newton_value_failure_is_wrapped_without_starting_range_recovery() {
    let mut calls = 0;
    let original = Error::MathDomain { value: -2.0 };
    let mut get = |_, _| {
        calls += 1;
        if calls == 1 {
            Ok(position(0.0, 0.5))
        } else {
            Err(original.clone())
        }
    };
    let error = planet(&mut get, MARS, 90.0, 1000.0, 0).unwrap_err();
    assert_eq!(calls, 2);
    match error {
        Error::PublicCrossing {
            message,
            cause: Some(cause),
            ..
        } => {
            assert_eq!(
                message,
                format!("Failed to calculate planet position during iteration: {original}")
            );
            assert_eq!(*cause, original);
        }
        other => panic!("wrong failure {other:?}"),
    }
}

#[test]
fn range_recovery_retains_the_fresh_second_range_error() {
    let mut calls = Vec::new();
    let mut get = |jd, flags| {
        calls.push((jd, flags));
        match calls.len() {
            1 => Ok(position(0.0, 0.5)),
            2 => Err(range(1180.0)),
            _ => Err(range(1000.0)),
        }
    };
    let error = planet(&mut get, MARS, 90.0, 1000.0, FLG_SWIEPH).unwrap_err();
    assert_eq!(error.primary_error(), &range(1000.0));
    assert!(std::error::Error::source(&error).is_none());

    assert_eq!(
        calls,
        vec![
            (1000.0, FLG_SWIEPH | FLG_SPEED),
            (1180.0, FLG_SWIEPH | FLG_SPEED),
            (1000.0, FLG_SWIEPH | FLG_SPEED)
        ]
    );
}

#[test]
fn recovery_value_error_propagates_instead_of_restoring_old_range() {
    let original = Error::MathDomain { value: -3.0 };
    let mut calls = 0;
    let mut get = |_, _| {
        calls += 1;
        match calls {
            1 => Ok(position(0.0, 0.5)),
            2 => Err(range(1180.0)),
            _ => Err(original.clone()),
        }
    };
    assert_eq!(
        planet(&mut get, MARS, 90.0, 1000.0, 0)
            .unwrap_err()
            .primary_error(),
        &original
    );
    assert_eq!(calls, 3);
}

#[test]
fn recovery_public_lookup_failure_restores_original_range() {
    let mut calls = 0;
    let mut get = |_, _| {
        calls += 1;
        match calls {
            1 => Ok(position(0.0, 0.5)),
            2 => Err(range(1180.0)),
            _ => Err(Error::UnknownBody { body_id: 99999 }),
        }
    };
    assert_eq!(
        planet(&mut get, MARS, 90.0, 1000.0, 0)
            .unwrap_err()
            .primary_error(),
        &range(1180.0)
    );
    assert_eq!(calls, 3);
}

#[test]
fn raw_scan_configuration_error_is_unwrapped() {
    let original = Error::ConfigurationMissing {
        missing_config: "observer_location",
        suggestion: "Call set_topo(lon, lat, alt) first",
    };
    let mut calls = 0;
    let mut get = |_, _| {
        calls += 1;
        if calls == 1 {
            Ok(position(0.0, -1.0))
        } else {
            Err(original.clone())
        }
    };
    assert_eq!(
        planet(&mut get, MARS, 10.0, 1000.0, 0).unwrap_err(),
        original
    );
    assert_eq!(calls, 2);
}

#[test]
fn entry_and_iteration_io_errors_remain_unwrapped() {
    for initial_failure in [true, false] {
        let original = io_failure();
        let mut calls = 0;
        let mut get = |_, _| {
            calls += 1;
            if initial_failure || calls > 1 {
                Err(original.clone())
            } else {
                Ok(position(0.0, 1.0))
            }
        };
        assert_eq!(
            luminary(&mut get, SUN, 10.0, 1000.0, 0, SearchDirection::Forward).unwrap_err(),
            original
        );
        assert_eq!(calls, if initial_failure { 1 } else { 2 });
    }
}

#[test]
fn newton_antipode_retains_positive_180_correction() {
    let mut times = Vec::new();
    let mut get = |jd, _| {
        times.push(jd);
        Ok(position(
            match times.len() {
                1 => 0.0,
                2 => 190.0,
                _ => 10.0,
            },
            1.0,
        ))
    };
    assert_eq!(
        luminary(&mut get, SUN, 10.0, 1000.0, 0, SearchDirection::Forward).unwrap(),
        1190.0
    );
    assert_eq!(times, vec![1000.0, 1010.0, 1190.0]);
    assert_eq!(signed(-180.0), 180.0);
    assert_eq!(modulo360(-0.0).to_bits(), 0.0_f64.to_bits());
    assert_eq!(normalize_target(-0.0).to_bits(), 0.0_f64.to_bits());
    assert_eq!(scan_delta(180.0, 0.0), -180.0);
}

#[test]
fn luminary_exhaustion_has_only_initial_plus_fifty_samples() {
    for body in [SUN, MOON] {
        let mut calls = 0;
        let mut get = |_, flags| {
            assert_eq!(flags, FLG_SWIEPH | FLG_SPEED);
            calls += 1;
            // A large speed prevents the one-degree residual leaving the envelope.
            Ok(position(if calls == 1 { 0.0 } else { 9.0 }, 1000.0))
        };
        let error = luminary(
            &mut get,
            body,
            10.0,
            1000.0,
            FLG_SWIEPH,
            SearchDirection::Forward,
        )
        .unwrap_err();
        assert_eq!(calls, 51);
        assert!(matches!(error, Error::PublicCrossing { cause: None, .. }));
    }
}

#[test]
fn generic_cap_stays_sixty_when_later_speed_requests_eighty() {
    let mut calls = 0;
    let stop = io_failure();
    let mut get = |jd, _| {
        calls += 1;
        if calls == 1 {
            return Ok(position(0.0, 0.5));
        }
        // The 61st callback is the last Newton sample. The next callback
        // reacquires start for the exhaustion forward bracket.
        if calls == 62 {
            assert_eq!(jd, 1000.0);
            return Err(stop.clone());
        }
        Ok(position(90.0 - 1e-6, 0.015))
    };
    assert_eq!(planet(&mut get, MARS, 90.0, 1000.0, 0).unwrap_err(), stop);
    assert_eq!(calls, 62);
}

#[test]
fn heliocentric_cap_has_no_exhaustion_read() {
    let mut calls = 0;
    let mut get = |_, flags| {
        calls += 1;
        assert_eq!(flags, FLG_SWIEPH | FLG_HELCTR | FLG_SPEED);
        Ok(position(
            if calls == 1 { 0.0 } else { 10.0 - 1e-6 },
            if calls == 1 { 1.0 } else { 0.015 },
        ))
    };
    let error = helio(
        &mut get,
        MARS,
        10.0,
        1000.0,
        FLG_SWIEPH,
        SearchDirection::Forward,
    )
    .unwrap_err();
    assert_eq!(calls, 61);
    assert_eq!(
        error.to_string(),
        "Maximum iterations reached in heliocentric crossing calculation"
    );
}

#[test]
fn backward_luminary_seed_uses_previous_cycle_at_start_root() {
    for (body, rate) in [(SUN, 0.9856), (MOON, 13.176)] {
        let mut times = Vec::new();
        let mut get = |jd, _| {
            times.push(jd);
            Ok(position(10.0, rate))
        };
        let expected = 1000.0 + -360.0 / rate;
        assert_eq!(
            luminary(&mut get, body, 10.0, 1000.0, 0, SearchDirection::Backward)
                .unwrap()
                .to_bits(),
            expected.to_bits()
        );
        assert_eq!(times, vec![1000.0, expected]);
    }
}

#[test]
fn node_imminent_root_consumes_component_four_and_retains_coordinates() {
    for direction in [SearchDirection::Forward, SearchDirection::Backward] {
        let dir = if direction == SearchDirection::Forward {
            1.0
        } else {
            -1.0
        };
        let root = 1000.0 + dir * 5e-7;
        let mut calls = Vec::new();
        let mut get = |jd, flags| {
            calls.push((jd, flags));
            // Component three deliberately disagrees with the latitude rate.
            Ok([27.0, jd - root, 1.0, 700.0, 1.0, 0.0])
        };
        let hit = node(&mut get, 1000.0, FLG_SWIEPH, direction).unwrap();
        assert_eq!(hit.time.to_bits(), root.to_bits());
        assert_eq!(hit.longitude.to_bits(), 27.0_f64.to_bits());
        assert_eq!(hit.latitude.to_bits(), 0.0_f64.to_bits());
        assert_eq!(
            calls,
            vec![
                (1000.0, FLG_SWIEPH | FLG_SPEED),
                (root, FLG_SWIEPH | FLG_SPEED)
            ]
        );
    }
}

#[test]
fn node_scan_error_has_entry_stage_and_no_later_sample() {
    let mut calls = 0;
    let original = range(1002.0);
    let mut get = |_, _| {
        calls += 1;
        if calls == 1 {
            Ok([1.0, 1.0, 1.0, 500.0, 1.0, 0.0])
        } else {
            Err(original.clone())
        }
    };
    match node(&mut get, 1000.0, 0, SearchDirection::Forward)
        .err()
        .unwrap()
    {
        Error::PublicCrossing {
            message,
            cause: Some(cause),
            ..
        } => {
            assert_eq!(
                message,
                format!("Failed to calculate Moon position: {original}")
            );
            assert_eq!(*cause, original);
        }
        other => panic!("wrong node scan stage {other:?}"),
    }
    assert_eq!(calls, 2);
}

#[test]
fn crossing_recovery_preserves_its_native_cause() {
    let cause = Error::ConfigurationMissing {
        missing_config: "observer_location",
        suggestion: "supply an observer",
    };
    let error = failure_with_cause("crossing refinement failed", cause.clone());
    assert_eq!(error.cause(), Some(&cause));
    assert_eq!(
        std::error::Error::source(&error).unwrap().to_string(),
        cause.to_string()
    );
}

#[test]
fn analytic_linear_longitude_crossings_respect_time_direction() {
    for direction in [SearchDirection::Forward, SearchDirection::Backward] {
        let mut sample = |jd: f64, _: i32| Ok(position((jd - 1000.0).rem_euclid(360.0), 1.0));
        let time = luminary(&mut sample, SUN, 10.0, 1000.0, 0, direction).unwrap();
        let expected = if direction == SearchDirection::Forward {
            1010.0
        } else {
            650.0
        };
        assert!((time - expected).abs() < 1e-8);
    }
}
