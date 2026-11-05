//! Provider-free station acquisition/order controls; no ephemeris fixtures.
use super::native;
use super::*;

fn state(v: f64) -> [f64; 6] {
    [0.0, 0.0, 1.0, v, 0.0, 0.0]
}

#[test]
fn solar_identities_consume_no_samples() {
    let mut never = |_, _| -> Result<[f64; 6], Error> { panic!("unused position") };
    assert!(!native::retrograde(&mut never, SUN, f64::NAN, -1).unwrap());
    assert_eq!(
        native::motion(&mut never, MOON, f64::NAN, -1).unwrap(),
        StationMotionKind::Direct
    );
    let e = native::find(&mut never, SUN, f64::NAN, StationRequest::Any, -1).unwrap_err();
    assert!(matches!(e, Error::BodyHasNoStations { body: SUN }));
}

#[test]
fn motion_reads_original_previous_next_with_speed() {
    let mut calls = Vec::new();
    let mut sample = |t: f64, flags| {
        calls.push((t, flags));
        Ok(state((2.0 - t) * 0.001))
    };
    assert_eq!(
        native::motion(&mut sample, MARS, 2.0, FLG_XYZ).unwrap(),
        StationMotionKind::StationaryRetrograde
    );
    assert_eq!(
        calls,
        vec![
            (2.0, FLG_XYZ | FLG_SPEED),
            (1.0, FLG_XYZ | FLG_SPEED),
            (3.0, FLG_XYZ | FLG_SPEED)
        ]
    );
}

#[test]
fn configuration_error_is_not_a_predicate_fold() {
    let mut sample = |_, _| {
        Err(Error::ConfigurationMissing {
            missing_config: "fixture",
            suggestion: "fixture",
        })
    };
    assert!(native::retrograde(&mut sample, MARS, 1.0, 0).is_err());
    assert!(native::motion(&mut sample, MARS, 1.0, 0).is_err());
}

#[test]
fn value_error_is_a_predicate_fold_but_not_search_recovery() {
    let mut sample = |_, _| Err(Error::NonFiniteArgument { value: f64::NAN });
    assert!(!native::retrograde(&mut sample, MARS, 1.0, 0).unwrap());
    assert_eq!(
        native::motion(&mut sample, MARS, 1.0, 0).unwrap(),
        StationMotionKind::Direct
    );
    assert!(matches!(
        native::find(&mut sample, MARS, 1.0, StationRequest::Any, 0),
        Err(Error::NonFiniteArgument { .. })
    ));
}

#[test]
fn finder_initial_sample_precedes_scan_sample_and_classification() {
    let mut calls = Vec::new();
    let mut sample = |t: f64, flags| {
        calls.push((t, flags));
        Ok(state(-t))
    };
    let (t, kind) = native::find(&mut sample, MARS, 0.0, StationRequest::Any, 2).unwrap();
    assert_eq!(t, 0.0);
    assert_eq!(kind, StationKind::StationaryRetrograde);
    assert_eq!(&calls[..2], &[(0.0, 258), (0.0, 258)]);
    assert_eq!(&calls[calls.len() - 2..], &[(-1.0, 258), (1.0, 258)]);
}

#[test]
fn full_attempt_catches_classification_error_and_advances_window() {
    let mut failed = false;
    let mut calls = Vec::new();
    let mut sample = |t: f64, flags| {
        calls.push((t, flags));
        if !failed && t == -1.0 {
            failed = true;
            return Err(Error::ConfigurationMissing {
                missing_config: "fixture",
                suggestion: "fixture",
            });
        }
        Ok(state(if t < 390.0 { -t } else { 390.0 - t }))
    };
    let (time, _) = native::find(&mut sample, MARS, 0.0, StationRequest::Any, 2).unwrap();
    assert_eq!(time, 390.0);
    assert!(calls.contains(&(390.0, 258)));
}

#[test]
fn zero_step_reads_both_positions_before_arithmetic_refusal() {
    let mut calls = Vec::new();
    let mut sample = |t, flags| {
        calls.push((t, flags));
        Ok(state(1.0))
    };
    assert!(matches!(
        native::velocity(&mut sample, 5.0, 2, -0.0),
        Err(Error::RationalDivisionByZero)
    ));
    assert_eq!(calls, vec![(5.0, 2), (5.0, 2)]);
}

#[test]
fn negative_velocity_step_and_wrap_preserve_original_flags() {
    let mut calls = Vec::new();
    let mut sample = |t: f64, flags| {
        calls.push((t, flags));
        Ok([
            if t == 6.0 { 1.0 } else { 359.0 },
            t,
            2.0 * t,
            0.0,
            0.0,
            0.0,
        ])
    };
    assert_eq!(
        native::velocity(&mut sample, 5.0, FLG_TRUEPOS, -1.0).unwrap(),
        (1.0, 1.0, 2.0)
    );
    assert_eq!(calls, vec![(6.0, FLG_TRUEPOS), (4.0, FLG_TRUEPOS)]);
}

#[test]
fn summary_double_failure_is_partial_with_first_owned_error() {
    let mut sample = |_, _| Ok(state(1.0));
    let value = native::info(&mut sample, MARS, 0.0, 2).unwrap();
    let StationInformation::Unavailable {
        currently_retrograde,
        velocity,
        error,
    } = value
    else {
        panic!("partial")
    };
    assert!(!currently_retrograde);
    assert_eq!(velocity, 1.0);
    assert!(matches!(
        error,
        Error::NoStationInInterval { planet: MARS, .. }
    ));
}

#[test]
fn summary_selected_longitude_does_not_add_speed() {
    let mut calls = Vec::new();
    let mut sample = |t: f64, flags| {
        calls.push((t, flags));
        Ok([t, 0.0, 1.0, -t, 0.0, 0.0])
    };
    let StationInformation::Complete(info) = native::info(&mut sample, MARS, 0.0, 2).unwrap()
    else {
        panic!("complete")
    };
    assert_eq!(info.jd_station, 0.0);
    assert_eq!(calls.last(), Some(&(0.0, 2)));
}

#[test]
fn aliases_unknown_and_case_preserve_source_request() {
    assert_eq!(station_request("rEtRoGrAdE"), StationRequest::Retrograde);
    assert_eq!(station_request("sD"), StationRequest::Direct);
    assert_eq!(station_request("unknown"), StationRequest::Any);
}

#[test]
fn public_solar_identities_do_not_resolve_native_main_configuration() {
    let session = Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
    let resolve = |_: &EphemerisContext, _| -> Result<EphemerisContext, Error> {
        panic!("unused configuration")
    };
    let value = session
        .station_query_with_position_context(
            StationQuery::Retrograde {
                body: SUN,
                jd_ut: f64::NAN,
                flags: -1,
            },
            &resolve,
        )
        .unwrap();
    assert_eq!(value, StationOutput::Retrograde(false));
}
