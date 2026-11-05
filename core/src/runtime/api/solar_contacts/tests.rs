//! Contact-specific acquisition/catch regressions; no astronomical fixtures.
use super::*;

#[test]
fn first_sample_owns_nonfinite_epoch_before_any_bracket_arithmetic() {
    let mut calls = Vec::new();
    let error = contact_source(
        &mut |time| {
            calls.push(time.to_bits());
            Err(Error::InvalidCalcMode {
                value: "marker".into(),
                source: "env",
            })
        },
        f64::NAN,
        false,
    )
    .unwrap_err();
    assert!(matches!(error, Error::InvalidCalcMode { .. }));
    assert_eq!(calls, [f64::NAN.to_bits()]);
}

#[test]
fn numeric_domain_failure_is_an_unresolved_contact_and_data_errors_propagate() {
    assert_eq!(
        contact_source(&mut |_| Err(Error::MathDomain { value: 2.0 }), 1.0, true).unwrap(),
        0.0
    );
    let missing = Error::UnknownBody { body_id: 777 };
    assert_eq!(
        contact_source(&mut |_| Err(missing.clone()), 1.0, true).unwrap_err(),
        missing
    );
}

#[test]
fn range_error_is_not_a_missing_contact() {
    let range = Error::EphemerisRange {
        requested_jd: 1.0,
        start_jd: 2.0,
        end_jd: 3.0,
        body: Some(MOON),
    };
    assert_eq!(
        contact_source(&mut |_| Err(range.clone()), 1.0, false).unwrap_err(),
        range
    );
}

#[test]
fn both_directions_reread_bracket_endpoints_at_the_original_stage() {
    for before in [true, false] {
        let sign = if before { -1.0 } else { 1.0 };
        let mut calls = Vec::new();
        let mut read = |time: f64| {
            calls.push(time);
            Ok((time - 10.0).abs() - 0.01)
        };
        let result = contact_source(&mut read, 10.0, before).unwrap();
        assert_eq!(calls[0], 10.0);
        assert_eq!(calls[1], 10.0 + sign * 0.02);
        let expected = if before {
            [10.0 - 0.02, 10.0]
        } else {
            [10.0, 10.0 + 0.02]
        };
        assert_eq!(&calls[2..4], &expected);
        assert!((result - (10.0 + sign * 0.01)).abs() < 1e-8);
    }
}

#[test]
fn a_caught_trial_failure_keeps_expanding_and_midpoint_failure_returns_zero() {
    let mut samples = Vec::new();
    let value = contact_source(
        &mut |time: f64| {
            samples.push(time);
            if samples.len() == 2 {
                return Err(Error::MathDomain { value: time });
            }
            Ok((time - 10.0).abs() - 0.01)
        },
        10.0,
        false,
    )
    .unwrap();
    assert_eq!(samples[2], 10.0 + 0.02 * 1.6);
    assert!((value - 10.01).abs() < 1e-8);
    let mut calls = 0;
    assert_eq!(
        contact_source(
            &mut |time: f64| {
                calls += 1;
                if calls == 5 {
                    return Err(Error::MathDomain { value: time });
                }
                Ok((time - 10.0).abs() - 0.01)
            },
            10.0,
            false
        )
        .unwrap(),
        0.0
    );
}
