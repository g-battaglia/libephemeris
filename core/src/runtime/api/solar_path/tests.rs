use super::*;
#[test]
fn invalid_sampling_step_precedes_main_resolver() {
    let session = Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
    let resolve = |_: &EphemerisContext, _: i32| -> Result<EphemerisContext, Error> {
        panic!("invalid step must precede main")
    };
    for edge in 0..3 {
        let result =
            session.solar_path_walk(2451545.0, 2451546.0, 0.0, FLG_SWIEPH, edge, Some(&resolve));
        assert!(matches!(result, Err(Error::InvalidEclipseSearch { .. })));
    }
}
#[test]
fn reversed_valid_window_is_empty_without_main() {
    let session = Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
    let resolve = |_: &EphemerisContext, _: i32| -> Result<EphemerisContext, Error> {
        panic!("reversed window must not acquire")
    };
    for edge in 0..3 {
        assert_eq!(
            session
                .solar_path_walk(2451546.0, 2451545.0, 1.0, FLG_SWIEPH, edge, Some(&resolve))
                .unwrap(),
            (vec![], vec![], vec![])
        );
    }
}
#[test]
fn width_catch_only_at_known_observer_family() {
    assert!(width_observer_failure(&Error::MathDomain {
        value: f64::INFINITY
    }));
    assert!(!width_observer_failure(&Error::UnknownBody {
        body_id: 999
    }));
    assert!(!width_observer_failure(&Error::ConfigurationMissing {
        missing_config: "fixture",
        suggestion: "fixture"
    }));
}

#[test]
fn besselian_cache_does_not_hold_borrow_during_resolver() {
    let session = Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
    let resolve = |_: &EphemerisContext, _: i32| -> Result<EphemerisContext, Error> {
        Err(Error::UnknownBody { body_id: 777 })
    };
    let native = NativeSolarPath::new(&session, FLG_SWIEPH, Some(&resolve));
    assert_eq!(
        native.besselian(2451545.0),
        Err(Error::UnknownBody { body_id: 777 })
    );
    assert!(native.cache.borrow_mut().is_empty());
}
