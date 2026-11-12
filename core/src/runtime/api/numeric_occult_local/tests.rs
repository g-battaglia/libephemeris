//! Public first-acquisition/owned-error regressions without ephemeris data.
use super::*;
#[test]
fn moon_refusal_precedes_unused_nonfinite_epoch_and_site() {
    let session = Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
    let error = session
        .lun_occult_when_loc(f64::NAN, MOON, [f64::NAN; 3], 2, 0)
        .unwrap_err();
    assert!(matches!(
        error.primary_error(),
        Error::InvalidEclipseSearch {
            entry: "lun_occult",
            ..
        }
    ));
}
#[test]
fn first_target_configuration_failure_precedes_unused_site_domain() {
    let session = Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
    let seen = std::cell::RefCell::new(Vec::new());
    let resolve = |_: &EphemerisContext, body| {
        seen.borrow_mut().push(body);
        Err(Error::ConfigurationMissing {
            missing_config: "actual_local_target_marker",
            suggestion: "retained",
        })
    };
    let error = session
        .lun_occult_when_loc_with_position_context(
            2451545.,
            VENUS,
            [f64::NAN, 100., 0.],
            2,
            ECL_ONE_TRY,
            &resolve,
        )
        .unwrap_err();
    assert_eq!(*seen.borrow(), vec![VENUS]);
    assert!(matches!(
        error,
        Error::ConfigurationMissing {
            missing_config: "actual_local_target_marker",
            ..
        }
    ));
    assert!(session.context().observer().is_none());
}

#[test]
fn finite_disc_radius_rejects_impossible_asin_geometry() {
    assert_eq!(checked_asin(0.0).unwrap(), 0.0);
    assert_eq!(checked_asin(1.0).unwrap(), 90.0);
    assert!(matches!(checked_asin(1.01), Err(Error::MathDomain { .. })));
}
