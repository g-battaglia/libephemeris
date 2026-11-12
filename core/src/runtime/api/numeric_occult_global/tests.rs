//! Native target/filter order and owned resolver contracts without any asset.
use super::*;

#[test]
fn global_moon_self_is_rejected_before_nonfinite_epoch_or_mask() {
    let session = Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
    let error = session
        .lun_occult_when_glob(f64::NAN, MOON, 2, ECL_PARTIAL | ECL_CENTRAL, 0)
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
fn impossible_mask_does_not_acquire_main_or_call_resolver() {
    let session = Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
    let called = std::cell::Cell::new(0);
    let resolve = |_: &EphemerisContext, _: i32| {
        called.set(called.get() + 1);
        Err(Error::UnknownBody { body_id: 999 })
    };
    let error = session
        .lun_occult_when_glob_with_position_context(f64::NAN, MARS, 2, ECL_ANNULAR, 0, &resolve)
        .unwrap_err();
    assert_eq!(called.get(), 0);
    assert!(matches!(
        error.primary_error(),
        Error::InvalidEclipseSearch { .. }
    ));
    assert!(session.context().spk_registry().is_empty());
}

#[test]
fn first_target_resolver_error_is_returned_without_a_no_event_substitute() {
    let session = Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
    let seen = std::cell::RefCell::new(Vec::new());
    let resolve = |_: &EphemerisContext, body: i32| {
        seen.borrow_mut().push(body);
        Err(Error::ConfigurationMissing {
            missing_config: "actual_resolver_marker",
            suggestion: "retained",
        })
    };
    let error = session
        .lun_occult_when_glob_with_position_context(2451545., VENUS, 2, 0, ECL_ONE_TRY, &resolve)
        .unwrap_err();
    assert_eq!(*seen.borrow(), vec![VENUS]);
    assert!(matches!(
        error,
        Error::ConfigurationMissing {
            missing_config: "actual_resolver_marker",
            ..
        }
    ));
}

/// A first coverage boundary consumes the real checked event walk, then the
/// same public wrapper used by the non-Native serving branch. No image or
/// provider is loaded; every later acquisition is an unconsumed marker.
#[test]
fn checked_first_boundary_exhausts_search_without_lunar_acquisition() {
    struct Boundary(std::cell::Cell<usize>);
    impl OccultGlobalEphemeris for Boundary {
        fn body_lonlat_deg(&self, jd: f64) -> Result<(f64, f64), Error> {
            self.0.set(self.0.get() + 1);
            Err(Error::LebOutOfRange {
                body_id: VENUS,
                jd,
                start: jd + 1.0,
                end: jd + 2.0,
            })
        }
        fn moon_lonlat_deg(&self, _: f64) -> Result<(f64, f64), Error> {
            panic!("Moon unconsumed")
        }
        fn sep_minus_radii_deg(&self, _: f64) -> Result<f64, Error> {
            panic!("objective unconsumed")
        }
        fn where_core(&self, _: f64) -> Result<WhereGeometry, Error> {
            panic!("geometry unconsumed")
        }
        fn ra_gap_deg(&self, _: f64) -> Result<f64, Error> {
            panic!("RA unconsumed")
        }
    }
    let boundary = Boundary(std::cell::Cell::new(0));
    let result = event::lun_occult_when_glob_via(
        &boundary,
        2451545.0,
        Some(VENUS),
        0,
        false,
        "3",
        false,
        false,
    );
    let raw = result.as_ref().unwrap_err();
    assert!(matches!(
        raw,
        Error::EclipseSearchExhausted {
            entry: "lun_occult_when_glob",
            ..
        }
    ));
    assert_eq!(boundary.0.get(), 1);
}
