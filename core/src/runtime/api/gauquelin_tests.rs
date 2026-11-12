//! Selection/ownership/catch regressions; original numeric acceptance separate.
use super::*;
use std::cell::RefCell;
#[test]
fn frame_wins_before_place_or_topocentric_site() {
    let error = star_sector_via(
        "missing",
        2451545.0,
        0,
        [999.0, 999.0, f64::NAN],
        0.0,
        0.0,
        FLG_TOPOCTR,
        &|_| {
            Err(Error::ConfigurationMissing {
                missing_config: "actual-frame",
                suggestion: "actual-frame",
            })
        },
        &|_, _| panic!("unselected place"),
        &|_| panic!("unselected event"),
    )
    .unwrap_err();
    assert!(matches!(
        error,
        Error::ConfigurationMissing {
            missing_config: "actual-frame",
            ..
        }
    ));
}
#[test]
fn scalar_star_uses_speed_and_method1_zeros_only_body_latitude() {
    let events = RefCell::new(Vec::new());
    let value = star_sector_via(
        "Sirius",
        2451545.0,
        1,
        [0.0, 0.0, f64::NAN],
        f64::NAN,
        f64::NAN,
        FLG_J2000,
        &|jd| {
            events.borrow_mut().push(("frame", jd.to_bits(), 0));
            Ok((150.0, 23.44))
        },
        &|jd, flags| {
            events.borrow_mut().push(("place", jd.to_bits(), flags));
            Ok((100.0, 1000.0))
        },
        &|_| panic!("method1 event feed"),
    )
    .unwrap();
    assert!(value.is_finite());
    assert_eq!(
        &events.borrow()[..],
        &[
            ("frame", 2451545.0_f64.to_bits(), 0),
            ("place", 2451545.0_f64.to_bits(), FLG_J2000 | FLG_SPEED)
        ]
    );
}
#[test]
fn method6_has_no_rise_or_set_feed() {
    let calls = RefCell::new(Vec::new());
    let error = star_sector_via(
        "Sirius",
        2451545.0,
        6,
        [0.0, 0.0, 0.0],
        0.0,
        0.0,
        0,
        &|_| {
            calls.borrow_mut().push("frame");
            Err(Error::ZeroVectorNormalize)
        },
        &|_, _| panic!("after failed frame"),
        &|_| panic!("method6 selected events"),
    )
    .unwrap_err();
    assert!(matches!(error, Error::ZeroVectorNormalize));
    assert_eq!(&calls.borrow()[..], &["frame"]);
}
#[test]
fn event_methods_keep_exact_query_flags_atmosphere_and_lookbacks() {
    for method in 2..=5 {
        let queries = RefCell::new(Vec::new());
        let value = star_sector_via(
            "Sirius",
            1000.0,
            method,
            [12.5, 48.0, 100.0],
            0.0,
            0.0,
            FLG_XYZ,
            &|_| panic!("event frame"),
            &|_, _| panic!("event scalar place"),
            &|query| {
                queries.borrow_mut().push(query);
                let index = queries.borrow().len();
                Ok((0, [1000.75, 1000.25, 999.75, 999.25][index - 1]))
            },
        )
        .unwrap();
        assert_eq!(value, 10.0);
        let queries = queries.borrow();
        assert_eq!(queries.len(), 4);
        assert_eq!(
            queries.iter().map(|q| q.jd_start_ut).collect::<Vec<_>>(),
            vec![1000.0, 1000.0, 999.25, 998.75]
        );
        let bits = (if method == 2 || method == 3 {
            BIT_DISC_CENTER
        } else {
            0
        }) | (if method == 2 || method == 4 {
            BIT_NO_REFRACTION
        } else {
            0
        });
        for (index, q) in queries.iter().enumerate() {
            assert_eq!(
                q.rsmi,
                (if index % 2 == 0 { CALC_RISE } else { CALC_SET }) | bits
            );
            assert_eq!(
                (q.atpress_mbar, q.attemp_c, q.flags),
                (1013.25, 15.0, FLG_XYZ)
            );
            assert_eq!((q.lon_deg_east, q.lat_deg, q.height_m), (12.5, 48.0, 100.0));
        }
    }
}

#[test]
fn ordinary_lookup_range_and_coordinate_failures_are_not_caught() {
    let errors = [
        Error::StarNotFound {
            star_id: "missing".into(),
            search_type: Some("name"),
        },
        Error::EphemerisRange {
            requested_jd: 1.0,
            start_jd: 2.0,
            end_jd: 3.0,
            body: None,
        },
        Error::InvalidObserver {
            field: "lat_deg",
            value: 99.0,
        },
    ];
    for original in errors {
        let error = star_sector_via(
            "missing",
            1000.0,
            4,
            [0.0, 99.0, 0.0],
            0.0,
            0.0,
            0,
            &|_| panic!(),
            &|_, _| panic!(),
            &|_| Err(original.clone()),
        )
        .unwrap_err();
        assert_eq!(error, original);
    }
}
#[test]
fn genuine_session_no_event_preparation_propagates_lookup_without_image() {
    let session = Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
    let error = session
        .star_gauquelin_sector(
            "not_a_catalog_star",
            2451545.0,
            2,
            [0.0, 0.0, 0.0],
            0.0,
            0.0,
            0,
        )
        .unwrap_err();
    assert!(matches!(error.primary_error(), Error::StarNotFound { .. }));
    assert!(session.context().observer().is_none());
}
#[test]
fn genuine_scalar_call_temporary_site_does_not_replace_retained_observer() {
    let observer = Observer::new(12.0, 34.0, 56.0).unwrap();
    let session = Session::from_context(EphemerisContext::new().with_observer(observer));
    let error = session
        .star_gauquelin_sector(
            "not_a_catalog_star",
            2451545.0,
            0,
            [3.0, 4.0, 5.0],
            0.0,
            0.0,
            FLG_TOPOCTR,
        )
        .unwrap_err();
    assert!(matches!(error.primary_error(), Error::StarNotFound { .. }));
    assert_eq!(session.context().observer(), Some(observer));
}

#[test]
fn configured_house_frame_first_clock_error_is_terminal() {
    let error = crate::time::public_sidereal::house_frame_at(2451545.0, &|_| {
        Err(Error::ConfigurationMissing {
            missing_config: "first-house-model-clock",
            suggestion: "first-house-model-clock",
        })
    })
    .unwrap_err();
    assert!(matches!(
        error,
        Error::ConfigurationMissing {
            missing_config: "first-house-model-clock",
            ..
        }
    ));
}
#[test]
fn configured_house_frame_modern_sidereal_reconsumes_original_ut_clock() {
    let calls = RefCell::new(Vec::new());
    let result = crate::time::public_sidereal::house_frame_at(2451545.0, &|ut| {
        calls.borrow_mut().push(ut.to_bits());
        Ok(0.0008)
    })
    .unwrap();
    assert_eq!(
        &calls.borrow()[..],
        &[2451545.0_f64.to_bits(), 2451545.0_f64.to_bits()]
    );
    assert!(result.0.is_finite() && result.1.is_finite());
}

#[test]
fn star_defining_observer_borrows_same_registry_and_keeps_original_site() {
    let original = Observer::new(12.0, 34.0, 56.0).unwrap();
    let selected = Observer::new(3.0, 4.0, 5.0).unwrap();
    let session = Session::from_context(
        EphemerisContext::new()
            .with_mode(ProviderMode::NativeJpl)
            .with_observer(original),
    );
    let pointer = session.context().spk_registry() as *const _;
    let error = session
        .star_ut_at_observer(
            "not_a_catalog_star",
            2451545.0,
            FLG_TOPOCTR | FLG_SPEED,
            selected,
        )
        .unwrap_err();
    assert!(matches!(error.primary_error(), Error::StarNotFound { .. }));
    assert_eq!(session.context().spk_registry() as *const _, pointer);
    assert_eq!(session.context().observer(), Some(original));
}
