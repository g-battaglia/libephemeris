//! Native house provider lifetime and explicit clock controls without SPK data.
use super::*;

#[test]
fn native_houses_provider_construction_does_not_consume_clock_site_or_alias() {
    let mut ctx = EphemerisContext::default().with_mode(ProviderMode::NativeJpl);
    ctx.set_sidereal_mode(SIDM_USER | SIDBIT_USER_UT, f64::MAX, 24.0);
    let session = Session::from_context(ctx);
    session
        .with_native_houses_providers((f64::NAN, f64::NAN, f64::NAN), |_, _, _, gauq| {
            assert!(gauq.rise_set_at.is_some());
            Ok(())
        })
        .unwrap();
    assert!(session.context().spk_registry().is_empty());
    assert_eq!(session.context().observer(), None);
}

#[test]
fn native_houses_defining_epoch_callbacks_keep_configured_model() {
    let mut ctx = EphemerisContext::default()
        .with_mode(ProviderMode::NativeJpl)
        .with_user_delta_t(Some(0.25));
    ctx.set_sidereal_mode(SIDM_USER | SIDBIT_USER_UT, 2440587.5, 24.0);
    let session = Session::from_context(ctx);
    session
        .with_native_houses_providers((f64::NAN, f64::NAN, f64::NAN), |_, _, sid, _| {
            assert_eq!((sid.user_t0_tt_at)()?, 2440587.75);
            assert_eq!(
                (sid.valens_epoch_tt_at)()?,
                crate::sidereal::VALENS_MOON_T0_UT + 0.25
            );
            Ok(())
        })
        .unwrap();
}

#[test]
fn native_houses_plain_and_extended_ayanamsha_keep_distinct_public_clocks() {
    let mut ctx = EphemerisContext::default()
        .with_mode(ProviderMode::NativeJpl)
        .with_user_delta_t(Some(0.25));
    ctx.set_sidereal_mode(SIDM_USER | SIDBIT_USER_UT, 2440587.5, 24.0);
    let session = Session::from_context(ctx);
    session
        .with_native_houses_providers((0.0, 0.0, 0.0), |selection, _, sid, _| {
            let delta = |_| Ok(0.25);
            let aya = SealedAyanamshaServing::new(selection, &delta, None);
            let epoch = 2451545.125;
            let entry_tt = crate::time::default_timescale_tt_from_ut1(epoch)?;
            let helper_ut = crate::time::default_timescale_ut1_from_tt(entry_tt)?;
            let expected_ex = crate::planetary_moons::moon_spk_ayanamsha(helper_ut, 0, &aya)?;
            let expected_mean = aya.mean_tt_with(entry_tt, false, false)?.0;
            assert_eq!(
                (sid.ayanamsha_ex_at)(epoch, 0)?.1.to_bits(),
                expected_ex.to_bits()
            );
            assert_eq!(
                (sid.mean_ut_plain_at)(epoch)?.to_bits(),
                expected_mean.to_bits()
            );
            // The configured quarter-day offset defines USER's epoch only.
            assert_ne!(expected_mean.to_bits(), aya.mean_ut(epoch)?.0.to_bits());
            Ok(())
        })
        .unwrap();
}

#[test]
fn native_houses_anchor_range_wrapper_uses_original_request_and_absent_body() {
    let registry = crate::spk::SpkRegistry::new();
    let provider =
        crate::positions::LocalSpkStateProvider::new(&registry, ProviderMode::NativeJpl).unwrap();
    let original = 2451545.125;
    let error = Error::SpkOutOfRange {
        target: 599,
        et: -1.0,
        start: 100.0,
        end: 200.0,
    };
    // With no main image this synthetic control uses the existing raw span
    // fallback. Actual source/native acceptance must bind the main DE image.
    assert_eq!(
        native_anchor_range_error(&provider, original, error),
        Error::EphemerisRange {
            requested_jd: original,
            start_jd: 2451545.0 + 100.0 / 86400.0,
            end_jd: 2451545.0 + 200.0 / 86400.0,
            body: None,
        }
    );
}

#[test]
fn native_houses_anchor_range_wrapper_keeps_missing_alias_identity() {
    let registry = crate::spk::SpkRegistry::new();
    let provider =
        crate::positions::LocalSpkStateProvider::new(&registry, ProviderMode::NativeJpl).unwrap();
    let error = Error::SpkUnregisteredBody { body: SUN };
    assert_eq!(
        native_anchor_range_error(&provider, 2451545.125, error.clone()),
        error
    );
}

#[test]
fn native_public_houses_do_not_consume_unselected_hostile_user_epoch() {
    let mut ctx = EphemerisContext::default()
        .with_mode(ProviderMode::NativeJpl)
        .with_user_delta_t(Some(0.0008));
    ctx.set_sidereal_mode(SIDM_USER | SIDBIT_USER_UT, f64::MAX, 24.0);
    let session = Session::from_context(ctx);
    let actual = session
        .houses_ex(2451545.125, 0, 52.52, 13.405, b'P')
        .unwrap();
    let expected = session
        .with_native_houses_providers((f64::NAN, f64::NAN, f64::NAN), |selection, jd, sid, _| {
            crate::houses::houses_ex_via_callbacks(
                2451545.125,
                52.52,
                13.405,
                'P',
                0,
                selection,
                jd,
                sid,
            )
        })
        .unwrap();
    assert_eq!(actual.cusps.as_slice(), expected.0.as_slice());
    assert_eq!(actual.angles, expected.1);
    assert_eq!(session.context().observer(), None);
    assert!(session.context().spk_registry().is_empty());
}

#[test]
fn native_houses_lazy_position_context_is_unconsumed_by_tropical_frames() {
    let mut ctx = EphemerisContext::default()
        .with_mode(ProviderMode::NativeJpl)
        .with_user_delta_t(Some(0.0008));
    ctx.set_sidereal_mode(SIDM_USER | SIDBIT_USER_UT, f64::MAX, 24.0);
    let session = Session::from_context(ctx);
    let resolve = |_: &EphemerisContext, _: i32| -> Result<EphemerisContext, Error> {
        panic!("ordinary tropical frames must not consume position configuration")
    };
    let actual = session
        .with_native_houses_providers_with_position_context(
            (f64::NAN, f64::NAN, f64::NAN),
            &resolve,
            |selection, jd, sid, gauq| {
                assert!(gauq.rise_set_at.is_some());
                crate::houses::houses_ex_via_callbacks(
                    2451545.125,
                    52.52,
                    13.405,
                    'P',
                    0,
                    selection,
                    jd,
                    sid,
                )
            },
        )
        .unwrap();
    let expected = session.houses(2451545.125, 52.52, 13.405, b'P').unwrap();
    assert_eq!(actual.0.as_slice(), expected.cusps.as_slice());
    assert_eq!(actual.1, expected.angles);
    assert!(session.context().spk_registry().is_empty());
    assert_eq!(session.context().observer(), None);
}

#[test]
fn native_houses_consumed_sun_and_body_resolve_before_position_clock() {
    let session = Session::from_context(
        EphemerisContext::default()
            .with_mode(ProviderMode::NativeJpl)
            .with_user_delta_t(Some(0.25)),
    );
    let calls = std::cell::RefCell::new(Vec::new());
    let marker = Error::ConfigurationMissing {
        missing_config: "test_position_configuration",
        suggestion: "exact resolver error must propagate",
    };
    let resolve = |ctx: &EphemerisContext, body| {
        assert_eq!(ctx, session.context());
        calls.borrow_mut().push(body);
        Err(marker.clone())
    };
    session
        .with_native_houses_providers_with_position_context(
            (f64::NAN, f64::NAN, f64::NAN),
            &resolve,
            |_, jd, _, gauq| {
                assert_eq!(
                    (jd.sun_declination_at)(f64::MAX, FLG_JPLEPH),
                    Err(marker.clone())
                );
                assert_eq!(
                    (gauq.body_ecliptic_at)(f64::MAX, MARS, FLG_SPEED),
                    Err(marker.clone())
                );
                assert_eq!(
                    (gauq.body_ecliptic_at)(f64::MAX, ECL_NUT, 0),
                    Err(marker.clone())
                );
                assert_eq!(
                    (gauq.body_ecliptic_at)(f64::MAX, 999, 0),
                    Err(marker.clone())
                );
                assert_eq!(
                    (gauq.body_ecliptic_at)(f64::MAX, EARTH, 0),
                    Err(marker.clone())
                );
                Ok(())
            },
        )
        .unwrap();
    // The boundary resolver sees exact public identities; its own standard-
    // body predicate decides which ones consume process configuration.
    assert_eq!(*calls.borrow(), [SUN, MARS, ECL_NUT, 999, EARTH]);
    assert_eq!(session.context().observer(), None);
}

#[test]
fn native_houses_consumed_site_precedes_position_configuration() {
    let session =
        Session::from_context(EphemerisContext::default().with_mode(ProviderMode::NativeJpl));
    let resolve = |_: &EphemerisContext, _: i32| -> Result<EphemerisContext, Error> {
        panic!("invalid consumed observer must fail before configuration")
    };
    session
        .with_native_houses_providers_with_position_context(
            (13.405, 91.0, 34.0),
            &resolve,
            |_, _, _, gauq| {
                assert!(matches!(
                    (gauq.body_ecliptic_at)(f64::MAX, MARS, FLG_TOPOCTR),
                    Err(Error::InvalidObserver { .. })
                ));
                Ok(())
            },
        )
        .unwrap();
    assert_eq!(session.context().observer(), None);
}

#[test]
fn native_houses_position_resolver_receives_only_temporary_observer_snapshot() {
    let caller = Observer::new(-70.0, -30.0, 100.0).unwrap();
    let site = Observer::new(13.405, 52.52, 34.0).unwrap();
    let session = Session::from_context(
        EphemerisContext::default()
            .with_mode(ProviderMode::NativeJpl)
            .with_observer(caller),
    );
    let seen = std::cell::RefCell::new(Vec::new());
    let marker = Error::ConfigurationMissing {
        missing_config: "test_observer_snapshot",
        suggestion: "no position sample is needed for this ownership control",
    };
    let resolve = |ctx: &EphemerisContext, body| {
        assert_eq!(body, MARS);
        seen.borrow_mut().push(ctx.observer());
        Err(marker.clone())
    };
    session
        .with_native_houses_providers_with_position_context(
            (13.405, 52.52, 34.0),
            &resolve,
            |_, _, _, gauq| {
                assert_eq!(
                    (gauq.body_ecliptic_at)(2451545.125, MARS, FLG_TOPOCTR),
                    Err(marker.clone())
                );
                assert_eq!(
                    (gauq.body_ecliptic_at)(2451545.125, MARS, 0),
                    Err(marker.clone())
                );
                Ok(())
            },
        )
        .unwrap();
    assert_eq!(*seen.borrow(), [Some(site), Some(caller)]);
    assert_eq!(session.context().observer(), Some(caller));
}

/// Original scalar event request, without any fabricated astronomical output.
fn native_sector_event_query(body_id: i32, modifiers: i32) -> crate::houses::GauquelinRiseSetQuery {
    crate::houses::GauquelinRiseSetQuery {
        jd_start_ut: 2451545.125,
        body_id,
        rsmi: CALC_RISE | modifiers,
        lon_deg_east: 13.405,
        lat_deg: 52.52,
        height_m: 34.0,
        atpress_mbar: 1013.25,
        attemp_c: 15.0,
        flags: FLG_JPLEPH | FLG_RADIANS,
    }
}

#[test]
fn native_gauquelin_event_feed_preserves_actual_coordinate_error_before_resolver() {
    let session =
        Session::from_context(EphemerisContext::default().with_mode(ProviderMode::NativeJpl));
    let resolve = |_: &EphemerisContext, _: i32| -> Result<EphemerisContext, Error> {
        panic!("invalid ordinary event observer precedes every consuming position call")
    };
    session
        .with_native_houses_providers_with_position_context(
            (13.405, 91.0, 34.0),
            &resolve,
            |_, _, _, gauq| {
                for method in 2..=5 {
                    assert_eq!(
                        crate::houses::gauquelin_sector_via(
                            2451545.125,
                            MARS,
                            method,
                            91.0,
                            13.405,
                            0,
                            34.0,
                            0.0,
                            0.0,
                            gauq,
                        ),
                        Err(Error::InvalidObserver {
                            field: "lat_deg",
                            value: 91.0
                        })
                    );
                }
                Ok(())
            },
        )
        .unwrap();
    assert_eq!(session.context().observer(), None);
}

#[test]
fn native_gauquelin_event_feed_consumes_temporary_site_and_keeps_configured_context() {
    let caller = Observer::new(-70.0, -30.0, 100.0).unwrap();
    let site = Observer::new(13.405, 52.52, 34.0).unwrap();
    let session = Session::from_context(
        EphemerisContext::default()
            .with_mode(ProviderMode::NativeJpl)
            .with_user_delta_t(Some(0.25))
            .with_observer(caller),
    );
    let calls = std::cell::Cell::new(0);
    let marker = Error::ConfigurationMissing {
        missing_config: "native_event_configuration",
        suggestion: "unchanged owning callback error",
    };
    let resolve = |ctx: &EphemerisContext, body| {
        calls.set(calls.get() + 1);
        assert_eq!(body, MARS);
        assert_eq!(ctx, &session.context().clone().with_observer(site));
        Err(marker.clone())
    };
    let query = native_sector_event_query(MARS, BIT_DISC_CENTER | BIT_NO_REFRACTION);
    assert_eq!(
        native_gauquelin_rise_set(&session, query, Some(&resolve)),
        Err(marker)
    );
    assert_eq!(calls.get(), 1);
    assert_eq!(session.context().observer(), Some(caller));
    assert!(session.context().spk_registry().is_empty());
}

#[test]
fn native_gauquelin_event_feed_preserves_first_probe_range_body_day_and_bounds() {
    let session =
        Session::from_context(EphemerisContext::default().with_mode(ProviderMode::NativeJpl));
    let marker = Error::EphemerisRange {
        requested_jd: 2451545.125,
        start_jd: 2450000.0,
        end_jd: 2451000.0,
        body: Some(999),
    };
    let calls = std::cell::Cell::new(0);
    let resolve = |_: &EphemerisContext, body| {
        assert_eq!(body, 999);
        calls.set(calls.get() + 1);
        Err(marker.clone())
    };
    assert_eq!(
        native_gauquelin_rise_set(
            &session,
            native_sector_event_query(999, BIT_NO_REFRACTION),
            Some(&resolve),
        ),
        Err(marker)
    );
    assert_eq!(calls.get(), 1);
}

#[test]
fn native_gauquelin_event_feed_undefined_earth_precedes_unused_config_and_clock() {
    let session =
        Session::from_context(EphemerisContext::default().with_mode(ProviderMode::NativeJpl));
    let resolve = |_: &EphemerisContext, _: i32| -> Result<EphemerisContext, Error> {
        panic!("Earth direction is undefined without an image or event clock")
    };
    let mut query = native_sector_event_query(EARTH, BIT_DISC_CENTER);
    query.jd_start_ut = f64::MAX;
    assert_eq!(
        native_gauquelin_rise_set(&session, query, Some(&resolve)),
        Err(Error::InvalidBody {
            body_id: EARTH,
            operation: "rise_trans"
        })
    );
}

#[test]
fn native_gauquelin_event_feed_forwards_consumed_atmosphere_without_shape_defaults() {
    let session =
        Session::from_context(EphemerisContext::default().with_mode(ProviderMode::NativeJpl));
    let resolve = |_: &EphemerisContext, _: i32| -> Result<EphemerisContext, Error> {
        panic!("nonfinite event input precedes position configuration")
    };
    let mut query = native_sector_event_query(MARS, BIT_DISC_CENTER);
    query.attemp_c = f64::INFINITY;
    assert_eq!(
        native_gauquelin_rise_set(&session, query, Some(&resolve)),
        Err(Error::NonFiniteArgument {
            value: f64::INFINITY
        })
    );
}
