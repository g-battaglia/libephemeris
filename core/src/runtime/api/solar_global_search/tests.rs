//! Native public search, physical acquisition and exact epoch-cache checks.
use super::*;

#[test]
fn mask_refusal_precedes_nonfinite_epoch_and_resolver() {
    let session = Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
    let called = std::cell::Cell::new(0);
    let resolve = |_: &EphemerisContext, _| {
        called.set(called.get() + 1);
        Err(Error::UnknownBody { body_id: 999 })
    };
    let error = session
        .sol_eclipse_when_glob_with_position_context(
            f64::NAN,
            0,
            crate::constants::ECL_CENTRAL | crate::constants::ECL_PARTIAL,
            false,
            &resolve,
        )
        .unwrap_err();
    assert!(matches!(
        error.primary_error(),
        Error::InvalidEclipseSearch { .. }
    ));
    assert_eq!(called.get(), 0);
}

#[test]
fn nonfinite_lunation_date_does_not_load_unused_binding_image() {
    let session = Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
    let called = std::cell::Cell::new(0);
    let resolve = |_: &EphemerisContext, _| {
        called.set(called.get() + 1);
        Err(Error::UnknownBody { body_id: 999 })
    };
    let error = session
        .sol_eclipse_when_glob_with_position_context(f64::NAN, 0, 0, false, &resolve)
        .unwrap_err();
    assert!(matches!(error.primary_error(),Error::NonFiniteArgument{value} if value.is_nan()));
    assert_eq!(called.get(), 0);
}

#[test]
fn exact_epoch_memoization_does_not_reuse_a_neighboring_besselian_sample() {
    let first = 2_451_545.0f64;
    let next = f64::from_bits(first.to_bits() + 1);
    let session = Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
    let reads = std::cell::Cell::new(0);
    let resolve = |_: &EphemerisContext, _| {
        reads.set(reads.get() + 1);
        Err(Error::UnknownBody { body_id: 999 })
    };
    let native = NativeSolarGlobal {
        session: &session,
        flags: FLG_SWIEPH,
        resolve: Some(&resolve),
        selected: SelectedMain::new(),
        besselian: std::cell::RefCell::new(std::collections::HashMap::from([(
            first.to_bits(),
            0.25,
        )])),
    };
    assert_eq!(native.besselian_x_er(first).unwrap(), 0.25);
    assert_eq!(reads.get(), 0);
    assert_eq!(
        native.besselian_x_er(next).unwrap_err(),
        Error::UnknownBody { body_id: 999 }
    );
    assert_eq!(reads.get(), 1);
}

#[test]
fn data_boundary_uses_the_primary_native_failure() {
    let range = Error::EphemerisRange {
        requested_jd: 1.0,
        start_jd: 2.0,
        end_jd: 3.0,
        body: Some(SUN),
    };
    assert!(solar::is_ephemeris_boundary(
        &range.clone().with_cause(Error::UnknownBody { body_id: 9 })
    ));
    assert!(!solar::is_ephemeris_boundary(
        &Error::UnknownBody { body_id: 9 }.with_cause(range)
    ));
}

#[test]
fn delayed_besselian_clock_keeps_callback_failure() {
    let mut called = 0;
    let marker = Error::UnknownBody { body_id: 999 };
    let error = solar::besselian_core_lazy((1.0, 0.0, 0.0), (0.002, 0.0001, 0.0), &mut || {
        called += 1;
        Err(marker.clone())
    })
    .unwrap_err();
    assert_eq!(error, marker);
    assert_eq!(called, 1);
}

struct RangeAtPosition(std::cell::Cell<usize>);
impl SolarGlobalEphemeris for RangeAtPosition {
    fn sun_moon_state(&self, _: f64) -> Result<SunMoonState, Error> {
        self.0.set(self.0.get() + 1);
        Err(Error::EphemerisRange {
            requested_jd: 1.,
            start_jd: 2.,
            end_jd: 3.,
            body: Some(SUN),
        })
    }
    fn mean_node_lon_deg(&self, _: f64) -> Result<f64, Error> {
        panic!("unused")
    }
    fn gamma_er(&self, _: f64) -> Result<f64, Error> {
        panic!("unused")
    }
    fn shadow_limits_er(&self, _: f64) -> Result<(f64, f64), Error> {
        panic!("unused")
    }
    fn where_gate(&self, _: f64) -> Result<(i32, f64), Error> {
        panic!("unused")
    }
    fn contact_residuals_km(&self, _: f64) -> Result<(f64, f64, f64), Error> {
        panic!("unused")
    }
    fn besselian_x_er(&self, _: f64) -> Result<f64, Error> {
        panic!("unused")
    }
    fn ground_core_width_km(&self, _: f64) -> Result<f64, Error> {
        panic!("unused")
    }
}
#[test]
fn source_preprobe_catches_only_forward_previous_conjunction() {
    for direction in [
        EclipseSearchDirection::Backward,
        EclipseSearchDirection::Bidirectional,
    ] {
        let eph = RangeAtPosition(std::cell::Cell::new(0));
        let error = solar::sol_eclipse_when_glob_source_via(&eph, 1., 0, direction).unwrap_err();
        assert!(matches!(error, Error::EphemerisRange { .. }));
        assert_eq!(eph.0.get(), 1);
    }
    let eph = RangeAtPosition(std::cell::Cell::new(0));
    let error =
        solar::sol_eclipse_when_glob_source_via(&eph, 1., 0, EclipseSearchDirection::Forward)
            .unwrap_err();
    assert!(matches!(error, Error::EclipseSearchExhausted { .. }));
    assert_eq!(eph.0.get(), 2);
}

#[test]
fn central_classifier_never_consumes_magnitude() {
    let state = SunMoonState::try_new(0., 0., 1., 0., 0., 0., 0.00257, 0.).unwrap();
    let calls = std::cell::Cell::new(0);
    let mut where_role = |magnitude: bool| {
        assert!(!magnitude, "unused local magnitude consumed");
        calls.set(calls.get() + 1);
        Ok((crate::constants::ECL_CENTRAL, 0.))
    };
    let result = crate::events::classify_solar_eclipse_consumed(
        &state,
        0.,
        0.5,
        -0.1,
        &mut where_role,
        &mut || Ok(false),
    )
    .unwrap();
    assert_eq!(
        result.eclipse_type,
        crate::constants::ECL_TOTAL | crate::constants::ECL_CENTRAL
    );
    assert_eq!(calls.get(), 1);
}
#[test]
fn near_miss_classifier_consumes_only_magnitude_role() {
    let state = SunMoonState::try_new(0., 0., 1., 0., 0., 0., 0.00257, 0.).unwrap();
    let calls = std::cell::Cell::new(0);
    let mut where_role = |magnitude: bool| {
        assert!(magnitude, "unrequested retflag consumed");
        calls.set(calls.get() + 1);
        Ok((0, 0.25))
    };
    let result = crate::events::classify_solar_eclipse_consumed(
        &state,
        1.499,
        0.5,
        -0.1,
        &mut where_role,
        &mut || panic!("unused hybrid"),
    )
    .unwrap();
    assert_eq!(
        result.eclipse_type,
        crate::constants::ECL_PARTIAL | crate::constants::ECL_NONCENTRAL
    );
    assert_eq!(calls.get(), 1);
}
#[test]
fn consumed_near_miss_error_is_terminal() {
    let state = SunMoonState::try_new(0., 0., 1., 0., 0., 0., 0.00257, 0.).unwrap();
    let error = crate::events::classify_solar_eclipse_consumed(
        &state,
        1.499,
        0.5,
        -0.1,
        &mut |role| {
            assert!(role);
            Err(Error::UnknownBody { body_id: 998 })
        },
        &mut || panic!("unused hybrid"),
    )
    .unwrap_err();
    assert_eq!(error, Error::UnknownBody { body_id: 998 });
}

/// Original first global search uses a geocentric observer with no surface
/// vector. The exact released source words use DE440, the Espenak-Meeus model,
/// automatic tidal acceleration, disabled IERS and the original 5e-6-day gate.
#[test]
#[ignore = "requires exact DE440 via LIBEPHEMERIS_TEST_JPL"]
fn original_first_global_solar_geocentric_observer() {
    use sha2::{Digest, Sha256};
    let path = std::env::var("LIBEPHEMERIS_TEST_JPL").expect("exact DE440 input");
    let bytes = std::fs::read(&path).expect("read exact DE440");
    assert_eq!(bytes.len(), 119799808);
    let digest: String = Sha256::digest(&bytes)
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect();
    assert_eq!(
        digest,
        "a4ce9bf9b3282becc9f4b2ac3cebe03a2ae7599981aabd7265fd8482fff7c4b5"
    );
    let mut registry = crate::spk::SpkRegistry::new();
    registry.register_image(&path, bytes).unwrap();
    crate::spk::register_all_naif_targets(&mut registry, &path).unwrap();
    let mut ctx = EphemerisContext::new()
        .with_mode(ProviderMode::NativeJpl)
        .with_network_policy(Some(crate::context::NetworkPolicy::Sealed))
        .with_auto_spk_download(Some(false))
        .with_strict_precision(Some(false))
        .with_delta_t_model(Some(crate::time::DeltaTModel::EspenakMeeus))
        .with_iers_delta_t(Some(false))
        .with_user_delta_t(None);
    *ctx.spk_registry_mut() = registry;
    let session = Session::from_context(ctx);
    let (bits, times) = session
        .sol_eclipse_when_glob_direction(
            f64::from_bits(0x4142c4bd570a3d71),
            2,
            0,
            EclipseSearchDirection::Forward,
        )
        .unwrap();
    assert_eq!(bits, 33);
    let source_words = [
        0x4142c4cb56d32406,
        0x4142c4cb54f10462,
        0x4142c4cb4869c1f8,
        0x4142c4cb654681a8,
        0x4142c4cb4df9c48c,
        0x4142c4cb5fb5103e,
        0x4142c4cb4dfa6b56,
        0x4142c4cb5fb3611a,
        0,
        0,
    ];
    for (slot, (actual, expected)) in times.into_iter().zip(source_words).enumerate() {
        let expected = f64::from_bits(expected);
        if expected == 0.0 {
            assert_eq!(actual.to_bits(), expected.to_bits(), "reserved slot {slot}");
        } else {
            assert!(
                (actual - expected).abs() <= 5e-6,
                "original slot {slot}: {actual:?} versus {expected:?}"
            );
        }
    }
}
