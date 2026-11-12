//! Synthetic callback/guard controls; no ephemeris image or source process.
use super::*;

fn request(rsmi: i32) -> RiseSetRequest {
    RiseSetRequest {
        jd_start_ut: 2_451_545.0,
        lon_deg_east: 12.5,
        lat_deg: 41.9,
        height_m: 0.0,
        atpress_mbar: 1_013.25,
        attemp_c: 15.0,
        horhgt_deg: 0.0,
        rsmi,
    }
}

#[test]
fn stars_are_centres_without_planet_disc_lift() {
    for bits in [0, BIT_DISC_CENTER, BIT_DISC_BOTTOM, BIT_FIXED_DISC_SIZE] {
        let p = star_event_parameters(&request(1 | bits)).unwrap();
        assert_eq!(p.disc, events::DiscPoint::Center);
        assert_eq!(p.radius.to_bits(), 0.0_f64.to_bits());
        assert_eq!(p.lift, 0.75);
    }
    assert_eq!(
        star_event_parameters(&request(1 | BIT_NO_REFRACTION))
            .unwrap()
            .lift,
        0.0
    );
}

#[test]
fn stellar_twilight_uses_source_precedence_and_unrefracted_centre() {
    for (bits, horizon) in [
        (
            BIT_CIVIL_TWILIGHT | BIT_NAUTIC_TWILIGHT | BIT_ASTRO_TWILIGHT,
            -6.0,
        ),
        (BIT_NAUTIC_TWILIGHT | BIT_ASTRO_TWILIGHT, -12.0),
        (BIT_ASTRO_TWILIGHT, -18.0),
    ] {
        let mut req = request(1 | bits);
        req.horhgt_deg = 3.0;
        let p = star_event_parameters(&req).unwrap();
        assert_eq!(p.horizon, horizon);
        assert_ne!(p.effective & BIT_NO_REFRACTION, 0);
        assert_ne!(p.effective & BIT_DISC_CENTER, 0);
        assert_eq!(p.lift.to_bits(), 0.0_f64.to_bits());
    }
}

#[test]
fn star_callback_preserves_ut_and_fixed_public_flags() {
    let ut = f64::from_bits(0x4142_b42c_8000_0001);
    for (internal, actual) in [
        (FLG_TOPOCTR | FLG_SPEED, FLG_SPEED),
        (FLG_EQUATORIAL | FLG_SPEED, FLG_EQUATORIAL | FLG_SPEED),
        (FLG_SPEED, FLG_SPEED),
    ] {
        let value = star_position_callback(ut, internal, &|seen, flags| {
            assert_eq!(seen.to_bits(), ut.to_bits());
            assert_eq!(flags, actual);
            Ok([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        })
        .unwrap();
        assert_eq!(value, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    }
}

#[test]
fn star_callback_propagates_actual_typed_refusal_without_retry() {
    let calls = std::sync::Mutex::new(0);
    let error = star_position_callback(2_451_545.0, FLG_SPEED, &|ut, _| {
        *calls.lock().unwrap() += 1;
        Err(Error::NonFiniteArgument { value: ut })
    })
    .unwrap_err();
    assert!(matches!(error, Error::NonFiniteArgument { value } if value == 2_451_545.0));
    assert_eq!(*calls.lock().unwrap(), 1);
}

#[test]
fn stellar_transit_consumes_equatorial_place_then_gast_only() {
    let mut req = request(4 | BIT_GEOCTR_NO_ECL_LAT);
    req.lon_deg_east = 0.0;
    let p = star_event_parameters(&req).unwrap();
    let order = std::sync::Mutex::new(Vec::new());
    let output = native_star_event_callbacks(
        -1,
        &req,
        &p,
        &|ut, internal| {
            star_position_callback(ut, internal, &|_, flags| {
                assert_eq!(flags, FLG_EQUATORIAL | FLG_SPEED);
                order.lock().unwrap().push("place");
                Ok([0.0; 6])
            })
        },
        &|_, _| panic!("transit consumes no altitude"),
        &|_, _| panic!("transit ignores Hindu projection"),
        &|_| {
            order.lock().unwrap().push("gast");
            Ok(0.0)
        },
    )
    .unwrap();
    assert_eq!(output.0, 0);
    assert_eq!(output.1[0], req.jd_start_ut);
    assert!(output.1[1..].iter().all(|x| x.to_bits() == 0));
    assert_eq!(&order.lock().unwrap()[..2], &["place", "gast"]);
}

#[test]
fn ordinary_star_circumpolar_precheck_never_requests_topocentre() {
    let mut req = request(1);
    req.lat_deg = 89.0;
    let p = star_event_parameters(&req).unwrap();
    let output = native_star_event_callbacks(
        -1,
        &req,
        &p,
        &|ut, internal| {
            star_position_callback(ut, internal, &|_, flags| {
                assert_eq!(flags, FLG_EQUATORIAL | FLG_SPEED);
                Ok([0.0, 80.0, 1.0, 0.0, 0.0, 0.0])
            })
        },
        &|_, _| panic!("circumpolar precheck returns before altitude"),
        &|_, _| panic!("ordinary precheck does not project"),
        &|_| panic!("rise/set does not consume transit GAST"),
    )
    .unwrap();
    assert_eq!(output.0, -2);
    assert!(output.1.iter().all(|x| x.to_bits() == 0));
}

#[test]
fn hindu_star_precheck_uses_speed_and_projected_declination() {
    let mut req = request(1 | BIT_GEOCTR_NO_ECL_LAT);
    req.lat_deg = 89.0;
    let p = star_event_parameters(&req).unwrap();
    let output = native_star_event_callbacks(
        -1,
        &req,
        &p,
        &|ut, internal| {
            star_position_callback(ut, internal, &|_, flags| {
                assert_eq!(flags, FLG_SPEED);
                Ok([30.0, 50.0, 1.0, 0.0, 0.0, 0.0])
            })
        },
        &|_, _| panic!("projected circumpolar precheck returns first"),
        &|_, place| {
            assert_eq!(place, [30.0, 0.0, 1.0]);
            Ok((0.0, 80.0))
        },
        &|_| panic!("rise/set does not consume transit GAST"),
    )
    .unwrap();
    assert_eq!(output.0, -2);
    assert!(output.1.iter().all(|x| x.to_bits() == 0));
}

#[test]
fn public_star_event_ignores_caller_flags_and_preserves_catalog_refusal() {
    let session = Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
    for flags in [0, FLG_XYZ | FLG_RADIANS, FLG_TOPOCTR, -1] {
        let result = session.star_rise_trans(
            "No Such Star xyz",
            2_451_545.0,
            1,
            [12.5, 41.9, 0.0],
            1_013.25,
            15.0,
            flags,
        );
        assert!(matches!(
            result.as_ref().err().map(Error::primary_error),
            Some(Error::StarNotFound { .. })
        ));
    }
    let result = session.star_rise_trans(
        "No Such Star xyz",
        f64::NAN,
        1,
        [12.5, 41.9, 0.0],
        1_013.25,
        15.0,
        0,
    );
    assert!(matches!(
        result.as_ref().err().map(Error::primary_error),
        Some(Error::StarNotFound { .. })
    ));
}

#[test]
fn pressure_power_refusal_precedes_catalog_and_owns_actual_operands() {
    let session = Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
    let error = session
        .star_rise_trans(
            "No Such Star xyz",
            2451545.0,
            CALC_MTRANSIT,
            [0.0, 0.0, -1e100],
            0.0,
            0.0,
            0,
        )
        .unwrap_err();
    assert!(matches!(
        error.primary_error(),
        Error::NumericalOverflow {
            operation: "barometric pressure",
            ..
        }
    ));
}

#[test]
fn pressure_infinite_input_nan_and_multiplication_are_distinct_stages() {
    assert!(star_pressure_from_altitude_m(f64::NAN).unwrap().is_nan());
    assert!(
        star_pressure_from_altitude_m(f64::NEG_INFINITY)
            .unwrap()
            .is_infinite()
    );
    assert_eq!(star_pressure_from_altitude_m(f64::INFINITY).unwrap(), 0.0);
    let base = 4e58_f64;
    let height = (1.0 - base) * 288.0 / 0.0065;
    let actual_base = 1.0 - 0.0065 * height / 288.0;
    assert!(actual_base.powf(5.255).is_finite());
    assert!(star_pressure_from_altitude_m(height).unwrap().is_infinite());
    assert_eq!(star_nonnegative(-0.0).to_bits(), (-0.0_f64).to_bits());
    assert!(star_nonnegative(f64::NAN).is_nan());
}

#[test]
fn star_transit_ignores_unused_ieee_fields_after_setup() {
    let mut req = request(CALC_MTRANSIT);
    req.lon_deg_east = 0.0;
    req.lat_deg = f64::NAN;
    req.attemp_c = f64::NAN;
    req.horhgt_deg = f64::NAN;
    let p = star_event_parameters(&req).unwrap();
    let out = native_star_event_callbacks(
        -1,
        &req,
        &p,
        &|_, _| Ok([0.0; 6]),
        &|_, _| panic!("unused latitude/altitude"),
        &|_, _| panic!("unused projected frame"),
        &|_| Ok(0.0),
    )
    .unwrap();
    assert_eq!(out.0, 0);
    assert_eq!(out.1[0], req.jd_start_ut);
    assert!(out.1[1..].iter().all(|word| word.to_bits() == 0));
}

#[test]
fn unknown_catalog_precedes_nan_epoch_on_all_star_event_branches() {
    let session = Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
    for event in [CALC_RISE, CALC_SET, CALC_MTRANSIT, CALC_ITRANSIT] {
        let result = session.star_rise_trans(
            "No Such Star xyz",
            f64::NAN,
            event,
            [0.0, f64::NAN, 0.0],
            1013.25,
            15.0,
            0,
        );
        assert!(matches!(
            result.as_ref().err().map(Error::primary_error),
            Some(Error::StarNotFound { .. })
        ));
    }
}

#[test]
fn checked_planetary_callback_admission_stays_before_acquisition() {
    let mut req = request(CALC_RISE);
    req.jd_start_ut = f64::NAN;
    let p = star_event_parameters(&req).unwrap();
    let result = super::super::rise::native_event_callbacks(
        MARS,
        &req,
        &p,
        &|_, _| panic!("checked callback consumed"),
        &|_, _| panic!("checked altitude consumed"),
        &|_, _| panic!("checked projection consumed"),
        &|_| panic!("checked GAST consumed"),
    );
    assert!(matches!(result, Err(Error::NonFiniteSearchValue { .. })));
}

#[test]
fn raw_transit_catalog_error_retains_identity_before_gast() {
    let error = Error::StarNotFound {
        star_id: "toy".into(),
        search_type: Some("name"),
    };
    let result =
        star_calculate_transit(f64::NAN, f64::NAN, false, &|_| Err(error.clone()), &|_| {
            panic!("GAST after catalog failure")
        });
    assert_eq!(result.unwrap_err(), error);
}
