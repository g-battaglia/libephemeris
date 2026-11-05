use super::*;

#[test]
fn acquired_aerosol_selection_is_not_rederived_from_the_scalar_interval() {
    assert_eq!(
        calc_aerosol_extinction_selected(f64::NAN, 0.3, 0.0, AerosolSelection::Range),
        Ok(3.912 / 0.3 - 0.106)
    );
    assert_eq!(
        calc_aerosol_extinction_selected(0.5, 10.0, 0.0, AerosolSelection::Humidity),
        Ok(0.04 + 0.18 * 0.5 + 0.12 * 0.5_f64.powi(2))
    );
    for range in [0.0, -0.0] {
        assert_eq!(
            calc_aerosol_extinction_selected(1e200, range, -1e7, AerosolSelection::Range),
            Err(SchaeferError::DivisionByZero)
        );
    }
    assert_eq!(
        calc_aerosol_extinction_selected(1e200, 10.0, -1e7, AerosolSelection::Humidity),
        Err(SchaeferError::PowerRange)
    );
    assert_eq!(
        calc_aerosol_extinction_selected(1e200, 0.3, -1e7, AerosolSelection::Direct),
        Ok(0.3)
    );
}

#[test]
fn rayleigh_sea_level_and_checked_exponential_domain() {
    assert_eq!(calc_rayleigh_extinction(1013.25, 0.0), Ok(0.1451));
    assert!(calc_rayleigh_extinction(1013.25, 8500.0).unwrap() < 0.1451);
    assert_eq!(
        calc_rayleigh_extinction(0.0, -1e7),
        Err(SchaeferError::MathRange)
    );
    assert!(
        calc_rayleigh_extinction(1.0, f64::NEG_INFINITY)
            .unwrap()
            .is_infinite()
    );
    assert!(
        calc_rayleigh_extinction(0.0, f64::NEG_INFINITY)
            .unwrap()
            .is_nan()
    );
    assert!(calc_rayleigh_extinction(-1.0, 0.0).unwrap() < 0.0);
}

#[test]
fn aerosol_direct_coefficient_skips_humidity_and_altitude_errors() {
    assert_eq!(calc_aerosol_extinction(1e308, 0.3, -1e7), Ok(0.3));
    assert_eq!(calc_aerosol_extinction(f64::NAN, 0.3, f64::NAN), Ok(0.3));
    assert_eq!(
        calc_aerosol_extinction(1e308, 0.0, 0.0),
        Err(SchaeferError::PowerRange)
    );
    assert_eq!(
        calc_aerosol_extinction(1e308, -1.0, -1e7),
        Err(SchaeferError::PowerRange)
    );
    assert_eq!(
        calc_aerosol_extinction(0.5, 1.0, -1e7),
        Err(SchaeferError::MathRange)
    );
    assert_eq!(
        calc_aerosol_extinction(0.5, 1000.0, 0.0).unwrap().to_bits(),
        0
    );
    assert!(
        calc_aerosol_extinction(0.5, 10.0, 1500.0).unwrap()
            < calc_aerosol_extinction(0.5, 10.0, 0.0).unwrap()
    );
    assert_eq!(
        calc_aerosol_extinction(-0.5, 0.0, 0.0),
        calc_aerosol_extinction(-0.5, f64::NAN, 0.0)
    );
}

#[test]
fn extinction_composition_keeps_unused_temperature_and_error_priority() {
    assert_eq!(calc_ozone_extinction(), 0.016);
    assert_eq!(calc_water_vapor_extinction(-1.0), -0.014);
    assert_eq!(
        calc_water_vapor_extinction(-0.0).to_bits(),
        (-0.0_f64).to_bits()
    );
    let ordinary = calc_total_extinction(1013.25, 15.0, 0.5, 0.0, 0.0);
    assert_eq!(
        ordinary,
        calc_total_extinction(1013.25, f64::NAN, 0.5, 0.0, 0.0)
    );
    assert_eq!(
        ordinary,
        calc_total_extinction(1013.25, f64::INFINITY, 0.5, 0.0, 0.0)
    );
    // Rayleigh is acquired before the humidity square can overflow.
    assert_eq!(
        calc_total_extinction(f64::NAN, 15.0, 1e308, 0.0, -1e7),
        Err(SchaeferError::MathRange)
    );
}

#[test]
fn airmass_clamps_and_extinction_signed_zero() {
    for altitude in [-100.0, -0.0, 0.0, f64::NEG_INFINITY] {
        assert_eq!(calc_airmass(altitude), 40.0);
    }
    for altitude in [90.0, 100.0, f64::INFINITY] {
        assert_eq!(calc_airmass(altitude), 1.0);
    }
    assert!(calc_airmass(f64::NAN).is_nan());
    assert!((1.0..=40.0).contains(&calc_airmass(f64::from_bits(1))));
    assert!(calc_airmass(1.0) > calc_airmass(45.0));
    assert_eq!(
        calc_extinction_magnitude(90.0, -0.0).to_bits(),
        (-0.0_f64).to_bits()
    );
    assert_eq!(calc_extinction_magnitude(0.0, 0.25), 10.0);
}

#[test]
fn twilight_preserves_shared_anchors_and_endpoint_discontinuities() {
    for boundary in [-6.0, -12.0] {
        let shallower = calc_twilight_brightness(boundary + 1e-10);
        let at = calc_twilight_brightness(boundary);
        let deeper = calc_twilight_brightness(boundary - 1e-10);
        assert!(shallower >= at && at >= deeper);
        assert!((shallower - deeper) / at < 1e-8);
    }
    assert_eq!(calc_twilight_brightness(-18.0).to_bits(), 0);
    assert!(calc_twilight_brightness(-18.0 + 1e-10) > 0.0);
    assert_eq!(calc_twilight_brightness(0.0), 1e10);
    assert!(calc_twilight_brightness(-1e-10) < 1e10);
    assert_eq!(calc_twilight_brightness(f64::INFINITY), 1e10);
    assert_eq!(calc_twilight_brightness(f64::NEG_INFINITY), 0.0);
    assert!(calc_twilight_brightness(f64::NAN).is_nan());
}

#[test]
fn moon_early_returns_and_infinite_phase_error() {
    assert_eq!(calc_moon_brightness(0.0, f64::INFINITY, f64::NAN), Ok(0.0));
    assert_eq!(calc_moon_brightness(45.0, 180.0, f64::NAN), Ok(0.0));
    assert_eq!(
        calc_moon_brightness(45.0, f64::INFINITY, 90.0),
        Err(SchaeferError::MathDomain)
    );
    assert_eq!(
        calc_moon_brightness(45.0, 0.0, -100.0),
        calc_moon_brightness(45.0, 0.0, 1.0)
    );
    assert_eq!(calc_moon_brightness(90.0, 0.0, 45.0), Ok(3000.0));
    assert!(
        calc_moon_brightness(45.0, 0.0, 20.0).unwrap()
            > calc_moon_brightness(45.0, 0.0, 45.0).unwrap()
    );
    assert!(calc_moon_brightness(45.0, f64::NAN, 90.0).unwrap().is_nan());
}

#[test]
fn zodiacal_domain_and_zero_factor_still_evaluates_latitude() {
    assert_eq!(calc_airglow_brightness(), 145.0);
    assert_eq!(calc_zodiacal_brightness(0.0, f64::NEG_INFINITY), Ok(0.0));
    assert_eq!(
        calc_zodiacal_brightness(0.0, f64::INFINITY),
        Err(SchaeferError::MathDomain)
    );
    assert!(calc_zodiacal_brightness(f64::NAN, 10.0).unwrap().is_nan());
    assert!(calc_zodiacal_brightness(0.0, 270.0).unwrap() < 0.0);
    assert_eq!(
        calc_zodiacal_brightness(30.0, 45.0),
        calc_zodiacal_brightness(-30.0, 45.0)
    );
    assert_eq!(calc_zodiacal_brightness(0.0, 90.0), Ok(150.0));
}

#[test]
fn total_sky_evaluates_zodiacal_at_zero_weight_and_moon_first() {
    assert_eq!(
        calc_total_sky_brightness(1.0, -90.0, 180.0, 90.0, 0.0, f64::INFINITY),
        Err(SchaeferError::MathDomain)
    );
    assert!(
        calc_total_sky_brightness(1.0, -90.0, 180.0, 90.0, f64::NAN, 10.0)
            .unwrap()
            .is_nan()
    );
    let at_12 = calc_total_sky_brightness(-12.0, -90.0, 180.0, 90.0, 0.0, 90.0).unwrap();
    let below_12 = calc_total_sky_brightness(-12.0 - 1e-10, -90.0, 180.0, 90.0, 0.0, 90.0).unwrap();
    assert!((at_12 - below_12).abs() < 1e-6);
    assert_eq!(
        calc_total_sky_brightness(-18.0, -90.0, 180.0, 90.0, 0.0, 90.0),
        Ok(295.0)
    );
}

#[test]
fn surface_brightness_sentinel_overflow_and_infinite_exponents() {
    assert_eq!(brightness_to_mag_arcsec2(-0.0), 22.0);
    assert_eq!(brightness_to_mag_arcsec2(-1.0), 22.0);
    assert!(brightness_to_mag_arcsec2(f64::NAN).is_nan());
    assert_eq!(brightness_to_mag_arcsec2(f64::INFINITY), f64::NEG_INFINITY);
    assert_eq!(
        mag_arcsec2_to_brightness(-1e6),
        Err(SchaeferError::PowerRange)
    );
    assert_eq!(
        mag_arcsec2_to_brightness(f64::NEG_INFINITY),
        Ok(f64::INFINITY)
    );
    assert_eq!(mag_arcsec2_to_brightness(f64::INFINITY), Ok(0.0));
    assert_eq!(mag_arcsec2_to_brightness(1e6).unwrap().to_bits(), 0);
    let brightness = mag_arcsec2_to_brightness(21.0).unwrap();
    assert!((brightness_to_mag_arcsec2(brightness) - 21.0).abs() < 1e-12);
}

#[test]
fn limiting_magnitude_nonpositive_and_unordered_inputs_follow_branches() {
    let baseline = calc_limiting_magnitude(100.0, 0.0, 90.0, 20.0, 1.0).unwrap();
    assert_eq!(baseline, 6.5);
    for snellen in [-1.0, 0.0, f64::NAN] {
        assert_eq!(
            calc_limiting_magnitude(100.0, 0.0, 90.0, 20.0, snellen),
            Ok(baseline)
        );
    }
    assert_eq!(
        calc_limiting_magnitude(f64::NAN, 0.0, 90.0, f64::NAN, 1.0),
        Ok(baseline)
    );
    assert!(
        calc_limiting_magnitude(100.0, 0.0, 90.0, 20.0, f64::from_bits(1))
            .unwrap()
            .is_infinite()
    );
    assert!(
        calc_limiting_magnitude(100.0, 0.0, f64::NAN, 20.0, 1.0)
            .unwrap()
            .is_nan()
    );
    assert!(calc_limiting_magnitude(100.0, 0.0, 90.0, 20.0, 2.0).unwrap() > baseline);
    assert!(calc_limiting_magnitude(100.0, 0.0, 90.0, 60.0, 1.0).unwrap() < baseline);
    assert_eq!(
        calc_limiting_magnitude(100.0, 0.0, 90.0, 20.0, f64::INFINITY),
        Err(SchaeferError::MathDomain)
    );
}

#[test]
fn visibility_early_returns_preserve_intrinsic_word_and_skip_bad_atmosphere() {
    let horizon = is_object_visible(
        -0.0,
        0.0,
        1.0,
        1013.25,
        15.0,
        1e308,
        0.0,
        -1e7,
        25.0,
        1.0,
        45.0,
        f64::INFINITY,
        90.0,
        true,
    )
    .unwrap();
    assert!(!horizon.0);
    assert_eq!(horizon.1.to_bits(), 0);
    assert_eq!(horizon.2.to_bits(), (-0.0_f64).to_bits());
    let day = is_object_visible(
        2.0,
        20.0,
        0.0,
        1013.25,
        15.0,
        1e308,
        0.0,
        -1e7,
        25.0,
        1.0,
        45.0,
        f64::INFINITY,
        90.0,
        true,
    )
    .unwrap();
    assert_eq!(day, (false, -5.0, 2.0));
}

#[test]
fn visibility_moon_toggle_and_intrinsic_comparison() {
    let included = is_object_visible(
        0.0,
        30.0,
        -20.0,
        1013.25,
        15.0,
        0.5,
        0.0,
        0.0,
        25.0,
        1.0,
        45.0,
        f64::INFINITY,
        90.0,
        true,
    );
    assert_eq!(included, Err(SchaeferError::MathDomain));
    let excluded = is_object_visible(
        5.7,
        30.0,
        -20.0,
        1013.25,
        15.0,
        0.5,
        0.0,
        0.0,
        25.0,
        1.0,
        45.0,
        f64::INFINITY,
        90.0,
        false,
    )
    .unwrap();
    assert!(excluded.0);
    assert!(excluded.2 > excluded.1); // Apparent magnitude is informational.
    assert!(
        !is_object_visible(
            -4.0, 30.0, -20.0, 1013.25, 15.0, 0.5, 0.0, 0.0, 25.0, 1.0, -90.0, 180.0, 6.9, false
        )
        .unwrap()
        .0
    );
}

#[test]
fn arcus_bins_are_half_open_and_unordered_input_uses_bright_fallback() {
    assert_eq!(get_arcus_visionis(-5.0), 5.0);
    assert_eq!(get_arcus_visionis(-2.0), 7.0);
    assert_eq!(get_arcus_visionis(0.0), 9.0);
    assert_eq!(get_arcus_visionis(1.5), 11.0);
    assert_eq!(get_arcus_visionis(3.0), 13.0);
    assert_eq!(get_arcus_visionis(4.5), 15.0);
    assert_eq!(get_arcus_visionis(6.0), 17.0);
    assert_eq!(get_arcus_visionis(f64::NAN), 4.0);
    assert_eq!(get_arcus_visionis(f64::NEG_INFINITY), 4.0);
    assert_eq!(get_arcus_visionis(f64::INFINITY), 17.0);
}

#[test]
fn optimal_sun_altitude_does_not_depend_on_morning_selector() {
    for (magnitude, value) in [
        (-2.1, -8.0),
        (-2.0, -9.0),
        (0.0, -10.0),
        (2.0, -11.0),
        (4.0, -12.0),
        (f64::NAN, -12.0),
    ] {
        assert_eq!(get_optimal_sun_altitude(magnitude, true), value);
        assert_eq!(get_optimal_sun_altitude(magnitude, false), value);
    }
}

#[test]
fn threshold_uses_ordered_floor_and_literal_sun_boundary_jumps() {
    assert_eq!(calc_heliacal_altitude_threshold(0.5, f64::NAN, 0.25), 1.0);
    assert_eq!(calc_heliacal_altitude_threshold(0.5, -15.0, f64::NAN), 1.0);
    assert_eq!(calc_heliacal_altitude_threshold(0.5, -15.0, -10.0), 1.0);
    assert!(
        calc_heliacal_altitude_threshold(0.5, -12.0, 0.25)
            > calc_heliacal_altitude_threshold(0.5, -12.0 - 1e-10, 0.25)
    );
    assert!(
        calc_heliacal_altitude_threshold(0.5, -9.0, 0.25)
            > calc_heliacal_altitude_threshold(0.5, -9.0 - 1e-10, 0.25)
    );
}

#[test]
fn condition_dictionary_keeps_intrinsic_margin_and_diagnostic_only_visibility() {
    let condition = get_visibility_conditions(-15.0, 30.0, 1.0, 1013.25, 0.5, 0.0).unwrap();
    assert_eq!(
        condition.magnitude_margin,
        condition.limiting_magnitude - 1.0
    );
    assert_eq!(condition.is_visible, 1.0 <= condition.limiting_magnitude);
    assert!(condition.apparent_magnitude > 1.0);
    let below_horizon = get_visibility_conditions(-20.0, -5.0, -20.0, 1013.25, 0.5, 0.0).unwrap();
    assert!(below_horizon.is_visible); // Dictionary has no horizon filter.
    assert!(
        get_visibility_conditions(-15.0, 30.0, f64::NAN, 1013.25, 0.5, 0.0)
            .unwrap()
            .magnitude_margin
            .is_nan()
    );
}
