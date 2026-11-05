// SPDX-License-Identifier: AGPL-3.0
// Copyright (c) 2025-2026 Giacomo Battaglia
//! Public-API assertion counterparts of the pinned Schaefer coverage tests.
//!
//! Values and properties come from `tests/test_cov100_schaefer.py` at v3.2.1;
//! source line mappings and omitted-argument defaults accompany this leaf.
//! These assertions supplement the separate exact 453-case acceptance bank.

use libephemeris_core::schaefer;

/// Scalar `pytest.approx(expected)` for the finite expectations in this leaf.
/// Its default asymmetric tolerance is max(1e-12, 1e-6*abs(expected)); equality
/// passes first. Nonfinite mismatches fail. No native-derived reference exists.
fn assert_source_approx(actual: f64, expected: f64) {
    if actual == expected {
        return;
    }
    let tolerance = 1e-12_f64.max(1e-6 * expected.abs());
    assert!(
        actual.is_finite() && expected.is_finite() && (actual - expected).abs() <= tolerance,
        "actual={actual:?}, source expectation={expected:?}, tolerance={tolerance:?}"
    );
}

// Source definition L24; assertions L27, L29.
#[test]
fn test_calc_rayleigh_extinction() {
    let sea = schaefer::calc_rayleigh_extinction(1013.25, 0.0).unwrap();
    assert_source_approx(sea, schaefer::RAYLEIGH_COEFF_SEA_LEVEL);
    let high = schaefer::calc_rayleigh_extinction(1013.25, 8500.0).unwrap();
    assert!(high < sea);
}

// Source definition L32; assertion L35; omitted altitude_m=0.
#[test]
fn test_calc_aerosol_extinction_met_range_ge_one() {
    let k = schaefer::calc_aerosol_extinction(0.5, 25.0, 0.0).unwrap();
    assert_source_approx(k, 0.0_f64.max(3.912 / 25.0 - 0.106));
}

// Source definition L38; assertion L41; omitted altitude_m=0.
#[test]
fn test_calc_aerosol_extinction_met_range_ge_one_clamped() {
    let k = schaefer::calc_aerosol_extinction(0.5, 1000.0, 0.0).unwrap();
    assert_eq!(k, 0.0);
}

// Source definition L44; assertion L47.
#[test]
fn test_calc_aerosol_extinction_direct_coeff() {
    let k = schaefer::calc_aerosol_extinction(0.5, 0.3, 0.0).unwrap();
    assert_source_approx(k, 0.3);
}

// Source definition L50; assertion L53.
#[test]
fn test_calc_aerosol_extinction_from_humidity() {
    let k = schaefer::calc_aerosol_extinction(1.0, 0.0, 0.0).unwrap();
    assert_source_approx(k, 0.04 + 0.18 + 0.12);
}

// Source definition L56; assertion L63.
#[test]
fn test_calc_aerosol_extinction_altitude_reduces() {
    let low = schaefer::calc_aerosol_extinction(0.5, 10.0, 0.0).unwrap();
    let high = schaefer::calc_aerosol_extinction(0.5, 10.0, 1500.0).unwrap();
    assert!(high < low);
}

// Source definition L66; assertions L70, L71.
#[test]
fn test_calc_aerosol_extinction_direct_coeff_not_altitude_scaled() {
    let sea = schaefer::calc_aerosol_extinction(0.5, 0.3, 0.0).unwrap();
    let high = schaefer::calc_aerosol_extinction(0.5, 0.3, 1500.0).unwrap();
    assert_eq!(sea, 0.3);
    assert_eq!(high, 0.3);
}

// Source definition L74; assertion L76.
#[test]
fn test_calc_ozone_extinction() {
    assert_eq!(schaefer::calc_ozone_extinction(), schaefer::OZONE_COEFF);
}

// Source definition L79; assertion L81-83.
#[test]
fn test_calc_water_vapor_extinction() {
    assert_source_approx(
        schaefer::calc_water_vapor_extinction(0.5),
        schaefer::WATER_VAPOR_COEFF * 0.5,
    );
}

// Source definition L86; assertion L95; all five omitted defaults materialized.
#[test]
fn test_calc_total_extinction() {
    let k = schaefer::calc_total_extinction(1013.25, 15.0, 0.5, 0.0, 0.0).unwrap();
    let expected = schaefer::calc_rayleigh_extinction(1013.25, 0.0).unwrap()
        + schaefer::calc_aerosol_extinction(0.5, 0.0, 0.0).unwrap()
        + schaefer::calc_ozone_extinction()
        + schaefer::calc_water_vapor_extinction(0.5);
    assert_source_approx(k, expected);
}

// Source definition L103; assertions L105, L106.
#[test]
fn test_calc_airmass_below_horizon() {
    assert_eq!(schaefer::calc_airmass(0.0), 40.0);
    assert_eq!(schaefer::calc_airmass(-5.0), 40.0);
}

// Source definition L109; assertions L111, L112.
#[test]
fn test_calc_airmass_zenith() {
    assert_eq!(schaefer::calc_airmass(90.0), 1.0);
    assert_eq!(schaefer::calc_airmass(95.0), 1.0);
}

// Source definition L115; assertion L118 (both chained comparisons).
#[test]
fn test_calc_airmass_normal() {
    let am = schaefer::calc_airmass(45.0);
    assert!(1.0 < am && am < 40.0);
}

// Source definition L121; assertion L124 (strict lower, inclusive upper).
#[test]
fn test_calc_airmass_near_horizon() {
    let am = schaefer::calc_airmass(0.001);
    assert!(1.0 < am && am <= 40.0);
}

// Source definition L132; assertion L135.
#[test]
fn test_calc_extinction_magnitude() {
    let mag = schaefer::calc_extinction_magnitude(45.0, 0.3);
    assert_source_approx(mag, 0.3 * schaefer::calc_airmass(45.0));
}

// Source definition L143; assertions L145, L146.
#[test]
fn test_calc_twilight_brightness_daytime() {
    assert_eq!(schaefer::calc_twilight_brightness(0.0), 1e10);
    assert_eq!(schaefer::calc_twilight_brightness(5.0), 1e10);
}

// Source definition L149; assertions L151, L152.
#[test]
fn test_calc_twilight_brightness_night() {
    assert_eq!(schaefer::calc_twilight_brightness(-18.0), 0.0);
    assert_eq!(schaefer::calc_twilight_brightness(-25.0), 0.0);
}

// Source definition L155; assertion L158.
#[test]
fn test_calc_twilight_brightness_civil() {
    let b = schaefer::calc_twilight_brightness(-3.0);
    assert_source_approx(b, 1e8 * (-3.0_f64 * 1.5).exp());
}

// Source definition L161; assertion L167.
#[test]
fn test_calc_twilight_brightness_nautical() {
    let b = schaefer::calc_twilight_brightness(-9.0);
    let b_civil_at_6 = 1e8 * (-6.0_f64 * 1.5).exp();
    let b_astro_at_12 = 1e3;
    let frac = (9.0 - 6.0) / 6.0;
    assert_source_approx(b, b_civil_at_6 * (b_astro_at_12 / b_civil_at_6).powf(frac));
}

// Source definition L170; assertion L173.
#[test]
fn test_calc_twilight_brightness_astronomical() {
    let b = schaefer::calc_twilight_brightness(-15.0);
    assert_source_approx(b, 1e3 * (-3.0_f64 * 0.4).exp());
}

// Source definition L181; assertion L183.
#[test]
fn test_calc_moon_brightness_below_horizon() {
    assert_eq!(
        schaefer::calc_moon_brightness(-5.0, 0.0, 45.0).unwrap(),
        0.0
    );
}

// Source definition L186; assertion L188.
#[test]
fn test_calc_moon_brightness_new_moon() {
    assert_eq!(
        schaefer::calc_moon_brightness(45.0, 180.0, 45.0).unwrap(),
        0.0
    );
}

// Source definition L191; assertion L194.
#[test]
fn test_calc_moon_brightness_close_distance() {
    let b = schaefer::calc_moon_brightness(45.0, 0.0, 0.0).unwrap();
    assert!(b > 0.0);
}

// Source definition L197; assertion L200.
#[test]
fn test_calc_moon_brightness_mid_distance() {
    let b = schaefer::calc_moon_brightness(45.0, 0.0, 20.0).unwrap();
    assert!(b > 0.0);
}

// Source definition L203; assertion L209.
#[test]
fn test_calc_moon_brightness_far_distance() {
    let b = schaefer::calc_moon_brightness(45.0, 0.0, 90.0).unwrap();
    let moon_airmass = schaefer::calc_airmass(45.0);
    let alt_factor = (-0.4 * (moon_airmass - 1.0)).exp();
    let expected = 3000.0 * 1.0 * alt_factor * 1.0;
    assert_source_approx(b, expected);
}

// Source definition L217; assertion L219.
#[test]
fn test_calc_airglow_brightness() {
    assert_eq!(
        schaefer::calc_airglow_brightness(),
        schaefer::AIRGLOW_BRIGHTNESS_NL
    );
}

// Source definition L222; assertion L224.
#[test]
fn test_calc_zodiacal_brightness_too_close_to_sun() {
    assert_eq!(schaefer::calc_zodiacal_brightness(0.0, 10.0).unwrap(), 0.0);
}

// Source definition L227; assertion L231.
#[test]
fn test_calc_zodiacal_brightness_peak() {
    let b = schaefer::calc_zodiacal_brightness(0.0, 45.0).unwrap();
    let expected = 150.0 * 0.0_f64.exp() * 45.0_f64.to_radians().sin();
    assert_source_approx(b, expected);
}

// Source definition L234; assertion L238.
#[test]
fn test_calc_zodiacal_brightness_anti_sun() {
    let b = schaefer::calc_zodiacal_brightness(0.0, 120.0).unwrap();
    let expected = 150.0 * 0.0_f64.exp() * 60.0_f64.to_radians().sin();
    assert_source_approx(b, expected);
}

// Source definition L241; assertion L245.
#[test]
fn test_calc_zodiacal_brightness_latitude_factor() {
    let near = schaefer::calc_zodiacal_brightness(0.0, 45.0).unwrap();
    let far = schaefer::calc_zodiacal_brightness(30.0, 45.0).unwrap();
    assert!(far < near);
}

// Source definition L253; assertion L274.
#[test]
fn test_calc_total_sky_brightness_with_zodiacal() {
    let b = schaefer::calc_total_sky_brightness(-15.0, -90.0, 180.0, 90.0, 0.0, 45.0).unwrap();
    let expected = schaefer::calc_twilight_brightness(-15.0)
        + 0.0
        + schaefer::calc_airglow_brightness()
        + 0.5 * schaefer::calc_zodiacal_brightness(0.0, 45.0).unwrap();
    assert_source_approx(b, expected);
}

// Source definition L277; assertion L289; omitted Moon/ecliptic defaults explicit.
#[test]
fn test_calc_total_sky_brightness_without_zodiacal() {
    let b = schaefer::calc_total_sky_brightness(-6.0, -90.0, 180.0, 90.0, 0.0, 45.0).unwrap();
    let expected =
        schaefer::calc_twilight_brightness(-6.0) + 0.0 + schaefer::calc_airglow_brightness() + 0.0;
    assert_source_approx(b, expected);
}

// Source definition L297; assertions L299, L300.
#[test]
fn test_brightness_to_mag_arcsec2_dark() {
    assert_eq!(schaefer::brightness_to_mag_arcsec2(0.0), 22.0);
    assert_eq!(schaefer::brightness_to_mag_arcsec2(-5.0), 22.0);
}

// Source definition L303; assertion L306.
#[test]
fn test_brightness_to_mag_arcsec2_positive() {
    let m = schaefer::brightness_to_mag_arcsec2(250.0);
    assert_source_approx(m, 26.33 - 2.5 * 250.0_f64.log10());
}

// Source definition L309; assertions L312, L313.
#[test]
fn test_mag_arcsec2_to_brightness_roundtrip() {
    let b = schaefer::mag_arcsec2_to_brightness(21.0).unwrap();
    assert_source_approx(b, 10.0_f64.powf((26.33 - 21.0) / 2.5));
    assert_source_approx(schaefer::brightness_to_mag_arcsec2(b), 21.0);
}

// Source definition L321; assertion L330.
#[test]
fn test_calc_limiting_magnitude_dark_young_normal() {
    let lim = schaefer::calc_limiting_magnitude(100.0, 0.0, 90.0, 20.0, 1.0).unwrap();
    assert_source_approx(lim, schaefer::IDEAL_NAKED_EYE_LIMIT_MAG);
}

// Source definition L333; assertion L349.
#[test]
fn test_calc_limiting_magnitude_bright_sky_old_poor_vision() {
    let lim = schaefer::calc_limiting_magnitude(2500.0, 0.3, 30.0, 60.0, 0.5).unwrap();
    let base = schaefer::IDEAL_NAKED_EYE_LIMIT_MAG;
    let sky_factor = 2.5 * (2500.0_f64 / 250.0).log10();
    let age_factor = schaefer::AGE_DEGRADATION_PER_DECADE * ((60.0 - 20.0) / 10.0);
    let snellen_factor = 2.5 * (1.0_f64 / 0.5).log10();
    let airmass = schaefer::calc_airmass(30.0);
    let extinction_loss = 0.3 * (airmass - 1.0);
    let expected = base - sky_factor - age_factor - snellen_factor - extinction_loss;
    assert_source_approx(lim, expected);
}

// Source definition L352; assertion L361.
#[test]
fn test_calc_limiting_magnitude_snellen_zero() {
    let lim = schaefer::calc_limiting_magnitude(100.0, 0.0, 90.0, 20.0, 0.0).unwrap();
    assert_source_approx(lim, schaefer::IDEAL_NAKED_EYE_LIMIT_MAG);
}

// Source definition L369; assertions L376-378; omitted source defaults explicit.
#[test]
fn test_is_object_visible_below_horizon() {
    let (vis, lim, app) = schaefer::is_object_visible(
        0.0, -1.0, -15.0, 1013.25, 15.0, 0.5, 0.0, 0.0, 25.0, 1.0, -90.0, 180.0, 90.0, true,
    )
    .unwrap();
    assert!(!vis);
    assert_eq!(lim, 0.0);
    assert_eq!(app, 0.0);
}

// Source definition L381; assertions L388-390; omitted source defaults explicit.
#[test]
fn test_is_object_visible_sun_up() {
    let (vis, lim, app) = schaefer::is_object_visible(
        0.0, 20.0, 1.0, 1013.25, 15.0, 0.5, 0.0, 0.0, 25.0, 1.0, -90.0, 180.0, 90.0, true,
    )
    .unwrap();
    assert!(!vis);
    assert_eq!(lim, -5.0);
    assert_eq!(app, 0.0);
}

// Source definition L393; assertions L403, L404; include_moon=True retained.
#[test]
fn test_is_object_visible_bright_object_with_moon() {
    let (vis, lim, app) = schaefer::is_object_visible(
        -4.0, 20.0, -15.0, 1013.25, 15.0, 0.5, 0.0, 0.0, 25.0, 1.0, -90.0, 180.0, 40.0, true,
    )
    .unwrap();
    assert!(vis);
    assert!(lim > app);
}

// Source definition L407; type assertions L416-418; include_moon=False retained.
#[test]
fn test_is_object_visible_normal_object_no_moon() {
    // Rust's statically checked tuple components are the source isinstance
    // counterparts. The real public call must still succeed at this input.
    let (_vis, _lim, _app): (bool, f64, f64) = schaefer::is_object_visible(
        0.0, 30.0, -15.0, 1013.25, 15.0, 0.5, 0.0, 0.0, 25.0, 1.0, -90.0, 180.0, 40.0, false,
    )
    .unwrap();
}

// Source definition L421; assertion L429; source default include_moon=True.
#[test]
fn test_is_object_visible_faint_object_elongation_fail() {
    let (vis, _lim, _app) = schaefer::is_object_visible(
        4.0, 30.0, -15.0, 1013.25, 15.0, 0.5, 0.0, 0.0, 25.0, 1.0, -90.0, 180.0, 5.0, true,
    )
    .unwrap();
    assert!(!vis);
}

// Source definition L432; assertion L440; source default include_moon=True.
#[test]
fn test_is_object_visible_bright_object_elongation_fail() {
    let (vis, _lim, _app) = schaefer::is_object_visible(
        -4.0, 30.0, -15.0, 1013.25, 15.0, 0.5, 0.0, 0.0, 25.0, 1.0, -90.0, 180.0, 1.0, true,
    )
    .unwrap();
    assert!(!vis);
}

// Source definition L448; assertion L450.
#[test]
fn test_get_arcus_visionis_table_match() {
    assert_eq!(schaefer::get_arcus_visionis(0.5), 9.0);
}

// Source definition L453; assertion L455.
#[test]
fn test_get_arcus_visionis_very_faint() {
    assert_eq!(schaefer::get_arcus_visionis(6.5), 17.0);
}

// Source definition L458; assertion L460.
#[test]
fn test_get_arcus_visionis_extremely_bright() {
    assert_eq!(schaefer::get_arcus_visionis(-10.0), 4.0);
}

// Source definition L478/assertion L480; parameter pair L471.
#[test]
fn test_get_optimal_sun_altitude_mag_neg_5() {
    assert_eq!(schaefer::get_optimal_sun_altitude(-5.0, true), -8.0);
}

// Source definition L478/assertion L480; parameter pair L472.
#[test]
fn test_get_optimal_sun_altitude_mag_neg_1() {
    assert_eq!(schaefer::get_optimal_sun_altitude(-1.0, true), -9.0);
}

// Source definition L478/assertion L480; parameter pair L473.
#[test]
fn test_get_optimal_sun_altitude_mag_1() {
    assert_eq!(schaefer::get_optimal_sun_altitude(1.0, true), -10.0);
}

// Source definition L478/assertion L480; parameter pair L474.
#[test]
fn test_get_optimal_sun_altitude_mag_3() {
    assert_eq!(schaefer::get_optimal_sun_altitude(3.0, true), -11.0);
}

// Source definition L478/assertion L480; parameter pair L475.
#[test]
fn test_get_optimal_sun_altitude_mag_5() {
    assert_eq!(schaefer::get_optimal_sun_altitude(5.0, true), -12.0);
}

// Source definition L488; assertion L493.
#[test]
fn test_calc_heliacal_altitude_threshold_deep_sun() {
    let th = schaefer::calc_heliacal_altitude_threshold(0.5, -15.0, 0.25);
    let base_av = schaefer::get_arcus_visionis(0.5);
    let expected = 1.0_f64.max(base_av * 1.0 * 0.8);
    assert_source_approx(th, expected);
}

// Source definition L496; assertion L501.
#[test]
fn test_calc_heliacal_altitude_threshold_mid_sun() {
    let th = schaefer::calc_heliacal_altitude_threshold(0.5, -10.0, 0.25);
    let base_av = schaefer::get_arcus_visionis(0.5);
    let expected = 1.0_f64.max(base_av * 1.0 * 1.0);
    assert_source_approx(th, expected);
}

// Source definition L504; assertion L510.
#[test]
fn test_calc_heliacal_altitude_threshold_shallow_sun() {
    let th = schaefer::calc_heliacal_altitude_threshold(0.5, -6.0, 0.25);
    let base_av = schaefer::get_arcus_visionis(0.5);
    let sun_factor = 1.2 + (-6.0 + 9.0) * 0.1;
    let expected = 1.0_f64.max(base_av * 1.0 * sun_factor);
    assert_source_approx(th, expected);
}

// Source definition L513; assertion L516.
#[test]
fn test_calc_heliacal_altitude_threshold_floor() {
    let th = schaefer::calc_heliacal_altitude_threshold(-4.0, -15.0, 0.0);
    assert!(th >= 1.0);
}

// Source definition L524; assertions L546, L547, L552-555.
#[test]
fn test_get_visibility_conditions_keys_and_values() {
    let conditions =
        schaefer::get_visibility_conditions(-15.0, 30.0, 1.0, 1013.25, 0.5, 0.0).unwrap();
    // An exhaustive pattern without '..' checks the actual thirteen-name
    // typed schema, the counterpart of the original dictionary key set.
    let schaefer::VisibilityConditions {
        extinction_total,
        extinction_rayleigh,
        extinction_aerosol,
        airmass,
        extinction_magnitude,
        sky_brightness_nl,
        sky_brightness_mag_arcsec2,
        limiting_magnitude,
        apparent_magnitude,
        magnitude_margin,
        is_visible,
        arcus_visionis,
        optimal_sun_altitude,
    } = conditions;
    let _ = (
        extinction_total,
        extinction_rayleigh,
        extinction_aerosol,
        airmass,
        extinction_magnitude,
        sky_brightness_nl,
        sky_brightness_mag_arcsec2,
        limiting_magnitude,
        apparent_magnitude,
        magnitude_margin,
        is_visible,
        arcus_visionis,
        optimal_sun_altitude,
    );
    let _: bool = is_visible;
    assert_source_approx(magnitude_margin, limiting_magnitude - 1.0);
    assert_eq!(is_visible, 1.0 <= limiting_magnitude);
}
