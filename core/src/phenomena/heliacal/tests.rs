use super::*;

fn standard_model() -> VislimModel {
    let atmo = VislimAtmosphere::standard(0.0, 41.9, 2451545.0).unwrap();
    VislimModel::new(&atmo, HeliacalObserver::naked_eye(), true)
}

#[test]
fn body_predicates_match_release_ids() {
    // Release `is_fixed_star`: body >= FIXSTAR_OFFSET (1000000).
    assert!(!is_fixed_star(FIXSTAR_OFFSET - 1));
    assert!(is_fixed_star(FIXSTAR_OFFSET));
    assert!(is_fixed_star(FIXSTAR_OFFSET + 4));
    for id in [i32::MIN, -1, 0, 1, SUN, MOON, MERCURY, VENUS, 999_999] {
        assert!(!is_fixed_star(id), "id={id}");
    }
    // Release `is_inner_planet` / INNER_PLANETS = {MERCURY, VENUS}.
    assert!((MERCURY, VENUS) == (2, 3));
    assert!(is_inner_planet(MERCURY));
    assert!(is_inner_planet(VENUS));
    for id in [
        0,
        1,
        4,
        5,
        9,
        10,
        15,
        FIXSTAR_OFFSET,
        FIXSTAR_OFFSET + 4,
        -1,
    ] {
        assert!(!is_inner_planet(id), "id={id}");
    }
}

#[test]
fn yallop_bands_and_boundaries() {
    assert_eq!(yallop_visibility_code(1.0), 1.0);
    assert_eq!(yallop_visibility_code(0.217), 1.0);
    assert_eq!(yallop_visibility_code(0.216), 2.0);
    assert_eq!(yallop_visibility_code(0.0), 2.0);
    assert_eq!(yallop_visibility_code(-0.014), 3.0);
    assert_eq!(yallop_visibility_code(-0.160), 4.0);
    assert_eq!(yallop_visibility_code(-0.232), 5.0);
    assert_eq!(yallop_visibility_code(-0.293), 6.0);
    assert_eq!(yallop_visibility_code(-2.0), 6.0);
}

#[test]
fn arc_of_light_known_values() {
    assert_eq!(arc_of_light(0.0, 0.0), 0.0);
    assert!((arc_of_light(90.0, 0.0) - 90.0).abs() < 1e-12);
    assert!((arc_of_light(0.0, 90.0) - 90.0).abs() < 1e-12);
    // cos is even and 360-periodic: DAZ sign and full turns change
    // nothing (up to 1 ulp of libm rounding, so tolerance, not bits).
    assert!((arc_of_light(10.0, 30.0) - arc_of_light(10.0, -30.0)).abs() < 1e-12);
    assert!((arc_of_light(10.0, 30.0) - arc_of_light(10.0, 390.0)).abs() < 1e-12);
    // Small-angle spot: the hypotenuse exceeds either leg.
    let arcl = arc_of_light(5.0, 5.0);
    assert!(arcl > 5.0 && arcl < 7.1, "arcl {arcl}");
}

#[test]
fn arc_of_light_preserves_unknown_geometry() {
    assert!(arc_of_light(f64::NAN, 0.0).is_nan());
    assert!(arc_of_light(0.0, f64::NAN).is_nan());
    assert_eq!(arc_of_light(0.0, 0.0), 0.0);
    assert!((arc_of_light(90.0, 0.0) - 90.0).abs() < 1e-12);
    assert!((arc_of_light(180.0, 0.0) - 180.0).abs() < 1e-12);
}

#[test]
fn humidity_bounds_preserve_nan_and_finite_limits() {
    let aerosol = |humidity| vl_extinction_components(15.0, humidity, 0.0, 0.0, 2451545.0, 0.0).1;
    assert!(aerosol(f64::NAN).is_nan());
    assert_eq!(aerosol(0.0), aerosol(1.0));
    assert_eq!(aerosol(100.0), aerosol(99.5));
}

#[test]
fn vl_month_cont_j2000_and_range() {
    // JD 2451545.0 is 2000-01-01 12:00 UT: January, day fraction 0.5/30.4.
    let m = vl_month_cont(2451545.0);
    assert!((m - 1.0).abs() < 0.05, "month {m}");
    for jd in [2451545.0, 2451600.0, 2451700.0, 2451800.0, 2451900.0] {
        let month = vl_month_cont(jd);
        assert!((1.0..13.0).contains(&month), "jd {jd} month {month}");
    }
}

#[test]
fn vl_extinction_defaults_and_pressure_branch() {
    let (kr, ka, ko, kw) = vl_extinction_components(15.0, 50.0, 0.0, 0.0, 2451545.0, 0.0);
    assert!((kr - 0.1066).abs() < 1e-15, "kr {kr}");
    assert!(ka > 0.05 && ka < 0.2, "ka {ka}");
    // Seasonal ozone at the equator in January sits below the 0.031
    // mean ( `3 + 0.4 * (LT*cos(RA) - cos(3LT))` over 3 arm).
    assert!(ko > 0.02 && ko < VL_OZ, "ko {ko}");
    assert!(kw > 0.0 && kw < 0.05, "kw {kw}");
    // Explicit site pressure encodes the column: altitude is then ignored.
    let (low, _, _, _) = vl_extinction_components(15.0, 50.0, 0.0, 0.0, 2451545.0, 900.0);
    let (high, _, _, _) = vl_extinction_components(15.0, 50.0, 5000.0, 0.0, 2451545.0, 900.0);
    assert_eq!(low.to_bits(), high.to_bits());
    assert!((low - 0.1066 * 900.0 / 1013.25).abs() < 1e-15);
    // Humidity clamps: 0% behaves as 1%, 100% as 99.5%.
    let (_, dry, _, _) = vl_extinction_components(15.0, 0.0, 0.0, 0.0, 2451545.0, 0.0);
    let (_, floor, _, _) = vl_extinction_components(15.0, 1.0, 0.0, 0.0, 2451545.0, 0.0);
    assert_eq!(dry.to_bits(), floor.to_bits());
}

#[test]
fn vl_airmass_horizon_and_zenith() {
    assert_eq!(vl_airmass(0.0), 40.0);
    assert_eq!(vl_airmass(-3.0), 40.0);
    let zenith = vl_airmass(90.0);
    assert!((zenith - 1.0).abs() < 1e-3, "zenith {zenith}");
    // Airmass grows as the object sinks.
    assert!(vl_airmass(10.0) > vl_airmass(30.0));
    assert!(vl_airmass(30.0) > vl_airmass(90.0));
    // Ozone shell guard: nearly-horizontal rays stay finite.
    assert!(vl_airmass_ozone(0.5).is_finite());
}

#[test]
fn hecht_threshold_regime_switch() {
    //  branch spelling, pinned bitwise against the formula.
    let bl = 2000.0_f64;
    let expected = 10_f64.powf(-8.350001) * (1.0 + (10_f64.powf(-5.9) * bl).sqrt()).powf(2.0);
    assert_eq!(vl_hecht_threshold(bl).to_bits(), expected.to_bits());
    let dark = 100.0_f64;
    let expected_dark = 10_f64.powf(-9.8) * (1.0 + (10_f64.powf(-1.9) * dark).sqrt()).powf(2.0);
    assert_eq!(vl_hecht_threshold(dark).to_bits(), expected_dark.to_bits());
    // Brighter skies need more contrast light.
    assert!(vl_hecht_threshold(5000.0) > vl_hecht_threshold(100.0));
}

#[test]
fn eye_pupil_reference_and_cap() {
    assert_eq!(
        vl_eye_pupil(VL_PUPIL_REF_AGE).to_bits(),
        (VL_PUPIL_D0 * (-(36.0_f64 * 36.0) / VL_PUPIL_K).exp()).to_bits()
    );
    // Below age 23 the pupil is capped at the age-23 value.
    assert_eq!(
        vl_eye_pupil(5.0).to_bits(),
        vl_eye_pupil(VL_PUPIL_AGE_MIN).to_bits()
    );
    assert!(vl_eye_pupil(70.0) < vl_eye_pupil(36.0));
}

#[test]
fn observer_optics_flag_gating() {
    let words = [25.0, 1.2, 1.0, 10.0, 50.0, 0.9];
    // Age and Snellen are always read; without the optical-params flag
    // the instrument words fall back to the naked-eye no-ops.
    let gated = parse_observer_optics(&words, 0).unwrap();
    assert_eq!(gated.age_years(), 25.0);
    assert_eq!(gated.snellen_ratio(), 1.2);
    assert!(!gated.binocular());
    assert_eq!(gated.telescope_mag(), 0.0);
    assert_eq!(gated.aperture_mm(), 0.0);
    assert_eq!(gated.transmission(), 1.0);
    // With the flag every word is read (binocular nonzero = true).
    let aided = parse_observer_optics(&words, crate::constants::HELFLAG_OPTICAL_PARAMS).unwrap();
    assert_eq!(aided.age_years(), 25.0);
    assert!(aided.binocular());
    assert_eq!(aided.aperture_mm(), 50.0);
    // Missing words fall back to the naked-eye defaults.
    let short = parse_observer_optics(&[], 0).unwrap();
    assert_eq!(short, HeliacalObserver::naked_eye());
    // Non-finite words are refused fail-closed.
    assert!(parse_observer_optics(&[f64::NAN, 1.0], 0).is_err());
}

#[test]
fn constructors_refuse_non_finite() {
    assert!(HeliacalObserver::new(f64::INFINITY, 1.0, false, 0.0, 0.0, 1.0).is_err());
    assert!(VislimAtmosphere::new(1013.25, 15.0, 50.0, 0.0, 0.0, 41.9, f64::NAN).is_err());
    // Finite extremes are admitted like the release (no range validation).
    assert!(VislimAtmosphere::new(-5.0, -273.0, 250.0, 0.0, -430.0, 95.0, 0.0).is_ok());
}

#[test]
fn seasonal_model_state_depends_on_requested_epoch_not_sample_history() {
    let start = 2_460_310.5;
    for pressure_scaled in [false, true] {
        let atmosphere =
            VislimAtmosphere::new(800.0, 5.0, 70.0, 0.25, 2000.0, 42.0, start).unwrap();
        let observer = HeliacalObserver::naked_eye();
        let mut reused = VislimModel::new(&atmosphere, observer, pressure_scaled);
        for step in [1.0, 0.25, -0.5, 2.0, 4.0, 0.25, 0.25] {
            let epoch = start + step;
            reused.update_season(epoch);
            let current =
                VislimAtmosphere::new(800.0, 5.0, 70.0, 0.25, 2000.0, 42.0, epoch).unwrap();
            let fresh = VislimModel::new(&current, observer, pressure_scaled);
            assert_eq!(reused, fresh);
            assert_eq!(
                reused.limiting_magnitude(-9.0, 30.0, 0.5, 20.0, 60.0, 50.0),
                fresh.limiting_magnitude(-9.0, 30.0, 0.5, 20.0, 60.0, 50.0)
            );
        }
    }
}

#[test]
fn optical_margin_retains_physical_units_and_geometric_exclusions() {
    let model = standard_model();
    let words = [10.0, -4.0, -8.0, 25.0, -90.0, 0.0, 180.0, 0.0];
    let margin = model.visibility_margin(words).unwrap().unwrap();
    assert!(margin > 0.0);
    assert_eq!(
        margin,
        model.limiting_magnitude(-8.0, -90.0, 0.0, 10.0, 25.0, 180.0) + 4.0
    );
    for (index, value) in [(0, -1.0), (2, 1.0), (3, 6.0)] {
        let mut excluded = words;
        excluded[index] = value;
        assert_eq!(model.visibility_margin(excluded).unwrap(), None);
    }
    let mut undefined = words;
    undefined[1] = f64::NAN;
    assert!(matches!(
        model.visibility_margin(undefined),
        Err(Error::NonFiniteArgument { .. })
    ));
}

#[test]
fn visibility_gates_and_margin() {
    let model = standard_model();
    // Below the horizon: never visible.
    assert!(!model.is_visible(-1.0, -4.0, -10.0, 30.0, -90.0, 0.0, 180.0, 0.0));
    // Sun above the horizon: never a heliacal visibility.
    assert!(!model.is_visible(20.0, -4.0, 5.0, 30.0, -90.0, 0.0, 180.0, 0.0));
    // Inside the 7-degree glare gate: never visible, however bright.
    assert!(!model.is_visible(20.0, -10.0, -10.0, 6.9, -90.0, 0.0, 180.0, 0.0));
    // Venus-bright at dusk clears the threshold at neutral margin...
    assert!(model.is_visible(10.0, -4.0, -8.0, 25.0, -90.0, 0.0, 180.0, 0.0));
    // ...but a large detection margin can still reject it.
    assert!(!model.is_visible(10.0, -4.0, -8.0, 25.0, -90.0, 0.0, 180.0, 99.0));
    // A faint star at the same geometry does not clear it.
    assert!(!model.is_visible(10.0, 6.0, -8.0, 25.0, -90.0, 0.0, 180.0, 0.0));
}

#[test]
fn vision_kind_defaults_and_overrides() {
    let model = standard_model();
    // Deep twilight reads scotopic without overrides...
    let dark = vislim_vision_kind(&model, -30.0, -90.0, 0.0, 45.0, 120.0, 180.0, 0);
    assert_eq!(dark, VisionKind::Scotopic);
    // ...the unforced verdict always follows the 1500 nL switch...
    let sky = model.sky_brightness_nl(-12.0, -90.0, 0.0, 30.0, 60.0, 180.0);
    let expected = if sky < SCOTOPIC_SKY_NL {
        VisionKind::Scotopic
    } else {
        VisionKind::Photopic
    };
    assert_eq!(
        vislim_vision_kind(&model, -12.0, -90.0, 0.0, 30.0, 60.0, 180.0, 0),
        expected
    );
    // ...unless photopic is forced, and scotopic wins when forced.
    assert_eq!(
        vislim_vision_kind(
            &model,
            -12.0,
            -90.0,
            0.0,
            30.0,
            60.0,
            180.0,
            HELFLAG_VISLIM_PHOTOPIC
        ),
        VisionKind::Photopic
    );
    assert_eq!(
        vislim_vision_kind(
            &model,
            -12.0,
            -90.0,
            0.0,
            30.0,
            60.0,
            180.0,
            HELFLAG_VISLIM_SCOTOPIC
        ),
        VisionKind::Scotopic
    );
    assert_eq!(VisionKind::Photopic.as_flag(), 0);
    assert_eq!(VisionKind::Scotopic.as_flag(), 1);
}

#[test]
fn body_event_and_side_gates() {
    assert_eq!(parse_heliacal_event(1), Some(HeliacalEvent::Rising));
    assert_eq!(parse_heliacal_event(4), Some(HeliacalEvent::MorningLast));
    assert_eq!(parse_heliacal_event(0), None);
    assert_eq!(parse_heliacal_event(5), None);
    // The Sun is never valid.
    assert!(validate_heliacal_body(SUN, HeliacalEvent::Rising).is_err());
    // The Moon has no rising/setting but owns the crescent events.
    assert!(validate_heliacal_body(MOON, HeliacalEvent::Rising).is_err());
    assert!(validate_heliacal_body(MOON, HeliacalEvent::EveningFirst).is_ok());
    // Stars rise and set but never open/close an apparition.
    assert!(validate_heliacal_body(FIXSTAR_OFFSET + 4, HeliacalEvent::Setting).is_ok());
    assert!(validate_heliacal_body(FIXSTAR_OFFSET + 4, HeliacalEvent::MorningLast).is_err());
    // Outer planets have no evening-first.
    assert!(validate_heliacal_body(4, HeliacalEvent::Rising).is_ok());
    assert!(validate_heliacal_body(4, HeliacalEvent::EveningFirst).is_err());
    assert!(validate_heliacal_body(VENUS, HeliacalEvent::MorningLast).is_ok());
    // Unknown ids refuse as unknown bodies, never as valid planets.
    assert!(matches!(
        validate_heliacal_body(99, HeliacalEvent::Rising),
        Err(Error::UnknownBody { body_id: 99 })
    ));
    // Morning events need the body west (negative diff), evening east.
    assert!(event_side_ok(-15.0, true));
    assert!(!event_side_ok(15.0, true));
    assert!(event_side_ok(15.0, false));
    assert!(!event_side_ok(-15.0, false));
    // Longitude-difference normalization matches positive-divisor remainder.
    assert!((sun_body_lon_diff_deg(10.0, 350.0) - 20.0).abs() < 1e-12);
    assert!((sun_body_lon_diff_deg(350.0, 10.0) + 20.0).abs() < 1e-12);
}

#[test]
fn rising_streak_needs_conjunction_gap() {
    // A mid-apparition flicker (visible before, short dip, no gap) is ignored.
    let flicker = [
        HeliacalDaySample {
            jd_day: 1.0,
            visible: true,
            jd_best: 1.1,
            elongation_deg: 40.0,
        },
        HeliacalDaySample {
            jd_day: 2.0,
            visible: false,
            jd_best: 0.0,
            elongation_deg: 41.0,
        },
        HeliacalDaySample {
            jd_day: 3.0,
            visible: false,
            jd_best: 0.0,
            elongation_deg: 42.0,
        },
        HeliacalDaySample {
            jd_day: 4.0,
            visible: true,
            jd_best: 4.1,
            elongation_deg: 43.0,
        },
    ];
    assert_eq!(find_heliacal_rising_day(&flicker, false, false), None);
    // A true emergence through the gap triggers on the first visible day.
    let mut emergence = vec![
        HeliacalDaySample {
            jd_day: 1.0,
            visible: false,
            jd_best: 0.0,
            elongation_deg: 4.0,
        },
        HeliacalDaySample {
            jd_day: 2.0,
            visible: false,
            jd_best: 0.0,
            elongation_deg: 5.0,
        },
        HeliacalDaySample {
            jd_day: 3.0,
            visible: false,
            jd_best: 0.0,
            elongation_deg: 6.0,
        },
        HeliacalDaySample {
            jd_day: 4.0,
            visible: false,
            jd_best: 0.0,
            elongation_deg: 7.0,
        },
        HeliacalDaySample {
            jd_day: 5.0,
            visible: false,
            jd_best: 0.0,
            elongation_deg: 8.0,
        },
    ];
    emergence.push(HeliacalDaySample {
        jd_day: 6.0,
        visible: true,
        jd_best: 6.1,
        elongation_deg: 12.0,
    });
    assert_eq!(find_heliacal_rising_day(&emergence, false, false), Some(5));
    // Setting mirrors the shape: last visible before the 5-day streak.
    let mut dusk = emergence.clone();
    dusk.push(HeliacalDaySample {
        jd_day: 7.0,
        visible: false,
        jd_best: 0.0,
        elongation_deg: 13.0,
    });
    dusk.push(HeliacalDaySample {
        jd_day: 8.0,
        visible: false,
        jd_best: 0.0,
        elongation_deg: 14.0,
    });
    dusk.push(HeliacalDaySample {
        jd_day: 9.0,
        visible: false,
        jd_best: 0.0,
        elongation_deg: 15.0,
    });
    dusk.push(HeliacalDaySample {
        jd_day: 10.0,
        visible: false,
        jd_best: 0.0,
        elongation_deg: 16.0,
    });
    dusk.push(HeliacalDaySample {
        jd_day: 11.0,
        visible: false,
        jd_best: 0.0,
        elongation_deg: 17.0,
    });
    assert_eq!(find_heliacal_setting_day(&dusk), Some(5));
}

#[test]
fn refine_converges_on_step_transition() {
    // Morning polarity: visibility turns ON with time; the refined
    // time brackets the step from above.
    let t0 = 2460000.25_f64;
    let morning = refine_heliacal_time(|jd| jd >= t0, t0 + 0.05, true).unwrap();
    assert!((morning - t0).abs() < 2e-6, "morning {morning}");
    // Evening polarity mirrors it.
    let evening = refine_heliacal_time(|jd| jd <= t0, t0 - 0.05, false).unwrap();
    assert!((evening - t0).abs() < 2e-6, "evening {evening}");
    // A non-finite anchor is refused, never bisected.
    assert!(refine_heliacal_time(|_| true, f64::NAN, true).is_err());
}
