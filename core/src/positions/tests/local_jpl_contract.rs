// SPDX-License-Identifier: AGPL-3.0
//! Focused tests for the parent position-routing module.

#[test]
fn direct_type21_retains_split_tdb_and_reads_final_fourth_sample() {
    let samples = std::cell::RefCell::new(Vec::new());
    let fraction = 0.125;
    let vectors = [
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [4.0, 0.0, 0.0],
    ];
    let value = super::direct_type21_relative_at(fraction, [0.0; 3], &|sample| {
        let index = samples.borrow().len();
        samples.borrow_mut().push(sample);
        Ok(vectors[index])
    })
    .unwrap();
    assert_eq!(value, vectors[3]);
    assert_eq!(
        *samples.borrow(),
        [
            fraction,
            fraction - 1.0 / 173.1446326846693,
            fraction - 2.0 / 173.1446326846693,
            fraction - 3.0 / 173.1446326846693
        ]
    );
}
#[test]
fn direct_minor_spk_reads_three_scalar_tt_samples_against_fixed_observer() {
    let epochs = std::cell::RefCell::new(Vec::new());
    let tt = 2451545.25;
    let observer = [0.25, -0.5, 0.75];
    let vectors = [[1.25, -0.5, 0.75], [2.25, -0.5, 0.75], [3.25, -0.5, 0.75]];
    let got = super::minor_spk_relative_at(tt, observer, &|epoch| {
        let index = epochs.borrow().len();
        epochs.borrow_mut().push(epoch);
        Ok(vectors[index])
    })
    .unwrap();
    assert_eq!(got, [3.0, 0.0, 0.0]);
    assert_eq!(
        *epochs.borrow(),
        [tt, tt - 1.0 / 173.144632674240, tt - 2.0 / 173.144632674240]
    );
}

#[test]
fn direct_minor_spk_read_failure_stops_the_original_sample_order() {
    let reads = std::cell::Cell::new(0);
    let result = super::minor_spk_relative_at(2451545.0, [0.0; 3], &|_| {
        reads.set(reads.get() + 1);
        Err(super::Error::SpkBodyNotFound { target: 2000001 })
    });
    assert!(matches!(
        result,
        Err(super::Error::SpkBodyNotFound { target: 2000001 })
    ));
    assert_eq!(reads.get(), 1);
}

#[test]
fn public_horizons_entries_use_owned_cache_and_native_only_for_refused_frames() {
    use crate::horizons::{HorizonsTimeType, StateVector};
    let ctx = crate::context::EphemerisContext::new()
        .with_mode(crate::context::ProviderMode::HorizonsRemote)
        .with_horizons_config(Default::default())
        .unwrap()
        .with_network_policy(Some(crate::context::NetworkPolicy::Sealed));
    let cfg = ctx.delta_t_config(crate::time::DeltaTModel::EspenakMeeus, false);
    let ut = super::FRAME_J2000;
    let tt = ut + super::deltat_with_config(ut, &cfg).unwrap();
    let client = ctx.horizons_client().unwrap();
    for (command, x) in [("499", 2.0), ("399", 1.0)] {
        client
            .get_or_fetch_with(command, tt, "@0", HorizonsTimeType::Tdb, || {
                StateVector::new(x, 0.0, 0.0, 0.0, 0.0, 0.0)
            })
            .unwrap();
    }
    let request = super::ServingRequest {
        ctx: &ctx,
        delta_t_model: crate::time::DeltaTModel::EspenakMeeus,
        iers_enabled: false,
        observer: None,
        tier_policy: super::TierPolicy::AttachOrder,
        skip_closed: false,
    };
    let (_, flags) = super::calc_horizons_ut(
        request,
        ut,
        crate::constants::MARS,
        crate::constants::FLG_TRUEPOS,
    )
    .unwrap();
    assert_eq!(
        flags & crate::constants::FLG_TRUEPOS,
        crate::constants::FLG_TRUEPOS
    );
    assert!(
        super::calc_horizons(
            request,
            tt,
            crate::constants::MARS,
            crate::constants::FLG_TRUEPOS
        )
        .is_ok()
    );
    assert!(matches!(
        super::calc_horizons_ut(
            request,
            ut,
            crate::constants::MARS,
            crate::constants::FLG_ICRS
        ),
        Err(super::Error::SpkUnregisteredBody { .. })
    ));
    // A different UT has no declared remote input: no source interpolation
    // or hidden native substitution can turn sealed transport into values.
    assert!(matches!(
        super::calc_horizons_ut(
            request,
            ut + 1.0,
            crate::constants::MARS,
            crate::constants::FLG_TRUEPOS
        ),
        Err(super::Error::NetworkSealed { .. })
    ));
}

#[test]
fn horizons_fallback_keeps_configuration_and_range_errors_distinct() {
    assert!(super::horizons_fallback_error(
        &super::Error::ZeroVectorNormalize
    ));
    assert!(super::horizons_fallback_error(
        &super::Error::SpkFileNotFound {
            path: "fixture".into()
        }
    ));
    assert!(!super::horizons_fallback_error(
        &super::Error::NetworkSealed {
            purpose: "fixture".into()
        }
    ));
    assert!(!super::horizons_fallback_error(
        &super::Error::ConfigurationMissing {
            missing_config: "fixture",
            suggestion: "fixture"
        }
    ));
    assert!(!super::horizons_fallback_error(
        &super::Error::EphemerisRange {
            requested_jd: 1.0,
            start_jd: 2.0,
            end_jd: 3.0,
            body: Some(1)
        }
    ));
}
#[test]
fn analytical_lunar_points_require_configured_main_image() {
    let registry = super::SpkRegistry::new();
    let leg = super::JplLeg::new(&registry);
    for body in [
        crate::constants::MEAN_NODE,
        crate::constants::MEAN_APOG,
        crate::constants::INTP_APOG,
        crate::constants::INTP_PERG,
    ] {
        // The public factory acquires the main DE before local formulas.
        // Mean models additionally validate TT in mean_lunar_apse; the
        // fitted interpolated formulas do not perform that range check.
        assert_eq!(
            super::local_jpl_lunar_position(
                &leg,
                2451545.0,
                None,
                body,
                crate::constants::FLG_NONUT,
                &|_| Ok(0.0),
            ),
            Err(super::Error::SpkUnregisteredBody {
                body: crate::constants::SUN
            })
        );
    }
}

#[test]
fn mean_and_fitted_apsides_use_main_sun_image_with_distinct_date_domains() {
    let mut registry = SpkRegistry::new();
    let target_offset = crate::spk::RECORD_BYTES + 3 * 8 + 2 * 40 + 16;
    let mut main = crate::spk::moon_test_kernel_bytes();
    main[target_offset..target_offset + 4].copy_from_slice(&10_i32.to_le_bytes());
    registry.register_image("main", main).unwrap();
    registry.register_body(SUN, "main", 10).unwrap();
    let mut moon = crate::spk::moon_test_kernel_bytes();
    moon[target_offset..target_offset + 4].copy_from_slice(&301_i32.to_le_bytes());
    for index in 0..3 {
        let descriptor = crate::spk::RECORD_BYTES + 3 * 8 + index * 40;
        moon[descriptor + 8..descriptor + 16].copy_from_slice(&100.0_f64.to_le_bytes());
    }
    registry.register_image("separate-moon", moon).unwrap();
    registry.register_body(MOON, "separate-moon", 301).unwrap();
    let leg = JplLeg::new(&registry);
    assert_eq!(
        local_jpl_main_window(&leg).unwrap(),
        (FRAME_J2000, FRAME_J2000 + 20_000_000.0 / 86400.0)
    );
    for body in [MEAN_NODE, MEAN_APOG] {
        // The independent Moon image cannot narrow a mean formula's main
        // file domain. No coefficient record is usable in these images.
        assert!(
            local_jpl_lunar_position(
                &leg,
                FRAME_J2000 + 10.0,
                None,
                body,
                FLG_NONUT,
                &|_| panic!("tropical model consumes no ayanamsha")
            )
            .is_ok()
        );
        let error =
            local_jpl_lunar_position(&leg, FRAME_J2000 - 1.0, None, body, FLG_NONUT, &|_| Ok(0.0))
                .unwrap_err();
        assert!(matches!(error, Error::EphemerisRange { requested_jd,
                body: Some(id), .. } if requested_jd == FRAME_J2000 - 1.0 && id == body));
    }
    for body in [INTP_APOG, INTP_PERG] {
        let (output, flags) =
            local_jpl_lunar_position(&leg, FRAME_J2000 - 1.0, None, body, FLG_NONUT, &|_| {
                panic!("fitted tropical point consumes no ayanamsha")
            })
            .unwrap();
        assert!(output.into_iter().all(f64::is_finite));
        assert_eq!(flags, FLG_NONUT);
    }
}

#[test]
fn standard_native_ut_range_precedes_time_and_keeps_original_epoch() {
    let mut registry = SpkRegistry::new();
    let mut main = crate::spk::moon_test_kernel_bytes();
    let target_offset = crate::spk::RECORD_BYTES + 3 * 8 + 2 * 40 + 16;
    main[target_offset..target_offset + 4].copy_from_slice(&4_i32.to_le_bytes());
    registry.register_image("main", main).unwrap();
    registry.register_body(MARS, "main", 4).unwrap();
    let mut context = crate::context::EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
    *context.spk_registry_mut() = registry;
    let request = ServingRequest {
        ctx: &context,
        delta_t_model: crate::time::DeltaTModel::EspenakMeeus,
        iers_enabled: false,
        observer: None,
        tier_policy: TierPolicy::AttachOrder,
        skip_closed: false,
    };
    for epoch in [-0.0_f64, 0.0, 1.0e308] {
        let error = calc_native_jpl_entry(request, epoch, MARS, 0, true).unwrap_err();
        let Error::EphemerisRange {
            requested_jd, body, ..
        } = error.primary_error()
        else {
            panic!("scalar public range must precede independent TT construction")
        };
        assert_eq!(requested_jd.to_bits(), epoch.to_bits());
        assert_eq!(*body, Some(MARS));
    }
}

#[test]
fn osculating_lunar_points_acquire_earth_before_main_envelope() {
    let mut registry = super::SpkRegistry::new();
    let mut image = crate::spk::moon_test_kernel_bytes();
    let target_offset = crate::spk::RECORD_BYTES + 3 * 8 + 2 * 40 + 16;
    image[target_offset..target_offset + 4].copy_from_slice(&301_i32.to_le_bytes());
    registry.register_image("main-moon", image).unwrap();
    registry.register_body(MOON, "main-moon", 301).unwrap();
    let mut main = crate::spk::moon_test_kernel_bytes();
    main[target_offset..target_offset + 4].copy_from_slice(&10_i32.to_le_bytes());
    registry.register_image("main-sun", main).unwrap();
    registry.register_body(SUN, "main-sun", 10).unwrap();
    let leg = super::JplLeg::new(&registry);
    let outside_main_tt = FRAME_J2000 - 1.0;
    for body in [TRUE_NODE, OSCU_APOG] {
        assert_eq!(
            super::local_jpl_lunar_position(
                &leg,
                outside_main_tt,
                None,
                body,
                FLG_NONUT,
                &|_| panic!("failed vector acquisition has no sidereal tail"),
            ),
            Err(Error::SpkUnregisteredBody { body: EARTH }),
        );
    }
    assert!(matches!(
        super::local_jpl_lunar_position(
            &leg,
            outside_main_tt,
            None,
            MEAN_NODE,
            FLG_NONUT,
            &|_| Ok(0.0)
        ),
        Err(Error::EphemerisRange {
            body: Some(MEAN_NODE),
            ..
        })
    ));
}

use super::*;

#[test]
fn native_classic_entry_clock_retains_original_ut_and_fresh_sample_inverse() {
    let original_ut: f64 = 2451545.125;
    let central_tt = crate::time::default_timescale_tt_from_ut1(original_ut).unwrap();
    assert_eq!(
        native_model_ut_at(central_tt, Some(original_ut))
            .unwrap()
            .to_bits(),
        original_ut.to_bits()
    );
    for sample_tt in [central_tt - 1.0, central_tt + 1.0] {
        let fresh_ut = crate::time::default_timescale_ut1_from_tt(sample_tt).unwrap();
        assert_eq!(
            native_model_ut_at(sample_tt, Some(original_ut))
                .unwrap()
                .to_bits(),
            fresh_ut.to_bits()
        );
        assert_ne!(fresh_ut.to_bits(), original_ut.to_bits());
    }
    assert_eq!(
        native_model_ut_at(central_tt, None).unwrap().to_bits(),
        crate::time::default_timescale_ut1_from_tt(central_tt)
            .unwrap()
            .to_bits()
    );
}

#[test]
fn native_public_ut_preserves_typed_model_tt_ranges_and_promotes_raw_ranges() {
    let context =
        crate::context::EphemerisContext::new().with_mode(crate::context::ProviderMode::NativeJpl);
    let request = ServingRequest {
        ctx: &context,
        delta_t_model: crate::time::DeltaTModel::EspenakMeeus,
        iers_enabled: false,
        observer: None,
        tier_policy: TierPolicy::AttachOrder,
        skip_closed: false,
    };
    let central_tt: f64 = 2451545.125;
    for original_ut in [0.0_f64, -0.0_f64, 2451545.125] {
        let result: Result<()> =
            with_native_jpl_entry_outputs(request, Some(original_ut), |_, _| {
                Err(Error::EphemerisRange {
                    requested_jd: central_tt + 1.0,
                    start_jd: 2287184.5,
                    end_jd: 2688976.5,
                    body: Some(NIBIRU),
                })
            });
        let Error::EphemerisRange {
            requested_jd,
            start_jd,
            end_jd,
            body,
        } = result.unwrap_err()
        else {
            panic!("public range shape lost")
        };
        assert_eq!(requested_jd.to_bits(), (central_tt + 1.0_f64).to_bits());
        assert_eq!(
            (start_jd, end_jd, body),
            (2287184.5, 2688976.5, Some(NIBIRU))
        );
    }
    let tt_error = Error::EphemerisRange {
        requested_jd: central_tt,
        start_jd: 2287184.5,
        end_jd: 2688976.5,
        body: Some(NIBIRU),
    };
    let registry = SpkRegistry::new();
    let leg = JplLeg::new(&registry);
    assert_eq!(
        local_jpl_public_range_error(&leg, -0.0, MEAN_NODE, tt_error.clone()),
        tt_error
    );
    let raw = Error::SpkOutOfRange {
        target: 10,
        et: -1.0,
        start: 0.0,
        end: 1.0,
    };
    let promoted = local_jpl_public_range_error(&leg, -0.0, NIBIRU, raw);
    let Error::EphemerisRange {
        requested_jd, body, ..
    } = promoted
    else {
        panic!("raw coverage failure was not promoted")
    };
    assert_eq!(requested_jd.to_bits(), (-0.0_f64).to_bits());
    assert_eq!(body, Some(NIBIRU));
    let unknown = Error::UnknownBody { body_id: NIBIRU };
    assert_eq!(
        local_jpl_public_range_error(&leg, -0.0, NIBIRU, unknown.clone()),
        unknown
    );
}

#[test]
fn native_classic_models_dispatch_tt_ut_before_unknown_target_resolution() {
    let context = crate::context::EphemerisContext::new()
        .with_mode(crate::context::ProviderMode::NativeJpl)
        .with_network_policy(Some(crate::context::NetworkPolicy::Sealed));
    let request = ServingRequest {
        ctx: &context,
        delta_t_model: crate::time::DeltaTModel::EspenakMeeus,
        iers_enabled: false,
        observer: None,
        tier_policy: TierPolicy::AttachOrder,
        skip_closed: false,
    };
    for body in CUPIDO..=NIBIRU {
        for is_ut in [false, true] {
            let result = with_local_jpl_calc_outputs(request, |engine, _| {
                if is_ut {
                    calc_ut(engine, FRAME_J2000, body, FLG_TRUEPOS)
                } else {
                    calc(engine, FRAME_J2000, body, FLG_TRUEPOS)
                }
            });
            assert_eq!(result, Err(Error::SpkUnregisteredBody { body: SUN }));
            let heliocentric = with_local_jpl_calc_outputs(request, |engine, _| {
                calc(
                    engine,
                    FRAME_J2000,
                    body,
                    FLG_HELCTR | FLG_J2000 | FLG_TRUEPOS,
                )
            });
            if body == NIBIRU {
                assert_eq!(heliocentric, Err(Error::UnknownBody { body_id: NIBIRU }));
            } else {
                assert!(heliocentric.unwrap().0.into_iter().all(f64::is_finite));
            }
        }
    }
}

#[test]
fn native_classic_range_uses_main_sun_image_envelope() {
    let mut registry = SpkRegistry::new();
    let mut image = crate::spk::moon_test_kernel_bytes();
    let target_offset = crate::spk::RECORD_BYTES + 3 * 8 + 16;
    image[target_offset..target_offset + 4].copy_from_slice(&10_i32.to_le_bytes());
    registry.register_image("classic-support", image).unwrap();
    registry.register_body(SUN, "classic-support", 10).unwrap();
    let kernel = registry.kernel("classic-support").unwrap();
    let spans: Vec<_> = kernel
        .index()
        .all_segments()
        .map(|s| s.descriptor)
        .collect();
    let start = spans
        .iter()
        .map(|s| s.start_et)
        .fold(f64::INFINITY, f64::min);
    let end = spans
        .iter()
        .map(|s| s.end_et)
        .fold(f64::NEG_INFINITY, f64::max);
    for body in CUPIDO..=NIBIRU {
        let error = local_jpl_public_range_error(
            &JplLeg::new(&registry),
            -1.5,
            body,
            Error::SpkOutOfRange {
                target: 10,
                et: -1.0,
                start: 0.0,
                end: 1.0,
            },
        );
        assert_eq!(
            error,
            Error::EphemerisRange {
                requested_jd: -1.5,
                start_jd: FRAME_J2000 + start / 86400.0,
                end_jd: FRAME_J2000 + end / 86400.0,
                body: Some(body),
            }
        );
    }
}

#[test]
fn retarded_split_preserves_fractional_light_time_and_validates_domain() {
    for (whole, fraction, delay, expected) in [
        (2451545.0_f64, 0.75_f64, 0.25_f64, 0.5_f64),
        (2415020.0, 0.125, 0.25, -0.125),
        (-1000.0, -0.125, 0.25, -0.375),
    ] {
        assert_eq!(
            retarded_tdb_split(whole, fraction, delay).unwrap(),
            (whole, expected)
        );
    }
    // A sub-JD-ulp delay is retained in the fractional component.
    let delay = 2.0_f64.powi(-40);
    let (whole, fraction) = retarded_tdb_split(2451545.0, 0.5, delay).unwrap();
    assert_eq!(whole, 2451545.0);
    assert_eq!(0.5 - fraction, delay);
    assert_eq!(2451545.5 - delay, 2451545.5);
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(retarded_tdb_split(value, 0.0, 0.0).is_err());
        assert!(retarded_tdb_split(2451545.0, value, 0.0).is_err());
        assert!(retarded_tdb_split(2451545.0, 0.0, value).is_err());
    }
    assert!(retarded_tdb_split(2451545.0, f64::MAX, -f64::MAX).is_err());
}

#[test]
fn local_spk_public_range_retains_central_tt_and_requested_body() {
    let registry = SpkRegistry::new();
    let leg = JplLeg::new(&registry);
    let tt = FRAME_J2000 + 0.25;
    assert_eq!(
        local_jpl_public_range_error(
            &leg,
            tt,
            MOON,
            Error::SpkOutOfRange {
                target: 399,
                et: -100.0,
                start: 0.0,
                end: 86400.0,
            }
        ),
        Error::EphemerisRange {
            requested_jd: tt,
            start_jd: FRAME_J2000,
            end_jd: FRAME_J2000 + 1.0,
            body: Some(MOON),
        }
    );
    let error = Error::SpkUnregisteredBody { body: MARS };
    assert_eq!(
        local_jpl_public_range_error(&leg, tt, MARS, error.clone()),
        error
    );
}

#[test]
fn local_spk_sidereal_uses_scalar_ut1_speed_samples() {
    use std::cell::RefCell;
    let samples = RefCell::new(Vec::new());
    let tt = FRAME_J2000 + 0.125;
    let ut = crate::time::default_timescale_ut1_from_tt(tt).unwrap();
    let h = 1.0 / 86400.0;
    let mut output = [120.0, 2.0, 1.0, 0.5, 0.25, 0.125];
    local_jpl_sidereal_word(&mut output, tt, None, FLG_SIDEREAL | FLG_SPEED, &|epoch| {
        samples.borrow_mut().push(epoch);
        Ok(20.0)
    })
    .unwrap();
    assert_eq!(*samples.borrow(), [ut, ut - h, ut + h]);
    assert_eq!(output, [100.0, 2.0, 1.0, 0.5, 0.25, 0.125]);
}

#[test]
fn native_sidereal_tail_retains_public_ut_and_reacquires_fresh_tt_samples() {
    use std::cell::RefCell;
    for original_ut in [-0.0_f64, FRAME_J2000 + 0.125] {
        let central_tt = crate::time::default_timescale_tt_from_ut1(original_ut).unwrap();
        for tt in [central_tt, central_tt + 1.0] {
            let expected_ut = if tt.to_bits() == central_tt.to_bits() {
                original_ut
            } else {
                crate::time::default_timescale_ut1_from_tt(tt).unwrap()
            };
            let samples = RefCell::new(Vec::new());
            let mut output = [120.0, 2.0, 1.0, 0.5, 0.25, 0.125];
            local_jpl_sidereal_word(
                &mut output,
                tt,
                Some(original_ut),
                FLG_SIDEREAL | FLG_SPEED,
                &|ut| {
                    samples.borrow_mut().push(ut.to_bits());
                    Ok(20.0)
                },
            )
            .unwrap();
            let h = 1.0 / 86400.0;
            assert_eq!(
                *samples.borrow(),
                [
                    expected_ut.to_bits(),
                    (expected_ut - h).to_bits(),
                    (expected_ut + h).to_bits()
                ]
            );
            assert_eq!(output, [100.0, 2.0, 1.0, 0.5, 0.25, 0.125]);
        }
    }
}

#[test]
fn local_spk_sidereal_skips_equatorial_and_origin_anchors() {
    for (mut output, flags) in [
        ([0.0; 6], FLG_SIDEREAL | FLG_SPEED),
        ([1.0; 6], FLG_SIDEREAL | FLG_EQUATORIAL | FLG_SPEED),
        ([1.0; 6], FLG_SPEED),
    ] {
        let original = output;
        local_jpl_sidereal_word(&mut output, FRAME_J2000, None, flags, &|_| {
            panic!("inapplicable ayanamsha channel read")
        })
        .unwrap();
        assert_eq!(output, original);
    }
}

#[test]
fn native_planet_range_rate_differentiates_the_converged_light_path() {
    use std::cell::Cell;
    let observer = [2.0, -0.5, 0.3];
    let observer_velocity = [0.003, -0.001, 0.002];
    let target = [7.0, 1.0, -0.2];
    let target_velocity = [0.01, 0.004, -0.003];
    let clock = 1.0000000003;
    let last_target = Cell::new(target);
    let state = |body, _, _| {
        Ok(if body == EARTH {
            (observer, observer_velocity)
        } else {
            (target, target_velocity)
        })
    };
    let retarded = |_, _, _, delay| {
        let p = std::array::from_fn(|i| target[i] - target_velocity[i] * delay);
        last_target.set(p);
        Ok((p, target_velocity))
    };
    for frame in [FLG_J2000, FLG_EQUATORIAL | FLG_J2000, FLG_NONUT] {
        let place = vector_place_from_state_routes_with_policy(
            &state,
            Some(&retarded),
            FRAME_J2000,
            0.0,
            MARS,
            frame | FLG_NOABERR | FLG_NOGDEFL,
            None,
            false,
            Some(clock),
            None,
            None,
        )
        .unwrap();
        let physical = place.range.unwrap();
        assert_eq!(place.coordinates[2], physical.distance);
        let ray: [f64; 3] = std::array::from_fn(|i| last_target.get()[i] - observer[i]);
        let unit = ray.map(|p| p / physical.distance);
        let numerator: f64 = (0..3)
            .map(|i| unit[i] * (target_velocity[i] - observer_velocity[i]))
            .sum();
        let target_projection: f64 = (0..3).map(|i| unit[i] * target_velocity[i]).sum();
        let exact = clock * numerator / (1.0 + target_projection / FRAME_LIGHT_SPEED_AU_DAY);
        let operation_bound = 64.0 * f64::EPSILON * vec3_dist(target_velocity);
        assert!((physical.rate - exact).abs() <= operation_bound);
    }
}

#[test]
fn registered_spk_retardation_has_its_own_observation_epoch_channel() {
    use std::cell::RefCell;
    let central = RefCell::new(Vec::new());
    let retarded = RefCell::new(Vec::new());
    let state = |id, w, f| {
        central.borrow_mut().push((id, w, f));
        Ok((
            if id == EARTH || id == SUN {
                [0.0; 3]
            } else {
                [0.0024, 0.0, 0.0]
            },
            [0.0; 3],
        ))
    };
    let emit = |id, w, f, delay| {
        retarded.borrow_mut().push((id, w, f, delay));
        Ok(([0.0024, 0.0, 0.0], [0.0; 3]))
    };
    let flags = FLG_J2000 | FLG_NOABERR | FLG_NOGDEFL;
    let _ =
        vector_place_from_state_routes(&state, Some(&emit), FRAME_J2000, 0.5, MOON, flags, None)
            .unwrap();
    assert_eq!(
        *central.borrow(),
        [(EARTH, FRAME_J2000, 0.5), (MOON, FRAME_J2000, 0.5)]
    );
    assert_eq!(
        *retarded.borrow(),
        [(MOON, FRAME_J2000, 0.5, 0.0024 / FRAME_LIGHT_SPEED_AU_DAY)]
    );
    for center_flags in [FLG_HELCTR, FLG_BARYCTR] {
        central.borrow_mut().clear();
        retarded.borrow_mut().clear();
        vector_place_from_state_routes(
            &state,
            Some(&emit),
            FRAME_J2000,
            0.5,
            MOON,
            flags | center_flags,
            None,
        )
        .unwrap();
        assert_eq!(retarded.borrow().len(), 3);
        assert!(
            retarded.borrow().iter().all(|&(id, w, f, delay)| id == MOON
                && w == FRAME_J2000
                && f == 0.5
                && delay > 0.0)
        );
        assert!(
            central
                .borrow()
                .iter()
                .all(|&(_, w, f)| w == FRAME_J2000 && f == 0.5)
        );
    }
    central.borrow_mut().clear();
    retarded.borrow_mut().clear();
    let _ = vector_place_from_state_routes(
        &state,
        Some(&emit),
        FRAME_J2000,
        0.5,
        MOON,
        flags | FLG_TRUEPOS,
        None,
    )
    .unwrap();
    assert!(
        retarded.borrow().is_empty(),
        "TRUEPOS never requests a retarded state"
    );
}

#[test]
fn named_sealed_factory_refuses_non_leb_modes_before_callback() {
    for mode in [ProviderMode::NativeJpl, ProviderMode::HorizonsRemote] {
        let ctx = EphemerisContext::new().with_mode(mode);
        let result: Result<()> =
            with_sealed_calc_engine(&ctx, DeltaTModel::EspenakMeeus, false, |_| {
                panic!("the named sealed factory must not construct a non-LEB callback")
            });
        assert!(matches!(
            result,
            Err(Error::ConfigurationMissing {
                missing_config: "non_leb_leg_wiring",
                ..
            })
        ));
    }
}

#[test]
fn vector_self_observation_never_reads_outside_derivative_epochs() {
    let never_read = |_: i32, _: f64, _: f64| -> Result<([f64; 3], [f64; 3])> {
        panic!("an exact-origin vector must not request a stored state")
    };
    for (body, flags) in [
        (EARTH, FLG_SWIEPH | FLG_SPEED),
        (EARTH, FLG_SWIEPH | FLG_SPEED | FLG_TOPOCTR | FLG_SIDEREAL),
    ] {
        let result =
            vector_position_from_states(&never_read, 2396758.5, 0.0, body, flags, None).unwrap();
        assert_eq!(result.map(f64::to_bits), [0; 6]);
    }
}

#[test]
fn vector_sun_heliocentric_reads_central_state_but_not_derivative_epochs() {
    let calls = std::cell::RefCell::new(Vec::new());
    let states = |id, whole, fraction| {
        calls.borrow_mut().push((id, whole, fraction));
        Ok(([0.1, 0.2, 0.3], [0.0; 3]))
    };
    let flags = FLG_SWIEPH | FLG_SPEED | FLG_HELCTR | FLG_J2000;
    let result = vector_position_from_states(&states, 2396758.5, 0.0, SUN, flags, None).unwrap();
    assert_eq!(result.map(f64::to_bits), [0; 6]);
    assert!(!calls.borrow().is_empty());
    assert!(
        calls
            .borrow()
            .iter()
            .all(|call| *call == (SUN, 2396758.5, 0.0))
    );
    let unavailable = |_: i32, _: f64, _: f64| -> Result<([f64; 3], [f64; 3])> {
        Err(Error::SpkOutOfRange {
            target: 10,
            et: -1.0,
            start: 0.0,
            end: 1.0,
        })
    };
    assert!(matches!(
        vector_position_from_states(&unavailable, 2396758.5, 0.0, SUN, flags, None),
        Err(Error::SpkOutOfRange { target: 10, .. })
    ));
}

#[test]
fn vector_speed_sampling_preserves_actual_fraction_and_split_epoch() {
    for (body, h) in [(MOON, 6.0 / 86400.0), (MARS, 1.0 / 86400.0)] {
        let calls = std::cell::RefCell::new(Vec::new());
        let states = |id, whole, fraction| {
            calls.borrow_mut().push((id, whole, fraction));
            Ok(([1.0, 0.2, 0.1], [0.0; 3]))
        };
        vector_position_from_states(
            &states,
            FRAME_J2000,
            0.2,
            body,
            FLG_SWIEPH | FLG_SPEED | FLG_TRUEPOS | FLG_BARYCTR | FLG_NONUT,
            None,
        )
        .unwrap();
        let calls = calls.borrow();
        assert_eq!(calls[0], (body, FRAME_J2000, 0.2));
        assert!(calls.iter().all(|&(id, whole, fraction)| {
            id == body && whole == FRAME_J2000 && fraction.is_finite()
        }));
        assert!(calls.iter().any(|&(_, _, fraction)| fraction == 0.2 - h));
        assert!(calls.iter().any(|&(_, _, fraction)| fraction == 0.2 + h));
        assert!(
            calls
                .iter()
                .any(|&(_, _, fraction)| (fraction - 0.2).abs() >= 1.0 / 64.0)
        );
    }
}

#[test]
fn local_spk_state_provider_preserves_mode_guard_and_missing_target() {
    let registry = SpkRegistry::new();
    assert!(matches!(
        LocalSpkStateProvider::new(&registry, ProviderMode::SealedLeb),
        Err(Error::ConfigurationMissing {
            missing_config: "jpl_source_access",
            ..
        })
    ));
    let provider = LocalSpkStateProvider::new(&registry, ProviderMode::NativeJpl).unwrap();
    assert_eq!(provider.provider_name(), "native-jpl/local-spk");
    assert!(!provider.has_body(MARS));
    assert!(matches!(
        provider.eval_raw(MARS, FRAME_J2000),
        Err(Error::SpkUnregisteredBody { body: MARS })
    ));
    let (result, events) = crate::trace::capture_events(|| {
        provider.eclipse_source_observed_vector(
            MARS,
            TwoPartEpoch {
                hi: FRAME_J2000,
                lo: 0.0,
            },
            [0.0; 3],
            [0.0; 3],
            None,
        )
    });
    assert!(matches!(
        result,
        Err(Error::SpkUnregisteredBody { body: MARS })
    ));
    assert!(
        events.is_empty(),
        "an unavailable SPK target was not consumed"
    );
    assert_eq!(registry.body_count(), 0);
}

#[test]
fn local_spk_construction_is_lazy_but_core_dispatch_requires_its_target() {
    let ctx = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
    let req = ServingRequest {
        ctx: &ctx,
        delta_t_model: DeltaTModel::EspenakMeeus,
        iers_enabled: false,
        observer: None,
        tier_policy: TierPolicy::AttachOrder,
        skip_closed: false,
    };
    with_shared_calc_outputs(req, |engine, cartesian| {
        // Callback construction is not a calculation. A core calculation
        // still fails with its actual missing target, never empty success.
        assert!(matches!(
            (engine.compute_tt)(FRAME_J2000, MARS, FLG_SWIEPH),
            Err(Error::SpkUnregisteredBody { body: MARS })
        ));
        for body in [
            MEAN_NODE, TRUE_NODE, MEAN_APOG, OSCU_APOG, INTP_APOG, INTP_PERG,
        ] {
            for center in [FLG_HELCTR, FLG_BARYCTR] {
                let flags = center | FLG_SPEED | FLG_SIDEREAL;
                assert_eq!(
                    (engine.compute_tt)(FRAME_J2000, body, flags)?,
                    ([0.0; 6], flags)
                );
            }
        }
        assert!((cartesian)(FRAME_J2000, MARS, FLG_SWIEPH)?.is_none());
        let (nutation, _) = calc(engine, FRAME_J2000, ECL_NUT, FLG_SWIEPH)?;
        assert!(nutation.iter().all(|component| component.is_finite()));
        Ok(())
    })
    .unwrap();
}

#[test]
fn native_empty_main_image_retains_public_sealed_acquisition_priority() {
    use crate::context::NetworkPolicy;

    for (tier, filename) in [
        (Tier::Base, "de440s.bsp"),
        (Tier::Medium, "de440.bsp"),
        (Tier::Extended, "de441.bsp"),
    ] {
        let ctx = EphemerisContext::new()
            .with_mode(ProviderMode::NativeJpl)
            .with_tier(tier)
            .with_network_policy(Some(NetworkPolicy::Sealed));
        let req = ServingRequest {
            ctx: &ctx,
            delta_t_model: DeltaTModel::EspenakMeeus,
            iers_enabled: false,
            observer: None,
            tier_policy: TierPolicy::AttachOrder,
            skip_closed: false,
        };
        for body in [SUN, MOON, MARS, EARTH] {
            for ut_entry in [false, true] {
                assert!(matches!(
                    calc_native_jpl_entry(req, FRAME_J2000, body, 0, ut_entry),
                    Err(Error::NetworkSealed { purpose })
                        if purpose == format!("JPL {filename} ephemeris download")
                ));
            }
        }
        assert!(matches!(
            calc_native_jpl_entry(req, FRAME_J2000, MARS, FLG_BARYCTR | FLG_MOSEPH, true),
            Err(Error::UnsupportedFlagCombination {
                entry: "calc_ut",
                ..
            })
        ));
        let error = calc_native_jpl_entry(req, f64::INFINITY, SUN, 0, true).unwrap_err();
        assert!(matches!(
            error.primary_error(),
            Error::NonFiniteArgument { .. }
        ));
        assert!(calc_native_jpl_entry(req, FRAME_J2000, ECL_NUT, 0, true).is_ok());
        assert!(matches!(
            calc_native_jpl_entry(req, FRAME_J2000, 999, 0, true),
            Err(Error::UnknownBody { body_id: 999 })
        ));
    }
}

#[test]
fn native_ut_central_time_retains_source_pair_without_scalar_reconstruction() {
    let ut = 2_460_400.25;
    let expected = crate::smh2016::timescale_parts_from_ut1(ut).unwrap();
    let tt = expected.0 + expected.1;
    let actual = native_position_central_time(tt, Some(ut)).unwrap();
    assert_eq!(actual.0.to_bits(), ut.to_bits());
    assert_eq!(actual.1.to_bits(), expected.1.to_bits());
    assert_ne!(actual.1.to_bits(), (tt - ut).to_bits());
    for epoch in [tt, tt + 0.25] {
        let fresh = native_position_central_time(epoch, None).unwrap();
        assert_eq!(fresh, (epoch.floor(), epoch - epoch.floor()));
    }
    let side = tt + 0.25;
    assert_eq!(
        native_position_central_time(side, Some(ut)).unwrap(),
        (side.floor(), side - side.floor())
    );
}

#[test]
fn low_level_vector_callback_retains_central_ut_pair_and_split_tt_sides() {
    let ut = 2_460_400.25;
    let expected = crate::smh2016::timescale_parts_from_ut1(ut).unwrap();
    let tt = expected.0 + expected.1;
    for (body, h) in [(MOON, 6.0_f64 / 86400.0), (MARS, 1.0_f64 / 86400.0)] {
        let samples = std::cell::RefCell::new(Vec::new());
        let place = |w: f64, f: f64| {
            let pair = if w.to_bits() == tt.to_bits() && f == 0.0 {
                native_position_central_time(tt, Some(ut))?
            } else {
                (w, f)
            };
            samples
                .borrow_mut()
                .push((pair.0.to_bits(), pair.1.to_bits()));
            Ok([30.0, 2.0, 1.0])
        };
        let result = vector_position_from_places(&place, tt, 0.0, body, FLG_SPEED).unwrap();
        assert_eq!(samples.borrow()[0], (ut.to_bits(), expected.1.to_bits()));
        assert!(
            samples
                .borrow()
                .iter()
                .skip(1)
                .all(|&(whole, _)| whole == tt.to_bits())
        );
        assert!(samples.borrow().contains(&(tt.to_bits(), (-h).to_bits())));
        assert!(samples.borrow().contains(&(tt.to_bits(), h.to_bits())));
        assert_eq!(&result[3..], &[0.0; 3]);
        samples.borrow_mut().clear();
        let no_speed = vector_position_from_places(&place, tt, 0.0, body, 0).unwrap();
        assert_eq!(
            *samples.borrow(),
            vec![(ut.to_bits(), expected.1.to_bits())]
        );
        assert_eq!(no_speed, [30.0, 2.0, 1.0, 0.0, 0.0, 0.0]);
    }
}

#[test]
fn local_spk_missing_registration_is_not_a_mode_wiring_refusal() {
    let registry = SpkRegistry::new();
    let leg = JplLeg::new(&registry);
    assert!(matches!(
        local_jpl_position(&leg, FRAME_J2000, FRAME_J2000, None, MARS, FLG_SWIEPH, None),
        Err(Error::SpkUnregisteredBody { body: MARS })
    ));
}

#[test]
fn local_spk_nonfinite_epoch_precedes_provider_access() {
    let registry = SpkRegistry::new();
    let leg = JplLeg::new(&registry);
    assert!(matches!(
        local_jpl_position(&leg, f64::NAN, f64::NAN, None, MARS, FLG_SWIEPH, None),
        Err(Error::NonFiniteArgument { .. })
    ));
}

#[test]
fn local_spk_equatorial_sidereal_does_not_require_ayanamsha_projection() {
    let registry = SpkRegistry::new();
    let leg = JplLeg::new(&registry);
    assert!(matches!(
        local_jpl_position(
            &leg,
            FRAME_J2000,
            FRAME_J2000,
            None,
            MARS,
            FLG_SWIEPH | FLG_EQUATORIAL | FLG_SIDEREAL,
            None
        ),
        Err(Error::SpkUnregisteredBody { body: MARS })
    ));
}

#[test]
fn local_spk_entry_never_unseals_leb_mode() {
    assert!(matches!(
        forbid_jpl_in_sealed_mode(ProviderMode::SealedLeb),
        Err(Error::ConfigurationMissing {
            missing_config: "jpl_source_access",
            ..
        })
    ));
    assert!(forbid_jpl_in_sealed_mode(ProviderMode::NativeJpl).is_ok());
}

#[test]
fn symbolic_center_state_requires_registered_de_support() {
    let registry = SpkRegistry::new();
    let leg = JplLeg::new(&registry);
    // Heliocentric symbolic models require their registered Sun/Earth support.
    assert_eq!(
        local_jpl_symbolic_center_state(
            &leg,
            FRAME_J2000,
            WHITE_MOON,
            FLG_HELCTR | FLG_SPEED,
            &|_| panic!("White Moon uses no Waldemath anchor")
        ),
        Err(Error::SpkUnregisteredBody { body: SUN })
    );
    // The internal symbolic helper consumes a previously acquired Time;
    // its IEEE epoch does not bypass the required Sun support read.
    assert_eq!(
        local_jpl_symbolic_center_state(&leg, f64::NAN, WHITE_MOON, FLG_HELCTR, &|_| panic!(
            "White Moon uses no Waldemath anchor"
        )),
        Err(Error::SpkUnregisteredBody { body: SUN })
    );
    // The outer public route still rejects a caller's nonfinite epoch.
    assert!(matches!(
        local_jpl_position(&leg, f64::NAN, f64::NAN, None, WHITE_MOON, FLG_HELCTR, None),
        Err(Error::NonFiniteArgument { .. })
    ));
}

#[test]
fn predicted_center_state_requires_registered_de_support() {
    let registry = SpkRegistry::new();
    let leg = JplLeg::new(&registry);
    // Pure heliocentric element propagation needs no kernel. Recentring to the
    // SSB or geocenter requires the registered Sun/Earth states.
    let helio = local_jpl_predicted_center_state(&leg, FRAME_J2000, VULCAN, FLG_HELCTR | FLG_SPEED)
        .unwrap();
    assert!((0.0..360.0).contains(&helio[0]));
    assert!(helio[1].is_finite() && helio[2] > 0.0);
    assert!(helio[3..].iter().all(|rate| rate.is_finite()));
    // Light time is real on this arm: the geometric place differs.
    let geometric =
        local_jpl_predicted_center_state(&leg, FRAME_J2000, VULCAN, FLG_HELCTR | FLG_TRUEPOS)
            .unwrap();
    assert_ne!(helio[0], geometric[0]);
    assert_eq!(
        local_jpl_predicted_center_state(&leg, FRAME_J2000, VULCAN, FLG_BARYCTR),
        Err(Error::SpkUnregisteredBody { body: SUN })
    );
    assert_eq!(
        local_jpl_predicted_center_state(&leg, FRAME_J2000, PROSERPINA, FLG_SPEED),
        Err(Error::SpkUnregisteredBody { body: SUN })
    );
    assert!(matches!(
        local_jpl_predicted_center_state(&leg, f64::NAN, VULCAN, FLG_HELCTR),
        Err(Error::NonFiniteArgument { .. })
    ));
}

#[test]
fn predicted_body_predicate_covers_vulcan_family_only() {
    use crate::horizons::is_predicted_hypothetical_body;
    // Weston Vulcan, Proserpina and the five classical predictions take
    // the heliocentric predicted arm; Nibiru stays model-less and the
    // Uranian/Transpluto/symbolic bodies keep their own arms.
    for body in [
        VULCAN,
        PROSERPINA,
        HARRINGTON,
        NEPTUNE_LEVERRIER,
        NEPTUNE_ADAMS,
        PLUTO_LOWELL,
        PLUTO_PICKERING,
    ] {
        assert!(is_predicted_hypothetical_body(body), "body {body}");
    }
    for body in [CUPIDO, ISIS, NIBIRU, WHITE_MOON, WALDEMATH] {
        assert!(!is_predicted_hypothetical_body(body), "body {body}");
    }
}

#[test]
fn keplerian_helio_needs_no_kernel_but_geocentric_needs_earth() {
    let registry = SpkRegistry::new();
    let leg = JplLeg::new(&registry);
    // Curated heliocentric Keplerian propagation needs no kernel. This route
    // is geometric; without SPEED its reported rates are zero.
    let helio = keplerian_position_at(&leg, FRAME_J2000, PHOLUS, FLG_HELCTR | FLG_SPEED).unwrap();
    assert!((0.0..360.0).contains(&helio.0));
    assert!(helio.1.is_finite() && helio.2 > 0.0);
    assert_eq!(
        helio,
        keplerian_position_at(
            &leg,
            FRAME_J2000,
            PHOLUS,
            FLG_HELCTR | FLG_SPEED | FLG_TRUEPOS
        )
        .unwrap()
    );
    let no_speed = local_jpl_keplerian_fallback(&leg, FRAME_J2000, PHOLUS, FLG_HELCTR).unwrap();
    assert_eq!(&no_speed[3..], &[0.0; 3]);
    assert_eq!(&no_speed[..3], &[helio.0, helio.1, helio.2]);
    // The shared observer boundary retains the registered-SPK wrapper's
    // complete raw model/stencil words; heliocentric acquisition never
    // touches either an injected observer or an empty kernel registry.
    for flags in [FLG_HELCTR, FLG_HELCTR | FLG_SPEED | FLG_J2000] {
        let injected = keplerian_fallback_with_earth(FRAME_J2000, PHOLUS, flags, &|_| {
            panic!("heliocentric samples need no Earth support")
        })
        .unwrap();
        let registered = local_jpl_keplerian_fallback(&leg, FRAME_J2000, PHOLUS, flags).unwrap();
        assert_eq!(injected.map(f64::to_bits), registered.map(f64::to_bits));
    }
    // The geocentric composition resolves Earth then Sun aliases before
    // reading the Sun observer state and the retarded Earth states.
    assert_eq!(
        keplerian_position_at(&leg, FRAME_J2000, PHOLUS, FLG_SPEED),
        Err(Error::SpkUnregisteredBody { body: EARTH })
    );
}

#[test]
fn numbered_asteroid_candidate_range_excludes_offset_and_star_sentinels() {
    for body in [CERES, CHIRON, AST_OFFSET + 90_123, FIXSTAR_OFFSET - 1] {
        assert!(local_jpl_minor_body_candidate(body));
    }
    for body in [AST_OFFSET, FIXSTAR_OFFSET, FIXSTAR_OFFSET + 1, MARS] {
        assert!(!local_jpl_minor_body_candidate(body));
    }
}

#[test]
fn numbered_asteroid_public_tt_ut_preserve_sealed_sbdb_refusal() {
    let context = crate::context::EphemerisContext::new()
        .with_mode(crate::context::ProviderMode::NativeJpl)
        .with_network_policy(Some(crate::context::NetworkPolicy::Sealed))
        .with_strict_precision(Some(false))
        .with_auto_spk_download(Some(false));
    let request = ServingRequest {
        ctx: &context,
        delta_t_model: crate::time::DeltaTModel::EspenakMeeus,
        iers_enabled: false,
        observer: None,
        tier_policy: TierPolicy::AttachOrder,
        skip_closed: false,
    };
    for body in [AST_OFFSET + 90_123, FIXSTAR_OFFSET - 1] {
        for ut in [false, true] {
            let error = with_local_jpl_calc_outputs(request, |engine, _| {
                if ut {
                    calc_ut(engine, FRAME_J2000, body, FLG_SPEED)
                } else {
                    calc(engine, FRAME_J2000, body, FLG_SPEED)
                }
            })
            .unwrap_err();
            assert_eq!(
                error,
                Error::NetworkSealed {
                    purpose: format!("SBDB orbital elements for asteroid {}", body - AST_OFFSET),
                }
            );
        }
    }
}

#[test]
fn numbered_asteroid_sbdb_failure_precedes_nutation_and_earth() {
    use std::cell::RefCell;
    let body = AST_OFFSET + 90_123;
    let seen = RefCell::new(Vec::new());
    let sealed = |number| {
        seen.borrow_mut().push(number);
        Err(Error::NetworkSealed {
            purpose: "injected SBDB refusal".into(),
        })
    };
    for epoch in [FRAME_J2000, f64::MAX] {
        assert!(matches!(
            numbered_asteroid_position_with_earth(epoch, body, 0, &sealed, &|_| {
                panic!("element refusal must precede Earth acquisition")
            }),
            Err(Error::NetworkSealed { .. })
        ));
    }
    assert_eq!(*seen.borrow(), [90_123, 90_123]);
    for failure in [false, true] {
        let unavailable = |_: i64| {
            if failure {
                Err(Error::UnknownBody { body_id: -99 })
            } else {
                Ok(None)
            }
        };
        assert_eq!(
            numbered_asteroid_position_with_earth(f64::MAX, body, 0, &unavailable, &|_| panic!(
                "an unavailable number must not read Earth"
            ),),
            Err(Error::UnknownBody { body_id: body }),
        );
    }
}

#[test]
fn numbered_asteroid_local_model_and_cached_elements_need_no_network() {
    let special = numbered_asteroid_position_with_earth(
        FRAME_J2000,
        AST_OFFSET + 2060,
        FLG_HELCTR | FLG_NONUT,
        &|_| panic!("special local number must precede SBDB"),
        &|_| panic!("heliocentric local number needs no Earth"),
    )
    .unwrap();
    let local = keplerian_position_with_earth(FRAME_J2000, CHIRON, FLG_HELCTR | FLG_NONUT, &|_| {
        panic!("heliocentric curated number needs no Earth")
    })
    .unwrap();
    assert_eq!(
        [special.0, special.1, special.2].map(f64::to_bits),
        [local.0, local.1, local.2].map(f64::to_bits),
    );
    let body = AST_OFFSET + 90_123;
    let registry = SpkRegistry::new();
    let cache = std::cell::RefCell::new(crate::minor_bodies::SbdbElementsCache::new());
    cache.borrow_mut().insert(
        90_123,
        crate::minor_bodies::minor_body_elements(CERES).unwrap(),
    );
    let elements = |number| {
        crate::minor_bodies::fetch_sbdb_elements(
            number,
            30.0,
            crate::context::EffectiveNetworkPolicy::Sealed,
            &mut cache.borrow_mut(),
        )
    };
    let output = local_jpl_numbered_asteroid_fallback(
        &JplLeg::new(&registry),
        FRAME_J2000,
        body,
        FLG_HELCTR | FLG_SPEED,
        &elements,
    )
    .unwrap();
    assert!(output.iter().all(|value| value.is_finite()));
    assert!(output[2] > 0.0);
    assert_eq!(cache.borrow().len(), 1);
    assert_eq!(
        numbered_asteroid_position_with_earth(FRAME_J2000, body, 0, &elements, &|_| {
            Err(Error::SpkUnregisteredBody { body: EARTH })
        }),
        Err(Error::SpkUnregisteredBody { body: EARTH }),
    );
}

#[test]
fn numbered_asteroid_registered_target_precedes_sealed_sbdb() {
    let body = AST_OFFSET + 90_123;
    // TT at J2000 maps to slightly negative TDB seconds. A day inside the
    // advertised kernel domain tests provider priority without crossing it.
    let epoch = FRAME_J2000 + 1.0;
    let mut context = crate::context::EphemerisContext::new()
        .with_mode(crate::context::ProviderMode::NativeJpl)
        .with_network_policy(Some(crate::context::NetworkPolicy::Sealed))
        .with_auto_spk_download(Some(false));
    context
        .spk_registry_mut()
        .register_image("numbered-test", crate::spk::moon_test_kernel_bytes())
        .unwrap();
    context
        .spk_registry_mut()
        .register_body(body, "numbered-test", 3)
        .unwrap();
    assert!(local_jpl_minor_body_servable(
        &JplLeg::new(context.spk_registry()),
        body,
        epoch,
        FLG_BARYCTR | FLG_TRUEPOS,
    ));
    let request = ServingRequest {
        ctx: &context,
        delta_t_model: crate::time::DeltaTModel::EspenakMeeus,
        iers_enabled: false,
        observer: None,
        tier_policy: TierPolicy::AttachOrder,
        skip_closed: false,
    };
    let result = with_local_jpl_calc_outputs(request, |engine, _| {
        calc(engine, epoch, body, FLG_BARYCTR | FLG_TRUEPOS)
    });
    // This summary-only image intentionally has no numerical directory.
    // Its provider failure proves the registered path was attempted; a
    // fallback would instead report the sealed SBDB acquisition refusal.
    assert!(
        matches!(result, Err(Error::SpkCorrupt { .. })),
        "{result:?}"
    );
}

#[test]
fn geometric_type21_admission_requires_the_speed_stencil_without_light_time() {
    let body = PHOLUS;
    let mut registry = SpkRegistry::new();
    registry
        .register_image("minor.bsp", crate::spk::geometric_minor_test_image())
        .unwrap();
    registry.register_body(body, "minor.bsp", 20005145).unwrap();
    let leg = JplLeg::new(&registry);
    let geometric = FLG_HELCTR | FLG_TRUEPOS;
    let first = FRAME_J2000 + 0.005;
    assert!(registered_spk_position_probe(&registry, first, body, geometric | FLG_SPEED).unwrap());
    assert!(matches!(
        registered_spk_heliocentric_probe(&registry, first, body),
        Err(Error::EphemerisRange { .. })
    ));
    assert!(local_jpl_minor_body_servable(
        &leg,
        body,
        first,
        geometric | FLG_SPEED
    ));
    assert!(!local_jpl_minor_body_servable(
        &leg,
        body,
        first,
        FLG_HELCTR | FLG_SPEED
    ));
    let near_start = FRAME_J2000 + 0.001;
    assert!(local_jpl_minor_body_servable(
        &leg, body, near_start, geometric
    ));
    assert!(!local_jpl_minor_body_servable(
        &leg,
        body,
        near_start,
        geometric | FLG_SPEED
    ));
    assert!(local_jpl_minor_body_servable(
        &leg,
        body,
        FRAME_J2000 + 0.995,
        geometric | FLG_SPEED
    ));
    assert!(!local_jpl_minor_body_servable(
        &leg,
        body,
        FRAME_J2000 + 1.001,
        geometric
    ));
}

#[test]
fn numbered_asteroid_support_aliases_precede_state_reads() {
    let mut registry = SpkRegistry::new();
    assert_eq!(
        earth_helio_true_of_date(&JplLeg::new(&registry), FRAME_J2000),
        Err(Error::SpkUnregisteredBody { body: EARTH }),
    );
    let mut image = crate::spk::moon_test_kernel_bytes();
    // The third hand-built descriptor is the SSB Earth link. Give it the
    // geocenter NAIF identity that the public Earth alias requires.
    let target_offset = crate::spk::RECORD_BYTES + (3 + 2 * crate::spk::SPK_SUMMARY_WORDS) * 8 + 16;
    image[target_offset..target_offset + 4].copy_from_slice(&399_i32.to_le_bytes());
    registry.register_image("support-test", image).unwrap();
    registry.register_body(EARTH, "support-test", 399).unwrap();
    // Both aliases must resolve before the summary-only Earth image's
    // missing numerical directory can be sampled.
    assert_eq!(
        earth_helio_true_of_date(&JplLeg::new(&registry), FRAME_J2000),
        Err(Error::SpkUnregisteredBody { body: SUN }),
    );
}

#[test]
fn numbered_asteroid_earth_retardation_requires_convergence() {
    let seen = std::cell::RefCell::new(Vec::new());
    let alternating = |whole, fraction| {
        let mut samples = seen.borrow_mut();
        samples.push((whole, fraction));
        Ok([if samples.len() % 2 == 1 { 2.0 } else { 1.0 }, 0.0, 0.0])
    };
    let error = retarded_earth_minus_sun_samples(
        FRAME_J2000,
        0.125,
        [0.0; 3],
        [1.0, 0.0, 0.0],
        &alternating,
    )
    .unwrap_err();
    assert!(matches!(
        error,
        Error::ConvergenceFailed {
            algorithm: "spk_light_time",
            iterations: 10,
            max_iterations: 10,
            tolerance: 1e-12,
            ..
        }
    ));
    assert_eq!(seen.borrow().len(), 10);
    let seen = seen.borrow();
    assert!(seen.iter().all(|(whole, fraction)| {
        *whole == FRAME_J2000 && *fraction < 0.125 && *fraction > 0.0
    }));
    let calls = std::cell::Cell::new(0);
    assert_eq!(
        retarded_earth_minus_sun_samples(FRAME_J2000, 0.125, [0.0; 3], [1.0, 0.0, 0.0], &|_, _| {
            calls.set(calls.get() + 1);
            Ok([1.0, 0.0, 0.0])
        },)
        .unwrap(),
        [1.0, 0.0, 0.0],
    );
    assert_eq!(calls.get(), 1);
    assert_eq!(
        retarded_earth_minus_sun_samples(FRAME_J2000, 0.125, [0.0; 3], [1.0, 0.0, 0.0], &|_, _| {
            Err(Error::SpkUnregisteredBody { body: EARTH })
        },),
        Err(Error::SpkUnregisteredBody { body: EARTH }),
    );
}

#[test]
fn earth_observation_reports_missing_earth_with_valid_sun_support() {
    let mut image = crate::spk::range_capture_chain_image(60.0, 60.0);
    // The constant Sun record and its directory really cover this epoch.
    let target_offset = crate::spk::RECORD_BYTES + 24 + 16;
    image[target_offset..target_offset + 4].copy_from_slice(&10_i32.to_le_bytes());
    let mut registry = SpkRegistry::new();
    registry.register_image("sun-support-test", image).unwrap();
    registry.register_body(SUN, "sun-support-test", 10).unwrap();
    let valid_sun = registry
        .eval_geometric_ssb_au_two_part_jd(SUN, FRAME_J2000, 30.0 / 86400.0)
        .unwrap();
    assert_eq!(
        valid_sun.0,
        [
            1.0 / 149_597_870.7,
            2.0 / 149_597_870.7,
            3.0 / 149_597_870.7
        ]
    );
    let leg = JplLeg::new(&registry);
    let epoch = FRAME_J2000 + 30.0 / 86400.0;
    assert_eq!(
        retarded_earth_minus_sun(&leg, epoch),
        Err(Error::SpkUnregisteredBody { body: EARTH })
    );
    assert_eq!(
        earth_helio_true_of_date(&leg, epoch),
        Err(Error::SpkUnregisteredBody { body: EARTH })
    );
}

#[test]
fn hypothetical_topo_site_gates_parallax_requests() {
    let site = Some((13.4, 52.5, 43.0));
    // Tropical geocentric topocentric requests serve with a site and
    // refuse without one, matching the reference missing-observer error.
    assert_eq!(hypothetical_topo_site(FLG_TOPOCTR, site).unwrap(), site);
    assert!(matches!(
        hypothetical_topo_site(FLG_TOPOCTR, None),
        Err(Error::ConfigurationMissing {
            missing_config: "observer_location",
            ..
        })
    ));
    assert!(matches!(
        hypothetical_topo_site(FLG_TOPOCTR | FLG_SPEED, None),
        Err(Error::ConfigurationMissing {
            missing_config: "observer_location",
            ..
        })
    ));
    // Parallax-inert requests never require a site: plain geocentric,
    // heliocentric/barycentric (with or without the topocentric bit).
    for flags in [0, FLG_HELCTR | FLG_TOPOCTR, FLG_BARYCTR | FLG_TOPOCTR] {
        assert_eq!(hypothetical_topo_site(flags, None).unwrap(), None);
    }
    // Topocentric parallax requires a site independently of the output frame.
    for flags in [
        FLG_TOPOCTR | FLG_J2000,
        FLG_TOPOCTR | FLG_SIDEREAL,
        FLG_TOPOCTR | FLG_EQUATORIAL,
    ] {
        assert_eq!(hypothetical_topo_site(flags, site).unwrap(), site);
        assert!(matches!(
            hypothetical_topo_site(flags, None),
            Err(Error::ConfigurationMissing {
                missing_config: "observer_location",
                ..
            })
        ));
    }
}

#[test]
fn hypothetical_topo_offset_has_earth_radius_scale() {
    // Berlin site: the offset magnitude is about one Earth radius in AU
    // and the rate is about the diurnal spin over one day; the
    // computation is deterministic and needs no kernel.
    let site = (13.4, 52.5, 43.0);
    let (pos, vel) = hypothetical_topo_offset_j2000(2460000.5, site).unwrap();
    let r = (pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]).sqrt();
    assert!((3.5e-5..5.0e-5).contains(&r), "offset radius {r}");
    let v = (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]).sqrt();
    assert!((1.0e-4..5.0e-4).contains(&v), "offset rate {v}");
    assert_eq!(
        hypothetical_topo_offset_j2000(2460000.5, site).unwrap().0,
        pos
    );
    assert!(matches!(
        hypothetical_topo_offset_j2000(f64::NAN, site),
        Err(Error::NonFiniteArgument { .. })
    ));
}

#[test]
fn symbolic_topo_serves_tropical_parallax_end_to_end() {
    use crate::context::{EphemerisContext, Observer, ProviderMode};
    // TRUEPOS keeps the symbolic construction kernel-free, so the whole
    // tropical path (symbolic core, ERFA nutation, site parallax) runs
    // on an empty registry through the public entry.
    let ctx = EphemerisContext::new()
        .with_mode(ProviderMode::NativeJpl)
        .with_observer(Observer::new(13.4, 52.5, 43.0).unwrap());
    let request = ServingRequest {
        ctx: &ctx,
        delta_t_model: crate::time::DeltaTModel::EspenakMeeus,
        iers_enabled: false,
        observer: None,
        tier_policy: TierPolicy::AttachOrder,
        skip_closed: false,
    };
    let tt = 2460000.5;
    let body = crate::constants::WHITE_MOON;
    let topo_flags = crate::constants::FLG_TRUEPOS | crate::constants::FLG_TOPOCTR;
    let (geo, _) = super::calc_horizons(request, tt, body, crate::constants::FLG_TRUEPOS).unwrap();
    let (topo, retflag) = super::calc_horizons(request, tt, body, topo_flags).unwrap();
    assert_eq!(
        retflag & crate::constants::FLG_TOPOCTR,
        crate::constants::FLG_TOPOCTR
    );
    for x in topo {
        assert!(x.is_finite());
    }
    // Parallax is real (arcseconds to arcminutes for Selena) and can
    // never reach a degree at these distances; an unapplied correction
    // would leave the words exactly equal instead.
    let mut dlon = (topo[0] - geo[0]).rem_euclid(360.0);
    if dlon > 180.0 {
        dlon -= 360.0;
    }
    assert!(
        dlon.abs() > 1e-6 && dlon.abs() < 1.0,
        "parallax shift {dlon}"
    );
    assert!(
        (topo[1] - geo[1]).abs() < 1.0 && (topo[2] - geo[2]).abs() < 1e-3,
        "topo {topo:?} vs geo {geo:?}"
    );
    // Without a configured site the same request is the typed
    // missing-observer refusal, not a kernel or fall-through error.
    let bare = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
    let bare_request = ServingRequest {
        ctx: &bare,
        delta_t_model: crate::time::DeltaTModel::EspenakMeeus,
        iers_enabled: false,
        observer: None,
        tier_policy: TierPolicy::AttachOrder,
        skip_closed: false,
    };
    assert!(matches!(
        super::calc_horizons(bare_request, tt, body, topo_flags),
        Err(Error::ConfigurationMissing {
            missing_config: "observer_location",
            ..
        })
    ));
}

#[test]
fn predicted_topo_gates_tropical_parallax_before_kernel_access() {
    use crate::context::{EphemerisContext, Observer, ProviderMode};
    // Tropical geocentric topocentric requests on the predicted arm
    // reach the parallax gate before any kernel acquisition: without a
    // site the typed missing-observer refusal owns them (previously the
    // fall-through unknown-body refusal, see the red-before bank probe).
    let bare = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
    let bare_request = ServingRequest {
        ctx: &bare,
        delta_t_model: crate::time::DeltaTModel::EspenakMeeus,
        iers_enabled: false,
        observer: None,
        tier_policy: TierPolicy::AttachOrder,
        skip_closed: false,
    };
    let tt = 2460000.5;
    let topo_flags = crate::constants::FLG_TOPOCTR | crate::constants::FLG_SPEED;
    assert!(matches!(
        super::calc_horizons(bare_request, tt, crate::constants::VULCAN, topo_flags),
        Err(Error::ConfigurationMissing {
            missing_config: "observer_location",
            ..
        })
    ));
    // A configured site does not replace missing Sun/Earth kernel support.
    let ctx = EphemerisContext::new()
        .with_mode(ProviderMode::NativeJpl)
        .with_observer(Observer::new(13.4, 52.5, 43.0).unwrap());
    let request = ServingRequest {
        ctx: &ctx,
        delta_t_model: crate::time::DeltaTModel::EspenakMeeus,
        iers_enabled: false,
        observer: None,
        tier_policy: TierPolicy::AttachOrder,
        skip_closed: false,
    };
    assert!(matches!(
        super::calc_horizons(request, tt, crate::constants::VULCAN, topo_flags),
        Err(Error::SpkUnregisteredBody { .. })
    ));
    // Topocentric parallax requires a site independently of the output frame.
    for flags in [
        topo_flags | crate::constants::FLG_J2000,
        topo_flags | crate::constants::FLG_SIDEREAL,
        topo_flags | crate::constants::FLG_EQUATORIAL,
    ] {
        assert!(matches!(
            super::calc_horizons(request, tt, crate::constants::VULCAN, flags),
            Err(Error::SpkUnregisteredBody { .. })
        ));
        assert!(matches!(
            super::calc_horizons(bare_request, tt, crate::constants::VULCAN, flags),
            Err(Error::ConfigurationMissing {
                missing_config: "observer_location",
                ..
            })
        ));
    }
}
