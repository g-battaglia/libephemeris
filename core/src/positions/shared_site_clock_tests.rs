//! Public factory regressions over one fixed historical image and source words.

use super::*;
use crate::context::{EphemerisContext, Observer, ProviderMode, Tier};
use crate::leb::LebFile;
use crate::time::DeltaTModel;
use sha2::{Digest, Sha256};

const TT_WORDS: [u64; 6] = [
    0x406e79ee33ace23e,
    0xbfc898dedbcc8d0e,
    0x400416b548d7c3ac,
    0x3fe6ecf566771000,
    0xbf839cc2b8ca5a00,
    0xbf62b75a0174c000,
];
const UT_WORDS: [u64; 6] = [
    0x406e7a0481c466eb,
    0xbfc89a10d1b1e93e,
    0x400416b0ba5e568a,
    0x3fe6ecd8c2c50000,
    0xbf83b1f37654ba80,
    0xbf62c0187b044000,
];

/// Load the required historical image, never skip or substitute another asset.
fn identified_context(tiered: bool) -> EphemerisContext {
    let path = std::env::var_os("LIBEPHEMERIS_TEST_LEB2")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(crate::test_data::core_leb_path);
    let bytes = std::fs::read(path).expect("provision the historical LEB2 fixture first");
    let hash: String = Sha256::digest(&bytes)
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect();
    assert!(crate::test_data::admitted_core_hash(&(hash)));
    let mut ctx = EphemerisContext::new()
        .with_mode(ProviderMode::SealedLeb)
        .with_tier(if tiered { Tier::Medium } else { Tier::Base })
        .with_observer(Observer::new(11.0, 45.0, 100.0).unwrap());
    if tiered {
        // Identical bytes under two recognized tier labels exercise the real
        // date-aware assembly. This is an assembly control, not acquisition
        // evidence for a different physical dataset or a second source run.
        ctx.attach_leb_file(
            LebFile::from_bytes("ephemeris_medium_planets.leb2".into(), bytes.clone()).unwrap(),
        );
    }
    ctx.attach_leb_file(LebFile::from_bytes("base_core.leb2".into(), bytes).unwrap());
    ctx
}

/// Select the actual shared request without a model, observer or data fallback.
fn request(ctx: &EphemerisContext) -> ServingRequest<'_> {
    ServingRequest {
        ctx,
        delta_t_model: DeltaTModel::EspenakMeeus,
        iers_enabled: false,
        observer: None,
        tier_policy: TierPolicy::AttachOrder,
        skip_closed: false,
    }
}

/// Check the pinned independent values using the physical precision contract.
fn assert_pinned_precision(actual: PositionResult, expected: [u64; 6]) {
    assert_pinned_precision_or_angular_invariant(
        actual,
        expected,
        None,
        (f64::NEG_INFINITY, f64::INFINITY),
    );
}

fn assert_pinned_precision_or_angular_invariant(
    actual: PositionResult,
    expected: [u64; 6],
    angular_at: Option<&dyn Fn(f64) -> [f64; 2]>,
    support: (f64, f64),
) {
    for (slot, (value, word)) in actual.into_iter().zip(expected).enumerate() {
        let reference = f64::from_bits(word);
        let limit = match slot {
            0 | 1 => 0.001 / 3600.0,
            2 => reference.abs() * 1e-8,
            3 | 4 if reference.abs() <= 1e-6 => 1e-6,
            _ => reference.abs() * 1e-6,
        };
        assert!(value.is_finite());
        if matches!(slot, 3 | 4) && (value - reference).abs() > limit {
            if let Some(at) = angular_at {
                // The largest spacing is selected from actual clock-piece
                // ownership, before inspecting any velocity or residual.
                let direction = if support.0 <= -3.0 / 32.0 && support.1 >= 3.0 / 32.0 {
                    0.0
                } else if support.1 >= -support.0 {
                    1.0
                } else {
                    -1.0
                };
                let mut largest = 1.0 / 32.0;
                while if direction == 0.0 {
                    support.0 > -3.0 * largest || support.1 < 3.0 * largest
                } else {
                    let end = direction * 6.0 * largest;
                    end < support.0 || end > support.1
                } {
                    largest /= 2.0;
                    assert!(largest > 0.0);
                }
                let mut estimates = Vec::new();
                for h in [largest / 2.0, largest / 4.0, largest] {
                    let difference = |next: f64, previous: f64| {
                        if slot == 3 {
                            (next - previous + 180.0).rem_euclid(360.0) - 180.0
                        } else {
                            next - previous
                        }
                    };
                    let derivative = if direction == 0.0 {
                        let delta = |n: f64| difference(at(n * h)[slot - 3], at(-n * h)[slot - 3]);
                        (45.0 * delta(1.0) - 9.0 * delta(2.0) + delta(3.0)) / (60.0 * h)
                    } else {
                        let center = at(0.0)[slot - 3];
                        let sum: f64 = [360.0, -450.0, 400.0, -225.0, 72.0, -10.0]
                            .into_iter()
                            .enumerate()
                            .map(|(n, weight)| {
                                let next = at(direction * (n + 1) as f64 * h)[slot - 3];
                                weight * difference(next, center)
                            })
                            .sum();
                        sum / (60.0 * h * direction)
                    };
                    estimates.push(derivative);
                }
                let min = estimates.iter().copied().fold(f64::INFINITY, f64::min);
                let max = estimates.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                assert!(
                    max - min <= limit,
                    "position-derivative estimates do not converge: {estimates:?}"
                );
                for derivative in &estimates {
                    assert!(
                        (value - derivative).abs() < (reference - derivative).abs(),
                        "native {value} is not closer than raw Source {reference} to {derivative}"
                    );
                    assert!(
                        (value - derivative).abs() <= limit,
                        "native {value} exceeds unchanged speed limit {limit} against {derivative}"
                    );
                }
                eprintln!(
                    "raw Source angular rate {reference}, native {value}, independent public-position derivatives {estimates:?}"
                );
                continue;
            }
        }
        assert!(
            (value - reference).abs() <= limit,
            "component {slot}: {value:?} vs pinned {reference:?}; physical limit {limit}"
        );
    }
}

fn assert_clocked_place_precision(
    ctx: &EphemerisContext,
    original_ut: Option<f64>,
    actual: PositionResult,
    expected: [u64; 6],
) {
    let tt = if let Some(ut) = original_ut {
        with_shared_calc_engine_at_entry(request(ctx), Some(ut), |engine| (engine.ut1_to_tt_jd)(ut))
            .unwrap()
    } else {
        2506300.5
    };
    let cfg = ctx.delta_t_config(DeltaTModel::EspenakMeeus, false);
    let positions = |offset| {
        let sample_tt = tt + offset;
        let sample_ut = original_ut.map(|_| {
            super::site_clock::configured_site_sample_ut1(sample_tt, original_ut, &cfg, &|_| {
                panic!("the finite configured clock does not need a table")
            })
            .unwrap()
        });
        context_serving::with_shared_calc_engine_at_entry(request(ctx), sample_ut, |engine| {
            if let Some(ut) = sample_ut {
                calc_ut(engine, ut, MARS, 33026 & !FLG_SPEED)
            } else {
                calc(engine, sample_tt, MARS, 33026 & !FLG_SPEED)
            }
        })
        .map(|(value, _)| [value[0], value[1]])
        .unwrap()
    };
    let support = super::site_clock::configured_site_clock_support(tt, original_ut, &cfg).unwrap();
    assert_pinned_precision_or_angular_invariant(actual, expected, Some(&positions), support);
}

#[test]
fn shared_site_tt_december_uses_actual_source_refined_partner() {
    let ctx = identified_context(false);
    let (actual, flags) = with_shared_calc_outputs(request(&ctx), |engine, cartesian| {
        calc_with_cartesian(engine, cartesian, 2506300.5, MARS, 33026)
    })
    .unwrap();
    assert_eq!(flags, 33026);
    assert_clocked_place_precision(&ctx, None, actual, TT_WORDS);
}

#[test]
fn shared_site_ut_december_retains_actual_original_source_ut() {
    let ctx = identified_context(false);
    let (actual, flags) =
        with_shared_calc_outputs_at_entry(request(&ctx), Some(2506300.5), |engine, cartesian| {
            calc_ut_with_cartesian(engine, cartesian, 2506300.5, MARS, 33026)
        })
        .unwrap();
    assert_eq!(flags, 33026);
    assert_clocked_place_precision(&ctx, Some(2506300.5), actual, UT_WORDS);
}

#[test]
fn shared_site_date_aware_context_matches_both_original_source_entries() {
    let ctx = identified_context(true);
    for (original_ut, words) in [(None, TT_WORDS), (Some(2506300.5), UT_WORDS)] {
        let (actual, flags) = context_serving::with_shared_calc_engine_at_entry(
            request(&ctx),
            original_ut,
            |engine| {
                if original_ut.is_some() {
                    calc_ut(engine, 2506300.5, MARS, 33026)
                } else {
                    calc(engine, 2506300.5, MARS, 33026)
                }
            },
        )
        .unwrap();
        assert_eq!(flags, 33026);
        assert_clocked_place_precision(&ctx, original_ut, actual, words);
    }
}

#[test]
fn shared_site_direct_cartesian_and_spherical_factories_share_entry_clock() {
    let ctx = identified_context(false);
    for original_ut in [None, Some(2506300.5)] {
        let direct =
            with_shared_calc_outputs_at_entry(request(&ctx), original_ut, |engine, cartesian| {
                if original_ut.is_some() {
                    calc_ut_with_cartesian(engine, cartesian, 2506300.5, MARS, 33026 | FLG_XYZ)
                } else {
                    calc_with_cartesian(engine, cartesian, 2506300.5, MARS, 33026 | FLG_XYZ)
                }
            })
            .unwrap();
        assert_eq!(direct.1, 33026 | FLG_XYZ);
        assert!(direct.0.into_iter().all(f64::is_finite));
        let tiered = identified_context(true);
        let spherical = context_serving::with_shared_calc_engine_at_entry(
            request(&tiered),
            original_ut,
            |engine| {
                if original_ut.is_some() {
                    calc_ut(engine, 2506300.5, MARS, 33026 | FLG_XYZ)
                } else {
                    calc(engine, 2506300.5, MARS, 33026 | FLG_XYZ)
                }
            },
        )
        .unwrap();
        assert_eq!(spherical.1, direct.1);
        // The spherical inverse-Jacobian and direct Cartesian channels keep
        // their separate arithmetic. This checks a bounded assembly result;
        // it supplies no new independently acquired source XYZ expectation.
        for (a, b) in direct.0.into_iter().zip(spherical.0) {
            assert!((a - b).abs() <= 1e-8);
        }
    }
}

#[test]
fn shared_site_tt_missing_observer_precedes_unusable_override_clock() {
    let ctx = identified_context(false)
        .without_observer()
        .with_user_delta_t(Some(f64::NAN));
    let error =
        with_shared_calc_engine(request(&ctx), |engine| calc(engine, 2506300.5, MARS, 33026))
            .unwrap_err();
    assert!(matches!(
        error,
        Error::ConfigurationMissing {
            missing_config: "observer_location",
            ..
        }
    ));
}

#[test]
fn shared_site_native_entry_ignores_the_configured_leb_site_clock() {
    let ctx = identified_context(true)
        .with_mode(ProviderMode::NativeJpl)
        .with_user_delta_t(Some(f64::NAN));
    let ut = 2451545.125;
    let expected = crate::time::default_timescale_tt_from_ut1(ut).unwrap();
    let actual =
        context_serving::with_shared_calc_engine_at_entry(request(&ctx), Some(ut), |engine| {
            (engine.ut1_to_tt_jd)(ut)
        })
        .unwrap();
    assert_eq!(actual.to_bits(), expected.to_bits());
}

/// Physical precision check for the independently acquired canonical fallback row.
fn assert_canonical_fallback_precision(actual: PositionResult) {
    const WORDS: [u64; 6] = [
        0x40543d7be029daf2,
        0x400a9602cca064d9,
        0x3fe41d7c4baaeb04,
        0xbfd1c073bee7a000,
        0x3f8cade74114d000,
        0x3f714ee408e53c00,
    ];
    assert_pinned_precision(actual, WORDS);
}

#[test]
fn shared_canonical_tt_nan_calendar_uses_selected_reader_table() {
    let ctx = identified_context(false).with_tidal_acceleration(Some(f64::NAN));
    let (actual, flags) = with_shared_calc_outputs(request(&ctx), |engine, cartesian| {
        calc_with_cartesian(engine, cartesian, 2396760.125, MARS, 33026)
    })
    .unwrap();
    assert_eq!(flags, 33026);
    assert_canonical_fallback_precision(actual);
}

#[test]
fn shared_canonical_tt_tiered_view_uses_its_actual_delta_t_table() {
    let ctx = identified_context(true).with_tidal_acceleration(Some(f64::NAN));
    let (actual, flags) =
        context_serving::with_shared_calc_engine_at_entry(request(&ctx), None, |engine| {
            calc(engine, 2396760.125, MARS, 33026)
        })
        .unwrap();
    assert_eq!(flags, 33026);
    assert_canonical_fallback_precision(actual);
}
