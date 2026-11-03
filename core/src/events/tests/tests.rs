use super::*;
use std::f64::consts::PI;

fn cfg() -> SearchConfig {
    SearchConfig::new(1e-12, 1e-13, 100).unwrap()
}

/// Read the released LEB2 core fixture, refusing unrelated coefficient files.
///
/// `LIBEPHEMERIS_TEST_LEB2` selects an explicit file; otherwise use
/// the current core in the configured test data directory. Admitted SHA-256
/// identities are checked below. Absence is an optional-data skip
/// (`None`), but a present file with the wrong SHA-256 is a fixture error,
/// not a numerical comparison: silent evaluation on unrelated coefficients previously
/// surfaced as misleading word-pin failures on these L28 pairs.
pub(super) fn leb2_core_fixture() -> Option<Vec<u8>> {
    use sha2::{Digest, Sha256};
    let path = std::env::var_os("LIBEPHEMERIS_TEST_LEB2")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(crate::test_data::core_leb_path);
    let bytes = match std::fs::read(&path) {
        Ok(bytes) => bytes,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => return None,
        Err(error) => panic!("cannot read fixture {}: {error}", path.display()),
    };
    assert!(
        crate::test_data::admitted_core_hash(
            &(Sha256::digest(&bytes)
                .iter()
                .map(|byte| format!("{byte:02x}"))
                .collect::<String>())
        ),
        "LEB2 core fixture requires the released coefficient file: {}",
        path.display()
    );
    Some(bytes)
}

#[test]
fn config_rejects_invalid_fields() {
    assert_eq!(
        SearchConfig::new(0.0, 1e-9, 10).unwrap_err(),
        Error::InvalidSearchTolerance { value: 0.0 }
    );
    assert!(matches!(
        SearchConfig::new(f64::NAN, 1e-9, 10).unwrap_err(),
        Error::InvalidSearchTolerance { value } if value.is_nan()
    ));
    assert_eq!(
        SearchConfig::new(1e-9, -1.0, 10).unwrap_err(),
        Error::InvalidSearchStep { value: -1.0 }
    );
    assert_eq!(
        SearchConfig::new(1e-9, 1e-9, 0).unwrap_err(),
        Error::InvalidSearchIterations { max_iter: 0 }
    );
    assert_eq!(
        NewtonConfig::new(1e-9, 0.0, 10).unwrap_err(),
        Error::InvalidSearchStep { value: 0.0 }
    );
}

#[test]
fn interval_gate_rejects_reversed_and_non_finite() {
    let mut f = |x: f64| x;
    assert_eq!(
        check_bracket(&mut f, 2.0, 1.0, cfg()).unwrap_err(),
        Error::InvalidSearchInterval { a: 2.0, b: 1.0 }
    );
    assert_eq!(
        check_bracket(&mut f, 1.0, 1.0, cfg()).unwrap_err(),
        Error::InvalidSearchInterval { a: 1.0, b: 1.0 }
    );
    assert!(matches!(
        check_bracket(&mut f, f64::NAN, 1.0, cfg()).unwrap_err(),
        Error::InvalidSearchInterval { a, b: 1.0 } if a.is_nan()
    ));
    assert_eq!(
        scan_first_bracket(&mut f, 0.0, 1.0, 0, 1e-9).unwrap_err(),
        Error::InvalidScanSamples { count: 0 }
    );
}

#[test]
fn non_finite_values_are_errors_not_outcomes() {
    let mut f = |x: f64| if x == 1.0 { f64::NAN } else { x - 2.0 };
    assert_eq!(
        check_bracket(&mut f, 1.0, 3.0, cfg()).unwrap_err(),
        Error::NonFiniteSearchValue { x: 1.0 }
    );
    // Interior non-finite values abort the refinement: endpoints are
    // finite with a sign change, but the first midpoint (1.5) is not.
    let mut g = |x: f64| if x == 1.5 { f64::INFINITY } else { x - 2.0 };
    assert_eq!(
        bisection(&mut g, 0.0, 3.0, cfg()).unwrap_err(),
        Error::NonFiniteSearchValue { x: 1.5 }
    );
    let mut h = |_: f64| (f64::NAN, 1.0);
    assert_eq!(
        newton(&mut h, 0.0, NewtonConfig::new(1e-12, 1e-12, 10).unwrap()).unwrap_err(),
        Error::NonFiniteSearchValue { x: 0.0 }
    );
}

#[test]
fn bisection_finds_sqrt2() {
    let mut f = |x: f64| x * x - 2.0;
    match bisection(&mut f, 1.0, 2.0, cfg()).unwrap() {
        RootOutcome::Converged { x, fx, .. } => {
            assert!((x - 2f64.sqrt()).abs() <= 1e-12);
            assert!(fx.abs() <= 1e-12);
        }
        other => panic!("expected convergence, got {other:?}"),
    }
}

#[test]
fn brent_finds_sqrt2_from_either_side() {
    for (a, b) in [(1.0, 2.0), (0.5, 4.0)] {
        let mut f = |x: f64| x * x - 2.0;
        match brent(&mut f, a, b, cfg()).unwrap() {
            RootOutcome::Converged { x, fx, .. } => {
                assert!((x - 2f64.sqrt()).abs() <= 1e-12, "bracket [{a}, {b}]");
                assert!(fx.abs() <= 1e-12);
            }
            other => panic!("expected convergence, got {other:?}"),
        }
    }
}

#[test]
fn endpoint_roots_need_no_iteration() {
    let mut f = |x: f64| x - 1.0;
    assert!(matches!(
        bisection(&mut f, 1.0, 2.0, cfg()).unwrap(),
        RootOutcome::EndpointRoot { x: 1.0, .. }
    ));
    assert!(matches!(
        brent(&mut f, 0.0, 1.0, cfg()).unwrap(),
        RootOutcome::EndpointRoot { x: 1.0, .. }
    ));
    assert!(matches!(
        check_bracket(&mut f, 0.0, 1.0, cfg()).unwrap(),
        BracketOutcome::EndpointRoot { .. }
    ));
}

#[test]
fn same_sign_interval_is_no_bracket() {
    // No sign change: x^2 + 1 has no real root anywhere.
    let mut f = |x: f64| x * x + 1.0;
    assert!(matches!(
        check_bracket(&mut f, -1.0, 1.0, cfg()).unwrap(),
        BracketOutcome::NoBracket { .. }
    ));
    assert!(matches!(
        brent(&mut f, -1.0, 1.0, cfg()).unwrap(),
        RootOutcome::NoBracket { .. }
    ));
}

#[test]
fn tangent_root_has_no_sign_change() {
    // Even-multiplicity (tangent) root: bracketed methods honestly
    // report NoBracket; Newton from nearby converges.
    let mut f = |x: f64| x * x;
    assert!(matches!(
        check_bracket(&mut f, -1.0, 1.0, cfg()).unwrap(),
        BracketOutcome::NoBracket { .. }
    ));
    let mut fdf = |x: f64| (x * x, 2.0 * x);
    match newton(&mut fdf, 0.5, NewtonConfig::new(1e-14, 1e-15, 100).unwrap()).unwrap() {
        RootOutcome::Converged { x, .. } => assert!(x.abs() <= 1e-6),
        other => panic!("expected Newton convergence, got {other:?}"),
    }
}

#[test]
fn newton_refuses_near_zero_derivative() {
    // Station analogue: starting exactly at the extremum of x^2 - 2
    // leaves df = 0, so Newton stalls instead of dividing by zero.
    let mut fdf = |x: f64| (x * x - 2.0, 2.0 * x);
    assert!(matches!(
        newton(&mut fdf, 0.0, NewtonConfig::new(1e-12, 1e-9, 50).unwrap()).unwrap(),
        RootOutcome::NonConverged { .. }
    ));
}

#[test]
fn exhausted_cap_reports_non_converged_with_best_estimate() {
    let mut f = |x: f64| x * x - 2.0;
    let tiny = SearchConfig::new(1e-300, 1e-300, 1).unwrap();
    match brent(&mut f, 1.0, 2.0, tiny).unwrap() {
        RootOutcome::NonConverged {
            best_x, best_fx, ..
        } => {
            assert!(best_x.is_finite() && best_fx.is_finite());
            assert!((1.0..=2.0).contains(&best_x));
        }
        other => panic!("expected non-convergence, got {other:?}"),
    }
}

#[test]
fn scan_isolates_first_sine_root() {
    let mut f = |x: f64| x.sin();
    match scan_first_bracket(&mut f, 2.0, 4.0 * PI, 16, 1e-12).unwrap() {
        ScanOutcome::Bracket { a, b, .. } => {
            assert!(a < PI && PI < b, "first bracket [{a}, {b}] misses pi");
            let mut g = |x: f64| x.sin();
            match brent(&mut g, a, b, cfg()).unwrap() {
                RootOutcome::Converged { x, .. } => {
                    assert!((x - PI).abs() <= 1e-9);
                }
                other => panic!("expected convergence, got {other:?}"),
            }
        }
        other => panic!("expected a bracket, got {other:?}"),
    }
}

#[test]
fn scan_reports_no_root_and_exact_sample_hits() {
    // sin is strictly positive on [0.5, 1.0]: genuine no-root scan.
    let mut f = |x: f64| x.sin();
    assert!(matches!(
        scan_first_bracket(&mut f, 0.5, 1.0, 8, 1e-12).unwrap(),
        ScanOutcome::NoRoot
    ));
    // A sample landing exactly on the root is an endpoint hit.
    let mut g = |x: f64| x - 0.75;
    assert!(matches!(
        scan_first_bracket(&mut g, 0.0, 1.0, 4, 1e-12).unwrap(),
        ScanOutcome::EndpointRoot { .. }
    ));
}

#[test]
fn repeated_root_bracket_still_converges() {
    // Odd-multiplicity repeated root (x-1)^3: sign change exists, so
    // Brent converges, more slowly than for a simple root.
    let mut f = |x: f64| (x - 1.0).powi(3);
    match brent(&mut f, 0.0, 2.0, cfg()).unwrap() {
        RootOutcome::Converged { x, fx, .. } => {
            assert!((x - 1.0).abs() <= 1e-6);
            assert!(fx.abs() <= 1e-12);
        }
        other => panic!("expected convergence, got {other:?}"),
    }
}

// --: longitude-crossing layer over analytic providers --------

/// Linear provider: `lon = lon0 + rate * (jd - jd0)`, exact speed.
fn linear(
    lon0: f64,
    rate: f64,
    jd0: f64,
) -> impl FnMut(f64) -> std::result::Result<(f64, f64), Error> {
    move |jd: f64| Ok((lon0 + rate * (jd - jd0), rate))
}

#[test]
fn residual_and_target_pins() {
    assert_eq!(wrapped_residual_deg(10.0, 10.0), 0.0);
    assert_eq!(wrapped_residual_deg(370.0, 10.0), 0.0);
    assert_eq!(wrapped_residual_deg(0.0, 359.0), 1.0);
    assert_eq!(wrapped_residual_deg(359.0, 0.0), -1.0);
    // Antipode sits at the ±180 discontinuity, never near 0.
    assert!((wrapped_residual_deg(190.0, 10.0).abs() - 180.0).abs() < 1e-9);
    assert_eq!(normalize_target_deg(720.0 + 45.0), 45.0);
    assert_eq!(adaptive_max_iterations(1.0), NR_MAX_ITER_PLANET);
    assert_eq!(adaptive_max_iterations(0.05), NR_MAX_ITER_SLOW);
    assert_eq!(adaptive_max_iterations(0.005), NR_MAX_ITER_SLOW);
    assert_eq!(adaptive_max_iterations(0.0005), NR_MAX_ITER_STATION);
    assert!(is_near_station(0.005));
    assert!(!is_near_station(0.5));
}

#[test]
fn seed_diff_direction_semantics() {
    // Forward: just-passed target maps near 360 (next cycle); on-target is now.
    assert!((seed_diff_deg(10.0, 9.0, SearchDirection::Forward, 1.0) - 359.0).abs() < 1e-12);
    assert_eq!(
        seed_diff_deg(10.0, 10.0, SearchDirection::Forward, 1.0),
        0.0
    );
    // Forward retrograde guard keeps the sign contract.
    assert_eq!(
        seed_diff_deg(10.0, 20.0, SearchDirection::Forward, -0.5),
        -350.0
    );
    // Backward: strictly-past — on-target resolves a full cycle back.
    assert_eq!(
        seed_diff_deg(10.0, 10.0, SearchDirection::Backward, 1.0),
        -360.0
    );
    // Backward just-after-event: small negative distance, not a full cycle.
    let d = seed_diff_deg(10.0, 9.999, SearchDirection::Backward, 1.0);
    assert!(d < 0.0 && d > -1.0, "got {d}");
}

#[test]
fn sun_crossing_linear_forward_and_backward() {
    let jd0 = 2_451_545.0;
    let rate = SUN_MEAN_RATE_DEG_DAY;
    // Forward: 10 deg ahead resolves in 10/rate days.
    let mut p = linear(100.0, rate, jd0);
    let hit = sun_crossing(&mut p, 110.0, jd0, SearchDirection::Forward).unwrap();
    assert!((hit - (jd0 + 10.0 / rate)).abs() < 1e-6, "got {hit}");
    // Wrap: start 359, target 1 → 2 deg ahead.
    let mut p = linear(359.0, 1.0, jd0);
    let hit = sun_crossing(&mut p, 1.0, jd0, SearchDirection::Forward).unwrap();
    assert!((hit - (jd0 + 2.0)).abs() < 1e-9, "got {hit}");
    // Backward strictly-past: on-target start goes a full cycle back.
    let mut p = linear(100.0, rate, jd0);
    let hit = sun_crossing(&mut p, 100.0, jd0, SearchDirection::Backward).unwrap();
    assert!((hit - (jd0 - 360.0 / rate)).abs() < 1e-6, "got {hit}");
}

#[test]
fn moon_crossing_fast_rate() {
    let jd0 = 2_451_545.0;
    let mut p = linear(0.0, MOON_MEAN_RATE_DEG_DAY, jd0);
    let hit = moon_crossing(&mut p, 90.0, jd0, SearchDirection::Forward).unwrap();
    assert!(
        (hit - (jd0 + 90.0 / MOON_MEAN_RATE_DEG_DAY)).abs() < 1e-6,
        "got {hit}"
    );
}

#[test]
fn moon_crossing_uses_release_point_one_floor_fix_c03() {
    // FIX-C03 regression: release `mooncross` substitutes the mean
    // rate (13.176 deg/day) for the Newton STEP once |speed| < 0.1
    // mid-iteration; the seed substitutes only on an exact-zero
    // start speed (verified at the pin: `if speed == 0` guards the
    // seed, `if abs(speed) < 0.1` guards the step). The old wiring
    // passed the shared 0.01 station threshold instead, so provider
    // speeds in [0.01, 0.1) took the raw-speed Newton path
    // mid-iteration.
    //
    // Two-phase probe: fast at start (seed lands ~0.76 d out), slow
    // (0.05 deg/day) everywhere else. With the 0.1 floor every step
    // divides by 13.176 while the true slope is 0.05, so the walk
    // advances ~0.75 d per step and trips the 31-day divergence
    // gate just past it. With the old 0.01 floor the first step
    // would divide by 0.05 (~199 d) and overshoot the gate on the
    // very first step instead.
    let jd0 = 2_451_545.0;
    let slow = |jd: f64| -> Result<(f64, f64)> {
        if jd == jd0 {
            Ok((80.0, MOON_MEAN_RATE_DEG_DAY))
        } else {
            Ok((80.0 + 0.05 * (jd - jd0), 0.05))
        }
    };
    let mut p = slow;
    let err = moon_crossing(&mut p, 90.0, jd0, SearchDirection::Forward)
        .expect_err("slow-phase probe must diverge");
    match err {
        Error::CrossingDiverged { start, estimate } => {
            assert_eq!(start, jd0, "start {start}");
            let jd_guess = jd0 + 10.0 / MOON_MEAN_RATE_DEG_DAY;
            let over = (estimate - jd_guess).abs();
            assert!(
                over > 31.0 && over < 40.0,
                "estimate {estimate} (over {over})"
            );
        }
        other => panic!("expected divergence, got {other:?}"),
    }
}

#[test]
fn newton_correction_steps_toward_target_fix_evt02() {
    // Stepped away from the root and diverged (exposed by the
    let jd0 = 2_451_545.0;
    let wobble = |jd: f64| -> Result<(f64, f64)> {
        let dt = jd - jd0;
        let lon = 100.0 + SUN_MEAN_RATE_DEG_DAY * dt + 0.05 * (dt / 3.0).sin();
        let speed = SUN_MEAN_RATE_DEG_DAY + 0.05 * (1.0 / 3.0) * (dt / 3.0).cos();
        Ok((lon, speed))
    };
    let mut p = wobble;
    let hit = sun_crossing(&mut p, 110.0, jd0, SearchDirection::Forward).unwrap();
    // About 10 days ahead (the wobble shifts it by < 0.1 d), forward.
    assert!(hit > jd0 && hit < jd0 + 30.0, "got {hit}");
    let q = wobble;
    let (lon_hit, _) = q(hit).unwrap();
    assert!(
        newton_residual_deg(lon_hit, 110.0).abs() < CROSSING_TOLERANCE_DEG,
        "lon {lon_hit}"
    );
}

#[test]
fn planet_crossing_direct_prograde() {
    // Mars-like direct motion: Newton seed lands in the right basin.
    let jd0 = 2_451_545.0;
    let mut p = linear(10.0, 0.5, jd0);
    let hit = planet_crossing(&mut p, 4, 20.0, jd0).unwrap();
    assert!((hit - (jd0 + 20.0)).abs() < 1e-6, "got {hit}");
}

#[test]
fn planet_crossing_retrograde_dip_prescan() {
    // Mercury-rate body starting retrograde with the target just behind:
    // Newton from a linear seed would project a full orbit ahead, but the
    // imminent retrograde dip reaches it within days.
    let jd0 = 2_451_545.0;
    let mut lon = move |jd: f64| {
        let t = jd - jd0;
        let lon = 10.0 + 0.5 * t - 8.0 * (2.0 * std::f64::consts::PI * t / 40.0).sin();
        let speed = 0.5
            - 8.0
                * (2.0 * std::f64::consts::PI / 40.0)
                * (2.0 * std::f64::consts::PI * t / 40.0).cos();
        Ok((lon, speed))
    };
    let hit = planet_crossing(&mut lon, 2, 9.5, jd0).unwrap();
    assert!(
        hit >= jd0 - 1e-6 && hit < jd0 + 60.0,
        "expected the imminent dip crossing, got {hit}"
    );
    let (lon_hit, _) = lon(hit).unwrap();
    assert!(wrapped_residual_deg(lon_hit, 9.5).abs() < CROSSING_TOLERANCE_DEG);
}

#[test]
fn planet_crossing_station_brent_path() {
    // Near-station start (speed 0.005 < 0.01): Brent over the bracketed
    // window refines the crossing 10 days out.
    let jd0 = 2_451_545.0;
    let mut p = linear(10.0, 0.005, jd0);
    let hit = planet_crossing(&mut p, 6, 10.05, jd0).unwrap();
    assert!((hit - (jd0 + 10.0)).abs() < 1e-4, "got {hit}");
}

// -- EVT-02H: heliocentric crossings over analytic providers --------

#[test]
fn helio_speed_table_pins() {
    assert_eq!(helio_typical_speed_deg_day(2), 4.09);
    assert_eq!(helio_typical_speed_deg_day(3), 1.60);
    assert_eq!(helio_typical_speed_deg_day(14), 0.986);
    assert_eq!(helio_typical_speed_deg_day(4), 0.524);
    assert_eq!(helio_typical_speed_deg_day(5), 0.083);
    assert_eq!(helio_typical_speed_deg_day(6), 0.034);
    assert_eq!(helio_typical_speed_deg_day(7), 0.012);
    assert_eq!(helio_typical_speed_deg_day(8), 0.006);
    assert_eq!(helio_typical_speed_deg_day(9), 0.004);
    assert_eq!(helio_typical_speed_deg_day(1000), 0.0);
    assert!(helio_speed_known(14));
    assert!(helio_speed_known(9));
    assert!(!helio_speed_known(1000));
    assert!(!helio_speed_known(10));
}

#[test]
fn helio_crossing_linear_forward() {
    // Mars-rate direct motion: Newton seed lands in the right basin.
    let jd0 = 2_451_545.0;
    let mut p = linear(10.0, 0.524, jd0);
    let hit = helio_crossing(&mut p, 4, 20.0, jd0, SearchDirection::Forward).unwrap();
    assert!((hit - (jd0 + 10.0 / 0.524)).abs() < 1e-6, "got {hit}");
    // Wrap: start 359, target 1 → 2 deg ahead.
    let mut p = linear(359.0, 1.60, jd0);
    let hit = helio_crossing(&mut p, 3, 1.0, jd0, SearchDirection::Forward).unwrap();
    assert!((hit - (jd0 + 2.0 / 1.60)).abs() < 1e-9, "got {hit}");
}

#[test]
fn helio_crossing_backward_strictly_past() {
    // On-target start resolves a full cycle back, never the start itself.
    let jd0 = 2_451_545.0;
    let rate = 0.524;
    let mut p = linear(100.0, rate, jd0);
    let hit = helio_crossing(&mut p, 4, 100.0, jd0, SearchDirection::Backward).unwrap();
    assert!((hit - (jd0 - 360.0 / rate)).abs() < 1e-6, "got {hit}");
}

#[test]
fn helio_crossing_out_of_table_bracketed() {
    // Body outside the calibrated table: deterministic
    // bracket-plus-Brent on the monotonic longitude.
    let jd0 = 2_451_545.0;
    let mut p = linear(10.0, 0.02, jd0);
    let hit = helio_crossing(&mut p, 1000, 20.0, jd0, SearchDirection::Forward).unwrap();
    assert!((hit - (jd0 + 10.0 / 0.02)).abs() < 1e-4, "got {hit}");
    // Backward from just past the target: the bracketed path finds
    // the adjacent strictly-past crossing (1 deg back at 0.02/d),
    // not a full cycle back.
    let mut p = linear(10.0, 0.02, jd0);
    let hit = helio_crossing(&mut p, 1000, 9.0, jd0, SearchDirection::Backward).unwrap();
    assert!((hit - (jd0 - 1.0 / 0.02)).abs() < 1e-4, "got {hit}");
}

#[test]
fn helio_crossing_near_station_brent_path() {
    // Near-station start (speed 0.005 < 0.01): Brent over the bracketed
    // window refines the crossing.
    let jd0 = 2_451_545.0;
    let mut p = linear(10.0, 0.005, jd0);
    let hit = helio_crossing(&mut p, 6, 10.05, jd0, SearchDirection::Forward).unwrap();
    assert!((hit - (jd0 + 10.0)).abs() < 1e-4, "got {hit}");
}

#[test]
fn helio_crossing_constant_provider_diverges() {
    // Frozen longitude off-target: the Newton walk escapes the
    // revolution-scaled guard.
    let mut p = |_: f64| Ok((10.0, 0.001));
    match helio_crossing(&mut p, 4, 20.0, 2_451_545.0, SearchDirection::Forward).unwrap_err() {
        Error::CrossingDiverged { .. } => {}
        other => panic!("expected divergence, got {other:?}"),
    }
}

#[test]
fn nearest_past_helio_repairs_wrong_basin() {
    // Newton landed a whole revolution too far back; the repair
    // returns the crossing adjacent to the start.
    let jd0 = 2_451_545.0;
    let rate = 0.524;
    let period = 360.0 / rate;
    let mut p = linear(10.0, rate, jd0);
    let jd_newton = jd0 - 10.0 / rate - period;
    let hit =
        nearest_past_helio_crossing(&mut p, 0.0, jd0, jd_newton, period, CROSSING_TOLERANCE_DEG)
            .unwrap();
    assert!((hit - (jd0 - 10.0 / rate)).abs() < 1e-4, "got {hit}");
}

#[test]
fn first_forward_helio_repairs_past_newton() {
    // Start 5 deg past the target: the first crossing ahead is the
    // next revolution, not the just-passed root.
    let jd0 = 2_451_545.0;
    let rate = 0.524;
    let period = 360.0 / rate;
    let mut p = linear(5.0, rate, jd0);
    let hit =
        first_forward_helio_crossing(&mut p, 0.0, jd0, period, CROSSING_TOLERANCE_DEG).unwrap();
    assert!((hit - (jd0 + 355.0 / rate)).abs() < 1e-4, "got {hit}");
}

#[test]
fn bracket_scan_rejects_antipodal_wrap() {
    // Residual oscillating ±1 deg around the antipode: sign changes jump
    // ~360 deg and must never bracket.
    let mut p = move |jd: f64| {
        let lon = 190.0 + (jd * 2.0 * std::f64::consts::PI).sin();
        Ok((lon, 1.0))
    };
    assert_eq!(
        find_crossing_bracket(&mut p, 10.0, 0.0, 5.0, 50, false).unwrap_err(),
        Error::NoCrossingInInterval {
            lo: 0.0,
            hi: 5.0,
            target: 10.0,
        }
    );
}

#[test]
fn constant_provider_diverges_honestly() {
    // No motion at all: no bracket anywhere, Newton cannot converge.
    let mut p = |_: f64| Ok((10.0, 0.0));
    match planet_crossing(&mut p, 4, 20.0, 2_451_545.0).unwrap_err() {
        Error::CrossingDiverged { .. } | Error::NoCrossingInInterval { .. } => {}
        other => panic!("expected a divergence outcome, got {other:?}"),
    }
}

#[test]
fn non_finite_provider_is_error() {
    let mut p = |_: f64| Ok((f64::NAN, 1.0));
    assert_eq!(
        sun_crossing(&mut p, 10.0, 2_451_545.0, SearchDirection::Forward).unwrap_err(),
        Error::NonFiniteSearchValue { x: 2_451_545.0 }
    );
}

#[test]
fn backward_exclusion_sizing() {
    // Slow bodies need a wider exclusion than the 1e-6 day floor.
    assert!(
        (backward_start_exclusion_days(0.004, CROSSING_TOLERANCE_DEG)
            - 10.0 * CROSSING_TOLERANCE_DEG / 0.004)
            .abs()
            < 1e-12
    );
    assert_eq!(
        backward_start_exclusion_days(0.0, CROSSING_TOLERANCE_DEG),
        1.0
    );
    assert_eq!(
        backward_start_exclusion_days(13.0, CROSSING_TOLERANCE_DEG),
        1e-6
    );
}

#[test]
fn adaptive_scan_finds_first_of_triple() {
    // Retrograde triple crossing: the adaptive scan returns the FIRST one.
    let jd0 = 2_451_545.0;
    let mut lon = move |jd: f64| {
        let t = jd - jd0;
        let lon = 0.5 * t - 6.0 * (2.0 * std::f64::consts::PI * t / 40.0).sin();
        let speed = 0.5
            - 6.0
                * (2.0 * std::f64::consts::PI / 40.0)
                * (2.0 * std::f64::consts::PI * t / 40.0).cos();
        Ok((lon, speed))
    };
    let first =
        first_forward_crossing_adaptive(&mut lon, 0.2, jd0, jd0 + 120.0, 2.5, 16.0, 2.0, 300.0)
            .unwrap()
            .expect("a crossing must exist");
    assert!(first >= jd0 - 1e-9, "got {first}");
    // No earlier crossing between start and the returned one.
    let earlier =
        first_forward_crossing_adaptive(&mut lon, 0.2, jd0, first - 0.6, 2.5, 16.0, 2.0, 300.0)
            .unwrap();
    assert!(earlier.is_none(), "a nearer crossing exists: {earlier:?}");
    let (lon_hit, _) = lon(first).unwrap();
    assert!(wrapped_residual_deg(lon_hit, 0.2).abs() < 1e-9);
}

// --: stations / retrograde / lunar nodes over analytic providers

/// Cosine speed: stations where `cos` crosses zero, alternating SR/SD.
fn cosine_speed(
    amplitude: f64,
    period_days: f64,
    jd0: f64,
) -> impl FnMut(f64) -> std::result::Result<f64, Error> {
    move |jd: f64| {
        let phase = 2.0 * std::f64::consts::PI * (jd - jd0) / period_days;
        Ok(amplitude * phase.cos())
    }
}

#[test]
fn station_predicates_and_body_guard() {
    assert!(is_retrograde_speed(-0.5));
    assert!(!is_retrograde_speed(0.0));
    assert!(!is_retrograde_speed(0.5));
    assert_eq!(
        classify_station(0.5, -0.5),
        StationKind::StationaryRetrograde
    );
    assert_eq!(classify_station(-0.5, 0.5), StationKind::StationaryDirect);
    // Release else-branch: degenerate trends are SD, preserved.
    assert_eq!(classify_station(-0.5, -0.5), StationKind::StationaryDirect);
    assert!(station_matches_request(
        StationKind::StationaryRetrograde,
        StationRequest::Any
    ));
    assert!(station_matches_request(
        StationKind::StationaryRetrograde,
        StationRequest::Retrograde
    ));
    assert!(!station_matches_request(
        StationKind::StationaryRetrograde,
        StationRequest::Direct
    ));
    assert_eq!(
        check_station_body(0).unwrap_err(),
        Error::BodyHasNoStations { body: 0 }
    );
    assert_eq!(
        check_station_body(1).unwrap_err(),
        Error::BodyHasNoStations { body: 1 }
    );
    assert!(check_station_body(4).is_ok());
    assert_eq!(
        find_station(
            &mut cosine_speed(0.3, 60.0, 0.0),
            0,
            0.0,
            StationRequest::Any
        )
        .unwrap_err(),
        Error::BodyHasNoStations { body: 0 }
    );
}

#[test]
fn brent_station_refines_linear_speed() {
    // v(t) = 0.1 * (t - 100): station exactly at t = 100.
    let mut v = |jd: f64| Ok(0.1 * (jd - 100.0));
    let hit = brent_station(
        &mut v,
        4,
        90.0,
        130.0,
        STATION_VELOCITY_TOLERANCE_DEG_DAY,
        100,
    )
    .unwrap();
    assert!((hit - 100.0).abs() < 1e-6, "got {hit}");
    assert!(checked_speed(&mut v, hit).unwrap().abs() < STATION_VELOCITY_TOLERANCE_DEG_DAY);
    // Unbracketed interval (same-sign speeds) is an honest no-station.
    let mut v = |jd: f64| Ok(0.1 * (jd - 100.0));
    assert_eq!(
        brent_station(
            &mut v,
            4,
            110.0,
            130.0,
            STATION_VELOCITY_TOLERANCE_DEG_DAY,
            100
        )
        .unwrap_err(),
        Error::NoStationInInterval {
            lo: 110.0,
            hi: 130.0,
            planet: 4,
        }
    );
    // Invalid configs are errors, never silent clamps.
    let mut v = |_: f64| Ok(0.0);
    assert_eq!(
        brent_station(&mut v, 4, 1.0, 1.0, 1e-6, 10).unwrap_err(),
        Error::InvalidSearchInterval { a: 1.0, b: 1.0 }
    );
    assert_eq!(
        brent_station(&mut v, 4, 0.0, 1.0, 0.0, 10).unwrap_err(),
        Error::InvalidSearchTolerance { value: 0.0 }
    );
}

#[test]
fn station_bracket_scan_isolates_first_sign_change() {
    let jd0 = 2_451_545.0;
    // cos period 60 d: stations at jd0+15 (SR), jd0+45 (SD),...
    let mut v = cosine_speed(0.3, 60.0, jd0);
    let (a, b) = find_station_bracket(&mut v, 4, jd0, jd0 + 60.0, 50).unwrap();
    assert!(a < jd0 + 15.0 && jd0 + 15.0 < b, "bracket [{a}, {b}]");
    // Constant prograde speed: no station anywhere in the window.
    let mut v = |_: f64| Ok(0.5);
    assert_eq!(
        find_station_bracket(&mut v, 4, jd0, jd0 + 100.0, 50).unwrap_err(),
        Error::NoStationInInterval {
            lo: jd0,
            hi: jd0 + 100.0,
            planet: 4,
        }
    );
}

#[test]
fn find_station_forward_with_type_filter() {
    let jd0 = 2_451_545.0;
    // Mars table (synodic 780 d, window 390 d) over a 60-day cosine
    // cycle: first station at jd0+15 is SR.
    let mut v = cosine_speed(0.3, 60.0, jd0);
    let (jd_sr, kind) = find_station(&mut v, 4, jd0, StationRequest::Any).unwrap();
    assert_eq!(kind, StationKind::StationaryRetrograde);
    assert!((jd_sr - (jd0 + 15.0)).abs() < 1e-4, "got {jd_sr}");
    // A Direct request skips the SR and lands on the SD at jd0+45.
    let mut v = cosine_speed(0.3, 60.0, jd0);
    let (jd_sd, kind) = find_station(&mut v, 4, jd0, StationRequest::Direct).unwrap();
    assert_eq!(kind, StationKind::StationaryDirect);
    assert!((jd_sd - (jd0 + 45.0)).abs() < 1e-4, "got {jd_sd}");
    // Constant speed: every window is empty, the full 4-window span
    // is reported honestly.
    let mut v = |_: f64| Ok(0.5);
    assert_eq!(
        find_station(&mut v, 4, jd0, StationRequest::Any).unwrap_err(),
        Error::NoStationInInterval {
            lo: jd0,
            hi: jd0 + 4.0 * 390.0,
            planet: 4,
        }
    );
}

#[test]
fn find_previous_station_is_strictly_past() {
    let jd0 = 2_451_545.0;
    // Stations at jd0+15 (SR), jd0+45 (SD): from jd0+50 the nearest
    // past station is the SD at jd0+45.
    let mut v = cosine_speed(0.3, 60.0, jd0);
    let (jd_past, kind) = find_previous_station(&mut v, 4, jd0 + 50.0).unwrap();
    assert_eq!(kind, StationKind::StationaryDirect);
    assert!((jd_past - (jd0 + 45.0)).abs() < 1e-4, "got {jd_past}");
    assert!(jd_past < jd0 + 50.0);
    // From exactly a station epoch the search still moves strictly
    // backward (the <= gate brackets at the endpoint, Brent returns
    // it only when it is the nearest past sign change).
    let mut v = cosine_speed(0.3, 60.0, jd0);
    let (jd_past, _) = find_previous_station(&mut v, 4, jd0 + 45.0).unwrap();
    assert!(jd_past <= jd0 + 45.0, "got {jd_past}");
}

#[test]
fn nearest_station_picks_closer_side() {
    let jd0 = 2_451_545.0;
    // Stations at jd0+15 (SR), jd0+45 (SD): from jd0+35 the past
    // station (jd0+15, 20 d away) loses to the future one (jd0+45,
    // 10 d away) → future wins, signed +. Motion between SR and SD
    // is retrograde, reported honestly.
    let mut v = cosine_speed(0.3, 60.0, jd0);
    let info = nearest_station(&mut v, 4, jd0 + 35.0).unwrap();
    assert_eq!(info.kind, StationKind::StationaryDirect);
    assert!((info.jd_station - (jd0 + 45.0)).abs() < 1e-4);
    assert!((info.days_to_station - 10.0).abs() < 1e-4);
    assert!(info.currently_retrograde);
    // From jd0+20 the past SR (5 d) beats the future SD (25 d).
    let mut v = cosine_speed(0.3, 60.0, jd0);
    let info = nearest_station(&mut v, 4, jd0 + 20.0).unwrap();
    assert_eq!(info.kind, StationKind::StationaryRetrograde);
    assert!(info.days_to_station < 0.0, "got {}", info.days_to_station);
    // Retrograde at the query epoch is reported honestly.
    let mut v = cosine_speed(0.3, 60.0, jd0);
    let info = nearest_station(&mut v, 4, jd0 + 30.0).unwrap();
    assert!(info.currently_retrograde, "v(jd0+30) must be negative");
}

#[test]
fn next_retrograde_returns_ordered_pair() {
    let jd0 = 2_451_545.0;
    // Direct at jd0 (cos = +0.3): next SR at jd0+15, next SD at jd0+45.
    let mut v = cosine_speed(0.3, 60.0, jd0);
    let (sr, sd) = next_retrograde(&mut v, 4, jd0).unwrap();
    assert!((sr - (jd0 + 15.0)).abs() < 1e-4, "got {sr}");
    assert!((sd - (jd0 + 45.0)).abs() < 1e-4, "got {sd}");
    assert!(sr < sd);
    // Currently retrograde (jd0+30): the returned pair is the NEXT
    // period, both strictly after the current one ends.
    let mut v = cosine_speed(0.3, 60.0, jd0);
    let (sr, sd) = next_retrograde(&mut v, 4, jd0 + 30.0).unwrap();
    assert!(sr > jd0 + 45.0, "got {sr}");
    assert!(sd > sr, "got ({sr}, {sd})");
    assert!((sr - (jd0 + 75.0)).abs() < 1.0, "got {sr}");
}

#[test]
fn station_central_velocity_handles_wrap() {
    // Longitude crossing 360→0 across the station epoch: the wrapped
    // difference is +2 deg over 2 days → +1 deg/day.
    let mut s = |jd: f64| Ok((359.0 + (jd - 100.0), 5.0, 1.0));
    let (v_lon, v_lat, v_dist) = station_central_velocity(&mut s, 100.0, 1.0).unwrap();
    assert!((v_lon - 1.0).abs() < 1e-12, "got {v_lon}");
    assert!(v_lat.abs() < 1e-12);
    assert!(v_dist.abs() < 1e-12);
    // Exact-wrap case: 359.5 → 0.5 is +1 deg, not −359.
    let mut s = |jd: f64| {
        let lon = if jd < 100.0 { 359.5 } else { 0.5 };
        Ok((lon, 0.0, 1.0))
    };
    let (v_lon, _, _) = station_central_velocity(&mut s, 100.0, 1.0).unwrap();
    assert!((v_lon - 0.5).abs() < 1e-12, "got {v_lon}");
    assert_eq!(
        station_central_velocity(&mut s, 100.0, 0.0).unwrap_err(),
        Error::InvalidSearchStep { value: 0.0 }
    );
    let mut bad = |_: f64| Ok((f64::NAN, 0.0, 1.0));
    assert_eq!(
        station_central_velocity(&mut bad, 100.0, 1.0).unwrap_err(),
        Error::NonFiniteSearchValue { x: 99.0 }
    );
}

/// Synthetic Moon: latitude sine (nodes every half period), longitude
/// linear, latitude speed analytic.
fn synthetic_moon(
    period_days: f64,
    jd0: f64,
) -> impl FnMut(f64) -> std::result::Result<(f64, f64, f64), Error> {
    move |jd: f64| {
        let phase = 2.0 * std::f64::consts::PI * (jd - jd0) / period_days;
        Ok((
            100.0 + 13.176 * (jd - jd0),
            5.0 * phase.sin(),
            5.0 * 2.0 * std::f64::consts::PI / period_days * phase.cos(),
        ))
    }
}

#[test]
fn moon_node_forward_finds_first_crossing() {
    let jd0 = 2_451_545.0;
    // Nodes at jd0, jd0+13.605,... (27.21-day period): from jd0+1
    // the first crossing is jd0+13.605.
    let mut moon = synthetic_moon(27.21, jd0);
    let hit = moon_node_crossing(&mut moon, jd0 + 1.0, SearchDirection::Forward).unwrap();
    assert!((hit.jd - (jd0 + 13.605)).abs() < 1e-6, "got {}", hit.jd);
    assert!(hit.latitude_deg.abs() < MOON_NODE_TOLERANCE_DEG);
    let mut moon = synthetic_moon(27.21, jd0);
    let (lon, _, _) = moon(hit.jd).unwrap();
    assert!((hit.longitude_deg - lon).abs() < 1e-9);
}

#[test]
fn moon_node_backward_is_strictly_past() {
    let jd0 = 2_451_545.0;
    // On-node start resolves a full half-month back, never itself.
    let mut moon = synthetic_moon(27.21, jd0);
    let hit = moon_node_crossing(&mut moon, jd0, SearchDirection::Backward).unwrap();
    assert!((hit.jd - (jd0 - 13.605)).abs() < 1e-6, "got {}", hit.jd);
    // Just after a node, backward finds that same node.
    let mut moon = synthetic_moon(27.21, jd0);
    let hit = moon_node_crossing(&mut moon, jd0 + 1.0, SearchDirection::Backward).unwrap();
    assert!((hit.jd - jd0).abs() < 1e-6, "got {}", hit.jd);
}

#[test]
fn moon_node_at_node_start_steps_forward() {
    let jd0 = 2_451_545.0;
    // Starting exactly on the node, forward must NOT bounce on its
    // own result: it steps a half nodal month to the next node.
    let mut moon = synthetic_moon(27.21, jd0);
    let hit = moon_node_crossing(&mut moon, jd0, SearchDirection::Forward).unwrap();
    assert!((hit.jd - (jd0 + 13.605)).abs() < 1e-6, "got {}", hit.jd);
}

#[test]
fn moon_node_reports_divergence_and_bad_inputs() {
    // Constant positive latitude with a negative rate: Newton walks
    // away from the start past the 30-day window and diverges honestly.
    let mut moon = |_: f64| Ok((100.0, 5.0, -1.0));
    match moon_node_crossing(&mut moon, 2_451_545.0, SearchDirection::Forward).unwrap_err() {
        Error::MoonNodeDiverged { .. } => {}
        other => panic!("expected divergence, got {other:?}"),
    }
    // Constant positive latitude with a positive rate: the direction
    // clamp keeps resetting the walk, so the cap exhausts instead —
    // the same `Maximum iterations` outcome the release raises here.
    let mut moon = |_: f64| Ok((100.0, 5.0, 1.0));
    match moon_node_crossing(&mut moon, 2_451_545.0, SearchDirection::Forward).unwrap_err() {
        Error::MoonNodeNotConverged { .. } => {}
        other => panic!("expected non-convergence, got {other:?}"),
    }
    // Non-finite provider values are errors, never crossings.
    let mut moon = |_: f64| Ok((100.0, f64::NAN, 1.0));
    assert_eq!(
        moon_node_crossing(&mut moon, 2_451_545.0, SearchDirection::Forward).unwrap_err(),
        Error::NonFiniteSearchValue { x: 2_451_545.0 }
    );
    assert!(matches!(
        moon_node_crossing(&mut synthetic_moon(27.21, 0.0), f64::NAN, SearchDirection::Forward)
            .unwrap_err(),
        Error::NonFiniteSearchValue { x } if x.is_nan()
    ));
}

// --: rise/set/transit over analytic providers

#[test]
fn event_kind_precedence_matches_release() {
    use crate::constants::{BIT_CIVIL_TWILIGHT, CALC_ITRANSIT, CALC_MTRANSIT};
    use crate::constants::{CALC_RISE, CALC_SET};
    assert_eq!(resolve_event_kind(CALC_RISE), RiseSetEvent::Rise);
    assert_eq!(resolve_event_kind(CALC_SET), RiseSetEvent::Set);
    assert_eq!(
        resolve_event_kind(CALC_MTRANSIT),
        RiseSetEvent::UpperTransit
    );
    assert_eq!(
        resolve_event_kind(CALC_ITRANSIT),
        RiseSetEvent::LowerTransit
    );
    // SET > RISE > ITRANSIT > MTRANSIT; modifier-only defaults to RISE.
    assert_eq!(resolve_event_kind(CALC_RISE | CALC_SET), RiseSetEvent::Set);
    assert_eq!(
        resolve_event_kind(CALC_MTRANSIT | CALC_ITRANSIT),
        RiseSetEvent::LowerTransit
    );
    assert_eq!(
        resolve_event_kind(CALC_RISE | CALC_MTRANSIT),
        RiseSetEvent::Rise
    );
    assert_eq!(resolve_event_kind(0), RiseSetEvent::Rise);
    assert_eq!(resolve_event_kind(BIT_CIVIL_TWILIGHT), RiseSetEvent::Rise);
}

#[test]
fn geopos_rejects_non_finite_components() {
    assert!(validate_geopos(12.5, 41.9, 0.0).is_ok());
    assert!(matches!(
        validate_geopos(f64::NAN, 41.9, 0.0).unwrap_err(),
        Error::InvalidObserver {
            field: "lon_deg",
            ..
        }
    ));
    assert_eq!(
        validate_geopos(12.5, f64::INFINITY, 0.0).unwrap_err(),
        Error::InvalidObserver {
            field: "lat_deg",
            value: f64::INFINITY,
        }
    );
    assert!(matches!(
        validate_geopos(12.5, 41.9, f64::NEG_INFINITY).unwrap_err(),
        Error::InvalidObserver { field: "alt_m", .. }
    ));
}

#[test]
fn pressure_estimate_pins() {
    assert_eq!(pressure_from_altitude_m(0.0), SEA_LEVEL_PRESSURE_HPA);
    let p1k = pressure_from_altitude_m(1000.0);
    assert!(
        (850.0..910.0).contains(&p1k),
        "1000 m pressure out of band: {p1k}"
    );
    assert!(pressure_from_altitude_m(100.0) > p1k);
    // Above ~44.3 km the base clamps: pressure ~0 instead of NaN.
    assert_eq!(pressure_from_altitude_m(50_000.0), 0.0);
    assert!(pressure_from_altitude_m(50_000.0).is_finite());
}

#[test]
fn twilight_override_rules() {
    use crate::constants::{BIT_ASTRO_TWILIGHT, BIT_CIVIL_TWILIGHT, BIT_NAUTIC_TWILIGHT};
    let civil = twilight_override(BIT_CIVIL_TWILIGHT, true).expect("sun twilight");
    assert_eq!(civil.horizon_alt_deg, CIVIL_TWILIGHT_ALT_DEG);
    assert!(civil.force_no_refraction && civil.force_disc_center);
    let nautic = twilight_override(BIT_NAUTIC_TWILIGHT, true).expect("sun twilight");
    assert_eq!(nautic.horizon_alt_deg, NAUTIC_TWILIGHT_ALT_DEG);
    let astro = twilight_override(BIT_ASTRO_TWILIGHT, true).expect("star twilight");
    assert_eq!(astro.horizon_alt_deg, ASTRO_TWILIGHT_ALT_DEG);
    // Civil has precedence when several bits are set.
    let both = twilight_override(BIT_CIVIL_TWILIGHT | BIT_ASTRO_TWILIGHT, true).unwrap();
    assert_eq!(both.horizon_alt_deg, CIVIL_TWILIGHT_ALT_DEG);
    // Planets and the Moon keep the ordinary horizon crossing.
    assert!(twilight_override(BIT_CIVIL_TWILIGHT, false).is_none());
    assert!(twilight_override(0, true).is_none());
}

#[test]
fn disc_selection_rules() {
    use crate::constants::{BIT_DISC_BOTTOM, BIT_DISC_CENTER};
    assert_eq!(select_disc_point(0, false), DiscPoint::UpperLimb);
    assert_eq!(select_disc_point(BIT_DISC_CENTER, false), DiscPoint::Center);
    assert_eq!(
        select_disc_point(BIT_DISC_BOTTOM, false),
        DiscPoint::LowerLimb
    );
    assert_eq!(
        select_disc_point(BIT_DISC_CENTER | BIT_DISC_BOTTOM, false),
        DiscPoint::Center
    );
    // Fixed stars always use the center, whatever the bits say.
    assert_eq!(select_disc_point(0, true), DiscPoint::Center);
    assert_eq!(select_disc_point(BIT_DISC_BOTTOM, true), DiscPoint::Center);
}

#[test]
fn disc_term_follows_release_asin_rule() {
    let term = disc_altitude_term_deg(DiscPoint::UpperLimb, 0.01, 1.0).unwrap();
    assert!((term - 0.01f64.asin().to_degrees()).abs() < 1e-15);
    let lower = disc_altitude_term_deg(DiscPoint::LowerLimb, 0.01, 1.0).unwrap();
    assert!((lower + term).abs() < 1e-15);
    assert_eq!(
        disc_altitude_term_deg(DiscPoint::Center, 0.01, 1.0).unwrap(),
        0.0
    );
    assert_eq!(
        disc_altitude_term_deg(DiscPoint::UpperLimb, 0.0, 1.0).unwrap(),
        0.0
    );
    // Clamped at the limb: R > d still yields exactly 90 deg, never NaN.
    let clamped = disc_altitude_term_deg(DiscPoint::UpperLimb, 2.0, 1.0).unwrap();
    assert!((clamped - 90.0).abs() < 1e-12);
    // A non-positive distance is an error, never a silent term.
    assert!(disc_altitude_term_deg(DiscPoint::UpperLimb, 0.01, 0.0).is_err());
}

#[test]
fn horizon_sentinel_resolves_past_the_dip() {
    assert_eq!(resolve_horizon_altitude_deg(5.0, -0.5), 5.0);
    assert_eq!(resolve_horizon_altitude_deg(-0.5, -0.5), -0.5);
    assert_eq!(
        resolve_horizon_altitude_deg(HORIZON_DIP_SENTINEL_DEG, -0.5),
        -0.5 + HORIZON_DIP_EPSILON_DEG
    );
    // Composition pieces: twilight overrides the horizon, the disc
    // term adds, refraction replaces the true altitude.
    let h = event_height_deg(0.5, None, 0.0);
    assert_eq!(h, 0.5);
    let h = event_height_deg(0.5, Some(1.0), CIVIL_TWILIGHT_ALT_DEG);
    assert_eq!(h, 7.0);
}

#[test]
fn tret_contract_is_ten_slots_with_jd_first() {
    let tret = make_tret(2_451_545.5);
    assert_eq!(tret.len(), 10);
    assert_eq!(tret[0], 2_451_545.5);
    assert!(tret[1..].iter().all(|&v| v == 0.0));
    // The circumpolar outcome carries 0.0 in slot 0.
    assert_eq!(make_tret(0.0)[0], 0.0);
}

#[test]
fn true_to_apparent_converges_and_clamps_to_dip() {
    // Constant refraction: the fixed point is true + 0.5 exactly.
    let mut r = |_: f64| Ok(0.5);
    let app = true_to_apparent_altitude(10.0, -1.0, &mut r, 0.0).unwrap();
    assert!((app - 10.5).abs() < 1e-9, "got {app}");
    // Refracted point below the dip: left unrefracted.
    let mut r = |_: f64| Ok(0.5);
    let app = true_to_apparent_altitude(-2.0, -1.0, &mut r, 0.0).unwrap();
    assert_eq!(app, -2.0);
    // Null refraction is the identity above the dip.
    let mut r = |_: f64| Ok(0.0);
    let app = true_to_apparent_altitude(3.0, -1.0, &mut r, 0.0).unwrap();
    assert!((app - 3.0).abs() < 1e-12);
    // Non-finite refraction is an error, never a silent altitude.
    let mut bad = |_: f64| Ok(f64::NAN);
    assert!(true_to_apparent_altitude(3.0, -1.0, &mut bad, 7.0).is_err());
    assert!(true_to_apparent_altitude(f64::NAN, -1.0, &mut r, 7.0).is_err());
}

#[test]
fn transit_finds_meridian_crossing() {
    let jd0 = 2_451_545.0;
    // Fixed star: RA = 0 h, GAST advancing 24 h/day -> transit at jd0.
    let mut ra = |_: f64| Ok((0.0, 10.0));
    let mut gast = move |jd: f64| Ok(24.0 * (jd - jd0));
    let sol = calculate_transit(jd0, 0.0, false, &mut ra, &mut gast).unwrap();
    assert!(sol.converged);
    assert!((sol.jd_ut - jd0).abs() < 1e-6, "got {}", sol.jd_ut);
    // Lower transit is half a day later under this GAST rate.
    let mut ra = |_: f64| Ok((0.0, 10.0));
    let mut gast = move |jd: f64| Ok(24.0 * (jd - jd0));
    let sol = calculate_transit(jd0, 0.0, true, &mut ra, &mut gast).unwrap();
    assert!(sol.converged);
    assert!((sol.jd_ut - (jd0 + 0.5)).abs() < 1e-6, "got {}", sol.jd_ut);
    // Searching from 6 h before still lands on jd0 (forward seed).
    let mut ra = |_: f64| Ok((0.0, 10.0));
    let mut gast = move |jd: f64| Ok(24.0 * (jd - jd0));
    let sol = calculate_transit(jd0 - 0.25, 0.0, false, &mut ra, &mut gast).unwrap();
    assert!(sol.converged);
    assert!((sol.jd_ut - jd0).abs() < 1e-6, "got {}", sol.jd_ut);
    // Non-finite GAST is an error, never a transit time.
    let mut ra = |_: f64| Ok((0.0, 10.0));
    let mut bad = |_: f64| Ok(f64::NAN);
    assert!(calculate_transit(jd0, 0.0, false, &mut ra, &mut bad).is_err());
}

/// Analytic diurnal altitude: sine with a 1-day period, ±30 deg.
fn sine_altitude(jd0: f64) -> impl FnMut(f64) -> std::result::Result<f64, Error> {
    move |jd: f64| Ok(30.0 * (2.0 * std::f64::consts::PI * (jd - jd0)).sin())
}

/// Flat declination provider for the circumpolar pre-check.
fn flat_dec(dec: f64) -> impl FnMut(f64) -> std::result::Result<(f64, f64), Error> {
    move |_: f64| Ok((0.0, dec))
}

#[test]
fn rise_set_finds_first_crossing_each_way() {
    let jd0 = 2_451_545.0;
    // Rising through 0 at jd0: start a quarter-day early.
    let out = calculate_rise_set(
        jd0 - 0.25,
        40.0,
        true,
        &mut sine_altitude(jd0),
        &mut flat_dec(0.0),
        0.0,
        0.0,
    )
    .unwrap();
    match out {
        RiseSetOutcome::Found { jd_ut } => {
            assert!((jd_ut - jd0).abs() < 1e-4, "got {jd_ut}")
        }
        RiseSetOutcome::Circumpolar => panic!("expected a rise"),
    }
    // Setting through 0 at jd0 + 0.5: start exactly at the rise.
    let out = calculate_rise_set(
        jd0,
        40.0,
        false,
        &mut sine_altitude(jd0),
        &mut flat_dec(0.0),
        0.0,
        0.0,
    )
    .unwrap();
    match out {
        RiseSetOutcome::Found { jd_ut } => {
            assert!((jd_ut - (jd0 + 0.5)).abs() < 1e-4, "got {jd_ut}")
        }
        RiseSetOutcome::Circumpolar => panic!("expected a set"),
    }
    // Polarity is direction-aware: a set search from below the horizon
    // skips the rise and lands on the next set.
    let out = calculate_rise_set(
        jd0 - 0.25,
        40.0,
        false,
        &mut sine_altitude(jd0),
        &mut flat_dec(0.0),
        0.0,
        0.0,
    )
    .unwrap();
    match out {
        RiseSetOutcome::Found { jd_ut } => {
            assert!((jd_ut - (jd0 + 0.5)).abs() < 1e-4, "got {jd_ut}")
        }
        RiseSetOutcome::Circumpolar => panic!("expected a set"),
    }
}

#[test]
fn circumpolar_verdicts_cover_both_branches() {
    let jd0 = 2_451_545.0;
    // Always above: the pre-check fires (dec 80 at lat 40: min_alt 30).
    let mut high = |_: f64| Ok(10.0);
    let out =
        calculate_rise_set(jd0, 40.0, true, &mut high, &mut flat_dec(80.0), 0.0, 0.0).unwrap();
    assert_eq!(out, RiseSetOutcome::Circumpolar);
    // Always below: max_alt is negative.
    let mut low = |_: f64| Ok(-10.0);
    let out =
        calculate_rise_set(jd0, 40.0, false, &mut low, &mut flat_dec(-80.0), 0.0, 0.0).unwrap();
    assert_eq!(out, RiseSetOutcome::Circumpolar);
    // Not circumpolar by declination, but no crossing in the window:
    // constant +5 deg height with equatorial geometry.
    let mut flat = |_: f64| Ok(5.0);
    let out = calculate_rise_set(jd0, 40.0, true, &mut flat, &mut flat_dec(0.0), 0.0, 0.0).unwrap();
    assert_eq!(out, RiseSetOutcome::Circumpolar);
}

#[test]
fn rise_set_rejects_bad_inputs() {
    let jd0 = 2_451_545.0;
    let mut h = sine_altitude(jd0);
    let mut rd = flat_dec(0.0);
    assert!(calculate_rise_set(f64::NAN, 40.0, true, &mut h, &mut rd, 0.0, 0.0).is_err());
    let mut h = sine_altitude(jd0);
    let mut rd = flat_dec(0.0);
    assert!(calculate_rise_set(jd0, f64::INFINITY, true, &mut h, &mut rd, 0.0, 0.0).is_err());
    // Non-finite heights are errors, never silent circumpolar verdicts.
    let mut bad = |_: f64| Ok(f64::NAN);
    let mut rd = flat_dec(0.0);
    assert!(calculate_rise_set(jd0, 40.0, true, &mut bad, &mut rd, 0.0, 0.0).is_err());
    // A non-finite start is an error for transits too.
    let mut ra = |_: f64| Ok((0.0, 0.0));
    let mut gast = |_: f64| Ok(0.0);
    assert!(calculate_transit(f64::NAN, 0.0, false, &mut ra, &mut gast).is_err());
}

// --: shadow-geometry data and fixed-state helpers ---------------
#[test]
fn safe_math_clamps_verbatim() {
    assert_eq!(safe_acos(2.0), 0.0);
    assert_eq!(safe_acos(-2.0), PI);
    assert_eq!(safe_acos(0.0), PI / 2.0);
    assert_eq!(safe_sqrt(-1e-12), 0.0);
    assert_eq!(safe_sqrt(4.0), 2.0);
}

#[test]
fn obscuration_branches() {
    // No overlap: at and past the tangent separation.
    assert_eq!(calculate_obscuration(0.26, 0.27, 0.53).unwrap(), 0.0);
    assert_eq!(calculate_obscuration(0.26, 0.27, 1.0).unwrap(), 0.0);
    // Concentric total (Moon larger): exactly 1.0.
    assert_eq!(calculate_obscuration(0.26, 0.27, 0.0).unwrap(), 1.0);
    // Concentric annular (Moon smaller): the disc-area ratio.
    let annular = calculate_obscuration(0.27, 0.26, 0.0).unwrap();
    assert_eq!(annular, (0.26 / 0.27_f64).powi(2));
    assert!(annular < 1.0);
    // Sub-lens separation takes the contained-disc arm, not the lens.
    let tiny = calculate_obscuration(0.26, 0.27, 1e-11).unwrap();
    assert_eq!(tiny, 1.0);
    // Genuine partial overlap: strictly inside (0, 1).
    let partial = calculate_obscuration(0.26, 0.27, 0.30).unwrap();
    assert!(partial > 0.0 && partial < 1.0);
    // Non-finite inputs are typed refusals, never NaN passthrough.
    assert!(calculate_obscuration(f64::NAN, 0.27, 0.3).is_err());
    assert!(calculate_obscuration(0.26, f64::INFINITY, 0.3).is_err());
}

#[test]
fn magnitude_branches() {
    // Central: full ratio.
    assert_eq!(calculate_eclipse_magnitude(0.0, 1.0).unwrap(), 1.0);
    // Past the limit: exactly zero.
    assert_eq!(calculate_eclipse_magnitude(1.6, 1.0).unwrap(), 0.0);
    // At the limit: the ramp reaches zero.
    assert_eq!(calculate_eclipse_magnitude(1.55, 1.0).unwrap(), 0.0);
    // Near-miss ramp: ((1.55 - 1.549) / 0.02) * 0.01 * ratio
    // (the 1.55 - 1.549 subtraction is itself inexact in binary).
    let ramp = calculate_eclipse_magnitude(1.549, 2.0).unwrap();
    assert_eq!(ramp, ((1.55 - 1.549) / 0.02) * 0.01 * 2.0);
    // Ceiling clamp.
    assert_eq!(calculate_eclipse_magnitude(0.0, 10.0).unwrap(), 1.5);
    // Custom-limit arm routes to the same core.
    assert_eq!(
        calculate_eclipse_magnitude_with_limit(0.0, 0.5, 1.0).unwrap(),
        0.5
    );
    assert!(calculate_eclipse_magnitude(f64::NAN, 1.0).is_err());
}

#[test]
fn shallow_and_near_miss_predicates() {
    assert!(is_shallow_eclipse(0.009));
    assert!(!is_shallow_eclipse(0.01));
    assert!(is_near_miss_eclipse(1.54));
    assert!(!is_near_miss_eclipse(1.0));
    assert!(is_near_miss_eclipse_with_limit(0.99, 1.0));
    assert!(!is_near_miss_eclipse_with_limit(0.5, 1.0));
}

#[test]
fn contact_time_sentinels_verbatim() {
    let jd_max = 2_460_409.26;
    assert_eq!(
        validate_contact_time(jd_max - 0.1, jd_max, 0.25),
        jd_max - 0.1
    );
    assert_eq!(validate_contact_time(0.0, jd_max, 0.25), 0.0);
    assert_eq!(validate_contact_time(-1.0, jd_max, 0.25), 0.0);
    assert_eq!(validate_contact_time(jd_max + 0.3, jd_max, 0.25), 0.0);
    // NaN passes through: every release comparison is false for NaN.
    assert!(validate_contact_time(f64::NAN, jd_max, 0.25).is_nan());
}

#[test]
fn eclipse_kind_parsing() {
    assert_eq!(parse_eclipse_kind("solar").unwrap(), EclipseKind::Solar);
    assert_eq!(parse_eclipse_kind("SOLAR").unwrap(), EclipseKind::Solar);
    assert_eq!(parse_eclipse_kind("Lunar").unwrap(), EclipseKind::Lunar);
    let err = parse_eclipse_kind("luner").unwrap_err();
    assert!(matches!(err, Error::InvalidEclipseType { .. }));
}

#[test]
fn cycle_rounding_is_bankers() {
    // Rounding uses ties to even at halfway values
    // for 0.5; the saros/inex counts need this spelling.
    assert_eq!(round_half_even(0.5), 0.0);
    assert_eq!(round_half_even(1.5), 2.0);
    assert_eq!(round_half_even(2.5), 2.0);
    assert_eq!(round_half_even(-0.5), 0.0);
}

#[test]
fn saros_tables_resolve_pinned_anchors() {
    // Solar anchor rows recover their exact (series, member).
    assert_eq!(saros_info(2460409.262, EclipseKind::Solar), (139.0, 30.0));
    assert_eq!(saros_number(2460409.262, EclipseKind::Solar), 139);
    // One mean cycle later: same series, next member.
    assert_eq!(
        saros_info(2460409.262 + SAROS_CYCLE_DAYS, EclipseKind::Solar),
        (139.0, 31.0)
    );
    // Lunar anchor row.
    assert_eq!(saros_info(2457117.0, EclipseKind::Lunar), (132.0, 30.0));
    // Far from any grid: the no-match sentinel pair.
    assert_eq!(
        saros_info(2450000.0, EclipseKind::Solar),
        (SAROS_NO_MATCH, SAROS_NO_MATCH)
    );
    assert_eq!(saros_number(2450000.0, EclipseKind::Solar), -99999999);
}

#[test]
fn inex_tables_resolve_pinned_anchors() {
    // Release docstring pins: 2024 solar inex 50, 2022 lunar inex 50.
    assert_eq!(inex_number(2460409.296, EclipseKind::Solar), 50);
    assert_eq!(inex_number(2459891.459, EclipseKind::Lunar), 50);
    // Table shape guards catch a truncated transcription.
    assert_eq!(SOLAR_SAROS_REFERENCES.len(), 93);
    assert_eq!(LUNAR_SAROS_REFERENCES.len(), 101);
    assert_eq!(SOLAR_INEX_REFERENCES.len(), 5);
    assert_eq!(LUNAR_INEX_REFERENCES.len(), 5);
}

#[test]
fn besselian_interpolation_matches_taylor_spelling() {
    let e = BesselianElements::new(
        2_460_409.25,
        0.3145,
        0.2731,
        7.5821,
        0.5436,
        -0.0047,
        359.9,
        0.5123,
        0.1456,
        0.0012,
        -0.0001,
        -0.0001,
        15.0041,
    )
    .unwrap();
    let later = interpolate_besselian_elements(&e, 2_460_409.25 + 0.5 / 24.0);
    assert_eq!(later.t0(), 2_460_409.25 + 0.5 / 24.0);
    // The Taylor step runs over the epoch-derived dt, which is itself
    // inexact in binary ((t0 + 0.5/24) - t0 is not exactly 0.5/24);
    // pin the association, not a hand-rounded decimal.
    let dt = (later.t0() - e.t0()) * 24.0;
    assert_eq!(later.x(), 0.3145 + 0.5123 * dt);
    assert_eq!(later.l2(), -0.0047 + -0.0001 * dt);
    assert_eq!(later.mu(), (359.9 + 15.0041 * dt).rem_euclid(360.0));
    // An exactly representable 6-hour step is fully bitwise.
    let six = interpolate_besselian_elements(&e, 2_460_409.25 + 0.25);
    assert_eq!(six.x(), 0.3145 + 0.5123 * 6.0);
    assert_eq!(six.mu(), (359.9_f64 + 15.0041 * 6.0).rem_euclid(360.0));
    // Rates carry over unchanged.
    assert_eq!(later.dx_dt(), 0.5123);
    assert_eq!(later.dmu_dt(), 15.0041);
    // Non-finite construction words are refused.
    assert!(
        BesselianElements::new(
            f64::NAN,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
        )
        .is_err()
    );
}
#[test]
fn besselian_core_on_axis_geometry() {
    // Moon near Earth on the Sun-Earth line: x = y = 0 by construction.
    // (A Moon parked at 1 AU would push sin_f past 1 — unphysical
    // input, NaN out on both sides; the sweep uses real states.)
    let core = besselian_core_via((1.0, 0.0, 0.0), (0.00257, 0.0, 0.0), 100.0).unwrap();
    assert_eq!(core.x(), 0.0);
    assert_eq!(core.y(), 0.0);
    assert_eq!(core.gamma(), 0.0);
    assert_eq!(core.d_deg(), 0.0);
    assert_eq!(core.mu_deg(), 100.0);
    assert_eq!(core.z_moon_er(), 0.00257 / ECL_REARTH_AU);
    assert!(core.l1() > 0.0);
    // Polar shadow axis: the ex = (1, 0, 0) rescue arm.
    let polar = besselian_core_via((0.0, 0.0, 1.0), (0.0, 0.0, 0.00257), 0.0).unwrap();
    assert_eq!(polar.x(), 0.0);
    assert_eq!(polar.y(), 0.0);
    assert_eq!(polar.d_deg(), 90.0);
    // Non-finite states are refused, never silent NaN cores.
    assert!(besselian_core_via((1.0, 0.0, f64::NAN), (0.9, 0.0, 0.0), 0.0).is_err());
}

#[test]
fn besselian_axis_norm_retains_square_residues() {
    use super::geometry::fundamental_plane_axis_norm as norm;
    // Rounded components of a unit axis. A platform hypot loses the square
    // residues here; source-identical vectors then acquire different x/y.
    let (x, y) = (
        f64::from_bits(0x3fee30d883b6a12f),
        f64::from_bits(0x3fd376e83e5e5e19),
    );
    // Pinned Python 3.12 norm, separately checked with a 180-digit decimal
    // square/sum/sqrt calculation. Python 3.9 has a different hypot word.
    let expected = f64::from_bits(0x3fefb87e11fb9479);
    for exponent in [-1000, -40, 0] {
        let scale = 2.0_f64.powi(exponent);
        for (a, b) in [(x, y), (y, x), (-x, y), (x, -y)] {
            assert_eq!(
                norm(a * scale, b * scale).to_bits(),
                (expected * scale).to_bits()
            );
        }
    }
    // Source-identical Cartesian states and GAST: the public core must use
    // the compensated norm, not retain an unused correct private helper.
    let w = f64::from_bits;
    let core = besselian_core_via(
        (
            w(0x3fee3c78cd6371c2),
            w(0x3fd37de8cbba05bc),
            w(0x3fc0e6a056bc7a63),
        ),
        (
            w(0x3f62a77e02a15495),
            w(0x3f47110df9527856),
            w(0x3f34dd4dac8dbc74),
        ),
        w(0x40583131c87b1836),
    )
    .unwrap();
    assert_eq!(core.x().to_bits(), 0xbfe5831272e72b0d);
    assert_eq!(core.y().to_bits(), 0x3fa086f28f4c26f4);
    assert_eq!(norm(-0.0, 0.0).to_bits(), 0);
    let tiny = f64::from_bits(1);
    assert_eq!(norm(tiny, tiny).to_bits(), 1);
    for bits in [
        0x000f_ffff_ffff_ffff,
        0x0010_0000_0000_0000,
        0x0010_0000_0000_0001,
        0x3fdf_ffff_ffff_ffff,
        0x3fe0_0000_0000_0000,
        0x3fe0_0000_0000_0001,
    ] {
        let value = f64::from_bits(bits);
        assert_eq!(norm(value, 0.0).to_bits(), bits);
    }
}

#[test]
fn besselian_rates_preserve_ordered_divisions() {
    // Analytic cancellation-free fixture, not astronomical golden data.
    // The one-minute quotient must round before conversion to hours;
    // a pre-rounded reciprocal changes this final word by one ULP.
    let numerator = f64::from_bits(0x3fb9999999999992);
    let make = |value| {
        BesselianCore::new(value, value, value, value, value, value, 0.01, 0.01, 1.0).unwrap()
    };
    let rates =
        besselian_rates_via(&make(0.0), &make(numerator), BESSELIAN_RATE_STEP_DAYS).unwrap();
    for value in [
        rates.dx_dt(),
        rates.dy_dt(),
        rates.dd_dt(),
        rates.dmu_dt(),
        rates.dl1_dt(),
        rates.dl2_dt(),
    ] {
        assert_eq!(value.to_bits(), 0x4007fffffffffff8);
    }
    assert!(matches!(
        besselian_rates_via(&make(0.0), &make(1.0), f64::from_bits(1)),
        Err(Error::MathDomain { .. })
    ));
}

#[test]
fn besselian_rates_fold_mu_only() {
    let mk =
        |mu: f64| BesselianCore::new(0.1, 0.2, 7.0, mu, 0.54, -0.005, 0.01, 0.01, 23000.0).unwrap();
    // Identical cores: zero rates.
    let same = besselian_rates_via(&mk(100.0), &mk(100.0), BESSELIAN_RATE_STEP_DAYS).unwrap();
    assert_eq!(same.dx_dt(), 0.0);
    assert_eq!(same.dmu_dt(), 0.0);
    // Forward wrap 359 -> 1 over the 1-minute stencil: +2 deg/h*60.
    let wrapped = besselian_rates_via(&mk(359.0), &mk(1.0), BESSELIAN_RATE_STEP_DAYS).unwrap();
    assert_eq!(wrapped.dmu_dt(), 60.0);
    // Backward wrap folds the other way.
    let back = besselian_rates_via(&mk(1.0), &mk(359.0), BESSELIAN_RATE_STEP_DAYS).unwrap();
    assert_eq!(back.dmu_dt(), -60.0);
    // Non-positive steps are refused.
    assert!(besselian_rates_via(&mk(1.0), &mk(2.0), 0.0).is_err());
    assert!(besselian_rates_via(&mk(1.0), &mk(2.0), f64::NAN).is_err());
}

// --: where/how fixed-state cores over analytic states --------
#[test]
fn vincenty_separation_pins() {
    assert_eq!(
        vincenty_separation_deg(10.0, 20.0, 10.0, 20.0).unwrap(),
        0.0
    );
    assert_eq!(vincenty_separation_deg(0.0, 0.0, 90.0, 0.0).unwrap(), 90.0);
    assert_eq!(
        vincenty_separation_deg(0.0, 0.0, 180.0, 0.0).unwrap(),
        180.0
    );
    // Wrap takes the short way: 359 -> 1 is 2 degrees (Vincenty
    // rounding, not bitwise — pin the association loosely).
    let wrapped = vincenty_separation_deg(359.0, 0.0, 1.0, 0.0).unwrap();
    assert!((wrapped - 2.0).abs() < 1e-9);
    // Pole: longitude differences collapse.
    let polar = vincenty_separation_deg(0.0, 89.0, 180.0, 89.0).unwrap();
    assert!((polar - 2.0).abs() < 1e-9);
    assert!(vincenty_separation_deg(f64::NAN, 0.0, 0.0, 0.0).is_err());
}

#[test]
fn where_core_syzygy_is_central_with_km_dcore() {
    // Moon between Earth and Sun on the x-axis: the shadow axis
    // pierces the globe (Sun-Moon distance ~1 AU keeps the cone
    // words physical; parking the Moon next to the Sun would push
    // sinf past 1 — unphysical input, NaN out on both sides).
    let w =
        eclipse_where_core_via((0.00257, 0.0, 0.0), (1.0, 0.0, 0.0), 100.0, ECL_RSUN_AU).unwrap();
    assert!(w.retc() & crate::constants::ECL_CENTRAL != 0);
    assert!(w.retc() & crate::constants::ECL_PARTIAL == 0);
    assert!(w.lon_deg() > -180.0 && w.lon_deg() <= 180.0);
    assert!(w.lat_deg() >= -90.0 && w.lat_deg() <= 90.0);
    let d = w.dcore();
    assert_eq!(d.len(), 10);
    assert_eq!((d[7], d[8], d[9]), (0.0, 0.0, 0.0));
    // Cone cosines are proper cosines; ground diameters are km-scale.
    assert!(d[5] > 0.0 && d[5] <= 1.0);
    assert!(d[6] > 0.0 && d[6] <= 1.0);
    assert!(d[0].abs() < 30000.0 && d[1] > 0.0);
}

#[test]
fn where_core_clean_miss_is_zero() {
    // Moon far off the Sun-Earth line: axis misses by ~65 radii.
    let w =
        eclipse_where_core_via((0.0, 0.00257, 0.001), (1.0, 0.0, 0.0), 100.0, ECL_RSUN_AU).unwrap();
    assert_eq!(w.retc(), 0);
    assert!(w.lon_deg().is_finite() && w.lat_deg().is_finite());
    assert!(w.dcore()[2] > 6378.0);
    // Non-finite inputs are typed refusals, never silent geometry.
    assert!(
        eclipse_where_core_via(
            (0.0, 0.00257, 0.001),
            (1.0, 0.0, 0.0),
            f64::NAN,
            ECL_RSUN_AU
        )
        .is_err()
    );
}

#[test]
fn sol_how_concentric_annular_contract() {
    // Concentric discs, Moon smaller: annular with the disc-area
    // ratio below 1 and the VISIBLE bit over a high Sun.
    let how = sol_how_core_via(
        (10.0, 20.0, 1.0),
        (10.0, 20.0, 0.00257),
        ECL_RSUN_AU,
        true,
        false,
        (180.0, 45.0, 45.5),
        0.0,
        (139.0, 30.0),
    )
    .unwrap();
    assert_eq!(
        how.retc() & crate::constants::ECL_ANNULAR,
        crate::constants::ECL_ANNULAR
    );
    let attr = how.attr();
    assert_eq!(attr.len(), 20);
    assert!(attr[2] > 0.9 && attr[2] < 1.0);
    assert_eq!(attr[2], (attr[1]).powi(2));
    assert_eq!(attr[7], 0.0);
    assert_eq!((attr[4], attr[5], attr[6]), (180.0, 45.0, 45.5));
    assert_eq!(attr[8], attr[1]);
    assert_eq!((attr[9], attr[10]), (139.0, 30.0));
    assert_eq!(attr[3], 0.0);
    assert!(how.retc() & crate::constants::ECL_VISIBLE != 0);
}

#[test]
fn sol_how_total_exceeds_one_solar_but_caps_occultation() {
    // Closer Moon: larger disc, total, obscuration above 1 on the
    // solar path (disc-area ratio, uncapped).
    let mk = |is_solar: bool, convention: bool| {
        sol_how_core_via(
            (10.0, 20.0, 1.0),
            (10.0, 20.0, 0.00240),
            ECL_RSUN_AU,
            is_solar,
            convention,
            (180.0, 45.0, 45.5),
            0.0,
            (139.0, 30.0),
        )
        .unwrap()
    };
    let total = mk(true, false);
    assert_eq!(
        total.retc() & crate::constants::ECL_TOTAL,
        crate::constants::ECL_TOTAL
    );
    assert!(total.attr()[2] > 1.0);
    // Occultation path caps the covered fraction at 1.0...
    let occ = mk(false, false);
    assert_eq!(
        occ.retc() & crate::constants::ECL_TOTAL,
        crate::constants::ECL_TOTAL
    );
    assert_eq!(occ.attr()[2], 1.0);
    assert_eq!(occ.attr()[8], 0.0);
    assert_eq!(occ.attr()[9], 0.0);
    // ... while the WHERE convention keeps the uncapped ratio.
    let wh = mk(false, true);
    assert_eq!(wh.attr()[2], total.attr()[2]);
}

#[test]
fn sol_how_no_overlap_sentinel_and_degenerate_radius() {
    // Clean miss: retc 0 with the 1.0 obscuration sentinel.
    let miss = sol_how_core_via(
        (10.0, 20.0, 1.0),
        (15.0, 20.0, 0.00257),
        ECL_RSUN_AU,
        true,
        false,
        (180.0, 45.0, 45.5),
        0.0,
        (139.0, 30.0),
    )
    .unwrap();
    assert_eq!(miss.retc(), 0);
    assert_eq!(miss.attr()[2], 1.0);
    // Zero body radius (star): magnitude 1.0, ratio 0, sentinel.
    let star = sol_how_core_via(
        (10.0, 20.0, 1.0),
        (10.0, 20.0, 0.00257),
        0.0,
        false,
        false,
        (180.0, 45.0, 45.5),
        0.0,
        (0.0, 0.0),
    )
    .unwrap();
    assert_eq!(
        (star.attr()[0], star.attr()[1], star.attr()[2]),
        (1.0, 0.0, 1.0)
    );
    // Non-finite words are refused.
    assert!(
        sol_how_core_via(
            (10.0, 20.0, 1.0),
            (10.0, 20.0, f64::NAN),
            ECL_RSUN_AU,
            true,
            false,
            (180.0, 45.0, 45.5),
            0.0,
            (139.0, 30.0),
        )
        .is_err()
    );
}

#[test]
fn sol_how_asin_domain_is_nan_not_raise() {
    // Unphysical radius/distance words (asin past 1): the release
    // `math.asin` raises; Rust yields NaN — but the NaN then fails
    // both `rsun > 0.0` guards, so the ratio branch stays
    // release-identical (0.0), the magnitude branch takes its
    // `else 1.0` arm, and the all-false classification lands on the
    // retc == 0 sentinel arm. Documented divergence, physical
    // states never trigger it.
    let how = sol_how_core_via(
        (10.0, 20.0, 1.0),
        (10.0, 20.0, 0.00257),
        2.0,
        true,
        false,
        (180.0, 45.0, 45.5),
        0.0,
        (139.0, 30.0),
    )
    .unwrap();
    assert_eq!(how.retc(), 0);
    assert_eq!(how.attr()[0], 1.0);
    assert_eq!(how.attr()[1], 0.0);
    assert_eq!(how.attr()[2], 1.0);
}

#[test]
fn sol_compositions_gate_and_zero_verbatim() {
    let w =
        eclipse_where_core_via((0.00257, 0.0, 0.0), (1.0, 0.0, 0.0), 100.0, ECL_RSUN_AU).unwrap();
    let how = sol_how_core_via(
        (10.0, 20.0, 1.0),
        (10.0, 20.0, 0.00240),
        ECL_RSUN_AU,
        true,
        false,
        (180.0, 45.0, 45.5),
        0.0,
        (139.0, 30.0),
    )
    .unwrap();
    // Where: global retc, 10-slot geopos, attr[3] from dcore[0].
    let (ret_where, geopos, attr_where) = sol_eclipse_where_parts(&w, &how);
    assert_eq!(ret_where, w.retc());
    assert_eq!(geopos.len(), 10);
    assert_eq!((geopos[0], geopos[1]), (w.lon_deg(), w.lat_deg()));
    assert!(geopos[2..].iter().all(|v| *v == 0.0));
    assert_eq!(attr_where[3], w.dcore()[0]);
    // How: local retc gains the global central bits.
    let (ret_how, attr_how) = sol_eclipse_how_parts(&how, &w);
    assert_eq!(
        ret_how & crate::constants::ECL_TOTAL,
        crate::constants::ECL_TOTAL
    );
    assert_eq!(
        ret_how & crate::constants::ECL_CENTRAL,
        w.retc() & crate::constants::ECL_CENTRAL
    );
    assert_eq!(attr_how[3], w.dcore()[0]);
    // Below the apparent horizon: ret 0 with slots zeroed, 4..7 kept.
    let low = sol_how_core_via(
        (10.0, 20.0, 1.0),
        (10.0, 20.0, 0.00240),
        ECL_RSUN_AU,
        true,
        false,
        (180.0, -3.0, -0.5),
        0.0,
        (139.0, 30.0),
    )
    .unwrap();
    let (ret_low, attr_low) = sol_eclipse_how_parts(&low, &w);
    assert_eq!(ret_low, 0);
    for i in [0, 1, 2, 3, 8, 9, 10] {
        assert_eq!(attr_low[i], 0.0);
    }
    assert_eq!((attr_low[4], attr_low[5], attr_low[6]), (180.0, -3.0, -0.5));
}

#[test]
fn lun_how_opposition_is_total_quadrature_is_none() {
    // Exact opposition: axis through the selenocenter, total with
    // a plausible umbral magnitude and saros fill.
    let total = lun_how_core_via((-0.00257, 0.0, 0.0), (1.0, 0.0, 0.0), (132.0, 30.0)).unwrap();
    assert_eq!(total.retc(), crate::constants::ECL_TOTAL);
    let attr = total.attr();
    assert!(attr[0] > 1.0 && attr[0] < 2.5);
    assert_eq!(attr[8], attr[0]);
    assert!(attr[1] > attr[0]);
    assert_eq!((attr[9], attr[10]), (132.0, 30.0));
    assert_eq!(attr[2], 0.0);
    assert_eq!(attr[3], 0.0);
    let d = total.dcore();
    assert_eq!(d.len(), 10);
    assert_eq!(d[0], 0.0);
    // AU-based cone cosines sit in slots 3..4 (slots 5..9 reserved).
    assert!((d[3] > 0.0) && (d[3] <= 1.0) && (d[4] > 0.0) && (d[4] <= 1.0));
    assert_eq!((d[5], d[6], d[7], d[8], d[9]), (0.0, 0.0, 0.0, 0.0, 0.0));
    // Quadrature: no eclipse, negative penumbral magnitude, slot 7
    // stays 0.
    let none = lun_how_core_via((0.0, 0.00257, 0.0), (1.0, 0.0, 0.0), (132.0, 30.0)).unwrap();
    assert_eq!(none.retc(), 0);
    assert!(none.attr()[1] < 0.0);
    assert_eq!(none.attr()[7], 0.0);
    assert_eq!(none.attr()[8], 0.0);
    assert!(lun_how_core_via((0.0, f64::NAN, 0.0), (1.0, 0.0, 0.0), (132.0, 30.0)).is_err());
}

#[test]
fn lun_magnitude_gamma_projections_match_pin_legs() {
    // Thin C7–C9 projections over the canonical core: opposition
    // serves the positive core channels; quadrature floors to 0.0.
    let moon = (-0.00257, 0.0, 0.0);
    let sun = (1.0, 0.0, 0.0);
    let saros = (132.0, 30.0);
    let core = lun_how_core_via(moon, sun, saros).unwrap();
    assert_eq!(
        lun_umbral_magnitude_via(moon, sun, saros).unwrap(),
        core.attr()[0].max(0.0)
    );
    assert_eq!(
        lun_penumbral_magnitude_via(moon, sun, saros).unwrap(),
        core.attr()[1].max(0.0)
    );
    assert!(core.attr()[0] > 0.0);
    // Exact opposition sits on the shadow axis: gamma is exactly 0.
    let eps = 23.43928_f64.to_radians();
    assert_eq!(lun_gamma_via(moon, sun, eps).unwrap(), 0.0);
    // Quadrature: no eclipse — magnitudes floor to 0.0, gamma lands
    // deep in the no-eclipse regime (> 1.5 Earth radii).
    let off = (0.0, 0.00257, 0.0);
    assert_eq!(lun_umbral_magnitude_via(off, sun, saros).unwrap(), 0.0);
    assert_eq!(lun_penumbral_magnitude_via(off, sun, saros).unwrap(), 0.0);
    assert!(lun_gamma_via(off, sun, eps).unwrap() > 1.5);
    // Non-finite legs refuse on all three arms.
    assert!(lun_umbral_magnitude_via((0.0, f64::NAN, 0.0), sun, saros).is_err());
    assert!(lun_penumbral_magnitude_via(moon, (f64::INFINITY, 0.0, 0.0), saros).is_err());
    assert!(lun_gamma_via(moon, sun, f64::NAN).is_err());
    assert!(lun_gamma_via((0.0, 0.0, f64::NAN), sun, eps).is_err());
}

#[test]
fn lun_composition_keeps_attrs_below_horizon() {
    let total = lun_how_core_via((-0.00257, 0.0, 0.0), (1.0, 0.0, 0.0), (132.0, 30.0)).unwrap();
    // Moon up: retc preserved with the horizontal triple filled in.
    let (ret_up, attr_up) = lun_eclipse_how_parts(&total, (200.0, 30.0, 30.5));
    assert_eq!(ret_up, crate::constants::ECL_TOTAL);
    assert_eq!((attr_up[4], attr_up[5], attr_up[6]), (200.0, 30.0, 30.5));
    assert_eq!(attr_up[0], total.attr()[0]);
    // Moon down: retc 0, but — unlike the solar path — the
    // attributes stay filled.
    let (ret_down, attr_down) = lun_eclipse_how_parts(&total, (200.0, -30.0, -1.0));
    assert_eq!(ret_down, 0);
    assert_eq!(attr_down[0], total.attr()[0]);
    assert_eq!(attr_down[6], -1.0);
}

// --: global eclipse searches and contact mapping -------------
use std::cell::Cell;

/// Circular coplanar Sun/Moon model with constant speeds: syzygies
/// are analytic, so the finders are pinned against exact answers.
fn circular_state(
    jd: f64,
    sun_lon0: f64,
    moon_lon0: f64,
    wsun: f64,
    wmoon: f64,
    jd0: f64,
) -> SunMoonState {
    SunMoonState::try_new(
        sun_lon0 + wsun * (jd - jd0),
        0.0,
        1.0,
        wsun,
        moon_lon0 + wmoon * (jd - jd0),
        0.0,
        0.00257,
        wmoon,
    )
    .unwrap()
}

#[test]
fn syzygy_finders_hit_analytic_conjunctions() {
    let (jd0, wsun, wmoon) = (2460000.0, 0.9856, 13.1764);
    let rel = wmoon - wsun;
    let mut state = |jd: f64| Ok(circular_state(jd, 100.0, 100.0, wsun, wmoon, jd0));
    // Conjunctions at jd0 + k * 360/rel (linear model: exact).
    let period = 360.0 / rel;
    let next = find_next_new_moon(&mut state, jd0 + 10.0).unwrap();
    assert!((next - (jd0 + 10.0 + (360.0 - rel * 10.0) / rel)).abs() < 1e-6);
    // The fastest lunation check: the returned instant is a
    // conjunction to the release gate.
    let at = state(next).unwrap();
    assert!(fold_180_deg(at.moon_lon_deg - at.sun_lon_deg).abs() < 1e-4);
    // Previous conjunction mirrors (elongation > 0 steps back).
    let prev = find_previous_new_moon(&mut state, jd0 + 10.0).unwrap();
    assert!((prev - jd0).abs() < 1e-6);
    // Opposition is half a synodic period off the conjunction grid.
    let full = find_next_full_moon(&mut state, jd0).unwrap();
    assert!((full - (jd0 + period / 2.0)).abs() < 1e-6);
    let prev_full = find_previous_full_moon(&mut state, jd0 + period).unwrap();
    assert!((prev_full - (jd0 + period / 2.0)).abs() < 1e-6);
    // Non-finite epochs are refused; provider misses propagate.
    assert!(find_next_new_moon(&mut state, f64::NAN).is_err());
    let mut missing = |_: f64| -> Result<SunMoonState> {
        Err(Error::EphemerisRange {
            requested_jd: 0.0,
            start_jd: 0.0,
            end_jd: 1.0,
            body: Some(10),
        })
    };
    assert!(find_next_new_moon(&mut missing, jd0).is_err());
}

#[test]
fn syzygy_speed_floor_never_divides_by_zero() {
    // Equal speeds: the elongation never changes, so the 0.1 floor
    // fires every step and the 20-step estimate returns finite.
    let mut state = |jd: f64| Ok(circular_state(jd, 100.0, 120.0, 1.0, 1.0, 2460000.0));
    let next = find_next_new_moon(&mut state, 2460000.0).unwrap();
    assert!(next.is_finite());
}

#[test]
fn node_distance_pins() {
    assert_eq!(moon_node_distance_deg(100.0, 100.0).unwrap(), 0.0);
    assert_eq!(moon_node_distance_deg(100.0, 280.0).unwrap(), 0.0);
    assert_eq!(moon_node_distance_deg(100.0, 190.0).unwrap(), 90.0);
    assert_eq!(moon_node_distance_deg(350.0, 10.0).unwrap(), 20.0);
    assert!(moon_node_distance_deg(f64::NAN, 0.0).is_err());
}

#[test]
fn golden_section_finds_parabola_minimum() {
    let mut f = |x: f64| Ok((x - 3.0).powi(2) + 5.0);
    let min = golden_section_min(&mut f, 0.0, 10.0, 1e-9, 60).unwrap();
    assert!((min - 3.0).abs() < 1e-6);
    assert!(golden_section_min(&mut f, 5.0, 5.0, 1e-9, 60).is_err());
    let mut failing = |_: f64| -> Result<f64> {
        Err(Error::EphemerisRange {
            requested_jd: 0.0,
            start_jd: 0.0,
            end_jd: 1.0,
            body: Some(10),
        })
    };
    assert!(golden_section_min(&mut failing, 0.0, 1.0, 1e-9, 60).is_err());
}

#[test]
fn golden_traced_matches_untraced_words() {
    // The traced core is the single loop behind `golden_section_min`:
    // same objective words must give the identical minimum word, and
    // the trace must replay the decisions faithfully.
    let mut f = |x: f64| Ok((x - 3.0).powi(2) + 5.0);
    let plain = golden_section_min(&mut f, 0.0, 10.0, 1e-9, 60).unwrap();
    let mut trace = Vec::new();
    let mut g = |x: f64| Ok((x - 3.0).powi(2) + 5.0);
    let traced = golden_section_min_traced(&mut g, 0.0, 10.0, 1e-9, 60, &mut trace).unwrap();
    assert_eq!(plain.to_bits(), traced.to_bits());
    assert!(!trace.is_empty() && trace.len() <= 60);
    for (k, row) in trace.iter().enumerate() {
        assert_eq!(row.level, k);
        assert_eq!(row.arm_a_taken, row.f_a < row.f_b);
        assert!(row.width_after.is_finite() && row.width_after > 0.0);
    }
    for pair in trace.windows(2) {
        assert!(pair[1].width_after < pair[0].width_after);
    }
    assert!(trace.last().unwrap().width_after < 1e-9);
}

#[test]
fn search_flag_helpers_match_release_tables() {
    use crate::constants::{
        ECL_ANNULAR, ECL_ANNULAR_TOTAL, ECL_CENTRAL, ECL_NONCENTRAL, ECL_PARTIAL, ECL_TOTAL,
    };
    assert!(sol_glob_accepts(0, ECL_TOTAL | ECL_CENTRAL));
    assert!(!sol_glob_accepts(ECL_CENTRAL, ECL_TOTAL | ECL_CENTRAL));
    assert!(sol_glob_accepts(ECL_TOTAL, ECL_TOTAL | ECL_CENTRAL));
    assert!(sol_glob_accepts(ECL_TOTAL, ECL_TOTAL | ECL_NONCENTRAL));
    assert!(!sol_glob_accepts(ECL_TOTAL, ECL_ANNULAR | ECL_CENTRAL));
    assert!(!sol_glob_accepts(
        ECL_TOTAL | ECL_ANNULAR,
        ECL_TOTAL | ECL_CENTRAL
    ));
    assert!(sol_glob_accepts(
        ECL_CENTRAL | ECL_TOTAL,
        ECL_TOTAL | ECL_CENTRAL
    ));
    assert!(!sol_glob_accepts(
        ECL_CENTRAL | ECL_TOTAL,
        ECL_TOTAL | ECL_NONCENTRAL
    ));
    assert!(sol_glob_accepts(
        ECL_NONCENTRAL | ECL_PARTIAL,
        ECL_PARTIAL | ECL_NONCENTRAL
    ));
    sol_glob_reject_impossible(0).unwrap();
    sol_glob_reject_impossible(ECL_TOTAL).unwrap();
    sol_glob_reject_impossible(ECL_CENTRAL | ECL_TOTAL).unwrap();
    sol_glob_reject_impossible(ECL_TOTAL | ECL_PARTIAL).unwrap();
    assert!(matches!(
        sol_glob_reject_impossible(ECL_CENTRAL | ECL_PARTIAL).unwrap_err(),
        Error::InvalidEclipseSearch { .. }
    ));
    assert!(matches!(
        sol_glob_reject_impossible(ECL_NONCENTRAL | ECL_ANNULAR_TOTAL).unwrap_err(),
        Error::InvalidEclipseSearch { .. }
    ));
    assert_eq!(
        parse_search_direction("forward").unwrap(),
        EclipseSearchDirection::Forward
    );
    assert_eq!(
        parse_search_direction("backward").unwrap(),
        EclipseSearchDirection::Backward
    );
    assert_eq!(
        parse_search_direction("bidirectional").unwrap(),
        EclipseSearchDirection::Bidirectional
    );
    assert!(parse_search_direction("sideways").is_err());
    assert!(coerce_backwards(BackwardsArg::Bool(true)).unwrap());
    assert!(!coerce_backwards(BackwardsArg::Int(0)).unwrap());
    assert!(coerce_backwards(BackwardsArg::Int(-2)).unwrap());
    assert!(coerce_backwards(BackwardsArg::Text("backward")).unwrap());
    assert!(!coerce_backwards(BackwardsArg::Text(" Forward ")).unwrap());
    assert!(!coerce_backwards(BackwardsArg::Text("")).unwrap());
    assert!(coerce_backwards(BackwardsArg::Text("sideways")).is_err());
    assert_eq!(strip_one_try_bit(ECL_ONE_TRY | 1), 1);
    assert_eq!(strip_one_try_bit(0), 0);
    assert_eq!(strip_one_try_bit(2048), 2048);
}

#[test]
fn gamma_dot_products_retain_addition_residuals() {
    assert_eq!(compensated_sum3([1e16, 1.0, -1e16]), 1.0);
    assert_eq!(compensated_sum3([1.0, 1e16, -1e16]), 1.0);
    // Matched apparent vectors for the 2023 hybrid maximum. Ordinary
    // dot-product summation loses the compensated low word here.
    let sun = (
        f64::from_bits(0x3febf8a4c922ba3f),
        f64::from_bits(0x3fdd0f442305c418),
        f64::from_bits(0x3fc93186c53301cd),
    );
    let moon = (
        f64::from_bits(0x3f61e80c14b6087d),
        f64::from_bits(0x3f52bd9f7ea61891),
        f64::from_bits(0x3f3f4cdcad5e5853),
    );
    assert_eq!(
        gamma_from_vectors(sun, moon).unwrap().to_bits(),
        0.39515565246468454_f64.to_bits(),
    );
}

#[test]
fn gamma_geometry_pins() {
    // Exact syzygy on the axis: zero shadow offset (bitwise).
    assert_eq!(
        gamma_from_vectors((1.0, 0.0, 0.0), (0.0025, 0.0, 0.0)).unwrap(),
        0.0
    );
    // 1e-4 AU off-axis at the Moon: ~2.35 Earth radii (band pin;
    // the release sweep pins the exact words).
    let g = gamma_from_vectors((1.0, 0.0, 0.0), (0.0025, 1e-4, 0.0)).unwrap();
    assert!(g > 2.3 && g < 2.4, "gamma = {g}");
    assert!(gamma_from_vectors((0.0, f64::NAN, 0.0), (1.0, 0.0, 0.0)).is_err());
    assert!(matches!(
        gamma_from_vectors((1.0, 0.0, 0.0), (1.0, 0.0, 0.0)).unwrap_err(),
        Error::ZeroVectorNormalize
    ));
}

#[test]
fn shadow_limit_signs_follow_geometry() {
    // Moon well inside the umbra vertex: negative (total) l2.
    let l2_near = umbra_limit_er(1.0, 0.0024).unwrap();
    assert!(l2_near < 0.0, "l2 = {l2_near}");
    // Moon past the vertex: positive (annular) l2.
    let l2_far = umbra_limit_er(1.0, 0.0027).unwrap();
    assert!(l2_far > 0.0, "l2 = {l2_far}");
    let l1 = penumbra_limit_er(1.0, 0.00257).unwrap();
    assert!(l1 > 0.0 && l1 > l2_near.abs() && l1 > l2_far);
    assert!(penumbra_limit_er(f64::INFINITY, 1.0).is_err());
}

fn solar_test_state() -> SunMoonState {
    SunMoonState::try_new(0.0, 0.0, 1.0, 0.9856, 0.5, 0.3, 0.0025, 13.1764).unwrap()
}

#[test]
fn solar_classification_covers_every_arm() {
    use crate::constants::{ECL_ANNULAR, ECL_ANNULAR_TOTAL, ECL_CENTRAL, ECL_NONCENTRAL};
    use crate::constants::{ECL_PARTIAL, ECL_TOTAL};
    let state = solar_test_state();
    let calls = Cell::new(0);
    let mut never = || -> Result<(i32, f64)> {
        calls.set(calls.get() + 1);
        panic!("where gate must not be consulted past the gamma limit");
    };
    let mut never_hybrid = || -> Result<bool> {
        panic!("hybrid probe must not run without a central core shadow");
    };
    // Beyond the penumbral limit: no eclipse, gate untouched.
    let none =
        classify_solar_eclipse(&state, 2.0, 0.53, -0.02, &mut never, &mut never_hybrid).unwrap();
    assert_eq!((none.eclipse_type, none.magnitude), (0, 0.0));
    assert_eq!(calls.get(), 0);
    assert_eq!(none.gamma, 2.0);
    // Near-miss with a negative where-magnitude: admitted by the
    // spherical test but rejected by the ellipsoidal gate.
    let mut gate_neg = || Ok((0, -0.1));
    let graze_out =
        classify_solar_eclipse(&state, 1.525, 0.53, -0.02, &mut gate_neg, &mut never_hybrid)
            .unwrap();
    assert_eq!((graze_out.eclipse_type, graze_out.magnitude), (0, 0.0));
    // Near-miss with overlap: shallow partial with the ramp magnitude.
    let mut gate_pos = || Ok((ECL_NONCENTRAL, 0.2));
    let graze =
        classify_solar_eclipse(&state, 1.525, 0.53, -0.02, &mut gate_pos, &mut never_hybrid)
            .unwrap();
    assert_eq!(graze.eclipse_type, ECL_PARTIAL | ECL_NONCENTRAL);
    let ratio = graze.moon_sun_ratio;
    let remaining = 1.53 - 1.525;
    let expected = ((remaining / NEAR_MISS_GAMMA_MARGIN) * SHALLOW_ECLIPSE_MAG_THRESHOLD * ratio)
        .clamp(0.0, ECLIPSE_MAGNITUDE_MAX);
    assert!((graze.magnitude - expected).abs() < 1e-12);
    // Central total / annular / hybrid via the where-centrality bit.
    for (l2, hybrid, expected) in [
        (-0.02, false, ECL_TOTAL | ECL_CENTRAL),
        (0.02, false, ECL_ANNULAR | ECL_CENTRAL),
        (0.005, true, ECL_ANNULAR_TOTAL | ECL_CENTRAL),
    ] {
        let mut gate = || Ok((ECL_CENTRAL, 0.5));
        let mut probe = || Ok(hybrid);
        let classified =
            classify_solar_eclipse(&state, 0.3, 0.53, l2, &mut gate, &mut probe).unwrap();
        assert_eq!(classified.eclipse_type, expected, "l2 = {l2}");
    }
    // Non-central core touch keeps the default partial label with
    // the non-central bit; a retc of 0 keeps plain partial.
    let mut gate_nc = || Ok((ECL_NONCENTRAL, 0.4));
    let mut probe = || Ok(false);
    let nc = classify_solar_eclipse(&state, 0.3, 0.53, -0.02, &mut gate_nc, &mut probe).unwrap();
    assert_eq!(nc.eclipse_type, ECL_TOTAL | ECL_NONCENTRAL);
    let mut gate_zero = || Ok((0, 0.4));
    let partial =
        classify_solar_eclipse(&state, 0.3, 0.53, -0.02, &mut gate_zero, &mut probe).unwrap();
    assert_eq!(partial.eclipse_type, ECL_PARTIAL | ECL_NONCENTRAL);
    assert!(
        classify_solar_eclipse(&state, f64::NAN, 0.53, 0.02, &mut gate_zero, &mut probe).is_err()
    );
}

#[test]
fn hybrid_probe_pins_fast_reject_and_sign_change() {
    let jd_max = 2460000.0;
    let calls = Cell::new(0);
    let mut contact = |_: f64, _: f64, _: bool| -> Result<f64> {
        calls.set(calls.get() + 1);
        panic!("fast reject must precede any contact solve");
    };
    let mut width = |_: f64| -> Result<f64> {
        panic!("fast reject must precede any width read");
    };
    assert!(!is_hybrid_solar_eclipse(jd_max, -0.01, &mut contact, &mut width).unwrap());
    assert!(!is_hybrid_solar_eclipse(jd_max, 0.01, &mut contact, &mut width).unwrap());
    assert_eq!(calls.get(), 0);
    // Contacts straddle the maximum; widths change sign along the path.
    let mut contact_ok = |jd: f64, _: f64, before: bool| -> Result<f64> {
        Ok(if before { jd - 0.05 } else { jd + 0.05 })
    };
    let mut width_cross = |jd: f64| -> Result<f64> { Ok((jd - jd_max) * 1.0e6) };
    assert!(is_hybrid_solar_eclipse(jd_max, 0.0, &mut contact_ok, &mut width_cross).unwrap());
    let mut width_pos = |_: f64| -> Result<f64> { Ok(1000.0) };
    assert!(!is_hybrid_solar_eclipse(jd_max, 0.0, &mut contact_ok, &mut width_pos).unwrap());
    // A coverage miss propagates; other read failures skip samples.
    let mut width_boundary = |_: f64| -> Result<f64> {
        Err(Error::EphemerisRange {
            requested_jd: 0.0,
            start_jd: 0.0,
            end_jd: 1.0,
            body: Some(10),
        })
    };
    assert!(is_hybrid_solar_eclipse(jd_max, 0.0, &mut contact_ok, &mut width_boundary).is_err());
    let mut width_skip = |_: f64| -> Result<f64> { Err(Error::NonFiniteSearchValue { x: 0.0 }) };
    assert!(!is_hybrid_solar_eclipse(jd_max, 0.0, &mut contact_ok, &mut width_skip).unwrap());
}

#[test]
fn gamma_contact_finder_pins_crossings_and_gates() {
    let mut gamma = |jd: f64| -> Result<f64> { Ok(0.1 + 50.0 * (jd - 7.0).powi(2)) };
    let before = find_contact_time_gamma(&mut gamma, 7.0, 0.2, true, 0.10).unwrap();
    let after = find_contact_time_gamma(&mut gamma, 7.0, 0.2, false, 0.10).unwrap();
    assert!((before - 6.955277).abs() < 1e-5, "before = {before}");
    assert!((after - 7.044723).abs() < 1e-5, "after = {after}");
    // Both window endpoints below the target: the 0.0 sentinel.
    let mut low = |jd: f64| -> Result<f64> { Ok(0.3 - 0.05 * (jd - 6.9)) };
    assert_eq!(
        find_contact_time_gamma(&mut low, 7.0, 0.35, true, 0.10).unwrap(),
        0.0
    );
}

#[test]
fn sign_bisection_pins_root_and_sentinels() {
    let mut f = |x: f64| -> Result<f64> { Ok(x * x - 2.0) };
    let root = bisect_sign_change(&mut f, 1.0, 2.0).unwrap();
    assert!((root - std::f64::consts::SQRT_2).abs() < 1e-8);
    let mut same_sign = |x: f64| -> Result<f64> { Ok(x * x + 1.0) };
    assert_eq!(bisect_sign_change(&mut same_sign, 0.0, 1.0).unwrap(), 0.0);
    let mut endpoint = |x: f64| -> Result<f64> { Ok(x - 1.0) };
    assert_eq!(bisect_sign_change(&mut endpoint, 1.0, 2.0).unwrap(), 1.0);
    let mut failing = |_: f64| -> Result<f64> {
        Err(Error::EphemerisRange {
            requested_jd: 0.0,
            start_jd: 0.0,
            end_jd: 1.0,
            body: Some(10),
        })
    };
    assert!(bisect_sign_change(&mut failing, 0.0, 1.0).is_err());
}

#[test]
fn besselian_contact_solve_pins_v_contacts() {
    // Half-width 0.2 fits the release 0.35-day expansion window.
    let mut v = |t: f64| -> Result<f64> { Ok((t - 5.0).abs() - 0.2) };
    assert!((solve_besselian_contact(&mut v, 5.0, true, 0.35).unwrap() - 4.8).abs() < 1e-8);
    assert!((solve_besselian_contact(&mut v, 5.0, false, 0.35).unwrap() - 5.2).abs() < 1e-8);
    let mut clear = |_: f64| -> Result<f64> { Ok(1.0) };
    assert_eq!(
        solve_besselian_contact(&mut clear, 5.0, true, 0.35).unwrap(),
        0.0
    );
    // Coverage misses propagate; a caught numerical-domain failure degrades.
    let mut boundary = |_: f64| -> Result<f64> {
        Err(Error::EphemerisRange {
            requested_jd: 0.0,
            start_jd: 0.0,
            end_jd: 1.0,
            body: Some(10),
        })
    };
    assert!(solve_besselian_contact(&mut boundary, 5.0, true, 0.35).is_err());
    let mut degraded = |_: f64| -> Result<f64> { Err(Error::MathDomain { value: -1.0 }) };
    assert_eq!(
        solve_besselian_contact(&mut degraded, 5.0, true, 0.35).unwrap(),
        0.0
    );
}

#[test]
fn besselian_global_range_corruption_and_unclassified_failures_propagate() {
    for primary in [
        Error::EphemerisRange {
            requested_jd: 100.0,
            start_jd: 0.0,
            end_jd: 99.0,
            body: Some(10),
        },
        Error::LebCorrupt {
            what: "contact fixture".into(),
            offset: 7,
        },
        Error::SpkCorrupt {
            what: "contact fixture".into(),
            offset: 9,
        },
    ] {
        let marker = primary;
        let mut residual = |_: f64| Err(marker.clone());
        assert_eq!(
            solve_besselian_contact(&mut residual, 100.0, true, 0.35).unwrap_err(),
            marker
        );
        let mut triple = |_: f64| Err(marker.clone());
        let mut x = |_: f64| Ok(0.0);
        assert_eq!(
            calculate_solar_phases_besselian(100.0, ECL_TOTAL, &mut triple, &mut x).unwrap_err(),
            marker
        );
        let mut triple = v_triple(100.0, 0.06, 0.02, 0.01);
        let mut x = |_: f64| Err(marker.clone());
        assert_eq!(
            calculate_solar_phases_besselian(100.0, ECL_TOTAL, &mut triple, &mut x).unwrap_err(),
            marker
        );
    }
    for marker in [
        Error::ConfigurationMissing {
            missing_config: "event-test",
            suggestion: "configure it",
        },
        Error::NonFiniteSearchValue { x: 100.0 },
    ] {
        let mut residual = |_: f64| Err(marker.clone());
        assert_eq!(
            solve_besselian_contact(&mut residual, 100.0, false, 0.35).unwrap_err(),
            marker
        );
    }
}

#[test]
fn besselian_contact_catches_at_initial_expansion_and_bisection_reads() {
    // A caught expansion sample is skipped; the next clearing sample still
    // brackets the same contact. Before and after walks use the same policy.
    for before in [true, false] {
        let mut calls = Vec::new();
        let mut residual = |jd: f64| {
            calls.push(jd);
            if calls.len() == 2 {
                Err(Error::MathDomain { value: -1.0 })
            } else {
                Ok((jd - 100.0).abs() - 0.03)
            }
        };
        let contact = solve_besselian_contact(&mut residual, 100.0, before, 0.35).unwrap();
        let expected = if before { 99.97 } else { 100.03 };
        assert!((contact - expected).abs() < 1e-8);
        assert_eq!(calls[1], if before { 99.98 } else { 100.02 });
        // This case clears at the first expansion point; read five is the
        // first bisection midpoint, where a caught failure returns the sentinel.
        let mut count = 0;
        let mut residual = |jd: f64| {
            count += 1;
            if count == 5 {
                Err(Error::MathDomain { value: -1.0 })
            } else {
                Ok((jd - 100.0).abs() - 0.01)
            }
        };
        assert_eq!(
            solve_besselian_contact(&mut residual, 100.0, before, 0.35).unwrap(),
            0.0
        );
        assert_eq!(count, 5);
    }
}

#[test]
fn besselian_residual_memo_retains_only_successful_triples() {
    let mut reads = 0;
    let mut residual = |_: f64| {
        reads += 1;
        if reads == 1 {
            Err(Error::MathDomain { value: -1.0 })
        } else {
            Ok((1.0, 2.0, 3.0))
        }
    };
    {
        let mut memo = MemoTriple {
            inner: &mut residual,
            memo: Vec::new(),
        };
        assert_eq!(memo.slot(0, 100.0).unwrap(), None);
        assert!(memo.memo.is_empty());
        assert_eq!(memo.slot(1, 100.0).unwrap(), Some(2.0));
        assert_eq!(memo.slot(2, 100.0).unwrap(), Some(3.0));
        assert_eq!(memo.memo.len(), 1);
    }
    assert_eq!(reads, 2);
}

/// V-shaped synthetic residuals around `tmax` with pen/umbra/axis
/// half-widths, mimicking a central total eclipse.
fn v_triple(
    tmax: f64,
    pen_half: f64,
    umb_half: f64,
    cen_half: f64,
) -> impl FnMut(f64) -> Result<(f64, f64, f64)> {
    move |jd: f64| {
        let d = (jd - tmax).abs();
        Ok((d - pen_half, d - umb_half, d - cen_half))
    }
}

#[test]
fn solar_phases_map_every_slot() {
    use crate::constants::{ECL_CENTRAL, ECL_NONCENTRAL, ECL_PARTIAL, ECL_TOTAL};
    let tmax = 100.0;
    let mut triple = v_triple(tmax, 0.06, 0.02, 0.01);
    let mut x = |jd: f64| -> Result<f64> { Ok(jd - tmax) };
    let tret = calculate_solar_phases_besselian(tmax, ECL_TOTAL | ECL_CENTRAL, &mut triple, &mut x)
        .unwrap();
    assert_eq!(tret[0], tmax);
    assert!((tret[1] - tmax).abs() < 1e-6, "noon = {}", tret[1]);
    assert!((tret[2] - (tmax - 0.06)).abs() < 1e-6);
    assert!((tret[3] - (tmax + 0.06)).abs() < 1e-6);
    assert!((tret[4] - (tmax - 0.02)).abs() < 1e-6);
    assert!((tret[5] - (tmax + 0.02)).abs() < 1e-6);
    assert!((tret[6] - (tmax - 0.01)).abs() < 1e-6);
    assert!((tret[7] - (tmax + 0.01)).abs() < 1e-6);
    assert_eq!((tret[8], tret[9]), (0.0, 0.0));
    // Non-central total: umbral contacts stay, center line goes.
    let mut triple = v_triple(tmax, 0.06, 0.02, 0.01);
    let mut x = |jd: f64| -> Result<f64> { Ok(jd - tmax) };
    let tret_nc =
        calculate_solar_phases_besselian(tmax, ECL_TOTAL | ECL_NONCENTRAL, &mut triple, &mut x)
            .unwrap();
    assert!(tret_nc[4] != 0.0 && tret_nc[5] != 0.0);
    assert_eq!((tret_nc[6], tret_nc[7]), (0.0, 0.0));
    // Partial without penumbral touch: the collapse rule fills
    // tret[4]/[5] with the maximum.
    let mut clear = |_: f64| -> Result<(f64, f64, f64)> { Ok((1.0, 1.0, 1.0)) };
    let mut x = |_: f64| -> Result<f64> { Ok(1.0) };
    let tret_collapsed =
        calculate_solar_phases_besselian(tmax, ECL_PARTIAL | ECL_NONCENTRAL, &mut clear, &mut x)
            .unwrap();
    assert_eq!((tret_collapsed[2], tret_collapsed[3]), (0.0, 0.0));
    assert_eq!((tret_collapsed[4], tret_collapsed[5]), (tmax, tmax));
    // Partial with a resolved penumbra: no collapse, no umbral slots.
    let mut triple = v_triple(tmax, 0.06, 0.02, 0.01);
    let mut x = |jd: f64| -> Result<f64> { Ok(jd - tmax) };
    let tret_p =
        calculate_solar_phases_besselian(tmax, ECL_PARTIAL | ECL_NONCENTRAL, &mut triple, &mut x)
            .unwrap();
    assert!(tret_p[2] != 0.0 && tret_p[3] != 0.0);
    assert_eq!(
        (tret_p[4], tret_p[5], tret_p[6], tret_p[7]),
        (0.0, 0.0, 0.0, 0.0)
    );
}

#[test]
fn lunar_classification_covers_every_arm() {
    use crate::constants::{ECL_PARTIAL, ECL_PENUMBRAL, ECL_TOTAL};
    let sun = (0.0, 0.0, 1.0);
    let mk = |moon_lat: f64| {
        SunMoonState::try_new(
            sun.0, sun.1, sun.2, 0.9856, 180.0, moon_lat, 0.00257, 13.1764,
        )
        .unwrap()
    };
    // Exact opposition: deep total.
    let total = classify_lunar_eclipse(&mk(0.0)).unwrap();
    assert_eq!(total.eclipse_type, ECL_TOTAL);
    assert!(total.magnitude_umbral >= 1.0);
    // One degree off-axis: penumbral-only.
    let pen = classify_lunar_eclipse(&mk(1.0)).unwrap();
    assert_eq!(pen.eclipse_type, ECL_PENUMBRAL);
    assert_eq!(pen.magnitude_umbral, 0.0);
    assert!(pen.magnitude_penumbral > 0.0);
    // Half a degree: partial umbral.
    let part = classify_lunar_eclipse(&mk(0.5)).unwrap();
    assert_eq!(part.eclipse_type, ECL_PARTIAL);
    assert!(part.magnitude_umbral > 0.0 && part.magnitude_umbral < 1.0);
    // Two degrees: clean miss.
    let none = classify_lunar_eclipse(&mk(2.0)).unwrap();
    assert_eq!(none.eclipse_type, 0);
    assert_eq!(
        (none.magnitude_umbral, none.magnitude_penumbral),
        (0.0, 0.0)
    );
    // Shallow penumbral: 0 < pen < 0.01 with the umbral slot zeroed.
    // The target distance is derived from the release constants so
    // the test pins the arm, not a hand constant.
    let pen_r = {
        let earth_sd = (LUNAR_EARTH_RADIUS_AU / 0.00257).atan().to_degrees();
        let sun_sd = (SUN_SEMIDIAMETER_ARCSEC / 3600.0) / 1.0;
        (earth_sd + sun_sd) * LUNAR_SHADOW_ENLARGEMENT
    };
    let moon_sd = (MOON_SEMIDIAMETER_ARCSEC / 3600.0) * (MOON_SEMIDIAMETER_REF_AU / 0.00257);
    let shallow_pen = classify_lunar_eclipse(&mk(pen_r + moon_sd - 0.005 * 2.0 * moon_sd)).unwrap();
    assert_eq!(shallow_pen.eclipse_type, ECL_PENUMBRAL);
    assert_eq!(shallow_pen.magnitude_umbral, 0.0);
    assert!(shallow_pen.magnitude_penumbral > 0.0 && shallow_pen.magnitude_penumbral < 0.01);
    // Degenerate near-zero lunar disc (finite but below the 1e-10
    // arm): classified no-eclipse with geometry.
    let flat = SunMoonState::try_new(0.0, 0.0, 1.0, 0.9856, 180.0, 0.0, 1.0e7, 13.1764).unwrap();
    let flat_class = classify_lunar_eclipse(&flat).unwrap();
    assert_eq!(flat_class.eclipse_type, 0);
    assert_eq!(flat_class.gamma, 0.0);
}

#[test]
fn selenocentric_depth_pins() {
    // Exact opposition on the axis: deep immersion (negative).
    let deep = selenocentric_depth_deg((-0.00257, 0.0, 0.0), (1.0, 0.0, 0.0)).unwrap();
    assert!(deep < -1.2 && deep > -1.3, "depth = {deep}");
    // Quadrature: far outside the shadow (positive).
    let out = selenocentric_depth_deg((0.0, 0.00257, 0.0), (1.0, 0.0, 0.0)).unwrap();
    assert!(out > 80.0, "depth = {out}");
    assert!(selenocentric_depth_deg((f64::NAN, 0.0, 0.0), (1.0, 0.0, 0.0)).is_err());
    assert!(matches!(
        selenocentric_depth_deg((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)).unwrap_err(),
        Error::ZeroVectorNormalize
    ));
}

/// Synthetic canonical shadow: fixed diameters/cosines with a
/// V-shaped axis offset around `tmax`, giving analytic contacts.
fn scripted_shadow(
    tmax: f64,
    rate_au_day: f64,
    retc: i32,
) -> impl FnMut(f64) -> Result<(i32, LunarShadowDcore)> {
    move |jd: f64| {
        Ok((
            retc,
            LunarShadowDcore::try_new(
                rate_au_day * (jd - tmax).abs(),
                4.324e-5,
                7.676e-5,
                1.0,
                1.0,
            )
            .unwrap(),
        ))
    }
}

#[test]
fn lunar_phase_times_map_every_slot() {
    use crate::constants::{ECL_PARTIAL, ECL_PENUMBRAL, ECL_TOTAL};
    let tmax = 200.0;
    // Analytic contacts from the scripted diameters and the exact
    // lunar radius (C/rate at rate 1e-3 AU/day).
    let c_pen = (7.676e-5 / 2.0 + ECL_RMOON_AU) / 1e-3;
    let c_par = (4.324e-5 / 2.0 + ECL_RMOON_AU) / 1e-3;
    let c_tot = (4.324e-5 / 2.0 - ECL_RMOON_AU) / 1e-3;
    let mut shadow = scripted_shadow(tmax, 1e-3, ECL_TOTAL);
    let tret = lun_eclipse_phase_times_via(tmax, ECL_TOTAL, &mut shadow).unwrap();
    assert_eq!((tret[0], tret[1], tret[8], tret[9]), (tmax, 0.0, 0.0, 0.0));
    assert!((tret[6] - (tmax - c_pen)).abs() < 2e-6, "P1 = {}", tret[6]);
    assert!((tret[7] - (tmax + c_pen)).abs() < 2e-6, "P4 = {}", tret[7]);
    assert!((tret[2] - (tmax - c_par)).abs() < 2e-6, "U1 = {}", tret[2]);
    assert!((tret[3] - (tmax + c_par)).abs() < 2e-6, "U4 = {}", tret[3]);
    assert!((tret[4] - (tmax - c_tot)).abs() < 2e-6, "U2 = {}", tret[4]);
    assert!((tret[5] - (tmax + c_tot)).abs() < 2e-6, "U3 = {}", tret[5]);
    // Penumbral-only: umbral/total slots stay zero.
    let mut shadow_pen = scripted_shadow(tmax, 1e-3, ECL_PENUMBRAL);
    let tret_pen = lun_eclipse_phase_times_via(tmax, ECL_PENUMBRAL, &mut shadow_pen).unwrap();
    assert!(tret_pen[6] != 0.0 && tret_pen[7] != 0.0);
    assert_eq!(
        (tret_pen[2], tret_pen[3], tret_pen[4], tret_pen[5]),
        (0.0, 0.0, 0.0, 0.0)
    );
    // Partial: total slots stay zero.
    let mut shadow_part = scripted_shadow(tmax, 1e-3, ECL_PARTIAL);
    let tret_part = lun_eclipse_phase_times_via(tmax, ECL_PARTIAL, &mut shadow_part).unwrap();
    assert!(tret_part[2] != 0.0 && tret_part[3] != 0.0);
    assert_eq!((tret_part[4], tret_part[5]), (0.0, 0.0));
    // Canonical contacts route through the same slots.
    let mut shadow = scripted_shadow(tmax, 1e-3, ECL_TOTAL);
    assert!((lunar_penumbral_contact_p1_via(tmax, &mut shadow).unwrap() - tret[6]).abs() < 1e-12);
    let mut shadow = scripted_shadow(tmax, 1e-3, ECL_TOTAL);
    assert!((lunar_penumbral_contact_p4_via(tmax, &mut shadow).unwrap() - tret[7]).abs() < 1e-12);
    let mut shadow = scripted_shadow(tmax, 1e-3, ECL_TOTAL);
    assert!((lunar_umbral_contact_u1_via(tmax, &mut shadow).unwrap() - tret[2]).abs() < 1e-12);
    let mut shadow = scripted_shadow(tmax, 1e-3, ECL_TOTAL);
    assert!((lunar_umbral_contact_u2_via(tmax, &mut shadow).unwrap() - tret[4]).abs() < 1e-12);
    let mut shadow = scripted_shadow(tmax, 1e-3, ECL_TOTAL);
    assert!((lunar_umbral_contact_u3_via(tmax, &mut shadow).unwrap() - tret[5]).abs() < 1e-12);
    let mut shadow = scripted_shadow(tmax, 1e-3, ECL_TOTAL);
    assert!((lunar_umbral_contact_u4_via(tmax, &mut shadow).unwrap() - tret[3]).abs() < 1e-12);
    // No eclipse at maximum: every contact is 0.0.
    let mut shadow_none = scripted_shadow(tmax, 1e-3, 0);
    assert_eq!(
        lunar_umbral_contact_u1_via(tmax, &mut shadow_none).unwrap(),
        0.0
    );
    // Out-of-range slot index is a typed refusal.
    let mut shadow = scripted_shadow(tmax, 1e-3, ECL_TOTAL);
    assert!(matches!(
        lun_contact_via_core(tmax, 10, &mut shadow).unwrap_err(),
        Error::InvalidEclipseSearch { .. }
    ));
    // Durations compose from the canonical contacts.
    let mut shadow = scripted_shadow(tmax, 1e-3, ECL_TOTAL);
    let total_min = lunar_total_duration_min_via(tmax, &mut shadow).unwrap();
    assert!(
        (total_min - 2.0 * c_tot * 24.0 * 60.0).abs() < 1e-3,
        "total = {total_min}"
    );
    let mut shadow = scripted_shadow(tmax, 1e-3, ECL_TOTAL);
    let umbral_min = lunar_umbral_duration_min_via(tmax, &mut shadow).unwrap();
    assert!(
        (umbral_min - 2.0 * c_par * 24.0 * 60.0).abs() < 1e-3,
        "umbral = {umbral_min}"
    );
    let mut shadow_pen = scripted_shadow(tmax, 1e-3, ECL_PENUMBRAL);
    assert_eq!(
        lunar_total_duration_min_via(tmax, &mut shadow_pen).unwrap(),
        0.0
    );
    assert_eq!(eclipse_phase_duration_min(0.0, 5.0), 0.0);
}

/// Scripted solar ephemeris: circular longitudes (genuine syzygy
/// finding) with eclipses scripted at absolute maxima `maxima`
/// (node prefilter, gamma parabola, central where-gate and V-shaped
/// residuals all keyed to the same instants). States outside
/// `[cover_lo, cover_hi]` are coverage misses, so exhaustion is
/// exercised genuinely.
struct ScriptedSolar {
    sun_lon0: f64,
    moon_lon0: f64,
    wsun: f64,
    wmoon: f64,
    jd0: f64,
    maxima: Vec<f64>,
    cover_lo: f64,
    cover_hi: f64,
    state_calls: Cell<usize>,
}

impl ScriptedSolar {
    fn moon_lon(&self, jd: f64) -> f64 {
        self.moon_lon0 + self.wmoon * (jd - self.jd0)
    }
    fn sun_lon(&self, jd: f64) -> f64 {
        self.sun_lon0 + self.wsun * (jd - self.jd0)
    }
    fn near_max(&self, jd: f64, window: f64) -> bool {
        self.maxima.iter().any(|m| (jd - m).abs() < window)
    }
    fn check_cover(&self, jd: f64) -> Result<()> {
        if jd < self.cover_lo || jd > self.cover_hi {
            return Err(Error::EphemerisRange {
                requested_jd: jd,
                start_jd: self.cover_lo,
                end_jd: self.cover_hi,
                body: Some(10),
            });
        }
        Ok(())
    }
}

impl SolarGlobalEphemeris for ScriptedSolar {
    fn sun_moon_state(&self, jd_ut: f64) -> Result<SunMoonState> {
        self.check_cover(jd_ut)?;
        self.state_calls.set(self.state_calls.get() + 1);
        SunMoonState::try_new(
            self.sun_lon(jd_ut),
            0.0,
            1.0,
            self.wsun,
            self.moon_lon(jd_ut),
            0.0,
            0.0025,
            self.wmoon,
        )
    }
    fn mean_node_lon_deg(&self, jd_ut: f64) -> Result<f64> {
        self.check_cover(jd_ut)?;
        // Near a scripted maximum the node sits 5 degrees off the
        // Moon (prefilter passes); elsewhere it sits 25+ degrees off
        // (prefilter skips, exercising the walk advance).
        Ok(if self.near_max(jd_ut, 16.0) {
            self.moon_lon(jd_ut) - 5.0
        } else {
            self.moon_lon(jd_ut) + 205.0
        })
    }
    fn gamma_er(&self, jd_ut: f64) -> Result<f64> {
        self.check_cover(jd_ut)?;
        let d = self
            .maxima
            .iter()
            .map(|m| (jd_ut - m).abs())
            .fold(f64::INFINITY, f64::min);
        Ok(0.2 + 200.0 * d * d)
    }
    fn shadow_limits_er(&self, jd_ut: f64) -> Result<(f64, f64)> {
        self.check_cover(jd_ut)?;
        Ok((0.53, -0.02))
    }
    fn where_gate(&self, jd_ut: f64) -> Result<(i32, f64)> {
        use crate::constants::ECL_CENTRAL;
        self.check_cover(jd_ut)?;
        Ok(if self.near_max(jd_ut, 1.0) {
            (ECL_CENTRAL, 0.5)
        } else {
            (0, -1.0)
        })
    }
    fn contact_residuals_km(&self, jd_ut: f64) -> Result<(f64, f64, f64)> {
        self.check_cover(jd_ut)?;
        let tmax = *self
            .maxima
            .iter()
            .min_by(|a, b| {
                (jd_ut - **a)
                    .abs()
                    .partial_cmp(&(jd_ut - **b).abs())
                    .unwrap()
            })
            .unwrap();
        let d = (jd_ut - tmax).abs() * 1.0e6;
        Ok((d - 0.06e6, d - 0.02e6, d - 0.01e6))
    }
    fn besselian_x_er(&self, jd_ut: f64) -> Result<f64> {
        self.check_cover(jd_ut)?;
        let tmax = *self
            .maxima
            .iter()
            .min_by(|a, b| {
                (jd_ut - **a)
                    .abs()
                    .partial_cmp(&(jd_ut - **b).abs())
                    .unwrap()
            })
            .unwrap();
        Ok(jd_ut - tmax)
    }
    fn ground_core_width_km(&self, jd_ut: f64) -> Result<f64> {
        self.check_cover(jd_ut)?;
        Ok(1000.0)
    }
}

#[test]
fn solar_walk_finds_filters_and_advances() {
    use crate::constants::{ECL_ANNULAR, ECL_CENTRAL, ECL_PARTIAL, ECL_TOTAL};
    // Circular model: conjunctions at jd0 + k * period.
    let (jd0, wsun, wmoon) = (2460000.0, 0.9856, 13.1764);
    let period = 360.0 / (wmoon - wsun);
    let t0 = jd0 + 2.0 * period + 0.03;
    let t1 = jd0 + 8.0 * period + 0.03;
    let mk = |maxima: Vec<f64>| ScriptedSolar {
        sun_lon0: 100.0,
        moon_lon0: 100.0,
        wsun,
        wmoon,
        jd0,
        maxima,
        cover_lo: jd0 - 400.0,
        cover_hi: jd0 + 400.0,
        state_calls: Cell::new(0),
    };
    // Forward search finds the first scripted maximum with the
    // full slot layout (golden-minimum + V contacts, all analytic).
    let eph = mk(vec![t0, t1]);
    let (retflag, tret) =
        sol_eclipse_when_glob_via(&eph, jd0 - 100.0, 0, EclipseSearchDirection::Forward).unwrap();
    assert_eq!(retflag, ECL_TOTAL | ECL_CENTRAL);
    assert!((tret[0] - t0).abs() < 1e-6, "tmax = {}", tret[0]);
    assert!((tret[2] - (t0 - 0.06)).abs() < 1e-6);
    assert!((tret[3] - (t0 + 0.06)).abs() < 1e-6);
    assert!((tret[4] - (t0 - 0.02)).abs() < 1e-6);
    assert!((tret[5] - (t0 + 0.02)).abs() < 1e-6);
    assert!((tret[1] - t0).abs() < 1e-6);
    // Restarting from the returned maximum advances (epoch margin).
    let (_, tret_next) =
        sol_eclipse_when_glob_via(&eph, tret[0], 0, EclipseSearchDirection::Forward).unwrap();
    assert!((tret_next[0] - t1).abs() < 1e-6);
    // Backward search mirrors.
    let (_, tret_back) =
        sol_eclipse_when_glob_via(&eph, t1 + 100.0, 0, EclipseSearchDirection::Backward).unwrap();
    assert!((tret_back[0] - t1).abs() < 1e-6);
    let (_, tret_prev) =
        sol_eclipse_when_glob_via(&eph, t1, 0, EclipseSearchDirection::Backward).unwrap();
    assert!((tret_prev[0] - t0).abs() < 1e-6);
    // Bidirectional search just after a maximum still finds the next.
    let (_, tret_bi) =
        sol_eclipse_when_glob_via(&eph, t0 + 5.0, 0, EclipseSearchDirection::Bidirectional)
            .unwrap();
    assert!((tret_bi[0] - t1).abs() < 1e-6);
    // A non-matching mask walks to the coverage edge: exhausted.
    let eph = mk(vec![t0, t1]);
    assert!(matches!(
        sol_eclipse_when_glob_via(
            &eph,
            jd0 - 100.0,
            ECL_ANNULAR,
            EclipseSearchDirection::Forward
        )
        .unwrap_err(),
        Error::EclipseSearchExhausted { .. }
    ));
    // An impossible mask is refused before any provider call.
    let eph = mk(vec![t0]);
    assert!(matches!(
        sol_eclipse_when_glob_via(
            &eph,
            jd0,
            ECL_CENTRAL | ECL_PARTIAL,
            EclipseSearchDirection::Forward
        )
        .unwrap_err(),
        Error::InvalidEclipseSearch { .. }
    ));
    assert_eq!(eph.state_calls.get(), 0);
}

/// Scripted lunar ephemeris: circular equatorial vectors (genuine
/// selenocentric-depth minima at oppositions) with the canonical
/// shadow scripted TOTAL near `maxima` and absent elsewhere.
struct ScriptedLunar {
    sun_lon0: f64,
    wsun: f64,
    wmoon: f64,
    jd0: f64,
    maxima: Vec<f64>,
    cover_lo: f64,
    cover_hi: f64,
}

impl ScriptedLunar {
    fn sun_lon(&self, jd: f64) -> f64 {
        (self.sun_lon0 + self.wsun * (jd - self.jd0)).to_radians()
    }
    fn moon_lon(&self, jd: f64) -> f64 {
        (self.sun_lon0 + 180.0 + self.wmoon * (jd - self.jd0)).to_radians()
            % (2.0 * std::f64::consts::PI)
    }
    fn near_max(&self, jd: f64, window: f64) -> bool {
        self.maxima.iter().any(|m| (jd - m).abs() < window)
    }
    fn check_cover(&self, jd: f64) -> Result<()> {
        if jd < self.cover_lo || jd > self.cover_hi {
            return Err(Error::EphemerisRange {
                requested_jd: jd,
                start_jd: self.cover_lo,
                end_jd: self.cover_hi,
                body: Some(10),
            });
        }
        Ok(())
    }
}

impl LunarGlobalEphemeris for ScriptedLunar {
    fn sun_moon_state(&self, jd_ut: f64) -> Result<SunMoonState> {
        self.check_cover(jd_ut)?;
        let sun_lon = self.sun_lon(jd_ut).to_degrees().rem_euclid(360.0);
        let moon_lon =
            (sun_lon + 180.0 + (self.wmoon - self.wsun) * (jd_ut - self.jd0)).rem_euclid(360.0);
        SunMoonState::try_new(
            sun_lon, 0.0, 1.0, self.wsun, moon_lon, 0.0, 0.00257, self.wmoon,
        )
    }
    fn mean_node_lon_deg(&self, jd_ut: f64) -> Result<f64> {
        self.check_cover(jd_ut)?;
        let sun_lon = self.sun_lon(jd_ut).to_degrees().rem_euclid(360.0);
        let moon_lon =
            (sun_lon + 180.0 + (self.wmoon - self.wsun) * (jd_ut - self.jd0)).rem_euclid(360.0);
        Ok(if self.near_max(jd_ut, 16.0) {
            moon_lon - 5.0
        } else {
            moon_lon + 205.0
        })
    }
    fn geocentric_xyz_au(&self, jd_ut: f64) -> Result<MoonSunXyzAu> {
        self.check_cover(jd_ut)?;
        let (sl, ml) = (self.sun_lon(jd_ut), self.moon_lon(jd_ut));
        Ok((
            (0.00257 * ml.cos(), 0.00257 * ml.sin(), 0.0),
            (sl.cos(), sl.sin(), 0.0),
        ))
    }
    fn lun_shadow(&self, jd_ut: f64) -> Result<(i32, LunarShadowDcore)> {
        use crate::constants::ECL_TOTAL;
        self.check_cover(jd_ut)?;
        if self.near_max(jd_ut, 1.0) {
            let tmax = self.maxima[0];
            Ok((
                ECL_TOTAL,
                LunarShadowDcore::try_new(
                    1e-3 * (jd_ut - tmax).abs(),
                    4.324e-5,
                    7.676e-5,
                    1.0,
                    1.0,
                )
                .unwrap(),
            ))
        } else {
            Ok((
                0,
                LunarShadowDcore::try_new(1.0, 4.324e-5, 7.676e-5, 1.0, 1.0).unwrap(),
            ))
        }
    }
}

#[test]
fn lunar_walk_finds_filters_and_advances() {
    use crate::constants::{ECL_PENUMBRAL, ECL_TOTAL};
    let (jd0, wsun, wmoon) = (2460000.0, 0.9856, 13.1764);
    let period = 360.0 / (wmoon - wsun);
    // Opposition grid: jd0 + k * period (the scripted state puts
    // elongation 180 deg at jd0 exactly).
    let u0 = jd0 + 2.0 * period;
    let mk = || ScriptedLunar {
        sun_lon0: 100.0,
        wsun,
        wmoon,
        jd0,
        maxima: vec![u0],
        cover_lo: jd0 - 400.0,
        cover_hi: jd0 + 400.0,
    };
    let eph = mk();
    let (retflag, tret) = lun_eclipse_when_via(&eph, jd0 - 100.0, 0, false).unwrap();
    assert_eq!(retflag, ECL_TOTAL);
    assert!((tret[0] - u0).abs() < 1e-6, "tmax = {}", tret[0]);
    assert!((tret[6] - (u0 - 0.05)).abs() < 1e-5, "P1 = {}", tret[6]);
    assert!((tret[2] - (u0 - 0.03324)).abs() < 1e-5, "U1 = {}", tret[2]);
    assert!((tret[4] - (u0 - 0.01)).abs() < 1e-5, "U2 = {}", tret[4]);
    // Restarting from the maximum advances past it: exhausted (only
    // one scripted eclipse inside coverage).
    assert!(matches!(
        lun_eclipse_when_via(&eph, tret[0], 0, false).unwrap_err(),
        Error::EclipseSearchExhausted { .. }
    ));
    // Backward search from after the maximum finds it.
    let (_, tret_back) = lun_eclipse_when_via(&eph, u0 + 100.0, 0, true).unwrap();
    assert!((tret_back[0] - u0).abs() < 1e-6);
    // A non-matching mask walks to the edge: exhausted.
    assert!(matches!(
        lun_eclipse_when_via(&eph, jd0 - 100.0, ECL_PENUMBRAL, false).unwrap_err(),
        Error::EclipseSearchExhausted { .. }
    ));
    // Annular-only lunar requests are refused up front.
    use crate::constants::ECL_ANNULAR;
    assert!(matches!(
        lun_eclipse_when_via(&eph, jd0, ECL_ANNULAR, false).unwrap_err(),
        Error::InvalidEclipseSearch { .. }
    ));
}

// ----: local eclipse search/visibility and maxima ----

/// Parabolic scripted separation: minimum `sep_min` at `jd_max`.
struct ParabolicSep {
    jd_max: f64,
    sep_min: f64,
    curve: f64,
}

impl ParabolicSep {
    fn sep(&self, jd: f64) -> Result<f64> {
        Ok(self.sep_min + self.curve * (jd - self.jd_max).powi(2))
    }
}

#[test]
fn local_max_time_finds_parabolic_minimum() {
    let model = ParabolicSep {
        jd_max: 2460409.28,
        sep_min: 0.05,
        curve: 200.0,
    };
    let mut sep = |jd: f64| model.sep(jd);
    let (jd_max, min_sep) = calc_local_eclipse_max_time_via(&mut sep, 2460409.30, 0.2).unwrap();
    assert!((jd_max - model.jd_max).abs() < 1e-6, "jd_max = {jd_max}");
    assert!((min_sep - 0.05).abs() < 1e-9, "min_sep = {min_sep}");
    // Bad window is refused, not silently searched.
    assert!(matches!(
        calc_local_eclipse_max_time_via(&mut sep, model.jd_max, 0.0).unwrap_err(),
        Error::InvalidSearchStep { .. }
    ));
}

#[test]
fn find_local_maximum_converges_in_served_window() {
    let model = ParabolicSep {
        jd_max: 2460409.28,
        sep_min: 0.0,
        curve: 500.0,
    };
    let mut sep = |jd: f64| model.sep(jd);
    let found =
        find_local_maximum_via(&mut sep, model.jd_max + 0.05, LOCAL_MAX_SEARCH_RANGE_DAYS).unwrap();
    assert!((found - model.jd_max).abs() < 1e-6, "found = {found}");
}

#[test]
fn local_contact_time_locates_linear_crossing() {
    // Separation falling through 0.5 deg at jd 2460000.1.
    let mut falling = |jd: f64| -> Result<f64> { Ok(1.0 - 5.0 * (jd - 2460000.0)) };
    let found = find_local_contact_time(&mut falling, 2460000.0, 2460000.2, 0.5, false).unwrap();
    assert!((found - 2460000.1).abs() < 1e-7, "found = {found}");
    let mut rising = |jd: f64| -> Result<f64> { Ok(5.0 * (jd - 2460000.0)) };
    let found = find_local_contact_time(&mut rising, 2460000.0, 2460000.2, 0.5, true).unwrap();
    assert!((found - 2460000.1).abs() < 1e-7, "found = {found}");
    // Unachievable target (constant separation above it) is the 0.0
    // sentinel on both arms, not a fabricated contact.
    let mut flat = |_: f64| -> Result<f64> { Ok(2.0) };
    assert_eq!(
        find_local_contact_time(&mut flat, 2460000.0, 2460000.2, 0.5, false).unwrap(),
        0.0
    );
    assert_eq!(
        find_local_contact_time(&mut flat, 2460000.0, 2460000.2, 0.5, true).unwrap(),
        0.0
    );
    // Inverted bracket is refused.
    assert!(matches!(
        find_local_contact_time(&mut flat, 2460000.2, 2460000.0, 0.5, false).unwrap_err(),
        Error::InvalidSearchInterval { .. }
    ));
}

#[test]
fn local_phases_partial_layout_and_ordering() {
    // Central-minimum parabola (sep 0.1 at max) with the Sun up and
    // mean distances: partial geometry, outer contacts only.
    let model = ParabolicSep {
        jd_max: 2460409.28,
        sep_min: 0.1,
        curve: 200.0,
    };
    let mut sep = |jd: f64| model.sep(jd);
    let mut altaz = |_: f64| -> Result<(f64, f64)> { Ok((30.0, 180.0)) };
    let mut dists = |_: f64| -> Result<(f64, f64)> { Ok((1.0, 0.002569)) };
    let phases =
        calculate_local_eclipse_phases_via(model.jd_max, &mut sep, &mut altaz, &mut dists).unwrap();
    let s = phases.as_array();
    assert!((s[0] - model.jd_max).abs() < 1e-6, "max = {}", s[0]);
    assert!(s[1] < s[0] && s[0] < s[4], "outer order: {s:?}");
    assert_eq!(s[2], 0.0, "no second contact for partial: {s:?}");
    assert_eq!(s[3], 0.0, "no third contact for partial: {s:?}");
    // Magnitude = (sum - sep) / sun_diameter with mean radii.
    let sun_r = (959.63 / 3600.0) / 1.0;
    let moon_r = (932.56 / 3600.0) * (0.002569 / 0.002569);
    let expected_mag = (sun_r + moon_r - 0.1) / (2.0 * sun_r);
    assert!((s[5] - expected_mag).abs() < 1e-9, "mag = {}", s[5]);
    assert!((s[6] - expected_mag).abs() < 1e-9, "fraction = {}", s[6]);
    assert!(s[7] > 0.0 && s[7] < 1.0, "obscuration = {}", s[7]);
    assert_eq!((s[8], s[9]), (180.0, 30.0));
}

#[test]
fn local_phases_zero_gates() {
    let model = ParabolicSep {
        jd_max: 2460409.28,
        sep_min: 0.0,
        curve: 200.0,
    };
    // Sun below the -1.0 deg gate at the global maximum.
    let mut sep = |jd: f64| model.sep(jd);
    let mut night = |_: f64| -> Result<(f64, f64)> { Ok((-5.0, 180.0)) };
    let mut dists = |_: f64| -> Result<(f64, f64)> { Ok((1.0, 0.002569)) };
    assert_eq!(
        calculate_local_eclipse_phases_via(model.jd_max, &mut sep, &mut night, &mut dists).unwrap(),
        LocalEclipsePhases::zeros()
    );
    // No disc overlap at a high Sun.
    let far = ParabolicSep {
        jd_max: model.jd_max,
        sep_min: 2.0,
        curve: 200.0,
    };
    let mut sep_far = |jd: f64| far.sep(jd);
    let mut day = |_: f64| -> Result<(f64, f64)> { Ok((30.0, 180.0)) };
    assert_eq!(
        calculate_local_eclipse_phases_via(model.jd_max, &mut sep_far, &mut day, &mut dists)
            .unwrap(),
        LocalEclipsePhases::zeros()
    );
}

#[test]
fn local_phases_central_total_has_inner_contacts() {
    // Exact central minimum with a closer (larger) Moon: total
    // geometry, so second/third contacts resolve inside first/fourth.
    let model = ParabolicSep {
        jd_max: 2460409.28,
        sep_min: 0.0,
        curve: 200.0,
    };
    let mut sep = |jd: f64| model.sep(jd);
    let mut altaz = |_: f64| -> Result<(f64, f64)> { Ok((45.0, 170.0)) };
    let mut dists = |_: f64| -> Result<(f64, f64)> { Ok((1.0, 0.0023)) };
    let s = calculate_local_eclipse_phases_via(model.jd_max, &mut sep, &mut altaz, &mut dists)
        .unwrap()
        .as_array();
    assert!(
        s[1] < s[2] && s[2] < s[0] && s[0] < s[3] && s[3] < s[4],
        "nesting: {s:?}"
    );
    assert!(s[5] >= 1.0, "total magnitude = {}", s[5]);
    assert_eq!(s[7], 1.0, "total obscuration = {}", s[7]);
}

#[test]
fn determine_local_type_branches() {
    use crate::constants::{
        ECL_ANNULAR, ECL_ANNULAR_TOTAL, ECL_CENTRAL, ECL_PARTIAL, ECL_TOTAL, ECL_VISIBLE,
    };
    assert_eq!(determine_local_eclipse_type(0.0, 1.05, true), 0);
    assert_eq!(
        determine_local_eclipse_type(0.5, 0.95, false),
        ECL_VISIBLE | ECL_PARTIAL
    );
    assert_eq!(
        determine_local_eclipse_type(1.1, 1.05, true),
        ECL_VISIBLE | ECL_CENTRAL | ECL_TOTAL
    );
    assert_eq!(
        determine_local_eclipse_type(0.98, 0.95, true),
        ECL_VISIBLE | ECL_CENTRAL | ECL_ANNULAR
    );
    // Near-unity Moon (0.99 < ratio < 1.0) is the hybrid flag arm.
    assert_eq!(
        determine_local_eclipse_type(1.0, 0.995, true),
        ECL_VISIBLE | ECL_CENTRAL | ECL_ANNULAR_TOTAL
    );
}

#[test]
fn magnitude_and_obscuration_at_loc_cores() {
    // Partial instant: sep 0.2 with mean distances.
    let mag = sol_eclipse_magnitude_at_loc_via(30.0, 0.2, 1.0, 0.002569).unwrap();
    let sun_r = (959.63 / 3600.0) / 1.0;
    let moon_r = 932.56 / 3600.0;
    assert!(
        (mag - (sun_r + moon_r - 0.2) / (2.0 * sun_r)).abs() < 1e-12,
        "mag = {mag}"
    );
    let obs = sol_eclipse_obscuration_at_loc_via(30.0, 0.2, 1.0, 0.002569).unwrap();
    assert!(obs > 0.0 && obs < mag, "obs = {obs}");
    // Gates: Sun down and no overlap both read 0.0.
    assert_eq!(
        sol_eclipse_magnitude_at_loc_via(-5.0, 0.2, 1.0, 0.002569).unwrap(),
        0.0
    );
    assert_eq!(
        sol_eclipse_magnitude_at_loc_via(30.0, 2.0, 1.0, 0.002569).unwrap(),
        0.0
    );
    assert_eq!(
        sol_eclipse_obscuration_at_loc_via(-5.0, 0.2, 1.0, 0.002569).unwrap(),
        0.0
    );
    assert_eq!(
        sol_eclipse_obscuration_at_loc_via(30.0, 2.0, 1.0, 0.002569).unwrap(),
        0.0
    );
    // Deep total instant with a close Moon: full area cover.
    assert_eq!(
        sol_eclipse_obscuration_at_loc_via(45.0, 0.001, 1.0, 0.0023).unwrap(),
        1.0
    );
}

/// Scripted served solar provider: one partial eclipse at `jd_e`
/// (parabolic separation, Sun up, circumpolar horizon, fixed
/// `how`/`where` words).
struct ScriptedSolarLoc {
    jd_e: f64,
    sep_min: f64,
    curve: f64,
    sun_up: bool,
}

impl ScriptedSolarLoc {
    fn boundary(&self, jd: f64) -> Error {
        Error::EphemerisRange {
            requested_jd: jd,
            start_jd: self.jd_e - 100.0,
            end_jd: self.jd_e + 100.0,
            body: None,
        }
    }
}

impl SolarLocalEphemeris for ScriptedSolarLoc {
    fn next_global_forward(&self, jd: f64) -> Result<(i32, [f64; 10])> {
        use crate::constants::ECL_PARTIAL;
        if jd <= self.jd_e {
            let mut tret = [0.0; 10];
            tret[0] = self.jd_e;
            Ok((ECL_PARTIAL, tret))
        } else {
            Err(self.boundary(jd))
        }
    }

    fn separation_deg(&self, jd: f64) -> Result<f64> {
        Ok(self.sep_min + self.curve * (jd - self.jd_e).powi(2))
    }

    fn radii_sep_deg(&self, jd: f64, inner_moon: bool) -> Result<(f64, f64, f64)> {
        let rsun = 0.2666;
        let mut rmoon = 0.2590;
        if inner_moon {
            rmoon *= 0.9991;
        }
        Ok((rsun, rmoon, self.separation_deg(jd)?))
    }

    fn sun_altaz_deg(&self, _jd: f64) -> Result<(f64, f64, f64)> {
        if self.sun_up {
            Ok((30.0, 180.0, 30.0))
        } else {
            Ok((-10.0, 180.0, -10.0))
        }
    }

    fn topo_diameters_deg(&self, _jd: f64) -> Result<(f64, f64)> {
        Ok((2.0 * 0.2590, 2.0 * 0.2666))
    }

    fn sun_horizon_crossings(&self, _jd_from: f64) -> Result<BodyHorizonCrossings> {
        Ok(BodyHorizonCrossings {
            rise_rc: -2,
            rise_jd: 0.0,
            set_rc: -2,
            set_jd: 0.0,
        })
    }

    fn how_core(&self, _jd: f64) -> Result<(i32, [f64; 20])> {
        use crate::constants::{ECL_PARTIAL, ECL_VISIBLE};
        Ok((ECL_PARTIAL | ECL_VISIBLE, [0.0; 20]))
    }

    fn where_core(&self, _jd: f64) -> Result<(i32, [f64; 10])> {
        use crate::constants::ECL_NONCENTRAL;
        let mut dcore = [0.0; 10];
        dcore[0] = 50.0;
        Ok((ECL_NONCENTRAL, dcore))
    }
}

#[test]
fn solar_loc_walk_finds_partial_and_advances() {
    use crate::constants::{
        ECL_1ST_VISIBLE, ECL_4TH_VISIBLE, ECL_MAX_VISIBLE, ECL_NONCENTRAL, ECL_PARTIAL, ECL_VISIBLE,
    };
    let eph = ScriptedSolarLoc {
        jd_e: 2460409.28,
        sep_min: 0.1,
        curve: 200.0,
        sun_up: true,
    };
    let (retflag, tret, attr) = sol_eclipse_when_loc_via(&eph, 2460409.28 - 10.0, false).unwrap();
    assert_eq!(
        retflag & (ECL_PARTIAL | ECL_VISIBLE | ECL_MAX_VISIBLE | ECL_1ST_VISIBLE | ECL_4TH_VISIBLE),
        ECL_PARTIAL | ECL_VISIBLE | ECL_MAX_VISIBLE | ECL_1ST_VISIBLE | ECL_4TH_VISIBLE,
    );
    assert_eq!(retflag & ECL_NONCENTRAL, ECL_NONCENTRAL);
    assert!((tret[0] - eph.jd_e).abs() < 1e-6, "max = {}", tret[0]);
    assert!(tret[1] < tret[0] && tret[0] < tret[4], "contacts: {tret:?}");
    assert_eq!(attr[3], 50.0);
    // Restarting from the reported maximum advances past the only
    // scripted eclipse: exhausted, never the same event twice.
    assert!(matches!(
        sol_eclipse_when_loc_via(&eph, tret[0], false).unwrap_err(),
        Error::EclipseSearchExhausted { .. }
    ));
    // Backward search from after the eclipse finds it.
    let (_, tret_back, _) = sol_eclipse_when_loc_via(&eph, eph.jd_e + 100.0, true).unwrap();
    assert!((tret_back[0] - eph.jd_e).abs() < 1e-6);
}

#[test]
fn solar_loc_walk_skips_night_lunation() {
    let eph = ScriptedSolarLoc {
        jd_e: 2460409.28,
        sep_min: 0.1,
        curve: 200.0,
        sun_up: false,
    };
    assert!(matches!(
        sol_eclipse_when_loc_via(&eph, 2460409.28 - 10.0, false).unwrap_err(),
        Error::EclipseSearchExhausted { .. }
    ));
}

/// Scripted served lunar provider: one total eclipse at `jd_max`
/// with penumbral contacts ±0.05 d and umbral contacts ±0.02 d.
struct ScriptedLunarLoc {
    jd_max: f64,
    moon_up: bool,
}

impl ScriptedLunarLoc {
    fn boundary(&self, jd: f64) -> Error {
        Error::EphemerisRange {
            requested_jd: jd,
            start_jd: self.jd_max - 100.0,
            end_jd: self.jd_max + 100.0,
            body: None,
        }
    }

    fn global_tret(&self) -> [f64; 10] {
        let mut tret = [0.0; 10];
        tret[0] = self.jd_max;
        tret[2] = self.jd_max - 0.02;
        tret[3] = self.jd_max + 0.02;
        tret[4] = self.jd_max - 0.01;
        tret[5] = self.jd_max + 0.01;
        tret[6] = self.jd_max - 0.05;
        tret[7] = self.jd_max + 0.05;
        tret
    }
}

impl LunarLocalEphemeris for ScriptedLunarLoc {
    fn next_global(&self, jd: f64, backwards: bool) -> Result<(i32, [f64; 10])> {
        use crate::constants::ECL_TOTAL;
        let hit = if backwards {
            jd >= self.jd_max
        } else {
            jd <= self.jd_max
        };
        if hit {
            Ok((ECL_TOTAL, self.global_tret()))
        } else {
            Err(self.boundary(jd))
        }
    }

    fn how_attr(&self, _jd: f64) -> Result<(i32, [f64; 20])> {
        use crate::constants::ECL_TOTAL;
        let mut attr = [0.0; 20];
        attr[6] = if self.moon_up { 30.0 } else { -10.0 };
        Ok((ECL_TOTAL, attr))
    }

    fn moon_horizon_crossings(&self, _jd_from: f64) -> Result<BodyHorizonCrossings> {
        Ok(BodyHorizonCrossings {
            rise_rc: -2,
            rise_jd: 0.0,
            set_rc: -2,
            set_jd: 0.0,
        })
    }
}

#[test]
fn lunar_loc_walk_finds_total_and_advances() {
    use crate::constants::{
        ECL_MAX_VISIBLE, ECL_PARTBEG_VISIBLE, ECL_PARTEND_VISIBLE, ECL_PENUMBBEG_VISIBLE,
        ECL_PENUMBEND_VISIBLE, ECL_TOTAL, ECL_TOTBEG_VISIBLE, ECL_TOTEND_VISIBLE, ECL_VISIBLE,
    };
    let eph = ScriptedLunarLoc {
        jd_max: 2460000.5,
        moon_up: true,
    };
    let (retflag, tret, _) = lun_eclipse_when_loc_via(&eph, 2460000.5 - 10.0, false).unwrap();
    assert_eq!(retflag & ECL_TOTAL, ECL_TOTAL);
    assert_eq!(retflag & ECL_VISIBLE, ECL_VISIBLE);
    for bit in [
        ECL_MAX_VISIBLE,
        ECL_PARTBEG_VISIBLE,
        ECL_PARTEND_VISIBLE,
        ECL_TOTBEG_VISIBLE,
        ECL_TOTEND_VISIBLE,
        ECL_PENUMBBEG_VISIBLE,
        ECL_PENUMBEND_VISIBLE,
    ] {
        assert_eq!(retflag & bit, bit, "missing bit {bit}: {retflag:#x}");
    }
    assert_eq!(tret[0], eph.jd_max);
    assert_eq!((tret[6], tret[7]), (eph.jd_max - 0.05, eph.jd_max + 0.05));
    // Re-entrancy advances past the only scripted eclipse.
    assert!(matches!(
        lun_eclipse_when_loc_via(&eph, tret[0], false).unwrap_err(),
        Error::EclipseSearchExhausted { .. }
    ));
    // Backward search from after the eclipse finds it.
    let (_, tret_back, _) = lun_eclipse_when_loc_via(&eph, eph.jd_max + 10.0, true).unwrap();
    assert_eq!(tret_back[0], eph.jd_max);
}

#[test]
fn lunar_loc_walk_skips_moon_down() {
    let eph = ScriptedLunarLoc {
        jd_max: 2460000.5,
        moon_up: false,
    };
    assert!(matches!(
        lun_eclipse_when_loc_via(&eph, 2460000.5 - 10.0, false).unwrap_err(),
        Error::EclipseSearchExhausted { .. }
    ));
}

#[test]
fn ecl_eph_flags_keeps_only_ephemeris_bits() {
    use crate::constants::{FLG_JPLEPH, FLG_MOSEPH, FLG_SWIEPH, FLG_TRUEPOS};
    assert_eq!(
        ecl_eph_flags(FLG_SWIEPH | FLG_TRUEPOS | 0x1234_0000_u32 as i32),
        FLG_SWIEPH
    );
    assert_eq!(
        ecl_eph_flags(FLG_JPLEPH | FLG_MOSEPH),
        FLG_JPLEPH | FLG_MOSEPH
    );
}

// -- analytic pins -------------------------------------------

#[test]
fn occ_body_radius_table_identities() {
    use crate::constants::{CERES, MOON, SUN, VENUS};
    // Star names are point sources.
    assert_eq!(occ_body_radius_au(None), 0.0);
    assert_eq!(occ_body_radius_au(Some(SUN)), ECL_RSUN_AU);
    assert_eq!(occ_body_radius_au(Some(MOON)), ECL_RMOON_AU);
    // Classical table entry: Venus 6051.8 km over 1 AU in km.
    assert_eq!(occ_body_radius_au(Some(VENUS)), 6051.8 / ECL_AU_KM);
    // Extra table entry: Ceres 469.7 km.
    assert_eq!(occ_body_radius_au(Some(CERES)), 469.7 / ECL_AU_KM);
    // Unknown ids are point sources.
    assert_eq!(occ_body_radius_au(Some(9999)), 0.0);
}

#[test]
fn planet_angular_radius_spelling_and_fallback() {
    use crate::constants::VENUS;
    //  asin(min(1, r/d)) spelling at 1 AU.
    let expected = (6051.8 / ECL_AU_KM).min(1.0).asin().to_degrees();
    assert_eq!(calc_planet_angular_radius_deg(VENUS, 1.0), expected);
    // Unknown bodies and non-positive distances: point-source fallback.
    assert_eq!(
        calc_planet_angular_radius_deg(9999, 1.0),
        OCC_POINT_SOURCE_RADIUS_DEG
    );
    assert_eq!(
        calc_planet_angular_radius_deg(VENUS, 0.0),
        OCC_POINT_SOURCE_RADIUS_DEG
    );
    assert_eq!(
        calc_planet_angular_radius_deg(VENUS, -2.0),
        OCC_POINT_SOURCE_RADIUS_DEG
    );
    // The extra-table bodies stay point sources on this path.
    assert_eq!(
        calc_planet_angular_radius_deg(crate::constants::CERES, 1.0),
        OCC_POINT_SOURCE_RADIUS_DEG
    );
}

#[test]
fn moon_self_and_pseudo_body_guards() {
    use crate::constants::{MARS, SUN};
    assert!(matches!(
        reject_moon_self_occultation(Some(MOON)).unwrap_err(),
        Error::InvalidEclipseSearch {
            entry: "lun_occult",
            ..
        }
    ));
    assert_eq!(
        reject_moon_self_occultation(Some(MOON))
            .unwrap_err()
            .to_string(),
        "lun_occult: lunar occultation of the Moon itself is undefined (body 1)."
    );
    reject_moon_self_occultation(Some(MARS)).unwrap();
    reject_moon_self_occultation(None).unwrap();
    assert_eq!(fold_pseudo_body_to_sun(-1), SUN);
    assert_eq!(fold_pseudo_body_to_sun(-1000), SUN);
    assert_eq!(fold_pseudo_body_to_sun(MARS), MARS);
    assert_eq!(fold_pseudo_body_to_sun(SUN), SUN);
}

#[test]
fn occultation_filter_expansion_and_refusals() {
    use crate::constants::{
        ECL_ANNULAR, ECL_ANNULAR_TOTAL, ECL_CENTRAL, ECL_NONCENTRAL, ECL_PARTIAL, ECL_TOTAL,
    };
    // Empty filter: everything possible; solar adds the annular bits.
    assert_eq!(
        normalize_occultation_filter(0, false, "3").unwrap(),
        ECL_TOTAL | ECL_PARTIAL | ECL_NONCENTRAL | ECL_CENTRAL
    );
    assert_eq!(
        normalize_occultation_filter(0, true, "0").unwrap(),
        ECL_TOTAL | ECL_PARTIAL | ECL_NONCENTRAL | ECL_CENTRAL | ECL_ANNULAR | ECL_ANNULAR_TOTAL
    );
    // Central partial cannot exist.
    assert!(matches!(
        normalize_occultation_filter(ECL_PARTIAL | ECL_CENTRAL, true, "0").unwrap_err(),
        Error::InvalidEclipseSearch { .. }
    ));
    // Annular-only for a non-solar body is refused with the body word.
    let err = normalize_occultation_filter(ECL_ANNULAR, false, "Regulus").unwrap_err();
    assert_eq!(
        err.to_string(),
        "lun_occult_when_glob: annular occultations do not exist for object Regulus"
    );
    // Combined annular bits are stripped for non-solar bodies.
    assert_eq!(
        normalize_occultation_filter(ECL_TOTAL | ECL_ANNULAR, false, "3").unwrap(),
        ECL_TOTAL | ECL_NONCENTRAL | ECL_CENTRAL
    );
    // Umbral classes imply both central variants; partial implies
    // non-central.
    assert_eq!(
        normalize_occultation_filter(ECL_TOTAL, false, "3").unwrap(),
        ECL_TOTAL | ECL_NONCENTRAL | ECL_CENTRAL
    );
    assert_eq!(
        normalize_occultation_filter(ECL_PARTIAL, false, "3").unwrap(),
        ECL_PARTIAL | ECL_NONCENTRAL
    );
}

#[test]
fn occultation_latitude_gates() {
    assert!(star_never_occulted(7.1));
    assert!(star_never_occulted(-8.0));
    assert!(!star_never_occulted(7.0));
    assert!(!star_never_occulted(0.0));
    assert!(conjunction_latitude_admits_occultation(1.0, -1.0));
    assert!(conjunction_latitude_admits_occultation(2.0, 0.0));
    assert!(!conjunction_latitude_admits_occultation(2.1, 0.0));
    for nonfinite in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(!conjunction_latitude_admits_occultation(nonfinite, 0.0));
        assert!(!conjunction_latitude_admits_occultation(0.0, nonfinite));
    }
}

#[test]
fn occultation_radii_star_and_finite_paths() {
    // Point source: body radius word is exactly 0.0.
    let (rb, rm) = occultation_radii_deg(0.0, 1.5, ECL_RMOON_AU, 0.0026).unwrap();
    assert_eq!(rb, 0.0);
    assert_eq!(rm, (ECL_RMOON_AU / 0.0026).asin().to_degrees());
    // Finite body spelling.
    let (rb2, _) = occultation_radii_deg(ECL_RSUN_AU, 1.0, ECL_RMOON_AU, 0.0026).unwrap();
    assert_eq!(rb2, (ECL_RSUN_AU / 1.0).asin().to_degrees());
    // Non-finite words and non-positive distances are refused.
    assert!(matches!(
        occultation_radii_deg(f64::NAN, 1.0, ECL_RMOON_AU, 0.0026).unwrap_err(),
        Error::NonFiniteArgument { .. }
    ));
    assert!(matches!(
        occultation_radii_deg(ECL_RSUN_AU, 0.0, ECL_RMOON_AU, 0.0026).unwrap_err(),
        Error::NonFiniteArgument { value } if value == 0.0
    ));
}

#[test]
fn occultation_phase_arms_finite_and_point_source() {
    use crate::constants::{ECL_ANNULAR, ECL_PARTIAL, ECL_TOTAL};
    // Finite body larger than the Moon: deep overlap is annular;
    // the total arm is unreachable there (total needs rbody < rmoon).
    assert_eq!(
        classify_occultation_phase(0.01, 0.30, 0.25).unwrap(),
        ECL_ANNULAR
    );
    assert_eq!(
        classify_occultation_phase(0.05, 0.20, 0.26).unwrap(),
        ECL_TOTAL
    );
    // Boundary (strict `<` arms): exactly on |rbody - rmoon| is partial.
    assert_eq!(
        classify_occultation_phase(0.06, 0.20, 0.26).unwrap(),
        ECL_PARTIAL
    );
    // Point-source star: annular arm unreachable, deep overlap total.
    assert_eq!(
        classify_occultation_phase(0.01, 0.0, 0.26).unwrap(),
        ECL_TOTAL
    );
    // Grazing overlap: partial.
    assert_eq!(
        classify_occultation_phase(0.40, 0.0, 0.26).unwrap(),
        ECL_PARTIAL
    );
    assert_eq!(
        classify_occultation_phase(0.40, 0.20, 0.26).unwrap(),
        ECL_PARTIAL
    );
    // In-progress gate mirrors the release skip rule.
    assert!(is_occultation_in_progress(0.26, 0.0, 0.26).unwrap());
    assert!(!is_occultation_in_progress(0.261, 0.0, 0.26).unwrap());
    assert!(matches!(
        classify_occultation_phase(f64::INFINITY, 0.0, 0.26).unwrap_err(),
        Error::NonFiniteArgument { .. }
    ));
}

#[test]
fn occultation_contact_residuals_and_star_coincidence() {
    // Outer residual: sep - (rb + rm); ~zero at external tangency
    // (0.25 + 0.26 need not be bit-identical to the 0.51 literal).
    assert!(
        occultation_outer_residual_deg(0.51, 0.25, 0.26)
            .unwrap()
            .abs()
            <= 1e-15
    );
    assert!((occultation_outer_residual_deg(0.40, 0.25, 0.26).unwrap() + 0.11).abs() < 1e-15);
    // Inner residual: |rb - rm| - sep; ~zero at internal tangency.
    assert!(
        occultation_inner_residual_deg(0.01, 0.25, 0.26)
            .unwrap()
            .abs()
            <= 1e-15
    );
    assert!((occultation_inner_residual_deg(0.0, 0.25, 0.26).unwrap() - 0.01).abs() < 1e-15);
    assert!(matches!(
        occultation_outer_residual_deg(f64::NAN, 0.25, 0.26).unwrap_err(),
        Error::NonFiniteArgument { .. }
    ));
    // Point-source coincidence, incl. the 0.0 no-contact case.
    assert_eq!(
        collapse_point_source_contacts(2460000.1, 2460000.2),
        (2460000.1, 2460000.2)
    );
    assert_eq!(collapse_point_source_contacts(0.0, 0.0), (0.0, 0.0));
}

#[test]
fn occult_where_composition_slots() {
    // Fixed where-core words with a total-geometry dcore.
    let moon = (0.9, 0.1, 0.0026);
    let sun = (0.0, 0.0, 1.0);
    let where_geom = eclipse_where_core_via(moon, sun, 100.0, 0.0).unwrap();
    let mut attr = [0.0; 20];
    attr[0] = 0.75;
    let out = compose_occult_where_outcome(&where_geom, attr);
    assert_eq!(out.retflag(), where_geom.retc());
    let geopos = out.geopos();
    assert_eq!(geopos[0], where_geom.lon_deg());
    assert_eq!(geopos[1], where_geom.lat_deg());
    assert!(geopos[2..].iter().all(|v| *v == 0.0));
    // attr[3] carries the core-shadow diameter, other slots preserved.
    assert_eq!(out.attr()[3], where_geom.dcore()[0]);
    assert_eq!(out.attr()[0], 0.75);
}

// conjunction seek, global/local occultation walks, contacts.
//
// All providers below are analytic scripts with known conjunctions,
// so every assertion pins the calculation order (Newton step,
// gates, windows, fallbacks, bit assembly) rather than an ephemeris.

/// Scripted global provider: fixed body lon/lat, Moon sweeping at a
/// constant rate, parabolic geocentric depth, scripted where-core and
/// right-ascension gap.
struct ScriptGlob {
    body_lon: f64,
    body_lat: f64,
    moon_rate: f64,
    moon_lat: f64,
    max_jd: f64,
    depth: f64,
    retc: i32,
    d2_base: f64,
    d2_curv: f64,
}

impl ScriptGlob {
    fn moon_lon(&self, jd: f64) -> f64 {
        self.body_lon + self.moon_rate * (jd - self.max_jd)
    }
    /// Nearest conjunction center: conjunctions sit where the sweep
    /// crosses `body_lon`, i.e. every `360/rate` days from `max_jd`.
    fn conj_center(&self, jd: f64) -> f64 {
        let period = 360.0 / self.moon_rate;
        let k = ((jd - self.max_jd) / period).round();
        self.max_jd + k * period
    }
}

impl OccultGlobalEphemeris for ScriptGlob {
    fn body_lonlat_deg(&self, jd: f64) -> Result<(f64, f64)> {
        let _ = jd;
        Ok((self.body_lon, self.body_lat))
    }
    fn moon_lonlat_deg(&self, jd: f64) -> Result<(f64, f64)> {
        Ok((self.moon_lon(jd), self.moon_lat))
    }
    fn sep_minus_radii_deg(&self, jd: f64) -> Result<f64> {
        let dt = jd - self.conj_center(jd);
        Ok(2.0 * dt * dt - self.depth)
    }
    fn where_core(&self, jd: f64) -> Result<WhereGeometry> {
        let dt = jd - self.conj_center(jd);
        let mut dcore = [0.0; 10];
        dcore[0] = 1.0;
        dcore[2] = self.d2_base + self.d2_curv * dt * dt;
        dcore[3] = 0.02;
        dcore[4] = 0.02;
        dcore[5] = 1e9;
        dcore[6] = 1e9;
        Ok(WhereGeometry {
            retc: self.retc,
            lon_deg: 10.0,
            lat_deg: 20.0,
            dcore,
        })
    }
    fn ra_gap_deg(&self, jd: f64) -> Result<f64> {
        Ok(wrap_ra_gap_deg(
            self.moon_rate * (jd - self.conj_center(jd)),
        ))
    }
}

#[test]
fn seek_converges_forward_and_backward() {
    // Body fixed at lon 100, Moon sweeping 13 deg/day through it.
    let mut body = |jd: f64| {
        let _ = jd;
        Ok((100.0, 1.0))
    };
    let mut moon = |jd: f64| Ok((90.0 + 13.0 * (jd - 1000.0), 0.5));
    let found = seek_moon_conjunction(&mut body, &mut moon, 1000.0, 1, false)
        .unwrap()
        .unwrap();
    // One Newton step of 10/13 d lands exactly on the conjunction.
    assert!((found.0 - (1000.0 + 10.0 / 13.0)).abs() < 1e-9);
    assert_eq!(found.1, 1.0);
    assert_eq!(found.2, 0.5);
    // Backward: the unfolded -350 deg gap steps back to the previous
    // conjunction; the folded residual is within tolerance.
    let mut body = |jd: f64| {
        let _ = jd;
        Ok((100.0, 1.0))
    };
    let mut moon = |jd: f64| Ok((90.0 + 13.0 * (jd - 1000.0), 0.5));
    let found = seek_moon_conjunction(&mut body, &mut moon, 1000.0, -1, false)
        .unwrap()
        .unwrap();
    assert!(found.0 < 1000.0);
    let gap = wrap_ra_gap_deg(100.0 - (90.0 + 13.0 * (found.0 - 1000.0)));
    assert!(gap.abs() <= OCC_CONJ_TOL_DEG);
}

#[test]
fn seek_star_gate_and_boundary_rules() {
    // Never-occultable star is refused, not walked past.
    let mut body = |jd: f64| {
        let _ = jd;
        Ok((100.0, 8.0))
    };
    let mut moon = |jd: f64| {
        let _ = jd;
        Ok((100.0, 0.0))
    };
    assert!(matches!(
        seek_moon_conjunction(&mut body, &mut moon, 1000.0, 1, true).unwrap_err(),
        Error::InvalidEclipseSearch { .. }
    ));
    // A first-evaluation coverage miss ends the walk (Ok(None))...
    let mut missing = |jd: f64| {
        let _ = jd;
        Err(Error::EphemerisRange {
            requested_jd: jd,
            start_jd: 0.0,
            end_jd: 900.0,
            body: Some(1),
        })
    };
    let mut moon = |jd: f64| {
        let _ = jd;
        Ok((100.0, 0.0))
    };
    assert!(
        seek_moon_conjunction(&mut missing, &mut moon, 1000.0, 1, false)
            .unwrap()
            .is_none()
    );
    // ... while an iteration failure propagates unchanged.
    let mut calls = 0;
    let mut flaky = |_jd: f64| {
        calls += 1;
        if calls == 1 {
            Ok((100.0, 1.0))
        } else {
            Err(Error::UnknownBody { body_id: 42 })
        }
    };
    let mut moon = |jd: f64| Ok((50.0 + 13.0 * (jd - 1000.0), 0.0));
    assert!(matches!(
        seek_moon_conjunction(&mut flaky, &mut moon, 1000.0, 1, false).unwrap_err(),
        Error::UnknownBody { body_id: 42 }
    ));
    // Non-finite epochs are refused.
    let mut body = |jd: f64| {
        let _ = jd;
        Ok((100.0, 1.0))
    };
    let mut moon = |jd: f64| {
        let _ = jd;
        Ok((100.0, 0.0))
    };
    assert!(matches!(
        seek_moon_conjunction(&mut body, &mut moon, f64::NAN, 1, false).unwrap_err(),
        Error::NonFiniteArgument { .. }
    ));
}

#[test]
fn wrap_ra_gap_pins() {
    assert_eq!(wrap_ra_gap_deg(10.0), 10.0);
    assert_eq!(wrap_ra_gap_deg(350.0), -10.0);
    assert_eq!(wrap_ra_gap_deg(180.0), 180.0);
    assert_eq!(wrap_ra_gap_deg(190.0), -170.0);
    assert_eq!(wrap_ra_gap_deg(-10.0), -10.0);
}

#[test]
fn glob_walk_hit_with_contacts_and_noon_slot() {
    use crate::constants::{ECL_CENTRAL, ECL_TOTAL};
    // Conjunction at JD 100, depth 0.01, TOTAL|CENTRAL geometry with
    // outer/total crossings inside the narrow windows.
    let eph = ScriptGlob {
        body_lon: 0.0,
        body_lat: 0.2,
        moon_rate: 13.0,
        moon_lat: 0.0,
        max_jd: 100.0,
        depth: 0.01,
        retc: ECL_TOTAL | ECL_CENTRAL,
        d2_base: 0.005,
        d2_curv: 2.0,
    };
    let (retflag, tret) =
        lun_occult_when_glob_via(&eph, 90.0, Some(5), 0, false, "Jupiter", false, false).unwrap();
    assert_eq!(retflag, ECL_TOTAL | ECL_CENTRAL);
    assert!((tret[0] - 100.0).abs() < 1e-6);
    // Outer crossings at d2 == ~0.01: |dt| == sqrt(0.005/2).
    assert!((tret[2] - (100.0 - 0.0500312)).abs() < 1e-4);
    assert!((tret[3] - (100.0 + 0.0500312)).abs() < 1e-4);
    // Totality crossings coincide here (same residual shape).
    assert!((tret[4] - tret[2]).abs() < 1e-9);
    assert!((tret[5] - tret[3]).abs() < 1e-9);
    // Central residual never crosses: the 0.0 sentinels stand.
    assert_eq!(tret[6], 0.0);
    assert_eq!(tret[7], 0.0);
    // The RA gap changes sign across begin/end: noon slot is set.
    assert!((tret[1] - 100.0).abs() < 1e-6);
    assert_eq!(tret[8], 0.0);
    assert_eq!(tret[9], 0.0);
}

#[test]
fn glob_walk_lat_gate_advances_and_one_try_suggests() {
    use crate::constants::{ECL_CENTRAL, ECL_TOTAL};
    // Moon latitude 5 deg before JD 110 (gate rejects the JD-100
    // conjunction), 0 afterwards (the JD-127.7 conjunction hits).
    struct LatGate {
        inner: ScriptGlob,
    }
    impl OccultGlobalEphemeris for LatGate {
        fn body_lonlat_deg(&self, jd: f64) -> Result<(f64, f64)> {
            self.inner.body_lonlat_deg(jd)
        }
        fn moon_lonlat_deg(&self, jd: f64) -> Result<(f64, f64)> {
            let (lon, _) = self.inner.moon_lonlat_deg(jd)?;
            Ok((lon, if jd < 110.0 { 5.0 } else { 0.0 }))
        }
        fn sep_minus_radii_deg(&self, jd: f64) -> Result<f64> {
            self.inner.sep_minus_radii_deg(jd)
        }
        fn where_core(&self, jd: f64) -> Result<WhereGeometry> {
            self.inner.where_core(jd)
        }
        fn ra_gap_deg(&self, jd: f64) -> Result<f64> {
            self.inner.ra_gap_deg(jd)
        }
    }
    let eph = LatGate {
        inner: ScriptGlob {
            body_lon: 0.0,
            body_lat: 0.2,
            moon_rate: 13.0,
            moon_lat: 0.0,
            max_jd: 100.0,
            depth: 0.01,
            retc: ECL_TOTAL | ECL_CENTRAL,
            d2_base: 0.005,
            d2_curv: 2.0,
        },
    };
    // Full walk skips JD 100 and reports the next conjunction.
    let (_, tret) =
        lun_occult_when_glob_via(&eph, 90.0, Some(5), 0, false, "Jupiter", false, false).unwrap();
    assert!((tret[0] - (100.0 + 360.0 / 13.0)).abs() < 1e-6);
    // One-try reports the miss with the continuation date.
    let (retflag, tret) =
        lun_occult_when_glob_via(&eph, 90.0, Some(5), 0, false, "Jupiter", false, true).unwrap();
    assert_eq!(retflag, 0);
    assert!((tret[0] - 101.0).abs() < 1e-6);
    assert!(tret[1..].iter().all(|v| *v == 0.0));
}

#[test]
fn glob_walk_epoch_margin_skips_and_exhaustion_reports() {
    use crate::constants::{ECL_CENTRAL, ECL_TOTAL};
    let eph = ScriptGlob {
        body_lon: 0.0,
        body_lat: 0.2,
        moon_rate: 13.0,
        moon_lat: 0.0,
        max_jd: 100.0,
        depth: 0.01,
        retc: ECL_TOTAL | ECL_CENTRAL,
        d2_base: 0.005,
        d2_curv: 2.0,
    };
    // Starting within the margin below the maximum skips to the next
    // conjunction rather than re-reporting it.
    let (_, tret) =
        lun_occult_when_glob_via(&eph, 99.999999, Some(5), 0, false, "Jupiter", false, false)
            .unwrap();
    assert!((tret[0] - (100.0 + 360.0 / 13.0)).abs() < 1e-6);
    // A permanently rejected geometry exhausts the walk honestly.
    // (body_lat 5 vs moon 0: every conjunction fails the 2 deg gate.)
    let never = ScriptGlob {
        body_lon: 0.0,
        body_lat: 5.0,
        moon_rate: 13.0,
        moon_lat: 0.0,
        max_jd: 100.0,
        depth: 0.01,
        retc: 0,
        d2_base: 0.005,
        d2_curv: 2.0,
    };
    assert!(matches!(
        lun_occult_when_glob_via(&never, 90.0, Some(5), 0, false, "Jupiter", false, false)
            .unwrap_err(),
        Error::EclipseSearchExhausted { .. }
    ));
    // The Moon itself and impossible filters never reach the walk.
    assert!(matches!(
        lun_occult_when_glob_via(&eph, 90.0, Some(1), 0, false, "Moon", false, false).unwrap_err(),
        Error::InvalidEclipseSearch { .. }
    ));
    assert!(matches!(
        lun_occult_when_glob_via(&eph, 90.0, Some(5), 8, false, "Jupiter", false, false)
            .unwrap_err(),
        Error::InvalidEclipseSearch { .. }
    ));
}

/// Scripted local provider: parabolic topocentric separation with
/// fixed radii, constant altitude, scripted rise/set and attr.
struct ScriptLoc {
    max_jd: f64,
    rbody: f64,
    rmoon: f64,
    rmoon_inner: f64,
    sep_base: f64,
    sep_curv: f64,
    alt: f64,
    rise: Option<f64>,
    set: Option<f64>,
    sun_rise: Option<f64>,
    sun_set: Option<f64>,
    retc: i32,
}

impl OccultLocalEphemeris for ScriptLoc {
    fn body_lonlat_deg(&self, jd: f64) -> Result<(f64, f64)> {
        let _ = jd;
        Ok((0.0, 0.2))
    }
    fn moon_lonlat_deg(&self, jd: f64) -> Result<(f64, f64)> {
        Ok((13.0 * (jd - self.max_jd), 0.0))
    }
    fn topo_radii_sep_deg(&self, jd: f64, inner: bool) -> Result<(f64, f64, f64)> {
        let dt = jd - self.max_jd;
        let rm = if inner { self.rmoon_inner } else { self.rmoon };
        Ok((self.rbody, rm, self.sep_base + self.sep_curv * dt * dt))
    }
    fn body_alt_deg(&self, jd: f64) -> Result<f64> {
        let _ = jd;
        Ok(self.alt)
    }
    fn local_how_attr(&self, jd: f64) -> Result<[f64; 20]> {
        let _ = jd;
        let mut attr = [0.5; 20];
        attr[0] = 1.5;
        attr[2] = 2.0;
        attr[6] = self.alt;
        Ok(attr)
    }
    fn where_core(&self, jd: f64) -> Result<WhereGeometry> {
        let _ = jd;
        let mut dcore = [0.0; 10];
        dcore[0] = -50.0;
        Ok(WhereGeometry {
            retc: self.retc,
            lon_deg: 10.0,
            lat_deg: 20.0,
            dcore,
        })
    }
    fn body_rise_set(&self, _jd_from: f64) -> Result<(Option<f64>, Option<f64>)> {
        Ok((self.rise, self.set))
    }
    fn sun_rise_set(&self, _jd: f64) -> Result<(Option<f64>, Option<f64>)> {
        Ok((self.sun_rise, self.sun_set))
    }
}

#[test]
fn loc_walk_planet_hit_with_bits_and_attr() {
    use crate::constants::{
        ECL_1ST_VISIBLE, ECL_2ND_VISIBLE, ECL_3RD_VISIBLE, ECL_4TH_VISIBLE, ECL_MAX_VISIBLE,
        ECL_NONCENTRAL, ECL_OCC_BEG_DAYLIGHT, ECL_OCC_END_DAYLIGHT, ECL_TOTAL, ECL_VISIBLE,
    };
    // Total geometry (min 0.005 < |0.25 - 0.26|), outer crossings at
    // +-0.1589 (wide fallback), inner at +-0.0474.
    let eph = ScriptLoc {
        max_jd: 100.0,
        rbody: 0.25,
        rmoon: 0.26,
        rmoon_inner: 0.20,
        sep_base: 0.005,
        sep_curv: 20.0,
        alt: 30.0,
        rise: Some(99.9),
        set: Some(100.1),
        sun_rise: Some(101.0),
        sun_set: Some(99.0),
        retc: ECL_TOTAL | ECL_NONCENTRAL,
    };
    let (retflag, tret, attr) = lun_occult_when_loc_via(&eph, 90.0, Some(5), false, false).unwrap();
    assert_eq!(
        retflag,
        ECL_TOTAL
            | ECL_VISIBLE
            | ECL_MAX_VISIBLE
            | ECL_1ST_VISIBLE
            | ECL_2ND_VISIBLE
            | ECL_3RD_VISIBLE
            | ECL_4TH_VISIBLE
            | ECL_OCC_BEG_DAYLIGHT
            | ECL_OCC_END_DAYLIGHT
            | ECL_NONCENTRAL
    );
    assert!((tret[0] - 100.0).abs() < 1e-6);
    assert!((tret[1] - (100.0 - 0.15892)).abs() < 1e-4);
    assert!((tret[2] - (100.0 - 0.04743)).abs() < 1e-4);
    assert!((tret[3] - (100.0 + 0.04743)).abs() < 1e-4);
    assert!((tret[4] - (100.0 + 0.15892)).abs() < 1e-4);
    assert_eq!(tret[5], 99.9);
    assert_eq!(tret[6], 100.1);
    assert_eq!(tret[7], 0.0);
    // Attr assembly: unity caps plus the dcore[0] overwrite.
    assert_eq!(attr[0], 1.0);
    assert_eq!(attr[2], 1.0);
    assert_eq!(attr[3], -50.0);
    assert_eq!(attr[6], 30.0);
}

#[test]
fn loc_walk_star_coincidence_and_invisible_miss() {
    use crate::constants::{ECL_NONCENTRAL, ECL_TOTAL};
    // Star (point source): outer contacts coincide with the inner.
    // (Steeper separation so the inner radius is crossed inside the
    // 1.5 h window: |0 - 0.20| - sep roots at +-0.057 d.)
    let star = ScriptLoc {
        max_jd: 100.0,
        rbody: 0.0,
        rmoon: 0.26,
        rmoon_inner: 0.20,
        sep_base: 0.005,
        sep_curv: 60.0,
        alt: 30.0,
        rise: Some(99.9),
        set: Some(100.1),
        sun_rise: Some(101.0),
        sun_set: Some(102.0),
        retc: ECL_TOTAL | ECL_NONCENTRAL,
    };
    let (retflag, tret, _) = lun_occult_when_loc_via(&star, 90.0, None, false, false).unwrap();
    assert!((retflag & ECL_TOTAL) != 0);
    assert_eq!(tret[1], tret[2]);
    assert_eq!(tret[4], tret[3]);
    assert_ne!(tret[1], 0.0);
    assert!((tret[1] - (100.0 - 0.05701)).abs() < 1e-4);
    // A body that never clears the horizon: one-try miss carries the
    // continuation date and a zero attr.
    // (alt drives both the visibility gate and attr[6] here.)
    let dark = ScriptLoc {
        max_jd: 100.0,
        rbody: 0.25,
        rmoon: 0.26,
        rmoon_inner: 0.20,
        sep_base: 0.005,
        sep_curv: 20.0,
        alt: -10.0,
        rise: Some(99.9),
        set: Some(100.1),
        sun_rise: Some(101.0),
        sun_set: Some(99.0),
        retc: ECL_TOTAL | ECL_NONCENTRAL,
    };
    let (retflag, tret, attr) = lun_occult_when_loc_via(&dark, 90.0, Some(5), false, true).unwrap();
    assert_eq!(retflag, 0);
    assert!((tret[0] - 101.0).abs() < 1e-4);
    assert!(attr.iter().all(|v| *v == 0.0));
}

#[test]
fn loc_attr_finish_caps_and_preserves_nan() {
    use crate::constants::ECL_NONCENTRAL;
    let mut dcore = [0.0; 10];
    dcore[0] = -80.0;
    let geom = WhereGeometry {
        retc: ECL_NONCENTRAL,
        lon_deg: 0.0,
        lat_deg: 0.0,
        dcore,
    };
    let mut attr = [0.25; 20];
    attr[0] = f64::NAN;
    attr[2] = 1.5;
    let (retflag, out) = finish_occult_local_attr(4, attr, &geom);
    assert_eq!(retflag, 4 | ECL_NONCENTRAL);
    //  NaN-preserving upper bound: NaN is preserved, 1.5 caps.
    assert!(out[0].is_nan());
    assert_eq!(out[2], 1.0);
    assert_eq!(out[3], -80.0);
}

// Planetary occultations over analytic providers
struct ScriptPlanet {
    rate_deg_day: f64,
    max_jd: f64,
    dec_off_deg: f64,
    occ_r_deg: f64,
    tgt_r_deg: f64,
    alt_deg: f64,
    az_deg: f64,
    fail_from: Option<f64>,
    fail_other: bool,
    /// Optional coverage edge gating ONLY the angular-radius
    /// acquisitions (coverage probe: the golden minimum samples
    /// separations but the refinement check also reads radii).
    fail_radii_from: Option<f64>,
    /// Optional coverage edge gating ONLY the altitude
    /// acquisitions (coverage probe: the global walk never reads
    /// altitudes, the local visibility scan does).
    fail_alt_from: Option<f64>,
}

impl ScriptPlanet {
    fn check(&self, jd: f64) -> Result<()> {
        if self.fail_other {
            return Err(Error::UnknownBody { body_id: 3 });
        }
        if let Some(edge) = self.fail_from {
            if jd >= edge {
                return Err(Error::EphemerisRange {
                    requested_jd: jd,
                    start_jd: 0.0,
                    end_jd: edge,
                    body: None,
                });
            }
        }
        Ok(())
    }
    fn base() -> Self {
        ScriptPlanet {
            rate_deg_day: 1.0,
            max_jd: 100.0,
            dec_off_deg: 0.005,
            occ_r_deg: 0.02,
            tgt_r_deg: 0.01,
            alt_deg: 30.0,
            az_deg: 180.0,
            fail_from: None,
            fail_other: false,
            fail_radii_from: None,
            fail_alt_from: None,
        }
    }
    fn range_miss(&self, jd: f64, edge: f64) -> Error {
        Error::EphemerisRange {
            requested_jd: jd,
            start_jd: 0.0,
            end_jd: edge,
            body: None,
        }
    }
}

impl PlanetOccultEphemeris for ScriptPlanet {
    fn occulting_ra_dec_deg(&self, jd: f64) -> Result<(f64, f64)> {
        self.check(jd)?;
        Ok((self.rate_deg_day * (jd - self.max_jd), 0.0))
    }
    fn target_ra_dec_deg(&self, jd: f64) -> Result<(f64, f64)> {
        self.check(jd)?;
        Ok((0.0, self.dec_off_deg))
    }
    fn occulting_radius_deg(&self, jd: f64) -> Result<f64> {
        self.check(jd)?;
        if let Some(edge) = self.fail_radii_from {
            if jd >= edge {
                return Err(self.range_miss(jd, edge));
            }
        }
        Ok(self.occ_r_deg)
    }
    fn target_radius_deg(&self, jd: f64) -> Result<f64> {
        self.check(jd)?;
        if let Some(edge) = self.fail_radii_from {
            if jd >= edge {
                return Err(self.range_miss(jd, edge));
            }
        }
        Ok(self.tgt_r_deg)
    }
}

impl PlanetOccultLocalEphemeris for ScriptPlanet {
    fn occulting_alt_az_deg(&self, jd: f64) -> Result<(f64, f64)> {
        self.check(jd)?;
        if let Some(edge) = self.fail_alt_from {
            if jd >= edge {
                return Err(self.range_miss(jd, edge));
            }
        }
        Ok((self.alt_deg, self.az_deg))
    }
    fn target_alt_az_deg(&self, jd: f64) -> Result<(f64, f64)> {
        self.check(jd)?;
        if let Some(edge) = self.fail_alt_from {
            if jd >= edge {
                return Err(self.range_miss(jd, edge));
            }
        }
        Ok((self.alt_deg, self.az_deg + 1.0))
    }
    fn local_radii_sep_deg(&self, jd: f64) -> Result<(f64, f64, f64)> {
        Ok((
            self.occulting_radius_deg(jd)?,
            self.target_radius_deg(jd)?,
            planet_occult_separation_via(self, jd)?,
        ))
    }
}

#[test]
fn planet_separation_identity_quadrant_wrap_and_poles() {
    assert_eq!(
        planet_angular_separation_deg(10.0, 20.0, 10.0, 20.0).unwrap(),
        0.0
    );
    assert!((planet_angular_separation_deg(0.0, 0.0, 90.0, 0.0).unwrap() - 90.0).abs() < 1e-12);
    assert!((planet_angular_separation_deg(0.0, 0.0, 180.0, 0.0).unwrap() - 180.0).abs() < 1e-12);
    assert!((planet_angular_separation_deg(359.0, 0.0, 1.0, 0.0).unwrap() - 2.0).abs() < 1e-9);
    assert!((planet_angular_separation_deg(0.0, 89.0, 180.0, 89.0).unwrap() - 2.0).abs() < 1e-9);
    // Coincident directions retain even sub-microarcsecond offsets.
    let tiny = planet_angular_separation_deg(0.0, 0.0, 1e-13, 0.0).unwrap();
    assert!(tiny.is_finite() && tiny < 1e-11);
    // Non-finite words are typed refusals (documented divergence).
    assert!(matches!(
        planet_angular_separation_deg(f64::NAN, 0.0, 0.0, 0.0),
        Err(Error::NonFiniteArgument { .. })
    ));
}

#[test]
fn planet_separation_resolves_small_equatorial_and_antipodal_offsets() {
    // Along the equator the central angle is exactly the longitude offset.
    // Powers of two are exact input degrees, independently of a sky model.
    // The same identity applies along a meridian. Test both ends of the
    // angular domain, where inversion of a rounded cosine loses precision.
    for exponent in [-21, -43, -35, -27, -14, -7, 0, 6] {
        let offset = 2.0_f64.powi(exponent);
        for (ra, dec, exact) in [
            (offset, 0.0, offset),
            (0.0, offset, offset),
            (180.0 - offset, 0.0, 180.0 - offset),
        ] {
            let actual = planet_angular_separation_deg(0.0, 0.0, ra, dec).unwrap();
            let roundoff = 64.0 * f64::EPSILON * exact.max(1.0);
            assert!(
                (actual - exact).abs() <= roundoff,
                "RA {ra}, Dec {dec}: expected {exact}, obtained {actual}"
            );
        }
    }
}

#[test]
fn planet_occultation_does_not_merge_resolved_disjoint_discs() {
    let script = ScriptPlanet {
        dec_off_deg: 0.0,
        occ_r_deg: 2.0_f64.powi(-27),
        tgt_r_deg: 2.0_f64.powi(-27),
        ..ScriptPlanet::base()
    };
    // Constant unit equatorial motion makes the center offset exactly
    // 2^-21 degrees, 32 times the combined radius. Rounding its cosine to
    // one used to fabricate overlap between these physically disjoint discs.
    let (overlap, separation, occulting_radius, target_radius) =
        check_planet_occultation_via(&script, 100.0 + 2.0_f64.powi(-21)).unwrap();
    assert!(!overlap);
    assert!(separation > occulting_radius + target_radius);
}

#[test]
fn planet_occult_threshold_and_classification_edges() {
    // The overlap rule is `<=`: exact tangency counts as a hit.
    assert!(planet_occult_threshold_met(0.03, 0.02, 0.01));
    assert!(!planet_occult_threshold_met(0.0300001, 0.02, 0.01));
    // The total rule is strict `<` on the absolute radius
    // difference: tangency of the inner contact is partial.
    use crate::constants::{ECL_PARTIAL, ECL_TOTAL, JUPITER, MOON, SUN, VENUS};
    assert_eq!(classify_planet_occultation(0.005, 0.02, 0.01), ECL_TOTAL);
    assert_eq!(classify_planet_occultation(0.01, 0.02, 0.01), ECL_PARTIAL);
    assert_eq!(classify_planet_occultation(0.029, 0.02, 0.01), ECL_PARTIAL);
    // Request validation arms with the release messages.
    assert!(validate_planet_occult_request(VENUS, Some(JUPITER), "").is_ok());
    assert!(validate_planet_occult_request(VENUS, None, "Regulus").is_ok());
    let missing = validate_planet_occult_request(VENUS, None, "").unwrap_err();
    assert!(format!("{missing:?}").contains("Either occulted_planet ID or starname"));
    let sun = validate_planet_occult_request(SUN, Some(JUPITER), "").unwrap_err();
    assert!(format!("{sun:?}").contains("Sun cannot be an occulting body"));
    let moon = validate_planet_occult_request(MOON, Some(JUPITER), "").unwrap_err();
    assert!(format!("{moon:?}").contains("Moon cannot be an occulting body"));
    let bad = validate_planet_occult_request(12345, Some(JUPITER), "").unwrap_err();
    assert!(format!("{bad:?}").contains("Invalid occulting planet ID: 12345"));
    let bad_tgt = validate_planet_occult_request(VENUS, Some(12345), "").unwrap_err();
    assert!(format!("{bad_tgt:?}").contains("Invalid occulted planet ID: 12345"));
    let same = validate_planet_occult_request(VENUS, Some(VENUS), "").unwrap_err();
    assert!(format!("{same:?}").contains("cannot be the same"));
}

#[test]
fn planet_occult_minimum_and_contact_chords() {
    use crate::constants::{JUPITER, VENUS};
    let script = ScriptPlanet::base();
    // Golden minimum recovers the conjunction epoch; the 0.005-deg
    // Dec offset is the minimum separation.
    let jd_max = planet_occult_minimum_via(&script, 99.5, 100.5).unwrap();
    assert!((jd_max - 100.0).abs() < 1e-6);
    let (hit, min_sep, occ_r, tgt_r) = check_planet_occultation_via(&script, jd_max).unwrap();
    assert!(hit);
    assert!((min_sep - 0.005).abs() < 1e-9);
    assert_eq!((occ_r, tgt_r), (0.02, 0.01));
    // Chord contacts: outer half-width sqrt(0.03² − 0.005²) at
    // 1 deg/day, inner half-width sqrt(0.01² − 0.005²); symmetric
    // about the maximum (dec = 0 keeps cos(dec) = 1, as in the
    // release unchanged-behavior arm).
    let (first, second, third, fourth) =
        planet_occult_contact_times_via(&script, jd_max, min_sep, occ_r, tgt_r).unwrap();
    let outer = (0.03f64 * 0.03 - 0.005 * 0.005).sqrt();
    let inner = (0.01f64 * 0.01 - 0.005 * 0.005).sqrt();
    assert!((jd_max - first - outer).abs() < 1e-6);
    assert!((fourth - jd_max - outer).abs() < 1e-6);
    assert!((jd_max - second - inner).abs() < 1e-6);
    assert!((third - jd_max - inner).abs() < 1e-6);
    // A minimum past the outer threshold yields zero sentinels.
    let (f0, s0, t0, fo0) =
        planet_occult_contact_times_via(&script, jd_max, 0.05, occ_r, tgt_r).unwrap();
    assert_eq!((f0, s0, t0, fo0), (0.0, 0.0, 0.0, 0.0));
    // Full global walk from JD 90: total hit with the release
    // 10-slot layout (slot 1 reserved zero).
    let (ecl_type, tret) =
        planet_occult_when_glob_via(&script, 90.0, VENUS, Some(JUPITER), "", 0).unwrap();
    use crate::constants::ECL_TOTAL;
    assert_eq!(ecl_type, ECL_TOTAL);
    assert!((tret[0] - 100.0).abs() < 1e-6);
    assert_eq!(tret[1], 0.0);
    assert!(tret[2] < tret[0] && tret[0] < tret[3]);
    assert!(tret[4] < tret[0] && tret[0] < tret[5]);
    assert!(tret[6..].iter().all(|v| *v == 0.0));
}

#[test]
fn event_positive_mod_preserves_divisor_sign_and_ieee_edges() {
    use super::geometry::positive_mod;
    for value in [0.0, -0.0, 360.0, -360.0, 720.0, -720.0] {
        assert_eq!(positive_mod(value, 360.0).to_bits(), 0.0_f64.to_bits());
    }
    for (value, expected) in [
        (-359.75, 0.25),
        (359.75, 359.75),
        (-0.25, 359.75),
        (720.25, 0.25),
    ] {
        assert_eq!(positive_mod(value, 360.0), expected);
    }
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(positive_mod(value, 360.0).is_nan());
    }
    // Floating addition may round a negative tiny remainder to the divisor;
    // the operation does not clamp that valid floating result to zero.
    assert_eq!(positive_mod(-f64::MIN_POSITIVE, 360.0), 360.0);
}

#[test]
fn planet_occult_contact_wraps_each_body_in_both_directions() {
    struct Pair {
        occ: [f64; 2],
        target: [f64; 2],
    }
    impl PlanetOccultEphemeris for Pair {
        fn occulting_ra_dec_deg(&self, jd: f64) -> Result<(f64, f64)> {
            Ok((self.occ[usize::from(jd > 2_451_545.0)], 0.0))
        }
        fn target_ra_dec_deg(&self, jd: f64) -> Result<(f64, f64)> {
            Ok((self.target[usize::from(jd > 2_451_545.0)], 0.0))
        }
        fn occulting_radius_deg(&self, _: f64) -> Result<f64> {
            Ok(0.1)
        }
        fn target_radius_deg(&self, _: f64) -> Result<f64> {
            Ok(0.01)
        }
    }
    // Literal contact words from the admitted Source wrap witnesses. The
    // reverse controls have the same speed magnitude and contact ordering.
    let expected = [
        0x4142_b42c_7a22_2222,
        0x4142_b42c_7b33_3333,
        0x4142_b42c_84cc_cccd,
        0x4142_b42c_85dd_ddde,
    ];
    for (occ, target) in [
        ([359.9, 0.1], [0.0, 0.0]),
        ([0.0, 0.0], [0.1, 359.9]),
        ([0.1, 359.9], [0.0, 0.0]),
        ([0.0, 0.0], [359.9, 0.1]),
    ] {
        let (a, b, c, d) =
            planet_occult_contact_times_via(&Pair { occ, target }, 2_451_545.0, 0.0, 0.1, 0.01)
                .unwrap();
        assert_eq!(
            [a.to_bits(), b.to_bits(), c.to_bits(), d.to_bits()],
            expected
        );
    }
    for occ in [[10.0, 10.25], [10.25, 10.0]] {
        let actual = planet_occult_contact_times_via(
            &Pair {
                occ,
                target: [0.0, 0.0],
            },
            2_451_545.0,
            0.0,
            0.1,
            0.01,
        )
        .unwrap();
        let outer = (0.1 + 0.01) / 3.0;
        let inner = (0.1_f64 - 0.01).abs() / 3.0;
        assert_eq!(
            actual,
            (
                2_451_545.0 - outer,
                2_451_545.0 - inner,
                2_451_545.0 + inner,
                2_451_545.0 + outer
            )
        );
    }
}

#[test]
fn planet_occult_glob_backward_and_partial() {
    use crate::constants::{ECL_PARTIAL, JUPITER, VENUS};
    // Backward search from JD 110 finds the same conjunction.
    let script = ScriptPlanet::base();
    let (_, tret) =
        planet_occult_when_glob_via(&script, 110.0, VENUS, Some(JUPITER), "", -1).unwrap();
    assert!((tret[0] - 100.0).abs() < 1e-6);
    // A 0.02-deg offset sits between the inner (0.01) and outer
    // (0.03) thresholds: partial with zero inner contacts.
    let partial = ScriptPlanet {
        dec_off_deg: 0.02,
        ..ScriptPlanet::base()
    };
    let (ecl_type, tret) =
        planet_occult_when_glob_via(&partial, 90.0, VENUS, Some(JUPITER), "", 0).unwrap();
    assert_eq!(ecl_type, ECL_PARTIAL);
    assert_ne!(tret[2], 0.0);
    assert_eq!((tret[4], tret[5]), (0.0, 0.0));
}

#[test]
fn planet_occult_glob_exhaustion_and_range_policy() {
    use crate::constants::{JUPITER, VENUS};
    // No conjunction ever: the walk runs the iteration bound and
    // raises the release exhaustion message.
    let quiet = ScriptPlanet {
        rate_deg_day: 0.0,
        dec_off_deg: 5.0,
        ..ScriptPlanet::base()
    };
    let err = planet_occult_when_glob_via(&quiet, 90.0, VENUS, Some(JUPITER), "", 0).unwrap_err();
    assert!(matches!(err, Error::EclipseSearchExhausted { .. }));
    assert!(format!("{err:?}").contains("No planetary occultation of planet 5 by planet 3"));
    // A provider range miss ends the walk fail-closed (exhaustion,
    // never a fabricated event); any other provider error
    // propagates untouched.
    let gapped = ScriptPlanet {
        fail_from: Some(95.0),
        ..ScriptPlanet::base()
    };
    let err = planet_occult_when_glob_via(&gapped, 90.0, VENUS, Some(JUPITER), "", 0).unwrap_err();
    assert!(matches!(err, Error::EclipseSearchExhausted { .. }));
    let broken = ScriptPlanet {
        fail_other: true,
        ..ScriptPlanet::base()
    };
    let err = planet_occult_when_glob_via(&broken, 90.0, VENUS, Some(JUPITER), "", 0).unwrap_err();
    assert!(matches!(err, Error::UnknownBody { .. }));
}

#[test]
fn planet_glob_refinement_range_miss_propagates_fix_c28() {
    // FIX-C28: the refinement check sits outside the release
    // `try` (only the main sample is guarded). A coverage miss on
    // the refinement check must propagate fail-closed — the walk
    // must NOT end it in exhaustion over an uncomputed minimum.
    // The radii-only edge keeps the golden-minimum separation
    // samples (RA/Dec) served while the refinement check at the
    // minimum reads a missing radius.
    use crate::constants::{JUPITER, VENUS};
    let gapped = ScriptPlanet {
        fail_radii_from: Some(99.99),
        ..ScriptPlanet::base()
    };
    let err = planet_occult_when_glob_via(&gapped, 90.0, VENUS, Some(JUPITER), "", 0).unwrap_err();
    assert!(matches!(err, Error::EphemerisRange { .. }));
}

#[test]
fn planet_loc_visibility_miss_breaks_to_exhaustion_fix_c28() {
    // FIX-C28: the release wraps the whole per-candidate local
    // body in `except (Error, RuntimeError)` ending the search. A
    // coverage miss on the visibility scan must break into
    // exhaustion — never propagate out of the local walk. The
    // altitude-only edge keeps the global re-search (RA/Dec and
    // radii) served so the candidate reaches the visibility scan.
    use crate::constants::{JUPITER, VENUS};
    let dark_gap = ScriptPlanet {
        fail_alt_from: Some(0.0),
        ..ScriptPlanet::base()
    };
    let err = planet_occult_when_loc_via(&dark_gap, 90.0, VENUS, Some(JUPITER), "", 41.9, 12.5)
        .unwrap_err();
    assert!(matches!(err, Error::EclipseSearchExhausted { .. }));
    assert!(format!("{err:?}").contains("visible from lat=41.9, lon=12.5"));
}

#[test]
fn planet_occult_loc_visible_hit_and_attr() {
    use crate::constants::{JUPITER, VENUS};
    let script = ScriptPlanet::base();
    let (retflag, times, attr) =
        planet_occult_when_loc_via(&script, 90.0, VENUS, Some(JUPITER), "", 41.9, 12.5).unwrap();
    use crate::constants::ECL_TOTAL;
    assert_eq!(retflag, ECL_TOTAL);
    // Reordered local layout: max, first, second, third, fourth.
    assert!((times[0] - 100.0).abs() < 1e-6);
    assert!(times[1] < times[0] && times[0] < times[4]);
    assert!(times[2] < times[0] && times[0] < times[3]);
    assert!(times[5..].iter().all(|v| *v == 0.0));
    // Attr composition: magnitude in (0, 1], ratio 2.0,
    // obscuration == magnitude, az/alt/apparent-alt/separation.
    assert!(attr[0] > 0.0 && attr[0] <= 1.0);
    assert_eq!(attr[1], 2.0);
    assert_eq!(attr[2], attr[0]);
    assert_eq!(attr[3], 0.0);
    assert_eq!(attr[4], 181.0);
    assert_eq!(attr[5], 30.0);
    assert!((attr[6] - 30.029099750246328).abs() < 1e-12);
    assert!((attr[7] - 0.005).abs() < 1e-9);
    assert!(attr[8..].iter().all(|v| *v == 0.0));
}

#[test]
fn planet_occult_loc_star_ratio_and_night_skip() {
    use crate::constants::VENUS;
    // Star target: the 0.0001-deg point source yields the 999.0
    // diameter ratio.
    let mut star = ScriptPlanet::base();
    star.tgt_r_deg = OCC_POINT_SOURCE_RADIUS_DEG;
    let (_, _, attr) =
        planet_occult_when_loc_via(&star, 90.0, VENUS, None, "Regulus", 41.9, 12.5).unwrap();
    assert_eq!(attr[1], PLANET_OCCULT_POINT_SOURCE_RATIO);
    // Both bodies below the −0.5-deg floor: the candidate is skipped
    // and, with a single conjunction on the script, the walk
    // exhausts with the release site message.
    let dark = ScriptPlanet {
        alt_deg: -10.0,
        ..ScriptPlanet::base()
    };
    let err =
        planet_occult_when_loc_via(&dark, 90.0, VENUS, None, "Regulus", 41.9, 12.5).unwrap_err();
    assert!(matches!(err, Error::EclipseSearchExhausted { .. }));
    assert!(format!("{err:?}").contains("visible from lat=41.9, lon=12.5"));
}

#[test]
fn planet_occult_apparent_altitude_gate() {
    // Past −1 deg the true altitude travels untouched; above it the independently computed
    // Saemundsson value holds to 1e-12.
    assert_eq!(planet_occult_apparent_altitude_deg(-5.0), -5.0);
    assert!((planet_occult_apparent_altitude_deg(45.0) - 45.016878460981225).abs() < 1e-12);
}

// -- tests ----------------------------------------------------
#[test]
fn sampling_step_validation_arms() {
    // Nominal minute step converts to days exactly.
    assert_eq!(
        eclipse_sampling_step_days(2460000.0, 2460001.0, 1.0, "probe").unwrap(),
        1.0 / 1440.0
    );
    // Non-finite window bounds are range errors naming the bad bound.
    let err = eclipse_sampling_step_days(f64::NAN, 2460001.0, 1.0, "probe").unwrap_err();
    assert!(matches!(
        err,
        Error::EphemerisRange { requested_jd, .. } if requested_jd.is_nan()
    ));
    assert!(matches!(
        eclipse_sampling_step_days(2460000.0, f64::INFINITY, 1.0, "probe"),
        Err(Error::EphemerisRange { .. })
    ));
    // NaN, zero, negative-zero and negative steps are refused before
    // any walk starts (NaN would otherwise slip the sign test).
    for bad in [f64::NAN, 0.0, -0.0, -2.5] {
        let err = eclipse_sampling_step_days(2460000.0, 2460001.0, bad, "probe").unwrap_err();
        assert!(matches!(err, Error::InvalidEclipseSearch { .. }), "{bad}");
    }
    // Resolution gate: at JD 2460000 consecutive doubles are
    // ~4.7e-10 days apart, so 1e-7 minutes cannot advance the
    // accumulator while 1 minute can.
    let err = eclipse_sampling_step_days(2460000.0, 2460001.0, 1e-7, "probe").unwrap_err();
    assert!(format!("{err:?}").contains("cannot advance"));
    // Exactly half the spacing is refused (round-half-even tie);
    // just above it is accepted.
    let pivot = 2460001.0f64;
    let spacing = pivot.next_up() - pivot;
    let tie_minutes = (spacing / 2.0) * MINUTES_PER_DAY;
    assert!(eclipse_sampling_step_days(pivot, pivot, tie_minutes, "probe").is_err());
    assert!(eclipse_sampling_step_days(pivot, pivot, tie_minutes.next_up(), "probe").is_ok());
}

#[test]
fn parse_geopos_keeps_release_messages() {
    // Short sequences are refused with the caller-specific message.
    let err =
        parse_geopos(&[1.0, 2.0], "sol_eclipse_when_loc", GEOPOS_SEQUENCE_MESSAGE).unwrap_err();
    assert!(matches!(err, Error::InvalidEclipseSearch { .. }));
    assert!(format!("{err:?}").contains("[longitude, latitude, altitude]"));
    let err = parse_geopos(&[1.0], "sol_eclipse_how_details", GEOPOS_ELEMENTS_MESSAGE).unwrap_err();
    assert!(format!("{err:?}").contains("[lon, lat, alt]"));
    // Order is (lon, lat, alt); extras are ignored like the release.
    assert_eq!(
        parse_geopos(
            &[-96.797, 32.7767, 0.0, 99.0],
            "probe",
            GEOPOS_SEQUENCE_MESSAGE
        )
        .unwrap(),
        (-96.797, 32.7767, 0.0)
    );
}

#[test]
fn path_longitude_normalization_pins() {
    assert_eq!(normalize_path_longitude_deg(190.0), -170.0);
    assert_eq!(normalize_path_longitude_deg(-190.0), 170.0);
    assert_eq!(normalize_path_longitude_deg(180.0), -180.0);
    assert_eq!(normalize_path_longitude_deg(0.0), -180.0 + 180.0);
    assert_eq!(normalize_path_longitude_deg(360.0), 0.0);
}

struct ScriptPath {
    retflag: i32,
    lon: f64,
    lat: f64,
    x: f64,
    y: f64,
}

impl EclipsePathEphemeris for ScriptPath {
    fn path_where(&self, _: f64) -> Result<(i32, f64, f64)> {
        Ok((self.retflag, self.lon, self.lat))
    }
    fn path_besselian_xy(&self, _: f64) -> Result<(f64, f64)> {
        Ok((self.x, self.y))
    }
}

#[test]
fn central_line_collects_central_only() {
    let central = ScriptPath {
        retflag: ECL_CENTRAL | ECL_TOTAL,
        lon: 10.0,
        lat: 20.0,
        x: 0.1,
        y: 0.2,
    };
    let (times, lats, lons) =
        eclipse_central_line_via(&central, 2460000.0, 2460000.125, 60.0).unwrap();
    // 3-hour window at 60-minute steps visits 4 instants.
    assert_eq!(times.len(), 4);
    assert!(lats.iter().all(|&v| v == 20.0));
    assert!(lons.iter().all(|&v| v == 10.0));
    // Partial-only eclipses (no ECL_CENTRAL bit) contribute nothing.
    let partial = ScriptPath {
        retflag: ECL_PARTIAL | ECL_NONCENTRAL,
        lon: 10.0,
        lat: 20.0,
        x: 0.1,
        y: 0.2,
    };
    let (times, _, _) = eclipse_central_line_via(&partial, 2460000.0, 2460000.125, 60.0).unwrap();
    assert!(times.is_empty());
    // A grazing gamma (>= 1.5) or a null island point is omitted too.
    let graze = ScriptPath {
        retflag: ECL_CENTRAL,
        lon: 10.0,
        lat: 20.0,
        x: 1.2,
        y: 1.2,
    };
    let (times, _, _) = eclipse_central_line_via(&graze, 2460000.0, 2460000.125, 60.0).unwrap();
    assert!(times.is_empty());
    // The step gate fires before any provider call.
    assert!(eclipse_central_line_via(&central, 2460000.0, 2460001.0, 0.0).is_err());
}

#[test]
fn limit_point_symmetry_and_gates() {
    // y = 0 with a small umbra: north/south mirror exactly.
    let north = eclipse_limit_point_deg(0.0, 0.0, 0.0, 0.0, 0.01, false).unwrap();
    let south = eclipse_limit_point_deg(0.0, 0.0, 0.0, 0.0, 0.01, true).unwrap();
    assert!((north.0 + south.0).abs() < 1e-12);
    assert_eq!(north.1, south.1);
    assert!(north.0 > 0.0 && north.0 < 90.0);
    // Zero umbra (point shadow) and a missed Earth omit the sample.
    assert!(eclipse_limit_point_deg(0.0, 0.0, 0.0, 0.0, 0.0, false).is_none());
    assert!(eclipse_limit_point_deg(2.0, 0.0, 0.0, 0.0, 0.01, false).is_none());
    assert!(eclipse_limit_point_deg(f64::NAN, 0.0, 0.0, 0.0, 0.01, false).is_none());
}

struct ScriptLimit {
    elements: (f64, f64, f64, f64, f64),
}

impl EclipseLimitEphemeris for ScriptLimit {
    fn limit_besselian(&self, _: f64) -> Result<(f64, f64, f64, f64, f64)> {
        Ok(self.elements)
    }
}

#[test]
fn limit_walks_share_one_core() {
    let hit = ScriptLimit {
        elements: (0.3, 0.2, 10.0, 50.0, -0.01),
    };
    let (nt, nlats, nlons) = eclipse_northern_limit_via(&hit, 2460000.0, 2460000.1, 30.0).unwrap();
    let (st, slats, slons) = eclipse_southern_limit_via(&hit, 2460000.0, 2460000.1, 30.0).unwrap();
    assert_eq!(nt.len(), st.len());
    assert!(!nt.is_empty());
    // Same instants, mirrored edges: northern latitudes sit above the
    // southern ones while longitudes stay near (the edge offset
    // shifts them slightly, ~0.15 deg here — not identical).
    for ((nlat, nlon), (slat, slon)) in nlats.iter().zip(&nlons).zip(slats.iter().zip(&slons)) {
        assert!(nlat > slat);
        assert!((nlon - slon).abs() < 1.0);
    }
    for v in nt.iter().chain(nlats.iter()).chain(nlons.iter()) {
        assert!(v.is_finite());
    }
    // Step errors name the calling sampler.
    let err = eclipse_northern_limit_via(&hit, 2460000.0, 2460001.0, -1.0).unwrap_err();
    assert!(format!("{err:?}").contains("calc_eclipse_northern_limit"));
    let err = eclipse_southern_limit_via(&hit, 2460000.0, 2460001.0, -1.0).unwrap_err();
    assert!(format!("{err:?}").contains("calc_eclipse_southern_limit"));
}

struct ScriptWidth {
    retflag: i32,
    sample: Option<PathWidthSample>,
    calls: std::cell::Cell<usize>,
    /// Coverage probe (a): fail `path_width_sample` with a coverage
    /// miss (`Some(true)`, LEB range) or a non-range provider
    /// error (`Some(false)`).
    sample_fail: Option<bool>,
}

impl EclipsePathWidthEphemeris for ScriptWidth {
    fn path_width_where(&self, _: f64) -> Result<(i32, f64, f64)> {
        Ok((self.retflag, 0.0, 0.0))
    }
    fn path_width_sample(&self, _: f64, _: f64, _: f64) -> Result<PathWidthSample> {
        self.calls.set(self.calls.get() + 1);
        if let Some(leb) = self.sample_fail {
            if leb {
                return Err(Error::LebOutOfRange {
                    body_id: 10,
                    jd: 2460000.0,
                    start: 0.0,
                    end: 2450000.0,
                });
            }
            return Err(Error::UnknownBody { body_id: 10 });
        }
        // (the script always carries a sample; missing-sample arms
        // return before this call)
        Ok(PathWidthSample {
            sun_dist_au: self.sample.as_ref().map_or(1.0, |s| s.sun_dist_au),
            moon_dist_au: self.sample.as_ref().map_or(0.0024, |s| s.moon_dist_au),
            sun_alt_true_deg: self.sample.as_ref().map_or(60.0, |s| s.sun_alt_true_deg),
            separation_deg: self.sample.as_ref().map_or(0.0, |s| s.separation_deg),
        })
    }
}

fn width_script(retflag: i32, sample: PathWidthSample) -> ScriptWidth {
    ScriptWidth {
        retflag,
        sample: Some(sample),
        calls: std::cell::Cell::new(0),
        sample_fail: None,
    }
}

#[test]
fn path_width_gate_and_geometry_arms() {
    // Non-central where-arm: 0.0 with no site sample taken.
    let idle = width_script(
        ECL_PARTIAL | ECL_NONCENTRAL,
        PathWidthSample {
            sun_dist_au: 1.0,
            moon_dist_au: 0.0024,
            sun_alt_true_deg: 60.0,
            separation_deg: 0.0,
        },
    );
    assert_eq!(
        eclipse_path_width_via(&idle, 2460000.0, None, None).unwrap(),
        0.0
    );
    assert_eq!(idle.calls.get(), 0);
    // Sun below the horizon: 0.0.
    let dark = width_script(
        ECL_CENTRAL,
        PathWidthSample {
            sun_dist_au: 1.0,
            moon_dist_au: 0.0024,
            sun_alt_true_deg: -5.0,
            separation_deg: 0.0,
        },
    );
    assert_eq!(
        eclipse_path_width_via(&dark, 2460000.0, None, None).unwrap(),
        0.0
    );
    // Discs apart: 0.0.
    let apart = width_script(
        ECL_CENTRAL,
        PathWidthSample {
            sun_dist_au: 1.0,
            moon_dist_au: 0.0024,
            sun_alt_true_deg: 60.0,
            separation_deg: 5.0,
        },
    );
    assert_eq!(
        eclipse_path_width_via(&apart, 2460000.0, None, None).unwrap(),
        0.0
    );
    // Total geometry at 60 deg altitude: a plausible path width.
    let total = width_script(
        ECL_CENTRAL,
        PathWidthSample {
            sun_dist_au: 1.0,
            moon_dist_au: 0.0024,
            sun_alt_true_deg: 60.0,
            separation_deg: 0.0,
        },
    );
    let width = eclipse_path_width_via(&total, 2460000.0, None, None).unwrap();
    assert!(width > 50.0 && width < 400.0, "{width}");
    // Grazing cone with a huge shadow: capped at 1500 km.
    let graze = width_script(
        ECL_CENTRAL,
        PathWidthSample {
            sun_dist_au: 1.0,
            moon_dist_au: 0.002,
            sun_alt_true_deg: 0.1,
            separation_deg: 0.0,
        },
    );
    assert_eq!(
        eclipse_path_width_via(&graze, 2460000.0, None, None).unwrap(),
        PATH_WIDTH_MAX_KM
    );
    // Explicit site skips the central gate (where arm must not run).
    let site = width_script(
        ECL_PARTIAL,
        PathWidthSample {
            sun_dist_au: 1.0,
            moon_dist_au: 0.0024,
            sun_alt_true_deg: 60.0,
            separation_deg: 0.0,
        },
    );
    assert!(eclipse_path_width_via(&site, 2460000.0, Some(32.0), Some(-96.0)).unwrap() > 0.0);
    assert_eq!(site.calls.get(), 1);
    // Non-finite epochs are refused.
    assert!(eclipse_path_width_via(&total, f64::NAN, None, None).is_err());
}

#[test]
fn path_width_sample_errors_follow_leaf_policy_fix_c31() {
    // FIX-C31(a): the release topo arm maps non-range acquisition
    // failures to 0.0 and re-raises LEB range errors only. A
    // coverage miss on the sample propagates fail-closed; any
    // other provider error degrades to 0.0, never a fabricated
    // width and never a propagated non-range error.
    let total = || {
        width_script(
            ECL_CENTRAL,
            PathWidthSample {
                sun_dist_au: 1.0,
                moon_dist_au: 0.0024,
                sun_alt_true_deg: 60.0,
                separation_deg: 0.0,
            },
        )
    };
    let mut gapped = total();
    gapped.sample_fail = Some(true);
    assert!(matches!(
        eclipse_path_width_via(&gapped, 2460000.0, None, None),
        Err(Error::LebOutOfRange { .. })
    ));
    let mut broken = total();
    broken.sample_fail = Some(false);
    assert_eq!(
        eclipse_path_width_via(&broken, 2460000.0, None, None).unwrap(),
        0.0
    );
}

struct ScriptDuration {
    retflag: i32,
    contacts: (f64, f64),
}

impl SolarDurationEphemeris for ScriptDuration {
    fn duration_where(&self, _: f64) -> Result<(i32, f64, f64)> {
        Ok((self.retflag, -96.0, 32.0))
    }
    fn duration_local_contacts(&self, _: f64, _: f64, _: f64) -> Result<(f64, f64)> {
        Ok(self.contacts)
    }
}

#[test]
fn solar_duration_minutes_and_zeros() {
    let central = ScriptDuration {
        retflag: ECL_CENTRAL | ECL_TOTAL,
        contacts: (2460000.0, 2460000.002),
    };
    // 0.002 days = 2.88 minutes of central phase (JD words near
    // 2.46e6 carry ~3.4e-7-minute representation noise).
    assert!((solar_eclipse_duration_min_via(&central, 2460000.0).unwrap() - 2.88).abs() < 1e-6);
    let partial = ScriptDuration {
        retflag: ECL_PARTIAL | ECL_NONCENTRAL,
        contacts: (2460000.0, 2460000.002),
    };
    assert_eq!(
        solar_eclipse_duration_min_via(&partial, 2460000.0).unwrap(),
        0.0
    );
    let missing = ScriptDuration {
        retflag: ECL_CENTRAL,
        contacts: (0.0, 2460000.002),
    };
    assert_eq!(
        solar_eclipse_duration_min_via(&missing, 2460000.0).unwrap(),
        0.0
    );
}

#[test]
fn max_time_dispatcher_arms() {
    let jd0 = 2460000.0;
    // Global arm: golden gamma minimum over a parabola.
    let mut gamma = |jd: f64| Ok((jd - jd0).powi(2) + 0.25);
    let mut sep = |_: f64| Ok(0.0);
    let (jd_max, gamma_at) =
        eclipse_max_time_via(&mut gamma, &mut sep, jd0 + 0.05, None, 0.125).unwrap();
    assert!((jd_max - jd0).abs() < 1e-6);
    assert!((gamma_at - 0.25).abs() < 1e-9);
    // Local arm delegates to the separation provider.
    let mut gamma_unused = |_: f64| -> Result<f64> { panic!("global arm not taken") };
    let mut vsep = |jd: f64| Ok((jd - jd0).abs() * 12.0);
    let (local_max, min_sep) = eclipse_max_time_via(
        &mut gamma_unused,
        &mut vsep,
        jd0 + 0.05,
        Some((32.0, -96.0)),
        0.125,
    )
    .unwrap();
    assert!((local_max - jd0).abs() < 1e-6);
    assert!(min_sep < 1e-6);
    // Degenerate windows and epochs are refused before any search.
    assert!(eclipse_max_time_via(&mut gamma, &mut sep, jd0, None, 0.0).is_err());
    assert!(eclipse_max_time_via(&mut gamma, &mut sep, f64::NAN, None, 0.125).is_err());
}

#[test]
fn kernel_bounds_min_max_and_empty() {
    assert_eq!(
        kernel_jd_bounds(&[(3.0, 4.0), (1.0, 9.0), (2.0, 5.0)]),
        Some((1.0, 9.0))
    );
    assert_eq!(kernel_jd_bounds(&[]), None);
    assert_eq!(kernel_jd_bounds(&[(5.0, 6.0)]), Some((5.0, 6.0)));
}

#[test]
fn position_angle_cardinal_pins() {
    // Moon due east of the Sun: 90 deg; due north: 0 deg.
    assert!((position_angle_deg(0.0, 0.0, 90.0, 0.0).unwrap() - 90.0).abs() < 1e-12);
    assert!((position_angle_deg(0.0, 0.0, 0.0, 45.0).unwrap()).abs() < 1e-12);
    // Moon 10 deg behind in RA on the equator: due west, 270 deg
    // (the negative atan2 folds through the +360 arm).
    assert!((position_angle_deg(10.0, 0.0, 0.0, 0.0).unwrap() - 270.0).abs() < 1e-9);
    assert!(position_angle_deg(f64::NAN, 0.0, 0.0, 0.0).is_err());
}

#[test]
fn details_pure_helper_arms() {
    // Disjoint discs: no obscuration.
    assert_eq!(details_max_obscuration(0.26, 0.27, 0.6), 0.0);
    // Larger-Moon containment: full cover.
    assert_eq!(details_max_obscuration(0.26, 0.27, 0.001), 1.0);
    // Larger-Sun containment: the disc-area ratio.
    let contained = details_max_obscuration(0.27, 0.26, 0.001);
    assert!((contained - (0.26f64 / 0.27).powi(2)).abs() < 1e-15);
    // Partial overlap lands strictly inside (0, 1).
    let partial = details_max_obscuration(0.26, 0.27, 0.3);
    assert!(partial > 0.0 && partial < 1.0);
    // Magnitude clamps at zero for a miss and follows overlap/diameter.
    assert_eq!(details_max_magnitude(0.52, 0.54, 1.0), 0.0);
    let got = details_max_magnitude(0.52, 0.54, 0.0);
    let want = (0.53f64 / 0.52).min(1.0 + 0.54 / 0.52);
    assert!((got - want).abs() < 1e-15);
    // Contact bisection finds a linear crossing point.
    let mut f = |jd: f64| Ok((jd - 2460000.0).abs() * 12.0);
    let c = details_contact_time(&mut f, 2459999.9, 2460000.0, 0.53, false).unwrap();
    assert!((c - (2460000.0 - 0.53 / 12.0)).abs() < 1e-6);
}

struct ScriptDetails {
    initial: (i32, [f64; 20]),
    jmax: f64,
    rate: f64,
    radii: (f64, f64),
    /// Coverage probe (b): fail `details_sun_alt_az_deg` from the
    /// second call on (the first call serves the maximum
    /// composition; contact tails start at call two) with a
    /// coverage miss (`fail_altaz_leb`) or a non-range provider
    /// error (`fail_altaz_other`).
    altaz_calls: std::cell::Cell<usize>,
    fail_altaz_leb: bool,
    fail_altaz_other: bool,
}

impl EclipseHowDetailsEphemeris for ScriptDetails {
    fn details_initial(&self, _: f64) -> Result<(i32, [f64; 20])> {
        Ok(self.initial)
    }
    fn details_separation_deg(&self, jd_ut: f64) -> Result<f64> {
        Ok((jd_ut - self.jmax).abs() * self.rate)
    }
    fn details_sun_alt_az_deg(&self, _: f64) -> Result<(f64, f64)> {
        let n = self.altaz_calls.get() + 1;
        self.altaz_calls.set(n);
        if n >= 2 {
            if self.fail_altaz_leb {
                return Err(Error::LebOutOfRange {
                    body_id: 10,
                    jd: self.jmax,
                    start: 0.0,
                    end: self.jmax - 1.0,
                });
            }
            if self.fail_altaz_other {
                return Err(Error::UnknownBody { body_id: 10 });
            }
        }
        Ok((60.0, 180.0))
    }
    fn details_disc_radii_deg(&self, _: f64) -> Result<(f64, f64)> {
        Ok(self.radii)
    }
    fn details_equatorial_deg(&self, _: f64) -> Result<(f64, f64, f64, f64)> {
        Ok((0.0, 0.0, 0.1, 0.1))
    }
}

#[test]
fn how_details_total_composition() {
    let jmax = 2460409.5;
    let mut attr = [0.0f64; 20];
    attr[0] = 1.0;
    attr[1] = 1.05;
    attr[2] = 1.0;
    attr[3] = 100.0;
    attr[4] = 180.0;
    attr[5] = 60.0;
    attr[7] = 0.0;
    let total = ScriptDetails {
        initial: (ECL_TOTAL | ECL_VISIBLE, attr),
        jmax,
        rate: 12.0,
        radii: (0.26, 0.27),
        altaz_calls: std::cell::Cell::new(0),
        fail_altaz_leb: false,
        fail_altaz_other: false,
    };
    let out = eclipse_how_details_via(&total, jmax, -96.797, 32.7767, 0.0).unwrap();
    // Instant slots come from the initial call; the maximum refines.
    assert_eq!(out.magnitude, 1.0);
    assert_eq!(out.obscuration_percent, 100.0);
    assert!((out.jd_max - jmax).abs() < 1e-6);
    // Central total: ordered contacts, full cover, re-classified flags.
    assert!(out.jd_c1 < out.jd_c2 && out.jd_c2 < out.jd_c3 && out.jd_c3 < out.jd_c4);
    assert!(out.jd_c1 > 0.0 && out.jd_c2 > 0.0 && out.jd_c3 > 0.0 && out.jd_c4 > 0.0);
    assert!(out.is_total && !out.is_partial && !out.is_annular);
    assert_eq!(out.max_obscuration, 1.0);
    assert!((out.ratio - 0.54 / 0.52).abs() < 1e-12);
    assert!(out.duration_partial_minutes > out.duration_total_minutes);
    assert!(out.duration_total_minutes > 0.0);
    for pa in [
        out.position_angle_c1,
        out.position_angle_c2,
        out.position_angle_c3,
        out.position_angle_c4,
    ] {
        assert!((0.0..360.0).contains(&pa));
    }
    assert_eq!(out.sun_alt_max, 60.0);
}

#[test]
fn how_details_contact_leb_miss_propagates_fix_c31() {
    // FIX-C31(b): the release re-raises LEB range errors out of
    // the per-contact blocks while other provider failures skip
    // the block tail. A coverage miss on a contact tail
    // propagates fail-closed; a non-range tail failure keeps the
    // solved contact with zeroed tail words.
    let jmax = 2460409.5;
    let mut attr = [0.0f64; 20];
    attr[0] = 1.0;
    attr[1] = 1.05;
    attr[2] = 1.0;
    attr[3] = 100.0;
    attr[4] = 180.0;
    attr[5] = 60.0;
    attr[7] = 0.0;
    let script = |leb: bool, other: bool| ScriptDetails {
        initial: (ECL_TOTAL | ECL_VISIBLE, attr),
        jmax,
        rate: 12.0,
        radii: (0.26, 0.27),
        altaz_calls: std::cell::Cell::new(0),
        fail_altaz_leb: leb,
        fail_altaz_other: other,
    };
    assert!(matches!(
        eclipse_how_details_via(&script(true, false), jmax, -96.797, 32.7767, 0.0),
        Err(Error::LebOutOfRange { .. })
    ));
    // UnknownBody is a public Error, outside the protected builtin hierarchy.
    // The original contact handler does not swallow this provider failure.
    assert!(matches!(
        eclipse_how_details_via(&script(false, true), jmax, -96.797, 32.7767, 0.0),
        Err(Error::UnknownBody { body_id: 10 })
    ));
}

#[test]
fn how_details_no_eclipse_returns_early() {
    // type 0: the initial assembly returns before any site provider
    // runs, so even failing providers are never touched.
    struct Silent;
    impl EclipseHowDetailsEphemeris for Silent {
        fn details_initial(&self, _: f64) -> Result<(i32, [f64; 20])> {
            Ok((0, [0.0; 20]))
        }
        fn details_separation_deg(&self, _: f64) -> Result<f64> {
            panic!("must not run without an eclipse")
        }
        fn details_sun_alt_az_deg(&self, _: f64) -> Result<(f64, f64)> {
            panic!("must not run without an eclipse")
        }
        fn details_disc_radii_deg(&self, _: f64) -> Result<(f64, f64)> {
            panic!("must not run without an eclipse")
        }
        fn details_equatorial_deg(&self, _: f64) -> Result<(f64, f64, f64, f64)> {
            panic!("must not run without an eclipse")
        }
    }
    let out = eclipse_how_details_via(&Silent, 2460000.0, 0.0, 0.0, 0.0).unwrap();
    assert_eq!(out.eclipse_type, 0);
    assert!(!out.is_visible);
    assert_eq!(out.jd_max, 0.0);
    assert_eq!(out.duration_partial_minutes, 0.0);
}

/// Named-star occultation with the independent bundled star clock.
/// Fixed-time circumstances use the pinned values; the selected maximum is
/// additionally tested against the derivative of the actual limb-overlap
/// geometry, because reducing propagation rounding can improve its epoch.
#[test]
fn lunloc_star_regulus_uses_bundled_star_clock() {
    let Some(bytes) = leb2_core_fixture() else {
        eprintln!("lunloc_star_regulus_uses_bundled_star_clock: SKIP (missing LEB2 core fixture)");
        return;
    };
    let index = crate::leb::FileIndex::decode(&bytes).expect("core LEB decodes");
    let cfg = crate::time::DeltaTConfig {
        model: crate::time::DeltaTModel::EspenakMeeus,
        user_override_days: None,
        tidal_acceleration: crate::constants::TIDAL_DEFAULT,
        iers_enabled: false,
        iers_lookup_seconds: None,
        smh2016_lookup_seconds: None,
    };
    let provider = crate::positions::SingleLeb::new(&index, &bytes);
    let nutation_at = |jd: f64| index.eval_nutation(&bytes, jd);
    let serving = crate::positions::SealedServing::new(&provider, &nutation_at);
    let deltat_at = |jd_ut: f64| crate::time::deltat_with_config(jd_ut, &cfg);
    let star_serving = crate::stars::SealedStarServing::new(
        &provider,
        &nutation_at,
        &deltat_at,
        crate::sidereal::SiderealSelection::unset(),
    )
    .with_default_timescale();
    let bound = SealedOccultLocalServing::try_new_star(
        &serving,
        &star_serving,
        &deltat_at,
        &nutation_at,
        crate::constants::FLG_SWIEPH,
        "Regulus",
        f64::from_bits(0xc058c886594af4f1),
        f64::from_bits(0x40336ebedfa43fe6),
        f64::from_bits(0x40a1800000000000),
    )
    .expect("Regulus local-occult serving binds");
    let (retflag, tret, attr) =
        lun_occult_when_loc_sealed(&bound, f64::from_bits(0x4142c412d6594af5), false, false)
            .expect("Regulus local-occult walk settles");
    assert_eq!(retflag, 16260, "occultation type travels back exact");
    let pin_tret: [u64; 10] = [
        0x4142c6f2443da0a4,
        0x4142c6f24212527c,
        0x4142c6f24212527c,
        0x4142c6f24671194c,
        0x4142c6f24671194c,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
    ];
    for (slot, (got, want)) in tret.iter().zip(pin_tret.iter()).enumerate() {
        let want = f64::from_bits(*want);
        if want == 0.0 {
            assert_eq!(*got, 0.0, "reserved event-time slot {slot}");
        } else {
            assert!(
                (got - want).abs() * 86400.0 <= 1.0,
                "event time {slot}: {got} vs pinned {want}"
            );
        }
    }
    let pin_attr: [u64; 20] = [
        0x3ff0000000000000,
        0x0000000000000000,
        0x3ff0000000000000,
        0xc0ab289999d44d53,
        0x4071b52fb3f419aa,
        0x40500ed493726773,
        0x40500f38284c9f2a,
        0x3fcc269ca5781491,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
    ];
    // Compare angles at the same original maximum, not at two different
    // event epochs. Source and native azimuth/altitudes agree exactly when
    // evaluated at either common time; their raw maximum epochs differ.
    let source_max = f64::from_bits(pin_tret[0]);
    let fixed_how = OccultLocalEphemeris::local_how_attr(&bound, source_max).unwrap();
    let fixed_where = OccultLocalEphemeris::where_core(&bound, source_max).unwrap();
    let (_, fixed_attr) = finish_occult_local_attr(retflag, fixed_how, &fixed_where);
    for (slot, (got, word)) in fixed_attr.into_iter().zip(pin_attr).enumerate() {
        let expected = f64::from_bits(word);
        let limit = match slot {
            3 => expected.abs() * 1e-8, // core width in km
            4..=7 => 0.001 / 3600.0,    // azimuth, altitudes, separation in degrees
            0..=2 | 8 => 1e-6,          // dimensionless geometry
            _ => 0.0,                   // reserved/discrete slots
        };
        assert!(
            (got - expected).abs() <= limit,
            "fixed-time attribute {slot}: {got} vs pinned {expected}; physical limit {limit}"
        );
        eprintln!(
            "Regulus raw event attr[{slot}] {} vs pinned {expected}; same-time {got}",
            attr[slot]
        );
    }
    // The maximum is the minimum of the outer-limb residual, not center
    // separation alone. Independent degree-six derivatives of actual
    // positions and radii check stationarity at three fixed dyadic spacings.
    let overlap = |jd| {
        let (body_radius, moon_radius, separation) =
            OccultLocalEphemeris::topo_radii_sep_deg(&bound, jd, false).unwrap();
        separation - body_radius - moon_radius
    };
    for h in [1.0 / 1024.0, 1.0 / 2048.0, 1.0 / 4096.0] {
        let correction_seconds = |t| {
            let center = overlap(t);
            let delta = |n: f64| overlap(t + n * h) - center;
            let first = (45.0 * (delta(1.0) - delta(-1.0)) - 9.0 * (delta(2.0) - delta(-2.0))
                + delta(3.0)
                - delta(-3.0))
                / (60.0 * h);
            let second = (270.0 * (delta(1.0) + delta(-1.0)) - 27.0 * (delta(2.0) + delta(-2.0))
                + 2.0 * (delta(3.0) + delta(-3.0)))
                / (180.0 * h * h);
            assert!(second > 0.0, "overlap must have a local minimum");
            -first / second * 86400.0
        };
        let native_error = correction_seconds(tret[0]).abs();
        let source_error = correction_seconds(source_max).abs();
        assert!(
            native_error <= 1.0,
            "maximum is outside the one-second event criterion"
        );
        assert!(
            native_error < source_error,
            "native maximum correction {native_error}s is not smaller than pinned {source_error}s"
        );
    }
}

#[test]
fn eclipse_details_retains_successful_geometry_when_position_angle_is_unresolved() {
    struct Observed;
    impl EclipseHowDetailsEphemeris for Observed {
        fn details_initial(&self, _: f64) -> Result<(i32, [f64; 20])> {
            Ok((ECL_TOTAL | ECL_VISIBLE, [0.0; 20]))
        }
        fn details_separation_deg(&self, jd: f64) -> Result<f64> {
            Ok((jd - 100.0).abs() * 12.0)
        }
        fn details_sun_alt_az_deg(&self, _: f64) -> Result<(f64, f64)> {
            Ok((60.0, 180.0))
        }
        fn details_disc_radii_deg(&self, _: f64) -> Result<(f64, f64)> {
            Ok((0.26, 0.27))
        }
        fn details_equatorial_deg(&self, _: f64) -> Result<(f64, f64, f64, f64)> {
            Err(Error::MathDomain { value: -1.0 })
        }
    }
    let out = eclipse_how_details_via(&Observed, 100.0, 0.0, 0.0, 0.0).unwrap();
    assert_eq!(out.sun_angular_radius, 0.26);
    assert_eq!(out.moon_angular_radius, 0.27);
    assert!((out.ratio - 0.27 / 0.26).abs() < 1e-15);
    assert_eq!(out.sun_alt_c1, 60.0);
    assert_eq!(out.sun_az_c1, 180.0);
    assert_eq!(out.sun_alt_c4, 60.0);
    assert_eq!(out.sun_az_c4, 180.0);
    // A caught position-angle error follows successful contact/horizon writes.
    assert!(out.jd_c1 > 0.0);
    assert_eq!(out.position_angle_c1, 0.0);
    // The larger lunar disc fully covers the Sun at maximum.
    assert_eq!(out.max_obscuration, 1.0);
}
