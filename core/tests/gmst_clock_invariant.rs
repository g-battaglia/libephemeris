//! Independent UT1 rotation contract for the sidereal origin (issue #50).
//!
//! A fixed celestial frame at TT cannot change Earth's UT1 spin rate. Advance
//! UT1 by one whole day at the same representable TT and boundary offsets:
//! only the nonintegral ERA rotation remains modulo one turn. This is not a
//! trajectory along a physical Delta-T history, and needs no ephemeris data.

use libephemeris_core::time::{
    SIDEREAL_MODERN_END_JD, SIDEREAL_MODERN_START_JD, SiderealDeltaT, mean_sidereal_time_deg,
    sidereal_continuity_offsets,
};

// PLAN's existing time tolerance, expressed as an angular displacement.
const ONE_MICROSECOND_SIDEREAL_ARCSEC: f64 = 0.000_015;
const ERA_FRACTIONAL_DEG_PER_UT1_DAY: f64 = 360.0 * 0.002_737_811_911_354_48;

fn rotation_errors_arcsec(epochs: &[f64]) -> Vec<(f64, f64)> {
    let delta = 64.0 / 86_400.0;
    let (lo, hi) = sidereal_continuity_offsets(delta, delta).unwrap();
    // Capitaine, Guinot & McCarthy (2000): ERA advances 1.00273781191135448
    // turns per UT1 day. Subtract the integral turn before binary64 evaluation.
    let expected_deg = ERA_FRACTIONAL_DEG_PER_UT1_DAY;
    epochs
        .iter()
        .map(|&ut| {
            let tt = ut + delta;
            let next_ut = ut + 1.0;
            let first = SiderealDeltaT::try_new(tt - ut, lo, hi).unwrap();
            let second = SiderealDeltaT::try_new(tt - next_ut, lo, hi).unwrap();
            assert_eq!(
                (ut + first.delta_t_days()).to_bits(),
                (next_ut + second.delta_t_days()).to_bits()
            );
            let a = mean_sidereal_time_deg(ut, &first).unwrap();
            let b = mean_sidereal_time_deg(next_ut, &second).unwrap();
            let error = ((b - a - expected_deg + 180.0).rem_euclid(360.0) - 180.0) * 3600.0;
            (ut, error)
        })
        .collect()
}

fn assert_rotation_contract(epochs: &[f64]) {
    let failures: Vec<_> = rotation_errors_arcsec(epochs)
        .into_iter()
        // NaN must fail too: an ordered comparison alone would silently pass it.
        .filter(|(_, error)| !error.is_finite() || error.abs() > ONE_MICROSECOND_SIDEREAL_ARCSEC)
        .collect();
    assert!(
        failures.is_empty(),
        "fixed-TT UT1 rotation failures: {failures:?}"
    );
}

#[test]
fn modern_mean_origin_preserves_independent_ut1_rotation() {
    assert_rotation_contract(&[2_396_759.5, 2_451_545.0, 2_469_805.5]);
}

#[test]
fn longterm_mean_origin_must_preserve_independent_ut1_rotation() {
    assert_rotation_contract(&[
        -3_100_015.5,
        SIDEREAL_MODERN_START_JD - 1.0,
        // These boundaries are positive finite numbers, so adjacent bit
        // patterns are their representable neighbours (also on Rust 1.86).
        // A one-day pair can cross the dispatcher: the physical contract must
        // hold across that transition too, not just inside a single branch.
        f64::from_bits(SIDEREAL_MODERN_START_JD.to_bits() - 1),
        SIDEREAL_MODERN_START_JD,
        f64::from_bits(SIDEREAL_MODERN_START_JD.to_bits() + 1),
        f64::from_bits(SIDEREAL_MODERN_END_JD.to_bits() - 1),
        SIDEREAL_MODERN_END_JD,
        f64::from_bits(SIDEREAL_MODERN_END_JD.to_bits() + 1),
        2_488_147.5,
        2_524_671.5,
        2_561_195.5,
        2_597_720.5,
        8_000_016.5,
    ]);
}
