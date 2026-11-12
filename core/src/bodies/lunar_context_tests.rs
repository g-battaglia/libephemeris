//! Owned synthetic lunar wiring regressions moved from the consumer.
use super::*;
use crate::constants::{
    FLG_BARYCTR, FLG_EQUATORIAL, FLG_HELCTR, FLG_SIDEREAL, FLG_SPEED, FLG_SWIEPH, MARS, MEAN_APOG,
    MEAN_NODE, OSCU_APOG, TRUE_NODE,
};
use crate::lunar::{no_nutation_frame, wants_lunar_helio_bary_zero};
/// Deterministic synthetic provider for configuration/order tests.
fn stub_dpsi(_: f64) -> LebResult<f64> {
    Ok(0.01)
}

/// Deterministic synthetic provider for configuration/order tests.
fn stub_nut_rate(_: f64, _: f64) -> LebResult<f64> {
    Ok(0.0)
}

/// Deterministic synthetic provider for configuration/order tests.
fn stub_moon_state(_: f64) -> LebResult<EclipticMoonState> {
    Ok(((0.0025, 0.0002, 0.0001), (-0.0002, 0.0165, 0.0007)))
}

#[test]
fn true_epoch_refusal_preserves_provider_order() {
    // the model shaping, data-free: without FLG_SPEED all three
    // speeds are exactly zero while the shaped positions are
    // untouched; a NaN epoch refuses before the Moon-state
    // provider is contacted.
    let speedy = true_lunar_state_with_providers(
        OSCU_APOG,
        2460000.5,
        FLG_SWIEPH | FLG_SPEED,
        &stub_moon_state,
        &stub_dpsi,
        &stub_nut_rate,
    )
    .expect("true lilith must serve");
    assert!(speedy[3].is_finite() && speedy[4].is_finite() && speedy[5].is_finite());
    let still = true_lunar_state_with_providers(
        OSCU_APOG,
        2460000.5,
        FLG_SWIEPH,
        &stub_moon_state,
        &stub_dpsi,
        &stub_nut_rate,
    )
    .expect("true lilith must serve");
    assert_eq!((still[3], still[4], still[5]), (0.0, 0.0, 0.0));
    assert_eq!(
        (still[0], still[1], still[2]),
        (speedy[0], speedy[1], speedy[2])
    );
    let touched = std::cell::Cell::new(false);
    let spy = |_: f64| {
        touched.set(true);
        stub_moon_state(0.0)
    };
    let error = true_lunar_state_with_providers(
        TRUE_NODE,
        f64::NAN,
        FLG_SWIEPH | FLG_SPEED,
        &spy,
        &stub_dpsi,
        &stub_nut_rate,
    )
    .expect_err("NaN epoch must refuse");
    assert!(
        format!("{error:?}").contains("NonFiniteArgument"),
        "wrong refusal: {error:?}"
    );
    // the calculation order: the central Moon-state acquisition precedes
    // the epoch gate, unlike the mean path whose range check runs
    // first.
    assert!(touched.get(), "provider runs before the epoch gate");
}

#[test]
fn true_central_matches_public_geometry_unnutated() {
    // Wiring identity, data-free: the true geometry consumes
    // true-ecliptic-of-date states, so nutation is already in the
    // raw values — unlike the mean path, the shaped central passes
    // the public raw geometry through unchanged (no dpsi shift).
    // The constant-provider stencil leaves exactly zero speeds
    // even under FLG_SPEED.
    use crate::lunar::{true_lilith_from_state, true_lunar_node_from_state};
    let (r, v) = stub_moon_state(2460000.5).expect("stub state");
    let (raw_node_lon, raw_node_lat, raw_node_dist) =
        true_lunar_node_from_state(r, v).expect("raw node geometry");
    assert_eq!(raw_node_lat, 0.0);
    let shaped = true_lunar_state_with_providers(
        TRUE_NODE,
        2460000.5,
        FLG_SWIEPH | FLG_SPEED,
        &stub_moon_state,
        &stub_dpsi,
        &stub_nut_rate,
    )
    .expect("true node must serve");
    assert_eq!(shaped[0], raw_node_lon);
    assert_eq!(shaped[1], raw_node_lat);
    assert_eq!(shaped[2], raw_node_dist);
    assert_eq!((shaped[3], shaped[4], shaped[5]), (0.0, 0.0, 0.0));
    let (raw_lil_lon, raw_lil_lat, raw_lil_dist) =
        true_lilith_from_state(r, v).expect("raw lilith geometry");
    let shaped = true_lunar_state_with_providers(
        OSCU_APOG,
        2460000.5,
        FLG_SWIEPH | FLG_SPEED,
        &stub_moon_state,
        &stub_dpsi,
        &stub_nut_rate,
    )
    .expect("true lilith must serve");
    assert_eq!(shaped[0], raw_lil_lon);
    assert_eq!(shaped[1], raw_lil_lat);
    assert_eq!(shaped[2], raw_lil_dist);
    assert_eq!((shaped[3], shaped[4], shaped[5]), (0.0, 0.0, 0.0));
}

#[test]
fn true_state_composition_rotates_the_geocentric_difference() {
    // Pure frame composition, data-free: the rotation preserves
    // vector lengths, actually rotates (the of-date frame is not
    // the identity), and swapping Moon/Earth negates the state.
    let jd_tt = 2460000.5_f64;
    let (r, v) = ecliptic_moon_state_of(
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        jd_tt,
        0.0,
        0.0,
    )
    .expect("frame composition must serve");
    let r_len = (r.0 * r.0 + r.1 * r.1 + r.2 * r.2).sqrt();
    let v_len = (v.0 * v.0 + v.1 * v.1 + v.2 * v.2).sqrt();
    assert!(
        (r_len - 1.0).abs() < 1e-12,
        "rotation preserves length: {r:?}"
    );
    assert!(
        (v_len - 1.0).abs() < 1e-12,
        "rotation preserves length: {v:?}"
    );
    let drift = ((r.0 - 1.0).powi(2) + r.1.powi(2) + r.2.powi(2)).sqrt();
    assert!(
        drift > 1e-6,
        "of-date frame must differ from identity: {r:?}"
    );
    let (r_swap, v_swap) = ecliptic_moon_state_of(
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        jd_tt,
        0.0,
        0.0,
    )
    .expect("swapped composition must serve");
    for (a, b) in [r_swap.0, r_swap.1, r_swap.2]
        .into_iter()
        .zip([r.0, r.1, r.2])
    {
        assert!((a + b).abs() < 1e-12, "swap negates: {r_swap:?} vs {r:?}");
    }
    for (a, b) in [v_swap.0, v_swap.1, v_swap.2]
        .into_iter()
        .zip([v.0, v.1, v.2])
    {
        assert!((a + b).abs() < 1e-12, "swap negates: {v_swap:?} vs {v:?}");
    }
    let error = ecliptic_moon_state_of(
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        f64::NAN,
        0.0,
        0.0,
    )
    .expect_err("NaN epoch must refuse");
    assert!(
        format!("{error:?}").contains("NonFiniteArgument"),
        "wrong refusal: {error:?}"
    );
}

#[test]
fn shaping_adds_nutation_and_gates_speeds() {
    // the model shaping, data-free: the constant stub dpsi lands on
    // the longitude (degrees, Euclidean modulo), and without
    // FLG_SPEED all three speeds are exactly zero.
    let dpsi_deg = 0.01_f64.to_degrees();
    let (lon_shifted, _, _, _, _, _) = mean_lunar_point_state(
        MEAN_NODE,
        2460000.5,
        FLG_SWIEPH | FLG_SPEED,
        LunarCoverage::Unenforced,
        &stub_dpsi,
        &stub_nut_rate,
    )
    .expect("mean node must serve");
    let (lon_plain, _, _, _, _, _) = mean_lunar_point_state(
        MEAN_NODE,
        2460000.5,
        FLG_SWIEPH | FLG_SPEED,
        LunarCoverage::Unenforced,
        &|_: f64| Ok(0.0),
        &stub_nut_rate,
    )
    .expect("mean node must serve");
    assert_eq!(lon_shifted, (lon_plain + dpsi_deg).rem_euclid(360.0));
    let speedy = mean_lunar_state_with_providers(
        MEAN_APOG,
        2460000.5,
        FLG_SWIEPH | FLG_SPEED,
        &stub_dpsi,
        &stub_nut_rate,
    )
    .expect("mean apogee must serve");
    assert_ne!(speedy[3], 0.0, "FLG_SPEED must carry a rate");
    let still = mean_lunar_state_with_providers(
        MEAN_APOG,
        2460000.5,
        FLG_SWIEPH,
        &stub_dpsi,
        &stub_nut_rate,
    )
    .expect("mean apogee must serve");
    assert_eq!((still[3], still[4], still[5]), (0.0, 0.0, 0.0));
    assert_eq!(still[0], speedy[0]);
    assert_eq!(still[1], speedy[1]);
}

#[test]
fn nutation_speed_stencil_uses_the_release_half_step() {
    // The dispatch drives the rate provider with the 0.001-day half step: a recording stub must observe exactly
    // that value on the FLG_SPEED mean path.
    use std::cell::Cell;
    let seen = Cell::new(-1.0_f64);
    let record = |_: f64, half_step: f64| {
        seen.set(half_step);
        Ok(0.0)
    };
    mean_lunar_state_with_providers(
        MEAN_NODE,
        2460000.5,
        FLG_SWIEPH | FLG_SPEED,
        &stub_dpsi,
        &record,
    )
    .expect("mean node must serve");
    assert_eq!(seen.get(), 0.001);
}

#[test]
fn non_finite_epoch_refuses_before_any_provider() {
    // Typed refusal, data-free: a NaN epoch never reaches the
    // nutation providers.
    let touched = std::cell::Cell::new(false);
    let spy = |_: f64| {
        touched.set(true);
        Ok(0.0)
    };
    let error = mean_lunar_state_with_providers(
        MEAN_NODE,
        f64::NAN,
        FLG_SWIEPH | FLG_SPEED,
        &spy,
        &stub_nut_rate,
    )
    .expect_err("NaN epoch must refuse");
    assert!(
        format!("{error:?}").contains("NonFiniteArgument"),
        "wrong refusal: {error:?}"
    );
    assert!(!touched.get(), "provider must not be contacted");
}

#[test]
fn helio_bary_sentinel_predicate_matches_release() {
    // the model public convention, data-free: the six lunar points
    // (hence both mean ids here) take the exact-zero sentinel under
    // HELCTR/BARYCTR and never otherwise; planets never do. The
    // adapter refuses these legs instead of surfacing zeros.
    for body in [MEAN_NODE, MEAN_APOG] {
        assert!(wants_lunar_helio_bary_zero(body, FLG_HELCTR));
        assert!(wants_lunar_helio_bary_zero(body, FLG_BARYCTR));
        assert!(!wants_lunar_helio_bary_zero(body, FLG_SWIEPH | FLG_SPEED));
    }
    assert!(!wants_lunar_helio_bary_zero(MARS, FLG_HELCTR));
    // The mean-ecliptic rule stays off on the adapter's flag sets:
    // the declination arm masks SIDEREAL out (never SIDEREAL +
    // EQUATORIAL together, which would stay mean-ecliptic), and no
    // NONUT/J2000 bits are ever set.
    assert!(!no_nutation_frame(FLG_SWIEPH | FLG_SPEED));
    assert!(!no_nutation_frame(FLG_SWIEPH | FLG_SPEED | FLG_EQUATORIAL));
    assert!(no_nutation_frame(
        FLG_SWIEPH | FLG_SPEED | FLG_SIDEREAL | FLG_EQUATORIAL
    ));
}
