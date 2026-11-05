//! Physical conjunction, contact geometry and native error domains.
use super::*;
use std::cell::Cell;

#[test]
fn conjunction_finds_both_directions_for_a_linear_lunar_orbit() {
    for direction in [-1, 1] {
        let result = source_conjunction(
            &mut |_| Ok((10.0, 0.0)),
            &mut |jd| Ok(((jd - 100.0) * OCC_MOON_RATE_DEG_DAY, 0.0)),
            100.0,
            direction,
        )
        .unwrap()
        .unwrap();
        let expected = 100.0 + if direction > 0 { 10.0 } else { -350.0 } / OCC_MOON_RATE_DEG_DAY;
        assert!((result.0 - expected).abs() < 1e-12);
        assert_eq!((result.1, result.2), (0.0, 0.0));
    }
}

#[test]
fn first_data_boundary_exhausts_search_but_iterative_failure_propagates() {
    let boundary = Error::EphemerisRange {
        requested_jd: 1.0,
        start_jd: 2.0,
        end_jd: 3.0,
        body: Some(VENUS),
    };
    assert!(
        source_conjunction(
            &mut |_| Err(boundary.clone()),
            &mut |_| panic!("no lunar sample outside coverage"),
            1.0,
            1
        )
        .unwrap()
        .is_none()
    );
    let reads = Cell::new(0);
    let error = source_conjunction(
        &mut |_| {
            reads.set(reads.get() + 1);
            if reads.get() == 1 {
                Ok((13.0, 0.0))
            } else {
                Err(boundary.clone())
            }
        },
        &mut |_| Ok((0.0, 0.0)),
        1.0,
        1,
    )
    .unwrap_err();
    assert_eq!(reads.get(), 2);
    assert_eq!(error, boundary);
}

#[test]
fn named_star_outside_lunar_latitude_band_is_impossible() {
    let error = source_conjunction_with_star(
        &mut |_| Ok((12.0, 8.0)),
        &mut |_| panic!("impossible latitude is resolved before lunar acquisition"),
        2451545.0,
        1,
        Some("star"),
    )
    .unwrap_err();
    assert!(matches!(error, Error::InvalidEclipseSearch { .. }));
}

#[test]
fn minimizer_uses_the_native_finite_interval_domain() {
    let minimum = golden_section_min(
        &mut |jd| Ok((jd - 100.0).powi(2)),
        99.0,
        101.0,
        GOLDEN_MIN_TOLERANCE_DAYS,
        GOLDEN_SECTION_ITERATIONS,
    )
    .unwrap();
    assert!((minimum - 100.0).abs() < 1e-6);
    assert!(
        golden_section_min(
            &mut |_| Ok(0.0),
            f64::NAN,
            f64::INFINITY,
            GOLDEN_MIN_TOLERANCE_DAYS,
            GOLDEN_SECTION_ITERATIONS
        )
        .is_err()
    );
    assert!(
        source_conjunction(
            &mut |_| Ok((0.0, 0.0)),
            &mut |_| Ok((0.0, 0.0)),
            f64::NAN,
            1
        )
        .is_err()
    );
}

#[test]
fn finite_disc_residual_matches_axis_and_right_angle_geometry() {
    let lunar_distance = 0.0025;
    let lunar_radius = (ECL_RMOON_AU / lunar_distance).asin().to_degrees();
    let target_radius = 0.001f64.asin().to_degrees();
    let aligned =
        numeric_occult_separation_source([1.0, 0.0, 0.0], [lunar_distance, 0.0, 0.0], 0.001)
            .unwrap();
    assert!((aligned + lunar_radius + target_radius).abs() < 1e-12);
    let square =
        numeric_occult_separation_source([1.0, 0.0, 0.0], [0.0, lunar_distance, 0.0], 0.001)
            .unwrap();
    assert!((square - (90.0 - lunar_radius - target_radius)).abs() < 1e-12);
}

#[test]
fn impossible_distance_radius_and_overflow_have_typed_errors() {
    assert!(matches!(
        numeric_occult_separation_source([0.0; 3], [1.0, 0.0, 0.0], 0.0),
        Err(Error::RationalDivisionByZero)
    ));
    assert!(matches!(
        numeric_occult_separation_source([1.0, 0.0, 0.0], [1e-12, 0.0, 0.0], 1.0),
        Err(Error::MathDomain { .. })
    ));
    for value in [f64::MAX, -f64::MAX] {
        assert!(matches!(
            finite_square(value),
            Err(Error::NumericalOverflow { .. })
        ));
    }
    for value in [f64::INFINITY, f64::NEG_INFINITY, f64::NAN] {
        assert!(matches!(
            finite_square(value),
            Err(Error::NonFiniteArgument { .. })
        ));
    }
}
