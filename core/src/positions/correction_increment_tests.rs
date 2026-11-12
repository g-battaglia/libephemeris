use super::*;

#[test]
fn transverse_lorentz_increment_retains_the_sub_ulp_longitudinal_term() {
    let beta = 2.0_f64.powi(-30);
    let corrected = apply_aberration_with_correction(
        [1.0, 0.0, 0.0],
        [0.0, beta * C_LIGHT_AU_DAY, 0.0],
        1.0 / C_LIGHT_AU_DAY,
        C_LIGHT_AU_DAY,
    )
    .unwrap();
    let exact_x = -(beta * beta) / (1.0 + (1.0 - beta * beta).sqrt());
    assert_eq!(corrected.position[0], 1.0);
    assert!(corrected.increment[0] < 0.0);
    assert!((corrected.increment[0] - exact_x).abs() <= 8.0 * f64::EPSILON * exact_x.abs());
    assert!((corrected.increment[1] - beta).abs() <= 8.0 * f64::EPSILON * beta);
    assert_eq!(corrected.increment[2], 0.0);
}

#[test]
fn ppn_increment_obeys_distance_scaling_even_when_added_ray_rounds_it_away() {
    let apply = |scale: f64| {
        apply_gravitational_deflection_with_lookup_correction(
            [scale, 0.2 * scale, 0.3 * scale],
            [0.1 * scale, -0.2 * scale, 0.3 * scale],
            FRAME_J2000,
            scale / C_LIGHT_AU_DAY,
            |body, _| Ok((body == SUN).then_some([3.0 * scale, 2.0 * scale, scale])),
        )
        .unwrap()
    };
    let ordinary = apply(1.0);
    let remote = apply(1.0e8);
    for axis in 0..3 {
        // Point-mass deflection angle scales as inverse distance. Scaling
        // every position together leaves its length-valued increment fixed.
        assert!(
            (ordinary.increment[axis] - remote.increment[axis]).abs()
                <= 64.0 * f64::EPSILON * ordinary.increment[axis].abs()
        );
    }
    assert!(remote.increment.iter().any(|value| *value != 0.0));
    assert_eq!(remote.position, [1.0e8, 0.2e8, 0.3e8]);
}

#[test]
fn combined_correction_retains_exact_existing_position_and_separate_increment() {
    let ray = [1.0, 0.2, 0.3];
    let observer = [0.1, -0.2, 0.3];
    let velocity = [0.002, -0.01, 0.003];
    let deflector = |body, _| Ok([3.0, 2.0 + body as f64, 1.0]);
    let lt = vec3_dist(ray) / C_LIGHT_AU_DAY;
    let ppn =
        apply_gravitational_deflection_with_correction(ray, observer, FRAME_J2000, lt, deflector)
            .unwrap();
    let lorentz =
        apply_aberration_with_correction(ppn.position, velocity, lt, C_LIGHT_AU_DAY).unwrap();
    let combined = apply_deflection_and_aberration_with_vector_route(
        ray,
        velocity,
        observer,
        FRAME_J2000,
        lt,
        MARS,
        0,
        deflector,
        None,
        false,
    )
    .unwrap();
    assert_eq!(
        combined.pos().map(f64::to_bits),
        lorentz.position.map(f64::to_bits)
    );
    assert_eq!(
        combined.correction_increment(),
        std::array::from_fn(|i| ppn.increment[i] + lorentz.increment[i])
    );
}
