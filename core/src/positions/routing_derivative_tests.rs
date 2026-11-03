use super::*;

#[test]
fn correction_partials_resolve_opposing_owned_clock_domains_without_crossing() {
    for (ephemeris_support, site_support) in [
        ((-1.0, 0.0), (0.0, 1.0)),
        ((-1.0, 2.0_f64.powi(-34)), (0.0, 1.0)),
    ] {
        let sample = |ephemeris: f64, site: f64| {
            assert!((ephemeris_support.0..=ephemeris_support.1).contains(&ephemeris));
            assert!((site_support.0..=site_support.1).contains(&site));
            // Every noncentral sample must vary exactly one physical clock.
            // The constant is large relative to the tiny common window, while
            // the cubic has an exact fourth-order directional derivative.
            assert!(ephemeris == 0.0 || site == 0.0);
            Ok(Some([
                1.0 + 0.002 * ephemeris - 0.003 * site + ephemeris.powi(3),
                -2.0 + 0.004 * ephemeris + 0.005 * site + site.powi(3),
                0.5 + 0.006 * ephemeris - 0.007 * site + ephemeris * site,
            ]))
        };
        let actual = correction_rates_from_independent_clocks(
            &sample,
            ephemeris_support,
            Some(site_support),
        )
        .unwrap();
        let expected = [-0.001, 0.009, -0.001];
        // Bound weighted binary64 subtraction at the original dyadic scale.
        let bound = 256.0 * f64::EPSILON * 2.0 * 64.0;
        for i in 0..3 {
            assert!((actual[i] - expected[i]).abs() <= bound);
        }
    }
}

#[test]
fn rotating_site_jacobian_resolves_diurnal_motion_on_one_sided_clock_pieces() {
    use reference_frames::*;
    let itrs = [3.0e-5, -1.0e-5, 2.0e-5];
    let omega = FRAME_TAU * 1.0027378119113546;
    let precession = 3.0e-7;
    for phase in [0.0, 0.7, 3.1, 6.2] {
        for support in [(-1.0, 1.0), (0.0, 1.0), (-1.0, 0.0)] {
            let orientation = |offset: f64| {
                if !(support.0..=support.1).contains(&offset) {
                    return Ok(None);
                }
                // A smooth forward clock with curvature and a slow moving
                // equinox have exact central derivatives. The diurnal sine
                // must not be approximated with a coarse one-sided stencil.
                let angle = phase + omega * (offset + 0.001 * offset.powi(2));
                Ok(Some((
                    rot_x_neg(-precession * offset),
                    (angle / FRAME_TAU * 24.0).rem_euclid(24.0),
                )))
            };
            let actual = rotating_site_velocity_tt(itrs, &orientation, support).unwrap();
            let q = matvec3(rot_z(phase), itrs).unwrap();
            let expected = [
                -omega * q[1],
                omega * q[0] + precession * q[2],
                -precession * q[1],
            ];
            // Matrix/angle roundoff amplified by the fixed 1/64-day scale;
            // the analytic synthetic trajectory has no unknown model error.
            let bound = 256.0 * f64::EPSILON * (1.0 + phase) * vec3_dist(itrs) * 64.0;
            for i in 0..3 {
                assert!(
                    (actual[i] - expected[i]).abs() <= bound,
                    "phase={phase} support={support:?} axis={i}: actual={} exact={} bound={bound}",
                    actual[i],
                    expected[i]
                );
            }
        }
    }
}

#[test]
fn topocentric_light_path_velocity_obeys_the_retarded_ray_jacobian() {
    let whole = 2415020.5;
    let site = (12.5, 41.9, 150.0);
    let target_velocity = [0.0001, 0.0005, -0.00003];
    let state = |body, w: f64, f: f64| {
        if body == EARTH {
            return Ok(([0.0; 3], [0.0; 3]));
        }
        assert_eq!(body, MOON);
        let offset = (w - whole) + f;
        Ok((
            std::array::from_fn(|i| [0.003, 0.0007, 0.0002][i] + target_velocity[i] * offset),
            target_velocity,
        ))
    };
    let observer_velocity = local_jpl_topo_velocity_tt(whole, 0.0, site).unwrap();
    let flags = FLG_TOPOCTR | FLG_NOGDEFL | FLG_NOABERR | FLG_J2000 | FLG_EQUATORIAL;
    let actual = vector_place_from_state_routes_with_policy(
        &state,
        None,
        whole,
        0.0,
        MOON,
        flags,
        Some(site),
        false,
        Some(1.0),
        Some(observer_velocity),
        None,
    )
    .unwrap();
    let radius = vec3_dist(actual.apparent_ray);
    let unit = actual.apparent_ray.map(|value| value / radius);
    let c = FRAME_LIGHT_SPEED_AU_DAY;
    let emission_rate = (c + compensated_dot3(unit, observer_velocity))
        / (c + compensated_dot3(unit, target_velocity));
    let expected = std::array::from_fn::<_, 3, _>(|i| {
        target_velocity[i] * emission_rate - observer_velocity[i]
    });
    for (i, value) in expected.into_iter().enumerate() {
        assert!(
            (actual.ray_velocity_tt[i] - value).abs() < 1e-15,
            "axis={i} actual={} implicit={value}",
            actual.ray_velocity_tt[i]
        );
    }
}

#[test]
fn affine_tiny_angular_rate_does_not_amplify_one_second_coordinate_roundoff() {
    let fraction = 0.375;
    let place = |whole: f64, part: f64| {
        assert_eq!(whole, FRAME_J2000);
        let t = part - fraction;
        Ok([22.0 + 0.000003 * t, 7.0 - 0.000004 * t, 1.0])
    };
    let result =
        vector_position_from_places(&place, FRAME_J2000, fraction, MARS, FLG_SPEED).unwrap();
    assert!((result[3] - 0.000003).abs() <= 0.000003 * 1e-6);
    assert!((result[4] + 0.000004).abs() <= 0.000004 * 1e-6);
}

#[test]
fn nearly_stationary_output_rate_includes_the_actual_moving_frame() {
    let fraction = 0.125;
    let place = |whole: f64, part: f64| {
        assert_eq!(whole, FRAME_J2000);
        let t = part - fraction;
        // Intrinsic drift and an independently rotating output frame almost
        // cancel. Differentiating frozen-frame states would lose this rate.
        let longitude = (22.0 + 0.050003 * t).to_radians();
        let latitude = (7.0 - 0.000004 * t).to_radians();
        let theta = (-0.05 * t).to_radians();
        let (sl, cl) = longitude.sin_cos();
        let (sb, cb) = latitude.sin_cos();
        let (st, ct) = theta.sin_cos();
        cartesian_to_spherical_deg([cb * (ct * cl - st * sl), cb * (st * cl + ct * sl), sb])
    };
    let result =
        vector_position_from_places(&place, FRAME_J2000, fraction, MARS, FLG_SPEED).unwrap();
    assert!((result[3] - (0.050003 - 0.05)).abs() <= 0.000003 * 1e-6);
    assert!((result[4] + 0.000004).abs() <= 0.000004 * 1e-6);
}

#[test]
fn actual_retarded_and_deflector_record_domains_remain_separate() {
    let domains = LocalJplReadDomains::default();
    domains.ordinary.borrow_mut().push(LocalJplReadDomain {
        body: 14,
        epoch: TwoPartEpoch {
            hi: FRAME_J2000,
            lo: 0.5,
        },
        offsets: (0.0, 1.0),
    });
    domains.ordinary.borrow_mut().push(LocalJplReadDomain {
        body: 5,
        epoch: TwoPartEpoch {
            hi: FRAME_J2000,
            lo: 0.49,
        },
        offsets: (-0.5, 0.001),
    });
    domains.retarded.borrow_mut().push(LocalJplReadDomain {
        body: 4,
        epoch: TwoPartEpoch {
            hi: FRAME_J2000,
            lo: 0.48,
        },
        offsets: (-0.5, 0.002),
    });
    let admission = LocalJplReadAdmission::new(&domains);
    assert!(
        admission
            .check(
                14,
                TwoPartEpoch {
                    hi: FRAME_J2000,
                    lo: 0.5005
                },
                false
            )
            .is_ok()
    );
    assert!(!admission.complete());
    let new_optional = LocalJplReadAdmission::new(&domains);
    assert!(
        new_optional
            .check(
                0,
                TwoPartEpoch {
                    hi: FRAME_J2000,
                    lo: 0.51
                },
                false
            )
            .is_err()
    );
    assert!(!new_optional.complete());
    let dropped_role = LocalJplReadAdmission::new(&domains);
    assert!(!dropped_role.complete());
    assert!(
        admission
            .check(
                4,
                TwoPartEpoch {
                    hi: FRAME_J2000,
                    lo: 0.4804
                },
                true
            )
            .is_ok()
    );
    assert!(
        admission
            .check(
                5,
                TwoPartEpoch {
                    hi: FRAME_J2000,
                    lo: 0.4915
                },
                false
            )
            .is_err()
    );
    assert!(admission.refused.get());
}

#[test]
fn stable_unavailable_deflector_does_not_consume_the_next_successful_role() {
    let domains = LocalJplReadDomains::default();
    domains.ordinary.borrow_mut().push(LocalJplReadDomain {
        body: 5,
        epoch: TwoPartEpoch {
            hi: FRAME_J2000,
            lo: 0.5,
        },
        offsets: (-1.0, 1.0),
    });
    let admission = LocalJplReadAdmission::new(&domains);
    admission.failed_read(
        0,
        TwoPartEpoch {
            hi: FRAME_J2000,
            lo: 0.51,
        },
        false,
    );
    assert!(!admission.refused.get());
    assert_eq!(admission.ordinary.get(), 0);
    assert!(
        admission
            .check(
                5,
                TwoPartEpoch {
                    hi: FRAME_J2000,
                    lo: 0.51
                },
                false
            )
            .is_ok()
    );
    assert!(admission.complete());
    let unavailable = LocalJplReadAdmission::new(&domains);
    unavailable.failed_read(
        5,
        TwoPartEpoch {
            hi: FRAME_J2000,
            lo: 2.0,
        },
        false,
    );
    assert!(unavailable.refused.get());
}

#[test]
fn shorter_light_time_iteration_still_requires_the_final_central_target_record() {
    let domains = LocalJplReadDomains::default();
    domains.retarded.borrow_mut().push(LocalJplReadDomain {
        body: 4,
        epoch: TwoPartEpoch {
            hi: FRAME_J2000,
            lo: 0.49,
        },
        offsets: (-0.49, 0.01_f64.next_down()),
    });
    domains.retarded.borrow_mut().push(LocalJplReadDomain {
        body: 4,
        epoch: TwoPartEpoch {
            hi: FRAME_J2000,
            lo: 0.51,
        },
        offsets: (-0.01, 0.49),
    });
    let earlier_record = LocalJplReadAdmission::new(&domains);
    assert!(
        earlier_record
            .check(
                4,
                TwoPartEpoch {
                    hi: FRAME_J2000,
                    lo: 0.49
                },
                true
            )
            .is_ok()
    );
    assert!(!earlier_record.complete());

    // Fewer fixed-point iterations are valid when the final emitted target
    // remains in the same owned polynomial as the central converged target.
    domains.retarded.borrow_mut()[0].offsets = (-0.49, 0.51);
    let same_final_record = LocalJplReadAdmission::new(&domains);
    assert!(
        same_final_record
            .check(
                4,
                TwoPartEpoch {
                    hi: FRAME_J2000,
                    lo: 0.52
                },
                true
            )
            .is_ok()
    );
    assert!(same_final_record.complete());
}
