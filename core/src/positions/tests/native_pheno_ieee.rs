// SPDX-License-Identifier: AGPL-3.0
//! Focused tests for the parent position-routing module.

use super::*;

#[test]
fn finite_singular_surface_and_superluminal_velocity_flow_to_sun_tuple() {
    let epoch = TwoPartEpoch {
        hi: 2451545.125,
        lo: 0.0,
    };
    let state = |body, _: TwoPartEpoch| {
        Ok(match body {
            SUN => [1.0, 0.0, 0.0],
            JUPITER => [0.2, 1.0, 0.0],
            SATURN => [0.3, 0.0, 1.0],
            _ => panic!("unexpected support"),
        })
    };
    let retarded = |_: TwoPartEpoch, _: f64| Ok([1.0, 0.0, 0.0]);
    for (surface, velocity) in [
        (Some([0.0; 3]), [0.0; 3]),
        (None, [2.0 * FRAME_LIGHT_SPEED_AU_DAY, 0.0, 0.0]),
    ] {
        let ray = native_pheno_observe(
            &state, &retarded, epoch, SUN, [0.0; 3], velocity, surface, false, true,
        )
        .unwrap();
        assert!(ray.into_iter().all(f64::is_nan));
        let tuple = crate::photometry::pheno_from_native_vectors(
            SUN,
            ray,
            [0.0; 3],
            [0.0; 3],
            epoch.value(),
            (0.0, 0.0),
            (0.0, 0.0),
            0.0,
        )
        .unwrap();
        assert!(tuple[3].is_nan());
        assert_eq!(tuple[4].to_bits(), (-26.86_f64).to_bits());
        assert!(
            tuple[..3]
                .iter()
                .chain(tuple[5..].iter())
                .all(|v| v.to_bits() == 0)
        );
    }
}

#[test]
fn raw_nonfinite_state_still_refuses_before_retardation_or_reduction() {
    let state = |_: i32, _: TwoPartEpoch| Ok([f64::INFINITY, 0.0, 0.0]);
    let retarded =
        |_: TwoPartEpoch, _: f64| -> Result<[f64; 3]> { panic!("invalid raw state cannot retard") };
    assert!(matches!(
        native_pheno_observe(
            &state,
            &retarded,
            TwoPartEpoch {
                hi: 2451545.125,
                lo: 0.0
            },
            SUN,
            [0.0; 3],
            [0.0; 3],
            None,
            false,
            true
        ),
        Err(Error::NonFiniteArgument { .. })
    ));
}
