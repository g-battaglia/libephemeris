// SPDX-License-Identifier: AGPL-3.0
//! Focused tests for the parent position-routing module.

use super::*;
use std::cell::RefCell;

#[test]
fn geometric_vector_retains_original_parts_and_skips_apparent_supports() {
    let epoch = TwoPartEpoch {
        hi: 2451545.125,
        lo: 0.00075,
    };
    let seen = RefCell::new(Vec::new());
    let state = |id, time: TwoPartEpoch| {
        seen.borrow_mut()
            .push((id, time.hi.to_bits(), time.lo.to_bits()));
        Ok([4.0, 6.0, 8.0])
    };
    let retard =
        |_: TwoPartEpoch, _: f64| -> Result<[f64; 3]> { panic!("geometric target cannot retard") };
    let result = native_pheno_observe(
        &state,
        &retard,
        epoch,
        MARS,
        [1.0, 2.0, 3.0],
        [0.0; 3],
        None,
        true,
        true,
    )
    .unwrap();
    assert_eq!(result, [3.0, 4.0, 5.0]);
    assert_eq!(
        &*seen.borrow(),
        &[(MARS, epoch.hi.to_bits(), epoch.lo.to_bits())]
    );
}

#[test]
fn astrometry_retards_from_original_time_without_deflector_reads() {
    let epoch = TwoPartEpoch {
        hi: 2451545.125,
        lo: 0.00075,
    };
    let state = |id, time: TwoPartEpoch| {
        assert_eq!(id, MARS);
        assert_eq!(time.hi.to_bits(), epoch.hi.to_bits());
        assert_eq!(time.lo.to_bits(), epoch.lo.to_bits());
        Ok([1.0, 0.0, 0.0])
    };
    let delays = RefCell::new(Vec::new());
    let retard = |time: TwoPartEpoch, delay| {
        assert_eq!(time.hi.to_bits(), epoch.hi.to_bits());
        assert_eq!(time.lo.to_bits(), epoch.lo.to_bits());
        delays.borrow_mut().push(delay);
        Ok([1.0, 0.0, 0.0])
    };
    let result = native_pheno_observe(
        &state, &retard, epoch, MARS, [0.0; 3], [0.0; 3], None, false, false,
    )
    .unwrap();
    assert_eq!(result, [1.0, 0.0, 0.0]);
    assert_eq!(&*delays.borrow(), &[1.0 / FRAME_LIGHT_SPEED_AU_DAY]);
}

#[test]
fn apparent_deflector_miss_is_strict_and_preserves_source_order() {
    let epoch = TwoPartEpoch {
        hi: 2451545.125,
        lo: 0.00075,
    };
    let reads = RefCell::new(Vec::new());
    let state = |id, time: TwoPartEpoch| {
        reads.borrow_mut().push((id, time));
        match id {
            MARS => Ok([1.0, 0.0, 0.0]),
            SUN => Ok([0.2, 1.0, 0.0]),
            JUPITER => Err(Error::SpkUnregisteredBody { body: JUPITER }),
            _ => panic!("later support read before missing Jupiter"),
        }
    };
    let retard = |_: TwoPartEpoch, _: f64| Ok([1.0, 0.0, 0.0]);
    assert!(matches!(
        native_pheno_observe(
            &state, &retard, epoch, MARS, [0.0; 3], [0.0; 3], None, false, true
        ),
        Err(Error::SpkUnregisteredBody { body: JUPITER })
    ));
    let reads = reads.borrow();
    assert_eq!(
        reads.iter().map(|(id, _)| *id).collect::<Vec<_>>(),
        vec![MARS, SUN, SUN, JUPITER]
    );
    // Closest approach stays on original TT whole, with an offset in the
    // original small fraction. A fresh scalar epoch would lose that pair.
    assert_eq!(reads[2].1.hi.to_bits(), epoch.hi.to_bits());
    assert!(reads[2].1.lo < epoch.lo);
    assert_eq!(reads[3].1.lo.to_bits(), epoch.lo.to_bits());
}

#[test]
fn native_vector_nonfinite_epoch_refuses_before_acquisition() {
    let state = |_: i32, _: TwoPartEpoch| -> Result<[f64; 3]> { panic!("invalid clock read") };
    let retard = |_: TwoPartEpoch, _: f64| -> Result<[f64; 3]> { panic!("invalid retardation") };
    assert!(matches!(
        native_pheno_observe(
            &state,
            &retard,
            TwoPartEpoch {
                hi: f64::NAN,
                lo: 0.0
            },
            MARS,
            [0.0; 3],
            [0.0; 3],
            None,
            false,
            false
        ),
        Err(Error::NonFiniteArgument { .. })
    ));
}
