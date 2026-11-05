// SPDX-License-Identifier: AGPL-3.0
//! Focused tests for the parent position-routing module.

use super::*;
use std::cell::RefCell;

#[test]
fn astropy_core_time_shared_apparent_tail_preserves_support_order_and_error() {
    let calls = RefCell::new(Vec::new());
    let epoch = TwoPartEpoch {
        hi: 2451545.0,
        lo: 0.125,
    };
    let state_at = |body, time: TwoPartEpoch| {
        calls
            .borrow_mut()
            .push((body, time.hi.to_bits(), time.lo.to_bits()));
        if body == SUN {
            Ok([10.0, 0.0, 0.0])
        } else {
            Err(Error::SpkUnregisteredBody { body })
        }
    };
    assert_eq!(
        native_apparent_ray(
            &state_at,
            epoch,
            [0.0, 1.0, 0.0],
            1.0 / FRAME_LIGHT_SPEED_AU_DAY,
            [0.0; 3],
            [0.0; 3],
            None
        ),
        Err(Error::SpkUnregisteredBody { body: JUPITER })
    );
    let ids: Vec<_> = calls.borrow().iter().map(|entry| entry.0).collect();
    assert_eq!(ids, vec![SUN, SUN, JUPITER]);
    assert_eq!(
        calls.borrow()[0],
        (SUN, epoch.hi.to_bits(), epoch.lo.to_bits())
    );
    assert_eq!(
        calls.borrow()[2],
        (JUPITER, epoch.hi.to_bits(), epoch.lo.to_bits())
    );
}
