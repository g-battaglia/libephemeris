use super::*;
use std::cell::RefCell;

#[test]
fn lunar_native_cartesian_order_and_speed_roles_are_distinct() {
    for speed in [false, true] {
        let reads = RefCell::new(Vec::new());
        let call = |jd: f64, body, flags| {
            reads.borrow_mut().push((jd.to_bits(), body, flags));
            Ok([body as f64, 2.0, 3.0, 4.0, 5.0, 6.0])
        };
        let flags = FLG_SWIEPH | FLG_SIDEREAL | FLG_TOPOCTR | FLG_SPEED;
        let pair = native_xyz_via(2451545.125, flags, speed, &call).unwrap();
        let expected = FLG_SWIEPH | FLG_EQUATORIAL | FLG_XYZ | if speed { FLG_SPEED } else { 0 };
        assert_eq!(
            *reads.borrow(),
            [
                (2451545.125_f64.to_bits(), MOON, expected),
                (2451545.125_f64.to_bits(), SUN, expected)
            ]
        );
        assert_eq!(pair, ((MOON as f64, 2.0, 3.0), (SUN as f64, 2.0, 3.0)));
    }
}

#[test]
fn lunar_native_first_position_failure_prevents_later_read() {
    let reads = RefCell::new(Vec::new());
    let call = |_: f64, body, _| {
        reads.borrow_mut().push(body);
        Err(Error::UnknownBody { body_id: body })
    };
    assert!(matches!(
        native_xyz_via(2451545.0, 2, true, &call),
        Err(Error::UnknownBody { body_id: MOON })
    ));
    assert_eq!(*reads.borrow(), [MOON]);
}
