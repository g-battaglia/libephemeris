use super::*;
use std::cell::RefCell;

#[test]
fn lunar_local_moon_lower_limb_rise_then_set_same_inputs() {
    let reads = RefCell::new(Vec::new());
    let call = |jd: f64, rsmi, site, flags| {
        reads.borrow_mut().push((jd.to_bits(), rsmi, site, flags));
        let mut times = [0.0; 10];
        times[0] = 100.0;
        Ok((-2, times))
    };
    let site = [10.0, 20.0, 300.0];
    let result = moon_crossings_via(99.699, site, 2, &call).unwrap();
    assert_eq!(
        *reads.borrow(),
        [
            (99.699_f64.to_bits(), CALC_RISE | BIT_DISC_BOTTOM, site, 2),
            (99.699_f64.to_bits(), CALC_SET | BIT_DISC_BOTTOM, site, 2),
        ]
    );
    assert_eq!(result.rise_rc, -2);
    assert_eq!(result.set_rc, -2);
}
#[test]
fn lunar_local_first_rise_failure_prevents_set() {
    let reads = RefCell::new(0);
    let error = Error::UnknownBody { body_id: MOON };
    let call = |_: f64, _, _, _| {
        *reads.borrow_mut() += 1;
        Err(error.clone())
    };
    assert_eq!(
        moon_crossings_via(100.0, [0.0; 3], 2, &call).unwrap_err(),
        error
    );
    assert_eq!(*reads.borrow(), 1);
}
#[test]
fn lunar_local_exhaustion_retains_the_native_search_domain() {
    let result = local_exhausted(2451545.0, [1.0, 2.0, 3.0]);
    assert!(matches!(
        result.primary_error(),
        Error::EclipseSearchExhausted {
            entry: "lun_eclipse_when_loc",
            ..
        }
    ));
}
