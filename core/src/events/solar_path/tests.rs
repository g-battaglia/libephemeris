use super::*;
use std::cell::RefCell;
struct Width {
    trace: RefCell<Vec<&'static str>>,
    central: bool,
    fail_alt: bool,
}
impl SolarPathWidthSource for Width {
    type Sample = ();
    fn width_where(&self, _: f64) -> Result<(i32, f64, f64)> {
        self.trace.borrow_mut().push("where");
        Ok((if self.central { ECL_CENTRAL } else { 0 }, 10.0, 20.0))
    }
    fn width_acquire(&self, _: f64, _: f64, _: f64) -> Result<Option<()>> {
        self.trace.borrow_mut().push("acquire");
        Ok(Some(()))
    }
    fn width_altitude(&self, _: &()) -> Result<f64> {
        self.trace.borrow_mut().push("altitude");
        if self.fail_alt {
            Err(Error::UnknownBody { body_id: 999 })
        } else {
            Ok(-1.0)
        }
    }
    fn width_distances(&self, _: &()) -> (f64, f64) {
        panic!("below horizon must not consume distances")
    }
    fn width_separation(&self, _: &()) -> f64 {
        panic!("below horizon must not consume separation")
    }
}
#[test]
fn noncentral_where_skips_observer() {
    let e = Width {
        trace: Default::default(),
        central: false,
        fail_alt: false,
    };
    assert_eq!(
        solar_path_width_source_via(&e, 2451545.0, None, None).unwrap(),
        0.0
    );
    assert_eq!(*e.trace.borrow(), ["where"]);
}
#[test]
fn explicit_site_ignores_where_and_horizon_skips_distances() {
    let e = Width {
        trace: Default::default(),
        central: false,
        fail_alt: false,
    };
    assert_eq!(
        solar_path_width_source_via(&e, 1.0, Some(2.0), Some(3.0)).unwrap(),
        0.0
    );
    assert_eq!(*e.trace.borrow(), ["acquire", "altitude"]);
}
#[test]
fn half_site_uses_where() {
    let e = Width {
        trace: Default::default(),
        central: true,
        fail_alt: false,
    };
    assert_eq!(
        solar_path_width_source_via(&e, 1.0, Some(2.0), None).unwrap(),
        0.0
    );
    assert_eq!(*e.trace.borrow(), ["where", "acquire", "altitude"]);
}
#[test]
fn altitude_failure_is_not_caught_as_zero() {
    let e = Width {
        trace: Default::default(),
        central: true,
        fail_alt: true,
    };
    assert_eq!(
        solar_path_width_source_via(&e, 1.0, None, None),
        Err(Error::UnknownBody { body_id: 999 })
    );
}
#[test]
fn undefined_declination_cannot_form_a_finite_limit_latitude() {
    let p = solar_limit_point_source(0.1, 0.2, f64::NAN, 10.0, 0.01, false)
        .unwrap()
        .unwrap();
    assert!(p.0.is_nan());
}
#[test]
fn limit_north_and_south_reuse_finite_geodetic_kernel() {
    for south in [false, true] {
        let old = eclipse_limit_point_deg(0.1, 0.2, 7.0, 100.0, -0.01, south);
        let new = solar_limit_point_source(0.1, 0.2, 7.0, 100.0, -0.01, south).unwrap();
        assert_eq!(old, new);
    }
}

#[test]
fn infinite_declination_has_a_native_domain_error() {
    let error = solar_limit_point_source(0.1, 0.2, f64::INFINITY, 10.0, 0.01, false).unwrap_err();
    assert!(matches!(error.primary_error(),Error::MathDomain { value } if value.is_infinite()));
}
