//! Native local-solar numerical domains and exhausted-search behavior.
use super::*;

#[test]
fn global_candidate_policy_uses_native_error_categories() {
    assert!(caught_global(&Error::NetworkSealed {
        purpose: "ephemeris download".into()
    }));
    assert!(caught_global(&Error::UnknownBody { body_id: 777 }));
    assert!(!caught_global(&Error::MathDomain { value: 2.0 }));
    assert!(!caught_global(&Error::NumericalOverflow {
        operation: "distance",
        message: "distance overflowed".into()
    }));
}

#[test]
fn exhausted_search_retains_an_actual_native_cause() {
    let cause = Error::UnknownBody { body_id: 777 };
    let error = local_exhausted(2451545.0, [12.0, 34.0, 0.0], Some(cause.clone()));
    assert!(matches!(
        error.primary_error(),
        Error::EclipseSearchExhausted {
            entry: "sol_eclipse_when_loc",
            ..
        }
    ));
    let source = std::error::Error::source(&error)
        .unwrap()
        .downcast_ref::<Error>()
        .unwrap();
    assert_eq!(source, &cause);
    assert!(
        std::error::Error::source(&local_exhausted(2451545.0, [12.0, 34.0, 0.0], None)).is_none()
    );
}

#[test]
fn angular_radius_uses_the_real_arcsine_domain() {
    assert_eq!(radius_asin(0.0).unwrap(), 0.0);
    assert!((radius_asin(0.5).unwrap() - 30.0).abs() < 1e-13);
    for ratio in [-2.0, 2.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(radius_asin(ratio), Err(Error::MathDomain { .. })));
    }
}

struct SealedGlobal {
    calls: std::cell::RefCell<Vec<f64>>,
    error: Error,
}

impl crate::events::SolarLocalEphemeris for SealedGlobal {
    fn next_global_forward(&self, jd: f64) -> Result<(i32, [f64; 10]), Error> {
        self.calls.borrow_mut().push(jd);
        Err(self.error.clone())
    }
    fn separation_deg(&self, _: f64) -> Result<f64, Error> {
        unreachable!()
    }
    fn radii_sep_deg(&self, _: f64, _: bool) -> Result<(f64, f64, f64), Error> {
        unreachable!()
    }
    fn sun_altaz_deg(&self, _: f64) -> Result<(f64, f64, f64), Error> {
        unreachable!()
    }
    fn topo_diameters_deg(&self, _: f64) -> Result<(f64, f64), Error> {
        unreachable!()
    }
    fn sun_horizon_crossings(&self, _: f64) -> Result<crate::events::BodyHorizonCrossings, Error> {
        unreachable!()
    }
    fn how_core(&self, _: f64) -> Result<(i32, [f64; 20]), Error> {
        unreachable!()
    }
    fn where_core(&self, _: f64) -> Result<(i32, [f64; 10]), Error> {
        unreachable!()
    }
}

#[test]
fn forward_global_failure_returns_exhaustion_with_its_native_cause() {
    let marker = Error::NetworkSealed {
        purpose: "ephemeris download".into(),
    };
    let eph = SealedGlobal {
        calls: Default::default(),
        error: marker.clone(),
    };
    let error = crate::events::sol_eclipse_when_loc_source_via(
        &eph,
        2451545.0,
        false,
        &caught_global,
        &|active| local_exhausted(2451545.0, [12.0, 34.0, 0.0], active),
    )
    .unwrap_err();
    assert_eq!(*eph.calls.borrow(), vec![2451545.0 - 0.35]);
    assert!(matches!(
        error.primary_error(),
        Error::EclipseSearchExhausted {
            entry: "sol_eclipse_when_loc",
            ..
        }
    ));
    assert_eq!(
        std::error::Error::source(&error)
            .unwrap()
            .downcast_ref::<Error>(),
        Some(&marker)
    );
}

#[test]
fn backward_search_remains_bounded_and_monotonic() {
    let eph = SealedGlobal {
        calls: Default::default(),
        error: Error::NetworkSealed {
            purpose: "ephemeris download".into(),
        },
    };
    let error = crate::events::sol_eclipse_when_loc_source_via(
        &eph,
        2451545.0,
        true,
        &caught_global,
        &|active| local_exhausted(2451545.0, [12.0, 34.0, 0.0], active),
    )
    .unwrap_err();
    assert_eq!(eph.calls.borrow().len(), 120);
    assert!(eph.calls.borrow().windows(2).all(|pair| pair[1] < pair[0]));
    assert!(matches!(
        error.primary_error(),
        Error::EclipseSearchExhausted {
            entry: "sol_eclipse_when_loc",
            ..
        }
    ));
    assert!(std::error::Error::source(&error).is_none());
}
