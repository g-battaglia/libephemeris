//! Provider-free callback/order/ownership controls; no ephemeris or source data.
use super::*;
use std::cell::RefCell;

fn how(app: f64) -> SolHowCircumstances {
    solar::sol_how_core_via(
        (0.0, 0.0, 1.0),
        (0.0, 0.0, 0.0025),
        solar::ECL_RSUN_AU,
        true,
        false,
        (10.0, 1.0, app),
        0.0,
        (1.0, 2.0),
    )
    .unwrap()
}
fn geometry() -> WhereGeometry {
    solar::eclipse_where_core_via((0.0025, 0.0, 0.0), (1.0, 0.0, 0.0), 0.0, solar::ECL_RSUN_AU)
        .unwrap()
}

#[test]
fn how_local_then_global_original_ut_site_stripped_flags() {
    let log = RefCell::new(Vec::new());
    let site = [12.0, 30.0, 40.0];
    let mut local = |ut, gp, f| {
        log.borrow_mut().push(("local", ut, gp, f));
        Ok(how(1.0))
    };
    let mut global = |ut, f| {
        log.borrow_mut().push(("global", ut, [0.0; 3], f));
        Ok(geometry())
    };
    solar_how_via(
        2451545.0,
        site,
        FLG_SWIEPH | FLG_TOPOCTR,
        &mut local,
        &mut global,
    )
    .unwrap();
    assert_eq!(
        *log.borrow(),
        vec![
            ("local", 2451545.0, site, FLG_SWIEPH),
            ("global", 2451545.0, [0.0; 3], FLG_SWIEPH)
        ]
    );
}
#[test]
fn where_global_then_actual_center_with_zero_height() {
    let log = RefCell::new(Vec::new());
    let geom = geometry();
    let mut global = |ut, f| {
        log.borrow_mut().push(("global", ut, [0.0; 3], f));
        Ok(geom)
    };
    let mut local = |ut, gp, f| {
        log.borrow_mut().push(("local", ut, gp, f));
        Ok(how(-1.0))
    };
    let (bits, gp, attr) = solar_where_via(2451545.0, FLG_SWIEPH, &mut global, &mut local).unwrap();
    assert_eq!(bits, geom.retc());
    assert_eq!(gp[0], geom.lon_deg());
    assert_eq!(gp[1], geom.lat_deg());
    assert_eq!(
        *log.borrow(),
        vec![
            ("global", 2451545.0, [0.0; 3], FLG_SWIEPH),
            (
                "local",
                2451545.0,
                [geom.lon_deg(), geom.lat_deg(), 0.0],
                FLG_SWIEPH
            )
        ]
    );
    assert_eq!(attr[6], -1.0);
    assert_eq!(attr[3], geom.dcore()[0]);
}
#[test]
fn local_failure_prevents_global_read_and_keeps_owned_error() {
    let error = Error::UnknownBody { body_id: 123 };
    let mut local = |_, _, _| Err(error.clone());
    let mut global = |_, _| panic!("unused global");
    assert_eq!(
        solar_how_via(1.0, [0.0; 3], 0, &mut local, &mut global),
        Err(error)
    );
}
#[test]
fn where_failure_prevents_local_read_and_keeps_owned_error() {
    let error = Error::UnknownBody { body_id: 124 };
    let mut global = |_, _| Err(error.clone());
    let mut local = |_, _, _| panic!("unused local");
    assert_eq!(solar_where_via(1.0, 0, &mut global, &mut local), Err(error));
}
#[test]
fn no_visible_how_clears_only_the_public_slots() {
    let original = how(-1.0);
    let mut local = |_, _, _| Ok(original);
    let geom = geometry();
    let mut global = |_, _| Ok(geom);
    let (bits, attr) = solar_how_via(1.0, [0.0; 3], 0, &mut local, &mut global).unwrap();
    assert_eq!(bits, 0);
    for i in [0, 1, 2, 3, 8, 9, 10] {
        assert_eq!(attr[i], 0.0);
    }
    for i in [4, 5, 6, 7] {
        assert_eq!(attr[i], original.attr()[i]);
    }
}
#[test]
fn lazy_local_horizontal_error_precedes_saros_callback() {
    let error = Error::UnknownBody { body_id: 125 };
    let mut horizontal = || Err(error.clone());
    let mut saros = || panic!("unused Saros");
    let got = solar::sol_how_core_lazy(
        (0.0, 0.0, 1.0),
        (0.0, 0.0, 0.0025),
        solar::ECL_RSUN_AU,
        true,
        false,
        &mut horizontal,
        0.0,
        &mut saros,
    );
    assert_eq!(got, Err(error));
}
#[test]
fn lazy_local_callbacks_and_algebra_equal_literal_kernel() {
    let order = RefCell::new(Vec::new());
    let mut horizontal = || {
        order.borrow_mut().push("horizontal");
        Ok((10.0, 1.0, 1.0))
    };
    let mut saros = || {
        order.borrow_mut().push("Saros");
        Ok((1.0, 2.0))
    };
    let result = solar::sol_how_core_lazy(
        (0.0, 0.0, 1.0),
        (0.0, 0.0, 0.0025),
        solar::ECL_RSUN_AU,
        true,
        false,
        &mut horizontal,
        0.0,
        &mut saros,
    )
    .unwrap();
    assert_eq!(result, how(1.0));
    assert_eq!(*order.borrow(), vec!["horizontal", "Saros"]);
}
#[test]
fn lazy_where_retains_geometry_and_one_sidereal_call() {
    let calls = RefCell::new(0);
    let mut gast = || {
        *calls.borrow_mut() += 1;
        Ok(0.0)
    };
    let result = solar::eclipse_where_core_lazy(
        (0.0025, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        &mut gast,
        solar::ECL_RSUN_AU,
    )
    .unwrap();
    assert_eq!(result, geometry());
    assert_eq!(*calls.borrow(), 1);
}
#[test]
fn native_how_consumes_main_configuration_before_unused_clock() {
    let session = Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
    let error = Error::UnknownBody { body_id: 126 };
    let calls = RefCell::new(Vec::new());
    let resolve = |_: &EphemerisContext, body| {
        calls.borrow_mut().push(body);
        Err(error.clone())
    };
    let result = session.sol_eclipse_how_with_position_context(f64::NAN, [0.0; 3], 0, &resolve);
    assert_eq!(result, Err(error));
    assert_eq!(*calls.borrow(), vec![SUN]);
}
