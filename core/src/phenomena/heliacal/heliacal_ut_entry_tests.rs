use super::*;
use crate::constants::VENUS;

#[test]
fn parse_heliacal_object_forms() {
    // Integer forms admit every non-Sun id (range gates live
    // downstream in `validate_heliacal_body`).
    assert_eq!(parse_heliacal_object("3").unwrap(), VislimBody::Planet(3));
    assert_eq!(
        parse_heliacal_object("  4 ").unwrap(),
        VislimBody::Planet(4)
    );
    assert_eq!(parse_heliacal_object("+5").unwrap(), VislimBody::Planet(5));
    assert_eq!(
        parse_heliacal_object("1").unwrap(),
        VislimBody::Planet(MOON)
    );
    assert_eq!(parse_heliacal_object("99").unwrap(), VislimBody::Planet(99));
    // Planet names (case/space insensitive like the release).
    assert_eq!(
        parse_heliacal_object("venus").unwrap(),
        VislimBody::Planet(3)
    );
    assert_eq!(
        parse_heliacal_object("  Mars ").unwrap(),
        VislimBody::Planet(4)
    );
    assert_eq!(
        parse_heliacal_object("MOON").unwrap(),
        VislimBody::Planet(MOON)
    );
    // A catalog star resolves through the fuzzy project helper.
    let sirius = crate::stars::resolve_star_name("Sirius").unwrap();
    assert_eq!(
        parse_heliacal_object("Sirius").unwrap(),
        VislimBody::Star(sirius)
    );
    // Out-of-range integers refuse downstream as unknown bodies.
    let req = HeliacalSearchRequest::new(2460000.5, 99, HeliacalEvent::Rising, 0);
    assert!(matches!(
        req.unwrap_err(),
        Error::UnknownBody { body_id: 99 }
    ));
    // The Moon rising gate refuses downstream (evening-first admits).
    assert!(HeliacalSearchRequest::new(2460000.5, MOON, HeliacalEvent::Rising, 0).is_err());
    assert!(HeliacalSearchRequest::new(2460000.5, MOON, HeliacalEvent::EveningFirst, 0).is_ok());
    let _ = VENUS;
}

#[test]
fn parse_heliacal_object_refusals() {
    // Empty name refuses.
    assert!(parse_heliacal_object("").is_err());
    // The Sun is never valid, in either spelling.
    assert!(matches!(
        parse_heliacal_object("0"),
        Err(Error::InvalidBody { body_id: SUN, .. })
    ));
    assert!(matches!(
        parse_heliacal_object("Sun"),
        Err(Error::InvalidBody { body_id: SUN, .. })
    ));
    // Unresolvable names refuse as star-not-found.
    assert!(matches!(
        parse_heliacal_object("Not A Star"),
        Err(Error::StarNotFound { .. })
    ));
}

#[test]
fn window_margin_arms() {
    let served = [1.5, 20.0, 250.0, -12.0, 260.0, -30.0, 100.0, 1.0, 0.0, 0.0];
    assert_eq!(heliacal_window_margin(0, &served), 0.5);
    // Below-horizon sentinel flag forces the flat miss margin.
    assert_eq!(heliacal_window_margin(-2, &served), -99.0);
    // Daylight (Sun at/above the horizon) forces the flat miss margin.
    let mut day = served;
    day[3] = 0.0;
    assert_eq!(heliacal_window_margin(0, &day), -99.0);
    day[3] = 5.0;
    assert_eq!(heliacal_window_margin(0, &day), -99.0);
}

#[test]
fn window_finds_parabola_optimum_and_crossings() {
    // Analytic twilight peak: margin 1 at t0, zero at t0 ± w.
    const T0: f64 = 2460000.55;
    const W: f64 = 0.05;
    let margin = |t: f64| 1.0 - ((t - T0) / W).powi(2);
    let (start, opt, end) = heliacal_visibility_window(2460000.5, margin).unwrap();
    assert!((opt - T0).abs() < 2e-5, "opt {opt}");
    assert!((start - (T0 - W)).abs() < 5e-5, "start {start}");
    assert!((end - (T0 + W)).abs() < 5e-5, "end {end}");
    assert!(start <= opt && opt <= end);
}

#[test]
fn window_degenerate_when_never_visible() {
    // A flat miss margin reports the event instant thrice.
    let (start, opt, end) = heliacal_visibility_window(2460000.5, |_| -99.0).unwrap();
    assert_eq!((start, opt, end), (2460000.5, 2460000.5, 2460000.5));
    // A non-finite anchor refuses fail-closed.
    assert!(heliacal_visibility_window(f64::NAN, |_| 1.0).is_err());
}
