use super::*;
use crate::constants::VENUS;
use std::rc::Rc;

/// Boxed word-provider set for tests (single concrete instantiation of
/// the sealed provider; no I/O, no ephemeris, no network).
type TestSealed<'a> = SealedHeliacalEphemeris<
    'a,
    (),
    Box<dyn Fn(f64) -> Result<f64>>,
    Box<dyn Fn(f64) -> Result<(f64, f64)>>,
    Box<dyn Fn(f64) -> (f64, f64)>,
    Box<dyn Fn(f64) -> (f64, f64)>,
    Box<dyn Fn(f64) -> f64>,
    Box<dyn Fn(f64) -> f64>,
    Box<dyn Fn(f64) -> (f64, f64, f64)>,
    Box<dyn Fn(f64) -> f64>,
>;

const JD0: f64 = 2460000.5;

// Clock and nutation fixture callbacks remain explicitly separate.
#[allow(clippy::type_complexity)]
fn healthy_epoch() -> (
    Box<dyn Fn(f64) -> Result<f64>>,
    Box<dyn Fn(f64) -> Result<(f64, f64)>>,
) {
    (Box::new(|_| Ok(0.0)), Box::new(|_| Ok((0.0, 0.0))))
}

fn const_altaz(v: (f64, f64)) -> Box<dyn Fn(f64) -> (f64, f64)> {
    Box::new(move |_| v)
}

fn const_word(v: f64) -> Box<dyn Fn(f64) -> f64> {
    Box::new(move |_| v)
}

fn const_moon(v: (f64, f64, f64)) -> Box<dyn Fn(f64) -> (f64, f64, f64)> {
    Box::new(move |_| v)
}

/// Build a sealed provider over scripted words; every provider call is
/// logged so acquisition order is assertable.
#[allow(clippy::too_many_arguments)]
fn build_sealed<'a>(
    serving: &'a (),
    request: HeliacalSearchRequest,
    log: Rc<RefCell<Vec<&'static str>>>,
    sun: impl Fn(f64) -> (f64, f64) + 'static,
    body: impl Fn(f64) -> (f64, f64) + 'static,
    elong: impl Fn(f64) -> f64 + 'static,
    mag: impl Fn(f64) -> f64 + 'static,
    moon: impl Fn(f64) -> (f64, f64, f64) + 'static,
    londiff: impl Fn(f64) -> f64 + 'static,
    delta_t: Box<dyn Fn(f64) -> Result<f64>>,
    nutation: Box<dyn Fn(f64) -> Result<(f64, f64)>>,
) -> TestSealed<'a> {
    let wrap = |name: &'static str, log: Rc<RefCell<Vec<&'static str>>>| {
        move || log.borrow_mut().push(name)
    };
    let sun_log = wrap("sun", log.clone());
    let body_log = wrap("body", log.clone());
    let elong_log = wrap("elong", log.clone());
    let mag_log = wrap("mag", log.clone());
    let moon_log = wrap("moon", log.clone());
    let lon_log = wrap("lon", log.clone());
    let dt_log = wrap("dt", log.clone());
    let nut_log = wrap("nut", log);
    let sun_word: Box<dyn Fn(f64) -> (f64, f64)> = Box::new(move |jd| {
        sun_log();
        sun(jd)
    });
    let body_word: Box<dyn Fn(f64) -> (f64, f64)> = Box::new(move |jd| {
        body_log();
        body(jd)
    });
    let elong_word: Box<dyn Fn(f64) -> f64> = Box::new(move |jd| {
        elong_log();
        elong(jd)
    });
    let mag_word: Box<dyn Fn(f64) -> f64> = Box::new(move |jd| {
        mag_log();
        mag(jd)
    });
    let moon_word: Box<dyn Fn(f64) -> (f64, f64, f64)> = Box::new(move |jd| {
        moon_log();
        moon(jd)
    });
    let lon_word: Box<dyn Fn(f64) -> f64> = Box::new(move |jd| {
        lon_log();
        londiff(jd)
    });
    let dt_word: Box<dyn Fn(f64) -> Result<f64>> = Box::new(move |jd| {
        dt_log();
        delta_t(jd)
    });
    let nut_word: Box<dyn Fn(f64) -> Result<(f64, f64)>> = Box::new(move |jd| {
        nut_log();
        nutation(jd)
    });
    SealedHeliacalEphemeris::new(
        serving,
        request,
        12.5,
        41.9,
        0.0,
        1013.25,
        15.0,
        0.5,
        0.0,
        &[],
        dt_word,
        nut_word,
        sun_word,
        body_word,
        elong_word,
        mag_word,
        moon_word,
        lon_word,
    )
    .unwrap()
}

/// Rising script: Sun climbs through twilight every morning, the body
/// steps up at `t0`, elongation crosses the conjunction gap at day 5.
/// Note the Sun word keys on the day fraction (`rem_euclid`), so a
/// refinement step must sit within the +/-0.1-day refinement window of
/// the first visible grid sample (grid hours ~19-25 UT sit at
/// fractions ~0.29-0.54 of the following JD number).
// The fixture preserves independent Sun and target altitude callback types.
#[allow(clippy::type_complexity)]
fn rising_script(t0: f64) -> (impl Fn(f64) -> (f64, f64), impl Fn(f64) -> (f64, f64)) {
    let sun = move |jd: f64| {
        let frac = jd.rem_euclid(1.0);
        (-20.0 + 30.0 * frac, 0.0)
    };
    let body = move |jd: f64| {
        if jd >= t0 { (10.0, 60.0) } else { (-5.0, 60.0) }
    };
    (sun, body)
}

/// Falling Sun through twilight every evening (for evening-polarity
/// twilight tests; the rising script never yields an evening centre).
fn falling_sun() -> impl Fn(f64) -> (f64, f64) {
    move |jd: f64| {
        let frac = jd.rem_euclid(1.0);
        (10.0 - 30.0 * frac, 0.0)
    }
}

fn rising_request() -> HeliacalSearchRequest {
    HeliacalSearchRequest::new(JD0, VENUS, HeliacalEvent::Rising, 0).unwrap()
}

#[test]
fn constructor_builds_entry_model_and_echoes_site() {
    let unit = ();
    let log = Rc::new(RefCell::new(Vec::new()));
    let (dt, nut) = healthy_epoch();
    let (sun, body) = rising_script(JD0 + 5.25);
    let provider = build_sealed(
        &unit,
        rising_request(),
        log,
        sun,
        body,
        |_| 20.0,
        |_| -4.0,
        |_| (-90.0, 0.0, 180.0),
        |_| -15.0,
        dt,
        nut,
    );
    // Site words echo through the accessors.
    assert_eq!(provider.lon_deg(), 12.5);
    assert_eq!(provider.lat_deg(), 41.9);
    assert_eq!(provider.altitude_m(), 0.0);
    assert_eq!(provider.pressure_mbar(), 1013.25);
    assert_eq!(provider.temperature_c(), 15.0);
    assert_eq!(provider.humidity_entry(), 0.5);
    assert_eq!(provider.met_range_km(), 0.0);
    assert_eq!(provider.request(), rising_request());
    // The model is the search altitude-form extinction over the
    // release humidity arm, seasoned at the search start.
    let atmo = heliacal_search_atmosphere(1013.25, 15.0, 0.5, 0.0, 41.9, JD0, 0.0).unwrap();
    let reference = heliacal_search_model(&atmo, &[], 0).unwrap();
    assert_eq!(provider.model.borrow().coeffs(), reference.coeffs());
    assert_eq!(provider.season_jd(), JD0);
    // The serving handle pins the channel the words come from.
    assert_eq!(provider.serving(), &());
}

#[test]
fn constructor_refuses_non_finite_site_words() {
    let unit = ();
    let words: [(&str, [f64; 7]); 7] = [
        ("lon", [f64::NAN, 41.9, 0.0, 1013.25, 15.0, 0.5, 0.0]),
        ("lat", [12.5, f64::INFINITY, 0.0, 1013.25, 15.0, 0.5, 0.0]),
        (
            "alt",
            [12.5, 41.9, f64::NEG_INFINITY, 1013.25, 15.0, 0.5, 0.0],
        ),
        ("pressure", [12.5, 41.9, 0.0, f64::NAN, 15.0, 0.5, 0.0]),
        ("temp", [12.5, 41.9, 0.0, 1013.25, f64::NAN, 0.5, 0.0]),
        ("humidity", [12.5, 41.9, 0.0, 1013.25, 15.0, f64::NAN, 0.0]),
        ("met", [12.5, 41.9, 0.0, 1013.25, 15.0, 0.5, f64::INFINITY]),
    ];
    for (name, site) in words {
        let (dt, nut) = healthy_epoch();
        let res: Result<TestSealed<'_>> = SealedHeliacalEphemeris::new(
            &unit,
            rising_request(),
            site[0],
            site[1],
            site[2],
            site[3],
            site[4],
            site[5],
            site[6],
            &[],
            dt,
            nut,
            const_altaz((0.0, 0.0)),
            const_altaz((0.0, 0.0)),
            const_word(20.0),
            const_word(0.0),
            const_moon((-90.0, 0.0, 180.0)),
            const_word(0.0),
        );
        let err = res.unwrap_err();
        assert!(
            matches!(err, Error::NonFiniteArgument { .. }),
            "{name}: {err:?}"
        );
        assert_eq!(
            err.category(),
            crate::errors::ErrorCategory::InputValidation
        );
    }
}

#[test]
fn constructor_keeps_humidity_and_optics_arms() {
    let unit = ();
    // Fraction and percent spellings build bitwise-identical models.
    let (dt, nut) = healthy_epoch();
    let frac: TestSealed<'_> = SealedHeliacalEphemeris::new(
        &unit,
        rising_request(),
        12.5,
        41.9,
        0.0,
        1013.25,
        15.0,
        0.5,
        0.0,
        &[],
        dt,
        nut,
        const_altaz((0.0, 0.0)),
        const_altaz((0.0, 0.0)),
        const_word(20.0),
        const_word(0.0),
        const_moon((-90.0, 0.0, 180.0)),
        const_word(0.0),
    )
    .unwrap();
    let (dt2, nut2) = healthy_epoch();
    let pct: TestSealed<'_> = SealedHeliacalEphemeris::new(
        &unit,
        rising_request(),
        12.5,
        41.9,
        0.0,
        1013.25,
        15.0,
        50.0,
        0.0,
        &[],
        dt2,
        nut2,
        const_altaz((0.0, 0.0)),
        const_altaz((0.0, 0.0)),
        const_word(20.0),
        const_word(0.0),
        const_moon((-90.0, 0.0, 180.0)),
        const_word(0.0),
    )
    .unwrap();
    assert_eq!(frac.model.borrow().coeffs(), pct.model.borrow().coeffs());
    // Instrument words without the optical-params flag stay naked-eye.
    assert!(!frac.model.borrow().observer().binocular());
    let (dt3, nut3) = healthy_epoch();
    let aided: TestSealed<'_> = SealedHeliacalEphemeris::new(
        &unit,
        HeliacalSearchRequest::new(
            JD0,
            VENUS,
            HeliacalEvent::Rising,
            crate::constants::HELFLAG_OPTICAL_PARAMS,
        )
        .unwrap(),
        12.5,
        41.9,
        0.0,
        1013.25,
        15.0,
        0.5,
        0.0,
        &[25.0, 1.2, 1.0, 10.0, 50.0, 0.9],
        dt3,
        nut3,
        const_altaz((0.0, 0.0)),
        const_altaz((0.0, 0.0)),
        const_word(20.0),
        const_word(0.0),
        const_moon((-90.0, 0.0, 180.0)),
        const_word(0.0),
    )
    .unwrap();
    assert!(aided.model.borrow().observer().binocular());
}

#[test]
fn day_scan_drives_search_to_known_crossing() {
    let unit = ();
    let log = Rc::new(RefCell::new(Vec::new()));
    // Step inside the day-5 grid span so the +/-0.1-day refinement
    // window covers it (grid samples fall at fractions ~0.29-0.54,
    // i.e. JD0+5.79.. JD0+6.04).
    let t0 = JD0 + 5.85;
    let (dt, nut) = healthy_epoch();
    let (sun, body) = rising_script(t0);
    let provider = build_sealed(
        &unit,
        rising_request(),
        log,
        sun,
        body,
        move |jd| {
            if jd < JD0 {
                40.0
            } else if jd < JD0 + 5.0 {
                5.0
            } else {
                12.0
            }
        },
        |_| -4.0,
        |_| (-90.0, 0.0, 180.0),
        |_| -15.0,
        dt,
        nut,
    );
    let found = provider.search().unwrap().unwrap();
    assert!((found - t0).abs() < 2e-6, "found {found} vs step {t0}");
}

#[test]
fn sun_acquired_before_body_every_instant() {
    let unit = ();
    let log = Rc::new(RefCell::new(Vec::new()));
    let t0 = JD0 + 5.25;
    let (dt, nut) = healthy_epoch();
    let (sun, body) = rising_script(t0);
    let provider = build_sealed(
        &unit,
        rising_request(),
        log.clone(),
        sun,
        body,
        |_| 12.0,
        |_| -4.0,
        |_| (-90.0, 0.0, 180.0),
        |_| -15.0,
        dt,
        nut,
    );
    // Full evaluation: epoch words, then Sun, then the body (the Sun
    // sits at fraction 0.3 here, deep in twilight).
    assert!(provider.refinement_visible(JD0 + 5.8));
    let calls = log.borrow();
    assert_eq!(&calls[0..4], &["dt", "nut", "sun", "body"]);
    drop(calls);
    log.borrow_mut().clear();
    // Twilight detector: the 25 hourly Sun-only reads precede any
    // body acquisition.
    let (visible, _) = provider.twilight_sample(JD0 + 5.0, true);
    assert!(visible);
    let calls = log.borrow();
    assert!(calls.len() > 30, "grid reads logged: {calls:?}");
    assert!(calls[2..27].iter().all(|c| *c == "sun"), "{calls:?}");
    assert_eq!(calls[29], "body", "{calls:?}");
}

#[test]
fn no_event_when_never_visible() {
    let unit = ();
    let log = Rc::new(RefCell::new(Vec::new()));
    let (dt, nut) = healthy_epoch();
    let (sun, _) = rising_script(JD0 + 5.25);
    let provider = build_sealed(
        &unit,
        rising_request(),
        log,
        sun,
        |_| (-5.0, 60.0),
        |_| 40.0,
        |_| -4.0,
        |_| (-90.0, 0.0, 180.0),
        |_| -15.0,
        dt,
        nut,
    );
    assert_eq!(provider.search().unwrap(), None);
}

#[test]
fn unserved_epoch_degrades_fail_closed() {
    let unit = ();
    // A refused Delta-T word misses every boundary: invisible
    // twilight, sentinel elongation, invisible refinement, no event.
    let log = Rc::new(RefCell::new(Vec::new()));
    let (sun, body) = rising_script(JD0 + 5.25);
    let provider = build_sealed(
        &unit,
        rising_request(),
        log,
        sun,
        body,
        |_| 12.0,
        |_| -4.0,
        |_| (-90.0, 0.0, 180.0),
        |_| -15.0,
        Box::new(|_| {
            Err(Error::DeltaTModelUnavailable {
                model: "sealed-stub",
            })
        }),
        Box::new(|_| Ok((0.0, 0.0))),
    );
    assert_eq!(provider.twilight_sample(JD0 + 5.0, true), (false, 0.0));
    assert_eq!(provider.elongation(JD0 + 5.0), ELONG_WORD_MISS_DEG);
    assert!(!provider.refinement_visible(JD0 + 5.3));
    assert_eq!(provider.search().unwrap(), None);
    // A non-finite Delta-T word degrades identically.
    let log2 = Rc::new(RefCell::new(Vec::new()));
    let (sun2, body2) = rising_script(JD0 + 5.25);
    let provider2 = build_sealed(
        &unit,
        rising_request(),
        log2,
        sun2,
        body2,
        |_| 12.0,
        |_| -4.0,
        |_| (-90.0, 0.0, 180.0),
        |_| -15.0,
        Box::new(|_| Ok(f64::NAN)),
        Box::new(|_| Ok((0.0, 0.0))),
    );
    assert_eq!(provider2.search().unwrap(), None);
}

#[test]
fn wrong_apparition_side_rejected() {
    let unit = ();
    // All-visible words on the evening side never count at morning.
    let log = Rc::new(RefCell::new(Vec::new()));
    let (dt, nut) = healthy_epoch();
    let (sun, body) = rising_script(JD0 - 1.0);
    let provider = build_sealed(
        &unit,
        rising_request(),
        log,
        sun,
        body,
        |_| 12.0,
        |_| -4.0,
        |_| (-90.0, 0.0, 180.0),
        |_| 15.0,
        dt,
        nut,
    );
    assert_eq!(provider.twilight_sample(JD0 + 5.0, true), (false, 0.0));
    // The same words on the matching evening side do count (evening
    // needs a falling Sun: the rising script never yields an evening
    // twilight centre).
    let log2 = Rc::new(RefCell::new(Vec::new()));
    let (dt2, nut2) = healthy_epoch();
    let provider2 = build_sealed(
        &unit,
        rising_request(),
        log2,
        falling_sun(),
        |_| (10.0, 60.0),
        |_| 12.0,
        |_| -4.0,
        |_| (-90.0, 0.0, 180.0),
        |_| 15.0,
        dt2,
        nut2,
    );
    assert!(provider2.twilight_sample(JD0 + 5.0, false).0);
}

#[test]
fn moonlight_counts_in_refinement_only() {
    let unit = ();
    // A risen full Moon beside the object: the moon-excluded twilight
    // detector still fires, but the full refinement refuses. Words are
    // calibrated to ~0.5-mag margins on each side (body 10 deg,
    // elongation word 60 deg, magnitude 1.5, Moon at 30 deg / full /
    // 1 deg separation).
    let log = Rc::new(RefCell::new(Vec::new()));
    let (dt, nut) = healthy_epoch();
    let (sun, body) = rising_script(JD0 - 1.0);
    let provider = build_sealed(
        &unit,
        rising_request(),
        log,
        sun,
        body,
        |_| 60.0,
        |_| 1.5,
        |_| (30.0, 1.0, 1.0),
        |_| -15.0,
        dt,
        nut,
    );
    assert!(provider.twilight_sample(JD0 + 5.0, true).0);
    assert!(!provider.refinement_visible(JD0 + 5.8));
    // Control: Moon set, same geometry refines visible.
    let log2 = Rc::new(RefCell::new(Vec::new()));
    let (dt2, nut2) = healthy_epoch();
    let (sun2, body2) = rising_script(JD0 - 1.0);
    let control = build_sealed(
        &unit,
        rising_request(),
        log2,
        sun2,
        body2,
        |_| 60.0,
        |_| 1.5,
        |_| (-90.0, 0.0, 180.0),
        |_| -15.0,
        dt2,
        nut2,
    );
    assert!(control.refinement_visible(JD0 + 5.8));
}

#[test]
fn elongation_word_passes_through_and_tt_epoch_agrees() {
    let unit = ();
    let log = Rc::new(RefCell::new(Vec::new()));
    let (dt, nut) = healthy_epoch();
    let (sun, body) = rising_script(JD0 + 5.25);
    let provider = build_sealed(
        &unit,
        rising_request(),
        log,
        sun,
        body,
        |_| 17.5,
        |_| -4.0,
        |_| (-90.0, 0.0, 180.0),
        |_| -15.0,
        dt,
        nut,
    );
    assert_eq!(provider.elongation(JD0 + 5.0).to_bits(), 17.5_f64.to_bits());
    // A non-finite elongation word degrades to the miss sentinel.
    let log2 = Rc::new(RefCell::new(Vec::new()));
    let (dt2, nut2) = healthy_epoch();
    let (sun2, body2) = rising_script(JD0 + 5.25);
    let broken = build_sealed(
        &unit,
        rising_request(),
        log2,
        sun2,
        body2,
        |_| f64::NAN,
        |_| -4.0,
        |_| (-90.0, 0.0, 180.0),
        |_| -15.0,
        dt2,
        nut2,
    );
    assert_eq!(broken.elongation(JD0 + 5.0), ELONG_WORD_MISS_DEG);
    // Zero Delta-T reproduces the UT epoch bitwise for the primary's
    // epoch-agreement check.
    assert_eq!(provider.tt_epoch(JD0).unwrap().to_bits(), JD0.to_bits());
    assert!(provider.tt_epoch(f64::NAN).is_none());
}
