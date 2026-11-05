use super::*;
use crate::constants::VENUS;
use std::rc::Rc;

/// Boxed-word instantiation of the binding (single concrete
/// instantiation; no I/O, no ephemeris, no network).
type TestBinding<'a> = PipelineAHeliacalBinding<
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
/// Body-word step inside the day-5 morning twilight window (found by
/// the search to the bisection tolerance once the stubs below place
/// the window; the elongation word steps here too).
const T0: f64 = JD0 + 4.65;

// Clock and nutation fixture callbacks remain explicitly separate.
#[allow(clippy::type_complexity)]
fn healthy_epoch() -> (
    Box<dyn Fn(f64) -> Result<f64>>,
    Box<dyn Fn(f64) -> Result<(f64, f64)>>,
) {
    (Box::new(|_| Ok(0.0)), Box::new(|_| Ok((0.0, 0.0))))
}

/// Analytically advancing GAST stub (24 h/day solar-rate sweep; a stub
/// rate only — it drives the diurnal Sun sweep through the reduction,
/// it states no physical GAST value).
fn stub_gast_value(jd: f64) -> f64 {
    (6.0 + 24.0 * (jd - JD0)).rem_euclid(24.0)
}

/// Fixed pre-reduced Sun horizontal stub `(alt_deg, az_deg)`: the
/// release `_leb_body_altaz` word for RA 6 h / Dec 0 deg over the
/// scripted GAST sweep, reduced eagerly through the shared
/// hour-angle law (the search dynamics are unchanged; only the
/// binding contract moved the reduction to the caller).
fn stub_sun_altaz() -> Box<dyn Fn(f64) -> (f64, f64)> {
    Box::new(move |jd_ut| altaz_from_radec(6.0, 0.0, 0.0, 0.0, stub_gast_value(jd_ut)))
}

/// Body horizontal stub stepping at `T0`: opposite the Sun before
/// (down through every morning twilight), 2 h west of the Sun after
/// (up through every morning twilight). Same scripted crossing as the
/// retired equatorial stub, pre-reduced eagerly.
fn stub_body_altaz_step() -> Box<dyn Fn(f64) -> (f64, f64)> {
    Box::new(move |jd_ut| {
        let (ra_hours, dec_deg) = if jd_ut < T0 { (12.0, 0.0) } else { (22.0, 0.0) };
        altaz_from_radec(ra_hours, dec_deg, 0.0, 0.0, stub_gast_value(jd_ut))
    })
}

/// Pheno-with-fallback elongation stub stepping with the body word, so
/// the lookback streak sits inside the conjunction gap and the
/// refinement sees a clear glare gate.
fn stub_elong_step() -> Box<dyn Fn(f64) -> f64> {
    Box::new(move |jd| if jd < T0 { 5.0 } else { 12.0 })
}

/// Build a binding over scripted words; every provider call is logged
/// so acquisition order is assertable.
#[allow(clippy::too_many_arguments)]
fn build_binding<'a>(
    serving: &'a (),
    log: Rc<RefCell<Vec<&'static str>>>,
    sun_altaz: Box<dyn Fn(f64) -> (f64, f64)>,
    body_altaz: Box<dyn Fn(f64) -> (f64, f64)>,
    elong: Box<dyn Fn(f64) -> f64>,
    mag: Box<dyn Fn(f64) -> f64>,
    moon: Box<dyn Fn(f64) -> (f64, f64, f64)>,
    lon: Box<dyn Fn(f64) -> f64>,
    delta_t: Box<dyn Fn(f64) -> Result<f64>>,
    nutation: Box<dyn Fn(f64) -> Result<(f64, f64)>>,
) -> TestBinding<'a> {
    let wrap = |name: &'static str, log: Rc<RefCell<Vec<&'static str>>>| {
        move || log.borrow_mut().push(name)
    };
    let sun_log = wrap("sun_altaz", log.clone());
    let body_log = wrap("body_altaz", log.clone());
    let elong_log = wrap("elong", log.clone());
    let mag_log = wrap("mag", log.clone());
    let moon_log = wrap("moon", log.clone());
    let lon_log = wrap("lon", log.clone());
    let dt_log = wrap("dt", log.clone());
    let nut_log = wrap("nut", log);
    let sun_word: Box<dyn Fn(f64) -> (f64, f64)> = Box::new(move |jd| {
        sun_log();
        sun_altaz(jd)
    });
    let body_word: Box<dyn Fn(f64) -> (f64, f64)> = Box::new(move |jd| {
        body_log();
        body_altaz(jd)
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
        lon(jd)
    });
    let dt_word: Box<dyn Fn(f64) -> Result<f64>> = Box::new(move |jd| {
        dt_log();
        delta_t(jd)
    });
    let nut_word: Box<dyn Fn(f64) -> Result<(f64, f64)>> = Box::new(move |jd| {
        nut_log();
        nutation(jd)
    });
    let request = HeliacalSearchRequest::new(JD0, VENUS, HeliacalEvent::Rising, 0).unwrap();
    PipelineAHeliacalBinding::new(
        serving,
        request,
        0.0,
        0.0,
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

fn healthy_binding<'a>(serving: &'a (), log: Rc<RefCell<Vec<&'static str>>>) -> TestBinding<'a> {
    let (dt, nut) = healthy_epoch();
    build_binding(
        serving,
        log,
        stub_sun_altaz(),
        stub_body_altaz_step(),
        stub_elong_step(),
        Box::new(|_| -4.0),
        Box::new(|_| (-90.0, 0.0, 180.0)),
        Box::new(|_| -15.0),
        dt,
        nut,
    )
}

#[test]
fn altaz_reduction_known_values() {
    // Equator, object on the meridian: zenith.
    let (alt, _) = altaz_from_radec(6.0, 0.0, 0.0, 0.0, 6.0);
    assert!((alt - 90.0).abs() < 1e-9, "alt {alt}");
    // Six sidereal hours off the meridian: horizon.
    let (alt, _) = altaz_from_radec(0.0, 0.0, 0.0, 0.0, 6.0);
    assert!(alt.abs() < 1e-9, "alt {alt}");
    // Northern observer, equatorial transit south of zenith:
    // altitude is the co-latitude, azimuth due south (180).
    let (alt, az) = altaz_from_radec(6.0, 0.0, 0.0, 41.9, 6.0);
    assert!((alt - (90.0 - 41.9)).abs() < 1e-9, "alt {alt}");
    assert!((az - 180.0).abs() < 1e-9, "az {az}");
    // Circumpolar transit north of zenith reads azimuth zero.
    let (alt, az) = altaz_from_radec(6.0, 60.0, 0.0, 41.9, 6.0);
    assert!((alt - (90.0 - (60.0 - 41.9))).abs() < 1e-9, "alt {alt}");
    assert!(az.abs() < 1e-9, "az {az}");
    // Azimuth always normalizes to [0, 360).
    for ra in [0.0, 6.0, 12.0, 18.0, 23.999] {
        let (_, az) = altaz_from_radec(ra, 20.0, -73.0, 41.0, 3.7);
        assert!((0.0..360.0).contains(&az), "ra {ra} az {az}");
    }
}

#[test]
fn unknown_equatorial_direction_does_not_become_zenith() {
    assert!(altaz_from_radec(f64::NAN, 0.0, 0.0, 0.0, 6.0).0.is_nan());
    assert!(altaz_from_radec(6.0, f64::NAN, 0.0, 0.0, 6.0).0.is_nan());
    let (altitude, _) = altaz_from_radec(6.0, 0.0, 0.0, 0.0, 6.0);
    assert!((altitude - 90.0).abs() < 1e-12);
}

#[test]
fn eph_flags_masking_is_verbatim() {
    // Only the JPL/SWISSEPH/MOSEPH selector bits survive; heliacal and
    // calculation bits (optical params, high precision, speed, topo,
    // truepos, equatorial) are stripped so they never reach pheno.
    let flags = crate::constants::FLG_SWIEPH
        | crate::constants::HELFLAG_OPTICAL_PARAMS
        | crate::constants::HELFLAG_HIGH_PRECISION
        | crate::constants::FLG_SPEED
        | crate::constants::FLG_TOPOCTR
        | crate::constants::FLG_EQUATORIAL
        | crate::constants::FLG_TRUEPOS;
    assert_eq!(heliacal_eph_flags(flags), crate::constants::FLG_SWIEPH);
    assert_eq!(
        heliacal_eph_flags(crate::constants::FLG_JPLEPH | crate::constants::FLG_MOSEPH),
        crate::constants::FLG_JPLEPH | crate::constants::FLG_MOSEPH
    );
    assert_eq!(heliacal_eph_flags(0), 0);
}

#[test]
fn constructor_validates_entry_words() {
    let unit = ();
    let (dt, nut) = healthy_epoch();
    let request = HeliacalSearchRequest::new(JD0, VENUS, HeliacalEvent::Rising, 0).unwrap();
    // Non-finite site words refuse fail-closed.
    let bad: Result<TestBinding<'_>> = PipelineAHeliacalBinding::new(
        &unit,
        request,
        f64::NAN,
        0.0,
        0.0,
        1013.25,
        15.0,
        0.5,
        0.0,
        &[],
        dt,
        nut,
        stub_sun_altaz(),
        stub_body_altaz_step(),
        stub_elong_step(),
        Box::new(|_| -4.0),
        Box::new(|_| (-90.0, 0.0, 180.0)),
        Box::new(|_| -15.0),
    );
    let err = bad.unwrap_err();
    assert!(matches!(err, Error::NonFiniteArgument { .. }), "{err:?}");
    // A non-finite humidity entry refuses through the reused
    // atmosphere arm (same error the search would surface).
    let (dt2, nut2) = healthy_epoch();
    let bad_hum: Result<TestBinding<'_>> = PipelineAHeliacalBinding::new(
        &unit,
        request,
        0.0,
        0.0,
        0.0,
        1013.25,
        15.0,
        f64::INFINITY,
        0.0,
        &[],
        dt2,
        nut2,
        stub_sun_altaz(),
        stub_body_altaz_step(),
        stub_elong_step(),
        Box::new(|_| -4.0),
        Box::new(|_| (-90.0, 0.0, 180.0)),
        Box::new(|_| -15.0),
    );
    assert!(
        matches!(bad_hum.unwrap_err(), Error::NonFiniteArgument { .. }),
        "humidity entry must refuse"
    );
    // Healthy words echo through the accessors with the serving pin.
    let log = Rc::new(RefCell::new(Vec::new()));
    let binding = healthy_binding(&unit, log);
    assert_eq!(binding.lon_deg(), 0.0);
    assert_eq!(binding.lat_deg(), 0.0);
    assert_eq!(binding.request(), request);
    assert_eq!(binding.serving(), &());
    assert_eq!(binding.tt_epoch(JD0).unwrap().to_bits(), JD0.to_bits());
}

#[test]
fn reduction_serves_pre_reduced_altaz_words() {
    let unit = ();
    let log = Rc::new(RefCell::new(Vec::new()));
    let binding = healthy_binding(&unit, log.clone());
    // One Sun read serves the single pre-reduced word (no equatorial /
    // GAST fan-out: the release `utils.azalt` ECL2HOR reduction lives
    // with the caller, exactly like `_leb_body_altaz`).
    let (alt, az) = binding.sun_altaz_at(JD0 + 4.7);
    assert!(alt.is_finite() && az.is_finite(), "sun ({alt}, {az})");
    assert_eq!(&log.borrow()[..], &["sun_altaz"]);
    log.borrow_mut().clear();
    // One body read follows the same shape on its own word.
    let (alt, az) = binding.body_altaz_at(JD0 + 4.7);
    assert!(alt.is_finite() && az.is_finite(), "body ({alt}, {az})");
    assert_eq!(&log.borrow()[..], &["body_altaz"]);
    // Sun-before-body: back-to-back reads keep the release
    // `_get_altitudes` order observable in the shared log.
    log.borrow_mut().clear();
    let _ = binding.sun_altaz_at(JD0 + 4.7);
    let _ = binding.body_altaz_at(JD0 + 4.7);
    let calls = log.borrow();
    let sun_pos = calls.iter().position(|c| *c == "sun_altaz").unwrap();
    let body_pos = calls.iter().position(|c| *c == "body_altaz").unwrap();
    assert!(sun_pos < body_pos, "{calls:?}");
}

#[test]
fn reduction_keys_altaz_words_by_ut_epoch() {
    use std::cell::Cell;
    let unit = ();
    let seen = Rc::new(Cell::new(f64::NAN));
    let seen_sun = seen.clone();
    // A nonzero Delta-T discriminates the UT key from the TT
    // instant: the pre-reduced closure observes the input UT epoch
    // (it resolves TT internally, exactly like `_leb_body_altaz`
    // threading `jd_ut` into `_topo_ecliptic` and `azalt`).
    let binding = build_binding(
        &unit,
        Rc::new(RefCell::new(Vec::new())),
        Box::new(move |jd| {
            seen_sun.set(jd);
            (10.0, 250.0)
        }),
        stub_body_altaz_step(),
        stub_elong_step(),
        Box::new(|_| -4.0),
        Box::new(|_| (-90.0, 0.0, 180.0)),
        Box::new(|_| -15.0),
        Box::new(|_| Ok(0.001)),
        healthy_epoch().1,
    );
    let jd_ut = JD0 + 4.7;
    let (alt, az) = binding.sun_altaz_at(jd_ut);
    assert!(alt.is_finite() && az.is_finite(), "sun ({alt}, {az})");
    assert_eq!(
        seen.get().to_bits(),
        jd_ut.to_bits(),
        "altaz word must be keyed by the UT epoch, not TT"
    );
}

#[test]
fn reduction_degrades_fail_closed_on_missing_words() {
    let unit = ();
    // A non-finite instant misses both words to (NaN, NaN).
    let log = Rc::new(RefCell::new(Vec::new()));
    let (dt, nut) = healthy_epoch();
    let binding = build_binding(
        &unit,
        log,
        stub_sun_altaz(),
        stub_body_altaz_step(),
        stub_elong_step(),
        Box::new(|_| -4.0),
        Box::new(|_| (-90.0, 0.0, 180.0)),
        Box::new(|_| -15.0),
        dt,
        nut,
    );
    for jd in [f64::NAN, f64::INFINITY] {
        let (alt, az) = binding.sun_altaz_at(jd);
        assert!(alt.is_nan() && az.is_nan(), "jd {jd}: ({alt}, {az})");
        let (alt, az) = binding.body_altaz_at(jd);
        assert!(alt.is_nan() && az.is_nan(), "jd {jd}: ({alt}, {az})");
    }
    // A refused Sun word misses only its own read.
    let log2 = Rc::new(RefCell::new(Vec::new()));
    let (dt2, nut2) = healthy_epoch();
    let broken = build_binding(
        &unit,
        log2,
        Box::new(|_| (f64::NAN, f64::NAN)),
        stub_body_altaz_step(),
        stub_elong_step(),
        Box::new(|_| -4.0),
        Box::new(|_| (-90.0, 0.0, 180.0)),
        Box::new(|_| -15.0),
        dt2,
        nut2,
    );
    let (alt, _) = broken.sun_altaz_at(JD0 + 4.7);
    assert!(alt.is_nan(), "refused Sun word must miss");
    assert!(broken.body_altaz_at(JD0 + 4.7).0.is_finite());
    // A refused body word misses only its own read.
    let log3 = Rc::new(RefCell::new(Vec::new()));
    let (dt3, nut3) = healthy_epoch();
    let body_broken = build_binding(
        &unit,
        log3,
        stub_sun_altaz(),
        Box::new(|_| (f64::NAN, f64::NAN)),
        stub_elong_step(),
        Box::new(|_| -4.0),
        Box::new(|_| (-90.0, 0.0, 180.0)),
        Box::new(|_| -15.0),
        dt3,
        nut3,
    );
    assert!(body_broken.sun_altaz_at(JD0 + 4.7).0.is_finite());
    assert!(body_broken.body_altaz_at(JD0 + 4.7).0.is_nan());
}

#[test]
fn binding_search_finds_scripted_crossing() {
    let unit = ();
    let log = Rc::new(RefCell::new(Vec::new()));
    let binding = healthy_binding(&unit, log);
    let found = binding.search().unwrap().unwrap();
    assert!((found - T0).abs() < 2e-6, "found {found} vs step {T0}");
}

#[test]
fn binding_search_sun_before_body_in_shared_log() {
    let unit = ();
    let log = Rc::new(RefCell::new(Vec::new()));
    let binding = healthy_binding(&unit, log.clone());
    let found = binding.search().unwrap().unwrap();
    assert!((found - T0).abs() < 2e-6, "found {found}");
    let calls = log.borrow();
    let sun_pos = calls.iter().position(|c| *c == "sun_altaz").unwrap();
    let body_pos = calls.iter().position(|c| *c == "body_altaz").unwrap();
    assert!(sun_pos < body_pos, "Sun word precedes body word: {calls:?}");
    // Every full instant reads the Sun word at least as often as the
    // body word (25 Sun-only centre reads precede any body read).
    assert!(
        calls.iter().filter(|c| **c == "sun_altaz").count()
            >= calls.iter().filter(|c| **c == "body_altaz").count(),
        "{calls:?}"
    );
}

#[test]
fn binding_search_missing_words_fail_closed() {
    let unit = ();
    // A refused Delta-T word misses every boundary: no event.
    let log = Rc::new(RefCell::new(Vec::new()));
    let (_, nut) = healthy_epoch();
    let no_epoch = build_binding(
        &unit,
        log,
        stub_sun_altaz(),
        stub_body_altaz_step(),
        stub_elong_step(),
        Box::new(|_| -4.0),
        Box::new(|_| (-90.0, 0.0, 180.0)),
        Box::new(|_| -15.0),
        Box::new(|_| {
            Err(Error::DeltaTModelUnavailable {
                model: "sealed-stub",
            })
        }),
        nut,
    );
    assert_eq!(no_epoch.search().unwrap(), None);
    // A refused Sun horizontal word blinds the twilight raster: the
    // search degrades to the release miss, never an error.
    let log2 = Rc::new(RefCell::new(Vec::new()));
    let (dt2, nut2) = healthy_epoch();
    let no_sun = build_binding(
        &unit,
        log2,
        Box::new(|_| (f64::NAN, f64::NAN)),
        stub_body_altaz_step(),
        stub_elong_step(),
        Box::new(|_| -4.0),
        Box::new(|_| (-90.0, 0.0, 180.0)),
        Box::new(|_| -15.0),
        dt2,
        nut2,
    );
    assert_eq!(no_sun.search().unwrap(), None);
}

// -- slice P2: `vis_limit_mag` entry layer --------------------

fn p2_model() -> VislimModel {
    let atmo = VislimAtmosphere::new(1013.25, 15.0, 40.0, 0.0, 0.0, 40.0, 2460000.0).unwrap();
    VislimModel::new(&atmo, HeliacalObserver::naked_eye(), true)
}

fn p2_geom() -> VislimGeometry {
    VislimGeometry {
        obj_alt_deg: 20.0,
        obj_az_deg: 250.0,
        obj_mag: 1.0,
        sun_alt_deg: -12.0,
        sun_az_deg: 260.0,
        moon_alt_deg: -30.0,
        moon_az_deg: 100.0,
        moon_phase: 0.1,
        sun_obj_angle_deg: 25.0,
        moon_obj_angle_deg: 120.0,
    }
}

#[test]
fn parse_vislim_object_forms() {
    use crate::constants::EARTH;
    // Integer forms.
    assert_eq!(parse_vislim_object("3").unwrap(), VislimBody::Planet(3));
    assert_eq!(parse_vislim_object("  4 ").unwrap(), VislimBody::Planet(4));
    assert_eq!(parse_vislim_object("+5").unwrap(), VislimBody::Planet(5));
    assert_eq!(
        parse_vislim_object("14").unwrap(),
        VislimBody::Planet(EARTH)
    );
    assert_eq!(parse_vislim_object("1").unwrap(), VislimBody::Planet(1));
    // Planet names (case/space insensitive like the release).
    assert_eq!(parse_vislim_object("venus").unwrap(), VislimBody::Planet(3));
    assert_eq!(
        parse_vislim_object("  Mars ").unwrap(),
        VislimBody::Planet(4)
    );
    assert_eq!(
        parse_vislim_object("MOON").unwrap(),
        VislimBody::Planet(MOON)
    );
    // A catalog star resolves to its id under the v2 grammar.
    let sirius = crate::stars::resolve_star2("Sirius").unwrap().star_id;
    assert_eq!(
        parse_vislim_object("Sirius").unwrap(),
        VislimBody::Star(sirius)
    );
}

#[test]
fn parse_vislim_object_refusals() {
    // Empty name refuses.
    assert!(parse_vislim_object("").is_err());
    // The Sun is never valid, in either spelling.
    assert!(matches!(
        parse_vislim_object("0"),
        Err(Error::InvalidBody { body_id: SUN, .. })
    ));
    assert!(matches!(
        parse_vislim_object("Sun"),
        Err(Error::InvalidBody { body_id: SUN, .. })
    ));
    // Parsed integers outside `_PLANET_MAP` refuse as unknown.
    assert!(matches!(
        parse_vislim_object("99"),
        Err(Error::UnknownBody { body_id: 99 })
    ));
    assert!(matches!(
        parse_vislim_object("-7"),
        Err(Error::UnknownBody { body_id: -7 })
    ));
    // Unresolvable names refuse as star-not-found.
    assert!(matches!(
        parse_vislim_object("Not A Star"),
        Err(Error::StarNotFound { .. })
    ));
}

#[test]
fn vis_limit_mag_below_horizon_sentinel() {
    let model = p2_model();
    let mut geom = p2_geom();
    geom.obj_alt_deg = -0.5;
    let (flag, dret) = vis_limit_mag_words(&model, 0, VislimBody::Planet(3), &geom);
    assert_eq!(flag, crate::constants::HELFLAG_BELOW_HORIZON);
    assert_eq!(dret[0], VISLIM_BELOW_HORIZON_MAG);
    assert!(dret[1..].iter().all(|w| *w == 0.0));
}

#[test]
fn vis_limit_mag_layout_and_flag_arms() {
    use crate::constants::{HELFLAG_BELOW_HORIZON, HELFLAG_VISLIM_DARK, HELFLAG_VISLIM_NOMOON};
    let model = p2_model();
    let geom = p2_geom();
    // Plain call: horizontal words echo, reserved slots stay zero.
    let (flag, dret) = vis_limit_mag_words(&model, 0, VislimBody::Planet(3), &geom);
    assert!(flag == 0 || flag == 1, "vision word {flag}");
    assert_ne!(flag, HELFLAG_BELOW_HORIZON);
    assert_eq!(dret[1], geom.obj_alt_deg);
    assert_eq!(dret[2], geom.obj_az_deg);
    assert_eq!(dret[3], geom.sun_alt_deg);
    assert_eq!(dret[4], geom.sun_az_deg);
    assert_eq!(dret[5], geom.moon_alt_deg);
    assert_eq!(dret[6], geom.moon_az_deg);
    assert_eq!(dret[7], geom.obj_mag);
    assert_eq!(dret[8], 0.0);
    assert_eq!(dret[9], 0.0);
    assert!(dret[0].is_finite(), "limiting magnitude {}", dret[0]);
    // Dark-sky arm forces the reported Sun altitude to −90°.
    let (_, dark) = vis_limit_mag_words(&model, HELFLAG_VISLIM_DARK, VislimBody::Planet(3), &geom);
    assert_eq!(dark[3], -90.0);
    // No-moon arm forces the reported Moon altitude to −90°.
    let (_, nomoon) =
        vis_limit_mag_words(&model, HELFLAG_VISLIM_NOMOON, VislimBody::Planet(3), &geom);
    assert_eq!(nomoon[5], -90.0);
    // The Moon object excludes its own light exactly like no-moon.
    let (_, moonobj) = vis_limit_mag_words(&model, 0, VislimBody::Planet(MOON), &geom);
    assert_eq!(moonobj[5], -90.0);
    assert_eq!(moonobj[0], nomoon[0]);
}
