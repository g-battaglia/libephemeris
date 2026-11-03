use super::*;
use crate::constants::MARS;

/// Analytic stub provider: boxed closures keep every test's timeline
/// explicit. No I/O, no ephemeris, no network.
struct StubProvider {
    twilight: Box<dyn Fn(f64, bool) -> (bool, f64)>,
    elong: Box<dyn Fn(f64) -> f64>,
    refine: Box<dyn Fn(f64) -> bool>,
}

impl HeliacalEphemeris for StubProvider {
    fn twilight_sample(&self, jd_day: f64, morning: bool) -> (bool, f64) {
        (self.twilight)(jd_day, morning)
    }
    fn elongation(&self, jd: f64) -> f64 {
        (self.elong)(jd)
    }
    fn refinement_visible(&self, jd: f64) -> bool {
        (self.refine)(jd)
    }
}

const JD0: f64 = 2460000.5;

#[test]
fn request_construction_refusals_and_gates() {
    // Non-finite start is refused fail-closed.
    assert!(HeliacalSearchRequest::new(JD0, VENUS, HeliacalEvent::Rising, 0).is_ok());
    assert!(HeliacalSearchRequest::new(f64::NAN, VENUS, HeliacalEvent::Rising, 0).is_err());
    // Invalid event codes reuse InvalidBody on the heliacal_ut entry
    // (dedicated variant stays primary-owned); they never alias a
    // neighbouring valid event.
    for code in [0, 5, 99, -1] {
        let err = HeliacalSearchRequest::from_code(JD0, VENUS, code, 0).unwrap_err();
        assert!(
            matches!(
                err,
                Error::InvalidBody {
                    body_id: 3,
                    operation: "heliacal_ut"
                }
            ),
            "code {code}: {err:?}"
        );
        assert_eq!(
            err.category(),
            crate::errors::ErrorCategory::InputValidation
        );
    }
    assert_eq!(
        HeliacalSearchRequest::from_code(JD0, VENUS, 1, 0)
            .unwrap()
            .event(),
        HeliacalEvent::Rising
    );
    // Slice-1 gates are reused, not duplicated: Sun, Moon rising,
    // outer-planet evening-first and unknown ids refuse here too.
    assert!(HeliacalSearchRequest::new(JD0, SUN, HeliacalEvent::Rising, 0).is_err());
    assert!(HeliacalSearchRequest::new(JD0, MOON, HeliacalEvent::Rising, 0).is_err());
    assert!(HeliacalSearchRequest::new(JD0, MOON, HeliacalEvent::EveningFirst, 0).is_ok());
    assert!(HeliacalSearchRequest::new(JD0, MARS, HeliacalEvent::EveningFirst, 0).is_err());
    assert!(HeliacalSearchRequest::new(JD0, 99, HeliacalEvent::Rising, 0).is_err());
    // Star split is derived, and the retflag contract spells codes/-1.
    let star =
        HeliacalSearchRequest::new(JD0, FIXSTAR_OFFSET + 7, HeliacalEvent::Setting, 0).unwrap();
    assert!(star.is_star());
    assert_eq!(heliacal_search_retflag(HeliacalEvent::Rising, true), 1);
    assert_eq!(heliacal_search_retflag(HeliacalEvent::MorningLast, true), 4);
    assert_eq!(heliacal_search_retflag(HeliacalEvent::Setting, false), -1);
}

#[test]
fn rising_scan_then_refine_wiring() {
    // Invisible days through the conjunction gap (elong 5), first
    // visibility on day 5: scan fires, refinement converges on the
    // stub step. The morning polarity is asserted inside the stub.
    let t0 = JD0 + 5.25;
    let stub = StubProvider {
        twilight: Box::new(|jd_day, morning| {
            assert!(morning, "rising scans the morning twilight");
            if jd_day >= JD0 + 5.0 {
                (true, jd_day + 0.25)
            } else {
                (false, 0.0)
            }
        }),
        elong: Box::new(|jd| if jd < JD0 + 5.0 { 5.0 } else { 12.0 }),
        refine: Box::new(move |jd| jd >= t0),
    };
    let req = HeliacalSearchRequest::new(JD0, VENUS, HeliacalEvent::Rising, 0).unwrap();
    let found = search_heliacal_event(&req, &stub).unwrap().unwrap();
    assert!((found - t0).abs() < 2e-6, "found {found} vs step {t0}");
    assert_eq!(heliacal_search_retflag(req.event(), true), 1);
}

#[test]
fn rising_ignores_mid_apparition_flicker() {
    // Visible on day 0, a 4-day dip at high elongation, visible again:
    // no 5-streak through a gap, so no event.
    let stub = StubProvider {
        twilight: Box::new(|jd_day, _| {
            let d = (jd_day - JD0).round() as i64;
            if d == 0 || d == 5 {
                (true, jd_day + 0.25)
            } else {
                (false, 0.0)
            }
        }),
        elong: Box::new(|_| 40.0),
        refine: Box::new(|_| true),
    };
    let req = HeliacalSearchRequest::new(JD0, MARS, HeliacalEvent::Rising, 0).unwrap();
    assert_eq!(search_heliacal_event(&req, &stub).unwrap(), None);
}

#[test]
fn forward_scan_evaluates_whole_batch_before_streak_decisions() {
    // Release `_HELIACAL_BATCH` order: the full 100-day batch runs
    // through the twilight sampler before the streak loop reads any
    // elongation word, so the sampler season advances to the batch
    // end. A day-by-day interleave would read each day's elongation
    // right after its twilight sample instead. The event sits in the
    // second batch so the first batch must complete first.
    use std::cell::RefCell;
    use std::rc::Rc;
    let log: Rc<RefCell<Vec<(char, i64)>>> = Rc::new(RefCell::new(Vec::new()));
    let event_day = 150_i64;
    let twilight = {
        let log = Rc::clone(&log);
        Box::new(move |jd_day: f64, _morning: bool| {
            let d = (jd_day - JD0).round() as i64;
            log.borrow_mut().push(('t', d));
            if d == event_day {
                (true, jd_day + 0.28)
            } else {
                (false, 0.0)
            }
        })
    };
    let elong = {
        let log = Rc::clone(&log);
        Box::new(move |jd: f64| {
            let d = (jd - JD0).round() as i64;
            log.borrow_mut().push(('e', d));
            5.0
        })
    };
    let stub = StubProvider {
        twilight,
        elong,
        refine: Box::new(|_| true),
    };
    let found = search_heliacal_rising(JD0, &stub, false, false)
        .expect("search runs")
        .expect("event found");
    assert!(
        found > JD0 + event_day as f64 - 1.0,
        "refined event near day {event_day}, got {found}"
    );
    let log = log.borrow();
    let twilight_99 = log
        .iter()
        .position(|entry| entry.0 == 't' && entry.1 == 99)
        .expect("day 99 twilight-sampled");
    let first_fwd_elong = log
        .iter()
        .position(|entry| entry.0 == 'e' && entry.1 >= 0)
        .expect("forward elongation read");
    assert!(
        twilight_99 < first_fwd_elong,
        "batch 0 completes before streak elongations"
    );
}

#[test]
fn setting_returns_last_visible_before_streak() {
    // Evenings 0-1 visible, then the 5-evening confirming streak: the
    // refined last visibility is returned; evening polarity asserted.
    let t0 = JD0 + 1.3;
    let stub = StubProvider {
        twilight: Box::new(|jd_day, morning| {
            assert!(!morning, "setting scans the evening twilight");
            if jd_day < JD0 + 2.0 {
                (true, jd_day + 0.3)
            } else {
                (false, 0.0)
            }
        }),
        elong: Box::new(|_| 30.0),
        refine: Box::new(move |jd| jd <= t0),
    };
    let req = HeliacalSearchRequest::new(JD0, MARS, HeliacalEvent::Setting, 0).unwrap();
    let found = search_heliacal_event(&req, &stub).unwrap().unwrap();
    assert!((found - t0).abs() < 2e-6, "found {found} vs step {t0}");
}

#[test]
fn setting_no_event_when_never_visible() {
    let stub = StubProvider {
        twilight: Box::new(|_, _| (false, 0.0)),
        elong: Box::new(|_| 45.0),
        refine: Box::new(|_| false),
    };
    let req = HeliacalSearchRequest::new(JD0, MARS, HeliacalEvent::Setting, 0).unwrap();
    assert_eq!(search_heliacal_event(&req, &stub).unwrap(), None);
    assert_eq!(heliacal_search_retflag(req.event(), false), -1);
}

#[test]
fn morning_last_confirms_through_gap() {
    // Visible on day 0, then a 3-morning streak dipping to elong 5:
    // the refined last visibility is returned.
    let t0 = JD0 + 0.2;
    let stub = StubProvider {
        twilight: Box::new(|jd_day, morning| {
            assert!(morning);
            if (jd_day - JD0).abs() < 0.5 {
                (true, JD0 + 0.2)
            } else {
                (false, 0.0)
            }
        }),
        elong: Box::new(|jd| if jd > JD0 + 0.5 { 5.0 } else { 20.0 }),
        refine: Box::new(move |jd| jd >= t0),
    };
    let req = HeliacalSearchRequest::new(JD0, VENUS, HeliacalEvent::MorningLast, 0).unwrap();
    let found = search_heliacal_event(&req, &stub).unwrap().unwrap();
    assert!((found - t0).abs() < 2e-6, "found {found} vs step {t0}");
}

#[test]
fn morning_last_trailing_day_is_raw() {
    // Visible on day 0, then invisible mornings that never dip inside
    // the gap: no confirmation, so the trailing day returns RAW
    // (unrefined, bitwise the provider word) — result layout.
    let stub = StubProvider {
        twilight: Box::new(|jd_day, _| {
            if (jd_day - JD0).abs() < 0.5 {
                (true, JD0 + 0.2)
            } else {
                (false, 0.0)
            }
        }),
        elong: Box::new(|_| 40.0),
        refine: Box::new(|_| true),
    };
    let req = HeliacalSearchRequest::new(JD0, VENUS, HeliacalEvent::MorningLast, 0).unwrap();
    assert_eq!(search_heliacal_event(&req, &stub).unwrap(), Some(JD0 + 0.2));
}

#[test]
fn evening_first_fast_path_in_conjunction() {
    // Elongation inside the gap with the 7-evening gate clear: the
    // next evening visibility IS the event (scan starts at day 1).
    let t0 = JD0 + 1.3;
    let stub = StubProvider {
        twilight: Box::new(|jd_day, morning| {
            assert!(!morning);
            if jd_day >= JD0 + 1.0 {
                (true, jd_day + 0.3)
            } else {
                (false, 0.0)
            }
        }),
        elong: Box::new(|_| 5.0),
        refine: Box::new(move |jd| jd <= t0),
    };
    let req = HeliacalSearchRequest::new(JD0, VENUS, HeliacalEvent::EveningFirst, 0).unwrap();
    let found = search_heliacal_event(&req, &stub).unwrap().unwrap();
    assert!((found - t0).abs() < 2e-6, "found {found} vs step {t0}");
}

#[test]
fn contract_gate_retries_when_event_predates_start() {
    // Day 0 classifies an already-past optimum (jd_best before the
    // start): the gate restarts one day later, where the stub timeline
    // is empty, so the outcome is None — not the stale event.
    let s = 2460000.7_f64;
    let stub = StubProvider {
        twilight: Box::new(move |jd_day, _| {
            if jd_day > s - 0.5 && jd_day < s + 0.5 {
                (true, s - 0.4)
            } else {
                (false, 0.0)
            }
        }),
        elong: Box::new(move |jd| if jd < s + 0.5 { 5.0 } else { 40.0 }),
        refine: Box::new(move |jd| jd >= s - 0.5),
    };
    let req = HeliacalSearchRequest::new(s, VENUS, HeliacalEvent::Rising, 0).unwrap();
    assert_eq!(search_heliacal_event(&req, &stub).unwrap(), None);
}

#[test]
fn search_construction_humidity_and_optics_chain() {
    // Release humidity arm: fractions scale to percent, larger words
    // are already percent; the search model keeps the altitude-form
    // extinction (pressure_scales_extinction = false).
    let frac = heliacal_search_atmosphere(1013.25, 15.0, 0.5, 0.0, 41.9, JD0, 0.0).unwrap();
    assert_eq!(frac.humidity_pct(), 50.0);
    let pct = heliacal_search_atmosphere(1013.25, 15.0, 50.0, 0.0, 41.9, JD0, 0.0).unwrap();
    assert_eq!(pct.humidity_pct(), 50.0);
    assert!(heliacal_search_atmosphere(1013.25, 15.0, f64::NAN, 0.0, 41.9, JD0, 0.0).is_err());
    let words = [25.0, 1.2, 1.0, 10.0, 50.0, 0.9];
    let aided =
        heliacal_search_model(&frac, &words, crate::constants::HELFLAG_OPTICAL_PARAMS).unwrap();
    assert!(aided.observer().binocular());
    let reference = VislimModel::new(&frac, aided.observer(), false);
    assert_eq!(aided.coeffs(), reference.coeffs());
}
