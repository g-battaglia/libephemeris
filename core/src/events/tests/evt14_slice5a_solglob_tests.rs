use super::*;
use crate::positions::StateProvider;
use std::cell::Cell;

/// Analytic global-solar stub: circular coplanar Sun/Moon anchored so
/// conjunctions fall exactly on `epoch + k * period`; only the listed
/// `maxima` carry an eclipse (node within the prefilter there, far
/// away at the other conjunctions, exercising the walk advance). The
/// gamma parabola, the linear contact residuals and the Besselian-x
/// line are closed-form scripts with distinct widths from the event-walk
/// model, so the expected slots are exact arithmetic.
/// Circular-model rates (deg/day) and epoch: conjunctions fall
/// exactly on `EPOCH + k * PERIOD`, so scripted maxima placed there
/// are true conjunctions the walk can refine.
const WSUN: f64 = 0.9856;
const WMOON: f64 = 13.1764;
const EPOCH: f64 = 2450000.0;
const PERIOD: f64 = 360.0 / (WMOON - WSUN);

struct StubSolGlob {
    epoch: f64,
    sun_lon0: f64,
    wsun: f64,
    wmoon: f64,
    maxima: Vec<f64>,
    state_calls: Cell<usize>,
}

impl StubSolGlob {
    fn sun_lon(&self, jd: f64) -> f64 {
        self.sun_lon0 + self.wsun * (jd - self.epoch)
    }
    fn moon_lon(&self, jd: f64) -> f64 {
        // Conjunctions exactly at `epoch + k * period` by construction.
        self.sun_lon(jd) + (self.wmoon - self.wsun) * (jd - self.epoch)
    }
    fn nearest_max(&self, jd: f64) -> f64 {
        let mut best = self.maxima[0];
        let mut dist = (jd - best).abs();
        for m in &self.maxima[1..] {
            let d = (jd - *m).abs();
            if d < dist {
                dist = d;
                best = *m;
            }
        }
        best
    }
    fn near_max(&self, jd: f64, window: f64) -> bool {
        self.maxima.iter().any(|m| (jd - m).abs() < window)
    }
}

impl SolarGlobalEphemeris for StubSolGlob {
    fn sun_moon_state(&self, jd_ut: f64) -> Result<SunMoonState> {
        self.state_calls.set(self.state_calls.get() + 1);
        SunMoonState::try_new(
            self.sun_lon(jd_ut),
            0.0,
            1.0,
            self.wsun,
            self.moon_lon(jd_ut),
            0.0,
            0.0025,
            self.wmoon,
        )
    }
    fn mean_node_lon_deg(&self, jd_ut: f64) -> Result<f64> {
        Ok(if self.near_max(jd_ut, 16.0) {
            self.moon_lon(jd_ut) - 5.0
        } else {
            self.moon_lon(jd_ut) + 205.0
        })
    }
    fn gamma_er(&self, jd_ut: f64) -> Result<f64> {
        let d = (jd_ut - self.nearest_max(jd_ut)).abs();
        Ok(0.15 + 150.0 * d * d)
    }
    fn shadow_limits_er(&self, _jd_ut: f64) -> Result<(f64, f64)> {
        Ok((0.53, -0.02))
    }
    fn where_gate(&self, jd_ut: f64) -> Result<(i32, f64)> {
        use crate::constants::ECL_CENTRAL;
        Ok(if self.near_max(jd_ut, 1.0) {
            (ECL_CENTRAL, 0.6)
        } else {
            (0, -1.0)
        })
    }
    fn contact_residuals_km(&self, jd_ut: f64) -> Result<(f64, f64, f64)> {
        let d = (jd_ut - self.nearest_max(jd_ut)).abs() * 1.0e6;
        Ok((d - 0.05e6, d - 0.02e6, d - 0.01e6))
    }
    fn besselian_x_er(&self, jd_ut: f64) -> Result<f64> {
        Ok(jd_ut - self.nearest_max(jd_ut))
    }
    fn ground_core_width_km(&self, _jd_ut: f64) -> Result<f64> {
        Ok(1000.0)
    }
}

fn stub(maxima: Vec<f64>) -> StubSolGlob {
    StubSolGlob {
        epoch: EPOCH,
        sun_lon0: 40.0,
        wsun: WSUN,
        wmoon: WMOON,
        maxima,
        state_calls: Cell::new(0),
    }
}

#[test]
fn solglob_walk_finds_filters_margin_and_mirrors() {
    use crate::constants::{ECL_ANNULAR, ECL_PARTIAL};
    use crate::constants::{ECL_CENTRAL, ECL_TOTAL};
    // Three consecutive conjunctions; the middle one is quiet (node
    // far), so the walk must skip it via the prefilter advance.
    // Conjunctions sit on `EPOCH + k * PERIOD` by construction, so
    // the scripted maxima below are true conjunctions.
    let c0 = EPOCH + 2.0 * PERIOD;
    let c1 = c0 + PERIOD;
    let c2 = c0 + 2.0 * PERIOD;
    let eph = stub(vec![c0, c2]);
    // Forward search finds the first scripted maximum with the full
    // 10-slot global layout (P1/P4 at +/-0.05 d, U1/U4 at +/-0.02 d,
    // center line at +/-0.01 d, noon at maximum, hybrid slots zero).
    let (retflag, tret) =
        sol_eclipse_when_glob_via(&eph, c0 - 100.0, 0, EclipseSearchDirection::Forward).unwrap();
    assert_eq!(retflag, ECL_TOTAL | ECL_CENTRAL);
    assert!((tret[0] - c0).abs() < 1e-6, "tmax = {}", tret[0]);
    assert!((tret[1] - c0).abs() < 1e-6, "noon = {}", tret[1]);
    assert!((tret[2] - (c0 - 0.05)).abs() < 1e-6);
    assert!((tret[3] - (c0 + 0.05)).abs() < 1e-6);
    assert!((tret[4] - (c0 - 0.02)).abs() < 1e-6);
    assert!((tret[5] - (c0 + 0.02)).abs() < 1e-6);
    assert!((tret[6] - (c0 - 0.01)).abs() < 1e-6);
    assert!((tret[7] - (c0 + 0.01)).abs() < 1e-6);
    assert_eq!(tret[8], 0.0);
    assert_eq!(tret[9], 0.0);
    // Restarting from the returned maximum advances past the quiet
    // conjunction to the second scripted maximum (epoch margin plus
    // prefilter skip in one step).
    let (_, tret_next) =
        sol_eclipse_when_glob_via(&eph, tret[0], 0, EclipseSearchDirection::Forward).unwrap();
    assert!((tret_next[0] - c2).abs() < 1e-6, "next = {}", tret_next[0]);
    assert!((tret_next[0] - c1).abs() > 1.0, "skips quiet c1");
    // Backward search mirrors: from past the last maximum back to it,
    // then from the maximum back past the quiet conjunction.
    let (_, tret_back) =
        sol_eclipse_when_glob_via(&eph, c2 + 100.0, 0, EclipseSearchDirection::Backward).unwrap();
    assert!((tret_back[0] - c2).abs() < 1e-6);
    let (_, tret_prev) =
        sol_eclipse_when_glob_via(&eph, c2, 0, EclipseSearchDirection::Backward).unwrap();
    assert!((tret_prev[0] - c0).abs() < 1e-6);
    // Bidirectional search just after a maximum still finds the next.
    let (_, tret_bi) =
        sol_eclipse_when_glob_via(&eph, c0 + 5.0, 0, EclipseSearchDirection::Bidirectional)
            .unwrap();
    assert!((tret_bi[0] - c2).abs() < 1e-6);
    // A non-matching mask walks the whole stub window: exhausted,
    // never a fabricated event.
    assert!(matches!(
        sol_eclipse_when_glob_via(
            &eph,
            c0 - 100.0,
            ECL_ANNULAR,
            EclipseSearchDirection::Forward
        )
        .unwrap_err(),
        Error::EclipseSearchExhausted { .. }
    ));
    // The impossible mask is refused before any provider call.
    let quiet_probe = stub(vec![c0]);
    assert!(matches!(
        sol_eclipse_when_glob_via(
            &quiet_probe,
            c0,
            ECL_CENTRAL | ECL_PARTIAL,
            EclipseSearchDirection::Forward
        )
        .unwrap_err(),
        Error::InvalidEclipseSearch { .. }
    ));
    assert_eq!(quiet_probe.state_calls.get(), 0);
}

#[test]
fn solglob_quiet_window_exhausts_without_fabrication() {
    // One scripted maximum two millennia in the past: every
    // conjunction the walk can reach fails the node prefilter, so
    // both directions must report exhaustion rather than an event.
    // This is the release no-event composition for a window with no
    // eclipse (the walk raises instead of returning a zeroed tret).
    // The start sits beyond the 24 800-lunation walk bound from the
    // lone maximum, so neither direction can reach it.
    let eph = stub(vec![EPOCH + 2.0 * PERIOD]);
    let far = EPOCH + 800_000.0;
    assert!(matches!(
        sol_eclipse_when_glob_via(&eph, far, 0, EclipseSearchDirection::Forward).unwrap_err(),
        Error::EclipseSearchExhausted { .. }
    ));
    assert!(matches!(
        sol_eclipse_when_glob_via(&eph, far, 0, EclipseSearchDirection::Backward).unwrap_err(),
        Error::EclipseSearchExhausted { .. }
    ));
}

/// Dummy provider that never serves anything, for the binding-level
/// gate tests (same shape as the slice-4 serving tests).
struct Unserved;
impl StateProvider for Unserved {
    fn provider_name(&self) -> &'static str {
        "evt14-5a-dummy/unserved"
    }
    fn has_body(&self, _body_id: i32) -> bool {
        false
    }
    fn body_coverage(&self, _body_id: i32) -> Option<(f64, f64)> {
        None
    }
    fn serving_coord_type(&self, _body_id: i32) -> Option<u32> {
        None
    }
    fn eval_raw(&self, body_id: i32, _jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        Err(Error::LebBodyNotFound { body_id })
    }
    fn eval_raw_two_part(
        &self,
        body_id: i32,
        _epoch: crate::types::TwoPartEpoch,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        Err(Error::LebBodyNotFound { body_id })
    }
    fn serving_entry_snapshot(&self, body_id: i32) -> Result<crate::positions::EntrySnapshot> {
        Err(Error::LebBodyNotFound { body_id })
    }
}

#[test]
fn solglob_sealed_binding_gates_and_fails_closed() {
    use crate::constants::{
        FLG_BARYCTR, FLG_HELCTR, FLG_J2000, FLG_NONUT, FLG_RADIANS, FLG_SIDEREAL, FLG_SWIEPH,
        FLG_XYZ,
    };
    let provider = Unserved;
    let nut = |_: f64| Ok((0.0, 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let dt = |_: f64| Ok(0.0008);
    // Selections outside the verified tropical basis refuse at
    // construction.
    for flags in [
        FLG_SWIEPH | FLG_SIDEREAL,
        FLG_SWIEPH | FLG_J2000,
        FLG_SWIEPH | FLG_NONUT,
        FLG_SWIEPH | FLG_XYZ,
        FLG_SWIEPH | FLG_RADIANS,
        FLG_SWIEPH | FLG_HELCTR,
        FLG_SWIEPH | FLG_BARYCTR,
    ] {
        assert!(matches!(
            SealedSolarGlobalServing::try_new(&serving, &dt, flags),
            Err(Error::UnsupportedFlagCombination { .. })
        ));
    }
    let bound = SealedSolarGlobalServing::try_new(&serving, &dt, FLG_SWIEPH)
        .expect("tropical basis serves");
    assert_eq!(bound.calc_flags(), FLG_SWIEPH);
    // Non-finite epochs refuse before any provider read.
    assert!(matches!(
        bound.jd_tt_for_ut(f64::NAN),
        Err(Error::NonFiniteArgument { .. })
    ));
    assert!(matches!(
        bound.sun_moon_state(f64::NAN),
        Err(Error::NonFiniteArgument { .. })
    ));
    assert!(matches!(
        sol_eclipse_when_glob_sealed(&bound, f64::NAN, 0, EclipseSearchDirection::Forward)
            .unwrap_err(),
        Error::NonFiniteArgument { .. }
    ));
    // The unserved provider fails closed with the typed miss (never
    // a fabricated event, never a mislabeled exhaustion). The
    // channel-less Sun read escalates to the unknown-body refusal
    // (same triple as the slice-4c serving tests).
    let err = sol_eclipse_when_glob_sealed(&bound, 2460400.0, 0, EclipseSearchDirection::Forward)
        .unwrap_err();
    assert!(
        matches!(
            err,
            Error::LebBodyNotFound { .. }
                | Error::EphemerisRange { .. }
                | Error::UnknownBody { .. }
        ),
        "unexpected: {err:?}"
    );
}
