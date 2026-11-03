use super::*;
use crate::positions::StateProvider;
use std::cell::Cell;

/// Analytic global-lunar stub: circular coplanar Sun/Moon anchored so
/// oppositions fall exactly on `epoch + k * period`; only the listed
/// `maxima` carry an eclipse (node within the prefilter there, far
/// away at the other oppositions, exercising the walk advance). The
/// circular equatorial vectors carry genuine selenocentric-depth
/// minima at the oppositions, and the scripted canonical shadow
/// (TOTAL near `maxima`, absent elsewhere) drives the same phase
/// slots as the walk tests (`P1/P4` at `±0.05 d`, `U1/U4` at
/// `±0.03324 d`, `U2/U3` at `±0.01 d`).
const WSUN: f64 = 0.9856;
const WMOON: f64 = 13.1764;
const EPOCH: f64 = 2460000.0;
const PERIOD: f64 = 360.0 / (WMOON - WSUN);

struct StubLunGlob {
    epoch: f64,
    sun_lon0: f64,
    wsun: f64,
    wmoon: f64,
    maxima: Vec<f64>,
    state_calls: Cell<usize>,
}

impl StubLunGlob {
    fn sun_lon_rad(&self, jd: f64) -> f64 {
        (self.sun_lon0 + self.wsun * (jd - self.epoch)).to_radians()
    }
    fn moon_lon_rad(&self, jd: f64) -> f64 {
        (self.sun_lon0 + 180.0 + self.wmoon * (jd - self.epoch)).to_radians()
    }
    fn sun_lon_deg(&self, jd: f64) -> f64 {
        self.sun_lon_rad(jd).to_degrees().rem_euclid(360.0)
    }
    fn moon_lon_deg(&self, jd: f64) -> f64 {
        (self.sun_lon_deg(jd) + 180.0 + (self.wmoon - self.wsun) * (jd - self.epoch))
            .rem_euclid(360.0)
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

impl LunarGlobalEphemeris for StubLunGlob {
    fn sun_moon_state(&self, jd_ut: f64) -> Result<SunMoonState> {
        self.state_calls.set(self.state_calls.get() + 1);
        SunMoonState::try_new(
            self.sun_lon_deg(jd_ut),
            0.0,
            1.0,
            self.wsun,
            self.moon_lon_deg(jd_ut),
            0.0,
            0.00257,
            self.wmoon,
        )
    }
    fn mean_node_lon_deg(&self, jd_ut: f64) -> Result<f64> {
        Ok(if self.near_max(jd_ut, 16.0) {
            self.moon_lon_deg(jd_ut) - 5.0
        } else {
            self.moon_lon_deg(jd_ut) + 205.0
        })
    }
    fn geocentric_xyz_au(&self, jd_ut: f64) -> Result<MoonSunXyzAu> {
        let (sl, ml) = (self.sun_lon_rad(jd_ut), self.moon_lon_rad(jd_ut));
        Ok((
            (0.00257 * ml.cos(), 0.00257 * ml.sin(), 0.0),
            (sl.cos(), sl.sin(), 0.0),
        ))
    }
    fn lun_shadow(&self, jd_ut: f64) -> Result<(i32, LunarShadowDcore)> {
        use crate::constants::ECL_TOTAL;
        if self.near_max(jd_ut, 1.0) {
            let tmax = self.nearest_max(jd_ut);
            Ok((
                ECL_TOTAL,
                LunarShadowDcore::try_new(
                    1e-3 * (jd_ut - tmax).abs(),
                    4.324e-5,
                    7.676e-5,
                    1.0,
                    1.0,
                )
                .unwrap(),
            ))
        } else {
            Ok((
                0,
                LunarShadowDcore::try_new(1.0, 4.324e-5, 7.676e-5, 1.0, 1.0).unwrap(),
            ))
        }
    }
}

fn stub(maxima: Vec<f64>) -> StubLunGlob {
    StubLunGlob {
        epoch: EPOCH,
        sun_lon0: 100.0,
        wsun: WSUN,
        wmoon: WMOON,
        maxima,
        state_calls: Cell::new(0),
    }
}

#[test]
fn lunglob_walk_finds_filters_margin_and_mirrors() {
    use crate::constants::ECL_PENUMBRAL;
    use crate::constants::ECL_TOTAL;
    // Three consecutive oppositions; the middle one is quiet (node
    // far), so the walk must skip it via the prefilter advance.
    // Oppositions sit on `EPOCH + k * PERIOD` by construction, so
    // the scripted maxima below are true oppositions.
    let u0 = EPOCH + 2.0 * PERIOD;
    let u1 = u0 + PERIOD;
    let u2 = u0 + 2.0 * PERIOD;
    let eph = stub(vec![u0, u2]);
    // Forward search finds the first scripted maximum with the full
    // 10-slot lunar layout (P1/P4 at +/-0.05 d, U1/U4 at +/-0.03324
    // d, U2/U3 at +/-0.01 d, reserved slots zero).
    let (retflag, tret) = lun_eclipse_when_via(&eph, u0 - 100.0, 0, false).unwrap();
    assert_eq!(retflag, ECL_TOTAL);
    assert!((tret[0] - u0).abs() < 1e-6, "tmax = {}", tret[0]);
    assert_eq!(tret[1], 0.0);
    assert!((tret[6] - (u0 - 0.05)).abs() < 1e-5, "P1 = {}", tret[6]);
    assert!((tret[7] - (u0 + 0.05)).abs() < 1e-5, "P4 = {}", tret[7]);
    assert!((tret[2] - (u0 - 0.03324)).abs() < 1e-5, "U1 = {}", tret[2]);
    assert!((tret[3] - (u0 + 0.03324)).abs() < 1e-5, "U4 = {}", tret[3]);
    assert!((tret[4] - (u0 - 0.01)).abs() < 1e-5, "U2 = {}", tret[4]);
    assert!((tret[5] - (u0 + 0.01)).abs() < 1e-5, "U3 = {}", tret[5]);
    assert_eq!(tret[8], 0.0);
    assert_eq!(tret[9], 0.0);
    // Restarting from the returned maximum advances past the quiet
    // opposition to the second scripted maximum (epoch margin plus
    // prefilter skip in one step).
    let (_, tret_next) = lun_eclipse_when_via(&eph, tret[0], 0, false).unwrap();
    assert!((tret_next[0] - u2).abs() < 1e-6, "next = {}", tret_next[0]);
    assert!((tret_next[0] - u1).abs() > 1.0, "skips quiet u1");
    // Backward search mirrors: from past the last maximum back to it,
    // then from the maximum back past the quiet opposition.
    let (_, tret_back) = lun_eclipse_when_via(&eph, u2 + 100.0, 0, true).unwrap();
    assert!((tret_back[0] - u2).abs() < 1e-6);
    let (_, tret_prev) = lun_eclipse_when_via(&eph, u2, 0, true).unwrap();
    assert!((tret_prev[0] - u0).abs() < 1e-6);
    // A non-matching mask walks the whole stub window: exhausted,
    // never a fabricated event.
    assert!(matches!(
        lun_eclipse_when_via(&eph, u0 - 100.0, ECL_PENUMBRAL, false).unwrap_err(),
        Error::EclipseSearchExhausted { .. }
    ));
    // The annular-only mask is refused before any provider call.
    let quiet_probe = stub(vec![u0]);
    assert!(matches!(
        lun_eclipse_when_via(&quiet_probe, u0, 8, false).unwrap_err(),
        Error::InvalidEclipseSearch { .. }
    ));
    assert_eq!(quiet_probe.state_calls.get(), 0);
}

#[test]
fn lunglob_quiet_window_exhausts_without_fabrication() {
    // One scripted maximum two millennia in the past: every
    // opposition the walk can reach fails the node prefilter, so
    // both directions must report exhaustion rather than an event.
    // This is the release no-event composition for a window with no
    // eclipse (the walk raises instead of returning a zeroed tret).
    // The start sits beyond the 24 800-lunation walk bound from the
    // lone maximum, so neither direction can reach it.
    let eph = stub(vec![EPOCH + 2.0 * PERIOD]);
    let far = EPOCH + 800_000.0;
    assert!(matches!(
        lun_eclipse_when_via(&eph, far, 0, false).unwrap_err(),
        Error::EclipseSearchExhausted { .. }
    ));
    assert!(matches!(
        lun_eclipse_when_via(&eph, far, 0, true).unwrap_err(),
        Error::EclipseSearchExhausted { .. }
    ));
}

/// Dummy provider that never serves anything, for the binding-level
/// gate tests (same shape as the slice-5a serving tests).
struct Unserved;
impl StateProvider for Unserved {
    fn provider_name(&self) -> &'static str {
        "evt14-5b-dummy/unserved"
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
fn lunglob_sealed_binding_gates_and_fails_closed() {
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
            SealedLunarGlobalServing::try_new(&serving, &dt, flags),
            Err(Error::UnsupportedFlagCombination { .. })
        ));
    }
    let bound = SealedLunarGlobalServing::try_new(&serving, &dt, FLG_SWIEPH)
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
        lun_eclipse_when_glob_sealed(&bound, f64::NAN, 0, false).unwrap_err(),
        Error::NonFiniteArgument { .. }
    ));
    // The unserved provider fails closed with the typed miss (never
    // a fabricated event, never a mislabeled exhaustion). The
    // channel-less reads escalate to the unknown-body refusal
    // (same triple as the slice-5a serving tests).
    let err = lun_eclipse_when_glob_sealed(&bound, 2460400.0, 0, false).unwrap_err();
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
