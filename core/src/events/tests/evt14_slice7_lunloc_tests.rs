use super::*;
use crate::positions::StateProvider;
use std::cell::Cell;

/// Analytic local-lunar stub: one total eclipse at `jd_max` with the
/// global `lun_eclipse_when` phase layout (maximum, umbral and
/// penumbral contacts) and a scripted apparent altitude. The
/// circumpolar Moon (`rc == -2`) leaves the geometric phase times
/// untouched, so the walk assembly — visibility bits, tret layout
/// and the final `how` composition — is directly observable.
struct StubLunLoc {
    jd_max: f64,
    moon_up: bool,
    /// Optional phase instant where the circumstance core reports a
    /// coverage miss (coverage regression probe; `None` serves all).
    miss_jd: Option<f64>,
}

impl StubLunLoc {
    fn boundary(&self, jd: f64) -> Error {
        Error::EphemerisRange {
            requested_jd: jd,
            start_jd: self.jd_max - 100.0,
            end_jd: self.jd_max + 100.0,
            body: None,
        }
    }

    fn global_tret(&self) -> [f64; 10] {
        let mut tret = [0.0; 10];
        tret[0] = self.jd_max;
        tret[2] = self.jd_max - 0.02;
        tret[3] = self.jd_max + 0.02;
        tret[4] = self.jd_max - 0.01;
        tret[5] = self.jd_max + 0.01;
        tret[6] = self.jd_max - 0.05;
        tret[7] = self.jd_max + 0.05;
        tret
    }
}

impl LunarLocalEphemeris for StubLunLoc {
    fn next_global(&self, jd: f64, backwards: bool) -> Result<(i32, [f64; 10])> {
        use crate::constants::ECL_TOTAL;
        let hit = if backwards {
            jd >= self.jd_max - WHEN_LOC_EPOCH_WINDOW_DAYS
        } else {
            jd <= self.jd_max + WHEN_LOC_EPOCH_WINDOW_DAYS
        };
        if hit {
            Ok((ECL_TOTAL, self.global_tret()))
        } else {
            Err(self.boundary(jd))
        }
    }

    fn how_attr(&self, jd: f64) -> Result<(i32, [f64; 20])> {
        use crate::constants::ECL_TOTAL;
        if self.miss_jd == Some(jd) {
            return Err(self.boundary(jd));
        }
        let mut attr = [0.0; 20];
        attr[6] = if self.moon_up { 30.0 } else { -10.0 };
        Ok((ECL_TOTAL, attr))
    }

    fn moon_horizon_crossings(&self, _jd_from: f64) -> Result<BodyHorizonCrossings> {
        Ok(BodyHorizonCrossings {
            rise_rc: -2,
            rise_jd: 0.0,
            set_rc: -2,
            set_jd: 0.0,
        })
    }
}

#[test]
fn lunloc_walk_finds_total_and_advances() {
    use crate::constants::{
        ECL_MAX_VISIBLE, ECL_PARTBEG_VISIBLE, ECL_PARTEND_VISIBLE, ECL_PENUMBBEG_VISIBLE,
        ECL_PENUMBEND_VISIBLE, ECL_TOTAL, ECL_TOTBEG_VISIBLE, ECL_TOTEND_VISIBLE, ECL_VISIBLE,
    };
    let eph = StubLunLoc {
        jd_max: 2460600.5,
        moon_up: true,
        miss_jd: None,
    };
    // Forward search finds the scripted total eclipse with the full
    // local layout (reported maximum plus the umbral/penumbral
    // contacts; moonrise/moonset stay 0.0 for the circumpolar Moon).
    let (retflag, tret, _) = lun_eclipse_when_loc_via(&eph, 2460600.5 - 10.0, false).unwrap();
    assert_eq!(retflag & ECL_TOTAL, ECL_TOTAL);
    assert_eq!(retflag & ECL_VISIBLE, ECL_VISIBLE);
    for bit in [
        ECL_MAX_VISIBLE,
        ECL_PARTBEG_VISIBLE,
        ECL_PARTEND_VISIBLE,
        ECL_TOTBEG_VISIBLE,
        ECL_TOTEND_VISIBLE,
        ECL_PENUMBBEG_VISIBLE,
        ECL_PENUMBEND_VISIBLE,
    ] {
        assert_eq!(retflag & bit, bit, "missing bit {bit}: {retflag:#x}");
    }
    assert_eq!(tret[0], eph.jd_max);
    assert_eq!((tret[6], tret[7]), (eph.jd_max - 0.05, eph.jd_max + 0.05));
    assert_eq!((tret[8], tret[9]), (0.0, 0.0));
    // Restarting from the reported maximum advances past the only
    // scripted eclipse: exhausted, never the same event twice.
    assert!(matches!(
        lun_eclipse_when_loc_via(&eph, tret[0], false).unwrap_err(),
        Error::EclipseSearchExhausted { .. }
    ));
    // Backward search from after the eclipse finds it.
    let (_, tret_back, _) = lun_eclipse_when_loc_via(&eph, eph.jd_max + 10.0, true).unwrap();
    assert_eq!(tret_back[0], eph.jd_max);
}

#[test]
fn lunloc_walk_skips_moon_down() {
    // A Moon that never clears the horizon leaves no visible phase:
    // exhausted, never a fabricated event.
    let eph = StubLunLoc {
        jd_max: 2460600.5,
        moon_up: false,
        miss_jd: None,
    };
    assert!(matches!(
        lun_eclipse_when_loc_via(&eph, 2460600.5 - 10.0, false).unwrap_err(),
        Error::EclipseSearchExhausted { .. }
    ));
}

#[test]
fn lunloc_phase_boundary_miss_propagates_fix_c21() {
    // FIX-C21: a coverage miss at one phase instant (the global
    // maximum itself is covered) propagates out of the walk,
    // fail-closed — the release calls the per-phase core with no
    // guard. It must NOT advance to another lunation, nor collapse
    // into walk exhaustion.
    let jd_max = 2460600.5;
    let eph = StubLunLoc {
        jd_max,
        moon_up: true,
        miss_jd: Some(jd_max + 0.05),
    };
    assert!(matches!(
        lun_eclipse_when_loc_via(&eph, jd_max - 10.0, false).unwrap_err(),
        Error::EphemerisRange { .. }
    ));
}

/// Call-counting provider that never serves anything, for the
/// binding-level gate tests (same shape as the slice-6 serving
/// tests, plus a read counter proving the mask refusal happens
/// before any provider call).
struct Unserved {
    reads: Cell<usize>,
}

impl StateProvider for Unserved {
    fn provider_name(&self) -> &'static str {
        "evt14-7-dummy/unserved"
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
        self.reads.set(self.reads.get() + 1);
        Err(Error::LebBodyNotFound { body_id })
    }
    fn eval_raw_two_part(
        &self,
        body_id: i32,
        _epoch: crate::types::TwoPartEpoch,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        self.reads.set(self.reads.get() + 1);
        Err(Error::LebBodyNotFound { body_id })
    }
    fn serving_entry_snapshot(&self, body_id: i32) -> Result<crate::positions::EntrySnapshot> {
        self.reads.set(self.reads.get() + 1);
        Err(Error::LebBodyNotFound { body_id })
    }
}

#[test]
fn lunloc_sealed_binding_gates_and_fails_closed() {
    use crate::constants::{
        FLG_BARYCTR, FLG_HELCTR, FLG_J2000, FLG_NONUT, FLG_RADIANS, FLG_SIDEREAL, FLG_SWIEPH,
        FLG_XYZ,
    };
    let provider = Unserved {
        reads: Cell::new(0),
    };
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
            SealedLunarLocalServing::try_new(&serving, &dt, flags, -99.1, 19.4, 2240.0),
            Err(Error::UnsupportedFlagCombination { .. })
        ));
    }
    // Non-finite observer components refuse as invalid-observer
    // .
    for (lon, lat, alt) in [
        (f64::NAN, 19.4, 2240.0),
        (-99.1, f64::INFINITY, 2240.0),
        (-99.1, 19.4, f64::NEG_INFINITY),
    ] {
        assert!(matches!(
            SealedLunarLocalServing::try_new(&serving, &dt, FLG_SWIEPH, lon, lat, alt),
            Err(Error::InvalidObserver { .. })
        ));
    }
    let bound = SealedLunarLocalServing::try_new(&serving, &dt, FLG_SWIEPH, -99.1, 19.4, 2240.0)
        .expect("tropical basis serves");
    assert_eq!(bound.calc_flags(), FLG_SWIEPH);
    assert_eq!(bound.observer_lon_lat_alt(), (-99.1, 19.4, 2240.0));
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
        lun_eclipse_when_loc_sealed(&bound, f64::NAN, false).unwrap_err(),
        Error::NonFiniteArgument { .. }
    ));
    // The annular-only mask is refused before any provider call
    // through the reused candidate step.
    assert!(matches!(
        lun_eclipse_when_via(&bound, 2460400.0, crate::constants::ECL_ANNULAR, false).unwrap_err(),
        Error::InvalidEclipseSearch { .. }
    ));
    assert_eq!(provider.reads.get(), 0);
    // The unserved provider fails closed (never a fabricated event):
    // the typed miss propagates — the lunar global walk
    // converts only ephemeris-boundary misses to exhaustion while
    // other provider failures return as-is (pre-existing
    // behavior, untouched by this slice).
    assert!(
        lun_eclipse_when_loc_sealed(&bound, 2460400.0, false).is_err(),
        "unserved provider must fail closed"
    );
}
