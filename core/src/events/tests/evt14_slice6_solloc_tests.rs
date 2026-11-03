use super::*;
use crate::positions::StateProvider;
use std::cell::Cell;

/// Analytic local-solar stub: one partial eclipse at `jd_e` with a
/// parabolic center separation, circular outer contacts and a
/// scripted global candidate. The global facet answers the same
/// maximum; the local facet
/// carries the overlap gate, the visibility bits and the
/// circumstance cores with distinct `attr[3]` / noncentral words so
/// the walk assembly is observable.
struct StubSolLoc {
    jd_e: f64,
    sep_min: f64,
    curve: f64,
    sun_up: bool,
}

impl StubSolLoc {
    fn boundary(&self, jd: f64) -> Error {
        Error::EphemerisRange {
            requested_jd: jd,
            start_jd: self.jd_e - 100.0,
            end_jd: self.jd_e + 100.0,
            body: None,
        }
    }

    fn global_tret(&self) -> [f64; 10] {
        let mut tret = [0.0; 10];
        tret[0] = self.jd_e;
        tret
    }

    fn separation_at(&self, jd: f64) -> f64 {
        self.sep_min + self.curve * (jd - self.jd_e).powi(2)
    }
}

impl SolarGlobalEphemeris for StubSolLoc {
    fn sun_moon_state(&self, _jd_ut: f64) -> Result<SunMoonState> {
        unreachable!("global facet unused: next_global_forward is scripted");
    }

    fn mean_node_lon_deg(&self, _jd_ut: f64) -> Result<f64> {
        unreachable!("global facet unused: next_global_forward is scripted");
    }

    fn gamma_er(&self, _jd_ut: f64) -> Result<f64> {
        unreachable!("global facet unused: next_global_forward is scripted");
    }

    fn shadow_limits_er(&self, _jd_ut: f64) -> Result<(f64, f64)> {
        unreachable!("global facet unused: next_global_forward is scripted");
    }

    fn where_gate(&self, _jd_ut: f64) -> Result<(i32, f64)> {
        unreachable!("global facet unused: next_global_forward is scripted");
    }

    fn contact_residuals_km(&self, _jd_ut: f64) -> Result<(f64, f64, f64)> {
        unreachable!("global facet unused: next_global_forward is scripted");
    }

    fn besselian_x_er(&self, _jd_ut: f64) -> Result<f64> {
        unreachable!("global facet unused: next_global_forward is scripted");
    }

    fn ground_core_width_km(&self, _jd_ut: f64) -> Result<f64> {
        unreachable!("global facet unused: next_global_forward is scripted");
    }
}

impl SolarLocalEphemeris for StubSolLoc {
    fn next_global_forward(&self, jd: f64) -> Result<(i32, [f64; 10])> {
        use crate::constants::ECL_PARTIAL;
        if jd <= self.jd_e + ECLIPSE_WHEN_EPOCH_MARGIN_DAYS {
            Ok((ECL_PARTIAL, self.global_tret()))
        } else {
            Err(self.boundary(jd))
        }
    }

    fn separation_deg(&self, jd: f64) -> Result<f64> {
        Ok(self.separation_at(jd))
    }

    fn radii_sep_deg(&self, jd: f64, inner_moon: bool) -> Result<(f64, f64, f64)> {
        let rsun = 0.2666;
        let mut rmoon = 0.2590;
        if inner_moon {
            rmoon *= 0.9991;
        }
        Ok((rsun, rmoon, self.separation_at(jd)))
    }

    fn sun_altaz_deg(&self, _jd: f64) -> Result<(f64, f64, f64)> {
        if self.sun_up {
            Ok((30.0, 180.0, 30.0))
        } else {
            Ok((-10.0, 180.0, -10.0))
        }
    }

    fn topo_diameters_deg(&self, _jd: f64) -> Result<(f64, f64)> {
        Ok((2.0 * 0.2590, 2.0 * 0.2666))
    }

    fn sun_horizon_crossings(&self, _jd_from: f64) -> Result<BodyHorizonCrossings> {
        Ok(BodyHorizonCrossings {
            rise_rc: -2,
            rise_jd: 0.0,
            set_rc: -2,
            set_jd: 0.0,
        })
    }

    fn how_core(&self, _jd: f64) -> Result<(i32, [f64; 20])> {
        use crate::constants::{ECL_PARTIAL, ECL_VISIBLE};
        Ok((ECL_PARTIAL | ECL_VISIBLE, [0.0; 20]))
    }

    fn where_core(&self, _jd: f64) -> Result<(i32, [f64; 10])> {
        use crate::constants::ECL_NONCENTRAL;
        let mut dcore = [0.0; 10];
        dcore[0] = 50.0;
        Ok((ECL_NONCENTRAL, dcore))
    }
}

#[test]
fn solloc_walk_finds_partial_and_advances() {
    use crate::constants::{
        ECL_1ST_VISIBLE, ECL_4TH_VISIBLE, ECL_MAX_VISIBLE, ECL_NONCENTRAL, ECL_PARTIAL, ECL_VISIBLE,
    };
    let eph = StubSolLoc {
        jd_e: 2460409.28,
        sep_min: 0.1,
        curve: 200.0,
        sun_up: true,
    };
    // Forward search finds the scripted partial eclipse with the
    // full local layout (maximum plus bracketing outer contacts,
    // the where-core width in `attr[3]`, the noncentral bit).
    let (retflag, tret, attr) = sol_eclipse_when_loc_via(&eph, 2460409.28 - 10.0, false).unwrap();
    assert_eq!(
        retflag & (ECL_PARTIAL | ECL_VISIBLE | ECL_MAX_VISIBLE | ECL_1ST_VISIBLE | ECL_4TH_VISIBLE),
        ECL_PARTIAL | ECL_VISIBLE | ECL_MAX_VISIBLE | ECL_1ST_VISIBLE | ECL_4TH_VISIBLE,
    );
    assert_eq!(retflag & ECL_NONCENTRAL, ECL_NONCENTRAL);
    assert!((tret[0] - eph.jd_e).abs() < 1e-6, "max = {}", tret[0]);
    assert!(tret[1] < tret[0] && tret[0] < tret[4], "contacts: {tret:?}");
    assert_eq!(attr[3], 50.0);
    // Restarting from the reported maximum advances past the only
    // scripted eclipse: exhausted, never the same event twice.
    assert!(matches!(
        sol_eclipse_when_loc_via(&eph, tret[0], false).unwrap_err(),
        Error::EclipseSearchExhausted { .. }
    ));
    // Backward search from after the eclipse finds it.
    let (_, tret_back, _) = sol_eclipse_when_loc_via(&eph, eph.jd_e + 100.0, true).unwrap();
    assert!((tret_back[0] - eph.jd_e).abs() < 1e-6);
}

#[test]
fn solloc_walk_skips_night_lunation() {
    // A Sun that never clears the horizon leaves no visible phase:
    // exhausted, never a fabricated event.
    let eph = StubSolLoc {
        jd_e: 2460409.28,
        sep_min: 0.1,
        curve: 200.0,
        sun_up: false,
    };
    assert!(matches!(
        sol_eclipse_when_loc_via(&eph, 2460409.28 - 10.0, false).unwrap_err(),
        Error::EclipseSearchExhausted { .. }
    ));
}

/// Call-counting provider that never serves anything, for the
/// binding-level gate tests (same shape as the slice-5a serving
/// tests, plus a read counter proving the mask refusal happens
/// before any provider call).
struct Unserved {
    reads: Cell<usize>,
}

impl StateProvider for Unserved {
    fn provider_name(&self) -> &'static str {
        "evt14-6-dummy/unserved"
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
fn solloc_sealed_binding_gates_and_fails_closed() {
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
            SealedSolarLocalServing::try_new(&serving, &dt, flags, 12.5, 41.9, 0.0),
            Err(Error::UnsupportedFlagCombination { .. })
        ));
    }
    // Non-finite observer components refuse as invalid-observer
    // .
    for (lon, lat, alt) in [
        (f64::NAN, 41.9, 0.0),
        (12.5, f64::INFINITY, 0.0),
        (12.5, 41.9, f64::NEG_INFINITY),
    ] {
        assert!(matches!(
            SealedSolarLocalServing::try_new(&serving, &dt, FLG_SWIEPH, lon, lat, alt),
            Err(Error::InvalidObserver { .. })
        ));
    }
    let bound = SealedSolarLocalServing::try_new(&serving, &dt, FLG_SWIEPH, 12.5, 41.9, 0.0)
        .expect("tropical basis serves");
    assert_eq!(bound.calc_flags(), FLG_SWIEPH);
    assert_eq!(bound.observer_lon_lat_alt(), (12.5, 41.9, 0.0));
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
        sol_eclipse_when_loc_sealed(&bound, f64::NAN, false).unwrap_err(),
        Error::NonFiniteArgument { .. }
    ));
    // The impossible global type mask (central-partial) is refused
    // before any provider call through the reused candidate step
    // .
    assert!(matches!(
        sol_eclipse_when_glob_via(
            &bound,
            2460400.0,
            crate::constants::ECL_CENTRAL | crate::constants::ECL_PARTIAL,
            EclipseSearchDirection::Forward,
        )
        .unwrap_err(),
        Error::InvalidEclipseSearch { .. }
    ));
    assert_eq!(provider.reads.get(), 0);
    // The unserved provider fails closed (never a fabricated event);
    // the local forward walk converts the coverage miss to the
    // documented exhaustion.
    assert!(matches!(
        sol_eclipse_when_loc_sealed(&bound, 2460400.0, false).unwrap_err(),
        Error::EclipseSearchExhausted { .. }
    ));
}
