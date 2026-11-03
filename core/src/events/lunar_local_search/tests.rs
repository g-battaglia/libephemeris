use super::*;
use std::cell::RefCell;

struct Probe {
    global: RefCell<Vec<f64>>,
    how: RefCell<Vec<f64>>,
    visible: bool,
    global_error_after: usize,
    phase_error: bool,
    horizon_error: bool,
    crossings: BodyHorizonCrossings,
}
fn marker() -> Error {
    Error::EphemerisRange {
        requested_jd: 100.0,
        start_jd: 101.0,
        end_jd: 102.0,
        body: Some(MOON),
    }
    .with_cause(Error::UnknownBody { body_id: 777 })
}
fn times() -> [f64; 10] {
    [100.0, 0.0, 99.8, 100.2, 99.9, 100.1, 99.7, 100.3, 0.0, 0.0]
}
fn probe() -> Probe {
    Probe {
        global: Default::default(),
        how: Default::default(),
        visible: true,
        global_error_after: usize::MAX,
        phase_error: false,
        horizon_error: false,
        crossings: BodyHorizonCrossings {
            rise_rc: -2,
            rise_jd: 0.0,
            set_rc: -2,
            set_jd: 0.0,
        },
    }
}
impl LunarLocalEphemeris for Probe {
    fn next_global(&self, jd: f64, _: bool) -> Result<(i32, [f64; 10])> {
        self.global.borrow_mut().push(jd);
        if self.global.borrow().len() > self.global_error_after {
            return Err(marker());
        }
        Ok((ECL_PARTIAL, times()))
    }
    fn how_attr(&self, jd: f64) -> Result<(i32, [f64; 20])> {
        self.how.borrow_mut().push(jd);
        if self.phase_error {
            return Err(marker());
        }
        let mut attr = [0.0; 20];
        attr[6] = if self.visible { 10.0 } else { -1.0 };
        Ok((ECL_PARTIAL, attr))
    }
    fn moon_horizon_crossings(&self, _: f64) -> Result<BodyHorizonCrossings> {
        if self.horizon_error {
            return Err(marker());
        }
        Ok(self.crossings)
    }
}
#[test]
fn lunar_local_each_candidate_has_epoch_window() {
    for backwards in [false, true] {
        let mut p = probe();
        p.visible = false;
        p.global_error_after = 1;
        assert_eq!(
            lun_eclipse_when_loc_source_via(&p, 95.0, backwards, &marker).unwrap_err(),
            marker()
        );
        let offset = if backwards { 0.35 } else { -0.35 };
        let advance = if backwards { -25.0 } else { 25.0 };
        assert_eq!(
            *p.global.borrow(),
            [95.0 + offset, 100.0 + advance + offset]
        );
    }
}
#[test]
fn lunar_local_phase_and_horizon_range_carriers_propagate() {
    for phase in [false, true] {
        let mut p = probe();
        p.phase_error = phase;
        p.horizon_error = !phase;
        assert_eq!(
            lun_eclipse_when_loc_source_via(&p, 95.0, false, &marker).unwrap_err(),
            marker()
        );
    }
}
#[test]
fn lunar_local_visibility_order_and_circumpolar_contacts() {
    let p = probe();
    let (_, actual, _) = lun_eclipse_when_loc_source_via(&p, 95.0, false, &marker).unwrap();
    assert_eq!(actual, times());
    assert_eq!(
        *p.how.borrow(),
        [100.3, 99.7, 100.1, 99.9, 100.2, 99.8, 100.0, 100.0]
    );
}
#[test]
fn lunar_local_rise_replaces_maximum_after_phase_clamping() {
    let mut p = probe();
    p.crossings = BodyHorizonCrossings {
        rise_rc: 0,
        rise_jd: 100.05,
        set_rc: 0,
        set_jd: 101.0,
    };
    let (_, actual, _) = lun_eclipse_when_loc_source_via(&p, 95.0, false, &marker).unwrap();
    assert_eq!(actual[0], 100.05);
    assert_eq!(actual[6], 0.0);
    assert_eq!(actual[2], 0.0);
    assert_eq!(actual[4], 0.0);
    assert_eq!(actual[8], 100.05);
    assert_eq!(p.how.borrow().last(), Some(&100.05));
}
#[test]
fn lunar_local_only_full_visibility_exhaustion_constructs_own_error() {
    let mut p = probe();
    p.visible = false;
    let error = Error::NetworkSealed {
        purpose: "own visibility exhaustion marker".into(),
    };
    assert_eq!(
        lun_eclipse_when_loc_source_via(&p, 95.0, false, &|| error.clone()).unwrap_err(),
        error
    );
    assert_eq!(p.global.borrow().len(), 250);
}
