use super::*;
use std::cell::RefCell;

struct Probe {
    reads: RefCell<Vec<&'static str>>,
    candidate_error: Error,
}
impl LunarGlobalEphemeris for Probe {
    fn sun_moon_state(&self, jd: f64) -> Result<SunMoonState> {
        self.reads.borrow_mut().push("fixed-state");
        SunMoonState::try_new(
            0.0,
            0.0,
            1.0,
            0.0,
            (240.0 + (jd - 2451545.0) * 12.0).rem_euclid(360.0),
            0.0,
            0.00257,
            12.0,
        )
    }
    fn candidate_moon_longitude_deg(&self, _: f64) -> Result<f64> {
        self.reads.borrow_mut().push("candidate-moon");
        Err(self.candidate_error.clone())
    }
    fn mean_node_lon_deg(&self, _: f64) -> Result<f64> {
        panic!("after failed Moon-only candidate")
    }
    fn geocentric_xyz_au(&self, _: f64) -> Result<MoonSunXyzAu> {
        panic!("after failed Moon-only candidate")
    }
    fn lun_shadow(&self, _: f64) -> Result<(i32, LunarShadowDcore)> {
        panic!("after failed Moon-only candidate")
    }
}
fn boundary() -> Error {
    Error::EphemerisRange {
        requested_jd: 2451540.0,
        start_jd: 2451541.0,
        end_jd: 2451550.0,
        body: Some(MOON),
    }
}
#[test]
fn lunar_forward_preprobe_classification_boundary_propagates() {
    let eph = Probe {
        reads: Default::default(),
        candidate_error: boundary(),
    };
    assert_eq!(
        lun_eclipse_when_source_via(&eph, 2451545.0, 0, false).unwrap_err(),
        boundary()
    );
    let reads = eph.reads.borrow();
    assert_eq!(reads.last(), Some(&"candidate-moon"));
    assert_eq!(
        reads
            .iter()
            .filter(|role| **role == "candidate-moon")
            .count(),
        1
    );
}
#[test]
fn lunar_main_walk_classification_boundary_becomes_exhaustion() {
    let eph = Probe {
        reads: Default::default(),
        candidate_error: boundary(),
    };
    assert!(matches!(
        lun_eclipse_when_source_via(&eph, 2451545.0, 0, true),
        Err(Error::EclipseSearchExhausted {
            entry: "lun_eclipse_when",
            ..
        })
    ));
    assert_eq!(
        eph.reads
            .borrow()
            .iter()
            .filter(|role| **role == "candidate-moon")
            .count(),
        1
    );
}
#[test]
fn lunar_annular_only_mask_precedes_provider_reads() {
    let eph = Probe {
        reads: Default::default(),
        candidate_error: boundary(),
    };
    assert!(matches!(
        lun_eclipse_when_source_via(&eph, 2451545.0, ECL_ANNULAR | ECL_CENTRAL, false),
        Err(Error::InvalidEclipseSearch {
            entry: "lun_eclipse_when",
            ..
        })
    ));
    assert!(eph.reads.borrow().is_empty());
}
#[test]
fn lunar_candidate_nonboundary_error_is_same_carrier() {
    let marker = Error::NetworkSealed {
        purpose: "actual acquisition".into(),
    }
    .with_cause(Error::UnknownBody { body_id: 778 });
    let eph = Probe {
        reads: Default::default(),
        candidate_error: marker.clone(),
    };
    assert_eq!(
        lun_eclipse_when_source_via(&eph, 2451545.0, 0, false).unwrap_err(),
        marker
    );
}

#[test]
fn lunar_annular_mask_precedes_nonfinite_epoch_rejection() {
    let eph = Probe {
        reads: Default::default(),
        candidate_error: boundary(),
    };
    assert!(matches!(
        lun_eclipse_when_source_via(&eph, f64::NAN, ECL_ANNULAR, false),
        Err(Error::InvalidEclipseSearch {
            entry: "lun_eclipse_when",
            ..
        })
    ));
    assert!(eph.reads.borrow().is_empty());
}
