use super::*;
use crate::sidereal::{DynamicAnchor, DynamicAnchorProvider};
use std::sync::{Arc, Mutex};

type Log = Arc<Mutex<Vec<String>>>;

struct NoProvider;
impl crate::positions::StateProvider for NoProvider {
    fn provider_name(&self) -> &'static str {
        "unused-live-anchor-provider"
    }
    fn has_body(&self, _: i32) -> bool {
        false
    }
    fn body_coverage(&self, _: i32) -> Option<(f64, f64)> {
        None
    }
    fn serving_coord_type(&self, _: i32) -> Option<u32> {
        None
    }
    fn eval_raw(&self, _: i32, _: f64) -> Result<(Vec<f64>, Vec<f64>), Error> {
        panic!("cache or constructor refusal reached provider")
    }
    fn eval_raw_two_part(
        &self,
        _: i32,
        _: crate::types::TwoPartEpoch,
    ) -> Result<(Vec<f64>, Vec<f64>), Error> {
        panic!("cache or constructor refusal reached split provider")
    }
    fn serving_entry_snapshot(&self, _: i32) -> Result<crate::positions::EntrySnapshot, Error> {
        panic!("cache or constructor refusal reached provider metadata")
    }
}

struct Rows {
    log: Log,
    cached: bool,
    refuse_lookup: bool,
    refuse_inputs: bool,
    refuse_publish: bool,
}
struct Row {
    log: Log,
    cached: bool,
    refuse_inputs: bool,
    refuse_publish: bool,
}
struct Sample {
    log: Log,
    refuse_inputs: bool,
    refuse_publish: bool,
}
fn marker() -> Error {
    Error::ConfigurationMissing {
        missing_config: "actual_anchor_callback",
        suggestion: "retained callback refusal",
    }
}
fn inputs() -> NativeStarInputs {
    NativeStarInputs {
        ra_hours: 12.345,
        dec_degrees: -12.5,
        ra_mas_per_year: -25.0,
        dec_mas_per_year: 5.0,
        parallax_mas: 0.0,
        radial_km_per_s: 0.0,
    }
}
impl AnchorRowAcquisition for Rows {
    fn acquire_row(&self, key: &str) -> Result<Arc<dyn AcquiredAnchorRow>, Error> {
        self.log.lock().unwrap().push(format!("lookup:{key}"));
        if self.refuse_lookup {
            return Err(marker());
        }
        Ok(Arc::new(Row {
            log: Arc::clone(&self.log),
            cached: self.cached,
            refuse_inputs: self.refuse_inputs,
            refuse_publish: self.refuse_publish,
        }))
    }
}
impl Drop for Row {
    fn drop(&mut self) {
        self.log.lock().unwrap().push("drop-row".into());
    }
}
impl AcquiredAnchorRow for Row {
    fn sample(&self, tt: f64, options: AnchorSampleOptions) -> Result<AnchorSample, Error> {
        self.log.lock().unwrap().push(format!(
            "sample:{tt}:{}:{}:{}:{}",
            options.nonut, options.geometric, options.noaberr, options.nogdefl
        ));
        assert_eq!(options.eps_true.to_bits(), 0.0f64.to_bits());
        if self.cached {
            return Ok(AnchorSample::Cached(123.456));
        }
        Ok(AnchorSample::Pending(Arc::new(Sample {
            log: Arc::clone(&self.log),
            refuse_inputs: self.refuse_inputs,
            refuse_publish: self.refuse_publish,
        })))
    }
    fn galactic_inputs(&self) -> Result<NativeStarInputs, Error> {
        self.log
            .lock()
            .unwrap()
            .push("galactic-four-keywords".into());
        if self.refuse_inputs {
            return Err(marker());
        }
        Ok(inputs())
    }
}
impl AcquiredAnchorSample for Sample {
    fn native_inputs(&self) -> Result<NativeStarInputs, Error> {
        self.log.lock().unwrap().push("six-keywords".into());
        if self.refuse_inputs {
            return Err(marker());
        }
        Ok(inputs())
    }
    fn publish(&self, longitude: f64) -> Result<(), Error> {
        assert!(longitude.is_finite());
        self.log.lock().unwrap().push("publish".into());
        if self.refuse_publish {
            Err(marker())
        } else {
            Ok(())
        }
    }
}
fn rows(log: &Log) -> Rows {
    Rows {
        log: Arc::clone(log),
        cached: false,
        refuse_lookup: false,
        refuse_inputs: false,
        refuse_publish: false,
    }
}

#[test]
fn actual_anchor_cache_hit_skips_constructor_provider_and_publication() {
    let log = Arc::new(Mutex::new(Vec::new()));
    let mut rows = rows(&log);
    rows.cached = true;
    let provider = NoProvider;
    let anchor = SealedAnchorProvider::new(&provider).with_live_rows(&rows);
    let got = anchor
        .anchor_longitude_deg(DynamicAnchor::Spica, 2451545.125, false, false)
        .unwrap();
    assert_eq!(got.to_bits(), 123.456f64.to_bits());
    assert_eq!(
        *log.lock().unwrap(),
        [
            "lookup:SPICA",
            "sample:2451545.125:true:false:false:false",
            "drop-row"
        ]
    );
}

#[test]
fn actual_anchor_reads_one_row_but_three_node_cache_keys() {
    let log = Arc::new(Mutex::new(Vec::new()));
    let rows = rows(&log);
    let provider = NoProvider;
    let ray = |_: ([f64; 3], [f64; 3]), tt, noaberr, nogdefl, role| {
        assert!(noaberr && nogdefl);
        assert_eq!(role, AnchorRayRole::EclipticLongitude);
        log.lock().unwrap().push(format!("ray:{tt}"));
        Ok([1.0, 2.0, 3.0])
    };
    let range = |_, error| error;
    let anchor = SealedAnchorProvider::new(&provider)
        .with_live_rows(&rows)
        .with_acquired_ray(&ray, &range);
    anchor
        .anchor_longitude_deg(DynamicAnchor::Spica, 2451545.125, true, true)
        .unwrap();
    let seen = log.lock().unwrap();
    assert_eq!(seen.iter().filter(|s| s.starts_with("lookup:")).count(), 1);
    assert_eq!(
        seen.iter()
            .filter(|s| s.starts_with("sample:"))
            .cloned()
            .collect::<Vec<_>>(),
        [
            "sample:2451536:true:false:true:true",
            "sample:2451544:true:false:true:true",
            "sample:2451552:true:false:true:true"
        ]
    );
    assert_eq!(
        seen.iter().filter(|s| s.as_str() == "six-keywords").count(),
        3
    );
    assert_eq!(seen.iter().filter(|s| s.as_str() == "publish").count(), 3);
    assert_eq!(seen.last().unwrap(), "drop-row");
}

#[test]
fn actual_anchor_refusal_never_substitutes_a_static_row() {
    for lookup in [true, false] {
        let log = Arc::new(Mutex::new(Vec::new()));
        let mut rows = rows(&log);
        rows.refuse_lookup = lookup;
        rows.refuse_inputs = !lookup;
        let provider = NoProvider;
        let anchor = SealedAnchorProvider::new(&provider).with_live_rows(&rows);
        assert_eq!(
            anchor.anchor_longitude_deg(DynamicAnchor::Revati, 2451545.0, false, false),
            Err(marker())
        );
        let seen = log.lock().unwrap();
        assert_eq!(seen[0], "lookup:REVATI");
        assert_eq!(
            seen.iter().filter(|s| s.as_str() == "six-keywords").count(),
            usize::from(!lookup)
        );
    }
}

#[test]
fn actual_anchor_publication_refusal_stops_later_nodes() {
    let log = Arc::new(Mutex::new(Vec::new()));
    let mut rows = rows(&log);
    rows.refuse_publish = true;
    let provider = NoProvider;
    let ray = |_, _, _, _, _| Ok([1.0, 2.0, 3.0]);
    let range = |_, error| error;
    let anchor = SealedAnchorProvider::new(&provider)
        .with_live_rows(&rows)
        .with_acquired_ray(&ray, &range);
    assert_eq!(
        anchor.anchor_longitude_deg(DynamicAnchor::Mula, 2451545.0, true, true),
        Err(marker())
    );
    assert_eq!(
        log.lock()
            .unwrap()
            .iter()
            .filter(|s| s.starts_with("sample:"))
            .count(),
        1
    );
}

#[test]
fn actual_geometric_pole_keeps_original_false_correction_cache_fields() {
    let log = Arc::new(Mutex::new(Vec::new()));
    let mut rows = rows(&log);
    rows.cached = true;
    let provider = NoProvider;
    let anchor = SealedAnchorProvider::new(&provider).with_live_rows(&rows);
    anchor
        .anchor_longitude_deg(DynamicAnchor::GalNorthPole, 2451545.0, true, true)
        .unwrap();
    let seen = log.lock().unwrap();
    assert!(
        seen.iter()
            .filter(|s| s.starts_with("sample:"))
            .all(|s| s.ends_with(":true:true:false:false"))
    );
}

#[test]
fn actual_galactic_ra_uses_four_keywords_without_ecliptic_cache() {
    let log = Arc::new(Mutex::new(Vec::new()));
    let rows = rows(&log);
    let provider = NoProvider;
    let ray = |_, _, _, _, role| {
        assert_eq!(role, AnchorRayRole::GalacticRightAscension);
        log.lock().unwrap().push("galactic-ray".into());
        Ok([1.0, 2.0, 3.0])
    };
    let range = |_, error| error;
    let anchor = SealedAnchorProvider::new(&provider)
        .with_live_rows(&rows)
        .with_acquired_ray(&ray, &range);
    anchor
        .anchor_longitude_deg(DynamicAnchor::MulaWilhelm, 2451545.0, true, true)
        .unwrap();
    assert_eq!(
        *log.lock().unwrap(),
        [
            "lookup:GAL_CENTER",
            "galactic-four-keywords",
            "galactic-ray",
            "drop-row"
        ]
    );
}
