//! Synthetic strict-acquisition and public request-clock controls.
//! No real file/provider is loaded; no numerical source fixtures are defined.
use super::*;
use crate::constants::{EARTH, J2000, JUPITER, SATURN, SUN};
use crate::positions::{EntrySnapshot, StateProvider};
use crate::sidereal::{DynamicAnchor, DynamicAnchorProvider};
use crate::types::TwoPartEpoch;
use std::cell::RefCell;

struct UnusedProvider;
impl StateProvider for UnusedProvider {
    fn provider_name(&self) -> &'static str {
        "unused-synthetic"
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
        panic!("native callback must own acquisition")
    }
    fn eval_raw_two_part(&self, _: i32, _: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>), Error> {
        panic!("native callback must own acquisition")
    }
    fn serving_entry_snapshot(&self, _: i32) -> Result<EntrySnapshot, Error> {
        panic!("native callback must own acquisition")
    }
}

fn raw_range(target: i32, jd: f64) -> Error {
    Error::SpkOutOfRange {
        target,
        et: (jd - J2000) * 86400.0,
        start: 100.0,
        end: 200.0,
    }
}

fn public_range(request_tt: f64, _: Error) -> Error {
    Error::EphemerisRange {
        requested_jd: request_tt,
        start_jd: J2000 - 100.0,
        end_jd: J2000 + 100.0,
        body: None,
    }
}

fn assert_public_request(error: Error, request: f64) {
    assert_eq!(error, public_range(request, raw_range(10, request)));
}

#[test]
fn native_anchor_every_observation_deflector_range_is_strict_with_raw_target_identity() {
    let request = J2000 + 0.125;
    for (failing_body, raw_target) in [(SUN, 10), (JUPITER, 599), (SATURN, 699)] {
        for noaberr in [false, true] {
            let reads = RefCell::new(Vec::new());
            let wraps = RefCell::new(Vec::new());
            let state_at = |body, jd| {
                reads.borrow_mut().push((body, jd));
                if body == failing_body {
                    return Err(raw_range(raw_target, jd));
                }
                Ok((
                    if body == EARTH {
                        [0.0; 3]
                    } else {
                        [2.0, 1.0, 0.5]
                    },
                    [0.0; 3],
                ))
            };
            let wrap = |tt, error: Error| {
                wraps.borrow_mut().push((tt, error.clone()));
                public_range(tt, error)
            };
            let anchors = SealedAnchorProvider::new(&UnusedProvider)
                .with_native_acquisition(&state_at, &wrap);
            let error = anchors
                .anchor_longitude_deg(DynamicAnchor::Spica, request, noaberr, false)
                .unwrap_err();
            assert_public_request(error, request);
            assert_eq!(*wraps.borrow(), [(request, raw_range(raw_target, request))]);
            assert_eq!(reads.borrow()[0], (EARTH, request));
            assert_eq!(reads.borrow().last(), Some(&(failing_body, request)));
        }
    }
}

#[test]
fn native_anchor_closest_approach_range_keeps_raw_sample_before_request_wrapping() {
    let request = J2000 + 0.125;
    let data = crate::sidereal::anchor_star("SPICA").unwrap();
    let direction = star_icrs_state_catalog_angles(&data, request).0;
    let sun_reads = RefCell::new(Vec::new());
    let wraps = RefCell::new(Vec::new());
    let state_at = |body, jd| {
        if body == EARTH {
            return Ok(([0.0; 3], [0.0; 3]));
        }
        assert_eq!(body, SUN);
        sun_reads.borrow_mut().push(jd);
        if sun_reads.borrow().len() == 2 {
            return Err(raw_range(10, jd));
        }
        // A synthetic deflector in the forward ray makes closest approach
        // strictly retarded; no golden astronomical value is asserted.
        Ok((direction.map(|value| value * 0.001), [0.0; 3]))
    };
    let wrap = |tt, error: Error| {
        wraps.borrow_mut().push((tt, error.clone()));
        public_range(tt, error)
    };
    let anchors =
        SealedAnchorProvider::new(&UnusedProvider).with_native_acquisition(&state_at, &wrap);
    assert_public_request(
        anchors
            .anchor_longitude_deg(DynamicAnchor::Spica, request, false, false)
            .unwrap_err(),
        request,
    );
    let samples = sun_reads.borrow();
    assert_eq!(samples.len(), 2);
    assert_eq!(samples[0], request);
    assert!(samples[1] < request);
    assert_eq!(*wraps.borrow(), [(request, raw_range(10, samples[1]))]);
}

#[test]
fn native_anchor_all_interpolation_node_failures_keep_original_request() {
    let request = J2000 + 0.125;
    let base = (request / 16.0).floor() * 16.0;
    let nodes = [base, base + 8.0, base + 16.0];
    for failed in 0..3 {
        let reads = RefCell::new(Vec::new());
        let wraps = RefCell::new(Vec::new());
        let state_at = |body, jd| {
            assert_eq!(body, EARTH);
            reads.borrow_mut().push(jd);
            if jd == nodes[failed] {
                return Err(raw_range(399, jd));
            }
            Ok(([0.0; 3], [0.0; 3]))
        };
        let wrap = |tt, error: Error| {
            wraps.borrow_mut().push((tt, error.clone()));
            public_range(tt, error)
        };
        let anchors =
            SealedAnchorProvider::new(&UnusedProvider).with_native_acquisition(&state_at, &wrap);
        assert_public_request(
            anchors
                .anchor_longitude_deg(DynamicAnchor::Spica, request, false, true)
                .unwrap_err(),
            request,
        );
        assert_eq!(*reads.borrow(), nodes[..=failed]);
        assert_eq!(*wraps.borrow(), [(request, raw_range(399, nodes[failed]))]);
    }
}

#[test]
fn native_anchor_geometric_and_nogdefl_paths_do_not_acquire_deflectors() {
    let request = J2000 + 0.125;
    for (noaberr, nogdefl) in [(true, true), (false, true)] {
        let reads = RefCell::new(Vec::new());
        let state_at = |body, jd| {
            assert_eq!(body, EARTH, "unused deflector was consumed");
            reads.borrow_mut().push(jd);
            Ok(([0.0; 3], [0.0; 3]))
        };
        let wrap = |_, _: Error| -> Error { panic!("no raw failure") };
        let anchors =
            SealedAnchorProvider::new(&UnusedProvider).with_native_acquisition(&state_at, &wrap);
        let longitude = anchors
            .anchor_longitude_deg(DynamicAnchor::Spica, request, noaberr, nogdefl)
            .unwrap();
        assert!(longitude.is_finite());
        assert_eq!(reads.borrow().len(), 3);
    }
}

#[test]
fn native_anchor_galactic_projection_wraps_earth_range_at_original_request() {
    let request = J2000 + 0.125;
    for anchor in [DynamicAnchor::GalCenter, DynamicAnchor::MulaWilhelm] {
        let samples = RefCell::new(Vec::new());
        let state_at = |body, jd| {
            assert_eq!(body, EARTH);
            samples.borrow_mut().push(jd);
            Err(raw_range(399, jd))
        };
        let anchors = SealedAnchorProvider::new(&UnusedProvider)
            .with_native_acquisition(&state_at, &public_range);
        assert_public_request(
            anchors
                .anchor_longitude_deg(anchor, request, false, false)
                .unwrap_err(),
            request,
        );
        assert_eq!(*samples.borrow(), [request]);
    }
}

#[test]
fn native_anchor_missing_alias_and_nonfinite_epoch_bypass_range_wrapper() {
    let state_at = |body, _| {
        if body == EARTH {
            Ok(([0.0; 3], [0.0; 3]))
        } else if body == SUN {
            Ok(([2.0, 1.0, 0.5], [0.0; 3]))
        } else {
            assert_eq!(body, JUPITER);
            Err(Error::SpkBodyNotFound { target: 599 })
        }
    };
    let wrap = |_, _: Error| -> Error { panic!("non-range failure was relabelled") };
    let anchors =
        SealedAnchorProvider::new(&UnusedProvider).with_native_acquisition(&state_at, &wrap);
    assert_eq!(
        anchors.anchor_longitude_deg(DynamicAnchor::Spica, J2000, false, false),
        Err(Error::SpkBodyNotFound { target: 599 })
    );
    assert!(matches!(
        anchors.anchor_longitude_deg(DynamicAnchor::Spica, f64::NAN, false, false),
        Err(Error::NonFiniteArgument { .. })
    ));
}

#[test]
fn native_anchor_opt_in_does_not_replace_legacy_optional_deflector_policy() {
    let legacy = SealedAnchorProvider::new(&UnusedProvider);
    assert!(legacy.native_state_at.is_none());
    assert!(legacy.range_error_at.is_none());
    let ray = [2.0, 3.0, 4.0];
    assert_eq!(
        crate::positions::apply_gravitational_deflection(ray, [0.0; 3], J2000, 0.25, |_, jd| Err(
            raw_range(10, jd)
        ))
        .unwrap(),
        ray
    );
}

#[test]
fn native_star_serving_routes_ray_nodes_and_wraps_only_at_original_request() {
    let request = J2000 + 0.125;
    let base = (request / 16.0).floor() * 16.0;
    let nodes = [base, base + 8.0, base + 16.0];
    for failed in 0..3 {
        let seen = RefCell::new(Vec::new());
        let ray = |_: &StarData, tt: f64, noaberr: bool, nogdefl: bool| {
            assert!(!noaberr && nogdefl);
            seen.borrow_mut().push(tt);
            if tt == nodes[failed] {
                return Err(raw_range(399, tt));
            }
            Ok([1.0, 2.0, 3.0])
        };
        let wraps = RefCell::new(Vec::new());
        let wrap = |tt, error: Error| {
            wraps.borrow_mut().push((tt, error.clone()));
            public_range(tt, error)
        };
        let nutation = |_| -> Result<(f64, f64), Error> { panic!("anchor frame owns nutation") };
        let delta = |_| -> Result<f64, Error> {
            panic!("configured Delta-T is not an anchor evaluation clock")
        };
        let serving = SealedStarServing::new(
            &UnusedProvider,
            &nutation,
            &delta,
            crate::sidereal::SiderealSelection::default(),
        )
        .with_default_timescale()
        .with_native_anchors(&ray, &wrap);
        assert_public_request(
            serving
                .anchor_provider()
                .anchor_longitude_deg(DynamicAnchor::Spica, request, false, true)
                .unwrap_err(),
            request,
        );
        assert_eq!(*seen.borrow(), nodes[..=failed]);
        assert_eq!(*wraps.borrow(), [(request, raw_range(399, nodes[failed]))]);
    }
}

#[test]
fn native_star_ray_seam_covers_exact_date_and_galactic_hour_circle() {
    let request = J2000 + 0.125;
    for anchor in [
        DynamicAnchor::Spica,
        DynamicAnchor::GalCenter,
        DynamicAnchor::MulaWilhelm,
    ] {
        let seen = RefCell::new(Vec::new());
        let ray = |_: &StarData, tt, _, _| {
            seen.borrow_mut().push(tt);
            Err(raw_range(599, tt))
        };
        let anchors =
            SealedAnchorProvider::new(&UnusedProvider).with_native_ray(&ray, &public_range);
        assert_public_request(
            anchors
                .anchor_longitude_deg(anchor, request, false, false)
                .unwrap_err(),
            request,
        );
        assert_eq!(*seen.borrow(), [request]);
    }
    let ray = |_: &StarData, _, _, _| Err(Error::SpkBodyNotFound { target: 599 });
    let wrap = |_, _: Error| -> Error { panic!("missing alias was promoted") };
    let anchors = SealedAnchorProvider::new(&UnusedProvider).with_native_ray(&ray, &wrap);
    assert_eq!(
        anchors.anchor_longitude_deg(DynamicAnchor::Spica, request, false, false),
        Err(Error::SpkBodyNotFound { target: 599 })
    );
}

#[test]
fn native_star_projection_keeps_unselected_defining_clocks_lazy() {
    use crate::constants::{
        FLG_EQUATORIAL, FLG_SIDEREAL, SIDBIT_ECL_T0, SIDBIT_USER_UT, SIDM_FAGAN_BRADLEY,
        SIDM_LAHIRI,
    };
    let star = resolve_star_ref("Regulus").unwrap();
    let raw =
        |_: i32, _: f64, _: bool, _: bool, _: bool, _: i32, _: Option<StarTopoSite>, _: bool| {
            Ok([10.0, 20.0, 1.0, 0.0, 0.0, 0.0])
        };
    let nutation = |_| -> Result<(f64, f64), Error> { panic!("fixed frame owns nutation") };
    let delta = |_| -> Result<f64, Error> { panic!("unused defining epoch consumed Delta-T") };
    for mode in [SIDM_LAHIRI, SIDM_FAGAN_BRADLEY] {
        let mut selection = crate::sidereal::SiderealSelection::default();
        selection.set_sidereal_mode(mode | SIDBIT_ECL_T0 | SIDBIT_USER_UT, 1e308, 0.0);
        let serving = SealedStarServing::new(&UnusedProvider, &nutation, &delta, selection)
            .with_default_timescale()
            .with_native_places(&raw, None);
        let (place, name, flags) = serving
            .fixstar_tt_wrapped(&star, J2000, FLG_EQUATORIAL | FLG_SIDEREAL, None)
            .unwrap();
        assert!(place.iter().all(|x| x.is_finite()));
        assert_eq!(name, star.canonical);
        assert_eq!(flags, FLG_EQUATORIAL | FLG_SIDEREAL);
    }
}

#[test]
fn native_star_projection_acquires_only_selected_defining_clock() {
    use crate::constants::{
        FLG_EQUATORIAL, FLG_SIDEREAL, SIDBIT_ECL_T0, SIDBIT_USER_UT, SIDM_USER, SIDM_VALENS_MOON,
    };
    let star = resolve_star_ref("Regulus").unwrap();
    let raw =
        |_: i32, _: f64, _: bool, _: bool, _: bool, _: i32, _: Option<StarTopoSite>, _: bool| {
            Ok([10.0, 20.0, 1.0, 0.0, 0.0, 0.0])
        };
    let nutation = |_| -> Result<(f64, f64), Error> { panic!("fixed frame owns nutation") };
    for (mode, stored, expected) in [
        (SIDM_USER, J2000, J2000),
        (SIDM_VALENS_MOON, 1e308, crate::sidereal::VALENS_MOON_T0_UT),
    ] {
        let reads = RefCell::new(Vec::new());
        let delta = |ut| {
            assert_eq!(ut, expected, "unselected defining epoch was consumed");
            reads.borrow_mut().push(ut);
            Ok(0.001)
        };
        let mut selection = crate::sidereal::SiderealSelection::default();
        selection.set_sidereal_mode(mode | SIDBIT_ECL_T0 | SIDBIT_USER_UT, stored, 0.0);
        let serving = SealedStarServing::new(&UnusedProvider, &nutation, &delta, selection)
            .with_default_timescale()
            .with_native_places(&raw, None);
        assert!(
            serving
                .fixstar_tt_wrapped(&star, J2000, FLG_EQUATORIAL | FLG_SIDEREAL, None)
                .is_ok()
        );
        assert_eq!(*reads.borrow(), [expected]);
    }
}

#[test]
fn native_star_projection_preserves_selected_model_failure() {
    use crate::constants::{
        FLG_EQUATORIAL, FLG_SIDEREAL, SIDBIT_ECL_T0, SIDBIT_USER_UT, SIDM_USER, SIDM_VALENS_MOON,
    };
    let star = resolve_star_ref("Regulus").unwrap();
    let raw =
        |_: i32, _: f64, _: bool, _: bool, _: bool, _: i32, _: Option<StarTopoSite>, _: bool| {
            Ok([10.0, 20.0, 1.0, 0.0, 0.0, 0.0])
        };
    let nutation = |_| -> Result<(f64, f64), Error> { panic!("fixed frame owns nutation") };
    let delta = |_| Ok(f64::INFINITY);
    for mode in [SIDM_USER, SIDM_VALENS_MOON] {
        let mut selection = crate::sidereal::SiderealSelection::default();
        selection.set_sidereal_mode(mode | SIDBIT_ECL_T0 | SIDBIT_USER_UT, J2000, 0.0);
        let serving = SealedStarServing::new(&UnusedProvider, &nutation, &delta, selection)
            .with_default_timescale()
            .with_native_places(&raw, None);
        assert!(matches!(
            serving.fixstar_tt_wrapped(&star, J2000, FLG_EQUATORIAL | FLG_SIDEREAL, None),
            Err(Error::NonFiniteArgument { value }) if value.is_infinite()
        ));
    }
}
