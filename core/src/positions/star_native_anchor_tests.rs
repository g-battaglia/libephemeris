//! Synthetic callback, Time-part and typed-error controls; no actual image is loaded.
use super::*;
use std::cell::RefCell;

fn raw_range(target: i32, sample: TwoPartEpoch) -> Error {
    Error::SpkOutOfRange {
        target,
        et: (sample.hi - crate::constants::J2000) * 86400.0 + sample.lo * 86400.0,
        start: 100.0,
        end: 200.0,
    }
}

#[test]
fn anchor_noaberr_uses_vector_difference_and_source_call_order() {
    let astro = [1.0, 2.0, 3.0];
    for (noaberr, nogdefl, wanted, order) in [
        (false, false, [20.0, 30.0, 40.0], vec![true]),
        (true, false, [16.0, 26.0, 36.0], vec![true, false]),
        (false, true, [5.0, 6.0, 7.0], vec![false]),
        (true, true, astro, vec![]),
    ] {
        let seen = RefCell::new(Vec::new());
        let got = native_anchor_correction_via(astro, noaberr, nogdefl, |default| {
            seen.borrow_mut().push(default);
            Ok(if default {
                [20.0, 30.0, 40.0]
            } else {
                [5.0, 6.0, 7.0]
            })
        })
        .unwrap();
        assert_eq!(got, wanted);
        assert_eq!(*seen.borrow(), order);
    }
    let raw = Error::SpkBodyNotFound { target: 599 };
    assert_eq!(
        native_anchor_correction_via(astro, true, false, |_| Err(raw.clone())),
        Err(raw)
    );
}

#[test]
fn native_anchor_uses_fresh_scalar_time_and_strict_original_supports() {
    let tt = crate::constants::J2000 + 0.125;
    let time = scalar_tt(tt).unwrap();
    let data = crate::sidereal::anchor_star("SPICA").unwrap();
    for (failed, target) in [
        (SUN, 10),
        (crate::constants::JUPITER, 599),
        (crate::constants::SATURN, 699),
    ] {
        let seen = RefCell::new(Vec::new());
        let at = |id, sample: TwoPartEpoch| {
            seen.borrow_mut().push((id, sample));
            if id == failed {
                return Err(raw_range(target, sample));
            }
            Ok((
                if id == EARTH {
                    [0.0; 3]
                } else {
                    [2.0, 1.0, 0.5]
                },
                [0.0; 3],
            ))
        };
        assert_eq!(
            native_anchor_ray_via(&at, &data, tt, false, false),
            Err(raw_range(target, time))
        );
        let reads = seen.borrow();
        assert_eq!(reads[0], (EARTH, time));
        assert_eq!(reads.last(), Some(&(failed, time)));
    }
    for toggles in [(true, true), (false, true)] {
        let reads = RefCell::new(Vec::new());
        let at = |id, sample: TwoPartEpoch| {
            assert_eq!(id, EARTH, "unused native deflector was consumed");
            reads.borrow_mut().push(sample);
            Ok(([0.0; 3], [0.0; 3]))
        };
        assert!(
            native_anchor_ray_via(&at, &data, tt, toggles.0, toggles.1)
                .unwrap()
                .iter()
                .all(|x| x.is_finite())
        );
        assert_eq!(*reads.borrow(), [time]);
    }
}

#[test]
fn native_anchor_closest_time_retains_split_subtraction_and_raw_error() {
    let tt = crate::constants::J2000 + 0.125;
    let data = crate::sidereal::anchor_star("SPICA").unwrap();
    let direction = crate::stars::star_icrs_initial_state_catalog_angles(&data).0;
    let reads = RefCell::new(Vec::new());
    let at = |id, time: TwoPartEpoch| {
        if id == EARTH {
            return Ok(([0.0; 3], [0.0; 3]));
        }
        assert_eq!(id, SUN);
        reads.borrow_mut().push(time);
        if reads.borrow().len() == 2 {
            return Err(raw_range(10, time));
        }
        Ok((direction.map(|v| v * 0.001), [0.0; 3]))
    };
    let got = native_anchor_ray_via(&at, &data, tt, true, false).unwrap_err();
    let times = reads.borrow();
    assert_eq!(times.len(), 2);
    assert_eq!(times[0], scalar_tt(tt).unwrap());
    assert!(times[1].value() < times[0].value());
    assert_eq!(got, raw_range(10, times[1]));
}

#[test]
fn main_sun_image_envelope_is_independent_of_observer_and_failed_support() {
    let mut registry = crate::spk::SpkRegistry::new();
    let mut main = crate::spk::moon_test_kernel_bytes();
    let descriptor = crate::spk::RECORD_BYTES + 3 * 8;
    main[descriptor + 16..descriptor + 20].copy_from_slice(&10_i32.to_le_bytes());
    registry.register_image("main-sun", main).unwrap();
    registry.register_body(SUN, "main-sun", 10).unwrap();
    let mut earth = crate::spk::moon_test_kernel_bytes();
    for index in 0..3 {
        let offset = descriptor + index * 40;
        earth[offset + 8..offset + 16].copy_from_slice(&100.0_f64.to_le_bytes());
    }
    earth[descriptor + 16..descriptor + 20].copy_from_slice(&399_i32.to_le_bytes());
    registry.register_image("separate-earth", earth).unwrap();
    registry
        .register_body(EARTH, "separate-earth", 399)
        .unwrap();
    let provider = LocalSpkStateProvider::new(&registry, ProviderMode::NativeJpl).unwrap();
    let request = -0.0;
    for center in [0, FLG_TOPOCTR, FLG_HELCTR, FLG_BARYCTR] {
        let raw = Error::SpkOutOfRange {
            target: 599,
            et: -1.0,
            start: 1.0,
            end: 2.0,
        };
        let got = provider.star_range_error(center, request, raw.clone());
        let expected = match provider.public_range_error(SUN, request, raw) {
            Error::EphemerisRange {
                requested_jd,
                start_jd,
                end_jd,
                ..
            } => Error::EphemerisRange {
                requested_jd,
                start_jd,
                end_jd,
                body: None,
            },
            other => panic!("unexpected: {other:?}"),
        };
        assert_eq!(got, expected);
        if let Error::EphemerisRange {
            requested_jd,
            end_jd,
            ..
        } = got
        {
            assert_eq!(requested_jd.to_bits(), request.to_bits());
            assert!(end_jd > crate::constants::J2000 + 1.0);
        }
    }
    let already_public = Error::EphemerisRange {
        requested_jd: 12.0,
        start_jd: 10.0,
        end_jd: 11.0,
        body: Some(EARTH),
    };
    assert_eq!(
        provider.star_range_error(0, request, already_public.clone()),
        already_public
    );
    let missing = Error::SpkBodyNotFound { target: 599 };
    assert_eq!(
        provider.star_range_error(0, request, missing.clone()),
        missing
    );
}

#[test]
fn fast_batch_per_star_range_does_not_use_shared_observer_wrapper() {
    let registry = crate::spk::SpkRegistry::new();
    let native = LocalSpkStateProvider::new(&registry, ProviderMode::NativeJpl).unwrap();
    let provider = ContextStarProvider {
        inner: ContextStarBackend::Native(&native),
    };
    let nutation = |_| Ok((0.0, 0.0));
    let delta = |_| Ok(0.0);
    let selection = crate::sidereal::SiderealSelection::default();
    let serving = SealedStarServing::new(&provider, &nutation, &delta, selection);
    let error = Error::SpkOutOfRange {
        target: 10,
        et: -1.0,
        start: 1.0,
        end: 2.0,
    };
    let at = |id, _| {
        if id == EARTH {
            Ok(([0.0; 3], [0.0; 3]))
        } else {
            Err(error.clone())
        }
    };
    let wrap = |_, _, _: Error| -> Error {
        panic!("per-star range was promoted as a shared-observer failure")
    };
    assert_eq!(
        native_batch_via(
            &at,
            &wrap,
            &serving,
            selection,
            &["Regulus"],
            2451545.125,
            2,
            None,
            false
        ),
        Err(error.clone())
    );
    assert_eq!(
        native_batch_via(
            &at,
            &wrap,
            &serving,
            selection,
            &["Regulus"],
            2451545.125,
            2,
            None,
            true
        )
        .unwrap(),
        vec![None]
    );
}
