//! Resource scope and physical invariants for model-only horizon conversions.

use super::*;
use crate::constants::{FLG_SPEED, FLG_SWIEPH, SUN};
use crate::engine_api::{Celestial, Horizontal};
use crate::events::{ECL2HOR, EQU2HOR, HOR2ECL, HOR2EQU};
use crate::leb::{IndexedLimits, LebReadPolicy};
use crate::time::{DeltaTModel, SIDEREAL_MODERN_END_JD, SIDEREAL_MODERN_START_JD};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

struct TestDirectory(PathBuf);

impl TestDirectory {
    fn new() -> Self {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let path = std::env::temp_dir().join(format!(
            "libephemeris-observer-model-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&path).unwrap();
        Self(path)
    }

    fn image(&self) -> PathBuf {
        self.0.join("observer.leb2")
    }

    fn install_image(&self) {
        std::fs::copy(crate::test_data::core_leb_path(), self.image()).unwrap();
    }
}

impl Drop for TestDirectory {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn context(path: &Path) -> EphemerisContext {
    let mut context = EphemerisContext::new()
        .with_mode(ProviderMode::SealedLeb)
        .with_iers_delta_t(Some(false))
        .with_iers_auto_download(Some(false));
    context.set_network_policy(Some(crate::context::NetworkPolicy::Sealed));
    context.set_leb_file(Some(path.to_path_buf()));
    context
}

fn horizontal(session: &Session) -> Horizontal {
    session
        .azalt(
            2451545.0, EQU2HOR, 12.5, 41.9, 250.0, 0.0, 15.0, 132.0, -12.0,
        )
        .unwrap()
}

fn horizontal_bits(value: Horizontal) -> [u64; 3] {
    [value.azimuth, value.true_altitude, value.apparent_altitude].map(f64::to_bits)
}

fn celestial_bits(value: Celestial) -> [u64; 2] {
    [value.longitude_or_ra, value.latitude_or_dec].map(f64::to_bits)
}

fn assert_unacquired(runtime: &Runtime, session: &Session) {
    let state = runtime.lock();
    assert!(state.reader.is_empty());
    assert!(state.reader_paths.is_none());
    assert!(state.reader_policy.is_none());
    assert!(state.reader_manifest_matches.is_empty());
    assert!(session.context().leb_files().is_empty());
}

fn assert_same_error(actual: Error, expected: &Error) {
    let actual = actual.primary_error();
    let expected = expected.primary_error();
    assert_eq!(
        std::mem::discriminant(actual),
        std::mem::discriminant(expected)
    );
    assert_eq!(actual.to_string(), expected.to_string());
}

#[test]
fn observer_models_preserve_configured_clocks_frames_and_atmosphere() {
    let directory = TestDirectory::new();
    directory.install_image();
    for (model, user_delta, iers) in [
        (DeltaTModel::Smh2016, None, false),
        (DeltaTModel::EspenakMeeus, None, false),
        (DeltaTModel::EspenakMeeus, Some(0.125), true),
        (DeltaTModel::Smh2016, None, true),
    ] {
        let mut configured = context(&directory.image())
            .with_delta_t_model(Some(model))
            .with_user_delta_t(user_delta)
            .with_tidal_acceleration(Some(-26.0))
            .with_iers_delta_t(Some(iers))
            .with_lapse_rate(Some(0.005));
        configured.attach_iers_rows(crate::iers::IersDataSnapshot::from_texts(
            None,
            None,
            Some("2000 1 1 63.0\n2000 2 1 66.1\n2000 3 1 71.9\n"),
        ));
        let without_iers = Session::from_context(configured.clone().with_iers_delta_t(Some(false)));
        let complete_runtime = Runtime::from_context(configured.clone());
        let complete = complete_runtime.shared_session().unwrap();
        let runtime = Runtime::from_context(configured);
        let geometry = runtime.model_session().unwrap();
        assert_unacquired(&runtime, &geometry);
        if let Some(days) = user_delta {
            assert_eq!(
                geometry.deltat(2451545.0).unwrap().to_bits(),
                days.to_bits()
            );
        } else if iers {
            assert_ne!(
                geometry.deltat(2451545.0).unwrap().to_bits(),
                without_iers.deltat(2451545.0).unwrap().to_bits()
            );
        }
        for jd in [
            SIDEREAL_MODERN_START_JD - 0.25,
            SIDEREAL_MODERN_START_JD,
            2451545.0,
            SIDEREAL_MODERN_END_JD,
            SIDEREAL_MODERN_END_JD + 0.25,
        ] {
            assert_eq!(
                geometry.deltat(jd).unwrap().to_bits(),
                complete.deltat(jd).unwrap().to_bits()
            );
            for flag in [ECL2HOR, EQU2HOR, 37] {
                for (lat, alt, pressure, temp) in
                    [(41.9, 250.0, 0.0, 15.0), (-35.0, -75.0, 900.0, -5.0)]
                {
                    let actual = geometry
                        .azalt(jd, flag, 12.5, lat, alt, pressure, temp, 132.0, -12.0)
                        .unwrap();
                    let expected = complete
                        .azalt(jd, flag, 12.5, lat, alt, pressure, temp, 132.0, -12.0)
                        .unwrap();
                    assert_eq!(horizontal_bits(actual), horizontal_bits(expected));
                    for reverse in [HOR2ECL, HOR2EQU, 37] {
                        assert_eq!(
                            celestial_bits(
                                geometry
                                    .azalt_rev(
                                        jd,
                                        reverse,
                                        12.5,
                                        lat,
                                        alt,
                                        actual.azimuth,
                                        actual.true_altitude
                                    )
                                    .unwrap()
                            ),
                            celestial_bits(
                                complete
                                    .azalt_rev(
                                        jd,
                                        reverse,
                                        12.5,
                                        lat,
                                        alt,
                                        expected.azimuth,
                                        expected.true_altitude
                                    )
                                    .unwrap()
                            )
                        );
                    }
                }
            }
        }
        assert_unacquired(&runtime, &geometry);
    }
}

#[test]
fn observer_models_defer_missing_and_corrupt_body_images() {
    let directory = TestDirectory::new();
    for corrupt in [false, true] {
        if corrupt {
            std::fs::write(directory.image(), b"invalid LEB image").unwrap();
        }
        let configured = context(&directory.image());
        let expected_error = Runtime::from_context(configured.clone())
            .shared_session()
            .unwrap_err();
        let runtime = Runtime::from_context(configured);
        let geometry = runtime.model_session().unwrap();
        assert!(horizontal(&geometry).azimuth.is_finite());
        assert_unacquired(&runtime, &geometry);
        assert_same_error(
            runtime
                .calc_ut(2451545.0, SUN, FLG_SWIEPH | FLG_SPEED)
                .unwrap_err(),
            &expected_error,
        );
        assert_same_error(runtime.leb_inventory().unwrap_err(), &expected_error);
        // Failed acquisition may retire the cached Session. The earlier
        // caller remains valid, and a new model snapshot keeps its semantics.
        let after_failure = runtime.model_session().unwrap();
        assert_eq!(
            horizontal_bits(horizontal(&after_failure)),
            horizontal_bits(horizontal(&geometry))
        );
        assert!(runtime.lock().reader.is_empty());
        assert!(after_failure.context().leb_files().is_empty());
        assert!(geometry.context().leb_files().is_empty());
    }
}

#[test]
fn observer_models_defer_valid_zero_budgets_until_body_acquisition() {
    let directory = TestDirectory::new();
    directory.install_image();
    let mut configured = context(&directory.image());
    let policy = LebReadPolicy::Indexed(IndexedLimits {
        metadata_bytes: 0,
        cache_bytes: 0,
        working_bytes: 0,
        ..IndexedLimits::default()
    });
    configured.set_leb_read_policy(Some(policy)).unwrap();
    let expected_error = Runtime::from_context(configured.clone())
        .shared_session()
        .unwrap_err();
    assert!(matches!(
        expected_error.primary_error(),
        Error::ResourceLimit { .. }
    ));
    let runtime = Runtime::from_context(configured);
    let geometry = runtime.model_session().unwrap();
    assert_eq!(geometry.context().leb_read_policy(), policy);
    assert!(horizontal(&geometry).azimuth.is_finite());
    assert_unacquired(&runtime, &geometry);
    assert_same_error(runtime.leb_inventory().unwrap_err(), &expected_error);
    assert_same_error(
        runtime
            .calc_ut(2451545.0, SUN, FLG_SWIEPH | FLG_SPEED)
            .unwrap_err(),
        &expected_error,
    );
}

#[test]
fn observer_models_still_reject_invalid_reader_configuration() {
    let sources = ConfigurationSources::load(
        std::collections::BTreeMap::from([(
            "LIBEPHEMERIS_LEB_READER".into(),
            "invalid-reader".into(),
        )]),
        None,
        None,
    );
    let runtime = Runtime::with_sources(sources);
    let expected = runtime.shared_session().unwrap_err();
    assert!(matches!(
        expected.primary_error(),
        Error::InvalidLebReader { .. }
    ));
    assert_same_error(runtime.model_session().unwrap_err(), &expected);
    assert!(runtime.lock().cached_session.is_none());
}

#[test]
fn observer_models_promote_once_and_reuse_the_retained_complete_snapshot() {
    let directory = TestDirectory::new();
    directory.install_image();
    let runtime = Runtime::from_context(context(&directory.image()));
    let geometry = runtime.model_session().unwrap();
    let expected = horizontal_bits(horizontal(&geometry));
    assert_unacquired(&runtime, &geometry);
    runtime
        .calc_ut(2451545.0, SUN, FLG_SWIEPH | FLG_SPEED)
        .unwrap();
    let inventory = runtime.leb_inventory().unwrap();
    let complete = runtime.shared_session().unwrap();
    assert_eq!(complete.context().leb_files().len(), 1);
    assert!(!Arc::ptr_eq(&geometry, &complete));
    std::fs::remove_file(directory.image()).unwrap();
    let retained = runtime.model_session().unwrap();
    assert!(Arc::ptr_eq(&complete, &retained));
    assert_eq!(horizontal_bits(horizontal(&retained)), expected);
    assert_eq!(runtime.leb_inventory().unwrap(), inventory);
    assert!(geometry.context().leb_files().is_empty());
    assert_eq!(horizontal_bits(horizontal(&geometry)), expected);

    runtime.set_leb_file(Some(directory.0.join("replacement-missing.leb2")));
    let refreshed = runtime.model_session().unwrap();
    assert!(!Arc::ptr_eq(&complete, &refreshed));
    assert_eq!(horizontal_bits(horizontal(&refreshed)), expected);
    assert!(matches!(
        runtime.leb_inventory(),
        Err(Error::LebFileNotFound { .. })
    ));
}

#[test]
fn observer_models_preserve_round_trip_and_consumed_domain_order() {
    let directory = TestDirectory::new();
    let runtime = Runtime::from_context(context(&directory.image()));
    let geometry = runtime.model_session().unwrap();
    for (lon, lat, ra, dec) in [(12.5, 41.9, 132.0, -12.0), (-75.0, -64.0, 359.0, 30.0)] {
        let horizon = geometry
            .azalt(2451545.0, EQU2HOR, lon, lat, 250.0, 1013.25, 15.0, ra, dec)
            .unwrap();
        let sky = geometry
            .azalt_rev(
                2451545.0,
                HOR2EQU,
                lon,
                lat,
                250.0,
                horizon.azimuth,
                horizon.true_altitude,
            )
            .unwrap();
        let ra_error = (sky.longitude_or_ra - ra + 180.0).rem_euclid(360.0) - 180.0;
        assert!(ra_error.abs() < 1e-9);
        assert!((sky.latitude_or_dec - dec).abs() < 1e-9);
    }
    // Reverse altitude is ignored physically, but its finite-input gate remains.
    assert!(matches!(
        geometry.azalt_rev(2451545.0, HOR2EQU, 12.5, 41.9, f64::INFINITY, 20.0, 30.0),
        Err(Error::NonFiniteArgument { value }) if value == f64::INFINITY
    ));
    assert!(matches!(
        geometry.azalt(2451545.0, EQU2HOR, 12.5, 41.9, 250.0, f64::INFINITY, 15.0, 132.0, -12.0),
        Err(Error::NonFiniteArgument { value }) if value == f64::INFINITY
    ));
    // The consumed model clock still precedes coordinate-domain validation.
    runtime.set_delta_t_userdef(Some(f64::INFINITY));
    let invalid_clock = runtime.model_session().unwrap();
    assert!(matches!(
        invalid_clock.azalt(2451545.0, EQU2HOR, f64::NAN, 41.9, 250.0, 0.0, 15.0, 132.0, -12.0),
        Err(Error::NonFiniteArgument { value }) if value == f64::INFINITY
    ));
    assert!(geometry.deltat(2451545.0).unwrap().is_finite());
    assert_unacquired(&runtime, &invalid_clock);
}
