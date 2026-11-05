//! Synthetic main-image ownership and lifecycle controls.
use super::*;
use crate::context::NetworkPolicy;
use std::cell::Cell;

#[test]
fn failed_acquisition_retries_and_only_success_is_shared() {
    let context = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
    let cache = NativeStarMainCache::default();
    let calls = Cell::new(0);
    let failure = Error::NetworkSealed {
        purpose: "synthetic main acquisition".into(),
    };
    let first = cache.snapshot_with(&context, None, |_| {
        calls.set(calls.get() + 1);
        Err(failure.clone())
    });
    assert_eq!(first.unwrap_err(), failure);
    assert!(cache.ready().is_none());
    let second = cache
        .snapshot_with(&context, None, |context| {
            calls.set(calls.get() + 1);
            Ok(context.clone())
        })
        .unwrap();
    let third = cache
        .snapshot_with(&context, None, |_| panic!("cached success must not reopen"))
        .unwrap();
    assert!(Arc::ptr_eq(&second, &third));
    assert_eq!(calls.get(), 2);
}

#[test]
fn supplied_planetary_registry_bypasses_configuration_and_acquisition() {
    let mut context = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
    let mut image = crate::spk::moon_test_kernel_bytes();
    image[1064..1068].copy_from_slice(&10_i32.to_le_bytes());
    image[1068..1072].copy_from_slice(&0_i32.to_le_bytes());
    context
        .spk_registry_mut()
        .register_image("synthetic-main", image)
        .unwrap();
    context
        .spk_registry_mut()
        .register_body(10, "synthetic-main", 10)
        .unwrap();
    let resolve = |_: &EphemerisContext, _| -> Result<EphemerisContext> {
        panic!("supplied registry must not consume configuration")
    };
    let cache = NativeStarMainCache::default();
    let selected = cache
        .snapshot_with(&context, Some(&resolve), |_| {
            panic!("supplied registry must not acquire")
        })
        .unwrap();
    assert!(!selected.spk_registry().is_empty());
    assert!(cache.ready().is_none());
}

#[test]
fn optional_images_do_not_replace_the_planetary_main() {
    let mut context = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
    context
        .spk_registry_mut()
        .register_image("optional-satellite", crate::spk::moon_test_kernel_bytes())
        .unwrap();
    let cache = NativeStarMainCache::default();
    let called = std::cell::Cell::new(false);
    let selected = cache
        .snapshot_with(&context, None, |current| {
            called.set(true);
            assert!(
                current
                    .spk_registry()
                    .kernel("optional-satellite")
                    .is_some()
            );
            Ok(current.clone())
        })
        .unwrap();
    assert!(called.get());
    assert!(
        selected
            .spk_registry()
            .kernel("optional-satellite")
            .is_some()
    );
}

#[test]
fn configuration_callback_is_unlocked_and_reentrant_success_wins() {
    let context = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
    let cache = NativeStarMainCache::default();
    let resolve = |context: &EphemerisContext, body| {
        assert_eq!(body, SUN);
        cache
            .snapshot_with(context, None, |context| Ok(context.clone()))
            .unwrap();
        Ok(context.clone())
    };
    let selected = cache
        .snapshot_with(&context, Some(&resolve), |_| {
            panic!("reentrant installed success wins")
        })
        .unwrap();
    assert!(Arc::ptr_eq(&selected, &cache.ready().unwrap()));
}

#[test]
fn configuration_failure_is_preserved_without_installing_an_owner() {
    let context = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
    let cache = NativeStarMainCache::default();
    let failure = Error::ConfigurationMissing {
        missing_config: "toy",
        suggestion: "toy",
    };
    let resolve = |_: &EphemerisContext, _| Err(failure.clone());
    assert_eq!(
        cache
            .snapshot_with(&context, Some(&resolve), |_| panic!("after refusal"))
            .unwrap_err(),
        failure
    );
    assert!(cache.ready().is_none());
    let replacement = |_: &EphemerisContext, _| Ok(EphemerisContext::new());
    assert!(matches!(
        cache.snapshot_with(&context, Some(&replacement), |_| panic!(
            "mode substitution"
        )),
        Err(Error::ConfigurationMissing {
            missing_config: "native_star_provider",
            ..
        })
    ));
}

#[test]
fn public_catalog_and_nonfinite_time_precede_binding_configuration() {
    let session = crate::engine_api::Session::from_context(
        EphemerisContext::new().with_mode(ProviderMode::NativeJpl),
    );
    let resolve = |_: &EphemerisContext, _| -> Result<EphemerisContext> {
        panic!("unconsumed configuration")
    };
    assert!(matches!(
        session.fixstar_ut_with_position_context("no-such-star-xyz", f64::NAN, 0, &resolve),
        Err(Error::StarNotFound { .. })
    ));
    assert!(matches!(
        session.fixstar_ut_with_position_context("Regulus", f64::NAN, 0, &resolve),
        Err(Error::NonFiniteArgument { .. })
    ));
}

#[test]
fn barycentric_public_place_still_consumes_main_selection() {
    let session = crate::engine_api::Session::from_context(
        EphemerisContext::new()
            .with_mode(ProviderMode::NativeJpl)
            .with_network_policy(Some(NetworkPolicy::Sealed)),
    );
    let failure = Error::ConfigurationMissing {
        missing_config: "toy main",
        suggestion: "toy",
    };
    let calls = Cell::new(0);
    let resolve = |_: &EphemerisContext, body| {
        calls.set(calls.get() + 1);
        assert_eq!(body, SUN);
        Err(failure.clone())
    };
    assert_eq!(
        session
            .fixstar_ut_with_position_context(
                "Regulus",
                2451545.125,
                crate::constants::FLG_BARYCTR,
                &resolve
            )
            .unwrap_err(),
        failure
    );
    assert_eq!(calls.get(), 1);
}

#[test]
fn automatic_local_barycentric_place_retains_owned_image_without_load_spk() {
    static ID: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let directory = std::env::temp_dir().join(format!(
        "lephem-star-main-{}-{}",
        std::process::id(),
        ID.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
    ));
    std::fs::create_dir(&directory).unwrap();
    struct Cleanup(std::path::PathBuf);
    impl Drop for Cleanup {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }
    let _cleanup = Cleanup(directory.clone());
    let mut image = crate::spk::moon_test_kernel_bytes();
    for (index, target) in [10i32, 399, 4].into_iter().enumerate() {
        let descriptor = 1048 + index * 40;
        image[descriptor + 16..descriptor + 20].copy_from_slice(&target.to_le_bytes());
        image[descriptor + 20..descriptor + 24].copy_from_slice(&0i32.to_le_bytes());
    }
    let path = directory.join("de440.bsp");
    std::fs::write(&path, image).unwrap();
    let context = EphemerisContext::new()
        .with_mode(ProviderMode::NativeJpl)
        .with_jpl_data_dir(Some(directory))
        .with_network_policy(Some(NetworkPolicy::Sealed));
    let session = crate::engine_api::Session::from_context(context.clone());
    let flags = crate::constants::FLG_BARYCTR | crate::constants::FLG_J2000;
    let first = session.fixstar_ut("Regulus", 2451545.125, flags).unwrap();
    assert!(context.spk_registry().is_empty());
    std::fs::remove_file(path).unwrap();
    let second = session
        .clone()
        .fixstar_ut("Regulus", 2451545.125, flags)
        .unwrap();
    assert_eq!(first, second);
    assert!(matches!(
        crate::engine_api::Session::from_context(context).fixstar_ut("Regulus", 2451545.125, flags),
        Err(Error::NetworkSealed { .. })
    ));
}
