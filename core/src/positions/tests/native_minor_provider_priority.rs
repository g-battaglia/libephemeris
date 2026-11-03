// SPDX-License-Identifier: AGPL-3.0
//! Focused tests for the parent position-routing module.

use super::*;
use crate::engine_api::MinorAssistAcquisition;
use std::sync::atomic::{AtomicUsize, Ordering};

fn constant_target(radius_au: f64) -> SpkRegistry {
    let mut image = crate::spk::moon_test_kernel_bytes();
    let descriptor = crate::spk::RECORD_BYTES + 3 * 8 + 2 * 40;
    let address = i32::from_le_bytes(image[descriptor + 32..descriptor + 36].try_into().unwrap());
    image[descriptor + 36..descriptor + 40].copy_from_slice(&(address + 8).to_le_bytes());
    let offset = (address as usize - 1) * 8;
    for (index, value) in [
        2_500_000.0,
        2_500_000.0,
        radius_au * 149_597_870.7,
        0.0,
        0.0,
        0.0,
        5_000_000.0,
        5.0,
        1.0,
    ]
    .into_iter()
    .enumerate()
    {
        image[offset + index * 8..offset + (index + 1) * 8].copy_from_slice(&value.to_le_bytes());
    }
    let mut registry = SpkRegistry::new();
    registry.register_image("constant-minor", image).unwrap();
    registry.register_body(CERES, "constant-minor", 3).unwrap();
    registry
}

struct Supplier {
    refreshed: SpkRegistry,
    refreshes: AtomicUsize,
}
impl MinorAssistAcquisition for Supplier {
    fn refreshed_spk(
        &self,
        _: i32,
        _: f64,
        _: &crate::context::EphemerisContext,
    ) -> Result<Option<SpkRegistry>> {
        self.refreshes.fetch_add(1, Ordering::SeqCst);
        Ok(Some(self.refreshed.clone()))
    }
    fn contains(&self, _: i32) -> Result<bool> {
        panic!("served SPK precedes ASSIST")
    }
    fn available(&self) -> Result<bool> {
        panic!("served SPK precedes ASSIST")
    }
    fn state_at(&self, _: i32, _: f64) -> Result<crate::nbody::PropagationResult> {
        panic!("served SPK precedes ASSIST")
    }
}

fn calculate(context: &crate::context::EphemerisContext, supplier: &Supplier) -> PositionResult {
    let request = ServingRequest {
        ctx: context,
        delta_t_model: crate::time::DeltaTModel::EspenakMeeus,
        iers_enabled: false,
        observer: None,
        tier_policy: TierPolicy::AttachOrder,
        skip_closed: false,
    };
    calc_native_jpl_entry_with_minor_assist(
        None,
        None,
        Some(supplier),
        request,
        FRAME_J2000 + 0.25,
        CERES,
        FLG_BARYCTR | FLG_TRUEPOS | FLG_J2000 | FLG_NONUT | FLG_XYZ | FLG_ICRS,
        false,
    )
    .unwrap()
    .0
}

#[test]
fn planetary_support_depends_on_the_actual_minor_origin_and_target() {
    let empty = SpkRegistry::new();
    let leg = JplLeg::new(&empty);
    assert!(!minor_requires_main_image(&leg, CERES, FRAME_J2000 + 0.25, FLG_HELCTR).unwrap());
    assert!(minor_requires_main_image(&leg, CERES, FRAME_J2000 + 0.25, 0).unwrap());
    let registered = constant_target(2.0);
    let leg = JplLeg::new(&registered);
    // The target is SSB-centered. Barycentric output is self-contained;
    // heliocentric/geocentric outputs need an actual Sun/Earth provider.
    assert!(!minor_requires_main_image(&leg, CERES, FRAME_J2000 + 0.25, FLG_BARYCTR).unwrap());
    assert!(minor_requires_main_image(&leg, CERES, FRAME_J2000 + 0.25, FLG_HELCTR).unwrap());
    assert!(minor_requires_main_image(&leg, CERES, FRAME_J2000 + 0.25, 0).unwrap());
}

#[test]
fn refreshed_spk_is_used_before_assist_and_preserves_the_original_snapshot() {
    let context = crate::context::EphemerisContext::new()
        .with_mode(ProviderMode::NativeJpl)
        .with_network_policy(Some(crate::context::NetworkPolicy::Sealed));
    let supplier = Supplier {
        refreshed: constant_target(2.0),
        refreshes: AtomicUsize::new(0),
    };
    let position = calculate(&context, &supplier);
    assert!((position[0] - 2.0).abs() < 1e-12, "{position:?}");
    assert_eq!(supplier.refreshes.load(Ordering::SeqCst), 1);
    assert!(context.spk_registry().is_empty());
}

#[test]
fn registered_spk_wins_without_provisioning_or_numerical_fallback() {
    let mut context = crate::context::EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
    *context.spk_registry_mut() = constant_target(2.0);
    let supplier = Supplier {
        refreshed: constant_target(4.0),
        refreshes: AtomicUsize::new(0),
    };
    let position = calculate(&context, &supplier);
    assert!((position[0] - 2.0).abs() < 1e-12, "{position:?}");
    assert_eq!(supplier.refreshes.load(Ordering::SeqCst), 0);
}
