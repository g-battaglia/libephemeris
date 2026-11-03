// SPDX-License-Identifier: AGPL-3.0
//! Focused tests for the parent position-routing module.

use super::*;
#[test]
fn selected_center_directory_failure_retains_its_actual_image() {
    let mut image = crate::spk::range_capture_chain_image(60.0, 600.0);
    // The target is valid at 200s. Only the center's real directory is
    // shorter than its descriptor; this test has one physical refusal.
    let descriptor_end = crate::spk::RECORD_BYTES + 24 + 8;
    image[descriptor_end..descriptor_end + 8].copy_from_slice(&600.0_f64.to_le_bytes());
    let mut registry = SpkRegistry::new();
    registry.register_image("owned-chain.bsp", image).unwrap();
    registry
        .register_body(EARTH, "owned-chain.bsp", 399)
        .unwrap();
    let leg = JplLeg::new(&registry);
    let directory = PublicModelDirectory {
        leg: &leg,
        raw: std::cell::RefCell::new(None),
    };
    let failure = directory
        .tdb(EARTH, 2451545.0, 200.0 / 86400.0)
        .unwrap_err();
    assert!(matches!(failure, Error::SpkOutOfRange { .. }));
    {
        let raw = directory.raw.borrow();
        let raw = raw.as_ref().unwrap();
        assert_eq!((raw.descriptor.target, raw.descriptor.center), (3, 0));
        assert_eq!(raw.segment_ordinal, 0);
    }
    let ordinary = registry
        .eval_geometric_ssb_au_two_part_jd(EARTH, 2451545.0, 15.0 / 86400.0)
        .unwrap();
    let acquired = directory.tdb(EARTH, 2451545.0, 15.0 / 86400.0).unwrap();
    assert_eq!(ordinary.0.map(f64::to_bits), acquired.0.map(f64::to_bits));
    assert_eq!(ordinary.1.map(f64::to_bits), acquired.1.map(f64::to_bits));
    assert!(directory.raw.borrow().is_none());
}

#[test]
fn actual_raw_payload_survives_promotion_overlay_and_context_update() {
    let mut image = crate::spk::moon_test_kernel_bytes();
    let offset = crate::spk::RECORD_BYTES + 3 * 8 + 2 * 40 + 16;
    image[offset..offset + 4].copy_from_slice(&10_i32.to_le_bytes());
    // The shared fixture is index-only. Give the selected Sun segment one
    // valid constant type2 record so decoding precedes the time refusal.
    let initial = i32::from_le_bytes(image[offset + 16..offset + 20].try_into().unwrap());
    image[offset + 20..offset + 24].copy_from_slice(&(initial + 8).to_le_bytes());
    // Keep the descriptor wider than its numerical directory so the
    // retained owner belongs to a genuinely selected, refused record.
    image[offset - 8..offset].copy_from_slice(&10_000_000.0_f64.to_le_bytes());
    let data = (initial as usize - 1) * 8;
    for (index, word) in [
        2_500_000.0_f64,
        2_500_000.0,
        1.0,
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
        image[data + index * 8..data + (index + 1) * 8].copy_from_slice(&word.to_le_bytes());
    }
    let mut registry = SpkRegistry::new();
    registry
        .register_image("/owned/summary.bsp", image)
        .unwrap();
    registry
        .register_body(SUN, "/owned/summary.bsp", 10)
        .unwrap();
    let leg = JplLeg::new(&registry);
    let directory = PublicModelDirectory {
        leg: &leg,
        raw: std::cell::RefCell::new(None),
    };
    let native = directory.state(SUN, 2451545.0 + 100.0, 0.125).unwrap_err();
    assert!(matches!(native, Error::SpkOutOfRange { .. }));
    let captured = directory.raw.borrow().as_ref().unwrap().clone();
    let error = directory.promote(2451545.0 + 100.0, CUPIDO, native);
    assert_eq!(error.backend_range().unwrap(), &captured);
    assert!(directory.raw.borrow().is_none());
    let error = error.with_cause(Error::SpkBodyNotFound { target: 999 });
    assert_eq!(error.backend_range().unwrap(), &captured);
    drop(registry);
    assert_eq!(error.backend_range().unwrap(), &captured);
}
#[test]
fn a_later_unregistered_read_cannot_reuse_an_ignored_range_payload() {
    let registry = SpkRegistry::new();
    let leg = JplLeg::new(&registry);
    let directory = PublicModelDirectory {
        leg: &leg,
        raw: std::cell::RefCell::new(None),
    };
    assert!(directory.tdb(999, 2451545.0, 0.0).is_err());
    assert!(directory.raw.borrow().is_none());
}
