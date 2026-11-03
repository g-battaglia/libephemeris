// SPDX-License-Identifier: AGPL-3.0
//! Focused tests for the parent position-routing module.

use super::*;
#[test]
fn selected_naif_target_is_independent_of_library_sun_alias() {
    let mut registry = SpkRegistry::new();
    let offset = crate::spk::RECORD_BYTES + 3 * 8 + 2 * 40 + 16;
    for (body, target, file) in [(SUN, 399_i32, "wrong-main"), (17, 10, "actual-main")] {
        let mut image = crate::spk::moon_test_kernel_bytes();
        image[offset..offset + 4].copy_from_slice(&target.to_le_bytes());
        registry.register_image(file, image).unwrap();
        registry.register_body(body, file, target).unwrap();
    }
    let selected = native_main_image_registration(&registry).unwrap();
    assert_eq!(
        (selected.body, selected.naif_id, selected.kernel.as_str()),
        (17, 10, "actual-main")
    );
    assert_eq!(
        local_jpl_target(&JplLeg::new(&registry), SUN).unwrap(),
        selected.body
    );
}
#[test]
fn same_target_order_and_missing_owner_follow_the_core_selector() {
    let mut registry = SpkRegistry::new();
    assert_eq!(
        native_main_image_registration(&registry),
        Err(Error::SpkUnregisteredBody { body: SUN })
    );
    let offset = crate::spk::RECORD_BYTES + 3 * 8 + 2 * 40 + 16;
    let mut image = crate::spk::moon_test_kernel_bytes();
    image[offset..offset + 4].copy_from_slice(&10_i32.to_le_bytes());
    registry.register_image("shared-main", image).unwrap();
    registry.register_body(17, "shared-main", 10).unwrap();
    registry.register_body(4, "shared-main", 10).unwrap();
    assert_eq!(native_main_image_registration(&registry).unwrap().body, 4);
}
