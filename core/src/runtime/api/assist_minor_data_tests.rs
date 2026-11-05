// SPDX-License-Identifier: AGPL-3.0
//! Explicit real-data derivative diagnostic for the held-state runtime model.
//! Higher-order samples use different, exactly representable TT grids. They
//! do not reintegrate side epochs or replace the original component budgets.
use super::*;
use sha2::{Digest, Sha256};
use std::io::Read;
use std::sync::Arc;

#[test]
#[ignore = "explicit authenticated DE440/N16 files required"]
fn geocentric_range_rates_agree_with_independent_sixth_order_grids() {
    let root = std::path::PathBuf::from(std::env::var("LIBEPHEMERIS_TEST_DATA_DIR").unwrap());
    for (name, expected) in [
        (
            "de440.bsp",
            "a4ce9bf9b3282becc9f4b2ac3cebe03a2ae7599981aabd7265fd8482fff7c4b5",
        ),
        (
            "assist/linux_p1550p2650.440",
            "29915576d0a6555766b99485ac3056ee415e86df4fce282611c31afb329ad062",
        ),
        (
            "assist/sb441-n16.bsp",
            "919d612ce3c72a78fc7158f9120156542d0f21e6b8b052e4c1339c759747fd90",
        ),
    ] {
        let mut image = std::fs::File::open(root.join(name)).unwrap();
        let mut hash = Sha256::new();
        let mut buffer = [0u8; 65536];
        loop {
            let count = image.read(&mut buffer).unwrap();
            if count == 0 {
                break;
            }
            hash.update(&buffer[..count]);
        }
        let actual: String = hash
            .finalize()
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect();
        assert_eq!(actual, expected, "{name}");
    }
    // Isolate exactly the declared primary image, not HOME companion files.
    let mut context = crate::context::EphemerisContext::new()
        .with_mode(crate::context::ProviderMode::NativeJpl)
        .with_ephemeris_file(Some("de440.bsp"))
        .with_jpl_data_dir(Some(root.join("absent-derivative-test-companions")))
        .with_network_policy(Some(crate::context::NetworkPolicy::Sealed));
    context.set_jpl_dir(Some(root.clone()));
    let main = super::super::MainImageLifecycle::new()
        .acquire(&context)
        .unwrap();
    let session = Session::from_context(main.native_provider_context().unwrap());
    let supplier = super::super::NativeMinorAssist::with_session(Arc::new(
        crate::nbody::AssistSession::new(crate::nbody::AssistEphemConfig {
            planets_file: Some(root.join("assist/linux_p1550p2650.440")),
            asteroids_file: Some(root.join("assist/sb441-n16.bsp")),
            data_dir: None,
        }),
    ));
    // Every curated body, including the observed small-rate body17066, uses
    // the same supplier and stencil. No body-specific runtime route is added.
    let bodies = [
        15, 16, 17, 18, 19, 20, 146199, 100377, 146108, 146472, 38978, 100482, 60000, 17066, 18405,
        20199, 235088, 30000, 109942, 10010, 10704, 10511, 10052, 10087, 10016, 10433, 11221,
        11566, 11685, 10080, 10055, 11181, 10944, 14179, 35143, 111955, 172173,
    ];
    for body in bodies {
        let central = supplier.state_at(body, J2000).unwrap();
        for modifiers in [0, FLG_TRUEPOS] {
            let flags = FLG_SPEED | modifiers;
            let actual = session
                .assist_minor_ecliptic(&supplier, body, J2000, flags)
                .unwrap()[5];
            let radius = |offset| {
                let tt = J2000 + offset;
                session
                    .assist_minor_position_from_state(&central.shifted(tt).unwrap(), tt, flags)
                    .unwrap()[2]
            };
            let old_step = 1.0 / 86400.0;
            let old = (radius(old_step) - radius(-old_step)) / (2.0 * old_step);
            let references = [1.0 / 32.0, 1.0 / 16.0, 1.0 / 8.0].map(|h| {
                (45.0 * (radius(h) - radius(-h)) - 9.0 * (radius(2.0 * h) - radius(-2.0 * h))
                    + radius(3.0 * h)
                    - radius(-3.0 * h))
                    / (60.0 * h)
            });
            for expected in references {
                let limit = expected.abs() * 1e-6;
                assert!(
                    (actual - expected).abs() <= limit,
                    "body={body}, flags={flags}: new={actual:.17e}, sixth={expected:.17e}, limit={limit:.17e}"
                );
            }
            eprintln!(
                "{}",
                serde_json::json!({"body":body,"flags":flags,
                "new_rate":actual,"old_one_second_rate":old,"sixth_order_rates":references})
            );
        }
    }
}
