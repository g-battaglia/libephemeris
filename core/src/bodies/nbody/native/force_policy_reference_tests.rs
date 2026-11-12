// SPDX-License-Identifier: AGPL-3.0
//! Explicit long-arc diagnostic; production always retains nonself sources.
use super::super::*;
use super::verified_force_data;
use serde_json::{Value, json};

fn words(value: &Value) -> [f64; 6] {
    let values = value.as_array().expect("six captured words");
    assert_eq!(values.len(), 6);
    std::array::from_fn(|index| {
        let value = f64::from_bits(values[index].as_u64().unwrap());
        assert!(value.is_finite());
        value
    })
}

fn poles(epoch: [f64; 2]) -> Result<HarmonicPoles, Error> {
    super::super::super::poles::PoleModel::IcrfTerrestrialIauSolar.at(epoch)
}

#[test]
#[ignore = "explicit exact-seed force-policy diagnosis; authenticated DE440/N16 required"]
fn transport_original_seeds_with_selected_and_major_only_datasets() {
    let input = std::env::var("LIBEPHEMERIS_FORCE_POLICY_INPUT").unwrap();
    let output =
        std::path::PathBuf::from(std::env::var("LIBEPHEMERIS_FORCE_POLICY_OUTPUT").unwrap());
    assert!(!output.exists(), "never overwrite a prior diagnosis");
    let packet: Value = serde_json::from_reader(std::fs::File::open(input).unwrap()).unwrap();
    assert_eq!(packet["schema"], "n16-force-policy-long-arc-v1");
    let requests = packet["requests"].as_array().unwrap();
    assert!(matches!(requests.len(), 5 | 7));
    let inventory = [
        (10016, 2451545.0),
        (10087, 2451545.0),
        (18, 2451545.0),
        (19, 2451545.0),
        (10087, 2469807.5),
        (10016, 2460310.5),
        (10704, 2460310.5),
    ];
    for (request, (body, endpoint)) in requests.iter().zip(inventory) {
        assert_eq!(request["body"].as_i64(), Some(body));
        assert_eq!(request["endpoint"].as_f64(), Some(endpoint));
    }
    let data = verified_force_data();
    let planet_path =
        std::path::PathBuf::from(std::env::var("LIBEPHEMERIS_TEST_DATA_DIR").unwrap())
            .join("assist/linux_p1550p2650.440");
    // This explicit diagnostic owner omits the asteroid image to reproduce
    // the observed v3 family policy. No runtime switch or fallback is added.
    let major_only = Arc::new(Ephemeris::open(&planet_path, None).unwrap());
    let mut results = Vec::new();
    for request in requests {
        let seed = words(&request["state_bits"]);
        let initial = State {
            position: seed[..3].try_into().unwrap(),
            velocity: seed[3..].try_into().unwrap(),
        };
        let epoch = request["epoch"].as_f64().unwrap();
        let endpoint = request["endpoint"].as_f64().unwrap();
        let target = i32::try_from(request["target"].as_i64().unwrap()).unwrap();
        let mut variants = Vec::new();
        for (name, owner) in [
            ("selected-n16-self-only", &data),
            ("major-only-diagnostic", &major_only),
        ] {
            let mut trajectory = ForcedTrajectory::new(
                Arc::clone(owner),
                initial,
                [epoch, 0.0],
                None,
                None,
                Some(target),
                poles,
            )
            .unwrap();
            let result = trajectory.advance(endpoint - epoch).unwrap();
            assert_eq!(trajectory.time(), endpoint - epoch);
            variants.push(json!({"model": name, "state_ssb": result.position.into_iter().chain(result.velocity).collect::<Vec<_>>()}));
        }
        results.push(json!({"request": request, "variants": variants}));
    }
    // Recheck all selected input hashes after the entire bounded experiment.
    verified_force_data();
    let file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(output)
        .unwrap();
    serde_json::to_writer_pretty(
        file,
        &json!({
            "schema": "n16-force-policy-long-arc-v1", "automatic_acceptance": false,
            "scope": "exact-old-seed force-family transport; not physical endpoint acceptance",
            "requests": results,
        }),
    )
    .unwrap();
}
