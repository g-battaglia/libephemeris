// SPDX-License-Identifier: AGPL-3.0
//! Explicit same-seed step refinement, with every physical force retained.
use super::super::*;
use super::verified_force_data;
use serde_json::{Value, json};

#[test]
#[ignore = "bounded real-file Pallas2100 refinement; authenticated DE440/N16 required"]
fn pallas_endpoint_with_predeclared_maximum_steps() {
    let input = std::env::var("LIBEPHEMERIS_FORCE_POLICY_INPUT").unwrap();
    let output =
        std::path::PathBuf::from(std::env::var("LIBEPHEMERIS_PALLAS_STEPS_OUTPUT").unwrap());
    assert!(!output.exists(), "never overwrite a prior step diagnosis");
    let input_bytes = std::fs::read(&input).unwrap();
    let packet: Value = serde_json::from_slice(&input_bytes).unwrap();
    assert_eq!(packet["schema"], "n16-force-policy-long-arc-v1");
    let rows = packet["requests"].as_array().unwrap();
    let matching: Vec<_> = rows.iter().filter(|row| row["body"] == 18).collect();
    assert_eq!(matching.len(), 1);
    let request = matching[0];
    assert_eq!(request["target"], 2_000_002);
    let words = request["state_bits"].as_array().unwrap();
    assert_eq!(words.len(), 6);
    let seed: [f64; 6] = std::array::from_fn(|i| f64::from_bits(words[i].as_u64().unwrap()));
    let initial = State {
        position: seed[..3].try_into().unwrap(),
        velocity: seed[3..].try_into().unwrap(),
    };
    let epoch = request["epoch"].as_f64().unwrap();
    let endpoint = 2488070.5;
    let data = verified_force_data();
    let mut results = Vec::new();
    for maximum_step in [32.0, 1.0, 0.25] {
        let mut trajectory = ForcedTrajectory::new(
            Arc::clone(&data),
            initial,
            [epoch, 0.0],
            None,
            None,
            Some(2_000_002),
            |at| super::super::super::poles::PoleModel::IcrfTerrestrialIauSolar.at(at),
        )
        .unwrap();
        // Same initial words, equations, poles, tolerance and controller.
        // The fixed cap sequence was declared before examining any endpoint.
        trajectory.integrator = Radau::new(
            initial,
            0.0,
            Settings {
                coefficient_tolerance: 1e-9,
                initial_step: 0.001,
                maximum_step,
                maximum_attempts: 20_000_000,
                maximum_iterations: 64,
            },
        )
        .unwrap();
        println!("Pallas2100 maximum step {maximum_step}: start");
        let result = trajectory.advance(endpoint - epoch).unwrap();
        assert_eq!(trajectory.time(), endpoint - epoch);
        results.push(json!({"maximum_step_days": maximum_step,
            "state_ssb_bits": result.position.into_iter().chain(result.velocity).map(f64::to_bits).collect::<Vec<_>>()}));
        println!("Pallas2100 maximum step {maximum_step}: complete");
    }
    verified_force_data();
    assert_eq!(std::fs::read(input).unwrap(), input_bytes);
    let file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(output)
        .unwrap();
    serde_json::to_writer_pretty(
        file,
        &json!({
            "scope": "same-seed Pallas2100 maximum-step refinement; no production change",
            "automatic_acceptance": false, "epoch": epoch, "endpoint": endpoint,
            "seed_ssb_bits": request["state_bits"], "coefficient_tolerance": 1e-9,
            "rows": results,
        }),
    )
    .unwrap();
}
