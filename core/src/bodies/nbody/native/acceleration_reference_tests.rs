// SPDX-License-Identifier: AGPL-3.0
//! Explicit real-data force sampling for an independent SPK derivative check.
//!
//! This ignored diagnostic reads exact SSB input words, evaluates the complete
//! production force without advancing a trajectory, and writes native quantities
//! only. It does not modify a force setting or declare orbital acceptance.

use super::super::super::newtonian;
use super::super::*;
use super::verified_force_data;
use serde_json::{Value, json};
use std::collections::BTreeSet;
use std::path::PathBuf;

fn input_words<const N: usize>(value: &Value) -> [f64; N] {
    let words = value.as_array().expect("exact binary64 input array");
    assert_eq!(words.len(), N, "complete input vector");
    std::array::from_fn(|index| {
        let value = f64::from_bits(words[index].as_u64().expect("binary64 bits"));
        assert!(value.is_finite(), "finite diagnostic input");
        value
    })
}

fn poles(epoch: [f64; 2]) -> Result<HarmonicPoles, Error> {
    super::super::super::poles::PoleModel::IcrfTerrestrialIauSolar.at(epoch)
}

#[test]
#[ignore = "explicit analytic-reference diagnostic on authenticated DE440/N16 data"]
fn sample_complete_force_at_independent_spk_states() {
    sample_forces("n16-analytic-acceleration-v1", 70);
}

#[test]
#[ignore = "explicit record-impulse diagnostic on authenticated DE440/N16 data"]
fn sample_complete_force_at_record_quadrature_nodes() {
    // Ten targets, seven records, and the fixed 8/16/32-point quadratures.
    // This is a force query at supplied states, never an integration trial.
    sample_forces("n16-record-impulse-v1", 10 * 7 * (8 + 16 + 32));
}

fn sample_forces(schema: &str, expected_requests: usize) {
    let input = PathBuf::from(
        std::env::var("LIBEPHEMERIS_ACCELERATION_INPUT").expect("explicit input packet path"),
    );
    let output = PathBuf::from(
        std::env::var("LIBEPHEMERIS_ACCELERATION_OUTPUT").expect("explicit output capture path"),
    );
    assert!(!output.exists(), "never overwrite an earlier diagnostic");
    let packet: Value = serde_json::from_reader(std::fs::File::open(input).unwrap()).unwrap();
    assert_eq!(packet["schema"], schema);
    let requests = packet["requests"].as_array().unwrap();
    assert_eq!(
        requests.len(),
        expected_requests,
        "complete diagnostic inventory"
    );
    let data = verified_force_data();
    let mut identifiers = BTreeSet::new();
    let mut results = Vec::with_capacity(requests.len());
    for request in requests {
        let id = request["id"].as_str().unwrap();
        assert!(identifiers.insert(id.to_owned()), "unique request identity");
        let body = request["body"].as_i64().unwrap();
        let mpc = match body {
            17 => 1,
            18 => 2,
            19 => 3,
            20 => 4,
            10010 => 10,
            10016 => 16,
            10052 => 52,
            10087 => 87,
            10511 => 511,
            10704 => 704,
            _ => panic!("not an original curated N16 target"),
        };
        let target = 2_000_000 + mpc;
        assert_eq!(request["target"], target);
        let epoch = input_words::<2>(&request["epoch_bits"]);
        let words = input_words::<6>(&request["state_bits"]);
        let initial = State {
            position: words[..3].try_into().unwrap(),
            velocity: words[3..].try_into().unwrap(),
        };
        let mut trajectory = ForcedTrajectory::new(
            Arc::clone(&data),
            initial,
            epoch,
            None,
            None,
            Some(target),
            poles,
        )
        .unwrap();
        let acceleration = trajectory.acceleration().unwrap();
        assert_eq!(trajectory.time(), 0.0, "force query must not advance time");
        for (before, after) in words.into_iter().zip(
            trajectory
                .state()
                .position
                .into_iter()
                .chain(trajectory.state().velocity),
        ) {
            assert_eq!(
                before.to_bits(),
                after.to_bits(),
                "force query preserves PV"
            );
        }
        let mut sources = Vec::new();
        data.frame(epoch, &mut sources).unwrap();
        let sun = sources
            .iter()
            .find(|source| source.point.body == 10)
            .unwrap();
        let stored_target = sources
            .iter()
            .find(|source| source.point.body == target)
            .unwrap();
        // The complete force above is always active. This separately reported
        // Newtonian subtotal is explanatory, never a replacement trajectory.
        let newton = newtonian::acceleration_from(
            initial.position,
            sources.iter().map(|source| source.point),
            Some(target),
        )
        .unwrap();
        results.push(json!({
            "request": request,
            "acceleration_au_day2": acceleration,
            "newtonian_au_day2": newton,
            "sun_position_au": sun.point.position,
            "sun_velocity_au_day": sun.velocity,
            "sun_acceleration_au_day2": sun.post_newtonian_acceleration.unwrap(),
            "stored_target_position_au": stored_target.point.position,
            "stored_target_velocity_au_day": stored_target.velocity,
            "source_count": sources.len(),
            "excluded_physical_target": target,
        }));
    }
    let capture = json!({
        "schema": schema,
        "pole_model": "IcrfTerrestrialIauSolar",
        "units": "SSB ICRF AU, AU/day, AU/day2; TDB days",
        "astronomical_unit_km": data.catalogue.astronomical_unit_km,
        "requests": results,
    });
    let file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(output)
        .unwrap();
    serde_json::to_writer_pretty(file, &capture).unwrap();
}
