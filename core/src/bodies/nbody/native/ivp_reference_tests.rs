// SPDX-License-Identifier: AGPL-3.0
//! Two original IVP paths, captured without another numerical/force route.
use super::super::super::radau_tables::NODES;
use super::super::*;
use super::verified_force_data;
use serde_json::{Value, json};

fn bits(state: State) -> [u64; 6] {
    std::array::from_fn(|i| {
        if i < 3 {
            state.position[i]
        } else {
            state.velocity[i - 3]
        }
        .to_bits()
    })
}

fn poles(epoch: [f64; 2]) -> Result<HarmonicPoles, Error> {
    super::super::super::poles::PoleModel::IcrfTerrestrialIauSolar.at(epoch)
}

#[test]
#[ignore = "explicit accepted-polynomial IVP reference input; authenticated DE440/N16 required"]
fn capture_original_2024_accepted_polynomials() {
    let input = std::env::var("LIBEPHEMERIS_IVP_INPUT").unwrap();
    let output = std::path::PathBuf::from(std::env::var("LIBEPHEMERIS_IVP_OUTPUT").unwrap());
    assert!(!output.exists(), "never overwrite a previous path capture");
    let packet: Value = serde_json::from_reader(std::fs::File::open(input).unwrap()).unwrap();
    assert_eq!(packet["schema"], "n16-force-policy-long-arc-v1");
    let requests = packet["requests"].as_array().unwrap();
    assert_eq!(requests.len(), 7);
    let selected = &requests[5..];
    let data = verified_force_data();
    let orientation = poles([2461000.5, 0.0]).unwrap();
    let normalized = |v: [f64; 3]| {
        let norm = v[0].hypot(v[1]).hypot(v[2]);
        v.map(|x| (x / norm).to_bits())
    };
    let mut records = Vec::new();
    for (request, (body, target)) in selected
        .iter()
        .zip([(10016, 2_000_016), (10704, 2_000_704)])
    {
        assert_eq!(request["body"].as_i64(), Some(body));
        assert_eq!(request["epoch"].as_f64(), Some(2461000.5));
        assert_eq!(request["endpoint"].as_f64(), Some(2460310.5));
        assert_eq!(request["target"].as_i64(), Some(target));
        let seed = request["state_bits"].as_array().unwrap();
        assert_eq!(seed.len(), 6);
        let values: [f64; 6] = std::array::from_fn(|i| f64::from_bits(seed[i].as_u64().unwrap()));
        let initial = State {
            position: values[..3].try_into().unwrap(),
            velocity: values[3..].try_into().unwrap(),
        }
        .checked()
        .unwrap();
        let exclude = i32::try_from(request["target"].as_i64().unwrap()).unwrap();
        let mut ordinary = ForcedTrajectory::new(
            Arc::clone(&data),
            initial,
            [2461000.5, 0.0],
            None,
            None,
            Some(exclude),
            poles,
        )
        .unwrap();
        let mut captured = ForcedTrajectory::new(
            Arc::clone(&data),
            initial,
            [2461000.5, 0.0],
            None,
            None,
            Some(exclude),
            poles,
        )
        .unwrap();
        captured.integrator.enable_step_capture();
        let ordinary_state = ordinary.advance(-690.0).unwrap();
        let captured_state = captured.advance(-690.0).unwrap();
        assert_eq!(bits(ordinary_state), bits(captured_state));
        let (steps, complete) = captured.integrator.step_capture().unwrap();
        assert!(
            complete && !steps.is_empty(),
            "complete accepted path required"
        );
        assert_eq!(steps.first().unwrap().start_time, 0.0);
        assert_eq!(bits(steps.first().unwrap().initial), bits(initial));
        assert_eq!(steps.last().unwrap().end_time, -690.0);
        assert_eq!(bits(steps.last().unwrap().endpoint), bits(captured_state));
        for pair in steps.windows(2) {
            assert_eq!(pair[0].end_time, pair[1].start_time);
            assert_eq!(bits(pair[0].endpoint), bits(pair[1].initial));
        }
        let path: Vec<_> = steps.iter().map(|step| json!({
            "start_time_bits": step.start_time.to_bits(),
            "end_time_bits": step.end_time.to_bits(),
            "step_bits": step.step.to_bits(),
            "initial_state_bits": bits(step.initial),
            "endpoint_state_bits": bits(step.endpoint),
            "acceleration_newton_bits": step.acceleration_coefficients.map(|row| row.map(f64::to_bits)),
        })).collect();
        records.push(json!({
            "request": request, "endpoint_state_bits": bits(captured_state),
            "untraced_endpoint_state_bits": bits(ordinary_state),
            "capture_complete": complete, "steps": path,
        }));
    }
    verified_force_data();
    let file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(output)
        .unwrap();
    serde_json::to_writer_pretty(file, &json!({
        "schema": "n16-accepted-polynomials-v1",
        "scope": "two original exact-seed2024 paths; inputs to a future uniform IVP bound, not trajectory acceptance",
        "automatic_acceptance": false, "global_error_bound_available": false,
        "node_bits": NODES.map(f64::to_bits), "records": records,
        "normalized_pole_bits": {"sun": normalized(orientation.sun), "earth": normalized(orientation.earth)},
    })).unwrap();
}
