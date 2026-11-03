// SPDX-License-Identifier: AGPL-3.0
//! Raw-state initialization against an independent official Horizons integration.
//! This explicit diagnostic does not infer a legacy pole from returned forces
//! or declare the JPL and selected force approximations identical.
use super::super::super::poles::PoleModel;
use super::*;

#[test]
#[ignore = "explicit independent JPL reference and authenticated DE441/N16 required"]
fn identical_horizons_ssb_seed_compares_independent_pole_models() {
    let path = std::env::var("LIBEPHEMERIS_HORIZONS_INPUT")
        .expect("provide independent JPL reference JSON");
    let reference: serde_json::Value =
        serde_json::from_reader(std::fs::File::open(path).unwrap()).unwrap();
    fn words(value: &serde_json::Value) -> State {
        let words = value.as_array().unwrap();
        assert_eq!(words.len(), 6);
        let state: [f64; 6] =
            std::array::from_fn(|axis| f64::from_bits(words[axis].as_u64().unwrap()));
        State {
            position: state[..3].try_into().unwrap(),
            velocity: state[3..].try_into().unwrap(),
        }
    }
    let rows = reference["rows"].as_array().unwrap();
    assert_eq!(rows.len(), 2);
    let start = rows[0]["epoch"].as_f64().unwrap();
    let end = rows[1]["epoch"].as_f64().unwrap();
    let initial = words(&rows[0]["state_bits"]);
    let expected = words(&rows[1]["state_bits"]);
    let data = verified_force_data_for("linux_m13000p17000.441");
    let policies: [(&str, PoleProvider); 2] = [
        ("PublishedDateIers", |epoch| {
            PoleModel::PublishedDateIers.at(epoch)
        }),
        ("IcrfTerrestrialIauSolar", |epoch| {
            PoleModel::IcrfTerrestrialIauSolar.at(epoch)
        }),
    ];
    for (label, poles) in policies {
        let mut trajectory = ForcedTrajectory::new(
            Arc::clone(&data),
            initial,
            [start, 0.0],
            None,
            None,
            None,
            poles,
        )
        .unwrap();
        let actual = trajectory.advance(end - start).unwrap();
        let position_error: [f64; 3] =
            std::array::from_fn(|axis| actual.position[axis] - expected.position[axis]);
        let velocity_error: [f64; 3] =
            std::array::from_fn(|axis| actual.velocity[axis] - expected.velocity[axis]);
        let position_norm = position_error[0]
            .hypot(position_error[1])
            .hypot(position_error[2]);
        let velocity_norm = velocity_error[0]
            .hypot(velocity_error[1])
            .hypot(velocity_error[2]);
        // The same fixed output-plane rotation is applied to each error;
        // a common selected Sun subtraction cancels exactly in the difference.
        let vertical_position = crate::nbody::equatorial_to_ecliptic(position_error)[2];
        let vertical_velocity = crate::nbody::equatorial_to_ecliptic(velocity_error)[2];
        eprintln!(
            "{}",
            serde_json::json!({"pole_model":label,
            "initial_bits":rows[0]["state_bits"],"epoch":start,"end":end,
            "native_state":[actual.position,actual.velocity],"reference_state":[expected.position,expected.velocity],
            "position_error":position_error,"velocity_error":velocity_error,
            "position_error_au":position_norm,"velocity_error_au_day":velocity_norm,
            "ecliptic_vertical_error":[vertical_position,vertical_velocity],
            "source_url":reference["source_url"]})
        );
    }
}
