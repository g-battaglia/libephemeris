// SPDX-License-Identifier: AGPL-3.0
//! Explicit bounded diagnosis on authenticated real force files, never CI data.
use super::*;
use sha2::{Digest, Sha256};
use std::io::Read;

#[path = "trajectory_horizons_tests.rs"]
mod horizons_tests;

#[path = "acceleration_reference_tests.rs"]
mod acceleration_reference_tests;

#[path = "force_policy_reference_tests.rs"]
mod force_policy_reference_tests;

#[path = "pallas_step_reference_tests.rs"]
mod pallas_step_reference_tests;

#[path = "ivp_reference_tests.rs"]
mod ivp_reference_tests;

fn verified_force_data() -> Arc<Ephemeris> {
    verified_force_data_for("linux_p1550p2650.440")
}

fn verified_force_data_for(planet_name: &str) -> Arc<Ephemeris> {
    let planet_hash = match planet_name {
        "linux_p1550p2650.440" => {
            "29915576d0a6555766b99485ac3056ee415e86df4fce282611c31afb329ad062"
        }
        "linux_m13000p17000.441" => {
            "476096486def4e41bfceb29aa27f50784da0bce318902bcf7b88caad058cd4da"
        }
        _ => panic!("explicit diagnostic requires a declared authenticated planet image"),
    };
    let root = std::path::PathBuf::from(
        std::env::var("LIBEPHEMERIS_TEST_DATA_DIR")
            .expect("set LIBEPHEMERIS_TEST_DATA_DIR for this explicit data diagnostic"),
    )
    .join("assist");
    for (name, expected) in [
        (planet_name, planet_hash),
        (
            "sb441-n16.bsp",
            "919d612ce3c72a78fc7158f9120156542d0f21e6b8b052e4c1339c759747fd90",
        ),
    ] {
        let mut file = std::fs::File::open(root.join(name)).unwrap();
        let mut digest = Sha256::new();
        let mut buffer = [0_u8; 65536];
        loop {
            let size = file.read(&mut buffer).unwrap();
            if size == 0 {
                break;
            }
            digest.update(&buffer[..size]);
        }
        let actual: String = digest
            .finalize()
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect();
        assert_eq!(actual, expected, "selected data identity");
    }
    Arc::new(Ephemeris::open(&root.join(planet_name), Some(&root.join("sb441-n16.bsp"))).unwrap())
}

#[test]
#[ignore = "explicit real-data diagnostic; requires authenticated DE440/N16 assets"]
fn near_terrestrial_forward_endpoint_has_bounded_force_work() {
    let data = verified_force_data();
    let elements = crate::minor_bodies::OrbitalElements::new(
        "near-terrestrial".into(),
        2451545.0,
        1.01,
        0.02,
        0.001,
        102.9,
        0.0,
        358.0,
        0.97,
    )
    .unwrap();
    // Reuse the actual supplied-element/public frame initialization. Nothing
    // is inferred from an already-propagated candidate or legacy endpoint.
    let heliocentric = crate::nbody::assist_element_state(&elements, 2451545.0).unwrap();
    let equatorial = crate::nbody::rotated(heliocentric, true);
    let sun = data.particle(0, [2451545.0, 0.0]).unwrap();
    let initial = State {
        position: std::array::from_fn(|axis| equatorial.position[axis] + sun.point.position[axis]),
        velocity: std::array::from_fn(|axis| equatorial.velocity[axis] + sun.velocity[axis]),
    };
    fn actual_poles(epoch: [f64; 2]) -> Result<HarmonicPoles, Error> {
        super::super::poles::PoleModel::IcrfTerrestrialIauSolar.at(epoch)
    }
    let mut trajectory = ForcedTrajectory::new(
        Arc::clone(&data),
        initial,
        [2451545.0, 0.0],
        None,
        None,
        None,
        actual_poles,
    )
    .unwrap();
    let mut candidate = trajectory.integrator.clone();
    let origin = trajectory.origin;
    let mut calls = 0_usize;
    // The cutoff bounds an investigation, not an accuracy tolerance or a new
    // public work/error contract. Exceeding it fails this diagnostic explicitly.
    let result = candidate.advance_with_error(30.0, &mut |relative, state| {
        calls += 1;
        if calls > 100_000 {
            return Err(ForceError::Numerical(NumericalError::StepLimit));
        }
        if calls % 10_000 == 0 {
            eprintln!("force_calls={calls} stage_relative_days={relative:.17e}");
        }
        let epoch = origin.add_offset(relative).map_err(ForceError::Frame)?;
        let epoch = [epoch.hi, epoch.lo];
        let orientation = actual_poles(epoch).map_err(ForceError::Frame)?;
        trajectory
            .evaluator
            .acceleration(epoch, *state, orientation)
    });
    let time = candidate.time();
    let earth = data.particle(3, [origin.hi, origin.lo + time]).unwrap();
    let separation = std::array::from_fn::<_, 3, _>(|axis| {
        candidate.state().position[axis] - earth.point.position[axis]
    });
    let distance = separation[0].hypot(separation[1]).hypot(separation[2]);
    eprintln!(
        "result={result:?} accepted_relative_days={time:.17e} force_calls={calls} earth_distance_au={distance:.17e}"
    );
    assert!(
        result.is_ok(),
        "bounded native forward-force diagnostic did not finish"
    );
    assert_eq!(candidate.time(), 30.0);
}

#[test]
#[ignore = "explicit force-family diagnostic on exact public-control input words"]
fn force_families_at_shared_control_state() {
    use super::super::forces::ForceModel;
    use super::super::poles::PoleModel;
    let input =
        std::env::var("LIBEPHEMERIS_FORCE_PROBE_INPUT").expect("provide diagnostic input JSON");
    let control: serde_json::Value =
        serde_json::from_reader(std::fs::File::open(input).unwrap()).unwrap();
    fn words<const N: usize>(value: &serde_json::Value) -> [f64; N] {
        let words = value.as_array().unwrap();
        assert_eq!(words.len(), N);
        std::array::from_fn(|axis| f64::from_bits(words[axis].as_u64().unwrap()))
    }
    let epoch = words::<2>(&control["epoch_bits"]);
    let state_words = words::<6>(&control["state_bits"]);
    let state = State {
        position: state_words[..3].try_into().unwrap(),
        velocity: state_words[3..].try_into().unwrap(),
    };
    let data = verified_force_data();
    let mut sources = Vec::new();
    data.frame(epoch, &mut sources).unwrap();
    let newtonian = ForceModel {
        relativity: None,
        zonals: Vec::new(),
        nongravitational: None,
    };
    let all_newtonian = newtonian.acceleration(state, &sources, None).unwrap();
    // Exact input words/selected files are shared, not an inferred pole or an
    // endpoint-fitted perturbation. Both pre-existing independently defined
    // poles are evaluated only as scientific diagnostic models here.
    for (label, pole_model) in [
        ("PublishedDateIers", PoleModel::PublishedDateIers),
        (
            "IcrfTerrestrialIauSolar",
            PoleModel::IcrfTerrestrialIauSolar,
        ),
    ] {
        let poles = pole_model.at(epoch).unwrap();
        for family in ["EARTH_HARMONICS", "SUN_HARMONICS", "GR_EIH", "TOTAL"] {
            let mut model = data.catalogue.force_model(poles, None).unwrap();
            match family {
                "EARTH_HARMONICS" => {
                    model.relativity = None;
                    model.zonals = vec![model.zonals[1]];
                }
                "SUN_HARMONICS" => {
                    model.relativity = None;
                    model.zonals = vec![model.zonals[0]];
                }
                "GR_EIH" => {
                    model.zonals.clear();
                }
                _ => {}
            }
            let total = model.acceleration(state, &sources, None).unwrap();
            let perturbation =
                std::array::from_fn::<_, 3, _>(|axis| total[axis] - all_newtonian[axis]);
            eprintln!(
                "{}",
                serde_json::json!({"pole_model":label,"family":family,"acceleration":perturbation,"earth_pole":poles.earth})
            );
        }
    }
    for family in ["SUN", "PLANETS", "ASTEROIDS"] {
        let selected: Vec<_> = sources
            .iter()
            .copied()
            .filter(|source| match family {
                "SUN" => source.point.body == 10,
                "PLANETS" => {
                    source.point.body != 10 && source.post_newtonian_acceleration.is_some()
                }
                _ => source.post_newtonian_acceleration.is_none(),
            })
            .collect();
        let acceleration = newtonian.acceleration(state, &selected, None).unwrap();
        eprintln!(
            "{}",
            serde_json::json!({"family":family,"acceleration":acceleration})
        );
    }
}
