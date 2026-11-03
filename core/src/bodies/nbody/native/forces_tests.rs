// SPDX-License-Identifier: AGPL-3.0
//! Prepared force-identity checks, not an executed qualification report.

use super::*;

fn sun() -> Source {
    Source {
        point: PointMass {
            body: 10,
            position: [0.0; 3],
            gravitational_parameter: 1.0,
        },
        velocity: [0.0; 3],
        post_newtonian_acceleration: Some([0.0; 3]),
    }
}
fn empty_model() -> ForceModel {
    ForceModel {
        relativity: None,
        zonals: Vec::new(),
        nongravitational: None,
    }
}
fn close(actual: [f64; 3], expected: [f64; 3]) {
    let scale = length(expected).unwrap();
    for axis in 0..3 {
        if scale == 0.0 {
            assert_eq!(actual[axis], expected[axis]);
        } else {
            assert!(
                (actual[axis] - expected[axis]).abs() <= 512.0 * f64::EPSILON * scale,
                "axis {axis}: {actual:?} != {expected:?}"
            );
        }
    }
}

#[test]
fn central_ppn_reduces_to_the_schwarzschild_test_particle_equation() {
    let source = sun();
    let state = State {
        position: [2.0, 0.3, -0.2],
        velocity: [0.01, 0.009, 0.003],
    };
    let c = 173.0;
    let radius = length(state.position).unwrap();
    let mu = source.point.gravitational_parameter;
    let speed_squared = dot(state.velocity, state.velocity);
    let radial_speed = dot(state.position, state.velocity);
    let factor = mu / (c * c * radius * radius * radius);
    let expected = std::array::from_fn(|axis| {
        factor
            * ((4.0 * mu / radius - speed_squared) * state.position[axis]
                + 4.0 * radial_speed * state.velocity[axis])
    });
    for outer in [RelativisticSources::Sun, RelativisticSources::MajorBodies] {
        let actual = Eih {
            speed_of_light: c,
            outer,
        }
        .acceleration(state, &[source], None)
        .unwrap();
        close(actual, expected);
    }
}

#[test]
fn ppn_potentials_outer_selection_and_source_acceleration_stay_distinct() {
    // Two prescribed stationary masses: Ui=1+2=3; U_sun=2/2=1;
    // U_other=1/2. Hence the Sun contributes +13/c^2 and the other
    // source -25/c^2. The asteroid must not enter these PPN sums.
    let mut primary = sun();
    let other = Source {
        point: PointMass {
            body: 5,
            position: [2.0, 0.0, 0.0],
            gravitational_parameter: 2.0,
        },
        ..sun()
    };
    let asteroid = Source {
        point: PointMass {
            body: 2_000_001,
            position: [3.0, 0.0, 0.0],
            gravitational_parameter: 0.4,
        },
        post_newtonian_acceleration: None,
        ..sun()
    };
    let state = State {
        position: [1.0, 0.0, 0.0],
        velocity: [0.0; 3],
    };
    let c = 173.0;
    let solar = Eih {
        speed_of_light: c,
        outer: RelativisticSources::Sun,
    };
    let major = Eih {
        outer: RelativisticSources::MajorBodies,
        ..solar
    };
    close(
        solar
            .acceleration(state, &[primary, other, asteroid], None)
            .unwrap(),
        [13.0 / (c * c), 0.0, 0.0],
    );
    close(
        major
            .acceleration(state, &[primary, other, asteroid], None)
            .unwrap(),
        [-12.0 / (c * c), 0.0, 0.0],
    );
    primary.velocity = [0.2, 0.0, 0.0];
    primary.post_newtonian_acceleration = Some([0.01, 0.0, 0.03]);
    // Radial 12.985 + velocity0.12 + source acceleration0.035;
    // normal acceleration is (7/2)*0.03. These literal contractions
    // exercise the terms lost by a heliocentric-only GR shortcut.
    close(
        solar
            .acceleration(state, &[primary, other, asteroid], None)
            .unwrap(),
        [13.14 / (c * c), 0.0, 0.105 / (c * c)],
    );
}

#[test]
fn zonal_accelerations_match_equatorial_and_polar_potential_gradients() {
    let source = sun();
    let radius = 2.0_f64;
    let reference = 0.1_f64;
    let coefficients = [0.02, -0.003, 0.0004];
    let zonals = Zonals::new(10, reference, [0.0, 0.0, 4.0], coefficients).unwrap();
    let q = reference / radius;
    let scale = source.point.gravitational_parameter / (radius * radius);
    close(
        zonals.acceleration([radius, 0.0, 0.0], source).unwrap(),
        [
            scale * (-1.5 * coefficients[0] * q.powi(2) + 15.0 / 8.0 * coefficients[2] * q.powi(4)),
            0.0,
            scale * 1.5 * coefficients[1] * q.powi(3),
        ],
    );
    close(
        zonals.acceleration([0.0, 0.0, radius], source).unwrap(),
        [
            0.0,
            0.0,
            scale
                * (3.0 * coefficients[0] * q.powi(2)
                    + 4.0 * coefficients[1] * q.powi(3)
                    + 5.0 * coefficients[2] * q.powi(4)),
        ],
    );
}

#[test]
fn zonal_parity_follows_its_legendre_degree() {
    let source = sun();
    let position = [1.2, -0.7, 0.4];
    for degree in 2..=4 {
        let mut coefficients = [0.0; 3];
        coefficients[degree - 2] = 0.01;
        let zonals = Zonals::new(10, 0.1, [0.1, -0.3, 0.9], coefficients).unwrap();
        let positive = zonals.acceleration(position, source).unwrap();
        let negative = zonals
            .acceleration(position.map(|value| -value), source)
            .unwrap();
        close(
            negative,
            positive.map(|value| if degree % 2 == 0 { -value } else { value }),
        );
    }
}

#[test]
fn radial_transverse_normal_parameters_retain_their_physical_orientation() {
    let state = State {
        position: [2.0, 0.0, 0.0],
        velocity: [0.4, 0.5, 0.0],
    };
    let parameters = [1e-8, -2e-8, 3e-8];
    let model = Nongravitational {
        parameters,
        law: DistanceLaw::InverseSquare,
    };
    close(
        model.acceleration(state, sun()).unwrap(),
        parameters.map(|value| value / 4.0),
    );
    let zero_plane = State {
        velocity: [0.4, 0.0, 0.0],
        ..state
    };
    assert!(matches!(
        model.acceleration(zero_plane, sun()),
        Err(NumericalError::InvalidOrbit(_))
    ));
    let radial = Nongravitational {
        parameters: [1e-8, 0.0, 0.0],
        ..model
    };
    close(
        radial.acceleration(zero_plane, sun()).unwrap(),
        [2.5e-9, 0.0, 0.0],
    );
    let zero = Nongravitational {
        parameters: [0.0; 3],
        ..model
    };
    assert_eq!(zero.acceleration(zero_plane, sun()).unwrap(), [0.0; 3]);
}

#[test]
fn logarithmic_marsden_law_is_normalized_without_an_overflowing_power() {
    // Toy exponents make the original expression independently elementary.
    let law = DistanceLaw::Marsden {
        alpha: 2.0,
        scale: 1.0,
        m: 2.0,
        n: 1.0,
        k: 1.0,
    };
    for radius in [0.01_f64, 0.5, 1.0, 2.0, 100.0] {
        let actual = law.value(radius).unwrap();
        let expected = 2.0 / (radius * radius * (1.0 + radius));
        assert!((actual - expected).abs() <= 64.0 * f64::EPSILON * expected);
    }
    assert_eq!(law.value(1e200).unwrap(), 0.0);
}

#[test]
fn self_exclusion_retains_every_other_newtonian_perturber() {
    let mut other = sun();
    other.point = PointMass {
        body: 5,
        position: [2.0, 0.0, 0.0],
        gravitational_parameter: 0.2,
    };
    let mut target = other;
    target.point = PointMass {
        body: 2_000_001,
        position: [1.0, 0.0, 0.0],
        gravitational_parameter: 0.001,
    };
    target.post_newtonian_acceleration = None;
    let sources = [sun(), other, target];
    let state = State {
        position: target.point.position,
        velocity: [0.0, 1.0, 0.0],
    };
    assert!(matches!(
        empty_model().acceleration(state, &sources, None),
        Err(NumericalError::SingularForce { body: 2_000_001 })
    ));
    close(
        empty_model()
            .acceleration(state, &sources, Some(2_000_001))
            .unwrap(),
        [-0.8, 0.0, 0.0],
    );
}

#[test]
fn full_force_composition_keeps_newton_harmonics_relativity_and_nongravitation() {
    let state = State {
        position: [2.0, 0.3, -0.1],
        velocity: [0.01, 0.02, 0.005],
    };
    let source = sun();
    let eih = Eih {
        speed_of_light: 173.0,
        outer: RelativisticSources::MajorBodies,
    };
    let zonals = Zonals::new(10, 0.1, [0.0, 0.0, 1.0], [0.001, -0.0001, 0.00001]).unwrap();
    let ng = Nongravitational {
        parameters: [1e-8, 2e-8, -3e-8],
        law: DistanceLaw::InverseSquare,
    };
    let model = ForceModel {
        relativity: Some(eih),
        zonals: vec![zonals],
        nongravitational: Some(ng),
    };
    let mut expected = newtonian::acceleration(state.position, &[source.point], None)
        .unwrap()
        .map(Sum::new);
    add(
        &mut expected,
        eih.acceleration(state, &[source], None).unwrap(),
    );
    add(
        &mut expected,
        zonals.acceleration(state.position, source).unwrap(),
    );
    add(&mut expected, ng.acceleration(state, source).unwrap());
    close(
        model.acceleration(state, &[source], None).unwrap(),
        expected.map(Sum::value),
    );
}

#[test]
fn missing_nonfinite_and_duplicate_sources_are_never_a_force_fallback() {
    let state = State {
        position: [1.0, 0.0, 0.0],
        velocity: [0.0, 1.0, 0.0],
    };
    let eih = Eih {
        speed_of_light: 173.0,
        outer: RelativisticSources::Sun,
    };
    let mut model = empty_model();
    model.relativity = Some(eih);
    assert_eq!(
        model.acceleration(state, &[], None),
        Err(NumericalError::MissingSource { body: 10 })
    );
    assert_eq!(
        model.acceleration(state, &[sun(), sun()], None),
        Err(NumericalError::DuplicateSource { body: 10 })
    );
    let mut invalid = sun();
    invalid.velocity[0] = f64::NAN;
    assert_eq!(
        model.acceleration(state, &[invalid], None),
        Err(NumericalError::NonFinite)
    );
    for c in [0.0, -1.0] {
        model.relativity = Some(Eih {
            speed_of_light: c,
            ..eih
        });
        assert_eq!(
            model.acceleration(state, &[sun()], None),
            Err(NumericalError::InvalidSettings)
        );
    }
    assert!(matches!(
        Zonals::new(10, 0.1, [0.0; 3], [0.001, 0.0, 0.0]),
        Err(NumericalError::InvalidSettings)
    ));
}
