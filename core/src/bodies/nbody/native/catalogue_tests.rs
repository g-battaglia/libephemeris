// SPDX-License-Identifier: AGPL-3.0
//! Prepared selected-parameter tests; no astronomical measurements are fixtures.

use super::*;

fn metadata(name: &str) -> Result<f64, DeError> {
    match name {
        "GMS" => Ok(1.0),
        "GMB" => Ok(8.0),
        "GM1" | "GM2" | "GM4" | "GM5" | "GM6" | "GM7" | "GM8" | "GM9" => Ok(0.5),
        "CLIGHT" => Ok(20.0),
        "ASUN" => Ok(4.0),
        "RE" => Ok(1.0),
        "J2SUN" => Ok(0.001),
        "J2E" => Ok(0.01),
        "J3E" => Ok(-0.002),
        "J4E" => Ok(-0.003),
        _ => Err(DeError::MissingConstant(name.into())),
    }
}
fn poles() -> HarmonicPoles {
    HarmonicPoles {
        sun: [0.0, 0.0, 1.0],
        earth: [0.0, 0.0, 1.0],
    }
}

#[test]
fn selected_metadata_preserves_mass_splitting_and_declared_physical_units() {
    let values = Catalogue::read(2.0, 3.0, metadata).unwrap();
    assert_eq!(MAJOR_TARGETS, [10, 1, 2, 399, 301, 4, 5, 6, 7, 8, 9]);
    assert_eq!(
        values.major_masses,
        [1.0, 0.5, 0.5, 6.0, 2.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
    );
    assert_eq!(values.major_masses[3] + values.major_masses[4], 8.0);
    assert_eq!(values.astronomical_unit_km, 2.0);
    assert_eq!(values.speed_of_light_au_day, 864000.0);
    assert_eq!(values.sun_radius, 2.0);
    assert_eq!(values.earth_radius, 0.5);
    assert_eq!(values.sun_j2, 0.001);
    assert_eq!(values.earth_j, [0.01, -0.002, -0.003]);
}

#[test]
fn observed_public_selection_keeps_all_families_and_the_inverse_square_law() {
    let catalogue = Catalogue::read(2.0, 3.0, metadata).unwrap();
    let model = catalogue
        .force_model(poles(), Some([1e-8, -2e-8, 3e-8]))
        .unwrap();
    assert!(matches!(
        model.relativity,
        Some(Eih {
            outer: RelativisticSources::Sun,
            ..
        })
    ));
    assert_eq!(model.zonals.len(), 2);
    let nongravitational = model.nongravitational.unwrap();
    assert_eq!(nongravitational.parameters, [1e-8, -2e-8, 3e-8]);
    assert!(matches!(nongravitational.law, DistanceLaw::InverseSquare));
    assert!(
        catalogue
            .force_model(poles(), None)
            .unwrap()
            .nongravitational
            .is_none()
    );
}

#[test]
fn missing_or_invalid_metadata_is_not_replaced_by_a_nominal_constant() {
    let absent = Catalogue::read(2.0, 3.0, |name| {
        if name == "J3E" {
            Err(DeError::MissingConstant(name.into()))
        } else {
            metadata(name)
        }
    });
    assert!(matches!(absent,Err(DeError::MissingConstant(name)) if name == "J3E"));
    for field in ["CLIGHT", "RE", "ASUN", "GM2"] {
        assert!(matches!(
            Catalogue::read(2.0, 3.0, |name| if name == field {
                Ok(-1.0)
            } else {
                metadata(name)
            }),
            Err(DeError::Malformed(_))
        ));
    }
    for bad in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(matches!(
            Catalogue::read(bad, 3.0, metadata),
            Err(DeError::Malformed(_))
        ));
        assert!(matches!(
            Catalogue::read(2.0, bad, metadata),
            Err(DeError::Malformed(_))
        ));
    }
}

#[test]
fn force_assembly_requires_both_explicit_valid_poles() {
    let catalogue = Catalogue::read(2.0, 3.0, metadata).unwrap();
    for invalid in [
        HarmonicPoles {
            sun: [0.0; 3],
            ..poles()
        },
        HarmonicPoles {
            earth: [0.0; 3],
            ..poles()
        },
    ] {
        assert!(matches!(
            catalogue.force_model(invalid, None),
            Err(NumericalError::InvalidSettings)
        ));
    }
    assert!(matches!(
        catalogue.force_model(poles(), Some([f64::NAN, 0.0, 0.0])),
        Err(NumericalError::NonFinite)
    ));
}

#[test]
fn actual_reader_with_missing_required_metadata_fails_instead_of_fabricating_a_model() {
    let fixture = super::super::classic_de::tests::Fixture::little();
    let (image, _temporary) = fixture.open();
    assert!(matches!(
        Catalogue::load(&image),
        Err(DeError::MissingConstant(_))
    ));
}
