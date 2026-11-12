//! Native target admission and retained-context checks.
use super::*;

#[test]
fn moon_self_rejection_precedes_nonfinite_epoch_and_empty_registry() {
    let session = Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
    let error = session.lun_occult_where(f64::NAN, MOON, 2).unwrap_err();
    assert!(matches!(
        error.primary_error(),
        Error::InvalidEclipseSearch {
            entry: "lun_occult",
            ..
        }
    ));
    assert!(session.context().spk_registry().is_empty());
}

#[test]
fn every_negative_numeric_identifier_folds_to_sun() {
    for body in [-1, -2, -100, i32::MIN] {
        assert_eq!(numeric_occultation_target(body).unwrap(), SUN);
    }
    for body in [SUN, VENUS, MARS, JUPITER, EARTH, CERES, 999999] {
        assert_eq!(numeric_occultation_target(body).unwrap(), body);
    }
}

#[test]
fn integer_catalog_alias_requires_the_explicit_stellar_route() {
    let body = 1_000_001;
    assert!(crate::stars::canonical_star_name(body).is_some());
    assert!(matches!(
        numeric_occultation_target(body),
        Err(Error::ConfigurationMissing {
            missing_config: "stellar_lunar_occultation_route",
            ..
        })
    ));
}

#[test]
fn classical_map_and_nonclassical_model_paths_stay_distinct() {
    for body in [
        SUN, MOON, MERCURY, VENUS, MARS, JUPITER, SATURN, URANUS, NEPTUNE, PLUTO, EARTH,
    ] {
        assert!(classical_numeric_target(body));
    }
    for body in [CHIRON, PHOLUS, CERES, PALLAS, JUNO, VESTA, 40, 999999] {
        assert!(!classical_numeric_target(body));
    }
}
