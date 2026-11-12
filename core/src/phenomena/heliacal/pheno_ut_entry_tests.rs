use super::*;
use crate::constants::{JUPITER, MARS, MOON, SATURN, SUN, VENUS};

#[test]
fn planet_names_map_case_insensitively() {
    // Release `_parse_object_name` planet table (upper-cased,
    // stripped); the Sun and Moon stay admitted for the phenomena
    // entry, unlike the search entries.
    assert_eq!(parse_pheno_object_name("Venus").unwrap(), VENUS);
    assert_eq!(parse_pheno_object_name("  mars ").unwrap(), MARS);
    assert_eq!(parse_pheno_object_name("JUPITER").unwrap(), JUPITER);
    assert_eq!(parse_pheno_object_name("sAtUrN").unwrap(), SATURN);
    assert_eq!(parse_pheno_object_name("Sun").unwrap(), SUN);
    assert_eq!(parse_pheno_object_name("Moon").unwrap(), MOON);
}

#[test]
fn star_names_resolve_and_unknowns_refuse() {
    // Catalog names resolve through the shared resolver.
    let sirius = parse_pheno_object_name("Sirius").unwrap();
    assert!(is_fixed_star(sirius));
    // Numeric spellings are NOT integer ids.
    assert!(parse_pheno_object_name("3").is_err());
    assert!(parse_pheno_object_name("1000001").is_err());
    // Unknown and empty names refuse `InvalidBody`.
    let err = parse_pheno_object_name("Xyzzy").unwrap_err();
    assert!(matches!(
        err,
        Error::InvalidBody {
            operation: "heliacal_pheno_ut",
            ..
        }
    ));
    assert!(parse_pheno_object_name("").is_err());
}

#[test]
fn geopos_length_defaults_and_finite_gate() {
    // Missing positions read 0.0.
    assert_eq!(
        parse_pheno_geopos(&[]).unwrap(),
        PhenoUtSite {
            lon_deg: 0.0,
            lat_deg: 0.0,
            altitude_m: 0.0,
        }
    );
    assert_eq!(
        parse_pheno_geopos(&[12.5]).unwrap(),
        PhenoUtSite {
            lon_deg: 12.5,
            lat_deg: 0.0,
            altitude_m: 0.0,
        }
    );
    let site = parse_pheno_geopos(&[12.5, 41.9, 50.0]).unwrap();
    assert_eq!(site.lon_deg.to_bits(), 12.5_f64.to_bits());
    // Extra words are ignored.
    let site = parse_pheno_geopos(&[12.5, 41.9, 50.0, 9.0]).unwrap();
    assert_eq!(site.altitude_m.to_bits(), 50.0_f64.to_bits());
    // Non-finite words refuse fail-closed (outside the swept
    // domain; the release would propagate NaN words).
    assert!(parse_pheno_geopos(&[f64::NAN, 0.0, 0.0]).is_err());
}

#[test]
fn atmo_length_defaults_and_percent_rule() {
    // Empty tuple reads the release defaults; humidity converts
    // percent-to-fraction for the acquisition core.
    let atmo = parse_pheno_atmo(&[]);
    assert_eq!(atmo.pressure_mbar.to_bits(), 1013.25_f64.to_bits());
    assert_eq!(atmo.temperature_c.to_bits(), 15.0_f64.to_bits());
    assert_eq!(atmo.humidity_frac.to_bits(), 0.4_f64.to_bits());
    assert_eq!(atmo.met_range_km.to_bits(), 0.0_f64.to_bits());
    // Non-positive pressure (including NaN, since `NaN > 0` is
    // false) falls back to the default, like the release.
    assert_eq!(
        parse_pheno_atmo(&[0.0, 15.0, 50.0, 0.0])
            .pressure_mbar
            .to_bits(),
        1013.25_f64.to_bits()
    );
    assert_eq!(
        parse_pheno_atmo(&[f64::NAN, 15.0, 50.0, 0.0])
            .pressure_mbar
            .to_bits(),
        1013.25_f64.to_bits()
    );
    let atmo = parse_pheno_atmo(&[1013.25, 15.0, 50.0, 10.0]);
    assert_eq!(atmo.humidity_frac.to_bits(), 0.5_f64.to_bits());
    assert_eq!(atmo.met_range_km.to_bits(), 10.0_f64.to_bits());
}
