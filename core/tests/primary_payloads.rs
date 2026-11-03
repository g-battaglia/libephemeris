//! Primary records preserve the native producer's bits and genuine error cause.
use libephemeris_core::{api as eph, errors, houses, native, refraction, utils};
use std::error::Error as _;

fn words(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (actual, expected) in actual.iter().zip(expected) {
        assert_eq!(actual.to_bits(), expected.to_bits());
    }
}

#[test]
fn pure_angle_and_rotation_records_preserve_ieee_words() {
    for value in [
        -720.5,
        -0.0,
        0.0,
        359.999,
        f64::INFINITY,
        f64::from_bits(0x7ff8_0000_0000_0042),
    ] {
        words(
            &[eph::degnorm(value).unwrap().value],
            &[utils::degnorm(value)],
        );
        words(
            &[eph::radnorm(value).unwrap().value],
            &[utils::radnorm(value)],
        );
        for second in [-13.25, -0.0, 0.0, 721.0] {
            words(
                &[eph::difdeg2n(value, second).unwrap().value],
                &[utils::difdeg2n(value, second)],
            );
            words(
                &[eph::difrad2n(value, second).unwrap().value],
                &[utils::difrad2n(value, second)],
            );
        }
    }
    for distance in [-0.0, 0.0, 1.234567890123] {
        let value = eph::cotrans(123.456, -45.678, distance, 23.439).unwrap();
        words(
            &value.values,
            &native::cotrans(123.456, -45.678, distance, 23.439).unwrap(),
        );
        let value =
            eph::cotrans_sp(123.456, -45.678, distance, 0.123, -0.456, -0.0, 23.439).unwrap();
        let (position, velocity) =
            native::cotrans_sp(123.456, -45.678, distance, 0.123, -0.456, -0.0, 23.439).unwrap();
        words(&value.position, &position);
        words(&value.velocity, &velocity);
    }
}

#[test]
fn position_and_house_records_retain_all_data_independent_channels() {
    let primary = eph::SessionBuilder::new().build().unwrap();
    let original = native::SessionBuilder::isolated().build().unwrap();
    for jd in [2415020.5, 2451545.0, 2460000.5] {
        for flags in [0, eph::FLG_SPEED, eph::FLG_RADIANS] {
            let value = primary.calc_ut(jd, eph::ECL_NUT, flags).unwrap();
            let expected = original.calc_ut(jd, native::ECL_NUT, flags).unwrap();
            words(&value.values, &expected.values);
            assert_eq!(value.returned_flags, expected.returned_flags);
            assert_eq!(value.diagnostic, expected.diagnostic);
            assert_eq!(value.messages, expected.messages);
            assert_eq!(
                value
                    .warnings
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>(),
                expected
                    .warnings
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
            );
        }
        for system in *b"PKWO" {
            let value = primary
                .houses_ex2(jd, eph::FLG_RADIANS, 41.9, 12.5, system)
                .unwrap();
            let expected = original
                .houses_ex2(jd, native::FLG_RADIANS, 41.9, 12.5, system)
                .unwrap();
            words(&value.cusps, &expected.cusps);
            words(&value.angles, &expected.angles);
            words(&value.cusp_speeds, &expected.cusp_speeds);
            words(&value.angle_speeds, &expected.angle_speeds);
            assert!(value.status.is_none());
            assert!(value.diagnostic.is_none());
        }
    }
}

#[test]
fn armc_preserves_twelve_or_thirty_six_cusps_without_padding() {
    let session = eph::SessionBuilder::new().build().unwrap();
    for system in *b"PG" {
        let value = session
            .houses_armc_ex2(123.456, 41.9, 23.439, system, -18.5)
            .unwrap();
        let (cusps, angles, speeds, angle_speeds) =
            houses::houses_armc_ex2(123.456, 41.9, 23.439, char::from(system), -18.5).unwrap();
        assert_eq!(value.cusps.len(), if system == b'G' { 36 } else { 12 });
        words(&value.cusps, &cusps);
        words(&value.angles, &angles);
        words(&value.cusp_speeds, &speeds);
        words(&value.angle_speeds, &angle_speeds);
    }
}

#[test]
fn refraction_records_preserve_every_detail_and_native_direction() {
    let session = eph::SessionBuilder::new().build().unwrap();
    for altitude in [-5.0, -0.0, 0.0, 0.5, 45.0, 90.0] {
        for flag in [0, 1, -7] {
            let value = session
                .refrac_extended(altitude, 1000.0, 1013.25, 10.0, 0.0065, flag)
                .unwrap();
            let (expected, (true_alt, apparent_alt, refraction, dip)) = if flag == 0 {
                refraction::refrac_extended_true_to_app(altitude, 1000.0, 1013.25, 10.0, 0.0065)
            } else {
                refraction::refrac_extended_app_to_true(altitude, 1000.0, 1013.25, 10.0, 0.0065)
            };
            words(&[value.value], &[expected]);
            words(&value.details, &[true_alt, apparent_alt, refraction, dip]);
            assert!(value.diagnostic.is_none());
        }
    }
}

#[test]
fn primary_errors_retain_the_original_structured_source() {
    let bits = 0x7ff8_0000_0000_0042;
    let error = eph::cotrans(f64::from_bits(bits), 0.0, 1.0, 23.439).unwrap_err();
    assert_eq!(error.kind(), eph::ErrorKind::InvalidInput);
    assert_eq!(error.operation(), Some("cotrans"));
    assert_eq!(error.native_scalar().unwrap().to_bits(), bits);
    assert!(error.native_status().is_none());
    let original = error
        .source()
        .unwrap()
        .downcast_ref::<errors::Error>()
        .unwrap();
    let errors::Error::NonFiniteArgument { value } = original else {
        panic!("native source lost")
    };
    assert_eq!(value.to_bits(), bits);
    assert_eq!(error.to_string(), original.to_string());
    assert_eq!(error.diagnostic(), Some(original.to_string().as_str()));
}
