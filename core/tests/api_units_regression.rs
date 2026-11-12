//! Same-engine flag/unit regressions; no external reference vectors.
use libephemeris_core::native as e;

/// Real LEB data is mandatory for this integration gate, never a passing skip.
fn session() -> Result<e::Session, e::Error> {
    let path = std::env::var_os("KERYKEION_LEB_FILE")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| {
            std::path::PathBuf::from(
                std::env::var_os("HOME").expect("HOME or explicit data required"),
            )
            .join(".libephemeris/leb/base_core.leb2")
        });
    let mut builder = e::SessionBuilder::isolated();
    builder.leb_file(&path)?.sid_mode(e::SIDM_LAHIRI, 0.0, 0.0);
    builder.build()
}

#[test]
fn speed3_keeps_sidereal_drift_without_changing_public_flag_echo() -> Result<(), e::Error> {
    let session = session()?;
    for jd in [2448087.8541666665, 2451545.0, 2460000.5] {
        for body in [e::SUN, e::MARS] {
            for is_ut in [false, true] {
                let calc = |flags| {
                    if is_ut {
                        session.calc_ut(jd, body, flags)
                    } else {
                        session.calc(jd, body, flags)
                    }
                };
                let drift = |speed| -> Result<f64, e::Error> {
                    let flags = e::FLG_SWIEPH | speed;
                    let tropical = calc(flags)?;
                    let sidereal = calc(flags | e::FLG_SIDEREAL)?;
                    assert_eq!(
                        sidereal.returned_flags & (e::FLG_SPEED | e::FLG_SPEED3),
                        speed
                    );
                    assert_eq!(sidereal.values[4..], tropical.values[4..]);
                    Ok(tropical.values[3] - sidereal.values[3])
                };
                assert_eq!(
                    drift(e::FLG_SPEED)?.to_bits(),
                    drift(e::FLG_SPEED3)?.to_bits()
                );
                assert_eq!(calc(e::FLG_SWIEPH | e::FLG_SIDEREAL)?.values[3..], [0.0; 3]);
            }
        }
    }
    Ok(())
}

#[test]
fn radians_shape_positions_not_degree_layer_house_rates() -> Result<(), e::Error> {
    let session = session()?;
    for jd in [2451545.0, 2460000.5] {
        for system in *b"PWI" {
            let degrees = session.houses_ex2(jd, 0, 41.9, 12.5, system)?;
            let radians = session.houses_ex2(jd, e::FLG_RADIANS, 41.9, 12.5, system)?;
            assert_eq!(
                degrees.cusp_speeds.map(f64::to_bits),
                radians.cusp_speeds.map(f64::to_bits)
            );
            assert_eq!(
                degrees.angle_speeds.map(f64::to_bits),
                radians.angle_speeds.map(f64::to_bits)
            );
            assert_eq!(degrees.cusps.map(f64::to_radians), radians.cusps);
            assert_eq!(degrees.angles.map(f64::to_radians), radians.angles);
        }
        let degrees = session.houses_gauquelin(jd, 0, 41.9, 12.5)?;
        let radians = session.houses_gauquelin(jd, e::FLG_RADIANS, 41.9, 12.5)?;
        assert_eq!(
            degrees.sector_speeds.map(f64::to_bits),
            radians.sector_speeds.map(f64::to_bits)
        );
        assert_eq!(
            degrees.angle_speeds.map(f64::to_bits),
            radians.angle_speeds.map(f64::to_bits)
        );
        assert_eq!(degrees.sectors.map(f64::to_radians), radians.sectors);
        assert_eq!(degrees.angles.map(f64::to_radians), radians.angles);
    }
    Ok(())
}

#[test]
fn lunar_center_sentinels_survive_ut_tt_and_output_shaping() -> Result<(), e::Error> {
    let session = session()?;
    let jd_ut = 2451545.0;
    let jd_tt = jd_ut + session.deltat(jd_ut)?;
    for body in [
        e::MEAN_NODE,
        e::TRUE_NODE,
        e::MEAN_APOG,
        e::OSCU_APOG,
        e::INTP_APOG,
        e::INTP_PERG,
    ] {
        for center in [e::FLG_HELCTR, e::FLG_BARYCTR] {
            for shape in [
                0,
                e::FLG_SPEED,
                e::FLG_SPEED3,
                e::FLG_EQUATORIAL,
                e::FLG_EQUATORIAL | e::FLG_J2000,
                e::FLG_XYZ | e::FLG_SPEED,
                e::FLG_RADIANS | e::FLG_SPEED,
            ] {
                let flags = e::FLG_SWIEPH | center | shape;
                for position in [
                    session.calc_ut(jd_ut, body, flags)?,
                    session.calc(jd_tt, body, flags)?,
                ] {
                    assert_eq!(position.values, [0.0; 6], "body={body}, flags={flags}");
                    assert_eq!(position.returned_flags & center, center);
                    assert_eq!(
                        position.returned_flags & (e::FLG_XYZ | e::FLG_RADIANS),
                        shape & (e::FLG_XYZ | e::FLG_RADIANS)
                    );
                }
            }
        }
        assert!(session.calc_ut(jd_ut, body, e::FLG_SWIEPH)?.distance() > 0.0);
    }
    Ok(())
}

#[test]
fn sidereal_zero_origins_preserve_values_and_returned_flags() -> Result<(), e::Error> {
    let session = session()?;
    let jd = 2451545.0;
    let tt = jd + session.deltat(jd)?;
    let cases = [
        (e::EARTH, 0),
        (e::EARTH, e::FLG_TOPOCTR),
        (e::SUN, e::FLG_HELCTR),
        (e::MEAN_NODE, e::FLG_HELCTR),
        (e::TRUE_NODE, e::FLG_BARYCTR),
        (e::MEAN_APOG, e::FLG_HELCTR),
        (e::OSCU_APOG, e::FLG_BARYCTR),
        (e::INTP_APOG, e::FLG_HELCTR),
        (e::INTP_PERG, e::FLG_BARYCTR),
    ];
    for (body, center) in cases {
        for speed in [0, e::FLG_SPEED, e::FLG_SPEED3] {
            let flags = e::FLG_SWIEPH | e::FLG_SIDEREAL | center | speed;
            for position in [
                session.calc_ut(jd, body, flags)?,
                session.calc(tt, body, flags)?,
            ] {
                assert_eq!(position.values, [0.0; 6], "body={body}, flags={flags}");
                assert_ne!(position.returned_flags & e::FLG_SIDEREAL, 0);
                assert_ne!(position.returned_flags & e::FLG_NONUT, 0);
                assert_eq!(
                    position.returned_flags & (e::FLG_SPEED | e::FLG_SPEED3),
                    speed
                );
            }
        }
    }
    Ok(())
}

#[test]
fn base_sidereal_echo_preserves_shape_and_nonlongitude_components() -> Result<(), e::Error> {
    let session = session()?;
    for body in [e::MARS, e::MEAN_NODE] {
        for speed in [0, e::FLG_SPEED, e::FLG_SPEED3] {
            let flags = e::FLG_SWIEPH | e::FLG_SIDEREAL | speed;
            for is_ut in [false, true] {
                let calc = |flags| {
                    if is_ut {
                        session.calc_ut(2451545.0, body, flags)
                    } else {
                        session.calc(2451545.0, body, flags)
                    }
                };
                let position = calc(flags)?;
                assert_eq!(position.returned_flags, flags | e::FLG_NONUT);
                assert_eq!(
                    position.returned_flags & (e::FLG_SPEED | e::FLG_SPEED3),
                    speed
                );
                // An implied NONUT echo does not alter the requested output
                // shape. The position offset and longitude-rate stencil are
                // provider operations, so public ayanamsha calls at separately
                // reconstructed epochs are not an exact position oracle.
                let tropical = calc(flags & !e::FLG_SIDEREAL)?;
                let without_speed = calc(e::FLG_SWIEPH | e::FLG_SIDEREAL)?;
                let words = position.values.map(f64::to_bits);
                let tropical_words = tropical.values.map(f64::to_bits);
                assert_eq!(words[..3], without_speed.values.map(f64::to_bits)[..3]);
                if speed == 0 {
                    assert_eq!(words[3..], [0; 3]);
                }
                assert_eq!(words[1..3], tropical_words[1..3]);
                assert_eq!(words[4..], tropical_words[4..]);
            }
        }
    }
    Ok(())
}

#[test]
fn twelve_cusp_ex2_refuses_gauquelin_before_data_acquisition() -> Result<(), e::Error> {
    let unbound = e::SessionBuilder::isolated().build()?;
    for system in *b"Gg" {
        assert!(matches!(
            unbound.houses_ex2(2451545.0, 0, 41.9, 12.5, system),
            Err(e::Error::UnsupportedFlagCombination { .. })
        ));
    }
    Ok(())
}
