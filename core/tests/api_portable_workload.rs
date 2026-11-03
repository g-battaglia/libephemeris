//! Standalone producer Session workload over explicit LEB data.
//! Consumer interoperability is tested through adapters, not mirrored producer
//! facades. This test has no Swiss dependency or cross-engine reference vectors.
use libephemeris_core::native as ephe;

fn workload(session: &ephe::Session) -> Result<(), ephe::Error> {
    let jd = ephe::julday(2000, 1, 1, 12.0, ephe::Calendar::Gregorian)?.value;
    let flags = ephe::FLG_SWIEPH;
    for shape in [ephe::FLG_RADIANS, ephe::FLG_XYZ, ephe::FLG_EQUATORIAL] {
        assert!(session.solcross_ut(45.0, jd, flags | shape).is_err());
        assert!(session.mooncross_node(jd, flags | shape).is_err());
        assert!(
            session
                .helio_cross(ephe::MARS, 45.0, jd, flags | shape, 1)
                .is_err()
        );
    }
    let hit = session.solcross_ut(45.0, jd, flags)?;
    let position = session.calc_ut(hit.time, ephe::SUN, flags)?;
    assert!(ephe::difdeg2n(position.values[0], 45.0).abs() < 1e-5);
    let hit = session.mooncross(45.0, jd, flags)?;
    let position = session.calc(hit.time, ephe::MOON, flags)?;
    assert!(ephe::difdeg2n(position.values[0], 45.0).abs() < 1e-5);
    assert!(session.solcross_ut(0.0, jd, flags)?.time >= jd);
    assert!(session.solcross(0.0, jd, flags)?.time >= jd);
    assert!(session.mooncross_ut(0.0, jd, flags)?.time >= jd);
    assert!(session.mooncross(0.0, jd, flags)?.time >= jd);
    assert!(session.mooncross_node_ut(jd, flags)?.longitude.is_finite());
    assert!(session.mooncross_node(jd, flags)?.latitude.is_finite());
    assert!(session.helio_cross_ut(ephe::MARS, 0.0, jd, flags, 1)?.time >= jd);
    assert!(session.helio_cross(ephe::MARS, 0.0, jd, flags, 1)?.time >= jd);
    let date = ephe::revjul(jd, ephe::Calendar::Gregorian)?;
    assert_eq!((date.year, date.month, date.day), (2000, 1, 1));
    for body in [
        ephe::SUN,
        ephe::MOON,
        ephe::MARS,
        ephe::MEAN_NODE,
        ephe::TRUE_NODE,
        ephe::MEAN_APOG,
        ephe::OSCU_APOG,
    ] {
        for flags in [
            ephe::FLG_SWIEPH | ephe::FLG_SPEED,
            ephe::FLG_SWIEPH | ephe::FLG_EQUATORIAL,
            ephe::FLG_SWIEPH | ephe::FLG_XYZ | ephe::FLG_SPEED,
        ] {
            let p = session.calc_ut(jd, body, flags)?;
            assert!(p.values.iter().all(|v| v.is_finite()));
            assert_eq!(p.returned_flags & ephe::FLG_XYZ, flags & ephe::FLG_XYZ);
        }
    }
    for body in [ephe::SUN, ephe::MOON, ephe::MARS] {
        let p = session.pheno_ut(jd, body, ephe::FLG_SWIEPH | ephe::FLG_SPEED)?;
        assert_eq!(p.values.len(), 20);
        assert!(p.values.iter().all(|v| v.is_finite()));
        assert!(
            session
                .pheno(jd + session.deltat(jd)?, body, ephe::FLG_SWIEPH)?
                .values
                .iter()
                .all(|v| v.is_finite())
        );
    }
    let tt = jd + session.deltat(jd)?;
    assert!(
        session
            .calc(tt, ephe::SUN, ephe::FLG_SWIEPH)?
            .values
            .iter()
            .all(|v| v.is_finite())
    );
    assert!(
        session
            .calc_pctr(tt, ephe::MOON, ephe::EARTH, ephe::FLG_SWIEPH)?
            .values
            .iter()
            .all(|v| v.is_finite())
    );
    for hsys in *b"PWI" {
        let h = session.houses_ex2(jd, 0, 41.9, 12.5, hsys)?;
        assert_eq!(h.cusps.len(), 12);
        assert_eq!(h.angles.len(), 8);
        assert!(
            h.cusp_speeds
                .iter()
                .chain(h.angle_speeds.iter())
                .all(|v| v.is_finite())
        );
        assert_eq!(session.houses(jd, 41.9, 12.5, hsys)?.cusps, h.cusps);
    }
    let g = session.houses_gauquelin(jd, 0, 41.9, 12.5)?;
    assert_eq!(g.sectors.len(), 36);
    assert!(g.sector_speeds.iter().all(|v| v.is_finite()));
    let horizontal = session.azalt(
        jd,
        ephe::EQU2HOR,
        12.5,
        41.9,
        0.0,
        1013.25,
        15.0,
        200.0,
        60.0,
    )?;
    let celestial = session.azalt_rev(
        jd,
        ephe::HOR2EQU,
        12.5,
        41.9,
        0.0,
        horizontal.azimuth,
        horizontal.true_altitude,
    )?;
    assert!((celestial.longitude_or_ra - 200.0).abs() < 1e-8);
    assert!((celestial.latitude_or_dec - 60.0).abs() < 1e-8);
    assert!(session.sidtime(jd)?.is_finite());
    assert!(session.sidtime0(jd, 23.4, 0.0)?.is_finite());
    assert!(session.get_ayanamsa_ut(jd)?.is_finite());
    assert!(session.get_ayanamsa(tt)?.is_finite());
    assert!(session.get_ayanamsa_ex_ut(jd, 0)?.value.is_finite());
    assert!(session.get_ayanamsa_ex(tt, 0)?.value.is_finite());
    assert_eq!(
        session
            .calc_ut(f64::NAN, ephe::SUN, 0)
            .unwrap_err()
            .category(),
        ephe::ErrorCategory::InputValidation
    );
    assert_eq!(
        ephe::julday(2000, 1, 1, f64::NAN, ephe::Calendar::Gregorian)
            .unwrap_err()
            .category(),
        ephe::ErrorCategory::InputValidation
    );
    let a = session
        .calc_ut(jd, ephe::SUN, ephe::FLG_SWIEPH | ephe::FLG_TOPOCTR)?
        .values;
    let mut variant = session.to_builder();
    variant
        .topo(-74.0, 40.7, 20.0)
        .sid_mode(ephe::SIDM_LAHIRI, 0.0, 0.0);
    let other = variant.build()?;
    assert!(
        other
            .calc_ut(jd, ephe::SUN, ephe::FLG_SWIEPH | ephe::FLG_SIDEREAL)?
            .values
            .iter()
            .all(|v| v.is_finite())
    );
    assert_eq!(
        a.map(f64::to_bits),
        session
            .calc_ut(jd, ephe::SUN, ephe::FLG_SWIEPH | ephe::FLG_TOPOCTR)?
            .values
            .map(f64::to_bits)
    );
    Ok(())
}

#[test]
fn portable_natal_surface() -> Result<(), ephe::Error> {
    let path = std::env::var_os("KERYKEION_LEB_FILE")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| {
            std::path::PathBuf::from(std::env::var_os("HOME").unwrap())
                .join(".libephemeris/leb/base_core.leb2")
        });
    let mut builder = ephe::SessionBuilder::isolated();
    builder.leb_file(&path)?.topo(12.5, 41.9, 0.0);
    workload(&builder.build()?)
}

#[test]
fn portable_coordinate_rotations() -> Result<(), ephe::Error> {
    for (longitude, latitude, distance) in [(0.0, 0.0, 1.0), (359.9, -40.0, 3.0), (45.0, 80.0, 2.0)]
    {
        let rotated = ephe::cotrans(longitude, latitude, distance, 23.44)?;
        let restored = ephe::cotrans(rotated[0], rotated[1], rotated[2], -23.44)?;
        assert!(ephe::difdeg2n(restored[0], longitude).abs() < 1e-10);
        assert!((restored[1] - latitude).abs() < 1e-10);
        assert_eq!(restored[2], distance);
        let (position, speed) =
            ephe::cotrans_sp(longitude, latitude, distance, 0.4, -0.2, 0.01, 23.44)?;
        assert!(position.into_iter().chain(speed).all(f64::is_finite));
        assert_eq!(position[2], distance);
        assert_eq!(speed[2], 0.01);
    }
    assert!(ephe::cotrans(f64::NAN, 0.0, 1.0, 23.44).is_err());
    assert!(ephe::cotrans_sp(0.0, 0.0, 1.0, f64::INFINITY, 0.0, 0.0, 23.44).is_err());
    Ok(())
}
