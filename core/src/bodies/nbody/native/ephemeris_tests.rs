// SPDX-License-Identifier: AGPL-3.0
//! Prepared common-frame/ownership checks using deliberately synthetic files.
use super::super::classic_de::tests::Fixture;
use super::*;
use std::sync::atomic::{AtomicUsize, Ordering};

struct Temporary(std::path::PathBuf);
impl Temporary {
    fn new(bytes: &[u8]) -> Self {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let path = std::env::temp_dir().join(format!(
            "libephemeris-force-frame-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::write(&path, bytes).unwrap();
        Self(path)
    }
}
impl Drop for Temporary {
    fn drop(&mut self) {
        std::fs::remove_file(&self.0).unwrap();
    }
}
fn same(actual: Source, expected: Source) {
    assert_eq!(actual.point.body, expected.point.body);
    assert_eq!(actual.point.position, expected.point.position);
    assert_eq!(
        actual.point.gravitational_parameter,
        expected.point.gravitational_parameter
    );
    assert_eq!(actual.velocity, expected.velocity);
    assert_eq!(
        actual.post_newtonian_acceleration,
        expected.post_newtonian_acceleration
    );
}

#[test]
fn planetary_frame_retains_identifiers_masses_accelerations_and_scalar_values() {
    let fixture = Fixture::force_ready();
    let (_reader, temporary) = fixture.open();
    let owner = Ephemeris::open(temporary.path(), None).unwrap();
    assert_eq!(owner.particle_count(), 11);
    assert_eq!(owner.coverage(), [(1000.0, 1064.0)]);
    let mut frame = Vec::new();
    owner.frame([1008.0, 0.0], &mut frame).unwrap();
    for (index, (source, target)) in frame.iter().copied().zip(MAJOR_TARGETS).enumerate() {
        assert_eq!(source.point.body, target);
        assert_eq!(
            source.point.gravitational_parameter,
            owner.catalogue.major_masses[index]
        );
        assert!(source.post_newtonian_acceleration.is_some());
        same(source, owner.particle(index, [1008.0, 0.0]).unwrap());
    }
    assert_eq!(frame[0].point.position, [48.5, 49.0, 49.5]);
}

#[test]
fn asteroid_mass_identity_and_same_epoch_declared_center_translation_are_explicit() {
    let fixture = Fixture::force_ready();
    let (_reader, planets) = fixture.open();
    let init = (1000.0 - J2000_JD_TDB) * SECONDS_PER_DAY;
    let mut bytes = crate::spk::force_test_image(
        crate::spk::DafEndian::Little,
        crate::spk::TYPE_CHEBYSHEV_POS,
        2,
        init,
        32.0 * SECONDS_PER_DAY,
    );
    // Reverse the two NAIF identities without moving coefficients. Mass
    // assignment must use MA0002/MA0001, not a fixed file-ordinal table.
    for (row, target) in [2_000_002_i32, 2_000_001].into_iter().enumerate() {
        let offset = 1024 + (3 + row * 5) * 8 + 16;
        bytes[offset..offset + 4].copy_from_slice(&target.to_le_bytes());
    }
    let asteroids = Temporary::new(&bytes);
    let owner = Ephemeris::open(planets.path(), Some(&asteroids.0)).unwrap();
    assert_eq!(owner.particle_count(), 13);
    assert_eq!(owner.coverage(), [(1000.0, 1064.0)]);
    let mut frame = Vec::with_capacity(owner.particle_count());
    owner.frame([1008.0, 0.0], &mut frame).unwrap();
    assert_eq!(frame[11].point.body, 2_000_002);
    assert_eq!(frame[11].point.gravitational_parameter, 0.002);
    assert_eq!(frame[12].point.body, 2_000_001);
    assert_eq!(frame[12].point.gravitational_parameter, 0.001);
    assert_eq!(frame[11].point.position[0], 48.8125);
    assert_eq!(frame[11].velocity, [0.125; 3]);
    assert!(frame[11].post_newtonian_acceleration.is_none());
    for (index, source) in frame.into_iter().enumerate() {
        same(source, owner.particle(index, [1008.0, 0.0]).unwrap());
    }
}

#[test]
fn curated_relative_seeds_share_bounded_force_reader_and_original_units() {
    for endian in [crate::spk::DafEndian::Little, crate::spk::DafEndian::Big] {
        for kind in [
            crate::spk::TYPE_CHEBYSHEV_POS,
            crate::spk::TYPE_CHEBYSHEV_POS_VEL,
        ] {
            let fixture = Fixture::force_ready();
            let (_reader, planets) = fixture.open();
            let bytes = crate::spk::force_test_image(
                endian,
                kind,
                2,
                (1000.0 - J2000_JD_TDB) * SECONDS_PER_DAY,
                32.0 * SECONDS_PER_DAY,
            );
            let asteroids = Temporary::new(&bytes);
            let owner = Ephemeris::open(planets.path(), Some(&asteroids.0)).unwrap();
            let mut original = crate::spk::SpkRegistry::new();
            let key = original.register_file(&asteroids.0).unwrap();
            original.register_body(2_000_001, &key, 2_000_001).unwrap();
            for epoch in [1000.0, 1008.0, 1032.0, 1064.0] {
                let actual = owner
                    .heliocentric_seed(2_000_001, [epoch, 0.0])
                    .unwrap()
                    .unwrap();
                let (expected, descriptor) = original
                    .eval_body_relative_two_part_jd(2_000_001, epoch, 0.0)
                    .unwrap();
                assert_eq!(descriptor.center, crate::spk::CENTER_SUN);
                assert_eq!(actual.position_au(), expected.position_au());
                assert_eq!(actual.velocity_au_per_day(), expected.velocity_au_per_day());
            }
            assert!(
                owner
                    .heliocentric_seed(2_000_099, [1008.0, 0.0])
                    .unwrap()
                    .is_none()
            );
            assert!(owner.heliocentric_seed(2_000_001, [999.0, 0.0]).is_err());
        }
    }
}

#[test]
fn curated_seed_does_not_assume_a_heliocentric_center_or_fabricate_missing_data() {
    let fixture = Fixture::force_ready();
    let (_reader, planets) = fixture.open();
    let absent = Ephemeris::open(planets.path(), None).unwrap();
    assert!(
        absent
            .heliocentric_seed(2_000_001, [1008.0, 0.0])
            .unwrap()
            .is_none()
    );
    let mut bytes = crate::spk::force_test_image(
        crate::spk::DafEndian::Little,
        crate::spk::TYPE_CHEBYSHEV_POS,
        1,
        (1000.0 - J2000_JD_TDB) * SECONDS_PER_DAY,
        32.0 * SECONDS_PER_DAY,
    );
    // The stored frame is valid for forcing, but explicit curated seeds
    // must preserve their separate heliocentric initialization contract.
    let center_offset = 1024 + 3 * 8 + 20;
    bytes[center_offset..center_offset + 4].copy_from_slice(&3_i32.to_le_bytes());
    let asteroids = Temporary::new(&bytes);
    let owner = Ephemeris::open(planets.path(), Some(&asteroids.0)).unwrap();
    assert!(matches!(
        owner.heliocentric_seed(2_000_001, [1008.0, 0.0]),
        Err(Error::SpkAmbiguousCenter { target: 2_000_001 })
    ));
}

#[test]
fn numbered_and_provisional_kbo_masses_use_their_selected_header_aliases() {
    // JPL SB441 Table1 identifies Eris2136199 and provisional2010 RF43
    // 3545742. Exact published GM/header word joins establish MA8001/MA8027;
    // these synthetic values test identity only, never fit an orbital result.
    for (target, label) in [(2_136_199_i32, "MA8001"), (3_545_742, "MA8027")] {
        let fixture = Fixture::force_ready().with_constant(label, 0.003);
        let (_reader, planets) = fixture.open();
        let mut bytes = crate::spk::force_test_image(
            crate::spk::DafEndian::Little,
            crate::spk::TYPE_CHEBYSHEV_POS,
            1,
            (1000.0 - J2000_JD_TDB) * SECONDS_PER_DAY,
            32.0 * SECONDS_PER_DAY,
        );
        let target_offset = 1024 + 3 * 8 + 16;
        bytes[target_offset..target_offset + 4].copy_from_slice(&target.to_le_bytes());
        let asteroids = Temporary::new(&bytes);
        let owner = Ephemeris::open(planets.path(), Some(&asteroids.0)).unwrap();
        let actual = owner.particle(11, [1008.0, 0.0]).unwrap();
        assert_eq!(actual.point.body, target);
        assert_eq!(actual.point.gravitational_parameter, 0.003);
    }
}

fn selected_asteroids(targets: &[i32]) -> Temporary {
    let mut bytes = crate::spk::force_test_image(
        crate::spk::DafEndian::Little,
        crate::spk::TYPE_CHEBYSHEV_POS,
        targets.len(),
        (1000.0 - J2000_JD_TDB) * SECONDS_PER_DAY,
        32.0 * SECONDS_PER_DAY,
    );
    for (index, target) in targets.iter().enumerate() {
        // DAF summary records hold 25 descriptors in this synthetic image.
        let offset = (1 + 2 * (index / 25)) * 1024 + (3 + (index % 25) * 5) * 8 + 16;
        bytes[offset..offset + 4].copy_from_slice(&target.to_le_bytes());
    }
    Temporary::new(&bytes)
}

#[test]
fn entire_sb441_kbo_inventory_uses_selected_masses_independent_of_record_order() {
    // Physical SPICE IDs in JPL Table 1's KBO order (MA8001 through MA8030).
    // Distinct exactly representable synthetic GMs discriminate every alias.
    let targets = [
        2136199, 2136108, 2090377, 2136472, 2050000, 2084522, 2090482, 2020000, 2055637, 2028978,
        2307261, 2174567, 3361580, 3308265, 2055565, 2145452, 2090568, 2208996, 2225088, 2019521,
        2120347, 2278361, 3525142, 2230965, 2042301, 2455502, 3545742, 2523639, 2528381, 3515022,
    ];
    let mut fixture = Fixture::force_ready();
    for index in 0..targets.len() {
        fixture =
            fixture.with_constant(&format!("MA{}", 8001 + index), (index + 1) as f64 / 1024.0);
    }
    let (_reader, planets) = fixture.open();
    for reverse in [false, true] {
        let mut order = targets;
        if reverse {
            order.reverse();
        }
        let asteroids = selected_asteroids(&order);
        let owner = Ephemeris::open(planets.path(), Some(&asteroids.0)).unwrap();
        assert_eq!(owner.particle_count(), 41);
        for (index, target) in order.into_iter().enumerate() {
            let source = owner.particle(11 + index, [1008.0, 0.0]).unwrap();
            let published_index = targets.iter().position(|value| *value == target).unwrap();
            assert_eq!(source.point.body, target);
            assert_eq!(
                source.point.gravitational_parameter,
                (published_index + 1) as f64 / 1024.0,
            );
        }
    }
}

#[test]
fn numbered_asteroids_cannot_borrow_reserved_kbo_mass_labels() {
    for number in 8001..=8030 {
        let fixture = Fixture::force_ready().with_constant(&format!("MA{number}"), 0.003);
        let (_reader, planets) = fixture.open();
        let asteroids = selected_asteroids(&[2_000_000 + number]);
        assert!(matches!(
            Ephemeris::open(planets.path(), Some(&asteroids.0)),
            Err(EphemerisError::De(DeError::UnsupportedTarget(target)))
                if target == 2_000_000 + number
        ));
    }
    let fixture = Fixture::force_ready().with_constant("MA8001", 0.003);
    let (_reader, planets) = fixture.open();
    let asteroids = selected_asteroids(&[2_136_199, 2_008_001]);
    assert!(matches!(
        Ephemeris::open(planets.path(), Some(&asteroids.0)),
        Err(EphemerisError::De(DeError::UnsupportedTarget(2_008_001)))
    ));
    for number in [8000, 8031] {
        let fixture = Fixture::force_ready().with_constant(&format!("MA{number}"), 0.005);
        let (_reader, planets) = fixture.open();
        let asteroids = selected_asteroids(&[2_000_000 + number]);
        let owner = Ephemeris::open(planets.path(), Some(&asteroids.0)).unwrap();
        assert_eq!(
            owner
                .particle(11, [1008.0, 0.0])
                .unwrap()
                .point
                .gravitational_parameter,
            0.005,
        );
    }
}

#[test]
fn missing_or_invalid_selected_kbo_mass_is_an_error() {
    let fixture = Fixture::force_ready();
    let (_reader, planets) = fixture.open();
    let asteroids = selected_asteroids(&[2_136_199]);
    assert!(matches!(
        Ephemeris::open(planets.path(), Some(&asteroids.0)),
        Err(EphemerisError::De(DeError::MissingConstant(name))) if name == "MA8001"
    ));
    let fixture = fixture.with_constant("MA8001", -0.003);
    let (_reader, planets) = fixture.open();
    assert!(matches!(
        Ephemeris::open(planets.path(), Some(&asteroids.0)),
        Err(EphemerisError::De(DeError::Malformed("asteroid GM")))
    ));
}

#[test]
fn coverage_intersections_keep_gaps_and_closed_single_epoch_contacts() {
    assert_eq!(
        intersection(
            &[(0.0, 5.0), (8.0, 10.0)],
            &[(2.0, 3.0), (4.0, 8.0), (9.0, 12.0)]
        ),
        [(2.0, 3.0), (4.0, 5.0), (8.0, 8.0), (9.0, 10.0)]
    );
    assert!(intersection(&[(0.0, 1.0)], &[(2.0, 3.0)]).is_empty());
}

#[test]
fn metadata_indices_and_uncovered_epochs_do_not_publish_a_partial_success() {
    let minimal = Fixture::little();
    let (_reader, temporary) = minimal.open();
    assert!(matches!(
        Ephemeris::open(temporary.path(), None),
        Err(EphemerisError::De(DeError::MissingConstant(_)))
    ));
    let fixture = Fixture::force_ready();
    let (_reader, temporary) = fixture.open();
    let owner = Ephemeris::open(temporary.path(), None).unwrap();
    assert!(matches!(
        owner.particle(11, [1008.0, 0.0]),
        Err(EphemerisError::InvalidParticle(11))
    ));
    let mut frame = Vec::new();
    assert!(matches!(
        owner.frame([999.0, 0.0], &mut frame),
        Err(EphemerisError::De(DeError::OutsideCoverage))
    ));
    let absent = temporary.path().with_extension("absent-spk");
    assert!(matches!(
        Ephemeris::open(temporary.path(), Some(&absent)),
        Err(EphemerisError::Spk(Error::SpkFileNotFound { .. }))
    ));
}
