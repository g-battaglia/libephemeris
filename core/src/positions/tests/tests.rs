use super::*;
use crate::leb::{
    HEADER_SIZE, LEB1_BODY_ENTRY_SIZE, LEB2_BODY_ENTRY_SIZE, MAGIC_LEB1, SECTION_DIR_SIZE,
};

fn push_u32(v: &mut Vec<u8>, x: u32) {
    v.extend_from_slice(&x.to_le_bytes());
}

fn push_i32(v: &mut Vec<u8>, x: i32) {
    v.extend_from_slice(&x.to_le_bytes());
}

fn push_u64(v: &mut Vec<u8>, x: u64) {
    v.extend_from_slice(&x.to_le_bytes());
}

fn push_f64(v: &mut Vec<u8>, x: f64) {
    v.extend_from_slice(&x.to_le_bytes());
}

fn pack_header(
    magic: &[u8; 4],
    version: u32,
    sections: u32,
    bodies: u32,
    jd_start: f64,
    jd_end: f64,
) -> Vec<u8> {
    let mut v = Vec::with_capacity(HEADER_SIZE);
    v.extend_from_slice(magic);
    push_u32(&mut v, version);
    push_u32(&mut v, sections);
    push_u32(&mut v, bodies);
    push_f64(&mut v, jd_start);
    push_f64(&mut v, jd_end);
    push_f64(&mut v, 0.0);
    push_u32(&mut v, 0);
    v.extend_from_slice(&[0u8; 20]);
    debug_assert_eq!(v.len(), HEADER_SIZE);
    v
}

fn pack_section(id: u32, offset: u64, size: u64) -> Vec<u8> {
    let mut v = Vec::with_capacity(SECTION_DIR_SIZE);
    push_u32(&mut v, id);
    push_u32(&mut v, 0);
    push_u64(&mut v, offset);
    push_u64(&mut v, size);
    v
}

fn pack_body(body_id: i32, coord: u32, data_offset: u64) -> Vec<u8> {
    // 2 segments over [0, 2], 1-day interval, degree 1, 3 components.
    let mut v = Vec::with_capacity(LEB1_BODY_ENTRY_SIZE);
    push_i32(&mut v, body_id);
    push_u32(&mut v, coord);
    push_u32(&mut v, 2);
    push_f64(&mut v, 0.0);
    push_f64(&mut v, 2.0);
    push_f64(&mut v, 1.0);
    push_u32(&mut v, 1);
    push_u32(&mut v, 3);
    push_u64(&mut v, data_offset);
    v
}

/// One-body LEB1 image. Blob per segment/component: `[c0, c1]` with
/// comp 0 constant 5, comp 1 slope 1, comp 2 zero — so at `jd = 0`
/// (`tau = -1`) the stored state is `pos = [5, -1, 0]`,
/// `vel = [0, 2, 0]` bit-exactly (velocity scale `2/interval`).
fn fixture(body_id: i32, coord: u32) -> Vec<u8> {
    let mut file = pack_header(&MAGIC_LEB1, 1, 1, 1, 0.0, 2.0);
    let body_at = HEADER_SIZE + SECTION_DIR_SIZE;
    let blob_at = body_at + LEB1_BODY_ENTRY_SIZE;
    file.extend_from_slice(&pack_section(0, body_at as u64, 52));
    file.extend_from_slice(&pack_body(body_id, coord, blob_at as u64));
    for _ in 0..2 {
        for comp in 0..3 {
            push_f64(&mut file, if comp == 0 { 5.0 } else { 0.0 });
            push_f64(&mut file, if comp == 1 { 1.0 } else { 0.0 });
        }
    }
    debug_assert_eq!(file.len(), blob_at + 96);
    file
}

fn composite_of(file: &[u8]) -> (FileIndex, CompositeIndex<'_>) {
    let index = FileIndex::decode(file).expect("fixture must decode");
    // Reborrow dance: decode borrows nothing (it copies metadata), so
    // a second decode for the composite is borrow-clean.
    let owned = FileIndex::decode(file).expect("fixture must decode");
    let composite = CompositeIndex::new(vec![(owned, file)]).expect("one reader must build");
    (index, composite)
}

#[test]
fn frame_classification_pins_all_channels() {
    assert_eq!(
        StateFrame::from_coord_type(COORD_ICRS_BARY).unwrap(),
        StateFrame::IcrsBarycentric
    );
    assert_eq!(
        StateFrame::from_coord_type(COORD_ICRS_BARY_SYSTEM).unwrap(),
        StateFrame::IcrsSystemBarycenter
    );
    assert_eq!(
        StateFrame::from_coord_type(COORD_ECLIPTIC).unwrap(),
        StateFrame::EclipticDirect
    );
    assert_eq!(
        StateFrame::from_coord_type(COORD_GEO_ECLIPTIC).unwrap(),
        StateFrame::GeocentricEcliptic
    );
    assert_eq!(
        StateFrame::from_coord_type(COORD_HELIO_ECL_RETIRED).unwrap(),
        StateFrame::HeliocentricEclipticRetired
    );
    assert!(matches!(
        StateFrame::from_coord_type(LEB2_BODY_ENTRY_SIZE as u32).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    assert!(StateFrame::IcrsBarycentric.is_icrs());
    assert!(StateFrame::IcrsSystemBarycenter.is_icrs());
    assert!(!StateFrame::EclipticDirect.is_icrs());
    assert_eq!(
        StateFrame::IcrsBarycentric.position_units(),
        ["AU", "AU", "AU"]
    );
    assert_eq!(
        StateFrame::EclipticDirect.position_units(),
        ["deg", "deg", "AU"]
    );
    assert_eq!(
        StateFrame::GeocentricEcliptic.velocity_units(),
        ["deg/day", "deg/day", "AU/day"]
    );
}

#[test]
fn icrs_state_passes_through_with_identity_metadata() {
    let file = fixture(0, COORD_ICRS_BARY);
    let (_index, composite) = composite_of(&file);
    let st = composite.geometric_state(0, 0.0).expect("in-range");
    assert_eq!(st.body_id(), 0);
    assert_eq!(st.jd_tt(), 0.0);
    assert_eq!(st.pos(), [5.0, -1.0, 0.0]);
    assert_eq!(st.vel(), [0.0, 2.0, 0.0]);
    assert_eq!(st.frame(), StateFrame::IcrsBarycentric);
    assert!(!st.is_planet_center());
    assert_eq!(st.tier(), None);
    assert!(composite.has_body(0));
    assert!(!composite.has_body(99));
    assert_eq!(composite.body_coverage(0), Some((0.0, 2.0)));
    assert_eq!(composite.serving_coord_type(0), Some(COORD_ICRS_BARY));
    assert_eq!(composite.provider_name(), "sealed-leb/composite");
}

#[test]
fn system_barycenter_retained_without_offset_provider() {
    let file = fixture(5, COORD_ICRS_BARY_SYSTEM);
    let (_index, composite) = composite_of(&file);
    let st = composite.geometric_state(5, 0.0).expect("in-range");
    assert_eq!(st.frame(), StateFrame::IcrsSystemBarycenter);
    // Explicit fallback: stored barycenter untouched, planet-center
    // flag clear.
    assert_eq!(st.pos(), [5.0, -1.0, 0.0]);
    assert_eq!(st.vel(), [0.0, 2.0, 0.0]);
    assert!(!st.is_planet_center());
}

#[test]
fn system_barycenter_offset_paths() {
    let file = fixture(5, COORD_ICRS_BARY_SYSTEM);
    let (_index, composite) = composite_of(&file);

    let st = composite
        .geometric_state_with_offset(5, 0.5, |_| Some([0.0; 3]))
        .expect("in-range");
    assert!(!st.is_planet_center());
    let retained = composite.geometric_state(5, 0.5).expect("in-range");
    assert_eq!(st.pos(), retained.pos());
    assert_eq!(st.vel(), retained.vel());
    // Constant nonzero offset: position shifts bit-exactly, and the
    // centered difference of a constant is zero so velocity is
    // bitwise unchanged; the planet-center flag is set.
    let st = composite
        .geometric_state_with_offset(5, 0.5, |_| Some([1.0, -2.0, 0.25]))
        .expect("in-range");
    assert!(st.is_planet_center());
    assert_eq!(st.pos()[0], retained.pos()[0] + 1.0);
    assert_eq!(st.pos()[1], retained.pos()[1] - 2.0);
    assert_eq!(st.pos()[2], retained.pos()[2] + 0.25);
    assert_eq!(st.vel(), retained.vel());
}

#[test]
fn ecliptic_channel_carries_frame_not_xyz() {
    let file = fixture(10, COORD_ECLIPTIC);
    let (_index, composite) = composite_of(&file);
    let st = composite.geometric_state(10, 0.0).expect("in-range");
    assert_eq!(st.frame(), StateFrame::EclipticDirect);
    assert_eq!(st.pos(), [5.0, -1.0, 0.0]);
    assert!(!st.is_planet_center());
}

#[test]
fn sealed_errors_name_the_requested_body() {
    let file = fixture(0, COORD_ICRS_BARY);
    let (_index, composite) = composite_of(&file);
    // Unserved body: no fallback, typed absence.
    assert_eq!(
        composite.geometric_state(3, 0.5).unwrap_err(),
        Error::LebBodyNotFound { body_id: 3 }
    );
    // Out-of-coverage epoch: range miss naming the REQUESTED body
    // (the sealed `_escalate_sealed_range_miss` concern for core
    // bodies is met structurally here: the error carries body 0,
    // never a support id).
    assert_eq!(
        composite.geometric_state(0, 7.5).unwrap_err(),
        Error::LebOutOfRange {
            body_id: 0,
            jd: 7.5,
            start: 0.0,
            end: 2.0,
        }
    );
    // Non-finite epoch is a typed refusal, never garbage.
    assert!(matches!(
        composite.geometric_state(0, f64::NAN).unwrap_err(),
        Error::NonFiniteArgument { .. }
    ));
    // Entry-shaped errors surface before evaluation: a missing
    // channel snapshot agrees with the eval miss.
    assert_eq!(
        composite.serving_entry_snapshot(3).unwrap_err(),
        Error::LebBodyNotFound { body_id: 3 }
    );
}

#[test]
fn two_part_agrees_with_scalar_and_keeps_metadata() {
    let file = fixture(0, COORD_ICRS_BARY);
    let (_index, composite) = composite_of(&file);
    let epoch = TwoPartEpoch::split(0.5).expect("finite splits");
    let a = composite.geometric_state(0, 0.5).expect("scalar");
    let b = composite
        .geometric_state_two_part(0, epoch)
        .expect("two-part");
    assert_eq!(b.body_id(), a.body_id());
    assert_eq!(b.frame(), a.frame());
    assert_eq!(b.is_planet_center(), a.is_planet_center());
    for k in 0..3 {
        assert!(
            (b.pos()[k] - a.pos()[k]).abs() <= 1e-12,
            "pos[{k}] diverged: {} vs {}",
            b.pos()[k],
            a.pos()[k]
        );
        assert!(
            (b.vel()[k] - a.vel()[k]).abs() <= 1e-12,
            "vel[{k}] diverged"
        );
    }
}

#[test]
fn ut_entry_applies_caller_delta_t() {
    let file = fixture(0, COORD_ICRS_BARY);
    let (_index, composite) = composite_of(&file);
    // jd_tt = 1.0 + 0.0008 exactly representable? 1.0008 is not;
    // compare against the TT path at the same collapsed argument.
    let jd_tt = 1.0 + 0.0008;
    let a = composite.geometric_state(0, jd_tt).expect("tt");
    let b = composite.geometric_state_ut(0, 1.0, 0.0008).expect("ut");
    assert_eq!(a.pos(), b.pos());
    assert_eq!(a.vel(), b.vel());
    assert_eq!(b.jd_tt(), jd_tt);
    assert!(matches!(
        composite
            .geometric_state_ut(0, f64::INFINITY, 0.0008)
            .unwrap_err(),
        Error::NonFiniteArgument { .. }
    ));
    assert!(matches!(
        composite.geometric_state_ut(0, 1.0, f64::NAN).unwrap_err(),
        Error::NonFiniteArgument { .. }
    ));
}

#[test]
fn tiered_provider_selects_covering_tier_and_reports_it() {
    // Base covers [0, 1], medium covers [1, 2]: both carry body 3
    // with the same blob so selection is observable via the tier tag.
    fn ranged_fixture(body_id: i32, start: f64, end: f64) -> Vec<u8> {
        let mut file = pack_header(&MAGIC_LEB1, 1, 1, 1, start, end);
        let body_at = HEADER_SIZE + SECTION_DIR_SIZE;
        let blob_at = body_at + LEB1_BODY_ENTRY_SIZE;
        file.extend_from_slice(&pack_section(0, body_at as u64, 52));
        let mut entry = Vec::with_capacity(LEB1_BODY_ENTRY_SIZE);
        push_i32(&mut entry, body_id);
        push_u32(&mut entry, COORD_ICRS_BARY);
        push_u32(&mut entry, 1);
        push_f64(&mut entry, start);
        push_f64(&mut entry, end);
        push_f64(&mut entry, end - start);
        push_u32(&mut entry, 1);
        push_u32(&mut entry, 3);
        push_u64(&mut entry, blob_at as u64);
        file.extend_from_slice(&entry);
        for comp in 0..3 {
            push_f64(&mut file, if comp == 0 { 5.0 } else { 0.0 });
            push_f64(&mut file, if comp == 1 { 1.0 } else { 0.0 });
        }
        file
    }
    let base_bytes = ranged_fixture(3, 0.0, 1.0);
    let med_bytes = ranged_fixture(3, 1.0, 2.0);
    let base_index = FileIndex::decode(&base_bytes).expect("base decodes");
    let med_index = FileIndex::decode(&med_bytes).expect("medium decodes");
    let base = CompositeIndex::new(vec![(base_index, &base_bytes[..])]).expect("base composite");
    let med = CompositeIndex::new(vec![(med_index, &med_bytes[..])]).expect("medium composite");
    let tiered =
        TieredIndex::new(vec![(LebTier::Base, base), (LebTier::Medium, med)]).expect("tiers build");
    // Date-selected tier wins per body: 0.25 -> base, 1.75 -> medium.
    let at_base = tiered.geometric_state(3, 0.25).expect("base covers");
    assert_eq!(at_base.tier(), Some(LebTier::Base));
    let at_med = tiered.geometric_state(3, 1.75).expect("medium covers");
    assert_eq!(at_med.tier(), Some(LebTier::Medium));
    assert_eq!(tiered.provider_name(), "sealed-leb/tiered");
    // Outside every tier: the highest-priority candidate's range
    // error surfaces for the sealed dispatcher (LEB-09) — sealed
    // means reported, not replaced.
    assert!(matches!(
        tiered.geometric_state(3, 9.0).unwrap_err(),
        Error::LebOutOfRange { body_id: 3, .. }
    ));
    assert_eq!(tiered.body_coverage(3), Some((0.0, 2.0)));
}

#[test]
fn single_provider_serves_monolith_directly() {
    let file = fixture(1, COORD_ICRS_BARY);
    let index = FileIndex::decode(&file).expect("fixture must decode");
    let single = SingleLeb::new(&index, &file);
    let st = single.geometric_state(1, 0.0).expect("in-range");
    assert_eq!(st.pos(), [5.0, -1.0, 0.0]);
    assert_eq!(st.vel(), [0.0, 2.0, 0.0]);
    assert_eq!(single.provider_name(), "sealed-leb/single");
    assert_eq!(
        single.geometric_state(2, 0.0).unwrap_err(),
        Error::LebBodyNotFound { body_id: 2 }
    );
}

#[test]
fn batch_states_match_repeated_scalars_with_per_epoch_tiers() {
    // LEB-12: one scalar `geometric_state` per epoch, in input
    // order, over a batch straddling endpoints, the interior seam
    // 1.0, adjacent float words, and a tier boundary.
    let file = fixture(0, COORD_ICRS_BARY);
    let index = FileIndex::decode(&file).expect("fixture must decode");
    let single = SingleLeb::new(&index, &file);
    let jds = [
        0.0,
        0.5,
        f64::next_down(1.0),
        1.0,
        f64::next_up(1.0),
        1.5,
        2.0,
    ];
    let batched = single.geometric_states(0, &jds).expect("batch");
    assert_eq!(batched.len(), jds.len());
    for (jd, st) in jds.iter().zip(batched.iter()) {
        assert_eq!(*st, single.geometric_state(0, *jd).expect("scalar"));
    }
    // Fail-fast: the out-of-range element reports the scalar error.
    assert_eq!(
        single.geometric_states(0, &[0.5, 9.0]).unwrap_err(),
        single.geometric_state(0, 9.0).unwrap_err()
    );
    // Two-part and UT batches loop their scalar entries too.
    let epochs: Vec<TwoPartEpoch> = jds
        .iter()
        .map(|jd| TwoPartEpoch::split(*jd).expect("finite splits"))
        .collect();
    let batched_two = single
        .geometric_states_two_part(0, &epochs)
        .expect("two-part batch");
    for (epoch, st) in epochs.iter().zip(batched_two.iter()) {
        assert_eq!(
            *st,
            single.geometric_state_two_part(0, *epoch).expect("scalar")
        );
    }
    let batched_ut = single
        .geometric_states_ut(0, &[0.5, 1.5], 0.0008)
        .expect("ut batch");
    assert_eq!(
        batched_ut[0],
        single.geometric_state_ut(0, 0.5, 0.0008).expect("scalar")
    );
    // Tier tags travel per element across a tiered batch.
    let base_bytes = {
        let mut file = pack_header(&MAGIC_LEB1, 1, 1, 1, 0.0, 1.0);
        let body_at = HEADER_SIZE + SECTION_DIR_SIZE;
        let blob_at = body_at + LEB1_BODY_ENTRY_SIZE;
        file.extend_from_slice(&pack_section(0, body_at as u64, 52));
        let mut entry = Vec::with_capacity(LEB1_BODY_ENTRY_SIZE);
        push_i32(&mut entry, 3);
        push_u32(&mut entry, COORD_ICRS_BARY);
        push_u32(&mut entry, 1);
        push_f64(&mut entry, 0.0);
        push_f64(&mut entry, 1.0);
        push_f64(&mut entry, 1.0);
        push_u32(&mut entry, 1);
        push_u32(&mut entry, 3);
        push_u64(&mut entry, blob_at as u64);
        file.extend_from_slice(&entry);
        for comp in 0..3 {
            push_f64(&mut file, if comp == 0 { 5.0 } else { 0.0 });
            push_f64(&mut file, if comp == 1 { 1.0 } else { 0.0 });
        }
        file
    };
    let med_bytes = {
        let mut file = pack_header(&MAGIC_LEB1, 1, 1, 1, 1.0, 2.0);
        let body_at = HEADER_SIZE + SECTION_DIR_SIZE;
        let blob_at = body_at + LEB1_BODY_ENTRY_SIZE;
        file.extend_from_slice(&pack_section(0, body_at as u64, 52));
        let mut entry = Vec::with_capacity(LEB1_BODY_ENTRY_SIZE);
        push_i32(&mut entry, 3);
        push_u32(&mut entry, COORD_ICRS_BARY);
        push_u32(&mut entry, 1);
        push_f64(&mut entry, 1.0);
        push_f64(&mut entry, 2.0);
        push_f64(&mut entry, 1.0);
        push_u32(&mut entry, 1);
        push_u32(&mut entry, 3);
        push_u64(&mut entry, blob_at as u64);
        file.extend_from_slice(&entry);
        for comp in 0..3 {
            push_f64(&mut file, if comp == 0 { 7.0 } else { 0.0 });
            push_f64(&mut file, if comp == 1 { 1.0 } else { 0.0 });
        }
        file
    };
    let base_index = FileIndex::decode(&base_bytes).expect("base decodes");
    let med_index = FileIndex::decode(&med_bytes).expect("medium decodes");
    let base = CompositeIndex::new(vec![(base_index, &base_bytes[..])]).expect("base composite");
    let med = CompositeIndex::new(vec![(med_index, &med_bytes[..])]).expect("medium composite");
    let tiered =
        TieredIndex::new(vec![(LebTier::Base, base), (LebTier::Medium, med)]).expect("tiers build");
    let tiered_batch = tiered
        .geometric_states(3, &[0.25, 1.75])
        .expect("tiered batch");
    assert_eq!(tiered_batch[0].tier(), Some(LebTier::Base));
    assert_eq!(tiered_batch[1].tier(), Some(LebTier::Medium));
    assert_eq!(tiered_batch[0].pos()[0], 5.0);
    assert_eq!(tiered_batch[1].pos()[0], 7.0);
}

mod center_composition {
    use super::*;
    use crate::constants::{
        FLG_BARYCTR, FLG_HELCTR, FLG_J2000, FLG_SPEED, FLG_TOPOCTR, MEAN_APOG, MEAN_NODE, MOON,
    };

    const TARGET_POS: [f64; 3] = [1.5, -2.5, 0.75];
    const TARGET_VEL: [f64; 3] = [0.01, 0.02, -0.005];
    const SUN_POS: [f64; 3] = [0.1, 0.05, -0.02];
    const SUN_VEL: [f64; 3] = [0.001, -0.002, 0.0005];
    const EARTH_POS: [f64; 3] = [-0.9, 0.3, 0.0001];
    const EARTH_VEL: [f64; 3] = [-0.017, -0.001, 0.0];

    fn compose(center: ObservationCenter) -> ObserverRelativeState {
        compose_relative_state(
            TARGET_POS, TARGET_VEL, SUN_POS, SUN_VEL, EARTH_POS, EARTH_VEL, center, None,
        )
        .expect("finite fixture must compose")
    }

    #[test]
    fn center_flag_priority_mirrors_release() {
        use ObservationCenter::*;
        // Single bits map directly.
        assert_eq!(observation_center(FLG_HELCTR), Heliocentric);
        assert_eq!(observation_center(FLG_BARYCTR), Barycentric);
        assert_eq!(observation_center(FLG_TOPOCTR), Topocentric);
        assert_eq!(observation_center(0), Geocentric);
        assert_eq!(observation_center(FLG_SPEED | FLG_J2000), Geocentric);

        assert_eq!(resolve_center_flags(FLG_HELCTR | FLG_BARYCTR), FLG_BARYCTR);
        assert_eq!(resolve_center_flags(FLG_HELCTR | FLG_TOPOCTR), FLG_TOPOCTR);
        assert_eq!(resolve_center_flags(FLG_BARYCTR | FLG_TOPOCTR), FLG_TOPOCTR);
        assert_eq!(
            resolve_center_flags(FLG_HELCTR | FLG_BARYCTR | FLG_TOPOCTR),
            FLG_TOPOCTR
        );
        // Resolution is idempotent and preserves unrelated bits.
        let mixed = FLG_HELCTR | FLG_BARYCTR | FLG_SPEED;
        let once = resolve_center_flags(mixed);
        assert_eq!(once, FLG_BARYCTR | FLG_SPEED);
        assert_eq!(resolve_center_flags(once), once);
        assert_eq!(observation_center(once), Barycentric);
    }

    #[test]
    fn geocentric_differences_earth_and_routes_its_velocity() {
        let st = compose(ObservationCenter::Geocentric);
        assert_eq!(
            st.geo_pos(),
            [
                TARGET_POS[0] - EARTH_POS[0],
                TARGET_POS[1] - EARTH_POS[1],
                TARGET_POS[2] - EARTH_POS[2],
            ]
        );
        assert_eq!(
            st.geo_vel(),
            [
                TARGET_VEL[0] - EARTH_VEL[0],
                TARGET_VEL[1] - EARTH_VEL[1],
                TARGET_VEL[2] - EARTH_VEL[2],
            ]
        );

        assert_eq!(st.aberration_vel(), EARTH_VEL);
        assert!(!st.is_origin());
    }

    #[test]
    fn heliocentric_differences_sun_and_zeroes_aberration_driver() {
        let st = compose(ObservationCenter::Heliocentric);
        assert_eq!(
            st.geo_pos(),
            [
                TARGET_POS[0] - SUN_POS[0],
                TARGET_POS[1] - SUN_POS[1],
                TARGET_POS[2] - SUN_POS[2],
            ]
        );
        assert_eq!(
            st.geo_vel(),
            [
                TARGET_VEL[0] - SUN_VEL[0],
                TARGET_VEL[1] - SUN_VEL[1],
                TARGET_VEL[2] - SUN_VEL[2],
            ]
        );

        assert_eq!(st.aberration_vel(), [0.0, 0.0, 0.0]);
    }

    #[test]
    fn barycentric_is_the_target_unchanged() {
        let st = compose(ObservationCenter::Barycentric);
        assert_eq!(st.geo_pos(), TARGET_POS);
        assert_eq!(st.geo_vel(), TARGET_VEL);
        assert_eq!(st.aberration_vel(), [0.0, 0.0, 0.0]);
    }

    #[test]
    fn topocentric_adds_the_explicit_offset_verbatim() {
        let topo_pos = [0.0001, -0.0002, 0.00003];
        let topo_vel = [1e-6, -2e-6, 3e-7];
        let st = compose_relative_state(
            TARGET_POS,
            TARGET_VEL,
            SUN_POS,
            SUN_VEL,
            EARTH_POS,
            EARTH_VEL,
            ObservationCenter::Topocentric,
            Some((topo_pos, topo_vel)),
        )
        .expect("finite");
        let obs = [
            EARTH_POS[0] + topo_pos[0],
            EARTH_POS[1] + topo_pos[1],
            EARTH_POS[2] + topo_pos[2],
        ];
        let obs_vel = [
            EARTH_VEL[0] + topo_vel[0],
            EARTH_VEL[1] + topo_vel[1],
            EARTH_VEL[2] + topo_vel[2],
        ];
        assert_eq!(
            st.geo_pos(),
            [
                TARGET_POS[0] - obs[0],
                TARGET_POS[1] - obs[1],
                TARGET_POS[2] - obs[2],
            ]
        );
        assert_eq!(
            st.geo_vel(),
            [
                TARGET_VEL[0] - obs_vel[0],
                TARGET_VEL[1] - obs_vel[1],
                TARGET_VEL[2] - obs_vel[2],
            ]
        );

        assert_eq!(st.aberration_vel(), obs_vel);
        // Without an explicit offset, topocentric == geocentric.
        let bare = compose_relative_state(
            TARGET_POS,
            TARGET_VEL,
            SUN_POS,
            SUN_VEL,
            EARTH_POS,
            EARTH_VEL,
            ObservationCenter::Topocentric,
            None,
        )
        .expect("finite");
        assert_eq!(bare, compose(ObservationCenter::Geocentric));
    }

    #[test]
    fn degenerate_and_sentinel_contracts_yield_exact_zeros() {
        // Sun heliocentric: exact zero in every slot.
        let st = compose_body_geometric(
            SUN, FLG_HELCTR, SUN_POS, SUN_VEL, SUN_POS, SUN_VEL, EARTH_POS, EARTH_VEL, None,
        )
        .expect("finite");
        assert_eq!(st, degenerate_origin());
        assert!(st.is_origin());
        // Earth geocentric: natural zero (same state both sides).
        let st = compose_body_geometric(
            EARTH, 0, EARTH_POS, EARTH_VEL, SUN_POS, SUN_VEL, EARTH_POS, EARTH_VEL, None,
        )
        .expect("finite");
        assert_eq!(st.geo_pos(), [0.0, 0.0, 0.0]);
        assert_eq!(st.geo_vel(), [0.0, 0.0, 0.0]);
        assert!(st.is_origin());
        // Earth topocentric: forced zero (not the -topo offset the
        // raw composition would give).
        let st = compose_body_geometric(
            EARTH,
            FLG_TOPOCTR,
            EARTH_POS,
            EARTH_VEL,
            SUN_POS,
            SUN_VEL,
            EARTH_POS,
            EARTH_VEL,
            Some(([0.0001, 0.0, 0.0], [1e-6, 0.0, 0.0])),
        )
        .expect("finite");
        assert_eq!(st, degenerate_origin());
        // Moon-derived points: zero under helio/bary, served under geo.
        for point in [
            MEAN_NODE, TRUE_NODE, MEAN_APOG, OSCU_APOG, INTP_APOG, INTP_PERG,
        ] {
            assert!(
                is_helio_bary_sentinel_point(point),
                "body {point} must be a sentinel point"
            );
            for flags in [FLG_HELCTR, FLG_BARYCTR] {
                let st = compose_body_geometric(
                    point, flags, TARGET_POS, TARGET_VEL, SUN_POS, SUN_VEL, EARTH_POS, EARTH_VEL,
                    None,
                )
                .expect("finite");
                assert_eq!(st, degenerate_origin(), "body {point} flags {flags}");
            }
            let st = compose_body_geometric(
                point, 0, TARGET_POS, TARGET_VEL, SUN_POS, SUN_VEL, EARTH_POS, EARTH_VEL, None,
            )
            .expect("finite");
            assert!(!st.is_origin());
        }
        // Non-sentinel bodies are never zeroed: Sun geocentric and
        // Earth heliocentric are ordinary differences.
        assert!(!is_helio_bary_sentinel_point(SUN));
        assert!(!is_helio_bary_sentinel_point(MOON));
        assert!(!is_helio_bary_sentinel_point(EARTH));
        assert!(!is_helio_bary_sentinel_point(PLUTO));
        let st = compose_body_geometric(
            SUN, 0, SUN_POS, SUN_VEL, SUN_POS, SUN_VEL, EARTH_POS, EARTH_VEL, None,
        )
        .expect("finite");
        assert!(!st.is_origin());
    }

    #[test]
    fn moon_geocentric_is_the_generic_difference() {
        let moon_pos = [0.0023, -0.0011, 0.0004];
        let moon_vel = [0.012, 0.004, -0.001];
        let st = compose_moon_geocentric(moon_pos, moon_vel, EARTH_POS, EARTH_VEL).expect("finite");
        assert_eq!(
            st.geo_pos(),
            [
                moon_pos[0] - EARTH_POS[0],
                moon_pos[1] - EARTH_POS[1],
                moon_pos[2] - EARTH_POS[2],
            ]
        );
        assert_eq!(st.aberration_vel(), EARTH_VEL);
    }

    #[test]
    fn non_finite_components_are_typed_refusals() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                compose_relative_state(
                    [bad, 0.0, 0.0],
                    TARGET_VEL,
                    SUN_POS,
                    SUN_VEL,
                    EARTH_POS,
                    EARTH_VEL,
                    ObservationCenter::Geocentric,
                    None,
                )
                .unwrap_err(),
                Error::NonFiniteArgument { .. }
            ));
            assert!(matches!(
                compose_relative_state(
                    TARGET_POS,
                    TARGET_VEL,
                    SUN_POS,
                    SUN_VEL,
                    EARTH_POS,
                    [0.0, bad, 0.0],
                    ObservationCenter::Barycentric,
                    None,
                )
                .unwrap_err(),
                Error::NonFiniteArgument { .. }
            ));
        }
    }
}

mod light_time {
    use super::*;
    use crate::constants::FLG_TRUEPOS;
    use std::cell::Cell;

    /// Constant-velocity target: `pos0 + vel*t`, analytic so the
    /// fixed-point residual can be checked against the propagation
    /// equation directly.
    fn linear_target(
        pos0: [f64; 3],
        vel: [f64; 3],
    ) -> impl Fn(f64) -> Result<([f64; 3], [f64; 3])> {
        move |jd| {
            Ok((
                [
                    pos0[0] + vel[0] * jd,
                    pos0[1] + vel[1] * jd,
                    pos0[2] + vel[2] * jd,
                ],
                vel,
            ))
        }
    }

    #[test]
    fn light_speed_constant_pins_release_value() {
        assert_eq!(C_LIGHT_AU_DAY.to_bits(), 0x4065A4A0D4B952CB);
        assert_eq!(
            C_LIGHT_AU_DAY.to_bits(),
            crate::spk::LIGHT_SPEED_AU_PER_DAY.to_bits()
        );
        assert_eq!(LIGHT_TIME_ITERATIONS, 3);
        assert!(wants_light_time(0));
        assert!(!wants_light_time(FLG_TRUEPOS));
        // Unrelated flag bits never disable light-time.
        assert!(wants_light_time(crate::constants::FLG_SPEED));
    }

    #[test]
    fn truepos_reads_once_and_returns_instantaneous_state() {
        let calls = Cell::new(0);
        let target = |jd: f64| {
            calls.set(calls.get() + 1);
            linear_target([1.0, 0.5, -0.25], [0.01, -0.02, 0.005])(jd)
        };
        let obs = [0.1, 0.05, -0.02];
        let obs_vel = [0.001, -0.002, 0.0005];
        let sol = solve_light_time(target, 2451545.0, obs, obs_vel, false).expect("finite");
        assert_eq!(calls.get(), 1);
        assert_eq!(sol.light_time_days(), 0.0);
        assert_eq!(sol.retarded_jd_tt(), 2451545.0);
        assert_eq!(
            sol.geo_pos(),
            [
                1.0 + 0.01 * 2451545.0 - obs[0],
                0.5 - 0.02 * 2451545.0 - obs[1],
                -0.25 + 0.005 * 2451545.0 - obs[2]
            ]
        );
        // No rate factor on the TRUEPOS path: plain difference.
        assert_eq!(
            sol.geo_vel(),
            [0.01 - obs_vel[0], -0.02 - obs_vel[1], 0.005 - obs_vel[2]]
        );
        assert!(!wants_light_time(FLG_TRUEPOS));
    }

    #[test]
    fn three_iterations_satisfy_the_propagation_equation() {
        // Earth-like rate target, observer at the origin: the
        // fixed-point loop must leave |ret-obs| = c*lt up to the
        // third-iterate residual, with exactly 1 + 3 evaluations.
        let calls = Cell::new(0);
        let vel = [0.0172, 0.0, 0.0];
        let target = |jd: f64| {
            calls.set(calls.get() + 1);
            linear_target([0.9, 0.4, 0.1], vel)(jd)
        };
        let sol =
            solve_light_time(target, 1.0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], true).expect("finite");
        assert_eq!(calls.get(), 1 + LIGHT_TIME_ITERATIONS);
        assert!(sol.light_time_days() > 0.0);
        assert_eq!(sol.retarded_jd_tt(), 1.0 - sol.light_time_days());
        // Retarded state matches the analytic line at the retarded epoch.
        let (rp, rv) = sol.retarded_state();
        let expect = linear_target([0.9, 0.4, 0.1], vel)(sol.retarded_jd_tt()).expect("finite");
        assert_eq!(rp, expect.0);
        assert_eq!(rv, expect.1);
        // Propagation residual after the fixed budget: sub-metre
        // (1e-12 AU ≈ 0.15 mm) for inner-system rates.
        let residual = (vec3_dist(sol.geo_pos()) - C_LIGHT_AU_DAY * sol.light_time_days()).abs();
        assert!(residual < 1e-12, "residual {residual:e} too large");
    }

    #[test]
    fn split_loop_settles_on_split_retarded_words() {
        use crate::positions::reference_frames::split_retarded_epoch_tt;
        let calls = Cell::new(0);
        let vel = [
            -0.016830471303132578,
            -0.003289765096742631,
            -0.001433528142177413,
        ];
        let target = |jd: f64| {
            calls.set(calls.get() + 1);
            linear_target(
                [-0.18623373246682284, 0.8829959128210167, 0.383310239881416],
                vel,
            )(jd)
        };
        let whole = 2451545.0;
        let tt_fraction = f64::from_bits(0x3f48_352d_6a6e_1050u64);
        let jd_tt = whole + tt_fraction;
        let obs = [0.9, -0.2, 0.05];
        let sol = solve_light_time_split(target, whole, tt_fraction, jd_tt, obs, [0.0; 3], true)
            .expect("finite");
        assert!(
            calls.get() <= 1 + SPLIT_LIGHT_TIME_ITERATIONS,
            "calls {}",
            calls.get()
        );
        assert!(sol.light_time_days() > 0.0);
        // The retarded epoch is exactly the split-preserving word for the
        // settled light time.
        let expect_jd =
            split_retarded_epoch_tt(whole, tt_fraction, sol.light_time_days()).expect("Sp");
        assert_eq!(sol.retarded_jd_tt().to_bits(), expect_jd.to_bits());
        // Retarded state matches the analytic line at that word.
        let (rp, rv) = sol.retarded_state();
        let expect = linear_target(
            [-0.18623373246682284, 0.8829959128210167, 0.383310239881416],
            vel,
        )(expect_jd)
        .expect("finite");
        assert_eq!(rp, expect.0);
        assert_eq!(rv, expect.1);

        let residual =
            (vec3_dist(sol.geo_pos()) - FRAME_LIGHT_SPEED_AU_DAY * sol.light_time_days()).abs();
        assert!(residual < 1e-12, "residual {residual:e} too large");
        // TRUEPOS path still reads once with zero light time.
        let solo = solve_light_time_split(
            linear_target([1.0, 0.0, 0.0], [0.0; 3]),
            whole,
            tt_fraction,
            jd_tt,
            obs,
            [0.0; 3],
            false,
        );
        assert!(solo.is_ok());
        // Non-finite splits are refused before any provider contact.
        assert!(matches!(
            solve_light_time_split(
                linear_target([1.0, 0.0, 0.0], [0.0; 3]),
                f64::NAN,
                tt_fraction,
                jd_tt,
                obs,
                [0.0; 3],
                true
            )
            .unwrap_err(),
            Error::NonFiniteArgument { .. }
        ));
    }

    #[test]
    fn zero_distance_breaks_before_any_retarded_read() {
        let calls = Cell::new(0);
        let target = |_jd: f64| {
            calls.set(calls.get() + 1);
            Ok(([0.3, -0.1, 0.2], [0.01, 0.02, -0.03]))
        };
        let sol = solve_light_time(target, 2460000.0, [0.3, -0.1, 0.2], [0.0, 0.0, 0.0], true)
            .expect("finite");

        assert_eq!(calls.get(), 1);
        assert_eq!(sol.light_time_days(), 0.0);
        assert_eq!(sol.geo_pos(), [0.0, 0.0, 0.0]);
        assert_eq!(sol.geo_vel(), [0.01, 0.02, -0.03]);
    }

    #[test]
    fn retardation_rate_matches_the_hand_derivative() {
        // Pure radial recession at speed u from distance d: with the
        // observer at rest, lt' = u / (c + u) and the corrected rate
        // is u*(1 - lt').
        let u = 0.02;
        let ret = [u, 0.0, 0.0];
        let geo = [1.3, 0.0, 0.0];
        let lt = 1.3 / C_LIGHT_AU_DAY;
        let ltdot = u / (C_LIGHT_AU_DAY + u);
        let got = light_time_corrected_velocity(ret, [0.0; 3], geo, lt, true);
        assert_eq!(got, [u * (1.0 - ltdot), 0.0, 0.0]);
        // Gates (verbatim): skipped, non-positive lt, or a zero
        // vector all yield the plain difference.
        assert_eq!(
            light_time_corrected_velocity(ret, [0.0; 3], geo, lt, false),
            [u, 0.0, 0.0]
        );
        assert_eq!(
            light_time_corrected_velocity(ret, [0.0; 3], geo, 0.0, true),
            [u, 0.0, 0.0]
        );
        assert_eq!(
            light_time_corrected_velocity(ret, [1.0, 2.0, 3.0], [0.0; 3], lt, true),
            [u - 1.0, -2.0, -3.0]
        );
    }

    #[test]
    fn provider_errors_and_non_finite_inputs_propagate() {
        let ok = |jd: f64| linear_target([1.0, 0.0, 0.0], [0.0, 0.0, 0.0])(jd);
        assert!(matches!(
            solve_light_time(ok, f64::NAN, [0.0; 3], [0.0; 3], true).unwrap_err(),
            Error::NonFiniteArgument { .. }
        ));
        assert!(matches!(
            solve_light_time(ok, 2460000.0, [f64::INFINITY, 0.0, 0.0], [0.0; 3], true).unwrap_err(),
            Error::NonFiniteArgument { .. }
        ));
        // A retarded-epoch provider failure (e.g. sealed range miss)
        // surfaces unchanged, never a substituted vector.
        let failing = |jd: f64| -> Result<([f64; 3], [f64; 3])> {
            if jd < 2460000.0 {
                Err(Error::LebOutOfRange {
                    body_id: 5,
                    jd,
                    start: 2460000.0,
                    end: 2461000.0,
                })
            } else {
                Ok(([2.0, 0.0, 0.0], [0.0, 0.0, 0.0]))
            }
        };
        assert_eq!(
            solve_light_time(failing, 2460000.0, [0.0; 3], [0.0; 3], true).unwrap_err(),
            Error::LebOutOfRange {
                body_id: 5,
                jd: 2460000.0 - 2.0 / C_LIGHT_AU_DAY,
                start: 2460000.0,
                end: 2461000.0,
            }
        );
    }

    #[test]
    fn retarded_provider_read_uses_the_same_composition() {
        let file = fixture(0, COORD_ICRS_BARY);
        let index = FileIndex::decode(&file).expect("fixture must decode");
        let single = SingleLeb::new(&index, &file);
        // Instantaneous distance ~5.2 AU -> lt ~0.03 d; jd 1.0
        // retards to ~0.97, inside [0, 2].
        let lt = vec3_dist([5.0, -1.0, 0.0]) / C_LIGHT_AU_DAY;
        let a = retarded_geometric_state(&single, 0, 1.0, lt, |_| None).expect("served");
        let b = single.geometric_state(0, 1.0 - lt).expect("scalar");
        assert_eq!(a, b);
        assert!(retarded_geometric_state(&single, 0, 0.005, lt, |_| None).is_err());
        assert!(matches!(
            retarded_geometric_state(&single, 0, f64::NAN, lt, |_| None).unwrap_err(),
            Error::NonFiniteArgument { .. }
        ));
        assert!(matches!(
            retarded_geometric_state(&single, 0, 1.0, f64::INFINITY, |_| None).unwrap_err(),
            Error::NonFiniteArgument { .. }
        ));
    }
}

mod aberration_deflection {
    use super::*;
    use crate::constants::{
        FLG_BARYCTR, FLG_HELCTR, FLG_NOABERR, FLG_NOGDEFL, FLG_SPEED, FLG_TRUEPOS, MOON,
    };

    /// All optional deflectors are absent, so the ray remains unchanged.
    fn unserved_lookup(_body_id: i32, _jd_tt: f64) -> Result<[f64; 3]> {
        Err(Error::LebBodyNotFound { body_id: 99 })
    }

    /// Synthetic AU geometry: Sun at the origin, Jupiter and Saturn on the +x axis.
    fn toy_lookup(body_id: i32, _jd_tt: f64) -> Result<[f64; 3]> {
        match body_id {
            SUN => Ok([0.0, 0.0, 0.0]),
            5 => Ok([5.2, 0.0, 0.0]),
            6 => Ok([9.5, 0.0, 0.0]),
            _ => Err(Error::LebBodyNotFound { body_id }),
        }
    }

    #[test]
    fn constants_pin_release_identities() {
        assert_eq!(GRAV_GS, 1.32712440017987e20);
        assert_eq!(SPEED_OF_LIGHT_M_S, 299792458.0);
        assert_eq!(AU_IN_M, 149_597_870_700.0);
        assert_eq!(DEFLECTORS, [(SUN, 1.0), (5, 1047.3486), (6, 3497.898)]);
        assert_eq!(DEFLECTION_EMAG_FLOOR_AU, 0.01);
        assert_eq!(DEFLECTION_EDOTP_LIMIT, 0.99999999999);
        assert_eq!(DEFLECTION_DLIM_SCALE, 1e-6);
        assert_eq!(ABERRATION_R_FLOOR, 1e-30);
    }

    #[test]
    fn gates_mirror_release_branches() {
        // Aberration: every BYPASS bit alone disables it; unrelated
        // bits (NOGDEFL, SPEED) never do.
        assert!(wants_aberration(0));
        assert!(wants_aberration(FLG_SPEED));
        assert!(wants_aberration(FLG_NOGDEFL));
        for flags in [FLG_NOABERR, FLG_HELCTR, FLG_BARYCTR, FLG_TRUEPOS] {
            assert!(!wants_aberration(flags), "flags {flags}");
            assert!(
                !wants_aberration(flags | FLG_SPEED),
                "flags {flags} + SPEED"
            );
        }
        // Deflection: every BYPASS bit alone disables it; the Moon and
        // a non-positive light time disable it; NOABERR never does.
        assert!(wants_deflection(0, 3, 0.005));
        assert!(wants_deflection(FLG_NOABERR, 3, 0.005));
        for flags in [FLG_NOGDEFL, FLG_HELCTR, FLG_BARYCTR, FLG_TRUEPOS] {
            assert!(!wants_deflection(flags, 3, 0.005), "flags {flags}");
        }
        assert!(!wants_deflection(0, MOON, 0.005));
        assert!(!wants_deflection(0, 3, 0.0));
        assert!(!wants_deflection(0, 3, -0.001));
    }

    #[test]
    fn aberration_degenerate_inputs_pass_through() {
        let geo = [1.2, -0.4, 0.3];
        let vel = [0.01, 0.02, -0.005];
        // Exactly-zero vector: no direction to aberrate.
        assert_eq!(
            apply_aberration([0.0, 0.0, 0.0], vel, 0.005).unwrap(),
            [0.0, 0.0, 0.0]
        );
        // Zero Earth velocity on the relativistic path: no aberration.
        assert_eq!(apply_aberration(geo, [0.0, 0.0, 0.0], 0.005).unwrap(), geo);
        // Zero Earth velocity on the Bradley path: identity too.
        assert_eq!(apply_aberration(geo, [0.0, 0.0, 0.0], 0.0).unwrap(), geo);
        // Non-finite inputs are typed refusals, never garbage.
        for bad in [f64::NAN, f64::INFINITY] {
            assert!(matches!(
                apply_aberration([bad, 0.0, 0.0], vel, 0.005).unwrap_err(),
                Error::NonFiniteArgument { .. }
            ));
            assert!(matches!(
                apply_aberration(geo, [0.0, bad, 0.0], 0.005).unwrap_err(),
                Error::NonFiniteArgument { .. }
            ));
            assert!(matches!(
                apply_aberration(geo, vel, bad).unwrap_err(),
                Error::NonFiniteArgument { .. }
            ));
        }
    }

    #[test]
    fn aberration_bradley_matches_hand_derivation() {
        // First-order path (lt = 0): transverse Earth velocity shifts
        // the apparent direction by v/c with the distance preserved.
        let v = 0.0172; // Earth-like, AU/day.
        let geo = [1.0, 0.0, 0.0];
        let out = apply_aberration(geo, [0.0, v, 0.0], 0.0).unwrap();
        // Independent hand derivation: a = u + w - (u.w)u, renormalized.
        let w = v / C_LIGHT_AU_DAY;
        let norm = (1.0 + w * w).sqrt();
        let expect = [1.0 / norm, w / norm, 0.0];
        assert_eq!(out, expect);
        // Distance preserved bit-exactly by the renormalization.
        assert_eq!(vec3_dist(out), vec3_dist(geo));
        // Radial velocity alone leaves a radial vector unchanged.
        let radial = apply_aberration(geo, [0.01, 0.0, 0.0], 0.0).unwrap();
        assert_eq!(radial, geo);
    }

    #[test]
    fn aberration_relativistic_shift_matches_bradley_angle() {
        // Consistent light time (lt = |geo|/c): the annual-aberration
        // angle agrees with v/c to first order (second-order terms are
        // O(beta^2) ~ 1e-8, far inside the 2% sandwich), and the
        // distance is preserved as a pure direction rotation.
        let v = 0.0172;
        let geo = [1.0, 0.0, 0.0];
        let lt = vec3_dist(geo) / C_LIGHT_AU_DAY;
        let out = apply_aberration(geo, [0.0, v, 0.0], lt).unwrap();
        let shift = vec3_dist([out[0] - geo[0], out[1] - geo[1], out[2] - geo[2]]);
        let expect_angle = v / C_LIGHT_AU_DAY;
        assert!(
            (shift - expect_angle).abs() / expect_angle < 0.02,
            "shift {shift:e} vs Bradley {expect_angle:e}"
        );
        let dist_err = (vec3_dist(out) - vec3_dist(geo)).abs() / vec3_dist(geo);
        assert!(dist_err < 1e-7, "distance error {dist_err:e} too large");
        // The shift is transverse: no radial component to first order.
        assert!((out[0] - geo[0]).abs() / shift < 0.02);
    }

    #[test]
    fn deflection_zero_vector_never_touches_the_reader() {
        let out = apply_gravitational_deflection(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            2460000.0,
            0.005,
            |_, _| panic!("no deflector read for a zero vector"),
        )
        .unwrap();
        assert_eq!(out, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn deflection_lookup_outcomes_match_release_policy() {
        let geo = [0.5, 1.0, -0.25];

        let out =
            apply_gravitational_deflection(geo, [1.0, 0.0, 0.0], 2460000.0, 0.005, unserved_lookup)
                .unwrap();
        assert_eq!(out, geo);
        // Same for an SPK range miss.
        let out =
            apply_gravitational_deflection(geo, [1.0, 0.0, 0.0], 2460000.0, 0.005, |_, jd| {
                Err(Error::SpkOutOfRange {
                    target: 0,
                    et: jd,
                    start: 0.0,
                    end: 1.0,
                })
            })
            .unwrap();
        assert_eq!(out, geo);
        // A fatal provisioning failure is NOT a skip: it propagates.
        assert!(matches!(
            apply_gravitational_deflection(geo, [1.0, 0.0, 0.0], 2460000.0, 0.005, |_, _| Err(
                Error::LebCorrupt {
                    what: "boom".to_string(),
                    offset: 0,
                }
            ),)
            .unwrap_err(),
            Error::LebCorrupt { .. }
        ));
        // Non-finite inputs are typed refusals.
        assert!(matches!(
            apply_gravitational_deflection(
                [f64::NAN, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                2460000.0,
                0.005,
                toy_lookup
            )
            .unwrap_err(),
            Error::NonFiniteArgument { .. }
        ));
    }

    #[test]
    fn deflection_skips_own_body_observer() {
        // Planet-centric observer AT the Sun's center: emag = 0 sits
        // below the 0.01 AU floor, so the singular far-field term is
        // skipped and the vector is untouched (Jupiter/Saturn unserved).
        let sun = [0.1, 0.05, -0.02];
        let out =
            apply_gravitational_deflection([1.0, 0.0, 0.0], sun, 2460000.0, 0.005, |id, _| {
                if id == SUN {
                    Ok(sun)
                } else {
                    Err(Error::LebBodyNotFound { body_id: id })
                }
            })
            .unwrap();
        assert_eq!(out, [1.0, 0.0, 0.0]);
    }

    #[test]
    fn deflection_skips_line_of_sight_targets() {
        // Target exactly behind the Sun as seen from [1,0,0]:
        // edotp = -1 trips the 0.99999999999 skip for every served
        // deflector on the axis, leaving the vector untouched.
        let out = apply_gravitational_deflection(
            [-3.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            2460000.0,
            0.02,
            toy_lookup,
        )
        .unwrap();
        assert_eq!(out, [-3.0, 0.0, 0.0]);
    }

    #[test]
    fn deflection_quadrature_order_of_magnitude() {
        // Observer at [1,0,0], Sun at the origin, target at 90 degrees
        // elongation (geo = [0,1,0]): the hand estimate is
        // 2GS/(c^2 * 1 AU)/1.7 ~ 1.2e-8 rad, so sandwich the shift in
        // [1e-9, 1e-7] rad — arcsecond-scale physics, not zero and not
        // arcminutes.
        let geo = [0.0, 1.0, 0.0];
        let out =
            apply_gravitational_deflection(geo, [1.0, 0.0, 0.0], 2460000.0, 0.006, toy_lookup)
                .unwrap();
        let shift = vec3_dist([out[0] - geo[0], out[1] - geo[1], out[2] - geo[2]]);
        assert!(
            (1e-9..1e-7).contains(&shift),
            "quadrature shift {shift:e} outside the physical sandwich"
        );
    }

    #[test]
    fn orchestrator_applies_gated_deflection_then_aberration() {
        let geo = [0.0, 1.0, 0.0];
        let earth_vel = [0.0, 0.0172, 0.0];
        let observer = [1.0, 0.0, 0.0];
        let lt = 0.006;
        // Full path: delta equals the direct deflection increment and
        // the reported vector equals aberration(deflected).
        let full = apply_deflection_and_aberration(
            geo, earth_vel, observer, 2460000.0, lt, 3, 0, toy_lookup, None,
        )
        .unwrap();
        let deflected =
            apply_gravitational_deflection(geo, observer, 2460000.0, lt, toy_lookup).unwrap();
        assert_eq!(
            full.deflection_delta(),
            [
                deflected[0] - geo[0],
                deflected[1] - geo[1],
                deflected[2] - geo[2]
            ]
        );
        assert_ne!(full.deflection_delta(), [0.0, 0.0, 0.0]);
        assert_eq!(
            full.pos(),
            apply_aberration(deflected, earth_vel, lt).unwrap()
        );
        // NOGDEFL: aberration only, zero delta.
        let no_defl = apply_deflection_and_aberration(
            geo,
            earth_vel,
            observer,
            2460000.0,
            lt,
            3,
            FLG_NOGDEFL,
            toy_lookup,
            None,
        )
        .unwrap();
        assert_eq!(no_defl.deflection_delta(), [0.0, 0.0, 0.0]);
        assert_eq!(no_defl.pos(), apply_aberration(geo, earth_vel, lt).unwrap());
        // NOABERR on the Moon: both corrections off, identity.
        let moon = apply_deflection_and_aberration(
            geo,
            earth_vel,
            observer,
            2460000.0,
            lt,
            MOON,
            FLG_NOABERR,
            toy_lookup,
            None,
        )
        .unwrap();
        assert_eq!(moon.pos(), geo);
        assert_eq!(moon.deflection_delta(), [0.0, 0.0, 0.0]);

        let site = [4.0e-5, 0.0, 3.0e-5];
        let moon_fast = apply_deflection_and_aberration(
            geo, earth_vel, observer, 2460000.0, lt, MOON, 0, toy_lookup, None,
        )
        .unwrap();
        assert_eq!(moon_fast.deflection_delta(), [0.0, 0.0, 0.0]);
        assert_eq!(
            moon_fast.pos(),
            apply_aberration(geo, earth_vel, lt).unwrap()
        );
        let moon_vector = apply_deflection_and_aberration(
            geo,
            earth_vel,
            observer,
            2460000.0,
            lt,
            MOON,
            0,
            toy_lookup,
            Some(site),
        )
        .unwrap();
        let deflected =
            apply_gravitational_deflection(geo, observer, 2460000.0, lt, toy_lookup).unwrap();
        let with_earth_manual = apply_aberration(
            apply_earth_deflection(deflected, site).unwrap(),
            earth_vel,
            lt,
        )
        .unwrap();
        assert_eq!(moon_vector.pos(), with_earth_manual);
        assert_ne!(moon_vector.deflection_delta(), [0.0, 0.0, 0.0]);
        // TRUEPOS keeps the geometric Moon vector exactly, site or not.
        let moon_tp = apply_deflection_and_aberration(
            geo,
            earth_vel,
            observer,
            2460000.0,
            lt,
            MOON,
            FLG_TRUEPOS,
            toy_lookup,
            Some(site),
        )
        .unwrap();
        assert_eq!(moon_tp.pos(), geo);
        assert_eq!(moon_tp.deflection_delta(), [0.0, 0.0, 0.0]);
        // TRUEPOS: both corrections off regardless of light time.
        let truepos = apply_deflection_and_aberration(
            geo,
            earth_vel,
            observer,
            2460000.0,
            lt,
            3,
            FLG_TRUEPOS,
            toy_lookup,
            None,
        )
        .unwrap();
        assert_eq!(truepos.pos(), geo);
        // HELCTR/BARYCTR: both corrections off.
        for flags in [FLG_HELCTR, FLG_BARYCTR] {
            let out = apply_deflection_and_aberration(
                geo, earth_vel, observer, 2460000.0, lt, 3, flags, toy_lookup, None,
            )
            .unwrap();
            assert_eq!(out.pos(), geo, "flags {flags}");
        }
    }

    #[test]
    fn orchestrator_applies_earth_term_on_surface_vectors() {
        let geo = [0.0, 1.0, 0.0];
        let earth_vel = [0.0, 0.0172, 0.0];
        let observer = [1.0, 0.0, 0.0];
        let site = [4.0e-5, 0.0, 3.0e-5];
        let lt = 0.006;
        let with_earth = apply_deflection_and_aberration(
            geo,
            earth_vel,
            observer,
            2460000.0,
            lt,
            3,
            0,
            toy_lookup,
            Some(site),
        )
        .unwrap();
        let plain = apply_deflection_and_aberration(
            geo, earth_vel, observer, 2460000.0, lt, 3, 0, toy_lookup, None,
        )
        .unwrap();

        assert_ne!(with_earth.pos(), plain.pos());
        let deflected =
            apply_gravitational_deflection(geo, observer, 2460000.0, lt, toy_lookup).unwrap();
        let with_earth_manual = apply_aberration(
            apply_earth_deflection(deflected, site).unwrap(),
            earth_vel,
            lt,
        )
        .unwrap();
        assert_eq!(with_earth.pos(), with_earth_manual);
        // The frozen delta carries the pre-aberration Earth increment
        // (aberration applies after, identically on both arms).
        let earth_pre = apply_earth_deflection(deflected, site).unwrap();
        let earth_delta = [
            earth_pre[0] - deflected[0],
            earth_pre[1] - deflected[1],
            earth_pre[2] - deflected[2],
        ];
        assert_eq!(
            with_earth.deflection_delta(),
            [
                plain.deflection_delta()[0] + earth_delta[0],
                plain.deflection_delta()[1] + earth_delta[1],
                plain.deflection_delta()[2] + earth_delta[2],
            ]
        );
        // TRUEPOS keeps the geometric vector exactly, site or not.
        let tp = apply_deflection_and_aberration(
            geo,
            earth_vel,
            observer,
            2460000.0,
            lt,
            3,
            FLG_TRUEPOS,
            toy_lookup,
            Some(site),
        )
        .unwrap();
        assert_eq!(tp.pos(), geo);
        assert_eq!(tp.deflection_delta(), [0.0, 0.0, 0.0]);
    }
}

mod routing {
    use super::*;
    use crate::constants::{MEAN_APOG, MEAN_NODE, OSCU_APOG, TRUE_NODE};
    use crate::context::EphemerisContext;

    #[test]
    fn tier_kernel_table_pins_release_filenames() {
        assert_eq!(tier_de_kernel(Tier::Base), "de440s.bsp");
        assert_eq!(tier_de_kernel(Tier::Medium), "de440.bsp");
        assert_eq!(tier_de_kernel(Tier::Extended), "de441.bsp");

        assert_eq!(Tier::default(), Tier::Medium);
        assert_eq!(
            crate::spk::spk_date_range_for_tier(Tier::Medium),
            ("1850-01-01", "2150-01-01")
        );
    }

    #[test]
    fn select_de_kernel_priority_mirrors_release() {
        let tier = Tier::Extended;
        // Full priority: env > explicit > TOML > tier default.
        assert_eq!(
            select_de_kernel(
                Some("explicit.bsp"),
                Some("env.bsp"),
                Some("toml.bsp"),
                tier
            ),
            "env.bsp"
        );
        assert_eq!(
            select_de_kernel(Some("explicit.bsp"), None, Some("toml.bsp"), tier),
            "explicit.bsp"
        );
        assert_eq!(
            select_de_kernel(None, None, Some("toml.bsp"), tier),
            "toml.bsp"
        );
        assert_eq!(select_de_kernel(None, None, None, tier), "de441.bsp");
        // Blank env/TOML spellings are unset, not filenames.
        assert_eq!(
            select_de_kernel(Some("explicit.bsp"), Some("  "), Some(""), tier),
            "explicit.bsp"
        );
        assert_eq!(
            select_de_kernel(None, Some(""), None, Tier::Base),
            "de440s.bsp"
        );
    }

    #[test]
    fn recognized_de_kernel_battery() {
        for good in [
            "de440.bsp",
            "de440s.bsp",
            "DE441.BSP",
            "De421s.Bsp",
            "de4.bsp",
            "/data/kernels/de440.bsp",
        ] {
            assert!(is_recognized_de_kernel(good), "{good} must be recognized");
        }
        for bad in [
            "bogus.bsp",
            "de440.lep",
            "de440bsp",
            "de.bsp",
            "des.bsp",
            "de440s.bsp ",
            "de440.bspx",
            "",
        ] {
            assert!(!is_recognized_de_kernel(bad), "{bad} must be refused");
        }
    }

    #[test]
    fn jpl_file_decision_matrix() {
        assert_eq!(decide_jpl_file(true, false), JplFileDecision::UseRequested);
        assert_eq!(decide_jpl_file(false, true), JplFileDecision::UseRequested);
        assert_eq!(decide_jpl_file(true, true), JplFileDecision::UseRequested);
        // Unrecognized + absent: tier default, never a download of the
        // requested name.
        assert_eq!(decide_jpl_file(false, false), JplFileDecision::TierDefault);
    }

    #[test]
    fn leb_body_classifiers_pin_release_sets() {
        // Core: SUN..=MEAN_APOG plus Earth and the fitted apsides.
        for core in [
            SUN, PLUTO, MEAN_NODE, TRUE_NODE, MEAN_APOG, EARTH, INTP_APOG, INTP_PERG,
        ] {
            assert!(is_leb_core_body(core), "body {core} must be core");
        }
        for non_core in [OSCU_APOG, 15, 22 + 1, -MEAN_NODE] {
            assert!(
                !is_leb_core_body(non_core),
                "body {non_core} must not be core"
            );
        }
        // Stored-channel only: planets plus Earth.
        for stored in [SUN, PLUTO, EARTH] {
            assert!(
                is_leb_stored_state_body(stored),
                "body {stored} must be stored-state"
            );
        }
        for non_stored in [MEAN_NODE, TRUE_NODE, MEAN_APOG, INTP_APOG, INTP_PERG] {
            assert!(
                !is_leb_stored_state_body(non_stored),
                "body {non_stored} must not be stored-state"
            );
        }
        // Model-served: mean node/apogee + fitted apsides.
        for model in [MEAN_NODE, MEAN_APOG, INTP_APOG, INTP_PERG] {
            assert!(
                is_leb_model_served_point(model),
                "body {model} must be model-served"
            );
        }
        for non_model in [SUN, TRUE_NODE, OSCU_APOG, EARTH] {
            assert!(
                !is_leb_model_served_point(non_model),
                "body {non_model} must not be model-served"
            );
        }
    }

    #[test]
    fn route_plans_mirror_release_gating() {
        use ProviderMode::*;
        let full = LegAvailability::new(true, true);
        let bare = LegAvailability::new(false, false);
        let net = EffectiveNetworkPolicy::Allow;
        // Sealed: exactly the LEB leg whatever is provisioned.
        assert_eq!(plan_route(SealedLeb, full, net).legs(), &[BackendLeg::Leb]);
        assert_eq!(plan_route(SealedLeb, bare, net).legs(), &[BackendLeg::Leb]);

        assert_eq!(
            plan_route(NativeJpl, full, net).legs(),
            &[BackendLeg::NativeJpl]
        );
        // Horizons: remote first, then the local fall-through.
        assert_eq!(
            plan_route(HorizonsRemote, full, net).legs(),
            &[BackendLeg::Horizons, BackendLeg::NativeJpl]
        );
        // Auto: LEB when provisioned, Horizons exactly when no local
        // JPL source exists, native JPL always last.
        assert_eq!(
            plan_route(Auto, full, net).legs(),
            &[BackendLeg::Leb, BackendLeg::NativeJpl]
        );
        assert_eq!(
            plan_route(Auto, LegAvailability::new(true, false), net).legs(),
            &[BackendLeg::Leb, BackendLeg::Horizons, BackendLeg::NativeJpl]
        );
        assert_eq!(
            plan_route(Auto, LegAvailability::new(false, true), net).legs(),
            &[BackendLeg::NativeJpl]
        );
        assert_eq!(
            plan_route(Auto, bare, net).legs(),
            &[BackendLeg::Horizons, BackendLeg::NativeJpl]
        );
        // The network gate fires at leg-attempt time, never by
        // omitting legs from the plan — and it travels in the plan.
        let sealed_plan = plan_route(Auto, bare, EffectiveNetworkPolicy::Sealed);
        assert_eq!(
            sealed_plan.legs(),
            &[BackendLeg::Horizons, BackendLeg::NativeJpl]
        );
        assert_eq!(sealed_plan.network_gate(), EffectiveNetworkPolicy::Sealed);
        assert!(LegAvailability::new(true, false).leb_ready());
        assert!(!LegAvailability::new(true, false).local_jpl_ready());
    }

    #[test]
    fn source_tags_name_the_backend_that_ran() {
        assert_eq!(BackendLeg::Leb.as_source_tag(), "LEB");
        assert_eq!(BackendLeg::NativeJpl.as_source_tag(), "NativeJpl");
        assert_eq!(BackendLeg::Horizons.as_source_tag(), "Horizons");
    }

    #[test]
    fn leb_miss_outcome_applies_sealed_policy() {
        // Inside declared coverage: served under every mode.
        for mode in [
            ProviderMode::SealedLeb,
            ProviderMode::Auto,
            ProviderMode::NativeJpl,
            ProviderMode::HorizonsRemote,
        ] {
            assert_eq!(
                leb_miss_outcome(mode, SUN, 1.5, Some((1.0, 2.0))),
                LebMissOutcome::Served
            );
            // Non-sealed modes always fall through outside coverage.
            assert_eq!(
                leb_miss_outcome(mode, SUN, 9.0, Some((1.0, 2.0))),
                if leb_miss_fallthrough(mode) {
                    LebMissOutcome::FallThrough
                } else {
                    LebMissOutcome::Fail(Error::EphemerisRange {
                        requested_jd: 9.0,
                        start_jd: 1.0,
                        end_jd: 2.0,
                        body: Some(SUN),
                    })
                }
            );
        }

        assert_eq!(
            leb_miss_outcome(ProviderMode::SealedLeb, OSCU_APOG, 9.0, Some((1.0, 2.0))),
            LebMissOutcome::FallThrough
        );
        // Sealed core body with no declared coverage but a stored
        // channel: declared missing body.
        assert_eq!(
            leb_miss_outcome(ProviderMode::SealedLeb, SUN, 9.0, None),
            LebMissOutcome::Fail(Error::UnknownBody { body_id: SUN })
        );
        // Sealed core body with a declared non-channel source and no
        // stored window: fall-through (e.g. the true node).
        assert_eq!(
            leb_miss_outcome(ProviderMode::SealedLeb, TRUE_NODE, 9.0, None),
            LebMissOutcome::FallThrough
        );
        assert!(!leb_miss_fallthrough(ProviderMode::SealedLeb));
        assert!(leb_miss_fallthrough(ProviderMode::Auto));
    }

    #[test]
    fn sealed_gates_refuse_with_typed_configuration_errors() {
        assert!(require_leb(true).is_ok());
        assert!(matches!(
            require_leb(false).unwrap_err(),
            Error::ConfigurationMissing { .. }
        ));
        assert!(forbid_jpl_in_sealed_mode(ProviderMode::Auto).is_ok());
        assert!(forbid_jpl_in_sealed_mode(ProviderMode::NativeJpl).is_ok());
        assert!(forbid_jpl_in_sealed_mode(ProviderMode::HorizonsRemote).is_ok());
        assert!(matches!(
            forbid_jpl_in_sealed_mode(ProviderMode::SealedLeb).unwrap_err(),
            Error::ConfigurationMissing { .. }
        ));
    }

    #[test]
    fn routing_reads_the_context_mode_and_stores() {
        // The planner consumes the context's resolved mode plus caller
        // availability: a fresh auto context with nothing bound plans
        // the remote-then-local chain, never a silent sealed answer.
        let ctx = EphemerisContext::new();
        let leg = JplLeg::new(ctx.spk_registry());
        assert!(!leg.has_local_source());
        assert!(ctx.leb_files().is_empty());
        let plan = plan_route(
            ctx.mode(),
            LegAvailability::new(!ctx.leb_files().is_empty(), leg.has_local_source()),
            ctx.effective_network_policy(),
        );
        assert_eq!(plan.legs(), &[BackendLeg::Horizons, BackendLeg::NativeJpl]);
    }

    /// Three-body ICRS fixture plus a stored nutation series (section
    /// 2 over JD 0..16, one 16-day segment, milliradian linear
    /// coefficients): the sealed wire sources nutation from the bound
    /// composite, so a series-free image takes the typed missing-series
    /// refusal instead of serving. Body blobs span JD 0..2.
    fn serving_fixture_nut(states: &[(i32, [f64; 3])]) -> Vec<u8> {
        use crate::leb::SECTION_NUTATION;
        let n = states.len();
        let mut file = pack_header(&MAGIC_LEB1, 1, 2, n as u32, 0.0, 16.0);
        let body_at = HEADER_SIZE + 2 * SECTION_DIR_SIZE;
        let blob_at = body_at + n * LEB1_BODY_ENTRY_SIZE;
        let mut nut = Vec::new();
        push_f64(&mut nut, 0.0);
        push_f64(&mut nut, 16.0);
        push_f64(&mut nut, 16.0);
        push_u32(&mut nut, 1);
        push_u32(&mut nut, 2);
        push_u32(&mut nut, 1);
        push_u32(&mut nut, 0);
        for c in [0.001, 0.0005, 0.002, -0.001] {
            push_f64(&mut nut, c);
        }
        let nut_at = blob_at + n * 96;
        file.extend_from_slice(&pack_section(
            0,
            body_at as u64,
            (n * LEB1_BODY_ENTRY_SIZE) as u64,
        ));
        file.extend_from_slice(&pack_section(
            SECTION_NUTATION,
            nut_at as u64,
            nut.len() as u64,
        ));
        for (i, (body_id, _)) in states.iter().enumerate() {
            file.extend_from_slice(&pack_body(
                *body_id,
                COORD_ICRS_BARY,
                (blob_at + i * 96) as u64,
            ));
        }
        for (_, pos) in states.iter() {
            for _ in 0..2 {
                for component in pos {
                    push_f64(&mut file, *component);
                    push_f64(&mut file, 0.0);
                }
            }
        }
        file.extend_from_slice(&nut);
        file
    }

    fn serving_fixture_with_nutation() -> Vec<u8> {
        serving_fixture_nut(&[
            (SUN, [0.0, 0.0, 0.0]),
            (EARTH, [1.0, 0.0, 0.0]),
            (MARS, [1.5, 0.1, 0.0]),
        ])
    }

    #[test]
    fn serv_wire_sealed_calc_roundtrip() {
        use crate::context::EphemerisContext;
        use crate::leb::LebFile;
        use crate::time::DeltaTModel;
        let leb = LebFile::from_bytes(
            "mem://serv-wire.leb1".to_string(),
            serving_fixture_with_nutation(),
        )
        .expect("fixture decodes");
        let mut ctx = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        ctx.attach_leb_file(leb);
        let (pos, retflag) =
            with_sealed_calc_engine(&ctx, DeltaTModel::EspenakMeeus, false, |engine| {
                calc(engine, 1.0, MARS, FLG_SWIEPH)
            })
            .expect("sealed wire serves");
        for v in pos {
            assert!(v.is_finite(), "finite core: {pos:?}");
        }
        assert!(pos[2] > 0.0, "positive distance: {pos:?}");
        assert!((0.0..360.0).contains(&pos[0]), "lon wrapped: {pos:?}");
        assert_ne!(retflag & FLG_SWIEPH, 0, "ephemeris echo: {retflag:#x}");
    }

    #[test]
    fn serv_wire_refuses_empty_store() {
        use crate::context::EphemerisContext;
        use crate::time::DeltaTModel;
        let ctx = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        let err = with_sealed_calc_engine(&ctx, DeltaTModel::EspenakMeeus, false, |_| Ok(()))
            .expect_err("no bound image must refuse");
        assert!(
            matches!(
                err,
                Error::ConfigurationMissing {
                    missing_config: "leb_reader",
                    ..
                }
            ),
            "unexpected refusal: {err:?}"
        );
    }

    #[test]
    fn serv_wire_refuses_non_leb_modes() {
        use crate::context::EphemerisContext;
        use crate::leb::LebFile;
        use crate::time::DeltaTModel;
        for mode in [ProviderMode::NativeJpl, ProviderMode::HorizonsRemote] {
            let leb = LebFile::from_bytes("mem://serv-wire.leb1".to_string(), serving_fixture())
                .expect("fixture decodes");
            let mut ctx = EphemerisContext::new().with_mode(mode);
            ctx.attach_leb_file(leb);
            let err = with_sealed_calc_engine(&ctx, DeltaTModel::EspenakMeeus, false, |_| Ok(()))
                .expect_err("non-LEB leg must refuse");
            assert!(
                matches!(
                    err,
                    Error::ConfigurationMissing {
                        missing_config: "non_leb_leg_wiring",
                        ..
                    }
                ),
                "unexpected refusal for {mode:?}: {err:?}"
            );
        }
    }

    #[test]
    fn serv_wire_auto_with_leb_serves() {
        use crate::context::EphemerisContext;
        use crate::leb::LebFile;
        use crate::time::DeltaTModel;
        let leb = LebFile::from_bytes(
            "mem://serv-wire.leb1".to_string(),
            serving_fixture_with_nutation(),
        )
        .expect("fixture decodes");
        let mut ctx = EphemerisContext::new();
        assert_eq!(ctx.mode(), ProviderMode::Auto);
        ctx.attach_leb_file(leb);
        let (pos, _) = with_sealed_calc_engine(&ctx, DeltaTModel::EspenakMeeus, false, |engine| {
            calc(engine, 1.0, MARS, FLG_SWIEPH | FLG_TRUEPOS)
        })
        .expect("auto with bound LEB serves");
        for v in pos {
            assert!(v.is_finite(), "finite core: {pos:?}");
        }
    }

    #[test]
    fn native_ecl_nut_ignores_body_model_override_and_stored_series() {
        use crate::context::EphemerisContext;
        use crate::leb::LebFile;
        use crate::time::DeltaTModel;
        let leb = LebFile::from_bytes(
            "mem://nutation-routing.leb1".to_string(),
            serving_fixture_with_nutation(),
        )
        .unwrap();
        let mut ctx = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        ctx.attach_leb_file(leb);
        let jd_ut = 2460232.249616082;
        let expected = ecl_nut_result_ut(jd_ut, FLG_SWIEPH).unwrap();
        for override_days in [None, Some(0.25)] {
            ctx.set_user_delta_t(override_days);
            let actual =
                with_sealed_calc_engine(&ctx, DeltaTModel::EspenakMeeus, false, |engine| {
                    calc_ut(engine, jd_ut, ECL_NUT, FLG_SWIEPH)
                })
                .unwrap();
            assert_eq!(actual.0.map(f64::to_bits), expected.0.map(f64::to_bits));
            assert_eq!(actual.1, expected.1);
        }
        // Explicit TT uses the native series too, irrespective of the
        // stored frame channel and ephemeris coverage of the fixture.
        let expected_tt = ecl_nut_result(jd_ut, FLG_SWIEPH, &crate::erfa::nut06a).unwrap();
        let actual_tt = with_sealed_calc_engine(&ctx, DeltaTModel::EspenakMeeus, false, |engine| {
            calc(engine, jd_ut, ECL_NUT, FLG_SWIEPH)
        })
        .unwrap();
        assert_eq!(
            actual_tt.0.map(f64::to_bits),
            expected_tt.0.map(f64::to_bits)
        );
    }

    #[test]
    fn serv1_stack_engine_matches_sealed_engine_bitwise() {
        use crate::context::EphemerisContext;
        use crate::leb::LebFile;
        use crate::time::DeltaTModel;
        let leb = LebFile::from_bytes(
            "mem://serv1.leb1".to_string(),
            serving_fixture_with_nutation(),
        )
        .expect("fixture decodes");
        let mut ctx = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        ctx.attach_leb_file(leb);
        let via_engine =
            with_sealed_calc_engine(&ctx, DeltaTModel::EspenakMeeus, false, |engine| {
                calc(engine, 1.0, MARS, FLG_SWIEPH)
            })
            .expect("sealed engine serves");
        let via_stack =
            with_sealed_serving_stack(&ctx, DeltaTModel::EspenakMeeus, false, false, |stack| {
                let ut1_to_tt_jd =
                    |tjdut: f64| Ok(tjdut + deltat_with_config(tjdut, stack.delta_cfg)?);
                let engine = CalcEngine {
                    ut1_to_tt_jd: &ut1_to_tt_jd,
                    compute_tt: &|jd, body, calc_flags| {
                        serving_compute_tt_escalated(
                            stack.mode,
                            stack.serving,
                            jd,
                            body,
                            calc_flags,
                        )
                    },
                    nutation_obliquity_ut: &|tjdut, nut_flags| {
                        let jd_tt = ut1_to_tt_jd(tjdut)?;
                        ecl_nut_result(jd_tt, nut_flags, stack.nutation_at)
                    },
                    nutation_obliquity_tt: &|jd_tt, nut_flags| {
                        ecl_nut_result(jd_tt, nut_flags, stack.nutation_at)
                    },
                    jpl_fallback_active: false,
                };
                calc(&engine, 1.0, MARS, FLG_SWIEPH)
            })
            .expect("stack engine serves");
        assert_eq!(via_engine.1, via_stack.1, "retflag words");
        assert_eq!(via_engine.0.len(), via_stack.0.len());
        for (a, b) in via_engine.0.iter().zip(via_stack.0.iter()) {
            assert_eq!(a.to_bits(), b.to_bits(), "bitwise position words");
        }
    }

    #[test]
    fn serv1_default_smh2016_serves_via_vendored_table() {
        use crate::context::EphemerisContext;
        use crate::leb::LebFile;
        use crate::time::DeltaTModel;
        let leb = LebFile::from_bytes(
            "mem://serv1.leb1".to_string(),
            serving_fixture_with_nutation(),
        )
        .expect("fixture decodes");
        let mut ctx = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        ctx.attach_leb_file(leb);
        let days = with_sealed_serving_stack(&ctx, DeltaTModel::Smh2016, false, false, |stack| {
            deltat_with_config(2451545.0, stack.delta_cfg)
        })
        .expect("default smh2016 serves through the vendored table");

        assert_eq!(days.to_bits(), (63.82896039639765 / 86400.0f64).to_bits());
    }

    #[test]
    fn serv1_shared_ctor_matches_legacy_engine_bitwise() {
        use crate::context::EphemerisContext;
        use crate::leb::LebFile;
        use crate::time::DeltaTModel;
        let leb = LebFile::from_bytes(
            "mem://serv1.leb1".to_string(),
            serving_fixture_with_nutation(),
        )
        .expect("fixture decodes");
        let mut ctx = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        ctx.attach_leb_file(leb);
        let via_engine =
            with_sealed_calc_engine(&ctx, DeltaTModel::EspenakMeeus, false, |engine| {
                calc_ut(engine, 0.2, MARS, FLG_SWIEPH)
            })
            .expect("legacy engine serves");
        let via_shared = with_shared_calc_engine(
            ServingRequest {
                ctx: &ctx,
                delta_t_model: DeltaTModel::EspenakMeeus,
                iers_enabled: false,
                observer: None,
                tier_policy: TierPolicy::AttachOrder,
                skip_closed: false,
            },
            |engine| calc_ut(engine, 0.2, MARS, FLG_SWIEPH),
        )
        .expect("shared constructor serves");
        assert_eq!(via_engine.1, via_shared.1, "retflag words");
        for (a, b) in via_engine.0.iter().zip(via_shared.0.iter()) {
            assert_eq!(a.to_bits(), b.to_bits(), "bitwise position words");
        }
    }

    #[test]
    fn serv1_tier_priority_orders_by_tier_token() {
        use crate::context::EphemerisContext;
        use crate::leb::LebFile;
        use crate::time::DeltaTModel;
        let medium_img = serving_fixture_nut(&[
            (SUN, [0.0, 0.0, 0.0]),
            (EARTH, [1.0, 0.0, 0.0]),
            (MARS, [7.0, 0.0, 0.0]),
        ]);
        let base_img = serving_fixture_nut(&[
            (SUN, [0.0, 0.0, 0.0]),
            (EARTH, [1.0, 0.0, 0.0]),
            (MARS, [1.5, 0.1, 0.0]),
        ]);
        // Baselines: each image alone.
        let solo = |label: &str, bytes: Vec<u8>| {
            let leb = LebFile::from_bytes(label.to_string(), bytes).expect("decodes");
            let mut ctx = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
            ctx.attach_leb_file(leb);
            with_shared_calc_engine(
                ServingRequest {
                    ctx: &ctx,
                    delta_t_model: DeltaTModel::EspenakMeeus,
                    iers_enabled: false,
                    observer: None,
                    tier_policy: TierPolicy::AttachOrder,
                    skip_closed: false,
                },
                |engine| calc_ut(engine, 0.2, MARS, FLG_SWIEPH),
            )
            .expect("solo serves")
            .0[0]
                .to_bits()
        };
        let medium_lon = solo("mem://medium_a.leb1", medium_img.clone());
        let base_lon = solo("mem://base_b.leb1", base_img.clone());
        assert_ne!(medium_lon, base_lon, "fixtures must discriminate");
        // Attach medium first: attach order serves medium, tier
        // priority serves base (`LebTier::PRIORITY` puts Base first).
        let mut ctx = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        ctx.attach_leb_file(
            LebFile::from_bytes("mem://medium_a.leb1".to_string(), medium_img).expect("decodes"),
        );
        ctx.attach_leb_file(
            LebFile::from_bytes("mem://base_b.leb1".to_string(), base_img).expect("decodes"),
        );
        let run = |policy| {
            with_shared_calc_engine(
                ServingRequest {
                    ctx: &ctx,
                    delta_t_model: DeltaTModel::EspenakMeeus,
                    iers_enabled: false,
                    observer: None,
                    tier_policy: policy,
                    skip_closed: false,
                },
                |engine| calc_ut(engine, 0.2, MARS, FLG_SWIEPH),
            )
            .expect("serves")
            .0[0]
                .to_bits()
        };
        assert_eq!(run(TierPolicy::AttachOrder), medium_lon);
        assert_eq!(run(TierPolicy::TierPriority), base_lon);
    }

    #[test]
    fn serv1_tier_priority_refuses_untokened_and_wide_sets() {
        use crate::context::EphemerisContext;
        use crate::context::Tier;
        use crate::leb::LebFile;
        use crate::time::DeltaTModel;
        // No tier token on the label: tier priority cannot be
        // established (grouped dispatch requires recognized tier labels).
        let leb = LebFile::from_bytes(
            "mem://serv1.leb1".to_string(),
            serving_fixture_with_nutation(),
        )
        .expect("fixture decodes");
        let mut ctx = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        ctx.attach_leb_file(leb);
        let err = with_shared_calc_engine(
            ServingRequest {
                ctx: &ctx,
                delta_t_model: DeltaTModel::EspenakMeeus,
                iers_enabled: false,
                observer: None,
                tier_policy: TierPolicy::TierPriority,
                skip_closed: false,
            },
            |_| Ok(()),
        )
        .expect_err("untokened labels refuse under tier priority");
        assert!(
            matches!(err, Error::ConfigurationMissing { .. }),
            "unexpected refusal: {err:?}"
        );
        // Medium image under a Base active tier: provisioning error,
        // never silent widening.
        let leb = LebFile::from_bytes(
            "mem://medium_a.leb1".to_string(),
            serving_fixture_with_nutation(),
        )
        .expect("fixture decodes");
        let mut narrow = EphemerisContext::new()
            .with_mode(ProviderMode::SealedLeb)
            .with_tier(Tier::Base);
        narrow.attach_leb_file(leb);
        let err = with_shared_calc_engine(
            ServingRequest {
                ctx: &narrow,
                delta_t_model: DeltaTModel::EspenakMeeus,
                iers_enabled: false,
                observer: None,
                tier_policy: TierPolicy::TierPriority,
                skip_closed: false,
            },
            |_| Ok(()),
        )
        .expect_err("wider installed set refuses");
        assert!(
            matches!(err, Error::LebCorrupt { .. }),
            "unexpected refusal: {err:?}"
        );
    }

    #[test]
    fn serv1_topo_binds_explicit_then_stored_observer() {
        use crate::context::{EphemerisContext, Observer};
        use crate::leb::LebFile;
        use crate::time::DeltaTModel;
        let leb = LebFile::from_bytes(
            "mem://serv1.leb1".to_string(),
            serving_fixture_with_nutation(),
        )
        .expect("fixture decodes");
        let mut ctx = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        ctx.attach_leb_file(leb);
        let topo = FLG_SWIEPH | FLG_TOPOCTR;
        let run = |ctx: &EphemerisContext, observer: Option<Observer>| {
            with_shared_calc_engine(
                ServingRequest {
                    ctx,
                    delta_t_model: DeltaTModel::EspenakMeeus,
                    iers_enabled: false,
                    observer,
                    tier_policy: TierPolicy::AttachOrder,
                    skip_closed: false,
                },
                |engine| calc_ut(engine, 0.2, MARS, topo),
            )
        };
        // No observer anywhere: the missing-observer refusal stays.
        let err = run(&ctx, None).expect_err("topo without observer refuses");
        assert!(
            matches!(err, Error::ConfigurationMissing { .. }),
            "unexpected refusal: {err:?}"
        );
        // Explicit observer serves finite words.
        let site_a = Observer::new(11.0, 45.0, 100.0).expect("valid site");
        let (pos_a, _) = run(&ctx, Some(site_a)).expect("explicit site serves");
        for v in &pos_a {
            assert!(v.is_finite(), "finite topo words: {pos_a:?}");
        }
        // A different explicit site answers different words.
        let site_b = Observer::new(-73.99, 40.7, 10.0).expect("valid site");
        let (pos_b, _) = run(&ctx, Some(site_b)).expect("second site serves");
        assert_ne!(pos_a[0].to_bits(), pos_b[0].to_bits());
        // The context-stored observer serves when no explicit one is
        // given, and the explicit one wins over it.
        let mut ctx_stored = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        ctx_stored.attach_leb_file(
            LebFile::from_bytes(
                "mem://serv1.leb1".to_string(),
                serving_fixture_with_nutation(),
            )
            .expect("fixture decodes"),
        );
        ctx_stored = ctx_stored.with_observer(site_b);
        let (pos_stored, _) = run(&ctx_stored, None).expect("stored observer serves");
        for v in &pos_stored {
            assert!(v.is_finite(), "finite topo words: {pos_stored:?}");
        }
        let (pos_over, _) = run(&ctx_stored, Some(site_a)).expect("explicit wins");
        for (a, b) in pos_over.iter().zip(pos_a.iter()) {
            assert_eq!(a.to_bits(), b.to_bits(), "explicit beats stored");
        }
    }

    #[test]
    fn serv1_topo_matches_explicit_site_path_bitwise() {
        use crate::context::{EphemerisContext, Observer};
        use crate::leb::LebFile;
        use crate::time::{DeltaTModel, deltat_with_config};
        let leb = LebFile::from_bytes(
            "mem://serv1.leb1".to_string(),
            serving_fixture_with_nutation(),
        )
        .expect("fixture decodes");
        let mut ctx = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        ctx.attach_leb_file(leb);
        let site = Observer::new(11.0, 45.0, 100.0).expect("valid site");
        let topo = FLG_SWIEPH | FLG_TOPOCTR;
        let (pos, retflag) = with_shared_calc_engine(
            ServingRequest {
                ctx: &ctx,
                delta_t_model: DeltaTModel::EspenakMeeus,
                iers_enabled: false,
                observer: Some(site),
                tier_policy: TierPolicy::AttachOrder,
                skip_closed: false,
            },
            |engine| calc_ut(engine, 0.2, MARS, topo),
        )
        .expect("topo serves");

        let cfg = ctx.delta_t_config(DeltaTModel::EspenakMeeus, false);
        let jd_tt = 0.2 + deltat_with_config(0.2, &cfg).expect("dt serves");
        let partner = jd_tt - deltat_with_config(jd_tt, &cfg).expect("dt serves");
        let direct =
            with_sealed_serving_stack(&ctx, DeltaTModel::EspenakMeeus, false, false, |stack| {
                let c2t_at = |tt: f64, ut1: f64| crate::erfa::c2t06a(tt, ut1, 0.0, 0.0);
                serving_compute_tt_topocentric(
                    stack.serving,
                    jd_tt,
                    partner,
                    MARS,
                    topo,
                    (site.lon_deg(), site.lat_deg(), site.height_m()),
                    &c2t_at,
                    false,
                    SiteEllipsoid::Wgs84,
                )
            })
            .expect("direct site path serves");
        for (a, b) in pos.iter().zip(direct.0.iter()) {
            assert_eq!(a.to_bits(), b.to_bits(), "bitwise topo words");
        }
        assert_eq!(direct.1, retflag, "retflag words");
    }

    #[test]
    fn serv1_native_jpl_does_not_substitute_bound_leb() {
        use crate::context::EphemerisContext;
        use crate::leb::LebFile;
        use crate::time::DeltaTModel;
        let leb = LebFile::from_bytes(
            "mem://serv1.leb1".to_string(),
            serving_fixture_with_nutation(),
        )
        .expect("fixture decodes");
        let mut ctx = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
        ctx.attach_leb_file(leb);
        let err = with_shared_calc_engine(
            ServingRequest {
                ctx: &ctx,
                delta_t_model: DeltaTModel::EspenakMeeus,
                iers_enabled: false,
                observer: None,
                tier_policy: TierPolicy::AttachOrder,
                skip_closed: false,
            },
            |engine| calc(engine, 2451545.0, MOON, FLG_SWIEPH),
        )
        .expect_err("NativeJpl dispatch requires its registered target despite a bound LEB");
        assert!(
            matches!(err, Error::SpkUnregisteredBody { body: MOON }),
            "unexpected refusal: {err:?}"
        );
    }

    #[test]
    fn serv1_closed_file_policy_legs() {
        use crate::context::EphemerisContext;
        use crate::leb::LebFile;
        use crate::time::DeltaTModel;
        let mut closed =
            LebFile::from_bytes("mem://serv1-closed.leb1".to_string(), serving_fixture())
                .expect("fixture decodes");
        closed.close();
        let mut ctx = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        ctx.attach_leb_file(closed);
        let err =
            with_sealed_serving_stack(&ctx, DeltaTModel::EspenakMeeus, false, false, |_| Ok(()))
                .expect_err("fail-fast must surface the closed reader");
        assert!(
            matches!(err, Error::LebReaderClosed),
            "unexpected refusal: {err:?}"
        );
        // Skip-closed counts the image as absent: the empty remainder
        // refuses via the composite (the bindings sealed-path arm).
        let err =
            with_sealed_serving_stack(&ctx, DeltaTModel::EspenakMeeus, false, true, |_| Ok(()))
                .expect_err("empty remainder must refuse");
        assert!(
            matches!(err, Error::LebCorrupt { .. }),
            "unexpected refusal: {err:?}"
        );
        // Skip-closed with a live remainder serves over it.
        let live = LebFile::from_bytes(
            "mem://serv1-live.leb1".to_string(),
            serving_fixture_with_nutation(),
        )
        .expect("fixture decodes");
        ctx.attach_leb_file(live);
        with_sealed_serving_stack(&ctx, DeltaTModel::EspenakMeeus, false, true, |stack| {
            assert_eq!(stack.composite.len(), 1, "only the live reader remains");
            let ut1_to_tt_jd = |tjdut: f64| Ok(tjdut + deltat_with_config(tjdut, stack.delta_cfg)?);
            let engine = CalcEngine {
                ut1_to_tt_jd: &ut1_to_tt_jd,
                compute_tt: &|jd, body, calc_flags| {
                    serving_compute_tt_escalated(stack.mode, stack.serving, jd, body, calc_flags)
                },
                nutation_obliquity_ut: &|tjdut, nut_flags| {
                    let jd_tt = ut1_to_tt_jd(tjdut)?;
                    ecl_nut_result(jd_tt, nut_flags, stack.nutation_at)
                },
                nutation_obliquity_tt: &|jd_tt, nut_flags| {
                    ecl_nut_result(jd_tt, nut_flags, stack.nutation_at)
                },
                jpl_fallback_active: false,
            };
            let (pos, _) = calc(&engine, 1.0, MARS, FLG_SWIEPH)?;
            for v in pos {
                assert!(v.is_finite(), "finite core: {pos:?}");
            }
            Ok(())
        })
        .expect("live remainder serves");
    }

    #[test]
    fn dt1doc_calc_ut_configured_subset_native() {
        use crate::context::EphemerisContext;
        use crate::leb::LebFile;
        use crate::time::DeltaTModel;
        // The hand fixture spans TT 0..2; the polynomial ΔT at these
        // ancient epochs is about +1.58 d, so the UT epoch 0.2 converts
        // to TT ≈ 1.78, inside coverage.
        let leb = LebFile::from_bytes(
            "mem://dt1doc.leb1".to_string(),
            serving_fixture_with_nutation(),
        )
        .expect("fixture decodes");
        let mut ctx = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        ctx.attach_leb_file(leb);
        let (pos, _) = with_sealed_calc_engine(&ctx, DeltaTModel::EspenakMeeus, false, |engine| {
            calc_ut(engine, 0.2, MARS, FLG_SWIEPH)
        })
        .expect("espenak_meeus UT serves");
        for v in pos {
            assert!(v.is_finite(), "finite core: {pos:?}");
        }
        // The default model serves through the vendored SMH-2016 table
        // (JD 0.2 lands on the long-term curve at about +1.60 d, so TT
        // ≈ 1.80 stays inside coverage, as with the polynomial leg).
        // The model travels the fully-default resolve path (no
        // programmatic/env/TOML selection), proving default-context
        // `calc_ut` end to end.
        use crate::context::{resolve_delta_t_model, resolve_iers_enabled};
        let model = resolve_delta_t_model(None, None, None);
        assert_eq!(model, DeltaTModel::Smh2016);
        let iers = resolve_iers_enabled(None, None, None);
        assert!(!iers);
        let (pos, _) = with_sealed_calc_engine(&ctx, model, iers, |engine| {
            calc_ut(engine, 0.2, MARS, FLG_SWIEPH)
        })
        .expect("default smh2016 UT serves");
        for v in pos {
            assert!(v.is_finite(), "finite core: {pos:?}");
        }
        // A user override wins verbatim even under the default model
        // (TT 1.0 + 0.0025 stays inside coverage).
        let leb = LebFile::from_bytes(
            "mem://dt1doc-override.leb1".to_string(),
            serving_fixture_with_nutation(),
        )
        .expect("fixture decodes");
        let mut ctx = EphemerisContext::new()
            .with_mode(ProviderMode::SealedLeb)
            .with_user_delta_t(Some(0.0025));
        ctx.attach_leb_file(leb);
        let (pos, _) = with_sealed_calc_engine(&ctx, DeltaTModel::Smh2016, false, |engine| {
            calc_ut(engine, 1.0, MARS, FLG_SWIEPH)
        })
        .expect("override UT serves");
        for v in pos {
            assert!(v.is_finite(), "finite core: {pos:?}");
        }
    }

    #[test]
    fn serv_wire_user_delta_t_shifts_ut() {
        use crate::context::EphemerisContext;
        use crate::leb::LebFile;
        use crate::time::DeltaTModel;
        let leb = LebFile::from_bytes("mem://serv-wire.leb1".to_string(), serving_fixture())
            .expect("fixture decodes");
        let mut ctx = EphemerisContext::new()
            .with_mode(ProviderMode::SealedLeb)
            .with_user_delta_t(Some(0.0025));
        ctx.attach_leb_file(leb);
        with_sealed_calc_engine(&ctx, DeltaTModel::EspenakMeeus, false, |engine| {
            let jd_tt = (engine.ut1_to_tt_jd)(2451545.0)?;
            assert_eq!(jd_tt, 2451545.0 + 0.0025, "verbatim override shift");
            Ok(())
        })
        .expect("override chain serves");
    }

    #[test]
    fn serv_wire_closed_file_refuses_typed() {
        use crate::context::EphemerisContext;
        use crate::leb::LebFile;
        use crate::time::DeltaTModel;
        let mut leb = LebFile::from_bytes("mem://serv-wire.leb1".to_string(), serving_fixture())
            .expect("fixture decodes");
        leb.close();
        let mut ctx = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        ctx.attach_leb_file(leb);
        let err = with_sealed_calc_engine(&ctx, DeltaTModel::EspenakMeeus, false, |_| Ok(()))
            .expect_err("closed image must refuse");
        assert!(
            matches!(err, Error::LebReaderClosed),
            "unexpected refusal: {err:?}"
        );
    }

    #[test]
    fn serv_wire_ecl_nut_without_stored_series_uses_native_series() {
        use crate::context::EphemerisContext;
        use crate::leb::LebFile;
        use crate::time::DeltaTModel;
        let leb = LebFile::from_bytes("mem://serv-wire.leb1".to_string(), serving_fixture())
            .expect("fixture decodes");
        let mut ctx = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        ctx.attach_leb_file(leb);
        let actual = with_sealed_calc_engine(&ctx, DeltaTModel::EspenakMeeus, false, |engine| {
            calc(engine, 1.0, ECL_NUT, FLG_SWIEPH)
        })
        .expect("ECL_NUT does not require stored nutation");
        let expected = ecl_nut_result(1.0, FLG_SWIEPH, &crate::erfa::nut06a).unwrap();
        assert_eq!(actual.0.map(f64::to_bits), expected.0.map(f64::to_bits));
        assert_eq!(actual.1, expected.1);
    }

    /// Minimal LE type-2 kernel: record 1 = file, record 2 = one
    /// summary (target 4, center 0, J2000, type 2 over et
    /// `[0, 86400]`), record 3 = one degree-0 record plus the
    /// trailing `[INIT, INTLEN, RSIZE, N]` directory.
    fn tiny_kernel_image() -> Vec<u8> {
        use crate::spk::{RECORD_BYTES, SPK_FILE_ID};
        let mut img = vec![0u8; 3 * RECORD_BYTES];
        img[0..8].copy_from_slice(SPK_FILE_ID);
        img[8..12].copy_from_slice(&2i32.to_le_bytes());
        img[12..16].copy_from_slice(&6i32.to_le_bytes());
        img[76..80].copy_from_slice(&2i32.to_le_bytes());
        img[80..84].copy_from_slice(&2i32.to_le_bytes());
        img[84..88].copy_from_slice(&385i32.to_le_bytes());
        img[88..96].copy_from_slice(b"LTL-IEEE");
        let put_f64 = |img: &mut [u8], off: usize, v: f64| {
            img[off..off + 8].copy_from_slice(&v.to_le_bytes());
        };
        let put_i32 = |img: &mut [u8], off: usize, v: i32| {
            img[off..off + 4].copy_from_slice(&v.to_le_bytes());
        };
        // Summary record: next/prev/nsum, then one 5-word summary.
        let s = RECORD_BYTES;
        put_f64(&mut img, s, 0.0);
        put_f64(&mut img, s + 8, 0.0);
        put_f64(&mut img, s + 16, 1.0);
        put_f64(&mut img, s + 24, 0.0);
        put_f64(&mut img, s + 32, 86400.0);
        // Integer slots in file order: target, center, frame, type,
        // initial address, final address (1-based words).
        put_i32(&mut img, s + 40, 4);
        put_i32(&mut img, s + 44, 0);
        put_i32(&mut img, s + 48, 1);
        put_i32(&mut img, s + 52, 2);
        put_i32(&mut img, s + 56, 257);
        put_i32(&mut img, s + 60, 265);
        // Element record: [mid, radius, cx, cy, cz, INIT, INTLEN,
        // RSIZE, N].
        let e = 2 * RECORD_BYTES;
        for (k, v) in [43200.0, 43200.0, 1.5, -2.5, 0.75, 0.0, 86400.0, 5.0, 1.0]
            .into_iter()
            .enumerate()
        {
            put_f64(&mut img, e + k * 8, v);
        }
        img
    }

    #[test]
    fn jpl_leg_serves_a_registered_chain_geometrically() {
        let registry = SpkRegistry::new();
        let leg = JplLeg::new(&registry);
        // Empty store: no source, no coverage, typed absence.
        assert!(!leg.has_local_source());
        assert!(!leg.covers_tdb(4, 2451545.0));
        assert!(!leg.covers_tdb(4, f64::NAN));
        assert!(matches!(
            leg.eval_geometric_ssb(4, 43200.0).unwrap_err(),
            Error::SpkUnregisteredBody { .. }
        ));
        // Register the tiny kernel: the leg answers geometrically.
        let mut registry = SpkRegistry::new();
        registry
            .register_image("tiny.bsp", tiny_kernel_image())
            .expect("fixture must decode");
        registry
            .register_body(4, "tiny.bsp", 4)
            .expect("chain must register");
        let leg = JplLeg::new(&registry);
        assert!(leg.has_local_source());
        let spans = registry.coverage_jd_tdb(4).expect("spans");
        assert_eq!(spans.len(), 1);
        let mid_jd = (spans[0].0 + spans[0].1) / 2.0;
        assert!(leg.covers_tdb(4, mid_jd));
        assert!(!leg.covers_tdb(4, spans[0].1 + 10.0));
        assert!(!leg.covers_tdb(99, mid_jd));
        let state = leg.eval_geometric_ssb(4, 43200.0).expect("covered");
        // Degree-0 type 2: constant position in km, zero velocity.
        assert_eq!(state.position_km, [1.5, -2.5, 0.75]);
        assert_eq!(state.velocity_km_s, [0.0, 0.0, 0.0]);
        // Outside stored coverage: the evaluator's typed range miss,
        // never invented data.
        assert!(matches!(
            leg.eval_geometric_ssb(4, 200_000.0).unwrap_err(),
            Error::SpkOutOfRange { .. }
        ));
    }

    #[test]
    fn minor_strict_gate_refuses_before_keplerian_fallback() {
        use crate::context::ProviderMode;

        let mut ctx = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
        ctx.set_strict_precision(Some(true));
        ctx.set_auto_spk_download(Some(false));
        let request = crate::positions::ServingRequest {
            ctx: &ctx,
            delta_t_model: crate::time::DeltaTModel::EspenakMeeus,
            iers_enabled: false,
            observer: None,
            tier_policy: crate::positions::TierPolicy::AttachOrder,
            skip_closed: false,
        };
        // CHIRON has curated elements but no kernel here: the original
        // strict-mode contract (register an SPK or enable downloads).
        assert!(matches!(
            crate::positions::calc_horizons(
                request,
                2451545.0,
                crate::constants::CHIRON,
                crate::constants::FLG_SPEED
            ),
            Err(Error::SpkRequired {
                body: crate::constants::CHIRON
            })
        ));

        assert!(!matches!(
            crate::positions::calc_horizons(request, 2451545.0, 99999, crate::constants::FLG_SPEED),
            Err(Error::SpkRequired { .. } | Error::SpkRegisteredOutOfCoverage { .. })
        ));
        // Registered to the tiny kernel (stored span J2000 +/- 0.5 d)
        // but requested far outside it: the registered-but-uncovered
        // variant, not the register-a-kernel one.
        let mut covered = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
        covered.set_strict_precision(Some(true));
        covered.set_auto_spk_download(Some(false));
        covered
            .spk_registry_mut()
            .register_image("tiny.bsp", tiny_kernel_image())
            .expect("fixture must decode");
        covered
            .spk_registry_mut()
            .register_body(crate::constants::CHIRON, "tiny.bsp", 4)
            .expect("chain must register");
        let registered = crate::positions::ServingRequest {
            ctx: &covered,
            delta_t_model: crate::time::DeltaTModel::EspenakMeeus,
            iers_enabled: false,
            observer: None,
            tier_policy: crate::positions::TierPolicy::AttachOrder,
            skip_closed: false,
        };
        assert!(matches!(
            crate::positions::calc_horizons(
                registered,
                2460000.5,
                crate::constants::CHIRON,
                crate::constants::FLG_SPEED
            ),
            Err(Error::SpkRegisteredOutOfCoverage {
                body: crate::constants::CHIRON,
                ..
            })
        ));
        // Inside the stored span the gate opens (a non-type-21 kernel
        // has no usable margins): whatever the kernel-free downstream
        // reports, it is not a strict-mode refusal.
        assert!(!matches!(
            crate::positions::calc_horizons(
                registered,
                2451545.5,
                crate::constants::CHIRON,
                crate::constants::FLG_SPEED
            ),
            Err(Error::SpkRequired { .. } | Error::SpkRegisteredOutOfCoverage { .. })
        ));

        let mut attempted = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
        attempted.set_strict_precision(Some(true));
        attempted.set_auto_spk_download(Some(true));
        let attempted_request = crate::positions::ServingRequest {
            ctx: &attempted,
            delta_t_model: crate::time::DeltaTModel::EspenakMeeus,
            iers_enabled: false,
            observer: None,
            tier_policy: crate::positions::TierPolicy::AttachOrder,
            skip_closed: false,
        };
        assert!(!matches!(
            crate::positions::calc_horizons(
                attempted_request,
                2451545.0,
                crate::constants::CHIRON,
                crate::constants::FLG_SPEED
            ),
            Err(Error::SpkRequired { .. } | Error::SpkRegisteredOutOfCoverage { .. })
        ));
    }
}

mod frames {
    use super::*;
    use crate::constants::{FLG_EQUATORIAL, FLG_ICRS, FLG_J2000, FLG_NONUT, FLG_SIDEREAL};

    /// Identity-frame bundle: every matrix is `I`, so each
    /// `frame_transform` branch reduces to its bare rotation and the
    /// branch choice itself is what the test pins.
    fn toy_frames() -> FrameMatrices {
        let ident = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        FrameMatrices::new(ident, ident, ident, ident, 0.41, 0.409).expect("toy bundle must build")
    }

    fn max_abs_diff(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> f64 {
        let mut m: f64 = 0.0;
        for i in 0..3 {
            for j in 0..3 {
                m = m.max((a[i][j] - b[i][j]).abs());
            }
        }
        m
    }

    fn is_orthogonal(m: [[f64; 3]; 3], tol: f64) -> bool {
        let mut mt = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                mt[i][j] = m[j][i];
            }
        }
        let prod = mat3_mul(m, mt).expect("finite");
        let ident = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        max_abs_diff(prod, ident) < tol
    }

    #[test]
    fn j2000_polynomial_pins_land_on_bare_constants() {
        assert_eq!(
            mean_obliquity_iau2006_deg(FRAME_J2000).unwrap(),
            84381.406 / 3600.0
        );
        assert_eq!(general_precession_deg(FRAME_J2000).unwrap(), 0.0);
        assert_eq!(
            general_precession_rate_deg_day(FRAME_J2000).unwrap(),
            5028.796195 / (3600.0 * 36525.0)
        );
    }

    #[test]
    fn fw_angles_at_j2000_match_bare_constants() {
        let (gamb, phib, psib, epsa) =
            iau2006_precession_angles(FRAME_J2000).expect("finite epoch");
        assert_eq!(gamb, -0.052928 * ARCSEC_TO_RAD);
        assert_eq!(phib, 84381.412819 * ARCSEC_TO_RAD);
        assert_eq!(psib, -0.041775 * ARCSEC_TO_RAD);
        assert_eq!(epsa, 84381.406 * ARCSEC_TO_RAD);
    }

    #[test]
    fn fw2m_is_a_rotation_at_j2000() {
        let (gamb, phib, psib, epsa) =
            iau2006_precession_angles(FRAME_J2000).expect("finite epoch");
        let r = fw2m(gamb, phib, psib, epsa).expect("finite angles");
        // Measured 3.4e-16: ordinary rounding over ~20 flops.
        assert!(is_orthogonal(r, 1e-15));
    }

    #[test]
    fn equatorial_ecliptic_rotation_round_trips() {
        let v = [0.3, -0.8, 0.5];
        let eps = 0.409092804018468;
        let ecl = rotate_equatorial_to_ecliptic(v, eps).expect("finite");
        let back = rotate_equatorial_to_ecliptic(ecl, -eps).expect("finite");
        for k in 0..3 {
            assert!((back[k] - v[k]).abs() < 2e-16, "k={k}");
        }
    }

    #[test]
    fn ecliptic_rotation_preserves_length_and_has_an_inverse() {
        // Two finite rotations exercise cancellation near the ecliptic plane.
        // The scalar bounds follow binary64 rotation error, independently of
        // a particular platform's evaluation order.
        let cases: [([u64; 3], u64, [u64; 3]); 2] = [
            (
                [0xbfd8cc2f6b4e5ca4, 0xc001332047be257d, 0xbfee8c0ff33f435e],
                0x3fda3e8603048b73,
                [0xbfd8cc2f6b4e5ca4, 0xc002d17d129a1b6a, 0xbf92c5176a71f520],
            ),
            (
                [0xbfe78881220fad9d, 0xbffa7bac87e788fd, 0xbfe655526ef0af1f],
                0x3fda15c2742fee72,
                [0xbfe78881220fad9d, 0xbffcbd7166811974, 0x3f8f71462dfcff40],
            ),
        ];
        for (v_bits, eps_bits, want_bits) in cases {
            let v = v_bits.map(f64::from_bits);
            let got = rotate_equatorial_to_ecliptic(v, f64::from_bits(eps_bits))
                .expect("finite rotation");
            let scale = v.iter().map(|x| x * x).sum::<f64>().sqrt();
            let inverse = rotate_equatorial_to_ecliptic(got, -f64::from_bits(eps_bits)).unwrap();
            for k in 0..3 {
                assert!(
                    (got[k] - f64::from_bits(want_bits[k])).abs() <= 16.0 * f64::EPSILON * scale
                );
                assert!((inverse[k] - v[k]).abs() <= 16.0 * f64::EPSILON * scale);
            }
            let norm = got.iter().map(|x| x * x).sum::<f64>().sqrt();
            assert!((norm - scale).abs() <= 16.0 * f64::EPSILON * scale);
        }
    }

    #[test]
    fn mat3_identity_laws_hold() {
        let ident = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let v = [1.0, -2.0, 0.5];
        assert_eq!(mat3_vec3(ident, v).unwrap(), v);
        let m = [[1.0, 2.0, 3.0], [0.0, -1.0, 4.0], [5.0, 0.0, -2.0]];
        assert_eq!(mat3_mul(m, ident).unwrap(), m);
        assert_eq!(mat3_mul(ident, m).unwrap(), m);
    }

    #[test]
    fn cartesian_spherical_covers_axes_wrap_and_origin() {
        assert_eq!(
            cartesian_to_spherical_deg([1.0, 0.0, 0.0]).unwrap(),
            [0.0, 0.0, 1.0]
        );
        assert_eq!(
            cartesian_to_spherical_deg([-1.0, 0.0, 0.0]).unwrap(),
            [180.0, 0.0, 1.0]
        );
        assert_eq!(
            cartesian_to_spherical_deg([0.0, 0.0, 0.0]).unwrap(),
            [0.0, 0.0, 0.0]
        );
        // Negative-y wraps forward into [0, 360).
        let s = cartesian_to_spherical_deg([0.0, -1.0, 0.0]).unwrap();
        assert_eq!(s, [270.0, 0.0, 1.0]);

        assert_eq!(
            cartesian_to_spherical_deg([-0.0, -0.0, -0.0]).unwrap(),
            [0.0, 0.0, 0.0]
        );
    }

    #[test]
    fn cotrans_normalizes_only_longitude_zero_sign() {
        for lon in [0.0, -0.0] {
            for lat in [0.0, -0.0] {
                for eps in [0.0, -0.0, 360.0, -360.0] {
                    let (longitude, latitude) = cotrans(lon, lat, eps).unwrap();
                    assert_eq!(longitude.to_bits(), 0.0_f64.to_bits());
                    assert_eq!(latitude, 0.0);
                }
            }
        }
        // Longitude canonicalization must not canonicalize computed latitude.
        assert_eq!(
            cotrans(0.0, -0.0, 0.0).unwrap().1.to_bits(),
            (-0.0_f64).to_bits()
        );
        assert_eq!(
            cotrans(-0.0, -0.0, 0.0).unwrap().1.to_bits(),
            0.0_f64.to_bits()
        );
        // Tiny angles and full-turn trig residuals retain the existing path:
        // no snap band or extra reduction of a rounded 360-degree result.
        for lon in [1e-14, -1e-14, 360.0, -360.0, 720.0, -720.0] {
            assert_ne!(cotrans(lon, 0.0, 0.0).unwrap().0, 0.0);
        }
    }

    #[test]
    fn cotrans_round_trips_through_both_signs() {
        let eps = 23.439279444444445;
        let (ra, dec) = cotrans(100.0, 20.0, eps).expect("finite");
        let (lon, lat) = cotrans(ra, dec, -eps).expect("finite");
        assert!((lon - 100.0).abs() < 1e-9, "lon={lon}");
        assert!((lat - 20.0).abs() < 1e-9, "lat={lat}");
    }

    #[test]
    fn nutation_matrix_is_identity_at_zero_nutation() {
        let eps = 0.409092804018468;
        let n = nutation_matrix(eps, 0.0, 0.0).expect("finite");
        let ident = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        assert!(max_abs_diff(n, ident) < 3e-16);
        // A sample nutation stays a rotation.
        let n2 = nutation_matrix(eps, 1.0e-6, 5.0e-7).expect("finite");
        assert!(is_orthogonal(n2, 1e-15));
    }

    #[test]
    fn nutation_matrix_preserves_erfa_rotation_words() {
        // PyERFA 2.0.1.5 numat outputs, with all arguments in radians.
        // Small off-diagonal terms expose cancellation and rotation order;
        // an orthogonality-only check cannot distinguish this defect.
        let cases = [
            (
                (0.4090926005957347, 1.2e-5, -4e-5),
                [
                    [0x3feffffffff61aba, 0xbee716d56d44e4b6, 0xbed405513556d5de],
                    [0x3ee716efaad0a3f8, 0x3fefffffff89b806, 0x3f04f8b4a1a8b7a4],
                    [0x3ed404d8274759c2, 0xbf04f8b66fe7963f, 0x3fefffffff907bb6],
                ],
            ),
            (
                (0.408987804, 7.133e-5, 1.516e-5),
                [
                    [0x3feffffffea25b81, 0xbf11281274ef3e7f, 0xbefdbe9808b2fe85],
                    [0x3f11280b1193e2ac, 0x3feffffffec9dca3, 0xbeefcb761ca353fd],
                    [0x3efdbedc37305d35, 0x3eefca76f3b2e4dc, 0x3fefffffffb8e89a],
                ],
            ),
            (
                (0.410042587, -8.252e-5, 3.119e-5),
                [
                    [0x3feffffffe2c0d1f, 0x3f13d6c25a1f1745, 0x3f013f49eb931bd7],
                    [0xbf13d6b0b9627a92, 0x3feffffffe339264, 0xbf005a6b027b1ada],
                    [0xbf013f9b06172d82, 0x3f005a157738c9e7, 0x3fefffffff72c6ce],
                ],
            ),
        ];
        for ((eps, dpsi, deps), expected) in cases {
            let actual = nutation_matrix(eps, dpsi, deps).unwrap();
            for (row, expected_row) in actual.iter().zip(expected) {
                for (value, expected_word) in row.iter().zip(expected_row) {
                    assert_eq!(
                        value.to_bits(),
                        expected_word,
                        "eps={eps}, dpsi={dpsi}, deps={deps}"
                    );
                }
            }
        }
    }

    #[test]
    fn vondrak_obliquity_at_j2000_matches_eps0() {
        let eps0_deg = 84381.406 / 3600.0;
        let pole_deg = vondrak_mean_obliquity_deg(FRAME_J2000).expect("finite");
        assert!((pole_deg - eps0_deg).abs() < 1e-9, "pole={pole_deg}");
        let series_deg = vondrak_mean_obliquity_series_rad(FRAME_J2000)
            .expect("finite")
            .to_degrees();
        assert!((series_deg - eps0_deg).abs() < 1e-9, "series={series_deg}");
    }

    #[test]
    fn vondrak_precession_at_j2000_is_near_identity() {
        let p = vondrak_precession_matrix(FRAME_J2000, false).expect("finite");
        let ident = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

        assert!(max_abs_diff(p, ident) < 1e-11);
        assert!(is_orthogonal(p, 1e-12));
    }

    #[test]
    fn vondrak_pn_with_zero_nutation_reproduces_precession() {
        let jd = 2460000.5;
        let (pn, eps_true) = vondrak_pn_matrix(jd, 0.0, 0.0, true).expect("finite");
        let p = vondrak_precession_matrix(jd, true).expect("finite");
        assert!(max_abs_diff(pn, p) < 1e-15);
        assert_eq!(eps_true, vondrak_mean_obliquity_rad(jd).unwrap());
        let (pn_free, _) = vondrak_pn_matrix(jd, 0.0, 0.0, false).expect("finite");
        let p_free = vondrak_precession_matrix(jd, false).expect("finite");
        assert!(max_abs_diff(pn_free, p_free) < 1e-15);
    }

    #[test]
    fn biased_precession_is_nobias_times_frame_bias() {
        let jd = 2370000.0;
        let biased = vondrak_precession_matrix(jd, true).expect("finite");
        let nobias = vondrak_precession_matrix(jd, false).expect("finite");
        // For orthogonal P, B = P F implies Pᵀ B = F. Check that identity
        // directly, independently of the production matrix-product helper.
        // gamma(12) bounds the short binary64 dot products and the input
        // rotations' accumulated rounding; this is far below 0.001 arcsec.
        let gamma12 = 12.0 * f64::EPSILON / (1.0 - 12.0 * f64::EPSILON);
        for (i, bias_row) in FRAME_BIAS_LTPB.iter().enumerate() {
            for (j, expected) in bias_row.iter().enumerate() {
                let recovered: f64 = (0..3).map(|k| nobias[k][i] * biased[k][j]).sum();
                assert!(
                    (recovered - expected).abs() <= gamma12,
                    "frame-bias entry ({i}, {j}): {recovered} vs {expected}"
                );
            }
        }
        for (index, rotation) in [nobias, biased].into_iter().enumerate() {
            for i in 0..3 {
                for j in 0..3 {
                    let inner: f64 = (0..3).map(|k| rotation[k][i] * rotation[k][j]).sum();
                    let expected = if index == 0 {
                        if i == j { 1.0 } else { 0.0 }
                    } else {
                        // The published first-order bias is not an exactly orthogonal
                        // rotation. Composition must preserve its Gram matrix.
                        (0..3)
                            .map(|k| FRAME_BIAS_LTPB[k][i] * FRAME_BIAS_LTPB[k][j])
                            .sum()
                    };
                    assert!((inner - expected).abs() <= gamma12);
                }
            }
        }
        // The bias itself is an ~8e-8 rotation, far above rounding.
        let ident = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let shift = max_abs_diff(biased, nobias);
        assert!(shift > 1e-9 && shift < 2e-7, "shift={shift:e}");
        assert!(max_abs_diff(FRAME_BIAS_J2000, ident) < 2e-7);
    }

    #[test]
    fn frame_transform_routes_every_release_branch() {
        let frames = toy_frames();
        let geo = [1.0, 0.0, 0.0];
        let jd = 2460000.5;
        // Default: true ecliptic of date with the true obliquity.
        let got = frame_transform(geo, jd, 0, false, frames).expect("finite");
        let want_vec = rotate_equatorial_to_ecliptic(geo, 0.41).expect("finite");
        assert_eq!(got, cartesian_to_spherical_deg(want_vec).unwrap());
        // Mean branch: same vector rotated by the mean obliquity instead.
        let got_mean = frame_transform(geo, jd, FLG_NONUT, false, frames).expect("finite");
        let want_mean = rotate_equatorial_to_ecliptic(geo, 0.409).expect("finite");
        assert_eq!(got_mean, cartesian_to_spherical_deg(want_mean).unwrap());
        // Equatorial J2000 without ICRS applies the frame bias; with
        // ICRS the vector passes through untouched.
        let biased =
            frame_transform(geo, jd, FLG_EQUATORIAL | FLG_J2000, true, frames).expect("finite");
        assert_eq!(biased, mat3_vec3(FRAME_BIAS_J2000, geo).unwrap());
        let plain = frame_transform(geo, jd, FLG_EQUATORIAL | FLG_J2000 | FLG_ICRS, true, frames)
            .expect("finite");
        assert_eq!(plain, geo);
        // Sidereal equator and plain mean equator both use P (identity
        // in the toy bundle).
        assert_eq!(
            frame_transform(geo, jd, FLG_EQUATORIAL | FLG_SIDEREAL, true, frames).expect("finite"),
            geo
        );
        assert_eq!(
            frame_transform(geo, jd, FLG_EQUATORIAL | FLG_NONUT, true, frames).expect("finite"),
            geo
        );
        // ICRS + default uses the bias-free true matrix (identity here)
        // with the true obliquity.
        let icrs = frame_transform(geo, jd, FLG_ICRS, true, frames).expect("finite");
        assert_eq!(icrs, want_vec);
    }

    #[test]
    fn frame_transform_rejects_nonfinite_inputs() {
        let frames = toy_frames();
        assert!(matches!(
            frame_transform([f64::NAN, 0.0, 0.0], 2460000.5, 0, false, frames).unwrap_err(),
            Error::NonFiniteArgument { .. }
        ));
        assert!(matches!(
            frame_transform([1.0, 0.0, 0.0], f64::INFINITY, 0, false, frames).unwrap_err(),
            Error::NonFiniteArgument { .. }
        ));
        assert!(matches!(
            FrameMatrices::new(
                [[f64::NAN, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                0.41,
                0.409,
            )
            .unwrap_err(),
            Error::NonFiniteArgument { .. }
        ));
    }

    #[test]
    fn j2000_native_matrices_are_rotations() {
        assert!(is_orthogonal(FRAME_BIAS_J2000, 1e-15));
        assert!(is_orthogonal(icrs_to_ecl_j2000_matrix(), 1e-15));
    }
}

mod topo {
    use super::*;
    use crate::constants::{EARTH, MOON, SUN};
    use crate::context::Observer;

    const IDENT: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    /// Exact `x/y` swap (orthogonal, binary-exact entries) to exercise
    /// the transpose path without inexact trigonometry.
    const SWAP_XY: [[f64; 3]; 3] = [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];

    fn max_abs_diff(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> f64 {
        let mut m: f64 = 0.0;
        for i in 0..3 {
            for j in 0..3 {
                m = m.max((a[i][j] - b[i][j]).abs());
            }
        }
        m
    }

    fn is_orthogonal(m: [[f64; 3]; 3], tol: f64) -> bool {
        let mut mt = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                mt[i][j] = m[j][i];
            }
        }
        let prod = mat3_mul(m, mt).expect("finite");
        max_abs_diff(prod, IDENT) < tol
    }

    fn vec_max_abs_diff(a: [f64; 3], b: [f64; 3]) -> f64 {
        (a[0] - b[0])
            .abs()
            .max((a[1] - b[1]).abs())
            .max((a[2] - b[2]).abs())
    }

    #[test]
    fn constants_pin_release_identities() {
        assert_eq!(WGS84_SEMI_MAJOR_M, 6_378_137.0);
        assert_eq!(WGS84_INV_FLATTENING, 298.257223563);
        assert_eq!(EARTH_OMEGA_RAD_S, 7.2921150e-5);
        assert_eq!(SECONDS_PER_DAY, 86400.0);
        assert_eq!(ERA_J2000_TURNS, 0.7790572732640);
        assert_eq!(ERA_RATE_TURNS_PER_UT1_DAY, 1.002_737_811_911_354_6);
        assert_eq!(AU_IN_M, 149_597_870_700.0);
    }

    #[test]
    fn geodetic_equator_pins_are_bitwise() {
        // lon 0, lat 0, h 0: cos/sin are exact, so x == a, y == z == 0.
        assert_eq!(
            geodetic_to_itrf_m(0.0, 0.0, 0.0).unwrap(),
            [6_378_137.0, 0.0, 0.0]
        );
        // Height shifts the equatorial radius exactly (values < 2^53).
        assert_eq!(
            geodetic_to_itrf_m(0.0, 0.0, 100.0).unwrap(),
            [6_378_237.0, 0.0, 0.0]
        );
        // Antimeridian: cos(π) == −1 exactly; sin(π) ~ 1.2e−16 leaks
        // sub-nanometre y, so only x is a bitwise pin.
        let anti = geodetic_to_itrf_m(std::f64::consts::PI, 0.0, 0.0).unwrap();
        assert_eq!(anti[0], -6_378_137.0);
        assert!(anti[1].abs() < 1e-9, "y = {}", anti[1]);
        assert_eq!(anti[2], 0.0);
    }

    #[test]
    fn geodetic_mid_latitude_pins_match_ephemeral_erfa() {
        // WGS84 reference values from ERFA gd2gc (SOFA-derived, BSD-licensed).
        // See `docs/sources.md` for attribution.
        let rome = geodetic_to_itrf_m(12.5_f64.to_radians(), 41.9_f64.to_radians(), 0.0).unwrap();
        assert_eq!(
            rome,
            [4641725.118266695, 1029045.6841753953, 4237343.073052207]
        );
        let nyc =
            geodetic_to_itrf_m((-74.0_f64).to_radians(), 40.7_f64.to_radians(), 10.0).unwrap();
        assert_eq!(
            nyc,
            [1334743.4315274705, -4654803.52193068, 4137229.2378602256]
        );
    }

    #[test]
    fn geodetic_poles_and_south_sign() {
        let north = geodetic_to_itrf_m(0.7, std::f64::consts::FRAC_PI_2, 0.0).unwrap();
        let south = geodetic_to_itrf_m(0.7, -std::f64::consts::FRAC_PI_2, 0.0).unwrap();
        // Polar radius (1 − e²)·N(90°) ≈ 6356752.314 m; z dominates,
        // x/y carry only the cos(π/2) ~ 6e−17 residue.
        assert!(
            (north[2] - 6_356_752.314_245_179).abs() < 1e-6,
            "z = {}",
            north[2]
        );
        assert!(north[0].abs() < 1e-9 && north[1].abs() < 1e-9);
        assert_eq!(south[2], -north[2]);
        // Longitude periodicity: 200° and −160° coincide to a few ulps
        // (the degree→radian scalings differ by ~1 ulp of 2π, amplified
        // by the ~6e6 m radius; 1e−8 m is ~10 ulps).
        let a = geodetic_to_itrf_m(200.0_f64.to_radians(), 0.1, 10.0).unwrap();
        let b = geodetic_to_itrf_m((-160.0_f64).to_radians(), 0.1, 10.0).unwrap();
        assert!(
            vec_max_abs_diff(a, b) < 1e-8,
            "diff = {}",
            vec_max_abs_diff(a, b)
        );
    }

    #[test]
    fn geodetic_rejects_nonfinite() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(geodetic_to_itrf_m(bad, 0.0, 0.0).is_err());
            assert!(geodetic_to_itrf_m(0.0, bad, 0.0).is_err());
            assert!(geodetic_to_itrf_m(0.0, 0.0, bad).is_err());
        }
    }

    #[test]
    fn era_pins_j2000_and_rate() {
        // ERA at J2000 is 0.7790572732640 turns = 280.46061837504°.
        let era0 = earth_rotation_angle_rad(FRAME_J2000).unwrap();
        let deg = era0 / std::f64::consts::TAU * 360.0;
        assert!((deg - 280.46061837504).abs() < 1e-9, "deg = {deg}");
        assert!((0.0..std::f64::consts::TAU).contains(&era0));
        // One mean solar day advances the angle by exactly the rate
        // (modulo 2π); the wrap lands back in range.
        let era1 = earth_rotation_angle_rad(FRAME_J2000 + 1.0).unwrap();
        let mut d = (era1 - era0) / std::f64::consts::TAU;
        d = d.rem_euclid(1.0);
        assert!((d - 0.00273781191135448).abs() < 1e-12, "d = {d}");
        assert!(!earth_rotation_angle_rad(f64::NAN).is_ok());
        assert!(!earth_rotation_angle_rad(f64::INFINITY).is_ok());
    }

    #[test]
    fn era_far_dates_match_reference_angles_and_local_rotation_rate() {
        // Independent ERFA orientation anchors retain their original values.
        // PLAN compares physical angles, independent of product association.
        let turn = std::f64::consts::TAU;
        let angle_limit = (0.001_f64 / 3600.0).to_radians();
        let principal =
            |angle: f64| (angle + std::f64::consts::PI).rem_euclid(turn) - std::f64::consts::PI;
        for (jd, expected) in [
            (1507900.5f64, 4.762058120203221f64),
            (2460600.0f64, 3.581064587166061f64),
        ] {
            let got = earth_rotation_angle_rad(jd).unwrap();
            assert!((0.0..turn).contains(&got));
            assert!(
                principal(got - expected).abs() <= angle_limit,
                "jd {jd}: {got} vs reference {expected}"
            );
            // Fractional-day rotation must retain the IAU UT1 spin rate,
            // including at distant epochs with many accumulated full turns.
            for step in [1.0 / 8.0, 1.0 / 1024.0] {
                let next = earth_rotation_angle_rad(jd + step).unwrap();
                let expected_advance = turn * 1.002_737_811_911_354_6 * step;
                let error = principal(next - got - expected_advance).abs();
                assert!(
                    error <= 128.0 * f64::EPSILON * turn,
                    "local rotation at {jd}, step {step}: angular error {error}"
                );
            }
        }
    }

    #[test]
    fn polar_motion_is_identity_at_zero_and_a_rotation() {
        assert_eq!(polar_motion_matrix(0.0).unwrap(), IDENT);
        let w =
            polar_motion_matrix(-47e-6 / 3600.0 * std::f64::consts::PI / 180.0).expect("finite s'");
        assert!(is_orthogonal(w, 1e-15));
        // R3(−s′) with s′ < 0 rotates by a positive micro-angle: the
        // off-diagonal has the active-layout sign.
        assert!(w[0][1] < 0.0 && w[1][0] > 0.0);
        assert!(!polar_motion_matrix(f64::NAN).is_ok());
    }

    #[test]
    fn c2t_assembly_orders_pom_theta_npbi() {
        // Identities compose to the identity (bitwise: zeros stay zero).
        assert_eq!(
            assemble_celestial_to_terrestrial(IDENT, 0.3, IDENT).unwrap(),
            rot_z_rad(0.3)
        );
        // Orthogonal inputs stay orthogonal through the assembly.
        let c2t = assemble_celestial_to_terrestrial(
            polar_motion_matrix(1e-9).unwrap(),
            1.234,
            FRAME_BIAS_J2000,
        )
        .unwrap();
        assert!(is_orthogonal(c2t, 1e-15));
        // The ERA slot is a genuine rotation: different angles differ.
        let other = assemble_celestial_to_terrestrial(IDENT, 2.345, IDENT).unwrap();
        assert!(max_abs_diff(c2t, other) > 1e-3);
        assert!(!assemble_celestial_to_terrestrial(IDENT, f64::NAN, IDENT).is_ok());
    }

    #[test]
    fn offset_with_identity_matrix_matches_itrf_over_au() {
        let off = topocentric_offset(30.0, 45.0, 100.0, IDENT, SiteEllipsoid::Wgs84).unwrap();
        let itrf = geodetic_to_itrf_m(30.0_f64.to_radians(), 45.0_f64.to_radians(), 100.0).unwrap();
        let expect_pos = [itrf[0] / AU_IN_M, itrf[1] / AU_IN_M, itrf[2] / AU_IN_M];
        assert_eq!(off.pos_au(), expect_pos);
        assert_eq!(off.pair().0, expect_pos);
        // ω × r structure: vx opposes y, vy follows x, vz is exact zero.
        let k = SECONDS_PER_DAY / AU_IN_M;
        assert_eq!(off.vel_au_day()[0], -EARTH_OMEGA_RAD_S * itrf[1] * k);
        assert_eq!(off.vel_au_day()[1], EARTH_OMEGA_RAD_S * itrf[0] * k);
        assert_eq!(off.vel_au_day()[2], 0.0);
        assert_eq!(off.pair().1, off.vel_au_day());
        // Norm preserved by the orthogonal map (1-ulp-class agreement).
        let norm_au = (off.pos_au()[0].powi(2) + off.pos_au()[1].powi(2) + off.pos_au()[2].powi(2))
            .sqrt()
            * AU_IN_M;
        let norm_m = (itrf[0].powi(2) + itrf[1].powi(2) + itrf[2].powi(2)).sqrt();
        assert!((norm_au - norm_m).abs() / norm_m < 1e-15);
    }

    #[test]
    fn offset_transpose_path_swaps_components() {
        let base = topocentric_offset(30.0, 45.0, 100.0, IDENT, SiteEllipsoid::Wgs84).unwrap();
        let swapped = topocentric_offset(30.0, 45.0, 100.0, SWAP_XY, SiteEllipsoid::Wgs84).unwrap();
        // SWAP_XY is symmetric, so it is its own transpose: x/y swap.
        assert_eq!(swapped.pos_au()[0], base.pos_au()[1]);
        assert_eq!(swapped.pos_au()[1], base.pos_au()[0]);
        assert_eq!(swapped.pos_au()[2], base.pos_au()[2]);
        assert_eq!(swapped.vel_au_day()[0], base.vel_au_day()[1]);
        assert_eq!(swapped.vel_au_day()[1], base.vel_au_day()[0]);
    }

    #[test]
    fn geodetic_general_site_words_are_pinned() {
        let itrf = geodetic_to_itrf_m(30.0_f64.to_radians(), 45.0_f64.to_radians(), 100.0).unwrap();
        let got = [itrf[0].to_bits(), itrf[1].to_bits(), itrf[2].to_bits()];
        assert_eq!(
            got,
            [
                0x414d_d96c_d9e2_b9b9,
                0x4141_3bc7_65ba_cfa9,
                0x4151_1e3e_c7a6_9c0a,
            ],
            "itrf words moved: {itrf:?}"
        );
    }

    #[test]
    fn site_ellipsoid_selects_published_constants() {
        assert_eq!(IERS2010_SEMI_MAJOR_M, 6_378_136.6);
        assert_eq!(IERS2010_INV_FLATTENING, 298.25642);
        assert_eq!(SiteEllipsoid::default(), SiteEllipsoid::Wgs84);
        assert_eq!(SiteEllipsoid::Wgs84.semi_major_m(), WGS84_SEMI_MAJOR_M);
        assert_eq!(SiteEllipsoid::Wgs84.inv_flattening(), WGS84_INV_FLATTENING);
        assert_eq!(
            SiteEllipsoid::Iers2010.semi_major_m(),
            IERS2010_SEMI_MAJOR_M
        );
        assert_eq!(
            SiteEllipsoid::Iers2010.inv_flattening(),
            IERS2010_INV_FLATTENING
        );
    }

    #[test]
    fn iers2010_itrf_satisfies_ellipsoid_identities() {
        // Closed-form identities on the IERS2010 ellipsoid: the equator
        // radius is `a + h`, the pole radius `a(1 − f) + h`, and the
        // default WGS84 entry point stays bitwise identical to the
        // explicit-ellipsoid spelling.
        let a = IERS2010_SEMI_MAJOR_M;
        let f = 1.0 / IERS2010_INV_FLATTENING;
        let eq = geodetic_to_itrf_m_in(0.3, 0.0, 100.0, SiteEllipsoid::Iers2010).unwrap();
        let eq_r = (eq[0].powi(2) + eq[1].powi(2)).sqrt();
        assert!((eq_r - (a + 100.0)).abs() / (a + 100.0) < 1e-15);
        assert!(eq[2].abs() < 1e-6);
        let pole = geodetic_to_itrf_m_in(
            0.3,
            std::f64::consts::FRAC_PI_2,
            100.0,
            SiteEllipsoid::Iers2010,
        )
        .unwrap();
        assert!(pole[0].abs() < 1e-6 && pole[1].abs() < 1e-6);
        assert!((pole[2] - ((1.0 - f) * a + 100.0)).abs() < 1e-6);
        let lon = 13.405_f64.to_radians();
        let lat = 52.52_f64.to_radians();
        let via_default = geodetic_to_itrf_m(lon, lat, 34.0).unwrap();
        let via_explicit = geodetic_to_itrf_m_in(lon, lat, 34.0, SiteEllipsoid::Wgs84).unwrap();
        assert_eq!(via_default, via_explicit);
        assert!(geodetic_to_itrf_m_in(f64::NAN, lat, 34.0, SiteEllipsoid::Iers2010).is_err());
    }

    #[test]
    fn iers2010_berlin_site_is_submetre_from_wgs84() {
        let lon = 13.405_f64.to_radians();
        let lat = 52.52_f64.to_radians();
        let wgs = geodetic_to_itrf_m_in(lon, lat, 34.0, SiteEllipsoid::Wgs84).unwrap();
        let grs = geodetic_to_itrf_m_in(lon, lat, 34.0, SiteEllipsoid::Iers2010).unwrap();
        let d = [grs[0] - wgs[0], grs[1] - wgs[1], grs[2] - wgs[2]];
        let norm = (d[0].powi(2) + d[1].powi(2) + d[2].powi(2)).sqrt();
        assert!(
            (0.30..0.60).contains(&norm),
            "berlin ellipsoid gap out of band: {d:?} ({norm} m)"
        );
        // Same ellipsoid through the offset path with an identity matrix.
        let off = topocentric_offset(13.405, 52.52, 34.0, IDENT, SiteEllipsoid::Iers2010).unwrap();
        let expect = [grs[0] / AU_IN_M, grs[1] / AU_IN_M, grs[2] / AU_IN_M];
        assert_eq!(off.pos_au(), expect);
        // The sampler carries the selection to its sample epochs.
        let s = TopoSampler::new(13.405, 52.52, 34.0, 0.25, SiteEllipsoid::Iers2010).unwrap();
        assert_eq!(s.ellipsoid(), SiteEllipsoid::Iers2010);
        assert_eq!(s.offset_at(2460000.5, IDENT).unwrap(), off);
    }

    #[test]
    // Direct matrix indices are the independent spelling this fixture is intended to check.
    #[allow(clippy::needless_range_loop)]
    fn offset_row_locals_match_direct_index_spelling() {
        // The row-local (`t2c_0/1/2`) grouping must be a bitwise no-op
        // versus direct `c2t[j][i]` indexing: same operands, same
        // left-to-right association, same single `/ AU` (pos) and
        // `* SECONDS_PER_DAY / AU` (vel) scaling.
        let c2t = rot_z_rad(0.7);
        let off = topocentric_offset(30.0, 45.0, 100.0, c2t, SiteEllipsoid::Wgs84).unwrap();
        let itrf = geodetic_to_itrf_m(30.0_f64.to_radians(), 45.0_f64.to_radians(), 100.0).unwrap();
        let (x_m, y_m, z_m) = (itrf[0], itrf[1], itrf[2]);
        let vx = -EARTH_OMEGA_RAD_S * y_m;
        let vy = EARTH_OMEGA_RAD_S * x_m;
        for i in 0..3 {
            let direct_pos = (c2t[0][i] * x_m + c2t[1][i] * y_m + c2t[2][i] * z_m) / AU_IN_M;
            let direct_vel =
                (c2t[0][i] * vx + c2t[1][i] * vy + c2t[2][i] * 0.0) * SECONDS_PER_DAY / AU_IN_M;
            assert_eq!(
                off.pos_au()[i].to_bits(),
                direct_pos.to_bits(),
                "pos row {i}"
            );
            assert_eq!(
                off.vel_au_day()[i].to_bits(),
                direct_vel.to_bits(),
                "vel row {i}"
            );
        }
    }

    #[test]
    fn offset_bundle_rejects_nonfinite() {
        assert!(TopocentricOffset::new([1.0, 0.0, 0.0], [0.0, f64::NAN, 0.0]).is_err());
        assert!(TopocentricOffset::new([f64::INFINITY, 0.0, 0.0], [0.0; 3]).is_err());
        assert!(topocentric_offset(f64::NAN, 0.0, 0.0, IDENT, SiteEllipsoid::Wgs84).is_err());
        assert!(topocentric_offset(0.0, 0.0, f64::INFINITY, IDENT, SiteEllipsoid::Wgs84).is_err());
        let mut bad = IDENT;
        bad[2][1] = f64::NAN;
        assert!(topocentric_offset(0.0, 0.0, 0.0, bad, SiteEllipsoid::Wgs84).is_err());
    }

    #[test]
    fn noop_set_matches_release_ids() {
        for id in [
            MEAN_NODE, TRUE_NODE, MEAN_APOG, OSCU_APOG, INTP_APOG, INTP_PERG,
        ] {
            assert!(is_topo_noop_point(id), "lunar point {id} must be a no-op");
        }
        for id in [SUN, MOON, EARTH, 1, 2, 3] {
            assert!(!is_topo_noop_point(id), "body {id} must not be a no-op");
        }
    }

    #[test]
    fn resolve_covers_flag_noop_override_and_missing_branches() {
        let obs = Observer::new(13.4, 52.5, 34.0).unwrap();
        // No flag: geocentric even with an observer and an override.
        assert_eq!(
            resolve_topo_request(0, SUN, Some((1.0, 2.0, 3.0)), Some(obs)),
            Ok(TopoPlan::Geocentric)
        );
        assert!(
            !resolve_topo_request(0, SUN, None, Some(obs))
                .unwrap()
                .uses_offset()
        );
        // Lunar points: no-op without any observer.
        assert_eq!(
            resolve_topo_request(FLG_TOPOCTR, MEAN_NODE, None, None),
            Ok(TopoPlan::TopoNoop)
        );
        assert!(
            !resolve_topo_request(FLG_TOPOCTR, OSCU_APOG, None, None)
                .unwrap()
                .uses_offset()
        );
        assert_eq!(
            resolve_topo_request(FLG_TOPOCTR, MEAN_NODE, None, None)
                .unwrap()
                .geopos(),
            None
        );
        // Explicit triple wins over the context observer.
        assert_eq!(
            resolve_topo_request(FLG_TOPOCTR, SUN, Some((1.0, 2.0, 3.0)), Some(obs)),
            Ok(TopoPlan::Topocentric {
                lon_deg: 1.0,
                lat_deg: 2.0,
                height_m: 3.0,
            })
        );
        // Context observer flows through unchanged.
        let plan = resolve_topo_request(FLG_TOPOCTR, MOON, None, Some(obs)).unwrap();
        assert!(plan.uses_offset());
        assert_eq!(plan.geopos(), Some((13.4, 52.5, 34.0)));
        let plan = resolve_topo_request(FLG_TOPOCTR, EARTH, None, Some(obs)).unwrap();
        assert_eq!(plan.geopos(), Some((13.4, 52.5, 34.0)));
        // Missing observer is a typed refusal, never silent geocentric.
        match resolve_topo_request(FLG_TOPOCTR, SUN, None, None) {
            Err(Error::ConfigurationMissing {
                missing_config,
                suggestion,
            }) => {
                assert_eq!(missing_config, "observer_location");
                assert!(suggestion.contains("FLG_TOPOCTR"));
            }
            other => panic!("expected ConfigurationMissing, got {other:?}"),
        }
    }

    #[test]
    fn sampler_freezes_the_central_scale_difference() {
        // Exactly representable epochs keep the frozen difference exact.
        assert_eq!(frozen_delta_t_days(2460000.5, 2460000.25).unwrap(), 0.25);
        assert!((frozen_delta_t_days(2460000.0, 2459999.99).unwrap() - 0.01).abs() < 1e-9);
        assert!(frozen_delta_t_days(f64::NAN, 0.0).is_err());
        let s = TopoSampler::new(13.4, 52.5, 34.0, 0.25, SiteEllipsoid::Wgs84).unwrap();
        assert_eq!(s.lon_deg(), 13.4);
        assert_eq!(s.lat_deg(), 52.5);
        assert_eq!(s.height_m(), 34.0);
        assert_eq!(s.delta_t_days(), 0.25);
        assert_eq!(s.ellipsoid(), SiteEllipsoid::Wgs84);
        // Every sample shifts by the same frozen difference.
        assert_eq!(s.ut1_at(2460000.5).unwrap(), 2460000.25);
        assert_eq!(s.ut1_at(2460001.5).unwrap(), 2460001.25);
        assert!(s.ut1_at(f64::INFINITY).is_err());
        assert!(TopoSampler::new(0.0, 0.0, 0.0, f64::NAN, SiteEllipsoid::Wgs84).is_err());
        // Same-epoch re-evaluation agrees with the direct offset path.
        let direct = topocentric_offset(13.4, 52.5, 34.0, IDENT, SiteEllipsoid::Wgs84).unwrap();
        assert_eq!(s.offset_at(2460000.5, IDENT).unwrap(), direct);
        assert!(s.offset_at(f64::NAN, IDENT).is_err());
    }
}

mod speeds {
    use super::*;
    use crate::constants::SUN;
    use crate::leb::CENTER_OFFSET_VELOCITY_STEP_DAYS;

    fn vec_close(a: [f64; 3], b: [f64; 3], tol: f64) -> bool {
        (a[0] - b[0])
            .abs()
            .max((a[1] - b[1]).abs())
            .max((a[2] - b[2]).abs())
            < tol
    }

    fn geocentric_center() -> StencilCenter {
        StencilCenter::new(
            [1.0, 0.5, -0.25],
            [0.01, -0.02, 0.005],
            [0.1, 0.2, 0.3],
            [0.001, 0.002, 0.003],
            [0.0; 3],
            5,
            false,
            false,
            None,
        )
        .expect("finite fixture")
    }

    #[test]
    fn step_pins_the_release_half_step() {
        assert_eq!(VEL_H_DAYS, 0.0005);

        assert_eq!(VEL_H_DAYS, CENTER_OFFSET_VELOCITY_STEP_DAYS);
    }

    #[test]
    fn nutation_rate_matches_analytic_slopes() {
        let rate = nutation_rate_deg_day(2460000.5, |jd| Ok(1e-5 + 2e-9 * jd)).unwrap();
        assert!(
            (rate - (2e-9_f64).to_degrees()).abs() < 1e-12,
            "rate = {rate}"
        );
        // Constant provider: zero rate, not noise.
        assert_eq!(nutation_rate_deg_day(2460000.5, |_| Ok(0.01)).unwrap(), 0.0);
        // Sine provider: central difference matches the analytic
        // derivative to O(h²) (here ~1e-16 deg/day; the provider is
        // written around the epoch so large-JD argument rounding
        // does not dominate the pin).
        let jd = 2460000.5;
        let rate = nutation_rate_deg_day(jd, |t| Ok(1e-6 * ((t - jd) * 0.01).sin())).unwrap();
        let expect = (1e-8_f64).to_degrees();
        assert!(
            (rate - expect).abs() < 1e-12,
            "rate = {rate}, expect = {expect}"
        );
        assert!(nutation_rate_deg_day(f64::NAN, |_| Ok(0.0)).is_err());
        assert!(
            nutation_rate_deg_day(2460000.5, |_| Err(Error::LebOutOfRange {
                body_id: 3,
                jd: 2460000.5,
                start: 0.0,
                end: 1.0,
            }))
            .is_err()
        );
    }

    #[test]
    fn extrapolation_and_curvature_follow_release_order() {
        let g = [1.0, 2.0, 3.0];
        let v = [0.1, -0.2, 0.05];
        // Zero offset reproduces the input bit-exactly.
        assert_eq!(extrapolate_geometric_vector(g, v, 0.0), g);
        // Verbatim association `(gv*s)` added component-wise.
        assert_eq!(extrapolate_geometric_vector(g, v, 0.5), [1.05, 1.9, 3.025]);
        // No curvature when the sample rides the linear path.
        let t0 = [0.01, 0.02, 0.03];
        let tv = [0.001, 0.002, 0.003];
        let tp = [
            t0[0] + tv[0] * 0.25,
            t0[1] + tv[1] * 0.25,
            t0[2] + tv[2] * 0.25,
        ];
        let e = extrapolate_geometric_vector(g, v, 0.25);
        assert_eq!(remove_topo_curvature(e, t0, tv, tp, 0.25), e);
        // A curved sample is removed exactly.
        let curved = [tp[0] + 1e-6, tp[1] - 2e-6, tp[2]];
        let fixed = remove_topo_curvature(e, t0, tv, curved, 0.25);
        assert_eq!(fixed, [e[0] - 1e-6, e[1] + 2e-6, e[2]]);
        // Sample light time normalizes by the extrapolated distance.
        assert_eq!(
            sample_light_time_days([3.0, 4.0, 0.0]),
            5.0 / C_LIGHT_AU_DAY
        );
    }

    #[test]
    fn sun_target_freezes_the_deflection_delta() {
        let mut center = geocentric_center();
        center = StencilCenter::new(
            center.geo_geom(),
            center.geo_vel(),
            center.observer(),
            center.observer_vel(),
            [1e-9, -2e-9, 3e-9],
            SUN,
            true,
            false,
            None,
        )
        .expect("finite fixture");
        let g = sample_apparent_vector(&center, VEL_H_DAYS, [0.0; 3], 2460000.5, |_, _, _, _| {
            panic!("Sun freeze must not re-evaluate deflection")
        })
        .unwrap();
        let expect = extrapolate_geometric_vector(center.geo_geom(), center.geo_vel(), VEL_H_DAYS);
        assert_eq!(g, [expect[0] + 1e-9, expect[1] - 2e-9, expect[2] + 3e-9]);
    }

    #[test]
    fn non_sun_sample_reevaluates_deflection_at_sample_geometry() {
        let center = StencilCenter::new(
            [1.0, 0.0, 0.0],
            [0.01, 0.0, 0.0],
            [0.1, 0.0, 0.0],
            [0.001, 0.0, 0.0],
            [0.0; 3],
            5,
            true,
            false,
            None,
        )
        .expect("finite fixture");
        let g = sample_apparent_vector(
            &center,
            -VEL_H_DAYS,
            [0.0; 3],
            2459999.0,
            |g, obs_s, jd_s, lt_s| {
                assert_eq!(jd_s, 2459999.0);
                assert_eq!(lt_s, sample_light_time_days(g));
                assert_eq!(
                    obs_s,
                    extrapolate_geometric_vector(
                        center.observer(),
                        center.observer_vel(),
                        -VEL_H_DAYS
                    )
                );
                Ok([g[0] + 1.0, g[1], g[2]])
            },
        )
        .unwrap();
        let expect = extrapolate_geometric_vector(center.geo_geom(), center.geo_vel(), -VEL_H_DAYS);
        assert_eq!(g, [expect[0] + 1.0, expect[1], expect[2]]);
    }

    #[test]
    fn aberration_uses_sample_earth_plus_topo_velocities() {
        let topo = StencilTopo::new(
            [0.01, 0.0, 0.0],
            [0.001, 0.0, 0.0],
            Some(([0.010001, 0.0, 0.0], [0.0011, 0.0002, 0.0])),
        )
        .expect("finite topo");
        let center = StencilCenter::new(
            [1.0, 0.5, 0.0],
            [0.0; 3],
            [0.1, 0.0, 0.0],
            [0.0; 3],
            [0.0; 3],
            5,
            false,
            true,
            Some(topo),
        )
        .expect("finite fixture");
        let earth_s = [0.01, 0.02, 0.03];
        let g = sample_apparent_vector(&center, VEL_H_DAYS, earth_s, 2460000.5, |_, _, _, _| {
            panic!("deflection is gated off here")
        })
        .unwrap();
        // Sample-pair velocity wins over the central one. The sample
        // position carries real diurnal curvature (5e-7 AU off the
        // linear path), so the expectation applies the same
        // curvature removal and re-times from the curved vector.
        let raw = extrapolate_geometric_vector(center.geo_geom(), center.geo_vel(), VEL_H_DAYS);
        let curved = remove_topo_curvature(
            raw,
            topo.central_pos(),
            topo.central_vel(),
            [0.010001, 0.0, 0.0],
            VEL_H_DAYS,
        );
        let expect = apply_aberration(
            curved,
            [earth_s[0] + 0.0011, earth_s[1] + 0.0002, earth_s[2]],
            sample_light_time_days(curved),
        )
        .unwrap();
        assert_eq!(g, expect);

        let topo_central_only =
            StencilTopo::new([0.01, 0.0, 0.0], [0.001, 0.0, 0.0], None).expect("finite topo");
        let center_fallback = StencilCenter::new(
            center.geo_geom(),
            center.geo_vel(),
            center.observer(),
            center.observer_vel(),
            [0.0; 3],
            5,
            false,
            true,
            Some(topo_central_only),
        )
        .expect("finite fixture");
        let g2 = sample_apparent_vector(
            &center_fallback,
            VEL_H_DAYS,
            earth_s,
            2460000.5,
            |_, _, _, _| panic!("deflection is gated off here"),
        )
        .unwrap();
        let expect2 = apply_aberration(
            extrapolate_geometric_vector(center.geo_geom(), center.geo_vel(), VEL_H_DAYS),
            [earth_s[0] + 0.001, earth_s[1], earth_s[2]],
            sample_light_time_days(extrapolate_geometric_vector(
                center.geo_geom(),
                center.geo_vel(),
                VEL_H_DAYS,
            )),
        )
        .unwrap();
        assert_eq!(g2, expect2);
    }

    #[test]
    fn differencing_pins_wrap_xyz_and_rate_scale() {
        let inv2h = 1.0 / (2.0 * VEL_H_DAYS);
        // Plain xyz branch.
        assert_eq!(
            difference_stencil_rate((1.0, 2.0, 3.0), (1.5, 1.0, 4.0), true).unwrap(),
            (0.5 * inv2h, -inv2h, 1.0 * inv2h)
        );
        // Forward cut crossing: 359.9 -> 0.1 is +0.2, not -359.8.
        let (dlon, _, _) =
            difference_stencil_rate((359.9, 0.0, 1.0), (0.1, 0.0, 1.0), false).unwrap();
        assert!((dlon - 0.2 * inv2h).abs() < 1e-9, "dlon = {dlon}");
        // Backward crossing mirrors it.
        let (dlon, _, _) =
            difference_stencil_rate((0.1, 0.0, 1.0), (359.9, 0.0, 1.0), false).unwrap();
        assert!((dlon + 0.2 * inv2h).abs() < 1e-9, "dlon = {dlon}");

        let (dlon, _, _) =
            difference_stencil_rate((180.0, 0.0, 1.0), (0.0, 0.0, 1.0), false).unwrap();
        assert_eq!(dlon, 180.0 * inv2h);
        // Lat/dist difference plainly on the spherical branch.
        let (_, dlat, ddist) =
            difference_stencil_rate((10.0, 20.0, 1.0), (10.5, 19.0, 1.25), false).unwrap();
        assert_eq!((dlat, ddist), (-inv2h, 0.25 * inv2h));
        assert!(difference_stencil_rate((f64::NAN, 0.0, 0.0), (0.0, 0.0, 0.0), true).is_err());
    }

    #[test]
    fn full_stencil_recovers_the_analytic_rate() {
        // No corrections, identity frame: the stencil must recover the
        // analytic input rate to discretization noise (~1e-9 here; the
        // exact-bit path is covered by the extrapolation pins above).
        let center = StencilCenter::new(
            [1.0, -2.0, 0.5],
            [0.03, 0.07, -0.11],
            [0.0; 3],
            [0.0; 3],
            [0.0; 3],
            5,
            false,
            false,
            None,
        )
        .expect("finite fixture");
        let jd_tt = 2460000.5;
        let (rx, ry, rz) = stencil_velocity(
            &center,
            jd_tt,
            [0.0; 3],
            [0.0; 3],
            true,
            |_, _, _, _| panic!("no deflection on this path"),
            |g, jd| {
                // Identity frame at the exact sample epochs.
                assert!(jd == jd_tt - VEL_H_DAYS || jd == jd_tt + VEL_H_DAYS);
                Ok((g[0], g[1], g[2]))
            },
        )
        .unwrap();
        assert!(
            vec_close([rx, ry, rz], center.geo_vel(), 1e-9),
            "rate = [{rx}, {ry}, {rz}]"
        );
        assert!(
            stencil_velocity(
                &center,
                f64::INFINITY,
                [0.0; 3],
                [0.0; 3],
                true,
                |g, _, _, _| Ok(g),
                |g, _| Ok((g[0], g[1], g[2])),
            )
            .is_err()
        );
    }

    #[test]
    fn constructors_reject_nonfinite_inputs() {
        assert!(StencilTopo::new([f64::NAN, 0.0, 0.0], [0.0; 3], None).is_err());
        assert!(
            StencilTopo::new(
                [0.0; 3],
                [0.0; 3],
                Some(([0.0; 3], [f64::INFINITY, 0.0, 0.0]))
            )
            .is_err()
        );
        assert!(
            StencilCenter::new(
                [f64::NAN, 0.0, 0.0],
                [0.0; 3],
                [0.0; 3],
                [0.0; 3],
                [0.0; 3],
                5,
                false,
                false,
                None,
            )
            .is_err()
        );
        let center = geocentric_center();
        assert!(
            sample_apparent_vector(&center, f64::NAN, [0.0; 3], 2460000.5, |g, _, _, _| Ok(g))
                .is_err()
        );
        assert!(
            sample_apparent_vector(
                &center,
                0.0,
                [f64::NAN, 0.0, 0.0],
                2460000.5,
                |g, _, _, _| Ok(g)
            )
            .is_err()
        );
    }
}

mod entry_tests {
    use super::*;
    use crate::constants::{
        FLG_CENTER_BODY, FLG_JPLHOR, FLG_MOSEPH, FLG_SPEED, FLG_SPEED3, MARS, MEAN_NODE, TRUE_NODE,
    };

    type TimeConvert = dyn Fn(f64) -> Result<f64>;
    type ComputeTt = dyn Fn(f64, i32, i32) -> Result<(PositionResult, i32)>;
    type NutationFn = dyn Fn(f64, i32) -> Result<(PositionResult, i32)>;

    fn stub_engine() -> (
        Box<TimeConvert>,
        Box<ComputeTt>,
        Box<NutationFn>,
        Box<NutationFn>,
    ) {
        // Core provider: identity spherical result (lon 100 + body,
        // lat 10, dist 1.5, zero rates), retflag = computation flags.
        let compute = Box::new(|_jd: f64, body: i32, calc_flags: i32| {
            Ok(([100.0 + body as f64, 10.0, 1.5, 0.0, 0.0, 0.0], calc_flags))
        }) as Box<dyn Fn(f64, i32, i32) -> Result<(PositionResult, i32)>>;
        let nut =
            Box::new(|_jd: f64, flags: i32| Ok(([23.44, 23.43, 0.001, 0.002, 0.0, 0.0], flags)))
                as Box<dyn Fn(f64, i32) -> Result<(PositionResult, i32)>>;
        let nut2 =
            Box::new(|_jd: f64, flags: i32| Ok(([23.44, 23.43, 0.001, 0.002, 0.0, 0.0], flags)))
                as Box<dyn Fn(f64, i32) -> Result<(PositionResult, i32)>>;
        let ut1_to_tt = Box::new(|jd: f64| Ok(jd + 0.0008)) as Box<dyn Fn(f64) -> Result<f64>>;
        (ut1_to_tt, compute, nut, nut2)
    }

    fn engine_with<'a>(
        ut1_to_tt: &'a dyn Fn(f64) -> Result<f64>,
        compute: &'a dyn Fn(f64, i32, i32) -> Result<(PositionResult, i32)>,
        nut_ut: &'a dyn Fn(f64, i32) -> Result<(PositionResult, i32)>,
        nut_tt: &'a dyn Fn(f64, i32) -> Result<(PositionResult, i32)>,
    ) -> CalcEngine<'a> {
        CalcEngine {
            ut1_to_tt_jd: ut1_to_tt,
            compute_tt: compute,
            nutation_obliquity_ut: nut_ut,
            nutation_obliquity_tt: nut_tt,
            jpl_fallback_active: false,
        }
    }

    #[test]
    fn cartesian_entry_preserves_words_flags_and_tt_acquisition() {
        let (ut, core, nut_ut, nut_tt) = stub_engine();
        let engine = engine_with(&*ut, &*core, &*nut_ut, &*nut_tt);
        let components = [-2.0, 0.75, -0.25, -0.5, 0.125, -0.0625];
        for format in [FLG_XYZ, FLG_XYZ | FLG_RADIANS] {
            let flags = FLG_SWIEPH | FLG_SPEED | format;
            let cartesian = |jd, body, calc_flags| {
                assert_eq!(jd, 100.0);
                assert_eq!(body, MARS);
                assert_eq!(calc_flags, FLG_SWIEPH | FLG_SPEED);
                Ok(Some((components, calc_flags)))
            };
            let actual = calc_with_cartesian(&engine, &cartesian, 100.0, MARS, flags).unwrap();
            assert_eq!(actual, (components, flags));
            let cartesian_ut = |jd, body, calc_flags| {
                assert_eq!(jd, 100.0 + 0.0008);
                assert_eq!(body, MARS);
                assert_eq!(calc_flags, FLG_SWIEPH | FLG_SPEED);
                Ok(Some((components, calc_flags)))
            };
            assert_eq!(
                calc_ut_with_cartesian(&engine, &cartesian_ut, 100.0, MARS, flags).unwrap(),
                actual
            );
        }
    }

    #[test]
    fn cartesian_entry_retains_spherical_fallback_and_special_dispatch() {
        let (ut, core, nut_ut, nut_tt) = stub_engine();
        let engine = engine_with(&*ut, &*core, &*nut_ut, &*nut_tt);
        let absent = |_, _, _| Ok(None);
        let flags = FLG_SWIEPH | FLG_SPEED | FLG_XYZ | FLG_RADIANS;
        assert_eq!(
            calc_with_cartesian(&engine, &absent, 100.0, MARS, flags).unwrap(),
            calc(&engine, 100.0, MARS, flags).unwrap(),
        );
        let unreachable = |_, _, _| panic!("direct provider must not enter special dispatch");
        for (body, request) in [
            (MARS, flags & !FLG_XYZ),
            (MARS, flags | FLG_SIDEREAL),
            (-MEAN_NODE, flags),
            (-TRUE_NODE, flags),
            (MEAN_APOG, (flags & !FLG_SPEED) | FLG_SPEED3),
            (ECL_NUT, flags),
        ] {
            assert_eq!(
                calc_with_cartesian(&engine, &unreachable, 100.0, body, request).unwrap(),
                calc(&engine, 100.0, body, request).unwrap(),
            );
        }
        let failed = |_, id, _| Err(Error::UnknownBody { body_id: id });
        assert_eq!(
            calc_with_cartesian(&engine, &failed, 100.0, MARS, flags),
            Err(Error::UnknownBody { body_id: MARS }),
        );
        assert!(matches!(
            calc_with_cartesian(&engine, &unreachable, f64::NAN, MARS, flags),
            Err(Error::NonFiniteArgument { .. }),
        ));
    }

    #[test]
    fn flag_normalization_battery() {
        // Empty request defaults to SWIEPH.
        assert_eq!(normalize_calc_flags(0, false), FLG_SWIEPH);
        // MOSEPH stripped, JPLEPH wins over SWIEPH, SPEED3 implies SPEED.
        let norm = normalize_calc_flags(FLG_MOSEPH | FLG_SPEED3, false);
        assert_eq!(norm & FLG_MOSEPH, 0);
        assert_eq!(norm & (FLG_SPEED | FLG_SPEED3), FLG_SPEED | FLG_SPEED3);
        // Explicit SPEED clears SPEED3.
        let norm = normalize_calc_flags(FLG_SPEED | FLG_SPEED3 | FLG_SWIEPH, false);
        assert_eq!(norm & FLG_SPEED3, 0);
        assert_ne!(norm & FLG_SPEED, 0);
        // JPLEPH priority keeps bits exclusive.
        let norm = normalize_calc_flags(FLG_JPLEPH | FLG_SWIEPH, false);
        assert_eq!(norm & (FLG_JPLEPH | FLG_SWIEPH), FLG_JPLEPH);
        // JPL fallback degrades the JPLEPH echo to SWIEPH.
        let norm = normalize_calc_flags(FLG_JPLEPH, true);
        assert_eq!(norm & (FLG_JPLEPH | FLG_SWIEPH), FLG_SWIEPH);
        // Center priority TOPOCTR > BARYCTR > HELCTR.
        let norm = normalize_calc_flags(FLG_HELCTR | FLG_BARYCTR | FLG_TOPOCTR, false);
        assert_eq!(norm & (FLG_HELCTR | FLG_BARYCTR | FLG_TOPOCTR), FLG_TOPOCTR);
        let norm = normalize_calc_flags(FLG_HELCTR | FLG_BARYCTR, false);
        assert_eq!(norm & (FLG_HELCTR | FLG_BARYCTR), FLG_BARYCTR);
    }

    #[test]
    fn consume_non_echoed_flags_battery() {
        // CENTER_BODY consumed for Sun..Mars, kept for nodes.
        assert_eq!(consume_non_echoed_flags(FLG_CENTER_BODY, MARS), 0);
        assert_eq!(
            consume_non_echoed_flags(FLG_CENTER_BODY, MEAN_NODE),
            FLG_CENTER_BODY
        );
        // JPLHOR bits always consumed.
        assert_eq!(consume_non_echoed_flags(FLG_JPLHOR, MEAN_NODE), 0);
        assert_eq!(
            consume_non_echoed_flags(FLG_JPLHOR_APPROX | FLG_SPEED, MOON),
            FLG_SPEED
        );
    }

    #[test]
    fn echo_and_implied_bits() {
        // SPEED3-only echo restores SPEED3, clears implied SPEED.
        assert_eq!(
            echo_request_bits(FLG_SPEED | FLG_SWIEPH, FLG_SPEED3),
            FLG_SPEED3 | FLG_SWIEPH
        );
        // Explicit SPEED keeps SPEED.
        assert_eq!(
            echo_request_bits(FLG_SPEED | FLG_SWIEPH, FLG_SPEED | FLG_SPEED3),
            FLG_SPEED | FLG_SWIEPH
        );
        // MOSEPH-only echo restores MOSEPH over the SWIEPH default.
        assert_eq!(echo_request_bits(FLG_SWIEPH, FLG_MOSEPH), FLG_MOSEPH);
        // Implied bits: J2000/SIDEREAL -> NONUT; helio/bary/truepos ->
        // NOGDEFL|NOABERR. Computation flags must not gain them.
        assert_eq!(implied_retflag_bits(FLG_J2000), FLG_NONUT);
        assert_eq!(implied_retflag_bits(FLG_SIDEREAL), FLG_NONUT);
        assert_eq!(implied_retflag_bits(FLG_HELCTR), FLG_NOGDEFL | FLG_NOABERR);
        assert_eq!(
            implied_retflag_bits(FLG_BARYCTR | FLG_TRUEPOS),
            FLG_NOGDEFL | FLG_NOABERR
        );
        assert_eq!(implied_retflag_bits(FLG_SPEED), 0);
        // TT echo strips the defaulted ephemeris bit only when the raw
        // request carried none.
        assert_eq!(
            calc_tt_epheflag_echo(FLG_SWIEPH | FLG_SPEED, FLG_SPEED),
            FLG_SPEED
        );
        assert_eq!(
            calc_tt_epheflag_echo(FLG_SWIEPH | FLG_SPEED, FLG_SWIEPH | FLG_SPEED),
            FLG_SWIEPH | FLG_SPEED
        );
        // Exclusive bit keeps MOSEPH (pctr/ECL_NUT convention).
        assert_eq!(
            exclusive_ephemeris_bit(FLG_MOSEPH, false) & FLG_MOSEPH,
            FLG_MOSEPH
        );
        assert_eq!(exclusive_ephemeris_bit(0, false) & FLG_SWIEPH, FLG_SWIEPH);
    }

    #[test]
    fn barycentric_moseph_refusal() {
        // BARYCTR + MOSEPH-only is refused, naming the entry point.
        let err = validate_barycentric_moseph(FLG_BARYCTR | FLG_MOSEPH, "calc_ut")
            .expect_err("must refuse");
        assert_eq!(
            err,
            Error::UnsupportedFlagCombination {
                entry: "calc_ut",
                flags: FLG_BARYCTR | FLG_MOSEPH,
            }
        );
        // TOPOCTR removes BARYCTR: allowed. JPLEPH/SWIEPH wins over
        // MOSEPH: allowed.
        assert!(
            validate_barycentric_moseph(FLG_BARYCTR | FLG_TOPOCTR | FLG_MOSEPH, "calc_ut").is_ok()
        );
        assert!(validate_barycentric_moseph(FLG_BARYCTR | FLG_MOSEPH | FLG_SWIEPH, "calc").is_ok());
        assert!(validate_barycentric_moseph(FLG_SPEED, "calc").is_ok());
    }

    #[test]
    fn output_shaping_xyz_priority_and_radians() {
        let sph = [180.0, 30.0, 2.0, 1.0, 0.5, 0.1];
        let xyz = apply_output_flags(sph, FLG_XYZ | FLG_RADIANS);
        // |pos| == dist and XYZ wins over RADIANS.
        let norm = (xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]).sqrt();
        assert!((norm - 2.0).abs() < 1e-12);
        let rad = apply_output_flags(sph, FLG_RADIANS);
        assert!((rad[0] - std::f64::consts::PI).abs() < 1e-15);
        assert_eq!(rad[2], 2.0);
        assert_eq!(rad[5], 0.1);
        // Zero rates give exactly zero velocity (verbatim branch).
        let still = spherical_to_cartesian_with_velocity(10.0, 20.0, 1.5, 0.0, 0.0, 0.0);
        assert_eq!([still[3], still[4], still[5]], [0.0, 0.0, 0.0]);
        // finalize echoes format bits plus implied bits.
        let (_, ret) = finalize_output_flags(sph, FLG_SPEED, FLG_SPEED | FLG_XYZ | FLG_J2000);
        assert_eq!(ret & (FLG_XYZ | FLG_NONUT), FLG_XYZ | FLG_NONUT);
    }

    #[test]
    fn south_node_antipode_in_three_representations() {
        let north = [100.0, 20.0, 1.5, 0.5, 0.1, 0.01];
        let south = south_node_from_north(north, FLG_SPEED);
        assert_eq!(south, [280.0, -20.0, 1.5, 0.5, -0.1, 0.01]);
        // Wrap across 0.
        let south = south_node_from_north([350.0, -5.0, 1.0, 0.0, 0.0, 0.0], FLG_SPEED);
        assert!((south[0] - 170.0).abs() < 1e-12);
        // Radians: half-turn is pi.
        let north_r = [1.0, 0.2, 1.5, 0.01, 0.02, 0.0];
        let south_r = south_node_from_north(north_r, FLG_RADIANS);
        assert!((south_r[0] - (1.0 + std::f64::consts::PI)).abs() < 1e-15);
        assert_eq!(south_r[1], -0.2);
        // Cartesian: full negation.
        let south_c = south_node_from_north([1.0, 2.0, 3.0, 0.1, 0.2, 0.3], FLG_XYZ);
        assert_eq!(south_c, [-1.0, -2.0, -3.0, -0.1, -0.2, -0.3]);
    }

    #[test]
    fn speed3_stencil_recovers_analytic_rate() {
        // Linear drift 2.5 deg/day around lon 359.9 exercises the wrap
        // unwrapping; the centered difference must recover the rate.
        let current = [359.9, 10.0, 1.0, 0.0, 0.0, 0.0];
        let previous = [359.9 - 2.5 * SPEED3_STEP_DAYS, 10.0, 1.0, 0.0, 0.0, 0.0];
        let following = [359.9 + 2.5 * SPEED3_STEP_DAYS, 10.0, 1.0, 0.0, 0.0, 0.0];
        let out = apply_central_speed3_stencil(current, previous, following, FLG_SPEED3);
        assert!((out[3] - 2.5).abs() < 1e-9);
        assert_eq!([out[0], out[1], out[2]], [359.9, 10.0, 1.0]);
        // Non-mean-point bodies pass through the rebuild unchanged.
        let same =
            rebuild_mean_point_speed3(current, MARS, FLG_SPEED3, |_| panic!("must not sample"))
                .expect("pass-through");
        assert_eq!(same, current);
        // Mean node rebuilds through the sampler (constant field ->
        // zero rate, position preserved).
        let rebuilt =
            rebuild_mean_point_speed3(current, MEAN_NODE, FLG_SPEED3 | FLG_SPEED, |_| Ok(current))
                .expect("rebuild");
        assert_eq!([rebuilt[0], rebuilt[1], rebuilt[2]], [359.9, 10.0, 1.0]);
        assert_eq!([rebuilt[3], rebuilt[4], rebuilt[5]], [0.0, 0.0, 0.0]);
    }

    #[test]
    fn calc_ut_entry_tail_over_stub_provider() {
        let (ut1_to_tt, compute, nut_ut, nut_tt) = stub_engine();
        let engine = engine_with(&ut1_to_tt, &compute, &nut_ut, &nut_tt);
        // Plain request: SWIEPH default echoed, position from provider
        // at the converted epoch.
        let (pos, ret) = calc_ut(&engine, 2451545.0, MARS, FLG_SPEED).expect("calc_ut");
        assert_eq!(pos[0], 100.0 + MARS as f64);
        assert_eq!(ret & FLG_SWIEPH, FLG_SWIEPH);
        assert_ne!(ret & FLG_SPEED, 0);
        // calc() with no ephemeris bit echoes none (TT convention).
        let (_, ret) = calc(&engine, 2451545.0, MARS, FLG_SPEED).expect("calc");
        assert_eq!(ret & (FLG_JPLEPH | FLG_SWIEPH | FLG_MOSEPH), 0);
        // J2000 echoes NONUT without changing the computation flags
        // the provider saw (provider echoes its input flags back).
        let (_, ret) = calc(&engine, 2451545.0, MARS, FLG_SPEED | FLG_J2000).expect("calc");
        assert_eq!(ret & FLG_NONUT, FLG_NONUT);
        // Degenerate topocentric Earth: exact zeros, TOPOCTR echoed.
        let (pos, ret) =
            calc(&engine, 2451545.0, EARTH, FLG_SPEED | FLG_TOPOCTR).expect("degenerate");
        assert_eq!(pos, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        assert_ne!(ret & FLG_TOPOCTR, 0);
        // Barycentric MOSEPH refused before any provider call.
        assert!(calc(&engine, 2451545.0, MARS, FLG_BARYCTR | FLG_MOSEPH).is_err());
        // Non-finite epoch refused.
        assert!(calc(&engine, f64::NAN, MARS, FLG_SPEED).is_err());
        assert!(calc_ut(&engine, f64::INFINITY, MARS, FLG_SPEED).is_err());
    }

    #[test]
    fn south_node_and_ecl_nut_branches() {
        let (ut1_to_tt, compute, nut_ut, nut_tt) = stub_engine();
        let engine = engine_with(&ut1_to_tt, &compute, &nut_ut, &nut_tt);
        // South node recurses through the north node then antipodes.
        let (pos, _) = calc(&engine, 2451545.0, -MEAN_NODE, FLG_SPEED).expect("south");
        assert!((pos[0] - (100.0 + MEAN_NODE as f64 + 180.0).rem_euclid(360.0)).abs() < 1e-9);
        assert_eq!(pos[1], -10.0);
        let (pos, _) = calc_ut(&engine, 2451545.0, -TRUE_NODE, FLG_SPEED).expect("south ut");
        assert!((pos[0] - (100.0 + TRUE_NODE as f64 + 180.0).rem_euclid(360.0)).abs() < 1e-9);
        // ECL_NUT returns the provider tuple with exclusive-bit echo:
        // MOSEPH-only request echoes MOSEPH, keeps no SWIEPH.
        let (pos, ret) = calc(&engine, 2451545.0, ECL_NUT, FLG_MOSEPH).expect("ecl_nut");
        assert_eq!(pos[0], 23.44);
        assert_eq!(ret & FLG_MOSEPH, FLG_MOSEPH);
        assert_eq!(ret & FLG_SWIEPH, 0);
        // calc_ut ECL_NUT with no ephemeris bit gains the SWIEPH default.
        let (_, ret) = calc_ut(&engine, 2451545.0, ECL_NUT, 0).expect("ecl_nut ut");
        assert_eq!(ret & FLG_SWIEPH, FLG_SWIEPH);
    }

    #[test]
    fn stub_orchestration_records_no_trace() {
        // Recording happens at the sealed serving legs, never in the
        // provider-agnostic orchestration: even with an active
        // session, stub-engine evaluations — including the south-node
        // transfer arm — leave the map empty.
        let (ut1_to_tt, compute, nut_ut, nut_tt) = stub_engine();
        let engine = engine_with(&ut1_to_tt, &compute, &nut_ut, &nut_tt);
        let token = crate::trace::start_tracing();
        let _ = calc_ut(&engine, 2451545.0, MARS, FLG_SPEED).expect("calc_ut");
        let _ = calc_ut(&engine, 2451545.0, -TRUE_NODE, FLG_SPEED).expect("south ut");
        let _ = calc(&engine, 2451545.0, -MEAN_NODE, FLG_SPEED).expect("south");
        assert!(crate::trace::get_trace_results().is_empty());
        assert!(crate::trace::reset_tracing(token));
    }

    #[test]
    fn pctr_entry_contracts() {
        // Target at [2, 0, 0] AU, observer at origin: TRUEPOS gives
        // lon 0, dist 2 with zero light time behind it.
        let state_at = |which: PctrEndpoint, _jd: f64| {
            Ok(match which {
                PctrEndpoint::Target => ([2.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
                PctrEndpoint::Observer => ([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
            })
        };
        // Ecliptic-of-date frame stub: identity rotation, spherical shaping.
        let frame_at = |p: [f64; 3], _jd: f64, _flags: i32| {
            cartesian_to_spherical_deg(p).map(|s| (s[0], s[1], s[2]))
        };
        let engine = PctrEngine {
            state_at: &state_at,
            apply_deflection: &|p, _, _, _| Ok(p),
            apply_aberration: &|p, _, _| Ok(p),
            frame_at: &frame_at,
            ayanamsha_deg: None,
            jpl_fallback_active: false,
        };
        // TRUEPOS vector core: exact difference, zero light time.
        let (p, lt) =
            pctr_apparent_vector(&engine, 2451545.0, FLG_TRUEPOS | FLG_SPEED).expect("truepos");
        assert_eq!(p, [2.0, 0.0, 0.0]);
        assert_eq!(lt, 0.0);
        // Light-time loop converges to dist/c for a static target.
        let (p, lt) = pctr_apparent_vector(&engine, 2451545.0, FLG_SPEED).expect("light time");
        assert_eq!(p, [2.0, 0.0, 0.0]);
        assert!((lt - 2.0 / C_AU_PER_DAY_PCTR).abs() < 1e-15);
        // Full entry: lon 0, dist 2, zero speeds without SPEED.
        let (pos, ret) = calc_pctr(&engine, 2451545.0, MARS, EARTH, FLG_SWIEPH).expect("pctr");
        assert!((pos[0] - 0.0).abs() < 1e-12);
        assert!((pos[2] - 2.0).abs() < 1e-12);
        assert_eq!([pos[3], pos[4], pos[5]], [0.0, 0.0, 0.0]);
        assert_eq!(ret & FLG_SWIEPH, FLG_SWIEPH);
        // SPEED3 alone is echoed but leaves velocity slots zero.
        let (pos, ret) = calc_pctr(&engine, 2451545.0, MARS, EARTH, FLG_SPEED3).expect("speed3");
        assert_eq!([pos[3], pos[4], pos[5]], [0.0, 0.0, 0.0]);
        assert_ne!(ret & FLG_SPEED3, 0);
        // HELCTR maps to the astrometric place: center bit stripped,
        // NOABERR|NOGDEFL echoed.
        let (_, ret) =
            calc_pctr(&engine, 2451545.0, MARS, EARTH, FLG_HELCTR | FLG_SPEED).expect("helctr");
        assert_eq!(ret & (FLG_HELCTR | FLG_BARYCTR), 0);
        assert_eq!(ret & (FLG_NOABERR | FLG_NOGDEFL), FLG_NOABERR | FLG_NOGDEFL);
        // Unsupported endpoints refused without touching providers.
        assert!(calc_pctr(&engine, 2451545.0, MEAN_NODE, EARTH, FLG_SPEED).is_err());
        assert!(calc_pctr(&engine, 2451545.0, MARS, MEAN_NODE, FLG_SPEED).is_err());
        assert!(calc_pctr(&engine, f64::NAN, MARS, EARTH, FLG_SPEED).is_err());
        // Sidereal without a provider refused (never silent tropical).
        assert!(calc_pctr(&engine, 2451545.0, MARS, EARTH, FLG_SIDEREAL | FLG_SPEED).is_err());
    }

    #[test]
    fn pctr_speed_stencil_wraps_and_rates() {
        // Rotating target: angle advances 1 deg/day; the stencil must
        // recover the rate and preserve the wrap fold.
        let state_at = |which: PctrEndpoint, jd: f64| {
            let angle = (jd * 1.0f64).to_radians();
            Ok(match which {
                PctrEndpoint::Target => ([angle.cos() * 2.0, angle.sin() * 2.0, 0.0], [0.0; 3]),
                PctrEndpoint::Observer => ([0.0; 3], [0.0; 3]),
            })
        };
        let frame_at = |p: [f64; 3], _jd: f64, _flags: i32| {
            cartesian_to_spherical_deg(p).map(|s| (s[0], s[1], s[2]))
        };
        let engine = PctrEngine {
            state_at: &state_at,
            apply_deflection: &|p, _, _, _| Ok(p),
            apply_aberration: &|p, _, _| Ok(p),
            frame_at: &frame_at,
            ayanamsha_deg: None,
            jpl_fallback_active: false,
        };

        let (pos, _) =
            calc_pctr(&engine, 100.0, MARS, EARTH, FLG_TRUEPOS | FLG_SPEED).expect("pctr speed");
        assert!((pos[3] - 1.0).abs() < 1e-9, "dlon = {}", pos[3]);
        assert!(pos[4].abs() < 1e-6, "dlat = {}", pos[4]);
        assert!(pos[5].abs() < 1e-6, "ddist = {}", pos[5]);
        // Half-step table: Moon six seconds, others one second.
        assert_eq!(pctr_speed_half_step_days(MOON), 6.0 / 86400.0);
        assert_eq!(pctr_speed_half_step_days(MARS), 1.0 / 86400.0);
        // pctr body set: planets plus Earth, never nodes.
        assert!(is_pctr_supported_body(EARTH));
        assert!(!is_pctr_supported_body(MEAN_NODE));
        assert!(!is_pctr_supported_body(99));
    }
}

#[test]
fn pos10_dispatch_table_predicates_and_names() {
    use crate::constants::{
        ADMETOS, ALGOL, ANGLE_OFFSET, APOLLON, ARABIC_OFFSET, CERES, CUPIDO, HADES, HARRINGTON,
        ISIS, JUNO, KRONOS, MERCURY, MOON_CHARON, MOON_TITAN, NIBIRU, PALLAS, PHOLUS,
        PLMOON_OFFSET, POSEIDON, PROSERPINA, REGULUS, SIRIUS, SPICA_STAR, VESTA, VULCAN, VULKANUS,
        WALDEMATH, WHITE_MOON, ZEUS,
    };
    // AST_OFFSET remap: the six built-in numbers resolve, the rest
    // pass through (verbatim `_AST_NUMBER_REMAP`).
    assert_eq!(remap_ast_offset(AST_OFFSET + 1), CERES);
    assert_eq!(remap_ast_offset(AST_OFFSET + 2), PALLAS);
    assert_eq!(remap_ast_offset(AST_OFFSET + 3), JUNO);
    assert_eq!(remap_ast_offset(AST_OFFSET + 4), VESTA);
    assert_eq!(remap_ast_offset(AST_OFFSET + 2060), CHIRON);
    assert_eq!(remap_ast_offset(AST_OFFSET + 5145), PHOLUS);
    assert_eq!(remap_ast_offset(AST_OFFSET + 5), AST_OFFSET + 5);
    assert_eq!(remap_ast_offset(CERES), CERES);
    assert_eq!(remap_ast_offset(MARS), MARS);
    // JPL-ephemeris predicate: exactly the `_PLANET_MAP` set.
    for id in [
        SUN, MOON, MERCURY, VENUS, MARS, JUPITER, SATURN, URANUS, NEPTUNE, PLUTO, EARTH,
    ] {
        assert!(body_uses_jpl_ephemeris(id), "jpl {id}");
    }
    for id in [MEAN_NODE, MEAN_APOG, CERES, CUPIDO, AST_OFFSET + 1] {
        assert!(!body_uses_jpl_ephemeris(id), "non-jpl {id}");
    }
    // Planetary-moon predicate: legacy aliases plus the canonical
    // range endpoints; neighbors of the range are not moons.
    assert!(is_planetary_moon_id(MOON_TITAN));
    assert!(is_planetary_moon_id(MOON_CHARON));
    assert!(is_planetary_moon_id(PLMOON_OFFSET + 501));
    assert!(is_planetary_moon_id(PLMOON_OFFSET + 999));
    assert!(!is_planetary_moon_id(PLMOON_OFFSET + 400));
    assert!(!is_planetary_moon_id(PLMOON_OFFSET + 1000));
    assert!(!is_planetary_moon_id(SATURN));
    // Names: verbatim table spot checks plus the number branch.
    assert_eq!(planet_name(SUN), "Sun");
    assert_eq!(planet_name(MEAN_NODE), "mean Node");
    assert_eq!(planet_name(OSCU_APOG), "osc. Apogee");
    assert_eq!(planet_name(VULKANUS), "Vulcanus");
    assert_eq!(planet_name(WHITE_MOON), "Selena/White Moon");
    assert_eq!(planet_name(NEPTUNE_LEVERRIER), "Leverrier (Neptune)");
    assert_eq!(planet_name(CERES), "Ceres");
    assert_eq!(planet_name(AST_OFFSET + 1), "Ceres");
    assert_eq!(planet_name(AST_OFFSET + 2060), "Chiron");
    assert_eq!(planet_name(99999), "");
    assert_eq!(planet_name(-MEAN_NODE), "");

    assert_eq!(classify_body(MARS), BodyKind::Planet);
    assert_eq!(classify_body(EARTH), BodyKind::Planet);
    assert_eq!(classify_body(MEAN_NODE), BodyKind::LunarNode);
    assert_eq!(classify_body(TRUE_NODE), BodyKind::LunarNode);
    assert_eq!(classify_body(MEAN_APOG), BodyKind::LunarApse);
    assert_eq!(classify_body(OSCU_APOG), BodyKind::LunarApse);
    assert_eq!(classify_body(INTP_APOG), BodyKind::InterpolatedApse);
    assert_eq!(classify_body(INTP_PERG), BodyKind::InterpolatedApse);
    assert_eq!(classify_body(CUPIDO), BodyKind::Uranian);
    assert_eq!(classify_body(POSEIDON), BodyKind::Uranian);
    assert_eq!(classify_body(HADES), BodyKind::Uranian);
    assert_eq!(classify_body(ZEUS), BodyKind::Uranian);
    assert_eq!(classify_body(KRONOS), BodyKind::Uranian);
    assert_eq!(classify_body(APOLLON), BodyKind::Uranian);
    assert_eq!(classify_body(ADMETOS), BodyKind::Uranian);
    assert_eq!(classify_body(ISIS), BodyKind::Transpluto);
    for id in [
        WHITE_MOON,
        VULCAN,
        PROSERPINA,
        WALDEMATH,
        NIBIRU,
        HARRINGTON,
        NEPTUNE_ADAMS,
        PLUTO_LOWELL,
        PLUTO_PICKERING,
    ] {
        assert_eq!(classify_body(id), BodyKind::FictHelio, "fict {id}");
    }
    assert_eq!(classify_body(MOON_TITAN), BodyKind::PlanetaryMoon);
    assert_eq!(classify_body(PLMOON_OFFSET + 501), BodyKind::PlanetaryMoon);
    assert_eq!(classify_body(CERES), BodyKind::MinorBody);
    assert_eq!(classify_body(PHOLUS), BodyKind::MinorBody);
    assert_eq!(classify_body(JUNO), BodyKind::MinorBody);
    assert_eq!(classify_body(VESTA), BodyKind::MinorBody);
    assert_eq!(classify_body(AST_OFFSET + 5), BodyKind::MinorBody);
    assert_eq!(classify_body(REGULUS), BodyKind::FixedStar);
    assert_eq!(classify_body(SPICA_STAR), BodyKind::FixedStar);
    assert_eq!(classify_body(ALGOL), BodyKind::FixedStar);
    assert_eq!(classify_body(SIRIUS), BodyKind::FixedStar);
    assert_eq!(classify_body(ANGLE_OFFSET + 1), BodyKind::Angle);
    assert_eq!(classify_body(ARABIC_OFFSET + 1), BodyKind::ArabicPart);
    assert_eq!(classify_body(ARABIC_OFFSET + 99), BodyKind::ArabicPart);
    assert_eq!(classify_body(ARABIC_OFFSET + 100), BodyKind::Unknown);

    assert_eq!(classify_body(99999), BodyKind::MinorBody);
    assert_eq!(classify_body(23), BodyKind::Unknown);
    // Every known family has an owning task; unknown has none.
    for kind in [
        BodyKind::Planet,
        BodyKind::LunarNode,
        BodyKind::LunarApse,
        BodyKind::InterpolatedApse,
        BodyKind::Uranian,
        BodyKind::Transpluto,
        BodyKind::FictHelio,
        BodyKind::PlanetaryMoon,
        BodyKind::MinorBody,
        BodyKind::FixedStar,
        BodyKind::Angle,
        BodyKind::ArabicPart,
    ] {
        assert!(!body_kind_owner(kind).starts_with("none"), "{kind:?}");
    }
    // Planet kind and the JPL predicate agree exactly.
    for id in -2..3_100_000 {
        assert_eq!(
            classify_body(id) == BodyKind::Planet,
            body_uses_jpl_ephemeris(id),
            "planet/predicate agreement at {id}"
        );
    }
}

#[test]
fn pos10_entry_remap_reaches_provider() {
    // `calc_ut` remaps AST_OFFSET+N before dispatch: both spellings
    // reach `compute_tt` as the dedicated body id.
    let seen = std::cell::RefCell::new(Vec::new());
    let compute_tt = |_jd: f64, body: i32, _flags: i32| {
        seen.borrow_mut().push(body);
        Ok(([1.0, 2.0, 3.0, 0.0, 0.0, 0.0], 0))
    };
    let engine = CalcEngine {
        ut1_to_tt_jd: &|jd| Ok(jd),
        compute_tt: &compute_tt,
        nutation_obliquity_ut: &|_, _| unreachable!(),
        nutation_obliquity_tt: &|_, _| unreachable!(),
        jpl_fallback_active: false,
    };
    let (pos_alias, _) =
        calc_ut(&engine, 2451545.0, AST_OFFSET + 1, FLG_SWIEPH).expect("alias spelling");
    let (pos_direct, _) = calc_ut(&engine, 2451545.0, CERES, FLG_SWIEPH).expect("direct spelling");
    assert_eq!(pos_alias, pos_direct);
    assert_eq!(
        *seen.borrow(),
        vec![CERES, CERES],
        "both spellings dispatch as CERES"
    );
}

mod batch_tests {
    use super::super::{
        CalcEngine, CalcRequest, PctrEngine, PctrRequest, calc, calc_many, calc_pctr,
        calc_pctr_many, calc_pctr_requests, calc_requests, calc_ut, calc_ut_many, calc_ut_requests,
    };
    use crate::constants::{CERES, ECL_NUT, FLG_SPEED, FLG_SWIEPH, MARS, MEAN_NODE, MOON, VENUS};
    use crate::errors::Error;

    #[allow(clippy::type_complexity)]
    fn calc_stub_engine() -> (
        Box<dyn Fn(f64) -> Result<f64, Error>>,
        Box<dyn Fn(f64, i32, i32) -> Result<(super::super::PositionResult, i32), Error>>,
        Box<dyn Fn(f64, i32) -> Result<(super::super::PositionResult, i32), Error>>,
        Box<dyn Fn(f64, i32) -> Result<(super::super::PositionResult, i32), Error>>,
    ) {
        // Epoch- and body-sensitive analytic core so batch/scalar
        // divergence cannot hide behind a constant: lon encodes both.
        let compute = Box::new(|jd: f64, body: i32, flags: i32| {
            Ok(([jd * 0.01 + body as f64, 10.0, 1.5, 0.5, 0.0, 0.0], flags))
        })
            as Box<dyn Fn(f64, i32, i32) -> Result<(super::super::PositionResult, i32), Error>>;
        let nut =
            Box::new(|_jd: f64, flags: i32| Ok(([23.44, 23.43, 0.001, 0.002, 0.0, 0.0], flags)))
                as Box<dyn Fn(f64, i32) -> Result<(super::super::PositionResult, i32), Error>>;
        let nut2 =
            Box::new(|_jd: f64, flags: i32| Ok(([23.44, 23.43, 0.001, 0.002, 0.0, 0.0], flags)))
                as Box<dyn Fn(f64, i32) -> Result<(super::super::PositionResult, i32), Error>>;
        let time = Box::new(|jd: f64| Ok(jd + 0.0008)) as Box<dyn Fn(f64) -> Result<f64, Error>>;
        (time, compute, nut, nut2)
    }

    fn engine_with<'a>(
        time: &'a dyn Fn(f64) -> Result<f64, Error>,
        compute: &'a dyn Fn(f64, i32, i32) -> Result<(super::super::PositionResult, i32), Error>,
        nut_ut: &'a dyn Fn(f64, i32) -> Result<(super::super::PositionResult, i32), Error>,
        nut_tt: &'a dyn Fn(f64, i32) -> Result<(super::super::PositionResult, i32), Error>,
    ) -> CalcEngine<'a> {
        CalcEngine {
            ut1_to_tt_jd: time,
            compute_tt: compute,
            nutation_obliquity_ut: nut_ut,
            nutation_obliquity_tt: nut_tt,
            jpl_fallback_active: false,
        }
    }

    #[test]
    fn homogeneous_ut_batch_matches_repeated_scalars() {
        let (time, compute, nut_ut, nut_tt) = calc_stub_engine();
        let engine = engine_with(&time, &compute, &nut_ut, &nut_tt);
        let jds = [2451545.0, 2451545.5, 2451546.0, 2460000.25];
        let batched = calc_ut_many(&engine, &jds, MARS, FLG_SWIEPH | FLG_SPEED).expect("batch");
        assert_eq!(batched.len(), jds.len());
        for (jd, row) in jds.iter().zip(batched.iter()) {
            let scalar = calc_ut(&engine, *jd, MARS, FLG_SWIEPH | FLG_SPEED).expect("scalar");
            assert_eq!(*row, scalar, "batch diverged from scalar at jd {jd}");
        }
        // Repeated batch calls are deterministic.
        let again = calc_ut_many(&engine, &jds, MARS, FLG_SWIEPH | FLG_SPEED).expect("batch");
        assert_eq!(batched, again);
        // Empty batch returns an empty vector, never an error.
        let empty: [f64; 0] = [];
        assert!(
            calc_ut_many(&engine, &empty, MARS, FLG_SWIEPH)
                .expect("empty")
                .is_empty()
        );
    }

    #[test]
    fn homogeneous_tt_batch_matches_repeated_scalars() {
        let (time, compute, nut_ut, nut_tt) = calc_stub_engine();
        let engine = engine_with(&time, &compute, &nut_ut, &nut_tt);
        let jds = [2451545.0, 2451545.5, 2451546.0];
        let batched = calc_many(&engine, &jds, VENUS, FLG_SWIEPH).expect("batch");
        for (jd, row) in jds.iter().zip(batched.iter()) {
            let scalar = calc(&engine, *jd, VENUS, FLG_SWIEPH).expect("scalar");
            assert_eq!(*row, scalar, "batch diverged from scalar at jd {jd}");
        }
    }

    #[test]
    fn heterogeneous_batch_covers_branches_per_request() {
        let (time, compute, nut_ut, nut_tt) = calc_stub_engine();
        let engine = engine_with(&time, &compute, &nut_ut, &nut_tt);
        // Mixed epochs, bodies, flags, alias spelling, south node,
        // and the ECL_NUT raw-flag branch — per-request semantics.
        let requests = [
            CalcRequest {
                jd: 2451545.0,
                body_id: MARS,
                flags: FLG_SWIEPH | FLG_SPEED,
            },
            CalcRequest {
                jd: 2451545.5,
                body_id: VENUS,
                flags: FLG_SWIEPH,
            },
            CalcRequest {
                jd: 2451546.0,
                body_id: super::super::AST_OFFSET + 1,
                flags: FLG_SWIEPH | FLG_SPEED,
            },
            CalcRequest {
                jd: 2451546.5,
                body_id: -MEAN_NODE,
                flags: FLG_SWIEPH | FLG_SPEED,
            },
            CalcRequest {
                jd: 2451547.0,
                body_id: ECL_NUT,
                flags: FLG_SWIEPH | FLG_SPEED,
            },
        ];
        let batched = calc_ut_requests(&engine, &requests).expect("batch");
        assert_eq!(batched.len(), requests.len());
        for (request, row) in requests.iter().zip(batched.iter()) {
            let scalar =
                calc_ut(&engine, request.jd, request.body_id, request.flags).expect("scalar");
            assert_eq!(
                *row, scalar,
                "heterogeneous batch diverged at {:?}",
                request
            );
        }
        // Alias spelling agrees with the direct id inside the batch.
        let alias = calc_ut_requests(
            &engine,
            &[CalcRequest {
                jd: 2451545.0,
                body_id: super::super::AST_OFFSET + 1,
                flags: FLG_SWIEPH,
            }],
        )
        .expect("alias batch");
        let direct = calc_ut_requests(
            &engine,
            &[CalcRequest {
                jd: 2451545.0,
                body_id: CERES,
                flags: FLG_SWIEPH,
            }],
        )
        .expect("direct batch");
        assert_eq!(alias, direct);
        // TT heterogeneous path carries the same per-request contract.
        let tt_batched = calc_requests(&engine, &requests).expect("tt batch");
        for (request, row) in requests.iter().zip(tt_batched.iter()) {
            let scalar = calc(&engine, request.jd, request.body_id, request.flags).expect("scalar");
            assert_eq!(*row, scalar, "TT batch diverged at {request:?}");
        }
        // Empty heterogeneous batch is Ok(empty).
        assert!(calc_ut_requests(&engine, &[]).expect("empty").is_empty());
    }

    #[test]
    fn batch_is_fail_fast_with_scalar_error_identity() {
        // Provider fails on its second invocation (the middle epoch);
        // the batch surfaces that same error (fail-fast) rather than
        // a synthesized batch error. A call counter is used because
        // `compute_tt` sees TT epochs (UT + 0.0008), so a UT-epoch
        // equality gate would never fire.
        let time = |jd: f64| Ok(jd + 0.0008);
        let calls = std::cell::Cell::new(0u32);
        let compute = |jd: f64, body: i32, flags: i32| {
            calls.set(calls.get() + 1);
            if calls.get() == 2 {
                return Err(Error::NonFiniteArgument { value: jd });
            }
            Ok(([jd * 0.01 + body as f64, 10.0, 1.5, 0.0, 0.0, 0.0], flags))
        };
        let nut = |_jd: f64, flags: i32| Ok(([23.44, 23.43, 0.001, 0.002, 0.0, 0.0], flags));
        let engine = CalcEngine {
            ut1_to_tt_jd: &time,
            compute_tt: &compute,
            nutation_obliquity_ut: &nut,
            nutation_obliquity_tt: &nut,
            jpl_fallback_active: false,
        };
        let jds = [2451545.0, 2451545.5, 2451546.0];
        let batch_err = calc_ut_many(&engine, &jds, MARS, FLG_SWIEPH).expect_err("must fail");
        // Re-arm the counter so the scalar comparison fails on the
        // same (second) invocation with the same TT epoch.
        calls.set(1);
        let scalar_err =
            calc_ut(&engine, 2451545.5, MARS, FLG_SWIEPH).expect_err("scalar must fail");
        assert_eq!(format!("{batch_err:?}"), format!("{scalar_err:?}"));
        // Heterogeneous fail-fast likewise mirrors the scalar error.
        let requests = [
            CalcRequest {
                jd: 2451545.0,
                body_id: MARS,
                flags: FLG_SWIEPH,
            },
            CalcRequest {
                jd: 2451545.5,
                body_id: MARS,
                flags: FLG_SWIEPH,
            },
        ];
        calls.set(0);
        let het_err = calc_ut_requests(&engine, &requests).expect_err("must fail");
        assert_eq!(format!("{het_err:?}"), format!("{scalar_err:?}"));
        // Non-finite epoch input fails before any provider call.
        let bad = calc_ut_many(&engine, &[f64::NAN], MARS, FLG_SWIEPH).expect_err("NaN");
        assert!(matches!(bad, Error::NonFiniteArgument { .. }));
    }

    #[test]
    fn pctr_batches_match_repeated_scalars() {
        use super::super::PctrEndpoint;
        let state_at = |endpoint: PctrEndpoint, jd: f64| {
            let base = match endpoint {
                PctrEndpoint::Target => 1.0 + jd * 1e-6,
                PctrEndpoint::Observer => 0.1,
            };
            Ok(([base, 0.2, 0.3], [0.01, 0.0, 0.0]))
        };
        let deflect = |p: [f64; 3], _o: [f64; 3], _jd: f64, _lt: f64| Ok(p);
        let aberrate = |p: [f64; 3], _v: [f64; 3], _lt: f64| Ok(p);
        let frame_at = |_p: [f64; 3], _jd: f64, _flags: i32| Ok((100.0, 10.0, 1.5));
        let engine = PctrEngine {
            state_at: &state_at,
            apply_deflection: &deflect,
            apply_aberration: &aberrate,
            frame_at: &frame_at,
            ayanamsha_deg: None,
            jpl_fallback_active: false,
        };
        let jds = [2451545.0, 2451545.5, 2451546.0];
        let batched =
            calc_pctr_many(&engine, &jds, MARS, super::super::EARTH, FLG_SWIEPH).expect("batch");
        for (jd, row) in jds.iter().zip(batched.iter()) {
            let scalar =
                calc_pctr(&engine, *jd, MARS, super::super::EARTH, FLG_SWIEPH).expect("scalar");
            assert_eq!(*row, scalar, "pctr batch diverged at jd {jd}");
        }
        let requests = [
            PctrRequest {
                jd_tt: 2451545.0,
                target: MARS,
                center: super::super::EARTH,
                flags: FLG_SWIEPH,
            },
            PctrRequest {
                jd_tt: 2451545.5,
                target: MOON,
                center: super::super::EARTH,
                flags: FLG_SWIEPH,
            },
        ];
        let het = calc_pctr_requests(&engine, &requests).expect("het batch");
        for (request, row) in requests.iter().zip(het.iter()) {
            let scalar = calc_pctr(
                &engine,
                request.jd_tt,
                request.target,
                request.center,
                request.flags,
            )
            .expect("scalar");
            assert_eq!(*row, scalar, "pctr het diverged at {request:?}");
        }
        assert!(
            calc_pctr_many(&engine, &[], MARS, super::super::EARTH, FLG_SWIEPH)
                .expect("empty")
                .is_empty()
        );
        // Unsupported endpoint fails fast with the scalar identity.
        let bad = [PctrRequest {
            jd_tt: 2451545.0,
            target: MEAN_NODE,
            center: super::super::EARTH,
            flags: FLG_SWIEPH,
        }];
        let batch_err = calc_pctr_requests(&engine, &bad).expect_err("must fail");
        let scalar_err = calc_pctr(
            &engine,
            2451545.0,
            MEAN_NODE,
            super::super::EARTH,
            FLG_SWIEPH,
        )
        .expect_err("scalar");
        assert_eq!(format!("{batch_err:?}"), format!("{scalar_err:?}"));
    }
}

fn fixture_icrs_const(states: &[(i32, [f64; 3])]) -> Vec<u8> {
    let n = states.len();
    let mut file = pack_header(&MAGIC_LEB1, 1, 1, n as u32, 0.0, 2.0);
    let body_at = HEADER_SIZE + SECTION_DIR_SIZE;
    let blob_at = body_at + n * LEB1_BODY_ENTRY_SIZE;
    file.extend_from_slice(&pack_section(
        0,
        body_at as u64,
        (n * LEB1_BODY_ENTRY_SIZE) as u64,
    ));
    for (i, (body_id, _)) in states.iter().enumerate() {
        file.extend_from_slice(&pack_body(
            *body_id,
            COORD_ICRS_BARY,
            (blob_at + i * 96) as u64,
        ));
    }
    for (_, pos) in states {
        for _ in 0..2 {
            for component in pos {
                push_f64(&mut file, *component);
                push_f64(&mut file, 0.0);
            }
        }
    }
    file
}

fn serving_fixture() -> Vec<u8> {
    fixture_icrs_const(&[
        (SUN, [0.0, 0.0, 0.0]),
        (EARTH, [1.0, 0.0, 0.0]),
        (MARS, [1.5, 0.1, 0.0]),
    ])
}

fn zero_nutation(_jd: f64) -> Result<(f64, f64)> {
    Ok((0.0, 0.0))
}

/// Run `f` over a bound sealed serving state on the shared
/// three-body fixture (borrowing stays inside, like
/// `composite_of`).
fn with_serving<R>(f: impl FnOnce(SealedServing<'_, CompositeIndex<'_>>) -> R) -> R {
    let file = serving_fixture();
    let owned = FileIndex::decode(&file).expect("fixture must decode");
    let composite = CompositeIndex::new(vec![(owned, &file[..])]).expect("one reader must build");
    f(SealedServing::new(&composite, &zero_nutation))
}

#[test]
fn direct_xyz_serving_differences_cartesian_samples() {
    with_serving(|serving| {
        let jd = 1.0;
        let flags = FLG_SWIEPH | FLG_SPEED | FLG_TRUEPOS | FLG_NOABERR | FLG_NOGDEFL;
        let geo = [0.5, 0.1, 0.0];
        let frame = |tt| {
            frame_transform(
                geo,
                tt,
                flags,
                true,
                frame_matrices_for_date(tt, 0.0, 0.0).unwrap(),
            )
            .unwrap()
        };
        let central = frame(jd);
        let minus = frame(jd - VEL_H_DAYS);
        let plus = frame(jd + VEL_H_DAYS);
        let legacy_rate = difference_stencil_rate(
            (minus[0], minus[1], minus[2]),
            (plus[0], plus[1], plus[2]),
            true,
        )
        .unwrap();
        let direct = serving_compute_xyz_tt(&serving, jd, MARS, flags).unwrap();
        assert_eq!(direct.0[..3], central);
        // Differentiate the known constant ICRS ray's rotating output frame
        // independently with degree-six centered interpolants. No ephemeris
        // speed or native SPEED result supplies these expectations.
        for h in [1.0 / 32.0, 1.0 / 64.0, 1.0 / 128.0] {
            let expected: [f64; 3] = std::array::from_fn(|component| {
                let delta = |n: f64| frame(jd + n * h)[component] - frame(jd - n * h)[component];
                (45.0 * delta(1.0) - 9.0 * delta(2.0) + delta(3.0)) / (60.0 * h)
            });
            for (component, reference) in expected.into_iter().enumerate() {
                assert!(
                    (direct.0[component + 3] - reference).abs() <= reference.abs() * 1e-6,
                    "Cartesian component {component}: {} vs independent {reference}",
                    direct.0[component + 3]
                );
            }
        }
        eprintln!(
            "legacy two-point Cartesian rate {legacy_rate:?}; complete-place rate {:?}",
            &direct.0[3..]
        );
    });
}

#[test]
fn serving_channel_gate_names_every_refusal() {
    let file = fixture_icrs_const(&[(MARS, [1.5, 0.1, 0.0])]);
    let owned = FileIndex::decode(&file).expect("fixture must decode");
    let composite = CompositeIndex::new(vec![(owned, &file[..])]).expect("one reader must build");
    // Fictitious runtime-model id: never LEB-served, even if stored.
    assert!(matches!(
        serving_channel_gate(&composite, CUPIDO).unwrap_err(),
        Error::UnknownBody { .. }
    ));
    // Absent body: the typed miss signal.
    assert!(matches!(
        serving_channel_gate(&composite, VENUS).unwrap_err(),
        Error::LebBodyNotFound { .. }
    ));
    // ICRS passes with its channel.
    assert_eq!(
        serving_channel_gate(&composite, MARS).unwrap(),
        StateFrame::IcrsBarycentric
    );
}

#[test]
fn serving_channel_gate_refuses_ecliptic_direct() {
    let file = fixture(MEAN_NODE, COORD_ECLIPTIC);
    let owned = FileIndex::decode(&file).expect("fixture must decode");
    let composite = CompositeIndex::new(vec![(owned, &file[..])]).expect("one reader must build");
    let err = serving_channel_gate(&composite, MEAN_NODE).unwrap_err();
    assert!(
        matches!(
            err,
            Error::SealedChannelNotServed {
                body_id,
                coord_type: 1
            } if body_id == MEAN_NODE
        ),
        "unexpected refusal: {err:?}"
    );
}

#[test]
fn serving_direct_cartesian_preserves_vector_words_and_refusals() {
    with_serving(|serving| {
        let flags = FLG_TRUEPOS | FLG_EQUATORIAL | FLG_J2000 | FLG_ICRS;
        // With no corrections or rotation, direct Cartesian output
        // must preserve target-minus-Earth exactly, not round-trip
        // through longitude/latitude/distance.
        let expected = [0.5, 0.1, 0.0];
        for speed in [0, FLG_SPEED] {
            let vector = serving_cartesian_tt(&serving, 1.0, MARS, flags | speed).unwrap();
            assert_eq!(vector.map(f64::to_bits), expected.map(f64::to_bits));
        }
        for (jd, body, extra) in [
            (f64::NAN, MARS, 0),
            (1.0, VENUS, 0),
            (10.0, MARS, 0),
            (1.0, MARS, FLG_TOPOCTR),
        ] {
            let direct = serving_cartesian_tt(&serving, jd, body, flags | extra).unwrap_err();
            let spherical = serving_compute_tt(&serving, jd, body, flags | extra).unwrap_err();
            assert_eq!(format!("{direct:?}"), format!("{spherical:?}"));
        }
    });
}

#[test]
fn serving_truepos_geocentric_matches_frame_routing() {
    with_serving(|serving| {
        let jd = 1.0;
        let flags = FLG_SWIEPH;
        let (out, retflag) = serving_compute_tt(&serving, jd, MARS, flags).expect("serve");
        assert_eq!(retflag, flags);
        assert!((0.0..360.0).contains(&out[0]), "lon wrapped: {out:?}");
        assert!(out[2] > 0.0, "positive distance: {out:?}");
        for v in out {
            assert!(v.is_finite(), "finite core: {out:?}");
        }
    });
}

#[test]
fn serving_geocentric_earth_returns_origin_zeros() {
    with_serving(|serving| {
        let (start, end) = serving.provider().body_coverage(EARTH).expect("coverage");
        for jd in [start, 1.0, end] {
            for flags in [
                FLG_SWIEPH,
                FLG_SWIEPH | FLG_SPEED,
                FLG_SWIEPH | FLG_EQUATORIAL | FLG_XYZ,
                FLG_SWIEPH | FLG_SPEED | FLG_EQUATORIAL | FLG_XYZ,
            ] {
                let (out, retflag) = serving_compute_tt(&serving, jd, EARTH, flags).expect("serve");
                assert_eq!(retflag, flags, "retflag echo for {flags}");
                assert_eq!(out.map(f64::to_bits), [0; 6], "Earth origin at {jd}");
            }
            let flags = FLG_SWIEPH | FLG_HELCTR | FLG_SPEED;
            let (out, retflag) = serving_compute_tt(&serving, jd, SUN, flags).expect("Sun origin");
            assert_eq!(retflag, flags);
            assert_eq!(out.map(f64::to_bits), [0; 6], "Sun origin at {jd}");
        }
        // Short-circuiting derived reads never waives central coverage.
        for jd in [start - 1.0, end + 1.0] {
            assert!(serving_compute_tt(&serving, jd, EARTH, FLG_SPEED).is_err());
        }
    });
}

#[test]
fn leb_vector_edge_preserves_reference_position_flags_and_radial_rate() {
    use crate::context::EphemerisContext;
    use crate::leb::LebFile;
    use crate::time::DeltaTModel;
    use sha2::{Digest, Sha256};
    // Reference words and source data must remain bound to the same image.
    let path = crate::test_data::published_core_path();
    let bytes = std::fs::read(&path)
        .expect("the published fixture bound to these reference words must exist");
    let hash: String = Sha256::digest(&bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect();
    assert_eq!(hash, crate::download::BASE_CORE_SHA256);
    let mut ctx = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
    ctx.attach_leb_file(LebFile::from_bytes(path.to_string_lossy().into_owned(), bytes).unwrap());
    // Independent pinned position/range-rate values at the first covered
    // polynomial. Angular rates are checked against the actual Cartesian
    // coefficient derivative in stored_vector_retry_tests; reproducing the
    // reference implementation's short coordinate-difference stencil is not
    // a physical invariant.
    let jd = f64::from_bits(0x4142492b400d1b71);
    let (actual, flags) =
        with_sealed_calc_engine(&ctx, DeltaTModel::EspenakMeeus, false, |engine| {
            calc(engine, jd, MOON, FLG_SWIEPH | FLG_TRUEPOS | FLG_SPEED)
        })
        .unwrap();
    assert_eq!(flags, 1810);
    for (index, word) in [
        (0, 0x4060cd3be870c448),
        (1, 0xbfed0336e16fcf0c),
        (2, 0x3f63d069a5d5c930),
        (5, 0x3efe3491454a5000),
    ] {
        assert!((actual[index] - f64::from_bits(word)).abs() <= 1e-8);
    }
    assert!(actual.into_iter().all(f64::is_finite));
}

#[test]
fn serving_truepos_pinned_against_manual_composition() {
    with_serving(|serving| {
        let jd = 1.0;
        let flags = FLG_SWIEPH | FLG_TRUEPOS;
        let (out, retflag) = serving_compute_tt(&serving, jd, MARS, flags).expect("serve");
        assert_eq!(retflag, flags);

        let rel = compose_relative_state(
            [1.5, 0.1, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            ObservationCenter::Geocentric,
            None,
        )
        .expect("compose");
        let frames = frame_matrices_for_date(jd, 0.0, 0.0).expect("frames");
        let expected = frame_transform(rel.geo_pos(), jd, flags, false, frames).expect("frame");
        assert_eq!(out[0], expected[0], "lon bit-identical under TRUEPOS");
        assert_eq!(out[1], expected[1], "lat bit-identical under TRUEPOS");
        assert_eq!(out[2], expected[2], "dist bit-identical under TRUEPOS");
        // No SPEED requested: exact zero rates.
        assert_eq!(&out[3..6], &[0.0, 0.0, 0.0]);
    });
}

#[test]
fn serving_degenerate_self_observations_are_exact_zeros() {
    with_serving(|serving| {
        let helio_sun =
            serving_compute_tt(&serving, 1.0, SUN, FLG_SWIEPH | FLG_HELCTR | FLG_TRUEPOS)
                .expect("sun helio");
        assert_eq!(helio_sun.0, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let geo_earth =
            serving_compute_tt(&serving, 1.0, EARTH, FLG_SWIEPH | FLG_TRUEPOS).expect("earth geo");
        assert_eq!(geo_earth.0, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    });
}

#[test]
fn serving_speed_path_is_finite_and_gated() {
    with_serving(|serving| {
        let (fast, _) =
            serving_compute_tt(&serving, 1.0, MARS, FLG_SWIEPH | FLG_SPEED).expect("speed");
        for v in fast {
            assert!(v.is_finite(), "finite speed core: {fast:?}");
        }

        let (slow, _) = serving_compute_tt(&serving, 1.0, MARS, FLG_SWIEPH).expect("no speed");
        assert_eq!(&slow[3..6], &[0.0, 0.0, 0.0]);
        assert_eq!(&slow[0..3], &fast[0..3], "position unaffected by SPEED");
    });
}

#[test]
fn serving_refusals_propagate_typed() {
    with_serving(|serving| {
        // Non-finite epoch.
        assert!(matches!(
            serving_compute_tt(&serving, f64::NAN, MARS, FLG_SWIEPH).unwrap_err(),
            Error::NonFiniteArgument { .. }
        ));
        // Outside the [0, 2] fixture coverage.
        assert!(matches!(
            serving_compute_tt(&serving, 5.0, MARS, FLG_SWIEPH).unwrap_err(),
            Error::LebOutOfRange { .. }
        ));

        assert!(matches!(
            serving_compute_tt(&serving, 1.0, MARS, FLG_SWIEPH | FLG_TOPOCTR).unwrap_err(),
            Error::ConfigurationMissing { .. }
        ));
        // Unknown body.
        assert!(matches!(
            serving_compute_tt(&serving, 1.0, CUPIDO, FLG_SWIEPH).unwrap_err(),
            Error::UnknownBody { .. }
        ));
    });
}

#[test]
fn escalated_sealed_core_miss_names_requested_body_with_its_window() {
    with_serving(|serving| {
        // Sealed core out-of-range: the provider miss becomes the
        // public EphemerisRange with the requested body's own window.
        let err =
            serving_compute_tt_escalated(ProviderMode::SealedLeb, &serving, 5.0, MARS, FLG_SWIEPH)
                .unwrap_err();
        assert_eq!(
            err,
            Error::EphemerisRange {
                requested_jd: 5.0,
                start_jd: 0.0,
                end_jd: 2.0,
                body: Some(MARS),
            }
        );
        // A served epoch is untouched by the escalation wrapper.
        assert!(
            serving_compute_tt_escalated(ProviderMode::SealedLeb, &serving, 1.0, MARS, FLG_SWIEPH)
                .is_ok()
        );

        let err = serving_compute_tt_escalated(ProviderMode::Auto, &serving, 5.0, MARS, FLG_SWIEPH)
            .unwrap_err();
        assert_eq!(
            err,
            Error::LebOutOfRange {
                body_id: MARS,
                jd: 5.0,
                start: 0.0,
                end: 2.0,
            }
        );
        // Sealed missing stored-state core body: declared missing
        // body, not a date-range error.
        let err = escalate_sealed_provider_miss(
            ProviderMode::SealedLeb,
            serving.provider,
            VENUS,
            1.0,
            Error::LebBodyNotFound { body_id: VENUS },
        );
        assert_eq!(err, Error::UnknownBody { body_id: VENUS });
        // Fictitious dispatch ids never escalate.
        let miss = Error::LebOutOfRange {
            body_id: CUPIDO,
            jd: 5.0,
            start: 0.0,
            end: 2.0,
        };
        let err = escalate_sealed_provider_miss(
            ProviderMode::SealedLeb,
            serving.provider,
            CUPIDO,
            5.0,
            miss,
        );
        assert_eq!(
            err,
            Error::LebOutOfRange {
                body_id: CUPIDO,
                jd: 5.0,
                start: 0.0,
                end: 2.0,
            }
        );
        // Non-finite epochs never escalate either.
        let miss = Error::LebOutOfRange {
            body_id: MARS,
            jd: f64::NAN,
            start: 0.0,
            end: 2.0,
        };
        let err = escalate_sealed_provider_miss(
            ProviderMode::SealedLeb,
            serving.provider,
            MARS,
            f64::NAN,
            miss,
        );
        assert!(matches!(err, Error::LebOutOfRange { .. }));
    });
}

#[test]
fn escalated_minor_miss_names_requested_body_only_past_support() {
    // Sun-only file over [0, 2]: the non-core CHIRON has no channel,
    // so its missed reads exercise the support rule.
    let file = fixture_icrs_const(&[(SUN, [0.0, 0.0, 0.0])]);
    let owned = FileIndex::decode(&file).expect("fixture must decode");
    let composite = CompositeIndex::new(vec![(owned, &file[..])]).expect("one reader must build");
    // Absent minor body while the Sun support covers the epoch: the
    // fall-through proceeds, so the miss propagates unchanged.
    let err = escalate_sealed_provider_miss(
        ProviderMode::SealedLeb,
        &composite,
        CHIRON,
        1.0,
        Error::LebBodyNotFound { body_id: CHIRON },
    );
    assert_eq!(err, Error::LebBodyNotFound { body_id: CHIRON });
    // Absent minor body past the Sun window: name the requested body
    // with the support window (never a "Body 0" mislabel).
    let err = escalate_sealed_provider_miss(
        ProviderMode::SealedLeb,
        &composite,
        CHIRON,
        5.0,
        Error::LebBodyNotFound { body_id: CHIRON },
    );
    assert_eq!(
        err,
        Error::EphemerisRange {
            requested_jd: 5.0,
            start_jd: 0.0,
            end_jd: 2.0,
            body: Some(CHIRON),
        }
    );
    // Own-window minor miss with support in range: fall through
    // unchanged.
    let err = escalate_sealed_provider_miss(
        ProviderMode::SealedLeb,
        &composite,
        CHIRON,
        1.5,
        Error::LebOutOfRange {
            body_id: CHIRON,
            jd: 1.5,
            start: 0.0,
            end: 1.0,
        },
    );
    assert_eq!(
        err,
        Error::LebOutOfRange {
            body_id: CHIRON,
            jd: 1.5,
            start: 0.0,
            end: 1.0,
        }
    );
    // Own-window minor miss past support too: escalate with the
    // requested window.
    let err = escalate_sealed_provider_miss(
        ProviderMode::SealedLeb,
        &composite,
        CHIRON,
        5.0,
        Error::LebOutOfRange {
            body_id: CHIRON,
            jd: 5.0,
            start: 0.0,
            end: 1.0,
        },
    );
    assert_eq!(
        err,
        Error::EphemerisRange {
            requested_jd: 5.0,
            start_jd: 0.0,
            end_jd: 1.0,
            body: Some(CHIRON),
        }
    );
    // A support-read miss naming another body is never relabeled.
    let err = escalate_sealed_provider_miss(
        ProviderMode::SealedLeb,
        &composite,
        CHIRON,
        5.0,
        Error::LebOutOfRange {
            body_id: SUN,
            jd: 5.0,
            start: 0.0,
            end: 2.0,
        },
    );
    assert_eq!(
        err,
        Error::LebOutOfRange {
            body_id: SUN,
            jd: 5.0,
            start: 0.0,
            end: 2.0,
        }
    );
}

#[test]
fn escalated_retired_minor_frame_miss_retains_requested_coverage() {
    let mut file = fixture_icrs_const(&[(SUN, [0.0, 0.0, 0.0]), (CHIRON, [2.0, 0.1, 0.0])]);
    let coord_at = HEADER_SIZE + SECTION_DIR_SIZE + LEB1_BODY_ENTRY_SIZE + 4;
    file[coord_at..coord_at + 4].copy_from_slice(&2_u32.to_le_bytes());
    let owned = FileIndex::decode(&file).expect("retired frame is valid metadata");
    let composite = CompositeIndex::new(vec![(owned, &file[..])]).unwrap();
    let serving = SealedServing::new(&composite, &zero_nutation);
    let refusal = Error::SealedChannelNotServed {
        body_id: CHIRON,
        coord_type: 2,
    };
    for flags in [
        FLG_SWIEPH,
        FLG_SWIEPH | FLG_ICRS,
        FLG_SWIEPH | FLG_HELCTR,
        FLG_SWIEPH | FLG_XYZ,
    ] {
        assert_eq!(
            serving_compute_tt_escalated(ProviderMode::SealedLeb, &serving, 5.0, CHIRON, flags)
                .unwrap_err(),
            Error::EphemerisRange {
                requested_jd: 5.0,
                start_jd: 0.0,
                end_jd: 2.0,
                body: Some(CHIRON),
            }
        );
        assert_eq!(
            serving_compute_tt_escalated(ProviderMode::SealedLeb, &serving, 1.0, CHIRON, flags)
                .unwrap_err(),
            refusal
        );
        assert_eq!(
            serving_compute_tt_escalated(ProviderMode::Auto, &serving, 5.0, CHIRON, flags)
                .unwrap_err(),
            refusal
        );
    }
}

const IDENTITY_C2T: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

fn identity_c2t(_jd_tt: f64, _jd_ut1: f64) -> Result<[[f64; 3]; 3]> {
    Ok(IDENTITY_C2T)
}

#[test]
fn topocentric_matches_manual_stage_composition_bitwise() {
    // Aberration/deflection gated off: the serving value must equal
    // the manual observer-selection + composition + frame chain
    // bit-exactly, both without and with light-time.
    with_serving(|serving| {
        let geopos = (13.4, 52.5, 34.0);
        let jd_tt = 1.0;
        let jd_ut1 = jd_tt - 0.0008;
        for truepos in [true, false] {
            let flags = FLG_SWIEPH
                | FLG_NOABERR
                | FLG_NOGDEFL
                | FLG_TOPOCTR
                | if truepos { FLG_TRUEPOS } else { 0 };
            let (out, retflag) = serving_compute_tt_topocentric(
                &serving,
                jd_tt,
                jd_ut1,
                MARS,
                flags,
                geopos,
                &identity_c2t,
                false,
                SiteEllipsoid::Wgs84,
            )
            .expect("topo serve");
            assert_eq!(retflag, flags);
            let topo = topocentric_offset(
                geopos.0,
                geopos.1,
                geopos.2,
                IDENTITY_C2T,
                SiteEllipsoid::Wgs84,
            )
            .expect("offset");
            let target = serving.provider.geometric_state(MARS, jd_tt).expect("mars");
            let sun = serving.provider.geometric_state(SUN, jd_tt).expect("sun");
            let earth = serving
                .provider
                .geometric_state(EARTH, jd_tt)
                .expect("earth");
            let rel = compose_relative_state(
                target.pos(),
                target.vel(),
                sun.pos(),
                sun.vel(),
                earth.pos(),
                earth.vel(),
                ObservationCenter::Topocentric,
                Some(topo.pair()),
            )
            .expect("compose");
            let geo = if truepos {
                rel.geo_pos()
            } else {
                let earth_state = serving
                    .provider
                    .geometric_state(EARTH, jd_tt)
                    .expect("earth");
                let topo2 = topocentric_offset(
                    geopos.0,
                    geopos.1,
                    geopos.2,
                    IDENTITY_C2T,
                    SiteEllipsoid::Wgs84,
                )
                .expect("offset");
                let (tp, _) = topo2.pair();
                let ep = earth_state.pos();
                let obs = [ep[0] + tp[0], ep[1] + tp[1], ep[2] + tp[2]];
                solve_light_time(
                    |jd| {
                        let state = serving
                            .provider
                            .geometric_state(MARS, jd)
                            .expect("retarded");
                        Ok((state.pos(), state.vel()))
                    },
                    jd_tt,
                    obs,
                    rel.aberration_vel(),
                    true,
                )
                .expect("lt")
                .geo_pos()
            };
            let frames = frame_matrices_for_date(jd_tt, 0.0, 0.0).expect("frames");
            let shaped = frame_transform(geo, jd_tt, flags, false, frames).expect("frame");
            assert_eq!([out[0], out[1], out[2]], [shaped[0], shaped[1], shaped[2]]);
            assert_eq!([out[3], out[4], out[5]], [0.0, 0.0, 0.0]);
        }
    });
}

#[test]
fn topocentric_samples_use_the_frozen_scale_difference() {
    use std::cell::RefCell;
    with_serving(|serving| {
        let geopos = (13.4, 52.5, 34.0);
        let jd_tt = 1.0;
        let jd_ut1 = jd_tt - 0.0008;
        let seen: RefCell<Vec<(f64, f64)>> = RefCell::new(Vec::new());
        let recording = |tt: f64, ut1: f64| -> Result<[[f64; 3]; 3]> {
            seen.borrow_mut().push((tt, ut1));
            Ok(IDENTITY_C2T)
        };
        let flags = FLG_SWIEPH | FLG_TOPOCTR | FLG_SPEED;
        let (out, _) = serving_compute_tt_topocentric(
            &serving,
            jd_tt,
            jd_ut1,
            MARS,
            flags,
            geopos,
            &recording,
            false,
            SiteEllipsoid::Wgs84,
        )
        .expect("topo speed");
        for v in out {
            assert!(v.is_finite(), "non-finite topo component");
        }
        // Every adaptive sample uses the same physical TT-minus-UT1
        // offset; its spacing and number do not define that clock relation.
        let seen = seen.borrow();
        assert!(seen.contains(&(jd_tt, jd_ut1)));
        assert!(seen.iter().any(|(tt, _)| *tt < jd_tt));
        assert!(seen.iter().any(|(tt, _)| *tt > jd_tt));
        let frozen_offset = jd_tt - jd_ut1;
        for &(tt, ut1) in seen.iter() {
            let epoch_scale = tt.abs() + ut1.abs() + jd_tt.abs() + jd_ut1.abs();
            let roundoff = f64::EPSILON * epoch_scale;
            assert!(((tt - ut1) - frozen_offset).abs() <= roundoff);
        }
        // Speed slots carry the diurnal rate; the position slots match
        // the no-speed computation exactly.
        let (slow, _) = serving_compute_tt_topocentric(
            &serving,
            jd_tt,
            jd_ut1,
            MARS,
            FLG_SWIEPH | FLG_TOPOCTR,
            geopos,
            &identity_c2t,
            false,
            SiteEllipsoid::Wgs84,
        )
        .expect("topo no speed");
        assert_eq!([out[0], out[1], out[2]], [slow[0], slow[1], slow[2]]);
        assert!([out[3], out[4], out[5]] != [0.0, 0.0, 0.0]);
        assert_eq!([slow[3], slow[4], slow[5]], [0.0, 0.0, 0.0]);
    });
}

#[test]
fn topocentric_delegates_without_flag_and_pins_earth_zero() {
    with_serving(|serving| {
        let geopos = (13.4, 52.5, 34.0);
        // No TOPOCTR bit: the site triple is ignored, geocentric value.
        let (topo, _) = serving_compute_tt_topocentric(
            &serving,
            1.0,
            1.0 - 0.0008,
            MARS,
            FLG_SWIEPH,
            geopos,
            &identity_c2t,
            false,
            SiteEllipsoid::Wgs84,
        )
        .expect("no-flag topo call");
        let (geo, _) = serving_compute_tt(&serving, 1.0, MARS, FLG_SWIEPH).expect("geocentric");
        assert_eq!(topo, geo);

        let (earth, _) = serving_compute_tt_topocentric(
            &serving,
            1.0,
            1.0 - 0.0008,
            EARTH,
            FLG_SWIEPH | FLG_TOPOCTR | FLG_TRUEPOS,
            geopos,
            &identity_c2t,
            false,
            SiteEllipsoid::Wgs84,
        )
        .expect("topo earth");
        let topo = topocentric_offset(
            geopos.0,
            geopos.1,
            geopos.2,
            IDENTITY_C2T,
            SiteEllipsoid::Wgs84,
        )
        .expect("offset");
        let (tp, _) = topo.pair();
        let frames = frame_matrices_for_date(1.0, 0.0, 0.0).expect("frames");
        let expect = frame_transform(
            [-tp[0], -tp[1], -tp[2]],
            1.0,
            FLG_SWIEPH | FLG_TOPOCTR | FLG_TRUEPOS,
            false,
            frames,
        )
        .expect("frame");
        assert!(
            [earth[0], earth[1], earth[2]] != [0.0, 0.0, 0.0],
            "fast path serves the offset, never the entry zero"
        );
        // Last-ulp association: the pipeline forms
        // `earth - (earth + topo)` (cancellation at the observer
        // scale) where the pin negates the offset directly.
        for k in 0..3 {
            assert!(
                (earth[k] - expect[k]).abs() < 1e-9,
                "slot {k}: {} vs {}",
                earth[k],
                expect[k]
            );
        }
        assert!(earth[2] > 0.0);
        // Non-finite UT1 partner is a typed refusal.
        assert!(matches!(
            serving_compute_tt_topocentric(
                &serving,
                1.0,
                f64::NAN,
                MARS,
                FLG_SWIEPH | FLG_TOPOCTR,
                geopos,
                &identity_c2t,
                false,
                SiteEllipsoid::Wgs84,
            )
            .unwrap_err(),
            Error::NonFiniteArgument { .. }
        ));
    });
}

#[test]
fn ecl_nut_tuple_layout_matches_release_contract() {
    // Hand-written nutation pair (radians): the tuple layout, the
    // true = mean + deps assembly, and the flag pass-through.
    let pair = |_jd: f64| -> Result<(f64, f64)> { Ok((0.001, 0.0005)) };
    let (out, ret) = ecl_nut_result(2451545.0, FLG_SWIEPH, &pair).expect("ecl_nut");
    let mean = vondrak_mean_obliquity_deg(2451545.0).expect("mean");
    assert_eq!(out[1], mean);
    assert_eq!(out[2], 0.001f64.to_degrees());
    assert_eq!(out[3], 0.0005f64.to_degrees());
    assert_eq!(out[0], mean + 0.0005f64.to_degrees());
    assert_eq!([out[4], out[5]], [0.0, 0.0]);
    assert_eq!(ret, FLG_SWIEPH);

    assert!((mean - OBLIQUITY_J2000_DEG).abs() < 5e-5, "mean = {mean}");
    // Refusals: non-finite epoch, pair, and provider error.
    assert!(matches!(
        ecl_nut_result(f64::NAN, 0, &pair).unwrap_err(),
        Error::NonFiniteArgument { .. }
    ));
    let bad = |_jd: f64| -> Result<(f64, f64)> { Ok((f64::INFINITY, 0.0)) };
    assert!(matches!(
        ecl_nut_result(2451545.0, 0, &bad).unwrap_err(),
        Error::NonFiniteArgument { .. }
    ));
    let failing = |_jd: f64| -> Result<(f64, f64)> { Err(Error::LebBodyNotFound { body_id: 99 }) };
    assert_eq!(
        ecl_nut_result(2451545.0, 0, &failing).unwrap_err(),
        Error::LebBodyNotFound { body_id: 99 }
    );
}

#[test]
fn observer_helper_matches_compose_selection() {
    let target = [1.5, 0.1, 0.0];
    let sun = [0.1, 0.0, 0.0];
    let earth = [1.0, 0.0, 0.0];
    let topo = Some(([0.001, 0.0, 0.0], [0.0, 0.0001, 0.0]));
    for center in [
        ObservationCenter::Heliocentric,
        ObservationCenter::Barycentric,
        ObservationCenter::Geocentric,
        ObservationCenter::Topocentric,
    ] {
        let (op, ov) = observer_state_for_center(center, sun, [0.0; 3], earth, [0.0; 3], topo);
        let rel = compose_relative_state(
            target, [0.0; 3], sun, [0.0; 3], earth, [0.0; 3], center, topo,
        )
        .expect("compose");
        assert_eq!(
            [target[0] - op[0], target[1] - op[1], target[2] - op[2]],
            rel.geo_pos(),
            "observer agrees on {center:?}"
        );
        let _ = ov;
    }
}
