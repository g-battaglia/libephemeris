use super::*;

mod inline_state;

#[path = "image_evaluation.rs"]
mod image_evaluation;
#[path = "mapping.rs"]
mod mapping;
#[path = "memory.rs"]
mod memory;
#[path = "ownership_index.rs"]
mod ownership_index;
mod positional;
mod prepared_tier_bounds;
#[path = "release_contracts.rs"]
mod release_contracts;
mod segment_domain;
#[path = "tiered_delta_t.rs"]
mod tiered_delta_t;

fn curvature_leb1_image() -> Vec<u8> {
    let offset = HEADER_SIZE + SECTION_DIR_SIZE + LEB1_BODY_ENTRY_SIZE;
    let mut file = pack_header(&MAGIC_LEB1, 1, 1, 1, 0.0, 4.0, 0.0, 0);
    file.extend(pack_section(
        SECTION_BODY_INDEX,
        (HEADER_SIZE + SECTION_DIR_SIZE) as u64,
        LEB1_BODY_ENTRY_SIZE as u64,
    ));
    file.extend(pack_body_leb1(
        0,
        COORD_ICRS_BARY,
        2,
        0.0,
        4.0,
        2.0,
        3,
        3,
        offset as u64,
    ));
    // x=T3, y=2*T2, z=constant on each two-day record. The right
    // record changes curvature so a join cannot borrow its predecessor.
    for multiplier in [1.0, 2.0] {
        for value in [
            0.0,
            0.0,
            0.0,
            multiplier,
            0.0,
            0.0,
            2.0 * multiplier,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
        ] {
            push_f64(&mut file, value);
        }
    }
    file
}

fn check_curvature_image(file: &[u8]) -> LebChunkCache {
    let index = FileIndex::decode(file).unwrap();
    let cache = LebChunkCache::default();
    let samples = [
        (TwoPartEpoch { hi: 0.0, lo: 0.0 }, [-24.0, 8.0, 0.0]),
        (TwoPartEpoch { hi: 1.0, lo: 0.25 }, [6.0, 8.0, 0.0]),
        (TwoPartEpoch { hi: 2.0, lo: 0.0 }, [-48.0, 16.0, 0.0]),
        (TwoPartEpoch { hi: 4.0, lo: 0.0 }, [48.0, 16.0, 0.0]),
    ];
    for (epoch, expected) in samples.into_iter().chain(samples) {
        assert_eq!(
            index
                .eval_raw_acceleration_two_part_source(LebSource::Slice(file), 0, epoch, &cache)
                .unwrap(),
            Some(expected)
        );
    }
    assert!(matches!(
        index.eval_raw_acceleration_two_part_source(
            LebSource::Slice(file),
            0,
            TwoPartEpoch { hi: 4.0, lo: 0.25 },
            &cache
        ),
        Err(Error::LebOutOfRange { .. })
    ));
    cache
}

#[test]
fn cartesian_acceleration_uses_the_selected_record_and_split_epoch() {
    check_curvature_image(&curvature_leb1_image());
}

#[test]
fn cartesian_acceleration_shares_compressed_records_and_chunk_ownership() {
    let leb1 = curvature_leb1_image();
    let raw = &leb1[HEADER_SIZE + SECTION_DIR_SIZE + LEB1_BODY_ENTRY_SIZE..];
    for version in [1, 2] {
        let mut file = pack_header(&MAGIC_LEB2, version, 1, 1, 0.0, 4.0, 0.0, 1);
        let body_at = file.len() + SECTION_DIR_SIZE;
        let data_at = body_at + LEB2_BODY_ENTRY_SIZE;
        file.extend(pack_section(
            SECTION_BODY_INDEX,
            body_at as u64,
            LEB2_BODY_ENTRY_SIZE as u64,
        ));
        let blob = if version == 1 {
            encode_frame(raw, 2, 3, 3)
        } else {
            let first = encode_frame(&raw[..96], 1, 3, 3);
            let second = encode_frame(&raw[96..], 1, 3, 3);
            let first_at = data_at + CHUNK_INDEX_HEADER_SIZE + 2 * CHUNK_ENTRY_SIZE;
            let mut blob = pack_chunk_header(2, 2.0);
            // Deliberately inaccurate advisory bound: segment ownership wins.
            blob.extend(pack_chunk_entry(
                0.0,
                2.5,
                first_at as u64,
                first.len() as u64,
                96,
                0,
                1,
            ));
            blob.extend(pack_chunk_entry(
                2.0,
                4.0,
                (first_at + first.len()) as u64,
                second.len() as u64,
                96,
                1,
                1,
            ));
            blob.extend(first);
            blob.extend(second);
            blob
        };
        file.extend(pack_body_leb2(
            0,
            COORD_ICRS_BARY,
            2,
            0.0,
            4.0,
            2.0,
            3,
            3,
            data_at as u64,
            blob.len() as u64,
            raw.len() as u64,
        ));
        file.extend(blob);
        let (hits, misses, entries) = check_curvature_image(&file).stats();
        assert_eq!(misses, u64::from(version));
        assert_eq!(entries, version as usize);
        assert_eq!(hits + misses, 8);
    }
}

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

/// 64-byte common header, field order `<4sIIIdddI20s`.
#[allow(clippy::too_many_arguments)]
fn pack_header(
    magic: &[u8; 4],
    version: u32,
    sections: u32,
    bodies: u32,
    jd_start: f64,
    jd_end: f64,
    generation_epoch: f64,
    flags: u32,
) -> Vec<u8> {
    let mut v = Vec::with_capacity(HEADER_SIZE);
    v.extend_from_slice(magic);
    push_u32(&mut v, version);
    push_u32(&mut v, sections);
    push_u32(&mut v, bodies);
    push_f64(&mut v, jd_start);
    push_f64(&mut v, jd_end);
    push_f64(&mut v, generation_epoch);
    push_u32(&mut v, flags);
    v.extend_from_slice(&[0u8; 20]);
    debug_assert_eq!(v.len(), HEADER_SIZE);
    v
}

/// 24-byte section directory entry (`<IIQQ`, reserved zero).
fn pack_section(id: u32, offset: u64, size: u64) -> Vec<u8> {
    let mut v = Vec::with_capacity(SECTION_DIR_SIZE);
    push_u32(&mut v, id);
    push_u32(&mut v, 0);
    push_u64(&mut v, offset);
    push_u64(&mut v, size);
    v
}

/// 52-byte LEB1 body entry (`<iIIdddIIQ`).
#[allow(clippy::too_many_arguments)]
fn pack_body_leb1(
    body_id: i32,
    coord: u32,
    segs: u32,
    jd_start: f64,
    jd_end: f64,
    interval: f64,
    degree: u32,
    comps: u32,
    data_offset: u64,
) -> Vec<u8> {
    let mut v = Vec::with_capacity(LEB1_BODY_ENTRY_SIZE);
    push_i32(&mut v, body_id);
    push_u32(&mut v, coord);
    push_u32(&mut v, segs);
    push_f64(&mut v, jd_start);
    push_f64(&mut v, jd_end);
    push_f64(&mut v, interval);
    push_u32(&mut v, degree);
    push_u32(&mut v, comps);
    push_u64(&mut v, data_offset);
    v
}

/// 68-byte LEB2 body entry (LEB1 fields plus blob sizes).
#[allow(clippy::too_many_arguments)]
fn pack_body_leb2(
    body_id: i32,
    coord: u32,
    segs: u32,
    jd_start: f64,
    jd_end: f64,
    interval: f64,
    degree: u32,
    comps: u32,
    data_offset: u64,
    compressed: u64,
    uncompressed: u64,
) -> Vec<u8> {
    let mut v = pack_body_leb1(
        body_id,
        coord,
        segs,
        jd_start,
        jd_end,
        interval,
        degree,
        comps,
        data_offset,
    );
    push_u64(&mut v, compressed);
    push_u64(&mut v, uncompressed);
    debug_assert_eq!(v.len(), LEB2_BODY_ENTRY_SIZE);
    v
}

/// F1: minimal one-body LEB1 — header, one section entry, one body
/// entry, and a 96-byte coefficient blob (`2·3·2·8`).
fn minimal_leb1() -> Vec<u8> {
    let mut file = pack_header(&MAGIC_LEB1, 1, 1, 1, 0.0, 2.0, 0.0, 0);
    let dir_at = file.len();
    let body_at = dir_at + SECTION_DIR_SIZE;
    let blob_at = body_at + LEB1_BODY_ENTRY_SIZE;
    file.extend_from_slice(&pack_section(0, body_at as u64, 52));
    file.extend_from_slice(&pack_body_leb1(
        0,
        0,
        2,
        0.0,
        2.0,
        1.0,
        1,
        3,
        blob_at as u64,
    ));
    // Hand-picked blob: constant channel T0 = 5, slope channel T1 = 1
    // per segment; LEB-02 only checks its bounds, LEB-06 evaluates it.
    for seg in 0..2 {
        for comp in 0..3 {
            push_f64(&mut file, if comp == 0 { 5.0 } else { 0.0 });
            push_f64(&mut file, if comp == 1 { 1.0 } else { 0.0 });
            let _ = seg;
        }
    }
    debug_assert_eq!(file.len(), blob_at + 96);
    file
}

#[test]
fn decode_minimal_leb1_file_index() {
    let file = minimal_leb1();
    let index = FileIndex::decode(&file).expect("valid fixture must decode");
    let header = index.header();
    assert_eq!(header.magic(), MAGIC_LEB1);
    assert_eq!(header.version(), 1);
    assert_eq!(header.section_count(), 1);
    assert_eq!(header.body_count(), 1);
    assert_eq!(header.jd_start(), 0.0);
    assert_eq!(header.jd_end(), 2.0);
    assert_eq!(header.format(), LebFormat::Leb1);
    assert_eq!(index.sections().len(), 1);
    assert_eq!(index.section(0).unwrap().size(), 52);
    assert!(index.section(6).is_none());
    assert_eq!(index.bodies().len(), 1);
    let body = &index.bodies()[0];
    assert_eq!(body.body_id(), 0);
    assert_eq!(body.coord_type(), COORD_ICRS_BARY);
    assert_eq!(body.segment_count(), 2);
    assert_eq!(body.jd_start(), 0.0);
    assert_eq!(body.jd_end(), 2.0);
    assert_eq!(body.interval_days(), 1.0);
    assert_eq!(body.degree(), 1);
    assert_eq!(body.components(), 3);
    assert_eq!(body.data_offset(), 140);
    assert_eq!(body.compressed_size(), None);
    assert_eq!(body.expected_uncompressed_size().unwrap(), 96);
}

#[test]
fn detect_classifies_formats_and_rejects_identity_failures() {
    let leb1 = pack_header(&MAGIC_LEB1, 1, 0, 0, 0.0, 1.0, 0.0, 0);
    assert_eq!(LebFormat::detect(&leb1).unwrap(), LebFormat::Leb1);
    let v1 = pack_header(&MAGIC_LEB2, 1, 0, 0, 0.0, 1.0, 0.0, 1);
    assert_eq!(LebFormat::detect(&v1).unwrap(), LebFormat::Leb2V1);
    let v2 = pack_header(&MAGIC_LEB2, 2, 0, 0, 0.0, 1.0, 0.0, 1);
    assert_eq!(LebFormat::detect(&v2).unwrap(), LebFormat::Leb2V2);
    // LEB1 ignores flags.
    let leb1_flags = pack_header(&MAGIC_LEB1, 1, 0, 0, 0.0, 1.0, 0.0, 0);
    assert_eq!(LebFormat::detect(&leb1_flags).unwrap(), LebFormat::Leb1);

    let bad = pack_header(b"BAD!", 1, 0, 0, 0.0, 1.0, 0.0, 0);
    assert_eq!(
        LebFormat::detect(&bad).unwrap_err(),
        Error::LebUnknownFormat { magic: *b"BAD!" }
    );
    let bad_version = pack_header(&MAGIC_LEB1, 2, 0, 0, 0.0, 1.0, 0.0, 0);
    assert_eq!(
        LebFormat::detect(&bad_version).unwrap_err(),
        Error::LebUnsupportedVersion { found: 2 }
    );
    let bad_v9 = pack_header(&MAGIC_LEB2, 9, 0, 0, 0.0, 1.0, 0.0, 1);
    assert_eq!(
        LebFormat::detect(&bad_v9).unwrap_err(),
        Error::LebUnsupportedVersion { found: 9 }
    );
    let no_comp = pack_header(&MAGIC_LEB2, 2, 0, 0, 0.0, 1.0, 0.0, 0);
    assert_eq!(
        LebFormat::detect(&no_comp).unwrap_err(),
        Error::LebUnsupportedCompression { found: 0 }
    );
    // 43 bytes cannot hold magic/version/flags.
    assert!(matches!(
        LebFormat::detect(&leb1[..43]).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
}

#[test]
fn malformed_headers_and_directories_are_corrupt() {
    let file = minimal_leb1();
    // Truncated header.
    assert!(matches!(
        FileIndex::decode(&file[..63]).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    // Declared directory entry missing.
    assert!(matches!(
        FileIndex::decode(&file[..HEADER_SIZE]).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    // Section payload past EOF.
    let mut bad = pack_header(&MAGIC_LEB1, 1, 1, 0, 0.0, 1.0, 0.0, 0);
    bad.extend_from_slice(&pack_section(0, 10_000, 52));
    assert!(matches!(
        FileIndex::decode(&bad).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    // Inflated body count trips the count guard (2^30 entries do not
    // fit in 52 section bytes).
    let mut inflated = pack_header(&MAGIC_LEB1, 1, 1, 1 << 30, 0.0, 1.0, 0.0, 0);
    inflated.extend_from_slice(&pack_section(0, 88, 52));
    assert!(matches!(
        FileIndex::decode(&inflated).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    // No section 0 at all.
    let mut no_index = pack_header(&MAGIC_LEB1, 1, 1, 1, 0.0, 1.0, 0.0, 0);
    no_index.extend_from_slice(&pack_section(4, 88, 64));
    no_index.extend_from_slice(&[0u8; 64]);
    assert!(matches!(
        FileIndex::decode(&no_index).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
}

fn minimal_leb1_with_body(body: Vec<u8>, body_count: u32) -> Vec<u8> {
    let mut file = pack_header(&MAGIC_LEB1, 1, 1, body_count, 0.0, 2.0, 0.0, 0);
    let dir_at = file.len();
    let body_at = dir_at + SECTION_DIR_SIZE;
    file.extend_from_slice(&pack_section(0, body_at as u64, body.len() as u64));
    file.extend_from_slice(&body);
    file.extend_from_slice(&[0u8; 96]);
    file
}

#[test]
fn malformed_body_entries_are_corrupt() {
    // Zero interval.
    let file = minimal_leb1_with_body(pack_body_leb1(0, 0, 2, 0.0, 2.0, 0.0, 1, 3, 140), 1);
    assert!(matches!(
        FileIndex::decode(&file).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    // Zero segment count.
    let file = minimal_leb1_with_body(pack_body_leb1(0, 0, 0, 0.0, 2.0, 1.0, 1, 3, 140), 1);
    assert!(matches!(
        FileIndex::decode(&file).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    // Data offset past EOF.
    let file = minimal_leb1_with_body(pack_body_leb1(0, 0, 2, 0.0, 2.0, 1.0, 1, 3, 1_000_000), 1);
    assert!(matches!(
        FileIndex::decode(&file).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    // Non-finite interval.
    let file = minimal_leb1_with_body(pack_body_leb1(0, 0, 2, 0.0, 2.0, f64::NAN, 1, 3, 140), 1);
    assert!(matches!(
        FileIndex::decode(&file).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    // Truncated body entry.
    let mut file = minimal_leb1();
    file.truncate(88 + 51);
    assert!(matches!(
        FileIndex::decode(&file).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    // Retired coord type 2 still decodes (calculation path decides).
    let file = minimal_leb1_with_body(
        pack_body_leb1(0, COORD_HELIO_ECL_RETIRED, 2, 0.0, 2.0, 1.0, 1, 3, 140),
        1,
    );
    let index = FileIndex::decode(&file).expect("retired type must decode");
    assert_eq!(index.bodies()[0].coord_type(), 2);
}

#[test]
fn leb2_body_entry_carries_blob_sizes() {
    let entry_bytes = pack_body_leb2(
        5,
        COORD_ICRS_BARY_SYSTEM,
        3425,
        0.0,
        100.0,
        31.99211678,
        13,
        3,
        500,
        1234,
        3425 * 14 * 3 * 8,
    );
    let entry = BodyEntry::decode(&entry_bytes, 0, 10_000, LebFormat::Leb2V2)
        .expect("valid LEB2 entry must decode");
    assert_eq!(entry.body_id(), 5);
    assert_eq!(entry.compressed_size(), Some(1234));
    assert_eq!(entry.uncompressed_size(), Some(3425 * 14 * 3 * 8));
    assert_eq!(
        entry.expected_uncompressed_size().unwrap(),
        3425 * 14 * 3 * 8
    );
}

#[test]
fn aux_headers_f4_decode_exactly() {
    // Nutation header: 1 segment, degree 1, 2 components, hand values.
    let mut nut = Vec::with_capacity(NUTATION_HEADER_SIZE);
    push_f64(&mut nut, 2400000.5);
    push_f64(&mut nut, 2400016.5);
    push_f64(&mut nut, 16.0);
    push_u32(&mut nut, 1);
    push_u32(&mut nut, 2);
    push_u32(&mut nut, 1);
    push_u32(&mut nut, 0);
    let header = NutationHeader::decode(&nut, 0).expect("valid nutation header");
    assert_eq!(header.jd_start(), 2400000.5);
    assert_eq!(header.jd_end(), 2400016.5);
    assert_eq!(header.interval_days(), 16.0);
    assert_eq!(header.degree(), 1);
    assert_eq!(header.components(), 2);
    assert_eq!(header.segment_count(), 1);
    assert!(header.is_present());
    // `--skip-aux` shape: segment_count == 0 decodes as absent.
    let mut absent = nut.clone();
    absent[32..36].copy_from_slice(&0u32.to_le_bytes());
    let header = NutationHeader::decode(&absent, 0).expect("absent table decodes");
    assert!(!header.is_present());
    assert!(matches!(
        NutationHeader::decode(&nut[..39], 0).unwrap_err(),
        Error::LebCorrupt { .. }
    ));

    // Delta-T header: 2 entries; capacity guard passes on 40 section
    // bytes and trips on an inflated count.
    let mut dt = Vec::with_capacity(DELTA_T_HEADER_SIZE);
    push_u32(&mut dt, 2);
    push_u32(&mut dt, 0);
    let header = DeltaTHeader::decode(&dt, 0).expect("valid Delta-T header");
    assert_eq!(header.n_entries(), 2);
    header
        .validate_capacity(8 + 2 * 16, 1000)
        .expect("exact fit passes");
    assert!(matches!(
        header.validate_capacity(8 + 2 * 16 - 1, 1000).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    let mut dt_bad = Vec::new();
    push_u32(&mut dt_bad, 1 << 30);
    push_u32(&mut dt_bad, 0);
    let bad = DeltaTHeader::decode(&dt_bad, 0).expect("count decodes");
    assert!(matches!(
        bad.validate_capacity(40, 1000).unwrap_err(),
        Error::LebCorrupt { .. }
    ));

    // One star entry, exact field round-trip.
    let mut star = Vec::with_capacity(STAR_ENTRY_SIZE);
    push_i32(&mut star, 12345);
    for x in [10.5, -20.25, 0.001, -0.002, 0.05, 11.0, 2.5] {
        push_f64(&mut star, x);
    }
    star.extend_from_slice(b"RSVD");
    let entry = StarEntry::decode(&star, 0).expect("valid star entry");
    assert_eq!(entry.star_id(), 12345);
    assert_eq!(entry.ra_j2000_deg(), 10.5);
    assert_eq!(entry.dec_j2000_deg(), -20.25);
    assert_eq!(entry.pm_ra_deg_per_yr(), 0.001);
    assert_eq!(entry.pm_dec_deg_per_yr(), -0.002);
    assert_eq!(entry.parallax_arcsec(), 0.05);
    assert_eq!(entry.rv_km_per_s(), 11.0);
    assert_eq!(entry.magnitude(), 2.5);
    assert!(matches!(
        StarEntry::decode(&star[..63], 0).unwrap_err(),
        Error::LebCorrupt { .. }
    ));

    // Star count is floor division: short tails are ignored.
    assert_eq!(star_count(92608), 1447);
    assert_eq!(star_count(92608 + 63), 1447);
    assert_eq!(star_count(63), 0);
}

#[test]
fn entry_count_guard_covers_exact_fit_shortfall_and_overflow() {
    validate_entry_count(2, 52, 104, "LEB body index", 88).expect("exact fit");
    assert!(matches!(
        validate_entry_count(2, 52, 103, "LEB body index", 88).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    assert!(matches!(
        validate_entry_count(u64::MAX, 2, u64::MAX, "LEB body index", 0).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    validate_entry_count(0, 52, 0, "LEB body index", 0).expect("empty index fits");
}

#[test]
fn monolithic_eval_matches_hand_computation_bitwise() {
    // F1 blob: segment-major, per segment comp0 = [5, 0] (constant 5),
    // comp1 = [0, 1] (T1 = tau), comp2 = [0, 0]. Degree-1 Clenshaw is
    // exact here, so every expectation below is bitwise, not tolerant:
    // value(comp0) = 5, value(comp1) = tau, d/dtau(comp1) = 1,
    // velocity scale 2/1.0 = 2.
    let file = minimal_leb1();
    let index = FileIndex::decode(&file).expect("valid fixture must decode");
    assert!(index.has_body(0));
    assert!(!index.has_body(99));
    assert_eq!(index.body_coverage(0), Some((0.0, 2.0)));
    assert_eq!(index.body_coverage(99), None);

    // jd = 0.25: seg 0, tau = 2*(0.25 - 0.5) = -0.5.
    let (pos, vel) = index.eval_body(&file, 0, 0.25).unwrap();
    assert_eq!(pos, vec![5.0, -0.5, 0.0]);
    assert_eq!(vel, vec![0.0, 2.0, 0.0]);
    // jd = 0.0: seg 0, tau clamped to -1.
    let (pos, vel) = index.eval_body(&file, 0, 0.0).unwrap();
    assert_eq!(pos, vec![5.0, -1.0, 0.0]);
    assert_eq!(vel, vec![0.0, 2.0, 0.0]);
    // jd = 1.0, the interior seam: int truncation owns it to seg 1,
    // tau = 2*(1.0 - 1.5) = -1.
    let (pos, vel) = index.eval_body(&file, 0, 1.0).unwrap();
    assert_eq!(pos, vec![5.0, -1.0, 0.0]);
    assert_eq!(vel, vec![0.0, 2.0, 0.0]);
    // jd = 2.0, the file end: clamped onto the last segment, tau = 1.
    let (pos, vel) = index.eval_body(&file, 0, 2.0).unwrap();
    assert_eq!(pos, vec![5.0, 1.0, 0.0]);
    assert_eq!(vel, vec![0.0, 2.0, 0.0]);
}

#[test]
fn eval_rejects_range_unknown_and_nonfinite() {
    let file = minimal_leb1();
    let index = FileIndex::decode(&file).expect("valid fixture must decode");
    for jd in [-0.5, -1e-9, 2.0 + 1e-9, 100.0] {
        assert_eq!(
            index.eval_body(&file, 0, jd).unwrap_err(),
            Error::LebOutOfRange {
                body_id: 0,
                jd,
                start: 0.0,
                end: 2.0,
            }
        );
    }
    assert_eq!(
        index.eval_body(&file, 99, 1.0).unwrap_err(),
        Error::LebBodyNotFound { body_id: 99 }
    );
    assert!(matches!(
        index.eval_body(&file, 0, f64::NAN).unwrap_err(),
        Error::NonFiniteArgument { .. }
    ));
    assert!(matches!(
        index.eval_body(&file, 0, f64::INFINITY).unwrap_err(),
        Error::NonFiniteArgument { .. }
    ));
}

/// Two ecliptic bodies plus one ICRS body, one segment each, degree 0:
/// ecliptic longitudes use positive-divisor remainder (negative included),
/// ICRS channels pass through, and degree-0 velocities are exactly 0.
fn wrap_fixture() -> Vec<u8> {
    let mut file = pack_header(&MAGIC_LEB1, 1, 1, 3, 0.0, 8.0, 0.0, 0);
    let dir_at = file.len();
    let entries_at = dir_at + SECTION_DIR_SIZE;
    let blob_at = entries_at + 3 * LEB1_BODY_ENTRY_SIZE;
    file.extend_from_slice(&pack_section(0, entries_at as u64, 3 * 52));
    file.extend_from_slice(&pack_body_leb1(
        10,
        1,
        1,
        0.0,
        8.0,
        8.0,
        0,
        3,
        blob_at as u64,
    ));
    file.extend_from_slice(&pack_body_leb1(
        11,
        1,
        1,
        0.0,
        8.0,
        8.0,
        0,
        3,
        (blob_at + 24) as u64,
    ));
    file.extend_from_slice(&pack_body_leb1(
        12,
        0,
        1,
        0.0,
        8.0,
        8.0,
        0,
        3,
        (blob_at + 48) as u64,
    ));
    for v in [370.0, -5.0, 1.5] {
        push_f64(&mut file, v);
    }
    for v in [-10.0, 45.0, 2.0] {
        push_f64(&mut file, v);
    }
    for v in [370.0, -5.0, 1.5] {
        push_f64(&mut file, v);
    }
    debug_assert_eq!(file.len(), blob_at + 72);
    file
}

#[test]
fn ecliptic_longitude_wraps_icrs_passes_through() {
    let file = wrap_fixture();
    let index = FileIndex::decode(&file).expect("valid fixture must decode");
    // tau = 0 at mid-segment; degree-0 values are the constants.
    let (pos, vel) = index.eval_body(&file, 10, 4.0).unwrap();
    assert_eq!(pos, vec![10.0, -5.0, 1.5]);
    assert_eq!(vel, vec![0.0, 0.0, 0.0]);
    // A negative longitude wraps to the corresponding positive angle.
    let (pos, _) = index.eval_body(&file, 11, 4.0).unwrap();
    assert_eq!(pos, vec![350.0, 45.0, 2.0]);
    // ICRS coord type: no wrap even above 360.
    let (pos, _) = index.eval_body(&file, 12, 4.0).unwrap();
    assert_eq!(pos, vec![370.0, -5.0, 1.5]);
}

#[test]
fn truncated_blob_is_corrupt_per_segment() {
    let mut file = minimal_leb1();
    // One byte short of the full 96-byte blob: seg 0 still fits,
    // seg 1 does not. Decode succeeds (only data_offset < EOF is a
    // decode-time check); evaluation enforces per-segment bounds.
    file.truncate(140 + 95);
    let index = FileIndex::decode(&file).expect("truncated blob still decodes");
    let (pos, _) = index.eval_body(&file, 0, 0.25).unwrap();
    assert_eq!(pos, vec![5.0, -0.5, 0.0]);
    assert!(matches!(
        index.eval_body(&file, 0, 1.5).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
}

#[test]
fn leb2_payload_refused_on_raw_path() {
    // Minimal LEB2 v2 image with a body entry; no blob bytes are even
    // needed — the raw path refuses before reading them.
    let mut file = pack_header(&MAGIC_LEB2, 2, 1, 1, 0.0, 2.0, 0.0, 1);
    file.extend_from_slice(&pack_section(0, 88, 68));
    file.extend_from_slice(&pack_body_leb2(0, 0, 2, 0.0, 2.0, 1.0, 1, 3, 156, 10, 96));
    // Padding so `data_offset` is inside the image: decode succeeds
    // and the raw-path gate is what fires.
    file.extend_from_slice(&[0u8; 96]);
    let index = FileIndex::decode(&file).expect("LEB2 index decodes");
    let err = index.eval_body(&file, 0, 1.0).unwrap_err();
    assert!(
        matches!(err, Error::LebCorrupt { .. }),
        "LEB2 on the raw path must be LebCorrupt, got {err:?}"
    );
    assert!(
        format!("{err}").contains("LEB-06"),
        "refusal must name the owning task, got {err}"
    );
}

fn push_f32(v: &mut Vec<u8>, x: f32) {
    v.extend_from_slice(&x.to_le_bytes());
}

/// Encode segment-major raw bytes through the generator-side pipeline
/// (coefficient-major reorder, byte shuffle, zstd), mirroring release
/// `leb_compression.compress_body` with full mantissa precision:
/// truncation is generator-side lossy policy and out of scope for the
/// reader contract, so fixtures keep all 52 bits.
fn encode_frame(segment_major: &[u8], segs: u32, deg: u32, comps: u32) -> Vec<u8> {
    let coeff_major = reorder_to_coeff_major(segment_major, segs, deg, comps)
        .expect("fixture geometry must be consistent");
    let shuffled = shuffle_bytes(&coeff_major).expect("fixture must shuffle");
    zstd::encode_all(std::io::Cursor::new(&shuffled), 19).expect("fixture must compress")
}

/// Raw segment-major fixture: 2 segments × 3 components × degree 1
/// (`2·3·2·8 = 96` bytes), the same hand values as the F1 LEB1 blob:
/// per segment, comp0 = `[5, 0]` (constant 5), comp1 = `[0, 1]`
/// (`T1`), comp2 = `[0, 0]`.
fn raw_two_segment_blob() -> Vec<u8> {
    let mut v = Vec::with_capacity(96);
    for _ in 0..2 {
        for comp in 0..3 {
            push_f64(&mut v, if comp == 0 { 5.0 } else { 0.0 });
            push_f64(&mut v, if comp == 1 { 1.0 } else { 0.0 });
        }
    }
    v
}

/// 16-byte chunk index header (`<IIfI`).
fn pack_chunk_header(chunk_count: u32, interval_days: f32) -> Vec<u8> {
    let mut v = Vec::with_capacity(CHUNK_INDEX_HEADER_SIZE);
    push_u32(&mut v, chunk_count);
    push_u32(&mut v, 0);
    push_f32(&mut v, interval_days);
    push_u32(&mut v, 0);
    debug_assert_eq!(v.len(), CHUNK_INDEX_HEADER_SIZE);
    v
}

/// 48-byte chunk index entry (`<ddQQQII`).
#[allow(clippy::too_many_arguments)]
fn pack_chunk_entry(
    jd_start: f64,
    jd_end: f64,
    blob_offset: u64,
    compressed: u64,
    uncompressed: u64,
    seg_start: u32,
    seg_count: u32,
) -> Vec<u8> {
    let mut v = Vec::with_capacity(CHUNK_ENTRY_SIZE);
    push_f64(&mut v, jd_start);
    push_f64(&mut v, jd_end);
    push_u64(&mut v, blob_offset);
    push_u64(&mut v, compressed);
    push_u64(&mut v, uncompressed);
    push_u32(&mut v, seg_start);
    push_u32(&mut v, seg_count);
    debug_assert_eq!(v.len(), CHUNK_ENTRY_SIZE);
    v
}

#[test]
fn shuffle_transposes_lanes_and_round_trips() {
    // Three words with distinct bytes: word j holds bytes j*8..j*8+8.
    let raw: Vec<u8> = (0u8..24).collect();
    let shuffled = shuffle_bytes(&raw).expect("multiple of 8 shuffles");
    // Lane `l` holds byte `l` of every word: raw[j*8+l].
    for lane in 0..8 {
        for j in 0..3 {
            assert_eq!(shuffled[lane * 3 + j], raw[j * 8 + lane]);
        }
    }
    assert_eq!(unshuffle_bytes(&shuffled).unwrap(), raw);
    // The other direction round-trips too.
    let back = shuffle_bytes(&unshuffle_bytes(&raw).unwrap()).unwrap();
    assert_eq!(back, raw);
    // Empty maps to empty; a partial word is corruption, not padding.
    assert_eq!(shuffle_bytes(&[]).unwrap(), Vec::<u8>::new());
    assert_eq!(unshuffle_bytes(&[]).unwrap(), Vec::<u8>::new());
    assert!(matches!(
        shuffle_bytes(&[0u8; 7]).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    assert!(matches!(
        unshuffle_bytes(&[0u8; 15]).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
}

#[test]
fn reorder_is_exact_bit_permutation_and_round_trips() {
    // segs = 2, comps = 3, deg = 1: 12 words, seg-major [s*6+k].
    let words: Vec<f64> = (0..12).map(|i| i as f64).collect();
    let bytes: Vec<u8> = words.iter().flat_map(|w| w.to_le_bytes()).collect();
    let cm = reorder_to_coeff_major(&bytes, 2, 1, 3).expect("fixture must reorder");
    let cm_words: Vec<f64> = cm
        .chunks_exact(8)
        .map(|c| f64::from_le_bytes(c.try_into().unwrap()))
        .collect();
    // Coeff-major words[k*2+s] = seg-major[s*6+k].
    assert_eq!(
        cm_words,
        vec![0.0, 6.0, 1.0, 7.0, 2.0, 8.0, 3.0, 9.0, 4.0, 10.0, 5.0, 11.0]
    );
    assert_eq!(reorder_to_segment_major(&cm, 2, 1, 3).unwrap(), bytes);
    // Non-finite words (truncated payloads stay binary64) permute
    // bit-exactly, never normalized or rejected here.
    let mut nan_bytes = bytes.clone();
    nan_bytes[0..8].copy_from_slice(&f64::NAN.to_le_bytes());
    nan_bytes[8..16].copy_from_slice(&f64::INFINITY.to_le_bytes());
    let round = reorder_to_segment_major(
        &reorder_to_coeff_major(&nan_bytes, 2, 1, 3).unwrap(),
        2,
        1,
        3,
    )
    .unwrap();
    assert_eq!(round, nan_bytes);
    // Length contract: one byte short or one word over is corruption.
    assert!(matches!(
        reorder_to_coeff_major(&bytes[..95], 2, 1, 3).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    let mut over = bytes.clone();
    over.extend_from_slice(&[0u8; 8]);
    assert!(matches!(
        reorder_to_segment_major(&over, 2, 1, 3).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
}

/// F2: a complete LEB2 v1 file image whose blob is a test-time zstd
/// frame decodes through the index entry back to the exact raw
/// segment-major bytes.
fn leb2_v1_image() -> (Vec<u8>, Vec<u8>) {
    let raw = raw_two_segment_blob();
    let blob = encode_frame(&raw, 2, 1, 3);
    let mut file = pack_header(&MAGIC_LEB2, 1, 2, 1, 0.0, 2.0, 0.0, 1);
    let body_at = file.len() + 2 * SECTION_DIR_SIZE;
    let blob_at = body_at + LEB2_BODY_ENTRY_SIZE;
    file.extend_from_slice(&pack_section(
        0,
        body_at as u64,
        LEB2_BODY_ENTRY_SIZE as u64,
    ));
    file.extend_from_slice(&pack_section(6, blob_at as u64, blob.len() as u64));
    file.extend_from_slice(&pack_body_leb2(
        0,
        0,
        2,
        0.0,
        2.0,
        1.0,
        1,
        3,
        blob_at as u64,
        blob.len() as u64,
        raw.len() as u64,
    ));
    file.extend_from_slice(&blob);
    debug_assert_eq!(file.len(), blob_at + blob.len());
    (file, raw)
}

#[test]
fn leb2_v1_monolith_decodes_to_raw_segment_major() {
    let (file, raw) = leb2_v1_image();
    let index = FileIndex::decode(&file).expect("LEB2 v1 image must decode");
    assert_eq!(index.header().format(), LebFormat::Leb2V1);
    let entry = index.find_body(0).expect("body 0 must be indexed");
    assert_eq!(entry.expected_uncompressed_size().unwrap(), 96);
    let compressed = read_blob(
        &file,
        entry.data_offset(),
        entry.compressed_size().unwrap(),
        "LEB2 test blob",
    )
    .expect("blob must slice");
    assert_eq!(decompress_monolith(compressed, entry).unwrap(), raw);
    // Stored-size mismatch twin: 95 stored vs 96 derived.
    let bad_bytes = pack_body_leb2(0, 0, 2, 0.0, 2.0, 1.0, 1, 3, 500, 10, 95);
    let bad = BodyEntry::decode(&bad_bytes, 0, 10_000, LebFormat::Leb2V1)
        .expect("sizes still decode structurally");
    assert!(matches!(
        decompress_monolith(compressed, &bad).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    // An LEB1 entry (no stored sizes) on the blob path is corruption.
    let leb1_bytes = pack_body_leb1(0, 0, 2, 0.0, 2.0, 1.0, 1, 3, 140);
    let leb1 =
        BodyEntry::decode(&leb1_bytes, 0, 10_000, LebFormat::Leb1).expect("LEB1 entry decodes");
    assert!(matches!(
        decompress_monolith(compressed, &leb1).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
}

/// F3: two chunks over the two F2 segments, each an independent
/// test-time frame; returns the index buffer plus the raw bytes.
fn two_chunk_buffer(chunk0_end: f64) -> (Vec<u8>, Vec<u8>) {
    let raw = raw_two_segment_blob();
    let blob0 = encode_frame(&raw[..48], 1, 1, 3);
    let blob1 = encode_frame(&raw[48..], 1, 1, 3);
    let mut buf = pack_chunk_header(2, 3652.5);
    let blob0_at = (CHUNK_INDEX_HEADER_SIZE + 2 * CHUNK_ENTRY_SIZE) as u64;
    let blob1_at = blob0_at + blob0.len() as u64;
    buf.extend_from_slice(&pack_chunk_entry(
        0.0,
        chunk0_end,
        blob0_at,
        blob0.len() as u64,
        48,
        0,
        1,
    ));
    buf.extend_from_slice(&pack_chunk_entry(
        1.0,
        2.0,
        blob1_at,
        blob1.len() as u64,
        48,
        1,
        1,
    ));
    buf.extend_from_slice(&blob0);
    buf.extend_from_slice(&blob1);
    (buf, raw)
}

#[test]
fn leb2_v2_chunks_decode_with_segment_ownership() {
    let (buf, raw) = two_chunk_buffer(1.0);
    let (header, chunks) = decode_chunk_index(&buf, 0, 0).expect("chunk index decodes");
    assert_eq!(header.chunk_count(), 2);
    assert_eq!(header.chunk_interval_days(), 3652.5f32);
    assert_eq!(chunks.len(), 2);
    assert_eq!(chunks[0].segment_start(), 0);
    assert_eq!(chunks[1].segment_start(), 1);
    // Binary-search hints.
    assert_eq!(find_chunk(&chunks, 0.25), Some(0));
    assert_eq!(find_chunk(&chunks, 1.5), Some(1));
    // Ownership confirms the hints.
    assert_eq!(select_chunk(&chunks, 0, 0, 0).unwrap(), 0);
    assert_eq!(select_chunk(&chunks, 1, 1, 0).unwrap(), 1);
    // Each chunk decodes to its own segment bytes exactly.
    for (i, chunk) in chunks.iter().enumerate() {
        let compressed = read_blob(&buf, chunk.blob_offset(), chunk.compressed_size(), "chunk")
            .expect("chunk blob must slice");
        let decoded = decompress_chunk(compressed, chunk, 1, 3, 0, i).unwrap();
        assert_eq!(decoded, raw[i * 48..(i + 1) * 48]);
    }
}

#[test]
fn chunk_lookup_trusts_segment_index_over_rounded_bounds() {
    // Boundary twin: chunk 0's JD bound is rounded half a ulp below
    // the true seam at jd 1.0, as write-time rounding can produce.
    // The release contract trusts the segment index, never the
    // rounded chunk boundary.
    let seam = 1.0f64;
    let rounded = f64::from_bits(seam.to_bits() - 1);
    assert!(rounded < seam);
    let (buf, _) = two_chunk_buffer(rounded);
    let (_, chunks) = decode_chunk_index(&buf, 0, 0).expect("index decodes");
    // jd exactly on the seam: truncation owns it to seg 1, and the
    // rounded bound routes the search to chunk 1, which covers it.
    assert_eq!(find_chunk(&chunks, seam), Some(1));
    assert_eq!(select_chunk(&chunks, 1, 1, 0).unwrap(), 1);
    // jd just below the seam: truncation owns it to seg 0, but the
    // rounded bound misroutes the search to chunk 1 — the ±1 probe
    // recovers chunk 0 instead of evaluating the wrong segment.
    assert_eq!(find_chunk(&chunks, rounded), Some(1));
    assert_eq!(select_chunk(&chunks, 0, 1, 0).unwrap(), 0);
}

#[test]
fn malformed_chunk_payloads_are_corrupt() {
    let (buf, _) = two_chunk_buffer(1.0);
    // Truncated index: header claims 2 chunks, one entry present.
    assert!(matches!(
        decode_chunk_index(&buf[..CHUNK_INDEX_HEADER_SIZE + CHUNK_ENTRY_SIZE], 0, 0).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    // Truncated header itself.
    assert!(matches!(
        decode_chunk_index(&buf[..10], 0, 0).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    // Inflated chunk count trips the count guard.
    let mut inflated = pack_chunk_header(u32::MAX, 3652.5);
    inflated.extend_from_slice(&buf[CHUNK_INDEX_HEADER_SIZE..]);
    assert!(matches!(
        decode_chunk_index(&inflated, 0, 0).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    // Index base past EOF.
    assert!(matches!(
        decode_chunk_index(&buf, buf.len() as u64 + 64, 0).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    let (_, chunks) = decode_chunk_index(&buf, 0, 0).expect("index decodes");
    // Blob offset past EOF.
    assert!(matches!(
        read_blob(&buf, u64::MAX, 10, "chunk").unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    assert!(matches!(
        read_blob(&buf, chunks[0].blob_offset(), u64::MAX, "chunk").unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    // Short frame: last byte of chunk 0's blob missing.
    let short: &[u8] = &buf[chunks[0].blob_offset() as usize
        ..chunks[0].blob_offset() as usize + chunks[0].compressed_size() as usize - 1];
    assert!(matches!(
        decompress_chunk(short, &chunks[0], 1, 3, 0, 0).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    // Stored-size mismatch twin: 47 stored vs 48 derived.
    let mut size_bad = pack_chunk_header(1, 3652.5);
    size_bad.extend_from_slice(&pack_chunk_entry(0.0, 1.0, 500, 10, 47, 0, 1));
    let (_, bad_chunks) = decode_chunk_index(&size_bad, 0, 7).expect("sizes decode");
    assert!(matches!(
        decompress_chunk(&[0u8; 10], &bad_chunks[0], 1, 3, 7, 0).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    // Non-frame bytes are not valid zstd.
    assert!(matches!(
        decompress_body(b"not a zstd frame", 96, 2, 1, 3).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    // Declared-size mismatch at the body level (95 vs 96 derived).
    let overlong_frame = encode_frame(&raw_two_segment_blob(), 2, 1, 3);
    assert!(matches!(
        decompress_body(&overlong_frame, 95, 2, 1, 3).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    // Segment covered by no chunk near the hint is corruption, never
    // a silent neighbouring evaluation: chunk 1 starts at seg 5.
    let mut gap = pack_chunk_header(2, 3652.5);
    gap.extend_from_slice(&pack_chunk_entry(0.0, 1.0, 500, 10, 48, 0, 1));
    gap.extend_from_slice(&pack_chunk_entry(1.0, 2.0, 600, 10, 48, 5, 1));
    let (_, gap_chunks) = decode_chunk_index(&gap, 0, 0).expect("gap index decodes");
    assert!(matches!(
        select_chunk(&gap_chunks, 1, 1, 0).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    // Out-of-range hint with no covering neighbour.
    assert!(matches!(
        select_chunk(&gap_chunks, 0, 7, 0).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    // Empty index: no hint, no selection.
    let empty: Vec<ChunkEntry> = Vec::new();
    assert_eq!(find_chunk(&empty, 0.5), None);
    assert!(matches!(
        select_chunk(&empty, 0, 0, 0).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
}

#[test]
fn scalar_selector_endpoints_and_seams_are_exact() {
    // Integer grid: jd_start 0, 2-day interval, 3 segments.
    let (seg, tau) = select_segment(0.0, 2.0, 3, 0.0).expect("start selects");
    assert_eq!(seg, 0);
    assert_eq!(tau.value(), -1.0);
    let (seg, tau) = select_segment(0.0, 2.0, 3, 6.0).expect("end selects");
    assert_eq!(seg, 2);
    assert_eq!(tau.value(), 1.0);
    // Interior seams belong to the right-hand segment at tau -1.
    for (seam, owner) in [(2.0, 1), (4.0, 2)] {
        let (seg, tau) = select_segment(0.0, 2.0, 3, seam).expect("seam selects");
        assert_eq!(seg, owner);
        assert_eq!(tau.value(), -1.0);
    }
    // Segment midpoints map to tau 0.
    for (mid, owner) in [(1.0, 0), (3.0, 1), (5.0, 2)] {
        let (seg, tau) = select_segment(0.0, 2.0, 3, mid).expect("mid selects");
        assert_eq!(seg, owner);
        assert_eq!(tau.value(), 0.0);
    }
}

#[test]
fn scalar_selector_covers_adjacent_float_words() {
    // Words around the interior seam 2.0: the word below stays in
    // segment 0 just under tau 1, the word above opens segment 1.
    let below = f64::next_down(2.0);
    let (seg, tau) = select_segment(0.0, 2.0, 3, below).expect("below selects");
    assert_eq!(seg, 0);
    assert!(tau.value() < 1.0 && tau.value() > 0.999_999_999_999);
    let above = f64::next_up(2.0);
    let (seg, tau) = select_segment(0.0, 2.0, 3, above).expect("above selects");
    assert_eq!(seg, 1);
    assert!(tau.value() > -1.0 && tau.value() < -0.999_999_999_999);
    // Outer edges: words just outside the table clamp to the edge
    // segment (callers range-check first; the clamp mirrors the
    // release `max(0, min(...))`).
    let (seg, _) = select_segment(0.0, 2.0, 3, f64::next_down(0.0)).expect("clamps low");
    assert_eq!(seg, 0);
    let (seg, tau) = select_segment(0.0, 2.0, 3, f64::next_up(6.0)).expect("clamps high");
    assert_eq!(seg, 2);
    assert_eq!(tau.value(), 1.0);
}

#[test]
fn scalar_selector_rejects_bad_inputs() {
    assert!(matches!(
        select_segment(0.0, 2.0, 3, f64::NAN).unwrap_err(),
        Error::NonFiniteArgument { .. }
    ));
    assert!(matches!(
        select_segment(0.0, 2.0, 3, f64::INFINITY).unwrap_err(),
        Error::NonFiniteArgument { .. }
    ));
    for interval in [0.0, -2.0, f64::NAN, f64::INFINITY] {
        assert!(
            matches!(
                select_segment(0.0, interval, 3, 1.0).unwrap_err(),
                Error::LebCorrupt { .. }
            ),
            "interval {interval} refused"
        );
    }
    assert!(matches!(
        select_segment(0.0, 2.0, 0, 1.0).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
}

/// Segment-major blob with distinct per-segment values (the body-evaluation layer):
/// seg0 comp0 = `[5, 0]` (constant 5), comp1 = `[0, 1]` (`T1`),
/// comp2 = `[0, 0]`; seg1 comp0 = `[7, 0]`, comp1 = `[1, 2]`
/// (`1 + 2*tau`), comp2 = `[0, 0]`. Degree-1 Clenshaw is exact, so
/// every expectation below is bitwise.
fn raw_two_segment_distinct_blob() -> Vec<u8> {
    let mut v = Vec::with_capacity(96);
    for seg in 0..2 {
        for comp in 0..3 {
            let (c0, c1) = match (seg, comp) {
                (0, 0) => (5.0, 0.0),
                (0, 1) => (0.0, 1.0),
                (1, 0) => (7.0, 0.0),
                (1, 1) => (1.0, 2.0),
                _ => (0.0, 0.0),
            };
            push_f64(&mut v, c0);
            push_f64(&mut v, c1);
        }
    }
    v
}

/// Complete LEB2 v1 image over the distinct two-segment blob: header,
/// sections 0 + 6, one body entry whose `data_offset` points at the
/// test-time zstd frame.
fn leb2_v1_distinct_image() -> Vec<u8> {
    let raw = raw_two_segment_distinct_blob();
    let blob = encode_frame(&raw, 2, 1, 3);
    let mut file = pack_header(&MAGIC_LEB2, 1, 2, 1, 0.0, 2.0, 0.0, 1);
    let body_at = file.len() + 2 * SECTION_DIR_SIZE;
    let blob_at = body_at + LEB2_BODY_ENTRY_SIZE;
    file.extend_from_slice(&pack_section(
        0,
        body_at as u64,
        LEB2_BODY_ENTRY_SIZE as u64,
    ));
    file.extend_from_slice(&pack_section(6, blob_at as u64, blob.len() as u64));
    file.extend_from_slice(&pack_body_leb2(
        0,
        0,
        2,
        0.0,
        2.0,
        1.0,
        1,
        3,
        blob_at as u64,
        blob.len() as u64,
        raw.len() as u64,
    ));
    file.extend_from_slice(&blob);
    file
}

/// Complete LEB2 v1 image over the plain two-segment blob (the same
/// hand values as the F1 LEB1 blob): the matched-format twin of
/// [`minimal_leb1`] for the verify/compare bridge tests.
fn leb2_v1_plain_image() -> Vec<u8> {
    let raw = raw_two_segment_blob();
    let blob = encode_frame(&raw, 2, 1, 3);
    let mut file = pack_header(&MAGIC_LEB2, 1, 2, 1, 0.0, 2.0, 0.0, 1);
    let body_at = file.len() + 2 * SECTION_DIR_SIZE;
    let blob_at = body_at + LEB2_BODY_ENTRY_SIZE;
    file.extend_from_slice(&pack_section(
        0,
        body_at as u64,
        LEB2_BODY_ENTRY_SIZE as u64,
    ));
    file.extend_from_slice(&pack_section(6, blob_at as u64, blob.len() as u64));
    file.extend_from_slice(&pack_body_leb2(
        0,
        0,
        2,
        0.0,
        2.0,
        1.0,
        1,
        3,
        blob_at as u64,
        blob.len() as u64,
        raw.len() as u64,
    ));
    file.extend_from_slice(&blob);
    file
}

/// Minimal LEB1 twin with a two-component body 0 over `[0, 2]`:
/// exercises the comparison component-count skip.
fn two_component_leb1() -> Vec<u8> {
    let mut file = pack_header(&MAGIC_LEB1, 1, 1, 1, 0.0, 2.0, 0.0, 0);
    let dir_at = file.len();
    let body_at = dir_at + SECTION_DIR_SIZE;
    let blob_at = body_at + LEB1_BODY_ENTRY_SIZE;
    file.extend_from_slice(&pack_section(0, body_at as u64, 52));
    file.extend_from_slice(&pack_body_leb1(
        0,
        0,
        2,
        0.0,
        2.0,
        1.0,
        1,
        2,
        blob_at as u64,
    ));
    for _ in 0..2 {
        for comp in 0..2 {
            push_f64(&mut file, if comp == 0 { 1.0 } else { 0.0 });
            push_f64(&mut file, 0.0);
        }
    }
    debug_assert_eq!(file.len(), blob_at + 64);
    file
}

/// Minimal LEB1 twin with body 0 over the disjoint interval
/// `[10, 12]`: exercises the comparison no-shared-coverage skip.
fn shifted_leb1() -> Vec<u8> {
    let mut file = pack_header(&MAGIC_LEB1, 1, 1, 1, 10.0, 12.0, 0.0, 0);
    let dir_at = file.len();
    let body_at = dir_at + SECTION_DIR_SIZE;
    let blob_at = body_at + LEB1_BODY_ENTRY_SIZE;
    file.extend_from_slice(&pack_section(0, body_at as u64, 52));
    file.extend_from_slice(&pack_body_leb1(
        0,
        0,
        2,
        10.0,
        12.0,
        1.0,
        1,
        3,
        blob_at as u64,
    ));
    for _ in 0..2 {
        for comp in 0..3 {
            push_f64(&mut file, if comp == 0 { 5.0 } else { 0.0 });
            push_f64(&mut file, if comp == 1 { 1.0 } else { 0.0 });
        }
    }
    debug_assert_eq!(file.len(), blob_at + 96);
    file
}

#[test]
fn verify_container_passes_on_minimal_leb1() {
    let report = verify_container(&minimal_leb1()).expect("fixture must verify");
    assert_eq!(report.format(), LebFormat::Leb1);
    assert_eq!(report.body_count(), 1);
    assert!(report.passed());
    assert!(report.failures().is_empty());
    let row = &report.bodies()[0];
    assert_eq!(row.body_id(), 0);
    assert_eq!((row.jd_start(), row.jd_end()), (0.0, 2.0));
    assert_eq!(row.probes(), 3);
    assert!(row.passed());
    assert_eq!(row.error(), None);
}

#[test]
fn verify_container_passes_on_matched_leb2_v1() {
    let report = verify_container(&leb2_v1_plain_image()).expect("fixture must verify");
    assert_eq!(report.format(), LebFormat::Leb2V1);
    assert_eq!(report.body_count(), 1);
    assert!(report.passed());
    assert_eq!(report.bodies()[0].probes(), 3);
}

#[test]
fn verify_container_refuses_undecodable_bytes() {
    assert!(verify_container(&[]).is_err());
    assert!(verify_container(&[0u8; 64]).is_err());
}

#[test]
fn compare_matched_leb1_leb2_pair_agrees_exactly() {
    // Same coefficients through both flavor paths and the same
    // Clenshaw kernel: full-mantissa fixtures agree bitwise, so the
    // worst observed error is exactly zero (not a tolerance).
    let report = compare_leb2_against_leb1(&leb2_v1_plain_image(), &minimal_leb1(), 9)
        .expect("matched pair must compare");
    assert_eq!(report.compared_count(), 1);
    assert!(report.skipped().is_empty());
    let row = &report.compared()[0];
    assert_eq!(row.body_id(), 0);
    assert_eq!(row.samples(), 9);
    assert_eq!(row.max_abs_err(), 0.0);
    assert_eq!(report.overall_max_abs_err(), Some(0.0));
}

#[test]
fn compare_clamps_zero_samples_to_a_midpoint_probe() {
    let report = compare_leb2_against_leb1(&leb2_v1_plain_image(), &minimal_leb1(), 0)
        .expect("zero samples clamp to one probe");
    assert_eq!(report.compared()[0].samples(), 1);
    assert_eq!(report.compared()[0].max_abs_err(), 0.0);
}

#[test]
fn compare_records_absent_body_skips() {
    // Bodies 10/11/12 live only in the candidate, body 0 only in the
    // reference: nothing compares, everything is named.
    let report = compare_leb2_against_leb1(&wrap_fixture(), &minimal_leb1(), 4)
        .expect("disjoint bodies must skip");
    assert_eq!(report.compared_count(), 0);
    assert_eq!(report.overall_max_abs_err(), None);
    let skipped: Vec<(i32, &str)> = report
        .skipped()
        .iter()
        .map(|row| (row.body_id(), row.reason()))
        .collect();
    assert_eq!(skipped.len(), 4);
    for (id, reason) in &skipped {
        assert!(reason.contains("absent"), "body {id}: {reason}");
    }
    assert!(skipped.iter().any(|(id, _)| *id == 0));
    assert!(skipped.iter().any(|(id, _)| *id == 12));
}

#[test]
fn compare_records_component_mismatch_skip() {
    let report = compare_leb2_against_leb1(&minimal_leb1(), &two_component_leb1(), 4)
        .expect("mismatched components must skip");
    assert_eq!(report.compared_count(), 0);
    assert_eq!(report.skipped().len(), 1);
    assert_eq!(report.skipped()[0].body_id(), 0);
    assert!(report.skipped()[0].reason().contains("component count"));
}

#[test]
fn compare_records_no_shared_coverage_skip() {
    let report = compare_leb2_against_leb1(&minimal_leb1(), &shifted_leb1(), 4)
        .expect("disjoint coverage must skip");
    assert_eq!(report.compared_count(), 0);
    assert_eq!(report.skipped().len(), 1);
    assert_eq!(report.skipped()[0].reason(), "no shared coverage");
}
/// one body entry whose `data_offset` points at a two-chunk index,
/// each chunk an independent test-time frame holding its own
/// segment. `chunk0_end` sets chunk 0's advisory JD bound (rounded
/// twins pin the segment-wins rule at eval level).
fn leb2_v2_distinct_image(chunk0_end: f64) -> Vec<u8> {
    leb2_v2_scaled_image(chunk0_end, 1.0)
}

/// Synthetic scaled coefficient twin for cross-image cache isolation.
fn leb2_v2_scaled_image(chunk0_end: f64, scale: f64) -> Vec<u8> {
    let raw: Vec<u8> = raw_two_segment_distinct_blob()
        .chunks_exact(8)
        .flat_map(|word| (f64::from_le_bytes(word.try_into().unwrap()) * scale).to_le_bytes())
        .collect();
    let blob0 = encode_frame(&raw[..48], 1, 1, 3);
    let blob1 = encode_frame(&raw[48..], 1, 1, 3);
    let index_len = CHUNK_INDEX_HEADER_SIZE + 2 * CHUNK_ENTRY_SIZE;
    let mut file = pack_header(&MAGIC_LEB2, 2, 1, 1, 0.0, 2.0, 0.0, 1);
    let body_at = file.len() + SECTION_DIR_SIZE;
    let index_at = body_at + LEB2_BODY_ENTRY_SIZE;
    let blob0_at = index_at + index_len;
    let blob1_at = blob0_at + blob0.len();
    file.extend_from_slice(&pack_section(
        0,
        body_at as u64,
        LEB2_BODY_ENTRY_SIZE as u64,
    ));
    file.extend_from_slice(&pack_body_leb2(
        0,
        0,
        2,
        0.0,
        2.0,
        1.0,
        1,
        3,
        index_at as u64,
        0,
        96,
    ));
    file.extend_from_slice(&pack_chunk_header(2, 3652.5));
    file.extend_from_slice(&pack_chunk_entry(
        0.0,
        chunk0_end,
        blob0_at as u64,
        blob0.len() as u64,
        48,
        0,
        1,
    ));
    file.extend_from_slice(&pack_chunk_entry(
        1.0,
        2.0,
        blob1_at as u64,
        blob1.len() as u64,
        48,
        1,
        1,
    ));
    file.extend_from_slice(&blob0);
    file.extend_from_slice(&blob1);
    debug_assert_eq!(file.len(), blob1_at + blob1.len());
    file
}

#[test]
fn leb2_v1_eval_matches_hand_values_bitwise() {
    // jd 0.25: seg 0, tau -0.5 -> [5, -0.5, 0], vel [0, 2, 0].
    // jd 1.0 (seam): seg 1, tau -1 -> [7, -1, 0], vel [0, 4, 0].
    // jd 1.5: seg 1, tau 0 -> [7, 1, 0], vel [0, 4, 0].
    // jd 2.0 (end): seg 1, tau 1 -> [7, 3, 0], vel [0, 4, 0].
    let file = leb2_v1_distinct_image();
    let index = FileIndex::decode(&file).expect("LEB2 v1 image must decode");
    assert_eq!(index.header().format(), LebFormat::Leb2V1);
    let (pos, vel) = index.eval_leb2_body(&file, 0, 0.25).unwrap();
    assert_eq!(pos, vec![5.0, -0.5, 0.0]);
    assert_eq!(vel, vec![0.0, 2.0, 0.0]);
    let (pos, vel) = index.eval_leb2_body(&file, 0, 1.0).unwrap();
    assert_eq!(pos, vec![7.0, -1.0, 0.0]);
    assert_eq!(vel, vec![0.0, 4.0, 0.0]);
    let (pos, vel) = index.eval_leb2_body(&file, 0, 1.5).unwrap();
    assert_eq!(pos, vec![7.0, 1.0, 0.0]);
    assert_eq!(vel, vec![0.0, 4.0, 0.0]);
    let (pos, vel) = index.eval_leb2_body(&file, 0, 2.0).unwrap();
    assert_eq!(pos, vec![7.0, 3.0, 0.0]);
    assert_eq!(vel, vec![0.0, 4.0, 0.0]);
    // Ecliptic wrap twin: same image with coord type 1 and a
    // negative longitude constant in seg 0 comp 0.
    let mut ecliptic = raw_two_segment_distinct_blob();
    ecliptic[0..8].copy_from_slice(&(-10.0f64).to_le_bytes());
    let blob = encode_frame(&ecliptic, 2, 1, 3);
    let mut file = pack_header(&MAGIC_LEB2, 1, 2, 1, 0.0, 2.0, 0.0, 1);
    let body_at = file.len() + 2 * SECTION_DIR_SIZE;
    let blob_at = body_at + LEB2_BODY_ENTRY_SIZE;
    file.extend_from_slice(&pack_section(
        0,
        body_at as u64,
        LEB2_BODY_ENTRY_SIZE as u64,
    ));
    file.extend_from_slice(&pack_section(6, blob_at as u64, blob.len() as u64));
    file.extend_from_slice(&pack_body_leb2(
        3,
        1,
        2,
        0.0,
        2.0,
        1.0,
        1,
        3,
        blob_at as u64,
        blob.len() as u64,
        96,
    ));
    file.extend_from_slice(&blob);
    let index = FileIndex::decode(&file).expect("ecliptic image must decode");
    let (pos, _) = index.eval_leb2_body(&file, 3, 0.25).unwrap();
    assert_eq!(pos[0], 350.0);
}

#[test]
fn leb2_v2_chunked_eval_matches_hand_values_bitwise() {
    let file = leb2_v2_distinct_image(1.0);
    let index = FileIndex::decode(&file).expect("LEB2 v2 image must decode");
    assert_eq!(index.header().format(), LebFormat::Leb2V2);
    // Chunk 0 serves seg 0, chunk 1 serves seg 1 — same hand values
    // as the v1 test above.
    let (pos, vel) = index.eval_leb2_body(&file, 0, 0.25).unwrap();
    assert_eq!(pos, vec![5.0, -0.5, 0.0]);
    assert_eq!(vel, vec![0.0, 2.0, 0.0]);
    let (pos, vel) = index.eval_leb2_body(&file, 0, 1.0).unwrap();
    assert_eq!(pos, vec![7.0, -1.0, 0.0]);
    assert_eq!(vel, vec![0.0, 4.0, 0.0]);
    let (pos, vel) = index.eval_leb2_body(&file, 0, 1.5).unwrap();
    assert_eq!(pos, vec![7.0, 1.0, 0.0]);
    assert_eq!(vel, vec![0.0, 4.0, 0.0]);
    // Rounded-bound twin at eval level: chunk 0's JD bound sits one
    // ulp below the seam. The epoch just below the seam still belongs
    // to seg 0, so chunk 0's constant-5 channel must come back —
    // never chunk 1's constant-7.
    let rounded = f64::from_bits(1.0f64.to_bits() - 1);
    let file = leb2_v2_distinct_image(rounded);
    let index = FileIndex::decode(&file).expect("rounded image must decode");
    let (pos, vel) = index.eval_leb2_body(&file, 0, rounded).unwrap();
    assert_eq!(pos[0], 5.0);
    assert_eq!(vel, vec![0.0, 2.0, 0.0]);
    assert!(pos[1] < 1.0 && pos[1] > 0.999_999_999_999);
}

#[test]
fn leb2_cached_eval_matches_direct_bitwise_with_hits() {
    // The decoded-chunk memo must serve the same words as
    // the direct decode on both LEB2 layouts, across chunk
    // boundaries and on repeats (the hit path).
    let cache = LebChunkCache::new();
    // v2 chunked image: chunk 0 serves seg 0, chunk 1 serves seg 1.
    let file = leb2_v2_distinct_image(1.0);
    let index = FileIndex::decode(&file).expect("LEB2 v2 image must decode");
    // Epochs spanning both chunks, each evaluated twice so the
    // second pass exercises the memo hit path.
    let jds = [0.25, 1.5, 0.25, 1.0, 1.5, 0.25];
    for &jd in &jds {
        let (pos, vel) = index
            .eval_leb2_body_cached(&file, 0, jd, &cache)
            .expect("cached eval serves");
        let (raw_pos, raw_vel) = index
            .eval_leb2_body(&file, 0, jd)
            .expect("direct eval serves");
        assert_eq!(pos, raw_pos, "cached pos words at jd {jd}");
        assert_eq!(vel, raw_vel, "cached vel words at jd {jd}");
    }
    // Two distinct chunks decoded (2 misses), four repeats served
    // from the memo (4 hits), two resident entries.
    assert_eq!(cache.stats(), (4, 2, 2));
    // Split-epoch twin shares the memo path (collapsed-word keying):
    // same words as the direct two-part eval, all memo hits.
    for &jd in &[0.25, 1.5] {
        let epoch = TwoPartEpoch::split(jd).unwrap();
        let (pos, vel) = index
            .eval_leb2_body_two_part_cached(&file, 0, epoch, &cache)
            .expect("cached two-part serves");
        let (raw_pos, raw_vel) = index
            .eval_leb2_body_two_part(&file, 0, epoch)
            .expect("direct two-part serves");
        assert_eq!(pos, raw_pos, "cached two-part pos words at jd {jd}");
        assert_eq!(vel, raw_vel, "cached two-part vel words at jd {jd}");
    }
    assert_eq!(cache.stats(), (6, 2, 2));
    // v1 monolith twin: same words, one memo entry.
    let cache_v1 = LebChunkCache::new();
    let file_v1 = leb2_v1_distinct_image();
    let index_v1 = FileIndex::decode(&file_v1).expect("LEB2 v1 image must decode");
    for &jd in &[0.25, 1.5, 0.25] {
        let (pos, vel) = index_v1
            .eval_leb2_body_cached(&file_v1, 0, jd, &cache_v1)
            .expect("cached v1 eval serves");
        let (raw_pos, raw_vel) = index_v1
            .eval_leb2_body(&file_v1, 0, jd)
            .expect("direct v1 serves");
        assert_eq!(pos, raw_pos, "cached v1 pos words at jd {jd}");
        assert_eq!(vel, raw_vel, "cached v1 vel words at jd {jd}");
    }
    assert_eq!(cache_v1.stats(), (2, 1, 1));
    // Error twin: an out-of-range epoch refuses before touching the
    // memo (no miss counted, nothing memoized).
    assert!(index.eval_leb2_body_cached(&file, 0, 5.0, &cache).is_err());
    assert_eq!(cache.stats(), (6, 2, 2));
}

#[test]
fn composite_memo_serves_direct_words_bitwise_with_hits() {
    // The per-constituent decoded-chunk memo must serve the same words
    // as the direct decode through the composite and tiered dispatch,
    // with LEB1 constituents bypassing the memo and shadowed readers
    // never leaking into served words.
    let leb2 = leb2_v2_distinct_image(1.0);
    let index_leb2 = FileIndex::decode(&leb2).expect("LEB2 image must decode");
    let leb1 = leb1_const_image(1, 9.0, 0.0, 4.0);
    let index_leb1 = FileIndex::decode(&leb1).expect("LEB1 image must decode");
    let composite = CompositeIndex::new(vec![
        (index_leb2, leb2.as_slice()),
        (index_leb1, leb1.as_slice()),
    ])
    .expect("composite must build");
    assert_eq!(composite.chunk_cache_stats(), (0, 0, 0));
    // Body 0 (LEB2, chunks 0 and 1): two passes so the second
    // exercises the memo hit path; every word bitwise vs direct.
    for &jd in &[0.25, 1.5, 0.25, 1.5] {
        let (pos, vel) = composite.eval_body_any(0, jd).expect("composite serves");
        let direct = FileIndex::decode(&leb2)
            .unwrap()
            .eval_leb2_body(&leb2, 0, jd)
            .expect("direct serves");
        assert_eq!((pos, vel), direct, "composite words at jd {jd}");
    }
    // Two distinct chunks decoded (2 misses), two repeats served
    // from the memo (2 hits), two resident entries.
    assert_eq!(composite.chunk_cache_stats(), (2, 2, 2));
    // Split-epoch twin shares the serving constituent's memo.
    for &jd in &[0.25, 1.5] {
        let epoch = TwoPartEpoch::split(jd).unwrap();
        let (pos, vel) = composite
            .eval_body_any_two_part(0, epoch)
            .expect("composite two-part serves");
        let direct = FileIndex::decode(&leb2)
            .unwrap()
            .eval_leb2_body_two_part(&leb2, 0, epoch)
            .expect("direct two-part serves");
        assert_eq!((pos, vel), direct, "composite two-part words at jd {jd}");
    }
    assert_eq!(composite.chunk_cache_stats(), (4, 2, 2));
    // Body 1 (LEB1) bypasses the memo: served words, stats untouched.
    assert_eq!(
        composite.eval_body_any(1, 1.0).unwrap(),
        (vec![9.0, 9.0, 9.0], vec![0.0, 0.0, 0.0])
    );
    assert_eq!(composite.chunk_cache_stats(), (4, 2, 2));
    // First-wins isolation: a leading LEB1 body-0 shadows the LEB2
    // reader, whose memo stays empty — shadowed bytes never leak.
    let shadow_leb1 = leb1_const_image(0, 5.0, 0.0, 2.0);
    let shadow_index = FileIndex::decode(&shadow_leb1).expect("shadow must decode");
    let shadowed = CompositeIndex::new(vec![
        (shadow_index, shadow_leb1.as_slice()),
        (FileIndex::decode(&leb2).unwrap(), leb2.as_slice()),
    ])
    .expect("shadowed composite must build");
    assert_eq!(
        shadowed.eval_body_any(0, 1.0).unwrap(),
        (vec![5.0, 5.0, 5.0], vec![0.0, 0.0, 0.0])
    );
    assert_eq!(shadowed.chunk_cache_stats(), (0, 0, 0));
    // Tiered dispatch inherits the memo through its composite.
    let base = CompositeIndex::new(vec![(FileIndex::decode(&leb2).unwrap(), leb2.as_slice())])
        .expect("tiered base must build");
    let tiered = TieredIndex::new(vec![(LebTier::Base, base)]).expect("tiered must build");
    for &jd in &[0.25, 0.25] {
        let (pos, vel) = tiered.eval_body_any(0, jd).expect("tiered serves");
        let direct = FileIndex::decode(&leb2)
            .unwrap()
            .eval_leb2_body(&leb2, 0, jd)
            .expect("direct serves");
        assert_eq!((pos, vel), direct, "tiered words at jd {jd}");
    }
    assert_eq!(tiered.chunk_cache_stats(), (1, 1, 1));
    // Error twin: an out-of-range epoch refuses before touching the
    // memo (no miss counted, nothing memoized).
    assert!(composite.eval_body_any(0, 5.0).is_err());
    assert_eq!(composite.chunk_cache_stats(), (4, 2, 2));
}

#[test]
fn leb2_eval_rejects_bad_inputs() {
    let file = leb2_v1_distinct_image();
    let index = FileIndex::decode(&file).expect("LEB2 v1 image must decode");
    for jd in [-0.5, 2.0 + 1e-9, 100.0] {
        assert_eq!(
            index.eval_leb2_body(&file, 0, jd).unwrap_err(),
            Error::LebOutOfRange {
                body_id: 0,
                jd,
                start: 0.0,
                end: 2.0,
            }
        );
    }
    assert_eq!(
        index.eval_leb2_body(&file, 99, 1.0).unwrap_err(),
        Error::LebBodyNotFound { body_id: 99 }
    );
    assert!(matches!(
        index.eval_leb2_body(&file, 0, f64::NAN).unwrap_err(),
        Error::NonFiniteArgument { .. }
    ));
    // LEB1 on the LEB2 path is a usage error, never silent raw reads.
    let leb1 = minimal_leb1();
    let leb1_index = FileIndex::decode(&leb1).expect("LEB1 fixture decodes");
    assert!(matches!(
        leb1_index.eval_leb2_body(&leb1, 0, 1.0).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    // Gap twin: chunk 1 starts at seg 5, so seg 1 near hint 1 is
    // uncovered — corruption, never a neighbouring evaluation.
    let raw = raw_two_segment_distinct_blob();
    let blob0 = encode_frame(&raw[..48], 1, 1, 3);
    let blob1 = encode_frame(&raw[48..], 1, 1, 3);
    let mut buf = pack_chunk_header(2, 3652.5);
    buf.extend_from_slice(&pack_chunk_entry(
        0.0,
        1.0,
        500,
        blob0.len() as u64,
        48,
        0,
        1,
    ));
    buf.extend_from_slice(&pack_chunk_entry(
        1.0,
        2.0,
        600,
        blob1.len() as u64,
        48,
        5,
        1,
    ));
    let mut file = pack_header(&MAGIC_LEB2, 2, 1, 1, 0.0, 2.0, 0.0, 1);
    let body_at = file.len() + SECTION_DIR_SIZE;
    let index_at = body_at + LEB2_BODY_ENTRY_SIZE;
    file.extend_from_slice(&pack_section(
        0,
        body_at as u64,
        LEB2_BODY_ENTRY_SIZE as u64,
    ));
    file.extend_from_slice(&pack_body_leb2(
        0,
        0,
        2,
        0.0,
        2.0,
        1.0,
        1,
        3,
        index_at as u64,
        0,
        96,
    ));
    file.extend_from_slice(&buf);
    // Pad so data_offset is inside the image.
    while file.len() < index_at {
        file.push(0);
    }
    let index = FileIndex::decode(&file).expect("gap image decodes");
    assert!(matches!(
        index.eval_leb2_body(&file, 0, 1.5).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    // Empty chunk index twin.
    let mut empty = pack_header(&MAGIC_LEB2, 2, 1, 1, 0.0, 2.0, 0.0, 1);
    let body_at = empty.len() + SECTION_DIR_SIZE;
    let index_at = body_at + LEB2_BODY_ENTRY_SIZE;
    empty.extend_from_slice(&pack_section(
        0,
        body_at as u64,
        LEB2_BODY_ENTRY_SIZE as u64,
    ));
    empty.extend_from_slice(&pack_body_leb2(
        0,
        0,
        2,
        0.0,
        2.0,
        1.0,
        1,
        3,
        index_at as u64,
        0,
        96,
    ));
    empty.extend_from_slice(&pack_chunk_header(0, 3652.5));
    let index = FileIndex::decode(&empty).expect("empty-chunk image decodes");
    assert!(matches!(
        index.eval_leb2_body(&empty, 0, 1.0).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
}

#[test]
fn degree_two_velocity_scale_is_exact() {
    // Single body, one segment, degree 2, one component [0, 0, 1]
    // (T2 = 2*tau^2 - 1), interval 2.0 d so the velocity scale is
    // exactly 1. At tau 0.5 (jd 1.5): value -0.5, d/dtau 2.0.
    let mut file = pack_header(&MAGIC_LEB1, 1, 1, 1, 0.0, 2.0, 0.0, 0);
    let body_at = file.len() + SECTION_DIR_SIZE;
    let blob_at = body_at + LEB1_BODY_ENTRY_SIZE;
    file.extend_from_slice(&pack_section(
        0,
        body_at as u64,
        LEB1_BODY_ENTRY_SIZE as u64,
    ));
    file.extend_from_slice(&pack_body_leb1(
        0,
        0,
        1,
        0.0,
        2.0,
        2.0,
        2,
        1,
        blob_at as u64,
    ));
    for c in [0.0, 0.0, 1.0] {
        push_f64(&mut file, c);
    }
    let index = FileIndex::decode(&file).expect("degree-2 fixture decodes");
    let (pos, vel) = index.eval_body(&file, 0, 1.5).unwrap();
    assert_eq!(pos, vec![-0.5]);
    assert_eq!(vel, vec![2.0]);
    let (pos, vel) = index.eval_body(&file, 0, 1.0).unwrap();
    assert_eq!(pos, vec![-1.0]);
    assert_eq!(vel, vec![0.0]);
}

/// Aux fixture: LEB1 header with an empty section 0 (required by the
/// index decoder) plus caller-supplied aux sections.
fn aux_image(extra_sections: Vec<(u32, Vec<u8>)>) -> Vec<u8> {
    let n = 1 + extra_sections.len() as u32;
    let mut file = pack_header(&MAGIC_LEB1, 1, n, 0, 0.0, 16.0, 0.0, 0);
    let mut at = file.len() + n as usize * SECTION_DIR_SIZE;
    let mut payloads = Vec::new();
    file.extend_from_slice(&pack_section(0, at as u64, 0));
    for (id, payload) in &extra_sections {
        file.extend_from_slice(&pack_section(*id, at as u64, payload.len() as u64));
        payloads.push(payload.clone());
        at += payload.len();
    }
    for payload in payloads {
        file.extend_from_slice(&payload);
    }
    file
}

fn nutation_payload(segs: u32, interval: f64) -> Vec<u8> {
    let mut v = Vec::with_capacity(NUTATION_HEADER_SIZE + 32);
    push_f64(&mut v, 0.0);
    push_f64(&mut v, 16.0);
    push_f64(&mut v, interval);
    push_u32(&mut v, 1);
    push_u32(&mut v, 2);
    push_u32(&mut v, segs);
    push_u32(&mut v, 0);
    // dpsi = [1.0, 0.5], deps = [2.0, -1.0]: at tau 0 -> (1, 2).
    for c in [1.0, 0.5, 2.0, -1.0] {
        push_f64(&mut v, c);
    }
    v
}

#[test]
fn nutation_eval_matches_hand_values() {
    let file = aux_image(vec![(SECTION_NUTATION, nutation_payload(1, 16.0))]);
    let index = FileIndex::decode(&file).expect("nutation image must decode");
    // Mid-segment: tau 0 -> the constant terms exactly.
    assert_eq!(index.eval_nutation(&file, 8.0).unwrap(), (1.0, 2.0));
    // Segment start: tau -1 -> (0.5, 3.0).
    assert_eq!(index.eval_nutation(&file, 0.0).unwrap(), (0.5, 3.0));
    assert_eq!(
        index.eval_nutation(&file, 17.0).unwrap_err(),
        Error::LebOutOfRange {
            body_id: -1,
            jd: 17.0,
            start: 0.0,
            end: 16.0,
        }
    );
    assert!(matches!(
        index.eval_nutation(&file, f64::NAN).unwrap_err(),
        Error::NonFiniteArgument { .. }
    ));
    // `--skip-aux` shape: segment_count 0 decodes as absent.
    let file = aux_image(vec![(SECTION_NUTATION, nutation_payload(0, 16.0))]);
    let index = FileIndex::decode(&file).expect("absent nutation decodes");
    assert!(matches!(
        index.eval_nutation(&file, 8.0).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    // Zero interval with a nonzero count is corruption, not division.
    let file = aux_image(vec![(SECTION_NUTATION, nutation_payload(1, 0.0))]);
    let index = FileIndex::decode(&file).expect("zero-interval header decodes");
    assert!(matches!(
        index.eval_nutation(&file, 8.0).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    // Missing section entirely.
    let file = aux_image(vec![]);
    let index = FileIndex::decode(&file).expect("bare image decodes");
    assert!(matches!(
        index.eval_nutation(&file, 8.0).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
}

fn delta_t_payload(pairs: &[(f64, f64)]) -> Vec<u8> {
    let mut v = Vec::with_capacity(DELTA_T_HEADER_SIZE + pairs.len() * DELTA_T_ENTRY_SIZE);
    push_u32(&mut v, pairs.len() as u32);
    push_u32(&mut v, 0);
    for (jd, dt) in pairs {
        push_f64(&mut v, *jd);
        push_f64(&mut v, *dt);
    }
    v
}

#[test]
fn delta_t_clamps_bisects_and_interpolates() {
    let file = aux_image(vec![(
        SECTION_DELTA_T,
        delta_t_payload(&[(0.0, 1.0), (10.0, 3.0)]),
    )]);
    let index = FileIndex::decode(&file).expect("Delta-T image must decode");
    assert_eq!(index.delta_t(&file, -5.0).unwrap(), 1.0);
    assert_eq!(index.delta_t(&file, 0.0).unwrap(), 1.0);
    assert_eq!(index.delta_t(&file, 5.0).unwrap(), 2.0);
    assert_eq!(index.delta_t(&file, 10.0).unwrap(), 3.0);
    assert_eq!(index.delta_t(&file, 20.0).unwrap(), 3.0);
    assert!(matches!(
        index.delta_t(&file, f64::NAN).unwrap_err(),
        Error::NonFiniteArgument { .. }
    ));
    // Empty table: absent data, never a silent zero.
    let file = aux_image(vec![(SECTION_DELTA_T, delta_t_payload(&[]))]);
    let index = FileIndex::decode(&file).expect("empty Delta-T decodes");
    assert!(matches!(
        index.delta_t(&file, 5.0).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    // Missing section entirely.
    let file = aux_image(vec![]);
    let index = FileIndex::decode(&file).expect("bare image decodes");
    assert!(matches!(
        index.delta_t(&file, 5.0).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
}

fn star_payload(ids: &[i32]) -> Vec<u8> {
    let mut v = Vec::new();
    for id in ids {
        push_i32(&mut v, *id);
        for x in [10.5, -20.25, 0.001, -0.002, 0.05, 11.0, 2.5] {
            push_f64(&mut v, x);
        }
        v.extend_from_slice(b"RSVD");
    }
    v
}

#[test]
fn star_lookup_finds_and_misses() {
    let file = aux_image(vec![(SECTION_STARS, star_payload(&[7, 9]))]);
    let index = FileIndex::decode(&file).expect("star image must decode");
    let star = index.get_star(&file, 7).expect("star 7 present");
    assert_eq!(star.star_id(), 7);
    assert_eq!(star.ra_j2000_deg(), 10.5);
    assert_eq!(star.magnitude(), 2.5);
    assert_eq!(index.get_star(&file, 9).unwrap().star_id(), 9);
    assert_eq!(
        index.get_star(&file, 8).unwrap_err(),
        Error::LebStarNotFound { star_id: 8 }
    );
    // Short tail ignored: 128 + 63 bytes still hold 2 stars.
    let mut long = file.clone();
    long.extend_from_slice(&[0u8; 63]);
    let index = FileIndex::decode(&long).expect("padded image decodes");
    assert_eq!(index.get_star(&long, 9).unwrap().star_id(), 9);
    // Missing section entirely.
    let file = aux_image(vec![]);
    let index = FileIndex::decode(&file).expect("bare image decodes");
    assert_eq!(
        index.get_star(&file, 7).unwrap_err(),
        Error::LebStarNotFound { star_id: 7 }
    );
}

#[test]
fn corrupt_reports_what_and_offset() {
    let err = Error::LebCorrupt {
        what: "LEB header".to_string(),
        offset: 63,
    };
    assert_eq!(format!("{err}"), "corrupt LEB LEB header at byte 63");
    let err = Error::LebOutOfRange {
        body_id: 3,
        jd: -1.0,
        start: 0.0,
        end: 2.0,
    };
    assert_eq!(format!("{err}"), "JD -1 outside range [0, 2] for body 3");
}

/// JD-scale one-body LEB1 for the split-epoch tests: two
/// 32-day segments over `[2451545, 2451609]` (seam and midpoints are
/// exact words), degree 1, one ICRS component; seg0 value is
/// `100 + 10*tau`, seg1 `200 - 20*tau`.
fn jd_scale_leb1() -> Vec<u8> {
    let (start, interval) = (2_451_545.0, 32.0);
    let mut file = pack_header(&MAGIC_LEB1, 1, 1, 1, start, start + 64.0, 0.0, 0);
    let dir_at = file.len();
    let body_at = dir_at + SECTION_DIR_SIZE;
    let blob_at = body_at + LEB1_BODY_ENTRY_SIZE;
    file.extend_from_slice(&pack_section(0, body_at as u64, 52));
    file.extend_from_slice(&pack_body_leb1(
        0,
        0,
        2,
        start,
        start + 64.0,
        interval,
        1,
        1,
        blob_at as u64,
    ));
    for c in [100.0, 10.0, 200.0, -20.0] {
        push_f64(&mut file, c);
    }
    debug_assert_eq!(file.len(), blob_at + 32);
    file
}

/// JD-scale nutation payload: one 32-day segment over
/// `[2451545, 2451577]`, degree 1, two components with the same
/// hand coefficients as [`nutation_payload`].
fn jd_nutation_payload() -> Vec<u8> {
    let mut v = Vec::with_capacity(NUTATION_HEADER_SIZE + 32);
    push_f64(&mut v, 2_451_545.0);
    push_f64(&mut v, 2_451_577.0);
    push_f64(&mut v, 32.0);
    push_u32(&mut v, 1);
    push_u32(&mut v, 2);
    push_u32(&mut v, 1);
    push_u32(&mut v, 0);
    // dpsi = [1.0, 0.5], deps = [2.0, -1.0]: at tau 0 -> (1, 2).
    for c in [1.0, 0.5, 2.0, -1.0] {
        push_f64(&mut v, c);
    }
    v
}

#[test]
fn two_part_selection_matches_scalar_on_collapsed_words() {
    // Ownership pins across a JD-scale grid: segment starts, the
    // interior seam (owned right), the global end (last segment),
    // adjacent float words around the seam, and mid-segment points.
    let seam: f64 = 2_451_577.0;
    let epochs = [
        2_451_545.0,
        2_451_550.25,
        2_451_561.0,
        seam.next_down(),
        seam,
        seam.next_up(),
        2_451_593.0,
        2_451_600.5,
        2_451_609.0,
    ];
    for v in epochs {
        let (s_scalar, tau_scalar) = select_segment(2_451_545.0, 32.0, 2, v).unwrap();
        let epoch = TwoPartEpoch::split(v).unwrap();
        assert_eq!(epoch.value().to_bits(), v.to_bits());
        let (s_split, tau_split) = select_segment_two_part(2_451_545.0, 32.0, 2, epoch).unwrap();
        assert_eq!(s_split, s_scalar, "ownership diverged at {v}");
        // Exact splits carry no tail: the two subtraction orders
        // agree to rounding noise. A priori bound: each displacement
        // is O(16) (ulp ~= 3.6e-15), tau scales by 2/32, so 4 eps of
        // tau carries several ulps of margin — fixed, not fitted.
        assert!(
            (tau_split.value() - tau_scalar.value()).abs() <= 4.0 * f64::EPSILON,
            "tau diverged at {v}"
        );
    }
    // Exact seam owned right; global end owned by the last segment.
    let (s, _) =
        select_segment_two_part(2_451_545.0, 32.0, 2, TwoPartEpoch::split(seam).unwrap()).unwrap();
    assert_eq!(s, 1);
    let (s, _) = select_segment_two_part(
        2_451_545.0,
        32.0,
        2,
        TwoPartEpoch::split(2_451_609.0).unwrap(),
    )
    .unwrap();
    assert_eq!(s, 1);
    // Degenerate grid is corruption, as in the scalar selector.
    assert!(matches!(
        select_segment_two_part(
            2_451_545.0,
            32.0,
            0,
            TwoPartEpoch::split(2_451_545.0).unwrap()
        )
        .unwrap_err(),
        Error::LebCorrupt { .. }
    ));
}

#[test]
fn two_part_tau_preserves_sub_ulp_tail() {
    let file = jd_scale_leb1();
    let index = FileIndex::decode(&file).expect("JD fixture must decode");
    // 1e-10 d is below half a JD ulp near 2451561 (ulp ~= 4.7e-10):
    // plain f64 addition loses it, the split epoch keeps it.
    assert_eq!(2_451_561.0_f64 + 1e-10, 2_451_561.0);
    let epoch = TwoPartEpoch::split(2_451_561.0)
        .unwrap()
        .add_offset(1e-10)
        .unwrap();
    let (seg, tau) = select_segment_two_part(2_451_545.0, 32.0, 2, epoch).unwrap();
    assert_eq!(seg, 0);
    // Expected tau = 2e-10/32 = 6.25e-12; the ordered difference is
    // exact through Sterbenz-adjacent magnitudes, so 2 ulps of the
    // expectation is generous, not fitted.
    let expect_tau = 6.25e-12;
    assert!(
        (tau.value() - expect_tau).abs() <= 2.0 * f64::EPSILON * expect_tau.abs(),
        "tail lost: tau = {}",
        tau.value()
    );
    // End to end: seg0 value is 100 + 10*tau; the scalar path on the
    // collapsed word gives exactly 100, the split path keeps the
    // 6.25e-11 tail in a distinct word.
    let (pos_split, vel_split) = index.eval_body_two_part(&file, 0, epoch).unwrap();
    let (pos_scalar, _) = index.eval_body(&file, 0, epoch.value()).unwrap();
    assert_eq!(pos_scalar[0], 100.0);
    assert_ne!(pos_split[0].to_bits(), 100.0f64.to_bits());
    // A priori bound: |c1| * tau noise + one Clenshaw rounding at
    // |v| ~= 100 (ulp(100) ~= 1.4e-14); 1e-12 carries ~70x margin.
    assert!(
        (pos_split[0] - (100.0 + 10.0 * expect_tau)).abs() <= 1e-12,
        "split value = {}",
        pos_split[0]
    );
    // Degree-1 x-derivative is exactly c1 = 10, so the velocity is
    // the exact word 10 * (2/32) = 0.625 on both paths.
    assert_eq!(vel_split[0], 0.625);
    let (_, vel_scalar) = index.eval_body(&file, 0, epoch.value()).unwrap();
    assert_eq!(vel_scalar[0], 0.625);
}

#[test]
fn retarded_epoch_steps_across_seam_with_re_split() {
    // The seam word is owned right; retard by 1e-9 d (above half a
    // JD ulp, so the collapsed word moves below the seam) and the
    // re-split retarded epoch owns left — exactly what the scalar
    // path sees on the collapsed word.
    let seam = TwoPartEpoch::split(2_451_577.0).unwrap();
    let retarded = retard_epoch(seam, 1e-9).unwrap();
    assert!(retarded.value() < 2_451_577.0);
    let (s_ret, _) = select_segment_two_part(2_451_545.0, 32.0, 2, retarded).unwrap();
    assert_eq!(s_ret, 0);
    let (s_scalar, _) = select_segment(2_451_545.0, 32.0, 2, retarded.value()).unwrap();
    assert_eq!(s_scalar, 0);
    // A sub-ulp retardation vanishes under f64 subtraction (pinned)
    // but survives in `lo`: the two-part diff recovers 1e-10 d with
    // a Sterbenz-exact leading difference (1e-18 is ~1e-8 relative,
    // generous).
    assert_eq!(2_451_577.0_f64 - 1e-10, 2_451_577.0);
    let kept = retard_epoch(seam, 1e-10).unwrap();
    assert_eq!(kept.value(), 2_451_577.0);
    assert!(
        (seam.diff(kept) - 1e-10).abs() <= 1e-18,
        "retarded diff = {}",
        seam.diff(kept)
    );
    // A negative delay advances the epoch; a non-finite delay is
    // refused rather than smuggled into the parts.
    let advanced = retard_epoch(seam, -1.0).unwrap();
    assert_eq!(advanced.value(), 2_451_578.0);
    assert!(matches!(
        retard_epoch(seam, f64::NAN).unwrap_err(),
        Error::NonFiniteArgument { .. }
    ));
}

#[test]
fn split_and_scalar_body_evals_agree_within_stated_bound() {
    let file = jd_scale_leb1();
    let index = FileIndex::decode(&file).expect("JD fixture must decode");
    // Sweep both segments including the seam and both ends. Exact
    // splits carry no tail, so this pins selection identity plus
    // subtraction-order noise. A priori bound: |c1| <= 20 times a
    // few-ulp tau difference (~1e-15) plus shared-kernel rounding at
    // |v| <= 220 (ulp ~= 2.8e-14); 1e-12 is all margin, not fitted.
    for v in [
        2_451_545.0,
        2_451_550.25,
        2_451_576.0,
        2_451_577.0,
        2_451_578.0,
        2_451_593.0,
        2_451_609.0,
    ] {
        let epoch = TwoPartEpoch::split(v).unwrap();
        let (p_split, v_split) = index.eval_body_two_part(&file, 0, epoch).unwrap();
        let (p_scalar, v_scalar) = index.eval_body(&file, 0, v).unwrap();
        assert!(
            (p_split[0] - p_scalar[0]).abs() <= 1e-12,
            "pos diverged at {v}"
        );
        assert!(
            (v_split[0] - v_scalar[0]).abs() <= 1e-12,
            "vel diverged at {v}"
        );
    }
    // Tailed splits (sub-ulp `lo`) select the same segment as the
    // scalar word and stay within the same bound of its answer.
    for lo in [1e-10, -1e-10, 2e-10] {
        let epoch = TwoPartEpoch::split(2_451_561.0)
            .unwrap()
            .add_offset(lo)
            .unwrap();
        let (p_split, _) = index.eval_body_two_part(&file, 0, epoch).unwrap();
        let (p_scalar, _) = index.eval_body(&file, 0, epoch.value()).unwrap();
        assert_eq!(p_scalar[0], 100.0);
        assert!(
            (p_split[0] - 100.0 - 10.0 * (2.0 * lo / 32.0)).abs() <= 1e-12,
            "tailed split diverged for lo = {lo}"
        );
    }
}

#[test]
fn nutation_two_part_matches_scalar_and_keeps_tail() {
    let file = aux_image(vec![(SECTION_NUTATION, jd_nutation_payload())]);
    let index = FileIndex::decode(&file).expect("JD nutation image must decode");
    // Segment midpoint: tau is exactly 0 on both paths.
    let mid = TwoPartEpoch::split(2_451_561.0).unwrap();
    assert_eq!(
        index.eval_nutation_two_part(&file, mid).unwrap(),
        (1.0, 2.0)
    );
    assert_eq!(index.eval_nutation(&file, 2_451_561.0).unwrap(), (1.0, 2.0));
    // Tailed split: scalar collapses to the midpoint (1, 2); the
    // split path keeps dpsi = 1 + 0.5*tau, deps = 2 - 1.0*tau with
    // tau = 6.25e-12. Bound: one Clenshaw rounding at |v| <= 2.
    let epoch = mid.add_offset(1e-10).unwrap();
    assert_eq!(
        index.eval_nutation(&file, epoch.value()).unwrap(),
        (1.0, 2.0)
    );
    let (dpsi, deps) = index.eval_nutation_two_part(&file, epoch).unwrap();
    assert_ne!(dpsi.to_bits(), 1.0f64.to_bits());
    assert!(
        (dpsi - (1.0 + 0.5 * 6.25e-12)).abs() <= 1e-15,
        "dpsi = {dpsi}"
    );
    assert!((deps - (2.0 - 6.25e-12)).abs() <= 1e-15, "deps = {deps}");
    // Range errors key on the collapsed word, same shape as scalar.
    assert_eq!(
        index
            .eval_nutation_two_part(&file, TwoPartEpoch::split(2_451_578.0).unwrap())
            .unwrap_err(),
        Error::LebOutOfRange {
            body_id: -1,
            jd: 2_451_578.0,
            start: 2_451_545.0,
            end: 2_451_577.0,
        }
    );
}

#[test]
fn leb2_two_part_matches_scalar_on_both_layouts() {
    // Monolith (v1) and chunked (v2) images over the same distinct
    // blob: the split path shares selection, hint, and size checks
    // with the scalar path, so agreement holds on both layouts.
    // Bound: hand coefficients are O(1..7); same derivation class
    // as the LEB1 sweep above with wider margin at 1e-12.
    for file in [leb2_v1_distinct_image(), leb2_v2_distinct_image(1.0)] {
        let index = FileIndex::decode(&file).expect("LEB2 image must decode");
        for v in [0.25, 1.0, 1.5, 2.0] {
            let epoch = TwoPartEpoch::split(v).unwrap();
            let (p_split, v_split) = index.eval_leb2_body_two_part(&file, 0, epoch).unwrap();
            let (p_scalar, v_scalar) = index.eval_leb2_body(&file, 0, v).unwrap();
            assert_eq!(p_split.len(), p_scalar.len());
            for (a, b) in p_split.iter().zip(p_scalar.iter()) {
                assert!((a - b).abs() <= 1e-12, "pos diverged at {v}");
            }
            for (a, b) in v_split.iter().zip(v_scalar.iter()) {
                assert!((a - b).abs() <= 1e-12, "vel diverged at {v}");
            }
        }
    }
}

#[test]
fn two_part_eval_rejects_bad_inputs() {
    let file = jd_scale_leb1();
    let index = FileIndex::decode(&file).expect("JD fixture must decode");
    let good = TwoPartEpoch::split(2_451_561.0).unwrap();
    // Hostile directly-constructed parts collapse non-finite.
    for hostile in [
        TwoPartEpoch {
            hi: f64::INFINITY,
            lo: 0.0,
        },
        TwoPartEpoch {
            hi: f64::NAN,
            lo: 0.0,
        },
        TwoPartEpoch {
            hi: 2_451_561.0,
            lo: f64::INFINITY,
        },
    ] {
        assert!(matches!(
            index.eval_body_two_part(&file, 0, hostile).unwrap_err(),
            Error::NonFiniteArgument { .. }
        ));
        assert!(matches!(
            index
                .eval_leb2_body_two_part(&file, 0, hostile)
                .unwrap_err(),
            // LEB1 image on the LEB2 path is corruption regardless.
            Error::NonFiniteArgument { .. } | Error::LebCorrupt { .. }
        ));
    }
    // Out of range keys on the collapsed word, same shape as scalar.
    assert_eq!(
        index
            .eval_body_two_part(&file, 0, TwoPartEpoch::split(2_451_610.0).unwrap())
            .unwrap_err(),
        Error::LebOutOfRange {
            body_id: 0,
            jd: 2_451_610.0,
            start: 2_451_545.0,
            end: 2_451_609.0,
        }
    );
    // Unknown body.
    assert!(matches!(
        index.eval_body_two_part(&file, 7, good).unwrap_err(),
        Error::LebBodyNotFound { .. }
    ));
    // Format gates hold in both directions on the two-part paths.
    assert!(matches!(
        index.eval_leb2_body_two_part(&file, 0, good).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    let (leb2_bytes, _) = leb2_v1_image();
    let leb2 = FileIndex::decode(&leb2_bytes).expect("LEB2 image must decode");
    let half = TwoPartEpoch::split(0.5).unwrap();
    assert!(matches!(
        leb2.eval_body_two_part(&leb2_bytes, 0, half).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
}

// LEB-08 fixtures: single-segment degree-0 LEB1 image with one body
// holding a constant value over [jd_start, jd_end]. Evaluates to
// pos = [value; 3], vel = [0; 3] everywhere inside the range.
fn leb1_const_image(body_id: i32, value: f64, jd_start: f64, jd_end: f64) -> Vec<u8> {
    let mut file = pack_header(&MAGIC_LEB1, 1, 1, 1, jd_start, jd_end, 0.0, 0);
    let dir_at = file.len();
    let body_at = dir_at + SECTION_DIR_SIZE;
    let blob_at = body_at + LEB1_BODY_ENTRY_SIZE;
    file.extend_from_slice(&pack_section(
        0,
        body_at as u64,
        LEB1_BODY_ENTRY_SIZE as u64,
    ));
    file.extend_from_slice(&pack_body_leb1(
        body_id,
        COORD_ICRS_BARY,
        1,
        jd_start,
        jd_end,
        jd_end - jd_start,
        0,
        3,
        blob_at as u64,
    ));
    for _ in 0..3 {
        push_f64(&mut file, value);
    }
    file
}

#[test]
fn group_and_tier_constants_match_release_partition() {
    // leb_groups.py pins,.
    assert_eq!(
        LEB1_GROUPS,
        ["planets", "asteroids", "exotics", "analytical"]
    );
    assert!(LEB1_COMPANION_GROUPS.is_empty());
    assert_eq!(LEB2_GROUPS, ["core", "asteroids", "exotics", "apogee"]);
    assert!(GROUP_SUFFIXES.contains(&"uranians"));
    assert_eq!(RETIRED_GROUP_SUFFIXES, ["uranians"]);
    assert_eq!(
        LebTier::PRIORITY,
        [LebTier::Base, LebTier::Medium, LebTier::Extended]
    );
    assert_eq!(LebTier::Base.name(), "base");
    assert_eq!(LebTier::parse("medium"), Some(LebTier::Medium));
    assert_eq!(LebTier::parse("core"), None);
}

#[test]
fn filename_policy_matches_composite_heuristics() {
    // Group-suffix classification on the extension-stripped name.
    assert_eq!(file_group_suffix("base_core.leb2"), Some("core"));
    assert_eq!(file_group_suffix("myset_base_apogee.leb"), Some("apogee"));
    assert_eq!(file_group_suffix("custom.leb"), None);
    assert_eq!(file_group_suffix("core.leb2"), Some("core"));
    // Tier extraction: exact prefix, custom-prefix token scan, bare
    // known-tier token, and untiered names.
    assert_eq!(file_tier("base_core.leb2"), Some(LebTier::Base));
    assert_eq!(file_tier("myset_base_core.leb2"), Some(LebTier::Base));
    assert_eq!(file_tier("ephemeris_base.leb"), Some(LebTier::Base));
    assert_eq!(file_tier("medium_apogee.leb2"), Some(LebTier::Medium));
    assert_eq!(file_tier("custom.leb"), None);
    // Retired companions classify but must never attach.
    assert!(is_retired_group_file("base_uranians.leb2"));
    assert!(!is_retired_group_file("base_core.leb2"));
    // Companion prefix: group files share it, merged files have none.
    assert_eq!(group_prefix("base_core.leb2"), Some("base".to_string()));
    assert_eq!(
        group_prefix("myset_base_core.leb2"),
        Some("myset_base".to_string())
    );
    assert_eq!(group_prefix("ephemeris_base.leb"), None);
    assert_eq!(group_prefix("custom.leb"), None);
    // Tier guard: agreement passes, mixture refuses, untiered is free.
    assert_eq!(
        check_single_tier(&["base_core.leb2", "base_apogee.leb2"]).unwrap(),
        Some(LebTier::Base)
    );
    assert!(matches!(
        check_single_tier(&["base_core.leb2", "medium_core.leb2"]).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    assert_eq!(
        check_single_tier(&["custom_a.leb", "custom_b.leb"]).unwrap(),
        None
    );
}

#[test]
fn composite_dispatch_is_first_wins_with_typed_misses() {
    let first = leb1_const_image(0, 5.0, 0.0, 2.0);
    let second = leb1_const_image(0, 7.0, 0.0, 4.0);
    let third = leb1_const_image(1, 9.0, 0.0, 4.0);
    let i0 = FileIndex::decode(&first).expect("first must decode");
    let i1 = FileIndex::decode(&second).expect("second must decode");
    let i2 = FileIndex::decode(&third).expect("third must decode");
    let composite = CompositeIndex::new(vec![
        (i0, first.as_slice()),
        (i1, second.as_slice()),
        (i2, third.as_slice()),
    ])
    .expect("composite must build");
    assert_eq!(composite.len(), 3);
    assert!(!composite.is_empty());
    // Body 0 comes from the first reader despite the overlap.
    assert_eq!(
        composite.eval_body_any(0, 1.0).unwrap(),
        (vec![5.0, 5.0, 5.0], vec![0.0, 0.0, 0.0])
    );
    // Body 1 only exists in the third reader.
    assert_eq!(
        composite.eval_body_any(1, 1.0).unwrap(),
        (vec![9.0, 9.0, 9.0], vec![0.0, 0.0, 0.0])
    );
    assert_eq!(composite.serving_index(0), Some(0));
    assert_eq!(composite.serving_index(1), Some(2));
    assert_eq!(composite.body_coverage(0), Some((0.0, 2.0)));
    assert_eq!(composite.jd_range(), (0.0, 4.0));
    assert!(composite.has_body(1));
    assert!(!composite.has_body(2));
    // Missing bodies — including missing support channels — name the id.
    assert_eq!(
        composite.eval_body_any(2, 1.0).unwrap_err(),
        Error::LebBodyNotFound { body_id: 2 }
    );
    // The serving entry exposes the channel for pipeline dispatch.
    assert_eq!(
        composite.serving_entry(0).unwrap().coord_type(),
        COORD_ICRS_BARY
    );
    assert!(composite.serving_entry(2).is_none());
    // Two-part path dispatches identically.
    let epoch = TwoPartEpoch::split(1.0).unwrap();
    assert_eq!(
        composite.eval_body_any_two_part(0, epoch).unwrap(),
        (vec![5.0, 5.0, 5.0], vec![0.0, 0.0, 0.0])
    );
    assert!(matches!(
        composite.eval_body_any_two_part(2, epoch).unwrap_err(),
        Error::LebBodyNotFound { .. }
    ));
    // Empty construction is an error, like the release ValueError.
    assert!(matches!(
        CompositeIndex::new(vec![]),
        Err(Error::LebCorrupt { .. })
    ));
}

#[test]
fn composite_routes_mixed_format_constituents() {
    let leb1 = minimal_leb1();
    let (leb2_bytes, _) = leb2_v1_image();
    let i_leb1 = FileIndex::decode(&leb1).expect("LEB1 must decode");
    let i_leb2 = FileIndex::decode(&leb2_bytes).expect("LEB2 must decode");
    // Body 0 exists in both: the LEB1 reader wins and serves the raw path.
    let composite = CompositeIndex::new(vec![
        (i_leb1, leb1.as_slice()),
        (i_leb2, leb2_bytes.as_slice()),
    ])
    .expect("mixed composite must build");
    assert_eq!(
        composite.eval_body_any(0, 0.25).unwrap(),
        FileIndex::decode(&leb1)
            .unwrap()
            .eval_body(&leb1, 0, 0.25)
            .unwrap()
    );
    // Reversed order: the LEB2 blob path serves body 0.
    let i_leb1 = FileIndex::decode(&leb1).expect("LEB1 must decode");
    let i_leb2 = FileIndex::decode(&leb2_bytes).expect("LEB2 must decode");
    let reversed = CompositeIndex::new(vec![
        (i_leb2, leb2_bytes.as_slice()),
        (i_leb1, leb1.as_slice()),
    ])
    .expect("reversed composite must build");
    assert_eq!(
        reversed.eval_body_any(0, 0.25).unwrap(),
        FileIndex::decode(&leb2_bytes)
            .unwrap()
            .eval_leb2_body(&leb2_bytes, 0, 0.25)
            .unwrap()
    );
}

#[test]
fn composite_aux_channels_come_from_first_carrier() {
    let bare = aux_image(vec![]);
    let nut = aux_image(vec![(SECTION_NUTATION, nutation_payload(1, 16.0))]);
    let i_bare = FileIndex::decode(&bare).expect("bare must decode");
    let i_nut = FileIndex::decode(&nut).expect("nutation image must decode");
    let composite = CompositeIndex::new(vec![(i_bare, bare.as_slice()), (i_nut, nut.as_slice())])
        .expect("composite must build");
    assert!(composite.has_nutation());
    // Mid-segment: tau 0 recovers the constant terms exactly.
    assert_eq!(composite.eval_nutation(8.0).unwrap(), (1.0, 2.0));
    // A `--skip-aux` zero-segment header does not shadow a later file.
    let skip = aux_image(vec![(SECTION_NUTATION, nutation_payload(0, 16.0))]);
    let i_skip = FileIndex::decode(&skip).expect("skip-aux must decode");
    let i_nut = FileIndex::decode(&nut).expect("nutation image must decode");
    let shadowed = CompositeIndex::new(vec![(i_skip, skip.as_slice()), (i_nut, nut.as_slice())])
        .expect("composite must build");
    assert_eq!(shadowed.eval_nutation(8.0).unwrap(), (1.0, 2.0));
    // No carrier at all: typed absence, never a panic or fallback.
    let i_bare = FileIndex::decode(&bare).expect("bare must decode");
    let empty =
        CompositeIndex::new(vec![(i_bare, bare.as_slice())]).expect("single composite must build");
    assert!(!empty.has_nutation());
    assert!(matches!(
        empty.eval_nutation(8.0).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    assert!(matches!(
        empty.delta_t(8.0).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    assert_eq!(
        empty.get_star(1).unwrap_err(),
        Error::LebStarNotFound { star_id: 1 }
    );
}

#[test]
fn composite_delta_t_and_stars_dispatch() {
    let dt = aux_image(vec![(
        SECTION_DELTA_T,
        delta_t_payload(&[(0.0, 1.0), (30.0, 2.0)]),
    )]);
    let stars = aux_image(vec![(SECTION_STARS, star_payload(&[7]))]);
    let i_dt = FileIndex::decode(&dt).expect("delta-t image must decode");
    let i_stars = FileIndex::decode(&stars).expect("star image must decode");
    let composite = CompositeIndex::new(vec![(i_dt, dt.as_slice()), (i_stars, stars.as_slice())])
        .expect("composite must build");
    // Linear interpolation between the two pairs.
    assert_eq!(composite.delta_t(15.0).unwrap(), 1.5);
    assert_eq!(composite.get_star(7).unwrap().star_id(), 7);
    assert_eq!(
        composite.get_star(8).unwrap_err(),
        Error::LebStarNotFound { star_id: 8 }
    );
}

#[test]
fn tiered_selection_prefers_highest_covering_tier() {
    let base_bytes = leb1_const_image(0, 5.0, 0.0, 2.0);
    let medium_bytes = leb1_const_image(0, 7.0, 1.0, 5.0);
    let base = CompositeIndex::new(vec![(
        FileIndex::decode(&base_bytes).expect("base must decode"),
        base_bytes.as_slice(),
    )])
    .expect("base composite must build");
    let medium = CompositeIndex::new(vec![(
        FileIndex::decode(&medium_bytes).expect("medium must decode"),
        medium_bytes.as_slice(),
    )])
    .expect("medium composite must build");
    // Installed out of order: priority still governs.
    let tiered = TieredIndex::new(vec![(LebTier::Medium, medium), (LebTier::Base, base)])
        .expect("tiered must build");
    assert_eq!(tiered.tier_at(0), Some(LebTier::Base));
    assert!(tiered.has_body(0));
    assert!(!tiered.has_body(3));
    // Overlap: base wins. Medium-only date: medium serves.
    assert_eq!(tiered.selected_tier(0, 1.5), Some(LebTier::Base));
    assert_eq!(tiered.selected_tier(0, 3.0), Some(LebTier::Medium));
    assert_eq!(tiered.selected_tier(0, 9.0), None);
    assert_eq!(
        tiered.eval_body_any(0, 1.5).unwrap(),
        (vec![5.0, 5.0, 5.0], vec![0.0, 0.0, 0.0])
    );
    assert_eq!(
        tiered.eval_body_any(0, 3.0).unwrap(),
        (vec![7.0, 7.0, 7.0], vec![0.0, 0.0, 0.0])
    );
    // Union coverage across tiers.
    assert_eq!(tiered.body_coverage(0), Some((0.0, 5.0)));
    assert_eq!(tiered.body_coverage(3), None);
    assert_eq!(tiered.jd_range(), (0.0, 5.0));
    // No tier covers: the highest-priority candidate evaluates so its
    // normal range error surfaces for the sealed dispatcher.
    assert_eq!(
        tiered.eval_body_any(0, 9.0).unwrap_err(),
        Error::LebOutOfRange {
            body_id: 0,
            jd: 9.0,
            start: 0.0,
            end: 2.0,
        }
    );
    // Unserved bodies name the id.
    assert_eq!(
        tiered.eval_body_any(3, 1.0).unwrap_err(),
        Error::LebBodyNotFound { body_id: 3 }
    );
    // Channel metadata comes from the highest-priority candidate.
    assert_eq!(
        tiered.serving_entry(0).unwrap().coord_type(),
        COORD_ICRS_BARY
    );
    // Construction guards: empty and duplicated tiers.
    let base_bytes = leb1_const_image(0, 5.0, 0.0, 2.0);
    let base = CompositeIndex::new(vec![(
        FileIndex::decode(&base_bytes).expect("base must decode"),
        base_bytes.as_slice(),
    )])
    .expect("base composite must build");
    assert!(matches!(
        TieredIndex::new(vec![]),
        Err(Error::LebCorrupt { .. })
    ));
    let base_bytes2 = leb1_const_image(0, 5.0, 0.0, 2.0);
    let base2 = CompositeIndex::new(vec![(
        FileIndex::decode(&base_bytes2).expect("base must decode"),
        base_bytes2.as_slice(),
    )])
    .expect("base composite must build");
    assert!(matches!(
        TieredIndex::new(vec![(LebTier::Base, base), (LebTier::Base, base2)]),
        Err(Error::LebCorrupt { .. })
    ));
}

#[test]
fn tiered_aux_dispatch_is_date_selected() {
    let nut = aux_image(vec![(SECTION_NUTATION, nutation_payload(1, 16.0))]);
    let bare = aux_image(vec![]);
    let nut_composite = CompositeIndex::new(vec![(
        FileIndex::decode(&nut).expect("nut must decode"),
        nut.as_slice(),
    )])
    .expect("nut composite must build");
    let bare_composite = CompositeIndex::new(vec![(
        FileIndex::decode(&bare).expect("bare must decode"),
        bare.as_slice(),
    )])
    .expect("bare composite must build");
    let tiered = TieredIndex::new(vec![
        (LebTier::Base, bare_composite),
        (LebTier::Medium, nut_composite),
    ])
    .expect("tiered must build");
    // Base range [0, 16] carries no nutation; the medium tier (same
    // file range here) serves it. A date outside every range fails.
    assert!(tiered.has_nutation());
    assert_eq!(tiered.eval_nutation(8.0).unwrap(), (1.0, 2.0));
    assert!(matches!(
        tiered.eval_nutation(100.0).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    // No nutation anywhere: typed absence.
    let bare2 = aux_image(vec![]);
    let bare_composite2 = CompositeIndex::new(vec![(
        FileIndex::decode(&bare2).expect("bare must decode"),
        bare2.as_slice(),
    )])
    .expect("bare composite must build");
    let bare_tiered =
        TieredIndex::new(vec![(LebTier::Base, bare_composite2)]).expect("tiered must build");
    assert!(!bare_tiered.has_nutation());
    assert!(matches!(
        bare_tiered.eval_nutation(8.0).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    assert!(matches!(
        bare_tiered.delta_t(8.0).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    assert_eq!(
        bare_tiered.get_star(1).unwrap_err(),
        Error::LebStarNotFound { star_id: 1 }
    );
}

#[test]
fn barycenter_channels_and_names_match_release() {
    assert_eq!(SYSTEM_BARY_BODY_IDS, [5, 6, 7, 8, 9]);
    assert!(is_system_barycenter_channel(COORD_ICRS_BARY_SYSTEM));
    assert!(!is_system_barycenter_channel(COORD_ICRS_BARY));
    assert!(!is_system_barycenter_channel(COORD_ECLIPTIC));
    assert_eq!(system_barycenter_name(5), Some("jupiter barycenter"));
    assert_eq!(system_barycenter_name(9), Some("pluto barycenter"));
    assert_eq!(system_barycenter_name(0), None);
    assert_eq!(CENTER_OFFSET_VELOCITY_STEP_DAYS, 0.0005);
}

#[test]
fn body_center_composition_applies_offsets_with_release_order() {
    let stored_pos = [1.0, 2.0, 4.0];
    let stored_vel = [0.25, -0.5, 1.0];
    // Non-system bodies pass through even with an available offset.
    let (pos, vel) = eval_body_center(stored_pos, stored_vel, false, 100.0, |_| {
        Some([10.0, 10.0, 10.0])
    });
    assert_eq!((pos, vel), (stored_pos, stored_vel));
    // No covering offset: the stored barycenter is retained explicitly.
    let (pos, vel) = eval_body_center(stored_pos, stored_vel, true, 100.0, |_| None);
    assert_eq!((pos, vel), (stored_pos, stored_vel));
    // Exact zero offset takes the release `center_pos == pos` fast path:
    // stored velocity returned bitwise unchanged.
    let (pos, vel) = eval_body_center(stored_pos, stored_vel, true, 100.0, |_| {
        Some([0.0, 0.0, 0.0])
    });
    assert_eq!((pos, vel), (stored_pos, stored_vel));
    // Linear offset rate [1, 2, -1] AU/day: position is the exact sum
    // (all operands exactly representable), velocity gains the
    // centered-difference rate. The numerator rounds each ~100-magnitude
    // side to at most half a ulp (~1.2e-14), and the 1/(2h) = 1000
    // scale lifts that to ~1.2e-11; the 1e-10 bound is derived, not fitted.
    let (pos, vel) = eval_body_center(stored_pos, stored_vel, true, 100.0, |jd| {
        Some([jd, 2.0 * jd, -jd])
    });
    assert_eq!(pos, [101.0, 202.0, -96.0]);
    for (got, want) in vel.iter().zip([1.25, 1.5, 0.0]) {
        assert!(
            (got - want).abs() <= 1e-10,
            "center velocity {got} too far from {want}"
        );
    }
    // Non-finite epochs skip the correction instead of propagating NaN.
    let (pos, vel) = eval_body_center(stored_pos, stored_vel, true, f64::NAN, |_| {
        Some([1.0, 1.0, 1.0])
    });
    assert_eq!((pos, vel), (stored_pos, stored_vel));
}

#[test]
fn batch_eval_matches_repeated_scalars_across_seam() {
    // LEB-12: the multi-epoch form shares the scalar path by
    // construction. Bounded grid over the F1 two-segment image:
    // both endpoints, the interior seam 1.0, adjacent float words
    // on both sides, and mid-segment points.
    let file = minimal_leb1();
    let index = FileIndex::decode(&file).expect("fixture must decode");
    let jds = [
        0.0,
        0.5,
        f64::next_down(1.0),
        1.0,
        f64::next_up(1.0),
        1.5,
        2.0,
    ];
    let batched = index
        .eval_body_many(&file, 0, &jds)
        .expect("in-range batch");
    assert_eq!(batched.len(), jds.len());
    for (jd, state) in jds.iter().zip(batched.iter()) {
        let scalar = index.eval_body(&file, 0, *jd).expect("scalar");
        assert_eq!(*state, scalar, "batch diverged from scalar at jd {jd}");
    }
    // Two-part batch agrees element-wise with repeated scalar
    // two-part calls on the same collapsed words.
    let epochs: Vec<TwoPartEpoch> = jds
        .iter()
        .map(|jd| TwoPartEpoch::split(*jd).expect("finite splits"))
        .collect();
    let batched_two = index
        .eval_body_two_part_many(&file, 0, &epochs)
        .expect("two-part batch");
    for (epoch, state) in epochs.iter().zip(batched_two.iter()) {
        let scalar = index.eval_body_two_part(&file, 0, *epoch).expect("scalar");
        assert_eq!(*state, scalar);
    }
}

#[test]
fn batch_eval_is_fail_fast_with_scalar_error_identity() {
    let file = minimal_leb1();
    let index = FileIndex::decode(&file).expect("fixture must decode");
    // The first out-of-range element aborts with the scalar error —
    // the release `_at` loop raises rather than skipping.
    let jds = [0.5, 1.5, 99.0, 0.25];
    assert_eq!(
        index.eval_body_many(&file, 0, &jds).unwrap_err(),
        index.eval_body(&file, 0, 99.0).unwrap_err()
    );
    // Non-finite element: the same typed refusal as the scalar path.
    assert!(matches!(
        index
            .eval_body_many(&file, 0, &[0.5, f64::NAN])
            .unwrap_err(),
        Error::NonFiniteArgument { .. }
    ));
    // Unserved body and empty input.
    assert_eq!(
        index.eval_body_many(&file, 7, &[0.5]).unwrap_err(),
        Error::LebBodyNotFound { body_id: 7 }
    );
    assert!(
        index
            .eval_body_many(&file, 0, &[])
            .expect("empty batch")
            .is_empty()
    );
}

#[test]
fn batch_eval_shares_scalar_semantics_over_composite_and_tiers() {
    // LEB-12 over dispatch: composite first-wins and per-epoch tier
    // selection still apply element-wise. Base covers [0, 1] with a
    // constant-5 channel, medium covers [1, 2] with constant-7.
    let base_bytes = leb1_const_image(0, 5.0, 0.0, 1.0);
    let med_bytes = leb1_const_image(0, 7.0, 1.0, 2.0);
    let base_index = FileIndex::decode(&base_bytes).expect("base decodes");
    let med_index = FileIndex::decode(&med_bytes).expect("medium decodes");
    let base = CompositeIndex::new(vec![(base_index, &base_bytes[..])]).expect("base composite");
    let med = CompositeIndex::new(vec![(med_index, &med_bytes[..])]).expect("medium composite");
    let tiered =
        TieredIndex::new(vec![(LebTier::Base, base), (LebTier::Medium, med)]).expect("tiers build");
    // Batch straddles the tier seam at 1.0 (owned by both ranges;
    // base wins there by priority).
    let jds = [0.25, f64::next_down(1.0), 1.0, 1.75];
    let batched = tiered.eval_body_any_many(0, &jds).expect("tiered batch");
    for (jd, state) in jds.iter().zip(batched.iter()) {
        let scalar = tiered.eval_body_any(0, *jd).expect("scalar");
        assert_eq!(*state, scalar, "tiered batch diverged at jd {jd}");
    }
    assert_eq!(batched[0].0[0], 5.0);
    assert_eq!(batched[2].0[0], 5.0);
    assert_eq!(batched[3].0[0], 7.0);
    // Composite batch over the merged image likewise loops scalars.
    let both = CompositeIndex::new(vec![
        (
            FileIndex::decode(&base_bytes).expect("base decodes"),
            &base_bytes[..],
        ),
        (
            FileIndex::decode(&med_bytes).expect("medium decodes"),
            &med_bytes[..],
        ),
    ])
    .expect("composite builds");
    let composite_batched = both
        .eval_body_any_many(0, &[0.25])
        .expect("composite batch");
    assert_eq!(
        composite_batched[0],
        both.eval_body_any(0, 0.25).expect("scalar")
    );
}

#[test]
fn coverage_tier_eligibility_is_cumulative() {
    // Release `_discover_reviewed_leb_tier_cores`: base admits only
    // base, medium admits base+medium, extended admits all three.
    assert_eq!(eligible_tiers(LebTier::Base), &[LebTier::Base]);
    assert_eq!(
        eligible_tiers(LebTier::Medium),
        &[LebTier::Base, LebTier::Medium]
    );
    assert_eq!(
        eligible_tiers(LebTier::Extended),
        &[LebTier::Base, LebTier::Medium, LebTier::Extended]
    );
    check_tier_subset(&[], LebTier::Base).expect("empty passes");
    check_tier_subset(&[LebTier::Base], LebTier::Base).expect("base under base");
    check_tier_subset(&[LebTier::Base, LebTier::Medium], LebTier::Medium)
        .expect("medium admits base");
    assert!(matches!(
        check_tier_subset(&[LebTier::Medium], LebTier::Base),
        Err(Error::LebCorrupt { .. })
    ));
    assert!(matches!(
        check_tier_subset(&[LebTier::Base, LebTier::Extended], LebTier::Medium),
        Err(Error::LebCorrupt { .. })
    ));
    // Index-level gate over the same rule.
    let base_bytes = leb1_const_image(0, 5.0, 0.0, 2.0);
    let base = CompositeIndex::new(vec![(
        FileIndex::decode(&base_bytes).expect("base must decode"),
        base_bytes.as_slice(),
    )])
    .expect("base composite must build");
    let tiered = TieredIndex::new(vec![(LebTier::Base, base)]).expect("tiered builds");
    tiered.check_eligible(LebTier::Base).expect("base eligible");
    tiered
        .check_eligible(LebTier::Medium)
        .expect("medium admits base");
    tiered
        .check_eligible(LebTier::Extended)
        .expect("wider active admits");
}

#[test]
fn fitted_model_windows_derive_from_generated_grid_metadata() {
    // Generated `lunar_apse_model` pins: START + STEP * (COUNT - 1).
    // Apogee: 2286820.5 + 10 * 40249 = 2689310.5.
    // Perigee: 2286820.5 + 2 * 201248 = 2689316.5.
    assert_eq!(APOGEE_MODEL_JD_END, 2689310.5);
    assert_eq!(PERIGEE_MODEL_JD_END, 2689316.5);
    assert_eq!(
        interpolated_apse_model_window(21),
        Some((2286820.5, 2689310.5))
    );
    assert_eq!(
        interpolated_apse_model_window(22),
        Some((2286820.5, 2689316.5))
    );
    assert_eq!(interpolated_apse_model_window(12), None);
    assert!(has_fitted_model_window(21));
    assert!(has_fitted_model_window(22));
    assert!(!has_fitted_model_window(12));
    // Body-id classifiers mirror the release planet sets.
    assert!(is_core_body(0) && is_core_body(12) && is_core_body(14));
    assert!(is_core_body(21) && is_core_body(22));
    assert!(!is_core_body(13) && !is_core_body(40) && !is_core_body(48));
    assert!(is_stored_state_body(0) && is_stored_state_body(9) && is_stored_state_body(14));
    assert!(!is_stored_state_body(10) && !is_stored_state_body(21));
    assert!(is_fictitious_body(40) && is_fictitious_body(58));
    assert!(!is_fictitious_body(39) && !is_fictitious_body(59));
    // Inclusive bounds; non-finite epochs never satisfy a bound.
    assert!(interval_contains(0.0, 2.0, 0.0));
    assert!(interval_contains(0.0, 2.0, 2.0));
    assert!(!interval_contains(0.0, 2.0, f64::NAN));
}

#[test]
fn tiered_coverage_for_reports_selected_or_union() {
    let base_bytes = leb1_const_image(0, 5.0, 0.0, 2.0);
    let medium_bytes = leb1_const_image(0, 7.0, 1.0, 5.0);
    let base = CompositeIndex::new(vec![(
        FileIndex::decode(&base_bytes).expect("base must decode"),
        base_bytes.as_slice(),
    )])
    .expect("base composite must build");
    let medium = CompositeIndex::new(vec![(
        FileIndex::decode(&medium_bytes).expect("medium must decode"),
        medium_bytes.as_slice(),
    )])
    .expect("medium composite must build");
    let tiered = TieredIndex::new(vec![(LebTier::Base, base), (LebTier::Medium, medium)])
        .expect("tiered must build");
    // Covered date: the selected tier's own interval.
    assert_eq!(
        tiered.coverage_for(0, Some(3.0)),
        Some(CoverageReport {
            start: 1.0,
            end: 5.0,
            selected: Some(1),
        })
    );
    assert_eq!(
        tiered.coverage_for(0, Some(0.5)),
        Some(CoverageReport {
            start: 0.0,
            end: 2.0,
            selected: Some(0),
        })
    );
    // Outside every stored interval: union, selected unset.
    assert_eq!(
        tiered.coverage_for(0, Some(9.0)),
        Some(CoverageReport {
            start: 0.0,
            end: 5.0,
            selected: None,
        })
    );
    // Date-blind: union as well.
    assert_eq!(
        tiered.coverage_for(0, None),
        Some(CoverageReport {
            start: 0.0,
            end: 5.0,
            selected: None,
        })
    );
    assert_eq!(tiered.coverage_for(3, Some(1.0)), None);
    // Fictitious ids never report coverage, even date-blind.
    assert_eq!(tiered.coverage_for(40, None), None);
    assert_eq!(tiered.coverage_for(40, Some(1.0)), None);
}

#[test]
fn sealed_outcome_matrix_matches_release_escalation() {
    // Core body in range: covered.
    assert_eq!(
        sealed_outcome(3, 1.0, Some((0.0, 2.0)), Some((0.0, 2.0))),
        SealedOutcome::Covered
    );
    // Core body out of range: refusal naming the requested body.
    assert_eq!(
        sealed_outcome(3, 9.0, Some((0.0, 2.0)), Some((0.0, 5.0))),
        SealedOutcome::RefusedRange {
            body_id: 3,
            jd: 9.0,
            start: 0.0,
            end: 2.0,
        }
    );
    // Stored-state core absent: missing-body, distinct from a range.
    assert_eq!(
        sealed_outcome(3, 1.0, None, Some((0.0, 5.0))),
        SealedOutcome::MissingCoreChannel { body_id: 3 }
    );
    // Derived core point absent (e.g. node without a stored channel):
    // no LEB miss to report; downstream model routing decides.
    assert_eq!(
        sealed_outcome(11, 1.0, None, Some((0.0, 5.0))),
        SealedOutcome::NoLebChannel { body_id: 11 }
    );
    // Curated companion out of its window while the Sun covers the
    // date: declared-local-model fall-through may proceed.
    assert_eq!(
        sealed_outcome(30, 9.0, Some((0.0, 2.0)), Some((0.0, 20.0))),
        SealedOutcome::LocalModelFallthrough
    );
    // Companion with Sun support also out: refusal names the
    // requested body with its own window.
    assert_eq!(
        sealed_outcome(30, 9.0, Some((0.0, 2.0)), Some((0.0, 5.0))),
        SealedOutcome::RefusedRange {
            body_id: 30,
            jd: 9.0,
            start: 0.0,
            end: 2.0,
        }
    );
    // Analytically modelled id with no channel of its own while
    // support is also out: refusal names the requested body with
    // the support window.
    assert_eq!(
        sealed_outcome(48, 9.0, None, Some((0.0, 5.0))),
        SealedOutcome::RefusedRange {
            body_id: 48,
            jd: 9.0,
            start: 0.0,
            end: 5.0,
        }
    );
    // Same id with support in range: not a range miss.
    assert_eq!(
        sealed_outcome(48, 3.0, None, Some((0.0, 5.0))),
        SealedOutcome::LocalModelFallthrough
    );
    // Fictitious ids carry no LEB channel, so the Sun-support rule
    // decides: support in range is fall-through, support out is a
    // refusal naming the requested body with the support window.
    assert_eq!(
        sealed_outcome(40, 3.0, None, Some((0.0, 5.0))),
        SealedOutcome::LocalModelFallthrough
    );
    assert_eq!(
        sealed_outcome(40, 9.0, None, Some((0.0, 5.0))),
        SealedOutcome::RefusedRange {
            body_id: 40,
            jd: 9.0,
            start: 0.0,
            end: 5.0,
        }
    );
    assert_eq!(
        sealed_outcome(40, 9.0, None, None),
        SealedOutcome::NoLebChannel { body_id: 40 }
    );
    // Fitted model window binds INTP_APOG even with stored coverage
    // claiming more: outside the grid is a refusal with the grid.
    assert_eq!(
        sealed_outcome(
            21,
            1000000.0,
            Some((0.0, 5000000.0)),
            Some((0.0, 5000000.0))
        ),
        SealedOutcome::RefusedRange {
            body_id: 21,
            jd: 1000000.0,
            start: 2286820.5,
            end: 2689310.5,
        }
    );
    // Refusals surface typed errors; routing dispositions carry none.
    assert_eq!(
        SealedOutcome::RefusedRange {
            body_id: 3,
            jd: 9.0,
            start: 0.0,
            end: 2.0,
        }
        .refusal_error(),
        Some(Error::LebOutOfRange {
            body_id: 3,
            jd: 9.0,
            start: 0.0,
            end: 2.0,
        })
    );
    assert_eq!(
        SealedOutcome::MissingCoreChannel { body_id: 3 }.refusal_error(),
        Some(Error::LebBodyNotFound { body_id: 3 })
    );
    assert_eq!(SealedOutcome::LocalModelFallthrough.refusal_error(), None);
    assert_eq!(
        SealedOutcome::NoLebChannel { body_id: 40 }.refusal_error(),
        None
    );
}

#[test]
fn sealed_precheck_agrees_with_index_unions() {
    let sun_bytes = leb1_const_image(0, 1.0, 0.0, 5.0);
    let companion_bytes = leb1_const_image(30, 2.0, 0.0, 2.0);
    // One tier holding both channels: unions are the stored windows.
    let both = CompositeIndex::new(vec![
        (
            FileIndex::decode(&sun_bytes).expect("sun must decode"),
            sun_bytes.as_slice(),
        ),
        (
            FileIndex::decode(&companion_bytes).expect("companion must decode"),
            companion_bytes.as_slice(),
        ),
    ])
    .expect("both composite must build");
    let tiered = TieredIndex::new(vec![(LebTier::Base, both)]).expect("tiered builds");
    assert_eq!(tiered.sealed_precheck(0, 4.0), SealedOutcome::Covered);
    assert_eq!(
        tiered.sealed_precheck(0, 9.0),
        SealedOutcome::RefusedRange {
            body_id: 0,
            jd: 9.0,
            start: 0.0,
            end: 5.0,
        }
    );
    // Companion outside its window with Sun support in range:
    // fall-through; with support out too: refusal.
    assert_eq!(
        tiered.sealed_precheck(30, 4.0),
        SealedOutcome::LocalModelFallthrough
    );
    assert_eq!(
        tiered.sealed_precheck(30, 9.0),
        SealedOutcome::RefusedRange {
            body_id: 30,
            jd: 9.0,
            start: 0.0,
            end: 2.0,
        }
    );
}

#[test]
fn batch_eval_shares_scalar_semantics_over_leb2_blobs() {
    // The LEB2 many-paths decode per element through the same chunk/
    // monolith pipeline as the scalar calls (bounded: one v1 and one
    // v2 batch, seam included).
    let v1 = leb2_v1_distinct_image();
    let v1_index = FileIndex::decode(&v1).expect("v1 must decode");
    let jds = [0.25, 1.0, 1.5, 2.0];
    let batched = v1_index
        .eval_leb2_body_many(&v1, 0, &jds)
        .expect("v1 batch");
    for (jd, state) in jds.iter().zip(batched.iter()) {
        assert_eq!(
            *state,
            v1_index.eval_leb2_body(&v1, 0, *jd).expect("scalar")
        );
    }
    let v2 = leb2_v2_distinct_image(1.0);
    let v2_index = FileIndex::decode(&v2).expect("v2 must decode");
    let batched = v2_index
        .eval_leb2_body_many(&v2, 0, &jds)
        .expect("v2 batch");
    for (jd, state) in jds.iter().zip(batched.iter()) {
        assert_eq!(
            *state,
            v2_index.eval_leb2_body(&v2, 0, *jd).expect("scalar")
        );
    }
}

#[test]
fn inventory_group_and_precision_rules_match_release() {
    // `_group_from_path`: stem rule, not the known-suffix table.
    assert_eq!(
        group_from_path(Some("/data/leb/base_core.leb2")),
        Some("core".to_string())
    );
    assert_eq!(
        group_from_path(Some("medium_apogee.leb2")),
        Some("apogee".to_string())
    );
    assert_eq!(
        group_from_path(Some("custom_myset.leb")),
        Some("myset".to_string())
    );
    assert_eq!(group_from_path(Some("nounderscore")), None);
    assert_eq!(group_from_path(Some("plain.leb")), None);
    assert_eq!(group_from_path(None), None);
    // `_leb_precision_class` branches over hand-picked witnesses.
    assert_eq!(
        precision_class(3, Some("base_core.leb2"), true),
        PrecisionClass::Ephemeris
    );
    for analytical in [10, 11, 12, 13, 21, 22, 40, 58] {
        assert_eq!(
            precision_class(analytical, Some("base_core.leb2"), true),
            PrecisionClass::Analytical,
            "body {analytical}"
        );
    }
    assert_eq!(
        precision_class(3, Some("base_core.leb2"), false),
        PrecisionClass::UnverifiedLocal
    );
    // N-body exotic 16 (Pholus): union, base-tier, extended-tier.
    assert_eq!(precision_class(16, None, true), PrecisionClass::Mixed);
    assert_eq!(
        precision_class(16, Some("base_exotics.leb2"), true),
        PrecisionClass::Ephemeris
    );
    assert_eq!(
        precision_class(16, Some("extended_exotics.leb2"), true),
        PrecisionClass::NumericalModel
    );
    // A perturber-set exotic (10010 Hygiea) is a plain ephemeris id.
    assert_eq!(
        precision_class(10010, Some("extended_asteroids.leb2"), true),
        PrecisionClass::Ephemeris
    );
    assert_eq!(PrecisionClass::Mixed.as_str(), "mixed");
}

#[test]
fn file_and_composite_coverage_records_label_the_serving_file() {
    let image = leb1_const_image(3, 1.0, 0.0, 5.0);
    let index = FileIndex::decode(&image).expect("must decode");
    // Single file: bounds, labels, closed-interval contains.
    let record = index
        .coverage_record(3, Some("/data/leb/base_core.leb2"), true)
        .expect("body 3 is indexed");
    assert_eq!(record.body_id(), 3);
    assert_eq!(record.source(), "LEB");
    assert_eq!(record.precision_class(), PrecisionClass::Ephemeris);
    assert_eq!((record.jd_start(), record.jd_end()), (0.0, 5.0));
    assert_eq!(record.data_file(), Some("/data/leb/base_core.leb2"));
    assert_eq!(record.group(), Some("core"));
    assert!(record.reviewed());
    assert!(record.contains(0.0) && record.contains(5.0));
    assert!(!record.contains(5.25));
    // Absent and fictitious bodies answer None, like the release.
    assert!(
        index
            .coverage_record(4, Some("base_core.leb2"), true)
            .is_none()
    );
    assert!(
        index
            .coverage_record(48, Some("base_core.leb2"), true)
            .is_none()
    );
    // Composite: first-wins serving file labels the record.
    let other = leb1_const_image(4, 2.0, 10.0, 12.0);
    let other_index = FileIndex::decode(&other).expect("must decode");
    let composite = CompositeIndex::new(vec![
        (index, image.as_slice()),
        (other_index, other.as_slice()),
    ])
    .expect("composite builds");
    let files = [
        ServingFile {
            path: Some("base_core.leb2"),
            reviewed: true,
        },
        ServingFile {
            path: Some("base_asteroids.leb2"),
            reviewed: true,
        },
    ];
    let record = composite
        .coverage_record(&files, 4)
        .expect("body 4 is served");
    assert_eq!((record.jd_start(), record.jd_end()), (10.0, 12.0));
    assert_eq!(record.data_file(), Some("base_asteroids.leb2"));
    assert_eq!(record.group(), Some("asteroids"));
    assert!(composite.coverage_record(&files, 48).is_none());
    // A short identity slice labels the record unlabeled, never fails.
    let record = composite
        .coverage_record(&files[..1], 4)
        .expect("short slice still resolves");
    assert_eq!(record.data_file(), None);
    assert_eq!(record.group(), None);
    assert!(!record.reviewed());
}

#[test]
fn tiered_coverage_reports_selected_interval_or_unlabeled_union() {
    let base_image = leb1_const_image(3, 1.0, 0.0, 5.0);
    let base_index = FileIndex::decode(&base_image).expect("must decode");
    let base = CompositeIndex::new(vec![(base_index, base_image.as_slice())]).expect("base builds");
    let medium_image = leb1_const_image(3, 1.0, 4.0, 9.0);
    let medium_index = FileIndex::decode(&medium_image).expect("must decode");
    let medium =
        CompositeIndex::new(vec![(medium_index, medium_image.as_slice())]).expect("medium builds");
    let tiered = TieredIndex::new(vec![(LebTier::Base, base), (LebTier::Medium, medium)])
        .expect("tiered builds");
    let file = |tier: LebTier, _pos: usize| {
        Some(ServingFile {
            path: match tier {
                LebTier::Base => Some("base_core.leb2"),
                LebTier::Medium => Some("medium_core.leb2"),
                LebTier::Extended => Some("extended_core.leb2"),
            },
            reviewed: true,
        })
    };
    // Covered dates report the selected tier's own interval + file.
    let record = tiered
        .coverage_record(3, Some(1.0), file)
        .expect("jd 1 is base-covered");
    assert_eq!((record.jd_start(), record.jd_end()), (0.0, 5.0));
    assert_eq!(record.data_file(), Some("base_core.leb2"));
    let record = tiered
        .coverage_record(3, Some(7.0), file)
        .expect("jd 7 is medium-covered");
    assert_eq!((record.jd_start(), record.jd_end()), (4.0, 9.0));
    assert_eq!(record.data_file(), Some("medium_core.leb2"));
    // Outside every interval — or dateless — the union stays unlabeled.
    let record = tiered
        .coverage_record(3, Some(20.0), file)
        .expect("union still reported");
    assert_eq!((record.jd_start(), record.jd_end()), (0.0, 9.0));
    assert_eq!(record.data_file(), None);
    assert!(record.reviewed());
    let record = tiered
        .coverage_record(3, None, file)
        .expect("dateless union reported");
    assert_eq!(record.data_file(), None);
    assert!(tiered.coverage_record(48, Some(1.0), file).is_none());
    assert!(tiered.coverage_record(5, Some(1.0), file).is_none());
}

#[test]
fn inspect_index_reports_identity_and_sorted_bodies() {
    let image = leb1_const_image(3, 1.0, 0.0, 5.0);
    let index = FileIndex::decode(&image).expect("must decode");
    let info = inspect_index("/data/leb/base_core.leb2", Some(4096), true, &index);
    assert_eq!(info.name(), "base_core.leb2");
    assert_eq!(info.path(), "/data/leb/base_core.leb2");
    assert_eq!(info.group(), Some("core"));
    assert_eq!(info.size_bytes(), Some(4096));
    assert!(info.reviewed());
    assert_eq!(info.body_count(), 1);
    assert_eq!(info.bodies()[0].body_id(), 3);
    assert_eq!(
        (info.bodies()[0].jd_start(), info.bodies()[0].jd_end()),
        (0.0, 5.0)
    );
    // Byte-level inspection derives the size from the image itself.
    let via_bytes = inspect_bytes("medium_apogee.leb2", &image, false).expect("must decode");
    assert_eq!(via_bytes.size_bytes(), Some(image.len() as u64));
    assert_eq!(via_bytes.group(), Some("apogee"));
    assert!(!via_bytes.reviewed());
    assert!(inspect_bytes("broken.leb", b"too short", false).is_err());
}

#[test]
fn runtime_requirements_expand_cumulative_tiers_over_manifest_pins() {
    let groups = ["core", "asteroids", "exotics", "apogee"];
    let pins = [
        ("base_core.leb2", "aa"),
        ("base_asteroids.leb2", "ab"),
        ("base_exotics.leb2", "ac"),
        ("base_apogee.leb2", "ad"),
        ("medium_core.leb2", "ba"),
        ("medium_asteroids.leb2", "bb"),
        ("medium_exotics.leb2", "bc"),
        ("medium_apogee.leb2", "bd"),
        ("extended_core.leb2", "ca"),
        ("extended_asteroids.leb2", "cb"),
        ("extended_exotics.leb2", "cc"),
        ("extended_apogee.leb2", "cd"),
    ];
    // Base tier: four files, tier-major order, joined data-dir paths.
    let base =
        runtime_data_requirements(LebTier::Base, &groups, &pins, "/data").expect("pins complete");
    assert_eq!(base.len(), 4);
    assert_eq!(base[0].name(), "base_core.leb2");
    assert_eq!(base[0].kind(), "leb2");
    assert_eq!(base[0].group(), "core");
    assert_eq!(base[0].path(), "/data/leb/base_core.leb2");
    assert_eq!(base[0].sha256(), "aa");
    assert_eq!(base[3].name(), "base_apogee.leb2");
    // Medium is cumulative (8), extended is cumulative (12).
    assert_eq!(
        runtime_data_requirements(LebTier::Medium, &groups, &pins, "/data")
            .expect("pins complete")
            .len(),
        8
    );
    let extended = runtime_data_requirements(LebTier::Extended, &groups, &pins, "/data")
        .expect("pins complete");
    assert_eq!(extended.len(), 12);
    assert_eq!(extended[11].name(), "extended_apogee.leb2");
    // A missing (or empty) pin is a packaging error, never partial output.
    let mut short = pins.to_vec();
    short.retain(|(name, _)| *name != "medium_exotics.leb2");
    assert!(matches!(
        runtime_data_requirements(LebTier::Medium, &groups, &short, "/data").unwrap_err(),
        Error::DataManifestMissing { name } if name == "medium_exotics.leb2"
    ));
    let mut emptied = pins;
    emptied[0] = ("base_core.leb2", "");
    assert!(matches!(
        runtime_data_requirements(LebTier::Base, &groups, &emptied, "/data").unwrap_err(),
        Error::DataManifestMissing { .. }
    ));
}

#[test]
fn tier_assets_report_ready_only_when_every_requirement_parses() {
    let groups = ["core", "asteroids"];
    let pins = [("base_core.leb2", "aa"), ("base_asteroids.leb2", "ab")];
    let core_image = leb1_const_image(3, 1.0, 0.0, 5.0);
    let asteroids_image = leb1_const_image(4, 2.0, 10.0, 12.0);
    // Direct rows: ready carries size + body count, absence is
    // missing, undecodable bytes are unreadable with a message.
    let requirements =
        runtime_data_requirements(LebTier::Base, &groups, &pins, "/data").expect("pins complete");
    let ready = check_asset(
        &requirements[0],
        Some(ObservedFile {
            name: "base_core.leb2",
            bytes: Some(&core_image),
            size_bytes: Some(core_image.len() as u64),
        }),
    );
    assert_eq!(ready.name(), "base_core.leb2");
    assert_eq!(ready.group(), "core");
    assert_eq!(ready.path(), "/data/leb/base_core.leb2");
    assert!(ready.is_ready());
    assert!(matches!(
        ready.state(),
        AssetState::Ready {
            size_bytes: _,
            body_count: 1
        }
    ));
    let missing = check_asset(&requirements[1], None);
    assert!(!missing.is_ready());
    assert_eq!(missing.state(), &AssetState::Missing);
    let mut broken = asteroids_image.clone();
    broken[0] = b'X';
    let unreadable = check_asset(
        &requirements[1],
        Some(ObservedFile {
            name: "base_asteroids.leb2",
            bytes: Some(&broken),
            size_bytes: None,
        }),
    );
    assert!(!unreadable.is_ready());
    assert!(matches!(
        unreadable.state(),
        AssetState::Unreadable { what } if !what.is_empty()
    ));
    // Full tier: one missing companion denies the pass and is named.
    let observed = [ObservedFile {
        name: "base_core.leb2",
        bytes: Some(&core_image),
        size_bytes: Some(core_image.len() as u64),
    }];
    let report = validate_tier_assets(LebTier::Base, &groups, &pins, "/data", &observed)
        .expect("pins complete");
    assert_eq!(report.tier(), LebTier::Base);
    assert!(!report.ready());
    assert_eq!(report.missing(), ["base_asteroids.leb2"]);
    assert!(report.unreadable().is_empty());
    // Retired companions left on disk never satisfy a requirement.
    let observed = [
        ObservedFile {
            name: "base_core.leb2",
            bytes: Some(&core_image),
            size_bytes: None,
        },
        ObservedFile {
            name: "base_asteroids.leb2",
            bytes: Some(&asteroids_image),
            size_bytes: None,
        },
        ObservedFile {
            name: "base_uranians.leb2",
            bytes: Some(&core_image),
            size_bytes: None,
        },
    ];
    let report = validate_tier_assets(LebTier::Base, &groups, &pins, "/data", &observed)
        .expect("pins complete");
    assert!(report.ready());
    assert_eq!(report.files().len(), 2);
    assert!(report.missing().is_empty());
    // A corrupt companion is reported unreadable, never ready.
    let observed = [
        ObservedFile {
            name: "base_core.leb2",
            bytes: Some(&core_image),
            size_bytes: None,
        },
        ObservedFile {
            name: "base_asteroids.leb2",
            bytes: Some(b"too short"),
            size_bytes: Some(9),
        },
    ];
    let report = validate_tier_assets(LebTier::Base, &groups, &pins, "/data", &observed)
        .expect("pins complete");
    assert!(!report.ready());
    assert!(report.missing().is_empty());
    assert_eq!(report.unreadable(), ["base_asteroids.leb2"]);
}

#[test]
fn pinned_sibling_rule_matches_release() {
    let manifest = [
        "base_core.leb2",
        "base_asteroids.leb2",
        "base_apogee.leb2",
        "base_exotics.leb2",
        "planet_centers_base.bsp",
    ];
    // Same-tier companions present: siblings available.
    assert!(has_pinned_siblings(
        "base_core.leb2",
        &["base_core.leb2", "base_asteroids.leb2"],
        &manifest,
    ));
    // The core itself never counts as its own sibling.
    assert!(!has_pinned_siblings(
        "base_core.leb2",
        &["base_core.leb2"],
        &manifest,
    ));
    // A companion on disk but outside the manifest does not qualify.
    assert!(!has_pinned_siblings(
        "base_core.leb2",
        &["base_core.leb2", "base_custom.leb2"],
        &manifest,
    ));
    // Non-LEB manifest entries (planet-center BSPs) do not qualify.
    assert!(!has_pinned_siblings(
        "planet_centers_base.bsp",
        &["planet_centers_base.bsp", "planet_centers_medium.bsp"],
        &manifest,
    ));
    // Stems without '_' and full paths follow the release rules.
    assert!(!has_pinned_siblings(
        "plain.leb2",
        &["plain.leb2"],
        &manifest
    ));
    assert!(has_pinned_siblings(
        "/data/leb/base_core.leb2",
        &["base_core.leb2", "base_apogee.leb2"],
        &manifest,
    ));
}

fn local_leb_dir() -> Option<std::path::PathBuf> {
    if let Ok(dir) = std::env::var("LIBEPHEMERIS_DATA_DIR") {
        let candidate = std::path::PathBuf::from(dir).join("leb");
        if candidate.is_dir() {
            return Some(candidate);
        }
    }
    if let Ok(home) = std::env::var("HOME") {
        let candidate = std::path::PathBuf::from(home)
            .join(".libephemeris")
            .join("leb");
        if candidate.is_dir() {
            return Some(candidate);
        }
    }
    None
}

#[test]
fn local_tier_assets_match_installed_distribution() {
    let Some(dir) = local_leb_dir() else {
        eprintln!("SKIP: no local LEB distribution installed");
        return;
    };
    let listing: Vec<String> = std::fs::read_dir(&dir)
        .expect("leb dir is readable")
        .filter_map(|entry| entry.ok())
        .map(|entry| entry.file_name().to_string_lossy().into_owned())
        .collect();
    let present: Vec<&str> = listing.iter().map(String::as_str).collect();
    // The 12 reviewed LEB2 names are the manifest under test.
    let mut manifest: Vec<String> = Vec::new();
    for tier in ["base", "medium", "extended"] {
        for group in LEB2_GROUPS {
            manifest.push(format!("{tier}_{group}.leb2"));
        }
    }
    let manifest_ref: Vec<&str> = manifest.iter().map(String::as_str).collect();
    // Base tier: every required file decodes with a non-empty index.
    let mut observed: Vec<ObservedFile> = Vec::new();
    let mut images: Vec<Vec<u8>> = Vec::new();
    for group in LEB2_GROUPS {
        let name = format!("base_{group}.leb2");
        let bytes = std::fs::read(dir.join(&name))
            .unwrap_or_else(|_| panic!("installed base file {name} is readable"));
        let index = FileIndex::decode(&bytes)
            .unwrap_or_else(|_| panic!("installed base file {name} parses"));
        assert!(
            !index.bodies().is_empty(),
            "installed base file {name} indexes bodies"
        );
        images.push(bytes);
    }
    let names: Vec<String> = LEB2_GROUPS
        .iter()
        .map(|group| format!("base_{group}.leb2"))
        .collect();
    for (name, bytes) in names.iter().zip(images.iter()) {
        observed.push(ObservedFile {
            name: name.as_str(),
            bytes: Some(bytes.as_slice()),
            size_bytes: Some(bytes.len() as u64),
        });
    }
    let groups: Vec<&str> = LEB2_GROUPS.to_vec();
    let pins: Vec<(&str, &str)> = manifest_ref
        .iter()
        .take(4)
        .map(|name| (*name, "pin-not-checked"))
        .collect();
    // Pins are non-empty so the requirement set builds; byte identity
    // against the manifest is install-time  work, not re-hashed
    // here over multi-GB assets.
    let report = validate_tier_assets(
        LebTier::Base,
        &groups,
        &pins,
        dir.to_string_lossy().as_ref(),
        &observed,
    )
    .expect("pins complete");
    assert!(report.ready(), "installed base tier validates ready");
    // Sibling companions are visible in the installed listing.
    assert!(has_pinned_siblings(
        "base_core.leb2",
        &present,
        &manifest_ref
    ));
    // Retired uranians companions are present on disk but required by
    // nothing and satisfy nothing.
    assert!(
        !manifest_ref.iter().any(|name| name.contains("uranians")),
        "retired companions stay outside the requirement set"
    );
    // Extended tier is absent here: the validator must report every
    // requirement missing, never a partial pass.
    if !present.contains(&"extended_core.leb2") {
        let empty: Vec<ObservedFile> = Vec::new();
        let extended_pins: Vec<(&str, &str)> = manifest_ref
            .iter()
            .map(|name| (*name, "pin-not-checked"))
            .collect();
        let report = validate_tier_assets(
            LebTier::Extended,
            &groups,
            &extended_pins,
            dir.to_string_lossy().as_ref(),
            &empty,
        )
        .expect("pins complete");
        assert!(!report.ready());
        assert_eq!(report.missing().len(), 12);
    }
}

// ------------------------------------------------------------------
// LEB-11: file/buffer ownership and warm/cool lifecycle.
// ------------------------------------------------------------------

/// Four-segment LEB1 builder: 1 body, degree 1, 3 components,
/// `seg_bytes = 2·3·8 = 48`, zero blob.
fn four_segment_leb1() -> (Vec<u8>, u64) {
    let mut file = pack_header(&MAGIC_LEB1, 1, 1, 1, 0.0, 4.0, 0.0, 0);
    let dir_at = file.len();
    let body_at = dir_at + SECTION_DIR_SIZE;
    let blob_at = body_at + LEB1_BODY_ENTRY_SIZE;
    file.extend_from_slice(&pack_section(0, body_at as u64, 52));
    file.extend_from_slice(&pack_body_leb1(
        0,
        0,
        4,
        0.0,
        4.0,
        1.0,
        1,
        3,
        blob_at as u64,
    ));
    file.extend_from_slice(&[0u8; 192]);
    (file, blob_at as u64)
}

#[test]
fn leb1_warm_ranges_follow_the_release_segment_rule() {
    let (file, blob_at) = four_segment_leb1();
    let index = FileIndex::decode(&file).expect("fixture must decode");
    // [0.5, 2.5]: first = int(0.5) = 0, last = int(2.5) = 2.
    assert_eq!(
        index.warm_ranges(&file, 0.5, 2.5).unwrap(),
        vec![(blob_at, 3 * 48)]
    );
    // Single-segment request pins one segment exactly.
    assert_eq!(
        index.warm_ranges(&file, 1.0, 1.0).unwrap(),
        vec![(blob_at + 48, 48)]
    );
    // Clamped tail: int(99.0) = 99 clamps to segment 3.
    assert_eq!(
        index.warm_ranges(&file, 3.0, 99.0).unwrap(),
        vec![(blob_at + 3 * 48, 48)]
    );
    // Disjoint range and reversed bounds advise nothing, like the
    // release loop skipping every body.
    assert!(index.warm_ranges(&file, 10.0, 12.0).unwrap().is_empty());
    assert!(index.warm_ranges(&file, 2.0, 1.0).unwrap().is_empty());
    // Advisory bounds must be finite: the release would silently
    // advise an unspecified subset instead.
    assert!(matches!(
        index.warm_ranges(&file, f64::NAN, 1.0).unwrap_err(),
        Error::NonFiniteArgument { .. }
    ));
    assert!(matches!(
        index.warm_ranges(&file, 0.0, f64::INFINITY).unwrap_err(),
        Error::NonFiniteArgument { .. }
    ));
}

#[test]
fn leb2_warm_ranges_report_blob_and_chunk_spans() {
    // v1: the whole body blob for any overlapping request.
    let (file, _) = leb2_v1_image();
    let index = FileIndex::decode(&file).expect("v1 image must decode");
    let entry = index.find_body(0).expect("body 0 indexed");
    assert_eq!(
        index.warm_ranges(&file, 0.0, 2.0).unwrap(),
        vec![(entry.data_offset(), entry.compressed_size().unwrap())]
    );
    assert!(index.warm_ranges(&file, 5.0, 6.0).unwrap().is_empty());
    // v2: one span per overlapping chunk, not the whole index.
    let file = leb2_v2_distinct_image(1.0);
    let index = FileIndex::decode(&file).expect("v2 image must decode");
    let (_, chunks) = decode_chunk_index(&file, index.find_body(0).unwrap().data_offset(), 0)
        .expect("chunk index decodes");
    assert_eq!(
        index.warm_ranges(&file, 0.0, 2.0).unwrap(),
        vec![
            (chunks[0].blob_offset(), chunks[0].compressed_size()),
            (chunks[1].blob_offset(), chunks[1].compressed_size()),
        ]
    );
    assert_eq!(
        index.warm_ranges(&file, 0.0, 0.5).unwrap(),
        vec![(chunks[0].blob_offset(), chunks[0].compressed_size())]
    );
}

#[test]
fn merge_ranges_pins_the_release_madvise_rule() {
    // Adjacent spans within a page merge; distant ones stay split.
    assert_eq!(
        merge_ranges(vec![(0, 100), (200, 100), (100_000, 50)], 1_000_000),
        vec![(0, 4096), (98_304, 4096)]
    );
    // Unsorted input merges identically.
    assert_eq!(
        merge_ranges(vec![(100_000, 50), (0, 100), (200, 100)], 1_000_000),
        vec![(0, 4096), (98_304, 4096)]
    );
    // Page alignment: start rounds down, end rounds up, end clamps
    // to the file length, empties drop.
    assert_eq!(merge_ranges(vec![(1, 1)], 100), vec![(0, 100)]);
    assert_eq!(merge_ranges(vec![(5000, 10)], 6000), vec![(4096, 1904)]);
    assert_eq!(merge_ranges(vec![(9000, 100)], 8000), Vec::new());
    assert!(merge_ranges(Vec::new(), 100).is_empty());
}

#[test]
fn leb_file_owns_eval_warm_and_close_lifecycle() {
    let (bytes, blob_at) = four_segment_leb1();
    let direct = FileIndex::decode(&bytes).expect("fixture must decode");
    let mut owned = LebFile::from_bytes("four.leb".to_string(), bytes.clone())
        .expect("owned image must decode");
    assert_eq!(owned.path(), "four.leb");
    assert!(!owned.is_closed());
    assert_eq!(owned.len(), bytes.len());
    assert_eq!(owned.format().unwrap(), LebFormat::Leb1);
    // Owned evaluation is bitwise-identical to the borrowed path.
    for jd in [0.0, 0.5, 1.0, 2.75, 4.0] {
        assert_eq!(
            owned.eval_body_any(0, jd).unwrap(),
            direct.eval_body(&bytes, 0, jd).unwrap(),
            "owned eval must equal borrowed eval at {jd}"
        );
    }
    let epoch = TwoPartEpoch::split(1.25).expect("epoch splits");
    assert_eq!(
        owned.eval_body_any_two_part(0, epoch).unwrap(),
        direct.eval_body_two_part(&bytes, 0, epoch).unwrap()
    );
    assert_eq!(
        owned.eval_body_any_many(0, &[0.5, 1.5]).unwrap(),
        vec![
            direct.eval_body(&bytes, 0, 0.5).unwrap(),
            direct.eval_body(&bytes, 0, 1.5).unwrap(),
        ]
    );
    // Warm advises the merged blob span: [0.5, 2.5] segments merge
    // with the page rule to one page-aligned range.
    let ranges = owned.warm_ranges(0.5, 2.5).unwrap();
    assert_eq!(
        ranges,
        merge_ranges(vec![(blob_at, 3 * 48)], bytes.len() as u64)
    );
    assert_eq!(owned.warm(0.5, 2.5).unwrap(), ranges.len());
    assert_eq!(owned.warm(10.0, 12.0).unwrap(), 0);
    owned.cool();
    // Partial reads serve exact slices; past-EOF is corruption.
    assert_eq!(
        owned.read_range(blob_at, 48).unwrap(),
        &bytes[blob_at as usize..][..48]
    );
    assert!(matches!(
        owned.read_range(bytes.len() as u64, 1).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    // Close drops the image: eval/warm fail closed, cool/close stay
    // safe and idempotent.
    owned.close();
    assert!(owned.is_closed());
    assert_eq!(owned.len(), 0);
    assert_eq!(
        owned.eval_body_any(0, 1.0).unwrap_err(),
        Error::LebReaderClosed
    );
    assert_eq!(owned.warm(0.0, 1.0).unwrap_err(), Error::LebReaderClosed);
    assert!(matches!(
        owned.read_range(0, 1).unwrap_err(),
        Error::LebReaderClosed
    ));
    owned.cool();
    owned.close();
    // Corrupt inputs never build a half-owned file: a truncated
    // image is corruption, and a full-size image with an unknown
    // magic keeps its identity variant.
    assert!(matches!(
        LebFile::from_bytes("short.leb".to_string(), vec![0u8; 10]).unwrap_err(),
        Error::LebCorrupt { .. }
    ));
    let mut bad_magic = vec![0u8; 64];
    bad_magic[0..4].copy_from_slice(b"BAD!");
    assert!(matches!(
        LebFile::from_bytes("bad.leb".to_string(), bad_magic).unwrap_err(),
        Error::LebUnknownFormat { .. }
    ));
}

#[test]
fn leb_file_open_reports_missing_files_and_detects_format() {
    let dir = std::env::temp_dir().join("libephemeris-leb11-probe");
    let missing = dir.join("no-such-file.leb");
    let _ = std::fs::remove_file(&missing);
    assert!(matches!(
        LebFile::open(&missing).unwrap_err(),
        Error::LebFileNotFound { .. }
    ));
    assert!(matches!(
        open_leb(&missing).unwrap_err(),
        Error::LebFileNotFound { .. }
    ));
    // Round-trip through a real file: format auto-detection serves
    // the LEB1 image with identical values.
    let _ = std::fs::create_dir_all(&dir);
    let path = dir.join("four.leb");
    let (bytes, _) = four_segment_leb1();
    std::fs::write(&path, &bytes).expect("scratch write");
    let owned = open_leb(&path).expect("scratch file opens");
    assert_eq!(owned.format().unwrap(), LebFormat::Leb1);
    let direct = FileIndex::decode(&bytes).expect("fixture must decode");
    assert_eq!(
        owned.eval_body_any(0, 1.5).unwrap(),
        direct.eval_body(&bytes, 0, 1.5).unwrap()
    );
    let _ = std::fs::remove_file(&path);
    let _ = std::fs::remove_dir(&dir);
}

#[test]
fn composite_and_tiered_warm_fan_out_per_image() {
    let (bytes_a, _) = four_segment_leb1();
    let (bytes_b, _) = four_segment_leb1();
    let index_a = FileIndex::decode(&bytes_a).expect("fixture A decodes");
    let index_b = FileIndex::decode(&bytes_b).expect("fixture B decodes");
    let composite = CompositeIndex::new(vec![
        (index_a.clone(), bytes_a.as_slice()),
        (index_b.clone(), bytes_b.as_slice()),
    ])
    .expect("composite builds");
    // One merged range per constituent image (offsets are per-file
    // address spaces and must never merge across files).
    assert_eq!(composite.warm(0.5, 2.5).unwrap(), vec![1, 1]);
    assert_eq!(composite.warm(10.0, 12.0).unwrap(), vec![0, 0]);
    composite.cool();
    let tiered = TieredIndex::new(vec![
        (
            LebTier::Base,
            CompositeIndex::new(vec![(index_a.clone(), bytes_a.as_slice())]).expect("base builds"),
        ),
        (
            LebTier::Medium,
            CompositeIndex::new(vec![(index_b.clone(), bytes_b.as_slice())])
                .expect("medium builds"),
        ),
    ])
    .expect("tiered builds");
    assert_eq!(tiered.warm(0.5, 2.5).unwrap(), vec![vec![1], vec![1]]);
    tiered.cool();
}
