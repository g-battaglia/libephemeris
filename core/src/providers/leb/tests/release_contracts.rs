//! Synthetic public LEB boundary and duplicate-directory contracts.

use super::*;

#[test]
fn delta_t_public_shape_guards_do_not_panic() {
    for (jds, vals) in [
        (&[0.0][..], &[][..]),
        (&[0.0][..], &[1.0, 2.0][..]),
        (&[0.0, 1.0][..], &[1.0][..]),
    ] {
        assert_eq!(
            eval_delta_t(jds, vals, 0.0).unwrap_err(),
            Error::LebCorrupt {
                what: "LEB Delta-T table (unequal lengths)".to_string(),
                offset: 0,
            }
        );
    }
    // The existing absent-table error retains priority, including NaN epochs.
    for vals in [&[][..], &[1.0][..]] {
        assert_eq!(
            eval_delta_t(&[], vals, f64::NAN).unwrap_err(),
            Error::LebCorrupt {
                what: "LEB Delta-T table (absent)".to_string(),
                offset: 0,
            }
        );
    }
}

#[test]
fn delta_t_public_nonfinite_epoch_and_singleton_guards() {
    for jd in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        match eval_delta_t(&[0.0], &[1.0], jd).unwrap_err() {
            Error::NonFiniteArgument { value } => assert_eq!(value.to_bits(), jd.to_bits()),
            other => panic!("expected epoch refusal, got {other:?}"),
        }
    }
    assert!(matches!(
        eval_delta_t(&[f64::NAN], &[1.0], 0.0),
        Err(Error::LebCorrupt { offset: 0, .. })
    ));
    // Shape refusal precedes the new epoch guard.
    assert!(matches!(
        eval_delta_t(&[0.0], &[], f64::NAN),
        Err(Error::LebCorrupt { offset: 0, .. })
    ));
}

#[test]
fn delta_t_valid_clamps_right_bisect_and_stored_words_remain() {
    for (jd, expected) in [
        (-5.0, 1.0),
        (0.0, 1.0),
        (5.0, 2.0),
        (10.0, 3.0),
        (20.0, 3.0),
    ] {
        assert_eq!(
            eval_delta_t(&[0.0, 10.0], &[1.0, 3.0], jd).unwrap(),
            expected
        );
    }
    // Endpoint clamp owns JD zero; interior search owns the last duplicate zero.
    assert_eq!(
        eval_delta_t(&[0.0, 0.0, 10.0], &[1.0, 2.0, 4.0], 0.0).unwrap(),
        1.0
    );
    assert_eq!(
        eval_delta_t(&[0.0, 0.0, 10.0], &[1.0, 2.0, 4.0], 5.0).unwrap(),
        3.0
    );
    for jd in [-f64::MAX, 0.0, f64::MAX] {
        assert_eq!(
            eval_delta_t(&[0.0], &[-0.0], jd).unwrap().to_bits(),
            (-0.0_f64).to_bits()
        );
    }
    assert!(
        eval_delta_t(&[0.0, 10.0], &[f64::NAN, 3.0], 5.0)
            .unwrap()
            .is_nan()
    );
}

#[test]
fn segment_selection_refuses_nonfinite_quotient_before_cast() {
    for (start, jd, expected) in [
        (-f64::MAX, f64::MAX, f64::INFINITY),
        (f64::MAX, -f64::MAX, f64::NEG_INFINITY),
    ] {
        assert_eq!(
            select_segment(start, 1.0, 1, jd).unwrap_err(),
            Error::NonFiniteArgument { value: expected }
        );
        assert!(matches!(
            select_segment_two_part(start, 1.0, 1, TwoPartEpoch::split(jd).unwrap()),
            Err(Error::NonFiniteArgument { .. })
        ));
    }
    assert!(matches!(
        select_segment(f64::NAN, 1.0, 1, 0.0),
        Err(Error::NonFiniteArgument { value }) if value.is_nan()
    ));
}

#[test]
fn segment_quotient_overflow_is_refused_inside_a_finite_body_range() {
    let body_at = HEADER_SIZE + SECTION_DIR_SIZE;
    let blob_at = body_at + LEB1_BODY_ENTRY_SIZE;
    let mut image = pack_header(&MAGIC_LEB1, 1, 1, 1, -f64::MAX, f64::MAX, 0.0, 0);
    image.extend_from_slice(&pack_section(
        SECTION_BODY_INDEX,
        body_at as u64,
        LEB1_BODY_ENTRY_SIZE as u64,
    ));
    image.extend_from_slice(&pack_body_leb1(
        0,
        0,
        1,
        -f64::MAX,
        f64::MAX,
        1.0,
        0,
        3,
        blob_at as u64,
    ));
    for _ in 0..3 {
        push_f64(&mut image, 1.0);
    }
    let index = FileIndex::decode(&image).unwrap();
    assert_eq!(index.body_coverage(0), Some((-f64::MAX, f64::MAX)));
    assert_eq!(
        index.eval_body(&image, 0, f64::MAX).unwrap_err(),
        Error::NonFiniteArgument {
            value: f64::INFINITY
        }
    );
}

#[test]
fn segment_selection_preserves_finite_quotients_and_infinite_tau_clamps() {
    for (jd, expected_tau) in [(f64::MAX, 1.0), (-f64::MAX, -1.0)] {
        let (segment, tau) = select_segment(0.0, 1.0, 2, jd).unwrap();
        assert_eq!(segment, usize::from(jd > 0.0));
        assert_eq!(tau.value(), expected_tau);
    }
    // Quotient zero is valid even when the ordered midpoint sum overflows.
    let (segment, tau) = select_segment(f64::MAX, f64::MAX, 1, f64::MAX).unwrap();
    assert_eq!((segment, tau.value()), (0, -1.0));
    for (jd, expected) in [(0.0, (0, -1.0)), (1.0, (1, -1.0)), (2.0, (1, 1.0))] {
        let (segment, tau) = select_segment(0.0, 1.0, 2, jd).unwrap();
        assert_eq!((segment, tau.value()), expected);
    }
    assert!(matches!(
        select_segment(0.0, 0.0, 0, f64::NAN),
        Err(Error::NonFiniteArgument { .. })
    ));
    assert!(matches!(
        select_segment(-f64::MAX, 0.0, 0, f64::MAX),
        Err(Error::LebCorrupt { .. })
    ));
}

/// Three physical entries: body zero, body five, then a replacement zero.
/// Degree-zero coefficients are input declarations, not native golden words.
fn duplicate_body_image(format: LebFormat) -> Vec<u8> {
    let is_leb1 = format == LebFormat::Leb1;
    let stride = if is_leb1 {
        LEB1_BODY_ENTRY_SIZE
    } else {
        LEB2_BODY_ENTRY_SIZE
    };
    let mut file = pack_header(
        if is_leb1 { &MAGIC_LEB1 } else { &MAGIC_LEB2 },
        if format == LebFormat::Leb2V2 { 2 } else { 1 },
        1,
        3,
        0.0,
        11.0,
        0.0,
        u32::from(!is_leb1),
    );
    let body_at = HEADER_SIZE + SECTION_DIR_SIZE;
    let mut blob_at = body_at + 3 * stride;
    file.extend_from_slice(&pack_section(
        SECTION_BODY_INDEX,
        body_at as u64,
        (3 * stride) as u64,
    ));
    let mut blobs = Vec::new();
    for (id, start, value) in [(0, 0.0, 7.0), (5, 2.0, 8.0), (0, 10.0, 9.0)] {
        let mut raw = Vec::new();
        for _ in 0..3 {
            push_f64(&mut raw, value);
        }
        if is_leb1 {
            file.extend_from_slice(&pack_body_leb1(
                id,
                0,
                1,
                start,
                start + 1.0,
                1.0,
                0,
                3,
                blob_at as u64,
            ));
            blob_at += raw.len();
            blobs.extend_from_slice(&raw);
        } else {
            let compressed = encode_frame(&raw, 1, 0, 3);
            let mut payload = Vec::new();
            if format == LebFormat::Leb2V2 {
                payload.extend_from_slice(&pack_chunk_header(1, 1.0));
                payload.extend_from_slice(&pack_chunk_entry(
                    start,
                    start + 1.0,
                    (blob_at + CHUNK_INDEX_HEADER_SIZE + CHUNK_ENTRY_SIZE) as u64,
                    compressed.len() as u64,
                    raw.len() as u64,
                    0,
                    1,
                ));
            }
            payload.extend_from_slice(&compressed);
            file.extend_from_slice(&pack_body_leb2(
                id,
                0,
                1,
                start,
                start + 1.0,
                1.0,
                0,
                3,
                blob_at as u64,
                compressed.len() as u64,
                raw.len() as u64,
            ));
            blob_at += payload.len();
            blobs.extend_from_slice(&payload);
        }
    }
    file.extend_from_slice(&blobs);
    file
}

#[test]
fn duplicate_body_last_record_owns_lookup_evaluation_and_inventory() {
    for format in [LebFormat::Leb1, LebFormat::Leb2V1, LebFormat::Leb2V2] {
        let image = duplicate_body_image(format);
        let index = FileIndex::decode(&image).unwrap();
        assert_eq!(index.header().body_count(), 3);
        assert_eq!(
            index
                .bodies()
                .iter()
                .map(BodyEntry::body_id)
                .collect::<Vec<_>>(),
            [0, 5]
        );
        assert_eq!(index.body_coverage(0), Some((10.0, 11.0)));
        assert_eq!(index.find_body(0).unwrap().jd_start(), 10.0);
        assert_eq!(index.body_coverage(99), None);
        assert_eq!(
            index.find_body(99).unwrap_err(),
            Error::LebBodyNotFound { body_id: 99 }
        );
        let composite = CompositeIndex::new(vec![(index.clone(), image.as_slice())]).unwrap();
        assert_eq!(
            composite.eval_body_any(0, 10.5).unwrap(),
            (vec![9.0; 3], vec![0.0; 3])
        );
        assert_eq!(
            composite.eval_body_any(5, 2.5).unwrap(),
            (vec![8.0; 3], vec![0.0; 3])
        );
        assert_eq!(
            composite.eval_body_any(0, 0.5).unwrap_err(),
            Error::LebOutOfRange {
                body_id: 0,
                jd: 0.5,
                start: 10.0,
                end: 11.0
            }
        );
        let info = inspect_bytes("base_core.leb2", &image, false).unwrap();
        assert_eq!(info.body_count(), 2);
        assert_eq!(
            info.bodies()
                .iter()
                .map(BodyRange::body_id)
                .collect::<Vec<_>>(),
            [0, 5]
        );
        assert_eq!(
            (info.bodies()[0].jd_start(), info.bodies()[0].jd_end()),
            (10.0, 11.0)
        );
    }
}

#[test]
fn duplicate_body_file_inspection_matches_memory_catalog() {
    use std::io::Write;
    for (serial, format) in [LebFormat::Leb1, LebFormat::Leb2V1, LebFormat::Leb2V2]
        .into_iter()
        .enumerate()
    {
        let image = duplicate_body_image(format);
        let path = std::env::temp_dir().join(format!(
            "liblephem-release-duplicate-{}-{serial}.leb",
            std::process::id()
        ));
        let mut file = std::fs::OpenOptions::new()
            .read(true)
            .write(true)
            .create_new(true)
            .open(&path)
            .unwrap();
        file.write_all(&image).unwrap();
        let result = FileIndex::read_from_file(&mut file);
        drop(file);
        std::fs::remove_file(path).unwrap();
        assert_eq!(result.unwrap(), FileIndex::decode(&image).unwrap());
    }
}

#[test]
fn duplicate_body_resolution_preserves_composite_and_tier_priority() {
    let first = duplicate_body_image(LebFormat::Leb1);
    let second = leb1_const_image(0, 99.0, 0.0, 20.0);
    let composite = CompositeIndex::new(vec![
        (FileIndex::decode(&first).unwrap(), first.as_slice()),
        (FileIndex::decode(&second).unwrap(), second.as_slice()),
    ])
    .unwrap();
    assert_eq!(composite.serving_index(0), Some(0));
    assert_eq!(
        composite.serving_entry(0),
        FileIndex::decode(&first).unwrap().find_body(0).ok()
    );
    assert_eq!(composite.body_coverage(0), Some((10.0, 11.0)));
    assert_eq!(composite.eval_body_any(0, 10.5).unwrap().0, vec![9.0; 3]);
    assert_eq!(
        composite.eval_body_any(0, 0.5).unwrap_err(),
        Error::LebOutOfRange {
            body_id: 0,
            jd: 0.5,
            start: 10.0,
            end: 11.0
        }
    );
    let base =
        CompositeIndex::new(vec![(FileIndex::decode(&first).unwrap(), first.as_slice())]).unwrap();
    let extended = CompositeIndex::new(vec![(
        FileIndex::decode(&second).unwrap(),
        second.as_slice(),
    )])
    .unwrap();
    let tiered =
        TieredIndex::new(vec![(LebTier::Extended, extended), (LebTier::Base, base)]).unwrap();
    assert_eq!(tiered.selected_tier(0, 10.5), Some(LebTier::Base));
    assert_eq!(tiered.eval_body_any(0, 10.5).unwrap().0, vec![9.0; 3]);
    assert_eq!(tiered.selected_tier(0, 0.5), Some(LebTier::Extended));
    assert_eq!(tiered.eval_body_any(0, 0.5).unwrap().0, vec![99.0; 3]);
}

#[test]
fn duplicate_star_last_complete_record_and_short_tail() {
    let mut payload = star_payload(&[42, 9, 42]);
    payload[4..12].copy_from_slice(&1.0_f64.to_le_bytes());
    payload[2 * STAR_ENTRY_SIZE + 4..2 * STAR_ENTRY_SIZE + 12]
        .copy_from_slice(&2.0_f64.to_le_bytes());
    payload.extend_from_slice(&[0xff; STAR_ENTRY_SIZE - 1]);
    let original = aux_image(vec![(SECTION_STARS, payload)]);
    for format in [LebFormat::Leb1, LebFormat::Leb2V1, LebFormat::Leb2V2] {
        let mut image = original.clone();
        if format != LebFormat::Leb1 {
            image[..4].copy_from_slice(&MAGIC_LEB2);
            image[40..44].copy_from_slice(&COMPRESSION_ZSTD_TRUNC_SHUFFLE.to_le_bytes());
            if format == LebFormat::Leb2V2 {
                image[4..8].copy_from_slice(&2_u32.to_le_bytes());
            }
        }
        let index = FileIndex::decode(&image).unwrap();
        assert_eq!(index.get_star(&image, 42).unwrap().ra_j2000_deg(), 2.0);
        assert_eq!(index.get_star(&image, 9).unwrap().ra_j2000_deg(), 10.5);
        assert_eq!(
            index.get_star(&image, 99).unwrap_err(),
            Error::LebStarNotFound { star_id: 99 }
        );
    }
}

#[test]
fn duplicate_star_resolution_preserves_first_catalog_carrier() {
    let mut first_payload = star_payload(&[42, 42]);
    first_payload[4..12].copy_from_slice(&1.0_f64.to_le_bytes());
    first_payload[STAR_ENTRY_SIZE + 4..STAR_ENTRY_SIZE + 12]
        .copy_from_slice(&2.0_f64.to_le_bytes());
    let first = aux_image(vec![(SECTION_STARS, first_payload)]);
    let mut second_payload = star_payload(&[42, 99]);
    second_payload[4..12].copy_from_slice(&9.0_f64.to_le_bytes());
    let second = aux_image(vec![(SECTION_STARS, second_payload)]);
    let composite = CompositeIndex::new(vec![
        (FileIndex::decode(&first).unwrap(), first.as_slice()),
        (FileIndex::decode(&second).unwrap(), second.as_slice()),
    ])
    .unwrap();
    assert_eq!(composite.get_star(42).unwrap().ra_j2000_deg(), 2.0);
    assert_eq!(
        composite.get_star(99).unwrap_err(),
        Error::LebStarNotFound { star_id: 99 }
    );
}
