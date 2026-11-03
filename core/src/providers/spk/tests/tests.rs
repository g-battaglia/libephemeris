use super::*;

// `SegSpec`, `put_f64`, `put_i32` and `build_image` live at module
// scope (test-only) so sibling-module registry tests can reuse the
// hand-built images through [`moon_test_kernel_bytes`].

fn one_seg() -> SegSpec {
    SegSpec {
        start: -1_000_000.0,
        end: 1_000_000.0,
        target: 499,
        center: 0,
        frame: 1,
        seg_type: 2,
        name: "MARS BARYCENTER",
    }
}

#[test]
fn source_empty_ssb_registration_preserves_strict_segment_guards_and_ownership() {
    let mut registry = SpkRegistry::new();
    registry
        .register_image("origin", const_image(499, 0, [1.0; 3]))
        .unwrap();
    assert!(matches!(
        registry.register_body(17, "origin", 0),
        Err(Error::SpkBodyNotFound { target: 0 })
    ));
    registry.register_body_or_ssb(17, "origin", 0).unwrap();
    let held = registry.clone();
    let mut captured = None;
    let zero = held
        .eval_geometric_ssb_au_two_part_jd_with_range(17, 2451545.0, 0.25, &mut captured)
        .unwrap();
    assert_eq!(zero, ([0.0; 3], [0.0; 3]));
    assert_eq!(
        held.eval_kernel_ssb_au_two_part_jd("origin", CENTER_SSB, 2451545.0, 0.25)
            .unwrap(),
        zero
    );
    assert!(captured.is_none());
    assert!(crate::positions::registered_spk_heliocentric_probe(&held, 2451545.0, 17).unwrap());
    assert!(matches!(
        held.eval_body_ssb(17, &PanicProvider, 0.0),
        Err(Error::SpkBodyNotFound { target: 0 })
    ));
    assert!(registry.unregister_body(17));
    assert_eq!(held.body_info(17).unwrap().naif_id, 0);
    registry
        .register_image("helio", const_image(499, 10, [1.0; 3]))
        .unwrap();
    assert!(matches!(
        registry.register_body_or_ssb(17, "helio", 0),
        Err(Error::SpkBodyNotFound { target: 0 })
    ));
    assert!(matches!(
        registry.register_body_or_ssb(17, "absent", 0),
        Err(Error::SpkFileNotFound { .. })
    ));
    assert!(matches!(
        registry.register_body_or_ssb(17, "origin", 777),
        Err(Error::SpkBodyNotFound { target: 777 })
    ));
}

#[test]
fn relative_state_cancels_identical_owned_centers_before_barycentric_roundoff() {
    let mut registry = SpkRegistry::new();
    for (path, target, center, value) in [
        ("parent", 3, 0, 100_000_000.0),
        ("earth", 399, 3, 1.0),
        ("moon", 301, 3, 2.0),
    ] {
        registry
            .register_image(path, const_image(target, center, [value; 3]))
            .unwrap();
        registry.register_body(target, path, target).unwrap();
    }
    let (position, velocity) = registry
        .eval_geometric_difference_au_two_part_jd(301, 399, J2000_JD_TDB, 0.0)
        .unwrap();
    // The shared barycenter cancels exactly: (parent+2)-(parent+1)=1 km.
    assert_eq!(position, [1.0 / AUNIT_TO_KM; 3]);
    assert_eq!(velocity, [0.0; 3]);
    assert!(matches!(
        registry.eval_geometric_difference_au_two_part_jd(301, 399, f64::NAN, 0.0),
        Err(Error::NonFiniteArgument { .. })
    ));
}

#[test]
fn relative_state_applies_observer_offset_after_shared_center_cancellation() {
    let mut registry = SpkRegistry::new();
    for (path, target, center, value) in [
        ("parent", 3, 0, 100_000_000.0),
        ("earth", 399, 3, 1.0),
        ("moon", 301, 3, 2.0),
    ] {
        registry
            .register_image(path, const_image(target, center, [value; 3]))
            .unwrap();
        registry.register_body(target, path, target).unwrap();
    }
    let offset = 0.0001;
    let (position, velocity) = registry
        .eval_geometric_difference_with_observer_au_two_part_jd(
            301,
            399,
            J2000_JD_TDB,
            0.0,
            Some(([offset; 3], [0.0; 3])),
        )
        .unwrap();
    let expected = 1.0 / AUNIT_TO_KM - offset;
    assert_eq!(position, [expected; 3]);
    assert_eq!(velocity, [0.0; 3]);
}

#[test]
fn au_chain_retains_replacement_and_fallback_validation() {
    let mut registry = SpkRegistry::new();
    registry
        .register_image("target", const_image(499, 0, [1.0; 3]))
        .unwrap();
    registry.register_body(4, "target", 499).unwrap();
    // Duplicate a valid descriptor, with an inconsistent declared center.
    let mut image = const_image(499, 0, [1.0; 3]);
    let start = RECORD_BYTES + 24;
    let descriptor = image[start..start + 40].to_vec();
    image[start + 40..start + 80].copy_from_slice(&descriptor);
    put_f64(&mut image, RECORD_BYTES + 16, 2.0, DafEndian::Little);
    put_i32(&mut image, start + 40 + 20, 10, DafEndian::Little);
    registry.register_image("target", image).unwrap();
    assert!(matches!(
        registry.eval_geometric_ssb_au_two_part_jd(4, J2000_JD_TDB, 0.0),
        Err(Error::SpkAmbiguousCenter { target: 499 })
    ));

    // The unregistered same-image center also needs whole-target validation.
    let mut image = const_image(499, 10, [1.0; 3]);
    let descriptor = image[start..start + 40].to_vec();
    for (i, center) in [(1, 0), (2, 3)] {
        let offset = start + i * 40;
        image[offset..offset + 40].copy_from_slice(&descriptor);
        put_i32(&mut image, offset + 16, 10, DafEndian::Little);
        put_i32(&mut image, offset + 20, center, DafEndian::Little);
    }
    put_f64(&mut image, RECORD_BYTES + 16, 3.0, DafEndian::Little);
    registry.register_image("target", image).unwrap();
    assert!(matches!(
        registry.eval_geometric_ssb_au_two_part_jd(4, J2000_JD_TDB, 0.0),
        Err(Error::SpkAmbiguousCenter { target: 10 })
    ));
}

#[test]
fn center_identity_never_confuses_library_and_naif_keys() {
    let mut registry = SpkRegistry::new();
    registry
        .register_image("parent", const_image(10, 0, [2.0; 3]))
        .unwrap();
    registry
        .register_image("child", const_image(499, 10, [3.0; 3]))
        .unwrap();
    registry
        .register_image("collision", const_image(599, 0, [1000.0; 3]))
        .unwrap();
    registry.register_body(0, "parent", 10).unwrap(); // ordinary SUN key
    registry.register_body(4, "child", 499).unwrap(); // ordinary MARS key
    registry.register_body(10, "collision", 599).unwrap(); // not a Sun registration
    let (pos, _) = registry
        .eval_geometric_ssb_au_two_part_jd(4, J2000_JD_TDB, 0.0)
        .unwrap();
    assert_eq!(pos, [2.0 / AUNIT_TO_KM + 3.0 / AUNIT_TO_KM; 3]);
}

#[test]
fn two_part_translated_epoch_preserves_fraction_and_analytic_rate() {
    // All powers-of-two and integer products below are exact. The fraction
    // is smaller than one ET ULP at this translated epoch, but remains a
    // source divmod offset: 2*(86400/2^36)/8 = 675/2^31.
    let jd = J2000_JD_TDB + 67_108_864.0;
    let init = 67_108_864.0 * SECONDS_PER_DAY;
    let rec = ChebRec {
        mid: init + 4.0,
        radius: 4.0,
        sets: vec![vec![1.0, 1.0]; 3],
    };
    let bytes = cheb_image(
        DafEndian::Little,
        TYPE_CHEBYSHEV_POS,
        301,
        0,
        init,
        init + 8.0,
        init,
        8.0,
        &[rec],
        None,
        "SPLIT FRACTION",
    );
    let segment = decode_only(&bytes, DafEndian::Little, 301).unwrap();
    let fraction = 1.0 / 68_719_476_736.0;
    assert_eq!(jd + fraction, jd);
    assert_eq!(
        eval_chebyshev(&segment, init).unwrap().position_km,
        [0.0; 3]
    );
    let state = eval_chebyshev_two_part_jd(&segment, jd, fraction).unwrap();
    assert_eq!(state.position_km, [675.0 / 2_147_483_648.0; 3]);
    assert_eq!(state.velocity_km_s, [0.25; 3]);
    let mut registry = SpkRegistry::new();
    registry.register_image("split", bytes).unwrap();
    registry.register_body(301, "split", 301).unwrap();
    assert_eq!(
        registry
            .eval_geometric_ssb_two_part_jd(301, jd, fraction)
            .unwrap(),
        state
    );
}

#[test]
fn two_part_exact_seams_neighbors_and_closed_coverage() {
    let jd = J2000_JD_TDB + 67_108_864.0;
    let init = 67_108_864.0 * SECONDS_PER_DAY;
    let interval = SECONDS_PER_DAY;
    let recs = [
        ChebRec {
            mid: init + interval / 2.0,
            radius: interval / 2.0,
            sets: vec![vec![10.0]; 3],
        },
        ChebRec {
            mid: init + 1.5 * interval,
            radius: interval / 2.0,
            sets: vec![vec![20.0]; 3],
        },
    ];
    let bytes = cheb_image(
        DafEndian::Little,
        TYPE_CHEBYSHEV_POS,
        499,
        0,
        init,
        init + 2.0 * interval,
        init,
        interval,
        &recs,
        None,
        "SPLIT SEAMS",
    );
    let index = SpkIndex::parse(&bytes).unwrap();
    let mut registry = SpkRegistry::new();
    registry
        .register_image("translated-seams", bytes.clone())
        .unwrap();
    registry
        .register_body(499, "translated-seams", 499)
        .unwrap();
    let native =
        |fraction| registry.eval_kernel_ssb_au_two_part_jd("translated-seams", 499, jd, fraction);
    for (fraction, expected) in [(0.0, 10.0), (1.0, 20.0), (2.0, 20.0)] {
        assert_eq!(native(fraction).unwrap().0, [expected / AUNIT_TO_KM; 3]);
    }
    // Directory acquisition must reject a negative split fraction even when
    // both collapsed JD and rounded ET would equal the translated start.
    let tiny_day = 1.0e-12;
    assert_eq!(jd - tiny_day, jd);
    assert_eq!(init - tiny_day * SECONDS_PER_DAY, init);
    assert!(matches!(
        native(-tiny_day),
        Err(Error::SpkOutOfRange { .. })
    ));
    for fraction in [f64::from_bits(2.0_f64.to_bits() + 1), 2.5] {
        assert!(matches!(native(fraction), Err(Error::SpkOutOfRange { .. })));
    }
    assert!(matches!(native(3.0), Err(Error::SpkOutOfRange { .. })));
    assert!(matches!(native(-0.5), Err(Error::SpkOutOfRange { .. })));
    assert!(matches!(
        registry.eval_geometric_ssb_au_two_part_jd(499, jd, 2.5),
        Err(Error::SpkOutOfRange { .. })
    ));
    let state =
        |fraction| eval_relative_two_part_jd(&bytes, DafEndian::Little, &index, 499, jd, fraction);
    assert_eq!(state(0.0).unwrap().0.position_km, [10.0; 3]);
    assert_eq!(
        state(f64::from_bits(1.0_f64.to_bits() - 1))
            .unwrap()
            .0
            .position_km,
        [10.0; 3]
    );
    assert_eq!(state(1.0).unwrap().0.position_km, [20.0; 3]);
    assert_eq!(
        state(f64::from_bits(1.0_f64.to_bits() + 1))
            .unwrap()
            .0
            .position_km,
        [20.0; 3]
    );
    assert_eq!(state(2.0).unwrap().0.position_km, [20.0; 3]);
    assert!(matches!(
        state(-f64::from_bits(1)),
        Err(Error::SpkOutOfRange { .. })
    ));
    assert!(matches!(
        state(f64::from_bits(2.0_f64.to_bits() + 1)),
        Err(Error::SpkOutOfRange { .. })
    ));
}

#[test]
fn registry_decoded_cache_reuses_storage_and_replacement_invalidates() {
    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<SpkRegistry>();
    let image = |value| {
        cheb_image(
            DafEndian::Little,
            TYPE_CHEBYSHEV_POS,
            499,
            0,
            0.0,
            8.0,
            0.0,
            8.0,
            &[ChebRec {
                mid: 4.0,
                radius: 4.0,
                sets: vec![vec![value]; 3],
            }],
            None,
            "CACHE",
        )
    };
    let mut registry = SpkRegistry::new();
    registry.register_image("cache", image(7.0)).unwrap();
    registry.register_body(499, "cache", 499).unwrap();
    let cache = Arc::clone(&registry.decoded["cache"]);
    assert!(cache.entries[0].get().is_none());
    let cold = registry.clone();
    let first = registry
        .eval_geometric_ssb_two_part_jd(499, J2000_JD_TDB, 0.0)
        .unwrap();
    assert_eq!(first.position_km, [7.0; 3]);
    let stored = cache.entries[0].get().unwrap() as *const _;
    for _ in 0..4 {
        assert_eq!(registry.eval_body_relative(499, 4.0).unwrap().0, first);
        assert_eq!(cache.entries[0].get().unwrap() as *const _, stored);
    }
    assert_eq!(registry, cold); // Cache state is not image identity.
    assert!(registry.register_image("cache", vec![0; 10]).is_err());
    assert!(Arc::ptr_eq(&cache, &registry.decoded["cache"]));
    registry.register_image("cache", image(9.0)).unwrap();
    assert!(!Arc::ptr_eq(&cache, &registry.decoded["cache"]));
    assert_eq!(
        registry.eval_body_relative(499, 4.0).unwrap().0.position_km,
        [9.0; 3]
    );
    registry.remove_kernel("cache");
    assert!(registry.decoded.is_empty());
}

#[test]
fn single_segment_little_endian_exact() {
    let img = build_image(DafEndian::Little, &[vec![one_seg()]]);
    let index = SpkIndex::parse(&img).expect("valid test image");
    assert_eq!(index.len(), 1);
    let seg = index.select(499, 0.0).unwrap().expect("covers J2000");
    let d = seg.descriptor;
    assert_eq!(d.start_et, -1_000_000.0);
    assert_eq!(d.end_et, 1_000_000.0);
    assert_eq!((d.target, d.center, d.frame, d.seg_type), (499, 0, 1, 2));
    assert!(d.initial_address >= 1 && d.final_address >= d.initial_address);
    assert_eq!(seg.name_str(), "MARS BARYCENTER");
    // Outside coverage is a clean miss, not an error.
    assert!(index.select(499, 2_000_000.0).unwrap().is_none());
    assert!(index.select(501, 0.0).unwrap().is_none());
    // Closed-interval endpoints are covered.
    assert!(index.select(499, -1_000_000.0).unwrap().is_some());
    assert!(index.select(499, 1_000_000.0).unwrap().is_some());
    assert_eq!(index.coverage(499), vec![(-1_000_000.0, 1_000_000.0)]);
    assert!(index.coverage(501).is_empty());
}

#[test]
fn big_endian_twin_decodes_identically() {
    let little = SpkIndex::parse(&build_image(DafEndian::Little, &[vec![one_seg()]]))
        .expect("little-endian image");
    let big = SpkIndex::parse(&build_image(DafEndian::Big, &[vec![one_seg()]]))
        .expect("big-endian image");
    assert_eq!(little, big);
}

#[test]
fn latest_segment_wins_on_overlap() {
    let early = SegSpec {
        start: 0.0,
        end: 100.0,
        target: 599,
        center: 0,
        frame: 1,
        seg_type: 2,
        name: "EARLY",
    };
    let late = SegSpec {
        start: 50.0,
        end: 150.0,
        target: 599,
        center: 0,
        frame: 1,
        seg_type: 3,
        name: "LATE",
    };
    let other = SegSpec {
        start: 0.0,
        end: 150.0,
        target: 699,
        center: 0,
        frame: 1,
        seg_type: 2,
        name: "OTHER",
    };
    let index = SpkIndex::parse(&build_image(DafEndian::Little, &[vec![early, late, other]]))
        .expect("valid test image");
    // Overlap: the later segment wins (NAIF "latest supersedes").
    assert_eq!(index.select(599, 75.0).unwrap().unwrap().name_str(), "LATE");
    // Exclusive regions keep their own segment.
    assert_eq!(
        index.select(599, 25.0).unwrap().unwrap().name_str(),
        "EARLY"
    );
    assert_eq!(
        index.select(599, 125.0).unwrap().unwrap().name_str(),
        "LATE"
    );
    // Shared endpoint (closed intervals both cover 50.0): later wins.
    assert_eq!(
        index
            .select(599, 50.0)
            .unwrap()
            .unwrap()
            .descriptor
            .seg_type,
        3
    );
    // Another target is unaffected by priority.
    assert_eq!(
        index.select(699, 75.0).unwrap().unwrap().name_str(),
        "OTHER"
    );
}

#[test]
fn empty_file_has_no_segments() {
    let img = build_image(DafEndian::Little, &[]);
    let index = SpkIndex::parse(&img).expect("empty image parses");
    assert!(index.is_empty());
    assert_eq!(index.len(), 0);
    assert!(index.select(499, 0.0).unwrap().is_none());
}

#[test]
fn empty_head_with_nonzero_tail_is_corrupt() {
    let mut img = build_image(DafEndian::Little, &[]);
    put_i32(&mut img, OFF_BWARD, 4, DafEndian::Little);
    assert!(matches!(
        SpkIndex::parse(&img),
        Err(Error::SpkCorrupt { .. })
    ));
}

#[test]
fn summary_chain_spans_two_records() {
    // 26 segments force a second summary record (25 fit in one).
    let mut group_a = Vec::new();
    let mut group_b = Vec::new();
    for i in 0..26 {
        let seg = SegSpec {
            start: f64::from(i) * 10.0,
            end: f64::from(i) * 10.0 + 10.0,
            target: 1000 + i,
            center: 0,
            frame: 1,
            seg_type: 2,
            name: "CHAIN",
        };
        if i < 25 {
            group_a.push(seg);
        } else {
            group_b.push(seg);
        }
    }
    let index = SpkIndex::parse(&build_image(DafEndian::Little, &[group_a, group_b]))
        .expect("chained image");
    assert_eq!(index.len(), 26);
    // A segment from the second record is reachable and file order
    // (priority) is preserved across the chain.
    let seg = index.select(1025, 255.0).unwrap().expect("covers");
    assert_eq!(seg.descriptor.target, 1025);
    assert_eq!(seg.name_str(), "CHAIN");
}

#[test]
fn malformed_images_are_rejected() {
    let good = build_image(DafEndian::Little, &[vec![one_seg()]]);
    // Helper: parse a mutated clone, expecting SpkCorrupt unless
    // noted.
    let corrupt = |img: &[u8]| {
        assert!(
            matches!(SpkIndex::parse(img), Err(Error::SpkCorrupt { .. })),
            "expected SpkCorrupt"
        );
    };

    // Bad identification word.
    let mut bad = good.clone();
    bad[0..8].copy_from_slice(b"NOT-A-DA");
    assert!(matches!(
        SpkIndex::parse(&bad),
        Err(Error::SpkUnknownFormat { .. })
    ));

    // CK-shaped summaries (ND = 2, NI = 2) are not SPK.
    let mut ck = good.clone();
    put_i32(&mut ck, OFF_NI, 2, DafEndian::Little);
    assert!(matches!(
        SpkIndex::parse(&ck),
        Err(Error::SpkUnsupportedNdNi { nd: 2, ni: 2 })
    ));

    // Unknown binary-format string.
    let mut fmt = good.clone();
    fmt[OFF_FORMAT..OFF_FORMAT + 8].copy_from_slice(b"VAX-FLOA");
    corrupt(&fmt);

    // Truncated file record.
    corrupt(&good[..1000]);

    // Summary head past the end of the image.
    let mut fwd = good.clone();
    put_i32(&mut fwd, OFF_FWARD, 99, DafEndian::Little);
    corrupt(&fwd);

    // Summary count above the record capacity.
    let mut nsum = good.clone();
    put_f64(&mut nsum, RECORD_BYTES + 16, 26.0, DafEndian::Little);
    corrupt(&nsum);

    // Non-integral control word.
    let mut next = good.clone();
    put_f64(&mut next, RECORD_BYTES, 2.5, DafEndian::Little);
    corrupt(&next);

    // NEXT self-loop: terminates as corruption, never hangs.
    let mut cyc = good.clone();
    put_f64(&mut cyc, RECORD_BYTES, 2.0, DafEndian::Little);
    corrupt(&cyc);

    // Non-finite segment epoch.
    let mut nan = good.clone();
    put_f64(&mut nan, RECORD_BYTES + 3 * 8, f64::NAN, DafEndian::Little);
    corrupt(&nan);

    // Inverted coverage.
    let mut inv = good.clone();
    put_f64(&mut inv, RECORD_BYTES + 3 * 8, 5.0, DafEndian::Little);
    put_f64(&mut inv, RECORD_BYTES + 4 * 8, -5.0, DafEndian::Little);
    corrupt(&inv);

    // Zero initial address (5th integer component).
    let mut addr = good.clone();
    put_i32(
        &mut addr,
        RECORD_BYTES + 3 * 8 + 16 + 16,
        0,
        DafEndian::Little,
    );
    corrupt(&addr);

    // Final address past the end of the image (6th component).
    let mut past = good.clone();
    put_i32(
        &mut past,
        RECORD_BYTES + 3 * 8 + 16 + 20,
        i32::MAX,
        DafEndian::Little,
    );
    corrupt(&past);

    // Missing name record (image ends after the summary record).
    corrupt(&good[..3 * RECORD_BYTES]);

    // Non-finite request epoch is a loud error, not a miss.
    let index = SpkIndex::parse(&good).expect("valid test image");
    assert!(matches!(
        index.select(499, f64::NAN),
        Err(Error::NonFiniteArgument { .. })
    ));
}

#[test]
fn coverage_merges_touching_but_not_disjoint() {
    let mk = |start: f64, end: f64| SegSpec {
        start,
        end,
        target: 5,
        center: 0,
        frame: 1,
        seg_type: 2,
        name: "COV",
    };
    let index = SpkIndex::parse(&build_image(
        DafEndian::Little,
        &[vec![
            mk(0.0, 10.0),
            mk(10.0, 20.0), // touches: merges
            mk(5.0, 15.0),  // overlaps: merges
            mk(100.0, 110.0),
        ]],
    ))
    .expect("valid test image");
    assert_eq!(index.coverage(5), vec![(0.0, 20.0), (100.0, 110.0)]);
}

#[test]
fn et_jd_conversions_are_exact_at_j2000() {
    assert_eq!(jd_tdb_to_et(J2000_JD_TDB).unwrap(), 0.0);
    assert_eq!(et_to_jd_tdb(0.0).unwrap(), J2000_JD_TDB);
    assert_eq!(et_to_jd_tdb(SECONDS_PER_DAY).unwrap(), J2000_JD_TDB + 1.0);
    assert_eq!(
        jd_tdb_to_et(J2000_JD_TDB + 0.5).unwrap(),
        0.5 * SECONDS_PER_DAY
    );
    assert!(matches!(
        jd_tdb_to_et(f64::INFINITY),
        Err(Error::NonFiniteArgument { .. })
    ));
    assert!(matches!(
        et_to_jd_tdb(f64::NAN),
        Err(Error::NonFiniteArgument { .. })
    ));
}

// -- Type 2/3 fixtures -------------------------------------------

struct ChebRec {
    mid: f64,
    radius: f64,
    sets: Vec<Vec<f64>>,
}

#[path = "subset_tests.rs"]
mod subset_tests;

/// Hand-built kernel image holding one type 2/3 segment: file
/// record, one summary record, one name record, then crafted
/// record data plus the 4-word directory. No real kernel bytes.
#[allow(clippy::too_many_arguments)]
fn cheb_image(
    endian: DafEndian,
    seg_type: i32,
    target: i32,
    center: i32,
    start_et: f64,
    end_et: f64,
    init: f64,
    intlen: f64,
    recs: &[ChebRec],
    dir_override: Option<[f64; 4]>,
    name: &str,
) -> Vec<u8> {
    let set_len = recs[0].sets[0].len();
    let rsize = 2 + recs[0].sets.len() * set_len;
    let data_words = recs.len() * rsize + 4;
    let n_records = 3 + data_words.div_ceil(WORDS_PER_RECORD);
    let mut img = vec![0u8; n_records * RECORD_BYTES];
    img[OFF_ID..OFF_ID + 8].copy_from_slice(SPK_FILE_ID);
    put_i32(&mut img, OFF_ND, SPK_ND, endian);
    put_i32(&mut img, OFF_NI, SPK_NI, endian);
    img[16..76].fill(b' ');
    let locfmt: &[u8; 8] = match endian {
        DafEndian::Little => b"LTL-IEEE",
        DafEndian::Big => b"BIG-IEEE",
    };
    img[OFF_FORMAT..OFF_FORMAT + 8].copy_from_slice(locfmt);
    put_i32(&mut img, OFF_FWARD, 2, endian);
    put_i32(&mut img, OFF_BWARD, 2, endian);

    let sbase = RECORD_BYTES;
    put_f64(&mut img, sbase + 16, 1.0, endian); // NSUM (NEXT/PREV stay 0)
    let off = sbase + 3 * 8;
    put_f64(&mut img, off, start_et, endian);
    put_f64(&mut img, off + 8, end_et, endian);
    let init_word = 3 * WORDS_PER_RECORD + 1;
    let vals = [
        target,
        center,
        1,
        seg_type,
        init_word as i32,
        (init_word + data_words - 1) as i32,
    ];
    for (k, v) in vals.iter().enumerate() {
        put_i32(&mut img, off + 16 + (k / 2) * 8 + (k % 2) * 4, *v, endian);
    }
    let nbase = 2 * RECORD_BYTES;
    img[nbase..nbase + SPK_NAME_BYTES].fill(b' ');
    let name_bytes = name.as_bytes();
    assert!(name_bytes.len() <= SPK_NAME_BYTES);
    img[nbase..nbase + name_bytes.len()].copy_from_slice(name_bytes);

    let dbase = 3 * RECORD_BYTES;
    let mut w = 0usize;
    let mut putw = |v: f64| {
        put_f64(&mut img, dbase + w * 8, v, endian);
        w += 1;
    };
    for rec in recs {
        putw(rec.mid);
        putw(rec.radius);
        for set in &rec.sets {
            assert_eq!(set.len(), set_len);
            for &c in set {
                putw(c);
            }
        }
    }
    for v in dir_override.unwrap_or([init, intlen, rsize as f64, recs.len() as f64]) {
        putw(v);
    }
    assert_eq!(w, data_words);
    img
}

fn decode_only(img: &[u8], endian: DafEndian, target: i32) -> Result<ChebSegment> {
    let index = SpkIndex::parse(img).expect("valid test image");
    let seg = index
        .segments_for(target)
        .next()
        .expect("target present")
        .clone();
    decode_chebyshev(img, &seg, endian)
}

#[test]
fn type2_linear_state_is_exact() {
    // Degree 1: X = 1000 + 50*tau with MID = 1000 s, RADIUS = 50 s.
    let rec = ChebRec {
        mid: 1000.0,
        radius: 50.0,
        sets: vec![vec![1000.0, 50.0], vec![0.0, -25.0], vec![300.0, 0.0]],
    };
    let img = cheb_image(
        DafEndian::Little,
        TYPE_CHEBYSHEV_POS,
        499,
        0,
        950.0,
        1050.0,
        950.0,
        100.0,
        &[rec],
        None,
        "LINEAR",
    );
    let seg = decode_only(&img, DafEndian::Little, 499).expect("decodes");
    // tau = 0: values exact, velocities are the scaled slopes.
    let s = eval_chebyshev(&seg, 1000.0).expect("evaluates");
    assert_eq!(s.position_km, [1000.0, 0.0, 300.0]);
    assert_eq!(s.velocity_km_s, [1.0, -0.5, 0.0]);
    // tau = 0.5 and tau = -1 (record edges): exact dyadic arithmetic.
    let s = eval_chebyshev(&seg, 1025.0).expect("evaluates");
    assert_eq!(s.position_km, [1025.0, -12.5, 300.0]);
    let s = eval_chebyshev(&seg, 950.0).expect("evaluates");
    assert_eq!(s.position_km, [950.0, 25.0, 300.0]);
    assert_eq!(s.velocity_km_s, [1.0, -0.5, 0.0]);
}

#[test]
fn type2_constant_has_zero_velocity() {
    let rec = ChebRec {
        mid: 0.0,
        radius: 64.0,
        sets: vec![vec![5.0], vec![-7.0], vec![11.0]],
    };
    let img = cheb_image(
        DafEndian::Little,
        TYPE_CHEBYSHEV_POS,
        5,
        0,
        -64.0,
        64.0,
        -64.0,
        128.0,
        &[rec],
        None,
        "CONST",
    );
    let seg = decode_only(&img, DafEndian::Little, 5).expect("decodes");
    let s = eval_chebyshev(&seg, 17.0).expect("evaluates");
    assert_eq!(s.position_km, [5.0, -7.0, 11.0]);
    assert_eq!(s.velocity_km_s, [0.0, 0.0, 0.0]);
}

#[test]
fn type2_quadratic_velocity_identity() {
    // X = T_2(tau) = 2*tau^2 - 1, RADIUS = 64 s: at tau = 0.5 the
    // value is -0.5 and dX/dt = 4*tau/RADIUS = 0.03125 km/s, all
    // exact in binary64.
    let rec = ChebRec {
        mid: 0.0,
        radius: 64.0,
        sets: vec![
            vec![0.0, 0.0, 1.0],
            vec![0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
        ],
    };
    let img = cheb_image(
        DafEndian::Little,
        TYPE_CHEBYSHEV_POS,
        5,
        0,
        -64.0,
        64.0,
        -64.0,
        128.0,
        &[rec],
        None,
        "QUAD",
    );
    let seg = decode_only(&img, DafEndian::Little, 5).expect("decodes");
    let s = eval_chebyshev(&seg, 32.0).expect("evaluates");
    assert_eq!(s.position_km[0], -0.5);
    assert_eq!(s.velocity_km_s[0], 0.03125);
}

#[test]
fn type3_velocity_sets_are_unscaled() {
    // RADIUS = 100 s (far from 1) proves the stored velocity sets
    // are used, with no 1/RADIUS scaling.
    let rec = ChebRec {
        mid: 0.0,
        radius: 100.0,
        sets: vec![
            vec![1.0],
            vec![2.0],
            vec![3.0],
            vec![4.0],
            vec![5.0],
            vec![6.0],
        ],
    };
    let img = cheb_image(
        DafEndian::Big,
        TYPE_CHEBYSHEV_POS_VEL,
        601,
        5,
        -100.0,
        100.0,
        -100.0,
        200.0,
        &[rec],
        None,
        "VEL",
    );
    let seg = decode_only(&img, DafEndian::Big, 601).expect("decodes");
    let s = eval_chebyshev(&seg, 33.0).expect("evaluates");
    assert_eq!(s.position_km, [1.0, 2.0, 3.0]);
    assert_eq!(s.velocity_km_s, [4.0, 5.0, 6.0]);
}

#[test]
fn record_selection_owns_actual_record_intervals_and_closed_endpoints() {
    let mk = |value: f64, midpoint: f64| ChebRec {
        mid: midpoint,
        radius: 50.0,
        sets: vec![vec![value], vec![0.0], vec![0.0]],
    };
    let img = cheb_image(
        DafEndian::Little,
        TYPE_CHEBYSHEV_POS,
        5,
        0,
        0.0,
        200.0,
        0.0,
        100.0,
        &[mk(10.0, 50.0), mk(20.0, 150.0)],
        None,
        "TWO",
    );
    let segment = decode_only(&img, DafEndian::Little, 5).expect("valid records");
    for (epoch, expected) in [
        (0.0, 10.0),
        (50.0, 10.0),
        (100.0, 20.0),
        (150.0, 20.0),
        (200.0, 20.0),
    ] {
        assert_eq!(
            eval_chebyshev(&segment, epoch).unwrap().position_km[0],
            expected
        );
    }
    for epoch in [-50.0, 250.0] {
        assert!(matches!(
            eval_chebyshev(&segment, epoch),
            Err(Error::SpkOutOfRange {
                start: 0.0,
                end: 200.0,
                ..
            })
        ));
    }
    assert!(matches!(
        eval_chebyshev(&segment, f64::NAN),
        Err(Error::NonFiniteArgument { .. })
    ));
}

#[test]
fn directory_and_record_defects_are_rejected() {
    let rec = || ChebRec {
        mid: 0.0,
        radius: 50.0,
        sets: vec![vec![1.0], vec![0.0], vec![0.0]],
    };
    let good = || {
        cheb_image(
            DafEndian::Little,
            TYPE_CHEBYSHEV_POS,
            5,
            0,
            -50.0,
            50.0,
            -50.0,
            100.0,
            &[rec()],
            None,
            "OK",
        )
    };
    // Wrong data type is a typed refusal, not corruption.
    let index = SpkIndex::parse(&good()).expect("valid test image");
    let seg = index.segments_for(5).next().unwrap().clone();
    let mut other = seg.clone();
    other.descriptor.seg_type = 9;
    assert!(matches!(
        decode_chebyshev(&good(), &other, DafEndian::Little),
        Err(Error::SpkUnsupportedType { seg_type: 9 })
    ));

    // RSIZE not divisible into 3 sets.
    let img = cheb_image(
        DafEndian::Little,
        TYPE_CHEBYSHEV_POS,
        5,
        0,
        -50.0,
        50.0,
        -50.0,
        100.0,
        &[rec()],
        Some([-50.0, 100.0, 7.0, 1.0]),
        "OK",
    );
    assert!(matches!(
        decode_only(&img, DafEndian::Little, 5),
        Err(Error::SpkCorrupt { .. })
    ));

    // N inconsistent with the addressed span.
    let img = cheb_image(
        DafEndian::Little,
        TYPE_CHEBYSHEV_POS,
        5,
        0,
        -50.0,
        50.0,
        -50.0,
        100.0,
        &[rec()],
        Some([-50.0, 100.0, 5.0, 2.0]),
        "OK",
    );
    assert!(matches!(
        decode_only(&img, DafEndian::Little, 5),
        Err(Error::SpkCorrupt { .. })
    ));

    // Non-positive record interval.
    let img = cheb_image(
        DafEndian::Little,
        TYPE_CHEBYSHEV_POS,
        5,
        0,
        -50.0,
        50.0,
        -50.0,
        0.0,
        &[rec()],
        None,
        "OK",
    );
    assert!(matches!(
        decode_only(&img, DafEndian::Little, 5),
        Err(Error::SpkCorrupt { .. })
    ));

    // Degenerate record radius fails at evaluation, not decode.
    let bad_rec = ChebRec {
        mid: 0.0,
        radius: 0.0,
        sets: vec![vec![1.0], vec![0.0], vec![0.0]],
    };
    let img = cheb_image(
        DafEndian::Little,
        TYPE_CHEBYSHEV_POS,
        5,
        0,
        -50.0,
        50.0,
        -50.0,
        100.0,
        &[bad_rec],
        None,
        "OK",
    );
    let seg = decode_only(&img, DafEndian::Little, 5).expect("decodes");
    assert!(matches!(
        eval_chebyshev(&seg, 0.0),
        Err(Error::SpkCorrupt { .. })
    ));
}

#[test]
fn km_au_converters_use_the_release_au() {
    let state = ChebState {
        position_km: [AUNIT_TO_KM, 0.0, -AUNIT_TO_KM],
        velocity_km_s: [0.0, 1.0, 0.0],
    };
    assert_eq!(state.position_au(), [1.0, 0.0, -1.0]);
    let v = state.velocity_au_per_day();
    assert_eq!(v[0], 0.0);
    assert!((v[1] - SECONDS_PER_DAY / AUNIT_TO_KM).abs() == 0.0);
}

fn local_asset(name: &str) -> Option<std::path::PathBuf> {
    if let Ok(dir) = std::env::var("LIBEPHEMERIS_DATA_DIR") {
        let p = std::path::PathBuf::from(dir).join(name);
        if p.exists() {
            return Some(p);
        }
    }
    if let Ok(home) = std::env::var("HOME") {
        let p = std::path::PathBuf::from(home)
            .join(".libephemeris")
            .join(name);
        if p.exists() {
            return Some(p);
        }
    }
    None
}

// -- Type 21 fixtures --------------------------------------------

struct EmdaRec {
    tl: f64,
    g: Vec<f64>,
    refpos: [f64; 3],
    refvel: [f64; 3],
    /// Difference table, `maxdim` rows × 3 columns.
    dt: Vec<[f64; 3]>,
    kqmax1: usize,
    kq: [usize; 3],
}

/// Hand-built kernel image holding one type 21 segment: file
/// record, one summary record, one name record, then `N` records,
/// the epoch table, the `N/100` directory, `MAXDIM`, and `N`.
/// `dir_override` replaces the trailing `[MAXDIM, N]` pair to
/// exercise directory defects. No real kernel bytes.
#[allow(clippy::too_many_arguments)]
fn emda_image(
    endian: DafEndian,
    target: i32,
    center: i32,
    start_et: f64,
    end_et: f64,
    maxdim: usize,
    recs: &[EmdaRec],
    epochs: &[f64],
    dir_override: Option<[f64; 2]>,
    name: &str,
) -> Vec<u8> {
    assert_eq!(recs.len(), epochs.len());
    let dlsize = 4 * maxdim + 11;
    let n = recs.len();
    let n_dir = n / 100;
    let data_words = n * dlsize + n + n_dir + 2;
    let n_records = 3 + data_words.div_ceil(WORDS_PER_RECORD);
    let mut img = vec![0u8; n_records * RECORD_BYTES];
    img[OFF_ID..OFF_ID + 8].copy_from_slice(SPK_FILE_ID);
    put_i32(&mut img, OFF_ND, SPK_ND, endian);
    put_i32(&mut img, OFF_NI, SPK_NI, endian);
    img[16..76].fill(b' ');
    let locfmt: &[u8; 8] = match endian {
        DafEndian::Little => b"LTL-IEEE",
        DafEndian::Big => b"BIG-IEEE",
    };
    img[OFF_FORMAT..OFF_FORMAT + 8].copy_from_slice(locfmt);
    put_i32(&mut img, OFF_FWARD, 2, endian);
    put_i32(&mut img, OFF_BWARD, 2, endian);

    let sbase = RECORD_BYTES;
    put_f64(&mut img, sbase + 16, 1.0, endian); // NSUM
    let off = sbase + 3 * 8;
    put_f64(&mut img, off, start_et, endian);
    put_f64(&mut img, off + 8, end_et, endian);
    let init_word = 3 * WORDS_PER_RECORD + 1;
    let vals = [
        target,
        center,
        1,
        TYPE_EMDA,
        init_word as i32,
        (init_word + data_words - 1) as i32,
    ];
    for (k, v) in vals.iter().enumerate() {
        put_i32(&mut img, off + 16 + (k / 2) * 8 + (k % 2) * 4, *v, endian);
    }
    let nbase = 2 * RECORD_BYTES;
    img[nbase..nbase + SPK_NAME_BYTES].fill(b' ');
    let name_bytes = name.as_bytes();
    assert!(name_bytes.len() <= SPK_NAME_BYTES);
    img[nbase..nbase + name_bytes.len()].copy_from_slice(name_bytes);

    let dbase = 3 * RECORD_BYTES;
    let mut w = 0usize;
    let mut putw = |v: f64| {
        put_f64(&mut img, dbase + w * 8, v, endian);
        w += 1;
    };
    for rec in recs {
        assert_eq!(rec.g.len(), maxdim);
        assert_eq!(rec.dt.len(), maxdim);
        putw(rec.tl);
        for &g in &rec.g {
            putw(g);
        }
        for c in 0..3 {
            putw(rec.refpos[c]);
            putw(rec.refvel[c]);
        }
        for c in 0..3 {
            for row in &rec.dt {
                putw(row[c]);
            }
        }
        putw(rec.kqmax1 as f64);
        for &k in &rec.kq {
            putw(k as f64);
        }
    }
    for &e in epochs {
        putw(e);
    }
    for i in 0..n_dir {
        putw(epochs[(i + 1) * 100 - 1]);
    }
    let [md, nn] = dir_override.unwrap_or([maxdim as f64, n as f64]);
    putw(md);
    putw(nn);
    assert_eq!(w, data_words);
    img
}

fn linear_emda_rec(tl: f64, maxdim: usize) -> EmdaRec {
    EmdaRec {
        tl,
        g: vec![10.0; maxdim],
        refpos: [100_000_000.0, 0.0, 0.0],
        refvel: [30.0, 0.0, 0.0],
        dt: vec![[0.0; 3]; maxdim],
        kqmax1: 3,
        kq: [0, 0, 0],
    }
}

fn decode_emda_only(img: &[u8], endian: DafEndian, target: i32) -> Result<EmdaSegment> {
    let index = SpkIndex::parse(img).expect("valid test image");
    let seg = index
        .segments_for(target)
        .next()
        .expect("target present")
        .clone();
    decode_emda(img, &seg, endian)
}

#[test]
fn emda_linear_extrapolation_is_exact() {
    // Zero difference table: the state is the reference state plus
    // linear motion, exactly (all dyadic values).
    let recs = vec![linear_emda_rec(1000.0, 2)];
    let img = emda_image(
        DafEndian::Little,
        20002060,
        10,
        0.0,
        2000.0,
        2,
        &recs,
        &[2000.0],
        None,
        "LINEAR",
    );
    let seg = decode_emda_only(&img, DafEndian::Little, 20002060).expect("decodes");
    // At the reference epoch the stored state returns.
    let s = eval_emda(&seg, 1000.0).expect("evaluates");
    assert_eq!(s.position_km, [100_000_000.0, 0.0, 0.0]);
    assert_eq!(s.velocity_km_s, [30.0, 0.0, 0.0]);
    // DELTA = 1000 s of linear motion, exact.
    let s = eval_emda(&seg, 2000.0).expect("evaluates");
    assert_eq!(s.position_km, [100_030_000.0, 0.0, 0.0]);
    assert_eq!(s.velocity_km_s, [30.0, 0.0, 0.0]);
}

#[test]
fn emda_quadratic_matches_hand_derivation() {
    // MAXDIM = 1, KQMAX1 = 3, DT = [[2, 0, 0]], G = [8], TL = 0,
    // REFPOS = 100, REFVEL = 3, evaluated at DELTA = 4. Tracing the
    // NAIF recurrence in 1-based indices (W starts [1, 1/2, 1/3],
    // FC[1] = 1, FC[2] = DELTA/G = 1/2, WC[1] = 1/2):
    //  position weights: W[3] <- FC[2]*W[2] - WC[1]*W[3] = 1/12;
    //  SUM = DT*W[2] = 2*(1/2) = 1; pos = 100 + 4*(3 + 4*1) = 128;
    //  velocity weights: W[2] <- FC[2]*W[1] - WC[1]*W[3];
    //  SUM = DT*W[1] = 2*1 = 2; vel = 3 + 4*2 = 11.
    // The read-back weights (W[2]/W[1]) are untouched exact halves,
    // so both results are exact in binary64.
    let rec = EmdaRec {
        tl: 0.0,
        g: vec![8.0],
        refpos: [100.0, 0.0, 0.0],
        refvel: [3.0, 1.0, -2.0],
        dt: vec![[2.0, 0.0, 0.0]],
        kqmax1: 3,
        kq: [1, 0, 0],
    };
    let img = emda_image(
        DafEndian::Little,
        20000001,
        10,
        0.0,
        100.0,
        1,
        &[rec],
        &[100.0],
        None,
        "QUAD",
    );
    let seg = decode_emda_only(&img, DafEndian::Little, 20000001).expect("decodes");
    let s = eval_emda(&seg, 4.0).expect("evaluates");
    assert_eq!(s.position_km[0], 128.0);
    assert_eq!(s.velocity_km_s[0], 11.0);
    // Zero-order components stay exactly linear.
    assert_eq!(s.position_km[1], 0.0 + 4.0 * 1.0);
    assert_eq!(s.position_km[2], 0.0 + 4.0 * -2.0);
    assert_eq!(s.velocity_km_s[1], 1.0);
    assert_eq!(s.velocity_km_s[2], -2.0);
}

#[test]
fn emda_record_selection_is_strictly_greater() {
    // Epochs [100, 200, 300]: a request exactly at an epoch moves
    // to the next record.
    let recs: Vec<EmdaRec> = [0.0, 100.0, 200.0]
        .iter()
        .map(|&tl| linear_emda_rec(tl, 1))
        .collect();
    let img = emda_image(
        DafEndian::Little,
        20000002,
        10,
        0.0,
        300.0,
        1,
        &recs,
        &[100.0, 200.0, 300.0],
        None,
        "SEL",
    );
    let seg = decode_emda_only(&img, DafEndian::Little, 20000002).expect("decodes");
    assert_eq!(emda_record_index(&seg, 50.0).unwrap(), 0);
    assert_eq!(emda_record_index(&seg, 100.0).unwrap(), 1);
    assert_eq!(emda_record_index(&seg, 199.999).unwrap(), 1);
    assert_eq!(emda_record_index(&seg, 200.0).unwrap(), 2);
    assert_eq!(emda_record_index(&seg, 300.0).unwrap(), 2);
    assert_eq!(emda_record_index(&seg, 1.0e9).unwrap(), 2);
    assert!(matches!(
        emda_record_index(&seg, f64::NAN),
        Err(Error::NonFiniteArgument { .. })
    ));
    assert!(matches!(
        eval_emda(&seg, f64::INFINITY),
        Err(Error::NonFiniteArgument { .. })
    ));
}

#[test]
fn emda_decode_defects_are_rejected() {
    let good = || {
        emda_image(
            DafEndian::Little,
            20000003,
            10,
            0.0,
            100.0,
            2,
            &[linear_emda_rec(0.0, 2)],
            &[100.0],
            None,
            "OK",
        )
    };
    // Wrong data type is a typed refusal, not corruption.
    let index = SpkIndex::parse(&good()).expect("valid test image");
    let seg = index.segments_for(20000003).next().unwrap().clone();
    let mut other = seg.clone();
    other.descriptor.seg_type = 2;
    assert!(matches!(
        decode_emda(&good(), &other, DafEndian::Little),
        Err(Error::SpkUnsupportedType { seg_type: 2 })
    ));

    // Directory defects: zero/oversize/non-integral MAXDIM or N.
    for dir in [
        [0.0, 1.0],
        [26.0, 1.0],
        [2.5, 1.0],
        [2.0, 0.0],
        [2.0, 1.5],
        [2.0, 2.0], // length mismatch
    ] {
        let img = emda_image(
            DafEndian::Little,
            20000003,
            10,
            0.0,
            100.0,
            2,
            &[linear_emda_rec(0.0, 2)],
            &[100.0],
            Some(dir),
            "OK",
        );
        assert!(
            matches!(
                decode_emda_only(&img, DafEndian::Little, 20000003),
                Err(Error::SpkCorrupt { .. })
            ),
            "dir {dir:?} must be rejected"
        );
    }

    // Decreasing epochs.
    let img = emda_image(
        DafEndian::Little,
        20000003,
        10,
        0.0,
        100.0,
        1,
        &[linear_emda_rec(0.0, 1), linear_emda_rec(50.0, 1)],
        &[100.0, 50.0],
        None,
        "OK",
    );
    assert!(matches!(
        decode_emda_only(&img, DafEndian::Little, 20000003),
        Err(Error::SpkCorrupt { .. })
    ));

    // Integration-order defects: KQMAX1 < 3, above MAXDIM + 1, or a
    // KQ entry above MAXDIM / non-integral.
    for (kqmax1, kq) in [
        (2, [0, 0, 0]),
        (4, [0, 0, 0]),
        (3, [2, 0, 0]),
        (3, [0, 0, 1]),
    ] {
        let mut rec = linear_emda_rec(0.0, 1);
        rec.kqmax1 = kqmax1;
        rec.kq = kq;
        // NOTE: kq = [0,0,1] with MAXDIM = 1 passes the range check
        // (1 <= min(1, 2)); it is exercised for shape, not rejection.
        let img = emda_image(
            DafEndian::Little,
            20000003,
            10,
            0.0,
            100.0,
            1,
            &[rec],
            &[100.0],
            None,
            "OK",
        );
        let result = decode_emda_only(&img, DafEndian::Little, 20000003);
        if kqmax1 == 3 && kq == [0, 0, 1] {
            assert!(result.is_ok());
        } else {
            assert!(
                matches!(result, Err(Error::SpkCorrupt { .. })),
                "kqmax1 {kqmax1} kq {kq:?} must be rejected"
            );
        }
    }
    // (3, [2,0,0]) above: KQ = 2 > MAXDIM = 1 is the rejection arm.

    // Non-integral KQ.
    let mut rec = linear_emda_rec(0.0, 1);
    rec.kq = [0, 0, 0];
    let mut img = emda_image(
        DafEndian::Little,
        20000003,
        10,
        0.0,
        100.0,
        1,
        &[rec],
        &[100.0],
        None,
        "OK",
    );
    // KQ[1] word: record base + 4*1+9 (0-based word 13 of the record).
    let dbase = 3 * RECORD_BYTES;
    put_f64(&mut img, dbase + 13 * 8, 1.5, DafEndian::Little);
    assert!(matches!(
        decode_emda_only(&img, DafEndian::Little, 20000003),
        Err(Error::SpkCorrupt { .. })
    ));

    // Zero stepsize in the used prefix fails at evaluation.
    let mut rec = linear_emda_rec(0.0, 1);
    rec.g = vec![0.0];
    let img = emda_image(
        DafEndian::Little,
        20000003,
        10,
        0.0,
        100.0,
        1,
        &[rec],
        &[100.0],
        None,
        "OK",
    );
    let seg = decode_emda_only(&img, DafEndian::Little, 20000003).expect("decodes");
    assert!(matches!(
        eval_emda(&seg, 50.0),
        Err(Error::SpkCorrupt { .. })
    ));

    // NaN coefficient fails at decode.
    let mut img = good();
    put_f64(
        &mut img,
        3 * RECORD_BYTES + 2 * 8,
        f64::NAN,
        DafEndian::Little,
    );
    assert!(matches!(
        decode_emda_only(&img, DafEndian::Little, 20000003),
        Err(Error::SpkCorrupt { .. })
    ));
}

#[test]
fn emda_directory_scales_to_hundreds_of_records() {
    // 150 records force one directory epoch (150/100 = 1).
    let recs: Vec<EmdaRec> = (0..150)
        .map(|i| linear_emda_rec(f64::from(i) * 10.0, 1))
        .collect();
    let epochs: Vec<f64> = (0..150).map(|i| f64::from(i) * 10.0 + 10.0).collect();
    let img = emda_image(
        DafEndian::Little,
        20000004,
        10,
        0.0,
        1500.0,
        1,
        &recs,
        &epochs,
        None,
        "BIG",
    );
    let seg = decode_emda_only(&img, DafEndian::Little, 20000004).expect("decodes");
    assert_eq!(emda_record_index(&seg, 1495.0).unwrap(), 149);
    // Record 149: TL = 1490, DELTA = 5 -> exact linear state.
    let s = eval_emda(&seg, 1495.0).expect("evaluates");
    assert_eq!(s.position_km, [100_000_150.0, 0.0, 0.0]);
    assert_eq!(s.velocity_km_s, [30.0, 0.0, 0.0]);
}

/// Exact reference-epoch identity on real Horizons data (gated on
/// the machine-local SPK cache): evaluating at a record's TL
/// reproduces its stored reference state bit-for-bit. Scans the
/// cache in sorted order for the first kernel with a type-21
/// segment; project-owned outputs only.
#[test]
fn emda_real_record_reference_identity() {
    let Some(cache) = local_asset("spk") else {
        eprintln!("SKIP: no local SPK cache; needs downloaded kernels");
        return;
    };
    let mut files: Vec<_> = std::fs::read_dir(&cache)
        .expect("cache readable")
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|x| x == "bsp"))
        .collect();
    files.sort();
    let mut found = false;
    for path in files {
        let Ok(bytes) = std::fs::read(&path) else {
            continue;
        };
        let Ok(index) = SpkIndex::parse(&bytes) else {
            continue;
        };
        let Some(seg) = index.segments_of_type(21).next().cloned() else {
            continue;
        };
        // DE-scale kernels ship little-endian; skip anything else
        // rather than guessing (the sweep pins its kernel anyway).
        let Ok(decoded) = decode_emda(&bytes, &seg, DafEndian::Little) else {
            continue;
        };
        let (pos, vel, tl) = emda_record_reference(&decoded, 0).expect("record 0");
        let state = eval_emda(&decoded, tl).expect("evaluates at TL");
        assert_eq!(state.position_km, pos);
        assert_eq!(state.velocity_km_s, vel);
        // And the last record too (different coefficients).
        let last = decoded.n_records - 1;
        let (pos, vel, tl) = emda_record_reference(&decoded, last).expect("last");
        let state = eval_emda(&decoded, tl).expect("evaluates at TL");
        assert_eq!(state.position_km, pos);
        assert_eq!(state.velocity_km_s, vel);
        found = true;
        break;
    }
    assert!(found, "no type-21 segment in the local SPK cache");
}

/// Sanity check (not a certified bound): real JPL DE data decoded
/// by this module puts the Earth-Moon barycenter ~1 AU from the
/// SSB and the Moon a few hundred thousand km from Earth. Bounds
/// are a priori physical margins, fixed before running.
#[test]
fn de440s_barycenter_distances_are_physical() {
    let Some(path) = local_asset("de440s.bsp") else {
        eprintln!("SKIP: de440s.bsp not installed; needs the data asset");
        return;
    };
    let bytes = std::fs::read(&path).expect("asset readable");
    let index = SpkIndex::parse(&bytes).expect("release asset parses");
    assert!(!index.is_empty());
    // DE440 ships Chebyshev segments only; anything else is a real
    // finding about the asset, reimplemented here rather than assumed.
    let mut types: Vec<i32> = index
        .segments_for(3)
        .chain(index.segments_for(399))
        .chain(index.segments_for(301))
        .map(|s| s.descriptor.seg_type)
        .collect();
    types.sort_unstable();
    types.dedup();
    assert!(
        types
            .iter()
            .all(|t| *t == TYPE_CHEBYSHEV_POS || *t == TYPE_CHEBYSHEV_POS_VEL),
        "unexpected segment types: {types:?}"
    );
    // First available body wins: Moon/Earth (geocentric), else
    // EMB/Earth (barycentric).
    let candidates: [(i32, i32, f64, f64); 3] = [
        (301, 3, 300_000.0, 500_000.0),
        (3, 0, 0.95 * AUNIT_TO_KM, 1.05 * AUNIT_TO_KM),
        (399, 0, 0.95 * AUNIT_TO_KM, 1.05 * AUNIT_TO_KM),
    ];
    let mut evaluated = false;
    for (target, center, lo, hi) in candidates {
        let Some(seg) = index
            .segments_for(target)
            .find(|s| s.descriptor.center == center)
        else {
            continue;
        };
        if !(seg.descriptor.start_et <= 0.0 && 0.0 <= seg.descriptor.end_et) {
            continue;
        }
        let endian = DafEndian::Little; // DE assets ship LTL-IEEE
        let decoded = decode_chebyshev(&bytes, seg, endian).expect("decodes");
        let state = eval_chebyshev(&decoded, 0.0).expect("evaluates at J2000");
        let r = state.position_km;
        let dist = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
        assert!(
            lo <= dist && dist <= hi,
            "target {target}: dist {dist} km outside [{lo}, {hi}]"
        );
        evaluated = true;
        break;
    }
    assert!(evaluated, "no candidate body segment covers J2000");
}

/// Provider-level check on the real DE asset (provider-profile tooling
/// scope): the NAIF-id
/// tooling registry reproduces the direct evaluation bit-exactly,
/// the Moon SSB chain composes exactly, and Sun/Moon distances sit
/// inside generous a priori physical margins. Bounds fixed before
/// running: the Sun stays within 2e6 km of the SSB (Jupiter alone
/// can displace it ~7.5e5 km; Saturn adds more), the Moon within
/// [300000, 500000] km of Earth (lunar range 356400–406700 km).
#[test]
fn de440s_provider_chain_is_exact_and_physical() {
    let Some(path) = local_asset("de440s.bsp") else {
        eprintln!("SKIP: de440s.bsp not installed; needs the data asset");
        return;
    };
    let path_str = path.to_string_lossy().into_owned();
    let mut registry = SpkRegistry::new();
    let bytes = std::fs::read(&path).expect("asset readable");
    registry
        .register_image(&path_str, bytes.clone())
        .expect("release asset parses");
    let registered =
        register_all_naif_targets(&mut registry, &path_str).expect("all targets resolve");
    for required in [3, 10, 301, 399] {
        assert!(
            registered.contains(&required),
            "de440s must serve NAIF {required}"
        );
    }
    let provider = NaifBaseProvider::new(&registry);
    for et in [0.0, 1_000_000.0] {
        // Registry path == direct path, bit-exactly.
        let (through_registry, _) = registry
            .eval_body_relative(301, et)
            .expect("Moon relative serves");
        let index = SpkIndex::parse(&bytes).expect("index");
        let (direct, _) = eval_relative(&bytes, DafEndian::Little, &index, 301, et)
            .expect("direct relative serves");
        assert_eq!(through_registry, direct, "registry adds no numerics");
        // Moon SSB == EMB SSB + Moon-relative, component-exact.
        let moon_ssb = registry
            .eval_body_ssb(301, &provider, et)
            .expect("chain serves");
        let (rel, _) = registry
            .eval_body_relative(301, et)
            .expect("relative serves");
        let emb_ssb = registry
            .eval_body_ssb(3, &provider, et)
            .expect("EMB serves");
        for i in 0..3 {
            assert_eq!(
                moon_ssb.position_km[i],
                emb_ssb.position_km[i] + rel.position_km[i],
                "composition is center-plus-relative"
            );
            assert_eq!(
                moon_ssb.velocity_km_s[i],
                emb_ssb.velocity_km_s[i] + rel.velocity_km_s[i],
                "composition is center-plus-relative"
            );
        }
    }
    // Physical margins at J2000 (independent evaluation, not a numerical comparison).
    let sun = registry
        .eval_body_ssb(10, &provider, 0.0)
        .expect("Sun serves");
    let sun_r = sun.position_km;
    let sun_dist = (sun_r[0] * sun_r[0] + sun_r[1] * sun_r[1] + sun_r[2] * sun_r[2]).sqrt();
    assert!(
        sun_dist < 2_000_000.0,
        "Sun-SSB distance {sun_dist} km exceeds the a priori margin"
    );
    let (moon_rel, _) = registry.eval_body_relative(301, 0.0).expect("Moon serves");
    let m = moon_rel.position_km;
    let moon_dist = (m[0] * m[0] + m[1] * m[1] + m[2] * m[2]).sqrt();
    assert!(
        (300_000.0..=500_000.0).contains(&moon_dist),
        "Moon-Earth distance {moon_dist} km outside the a priori margin"
    );
}

// -- chain fixtures -----------------------------------------

struct ConstProvider {
    state: ChebState,
}

impl CenterProvider for ConstProvider {
    fn center_state(&self, _center: i32, et: f64) -> Result<ChebState> {
        if !et.is_finite() {
            return Err(Error::NonFiniteArgument { value: et });
        }
        Ok(self.state)
    }
}

struct PanicProvider;

impl CenterProvider for PanicProvider {
    fn center_state(&self, _center: i32, _et: f64) -> Result<ChebState> {
        panic!("barycentric chain must not consult a provider");
    }
}

/// Constant type-2 state image: one degree-0 segment, so the
/// relative state is exactly `pos` with zero velocity.
fn const_image(target: i32, center: i32, pos: [f64; 3]) -> Vec<u8> {
    let rec = ChebRec {
        mid: 0.0,
        radius: 64.0,
        sets: vec![vec![pos[0]], vec![pos[1]], vec![pos[2]]],
    };
    cheb_image(
        DafEndian::Little,
        TYPE_CHEBYSHEV_POS,
        target,
        center,
        -64.0,
        64.0,
        -64.0,
        128.0,
        &[rec],
        None,
        "CONST",
    )
}

#[test]
fn geometric_source_composes_explicit_images_and_retains_them_after_clear() {
    let mut registry = SpkRegistry::new();
    for (path, target, center, position) in [
        ("minor", 2_000_001, 10, [1.0, 2.0, 3.0]),
        ("main", 10, 0, [4.0, 5.0, 6.0]),
    ] {
        registry
            .register_image(path, const_image(target, center, position))
            .unwrap();
        registry.register_body(target, path, target).unwrap();
    }
    let source = GeometricSpkSource::from_kernel(registry.clone(), "minor", 2_000_001).unwrap();
    registry.clear();
    assert_eq!(source.target(), 2_000_001);
    assert_eq!(source.paths(), &["minor", "main"]);
    assert_eq!(
        source.coverage_tdb(),
        &[(J2000_JD_TDB - 64.0 / 86400.0, J2000_JD_TDB + 64.0 / 86400.0)]
    );
    let (position, velocity) = source.state_tdb(J2000_JD_TDB).unwrap();
    for (actual, kilometers) in position.into_iter().zip([5.0, 7.0, 9.0]) {
        assert!((actual - kilometers / AUNIT_TO_KM).abs() < 1e-22);
    }
    assert_eq!(velocity, [0.0; 3]);
    assert_eq!(
        source.states_tdb(&[J2000_JD_TDB]).unwrap(),
        vec![[position[0], position[1], position[2], 0.0, 0.0, 0.0]]
    );
    assert!(source.states_tdb(&[]).unwrap().is_empty());
    for epoch in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            source.state_tdb(epoch),
            Err(Error::NonFiniteArgument { .. })
        ));
    }
}

#[test]
fn geometric_source_intersects_support_gaps_and_rejects_missing_targets() {
    let mut registry = SpkRegistry::new();
    registry
        .register_image("primary", const_image(499, 10, [1.0; 3]))
        .unwrap();
    let gap = |start, end| SegSpec {
        start,
        end,
        target: 10,
        center: 0,
        frame: 1,
        seg_type: 2,
        name: "CENTER GAP",
    };
    registry
        .register_image(
            "support",
            build_image(
                DafEndian::Little,
                &[vec![gap(-64.0, -16.0), gap(16.0, 64.0)]],
            ),
        )
        .unwrap();
    registry.register_body(10, "support", 10).unwrap();
    let source = GeometricSpkSource::from_kernel(registry.clone(), "primary", 499).unwrap();
    assert_eq!(
        source.coverage_tdb(),
        &[
            (J2000_JD_TDB - 64.0 / 86400.0, J2000_JD_TDB - 16.0 / 86400.0),
            (J2000_JD_TDB + 16.0 / 86400.0, J2000_JD_TDB + 64.0 / 86400.0),
        ]
    );
    assert!(matches!(
        source.state_tdb(J2000_JD_TDB),
        Err(Error::SpkOutOfRange { target: 10, .. })
    ));
    assert!(matches!(
        GeometricSpkSource::from_kernel(registry, "primary", 999),
        Err(Error::SpkBodyNotFound { target: 999 })
    ));
}

#[cfg(feature = "nbody")]
#[test]
fn nbody_seed_uses_primary_heliocentric_coordinates_without_a_supporting_sun() {
    let mut registry = SpkRegistry::new();
    for (path, target, center, position) in [
        ("seed", 2_000_001, 10, [1.0, 2.0, 3.0]),
        ("sun", 10, 0, [40.0, 50.0, 60.0]),
    ] {
        registry
            .register_image(path, const_image(target, center, position))
            .unwrap();
        registry.register_body(target, path, target).unwrap();
    }
    let seed = GeometricSpkSource::from_kernel(registry.clone(), "seed", 2_000_001).unwrap();
    let (position, velocity) = seed.heliocentric_seed_tdb(J2000_JD_TDB).unwrap();
    assert_eq!(position, [1.0, 2.0, 3.0].map(|v| v / AUNIT_TO_KM));
    assert_eq!(velocity, [0.0; 3]);
    let barycentric = GeometricSpkSource::from_kernel(registry, "sun", 10).unwrap();
    assert!(matches!(
        barycentric.heliocentric_seed_tdb(J2000_JD_TDB),
        Err(Error::SpkAmbiguousCenter { target: 10 })
    ));
    assert!(seed.heliocentric_seed_tdb(f64::NAN).is_err());
    assert!(seed.heliocentric_seed_tdb(J2000_JD_TDB + 1.0).is_err());
}

#[test]
fn geometric_source_never_replaces_primary_target_with_a_wider_image() {
    let mut registry = SpkRegistry::new();
    let mut primary = const_image(499, 0, [1.0; 3]);
    put_f64(&mut primary, RECORD_BYTES + 24, -16.0, DafEndian::Little);
    put_f64(&mut primary, RECORD_BYTES + 32, 16.0, DafEndian::Little);
    registry.register_image("primary", primary).unwrap();
    registry
        .register_image("alternate", const_image(499, 0, [1000.0; 3]))
        .unwrap();
    registry.register_body(499, "alternate", 499).unwrap();
    let source = GeometricSpkSource::from_kernel(registry, "primary", 499).unwrap();
    assert_eq!(source.paths(), &["primary"]);
    let (position, _) = source.state_tdb(J2000_JD_TDB).unwrap();
    assert_eq!(position, [1.0 / AUNIT_TO_KM; 3]);
    assert!(matches!(
        source.state_tdb(J2000_JD_TDB + 32.0 / 86400.0),
        Err(Error::SpkOutOfRange { target: 499, .. })
    ));
}

#[test]
fn chain_single_center_resolves() {
    let pair = |center| SegSpec {
        start: -100.0,
        end: 100.0,
        target: 499,
        center,
        frame: 1,
        seg_type: 2,
        name: "SEG",
    };
    let index = SpkIndex::parse(&build_image(DafEndian::Little, &[vec![pair(10), pair(10)]]))
        .expect("valid test image");
    assert_eq!(declared_center(&index, 499).unwrap(), 10);
    assert_eq!(check_frame(&index, 499).unwrap(), FRAME_J2000);
    assert_eq!(
        resolve_chain(&index, 499).unwrap(),
        TargetChain {
            target: 499,
            center: 10,
            frame: FRAME_J2000,
        }
    );
    assert_eq!(
        classify_center(0),
        CenterKind::Barycenter,
        "center 0 is the barycenter itself"
    );
    assert_eq!(classify_center(10), CenterKind::Sun);
    assert_eq!(classify_center(399), CenterKind::Other(399));
}

#[test]
fn chain_mixed_centers_rejected_at_registration_and_eval() {
    let seg = |center: i32, name: &'static str| SegSpec {
        start: -100.0,
        end: 100.0,
        target: 2002060,
        center,
        frame: 1,
        seg_type: 2,
        name,
    };
    let img = build_image(DafEndian::Little, &[vec![seg(0, "A"), seg(10, "B")]]);
    let index = SpkIndex::parse(&img).expect("valid test image");
    assert!(matches!(
        declared_center(&index, 2002060),
        Err(Error::SpkAmbiguousCenter { target: 2002060 })
    ));
    assert!(matches!(
        resolve_chain(&index, 2002060),
        Err(Error::SpkAmbiguousCenter { target: 2002060 })
    ));
    // The runtime gate fires before any segment is served.
    let provider = ConstProvider {
        state: ChebState {
            position_km: [0.0; 3],
            velocity_km_s: [0.0; 3],
        },
    };
    assert!(matches!(
        eval_ssb(&img, DafEndian::Little, &index, 2002060, &provider, 0.0),
        Err(Error::SpkAmbiguousCenter { target: 2002060 })
    ));
}

#[test]
fn chain_unknown_target_is_not_found() {
    let index = SpkIndex::parse(&build_image(DafEndian::Little, &[vec![one_seg()]]))
        .expect("valid test image");
    assert!(matches!(
        declared_center(&index, 501),
        Err(Error::SpkBodyNotFound { target: 501 })
    ));
    assert!(matches!(
        check_frame(&index, 501),
        Err(Error::SpkBodyNotFound { target: 501 })
    ));
    assert!(matches!(
        resolve_chain(&index, 501),
        Err(Error::SpkBodyNotFound { target: 501 })
    ));
    let provider = PanicProvider;
    assert!(matches!(
        eval_ssb(
            &build_image(DafEndian::Little, &[vec![one_seg()]]),
            DafEndian::Little,
            &index,
            501,
            &provider,
            0.0
        ),
        Err(Error::SpkBodyNotFound { target: 501 })
    ));
}

#[test]
fn chain_unsupported_frame_rejected() {
    let seg = |frame: i32, name: &'static str| SegSpec {
        start: -100.0,
        end: 100.0,
        target: 499,
        center: 0,
        frame,
        seg_type: 2,
        name,
    };
    // First offending segment in file order names the error.
    let img = build_image(DafEndian::Little, &[vec![seg(1, "OK"), seg(17, "BAD")]]);
    let index = SpkIndex::parse(&img).expect("valid test image");
    assert!(matches!(
        check_frame(&index, 499),
        Err(Error::SpkUnsupportedFrame { frame: 17 })
    ));
    assert!(matches!(
        resolve_chain(&index, 499),
        Err(Error::SpkUnsupportedFrame { frame: 17 })
    ));
    assert!(matches!(
        eval_ssb(&img, DafEndian::Little, &index, 499, &PanicProvider, 0.0),
        Err(Error::SpkUnsupportedFrame { frame: 17 })
    ));
}

#[test]
fn chain_ssb_center_needs_no_provider() {
    let img = const_image(499, 0, [5.0, -7.0, 11.0]);
    let index = SpkIndex::parse(&img).expect("valid test image");
    assert_eq!(
        resolve_chain(&index, 499).unwrap(),
        TargetChain {
            target: 499,
            center: 0,
            frame: FRAME_J2000,
        }
    );
    // A panicking provider proves the barycentric path never calls it.
    let served =
        eval_ssb(&img, DafEndian::Little, &index, 499, &PanicProvider, 17.0).expect("evaluates");
    let (relative, descriptor) =
        eval_relative(&img, DafEndian::Little, &index, 499, 17.0).expect("relative");
    assert_eq!(descriptor.center, CENTER_SSB);
    assert_eq!(served, relative, "SSB state is the relative state");
    assert_eq!(served.position_km, [5.0, -7.0, 11.0]);
    assert_eq!(served.velocity_km_s, [0.0, 0.0, 0.0]);
}

#[test]
fn chain_heliocentric_composition_is_exact() {
    let img = const_image(2002060, 10, [1.0, 2.0, 3.0]);
    let index = SpkIndex::parse(&img).expect("valid test image");
    let provider = ConstProvider {
        state: ChebState {
            position_km: [10.0, 20.0, 30.0],
            velocity_km_s: [1.0, 2.0, 3.0],
        },
    };
    let served =
        eval_ssb(&img, DafEndian::Little, &index, 2002060, &provider, -7.0).expect("evaluates");
    // Center-plus-relative, matching the release operand order.
    assert_eq!(served.position_km, [11.0, 22.0, 33.0]);
    assert_eq!(served.velocity_km_s, [1.0, 2.0, 3.0]);
}

#[test]
fn chain_composes_two_kernel_images_end_to_end() {
    // Target kernel: body rel. Sun; center kernel: Sun rel. SSB.
    let tgt = const_image(2002060, 10, [1_000.0, 0.0, 0.0]);
    let sun = const_image(10, 0, [150_000_000.0, 0.0, 0.0]);
    let tgt_index = SpkIndex::parse(&tgt).expect("target image");
    let sun_index = SpkIndex::parse(&sun).expect("center image");
    struct ImageProvider<'a> {
        bytes: &'a [u8],
        index: SpkIndex,
    }
    impl CenterProvider for ImageProvider<'_> {
        fn center_state(&self, center: i32, et: f64) -> Result<ChebState> {
            // Centers resolve against their own kernel's SSB chain.
            let (state, descriptor) =
                eval_relative(self.bytes, DafEndian::Little, &self.index, center, et)?;
            assert_eq!(descriptor.center, CENTER_SSB);
            Ok(state)
        }
    }
    let provider = ImageProvider {
        bytes: &sun,
        index: sun_index,
    };
    let served =
        eval_ssb(&tgt, DafEndian::Little, &tgt_index, 2002060, &provider, 3.0).expect("evaluates");
    assert_eq!(served.position_km, [150_001_000.0, 0.0, 0.0]);
    assert_eq!(served.velocity_km_s, [0.0, 0.0, 0.0]);
}

#[test]
fn naif_tooling_registers_all_targets_and_serves_two_level_chain() {
    // Moon-like two-level chain inside one registry: 301 rel. 3,
    // 3 rel. SSB. All bodies keyed by their NAIF id.
    let mut registry = SpkRegistry::new();
    registry
        .register_image("moon.bsp", const_image(301, 3, [1.0, 2.0, 3.0]))
        .expect("image loads");
    registry
        .register_image("emb.bsp", const_image(3, 0, [10.0, 20.0, 30.0]))
        .expect("image loads");
    assert_eq!(
        register_all_naif_targets(&mut registry, "moon.bsp"),
        Ok(vec![301])
    );
    assert_eq!(
        register_all_naif_targets(&mut registry, "emb.bsp"),
        Ok(vec![3])
    );
    let provider = NaifBaseProvider::new(&registry);
    let served = registry
        .eval_body_ssb(301, &provider, 3.0)
        .expect("chain serves");
    // Exact center-plus-relative composition through the registry path.
    assert_eq!(served.position_km, [11.0, 22.0, 33.0]);
    assert_eq!(served.velocity_km_s, [0.0, 0.0, 0.0]);
    // The registry relative path is bit-identical to the direct one.
    let (through_registry, _) = registry.eval_body_relative(301, 3.0).expect("relative");
    let bytes = registry
        .kernel("moon.bsp")
        .expect("kernel")
        .bytes()
        .to_vec();
    let index = SpkIndex::parse(&bytes).expect("index");
    let (direct, _) =
        eval_relative(&bytes, DafEndian::Little, &index, 301, 3.0).expect("direct relative");
    assert_eq!(through_registry, direct);
}

#[test]
fn naif_tooling_errors_are_typed() {
    // Unknown kernel key.
    let mut registry = SpkRegistry::new();
    assert!(matches!(
        register_all_naif_targets(&mut registry, "missing.bsp"),
        Err(Error::SpkFileNotFound { .. })
    ));
    // Mixed-center image is refused, not half-served.
    let seg = |center: i32| SegSpec {
        start: -100.0,
        end: 100.0,
        target: 7,
        center,
        frame: 1,
        seg_type: 2,
        name: "SEG",
    };
    registry
        .register_image(
            "mixed.bsp",
            build_image(DafEndian::Little, &[vec![seg(0), seg(10)]]),
        )
        .expect("image loads");
    assert!(matches!(
        register_all_naif_targets(&mut registry, "mixed.bsp"),
        Err(Error::SpkAmbiguousCenter { target: 7 })
    ));
    // Dangling center: registered target whose center has no body.
    registry
        .register_image("dangle.bsp", const_image(8, 9, [1.0, 0.0, 0.0]))
        .expect("image loads");
    register_all_naif_targets(&mut registry, "dangle.bsp").expect("registers 8");
    // A center cycle fails with a typed error, not a stack overflow.
    registry
        .register_image("cyc-a.bsp", const_image(5, 6, [1.0, 0.0, 0.0]))
        .expect("image loads");
    registry
        .register_image("cyc-b.bsp", const_image(6, 5, [0.0, 1.0, 0.0]))
        .expect("image loads");
    register_all_naif_targets(&mut registry, "cyc-a.bsp").expect("registers 5");
    register_all_naif_targets(&mut registry, "cyc-b.bsp").expect("registers 6");
    // All registrations done: one provider for the evaluation checks.
    let provider = NaifBaseProvider::new(&registry);
    assert!(matches!(
        registry.eval_body_ssb(8, &provider, 3.0),
        Err(Error::SpkUnregisteredBody { body: 9 })
    ));
    assert!(matches!(
        registry.eval_body_ssb(5, &provider, 3.0),
        Err(Error::SpkCorrupt { .. })
    ));
}

#[test]
fn chain_serves_chebyshev_velocity_sets_and_emda() {
    // Type 3: stored velocity sets pass through unscaled.
    let rec = ChebRec {
        mid: 0.0,
        radius: 100.0,
        sets: vec![
            vec![1.0],
            vec![2.0],
            vec![3.0],
            vec![4.0],
            vec![5.0],
            vec![6.0],
        ],
    };
    let t3 = cheb_image(
        DafEndian::Little,
        TYPE_CHEBYSHEV_POS_VEL,
        5,
        0,
        -50.0,
        50.0,
        -50.0,
        100.0,
        &[rec],
        None,
        "T3",
    );
    let t3_index = SpkIndex::parse(&t3).expect("type-3 image");
    let s = eval_ssb(&t3, DafEndian::Little, &t3_index, 5, &PanicProvider, 0.0).expect("evaluates");
    assert_eq!(s.position_km, [1.0, 2.0, 3.0]);
    assert_eq!(s.velocity_km_s, [4.0, 5.0, 6.0]);
    // Type 21: linear record evaluates through the chain untouched.
    let recs = vec![linear_emda_rec(1000.0, 2)];
    let e21 = emda_image(
        DafEndian::Little,
        20002060,
        0,
        0.0,
        2000.0,
        2,
        &recs,
        &[2000.0],
        None,
        "E21",
    );
    let e21_index = SpkIndex::parse(&e21).expect("type-21 image");
    let s = eval_ssb(
        &e21,
        DafEndian::Little,
        &e21_index,
        20002060,
        &PanicProvider,
        1000.0,
    )
    .expect("evaluates");
    assert_eq!(s.position_km, [100_000_000.0, 0.0, 0.0]);
    assert_eq!(s.velocity_km_s, [30.0, 0.0, 0.0]);
}

#[test]
fn chain_gap_reports_merged_span() {
    // Two disjoint descriptor-only segments: the gap epoch must
    // report the merged stored span, never extrapolate.
    let seg = |start: f64, end: f64, name: &'static str| SegSpec {
        start,
        end,
        target: 499,
        center: 0,
        frame: 1,
        seg_type: 2,
        name,
    };
    let img = build_image(
        DafEndian::Little,
        &[vec![seg(0.0, 10.0, "A"), seg(20.0, 30.0, "B")]],
    );
    let index = SpkIndex::parse(&img).expect("valid test image");
    assert!(matches!(
        eval_ssb(&img, DafEndian::Little, &index, 499, &PanicProvider, 15.0),
        Err(Error::SpkOutOfRange {
            target: 499,
            et: 15.0,
            start: 0.0,
            end: 30.0,
        })
    ));
    // Non-finite epochs are rejected before selection.
    assert!(matches!(
        eval_ssb(
            &img,
            DafEndian::Little,
            &index,
            499,
            &PanicProvider,
            f64::NAN
        ),
        Err(Error::NonFiniteArgument { .. })
    ));
}

#[test]
fn chain_unsupported_type_surfaces_through_eval() {
    let img = build_image(
        DafEndian::Little,
        &[vec![SegSpec {
            start: -100.0,
            end: 100.0,
            target: 499,
            center: 0,
            frame: 1,
            seg_type: 9,
            name: "T9",
        }]],
    );
    let index = SpkIndex::parse(&img).expect("valid test image");
    // Registration sees a consistent chain; evaluation refuses the type.
    assert!(resolve_chain(&index, 499).is_ok());
    assert!(matches!(
        eval_ssb(&img, DafEndian::Little, &index, 499, &PanicProvider, 0.0),
        Err(Error::SpkUnsupportedType { seg_type: 9 })
    ));
}

// -- Registration, selection, and lifetime  ------------------

/// Registry holding one barycentric constant image for `target`,
/// with `body` registered onto it.
fn registry_with_const(body: i32, target: i32, pos: [f64; 3]) -> (SpkRegistry, String) {
    let mut reg = SpkRegistry::new();
    let key = format!("const-{target}.bsp");
    reg.register_image(&key, const_image(target, CENTER_SSB, pos))
        .expect("const image registers");
    reg.register_body(body, &key, target)
        .expect("body registers");
    (reg, key)
}

#[test]
fn registry_body_lifecycle_and_list() {
    let mut reg = SpkRegistry::new();
    assert!(reg.is_empty());
    assert_eq!((reg.kernel_count(), reg.body_count()), (0, 0));
    assert!(reg.kernel_paths().is_empty());
    assert!(reg.registered_bodies().is_empty());

    reg.register_image("a.bsp", const_image(499, CENTER_SSB, [5.0, -7.0, 11.0]))
        .expect("image registers");
    assert!(!reg.is_empty());
    let kernel = reg.kernel("a.bsp").expect("stored");
    assert_eq!(kernel.path(), "a.bsp");
    assert_eq!(kernel.endian(), DafEndian::Little);
    assert_eq!(kernel.index().len(), 1);
    assert!(format!("{kernel:?}").contains("a.bsp"));

    reg.register_body(2060, "a.bsp", 499)
        .expect("body registers");
    assert_eq!(reg.body_count(), 1);
    assert_eq!(
        reg.body_info(2060),
        Some(BodyRegistration {
            body: 2060,
            kernel: "a.bsp".to_string(),
            naif_id: 499,
        })
    );
    assert_eq!(reg.body_info(9999), None);
    assert_eq!(reg.registered_bodies().len(), 1);

    // Barycentric target: the SSB path serves without a provider.
    let served = reg
        .eval_body_ssb(2060, &PanicProvider, 0.0)
        .expect("evaluates");
    assert_eq!(served.position_km, [5.0, -7.0, 11.0]);
    assert_eq!(served.velocity_km_s, [0.0, 0.0, 0.0]);
    let (relative, descriptor) = reg.eval_body_relative(2060, 0.0).expect("evaluates");
    assert_eq!(relative, served);
    assert_eq!((descriptor.target, descriptor.center), (499, CENTER_SSB));

    // Stored coverage in both native scales.
    assert_eq!(reg.coverage_et(2060).unwrap(), vec![(-64.0, 64.0)]);
    assert_eq!(
        reg.coverage_jd_tdb(2060).unwrap(),
        vec![(et_to_jd_tdb(-64.0).unwrap(), et_to_jd_tdb(64.0).unwrap())]
    );

    assert!(reg.unregister_body(2060));
    assert!(!reg.unregister_body(2060));
    assert_eq!(reg.body_count(), 0);
    assert!(matches!(
        reg.eval_body_ssb(2060, &PanicProvider, 0.0),
        Err(Error::SpkUnregisteredBody { body: 2060 })
    ));
    assert!(matches!(
        reg.coverage_et(2060),
        Err(Error::SpkUnregisteredBody { body: 2060 })
    ));
}

#[test]
fn registry_last_registration_wins() {
    let mut reg = SpkRegistry::new();
    reg.register_image("a.bsp", const_image(499, CENTER_SSB, [5.0, -7.0, 11.0]))
        .expect("image A registers");
    reg.register_image("b.bsp", const_image(499, CENTER_SSB, [1.0, 2.0, 3.0]))
        .expect("image B registers");
    reg.register_body(2060, "a.bsp", 499)
        .expect("first mapping");
    assert_eq!(
        reg.eval_body_ssb(2060, &PanicProvider, 0.0)
            .expect("evaluates")
            .position_km,
        [5.0, -7.0, 11.0]
    );
    // Re-registering the same body replaces the mapping outright:
    // no fallback to the earlier kernel.
    reg.register_body(2060, "b.bsp", 499)
        .expect("second mapping");
    assert_eq!(reg.body_info(2060).expect("info").kernel, "b.bsp");
    assert_eq!(
        reg.eval_body_ssb(2060, &PanicProvider, 0.0)
            .expect("evaluates")
            .position_km,
        [1.0, 2.0, 3.0]
    );
    // Replacing the image serves the new bytes under the same key;
    // the dependent registration needs no update.
    reg.register_image("b.bsp", const_image(499, CENTER_SSB, [9.0, 9.0, 9.0]))
        .expect("image replaces");
    assert_eq!(
        reg.eval_body_ssb(2060, &PanicProvider, 0.0)
            .expect("evaluates")
            .position_km,
        [9.0, 9.0, 9.0]
    );
}

#[test]
fn registry_rejects_unusable_kernels_at_registration() {
    let mut reg = SpkRegistry::new();
    // Not an SPK image at all.
    assert!(matches!(
        reg.register_image("bad.bsp", vec![0u8; RECORD_BYTES]),
        Err(Error::SpkUnknownFormat { .. })
    ));
    // Truncated file record.
    assert!(matches!(
        reg.register_image("short.bsp", vec![0u8; 100]),
        Err(Error::SpkCorrupt { .. })
    ));
    // An empty registry key can never resolve to a file.
    let valid = const_image(499, CENTER_SSB, [1.0, 2.0, 3.0]);
    assert!(matches!(
        reg.register_image("", valid),
        Err(Error::SpkFileNotFound { .. })
    ));
    // Nothing usable was stored by the failed attempts.
    assert!(reg.is_empty());

    reg.register_image("a.bsp", const_image(499, CENTER_SSB, [5.0, -7.0, 11.0]))
        .expect("image registers");
    // Unknown kernel key names the missing file.
    assert!(matches!(
        reg.register_body(2060, "missing.bsp", 499),
        Err(Error::SpkFileNotFound { .. })
    ));
    // NAIF id with no segment in this kernel.
    assert!(matches!(
        reg.register_body(2060, "a.bsp", 777),
        Err(Error::SpkBodyNotFound { target: 777 })
    ));
    assert_eq!(reg.body_count(), 0);

    // Mixed centers fail at registration, never per evaluation.
    let pair = |center: i32, name: &'static str| SegSpec {
        start: -100.0,
        end: 100.0,
        target: 499,
        center,
        frame: 1,
        seg_type: 2,
        name,
    };
    reg.register_image(
        "mixed.bsp",
        build_image(DafEndian::Little, &[vec![pair(0, "A"), pair(10, "B")]]),
    )
    .expect("summaries parse");
    assert!(matches!(
        reg.register_body(2060, "mixed.bsp", 499),
        Err(Error::SpkAmbiguousCenter { target: 499 })
    ));
    // Non-J2000 frame fails at registration too.
    reg.register_image(
        "fr17.bsp",
        build_image(
            DafEndian::Little,
            &[vec![SegSpec {
                start: -100.0,
                end: 100.0,
                target: 499,
                center: 0,
                frame: 17,
                seg_type: 2,
                name: "F17",
            }]],
        ),
    )
    .expect("summaries parse");
    assert!(matches!(
        reg.register_body(2060, "fr17.bsp", 499),
        Err(Error::SpkUnsupportedFrame { frame: 17 })
    ));
    // An unimplemented data type registers (chain is consistent)
    // but evaluation refuses it instead of inventing data.
    reg.register_image(
        "t9.bsp",
        build_image(
            DafEndian::Little,
            &[vec![SegSpec {
                start: -100.0,
                end: 100.0,
                target: 499,
                center: 0,
                frame: 1,
                seg_type: 9,
                name: "T9",
            }]],
        ),
    )
    .expect("summaries parse");
    reg.register_body(2060, "t9.bsp", 499)
        .expect("chain is consistent");
    assert!(matches!(
        reg.eval_body_ssb(2060, &PanicProvider, 0.0),
        Err(Error::SpkUnsupportedType { seg_type: 9 })
    ));
    reg.unregister_body(2060);
}

#[test]
fn registry_eval_errors_are_typed() {
    let (reg, _) = registry_with_const(2060, 499, [5.0, -7.0, 11.0]);
    // No mapping for this body: never a silent fallback.
    assert!(matches!(
        reg.eval_body_ssb(9999, &PanicProvider, 0.0),
        Err(Error::SpkUnregisteredBody { body: 9999 })
    ));
    // Outside the stored span: the merged span bounds the error.
    assert!(matches!(
        reg.eval_body_ssb(2060, &PanicProvider, 1000.0),
        Err(Error::SpkOutOfRange {
            target: 499,
            et: 1000.0,
            start: -64.0,
            end: 64.0,
        })
    ));
    // Non-finite epochs are rejected before selection.
    assert!(matches!(
        reg.eval_body_relative(2060, f64::NAN),
        Err(Error::NonFiniteArgument { .. })
    ));

    // A center the provider cannot resolve surfaces the provider's
    // typed error, never invented data.
    struct FailingProvider;
    impl CenterProvider for FailingProvider {
        fn center_state(&self, center: i32, _et: f64) -> Result<ChebState> {
            Err(Error::SpkBodyNotFound { target: center })
        }
    }
    let mut reg2 = SpkRegistry::new();
    reg2.register_image("helio.bsp", const_image(600, 10, [7.0, 0.0, 0.0]))
        .expect("image registers");
    reg2.register_body(6000, "helio.bsp", 600)
        .expect("body registers");
    assert!(matches!(
        reg2.eval_body_ssb(6000, &FailingProvider, 0.0),
        Err(Error::SpkBodyNotFound { target: 10 })
    ));
}

#[test]
fn registry_remove_and_clear_cascade() {
    let mut reg = SpkRegistry::new();
    reg.register_image("a.bsp", const_image(499, CENTER_SSB, [5.0, -7.0, 11.0]))
        .expect("image A registers");
    reg.register_image("b.bsp", const_image(501, CENTER_SSB, [1.0, 2.0, 3.0]))
        .expect("image B registers");
    reg.register_body(1, "a.bsp", 499).expect("body 1");
    reg.register_body(2, "b.bsp", 501).expect("body 2");
    assert_eq!(reg.kernel_paths(), vec!["a.bsp", "b.bsp"]);

    // Unknown key: no-op.
    assert!(!reg.remove_kernel("missing.bsp"));
    assert_eq!((reg.kernel_count(), reg.body_count()), (2, 2));
    // Dropping an image cascades to its dependent registrations;
    // the surviving kernel keeps serving.
    assert!(reg.remove_kernel("a.bsp"));
    assert_eq!((reg.kernel_count(), reg.body_count()), (1, 1));
    assert_eq!(reg.kernel_paths(), vec!["b.bsp"]);
    assert!(matches!(
        reg.eval_body_ssb(1, &PanicProvider, 0.0),
        Err(Error::SpkUnregisteredBody { body: 1 })
    ));
    assert_eq!(
        reg.eval_body_ssb(2, &PanicProvider, 0.0)
            .expect("survivor evaluates")
            .position_km,
        [1.0, 2.0, 3.0]
    );
    // Full teardown: images and registrations go as one unit.
    reg.clear();
    assert!(reg.is_empty());
    assert_eq!(reg.body_count(), 0);
    assert!(reg.registered_bodies().is_empty());
}

fn discovery_seg(target: i32, start: f64, end: f64, frame: i32) -> SegSpec {
    SegSpec {
        start,
        end,
        target,
        center: CENTER_SSB,
        frame,
        seg_type: 2,
        name: "DISCOVERY PROBE",
    }
}

fn discovery_image(target: i32, start: f64, end: f64, frame: i32) -> Vec<u8> {
    build_image(
        DafEndian::Little,
        &[vec![discovery_seg(target, start, end, frame)]],
    )
}

#[test]
fn discover_local_spks_registers_matches_and_names() {
    // Scratch directory: one minor body per NAIF-convention arm
    // plus skipped files.
    let dir = std::env::temp_dir().join(format!("libeph-spkdisc-{}-a", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("scratch dir");
    // Chiron under the legacy 2000000+N convention.
    std::fs::write(
        dir.join("a_chiron.bsp"),
        discovery_image(2_002_060, 0.0, 200.0 * 86_400.0, FRAME_J2000),
    )
    .expect("write chiron");
    // Pholus with an unserved frame: decode succeeds, registration
    // fails, so the status is the `"error:..."` leg.
    std::fs::write(
        dir.join("b_pholus.bsp"),
        discovery_image(2_005_145, 0.0, 200.0 * 86_400.0, 17),
    )
    .expect("write pholus");
    // Chariklo under the Horizons 20000000+N alias arm.
    std::fs::write(
        dir.join("c_chariklo.bsp"),
        discovery_image(20_010_199, 0.0, 200.0 * 86_400.0, FRAME_J2000),
    )
    .expect("write chariklo");
    // Varuna exercises the exotic-registry display-name fallback
    // (no hard-coded `SPK_BODY_NAMES` row).
    std::fs::write(
        dir.join("d_varuna.bsp"),
        discovery_image(2_020_000, 0.0, 200.0 * 86_400.0, FRAME_J2000),
    )
    .expect("write varuna");
    // Bennu exercises the -Horizons-id fallback (neither
    // hard-coded nor exotic name).
    std::fs::write(
        dir.join("e_bennu.bsp"),
        discovery_image(2_101_955, 0.0, 200.0 * 86_400.0, FRAME_J2000),
    )
    .expect("write bennu");
    // Undecodable image and non-.bsp names are skipped silently.
    std::fs::write(dir.join("junk.bsp"), b"not a kernel").expect("write junk");
    std::fs::write(dir.join("notes.txt"), b"ignore me").expect("write notes");

    let mut reg = SpkRegistry::new();
    let entries = discover_local_spks(&mut reg, &dir).expect("discovers");
    assert_eq!(
        entries,
        vec![
            SpkDiscoveryEntry {
                body_name: "Chiron".to_string(),
                status: "registered".to_string(),
            },
            SpkDiscoveryEntry {
                body_name: "Pholus".to_string(),
                status: "error: unsupported SPK reference frame: 17".to_string(),
            },
            SpkDiscoveryEntry {
                body_name: "Chariklo".to_string(),
                status: "registered".to_string(),
            },
            SpkDiscoveryEntry {
                body_name: "Varuna".to_string(),
                status: "registered".to_string(),
            },
            SpkDiscoveryEntry {
                body_name: "Bennu;".to_string(),
                status: "registered".to_string(),
            },
        ]
    );
    // Only the success legs install body mappings (Chiron 15,
    // Chariklo 20199, Varuna 30000, Bennu 111955); Pholus 16 stays
    // unmapped. Probe decodes stay side-effect-free: the junk and
    // text files never lodge images.
    assert!(reg.body_info(15).is_some());
    assert!(reg.body_info(16).is_none());
    assert!(reg.body_info(20199).is_some());
    assert!(reg.body_info(30000).is_some());
    assert!(reg.body_info(111955).is_some());
    assert_eq!(reg.body_count(), 4);
    assert_eq!(reg.kernel_count(), 5);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn discover_local_spks_reregisters_wider_spans_only() {
    let dir = std::env::temp_dir().join(format!("libeph-spkdisc-{}-b", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("scratch dir");
    std::fs::write(
        dir.join("n1_narrow.bsp"),
        discovery_image(2_002_060, 0.0, 100.0 * 86_400.0, FRAME_J2000),
    )
    .expect("write narrow");
    let mut reg = SpkRegistry::new();
    let first = discover_local_spks(&mut reg, &dir).expect("first discovers");
    assert_eq!(first.len(), 1);
    assert_eq!(first[0].status, "registered");
    let narrow_key = reg.body_info(15).expect("chiron mapped").kernel;

    // Repeat visit with no wider span: the name keeps its
    // first-seen outcome through the no-overwrite leg.
    let second = discover_local_spks(&mut reg, &dir).expect("second discovers");
    assert_eq!(second.len(), 1);
    assert_eq!(second[0].status, "already_registered");
    assert_eq!(reg.body_info(15).expect("still mapped").kernel, narrow_key);

    // A wider span re-registers and takes over the mapping.
    std::fs::write(
        dir.join("n2_wide.bsp"),
        discovery_image(2_002_060, 0.0, 1000.0 * 86_400.0, FRAME_J2000),
    )
    .expect("write wide");
    let third = discover_local_spks(&mut reg, &dir).expect("third discovers");
    assert_eq!(third.len(), 1);
    assert_eq!(third[0].status, "registered");
    let wide_key = reg.body_info(15).expect("remapped").kernel;
    assert_ne!(wide_key, narrow_key);

    // An equal-or-narrower span never displaces the wider mapping,
    // and the repeat status stays put (strict `>` comparison).
    std::fs::write(
        dir.join("n3_equal.bsp"),
        discovery_image(2_002_060, 0.0, 1000.0 * 86_400.0, FRAME_J2000),
    )
    .expect("write equal");
    let fourth = discover_local_spks(&mut reg, &dir).expect("fourth discovers");
    assert_eq!(fourth.len(), 1);
    assert_eq!(fourth[0].status, "already_registered");
    assert_eq!(reg.body_info(15).expect("wide kept").kernel, wide_key);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn discover_local_spks_empty_legs() {
    let mut reg = SpkRegistry::new();
    // Not a directory: no entries, never an error.
    let missing =
        std::env::temp_dir().join(format!("libeph-spkdisc-{}-absent", std::process::id()));
    let _ = std::fs::remove_dir_all(&missing);
    assert!(
        discover_local_spks(&mut reg, &missing)
            .expect("absent ok")
            .is_empty()
    );
    // A file is not a directory either.
    let file = std::env::temp_dir().join(format!("libeph-spkdisc-{}-file", std::process::id()));
    std::fs::write(&file, b"x").expect("scratch file");
    assert!(
        discover_local_spks(&mut reg, &file)
            .expect("file ok")
            .is_empty()
    );
    let _ = std::fs::remove_file(&file);
    // An empty directory discovers nothing and stores nothing.
    let empty = std::env::temp_dir().join(format!("libeph-spkdisc-{}-empty", std::process::id()));
    let _ = std::fs::remove_dir_all(&empty);
    std::fs::create_dir_all(&empty).expect("scratch dir");
    assert!(
        discover_local_spks(&mut reg, &empty)
            .expect("empty ok")
            .is_empty()
    );
    assert!(reg.is_empty());
    let _ = std::fs::remove_dir_all(&empty);
}

#[test]
fn registry_file_round_trip_and_path_resolution() {
    // Scratch kernel file under a unique temp directory.
    let dir = std::env::temp_dir().join(format!("libeph-spk06-{}", std::process::id()));
    let data_dir = dir.join("data");
    let spk_dir = data_dir.join("spk");
    std::fs::create_dir_all(&spk_dir).expect("scratch dirs");
    let direct = data_dir.join("direct.bsp");
    let nested = spk_dir.join("nested.bsp");
    std::fs::write(&direct, const_image(499, CENTER_SSB, [5.0, -7.0, 11.0])).expect("write direct");
    std::fs::write(&nested, const_image(501, CENTER_SSB, [1.0, 2.0, 3.0])).expect("write nested");

    // Candidate order: `<data>/spk/<name>` before `<data>/<name>`.
    let dir_s = data_dir.to_string_lossy().into_owned();
    assert_eq!(
        candidate_spk_paths("x.bsp", Some(&dir_s)),
        vec![format!("{dir_s}/spk/x.bsp"), format!("{dir_s}/x.bsp"),]
    );
    // Absolute requests pass through untouched.
    assert_eq!(
        candidate_spk_paths("/abs/x.bsp", Some(&dir_s)),
        vec!["/abs/x.bsp"]
    );
    // No data directory: no candidates, never an error.
    assert!(candidate_spk_paths("x.bsp", None).is_empty());

    // Resolution finds the direct file via the fallback candidate
    // and the nested file via the first candidate.
    assert_eq!(
        resolve_spk_path("direct.bsp", Some(&dir_s)),
        Some(direct.to_string_lossy().into_owned())
    );
    assert_eq!(
        resolve_spk_path("nested.bsp", Some(&dir_s)),
        Some(nested.to_string_lossy().into_owned())
    );
    assert_eq!(resolve_spk_path("absent.bsp", Some(&dir_s)), None);

    // File registration reads, parses eagerly, and serves.
    let mut reg = SpkRegistry::new();
    let key = reg.register_file(&nested).expect("file registers");
    assert_eq!(key, nested.to_string_lossy());
    reg.register_body(2060, &key, 501).expect("body registers");
    assert_eq!(
        reg.eval_body_ssb(2060, &PanicProvider, 0.0)
            .expect("evaluates")
            .position_km,
        [1.0, 2.0, 3.0]
    );
    // A missing file names itself in the error.
    let missing = data_dir.join("absent.bsp");
    assert!(matches!(
        reg.register_file(&missing),
        Err(Error::SpkFileNotFound { .. })
    ));

    // Changing the configured data directory does not disturb
    // installed images: owned bytes keep serving.
    assert_eq!(resolve_spk_path("nested.bsp", Some("/elsewhere")), None);
    assert_eq!(
        reg.eval_body_ssb(2060, &PanicProvider, 0.0)
            .expect("still evaluates")
            .position_km,
        [1.0, 2.0, 3.0]
    );

    std::fs::remove_dir_all(&dir).expect("scratch cleanup");
}

#[test]
fn spk_cache_dir_resolver_follows_precedence() {
    // Programmatic > env > TOML > default-None.
    assert_eq!(resolve_spk_cache_dir(None, None, None), None);
    assert_eq!(resolve_spk_cache_dir(Some(""), Some(""), Some("")), None);
    assert_eq!(
        resolve_spk_cache_dir(Some("/prog"), Some("/env"), Some("/toml")),
        Some("/prog".to_string())
    );
    assert_eq!(
        resolve_spk_cache_dir(None, Some("/env"), Some("/toml")),
        Some("/env".to_string())
    );
    assert_eq!(
        resolve_spk_cache_dir(None, None, Some("/toml")),
        Some("/toml".to_string())
    );
    // An empty override level falls through to the next one.
    assert_eq!(
        resolve_spk_cache_dir(Some(""), Some("/env"), Some("/toml")),
        Some("/env".to_string())
    );
}

#[test]
fn icrs_rotation_pins_release_order() {
    assert_eq!(OBLIQUITY_J2000_DEG, 23.4392911);
    // The x axis passes through untouched.
    assert_eq!(icrs_to_ecliptic_j2000([3.0, 0.0, 0.0]), [3.0, 0.0, 0.0]);
    let (se, ce) = OBLIQUITY_J2000_DEG.to_radians().sin_cos();
    // y -> (ce, -se), z -> (se, ce): the release component order.
    assert_eq!(icrs_to_ecliptic_j2000([0.0, 1.0, 0.0]), [0.0, ce, -se]);
    assert_eq!(icrs_to_ecliptic_j2000([0.0, 0.0, 1.0]), [0.0, se, ce]);
    // A general vector keeps its norm within 2 ulps (pure rotation).
    let v = [1.0, -2.5, 4.25];
    let w = icrs_to_ecliptic_j2000(v);
    let (n0, n1) = (
        (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt(),
        (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt(),
    );
    let ulp = f64::EPSILON * n0;
    let diff = (n1 - n0).abs();
    assert!(
        diff <= 2.0 * ulp,
        "rotation changed norm by {diff} (> 2 ulp {ulp})"
    );
}

#[test]
fn chain_real_assets_use_j2000_single_center_chains() {
    // Structural pin on the installed release assets (read-only):
    // every stored segment is J2000 and every target resolves one
    // center — the posture eval_ssb relies on.
    let mut checked = 0;
    for name in ["de440s.bsp", "10199_199001_201001.bsp"] {
        let Some(path) = local_asset(name) else {
            eprintln!("SKIP: {name} not installed; needs the data asset");
            continue;
        };
        let bytes = std::fs::read(&path).expect("asset readable");
        let index = SpkIndex::parse(&bytes).expect("release asset parses");
        assert!(!index.is_empty());
        let mut targets: Vec<i32> = index.all_segments().map(|s| s.descriptor.target).collect();
        targets.sort_unstable();
        targets.dedup();
        for target in targets {
            for seg in index.segments_for(target) {
                assert_eq!(
                    seg.descriptor.frame, FRAME_J2000,
                    "{name} target {target} names a non-J2000 frame"
                );
            }
            resolve_chain(&index, target).expect("single-center J2000 chain");
            checked += 1;
        }
    }
    // The Horizons Chariklo kernel of the sweep declares the
    // Sun (center 10) for its target, matching that evidence. It
    // lives in the `spk` cache subdirectory when installed.
    let chariklo = local_asset("10199_199001_201001.bsp").or_else(|| {
        std::env::var("HOME").ok().and_then(|home| {
            let p = std::path::PathBuf::from(home)
                .join(".libephemeris")
                .join("spk")
                .join("10199_199001_201001.bsp");
            p.exists().then_some(p)
        })
    });
    if let Some(path) = chariklo {
        let bytes = std::fs::read(&path).expect("asset readable");
        let index = SpkIndex::parse(&bytes).expect("parses");
        assert_eq!(declared_center(&index, 20010199).unwrap(), CENTER_SUN);
        eprintln!("chariklo center pin checked at {}", path.display());
    } else {
        eprintln!("SKIP: chariklo kernel not installed; needs the data asset");
    }
    eprintln!("chain asset pins checked for {checked} targets");
}

// ---: source-helper and auto-resolution pins -----------------
// Expectations below are hand-derived from the release source read
// (docstring examples and branch rules in `spk.py`/`spk_auto.py` at
// the pin); the release-agreement sweep over these helpers is
// uncommitted scratch evidence (see the evidence log).

#[test]
fn extract_asteroid_number_docstring_cases() {
    assert_eq!(extract_asteroid_number("2060"), Some(2060));
    assert_eq!(extract_asteroid_number("2060 Chiron"), Some(2060));
    assert_eq!(extract_asteroid_number("(136199) Eris"), Some(136199));
    assert_eq!(extract_asteroid_number("Chiron"), None);
    assert_eq!(extract_asteroid_number("2003 UB313"), None);
    assert_eq!(extract_asteroid_number("1P/Halley"), None);
    // Adjacent branches from the same source rules.
    assert_eq!(extract_asteroid_number("(2060)"), Some(2060));
    assert_eq!(extract_asteroid_number("1998 KY26"), None);
    assert_eq!(extract_asteroid_number("C/1995 O1"), None);
    assert_eq!(extract_asteroid_number("2014 MU69"), None);
    // "2003 Harding" is not provisional (lowercase continuation) so
    // the numbered branch claims the leading number, as release does.
    assert_eq!(extract_asteroid_number("2003 Harding"), Some(2003));
    assert_eq!(extract_asteroid_number("  2060 Chiron  "), Some(2060));
    assert_eq!(extract_asteroid_number(""), None);
    assert_eq!(extract_asteroid_number("(136199)"), Some(136199));
    // Six-digit parenthesized numbers are not four-digit provisionals.
    assert_eq!(extract_asteroid_number("(2003) UB313"), None);
    // Unbounded digit runs cannot become NAIF ids downstream.
    assert_eq!(extract_asteroid_number("99999999999999999999"), None);
}

#[test]
fn deduce_naif_id_offsets() {
    assert_eq!(deduce_naif_id("2060 Chiron", None), Some(2_002_060));
    assert_eq!(deduce_naif_id("anything", Some(2060)), Some(2_002_060));
    assert_eq!(deduce_naif_id("Chiron", None), None);
    assert_eq!(deduce_naif_id("2003 UB313", None), None);
    assert_eq!(deduce_naif_id("x", Some(i64::MAX)), None);
}

#[test]
fn naif_convention_selection_prefers_horizons() {
    assert_eq!(
        find_naif_id_for_asteroid(&[2_002_060, 20_002_060], 2060),
        Some(20_002_060)
    );
    assert_eq!(
        find_naif_id_for_asteroid(&[2_002_060], 2060),
        Some(2_002_060)
    );
    assert_eq!(find_naif_id_for_asteroid(&[42], 2060), None);
    assert_eq!(find_naif_id_for_asteroid(&[], 2060), None);
}

#[test]
fn detect_naif_id_prefers_single_small_body() {
    assert_eq!(
        detect_naif_id_from_targets(&[2_002_060, 20_002_060], Some(2_002_060)),
        Some(2_002_060)
    );
    assert_eq!(
        detect_naif_id_from_targets(&[399, 20_002_060], None),
        Some(20_002_060)
    );
    // Two segments sharing one target still count twice release-side.
    assert_eq!(
        detect_naif_id_from_targets(&[20_010_199, 20_010_199], None),
        None
    );
    assert_eq!(detect_naif_id_from_targets(&[399, 301], None), None);
    assert_eq!(detect_naif_id_from_targets(&[], None), None);
    assert_eq!(detect_naif_id_from_targets(&[], Some(5)), None);
}

#[test]
fn sanitize_filename_shared_rule() {
    // The docstring pin: writer and lookup agree on "Ceres;".
    assert_eq!(sanitize_filename("Ceres;"), "ceres");
    assert_eq!(sanitize_filename("2060 Chiron"), "2060_chiron");
    assert_eq!(sanitize_filename("(136199) Eris"), "136199_eris");
    assert_eq!(sanitize_filename("2060"), "2060");
    assert_eq!(sanitize_filename("a--b"), "a--b");
    assert_eq!(sanitize_filename(""), "body");
    assert_eq!(sanitize_filename("___"), "body");
    assert_eq!(sanitize_filename("a  b"), "a_b");
}

#[test]
fn auto_spk_stem_keeps_runs() {
    // The AutoSpkConfig rule does not collapse or strip.
    assert_eq!(auto_spk_cache_stem("Ceres;"), "ceres_");
    assert_eq!(auto_spk_cache_stem("2060 Chiron"), "2060_chiron");
}

#[test]
fn horizons_url_shape_and_escaping() {
    let req = HorizonsSpkRequest::new("2060", "2000-01-01", "2100-01-01");
    assert_eq!(
        build_horizons_url(&req),
        "https://ssd.jpl.nasa.gov/api/horizons.api?format=json&COMMAND='2060'&OBJ_DATA=NO&MAKE_EPHEM=YES&EPHEM_TYPE=SPK&CENTER='500@0'&START_TIME='2000-01-01'&STOP_TIME='2100-01-01'"
    );
    // The semicolon bypass syntax must arrive encoded.
    let named = HorizonsSpkRequest::new("Ceres;", "2000-01-01", "2100-01-01");
    let url = build_horizons_url(&named);
    assert!(url.contains("COMMAND='Ceres%3B'"), "{url}");
    assert!(url.contains("CENTER='500@0'"), "{url}");
    // Spaces become `+`, never `%20`.
    let spaced = HorizonsSpkRequest::new("2060 Chiron", "2000-01-01", "2100-01-01");
    let url = build_horizons_url(&spaced);
    assert!(url.contains("COMMAND='2060+Chiron'"), "{url}");
}

#[test]
fn base64_standard_vectors_and_rejections() {
    assert_eq!(decode_base64_standard("Zg==").unwrap(), b"f");
    assert_eq!(decode_base64_standard("Zm9v").unwrap(), b"foo");
    assert_eq!(decode_base64_standard("Zm9vYmFy").unwrap(), b"foobar");
    assert_eq!(decode_base64_standard("QUJD").unwrap(), b"ABC");
    // Empty, bad length, bad alphabet, stray whitespace.
    assert!(decode_base64_standard("").is_none());
    assert!(decode_base64_standard("Zg=").is_none());
    assert!(decode_base64_standard("Zm9v\n").is_none());
    assert!(decode_base64_standard("Zm!v").is_none());
    // Padding only at the very end, at most two chars.
    assert!(decode_base64_standard("====").is_none());
    assert!(decode_base64_standard("Zm==Zm9v").is_none());
    assert!(decode_base64_standard("Zg==Zm9v").is_none());
    assert!(decode_base64_standard("Zm9vYmE=").unwrap() == b"fooba");
}

#[test]
fn horizons_base64_matches_ascii_discard_and_padding_rules() {
    // Primitive oracle: CPython 3.12.12 base64.b64decode, validate=False.
    // The pinned download_spk calls that default primitive directly.
    let accepted: &[(&str, &[u8])] = &[
        ("", b""),
        ("!!", b""),
        ("====", b""),
        ("===QUJD", b"ABC"),
        ("Zm9v\n", b"foo"),
        ("Q U!J\r\nD\t", b"ABC"),
        ("\0QUJD\x7f", b"ABC"),
        ("Zg===", b"f"),
        ("Zg==Zm9v", b"f"),
        ("Zm==Zm9v", b"f"),
        ("Z===g==", b"f"),
        ("Z=g==", b"f"),
        ("Zg=!\n=", b"f"),
        ("Zg=!=x", b"f"),
        ("Zm9=v", b"fo"),
        ("Zm9v=", b"foo"),
        ("Zm9v====", b"foo"),
        ("Zm9v=YmFy", b"foobar"),
        ("Zm9v==YmFy", b"foobar"),
        ("AA=AA", &[0, 0, 0]),
        ("AAA=A", &[0, 0]),
        ("AAAA=AAAA", &[0, 0, 0, 0, 0, 0]),
        ("AA===A", &[0]),
    ];
    for &(input, expected) in accepted {
        assert_eq!(
            super::decode::decode_horizons_base64(input).as_deref(),
            Some(expected),
            "{input:?}"
        );
    }
    for input in [
        "Z", "Zg", "Zg=", "Zg=Z", "A===", "A===AAAA", "éQUJD", "QUJDé", "Zg==é",
    ] {
        assert!(
            super::decode::decode_horizons_base64(input).is_none(),
            "{input:?}"
        );
    }
}

#[test]
fn horizons_response_envelope_branches() {
    // API refusal carries the server text.
    let err = parse_horizons_spk_response(r#"{"error": "unknown target"}"#, "Nope").unwrap_err();
    assert!(matches!(err, Error::HorizonsApiError { .. }));
    assert!(err.to_string().contains("unknown target"));
    // Payload decodes through the envelope.
    let payload = parse_horizons_spk_response(r#"{"spk_file_id": "7", "spk": "QUJD"}"#, "2060")
        .expect("payload");
    assert_eq!(payload, b"ABC");
    // Body-not-found answers arrive as `result` without a payload.
    let err = parse_horizons_spk_response(r#"{"result": "no such body"}"#, "Nope").unwrap_err();
    assert!(matches!(err, Error::HorizonsApiError { .. }));
    // Invalid padding and shapeless bodies are malformed, never silent.
    assert!(parse_horizons_spk_response(r#"{"spk_file_id": "7", "spk": "A"}"#, "2060").is_err());
    assert!(parse_horizons_spk_response(r#"{"hello": 1}"#, "2060").is_err());
    assert!(parse_horizons_spk_response("not json", "2060").is_err());
}

#[test]
fn horizons_response_accepts_wrapped_payload_before_separate_bsp_validation() {
    let payload = parse_horizons_spk_response(
        r#"{"spk_file_id":"fixture","spk":"\tREFG\r\nL1NQ!SyA=\n"}"#,
        "fixture",
    )
    .unwrap();
    assert_eq!(payload, SPK_FILE_ID);
    // Decoding accepts discarded input; the independent BSP structural
    // validator still refuses an empty file or a bare eight-byte header.
    let empty =
        parse_horizons_spk_response(r#"{"spk_file_id":"fixture","spk":"!!"}"#, "fixture").unwrap();
    assert!(empty.is_empty());
    assert!(FileRecord::parse(&empty).is_err());
    assert!(FileRecord::parse(&payload).is_err());
    assert!(
        parse_horizons_spk_response(r#"{"spk_file_id":"fixture","spk":"QUJD\u00e9"}"#, "fixture",)
            .is_err()
    );
}

#[test]
fn horizons_wrapped_synthetic_bsp_preserves_all_kernel_bytes() {
    // Independent direction: encode the existing hand-built DAF fixture,
    // with MIME-style wrapping, then consume the real response parser.
    const ALPHABET: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let image = build_image(DafEndian::Little, &[vec![one_seg()]]);
    let mut encoded = String::new();
    for bytes in image.chunks(3) {
        let word = (u32::from(bytes[0]) << 16)
            | (u32::from(bytes.get(1).copied().unwrap_or(0)) << 8)
            | u32::from(bytes.get(2).copied().unwrap_or(0));
        encoded.push(ALPHABET[(word >> 18) as usize] as char);
        encoded.push(ALPHABET[((word >> 12) & 63) as usize] as char);
        encoded.push(if bytes.len() > 1 {
            ALPHABET[((word >> 6) & 63) as usize] as char
        } else {
            '='
        });
        encoded.push(if bytes.len() > 2 {
            ALPHABET[(word & 63) as usize] as char
        } else {
            '='
        });
    }
    let wrapped = encoded
        .as_bytes()
        .chunks(76)
        .map(|line| std::str::from_utf8(line).unwrap())
        .collect::<Vec<_>>()
        .join(r"\r\n");
    let response = format!(r#"{{"spk_file_id":"fixture","spk":"{wrapped}"}}"#);
    let decoded = parse_horizons_spk_response(&response, "fixture").unwrap();
    assert_eq!(decoded, image);
    let index = SpkIndex::parse(&decoded).unwrap();
    assert_eq!(index.segments_for(499).count(), 1);
}

#[test]
fn horizons_download_keeps_valid_cached_kernel() {
    use crate::context::EffectiveNetworkPolicy;
    use crate::net::TransportConfig;
    let dir = std::env::temp_dir().join(format!("libeph-spkdl-{}", std::process::id()));
    let _ = std::fs::create_dir_all(&dir);
    let dest = dir.join("2060_2458849_2462502.bsp");
    std::fs::write(&dest, build_image(DafEndian::Little, &[vec![one_seg()]])).expect("seed kernel");
    // Sealed policy would refuse socket work, but the keep leg must
    // answer first without touching the network.
    let fetched = download_horizons_spk_to(
        "2060",
        "2020-01-01",
        "2030-01-01",
        "500@0",
        &dest,
        false,
        EffectiveNetworkPolicy::Sealed,
        &TransportConfig::default(),
    )
    .expect("keep leg");
    assert!(!fetched);
    assert!(dest.exists());
    let _ = std::fs::remove_file(&dest);
    // A corrupt cached file under a sealed policy surfaces the
    // refusal (never a fetch attempt, never a deletion).
    std::fs::write(&dest, b"truncated").expect("seed corrupt");
    let err = download_horizons_spk_to(
        "2060",
        "2020-01-01",
        "2030-01-01",
        "500@0",
        &dest,
        false,
        EffectiveNetworkPolicy::Sealed,
        &TransportConfig::default(),
    )
    .unwrap_err();
    assert!(matches!(err, Error::NetworkSealed { .. }));
    assert!(dest.exists());
    let _ = std::fs::remove_file(&dest);
}

#[test]
fn coverage_tolerance_bands() {
    assert!(spk_coverage_covers_request(100.0, 200.0, 100.0, 200.0));
    // Start may sit up to 1.5 d inside the request...
    assert!(spk_coverage_covers_request(101.5, 200.0, 100.0, 200.0));
    assert!(!spk_coverage_covers_request(
        101.5 + 1e-9,
        200.0,
        100.0,
        200.0
    ));
    // ...and the end up to 0.5 d short.
    assert!(spk_coverage_covers_request(100.0, 199.5, 100.0, 200.0));
    assert!(!spk_coverage_covers_request(
        100.0,
        199.5 - 1e-9,
        100.0,
        200.0
    ));
}

#[test]
fn type21_reach_table() {
    assert_eq!(type21_earth_reach_au(0), Some(1.1));
    assert_eq!(type21_earth_reach_au(10), Some(1.1));
    assert_eq!(type21_earth_reach_au(3), Some(2.2));
    assert_eq!(type21_earth_reach_au(399), Some(2.2));
    assert_eq!(type21_earth_reach_au(599), Some(5.5 + 1.1));
    assert_eq!(type21_earth_reach_au(-1), None);
    assert_eq!(type21_earth_reach_au(12345), None);
    assert!(type21_center_reach_ok(10));
    assert!(!type21_center_reach_ok(-7));
    // Usable-span arithmetic core: start loses light-time plus the
    // stencil, the end loses only the stencil.
    let (start, end) = type21_usable_span(100.0, 200.0, 5.0, 1.1, 0.01);
    assert!((start - (100.0 + 6.1 / LIGHT_SPEED_AU_PER_DAY + 0.01)).abs() < 1e-15);
    assert!((end - 199.99).abs() < 1e-15);
}

#[test]
fn usable_coverage_query_and_stencil_selectors() {
    // Request stencil: zero
    // without SPEED, the 1 s body step with SPEED, the 300 s geometric
    // asteroid step under TRUEPOS; the default is their maximum.
    assert_eq!(type21_request_stencil_days(0), 0.0);
    assert_eq!(
        type21_request_stencil_days(crate::constants::FLG_SPEED),
        1.0 / 86400.0
    );
    assert_eq!(
        type21_request_stencil_days(crate::constants::FLG_SPEED | crate::constants::FLG_TRUEPOS),
        300.0 / 86400.0
    );
    assert_eq!(type21_default_stencil_days(), 300.0 / 86400.0);
    // One linear EMDA segment, center Sun, one day of stored span: the
    // query extracts the stored span, resolves the reach, probes the
    // center-relative distance a micro-day inside the start and applies
    // the arithmetic core. The probe advances the linear record by
    // 0.0864 s at 30 km/s from [1e8, 0, 0] km.
    let recs = [linear_emda_rec(0.0, 4)];
    let img = emda_image(
        DafEndian::Little,
        20002060,
        10,
        0.0,
        86400.0,
        4,
        &recs,
        &[0.0],
        None,
        "CHIRON",
    );
    let mut registry = SpkRegistry::new();
    registry.register_image("c.bsp", img).expect("image loads");
    registry
        .register_body(15, "c.bsp", 20002060)
        .expect("registers");
    let probe_x_km = 100_000_000.0 + 30.0 * 0.0864;
    let distance_au = probe_x_km / AUNIT_TO_KM;
    let stencil = 300.0 / 86400.0;
    let expected_start = 2_451_545.0 + (distance_au + 1.1) / LIGHT_SPEED_AU_PER_DAY + stencil;
    let expected_end = 2_451_546.0 - stencil;
    let (start, end) = registry
        .usable_coverage_jd_tdb(15, stencil)
        .expect("usable span");
    assert!(
        (start - expected_start).abs() < 1e-9,
        "{start} vs {expected_start}"
    );
    assert!(
        (end - expected_end).abs() < 1e-12,
        "{end} vs {expected_end}"
    );
    // Zero stencil widens both ends by exactly the stencil.
    let (start0, end0) = registry
        .usable_coverage_jd_tdb(15, 0.0)
        .expect("zero-stencil span");
    assert!((start0 - (expected_start - stencil)).abs() < 1e-9);
    assert!((end0 - (expected_end + stencil)).abs() < 1e-12);
    // Unregistered bodies and non-type-21 registrations have no usable
    // span; callers fall back to the stored span.
    assert_eq!(registry.usable_coverage_jd_tdb(16, stencil), None);
    let mut legacy = SpkRegistry::new();
    legacy
        .register_image("t.bsp", const_image(8, 9, [1.0, 0.0, 0.0]))
        .expect("image loads");
    legacy.register_body(15, "t.bsp", 8).expect("registers");
    assert_eq!(legacy.usable_coverage_jd_tdb(15, stencil), None);
}

/// One day of linear heliocentric type-21 motion for provider-admission tests.
pub(crate) fn geometric_minor_test_image() -> Vec<u8> {
    emda_image(
        DafEndian::Little,
        20005145,
        10,
        0.0,
        86400.0,
        4,
        &[linear_emda_rec(0.0, 4)],
        &[0.0],
        None,
        "LINEAR MINOR",
    )
}

#[test]
fn jd_iso_round_trip() {
    // Docstring pins from `_iso_to_jd`.
    assert_eq!(iso_to_jd("2000-01-01"), Some(2_451_544.5));
    assert_eq!(iso_to_jd("1550-01-01"), Some(2_287_185.5));
    assert_eq!(jd_to_iso_date(2_451_544.5), "2000-01-01");
    assert_eq!(jd_to_iso_date(2_287_185.5), "1550-01-01");
    // Sub-day fractions select the containing civil date.
    assert_eq!(jd_to_iso_date(2_451_544.9), "2000-01-01");
    assert_eq!(jd_to_iso_date(2_451_545.5), "2000-01-02");
    // BCE round-trip through the leading-dash spelling.
    let ancient = jd_to_iso_date(1_000_000.0);
    assert!(ancient.starts_with('-'), "{ancient}");
    let back = iso_to_jd(&ancient).expect("BCE spelling parses");
    assert!((jd_to_iso_date(back) == ancient), "{ancient} vs {back}");
    // Malformed inputs are `None`, never a fabricated day.
    assert_eq!(iso_to_jd("2000-01"), None);
    assert_eq!(iso_to_jd("2000-01-01-02"), None);
    assert_eq!(iso_to_jd("2000-xx-01"), None);
    assert_eq!(iso_to_jd(""), None);
    // Day-granular cache filenames truncate toward zero.
    assert_eq!(
        generate_spk_cache_filename("2060", 2_458_849.7, 2_462_502.2),
        "2060_2458849_2462502.bsp"
    );
    assert_eq!(
        generate_spk_cache_filename("Ceres;", 2_458_849.0, 2_462_502.0),
        "ceres_2458849_2462502.bsp"
    );
}

#[test]
fn covering_lookup_prefix_and_gates() {
    let candidates = [
        ("2060_2458849_2462502.bsp", 100.0, 300.0),
        ("2060_2451000_2452000.bsp", 10.0, 20.0),
        ("ceres_2458849_2462502.bsp", 100.0, 300.0),
        ("2060_notes.txt", 100.0, 300.0),
    ];
    let ok = |_: &str| true;
    assert_eq!(
        find_covering_spk("2060", 150.0, 250.0, &candidates, ok),
        Some("2060_2458849_2462502.bsp")
    );
    // Wrong body prefix never matches.
    assert_eq!(find_covering_spk("10", 150.0, 250.0, &candidates, ok), None);
    // Out-of-range requests miss; the integrity gate can refuse.
    assert_eq!(
        find_covering_spk("2060", 150.0, 400.0, &candidates, ok),
        None
    );
    let refuse = |_: &str| false;
    assert_eq!(
        find_covering_spk("2060", 150.0, 250.0, &candidates, refuse),
        None
    );
}

#[test]
fn tier_windows_are_cumulative() {
    use crate::context::Tier;
    assert_eq!(
        spk_date_range_for_tier(Tier::Base),
        ("1850-01-01", "2150-01-01")
    );
    // Medium unions base: not its narrower literal.
    assert_eq!(
        spk_date_range_for_tier(Tier::Medium),
        ("1850-01-01", "2150-01-01")
    );
    assert_eq!(
        spk_date_range_for_tier(Tier::Extended),
        ("1600-01-01", "2500-01-01")
    );
    assert_eq!(
        spk_date_range_for_tier_name("extended").unwrap(),
        ("1600-01-01", "2500-01-01")
    );
    assert!(matches!(
        spk_date_range_for_tier_name("ultra"),
        Err(Error::InvalidTier { .. })
    ));
}

#[test]
fn auto_spk_policy_precedence() {
    // Programmatic wins over everything.
    assert!(!resolve_auto_spk_download(
        Some(false),
        Some("1"),
        Some(true)
    ));
    assert!(resolve_auto_spk_download(
        Some(true),
        Some("0"),
        Some(false)
    ));
    // Environment truth tables, case-insensitive and trimmed.
    assert!(!resolve_auto_spk_download(None, Some("0"), None));
    assert!(!resolve_auto_spk_download(None, Some(" False "), None));
    assert!(!resolve_auto_spk_download(None, Some("OFF"), None));
    assert!(resolve_auto_spk_download(None, Some("1"), None));
    assert!(resolve_auto_spk_download(None, Some(" yes "), None));
    assert!(resolve_auto_spk_download(None, Some("ENABLED"), None));
    // Unknown env values fall through to TOML, then the default.
    assert!(!resolve_auto_spk_download(None, Some("maybe"), Some(false)));
    assert!(resolve_auto_spk_download(None, Some("maybe"), None));
    assert!(resolve_auto_spk_download(None, None, None));
}

#[test]
fn major_asteroid_table_and_gates() {
    use crate::constants::{BENNU, CERES, CHIRON, ERIS, NAIF_CERES, NAIF_CHIRON};
    assert_eq!(
        get_major_asteroid_info(CERES),
        Some((1, "Ceres;", NAIF_CERES, "Ceres"))
    );
    assert_eq!(
        get_major_asteroid_info(CHIRON),
        Some((2060, "2060", NAIF_CHIRON, "Chiron"))
    );
    assert_eq!(get_major_asteroid_info(ERIS), None);
    assert!(is_major_asteroid(CERES));
    assert!(!is_major_asteroid(ERIS));
    assert!(is_spk_downloadable(CERES));
    assert!(is_spk_downloadable(CHIRON));
    // Bennu has map entries but Horizons blocks generation.
    assert!(!is_spk_downloadable(BENNU));
    assert!(!is_spk_downloadable(-999));
    // Backward-compatibility constant.
    assert!(!is_in_jpl_major_body_index(CERES));
}

#[test]
fn spk_diagnostic_names() {
    use crate::constants::{CERES, CHIRON, PHOLUS};
    assert_eq!(spk_body_name(CHIRON), Some("Chiron"));
    assert_eq!(spk_body_name(PHOLUS), Some("Pholus"));
    assert_eq!(spk_body_name(CERES), Some("Ceres"));
    assert_eq!(spk_body_name(20199), Some("Chariklo"));
    assert_eq!(spk_body_name(999_999), None);
}

#[test]
fn auto_spk_config_descriptor() {
    use crate::constants::CHIRON;
    let cfg = AutoSpkConfig::new(CHIRON, "2060");
    assert_eq!(cfg.ipl, CHIRON);
    assert_eq!(cfg.start, "2000-01-01");
    assert_eq!(cfg.end, "2100-01-01");
    assert_eq!(cfg.naif_id, None);
    assert!(cfg.enabled);
    assert_eq!(cfg.config_string(), "2060_2000-01-01_2100-01-01");
    assert_eq!(cfg.cache_filename_stem(), "2060");
    let named = AutoSpkConfig::new(17, "Ceres;");
    assert_eq!(named.cache_filename_stem(), "ceres_");
}

// API-07 U16h: tier download/ensure composition regression.
#[test]
fn u16h_download_filename_matches_release_convention() {
    assert_eq!(
        spk_download_filename("Chiron", "1850-01-01", "2150-01-01"),
        "chiron_185001_215001.bsp"
    );
    // The release `Ceres;` syntax sanitizes to `ceres` (shared
    // `sanitize_filename` kernel, swept under ).
    assert_eq!(
        spk_download_filename("Ceres;", "1600-01-01", "2500-01-01"),
        "ceres_160001_250001.bsp"
    );
    // Short labels keep their available characters.
    assert_eq!(
        spk_download_filename("2060", "2000", "21"),
        "2060_2000_21.bsp"
    );
}

#[test]
fn u16h_ensure_plan_is_pure_and_tier_gated() {
    use crate::constants::{SPK_AUTO_DOWNLOAD_BLOCKED, SPK_BODY_NAME_MAP};
    let plan = ensure_all_plan("base").expect("base tier plans");
    assert_eq!(plan.len(), SPK_BODY_NAME_MAP.len());
    let (start, end) = spk_date_range_for_tier_name("base").expect("base window");
    for row in &plan {
        let map = SPK_BODY_NAME_MAP
            .iter()
            .find(|(ipl, _, _)| *ipl == row.ipl)
            .expect("plan row maps");
        assert_eq!((row.horizons_id, row.naif_id), (map.1, map.2));
        assert!(row.filename.ends_with(".bsp"));
        assert!(row.filename.contains(&start.replace('-', "")[..6]));
        assert!(row.filename.contains(&end.replace('-', "")[..6]));
    }
    let blocked = plan.iter().filter(|row| row.blocked.is_some()).count();
    assert_eq!(blocked, SPK_AUTO_DOWNLOAD_BLOCKED.len());
    assert!(matches!(
        ensure_all_plan("nope"),
        Err(Error::InvalidTier { .. })
    ));
}

#[test]
fn u16h_request_labels_parse_to_ordered_epochs() {
    let a = request_label_to_et("1850-01-01").expect("base start parses");
    let b = request_label_to_et("2150-01-01").expect("base end parses");
    assert!(a.is_finite() && b.is_finite() && a < b);
    assert!(request_label_to_et("not-a-date").is_none());
    assert!(request_label_to_et("1850-01").is_none());
}

#[test]
fn u16h_ensure_sweep_reports_without_network() {
    use crate::net::TransportConfig;
    // Scratch dirs: empty data dir (tier kernel missing counts one
    // error, like the release `error:` arm) and cache dir.
    let root = std::env::temp_dir().join(format!("libephemeris-spk-u16h-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&root);
    let data_dir = root.join("data");
    let cache_dir = root.join("spk");
    std::fs::create_dir_all(&cache_dir).expect("cache dir");
    // Register one mapped body from a hand-built kernel so its row
    // reports `cached` with no socket work.
    let ipl = crate::constants::SPK_BODY_NAME_MAP[0].0;
    let naif = crate::constants::SPK_BODY_NAME_MAP[0].2;
    let img = discovery_image(naif, -1.0e12, 1.0e12, 1);
    let kernel_path = cache_dir.join("seed.bsp");
    std::fs::write(&kernel_path, &img).expect("seed kernel");
    let mut registry = SpkRegistry::new();
    let key = registry
        .register_file(&kernel_path)
        .expect("seed registers");
    registry
        .register_body(ipl, &key, naif)
        .expect("seed body registers");
    let config = TransportConfig::default();
    // Sealed policy: blocked + cached rows report without socket
    // work; the first downloadable row surfaces the refusal.
    let report = ensure_all_ephemerides_to(
        "base",
        &data_dir,
        &cache_dir,
        &mut registry,
        false,
        crate::context::EffectiveNetworkPolicy::Sealed,
        &config,
    );
    assert!(matches!(report, Err(Error::NetworkSealed { .. })));
    let _ = std::fs::remove_dir_all(&root);
}

#[test]
fn u16h_download_register_keeps_valid_file_without_network() {
    use crate::net::TransportConfig;
    let root =
        std::env::temp_dir().join(format!("libephemeris-spk-u16h-keep-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&root);
    let cache_dir = root.join("spk");
    std::fs::create_dir_all(&cache_dir).expect("cache dir");
    // Numeric body: the file-first chain matches the Horizons
    // convention id inside the kernel's own targets.
    let body = "2060";
    let (start, end) = ("1850-01-01", "2150-01-01");
    let req_start = request_label_to_et(start).expect("start parses");
    let req_end = request_label_to_et(end).expect("end parses");
    let img = discovery_image(20_002_060, req_start - 10.0, req_end + 10.0, 1);
    let dest = spk_download_dest(&cache_dir, body, start, end);
    std::fs::write(&dest, &img).expect("kept kernel");
    let mut registry = SpkRegistry::new();
    let config = TransportConfig::default();
    let out = download_and_register_spk_to(
        &mut registry,
        &cache_dir,
        body,
        2060,
        start,
        end,
        None,
        "500@0",
        false,
        crate::context::EffectiveNetworkPolicy::Sealed,
        &config,
    )
    .expect("kept file registers with no socket work");
    assert_eq!(out, dest.to_string_lossy());
    let info = registry.body_info(2060).expect("body registered");
    assert_eq!(info.naif_id, 20_002_060);
    let _ = std::fs::remove_dir_all(&root);
}

#[test]
fn cached_download_margins_use_days_and_consume_no_policy() {
    let root = std::env::temp_dir().join(format!("lephem-spk-day-margins-{}", std::process::id()));
    let cache_dir = root.join("spk");
    std::fs::create_dir_all(&cache_dir).unwrap();
    let body = "2060 Chiron";
    let (start, end) = ("2000-01-01", "2000-02-01");
    // The admitted cache contract permits 1.5 days at the start and 0.5 at
    // the end. These independently specified descriptor endpoints are within
    // those margins; interpreting the margins as seconds would refetch.
    let start_jd = crate::time::julday(2000, 1, 1, 0.0, crate::time::Calendar::Gregorian).unwrap();
    let end_jd = crate::time::julday(2000, 2, 1, 0.0, crate::time::Calendar::Gregorian).unwrap();
    let image = discovery_image(
        20_002_060,
        (start_jd + 1.0 - 2451545.0) * 86400.0,
        (end_jd - 0.25 - 2451545.0) * 86400.0,
        1,
    );
    let dest = spk_download_dest(&cache_dir, body, start, end);
    std::fs::write(&dest, image).unwrap();
    let mut registry = SpkRegistry::new();
    let result = download_and_register_spk_with_policy(
        &mut registry,
        &cache_dir,
        body,
        15,
        start,
        end,
        None,
        "500@0",
        false,
        &|| panic!("covering cache must not read policy"),
        &crate::net::TransportConfig::default(),
    )
    .unwrap();
    assert_eq!(result, dest.to_string_lossy());
    assert_eq!(registry.body_info(15).unwrap().naif_id, 20_002_060);
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn cached_coverage_refetch_precedes_missing_naif_deduction() {
    let root =
        std::env::temp_dir().join(format!("lephem-spk-refetch-order-{}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    let body = "unidentified";
    let (start, end) = ("1850-01-01", "2150-01-01");
    let dest = spk_download_dest(&root, body, start, end);
    std::fs::write(&dest, range_capture_chain_image(60.0, 60.0)).unwrap();
    let count = std::cell::Cell::new(0);
    let policy = || {
        count.set(count.get() + 1);
        Ok(crate::context::EffectiveNetworkPolicy::Sealed)
    };
    let error = download_and_register_spk_with_policy(
        &mut SpkRegistry::new(),
        &root,
        body,
        15,
        start,
        end,
        None,
        "500@0",
        false,
        &policy,
        &crate::net::TransportConfig::default(),
    )
    .unwrap_err();
    assert!(matches!(error, Error::NetworkSealed { .. }));
    assert_eq!(count.get(), 1);
    std::fs::remove_dir_all(root).unwrap();
}

/// Two readable scalar type2 links for acquisition-order and owner tests.
/// Both segments use one constant record; intervals are TDB seconds from J2000.
pub(crate) fn range_capture_chain_image(center_interval: f64, target_interval: f64) -> Vec<u8> {
    let rec = |value, interval| ChebRec {
        mid: interval / 2.0,
        radius: interval / 2.0,
        sets: vec![vec![value], vec![value + 1.0], vec![value + 2.0]],
    };
    let first = cheb_image(
        DafEndian::Little,
        2,
        3,
        0,
        0.0,
        center_interval,
        0.0,
        center_interval,
        &[rec(1.0, center_interval)],
        None,
        "CENTER",
    );
    let second = cheb_image(
        DafEndian::Little,
        2,
        399,
        3,
        0.0,
        target_interval,
        0.0,
        target_interval,
        &[rec(4.0, target_interval)],
        None,
        "TARGET",
    );
    let mut image = first;
    let summary = RECORD_BYTES;
    put_f64(&mut image, summary + 16, 2.0, DafEndian::Little);
    let first_word = 3 * WORDS_PER_RECORD + 1;
    let data_words = 9;
    let off = summary + 24 + SPK_NAME_BYTES;
    image[off..off + SPK_NAME_BYTES]
        .copy_from_slice(&second[summary + 24..summary + 24 + SPK_NAME_BYTES]);
    put_i32(
        &mut image,
        off + 32,
        (first_word + data_words) as i32,
        DafEndian::Little,
    );
    put_i32(
        &mut image,
        off + 36,
        (first_word + data_words * 2 - 1) as i32,
        DafEndian::Little,
    );
    image[2 * RECORD_BYTES + SPK_NAME_BYTES..2 * RECORD_BYTES + SPK_NAME_BYTES * 2]
        .copy_from_slice(&second[2 * RECORD_BYTES..2 * RECORD_BYTES + SPK_NAME_BYTES]);
    image[3 * RECORD_BYTES + data_words * 8..3 * RECORD_BYTES + data_words * 16]
        .copy_from_slice(&second[3 * RECORD_BYTES..3 * RECORD_BYTES + data_words * 8]);
    image
}

#[test]
fn registry_snapshots_share_images_and_survive_replacement_and_clear() {
    let original = const_image(10, 0, [1.0, 2.0, 3.0]);
    let mut registry = SpkRegistry::new();
    registry
        .register_image("shared.bsp", original.clone())
        .unwrap();
    registry.register_body(10, "shared.bsp", 10).unwrap();
    let retained = registry.clone();
    let held = &retained.kernels["shared.bsp"];
    let installed = &registry.kernels["shared.bsp"];
    assert!(Arc::ptr_eq(&held.bytes, &installed.bytes));
    assert!(Arc::ptr_eq(&held.index, &installed.index));
    assert!(Arc::ptr_eq(
        &retained.decoded["shared.bsp"],
        &registry.decoded["shared.bsp"]
    ));
    let before = retained.eval_body_relative(10, 0.0).unwrap();

    registry
        .register_image("shared.bsp", const_image(10, 0, [4.0, 5.0, 6.0]))
        .unwrap();
    assert!(!Arc::ptr_eq(
        &held.bytes,
        &registry.kernels["shared.bsp"].bytes
    ));
    assert!(!Arc::ptr_eq(
        &held.index,
        &registry.kernels["shared.bsp"].index
    ));
    assert_eq!(held.bytes(), original);
    assert_eq!(held.index().len(), 1);
    assert_eq!(
        retained.body_info(10),
        Some(BodyRegistration {
            body: 10,
            kernel: "shared.bsp".into(),
            naif_id: 10,
        })
    );
    registry.clear();
    assert_eq!(held.bytes(), original);
    assert_eq!(held.index().len(), 1);
    assert_eq!(retained.body_count(), 1);
    assert_eq!(retained.eval_body_relative(10, 0.0).unwrap(), before);
}

/// A small owned file for explicit-registration resource controls. Its unique
/// path prevents concurrent tests from publishing or replacing one another's
/// images, and cleanup never touches a caller-supplied asset.
#[cfg(unix)]
struct RegisteredSpkTestFile(std::path::PathBuf);

#[cfg(unix)]
impl RegisteredSpkTestFile {
    fn new(bytes: &[u8]) -> Self {
        use std::sync::atomic::{AtomicU64, Ordering};
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "lephem-explicit-spk-{}-{nonce}-{}.bsp",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::write(&path, bytes).unwrap();
        Self(path)
    }
}

#[cfg(unix)]
impl Drop for RegisteredSpkTestFile {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

/// Repeated explicit loads share only immutable storage. Each load still owns
/// a distinct actual descriptor/cursor, and mutable body maps stay independent.
#[cfg(unix)]
#[test]
fn explicit_file_registration_shares_storage_but_keeps_files_and_maps_independent() {
    use std::io::{Seek, SeekFrom};
    use std::os::fd::AsRawFd;

    let bytes = const_image(10, 0, [1.0, 2.0, 3.0]);
    let path = RegisteredSpkTestFile::new(&bytes);
    let mut first = SpkRegistry::new();
    let key = first.register_file(&path.0).unwrap();
    first.register_body(42, &key, 10).unwrap();
    let expected = first.eval_body_relative(42, 0.0).unwrap();
    let mut first_file = first
        .kernel(&key)
        .unwrap()
        .retained_source_file()
        .unwrap()
        .unwrap();
    first_file.seek(SeekFrom::Start(17)).unwrap();

    let mut second = SpkRegistry::new();
    assert_eq!(second.register_file(&path.0).unwrap(), key);
    let left = first.kernel(&key).unwrap();
    let right = second.kernel(&key).unwrap();
    assert!(Arc::ptr_eq(&left.bytes, &right.bytes));
    assert!(Arc::ptr_eq(&left.index, &right.index));
    assert!(Arc::ptr_eq(&first.decoded[&key], &second.decoded[&key]));
    assert!(!Arc::ptr_eq(
        left.source_file.as_ref().unwrap(),
        right.source_file.as_ref().unwrap()
    ));
    assert_ne!(
        left.source_file.as_ref().unwrap().as_raw_fd(),
        right.source_file.as_ref().unwrap().as_raw_fd()
    );
    let mut second_file = right.retained_source_file().unwrap().unwrap();
    assert_eq!(second_file.stream_position().unwrap(), bytes.len() as u64);
    assert_eq!(first_file.stream_position().unwrap(), 17);
    second_file.seek(SeekFrom::Start(29)).unwrap();
    assert_eq!(first_file.stream_position().unwrap(), 17);
    assert_eq!(second.body_count(), 0);
    second.register_body(43, &key, 10).unwrap();
    assert!(second.body_info(42).is_none());
    assert!(first.body_info(43).is_none());
    assert_eq!(second.eval_body_relative(43, 0.0).unwrap(), expected);

    // Explicit caller bytes remain a separate admission with new storage and
    // no file capability; they do not silently enter the file-hint channel.
    let retained = second.clone();
    second.register_image(&key, bytes).unwrap();
    assert!(!Arc::ptr_eq(
        &retained.kernel(&key).unwrap().bytes,
        &second.kernel(&key).unwrap().bytes
    ));
    assert!(second.kernel(&key).unwrap().source_file.is_none());
    assert_eq!(first.eval_body_relative(42, 0.0).unwrap(), expected);
}

/// A deliberately colliding private hint proves that stamps are never content
/// authority. A changed byte beyond the 64 KiB comparison buffer, followed by
/// changed coefficients, must produce newly decoded actual bytes even when the
/// hint names the current metadata and the old mtime has been restored.
#[cfg(unix)]
#[test]
fn explicit_file_registration_compares_all_bytes_even_with_matching_hint() {
    use super::file_image_cache::{self, FileImageKey};

    let mut original = const_image(10, 0, [1.0, 2.0, 3.0]);
    original.resize(128 * 1024, 0);
    let path = RegisteredSpkTestFile::new(&original);
    let modified = std::fs::metadata(&path.0).unwrap().modified().unwrap();
    let mut source = SpkRegistry::new();
    let key = source.register_file(&path.0).unwrap();
    source.register_body(10, &key, 10).unwrap();
    let expected = source.eval_body_relative(10, 0.0).unwrap();

    let mut changed = original.clone();
    *changed.last_mut().unwrap() = 7;
    std::fs::write(&path.0, &changed).unwrap();
    let file = std::fs::File::options().write(true).open(&path.0).unwrap();
    file.set_times(std::fs::FileTimes::new().set_modified(modified))
        .unwrap();
    let current_hint = FileImageKey::acquired(&key, &file).unwrap();
    // Force the exact metadata collision that a coarse or misleading host
    // could present, without changing the production comparison/admission path.
    file_image_cache::remember(current_hint, &source);
    drop(file);
    let mut destination = SpkRegistry::new();
    destination.register_file(&path.0).unwrap();
    assert_eq!(destination.kernel(&key).unwrap().bytes(), changed);
    assert!(!Arc::ptr_eq(
        &source.kernel(&key).unwrap().bytes,
        &destination.kernel(&key).unwrap().bytes
    ));
    destination.register_body(10, &key, 10).unwrap();
    assert_eq!(destination.eval_body_relative(10, 0.0).unwrap(), expected);

    let mut changed = const_image(10, 0, [4.0, 5.0, 6.0]);
    changed.resize(original.len(), 0);
    std::fs::write(&path.0, &changed).unwrap();
    let file = std::fs::File::options().write(true).open(&path.0).unwrap();
    file.set_times(std::fs::FileTimes::new().set_modified(modified))
        .unwrap();
    file_image_cache::remember(FileImageKey::acquired(&key, &file).unwrap(), &source);
    drop(file);
    destination.register_file(&path.0).unwrap();
    assert_eq!(destination.kernel(&key).unwrap().bytes(), changed);
    assert_eq!(
        destination
            .eval_body_relative(10, 0.0)
            .unwrap()
            .0
            .position_km,
        [4.0, 5.0, 6.0]
    );
    assert_eq!(source.eval_body_relative(10, 0.0).unwrap(), expected);
}

/// Replacing or removing the pathname must still acquire/refuse that pathname;
/// previously retained descriptors and snapshots continue to serve their bytes.
#[cfg(unix)]
#[test]
fn explicit_file_registration_replacement_removal_and_refusal_preserve_snapshots() {
    use std::io::{Read, Seek, SeekFrom};

    let original = const_image(10, 0, [1.0, 2.0, 3.0]);
    let changed = const_image(10, 0, [4.0, 5.0, 6.0]);
    let path = RegisteredSpkTestFile::new(&original);
    let mut current = SpkRegistry::new();
    let key = current.register_file(&path.0).unwrap();
    current.register_body(10, &key, 10).unwrap();
    let snapshot = current.clone();
    let replacement = RegisteredSpkTestFile::new(&changed);
    std::fs::rename(&replacement.0, &path.0).unwrap();
    current.register_file(&path.0).unwrap();
    assert!(!Arc::ptr_eq(
        &snapshot.kernel(&key).unwrap().bytes,
        &current.kernel(&key).unwrap().bytes
    ));
    assert_eq!(
        current.eval_body_relative(10, 0.0).unwrap().0.position_km,
        [4.0, 5.0, 6.0]
    );
    assert_eq!(
        snapshot.eval_body_relative(10, 0.0).unwrap().0.position_km,
        [1.0, 2.0, 3.0]
    );
    let mut old_file = snapshot
        .kernel(&key)
        .unwrap()
        .retained_source_file()
        .unwrap()
        .unwrap();
    old_file.seek(SeekFrom::Start(0)).unwrap();
    let mut acquired = Vec::new();
    old_file.read_to_end(&mut acquired).unwrap();
    assert_eq!(acquired, original);

    std::fs::remove_file(&path.0).unwrap();
    let before = current.clone();
    assert!(matches!(
        current.register_file(&path.0),
        Err(Error::SpkFileNotFound { .. })
    ));
    assert_eq!(current, before);
    std::fs::write(&path.0, [0u8; RECORD_BYTES]).unwrap();
    assert!(matches!(
        current.register_file(&path.0),
        Err(Error::SpkUnknownFormat { .. })
    ));
    assert_eq!(current, before);
}

/// The weak cache cannot root a descriptor, decoded cache, snapshot or the last
/// held buffer. A standalone range owner can reuse still-live immutable bytes
/// even after the old registry's coefficient cells have expired.
#[cfg(unix)]
#[test]
fn explicit_file_registration_last_buffer_and_snapshot_expire_without_cache_roots() {
    use super::file_image_cache::{self, FileImageKey};

    let bytes = const_image(10, 0, [1.0, 2.0, 3.0]);
    let path = RegisteredSpkTestFile::new(&bytes);
    let mut first = SpkRegistry::new();
    let key = first.register_file(&path.0).unwrap();
    let kernel = first.kernel(&key).unwrap();
    let hint = FileImageKey::acquired(&key, kernel.source_file.as_ref().unwrap()).unwrap();
    let held = RangeImage::acquired(kernel);
    let image = Arc::downgrade(&kernel.bytes);
    let index = Arc::downgrade(&kernel.index);
    let first_file = Arc::downgrade(kernel.source_file.as_ref().unwrap());
    let old_decoded = Arc::downgrade(&first.decoded[&key]);
    first.clear();
    assert!(old_decoded.upgrade().is_none());
    let mut second = SpkRegistry::new();
    second.register_file(&path.0).unwrap();
    assert!(Arc::ptr_eq(
        &held.retained_kernel().bytes,
        &second.kernel(&key).unwrap().bytes
    ));
    second.register_body(10, &key, 10).unwrap();
    let second_file = Arc::downgrade(second.kernel(&key).unwrap().source_file.as_ref().unwrap());
    let decoded = Arc::downgrade(&second.decoded[&key]);
    let snapshot = second.clone();
    second.clear();
    std::fs::remove_file(&path.0).unwrap();
    assert_eq!(
        snapshot.eval_body_relative(10, 0.0).unwrap().0.position_km,
        [1.0, 2.0, 3.0]
    );
    assert_eq!(held.bytes(), bytes);
    drop(snapshot);
    assert!(second_file.upgrade().is_none());
    assert!(decoded.upgrade().is_none());
    assert!(image.upgrade().is_some());
    assert!(first_file.upgrade().is_some());
    drop(held);
    assert!(image.upgrade().is_none());
    assert!(index.upgrade().is_none());
    assert!(first_file.upgrade().is_none());
    assert!(file_image_cache::retained(&hint).is_none());
}

/// A cache entry may remain after the final kernel, but it owns only a weak
/// small carrier. In particular, it adds no direct weak unsized-byte owner:
/// the sole weak byte/index fields live in the carrier and die with the last
/// kernel. No test-created weak byte handle can itself reserve the allocation.
#[cfg(unix)]
#[test]
fn explicit_file_registration_dead_hint_keeps_no_unsized_byte_weak_owner() {
    use super::file_image_cache::{FileImageCache, FileImageKey};

    let bytes = const_image(10, 0, [1.0, 2.0, 3.0]);
    let path = RegisteredSpkTestFile::new(&bytes);
    let key = path.0.to_string_lossy().into_owned();
    let mut source = SpkRegistry::new();
    source.retain_owned_kernel(
        SpkKernel::decode_file_image(&key, bytes.clone(), std::fs::File::open(&path.0).unwrap())
            .unwrap(),
    );
    let kernel = source.kernel(&key).unwrap();
    let hint = FileImageKey::acquired(&key, kernel.source_file.as_ref().unwrap()).unwrap();
    assert_eq!(Arc::strong_count(&kernel.bytes), 1);
    assert_eq!(Arc::weak_count(&kernel.bytes), 1);
    assert_eq!(Arc::weak_count(&kernel.index), 1);
    let owner = Arc::downgrade(kernel.file_image_owner.as_ref().unwrap());
    let mut cache = FileImageCache::default();
    cache.remember(hint.clone(), &source);
    assert_eq!(Arc::strong_count(&kernel.bytes), 1);
    assert_eq!(Arc::weak_count(&kernel.bytes), 1);
    assert_eq!(Arc::weak_count(&kernel.index), 1);
    source.clear();
    assert_eq!(cache.len(), 1);
    assert!(owner.upgrade().is_none());
    assert!(cache.retained(&hint).is_none());
    assert_eq!(cache.len(), 0);

    source.register_image(&key, bytes).unwrap();
    let image = source.kernel(&key).unwrap();
    assert!(image.file_image_owner.is_none());
    assert_eq!(Arc::weak_count(&image.bytes), 0);
}

/// Hints are bounded even while every image remains owned, and expired images
/// disappear at the next lookup. The local cache avoids global test ordering.
#[cfg(unix)]
#[test]
fn explicit_file_registration_weak_bookkeeping_is_bounded_and_expires() {
    use super::file_image_cache::{FileImageCache, FileImageKey};

    let bytes = const_image(10, 0, [1.0, 2.0, 3.0]);
    let path = RegisteredSpkTestFile::new(&bytes);
    let file = std::fs::File::open(&path.0).unwrap();
    let mut cache = FileImageCache::default();
    let mut owners = Vec::new();
    let mut keys = Vec::new();
    for ordinal in 0..129 {
        let key = format!("bounded-image-{ordinal}");
        let mut source = SpkRegistry::new();
        source.retain_owned_kernel(
            SpkKernel::decode_file_image(
                &key,
                bytes.clone(),
                std::fs::File::open(&path.0).unwrap(),
            )
            .unwrap(),
        );
        let hint = FileImageKey::acquired(&key, &file).unwrap();
        cache.remember(hint.clone(), &source);
        owners.push(source);
        keys.push(hint);
    }
    assert_eq!(cache.len(), 128);
    assert!(cache.retained(&keys[0]).is_none());
    assert!(cache.retained(&keys[128]).is_some());
    drop(owners);
    assert!(cache.retained(&keys[128]).is_none());
    assert_eq!(cache.len(), 0);
}

/// Observe physical allocations, not just the registry's logical removal.
/// A real file-backed constant segment is sampled before and after retirement;
/// weak handles never retain its bytes, index, source file or decoded cache.
#[test]
fn native_spk_final_lease_expires_actual_image_index_file_and_cache() {
    struct RemoveFile(std::path::PathBuf);
    impl Drop for RemoveFile {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }

    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path = RemoveFile(std::env::temp_dir().join(format!(
        "lephem-spk-final-lease-{}-{nonce}.bsp",
        std::process::id()
    )));
    let original = const_image(10, 0, [1.0, 2.0, 3.0]);
    std::fs::write(&path.0, &original).unwrap();
    let mut current = SpkRegistry::new();
    let key = current.register_file(&path.0).unwrap();
    current.register_body(10, &key, 10).unwrap();
    let before = current.eval_body_relative(10, 0.0).unwrap();
    let installed = current.kernel(&key).unwrap();
    let bytes = Arc::downgrade(&installed.bytes);
    let index = Arc::downgrade(&installed.index);
    let file = Arc::downgrade(installed.source_file.as_ref().unwrap());
    let decoded = Arc::downgrade(&current.decoded[&key]);
    let retained_kernel = installed.clone();
    let snapshot = Arc::new(current.clone());
    let active = Arc::clone(&snapshot);

    // Current publication replacement and removal cannot consume any lease
    // captured by the independent immutable snapshot or kernel clone.
    current
        .register_image(&key, const_image(10, 0, [4.0, 5.0, 6.0]))
        .unwrap();
    assert!(!Arc::ptr_eq(
        &retained_kernel.bytes,
        &current.kernel(&key).unwrap().bytes
    ));
    current.clear();
    assert!(current.is_empty());
    assert_eq!(active.eval_body_relative(10, 0.0).unwrap(), before);
    assert_eq!(retained_kernel.bytes(), original);
    drop(retained_kernel);
    drop(snapshot);
    assert_eq!(bytes.strong_count(), 1);
    assert_eq!(index.strong_count(), 1);
    assert_eq!(file.strong_count(), 1);
    assert_eq!(decoded.strong_count(), 1);

    // This is a genuine last kernel lease. The decoded registry cache has a
    // shorter lifetime and must expire independently when its owner retires.
    let final_kernel = active.kernel(&key).unwrap().clone();
    drop(active);
    assert!(decoded.upgrade().is_none());
    assert_eq!(bytes.strong_count(), 1);
    assert_eq!(index.strong_count(), 1);
    assert_eq!(file.strong_count(), 1);
    assert_eq!(final_kernel.bytes(), original);
    assert_eq!(final_kernel.index().len(), 1);
    assert_eq!(
        file.upgrade().unwrap().metadata().unwrap().len(),
        original.len() as u64
    );
    drop(final_kernel);
    assert!(bytes.upgrade().is_none());
    assert!(index.upgrade().is_none());
    assert!(file.upgrade().is_none());
    assert!(decoded.upgrade().is_none());
}

#[test]
fn native_directory_refusal_retains_the_actual_immutable_image() {
    let mut image = const_image(499, 0, [1.0; 3]);
    // The descriptor is wider than its stored directory: a real directory
    // refusal may retain that descriptor, but cannot extrapolate its record.
    put_f64(&mut image, RECORD_BYTES + 24, -128.0, DafEndian::Little);
    put_f64(&mut image, RECORD_BYTES + 32, 128.0, DafEndian::Little);
    let mut registry = SpkRegistry::new();
    registry
        .register_image("owned-range", image.clone())
        .unwrap();
    registry.register_body(499, "owned-range", 499).unwrap();
    let mut captured = None;
    assert!(matches!(
        registry.eval_geometric_ssb_au_two_part_jd_with_range(
            499,
            J2000_JD_TDB,
            100.0 / SECONDS_PER_DAY,
            &mut captured
        ),
        Err(Error::SpkOutOfRange { target: 499, .. })
    ));
    let raw = captured.take().unwrap();
    assert_eq!(raw.descriptor.target, 499);
    assert_eq!(raw.record_count, 1);
    assert_eq!(raw.image.bytes(), image);
    registry.clear();
    assert_eq!(raw.image.bytes(), image);
    assert!(
        registry
            .eval_geometric_ssb_au_two_part_jd_with_range(499, J2000_JD_TDB, 0.0, &mut captured)
            .is_err()
    );
    assert!(captured.is_none());
}

/// Digest reuse follows the same exact-file admission as immutable byte reuse.
/// A deliberately colliding hint cannot make changed bytes inherit a digest.
#[cfg(unix)]
#[test]
fn file_image_digest_reuses_only_fully_compared_bytes_and_retires_with_owner() {
    use super::file_image_cache::{self, FileImageKey};
    use sha2::{Digest, Sha256};

    let mut original = const_image(10, 0, [1.0, 2.0, 3.0]);
    original.resize(128 * 1024, 0);
    let path = RegisteredSpkTestFile::new(&original);
    let mut first = SpkRegistry::new();
    let key = first.register_file(&path.0).unwrap();
    let digest = first.kernel(&key).unwrap().image_sha256();
    let expected: [u8; 32] = Sha256::digest(&original).into();
    assert_eq!(digest, expected);
    let snapshot = first.clone();
    let owner = Arc::downgrade(
        first
            .kernel(&key)
            .unwrap()
            .file_image_owner
            .as_ref()
            .unwrap(),
    );
    let mut repeated = SpkRegistry::new();
    repeated.register_file(&path.0).unwrap();
    assert!(Arc::ptr_eq(
        first
            .kernel(&key)
            .unwrap()
            .file_image_owner
            .as_ref()
            .unwrap(),
        repeated
            .kernel(&key)
            .unwrap()
            .file_image_owner
            .as_ref()
            .unwrap()
    ));
    assert_eq!(repeated.kernel(&key).unwrap().image_sha256(), digest);

    // Change a byte beyond the comparison scratch length. Restore mtime and
    // force a current-metadata hint to exercise full comparison, not stamp luck.
    let modified = std::fs::metadata(&path.0).unwrap().modified().unwrap();
    let mut changed = original.clone();
    *changed.last_mut().unwrap() = 7;
    std::fs::write(&path.0, &changed).unwrap();
    let file = std::fs::File::options().write(true).open(&path.0).unwrap();
    file.set_times(std::fs::FileTimes::new().set_modified(modified))
        .unwrap();
    file_image_cache::remember(FileImageKey::acquired(&key, &file).unwrap(), &first);
    drop(file);
    repeated.register_file(&path.0).unwrap();
    let expected_changed: [u8; 32] = Sha256::digest(&changed).into();
    assert_ne!(expected_changed, digest);
    assert_eq!(
        repeated.kernel(&key).unwrap().image_sha256(),
        expected_changed
    );
    assert!(!Arc::ptr_eq(
        first
            .kernel(&key)
            .unwrap()
            .file_image_owner
            .as_ref()
            .unwrap(),
        repeated
            .kernel(&key)
            .unwrap()
            .file_image_owner
            .as_ref()
            .unwrap()
    ));
    assert_eq!(snapshot.kernel(&key).unwrap().image_sha256(), digest);
    first.clear();
    assert!(owner.upgrade().is_some());
    drop(snapshot);
    assert!(owner.upgrade().is_none());
}

/// Atomic replacement, absence and corrupt replacement still acquire/refuse
/// the current path before a digest can be reused; old snapshots stay valid.
#[cfg(unix)]
#[test]
fn file_image_digest_atomic_replacement_and_refusal_keep_old_snapshot_independent() {
    use sha2::{Digest, Sha256};

    let original = const_image(10, 0, [1.0, 2.0, 3.0]);
    let changed = const_image(10, 0, [4.0, 5.0, 6.0]);
    let path = RegisteredSpkTestFile::new(&original);
    let mut current = SpkRegistry::new();
    let key = current.register_file(&path.0).unwrap();
    let snapshot = current.clone();
    let old_digest = snapshot.kernel(&key).unwrap().image_sha256();
    let replacement = RegisteredSpkTestFile::new(&changed);
    std::fs::rename(&replacement.0, &path.0).unwrap();
    current.register_file(&path.0).unwrap();
    let expected: [u8; 32] = Sha256::digest(&changed).into();
    assert_eq!(current.kernel(&key).unwrap().image_sha256(), expected);
    assert_ne!(expected, old_digest);
    assert_eq!(snapshot.kernel(&key).unwrap().image_sha256(), old_digest);
    let accepted = current.clone();
    std::fs::remove_file(&path.0).unwrap();
    assert!(matches!(
        current.register_file(&path.0),
        Err(Error::SpkFileNotFound { .. })
    ));
    assert_eq!(current, accepted);
    std::fs::write(&path.0, [0u8; RECORD_BYTES]).unwrap();
    assert!(matches!(
        current.register_file(&path.0),
        Err(Error::SpkUnknownFormat { .. })
    ));
    assert_eq!(current, accepted);
    assert_eq!(current.kernel(&key).unwrap().image_sha256(), expected);
    assert_eq!(snapshot.kernel(&key).unwrap().image_sha256(), old_digest);
}
