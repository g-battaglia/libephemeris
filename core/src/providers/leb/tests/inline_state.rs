//! Inline geometric-state buffers retain record selection, words and refusals.

use super::*;
use crate::positions::context_serving::BoundProvider;
use crate::positions::{SingleLeb, StateProvider};
use std::sync::atomic::{AtomicU64, Ordering};

struct TemporaryImage(std::path::PathBuf);
impl TemporaryImage {
    fn new(bytes: &[u8]) -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "leb-inline-state-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::write(&path, bytes).unwrap();
        Self(path)
    }
}
impl Drop for TemporaryImage {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

fn words(values: &[f64]) -> Vec<u64> {
    values.iter().map(|v| v.to_bits()).collect()
}

fn check_provider(provider: &impl StateProvider, body: i32, epochs: &[TwoPartEpoch]) {
    for &epoch in epochs {
        let jd = epoch.value();
        let raw = provider.eval_raw(body, jd).unwrap();
        let cartesian = provider
            .serving_coord_type(body)
            .is_some_and(|coord| matches!(coord, COORD_ICRS_BARY | COORD_ICRS_BARY_SYSTEM));
        if let Some(inline) = provider.eval_raw_inline(body, jd).unwrap() {
            assert!(cartesian);
            assert_eq!(words(&inline.0), words(&raw.0));
            assert_eq!(words(&inline.1), words(&raw.1));
        } else {
            assert!(!cartesian);
        }
        let state = provider
            .geometric_state_with_offset(body, jd, |_| None)
            .unwrap();
        assert_eq!(words(&state.pos()), words(&raw.0));
        assert_eq!(words(&state.vel()), words(&raw.1));
        assert_eq!(state.tier(), provider.selected_tier(body, jd));
        assert!(!state.is_planet_center());
        let raw = provider.eval_raw_two_part(body, epoch).unwrap();
        if let Some(inline) = provider.eval_raw_inline_two_part(body, epoch).unwrap() {
            assert!(cartesian);
            assert_eq!(words(&inline.0), words(&raw.0));
            assert_eq!(words(&inline.1), words(&raw.1));
        } else {
            assert!(!cartesian);
        }
        let state = provider
            .geometric_state_two_part_with_offset(body, epoch, |_| None)
            .unwrap();
        assert_eq!(words(&state.pos()), words(&raw.0));
        assert_eq!(words(&state.vel()), words(&raw.1));
        assert_eq!(state.tier(), provider.selected_tier(body, jd));
    }
    for jd in [-1.0, 3.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let raw = provider.eval_raw(body, jd).unwrap_err();
        let inline = if provider
            .serving_coord_type(body)
            .is_some_and(|coord| matches!(coord, COORD_ICRS_BARY | COORD_ICRS_BARY_SYSTEM))
        {
            provider.eval_raw_inline(body, jd).unwrap_err()
        } else {
            assert!(provider.eval_raw_inline(body, jd).unwrap().is_none());
            provider.geometric_state(body, jd).unwrap_err()
        };
        assert_eq!(format!("{raw:?}"), format!("{inline:?}"));
        let epoch = TwoPartEpoch { hi: jd, lo: 0.0 };
        let raw = provider.eval_raw_two_part(body, epoch).unwrap_err();
        let inline = if provider
            .serving_coord_type(body)
            .is_some_and(|coord| matches!(coord, COORD_ICRS_BARY | COORD_ICRS_BARY_SYSTEM))
        {
            provider.eval_raw_inline_two_part(body, epoch).unwrap_err()
        } else {
            assert!(
                provider
                    .eval_raw_inline_two_part(body, epoch)
                    .unwrap()
                    .is_none()
            );
            provider.geometric_state_two_part(body, epoch).unwrap_err()
        };
        assert_eq!(format!("{raw:?}"), format!("{inline:?}"));
    }
    for jd in [f64::NAN, 0.5] {
        assert!(matches!(
            provider.eval_raw_inline(12345, jd),
            Err(Error::LebBodyNotFound { body_id: 12345 })
        ));
    }
}

#[test]
fn inline_state_formats_storage_composition_and_split_record_words() {
    let epochs = [
        TwoPartEpoch { hi: 0.0, lo: 0.0 },
        TwoPartEpoch {
            hi: 0.25,
            lo: 1e-18,
        },
        TwoPartEpoch {
            hi: 1.0_f64.next_down(),
            lo: 1e-18,
        },
        TwoPartEpoch {
            hi: 1.0,
            lo: -1e-18,
        },
        TwoPartEpoch { hi: 1.0, lo: 0.0 },
        TwoPartEpoch { hi: 1.0, lo: 1e-18 },
        TwoPartEpoch {
            hi: 1.0_f64.next_up(),
            lo: -1e-18,
        },
        TwoPartEpoch {
            hi: 1.75,
            lo: -1e-18,
        },
        TwoPartEpoch { hi: 2.0, lo: 0.0 },
    ];
    for image in [
        minimal_leb1(),
        leb2_v1_distinct_image(),
        leb2_v2_distinct_image(1.0),
        leb2_v2_distinct_image(1.0_f64.next_up()),
        leb2_v2_distinct_image(1.0_f64.next_down()),
    ] {
        for coord in [
            COORD_ICRS_BARY,
            COORD_ICRS_BARY_SYSTEM,
            COORD_ECLIPTIC,
            COORD_GEO_ECLIPTIC,
            COORD_HELIO_ECL_RETIRED,
        ] {
            let mut image = image.clone();
            let decoded = FileIndex::decode(&image).unwrap();
            let body_offset = decoded
                .sections()
                .iter()
                .find(|s| s.section_id() == SECTION_BODY_INDEX)
                .unwrap()
                .offset() as usize;
            image[body_offset + 4..body_offset + 8].copy_from_slice(&coord.to_le_bytes());
            let temporary = TemporaryImage::new(&image);
            let owned = LebFile::from_bytes("base_core.leb2".into(), image).unwrap();
            let indexed = LebFile::open_indexed(&temporary.0, IndexedLimits::default()).unwrap();
            for file in [&owned, &indexed] {
                check_provider(&SingleLeb::from_file(file).unwrap(), 0, &epochs);
                check_provider(&CompositeIndex::from_files([file]).unwrap(), 0, &epochs);
                check_provider(
                    &BoundProvider::Composite(CompositeIndex::from_files([file]).unwrap()),
                    0,
                    &epochs,
                );
                let tiered = TieredIndex::new(vec![(
                    LebTier::Medium,
                    CompositeIndex::from_files([file]).unwrap(),
                )])
                .unwrap();
                check_provider(&tiered, 0, &epochs);
            }
        }
    }
}

#[test]
fn inline_state_tier_and_composite_precedence_follow_actual_epoch() {
    let broad = leb1_const_image(0, 2.0, 0.0, 2.0);
    let narrow = leb1_const_image(0, 1.0, 0.5, 1.5);
    let first = LebFile::from_bytes("base_core.leb".into(), narrow).unwrap();
    let second = LebFile::from_bytes("medium_core.leb".into(), broad).unwrap();
    let tiered = TieredIndex::new(vec![
        (
            LebTier::Medium,
            CompositeIndex::from_files([&second]).unwrap(),
        ),
        (LebTier::Base, CompositeIndex::from_files([&first]).unwrap()),
    ])
    .unwrap();
    for (jd, value, tier) in [
        (0.0, 2.0, LebTier::Medium),
        (0.5, 1.0, LebTier::Base),
        (1.0, 1.0, LebTier::Base),
        (1.5, 1.0, LebTier::Base),
        (2.0, 2.0, LebTier::Medium),
    ] {
        let state = tiered.geometric_state(0, jd).unwrap();
        assert_eq!(state.pos(), [value; 3]);
        assert_eq!(state.tier(), Some(tier));
        assert_eq!(
            tiered
                .geometric_state_two_part(0, TwoPartEpoch { hi: jd, lo: 1e-18 })
                .unwrap()
                .pos(),
            [value; 3]
        );
    }
    let composite = CompositeIndex::from_files([&first, &second]).unwrap();
    assert_eq!(
        composite.eval_raw_inline(0, 1.0).unwrap().unwrap().0,
        [1.0; 3]
    );
    assert!(matches!(
        composite.eval_raw_inline(0, 0.0),
        Err(Error::LebOutOfRange { .. })
    ));
}

#[test]
fn inline_state_corrupt_payload_and_noncartesian_count_keep_existing_errors() {
    for mut image in [
        minimal_leb1(),
        leb2_v1_distinct_image(),
        leb2_v2_distinct_image(1.0),
    ] {
        let index = FileIndex::decode(&image).unwrap();
        let offset = index.find_body(0).unwrap().data_offset() as usize;
        if index.header().format() == LebFormat::Leb1 {
            image[offset..offset + 8].copy_from_slice(&f64::NAN.to_le_bytes());
        } else if index.header().format() == LebFormat::Leb2V1 {
            image[offset..offset + 4].fill(0);
        } else {
            // The directory remains structurally valid; its selected zstd frame fails.
            let blob_offset =
                u64::from_le_bytes(image[offset + 32..offset + 40].try_into().unwrap()) as usize;
            image[blob_offset..blob_offset + 4].fill(0);
        }
        let file = LebFile::from_bytes("corrupt".into(), image).unwrap();
        let provider = SingleLeb::from_file(&file).unwrap();
        let raw = provider.eval_raw(0, 0.25).unwrap_err();
        let inline = provider.eval_raw_inline(0, 0.25).unwrap_err();
        assert_eq!(format!("{raw:?}"), format!("{inline:?}"));
    }
    let offset = HEADER_SIZE + SECTION_DIR_SIZE + LEB1_BODY_ENTRY_SIZE;
    let mut image = pack_header(&MAGIC_LEB1, 1, 1, 1, 0.0, 2.0, 0.0, 0);
    image.extend(pack_section(
        0,
        (HEADER_SIZE + SECTION_DIR_SIZE) as u64,
        LEB1_BODY_ENTRY_SIZE as u64,
    ));
    image.extend(pack_body_leb1(
        0,
        COORD_ICRS_BARY,
        1,
        0.0,
        2.0,
        2.0,
        0,
        2,
        offset as u64,
    ));
    image.extend(1.0_f64.to_le_bytes());
    image.extend(2.0_f64.to_le_bytes());
    let file = LebFile::from_bytes("two-components".into(), image).unwrap();
    let provider = SingleLeb::from_file(&file).unwrap();
    assert!(provider.eval_raw_inline(0, 1.0).unwrap().is_none());
    assert!(
        matches!(provider.geometric_state(0, 1.0), Err(Error::LebCorrupt { what, offset: 0 }) if what == "LEB stored state (body 0)")
    );
    assert!(matches!(
        provider.geometric_state(0, f64::NAN),
        Err(Error::NonFiniteArgument { .. })
    ));
}

/// An arbitrary provider keeps its historical raw-route ordering by default.
struct RecordingProvider<'a> {
    inner: SingleLeb<'a>,
    calls: std::cell::RefCell<Vec<&'static str>>,
    malformed_shape: bool,
    raw_failure: bool,
}
impl StateProvider for RecordingProvider<'_> {
    fn provider_name(&self) -> &'static str {
        "test/recording-vector"
    }
    fn has_body(&self, body: i32) -> bool {
        self.inner.has_body(body)
    }
    fn body_coverage(&self, body: i32) -> Option<(f64, f64)> {
        self.inner.body_coverage(body)
    }
    fn serving_coord_type(&self, body: i32) -> Option<u32> {
        self.inner.serving_coord_type(body)
    }
    fn serving_entry_snapshot(&self, body: i32) -> Result<crate::positions::EntrySnapshot> {
        self.calls.borrow_mut().push("entry");
        self.inner.serving_entry_snapshot(body)
    }
    fn eval_raw(&self, body: i32, jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        self.calls.borrow_mut().push("raw");
        if self.raw_failure {
            return Err(Error::LebCorrupt {
                what: "recorded raw refusal".into(),
                offset: 73,
            });
        }
        let (mut pos, vel) = self.inner.eval_raw(body, jd)?;
        if self.malformed_shape {
            pos.pop();
        }
        Ok((pos, vel))
    }
    fn eval_raw_two_part(&self, body: i32, epoch: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)> {
        self.calls.borrow_mut().push("raw-split");
        self.inner.eval_raw_two_part(body, epoch)
    }
    fn selected_tier(&self, _: i32, _: f64) -> Option<LebTier> {
        self.calls.borrow_mut().push("tier");
        Some(LebTier::Medium)
    }
}

#[test]
fn inline_state_default_provider_keeps_call_error_and_center_order() {
    let mut bytes = minimal_leb1();
    bytes[HEADER_SIZE + SECTION_DIR_SIZE + 4..HEADER_SIZE + SECTION_DIR_SIZE + 8]
        .copy_from_slice(&COORD_ICRS_BARY_SYSTEM.to_le_bytes());
    let file = LebFile::from_bytes("system-center".into(), bytes).unwrap();
    let mut provider = RecordingProvider {
        inner: SingleLeb::from_file(&file).unwrap(),
        calls: std::cell::RefCell::new(Vec::new()),
        malformed_shape: false,
        raw_failure: false,
    };
    assert!(provider.eval_raw_inline(0, 1.0).unwrap().is_none());
    assert!(
        provider
            .eval_raw_inline_two_part(0, TwoPartEpoch { hi: 1.0, lo: 0.0 })
            .unwrap()
            .is_none()
    );
    assert!(provider.calls.borrow().is_empty());
    let raw = provider.inner.eval_raw(0, 1.0).unwrap();
    let offset = |jd: f64| Some([0.25 + 0.125 * jd, -0.5, 1.0]);
    let expected = eval_body_center(
        raw.0.try_into().unwrap(),
        raw.1.try_into().unwrap(),
        true,
        1.0,
        offset,
    );
    let state = provider
        .geometric_state_with_offset(0, 1.0, |jd| {
            provider.calls.borrow_mut().push("offset");
            offset(jd)
        })
        .unwrap();
    assert_eq!(words(&state.pos()), words(&expected.0));
    assert_eq!(words(&state.vel()), words(&expected.1));
    assert_eq!(
        *provider.calls.borrow(),
        ["entry", "raw", "tier", "offset", "offset", "offset"]
    );
    assert!(state.is_planet_center());
    provider.calls.borrow_mut().clear();
    let state = provider
        .geometric_state_two_part_with_offset(0, TwoPartEpoch { hi: 1.0, lo: 0.0 }, offset)
        .unwrap();
    assert_eq!(words(&state.pos()), words(&expected.0));
    assert_eq!(words(&state.vel()), words(&expected.1));
    assert_eq!(*provider.calls.borrow(), ["entry", "raw-split", "tier"]);
    provider.calls.borrow_mut().clear();
    provider.malformed_shape = true;
    assert!(matches!(
        provider.geometric_state_with_offset(0, 1.0, |_| panic!(
            "invalid shape cannot consume offsets"
        )),
        Err(Error::LebCorrupt { offset: 0, .. })
    ));
    assert_eq!(*provider.calls.borrow(), ["entry", "raw", "tier"]);
    provider.calls.borrow_mut().clear();
    provider.raw_failure = true;
    assert!(matches!(
        provider.geometric_state_with_offset(0, 1.0, |_| None),
        Err(Error::LebCorrupt { offset: 73, .. })
    ));
    assert_eq!(*provider.calls.borrow(), ["entry", "raw"]);
}

fn check_scalar_split_reuse(provider: &impl StateProvider, epochs: &[f64]) {
    for &jd in epochs {
        let epoch = TwoPartEpoch::split(jd).unwrap();
        assert!(provider.can_reuse_scalar_state(0, jd, epoch));
        let scalar = provider.geometric_state(0, jd).unwrap();
        let split = provider.geometric_state_two_part(0, epoch).unwrap();
        assert_eq!(words(&scalar.pos()), words(&split.pos()));
        assert_eq!(words(&scalar.vel()), words(&split.vel()));
        assert_eq!(scalar.body_id(), split.body_id());
        assert_eq!(scalar.jd_tt().to_bits(), split.jd_tt().to_bits());
        assert_eq!(scalar.frame(), split.frame());
        assert_eq!(scalar.tier(), split.tier());
        assert_eq!(scalar.is_planet_center(), split.is_planet_center());
        assert!(!provider.can_reuse_scalar_state(0, jd, TwoPartEpoch { hi: jd, lo: 1e-20 }));
        assert!(!provider.can_reuse_scalar_state(12345, jd, epoch));
    }
}

#[test]
fn scalar_split_reuse_preserves_whole_states_across_formats_readers_and_budgets() {
    let epochs = [
        0.0,
        0.25,
        0.5,
        1.0_f64.next_down(),
        1.0,
        1.0_f64.next_up(),
        1.75,
        2.0,
    ];
    for image in [
        minimal_leb1(),
        leb2_v1_distinct_image(),
        leb2_v2_distinct_image(1.0),
        leb2_v2_distinct_image(1.0_f64.next_up()),
        leb2_v2_distinct_image(1.0_f64.next_down()),
    ] {
        for coord in [COORD_ICRS_BARY, COORD_ICRS_BARY_SYSTEM] {
            let mut image = image.clone();
            let body_offset = FileIndex::decode(&image)
                .unwrap()
                .sections()
                .iter()
                .find(|s| s.section_id() == SECTION_BODY_INDEX)
                .unwrap()
                .offset() as usize;
            image[body_offset + 4..body_offset + 8].copy_from_slice(&coord.to_le_bytes());
            let temporary = TemporaryImage::new(&image);
            for cache_bytes in [0, 4096] {
                let owned = LebFile::from_bytes_with_cache_budget(
                    "base_core.leb2".into(),
                    image.clone(),
                    cache_bytes,
                )
                .unwrap();
                let indexed = LebFile::open_indexed(
                    &temporary.0,
                    IndexedLimits {
                        metadata_bytes: 4096,
                        cache_bytes,
                        // The synthetic encoder omits content size and advertises
                        // an 8 MiB window even for these small coefficient rows.
                        working_bytes: 9 * 1024 * 1024,
                        compressed_bytes: 4096,
                        frame_bytes: 4096,
                        window_log_max: 23,
                    },
                )
                .unwrap();
                for file in [&owned, &indexed] {
                    check_scalar_split_reuse(&SingleLeb::from_file(file).unwrap(), &epochs);
                    check_scalar_split_reuse(&CompositeIndex::from_files([file]).unwrap(), &epochs);
                    check_scalar_split_reuse(
                        &BoundProvider::Composite(CompositeIndex::from_files([file]).unwrap()),
                        &epochs,
                    );
                    check_scalar_split_reuse(
                        &TieredIndex::new(vec![(
                            LebTier::Medium,
                            CompositeIndex::from_files([file]).unwrap(),
                        )])
                        .unwrap(),
                        &epochs,
                    );
                    if cache_bytes == 0 {
                        assert_eq!(file.chunk_cache().unwrap().stats().2, 0);
                    }
                }
            }
        }
    }
}

#[test]
fn scalar_split_reuse_selects_actual_attach_order_and_date_tier() {
    let broad =
        LebFile::from_bytes("medium_core.leb".into(), leb1_const_image(0, 2.0, 0.0, 2.0)).unwrap();
    let narrow =
        LebFile::from_bytes("base_core.leb".into(), leb1_const_image(0, 1.0, 0.5, 1.5)).unwrap();
    let tiered = TieredIndex::new(vec![
        (
            LebTier::Medium,
            CompositeIndex::from_files([&broad]).unwrap(),
        ),
        (
            LebTier::Base,
            CompositeIndex::from_files([&narrow]).unwrap(),
        ),
    ])
    .unwrap();
    check_scalar_split_reuse(&tiered, &[0.0, 0.5, 1.0, 1.5, 2.0]);
    let composite = CompositeIndex::from_files([&narrow, &broad]).unwrap();
    check_scalar_split_reuse(&composite, &[0.5, 1.0, 1.5]);
    // Attach order owns the body even where a later image would cover TT.
    assert!(!composite.can_reuse_scalar_state(0, 0.0, TwoPartEpoch::split(0.0).unwrap()));
    assert!(matches!(
        composite.geometric_state_two_part(0, TwoPartEpoch::split(0.0).unwrap()),
        Err(Error::LebOutOfRange { .. })
    ));
}
