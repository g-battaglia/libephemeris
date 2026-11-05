//! Shared-snapshot and retained-cache memory regressions using synthetic data.
use super::*;
use crate::context::{EphemerisContext, ProviderMode, Tier};
use crate::positions::{ServingRequest, TierPolicy, with_shared_calc_engine};
use crate::time::DeltaTModel;

#[test]
fn retained_dispatch_shares_images_and_releases_its_last_owner() {
    let mut file =
        LebFile::from_bytes("base_core.leb2".into(), leb2_v2_distinct_image(1.0)).unwrap();
    let image = Arc::downgrade(file.image.as_ref().unwrap());
    let mut context = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
    context.attach_leb_file(file.clone());
    let index = context
        .leb_composite(crate::positions::TierPolicy::AttachOrder, false)
        .unwrap();
    let expected = index.eval_body_any(0, 0.25).unwrap();
    assert!(std::ptr::eq(index.entries.get(0).0, file.index().unwrap()));
    assert!(Arc::ptr_eq(&index.caches[0], file.chunk_cache().unwrap()));
    let mut other = context.clone();
    assert!(std::ptr::eq(
        index,
        other
            .leb_composite(crate::positions::TierPolicy::AttachOrder, false)
            .unwrap()
    ));
    file.close();
    context.close();
    assert!(image.upgrade().is_some());
    assert_eq!(
        other
            .leb_composite(crate::positions::TierPolicy::AttachOrder, false)
            .unwrap()
            .eval_body_any(0, 0.25)
            .unwrap(),
        expected
    );
    other.close();
    assert!(image.upgrade().is_none());
}

#[test]
fn retained_dispatch_rebuilds_after_relevant_edits_and_preserves_strict_errors() {
    use crate::positions::TierPolicy;
    let first = LebFile::from_bytes("base_core.leb2".into(), leb2_v2_distinct_image(1.0)).unwrap();
    let second =
        LebFile::from_bytes("base_core.leb2".into(), leb2_v2_scaled_image(1.0, 2.0)).unwrap();
    let mut context = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
    context.attach_leb_file(first);
    for edit in 0..4 {
        let snapshot = context.clone();
        let before = snapshot
            .leb_composite(TierPolicy::AttachOrder, false)
            .unwrap();
        match edit {
            0 => context.set_tier_name("medium").unwrap(),
            1 => context.set_mode_name("auto").unwrap(),
            2 => context.spk_registry_mut().clear(),
            _ => context.attach_leb_file(second.clone()),
        }
        let after = context
            .leb_composite(TierPolicy::AttachOrder, false)
            .unwrap();
        assert!(!std::ptr::eq(before, after));
        assert_eq!(
            before.eval_body_any(0, 0.25).unwrap(),
            after.eval_body_any(0, 0.25).unwrap()
        );
    }
    let snapshot = context.clone();
    let before = snapshot.leb_provider(false).unwrap();
    let mut closed = second.clone();
    closed.close();
    context.attach_leb_file(closed);
    assert!(context.leb_provider(true).is_ok());
    assert!(matches!(
        context.leb_provider(false),
        Err(Error::LebReaderClosed)
    ));
    assert!(context.leb_composite(TierPolicy::AttachOrder, true).is_ok());
    assert!(matches!(
        context.leb_composite(TierPolicy::AttachOrder, false),
        Err(Error::LebReaderClosed)
    ));
    context.release_leb_images();
    assert!(context.leb_provider(false).is_err());
    context.attach_leb_file(second);
    assert!(context.leb_provider(false).is_ok());
    assert!(crate::positions::StateProvider::has_body(before, 0));
}

#[test]
fn retained_dispatch_keeps_initial_tier_validation_order() {
    let unknown = LebFile::from_bytes("custom.leb2".into(), leb2_v2_distinct_image(1.0)).unwrap();
    let mut closed = unknown.clone();
    closed.close();
    let mut context = EphemerisContext::new();
    context.attach_leb_file(unknown);
    context.attach_leb_file(closed);
    assert!(matches!(
        context.leb_composite(crate::positions::TierPolicy::TierPriority, false),
        Err(Error::ConfigurationMissing {
            missing_config: "tier_priority_dispatch",
            ..
        })
    ));
}

#[test]
fn clone_shares_allocation_close_releases_last_owner() {
    fn send_sync<T: Send + Sync>() {}
    send_sync::<LebFile>();
    send_sync::<EphemerisContext>();
    let bytes = leb2_v2_distinct_image(1.0);
    let pointer = bytes.as_ptr();
    let mut file = LebFile::from_bytes("base_core.leb2".into(), bytes).unwrap();
    assert_eq!(file.bytes().unwrap().as_ptr(), pointer);
    let weak = Arc::downgrade(file.image.as_ref().unwrap());
    let clone = file.clone();
    assert_eq!(clone.bytes().unwrap().as_ptr(), pointer);
    assert!(std::ptr::eq(file.index().unwrap(), clone.index().unwrap()));
    assert_eq!(file.retained_chunk_index_bytes(), 0);
    file.eval_body_any(0, 0.25).unwrap();
    let metadata = file.retained_chunk_index_bytes();
    assert!(metadata > 0);
    assert_eq!(clone.retained_chunk_index_bytes(), metadata);
    assert_eq!(file, clone); // cache mutation is not semantic state
    let mut ctx = EphemerisContext::new();
    ctx.attach_leb_file(clone);
    let mut other = ctx.clone();
    assert_eq!(other.leb_files()[0].bytes().unwrap().as_ptr(), pointer);
    file.close();
    file.close();
    assert!(file.index().is_err());
    assert_eq!(file.retained_chunk_index_bytes(), 0);
    ctx.reset_session();
    assert_eq!(ctx.leb_files()[0].retained_chunk_index_bytes(), metadata);
    assert_eq!(ctx.leb_files()[0].chunk_cache_stats().unwrap(), (0, 1, 1));
    ctx.close();
    assert!(weak.upgrade().is_some());
    other.close();
    assert!(weak.upgrade().is_none());
}

#[test]
fn persistent_cache_survives_provider_rebuild_and_keeps_images_isolated() {
    let first = LebFile::from_bytes("base_core.leb2".into(), leb2_v2_distinct_image(1.0)).unwrap();
    let second =
        LebFile::from_bytes("base_core.leb2".into(), leb2_v2_scaled_image(1.0, 2.0)).unwrap();
    let direct = first
        .index()
        .unwrap()
        .eval_leb2_body(first.bytes().unwrap(), 0, 0.25)
        .unwrap();
    for _ in 0..3 {
        let composite = CompositeIndex::from_files([&first]).unwrap();
        assert!(std::ptr::eq(
            composite.entries.get(0).0,
            first.index().unwrap()
        ));
        assert!(Arc::ptr_eq(
            &composite.caches[0],
            first.chunk_cache().unwrap()
        ));
        assert_eq!(composite.eval_body_any(0, 0.25).unwrap(), direct);
    }
    assert_eq!(first.chunk_cache_stats().unwrap(), (2, 1, 1));
    let metadata = first.retained_chunk_index_bytes();
    assert!(metadata > 0);
    assert_eq!(second.retained_chunk_index_bytes(), 0);
    let single = crate::positions::SingleLeb::from_file(&first).unwrap();
    assert_eq!(
        crate::positions::StateProvider::eval_raw(&single, 0, 0.25).unwrap(),
        direct
    );
    assert_eq!(first.chunk_cache_stats().unwrap(), (3, 1, 1));
    assert_eq!(second.chunk_cache_stats().unwrap(), (0, 0, 0));
    assert_ne!(second.eval_body_any(0, 0.25).unwrap(), direct);
    let epoch = TwoPartEpoch::split(0.25).unwrap();
    assert_eq!(first.eval_body_any_two_part(0, epoch).unwrap(), direct);
    std::thread::scope(|scope| {
        for _ in 0..8 {
            let first = &first;
            let direct = &direct;
            scope.spawn(move || {
                for _ in 0..20 {
                    assert_eq!(&first.eval_body_any(0, 0.25).unwrap(), direct);
                }
            });
        }
    });
    let (resident, limit) = first.chunk_cache_memory().unwrap();
    assert!(resident > 0 && resident <= limit);
    first.cool();
    first.cool();
    assert_eq!(first.chunk_cache_memory().unwrap().0, 0);
    assert_eq!(first.retained_chunk_index_bytes(), metadata);
    assert_eq!(first.eval_body_any(0, 0.25).unwrap(), direct);
}

#[test]
fn byte_budget_eviction_oversize_and_disabled_retention_preserve_words() {
    let bytes = leb2_v2_distinct_image(1.0);
    let index = FileIndex::decode(&bytes).unwrap();
    let chunk_bytes = 48; // fixture: 1 segment × 3 components × 2 words × 8 bytes
    for budget in [0, 1, chunk_bytes, DEFAULT_CHUNK_CACHE_BYTES] {
        let file =
            LebFile::from_bytes_with_cache_budget("base_core.leb2".into(), bytes.clone(), budget)
                .unwrap();
        for jd in [
            0.0,
            0.25,
            1.0_f64.next_down(),
            1.0,
            1.0_f64.next_up(),
            1.5,
            2.0,
            0.25,
        ] {
            let expected = index.eval_leb2_body(&bytes, 0, jd).unwrap();
            let got = file.eval_body_any(0, jd).unwrap();
            assert_eq!(
                got.0.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                expected.0.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
            );
            assert_eq!(
                got.1.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                expected.1.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
            );
            assert!(file.chunk_cache_memory().unwrap().0 <= budget);
        }
        if budget < chunk_bytes {
            assert_eq!(file.chunk_cache_stats().unwrap().2, 0);
        }
        if budget == chunk_bytes {
            assert_eq!(file.chunk_cache_stats().unwrap().2, 1);
        }
        let metadata = file.retained_chunk_index_bytes();
        assert!(metadata > 0);
        assert_eq!(file.chunk_cache().unwrap().chunk_index_test_stats(), (1, 1));
        let before = file.chunk_cache_stats().unwrap();
        assert!(file.eval_body_any(0, 5.0).is_err());
        assert_eq!(file.chunk_cache_stats().unwrap(), before);
    }
}

#[test]
fn concurrent_first_directory_parse_is_shared_and_capacity_is_counted() {
    let bytes = leb2_v2_distinct_image(1.0);
    let index = FileIndex::decode(&bytes).unwrap();
    let (_, directory) =
        decode_chunk_index(&bytes, index.find_body(0).unwrap().data_offset(), 0).unwrap();
    let expected_bytes = directory.capacity() * std::mem::size_of::<ChunkEntry>();
    // Accounting follows actual capacity, whether parsing reserves exactly
    // the validated count or the allocator retains spare entries.
    assert!(directory.capacity() >= directory.len());
    let file = LebFile::from_bytes("base_core.leb2".into(), bytes).unwrap();
    let barrier = std::sync::Barrier::new(8);
    std::thread::scope(|scope| {
        for _ in 0..8 {
            let file = &file;
            let barrier = &barrier;
            scope.spawn(move || {
                barrier.wait();
                for jd in [0.25, 1.25] {
                    let expected = file
                        .index()
                        .unwrap()
                        .eval_leb2_body(file.bytes().unwrap(), 0, jd)
                        .unwrap();
                    assert_eq!(file.eval_body_any(0, jd).unwrap(), expected);
                }
            });
        }
    });
    assert_eq!(file.chunk_cache().unwrap().chunk_index_test_stats(), (1, 1));
    assert_eq!(file.retained_chunk_index_bytes(), expected_bytes);
    let debug = format!("{:?}", file.chunk_cache().unwrap());
    assert!(debug.contains("chunk_index_bytes"));
    assert!(!debug.contains("ChunkEntry") && !debug.contains("entries:"));
}

#[test]
fn directory_errors_remain_lazy_and_failed_parses_are_not_retained() {
    for corruption in [0, 1, 2] {
        let mut bytes = leb2_v2_distinct_image(1.0);
        let index = FileIndex::decode(&bytes).unwrap();
        let entry = index.find_body(0).unwrap();
        let offset = entry.data_offset() as usize;
        match corruption {
            0 => bytes.truncate(offset + CHUNK_INDEX_HEADER_SIZE - 1),
            1 => bytes[offset..offset + 4].copy_from_slice(&u32::MAX.to_le_bytes()),
            _ => bytes[offset..offset + 4].copy_from_slice(&0_u32.to_le_bytes()),
        }
        let file = LebFile::from_bytes("base_core.leb2".into(), bytes).unwrap();
        for _ in 0..2 {
            let expected = format!(
                "{:?}",
                file.index()
                    .unwrap()
                    .eval_leb2_body(file.bytes().unwrap(), 0, 0.25)
                    .unwrap_err()
            );
            let actual = format!("{:?}", file.eval_body_any(0, 0.25).unwrap_err());
            assert_eq!(actual, expected);
            assert_eq!(file.retained_chunk_index_bytes(), 0);
            assert_eq!(file.chunk_cache_stats().unwrap(), (0, 0, 0));
        }
        assert_eq!(file.chunk_cache().unwrap().chunk_index_test_stats(), (2, 0));
    }
    // A valid directory cannot make a corrupt compressed payload admissible.
    let mut bytes = leb2_v2_distinct_image(1.0);
    let index = FileIndex::decode(&bytes).unwrap();
    let (_, chunks) =
        decode_chunk_index(&bytes, index.find_body(0).unwrap().data_offset(), 0).unwrap();
    bytes[chunks[0].blob_offset() as usize] ^= 0xff;
    let file = LebFile::from_bytes("base_core.leb2".into(), bytes).unwrap();
    let expected = format!(
        "{:?}",
        file.index()
            .unwrap()
            .eval_leb2_body(file.bytes().unwrap(), 0, 0.25)
            .unwrap_err()
    );
    assert_eq!(
        format!("{:?}", file.eval_body_any(0, 0.25).unwrap_err()),
        expected
    );
    assert!(file.retained_chunk_index_bytes() > 0);
    assert_eq!(file.chunk_cache_stats().unwrap(), (0, 1, 0));
}

#[test]
fn nonchunked_formats_retain_no_chunk_directory() {
    for bytes in [leb1_const_image(0, 1.0, 0.0, 2.0), leb2_v1_distinct_image()] {
        let file = LebFile::from_bytes("opaque.dat".into(), bytes).unwrap();
        file.eval_body_any(0, 0.25).unwrap();
        assert_eq!(file.retained_chunk_index_bytes(), 0);
        assert_eq!(file.chunk_cache().unwrap().chunk_index_test_stats(), (0, 0));
    }
}

#[test]
fn original_names_use_canonical_tier_and_sdk_file_views() {
    for tier in [Tier::Base, Tier::Medium, Tier::Extended] {
        for name in [
            format!("ephemeris_{}.leb", tier.as_name()),
            format!("ephemeris_{}_planets.leb", tier.as_name()),
            format!("{}_core.leb2", tier.as_name()),
        ] {
            let mut ctx = EphemerisContext::new()
                .with_mode(ProviderMode::SealedLeb)
                .with_tier(tier);
            ctx.attach_leb_file(
                LebFile::from_bytes(name, leb1_const_image(0, 1.0, 0.0, 2.0)).unwrap(),
            );
            let req = ServingRequest {
                ctx: &ctx,
                delta_t_model: DeltaTModel::Smh2016,
                iers_enabled: false,
                observer: None,
                tier_policy: TierPolicy::TierPriority,
                skip_closed: false,
            };
            with_shared_calc_engine(req, |_| Ok(())).unwrap();
        }
    }
}

#[test]
fn selected_segment_borrows_the_decoded_row_and_preserves_bounds() {
    let bytes = leb2_v2_distinct_image(1.0);
    let index = FileIndex::decode(&bytes).unwrap();
    let entry = index.find_body(0).unwrap();
    let per_segment = entry.components() as usize * (entry.degree() as usize + 1);
    let words: Vec<_> = (0..per_segment * entry.segment_count() as usize)
        .map(|word| word as f64)
        .collect();
    let row = super::super::index::slice_one_segment(&words, entry, 1, 0).unwrap();
    assert_eq!(row.len(), per_segment);
    assert_eq!(row.as_ptr(), words[per_segment..].as_ptr());
    assert!(
        super::super::index::slice_one_segment(&words, entry, entry.segment_count() as usize, 0)
            .is_err()
    );
    assert!(
        super::super::index::slice_one_segment(&words[..words.len() - 1], entry, 0, 0).is_err()
    );
}
