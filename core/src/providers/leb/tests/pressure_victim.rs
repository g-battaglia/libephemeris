use super::*;

fn entry(body: i32, width: u32) -> BodyEntry {
    BodyEntry {
        body_id: body,
        coord_type: COORD_ICRS_BARY,
        segment_count: 512,
        jd_start: 0.0,
        jd_end: 512.0,
        interval_days: 1.0,
        degree: 0,
        components: width,
        data_offset: 0,
        compressed_size: None,
        uncompressed_size: None,
    }
}

fn index(width: u32) -> FileIndex {
    let fixture = crate::test_data::core_bytes();
    let mut index = FileIndex::decode(fixture).unwrap();
    index.bodies = (0..20).map(|body| entry(body, width)).collect();
    index
}

fn raw(width: usize) -> TrackedVec<f64> {
    let patterns = [0x8000_0000_0000_0000, 0x7ff8_0123_4567_89ab, 1];
    TrackedVec {
        values: (0..512 * width)
            .map(|word| f64::from_bits(patterns[(word % width) % patterns.len()]))
            .collect(),
        _allocation: None,
    }
}

fn obtain(cache: &LebChunkCache, index: &FileIndex, body: i32, width: usize) -> CachedStorage {
    cache
        .get_or_decode_storage(CacheKey::Body(body, 0), Some(index), || Ok(raw(width)))
        .unwrap()
}

fn is_packed(cache: &LebChunkCache, body: i32) -> bool {
    cache
        .lock_inner()
        .entries
        .get(&CacheKey::Body(body, 0))
        .is_some_and(|entry| matches!(entry.storage, CachedStorage::Packed(_)))
}

#[test]
fn victim_only_admission_keeps_dense_hits_and_is_not_sticky() {
    let index = index(3);
    let bytes = 512 * 3 * 8;
    let cache = LebChunkCache::with_image_byte_budget(2 * bytes + bytes / 2);
    obtain(&cache, &index, 0, 3);
    obtain(&cache, &index, 1, 3);
    assert!(!is_packed(&cache, 0) && !is_packed(&cache, 1));
    obtain(&cache, &index, 2, 3);
    assert!(is_packed(&cache, 0));
    assert!(!is_packed(&cache, 1) && !is_packed(&cache, 2));
    let inner = cache.lock_inner();
    assert_eq!(inner.entries.len(), 3);
    assert_eq!(
        inner.resident_bytes,
        inner
            .entries
            .values()
            .map(|entry| entry.storage.byte_capacity())
            .sum::<usize>()
    );
    assert!(inner.resident_bytes <= inner.max_bytes);
    drop(inner);
    // The first LRU is already packed. Do not search/pack the second entry:
    // ordinary eviction removes both oldest entries before retaining this raw one.
    obtain(&cache, &index, 3, 3);
    let inner = cache.lock_inner();
    assert!(!inner.entries.contains_key(&CacheKey::Body(0, 0)));
    assert!(!inner.entries.contains_key(&CacheKey::Body(1, 0)));
    assert!(
        inner
            .entries
            .values()
            .all(|entry| matches!(entry.storage, CachedStorage::Dense(_)))
    );
}

#[test]
fn packed_rows_keep_64_word_stack_bits_errors_and_owner_release() {
    let index = index(64);
    let bytes = 512 * 64 * 8;
    let cache = LebChunkCache::with_image_byte_budget(2 * bytes + bytes / 2);
    obtain(&cache, &index, 0, 64);
    obtain(&cache, &index, 1, 64);
    obtain(&cache, &index, 2, 64);
    let stored = cache
        .get_or_decode_storage(CacheKey::Body(0, 0), Some(&index), || {
            panic!("packed hit decoded")
        })
        .unwrap();
    let CachedStorage::Packed(owner) = &stored else {
        panic!("first victim was not compacted")
    };
    let weak = Arc::downgrade(owner);
    let expected = raw(64);
    stored
        .with_segment(&index.bodies[0], 511, 512, 0, |row| {
            assert_eq!(
                row.iter().map(|value| value.to_bits()).collect::<Vec<_>>(),
                expected.values[511 * 64..]
                    .iter()
                    .map(|value| value.to_bits())
                    .collect::<Vec<_>>()
            );
            Ok(())
        })
        .unwrap();
    let dense_error = slice_one_segment_in(&expected, &index.bodies[0], 512, 512, 0).unwrap_err();
    let packed_error = stored
        .with_segment(&index.bodies[0], 512, 512, 0, |_| Ok(()))
        .unwrap_err();
    assert_eq!(dense_error.to_string(), packed_error.to_string());
    cache.set_byte_budget(0);
    assert_eq!(cache.memory_usage(), (0, 0));
    assert!(weak.upgrade().is_some());
    drop(stored);
    assert!(weak.upgrade().is_none());
    println!(
        "PRESSURE_VICTIM_LAYOUT entry={} storage={} owner={} owner_fee={} inline_row_bytes={}",
        std::mem::size_of::<CachedWords>(),
        std::mem::size_of::<CachedStorage>(),
        std::mem::size_of::<PackedBody>(),
        owner_fee(32),
        std::mem::size_of::<[f64; 64]>()
    );
}

#[test]
fn count_table_public_and_large_row_pressure_keep_original_dense_eviction() {
    let index = index(3);
    let bytes = 512 * 3 * 8;
    let public = LebChunkCache::with_byte_budget(2 * bytes + bytes / 2);
    for body in 0..3 {
        obtain(&public, &index, body, 3);
    }
    assert_eq!(public.stats().2, 2);
    assert!(!is_packed(&public, 1));
    let counted = LebChunkCache::with_image_byte_budget(1024 * 1024);
    for body in 0..17 {
        obtain(&counted, &index, body, 3);
    }
    assert_eq!(counted.stats().2, MAX_CACHED_CHUNKS);
    assert!(
        counted
            .lock_inner()
            .entries
            .values()
            .all(|entry| matches!(entry.storage, CachedStorage::Dense(_)))
    );
    let tables = LebChunkCache::with_image_byte_budget(2 * bytes + bytes / 2);
    tables
        .get_or_decode_key(CacheKey::Table { offset: 0, bytes }, || Ok(raw(3)))
        .unwrap();
    obtain(&tables, &index, 0, 3);
    obtain(&tables, &index, 1, 3);
    assert!(
        tables
            .lock_inner()
            .entries
            .values()
            .all(|entry| matches!(entry.storage, CachedStorage::Dense(_)))
    );
    let wide_index = self::index(65);
    let wide_bytes = 512 * 65 * 8;
    let wide = LebChunkCache::with_image_byte_budget(2 * wide_bytes + wide_bytes / 2);
    for body in 0..3 {
        obtain(&wide, &wide_index, body, 65);
    }
    assert_eq!(wide.stats().2, 2);
    assert!(
        wide.lock_inner()
            .entries
            .values()
            .all(|entry| matches!(entry.storage, CachedStorage::Dense(_)))
    );
}

#[test]
fn failed_decode_and_indexed_leases_never_enter_optional_compaction() {
    let index = index(3);
    let bytes = 512 * 3 * 8;
    let cache = LebChunkCache::with_image_byte_budget(2 * bytes + bytes / 2);
    obtain(&cache, &index, 0, 3);
    obtain(&cache, &index, 1, 3);
    let before = cache.memory_usage();
    let result = cache.get_or_decode_storage(CacheKey::Body(2, 0), Some(&index), || {
        Err(Error::LebCorrupt {
            what: "test refusal".into(),
            offset: 7,
        })
    });
    assert!(matches!(result, Err(Error::LebCorrupt { offset: 7, .. })));
    assert_eq!(cache.memory_usage(), before);
    assert!(!is_packed(&cache, 0));
    let limits = IndexedLimits {
        cache_bytes: 2 * bytes + bytes / 2,
        working_bytes: 8 * bytes,
        ..IndexedLimits::default()
    };
    let resources = IndexedResources::new(limits);
    let indexed = LebChunkCache::with_indexed_resources(Arc::clone(&resources));
    let held = indexed
        .get_or_decode_storage(CacheKey::Body(0, 0), Some(&index), || {
            let mut data = raw(3);
            data._allocation = Some(resources.working.reserve(bytes)?);
            Ok(data)
        })
        .unwrap();
    for body in 1..3 {
        indexed
            .get_or_decode_storage(CacheKey::Body(body, 0), Some(&index), || {
                let mut data = raw(3);
                data._allocation = Some(resources.working.reserve(bytes)?);
                Ok(data)
            })
            .unwrap();
    }
    assert_eq!(indexed.stats().2, 2);
    assert!(
        indexed
            .lock_inner()
            .entries
            .values()
            .all(|entry| matches!(entry.storage, CachedStorage::Dense(_)))
    );
    indexed.clear();
    assert_eq!(resources.retained.used(), 0);
    assert_eq!(resources.working.used(), bytes);
    drop(held);
    assert_eq!(resources.working.used(), 0);
}

#[test]
fn concurrent_dense_borrower_survives_victim_conversion_and_clear() {
    let index = index(3);
    let bytes = 512 * 3 * 8;
    let cache = LebChunkCache::with_image_byte_budget(2 * bytes + bytes / 2);
    let CachedStorage::Dense(held) = obtain(&cache, &index, 0, 3) else {
        unreachable!()
    };
    let weak = Arc::downgrade(&held);
    obtain(&cache, &index, 1, 3);
    let (ready_tx, ready_rx) = std::sync::mpsc::channel();
    let (release_tx, release_rx) = std::sync::mpsc::channel();
    let borrower = std::thread::spawn(move || {
        ready_tx.send(()).unwrap();
        let _ = release_rx.recv();
        assert_eq!(held.values[0].to_bits(), 0x8000_0000_0000_0000);
        assert_eq!(held.values[1].to_bits(), 0x7ff8_0123_4567_89ab);
    });
    ready_rx.recv().unwrap();
    obtain(&cache, &index, 2, 3);
    assert!(is_packed(&cache, 0));
    cache.clear();
    assert!(weak.upgrade().is_some());
    release_tx.send(()).unwrap();
    borrower.join().unwrap();
    assert!(weak.upgrade().is_none());
    assert_eq!(cache.memory_usage().0, 0);
}
