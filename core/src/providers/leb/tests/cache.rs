//! Finite cache retention/concurrency risks independent of astronomy fixtures.
use super::*;

#[test]
fn capacity_and_entry_limits_account_allocations_not_lengths() {
    let cache = LebChunkCache::with_byte_budget(1024);
    for id in 0..17 {
        cache.get_or_decode((id, 0), || Ok(vec![1.0])).unwrap();
        let (bytes, limit) = cache.memory_usage();
        assert!(bytes <= limit);
        assert!(cache.stats().2 <= 16);
    }
    assert_eq!(cache.stats().2, 16);
    cache.clear();
    cache.set_byte_budget(64);
    let result = cache
        .get_or_decode((1, 0), || {
            let mut words = Vec::with_capacity(100);
            words.push(1.0);
            Ok(words)
        })
        .unwrap();
    assert_eq!(result[0], 1.0);
    assert_eq!(cache.memory_usage(), (0, 64));
    assert_eq!(cache.stats().2, 0); // length 8 bytes, allocation 800: do not retain
}

#[test]
fn auxiliary_walk_preserves_recent_body_coefficients_with_fixed_capacity() {
    let cache = LebChunkCache::new();
    let body = cache.get_or_decode((0, 0), || Ok(vec![2.0; 8])).unwrap();
    // Real indexed walks interleave body chunks with many small table ranges.
    // A table miss must not discard the frequently consumed body coefficients.
    for offset in 0..100 {
        let hit = cache
            .get_or_decode((0, 0), || panic!("hot body was evicted"))
            .unwrap();
        assert!(Arc::ptr_eq(&body, &hit));
        cache
            .get_or_decode_key(CacheKey::Table { offset, bytes: 8 }, || {
                Ok(TrackedVec {
                    values: vec![offset as f64],
                    _allocation: None,
                })
            })
            .unwrap();
        assert!(cache.stats().2 <= MAX_CACHED_CHUNKS);
    }
    assert_eq!(cache.memory_usage().0, 64 + 15 * 8);
}

#[test]
fn byte_pressure_evicts_only_needed_words_and_survives_access_clock_rollover() {
    let cache = LebChunkCache::with_byte_budget(64);
    for id in 0..4 {
        cache
            .get_or_decode((id, 0), || Ok(vec![id as f64; 2]))
            .unwrap();
    }
    cache.lock_inner().access_clock = u64::MAX;
    // Refresh one old entry across rollover, leaving 1/2/3 older than 0.
    cache.get_or_decode((0, 0), || unreachable!()).unwrap();
    let held = cache.get_or_decode((3, 0), || unreachable!()).unwrap();
    cache.get_or_decode((4, 0), || Ok(vec![4.0; 4])).unwrap();
    let inner = cache.lock_inner();
    assert_eq!(inner.resident_bytes, 64);
    for id in [0, 3, 4] {
        assert!(inner.entries.contains_key(&CacheKey::Body(id, 0)));
    }
    for id in [1, 2] {
        assert!(!inner.entries.contains_key(&CacheKey::Body(id, 0)));
    }
    assert_eq!(held.as_slice(), &[3.0, 3.0]);
}

#[test]
fn failed_decode_is_never_cached() {
    let cache = LebChunkCache::new();
    for _ in 0..2 {
        assert!(
            cache
                .get_or_decode((0, 0), || Err(Error::LebCorrupt {
                    what: "synthetic decode failure".into(),
                    offset: 0
                }))
                .is_err()
        );
    }
    assert_eq!(cache.stats(), (0, 2, 0));
    assert_eq!(cache.memory_usage().0, 0);
}

#[test]
fn concurrent_same_key_retains_one_allocation_and_shrink_wins_inflight_decode() {
    let cache = LebChunkCache::with_byte_budget(1024);
    let barrier = std::sync::Barrier::new(8);
    std::thread::scope(|scope| {
        for _ in 0..8 {
            let cache = &cache;
            let barrier = &barrier;
            scope.spawn(move || {
                let words = cache
                    .get_or_decode((0, 0), || {
                        barrier.wait();
                        Ok(vec![2.0; 8])
                    })
                    .unwrap();
                assert_eq!(words.as_slice(), &[2.0; 8]);
            });
        }
    });
    assert_eq!(cache.stats(), (0, 8, 1));
    assert_eq!(cache.memory_usage(), (64, 1024));
    cache.clear();
    let (tx, rx) = std::sync::mpsc::channel();
    let resume = std::sync::Barrier::new(2);
    std::thread::scope(|scope| {
        scope.spawn(|| {
            cache
                .get_or_decode((1, 0), || {
                    tx.send(()).unwrap();
                    resume.wait();
                    Ok(vec![3.0; 8])
                })
                .unwrap();
        });
        rx.recv().unwrap();
        cache.set_byte_budget(0);
        resume.wait();
    });
    assert_eq!(cache.memory_usage(), (0, 0));
    assert_eq!(cache.stats().2, 0);
}

#[test]
fn indexed_eviction_releases_retention_but_not_borrowed_words() {
    let resources = IndexedResources::new(IndexedLimits::default());
    let cache = LebChunkCache::with_indexed_resources(Arc::clone(&resources));
    let words = cache
        .get_or_decode_tracked((0, 0), || {
            let allocation = resources.working.reserve(64)?;
            Ok(TrackedVec {
                values: vec![2.0; 8],
                _allocation: Some(allocation),
            })
        })
        .unwrap();
    assert_eq!(resources.retained.used(), 64);
    assert_eq!(resources.working.used(), 64);
    cache.set_byte_budget(0);
    assert_eq!(resources.retained.used(), 0);
    assert_eq!(resources.working.used(), 64);
    assert_eq!(words.as_slice(), &[2.0; 8]);
    drop(words);
    assert_eq!(resources.working.used(), 0);
}

#[test]
fn indexed_store_serializes_cold_decode_and_rechecks_the_cache() {
    use std::sync::atomic::{AtomicUsize, Ordering};
    let resources = IndexedResources::new(IndexedLimits::default());
    let cache = LebChunkCache::with_indexed_resources(Arc::clone(&resources));
    let barrier = std::sync::Barrier::new(8);
    let calls = AtomicUsize::new(0);
    std::thread::scope(|scope| {
        for _ in 0..8 {
            scope.spawn(|| {
                barrier.wait();
                let words = cache
                    .get_or_decode_tracked((0, 0), || {
                        calls.fetch_add(1, Ordering::Relaxed);
                        let allocation = resources.working.reserve(64)?;
                        Ok(TrackedVec {
                            values: vec![2.0; 8],
                            _allocation: Some(allocation),
                        })
                    })
                    .unwrap();
                assert_eq!(words.as_slice(), &[2.0; 8]);
            });
        }
    });
    assert_eq!(calls.load(Ordering::Relaxed), 1);
    assert_eq!(cache.stats(), (7, 1, 1));
    cache.clear();
    assert_eq!(resources.working.used(), 0);
}

#[test]
fn shared_retention_limit_never_hides_live_uncached_allocations() {
    let resources = IndexedResources::new(IndexedLimits {
        cache_bytes: 64,
        ..IndexedLimits::default()
    });
    let left = LebChunkCache::with_indexed_resources(Arc::clone(&resources));
    let right = LebChunkCache::with_indexed_resources(Arc::clone(&resources));
    let decode = |value| {
        let allocation = resources.working.reserve(64)?;
        Ok(TrackedVec {
            values: vec![value; 8],
            _allocation: Some(allocation),
        })
    };
    let a = left.get_or_decode_tracked((0, 0), || decode(1.0)).unwrap();
    let b = right.get_or_decode_tracked((0, 0), || decode(2.0)).unwrap();
    assert_eq!(resources.retained.used(), 64);
    assert_eq!(resources.working.used(), 128);
    assert_eq!(right.stats().2, 0);
    assert_eq!(a[0], 1.0);
    assert_eq!(b[0], 2.0);
    drop(b);
    assert_eq!(resources.working.used(), 64);
    left.clear();
    drop(a);
    assert_eq!(resources.working.used(), 0);
}

#[test]
fn shared_pressure_evicts_only_enough_local_entries_and_keeps_live_leases() {
    let resources = IndexedResources::new(IndexedLimits {
        cache_bytes: 64,
        ..IndexedLimits::default()
    });
    let left = LebChunkCache::with_indexed_resources(Arc::clone(&resources));
    let right = LebChunkCache::with_indexed_resources(Arc::clone(&resources));
    let decode = |length| {
        let allocation = resources.working.reserve(length * 8)?;
        Ok(TrackedVec {
            values: vec![1.0; length],
            _allocation: Some(allocation),
        })
    };
    let held: Vec<_> = (0..3)
        .map(|id| left.get_or_decode_tracked((id, 0), || decode(2)).unwrap())
        .collect();
    right.get_or_decode_tracked((0, 0), || decode(2)).unwrap();
    left.get_or_decode_tracked((3, 0), || decode(4)).unwrap();
    assert_eq!(left.stats().2, 2);
    assert_eq!(right.stats().2, 1);
    assert_eq!(resources.retained.used(), 64);
    assert_eq!(resources.working.used(), 96);
    left.get_or_decode_tracked((2, 0), || panic!("recent retained entry was evicted"))
        .unwrap();
    left.clear();
    right.clear();
    assert_eq!(resources.retained.used(), 0);
    assert_eq!(resources.working.used(), 48);
    drop(held);
    assert_eq!(resources.working.used(), 0);
}
