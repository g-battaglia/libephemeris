//! Process-wide retention for runtime-owned LEB images. Each shard has its
//! own least-recently-used queue and byte quota. No lock spans a decode.

use super::{TrackedVec, index::CacheKey};
use std::collections::{BTreeMap, HashMap};
use std::hash::{DefaultHasher, Hash, Hasher};
use std::num::NonZeroU64;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering::Relaxed};
use std::sync::{Arc, Mutex, OnceLock};

const SHARDS: usize = 8;
type Key = (NonZeroU64, CacheKey);
type Words = Arc<TrackedVec<f64>>;

#[derive(Debug)]
struct Entry {
    words: Words,
    access: u64,
}

#[derive(Debug, Default)]
struct Shard {
    entries: HashMap<Key, Entry>,
    age: BTreeMap<u64, Key>,
    clock: u64,
    bytes: usize,
    limit: usize,
}

fn size(words: &Words) -> usize {
    words.capacity() * size_of::<f64>()
}

impl Shard {
    fn tick(&mut self) -> u64 {
        if self.clock == u64::MAX {
            let order = std::mem::take(&mut self.age);
            self.clock = 0;
            for (_, key) in order {
                self.clock += 1;
                self.entries.get_mut(&key).expect("LRU entry").access = self.clock;
                self.age.insert(self.clock, key);
            }
        }
        self.clock += 1;
        self.clock
    }

    fn hit(&mut self, key: Key) -> Option<Words> {
        let access = self.tick();
        let entry = self.entries.get_mut(&key)?;
        self.age.remove(&entry.access);
        entry.access = access;
        self.age.insert(access, key);
        Some(Arc::clone(&entry.words))
    }

    fn remove(&mut self, key: &Key) {
        if let Some(entry) = self.entries.remove(key) {
            self.age.remove(&entry.access);
            self.bytes -= size(&entry.words);
        }
    }

    fn make_room(&mut self, extra: usize) {
        while self.bytes > self.limit - extra {
            let (&_, &key) = self.age.first_key_value().expect("resident LRU entry");
            self.remove(&key);
        }
    }
}

#[derive(Debug)]
pub(crate) struct ProcessCache {
    shards: [Mutex<Shard>; SHARDS],
    // Resizing and its diagnostics are rare. Cache hits only touch one shard.
    budget: Mutex<usize>,
}

impl ProcessCache {
    fn new(bytes: usize) -> Arc<Self> {
        Arc::new(Self {
            shards: std::array::from_fn(|i| {
                Mutex::new(Shard {
                    limit: Self::quota(bytes, i),
                    ..Shard::default()
                })
            }),
            budget: Mutex::new(bytes),
        })
    }

    fn quota(bytes: usize, shard: usize) -> usize {
        bytes / SHARDS + usize::from(shard < bytes % SHARDS)
    }

    fn shard(&self, key: Key) -> &Mutex<Shard> {
        let mut hash = DefaultHasher::new();
        key.hash(&mut hash);
        &self.shards[hash.finish() as usize % SHARDS]
    }

    fn resize(&self, bytes: usize) {
        let mut budget = self.budget.lock().unwrap_or_else(|e| e.into_inner());
        if *budget == bytes {
            return;
        }
        for (i, shard) in self.shards.iter().enumerate() {
            let mut shard = shard.lock().unwrap_or_else(|e| e.into_inner());
            shard.limit = Self::quota(bytes, i);
            shard.make_room(0);
        }
        *budget = bytes;
    }

    fn clear_image(&self, image: NonZeroU64) {
        for shard in &self.shards {
            let mut shard = shard.lock().unwrap_or_else(|e| e.into_inner());
            let keys: Vec<_> = shard
                .entries
                .keys()
                .filter(|(owner, _)| *owner == image)
                .copied()
                .collect();
            for key in keys {
                shard.remove(&key);
            }
        }
    }

    fn usage(&self, image: Option<NonZeroU64>) -> (usize, usize, usize) {
        let budget = self.budget.lock().unwrap_or_else(|e| e.into_inner());
        let (mut bytes, mut entries) = (0, 0);
        for shard in &self.shards {
            let shard = shard.lock().unwrap_or_else(|e| e.into_inner());
            if let Some(image) = image {
                for ((owner, _), entry) in &shard.entries {
                    if *owner == image {
                        bytes += size(&entry.words);
                        entries += 1;
                    }
                }
            } else {
                bytes += shard.bytes;
                entries += shard.entries.len();
            }
        }
        (bytes, *budget, entries)
    }
}

static PROCESS_CACHE: OnceLock<Arc<ProcessCache>> = OnceLock::new();

pub(crate) fn configured_cache(bytes: Option<usize>) -> Option<Arc<ProcessCache>> {
    match bytes {
        Some(bytes) => Some(Arc::clone(
            PROCESS_CACHE.get_or_init(|| ProcessCache::new(bytes)),
        )),
        None => PROCESS_CACHE.get().cloned(),
    }
}

pub(crate) fn set_process_cache_bytes(bytes: usize) {
    PROCESS_CACHE
        .get_or_init(|| ProcessCache::new(bytes))
        .resize(bytes);
}

pub(crate) fn process_cache_usage() -> Option<(usize, usize, usize)> {
    PROCESS_CACHE.get().map(|cache| cache.usage(None))
}

/// A lease tied to the immutable image, with per-image hit/miss counters.
#[derive(Debug)]
pub(super) struct ManagedCache {
    pool: Arc<ProcessCache>,
    image: NonZeroU64,
    enabled: AtomicBool,
    hits: AtomicU64,
    misses: AtomicU64,
}

impl ManagedCache {
    pub(super) fn new(pool: Arc<ProcessCache>, image: NonZeroU64) -> Self {
        Self {
            pool,
            image,
            enabled: AtomicBool::new(true),
            hits: AtomicU64::new(0),
            misses: AtomicU64::new(0),
        }
    }

    pub(super) fn enabled(&self) -> bool {
        self.enabled.load(std::sync::atomic::Ordering::Acquire)
    }

    pub(super) fn disable(&self) {
        self.enabled
            .store(false, std::sync::atomic::Ordering::Release);
        self.clear();
    }

    pub(super) fn clear(&self) {
        self.pool.clear_image(self.image);
    }

    pub(super) fn memory_usage(&self) -> (usize, usize) {
        let (bytes, limit, _) = self.pool.usage(Some(self.image));
        (bytes, limit)
    }

    pub(super) fn stats(&self) -> (u64, u64, usize) {
        let (_, _, entries) = self.pool.usage(Some(self.image));
        (self.hits.load(Relaxed), self.misses.load(Relaxed), entries)
    }

    pub(super) fn get_or_decode(
        &self,
        key: CacheKey,
        decode: impl FnOnce() -> super::Result<TrackedVec<f64>>,
    ) -> super::Result<Words> {
        let key = (self.image, key);
        let shard = self.pool.shard(key);
        if let Some(words) = shard.lock().unwrap_or_else(|e| e.into_inner()).hit(key) {
            self.hits.fetch_add(1, Relaxed);
            return Ok(words);
        }
        self.misses.fetch_add(1, Relaxed);
        let words = Arc::new(decode()?);
        let bytes = size(&words);
        let mut shard = shard.lock().unwrap_or_else(|e| e.into_inner());
        if !self.enabled() || bytes == 0 || bytes > shard.limit {
            return Ok(words);
        }
        if let Some(existing) = shard.hit(key) {
            return Ok(existing);
        }
        shard.make_room(bytes);
        let access = shard.tick();
        shard.age.insert(access, key);
        shard.entries.insert(
            key,
            Entry {
                words: Arc::clone(&words),
                access,
            },
        );
        shard.bytes += bytes;
        Ok(words)
    }
}

impl Drop for ManagedCache {
    fn drop(&mut self) {
        self.clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn images_share_retention_and_resize_keeps_live_borrows_valid() {
        let pool = ProcessCache::new(8 * SHARDS);
        let id = |i| NonZeroU64::new(i).unwrap();
        let first = ManagedCache::new(Arc::clone(&pool), id(1));
        let second = ManagedCache::new(Arc::clone(&pool), id(2));
        let k = CacheKey::Body(4, 1);
        let shard = pool.shard((id(1), k));
        let j = (0..)
            .map(|i| CacheKey::Body(5, i))
            .find(|key| std::ptr::eq(pool.shard((id(2), *key)), shard))
            .unwrap();
        let decode = |value| {
            Ok(TrackedVec {
                values: vec![value],
                _allocation: None,
            })
        };
        let old = first.get_or_decode(k, || decode(-0.0)).unwrap();
        assert_eq!(first.stats(), (0, 1, 1));
        let good = first.get_or_decode(k, || panic!("expected a hit")).unwrap();
        assert!(Arc::ptr_eq(&good, &old));
        second.get_or_decode(j, || decode(2.0)).unwrap();
        assert_eq!(first.stats(), (1, 1, 0));
        assert_eq!(pool.usage(None), (8, 8 * SHARDS, 1));
        pool.resize(0);
        assert_eq!(pool.usage(None), (0, 0, 0));
        assert_eq!(old[0].to_bits(), (-0.0_f64).to_bits());
        second.get_or_decode(j, || decode(3.0)).unwrap();
        assert_eq!(pool.usage(None).0, 0);
        pool.resize(16 * SHARDS);
        second.get_or_decode(j, || decode(4.0)).unwrap();
        drop(second);
        assert_eq!(pool.usage(None).0, 0);
    }
}
