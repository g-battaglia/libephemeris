//! Small record cache shared by the lazily admitted segments of one kernel.

use super::*;
use std::collections::VecDeque;
use std::sync::Mutex;

const SHARDS: usize = 8;
const ENTRIES_PER_SHARD: usize = 16;
const BYTES_PER_SHARD: usize = 128 * 1024;
type Key = (usize, usize);
type Record = Arc<[f64]>;

#[derive(Default)]
struct Shard {
    records: VecDeque<(Key, Record)>,
    bytes: usize,
}

impl Shard {
    fn get(&mut self, key: Key) -> Option<Record> {
        let index = self.records.iter().position(|(found, _)| *found == key)?;
        let (_, record) = self.records.remove(index).expect("record in shard");
        self.records.push_back((key, Arc::clone(&record)));
        Some(record)
    }
}

#[derive(Default)]
pub(super) struct RecordCache {
    shards: [Mutex<Shard>; SHARDS],
}

impl RecordCache {
    pub(super) fn acquire(
        &self,
        offset: usize,
        length: usize,
        read: impl FnOnce() -> Result<Record>,
    ) -> Result<Record> {
        let key = (offset, length);
        let lock = &self.shards[((offset / 8).wrapping_mul(0x9e37_79b9) >> 12) % SHARDS];
        if let Some(record) = lock.lock().unwrap_or_else(|e| e.into_inner()).get(key) {
            return Ok(record);
        }
        let record = read()?;
        let size = std::mem::size_of_val(record.as_ref());
        if size > BYTES_PER_SHARD {
            return Ok(record);
        }
        let mut shard = lock.lock().unwrap_or_else(|e| e.into_inner());
        if let Some(existing) = shard.get(key) {
            return Ok(existing);
        }
        while shard.bytes + size > BYTES_PER_SHARD || shard.records.len() >= ENTRIES_PER_SHARD {
            let (_, old) = shard.records.pop_front().expect("resident record");
            shard.bytes -= std::mem::size_of_val(old.as_ref());
        }
        shard.records.push_back((key, Arc::clone(&record)));
        shard.bytes += size;
        Ok(record)
    }
}

#[derive(Clone)]
pub(super) struct ImageRecords {
    bytes: Arc<Vec<u8>>,
    endian: DafEndian,
    cache: Arc<RecordCache>,
}

impl ImageRecords {
    pub(super) fn new(kernel: &SpkKernel, cache: Arc<RecordCache>) -> Self {
        Self {
            bytes: Arc::clone(&kernel.bytes),
            endian: kernel.endian,
            cache,
        }
    }

    pub(super) fn contains(&self, offset: usize, words: usize) -> bool {
        words
            .checked_mul(8)
            .and_then(|size| offset.checked_add(size))
            .is_some_and(|end| end <= self.bytes.len())
    }

    pub(super) fn record(&self, offset: usize, words: usize) -> Result<Record> {
        if !self.contains(offset, words) {
            return Err(Error::SpkCorrupt {
                what: "SPK Chebyshev record".into(),
                offset,
            });
        }
        self.cache.acquire(offset, words, || {
            Ok(self.bytes[offset..offset + words * 8]
                .chunks_exact(8)
                .map(|word| self.endian.read_f64(word))
                .collect())
        })
    }
}

impl std::fmt::Debug for ImageRecords {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ImageRecords")
            .field("image_bytes", &self.bytes.len())
            .finish()
    }
}

impl PartialEq for ImageRecords {
    fn eq(&self, other: &Self) -> bool {
        self.endian == other.endian
            && (Arc::ptr_eq(&self.bytes, &other.bytes) || self.bytes == other.bytes)
    }
}
