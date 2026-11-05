//! Optional lossless compaction of the first byte-pressure LRU victim only.

use super::*;
use crate::leb::packing::PackedWords;

const INLINE_ROW_WORDS: usize = 64;

#[derive(Debug, Clone)]
pub(super) enum CachedStorage {
    Dense(Arc<TrackedVec<f64>>),
    Packed(Arc<PackedBody>),
}

#[derive(Debug)]
pub(super) struct PackedBody {
    words: PackedWords,
    charged_bytes: usize,
}

impl CachedStorage {
    pub(super) fn byte_capacity(&self) -> usize {
        match self {
            Self::Dense(words) => words.capacity() * std::mem::size_of::<f64>(),
            Self::Packed(words) => words.charged_bytes,
        }
    }

    fn with_segment<R>(
        &self,
        entry: &BodyEntry,
        local_idx: usize,
        buffer_segments: u32,
        body_id: i32,
        consume: impl FnOnce(&[f64]) -> Result<R>,
    ) -> Result<R> {
        match self {
            Self::Dense(words) => consume(slice_one_segment_in(
                words,
                entry,
                local_idx,
                buffer_segments,
                body_id,
            )?),
            Self::Packed(owner) => {
                let range = segment_word_range(
                    owner.words.row_count() * owner.words.row_width(),
                    entry,
                    local_idx,
                    buffer_segments,
                    body_id,
                )?;
                let mut row = [0.0; INLINE_ROW_WORDS];
                let selected = row
                    .get_mut(..range.len())
                    .ok_or_else(|| Error::LebCorrupt {
                        what: format!("LEB2 coefficient blob (body {body_id})"),
                        offset: 0,
                    })?;
                owner.words.decode_row(local_idx, selected)?;
                consume(selected)
            }
        }
    }
}

/// Conservative additional shared-owner and common entry-metadata charge.
/// The original dense map/Arc overhead remains under its documented exclusion;
/// each compact body pays for its owner and a conservative capacity charge for
/// the new storage tag across the common map, including possible next growth.
fn owner_fee(index_capacity: usize) -> usize {
    #[allow(dead_code)]
    struct DenseEntryLayout {
        words: Arc<TrackedVec<f64>>,
        last_access: u64,
        retention: Option<resources::Reservation>,
    }
    let extra =
        std::mem::size_of::<CachedWords>().saturating_sub(std::mem::size_of::<DenseEntryLayout>());
    std::mem::size_of::<PackedBody>()
        + 4 * std::mem::size_of::<usize>()
        + std::mem::align_of::<PackedBody>()
        + 2 * (index_capacity + MAX_CACHED_CHUNKS) * extra
}

impl LebChunkCache {
    /// Dense hits still take one ordinary cache lock and clone their raw Arc;
    /// they do not inspect a codec, reconstruct a row or allocate heap storage.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn with_cached_body_segment<R>(
        &self,
        key: (i32, u32),
        index: &FileIndex,
        entry: &BodyEntry,
        local_idx: usize,
        buffer_segments: u32,
        decode: impl FnOnce() -> Result<TrackedVec<f64>>,
        consume: impl FnOnce(&[f64]) -> Result<R>,
    ) -> Result<R> {
        self.get_or_decode_storage(CacheKey::Body(key.0, key.1), Some(index), decode)?
            .with_segment(entry, local_idx, buffer_segments, key.0, consume)
    }

    /// Called only for imminent byte eviction with a spare entry slot.
    /// Packing is pure internal work under the cold cache lock; user decoders
    /// and row callbacks still run outside it. The lock makes clear/shrink,
    /// concurrent hits and victim replacement atomic with this optional attempt.
    pub(super) fn compact_first_victim(
        &self,
        inner: &mut ChunkCacheInner,
        index: &FileIndex,
        incoming_bytes: usize,
    ) {
        // Public arbitrary-source caches and every Indexed cache retain their
        // original storage and live/retained reservation policy.
        if self.resources.is_some() || self.word_image_identity.is_none() {
            return;
        }
        let Some(key) = inner
            .entries
            .iter()
            .max_by_key(|(_, entry)| inner.access_clock.wrapping_sub(entry.last_access))
            .map(|(key, _)| *key)
        else {
            return;
        };
        let CacheKey::Body(body, _) = key else { return };
        let Some(victim) = inner.entries.get(&key) else {
            return;
        };
        let CachedStorage::Dense(raw) = &victim.storage else {
            return;
        };
        // A tracked allocation belongs to an Indexed lease. Never shrink,
        // transfer or replace a live lease through this contiguous-image mechanism.
        if raw._allocation.is_some() {
            return;
        }
        let Ok(entry) = index.find_body(body) else {
            return;
        };
        let Some(width) = (entry.degree() as usize)
            .checked_add(1)
            .and_then(|degree| degree.checked_mul(entry.components() as usize))
            .filter(|width| (1..=INLINE_ROW_WORDS).contains(width))
        else {
            return;
        };
        let old_bytes = raw.capacity() * std::mem::size_of::<f64>();
        let Ok(Some(words)) = PackedWords::pack(raw, width) else {
            return;
        };
        let Some(charged_bytes) = words
            .heap_bytes()
            .checked_add(owner_fee(inner.entries.capacity()))
        else {
            return;
        };
        let Some(total) = inner
            .resident_bytes
            .checked_sub(old_bytes)
            .and_then(|bytes| bytes.checked_add(charged_bytes))
            .and_then(|bytes| bytes.checked_add(incoming_bytes))
        else {
            return;
        };
        // Decline without changing the LRU victim or looking at another entry.
        // The ordinary admission loop below will perform its original eviction.
        if charged_bytes >= old_bytes || total > inner.max_bytes {
            return;
        }
        let packed = Arc::new(PackedBody {
            words,
            charged_bytes,
        });
        let Some(victim) = inner.entries.get_mut(&key) else {
            return;
        };
        victim.storage = CachedStorage::Packed(packed);
        // last_access deliberately does not change: this remains the same LRU
        // entry. Outstanding dense borrowers retain their original raw Arc.
        inner.resident_bytes = inner.resident_bytes - old_bytes + charged_bytes;
    }
}

#[cfg(test)]
#[path = "tests/pressure_victim.rs"]
mod tests;
