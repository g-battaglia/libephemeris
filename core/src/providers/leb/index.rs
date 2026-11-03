//! Index implementation for the parent module.
//! Public entry points and types are exposed through the parent API.

use super::*;
use std::sync::{Arc, Mutex, MutexGuard};

#[path = "pressure_victim.rs"]
mod pressure_victim;
use pressure_victim::CachedStorage;

/// Bound on memoized decoded chunks per [`LebChunkCache`].
///
/// Body chunks, row windows and auxiliary tables share this bound. When admission needs room,
/// evict the least recently accessed entry, preserving other useful chunks.
/// The bound caps retained coefficient entries, never correctness: every entry is
/// a pure memo of the owner's borrowed image.
const MAX_CACHED_CHUNKS: usize = 16;

/// Maximum coefficient allocation retained for one contiguous-image window.
const MAX_COEFFICIENT_WINDOW_BYTES: usize = 64 * 1024;

/// A complete row must fit; wider/empty layouts keep their original path.
fn coefficient_window(
    entry: &BodyEntry,
    local: usize,
    segments: u32,
) -> Option<std::ops::Range<u32>> {
    let row_bytes = (entry.degree() as usize)
        .checked_add(1)?
        .checked_mul(entry.components() as usize)?
        .checked_mul(std::mem::size_of::<f64>())?;
    let rows = MAX_COEFFICIENT_WINDOW_BYTES.checked_div(row_bytes)?;
    if rows == 0 || local >= segments as usize {
        return None;
    }
    let first = local / rows * rows;
    let end = first + rows.min(segments as usize - first);
    Some(u32::try_from(first).ok()?..u32::try_from(end).ok()?)
}

/// Default retained decoded-coefficient budget per image: 1 MiB.
///
/// Image/context clones share this budget; it is not multiplied by request
/// count. Contiguous-image windows have a separate 64 KiB entry bound;
/// other oversized entries still evaluate without retention. Whole-file bytes,
/// zstd workspace and concurrently borrowed/decoding chunks are additional
/// memory, so this is a cache bound, not a total-RSS promise.
pub const DEFAULT_CHUNK_CACHE_BYTES: usize = 1024 * 1024;

/// Decoded LEB2 chunks and directories retained for one immutable image.
///
/// Entries are keyed by body and chunk index, with an additional segment
/// offset for bounded windows; `u32::MAX` identifies a monolithic LEB2 v1 blob.
/// Raw-index providers own private caches;
/// file-backed providers share the corresponding [`LebFile`] cache. A path
/// replacement creates another owner and never reuses an old image's words.
///
/// A caller-created cache must be dedicated to one immutable container and
/// its matching [`FileIndex`] for its entire lifetime. Raw cached evaluators
/// do not hash or validate image identity before a hit. Reusing a cache with
/// different bytes can return another image's words or directory offsets.
/// [`Self::clear`] and a zero byte budget do **not** make that reuse valid:
/// both retain directories. Construct a new cache for a different image.
///
/// Directory metadata is initialized lazily once per body and is additional
/// to the retained coefficient budget. Decompression runs outside the
/// coefficient-cache lock; no callbacks execute under either lock. Entry count
/// and the configured byte budget bound retained coefficients, not total RAM;
/// oversized chunks can be evaluated without retention.
pub struct LebChunkCache {
    inner: Mutex<ChunkCacheInner>,
    chunk_indexes: Mutex<ChunkIndexCacheInner>,
    resources: Option<Arc<IndexedResources>>,
    // Additional per-cache metadata: one optional nonzero u64 (8 bytes).
    // Only private immutable-image constructors enable numerical reuse.
    word_image_identity: Option<std::num::NonZeroU64>,
    // Only owned/static/shared contiguous images use partial decoded windows.
    // Indexed leases and mmap acquisition retain their original decode route.
    segment_windows: bool,
}

/// Immutable per-body directories paired with the cache owner's snapshot.
/// Metadata is bounded by the requested bodies and their stored directories,
/// not by request count; it survives coefficient eviction and cooling.
#[derive(Default)]
struct ChunkIndexCacheInner {
    entries: HashMap<i32, Arc<TrackedVec<ChunkEntry>>>,
    /// Directory vector capacities, excluding map/Arc/allocator overhead.
    resident_bytes: usize,
    #[cfg(test)]
    parse_calls: usize,
}

impl std::fmt::Debug for LebChunkCache {
    /// Report aggregate residency only; never print coefficient or directory data.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let (resident, limit) = self.memory_usage();
        f.debug_struct("LebChunkCache")
            .field("coefficient_bytes", &resident)
            .field("coefficient_limit", &limit)
            .field("coefficient_stats", &self.stats())
            .field("chunk_index_bytes", &self.retained_chunk_index_bytes())
            .finish_non_exhaustive()
    }
}

/// Mutable state behind [`LebChunkCache`]: decoded words plus
/// hit/miss counters for measurement.
#[derive(Debug, Default)]
struct ChunkCacheInner {
    /// Decoded body/chunk/window words and auxiliary-table words.
    entries: HashMap<CacheKey, CachedWords>,
    /// Sum of charged retained storage; historical dense map overhead is excluded.
    resident_bytes: usize,
    /// Maximum retained coefficient allocation; zero disables retention.
    max_bytes: usize,
    /// Cache hits (the decode pipeline was skipped).
    hits: u64,
    /// Decode attempts on cache miss; retention is optional.
    misses: u64,
    /// Wrapping access sequence; used only to choose a replacement victim.
    access_clock: u64,
}

/// Body/chunk keys and raw-table byte ranges are distinct namespaces.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(super) enum CacheKey {
    Body(i32, u32),
    Window {
        body: i32,
        chunk: u32,
        first_segment: u32,
    },
    Table {
        offset: usize,
        bytes: usize,
    },
}

#[derive(Debug)]
struct CachedWords {
    storage: CachedStorage,
    last_access: u64,
    // Retention ends on eviction, while the words' live lease follows readers.
    _retention: Option<resources::Reservation>,
}

impl ChunkCacheInner {
    fn next_access(&mut self) -> u64 {
        self.access_clock = self.access_clock.wrapping_add(1);
        self.access_clock
    }

    fn get(&mut self, key: &CacheKey) -> Option<CachedStorage> {
        let access = self.next_access();
        let entry = self.entries.get_mut(key)?;
        entry.last_access = access;
        Some(entry.storage.clone())
    }

    /// Drop retention credit only; active readers keep their live allocation.
    /// Modular ages handle clock rollover without changing stored coefficients.
    fn evict_oldest(&mut self) -> bool {
        let Some(key) = self
            .entries
            .iter()
            .max_by_key(|(_, entry)| self.access_clock.wrapping_sub(entry.last_access))
            .map(|(key, _)| *key)
        else {
            return false;
        };
        if let Some(entry) = self.entries.remove(&key) {
            self.resident_bytes -= entry.storage.byte_capacity();
        }
        true
    }
}

impl LebChunkCache {
    pub(super) fn with_indexed_resources(resources: Arc<IndexedResources>) -> Self {
        let mut cache = Self::with_image_byte_budget(resources.limits.cache_bytes);
        cache.resources = Some(resources);
        cache
    }

    /// Empty cache with [`DEFAULT_CHUNK_CACHE_BYTES`] retained-coefficient budget.
    /// Dedicate it to one immutable image as described by [`LebChunkCache`].
    pub fn new() -> Self {
        Self::with_byte_budget(DEFAULT_CHUNK_CACHE_BYTES)
    }

    /// Construct with a retained decoded-byte limit. Zero disables retention;
    /// oversize chunks are decoded normally on every request, never refused.
    /// The entry count remains bounded at 16. Decode scratch space and live
    /// callers' references are not included in retained-byte accounting.
    pub fn with_byte_budget(max_bytes: usize) -> Self {
        Self {
            inner: Mutex::new(ChunkCacheInner {
                max_bytes,
                ..ChunkCacheInner::default()
            }),
            chunk_indexes: Mutex::new(ChunkIndexCacheInner::default()),
            resources: None,
            word_image_identity: None,
            segment_windows: false,
        }
    }

    /// Pair this private cache with one immutable native index/source owner.
    /// Public arbitrary-source caches remain ineligible, including after moves.
    pub(super) fn with_image_byte_budget(max_bytes: usize) -> Self {
        let mut cache = Self::with_byte_budget(max_bytes);
        cache.word_image_identity = crate::pure_cache::new_word_image_identity();
        cache
    }

    pub(super) fn with_contiguous_window_byte_budget(max_bytes: usize) -> Self {
        let mut cache = Self::with_image_byte_budget(max_bytes);
        cache.segment_windows = true;
        cache
    }

    pub(super) fn word_image_identity(&self) -> Option<std::num::NonZeroU64> {
        self.word_image_identity
    }

    /// Retained directory vector capacities, excluding map/Arc/allocator overhead.
    /// No parsing, file reads or hashing occur; coefficient budgets are independent.
    pub(crate) fn retained_chunk_index_bytes(&self) -> usize {
        self.lock_chunk_indexes().resident_bytes
    }

    /// Lock directory metadata without holding the coefficient guard.
    /// As with the coefficient guard, recover poisoning rather than invalidate
    /// an otherwise usable immutable snapshot.
    fn lock_chunk_indexes(&self) -> MutexGuard<'_, ChunkIndexCacheInner> {
        self.chunk_indexes
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
    }

    /// Reuse a body's immutable directory, preserving lazy decoder errors.
    /// Only directory parsing occurs under this guard, never zstd or astronomy.
    /// Failures publish nothing; shared concurrent readers parse successfully once.
    fn chunk_index(
        &self,
        bytes: LebSource<'_>,
        entry: &BodyEntry,
        body_id: i32,
    ) -> Result<Arc<TrackedVec<ChunkEntry>>> {
        let mut inner = self.lock_chunk_indexes();
        if let Some(chunks) = inner.entries.get(&body_id) {
            return Ok(Arc::clone(chunks));
        }
        if let LebSource::Positional(file) = bytes {
            let header_bytes = bytes.blob(
                entry.data_offset,
                CHUNK_INDEX_HEADER_SIZE as u64,
                "LEB chunk header",
            )?;
            let header = ChunkIndexHeader::decode(&header_bytes, 0)?;
            let additional = (header.chunk_count() as usize)
                .checked_mul(std::mem::size_of::<ChunkEntry>())
                .and_then(|count| count.checked_add(inner.resident_bytes))
                .and_then(|count| count.checked_add(file.main_index_bytes))
                .ok_or_else(|| Error::LebCorrupt {
                    what: "LEB directory allocation".into(),
                    offset: 0,
                })?;
            // The directory lock serializes admission and publication; no
            // coefficient decode or astronomy callback runs under this guard.
            access::limit("LEB index metadata", additional, file.limits.metadata_bytes)?;
        }
        #[cfg(test)]
        {
            inner.parse_calls += 1;
        }
        let (_, chunks) = decode::decode_chunk_index_source(bytes, entry.data_offset, body_id)?;
        if chunks.is_empty() {
            return Err(Error::LebCorrupt {
                what: format!("LEB2 chunk index (body {body_id})"),
                offset: 0,
            });
        }
        let size = chunks
            .capacity()
            .saturating_mul(std::mem::size_of::<ChunkEntry>());
        let chunks = Arc::new(chunks);
        inner.resident_bytes = inner.resident_bytes.saturating_add(size);
        inner.entries.insert(body_id, Arc::clone(&chunks));
        Ok(chunks)
    }

    /// Test-only initialization evidence; never changes public coefficient counters.
    #[cfg(test)]
    pub(super) fn chunk_index_test_stats(&self) -> (usize, usize) {
        let inner = self.lock_chunk_indexes();
        (inner.parse_calls, inner.entries.len())
    }

    /// `(resident_coefficient_bytes, byte_limit)` for memory diagnostics.
    /// Dense entries count vector capacity. Compact entries additionally charge
    /// their shared owner and conservative extra map metadata; historical dense
    /// map/Arc overhead remains excluded.
    pub fn memory_usage(&self) -> (usize, usize) {
        let inner = self.lock_inner();
        (inner.resident_bytes, inner.max_bytes)
    }

    /// Change the retained-byte budget. Shrinking below current residency
    /// releases retained words; zero disables retention. All shared handles
    /// see the same limit. In-flight decodes recheck it under the same lock
    /// before insertion, so a concurrent shrink cannot violate the new bound.
    pub fn set_byte_budget(&self, max_bytes: usize) {
        let mut inner = self.lock_inner();
        inner.max_bytes = max_bytes;
        if inner.resident_bytes > max_bytes || max_bytes == 0 {
            inner.entries.clear();
            inner.resident_bytes = 0;
        }
    }

    /// Release retained decoded coefficients, without closing the image or resetting
    /// counters. Immutable directory metadata remains; in-flight coefficient readers
    /// keep their own references until they finish and decodes may repopulate
    /// entries. This is not an image reset: never reuse the cache with another
    /// container after clearing.
    pub fn clear(&self) {
        let mut inner = self.lock_inner();
        inner.entries.clear();
        inner.resident_bytes = 0;
    }

    /// `(hits, misses, resident_entries)` for measurement.
    pub fn stats(&self) -> (u64, u64, usize) {
        let inner = self.lock_inner();
        (inner.hits, inner.misses, inner.entries.len())
    }

    /// Lock the memo, recovering the guard if a foreign panic
    /// poisoned it. Poisoning from this module is unreachable: the
    /// guarded sections never run caller code. Optional immutable-image compaction
    /// uses only validated internal storage; recovery keeps a live provider usable.
    fn lock_inner(&self) -> MutexGuard<'_, ChunkCacheInner> {
        self.inner
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
    }

    /// Return the memoized words for `key`, decoding (outside the
    /// lock) and memoizing them on miss. Decoding failures propagate
    /// and memoize nothing; the miss is still counted.
    #[cfg(test)]
    fn get_or_decode(
        &self,
        key: (i32, u32),
        decode: impl FnOnce() -> Result<Vec<f64>>,
    ) -> Result<Arc<TrackedVec<f64>>> {
        self.get_or_decode_tracked(key, || {
            Ok(TrackedVec {
                values: decode()?,
                _allocation: None,
            })
        })
    }

    #[cfg(test)]
    fn get_or_decode_tracked(
        &self,
        key: (i32, u32),
        decode: impl FnOnce() -> Result<TrackedVec<f64>>,
    ) -> Result<Arc<TrackedVec<f64>>> {
        self.get_or_decode_key(CacheKey::Body(key.0, key.1), decode)
    }

    pub(super) fn get_or_decode_key(
        &self,
        key: CacheKey,
        decode: impl FnOnce() -> Result<TrackedVec<f64>>,
    ) -> Result<Arc<TrackedVec<f64>>> {
        match self.get_or_decode_storage(key, None, decode)? {
            CachedStorage::Dense(words) => Ok(words),
            CachedStorage::Packed(_) => unreachable!("raw cache entry remains dense"),
        }
    }

    fn get_or_decode_storage(
        &self,
        key: CacheKey,
        index: Option<&FileIndex>,
        decode: impl FnOnce() -> Result<TrackedVec<f64>>,
    ) -> Result<CachedStorage> {
        if let Some(hit) = self.hit(&key) {
            if index.is_some() || matches!(hit, CachedStorage::Dense(_)) {
                return Ok(hit);
            }
        }
        // No coefficient-table lock is held while waiting or decoding. Raw
        // slice caches retain their independent parallel decode behavior.
        let _decode_guard = self.resources.as_ref().map(|resources| {
            resources
                .decode_gate
                .lock()
                .unwrap_or_else(|error| error.into_inner())
        });
        if self.resources.is_some() {
            if let Some(hit) = self.hit(&key) {
                if index.is_some() || matches!(hit, CachedStorage::Dense(_)) {
                    return Ok(hit);
                }
            }
        }
        {
            let mut inner = self.lock_inner();
            inner.misses = inner.misses.saturating_add(1);
        }
        let words = decode()?;
        let size = words.capacity().saturating_mul(std::mem::size_of::<f64>());
        let shared = CachedStorage::Dense(Arc::new(words));
        let mut inner = self.lock_inner();
        if inner.max_bytes == 0
            || size > inner.max_bytes
            || (matches!(key, CacheKey::Window { .. }) && size > MAX_COEFFICIENT_WINDOW_BYTES)
        {
            return Ok(shared);
        }
        // Another reader may have decoded the same immutable chunk while we
        // were outside the lock. Do not count/retain its allocation twice.
        if let Some(existing) = inner.get(&key) {
            // Raw/Table callers keep a complete dense owner even if a private
            // native body callback previously compacted this immutable key.
            return Ok(
                if index.is_none() && matches!(existing, CachedStorage::Packed(_)) {
                    shared
                } else {
                    existing
                },
            );
        }
        if inner.entries.len() < MAX_CACHED_CHUNKS
            && size > inner.max_bytes.saturating_sub(inner.resident_bytes)
        {
            if let Some(index) = index {
                self.compact_first_victim(&mut inner, index, size);
            }
        }
        while inner.entries.len() >= MAX_CACHED_CHUNKS
            || size > inner.max_bytes.saturating_sub(inner.resident_bytes)
        {
            inner.evict_oldest();
        }
        let retention = if let Some(resources) = &self.resources {
            loop {
                match resources.retained.reserve(size) {
                    Ok(lease) => break Some(lease),
                    Err(_) if inner.evict_oldest() => continue,
                    // Cache retention is optional. Live words remain
                    // accounted even if another image occupies the cache.
                    Err(_) => return Ok(shared),
                }
            }
        } else {
            None
        };
        inner.resident_bytes += size;
        let last_access = inner.next_access();
        inner.entries.insert(
            key,
            CachedWords {
                storage: shared.clone(),
                last_access,
                _retention: retention,
            },
        );
        Ok(shared)
    }

    /// Hit-path probe: clones the resident `Arc` and counts the hit
    /// under one guard. Returns `None` on a miss (counted by the
    /// caller in [`LebChunkCache::get_or_decode_storage`]).
    fn hit(&self, key: &CacheKey) -> Option<CachedStorage> {
        let mut inner = self.lock_inner();
        let hit = inner.get(key)?;
        inner.hits = inner.hits.saturating_add(1);
        Some(hit)
    }
}

#[cfg(test)]
#[path = "tests/cache.rs"]
mod cache_tests;

#[cfg(test)]
#[path = "tests/window_cache.rs"]
mod window_cache_tests;

#[cfg(test)]
mod word_image_identity_tests {
    use super::*;

    #[test]
    fn public_cache_constructors_do_not_license_arbitrary_index_source_pairs() {
        assert!(LebChunkCache::new().word_image_identity().is_none());
        assert!(
            LebChunkCache::with_byte_budget(0)
                .word_image_identity()
                .is_none()
        );
        assert!(
            LebChunkCache::with_byte_budget(usize::MAX)
                .word_image_identity()
                .is_none()
        );
    }

    #[test]
    fn image_cache_tokens_survive_clone_move_and_try_unwrap_without_address_keys() {
        let first = Arc::new(LebChunkCache::with_image_byte_budget(0));
        let identity = first.word_image_identity().unwrap();
        let clone = Arc::clone(&first);
        assert_eq!(clone.word_image_identity(), Some(identity));
        drop(first);
        let moved = Arc::try_unwrap(clone).unwrap();
        assert_eq!(moved.word_image_identity(), Some(identity));
        assert_ne!(
            LebChunkCache::with_image_byte_budget(0).word_image_identity(),
            Some(identity)
        );
        assert_eq!(std::mem::size_of::<Option<std::num::NonZeroU64>>(), 8);
        println!(
            "LebChunkCache bytes with identity: {}; explicit per-cache identity bytes: {}",
            std::mem::size_of::<LebChunkCache>(),
            std::mem::size_of::<Option<std::num::NonZeroU64>>()
        );
    }
}

#[cfg(test)]
#[path = "tests/scalar_split_input_keys.rs"]
mod scalar_split_input_keys_tests;

impl Default for LebChunkCache {
    /// Empty memo (same as [`LebChunkCache::new`]).
    fn default() -> Self {
        Self::new()
    }
}

/// Resolve a physical body table to one entry per identifier.
///
/// Later records replace earlier values in place, retaining each identifier's
/// first insertion position. The common header still records the physical
/// entry count; evaluation and metadata views consume this unique catalog.
fn resolve_body_entries(entries: Vec<BodyEntry>) -> Vec<BodyEntry> {
    let mut slots = HashMap::new();
    // The validated physical count bounds the final unique table as well.
    // Exact reservation makes allocation admission independent of Vec growth.
    let mut bodies = Vec::with_capacity(entries.len());
    for entry in entries {
        let slot = *slots.entry(entry.body_id()).or_insert(bodies.len());
        if slot == bodies.len() {
            bodies.push(entry);
        } else {
            bodies[slot] = entry;
        }
    }
    bodies
}

impl FileIndex {
    /// Decode a complete container image: header, directory, then the
    /// body index in section 0 guarded by [`validate_entry_count`]. A
    /// missing section 0 is [`Error::LebCorrupt`], even for a zero-body image.
    /// Every physical body record is decoded and validated. Repeated body ids
    /// retain their last record at the first id's insertion position.
    pub fn decode(bytes: &[u8]) -> Result<Self> {
        let header = FileHeader::decode(bytes)?;
        let sections = decode_section_table(bytes, header.section_count, bytes.len())?;
        let body_section = sections
            .iter()
            .find(|s| s.section_id == SECTION_BODY_INDEX)
            .ok_or_else(|| Error::LebCorrupt {
                what: "LEB body index".to_string(),
                offset: HEADER_SIZE,
            })?;
        validate_entry_count(
            u64::from(header.body_count),
            header.format.body_entry_size() as u64,
            body_section.size,
            "LEB body index",
            body_section.offset as usize,
        )?;
        let mut bodies = Vec::new();
        for i in 0..header.body_count {
            let at = (body_section.offset as usize)
                .checked_add((i as usize) * header.format.body_entry_size())
                .ok_or_else(|| Error::LebCorrupt {
                    what: "LEB body index".to_string(),
                    offset: body_section.offset as usize,
                })?;
            bodies.push(BodyEntry::decode(bytes, at, bytes.len(), header.format)?);
        }
        Ok(Self {
            header,
            sections,
            bodies: resolve_body_entries(bodies),
        })
    }

    /// Inspect only the header, directory and body entries on an open regular
    /// file. Coefficient payloads and auxiliary tables are not read. Uses the
    /// same field decoders/bounds as [`Self::decode`], with no whole-file buffer.
    /// Changes the descriptor cursor; the caller must keep data stable while
    /// inspecting if it needs a coherent inventory (this does not make a snapshot).
    /// Repeated body ids use the same last-record and first-position resolution
    /// as [`Self::decode`]; header counts continue to describe physical records.
    pub fn read_from_file(file: &mut std::fs::File) -> Result<Self> {
        use std::io::{Read, Seek, SeekFrom};
        let io_error = |error: std::io::Error| Error::LebCorrupt {
            what: format!("LEB index read: {error}"),
            offset: 0,
        };
        let metadata = file.metadata().map_err(io_error)?;
        if !metadata.is_file() {
            return Err(Error::LebCorrupt {
                what: "LEB index regular file".into(),
                offset: 0,
            });
        }
        let file_len = usize::try_from(metadata.len()).map_err(|_| Error::LebCorrupt {
            what: "LEB index file length".into(),
            offset: 0,
        })?;
        file.seek(SeekFrom::Start(0)).map_err(io_error)?;
        let mut prefix = [0; HEADER_SIZE];
        file.read_exact(&mut prefix).map_err(io_error)?;
        let header = FileHeader::decode(&prefix)?;
        validate_entry_count(
            u64::from(header.section_count),
            SECTION_DIR_SIZE as u64,
            file_len.saturating_sub(HEADER_SIZE) as u64,
            "LEB section directory",
            HEADER_SIZE,
        )?;
        let mut sections = Vec::with_capacity(header.section_count as usize);
        let mut section_buffer = [0; SECTION_DIR_SIZE];
        for i in 0..header.section_count {
            file.read_exact(&mut section_buffer).map_err(io_error)?;
            let at = HEADER_SIZE + i as usize * SECTION_DIR_SIZE;
            let entry = SectionEntry::decode(&section_buffer, 0).map_err(|error| match error {
                Error::LebCorrupt { what, offset } => Error::LebCorrupt {
                    what,
                    offset: at.saturating_add(offset),
                },
                other => other,
            })?;
            super::decode::insert_checked_section(&mut sections, entry, at, file_len as u64)?;
        }
        let body_section = sections
            .iter()
            .find(|s| s.section_id == SECTION_BODY_INDEX)
            .ok_or_else(|| Error::LebCorrupt {
                what: "LEB body index".into(),
                offset: HEADER_SIZE,
            })?;
        let stride = header.format.body_entry_size();
        validate_entry_count(
            u64::from(header.body_count),
            stride as u64,
            body_section.size,
            "LEB body index",
            body_section.offset as usize,
        )?;
        file.seek(SeekFrom::Start(body_section.offset))
            .map_err(io_error)?;
        let mut buffer = vec![0; stride];
        let mut bodies = Vec::with_capacity(header.body_count as usize);
        for i in 0..header.body_count {
            file.read_exact(&mut buffer).map_err(io_error)?;
            let at = body_section.offset as usize + i as usize * stride;
            let entry = BodyEntry::decode(&buffer, 0, file_len, header.format).map_err(
                |error| match error {
                    Error::LebCorrupt { what, offset } => Error::LebCorrupt {
                        what,
                        offset: at.saturating_add(offset),
                    },
                    other => other,
                },
            )?;
            bodies.push(entry);
        }
        Ok(Self {
            header,
            sections,
            bodies: resolve_body_entries(bodies),
        })
    }

    /// Decoded common header.
    pub fn header(&self) -> &FileHeader {
        &self.header
    }
    /// Decoded section directory (duplicate ids resolved last-wins).
    pub fn sections(&self) -> &[SectionEntry] {
        &self.sections
    }
    /// Unique body index entries in first-id insertion order, using the last
    /// physical record for each repeated id. [`FileHeader::body_count`] retains
    /// the physical table count and can therefore exceed this slice's length.
    pub fn bodies(&self) -> &[BodyEntry] {
        &self.bodies
    }
    /// Section payload by id, if present.
    pub fn section(&self, id: u32) -> Option<&SectionEntry> {
        self.sections.iter().find(|s| s.section_id == id)
    }

    /// True when `body_id` has an index entry in this file.
    pub fn has_body(&self, body_id: i32) -> bool {
        self.bodies.iter().any(|b| b.body_id == body_id)
    }

    /// Stored per-body JD coverage, or `None` when the body is absent.
    /// Per-body ranges own evaluation; the file-level range never does.
    /// A repeated id's last physical record owns its range.
    pub fn body_coverage(&self, body_id: i32) -> Option<(f64, f64)> {
        self.bodies
            .iter()
            .find(|b| b.body_id == body_id)
            .map(|b| (b.jd_start, b.jd_end))
    }

    /// Last physical index entry for `body_id`, or [`Error::LebBodyNotFound`].
    pub fn find_body(&self, body_id: i32) -> Result<&BodyEntry> {
        let Some(entry) = self.bodies.iter().find(|b| b.body_id == body_id) else {
            return Err(Error::LebBodyNotFound { body_id });
        };
        Ok(entry)
    }

    /// Containing polynomial's TT bounds for a body, or `None` if unserved.
    /// Interior joins follow the evaluation selector's right ownership.
    /// Invalid segment metadata returns its native reader error.
    pub fn body_smooth_interval(&self, body_id: i32, jd: f64) -> Result<Option<(f64, f64)>> {
        match self.find_body(body_id) {
            Ok(entry) => entry.smooth_interval(jd),
            Err(Error::LebBodyNotFound { .. }) => Ok(None),
            Err(error) => Err(error),
        }
    }

    /// Prove scalar/split argument identity using this actual selected record.
    pub(crate) fn scalar_state_matches_split(
        &self,
        body: i32,
        jd: f64,
        epoch: TwoPartEpoch,
    ) -> bool {
        self.find_body(body)
            .is_ok_and(|entry| scalar_and_split_share_argument(entry, jd, epoch))
    }

    /// Evaluate three stored ICRS components without heap result buffers. The
    /// ordinary scalar path still owns format/domain/record/error decisions.
    pub(crate) fn eval_state_inline_source(
        &self,
        bytes: LebSource<'_>,
        body: i32,
        jd: f64,
        cache: &LebChunkCache,
    ) -> Result<Option<CartesianState>> {
        let entry = self.find_body(body)?;
        if entry.components() != 3
            || !matches!(entry.coord_type(), COORD_ICRS_BARY | COORD_ICRS_BARY_SYSTEM)
        {
            return Ok(None);
        }
        let state = if self.header.format == LebFormat::Leb1 {
            self.eval_body_source_with(bytes, body, jd, eval_one_segment_inline)
        } else {
            self.eval_leb2_body_cached_source_with(bytes, body, jd, cache, eval_one_segment_inline)
        }?;
        Ok(Some(state))
    }

    /// Split-epoch form of the same stored-component operation. Both the
    /// collapsed-word record ownership and the two-part normalized time remain
    /// those of the existing reader path.
    pub(crate) fn eval_state_inline_two_part_source(
        &self,
        bytes: LebSource<'_>,
        body: i32,
        epoch: TwoPartEpoch,
        cache: &LebChunkCache,
    ) -> Result<Option<CartesianState>> {
        let entry = self.find_body(body)?;
        if entry.components() != 3
            || !matches!(entry.coord_type(), COORD_ICRS_BARY | COORD_ICRS_BARY_SYSTEM)
        {
            return Ok(None);
        }
        let state = if self.header.format == LebFormat::Leb1 {
            self.eval_body_two_part_source_with(bytes, body, epoch, eval_one_segment_inline)
        } else {
            self.eval_leb2_body_two_part_cached_source_with(
                bytes,
                body,
                epoch,
                cache,
                eval_one_segment_inline,
            )
        }?;
        Ok(Some(state))
    }

    /// Evaluate a monolithic LEB1 body at Julian Day TT.
    ///
    /// The selected Chebyshev segment supplies position and analytic velocity;
    /// the latter is scaled by `2/interval_days`. ICRS components are AU and
    /// AU/day. Ecliptic components are degrees/degrees/AU and their daily rates;
    /// longitude is normalized to `[0,360)`.
    ///
    /// Missing bodies, epochs outside the body's stored interval, malformed
    /// coefficient ranges and compressed payloads on this raw path return typed
    /// errors. No polynomial is extrapolated beyond its declared body domain.
    pub(crate) fn eval_body_source(
        &self,
        bytes: LebSource<'_>,
        body_id: i32,
        jd: f64,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        self.eval_body_source_with(bytes, body_id, jd, eval_one_segment)
    }

    fn eval_body_source_with<R>(
        &self,
        bytes: LebSource<'_>,
        body_id: i32,
        jd: f64,
        evaluate: impl FnOnce(&BodyEntry, &[f64], f64) -> Result<R>,
    ) -> Result<R> {
        if self.header.format != LebFormat::Leb1 {
            return Err(Error::LebCorrupt {
                what: "LEB2 compressed body (evaluation over decoded blobs is LEB-06 work)"
                    .to_string(),
                offset: 0,
            });
        }
        let entry = self.find_body(body_id)?;
        if !jd.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd });
        }
        if jd < entry.jd_start || jd > entry.jd_end {
            return Err(Error::LebOutOfRange {
                body_id,
                jd,
                start: entry.jd_start,
                end: entry.jd_end,
            });
        }
        // Defensive recheck: decode already rejects these, but evaluation
        // must never divide by zero or index an empty table.
        if entry.segment_count == 0
            || !entry.interval_days.is_finite()
            || entry.interval_days <= 0.0
        {
            return Err(Error::LebCorrupt {
                what: format!("LEB body index (body {body_id})"),
                offset: usize::try_from(entry.data_offset).unwrap_or(usize::MAX),
            });
        }
        let (seg_idx, tau_norm) =
            select_segment(entry.jd_start, entry.interval_days, entry.segment_count, jd)?;
        let coefficients = read_segment_coefficients(bytes, entry, seg_idx)?;
        evaluate(entry, &coefficients, tau_norm.value())
    }

    /// Evaluate a body's position and velocity at the two-part epoch
    /// `epoch` (TT) from a monolithic (LEB1 raw-coefficient) container
    /// image — the split-epoch form of [`FileIndex::eval_body`].
    ///
    /// The release has no two-part LEB API, so this form is Rust-native
    /// (like the timescale two-part scale converters): segment ownership is defined
    /// by the collapsed word, hence bitwise-identical to the scalar path
    /// on `epoch.value` by construction, while `tau` keeps the
    /// sub-ulp tail via [`select_segment_two_part`]. Errors mirror the
    /// scalar path; a directly constructed non-finite part is
    /// [`Error::NonFiniteArgument`].
    pub(crate) fn eval_body_two_part_source(
        &self,
        bytes: LebSource<'_>,
        body_id: i32,
        epoch: TwoPartEpoch,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        self.eval_body_two_part_source_with(bytes, body_id, epoch, eval_one_segment)
    }

    fn eval_body_two_part_source_with<R>(
        &self,
        bytes: LebSource<'_>,
        body_id: i32,
        epoch: TwoPartEpoch,
        evaluate: impl FnOnce(&BodyEntry, &[f64], f64) -> Result<R>,
    ) -> Result<R> {
        if self.header.format != LebFormat::Leb1 {
            return Err(Error::LebCorrupt {
                what: "LEB2 compressed body (evaluation over decoded blobs is LEB-06 work)"
                    .to_string(),
                offset: 0,
            });
        }
        let entry = self.find_body(body_id)?;
        let jd = epoch.value();
        if !jd.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd });
        }
        if jd < entry.jd_start || jd > entry.jd_end {
            return Err(Error::LebOutOfRange {
                body_id,
                jd,
                start: entry.jd_start,
                end: entry.jd_end,
            });
        }
        // Defensive recheck: decode already rejects these, but evaluation
        // must never divide by zero or index an empty table.
        if entry.segment_count == 0
            || !entry.interval_days.is_finite()
            || entry.interval_days <= 0.0
        {
            return Err(Error::LebCorrupt {
                what: format!("LEB body index (body {body_id})"),
                offset: usize::try_from(entry.data_offset).unwrap_or(usize::MAX),
            });
        }
        let (seg_idx, tau_norm) = select_segment_two_part(
            entry.jd_start,
            entry.interval_days,
            entry.segment_count,
            epoch,
        )?;
        let coefficients = read_segment_coefficients(bytes, entry, seg_idx)?;
        evaluate(entry, &coefficients, tau_norm.value())
    }

    /// Evaluate a body's position and velocity at each Julian Day in
    /// `jds` (TT) from a monolithic (LEB1 raw-coefficient) container
    /// image — the multi-epoch form of [`FileIndex::eval_body`].
    ///
    /// for the loop contract: a
    /// vector time is served by one scalar evaluation per element, in
    /// input order, through the same evaluator (same segment ownership,
    /// `tau` map, and errors — the release lets the first miss raise).
    /// The Rust form returns the states in input order and fails fast
    /// with the first element's scalar error; an empty input yields an
    /// empty vector. No new numerical path is introduced here.
    pub fn eval_body_many(
        &self,
        bytes: &[u8],
        body_id: i32,
        jds: &[f64],
    ) -> Result<Vec<(Vec<f64>, Vec<f64>)>> {
        batch_eval(jds, |jd| self.eval_body(bytes, body_id, jd))
    }

    /// Multi-epoch form of [`FileIndex::eval_body_two_part`]: one
    /// split-epoch scalar evaluation per element, in input order,
    /// failing fast with the first element's scalar error.
    pub fn eval_body_two_part_many(
        &self,
        bytes: &[u8],
        body_id: i32,
        epochs: &[TwoPartEpoch],
    ) -> Result<Vec<(Vec<f64>, Vec<f64>)>> {
        batch_eval_epochs(epochs, |epoch| {
            self.eval_body_two_part(bytes, body_id, epoch)
        })
    }

    /// Evaluate a body's position and velocity at Julian Day `jd` (TT)
    /// from a compressed (LEB2 v1 monolithic or v2 chunked) container
    /// image, decoding the owning blob through the chunked-blob decode pipeline.
    ///
    /// Range check, corruption
    /// guard, global segment selection ([`select_segment`]), chunk hint
    /// ([`find_chunk`]) with the segment-wins ownership rule
    /// ([`select_chunk`]), `tau` map + clamp, per-component Clenshaw
    /// value/derivative with the `2 / interval_days` velocity scale,
    /// and ecliptic `pos[0] % 360.0` wrap are all the release order —
    /// the only structural difference is that decompression is a pure
    /// function of the caller's byte slice (no memo cache; the file-lifecycle layer owns
    /// the lifecycle) rather than a lock-guarded method cache.
    ///
    /// Errors mirror [`FileIndex::eval_body`]: [`Error::LebBodyNotFound`],
    /// [`Error::LebOutOfRange`], [`Error::LebCorrupt`] (empty chunk
    /// index, uncovered segment, stored-vs-derived size mismatch, zstd
    /// failure, blob truncation), [`Error::NonFiniteArgument`]. Calling
    /// this on an LEB1 image is [`Error::LebCorrupt`] (raw coefficients
    /// belong to [`FileIndex::eval_body`]).
    pub(crate) fn eval_leb2_body_source(
        &self,
        bytes: LebSource<'_>,
        body_id: i32,
        jd: f64,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        if self.header.format == LebFormat::Leb1 {
            return Err(Error::LebCorrupt {
                what: "LEB1 raw body (use eval_body for monolithic coefficients)".to_string(),
                offset: 0,
            });
        }
        let entry = self.find_body(body_id)?;
        if !jd.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd });
        }
        if jd < entry.jd_start || jd > entry.jd_end {
            return Err(Error::LebOutOfRange {
                body_id,
                jd,
                start: entry.jd_start,
                end: entry.jd_end,
            });
        }
        if entry.segment_count == 0
            || !entry.interval_days.is_finite()
            || entry.interval_days <= 0.0
        {
            return Err(Error::LebCorrupt {
                what: format!("LEB2 body index (body {body_id})"),
                offset: usize::try_from(entry.data_offset).unwrap_or(usize::MAX),
            });
        }
        let (seg_idx, tau_norm) =
            select_segment(entry.jd_start, entry.interval_days, entry.segment_count, jd)?;
        let tau = tau_norm.value();
        self.leb2_eval_at(bytes, body_id, entry, seg_idx, jd, tau, None)
    }

    /// Evaluate a body's position and velocity at Julian Day `jd` (TT)
    /// serving immutable chunk directories and the chunk-blob decode from `cache`.
    ///
    /// Bitwise-identical to [`FileIndex::eval_leb2_body`]: guards,
    /// segment ownership, `tau` and the Clenshaw tail are the scalar
    /// path's; only immutable directories and blob decode consult `cache` (see
    /// [`LebChunkCache`] for the key soundness and bound). Errors
    /// mirror [`FileIndex::eval_leb2_body`].
    pub(crate) fn eval_leb2_body_cached_source(
        &self,
        bytes: LebSource<'_>,
        body_id: i32,
        jd: f64,
        cache: &LebChunkCache,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        self.eval_leb2_body_cached_source_with(bytes, body_id, jd, cache, eval_one_segment)
    }

    fn eval_leb2_body_cached_source_with<R>(
        &self,
        bytes: LebSource<'_>,
        body_id: i32,
        jd: f64,
        cache: &LebChunkCache,
        evaluate: impl FnOnce(&BodyEntry, &[f64], f64) -> Result<R>,
    ) -> Result<R> {
        if self.header.format == LebFormat::Leb1 {
            return Err(Error::LebCorrupt {
                what: "LEB1 raw body (use eval_body for monolithic coefficients)".to_string(),
                offset: 0,
            });
        }
        let entry = self.find_body(body_id)?;
        if !jd.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd });
        }
        if jd < entry.jd_start || jd > entry.jd_end {
            return Err(Error::LebOutOfRange {
                body_id,
                jd,
                start: entry.jd_start,
                end: entry.jd_end,
            });
        }
        if entry.segment_count == 0
            || !entry.interval_days.is_finite()
            || entry.interval_days <= 0.0
        {
            return Err(Error::LebCorrupt {
                what: format!("LEB2 body index (body {body_id})"),
                offset: usize::try_from(entry.data_offset).unwrap_or(usize::MAX),
            });
        }
        let (seg_idx, tau_norm) =
            select_segment(entry.jd_start, entry.interval_days, entry.segment_count, jd)?;
        let tau = tau_norm.value();
        self.with_leb2_segment_coefficients(
            bytes,
            body_id,
            entry,
            seg_idx,
            jd,
            Some(cache),
            |coefficients| evaluate(entry, coefficients, tau),
        )
    }

    /// Evaluate a body's position and velocity at the two-part epoch
    /// `epoch` (TT) from a compressed (LEB2 v1 monolithic or v2 chunked)
    /// container image — the split-epoch form of
    /// [`FileIndex::eval_leb2_body`].
    ///
    /// Segment ownership is defined by the collapsed word (identical to
    /// the scalar path on `epoch.value` by construction, including the
    /// segment-wins chunk rule, whose hint also keys on the collapsed
    /// word); `tau` keeps the sub-ulp tail via
    /// [`select_segment_two_part`]. Errors mirror the scalar path.
    pub(crate) fn eval_leb2_body_two_part_source(
        &self,
        bytes: LebSource<'_>,
        body_id: i32,
        epoch: TwoPartEpoch,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        if self.header.format == LebFormat::Leb1 {
            return Err(Error::LebCorrupt {
                what: "LEB1 raw body (use eval_body for monolithic coefficients)".to_string(),
                offset: 0,
            });
        }
        let entry = self.find_body(body_id)?;
        let jd = epoch.value();
        if !jd.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd });
        }
        if jd < entry.jd_start || jd > entry.jd_end {
            return Err(Error::LebOutOfRange {
                body_id,
                jd,
                start: entry.jd_start,
                end: entry.jd_end,
            });
        }
        if entry.segment_count == 0
            || !entry.interval_days.is_finite()
            || entry.interval_days <= 0.0
        {
            return Err(Error::LebCorrupt {
                what: format!("LEB2 body index (body {body_id})"),
                offset: usize::try_from(entry.data_offset).unwrap_or(usize::MAX),
            });
        }
        let (seg_idx, tau_norm) = select_segment_two_part(
            entry.jd_start,
            entry.interval_days,
            entry.segment_count,
            epoch,
        )?;
        let tau = tau_norm.value();
        self.leb2_eval_at(bytes, body_id, entry, seg_idx, jd, tau, None)
    }

    /// Evaluate a body's position and velocity at the two-part epoch
    /// `epoch` (TT) serving the chunk-blob decode from `cache`.
    ///
    /// Bitwise-identical to [`FileIndex::eval_leb2_body_two_part`]:
    /// segment ownership is defined by the collapsed word and `tau`
    /// keeps the sub-ulp tail; only the blob decode consults `cache`
    /// (see [`LebChunkCache`] for the key soundness and bound). Errors
    /// mirror [`FileIndex::eval_leb2_body_two_part`].
    pub(crate) fn eval_leb2_body_two_part_cached_source(
        &self,
        bytes: LebSource<'_>,
        body_id: i32,
        epoch: TwoPartEpoch,
        cache: &LebChunkCache,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        self.eval_leb2_body_two_part_cached_source_with(
            bytes,
            body_id,
            epoch,
            cache,
            eval_one_segment,
        )
    }

    fn eval_leb2_body_two_part_cached_source_with<R>(
        &self,
        bytes: LebSource<'_>,
        body_id: i32,
        epoch: TwoPartEpoch,
        cache: &LebChunkCache,
        evaluate: impl FnOnce(&BodyEntry, &[f64], f64) -> Result<R>,
    ) -> Result<R> {
        let (entry, jd, seg_idx, tau) = self.leb2_two_part_record(body_id, epoch)?;
        self.with_leb2_segment_coefficients(
            bytes,
            body_id,
            entry,
            seg_idx,
            jd,
            Some(cache),
            |coefficients| evaluate(entry, coefficients, tau),
        )
    }

    /// Ordinary split-input admission shared by row acquisition and its native
    /// borrowed-row consumer. Keep the scalar normalization independent of the
    /// Word recurrence; both select the same immutable record from these words.
    fn leb2_two_part_record(
        &self,
        body_id: i32,
        epoch: TwoPartEpoch,
    ) -> Result<(&BodyEntry, f64, usize, f64)> {
        if self.header.format == LebFormat::Leb1 {
            return Err(Error::LebCorrupt {
                what: "LEB1 raw body (use eval_body for monolithic coefficients)".to_string(),
                offset: 0,
            });
        }
        let entry = self.find_body(body_id)?;
        let jd = epoch.value();
        if !jd.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd });
        }
        if jd < entry.jd_start || jd > entry.jd_end {
            return Err(Error::LebOutOfRange {
                body_id,
                jd,
                start: entry.jd_start,
                end: entry.jd_end,
            });
        }
        if entry.segment_count == 0
            || !entry.interval_days.is_finite()
            || entry.interval_days <= 0.0
        {
            return Err(Error::LebCorrupt {
                what: format!("LEB2 body index (body {body_id})"),
                offset: usize::try_from(entry.data_offset).unwrap_or(usize::MAX),
            });
        }
        let (seg_idx, tau_norm) = select_segment_two_part(
            entry.jd_start,
            entry.interval_days,
            entry.segment_count,
            epoch,
        )?;
        let tau = tau_norm.value();
        Ok((entry, jd, seg_idx, tau))
    }

    /// Evaluate the ordinary scalar state at its original native call position.
    /// The row must come from this index/body and these exact split words through
    /// `with_native_word_parts_two_part_source`. That private synchronous borrow
    /// preserves its selected coefficients without retaining another owner.
    pub(crate) fn eval_borrowed_state_inline_two_part(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
        coefficients: &[f64],
    ) -> Result<Option<CartesianState>> {
        let entry = self.find_body(body)?;
        if entry.components() != 3
            || !matches!(entry.coord_type(), COORD_ICRS_BARY | COORD_ICRS_BARY_SYSTEM)
        {
            return Ok(None);
        }
        let (entry, _, _, tau) = self.leb2_two_part_record(body, epoch)?;
        eval_one_segment_inline(entry, coefficients, tau).map(Some)
    }

    /// Evaluate a body's position and velocity at each Julian Day in
    /// `jds` (TT) from a compressed (LEB2 v1 monolithic or v2 chunked)
    /// container image — the multi-epoch form of
    /// [`FileIndex::eval_leb2_body`]. Same per-element scalar loop and
    /// fail-fast contract as [`FileIndex::eval_body_many`].
    pub fn eval_leb2_body_many(
        &self,
        bytes: &[u8],
        body_id: i32,
        jds: &[f64],
    ) -> Result<Vec<(Vec<f64>, Vec<f64>)>> {
        batch_eval(jds, |jd| self.eval_leb2_body(bytes, body_id, jd))
    }

    /// Multi-epoch form of [`FileIndex::eval_leb2_body_two_part`]: one
    /// split-epoch scalar evaluation per element, in input order,
    /// failing fast with the first element's scalar error.
    pub fn eval_leb2_body_two_part_many(
        &self,
        bytes: &[u8],
        body_id: i32,
        epochs: &[TwoPartEpoch],
    ) -> Result<Vec<(Vec<f64>, Vec<f64>)>> {
        batch_eval_epochs(epochs, |epoch| {
            self.eval_leb2_body_two_part(bytes, body_id, epoch)
        })
    }

    /// Evaluate one selected LEB2 segment at normalized epoch `tau`:
    /// monolith/chunk blob decode plus the shared Chebyshev tail, with
    /// the chunk-search hint `hint_jd`.
    /// Used by both the scalar ([`FileIndex::eval_leb2_body`]) and the
    /// two-part ([`FileIndex::eval_leb2_body_two_part`]) paths, and by
    /// their cached twins.
    ///
    /// `cache` restores the release-side decoded-blob memo without
    /// changing any word: with `Some`, the blob decode (zstd, reorder,
    /// `f64` conversion) is served from the caller's
    /// [`LebChunkCache`] keyed `(body_id, chunk_index)` — `u32::MAX`
    /// for LEB2 v1 monoliths — while chunk ownership, segment slicing
    /// and the Clenshaw tail run unchanged. A chunk index that cannot
    /// key the memo falls back to the direct decode, never a new
    /// failure mode. `None` keeps the historical direct decode.
    #[allow(clippy::too_many_arguments)]
    fn leb2_eval_at(
        &self,
        bytes: LebSource<'_>,
        body_id: i32,
        entry: &BodyEntry,
        seg_idx: usize,
        hint_jd: f64,
        tau: f64,
        cache: Option<&LebChunkCache>,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        self.with_leb2_segment_coefficients(
            bytes,
            body_id,
            entry,
            seg_idx,
            hint_jd,
            cache,
            |coefficients| eval_one_segment(entry, coefficients, tau),
        )
    }

    /// Decode the actual selected polynomial, sharing the ordinary image cache.
    pub(super) fn leb2_segment_coefficients(
        &self,
        bytes: LebSource<'_>,
        body_id: i32,
        entry: &BodyEntry,
        seg_idx: usize,
        hint_jd: f64,
        cache: Option<&LebChunkCache>,
    ) -> Result<CoefficientWords> {
        self.with_leb2_segment_coefficients(
            bytes,
            body_id,
            entry,
            seg_idx,
            hint_jd,
            cache,
            |coefficients| bytes.copy_words(coefficients),
        )
    }

    /// Consume the selected row while its immutable decoded owner is retained.
    /// The callback cannot return a borrow tied to this temporary row. Ordinary
    /// scalar evaluation avoids a row copy; relative polynomial operations use
    /// the owned projection above when coefficients must outlive this scope.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn with_leb2_segment_coefficients<R>(
        &self,
        bytes: LebSource<'_>,
        body_id: i32,
        entry: &BodyEntry,
        seg_idx: usize,
        hint_jd: f64,
        cache: Option<&LebChunkCache>,
        consume: impl FnOnce(&[f64]) -> Result<R>,
    ) -> Result<R> {
        let seg_u32 = u32::try_from(seg_idx).map_err(|_| Error::LebCorrupt {
            what: format!("LEB2 chunk index (body {body_id})"),
            offset: 0,
        })?;
        match self.header.format {
            LebFormat::Leb1 => unreachable!("gated above"),
            LebFormat::Leb2V1 => {
                let decode = || {
                    let what = format!("LEB2 body entry {body_id}");
                    let size = entry.compressed_size().ok_or_else(|| Error::LebCorrupt {
                        what: what.clone(),
                        offset: 0,
                    })?;
                    bytes.check_frame(entry.expected_uncompressed_size()?)?;
                    let allocation =
                        bytes.reserve_decode(entry.expected_uncompressed_size()?, size)?;
                    let compressed = bytes.blob(entry.data_offset, size, &what)?;
                    let values = reader::decompress_body_words_limited(
                        &compressed,
                        entry.uncompressed_size().ok_or_else(|| Error::LebCorrupt {
                            what: what.clone(),
                            offset: 0,
                        })?,
                        entry.segment_count(),
                        entry.degree(),
                        entry.components(),
                        bytes.window_limit(),
                    )?;
                    let capacity = values.capacity() * std::mem::size_of::<f64>();
                    drop(compressed);
                    let allocation = allocation.map(|lease| lease.retain(capacity)).transpose()?;
                    Ok(TrackedVec {
                        values,
                        _allocation: allocation,
                    })
                };
                match cache {
                    Some(memo) => memo.with_cached_body_segment(
                        (body_id, u32::MAX),
                        self,
                        entry,
                        seg_idx,
                        entry.segment_count(),
                        decode,
                        consume,
                    ),
                    None => consume(slice_one_segment(&decode()?, entry, seg_idx, body_id)?),
                }
            }
            LebFormat::Leb2V2 => {
                let (chunk, chunk_idx) =
                    self.owning_chunk(bytes, body_id, entry, hint_jd, seg_u32, cache)?;
                let chunk = &chunk;
                let local = seg_idx
                    .checked_sub(chunk.segment_start() as usize)
                    .ok_or_else(|| Error::LebCorrupt {
                        what: format!("LEB2 chunk index (body {body_id})"),
                        offset: 0,
                    })?;
                let seg_count = chunk.segment_count();
                if let (Some(memo), LebSource::Slice(_), Ok(key_idx), Some(window)) = (
                    cache.filter(|memo| memo.segment_windows),
                    bytes,
                    u32::try_from(chunk_idx),
                    coefficient_window(entry, local, seg_count),
                ) {
                    let first = window.start;
                    let window_segments = window.end - first;
                    let words = memo.get_or_decode_key(
                        CacheKey::Window {
                            body: body_id,
                            chunk: key_idx,
                            first_segment: first,
                        },
                        || {
                            let what = format!("LEB2 chunk (body {body_id}, chunk {chunk_idx})");
                            let expected = reader::checked_chunk_uncompressed_size(
                                chunk,
                                entry.degree(),
                                entry.components(),
                                &what,
                            )?;
                            bytes.check_frame(expected)?;
                            let compressed =
                                bytes.blob(chunk.blob_offset(), chunk.compressed_size(), &what)?;
                            // Validate the complete frame before retaining only the
                            // selected rows. Slice images have no indexed leases.
                            let values = reader::decompress_body_window_words_limited(
                                &compressed,
                                expected,
                                seg_count,
                                entry.degree(),
                                entry.components(),
                                window,
                                bytes.window_limit(),
                            )?;
                            Ok(TrackedVec {
                                values,
                                _allocation: None,
                            })
                        },
                    )?;
                    return consume(slice_one_segment_in(
                        &words,
                        entry,
                        local - first as usize,
                        window_segments,
                        body_id,
                    )?);
                }
                let decode = || {
                    let what = format!("LEB2 chunk (body {body_id}, chunk {chunk_idx})");
                    let expected = reader::checked_chunk_uncompressed_size(
                        chunk,
                        entry.degree(),
                        entry.components(),
                        &what,
                    )?;
                    bytes.check_frame(expected)?;
                    let allocation = bytes.reserve_decode(expected, chunk.compressed_size())?;
                    let compressed =
                        bytes.blob(chunk.blob_offset(), chunk.compressed_size(), &what)?;
                    let values = reader::decompress_body_words_limited(
                        &compressed,
                        expected,
                        seg_count,
                        entry.degree(),
                        entry.components(),
                        bytes.window_limit(),
                    )?;
                    let capacity = values.capacity() * std::mem::size_of::<f64>();
                    drop(compressed);
                    let allocation = allocation.map(|lease| lease.retain(capacity)).transpose()?;
                    Ok(TrackedVec {
                        values,
                        _allocation: allocation,
                    })
                };
                match (cache, u32::try_from(chunk_idx)) {
                    (Some(memo), Ok(key_idx)) => memo.with_cached_body_segment(
                        (body_id, key_idx),
                        self,
                        entry,
                        local,
                        seg_count,
                        decode,
                        consume,
                    ),
                    _ => consume(slice_one_segment_in(
                        &decode()?,
                        entry,
                        local,
                        seg_count,
                        body_id,
                    )?),
                }
            }
        }
    }

    /// Locate the chunk owning global segment `seg_idx` at epoch `jd`:
    /// reuse the body's cached directory or decode it directly, hint with [`find_chunk`], confirm
    /// with [`select_chunk`]. An empty index is [`Error::LebCorrupt`]
    fn owning_chunk(
        &self,
        bytes: LebSource<'_>,
        body_id: i32,
        entry: &BodyEntry,
        jd: f64,
        seg_idx: u32,
        cache: Option<&LebChunkCache>,
    ) -> Result<(ChunkEntry, usize)> {
        // Uncached callers retain their original vector-only path; only the
        // cached path adds shared ownership. Both use identical selection logic.
        let decoded;
        let shared;
        let chunks = match cache {
            Some(memo) => {
                shared = memo.chunk_index(bytes, entry, body_id)?;
                shared.as_slice()
            }
            None => {
                decoded = decode::decode_chunk_index_source(bytes, entry.data_offset, body_id)?.1;
                decoded.as_slice()
            }
        };
        if chunks.is_empty() {
            return Err(Error::LebCorrupt {
                what: format!("LEB2 chunk index (body {body_id})"),
                offset: 0,
            });
        }
        let hint = find_chunk(chunks, jd).ok_or_else(|| Error::LebCorrupt {
            what: format!("LEB2 chunk index (body {body_id})"),
            offset: 0,
        })?;
        let owned = select_chunk(chunks, seg_idx, hint, body_id)?;
        Ok((chunks[owned].clone(), owned))
    }

    /// Evaluate stored nutation at Julian Day TT, returning `(dpsi,deps)` in radians.
    /// The selected polynomial uses Clenshaw evaluation. Missing/empty tables
    /// or zero segment intervals are corruption errors; non-finite epochs and
    /// stored-domain misses retain their distinct typed errors.
    pub(crate) fn eval_nutation_source(&self, bytes: LebSource<'_>, jd: f64) -> Result<(f64, f64)> {
        if !jd.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd });
        }
        let (base, header) = self.nutation_table(bytes)?;
        if jd < header.jd_start || jd > header.jd_end {
            return Err(Error::LebOutOfRange {
                body_id: -1,
                jd,
                start: header.jd_start,
                end: header.jd_end,
            });
        }
        let (seg_idx, tau_norm) = select_segment(
            header.jd_start,
            header.interval_days,
            header.segment_count,
            jd,
        )?;
        self.nutation_eval_at(bytes, base, &header, seg_idx, tau_norm.value())
    }

    /// Containing nutation polynomial's TT interval, in Julian days.
    /// Missing data or an uncovered epoch returns `None`; malformed metadata fails.
    /// Bounds follow the actual binary64 selector, including rounded joins.
    pub(crate) fn nutation_smooth_interval_source(
        &self,
        bytes: LebSource<'_>,
        jd: f64,
    ) -> Result<Option<(f64, f64)>> {
        if self.section(SECTION_NUTATION).is_none() {
            return Ok(None);
        }
        let (_, header) = self.nutation_table(bytes)?;
        if !(header.jd_start <= jd && jd <= header.jd_end) {
            return Ok(None);
        }
        super::segment_domain::smooth_segment_interval(
            header.jd_start,
            header.jd_end,
            header.interval_days,
            header.segment_count,
            jd,
        )
        .map(Some)
    }

    /// Exact TT derivative of the selected nutation polynomial, radians/day.
    /// Missing data or uncovered epochs return `None`; corrupt metadata fails.
    /// At a join the derivative uses the same binary64-selected record as its value.
    pub(crate) fn nutation_rate_source(
        &self,
        bytes: LebSource<'_>,
        jd: f64,
    ) -> Result<Option<(f64, f64)>> {
        if self.section(SECTION_NUTATION).is_none() {
            return Ok(None);
        }
        let (base, header) = self.nutation_table(bytes)?;
        if !(header.jd_start <= jd && jd <= header.jd_end) {
            return Ok(None);
        }
        let (segment, tau) = select_segment(
            header.jd_start,
            header.interval_days,
            header.segment_count,
            jd,
        )?;
        let data_base =
            base.checked_add(NUTATION_HEADER_SIZE)
                .ok_or_else(|| Error::LebCorrupt {
                    what: "LEB nutation table".to_owned(),
                    offset: base,
                })?;
        with_table_words(
            bytes,
            data_base,
            &header,
            segment,
            "LEB nutation table",
            |words| {
                let length = header.degree as usize + 1;
                if header.components != 2 || words.len() != 2 * length {
                    return Err(Error::LebCorrupt {
                        what: "LEB nutation table".to_owned(),
                        offset: base,
                    });
                }
                let scale = 2.0 / header.interval_days;
                // Preserve the conservative one-row working admission before
                // either checked derivative, including the bounded stack route.
                let _derivative_allocation =
                    bytes.reserve_work(length * std::mem::size_of::<f64>())?;
                let derivative = |coefficients: &[f64]| -> Result<f64> {
                    // Capacity is a storage choice, never a polynomial degree
                    // limit. Keep admission above and checked value evaluation
                    // below in their original order for both storage routes.
                    let mut inline = [0.0; 64];
                    let count = coefficients.len().saturating_sub(1).max(1);
                    let value = if count <= inline.len() {
                        let row = &mut inline[..count];
                        crate::math::chebyshev::chebyshev_deriv_coeffs_into(coefficients, row);
                        chebyshev_value(row, tau.value())?
                    } else {
                        chebyshev_value(&chebyshev_deriv_coeffs(coefficients), tau.value())?
                    };
                    Ok(value * scale)
                };
                Ok(Some((
                    derivative(&words[..length])?,
                    derivative(&words[length..])?,
                )))
            },
        )
    }

    /// Evaluate the nutation angles at the two-part epoch `epoch` (TT) —
    /// the split-epoch form of [`FileIndex::eval_nutation`].
    ///
    /// Segment ownership is defined by the collapsed word (identical to
    /// the scalar path on `epoch.value` by construction); `tau` keeps
    /// the sub-ulp tail via [`select_segment_two_part`]. Errors mirror
    /// the scalar path.
    pub(crate) fn eval_nutation_two_part_source(
        &self,
        bytes: LebSource<'_>,
        epoch: TwoPartEpoch,
    ) -> Result<(f64, f64)> {
        let jd = epoch.value();
        if !jd.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd });
        }
        let (base, header) = self.nutation_table(bytes)?;
        if jd < header.jd_start || jd > header.jd_end {
            return Err(Error::LebOutOfRange {
                body_id: -1,
                jd,
                start: header.jd_start,
                end: header.jd_end,
            });
        }
        let (seg_idx, tau_norm) = select_segment_two_part(
            header.jd_start,
            header.interval_days,
            header.segment_count,
            epoch,
        )?;
        self.nutation_eval_at(bytes, base, &header, seg_idx, tau_norm.value())
    }

    /// Decode and validate the nutation table header: missing section or
    /// absent table (`segment_count == 0`) and zero/non-finite interval
    /// are [`Error::LebCorrupt`], matching the release `eval_nutation`
    /// guard order. Returns the section base plus the header.
    fn nutation_table(&self, bytes: LebSource<'_>) -> Result<(usize, NutationHeader)> {
        let section = self
            .section(SECTION_NUTATION)
            .ok_or_else(|| Error::LebCorrupt {
                what: "LEB nutation table (absent)".to_string(),
                offset: 0,
            })?;
        let base = usize::try_from(section.offset()).map_err(|_| Error::LebCorrupt {
            what: "LEB nutation table".to_string(),
            offset: usize::MAX,
        })?;
        let header = bytes.nutation_header(base)?;
        if !header.is_present() {
            return Err(Error::LebCorrupt {
                what: "LEB nutation table (absent)".to_string(),
                offset: base,
            });
        }
        if !header.interval_days.is_finite() || header.interval_days <= 0.0 {
            return Err(Error::LebCorrupt {
                what: "LEB nutation table".to_string(),
                offset: base,
            });
        }
        Ok((base, header))
    }

    /// Evaluate one selected nutation segment at normalized epoch `tau`:
    /// coefficient read plus the value-only Clenshaw tail. Used by both
    /// the scalar ([`FileIndex::eval_nutation`]) and the two-part
    /// ([`FileIndex::eval_nutation_two_part`]) paths.
    fn nutation_eval_at(
        &self,
        bytes: LebSource<'_>,
        base: usize,
        header: &NutationHeader,
        seg_idx: usize,
        tau: f64,
    ) -> Result<(f64, f64)> {
        let data_base =
            base.checked_add(NUTATION_HEADER_SIZE)
                .ok_or_else(|| Error::LebCorrupt {
                    what: "LEB nutation table".to_string(),
                    offset: base,
                })?;
        with_table_words(
            bytes,
            data_base,
            header,
            seg_idx,
            "LEB nutation table",
            |words| eval_nutation_words(header, words, tau),
        )
    }

    /// Read stored Delta-T (TT−UT1) in days at Julian Day `jd`.
    /// Interior table values are linearly interpolated, with the last equal
    /// key selected. Outside values use the nearest endpoint. An empty/missing
    /// table is corrupt; non-finite epochs are refused. A zero-span pair uses
    /// its left value.
    pub(crate) fn delta_t_source(&self, bytes: LebSource<'_>, jd: f64) -> Result<f64> {
        if !jd.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd });
        }
        let section = self
            .section(SECTION_DELTA_T)
            .ok_or_else(|| Error::LebCorrupt {
                what: "LEB Delta-T table (absent)".to_string(),
                offset: 0,
            })?;
        let base = usize::try_from(section.offset()).map_err(|_| Error::LebCorrupt {
            what: "LEB Delta-T table".to_string(),
            offset: usize::MAX,
        })?;
        let header = bytes.delta_t_header(base).map_err(|error| match error {
            Error::LebCorrupt { .. } => Error::LebCorrupt {
                what: "LEB Delta-T table".into(),
                offset: base,
            },
            other => other,
        })?;
        header.validate_capacity(section.size(), base)?;
        let data_base = base
            .checked_add(DELTA_T_HEADER_SIZE)
            .ok_or_else(|| Error::LebCorrupt {
                what: "LEB Delta-T table".to_string(),
                offset: base,
            })?;
        if let LebSource::Slice(_) = bytes {
            let (jds, vals) = decode_delta_t_table(bytes, data_base, header.n_entries(), base)?;
            return eval_delta_t(&jds, &vals, jd);
        }
        evaluation::eval_delta_t_rows(header.n_entries() as usize, jd, |index| {
            let at = index
                .checked_mul(DELTA_T_ENTRY_SIZE)
                .and_then(|offset| data_base.checked_add(offset))
                .ok_or_else(|| Error::LebCorrupt {
                    what: "LEB Delta-T row".into(),
                    offset: data_base,
                })?;
            let row = bytes.range(at, DELTA_T_ENTRY_SIZE, "LEB Delta-T row")?;
            Ok((
                read_f64(&row, 0, "LEB Delta-T row", at)?,
                read_f64(&row, 8, "LEB Delta-T row", at)?,
            ))
        })
    }

    /// Look up a fixed-star record in this container image by internal id.
    ///
    /// The returned [`StarEntry`] contains the J2000 position, proper motion,
    /// parallax, radial velocity and magnitude in its documented units. All
    /// complete catalog records are decoded before return; repeated ids use
    /// the last complete record. A short trailing record is ignored according
    /// to [`star_count`]. A missing catalog or absent id is
    /// [`Error::LebStarNotFound`]; a malformed complete record is
    /// [`Error::LebCorrupt`]. No catalog storage is retained or mutated.
    pub(crate) fn get_star_source(&self, bytes: LebSource<'_>, star_id: i32) -> Result<StarEntry> {
        let section = self
            .section(SECTION_STARS)
            .ok_or(Error::LebStarNotFound { star_id })?;
        let base = usize::try_from(section.offset()).unwrap_or(usize::MAX);
        let count = star_count(section.size());
        let mut at = base;
        let mut found = None;
        for _ in 0..count {
            let record = bytes.range(at, STAR_ENTRY_SIZE, "LEB star catalog")?;
            let entry = StarEntry::decode(&record, 0).map_err(|_| Error::LebCorrupt {
                what: "LEB star catalog".to_string(),
                offset: at,
            })?;
            if entry.star_id() == star_id {
                found = Some(entry);
            }
            at = at
                .checked_add(STAR_ENTRY_SIZE)
                .ok_or(Error::LebStarNotFound { star_id })?;
        }
        found.ok_or(Error::LebStarNotFound { star_id })
    }
}

/// Borrow global segment `seg_idx` from a whole-body decoded word
/// buffer (LEB2 v1 monolith layout: every segment present).
pub(super) fn slice_one_segment<'a>(
    words: &'a [f64],
    entry: &BodyEntry,
    seg_idx: usize,
    body_id: i32,
) -> Result<&'a [f64]> {
    slice_one_segment_in(words, entry, seg_idx, entry.segment_count(), body_id)
}

/// Borrow segment `local_idx` from a decoded word buffer holding
/// `buffer_segments` segments (a v2 chunk's local table or a whole
/// monolith). Out-of-bounds slicing is [`Error::LebCorrupt`], never a
/// panic.
pub(super) fn slice_one_segment_in<'a>(
    words: &'a [f64],
    entry: &BodyEntry,
    local_idx: usize,
    buffer_segments: u32,
    body_id: i32,
) -> Result<&'a [f64]> {
    let range = segment_word_range(words.len(), entry, local_idx, buffer_segments, body_id)?;
    Ok(&words[range])
}

#[inline]
fn segment_word_range(
    word_count: usize,
    entry: &BodyEntry,
    local_idx: usize,
    buffer_segments: u32,
    body_id: i32,
) -> Result<std::ops::Range<usize>> {
    let invalid = || Error::LebCorrupt {
        what: format!("LEB2 coefficient blob (body {body_id})"),
        offset: 0,
    };
    let deg1 = (entry.degree() as usize)
        .checked_add(1)
        .ok_or_else(invalid)?;
    let per_seg = (entry.components() as usize)
        .checked_mul(deg1)
        .ok_or_else(invalid)?;
    let expect = (buffer_segments as usize)
        .checked_mul(per_seg)
        .ok_or_else(invalid)?;
    if word_count != expect || local_idx >= buffer_segments as usize {
        return Err(invalid());
    }
    let base = local_idx * per_seg;
    Ok(base..base + per_seg)
}

/// Conservative proof that the scalar and split kernels consume the same
/// coefficients and binary64 normalized argument. This checks representation,
/// not sampled coordinate agreement, and performs no coefficient acquisition.
fn scalar_and_split_share_argument(entry: &BodyEntry, jd: f64, epoch: TwoPartEpoch) -> bool {
    if entry.components != 3
        || !matches!(entry.coord_type, COORD_ICRS_BARY | COORD_ICRS_BARY_SYSTEM)
        || jd < entry.jd_start
        || jd > entry.jd_end
    {
        return false;
    }
    let Ok(canonical) = TwoPartEpoch::split(jd) else {
        return false;
    };
    if epoch.hi.to_bits() != canonical.hi.to_bits()
        || epoch.lo.to_bits() != canonical.lo.to_bits()
        || epoch.value().to_bits() != jd.to_bits()
    {
        return false;
    }
    let Ok((scalar_segment, scalar_tau)) =
        select_segment(entry.jd_start, entry.interval_days, entry.segment_count, jd)
    else {
        return false;
    };
    let Ok((split_segment, split_tau)) = select_segment_two_part(
        entry.jd_start,
        entry.interval_days,
        entry.segment_count,
        epoch,
    ) else {
        return false;
    };
    // These are the actual polynomial input keys, not evaluated coordinates.
    // Both routes use this same entry's interval scale and coefficient row.
    // Compare tau by bits so signed zero and one-word changes remain distinct.
    scalar_segment == split_segment && scalar_tau.value().to_bits() == split_tau.value().to_bits()
}

/// Select a segment and normalized epoch on a uniform Julian-Day grid.
///
/// `jd_start` is the first segment's start, `interval_days` is the positive
/// segment width in days, and `segment_count` is the nonzero number of segments.
/// `jd` is a finite scalar Julian Day in the same time scale as the grid (TT
/// for LEB body and nutation tables). This helper does not enforce coverage;
/// evaluation callers first check the stored body's or nutation table's range.
///
/// Rounded nominal starts `jd_start + segment * interval_days` define
/// half-open records: an interior join belongs to the following record.
/// The quotient `(jd - jd_start) / interval_days` supplies an index estimate;
/// a bounded integer search corrects it when rounded subtraction or division
/// disagrees with those starts. The first endpoint remains in record zero,
/// including when multiple starts round to that same epoch. Outside the grid,
/// selection clamps to its first or last record. Its midpoint is computed in order
/// as `seg_start = jd_start + segment * interval_days`, then
/// `seg_mid = seg_start + 0.5 * interval_days`. The dimensionless normalized
/// epoch is `tau = 2 * (jd - seg_mid) / interval_days`, clamped to `[-1, 1]`.
/// The final endpoint is clamped to the last segment. The result is a zero-based
/// segment index and a checked [`NormalizedTime`] for polynomial evaluation.
///
/// A non-finite requested epoch or derived segment quotient is
/// [`Error::NonFiniteArgument`], with quotient refusal before integer conversion.
/// A zero segment count or non-finite/non-positive interval is
/// [`Error::LebCorrupt`]. Finite quotients of any magnitude retain the index
/// clamp. Intermediate infinite `tau` is permitted and clamps to its endpoint;
/// a resulting NaN is refused by [`NormalizedTime::new`].
pub fn select_segment(
    jd_start: f64,
    interval_days: f64,
    segment_count: u32,
    jd: f64,
) -> Result<(usize, NormalizedTime)> {
    if !jd.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd });
    }
    if segment_count == 0 || !interval_days.is_finite() || interval_days <= 0.0 {
        return Err(Error::LebCorrupt {
            what: "LEB segment table (empty or non-positive interval)".to_string(),
            offset: 0,
        });
    }
    let quotient = (jd - jd_start) / interval_days;
    if !quotient.is_finite() {
        return Err(Error::NonFiniteArgument { value: quotient });
    }
    let quotient = quotient.trunc();
    let last = (segment_count - 1) as i64;
    let mut seg = (quotient as i64).clamp(0, last) as usize;
    if jd <= jd_start {
        seg = 0;
    } else {
        let before_start = jd < jd_start + seg as f64 * interval_days;
        let after_end = seg < last as usize && jd >= jd_start + (seg + 1) as f64 * interval_days;
        if before_start || after_end {
            // Multiplication and division are not inverse in binary64.
            // Search the monotone nominal starts rather than walking ulps or
            // assigning a query-dependent record. At most 32 index halvings
            // suffice for the u32 segment count, including duplicate starts.
            let mut lower = 0_u32;
            let mut upper = segment_count - 1;
            while lower < upper {
                let middle = lower + (upper - lower).div_ceil(2);
                if jd_start + f64::from(middle) * interval_days <= jd {
                    lower = middle;
                } else {
                    upper = middle - 1;
                }
            }
            seg = lower as usize;
        }
    }
    let seg_start = jd_start + seg as f64 * interval_days;
    let seg_mid = seg_start + 0.5 * interval_days;
    // Finite inputs can still overflow the ordered midpoint or tau arithmetic.
    // Infinite tau clamps to an endpoint; NormalizedTime rejects a NaN result.
    let tau = (2.0 * (jd - seg_mid) / interval_days).clamp(-1.0, 1.0);
    Ok((seg, NormalizedTime::new(tau)?))
}

/// Canonical segment selection plus normalized-epoch construction for a
/// two-part epoch, shared by the monolithic (LEB1), decoded-blob (LEB2)
/// and nutation evaluation paths.
///
/// Segment ownership is defined by the collapsed word: selection runs the
/// scalar [`select_segment`] on `epoch.value`, hence agrees with the
/// scalar path on that word by construction — including rounded interior
/// seams, the global end (last segment), and the
/// below-start/above-end edge clamps (callers range-check first, as in
/// the scalar paths). What the split path adds is precision: `tau`
/// comes from the ordered two-part displacement
/// `(hi - seg_mid) + lo` ([`TwoPartEpoch::diff`]), so a fractional tail
/// below half a Julian-day ulp (~2.4e-10 d near J2000) survives instead
/// of collapsing before the segment midpoint subtraction.
///
/// Documented boundary note (same class as the timescale sub-ulp note):
/// an epoch whose exact real value sits on the far side of a seam from
/// its collapsed word is owned by the collapsed word's segment. Callers
/// that step across a seam (e.g. retarded arguments) must re-split at
/// the stepped epoch — which [`retard_epoch`] does — rather than reuse
/// a stale split.
///
/// Errors: [`Error::NonFiniteArgument`] for a non-finite collapsed word,
/// derived segment quotient or ordered displacement (only reachable with
/// directly constructed hostile parts or a hostile grid far outside
/// any ephemeris domain); [`Error::LebCorrupt`] for a zero
/// `segment_count` or a non-finite/non-positive `interval_days`, as in
/// [`select_segment`].
pub fn select_segment_two_part(
    jd_start: f64,
    interval_days: f64,
    segment_count: u32,
    epoch: TwoPartEpoch,
) -> Result<(usize, NormalizedTime)> {
    let collapsed = epoch.value();
    let (seg, _) = select_segment(jd_start, interval_days, segment_count, collapsed)?;
    // Polynomial midpoint, recomputed exactly as in
    // `select_segment` so the displacement below anchors at the same
    // word the scalar path subtracts.
    let seg_start = jd_start + seg as f64 * interval_days;
    let seg_mid = seg_start + 0.5 * interval_days;
    // `select_segment` already established a finite positive interval;
    // a non-finite midpoint can result from a non-finite grid input or
    // overflow in the ordered sum of finite grid words. Both are surfaced
    // as a typed error rather than a panic.
    let mid = TwoPartEpoch::split(seg_mid)?;
    let dt = epoch.diff(mid);
    if !dt.is_finite() {
        return Err(Error::NonFiniteArgument { value: dt });
    }
    // Finite `dt` and finite positive `interval_days`: `clamp` cannot
    // hit its NaN-panic arm, matching the release `tau` clamp.
    let tau = (2.0 * dt / interval_days).clamp(-1.0, 1.0);
    NormalizedTime::new(tau).map(|t| (seg, t))
}

impl ChunkIndexHeader {
    /// Decode the header at `offset` (pure structural decode; the two
    /// reserved words are skipped, matching the release unpacker).
    pub fn decode(bytes: &[u8], offset: usize) -> Result<Self> {
        const WHAT: &str = "LEB2 chunk index header";
        let chunk = slice_at(bytes, offset, CHUNK_INDEX_HEADER_SIZE, WHAT)?;
        Ok(Self {
            chunk_count: read_u32(chunk, 0, WHAT, offset)?,
            chunk_interval_days: read_f32(chunk, 8, WHAT, offset)?,
        })
    }

    /// Declared chunk entry count.
    pub fn chunk_count(&self) -> u32 {
        self.chunk_count
    }
    /// Nominal chunk interval in days (`f32`, advisory only).
    pub fn chunk_interval_days(&self) -> f32 {
        self.chunk_interval_days
    }
}

/// Confirm that global segment `seg_idx` is owned by the hinted chunk,
/// probing the immediate neighbours before reporting corruption.
///
/// In the half-ulp window where a
/// rounded chunk JD bound disagrees with the segment arithmetic, the
/// segment index wins: a hint that misses is retried at `hint - 1` and
/// `hint + 1`, and only a miss in all three is [`Error::LebCorrupt`]
/// (never a silent neighbouring-segment evaluation).
/// Guard failures the release raises as `ValueError` surface here as [`Error::LebCorrupt`] (typed form, same failure).
pub fn select_chunk(
    chunks: &[ChunkEntry],
    seg_idx: u32,
    hint: usize,
    body_id: i32,
) -> Result<usize> {
    if hint < chunks.len() && chunks[hint].covers(seg_idx) {
        return Ok(hint);
    }
    for a in [hint.checked_sub(1), hint.checked_add(1)]
        .into_iter()
        .flatten()
    {
        if a < chunks.len() && chunks[a].covers(seg_idx) {
            return Ok(a);
        }
    }
    Err(Error::LebCorrupt {
        what: format!("LEB2 chunk index (body {body_id})"),
        offset: 0,
    })
}

/// Reorder coefficient bytes from coefficient-major back to
/// segment-major layout (the reader direction of release
/// `leb_compression.reorder_segment_major`, applied after unshuffling).
/// Exact bit-pattern permutation; same length contract as
/// [`reorder_to_coeff_major`].
pub fn reorder_to_segment_major(
    data: &[u8],
    segment_count: u32,
    degree: u32,
    components: u32,
) -> Result<Vec<u8>> {
    const WHAT: &str = "LEB2 coefficient layout";
    reorder_words(data, segment_count, degree, components, WHAT, false)
}

impl LebTier {
    /// Selection priority, highest first.
    pub const PRIORITY: [LebTier; 3] = [LebTier::Base, LebTier::Medium, LebTier::Extended];
    /// Canonical tier token as encoded in `{tier}_{group}` file names.
    pub fn name(self) -> &'static str {
        match self {
            LebTier::Base => "base",
            LebTier::Medium => "medium",
            LebTier::Extended => "extended",
        }
    }
    /// Parse a tier token; anything else is not a tier (`None`), never an
    /// error — files without a recognizable tier are unconstrained by the
    /// tier guard.
    pub fn parse(token: &str) -> Option<LebTier> {
        match token {
            "base" => Some(LebTier::Base),
            "medium" => Some(LebTier::Medium),
            "extended" => Some(LebTier::Extended),
            _ => None,
        }
    }
}

/// Tier token encoded in a container file name, or `None` when no
/// recognizable tier is present. Mirrors release `_file_tier`: a
/// `{tier}_{group}` prefix that is itself a known tier wins; otherwise
/// (custom multi-token prefixes such as `myset_base`) the tokens are
/// scanned and the first known tier token is reported.
pub fn file_tier(file_name: &str) -> Option<LebTier> {
    let stem = if file_name.contains('.') {
        file_name
            .rsplit_once('.')
            .map(|(s, _)| s)
            .unwrap_or(file_name)
    } else {
        file_name
    };
    let parts: Vec<&str> = stem.split('_').collect();
    if parts.len() >= 2 && file_group_suffix(file_name).is_some() {
        let prefix = parts[..parts.len() - 1].join("_");
        if let Some(tier) = LebTier::parse(&prefix) {
            return Some(tier);
        }
    }
    parts.iter().find_map(|token| LebTier::parse(token))
}

impl FileIndex {
    /// Presence probe for usable nutation data: the nutation section is
    /// present, its header decodes, the table is marked present
    /// (`segment_count > 0`, the release `has_nutation` rule — a
    /// `--skip-aux` header with zero segments must not shadow a later
    /// file), and the interval is finite positive. Never errors: any
    /// decode problem reports absence, and evaluation owns validation.
    pub fn has_nutation_data(&self, bytes: &[u8]) -> bool {
        self.has_nutation_data_source(LebSource::Slice(bytes))
            .unwrap_or(false)
    }

    /// Fallible acquisition preserves positional I/O/resource errors in composites.
    pub(super) fn has_nutation_data_source(&self, bytes: LebSource<'_>) -> Result<bool> {
        let Some(section) = self.section(SECTION_NUTATION) else {
            return Ok(false);
        };
        let Ok(base) = usize::try_from(section.offset()) else {
            return Ok(false);
        };
        let header = match bytes.nutation_header(base) {
            Ok(header) => header,
            Err(Error::LebCorrupt { .. }) => return Ok(false),
            Err(error) => return Err(error),
        };
        Ok(header.is_present() && header.interval_days.is_finite() && header.interval_days > 0.0)
    }

    /// Presence probe for usable Delta-T data: the section is present,
    /// its header decodes, and at least one pair is stored. Never
    /// errors; evaluation owns validation.
    /// The release selects on the decoded non-empty `_delta_t_jds` table; probing `n_entries > 0` here is the typed equivalent.
    pub fn has_delta_t_data(&self, bytes: &[u8]) -> bool {
        self.has_delta_t_data_source(LebSource::Slice(bytes))
            .unwrap_or(false)
    }

    /// Presence is tolerant of malformed header fields, not failed file reads.
    fn has_delta_t_data_source(&self, bytes: LebSource<'_>) -> Result<bool> {
        let Some(section) = self.section(SECTION_DELTA_T) else {
            return Ok(false);
        };
        let Ok(base) = usize::try_from(section.offset()) else {
            return Ok(false);
        };
        let header = match bytes.delta_t_header(base) {
            Ok(header) => header,
            Err(Error::LebCorrupt { .. }) => return Ok(false),
            Err(error) => return Err(error),
        };
        Ok(header.n_entries() > 0)
    }

    /// Presence probe for a usable star catalog: the section is present
    /// and holds at least one complete 64-byte entry. Never errors; evaluation owns
    /// validation.
    pub fn has_star_data(&self, bytes: &[u8]) -> bool {
        let _ = bytes;
        self.section(SECTION_STARS)
            .is_some_and(|section| star_count(section.size()) > 0)
    }
}

impl CompositeSources<'_> {
    fn len(&self) -> usize {
        match self {
            Self::Raw(entries) => entries.len(),
            Self::Shared(files) => files.len(),
        }
    }

    pub(super) fn get(&self, index: usize) -> (&FileIndex, LebSource<'_>) {
        match self {
            Self::Raw(entries) => (&entries[index].0, entries[index].1),
            Self::Shared(files) => {
                // Construction checks every handle. These private immutable
                // clones cannot be closed by another owner of the same image.
                let image = files[index].image.as_ref().expect("validated LEB owner");
                (&image.index, image.bytes.source())
            }
        }
    }

    fn iter(&self) -> impl Iterator<Item = (&FileIndex, LebSource<'_>)> {
        (0..self.len()).map(|index| self.get(index))
    }

    fn files(&self) -> &[LebFile] {
        match self {
            Self::Raw(_) => &[],
            Self::Shared(files) => files,
        }
    }
}

impl<'a> CompositeIndex<'a> {
    /// Build the first-wins dispatch over at least one decoded reader.
    /// Each index must describe its paired immutable byte slice. Slices remain
    /// borrowed for `'a`; private image-specific caches are created here.
    /// An empty list is [`Error::LebCorrupt`].
    pub fn new(entries: Vec<(FileIndex, &'a [u8])>) -> Result<Self> {
        let caches = (0..entries.len())
            .map(|_| Arc::new(LebChunkCache::new()))
            .collect();
        let entries = entries
            .into_iter()
            .map(|(index, bytes)| (index, LebSource::Slice(bytes)))
            .collect();
        Self::from_parts(CompositeSources::Raw(entries), caches)
    }

    fn from_parts(entries: CompositeSources<'a>, caches: Vec<Arc<LebChunkCache>>) -> Result<Self> {
        if entries.len() == 0 {
            return Err(Error::LebCorrupt {
                what: "LEB composite (no readers)".to_string(),
                offset: 0,
            });
        }
        let mut body_map: HashMap<i32, BodyLocation> = HashMap::new();
        for (source_position, (index, _)) in entries.iter().enumerate() {
            for (body_row, body) in index.bodies().iter().enumerate() {
                body_map.entry(body.body_id()).or_insert(BodyLocation {
                    source_position,
                    body_row,
                });
            }
        }
        let builtin_bodies = Box::new(std::array::from_fn(|id| {
            body_map.get(&(id as i32)).copied()
        }));
        let mut nutation_at = None;
        let mut delta_t_at = None;
        let mut stars_at = None;
        for (i, (index, bytes)) in entries.iter().enumerate() {
            if nutation_at.is_none() && index.has_nutation_data_source(bytes)? {
                nutation_at = Some(i);
            }
            if delta_t_at.is_none() && index.has_delta_t_data_source(bytes)? {
                delta_t_at = Some(i);
            }
            if stars_at.is_none()
                && index
                    .section(SECTION_STARS)
                    .is_some_and(|section| star_count(section.size()) > 0)
            {
                stars_at = Some(i);
            }
        }
        Ok(Self {
            entries,
            caches,
            body_map,
            builtin_bodies,
            nutation_at,
            delta_t_at,
            stars_at,
        })
    }

    /// Bind immutable file handles, sharing their decoded caches across
    /// provider rebuilds and context clones. Closed handles fail; empty input
    /// has the same error as [`Self::new`]. Neither image bytes nor decoded
    /// indexes are copied; raw-index providers retain their separate owned API.
    /// The private image/cache pairing prevents cross-file chunk-key aliasing.
    pub fn from_files(files: impl IntoIterator<Item = &'a LebFile>) -> Result<Self> {
        Self::from_shared_files(files)
    }

    /// Retain immutable image owners without a borrow from the supplied handles.
    /// Images, indexes and bounded caches are shared; no file data are copied.
    /// Closed handles and malformed metadata retain ordinary reader errors.
    pub(crate) fn from_shared_files<'f>(
        files: impl IntoIterator<Item = &'f LebFile>,
    ) -> Result<CompositeIndex<'static>> {
        let files: Vec<_> = files.into_iter().cloned().collect();
        let caches = files
            .iter()
            .map(|file| file.chunk_cache().cloned())
            .collect::<Result<_>>()?;
        CompositeIndex::from_parts(CompositeSources::Shared(files), caches)
    }

    /// Number of constituent readers.
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    /// Always false after successful construction; kept so the
    /// `len`/`is_empty` pair stays complete.
    pub fn is_empty(&self) -> bool {
        self.entries.len() == 0
    }

    /// Body ids served by this composite, in unordered `HashMap` order (first-wins applies to serving dispatch, not to this listing).
    pub fn body_ids(&self) -> Vec<i32> {
        self.body_map.keys().copied().collect()
    }

    /// True when any constituent indexes `body_id`.
    pub fn has_body(&self, body_id: i32) -> bool {
        self.body_location(body_id).is_some()
    }

    /// Constituent position serving `body_id` (first-wins), for
    /// inventory consumers.
    pub fn serving_index(&self, body_id: i32) -> Option<usize> {
        self.body_location(body_id)
            .map(|location| location.source_position)
    }

    /// Look up immutable dispatch without hashing the fixed built-in vocabulary.
    /// Missing built-ins remain absent; all other IDs retain randomized lookup.
    fn body_location(&self, body_id: i32) -> Option<BodyLocation> {
        match usize::try_from(body_id)
            .ok()
            .filter(|id| *id < BUILTIN_BODY_COUNT)
        {
            Some(id) => self.builtin_bodies[id],
            None => self.body_map.get(&body_id).copied(),
        }
    }

    /// Decoded index entry actually serving `body_id` — the entry whose
    /// `coord_type` classifies the channel (planet center versus system
    /// barycenter) for the calculation pipeline.
    pub fn serving_entry(&self, body_id: i32) -> Option<&BodyEntry> {
        let location = self.body_location(body_id)?;
        self.entries
            .get(location.source_position)
            .0
            .bodies()
            .get(location.body_row)
    }

    /// Stored per-body JD coverage from the serving reader, or `None`
    /// when no constituent carries the body.
    pub fn body_coverage(&self, body_id: i32) -> Option<(f64, f64)> {
        self.serving_entry(body_id)
            .map(|entry| (entry.jd_start(), entry.jd_end()))
    }

    /// TT interval of the body's actual serving polynomial.
    /// Missing bodies/coverage return `None`; corrupt segment metadata fails.
    pub fn body_smooth_interval(&self, body_id: i32, jd: f64) -> Result<Option<(f64, f64)>> {
        let Some(index) = self.serving_index(body_id) else {
            return Ok(None);
        };
        self.entries.get(index).0.body_smooth_interval(body_id, jd)
    }

    /// TT interval of the selected nutation polynomial, or `None` if absent.
    /// Malformed selected metadata returns its reader error.
    pub fn nutation_smooth_interval(&self, jd: f64) -> Result<Option<(f64, f64)>> {
        let Some(index) = self.nutation_at else {
            return Ok(None);
        };
        let (index, bytes) = self.entries.get(index);
        index.nutation_smooth_interval_source(bytes, jd)
    }

    /// Aggregate file-header Julian-day bounds in constituent order.
    ///
    /// Construction guarantees at least one reader. Its stored bounds seed
    /// the result; later starts replace only under `<`, and later ends only
    /// under `>`. Equal bounds retain the first word, including signed zero.
    /// A first NaN remains selected, while later NaNs cannot replace a selected
    /// value. Header bounds are carried; per-body coverage separately
    /// governs position evaluation.
    pub fn jd_range(&self) -> (f64, f64) {
        let first = self.entries.get(0).0.header();
        let mut start = first.jd_start();
        let mut end = first.jd_end();
        for (index, _) in self.entries.iter().skip(1) {
            let candidate_start = index.header().jd_start();
            let candidate_end = index.header().jd_end();
            if candidate_start < start {
                start = candidate_start;
            }
            if candidate_end > end {
                end = candidate_end;
            }
        }
        (start, end)
    }

    /// True when any constituent carries usable nutation data.
    pub fn has_nutation(&self) -> bool {
        self.nutation_at.is_some()
    }

    /// Prove identity on the same attach-order image as ordinary evaluation.
    pub(crate) fn scalar_state_matches_split(
        &self,
        body: i32,
        jd: f64,
        epoch: TwoPartEpoch,
    ) -> bool {
        self.serving_index(body).is_some_and(|image| {
            self.entries
                .get(image)
                .0
                .scalar_state_matches_split(body, jd, epoch)
        })
    }

    /// Inline-storage projection of the ordinary scalar stored-state route.
    pub(crate) fn eval_state_inline(&self, body: i32, jd: f64) -> Result<Option<CartesianState>> {
        let Some(image) = self.serving_index(body) else {
            return Err(Error::LebBodyNotFound { body_id: body });
        };
        self.evaluation(image).eval_inline(body, jd)
    }

    /// Inline-storage projection retaining ordinary split-epoch selection.
    pub(crate) fn eval_state_inline_two_part(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
    ) -> Result<Option<CartesianState>> {
        let Some(image) = self.serving_index(body) else {
            return Err(Error::LebBodyNotFound { body_id: body });
        };
        self.evaluation(image).eval_inline_two_part(body, epoch)
    }

    /// Evaluate a body's position and velocity at TT Julian Day `jd`.
    ///
    /// Dispatches to the selected reader's monolithic LEB1 coefficients or
    /// bounded decoded LEB2 cache. Values retain their stored channel units;
    /// rates are per TT day. Missing bodies return [`Error::LebBodyNotFound`];
    /// range and corruption errors propagate from the selected reader.
    pub fn eval_body_any(&self, body_id: i32, jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        let i = self
            .serving_index(body_id)
            .ok_or(Error::LebBodyNotFound { body_id })?;
        self.evaluation(i).eval_any(body_id, jd)
    }

    /// Two-part-epoch form of [`CompositeIndex::eval_body_any`]: segment
    /// ownership keys on the collapsed word, `tau` keeps the sub-ulp
    /// tail (the split-epoch layer). The LEB2 memo path mirrors the
    /// scalar form through the serving constituent's handle.
    pub fn eval_body_any_two_part(
        &self,
        body_id: i32,
        epoch: TwoPartEpoch,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        let i = self
            .serving_index(body_id)
            .ok_or(Error::LebBodyNotFound { body_id })?;
        self.evaluation(i).eval_any_two_part(body_id, epoch)
    }

    /// `(hits, misses, resident_entries)` summed over the constituent
    /// memos, for measurement. LEB1-only composites always report zeros.
    pub fn chunk_cache_stats(&self) -> (u64, u64, usize) {
        self.caches.iter().fold((0, 0, 0), |(h, m, e), c| {
            let (ch, cm, ce) = c.stats();
            (h + ch, m + cm, e + ce)
        })
    }

    /// Multi-epoch form of [`CompositeIndex::eval_body_any`]: one
    /// scalar dispatch per element, in input order, failing fast with
    /// the first element's scalar error (the multi-epoch contract; same loop contract as
    /// [`FileIndex::eval_body_many`]).
    pub fn eval_body_any_many(
        &self,
        body_id: i32,
        jds: &[f64],
    ) -> Result<Vec<(Vec<f64>, Vec<f64>)>> {
        batch_eval(jds, |jd| self.eval_body_any(body_id, jd))
    }

    /// Multi-epoch form of [`CompositeIndex::eval_body_any_two_part`].
    pub fn eval_body_any_two_part_many(
        &self,
        body_id: i32,
        epochs: &[TwoPartEpoch],
    ) -> Result<Vec<(Vec<f64>, Vec<f64>)>> {
        batch_eval_epochs(epochs, |epoch| self.eval_body_any_two_part(body_id, epoch))
    }

    /// Nutation pair from the first constituent carrying usable data.
    /// No carrier is [`Error::LebCorrupt`].
    pub fn eval_nutation(&self, jd: f64) -> Result<(f64, f64)> {
        let i = self.nutation_at.ok_or_else(|| Error::LebCorrupt {
            what: "LEB composite nutation (absent)".to_string(),
            offset: 0,
        })?;
        self.evaluation(i).nutation(jd)
    }

    /// Exact selected nutation-polynomial derivative in radians per TT day.
    /// Missing data returns `None`; malformed selected data returns its reader error.
    pub fn nutation_rate(&self, jd: f64) -> Result<Option<(f64, f64)>> {
        let Some(position) = self.nutation_at else {
            return Ok(None);
        };
        self.evaluation(position).nutation_rate(jd)
    }

    /// Delta-T from the first constituent carrying a usable table. No
    /// carrier is [`Error::LebCorrupt`].
    pub fn delta_t(&self, jd: f64) -> Result<f64> {
        let i = self.delta_t_at.ok_or_else(|| Error::LebCorrupt {
            what: "LEB composite Delta-T (absent)".to_string(),
            offset: 0,
        })?;
        self.evaluation(i).delta_t(jd)
    }

    /// Return a selected Delta-T value in days, or `None` for table absence.
    ///
    /// A selected immutable carrier evaluates at Julian Day `jd` through
    /// [`Self::delta_t`], preserving every value and error. Without a carrier,
    /// each present table header must decode and its declared count must fit
    /// the section before absence is returned. This distinguishes genuine
    /// missing/empty tables from corrupt headers hidden by the non-throwing
    /// construction probe. It neither selects a different constituent nor
    /// validates, sorts or changes stored epoch/value words.
    ///
    /// Malformed metadata is [`Error::LebCorrupt`]; a selected carrier also
    /// retains its strict [`Error::NonFiniteArgument`] request guard. No
    /// evaluation error is converted to `None`.
    fn delta_t_for_tier(&self, jd: f64) -> Result<Option<f64>> {
        if self.delta_t_at.is_some() {
            return self.delta_t(jd).map(Some);
        }
        for (index, bytes) in self.entries.iter() {
            if let Some(section) = index.section(SECTION_DELTA_T) {
                let base = usize::try_from(section.offset()).map_err(|_| Error::LebCorrupt {
                    what: "LEB Delta-T table".to_string(),
                    offset: usize::MAX,
                })?;
                let header = bytes.delta_t_header(base)?;
                header.validate_capacity(section.size(), base)?;
            }
        }
        Ok(None)
    }

    /// Star entry from the first constituent carrying a usable catalog.
    /// No carrier is [`Error::LebStarNotFound`].
    pub fn get_star(&self, star_id: i32) -> Result<StarEntry> {
        let i = self.stars_at.ok_or(Error::LebStarNotFound { star_id })?;
        let (index, bytes) = self.entries.get(i);
        index.get_star_source(bytes, star_id)
    }
}

impl<'a> TieredIndex<'a> {
    /// Build the tiered selector. Tiers are reordered into priority
    /// order; an empty list and a duplicated tier are
    /// [`Error::LebCorrupt`].
    /// A duplicated tier surfaces here as [`Error::LebCorrupt`]; it has no release counterpart (mapping keys are unique), while the release's unknown-tier `ValueError` is unrepresentable here (`LebTier` is typed).
    pub fn new(tiers: Vec<(LebTier, CompositeIndex<'a>)>) -> Result<Self> {
        if tiers.is_empty() {
            return Err(Error::LebCorrupt {
                what: "LEB tiered index (no tiers)".to_string(),
                offset: 0,
            });
        }
        let mut ordered: Vec<(LebTier, CompositeIndex<'a>)> = Vec::with_capacity(tiers.len());
        for (tier, composite) in tiers {
            if ordered.iter().any(|(installed, _)| *installed == tier) {
                return Err(Error::LebCorrupt {
                    what: format!("LEB tiered index (duplicate tier {})", tier.name()),
                    offset: 0,
                });
            }
            ordered.push((tier, composite));
        }
        // Priority order governs selection, never installation order.
        ordered.sort_by_key(|(tier, _)| {
            LebTier::PRIORITY
                .iter()
                .position(|t| *t == *tier)
                .unwrap_or(LebTier::PRIORITY.len())
        });
        let mut candidates: HashMap<i32, Vec<TierCandidate>> = HashMap::new();
        for (pos, (_, composite)) in ordered.iter().enumerate() {
            for body_id in composite.body_ids() {
                candidates.entry(body_id).or_default().push(TierCandidate {
                    position: pos,
                    bounds: composite.body_coverage(body_id),
                });
            }
        }
        Ok(Self {
            tiers: ordered,
            candidates,
        })
    }

    /// Tier of the constituent at priority position `pos`.
    pub fn tier_at(&self, pos: usize) -> Option<LebTier> {
        self.tiers.get(pos).map(|(tier, _)| *tier)
    }

    /// True when any installed tier serves `body_id`.
    pub fn has_body(&self, body_id: i32) -> bool {
        self.candidates.contains_key(&body_id)
    }

    /// Highest-priority tier position covering `body_id` at `jd`
    /// (inclusive stored bounds, like the release `selected_body_reader`).
    /// A non-finite `jd` never satisfies a bound, so it selects nothing.
    pub fn selected(&self, body_id: i32, jd: f64) -> Option<usize> {
        self.candidates.get(&body_id)?.iter().find_map(|candidate| {
            candidate
                .bounds
                .filter(|(start, end)| *start <= jd && jd <= *end)
                .map(|_| candidate.position)
        })
    }

    /// Tier selected for `body_id` at `jd`, if any tier covers it.
    pub fn selected_tier(&self, body_id: i32, jd: f64) -> Option<LebTier> {
        self.selected(body_id, jd).and_then(|pos| self.tier_at(pos))
    }

    /// Serving tier position for `body_id`, optionally selected for a
    /// target date (`None` returns the highest-priority candidate —
    /// the release `body_reader(body_id, jd=None)` rule). The entry is
    /// used for coordinate-channel metadata only; numerical evaluation
    /// still goes through [`TieredIndex::eval_body_any`] date selection.
    pub fn serving(&self, body_id: i32, jd: Option<f64>) -> Option<usize> {
        match jd {
            Some(jd) => self.selected(body_id, jd),
            None => self.candidates.get(&body_id)?.first().map(|c| c.position),
        }
    }

    /// Decoded index entry serving `body_id` at the highest-priority
    /// candidate (channel metadata only; see [`TieredIndex::serving`]).
    pub fn serving_entry(&self, body_id: i32) -> Option<&BodyEntry> {
        let pos = self.serving(body_id, None)?;
        self.tiers[pos].1.serving_entry(body_id)
    }

    /// Union of the stored per-body intervals across all installed tiers,
    /// or `None` when no tier serves the body.
    pub fn body_coverage(&self, body_id: i32) -> Option<(f64, f64)> {
        let positions = self.candidates.get(&body_id)?;
        let mut start = f64::INFINITY;
        let mut end = f64::NEG_INFINITY;
        for candidate in positions {
            if let Some((s, e)) = candidate.bounds {
                start = start.min(s);
                end = end.max(e);
            }
        }
        if start.is_finite() && end.is_finite() {
            Some((start, end))
        } else {
            None
        }
    }

    /// TT interval of the polynomial in the actual date-selected tier.
    /// Missing bodies/coverage return `None`; invalid metadata returns an error.
    pub fn body_smooth_interval(&self, body_id: i32, jd: f64) -> Result<Option<(f64, f64)>> {
        let Some(index) = self.selected(body_id, jd) else {
            return Ok(None);
        };
        let Some((mut start, mut end)) = self.tiers[index].1.body_smooth_interval(body_id, jd)?
        else {
            return Ok(None);
        };
        if let Some(candidates) = self.candidates.get(&body_id) {
            for candidate in candidates
                .iter()
                .take_while(|candidate| candidate.position != index)
            {
                if let Some((lower, upper)) = candidate.bounds {
                    if upper < jd {
                        start = start.max(upper.next_up());
                    }
                    if lower > jd {
                        end = end.min(lower.next_down());
                    }
                }
            }
        }
        Ok(Some((start, end)))
    }

    /// TT interval of the date-selected nutation polynomial, or `None`.
    /// The tier and carrier rules are identical to nutation evaluation.
    pub fn nutation_smooth_interval(&self, jd: f64) -> Result<Option<(f64, f64)>> {
        for (_, composite) in &self.tiers {
            let (start, end) = composite.jd_range();
            if start <= jd && jd <= end && composite.has_nutation() {
                return composite.nutation_smooth_interval(jd);
            }
        }
        Ok(None)
    }

    /// Aggregate file-header Julian-day bounds in installed tier priority.
    ///
    /// Construction guarantees at least one tier and orders it before lookup.
    /// The first composite range seeds strict `<` / `>` reductions: equal
    /// bounds retain its words, first NaNs persist, and later NaNs are ignored.
    /// Per-body serving intervals remain separate from these header bounds.
    pub fn jd_range(&self) -> (f64, f64) {
        let (mut start, mut end) = self.tiers[0].1.jd_range();
        for (_, composite) in self.tiers.iter().skip(1) {
            let (candidate_start, candidate_end) = composite.jd_range();
            if candidate_start < start {
                start = candidate_start;
            }
            if candidate_end > end {
                end = candidate_end;
            }
        }
        (start, end)
    }

    /// True when any installed tier carries usable nutation data.
    pub fn has_nutation(&self) -> bool {
        self.tiers.iter().any(|(_, c)| c.has_nutation())
    }

    /// Prove identity on the actual date-selected tier, never its metadata-only
    /// highest-priority entry. Missing coverage retains the ordinary error path.
    pub(crate) fn scalar_state_matches_split(
        &self,
        body: i32,
        jd: f64,
        epoch: TwoPartEpoch,
    ) -> bool {
        self.selected(body, jd).is_some_and(|selected| {
            self.tiers[selected]
                .1
                .scalar_state_matches_split(body, jd, epoch)
        })
    }

    /// Inline-storage projection of the ordinary scalar stored-state route.
    pub(crate) fn eval_state_inline(&self, body: i32, jd: f64) -> Result<Option<CartesianState>> {
        let candidates = self
            .candidates
            .get(&body)
            .ok_or(Error::LebBodyNotFound { body_id: body })?;
        let selected = self.selected(body, jd).unwrap_or(candidates[0].position);
        self.tiers[selected].1.eval_state_inline(body, jd)
    }

    /// Inline-storage projection retaining ordinary split-epoch selection.
    pub(crate) fn eval_state_inline_two_part(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
    ) -> Result<Option<CartesianState>> {
        let candidates = self
            .candidates
            .get(&body)
            .ok_or(Error::LebBodyNotFound { body_id: body })?;
        let selected = self
            .selected(body, epoch.value())
            .unwrap_or(candidates[0].position);
        self.tiers[selected]
            .1
            .eval_state_inline_two_part(body, epoch)
    }

    /// Evaluate a body's position and velocity at TT Julian Day `jd`.
    ///
    /// Uses the highest-priority tier covering the epoch. Values retain their
    /// stored channel units and rates are per TT day. An unserved body returns
    /// [`Error::LebBodyNotFound`]. If no tier covers the epoch, the ordinary
    /// evaluation on the highest-priority candidate returns its range error.
    pub fn eval_body_any(&self, body_id: i32, jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        let positions = self
            .candidates
            .get(&body_id)
            .ok_or(Error::LebBodyNotFound { body_id })?;
        let pos = self.selected(body_id, jd).unwrap_or(positions[0].position);
        self.tiers[pos].1.eval_body_any(body_id, jd)
    }

    /// Two-part-epoch form of [`TieredIndex::eval_body_any`].
    pub fn eval_body_any_two_part(
        &self,
        body_id: i32,
        epoch: TwoPartEpoch,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        let positions = self
            .candidates
            .get(&body_id)
            .ok_or(Error::LebBodyNotFound { body_id })?;
        let pos = self
            .selected(body_id, epoch.value())
            .unwrap_or(positions[0].position);
        self.tiers[pos].1.eval_body_any_two_part(body_id, epoch)
    }

    /// Multi-epoch form of [`TieredIndex::eval_body_any`]: per-epoch
    /// tier selection still applies element-wise (each epoch is served
    /// by its own highest-priority covering tier), failing fast with
    /// the first element's scalar error (the multi-epoch contract).
    pub fn eval_body_any_many(
        &self,
        body_id: i32,
        jds: &[f64],
    ) -> Result<Vec<(Vec<f64>, Vec<f64>)>> {
        batch_eval(jds, |jd| self.eval_body_any(body_id, jd))
    }

    /// Multi-epoch form of [`TieredIndex::eval_body_any_two_part`].
    pub fn eval_body_any_two_part_many(
        &self,
        body_id: i32,
        epochs: &[TwoPartEpoch],
    ) -> Result<Vec<(Vec<f64>, Vec<f64>)>> {
        batch_eval_epochs(epochs, |epoch| self.eval_body_any_two_part(body_id, epoch))
    }

    /// `(hits, misses, resident_entries)` summed over the installed
    /// tier composites' constituent memos, for measurement.
    pub fn chunk_cache_stats(&self) -> (u64, u64, usize) {
        self.tiers
            .iter()
            .map(|(_, c)| c.chunk_cache_stats())
            .fold((0, 0, 0), |(h, m, e), (ch, cm, ce)| {
                (h + ch, m + cm, e + ce)
            })
    }

    /// Nutation pair from the highest-priority tier whose file range
    /// covers `jd` and carries usable data. No covering carrier is
    /// [`Error::LebCorrupt`].
    pub fn eval_nutation(&self, jd: f64) -> Result<(f64, f64)> {
        for (_, composite) in &self.tiers {
            let (start, end) = composite.jd_range();
            if start <= jd && jd <= end && composite.has_nutation() {
                return composite.eval_nutation(jd);
            }
        }
        Err(Error::LebCorrupt {
            what: format!("LEB tiered nutation (no coverage at JD {jd})"),
            offset: 0,
        })
    }

    /// Exact nutation derivative from the date-selected tier, radians/TT day.
    /// Missing coverage/data returns `None`; selected corruption is an error.
    pub fn nutation_rate(&self, jd: f64) -> Result<Option<(f64, f64)>> {
        for (_, composite) in &self.tiers {
            let (start, end) = composite.jd_range();
            if start <= jd && jd <= end && composite.has_nutation() {
                return composite.nutation_rate(jd);
            }
        }
        Ok(None)
    }

    /// Look up Delta-T (TT minus UT1), in days, at Julian Day `jd`.
    ///
    /// Installed tier priority and inclusive composite file-header bounds
    /// select the first covering tier with a non-empty table carrier. A
    /// covering tier without a carrier is skipped only after present table
    /// headers pass structural checks; exhaustion is [`Error::LebCorrupt`].
    /// Selected table evaluation keeps its clamping,
    /// right-biased search, interpolation order and stored IEEE words.
    ///
    /// Once a carrier is selected, its result is returned unchanged. Table
    /// corruption and [`Error::NonFiniteArgument`] must not select another
    /// tier. The immutable concrete readers' recoverable lookup failure is
    /// table absence, which the composite records separately from evaluation
    /// errors. The shared [`Error::LebCorrupt`] variant also represents
    /// malformed selected tables, so it cannot classify absence by itself.
    /// NaN `jd` selects no tier; an infinity admitted by stored header bounds
    /// retains the selected reader's strict non-finite-input refusal.
    pub fn delta_t(&self, jd: f64) -> Result<f64> {
        for (_, composite) in &self.tiers {
            let (start, end) = composite.jd_range();
            if start <= jd && jd <= end {
                if let Some(value) = composite.delta_t_for_tier(jd)? {
                    return Ok(value);
                }
            }
        }
        Err(Error::LebCorrupt {
            what: format!("LEB tiered Delta-T (no coverage at JD {jd})"),
            offset: 0,
        })
    }

    /// Star entry from the first installed tier carrying it. Exhaustion is [`Error::LebStarNotFound`].
    pub fn get_star(&self, star_id: i32) -> Result<StarEntry> {
        for (_, composite) in &self.tiers {
            match composite.get_star(star_id) {
                Ok(entry) => return Ok(entry),
                Err(Error::LebStarNotFound { .. }) => continue,
                Err(other) => return Err(other),
            }
        }
        Err(Error::LebStarNotFound { star_id })
    }
}

/// Tiers eligible under the configured active tier: `base` admits only
/// base, `medium` admits base+medium, `extended` admits all three. Source:
/// release `state._discover_reviewed_leb_tier_cores` (`priority` slice up
/// to and including the active tier).
pub fn eligible_tiers(active: LebTier) -> &'static [LebTier] {
    match active {
        LebTier::Base => &[LebTier::Base],
        LebTier::Medium => &[LebTier::Base, LebTier::Medium],
        LebTier::Extended => &[LebTier::Base, LebTier::Medium, LebTier::Extended],
    }
}

impl CoverageReport {
    /// True when `jd` falls inside the reported stored interval
    /// (inclusive bounds).
    pub fn contains(&self, jd: f64) -> bool {
        interval_contains(self.start, self.end, jd)
    }
}

impl<'a> CompositeIndex<'a> {
    /// Coverage for `body_id` from the serving (first-wins) reader.
    /// Fictitious ids always answer `None`: no file provides coverage
    /// for them, even one carrying retired channels.
    pub fn coverage_for(&self, body_id: i32) -> Option<CoverageReport> {
        if is_fictitious_body(body_id) {
            return None;
        }
        let (start, end) = self.body_coverage(body_id)?;
        Some(CoverageReport {
            start,
            end,
            selected: self.serving_index(body_id),
        })
    }
}

impl<'a> TieredIndex<'a> {
    /// Coverage for `body_id`, optionally resolved for a target date.
    /// A covered date reports the selected tier's own interval with its
    /// priority position; a date outside every stored interval (or no
    /// date at all) reports the union interval with `selected` unset.
    /// Fictitious ids always answer `None`.
    pub fn coverage_for(&self, body_id: i32, jd: Option<f64>) -> Option<CoverageReport> {
        if is_fictitious_body(body_id) {
            return None;
        }
        if let Some(jd) = jd {
            if let Some(pos) = self.selected(body_id, jd) {
                let (start, end) = self.tiers[pos].1.body_coverage(body_id)?;
                return Some(CoverageReport {
                    start,
                    end,
                    selected: Some(pos),
                });
            }
        }
        let (start, end) = self.body_coverage(body_id)?;
        Some(CoverageReport {
            start,
            end,
            selected: None,
        })
    }

    /// Check this index's installed tiers against the active
    /// configuration tier (see [`check_tier_subset`]).
    pub fn check_eligible(&self, active: LebTier) -> Result<()> {
        let installed: Vec<LebTier> = self.tiers.iter().map(|(tier, _)| *tier).collect();
        check_tier_subset(&installed, active)
    }

    /// Sealed-mode disposition of `body_id` at `jd` over this index's
    /// union intervals (see [`sealed_outcome`]). Evaluation itself is
    /// unchanged: covered dates evaluate, and every miss surfaces the
    /// disposition's typed error (or the evaluator's own
    /// [`Error::LebBodyNotFound`] when no channel exists) — never a
    /// silent substitution.
    pub fn sealed_precheck(&self, body_id: i32, jd: f64) -> SealedOutcome {
        sealed_outcome(
            body_id,
            jd,
            self.body_coverage(body_id),
            self.body_coverage(SUN_BODY_ID),
        )
    }
}

impl BodyCoverage {
    /// Assemble a coverage record over an already-resolved stored
    /// interval. Fictitious ids answer `None`; the caller
    /// supplies the serving file identity explicitly (see the data-inventory
    /// module notes).
    pub fn new(
        body_id: i32,
        bounds: (f64, f64),
        data_file: Option<&str>,
        reviewed: bool,
    ) -> Option<Self> {
        if is_fictitious_body(body_id) {
            return None;
        }
        let class = precision_class(body_id, data_file, reviewed);
        Some(Self {
            body_id,
            precision_class: class,
            jd_start: bounds.0,
            jd_end: bounds.1,
            data_file: data_file.map(str::to_string),
            group: group_from_path(data_file),
            reviewed,
        })
    }

    /// Body identifier in the public numbering.
    pub fn body_id(&self) -> i32 {
        self.body_id
    }
    /// Numerical origin of the serving artifact (see [`PrecisionClass`]).
    pub fn precision_class(&self) -> PrecisionClass {
        self.precision_class
    }
    /// Calculation source of every record built here.
    pub fn source(&self) -> &'static str {
        "LEB"
    }
    /// Start of the stored coefficient interval (Julian Day).
    pub fn jd_start(&self) -> f64 {
        self.jd_start
    }
    /// End of the stored coefficient interval (Julian Day).
    pub fn jd_end(&self) -> f64 {
        self.jd_end
    }
    /// Serving file path, or `None` for a date-less tiered union (the
    /// release refuses to mislabel the priority fallback file as the
    /// source of an out-of-range result).
    pub fn data_file(&self) -> Option<&str> {
        self.data_file.as_deref()
    }
    /// Canonical group inferred from [`BodyCoverage::data_file`].
    pub fn group(&self) -> Option<&str> {
        self.group.as_deref()
    }
    /// Whether the serving file is a reviewed manifest artifact.
    pub fn reviewed(&self) -> bool {
        self.reviewed
    }
    /// Whether `jd` is inside this body's closed stored interval (the
    /// release `BodyCoverage.contains` rule).
    pub fn contains(&self, jd: f64) -> bool {
        self.jd_start <= jd && jd <= self.jd_end
    }
}

/// Build the [`FileInfo`] inventory row for one decoded container. Ports
/// release `_reader_file_info` over explicit inputs: `path` is the full
/// container path, `size_bytes` the caller-observed length (`None` when
/// stating the file failed), `reviewed` the manifest-verified flag, and
/// `index` the decoded header/body directory.
pub fn inspect_index(
    path: &str,
    size_bytes: Option<u64>,
    reviewed: bool,
    index: &FileIndex,
) -> FileInfo {
    let name = path.rsplit('/').next().unwrap_or(path).to_string();
    let mut ids: Vec<i32> = index.bodies().iter().map(BodyEntry::body_id).collect();
    ids.sort_unstable();
    let mut bodies = Vec::with_capacity(ids.len());
    for body_id in ids {
        if let Some((start, end)) = index.body_coverage(body_id) {
            bodies.push(BodyRange {
                body_id,
                jd_start: start,
                jd_end: end,
            });
        }
    }
    FileInfo {
        name,
        group: group_from_path(Some(path)),
        path: path.to_string(),
        size_bytes,
        reviewed,
        bodies,
    }
}

impl FileIndex {
    /// Coverage record for `body_id` from this container (the single-file
    /// shape of release `get_reader_body_coverage`): fictitious ids
    /// answer `None`, and `data_file` / `reviewed` label the record like
    /// the release `.path` / `._manifest_verified` attributes.
    pub fn coverage_record(
        &self,
        body_id: i32,
        data_file: Option<&str>,
        reviewed: bool,
    ) -> Option<BodyCoverage> {
        let bounds = self.body_coverage(body_id)?;
        BodyCoverage::new(body_id, bounds, data_file, reviewed)
    }
}

impl<'a> CompositeIndex<'a> {
    /// Coverage record for `body_id` from the first-wins serving
    /// constituent (the composite shape of release
    /// `get_reader_body_coverage`): `files` runs parallel to the
    /// construction order and names each constituent's serving file; a
    /// short slice labels the record unlabeled rather than failing.
    /// Fictitious ids answer `None`.
    pub fn coverage_record(&self, files: &[ServingFile<'_>], body_id: i32) -> Option<BodyCoverage> {
        if is_fictitious_body(body_id) {
            return None;
        }
        let bounds = self.body_coverage(body_id)?;
        let pos = self.serving_index(body_id)?;
        let serving = files.get(pos).copied().unwrap_or(ServingFile::unlabeled());
        BodyCoverage::new(body_id, bounds, serving.path, serving.reviewed)
    }
}

impl<'a> TieredIndex<'a> {
    /// Coverage record for `body_id`, optionally resolved for a target
    /// date (the tiered shape of release `get_reader_body_coverage`). A
    /// covered date reports the selected tier's own interval with the
    /// serving file labeled; a date outside every stored interval — or no
    /// date at all — reports the union interval unlabeled (`data_file`
    /// `None`), never mislabeling the priority fallback file as the
    /// source. `file` resolves `(tier, constituent position)` to serving
    /// file identity; unknown entries resolve unlabeled. Fictitious ids
    /// answer `None`.
    pub fn coverage_record<'s, F>(
        &self,
        body_id: i32,
        jd: Option<f64>,
        mut file: F,
    ) -> Option<BodyCoverage>
    where
        F: FnMut(LebTier, usize) -> Option<ServingFile<'s>>,
    {
        if is_fictitious_body(body_id) {
            return None;
        }
        if let Some(jd) = jd {
            if let Some(pos) = self.selected(body_id, jd) {
                let (tier, composite) = &self.tiers[pos];
                let bounds = composite.body_coverage(body_id)?;
                let file_pos = composite.serving_index(body_id).unwrap_or(0);
                let serving = file(*tier, file_pos).unwrap_or(ServingFile::unlabeled());
                return BodyCoverage::new(body_id, bounds, serving.path, serving.reviewed);
            }
        }
        let bounds = self.body_coverage(body_id)?;
        // Union case: keep the highest-priority candidate's reviewed flag
        //  while leaving the file unlabeled.
        let reviewed = self
            .serving(body_id, None)
            .and_then(|pos| {
                let (tier, composite) = &self.tiers[pos];
                let file_pos = composite.serving_index(body_id).unwrap_or(0);
                file(*tier, file_pos)
            })
            .is_some_and(|serving| serving.reviewed);
        BodyCoverage::new(body_id, bounds, None, reviewed)
    }
}

impl TierAssets {
    /// Active tier this report covers.
    pub fn tier(&self) -> LebTier {
        self.tier
    }
    /// Per-file rows in tier-major, group-minor requirement order.
    pub fn files(&self) -> &[AssetStatus] {
        &self.files
    }
    /// True only when every required asset is ready.
    pub fn ready(&self) -> bool {
        self.files.iter().all(AssetStatus::is_ready)
    }
    /// Names of the required assets with no observed bytes.
    pub fn missing(&self) -> Vec<&str> {
        self.files
            .iter()
            .filter(|file| matches!(file.state, AssetState::Missing))
            .map(|file| file.name())
            .collect()
    }
    /// Names of the required assets whose bytes failed to decode.
    pub fn unreadable(&self) -> Vec<&str> {
        self.files
            .iter()
            .filter(|file| matches!(file.state, AssetState::Unreadable { .. }))
            .map(|file| file.name())
            .collect()
    }
}

impl FileIndex {
    /// Unmerged byte ranges overlapping the supplied TT Julian-day bounds.
    ///
    /// LEB1 uses every overlapping body's stored segment grid. For degree
    /// `d`, component count `c`, and segment count `n`, one segment occupies
    /// `(d + 1) * c * 8` bytes. The first index is
    /// `max(0, trunc((jd_start - body_start) / interval))`; the last is
    /// `min(n - 1, trunc((jd_end - body_start) / interval))`. Truncation is
    /// toward zero. A body contributes nothing when its first index exceeds
    /// its last. Bounds are used in supplied order, without a reversal guard.
    /// LEB2 v1 contributes overlapping body blobs; v2 contributes overlapping
    /// chunk blobs. Ranges remain per-image; use [`merge_ranges`] before advice.
    ///
    /// Request bounds and both LEB1 segment quotients must be finite;
    /// otherwise [`Error::NonFiniteArgument`] identifies the failing value.
    /// Bounds are checked before format dispatch. Each body's overlap and
    /// checked segment-byte arithmetic precede its first quotient, which is
    /// checked and converted before evaluating the last quotient. No ranges
    /// escape an error, so callers perform no page advice for that request.
    /// Byte-count/offset overflow is [`Error::LebCorrupt`]. Malformed LEB2
    /// chunk directories are skipped by this advisory operation.
    pub(crate) fn warm_ranges_source(
        &self,
        bytes: LebSource<'_>,
        jd_start: f64,
        jd_end: f64,
    ) -> Result<Vec<(u64, u64)>> {
        if !jd_start.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_start });
        }
        if !jd_end.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_end });
        }
        let mut ranges = Vec::new();
        match self.header.format {
            LebFormat::Leb1 => {
                for entry in &self.bodies {
                    if entry.jd_end < jd_start || entry.jd_start > jd_end {
                        continue;
                    }
                    let seg_bytes = u64::from(entry.degree)
                        .checked_add(1)
                        .and_then(|v| v.checked_mul(u64::from(entry.components)))
                        .and_then(|v| v.checked_mul(8))
                        .ok_or_else(|| Error::LebCorrupt {
                            what: format!("LEB warm range (body {})", entry.body_id),
                            offset: usize::try_from(entry.data_offset).unwrap_or(usize::MAX),
                        })?;
                    let count = entry.segment_count as i64;
                    let first_quotient = (jd_start - entry.jd_start) / entry.interval_days;
                    if !first_quotient.is_finite() {
                        return Err(Error::NonFiniteArgument {
                            value: first_quotient,
                        });
                    }
                    let first = first_quotient.trunc() as i64;
                    let last_quotient = (jd_end - entry.jd_start) / entry.interval_days;
                    if !last_quotient.is_finite() {
                        return Err(Error::NonFiniteArgument {
                            value: last_quotient,
                        });
                    }
                    let last = last_quotient.trunc() as i64;
                    let first = first.max(0);
                    let last = last.min(count - 1);
                    if first > last {
                        continue;
                    }
                    let offset = entry
                        .data_offset
                        .checked_add((first as u64).checked_mul(seg_bytes).ok_or_else(|| {
                            Error::LebCorrupt {
                                what: format!("LEB warm range (body {})", entry.body_id),
                                offset: usize::try_from(entry.data_offset).unwrap_or(usize::MAX),
                            }
                        })?)
                        .ok_or_else(|| Error::LebCorrupt {
                            what: format!("LEB warm range (body {})", entry.body_id),
                            offset: usize::try_from(entry.data_offset).unwrap_or(usize::MAX),
                        })?;
                    let size = ((last - first + 1) as u64)
                        .checked_mul(seg_bytes)
                        .ok_or_else(|| Error::LebCorrupt {
                            what: format!("LEB warm range (body {})", entry.body_id),
                            offset: usize::try_from(entry.data_offset).unwrap_or(usize::MAX),
                        })?;
                    ranges.push((offset, size));
                }
            }
            LebFormat::Leb2V1 => {
                for entry in &self.bodies {
                    if entry.jd_end < jd_start || entry.jd_start > jd_end {
                        continue;
                    }
                    let Some(compressed) = entry.compressed_size else {
                        return Err(Error::LebCorrupt {
                            what: format!("LEB2 warm range (body {})", entry.body_id),
                            offset: usize::try_from(entry.data_offset).unwrap_or(usize::MAX),
                        });
                    };
                    ranges.push((entry.data_offset, compressed));
                }
            }
            LebFormat::Leb2V2 => {
                for entry in &self.bodies {
                    let chunks = match decode::decode_chunk_index_source(
                        bytes,
                        entry.data_offset,
                        entry.body_id,
                    ) {
                        Ok((_, chunks)) => chunks,
                        Err(
                            error @ (Error::ResourceLimit { .. } | Error::EphemerisLoaderIo { .. }),
                        ) => return Err(error),
                        // Malformed unrelated directories cannot improve
                        // advisory warming. Explicit resource/I/O refusals,
                        // however, must not look like absent metadata.
                        Err(_) => continue,
                    };
                    for chunk in chunks.iter() {
                        if chunk.jd_end < jd_start || chunk.jd_start > jd_end {
                            continue;
                        }
                        ranges.push((chunk.blob_offset, chunk.compressed_size));
                    }
                }
            }
        }
        Ok(ranges)
    }
}

impl<'a> CompositeIndex<'a> {
    /// Prefault every constituent image over `[jd_start, jd_end]`,
    /// returning the merged range count per constituent in order (the
    /// release composite `warm` fans out to one `warm` call per reader;
    /// ranges stay per-image because offsets are per-file address
    /// spaces and must never merge across files).
    pub fn warm(&self, jd_start: f64, jd_end: f64) -> Result<Vec<usize>> {
        if !self.entries.files().is_empty() {
            return self
                .entries
                .files()
                .iter()
                .map(|file| file.warm(jd_start, jd_end))
                .collect();
        }
        self.entries
            .iter()
            .map(|(index, bytes)| {
                let LebSource::Slice(bytes) = bytes else {
                    unreachable!("raw views use contiguous slices")
                };
                let ranges = index.warm_ranges(bytes, jd_start, jd_end)?;
                let merged = merge_ranges(ranges, bytes.len() as u64);
                prefault(bytes, &merged);
                Ok(merged.len())
            })
            .collect()
    }

    /// Release coefficients and advise mapped-page release for file-backed views.
    /// Raw borrowed slices have no OS ownership, so are never passed to madvise.
    pub fn cool(&self) {
        for cache in &self.caches {
            cache.clear();
        }
        for file in self.entries.files() {
            file.release_file_pages();
        }
    }
}

impl<'a> TieredIndex<'a> {
    /// Prefault every tier composite over `[jd_start, jd_end]`,
    /// returning per-tier per-constituent merged range counts in
    /// priority order.
    pub fn warm(&self, jd_start: f64, jd_end: f64) -> Result<Vec<Vec<usize>>> {
        self.tiers
            .iter()
            .map(|(_, composite)| composite.warm(jd_start, jd_end))
            .collect()
    }

    /// Release each tier's caches and mapped pages (see [`CompositeIndex::cool`]).
    pub fn cool(&self) {
        for (_, composite) in &self.tiers {
            composite.cool();
        }
    }
}

#[cfg(test)]
mod nutation_derivative_storage_tests {
    use super::*;

    fn fixture(degree: u32) -> (FileIndex, Vec<u8>, Vec<f64>) {
        let mut bytes = Vec::new();
        for value in [0.0_f64, 16.0, 16.0] {
            bytes.extend(value.to_le_bytes());
        }
        for value in [degree, 2, 1, 0] {
            bytes.extend(value.to_le_bytes());
        }
        let coefficients: Vec<_> = (0..=degree)
            .map(|i| {
                if i == 0 {
                    f64::NAN
                } else {
                    (i as f64 - 3.0) * 2.0_f64.powi(-20)
                }
            })
            .collect();
        for _ in 0..2 {
            for coefficient in &coefficients {
                bytes.extend(coefficient.to_le_bytes());
            }
        }
        let index = FileIndex {
            header: FileHeader {
                magic: MAGIC_LEB1,
                version: VERSION_LEB1,
                section_count: 1,
                body_count: 0,
                jd_start: 0.0,
                jd_end: 16.0,
                generation_epoch: 0.0,
                flags: 0,
                format: LebFormat::Leb1,
            },
            sections: vec![SectionEntry {
                section_id: SECTION_NUTATION,
                offset: 0,
                size: bytes.len() as u64,
            }],
            bodies: vec![],
        };
        (index, bytes, coefficients)
    }

    #[test]
    fn inline_and_heap_nutation_derivatives_keep_the_original_words() {
        for degree in [0, 1, 16, 63, 64, 65] {
            let (index, bytes, coefficients) = fixture(degree);
            for (jd, tau) in [(0.0, -1.0), (8.0, 0.0), (16.0, 1.0)] {
                let expected =
                    chebyshev_value(&chebyshev_deriv_coeffs(&coefficients), tau).unwrap() * 0.125;
                let actual = index
                    .nutation_rate_source(LebSource::Slice(&bytes), jd)
                    .unwrap()
                    .unwrap();
                assert_eq!(actual.0.to_bits(), expected.to_bits());
                assert_eq!(actual.1.to_bits(), expected.to_bits());
            }
        }
    }

    struct Temporary(std::path::PathBuf);
    impl Drop for Temporary {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }

    #[test]
    fn derivative_storage_keeps_working_admission_before_coefficient_checks() {
        use std::sync::atomic::{AtomicU64, Ordering};
        static NEXT: AtomicU64 = AtomicU64::new(0);
        for degree in [1, 64, 65] {
            let (index, mut bytes, _) = fixture(degree);
            bytes[NUTATION_HEADER_SIZE + 8..NUTATION_HEADER_SIZE + 16]
                .copy_from_slice(&f64::INFINITY.to_le_bytes());
            let path = Temporary(std::env::temp_dir().join(format!(
                "nutation-derivative-storage-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            )));
            std::fs::write(&path.0, &bytes).unwrap();
            let working_bytes = 1 << 20;
            let resources = resources::IndexedResources::new(IndexedLimits {
                working_bytes,
                ..Default::default()
            });
            let file = access::PositionalFile::open(&path.0, resources.clone()).unwrap();
            let source = LebSource::Positional(&file);
            let ordinary = index.nutation_rate_source(source, 8.0).unwrap_err();
            assert!(matches!(ordinary, Error::NonFiniteArgument { value } if value.is_infinite()));
            let resident_working = resources.working.used();
            // The row is now cached. Keep exactly one byte too little for the
            // original reservation, including when the derivative uses stack.
            let required = (degree as usize + 1) * std::mem::size_of::<f64>();
            let held = resources
                .working
                .reserve(working_bytes - resident_working - required + 1)
                .unwrap();
            assert!(matches!(
                index.nutation_rate_source(source, 8.0),
                Err(Error::ResourceLimit { required, limit, .. })
                    if required == working_bytes + 1 && limit == working_bytes
            ));
            drop(held);
            assert_eq!(resources.working.used(), resident_working);
            let restored = index.nutation_rate_source(source, 8.0).unwrap_err();
            assert_eq!(format!("{restored:?}"), format!("{ordinary:?}"));
            assert_eq!(resources.working.used(), resident_working);
        }
    }
}
