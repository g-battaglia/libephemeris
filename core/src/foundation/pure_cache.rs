//! Bounded reuse of exact-input numerical results, without resource ownership.
//!
//! Values are specified by their binary64 operands or an immutable image token
//! and exact selected-record/time operands. Admission, callbacks, acquisition
//! and tracing remain with callers. Fixed numerical blocks share one process
//! budget; no image, Session, callback or request context is retained.

use crate::positions::{CartesianWordParts, FrameMatrices};
use std::num::NonZeroU64;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Mutex, OnceLock};

/// Exact native polynomial identity, without a pointer or resource owner.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) struct RawWordStateKey {
    image: NonZeroU64,
    body: i32,
    segment: u32,
    tt: [u64; 2],
}

impl RawWordStateKey {
    pub(crate) fn new(image: NonZeroU64, body: i32, segment: u32, tt: [u64; 2]) -> Self {
        Self {
            image,
            body,
            segment,
            tt,
        }
    }

    fn words(self) -> [u64; 4] {
        [
            self.image.get(),
            u64::from(self.body as u32) | (u64::from(self.segment) << 32),
            self.tt[0],
            self.tt[1],
        ]
    }
}

// Tokens never wrap or identify an address. Exhaustion only disables reuse.
static NEXT_WORD_IMAGE: AtomicU64 = AtomicU64::new(1);

fn issue_word_image_identity(next: &AtomicU64) -> Option<NonZeroU64> {
    next.fetch_update(Ordering::Relaxed, Ordering::Relaxed, |value| {
        if value == 0 {
            None
        } else {
            value.checked_add(1)
        }
    })
    .ok()
    .and_then(NonZeroU64::new)
}

/// Allocate only when an immutable native image binds its own cache.
pub(crate) fn new_word_image_identity() -> Option<NonZeroU64> {
    issue_word_image_identity(&NEXT_WORD_IMAGE)
}

/// The ordinary twelve-cusp layout, without an owned output allocation.
#[derive(Clone, Copy)]
pub(crate) struct HouseGeometry {
    pub(crate) cusps: [f64; 12],
    pub(crate) angles: [f64; 8],
}

/// Row-major output-frame matrix and its derivative per TT day.
pub(crate) type FrameRate = ([[f64; 3]; 3], [[f64; 3]; 3]);

/// Complete dimensionless matrix words and their rates per TT day.
#[derive(Clone, Copy)]
pub(crate) struct AnalyticFrame {
    pub(crate) matrix: [[[f64; 2]; 3]; 3],
    pub(crate) rate: [[[f64; 2]; 3]; 3],
}

const BLOCK_WORDS: usize = 16;
const BLOCKS: usize = 432;
const ENTRIES: usize = 384;
const BUCKETS: usize = 1024;
const MASK: usize = BUCKETS - 1;
const DOMAINS: usize = 8;
const NONE: u16 = u16::MAX;

/// Each domain has one fixed codec and the original independent entry limit.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[repr(u8)]
enum Kind {
    Empty,
    Frame,
    OwnedApparent,
    RawWord,
    FrameRate,
    AnalyticFrame,
    Houses,
    ReferenceNutation,
    ErfaNutation,
}

impl Kind {
    fn index(self) -> usize {
        self as usize - 1
    }

    fn key_words(self) -> usize {
        match self {
            Self::Frame => 3,
            Self::OwnedApparent => 9,
            Self::RawWord => 4,
            Self::FrameRate | Self::AnalyticFrame => 6,
            Self::Houses => 5,
            Self::ReferenceNutation | Self::ErfaNutation => 1,
            Self::Empty => unreachable!("empty numerical descriptor"),
        }
    }

    fn payload_words(self) -> usize {
        match self {
            Self::Frame => 49,
            Self::OwnedApparent => 15,
            Self::RawWord => 12,
            Self::FrameRate => 18,
            Self::AnalyticFrame => 36,
            Self::Houses => 20,
            Self::ReferenceNutation | Self::ErfaNutation => 2,
            Self::Empty => unreachable!("empty numerical descriptor"),
        }
    }

    fn blocks(self) -> usize {
        (self.key_words() + self.payload_words()).div_ceil(BLOCK_WORDS)
    }

    fn limit(self) -> u16 {
        match self {
            Self::Frame => 31,
            Self::RawWord => 64,
            Self::FrameRate | Self::AnalyticFrame => 16,
            Self::Houses | Self::ReferenceNutation | Self::ErfaNutation => 128,
            // Only this domain may use otherwise unoccupied global capacity.
            Self::OwnedApparent => ENTRIES as u16,
            Self::Empty => unreachable!("empty numerical descriptor"),
        }
    }
}

/// Links are replacement metadata, not resource references or access history.
/// A 32-bit hash is only an index hint: the complete kind/key decides equality.
#[derive(Clone, Copy)]
struct Descriptor {
    hash: u32,
    first: u16,
    global_prev: u16,
    global_next: u16,
    domain_prev: u16,
    domain_next: u16,
    kind: Kind,
}

impl Descriptor {
    const EMPTY: Self = Self {
        hash: 0,
        first: NONE,
        global_prev: NONE,
        global_next: NONE,
        domain_prev: NONE,
        domain_next: NONE,
        kind: Kind::Empty,
    };
}

struct PureCache {
    blocks: [[u64; BLOCK_WORDS]; BLOCKS],
    // Entries can span nonadjacent free blocks; fragmentation never requires
    // an extra eviction or a temporary payload allocation.
    block_next: [u16; BLOCKS],
    entries: [Descriptor; ENTRIES],
    buckets: [u16; BUCKETS], // descriptor index + 1; zero is empty
    domain_oldest: [u16; DOMAINS],
    domain_newest: [u16; DOMAINS],
    domain_count: [u16; DOMAINS],
    oldest: u16,
    newest: u16,
    free_entry: u16,
    free_block: u16,
    free_blocks: u16,
    len: u16,
}

impl PureCache {
    const fn new() -> Self {
        let mut entries = [Descriptor::EMPTY; ENTRIES];
        let mut index = 0;
        while index < ENTRIES {
            entries[index].global_next = if index + 1 < ENTRIES {
                (index + 1) as u16
            } else {
                NONE
            };
            index += 1;
        }
        let mut block_next = [NONE; BLOCKS];
        index = 0;
        while index + 1 < BLOCKS {
            block_next[index] = (index + 1) as u16;
            index += 1;
        }
        Self {
            blocks: [[0; BLOCK_WORDS]; BLOCKS],
            block_next,
            entries,
            buckets: [0; BUCKETS],
            domain_oldest: [NONE; DOMAINS],
            domain_newest: [NONE; DOMAINS],
            domain_count: [0; DOMAINS],
            oldest: NONE,
            newest: NONE,
            free_entry: 0,
            free_block: 0,
            free_blocks: BLOCKS as u16,
            len: 0,
        }
    }

    fn hash(kind: Kind, key: &[u64]) -> u32 {
        let mut hash = (kind as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15);
        for &word in key {
            hash ^= word.wrapping_add(0x9e37_79b9_7f4a_7c15);
            hash = hash.rotate_left(27).wrapping_mul(0xff51_afd7_ed55_8ccd);
        }
        hash ^= hash >> 33;
        hash = hash.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
        hash ^= hash >> 33;
        // Hash truncation affects collisions only, never key or payload words.
        (hash ^ (hash >> 32)) as u32
    }

    fn home(hash: u32) -> usize {
        hash as usize & MASK
    }

    fn find_bucket(&self, kind: Kind, key: &[u64], hash: u32) -> Option<usize> {
        let mut bucket = Self::home(hash);
        for _ in 0..BUCKETS {
            let slot = usize::from(self.buckets[bucket].checked_sub(1)?);
            let entry = self.entries[slot];
            // Every key fits wholly inside the first block.
            if entry.kind == kind
                && entry.hash == hash
                && self.blocks[usize::from(entry.first)][..key.len()] == *key
            {
                return Some(bucket);
            }
            bucket = (bucket + 1) & MASK;
        }
        unreachable!("bounded numerical descriptors cannot fill the index")
    }

    fn get(&self, kind: Kind, key: &[u64], payload: &mut [u64]) -> bool {
        assert_eq!(key.len(), kind.key_words());
        assert_eq!(payload.len(), kind.payload_words());
        let Some(bucket) = self.find_bucket(kind, key, Self::hash(kind, key)) else {
            return false;
        };
        let entry = self.entries[usize::from(self.buckets[bucket] - 1)];
        let mut block = entry.first;
        let mut word = kind.key_words();
        for output in payload {
            if word == BLOCK_WORDS {
                block = self.block_next[usize::from(block)];
                word = 0;
            }
            *output = self.blocks[usize::from(block)][word];
            word += 1;
        }
        true
    }

    fn remove_bucket(&mut self, mut hole: usize) {
        self.buckets[hole] = 0;
        let mut current = (hole + 1) & MASK;
        for _ in 0..BUCKETS {
            let Some(slot) = self.buckets[current].checked_sub(1) else {
                return;
            };
            let home = Self::home(self.entries[usize::from(slot)].hash);
            if (hole.wrapping_sub(home) & MASK) < (current.wrapping_sub(home) & MASK) {
                self.buckets[hole] = self.buckets[current];
                self.buckets[current] = 0;
                hole = current;
            }
            current = (current + 1) & MASK;
        }
        unreachable!("bounded numerical descriptors cannot fill the index")
    }

    fn evict(&mut self, slot: u16) {
        let entry = self.entries[usize::from(slot)];
        let domain = entry.kind.index();
        let mut bucket = Self::home(entry.hash);
        while self.buckets[bucket] != slot + 1 {
            bucket = (bucket + 1) & MASK;
        }
        self.remove_bucket(bucket);
        if entry.global_prev == NONE {
            self.oldest = entry.global_next;
        } else {
            self.entries[usize::from(entry.global_prev)].global_next = entry.global_next;
        }
        if entry.global_next == NONE {
            self.newest = entry.global_prev;
        } else {
            self.entries[usize::from(entry.global_next)].global_prev = entry.global_prev;
        }
        if entry.domain_prev == NONE {
            self.domain_oldest[domain] = entry.domain_next;
        } else {
            self.entries[usize::from(entry.domain_prev)].domain_next = entry.domain_next;
        }
        if entry.domain_next == NONE {
            self.domain_newest[domain] = entry.domain_prev;
        } else {
            self.entries[usize::from(entry.domain_next)].domain_prev = entry.domain_prev;
        }
        let mut block = entry.first;
        for _ in 0..entry.kind.blocks() {
            let next = self.block_next[usize::from(block)];
            self.block_next[usize::from(block)] = self.free_block;
            self.free_block = block;
            block = next;
        }
        self.free_blocks += entry.kind.blocks() as u16;
        self.domain_count[domain] -= 1;
        self.len -= 1;
        self.entries[usize::from(slot)] = Descriptor {
            global_next: self.free_entry,
            ..Descriptor::EMPTY
        };
        self.free_entry = slot;
    }

    fn insert(&mut self, kind: Kind, key: &[u64], payload: &[u64]) {
        assert_eq!(key.len(), kind.key_words());
        assert_eq!(payload.len(), kind.payload_words());
        let hash = Self::hash(kind, key);
        if self.find_bucket(kind, key, hash).is_some() {
            return;
        }
        let domain = kind.index();
        if self.domain_count[domain] == kind.limit() {
            self.evict(self.domain_oldest[domain]);
        }
        while usize::from(self.free_blocks) < kind.blocks() || usize::from(self.len) == ENTRIES {
            self.evict(self.oldest);
        }
        let slot = self.free_entry;
        self.free_entry = self.entries[usize::from(slot)].global_next;
        let first = self.free_block;
        let mut block = first;
        for index in 0..kind.blocks() {
            let next = self.block_next[usize::from(block)];
            self.free_block = next;
            if index + 1 == kind.blocks() {
                self.block_next[usize::from(block)] = NONE;
            } else {
                block = next;
            }
        }
        block = first;
        let mut word = 0;
        for &value in key.iter().chain(payload) {
            if word == BLOCK_WORDS {
                block = self.block_next[usize::from(block)];
                word = 0;
            }
            self.blocks[usize::from(block)][word] = value;
            word += 1;
        }
        self.free_blocks -= kind.blocks() as u16;
        self.entries[usize::from(slot)] = Descriptor {
            hash,
            first,
            kind,
            global_prev: self.newest,
            global_next: NONE,
            domain_prev: self.domain_newest[domain],
            domain_next: NONE,
        };
        if self.newest == NONE {
            self.oldest = slot;
        } else {
            self.entries[usize::from(self.newest)].global_next = slot;
        }
        if self.domain_newest[domain] == NONE {
            self.domain_oldest[domain] = slot;
        } else {
            self.entries[usize::from(self.domain_newest[domain])].domain_next = slot;
        }
        self.newest = slot;
        self.domain_newest[domain] = slot;
        self.domain_count[domain] += 1;
        self.len += 1;
        let mut bucket = Self::home(hash);
        while self.buckets[bucket] != 0 {
            bucket = (bucket + 1) & MASK;
        }
        self.buckets[bucket] = slot + 1;
    }
}

static CACHE: OnceLock<Mutex<Box<PureCache>>> = {
    // Keep the complete retained-state budget a compile-time invariant.
    assert!(RETAINED_BYTES <= 64 * 1024);
    OnceLock::new()
};
static NEXT_OWNED_APPARENT_ID: AtomicU64 = AtomicU64::new(1);
const RETAINED_BYTES: usize = std::mem::size_of::<PureCache>()
    + std::mem::size_of::<OnceLock<Mutex<Box<PureCache>>>>()
    + 2 * std::mem::size_of::<AtomicU64>();

fn owned_apparent_identity(counter: &AtomicU64) -> Option<NonZeroU64> {
    counter
        .fetch_update(Ordering::Relaxed, Ordering::Relaxed, |value| {
            value.checked_add(1)
        })
        .ok()
        .and_then(NonZeroU64::new)
}

/// New immutable Session producer. A clone may share it; rebuilding may not.
pub(crate) fn new_owned_apparent_identity() -> Option<NonZeroU64> {
    owned_apparent_identity(&NEXT_OWNED_APPARENT_ID)
}

fn with_cache<T>(f: impl FnOnce(&mut PureCache) -> T) -> T {
    let cache = CACHE.get_or_init(|| Mutex::new(Box::new(PureCache::new())));
    let mut guard = cache
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    f(&mut guard)
}

/// Query only after concrete Owned Session and current request admission.
pub(crate) fn owned_apparent(key: [u64; 9]) -> Option<[f64; 15]> {
    let mut words = [0; 15];
    with_cache(|cache| cache.get(Kind::OwnedApparent, &key, &mut words))
        .then(|| words.map(f64::from_bits))
}

/// Retain a completed physical state after the original core trace succeeds.
pub(crate) fn retain_owned_apparent(key: [u64; 9], value: [f64; 15]) {
    let words = value.map(f64::to_bits);
    with_cache(|cache| cache.insert(Kind::OwnedApparent, &key, &words));
}

pub(crate) fn frame(key: [u64; 3]) -> Option<FrameMatrices> {
    let mut words = [0; 49];
    with_cache(|cache| cache.get(Kind::Frame, &key, &mut words))
        .then(|| FrameMatrices::from_cache_words(words))
}

pub(crate) fn retain_frame(key: [u64; 3], value: FrameMatrices) {
    let words = value.cache_words();
    with_cache(|cache| cache.insert(Kind::Frame, &key, &words));
}

fn raw_words(value: CartesianWordParts) -> [u64; 12] {
    std::array::from_fn(|index| {
        if index < 6 {
            value.0[index / 2][index % 2].to_bits()
        } else {
            value.1[(index - 6) / 2][index % 2].to_bits()
        }
    })
}

fn raw_from_words(words: [u64; 12]) -> CartesianWordParts {
    let component = |offset| {
        std::array::from_fn(|axis| {
            [
                f64::from_bits(words[offset + 2 * axis]),
                f64::from_bits(words[offset + 2 * axis + 1]),
            ]
        })
    };
    (component(0), component(6))
}

/// Lookup after the original coefficient acquisition and admission; no owner kept.
pub(crate) fn raw_word_state(key: RawWordStateKey) -> Option<CartesianWordParts> {
    let mut words = [0; 12];
    with_cache(|cache| cache.get(Kind::RawWord, &key.words(), &mut words))
        .then(|| raw_from_words(words))
}

/// Retain only a successful finite raw state of the selected native polynomial.
pub(crate) fn retain_raw_word_state(key: RawWordStateKey, value: CartesianWordParts) {
    let words = raw_words(value);
    with_cache(|cache| cache.insert(Kind::RawWord, &key.words(), &words));
}

fn frame_rate_words(value: FrameRate) -> [u64; 18] {
    std::array::from_fn(|index| {
        if index < 9 {
            value.0[index / 3][index % 3].to_bits()
        } else {
            value.1[(index - 9) / 3][index % 3].to_bits()
        }
    })
}

fn frame_rate_from_words(words: [u64; 18]) -> FrameRate {
    let matrix = |offset| {
        std::array::from_fn(|row| {
            std::array::from_fn(|column| f64::from_bits(words[offset + 3 * row + column]))
        })
    };
    (matrix(0), matrix(9))
}

/// Reuse a matrix/rate only for identical original numerical operands.
pub(crate) fn frame_rate(key: [u64; 6]) -> Option<FrameRate> {
    let mut words = [0; 18];
    with_cache(|cache| cache.get(Kind::FrameRate, &key, &mut words))
        .then(|| frame_rate_from_words(words))
}

/// Retain one successful pure matrix/rate result; failures remain outside.
pub(crate) fn retain_frame_rate(key: [u64; 6], value: FrameRate) {
    let words = frame_rate_words(value);
    with_cache(|cache| cache.insert(Kind::FrameRate, &key, &words));
}

fn analytic_words(value: AnalyticFrame) -> [u64; 36] {
    std::array::from_fn(|index| {
        let word = index % 18;
        let component = if index < 18 { value.matrix } else { value.rate };
        component[word / 6][(word / 2) % 3][word % 2].to_bits()
    })
}

fn analytic_from_words(words: [u64; 36]) -> AnalyticFrame {
    let matrix = |offset| {
        std::array::from_fn(|row| {
            std::array::from_fn(|column| {
                std::array::from_fn(|lane| {
                    f64::from_bits(words[offset + 6 * row + 2 * column + lane])
                })
            })
        })
    };
    AnalyticFrame {
        matrix: matrix(0),
        rate: matrix(18),
    }
}

/// Sampled and analytic realizations remain distinct even for equal keys.
pub(crate) fn analytic_frame(key: [u64; 6]) -> Option<AnalyticFrame> {
    let mut words = [0; 36];
    with_cache(|cache| cache.get(Kind::AnalyticFrame, &key, &mut words))
        .then(|| analytic_from_words(words))
}

/// Computation, custom callbacks and admission remain outside the mutex.
pub(crate) fn retain_analytic_frame(key: [u64; 6], value: AnalyticFrame) {
    let words = analytic_words(value);
    with_cache(|cache| cache.insert(Kind::AnalyticFrame, &key, &words));
}

pub(crate) fn houses(key: [u64; 5]) -> Option<HouseGeometry> {
    let mut words = [0; 20];
    with_cache(|cache| cache.get(Kind::Houses, &key, &mut words)).then(|| HouseGeometry {
        cusps: std::array::from_fn(|i| f64::from_bits(words[i])),
        angles: std::array::from_fn(|i| f64::from_bits(words[i + 12])),
    })
}

pub(crate) fn retain_houses(key: [u64; 5], value: HouseGeometry) {
    let words: [u64; 20] = std::array::from_fn(|i| {
        if i < 12 {
            value.cusps[i].to_bits()
        } else {
            value.angles[i - 12].to_bits()
        }
    });
    with_cache(|cache| cache.insert(Kind::Houses, &key, &words));
}

pub(crate) fn reference_nutation(key: u64) -> Option<(f64, f64)> {
    let mut words = [0; 2];
    with_cache(|cache| cache.get(Kind::ReferenceNutation, &[key], &mut words))
        .then(|| (f64::from_bits(words[0]), f64::from_bits(words[1])))
}

pub(crate) fn retain_reference_nutation(key: u64, value: (f64, f64)) {
    let words = [value.0.to_bits(), value.1.to_bits()];
    with_cache(|cache| cache.insert(Kind::ReferenceNutation, &[key], &words));
}

pub(crate) fn erfa_nutation(key: u64) -> Option<(f64, f64)> {
    let mut words = [0; 2];
    with_cache(|cache| cache.get(Kind::ErfaNutation, &[key], &mut words))
        .then(|| (f64::from_bits(words[0]), f64::from_bits(words[1])))
}

pub(crate) fn retain_erfa_nutation(key: u64, value: (f64, f64)) {
    let words = [value.0.to_bits(), value.1.to_bits()];
    with_cache(|cache| cache.insert(Kind::ErfaNutation, &[key], &words));
}

#[cfg(test)]
mod test_support {
    use super::*;

    pub(super) const IEEE: [u64; 12] = [
        0,
        1 << 63,
        1,
        u64::MAX,
        0x7ff8_0000_0000_1234,
        0xfff8_0000_0000_5678,
        0x7ff0_0000_0000_0000,
        0xfff0_0000_0000_0000,
        0x0010_0000_0000_0000,
        0x000f_ffff_ffff_ffff,
        0x3ff0_0000_0000_0000,
        0xbff0_0000_0000_0000,
    ];

    pub(super) fn key(image: u64, body: i32, segment: u32, tt: [u64; 2]) -> RawWordStateKey {
        RawWordStateKey::new(NonZeroU64::new(image).unwrap(), body, segment, tt)
    }

    pub(super) fn get(cache: &PureCache, kind: Kind, key: &[u64]) -> Option<Vec<u64>> {
        let mut words = vec![0; kind.payload_words()];
        cache.get(kind, key, &mut words).then_some(words)
    }

    pub(super) fn put(cache: &mut PureCache, kind: Kind, key: &[u64], word: u64) {
        cache.insert(kind, key, &vec![word; kind.payload_words()]);
    }

    pub(super) fn valid(cache: &PureCache) {
        let mut seen_entries = [false; ENTRIES];
        let mut seen_blocks = [false; BLOCKS];
        let mut counts = [0_u16; DOMAINS];
        let mut previous = NONE;
        let mut slot = cache.oldest;
        while slot != NONE {
            let index = usize::from(slot);
            assert!(!seen_entries[index]);
            seen_entries[index] = true;
            let entry = cache.entries[index];
            assert_ne!(entry.kind, Kind::Empty);
            assert_eq!(entry.global_prev, previous);
            counts[entry.kind.index()] += 1;
            let mut block = entry.first;
            for _ in 0..entry.kind.blocks() {
                let index = usize::from(block);
                assert!(!seen_blocks[index]);
                seen_blocks[index] = true;
                block = cache.block_next[index];
            }
            assert_eq!(block, NONE);
            previous = slot;
            slot = entry.global_next;
        }
        assert_eq!(previous, cache.newest);
        assert_eq!(counts, cache.domain_count);
        assert_eq!(counts.iter().sum::<u16>(), cache.len);
        for kind in [
            Kind::Frame,
            Kind::OwnedApparent,
            Kind::RawWord,
            Kind::FrameRate,
            Kind::AnalyticFrame,
            Kind::Houses,
            Kind::ReferenceNutation,
            Kind::ErfaNutation,
        ] {
            let domain = kind.index();
            assert!(counts[domain] <= kind.limit());
            let mut previous = NONE;
            let mut slot = cache.domain_oldest[domain];
            let mut count = 0;
            while slot != NONE {
                let entry = cache.entries[usize::from(slot)];
                assert_eq!(entry.kind, kind);
                assert_eq!(entry.domain_prev, previous);
                assert!(seen_entries[usize::from(slot)]);
                previous = slot;
                slot = entry.domain_next;
                count += 1;
            }
            assert_eq!(previous, cache.domain_newest[domain]);
            assert_eq!(count, counts[domain]);
        }
        let mut block = cache.free_block;
        let mut free = 0;
        while block != NONE {
            let index = usize::from(block);
            assert!(!seen_blocks[index]);
            seen_blocks[index] = true;
            block = cache.block_next[index];
            free += 1;
        }
        assert_eq!(free, cache.free_blocks);
        assert!(seen_blocks.into_iter().all(|seen| seen));
        let mut slot = cache.free_entry;
        let mut free = 0;
        while slot != NONE {
            let index = usize::from(slot);
            assert!(!seen_entries[index]);
            seen_entries[index] = true;
            assert_eq!(cache.entries[index].kind, Kind::Empty);
            slot = cache.entries[index].global_next;
            free += 1;
        }
        assert_eq!(free + usize::from(cache.len), ENTRIES);
        assert!(seen_entries.into_iter().all(|seen| seen));
        let mut indexed = [false; ENTRIES];
        for &bucket in &cache.buckets {
            if let Some(slot) = bucket.checked_sub(1) {
                assert!(!indexed[usize::from(slot)]);
                indexed[usize::from(slot)] = true;
                assert_ne!(cache.entries[usize::from(slot)].kind, Kind::Empty);
            }
        }
        assert_eq!(
            indexed.into_iter().filter(|&seen| seen).count(),
            usize::from(cache.len)
        );
    }
}

#[cfg(test)]
mod tests {
    use super::test_support::*;
    use super::*;

    #[test]
    fn replacement_and_racing_insertions_keep_the_fixed_capacity() {
        let mut cache = PureCache::new();
        for word in 0..64 {
            put(&mut cache, Kind::RawWord, &[1, 0, word, 0], word);
        }
        let oldest = cache.oldest;
        let free = cache.free_blocks;
        put(&mut cache, Kind::RawWord, &[1, 0, 0, 0], 999);
        assert_eq!(cache.oldest, oldest);
        assert_eq!(cache.free_blocks, free);
        assert_eq!(get(&cache, Kind::RawWord, &[1, 0, 0, 0]), Some(vec![0; 12]));
        for word in 64..10000 {
            put(&mut cache, Kind::RawWord, &[1, 0, word, 0], word);
        }
        assert_eq!(cache.domain_count[Kind::RawWord.index()], 64);
        assert!(get(&cache, Kind::RawWord, &[1, 0, 9935, 0]).is_none());
        for word in 9936..10000 {
            assert_eq!(
                get(&cache, Kind::RawWord, &[1, 0, word, 0]),
                Some(vec![word; 12])
            );
        }
        valid(&cache);
    }

    #[test]
    fn exact_keys_preserve_signed_zero_and_single_word_changes() {
        let mut cache = PureCache::new();
        let key = [1.0_f64.to_bits(), 0.0_f64.to_bits(), (-0.0_f64).to_bits()];
        let mut payload = [0; 49];
        payload[0] = 7;
        cache.insert(Kind::Frame, &key, &payload);
        assert_eq!(get(&cache, Kind::Frame, &key), Some(payload.to_vec()));
        for index in 0..key.len() {
            let mut different = key;
            different[index] ^= 1;
            assert!(get(&cache, Kind::Frame, &different).is_none());
        }
        println!("retained pure-cache bytes: {RETAINED_BYTES}");
    }

    #[test]
    fn frame_rate_keys_keep_every_operand_and_signed_zero() {
        let mut cache = PureCache::new();
        let key = [2451545.0_f64.to_bits(), 258, 0, (-0.0_f64).to_bits(), 1, 2];
        put(&mut cache, Kind::FrameRate, &key, 7);
        for index in 0..key.len() {
            let mut changed = key;
            changed[index] ^= 1;
            assert!(get(&cache, Kind::FrameRate, &changed).is_none());
        }
        let mut changed = key;
        changed[3] = 0.0_f64.to_bits();
        assert!(get(&cache, Kind::FrameRate, &changed).is_none());
        assert_eq!(get(&cache, Kind::FrameRate, &key), Some(vec![7; 18]));
    }
}

#[cfg(test)]
mod analytic_retention_tests {
    use super::test_support::*;
    use super::*;

    #[test]
    fn analytic_retention_has_exact_keys_and_separate_numerical_models() {
        let mut cache = PureCache::new();
        let key = [2451545.0_f64.to_bits(), 258, 0, (-0.0_f64).to_bits(), 1, 2];
        let analytic = analytic_from_words(std::array::from_fn(|i| IEEE[i % IEEE.len()]));
        let sampled = frame_rate_from_words(std::array::from_fn(|i| IEEE[(i + 3) % IEEE.len()]));
        cache.insert(Kind::AnalyticFrame, &key, &analytic_words(analytic));
        cache.insert(Kind::FrameRate, &key, &frame_rate_words(sampled));
        assert_eq!(
            get(&cache, Kind::AnalyticFrame, &key),
            Some(analytic_words(analytic).to_vec())
        );
        assert_eq!(
            get(&cache, Kind::FrameRate, &key),
            Some(frame_rate_words(sampled).to_vec())
        );
        for index in 0..key.len() {
            let mut different = key;
            different[index] ^= 1;
            assert!(get(&cache, Kind::AnalyticFrame, &different).is_none());
        }
        for index in [2, 3] {
            let mut different = key;
            different[index] ^= 1_u64 << 63;
            assert!(get(&cache, Kind::AnalyticFrame, &different).is_none());
        }
        assert_eq!(
            analytic_words(analytic_from_words(analytic_words(analytic))),
            analytic_words(analytic)
        );
        println!(
            "analytic record words/blocks/global bytes: {}/{}/{RETAINED_BYTES}",
            Kind::AnalyticFrame.key_words() + Kind::AnalyticFrame.payload_words(),
            Kind::AnalyticFrame.blocks()
        );
    }

    #[test]
    fn analytic_retention_eviction_and_duplicate_completion_are_bounded() {
        let mut cache = PureCache::new();
        for word in 0..16 {
            put(&mut cache, Kind::AnalyticFrame, &[word; 6], word);
        }
        let oldest = cache.oldest;
        let newest = cache.newest;
        put(&mut cache, Kind::AnalyticFrame, &[4; 6], 999);
        assert_eq!((cache.oldest, cache.newest), (oldest, newest));
        assert_eq!(get(&cache, Kind::AnalyticFrame, &[4; 6]), Some(vec![4; 36]));
        for word in 16..1000 {
            put(&mut cache, Kind::AnalyticFrame, &[word; 6], word);
        }
        assert_eq!(cache.domain_count[Kind::AnalyticFrame.index()], 16);
        assert!(get(&cache, Kind::AnalyticFrame, &[983; 6]).is_none());
        for word in 984..1000 {
            assert_eq!(
                get(&cache, Kind::AnalyticFrame, &[word; 6]),
                Some(vec![word; 36])
            );
        }
        valid(&cache);
    }
}

#[cfg(test)]
mod raw_word_retention_tests {
    use super::test_support::*;
    use super::*;

    #[test]
    fn indexed_words_preserve_colliding_wraparound_chains_and_fifo() {
        let keys: Vec<_> = (0..)
            .map(|word| key(1, -3, 7, [word, 0]))
            .filter(|&key| PureCache::home(PureCache::hash(Kind::RawWord, &key.words())) == MASK)
            .take(80)
            .collect();
        let mut cache = PureCache::new();
        let mut original_fifo = std::collections::VecDeque::new();
        for (index, &key) in keys.iter().enumerate() {
            let words = [index as u64; 12];
            cache.insert(Kind::RawWord, &key.words(), &words);
            if original_fifo.len() == 64 {
                original_fifo.pop_front();
            }
            original_fifo.push_back((key, words));
            let order = (cache.oldest, cache.newest, cache.free_blocks);
            cache.insert(Kind::RawWord, &key.words(), &[u64::MAX; 12]);
            assert_eq!((cache.oldest, cache.newest, cache.free_blocks), order);
            for &probe in &keys {
                assert_eq!(
                    get(&cache, Kind::RawWord, &probe.words()),
                    original_fifo
                        .iter()
                        .find(|(key, _)| *key == probe)
                        .map(|(_, words)| words.to_vec())
                );
            }
            valid(&cache);
        }
    }

    #[test]
    fn indexed_words_keep_full_keys_and_all_payload_bits_through_replacement() {
        let mut state = 0x79b9_7f4a_7c15_9e37_u64;
        let mut word = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let keys: [RawWordStateKey; 128] = std::array::from_fn(|_| {
            key(word() | 1, word() as i32, word() as u32, [word(), word()])
        });
        let mut cache = PureCache::new();
        let mut original_fifo = std::collections::VecDeque::new();
        for operation in 0..4096 {
            let key = keys[word() as usize % keys.len()];
            let words: [u64; 12] = std::array::from_fn(|_| word());
            cache.insert(Kind::RawWord, &key.words(), &words);
            if !original_fifo.iter().any(|(old, _)| *old == key) {
                if original_fifo.len() == 64 {
                    original_fifo.pop_front();
                }
                original_fifo.push_back((key, words));
            }
            if operation % 16 == 0 {
                for &probe in &keys {
                    assert_eq!(
                        get(&cache, Kind::RawWord, &probe.words()),
                        original_fifo
                            .iter()
                            .find(|(key, _)| *key == probe)
                            .map(|(_, words)| words.to_vec())
                    );
                }
                valid(&cache);
            }
        }
        let original = key(1, -3, 7, [2451545.0_f64.to_bits(), 0]);
        put(
            &mut cache,
            Kind::RawWord,
            &original.words(),
            0x7ff8_0000_0000_1234,
        );
        assert_eq!(
            get(&cache, Kind::RawWord, &original.words()),
            Some(vec![0x7ff8_0000_0000_1234; 12])
        );
        for changed in [
            key(2, -3, 7, original.tt),
            key(1, 3, 7, original.tt),
            key(1, -3, 8, original.tt),
            key(1, -3, 7, [original.tt[0] ^ 1, 0]),
            key(1, -3, 7, [original.tt[0], 1]),
            key(1, -3, 7, [original.tt[0], (-0.0_f64).to_bits()]),
        ] {
            assert!(get(&cache, Kind::RawWord, &changed.words()).is_none());
        }
    }

    #[test]
    fn shared_pool_charges_all_blocks_descriptors_and_index_inside_the_global_cap() {
        assert_eq!(std::mem::size_of::<Descriptor>(), 16);
        assert_eq!(std::mem::size_of::<[[u64; BLOCK_WORDS]; BLOCKS]>(), 55_296);
        assert_eq!(std::mem::size_of::<[Descriptor; ENTRIES]>(), 6_144);
        assert_eq!(std::mem::size_of::<[u16; BLOCKS]>(), 864);
        assert_eq!(std::mem::size_of::<[u16; BUCKETS]>(), 2_048);
        const { assert!(RETAINED_BYTES <= 64 * 1024) };
        println!("shared numerical pool retained bytes: {RETAINED_BYTES}");
    }

    #[test]
    fn compact_word_key_keeps_identity_body_record_and_both_original_time_words() {
        let original = key(1, -3, 7, [2451545.0_f64.to_bits(), 0.0_f64.to_bits()]);
        let mut cache = PureCache::new();
        put(&mut cache, Kind::RawWord, &original.words(), 1);
        assert_eq!(
            get(&cache, Kind::RawWord, &original.words()),
            Some(vec![1; 12])
        );
        for different in [
            key(2, -3, 7, original.tt),
            key(1, 3, 7, original.tt),
            key(1, -3, 8, original.tt),
            key(1, -3, 7, [original.tt[0] ^ 1, original.tt[1]]),
            key(1, -3, 7, [original.tt[0], original.tt[1] ^ 1]),
            key(1, -3, 7, [original.tt[0], (-0.0_f64).to_bits()]),
        ] {
            assert!(get(&cache, Kind::RawWord, &different.words()).is_none());
        }
    }

    #[test]
    fn compact_word_table_copies_all_twelve_ieee_words_without_arithmetic() {
        let mut cache = PureCache::new();
        let identity = key(1, 0, 0, [0, 0]);
        let value = raw_from_words(IEEE);
        cache.insert(Kind::RawWord, &identity.words(), &raw_words(value));
        let mut words = [0; 12];
        assert!(cache.get(Kind::RawWord, &identity.words(), &mut words));
        assert_eq!(raw_words(raw_from_words(words)), IEEE);
    }

    #[test]
    fn compact_word_capacity_and_duplicate_completion_are_bounded() {
        let mut cache = PureCache::new();
        for word in 0..64 {
            put(
                &mut cache,
                Kind::RawWord,
                &key(1, 0, 0, [word, 0]).words(),
                word,
            );
        }
        let order = (cache.oldest, cache.newest, cache.free_blocks);
        put(
            &mut cache,
            Kind::RawWord,
            &key(1, 0, 0, [3, 0]).words(),
            999,
        );
        assert_eq!((cache.oldest, cache.newest, cache.free_blocks), order);
        assert_eq!(
            get(&cache, Kind::RawWord, &key(1, 0, 0, [3, 0]).words()),
            Some(vec![3; 12])
        );
        for word in 64..1000 {
            put(
                &mut cache,
                Kind::RawWord,
                &key(1, 0, 0, [word, 0]).words(),
                word,
            );
        }
        assert_eq!(cache.domain_count[Kind::RawWord.index()], 64);
        for word in 936..1000 {
            assert_eq!(
                get(&cache, Kind::RawWord, &key(1, 0, 0, [word, 0]).words()),
                Some(vec![word; 12])
            );
        }
        valid(&cache);
    }

    #[test]
    fn word_image_identity_exhaustion_never_wraps_or_reuses_a_token() {
        let counter = AtomicU64::new(u64::MAX - 1);
        assert_eq!(
            issue_word_image_identity(&counter).unwrap().get(),
            u64::MAX - 1
        );
        for _ in 0..4 {
            assert!(issue_word_image_identity(&counter).is_none());
        }
        assert_eq!(counter.load(Ordering::Relaxed), u64::MAX);
        assert!(issue_word_image_identity(&AtomicU64::new(0)).is_none());
    }

    #[test]
    fn concurrent_identity_issue_and_fixed_table_access_do_not_alias() {
        let counter = std::sync::Arc::new(AtomicU64::new(1));
        let mut workers = Vec::new();
        for _ in 0..8 {
            let counter = std::sync::Arc::clone(&counter);
            workers.push(std::thread::spawn(move || {
                (0..64)
                    .map(|_| issue_word_image_identity(&counter).unwrap().get())
                    .collect::<Vec<_>>()
            }));
        }
        let mut identities = workers
            .into_iter()
            .flat_map(|worker| worker.join().unwrap())
            .collect::<Vec<_>>();
        identities.sort_unstable();
        assert_eq!(identities, (1..=512).collect::<Vec<_>>());
        let cache = std::sync::Arc::new(Mutex::new(PureCache::new()));
        let mut workers = Vec::new();
        for lane in 0..8 {
            let cache = std::sync::Arc::clone(&cache);
            workers.push(std::thread::spawn(move || {
                for index in 0..32 {
                    let id = lane * 32 + index;
                    let identity = key(id + 1, lane as i32, index as u32, [id, 0]);
                    let mut cache = cache.lock().unwrap();
                    put(&mut cache, Kind::RawWord, &identity.words(), id);
                    assert_eq!(
                        get(&cache, Kind::RawWord, &identity.words()),
                        Some(vec![id; 12])
                    );
                }
            }));
        }
        for worker in workers {
            worker.join().unwrap();
        }
        let cache = cache.lock().unwrap();
        assert_eq!(cache.domain_count[Kind::RawWord.index()], 64);
        let mut slot = cache.domain_oldest[Kind::RawWord.index()];
        while slot != NONE {
            let entry = cache.entries[usize::from(slot)];
            let key = &cache.blocks[usize::from(entry.first)][..4];
            assert_eq!(get(&cache, Kind::RawWord, key), Some(vec![key[2]; 12]));
            slot = entry.domain_next;
        }
        valid(&cache);
    }

    #[test]
    fn shared_pool_preserves_each_original_domain_maximum() {
        for (kind, limit) in [
            (Kind::Frame, 31),
            (Kind::RawWord, 64),
            (Kind::FrameRate, 16),
            (Kind::AnalyticFrame, 16),
            (Kind::Houses, 128),
            (Kind::ReferenceNutation, 128),
            (Kind::ErfaNutation, 128),
        ] {
            let mut cache = PureCache::new();
            assert_eq!(kind.limit(), limit);
            for word in 0..u64::from(limit) {
                put(&mut cache, kind, &vec![word; kind.key_words()], word);
            }
            assert_eq!(cache.domain_count[kind.index()], limit);
            put(&mut cache, kind, &vec![99_999; kind.key_words()], 99_999);
            assert_eq!(cache.domain_count[kind.index()], limit);
            assert!(get(&cache, kind, &vec![0; kind.key_words()]).is_none());
            assert_eq!(
                get(&cache, kind, &vec![99_999; kind.key_words()]),
                Some(vec![99_999; kind.payload_words()])
            );
            valid(&cache);
        }
    }
}

#[cfg(test)]
mod owned_apparent_cache_tests {
    use super::test_support::*;
    use super::*;

    #[test]
    fn owned_producer_ids_never_wrap_or_publish_zero() {
        let counter = AtomicU64::new(1);
        assert_eq!(
            owned_apparent_identity(&counter).map(|id| id.get()),
            Some(1)
        );
        assert_eq!(
            owned_apparent_identity(&counter).map(|id| id.get()),
            Some(2)
        );
        let exhausted = AtomicU64::new(u64::MAX);
        assert_eq!(owned_apparent_identity(&exhausted), None);
        assert_eq!(owned_apparent_identity(&exhausted), None);
        assert_eq!(
            exhausted.load(std::sync::atomic::Ordering::Relaxed),
            u64::MAX
        );
    }

    #[test]
    fn owned_entry_preserves_all_key_and_payload_words() {
        let mut cache = PureCache::new();
        let key: [u64; 9] = std::array::from_fn(|i| 0x8000_0000_0000_0000 | i as u64);
        let words: [u64; 15] = std::array::from_fn(|i| IEEE[i % IEEE.len()]);
        cache.insert(Kind::OwnedApparent, &key, &words);
        assert_eq!(get(&cache, Kind::OwnedApparent, &key), Some(words.to_vec()));
        for index in 0..key.len() {
            let mut other = key;
            other[index] ^= 1;
            assert!(get(&cache, Kind::OwnedApparent, &other).is_none());
        }
        cache.insert(Kind::OwnedApparent, &[0; 9], &[1; 15]);
        // Capacity is deliberately shared now, rather than a one-entry contract.
        assert_eq!(get(&cache, Kind::OwnedApparent, &key), Some(words.to_vec()));
    }

    #[test]
    fn owned_entries_and_counter_stay_inside_original_combined_cap() {
        let mut cache = PureCache::new();
        for word in 0..100 {
            put(&mut cache, Kind::OwnedApparent, &[word; 9], word);
        }
        assert_eq!(cache.domain_count[Kind::OwnedApparent.index()], 100);
        for word in 0..100 {
            assert_eq!(
                get(&cache, Kind::OwnedApparent, &[word; 9]),
                Some(vec![word; 15])
            );
        }
        const { assert!(RETAINED_BYTES <= 64 * 1024) };
        valid(&cache);
    }
}

#[cfg(test)]
mod shared_pool_tests {
    use super::test_support::*;
    use super::*;

    #[test]
    fn linked_blocks_reuse_fragmented_space_and_global_fifo_without_extra_eviction() {
        let mut cache = PureCache::new();
        for word in 0..64 {
            put(&mut cache, Kind::RawWord, &[1, 0, word, 0], word);
        }
        for word in 0..31 {
            put(&mut cache, Kind::Frame, &[word; 3], word);
        }
        for word in 0..100 {
            put(&mut cache, Kind::OwnedApparent, &[word; 9], word);
        }
        for word in (0..64).step_by(2) {
            let key = [1, 0, word, 0];
            let bucket = cache
                .find_bucket(Kind::RawWord, &key, PureCache::hash(Kind::RawWord, &key))
                .unwrap();
            cache.evict(cache.buckets[bucket] - 1);
        }
        for word in 0..16 {
            put(&mut cache, Kind::AnalyticFrame, &[word; 6], word);
        }
        for word in 0..14 {
            put(&mut cache, Kind::FrameRate, &[word; 6], word);
        }
        assert_eq!(cache.free_blocks, 0);
        for word in (1..64).step_by(2) {
            assert_eq!(
                get(&cache, Kind::RawWord, &[1, 0, word, 0]),
                Some(vec![word; 12])
            );
        }
        for word in 0..100 {
            assert!(get(&cache, Kind::OwnedApparent, &[word; 9]).is_some());
        }
        valid(&cache);
        put(&mut cache, Kind::OwnedApparent, &[100; 9], 100);
        // Exactly the two oldest one-block entries fund this two-block record.
        assert!(get(&cache, Kind::RawWord, &[1, 0, 1, 0]).is_none());
        assert!(get(&cache, Kind::RawWord, &[1, 0, 3, 0]).is_none());
        assert!(get(&cache, Kind::RawWord, &[1, 0, 5, 0]).is_some());
        assert_eq!(cache.free_blocks, 0);
        valid(&cache);
    }

    #[test]
    fn domain_limits_prevent_intermediate_raw_churn_from_flooding_spare_capacity() {
        let mut cache = PureCache::new();
        // Synthetic cardinalities demonstrate storage, not measured native hits.
        for word in 0..100 {
            put(&mut cache, Kind::OwnedApparent, &[word; 9], word);
        }
        for word in 0..10 {
            put(&mut cache, Kind::Frame, &[word; 3], word);
            put(&mut cache, Kind::AnalyticFrame, &[word; 6], word);
            put(&mut cache, Kind::ReferenceNutation, &[word], word);
            put(&mut cache, Kind::ErfaNutation, &[word], word);
        }
        for word in 0..4096 {
            put(&mut cache, Kind::RawWord, &[1, 0, word, 0], word);
        }
        assert_eq!(cache.free_blocks, (BLOCKS - 354) as u16);
        assert_eq!(cache.len, 204);
        assert_eq!(cache.domain_count[Kind::RawWord.index()], 64);
        for word in 0..100 {
            assert_eq!(
                get(&cache, Kind::OwnedApparent, &[word; 9]),
                Some(vec![word; 15])
            );
        }
        for word in 0..10 {
            assert!(get(&cache, Kind::Frame, &[word; 3]).is_some());
            assert!(get(&cache, Kind::AnalyticFrame, &[word; 6]).is_some());
        }
        valid(&cache);
    }

    #[test]
    fn typed_source_frame_and_optional_native_fields_roundtrip_all_active_ieee_bits() {
        let mut cache = PureCache::new();
        let words: [u64; 12] = std::array::from_fn(|i| IEEE[i]);
        assert_eq!(raw_words(raw_from_words(words)), words);
        for some in [false, true] {
            let mut words: [u64; 49] = std::array::from_fn(|i| IEEE[i % IEEE.len()]);
            words[38] = u64::from(some);
            if !some {
                words[39..].fill(0);
            }
            let frame = FrameMatrices::from_cache_words(words);
            assert_eq!(frame.cache_words(), words);
            cache.insert(Kind::Frame, &[u64::from(some), 0, 0], &frame.cache_words());
            let mut restored = [0; 49];
            assert!(cache.get(Kind::Frame, &[u64::from(some), 0, 0], &mut restored));
            assert_eq!(
                FrameMatrices::from_cache_words(restored).cache_words(),
                words
            );
        }
        let analytic: [u64; 36] = std::array::from_fn(|i| IEEE[i % IEEE.len()]);
        assert_eq!(analytic_words(analytic_from_words(analytic)), analytic);
        let rate: [u64; 18] = std::array::from_fn(|i| IEEE[i % IEEE.len()]);
        assert_eq!(frame_rate_words(frame_rate_from_words(rate)), rate);
        for kind in [Kind::Houses, Kind::ReferenceNutation, Kind::ErfaNutation] {
            let payload: Vec<_> = (0..kind.payload_words())
                .map(|i| IEEE[i % IEEE.len()])
                .collect();
            let key = vec![0x8000_0000_0000_0000; kind.key_words()];
            cache.insert(kind, &key, &payload);
            assert_eq!(get(&cache, kind, &key), Some(payload));
        }
        valid(&cache);
    }

    #[test]
    fn cross_domain_global_pressure_duplicate_completions_and_repeated_wrap_preserve_words() {
        let mut cache = PureCache::new();
        for word in 0..2000_u64 {
            let kind = match word % 4 {
                0 => Kind::OwnedApparent,
                1 => Kind::Houses,
                2 => Kind::FrameRate,
                _ => Kind::AnalyticFrame,
            };
            let key = vec![word; kind.key_words()];
            let words: Vec<_> = (0..kind.payload_words())
                .map(|i| IEEE[i % IEEE.len()] ^ word)
                .collect();
            cache.insert(kind, &key, &words);
            assert_eq!(get(&cache, kind, &key), Some(words.clone()));
            let state = (
                cache.oldest,
                cache.newest,
                cache.free_blocks,
                cache.len,
                cache.domain_count,
            );
            cache.insert(kind, &key, &vec![0; kind.payload_words()]);
            assert_eq!(
                (
                    cache.oldest,
                    cache.newest,
                    cache.free_blocks,
                    cache.len,
                    cache.domain_count
                ),
                state
            );
            assert_eq!(get(&cache, kind, &key), Some(words));
            if word % 17 == 0 {
                valid(&cache);
            }
        }
        valid(&cache);
    }

    #[test]
    fn concurrent_cross_domain_completions_keep_the_first_complete_words() {
        let cache = std::sync::Arc::new(Mutex::new(PureCache::new()));
        let start = std::sync::Arc::new(std::sync::Barrier::new(8));
        let mut workers = Vec::new();
        for &payload_word in IEEE.iter().take(8) {
            let cache = std::sync::Arc::clone(&cache);
            let start = std::sync::Arc::clone(&start);
            workers.push(std::thread::spawn(move || {
                start.wait();
                let mut observed = Vec::new();
                for word in 0..4 {
                    for kind in [Kind::OwnedApparent, Kind::Houses, Kind::FrameRate] {
                        // Each producer finishes its pure payload before locking.
                        let key = vec![word; kind.key_words()];
                        let payload = vec![payload_word; kind.payload_words()];
                        let mut cache = cache.lock().unwrap();
                        cache.insert(kind, &key, &payload);
                        let kept = get(&cache, kind, &key).unwrap();
                        assert!(kept.iter().all(|&value| value == kept[0]));
                        observed.push((kind, key, kept));
                    }
                }
                observed
            }));
        }
        let observed: Vec<_> = workers
            .into_iter()
            .flat_map(|worker| worker.join().unwrap())
            .collect();
        let cache = cache.lock().unwrap();
        for (kind, key, payload) in observed {
            assert_eq!(get(&cache, kind, &key), Some(payload));
        }
        assert_eq!(cache.len, 12);
        assert_eq!(cache.free_blocks, (BLOCKS - 24) as u16);
        valid(&cache);
    }
}
