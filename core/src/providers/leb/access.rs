//! Bounded range access shared by contiguous and positional LEB images.
//!
//! This layer owns byte acquisition only. Segment ownership, polynomial
//! evaluation and astronomical units remain in the existing reader kernels.

use super::*;
use std::{
    borrow::Cow,
    fs::{File, Metadata},
    io,
    path::Path,
};

/// Allocation limits for a shared positional LEB domain, in bytes except the
/// window log. Limits constrain metadata, retained/live coefficient capacities
/// and decode reservations, not total RSS. OS cache, allocator/decoder overhead,
/// fixed-size header reads and caller-owned return values are additional.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IndexedLimits {
    /// Shared parsed metadata and cached directory capacity, excluding map and
    /// allocator bookkeeping. Temporary directory bytes count while parsing.
    pub metadata_bytes: usize,
    /// Shared retained coefficient budget; zero disables cache retention.
    pub cache_bytes: usize,
    /// Shared live coefficient/decode reservation budget, including advertised
    /// decoder-window capacity. OS/allocator overhead and return vectors are extra.
    pub working_bytes: usize,
    /// Maximum single compressed range materialized before decode.
    pub compressed_bytes: usize,
    /// Maximum decoded coefficient frame, including legacy monolithic frames.
    pub frame_bytes: usize,
    /// Largest zstd window admitted, expressed as log2(bytes).
    pub window_log_max: u32,
}

impl Default for IndexedLimits {
    fn default() -> Self {
        Self {
            metadata_bytes: 8 * 1024 * 1024,
            cache_bytes: 8 * 1024 * 1024,
            working_bytes: 64 * 1024 * 1024,
            compressed_bytes: 4 * 1024 * 1024,
            frame_bytes: 8 * 1024 * 1024,
            window_log_max: 23,
        }
    }
}

/// Storage policy for future LEB file acquisitions. Existing explicitly
/// attached images and already acquired Session snapshots are never reopened.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum LebReadPolicy {
    /// Whole-file immutable snapshot; safe for arbitrary regular input files.
    #[default]
    Owned,
    /// Positional reads under the immutable-inode and resource contracts of
    /// [`LebFile::open_indexed`]. No automatic owned-buffer fallback is allowed.
    Indexed(IndexedLimits),
}

impl LebReadPolicy {
    /// Open an image with this policy; units and errors are those of the
    /// corresponding `LebFile` constructor. No discovery or network is used.
    pub fn open(self, path: &Path) -> Result<LebFile> {
        match self {
            Self::Owned => LebFile::open(path),
            Self::Indexed(limits) => LebFile::open_indexed(path, limits),
        }
    }

    /// Validate resource settings without opening a file or changing state.
    /// Zero byte limits are permitted and refuse work at acquisition/decode.
    pub fn validate(self) -> Result<()> {
        if let Self::Indexed(limits) = self {
            if !(10..=30).contains(&limits.window_log_max) {
                return Err(Error::InvalidLebReader {
                    value: "window_log_max must be in 10..=30".into(),
                });
            }
        }
        Ok(())
    }
}

/// Descriptor identity is retained across atomic pathname replacement.
/// Size/mtime checks diagnose some in-place writes, but do not prove immutability.
pub(crate) struct PositionalFile {
    file: File,
    identity: Metadata,
    path: String,
    pub(crate) limits: IndexedLimits,
    pub(super) digest_gate: std::sync::Mutex<()>,
    pub(super) main_index_bytes: usize,
    pub(crate) resources: Arc<IndexedResources>,
    _index_allocation: Option<resources::Reservation>,
    pub(super) cache: Arc<LebChunkCache>,
    nutation_header: std::sync::OnceLock<NutationHeader>,
    delta_t_header: std::sync::OnceLock<DeltaTHeader>,
}

impl PositionalFile {
    pub(super) fn open(path: &Path, resources: Arc<IndexedResources>) -> Result<Self> {
        let limits = resources.limits;
        LebReadPolicy::Indexed(limits).validate()?;
        let file = super::open_regular_file(path).map_err(|error| {
            if error.kind() == io::ErrorKind::NotFound {
                Error::LebFileNotFound {
                    path: path.display().to_string(),
                }
            } else {
                io_error(path.to_string_lossy().as_ref(), error)
            }
        })?;
        let identity = file
            .metadata()
            .map_err(|error| io_error(&path.display().to_string(), error))?;
        if identity.len() > isize::MAX as u64 {
            return Err(corrupt("file length", usize::MAX));
        }
        let cache = Arc::new(LebChunkCache::with_indexed_resources(Arc::clone(
            &resources,
        )));
        Ok(Self {
            file,
            identity,
            path: path.display().to_string(),
            limits,
            digest_gate: std::sync::Mutex::new(()),
            main_index_bytes: 0,
            resources,
            _index_allocation: None,
            cache,
            nutation_header: std::sync::OnceLock::new(),
            delta_t_header: std::sync::OnceLock::new(),
        })
    }

    pub(super) fn len(&self) -> usize {
        self.identity.len() as usize
    }

    /// Identity of the descriptor admitted by the indexed reader. Offline
    /// writers use it to protect the acquired source inode from output aliases.
    pub(super) fn acquired_metadata(&self) -> &Metadata {
        &self.identity
    }

    pub(super) fn validate_identity(&self) -> Result<()> {
        let now = self
            .file
            .metadata()
            .map_err(|error| io_error(&self.path, error))?;
        if now.len() != self.identity.len() || now.modified().ok() != self.identity.modified().ok()
        {
            return Err(corrupt("immutable file changed", 0));
        }
        Ok(())
    }

    pub(super) fn read_exact(&self, offset: usize, mut output: &mut [u8]) -> Result<()> {
        if offset
            .checked_add(output.len())
            .is_none_or(|end| end > self.len())
        {
            return Err(corrupt("positional range", offset));
        }
        self.validate_identity()?;
        let mut at = offset as u64;
        while !output.is_empty() {
            #[cfg(unix)]
            let result = {
                use std::os::unix::fs::FileExt;
                self.file.read_at(output, at)
            };
            #[cfg(windows)]
            let result = {
                use std::os::windows::fs::FileExt;
                self.file.seek_read(output, at)
            };
            #[cfg(not(any(unix, windows)))]
            let result: io::Result<usize> = Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "positional LEB reads unavailable",
            ));
            match result {
                Ok(0) => return Err(corrupt("truncated positional read", at as usize)),
                Ok(count) => {
                    at += count as u64;
                    output = &mut output[count..];
                }
                Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
                Err(error) => return Err(io_error(&self.path, error)),
            }
        }
        self.validate_identity()
    }

    pub(super) fn read_index(&mut self) -> Result<FileIndex> {
        let mut header_bytes = [0; HEADER_SIZE];
        self.read_exact(0, &mut header_bytes)?;
        let header = FileHeader::decode(&header_bytes)?;
        // Bound decoded capacities before the existing index reader allocates.
        let bytes = (header.body_count() as usize)
            .checked_mul(std::mem::size_of::<BodyEntry>())
            .and_then(|count| {
                (header.section_count() as usize)
                    .checked_mul(std::mem::size_of::<SectionEntry>())
                    .and_then(|sections| count.checked_add(sections))
            })
            .and_then(|count| count.checked_mul(2))
            .ok_or_else(|| corrupt("main index allocation", 0))?;
        limit("LEB main index", bytes, self.limits.metadata_bytes)?;
        let allocation = self.resources.metadata.reserve(bytes)?;
        let mut descriptor = self
            .file
            .try_clone()
            .map_err(|error| io_error(&self.path, error))?;
        // The descriptor is not published until this cursor-based parse ends.
        let index = FileIndex::read_from_file(&mut descriptor)?;
        self.main_index_bytes = index.bodies.capacity() * std::mem::size_of::<BodyEntry>()
            + index.sections.capacity() * std::mem::size_of::<SectionEntry>();
        limit(
            "LEB main index",
            self.main_index_bytes,
            self.limits.metadata_bytes,
        )?;
        self._index_allocation = Some(allocation.retain(self.main_index_bytes)?);
        self.validate_identity()?;
        Ok(index)
    }

    pub(super) fn sha256(&self) -> Result<String> {
        use sha2::{Digest, Sha256};
        let mut hasher = Sha256::new();
        let mut buffer = vec![0; 64 * 1024];
        let mut offset = 0;
        while offset < self.len() {
            let count = buffer.len().min(self.len() - offset);
            self.read_exact(offset, &mut buffer[..count])?;
            hasher.update(&buffer[..count]);
            offset += count;
        }
        self.validate_identity()?;
        Ok(crate::download::hex_of(&hasher.finalize()))
    }
}

/// A copyable borrowed view; owned images retain their zero-copy slice path.
#[derive(Clone, Copy)]
pub(crate) enum LebSource<'a> {
    Slice(&'a [u8]),
    Positional(&'a PositionalFile),
}

impl<'a> LebSource<'a> {
    pub(crate) fn len(self) -> usize {
        match self {
            Self::Slice(bytes) => bytes.len(),
            Self::Positional(file) => file.len(),
        }
    }

    pub(super) fn range(self, offset: usize, count: usize, what: &str) -> Result<Cow<'a, [u8]>> {
        match self {
            Self::Slice(bytes) => Ok(Cow::Borrowed(reader::slice_at(bytes, offset, count, what)?)),
            Self::Positional(file) => {
                limit("LEB positional range", count, file.limits.compressed_bytes)?;
                if offset.checked_add(count).is_none_or(|end| end > file.len()) {
                    return Err(corrupt(what, offset));
                }
                let mut bytes = vec![0; count];
                file.read_exact(offset, &mut bytes)?;
                Ok(Cow::Owned(bytes))
            }
        }
    }

    pub(super) fn nutation_header(self, base: usize) -> Result<NutationHeader> {
        if let Self::Positional(file) = self {
            if let Some(header) = file.nutation_header.get() {
                return Ok(header.clone());
            }
        }
        let data = self.range(base, NUTATION_HEADER_SIZE, "LEB nutation header")?;
        let header = NutationHeader::decode(&data, 0)?;
        if let Self::Positional(file) = self {
            let _ = file.nutation_header.set(header.clone());
        }
        Ok(header)
    }

    pub(super) fn delta_t_header(self, base: usize) -> Result<DeltaTHeader> {
        if let Self::Positional(file) = self {
            if let Some(header) = file.delta_t_header.get() {
                return Ok(header.clone());
            }
        }
        let data = self.range(base, DELTA_T_HEADER_SIZE, "LEB Delta-T header")?;
        let header = DeltaTHeader::decode(&data, 0)?;
        if let Self::Positional(file) = self {
            let _ = file.delta_t_header.set(header.clone());
        }
        Ok(header)
    }

    pub(super) fn check_frame(self, bytes: u64) -> Result<()> {
        if let Self::Positional(file) = self {
            let bytes = usize::try_from(bytes).map_err(|_| corrupt("decoded frame size", 0))?;
            limit("LEB decoded frame", bytes, file.limits.frame_bytes)?;
        }
        Ok(())
    }

    /// Reserve scratch before reading compressed bytes. Four frame capacities
    /// cover decoded/shuffled/reordered/word buffers; the gate permits one
    /// decoder window per shared store. Only final word capacity remains leased.
    pub(super) fn reserve_decode(
        self,
        frame: u64,
        compressed: u64,
    ) -> Result<Option<resources::Reservation>> {
        let Self::Positional(file) = self else {
            return Ok(None);
        };
        let bytes = frame
            .checked_mul(4)
            .and_then(|count| count.checked_add(compressed))
            .and_then(|count| count.checked_add(1u64 << file.limits.window_log_max))
            .and_then(|count| count.checked_add(1))
            .and_then(|count| usize::try_from(count).ok())
            .ok_or_else(|| corrupt("decode reservation", 0))?;
        Ok(Some(file.resources.working.reserve(bytes)?))
    }

    pub(super) fn reserve_work(self, bytes: usize) -> Result<Option<resources::Reservation>> {
        match self {
            Self::Slice(_) => Ok(None),
            Self::Positional(file) => Ok(Some(file.resources.working.reserve(bytes)?)),
        }
    }

    pub(super) fn copy_words(self, words: &[f64]) -> Result<CoefficientWords> {
        let bytes = words
            .len()
            .checked_mul(std::mem::size_of::<f64>())
            .ok_or_else(|| corrupt("coefficient row size", 0))?;
        let allocation = self.reserve_work(bytes)?;
        Ok(CoefficientWords::Local(TrackedVec {
            values: words.to_vec(),
            _allocation: allocation,
        }))
    }

    pub(super) fn read_words(
        self,
        offset: usize,
        count: usize,
        what: &str,
    ) -> Result<CoefficientWords> {
        self.check_frame(count as u64)?;
        match self {
            Self::Slice(_) => self
                .read_words_uncached(offset, count, what)
                .map(CoefficientWords::Local),
            Self::Positional(file) => file
                .cache
                .get_or_decode_key(
                    index::CacheKey::Table {
                        offset,
                        bytes: count,
                    },
                    || self.read_words_uncached(offset, count, what),
                )
                .map(CoefficientWords::Shared),
        }
    }

    /// Consume a short contiguous row without allocating coefficient storage.
    /// Indexed reads retain their ordinary shared cache and allocation leases.
    /// The stack capacity is only a storage choice: larger rows use the same
    /// checked heap reader, with no restriction on the file's polynomial degree.
    pub(super) fn with_words<R>(
        self,
        offset: usize,
        count: usize,
        what: &str,
        consume: impl FnOnce(&[f64]) -> Result<R>,
    ) -> Result<R> {
        const INLINE_WORDS: usize = 64;
        if matches!(self, Self::Slice(_)) && count <= INLINE_WORDS * 8 && count % 8 == 0 {
            let data = self.range(offset, count, what)?;
            let mut words = [0.0; INLINE_WORDS];
            for (word, encoded) in words.iter_mut().zip(data.chunks_exact(8)) {
                *word = f64::from_le_bytes(encoded.try_into().expect("eight-byte word"));
            }
            return consume(&words[..count / 8]);
        }
        let words = self.read_words(offset, count, what)?;
        consume(&words)
    }

    fn read_words_uncached(
        self,
        offset: usize,
        count: usize,
        what: &str,
    ) -> Result<TrackedVec<f64>> {
        let temporary = count.checked_mul(2).ok_or_else(|| corrupt(what, offset))?;
        let allocation = self.reserve_work(temporary)?;
        let data = self.range(offset, count, what)?;
        let values = reader::f64_words(&data, what, offset)?;
        drop(data);
        let capacity = values.capacity() * std::mem::size_of::<f64>();
        let allocation = allocation.map(|lease| lease.retain(capacity)).transpose()?;
        Ok(TrackedVec {
            values,
            _allocation: allocation,
        })
    }

    pub(super) fn window_limit(self) -> Option<u32> {
        match self {
            Self::Slice(_) => None,
            Self::Positional(file) => Some(file.limits.window_log_max),
        }
    }

    pub(super) fn blob(self, offset: u64, count: u64, what: &str) -> Result<Cow<'a, [u8]>> {
        let offset = usize::try_from(offset).map_err(|_| corrupt(what, usize::MAX))?;
        let count = usize::try_from(count).map_err(|_| corrupt(what, offset))?;
        self.range(offset, count, what)
    }
}

pub(super) fn limit(resource: &'static str, required: usize, maximum: usize) -> Result<()> {
    if required > maximum {
        Err(Error::ResourceLimit {
            resource,
            required,
            limit: maximum,
        })
    } else {
        Ok(())
    }
}

/// Inspect the advertised window before the decoder can allocate it. Invalid
/// frame headers remain the ordinary decoder's responsibility.
pub(super) fn check_zstd_window(compressed: &[u8], log_max: u32) -> Result<()> {
    if compressed.len() < 6 || compressed[..4] != [0x28, 0xb5, 0x2f, 0xfd] {
        return Ok(());
    }
    let window = if compressed[4] & 0x20 == 0 {
        let descriptor = compressed[5];
        let base = 1u64 << (10 + (descriptor >> 3));
        base + (base / 8) * u64::from(descriptor & 7)
    } else {
        match zstd::zstd_safe::get_frame_content_size(compressed) {
            Ok(Some(size)) => size,
            _ => return Ok(()),
        }
    };
    let maximum = 1usize
        .checked_shl(log_max)
        .ok_or_else(|| corrupt("zstd window limit", 0))?;
    limit(
        "LEB zstd window",
        usize::try_from(window).unwrap_or(usize::MAX),
        maximum,
    )
}

fn io_error(path: &str, error: io::Error) -> Error {
    Error::EphemerisLoaderIo {
        path: path.into(),
        kind: error.kind(),
        errno: error.raw_os_error(),
        message: error.to_string(),
    }
}

fn corrupt(what: &str, offset: usize) -> Error {
    Error::LebCorrupt {
        what: what.into(),
        offset,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scoped_contiguous_words_preserve_bits_across_storage_capacities() {
        let patterns = [
            0,
            1 << 63,
            1,
            0x3ff0_0000_0000_0000,
            0x7ff0_0000_0000_0000,
            0x7ff0_0000_0000_0001,
            0xfff8_0123_4567_89ab,
            u64::MAX,
        ];
        for count in [0, 1, 34, 64, 65, 4096] {
            let expected: Vec<u64> = patterns.into_iter().cycle().take(count).collect();
            // An unaligned start must work without a native f64 pointer cast.
            let mut encoded = vec![0xff];
            encoded.extend(expected.iter().flat_map(|word| word.to_le_bytes()));
            let source = LebSource::Slice(&encoded);
            source
                .with_words(1, count * 8, "test words", |words| {
                    assert_eq!(
                        words.iter().map(|word| word.to_bits()).collect::<Vec<_>>(),
                        expected
                    );
                    Ok(())
                })
                .unwrap();
        }
    }

    #[test]
    fn scoped_contiguous_words_preserve_range_and_word_errors() {
        let encoded = [0; 24];
        let source = LebSource::Slice(&encoded);
        for (offset, count) in [(0, 7), (23, 8), (24, 1), (usize::MAX, 8)] {
            let original = source
                .read_words(offset, count, "test words")
                .err()
                .unwrap();
            let scoped = source
                .with_words(offset, count, "test words", |_| -> Result<()> {
                    panic!("invalid input reached the coefficient consumer")
                })
                .unwrap_err();
            assert_eq!(format!("{original:?}"), format!("{scoped:?}"));
        }
        let error = source
            .with_words(0, 8, "test words", |_| -> Result<()> {
                Err(Error::NonFiniteArgument {
                    value: f64::INFINITY,
                })
            })
            .unwrap_err();
        assert!(matches!(error, Error::NonFiniteArgument { value } if value == f64::INFINITY));
    }

    fn fixture() -> std::path::PathBuf {
        crate::test_data::core_leb_path()
    }

    #[test]
    fn indexed_body_and_auxiliary_words_match_owned_image() {
        let owned = LebFile::open(&fixture()).unwrap();
        let indexed = LebFile::open_indexed(&fixture(), IndexedLimits::default()).unwrap();
        assert_eq!(indexed.allocated_image_bytes(), 0);
        assert_eq!(indexed.mapped_file_bytes(), 0);
        assert_eq!(indexed.len(), owned.len());
        assert!(matches!(indexed.bytes(), Err(Error::LebBytesUnavailable)));
        assert!(matches!(
            indexed.read_range(0, 8),
            Err(Error::LebBytesUnavailable)
        ));
        for body in owned.index().unwrap().bodies() {
            for fraction in [0.0, 0.001, 0.25, 0.5, 0.999, 1.0] {
                let jd = body.jd_start() + fraction * (body.jd_end() - body.jd_start());
                let words = |value: (Vec<f64>, Vec<f64>)| {
                    (
                        value.0.into_iter().map(f64::to_bits).collect::<Vec<_>>(),
                        value.1.into_iter().map(f64::to_bits).collect::<Vec<_>>(),
                    )
                };
                assert_eq!(
                    words(owned.eval_body_any(body.body_id(), jd).unwrap()),
                    words(indexed.eval_body_any(body.body_id(), jd).unwrap())
                );
                let epoch = TwoPartEpoch { hi: jd, lo: 1e-12 };
                assert_eq!(
                    words(owned.eval_body_any_two_part(body.body_id(), epoch).unwrap()),
                    words(
                        indexed
                            .eval_body_any_two_part(body.body_id(), epoch)
                            .unwrap()
                    )
                );
            }
        }
        for jd in [2400000.0, 2451545.0, 2500000.0] {
            assert_eq!(
                owned.eval_nutation(jd).unwrap(),
                indexed.eval_nutation(jd).unwrap()
            );
            assert_eq!(
                owned.delta_t(jd).unwrap().to_bits(),
                indexed.delta_t(jd).unwrap().to_bits()
            );
            assert_eq!(
                owned
                    .index()
                    .unwrap()
                    .nutation_rate(owned.bytes().unwrap(), jd)
                    .unwrap(),
                indexed
                    .index()
                    .unwrap()
                    .nutation_rate_source(indexed.source().unwrap(), jd)
                    .unwrap()
            );
        }
        assert_eq!(owned.sha256().unwrap(), indexed.sha256().unwrap());
        assert!(indexed == indexed.clone());
        assert!(owned != indexed);
        assert!(LebFile::open_indexed(&fixture(), IndexedLimits::default()).unwrap() != indexed);
    }

    #[test]
    fn zstd_window_admission_handles_descriptor_and_single_segment_frames() {
        let descriptor = [0x28, 0xb5, 0x2f, 0xfd, 0, 1]; // 1,152-byte window
        assert!(matches!(
            check_zstd_window(&descriptor, 10),
            Err(Error::ResourceLimit {
                resource: "LEB zstd window",
                required: 1152,
                limit: 1024
            })
        ));
        check_zstd_window(&descriptor, 11).unwrap();
        let mut single = vec![0x28, 0xb5, 0x2f, 0xfd, 0xe0];
        single.extend_from_slice(&4096u64.to_le_bytes());
        assert!(matches!(
            check_zstd_window(&single, 11),
            Err(Error::ResourceLimit {
                required: 4096,
                limit: 2048,
                ..
            })
        ));
        check_zstd_window(&single, 12).unwrap();
        check_zstd_window(b"bad", 10).unwrap(); // malformed frames belong to zstd
    }

    #[test]
    fn positional_limits_refuse_main_index_and_frames_before_decode() {
        assert!(matches!(
            LebFile::open_indexed(
                &fixture(),
                IndexedLimits {
                    metadata_bytes: 0,
                    ..IndexedLimits::default()
                }
            ),
            Err(Error::ResourceLimit {
                resource: "LEB main index",
                ..
            })
        ));
        let indexed = LebFile::open_indexed(
            &fixture(),
            IndexedLimits {
                frame_bytes: 0,
                ..IndexedLimits::default()
            },
        )
        .unwrap();
        assert!(matches!(
            indexed.eval_body_any(0, 2451545.0),
            Err(Error::ResourceLimit {
                resource: "LEB decoded frame",
                ..
            })
        ));
    }
}
