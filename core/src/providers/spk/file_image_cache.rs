//! Bounded weak hints for immutable images admitted from real files.
//!
//! A hint never admits a file: registration still opens its pathname, compares
//! every byte through EOF, and verifies the acquired descriptor's identity on
//! both sides of that comparison. Only bytes, the summary index and lazy
//! coefficient storage can be shared. No file, cursor, body map, configuration
//! or foreign-runtime object is retained by this process-wide bookkeeping.

use super::{DafEndian, DecodedSpkCache, SpkIndex, SpkKernel, SpkRegistry};
use std::collections::VecDeque;
use std::fs::File;
use std::io::{self, Read};
use std::sync::{Arc, Mutex, OnceLock, Weak};

/// Maximum number of weak hints, including hints whose images are still live.
const MAX_FILE_IMAGE_HINTS: usize = 128;

/// Acquired regular-file metadata, used only to locate a comparison candidate.
/// Exact path spelling is the registry identity; device/inode, byte length and
/// nanosecond modification/change clocks describe the opened descriptor.
/// Coarse, restored or colliding metadata never establishes content equality.
#[derive(Clone, PartialEq, Eq)]
pub(super) struct FileImageKey {
    path: String,
    device: u64,
    inode: u64,
    length: u64,
    modified: (i64, i64),
    changed: (i64, i64),
}

impl FileImageKey {
    /// Read the descriptor without resolving its pathname. Metadata failures,
    /// nonregular files and unsupported hosts disable hints rather than adding
    /// a refusal to the ordinary open/read/decode contract.
    pub(super) fn acquired(path: &str, file: &File) -> Option<Self> {
        #[cfg(unix)]
        {
            use std::os::unix::fs::MetadataExt;
            let metadata = file.metadata().ok()?;
            if !metadata.is_file() {
                return None;
            }
            Some(Self {
                path: path.to_owned(),
                device: metadata.dev(),
                inode: metadata.ino(),
                length: metadata.len(),
                modified: (metadata.mtime(), metadata.mtime_nsec()),
                changed: (metadata.ctime(), metadata.ctime_nsec()),
            })
        }
        #[cfg(not(unix))]
        {
            let _ = (path, file);
            None
        }
    }
}

/// Small lifetime carrier owned strongly only by actual file-backed kernels.
/// The cache holds a weak carrier, rather than a direct weak unsized byte Arc:
/// its weak byte/index fields are dropped with the last kernel, allowing the
/// large byte allocation itself to retire even before the hint is pruned.
pub(super) struct FileImageOwner {
    bytes: Weak<[u8]>,
    index: Weak<SpkIndex>,
    endian: DafEndian,
    // Fixed metadata on this existing owner; no bytes or descriptor is retained.
    sha256: OnceLock<[u8; 32]>,
}

impl FileImageOwner {
    /// Describe existing immutable storage without adding strong image owners.
    pub(super) fn acquired(kernel: &SpkKernel) -> Self {
        Self {
            bytes: Arc::downgrade(&kernel.bytes),
            index: Arc::downgrade(&kernel.index),
            endian: kernel.endian,
            sha256: OnceLock::new(),
        }
    }

    /// Hash only this owner's immutable image, after its ordinary admission.
    /// The owning kernel supplies the same bytes described by this carrier.
    /// Publication shares the fixed digest with clones and validated file hints;
    /// it never admits a pathname or replaces fresh open/full-byte comparison.
    pub(super) fn sha256(&self, bytes: &[u8]) -> &[u8; 32] {
        use sha2::{Digest, Sha256};
        self.sha256.get_or_init(|| Sha256::digest(bytes).into())
    }
}

/// Weak lifetime carrier and independently weak lazy coefficient storage.
/// Coefficient cells can expire before bytes held by a standalone kernel or
/// captured range; such an owner permits byte/index reuse with fresh cold cells.
struct WeakFileImage {
    owner: Weak<FileImageOwner>,
    decoded: Weak<DecodedSpkCache>,
}

impl WeakFileImage {
    /// True only while real native owners retain both image and summary index.
    fn live(&self) -> bool {
        self.owner.strong_count() != 0
    }

    /// Upgrade immutable storage into a temporary map-free source registry.
    /// It owns no acquisition descriptor. The destination receives its newly
    /// opened file only through `retain_kernel_from_acquired_file`.
    fn retained(&self, key: &str) -> Option<SpkRegistry> {
        let owner = self.owner.upgrade()?;
        let bytes = owner.bytes.upgrade()?;
        let index = owner.index.upgrade()?;
        let decoded = self.decoded.upgrade().unwrap_or_else(|| {
            Arc::new(DecodedSpkCache {
                entries: (0..index.len()).map(|_| OnceLock::new()).collect(),
            })
        });
        let mut source = SpkRegistry::new();
        source.kernels.insert(
            key.to_owned(),
            SpkKernel {
                path: key.to_owned(),
                bytes,
                endian: owner.endian,
                index,
                source_file: None,
                file_image_owner: Some(owner),
            },
        );
        source.decoded.insert(key.to_owned(), decoded);
        Some(source)
    }
}

/// Least-recently-used weak bookkeeping, bounded independently of owner lives.
/// Locks protect only this deque; no open, read, seek or decode occurs locked.
#[derive(Default)]
pub(super) struct FileImageCache {
    entries: VecDeque<(FileImageKey, WeakFileImage)>,
}

impl FileImageCache {
    /// Upgrade a candidate, pruning dead hints and moving a live hit to the end.
    pub(super) fn retained(&mut self, key: &FileImageKey) -> Option<SpkRegistry> {
        self.entries.retain(|(_, image)| image.live());
        let ordinal = self.entries.iter().position(|(stored, _)| stored == key)?;
        let entry = self.entries.remove(ordinal)?;
        let retained = entry.1.retained(&key.path);
        if retained.is_some() {
            self.entries.push_back(entry);
        }
        retained
    }

    /// Store weak allocation handles after successful stable-file admission.
    /// Dead hints retire first; live oldest hints can be evicted at the bound
    /// without changing any registry, file or active immutable owner.
    pub(super) fn remember(&mut self, key: FileImageKey, source: &SpkRegistry) {
        let Some(kernel) = source.kernel(&key.path) else {
            return;
        };
        let Some(decoded) = source.decoded.get(&key.path) else {
            return;
        };
        let Some(owner) = &kernel.file_image_owner else {
            return;
        };
        self.entries
            .retain(|(stored, image)| stored != &key && image.live());
        if self.entries.len() == MAX_FILE_IMAGE_HINTS {
            self.entries.pop_front();
        }
        self.entries.push_back((
            key,
            WeakFileImage {
                owner: Arc::downgrade(owner),
                decoded: Arc::downgrade(decoded),
            },
        ));
    }

    /// Number of retained bookkeeping entries, for bounded-ownership controls.
    #[cfg(test)]
    pub(super) fn len(&self) -> usize {
        self.entries.len()
    }
}

/// Process-wide native weak hints; no immutable image or file is rooted here.
fn images() -> &'static Mutex<FileImageCache> {
    static IMAGES: OnceLock<Mutex<FileImageCache>> = OnceLock::new();
    IMAGES.get_or_init(|| Mutex::new(FileImageCache::default()))
}

/// Acquire a candidate under a short bookkeeping lock, released before I/O.
pub(super) fn retained(key: &FileImageKey) -> Option<SpkRegistry> {
    images()
        .lock()
        .unwrap_or_else(|error| error.into_inner())
        .retained(key)
}

/// Publish weak storage after successful decode or complete byte comparison.
pub(super) fn remember(key: FileImageKey, source: &SpkRegistry) {
    images()
        .lock()
        .unwrap_or_else(|error| error.into_inner())
        .remember(key, source);
}

/// Compare the actual descriptor stream with immutable bytes using 64 KiB of
/// scratch storage. A match consumes EOF and checks exact length. A mismatch
/// returns before admission so the caller can seek back and read/decode the
/// changed image. Interrupted reads retry; all other I/O errors remain errors.
pub(super) fn file_matches_image(file: &mut File, expected: &[u8]) -> io::Result<bool> {
    let mut buffer = [0u8; 64 * 1024];
    let mut consumed = 0;
    loop {
        let count = match file.read(&mut buffer) {
            Err(error) if error.kind() == io::ErrorKind::Interrupted => continue,
            result => result?,
        };
        if count == 0 {
            return Ok(consumed == expected.len());
        }
        if expected.get(consumed..consumed + count) != Some(&buffer[..count]) {
            return Ok(false);
        }
        consumed += count;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn immutable_digest_is_shared_by_concurrent_readers_and_has_fixed_owner_cost() {
        use sha2::{Digest, Sha256};
        use std::sync::Barrier;

        // The existing synthetic SPK helper supplies valid bytes, not a file or
        // an astronomical reference. Byte-backed kernels have no file carrier.
        let bytes = crate::spk::moon_test_kernel_bytes();
        let expected: [u8; 32] = Sha256::digest(&bytes).into();
        let kernel = SpkKernel::decode("digest-synthetic", bytes).unwrap();
        assert!(kernel.file_image_owner.is_none());
        assert_eq!(kernel.image_sha256(), expected);
        let owner = FileImageOwner::acquired(&kernel);
        assert!(owner.sha256.get().is_none());
        let barrier = Barrier::new(8);
        std::thread::scope(|scope| {
            let handles: Vec<_> = (0..8)
                .map(|_| {
                    scope.spawn(|| {
                        barrier.wait();
                        owner.sha256(kernel.bytes())
                    })
                })
                .collect();
            let results: Vec<_> = handles
                .into_iter()
                .map(|handle| handle.join().unwrap())
                .collect();
            let published = owner.sha256(kernel.bytes());
            for result in results {
                assert_eq!(*result, expected);
                assert!(std::ptr::eq(result, published));
            }
        });
        assert_eq!(owner.sha256.get(), Some(&expected));
        assert!(std::ptr::eq(
            owner.sha256(kernel.bytes()),
            owner.sha256(kernel.bytes())
        ));
        // Report the compiled storage delta rather than predicting padding.
        struct PreviousFileImageOwner {
            _bytes: Weak<[u8]>,
            _index: Weak<SpkIndex>,
            _endian: DafEndian,
        }
        let previous = std::mem::size_of::<PreviousFileImageOwner>();
        let current = std::mem::size_of::<FileImageOwner>();
        let digest = std::mem::size_of::<OnceLock<[u8; 32]>>();
        assert!(current >= previous);
        assert!(current <= previous + digest + std::mem::align_of::<FileImageOwner>());
        println!(
            "SPK_IMAGE_DIGEST_OWNER_BYTES previous={previous} current={current} delta={} digest={digest}",
            current - previous
        );
    }
}
