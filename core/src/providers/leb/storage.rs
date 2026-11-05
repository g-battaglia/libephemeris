//! Owned bytes and caller-guaranteed immutable file maps. All OS advice is
//! best effort: residency and filesystem cache are not an application RAM cap.
use super::{Error, Result};
use std::fs::{File, OpenOptions};
use std::path::Path;

/// Backing and descriptor stay with one image; no pathname reopening for advice.
pub(super) enum ImageStorage {
    Owned(Vec<u8>),
    /// Immutable executable data, shared without a second heap image.
    Static(&'static [u8]),
    /// Shared immutable snapshot retained across strict and Source admission.
    SharedOwned(std::sync::Arc<Vec<u8>>),
    Positional(Box<super::access::PositionalFile>),
    Mapped {
        map: memmap2::Mmap,
        file: File,
    },
}

/// Open read-only without blocking on a FIFO, then reject non-regular files.
/// Symlinks to regular files are allowed; identity belongs to the open inode.
pub(crate) fn open_regular_file(path: &Path) -> std::io::Result<File> {
    let mut options = OpenOptions::new();
    options.read(true);
    #[cfg(unix)]
    {
        use std::os::unix::fs::OpenOptionsExt;
        options.custom_flags(rustix::fs::OFlags::NONBLOCK.bits() as i32);
    }
    let file = options.open(path)?;
    if !file.metadata()?.is_file() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "LEB data must be a regular file",
        ));
    }
    Ok(file)
}

impl ImageStorage {
    /// Map an immutable regular file without pre-faulting its pages.
    ///
    /// # Safety
    /// The caller must prevent writes/truncation of this inode until every
    /// image owner and borrowed view has gone away. Read-only permissions or
    /// MAP_PRIVATE alone do not establish that invariant.
    pub(super) unsafe fn open_mapped(path: &Path) -> Result<Self> {
        let file = open_regular_file(path).map_err(|error| {
            if error.kind() == std::io::ErrorKind::NotFound {
                Error::LebFileNotFound {
                    path: path.to_string_lossy().into_owned(),
                }
            } else {
                Error::LebCorrupt {
                    what: format!("LEB file open: {error}"),
                    offset: 0,
                }
            }
        })?;
        let len = file
            .metadata()
            .map_err(|error| Error::LebCorrupt {
                what: format!("LEB file metadata: {error}"),
                offset: 0,
            })?
            .len();
        if len == 0 || len > isize::MAX as u64 {
            return Err(Error::LebCorrupt {
                what: "LEB mapping length".into(),
                offset: 0,
            });
        }
        // SAFETY: open_mapped's caller guarantees a stable immutable inode for
        // the mapping lifetime. The descriptor is read-only and regular; no
        // mutable map, unmap, remap or extended-lifetime slice is exposed.
        let map = unsafe { memmap2::MmapOptions::new().map(&file) }.map_err(|error| {
            Error::LebCorrupt {
                what: format!("LEB file mapping: {error}"),
                offset: 0,
            }
        })?;
        Ok(Self::Mapped { map, file })
    }

    pub(super) fn source(&self) -> super::LebSource<'_> {
        match self {
            Self::Owned(bytes) => super::LebSource::Slice(bytes),
            Self::Static(bytes) => super::LebSource::Slice(bytes),
            Self::SharedOwned(bytes) => super::LebSource::Slice(bytes),
            Self::Mapped { map, .. } => super::LebSource::Slice(map),
            Self::Positional(file) => super::LebSource::Positional(file),
        }
    }

    pub(super) fn bytes(&self) -> Result<&[u8]> {
        match self.source() {
            super::LebSource::Slice(bytes) => Ok(bytes),
            super::LebSource::Positional(_) => Err(Error::LebBytesUnavailable),
        }
    }
    pub(super) fn owned_capacity(&self) -> usize {
        match self {
            Self::Owned(bytes) => bytes.capacity(),
            Self::SharedOwned(bytes) => bytes.capacity(),
            Self::Mapped { .. } | Self::Static(_) | Self::Positional(_) => 0,
        }
    }
    pub(super) fn mapped_len(&self) -> usize {
        match self {
            Self::Owned(_) | Self::SharedOwned(_) | Self::Static(_) | Self::Positional(_) => 0,
            Self::Mapped { map, .. } => map.len(),
        }
    }

    /// Advisory prefetch, not eager evaluation or whole-file warming.
    pub(super) fn warm(&self, ranges: &[(u64, u64)]) {
        match self {
            Self::Owned(bytes) => super::prefault(bytes, ranges),
            Self::Static(bytes) => super::prefault(bytes, ranges),
            Self::SharedOwned(bytes) => super::prefault(bytes, ranges),
            Self::Positional(_) => {}
            Self::Mapped { map, .. } => {
                #[cfg(unix)]
                for &(offset, size) in ranges {
                    if let Some((start, len)) =
                        aligned_range(offset, size, map.len(), rustix::param::page_size())
                    {
                        let _ = map.advise_range(memmap2::Advice::WillNeed, start, len);
                    }
                }
                #[cfg(not(unix))]
                let _ = (map, ranges); // Like Python: no advice where unavailable.
            }
        }
    }

    /// Drop process page references, then ask Linux to evict clean file cache.
    /// Existing borrows remain valid; re-faults read the identical immutable data.
    /// Not guaranteed to release every page, especially with other active readers.
    pub(super) fn release_pages(&self) {
        if let Self::Positional(file) = self {
            file.release_pages();
        }
        if let Self::Mapped { map, file } = self {
            #[cfg(unix)]
            {
                // SAFETY: only a read-only file map of an immutable inode. No
                // dirty private/anonymous pages can be lost or zero-filled.
                let _ = unsafe { map.unchecked_advise(memmap2::UncheckedAdvice::DontNeed) };
            }
            #[cfg(not(unix))]
            let _ = map;
            release_file_cache(file);
        }
    }
}

/// Page-aligned bounded interval, with overflow-safe clamping at file end.
#[cfg(unix)]
fn aligned_range(
    offset: u64,
    size: u64,
    file_len: usize,
    page_size: usize,
) -> Option<(usize, usize)> {
    if size == 0 || offset >= file_len as u64 || page_size == 0 {
        return None;
    }
    let start = offset as usize / page_size * page_size;
    let end = offset.saturating_add(size).min(file_len as u64) as usize;
    let end = end.saturating_add(page_size - 1) / page_size * page_size;
    let end = end.min(file_len);
    Some((start, end - start))
}

/// Native counterpart of POSIX_FADV_DONTNEED. Zero length means through EOF.
/// No-op outside Linux, on unsupported operating systems.
pub(crate) fn release_file_cache(file: &File) {
    #[cfg(target_os = "linux")]
    let _ = rustix::fs::fadvise(file, 0, None, rustix::fs::Advice::DontNeed);
    #[cfg(not(target_os = "linux"))]
    let _ = file;
}

/// Retire only fully consumed pages during a sequential integrity scan.
/// Advice is best effort and never changes the digest or read error semantics.
/// Ordinary coefficient reads do not use this helper, preserving their hot cache.
pub(crate) struct HashCacheAdvice<'a> {
    file: &'a File,
    consumed: u64,
    released: u64,
}

impl<'a> HashCacheAdvice<'a> {
    pub(crate) fn new(file: &'a File) -> Self {
        Self {
            file,
            consumed: 0,
            released: 0,
        }
    }

    pub(crate) fn consumed(&mut self, offset: u64) {
        self.consumed = offset;
        if self.consumed.saturating_sub(self.released) >= 8 * 1024 * 1024 {
            self.release();
        }
    }

    fn release(&mut self) {
        #[cfg(target_os = "linux")]
        {
            let page = rustix::param::page_size() as u64;
            let end = self.consumed / page * page;
            if let Some(len) = std::num::NonZeroU64::new(end.saturating_sub(self.released)) {
                let _ = rustix::fs::fadvise(
                    self.file,
                    self.released,
                    Some(len),
                    rustix::fs::Advice::DontNeed,
                );
                self.released = end;
            }
        }
        #[cfg(not(target_os = "linux"))]
        {
            let _ = self.file;
            self.released = self.consumed;
        }
    }
}

impl Drop for HashCacheAdvice<'_> {
    fn drop(&mut self) {
        // Also release the completed prefix on a failed read, without touching
        // an unread range or the unaligned tail of the final page.
        self.release();
    }
}

/// Best-effort release of filesystem cache below an explicitly supplied data
/// directory. No file contents are read, changed or deleted. Directory symlinks
/// are not descended (no cycles); file symlinks to regular files follow stat/open
/// semantics. Per-entry errors are ignored. No-op where fadvise is unavailable.
/// This does not clear coefficient caches owned by live readers.
pub fn release_data_cache(directory: &Path) {
    #[cfg(target_os = "linux")]
    {
        let mut pending = vec![directory.to_path_buf()];
        while let Some(dir) = pending.pop() {
            let Ok(entries) = std::fs::read_dir(dir) else {
                continue;
            };
            for entry in entries.flatten() {
                let path = entry.path();
                let Ok(kind) = entry.file_type() else {
                    continue;
                };
                if kind.is_dir() {
                    pending.push(path);
                } else if let Ok(file) = open_regular_file(&path) {
                    release_file_cache(&file);
                }
            }
        }
    }
    #[cfg(not(target_os = "linux"))]
    let _ = directory;
}

#[cfg(all(test, unix))]
mod tests {
    use super::*;
    #[test]
    fn advice_ranges_use_host_page_alignment_without_overflow() {
        for page in [4096, 16384] {
            assert_eq!(aligned_range(1, 2, 100000, page), Some((0, page)));
            assert_eq!(
                aligned_range(99999, u64::MAX, 100000, page),
                Some((98304, 1696))
            );
            assert_eq!(aligned_range(u64::MAX, 10, 100000, page), None);
            assert_eq!(aligned_range(0, 0, 100000, page), None);
        }
    }
}
