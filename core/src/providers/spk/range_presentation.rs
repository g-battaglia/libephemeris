//! Owned segment metadata for a scalar SPK coverage failure.
//!
//! Capture occurs after actual descriptor selection and directory decoding,
//! before the typed error leaves its provider. No dates, polynomial words or
//! backend attributes are inferred from a diagnostic string. The image handle
//! keeps the selected immutable bytes alive independently of registry mutation.

use super::{ChebSegment, Descriptor, Segment, SpkKernel};
use sha2::{Digest, Sha256};
use std::sync::Arc;

/// An acquired immutable image, with a content identity computed from its bytes.
/// The path is the registry identity, not proof of a later file's contents.
#[derive(Clone)]
pub struct RangeImage {
    /// Registry identity at acquisition; the pathname is not reopened for inspection.
    pub path: String,
    /// SHA-256 of the acquired bytes as 32 raw bytes, computed at acquisition.
    pub sha256: [u8; 32],
    bytes: Arc<[u8]>,
    kernel: SpkKernel,
}

/// Diagnostic identity omits large image contents and coefficient bytes.
impl std::fmt::Debug for RangeImage {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("RangeImage")
            .field("path", &self.path)
            .field("sha256", &hex(&self.sha256))
            .field("bytes", &self.bytes.len())
            .finish()
    }
}

/// Equality retains the same immutable image owner, not a bytewise comparison
/// of a potentially large resource. Path/digest/length remain observable facts;
/// cloning a captured result shares the owner and never copies image contents.
impl PartialEq for RangeImage {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.bytes, &other.bytes)
            && self.path == other.path
            && self.sha256 == other.sha256
            && self.bytes.len() == other.bytes.len()
    }
}

impl RangeImage {
    /// Duplicate the retained acquisition descriptor, if one was supplied.
    /// No path is reopened; in-memory images have no file capability. The
    /// duplicate owns its lifetime and can outlive registry replacement, but
    /// may share a cursor with other duplicates. Access permissions are those
    /// of the original descriptor; see [`SpkKernel::decode_file_image`].
    /// Numerical evaluation uses immutable bytes, not this cursor or file.
    pub fn retained_source_file(&self) -> std::io::Result<Option<std::fs::File>> {
        self.kernel.retained_source_file()
    }
    /// Own the actual installed bytes and their independently computed identity.
    pub(crate) fn acquired(kernel: &SpkKernel) -> Self {
        Self {
            path: kernel.path.clone(),
            sha256: Sha256::digest(kernel.bytes()).into(),
            bytes: Arc::clone(&kernel.bytes),
            kernel: kernel.clone(),
        }
    }

    /// Actual captured bytes. No registry, file handle or provider borrow escapes.
    pub fn bytes(&self) -> &[u8] {
        &self.bytes
    }

    /// Retain this same acquired image and decoded index without reopening or
    /// copying its bytes. Callers can inspect the selected
    /// physical segment independently of later registry/cache replacement.
    pub fn retained_kernel(&self) -> SpkKernel {
        self.kernel.clone()
    }
}

impl SpkKernel {
    /// Decode an owned byte image and retain the supplied acquisition file.
    ///
    /// The caller supplies bytes consumed from `file`; no pathname is opened
    /// and their correspondence with current file contents is not verified here.
    /// Numerical reads use the immutable bytes, independently of file mutation
    /// or cursor changes. Structural errors follow [`Self::decode`].
    ///
    /// File-based registry acquisition supplies a read-only file. This public
    /// constructor does not enforce that access mode; later descriptor clones
    /// inherit the supplied file's capabilities and may share its cursor.
    pub fn decode_file_image(
        path: &str,
        bytes: Vec<u8>,
        file: std::fs::File,
    ) -> crate::errors::Result<Self> {
        let mut kernel = Self::decode(path, bytes)?;
        kernel.file_image_owner = Some(Arc::new(
            super::file_image_cache::FileImageOwner::acquired(&kernel),
        ));
        kernel.source_file = Some(Arc::new(file));
        Ok(kernel)
    }

    /// Duplicate the retained acquisition file, or return `Ok(None)` if absent.
    /// No current pathname is resolved or opened. The returned handle owns its
    /// lifetime, can share the original cursor, and retains the access mode
    /// supplied at acquisition. Duplication failures return ordinary I/O errors.
    /// Its current contents need not equal the separately retained byte snapshot.
    pub fn retained_source_file(&self) -> std::io::Result<Option<std::fs::File>> {
        self.source_file
            .as_ref()
            .map(|file| file.try_clone())
            .transpose()
    }
}

/// Selected descriptor and decoded directory of one actual scalar refusal.
/// Epoch and directory floating-point values are retained as binary64 bit
/// patterns, decoded with [`f64::from_bits`]. The two epoch parts add to a TDB
/// Julian date; the directory uses seconds from J2000 and seconds per record.
/// Type 2/3 Chebyshev directories supply these fields. Retaining this diagnostic
/// keeps the entire shared image alive, not just the selected segment.
#[derive(Clone, Debug, PartialEq)]
pub struct OwnedSpkRange {
    /// Shared owner of the bytes actually used for descriptor and directory reads.
    pub image: RangeImage,
    /// File-order ordinal, required even when descriptors compare equal.
    pub segment_ordinal: usize,
    /// Selected physical DAF descriptor, before directory coverage was checked.
    pub descriptor: Descriptor,
    /// The DAF summary name after removing edge HT/LF/VT/FF/CR/space bytes.
    pub source: Vec<u8>,
    /// Binary64 bits of the first additive TDB Julian-date part; not an integer day.
    pub tdb_whole: u64,
    /// Binary64 bits of the second additive TDB Julian-date part.
    pub tdb_fraction: u64,
    /// Binary64 bits of the directory's initial epoch, in seconds from J2000 TDB.
    pub init_et: u64,
    /// Binary64 bits of the duration of each directory record, in seconds.
    pub interval_seconds: u64,
    /// Number of 64-bit words per record, including midpoint and half-width.
    pub record_words: usize,
    /// Number of Chebyshev records declared by the selected segment directory.
    pub record_count: usize,
    /// Coefficients per component: `(record_words - 2) / component_count`.
    pub coefficient_count: usize,
    /// Three position components for type 2; six position/velocity components for type 3.
    pub component_count: usize,
}

impl OwnedSpkRange {
    /// Retain the actual selected image/segment only after its directory refused.
    pub(super) fn acquired(
        kernel: &SpkKernel,
        selected: &Segment,
        decoded: &ChebSegment,
        jd: f64,
        fraction: f64,
    ) -> Self {
        let components = if decoded.seg_type == 2 { 3 } else { 6 };
        let segment_ordinal = kernel
            .index
            .segments
            .iter()
            .position(|segment| std::ptr::eq(segment, selected))
            .expect("selected segment belongs to the acquired kernel");
        let source = strip_daf_summary_name(&selected.name).to_vec();
        Self {
            image: RangeImage::acquired(kernel),
            segment_ordinal,
            descriptor: selected.descriptor,
            source,
            tdb_whole: jd.to_bits(),
            tdb_fraction: fraction.to_bits(),
            init_et: decoded.init_et.to_bits(),
            interval_seconds: decoded.int_len.to_bits(),
            record_words: decoded.record_words,
            record_count: decoded.n_records,
            coefficient_count: (decoded.record_words - 2) / components,
            component_count: components,
        }
    }
}

/// Remove edge HT, LF, VT, FF, CR and space bytes from a DAF summary name.
/// Interior bytes and every other byte remain unchanged; this performs no
/// Unicode decoding. The returned slice borrows the input, including an empty
/// slice for an empty name or a name containing only the six removable bytes.
fn strip_daf_summary_name(bytes: &[u8]) -> &[u8] {
    let is_edge_space = |byte| matches!(byte, b'\t' | b'\n' | 0x0b | 0x0c | b'\r' | b' ');
    let mut start = 0;
    let mut end = bytes.len();
    while start < end && is_edge_space(bytes[start]) {
        start += 1;
    }
    while end > start && is_edge_space(bytes[end - 1]) {
        end -= 1;
    }
    &bytes[start..end]
}

/// Lowercase hex encodes actual bytes, including opaque DAF names.
fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|byte| format!("{byte:02x}")).collect()
}

#[cfg(test)]
mod tests {
    use super::strip_daf_summary_name;

    #[test]
    fn daf_name_edges_remove_exactly_six_byte_values() {
        for byte in 0..=u8::MAX {
            let name = [byte, b'A', byte];
            let removable = matches!(byte, 0x09..=0x0d | 0x20);
            let expected: &[u8] = if removable { b"A" } else { &name };
            assert_eq!(strip_daf_summary_name(&name), expected, "byte {byte:02x}");
        }
        assert!(strip_daf_summary_name(b"").is_empty());
        assert!(strip_daf_summary_name(b"\t\n\x0b\x0c\r ").is_empty());
    }

    #[test]
    fn daf_name_interior_and_nonspace_bytes_are_preserved() {
        let interior = b"A\t\n\x0b\x0c\r \x00\x85\xa0\xffZ";
        let mut padded = b"\t\n\x0b\x0c\r ".to_vec();
        padded.extend_from_slice(interior);
        padded.extend_from_slice(b" \r\x0c\x0b\n\t");
        assert_eq!(strip_daf_summary_name(&padded), interior);
        assert_eq!(strip_daf_summary_name(interior), interior);
    }
}
