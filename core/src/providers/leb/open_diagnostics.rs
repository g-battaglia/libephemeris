//! Reader-constructor diagnostics over one retained image snapshot.
//!
//! This controller inspects metadata only, in the released LEB1/LEB2 parse
//! order. It does not decode coefficients, validate numerical table domains,
//! relax LebFile invariants or fabricate a metadata failure for a native-only
//! safety refusal. Offsets/counts retain their exact integer values.

/// Container reader selected by the preceding four-byte magic probe.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReaderKind {
    /// Uncompressed LEB1 container.
    Leb1,
    /// Chunk-compressed LEB2 container.
    Leb2,
}

/// First metadata read or consistency check that refuses the acquired image.
/// This diagnostic does not certify coefficient values or numerical coverage.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ReaderParseDiagnostic {
    /// Mmap precedes all parsing and refuses an empty selected image.
    EmptyImage,
    /// Header magic differs from the selected reader class.
    InvalidMagic {
        /// Four bytes read from the start of the selected image.
        actual: [u8; 4],
    },
    /// Actual version field, before compression or directory acquisition.
    UnsupportedVersion {
        /// Version value read from the container header.
        found: u32,
    },
    /// Actual LEB2 header compression flags.
    UnsupportedCompression {
        /// Compression flags read from the LEB2 header.
        found: u32,
    },
    /// Fixed-size metadata read extending beyond the acquired image.
    Unpack {
        /// Wire-layout notation identifying the requested fixed-size record.
        format: &'static str,
        /// Required record length, in bytes.
        size: u32,
        /// Requested byte offset from the start of the image.
        offset: u128,
        /// Total length of the acquired image, in bytes.
        image_len: usize,
    },
    /// Metadata offset exceeds the platform's signed indexing range.
    OffsetOverflow {
        /// Wire-layout notation of the record whose read was refused.
        format: &'static str,
        /// Requested byte offset, retained without narrowing.
        offset: u128,
    },
    /// Declared records exceed the space available in their section.
    EntryCount {
        /// Number of entries declared by the section or index.
        count: u32,
        /// Required length of each entry, in bytes.
        entry_size: u32,
        /// Bytes remaining after fixed headers; negative for undersized sections.
        available: i128,
        /// Human-readable section or index label, including a body where applicable.
        what: String,
    },
}

const HEADER: (&str, u32) = ("<4sIIIdddI20s", 64);
const SECTION: (&str, u32) = ("<IIQQ", 24);
const BODY: (&str, u32) = ("<iIIdddIIQ", 52);
const COMPRESSED_BODY: (&str, u32) = ("<iIIdddIIQQQ", 68);
const CHUNK_HEADER: (&str, u32) = ("<IIfI", 16);
const CHUNK: (&str, u32) = ("<ddQQQII", 48);
const NUTATION: (&str, u32) = ("<dddIIII", 40);
const DELTA_HEADER: (&str, u32) = ("<II", 8);
const DELTA_ENTRY: (&str, u32) = ("<dd", 16);
const STAR: (&str, u32) = ("<iddddddd4s", 64);

fn read<'a>(
    bytes: &'a [u8],
    offset: u128,
    spec: (&'static str, u32),
) -> Result<&'a [u8], ReaderParseDiagnostic> {
    if offset > isize::MAX as u128 {
        return Err(ReaderParseDiagnostic::OffsetOverflow {
            format: spec.0,
            offset,
        });
    }
    let end = offset + u128::from(spec.1);
    if end > bytes.len() as u128 {
        return Err(ReaderParseDiagnostic::Unpack {
            format: spec.0,
            size: spec.1,
            offset,
            image_len: bytes.len(),
        });
    }
    Ok(&bytes[offset as usize..end as usize])
}
fn u32_at(bytes: &[u8], offset: usize) -> u32 {
    u32::from_le_bytes(
        bytes[offset..offset + 4]
            .try_into()
            .expect("checked metadata read"),
    )
}
fn i32_at(bytes: &[u8], offset: usize) -> i32 {
    i32::from_le_bytes(
        bytes[offset..offset + 4]
            .try_into()
            .expect("checked metadata read"),
    )
}
fn u64_at(bytes: &[u8], offset: usize) -> u64 {
    u64::from_le_bytes(
        bytes[offset..offset + 8]
            .try_into()
            .expect("checked metadata read"),
    )
}
fn count(
    count: u32,
    size: u32,
    available: i128,
    what: String,
) -> Result<(), ReaderParseDiagnostic> {
    if i128::from(count) * i128::from(size) > available {
        return Err(ReaderParseDiagnostic::EntryCount {
            count,
            entry_size: size,
            available,
            what,
        });
    }
    Ok(())
}

fn parse(bytes: &[u8], kind: ReaderKind) -> Result<(), ReaderParseDiagnostic> {
    if bytes.is_empty() {
        return Err(ReaderParseDiagnostic::EmptyImage);
    }
    let header = read(bytes, 0, HEADER)?;
    let magic: [u8; 4] = header[..4].try_into().expect("four header bytes");
    let expected = match kind {
        ReaderKind::Leb1 => *b"LEB1",
        ReaderKind::Leb2 => *b"LEB2",
    };
    if magic != expected {
        return Err(ReaderParseDiagnostic::InvalidMagic { actual: magic });
    }
    let version = u32_at(header, 4);
    if (kind == ReaderKind::Leb1 && version != 1)
        || (kind == ReaderKind::Leb2 && !matches!(version, 1 | 2))
    {
        return Err(ReaderParseDiagnostic::UnsupportedVersion { found: version });
    }
    if kind == ReaderKind::Leb2 {
        let found = u32_at(header, 40);
        if found != 1 {
            return Err(ReaderParseDiagnostic::UnsupportedCompression { found });
        }
    }
    // Source dictionaries overwrite a repeated section while preserving its
    // original insertion location; only keyed selections matter below.
    let mut sections = Vec::<(u32, u64, u64)>::new();
    for index in 0..u32_at(header, 8) {
        let row = read(bytes, 64 + u128::from(index) * 24, SECTION)?;
        let value = (u32_at(row, 0), u64_at(row, 8), u64_at(row, 16));
        if let Some(old) = sections.iter_mut().find(|old| old.0 == value.0) {
            *old = value;
        } else {
            sections.push(value);
        }
    }
    let section = |id| {
        sections
            .iter()
            .find(|entry| entry.0 == id)
            .map(|entry| (entry.1, entry.2))
    };
    let mut bodies = Vec::<(i32, u64)>::new();
    if let Some((offset, size)) = section(0) {
        let body_count = u32_at(header, 12);
        let (spec, label) = match kind {
            ReaderKind::Leb1 => (BODY, "LEB body index"),
            ReaderKind::Leb2 => (COMPRESSED_BODY, "LEB2 body index"),
        };
        count(body_count, spec.1, i128::from(size), label.to_owned())?;
        for index in 0..body_count {
            let row = read(
                bytes,
                u128::from(offset) + u128::from(index) * u128::from(spec.1),
                spec,
            )?;
            let value = (i32_at(row, 0), u64_at(row, 44));
            if let Some(old) = bodies.iter_mut().find(|old| old.0 == value.0) {
                *old = value;
            } else {
                bodies.push(value);
            }
        }
    }
    if kind == ReaderKind::Leb2 && version >= 2 {
        for (body, offset) in bodies {
            let row = read(bytes, u128::from(offset), CHUNK_HEADER)?;
            let chunk_count = u32_at(row, 0);
            let first = u128::from(offset) + 16;
            count(
                chunk_count,
                48,
                bytes.len() as i128 - first as i128,
                format!("LEB2 chunk index (body {body})"),
            )?;
            for index in 0..chunk_count {
                read(bytes, first + u128::from(index) * 48, CHUNK)?;
            }
        }
    }
    if let Some((offset, _)) = section(2) {
        read(bytes, u128::from(offset), NUTATION)?;
    }
    if let Some((offset, size)) = section(3) {
        let row = read(bytes, u128::from(offset), DELTA_HEADER)?;
        let entries = u32_at(row, 0);
        let label = match kind {
            ReaderKind::Leb1 => "LEB Delta-T table",
            ReaderKind::Leb2 => "LEB2 Delta-T table",
        };
        count(entries, 16, i128::from(size) - 8, label.to_owned())?;
        for index in 0..entries {
            read(
                bytes,
                u128::from(offset) + 8 + u128::from(index) * 16,
                DELTA_ENTRY,
            )?;
        }
    }
    if let Some((offset, size)) = section(4) {
        for index in 0..size / 64 {
            // The fourth-byte field is reserved. Source read_star_entry
            // unpacks but never decodes or publishes it.
            read(bytes, u128::from(offset) + u128::from(index) * 64, STAR)?;
        }
    }
    Ok(())
}

/// First actual constructor diagnostic, or None when ordinary Source metadata
/// parsing succeeds. The caller retains this same acquired image/path. None
/// does not waive strict native validation: unsupported numerical metadata
/// requires separate lazy numerical admission in the Source-facing owner.
pub fn source_constructor_failure(bytes: &[u8], kind: ReaderKind) -> Option<ReaderParseDiagnostic> {
    parse(bytes, kind).err()
}

#[cfg(test)]
mod tests {
    use super::*;
    fn header(kind: ReaderKind, sections: u32, bodies: u32) -> Vec<u8> {
        let mut bytes = vec![0; 64];
        bytes[..4].copy_from_slice(match kind {
            ReaderKind::Leb1 => b"LEB1",
            ReaderKind::Leb2 => b"LEB2",
        });
        bytes[4..8].copy_from_slice(&1u32.to_le_bytes());
        bytes[8..12].copy_from_slice(&sections.to_le_bytes());
        bytes[12..16].copy_from_slice(&bodies.to_le_bytes());
        bytes[40..44].copy_from_slice(&1u32.to_le_bytes());
        bytes
    }
    #[test]
    fn directory_truncation_precedes_body_fields() {
        let bytes = header(ReaderKind::Leb1, 1, 100);
        assert_eq!(
            source_constructor_failure(&bytes, ReaderKind::Leb1),
            Some(ReaderParseDiagnostic::Unpack {
                format: SECTION.0,
                size: 24,
                offset: 64,
                image_len: 64
            })
        );
    }
    #[test]
    fn absent_body_section_is_not_an_invented_source_failure() {
        let bytes = header(ReaderKind::Leb1, 0, 10);
        assert_eq!(source_constructor_failure(&bytes, ReaderKind::Leb1), None);
    }
    #[test]
    fn declared_count_preserves_negative_delta_available() {
        let mut bytes = header(ReaderKind::Leb2, 1, 0);
        bytes.resize(96, 0);
        bytes[64..68].copy_from_slice(&3u32.to_le_bytes());
        bytes[72..80].copy_from_slice(&88u64.to_le_bytes());
        assert_eq!(
            source_constructor_failure(&bytes, ReaderKind::Leb2),
            Some(ReaderParseDiagnostic::EntryCount {
                count: 0,
                entry_size: 16,
                available: -8,
                what: "LEB2 Delta-T table".into()
            })
        );
    }
}
