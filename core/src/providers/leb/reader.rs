//! LEB header decoding, coefficient transforms, inspection and image-backed readers.
//!
//! Body and nutation evaluation use TT Julian days. Constructors describe the
//! backing ownership and resource limits; payloads are also validated when read.

use super::*;

impl LebFormat {
    /// Bytes needed to classify a header: magic, version, and flags live
    /// at offsets `0..4`, `4..8`, and `40..44`.
    const DETECT_LEN: usize = 44;

    /// Classify the first 44 bytes by magic, version and compression flags.
    ///
    /// Unknown magic, unsupported versions and unsupported LEB2 compression
    /// return their distinct format errors. A shorter prefix returns
    /// [`Error::LebCorrupt`]. This does not validate the complete container.
    pub fn detect(bytes: &[u8]) -> Result<Self> {
        let head = slice_at(bytes, 0, Self::DETECT_LEN, "LEB header")?;
        let magic: [u8; 4] = head[0..4].try_into().map_err(|_| Error::LebCorrupt {
            what: "LEB header".to_string(),
            offset: 0,
        })?;
        let version = u32::from_le_bytes(head[4..8].try_into().map_err(|_| Error::LebCorrupt {
            what: "LEB header".to_string(),
            offset: 4,
        })?);
        let flags = u32::from_le_bytes(head[40..44].try_into().map_err(|_| Error::LebCorrupt {
            what: "LEB header".to_string(),
            offset: 40,
        })?);
        Self::classify(magic, version, flags)
    }

    /// Classify decoded magic, version and flags without reading a file.
    /// The version is checked before compression flags. LEB1 ignores `flags`;
    /// LEB2 requires `COMPRESSION_ZSTD_TRUNC_SHUFFLE`. Invalid identity fields
    /// return their distinct format errors.
    pub fn classify(magic: [u8; 4], version: u32, flags: u32) -> Result<Self> {
        if magic == MAGIC_LEB1 {
            if version == VERSION_LEB1 {
                Ok(Self::Leb1)
            } else {
                Err(Error::LebUnsupportedVersion { found: version })
            }
        } else if magic == MAGIC_LEB2 {
            if version != VERSION_LEB2_V1 && version != VERSION_LEB2_V2 {
                Err(Error::LebUnsupportedVersion { found: version })
            } else if flags != COMPRESSION_ZSTD_TRUNC_SHUFFLE {
                Err(Error::LebUnsupportedCompression { found: flags })
            } else if version == VERSION_LEB2_V1 {
                Ok(Self::Leb2V1)
            } else {
                Ok(Self::Leb2V2)
            }
        } else {
            Err(Error::LebUnknownFormat { magic })
        }
    }

    /// Body index entry size for this flavor (52 for LEB1, 68 for LEB2).
    pub fn body_entry_size(self) -> usize {
        match self {
            Self::Leb1 => LEB1_BODY_ENTRY_SIZE,
            Self::Leb2V1 | Self::Leb2V2 => LEB2_BODY_ENTRY_SIZE,
        }
    }
}

/// Borrow `len` bytes at `offset`, mapping truncation and address
/// overflow to [`Error::LebCorrupt`] with the failing structure label.
pub(super) fn slice_at<'a>(
    bytes: &'a [u8],
    offset: usize,
    len: usize,
    what: &str,
) -> Result<&'a [u8]> {
    let end = offset.checked_add(len).ok_or_else(|| Error::LebCorrupt {
        what: what.to_string(),
        offset,
    })?;
    bytes.get(offset..end).ok_or_else(|| Error::LebCorrupt {
        what: what.to_string(),
        offset,
    })
}

/// Copy exactly 4 bytes out of a checked slice prefix.
pub(super) fn bytes4(chunk: &[u8], what: &str, offset: usize) -> Result<[u8; 4]> {
    chunk
        .get(0..4)
        .ok_or_else(|| Error::LebCorrupt {
            what: what.to_string(),
            offset,
        })?
        .try_into()
        .map_err(|_| Error::LebCorrupt {
            what: what.to_string(),
            offset,
        })
}

impl FileHeader {
    /// Decode the 64-byte little-endian header and classify its format.
    ///
    /// Truncation returns [`Error::LebCorrupt`]; identity failures retain their
    /// specific format errors. Reserved bytes and numerical date fields are
    /// not validated here. File-level dates are metadata; body entries govern
    /// evaluation coverage.
    pub fn decode(bytes: &[u8]) -> Result<Self> {
        let chunk = slice_at(bytes, 0, HEADER_SIZE, "LEB header")?;
        let magic = bytes4(chunk, "LEB header", 0)?;
        let version = read_u32(chunk, 4, "LEB header", 0)?;
        let section_count = read_u32(chunk, 8, "LEB header", 0)?;
        let body_count = read_u32(chunk, 12, "LEB header", 0)?;
        let jd_start = read_f64(chunk, 16, "LEB header", 0)?;
        let jd_end = read_f64(chunk, 24, "LEB header", 0)?;
        let generation_epoch = read_f64(chunk, 32, "LEB header", 0)?;
        let flags = read_u32(chunk, 40, "LEB header", 0)?;
        let format = LebFormat::classify(magic, version, flags)?;
        Ok(Self {
            magic,
            version,
            section_count,
            body_count,
            jd_start,
            jd_end,
            generation_epoch,
            flags,
            format,
        })
    }

    /// Raw magic bytes.
    pub fn magic(&self) -> [u8; 4] {
        self.magic
    }
    /// Raw version field.
    pub fn version(&self) -> u32 {
        self.version
    }
    /// Declared number of section directory entries.
    pub fn section_count(&self) -> u32 {
        self.section_count
    }
    /// Declared number of body index entries.
    pub fn body_count(&self) -> u32 {
        self.body_count
    }
    /// File-level coverage start (Julian Day; informational only).
    pub fn jd_start(&self) -> f64 {
        self.jd_start
    }
    /// File-level coverage end (Julian Day; informational only).
    pub fn jd_end(&self) -> f64 {
        self.jd_end
    }
    /// Generator epoch stamp (Julian Day; informational only).
    pub fn generation_epoch(&self) -> f64 {
        self.generation_epoch
    }
    /// Raw compression-scheme flags (meaningful for LEB2 only).
    pub fn flags(&self) -> u32 {
        self.flags
    }
    /// Classified container flavor.
    pub fn format(&self) -> LebFormat {
        self.format
    }
}

impl SectionEntry {
    /// Decode a little-endian section-directory entry at byte `offset`.
    /// The four reserved bytes are ignored. Truncation returns
    /// [`Error::LebCorrupt`]; payload bounds are checked by the index/section reader.
    pub fn decode(bytes: &[u8], offset: usize) -> Result<Self> {
        let chunk = slice_at(bytes, offset, SECTION_DIR_SIZE, "LEB section directory")?;
        Ok(Self {
            section_id: read_u32(chunk, 0, "LEB section directory", offset)?,
            offset: read_u64(chunk, 8, "LEB section directory", offset)?,
            size: read_u64(chunk, 16, "LEB section directory", offset)?,
        })
    }

    /// Section id (0 body index, 2 nutation, 3 Delta-T, 4 stars, 6 LEB2
    /// blobs; 1 raw Chebyshev, 5 reserved).
    pub fn section_id(&self) -> u32 {
        self.section_id
    }
    /// Absolute byte offset of the section payload.
    pub fn offset(&self) -> u64 {
        self.offset
    }
    /// Section payload size in bytes.
    pub fn size(&self) -> u64 {
        self.size
    }
}

impl BodyEntry {
    /// Decode a body entry at byte `offset` for the supplied container format.
    ///
    /// Truncation, zero segments or components, non-finite JD bounds, a
    /// non-finite or nonpositive interval, or a data offset at or beyond
    /// `file_len` return [`Error::LebCorrupt`]. Coordinate codes are retained
    /// without checking calculation-route support. Payload extents and
    /// coefficients are validated by their readers.
    pub fn decode(bytes: &[u8], offset: usize, file_len: usize, format: LebFormat) -> Result<Self> {
        let what = match format {
            LebFormat::Leb1 => "LEB body index",
            LebFormat::Leb2V1 | LebFormat::Leb2V2 => "LEB2 body index",
        };
        let chunk = slice_at(bytes, offset, format.body_entry_size(), what)?;
        let body_id = read_i32(chunk, 0, what, offset)?;
        let coord_type = read_u32(chunk, 4, what, offset)?;
        let segment_count = read_u32(chunk, 8, what, offset)?;
        let jd_start = read_f64(chunk, 12, what, offset)?;
        let jd_end = read_f64(chunk, 20, what, offset)?;
        let interval_days = read_f64(chunk, 28, what, offset)?;
        let degree = read_u32(chunk, 36, what, offset)?;
        let components = read_u32(chunk, 40, what, offset)?;
        let data_offset = read_u64(chunk, 44, what, offset)?;
        let (compressed_size, uncompressed_size) = match format {
            LebFormat::Leb1 => (None, None),
            LebFormat::Leb2V1 | LebFormat::Leb2V2 => (
                Some(read_u64(chunk, 52, what, offset)?),
                Some(read_u64(chunk, 60, what, offset)?),
            ),
        };
        let entry = Self {
            body_id,
            coord_type,
            segment_count,
            jd_start,
            jd_end,
            interval_days,
            degree,
            components,
            data_offset,
            compressed_size,
            uncompressed_size,
        };
        entry.validate_table(what, offset, file_len)?;
        Ok(entry)
    }

    fn validate_table(&self, what: &str, offset: usize, file_len: usize) -> Result<()> {
        let labeled = format!("{what} (body {})", self.body_id);
        if self.segment_count == 0
            || !self.interval_days.is_finite()
            || self.interval_days <= 0.0
            || !self.jd_start.is_finite()
            || !self.jd_end.is_finite()
            || self.components == 0
            || self.data_offset >= file_len as u64
        {
            return Err(Error::LebCorrupt {
                what: labeled,
                offset,
            });
        }
        Ok(())
    }

    /// Expected decoded coefficient bytes:
    /// `segment_count * (degree + 1) * components * 8`.
    /// Overflow returns [`Error::LebCorrupt`]. Decoders compare this geometry
    /// with the stored uncompressed size before decoding.
    pub fn expected_uncompressed_size(&self) -> Result<u64> {
        let what = format!("LEB body index (body {})", self.body_id);
        let offset = usize::try_from(self.data_offset).unwrap_or(usize::MAX);
        let degree_plus_one = u64::from(self.degree)
            .checked_add(1)
            .ok_or(Error::LebCorrupt {
                what: what.clone(),
                offset,
            })?;
        u64::from(self.segment_count)
            .checked_mul(degree_plus_one)
            .and_then(|v| v.checked_mul(u64::from(self.components)))
            .and_then(|v| v.checked_mul(8))
            .ok_or(Error::LebCorrupt { what, offset })
    }

    /// Smooth polynomial interval containing TT Julian day `jd`.
    ///
    /// Uses the same binary64 segment selector as coefficient evaluation.
    /// Interior nominal joins belong to the following record. The bounds
    /// include exactly the representable epochs selecting this record, clipped
    /// to body coverage and the selector's finite-quotient domain. They are
    /// inclusive TT days; a polynomial join is never a smooth interior point.
    /// Returns `None` outside coverage; invalid segment metadata returns an error.
    pub fn smooth_interval(&self, jd: f64) -> Result<Option<(f64, f64)>> {
        if !(self.jd_start <= jd && jd <= self.jd_end) {
            return Ok(None);
        }
        super::segment_domain::smooth_segment_interval(
            self.jd_start,
            self.jd_end,
            self.interval_days,
            self.segment_count,
            jd,
        )
        .map(Some)
    }

    /// Body identifier in the library's public numbering.
    pub fn body_id(&self) -> i32 {
        self.body_id
    }
    /// Stored coordinate code; decoding alone does not establish support.
    pub fn coord_type(&self) -> u32 {
        self.coord_type
    }
    /// Number of Chebyshev segments in this body.
    pub fn segment_count(&self) -> u32 {
        self.segment_count
    }
    /// Inclusive per-body coverage start, TT Julian day.
    pub fn jd_start(&self) -> f64 {
        self.jd_start
    }
    /// Inclusive per-body coverage end, TT Julian day.
    pub fn jd_end(&self) -> f64 {
        self.jd_end
    }
    /// Segment length in days, used (never a nominal constant).
    pub fn interval_days(&self) -> f64 {
        self.interval_days
    }
    /// Chebyshev degree (coefficients per component: `degree + 1`).
    pub fn degree(&self) -> u32 {
        self.degree
    }
    /// Declared number of coordinate channels per segment.
    pub fn components(&self) -> u32 {
        self.components
    }
    /// Absolute offset of the coefficient blob / chunk index.
    pub fn data_offset(&self) -> u64 {
        self.data_offset
    }
    /// Stored compressed blob size (`None` for LEB1).
    pub fn compressed_size(&self) -> Option<u64> {
        self.compressed_size
    }
    /// Stored decoded size (`None` for LEB1).
    pub fn uncompressed_size(&self) -> Option<u64> {
        self.uncompressed_size
    }
}

/// True for coordinate types whose longitude channel wraps to `[0, 360)`
/// on evaluation: ecliptic (1), retired heliocentric-ecliptic (2), and
/// reserved geo-ecliptic (3). ICRS types (0, 4) pass through untouched.
pub(super) fn is_ecliptic_coord(coord_type: u32) -> bool {
    coord_type == COORD_ECLIPTIC
        || coord_type == COORD_HELIO_ECL_RETIRED
        || coord_type == COORD_GEO_ECLIPTIC
}

/// Decode little-endian binary64 words from segment-major bytes.
/// Length must be divisible by eight or [`Error::LebCorrupt`] is returned.
/// Bit patterns are read without floating-point arithmetic.
pub(super) fn f64_words(data: &[u8], what: &str, offset: usize) -> Result<Vec<f64>> {
    if data.len() % 8 != 0 {
        return Err(Error::LebCorrupt {
            what: what.to_string(),
            offset,
        });
    }
    let mut out = Vec::with_capacity(data.len() / 8);
    for chunk in data.chunks_exact(8) {
        let word: [u8; 8] = chunk.try_into().map_err(|_| Error::LebCorrupt {
            what: what.to_string(),
            offset,
        })?;
        out.push(f64::from_le_bytes(word));
    }
    Ok(out)
}

/// Shift a two-part TT Julian day earlier by `delay_days`.
///
/// Positive delays move the epoch earlier; negative delays advance it.
/// [`TwoPartEpoch::add_offset`] retains the low-word contribution and owns
/// validation of operands and result. Its errors propagate unchanged.
pub fn retard_epoch(epoch: TwoPartEpoch, delay_days: f64) -> Result<TwoPartEpoch> {
    epoch.add_offset(-delay_days)
}

impl NutationHeader {
    /// Decode the header at `offset` (pure structural decode).
    pub fn decode(bytes: &[u8], offset: usize) -> Result<Self> {
        const WHAT: &str = "LEB nutation header";
        let chunk = slice_at(bytes, offset, NUTATION_HEADER_SIZE, WHAT)?;
        Ok(Self {
            jd_start: read_f64(chunk, 0, WHAT, offset)?,
            jd_end: read_f64(chunk, 8, WHAT, offset)?,
            interval_days: read_f64(chunk, 16, WHAT, offset)?,
            degree: read_u32(chunk, 24, WHAT, offset)?,
            components: read_u32(chunk, 28, WHAT, offset)?,
            segment_count: read_u32(chunk, 32, WHAT, offset)?,
        })
    }

    /// Coverage start (Julian Day, TT at evaluation).
    pub fn jd_start(&self) -> f64 {
        self.jd_start
    }
    /// Coverage end (Julian Day, TT at evaluation).
    pub fn jd_end(&self) -> f64 {
        self.jd_end
    }
    /// Segment length in days.
    pub fn interval_days(&self) -> f64 {
        self.interval_days
    }
    /// Chebyshev degree of the `(dpsi, deps)` channels.
    pub fn degree(&self) -> u32 {
        self.degree
    }
    /// Channel count (2: `dpsi`, `deps`, radians).
    pub fn components(&self) -> u32 {
        self.components
    }
    /// Segment count; zero means the table is absent.
    pub fn segment_count(&self) -> u32 {
        self.segment_count
    }
    /// False exactly when `segment_count == 0`.
    pub fn is_present(&self) -> bool {
        self.segment_count > 0
    }
}

impl DeltaTHeader {
    /// Decode the header at `offset` (pure structural decode).
    pub fn decode(bytes: &[u8], offset: usize) -> Result<Self> {
        const WHAT: &str = "LEB Delta-T header";
        let chunk = slice_at(bytes, offset, DELTA_T_HEADER_SIZE, WHAT)?;
        Ok(Self {
            n_entries: read_u32(chunk, 0, WHAT, offset)?,
        })
    }

    /// Declared `(jd, delta_t_days)` pair count.
    pub fn n_entries(&self) -> u32 {
        self.n_entries
    }

    /// Guard the declared pair count against the section payload size
    /// (the 8 header bytes excluded), via [`validate_entry_count`].
    pub fn validate_capacity(&self, section_size: u64, section_offset: usize) -> Result<()> {
        let available = section_size
            .checked_sub(DELTA_T_HEADER_SIZE as u64)
            .ok_or_else(|| Error::LebCorrupt {
                what: "LEB Delta-T table".to_string(),
                offset: section_offset,
            })?;
        validate_entry_count(
            u64::from(self.n_entries),
            DELTA_T_ENTRY_SIZE as u64,
            available,
            "LEB Delta-T table",
            section_offset,
        )
    }
}

impl StarEntry {
    /// Decode one entry at `offset` (pure structural decode).
    pub fn decode(bytes: &[u8], offset: usize) -> Result<Self> {
        const WHAT: &str = "LEB star catalog";
        let chunk = slice_at(bytes, offset, STAR_ENTRY_SIZE, WHAT)?;
        Ok(Self {
            star_id: read_i32(chunk, 0, WHAT, offset)?,
            ra_j2000_deg: read_f64(chunk, 4, WHAT, offset)?,
            dec_j2000_deg: read_f64(chunk, 12, WHAT, offset)?,
            pm_ra_deg_per_yr: read_f64(chunk, 20, WHAT, offset)?,
            pm_dec_deg_per_yr: read_f64(chunk, 28, WHAT, offset)?,
            parallax_arcsec: read_f64(chunk, 36, WHAT, offset)?,
            rv_km_per_s: read_f64(chunk, 44, WHAT, offset)?,
            magnitude: read_f64(chunk, 52, WHAT, offset)?,
        })
    }

    /// Catalog identifier.
    pub fn star_id(&self) -> i32 {
        self.star_id
    }
    /// Right ascension at J2000, degrees.
    pub fn ra_j2000_deg(&self) -> f64 {
        self.ra_j2000_deg
    }
    /// Declination at J2000, degrees.
    pub fn dec_j2000_deg(&self) -> f64 {
        self.dec_j2000_deg
    }
    /// Proper motion in RA, deg/yr (includes the cos(dec) factor).
    pub fn pm_ra_deg_per_yr(&self) -> f64 {
        self.pm_ra_deg_per_yr
    }
    /// Proper motion in Dec, deg/yr.
    pub fn pm_dec_deg_per_yr(&self) -> f64 {
        self.pm_dec_deg_per_yr
    }
    /// Parallax, arcsec.
    pub fn parallax_arcsec(&self) -> f64 {
        self.parallax_arcsec
    }
    /// Radial velocity, km/s.
    pub fn rv_km_per_s(&self) -> f64 {
        self.rv_km_per_s
    }
    /// Catalog magnitude.
    pub fn magnitude(&self) -> f64 {
        self.magnitude
    }
}

/// Number of complete 64-byte star records in a section payload.
/// Any trailing partial record is ignored.
pub fn star_count(section_size: u64) -> u64 {
    section_size / STAR_ENTRY_SIZE as u64
}

/// Compute `segment_count * (degree + 1) * components * 8` with checked
/// arithmetic. Overflow returns [`Error::LebCorrupt`]. This is a geometry
/// calculation, not a limit on total reader memory.
pub(super) fn coeff_bytes(
    segment_count: u32,
    degree: u32,
    components: u32,
    what: &str,
    offset: usize,
) -> Result<u64> {
    let deg1 = u64::from(degree)
        .checked_add(1)
        .ok_or_else(|| Error::LebCorrupt {
            what: what.to_string(),
            offset,
        })?;
    u64::from(segment_count)
        .checked_mul(deg1)
        .and_then(|v| v.checked_mul(u64::from(components)))
        .and_then(|v| v.checked_mul(8))
        .ok_or_else(|| Error::LebCorrupt {
            what: what.to_string(),
            offset,
        })
}

/// Check a chunk's declared size against its coefficient geometry.
/// Overflow and mismatch retain the caller's error label at offset zero.
pub(super) fn checked_chunk_uncompressed_size(
    chunk: &ChunkEntry,
    degree: u32,
    components: u32,
    what: &str,
) -> Result<u64> {
    let expected = coeff_bytes(chunk.segment_count, degree, components, what, 0)?;
    if chunk.uncompressed_size != expected {
        return Err(Error::LebCorrupt {
            what: what.to_string(),
            offset: 0,
        });
    }
    Ok(expected)
}

impl ChunkEntry {
    /// Decode one entry at `offset` (pure structural decode; per-chunk
    /// size/ownership checks happen at use).
    pub fn decode(bytes: &[u8], offset: usize) -> Result<Self> {
        const WHAT: &str = "LEB2 chunk index";
        let chunk = slice_at(bytes, offset, CHUNK_ENTRY_SIZE, WHAT)?;
        Ok(Self {
            jd_start: read_f64(chunk, 0, WHAT, offset)?,
            jd_end: read_f64(chunk, 8, WHAT, offset)?,
            blob_offset: read_u64(chunk, 16, WHAT, offset)?,
            compressed_size: read_u64(chunk, 24, WHAT, offset)?,
            uncompressed_size: read_u64(chunk, 32, WHAT, offset)?,
            segment_start: read_u32(chunk, 40, WHAT, offset)?,
            segment_count: read_u32(chunk, 44, WHAT, offset)?,
        })
    }

    /// Advisory chunk coverage start (Julian Day).
    pub fn jd_start(&self) -> f64 {
        self.jd_start
    }
    /// Advisory chunk coverage end (Julian Day).
    pub fn jd_end(&self) -> f64 {
        self.jd_end
    }
    /// Absolute byte offset of this chunk's zstd blob.
    pub fn blob_offset(&self) -> u64 {
        self.blob_offset
    }
    /// Stored compressed blob size in bytes.
    pub fn compressed_size(&self) -> u64 {
        self.compressed_size
    }
    /// Stored decoded size in bytes (checked against the derived
    /// geometry before any allocation).
    pub fn uncompressed_size(&self) -> u64 {
        self.uncompressed_size
    }
    /// First global segment index covered by this chunk.
    pub fn segment_start(&self) -> u32 {
        self.segment_start
    }
    /// Number of global segments covered by this chunk.
    pub fn segment_count(&self) -> u32 {
        self.segment_count
    }

    /// True when global segment `seg_idx` lies in
    /// `[segment_start, segment_start + segment_count)`, in `u64` so a
    /// saturated `u32` range can never wrap into a false positive.
    pub(super) fn covers(&self, seg_idx: u32) -> bool {
        let seg = u64::from(seg_idx);
        let start = u64::from(self.segment_start);
        seg >= start && seg - start < u64::from(self.segment_count)
    }
}

/// Find a candidate chunk by binary search over TT `jd_end` bounds.
///
/// For an index sorted by end time, the search advances past bounds with
/// `jd >= jd_end`. An empty index returns `None`. This is a hint only:
/// rounded chunk dates do not establish segment ownership. Confirm the
/// candidate with [`select_chunk`].
pub fn find_chunk(chunks: &[ChunkEntry], jd: f64) -> Option<usize> {
    if chunks.is_empty() {
        return None;
    }
    let mut lo = 0usize;
    let mut hi = chunks.len() - 1;
    while lo < hi {
        let mid = (lo + hi) / 2;
        if jd >= chunks[mid].jd_end {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    Some(lo)
}

/// Transpose `N` eight-byte words into eight byte lanes.
///
/// Output byte `lane * N + j` contains input byte `j * 8 + lane`.
/// All bit patterns are preserved without floating-point arithmetic. Empty
/// input is accepted; a length not divisible by eight returns
/// [`Error::LebCorrupt`].
pub fn shuffle_bytes(data: &[u8]) -> Result<Vec<u8>> {
    const WHAT: &str = "LEB2 shuffle frame";
    if data.len() % SHUFFLE_LANES != 0 {
        return Err(Error::LebCorrupt {
            what: WHAT.to_string(),
            offset: 0,
        });
    }
    let n = data.len() / SHUFFLE_LANES;
    let mut out = vec![0u8; data.len()];
    for (j, row) in data.chunks_exact(SHUFFLE_LANES).enumerate() {
        for lane in 0..SHUFFLE_LANES {
            out[lane * n + j] = row[lane];
        }
    }
    Ok(out)
}

/// Invert [`shuffle_bytes`], restoring consecutive eight-byte words.
/// All bit patterns are preserved. Empty input is accepted; a length not
/// divisible by eight returns [`Error::LebCorrupt`].
pub fn unshuffle_bytes(data: &[u8]) -> Result<Vec<u8>> {
    const WHAT: &str = "LEB2 shuffle frame";
    if data.len() % SHUFFLE_LANES != 0 {
        return Err(Error::LebCorrupt {
            what: WHAT.to_string(),
            offset: 0,
        });
    }
    let n = data.len() / SHUFFLE_LANES;
    let mut out = vec![0u8; data.len()];
    for j in 0..n {
        for lane in 0..SHUFFLE_LANES {
            out[j * SHUFFLE_LANES + lane] = data[lane * n + j];
        }
    }
    Ok(out)
}

/// Reorder coefficient words from segment-major to coefficient-major:
/// `[segment * coefficients_per_segment + k]` becomes
/// `[k * segment_count + segment]`.
///
/// Words are copied as bit patterns, including non-finite values and signed
/// zeros. Input length must equal
/// `segment_count * (degree + 1) * components * 8`; a mismatch or arithmetic
/// overflow returns [`Error::LebCorrupt`].
pub fn reorder_to_coeff_major(
    data: &[u8],
    segment_count: u32,
    degree: u32,
    components: u32,
) -> Result<Vec<u8>> {
    const WHAT: &str = "LEB2 coefficient layout";
    reorder_words(data, segment_count, degree, components, WHAT, true)
}

/// Shared bit-pattern transpose for both reorder directions.
pub(super) fn reorder_words(
    data: &[u8],
    segment_count: u32,
    degree: u32,
    components: u32,
    what: &str,
    to_coeff_major: bool,
) -> Result<Vec<u8>> {
    let expected = coeff_bytes(segment_count, degree, components, what, 0)?;
    if data.len() as u64 != expected {
        return Err(Error::LebCorrupt {
            what: what.to_string(),
            offset: 0,
        });
    }
    let segs = segment_count as usize;
    let deg1 = degree as usize + 1;
    let comps = components as usize;
    let n_coeffs = comps * deg1;
    let mut words = Vec::with_capacity(segs * n_coeffs);
    for chunk in data.chunks_exact(8) {
        let word: [u8; 8] = chunk.try_into().map_err(|_| Error::LebCorrupt {
            what: what.to_string(),
            offset: 0,
        })?;
        words.push(u64::from_le_bytes(word));
    }
    let mut out = vec![0u64; words.len()];
    for s in 0..segs {
        for k in 0..n_coeffs {
            let (from, to) = if to_coeff_major {
                (s * n_coeffs + k, k * segs + s)
            } else {
                (k * segs + s, s * n_coeffs + k)
            };
            out[to] = words[from];
        }
    }
    let mut bytes = Vec::with_capacity(data.len());
    for w in out {
        bytes.extend_from_slice(&w.to_le_bytes());
    }
    Ok(bytes)
}

/// Decompress an LEB2 payload to segment-major coefficient bytes.
///
/// The pipeline is zstd decoding, byte-lane unshuffle and coefficient
/// permutation. Stored binary64 words are preserved; mantissa bits already
/// discarded by a generator cannot be recovered.
///
/// `uncompressed_size` must equal
/// `segment_count * (degree + 1) * components * 8`. Overflow, zstd failures
/// and decoded-length mismatches return [`Error::LebCorrupt`]. At most one
/// extra output byte is read to detect overlong payloads. This is not a
/// total allocation or decompressor-workspace cap; indexed readers apply
/// their separate resource and window limits.
pub fn decompress_body(
    compressed: &[u8],
    uncompressed_size: u64,
    segment_count: u32,
    degree: u32,
    components: u32,
) -> Result<Vec<u8>> {
    decompress_body_limited(
        compressed,
        uncompressed_size,
        segment_count,
        degree,
        components,
        None,
    )
}

/// The same lossless decode with an optional positional-reader window policy.
pub(super) fn decompress_body_limited(
    compressed: &[u8],
    uncompressed_size: u64,
    segment_count: u32,
    degree: u32,
    components: u32,
    window_log_max: Option<u32>,
) -> Result<Vec<u8>> {
    let decoded = decompress_shuffled_limited(
        compressed,
        uncompressed_size,
        segment_count,
        degree,
        components,
        window_log_max,
    )?;
    let mut bytes = Vec::with_capacity(decoded.len());
    for_each_segment_word(&decoded, segment_count, |word| {
        bytes.extend_from_slice(&word)
    });
    Ok(bytes)
}

/// Decode directly to the segment-major binary64 owner consumed by readers.
/// Byte-lane unshuffle, coefficient permutation and little-endian word
/// reconstruction share one pass; there are no intermediate whole-frame
/// unshuffled/reordered buffers. Every binary64 bit pattern is retained.
pub(super) fn decompress_body_words_limited(
    compressed: &[u8],
    uncompressed_size: u64,
    segment_count: u32,
    degree: u32,
    components: u32,
    window_log_max: Option<u32>,
) -> Result<Vec<f64>> {
    let decoded = decompress_shuffled_limited(
        compressed,
        uncompressed_size,
        segment_count,
        degree,
        components,
        window_log_max,
    )?;
    let mut words = Vec::with_capacity(decoded.len() / SHUFFLE_LANES);
    for_each_segment_word(&decoded, segment_count, |word| {
        words.push(f64::from_le_bytes(word))
    });
    Ok(words)
}

/// Validate the same complete frame, then reconstruct only complete selected
/// rows. The retained window is ordinary segment-major binary64 storage;
/// omitted rows still participate in frame length/checksum validation.
pub(super) fn decompress_body_window_words_limited(
    compressed: &[u8],
    uncompressed_size: u64,
    segment_count: u32,
    degree: u32,
    components: u32,
    window: std::ops::Range<u32>,
    window_log_max: Option<u32>,
) -> Result<Vec<f64>> {
    let decoded = decompress_shuffled_limited(
        compressed,
        uncompressed_size,
        segment_count,
        degree,
        components,
        window_log_max,
    )?;
    if window.start > window.end || window.end > segment_count {
        return Err(Error::LebCorrupt {
            what: "LEB2 coefficient window".into(),
            offset: 0,
        });
    }
    let word_count = usize::try_from(
        coeff_bytes(
            window.end - window.start,
            degree,
            components,
            "LEB2 coefficient window",
            0,
        )? / SHUFFLE_LANES as u64,
    )
    .map_err(|_| Error::LebCorrupt {
        what: "LEB2 coefficient window".into(),
        offset: 0,
    })?;
    let mut words = Vec::with_capacity(word_count);
    for_each_segment_word_range(&decoded, segment_count, window, |word| {
        words.push(f64::from_le_bytes(word))
    });
    Ok(words)
}

/// Visit validated shuffled words in their final segment-major order.
/// The shared decoder has proved exact geometry before this helper is called:
/// `decoded.len() = segments * coefficients_per_segment * 8`. Thus every
/// lane index is in bounds, including zero-segment/zero-component layouts.
/// Byte arrays carry payload bits without any floating-point arithmetic.
fn for_each_segment_word(
    decoded: &[u8],
    segment_count: u32,
    consume: impl FnMut([u8; SHUFFLE_LANES]),
) {
    for_each_segment_word_range(decoded, segment_count, 0..segment_count, consume);
}

fn for_each_segment_word_range(
    decoded: &[u8],
    segment_count: u32,
    window: std::ops::Range<u32>,
    mut consume: impl FnMut([u8; SHUFFLE_LANES]),
) {
    let words = decoded.len() / SHUFFLE_LANES;
    let segments = segment_count as usize;
    if words == 0 {
        return;
    }
    let coefficients_per_segment = words / segments;
    for segment in window.start as usize..window.end as usize {
        for coefficient in 0..coefficients_per_segment {
            let source_word = coefficient * segments + segment;
            let word = std::array::from_fn(|lane| decoded[lane * words + source_word]);
            consume(word);
        }
    }
}

/// Decode only the shuffled coefficient-major payload. Geometry, streaming
/// output and positional window admission belong to this single shared path.
fn decompress_shuffled_limited(
    compressed: &[u8],
    uncompressed_size: u64,
    segment_count: u32,
    degree: u32,
    components: u32,
    window_log_max: Option<u32>,
) -> Result<Vec<u8>> {
    const WHAT: &str = "LEB2 coefficient blob";
    let expected = coeff_bytes(segment_count, degree, components, WHAT, 0)?;
    if uncompressed_size != expected {
        return Err(Error::LebCorrupt {
            what: WHAT.to_string(),
            offset: 0,
        });
    }
    let expected_len = usize::try_from(expected).map_err(|_| Error::LebCorrupt {
        what: WHAT.to_string(),
        offset: 0,
    })?;
    let limit = expected_len
        .checked_add(1)
        .ok_or_else(|| Error::LebCorrupt {
            what: WHAT.to_string(),
            offset: 0,
        })?;
    // Bounded decode: at most `expected + 1` bytes are ever materialized,
    // so a crafted frame cannot force an unbounded allocation. One byte
    // over the budget is decoded deliberately to distinguish "exact" from
    // "over-long" instead of truncating silently.
    if let Some(window) = window_log_max {
        access::check_zstd_window(compressed, window)?;
    }
    // A complete single frame can decode directly into its final shuffled
    // buffer, without the streaming decoder's history/output copies. Owned
    // callers admit at most 1 MiB eagerly even when geometry is forged, and
    // require a similarly small window. Larger windows retain the streaming
    // decoder's ordinary default-window refusal and error classification;
    // indexed callers already reserved their explicitly bounded frame/workspace.
    // Keep concatenated/skippable sequences on the ordinary streaming path so
    // every later frame still receives the decoder's window enforcement.
    const DIRECT_WINDOW_LOG: u32 = 20;
    let direct_frame = window_log_max.is_some()
        || (limit <= 1 << DIRECT_WINDOW_LOG
            && access::check_zstd_window(compressed, DIRECT_WINDOW_LOG).is_ok());
    if direct_frame
        && zstd::zstd_safe::find_frame_compressed_size(compressed)
            .is_ok_and(|frame_size| frame_size == compressed.len())
    {
        let mut decoder = zstd::bulk::Decompressor::new().map_err(|_| Error::LebCorrupt {
            what: "LEB2 zstd frame".into(),
            offset: 0,
        })?;
        if let Some(window) = window_log_max {
            decoder
                .window_log_max(window)
                .map_err(|_| Error::LebCorrupt {
                    what: "LEB2 zstd window".into(),
                    offset: 0,
                })?;
        }
        let mut decoded = Vec::with_capacity(limit);
        decoder
            .decompress_to_buffer(compressed, &mut decoded)
            .map_err(|_| Error::LebCorrupt {
                what: "LEB2 zstd frame".into(),
                offset: 0,
            })?;
        if decoded.len() != expected_len {
            return Err(Error::LebCorrupt {
                what: "LEB2 zstd frame".into(),
                offset: 0,
            });
        }
        return Ok(decoded);
    }
    let mut decoder = zstd::Decoder::new(compressed).map_err(|_| Error::LebCorrupt {
        what: "LEB2 zstd frame".to_string(),
        offset: 0,
    })?;
    if let Some(window) = window_log_max {
        decoder
            .window_log_max(window)
            .map_err(|_| Error::LebCorrupt {
                what: "LEB2 zstd window".into(),
                offset: 0,
            })?;
    }
    // On this streaming path, only an admitted positional frame may reserve
    // its complete output. Larger unrestricted slice inputs keep incremental
    // growth, so forged geometry cannot trigger a huge eager allocation.
    let mut decoded = if window_log_max.is_some() {
        Vec::with_capacity(limit)
    } else {
        Vec::new()
    };
    {
        use std::io::Read as _;
        let mut capped = (&mut decoder).take(limit as u64);
        capped
            .read_to_end(&mut decoded)
            .map_err(|_| Error::LebCorrupt {
                what: "LEB2 zstd frame".to_string(),
                offset: 0,
            })?;
    }
    if decoded.len() != expected_len {
        return Err(Error::LebCorrupt {
            what: "LEB2 zstd frame".to_string(),
            offset: 0,
        });
    }
    Ok(decoded)
}

/// Decompress an LEB2 v2 chunk from its already-sliced blob.
///
/// The stored size must match the chunk segment count, degree and component
/// count, or [`Error::LebCorrupt`] is returned. Further errors follow
/// [`decompress_body`]. Chunk selection, file offsets and retention belong
/// to the caller.
pub fn decompress_chunk(
    compressed: &[u8],
    chunk: &ChunkEntry,
    degree: u32,
    components: u32,
    body_id: i32,
    chunk_idx: usize,
) -> Result<Vec<u8>> {
    let what = format!("LEB2 chunk (body {body_id}, chunk {chunk_idx})");
    checked_chunk_uncompressed_size(chunk, degree, components, &what)?;
    decompress_body(
        compressed,
        chunk.uncompressed_size,
        chunk.segment_count,
        degree,
        components,
    )
}

/// Decompress a monolithic LEB2 payload using its body-entry geometry.
///
/// A missing stored size, including an LEB1 entry, or a size inconsistent
/// with body geometry returns [`Error::LebCorrupt`]. Further errors follow
/// [`decompress_body`]. The caller slices the blob and manages retention.
pub fn decompress_monolith(compressed: &[u8], entry: &BodyEntry) -> Result<Vec<u8>> {
    let what = format!("LEB2 body entry {}", entry.body_id());
    let stored = entry.uncompressed_size().ok_or_else(|| Error::LebCorrupt {
        what: what.clone(),
        offset: 0,
    })?;
    let expected = entry.expected_uncompressed_size()?;
    if stored != expected {
        return Err(Error::LebCorrupt { what, offset: 0 });
    }
    decompress_body(
        compressed,
        stored,
        entry.segment_count(),
        entry.degree(),
        entry.components(),
    )
}

/// Known group suffix of an extension-stripped file name.
/// For example, `base_core.leb2` yields `Some("core")`. The prefix is not
/// validated as a tier name; unknown suffixes return `None`.
pub fn file_group_suffix(file_name: &str) -> Option<&str> {
    // The release `rsplit(".", 1)[0]` keeps the full name for dotless
    // inputs; `rsplit_once` returns `None` there, so fall back to the
    // whole name in that case.
    let stem = if file_name.contains('.') {
        file_name
            .rsplit_once('.')
            .map(|(s, _)| s)
            .unwrap_or(file_name)
    } else {
        file_name
    };
    let last = stem.rsplit('_').next().unwrap_or(stem);
    GROUP_SUFFIXES.iter().find(|g| **g == last).copied()
}

/// Whether the file name carries the retired `uranians` group suffix.
/// This classifies a name without opening or attaching the file.
pub fn is_retired_group_file(file_name: &str) -> bool {
    file_group_suffix(file_name).is_some_and(|suffix| RETIRED_GROUP_SUFFIXES.contains(&suffix))
}

/// Prefix preceding a known group suffix, such as `base` in
/// `base_core.leb2`. Companion candidates share this prefix and have a known
/// group suffix. Other names, including `ephemeris_base.leb`, return `None`;
/// this function does not inspect their inventories.
pub fn group_prefix(file_name: &str) -> Option<String> {
    file_group_suffix(file_name)?;
    let stem = if file_name.contains('.') {
        file_name
            .rsplit_once('.')
            .map(|(s, _)| s)
            .unwrap_or(file_name)
    } else {
        file_name
    };
    let mut parts: Vec<&str> = stem.split('_').collect();
    parts.pop();
    Some(parts.join("_"))
}

/// System-barycenter name for a body needing a center-of-body offset.
/// Returns `None` for identifiers absent from the configured mapping.
pub fn system_barycenter_name(body_id: i32) -> Option<&'static str> {
    match body_id {
        5 => Some("jupiter barycenter"),
        6 => Some("saturn barycenter"),
        7 => Some("uranus barycenter"),
        8 => Some("neptune barycenter"),
        9 => Some("pluto barycenter"),
        _ => None,
    }
}

/// True for stored system-barycenter channels
/// (`COORD_ICRS_BARY_SYSTEM`); planet-center (`COORD_ICRS_BARY`) and
/// ecliptic channels need no center-of-body correction.
pub fn is_system_barycenter_channel(coord_type: u32) -> bool {
    coord_type == COORD_ICRS_BARY_SYSTEM
}

/// Add a barycenter-to-center offset component by component.
/// Both vectors are ICRS Cartesian positions in AU. This helper does not
/// validate finite values or source identity.
pub fn apply_center_offset(barycenter_au: [f64; 3], offset_au: [f64; 3]) -> [f64; 3] {
    [
        barycenter_au[0] + offset_au[0],
        barycenter_au[1] + offset_au[1],
        barycenter_au[2] + offset_au[2],
    ]
}

/// Keep only finite offsets: a non-finite provider result is treated as
/// unavailable.
pub(super) fn finite_offset(offset: [f64; 3]) -> Option<[f64; 3]> {
    if offset.iter().all(|v| v.is_finite()) {
        Some(offset)
    } else {
        None
    }
}

/// Whether `body_id` belongs to the analytical fictitious-body range `40..=58`.
pub fn is_fictitious_body(body_id: i32) -> bool {
    (FICTITIOUS_BODY_MIN..=FICTITIOUS_BODY_MAX).contains(&body_id)
}

/// Whether `body_id` is subject to sealed core-coverage policy:
/// `0..=12`, Earth (14), and the interpolated apsides (21, 22).
/// This classifies an identifier, not the presence of a stored channel.
pub fn is_core_body(body_id: i32) -> bool {
    (0..=12).contains(&body_id) || matches!(body_id, 14 | 21 | 22)
}

/// Whether a core body requires its own stored state channel:
/// `0..=9` and Earth (14). See [`sealed_outcome`] for missing-channel policy.
pub fn is_stored_state_body(body_id: i32) -> bool {
    (0..=9).contains(&body_id) || body_id == 14
}

/// Whether a body uses a packaged fitted-model window: the interpolated
/// apsides (21, 22), independently of installed LEB groups.
pub fn has_fitted_model_window(body_id: i32) -> bool {
    matches!(body_id, 21 | 22)
}

/// Inclusive comparison `start <= jd && jd <= end`.
/// Non-finite `jd` cannot satisfy finite bounds. No bounds validation or
/// time-scale conversion is performed.
pub fn interval_contains(start: f64, end: f64, jd: f64) -> bool {
    start <= jd && jd <= end
}

impl SealedOutcome {
    /// The typed error this disposition must surface, if any. `Covered`,
    /// `LocalModelFallthrough`, and `NoLebChannel` carry no error: the
    /// caller evaluates or routes onward.
    pub fn refusal_error(&self) -> Option<Error> {
        match *self {
            SealedOutcome::RefusedRange {
                body_id,
                jd,
                start,
                end,
            } => Some(Error::LebOutOfRange {
                body_id,
                jd,
                start,
                end,
            }),
            SealedOutcome::MissingCoreChannel { body_id } => {
                Some(Error::LebBodyNotFound { body_id })
            }
            _ => None,
        }
    }
}

/// Classify sealed coverage from the requested body's stored interval
/// and the Sun's support interval, in TT Julian days.
///
/// Packaged fitted-model windows are checked first. Body coverage is
/// inclusive. Uncovered core bodies are refused; missing stored-state core
/// channels have a distinct missing-channel outcome. For non-core bodies,
/// Sun support can permit the declared local-model route. Pass `None` when
/// no stored interval exists.
///
/// This performs no file reads, model evaluation or independent validation
/// of the supplied coverage metadata.
pub fn sealed_outcome(
    body_id: i32,
    jd: f64,
    body_cov: Option<(f64, f64)>,
    sun_cov: Option<(f64, f64)>,
) -> SealedOutcome {
    if has_fitted_model_window(body_id) {
        if let Some((start, end)) = interpolated_apse_model_window(body_id) {
            if !interval_contains(start, end, jd) {
                return SealedOutcome::RefusedRange {
                    body_id,
                    jd,
                    start,
                    end,
                };
            }
        }
    }
    match body_cov {
        Some((start, end)) => {
            if interval_contains(start, end, jd) {
                SealedOutcome::Covered
            } else if is_core_body(body_id) {
                SealedOutcome::RefusedRange {
                    body_id,
                    jd,
                    start,
                    end,
                }
            } else {
                match sun_cov {
                    Some((s_start, s_end)) if interval_contains(s_start, s_end, jd) => {
                        SealedOutcome::LocalModelFallthrough
                    }
                    _ => SealedOutcome::RefusedRange {
                        body_id,
                        jd,
                        start,
                        end,
                    },
                }
            }
        }
        None => {
            if is_core_body(body_id) && is_stored_state_body(body_id) {
                SealedOutcome::MissingCoreChannel { body_id }
            } else if is_core_body(body_id) {
                SealedOutcome::NoLebChannel { body_id }
            } else {
                match sun_cov {
                    Some((s_start, s_end)) if interval_contains(s_start, s_end, jd) => {
                        SealedOutcome::LocalModelFallthrough
                    }
                    Some((s_start, s_end)) => SealedOutcome::RefusedRange {
                        body_id,
                        jd,
                        start: s_start,
                        end: s_end,
                    },
                    None => SealedOutcome::NoLebChannel { body_id },
                }
            }
        }
    }
}

impl PrecisionClass {
    /// String representation of this provenance class.
    pub fn as_str(self) -> &'static str {
        match self {
            PrecisionClass::Ephemeris => "ephemeris",
            PrecisionClass::Analytical => "analytical",
            PrecisionClass::NumericalModel => "numerical-model",
            PrecisionClass::Mixed => "mixed",
            PrecisionClass::UnverifiedLocal => "unverified-local",
        }
    }
}

/// True for the N-body extended-tier trajectory ids (see
/// [`NBODY_MODEL_BODY_IDS`]).
pub fn is_nbody_model_body(body_id: i32) -> bool {
    NBODY_MODEL_BODY_IDS.contains(&body_id)
}

/// Classify provenance from the body identifier, file-name tokens and
/// caller-supplied `reviewed` flag.
///
/// Unreviewed inputs are `UnverifiedLocal`; analytical identifiers are
/// `Analytical`. An N-body identifier without a file is `Mixed`, or
/// `NumericalModel` when the file stem contains `extended`. Other reviewed
/// inputs are `Ephemeris`. No file verification or error bound is computed.
pub fn precision_class(body_id: i32, file_name: Option<&str>, reviewed: bool) -> PrecisionClass {
    if !reviewed {
        return PrecisionClass::UnverifiedLocal;
    }
    if ANALYTICAL_BODY_IDS.contains(&body_id) || is_fictitious_body(body_id) {
        return PrecisionClass::Analytical;
    }
    if is_nbody_model_body(body_id) {
        let Some(name) = file_name else {
            return PrecisionClass::Mixed;
        };
        if file_stem_tokens(name).contains(&"extended") {
            return PrecisionClass::NumericalModel;
        }
    }
    PrecisionClass::Ephemeris
}

/// Stem tokens of a file path: the base name after the last `/` with the
/// last `.suffix` stripped, split on `'_'`.
pub(super) fn file_stem_tokens(path: &str) -> Vec<&str> {
    let base = path.rsplit('/').next().unwrap_or(path);
    let stem = if base.contains('.') {
        base.rsplit_once('.').map(|(s, _)| s).unwrap_or(base)
    } else {
        base
    };
    stem.split('_').collect()
}

/// Canonical LEB group name inferred from a manifest-style file path
/// (`Some("core")` for `.../base_core.leb2`): `None` for a missing path or a stem without
/// `'_'`, otherwise the text after the last `'_'`. Unlike
/// [`file_group_suffix`] this reports any trailing token, not only the
/// known distribution suffixes.
pub fn group_from_path(path: Option<&str>) -> Option<String> {
    let path = path?;
    let base = path.rsplit('/').next().unwrap_or(path);
    let stem = if base.contains('.') {
        base.rsplit_once('.').map(|(s, _)| s).unwrap_or(base)
    } else {
        base
    };
    if !stem.contains('_') {
        return None;
    }
    stem.rsplit('_').next().map(str::to_string)
}

impl<'s> ServingFile<'s> {
    /// Unreviewed file with no path label (the tiered-union case).
    pub fn unlabeled() -> Self {
        Self {
            path: None,
            reviewed: false,
        }
    }
}

impl BodyRange {
    /// Body identifier in the public numbering.
    pub fn body_id(&self) -> i32 {
        self.body_id
    }
    /// Start of the stored interval (Julian Day).
    pub fn jd_start(&self) -> f64 {
        self.jd_start
    }
    /// End of the stored interval (Julian Day).
    pub fn jd_end(&self) -> f64 {
        self.jd_end
    }
}

impl FileInfo {
    /// Base file name.
    pub fn name(&self) -> &str {
        &self.name
    }
    /// Caller-supplied container path or image label.
    pub fn path(&self) -> &str {
        &self.path
    }
    /// Canonical group inferred from the path (see [`group_from_path`]).
    pub fn group(&self) -> Option<&str> {
        self.group.as_deref()
    }
    /// Caller-observed file length in bytes, if known.
    pub fn size_bytes(&self) -> Option<u64> {
        self.size_bytes
    }
    /// The caller-supplied reviewed flag; this accessor performs no verification.
    pub fn reviewed(&self) -> bool {
        self.reviewed
    }
    /// Per-body stored intervals in id order.
    pub fn bodies(&self) -> &[BodyRange] {
        &self.bodies
    }
    /// Number of indexed bodies.
    pub fn body_count(&self) -> usize {
        self.bodies.len()
    }
}

/// Decode a container index and return its inventory metadata.
/// `path` labels the row; `reviewed` is supplied by the caller, not verified
/// here. Errors propagate from [`FileIndex::decode`]. Coefficient payloads
/// are not evaluated and scientific accuracy is not checked.
pub fn inspect_bytes(path: &str, bytes: &[u8], reviewed: bool) -> Result<FileInfo> {
    let index = FileIndex::decode(bytes)?;
    let size = u64::try_from(bytes.len()).unwrap_or(u64::MAX);
    Ok(inspect_index(path, Some(size), reviewed, &index))
}

/// Inspect a regular file's header/directory/body metadata without mapping or
/// loading its coefficient payloads. Safe for generic paths; the caller must
/// keep the file stable during inspection for a coherent row. Always unverified:
/// no hash, immutable snapshot, or scientific acceptance is implied.
pub fn inspect_file(path: &std::path::Path) -> Result<FileInfo> {
    let mut file = storage::open_regular_file(path).map_err(|error| {
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
    let size = file.metadata().ok().map(|m| m.len());
    let index = FileIndex::read_from_file(&mut file)?;
    Ok(inspect_index(&path.to_string_lossy(), size, false, &index))
}

impl LebBodyVerify {
    /// Indexed body this row covers.
    pub fn body_id(&self) -> i32 {
        self.body_id
    }
    /// Start of the stored interval (Julian Day).
    pub fn jd_start(&self) -> f64 {
        self.jd_start
    }
    /// End of the stored interval (Julian Day).
    pub fn jd_end(&self) -> f64 {
        self.jd_end
    }
    /// Number of evaluation probes executed for this body (0 on an
    /// index-level refusal, before any probe ran).
    pub fn probes(&self) -> u32 {
        self.probes
    }
    /// First probe failure detail, or `None` on a passing body.
    pub fn error(&self) -> Option<&str> {
        self.error.as_deref()
    }
    /// True only when every probe for this body evaluated.
    pub fn passed(&self) -> bool {
        self.error.is_none()
    }
}

impl LebVerifyReport {
    /// Container flavor the report covers.
    pub fn format(&self) -> LebFormat {
        self.format
    }
    /// Per-body rows in id order.
    pub fn bodies(&self) -> &[LebBodyVerify] {
        &self.bodies
    }
    /// Number of indexed bodies.
    pub fn body_count(&self) -> usize {
        self.bodies.len()
    }
    /// Ids of the bodies with a probe failure, in id order.
    pub fn failures(&self) -> Vec<i32> {
        self.bodies
            .iter()
            .filter(|row| !row.passed())
            .map(LebBodyVerify::body_id)
            .collect()
    }
    /// Whether the report contains at least one body and every body passed.
    /// An empty report returns `false`.
    pub fn passed(&self) -> bool {
        !self.bodies.is_empty() && self.bodies.iter().all(LebBodyVerify::passed)
    }
}

/// Deterministic probe epochs for one stored interval: both endpoints
/// (segment-clamp boundaries, where segment-selection seam ownership applies) plus
/// the midpoint. Degenerate (zero-width) intervals probe a single epoch.
pub(super) fn verify_probe_epochs(jd_start: f64, jd_end: f64) -> Vec<f64> {
    if jd_end <= jd_start {
        return vec![jd_start];
    }
    let mid = jd_start + (jd_end - jd_start) / 2.0;
    if mid <= jd_start || mid >= jd_end {
        return vec![jd_start, jd_end];
    }
    vec![jd_start, mid, jd_end]
}

/// Decode a container index and probe each body at its stored TT interval
/// endpoints and midpoint.
///
/// Index-decoding failures return an error. Once the index is admitted, a
/// failing probe is recorded for that body and other bodies are still
/// checked. These are finite structural/evaluation probes, not a proof of
/// every coefficient or astronomical accuracy.
pub fn verify_container(bytes: &[u8]) -> Result<LebVerifyReport> {
    let index = FileIndex::decode(bytes)?;
    let format = index.header().format();
    let mut ids: Vec<i32> = index.bodies().iter().map(BodyEntry::body_id).collect();
    ids.sort_unstable();
    let mut bodies = Vec::with_capacity(ids.len());
    for body_id in ids {
        let entry = index.find_body(body_id)?;
        let (jd_start, jd_end) = (entry.jd_start(), entry.jd_end());
        let mut row = LebBodyVerify {
            body_id,
            jd_start,
            jd_end,
            probes: 0,
            error: None,
        };
        if entry.segment_count() == 0
            || !entry.interval_days().is_finite()
            || entry.interval_days() <= 0.0
        {
            row.error = Some(format!(
                "body {body_id}: index entry has segment_count={} interval_days={}",
                entry.segment_count(),
                entry.interval_days()
            ));
            bodies.push(row);
            continue;
        }
        for jd in verify_probe_epochs(jd_start, jd_end) {
            match eval_verify_probe(&index, bytes, body_id, jd) {
                Ok(_) => row.probes += 1,
                Err(error) => {
                    row.error = Some(format!("body {body_id} probe jd={jd}: {error}"));
                    break;
                }
            }
        }
        bodies.push(row);
    }
    Ok(LebVerifyReport { format, bodies })
}

impl Leb2BodyComparison {
    /// Compared body id.
    pub fn body_id(&self) -> i32 {
        self.body_id
    }
    /// Number of shared-coverage samples evaluated.
    pub fn samples(&self) -> u32 {
        self.samples
    }
    /// Worst absolute position/velocity component difference observed.
    pub fn max_abs_err(&self) -> f64 {
        self.max_abs_err
    }
}

impl Leb2BodySkipped {
    /// Skipped body id.
    pub fn body_id(&self) -> i32 {
        self.body_id
    }
    /// Why the body was not compared.
    pub fn reason(&self) -> &str {
        &self.reason
    }
}

impl Leb2CompareReport {
    /// Per-body comparison rows in id order.
    pub fn compared(&self) -> &[Leb2BodyComparison] {
        &self.compared
    }
    /// Per-body skip rows in id order.
    pub fn skipped(&self) -> &[Leb2BodySkipped] {
        &self.skipped
    }
    /// Number of bodies compared.
    pub fn compared_count(&self) -> usize {
        self.compared.len()
    }
    /// Worst per-body error over all compared bodies, or `None` when no
    /// body was compared.
    pub fn overall_max_abs_err(&self) -> Option<f64> {
        self.compared
            .iter()
            .map(Leb2BodyComparison::max_abs_err)
            .fold(None, |worst: Option<f64>, err| {
                Some(worst.map_or(err, |w: f64| w.max(err)))
            })
    }
}

/// Deterministic uniform sample grid over `[lo, hi]`, endpoints
/// included. A single sample probes the midpoint; otherwise the grid is
/// `lo + (hi - lo) * i / (n - 1)`.
pub(super) fn compare_sample_grid(lo: f64, hi: f64, n: usize) -> Vec<f64> {
    if n <= 1 || hi <= lo {
        return vec![lo + (hi - lo) / 2.0];
    }
    (0..n)
        .map(|i| lo + (hi - lo) * (i as f64) / ((n - 1) as f64))
        .collect()
}

/// Compare two containers body by body over shared TT coverage.
///
/// Each image uses its own raw or compressed evaluation route; either may
/// use any admitted LEB format. `samples_per_body` is clamped to
/// `1..=LEB2_COMPARE_MAX_SAMPLES`. Absent bodies, unequal component counts
/// and disjoint coverage are recorded as skips. Decode/evaluation errors
/// propagate.
///
/// Stored values and TT-day rates are compared without coordinate or unit
/// conversion; the caller must choose compatible conventions. Reported
/// maxima combine unscaled component differences. No tolerance or
/// acceptance verdict is applied.
pub fn compare_leb2_against_leb1(
    candidate: &[u8],
    reference: &[u8],
    samples_per_body: usize,
) -> Result<Leb2CompareReport> {
    let cand_index = FileIndex::decode(candidate)?;
    let ref_index = FileIndex::decode(reference)?;
    let n = samples_per_body.clamp(1, LEB2_COMPARE_MAX_SAMPLES);
    let mut cand_ids: Vec<i32> = cand_index.bodies().iter().map(BodyEntry::body_id).collect();
    cand_ids.sort_unstable();
    let mut ref_ids: Vec<i32> = ref_index.bodies().iter().map(BodyEntry::body_id).collect();
    ref_ids.sort_unstable();
    let mut compared = Vec::new();
    let mut skipped = Vec::new();
    for body_id in cand_ids.iter().chain(ref_ids.iter()) {
        let body_id = *body_id;
        if compared
            .iter()
            .any(|row: &Leb2BodyComparison| row.body_id == body_id)
            || skipped
                .iter()
                .any(|row: &Leb2BodySkipped| row.body_id == body_id)
        {
            continue;
        }
        let (cand_entry, ref_entry) =
            match (cand_index.find_body(body_id), ref_index.find_body(body_id)) {
                (Ok(cand), Ok(reference)) => (cand, reference),
                (Err(_), _) => {
                    skipped.push(Leb2BodySkipped {
                        body_id,
                        reason: "absent from candidate".to_string(),
                    });
                    continue;
                }
                (_, Err(_)) => {
                    skipped.push(Leb2BodySkipped {
                        body_id,
                        reason: "absent from reference".to_string(),
                    });
                    continue;
                }
            };
        if cand_entry.components() != ref_entry.components() {
            skipped.push(Leb2BodySkipped {
                body_id,
                reason: format!(
                    "component count differs (candidate {} vs reference {})",
                    cand_entry.components(),
                    ref_entry.components()
                ),
            });
            continue;
        }
        let lo = cand_entry.jd_start().max(ref_entry.jd_start());
        let hi = cand_entry.jd_end().min(ref_entry.jd_end());
        if hi < lo {
            skipped.push(Leb2BodySkipped {
                body_id,
                reason: "no shared coverage".to_string(),
            });
            continue;
        }
        let mut worst = 0.0f64;
        let mut samples = 0u32;
        for jd in compare_sample_grid(lo, hi, n) {
            let (cand_pos, cand_vel) = eval_verify_probe(&cand_index, candidate, body_id, jd)?;
            let (ref_pos, ref_vel) = eval_verify_probe(&ref_index, reference, body_id, jd)?;
            for (cand, reference) in cand_pos
                .iter()
                .chain(cand_vel.iter())
                .zip(ref_pos.iter().chain(ref_vel.iter()))
            {
                worst = worst.max((cand - reference).abs());
            }
            samples += 1;
        }
        compared.push(Leb2BodyComparison {
            body_id,
            samples,
            max_abs_err: worst,
        });
    }
    compared.sort_by_key(Leb2BodyComparison::body_id);
    skipped.sort_by_key(Leb2BodySkipped::body_id);
    Ok(Leb2CompareReport { compared, skipped })
}

impl RuntimeDataRequirement {
    /// Manifest file name (`{tier}_{group}.leb2`).
    pub fn name(&self) -> &str {
        &self.name
    }
    /// Artifact kind: always `"leb2"`.
    pub fn kind(&self) -> &'static str {
        "leb2"
    }
    /// Distribution group (`LEB2_GROUPS` entry).
    pub fn group(&self) -> &str {
        &self.group
    }
    /// Expected path under `data_dir/leb`; a relative root stays relative.
    pub fn path(&self) -> &str {
        &self.path
    }
    /// Reviewed manifest SHA-256 pin.
    pub fn sha256(&self) -> &str {
        &self.sha256
    }
}

/// Build the required file list for sealed operation at the active tier.
///
/// `groups` supplies partitions; `pins` supplies names and SHA-256 strings.
/// A missing or empty digest returns [`Error::DataManifestMissing`] instead
/// of a partial list. The data root is used as supplied. Output follows
/// cumulative [`eligible_tiers`], then the supplied group order. This
/// constructs requirements; it does not inspect files or verify hashes.
pub fn runtime_data_requirements(
    active: LebTier,
    groups: &[&str],
    pins: &[(&str, &str)],
    data_dir: &str,
) -> Result<Vec<RuntimeDataRequirement>> {
    let mut requirements = Vec::new();
    for tier in eligible_tiers(active) {
        for group in groups {
            let name = format!("{}_{group}.leb2", tier.name());
            let sha256 = pins
                .iter()
                .find_map(|(file, sha)| (*file == name).then_some(*sha))
                .filter(|sha| !sha.is_empty());
            let Some(sha256) = sha256 else {
                return Err(Error::DataManifestMissing { name });
            };
            requirements.push(RuntimeDataRequirement {
                path: format!("{data_dir}/leb/{name}"),
                name,
                group: (*group).to_string(),
                sha256: sha256.to_string(),
            });
        }
    }
    Ok(requirements)
}

impl AssetStatus {
    /// Manifest file name (`{tier}_{group}.leb2`).
    pub fn name(&self) -> &str {
        &self.name
    }
    /// Distribution group (`LEB2_GROUPS` entry).
    pub fn group(&self) -> &str {
        &self.group
    }
    /// Expected path under `data_dir/leb`; a relative root stays relative.
    pub fn path(&self) -> &str {
        &self.path
    }
    /// Validation outcome (see [`AssetState`]).
    pub fn state(&self) -> &AssetState {
        &self.state
    }
    /// True only for present, header-parseable assets.
    pub fn is_ready(&self) -> bool {
        matches!(self.state, AssetState::Ready { .. })
    }
}

/// Whether a core file name has a listed sibling companion.
///
/// The prefix is the stem before its final underscore. A sibling must be a
/// different `.leb` or `.leb2` name with that prefix, present in both
/// `manifest_names` and `dir_files`. A missing or empty prefix returns
/// `false`. This checks names only: no files are opened or digests verified.
pub fn has_pinned_siblings(core_name: &str, dir_files: &[&str], manifest_names: &[&str]) -> bool {
    let base = core_name.rsplit('/').next().unwrap_or(core_name);
    let no_ext = base.rsplit_once('.').map(|(stem, _)| stem).unwrap_or(base);
    let Some((stem_prefix, _)) = no_ext.rsplit_once('_') else {
        return false;
    };
    if stem_prefix.is_empty() {
        return false;
    }
    let want = format!("{stem_prefix}_");
    manifest_names.iter().any(|manifest| {
        *manifest != base
            && (manifest.ends_with(".leb") || manifest.ends_with(".leb2"))
            && manifest.starts_with(&want)
            && dir_files.contains(manifest)
    })
}

/// Sort and merge advisory byte ranges within one image.
///
/// Ranges separated by at most `ADVISE_PAGE_SIZE` are combined, then
/// page-aligned and clamped to `file_len`. Empty resulting ranges are
/// discarded. Saturating arithmetic handles out-of-range inputs without an
/// error. Merge each file separately because offsets are per-image.
pub fn merge_ranges(mut ranges: Vec<(u64, u64)>, file_len: u64) -> Vec<(u64, u64)> {
    if ranges.is_empty() {
        return ranges;
    }
    ranges.sort();
    let mut merged: Vec<(u64, u64)> = Vec::with_capacity(ranges.len());
    for (offset, size) in ranges {
        match merged.last_mut() {
            Some((prev_off, prev_size))
                if offset
                    <= prev_off
                        .saturating_add(*prev_size)
                        .saturating_add(ADVISE_PAGE_SIZE) =>
            {
                let end = offset
                    .saturating_add(size)
                    .max(prev_off.saturating_add(*prev_size));
                *prev_size = end.saturating_sub(*prev_off);
            }
            _ => merged.push((offset, size)),
        }
    }
    let mut out = Vec::with_capacity(merged.len());
    for (offset, size) in merged {
        let start_page = offset & !(ADVISE_PAGE_SIZE - 1);
        let end_page = offset
            .saturating_add(size)
            .saturating_add(ADVISE_PAGE_SIZE - 1)
            & !(ADVISE_PAGE_SIZE - 1);
        let end = end_page.min(file_len);
        if end > start_page {
            out.push((start_page, end - start_page));
        }
    }
    out
}

/// Touch one byte at each advisory page step within the supplied ranges.
/// Ranges are clamped to the slice. Values are not modified; this does not
/// guarantee that pages remain resident.
pub(super) fn prefault(bytes: &[u8], ranges: &[(u64, u64)]) {
    let len = bytes.len() as u64;
    for (offset, size) in ranges {
        let start = (*offset).min(len) as usize;
        let end = offset.saturating_add(*size).min(len) as usize;
        let mut i = start;
        while i < end {
            std::hint::black_box(bytes[i]);
            // `ADVISE_PAGE_SIZE >= 1`, so this strictly advances while
            // `i < end` and the loop always terminates.
            i = i.saturating_add(ADVISE_PAGE_SIZE as usize).min(end);
        }
    }
}

impl LebFile {
    /// Decode the container index and take ownership of the image vector.
    ///
    /// `path` is only a label. The buffer is moved without copying; clones share
    /// the backing, index and coefficient cache. Format/structural errors
    /// return no handle. Payload validation also occurs when coefficients are
    /// read or evaluated.
    pub fn from_bytes(path: String, bytes: Vec<u8>) -> Result<Self> {
        Self::from_bytes_with_cache_budget(path, bytes, DEFAULT_CHUNK_CACHE_BYTES)
    }

    /// Decode a process-lifetime immutable byte slice without copying it.
    ///
    /// Structural admission follows [`Self::from_bytes`]. Clones share the index
    /// and coefficient cache. These owned resources are released when their
    /// last owner drops; the static bytes retain their original lifetime.
    pub fn from_static_bytes(path: String, bytes: &'static [u8]) -> Result<Self> {
        let index = FileIndex::decode(bytes)?;
        Ok(Self {
            path,
            image: Some(Arc::new(LebImage {
                bytes: storage::ImageStorage::Static(bytes),
                index,
                cache: Arc::new(LebChunkCache::with_contiguous_window_byte_budget(
                    DEFAULT_CHUNK_CACHE_BYTES,
                )),
                sha256: std::sync::OnceLock::new(),
            })),
        })
    }

    /// Decode a shared immutable snapshot without copying its byte buffer.
    ///
    /// The caller may retain its `Arc` after a parse failure. Each construction
    /// creates an index and coefficient cache; clones of that reader share them.
    /// Structural admission follows [`Self::from_bytes`].
    pub fn from_shared_bytes(path: String, bytes: Arc<Vec<u8>>) -> Result<Self> {
        let index = FileIndex::decode(&bytes)?;
        Ok(Self {
            path,
            image: Some(Arc::new(LebImage {
                bytes: storage::ImageStorage::SharedOwned(bytes),
                index,
                cache: Arc::new(LebChunkCache::with_contiguous_window_byte_budget(
                    DEFAULT_CHUNK_CACHE_BYTES,
                )),
                sha256: std::sync::OnceLock::new(),
            })),
        })
    }

    /// Decode an owned snapshot with a retained-coefficient budget in bytes.
    ///
    /// The buffer is moved without copying. Image storage and decoding workspace
    /// are additional to this budget. Zero disables retention without
    /// prohibiting evaluation. Structural errors follow [`Self::from_bytes`].
    pub fn from_bytes_with_cache_budget(
        path: String,
        bytes: Vec<u8>,
        max_bytes: usize,
    ) -> Result<Self> {
        let index = FileIndex::decode(&bytes)?;
        Ok(Self {
            path,
            image: Some(Arc::new(LebImage {
                bytes: storage::ImageStorage::Owned(bytes),
                index,
                cache: Arc::new(LebChunkCache::with_contiguous_window_byte_budget(max_bytes)),
                sha256: std::sync::OnceLock::new(),
            })),
        })
    }

    /// SHA-256 of the acquired image, cached after a successful computation.
    ///
    /// Clones share the digest. It hashes contiguous backing or streams the
    /// acquired positional descriptor, never a reopened pathname. Failed hashes
    /// are not cached. Mapped/indexed files still require their constructor's
    /// immutable-inode contract; a digest does not establish it. An indexed
    /// cached digest performs no new metadata check; use
    /// [`Self::validate_identity`] for that probe. Closed handles return
    /// [`Error::LebReaderClosed`].
    pub fn sha256(&self) -> Result<&str> {
        let image = self.image.as_ref().ok_or(Error::LebReaderClosed)?;
        if let LebSource::Positional(file) = image.bytes.source() {
            let _guard = file
                .digest_gate
                .lock()
                .unwrap_or_else(|error| error.into_inner());
            if image.sha256.get().is_none() {
                let digest = file.sha256()?;
                let _ = image.sha256.set(digest);
            }
            // The successful hash above initialized the cell; errors publish nothing.
            return Ok(image
                .sha256
                .get()
                .expect("successful image digest")
                .as_str());
        }
        Ok(image
            .sha256
            .get_or_init(|| {
                use sha2::{Digest, Sha256};
                let LebSource::Slice(bytes) = image.bytes.source() else {
                    unreachable!()
                };
                crate::download::hex_of(&Sha256::digest(bytes))
            })
            .as_str())
    }

    /// Borrow the cache paired with this immutable image. Kept crate-local so
    /// callers cannot accidentally pair a persistent memo with another image.
    pub(crate) fn chunk_cache(&self) -> Result<&Arc<LebChunkCache>> {
        Ok(&self.image.as_ref().ok_or(Error::LebReaderClosed)?.cache)
    }

    /// Shared indexed `(metadata, retained coefficients, working reservations)`
    /// in bytes, or `None` for contiguous/closed images. Clones and images
    /// acquired together by a Runtime/Session share them; independently opened
    /// stores do not. Never sum a shared domain once per image. Working
    /// reservations include live decoded words and worst-case cold scratch,
    /// not RSS, OS cache, allocator overhead or caller-owned result vectors.
    pub fn indexed_resource_usage(&self) -> Option<(usize, usize, usize)> {
        let LebSource::Positional(file) = self.image.as_ref()?.bytes.source() else {
            return None;
        };
        Some((
            file.resources.metadata.used(),
            file.resources.retained.used(),
            file.resources.working.used(),
        ))
    }

    /// `(hits, misses, resident_entries)` shared by all clones/providers of
    /// this image, or a closed-reader error. No file reads or hashing occur.
    pub fn chunk_cache_stats(&self) -> Result<(u64, u64, usize)> {
        Ok(self.chunk_cache()?.stats())
    }

    /// Cache-accounted retained bytes and configured byte limit for this image.
    /// Image bytes, other metadata and active decoding workspace are additional.
    /// These are not total-allocation or resident-memory measurements. Closed
    /// handles return [`Error::LebReaderClosed`].
    pub fn chunk_cache_memory(&self) -> Result<(usize, usize)> {
        Ok(self.chunk_cache()?.memory_usage())
    }

    /// Retained LEB2 v2 directory vector capacities in bytes (zero when closed).
    /// Shared clones/providers observe the same metadata; no I/O, parsing or
    /// hashing occurs. This is additional to image bytes and the decoded
    /// coefficient budget, and excludes map/Arc/allocator overhead. Directories
    /// survive cooling and budget-zero requests; the last image owner frees them.
    pub fn retained_chunk_index_bytes(&self) -> usize {
        self.image
            .as_ref()
            .map_or(0, |image| image.cache.retained_chunk_index_bytes())
    }

    /// Set this image's retained decoded-byte budget. All clones/providers
    /// share the limit; shrinking releases excess retained entries. Zero
    /// disables retention without changing results. An indexed image also
    /// respects its shared domain's fixed retention ceiling. Closed handles fail.
    pub fn set_chunk_cache_byte_budget(&self, max_bytes: usize) -> Result<()> {
        self.chunk_cache()?.set_byte_budget(max_bytes);
        Ok(())
    }

    /// Read all file bytes into an owned snapshot and decode its index.
    ///
    /// Later source modification, truncation or pathname replacement does not
    /// alter the snapshot. Concurrent writes during the read are not made
    /// atomic. Missing files return [`Error::LebFileNotFound`]; other read
    /// failures return [`Error::LebCorrupt`] at offset zero. Structural errors
    /// follow [`Self::from_bytes`]. Vector capacity is reported by
    /// [`Self::allocated_image_bytes`], not guaranteed to equal file length.
    pub fn open(path: &std::path::Path) -> Result<Self> {
        let key = path.to_string_lossy().into_owned();
        let bytes = std::fs::read(path).map_err(|e| {
            if e.kind() == std::io::ErrorKind::NotFound {
                Error::LebFileNotFound { path: key.clone() }
            } else {
                Error::LebCorrupt {
                    what: "LEB file read".to_string(),
                    offset: 0,
                }
            }
        })?;
        Self::from_bytes(key, bytes)
    }

    /// Open a regular immutable file as a read-only mapping and decode its index.
    ///
    /// Mapping avoids a whole-file heap copy; parsing still reads metadata.
    /// Mapping/structural failures refuse the open without an owned fallback.
    /// Clones share backing, descriptor, index and coefficient cache.
    ///
    /// # Safety
    /// The acquired inode must not be modified or truncated until all image
    /// owners and borrowed views are dropped. Read-only permissions or a digest
    /// do not guarantee this. Atomic pathname replacement is permitted when it
    /// leaves the acquired inode unchanged. Use [`Self::open`] to own a copy
    /// whose bytes survive later source mutations.
    pub unsafe fn open_mapped(path: &std::path::Path) -> Result<Self> {
        // SAFETY: the public caller supplies the immutable-inode lifetime contract.
        let bytes = unsafe { storage::ImageStorage::open_mapped(path) }?;
        let index = FileIndex::decode(bytes.bytes()?)?;
        Ok(Self {
            path: path.to_string_lossy().into_owned(),
            image: Some(Arc::new(LebImage {
                bytes,
                index,
                cache: Arc::new(LebChunkCache::with_image_byte_budget(
                    DEFAULT_CHUNK_CACHE_BYTES,
                )),
                sha256: std::sync::OnceLock::new(),
            })),
        })
    }

    /// Open a regular immutable LEB file using bounded positional reads.
    ///
    /// Metadata, retained coefficients, decoding workspace and zstd windows are
    /// subject to `limits`, not a total process-memory cap. Invalid limits,
    /// format, I/O and resource failures propagate without an owned/mapped
    /// fallback. [`Self::bytes`] and [`Self::read_range`] cannot borrow bytes
    /// from this backing.
    ///
    /// The acquired inode must not be rewritten or truncated while any image
    /// owner is live. Atomic replacement may publish a different inode.
    /// Size/mtime checks detect some mutations but do not prove immutability;
    /// unlike [`Self::open`], this is not an arbitrary-file snapshot.
    pub fn open_indexed(path: &std::path::Path, limits: IndexedLimits) -> Result<Self> {
        Self::open_indexed_with_resources(path, IndexedResources::new(limits))
    }

    /// Attach another immutable image to an existing indexed allocation domain.
    pub(crate) fn open_indexed_with_resources(
        path: &std::path::Path,
        resources: Arc<IndexedResources>,
    ) -> Result<Self> {
        let mut file = access::PositionalFile::open(path, Arc::clone(&resources))?;
        let index = file.read_index()?;
        let cache = Arc::clone(&file.cache);
        Ok(Self {
            path: path.to_string_lossy().into_owned(),
            image: Some(Arc::new(LebImage {
                bytes: storage::ImageStorage::Positional(Box::new(file)),
                index,
                cache,
                sha256: std::sync::OnceLock::new(),
            })),
        })
    }

    /// Whether this live handle owns immutable heap bytes, without mappings or
    /// positional admission. Session reuse must also establish native dispatch.
    pub(crate) fn is_owned_image(&self) -> bool {
        self.image.as_ref().is_some_and(|image| {
            matches!(
                &image.bytes,
                storage::ImageStorage::Owned(_) | storage::ImageStorage::SharedOwned(_)
            )
        })
    }

    /// Path or image label supplied when this handle was constructed.
    pub fn path(&self) -> &str {
        &self.path
    }

    /// True after [`LebFile::close`]; `eval`/`warm` then fail with
    /// [`Error::LebReaderClosed`] while `cool`/`close` stay safe.
    pub fn is_closed(&self) -> bool {
        self.image.is_none()
    }

    /// Logical image length for any backing, zero when closed. Not an RSS metric.
    pub fn len(&self) -> usize {
        self.image
            .as_ref()
            .map_or(0, |image| image.bytes.source().len())
    }

    /// Retained image-vector capacity (zero for static/mapped/indexed/closed),
    /// excluding index/cache metadata. Shared clones count the same allocation.
    /// This can exceed [`Self::len`] for caller-provided overallocated vectors.
    pub fn allocated_image_bytes(&self) -> usize {
        self.image
            .as_ref()
            .map_or(0, |image| image.bytes.owned_capacity())
    }

    /// Mapped address-space length, zero for non-mapped or closed images.
    /// This is not resident RAM or filesystem-cache usage; the query touches
    /// no image pages.
    pub fn mapped_file_bytes(&self) -> usize {
        self.image
            .as_ref()
            .map_or(0, |image| image.bytes.mapped_len())
    }

    /// Whether [`Self::len`] is zero. Closed handles also return `true`;
    /// use [`Self::is_closed`] to distinguish them.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Decoded container index, or [`Error::LebReaderClosed`].
    pub fn index(&self) -> Result<&FileIndex> {
        Ok(&self.image.as_ref().ok_or(Error::LebReaderClosed)?.index)
    }

    /// Borrowed contiguous image bytes, or [`Error::LebReaderClosed`]. Indexed
    /// images return [`Error::LebBytesUnavailable`], never a whole-file copy.
    /// The borrow cannot outlive its owner; taking it does not fault every page.
    pub fn bytes(&self) -> Result<&[u8]> {
        self.image
            .as_ref()
            .ok_or(Error::LebReaderClosed)?
            .bytes
            .bytes()
    }

    /// Probe the acquired positional descriptor for detectable size/mtime
    /// changes, even when decoded words or a digest are already cached. This
    /// does not prove immutability. Contiguous images need no positional probe;
    /// closed handles return [`Error::LebReaderClosed`].
    pub fn validate_identity(&self) -> Result<()> {
        match self.source()? {
            LebSource::Slice(_) => Ok(()),
            LebSource::Positional(file) => file.validate_identity(),
        }
    }

    /// Internal byte acquisition view, paired with this image's index/cache.
    pub(crate) fn source(&self) -> Result<LebSource<'_>> {
        Ok(self
            .image
            .as_ref()
            .ok_or(Error::LebReaderClosed)?
            .bytes
            .source())
    }

    /// Container flavor, or [`Error::LebReaderClosed`].
    pub fn format(&self) -> Result<LebFormat> {
        Ok(self.index()?.header.format)
    }

    /// Borrow `size` bytes at absolute byte `offset` from contiguous backing.
    /// Out-of-bounds/unrepresentable ranges return [`Error::LebCorrupt`];
    /// indexed backing returns [`Error::LebBytesUnavailable`]. A closed handle
    /// returns [`Error::LebReaderClosed`]. No image copy is made.
    pub fn read_range(&self, offset: u64, size: u64) -> Result<&[u8]> {
        let bytes = self.bytes()?;
        read_blob(bytes, offset, size, "LEB partial read")
    }

    /// True when `body_id` has an index entry, or
    /// [`Error::LebReaderClosed`].
    pub fn has_body(&self, body_id: i32) -> Result<bool> {
        Ok(self.index()?.has_body(body_id))
    }

    /// Inclusive stored body coverage in TT Julian days, or `None` for an
    /// absent body. Closed handles return [`Error::LebReaderClosed`].
    pub fn body_coverage(&self, body_id: i32) -> Result<Option<(f64, f64)>> {
        Ok(self.index()?.body_coverage(body_id))
    }

    /// The file header's declared start and end TT Julian days.
    /// This is not recomputed from body intervals and does not establish every
    /// body's coverage or all support-channel availability. Closed handles return
    /// [`Error::LebReaderClosed`].
    pub fn jd_range(&self) -> Result<(f64, f64)> {
        let index = self.index()?;
        Ok((index.header().jd_start(), index.header().jd_end()))
    }

    /// Evaluate an LEB1 body's stored coordinates and TT-day rates at TT JD `jd`.
    ///
    /// Coordinate units, body coverage and errors follow [`FileIndex::eval_body`].
    /// This raw reader refuses non-LEB1 formats, but is not limited to owned
    /// backing. Closed handles return [`Error::LebReaderClosed`].
    pub fn eval_body(&self, body_id: i32, jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        self.evaluation()?.eval_leb1(body_id, jd)
    }

    /// Two-part-epoch form of [`LebFile::eval_body`] (the split-epoch layer).
    pub fn eval_body_two_part(
        &self,
        body_id: i32,
        epoch: TwoPartEpoch,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        self.evaluation()?.eval_leb1_two_part(body_id, epoch)
    }

    /// Evaluate an LEB2 body's stored coordinates and TT-day rates at TT JD `jd`.
    /// Coordinate units, coverage and decoding/numerical errors follow
    /// [`FileIndex::eval_leb2_body`]. The shared coefficient cache is used; this
    /// is not limited to owned backing. Closed handles return
    /// [`Error::LebReaderClosed`].
    pub fn eval_leb2_body(&self, body_id: i32, jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        self.evaluation()?.eval_leb2(body_id, jd)
    }

    /// Two-part-epoch form of [`LebFile::eval_leb2_body`] (the split-epoch layer).
    pub fn eval_leb2_body_two_part(
        &self,
        body_id: i32,
        epoch: TwoPartEpoch,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        self.evaluation()?.eval_leb2_two_part(body_id, epoch)
    }

    /// Evaluate stored coordinates and TT-day rates using the detected format.
    /// LEB1 uses [`Self::eval_body`]; LEB2 uses [`Self::eval_leb2_body`] and the
    /// shared coefficient cache. Their units, coverage and errors apply; closed
    /// handles return [`Error::LebReaderClosed`].
    pub fn eval_body_any(&self, body_id: i32, jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        self.evaluation()?.eval_any(body_id, jd)
    }

    /// Two-part-epoch form of [`LebFile::eval_body_any`].
    pub fn eval_body_any_two_part(
        &self,
        body_id: i32,
        epoch: TwoPartEpoch,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        self.evaluation()?.eval_any_two_part(body_id, epoch)
    }

    /// Evaluate [`Self::eval_body_any`] in input order, with one state per epoch.
    /// The first failing epoch returns its scalar-call error rather than a
    /// partial result vector.
    pub fn eval_body_any_many(
        &self,
        body_id: i32,
        jds: &[f64],
    ) -> Result<Vec<(Vec<f64>, Vec<f64>)>> {
        batch_eval(jds, |jd| self.eval_body_any(body_id, jd))
    }

    /// Multi-epoch form of [`LebFile::eval_body_any_two_part`].
    pub fn eval_body_any_two_part_many(
        &self,
        body_id: i32,
        epochs: &[TwoPartEpoch],
    ) -> Result<Vec<(Vec<f64>, Vec<f64>)>> {
        batch_eval_epochs(epochs, |epoch| self.eval_body_any_two_part(body_id, epoch))
    }

    /// Whether this image declares usable stored nutation data.
    ///
    /// A missing or undecodable header, zero segments, or a non-positive or
    /// non-finite interval returns `false`. This probes the header only: it does
    /// not validate coefficient payloads or establish coverage for an epoch.
    /// Positional I/O/resource failures propagate, and a closed handle returns
    /// [`Error::LebReaderClosed`]. No analytical nutation model is substituted.
    pub fn has_nutation(&self) -> Result<bool> {
        self.evaluation()?.has_nutation()
    }

    /// Evaluate stored `(dpsi, deps)` nutation angles in radians at TT JD `jd`.
    /// Missing/malformed data, non-finite or uncovered epochs, backing failures
    /// and closed handles return errors. See [`FileIndex::eval_nutation`].
    pub fn eval_nutation(&self, jd: f64) -> Result<(f64, f64)> {
        self.evaluation()?.nutation(jd)
    }

    /// Interpolate stored TT−UT1 in days at the table coordinate `jd`.
    /// No clock conversion or freshness test is performed. Outside the table,
    /// the nearest endpoint is used, not a predictive model. Nonfinite epochs,
    /// missing/malformed data, backing failures and closed handles return errors.
    /// See [`FileIndex::delta_t`] and [`eval_delta_t`] for the raw table contract.
    pub fn delta_t(&self, jd: f64) -> Result<f64> {
        self.evaluation()?.delta_t(jd)
    }

    /// Read a stored catalog entry by `star_id`. Missing catalogs or identifiers
    /// return [`Error::LebStarNotFound`]; malformed records and backing failures
    /// propagate. Closed handles return [`Error::LebReaderClosed`].
    /// See [`FileIndex::get_star`].
    pub fn get_star(&self, star_id: i32) -> Result<StarEntry> {
        self.evaluation()?.star(star_id)
    }

    /// Merged advisory `(byte_offset, byte_count)` ranges for the TT-day window
    /// `[jd_start, jd_end]`. Metadata/backing errors propagate; closed handles
    /// return [`Error::LebReaderClosed`]. This does not prefetch data.
    pub fn warm_ranges(&self, jd_start: f64, jd_end: f64) -> Result<Vec<(u64, u64)>> {
        let (index, bytes) = (self.index()?, self.source()?);
        let ranges = index.warm_ranges_source(bytes, jd_start, jd_end)?;
        Ok(merge_ranges(ranges, bytes.len() as u64))
    }

    /// Compute advisory ranges for the TT-day window `[jd_start, jd_end]` and
    /// return their merged count.
    ///
    /// Heap/static backing touches pages; mapped backing uses best-effort OS
    /// advice where supported. Indexed backing does not prefetch coefficients.
    /// This does not evaluate coefficients or guarantee residency. Metadata or
    /// backing errors propagate; closed handles return [`Error::LebReaderClosed`].
    pub fn warm(&self, jd_start: f64, jd_end: f64) -> Result<usize> {
        let (index, bytes) = (self.index()?, self.source()?);
        let ranges = merge_ranges(
            index.warm_ranges_source(bytes, jd_start, jd_end)?,
            bytes.len() as u64,
        );
        self.image
            .as_ref()
            .ok_or(Error::LebReaderClosed)?
            .bytes
            .warm(&ranges);
        Ok(ranges.len())
    }

    /// Clear shared retained coefficients, including for other clones. Active
    /// evaluations keep their borrowed words; later queries may decode again.
    /// Retained chunk directories and owned image bytes remain.
    ///
    /// Mapped and indexed backing receive best-effort page-release advice without
    /// invalidating the mapping. No immediate RSS or filesystem-cache reduction
    /// is guaranteed. Closed handles and repeated cooling are safe.
    pub fn cool(&self) {
        if let Some(image) = &self.image {
            image.cache.clear();
            image.bytes.release_pages();
        }
    }

    /// Advise release of mapped/indexed clean file-cache pages where
    /// supported, without clearing retained coefficients. Other backings are
    /// unchanged. Safe on closed handles; other readers may retain or refault
    /// pages. Use [`Self::cool`] to clear coefficients as well.
    pub fn release_file_pages(&self) {
        if let Some(image) = &self.image {
            image.bytes.release_pages();
        }
    }

    /// Release this handle's reference to the shared image. Other clones remain
    /// usable; backing, index and cache follow their remaining owners. Static
    /// bytes or a caller-retained shared buffer can outlive every handle.
    /// Repeated close is a no-op.
    pub fn close(&mut self) {
        self.image = None;
    }
}

/// Open an owned snapshot with header-based format detection.
/// The extension does not select a format. File-reading and structural
/// errors follow [`LebFile::open`]; payload validation also occurs when
/// coefficients are evaluated.
pub fn open_leb(path: &std::path::Path) -> Result<LebFile> {
    LebFile::open(path)
}

#[cfg(test)]
mod digest_tests {
    use super::*;

    #[test]
    fn concurrent_clones_share_one_digest_and_new_images_do_not() {
        let path = crate::test_data::core_leb_path();
        let bytes = std::fs::read(path).unwrap();
        let file =
            LebFile::from_bytes("/not/reopened/base_core.leb2".into(), bytes.clone()).unwrap();
        assert!(file.image.as_ref().unwrap().sha256.get().is_none());
        let addresses = std::thread::scope(|scope| {
            let calls: Vec<_> = (0..8)
                .map(|_| {
                    scope.spawn(|| {
                        let clone = file.clone();
                        let digest = clone.sha256().unwrap();
                        assert!(crate::test_data::admitted_core_hash(digest));
                        digest.as_ptr() as usize
                    })
                })
                .collect();
            calls
                .into_iter()
                .map(|call| call.join().unwrap())
                .collect::<Vec<_>>()
        });
        assert!(addresses.iter().all(|address| *address == addresses[0]));
        let mut changed = bytes;
        changed.push(0);
        let other = LebFile::from_bytes(file.path().into(), changed).unwrap();
        assert!(other.image.as_ref().unwrap().sha256.get().is_none());
        assert_ne!(other.sha256().unwrap(), file.sha256().unwrap());
        assert_eq!(
            other.eval_body_any(0, 2451545.0).unwrap(),
            file.eval_body_any(0, 2451545.0).unwrap()
        );
        let mut closed = file.clone();
        closed.close();
        assert!(matches!(closed.sha256(), Err(Error::LebReaderClosed)));
        assert_eq!(file.sha256().unwrap().as_ptr() as usize, addresses[0]);
    }
}

#[cfg(test)]
#[path = "tests/fused_decode.rs"]
mod fused_decode_tests;

#[cfg(test)]
#[path = "tests/manifest_identity.rs"]
mod manifest_identity_tests;

#[cfg(test)]
mod owned_session_image_tests {
    use super::*;

    #[test]
    fn owned_capability_uses_backing_and_close_is_handle_local() {
        let bytes = crate::test_data::core_bytes();
        let mut owned = LebFile::from_bytes("base_core.leb2".into(), bytes.to_vec()).unwrap();
        let held = owned.clone();
        assert!(owned.is_owned_image());
        owned.close();
        assert!(!owned.is_owned_image());
        assert!(held.is_owned_image());
        assert!(held.index().is_ok());
        let shared = LebFile::from_shared_bytes(
            "base_core.leb2".into(),
            std::sync::Arc::new(bytes.to_vec()),
        )
        .unwrap();
        assert!(shared.is_owned_image());
        let static_image = LebFile::from_static_bytes("base_core.leb2".into(), bytes).unwrap();
        // A generic Slice view does not qualify: backing kind is authoritative.
        assert!(!static_image.is_owned_image());
    }
}

#[cfg(test)]
mod nonowned_session_image_tests {
    use super::*;

    #[test]
    fn indexed_and_mapped_images_never_supply_owned_session_capability() {
        let bytes = crate::test_data::core_bytes();
        let path = crate::test_data::core_leb_path();
        let indexed = LebFile::open_indexed(&path, IndexedLimits::default()).unwrap();
        assert!(!indexed.is_owned_image());
        assert!(indexed.indexed_resource_usage().is_some());
        // SAFETY: this test maps its archive's immutable read-only fixture;
        // no test writes, truncates or replaces it while the mapping lives.
        // The reader validates length and owns the map/file lifetime. Every
        // alias is read-only and no borrowed slice outlives the mapped handle.
        let mapped = unsafe { LebFile::open_mapped(&path) }.unwrap();
        assert!(!mapped.is_owned_image());
        assert_eq!(mapped.bytes().unwrap(), bytes);
        drop(mapped);
        drop(indexed);
    }
}
