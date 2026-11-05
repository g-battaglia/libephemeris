//! Index implementation for the parent module.
//! Public entry points and types are exposed through the parent API.

use super::*;

impl Segment {
    /// Segment name as text: trailing spaces and NULs trimmed,
    /// non-UTF-8 bytes replaced (names are advisory only).
    pub fn name_str(&self) -> String {
        String::from_utf8_lossy(&self.name)
            .trim_end_matches([' ', '\0'])
            .to_string()
    }
}

impl SpkIndex {
    /// Decode every segment descriptor of the kernel image `bytes`.
    ///
    /// Walks the summary-record chain from `FWARD`, decoding at most
    /// [`MAX_SUMMARIES_PER_RECORD`] descriptors per record and reading
    /// each paired name record. A `FWARD` of 0 means no segments.
    pub fn parse(bytes: &[u8]) -> Result<Self> {
        let header = FileRecord::parse(bytes)?;
        Self::parse_directory(header, bytes.len(), usize::MAX, |record| {
            let base = record_offset(record, bytes)?;
            Ok(bytes[base..base + RECORD_BYTES]
                .try_into()
                .expect("validated DAF record"))
        })
    }

    /// Read only the file record and paired summary/name records from one
    /// acquired file. Coefficient payloads are not copied or evaluated; their
    /// declared addresses retain the same actual-file-length validation as an
    /// in-memory index. Discovery opens a full owned image only on registration.
    pub(crate) fn read_directory(path: &std::path::Path) -> Result<Self> {
        let mut file = std::fs::File::open(path).map_err(|_| Error::SpkCorrupt {
            what: "SPK directory open".to_owned(),
            offset: 0,
        })?;
        let length = file
            .metadata()
            .map_err(|_| Error::SpkCorrupt {
                what: "SPK directory length".to_owned(),
                offset: 0,
            })?
            .len();
        let length = usize::try_from(length).map_err(|_| Error::SpkCorrupt {
            what: "SPK directory length".to_owned(),
            offset: 0,
        })?;
        Self::read_directory_from(&mut file, length)
    }

    /// Decode through bounded record reads, sharing every descriptor and
    /// directory guard with byte-image parsing. `length` is the acquired file
    /// length; read/seek failures remain corrupt-directory refusals.
    pub(super) fn read_directory_from(
        reader: &mut (impl std::io::Read + std::io::Seek),
        length: usize,
    ) -> Result<Self> {
        Self::read_directory_bounded_from(reader, length, usize::MAX)
    }

    /// The same directory decoder with an explicit metadata admission bound.
    /// Existing image/discovery callers retain their original unlimited policy.
    pub(super) fn read_directory_bounded_from(
        reader: &mut (impl std::io::Read + std::io::Seek),
        length: usize,
        maximum_segments: usize,
    ) -> Result<Self> {
        use std::io::SeekFrom;
        let mut read = |record| {
            let base = super::reader::record_offset_for_length(record, length)?;
            let failed = |_| Error::SpkCorrupt {
                what: "SPK directory read".to_owned(),
                offset: base,
            };
            reader.seek(SeekFrom::Start(base as u64)).map_err(failed)?;
            let mut bytes = [0u8; RECORD_BYTES];
            reader.read_exact(&mut bytes).map_err(failed)?;
            Ok(bytes)
        };
        let first = read(1)?;
        let header = FileRecord::parse(&first)?;
        Self::parse_directory(header, length, maximum_segments, read)
    }

    /// Walk the shared strict summary chain over bounded physical records.
    fn parse_directory(
        header: FileRecord,
        length: usize,
        maximum_segments: usize,
        mut read_record: impl FnMut(u64) -> Result<[u8; RECORD_BYTES]>,
    ) -> Result<Self> {
        let endian = header.endian;
        if header.first_summary_record == 0 {
            if header.last_summary_record != 0 {
                return Err(Error::SpkCorrupt {
                    what: "DAF summary chain (empty head, nonzero tail)".to_string(),
                    offset: OFF_BWARD,
                });
            }
            return Ok(Self {
                segments: Vec::new(),
            });
        }
        let total_words = length / 8;
        // Step budget: at most one visit per physical record, so a
        // corrupted NEXT cycle terminates instead of looping forever.
        let mut remaining = length / RECORD_BYTES + 1;
        let mut record = header.first_summary_record;
        let mut segments = Vec::new();
        loop {
            if remaining == 0 {
                return Err(Error::SpkCorrupt {
                    what: "DAF summary chain (cycle)".to_string(),
                    offset: super::reader::record_offset_for_length(record, length)?,
                });
            }
            remaining -= 1;
            let base = super::reader::record_offset_for_length(record, length)?;
            let summary = read_record(record)?;
            // Control words are stored as doubles; they must be exact
            // non-negative integers.
            let control = |word: usize| -> Result<u64> {
                let value = endian.read_f64(&summary[word * 8..word * 8 + 8]);
                if value < 0.0 || value.fract() != 0.0 || value > u32::MAX as f64 {
                    return Err(Error::SpkCorrupt {
                        what: "DAF summary control word".to_string(),
                        offset: base + word * 8,
                    });
                }
                Ok(value as u64)
            };
            let next = control(0)?;
            let nsum = control(2)?;
            if nsum as usize > MAX_SUMMARIES_PER_RECORD {
                return Err(Error::SpkCorrupt {
                    what: "DAF summary count".to_string(),
                    offset: base + 16,
                });
            }
            if nsum as usize > maximum_segments.saturating_sub(segments.len()) {
                return Err(Error::SpkCorrupt {
                    what: "SPK bounded directory metadata limit".into(),
                    offset: base + 16,
                });
            }
            let _name_base =
                super::reader::record_offset_for_length(record + 1, length).map_err(|_| {
                    Error::SpkCorrupt {
                        what: "DAF name record".to_string(),
                        offset: base,
                    }
                })?;
            let names = read_record(record + 1)?;
            for i in 0..nsum as usize {
                let local = (3 + i * SPK_SUMMARY_WORDS) * 8;
                let off = base + local;
                let start_et = endian.read_f64(&summary[local..local + 8]);
                let end_et = endian.read_f64(&summary[local + 8..local + 16]);
                if !start_et.is_finite() || !end_et.is_finite() || start_et > end_et {
                    return Err(Error::SpkCorrupt {
                        what: "SPK segment epochs".to_string(),
                        offset: off,
                    });
                }
                let mut ints = [0i32; 6];
                for (k, slot) in ints.iter_mut().enumerate() {
                    // Integer pairs are equivalenced into doubles in
                    // file order: low half first in memory.
                    let w = local + 16 + (k / 2) * 8 + (k % 2) * 4;
                    *slot = endian.read_i32(&summary[w..w + 4]);
                }
                let [target, center, frame, seg_type, init, fin] = ints;
                if init < 1
                    || fin < init
                    || (fin as u64) > total_words as u64
                    || init as u64 > total_words as u64
                {
                    return Err(Error::SpkCorrupt {
                        what: "SPK segment addresses".to_string(),
                        offset: off + 16,
                    });
                }
                let name_off = i * SPK_NAME_BYTES;
                let mut name = [0u8; SPK_NAME_BYTES];
                name.copy_from_slice(&names[name_off..name_off + SPK_NAME_BYTES]);
                segments.push(Segment {
                    descriptor: Descriptor {
                        start_et,
                        end_et,
                        target,
                        center,
                        frame,
                        seg_type,
                        initial_address: init,
                        final_address: fin,
                    },
                    name,
                });
            }
            if next == 0 {
                break;
            }
            record = next;
        }
        Ok(Self { segments })
    }
    /// Segments for `target` in file (priority) order.
    pub fn segments_for(&self, target: i32) -> impl Iterator<Item = &Segment> {
        self.segments
            .iter()
            .filter(move |s| s.descriptor.target == target)
    }

    /// Every segment in the index, in file order (e.g. for inventory
    /// and registration-time chain validation over all targets).
    pub fn all_segments(&self) -> impl Iterator<Item = &Segment> {
        self.segments.iter()
    }

    /// Segments of one data type in file order (e.g. all type 21
    /// segments of a mixed kernel).
    pub fn segments_of_type(&self, seg_type: i32) -> impl Iterator<Item = &Segment> {
        self.segments
            .iter()
            .filter(move |s| s.descriptor.seg_type == seg_type)
    }

    /// Highest-priority segment covering `et` for `target`: the last
    /// segment in file order whose closed coverage interval contains
    /// `et`. Returns `Ok(None)` when no segment covers the epoch (the
    /// release answers that miss with `None`/range error downstream,
    /// never with invented data).
    pub fn select(&self, target: i32, et: f64) -> Result<Option<&Segment>> {
        if !et.is_finite() {
            return Err(Error::NonFiniteArgument { value: et });
        }
        Ok(self
            .segments_for(target)
            .filter(|s| s.descriptor.covers(et))
            .last())
    }

    /// Merged closed coverage intervals for `target`, ascending.
    /// Touching intervals (`end >= next start`) merge; disjoint ones
    /// stay separate, mirroring the `spkcov` window contract.
    pub fn coverage(&self, target: i32) -> Vec<(f64, f64)> {
        let mut spans: Vec<(f64, f64)> = self
            .segments_for(target)
            .map(|s| (s.descriptor.start_et, s.descriptor.end_et))
            .collect();
        spans.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
        let mut merged: Vec<(f64, f64)> = Vec::with_capacity(spans.len());
        for (start, end) in spans {
            if let Some(last) = merged.last_mut() {
                if start <= last.1 {
                    last.1 = last.1.max(end);
                    continue;
                }
            }
            merged.push((start, end));
        }
        merged
    }

    /// Number of segments in the index (all targets).
    pub fn len(&self) -> usize {
        self.segments.len()
    }

    /// True when the kernel holds no segments.
    pub fn is_empty(&self) -> bool {
        self.segments.is_empty()
    }
}

/// Registration-time validation for `target`: single declared center
/// plus served reference frame.
pub fn resolve_chain(index: &SpkIndex, target: i32) -> Result<TargetChain> {
    let center = declared_center(index, target)?;
    let frame = check_frame(index, target)?;
    Ok(TargetChain {
        target,
        center,
        frame,
    })
}

impl SpkKernel {
    /// Decode a kernel image under the registry key `path`.
    ///
    /// Parsing is eager: a corrupt or non-SPK image fails here
    /// ([`Error::SpkUnknownFormat`] / [`Error::SpkUnsupportedNdNi`] /
    /// [`Error::SpkCorrupt`]), never at the first evaluation. An empty
    /// registry key is [`Error::SpkFileNotFound`].
    pub fn decode(path: &str, bytes: Vec<u8>) -> Result<Self> {
        if path.is_empty() {
            return Err(Error::SpkFileNotFound {
                path: path.to_string(),
            });
        }
        let endian = FileRecord::parse(&bytes)?.endian;
        let index = SpkIndex::parse(&bytes)?;
        Ok(Self {
            path: path.to_string(),
            bytes: bytes.into(),
            endian,
            index: Arc::new(index),
            source_file: None,
            file_image_owner: None,
        })
    }

    /// Registry key this kernel was stored under (the exact string
    /// passed at registration; no canonicalization is applied).
    pub fn path(&self) -> &str {
        &self.path
    }

    /// Owned kernel bytes (evaluation reads segment words from these).
    pub fn bytes(&self) -> &[u8] {
        &self.bytes
    }

    /// SHA-256 of the immutable, already admitted image bytes.
    ///
    /// File-backed clones and fully compared file hints share one fixed digest
    /// on their existing owner. Caller-byte kernels calculate normally. This
    /// method performs no file admission, I/O, registration or numeric work.
    pub(crate) fn image_sha256(&self) -> [u8; 32] {
        use sha2::{Digest, Sha256};
        match &self.file_image_owner {
            Some(owner) => *owner.sha256(&self.bytes),
            None => Sha256::digest(&self.bytes).into(),
        }
    }

    /// Binary format declared by the image's file record.
    pub fn endian(&self) -> DafEndian {
        self.endian
    }

    /// Decoded summary index (file-order segments).
    pub fn index(&self) -> &SpkIndex {
        &self.index
    }
}

impl std::fmt::Debug for SpkKernel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SpkKernel")
            .field("path", &self.path)
            .field("endian", &self.endian)
            .field("segments", &self.index.len())
            .finish_non_exhaustive()
    }
}

/// Effective SPK cache directory from explicit inputs: the programmatic override first,
/// then `LIBEPHEMERIS_SPK_DIR`, then the TOML `spk_dir` key, else `None`
/// (the caller falls back to its default cache location).
///
/// Empty strings are unset; whitespace-only strings are retained. Values are
/// copied as supplied without environment/TOML acquisition, home expansion,
/// path normalization, existence checks or directory creation. Download policy
/// and the eventual default path belong to the caller.
pub fn resolve_spk_cache_dir(
    explicit: Option<&str>,
    env: Option<&str>,
    toml: Option<&str>,
) -> Option<String> {
    for dir in [explicit, env, toml].into_iter().flatten() {
        if !dir.is_empty() {
            return Some(dir.to_string());
        }
    }
    None
}

/// Body stem rule of `AutoSpkConfig.get_cache_filename`: lowercase with non-alphanumeric (other than `-_`)
/// characters mapped to `_`, with *no* run-collapsing or stripping —
/// deliberately different from [`sanitize_filename`]. The trailing
/// `{md5(config)[:8]}.bsp` scheme needs an MD5 dependency and is
/// completed in the download layer; [`AutoSpkConfig::config_string`] exposes the
/// exact hash preimage.
pub fn auto_spk_cache_stem(body_id: &str) -> String {
    body_id
        .to_lowercase()
        .chars()
        .map(|c| {
            if c.is_alphanumeric() || c == '-' || c == '_' {
                c
            } else {
                '_'
            }
        })
        .collect()
}

/// Acquisition eligibility comparison for supplied Julian-day envelopes.
///
/// All endpoints must use the same coordinate convention. Accepts
/// `cov_start <= req_start + 1.5 days` and
/// `cov_end >= req_end - 0.5 days`, allowing the established light-time/stencil
/// and service date-rounding margins. This is not a strict guarantee of
/// complete requested-epoch, target or center-chain coverage. It neither reads
/// a file nor validates interval ordering/finiteness; comparisons govern the
/// result for unchecked inputs.
pub fn spk_coverage_covers_request(
    cov_start: f64,
    cov_end: f64,
    req_start: f64,
    req_end: f64,
) -> bool {
    cov_start <= req_start + SPK_COVERAGE_START_TOLERANCE_DAYS
        && cov_end >= req_end - SPK_COVERAGE_END_TOLERANCE_DAYS
}

/// Generate an SPK cache filename:
/// `{sanitized}_{truncated start}_{truncated end}.bsp`, where the JDs
/// truncate toward zero (`int` semantics — day-granular lookup needs
/// no more).
pub fn generate_spk_cache_filename(body_id: &str, jd_start: f64, jd_end: f64) -> String {
    format!(
        "{}_{}_{}.bsp",
        sanitize_filename(body_id),
        jd_start.trunc() as i64,
        jd_end.trunc() as i64
    )
}

/// SPK date window for a precision tier: tiers are cumulative
/// (`base` ⊂ `medium` ⊂ `extended`), so the window is the union with
/// every lower-priority tier's window — a fresh `medium` install
/// otherwise fetches a *narrower* asteroid-SPK window than `base` and
/// silently drops to the Keplerian fallback at e.g. 1870.
pub fn spk_date_range_for_tier(tier: crate::context::Tier) -> (&'static str, &'static str) {
    use crate::context::Tier;
    match tier {
        Tier::Base => ("1850-01-01", "2150-01-01"),
        Tier::Medium => ("1850-01-01", "2150-01-01"),
        Tier::Extended => ("1600-01-01", "2500-01-01"),
    }
}

/// SPK date window for a tier name: unknown names are
/// [`Error::InvalidTier`], mirroring the release `ValueError`.
pub fn spk_date_range_for_tier_name(name: &str) -> Result<(&'static str, &'static str)> {
    match crate::context::Tier::from_name(name) {
        Some(tier) => Ok(spk_date_range_for_tier(tier)),
        None => Err(Error::InvalidTier {
            value: name.to_string(),
        }),
    }
}

/// Whether a body is in JPL's major-body index: always `false` — the
/// `"Name;"` syntax bypasses the restriction, so this is kept only for
/// backward compatibility.
pub fn is_in_jpl_major_body_index(_body_id: i32) -> bool {
    false
}

/// Hand-built chain-complete satellite-kernel image for planetary-moon
/// registry tests: Io (NAIF 501, center Jupiter barycenter 5) plus the
/// center-5 segment itself (center SSB 0), so the SSB chain resolves
/// entirely inside the image. Index-only like
/// [`moon_test_kernel_bytes`] (addresses point at the zeroed element
/// record: registration and chain resolution work, coefficient
/// evaluation does not).
#[cfg(test)]
pub(crate) fn moon_chain_test_kernel_bytes() -> Vec<u8> {
    let seg = |start: f64, end: f64, target: i32, center: i32, name: &'static str| SegSpec {
        start,
        end,
        target,
        center,
        frame: 1,
        seg_type: 2,
        name,
    };
    build_image(
        DafEndian::Little,
        &[vec![
            seg(0.0, 10_000_000.0, 501, 5, "IO"),
            seg(0.0, 10_000_000.0, 5, 0, "JUPBARY"),
        ]],
    )
}

/// Hand-built satellite-kernel image for planetary-moon registry tests: Io
/// (NAIF 501, center Jupiter barycenter 5) and Titan (NAIF 606,
/// center Saturn barycenter 6) type-2 J2000 segments plus one
/// non-satellite segment (NAIF 3, excluded from auto-detect).
#[cfg(test)]
pub(crate) fn moon_test_kernel_bytes() -> Vec<u8> {
    let seg = |start: f64, end: f64, target: i32, center: i32, name: &'static str| SegSpec {
        start,
        end,
        target,
        center,
        frame: 1,
        seg_type: 2,
        name,
    };
    build_image(
        DafEndian::Little,
        &[vec![
            seg(0.0, 10_000_000.0, 501, 5, "IO"),
            seg(0.0, 20_000_000.0, 606, 6, "TITAN"),
            seg(0.0, 5_000_000.0, 3, 0, "EARTH"),
        ]],
    )
}

#[cfg(test)]
mod directory_read_tests {
    use super::*;
    use std::collections::BTreeMap;
    use std::io::{Read, Seek, SeekFrom};

    /// A file that exposes directory records only. Any coefficient read fails,
    /// including reads across a large sparse interval between header and index.
    struct DirectoryOnly {
        records: BTreeMap<u64, [u8; RECORD_BYTES]>,
        position: u64,
        bytes_read: usize,
    }
    impl Read for DirectoryOnly {
        fn read(&mut self, buffer: &mut [u8]) -> std::io::Result<usize> {
            let record = self.position / RECORD_BYTES as u64;
            let within = (self.position % RECORD_BYTES as u64) as usize;
            let bytes = self
                .records
                .get(&record)
                .ok_or_else(|| std::io::Error::other("coefficient payload must not be read"))?;
            let count = buffer.len().min(RECORD_BYTES - within);
            buffer[..count].copy_from_slice(&bytes[within..within + count]);
            self.position += count as u64;
            self.bytes_read += count;
            Ok(count)
        }
    }
    impl Seek for DirectoryOnly {
        fn seek(&mut self, seek: SeekFrom) -> std::io::Result<u64> {
            let SeekFrom::Start(position) = seek else {
                return Err(std::io::Error::other("only acquired absolute record seeks"));
            };
            self.position = position;
            Ok(position)
        }
    }
    fn image(endian: DafEndian) -> Vec<u8> {
        super::super::reader::build_image(
            endian,
            &[vec![SegSpec {
                start: 0.0,
                end: 86400.0,
                target: 499,
                center: 0,
                frame: 1,
                seg_type: 2,
                name: "DIRECTORY ONLY",
            }]],
        )
    }

    #[test]
    fn directory_probe_reads_only_three_records_even_with_distant_summaries() {
        for endian in [DafEndian::Little, DafEndian::Big] {
            let bytes = image(endian);
            let expected = SpkIndex::parse(&bytes).unwrap();
            let mut header: [u8; RECORD_BYTES] = bytes[..RECORD_BYTES].try_into().unwrap();
            let distant_record = 1_048_576;
            super::super::reader::put_i32(&mut header, OFF_FWARD, distant_record, endian);
            super::super::reader::put_i32(&mut header, OFF_BWARD, distant_record, endian);
            let mut source = DirectoryOnly {
                records: BTreeMap::from([
                    (0, header),
                    (
                        (distant_record - 1) as u64,
                        bytes[RECORD_BYTES..2 * RECORD_BYTES].try_into().unwrap(),
                    ),
                    (
                        distant_record as u64,
                        bytes[2 * RECORD_BYTES..3 * RECORD_BYTES]
                            .try_into()
                            .unwrap(),
                    ),
                ]),
                position: 0,
                bytes_read: 0,
            };
            let index = SpkIndex::read_directory_from(
                &mut source,
                (distant_record as usize + 1) * RECORD_BYTES,
            )
            .unwrap();
            assert_eq!(
                index.all_segments().collect::<Vec<_>>(),
                expected.all_segments().collect::<Vec<_>>()
            );
            assert_eq!(source.bytes_read, 3 * RECORD_BYTES);
        }
    }

    #[test]
    fn directory_probe_preserves_byte_parser_corruption_and_address_checks() {
        for endian in [DafEndian::Little, DafEndian::Big] {
            let original = image(endian);
            let mut variants = Vec::new();
            let mut count = original.clone();
            super::super::reader::put_f64(&mut count, RECORD_BYTES + 16, 26.0, endian);
            variants.push(count);
            let mut epoch = original.clone();
            super::super::reader::put_f64(&mut epoch, RECORD_BYTES + 24, f64::NAN, endian);
            variants.push(epoch);
            let mut address = original.clone();
            super::super::reader::put_i32(
                &mut address,
                RECORD_BYTES + 24 + 16 + 5 * 4,
                100_000,
                endian,
            );
            variants.push(address);
            variants.push(original[..2 * RECORD_BYTES].to_vec());
            for bytes in variants {
                let expected = SpkIndex::parse(&bytes).unwrap_err().to_string();
                let mut reader = std::io::Cursor::new(&bytes);
                let actual = SpkIndex::read_directory_from(&mut reader, bytes.len()).unwrap_err();
                assert_eq!(actual.to_string(), expected);
            }
        }
    }
}
