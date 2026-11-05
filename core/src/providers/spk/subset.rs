// SPDX-License-Identifier: AGPL-3.0
//! Streamed, lossless extraction of explicit J2000 Chebyshev SPK segments.

use super::{DafEndian, Descriptor, FileRecord, RECORD_BYTES, Segment, SpkIndex};
use sha2::{Digest, Sha256};
use std::fs::{File, OpenOptions};
use std::io::{self, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

/// One required source file and exact NAIF target/center pair.
///
/// Every intersecting type-2/3 J2000 segment is retained in source order.
/// Other targets are excluded; unsupported types or frames for this pair fail.
#[derive(Clone, Debug)]
pub struct SpkSubsetInput {
    /// Existing local NASA/JPL SPK; never modified or downloaded.
    pub path: PathBuf,
    /// Exact target NAIF identifier.
    pub target: i32,
    /// Exact center NAIF identifier.
    pub center: i32,
    /// Optional expected whole-file SHA-256. A mismatch prevents publication.
    pub expected_sha256: Option<String>,
}

/// Actual input identity and coverage of one successful SPK extraction.
#[derive(Debug)]
pub struct SpkSubsetReport {
    /// Source path and SHA-256, checked before and after extraction.
    pub sources: Vec<(PathBuf, String)>,
    /// Output segment descriptors in file priority order, including actual ET bounds.
    pub segments: Vec<Descriptor>,
    /// Published file length in bytes, including DAF record padding.
    pub bytes: u64,
}

struct Source {
    file: File,
    path: PathBuf,
    endian: DafEndian,
    index: SpkIndex,
    hash: String,
}

struct PlannedSegment {
    source: usize,
    segment: Segment,
    offset: u64,
    words: usize,
    trailer: [f64; 4],
}

fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

fn hash_file(file: &mut File) -> io::Result<String> {
    file.seek(SeekFrom::Start(0))?;
    let mut hasher = Sha256::new();
    let mut buffer = vec![0; 1024 * 1024];
    loop {
        let n = file.read(&mut buffer)?;
        if n == 0 {
            break;
        }
        hasher.update(&buffer[..n]);
    }
    Ok(crate::download::hex_of(&hasher.finalize()))
}

fn protect_source(output: &Path, source: &Source) -> io::Result<()> {
    let destination = match std::fs::metadata(output) {
        Ok(metadata) => metadata,
        Err(error) if error.kind() == io::ErrorKind::NotFound => return Ok(()),
        Err(error) => return Err(error),
    };
    let same_path = std::fs::canonicalize(output)? == std::fs::canonicalize(&source.path)?;
    #[cfg(unix)]
    let same_file = {
        use std::os::unix::fs::MetadataExt;
        let original = source.file.metadata()?;
        destination.dev() == original.dev() && destination.ino() == original.ino()
    };
    #[cfg(not(unix))]
    let same_file = {
        let _ = destination;
        false
    };
    if same_path || same_file {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "output aliases a source SPK",
        ));
    }
    Ok(())
}

fn plan_segment(
    source: &mut Source,
    mut segment: Segment,
    clip: Option<(f64, f64)>,
) -> io::Result<Option<PlannedSegment>> {
    let descriptor = &mut segment.descriptor;
    if descriptor.frame != 1 || !matches!(descriptor.seg_type, 2 | 3) {
        return Err(invalid(
            "SPK extraction requires type 2/3 and J2000 frame 1",
        ));
    }
    if let Some((start, end)) = clip {
        descriptor.start_et = descriptor.start_et.max(start);
        descriptor.end_et = descriptor.end_et.min(end);
    }
    if descriptor.start_et >= descriptor.end_et {
        return Ok(None);
    }
    let offset = (descriptor.initial_address as u64 - 1) * 8;
    let length = (descriptor.final_address - descriptor.initial_address + 1) as usize;
    if length < 4 {
        return Err(invalid(
            "SPK Chebyshev segment is shorter than its directory",
        ));
    }
    source
        .file
        .seek(SeekFrom::Start(offset + (length as u64 - 4) * 8))?;
    let mut raw = [0; 32];
    source.file.read_exact(&mut raw)?;
    let trailer = std::array::from_fn(|i| source.endian.read_f64(&raw[i * 8..i * 8 + 8]));
    let directory =
        super::decode::chebyshev_directory(descriptor.seg_type, length, offset as usize, trailer)
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
    let end = directory.init_et + directory.n_records as f64 * directory.int_len;
    if !end.is_finite() || descriptor.start_et < directory.init_et || descriptor.end_et > end {
        return Err(invalid("SPK descriptor exceeds its coefficient directory"));
    }
    // Keep the neighbouring record on either edge. This protects support at
    // representable joins when division rounds onto an adjacent record. Only
    // the descriptor is clipped; every retained coefficient word is unchanged.
    let first = (((descriptor.start_et - directory.init_et) / directory.int_len).floor() as usize)
        .saturating_sub(1)
        .min(directory.n_records - 1);
    let last = (((descriptor.end_et - directory.init_et) / directory.int_len).floor() as usize)
        .saturating_add(1)
        .min(directory.n_records - 1);
    let count = last - first + 1;
    Ok(Some(PlannedSegment {
        source: 0,
        segment,
        offset: offset + (first * directory.record_words) as u64 * 8,
        words: count * directory.record_words,
        trailer: [
            directory.init_et + first as f64 * directory.int_len,
            directory.int_len,
            directory.record_words as f64,
            count as f64,
        ],
    }))
}

fn check_coverage(plans: &[PlannedSegment], required: &[(i32, i32, f64, f64)]) -> io::Result<()> {
    for &(target, center, start, end) in required {
        if !start.is_finite() || !end.is_finite() || start >= end {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "invalid required ET interval",
            ));
        }
        let mut spans: Vec<_> = plans
            .iter()
            .map(|p| p.segment.descriptor)
            .filter(|d| d.target == target && d.center == center)
            .map(|d| (d.start_et, d.end_et))
            .collect();
        spans.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut covered = start;
        for (low, high) in spans {
            if low > covered {
                break;
            }
            covered = covered.max(high);
        }
        if covered < end {
            return Err(invalid(format!(
                "SPK target {target} relative to {center} has incomplete required coverage"
            )));
        }
    }
    Ok(())
}

struct Temporary(PathBuf);
impl Drop for Temporary {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

fn create_temporary(output: &Path) -> io::Result<(Temporary, File)> {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let name = output
        .file_name()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "output needs a filename"))?;
    for _ in 0..64 {
        let path = output.with_file_name(format!(
            ".{}-{}-{}.spk-tmp",
            name.to_string_lossy(),
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        match OpenOptions::new()
            .read(true)
            .write(true)
            .create_new(true)
            .open(&path)
        {
            Ok(file) => return Ok((Temporary(path), file)),
            Err(error) if error.kind() == io::ErrorKind::AlreadyExists => continue,
            Err(error) => return Err(error),
        }
    }
    Err(io::Error::new(
        io::ErrorKind::AlreadyExists,
        "cannot allocate an exclusive SPK temporary file",
    ))
}

fn write_directory(out: &mut File, plans: &mut [PlannedSegment]) -> io::Result<u64> {
    let groups = plans.len().div_ceil(super::MAX_SUMMARIES_PER_RECORD);
    let mut address = (1 + groups * 2) * super::WORDS_PER_RECORD + 1;
    for plan in plans.iter_mut() {
        plan.segment.descriptor.initial_address =
            i32::try_from(address).map_err(|_| invalid("SPK exceeds DAF address range"))?;
        address = address
            .checked_add(plan.words + 4)
            .ok_or_else(|| invalid("SPK address overflow"))?;
        plan.segment.descriptor.final_address =
            i32::try_from(address - 1).map_err(|_| invalid("SPK exceeds DAF address range"))?;
    }
    let free = i32::try_from(address).map_err(|_| invalid("SPK exceeds DAF address range"))?;
    let mut header = [0; RECORD_BYTES];
    header[..8].copy_from_slice(super::SPK_FILE_ID);
    header[8..12].copy_from_slice(&super::SPK_ND.to_le_bytes());
    header[12..16].copy_from_slice(&super::SPK_NI.to_le_bytes());
    header[16..76].fill(b' ');
    let name = b"LibEphemeris J2000 SPK subset";
    header[16..16 + name.len()].copy_from_slice(name);
    header[76..80].copy_from_slice(&2_i32.to_le_bytes());
    header[80..84].copy_from_slice(&((groups * 2) as i32).to_le_bytes());
    header[84..88].copy_from_slice(&free.to_le_bytes());
    header[88..96].copy_from_slice(b"LTL-IEEE");
    header[699..727].copy_from_slice(b"FTPSTR:\r:\n:\r\n:\r\0:\x81:\x10\xce:ENDFTP");
    out.write_all(&header)?;
    for (group, chunk) in plans.chunks(super::MAX_SUMMARIES_PER_RECORD).enumerate() {
        let mut summary = [0; RECORD_BYTES];
        let mut names = [b' '; RECORD_BYTES];
        let controls = [
            if group + 1 < groups {
                (4 + group * 2) as f64
            } else {
                0.0
            },
            if group > 0 { (group * 2) as f64 } else { 0.0 },
            chunk.len() as f64,
        ];
        for (index, value) in controls.into_iter().enumerate() {
            summary[index * 8..index * 8 + 8].copy_from_slice(&value.to_le_bytes());
        }
        for (index, plan) in chunk.iter().enumerate() {
            let d = plan.segment.descriptor;
            let start = 24 + index * 40;
            summary[start..start + 8].copy_from_slice(&d.start_et.to_le_bytes());
            summary[start + 8..start + 16].copy_from_slice(&d.end_et.to_le_bytes());
            for (i, value) in [
                d.target,
                d.center,
                d.frame,
                d.seg_type,
                d.initial_address,
                d.final_address,
            ]
            .into_iter()
            .enumerate()
            {
                summary[start + 16 + i * 4..start + 20 + i * 4]
                    .copy_from_slice(&value.to_le_bytes());
            }
            names[index * 40..(index + 1) * 40].copy_from_slice(&plan.segment.name);
        }
        out.write_all(&summary)?;
        out.write_all(&names)?;
    }
    Ok((address as u64 - 1) * 8)
}

/// Extract exact target/center pairs into an atomic little-endian DAF/SPK.
///
/// `clip_et` is an optional `(start, end)` clipping window in TDB seconds past
/// J2000. Each `required_et` row is `(NAIF target, NAIF center, start, end)` in
/// the same units; it requires descriptor coverage, not an SSB center chain.
/// Clipping retains
/// complete type-2/3 coefficient records, including edge neighbours; there is
/// no resampling, numerical fit or frame conversion. Input order preserves
/// overlapping-segment priority. Every source must contribute, and required
/// coverage checks the union without bridging gaps. SHA-256 is checked before
/// and after copying each source, with bounded streaming memory. Keep source
/// files stable during extraction: hash comparisons detect changes but do not
/// create an atomic filesystem snapshot. The output parent must already exist.
///
/// Errors: ordinary I/O failures; `InvalidInput` for invalid intervals or a
/// source/output alias; `InvalidData` for malformed, unsupported or incomplete
/// kernels. `overwrite=false` publishes without replacing an existing path,
/// including one created during extraction. Every failed attempt preserves the
/// destination and removes only its exclusively owned temporary file.
///
/// Format: [NASA/NAIF DAF](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html)
/// and SPK type-2/3 record directories. Other SPK types are not accepted here.
pub fn write_spk_subset(
    inputs: &[SpkSubsetInput],
    output: &Path,
    clip_et: Option<(f64, f64)>,
    required_et: &[(i32, i32, f64, f64)],
    overwrite: bool,
) -> io::Result<SpkSubsetReport> {
    if inputs.is_empty() || clip_et.is_some_and(|(a, b)| !a.is_finite() || !b.is_finite() || a >= b)
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid SPK subset request",
        ));
    }
    if !overwrite && output.symlink_metadata().is_ok() {
        return Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            "SPK output exists; use --force to replace it",
        ));
    }
    let mut sources: Vec<Source> = Vec::new();
    let mut plans = Vec::new();
    for input in inputs {
        let path = std::fs::canonicalize(&input.path)?;
        let source_index = match sources.iter().position(|source| source.path == path) {
            Some(index) => index,
            None => {
                let mut file = File::open(&path)?;
                let length = usize::try_from(file.metadata()?.len())
                    .map_err(|_| invalid("SPK is too large for this platform"))?;
                let mut header = [0; RECORD_BYTES];
                file.read_exact(&mut header)?;
                let endian = FileRecord::parse(&header)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?
                    .endian;
                let index = SpkIndex::read_directory_from(&mut file, length)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                let hash = hash_file(&mut file)?;
                sources.push(Source {
                    file,
                    path,
                    endian,
                    index,
                    hash,
                });
                sources.len() - 1
            }
        };
        let source = &mut sources[source_index];
        if input
            .expected_sha256
            .as_ref()
            .is_some_and(|expected| !source.hash.eq_ignore_ascii_case(expected))
        {
            return Err(invalid(format!(
                "unexpected SHA-256 for source {}",
                input.path.display()
            )));
        }
        protect_source(output, source)?;
        let selected: Vec<_> = source
            .index
            .segments_for(input.target)
            .filter(|s| s.descriptor.center == input.center)
            .cloned()
            .collect();
        let before = plans.len();
        for segment in selected {
            if let Some(mut plan) = plan_segment(source, segment, clip_et)? {
                plan.source = source_index;
                plans.push(plan);
            }
        }
        if plans.len() == before {
            return Err(invalid(format!(
                "source {} has no requested target/center coverage",
                input.path.display()
            )));
        }
    }
    check_coverage(&plans, required_et)?;
    let (temporary, mut out) = create_temporary(output)?;
    let length = write_directory(&mut out, &mut plans)?;
    let mut buffer = vec![0; 1024 * 1024];
    for plan in &plans {
        let source = &mut sources[plan.source];
        source.file.seek(SeekFrom::Start(plan.offset))?;
        let mut left = plan.words * 8;
        while left > 0 {
            let count = left.min(buffer.len());
            source.file.read_exact(&mut buffer[..count])?;
            if source.endian == DafEndian::Big {
                for word in buffer[..count].chunks_exact_mut(8) {
                    word.reverse();
                }
            }
            let first_word = plan.words - left / 8;
            let record_words = plan.trailer[2] as usize;
            for (index, word) in buffer[..count].chunks_exact(8).enumerate() {
                let value = DafEndian::Little.read_f64(word);
                if !value.is_finite() || ((first_word + index) % record_words == 1 && value <= 0.0)
                {
                    return Err(invalid(
                        "nonfinite SPK coefficient or invalid record radius",
                    ));
                }
            }
            out.write_all(&buffer[..count])?;
            left -= count;
        }
        for value in plan.trailer {
            out.write_all(&value.to_le_bytes())?;
        }
    }
    let padded = length.div_ceil(RECORD_BYTES as u64) * RECORD_BYTES as u64;
    out.set_len(padded)?;
    out.flush()?;
    out.sync_all()?;
    let output_length = usize::try_from(padded)
        .map_err(|_| invalid("SPK output is too large for this platform"))?;
    SpkIndex::read_directory_from(&mut out, output_length)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    drop(out);
    for source in &mut sources {
        if hash_file(&mut source.file)? != source.hash {
            return Err(invalid("SPK source changed during extraction"));
        }
        protect_source(output, source)?;
    }
    if overwrite {
        std::fs::rename(&temporary.0, output)?;
    } else {
        std::fs::hard_link(&temporary.0, output)?;
    }
    Ok(SpkSubsetReport {
        sources: sources.into_iter().map(|s| (s.path, s.hash)).collect(),
        segments: plans.into_iter().map(|p| p.segment.descriptor).collect(),
        bytes: padded,
    })
}
