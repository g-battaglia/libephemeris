//! Streamed, lossless LEB1-to-LEB2 v2 conversion and checked publication.

use super::*;
use std::collections::BTreeSet;
use std::fs::{self, File, Metadata, OpenOptions};
use std::io::{Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

const COPY_BLOCK: usize = 64 * 1024;

/// Selection and bounded chunk geometry for lossless LEB2 v2 output.
/// No coefficient mantissa bits are discarded.
#[derive(Clone, Debug)]
pub struct Leb2ConvertOptions {
    /// Exact body identifiers, or all source bodies with `None`. Source order
    /// is retained; duplicate selections and absent identifiers are errors.
    pub body_ids: Option<Vec<i32>>,
    /// Copy every non-coefficient auxiliary section verbatim, including opaque
    /// section ids. Conflicting/duplicate reserved sections are refused.
    pub include_auxiliary: bool,
    /// Atomically replace an existing output. Input aliases are always refused.
    pub overwrite: bool,
    /// Nominal TT-grid chunk span in days, positive and finite both as `f64`
    /// and when stored as `f32`. A segment is never split; a single wider
    /// segment is retained if its coefficient bytes fit.
    pub chunk_interval_days: f64,
    /// Maximum decoded coefficient bytes per chunk, at least eight.
    /// A single larger segment is refused before reading its coefficients.
    /// Compression/decode scratch and index metadata are additional.
    pub max_chunk_bytes: usize,
}

impl Default for Leb2ConvertOptions {
    fn default() -> Self {
        Self {
            body_ids: None,
            include_auxiliary: true,
            overwrite: false,
            chunk_interval_days: 3652.5,
            max_chunk_bytes: 1024 * 1024,
        }
    }
}

/// Actual publication result; all reported chunks were decoded and compared
/// byte-for-byte with their original segment-major binary64 coefficient words.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Leb2ConversionReport {
    /// Selected bodies, in original source index order.
    pub body_ids: Vec<i32>,
    /// Actual published file length in bytes, including all directories.
    pub bytes_written: u64,
    /// Number of independently compressed/verified coefficient chunks.
    pub chunks: u64,
}

struct Plan {
    entry: BodyEntry,
    segment_bytes: usize,
    chunk_segments: u32,
    chunks: u32,
}

struct OutputSection {
    id: u32,
    reserved: [u8; 4],
    offset: u64,
    size: u64,
}

struct Temporary(PathBuf);
impl Drop for Temporary {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.0);
    }
}

fn corrupt(what: &str) -> Error {
    Error::LebCorrupt {
        what: format!("LEB conversion: {what}"),
        offset: 0,
    }
}

fn io_error(path: &Path, error: std::io::Error) -> Error {
    Error::EphemerisLoaderIo {
        path: path.display().to_string(),
        kind: error.kind(),
        errno: error.raw_os_error(),
        message: error.to_string(),
    }
}

fn invalid(message: &str) -> Error {
    Error::EphemerisLoaderValue {
        message: format!("LEB conversion: {message}"),
    }
}

fn add(a: u64, b: u64) -> Result<u64> {
    a.checked_add(b)
        .ok_or_else(|| corrupt("byte offset overflow"))
}

fn length(bytes: u64) -> Result<usize> {
    usize::try_from(bytes).map_err(|_| corrupt("range exceeds platform address space"))
}

fn same_inode(a: &Metadata, b: &Metadata) -> bool {
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        a.dev() == b.dev() && a.ino() == b.ino()
    }
    #[cfg(not(unix))]
    {
        let _ = (a, b);
        false
    }
}

fn protect_input(source: &Path, identity: &Metadata, output: &Path) -> Result<()> {
    if source == output {
        return Err(invalid("output aliases the input"));
    }
    match fs::metadata(output) {
        Ok(other) => {
            let canonical = fs::canonicalize(output).map_err(|e| io_error(output, e))?;
            if canonical == source || same_inode(identity, &other) {
                return Err(invalid("output aliases the input"));
            }
        }
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
        Err(error) => return Err(io_error(output, error)),
    }
    Ok(())
}

fn check_source(input: &LebFile) -> Result<()> {
    match input.source()? {
        LebSource::Positional(file) => file.validate_identity(),
        LebSource::Slice(_) => Err(corrupt("conversion requires an indexed source")),
    }
}

struct Target {
    source: PathBuf,
    output: PathBuf,
}

/// Resolve existing ancestors without creating the planned output directory.
fn output_parent(parent: &Path) -> Result<PathBuf> {
    let mut cursor = std::path::absolute(parent).map_err(|e| io_error(parent, e))?;
    let mut missing = Vec::new();
    loop {
        match fs::canonicalize(&cursor) {
            Ok(mut resolved) => {
                for part in missing.into_iter().rev() {
                    if part == std::ffi::OsStr::new("..") {
                        resolved.pop();
                    } else if part != std::ffi::OsStr::new(".") {
                        resolved.push(part);
                    }
                }
                return Ok(resolved);
            }
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
                // A dangling symlink is not an absent directory to create.
                match fs::symlink_metadata(&cursor) {
                    Ok(_) => return Err(io_error(&cursor, error)),
                    Err(other) if other.kind() == std::io::ErrorKind::NotFound => {}
                    Err(other) => return Err(io_error(&cursor, other)),
                }
                let part = match cursor.components().next_back() {
                    Some(std::path::Component::Normal(part)) => part.to_owned(),
                    Some(std::path::Component::ParentDir) => std::ffi::OsString::from(".."),
                    _ => return Err(io_error(&cursor, error)),
                };
                missing.push(part);
                if !cursor.pop() {
                    return Err(io_error(parent, error));
                }
            }
            Err(error) => return Err(io_error(&cursor, error)),
        }
    }
}

fn target(source: &Path, output: &Path, overwrite: bool) -> Result<Target> {
    let source = fs::canonicalize(source).map_err(|e| io_error(source, e))?;
    let identity = fs::metadata(&source).map_err(|e| io_error(&source, e))?;
    if !identity.is_file() {
        return Err(invalid("input must be a regular file"));
    }
    let parent = output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let parent = output_parent(parent)?;
    let output = parent.join(
        output
            .file_name()
            .ok_or_else(|| invalid("output needs a filename"))?,
    );
    protect_input(&source, &identity, &output)?;
    match fs::symlink_metadata(&output) {
        Ok(_) if !overwrite => {
            return Err(io_error(
                &output,
                std::io::Error::new(
                    std::io::ErrorKind::AlreadyExists,
                    "LEB output exists; use overwrite explicitly",
                ),
            ));
        }
        Ok(_) => {}
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
        Err(error) => return Err(io_error(&output, error)),
    }
    Ok(Target { source, output })
}

/// Preflight one output without creating or deleting any filesystem entry.
/// The converter independently reacquires/rechecks its actual input and target
/// before publication; this preliminary check is not a reservation.
pub(crate) fn validate_leb2_conversion_target(
    source: &Path,
    output: &Path,
    overwrite: bool,
) -> Result<()> {
    target(source, output, overwrite).map(|_| ())
}

fn publish(temporary: &Path, output: &Path, overwrite: bool) -> Result<()> {
    if overwrite {
        fs::rename(temporary, output).map_err(|e| io_error(output, e))
    } else {
        fs::hard_link(temporary, output).map_err(|e| io_error(output, e))
    }
}

fn temporary(parent: &Path) -> Result<(Temporary, File)> {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    for _ in 0..64 {
        let path = parent.join(format!(
            ".leb2-convert-{}-{}.tmp",
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
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => return Err(io_error(&path, error)),
        }
    }
    Err(invalid("could not create an exclusive temporary output"))
}

fn plans(input: &LebFile, options: &Leb2ConvertOptions) -> Result<Vec<Plan>> {
    let index = input.index()?;
    if index.header.format != LebFormat::Leb1 {
        return Err(invalid("input must be LEB1 version 1"));
    }
    if index.header.body_count as usize != index.bodies.len()
        || index.header.section_count as usize != index.sections.len()
        || index.section(SECTION_COMPRESSED_CHEBYSHEV).is_some()
    {
        return Err(corrupt(
            "duplicate body/section or conflicting compressed section",
        ));
    }
    let requested = options
        .body_ids
        .as_ref()
        .map(|ids| ids.iter().copied().collect::<BTreeSet<_>>());
    if let (Some(ids), Some(unique)) = (&options.body_ids, &requested) {
        if ids.len() != unique.len() {
            return Err(invalid("duplicate body selection"));
        }
        for id in ids {
            if !index.has_body(*id) {
                return Err(Error::LebBodyNotFound { body_id: *id });
            }
        }
    }
    let mut out = Vec::new();
    for entry in &index.bodies {
        if requested
            .as_ref()
            .is_some_and(|ids| !ids.contains(&entry.body_id))
        {
            continue;
        }
        let total = entry.expected_uncompressed_size()?;
        let segment_bytes = length(total / u64::from(entry.segment_count))?;
        if segment_bytes > options.max_chunk_bytes {
            return Err(Error::ResourceLimit {
                resource: "LEB conversion segment bytes",
                required: segment_bytes,
                limit: options.max_chunk_bytes,
            });
        }
        let last_start = entry.jd_start + f64::from(entry.segment_count - 1) * entry.interval_days;
        let end = entry.jd_start + f64::from(entry.segment_count) * entry.interval_days;
        if !end.is_finite() || entry.jd_end <= last_start || entry.jd_end > end {
            return Err(corrupt("body domain does not match its stored segments"));
        }
        let data_end = add(entry.data_offset, total)?;
        if data_end > input.len() as u64 {
            return Err(corrupt("truncated coefficient range"));
        }
        let prefix = add(
            HEADER_SIZE as u64,
            u64::from(index.header.section_count) * SECTION_DIR_SIZE as u64,
        )?;
        if entry.data_offset < prefix {
            return Err(corrupt("coefficients overlap the file header/directory"));
        }
        for section in &index.sections {
            let section_end = add(section.offset, section.size)?;
            if section.section_id != SECTION_CHEBYSHEV
                && section.size != 0
                && entry.data_offset < section_end
                && section.offset < data_end
            {
                return Err(corrupt("coefficients overlap an index/auxiliary section"));
            }
        }
        if let Some(raw) = index.section(SECTION_CHEBYSHEV) {
            if entry.data_offset < raw.offset || data_end > add(raw.offset, raw.size)? {
                return Err(corrupt("coefficients exceed the raw coefficient section"));
            }
        }
        let by_time = (options.chunk_interval_days / entry.interval_days)
            .floor()
            .max(1.0)
            .min(f64::from(entry.segment_count)) as u32;
        let by_bytes =
            (options.max_chunk_bytes / segment_bytes).min(entry.segment_count as usize) as u32;
        let chunk_segments = by_time.min(by_bytes);
        out.push(Plan {
            entry: entry.clone(),
            segment_bytes,
            chunk_segments,
            chunks: entry.segment_count.div_ceil(chunk_segments),
        });
    }
    Ok(out)
}

fn write_at(file: &mut File, path: &Path, offset: u64, bytes: &[u8]) -> Result<()> {
    file.seek(SeekFrom::Start(offset))
        .map_err(|e| io_error(path, e))?;
    file.write_all(bytes).map_err(|e| io_error(path, e))
}

fn position(file: &mut File, path: &Path) -> Result<u64> {
    file.stream_position().map_err(|e| io_error(path, e))
}

fn raw_chunk(input: &LebFile, plan: &Plan, first: u32, count: u32) -> Result<Vec<u8>> {
    let offset = add(
        plan.entry.data_offset,
        u64::from(first) * plan.segment_bytes as u64,
    )?;
    let bytes = input
        .source()?
        .range(
            length(offset)?,
            count as usize * plan.segment_bytes,
            "LEB conversion coefficients",
        )?
        .into_owned();
    for word in bytes.chunks_exact(8) {
        let value = f64::from_le_bytes(word.try_into().expect("eight-byte coefficient"));
        if !value.is_finite() {
            return Err(corrupt("nonfinite coefficient"));
        }
    }
    Ok(bytes)
}

fn write_body(
    file: &mut File,
    path: &Path,
    input: &LebFile,
    plan: &Plan,
    options: &Leb2ConvertOptions,
    max_read: &mut usize,
) -> Result<BodyEntry> {
    let base = position(file, path)?;
    let mut header = [0; CHUNK_INDEX_HEADER_SIZE];
    header[..4].copy_from_slice(&plan.chunks.to_le_bytes());
    header[8..12].copy_from_slice(&(options.chunk_interval_days as f32).to_le_bytes());
    file.write_all(&header).map_err(|e| io_error(path, e))?;
    let directory_bytes = u64::from(plan.chunks) * CHUNK_ENTRY_SIZE as u64;
    *max_read = (*max_read).max(length(directory_bytes)?);
    let mut end = add(add(base, CHUNK_INDEX_HEADER_SIZE as u64)?, directory_bytes)?;
    file.set_len(end).map_err(|e| io_error(path, e))?;
    for ordinal in 0..plan.chunks {
        let first = ordinal * plan.chunk_segments;
        let count = (plan.entry.segment_count - first).min(plan.chunk_segments);
        let raw = raw_chunk(input, plan, first, count)?;
        let reordered =
            reorder_to_coeff_major(&raw, count, plan.entry.degree, plan.entry.components)?;
        let shuffled = shuffle_bytes(&reordered)?;
        let compressed = zstd::encode_all(shuffled.as_slice(), 3).map_err(|e| io_error(path, e))?;
        *max_read = (*max_read).max(compressed.len());
        write_at(file, path, end, &compressed)?;
        let start_jd = plan.entry.jd_start + f64::from(first) * plan.entry.interval_days;
        let end_jd = if first + count == plan.entry.segment_count {
            plan.entry.jd_end
        } else {
            plan.entry.jd_start + f64::from(first + count) * plan.entry.interval_days
        };
        let mut record = [0; CHUNK_ENTRY_SIZE];
        record[..8].copy_from_slice(&start_jd.to_le_bytes());
        record[8..16].copy_from_slice(&end_jd.to_le_bytes());
        record[16..24].copy_from_slice(&end.to_le_bytes());
        record[24..32].copy_from_slice(&(compressed.len() as u64).to_le_bytes());
        record[32..40].copy_from_slice(&(raw.len() as u64).to_le_bytes());
        record[40..44].copy_from_slice(&first.to_le_bytes());
        record[44..48].copy_from_slice(&count.to_le_bytes());
        let at = add(
            add(base, CHUNK_INDEX_HEADER_SIZE as u64)?,
            u64::from(ordinal) * CHUNK_ENTRY_SIZE as u64,
        )?;
        write_at(file, path, at, &record)?;
        end = add(end, compressed.len() as u64)?;
    }
    file.seek(SeekFrom::Start(end))
        .map_err(|e| io_error(path, e))?;
    let mut entry = plan.entry.clone();
    entry.data_offset = base;
    entry.compressed_size = Some(end - base);
    entry.uncompressed_size = Some(entry.expected_uncompressed_size()?);
    Ok(entry)
}

fn body_bytes(entry: &BodyEntry) -> Vec<u8> {
    let mut out = Vec::with_capacity(LEB2_BODY_ENTRY_SIZE);
    out.extend_from_slice(&entry.body_id.to_le_bytes());
    out.extend_from_slice(&entry.coord_type.to_le_bytes());
    out.extend_from_slice(&entry.segment_count.to_le_bytes());
    out.extend_from_slice(&entry.jd_start.to_le_bytes());
    out.extend_from_slice(&entry.jd_end.to_le_bytes());
    out.extend_from_slice(&entry.interval_days.to_le_bytes());
    out.extend_from_slice(&entry.degree.to_le_bytes());
    out.extend_from_slice(&entry.components.to_le_bytes());
    out.extend_from_slice(&entry.data_offset.to_le_bytes());
    out.extend_from_slice(&entry.compressed_size.expect("converted body").to_le_bytes());
    out.extend_from_slice(
        &entry
            .uncompressed_size
            .expect("converted body")
            .to_le_bytes(),
    );
    out
}

fn copy_auxiliary(
    input: &LebFile,
    output: &mut File,
    path: &Path,
    section: &SectionEntry,
) -> Result<()> {
    let mut copied = 0;
    while copied < section.size {
        let count = (section.size - copied).min(COPY_BLOCK as u64);
        let bytes = input.source()?.range(
            length(add(section.offset, copied)?)?,
            length(count)?,
            "LEB conversion auxiliary section",
        )?;
        output.write_all(&bytes).map_err(|e| io_error(path, e))?;
        copied += count;
    }
    Ok(())
}

fn verify_output(
    input: &LebFile,
    path: &Path,
    plans: &[Plan],
    sections: &[OutputSection],
    limits: IndexedLimits,
) -> Result<()> {
    let output = LebFile::open_indexed(path, limits)?;
    let index = output.index()?;
    if index.header.format != LebFormat::Leb2V2 || index.bodies.len() != plans.len() {
        return Err(corrupt("output format/body count"));
    }
    for plan in plans {
        let entry = index
            .bodies
            .iter()
            .find(|entry| entry.body_id == plan.entry.body_id)
            .ok_or_else(|| corrupt("output body missing"))?;
        let (_, chunks) =
            decode::decode_chunk_index_source(output.source()?, entry.data_offset, entry.body_id)?;
        if chunks.values.len() != plan.chunks as usize {
            return Err(corrupt("output chunk count"));
        }
        let mut first = 0;
        for (ordinal, chunk) in chunks.values.iter().enumerate() {
            if chunk.segment_start != first || chunk.uncompressed_size > limits.frame_bytes as u64 {
                return Err(corrupt("output chunk geometry"));
            }
            let compressed = output.source()?.range(
                length(chunk.blob_offset)?,
                length(chunk.compressed_size)?,
                "LEB conversion verification",
            )?;
            let decoded = decompress_chunk(
                &compressed,
                chunk,
                entry.degree,
                entry.components,
                entry.body_id,
                ordinal,
            )?;
            if decoded != raw_chunk(input, plan, first, chunk.segment_count)? {
                return Err(corrupt("output coefficient words differ"));
            }
            first += chunk.segment_count;
        }
        if first != plan.entry.segment_count {
            return Err(corrupt("output does not cover every segment"));
        }
    }
    for section in sections
        .iter()
        .filter(|s| s.id != SECTION_BODY_INDEX && s.id != SECTION_COMPRESSED_CHEBYSHEV)
    {
        let old = input
            .index()?
            .section(section.id)
            .ok_or_else(|| corrupt("auxiliary missing"))?;
        let mut at = 0;
        while at < section.size {
            let count = (section.size - at).min(COPY_BLOCK as u64);
            let a = input.source()?.range(
                length(add(old.offset, at)?)?,
                length(count)?,
                "LEB auxiliary verification",
            )?;
            let b = output.source()?.range(
                length(add(section.offset, at)?)?,
                length(count)?,
                "LEB auxiliary verification",
            )?;
            if a != b {
                return Err(corrupt("auxiliary bytes differ"));
            }
            at += count;
        }
    }
    Ok(())
}

/// Convert an immutable regular LEB1 file to lossless chunked LEB2 v2.
///
/// The caller must keep the acquired source inode immutable for the whole
/// operation (the indexed-reader contract). Atomic pathname replacement after
/// acquisition is permitted; both readers retain the original inode. This is
/// not an arbitrary-file snapshot. Each decoded chunk
/// is bounded by `max_chunk_bytes`; reorder/shuffle/zstd/verification scratch
/// adds a constant number of chunk buffers and a decoder window. Parsed index
/// and one body's chunk directory use the ordinary indexed metadata budget.
/// Auxiliary sections are copied/compared in blocks of at most 64 KiB. No full-file
/// coefficient buffer, mantissa truncation or astronomical evaluation is used.
/// Missing output parents are created only after input/selection planning.
///
/// Invalid options, ambiguous/overlapping tables, nonfinite coefficients,
/// truncated data, absent selected bodies, resource limits and native codec/I/O
/// errors prevent publication. Existing output survives failure; default
/// publication uses an exclusive hard link, so a raced-in output is not replaced.
/// `overwrite` uses same-directory atomic rename. Canonical/hard-link input
/// aliases are always refused. Only the newly created temporary is cleaned up.
pub fn convert_leb1_to_leb2(
    source: &Path,
    output: &Path,
    options: &Leb2ConvertOptions,
) -> Result<Leb2ConversionReport> {
    if !options.chunk_interval_days.is_finite()
        || options.chunk_interval_days <= 0.0
        || !(options.chunk_interval_days as f32).is_finite()
        || options.chunk_interval_days as f32 <= 0.0
        || options.max_chunk_bytes < 8
    {
        return Err(invalid(
            "positive finite f32 chunk interval and at least8 chunk bytes required",
        ));
    }
    let Target { source, output } = target(source, output, options.overwrite)?;
    let parent = output.parent().expect("canonical output parent");
    let mut limits = IndexedLimits {
        cache_bytes: 0,
        compressed_bytes: options.max_chunk_bytes.max(COPY_BLOCK),
        frame_bytes: options.max_chunk_bytes,
        ..IndexedLimits::default()
    };
    let input = LebFile::open_indexed(&source, limits)?;
    // Alias protection uses metadata of the very same descriptor supplying
    // index and coefficient reads, not a second pathname acquisition.
    let LebSource::Positional(acquired) = input.source()? else {
        return Err(corrupt("conversion requires an indexed source"));
    };
    let identity = acquired.acquired_metadata().clone();
    check_source(&input)?;
    protect_input(&source, &identity, &output)?;
    let plans = plans(&input, options)?;
    let auxiliaries: Vec<_> = input
        .index()?
        .sections
        .iter()
        .filter(|section| {
            options.include_auxiliary
                && !matches!(section.section_id, SECTION_BODY_INDEX | SECTION_CHEBYSHEV)
        })
        .cloned()
        .collect();
    let section_count =
        u32::try_from(auxiliaries.len() + 2).map_err(|_| corrupt("section count overflow"))?;
    let index_size = plans
        .len()
        .checked_mul(LEB2_BODY_ENTRY_SIZE)
        .ok_or_else(|| corrupt("body index size overflow"))? as u64;
    let index_offset = add(
        HEADER_SIZE as u64,
        u64::from(section_count) * SECTION_DIR_SIZE as u64,
    )?;
    let data_offset = add(index_offset, index_size)?;
    // Batch callers can preflight every destination before creating any parent.
    fs::create_dir_all(parent).map_err(|e| io_error(parent, e))?;
    let (temporary, mut file) = temporary(parent)?;
    let path = &temporary.0;
    file.set_len(data_offset).map_err(|e| io_error(path, e))?;
    file.seek(SeekFrom::Start(data_offset))
        .map_err(|e| io_error(path, e))?;
    let mut max_read = COPY_BLOCK.max(options.max_chunk_bytes);
    for (ordinal, plan) in plans.iter().enumerate() {
        let entry = write_body(&mut file, path, &input, plan, options, &mut max_read)?;
        let end = position(&mut file, path)?;
        write_at(
            &mut file,
            path,
            add(index_offset, ordinal as u64 * LEB2_BODY_ENTRY_SIZE as u64)?,
            &body_bytes(&entry),
        )?;
        file.seek(SeekFrom::Start(end))
            .map_err(|e| io_error(path, e))?;
    }
    let compressed_end = position(&mut file, path)?;
    let mut sections = vec![
        OutputSection {
            id: SECTION_BODY_INDEX,
            reserved: [0; 4],
            offset: index_offset,
            size: index_size,
        },
        OutputSection {
            id: SECTION_COMPRESSED_CHEBYSHEV,
            reserved: [0; 4],
            offset: data_offset,
            size: compressed_end - data_offset,
        },
    ];
    for section in auxiliaries {
        let at = position(&mut file, path)?;
        copy_auxiliary(&input, &mut file, path, &section)?;
        let ordinal = input
            .index()?
            .sections
            .iter()
            .position(|s| s.section_id == section.section_id)
            .ok_or_else(|| corrupt("auxiliary directory missing"))?;
        let raw = input.source()?.range(
            HEADER_SIZE + ordinal * SECTION_DIR_SIZE,
            SECTION_DIR_SIZE,
            "LEB conversion section metadata",
        )?;
        sections.push(OutputSection {
            id: section.section_id,
            reserved: raw[4..8].try_into().expect("reserved section word"),
            offset: at,
            size: section.size,
        });
    }
    let bytes_written = position(&mut file, path)?;
    let mut header = input
        .source()?
        .range(0, HEADER_SIZE, "LEB conversion header")?
        .into_owned();
    header[..4].copy_from_slice(&MAGIC_LEB2);
    header[4..8].copy_from_slice(&VERSION_LEB2_V2.to_le_bytes());
    header[8..12].copy_from_slice(&section_count.to_le_bytes());
    header[12..16].copy_from_slice(&(plans.len() as u32).to_le_bytes());
    header[40..44].copy_from_slice(&COMPRESSION_ZSTD_TRUNC_SHUFFLE.to_le_bytes());
    write_at(&mut file, path, 0, &header)?;
    for (ordinal, section) in sections.iter().enumerate() {
        let mut raw = Vec::with_capacity(SECTION_DIR_SIZE);
        raw.extend_from_slice(&section.id.to_le_bytes());
        raw.extend_from_slice(&section.reserved);
        raw.extend_from_slice(&section.offset.to_le_bytes());
        raw.extend_from_slice(&section.size.to_le_bytes());
        write_at(
            &mut file,
            path,
            HEADER_SIZE as u64 + ordinal as u64 * SECTION_DIR_SIZE as u64,
            &raw,
        )?;
    }
    file.sync_all().map_err(|e| io_error(path, e))?;
    limits.compressed_bytes = max_read;
    verify_output(&input, path, &plans, &sections, limits)?;
    check_source(&input)?;
    protect_input(&source, &identity, &output)?;
    drop(file);
    publish(path, &output, options.overwrite)?;
    Ok(Leb2ConversionReport {
        body_ids: plans.iter().map(|p| p.entry.body_id).collect(),
        bytes_written,
        chunks: plans.iter().map(|p| u64::from(p.chunks)).sum(),
    })
}

#[cfg(test)]
#[path = "tests/conversion.rs"]
mod tests;
