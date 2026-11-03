// SPDX-License-Identifier: AGPL-3.0
//! Retained bounded-record access for immutable SPK2/3 force files.
//!
//! Directory decoding and coefficient mathematics are shared with the ordinary
//! owned-image provider. Only storage/record acquisition differs. No complete
//! coefficient image or growing epoch/result cache is retained. Callers own
//! authentication and immutability, and compose the declared center explicitly.

use super::decode::{ChebyshevDirectory, chebyshev_directory};
use super::evaluation::{
    eval_chebyshev_record_values, jd_parts_seconds, select_two_part_jd, split_seconds_covers,
};
use super::*;
use crate::types::TwoPartEpoch;
use std::collections::{BTreeSet, VecDeque};
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;
use std::sync::{Arc, Mutex};

const MAXIMUM_SEGMENTS: usize = 65536;
const MAXIMUM_RECORD_WORDS: usize = 1 << 20;
const CACHE_BYTES: usize = 4 * 1024 * 1024;
const CACHE_ENTRIES: usize = 128;

type CachedRecord = ((usize, usize), Arc<[f64]>);
struct RecordFile {
    file: File,
    records: VecDeque<CachedRecord>,
    /// Retained coefficient bytes, not allocator/RSS or in-flight old owners.
    bytes: usize,
}

/// Acquired immutable force-file owner. State evaluation returns native km/km/s
/// and the actual descriptor center; it never assumes that center is the Sun.
pub(crate) struct ChebyshevFile {
    index: SpkIndex,
    directories: Vec<ChebyshevDirectory>,
    endian: DafEndian,
    file: Mutex<RecordFile>,
}

impl ChebyshevFile {
    /// Acquire this exact local file and admit bounded J2000 type2/3 metadata.
    /// I/O, invalid directory/shape, unsupported frame/type and admission-limit
    /// errors are explicit. No replacement image, network request or mmap occurs.
    pub(crate) fn open(path: &Path) -> Result<Self> {
        let mut file = File::open(path).map_err(|error| {
            if error.kind() == std::io::ErrorKind::NotFound {
                Error::SpkFileNotFound {
                    path: path.to_string_lossy().into_owned(),
                }
            } else {
                corrupt("SPK force-file open", 0)
            }
        })?;
        let length = file
            .metadata()
            .map_err(|_| corrupt("SPK force-file length", 0))?
            .len();
        let length = usize::try_from(length).map_err(|_| corrupt("SPK force-file length", 0))?;
        let mut first = [0; RECORD_BYTES];
        file.read_exact(&mut first)
            .map_err(|_| corrupt("SPK force-file header", 0))?;
        let endian = FileRecord::parse(&first)?.endian;
        let index = SpkIndex::read_directory_bounded_from(&mut file, length, MAXIMUM_SEGMENTS)?;
        let mut directories = Vec::with_capacity(index.len());
        for segment in index.all_segments() {
            let descriptor = segment.descriptor;
            if descriptor.frame != FRAME_J2000 {
                return Err(Error::SpkUnsupportedFrame {
                    frame: descriptor.frame,
                });
            }
            if !matches!(
                descriptor.seg_type,
                TYPE_CHEBYSHEV_POS | TYPE_CHEBYSHEV_POS_VEL
            ) {
                return Err(Error::SpkUnsupportedType {
                    seg_type: descriptor.seg_type,
                });
            }
            let words = (descriptor.final_address - descriptor.initial_address + 1) as usize;
            if words < 4 {
                return Err(corrupt("SPK force segment length", 0));
            }
            let offset = (descriptor.initial_address as usize - 1) * 8;
            let trailer_offset = (descriptor.final_address as usize - 4) * 8;
            let mut trailer = [0; 32];
            file.seek(SeekFrom::Start(trailer_offset as u64))
                .map_err(|_| corrupt("SPK force directory seek", trailer_offset))?;
            file.read_exact(&mut trailer)
                .map_err(|_| corrupt("SPK force directory read", trailer_offset))?;
            let directory = chebyshev_directory(
                descriptor.seg_type,
                words,
                offset,
                std::array::from_fn(|i| endian.read_f64(&trailer[i * 8..i * 8 + 8])),
            )?;
            if directory.record_words > MAXIMUM_RECORD_WORDS
                || !(directory.init_et + directory.n_records as f64 * directory.int_len).is_finite()
            {
                return Err(corrupt("SPK force record/domain limit", trailer_offset));
            }
            directories.push(directory);
        }
        Ok(Self {
            index,
            directories,
            endian,
            file: Mutex::new(RecordFile {
                file,
                records: VecDeque::new(),
                bytes: 0,
            }),
        })
    }

    /// Unique targets in first file occurrence order; descriptor priority is
    /// still last-covering-segment-first during an actual state evaluation.
    pub(crate) fn targets(&self) -> Vec<i32> {
        let mut seen = BTreeSet::new();
        self.index
            .all_segments()
            .filter_map(|segment| {
                let target = segment.descriptor.target;
                seen.insert(target).then_some(target)
            })
            .collect()
    }

    /// Actual merged descriptor coverage in TDB seconds from J2000.
    pub(crate) fn coverage(&self, target: i32) -> Vec<(f64, f64)> {
        self.index.coverage(target)
    }

    /// Evaluate the last covering descriptor at the two-part TDB Julian day.
    /// Return its declared center and geometric km/km-s state; no center/frame
    /// conversion occurs. Missing targets, uncovered epochs, bad coefficients,
    /// I/O and nonfinite arithmetic propagate as typed SPK errors.
    pub(crate) fn state(&self, target: i32, jd: f64, fraction: f64) -> Result<(i32, ChebState)> {
        let segment = select_two_part_jd(&self.index, target, jd, fraction)?;
        let ordinal = self
            .index
            .all_segments()
            .position(|candidate| std::ptr::eq(candidate, segment))
            .expect("selected segment belongs to this retained index");
        let directory = &self.directories[ordinal];
        let (_, _, hi, lo) = jd_parts_seconds(jd, fraction)?;
        let time = TwoPartEpoch { hi, lo };
        let end = directory.init_et + directory.n_records as f64 * directory.int_len;
        if !split_seconds_covers(hi, lo, directory.init_et, end) {
            return Err(Error::SpkOutOfRange {
                target,
                et: hi + lo,
                start: directory.init_et,
                end,
            });
        }
        let relative = time.add_offset(-directory.init_et)?;
        let approximate = (relative.hi / directory.int_len).floor();
        let mut record = approximate.max(0.0).min((directory.n_records - 1) as f64) as usize;
        // Select with both epoch words. A tiny negative residual at a seam
        // must remain on the preceding record, even when its scalar offset
        // rounds up to the complete interval length.
        let start = directory.init_et + record as f64 * directory.int_len;
        if negative(time.add_offset(-start)?) && record > 0 {
            record -= 1;
        }
        let next = directory.init_et + (record + 1) as f64 * directory.int_len;
        if !negative(time.add_offset(-next)?) && record + 1 < directory.n_records {
            record += 1;
        }
        let words = self.record(ordinal, record)?;
        let midpoint = words[0];
        let radius = words[1];
        if !midpoint.is_finite() || !radius.is_finite() || radius <= 0.0 {
            return Err(corrupt("SPK force record domain", 0));
        }
        let delta = time.add_offset(-midpoint)?;
        if negative(delta.add_offset(radius)?) || positive(delta.add_offset(-radius)?) {
            return Err(Error::SpkOutOfRange {
                target,
                et: hi + lo,
                start: midpoint - radius,
                end: midpoint + radius,
            });
        }
        let tau = delta.hi / radius + delta.lo / radius;
        Ok((
            segment.descriptor.center,
            eval_chebyshev_record_values(segment.descriptor.seg_type, &words, tau, radius)?,
        ))
    }

    fn record(&self, segment: usize, record: usize) -> Result<Arc<[f64]>> {
        let key = (segment, record);
        let mut owner = self
            .file
            .lock()
            .map_err(|_| corrupt("SPK force-file owner poisoned", 0))?;
        if let Some(position) = owner
            .records
            .iter()
            .position(|(candidate, _)| *candidate == key)
        {
            let entry = owner
                .records
                .remove(position)
                .expect("located cached record");
            let result = Arc::clone(&entry.1);
            owner.records.push_back(entry);
            return Ok(result);
        }
        let descriptor = self
            .index
            .all_segments()
            .nth(segment)
            .expect("admitted segment ordinal")
            .descriptor;
        let directory = &self.directories[segment];
        if record >= directory.n_records {
            return Err(corrupt("SPK force record ordinal", 0));
        }
        let offset =
            (descriptor.initial_address as usize - 1 + record * directory.record_words) * 8;
        let size = directory.record_words * 8;
        let mut bytes = vec![0; size];
        owner
            .file
            .seek(SeekFrom::Start(offset as u64))
            .map_err(|_| corrupt("SPK force record seek", offset))?;
        owner
            .file
            .read_exact(&mut bytes)
            .map_err(|_| corrupt("SPK force record read", offset))?;
        let values: Vec<f64> = bytes
            .chunks_exact(8)
            .map(|word| self.endian.read_f64(word))
            .collect();
        let values: Arc<[f64]> = values.into();
        // A single admitted large record still works without being retained
        // beyond this borrow. Storage admission never substitutes mathematics.
        if size <= CACHE_BYTES {
            while owner.bytes + size > CACHE_BYTES || owner.records.len() >= CACHE_ENTRIES {
                let (_, old) = owner
                    .records
                    .pop_front()
                    .expect("nonempty bounded cache exceeds admission");
                owner.bytes -= old.len() * 8;
            }
            owner.bytes += size;
            owner.records.push_back((key, Arc::clone(&values)));
        }
        Ok(values)
    }
}

fn negative(value: TwoPartEpoch) -> bool {
    value.hi < 0.0 || (value.hi == 0.0 && value.lo < 0.0)
}
fn positive(value: TwoPartEpoch) -> bool {
    value.hi > 0.0 || (value.hi == 0.0 && value.lo > 0.0)
}
fn corrupt(what: &str, offset: usize) -> Error {
    Error::SpkCorrupt {
        what: what.into(),
        offset,
    }
}

#[cfg(test)]
#[path = "bounded_tests.rs"]
pub(super) mod tests;
