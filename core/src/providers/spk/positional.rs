//! Bounded record reads for immutable SPK images.
//!
//! The caller keeps the acquired file immutable. One descriptor index and an
//! open file are retained; coefficient payloads are read only for the selected
//! record. Owned and positional paths share record selection and evaluation.

use super::decode::chebyshev_directory;
use super::evaluation::{
    chebyshev_record_coordinates, eval_selected_chebyshev_record, select_two_part_jd,
};
use super::{ChebState, DafEndian, Descriptor, Error, FileRecord, Result, SpkIndex};
use std::fs::{File, Metadata};
use std::io::{Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::{Mutex, OnceLock};

#[cfg(feature = "nbody")]
pub(crate) struct PositionalChebyshev {
    file: File,
    endian: DafEndian,
    index: SpkIndex,
}

#[cfg(feature = "nbody")]
impl PositionalChebyshev {
    /// Acquire one immutable regular SPK without copying its payload.
    pub(crate) fn open(path: &Path) -> Result<Self> {
        let mut file = open(path)?;
        let identity = metadata(&file)?;
        let (endian, index) = read_index(&mut file, &identity)?;
        Ok(Self {
            file,
            endian,
            index,
        })
    }

    pub(crate) fn contains(&self, target: i32) -> bool {
        self.index.segments_for(target).next().is_some()
    }

    /// Only a covered Type2 record is eligible for the private trajectory.
    /// Other types retain the caller's existing integration policy. Record or
    /// coefficient failures in an eligible Type2 source must not decline it.
    pub(crate) fn trajectory_covers(&self, target: i32, jd: f64, fraction: f64) -> Result<bool> {
        super::resolve_chain(&self.index, target)?;
        match select_two_part_jd(&self.index, target, jd, fraction) {
            Ok(segment) => Ok(segment.descriptor.seg_type == super::TYPE_CHEBYSHEV_POS),
            Err(Error::SpkOutOfRange { .. }) => Ok(false),
            Err(error) => Err(error),
        }
    }

    /// A trajectory requires dP/dTDB. Type2 derives that from its position
    /// polynomial; Type3's separately stored velocity need not be its tangent.
    /// Keep the private trajectory policy explicit until both quantities can
    /// be carried separately through an owned sampler.
    pub(crate) fn trajectory_state(
        &mut self,
        target: i32,
        jd: f64,
        fraction: f64,
    ) -> Result<(ChebState, Descriptor)> {
        let (state, descriptor) = self.relative_state(target, jd, fraction)?;
        if descriptor.seg_type != super::TYPE_CHEBYSHEV_POS {
            return Err(Error::SpkUnsupportedType {
                seg_type: descriptor.seg_type,
            });
        }
        Ok((state, descriptor))
    }

    /// Center-relative J2000 state at split TDB JD, in km and km/s.
    /// Missing coverage, unsupported metadata and read errors are explicit.
    pub(crate) fn relative_state(
        &mut self,
        target: i32,
        jd: f64,
        fraction: f64,
    ) -> Result<(ChebState, Descriptor)> {
        super::resolve_chain(&self.index, target)?;
        relative_state(&self.index, target, jd, fraction, |offset, count| {
            read_words::<Vec<f64>>(&mut self.file, self.endian, offset, count)
        })
    }
}

fn open(path: &Path) -> Result<File> {
    crate::leb::open_regular_file(path).map_err(|error| {
        if error.kind() == std::io::ErrorKind::NotFound {
            Error::SpkFileNotFound {
                path: path.display().to_string(),
            }
        } else {
            corrupt("SPK file read", 0)
        }
    })
}

fn metadata(file: &File) -> Result<Metadata> {
    file.metadata().map_err(|_| corrupt("SPK file length", 0))
}

fn read_index(file: &mut File, identity: &Metadata) -> Result<(DafEndian, SpkIndex)> {
    let length = usize::try_from(identity.len()).map_err(|_| corrupt("SPK file length", 0))?;
    let mut header = [0; super::RECORD_BYTES];
    file.seek(SeekFrom::Start(0))
        .and_then(|_| file.read_exact(&mut header))
        .map_err(|_| corrupt("SPK file record", 0))?;
    let endian = FileRecord::parse(&header)?.endian;
    let index = SpkIndex::read_directory_from(file, length)?;
    Ok((endian, index))
}

fn relative_state<W: AsRef<[f64]>>(
    index: &SpkIndex,
    target: i32,
    jd: f64,
    fraction: f64,
    mut read: impl FnMut(usize, usize) -> Result<W>,
) -> Result<(ChebState, Descriptor)> {
    let descriptor = select_two_part_jd(index, target, jd, fraction)?.descriptor;
    if descriptor.frame != super::FRAME_J2000 {
        return Err(Error::SpkUnsupportedFrame {
            frame: descriptor.frame,
        });
    }
    if !matches!(
        descriptor.seg_type,
        super::TYPE_CHEBYSHEV_POS | super::TYPE_CHEBYSHEV_POS_VEL
    ) {
        return Err(Error::SpkUnsupportedType {
            seg_type: descriptor.seg_type,
        });
    }
    let start = (descriptor.initial_address as usize - 1) * 8;
    let length = (descriptor.final_address - descriptor.initial_address + 1) as usize;
    if length < 4 {
        return Err(corrupt("SPK segment length", start));
    }
    let trailer = read(start + (length - 4) * 8, 4)?;
    let directory = chebyshev_directory(
        descriptor.seg_type,
        length,
        start,
        trailer.as_ref().try_into().expect("four directory words"),
    )?;
    let end = directory.init_et + directory.n_records as f64 * directory.int_len;
    if !end.is_finite() {
        return Err(corrupt("SPK Chebyshev directory", start));
    }
    let (index, offset, _interior, hi, lo) = chebyshev_record_coordinates(
        directory.init_et,
        directory.int_len,
        directory.n_records,
        jd,
        fraction,
    )?;
    let record_offset = start + index * directory.record_words * 8;
    let record = read(record_offset, directory.record_words)?;
    let state = eval_selected_chebyshev_record(
        &directory,
        descriptor.seg_type,
        record.as_ref(),
        index,
        offset,
        record_offset,
        Some((hi, lo)),
    )?
    .0;
    Ok((state, descriptor))
}

fn read_words<W: FromIterator<f64>>(
    file: &mut (impl Read + Seek),
    endian: DafEndian,
    offset: usize,
    count: usize,
) -> Result<W> {
    let length = count
        .checked_mul(8)
        .ok_or_else(|| corrupt("SPK record length", offset))?;
    let mut bytes = vec![0; length];
    file.seek(SeekFrom::Start(offset as u64))
        .map_err(|_| corrupt("SPK record seek", offset))?;
    file.read_exact(&mut bytes)
        .map_err(|_| corrupt("SPK record read", offset))?;
    Ok(bytes
        .chunks_exact(8)
        .map(|word| endian.read_f64(word))
        .collect())
}

/// Lazily opens the first admissible pinned planet-center image. The opened
/// descriptor remains authoritative across pathname replacement.
#[derive(Debug)]
pub(crate) struct PositionalCenters {
    candidates: Vec<(PathBuf, &'static str)>,
    owner: OnceLock<Option<VerifiedImage>>,
}

impl PositionalCenters {
    pub(crate) fn new(candidates: Vec<(PathBuf, &'static str)>) -> Self {
        Self {
            candidates,
            owner: OnceLock::new(),
        }
    }

    pub(crate) fn relative_state(
        &self,
        target: i32,
        jd: f64,
        fraction: f64,
    ) -> Result<(ChebState, Descriptor)> {
        let image = self
            .owner
            .get_or_init(|| {
                self.candidates.iter().find_map(|(path, pin)| {
                    if !path.is_file() {
                        return None;
                    }
                    match VerifiedImage::open(path, pin) {
                        Ok(image) => Some(image),
                        Err(error) => {
                            log::warn!("Ignoring planet-center image {}: {error}", path.display());
                            None
                        }
                    }
                })
            })
            .as_ref()
            .ok_or(Error::SpkBodyNotFound { target })?;
        image.relative_state(target, jd, fraction)
    }
}

/// Pinned data directories replace files by rename. Size and modification
/// time guards detect ordinary in-place writes even on cached record reads.
struct VerifiedImage {
    file: File,
    cursor: Mutex<()>,
    identity: Metadata,
    endian: DafEndian,
    index: SpkIndex,
    cache: super::records::RecordCache,
}

impl std::fmt::Debug for VerifiedImage {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("VerifiedImage")
            .field("bytes", &self.identity.len())
            .field("descriptors", &self.index.len())
            .finish()
    }
}

impl VerifiedImage {
    fn open(path: &Path, pin: &str) -> Result<Self> {
        use sha2::{Digest, Sha256};

        let mut file = open(path)?;
        let identity = metadata(&file)?;
        let mut hash = Sha256::new();
        let mut buffer = vec![0; 64 * 1024];
        let mut remaining = identity.len();
        while remaining > 0 {
            let count = remaining.min(buffer.len() as u64) as usize;
            file.read_exact(&mut buffer[..count])
                .map_err(|_| corrupt("SPK digest read", 0))?;
            hash.update(&buffer[..count]);
            remaining -= count as u64;
        }
        validate_identity(&file, &identity)?;
        if crate::download::hex_of(&hash.finalize()) != pin {
            return Err(corrupt("SPK digest does not match the pinned image", 0));
        }
        let (endian, index) = read_index(&mut file, &identity)?;
        validate_identity(&file, &identity)?;
        for (target, parent) in super::AUTOMATIC_PLANET_CENTERS {
            if index
                .segments_for(target)
                .any(|segment| segment.descriptor.center != parent)
            {
                return Err(corrupt(
                    "planet-center segment has the wrong system barycenter",
                    0,
                ));
            }
        }
        Ok(Self {
            file,
            cursor: Mutex::new(()),
            identity,
            endian,
            index,
            cache: Default::default(),
        })
    }

    fn relative_state(
        &self,
        target: i32,
        jd: f64,
        fraction: f64,
    ) -> Result<(ChebState, Descriptor)> {
        validate_identity(&self.file, &self.identity)?;
        let result = relative_state(&self.index, target, jd, fraction, |offset, count| {
            self.cache.acquire(offset, count, || {
                let _cursor = self
                    .cursor
                    .lock()
                    .unwrap_or_else(|error| error.into_inner());
                read_words(&mut &self.file, self.endian, offset, count)
            })
        });
        validate_identity(&self.file, &self.identity)?;
        result
    }
}

fn validate_identity(file: &File, expected: &Metadata) -> Result<()> {
    let actual = metadata(file)?;
    if actual.len() != expected.len() || actual.modified().ok() != expected.modified().ok() {
        return Err(corrupt("immutable SPK file changed", 0));
    }
    Ok(())
}

fn corrupt(what: &str, offset: usize) -> Error {
    Error::SpkCorrupt {
        what: what.to_owned(),
        offset,
    }
}

#[cfg(all(test, feature = "nbody"))]
mod tests {
    use super::super::{SegSpec, SpkKernel};
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};

    fn image(kind: i32, endian: DafEndian) -> Vec<u8> {
        let mut bytes = super::super::build_image(
            endian,
            &[vec![SegSpec {
                start: -172_800.0,
                end: 172_800.0,
                target: 2000002,
                center: 10,
                frame: 1,
                seg_type: kind,
                name: "SYNTHETIC POSITIONAL RECORDS",
            }]],
        );
        let components = if kind == 2 { 3 } else { 6 };
        let record_words = 2 + 3 * components;
        let put_i32 = super::super::reader::put_i32;
        let put_f64 = super::super::reader::put_f64;
        put_i32(&mut bytes, 1080, 385, endian);
        put_i32(&mut bytes, 1084, 384 + 2 * record_words + 4, endian);
        let mut values = Vec::new();
        for midpoint in [-86_400.0, 86_400.0] {
            values.extend([midpoint, 86_400.0]);
            for axis in 0..components {
                values.extend([axis as f64 + 1.0, 0.25, 0.125]);
            }
        }
        values.extend([-172_800.0, 172_800.0, record_words as f64, 2.0]);
        for (i, value) in values.into_iter().enumerate() {
            put_f64(&mut bytes, 3072 + i * 8, value, endian);
        }
        bytes
    }

    fn save(bytes: &[u8]) -> std::path::PathBuf {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let path = std::env::temp_dir().join(format!(
            "libeph-positional-{}-{}.bsp",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::write(&path, bytes).unwrap();
        path
    }

    #[test]
    fn verified_reads_reject_mutation_after_cache_hits() {
        use sha2::{Digest, Sha256};
        let bytes = image(2, DafEndian::Little);
        let path = save(&bytes);
        assert!(VerifiedImage::open(&path, "wrong pin").is_err());
        let pin = crate::download::hex_of(&Sha256::digest(&bytes));
        let image = VerifiedImage::open(&path, &pin).unwrap();
        let expected = image.relative_state(2000002, 2451545.0, 0.5).unwrap();
        assert_eq!(
            image.relative_state(2000002, 2451545.0, 0.5).unwrap(),
            expected
        );
        std::fs::OpenOptions::new()
            .write(true)
            .open(&path)
            .unwrap()
            .set_len(3080)
            .unwrap();
        assert!(matches!(
            image.relative_state(2000002, 2451545.0, 0.5),
            Err(Error::SpkCorrupt { what, .. }) if what == "immutable SPK file changed"
        ));
        std::fs::remove_file(path).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn verified_reads_retain_the_acquired_file_after_replacement() {
        use sha2::{Digest, Sha256};
        let bytes = image(2, DafEndian::Little);
        let path = save(&bytes);
        let pin = crate::download::hex_of(&Sha256::digest(&bytes));
        let image = VerifiedImage::open(&path, &pin).unwrap();
        let expected = image.relative_state(2000002, 2451545.0, -0.5).unwrap();
        let replacement = save(&vec![0; bytes.len()]);
        std::fs::rename(replacement, &path).unwrap();
        assert!(VerifiedImage::open(&path, &pin).is_err());
        assert_eq!(
            image.relative_state(2000002, 2451545.0, -0.5).unwrap(),
            expected
        );
        assert!(image.relative_state(2000002, 2451545.0, 0.5).is_ok());
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn positional_records_match_owned_type2_type3_words_on_both_endiannesses() {
        for endian in [DafEndian::Little, DafEndian::Big] {
            for kind in [2, 3] {
                let bytes = image(kind, endian);
                let path = save(&bytes);
                let mut reader = PositionalChebyshev::open(&path).unwrap();
                let kernel = SpkKernel::decode("owned", bytes).unwrap();
                for fraction in [
                    -2.0_f64,
                    (-2.0_f64).next_up(),
                    -1.0,
                    -1e-12,
                    0.0,
                    1e-12,
                    1.0,
                    2.0,
                ] {
                    let segment =
                        select_two_part_jd(kernel.index(), 2000002, 2451545.0, fraction).unwrap();
                    let decoded =
                        super::super::decode_chebyshev(kernel.bytes(), segment, endian).unwrap();
                    let expected =
                        super::super::eval_chebyshev_two_part_jd(&decoded, 2451545.0, fraction)
                            .unwrap();
                    let actual = reader.relative_state(2000002, 2451545.0, fraction).unwrap();
                    assert_eq!(actual, (expected, segment.descriptor));
                }
                assert!(!reader.contains(1));
                assert!(
                    !reader
                        .trajectory_covers(2000002, 2451545.0, 2.0_f64.next_up())
                        .unwrap()
                );
                assert!(matches!(
                    reader.relative_state(2000002, 2451545.0, 3.0),
                    Err(Error::SpkOutOfRange { .. })
                ));
                std::fs::remove_file(path).unwrap();
            }
        }
    }

    #[test]
    fn stored_velocity_is_not_substituted_for_a_type3_trajectory_tangent() {
        let path = save(&image(3, DafEndian::Little));
        let mut reader = PositionalChebyshev::open(&path).unwrap();
        assert!(reader.relative_state(2000002, 2451545.0, 0.5).is_ok());
        assert!(
            !reader.trajectory_covers(2000002, 2451545.0, 0.5).unwrap(),
            "Type3 must retain the existing integration route, not introduce a new public error"
        );
        assert!(matches!(
            reader.trajectory_state(2000002, 2451545.0, 0.5),
            Err(Error::SpkUnsupportedType { seg_type: 3 })
        ));
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn positional_reader_keeps_a_sparse_large_file_without_copying_its_image() {
        let path = save(&image(2, DafEndian::Little));
        let file = std::fs::OpenOptions::new().write(true).open(&path).unwrap();
        file.set_len(616 * 1024 * 1024).unwrap();
        let mut reader = PositionalChebyshev::open(&path).unwrap();
        assert_eq!(reader.index.len(), 1);
        assert!(reader.relative_state(2000002, 2451545.0, 0.5).is_ok());
        // There is no whole-image field; only the file and small index remain.
        assert!(std::mem::size_of_val(&reader) < 1024);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn metadata_and_record_failures_are_distinct_from_descriptor_misses() {
        let mut bytes = image(2, DafEndian::Little);
        let path = save(&bytes);
        let mut reader = PositionalChebyshev::open(&path).unwrap();
        assert!(reader.trajectory_covers(2000002, 2451545.0, 0.5).unwrap());
        std::fs::OpenOptions::new()
            .write(true)
            .open(&path)
            .unwrap()
            .set_len(3080)
            .unwrap();
        assert!(matches!(
            reader.relative_state(2000002, 2451545.0, 0.5),
            Err(Error::SpkCorrupt { .. })
        ));
        std::fs::remove_file(path).unwrap();
        for (offset, value) in [(1072, 17), (1076, 21)] {
            super::super::reader::put_i32(&mut bytes, offset, value, DafEndian::Little);
            let path = save(&bytes);
            let mut reader = PositionalChebyshev::open(&path).unwrap();
            assert!(reader.relative_state(2000002, 2451545.0, 0.0).is_err());
            std::fs::remove_file(path).unwrap();
            bytes = image(2, DafEndian::Little);
        }
    }
}
