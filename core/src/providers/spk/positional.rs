//! Bounded record reads for the retained ASSIST asteroid image.
//!
//! The caller keeps the acquired file immutable. One descriptor index and an
//! open file are retained; coefficient payloads are read only for the selected
//! record. Owned and positional paths share record selection and evaluation.

use super::decode::chebyshev_directory;
use super::evaluation::{
    chebyshev_record_coordinates, eval_selected_chebyshev_record, select_two_part_jd,
};
use super::{ChebState, DafEndian, Descriptor, Error, FileRecord, Result, SpkIndex};
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

pub(crate) struct PositionalChebyshev {
    file: File,
    endian: DafEndian,
    index: SpkIndex,
}

impl PositionalChebyshev {
    /// Acquire one immutable regular SPK without copying its payload.
    pub(crate) fn open(path: &Path) -> Result<Self> {
        let mut file = File::open(path).map_err(|error| {
            if error.kind() == std::io::ErrorKind::NotFound {
                Error::SpkFileNotFound {
                    path: path.display().to_string(),
                }
            } else {
                corrupt("SPK file read", 0)
            }
        })?;
        let length = file.metadata().map_err(|_| corrupt("SPK file length", 0))?;
        if !length.is_file() {
            return Err(corrupt("SPK source must be a regular file", 0));
        }
        let length = usize::try_from(length.len()).map_err(|_| corrupt("SPK file length", 0))?;
        let mut header = [0; super::RECORD_BYTES];
        file.read_exact(&mut header)
            .map_err(|_| corrupt("SPK file record", 0))?;
        let endian = FileRecord::parse(&header)?.endian;
        let index = SpkIndex::read_directory_from(&mut file, length)?;
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
        let descriptor = select_two_part_jd(&self.index, target, jd, fraction)?.descriptor;
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
        let trailer = self.words(start + (length - 4) * 8, 4)?;
        let directory = chebyshev_directory(
            descriptor.seg_type,
            length,
            start,
            trailer.try_into().expect("four directory words"),
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
        let record = self.words(record_offset, directory.record_words)?;
        let state = eval_selected_chebyshev_record(
            &directory,
            descriptor.seg_type,
            &record,
            index,
            offset,
            record_offset,
            Some((hi, lo)),
        )?
        .0;
        Ok((state, descriptor))
    }

    fn words(&mut self, offset: usize, count: usize) -> Result<Vec<f64>> {
        let length = count
            .checked_mul(8)
            .ok_or_else(|| corrupt("SPK record length", offset))?;
        let mut bytes = vec![0; length];
        self.file
            .seek(SeekFrom::Start(offset as u64))
            .map_err(|_| corrupt("SPK record seek", offset))?;
        self.file
            .read_exact(&mut bytes)
            .map_err(|_| corrupt("SPK record read", offset))?;
        Ok(bytes
            .chunks_exact(8)
            .map(|word| self.endian.read_f64(word))
            .collect())
    }
}

fn corrupt(what: &str, offset: usize) -> Error {
    Error::SpkCorrupt {
        what: what.to_owned(),
        offset,
    }
}

#[cfg(test)]
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
