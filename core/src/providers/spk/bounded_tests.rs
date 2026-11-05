// SPDX-License-Identifier: AGPL-3.0
//! Prepared bounded-storage checks; no astronomical coefficients are fixtures.
use super::*;
use std::sync::atomic::{AtomicUsize, Ordering};

struct Temporary(std::path::PathBuf);
impl Temporary {
    fn new(bytes: &[u8]) -> Self {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let path = std::env::temp_dir().join(format!(
            "libephemeris-force-spk-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::write(&path, bytes).unwrap();
        Self(path)
    }
}
impl Drop for Temporary {
    fn drop(&mut self) {
        std::fs::remove_file(&self.0).unwrap();
    }
}

pub(crate) fn image(
    endian: DafEndian,
    seg_type: i32,
    targets: usize,
    init: f64,
    interval: f64,
) -> Vec<u8> {
    let specs: Vec<_> = (0..targets)
        .map(|i| SegSpec {
            start: init,
            end: init + 2.0 * interval,
            target: 2_000_001 + i as i32,
            center: 10,
            frame: FRAME_J2000,
            seg_type,
            name: "SYNTHETIC FORCE RECORD",
        })
        .collect();
    let groups: Vec<_> = specs.chunks(25).map(|group| group.to_vec()).collect();
    let mut bytes = build_image(endian, &groups);
    for i in 0..targets {
        let components = if seg_type == TYPE_CHEBYSHEV_POS_VEL {
            6
        } else {
            3
        };
        let rsize = 2 + components * 3;
        let length = 2 * rsize + 4;
        let offset = bytes.len();
        bytes.resize(offset + length * 8, 0);
        let summary = (1 + 2 * (i / 25)) * RECORD_BYTES + (3 + (i % 25) * SPK_SUMMARY_WORDS) * 8;
        put_i32(&mut bytes, summary + 32, (offset / 8 + 1) as i32, endian);
        put_i32(
            &mut bytes,
            summary + 36,
            (offset / 8 + length) as i32,
            endian,
        );
        for record in 0..2 {
            let base = offset + record * rsize * 8;
            put_f64(
                &mut bytes,
                base,
                init + (record as f64 + 0.5) * interval,
                endian,
            );
            put_f64(&mut bytes, base + 8, interval / 2.0, endian);
            for component in 0..components {
                let set = base + (2 + component * 3) * 8;
                put_f64(
                    &mut bytes,
                    set,
                    (i * 20 + record * 10 + component + 1) as f64,
                    endian,
                );
                put_f64(&mut bytes, set + 8, 0.5, endian);
                put_f64(&mut bytes, set + 16, 0.25, endian);
            }
        }
        let tail = offset + 2 * rsize * 8;
        for (j, value) in [init, interval, rsize as f64, 2.0].into_iter().enumerate() {
            put_f64(&mut bytes, tail + j * 8, value, endian);
        }
    }
    bytes.resize(bytes.len().div_ceil(RECORD_BYTES) * RECORD_BYTES, 0);
    bytes
}

#[test]
fn bounded_and_owned_storage_share_bit_identical_coefficient_mathematics() {
    for endian in [DafEndian::Little, DafEndian::Big] {
        for seg_type in [TYPE_CHEBYSHEV_POS, TYPE_CHEBYSHEV_POS_VEL] {
            let bytes = image(endian, seg_type, 2, 0.0, 2.0);
            let temporary = Temporary::new(&bytes);
            let bounded = crate::spk::ChebyshevFile::open(&temporary.0).unwrap();
            let index = SpkIndex::parse(&bytes).unwrap();
            assert_eq!(bounded.targets(), [2_000_001, 2_000_002]);
            assert_eq!(bounded.coverage(2_000_001), [(0.0, 4.0)]);
            for target in bounded.targets() {
                for seconds in [0.0_f64, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0] {
                    let fraction = seconds / SECONDS_PER_DAY;
                    let segment =
                        select_two_part_jd(&index, target, J2000_JD_TDB, fraction).unwrap();
                    let decoded = decode_chebyshev(&bytes, segment, endian).unwrap();
                    let owned =
                        eval_chebyshev_two_part_jd(&decoded, J2000_JD_TDB, fraction).unwrap();
                    let (center, actual) = bounded.state(target, J2000_JD_TDB, fraction).unwrap();
                    assert_eq!(center, 10);
                    assert_eq!(actual, owned);
                }
            }
        }
    }
}

#[test]
fn tiny_split_epoch_controls_record_ownership_and_closed_global_coverage() {
    let bytes = image(DafEndian::Little, TYPE_CHEBYSHEV_POS, 1, -172800.0, 86400.0);
    let temporary = Temporary::new(&bytes);
    let bounded = ChebyshevFile::open(&temporary.0).unwrap();
    let tiny = 2.0_f64.powi(-80);
    assert_eq!(
        bounded
            .state(2_000_001, J2000_JD_TDB - 1.0, -tiny)
            .unwrap()
            .1
            .position_km[0],
        1.75
    );
    for fraction in [0.0, tiny] {
        assert_eq!(
            bounded
                .state(2_000_001, J2000_JD_TDB - 1.0, fraction)
                .unwrap()
                .1
                .position_km[0],
            10.75
        );
    }
    assert!(matches!(
        bounded.state(2_000_001, J2000_JD_TDB, tiny),
        Err(Error::SpkOutOfRange { .. })
    ));
    assert!(matches!(
        bounded.state(2_000_001, J2000_JD_TDB - 2.0, -tiny),
        Err(Error::SpkOutOfRange { .. })
    ));
}

#[test]
fn coefficient_cache_has_bounded_ownership_and_does_not_retain_evicted_records() {
    let bytes = image(DafEndian::Little, TYPE_CHEBYSHEV_POS, 129, 0.0, 2.0);
    let temporary = Temporary::new(&bytes);
    let bounded = ChebyshevFile::open(&temporary.0).unwrap();
    let first = bounded.record(0, 0).unwrap();
    let weak = Arc::downgrade(&first);
    assert!(Arc::ptr_eq(&first, &bounded.record(0, 0).unwrap()));
    for segment in 1..129 {
        bounded.record(segment, 0).unwrap();
    }
    {
        let owner = bounded.file.lock().unwrap();
        assert_eq!(owner.records.len(), CACHE_ENTRIES);
        assert!(owner.bytes <= CACHE_BYTES);
        assert!(owner.records.iter().all(|(key, _)| *key != (0, 0)));
    }
    assert!(weak.upgrade().is_some());
    drop(first);
    assert!(weak.upgrade().is_none());
}

#[test]
fn bounded_metadata_admission_and_unserved_inputs_are_explicit_errors() {
    let bytes = image(DafEndian::Little, TYPE_CHEBYSHEV_POS, 2, 0.0, 2.0);
    let mut cursor = std::io::Cursor::new(&bytes);
    assert!(matches!(
        SpkIndex::read_directory_bounded_from(&mut cursor, bytes.len(), 1),
        Err(Error::SpkCorrupt { .. })
    ));
    let temporary = Temporary::new(&bytes);
    let bounded = ChebyshevFile::open(&temporary.0).unwrap();
    assert!(bounded.coverage(999).is_empty());
    assert!(matches!(
        bounded.state(999, J2000_JD_TDB, 0.0),
        Err(Error::SpkBodyNotFound { target: 999 })
    ));
    assert!(matches!(
        bounded.state(2_000_001, f64::NAN, 0.0),
        Err(Error::NonFiniteArgument { .. })
    ));
    let unsupported = image(DafEndian::Little, TYPE_EMDA, 1, 0.0, 2.0);
    let other = Temporary::new(&unsupported);
    assert!(matches!(
        ChebyshevFile::open(&other.0),
        Err(Error::SpkUnsupportedType {
            seg_type: TYPE_EMDA
        })
    ));
}

#[test]
fn failed_record_io_leaves_the_cache_intact_and_can_be_retried() {
    let bytes = image(DafEndian::Little, TYPE_CHEBYSHEV_POS, 2, 0.0, 2.0);
    let temporary = Temporary::new(&bytes);
    let bounded = ChebyshevFile::open(&temporary.0).unwrap();
    let retained = bounded.record(0, 0).unwrap();
    let expected = ChebyshevFile::open(&temporary.0)
        .unwrap()
        .state(2_000_002, J2000_JD_TDB, 0.5 / SECONDS_PER_DAY)
        .unwrap();

    // Inject a real read failure without modifying any bytes or truncating an
    // open image. This synthetic reader uses ordinary File I/O, never mmap.
    // Only the test can replace its private handle with a write-only handle.
    let write_only = std::fs::OpenOptions::new()
        .write(true)
        .open(&temporary.0)
        .unwrap();
    let readable = {
        let mut owner = bounded.file.lock().unwrap();
        std::mem::replace(&mut owner.file, write_only)
    };
    assert!(matches!(
        bounded.state(2_000_002, J2000_JD_TDB, 0.5 / SECONDS_PER_DAY),
        Err(Error::SpkCorrupt { .. })
    ));
    {
        let mut owner = bounded.file.lock().unwrap();
        assert_eq!(owner.records.len(), 1);
        assert_eq!(owner.records[0].0, (0, 0));
        assert!(Arc::ptr_eq(&owner.records[0].1, &retained));
        assert_eq!(owner.bytes, retained.len() * 8);
        owner.file = readable;
    }
    assert_eq!(
        bounded
            .state(2_000_002, J2000_JD_TDB, 0.5 / SECONDS_PER_DAY)
            .unwrap(),
        expected
    );
    assert!(Arc::ptr_eq(&bounded.record(0, 0).unwrap(), &retained));
}

#[test]
fn independent_worker_evaluations_share_only_immutable_coefficient_records() {
    let bytes = image(DafEndian::Big, TYPE_CHEBYSHEV_POS, 2, 0.0, 2.0);
    let temporary = Temporary::new(&bytes);
    let bounded = ChebyshevFile::open(&temporary.0).unwrap();
    let reference = bounded
        .state(2_000_002, J2000_JD_TDB, 1.5 / SECONDS_PER_DAY)
        .unwrap();
    std::thread::scope(|scope| {
        for _ in 0..8 {
            let bounded = &bounded;
            scope.spawn(move || {
                for target in [2_000_001, 2_000_002] {
                    bounded
                        .state(target, J2000_JD_TDB, 0.5 / SECONDS_PER_DAY)
                        .unwrap();
                }
                assert_eq!(
                    bounded
                        .state(2_000_002, J2000_JD_TDB, 1.5 / SECONDS_PER_DAY)
                        .unwrap(),
                    reference
                );
            });
        }
    });
}
