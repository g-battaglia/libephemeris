// SPDX-License-Identifier: AGPL-3.0
//! Synthetic format tests contain no astronomical measurements or vendor data.

use super::*;
use std::sync::atomic::{AtomicUsize, Ordering};

pub(in super::super) struct Fixture {
    bytes: Vec<u8>,
    words: usize,
    order: Endian,
}

impl Fixture {
    /// Shared minimal format image for other native data-boundary refusal tests.
    pub(in super::super) fn little() -> Self {
        Self::new(Endian::Little, 440)
    }

    /// Complete, deliberately non-astronomical metadata for force-owner tests.
    #[cfg(feature = "nbody")]
    pub(in super::super) fn force_ready() -> Self {
        let mut fixture = Self::little();
        let constants = [
            ("DENUM", 440.0),
            ("AU", 2.0),
            ("EMRAT", 3.0),
            ("GMS", 1.0),
            ("J2SUN", 0.001),
            ("GMB", 8.0),
            ("GM1", 0.5),
            ("GM2", 0.5),
            ("GM4", 0.5),
            ("GM5", 0.5),
            ("GM6", 0.5),
            ("GM7", 0.5),
            ("GM8", 0.5),
            ("GM9", 0.5),
            ("CLIGHT", 20.0),
            ("ASUN", 4.0),
            ("RE", 1.0),
            ("J2E", 0.01),
            ("J3E", -0.002),
            ("J4E", -0.003),
            ("MA0001", 0.001),
            ("MA0002", 0.002),
            ("MA0004", 0.004),
        ];
        fixture.integer(2676, constants.len());
        for (index, (name, value)) in constants.into_iter().enumerate() {
            let offset = TITLE_BYTES + index * NAME_BYTES;
            fixture.bytes[offset..offset + NAME_BYTES].fill(b' ');
            fixture.bytes[offset..offset + name.len()].copy_from_slice(name.as_bytes());
            fixture.scalar(fixture.words * 8 + index * 8, value);
        }
        fixture
    }

    /// Change only the synthetic Sun's declared GM (AU³/day²). All stored
    /// trajectories and other constants remain identical, so a force-state
    /// difference cannot be attributed to a different SPK or initial state.
    #[cfg(feature = "nbody")]
    pub(in super::super) fn with_solar_gm(mut self, gm: f64) -> Self {
        self.scalar(self.words * 8 + 3 * 8, gm);
        self
    }

    /// Append one explicitly named synthetic header constant. The DE name
    /// field is six bytes; no astronomical value or external file is copied.
    #[cfg(feature = "nbody")]
    pub(in super::super) fn with_constant(mut self, name: &str, value: f64) -> Self {
        assert!(!name.is_empty() && name.len() <= NAME_BYTES && name.is_ascii());
        let index = self.order.integer(&self.bytes, 2676).unwrap();
        assert!(index < 400 && index < self.words);
        let offset = TITLE_BYTES + index * NAME_BYTES;
        self.bytes[offset..offset + NAME_BYTES].fill(b' ');
        self.bytes[offset..offset + name.len()].copy_from_slice(name.as_bytes());
        self.scalar(self.words * 8 + index * 8, value);
        self.integer(2676, index + 1);
        self
    }

    fn new(order: Endian, version: usize) -> Self {
        // Two records, every pointer family populated. Degree-two coefficients
        // are embedded in longer rows so the header fits its own record.
        let mut starts = [0; 15];
        let mut words = 2;
        for (slot, start) in starts.iter_mut().enumerate() {
            *start = words;
            let components = match slot {
                11 => 2,
                14 => 1,
                _ => 3,
            };
            words += components * 10 * 2;
        }
        let mut result = Self {
            bytes: vec![0; 4 * words * 8],
            words,
            order,
        };
        result.scalar(2652, 1000.0);
        result.scalar(2660, 1064.0);
        result.scalar(2668, 32.0);
        result.integer(2676, 5);
        result.scalar(2680, 2.0);
        result.scalar(2688, 3.0);
        result.integer(2840, version);
        for (index, name) in ["DENUM ", "AU    ", "EMRAT ", "GMS   ", "J2SUN "]
            .into_iter()
            .enumerate()
        {
            let start = TITLE_BYTES + index * NAME_BYTES;
            result.bytes[start..start + NAME_BYTES].copy_from_slice(name.as_bytes());
        }
        for (index, value) in [version as f64, 2.0, 3.0, 1.0, 0.0].into_iter().enumerate() {
            result.scalar(words * 8 + index * 8, value);
        }
        for (slot, start) in starts.into_iter().enumerate() {
            let offset = match slot {
                0..=11 => 2696 + slot * 12,
                12 => 2844,
                _ => FIXED_HEADER_BYTES + (slot - 13) * 12,
            };
            result.integer(offset, start + 1);
            result.integer(offset + 4, 10);
            result.integer(offset + 8, 2);
            let components = match slot {
                11 => 2,
                14 => 1,
                _ => 3,
            };
            for record in 0..2 {
                let base = (record + 2) * words;
                result.scalar(base * 8, 1000.0 + 32.0 * record as f64);
                result.scalar((base + 1) * 8, 1032.0 + 32.0 * record as f64);
                for subinterval in 0..2 {
                    for component in 0..components {
                        let offset =
                            (base + start + (subinterval * components + component) * 10) * 8;
                        result.scalar(offset, (10 * slot + component) as f64);
                        result.scalar(offset + 8, 2.0);
                        result.scalar(offset + 16, 3.0);
                    }
                }
            }
        }
        result
    }

    fn scalar(&mut self, offset: usize, value: f64) {
        let bytes = match self.order {
            Endian::Little => value.to_le_bytes(),
            Endian::Big => value.to_be_bytes(),
        };
        self.bytes[offset..offset + 8].copy_from_slice(&bytes);
    }

    fn integer(&mut self, offset: usize, value: usize) {
        let bytes = match self.order {
            Endian::Little => (value as i32).to_le_bytes(),
            Endian::Big => (value as i32).to_be_bytes(),
        };
        self.bytes[offset..offset + 4].copy_from_slice(&bytes);
    }

    fn header(&self) -> Result<Header, DeError> {
        Header::parse(&self.bytes, self.bytes.len() as u64)
    }

    pub(in super::super) fn open(&self) -> (ClassicDe, Temporary) {
        let temporary = Temporary::new(&self.bytes);
        (ClassicDe::open(&temporary.0).unwrap(), temporary)
    }
}

pub(in super::super) struct Temporary(std::path::PathBuf);

impl Temporary {
    #[cfg(feature = "nbody")]
    pub(in super::super) fn path(&self) -> &std::path::Path {
        &self.0
    }

    pub(in super::super) fn new(bytes: &[u8]) -> Self {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let path = std::env::temp_dir().join(format!(
            "libephemeris-native-de-{}-{}",
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

#[test]
fn both_versions_and_endian_orders_have_identical_states_and_constants() {
    for version in [440, 441] {
        for order in [Endian::Little, Endian::Big] {
            let fixture = Fixture::new(order, version);
            let (reader, _temporary) = fixture.open();
            assert_eq!(reader.header.version, version);
            assert_eq!(reader.header.record_words, fixture.words);
            assert_eq!(reader.constant("GMS").unwrap(), 1.0);
            assert!(matches!(
                reader.constant("absent"),
                Err(DeError::MissingConstant(_))
            ));
            // t=1008 is x=0 in the first 16-day subinterval. Its polynomial is
            // c0 + 2*T1(x) + 3*T2(x), with AU=2 and dx/dt=1/8 day^-1.
            let sun = reader.state(10, [1008.0, 0.0]).unwrap();
            assert_eq!(sun.state.position, [48.5, 49.0, 49.5]);
            assert_eq!(sun.state.velocity, [0.125; 3]);
            assert_eq!(sun.acceleration, [0.09375; 3]);
            assert_eq!(
                reader.state(0, [1008.0, 0.0]).unwrap().state.position,
                [0.0; 3]
            );
            assert!(matches!(
                reader.state(999, [1008.0, 0.0]),
                Err(DeError::UnsupportedTarget(999))
            ));
        }
    }
}

#[test]
fn batch_frame_preserves_order_duplicates_and_the_scalar_mathematical_route() {
    let fixture = Fixture::new(Endian::Little, 441);
    let (reader, _temporary) = fixture.open();
    let targets = [10, 399, 301, 0, 10, 5, 3];
    let frame = reader
        .states(targets, [1032.0, -2.0_f64.powi(-80)])
        .unwrap();
    for (body, batched) in targets.into_iter().zip(frame) {
        let scalar = reader.state(body, [1032.0, -2.0_f64.powi(-80)]).unwrap();
        assert_eq!(scalar.state, batched.state);
        assert_eq!(scalar.acceleration, batched.acceleration);
    }
    assert_eq!(frame[0].state, frame[4].state);
    assert_eq!(reader.file.lock().unwrap().current.as_ref().unwrap().0, 0);
}

#[test]
fn only_one_record_is_retained_and_old_images_follow_borrowing_lifetimes() {
    let fixture = Fixture::new(Endian::Big, 440);
    let (reader, _temporary) = fixture.open();
    assert!(reader.file.lock().unwrap().current.is_none());
    let first = reader.record(0).unwrap();
    let old = Arc::downgrade(&first);
    assert!(Arc::ptr_eq(&first, &reader.record(0).unwrap()));
    let second = reader.record(1).unwrap();
    assert!(!Arc::ptr_eq(&first, &second));
    assert_eq!(reader.file.lock().unwrap().current.as_ref().unwrap().0, 1);
    assert_eq!(Arc::strong_count(&first), 1);
    drop(first);
    assert!(old.upgrade().is_none());
    assert!(Arc::ptr_eq(&second, &reader.record(1).unwrap()));
}

#[test]
fn zero_or_empty_frames_validate_coverage_without_acquiring_a_record() {
    let fixture = Fixture::new(Endian::Little, 441);
    let (reader, _temporary) = fixture.open();
    assert_eq!(reader.states([], [1000.0, 0.0]).unwrap().len(), 0);
    let frame = reader.states([0, 0], [1064.0, 0.0]).unwrap();
    assert!(
        frame
            .into_iter()
            .all(|value| value.state.position == [0.0; 3]
                && value.state.velocity == [0.0; 3]
                && value.acceleration == [0.0; 3])
    );
    assert!(reader.file.lock().unwrap().current.is_none());
    assert!(matches!(
        reader.states([], [1064.0, 2.0_f64.powi(-80)]),
        Err(DeError::OutsideCoverage)
    ));
    assert!(matches!(
        reader.states([0], [1000.0, -2.0_f64.powi(-80)]),
        Err(DeError::OutsideCoverage)
    ));
}

#[test]
fn selected_numeric_coverage_is_not_replaced_by_the_text_title() {
    let mut fixture = Fixture::new(Endian::Little, 441);
    fixture.bytes[..TITLE_BYTES].fill(b'9');
    let header = fixture.header().unwrap();
    assert_eq!(
        (header.start, header.end, header.step),
        (1000.0, 1064.0, 32.0)
    );
}

#[test]
fn split_epoch_retains_sub_ulp_offsets_at_record_boundaries() {
    let fixture = Fixture::new(Endian::Little, 440);
    let header = fixture.header().unwrap();
    let small = 2.0_f64.powi(-45);
    assert_eq!(1032.0 + small, 1032.0);
    assert_eq!(
        header.record_index([1032.0, small]).unwrap(),
        (1, [small, 0.0])
    );
    let (previous, delta) = header.record_index([1032.0, -small]).unwrap();
    assert_eq!(previous, 0);
    assert!(delta[0] + delta[1] < 32.0);
    assert_eq!(
        header.record_index([1064.0, 0.0]).unwrap(),
        (1, [32.0, 0.0])
    );
    assert!(matches!(
        header.record_index([1064.0, small]),
        Err(DeError::OutsideCoverage)
    ));
    assert!(matches!(
        header.record_index([1000.0, -small]),
        Err(DeError::OutsideCoverage)
    ));
    assert!(matches!(
        header.record_index([f64::NAN, 0.0]),
        Err(DeError::Malformed(_))
    ));
}

#[test]
fn offsets_smaller_than_a_record_ulp_still_control_ownership_and_domain() {
    let fixture = Fixture::new(Endian::Little, 440);
    let header = fixture.header().unwrap();
    let tiny = 2.0_f64.powi(-80);
    assert_eq!(32.0 - tiny, 32.0);
    assert_eq!(header.record_index([1032.0, -tiny]).unwrap().0, 0);
    assert!(matches!(
        header.record_index([1064.0, tiny]),
        Err(DeError::OutsideCoverage)
    ));
    assert!(matches!(
        header.record_index([1000.0, -tiny]),
        Err(DeError::OutsideCoverage)
    ));
}

#[test]
fn tiny_offsets_also_select_the_correct_polynomial_subinterval() {
    let mut fixture = Fixture::new(Endian::Little, 440);
    let header = fixture.header().unwrap();
    let sun = header.blocks[10];
    let next_constant = (2 * fixture.words + sun.start + 3 * sun.coefficients) * 8;
    fixture.scalar(next_constant, 1000.0);
    let (reader, _temporary) = fixture.open();
    let tiny = 2.0_f64.powi(-80);
    // A deliberately discontinuous synthetic row makes ownership observable;
    // real ephemeris fits are smooth and must not be used to hide this defect.
    let before = reader.state(10, [1016.0, -tiny]).unwrap();
    let at = reader.state(10, [1016.0, 0.0]).unwrap();
    assert_eq!(before.state.position[0], (100.0 + 2.0 + 3.0) / 2.0);
    assert_eq!(at.state.position[0], (1000.0 - 2.0 + 3.0) / 2.0);
}

#[test]
fn earth_and_moon_are_reconstructed_from_emb_and_geocentric_moon() {
    let fixture = Fixture::new(Endian::Little, 440);
    let (reader, _temporary) = fixture.open();
    let emb = reader.state(3, [1008.0, 0.0]).unwrap();
    let earth = reader.state(399, [1008.0, 0.0]).unwrap();
    let moon = reader.state(301, [1008.0, 0.0]).unwrap();
    for axis in 0..3 {
        assert_eq!(
            (3.0 * earth.state.position[axis] + moon.state.position[axis]) / 4.0,
            emb.state.position[axis]
        );
        assert_eq!(
            (3.0 * earth.state.velocity[axis] + moon.state.velocity[axis]) / 4.0,
            emb.state.velocity[axis]
        );
        assert_eq!(
            (3.0 * earth.acceleration[axis] + moon.acceleration[axis]) / 4.0,
            emb.acceleration[axis]
        );
    }
}

#[test]
fn corrupt_shapes_overlaps_versions_and_lengths_are_refused() {
    for (offset, value) in [
        (2840, 999),
        (2676, LIMIT_CONSTANTS + 1),
        (2696, 1),
        (2700, 0),
        (2708, 3),
    ] {
        let mut fixture = Fixture::new(Endian::Little, 440);
        fixture.integer(offset, value);
        assert!(fixture.header().is_err(), "offset {offset}, value {value}");
    }
    let mut fixture = Fixture::new(Endian::Little, 440);
    fixture.bytes.pop();
    assert!(matches!(fixture.header(), Err(DeError::Malformed(_))));
    let mut fixture = Fixture::new(Endian::Little, 440);
    fixture.scalar(2660, 1063.0);
    assert!(matches!(fixture.header(), Err(DeError::Malformed(_))));
    fixture.scalar(2660, f64::NAN);
    assert!(matches!(fixture.header(), Err(DeError::Malformed(_))));
}

#[test]
fn inconsistent_constant_identity_and_changed_record_bounds_are_refused() {
    let mut fixture = Fixture::new(Endian::Little, 440);
    fixture.scalar(fixture.words * 8 + 8, 3.0);
    let temporary = Temporary::new(&fixture.bytes);
    assert!(matches!(
        ClassicDe::open(&temporary.0),
        Err(DeError::Malformed(_))
    ));
    let mut fixture = Fixture::new(Endian::Little, 440);
    fixture.scalar(2 * fixture.words * 8, 999.0);
    let (reader, _temporary) = fixture.open();
    assert!(matches!(
        reader.state(10, [1008.0, 0.0]),
        Err(DeError::Malformed(_))
    ));
}

#[test]
fn retained_file_owner_is_thread_safe_without_duplicate_images() {
    let fixture = Fixture::new(Endian::Big, 441);
    let (reader, _temporary) = fixture.open();
    let reader = std::sync::Arc::new(reader);
    let expected = reader.state(10, [1008.0, 0.0]).unwrap();
    let threads: Vec<_> = (0..8)
        .map(|_| {
            let shared = std::sync::Arc::clone(&reader);
            std::thread::spawn(move || shared.state(10, [1008.0, 0.0]).unwrap())
        })
        .collect();
    for thread in threads {
        let value = thread.join().unwrap();
        assert_eq!(value.state, expected.state);
        assert_eq!(value.acceleration, expected.acceleration);
    }
}
