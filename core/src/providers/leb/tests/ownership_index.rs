//! Ordered file-header bounds and pre-advice segment arithmetic guards.

use super::*;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

/// One finite-grid LEB1 body: degree zero, three components, one segment,
/// body bounds [-MAX, MAX], interval one day, coefficient bytes at offset164.
const WARM_OVERFLOW_IMAGE: [u8; 188] = [
    0x4c, 0x45, 0x42, 0x31, 0x01, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xef, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xef, 0x7f,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x70, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x34, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0xa4, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x18, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xef, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xef, 0x7f, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0xf0, 0x3f, 0x00, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0xa4, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x40, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x08, 0x40,
];

/// Carry arbitrary file-header words while keeping the numerical body grid finite.
fn header_image(start: f64, end: f64) -> Vec<u8> {
    let mut image = leb1_const_image(0, 1.0, 0.0, 1.0);
    image[16..24].copy_from_slice(&start.to_le_bytes());
    image[24..32].copy_from_slice(&end.to_le_bytes());
    image
}

/// Read range words from genuine decoded constituent indexes in the given order.
fn composite_range(images: &[Vec<u8>]) -> (f64, f64) {
    CompositeIndex::new(
        images
            .iter()
            .map(|image| (FileIndex::decode(image).unwrap(), image.as_slice()))
            .collect(),
    )
    .unwrap()
    .jd_range()
}

/// Build genuine tier composites; constructor priority supersedes input order.
fn tiered_range(images: &[(LebTier, Vec<u8>)]) -> (f64, f64) {
    TieredIndex::new(
        images
            .iter()
            .map(|(tier, image)| {
                (
                    *tier,
                    CompositeIndex::new(vec![(
                        FileIndex::decode(image).unwrap(),
                        image.as_slice(),
                    )])
                    .unwrap(),
                )
            })
            .collect(),
    )
    .unwrap()
    .jd_range()
}

/// Compare payloads and zero signs without NaN equality or numeric normalization.
fn assert_range_words(actual: (f64, f64), expected: (f64, f64)) {
    assert_eq!(actual.0.to_bits(), expected.0.to_bits());
    assert_eq!(actual.1.to_bits(), expected.1.to_bits());
}

#[test]
fn leb_index_ownership_composite_first_later_single_and_all_nan_bounds() {
    let a = f64::from_bits(0x7ff8_0000_0000_0011);
    let b = f64::from_bits(0xfff8_0000_0000_0022);
    for (bounds, expected) in [
        (vec![(a, b), (1000.0, 1100.0)], (a, b)),
        (vec![(1000.0, 1100.0), (a, b)], (1000.0, 1100.0)),
        (vec![(a, b), (b, a)], (a, b)),
        (vec![(a, b)], (a, b)),
        (vec![(a, 5.0), (-3.0, b), (1.0, 10.0)], (a, 10.0)),
        (
            vec![(f64::INFINITY, f64::NEG_INFINITY), (-3.0, 7.0)],
            (-3.0, 7.0),
        ),
    ] {
        let images: Vec<_> = bounds
            .into_iter()
            .map(|(s, e)| header_image(s, e))
            .collect();
        assert_range_words(composite_range(&images), expected);
    }
}

#[test]
fn leb_index_ownership_composite_signed_zero_ties_retain_first_words() {
    for first in [(0.0, -0.0), (-0.0, 0.0)] {
        let images = [
            header_image(first.0, first.1),
            header_image(-first.0, -first.1),
        ];
        assert_range_words(composite_range(&images), first);
    }
}

#[test]
fn leb_index_ownership_tiered_first_later_single_and_all_nan_bounds() {
    let a = f64::from_bits(0x7ff8_0000_0000_0011);
    let b = f64::from_bits(0xfff8_0000_0000_0022);
    for (base, medium, expected) in [
        ((a, b), (1000.0, 1100.0), (a, b)),
        ((1000.0, 1100.0), (a, b), (1000.0, 1100.0)),
        ((a, b), (b, a), (a, b)),
    ] {
        // Installation is reversed deliberately; base owns the first range.
        let images = [
            (LebTier::Medium, header_image(medium.0, medium.1)),
            (LebTier::Base, header_image(base.0, base.1)),
        ];
        assert_range_words(tiered_range(&images), expected);
    }
    assert_range_words(
        tiered_range(&[(LebTier::Extended, header_image(a, b))]),
        (a, b),
    );
}

#[test]
fn leb_index_ownership_tiered_signed_zero_uses_priority_not_installation() {
    for base in [(0.0, -0.0), (-0.0, 0.0)] {
        let images = [
            (LebTier::Extended, header_image(-base.0, -base.1)),
            (LebTier::Base, header_image(base.0, base.1)),
        ];
        assert_range_words(tiered_range(&images), base);
    }
}

#[test]
fn leb_index_ownership_nan_composite_range_defers_aux_to_covering_tier() {
    let mut first = aux_image(vec![(
        SECTION_DELTA_T,
        delta_t_payload(&[(0.0, 1.0), (16.0, 1.0)]),
    )]);
    first[16..24].copy_from_slice(&f64::NAN.to_le_bytes());
    first[24..32].copy_from_slice(&f64::NAN.to_le_bytes());
    let second = header_image(0.0, 16.0);
    let medium = aux_image(vec![(
        SECTION_DELTA_T,
        delta_t_payload(&[(0.0, 2.0), (16.0, 2.0)]),
    )]);
    let base = CompositeIndex::new(vec![
        (FileIndex::decode(&first).unwrap(), first.as_slice()),
        (FileIndex::decode(&second).unwrap(), second.as_slice()),
    ])
    .unwrap();
    let medium = CompositeIndex::new(vec![(
        FileIndex::decode(&medium).unwrap(),
        medium.as_slice(),
    )])
    .unwrap();
    let tiered = TieredIndex::new(vec![(LebTier::Medium, medium), (LebTier::Base, base)]).unwrap();
    assert_eq!(tiered.delta_t(1.0).unwrap(), 2.0);
}

#[test]
fn leb_index_ownership_warm_rejects_first_and_last_positive_overflow_quotients() {
    let index = FileIndex::decode(&WARM_OVERFLOW_IMAGE).unwrap();
    assert_eq!(index.body_coverage(0), Some((-f64::MAX, f64::MAX)));
    // First quotient overflows; second input isolates last-quotient overflow.
    for start in [f64::MAX, -f64::MAX] {
        assert_eq!(
            index
                .warm_ranges(&WARM_OVERFLOW_IMAGE, start, f64::MAX)
                .unwrap_err(),
            Error::NonFiniteArgument {
                value: f64::INFINITY
            },
        );
    }
}

#[test]
fn leb_index_ownership_warm_rejects_negative_first_before_finite_last() {
    let mut image = WARM_OVERFLOW_IMAGE;
    image[16..24].copy_from_slice(&f64::MAX.to_le_bytes());
    image[124..132].copy_from_slice(&f64::MAX.to_le_bytes());
    let index = FileIndex::decode(&image).unwrap();
    assert_eq!(
        index.warm_ranges(&image, -f64::MAX, f64::MAX).unwrap_err(),
        Error::NonFiniteArgument {
            value: f64::NEG_INFINITY
        },
    );
}

#[test]
fn leb_index_ownership_warm_preserves_overlap_input_and_finite_segment_rules() {
    let mut image = WARM_OVERFLOW_IMAGE;
    for offset in [16, 124] {
        image[offset..offset + 8].copy_from_slice(&1000.0_f64.to_le_bytes());
    }
    for offset in [24, 132] {
        image[offset..offset + 8].copy_from_slice(&1100.0_f64.to_le_bytes());
    }
    image[140..148].copy_from_slice(&100.0_f64.to_le_bytes());
    let index = FileIndex::decode(&image).unwrap();
    assert_eq!(
        index.warm_ranges(&image, 1001.0, 1002.0).unwrap(),
        vec![(164, 24)]
    );
    // Reversed bounds in the same segment remain a range, without a new guard.
    assert_eq!(
        index.warm_ranges(&image, 1002.0, 1001.0).unwrap(),
        vec![(164, 24)]
    );
    assert!(
        index
            .warm_ranges(&image, f64::MAX, f64::MAX)
            .unwrap()
            .is_empty()
    );
    for (start, end, expected) in [
        (f64::NAN, 1001.0, f64::NAN),
        (1001.0, f64::INFINITY, f64::INFINITY),
    ] {
        match index.warm_ranges(&image, start, end).unwrap_err() {
            Error::NonFiniteArgument { value } => assert_eq!(value.to_bits(), expected.to_bits()),
            other => panic!("expected bound refusal, got {other:?}"),
        }
    }
}

#[test]
fn leb_index_ownership_owned_warm_errors_before_page_touch_and_closed_first() {
    let mut file =
        LebFile::from_bytes("warm-overflow.leb".into(), WARM_OVERFLOW_IMAGE.to_vec()).unwrap();
    for error in [
        file.warm_ranges(f64::MAX, f64::MAX).unwrap_err(),
        file.warm(f64::MAX, f64::MAX).unwrap_err(),
    ] {
        assert_eq!(
            error,
            Error::NonFiniteArgument {
                value: f64::INFINITY
            }
        );
    }
    file.close();
    assert_eq!(
        file.warm(f64::NAN, f64::INFINITY).unwrap_err(),
        Error::LebReaderClosed
    );
    assert_eq!(
        file.warm_ranges(f64::MAX, f64::MAX).unwrap_err(),
        Error::LebReaderClosed
    );
}

/// Own one regular fixture inode until all borrowed/mapped handles have dropped.
struct Scratch(PathBuf);
impl Scratch {
    fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "leb-index-ownership-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
}
impl Drop for Scratch {
    fn drop(&mut self) {
        std::fs::remove_dir_all(&self.0).unwrap();
    }
}

#[test]
fn leb_index_ownership_mapped_warm_errors_before_advice_and_closed_first() {
    let scratch = Scratch::new();
    let path = scratch.0.join("warm-overflow.leb");
    std::fs::write(&path, WARM_OVERFLOW_IMAGE).unwrap();
    // SAFETY: this test owns the regular inode, never changes its bytes while
    // a map/clone is alive, and drops the mapped handle before scratch cleanup.
    let mut file = unsafe { LebFile::open_mapped(&path) }.unwrap();
    assert_eq!(
        file.warm(f64::MAX, f64::MAX).unwrap_err(),
        Error::NonFiniteArgument {
            value: f64::INFINITY
        }
    );
    assert_eq!(
        file.warm_ranges(f64::MAX, f64::MAX).unwrap_err(),
        Error::NonFiniteArgument {
            value: f64::INFINITY
        }
    );
    file.close();
    assert_eq!(
        file.warm(f64::NAN, f64::INFINITY).unwrap_err(),
        Error::LebReaderClosed
    );
}
