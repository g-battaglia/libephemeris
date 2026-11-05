//! Image-affinity regressions for internal owner-bound evaluation views.

use super::*;
use crate::positions::{SingleLeb, StateProvider};

fn bits((position, velocity): (Vec<f64>, Vec<f64>)) -> (Vec<u64>, Vec<u64>) {
    (
        position.into_iter().map(f64::to_bits).collect(),
        velocity.into_iter().map(f64::to_bits).collect(),
    )
}

#[test]
fn file_single_and_composite_views_cannot_mix_same_key_images() {
    for budget in [0, 1, 48, DEFAULT_CHUNK_CACHE_BYTES] {
        let first = LebFile::from_bytes_with_cache_budget(
            "same-label.leb2".into(),
            leb2_v2_scaled_image(1.0, 1.0),
            budget,
        )
        .unwrap();
        let second = LebFile::from_bytes_with_cache_budget(
            "same-label.leb2".into(),
            leb2_v2_scaled_image(1.0, 2.0),
            budget,
        )
        .unwrap();
        let first_single = SingleLeb::from_file(&first).unwrap();
        let second_single = SingleLeb::from_file(&second).unwrap();
        let first_composite = CompositeIndex::from_files([&first]).unwrap();
        let second_composite = CompositeIndex::from_files([&second]).unwrap();
        let first_count = Arc::strong_count(first.chunk_cache().unwrap());
        let second_count = Arc::strong_count(second.chunk_cache().unwrap());
        for _ in 0..3 {
            for jd in [0.25, 1.0, 1.75] {
                let epoch = TwoPartEpoch {
                    hi: jd,
                    lo: 2.0_f64.powi(-55),
                };
                for (file, single, composite) in [
                    (&first, &first_single, &first_composite),
                    (&second, &second_single, &second_composite),
                ] {
                    let expected = bits(
                        file.index()
                            .unwrap()
                            .eval_leb2_body(file.bytes().unwrap(), 0, jd)
                            .unwrap(),
                    );
                    let expected_split = bits(
                        file.index()
                            .unwrap()
                            .eval_leb2_body_two_part(file.bytes().unwrap(), 0, epoch)
                            .unwrap(),
                    );
                    let view = file.evaluation().unwrap();
                    assert_eq!(bits(view.eval_any(0, jd).unwrap()), expected);
                    assert_eq!(
                        bits(StateProvider::eval_raw(single, 0, jd).unwrap()),
                        expected
                    );
                    assert_eq!(bits(composite.eval_body_any(0, jd).unwrap()), expected);
                    assert_eq!(
                        bits(view.eval_any_two_part(0, epoch).unwrap()),
                        expected_split
                    );
                    assert_eq!(
                        bits(StateProvider::eval_raw_two_part(single, 0, epoch).unwrap()),
                        expected_split
                    );
                    assert_eq!(
                        bits(composite.eval_body_any_two_part(0, epoch).unwrap()),
                        expected_split
                    );
                }
                assert_ne!(
                    bits(first.eval_body_any(0, jd).unwrap()),
                    bits(second.eval_body_any(0, jd).unwrap())
                );
            }
            // Clearing retained words does not reassociate a directory/cache.
            first.cool();
            second.cool();
        }
        assert_eq!(Arc::strong_count(first.chunk_cache().unwrap()), first_count);
        assert_eq!(
            Arc::strong_count(second.chunk_cache().unwrap()),
            second_count
        );
    }
}

#[test]
fn cloned_owner_survives_another_handle_close_and_closed_admission_stays_first() {
    let mut file =
        LebFile::from_bytes("label.leb2".into(), leb2_v2_scaled_image(1.0, 1.0)).unwrap();
    let clone = file.clone();
    let expected = bits(clone.eval_body_any(0, 0.25).unwrap());
    file.close();
    assert!(matches!(file.evaluation(), Err(Error::LebReaderClosed)));
    assert!(matches!(
        file.eval_body_any(-1, f64::NAN),
        Err(Error::LebReaderClosed)
    ));
    assert_eq!(
        bits(clone.evaluation().unwrap().eval_any(0, 0.25).unwrap()),
        expected
    );
}

#[test]
fn grouping_the_opaque_provider_binding_preserves_its_storage_size() {
    // This is an ownership grouping, not a larger retained owner/cache. These
    // are the exact pre-grouping storage components; no public repr is added.
    struct OriginalStorage<'a> {
        _index: &'a FileIndex,
        _source: LebSource<'a>,
        _cache: Arc<LebChunkCache>,
    }
    assert_eq!(
        std::mem::size_of::<SingleLeb<'_>>(),
        std::mem::size_of::<OriginalStorage<'_>>()
    );
    assert_eq!(
        std::mem::align_of::<SingleLeb<'_>>(),
        std::mem::align_of::<OriginalStorage<'_>>()
    );
}
