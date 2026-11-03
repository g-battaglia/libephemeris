//! All storage modes use the same ownership and polynomial kernels.
use super::*;
use crate::positions::{SingleLeb, StateProvider};
use std::{
    path::PathBuf,
    sync::atomic::{AtomicUsize, Ordering},
};

struct TemporaryImage(PathBuf);
impl TemporaryImage {
    fn new(bytes: &[u8]) -> Self {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let path = std::env::temp_dir().join(format!(
            "leb-positional-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::write(&path, bytes).unwrap();
        Self(path)
    }
    fn open(&self) -> LebFile {
        LebFile::open_indexed(&self.0, IndexedLimits::default()).unwrap()
    }
}
impl Drop for TemporaryImage {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

fn bits(value: (Vec<f64>, Vec<f64>)) -> (Vec<u64>, Vec<u64>) {
    (
        value.0.into_iter().map(f64::to_bits).collect(),
        value.1.into_iter().map(f64::to_bits).collect(),
    )
}

#[test]
fn every_format_shares_scalar_split_provider_and_composite_words() {
    for bytes in [
        minimal_leb1(),
        leb2_v1_distinct_image(),
        leb2_v2_distinct_image(1.0),
    ] {
        let temporary = TemporaryImage::new(&bytes);
        let indexed = temporary.open();
        let owned = LebFile::from_bytes("test".into(), bytes).unwrap();
        let provider = SingleLeb::from_file(&indexed).unwrap();
        let indexed_composite = CompositeIndex::from_files([&indexed]).unwrap();
        let owned_composite = CompositeIndex::from_files([&owned]).unwrap();
        for body in owned.index().unwrap().bodies() {
            for fraction in [0.0, 0.25, 0.5, 0.75, 1.0] {
                let jd = body.jd_start() + fraction * (body.jd_end() - body.jd_start());
                let expected = bits(owned.eval_body_any(body.body_id(), jd).unwrap());
                assert_eq!(
                    bits(indexed.eval_body_any(body.body_id(), jd).unwrap()),
                    expected
                );
                assert_eq!(
                    bits(provider.eval_raw(body.body_id(), jd).unwrap()),
                    expected
                );
                assert_eq!(
                    bits(indexed_composite.eval_body_any(body.body_id(), jd).unwrap()),
                    expected
                );
                let epoch = TwoPartEpoch {
                    hi: jd,
                    lo: if fraction == 1.0 { -1e-18 } else { 1e-18 },
                };
                assert_eq!(
                    bits(
                        indexed_composite
                            .eval_body_any_two_part(body.body_id(), epoch)
                            .unwrap()
                    ),
                    bits(
                        owned_composite
                            .eval_body_any_two_part(body.body_id(), epoch)
                            .unwrap()
                    )
                );
            }
        }
        assert!(indexed.warm(0.0, 1.0).is_ok());
        assert!(indexed.warm(f64::NAN, 1.0).is_err());
        let mut closed = indexed.clone();
        closed.close();
        assert!(matches!(closed.index(), Err(Error::LebReaderClosed)));
        assert!(!indexed.is_closed());
    }
}

#[test]
fn atomic_replacement_preserves_old_descriptor_and_new_open_has_new_words() {
    let temporary = TemporaryImage::new(&leb2_v2_scaled_image(1.0, 1.0));
    let old = temporary.open();
    let expected = bits(old.eval_body_any(0, 0.25).unwrap());
    let replacement = TemporaryImage::new(&leb2_v2_scaled_image(1.0, 2.0));
    std::fs::rename(&replacement.0, &temporary.0).unwrap();
    let new = temporary.open();
    old.cool();
    assert_eq!(bits(old.eval_body_any(0, 0.25).unwrap()), expected);
    assert_ne!(bits(new.eval_body_any(0, 0.25).unwrap()), expected);
    assert_ne!(old.sha256().unwrap(), new.sha256().unwrap());
}

#[test]
fn truncation_is_a_reader_error_and_failed_digest_is_not_cached() {
    let temporary = TemporaryImage::new(&leb2_v2_distinct_image(1.0));
    let file = temporary.open();
    std::fs::OpenOptions::new()
        .write(true)
        .open(&temporary.0)
        .unwrap()
        .set_len(64)
        .unwrap();
    assert!(matches!(
        file.eval_body_any(0, 0.25),
        Err(Error::LebCorrupt { .. })
    ));
    assert!(matches!(file.sha256(), Err(Error::LebCorrupt { .. })));
    assert!(matches!(file.sha256(), Err(Error::LebCorrupt { .. })));
}

#[test]
fn eight_threads_share_cached_words_without_mixing_images() {
    let a = TemporaryImage::new(&leb2_v2_scaled_image(1.0, 1.0));
    let b = TemporaryImage::new(&leb2_v2_scaled_image(1.0, 2.0));
    let images = [a.open(), b.open()];
    let expected = images
        .each_ref()
        .map(|file| bits(file.eval_body_any(0, 0.25).unwrap()));
    for file in &images {
        file.cool();
    }
    std::thread::scope(|scope| {
        for _ in 0..8 {
            scope.spawn(|| {
                for _ in 0..32 {
                    for (file, expected) in images.iter().zip(&expected) {
                        assert_eq!(bits(file.eval_body_any(0, 0.25).unwrap()), *expected);
                    }
                }
            });
        }
    });
}

#[test]
fn advisory_warm_propagates_resource_refusal_and_can_retry() {
    let temporary = TemporaryImage::new(&leb2_v2_distinct_image(1.0));
    let file = temporary.open();
    let LebSource::Positional(source) = file.source().unwrap() else {
        unreachable!()
    };
    let resources = &source.resources;
    // Occupy the remaining shared metadata allowance as another image could.
    let held = resources
        .metadata
        .reserve(resources.limits.metadata_bytes - resources.metadata.used())
        .unwrap();
    assert!(matches!(
        file.warm(0.0, 1.0),
        Err(Error::ResourceLimit { .. })
    ));
    drop(held);
    assert!(file.warm(0.0, 1.0).is_ok());
}
