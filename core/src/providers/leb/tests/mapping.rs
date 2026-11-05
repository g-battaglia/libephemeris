//! Immutable mapped-file contract. Only owned scratch fixtures are written;
//! an inode is never changed while a live map refers to it.
use super::*;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

struct Scratch(PathBuf);
impl Scratch {
    fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "leb-mmap-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
    fn image(&self, bytes: &[u8]) -> LebFile {
        let path = self.0.join("opaque.bin");
        std::fs::write(&path, bytes).unwrap();
        // SAFETY: this fixture is not rewritten until all returned maps drop.
        unsafe { LebFile::open_mapped(&path) }.unwrap()
    }
}
impl Drop for Scratch {
    fn drop(&mut self) {
        std::fs::remove_dir_all(&self.0).unwrap();
    }
}
fn bits(result: (Vec<f64>, Vec<f64>)) -> Vec<u64> {
    result
        .0
        .into_iter()
        .chain(result.1)
        .map(f64::to_bits)
        .collect()
}

#[test]
fn mapped_formats_match_owned_words_errors_and_budget_at_seams() {
    for bytes in [
        leb1_const_image(0, 1.0, 0.0, 2.0),
        leb2_v1_distinct_image(),
        leb2_v2_distinct_image(1.0),
    ] {
        let scratch = Scratch::new();
        let mapped = scratch.image(&bytes);
        let owned = LebFile::from_bytes(mapped.path().into(), bytes.clone()).unwrap();
        assert_eq!(mapped.allocated_image_bytes(), 0);
        assert_eq!(mapped.mapped_file_bytes(), bytes.len());
        assert_eq!(mapped.bytes().unwrap(), owned.bytes().unwrap());
        assert_eq!(mapped, owned);
        for budget in [0, 1, 48, DEFAULT_CHUNK_CACHE_BYTES] {
            mapped.set_chunk_cache_byte_budget(budget).unwrap();
            for jd in [0.0, 0.25, 1_f64.next_down(), 1.0, 1_f64.next_up(), 1.5, 2.0] {
                assert_eq!(
                    bits(mapped.eval_body_any(0, jd).unwrap()),
                    bits(owned.eval_body_any(0, jd).unwrap())
                );
                let epoch = TwoPartEpoch::split(jd).unwrap();
                assert_eq!(
                    bits(mapped.eval_body_any_two_part(0, epoch).unwrap()),
                    bits(owned.eval_body_any_two_part(0, epoch).unwrap())
                );
            }
            assert!(mapped.chunk_cache_memory().unwrap().0 <= budget);
        }
        for (body, jd) in [(0, 9.0), (999, 0.25)] {
            assert_eq!(
                mapped.eval_body_any(body, jd).unwrap_err(),
                owned.eval_body_any(body, jd).unwrap_err()
            );
        }
        assert_eq!(
            mapped.warm(0.0, 1.0).unwrap(),
            owned.warm(0.0, 1.0).unwrap()
        );
        assert!(
            matches!(mapped.warm(f64::NAN,1.0),Err(Error::NonFiniteArgument{value}) if value.is_nan())
        );
        assert!(
            matches!(owned.warm(f64::NAN,1.0),Err(Error::NonFiniteArgument{value}) if value.is_nan())
        );
        mapped.release_file_pages();
        mapped.cool();
        assert_eq!(
            bits(mapped.eval_body_any(0, 0.25).unwrap()),
            bits(owned.eval_body_any(0, 0.25).unwrap())
        );
        assert!(!format!("{mapped:?}").contains("Owned("));
    }
}

#[test]
fn mapped_clones_composites_cooling_and_last_owner_are_independent_of_handle_close() {
    let scratch = Scratch::new();
    let mut mapped = scratch.image(&leb2_v2_distinct_image(1.0));
    let weak = Arc::downgrade(mapped.image.as_ref().unwrap());
    let clone = mapped.clone();
    let pointer = mapped.bytes().unwrap().as_ptr();
    assert_eq!(clone.bytes().unwrap().as_ptr(), pointer);
    let expected = bits(mapped.eval_body_any(0, 0.25).unwrap());
    mapped.close();
    mapped.cool();
    assert_eq!(mapped.mapped_file_bytes(), 0);
    assert_eq!(mapped.warm(0.0, 1.0).unwrap_err(), Error::LebReaderClosed);
    assert!(weak.upgrade().is_some());
    {
        let composite = CompositeIndex::from_files([&clone]).unwrap();
        composite.warm(0.0, 1.0).unwrap();
        composite.cool();
        assert_eq!(clone.chunk_cache_memory().unwrap().0, 0);
        assert!(clone.retained_chunk_index_bytes() > 0);
        std::thread::scope(|scope| {
            for i in 0..4 {
                let clone = &clone;
                let expected = &expected;
                scope.spawn(move || {
                    for _ in 0..20 {
                        if i == 0 {
                            clone.cool();
                        }
                        assert_eq!(bits(clone.eval_body_any(0, 0.25).unwrap()), *expected);
                    }
                });
            }
        });
    }
    drop(clone);
    assert!(weak.upgrade().is_none());
}

#[test]
fn mapped_validation_is_lazy_and_invalid_open_does_not_leave_a_handle() {
    let scratch = Scratch::new();
    // SAFETY: all paths here either fail open or remain untouched while mapped.
    assert!(unsafe { LebFile::open_mapped(&scratch.0) }.is_err());
    assert!(matches!(
        unsafe { LebFile::open_mapped(&scratch.0.join("missing")) },
        Err(Error::LebFileNotFound { .. })
    ));
    for bytes in [Vec::new(), b"short".to_vec(), vec![0; 64]] {
        assert!(std::fs::write(scratch.0.join("bad"), bytes).is_ok());
        assert!(unsafe { LebFile::open_mapped(&scratch.0.join("bad")) }.is_err());
    }
    let mut bytes = leb2_v2_distinct_image(1.0);
    let idx = FileIndex::decode(&bytes).unwrap();
    let offset = idx.find_body(0).unwrap().data_offset() as usize;
    bytes[offset..offset + 4].copy_from_slice(&u32::MAX.to_le_bytes());
    let map = scratch.image(&bytes);
    assert_eq!(map.retained_chunk_index_bytes(), 0);
    assert!(map.eval_body_any(0, 0.25).is_err());
    assert_eq!(map.retained_chunk_index_bytes(), 0);
}

#[test]
fn metadata_inspection_is_bounded_even_for_a_large_sparse_image() {
    use std::io::{Seek, Write};
    let scratch = Scratch::new();
    for bytes in [
        leb1_const_image(0, 1.0, 0.0, 2.0),
        leb2_v1_distinct_image(),
        leb2_v2_distinct_image(1.0),
    ] {
        let path = scratch.0.join("sparse");
        let mut file = std::fs::OpenOptions::new()
            .create(true)
            .truncate(true)
            .write(true)
            .read(true)
            .open(path)
            .unwrap();
        file.write_all(&bytes).unwrap();
        file.set_len(1 << 30).unwrap(); // Before any mapping; holes need no dataset copy.
        assert_eq!(
            FileIndex::read_from_file(&mut file).unwrap(),
            FileIndex::decode(&bytes).unwrap()
        );
        assert!(file.stream_position().unwrap() < 4096);
    }
}

#[test]
fn bounded_index_matches_slice_decoder_for_malformed_metadata() {
    let scratch = Scratch::new();
    let path = scratch.0.join("metadata");
    for original in [
        leb1_const_image(0, 1.0, 0.0, 2.0),
        leb2_v1_distinct_image(),
        leb2_v2_distinct_image(1.0),
    ] {
        let index = FileIndex::decode(&original).unwrap();
        let metadata_end = index.bodies()[0].data_offset() as usize;
        // Mutate every metadata byte, including counts, offsets, sizes, magic,
        // compression and floating-point fields. Never map these mutable fixtures.
        for at in 0..metadata_end {
            for word in [0, 0xff] {
                let mut bytes = original.clone();
                bytes[at] = word;
                std::fs::write(&path, &bytes).unwrap();
                let mut file = std::fs::File::open(&path).unwrap();
                let expected = FileIndex::decode(&bytes);
                let actual = FileIndex::read_from_file(&mut file);
                assert_eq!(actual, expected, "metadata offset={at}, word={word}");
            }
        }
        for len in [0, 4, HEADER_SIZE - 1, HEADER_SIZE, metadata_end - 1] {
            std::fs::write(&path, &original[..len]).unwrap();
            let mut file = std::fs::File::open(&path).unwrap();
            assert!(FileIndex::decode(&original[..len]).is_err());
            assert!(FileIndex::read_from_file(&mut file).is_err());
        }
    }
}

#[cfg(unix)]
#[test]
fn special_files_refuse_without_blocking_and_regular_symlinks_work() {
    use std::os::unix::fs::symlink;
    let scratch = Scratch::new();
    let fifo = scratch.0.join("fifo");
    assert!(
        std::process::Command::new("mkfifo")
            .arg(&fifo)
            .status()
            .unwrap()
            .success()
    );
    assert!(inspect_file(&fifo).is_err());
    // SAFETY: non-regular file must be refused before mapping; no fixture writes.
    assert!(unsafe { LebFile::open_mapped(&fifo) }.is_err());
    let original = scratch.image(&leb1_const_image(0, 1.0, 0.0, 2.0));
    let link = scratch.0.join("link");
    symlink(scratch.0.join("opaque.bin"), &link).unwrap();
    // SAFETY: symlink and its regular target stay immutable while owners exist.
    let linked = unsafe { LebFile::open_mapped(&link) }.unwrap();
    assert_eq!(linked.bytes().unwrap(), original.bytes().unwrap());
    assert!(inspect_file(&link).is_ok());
}

#[cfg(unix)]
#[test]
fn repeated_map_clone_drop_does_not_leak_file_descriptors() {
    use std::os::unix::fs::MetadataExt;
    let scratch = Scratch::new();
    let file = scratch.image(&leb1_const_image(1, 0.0, 0.0, 2.0));
    let metadata = std::fs::metadata(scratch.0.join("opaque.bin")).unwrap();
    let descriptors = || {
        std::fs::read_dir("/dev/fd")
            .unwrap()
            .flatten()
            .filter_map(|entry| std::fs::metadata(entry.path()).ok())
            .filter(|candidate| {
                candidate.dev() == metadata.dev() && candidate.ino() == metadata.ino()
            })
            .count()
    };
    let before = descriptors();
    for _ in 0..32 {
        // SAFETY: the scratch inode is unchanged until all owners below drop.
        let mapped = unsafe { LebFile::open_mapped(&scratch.0.join("opaque.bin")) }.unwrap();
        let clone = mapped.clone();
        assert_eq!(
            mapped.bytes().unwrap().as_ptr(),
            clone.bytes().unwrap().as_ptr()
        );
        drop(mapped);
        assert_eq!(clone.bytes().unwrap().len(), file.bytes().unwrap().len());
        drop(clone);
    }
    assert!(
        descriptors() <= before + 1,
        "descriptor count grew from {before} to {}",
        descriptors()
    );
}

#[cfg(unix)]
#[test]
fn atomic_replacement_keeps_old_inode_and_new_cache_generation_separate() {
    let scratch = Scratch::new();
    let old = scratch.image(&leb2_v2_distinct_image(1.0));
    let original = bits(old.eval_body_any(0, 0.25).unwrap());
    let replacement = scratch.0.join("next");
    std::fs::write(&replacement, leb2_v2_scaled_image(1.0, 2.0)).unwrap();
    std::fs::rename(replacement, scratch.0.join("opaque.bin")).unwrap();
    // SAFETY: rename changed the pathname, not the old inode; the new inode is
    // also never modified while mapped. No truncate-in-place test is valid here.
    let new = unsafe { LebFile::open_mapped(&scratch.0.join("opaque.bin")) }.unwrap();
    assert_eq!(new.retained_chunk_index_bytes(), 0);
    assert_ne!(bits(new.eval_body_any(0, 0.25).unwrap()), original);
    old.release_file_pages();
    assert_eq!(bits(old.eval_body_any(0, 0.25).unwrap()), original);
    assert_ne!(old.bytes().unwrap().as_ptr(), new.bytes().unwrap().as_ptr());
}
