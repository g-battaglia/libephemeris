//! Contiguous-image coefficient windows, independent of astronomical data.
use super::*;

fn shuffled(words: &[u64], segments: u32, degree: u32, components: u32) -> Vec<u8> {
    let width = (degree as usize + 1) * components as usize;
    assert_eq!(words.len(), segments as usize * width);
    let mut bytes = Vec::with_capacity(words.len() * 8);
    for lane in 0..8 {
        for coefficient in 0..width {
            for segment in 0..segments as usize {
                bytes.push(words[segment * width + coefficient].to_le_bytes()[lane]);
            }
        }
    }
    bytes
}

fn image(segments: u32, degree: u32, components: u32, words: &[u64]) -> Vec<u8> {
    fn put<const N: usize>(bytes: &mut [u8], offset: usize, value: [u8; N]) {
        bytes[offset..offset + N].copy_from_slice(&value);
    }
    let frame =
        zstd::encode_all(shuffled(words, segments, degree, components).as_slice(), 1).unwrap();
    let body = HEADER_SIZE + SECTION_DIR_SIZE;
    let directory = body + LEB2_BODY_ENTRY_SIZE;
    let chunk = directory + CHUNK_INDEX_HEADER_SIZE;
    let payload = chunk + CHUNK_ENTRY_SIZE;
    let mut bytes = vec![0; payload];
    put(&mut bytes, 0, MAGIC_LEB2);
    put(&mut bytes, 4, 2_u32.to_le_bytes());
    put(&mut bytes, 8, 1_u32.to_le_bytes());
    put(&mut bytes, 12, 1_u32.to_le_bytes());
    put(&mut bytes, 24, f64::from(segments).to_le_bytes());
    put(&mut bytes, 40, 1_u32.to_le_bytes());
    put(&mut bytes, HEADER_SIZE + 8, (body as u64).to_le_bytes());
    put(
        &mut bytes,
        HEADER_SIZE + 16,
        (LEB2_BODY_ENTRY_SIZE as u64).to_le_bytes(),
    );
    put(&mut bytes, body + 8, segments.to_le_bytes());
    put(&mut bytes, body + 20, f64::from(segments).to_le_bytes());
    put(&mut bytes, body + 28, 1.0_f64.to_le_bytes());
    put(&mut bytes, body + 36, degree.to_le_bytes());
    put(&mut bytes, body + 40, components.to_le_bytes());
    put(&mut bytes, body + 44, (directory as u64).to_le_bytes());
    put(
        &mut bytes,
        body + 60,
        (words.len() as u64 * 8).to_le_bytes(),
    );
    put(&mut bytes, directory, 1_u32.to_le_bytes());
    put(&mut bytes, directory + 8, (segments as f32).to_le_bytes());
    put(&mut bytes, chunk + 8, f64::from(segments).to_le_bytes());
    put(&mut bytes, chunk + 16, (payload as u64).to_le_bytes());
    put(&mut bytes, chunk + 24, (frame.len() as u64).to_le_bytes());
    put(
        &mut bytes,
        chunk + 32,
        (words.len() as u64 * 8).to_le_bytes(),
    );
    put(&mut bytes, chunk + 44, segments.to_le_bytes());
    bytes.extend(frame);
    bytes
}

fn bits(words: &[f64]) -> Vec<u64> {
    words.iter().map(|word| word.to_bits()).collect()
}

#[test]
fn windows_reconstruct_exact_words_ranges_and_partial_tails() {
    let patterns = [
        0,
        0x8000_0000_0000_0000,
        1,
        0x8000_0000_0000_0001,
        0x7ff0_0000_0000_0000,
        0xfff0_0000_0000_0000,
        0x7ff8_1234_5678_9abc,
        0x7ff0_0000_0000_0001,
        0x0123_4567_89ab_cdef,
    ];
    for (segments, degree, components) in [(1, 0, 1), (7, 5, 3), (23, 15, 5)] {
        let width = (degree as usize + 1) * components as usize;
        let words: Vec<_> = (0..segments as usize * width)
            .map(|word| patterns[word % patterns.len()] ^ (word as u64 >> 4))
            .collect();
        let frame =
            zstd::encode_all(shuffled(&words, segments, degree, components).as_slice(), 1).unwrap();
        for first in 0..segments {
            let end = (first + 4).min(segments);
            let selected = reader::decompress_body_window_words_limited(
                &frame,
                (words.len() * 8) as u64,
                segments,
                degree,
                components,
                first..end,
                None,
            )
            .unwrap();
            assert_eq!(
                bits(&selected),
                words[first as usize * width..end as usize * width]
            );
        }
    }
    let file = image(5000, 1, 3, &vec![0; 5000 * 6]);
    let index = FileIndex::decode(&file).unwrap();
    let entry = &index.bodies[0];
    assert_eq!(coefficient_window(entry, 1364, 5000), Some(0..1365));
    assert_eq!(coefficient_window(entry, 1365, 5000), Some(1365..2730));
    assert_eq!(coefficient_window(entry, 4999, 5000), Some(4095..5000));
    assert_eq!(coefficient_window(entry, 5000, 5000), None);
    let mut wide = entry.clone();
    wide.components = 8193;
    assert_eq!(coefficient_window(&wide, 0, 5000), None);
    wide.components = 0;
    assert_eq!(coefficient_window(&wide, 0, 5000), None);
}

#[test]
fn windows_validate_the_complete_frame_before_selecting_rows() {
    use std::io::Write as _;
    let words: Vec<_> = (0..24_u64).collect();
    let payload = shuffled(&words, 4, 1, 3);
    let mut encoder = zstd::stream::Encoder::new(Vec::new(), 1).unwrap();
    encoder.include_checksum(true).unwrap();
    encoder.write_all(&payload).unwrap();
    let frame = encoder.finish().unwrap();
    let mut corrupt = frame.clone();
    *corrupt.last_mut().unwrap() ^= 1;
    let short = zstd::encode_all(&payload[..payload.len() - 1], 1).unwrap();
    let long = zstd::encode_all([payload.as_slice(), &[1]].concat().as_slice(), 1).unwrap();
    for invalid in [corrupt.as_slice(), &frame[..frame.len() - 1], &short, &long] {
        let full = reader::decompress_body_words_limited(invalid, 192, 4, 1, 3, None).unwrap_err();
        let window =
            reader::decompress_body_window_words_limited(invalid, 192, 4, 1, 3, 0..1, None)
                .unwrap_err();
        assert_eq!(full.to_string(), window.to_string());
    }
    for (first, end) in [(3, 2), (0, 5)] {
        assert!(matches!(
            reader::decompress_body_window_words_limited(&frame, 192, 4, 1, 3, first..end, None),
            Err(Error::LebCorrupt { .. })
        ));
    }
    let descriptor = [0x28, 0xb5, 0x2f, 0xfd, 0, 1];
    assert!(matches!(
        reader::decompress_body_window_words_limited(&descriptor, 192, 4, 1, 3, 0..1, Some(10),),
        Err(Error::ResourceLimit {
            resource: "LEB zstd window",
            required: 1152,
            limit: 1024,
        })
    ));
}

#[test]
fn contiguous_windows_share_hits_without_changing_rows_or_scalar_results() {
    let words: Vec<_> = (0..30_000_u32)
        .map(|word| (f64::from(word) / 32.0).to_bits())
        .collect();
    let bytes = image(5000, 1, 3, &words);
    let reference = FileIndex::decode(&bytes).unwrap();
    let file = LebFile::from_bytes("window.leb2".into(), bytes.clone()).unwrap();
    let cache = file.chunk_cache().unwrap();
    let token = cache.word_image_identity();
    for segment in [0, 1364, 1365, 2730, 4095, 4999] {
        let jd = f64::from(segment) + 0.25;
        let expected = reference.eval_leb2_body(&bytes, 0, jd).unwrap();
        let actual = file.eval_leb2_body(0, jd).unwrap();
        assert_eq!(bits(&actual.0), bits(&expected.0));
        assert_eq!(bits(&actual.1), bits(&expected.1));
    }
    assert_eq!(cache.stats(), (2, 4, 4));
    assert_eq!(cache.word_image_identity(), token);
    {
        let inner = cache.lock_inner();
        assert!(
            inner
                .entries
                .keys()
                .all(|key| matches!(key, CacheKey::Window { .. }))
        );
        assert!(inner.entries.values().all(|entry| {
            entry.storage.byte_capacity() <= MAX_COEFFICIENT_WINDOW_BYTES
                && matches!(entry.storage, CachedStorage::Dense(_))
        }));
    }
    let shared = LebFile::from_shared_bytes("shared.leb2".into(), Arc::new(bytes.clone())).unwrap();
    assert!(shared.chunk_cache().unwrap().segment_windows);
    assert_eq!(
        shared.eval_leb2_body(0, 5000.0).unwrap(),
        reference.eval_leb2_body(&bytes, 0, 5000.0).unwrap()
    );
    let static_bytes = Box::leak(bytes.clone().into_boxed_slice());
    let static_file = LebFile::from_static_bytes("static.leb2".into(), static_bytes).unwrap();
    assert!(static_file.chunk_cache().unwrap().segment_windows);
    assert_eq!(
        static_file.eval_leb2_body(0, 4999.25).unwrap(),
        shared.eval_leb2_body(0, 4999.25).unwrap()
    );
    let zero = LebFile::from_bytes_with_cache_budget("zero.leb2".into(), bytes.clone(), 0).unwrap();
    assert_eq!(
        zero.eval_leb2_body(0, 1.25).unwrap(),
        file.eval_leb2_body(0, 1.25).unwrap()
    );
    assert_eq!(zero.chunk_cache_stats().unwrap(), (0, 1, 0));
    assert_eq!(zero.chunk_cache_memory().unwrap(), (0, 0));
    // Public raw caches retain the complete original body key and owner.
    let raw = LebChunkCache::new();
    reference
        .with_leb2_segment_coefficients(
            LebSource::Slice(&bytes),
            0,
            &reference.bodies[0],
            0,
            0.25,
            Some(&raw),
            |row| {
                assert_eq!(bits(row), words[..6]);
                Ok(())
            },
        )
        .unwrap();
    assert!(!raw.segment_windows);
    assert!(raw.lock_inner().entries.contains_key(&CacheKey::Body(0, 0)));
    assert!(!LebChunkCache::with_image_byte_budget(1024).segment_windows);
    assert!(
        !LebChunkCache::with_indexed_resources(IndexedResources::new(IndexedLimits::default()))
            .segment_windows
    );
}

#[test]
fn window_entries_keep_existing_concurrency_lru_and_live_borrower_bounds() {
    let cache = LebChunkCache::with_contiguous_window_byte_budget(1024);
    let barrier = std::sync::Barrier::new(8);
    let key = CacheKey::Window {
        body: 1,
        chunk: 2,
        first_segment: 3,
    };
    std::thread::scope(|scope| {
        for _ in 0..8 {
            let cache = &cache;
            let barrier = &barrier;
            scope.spawn(move || {
                let owner = cache
                    .get_or_decode_key(key, || {
                        barrier.wait();
                        Ok(TrackedVec {
                            values: vec![-0.0; 8],
                            _allocation: None,
                        })
                    })
                    .unwrap();
                assert_eq!(bits(&owner), vec![(-0.0_f64).to_bits(); 8]);
            });
        }
    });
    assert_eq!(cache.stats(), (0, 8, 1));
    assert_eq!(cache.memory_usage(), (64, 1024));
    let held = cache
        .get_or_decode_key(key, || panic!("resident window decoded"))
        .unwrap();
    cache.set_byte_budget(0);
    assert_eq!(cache.memory_usage(), (0, 0));
    assert_eq!(bits(&held), vec![(-0.0_f64).to_bits(); 8]);
    cache.set_byte_budget(DEFAULT_CHUNK_CACHE_BYTES);
    for first_segment in 0..17 {
        cache
            .get_or_decode_key(
                CacheKey::Window {
                    body: 1,
                    chunk: 2,
                    first_segment,
                },
                || {
                    Ok(TrackedVec {
                        values: vec![1.0],
                        _allocation: None,
                    })
                },
            )
            .unwrap();
    }
    assert_eq!(cache.stats().2, MAX_CACHED_CHUNKS);
    let before = cache.memory_usage();
    let oversize = cache
        .get_or_decode_key(
            CacheKey::Window {
                body: 1,
                chunk: 2,
                first_segment: u32::MAX,
            },
            || {
                Ok(TrackedVec {
                    values: vec![1.0; MAX_COEFFICIENT_WINDOW_BYTES / 8 + 1],
                    _allocation: None,
                })
            },
        )
        .unwrap();
    assert_eq!(oversize.len(), MAX_COEFFICIENT_WINDOW_BYTES / 8 + 1);
    assert_eq!(cache.memory_usage(), before);
    assert_eq!(cache.stats().2, MAX_CACHED_CHUNKS);
}
