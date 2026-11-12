//! Bit-pattern and admission checks for the fused reader reconstruction.
use super::*;

/// Build coefficient-major byte lanes independently from the reader visitor.
/// The raw segment-major words are the oracle, including arbitrary IEEE bits.
fn shuffled_words(words: &[u64], segments: u32, degree: u32, components: u32) -> Vec<u8> {
    let per_segment = (degree as usize + 1) * components as usize;
    assert_eq!(words.len(), segments as usize * per_segment);
    let mut shuffled = Vec::with_capacity(words.len() * SHUFFLE_LANES);
    for lane in 0..SHUFFLE_LANES {
        for coefficient in 0..per_segment {
            for segment in 0..segments as usize {
                shuffled.push(words[segment * per_segment + coefficient].to_le_bytes()[lane]);
            }
        }
    }
    shuffled
}

fn encode_words(words: &[u64], segments: u32, degree: u32, components: u32) -> Vec<u8> {
    zstd::encode_all(
        shuffled_words(words, segments, degree, components).as_slice(),
        1,
    )
    .unwrap()
}

fn encoded_frame(payload: &[u8], content_size: bool, checksum: bool) -> Vec<u8> {
    use std::io::Write as _;
    let mut encoder = zstd::stream::Encoder::new(Vec::new(), 1).unwrap();
    encoder.window_log(10).unwrap();
    encoder.include_contentsize(content_size).unwrap();
    encoder
        .set_pledged_src_size(content_size.then_some(payload.len() as u64))
        .unwrap();
    encoder.include_checksum(checksum).unwrap();
    encoder.write_all(payload).unwrap();
    encoder.finish().unwrap()
}

/// One final RAW block with no content-size field. The advertised window is
/// intentionally independent of the small payload, as permitted by zstd.
fn raw_frame(payload: &[u8], window_log: u8) -> Vec<u8> {
    assert!(payload.len() < 1 << 17);
    assert!((10..=31).contains(&window_log));
    let mut frame = vec![0x28, 0xb5, 0x2f, 0xfd, 0, (window_log - 10) << 3];
    let block_header = ((payload.len() as u32) << 3) | 1;
    frame.extend_from_slice(&block_header.to_le_bytes()[..3]);
    frame.extend_from_slice(payload);
    frame
}

fn skippable_frame(payload: &[u8]) -> Vec<u8> {
    let mut frame = 0x184d_2a50_u32.to_le_bytes().to_vec();
    frame.extend_from_slice(&(payload.len() as u32).to_le_bytes());
    frame.extend_from_slice(payload);
    frame
}

fn decoded_bits(
    frame: &[u8],
    size: u64,
    segments: u32,
    degree: u32,
    components: u32,
    window: Option<u32>,
) -> Result<Vec<u64>> {
    decompress_body_words_limited(frame, size, segments, degree, components, window)
        .map(|words| words.into_iter().map(f64::to_bits).collect())
}

#[test]
fn fused_decode_preserves_full_words_and_independent_layouts() {
    // Include both zero signs, infinities, quiet/signaling NaN payloads,
    // subnormals, and deliberately asymmetric little-endian byte patterns.
    let patterns = [
        0x0000_0000_0000_0000,
        0x8000_0000_0000_0000,
        0x7ff0_0000_0000_0000,
        0xfff0_0000_0000_0000,
        0x7ff8_1234_5678_9abc,
        0xfff8_fedc_ba98_7654,
        0x7ff0_0000_0000_0001,
        0x0000_0000_0000_0001,
        0x8000_0000_0000_0001,
        0x0010_0000_0000_0000,
        0x0123_4567_89ab_cdef,
        0xfedc_ba98_7654_3210,
    ];
    for segments in [0, 1, 2, 7, 31] {
        for degree in [0, 1, 5, 15] {
            for components in [0, 1, 2, 3, 5] {
                let count = segments as usize * (degree as usize + 1) * components as usize;
                let words: Vec<_> = (0..count)
                    .map(|index| {
                        if index < patterns.len() {
                            patterns[index]
                        } else {
                            // Vary every word so an incorrect permutation
                            // cannot hide behind repeating special values.
                            (index as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15)
                        }
                    })
                    .collect();
                let frame = encode_words(&words, segments, degree, components);
                let size = (count * SHUFFLE_LANES) as u64;
                assert_eq!(
                    decoded_bits(&frame, size, segments, degree, components, None).unwrap(),
                    words,
                    "segments={segments}, degree={degree}, components={components}"
                );
                let expected_bytes: Vec<_> =
                    words.iter().flat_map(|word| word.to_le_bytes()).collect();
                assert_eq!(
                    decompress_body(&frame, size, segments, degree, components).unwrap(),
                    expected_bytes
                );
            }
        }
    }
}

#[test]
fn fused_decode_preserves_all_byte_values_in_every_lane() {
    // 2,048 words put each possible byte at each lane offset.
    let words: Vec<_> = (0..8)
        .flat_map(|lane| (0..=255_u64).map(move |byte| byte << (lane * 8)))
        .collect();
    let frame = encode_words(&words, 32, 15, 4);
    assert_eq!(
        decoded_bits(&frame, (words.len() * 8) as u64, 32, 15, 4, None).unwrap(),
        words
    );
}

#[test]
fn fused_decode_rejects_short_overlong_and_corrupt_frames() {
    for length in [0, 1, 7, 9, 16, 1 << 20] {
        let frame = zstd::encode_all(vec![0xa5; length].as_slice(), 1).unwrap();
        assert!(matches!(
            decoded_bits(&frame, 8, 1, 0, 1, None),
            Err(Error::LebCorrupt { .. })
        ));
        assert!(matches!(
            decompress_body(&frame, 8, 1, 0, 1),
            Err(Error::LebCorrupt { .. })
        ));
    }
    let frame = encode_words(&[0x0123_4567_89ab_cdef], 1, 0, 1);
    for length in [0, 1, frame.len() - 1] {
        assert!(matches!(
            decoded_bits(&frame[..length], 8, 1, 0, 1, None),
            Err(Error::LebCorrupt { .. })
        ));
    }
    assert!(matches!(
        decoded_bits(b"invalid zstd", 8, 1, 0, 1, None),
        Err(Error::LebCorrupt { .. })
    ));
    // Geometry is checked before malformed compressed bytes are consumed.
    for size in [0, 7, 9, u64::MAX] {
        assert!(matches!(
            decoded_bits(b"invalid zstd", size, 1, 0, 1, None),
            Err(Error::LebCorrupt { .. })
        ));
    }
    assert!(matches!(
        decoded_bits(b"invalid zstd", 0, u32::MAX, u32::MAX, u32::MAX, None),
        Err(Error::LebCorrupt { .. })
    ));
}

#[test]
fn fused_decode_retains_window_limits_before_frame_allocation() {
    // A valid zstd prefix advertises 1,152 bytes before any block parsing.
    let descriptor = [0x28, 0xb5, 0x2f, 0xfd, 0, 1];
    for error in [
        decoded_bits(&descriptor, 8, 1, 0, 1, Some(10)).unwrap_err(),
        decompress_body_limited(&descriptor, 8, 1, 0, 1, Some(10)).unwrap_err(),
    ] {
        assert!(matches!(
            error,
            Error::ResourceLimit {
                resource: "LEB zstd window",
                required: 1152,
                limit: 1024,
            }
        ));
    }
    // A single-segment frame's content size is its decoder window.
    let mut single = vec![0x28, 0xb5, 0x2f, 0xfd, 0xe0];
    single.extend_from_slice(&4096_u64.to_le_bytes());
    assert!(matches!(
        decoded_bits(&single, 8, 1, 0, 1, Some(11)),
        Err(Error::ResourceLimit {
            resource: "LEB zstd window",
            required: 4096,
            limit: 2048
        })
    ));
    let words: Vec<_> = (0..512_u64)
        .map(|word| word.wrapping_mul(0x0102_0304_0506_0708))
        .collect();
    let frame = encode_words(&words, 8, 15, 4);
    assert_eq!(
        decoded_bits(&frame, 4096, 8, 15, 4, Some(23)).unwrap(),
        words
    );
    // A malformed prefix stays a decoder error, not a resource refusal.
    assert!(matches!(
        decoded_bits(b"invalid zstd", 8, 1, 0, 1, Some(10)),
        Err(Error::LebCorrupt { .. })
    ));
}

#[test]
fn fused_decode_accepts_complete_frames_with_optional_size_and_checksum() {
    let words: Vec<_> = (0..12_u64)
        .map(|word| word.wrapping_mul(0x9e37_79b9_7f4a_7c15))
        .collect();
    let shuffled = shuffled_words(&words, 2, 1, 3);
    let expected: Vec<_> = words.iter().flat_map(|word| word.to_le_bytes()).collect();
    for content_size in [false, true] {
        for checksum in [false, true] {
            let frame = encoded_frame(&shuffled, content_size, checksum);
            assert_eq!(
                zstd::zstd_safe::get_frame_content_size(&frame).unwrap(),
                content_size.then_some(shuffled.len() as u64)
            );
            for window in [None, Some(10)] {
                assert_eq!(decoded_bits(&frame, 96, 2, 1, 3, window).unwrap(), words);
                assert_eq!(
                    decompress_body_limited(&frame, 96, 2, 1, 3, window).unwrap(),
                    expected
                );
            }
        }
    }
}

#[test]
fn fused_decode_accepts_concatenated_empty_and_skippable_frames() {
    let words = [0x7ff8_0123_4567_89ab, 0x8000_0000_0000_0000, 1, u64::MAX];
    let shuffled = shuffled_words(&words, 2, 0, 2);
    // A frame boundary may split a byte lane; framing does not alter layout.
    let first = encoded_frame(&shuffled[..13], false, true);
    let second = encoded_frame(&shuffled[13..], true, true);
    let empty = encoded_frame(&[], true, true);
    let skip = skippable_frame(b"opaque metadata");
    let expected: Vec<_> = words.iter().flat_map(|word| word.to_le_bytes()).collect();
    for frames in [
        vec![first.as_slice(), second.as_slice()],
        vec![skip.as_slice(), first.as_slice(), second.as_slice()],
        vec![first.as_slice(), skip.as_slice(), second.as_slice()],
        vec![first.as_slice(), second.as_slice(), skip.as_slice()],
        vec![
            empty.as_slice(),
            first.as_slice(),
            second.as_slice(),
            empty.as_slice(),
        ],
    ] {
        let compressed = frames.concat();
        for window in [None, Some(10)] {
            assert_eq!(
                decoded_bits(&compressed, 32, 2, 0, 2, window).unwrap(),
                words
            );
            assert_eq!(
                decompress_body_limited(&compressed, 32, 2, 0, 2, window).unwrap(),
                expected
            );
        }
    }
    for compressed in [empty, skip, skippable_frame(&[])] {
        for window in [None, Some(10)] {
            assert!(
                decoded_bits(&compressed, 0, 0, 0, 1, window)
                    .unwrap()
                    .is_empty()
            );
            assert!(
                decompress_body_limited(&compressed, 0, 0, 0, 1, window)
                    .unwrap()
                    .is_empty()
            );
        }
    }
}

#[test]
fn fused_decode_rejects_checksum_damage_and_unfinished_trailing_frames() {
    let shuffled = shuffled_words(&[0x0123_4567_89ab_cdef], 1, 0, 1);
    let valid = encoded_frame(&shuffled, false, true);
    let mut damaged = valid.clone();
    *damaged.last_mut().unwrap() ^= 1;
    let mut truncated_skip = skippable_frame(b"metadata");
    truncated_skip.pop();
    let mut corrupt = vec![damaged];
    for tail in [vec![0], vec![0x28, 0xb5, 0x2f, 0xfd], truncated_skip] {
        let mut frame = valid.clone();
        frame.extend(tail);
        corrupt.push(frame);
    }
    for frame in corrupt {
        for window in [None, Some(10)] {
            assert!(matches!(
                decoded_bits(&frame, 8, 1, 0, 1, window),
                Err(Error::LebCorrupt { .. })
            ));
            assert!(matches!(
                decompress_body_limited(&frame, 8, 1, 0, 1, window),
                Err(Error::LebCorrupt { .. })
            ));
        }
    }
}

#[test]
fn fused_decode_preserves_window_refusals_in_owned_and_later_frames() {
    // Leaving the direct path is not a refusal: ordinary owned windows larger
    // than 1 MiB remain valid under the original streaming contract.
    let admitted = raw_frame(&[0xa5; 8], 21);
    assert_eq!(
        decoded_bits(&admitted, 8, 1, 0, 1, None).unwrap(),
        [0xa5a5_a5a5_a5a5_a5a5]
    );
    // The original unrestricted streaming decoder refuses windows above its
    // default 128 MiB limit, even if this complete RAW frame has only one word.
    let wide = raw_frame(&[0xa5; 8], 28);
    assert!(matches!(
        decoded_bits(&wide, 8, 1, 0, 1, None),
        Err(Error::LebCorrupt { .. })
    ));
    assert!(matches!(
        decompress_body(&wide, 8, 1, 0, 1),
        Err(Error::LebCorrupt { .. })
    ));

    // A valid first frame (including a skippable one) must not hide a later
    // frame whose window exceeds the positional reader's explicit limit.
    let later = raw_frame(&[0x5a; 8], 11);
    for (first, size) in [
        (raw_frame(&[0xa5; 8], 10), 16),
        (skippable_frame(b"metadata"), 8),
    ] {
        let compressed = [first, later.clone()].concat();
        assert!(matches!(
            decoded_bits(&compressed, size, (size / 8) as u32, 0, 1, Some(10)),
            Err(Error::LebCorrupt { .. })
        ));
        assert!(matches!(
            decompress_body_limited(&compressed, size, (size / 8) as u32, 0, 1, Some(10)),
            Err(Error::LebCorrupt { .. })
        ));
    }
}

#[test]
fn fused_decode_preserves_owned_output_at_the_one_mib_boundary() {
    for size in [(1 << 20) - 8, 1 << 20, (1 << 20) + 8] {
        let segments = size / 8;
        let words: Vec<_> = (0..u64::from(segments))
            .map(|word| word.wrapping_mul(0x9e37_79b9_7f4a_7c15))
            .collect();
        let frame = encoded_frame(&shuffled_words(&words, segments, 0, 1), false, true);
        assert_eq!(
            decoded_bits(&frame, u64::from(size), segments, 0, 1, None).unwrap(),
            words
        );
        let expected: Vec<_> = words.iter().flat_map(|word| word.to_le_bytes()).collect();
        assert_eq!(
            decompress_body(&frame, u64::from(size), segments, 0, 1).unwrap(),
            expected
        );
        // Claimed geometry alone cannot turn a short complete frame into a
        // valid large coefficient block on either side of the owned boundary.
        let short = encoded_frame(&[0; 8], true, false);
        assert!(matches!(
            decoded_bits(&short, u64::from(size), segments, 0, 1, None),
            Err(Error::LebCorrupt { .. })
        ));
    }
}
