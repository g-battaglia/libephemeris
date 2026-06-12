"""Tests for LEB2 compression primitives."""

import numpy as np

from libephemeris.leb_compression import (
    compress_body,
    compute_mantissa_bits,
    decompress_body,
    reorder_coeff_major,
    reorder_segment_major,
    shuffle_bytes,
    truncate_mantissa,
    unshuffle_bytes,
)


class TestShuffleBytes:
    def test_round_trip(self):
        data = np.arange(24, dtype=np.float64).tobytes()
        shuffled = shuffle_bytes(data)
        assert shuffled != data
        unshuffled = unshuffle_bytes(shuffled)
        assert unshuffled == data

    def test_round_trip_random(self):
        rng = np.random.default_rng(42)
        data = rng.random(100).tobytes()
        assert unshuffle_bytes(shuffle_bytes(data)) == data

    def test_empty(self):
        assert shuffle_bytes(b"") == b""
        assert unshuffle_bytes(b"") == b""

    def test_single_float(self):
        data = np.array([3.14], dtype=np.float64).tobytes()
        assert unshuffle_bytes(shuffle_bytes(data)) == data


class TestReorder:
    def test_round_trip(self):
        rng = np.random.default_rng(42)
        segments, degree, components = 100, 13, 3
        data = (
            rng.random(segments * components * (degree + 1))
            .astype(np.float64)
            .tobytes()
        )
        reordered = reorder_coeff_major(data, segments, degree, components)
        assert reordered != data
        restored = reorder_segment_major(reordered, segments, degree, components)
        assert restored == data

    def test_small(self):
        segments, degree, components = 2, 1, 1
        arr = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64)
        data = arr.tobytes()
        reordered = reorder_coeff_major(data, segments, degree, components)
        # Coefficient-major: [c0_seg0, c0_seg1, c1_seg0, c1_seg1]
        expected = np.array([1.0, 3.0, 2.0, 4.0], dtype=np.float64).tobytes()
        assert reordered == expected


class TestTruncateMantissa:
    def test_full_precision(self):
        arr = np.array([1.23456789012345], dtype=np.float64)
        result = truncate_mantissa(arr, 52)
        np.testing.assert_array_equal(result, arr)

    def test_zero_bits(self):
        arr = np.array([1.23456789012345, -5.678], dtype=np.float64)
        result = truncate_mantissa(arr, 0)
        np.testing.assert_array_equal(result, np.zeros(2))

    def test_truncation_reduces_precision(self):
        arr = np.array([1.23456789012345], dtype=np.float64)
        result = truncate_mantissa(arr, 10)
        assert result[0] != arr[0]
        assert abs(result[0] - arr[0]) < 0.001  # still close

    def test_preserves_sign(self):
        arr = np.array([-3.14, 2.71], dtype=np.float64)
        result = truncate_mantissa(arr, 20)
        assert result[0] < 0
        assert result[1] > 0


class TestComputeMantissaBits:
    def test_zero_coefficients(self):
        coeffs = np.zeros((10, 3, 5))
        bits = compute_mantissa_bits(coeffs)
        assert bits == [0, 0, 0, 0, 0]

    def test_large_c0_small_high_order(self):
        rng = np.random.default_rng(42)
        coeffs = np.zeros((10, 3, 4))
        coeffs[:, :, 0] = 1.0  # large c0
        coeffs[:, :, 1] = 1e-3  # medium c1
        coeffs[:, :, 2] = 1e-6  # small c2
        coeffs[:, :, 3] = 1e-12  # tiny c3
        bits = compute_mantissa_bits(coeffs)
        assert bits[0] > bits[1] > bits[2] > bits[3]
        assert bits[0] > 20  # c0 needs many bits

    def test_below_target_zeroed(self):
        coeffs = np.zeros((10, 3, 3))
        coeffs[:, :, 0] = 1.0
        coeffs[:, :, 1] = 1e-12  # below 5e-9 target
        coeffs[:, :, 2] = 1e-15  # way below
        bits = compute_mantissa_bits(coeffs)
        assert bits[0] > 0
        assert bits[1] == 0
        assert bits[2] == 0


class TestCompressDecompress:
    def test_round_trip(self):
        rng = np.random.default_rng(42)
        segments, degree, components = 50, 13, 3
        deg1 = degree + 1
        n_coeffs = components * deg1
        raw = rng.random(segments * n_coeffs).astype(np.float64)

        # Make realistic: high-order coefficients small
        arr = raw.reshape(segments, components, deg1)
        for k in range(deg1):
            arr[:, :, k] *= 10.0 ** (-k)

        raw_bytes = arr.reshape(segments, n_coeffs).tobytes()
        bits = compute_mantissa_bits(arr)

        compressed = compress_body(raw_bytes, segments, degree, components, bits)
        decompressed = decompress_body(
            compressed, len(raw_bytes), segments, degree, components
        )

        assert len(decompressed) == len(raw_bytes)
        # Decompressed should be close to original (lossy)
        original = np.frombuffer(raw_bytes, dtype=np.float64)
        restored = np.frombuffer(decompressed, dtype=np.float64)
        np.testing.assert_allclose(restored, original, atol=1e-5)

    def test_compression_ratio(self):
        rng = np.random.default_rng(42)
        segments, degree, components = 200, 13, 3
        deg1 = degree + 1
        n_coeffs = components * deg1
        arr = np.zeros((segments, components, deg1))
        arr[:, :, 0] = rng.random((segments, components))  # c0 ~ 1.0
        arr[:, :, 1] = rng.random((segments, components)) * 0.01
        # c2-c13 are effectively zero

        raw_bytes = arr.reshape(segments, n_coeffs).tobytes()
        bits = compute_mantissa_bits(arr)
        compressed = compress_body(raw_bytes, segments, degree, components, bits)

        ratio = len(raw_bytes) / len(compressed)
        assert ratio > 3.0  # should compress well with mostly-zero high orders
