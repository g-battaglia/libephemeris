# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
LEB2 compression primitives.

Implements error-bounded lossy compression for Chebyshev polynomial coefficients:

1. Mantissa truncation — zeros out unneeded mantissa bits per coefficient order,
   keeping only the bits required for the target precision (0.001" = 5e-9 AU).
2. Coefficient-major reorder — transposes (segments, coeffs) to (coeffs, segments)
   so that same-order coefficients are contiguous (better compression).
3. Byte shuffle — transposes byte lanes across float64 values so exponent bytes
   cluster together (standard technique from HDF5/Blosc).
4. zstd compression — the shuffled data with zeroed mantissa bits compresses well.

The pipeline is fully reversible: decompress -> unshuffle -> inverse reorder.
Truncated floats are valid float64 values — no special handling needed at eval time.
"""

from __future__ import annotations

import threading

import numpy as np
import zstandard as zstd

from .exceptions import LEBCorruptionError
from .exotic_bodies import target_au_map as _exotic_target_au

# Default target: 0.001 arcsecond expressed in AU
# 0.001" = 4.848e-9 radians; at 1 AU distance = 4.848e-9 AU positional error
DEFAULT_TARGET_AU = 5e-9

# Per-body precision targets based on minimum geocentric distance.
# The angular error from a position error delta_AU at distance d_AU is:
#   delta_arcsec = delta_AU / d_AU * 206265
# Bodies closer to Earth need tighter AU precision for the same angular precision.
# Moon/Earth positions also feed into the pipeline's light-time, deflection,
# and aberration corrections for ALL other bodies, further amplifying errors.
BODY_TARGET_AU: dict[int, float] = {
    1: 1e-12,  # Moon    — d_geo ~ 0.002 AU, need <0.001" through full pipeline
    14: 1e-12,  # Earth   — used in corrections for every other body
    0: 1e-10,  # Sun     — deflector for all bodies, high leverage
    2: 1e-10,  # Mercury — d_geo ~ 0.55 AU, fast orbit amplifies errors
    3: 1e-10,  # Venus   — d_geo ~ 0.27 AU, closest inner planet to Earth
    4: 1e-10,  # Mars    — d_geo ~ 0.37 AU
}

# Exotic minor bodies with tighter-than-default targets (close-approach NEAs).
# Registry is the single source of truth: libephemeris.exotic_bodies.
BODY_TARGET_AU.update(_exotic_target_au())


# ---------------------------------------------------------------------------
# Mantissa truncation
# ---------------------------------------------------------------------------


def compute_mantissa_bits(
    coeffs: np.ndarray,
    target_precision: float = DEFAULT_TARGET_AU,
) -> list[int]:
    """Compute the minimum mantissa bits needed per coefficient order.

    Args:
        coeffs: Array of shape (segment_count, components, degree+1) with float64.
        target_precision: Maximum acceptable absolute error per coefficient (AU).

    Returns:
        List of length degree+1 with the number of mantissa bits to keep for
        each coefficient order. 0 means the coefficient can be zeroed entirely.
    """
    deg1 = coeffs.shape[2]
    bits = []
    for k in range(deg1):
        mx = float(np.max(np.abs(coeffs[:, :, k])))
        if mx < target_precision or mx < 1e-30:
            bits.append(0)
        else:
            rel = target_precision / mx
            if rel >= 1.0:
                bits.append(1)
            else:
                b = max(1, int(np.ceil(-np.log2(rel))))
                bits.append(min(b, 52))
    return bits


def truncate_mantissa(arr: np.ndarray, keep_bits: int) -> np.ndarray:
    """Zero out mantissa bits beyond `keep_bits` most significant ones.

    Args:
        arr: Float64 array (any shape).
        keep_bits: Number of mantissa bits to preserve (0-52).
            0 = zero the entire value. 52 = keep all bits (no-op).

    Returns:
        New float64 array with truncated mantissa.
    """
    if keep_bits >= 52:
        return arr.copy()
    if keep_bits <= 0:
        return np.zeros_like(arr)
    shift = np.uint64(52 - keep_bits)
    mask = np.uint64(np.iinfo(np.uint64).max) << shift
    u = arr.view(np.uint64).copy()
    u &= mask
    return u.view(np.float64)


def apply_truncation(
    coeffs: np.ndarray,
    bits_per_order: list[int],
) -> np.ndarray:
    """Apply mantissa truncation to a coefficient array.

    Args:
        coeffs: Shape (segment_count, components, degree+1), float64.
        bits_per_order: Bits to keep per coefficient order (from compute_mantissa_bits).

    Returns:
        New array with truncated mantissa values.
    """
    result = coeffs.copy()
    for k, bits in enumerate(bits_per_order):
        result[:, :, k] = truncate_mantissa(coeffs[:, :, k], bits)
    return result


# ---------------------------------------------------------------------------
# Coefficient-major reorder
# ---------------------------------------------------------------------------


def reorder_coeff_major(
    data: bytes,
    segment_count: int,
    degree: int,
    components: int,
) -> bytes:
    """Reorder from segment-major to coefficient-major layout.

    Input: [seg0: c0_x, c1_x, ..., c0_y, ...] [seg1: ...] ...
    Output: [c0_x_seg0, c0_x_seg1, ...] [c0_y_seg0, ...] [c1_x_seg0, ...]
    """
    deg1 = degree + 1
    n_coeffs = components * deg1
    arr = np.frombuffer(data, dtype=np.float64).reshape(segment_count, n_coeffs)
    return np.ascontiguousarray(arr.T).tobytes()


def reorder_segment_major(
    data: bytes,
    segment_count: int,
    degree: int,
    components: int,
) -> bytes:
    """Inverse of reorder_coeff_major: coefficient-major -> segment-major."""
    deg1 = degree + 1
    n_coeffs = components * deg1
    arr = np.frombuffer(data, dtype=np.float64).reshape(n_coeffs, segment_count)
    return np.ascontiguousarray(arr.T).tobytes()


# ---------------------------------------------------------------------------
# Byte shuffle (HDF5/Blosc-style)
# ---------------------------------------------------------------------------


def shuffle_bytes(data: bytes, element_size: int = 8) -> bytes:
    """Transpose byte lanes: (N_floats x 8 bytes) -> (8 lanes x N bytes)."""
    arr = np.frombuffer(data, dtype=np.uint8).reshape(-1, element_size)
    return np.ascontiguousarray(arr.T).tobytes()


def unshuffle_bytes(data: bytes, element_size: int = 8) -> bytes:
    """Inverse transpose: (8 lanes x N bytes) -> (N_floats x 8 bytes)."""
    n_elements = len(data) // element_size
    arr = np.frombuffer(data, dtype=np.uint8).reshape(element_size, n_elements)
    return np.ascontiguousarray(arr.T).tobytes()


# ---------------------------------------------------------------------------
# Full pipeline
# ---------------------------------------------------------------------------

_COMPRESSOR = zstd.ZstdCompressor(level=19)
_DECOMPRESSOR = zstd.ZstdDecompressor()
# python-zstandard forbids concurrent use of one (de)compressor context;
# the decompressor runs on the process-wide LEB reader from any thread
# calling calc_ut, so it must be serialized.
_DECOMPRESSOR_LOCK = threading.Lock()


def compress_body(
    raw_coeffs: bytes,
    segment_count: int,
    degree: int,
    components: int,
    bits_per_order: list[int],
) -> bytes:
    """Compress a body's Chebyshev coefficients.

    Pipeline: truncation -> coeff-major reorder -> byte shuffle -> zstd.

    Args:
        raw_coeffs: Raw coefficient bytes in segment-major layout.
        segment_count: Number of time segments.
        degree: Chebyshev polynomial degree.
        components: Number of components (always 3).
        bits_per_order: Mantissa bits to keep per coefficient order.

    Returns:
        Compressed bytes blob.
    """
    deg1 = degree + 1
    n_coeffs = components * deg1
    arr = np.frombuffer(raw_coeffs, dtype=np.float64).reshape(
        segment_count, components, deg1
    )
    truncated = apply_truncation(arr, bits_per_order)
    flat = truncated.reshape(segment_count, n_coeffs)
    reordered = np.ascontiguousarray(flat.T).tobytes()
    shuffled = shuffle_bytes(reordered)
    return _COMPRESSOR.compress(shuffled)


def compress_body_chunked(
    raw_coeffs: bytes,
    segment_count: int,
    degree: int,
    components: int,
    bits_per_order: list[int],
    chunk_segments: int,
) -> list[tuple[bytes, int]]:
    """Compress a body's coefficients in temporal chunks.

    Splits the segment-major coefficient data into chunks of
    ``chunk_segments`` segments each, compresses each chunk
    independently.  This enables per-chunk decompression at read time.

    Args:
        raw_coeffs: Raw coefficient bytes in segment-major layout.
        segment_count: Total number of time segments.
        degree: Chebyshev polynomial degree.
        components: Number of components (always 3).
        bits_per_order: Mantissa bits to keep per coefficient order.
        chunk_segments: Number of segments per chunk.

    Returns:
        List of (compressed_bytes, uncompressed_size) tuples, one per chunk.
    """
    deg1 = degree + 1
    n_coeffs = components * deg1
    seg_bytes = n_coeffs * 8
    chunks = []
    offset = 0
    remaining = segment_count

    while remaining > 0:
        n = min(chunk_segments, remaining)
        chunk_raw = raw_coeffs[offset: offset + n * seg_bytes]
        compressed = compress_body(chunk_raw, n, degree, components, bits_per_order)
        chunks.append((compressed, n * seg_bytes))
        offset += n * seg_bytes
        remaining -= n

    return chunks


def decompress_body(
    compressed: bytes,
    uncompressed_size: int,
    segment_count: int,
    degree: int,
    components: int,
) -> bytes:
    """Decompress a body's Chebyshev coefficients.

    Pipeline: zstd -> unshuffle -> segment-major reorder.

    Returns:
        Raw coefficient bytes in segment-major layout (same as LEB1).

    Raises:
        ValueError: If the blob is not valid zstd data (corrupted file).
    """
    try:
        with _DECOMPRESSOR_LOCK:
            decompressed = _DECOMPRESSOR.decompress(
                compressed, max_output_size=uncompressed_size
            )
    except zstd.ZstdError as exc:
        raise LEBCorruptionError(
            f"Corrupted LEB2 data (possibly truncated file): {exc}"
        ) from exc
    if len(decompressed) != uncompressed_size:
        # A valid-but-short frame would otherwise surface as a cryptic
        # numpy reshape error in reorder_segment_major.
        raise LEBCorruptionError(
            f"Corrupted LEB2 data: decompressed to {len(decompressed)} bytes, "
            f"expected {uncompressed_size}"
        )
    unshuffled = unshuffle_bytes(decompressed)
    return reorder_segment_major(unshuffled, segment_count, degree, components)
