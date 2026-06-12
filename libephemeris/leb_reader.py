# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Memory-mapped reader for .leb binary ephemeris files.

This module provides the LEBReader class which opens a .leb file via mmap
and evaluates precomputed Chebyshev polynomials for fast ephemeris lookups.

All evaluation uses pure Python scalar math (no numpy) for optimal
single-point performance (~1.5us per Clenshaw evaluation).
"""

from __future__ import annotations

import mmap
import os
import struct
from bisect import bisect_right
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple, Union

if TYPE_CHECKING:
    from .leb2_reader import LEB2Reader

from .leb_format import (
    BODY_ENTRY_SIZE,
    COORD_ECLIPTIC,
    COORD_GEO_ECLIPTIC,
    COORD_HELIO_ECL,
    DELTA_T_ENTRY_FMT,
    DELTA_T_ENTRY_SIZE,
    DELTA_T_HEADER_FMT,
    DELTA_T_HEADER_SIZE,
    HEADER_SIZE,
    LEB2_MAGIC,
    MAGIC,
    NUTATION_HEADER_SIZE,
    SECTION_BODY_INDEX,
    SECTION_DELTA_T,
    SECTION_DIR_SIZE,
    SECTION_NUTATION,
    SECTION_STARS,
    STAR_ENTRY_SIZE,
    VERSION,
    BodyEntry,
    NutationHeader,
    SectionEntry,
    StarEntry,
    read_body_entry,
    read_header,
    read_nutation_header,
    read_section_dir,
    read_star_entry,
    _madvise_ranges,
    _madvise_dontneed,
)


# =============================================================================
# CLENSHAW ALGORITHM (Pure Python)
# =============================================================================


def _clenshaw(coeffs: tuple, tau: float) -> float:
    """Evaluate a Chebyshev series at tau in [-1, 1] using Clenshaw algorithm.

    Args:
        coeffs: Chebyshev coefficients (c0, c1, ..., cN).
        tau: Evaluation point in [-1, 1].

    Returns:
        The value of the Chebyshev series at tau.
    """
    n = len(coeffs) - 1
    if n == 0:
        return coeffs[0]

    b_k1 = 0.0  # b_{k+1}
    b_k2 = 0.0  # b_{k+2}
    two_tau = 2.0 * tau
    for k in range(n, 0, -1):
        b_k = coeffs[k] + two_tau * b_k1 - b_k2
        b_k2 = b_k1
        b_k1 = b_k
    return coeffs[0] + tau * b_k1 - b_k2


def _deriv_coeffs(coeffs: tuple) -> tuple:
    """Compute derivative Chebyshev coefficients using the standard recurrence.

    Given f(x) = sum_{k=0}^{n} c_k T_k(x), the derivative f'(x) = sum_{k=0}^{n-1} d_k T_k(x)
    where the d_k satisfy:
        d_{n-1} = 2*n * c_n
        d_k = d_{k+2} + 2*(k+1) * c_{k+1}    for k = n-2, ..., 1
        d_0 = d_2/2 + c_1

    Args:
        coeffs: Chebyshev coefficients (c0, c1, ..., cN).

    Returns:
        Derivative coefficients (d0, d1, ..., d_{n-1}).
    """
    n = len(coeffs) - 1
    if n == 0:
        return (0.0,)
    if n == 1:
        return (coeffs[1],)

    d = [0.0] * n  # d[0]..d[n-1]
    d[n - 1] = 2.0 * n * coeffs[n]
    for k in range(n - 2, 0, -1):
        d_k2 = d[k + 2] if k + 2 < n else 0.0
        d[k] = d_k2 + 2.0 * (k + 1) * coeffs[k + 1]
    # d[0] special case (T_0 coefficient halving)
    d_2 = d[2] if n >= 3 else 0.0
    d[0] = d_2 / 2.0 + coeffs[1]

    return tuple(d)


def _clenshaw_derivative(coeffs: tuple, tau: float) -> float:
    """Evaluate the derivative of a Chebyshev series at tau in [-1, 1].

    Computes derivative coefficients via the standard recurrence, then
    evaluates them using Clenshaw.

    Args:
        coeffs: Chebyshev coefficients (c0, c1, ..., cN).
        tau: Evaluation point in [-1, 1].

    Returns:
        d(value)/d(tau) -- derivative with respect to the normalized variable.
    """
    return _clenshaw(_deriv_coeffs(coeffs), tau)


def _clenshaw_with_derivative(coeffs: tuple, tau: float) -> Tuple[float, float]:
    """Evaluate Chebyshev series and its derivative simultaneously at tau.

    Args:
        coeffs: Chebyshev coefficients (c0, c1, ..., cN).
        tau: Evaluation point in [-1, 1].

    Returns:
        (value, derivative) where derivative is d(value)/d(tau).
    """
    value = _clenshaw(coeffs, tau)
    derivative = _clenshaw(_deriv_coeffs(coeffs), tau)
    return value, derivative


# =============================================================================
# LEB READER
# =============================================================================


class LEBReader:
    """Memory-mapped reader for .leb binary ephemeris files.

    Usage:
        reader = LEBReader("/path/to/ephemeris.leb")
        pos, vel = reader.eval_body(SUN, jd_tt)
        reader.close()

    Or as context manager:
        with LEBReader("/path/to/ephemeris.leb") as reader:
            pos, vel = reader.eval_body(SUN, jd_tt)
    """

    def __init__(self, path: str) -> None:
        """Open and parse a .leb file.

        Args:
            path: Path to the .leb file.

        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError: If the file is not a valid .leb file.
        """
        if not os.path.exists(path):
            raise FileNotFoundError(f"LEB file not found: {path}")

        self._path = path
        self._file = open(path, "rb")
        self._mm = mmap.mmap(self._file.fileno(), 0, access=mmap.ACCESS_READ)
        self._eval_cache: Dict[
            Tuple[int, float],
            Tuple[Tuple[float, float, float], Tuple[float, float, float]],
        ] = {}

        try:
            self._parse()
        except (struct.error, IndexError, KeyError) as exc:
            # Truncated/corrupted files surface as struct/index/key errors
            # deep in the unpackers; normalize to ValueError so callers
            # (validation, Skyfield fallback) see one failure type.
            self.close()
            raise ValueError(f"Corrupted LEB file {path!r}: {exc}") from exc
        except (OSError, ValueError):
            # Clean up mmap and file handle on parse failure
            self.close()
            raise

    def _parse(self) -> None:
        """Parse the .leb file header, sections, and metadata.

        Raises:
            ValueError: If the file is not a valid .leb file.
        """
        # Parse header
        self._header = read_header(self._mm, 0)
        if self._header.magic != MAGIC:
            raise ValueError(
                f"Invalid LEB magic: {self._header.magic!r} (expected {MAGIC!r})"
            )
        if self._header.version != VERSION:
            raise ValueError(
                f"Unsupported LEB version: {self._header.version} (expected {VERSION})"
            )

        # Parse section directory
        self._sections: Dict[int, SectionEntry] = {}
        for i in range(self._header.section_count):
            offset = HEADER_SIZE + i * SECTION_DIR_SIZE
            sec = read_section_dir(self._mm, offset)
            self._sections[sec.section_id] = sec

        # Parse body index
        self._bodies: Dict[int, BodyEntry] = {}
        if SECTION_BODY_INDEX in self._sections:
            sec = self._sections[SECTION_BODY_INDEX]
            for i in range(self._header.body_count):
                offset = sec.offset + i * BODY_ENTRY_SIZE
                entry = read_body_entry(self._mm, offset)
                self._bodies[entry.body_id] = entry

        # Parse nutation header
        self._nutation: Optional[NutationHeader] = None
        self._nutation_data_offset: int = 0
        if SECTION_NUTATION in self._sections:
            sec = self._sections[SECTION_NUTATION]
            self._nutation = read_nutation_header(self._mm, sec.offset)
            self._nutation_data_offset = sec.offset + NUTATION_HEADER_SIZE

        # Parse Delta-T table
        self._delta_t_jds: List[float] = []
        self._delta_t_vals: List[float] = []
        if SECTION_DELTA_T in self._sections:
            sec = self._sections[SECTION_DELTA_T]
            n_entries, _ = struct.unpack_from(DELTA_T_HEADER_FMT, self._mm, sec.offset)
            data_offset = sec.offset + DELTA_T_HEADER_SIZE
            for i in range(n_entries):
                off = data_offset + i * DELTA_T_ENTRY_SIZE
                jd, dt = struct.unpack_from(DELTA_T_ENTRY_FMT, self._mm, off)
                self._delta_t_jds.append(jd)
                self._delta_t_vals.append(dt)

        # Parse star catalog
        self._stars: Dict[int, StarEntry] = {}
        if SECTION_STARS in self._sections:
            sec = self._sections[SECTION_STARS]
            n_stars = sec.size // STAR_ENTRY_SIZE
            for i in range(n_stars):
                offset = sec.offset + i * STAR_ENTRY_SIZE
                star = read_star_entry(self._mm, offset)
                self._stars[star.star_id] = star

    def __enter__(self) -> "LEBReader":
        return self

    def __exit__(self, *args) -> None:
        self.close()

    @property
    def path(self) -> str:
        """Return the file path of this .leb file."""
        return self._path

    @property
    def jd_range(self) -> Tuple[float, float]:
        """Return (jd_start, jd_end) covered by this file."""
        return (self._header.jd_start, self._header.jd_end)

    def warm(self, jd_start: float, jd_end: float) -> None:
        """Pre-fault mmap pages for segments covering [jd_start, jd_end].

        Computes which Chebyshev segments overlap the requested JD range for
        each body and calls ``madvise(MADV_WILLNEED)`` on the corresponding
        byte ranges.

        Args:
            jd_start: Start of the Julian Day range to pre-fault.
            jd_end: End of the Julian Day range to pre-fault.
        """
        ranges: list[tuple[int, int]] = []
        for entry in self._bodies.values():
            if entry.jd_end < jd_start or entry.jd_start > jd_end:
                continue
            seg_bytes = (entry.degree + 1) * entry.components * 8
            first_seg = max(0, int((jd_start - entry.jd_start) / entry.interval_days))
            last_seg = min(
                entry.segment_count - 1,
                int((jd_end - entry.jd_start) / entry.interval_days),
            )
            if first_seg > last_seg:
                continue
            offset = entry.data_offset + first_seg * seg_bytes
            size = (last_seg - first_seg + 1) * seg_bytes
            ranges.append((offset, size))

        _madvise_ranges(self._mm, ranges)

    def cool(self) -> None:
        """Advise the kernel that mmap pages can be reclaimed.

        Idempotent and safe to call on closed readers.
        """
        if self._mm is not None:
            _madvise_dontneed(self._mm)

    def has_body(self, body_id: int) -> bool:
        """Check if a body is available in this .leb file."""
        return body_id in self._bodies

    def eval_body(
        self, body_id: int, jd: float
    ) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
        """Evaluate a body's position and velocity at a given Julian Day.

        Args:
            body_id: SE_* body constant (e.g., SUN=0, MOON=1).
            jd: Julian Day in TT (Terrestrial Time).

        Returns:
            ((pos_0, pos_1, pos_2), (vel_0, vel_1, vel_2))

            For ICRS_BARY bodies: ((x, y, z), (vx, vy, vz)) in AU and AU/day.
            For ECLIPTIC bodies: ((lon, lat, dist), (dlon, dlat, ddist)).
            For HELIO_ECL bodies: same as ECLIPTIC.

            Longitude values are already wrapped to [0, 360) for ecliptic bodies.
            Velocity is the analytical Chebyshev derivative.

        Raises:
            KeyError: If body_id is not in this .leb file.
            ValueError: If jd is outside the body's coverage range.
        """
        # Check instance-level eval cache first.  During a multi-body chart
        # calculation at the same jd_tt the observer (Earth) and gravitational
        # deflectors (Sun, Jupiter, Saturn) are re-evaluated for every planet.
        # Caching avoids ~40 redundant Chebyshev evaluations per chart.
        _cache_key = (body_id, jd)
        cached = self._eval_cache.get(_cache_key)
        if cached is not None:
            return cached

        if body_id not in self._bodies:
            raise KeyError(f"Body {body_id} not in LEB file")

        body = self._bodies[body_id]

        # Check range
        if jd < body.jd_start or jd > body.jd_end:
            raise ValueError(
                f"JD {jd} outside range [{body.jd_start}, {body.jd_end}] "
                f"for body {body_id}"
            )

        # Corrupted-header guard: a zero interval or empty segment table
        # would otherwise surface as ZeroDivisionError / garbage reads.
        if body.interval_days <= 0.0 or body.segment_count <= 0:
            raise ValueError(
                f"Corrupted LEB body entry {body_id}: interval_days="
                f"{body.interval_days}, segment_count={body.segment_count}"
            )

        # O(1) segment lookup
        seg_idx = int((jd - body.jd_start) / body.interval_days)
        seg_idx = max(0, min(seg_idx, body.segment_count - 1))

        # Compute tau (map jd to [-1, 1] within segment)
        seg_start = body.jd_start + seg_idx * body.interval_days
        seg_mid = seg_start + 0.5 * body.interval_days
        tau = 2.0 * (jd - seg_mid) / body.interval_days

        # Clamp tau to [-1, 1] for safety
        if tau > 1.0:
            tau = 1.0
        elif tau < -1.0:
            tau = -1.0

        # Read coefficients from mmap (zero-copy)
        deg1 = body.degree + 1
        n_coeffs = body.components * deg1
        byte_offset = body.data_offset + seg_idx * n_coeffs * 8
        coeffs = struct.unpack_from(f"<{n_coeffs}d", self._mm, byte_offset)

        # Evaluate each component via Clenshaw
        pos = []
        vel = []
        scale = (
            2.0 / body.interval_days
        )  # derivative scaling: d/d(jd) = d/d(tau) * 2/interval

        for c in range(body.components):
            comp_coeffs = coeffs[c * deg1 : (c + 1) * deg1]
            val, deriv = _clenshaw_with_derivative(comp_coeffs, tau)
            pos.append(val)
            vel.append(deriv * scale)

        # Wrap longitude for ecliptic-frame bodies (COORD_GEO_ECLIPTIC is
        # reserved/unused but included for format completeness)
        if body.coord_type in (COORD_ECLIPTIC, COORD_HELIO_ECL, COORD_GEO_ECLIPTIC):
            pos[0] = pos[0] % 360.0

        result = tuple(pos), tuple(vel)  # type: ignore[return-value]

        # Store in cache (bounded: clear when exceeding 256 entries to
        # prevent unbounded growth during ephemeris scans)
        if len(self._eval_cache) > 256:
            self._eval_cache.clear()
        self._eval_cache[_cache_key] = result  # type: ignore[assignment]

        return result  # type: ignore[return-value]

    def has_nutation(self) -> bool:
        """Return True if this LEB file contains usable nutation data.

        A nutation section written with ``--skip-aux`` carries a header with
        ``segment_count == 0`` and no coefficients; treating it as present
        would make eval_nutation read the adjacent section's bytes as
        coefficients, silently corrupting every position.
        """
        return self._nutation is not None and self._nutation.segment_count > 0

    def eval_nutation(self, jd_tt: float) -> Tuple[float, float]:
        """Evaluate nutation angles at a given Julian Day.

        Args:
            jd_tt: Julian Day in TT.

        Returns:
            (dpsi, deps) in radians (IAU 2006/2000A).

        Raises:
            ValueError: If nutation data is not available or JD out of range.
        """
        if self._nutation is None or self._nutation.segment_count <= 0:
            raise ValueError("No nutation data in this LEB file")

        nut = self._nutation

        if nut.interval_days <= 0.0:
            raise ValueError("Corrupted LEB nutation header: zero interval")

        if jd_tt < nut.jd_start or jd_tt > nut.jd_end:
            raise ValueError(
                f"JD {jd_tt} outside nutation range [{nut.jd_start}, {nut.jd_end}]"
            )

        # O(1) segment lookup
        seg_idx = int((jd_tt - nut.jd_start) / nut.interval_days)
        seg_idx = max(0, min(seg_idx, nut.segment_count - 1))

        # Compute tau
        seg_start = nut.jd_start + seg_idx * nut.interval_days
        seg_mid = seg_start + 0.5 * nut.interval_days
        tau = 2.0 * (jd_tt - seg_mid) / nut.interval_days
        if tau > 1.0:
            tau = 1.0
        elif tau < -1.0:
            tau = -1.0

        # Read coefficients
        deg1 = nut.degree + 1
        n_coeffs = nut.components * deg1
        seg_size = n_coeffs * 8
        byte_offset = self._nutation_data_offset + seg_idx * seg_size
        coeffs = struct.unpack_from(f"<{n_coeffs}d", self._mm, byte_offset)

        # Evaluate dpsi and deps
        dpsi = _clenshaw(coeffs[0:deg1], tau)
        deps = _clenshaw(coeffs[deg1 : 2 * deg1], tau)

        return dpsi, deps

    def delta_t(self, jd: float) -> float:
        """Get Delta-T (TT - UT1) at a given Julian Day.

        Uses linear interpolation on the sparse table (sufficient for the
        30-day sample spacing).

        Args:
            jd: Julian Day (UT or TT -- the difference is negligible for lookup).

        Returns:
            Delta-T in days.

        Raises:
            ValueError: If no Delta-T data is available.
        """
        if not self._delta_t_jds:
            raise ValueError("No Delta-T data in this LEB file")

        jds = self._delta_t_jds
        vals = self._delta_t_vals
        n = len(jds)

        # Clamp to range
        if jd <= jds[0]:
            return vals[0]
        if jd >= jds[-1]:
            return vals[-1]

        # Binary search for the interval
        idx = bisect_right(jds, jd) - 1
        idx = max(0, min(idx, n - 2))

        # Linear interpolation (fast, sufficient for 30-day spacing)
        span = jds[idx + 1] - jds[idx]
        if span == 0.0:
            return vals[idx]
        t = (jd - jds[idx]) / span
        return vals[idx] + t * (vals[idx + 1] - vals[idx])

    def get_star(self, star_id: int) -> StarEntry:
        """Look up a fixed star record.

        Args:
            star_id: Internal star ID.

        Returns:
            StarEntry with J2000 position, proper motion, etc.

        Raises:
            KeyError: If star_id is not in the catalog.
        """
        if star_id not in self._stars:
            raise KeyError(f"Star {star_id} not in LEB catalog")
        return self._stars[star_id]

    def close(self) -> None:
        """Close the memory-mapped file and release resources."""
        self._eval_cache.clear()
        if self._mm is not None:
            try:
                self._mm.close()
            except (OSError, ValueError, KeyError):
                pass
            self._mm = None  # type: ignore[assignment]
        if self._file is not None:
            try:
                self._file.close()
            except (OSError, ValueError, KeyError):
                pass
            self._file = None  # type: ignore[assignment]


def open_leb(path: str) -> Union["LEBReader", "LEB2Reader"]:
    """Open a .leb or .leb2 file, auto-detecting LEB1 vs LEB2 format.

    Args:
        path: Path to the .leb or .leb2 file.

    Returns:
        LEBReader for LEB1 files, LEB2Reader for LEB2 files.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file has an unknown magic.
    """
    with open(path, "rb") as f:
        magic = f.read(4)
    if magic == MAGIC:
        return LEBReader(path)
    elif magic == LEB2_MAGIC:
        from .leb2_reader import LEB2Reader

        return LEB2Reader(path)
    else:
        raise ValueError(f"Unknown LEB format magic: {magic!r}")
