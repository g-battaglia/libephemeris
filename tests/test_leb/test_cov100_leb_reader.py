# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Branch-coverage tests for ``libephemeris.leb_reader``.

These tests target the remaining uncovered lines and branch arcs in the
LEB reader: parse-time error normalization, absent optional sections,
tau clamping in ``eval_body``/``eval_nutation``, corrupted-header guards,
the ``warm`` empty-overlap case, the Delta-T no-data path, and the
exception-swallowing close path.

All tests build small synthetic ``.leb`` buffers in ``tmp_path`` (no
network, no external fixtures required), then exercise edge branches by
constructing or lightly mutating in-memory state.  Mutating the live
module is deliberately avoided; only per-reader instance state is touched.
"""

from __future__ import annotations

import struct
from typing import Any

import pytest

from libephemeris.constants import SUN
from libephemeris.leb_format import (
    BODY_ENTRY_SIZE,
    DELTA_T_ENTRY_FMT,
    DELTA_T_ENTRY_SIZE,
    DELTA_T_HEADER_FMT,
    DELTA_T_HEADER_SIZE,
    HEADER_SIZE,
    MAGIC,
    NUTATION_HEADER_SIZE,
    SECTION_BODY_INDEX,
    SECTION_CHEBYSHEV,
    SECTION_DELTA_T,
    SECTION_DIR_SIZE,
    SECTION_NUTATION,
    VERSION,
    BodyEntry,
    FileHeader,
    NutationHeader,
    SectionEntry,
    segment_byte_size,
    write_body_entry,
    write_header,
    write_nutation_header,
    write_section_dir,
)
from libephemeris.leb_reader import LEBReader

JD0 = 2451545.0


# =============================================================================
# SYNTHETIC FILE BUILDERS
# =============================================================================


def _write_header_only(path, *, version: int = VERSION, section_count: int = 0) -> str:
    """Write a header-only LEB file (no section directory entries).

    Args:
        path: Destination path object.
        version: LEB format version to stamp in the header.
        section_count: Number of section directory entries to *claim* (the
            body of those entries is omitted, which truncates the file for
            section_count > 0).

    Returns:
        The string path to the written file.
    """
    buf = bytearray(HEADER_SIZE)
    write_header(
        buf,
        FileHeader(
            magic=MAGIC,
            version=version,
            section_count=section_count,
            body_count=0,
            jd_start=JD0,
            jd_end=JD0 + 100.0,
            generation_epoch=0.0,
            flags=0,
        ),
    )
    path.write_bytes(bytes(buf))
    return str(path)


def _build_single_body(
    path,
    *,
    n_seg: int = 3,
    interval: float = 10.0,
    degree: int = 2,
    components: int = 3,
    with_nutation: bool = False,
    with_delta_t: bool = False,
) -> str:
    """Build a small, valid LEB file with one body (SUN).

    Optionally appends a nutation section and/or a Delta-T section so the
    reader exposes mutable ``_nutation`` / ``_delta_t_*`` state.

    Args:
        path: Destination path object.
        n_seg: Number of Chebyshev segments for the body.
        interval: Segment interval in days.
        degree: Chebyshev polynomial degree.
        components: Number of components per segment (3).
        with_nutation: Append a nutation section.
        with_delta_t: Append a Delta-T section.

    Returns:
        The string path to the written file.
    """
    seg_bytes = segment_byte_size(degree, components)
    sections: list[tuple[int, int, int]] = []  # (section_id, offset, size)

    n_sections = 2 + int(with_nutation) + int(with_delta_t)
    sec_dir_off = HEADER_SIZE
    body_idx_off = sec_dir_off + n_sections * SECTION_DIR_SIZE
    cheby_off = body_idx_off + BODY_ENTRY_SIZE
    cheby_size = n_seg * seg_bytes
    cursor = cheby_off + cheby_size

    sections.append((SECTION_BODY_INDEX, body_idx_off, BODY_ENTRY_SIZE))
    sections.append((SECTION_CHEBYSHEV, cheby_off, cheby_size))

    nut_off = nut_data_off = 0
    nut_nseg = nut_deg = nut_comps = 0
    nut_interval = 20.0
    if with_nutation:
        nut_deg, nut_comps, nut_nseg = 2, 2, 2
        nut_data_size = nut_nseg * nut_comps * (nut_deg + 1) * 8
        nut_off = cursor
        nut_data_off = nut_off + NUTATION_HEADER_SIZE
        nut_size = NUTATION_HEADER_SIZE + nut_data_size
        sections.append((SECTION_NUTATION, nut_off, nut_size))
        cursor += nut_size

    dt_off = 0
    dt_entries: list[tuple[float, float]] = []
    if with_delta_t:
        dt_entries = [(JD0, 0.0008), (JD0 + 30.0, 0.0009), (JD0 + 60.0, 0.0010)]
        dt_off = cursor
        dt_size = DELTA_T_HEADER_SIZE + len(dt_entries) * DELTA_T_ENTRY_SIZE
        sections.append((SECTION_DELTA_T, dt_off, dt_size))
        cursor += dt_size

    buf = bytearray(cursor)
    write_header(
        buf,
        FileHeader(
            magic=MAGIC,
            version=VERSION,
            section_count=n_sections,
            body_count=1,
            jd_start=JD0,
            jd_end=JD0 + n_seg * interval,
            generation_epoch=0.0,
            flags=0,
        ),
    )
    for i, (sid, off, size) in enumerate(sections):
        write_section_dir(
            buf, sec_dir_off + i * SECTION_DIR_SIZE, SectionEntry(sid, off, size)
        )

    write_body_entry(
        buf,
        body_idx_off,
        BodyEntry(
            body_id=SUN,
            coord_type=0,
            segment_count=n_seg,
            jd_start=JD0,
            jd_end=JD0 + n_seg * interval,
            interval_days=interval,
            degree=degree,
            components=components,
            data_offset=cheby_off,
        ),
    )
    for s in range(n_seg):
        for ci in range(components * (degree + 1)):
            struct.pack_into(
                "<d", buf, cheby_off + s * seg_bytes + ci * 8, float(ci + 1)
            )

    if with_nutation:
        write_nutation_header(
            buf,
            nut_off,
            NutationHeader(
                jd_start=JD0,
                jd_end=JD0 + nut_nseg * nut_interval,
                interval_days=nut_interval,
                degree=nut_deg,
                components=nut_comps,
                segment_count=nut_nseg,
                reserved=0,
            ),
        )
        stride = nut_comps * (nut_deg + 1) * 8
        for s in range(nut_nseg):
            for ci in range(nut_comps * (nut_deg + 1)):
                struct.pack_into(
                    "<d", buf, nut_data_off + s * stride + ci * 8, 1e-5 * (ci + 1)
                )

    if with_delta_t:
        struct.pack_into(DELTA_T_HEADER_FMT, buf, dt_off, len(dt_entries), 0)
        data_off = dt_off + DELTA_T_HEADER_SIZE
        for i, (jd, dt) in enumerate(dt_entries):
            struct.pack_into(
                DELTA_T_ENTRY_FMT, buf, data_off + i * DELTA_T_ENTRY_SIZE, jd, dt
            )

    path.write_bytes(bytes(buf))
    return str(path)


class _JDStartTrap:
    """Proxy over a body/nutation entry that flips ``jd_start`` per access.

    The first ``jd_start`` read (the range check) returns a *low* value so
    the in-range guard passes; subsequent reads (the segment/tau math)
    return a *high* value so the computed tau falls below -1.0, forcing the
    ``tau < -1.0`` clamp branch.  Every other attribute delegates to the
    wrapped entry unchanged.
    """

    def __init__(self, base: Any, delta: float = 1000.0) -> None:
        object.__setattr__(self, "_base", base)
        object.__setattr__(self, "_delta", delta)
        object.__setattr__(self, "_n", 0)

    @property
    def jd_start(self) -> float:  # type: ignore[override]
        base = object.__getattribute__(self, "_base")
        delta = object.__getattribute__(self, "_delta")
        n = object.__getattribute__(self, "_n") + 1
        object.__setattr__(self, "_n", n)
        return base.jd_start - delta if n == 1 else base.jd_start + delta

    def __getattr__(self, name: str) -> Any:
        return getattr(object.__getattribute__(self, "_base"), name)


# =============================================================================
# __init__ / _parse error normalization
# =============================================================================


class TestParseErrors:
    """Cover parse-time error normalization and version/section branches."""

    @pytest.mark.unit
    def test_truncated_file_raises_value_error(self, tmp_path):
        """A header claiming sections but lacking their bytes -> ValueError.

        Covers the ``except (struct.error, IndexError, KeyError)`` arm that
        normalizes truncated/corrupted files into a single ValueError type
        (lines 196-197).
        """
        path = _write_header_only(tmp_path / "trunc.leb", section_count=2)
        with pytest.raises(ValueError, match="Corrupted LEB file"):
            LEBReader(path)

    @pytest.mark.unit
    def test_unsupported_version_raises(self, tmp_path):
        """Correct magic but wrong version -> ValueError (line 216)."""
        path = _write_header_only(tmp_path / "ver.leb", version=VERSION + 99)
        with pytest.raises(ValueError, match="Unsupported LEB version"):
            LEBReader(path)

    @pytest.mark.unit
    def test_no_optional_sections(self, tmp_path):
        """A header-only file has no body/nutation/delta-t/stars sections.

        Covers the *false* arc of each optional-section guard:
        229->237 (no body index), 239->245 (no nutation),
        247->258 (no delta-t), 259->exit (no stars).
        """
        path = _write_header_only(tmp_path / "empty.leb", section_count=0)
        reader = LEBReader(path)
        try:
            assert reader._bodies == {}
            assert reader._nutation is None
            assert reader._delta_t_jds == []
            assert reader._stars == {}
        finally:
            reader.close()


# =============================================================================
# warm() empty-overlap branch
# =============================================================================


class TestWarmEmptyOverlap:
    """Cover the warm() first_seg > last_seg continue branch (line 305)."""

    @pytest.mark.unit
    def test_warm_first_seg_past_last_seg(self, tmp_path):
        """A range overlapping a clamped single-segment body skips warming."""
        path = _build_single_body(tmp_path / "warm.leb", n_seg=3, interval=10.0)
        reader = LEBReader(path)
        try:
            body = reader._bodies[SUN]
            # Clamp last_seg to 0 while requesting a window whose first_seg
            # lands at segment index 2 -> first_seg (2) > last_seg (0).
            body.segment_count = 1
            start = body.jd_start + 2 * body.interval_days
            reader.warm(start, start + 0.1)  # must not raise
        finally:
            reader.close()


# =============================================================================
# eval_body branches
# =============================================================================


class TestEvalBodyBranches:
    """Cover eval_body corrupted-header guard and tau clamps."""

    @pytest.mark.unit
    def test_corrupted_body_entry_guard(self, tmp_path):
        """Zero interval -> ValueError corrupted-header guard (line 371)."""
        path = _build_single_body(tmp_path / "corrupt.leb")
        reader = LEBReader(path)
        try:
            reader._bodies[SUN].interval_days = 0.0
            with pytest.raises(ValueError, match="Corrupted LEB body entry"):
                reader.eval_body(SUN, JD0 + 5.0)
        finally:
            reader.close()

    @pytest.mark.unit
    def test_tau_clamped_high(self, tmp_path):
        """JD beyond the last segment forces tau > 1.0 clamp (line 387)."""
        path = _build_single_body(tmp_path / "tauhi.leb", n_seg=3, interval=10.0)
        reader = LEBReader(path)
        try:
            body = reader._bodies[SUN]
            # Extend the declared end past the real segment coverage so a
            # query near jd_end clamps to the last segment with tau > 1.
            body.jd_end = (
                body.jd_start + body.segment_count * body.interval_days + 100.0
            )
            pos, vel = reader.eval_body(SUN, body.jd_end - 1.0)
            assert len(pos) == 3 and len(vel) == 3
        finally:
            reader.close()

    @pytest.mark.unit
    def test_tau_clamped_low(self, tmp_path):
        """A jd_start trap forces tau < -1.0 clamp (line 389)."""
        path = _build_single_body(tmp_path / "taulo.leb", n_seg=3, interval=10.0)
        reader = LEBReader(path)
        try:
            real = reader._bodies[SUN]
            reader._bodies[SUN] = _JDStartTrap(real)
            # Range check sees low jd_start (passes); tau math sees high
            # jd_start (segment origin above jd) -> tau < -1.0.
            pos, vel = reader.eval_body(SUN, real.jd_start)
            assert len(pos) == 3 and len(vel) == 3
        finally:
            reader._bodies[SUN] = real
            reader.close()


# =============================================================================
# eval_nutation branches
# =============================================================================


class TestEvalNutationBranches:
    """Cover eval_nutation error guards and tau clamps."""

    @pytest.mark.unit
    def test_no_nutation_data(self, tmp_path):
        """Missing nutation header -> ValueError (line 448)."""
        path = _build_single_body(tmp_path / "nonut.leb", with_nutation=True)
        reader = LEBReader(path)
        try:
            reader._nutation = None
            with pytest.raises(ValueError, match="No nutation data"):
                reader.eval_nutation(JD0 + 5.0)
        finally:
            reader.close()

    @pytest.mark.unit
    def test_zero_interval_nutation(self, tmp_path):
        """Zero nutation interval -> corrupted-header ValueError (line 453)."""
        path = _build_single_body(tmp_path / "nutzero.leb", with_nutation=True)
        reader = LEBReader(path)
        try:
            reader._nutation.interval_days = 0.0
            with pytest.raises(ValueError, match="zero interval"):
                reader.eval_nutation(JD0 + 5.0)
        finally:
            reader.close()

    @pytest.mark.unit
    def test_jd_out_of_nutation_range(self, tmp_path):
        """JD before nutation start -> ValueError out-of-range (line 456)."""
        path = _build_single_body(tmp_path / "nutrange.leb", with_nutation=True)
        reader = LEBReader(path)
        try:
            nut = reader._nutation
            with pytest.raises(ValueError, match="outside nutation range"):
                reader.eval_nutation(nut.jd_start - 5.0)
        finally:
            reader.close()

    @pytest.mark.unit
    def test_nutation_tau_clamped_high(self, tmp_path):
        """JD beyond last nutation segment forces tau > 1.0 (line 469)."""
        path = _build_single_body(tmp_path / "nuttauhi.leb", with_nutation=True)
        reader = LEBReader(path)
        try:
            nut = reader._nutation
            nut.jd_end = (
                nut.jd_start + nut.segment_count * nut.interval_days + 100.0
            )
            dpsi, deps = reader.eval_nutation(nut.jd_end - 1.0)
            assert isinstance(dpsi, float) and isinstance(deps, float)
        finally:
            reader.close()

    @pytest.mark.unit
    def test_nutation_tau_clamped_low(self, tmp_path):
        """A jd_start trap forces tau < -1.0 in eval_nutation (line 471)."""
        path = _build_single_body(tmp_path / "nuttaulo.leb", with_nutation=True)
        reader = LEBReader(path)
        try:
            real = reader._nutation
            reader._nutation = _JDStartTrap(real)
            dpsi, deps = reader.eval_nutation(real.jd_start)
            assert isinstance(dpsi, float) and isinstance(deps, float)
        finally:
            reader._nutation = real
            reader.close()


# =============================================================================
# delta_t no-data branch
# =============================================================================


class TestDeltaTNoData:
    """Cover delta_t empty-table guard (line 502)."""

    @pytest.mark.unit
    def test_no_delta_t_data(self, tmp_path):
        """An empty Delta-T table -> ValueError (line 502)."""
        path = _build_single_body(tmp_path / "nodt.leb", with_delta_t=True)
        reader = LEBReader(path)
        try:
            reader._delta_t_jds = []
            with pytest.raises(ValueError, match="No Delta-T data"):
                reader.delta_t(JD0 + 5.0)
        finally:
            reader.close()


# =============================================================================
# close() exception swallowing
# =============================================================================


class TestCloseSwallowsErrors:
    """Cover close() exception-swallowing arms (lines 547-548, 553-554)."""

    @pytest.mark.unit
    def test_close_ignores_mm_and_file_errors(self, tmp_path):
        """close() swallows errors from both mmap.close and file.close."""
        path = _build_single_body(tmp_path / "close.leb")
        reader = LEBReader(path)

        class _BadMM:
            def close(self) -> None:
                raise OSError("mm close failed")

        class _BadFile:
            def close(self) -> None:
                raise ValueError("file close failed")

        # Release the real handles first so the OS resources are freed,
        # then swap in raising stand-ins to exercise the except arms.
        reader._mm.close()
        reader._file.close()
        reader._mm = _BadMM()  # type: ignore[assignment]
        reader._file = _BadFile()  # type: ignore[assignment]

        reader.close()  # must not raise

        assert reader._mm is None
        assert reader._file is None
