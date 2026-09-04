# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Regression tests for issue #49: no silent type-21 degradation.

Two compounding defects made a type-21 body answer from the Keplerian
approximation while every signal available to the caller said otherwise.

**Boundary semantics disagreed.** ``calc_spk_body_position`` gated the epoch
end-inclusively against the advertised coverage; the vendored reader's upper
bound is exclusive, and the apparent-place pipeline retards the epoch by the
one-way light time before evaluating it. So two windows *inside* the
advertised coverage could not be evaluated: the last instant, and the first
light-time of the span (~0.11 d for Chiron, over a day at Sedna-class
distances). Both raised inside the reader, were caught, returned None, and
were answered from the Keplerian model — logged at debug level only.

**Registration validated nothing.** The type-21 branch skipped the target
check on the premise that the reader could not list its segments. It can:
``.target``, ``.center`` and ``.data_type`` are all exposed. A mismatched
NAIF ID, or a kernel with a center other than the hardcoded Sun, registered
cleanly and then made *every* evaluation fail into the same silent fallback.

The doubles below reproduce the reader's real contract — exclusive upper
bound, and ``ValueError('Invalid Target and/or Center')`` for an unmatched
pair — so the tests prove the wrapper agrees with it.
"""

from __future__ import annotations

import math
from unittest.mock import patch

import pytest

from libephemeris import spk, state
from libephemeris.constants import CHIRON
from libephemeris.exceptions import EphemerisRangeError

pytestmark = pytest.mark.usefixtures("allow_mocked_network")

_AU_KM = 149597870.7
_C_AU_DAY = 173.1446326846693

_NAIF = 20002060
_START_JD = 2458849.5  # 2020-01-01
_END_JD = 2462502.5  # 2030-01-01
_HELIO_AU = 18.8  # Chiron-like, ~0.11 d of one-way light time


class _Segment:
    def __init__(self, target, center, start_jd, end_jd, data_type=21):
        self.target = target
        self.center = center
        self.start_jd = start_jd
        self.end_jd = end_jd
        self.data_type = data_type


class _Kernel:
    """A type-21 reader double with the vendored reader's real semantics."""

    def __init__(self, target=_NAIF, center=10):
        self.segments = [_Segment(target, center, _START_JD, _END_JD)]
        self.calls: list[tuple[int, int, float]] = []

    def compute_type21(self, center, target, jd1, jd2=0.0):
        self.calls.append((center, target, jd1))
        matched = [
            seg
            for seg in self.segments
            if seg.target == target and seg.center == center
        ]
        if not matched:
            raise ValueError("Invalid Target and/or Center")
        # The upper bound is exclusive, exactly as get_MDA_record has it.
        if jd1 < matched[0].start_jd or jd1 >= matched[-1].end_jd:
            raise ValueError("Invalid Time to evaluate")
        # A fixed heliocentric distance is enough: only its magnitude feeds
        # the light-time margin under test.
        return ([_HELIO_AU * _AU_KM, 0.0, 0.0], [0.0, 0.001, 0.0])


@pytest.fixture
def registered_type21(tmp_path):
    """Register a type-21 double for CHIRON and restore the map afterwards."""
    path = tmp_path / "chiron.bsp"
    path.write_bytes(b"\x00" * 16)
    kernel = _Kernel()
    saved = dict(state._SPK_BODY_MAP)
    state._SPK_BODY_MAP[CHIRON] = (str(path), _NAIF)
    with (
        patch("libephemeris.spk._detect_spk_type", return_value=21),
        patch("libephemeris.spk._load_type21_kernel", return_value=kernel),
        patch(
            "libephemeris.spk.get_spk_coverage",
            return_value=(_START_JD, _END_JD),
        ),
    ):
        yield kernel
    state._SPK_BODY_MAP.clear()
    state._SPK_BODY_MAP.update(saved)


def _expected_margin() -> float:
    return max(0.2, (_HELIO_AU + 1.0) / _C_AU_DAY * 1.2)


# ---------------------------------------------------------------------------
# The advertised coverage is not the usable coverage
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_the_usable_span_excludes_both_light_time_edges():
    kernel = _Kernel()
    usable_start, usable_end = spk._type21_effective_coverage(
        kernel, _NAIF, 10, "unused.bsp", _START_JD + 100.0
    )
    margin = _expected_margin()
    assert usable_start == pytest.approx(_START_JD + margin)
    assert usable_end == pytest.approx(_END_JD - margin)
    # Chiron's light time is about 0.11 d; the floor keeps the margin at 0.2.
    assert margin >= _HELIO_AU / _C_AU_DAY


@pytest.mark.unit
def test_exact_end_boundary_raises_instead_of_degrading(registered_type21):
    """The reader's upper bound is exclusive; the gate was inclusive."""
    ts = state.get_timescale()
    t = ts.tdb_jd(_END_JD)
    with pytest.raises(EphemerisRangeError):
        spk.calc_spk_body_position(t, CHIRON, 0)


@pytest.mark.unit
@pytest.mark.parametrize("offset", [0.001, 0.05, 0.1])
def test_light_time_dead_zone_at_the_start_raises(registered_type21, offset):
    """Inside the advertised span, but the retarded epoch is not."""
    ts = state.get_timescale()
    t = ts.tdb_jd(_START_JD + offset)
    with pytest.raises(EphemerisRangeError):
        spk.calc_spk_body_position(t, CHIRON, 0)


@pytest.mark.unit
def test_the_error_reports_the_bounds_it_enforced(registered_type21):
    ts = state.get_timescale()
    t = ts.tdb_jd(_START_JD + 0.05)
    with pytest.raises(EphemerisRangeError) as excinfo:
        spk.calc_spk_body_position(t, CHIRON, 0)
    error = excinfo.value
    assert error.body_id == CHIRON
    margin = _expected_margin()
    assert error.start_jd == pytest.approx(_START_JD + margin)
    assert error.end_jd == pytest.approx(_END_JD - margin)


@pytest.mark.unit
@pytest.mark.parametrize("offset", [1.0, 100.0, 1000.0])
def test_the_interior_still_computes(registered_type21, offset):
    """Narrowing the span must not cost any epoch that always worked."""
    ts = state.get_timescale()
    t = ts.tdb_jd(_START_JD + offset)
    result = spk.calc_spk_body_position(t, CHIRON, 0)
    assert result is not None
    assert 0.0 <= result[0] < 360.0


@pytest.mark.unit
def test_outside_the_advertised_span_still_raises(registered_type21):
    """The pre-existing contract is unchanged."""
    ts = state.get_timescale()
    for jd in (_START_JD - 10.0, _END_JD + 10.0):
        with pytest.raises(EphemerisRangeError):
            spk.calc_spk_body_position(ts.tdb_jd(jd), CHIRON, 0)


# ---------------------------------------------------------------------------
# Both routes must agree on where the kernel can be evaluated
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize("offset", [0.001, 0.05, 0.1])
def test_the_skyfield_route_declines_exactly_where_the_direct_one_refuses(
    registered_type21, offset
):
    """The two routes shared no code and disagreed at the edges."""
    assert spk.get_spk_type21_target(CHIRON, _START_JD + offset) is None


@pytest.mark.unit
def test_the_skyfield_route_serves_the_interior(registered_type21):
    assert spk.get_spk_type21_target(CHIRON, _START_JD + 100.0) is not None


# ---------------------------------------------------------------------------
# The center is read from the kernel, not assumed
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_the_segment_center_is_read_from_the_kernel():
    assert spk._type21_center(_Kernel(center=10), _NAIF) == 10
    assert spk._type21_center(_Kernel(center=0), _NAIF) == 0
    # No matching segment: the documented Horizons default.
    assert spk._type21_center(_Kernel(target=999), _NAIF) == 10


@pytest.mark.unit
def test_a_non_solar_center_is_used_rather_than_the_hardcoded_sun(tmp_path):
    """Every call site passed 10; a barycentric kernel failed every call."""
    path = tmp_path / "barycentric.bsp"
    path.write_bytes(b"\x00" * 16)
    kernel = _Kernel(center=0)
    saved = dict(state._SPK_BODY_MAP)
    state._SPK_BODY_MAP[CHIRON] = (str(path), _NAIF)
    try:
        with (
            patch("libephemeris.spk._detect_spk_type", return_value=21),
            patch("libephemeris.spk._load_type21_kernel", return_value=kernel),
            patch(
                "libephemeris.spk.get_spk_coverage",
                return_value=(_START_JD, _END_JD),
            ),
        ):
            ts = state.get_timescale()
            result = spk.calc_spk_body_position(ts.tdb_jd(_START_JD + 100.0), CHIRON, 0)
        assert result is not None
        assert kernel.calls, "the kernel was never evaluated"
        assert {center for center, _target, _jd in kernel.calls} == {0}
    finally:
        state._SPK_BODY_MAP.clear()
        state._SPK_BODY_MAP.update(saved)


# ---------------------------------------------------------------------------
# Registration validates the target
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_registration_refuses_a_naif_id_the_kernel_does_not_carry(tmp_path):
    path = tmp_path / "other.bsp"
    path.write_bytes(b"\x00" * 16)
    saved = dict(state._SPK_BODY_MAP)
    state._SPK_BODY_MAP.pop(CHIRON, None)
    try:
        with (
            patch("libephemeris.spk._detect_spk_type", return_value=21),
            patch(
                "libephemeris.spk._load_type21_kernel",
                return_value=_Kernel(target=2000001),
            ),
        ):
            with pytest.raises(ValueError, match="no type 21 segment"):
                spk.register_spk_body(CHIRON, str(path), _NAIF)
        assert CHIRON not in state._SPK_BODY_MAP
    finally:
        state._SPK_BODY_MAP.clear()
        state._SPK_BODY_MAP.update(saved)


@pytest.mark.unit
def test_registration_error_lists_what_the_kernel_does_carry(tmp_path):
    path = tmp_path / "other.bsp"
    path.write_bytes(b"\x00" * 16)
    saved = dict(state._SPK_BODY_MAP)
    try:
        with (
            patch("libephemeris.spk._detect_spk_type", return_value=21),
            patch(
                "libephemeris.spk._load_type21_kernel",
                return_value=_Kernel(target=2000001),
            ),
        ):
            with pytest.raises(ValueError) as excinfo:
                spk.register_spk_body(CHIRON, str(path), _NAIF)
        assert "2000001" in str(excinfo.value)
    finally:
        state._SPK_BODY_MAP.clear()
        state._SPK_BODY_MAP.update(saved)


@pytest.mark.unit
def test_a_reader_exposing_no_segments_is_not_refused(tmp_path):
    """Silence from the reader is not evidence that the target is absent.

    Refusing every id because a reader reports no summaries would be a
    worse failure than the missing check it replaces.
    """

    class _Opaque:
        segments: list = []

        def compute_type21(self, center, target, jd1, jd2=0.0):  # pragma: no cover
            return ([0.0, 0.0, 0.0], [0.0, 0.0, 0.0])

    path = tmp_path / "opaque.bsp"
    path.write_bytes(b"\x00" * 16)
    saved = dict(state._SPK_BODY_MAP)
    try:
        with (
            patch("libephemeris.spk._detect_spk_type", return_value=21),
            patch("libephemeris.spk._load_type21_kernel", return_value=_Opaque()),
        ):
            spk.register_spk_body(CHIRON, str(path), _NAIF)
        assert CHIRON in state._SPK_BODY_MAP
    finally:
        state._SPK_BODY_MAP.clear()
        state._SPK_BODY_MAP.update(saved)


@pytest.mark.unit
def test_registration_accepts_a_matching_target(tmp_path):
    path = tmp_path / "chiron.bsp"
    path.write_bytes(b"\x00" * 16)
    saved = dict(state._SPK_BODY_MAP)
    try:
        with (
            patch("libephemeris.spk._detect_spk_type", return_value=21),
            patch("libephemeris.spk._load_type21_kernel", return_value=_Kernel()),
        ):
            spk.register_spk_body(CHIRON, str(path), _NAIF)
        assert state._SPK_BODY_MAP[CHIRON] == (str(path), _NAIF)
    finally:
        state._SPK_BODY_MAP.clear()
        state._SPK_BODY_MAP.update(saved)


@pytest.mark.unit
def test_the_advertised_coverage_comes_from_the_targets_own_segments():
    """Reading every segment would advertise a span the target lacks."""
    kernel = _Kernel()
    kernel.segments.append(_Segment(999, 10, 0.0, 9e6))
    assert spk._type21_coverage(kernel, _NAIF) == (_START_JD, _END_JD)


@pytest.mark.unit
def test_margin_scales_with_distance():
    """A margin sized for Chiron is far too small for a Sedna-class body."""

    class _Distant(_Kernel):
        def compute_type21(self, center, target, jd1, jd2=0.0):
            super().compute_type21(center, target, jd1)
            return ([500.0 * _AU_KM, 0.0, 0.0], [0.0, 0.0, 0.0])

    near = spk._type21_light_time_margin(_Kernel(), _NAIF, 10, _START_JD + 100.0)
    far = spk._type21_light_time_margin(_Distant(), _NAIF, 10, _START_JD + 100.0)
    assert far > near
    assert far == pytest.approx((500.0 + 1.0) / _C_AU_DAY * 1.2)
    assert not math.isclose(near, far)
