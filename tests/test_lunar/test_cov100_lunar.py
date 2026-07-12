# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Branch-coverage completion tests for :mod:`libephemeris.lunar`.

These tests target the residual uncovered lines and branch arcs of the
lunar node / apogee / perigee module that the main suites do not exercise:

- Import-time ``ImportError`` fallbacks for the optional correction tables.
- Edge / taper / out-of-range branches of the correction interpolators.
- Ephemeris-range discovery fallbacks (segment failures, outer ``except``).
- Boundary-window sampling fallbacks for the osculating apogee and perigee
  samplers (constrained windows, single-sample windows, failed samples,
  zero valid samples, fewer-than-two valid samples).
- The degenerate (circular-orbit) distance branch of the true lunar node.
- The cubic-spline singular-pivot guard and the linear-system pivot swap /
  singular-matrix return.
- The correction-table skip branches of the interpolated apsides and the
  longitude-wrap branches of ``compare_true_lilith_methods``.

The tests favour deterministic inputs and monkeypatch the module's own
helper seams; no network access is required.
"""

from __future__ import annotations

import importlib.util
import math
import sys

import pytest

from libephemeris import lunar
from libephemeris.exceptions import EphemerisRangeError


# ---------------------------------------------------------------------------
# Import-time optional-dependency guards (lines 161-166, 181-190)
# ---------------------------------------------------------------------------


def _exec_lunar_without(deps: tuple[str, ...]):
    """Exec a throwaway copy of lunar.py with given submodules absent.

    Loads the live module source into a fresh module object with the named
    optional dependencies forced to ``None`` in :data:`sys.modules`, so the
    ``try/except ImportError`` blocks take the fallback path. The live
    module is never reloaded.
    """
    spec = importlib.util.spec_from_file_location(
        "libephemeris._lunar_cov_throwaway", lunar.__file__
    )
    module = importlib.util.module_from_spec(spec)
    module.__package__ = "libephemeris"

    saved: dict[str, object] = {}
    for name in deps:
        saved[name] = sys.modules.get(name)
        sys.modules[name] = None  # type: ignore[assignment]
    try:
        spec.loader.exec_module(module)
    finally:
        for name, value in saved.items():
            if value is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = value
    return module


def test_import_falls_back_when_perigee_corrections_missing():
    """Missing ``lunar_corrections`` selects the ImportError fallback."""
    module = _exec_lunar_without(("libephemeris.lunar_corrections",))

    assert module._PERIGEE_CORRECTIONS_AVAILABLE is False
    assert module.PERIGEE_PERTURBATION_CORRECTIONS == ()
    assert module.PERIGEE_CORRECTION_START_YEAR == 0
    assert module.PERIGEE_CORRECTION_END_YEAR == 0
    assert module.PERIGEE_CORRECTION_STEP_YEARS == 2


def test_import_falls_back_when_apse_corrections_missing():
    """Missing ``lunar_apse_corrections`` selects the ImportError fallback."""
    module = _exec_lunar_without(("libephemeris.lunar_apse_corrections",))

    assert module._APSE_CORRECTIONS_AVAILABLE is False
    assert module.APOGEE_CT_JD_START == 0.0
    assert module.APOGEE_CT_STEP_DAYS == 10.0
    assert module.APOGEE_CT_COUNT == 0
    assert module.APOGEE_CORRECTIONS == ()
    assert module.PERIGEE_CT_JD_START == 0.0
    assert module.PERIGEE_CT_STEP_DAYS == 2.0
    assert module.PERIGEE_CT_COUNT == 0
    assert module.PERIGEE_CORRECTIONS == ()


# ---------------------------------------------------------------------------
# _interpolate_perigee_correction (lines 241, 246-249, 251-254, 262)
# ---------------------------------------------------------------------------


def _jd_for_year(year: float) -> float:
    """Inverse of lunar._jd_to_year for a small synthetic table."""
    return (year - 2000.0) * 365.25 + 2451545.0


@pytest.fixture
def synthetic_perigee_table(monkeypatch):
    """Install a tiny deterministic perigee correction table."""
    monkeypatch.setattr(lunar, "_PERIGEE_CORRECTIONS_AVAILABLE", True)
    monkeypatch.setattr(lunar, "PERIGEE_PERTURBATION_CORRECTIONS", (10.0, 20.0, 30.0))
    monkeypatch.setattr(lunar, "PERIGEE_CORRECTION_START_YEAR", 2000)
    monkeypatch.setattr(lunar, "PERIGEE_CORRECTION_END_YEAR", 2010)
    monkeypatch.setattr(lunar, "PERIGEE_CORRECTION_STEP_YEARS", 2)


def test_perigee_correction_unavailable(monkeypatch):
    """Line 241: returns 0.0 when corrections are unavailable."""
    monkeypatch.setattr(lunar, "_PERIGEE_CORRECTIONS_AVAILABLE", False)
    assert lunar._interpolate_perigee_correction(2451545.0) == 0.0


def test_perigee_correction_empty_table(monkeypatch):
    """Line 241: returns 0.0 when the table is empty."""
    monkeypatch.setattr(lunar, "_PERIGEE_CORRECTIONS_AVAILABLE", True)
    monkeypatch.setattr(lunar, "PERIGEE_PERTURBATION_CORRECTIONS", ())
    assert lunar._interpolate_perigee_correction(2451545.0) == 0.0


def test_perigee_correction_below_start_within_taper(synthetic_perigee_table):
    """Line 249: below start but inside the one-year taper returns a scaled
    edge value."""
    value = lunar._interpolate_perigee_correction(_jd_for_year(1999.5))
    assert value == pytest.approx(5.0)  # corr[0] * (1 - 0.5)


def test_perigee_correction_below_start_beyond_taper(synthetic_perigee_table):
    """Line 248: well below start tapers to 0.0."""
    assert lunar._interpolate_perigee_correction(_jd_for_year(1990.0)) == 0.0


def test_perigee_correction_above_end_within_taper(synthetic_perigee_table):
    """Line 254: above end but inside the taper returns a scaled edge value."""
    value = lunar._interpolate_perigee_correction(_jd_for_year(2010.5))
    assert value == pytest.approx(15.0)  # corr[-1] * (1 - 0.5)


def test_perigee_correction_above_end_beyond_taper(synthetic_perigee_table):
    """Line 253: well above end tapers to 0.0."""
    assert lunar._interpolate_perigee_correction(_jd_for_year(2020.0)) == 0.0


def test_perigee_correction_idx_at_or_past_last(synthetic_perigee_table):
    """Line 262: in range but index runs past the last entry returns corr[-1]."""
    # year=2008 -> idx_float=4 -> idx_low=4 >= len-1 (==2)
    assert lunar._interpolate_perigee_correction(_jd_for_year(2008.0)) == 30.0


def test_perigee_correction_normal_interpolation(synthetic_perigee_table):
    """Sanity: an in-range interior date interpolates linearly."""
    # year=2003 -> idx_float=1.5 -> between corr[1]=20 and corr[2]=30 -> 25
    assert lunar._interpolate_perigee_correction(_jd_for_year(2003.0)) == pytest.approx(
        25.0
    )


# ---------------------------------------------------------------------------
# _interpolate_apse_correction (lines 294, 299-303, 305-309)
# ---------------------------------------------------------------------------

_APSE_CORR = (100.0, 200.0, 300.0)
_APSE_START = 2451545.0
_APSE_STEP = 10.0
_APSE_N = 3


def test_apse_correction_empty():
    """Line 294: empty corrections / zero count returns 0.0."""
    assert (
        lunar._interpolate_apse_correction(_APSE_START, _APSE_START, 10.0, (), 0) == 0.0
    )


def test_apse_correction_before_start_within_taper():
    """Line 303: before the table but inside the taper returns a scaled value."""
    value = lunar._interpolate_apse_correction(
        _APSE_START - 5.0, _APSE_START, _APSE_STEP, _APSE_CORR, _APSE_N
    )
    expected = 100.0 * (1.0 - 5.0 / lunar._APSE_EDGE_TAPER_DAYS)
    assert value == pytest.approx(expected)


def test_apse_correction_before_start_beyond_taper():
    """Line 302: well before the table tapers to 0.0."""
    value = lunar._interpolate_apse_correction(
        _APSE_START - 400.0, _APSE_START, _APSE_STEP, _APSE_CORR, _APSE_N
    )
    assert value == 0.0


def test_apse_correction_after_end_within_taper():
    """Line 309: after the table but inside the taper returns a scaled value."""
    # idx_float = 2.5 -> dist_days = 0.5 * 10 = 5
    value = lunar._interpolate_apse_correction(
        _APSE_START + 25.0, _APSE_START, _APSE_STEP, _APSE_CORR, _APSE_N
    )
    expected = 300.0 * (1.0 - 5.0 / lunar._APSE_EDGE_TAPER_DAYS)
    assert value == pytest.approx(expected)


def test_apse_correction_after_end_beyond_taper():
    """Line 308: well after the table tapers to 0.0."""
    value = lunar._interpolate_apse_correction(
        _APSE_START + 500.0, _APSE_START, _APSE_STEP, _APSE_CORR, _APSE_N
    )
    assert value == 0.0


def test_apse_correction_interior_interpolation():
    """Sanity: interior date interpolates linearly."""
    value = lunar._interpolate_apse_correction(
        _APSE_START + 15.0, _APSE_START, _APSE_STEP, _APSE_CORR, _APSE_N
    )
    assert value == pytest.approx(250.0)


# ---------------------------------------------------------------------------
# _get_ephemeris_range fallbacks (lines 2381-2382, 2385, 2388-2389)
# ---------------------------------------------------------------------------

_DEFAULT_RANGE = (2415020.0, 2471184.0)


def test_get_ephemeris_range_all_segments_fail(monkeypatch):
    """Lines 2381-2382 + 2385: every segment raises -> min stays inf ->
    default range."""

    class _FailingSegment:
        def time_range(self, ts):
            raise AttributeError("no time_range")

    class _FakePlanets:
        segments = [_FailingSegment(), _FailingSegment()]

    monkeypatch.setattr(lunar, "get_planets", lambda: _FakePlanets())
    assert lunar._get_ephemeris_range() == _DEFAULT_RANGE


def test_get_ephemeris_range_outer_except(monkeypatch):
    """Lines 2388-2389: accessing ``.segments`` raises -> outer except ->
    default range."""

    class _NoSegments:
        pass

    monkeypatch.setattr(lunar, "get_planets", lambda: _NoSegments())
    assert lunar._get_ephemeris_range() == _DEFAULT_RANGE


# ---------------------------------------------------------------------------
# Osculating apogee/perigee boundary sampling
# Apogee: lines 2443-2444, 2456, 2467, 2491-2493, 2498-2513, 2518-2519
# Perigee: lines 2829-2872, 2891-2893, 2898-2913, 2918-2919
# ---------------------------------------------------------------------------


@pytest.fixture
def stub_apogee_calc(monkeypatch):
    """Replace calc_true_lilith with a deterministic stub; return its setter."""

    def _set(func):
        monkeypatch.setattr(lunar, "calc_true_lilith", func)

    return _set


@pytest.fixture
def stub_perigee_calc(monkeypatch):
    """Replace calc_osculating_perigee with a deterministic stub."""

    def _set(func):
        monkeypatch.setattr(lunar, "calc_osculating_perigee", func)

    return _set


def _ok(jd):
    """A well-behaved sample value."""
    return (jd % 360.0, 0.0, 0.0025)


def test_apogee_near_end_window(monkeypatch, stub_apogee_calc):
    """Lines 2443-2444: target near the ephemeris end shifts window backward."""
    monkeypatch.setattr(lunar, "_get_ephemeris_range", lambda: (0.0, 1000.0))
    stub_apogee_calc(_ok)
    times, lons, lats, eccs, idx = lunar._sample_osculating_apogee_with_fallback(
        950.0, 100.0, 9
    )
    assert len(times) == 9
    assert max(times) <= 1000.0


def test_apogee_constrained_window_reduces_samples(monkeypatch, stub_apogee_calc):
    """Line 2456: a window narrower than the half-window reduces sample count."""
    monkeypatch.setattr(lunar, "_get_ephemeris_range", lambda: (0.0, 50.0))
    stub_apogee_calc(_ok)
    times, *_ = lunar._sample_osculating_apogee_with_fallback(25.0, 100.0, 9)
    assert len(times) == 3


def test_apogee_single_sample_window(monkeypatch, stub_apogee_calc):
    """Line 2467: a single requested sample at a boundary returns just jd_tt."""
    monkeypatch.setattr(lunar, "_get_ephemeris_range", lambda: (0.0, 1000.0))
    stub_apogee_calc(_ok)
    times, *_ = lunar._sample_osculating_apogee_with_fallback(50.0, 100.0, 1)
    assert times == [50.0]


def test_apogee_some_samples_fail(monkeypatch, stub_apogee_calc):
    """Lines 2491-2493 + 2498-2505: some samples raise but >=2 survive."""
    monkeypatch.setattr(lunar, "_get_ephemeris_range", lambda: (0.0, 1000.0))

    def _calc(jd):
        if jd < 450.0:
            raise EphemerisRangeError("out of range")
        return _ok(jd)

    stub_apogee_calc(_calc)
    times, lons, lats, eccs, idx = lunar._sample_osculating_apogee_with_fallback(
        500.0, 100.0, 9
    )
    assert 2 <= len(times) < 9
    assert all(t >= 450.0 for t in times)


def test_apogee_one_valid_sample(monkeypatch, stub_apogee_calc):
    """Lines 2518-2519: exactly one valid sample falls back to the target."""
    monkeypatch.setattr(lunar, "_get_ephemeris_range", lambda: (0.0, 1000.0))

    def _calc(jd):
        if abs(jd - 500.0) > 1.0:
            raise EphemerisRangeError("out of range")
        return _ok(jd)

    stub_apogee_calc(_calc)
    times, lons, lats, eccs, idx = lunar._sample_osculating_apogee_with_fallback(
        500.0, 100.0, 9
    )
    assert times == [500.0]
    assert idx == 0


def test_apogee_no_valid_samples_target_ok(monkeypatch, stub_apogee_calc):
    """Lines 2507-2510: all samples fail but the target date itself succeeds."""
    monkeypatch.setattr(lunar, "_get_ephemeris_range", lambda: (0.0, 1000.0))

    state = {"calls": 0}

    def _calc(jd):
        state["calls"] += 1
        if state["calls"] <= 9:  # the nine sample calls all fail
            raise EphemerisRangeError("out of range")
        return (123.0, 0.0, 0.0025)

    stub_apogee_calc(_calc)
    times, lons, lats, eccs, idx = lunar._sample_osculating_apogee_with_fallback(
        500.0, 100.0, 9
    )
    assert times == [500.0]
    assert lons == [123.0]
    assert idx == 0


def test_apogee_no_valid_samples_target_fails(monkeypatch, stub_apogee_calc):
    """Lines 2511-2513: all samples and the target fail -> re-raise."""
    monkeypatch.setattr(lunar, "_get_ephemeris_range", lambda: (0.0, 1000.0))

    def _calc(jd):
        raise EphemerisRangeError("out of range")

    stub_apogee_calc(_calc)
    with pytest.raises(EphemerisRangeError):
        lunar._sample_osculating_apogee_with_fallback(500.0, 100.0, 9)


def _perig_ok(jd):
    return (jd % 360.0, 0.0, 0.0024)


def test_perigee_near_end_window(monkeypatch, stub_perigee_calc):
    """Perigee mirror: target near the end shifts window backward."""
    monkeypatch.setattr(lunar, "_get_ephemeris_range", lambda: (0.0, 1000.0))
    stub_perigee_calc(_perig_ok)
    times, *_ = lunar._sample_osculating_perigee_with_fallback(950.0, 100.0, 9)
    assert len(times) == 9
    assert max(times) <= 1000.0


def test_perigee_near_start_window(monkeypatch, stub_perigee_calc):
    """Perigee mirror: target near the start shifts window forward."""
    monkeypatch.setattr(lunar, "_get_ephemeris_range", lambda: (0.0, 1000.0))
    stub_perigee_calc(_perig_ok)
    times, *_ = lunar._sample_osculating_perigee_with_fallback(50.0, 100.0, 9)
    assert len(times) == 9
    assert min(times) >= 0.0


def test_perigee_constrained_window(monkeypatch, stub_perigee_calc):
    """Perigee mirror: narrow window reduces sample count."""
    monkeypatch.setattr(lunar, "_get_ephemeris_range", lambda: (0.0, 50.0))
    stub_perigee_calc(_perig_ok)
    times, *_ = lunar._sample_osculating_perigee_with_fallback(25.0, 100.0, 9)
    assert len(times) == 3


def test_perigee_single_sample_window(monkeypatch, stub_perigee_calc):
    """Perigee mirror: single requested sample returns just jd_tt."""
    monkeypatch.setattr(lunar, "_get_ephemeris_range", lambda: (0.0, 1000.0))
    stub_perigee_calc(_perig_ok)
    times, *_ = lunar._sample_osculating_perigee_with_fallback(50.0, 100.0, 1)
    assert times == [50.0]


def test_perigee_some_samples_fail(monkeypatch, stub_perigee_calc):
    """Perigee mirror: lines 2891-2893 + 2898-2905."""
    monkeypatch.setattr(lunar, "_get_ephemeris_range", lambda: (0.0, 1000.0))

    def _calc(jd):
        if jd < 450.0:
            raise EphemerisRangeError("out of range")
        return _perig_ok(jd)

    stub_perigee_calc(_calc)
    times, *_ = lunar._sample_osculating_perigee_with_fallback(500.0, 100.0, 9)
    assert 2 <= len(times) < 9
    assert all(t >= 450.0 for t in times)


def test_perigee_one_valid_sample(monkeypatch, stub_perigee_calc):
    """Perigee mirror: lines 2918-2919."""
    monkeypatch.setattr(lunar, "_get_ephemeris_range", lambda: (0.0, 1000.0))

    def _calc(jd):
        if abs(jd - 500.0) > 1.0:
            raise EphemerisRangeError("out of range")
        return _perig_ok(jd)

    stub_perigee_calc(_calc)
    times, *_ = lunar._sample_osculating_perigee_with_fallback(500.0, 100.0, 9)
    assert times == [500.0]


def test_perigee_no_valid_samples_target_ok(monkeypatch, stub_perigee_calc):
    """Perigee mirror: lines 2908-2910 (all fail, target ok)."""
    monkeypatch.setattr(lunar, "_get_ephemeris_range", lambda: (0.0, 1000.0))

    state = {"calls": 0}

    def _calc(jd):
        state["calls"] += 1
        if state["calls"] <= 9:
            raise EphemerisRangeError("out of range")
        return (200.0, 0.0, 0.0024)

    stub_perigee_calc(_calc)
    times, lons, *_ = lunar._sample_osculating_perigee_with_fallback(500.0, 100.0, 9)
    assert times == [500.0]
    assert lons == [200.0]


def test_perigee_no_valid_samples_target_fails(monkeypatch, stub_perigee_calc):
    """Perigee mirror: lines 2911-2913 (all fail, target fails -> re-raise)."""
    monkeypatch.setattr(lunar, "_get_ephemeris_range", lambda: (0.0, 1000.0))

    def _calc(jd):
        raise EphemerisRangeError("out of range")

    stub_perigee_calc(_calc)
    with pytest.raises(EphemerisRangeError):
        lunar._sample_osculating_perigee_with_fallback(500.0, 100.0, 9)


# ---------------------------------------------------------------------------
# Degenerate (circular-orbit) node distance (line 2016)
# ---------------------------------------------------------------------------


def test_true_node_degenerate_distance(monkeypatch):
    """Line 2016: a circular orbit (n_mag == 0) uses the mean distance."""
    gm = lunar.GM_EARTH_MOON_AU3_DAY2
    d = 384400.0 / 149597870.7
    v_circular = math.sqrt(gm / d)

    class _Vector:
        def __init__(self, values):
            self._values = values

        @property
        def au(self):
            return self._values

        @property
        def au_per_d(self):
            return self._values

    class _Position:
        def frame_xyz(self, frame):
            return _Vector([d, 0.0, 0.0])

        def frame_xyz_and_velocity(self, frame):
            # Velocity perpendicular to r, in the xy-plane: angular momentum
            # is purely along +z, so n_mag = sqrt(h_x^2 + h_y^2) == 0.
            return (_Vector([d, 0.0, 0.0]), _Vector([0.0, v_circular, 0.0]))

    class _Difference:
        def at(self, t):
            return _Position()

    class _Body:
        def __sub__(self, other):
            return _Difference()

    fake_planets = {"earth": _Body(), "moon": _Body()}

    from libephemeris import cache as cache_module

    monkeypatch.setattr(lunar, "get_planets", lambda: fake_planets)
    monkeypatch.setattr(cache_module, "get_cached_time_tt", lambda jd: object())

    lon, lat, dist = lunar.calc_true_lunar_node(2451545.0)

    assert lat == 0.0
    assert dist == pytest.approx(d)


# ---------------------------------------------------------------------------
# Cubic spline singular-pivot guard (line 2735)
# ---------------------------------------------------------------------------


def test_cubic_spline_singular_diagonal_guard():
    """Line 2735: a degenerate node set forces a near-zero diagonal pivot.

    Passing x with x[2] == x[0] makes diag[1] == 0, exercising the
    division-by-zero guard. The result must still be a finite number.
    """
    x = [1.0, 2.0, 1.0, 3.0]
    y = [0.0, 1.0, 0.0, 1.0]
    result = lunar._cubic_spline_interpolate(x, y, 1.5)
    assert math.isfinite(result)


def test_cubic_spline_trivial_returns_first():
    """Sanity: too few points returns y[0]."""
    assert lunar._cubic_spline_interpolate([5.0], [7.0], 5.0) == 7.0
    assert lunar._cubic_spline_interpolate([], [], 0.0) == 0.0


def test_cubic_spline_evaluation_branches():
    """Sanity: left-clamp, right-clamp, and interior evaluation."""
    x = [0.0, 1.0, 2.0, 3.0]
    y = [0.0, 1.0, 4.0, 9.0]
    assert math.isfinite(lunar._cubic_spline_interpolate(x, y, -1.0))
    assert math.isfinite(lunar._cubic_spline_interpolate(x, y, 5.0))
    assert math.isfinite(lunar._cubic_spline_interpolate(x, y, 1.5))


def test_cubic_spline_interval_search_exhausts(monkeypatch):
    """Arc 2768->2772: the interval search loop exits without a break.

    A NaN evaluation point passes the ``x[0] < x_eval < x[n]`` clamp checks
    (all NaN comparisons are False, so neither clamp branch is taken) and
    then matches no interval, so the binary-search loop runs to exhaustion
    and falls through to the post-loop ``dx`` assignment with ``idx == 0``.
    """
    x = [0.0, 1.0, 2.0, 3.0]
    y = [0.0, 1.0, 4.0, 9.0]
    result = lunar._cubic_spline_interpolate(x, y, float("nan"))
    assert math.isnan(result)


# ---------------------------------------------------------------------------
# Linear-system solver (lines 3046, 3054)
# ---------------------------------------------------------------------------


def test_solve_linear_system_pivot_swap():
    """Line 3046: a larger pivot in a later row triggers a row swap."""
    A = [[0.0, 1.0], [2.0, 1.0]]
    b = [1.0, 3.0]
    x = lunar._solve_linear_system(A, b)
    assert x == pytest.approx([1.0, 1.0])


def test_solve_linear_system_singular_returns_zeros():
    """Line 3054: a singular matrix returns the zero vector."""
    A = [[0.0, 0.0], [0.0, 0.0]]
    b = [1.0, 2.0]
    assert lunar._solve_linear_system(A, b) == [0.0, 0.0]


# ---------------------------------------------------------------------------
# Interpolated apside correction-skip branches (arcs 2636->2646, 2999->3009)
# ---------------------------------------------------------------------------


def test_interpolated_apogee_skips_correction_when_table_empty(monkeypatch):
    """Arc 2636->2646: empty apogee correction table skips the residual block."""
    monkeypatch.setattr(lunar, "APOGEE_CORRECTIONS", ())
    lon, lat, dist = lunar.calc_interpolated_apogee(2451545.0)
    assert 0.0 <= lon < 360.0
    assert dist == pytest.approx(0.0027099)


def test_interpolated_perigee_skips_correction_when_table_empty(monkeypatch):
    """Arc 2999->3009: empty perigee correction table skips the residual block."""
    monkeypatch.setattr(lunar, "PERIGEE_CORRECTIONS", ())
    lon, lat, dist = lunar.calc_interpolated_perigee(2451545.0)
    assert 0.0 <= lon < 360.0
    assert dist == pytest.approx(0.0024222)


# ---------------------------------------------------------------------------
# compare_true_lilith_methods longitude wrap branches (lines 3099-3114)
# ---------------------------------------------------------------------------


def test_compare_true_lilith_methods_identical():
    """Lines 3099-3105 + 3111-3114: identical methods give zero differences."""
    result = lunar.compare_true_lilith_methods(2451545.0)
    assert result["lon_diff"] == pytest.approx(0.0)
    assert result["lat_diff"] == pytest.approx(0.0)
    assert result["e_diff"] == pytest.approx(0.0)


def test_compare_true_lilith_methods_wrap_positive(monkeypatch):
    """Lines 3106-3107: a raw difference > 180 deg wraps down by 360.

    The orbital-elements longitude is offset so the raw difference
    (ev_lon - oe_lon) is exactly +200, forcing the ``> 180`` branch.
    """
    base = lunar.calc_true_lilith(2451545.0)
    monkeypatch.setattr(
        lunar,
        "calc_true_lilith_orbital_elements",
        lambda jd: (base[0] - 200.0, base[1], base[2]),
    )
    result = lunar.compare_true_lilith_methods(2451545.0)
    assert result["lon_diff"] == pytest.approx(-160.0)


def test_compare_true_lilith_methods_wrap_negative(monkeypatch):
    """Lines 3108-3109: a raw difference < -180 deg wraps up by 360.

    The orbital-elements longitude is offset so the raw difference
    (ev_lon - oe_lon) is exactly -200, forcing the ``< -180`` branch.
    """
    base = lunar.calc_true_lilith(2451545.0)
    monkeypatch.setattr(
        lunar,
        "calc_true_lilith_orbital_elements",
        lambda jd: (base[0] + 200.0, base[1], base[2]),
    )
    result = lunar.compare_true_lilith_methods(2451545.0)
    assert result["lon_diff"] == pytest.approx(160.0)
