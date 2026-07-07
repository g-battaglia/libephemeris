"""Regression tests for the v3 crossing fixes.

Covers two independent bugs:

* ``helio_cross_ut`` / ``helio_cross`` backward search jumping a whole
  revolution too far back near a planet's aphelion, where the instantaneous
  speed runs well below the mean (Mercury ~2.75 vs ~4.09 deg/day) and Newton
  converged into the wrong basin.
* ``cross_ut`` diverging on valid but long-period bodies outside the
  typical-speed table (Mean Node, ~18.6 yr and perpetually retrograde;
  Mean Apogee, ~8.85 yr), whose genuine crossings lie thousands of days out,
  well beyond the generic ~800 day horizon.
"""

from __future__ import annotations

import random

import pytest

from libephemeris.constants import (
    FLG_HELCTR,
    FLG_SPEED,
    FLG_SWIEPH,
    MEAN_APOG,
    MEAN_NODE,
    MERCURY,
)
from libephemeris.crossing import cross_ut, helio_cross, helio_cross_ut
from libephemeris.planets import calc, calc_ut
from libephemeris.time_utils import julday

# Heliocentric orbital period of Mercury (days). The nearest strictly-past
# crossing of a fixed longitude is always within one period, so an answer with
# |offset| > this signals a revolution was skipped.
MERCURY_HELIO_PERIOD = 87.97


def _wrap(a: float, b: float) -> float:
    """Signed smallest angular difference a - b, wrapped to (-180, 180]."""
    return (a - b + 180.0) % 360.0 - 180.0


class TestHelioCrossBackwardNoSkippedRevolution:
    """C1: backward helio search must return the nearest strictly-past crossing."""

    @pytest.mark.unit
    def test_mercury_backward_exact_repro(self):
        # Exact reported case: the buggy result was -175.79 days (two Mercury
        # revolutions back); the correct nearest-past crossing is -87.82 days.
        start = 2451980.089
        result = helio_cross_ut(MERCURY, 241.0510, start, 0, True)
        offset = result - start
        assert offset == pytest.approx(-87.82, abs=1.0), (
            f"expected ~-87.82 d, got {offset:.4f} d "
            "(a skipped revolution lands near -175.79 d)"
        )
        # Longitude at the crossing equals the target.
        lon = calc_ut(result, MERCURY, FLG_SWIEPH | FLG_HELCTR | FLG_SPEED)[0][0]
        assert abs(_wrap(lon, 241.0510)) * 3600.0 < 1.0

    @pytest.mark.unit
    def test_mercury_backward_tt_twin_exact_repro(self):
        start = 2451980.089
        result = helio_cross(MERCURY, 241.0510, start, 0, True)
        offset = result - start
        assert offset == pytest.approx(-87.82, abs=1.0)
        lon = calc(result, MERCURY, FLG_SWIEPH | FLG_HELCTR | FLG_SPEED)[0][0]
        assert abs(_wrap(lon, 241.0510)) * 3600.0 < 1.0

    @pytest.mark.unit
    def test_mercury_backward_fuzz_no_skipped_revolution(self):
        # 20 seeded trials: every backward crossing must be strictly in the
        # past, land on the target longitude, and lie within one orbital period
        # of the start (the invariant a skipped revolution violates).
        rng = random.Random(20240707)
        base = julday(2000, 1, 1)
        for _ in range(20):
            start = base + rng.uniform(-3000.0, 3000.0)
            target = rng.uniform(0.0, 360.0)
            result = helio_cross_ut(MERCURY, target, start, 0, True)
            offset = result - start
            assert result <= start + 1e-6, f"not in the past: offset={offset:.3f}"
            assert abs(offset) <= MERCURY_HELIO_PERIOD + 2.0, (
                f"skipped a revolution: target={target:.2f}, offset={offset:.3f}"
            )
            lon = calc_ut(result, MERCURY, FLG_SWIEPH | FLG_HELCTR | FLG_SPEED)[0][0]
            assert abs(_wrap(lon, target)) * 3600.0 < 1.0

    @pytest.mark.unit
    def test_mercury_backward_tt_fuzz_no_skipped_revolution(self):
        rng = random.Random(999001)
        base = julday(2000, 1, 1)
        for _ in range(20):
            start = base + rng.uniform(-3000.0, 3000.0)
            target = rng.uniform(0.0, 360.0)
            result = helio_cross(MERCURY, target, start, 0, True)
            offset = result - start
            assert result <= start + 1e-6
            assert abs(offset) <= MERCURY_HELIO_PERIOD + 2.0
            lon = calc(result, MERCURY, FLG_SWIEPH | FLG_HELCTR | FLG_SPEED)[0][0]
            assert abs(_wrap(lon, target)) * 3600.0 < 1.0


class TestCrossLongPeriodBodies:
    """C3: cross_ut must resolve crossings of slow off-table bodies."""

    @pytest.mark.unit
    def test_mean_node_far_crossing(self):
        # Mean Node is perpetually retrograde (~18.6 yr period). Reaching 123 deg
        # from 2024-01-01 takes ~4869 days; the old ~800 day horizon reported
        # divergence.
        start = julday(2024, 1, 1)
        result = cross_ut(MEAN_NODE, 123.0, start)
        assert result >= start - 1e-6
        assert (result - start) == pytest.approx(4869.30, abs=2.0)
        lon = calc_ut(result, MEAN_NODE, FLG_SWIEPH)[0][0]
        assert abs(_wrap(lon, 123.0)) * 3600.0 < 1.0

    @pytest.mark.unit
    def test_mean_apogee_far_crossing(self):
        # Mean Apogee (~8.85 yr period, prograde) reaching 300 deg from
        # 2024-01-01 takes ~1257 days, again beyond the generic horizon.
        start = julday(2024, 1, 1)
        result = cross_ut(MEAN_APOG, 300.0, start)
        assert result >= start - 1e-6
        assert (result - start) == pytest.approx(1256.92, abs=2.0)
        lon = calc_ut(result, MEAN_APOG, FLG_SWIEPH)[0][0]
        assert abs(_wrap(lon, 300.0)) * 3600.0 < 1.0

    @pytest.mark.unit
    def test_mean_node_near_crossing_unchanged(self):
        # A nearby target (~5 deg retrograde ahead, reached in ~94 days) already
        # worked before the fix and must be unchanged.
        start = julday(2024, 1, 1)
        node_lon = calc_ut(start, MEAN_NODE, FLG_SWIEPH)[0][0]
        target = (node_lon - 5.0) % 360.0
        result = cross_ut(MEAN_NODE, target, start)
        assert result >= start - 1e-6
        assert (result - start) == pytest.approx(94.42, abs=2.0)
        lon = calc_ut(result, MEAN_NODE, FLG_SWIEPH)[0][0]
        assert abs(_wrap(lon, target)) * 3600.0 < 1.0
