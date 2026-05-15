"""
Delta T tests across centuries.

Verifies that libephemeris deltat() returns physically plausible
values across the full date range, with correct sign and magnitude.
"""

from __future__ import annotations

import math
import random

import pytest

import libephemeris as swe


def _jd_from_year(year: float) -> float:
    """Approximate JD for a given year (Jan 1 noon)."""
    return 2451545.0 + (year - 2000.0) * 365.25


class TestDeltaTBasic:
    """Basic Delta T functionality."""

    @pytest.mark.unit
    def test_deltat_returns_float(self):
        """deltat returns a float value."""
        jd = 2451545.0
        dt = swe.deltat(jd)
        assert isinstance(dt, float)
        assert math.isfinite(dt)

    @pytest.mark.unit
    def test_deltat_at_j2000(self):
        """Delta T at J2000 should be approximately 63.8 seconds."""
        jd = 2451545.0
        dt = swe.deltat(jd)
        dt_seconds = dt * 86400.0
        assert 60 < dt_seconds < 70, (
            f"Delta T at J2000 = {dt_seconds:.1f}s, expected ~63.8s"
        )

    @pytest.mark.unit
    def test_deltat_positive_modern(self):
        """Delta T should be positive for modern dates after 1902 (TT > UT)."""
        for year in [1950, 2000, 2024]:
            jd = _jd_from_year(year)
            dt = swe.deltat(jd)
            assert dt > 0, f"Delta T negative at year {year}: {dt}"


class TestDeltaTAcrossCenturies:
    """Test Delta T values across the full date range."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "year,min_sec,max_sec",
        [
            (1600, 100, 160),  # ~120s
            (1700, 5, 25),  # ~21s
            (1800, 10, 20),  # ~13s
            (1900, -5, 5),  # ~-2.7s
            (1950, 25, 35),  # ~29s
            (2000, 60, 70),  # ~63.8s
            (2020, 68, 72),  # ~69s
        ],
    )
    def test_deltat_approximate_value(self, year: int, min_sec: float, max_sec: float):
        """Delta T at known epochs should be within expected range."""
        jd = _jd_from_year(year)
        dt_seconds = swe.deltat(jd) * 86400.0
        assert min_sec <= dt_seconds <= max_sec, (
            f"Delta T at {year} = {dt_seconds:.1f}s, expected [{min_sec}, {max_sec}]"
        )

    @pytest.mark.unit
    def test_deltat_monotonic_increase_recent(self):
        """Delta T should generally increase from 1950 to present."""
        prev_dt = None
        for year in range(1960, 2025, 5):
            jd = _jd_from_year(year)
            dt = swe.deltat(jd)
            if prev_dt is not None:
                assert dt >= prev_dt - 1e-6, (
                    f"Delta T decreased from {year - 5} to {year}"
                )
            prev_dt = dt


class TestDeltaTFinite:
    """Test Delta T is finite across the full range."""

    @pytest.mark.unit
    @pytest.mark.parametrize("year", list(range(1550, 2650, 50)))
    def test_deltat_finite_at_century(self, year: int):
        """Delta T is finite at each 50-year interval."""
        jd = _jd_from_year(year)
        dt = swe.deltat(jd)
        assert math.isfinite(dt), f"Delta T not finite at year {year}"
        # Should be a reasonable number of days (< 1 day)
        assert abs(dt) < 1.0, f"Delta T = {dt} days at year {year} too large"

    @pytest.mark.unit
    @pytest.mark.parametrize("year", list(range(1600, 2500, 10)))
    def test_deltat_smooth(self, year: int):
        """Delta T should change smoothly (no huge jumps between decades)."""
        jd1 = _jd_from_year(year)
        jd2 = _jd_from_year(year + 10)
        dt1 = swe.deltat(jd1) * 86400.0
        dt2 = swe.deltat(jd2) * 86400.0
        diff = abs(dt2 - dt1)
        # Over 10 years, Delta T change grows with distance from present.
        # Far-future extrapolation can change ~45s/decade. Use 60s limit.
        assert diff < 60.0, f"Delta T jump {diff:.2f}s between {year} and {year + 10}"


class TestDeltaTHighVolume:
    """High-volume random date Delta T tests."""

    @pytest.mark.unit
    @pytest.mark.parametrize("idx", range(200))
    def test_deltat_200_random_dates(self, idx: int):
        """Delta T is finite at 200 random dates."""
        rng = random.Random(idx + 8888)
        year = rng.uniform(1550, 2650)
        jd = _jd_from_year(year)
        dt = swe.deltat(jd)
        assert math.isfinite(dt), f"Delta T not finite at year {year:.1f}"
        # Delta T can be up to ~2000s (~0.023 days) at far-future/past dates
        assert abs(dt) < 0.03, (
            f"Delta T = {dt:.6f} days ({dt * 86400:.1f}s) at year {year:.1f}"
        )
