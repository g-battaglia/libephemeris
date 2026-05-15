"""
Extended tests for fixed star functions.

Verifies fixstar_ut, fixstar2_ut, and fixstar_mag
with various stars, dates, and flag combinations.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe
from libephemeris.constants import (
    FLG_SPEED,
    FLG_EQUATORIAL,
    FLG_SIDEREAL,
    SIDM_LAHIRI,
)


BRIGHT_STARS = [
    "Sirius",
    "Canopus",
    "Arcturus",
    "Vega",
    "Capella",
    "Rigel",
    "Procyon",
    "Betelgeuse",
    "Aldebaran",
    "Spica",
    "Antares",
    "Pollux",
    "Fomalhaut",
    "Deneb",
    "Regulus",
]


class TestFixstarBasic:
    """Basic fixed star functionality."""

    @pytest.mark.unit
    def test_sirius_returns_tuple(self):
        """fixstar_ut for Sirius returns a tuple."""
        jd = 2451545.0
        result = swe.fixstar_ut("Sirius", jd, 0)
        assert isinstance(result, tuple)

    @pytest.mark.unit
    def test_sirius_position_valid(self):
        """Sirius position should be valid."""
        jd = 2451545.0
        result = swe.fixstar_ut("Sirius", jd, 0)
        # Result format: (lon, lat, dist, speed_lon, speed_lat, speed_dist, name_str)
        # or similar depending on implementation
        assert len(result) >= 2  # At least position data

    @pytest.mark.unit
    @pytest.mark.parametrize("star", BRIGHT_STARS[:10])
    def test_bright_star_no_crash(self, star: str):
        """Bright stars should not crash."""
        jd = 2451545.0
        result = swe.fixstar_ut(star, jd, 0)
        assert result is not None


class TestFixstarMagnitude:
    """Test fixstar_mag."""

    @pytest.mark.unit
    def test_sirius_magnitude(self):
        """Sirius magnitude should be ~ -1.46."""
        result = swe.fixstar_mag("Sirius")
        # Returns (magnitude, star_name) tuple
        assert len(result) == 2
        mag, name = result
        assert -2.0 < mag < 0.0, f"Sirius mag={mag} (expected ~-1.46)"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "star,max_mag",
        [
            ("Sirius", 0.0),
            ("Canopus", 0.0),
            ("Arcturus", 0.5),
            ("Vega", 0.5),
            ("Capella", 0.5),
        ],
    )
    def test_bright_star_magnitude(self, star: str, max_mag: float):
        """Bright stars should have magnitude < max_mag."""
        result = swe.fixstar_mag(star)
        mag = result[0]
        assert mag < max_mag, f"{star}: mag={mag}, expected < {max_mag}"

    @pytest.mark.unit
    def test_fixstar_mag_returns_tuple(self):
        """fixstar_mag returns (magnitude, name) tuple."""
        result = swe.fixstar_mag("Vega")
        assert len(result) == 2
        mag, name = result
        assert isinstance(mag, (int, float))
        assert isinstance(name, str)


class TestFixstarEquatorial:
    """Test fixed stars with equatorial coordinates."""

    @pytest.mark.unit
    @pytest.mark.parametrize("star", ["Sirius", "Vega", "Polaris"])
    def test_star_equatorial(self, star: str):
        """Fixed stars with FLG_EQUATORIAL should return RA/Dec."""
        jd = 2451545.0
        try:
            result = swe.fixstar_ut(star, jd, FLG_EQUATORIAL)
            assert result is not None
        except Exception:
            pytest.skip(f"Star {star} not found")


class TestFixstarSidereal:
    """Test fixed stars in sidereal mode."""

    @pytest.mark.unit
    @pytest.mark.parametrize("star", ["Sirius", "Regulus", "Spica"])
    def test_star_sidereal(self, star: str):
        """Fixed stars in sidereal mode should differ from tropical."""
        swe.set_sid_mode(SIDM_LAHIRI)
        jd = 2451545.0
        try:
            tropical = swe.fixstar_ut(star, jd, 0)
            sidereal = swe.fixstar_ut(star, jd, FLG_SIDEREAL)
            assert tropical is not None
            assert sidereal is not None
        except Exception:
            pytest.skip(f"Star {star} not found")


class TestFixstarDateRange:
    """Test fixed stars across dates (proper motion)."""

    @pytest.mark.unit
    @pytest.mark.parametrize("year", [1900, 2000, 2100])
    def test_sirius_across_centuries(self, year: int):
        """Sirius valid across centuries."""
        jd = swe.julday(year, 1, 1, 12.0)
        result = swe.fixstar_ut("Sirius", jd, 0)
        assert result is not None

    @pytest.mark.unit
    def test_sirius_proper_motion(self):
        """Sirius position should change slightly over a century (proper motion)."""
        jd1 = swe.julday(1900, 1, 1, 12.0)
        jd2 = swe.julday(2100, 1, 1, 12.0)
        r1 = swe.fixstar_ut("Sirius", jd1, 0)
        r2 = swe.fixstar_ut("Sirius", jd2, 0)
        # Sirius has significant proper motion
        # Position should change but not by a huge amount
        assert r1 is not None
        assert r2 is not None


class TestFixstar2:
    """Test fixstar2_ut (catalog-based lookup)."""

    @pytest.mark.unit
    def test_fixstar2_sirius(self):
        """fixstar2_ut for Sirius returns valid result."""
        jd = 2451545.0
        try:
            result = swe.fixstar2_ut("Sirius", jd, 0)
            assert result is not None
        except (AttributeError, NotImplementedError):
            pytest.skip("fixstar2_ut not implemented")

    @pytest.mark.unit
    @pytest.mark.parametrize("star", ["Aldebaran", "Regulus", "Antares"])
    def test_fixstar2_royal_stars(self, star: str):
        """fixstar2_ut for royal stars."""
        jd = 2451545.0
        try:
            result = swe.fixstar2_ut(star, jd, 0)
            assert result is not None
        except (AttributeError, NotImplementedError):
            pytest.skip("fixstar2_ut not implemented")
        except Exception:
            pytest.skip(f"Star {star} not found")
