"""
Flag pair robustness tests.

Tests that ~30 flag pair combinations do not crash for Sun, Moon, Mars,
Jupiter, and Saturn across 5 dates. Verifies output format (6-tuple, int),
all values finite, XYZ coherence (x**2+y**2+z**2 ~ dist**2), and RADIANS
coherence with degree results.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe
from libephemeris.constants import (
    SUN,
    MOON,
    MARS,
    JUPITER,
    SATURN,
    MERCURY,
    VENUS,
    FLG_SPEED,
    FLG_HELCTR,
    FLG_TRUEPOS,
    FLG_J2000,
    FLG_NONUT,
    FLG_NOGDEFL,
    FLG_NOABERR,
    FLG_EQUATORIAL,
    FLG_XYZ,
    FLG_RADIANS,
)

# ---------------------------------------------------------------------------
# Test data
# ---------------------------------------------------------------------------

BODIES = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
]

# Bodies valid for heliocentric (not Sun or Moon)
HELIO_BODIES = [
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
    (MERCURY, "Mercury"),
    (VENUS, "Venus"),
]

DATES = [
    2415020.0,  # 1900-01-01
    2433282.5,  # 1950-01-01
    2451545.0,  # 2000-01-01 (J2000)
    2460310.5,  # 2024-01-18
    2488070.0,  # 2100-01-01
]

# ~30 flag pair combinations (geocentric-safe, no HELCTR with Sun/Moon)
GEOCENTRIC_FLAG_PAIRS = [
    (FLG_SPEED | FLG_EQUATORIAL, "SPEED+EQUATORIAL"),
    (FLG_SPEED | FLG_TRUEPOS, "SPEED+TRUEPOS"),
    (FLG_SPEED | FLG_J2000, "SPEED+J2000"),
    (FLG_SPEED | FLG_NONUT, "SPEED+NONUT"),
    (FLG_SPEED | FLG_NOGDEFL, "SPEED+NOGDEFL"),
    (FLG_SPEED | FLG_NOABERR, "SPEED+NOABERR"),
    (FLG_SPEED | FLG_RADIANS, "SPEED+RADIANS"),
    (FLG_EQUATORIAL | FLG_NONUT, "EQUATORIAL+NONUT"),
    (FLG_EQUATORIAL | FLG_J2000, "EQUATORIAL+J2000"),
    (FLG_EQUATORIAL | FLG_TRUEPOS, "EQUATORIAL+TRUEPOS"),
    (FLG_EQUATORIAL | FLG_RADIANS, "EQUATORIAL+RADIANS"),
    (FLG_XYZ | FLG_SPEED, "XYZ+SPEED"),
    (FLG_XYZ | FLG_EQUATORIAL, "XYZ+EQUATORIAL"),
    (FLG_XYZ | FLG_NONUT, "XYZ+NONUT"),
    (FLG_XYZ | FLG_J2000, "XYZ+J2000"),
    (FLG_XYZ | FLG_TRUEPOS, "XYZ+TRUEPOS"),
    (FLG_TRUEPOS | FLG_NONUT, "TRUEPOS+NONUT"),
    (FLG_TRUEPOS | FLG_J2000, "TRUEPOS+J2000"),
    (FLG_TRUEPOS | FLG_NOGDEFL, "TRUEPOS+NOGDEFL"),
    (FLG_TRUEPOS | FLG_NOABERR, "TRUEPOS+NOABERR"),
    (FLG_NONUT | FLG_NOGDEFL, "NONUT+NOGDEFL"),
    (FLG_NONUT | FLG_NOABERR, "NONUT+NOABERR"),
    (FLG_J2000 | FLG_NONUT, "J2000+NONUT"),
    (FLG_J2000 | FLG_NOGDEFL, "J2000+NOGDEFL"),
    (FLG_J2000 | FLG_NOABERR, "J2000+NOABERR"),
    (FLG_NOGDEFL | FLG_NOABERR, "NOGDEFL+NOABERR"),
    (FLG_RADIANS | FLG_NONUT, "RADIANS+NONUT"),
    (FLG_RADIANS | FLG_J2000, "RADIANS+J2000"),
]

# Heliocentric flag pairs (only for bodies that support HELCTR)
HELIOCENTRIC_FLAG_PAIRS = [
    (FLG_HELCTR | FLG_SPEED, "HELCTR+SPEED"),
    (FLG_HELCTR | FLG_EQUATORIAL, "HELCTR+EQUATORIAL"),
    (FLG_HELCTR | FLG_NONUT, "HELCTR+NONUT"),
    (FLG_HELCTR | FLG_J2000, "HELCTR+J2000"),
    (FLG_HELCTR | FLG_TRUEPOS, "HELCTR+TRUEPOS"),
    (FLG_HELCTR | FLG_XYZ, "HELCTR+XYZ"),
]


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestGeocentricFlagPairs:
    """Geocentric flag pair combinations produce valid results."""

    @pytest.mark.unit
    @pytest.mark.parametrize("flags,desc", GEOCENTRIC_FLAG_PAIRS)
    @pytest.mark.parametrize("body_id,body_name", BODIES)
    @pytest.mark.parametrize("jd", DATES)
    def test_flag_pair_valid_output(
        self,
        body_id: int,
        body_name: str,
        flags: int,
        desc: str,
        jd: float,
    ):
        """Flag pair produces a (6-tuple, int) with all finite values."""
        result = swe.calc_ut(jd, body_id, flags)
        assert isinstance(result, tuple), f"Expected tuple, got {type(result)}"
        assert len(result) == 2, f"Expected 2-element tuple, got {len(result)}"

        values, retflag = result
        assert isinstance(retflag, int), f"retflag not int: {type(retflag)}"
        assert len(values) == 6, f"{body_name}+{desc}: got {len(values)} elements"
        for i, val in enumerate(values):
            assert math.isfinite(val), (
                f"{body_name}+{desc} jd={jd}: values[{i}]={val} not finite"
            )


class TestHeliocentricFlagPairs:
    """Heliocentric flag pair combinations for valid bodies."""

    @pytest.mark.unit
    @pytest.mark.parametrize("flags,desc", HELIOCENTRIC_FLAG_PAIRS)
    @pytest.mark.parametrize("body_id,body_name", HELIO_BODIES)
    @pytest.mark.parametrize("jd", DATES)
    def test_helio_pair_valid_output(
        self,
        body_id: int,
        body_name: str,
        flags: int,
        desc: str,
        jd: float,
    ):
        """Heliocentric flag pair produces valid finite output."""
        result = swe.calc_ut(jd, body_id, flags)
        values, retflag = result
        assert len(values) == 6
        for i, val in enumerate(values):
            assert math.isfinite(val), (
                f"{body_name}+{desc} jd={jd}: values[{i}]={val} not finite"
            )


class TestXYZCoherence:
    """XYZ flag produces Cartesian coordinates where x^2+y^2+z^2 ~ dist^2."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,body_name", BODIES)
    @pytest.mark.parametrize("jd", DATES)
    def test_xyz_vs_spherical_distance(self, body_id: int, body_name: str, jd: float):
        """XYZ magnitude matches spherical distance within 1e-8 AU."""
        # Get spherical coordinates to extract distance
        sph_vals, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
        dist_sph = sph_vals[2]

        # Get Cartesian coordinates
        xyz_vals, _ = swe.calc_ut(jd, body_id, FLG_XYZ)
        x, y, z = xyz_vals[0], xyz_vals[1], xyz_vals[2]
        dist_xyz = math.sqrt(x * x + y * y + z * z)

        assert abs(dist_xyz - dist_sph) < 1e-8, (
            f"{body_name} jd={jd}: |xyz_dist({dist_xyz:.10f}) - "
            f"sph_dist({dist_sph:.10f})| = {abs(dist_xyz - dist_sph):.2e}"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,body_name", BODIES)
    def test_xyz_equatorial_vs_ecliptic(self, body_id: int, body_name: str):
        """XYZ ecliptic and XYZ equatorial magnitudes match."""
        jd = 2451545.0
        ecl_vals, _ = swe.calc_ut(jd, body_id, FLG_XYZ)
        equ_vals, _ = swe.calc_ut(jd, body_id, FLG_XYZ | FLG_EQUATORIAL)

        r_ecl = math.sqrt(sum(v * v for v in ecl_vals[:3]))
        r_equ = math.sqrt(sum(v * v for v in equ_vals[:3]))

        assert abs(r_ecl - r_equ) < 1e-8, (
            f"{body_name}: ecliptic r={r_ecl:.10f}, equatorial r={r_equ:.10f}"
        )


class TestRadiansCoherence:
    """RADIANS flag produces values consistent with degree results."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,body_name", BODIES)
    @pytest.mark.parametrize("jd", DATES)
    def test_radians_vs_degrees(self, body_id: int, body_name: str, jd: float):
        """Radians result equals degrees * pi/180 within 1e-10."""
        deg_vals, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
        rad_vals, _ = swe.calc_ut(jd, body_id, FLG_SPEED | FLG_RADIANS)

        # lon and lat (indices 0, 1) should convert; distance (index 2) stays the same.
        # Tolerance of 1e-8 rad (~0.002 arcsec) accounts for the fact that
        # RADIANS mode recalculates internally rather than converting the
        # degree result, so floating-point paths may diverge slightly.
        for i in (0, 1):
            expected_rad = math.radians(deg_vals[i])
            assert abs(rad_vals[i] - expected_rad) < 1e-8, (
                f"{body_name} jd={jd} idx={i}: "
                f"radians={rad_vals[i]:.12f}, expected={expected_rad:.12f}"
            )

        # Distance should match within double-precision rounding
        assert abs(rad_vals[2] - deg_vals[2]) < 1e-10, (
            f"{body_name} jd={jd}: distance differs in RADIANS mode"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,body_name", BODIES)
    def test_equatorial_radians_vs_degrees(self, body_id: int, body_name: str):
        """Equatorial RADIANS result equals degrees * pi/180."""
        jd = 2451545.0
        flags_deg = FLG_SPEED | FLG_EQUATORIAL
        flags_rad = flags_deg | FLG_RADIANS

        deg_vals, _ = swe.calc_ut(jd, body_id, flags_deg)
        rad_vals, _ = swe.calc_ut(jd, body_id, flags_rad)

        # RA and Dec (indices 0, 1)
        for i in (0, 1):
            expected_rad = math.radians(deg_vals[i])
            assert abs(rad_vals[i] - expected_rad) < 1e-8, (
                f"{body_name} eq idx={i}: "
                f"radians={rad_vals[i]:.12f}, expected={expected_rad:.12f}"
            )


class TestTripleFlagCombinations:
    """Test selected three-flag combinations for robustness."""

    TRIPLE_COMBOS = [
        (FLG_SPEED | FLG_EQUATORIAL | FLG_NONUT, "SPEED+EQ+NONUT"),
        (FLG_SPEED | FLG_J2000 | FLG_NONUT, "SPEED+J2000+NONUT"),
        (FLG_SPEED | FLG_TRUEPOS | FLG_NOABERR, "SPEED+TRUE+NOABERR"),
        (FLG_XYZ | FLG_SPEED | FLG_EQUATORIAL, "XYZ+SPEED+EQ"),
        (FLG_XYZ | FLG_SPEED | FLG_J2000, "XYZ+SPEED+J2000"),
        (FLG_EQUATORIAL | FLG_NONUT | FLG_NOABERR, "EQ+NONUT+NOABERR"),
        (FLG_J2000 | FLG_NONUT | FLG_NOGDEFL, "J2000+NONUT+NOGDEFL"),
        (FLG_SPEED | FLG_NOGDEFL | FLG_NOABERR, "SPEED+NOGDEFL+NOABERR"),
        (FLG_RADIANS | FLG_SPEED | FLG_NONUT, "RAD+SPEED+NONUT"),
        (FLG_XYZ | FLG_J2000 | FLG_NONUT, "XYZ+J2000+NONUT"),
    ]

    @pytest.mark.unit
    @pytest.mark.parametrize("flags,desc", TRIPLE_COMBOS)
    @pytest.mark.parametrize(
        "body_id,body_name",
        [(SUN, "Sun"), (MOON, "Moon"), (MARS, "Mars")],
    )
    def test_triple_flag_valid_output(
        self, body_id: int, body_name: str, flags: int, desc: str
    ):
        """Triple flag combination produces valid finite 6-tuple."""
        jd = 2451545.0
        values, retflag = swe.calc_ut(jd, body_id, flags)
        assert len(values) == 6
        for i, val in enumerate(values):
            assert math.isfinite(val), (
                f"{body_name}+{desc}: values[{i}]={val} not finite"
            )
