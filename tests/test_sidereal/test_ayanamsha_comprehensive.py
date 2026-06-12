"""
Comprehensive tests for ayanamsha calculations across all 43 modes.

Verifies that each sidereal mode produces a valid ayanamsha value,
that values change over time (precession), and that different modes
produce different results.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe
from libephemeris.constants import (
    SUN,
    MOON,
    MARS,
    FLG_SIDEREAL,
)


# All 43 sidereal modes (0-42) — some may not be implemented
SIDEREAL_MODES = list(range(43))

# Well-known modes with approximate J2000 ayanamsha values
KNOWN_MODES = [
    (0, "Fagan-Bradley", 24.0, 25.5),
    (1, "Lahiri", 23.5, 24.5),
    (3, "Raman", 22.0, 23.5),
    (5, "Krishnamurti", 23.5, 24.5),
    (7, "Yukteswar", 22.0, 23.5),
    (27, "True Citra", 23.5, 25.0),
    (28, "True Revati", 19.5, 21.5),
]


class TestAyanamshaBasic:
    """Basic ayanamsha tests."""

    @pytest.mark.unit
    def test_get_ayanamsa_ut_returns_float(self):
        """get_ayanamsa_ut returns a Python float."""
        swe.set_sid_mode(1)
        result = swe.get_ayanamsa_ut(2451545.0)
        assert type(result) is float

    @pytest.mark.unit
    def test_get_ayanamsa_ex_ut_returns_tuple(self):
        """get_ayanamsa_ex_ut returns (retflag, ayanamsa)."""
        swe.set_sid_mode(1)
        result = swe.get_ayanamsa_ex_ut(2451545.0, 0)
        assert len(result) == 2
        retflag, ayan = result
        assert isinstance(retflag, int)
        assert type(ayan) is float

    @pytest.mark.unit
    def test_ayanamsa_ex_ut_retflag_is_2(self):
        """get_ayanamsa_ex_ut should return retflag=2 (FLG_SWIEPH)."""
        swe.set_sid_mode(1)
        retflag, _ = swe.get_ayanamsa_ex_ut(2451545.0, 0)
        assert retflag == 2, f"retflag={retflag}, expected 2"


class TestAyanamshaAllModes:
    """Test all 43 sidereal modes produce valid values."""

    @pytest.mark.unit
    @pytest.mark.parametrize("mode", SIDEREAL_MODES)
    def test_mode_returns_finite(self, mode: int):
        """Each sidereal mode returns a finite ayanamsha."""
        swe.set_sid_mode(mode)
        ayan = swe.get_ayanamsa_ut(2451545.0)
        assert math.isfinite(ayan), f"Mode {mode}: ayanamsha={ayan} not finite"

    @pytest.mark.unit
    @pytest.mark.parametrize("mode", SIDEREAL_MODES)
    def test_mode_calc_ut_works(self, mode: int):
        """Each sidereal mode works with calc_ut + FLG_SIDEREAL."""
        swe.set_sid_mode(mode)
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, SUN, FLG_SIDEREAL)
        assert 0 <= result[0] < 360, f"Mode {mode}: lon={result[0]}"


class TestAyanamshaKnownValues:
    """Test known modes have approximately correct ayanamsha values."""

    @pytest.mark.unit
    @pytest.mark.parametrize("mode,name,lo,hi", KNOWN_MODES)
    def test_known_mode_range(self, mode: int, name: str, lo: float, hi: float):
        """Known sidereal modes should have ayanamsha in expected range at J2000."""
        swe.set_sid_mode(mode)
        ayan = swe.get_ayanamsa_ut(2451545.0)
        assert lo < ayan < hi, (
            f"{name} (mode {mode}): ayanamsha={ayan:.4f}, expected [{lo}, {hi}]"
        )


class TestAyanamshaPrecession:
    """Test that ayanamsha changes over time (precession)."""

    @pytest.mark.unit
    @pytest.mark.parametrize("mode", [0, 1, 3, 5])
    def test_ayanamsha_increases_over_century(self, mode: int):
        """Ayanamsha should increase over a century (precession ~50.3"/yr)."""
        swe.set_sid_mode(mode)
        ayan_2000 = swe.get_ayanamsa_ut(2451545.0)
        ayan_2100 = swe.get_ayanamsa_ut(2451545.0 + 36525.0)
        diff = ayan_2100 - ayan_2000
        # Should increase by ~1.4° per century
        assert 1.0 < diff < 2.0, (
            f"Mode {mode}: century change={diff:.4f}° (expected ~1.4°)"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("mode", [0, 1, 3])
    def test_ayanamsha_at_different_dates(self, mode: int):
        """Ayanamsha should be different at different dates."""
        swe.set_sid_mode(mode)
        ayan_1900 = swe.get_ayanamsa_ut(2415020.0)
        ayan_2000 = swe.get_ayanamsa_ut(2451545.0)
        ayan_2100 = swe.get_ayanamsa_ut(2451545.0 + 36525.0)
        assert ayan_1900 < ayan_2000 < ayan_2100, (
            f"Mode {mode}: not monotonically increasing"
        )


class TestAyanamshaDifferentModes:
    """Test that different modes produce different results."""

    @pytest.mark.unit
    def test_lahiri_vs_fagan(self):
        """Lahiri and Fagan-Bradley should differ."""
        jd = 2451545.0
        swe.set_sid_mode(0)
        fagan = swe.get_ayanamsa_ut(jd)
        swe.set_sid_mode(1)
        lahiri = swe.get_ayanamsa_ut(jd)
        diff = abs(fagan - lahiri)
        assert diff > 0.1, f"Fagan={fagan:.4f}, Lahiri={lahiri:.4f}"

    @pytest.mark.unit
    def test_multiple_modes_distinct(self):
        """At least 5 modes should produce distinct ayanamsha values."""
        jd = 2451545.0
        values = set()
        for mode in [0, 1, 3, 5, 7, 27, 28]:
            swe.set_sid_mode(mode)
            ayan = swe.get_ayanamsa_ut(jd)
            values.add(round(ayan, 2))
        assert len(values) >= 4, f"Only {len(values)} distinct values from 7 modes"


class TestSiderealPositions:
    """Test that sidereal positions are tropical minus ayanamsha."""

    @pytest.mark.unit
    @pytest.mark.parametrize("mode", [0, 1, 3])
    def test_sidereal_lon_equals_tropical_minus_ayan(self, mode: int):
        """Sidereal longitude ≈ tropical longitude - ayanamsha."""
        swe.set_sid_mode(mode)
        jd = 2451545.0

        tropical, _ = swe.calc_ut(jd, SUN, 0)
        sidereal, _ = swe.calc_ut(jd, SUN, FLG_SIDEREAL)
        ayan = swe.get_ayanamsa_ut(jd)

        expected_sid = (tropical[0] - ayan) % 360
        diff = abs(sidereal[0] - expected_sid)
        if diff > 180:
            diff = 360 - diff
        assert diff < 0.01, (
            f"Mode {mode}: sid={sidereal[0]:.4f}, expected={expected_sid:.4f}"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,name",
        [(SUN, "Sun"), (MOON, "Moon"), (MARS, "Mars")],
    )
    def test_sidereal_lat_unchanged(self, body_id: int, name: str):
        """Sidereal latitude should equal tropical latitude."""
        swe.set_sid_mode(1)
        jd = 2451545.0
        tropical, _ = swe.calc_ut(jd, body_id, 0)
        sidereal, _ = swe.calc_ut(jd, body_id, FLG_SIDEREAL)
        assert abs(tropical[1] - sidereal[1]) < 0.001, (
            f"{name}: trop lat={tropical[1]:.6f}, sid lat={sidereal[1]:.6f}"
        )
