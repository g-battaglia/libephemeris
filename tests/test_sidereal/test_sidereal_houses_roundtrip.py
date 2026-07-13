"""
Tests for sidereal houses round-trip and consistency.

Verifies that sidereal house cusps = tropical cusps - ayanamsha,
and that sidereal houses_ex is consistent across modes.
"""

from __future__ import annotations


import pytest

import libephemeris as swe
from libephemeris.constants import (
    FLG_SIDEREAL,
    SIDM_J1900,
    SIDM_J2000,
)


HOUSE_SYSTEMS_FOR_SIDEREAL = [
    (ord("P"), "Placidus"),
    (ord("K"), "Koch"),
    (ord("R"), "Regiomontanus"),
    (ord("E"), "Equal"),
    (ord("W"), "Whole Sign"),
    (ord("O"), "Porphyry"),
    (ord("C"), "Campanus"),
]


class TestSiderealHousesRoundTrip:
    """Test sidereal house cusp = tropical - ayanamsha."""

    @pytest.mark.unit
    @pytest.mark.parametrize("hsys,name", HOUSE_SYSTEMS_FOR_SIDEREAL)
    def test_cusps_differ_by_ayanamsha(self, hsys: int, name: str):
        """Sidereal cusps should differ from tropical by ~ayanamsha."""
        jd = 2451545.0
        lat, lon = 41.9, 12.5

        swe.set_sid_mode(SIDM_J2000)

        cusps_trop, _ = swe.houses_ex(jd, lat, lon, hsys, 0)
        cusps_sid, _ = swe.houses_ex(jd, lat, lon, hsys, FLG_SIDEREAL)

        ayan = swe.get_ayanamsa_ut(jd)

        # For most systems, cusp_sid ≈ cusp_trop - ayan (mod 360)
        # Whole Sign is special — recalculated from sidereal ASC
        diffs = []
        for i in range(12):
            expected_sid = (cusps_trop[i] - ayan) % 360
            actual_diff = abs(cusps_sid[i] - expected_sid)
            if actual_diff > 180:
                actual_diff = 360 - actual_diff
            diffs.append(actual_diff)

        avg_diff = sum(diffs) / len(diffs)
        # Most systems should have avg diff < 1 deg
        # Whole Sign may differ more (cusps snap to sign boundaries)
        if name == "Whole Sign":
            assert avg_diff < 30, f"{name}: avg cusp diff {avg_diff:.2f} deg"
        else:
            assert avg_diff < 1.0, f"{name}: avg cusp diff {avg_diff:.2f} deg"

    @pytest.mark.unit
    def test_sidereal_cusps_in_range(self):
        """All sidereal cusps should be [0, 360)."""
        jd = 2451545.0
        swe.set_sid_mode(SIDM_J2000)
        cusps, _ = swe.houses_ex(jd, 41.9, 12.5, ord("P"), FLG_SIDEREAL)
        for i, c in enumerate(cusps[:12]):
            assert 0 <= c < 360, f"Sidereal cusp {i + 1}: {c} deg"


class TestSiderealModesHouses:
    """Test independently defined frames produce different house cusps."""

    @pytest.mark.unit
    def test_j2000_vs_j1900_cusps_differ(self):
        """The J2000 and J1900 frames should produce different house cusps."""
        jd = 2451545.0

        swe.set_sid_mode(SIDM_J2000)
        cusps_j2000, _ = swe.houses_ex(jd, 41.9, 12.5, ord("P"), FLG_SIDEREAL)

        swe.set_sid_mode(SIDM_J1900)
        cusps_j1900, _ = swe.houses_ex(jd, 41.9, 12.5, ord("P"), FLG_SIDEREAL)

        diffs = []
        for i in range(12):
            diff = abs(cusps_j2000[i] - cusps_j1900[i])
            if diff > 180:
                diff = 360 - diff
            diffs.append(diff)

        avg_diff = sum(diffs) / len(diffs)
        assert avg_diff > 1.0, f"Avg diff: {avg_diff:.2f} deg"


class TestSiderealHousesAcrossDates:
    """Test sidereal houses across dates."""

    @pytest.mark.unit
    @pytest.mark.parametrize("year", [1900, 1950, 2000, 2050, 2100])
    def test_sidereal_placidus_across_years(self, year: int):
        """Sidereal Placidus valid across years."""
        jd = swe.julday(year, 6, 21, 12.0)
        swe.set_sid_mode(SIDM_J2000)
        cusps, ascmc = swe.houses_ex(jd, 41.9, 12.5, ord("P"), FLG_SIDEREAL)
        assert len(cusps) >= 12
        for i, c in enumerate(cusps[:12]):
            assert 0 <= c < 360, f"Year {year}: cusp {i + 1} = {c}"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "lat,lon,name",
        [
            (41.9, 12.5, "Rome"),
            (0.0, 0.0, "Equator"),
            (-33.9, 151.2, "Sydney"),
            (55.7, 37.6, "Moscow"),
        ],
    )
    def test_sidereal_various_locations(self, lat: float, lon: float, name: str):
        """Sidereal houses valid at various locations."""
        jd = 2451545.0
        swe.set_sid_mode(SIDM_J2000)
        cusps, _ = swe.houses_ex(jd, lat, lon, ord("R"), FLG_SIDEREAL)
        assert len(cusps) >= 12
        for c in cusps[:12]:
            assert 0 <= c < 360, f"{name}: cusp out of range"
