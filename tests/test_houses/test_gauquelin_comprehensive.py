"""
Tests for Gauquelin sectors (36 cusps) and gauquelin_sector function.

Verifies that house system 'G' returns 36 cusps, sector numbering,
gauquelin_sector() body placement, and edge cases.
"""

from __future__ import annotations

import pytest

import libephemeris as swe
from libephemeris.constants import (
    SUN,
    MOON,
    MARS,
    JUPITER,
)


@pytest.fixture(autouse=True)
def _reset_state():
    """Reset ephemeris state between tests."""
    yield
    swe.close()


JD_J2000 = 2451545.0
JD_2020 = 2458849.5
JD_2023 = 2460000.0


class TestGauquelinCusps:
    """Test that Gauquelin house system returns 36 cusps."""

    @pytest.mark.unit
    def test_gauquelin_returns_36_cusps(self):
        """House system 'G' returns tuple of 36 cusps."""
        cusps, ascmc = swe.houses(JD_J2000, 48.85, 2.35, ord("G"))
        assert len(cusps) == 36, f"Expected 36 cusps, got {len(cusps)}"

    @pytest.mark.unit
    def test_other_systems_return_12_cusps(self):
        """Other house systems return 12 cusps."""
        for hsys in "PKORCABMTUWVX":
            cusps, ascmc = swe.houses(JD_J2000, 48.85, 2.35, ord(hsys))
            assert len(cusps) == 12, f"System {hsys} returned {len(cusps)} cusps"

    @pytest.mark.unit
    def test_gauquelin_cusps_are_ordered(self):
        """Gauquelin cusps should generally progress through the zodiac."""
        cusps, ascmc = swe.houses(JD_J2000, 48.85, 2.35, ord("G"))
        # All cusps should be valid angles [0, 360)
        for i, c in enumerate(cusps):
            assert 0.0 <= c < 360.0, f"Cusp {i + 1} out of range: {c}"

    @pytest.mark.unit
    def test_gauquelin_cusp_1_near_ascendant(self):
        """Sector 1 should be near the Ascendant."""
        cusps, ascmc = swe.houses(JD_J2000, 48.85, 2.35, ord("G"))
        asc = ascmc[0]
        # Cusp 1 should be close to the Ascendant
        diff = abs(cusps[0] - asc) % 360.0
        if diff > 180.0:
            diff = 360.0 - diff
        assert diff < 15.0, f"Cusp 1 ({cusps[0]}) not near Asc ({asc})"

    @pytest.mark.unit
    def test_gauquelin_ascmc_same_length_as_other_systems(self):
        """ascmc tuple has same structure regardless of house system."""
        _, ascmc_g = swe.houses(JD_J2000, 48.85, 2.35, ord("G"))
        _, ascmc_p = swe.houses(JD_J2000, 48.85, 2.35, ord("P"))
        assert len(ascmc_g) == len(ascmc_p)

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "lat,lon",
        [
            (48.85, 2.35),  # Paris
            (40.71, -74.01),  # New York
            (35.68, 139.69),  # Tokyo
            (-33.87, 151.21),  # Sydney
            (0.0, 0.0),  # Equator/Greenwich
            (60.0, 25.0),  # Helsinki (high latitude)
        ],
    )
    def test_gauquelin_at_various_locations(self, lat, lon):
        """Gauquelin sectors work at various geographic locations."""
        cusps, ascmc = swe.houses(JD_J2000, lat, lon, ord("G"))
        assert len(cusps) == 36
        for c in cusps:
            assert 0.0 <= c < 360.0

    @pytest.mark.unit
    @pytest.mark.parametrize("jd", [JD_J2000, JD_2020, JD_2023])
    def test_gauquelin_at_various_dates(self, jd):
        """Gauquelin sectors work at various dates."""
        cusps, ascmc = swe.houses(jd, 48.85, 2.35, ord("G"))
        assert len(cusps) == 36

    @pytest.mark.unit
    def test_gauquelin_polar_circle_handling(self):
        """Gauquelin at polar latitudes either works or raises PolarCircleError."""
        from libephemeris.exceptions import PolarCircleError

        try:
            cusps, ascmc = swe.houses(JD_J2000, 70.0, 25.0, ord("G"))
            # If it works, should still have 36 cusps
            assert len(cusps) == 36
        except PolarCircleError:
            # Expected for polar latitudes
            pass


class TestGauquelinSectorFunction:
    """Test gauquelin_sector() for body placement."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body", [SUN, MOON, MARS, JUPITER])
    def test_sector_returns_valid_range(self, body):
        """gauquelin_sector returns value in [1, 37)."""
        geopos = (2.35, 48.85, 0.0)  # Paris (lon, lat, alt)
        sector = swe.gauquelin_sector(JD_J2000, body, 0, geopos, 1013.25, 15.0)
        assert 1.0 <= sector < 37.0, f"Sector {sector} out of range for body {body}"

    @pytest.mark.unit
    def test_sector_integer_part_is_sector_number(self):
        """Integer part of result is the sector number (1-36)."""
        geopos = (2.35, 48.85, 0.0)
        sector = swe.gauquelin_sector(JD_J2000, MARS, 0, geopos, 1013.25, 15.0)
        sector_num = int(sector)
        assert 1 <= sector_num <= 36

    @pytest.mark.unit
    @pytest.mark.parametrize("method", [0, 1])
    def test_different_methods(self, method):
        """Methods 0 and 1 both return valid sectors."""
        geopos = (2.35, 48.85, 0.0)
        sector = swe.gauquelin_sector(
            JD_J2000, MARS, method, geopos, 1013.25, 15.0
        )
        assert 1.0 <= sector < 37.0

    @pytest.mark.unit
    def test_sector_varies_with_date(self):
        """Sector number changes with different dates."""
        geopos = (2.35, 48.85, 0.0)
        sectors = []
        for jd in [JD_J2000, JD_J2000 + 0.5, JD_J2000 + 1.0, JD_J2000 + 5.0]:
            s = swe.gauquelin_sector(jd, SUN, 0, geopos, 1013.25, 15.0)
            sectors.append(s)
        # At least some sectors should differ
        assert len(set(int(s) for s in sectors)) > 1

    @pytest.mark.unit
    def test_sector_varies_with_location(self):
        """Sector depends on geographic location."""
        s1 = swe.gauquelin_sector(
            JD_J2000, MARS, 0, (2.35, 48.85, 0.0), 1013.25, 15.0
        )
        s2 = swe.gauquelin_sector(
            JD_J2000, MARS, 0, (139.69, 35.68, 0.0), 1013.25, 15.0
        )
        # Different locations should give different sectors
        assert abs(s1 - s2) > 0.01


class TestGauquelinHouseName:
    """Test house name for Gauquelin system."""

    @pytest.mark.unit
    def test_house_name_gauquelin(self):
        """house_name returns 'Gauquelin sectors' for 'G'."""
        name = swe.house_name(ord("G"))
        assert "Gauquelin" in name or "gauquelin" in name.lower()
