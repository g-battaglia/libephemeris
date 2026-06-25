"""
Tests for solar and lunar eclipse functions.

Verifies that eclipse search functions return valid results,
correct eclipse types, and plausible timing values.
"""

from __future__ import annotations


import pytest

import libephemeris as swe
from libephemeris.constants import (
    ECL_TOTAL,
    ECL_ANNULAR,
    ECL_PARTIAL,
    ECL_PENUMBRAL,
)


# Known solar eclipses for validation
KNOWN_SOLAR_ECLIPSES = [
    # (search_jd, expected_year, expected_month, eclipse_type_description)
    (2451545.0, 2000, 2, "after J2000"),  # First eclipse after J2000
    (2458849.5, 2020, 1, "after 2020-Jan"),  # First eclipse after 2020
    (2460310.5, 2024, 1, "after 2024-Jan"),  # First eclipse after 2024
]

KNOWN_LUNAR_ECLIPSES = [
    (2451545.0, 2000, 1, "after J2000"),
    (2458849.5, 2020, 1, "after 2020-Jan"),
]


class TestSolarEclipseGlobal:
    """Tests for sol_eclipse_when_glob."""

    @pytest.mark.unit
    def test_sol_eclipse_returns_valid_format(self):
        """sol_eclipse_when_glob returns (retflag, tret_10_tuple)."""
        jd = 2451545.0
        result = swe.sol_eclipse_when_glob(jd)
        assert len(result) == 2, f"Expected 2 return values, got {len(result)}"
        retflag, tret = result
        assert isinstance(retflag, int)
        assert len(tret) == 10, f"tret should have 10 elements, got {len(tret)}"

    @pytest.mark.unit
    def test_sol_eclipse_retflag_has_type(self):
        """Return flag should indicate an eclipse type."""
        jd = 2451545.0
        retflag, tret = swe.sol_eclipse_when_glob(jd)
        # Should have at least one type bit set
        type_bits = ECL_TOTAL | ECL_ANNULAR | ECL_PARTIAL
        assert retflag & type_bits, f"retflag {retflag} has no eclipse type bits"

    @pytest.mark.unit
    def test_sol_eclipse_time_after_search(self):
        """Eclipse maximum time should be after search start."""
        jd = 2451545.0
        retflag, tret = swe.sol_eclipse_when_glob(jd)
        eclipse_max = tret[0]
        assert eclipse_max > jd, (
            f"Eclipse max {eclipse_max} not after search start {jd}"
        )

    @pytest.mark.unit
    def test_sol_eclipse_time_reasonable(self):
        """Eclipse should be found within ~1 year of search start."""
        jd = 2451545.0
        retflag, tret = swe.sol_eclipse_when_glob(jd)
        eclipse_max = tret[0]
        # Solar eclipses happen ~2-5 times per year
        assert eclipse_max - jd < 365.25, (
            f"Eclipse {eclipse_max - jd:.1f} days after search start"
        )

    @pytest.mark.unit
    def test_sol_eclipse_sequential_search(self):
        """Searching from eclipse time+1 should find the next eclipse."""
        jd = 2451545.0
        retflag1, tret1 = swe.sol_eclipse_when_glob(jd)
        retflag2, tret2 = swe.sol_eclipse_when_glob(tret1[0] + 1.0)
        assert tret2[0] > tret1[0], "Second eclipse should be after first"
        # Gap between successive eclipses should be < 1 year
        gap = tret2[0] - tret1[0]
        assert gap < 365.25, f"Gap between eclipses: {gap:.1f} days"

    @pytest.mark.unit
    def test_sol_eclipse_find_5_eclipses(self):
        """Find 5 successive solar eclipses starting from J2000."""
        jd = 2451545.0
        eclipses = []
        for _ in range(5):
            retflag, tret = swe.sol_eclipse_when_glob(jd)
            eclipses.append((tret[0], retflag))
            jd = tret[0] + 1.0

        # Verify ascending order
        for i in range(4):
            assert eclipses[i + 1][0] > eclipses[i][0]

        # Verify all have type bits
        type_bits = ECL_TOTAL | ECL_ANNULAR | ECL_PARTIAL
        for t, flag in eclipses:
            assert flag & type_bits, f"Eclipse at JD {t} has no type"


class TestLunarEclipse:
    """Tests for lun_eclipse_when."""

    @pytest.mark.unit
    def test_lun_eclipse_returns_valid_format(self):
        """lun_eclipse_when returns (retflag, tret_10_tuple)."""
        jd = 2451545.0
        result = swe.lun_eclipse_when(jd)
        assert len(result) == 2
        retflag, tret = result
        assert isinstance(retflag, int)
        assert len(tret) == 10

    @pytest.mark.unit
    def test_lun_eclipse_retflag_has_type(self):
        """Return flag should indicate a lunar eclipse type."""
        jd = 2451545.0
        retflag, tret = swe.lun_eclipse_when(jd)
        type_bits = ECL_TOTAL | ECL_PARTIAL | ECL_PENUMBRAL
        assert retflag & type_bits, f"retflag {retflag} has no lunar eclipse type bits"

    @pytest.mark.unit
    def test_lun_eclipse_time_after_search(self):
        """Eclipse time should be after search start."""
        jd = 2451545.0
        retflag, tret = swe.lun_eclipse_when(jd)
        assert tret[0] > jd

    @pytest.mark.unit
    def test_lun_eclipse_find_5_eclipses(self):
        """Find 5 successive lunar eclipses."""
        jd = 2451545.0
        eclipses = []
        for _ in range(5):
            retflag, tret = swe.lun_eclipse_when(jd)
            eclipses.append((tret[0], retflag))
            jd = tret[0] + 1.0

        for i in range(4):
            assert eclipses[i + 1][0] > eclipses[i][0]


class TestSolarEclipseLocal:
    """Tests for sol_eclipse_when_loc."""

    @pytest.mark.unit
    def test_sol_eclipse_loc_returns_valid(self):
        """sol_eclipse_when_loc returns valid results."""
        jd = 2451545.0
        # sol_eclipse_when_loc takes geopos as (lon, lat, alt)
        geopos = (12.5, 41.9, 0.0)  # Rome
        result = swe.sol_eclipse_when_loc(jd, geopos)
        assert len(result) >= 2
        # First element is retflag or similar

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "lat,lon,name",
        [
            (41.9, 12.5, "Rome"),
            (40.7, -74.0, "New York"),
            (35.7, 139.7, "Tokyo"),
            (-33.9, 151.2, "Sydney"),
            (0.0, 0.0, "Equator"),
        ],
    )
    def test_sol_eclipse_loc_various_locations(self, lat: float, lon: float, name: str):
        """Eclipse search works at various locations."""
        jd = 2451545.0
        # sol_eclipse_when_loc takes geopos as (lon, lat, alt)
        geopos = (lon, lat, 0.0)
        result = swe.sol_eclipse_when_loc(jd, geopos)
        assert len(result) >= 2, f"{name}: unexpected result length"


class TestEclipseTimingConsistency:
    """Test that eclipse timing values are internally consistent."""

    @pytest.mark.unit
    def test_sol_eclipse_tret_ordering(self):
        """Solar eclipse tret values should have reasonable ordering."""
        jd = 2451545.0
        retflag, tret = swe.sol_eclipse_when_glob(jd)
        # tret[0] = max eclipse, tret[2] = begin, tret[3] = end
        if tret[2] > 0 and tret[3] > 0:
            assert tret[2] < tret[0] < tret[3], (
                f"Eclipse ordering: begin={tret[2]} max={tret[0]} end={tret[3]}"
            )

    @pytest.mark.unit
    def test_lun_eclipse_tret_ordering(self):
        """Lunar eclipse tret values should have reasonable ordering."""
        jd = 2451545.0
        retflag, tret = swe.lun_eclipse_when(jd)
        if tret[2] > 0 and tret[3] > 0:
            assert tret[2] < tret[0] < tret[3], (
                f"Eclipse ordering: begin={tret[2]} max={tret[0]} end={tret[3]}"
            )

    @pytest.mark.unit
    def test_eclipse_duration_reasonable(self):
        """Eclipse duration should be less than 12 hours."""
        jd = 2451545.0
        retflag, tret = swe.sol_eclipse_when_glob(jd)
        if tret[2] > 0 and tret[3] > 0:
            duration_hours = (tret[3] - tret[2]) * 24
            assert duration_hours < 12, (
                f"Eclipse duration {duration_hours:.1f} hours too long"
            )
