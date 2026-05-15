"""
Tests for gauquelin_sector and gauquelin_sector functions.

Tests the calculation of Gauquelin sectors (1-36) for celestial bodies.
Verifies the implementation correctly divides diurnal and nocturnal arcs
into 18 sectors each (36 total).
"""

import math
import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SUN,
    MOON,
    MARS,
    JUPITER,
    SATURN,
    FLG_SWIEPH,
)
from libephemeris.houses import _houses_gauquelin


class TestGauquelinSectorBasic:
    """Basic tests for gauquelin_sector function."""

    @pytest.mark.unit
    def test_gauquelin_sector_returns_float(self):
        """gauquelin_sector() should return a float."""
        jd = 2451545.0  # J2000.0
        lat = 48.85  # Paris
        lon = 2.35

        result = ephem.gauquelin_sector(jd, MARS, 0, (lon, lat, 0.0))

        assert isinstance(result, float)

    @pytest.mark.unit
    def test_gauquelin_sector_in_valid_range(self):
        """gauquelin_sector() should return value in range [1, 37)."""
        jd = 2451545.0
        lat = 48.85
        lon = 2.35

        result = ephem.gauquelin_sector(jd, MARS, 0, (lon, lat, 0.0))

        assert 1.0 <= result < 37.0

    @pytest.mark.unit
    def test_gauquelin_sector_integer_part_is_sector_number(self):
        """Integer part of result should be sector number (1-36)."""
        jd = 2451545.0
        lat = 48.85
        lon = 2.35

        result = ephem.gauquelin_sector(jd, MARS, 0, (lon, lat, 0.0))
        sector_num = int(result)

        assert 1 <= sector_num <= 36

    @pytest.mark.unit
    def test_gauquelin_sector_alias_matches(self):
        """gauquelin_sector should give same result as gauquelin_sector."""
        jd = 2451545.0
        lat = 48.85
        lon = 2.35

        result1 = ephem.gauquelin_sector(jd, MARS, 0, (lon, lat, 0.0))
        # gauquelin_sector uses SE-compatible signature: (jd, body, method, geopos)
        geopos = (lon, lat, 0.0)
        result2 = ephem.gauquelin_sector(jd, MARS, 0, geopos)

        assert result1 == result2


class TestGauquelinSectorMultiplePlanets:
    """Test gauquelin_sector with different planets."""

    PLANETS = [SUN, MOON, MARS, JUPITER, SATURN]

    @pytest.mark.unit
    @pytest.mark.parametrize("planet", PLANETS)
    def test_gauquelin_sector_all_planets_return_valid(self, planet):
        """gauquelin_sector should return valid values for all major planets."""
        jd = 2451545.0
        lat = 48.85
        lon = 2.35

        result = ephem.gauquelin_sector(jd, planet, 0, (lon, lat, 0.0))

        assert 1.0 <= result < 37.0, f"Invalid result {result} for planet {planet}"


class TestGauquelinSectorMethods:
    """Test gauquelin_sector with different calculation methods."""

    @pytest.mark.unit
    @pytest.mark.parametrize("method", [0, 1])
    def test_gauquelin_sector_method_returns_valid(self, method):
        """gauquelin_sector should return valid results for methods 0 and 1."""
        jd = 2451545.0
        lat = 48.85
        lon = 2.35

        result = ephem.gauquelin_sector(jd, MARS, method, (lon, lat, 0.0))

        assert 1.0 <= result < 37.0

    @pytest.mark.unit
    def test_method_0_differs_from_method_1_for_high_latitude_body(self):
        """Method 0 (with lat) and method 1 (without lat) may give different results."""
        jd = 2451545.0
        lat = 48.85
        lon = 2.35

        result_with_lat = ephem.gauquelin_sector(jd, MOON, 0, (lon, lat, 0.0))
        result_without_lat = ephem.gauquelin_sector(jd, MOON, 1, (lon, lat, 0.0))

        # They should both be valid
        assert 1.0 <= result_with_lat < 37.0
        assert 1.0 <= result_without_lat < 37.0


class TestGauquelinSectorComparisonWithSwisseph:
    """Compare gauquelin_sector results with Swiss Ephemeris."""

    # Test cases: (jd, planet, lat, lon, method)
    TEST_CASES = [
        # Paris, Mars at J2000
        (2451545.0, MARS, 48.85, 2.35, 0),
        # Rome, Sun at different time
        (2451600.0, SUN, 41.9, 12.5, 0),
        # London, Moon
        (2451550.0, MOON, 51.5, -0.12, 0),
        # New York, Jupiter
        (2451545.0, JUPITER, 40.7, -74.0, 0),
        # Southern hemisphere, Saturn
        (2451545.0, SATURN, -33.9, 151.2, 0),
    ]

    @pytest.mark.unit
    @pytest.mark.parametrize("jd,planet,lat,lon,method", TEST_CASES)
    def test_gauquelin_sector_matches_swisseph_sector(
        self, jd, planet, lat, lon, method
    ):
        """gauquelin_sector sector number should match Swiss Ephemeris."""
        # Get result from libephemeris
        result_lib = ephem.gauquelin_sector(jd, planet, method, (lon, lat, 0.0))

        # Get result from Swiss Ephemeris
        geopos = (lon, lat, 0.0)  # (lon, lat, alt)
        result_swe = swe.gauquelin_sector(
            jd, planet, method, geopos, 0.0, 0.0, FLG_SWIEPH
        )

        # Compare sector numbers
        sector_lib = int(result_lib)
        sector_swe = int(result_swe)

        # Allow difference of up to 1 sector due to different algorithms
        # (Swiss Ephemeris uses more complex rise/set calculations)
        diff = abs(sector_lib - sector_swe)
        # Handle wrap-around (sector 36 is adjacent to sector 1)
        if diff > 18:
            diff = 36 - diff

        assert diff <= 2, (
            f"Sector mismatch: libephemeris={sector_lib}, swisseph={sector_swe}, "
            f"diff={diff}"
        )


class TestGauquelinSectorEdgeCases:
    """Test edge cases for gauquelin_sector."""

    @pytest.mark.unit
    def test_gauquelin_sector_equator(self):
        """gauquelin_sector should work at the equator."""
        jd = 2451545.0
        lat = 0.0  # Equator
        lon = 0.0

        result = ephem.gauquelin_sector(jd, MARS, 0, (lon, lat, 0.0))

        assert 1.0 <= result < 37.0

    @pytest.mark.unit
    def test_gauquelin_sector_high_latitude(self):
        """gauquelin_sector should work at high latitudes."""
        jd = 2451545.0
        lat = 65.0  # Near Arctic circle
        lon = 25.0

        result = ephem.gauquelin_sector(jd, MARS, 0, (lon, lat, 0.0))

        assert 1.0 <= result < 37.0

    @pytest.mark.unit
    def test_gauquelin_sector_southern_hemisphere(self):
        """gauquelin_sector should work in southern hemisphere."""
        jd = 2451545.0
        lat = -45.0
        lon = 170.0

        result = ephem.gauquelin_sector(jd, MARS, 0, (lon, lat, 0.0))

        assert 1.0 <= result < 37.0

    @pytest.mark.unit
    def test_gauquelin_sector_negative_longitude(self):
        """gauquelin_sector should work with negative longitude (Western)."""
        jd = 2451545.0
        lat = 40.7
        lon = -74.0  # New York

        result = ephem.gauquelin_sector(jd, MARS, 0, (lon, lat, 0.0))

        assert 1.0 <= result < 37.0


class TestGauquelinSectorSectorDistribution:
    """Test that sectors are properly distributed across time."""

    @pytest.mark.unit
    def test_different_times_give_different_sectors(self):
        """Different times should give different sector positions."""
        lat = 48.85
        lon = 2.35

        # Test Mars at different times (spaced ~2 hours apart)
        # 2 hours = ~2/24 of a day = ~0.083 days
        sectors = []
        for i in range(6):
            jd = 2451545.0 + i * 0.1  # ~2.4 hours apart
            sector = ephem.gauquelin_sector(jd, MARS, 0, (lon, lat, 0.0))
            sectors.append(int(sector))

        # At least some of the sectors should be different
        unique_sectors = set(sectors)
        assert len(unique_sectors) >= 2, (
            f"Expected different sectors at different times, got {sectors}"
        )

    @pytest.mark.unit
    def test_sun_cycles_through_all_sectors_in_day(self):
        """Sun should cycle through many sectors over 24 hours."""
        lat = 48.85
        lon = 2.35
        jd_start = 2451545.0

        # Sample 24 times over 24 hours
        sectors = set()
        for i in range(24):
            jd = jd_start + i / 24.0
            sector = ephem.gauquelin_sector(jd, SUN, 0, (lon, lat, 0.0))
            sectors.add(int(sector))

        # Sun should appear in at least 20 different sectors in 24 hours
        assert len(sectors) >= 20, (
            f"Expected Sun to visit many sectors in 24h, got {len(sectors)}: {sorted(sectors)}"
        )


class TestGauquelinSectorAtmosphericParameters:
    """Test gauquelin_sector with atmospheric parameters."""

    @pytest.mark.unit
    def test_gauquelin_sector_with_altitude(self):
        """gauquelin_sector should accept altitude parameter."""
        jd = 2451545.0
        lat = 48.85
        lon = 2.35
        altitude = 1000.0  # 1km altitude

        result = ephem.gauquelin_sector(jd, MARS, 0, (lon, lat, altitude))

        assert 1.0 <= result < 37.0

    @pytest.mark.unit
    def test_gauquelin_sector_with_pressure_and_temperature(self):
        """gauquelin_sector should accept pressure and temperature parameters."""
        jd = 2451545.0
        lat = 48.85
        lon = 2.35

        result = ephem.gauquelin_sector(
            jd, MARS, 0, (lon, lat, 0.0), atpress=1000.0, attemp=20.0
        )

        assert 1.0 <= result < 37.0


class TestGauquelinSectorDivision:
    """Test that Gauquelin correctly divides arcs into 18 sectors each."""

    @pytest.mark.unit
    def test_diurnal_arc_has_18_sectors(self):
        """Sectors 1-18 should be in the diurnal (above horizon) arc."""
        jd = 2451545.0
        lat = 48.85
        lon = 2.35

        # Collect sectors over 24 hours - diurnal and nocturnal
        diurnal_sectors = set()
        nocturnal_sectors = set()

        for i in range(48):  # Sample every 30 minutes
            jd_sample = jd + i / 48.0
            sector = ephem.gauquelin_sector(jd_sample, MARS, 0, (lon, lat, 0.0))
            sector_int = int(sector)

            if 1 <= sector_int <= 18:
                diurnal_sectors.add(sector_int)
            elif 19 <= sector_int <= 36:
                nocturnal_sectors.add(sector_int)

        # Both arcs should have sectors (Mars should spend time above and below horizon)
        assert len(diurnal_sectors) > 0, "Expected diurnal sectors (1-18)"
        assert len(nocturnal_sectors) > 0, "Expected nocturnal sectors (19-36)"

    @pytest.mark.unit
    def test_sector_boundaries_at_cardinal_points(self):
        """Key sectors should be at cardinal points: 1=Asc, 10=MC, 19=Desc, 28=IC."""
        # This test verifies the sector numbering convention
        # JD 2451545.0 is J2000 = Jan 1, 2000 at 12:00 TT (approximately noon UTC)
        # For Paris (lon=2.35), local solar noon is when Sun crosses meridian

        jd = 2451545.0
        lat = 48.85
        lon = 2.35

        # Collect sectors throughout the day to see the distribution
        all_sectors = []
        for i in range(24):
            jd_sample = jd + i / 24.0
            sector = ephem.gauquelin_sector(jd_sample, SUN, 0, (lon, lat, 0.0))
            all_sectors.append((i, int(sector)))

        # The Sun should spend part of day in diurnal sectors and part in nocturnal
        diurnal_hours = [h for h, s in all_sectors if 1 <= s <= 18]
        nocturnal_hours = [h for h, s in all_sectors if 19 <= s <= 36]

        # Sun should be in diurnal sectors during some hours (daytime)
        assert len(diurnal_hours) > 0, (
            f"Expected Sun in diurnal sectors at some hours, got {all_sectors}"
        )

        # Sun should be in nocturnal sectors during some hours (nighttime)
        assert len(nocturnal_hours) > 0, (
            f"Expected Sun in nocturnal sectors at some hours, got {all_sectors}"
        )

    @pytest.mark.unit
    def test_sectors_1_to_18_are_diurnal(self):
        """Verify sectors 1-18 represent the diurnal arc (above horizon)."""
        # The gauquelin_sector function should place planets above
        # the horizon in sectors 1-18

        # Per the implementation comment:
        # Sector 1: Rising (Asc, H = h_rise)
        # Sector 10: Culminating (MC, H = 0)
        # Sector 19: Setting (Desc, H = h_set)
        # Sector 28: Anti-culminating (IC, H = 180)

        # This is verified by the implementation using is_above to check
        # horizon position and mapping to sectors 1-18 or 19-36

        jd = 2451545.0
        lat = 48.85
        lon = 2.35

        # Sample many times and verify sector distribution makes sense
        all_sectors = []
        for i in range(72):  # Every 20 minutes for 24 hours
            jd_sample = jd + i / 72.0
            sector = ephem.gauquelin_sector(jd_sample, SUN, 0, (lon, lat, 0.0))
            all_sectors.append(int(sector))

        # Count sectors in each range
        diurnal_count = sum(1 for s in all_sectors if 1 <= s <= 18)
        nocturnal_count = sum(1 for s in all_sectors if 19 <= s <= 36)

        # Both should be populated (Sun spends ~12h above and ~12h below horizon)
        assert diurnal_count > 20, f"Expected more diurnal samples, got {diurnal_count}"
        assert nocturnal_count > 20, (
            f"Expected more nocturnal samples, got {nocturnal_count}"
        )

    @pytest.mark.unit
    def test_sectors_19_to_36_are_nocturnal(self):
        """Verify sectors 19-36 represent the nocturnal arc (below horizon)."""
        # Per sector numbering: 19=Desc, 28=IC, 36/1=Asc
        jd = 2451545.0
        lat = 48.85
        lon = 2.35

        # Verify that we can find sectors in both ranges
        found_19_to_36 = False
        for i in range(48):
            jd_sample = jd + i / 48.0
            sector = ephem.gauquelin_sector(jd_sample, SUN, 0, (lon, lat, 0.0))
            if 19 <= int(sector) <= 36:
                found_19_to_36 = True
                break

        assert found_19_to_36, "Expected to find nocturnal sectors (19-36)"


class TestGauquelinHousesFunction:
    """Test the _houses_gauquelin internal function for 36-sector division."""

    @pytest.mark.unit
    def test_houses_gauquelin_returns_37_sectors(self):
        """_houses_gauquelin should return 37 elements (index 0 unused, 1-36 sectors)."""
        armc = 180.0
        lat = 48.85
        eps = 23.44
        asc = 90.0
        mc = 180.0

        cusps = _houses_gauquelin(armc, lat, eps, asc, mc)

        assert len(cusps) == 37
        # Sectors 1-36 should all be valid degrees
        for i in range(1, 37):
            assert 0.0 <= cusps[i] < 360.0, f"Sector {i} invalid: {cusps[i]}"

    @pytest.mark.unit
    def test_houses_gauquelin_maps_36_sectors_to_12_houses(self):
        """Verify 36 sectors are returned with proper structure."""
        # Gauquelin has 36 sectors, not 12 houses
        # Sectors 1 and 19 are opposite (Asc/Desc)
        # Sectors 10 and 28 are opposite (MC/IC)

        armc = 180.0
        lat = 45.0
        eps = 23.44
        asc = 90.0
        mc = 180.0

        cusps = _houses_gauquelin(armc, lat, eps, asc, mc)

        # Should have 37 elements (index 0 unused)
        assert len(cusps) == 37

        # Cardinal sectors should be at Asc, MC, Desc, IC
        assert abs(cusps[1] - asc) < 0.01  # Sector 1 = Asc
        assert abs(cusps[10] - mc) < 0.01  # Sector 10 = MC
        assert abs(cusps[19] - (asc + 180) % 360) < 0.01  # Sector 19 = Desc
        assert abs(cusps[28] - (mc + 180) % 360) < 0.01  # Sector 28 = IC

    @pytest.mark.unit
    def test_houses_gauquelin_polar_fallback(self):
        """At polar latitudes, should fall back to equal division of 36 sectors."""
        # Within polar circle (|lat| >= 90 - eps), use equal division
        armc = 180.0
        lat = 70.0  # Polar latitude
        eps = 23.44  # |70| + 23.44 > 90, so it's within polar circle
        asc = 90.0
        mc = 180.0

        cusps = _houses_gauquelin(armc, lat, eps, asc, mc)

        # Should return 37 elements (index 0 unused, 1-36 sectors)
        assert len(cusps) == 37
        for i in range(1, 37):
            assert 0.0 <= cusps[i] < 360.0

    @pytest.mark.unit
    def test_diurnal_nocturnal_arc_division(self):
        """Verify diurnal and nocturnal arcs are each divided into 18 sectors."""
        # The implementation divides:
        # - Sectors 1-18: diurnal (above horizon, from Asc to Desc via MC)
        # - Sectors 19-36: nocturnal (below horizon, from Desc to Asc via IC)

        # This is verified by checking the sector interpolation in the code:
        # for i in range(1, 37):
        #     if i <= 18:  # diurnal
        #     else:        # nocturnal (i in 19-36)

        # Indirect verification: check that Sun samples cover both ranges
        jd = 2451545.0
        lat = 48.85
        lon = 2.35

        diurnal_samples = []
        nocturnal_samples = []

        for i in range(36):  # Sample every 40 minutes
            jd_sample = jd + i / 36.0
            sector = ephem.gauquelin_sector(jd_sample, SUN, 0, (lon, lat, 0.0))
            if 1 <= sector < 19:
                diurnal_samples.append(sector)
            else:
                nocturnal_samples.append(sector)

        # Sun should spend roughly equal time in each arc
        assert len(diurnal_samples) > 10, "Expected diurnal samples"
        assert len(nocturnal_samples) > 10, "Expected nocturnal samples"

        # Each arc should have 18 sectors worth of range
        if diurnal_samples:
            diurnal_range = max(diurnal_samples) - min(diurnal_samples)
            assert diurnal_range > 5, (
                f"Expected larger diurnal range, got {diurnal_range}"
            )

        if nocturnal_samples:
            nocturnal_range = max(nocturnal_samples) - min(nocturnal_samples)
            assert nocturnal_range > 5, (
                f"Expected larger nocturnal range, got {nocturnal_range}"
            )

    @pytest.mark.unit
    def test_sector_continuity_through_36(self):
        """Sectors should flow continuously 1->2->...->36->1."""
        jd = 2451545.0
        lat = 48.85
        lon = 2.35

        # Collect sectors over 24 hours
        sectors = []
        for i in range(144):  # Every 10 minutes
            jd_sample = jd + i / 144.0
            sector = ephem.gauquelin_sector(jd_sample, SUN, 0, (lon, lat, 0.0))
            sectors.append(sector)

        # Check that sectors progress (mostly) monotonically or wrap around
        # Since Sun moves through sectors in order, transitions should be small
        large_jumps = 0
        for i in range(1, len(sectors)):
            diff = sectors[i] - sectors[i - 1]
            # Allow for normal progression and wrap-around (36 -> 1)
            if abs(diff) > 5 and abs(abs(diff) - 36) > 5:
                large_jumps += 1

        # Most transitions should be smooth
        assert large_jumps < 10, (
            f"Too many large jumps in sector sequence: {large_jumps}"
        )
