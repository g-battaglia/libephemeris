"""
Tests for lunar eclipse location-based calculations in libephemeris.

Tests the lun_eclipse_when_loc function which finds lunar eclipses
visible from a specific geographic location.

Reference data from NASA Eclipse website:
https://eclipse.gsfc.nasa.gov/lunar.html

Times tuple layout (10 elements, pyswisseph-compatible):
    [0]: Time of maximum eclipse
    [1]: Reserved
    [2]: Time of partial eclipse beginning (Moon enters umbra)
    [3]: Time of partial eclipse ending (Moon leaves umbra)
    [4]: Time of totality beginning
    [5]: Time of totality ending
    [6]: Time of penumbral eclipse beginning
    [7]: Time of penumbral eclipse ending
    [8]: Time of moonrise (if Moon rises during eclipse)
    [9]: Time of moonset (if Moon sets during eclipse)

Attr tuple layout (20 elements, pyswisseph-compatible):
    [0]: Umbral magnitude
    [1]: Penumbral magnitude
    [2]: Reserved
    [3]: Reserved
    [4]: Azimuth of Moon at maximum (degrees)
    [5]: True altitude of Moon at maximum (degrees)
    [6]: Apparent altitude of Moon (with refraction)
    [7]: Distance from opposition (degrees)
    [8]: Umbral magnitude (equals [0])
    [9]: Saros series number
    [10]: Saros series member number
    [11-19]: Reserved
"""

from libephemeris import (
    julday,
    revjul,
    lun_eclipse_when_loc,
    ECL_TOTAL,
    ECL_VISIBLE,
    ECL_MAX_VISIBLE,
    ECL_1ST_VISIBLE,
    ECL_2ND_VISIBLE,
    ECL_3RD_VISIBLE,
    ECL_4TH_VISIBLE,
)


class TestLunEclipseWhenLoc:
    """Test suite for lun_eclipse_when_loc function."""

    def test_finds_lunar_eclipse_from_rome(self):
        """Test that function finds a lunar eclipse visible from Rome."""
        # Start from Jan 1, 2024
        jd_start = julday(2024, 1, 1, 0)
        rome_lat, rome_lon = 41.9028, 12.4964

        ecl_type, times, attr = lun_eclipse_when_loc(
            jd_start, (rome_lon, rome_lat, 0.0)
        )

        # Should find an eclipse
        assert ecl_type != 0
        assert times[0] > jd_start  # Maximum should be after start
        assert times[0] > 0  # Should have valid maximum time
        # Should be marked as visible
        assert ecl_type & ECL_VISIBLE

    def test_returns_correct_tuple_structure(self):
        """Test that return values have correct structure."""
        jd_start = julday(2024, 1, 1, 0)
        london_lat, london_lon = 51.5074, -0.1278

        ecl_type, times, attr = lun_eclipse_when_loc(
            jd_start, (london_lon, london_lat, 0.0)
        )

        # Should return 10-element times tuple
        assert len(times) == 10
        # All elements should be floats
        assert all(isinstance(t, float) for t in times)

        # Should return 20-element attr tuple (pyswisseph layout)
        assert len(attr) == 20

        # Eclipse type should be int
        assert isinstance(ecl_type, int)

    def test_attributes_have_valid_values(self):
        """Test that eclipse attributes have reasonable values."""
        jd_start = julday(2024, 1, 1, 0)
        tokyo_lat, tokyo_lon = 35.6762, 139.6503

        ecl_type, times, attr = lun_eclipse_when_loc(
            jd_start, (tokyo_lon, tokyo_lat, 0.0)
        )

        # Umbral magnitude should be between 0 and ~2 for visible eclipses
        assert 0 <= attr[0] <= 2.5

        # Penumbral magnitude should be positive
        assert attr[1] >= 0

        # Azimuth should be 0-360 (pyswisseph index [4])
        assert 0 <= attr[4] <= 360

        # Altitude could be any value but should be reasonable (pyswisseph index [5])
        assert -90 <= attr[5] <= 90

    def test_phase_times_ordering(self):
        """Test that phase times are in correct chronological order."""
        jd_start = julday(2024, 1, 1, 0)
        sydney_lat, sydney_lon = -33.8688, 151.2093

        ecl_type, times, attr = lun_eclipse_when_loc(
            jd_start, (sydney_lon, sydney_lat, 0.0)
        )

        # Penumbral begin should be first (if present)
        if times[6] > 0:
            assert times[6] < times[0]  # Penumbral begin < maximum

        # Penumbral end should be last (if present)
        if times[7] > 0:
            assert times[7] > times[0]  # Penumbral end > maximum

        # Partial phases should be between penumbral phases
        if times[2] > 0 and times[6] > 0:
            assert times[2] >= times[6]  # Partial begin >= penumbral begin
        if times[3] > 0 and times[7] > 0:
            assert times[3] <= times[7]  # Partial end <= penumbral end

    def test_different_locations_may_see_different_eclipses(self):
        """Test that different locations might see different eclipses."""
        jd_start = julday(2024, 1, 1, 0)

        # Two locations on opposite sides of the Earth
        new_york_lat, new_york_lon = 40.7128, -74.0060
        beijing_lat, beijing_lon = 39.9042, 116.4074

        ecl_type_ny, times_ny, attr_ny = lun_eclipse_when_loc(
            jd_start, (new_york_lon, new_york_lat, 0.0)
        )
        ecl_type_bj, times_bj, attr_bj = lun_eclipse_when_loc(
            jd_start, (beijing_lon, beijing_lat, 0.0)
        )

        # Both should find eclipses (may be same or different)
        assert ecl_type_ny & ECL_VISIBLE
        assert ecl_type_bj & ECL_VISIBLE

    def test_visibility_flags_are_set(self):
        """Test that visibility flags are properly set."""
        jd_start = julday(2024, 1, 1, 0)
        berlin_lat, berlin_lon = 52.5200, 13.4050

        ecl_type, times, attr = lun_eclipse_when_loc(
            jd_start, (berlin_lon, berlin_lat, 0.0)
        )

        # ECL_VISIBLE should always be set for returned eclipses
        assert ecl_type & ECL_VISIBLE

        # At least some visibility flags should be set. The lunar bits
        # span MAX plus the six phase bits including the penumbral ones
        # (ECL_PENUMBBEG/PENUMBEND at 8192/16384).
        from libephemeris import (
            ECL_PENUMBBEG_VISIBLE,
            ECL_PENUMBEND_VISIBLE,
        )

        visibility_flags = (
            ECL_MAX_VISIBLE
            | ECL_1ST_VISIBLE
            | ECL_2ND_VISIBLE
            | ECL_3RD_VISIBLE
            | ECL_4TH_VISIBLE
            | ECL_PENUMBBEG_VISIBLE
            | ECL_PENUMBEND_VISIBLE
        )
        # At least one phase should be visible
        assert ecl_type & visibility_flags

    def test_moon_altitude_is_valid(self):
        """Test that Moon altitude is reasonable for visible eclipse."""
        jd_start = julday(2024, 1, 1, 0)
        paris_lat, paris_lon = 48.8566, 2.3522

        ecl_type, times, attr = lun_eclipse_when_loc(
            jd_start, (paris_lon, paris_lat, 0.0)
        )

        moon_alt = attr[5]  # True altitude (pyswisseph index)

        # For a visible eclipse, Moon should be above horizon at some point
        # The altitude at maximum might be below horizon if eclipse is partially visible
        # during moonrise/moonset, but it should still be reasonable
        assert -90 <= moon_alt <= 90

    def test_swe_alias(self):
        """Test that lun_eclipse_when_loc matches lun_eclipse_when_loc."""
        jd_start = julday(2024, 1, 1, 0)
        cape_town_lat, cape_town_lon = -33.9249, 18.4241

        ecl_type1, times1, attr1 = lun_eclipse_when_loc(
            jd_start, (cape_town_lon, cape_town_lat, 0.0)
        )
        # swe_ version takes geopos tuple: [lon, lat, alt]
        geopos = [cape_town_lon, cape_town_lat, 0]
        ecl_type2, times2, attr2 = lun_eclipse_when_loc(jd_start, geopos)

        assert times1 == times2
        assert attr1 == attr2
        assert ecl_type1 == ecl_type2

    def test_sequential_visible_eclipses(self):
        """Test finding multiple sequential visible lunar eclipses."""
        jd = julday(2024, 1, 1, 0)
        mumbai_lat, mumbai_lon = 19.0760, 72.8777
        eclipses = []

        # Find 3 sequential visible eclipses
        for _ in range(3):
            ecl_type, times, attr = lun_eclipse_when_loc(
                jd, (mumbai_lon, mumbai_lat, 0.0)
            )
            eclipses.append((times[0], ecl_type))
            jd = times[0] + 1  # Start after this eclipse

        # Each eclipse should be later than the previous
        for i in range(1, len(eclipses)):
            assert eclipses[i][0] > eclipses[i - 1][0]

    def test_known_eclipse_may_2022_visible_from_americas(self):
        """Test the total lunar eclipse of May 16, 2022 visible from Americas.

        Reference: NASA - Total Lunar Eclipse of May 16, 2022
        Maximum eclipse: approximately 04:12 UTC
        Visible from North and South America, Europe, Africa
        """
        # Start searching from May 1, 2022
        jd_start = julday(2022, 5, 1, 0)
        # Rio de Janeiro - should see this eclipse
        rio_lat, rio_lon = -22.9068, -43.1729

        ecl_type, times, attr = lun_eclipse_when_loc(jd_start, (rio_lon, rio_lat, 0.0))

        # Should be a total eclipse
        assert ecl_type & ECL_TOTAL
        assert ecl_type & ECL_VISIBLE

        # Check maximum is on May 16, 2022
        year, month, day, hour = revjul(times[0])
        assert year == 2022
        assert month == 5
        assert day == 16

    def test_altitude_parameter(self):
        """Test that altitude parameter is accepted and doesn't cause errors."""
        jd_start = julday(2024, 1, 1, 0)
        denver_lat, denver_lon = 39.7392, -104.9903
        denver_altitude = 1609  # 1 mile high

        ecl_type, times, attr = lun_eclipse_when_loc(
            jd_start, (denver_lon, denver_lat, denver_altitude)
        )

        # Should find an eclipse
        assert ecl_type & ECL_VISIBLE
        assert times[0] > jd_start

    def test_southern_hemisphere(self):
        """Test finding eclipses from southern hemisphere."""
        jd_start = julday(2024, 1, 1, 0)
        # Melbourne, Australia
        melbourne_lat, melbourne_lon = -37.8136, 144.9631

        ecl_type, times, attr = lun_eclipse_when_loc(
            jd_start, (melbourne_lon, melbourne_lat, 0.0)
        )

        # Should find an eclipse
        assert ecl_type & ECL_VISIBLE
        assert times[0] > jd_start

    def test_high_latitude_location(self):
        """Test finding eclipses from high latitude location."""
        jd_start = julday(2024, 1, 1, 0)
        # Reykjavik, Iceland
        reykjavik_lat, reykjavik_lon = 64.1466, -21.9426

        ecl_type, times, attr = lun_eclipse_when_loc(
            jd_start, (reykjavik_lon, reykjavik_lat, 0.0)
        )

        # Should find an eclipse (might take longer due to polar location)
        assert ecl_type & ECL_VISIBLE
        assert times[0] > jd_start


class TestLunEclipseWhenLocEdgeCases:
    """Test edge cases for location-based lunar eclipse calculations."""

    def test_equatorial_location(self):
        """Test finding eclipses from equatorial location."""
        jd_start = julday(2024, 1, 1, 0)
        # Quito, Ecuador (nearly on the equator)
        quito_lat, quito_lon = -0.1807, -78.4678

        ecl_type, times, attr = lun_eclipse_when_loc(
            jd_start, (quito_lon, quito_lat, 0.0)
        )

        assert ecl_type & ECL_VISIBLE
        assert times[0] > jd_start

    def test_prime_meridian(self):
        """Test finding eclipses from location on prime meridian."""
        jd_start = julday(2024, 1, 1, 0)
        # Greenwich, UK
        greenwich_lat, greenwich_lon = 51.4772, 0.0

        ecl_type, times, attr = lun_eclipse_when_loc(
            jd_start, (greenwich_lon, greenwich_lat, 0.0)
        )

        assert ecl_type & ECL_VISIBLE
        assert times[0] > jd_start

    def test_international_date_line(self):
        """Test finding eclipses from location near date line."""
        jd_start = julday(2024, 1, 1, 0)
        # Fiji
        fiji_lat, fiji_lon = -18.1416, 178.4419

        ecl_type, times, attr = lun_eclipse_when_loc(
            jd_start, (fiji_lon, fiji_lat, 0.0)
        )

        assert ecl_type & ECL_VISIBLE
        assert times[0] > jd_start

    def test_early_20th_century(self):
        """Test finding eclipse in early 20th century."""
        # 1920 CE (within ephemeris range)
        jd_start = julday(1920, 1, 1, 0)
        london_lat, london_lon = 51.5074, -0.1278

        ecl_type, times, attr = lun_eclipse_when_loc(
            jd_start, (london_lon, london_lat, 0.0)
        )

        assert ecl_type & ECL_VISIBLE
        assert times[0] > jd_start

    def test_mid_21st_century(self):
        """Test finding eclipse in mid 21st century."""
        # 2050 CE (within ephemeris range)
        jd_start = julday(2050, 1, 1, 0)
        new_york_lat, new_york_lon = 40.7128, -74.0060

        ecl_type, times, attr = lun_eclipse_when_loc(
            jd_start, (new_york_lon, new_york_lat, 0.0)
        )

        assert ecl_type & ECL_VISIBLE
        assert times[0] > jd_start

    def test_moonrise_moonset_during_eclipse(self):
        """Test that moonrise/moonset times are calculated when applicable."""
        # This test verifies the function doesn't crash when calculating
        # moonrise/moonset times. The actual times depend on specific eclipse
        # geometry and location.
        jd_start = julday(2024, 1, 1, 0)
        moscow_lat, moscow_lon = 55.7558, 37.6173

        ecl_type, times, attr = lun_eclipse_when_loc(
            jd_start, (moscow_lon, moscow_lat, 0.0)
        )

        # times[8] is moonrise, times[9] is moonset (pyswisseph layout)
        # These may be 0 if Moon doesn't rise/set during eclipse
        assert isinstance(times[8], float)
        assert isinstance(times[9], float)
        # If set, they should be within a reasonable window around the eclipse
        if times[8] > 0:
            # Moonrise should be near the eclipse (within 1 day of maximum)
            assert abs(times[8] - times[0]) < 1.0
        if times[9] > 0:
            # Moonset should be near the eclipse (within 1 day of maximum)
            assert abs(times[9] - times[0]) < 1.0
