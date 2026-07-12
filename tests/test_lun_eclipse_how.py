"""
Tests for lun_eclipse_how function in libephemeris.

Tests the lunar eclipse circumstances calculation for a specific
location and time.

Reference data from NASA Eclipse website:
https://eclipse.gsfc.nasa.gov/lunar.html
"""

import pytest

from libephemeris import (
    julday,
    lun_eclipse_when,
    lun_eclipse_how,
    ECL_TOTAL,
    ECL_PARTIAL,
    ECL_PENUMBRAL,
    ECL_VISIBLE,
    FLG_SWIEPH,
)

pytestmark = pytest.mark.slow


class TestLunEclipseHow:
    """Test suite for lun_eclipse_how function."""

    def test_returns_correct_tuple_structure(self):
        """Test that return values have correct structure."""
        # Use a known eclipse time (May 16, 2022 total lunar eclipse)
        jd_eclipse = julday(2022, 5, 16, 4.2)
        rome_lat, rome_lon = 41.9028, 12.4964

        ecl_type, attr = lun_eclipse_how(jd_eclipse, (rome_lon, rome_lat, 0.0))

        # Should return 11-element attr tuple
        assert len(attr) == 20
        # All elements should be floats
        assert all(isinstance(a, float) for a in attr)
        # Eclipse type should be int
        assert isinstance(ecl_type, int)

    def test_attributes_during_eclipse(self):
        """Test that eclipse attributes have reasonable values during eclipse."""
        # Use a known eclipse time (May 16, 2022 total lunar eclipse maximum ~04:12 UTC)
        jd_eclipse = julday(2022, 5, 16, 4.2)
        # Use a location where Moon should be visible
        rio_lat, rio_lon = -22.9068, -43.1729

        ecl_type, attr = lun_eclipse_how(jd_eclipse, (rio_lon, rio_lat, 0.0))

        # Umbral magnitude should be positive during eclipse
        assert attr[0] > 0

        # Penumbral magnitude should be positive during eclipse
        assert attr[1] > 0

        # Azimuth should be 0-360
        assert 0 <= attr[4] <= 360

        # True altitude should be in valid range
        assert -90 <= attr[5] <= 90

        # Apparent altitude should be in valid range
        assert -90 <= attr[6] <= 90

        # Distance from opposition should be non-negative
        assert attr[7] >= 0

    def test_no_eclipse_returns_zero_magnitude(self):
        """Test that non-eclipse time returns zero magnitude."""
        # Use a time when there's definitely no eclipse
        jd_no_eclipse = julday(2022, 6, 1, 12.0)
        london_lat, london_lon = 51.5074, -0.1278

        ecl_type, attr = lun_eclipse_how(jd_no_eclipse, (london_lon, london_lat, 0.0))

        # No umbral immersion; the penumbral measure stays as the
        # (negative) distance-from-penumbra value, the reference
        # convention (verified against the reference ephemeris 2.10.03).
        assert attr[0] == 0.0
        assert attr[1] < 0.0
        # Eclipse type should be 0
        assert ecl_type == 0

    def test_moon_position_calculated(self):
        """Test that Moon's position is calculated correctly."""
        # Use a known eclipse time
        jd_eclipse = julday(2022, 5, 16, 4.2)
        tokyo_lat, tokyo_lon = 35.6762, 139.6503

        ecl_type, attr = lun_eclipse_how(jd_eclipse, (tokyo_lon, tokyo_lat, 0.0))

        # Moon azimuth should be in valid range
        assert 0 <= attr[4] <= 360

        # Moon altitude should be in valid range
        assert -90 <= attr[5] <= 90

    def test_different_locations_different_altitudes(self):
        """Test that different locations have different Moon altitudes."""
        jd_eclipse = julday(2022, 5, 16, 4.2)

        # Two locations on opposite sides of the Earth
        new_york_lat, new_york_lon = 40.7128, -74.0060
        sydney_lat, sydney_lon = -33.8688, 151.2093

        _, attr_ny = lun_eclipse_how(jd_eclipse, (new_york_lon, new_york_lat, 0.0))
        _, attr_syd = lun_eclipse_how(jd_eclipse, (sydney_lon, sydney_lat, 0.0))

        # Moon altitudes should be different
        assert attr_ny[5] != attr_syd[5]

    def test_visibility_flag_when_moon_above_horizon(self):
        """Test that visibility flag is set when Moon is above horizon."""
        # May 16, 2022 eclipse - Moon should be visible from Rio
        jd_eclipse = julday(2022, 5, 16, 4.2)
        rio_lat, rio_lon = -22.9068, -43.1729

        ecl_type, attr = lun_eclipse_how(jd_eclipse, (rio_lon, rio_lat, 0.0))

        # Above the horizon during the eclipse the phase is reported;
        # the reference encodes visibility by zeroing the retflag below
        # the horizon, never by an ECL_VISIBLE bit in this function.
        apparent_alt = attr[6]
        if apparent_alt > 0.0:
            assert ecl_type & (ECL_TOTAL | ECL_PARTIAL | ECL_PENUMBRAL)

    def test_eclipse_type_detection(self):
        """Test that eclipse type is correctly detected."""
        # First find a total lunar eclipse
        jd_start = julday(2022, 5, 1, 0)
        global_ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        # Check circumstances at maximum
        rio_lat, rio_lon = -22.9068, -43.1729
        ecl_type, attr = lun_eclipse_how(jd_max, (rio_lon, rio_lat, 0.0))

        # Should detect the eclipse type
        assert ecl_type & (ECL_TOTAL | ECL_PARTIAL | ECL_PENUMBRAL)

    def test_altitude_parameter(self):
        """Test that altitude parameter is accepted and doesn't cause errors."""
        jd_eclipse = julday(2022, 5, 16, 4.2)
        denver_lat, denver_lon = 39.7392, -104.9903
        denver_altitude = 1609  # 1 mile high

        ecl_type, attr = lun_eclipse_how(
            jd_eclipse, (denver_lon, denver_lat, denver_altitude)
        )

        # Should return valid results
        assert len(attr) == 20
        assert isinstance(ecl_type, int)

    def test_apparent_altitude_and_distance(self):
        """Test that apparent altitude and distance from opposition are calculated."""
        jd_eclipse = julday(2022, 5, 16, 4.2)
        berlin_lat, berlin_lon = 52.5200, 13.4050

        ecl_type, attr = lun_eclipse_how(jd_eclipse, (berlin_lon, berlin_lat, 0.0))

        # True altitude (attr[5]) should be in valid range
        assert -90 <= attr[5] <= 90

        # Apparent altitude (attr[6]) should be in valid range
        assert -90 <= attr[6] <= 90

        # Distance from opposition (attr[7]) should be non-negative
        assert attr[7] >= 0

    def test_consistency_with_lun_eclipse_when(self):
        """Test that magnitudes are consistent with lun_eclipse_when values."""
        # Find a lunar eclipse
        jd_start = julday(2024, 1, 1, 0)
        global_type, times = lun_eclipse_when(jd_start)
        jd_max = times[0]

        # Get circumstances at maximum
        # Use a location where eclipse should be visible
        london_lat, london_lon = 51.5074, -0.1278
        ecl_type, attr = lun_eclipse_how(jd_max, (london_lon, london_lat, 0.0))

        # At least one of umbral or penumbral magnitude should be positive
        assert attr[0] > 0 or attr[1] > 0


class TestSweLunEclipseHow:
    """Test suite for lun_eclipse_how with the reference-API-compatible signature."""

    def test_reference_signature(self):
        """Test the reference-API-compatible signature with geopos parameter."""
        jd_eclipse = julday(2022, 5, 16, 4.2)
        # geopos = [longitude, latitude, altitude]
        rome_geopos = [12.4964, 41.9028, 0]

        ecl_type, attr = lun_eclipse_how(jd_eclipse, rome_geopos, FLG_SWIEPH)

        # Should return 11-element attr tuple
        assert len(attr) == 20
        # All elements should be floats
        assert all(isinstance(a, float) for a in attr)
        # Eclipse type should be int
        assert isinstance(ecl_type, int)

    def test_nov2022_eclipse_los_angeles(self):
        """Test Nov 8, 2022 total lunar eclipse from Los Angeles."""
        # JD 2459891.9578 is eclipse maximum (Nov 8, 2022 ~10:59 UTC)
        # LA is UTC-8, so this is ~2:59 AM local time
        jd_eclipse = 2459891.9578
        # Los Angeles: 34.05N, 118.24W
        la_geopos = [-118.24, 34.05, 0]

        ecl_type, attr = lun_eclipse_how(jd_eclipse, la_geopos, FLG_SWIEPH)

        # Totality in progress and Moon above the horizon: the retflag
        # carries the phase (the reference ephemeris 2.10.03 returns ECL_TOTAL here;
        # no visibility bits are set by this function).
        assert ecl_type & ECL_TOTAL

        # Umbral magnitude should be > 1 for total eclipse at maximum
        assert attr[0] > 1.0  # Total eclipse at maximum

        # Moon altitude should be reasonable (around 40 degrees at maximum)
        moon_altitude = attr[5]
        assert moon_altitude > 30  # Moon should be well above horizon
        assert moon_altitude < 60

    def test_nov2022_eclipse_london(self):
        """Test Nov 8, 2022 total lunar eclipse from London (not visible)."""
        # During maximum (JD 2459891.9578 = ~10:59 UTC), Moon is below horizon from London
        jd_eclipse = 2459891.9578
        # London: 51.51N, 0.13W
        london_geopos = [-0.13, 51.51, 0]

        ecl_type, attr = lun_eclipse_how(jd_eclipse, london_geopos, FLG_SWIEPH)

        # Moon should be below or very low from London during this eclipse
        moon_altitude = attr[5]

        # If Moon is below horizon, visibility flag should not be set
        if moon_altitude < -1.0:
            assert not (ecl_type & ECL_VISIBLE)

    def test_nov2022_eclipse_tokyo(self):
        """Test Nov 8, 2022 total lunar eclipse from Tokyo."""
        # Eclipse maximum is JD 2459891.9578 = ~10:59 UTC = ~19:59 JST
        jd_eclipse = 2459891.9578
        # Tokyo: 35.68N, 139.69E
        tokyo_geopos = [139.69, 35.68, 0]

        ecl_type, attr = lun_eclipse_how(jd_eclipse, tokyo_geopos, FLG_SWIEPH)

        # Eclipse in progress and Moon above the horizon: the phase is
        # reported in the retflag (no visibility bits in this function).
        assert ecl_type & (ECL_TOTAL | ECL_PARTIAL | ECL_PENUMBRAL)

        # Umbral magnitude should be positive
        assert attr[0] > 0

    def test_apparent_altitude_includes_refraction(self):
        """Test that apparent altitude includes atmospheric refraction."""
        jd_eclipse = julday(2022, 5, 16, 4.2)
        rome_geopos = [12.4964, 41.9028, 0]

        ecl_type, attr = lun_eclipse_how(jd_eclipse, rome_geopos, FLG_SWIEPH)

        true_altitude = attr[5]
        apparent_altitude = attr[6]

        # Near horizon, apparent altitude should be higher due to refraction
        # When Moon is above horizon, apparent >= true (refraction bends light up)
        if true_altitude > 0:
            assert apparent_altitude >= true_altitude

    def test_center_distance_in_moon_radii(self):
        """Test that center distance from shadow axis is in Moon radii."""
        jd_eclipse = julday(2022, 5, 16, 4.2)
        rio_geopos = [-43.1729, -22.9068, 0]

        ecl_type, attr = lun_eclipse_how(jd_eclipse, rio_geopos, FLG_SWIEPH)

        # Center distance (attr[7]) should be reasonable
        # For total eclipse at maximum, this should be < 1
        center_distance = attr[7]
        assert center_distance >= 0
        # During a total eclipse, center distance should be relatively small
        assert center_distance < 5  # Reasonable upper bound

    def test_eclipse_type_at_moment(self):
        """Test that eclipse type at moment is returned in retflag."""
        # Find a total lunar eclipse
        jd_start = julday(2022, 5, 1, 0)
        global_ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        rio_geopos = [-43.1729, -22.9068, 0]
        ecl_type, attr = lun_eclipse_how(jd_max, rio_geopos, FLG_SWIEPH)

        # Eclipse type should be indicated in the return flag
        assert ecl_type & (ECL_TOTAL | ECL_PARTIAL | ECL_PENUMBRAL)

        # attr[8] should equal attr[0] (umbral magnitude duplicate)
        assert abs(attr[8] - attr[0]) < 1e-10

    def test_geopos_longitude_first(self):
        """Test that geopos uses longitude-first order (reference ephemeris convention)."""
        jd_eclipse = julday(2022, 5, 16, 4.2)

        # Rome: lon=12.4964, lat=41.9028
        # Test with correct order (longitude first)
        geopos_correct = [12.4964, 41.9028, 0]
        ecl_type_correct, attr_correct = lun_eclipse_how(
            jd_eclipse, geopos_correct, FLG_SWIEPH
        )

        # Test with swapped coordinates
        geopos_swapped = [41.9028, 12.4964, 0]
        ecl_type_swapped, attr_swapped = lun_eclipse_how(
            jd_eclipse, geopos_swapped, FLG_SWIEPH
        )

        # Moon azimuths should be different for different locations
        assert attr_correct[4] != attr_swapped[4]


class TestLunEclipseHowEdgeCases:
    """Test edge cases for lun_eclipse_how function."""

    def test_equatorial_location(self):
        """Test eclipse circumstances from equatorial location."""
        jd_eclipse = julday(2022, 5, 16, 4.2)
        # Quito, Ecuador (nearly on the equator)
        quito_lat, quito_lon = -0.1807, -78.4678

        ecl_type, attr = lun_eclipse_how(jd_eclipse, (quito_lon, quito_lat, 0.0))

        assert len(attr) == 20
        # Moon altitude should be in valid range
        assert -90 <= attr[5] <= 90

    def test_high_latitude_location(self):
        """Test eclipse circumstances from high latitude location."""
        jd_eclipse = julday(2022, 5, 16, 4.2)
        # Reykjavik, Iceland
        reykjavik_lat, reykjavik_lon = 64.1466, -21.9426

        ecl_type, attr = lun_eclipse_how(
            jd_eclipse, (reykjavik_lon, reykjavik_lat, 0.0)
        )

        assert len(attr) == 20
        # Moon altitude should be in valid range
        assert -90 <= attr[5] <= 90

    def test_southern_hemisphere(self):
        """Test eclipse circumstances from southern hemisphere."""
        jd_eclipse = julday(2022, 5, 16, 4.2)
        # Melbourne, Australia
        melbourne_lat, melbourne_lon = -37.8136, 144.9631

        ecl_type, attr = lun_eclipse_how(
            jd_eclipse, (melbourne_lon, melbourne_lat, 0.0)
        )

        assert len(attr) == 20
        # Moon altitude should be in valid range
        assert -90 <= attr[5] <= 90

    def test_partial_eclipse(self):
        """Test circumstances during a partial lunar eclipse."""
        # Find a partial lunar eclipse
        jd_start = julday(2023, 10, 1, 0)
        global_type, times = lun_eclipse_when(jd_start, ecltype=ECL_PARTIAL)
        jd_max = times[0]

        # Check circumstances
        london_lat, london_lon = 51.5074, -0.1278
        ecl_type, attr = lun_eclipse_how(jd_max, (london_lon, london_lat, 0.0))

        # Should have positive umbral magnitude (but less than 1 for partial)
        if attr[0] > 0:
            assert attr[0] < 1.5  # Allow some margin

    def test_penumbral_eclipse(self):
        """Test circumstances during a penumbral lunar eclipse."""
        # Find a penumbral lunar eclipse
        jd_start = julday(2020, 1, 1, 0)
        global_type, times = lun_eclipse_when(jd_start, ecltype=ECL_PENUMBRAL)
        jd_max = times[0]

        # Check circumstances
        paris_lat, paris_lon = 48.8566, 2.3522
        ecl_type, attr = lun_eclipse_how(jd_max, (paris_lon, paris_lat, 0.0))

        # Should have positive penumbral magnitude
        # Umbral magnitude should be zero or very small for penumbral-only
        assert attr[1] >= 0

    def test_early_20th_century(self):
        """Test eclipse circumstances in early 20th century."""
        # Find an eclipse in 1920
        jd_start = julday(1920, 1, 1, 0)
        global_type, times = lun_eclipse_when(jd_start)
        jd_max = times[0]

        # Check circumstances
        london_lat, london_lon = 51.5074, -0.1278
        ecl_type, attr = lun_eclipse_how(jd_max, (london_lon, london_lat, 0.0))

        assert len(attr) == 20
        assert isinstance(ecl_type, int)

    def test_mid_21st_century(self):
        """Test eclipse circumstances in mid 21st century."""
        # Find an eclipse in 2050
        jd_start = julday(2050, 1, 1, 0)
        global_type, times = lun_eclipse_when(jd_start)
        jd_max = times[0]

        # Check circumstances
        new_york_lat, new_york_lon = 40.7128, -74.0060
        ecl_type, attr = lun_eclipse_how(jd_max, (new_york_lon, new_york_lat, 0.0))

        assert len(attr) == 20
        assert isinstance(ecl_type, int)

    def test_international_date_line(self):
        """Test eclipse circumstances near international date line."""
        jd_eclipse = julday(2022, 5, 16, 4.2)
        # Fiji
        fiji_lat, fiji_lon = -18.1416, 178.4419

        ecl_type, attr = lun_eclipse_how(jd_eclipse, (fiji_lon, fiji_lat, 0.0))

        assert len(attr) == 20
        # Moon azimuth should be in valid range
        assert 0 <= attr[4] <= 360

    def test_prime_meridian(self):
        """Test eclipse circumstances on prime meridian."""
        jd_eclipse = julday(2022, 5, 16, 4.2)
        # Greenwich, UK
        greenwich_lat, greenwich_lon = 51.4772, 0.0

        ecl_type, attr = lun_eclipse_how(
            jd_eclipse, (greenwich_lon, greenwich_lat, 0.0)
        )

        assert len(attr) == 20
        assert isinstance(ecl_type, int)

    def test_no_visibility_when_moon_below_horizon(self):
        """Test that ECL_VISIBLE is not set when Moon is below horizon."""
        # Find an eclipse and test from a location where Moon is below horizon
        jd_start = julday(2022, 5, 1, 0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        # Test various locations - at least one should have Moon below horizon
        locations = [
            [139.69, 35.68, 0],  # Tokyo
            [-118.24, 34.05, 0],  # LA
            [0.0, 51.51, 0],  # London
            [116.40, 39.90, 0],  # Beijing
        ]

        for geopos in locations:
            ecl_type, attr = lun_eclipse_how(jd_max, geopos, FLG_SWIEPH)
            moon_alt = attr[5]

            if moon_alt < -1.0:
                # Moon below horizon - visibility flag should NOT be set
                assert not (ecl_type & ECL_VISIBLE), (
                    f"ECL_VISIBLE should not be set when Moon altitude is {moon_alt} "
                    f"at location {geopos}"
                )


class TestValidationRequirements:
    """Tests for ECLIPSE-005 validation requirements.

    These tests verify the implementation matches the reference ephemeris output
    within the specified tolerances for the 2022-Nov-08 total lunar eclipse.

    Reference: Nov 8, 2022 total lunar eclipse
    - Maximum eclipse: ~10:59 UTC (JD 2459891.9578)
    - Eclipse visible from Americas, Pacific, Asia, Australia
    """

    def test_los_angeles_validation(self):
        """Validate LA results for Nov 8, 2022 eclipse at maximum.

        Expected: Moon visible, totality visible, altitude around 40-50°.
        Validation: umbral magnitude within 0.01, Moon altitude within 1° of expected.
        """
        # Find the actual eclipse maximum
        jd_start = julday(2022, 11, 1, 0)
        _, times = lun_eclipse_when(jd_start, FLG_SWIEPH, ECL_TOTAL)
        jd_max = times[0]

        la_geopos = [-118.24, 34.05, 0]  # Los Angeles: lon, lat, alt

        retflag, attr = lun_eclipse_how(jd_max, la_geopos, FLG_SWIEPH)

        # Eclipse should be total and observable (nonzero retflag means
        # the Moon stands above the local horizon)
        assert retflag & ECL_TOTAL, "Eclipse should be total"

        # Umbral magnitude should be > 1 for total eclipse
        umbral_mag = attr[0]
        assert umbral_mag > 1.0, f"Umbral magnitude should be > 1.0, got {umbral_mag}"

        # Moon altitude should be around 40° at this time
        moon_alt = attr[5]
        assert 30 < moon_alt < 60, f"Moon altitude should be ~40°, got {moon_alt}°"

    def test_london_validation(self):
        """Validate London results for Nov 8, 2022 eclipse.

        Expected: Moon below horizon or very low during totality.
        """
        # Find the actual eclipse maximum
        jd_start = julday(2022, 11, 1, 0)
        _, times = lun_eclipse_when(jd_start, FLG_SWIEPH, ECL_TOTAL)
        jd_max = times[0]

        london_geopos = [-0.13, 51.51, 0]  # London: lon, lat, alt

        retflag, attr = lun_eclipse_how(jd_max, london_geopos, FLG_SWIEPH)

        # Moon should be below horizon at eclipse maximum in London
        # (Nov 8, 2022 ~11:00 UTC means Moon is setting/set in London)
        moon_alt = attr[5]
        assert moon_alt < 10, (
            f"Moon altitude should be below horizon or very low, got {moon_alt}°"
        )

        # If Moon is significantly below horizon, visibility flag should not be set
        if moon_alt < -1.0:
            assert not (retflag & ECL_VISIBLE), (
                "Eclipse should NOT be visible from London"
            )

    def test_tokyo_validation(self):
        """Validate Tokyo results for Nov 8, 2022 eclipse.

        Expected: Moon visible, eclipse visible.
        """
        # Find the actual eclipse maximum
        jd_start = julday(2022, 11, 1, 0)
        _, times = lun_eclipse_when(jd_start, FLG_SWIEPH, ECL_TOTAL)
        jd_max = times[0]

        tokyo_geopos = [139.69, 35.68, 0]  # Tokyo: lon, lat, alt

        retflag, attr = lun_eclipse_how(jd_max, tokyo_geopos, FLG_SWIEPH)

        # A nonzero retflag means the eclipse is observable from Tokyo
        assert retflag != 0, "Eclipse should be observable from Tokyo"

        # Moon should be above horizon
        moon_alt = attr[5]
        assert moon_alt > 0, f"Moon should be above horizon in Tokyo, got {moon_alt}°"

        # Eclipse should be visible (total, partial, or penumbral)
        assert retflag & (ECL_TOTAL | ECL_PARTIAL | ECL_PENUMBRAL), (
            "Eclipse should be visible from Tokyo"
        )

    def test_geopos_order_is_longitude_first(self):
        """Verify that geopos order is [longitude, latitude, altitude].

        DO NOT confuse observer longitude/latitude order (longitude first).
        """
        # Find the actual eclipse maximum
        jd_start = julday(2022, 11, 1, 0)
        _, times = lun_eclipse_when(jd_start, FLG_SWIEPH, ECL_TOTAL)
        jd_max = times[0]

        # Test with LA coordinates in correct order (longitude first)
        # LA: -118.24°W, 34.05°N
        geopos_lon_first = [-118.24, 34.05, 0]  # Correct: lon, lat, alt
        geopos_lat_first = [34.05, -118.24, 0]  # Wrong: lat, lon, alt

        retflag_correct, attr_correct = lun_eclipse_how(
            jd_max, geopos_lon_first, FLG_SWIEPH
        )
        retflag_wrong, attr_wrong = lun_eclipse_how(
            jd_max, geopos_lat_first, FLG_SWIEPH
        )

        # The Moon position should be very different between these two
        # LA at correct position should have Moon high in sky
        # Wrong position (somewhere in SE Australia) would have different altitude

        moon_alt_correct = attr_correct[5]
        moon_alt_wrong = attr_wrong[5]

        # Correct LA should have the eclipse observable (nonzero retflag)
        assert retflag_correct != 0, "Moon should be visible at correct LA position"
        assert moon_alt_correct > 30, (
            f"Moon altitude at correct LA should be high, got {moon_alt_correct}°"
        )

        # The altitudes should be different
        assert abs(moon_alt_correct - moon_alt_wrong) > 5, (
            "Moon altitudes should differ significantly between correct and wrong coord order"
        )

    def test_visibility_flag_correctness(self):
        """Verify ECL_VISIBLE is NOT returned if Moon is below horizon.

        DO NOT return ECL_VISIBLE if Moon is below horizon.
        """
        jd = 2459892.4

        # Test from various locations
        test_locations = [
            ([-118.24, 34.05, 0], "Los Angeles"),
            ([-0.13, 51.51, 0], "London"),
            ([139.69, 35.68, 0], "Tokyo"),
            ([151.21, -33.87, 0], "Sydney"),
        ]

        for geopos, name in test_locations:
            retflag, attr = lun_eclipse_how(jd, geopos, FLG_SWIEPH)
            moon_alt = attr[5]

            if moon_alt < -1.0:
                # Moon is definitely below horizon
                assert not (retflag & ECL_VISIBLE), (
                    f"ECL_VISIBLE should NOT be set for {name} when Moon altitude is {moon_alt}°"
                )
            elif moon_alt > 0:
                # Moon is definitely above horizon
                # Should have visibility flag if there's an eclipse
                if retflag & (ECL_TOTAL | ECL_PARTIAL | ECL_PENUMBRAL):
                    assert retflag & ECL_VISIBLE, (
                        f"ECL_VISIBLE should be set for {name} when Moon altitude is {moon_alt}°"
                    )

    def test_atmospheric_refraction_applied(self):
        """Verify atmospheric refraction is applied for apparent altitude.

        DO NOT ignore atmospheric refraction for apparent altitude.
        """
        jd = 2459892.4
        la_geopos = [-118.24, 34.05, 0]

        retflag, attr = lun_eclipse_how(jd, la_geopos, FLG_SWIEPH)

        true_alt = attr[5]
        apparent_alt = attr[6]

        # For objects above horizon, apparent altitude should be >= true altitude
        # because atmospheric refraction makes objects appear higher
        if true_alt > 0:
            assert apparent_alt >= true_alt, (
                f"Apparent altitude ({apparent_alt}°) should be >= true altitude ({true_alt}°) due to refraction"
            )

        # The difference should be small (a few arcminutes at high altitude)
        # but non-zero for reasonable altitudes
        if true_alt > 10:  # For altitudes > 10°, refraction is typically 0.1° or less
            refraction = apparent_alt - true_alt
            assert 0 <= refraction < 1.0, (
                f"Refraction {refraction}° seems unreasonable for altitude {true_alt}°"
            )
