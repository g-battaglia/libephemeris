"""
Tests for the azalt_rev coordinate transformation function.

Tests verify that azalt_rev correctly transforms horizontal (azimuth/altitude)
coordinates to equatorial/ecliptic coordinates, matching pyswisseph's
swe.azalt_rev() behavior.

pyswisseph signature: azalt_rev(tjdut, flag, geopos, azimuth, true_altitude)
  - geopos: (longitude, latitude, altitude_meters)
  - Returns: (x1, x2) depending on flag
"""

import pytest
import swisseph as swe
import libephemeris as ephem


class TestAzaltRevBasic:
    """Basic functionality tests for azalt_rev."""

    def test_azalt_rev_exported(self):
        """Test that azalt_rev is exported from the package."""
        assert hasattr(ephem, "azalt_rev")
        assert callable(ephem.azalt_rev)

    def test_azalt_rev_constants_exported(self):
        """Test that azalt_rev constants are exported."""
        assert hasattr(ephem, "HOR2ECL")
        assert hasattr(ephem, "HOR2EQU")
        assert ephem.HOR2ECL == 0
        assert ephem.HOR2EQU == 1

    def test_azalt_rev_returns_tuple(self):
        """Test that azalt_rev returns a tuple of two elements."""
        jd = ephem.julday(2024, 6, 15, 12.0)
        geopos = (12.5, 41.9, 0)  # lon, lat, alt
        result = ephem.azalt_rev(jd, ephem.HOR2EQU, geopos, 90.0, 45.0)

        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_azalt_rev_ra_in_range(self):
        """Test that RA is in [0, 360) range."""
        jd = ephem.julday(2024, 6, 15, 12.0)
        geopos = (12.5, 41.9, 0)
        azimuths = [0.0, 90.0, 180.0, 270.0]

        for az in azimuths:
            ra, dec = ephem.azalt_rev(jd, ephem.HOR2EQU, geopos, az, 30.0)
            assert 0.0 <= ra < 360.0, f"RA {ra} out of range for azimuth {az}"

    def test_azalt_rev_dec_in_range(self):
        """Test that declination is in valid range (-90, 90)."""
        jd = ephem.julday(2024, 6, 15, 12.0)
        geopos = (12.5, 41.9, 0)
        ra, dec = ephem.azalt_rev(jd, ephem.HOR2EQU, geopos, 90.0, 45.0)

        assert -90.0 <= dec <= 90.0

    def test_azalt_rev_ecl_lon_in_range(self):
        """Test that ecliptic longitude is in [0, 360) range."""
        jd = ephem.julday(2024, 6, 15, 12.0)
        geopos = (12.5, 41.9, 0)
        ecl_lon, ecl_lat = ephem.azalt_rev(jd, ephem.HOR2ECL, geopos, 90.0, 45.0)

        assert 0.0 <= ecl_lon < 360.0


class TestAzaltRevVsSwisseph:
    """Comparison tests with pyswisseph's swe.azalt_rev()."""

    @pytest.mark.parametrize(
        "azimuth,altitude",
        [
            (0.0, 30.0),  # South
            (90.0, 45.0),  # West
            (180.0, 60.0),  # North
            (270.0, 30.0),  # East
            (45.0, 20.0),
            (135.0, 50.0),
            (225.0, 70.0),
            (315.0, 15.0),
        ],
    )
    def test_azalt_rev_equatorial_matches_swisseph(self, azimuth, altitude):
        """Test horizontal to equatorial transformation matches pyswisseph."""
        jd = ephem.julday(2024, 6, 15, 12.0)
        geopos = (12.5, 41.9, 100.0)  # lon, lat, alt_m

        result_lib = ephem.azalt_rev(jd, ephem.HOR2EQU, geopos, azimuth, altitude)
        result_swe = swe.azalt_rev(jd, swe.HOR2EQU, geopos, azimuth, altitude)

        # Compare RA (handle wrap-around)
        ra_diff = abs(result_lib[0] - result_swe[0])
        if ra_diff > 180:
            ra_diff = 360 - ra_diff
        assert ra_diff < 0.1, (
            f"RA mismatch: lib={result_lib[0]:.4f}, swe={result_swe[0]:.4f}"
        )

        # Compare Declination
        assert abs(result_lib[1] - result_swe[1]) < 0.1, (
            f"Dec mismatch: lib={result_lib[1]:.4f}, swe={result_swe[1]:.4f}"
        )

    @pytest.mark.parametrize(
        "azimuth,altitude",
        [
            (0.0, 30.0),
            (90.0, 45.0),
            (180.0, 60.0),
            (45.0, 20.0),
        ],
    )
    def test_azalt_rev_ecliptic_matches_swisseph(self, azimuth, altitude):
        """Test horizontal to ecliptic transformation matches pyswisseph."""
        jd = ephem.julday(2024, 6, 15, 12.0)
        geopos = (12.5, 41.9, 0.0)  # lon, lat, alt_m

        result_lib = ephem.azalt_rev(jd, ephem.HOR2ECL, geopos, azimuth, altitude)
        result_swe = swe.azalt_rev(jd, swe.HOR2ECL, geopos, azimuth, altitude)

        # Compare ecliptic longitude (handle wrap-around)
        lon_diff = abs(result_lib[0] - result_swe[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff
        assert lon_diff < 0.5, (
            f"Ecliptic lon mismatch: lib={result_lib[0]:.4f}, swe={result_swe[0]:.4f}"
        )

        # Compare ecliptic latitude
        assert abs(result_lib[1] - result_swe[1]) < 0.5, (
            f"Ecliptic lat mismatch: lib={result_lib[1]:.4f}, swe={result_swe[1]:.4f}"
        )


class TestAzaltRevRoundTrip:
    """Test that azalt and azalt_rev are consistent inverses."""

    @pytest.mark.parametrize(
        "ra,dec",
        [
            (0.0, 0.0),
            (90.0, 23.5),
            (180.0, -15.0),
            (270.0, 45.0),
            (45.0, 30.0),
            (120.0, -30.0),
        ],
    )
    def test_azalt_azalt_rev_roundtrip(self, ra, dec):
        """Test roundtrip: equatorial -> horizontal -> equatorial."""
        jd = ephem.julday(2024, 6, 15, 12.0)
        geopos = (12.5, 41.9, 0.0)  # lon, lat, alt

        # Convert equatorial to horizontal
        coord = (ra, dec, 1.0)
        az, alt_true, _ = ephem.azalt(jd, ephem.EQU2HOR, geopos, 0.0, 15.0, coord)

        # Convert back to equatorial
        ra_back, dec_back = ephem.azalt_rev(jd, ephem.HOR2EQU, geopos, az, alt_true)

        # Compare RA (handle wrap-around)
        ra_diff = abs(ra - ra_back)
        if ra_diff > 180:
            ra_diff = 360 - ra_diff
        assert ra_diff < 0.01, (
            f"RA roundtrip error: original={ra:.4f}, recovered={ra_back:.4f}"
        )

        # Compare Declination
        assert abs(dec - dec_back) < 0.01, (
            f"Dec roundtrip error: original={dec:.4f}, recovered={dec_back:.4f}"
        )


class TestAzaltRevDifferentLocations:
    """Test azalt_rev at different geographic locations."""

    @pytest.mark.parametrize(
        "lat,lon",
        [
            (0.0, 0.0),  # Equator at prime meridian
            (51.5, -0.12),  # London
            (40.7, -74.0),  # New York
            (-33.9, 151.2),  # Sydney
            (35.7, 139.7),  # Tokyo
            (64.1, -21.9),  # Reykjavik (high latitude)
        ],
    )
    def test_azalt_rev_various_locations(self, lat, lon):
        """Test azalt_rev at various geographic locations."""
        jd = ephem.julday(2024, 6, 21, 12.0)  # Summer solstice noon
        geopos = (lon, lat, 0)  # Note: pyswisseph order is (lon, lat, alt)
        azimuth, altitude = 90.0, 45.0

        ra, dec = ephem.azalt_rev(jd, ephem.HOR2EQU, geopos, azimuth, altitude)

        # Basic sanity checks
        assert 0.0 <= ra < 360.0
        assert -90.0 <= dec <= 90.0


class TestAzaltRevEdgeCases:
    """Test azalt_rev with edge cases."""

    def test_azalt_rev_zenith(self):
        """Test object at zenith (altitude = 90)."""
        jd = ephem.julday(2024, 6, 15, 12.0)
        geopos = (12.5, 41.9, 0)

        ra, dec = ephem.azalt_rev(jd, ephem.HOR2EQU, geopos, 0.0, 90.0)

        # At zenith, declination should approximately equal latitude
        assert abs(dec - 41.9) < 0.1, (
            f"At zenith, Dec ({dec:.4f}) should equal lat ({41.9:.4f})"
        )

    def test_azalt_rev_horizon(self):
        """Test object at horizon (altitude = 0)."""
        jd = ephem.julday(2024, 6, 15, 12.0)
        geopos = (12.5, 41.9, 0)

        ra, dec = ephem.azalt_rev(jd, ephem.HOR2EQU, geopos, 180.0, 0.0)

        # Should return valid coordinates
        assert 0.0 <= ra < 360.0
        assert -90.0 <= dec <= 90.0

    def test_azalt_rev_below_horizon(self):
        """Test object below horizon (negative altitude)."""
        jd = ephem.julday(2024, 6, 15, 12.0)
        geopos = (12.5, 41.9, 0)

        # Negative altitude (below horizon)
        ra, dec = ephem.azalt_rev(jd, ephem.HOR2EQU, geopos, 90.0, -10.0)

        # Should still return valid coordinates
        assert 0.0 <= ra < 360.0
        assert -90.0 <= dec <= 90.0


class TestAzaltRevDifferentTimes:
    """Test azalt_rev at different times."""

    @pytest.mark.parametrize(
        "hour",
        [0.0, 6.0, 12.0, 18.0, 23.99],
    )
    def test_azalt_rev_different_hours(self, hour):
        """Test azalt_rev at different hours of the day."""
        jd = ephem.julday(2024, 6, 15, hour)
        geopos = (12.5, 41.9, 0)
        azimuth, altitude = 90.0, 30.0

        ra, dec = ephem.azalt_rev(jd, ephem.HOR2EQU, geopos, azimuth, altitude)

        # RA should change with time (due to Earth's rotation)
        # Just verify valid output
        assert 0.0 <= ra < 360.0
        assert -90.0 <= dec <= 90.0
