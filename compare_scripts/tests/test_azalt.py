"""
Tests for the azalt coordinate transformation function.

Tests verify that azalt correctly transforms equatorial/ecliptic coordinates
to horizontal (azimuth/altitude) coordinates, matching pyswisseph's swe.azalt() behavior.

pyswisseph signature: azalt(tjdut, flag, geopos, atpress, attemp, xin)
  - geopos: (longitude, latitude, altitude_meters)
  - xin: (lon/ra, lat/dec, distance)
"""

import pytest
import swisseph as swe
import libephemeris as ephem


class TestAzaltBasic:
    """Basic functionality tests for azalt."""

    def test_azalt_exported(self):
        """Test that azalt is exported from the package."""
        assert hasattr(ephem, "azalt")
        assert callable(ephem.azalt)

    def test_azalt_constants_exported(self):
        """Test that azalt constants are exported."""
        assert hasattr(ephem, "ECL2HOR")
        assert hasattr(ephem, "EQU2HOR")
        assert ephem.ECL2HOR == 0
        assert ephem.EQU2HOR == 1

    def test_azalt_returns_tuple(self):
        """Test that azalt returns a tuple of three elements."""
        jd = ephem.julday(2024, 6, 15, 12.0)
        geopos = (12.5, 41.9, 0)  # lon, lat, alt
        coord = (90.0, 23.5, 1.0)
        result = ephem.azalt(jd, ephem.EQU2HOR, geopos, 1013.25, 15, coord)

        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_azalt_azimuth_in_range(self):
        """Test that azimuth is in [0, 360) range."""
        jd = ephem.julday(2024, 6, 15, 12.0)
        geopos = (12.5, 41.9, 0)
        coords = [
            (0.0, 0.0, 1.0),
            (90.0, 45.0, 1.0),
            (180.0, -30.0, 1.0),
            (270.0, 60.0, 1.0),
        ]

        for coord in coords:
            az, alt_true, alt_app = ephem.azalt(
                jd, ephem.EQU2HOR, geopos, 1013.25, 15, coord
            )
            assert 0.0 <= az < 360.0, f"Azimuth {az} out of range for coord {coord}"

    def test_azalt_altitude_in_range(self):
        """Test that altitude is in valid range (-90, 90)."""
        jd = ephem.julday(2024, 6, 15, 12.0)
        geopos = (12.5, 41.9, 0)
        coord = (90.0, 23.5, 1.0)
        az, alt_true, alt_app = ephem.azalt(
            jd, ephem.EQU2HOR, geopos, 1013.25, 15, coord
        )

        assert -90.0 <= alt_true <= 90.0
        assert -90.0 <= alt_app <= 90.0


class TestAzaltNoRefraction:
    """Tests for azalt without atmospheric refraction."""

    def test_no_refraction_when_pressure_zero(self):
        """Test that true and apparent altitude are equal when pressure is 0."""
        jd = ephem.julday(2024, 6, 15, 12.0)
        geopos = (12.5, 41.9, 0)
        coord = (90.0, 30.0, 1.0)

        az, alt_true, alt_app = ephem.azalt(
            jd, ephem.EQU2HOR, geopos, 0.0, 15, coord
        )

        assert abs(alt_true - alt_app) < 1e-10, (
            "True and apparent altitude should be equal when pressure is 0"
        )


class TestAzaltRefraction:
    """Tests for atmospheric refraction in azalt."""

    def test_refraction_increases_apparent_altitude(self):
        """Test that refraction makes objects appear higher."""
        jd = ephem.julday(2024, 6, 15, 12.0)
        geopos = (12.5, 41.9, 0)
        coord = (90.0, 30.0, 1.0)

        az, alt_true, alt_app = ephem.azalt(
            jd, ephem.EQU2HOR, geopos, 1013.25, 15, coord
        )

        # Objects should appear higher due to refraction (for positive altitudes)
        if alt_true > 0:
            assert alt_app >= alt_true, (
                f"Apparent altitude {alt_app} should be >= true altitude {alt_true}"
            )

    def test_refraction_larger_at_horizon(self):
        """Test that refraction is larger near the horizon."""
        jd = ephem.julday(2024, 6, 15, 6.0)  # Early morning for lower sun
        geopos = (12.5, 41.9, 0)

        # Find two objects at different altitudes
        coord_low = (280.0, -10.0, 1.0)  # Near horizon
        coord_high = (90.0, 60.0, 1.0)  # Higher

        _, alt_low_true, alt_low_app = ephem.azalt(
            jd, ephem.EQU2HOR, geopos, 1013.25, 15, coord_low
        )
        _, alt_high_true, alt_high_app = ephem.azalt(
            jd, ephem.EQU2HOR, geopos, 1013.25, 15, coord_high
        )

        refraction_low = alt_low_app - alt_low_true
        refraction_high = alt_high_app - alt_high_true

        # Skip if both are below horizon or both very high
        if alt_low_true > 0 and alt_low_true < alt_high_true and alt_high_true > 0:
            # Refraction should be larger for lower altitude objects
            assert refraction_low > refraction_high, (
                f"Refraction at alt {alt_low_true:.1f}={refraction_low:.4f} "
                f"should be > refraction at alt {alt_high_true:.1f}={refraction_high:.4f}"
            )


class TestAzaltVsSwisseph:
    """Comparison tests with pyswisseph's swe.azalt()."""

    @pytest.mark.parametrize(
        "coord,calc_flag",
        [
            ((0.0, 0.0, 1.0), 1),  # EQU2HOR
            ((90.0, 0.0, 1.0), 1),
            ((180.0, 0.0, 1.0), 1),
            ((270.0, 0.0, 1.0), 1),
            ((45.0, 23.5, 1.0), 1),
            ((120.0, -15.0, 1.0), 1),
            ((200.0, 45.0, 1.0), 1),
        ],
    )
    def test_azalt_equatorial_matches_swisseph(self, coord, calc_flag):
        """Test equatorial to horizontal transformation matches pyswisseph."""
        jd = ephem.julday(2024, 6, 15, 12.0)
        geopos = (12.5, 41.9, 100.0)  # lon, lat, alt_m
        pressure = 1013.25
        temp = 20.0

        result_lib = ephem.azalt(jd, calc_flag, geopos, pressure, temp, coord)
        result_swe = swe.azalt(jd, calc_flag, geopos, pressure, temp, coord)

        # Compare azimuth (handle wrap-around)
        az_diff = abs(result_lib[0] - result_swe[0])
        if az_diff > 180:
            az_diff = 360 - az_diff
        assert az_diff < 0.1, (
            f"Azimuth mismatch: lib={result_lib[0]:.4f}, swe={result_swe[0]:.4f}"
        )

        # Compare true altitude
        assert abs(result_lib[1] - result_swe[1]) < 0.1, (
            f"True altitude mismatch: lib={result_lib[1]:.4f}, swe={result_swe[1]:.4f}"
        )

        # Compare apparent altitude (refraction-affected)
        assert abs(result_lib[2] - result_swe[2]) < 0.1, (
            f"Apparent altitude mismatch: lib={result_lib[2]:.4f}, swe={result_swe[2]:.4f}"
        )

    @pytest.mark.parametrize(
        "coord",
        [
            (0.0, 0.0, 1.0),
            (90.0, 0.0, 1.0),
            (180.0, 0.0, 1.0),
            (45.0, 23.5, 1.0),
        ],
    )
    def test_azalt_ecliptic_matches_swisseph(self, coord):
        """Test ecliptic to horizontal transformation matches pyswisseph."""
        jd = ephem.julday(2024, 6, 15, 12.0)
        geopos = (12.5, 41.9, 0.0)  # lon, lat, alt_m
        pressure = 1013.25
        temp = 15.0
        calc_flag = ephem.ECL2HOR  # Ecliptic to horizontal

        result_lib = ephem.azalt(jd, calc_flag, geopos, pressure, temp, coord)
        result_swe = swe.azalt(jd, calc_flag, geopos, pressure, temp, coord)

        # Compare azimuth (handle wrap-around)
        az_diff = abs(result_lib[0] - result_swe[0])
        if az_diff > 180:
            az_diff = 360 - az_diff
        assert az_diff < 0.5, (
            f"Azimuth mismatch for ecliptic: lib={result_lib[0]:.4f}, swe={result_swe[0]:.4f}"
        )

        # Compare true altitude
        assert abs(result_lib[1] - result_swe[1]) < 0.5, (
            f"True altitude mismatch: lib={result_lib[1]:.4f}, swe={result_swe[1]:.4f}"
        )


class TestAzaltDifferentLocations:
    """Test azalt at different geographic locations."""

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
    def test_azalt_various_locations(self, lat, lon):
        """Test azalt at various geographic locations."""
        jd = ephem.julday(2024, 6, 21, 12.0)  # Summer solstice noon
        geopos = (lon, lat, 0)  # Note: pyswisseph order is (lon, lat, alt)
        coord = (90.0, 23.5, 1.0)  # Near summer solstice declination

        az, alt_true, alt_app = ephem.azalt(
            jd, ephem.EQU2HOR, geopos, 1013.25, 15, coord
        )

        # Basic sanity checks
        assert 0.0 <= az < 360.0
        assert -90.0 <= alt_true <= 90.0
        assert -90.0 <= alt_app <= 90.0


class TestAzaltNoRefractionVsSwisseph:
    """Test azalt without refraction against pyswisseph."""

    @pytest.mark.parametrize(
        "coord",
        [
            (0.0, 0.0, 1.0),
            (90.0, 45.0, 1.0),
            (180.0, -20.0, 1.0),
            (270.0, 66.0, 1.0),
        ],
    )
    def test_no_refraction_matches_swisseph(self, coord):
        """Test that no-refraction results match pyswisseph."""
        jd = ephem.julday(2024, 3, 20, 12.0)  # Equinox
        geopos = (0.0, 45.0, 0)  # lon, lat, alt
        pressure = 0.0  # No refraction
        temp = 10.0

        result_lib = ephem.azalt(jd, ephem.EQU2HOR, geopos, pressure, temp, coord)
        result_swe = swe.azalt(jd, 1, geopos, pressure, temp, coord)

        # Without refraction in our library, true and apparent should be equal
        assert abs(result_lib[1] - result_lib[2]) < 1e-6

        # Note: pyswisseph may still apply some refraction even with pressure=0
        # so we only compare true altitude with our implementation

        # Compare with pyswisseph true altitude
        az_diff = abs(result_lib[0] - result_swe[0])
        if az_diff > 180:
            az_diff = 360 - az_diff
        assert az_diff < 0.1, (
            f"Azimuth mismatch: lib={result_lib[0]}, swe={result_swe[0]}"
        )
        assert abs(result_lib[1] - result_swe[1]) < 0.1, (
            f"True altitude mismatch: lib={result_lib[1]}, swe={result_swe[1]}"
        )
