"""
Tests for lowercase 'i' house system (Sunshine alternative).

Verifies that lowercase 'i' produces identical results to uppercase 'I'
(Sunshine/Makransky house system).
"""

import pytest
import swisseph as swe
import libephemeris as ephem


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# Test locations with names for clarity
TEST_LOCATIONS = [
    (
        "johnny_depp",
        2439606.15625,
        37.7749,
        -122.4194,
    ),  # 1963-06-09 03:45:00, San Francisco
    (
        "paul_mccartney",
        2431259.4375,
        53.4084,
        -2.9916,
    ),  # 1942-06-18 14:30:00, Liverpool
    (
        "john_lennon",
        2428867.7916666665,
        53.4084,
        -2.9916,
    ),  # 1940-10-09 07:00:00, Liverpool
]

# Standard test locations
STANDARD_LOCATIONS = [
    ("Rome", 41.9028, 12.4964),
    ("New York", 40.7128, -74.0060),
    ("Sydney", -33.8688, 151.2093),
    ("Equator", 0.0, 0.0),
]

TEST_DATES = [
    (2000, 1, 1, 12.0, "J2000"),
    (2024, 6, 15, 0.0, "Summer 2024"),
    (1980, 5, 20, 14.5, "Past"),
]

TOLERANCE = 0.0001  # 0.36 arcseconds


class TestSunshineLowercase:
    """Tests for lowercase 'i' (Sunshine alternative) house system."""

    @pytest.mark.parametrize("name,jd,lat,lon", TEST_LOCATIONS)
    def test_sunshine_lowercase_vs_uppercase(self, name, jd, lat, lon):
        """Test that lowercase 'i' produces identical results to uppercase 'I'."""
        # Calculate with uppercase 'I'
        cusps_I, ascmc_I = ephem.houses(jd, lat, lon, ord("I"))

        # Calculate with lowercase 'i'
        cusps_i, ascmc_i = ephem.houses(jd, lat, lon, ord("i"))

        # Verify cusps are identical
        for house_num in range(12):
            diff = angular_diff(cusps_I[house_num], cusps_i[house_num])
            assert diff < TOLERANCE, (
                f"{name}: House {house_num + 1} differs between 'I' and 'i': "
                f"'I'={cusps_I[house_num]:.6f}°, 'i'={cusps_i[house_num]:.6f}°, diff={diff:.6f}°"
            )

        # Verify ASCMC values are identical
        for i in range(min(len(ascmc_I), len(ascmc_i))):
            diff = angular_diff(ascmc_I[i], ascmc_i[i])
            assert diff < TOLERANCE, (
                f"{name}: ASCMC[{i}] differs between 'I' and 'i': "
                f"'I'={ascmc_I[i]:.6f}°, 'i'={ascmc_i[i]:.6f}°, diff={diff:.6f}°"
            )

    @pytest.mark.parametrize("name,jd,lat,lon", TEST_LOCATIONS)
    def test_sunshine_lowercase_vs_swisseph(self, name, jd, lat, lon):
        """Test that lowercase 'i' matches Swiss Ephemeris 'I' system."""
        # Swiss Ephemeris (uppercase 'I')
        cusps_swe, ascmc_swe = swe.houses(jd, lat, lon, b"I")

        # libephemeris with lowercase 'i'
        cusps_py, ascmc_py = ephem.houses(jd, lat, lon, ord("i"))

        # Verify cusps match
        max_diff = 0.0
        for house_num in range(12):
            diff = angular_diff(cusps_swe[house_num], cusps_py[house_num])
            max_diff = max(max_diff, diff)

        # Tolerance is slightly relaxed for comparison with Swiss Ephemeris
        tolerance = 0.001  # ~3.6 arcseconds
        assert max_diff < tolerance, (
            f"{name}: Max cusp diff {max_diff:.6f}° exceeds tolerance "
            f"when comparing lowercase 'i' with Swiss Ephemeris 'I'"
        )

    @pytest.mark.parametrize("name,lat,lon", STANDARD_LOCATIONS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_sunshine_both_cases_match(
        self, name, lat, lon, year, month, day, hour, date_desc
    ):
        """Test that both 'I' and 'i' work correctly at various locations/dates."""
        jd = swe.julday(year, month, day, hour)

        # Both should produce the same results
        cusps_I, ascmc_I = ephem.houses(jd, lat, lon, ord("I"))
        cusps_i, ascmc_i = ephem.houses(jd, lat, lon, ord("i"))

        # Verify identical results
        for house_num in range(12):
            diff = angular_diff(cusps_I[house_num], cusps_i[house_num])
            assert diff < TOLERANCE, (
                f"House {house_num + 1} at {name} ({date_desc}): "
                f"'I'={cusps_I[house_num]:.6f}°, 'i'={cusps_i[house_num]:.6f}°, diff={diff:.6f}°"
            )


class TestSunshineName:
    """Tests for house_name with Sunshine systems."""

    def test_uppercase_I_name(self):
        """Test that 'I' returns 'Sunshine' name."""
        name = ephem.house_name(ord("I"))
        assert name == "Sunshine", f"Expected 'Sunshine' for 'I', got '{name}'"

    def test_lowercase_i_name(self):
        """Test that 'i' returns 'Sunshine/alt.' name."""
        name = ephem.house_name(ord("i"))
        assert name == "Sunshine/alt.", (
            f"Expected 'Sunshine/alt.' for 'i', got '{name}'"
        )


class TestSunshineNotFallingBackToPlacidus:
    """Verify lowercase 'i' doesn't incorrectly fall back to Placidus."""

    @pytest.mark.parametrize("name,jd,lat,lon", TEST_LOCATIONS)
    def test_sunshine_different_from_placidus(self, name, jd, lat, lon):
        """Verify Sunshine 'i' produces different results from Placidus (the old default)."""
        # Calculate with lowercase 'i' (should be Sunshine)
        cusps_i, _ = ephem.houses(jd, lat, lon, ord("i"))

        # Calculate with Placidus
        cusps_P, _ = ephem.houses(jd, lat, lon, ord("P"))

        # They should be DIFFERENT (if they're the same, lowercase was falling back to Placidus)
        total_diff = sum(angular_diff(cusps_i[i], cusps_P[i]) for i in range(12))

        # Sunshine and Placidus should differ significantly (usually 6-10+ degrees total)
        assert total_diff > 1.0, (
            f"{name}: Sunshine 'i' is too similar to Placidus! "
            f"Total difference is only {total_diff:.4f}°. "
            f"This suggests 'i' is incorrectly falling back to Placidus."
        )

    @pytest.mark.parametrize("name,jd,lat,lon", TEST_LOCATIONS)
    def test_sunshine_matches_uppercase(self, name, jd, lat, lon):
        """Verify lowercase 'i' matches uppercase 'I' exactly."""
        cusps_I, _ = ephem.houses(jd, lat, lon, ord("I"))
        cusps_i, _ = ephem.houses(jd, lat, lon, ord("i"))

        # They should be IDENTICAL
        total_diff = sum(angular_diff(cusps_I[i], cusps_i[i]) for i in range(12))

        assert total_diff < 0.001, (
            f"{name}: Sunshine 'I' and 'i' differ by {total_diff:.6f}° total! "
            f"They should be identical."
        )
