"""
Tests for heliacal_ut function in libephemeris.

Tests the calculation of heliacal rising and setting events for celestial bodies.

Heliacal events are the first/last visibility of a celestial body at dawn/dusk.
These were crucial for ancient calendars (e.g., Egyptian calendar based on Sirius).

Reference data based on astronomical calculations and historical records.
"""

import pytest

pytestmark = pytest.mark.slow

from libephemeris import (
    julday,
    revjul,
    heliacal_ut,
    heliacal_pheno_ut,
    Error,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    HELIACAL_RISING,
    HELIACAL_SETTING,
    EVENING_FIRST,
    MORNING_LAST,
)


class TestHeliacalBasic:
    """Basic tests for heliacal_ut function."""

    def test_venus_heliacal_rising_returns_valid_result(self):
        """Test that Venus heliacal rising returns a valid Julian Day."""
        # Start from January 1, 2024
        jd_start = julday(2024, 1, 1, 0)
        # Rome, Italy
        geopos = (12.4964, 41.9028, 0.0)

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Venus",
            HELIACAL_RISING,
        )

        # Venus synodic period is ~584 days, so the next heliacal rising
        # after inferior conjunction may be up to ~500 days away
        assert result[0] > jd_start
        assert result[0] < jd_start + 600

    def test_mercury_heliacal_rising_returns_valid_result(self):
        """Test that Mercury heliacal rising returns a valid Julian Day."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Mercury",
            HELIACAL_RISING,
        )

        # Mercury has ~3 synodic periods per year, so should find one quickly
        assert result[0] > jd_start
        assert result[0] < jd_start + 120  # Within ~4 months

    def test_jupiter_heliacal_rising(self):
        """Test Jupiter heliacal rising."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (-0.1278, 51.5074, 0.0)  # London

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Jupiter",
            HELIACAL_RISING,
        )

        assert result[0] > jd_start
        assert result[0] < jd_start + 400  # Within a bit more than a year

    def test_saturn_heliacal_rising(self):
        """Test Saturn heliacal rising."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (-74.0060, 40.7128, 0.0)  # New York

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Saturn",
            HELIACAL_RISING,
        )

        # Saturn may not have heliacal rising within search window for all start dates
        # Just verify it returns a valid response
        if result[0] > 0:
            assert result[0] > jd_start
            assert result[0] < jd_start + 400

    def test_mars_heliacal_rising(self):
        """Test Mars heliacal rising."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (139.6503, 35.6762, 0.0)  # Tokyo

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Mars",
            HELIACAL_RISING,
        )

        # Mars has ~780 day synodic period, may not have event in window
        if result[0] > 0:
            assert result[0] > jd_start
            assert result[0] < jd_start + 780


class TestHeliacalEventTypes:
    """Test different heliacal event types."""

    def test_heliacal_setting(self):
        """Test heliacal setting (evening last visibility)."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Venus",
            HELIACAL_SETTING,
        )

        assert result[0] > jd_start

    def test_evening_first(self):
        """Test evening first visibility (after superior conjunction)."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Venus",
            EVENING_FIRST,
        )

        # May or may not find this event depending on current position
        # Just verify it returns something valid
        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_morning_last(self):
        """Test morning last visibility (before superior conjunction)."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Mercury",
            MORNING_LAST,
        )

        # May or may not find this event
        assert isinstance(result, tuple)
        assert len(result) == 3


class TestHeliacalValidation:
    """Test input validation for heliacal_ut."""

    def test_sun_raises_error(self):
        """Test that Sun raises ValueError."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)

        with pytest.raises(Error, match="Sun"):
            heliacal_ut(
                jd_start,
                geopos,
                (1013.25, 15.0, 50.0, 0.0),
                (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                "Sun",
                HELIACAL_RISING,
            )

    def test_moon_raises_error(self):
        """Test that Moon raises ValueError."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)

        with pytest.raises(Error, match="Moon"):
            heliacal_ut(
                jd_start,
                geopos,
                (1013.25, 15.0, 50.0, 0.0),
                (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                "Moon",
                HELIACAL_RISING,
            )

    def test_invalid_event_type_raises_error(self):
        """Test that invalid event type raises ValueError."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)

        with pytest.raises(Error, match="Invalid event_type"):
            heliacal_ut(
                jd_start,
                geopos,
                (1013.25, 15.0, 50.0, 0.0),
                (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                "Venus",
                99,
            )


class TestHeliacalLocations:
    """Test heliacal calculations at various locations."""

    def test_northern_hemisphere(self):
        """Test at northern latitude."""
        jd_start = julday(2024, 6, 1, 0)  # Summer
        geopos = (24.9384, 60.1699, 0.0)  # Helsinki

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Jupiter",
            HELIACAL_RISING,
        )

        # High latitude summer nights are short, may affect visibility
        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_southern_hemisphere(self):
        """Test at southern latitude."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (151.2093, -33.8688, 0.0)  # Sydney

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Venus",
            HELIACAL_RISING,
        )

        # In the southern hemisphere, twilight times differ
        # Allow for either success or not-found
        if result[0] > 0:
            assert result[0] > jd_start

    def test_equatorial_location(self):
        """Test at equatorial location."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (32.5825, 0.3476, 0.0)  # Kampala, Uganda

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Saturn",
            HELIACAL_RISING,
        )

        # Equatorial locations have consistent twilight times
        if result[0] > 0:
            assert result[0] > jd_start


class TestHeliacalAtmosphericConditions:
    """Test heliacal calculations with different atmospheric conditions."""

    def test_with_altitude(self):
        """Test calculation with observer at altitude."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 2000.0)  # Rome, 2000 meters

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Venus",
            HELIACAL_RISING,
        )

        assert result[0] > jd_start

    def test_with_pressure(self):
        """Test calculation with custom pressure."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        result = heliacal_ut(
            jd_start,
            geopos,
            (1000.0, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Venus",
            HELIACAL_RISING,
        )

        assert result[0] > jd_start

    def test_with_humidity(self):
        """Test calculation with different humidity."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 80.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Venus",
            HELIACAL_RISING,
        )

        assert result[0] > jd_start


class TestSweHeliacalUt:
    """Test heliacal_ut reference-API-compatible API."""

    def test_swe_heliacal_ut_basic_call(self):
        """Test basic heliacal_ut call with array parameters."""
        jd_start = julday(2024, 1, 1, 0)
        # Geographic position: Rome (lon, lat, altitude)
        geopos = (12.4964, 41.9028, 0.0)
        # Atmospheric conditions: pressure, temp, humidity, extinction
        datm = (1013.25, 15.0, 40.0, 0.0)
        # Observer: age, Snellen ratio, and optical params
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        result = heliacal_ut(jd_start, geopos, datm, dobs, "Venus", HELIACAL_RISING)

        # Should return a tuple of 3 floats (jd1, jd2, jd3)
        assert isinstance(result, tuple)
        assert len(result) == 3
        # Venus synodic period is ~584 days, so next heliacal rising
        # after inferior conjunction may be up to ~500 days away
        if result[0] > 0:
            assert result[0] > jd_start
            assert result[0] < jd_start + 600

    def test_swe_heliacal_ut_with_planet_name_mercury(self):
        """Test heliacal_ut with Mercury."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        result = heliacal_ut(jd_start, geopos, datm, dobs, "Mercury", HELIACAL_RISING)

        assert isinstance(result, tuple)
        assert len(result) == 3
        # Mercury has ~3 synodic periods per year
        if result[0] > 0:
            assert result[0] > jd_start
            assert result[0] < jd_start + 120

    def test_swe_heliacal_ut_with_planet_name_mars(self):
        """Test heliacal_ut with Mars."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        result = heliacal_ut(jd_start, geopos, datm, dobs, "Mars", HELIACAL_RISING)

        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_swe_heliacal_ut_with_planet_name_jupiter(self):
        """Test heliacal_ut with Jupiter."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (-0.1278, 51.5074, 0.0)  # London
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        result = heliacal_ut(jd_start, geopos, datm, dobs, "Jupiter", HELIACAL_RISING)

        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_swe_heliacal_ut_with_planet_name_saturn(self):
        """Test heliacal_ut with Saturn."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        result = heliacal_ut(jd_start, geopos, datm, dobs, "Saturn", HELIACAL_RISING)

        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_swe_heliacal_ut_heliacal_setting(self):
        """Test heliacal_ut for heliacal setting event."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        result = heliacal_ut(jd_start, geopos, datm, dobs, "Venus", HELIACAL_SETTING)

        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_swe_heliacal_ut_evening_first(self):
        """Test heliacal_ut for evening first event."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        result = heliacal_ut(jd_start, geopos, datm, dobs, "Venus", EVENING_FIRST)

        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_swe_heliacal_ut_morning_last(self):
        """Test heliacal_ut for morning last event."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        result = heliacal_ut(jd_start, geopos, datm, dobs, "Mercury", MORNING_LAST)

        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_swe_heliacal_ut_sun_raises_error(self):
        """Test that Sun raises ValueError."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        with pytest.raises(Error, match="Sun"):
            heliacal_ut(jd_start, geopos, datm, dobs, "Sun", HELIACAL_RISING)

    def test_swe_heliacal_ut_moon_raises_error(self):
        """Test that Moon raises ValueError."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        with pytest.raises(Error, match="Moon"):
            heliacal_ut(jd_start, geopos, datm, dobs, "Moon", HELIACAL_RISING)

    def test_swe_heliacal_ut_invalid_event_type(self):
        """Test that invalid event type raises ValueError."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        with pytest.raises(Error, match="Invalid event_type"):
            heliacal_ut(jd_start, geopos, datm, dobs, "Venus", 99)

    def test_swe_heliacal_ut_invalid_object_name(self):
        """Test that unknown object name raises ValueError."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        with pytest.raises(Error, match="not recognized"):
            heliacal_ut(jd_start, geopos, datm, dobs, "InvalidPlanet", HELIACAL_RISING)

    def test_swe_heliacal_ut_case_insensitive(self):
        """Test that planet name matching is case insensitive."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        # Should work with lowercase
        result1 = heliacal_ut(jd_start, geopos, datm, dobs, "venus", HELIACAL_RISING)
        # Should work with uppercase
        result2 = heliacal_ut(jd_start, geopos, datm, dobs, "VENUS", HELIACAL_RISING)
        # Should work with mixed case
        result3 = heliacal_ut(jd_start, geopos, datm, dobs, "Venus", HELIACAL_RISING)

        # All should return same result
        assert result1[0] == result2[0] == result3[0]

    def test_swe_heliacal_ut_default_atmospheric(self):
        """Test that zero atmospheric values get defaults."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        # All zeros should get defaults
        datm = (0, 0, 0, 0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        result = heliacal_ut(jd_start, geopos, datm, dobs, "Venus", HELIACAL_RISING)

        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_swe_heliacal_ut_with_altitude(self):
        """Test heliacal_ut with observer at altitude."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 2000.0)  # 2000m altitude
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        result = heliacal_ut(jd_start, geopos, datm, dobs, "Venus", HELIACAL_RISING)

        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_swe_heliacal_ut_southern_hemisphere(self):
        """Test heliacal_ut in southern hemisphere."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (151.2093, -33.8688, 0.0)  # Sydney
        datm = (1013.25, 25.0, 60.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        result = heliacal_ut(jd_start, geopos, datm, dobs, "Venus", HELIACAL_RISING)

        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_swe_heliacal_ut_with_integer_planet_id(self):
        """Test heliacal_ut with an integer planet ID."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        # Venus is planet ID 3
        result = heliacal_ut(jd_start, geopos, datm, dobs, 3, HELIACAL_RISING)

        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_swe_heliacal_ut_numeric_string_not_planet_id(self):
        """Numeric object names are not planet IDs on this API surface."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        with pytest.raises(Error):
            heliacal_ut(jd_start, geopos, datm, dobs, "3", HELIACAL_RISING)

    def test_swe_heliacal_ut_visibility_window_ordering(self):
        """Detailed heliacal windows keep start <= optimum <= end."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        jd1, jd2, jd3 = heliacal_ut(
            jd_start, geopos, datm, dobs, "Venus", HELIACAL_RISING
        )

        assert jd1 > jd_start
        assert jd1 <= jd2 <= jd3
        assert jd3 > jd1


class TestHeliacalDateValidation:
    """Test that heliacal events return sensible dates."""

    def test_event_date_is_reasonable(self):
        """Test that returned date is within expected range."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Venus",
            HELIACAL_RISING,
        )

        if result[0] > 0:
            year, month, day, hour = revjul(result[0])
            # Event should be in 2024 or 2025
            assert year in (2024, 2025)
            # Month and day should be valid
            assert 1 <= month <= 12
            assert 1 <= day <= 31


# =============================================================================
# HELIACAL_PHENO_UT TESTS
# =============================================================================


class TestHeliacalPhenoBasic:
    """Basic tests for heliacal_pheno_ut function."""

    def test_venus_pheno_returns_valid_result(self):
        """Test that Venus heliacal phenomena returns a valid result tuple."""
        jd = julday(2024, 1, 1, 12)  # Noon
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        dret = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Venus",
            HELIACAL_RISING,
        )

        # Should return a tuple of 50 floats
        assert isinstance(dret, tuple)
        assert len(dret) == 50
        assert all(isinstance(x, float) for x in dret)

    def test_mercury_pheno_returns_valid_result(self):
        """Test that Mercury heliacal phenomena returns a valid result."""
        jd = julday(2024, 3, 15, 5)  # Morning
        geopos = (-0.1278, 51.5074, 0.0)  # London

        dret = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Mercury",
            HELIACAL_RISING,
        )

        assert isinstance(dret, tuple)
        assert len(dret) == 50

    def test_jupiter_pheno_returns_valid_result(self):
        """Test that Jupiter heliacal phenomena returns a valid result."""
        jd = julday(2024, 6, 1, 4)  # Early morning
        geopos = (-74.0060, 40.7128, 0.0)  # New York

        dret = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Jupiter",
            HELIACAL_RISING,
        )

        assert isinstance(dret, tuple)
        assert len(dret) == 50

    def test_mars_pheno_evening(self):
        """Test Mars heliacal phenomena for evening setting."""
        jd = julday(2024, 7, 20, 20)  # Evening
        geopos = (139.6503, 35.6762, 0.0)  # Tokyo

        dret = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Mars",
            HELIACAL_SETTING,
        )

        assert isinstance(dret, tuple)
        assert len(dret) == 50


class TestHeliacalPhenoValues:
    """Test that heliacal_pheno_ut returns sensible values."""

    def test_altitude_values_are_reasonable(self):
        """Test that altitude values are within valid range."""
        jd = julday(2024, 1, 15, 6)  # Dawn
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        dret = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Venus",
            HELIACAL_RISING,
        )

        # Altitudes should be in -90 to +90 range
        alt_o = dret[0]  # Object topocentric altitude
        app_alt_o = dret[1]  # Object apparent altitude
        geo_alt_o = dret[2]  # Object geocentric altitude
        alt_s = dret[4]  # Sun altitude

        assert -90 <= alt_o <= 90
        assert -90 <= app_alt_o <= 90
        assert -90 <= geo_alt_o <= 90
        assert -90 <= alt_s <= 90

    def test_azimuth_values_are_reasonable(self):
        """Test that azimuth values are within valid range."""
        jd = julday(2024, 1, 15, 6)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        dret = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Venus",
            HELIACAL_RISING,
        )

        azi_o = dret[3]  # Object azimuth
        azi_s = dret[5]  # Sun azimuth

        assert 0 <= azi_o <= 360
        assert 0 <= azi_s <= 360

    def test_arcus_visionis_values(self):
        """Test that arcus visionis values are calculated."""
        jd = julday(2024, 1, 15, 6)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        dret = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Venus",
            HELIACAL_RISING,
        )

        tav_act = dret[6]  # Topocentric arcus visionis
        arcv_act = dret[7]  # Geocentric arcus visionis

        # TAV is the altitude difference between body and Sun
        # Should be a reasonable value (typically -180 to +180)
        assert -180 <= tav_act <= 180
        assert -180 <= arcv_act <= 180

    def test_extinction_coefficient_is_positive(self):
        """Test that extinction coefficient is a sensible positive value."""
        jd = julday(2024, 1, 15, 6)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        dret = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Venus",
            HELIACAL_RISING,
        )

        k_act = dret[10]  # Extinction coefficient

        # Extinction coefficient typically 0.1 to 0.6
        assert 0 < k_act < 1.0

    def test_magnitude_for_venus(self):
        """Test that magnitude is calculated for Venus."""
        jd = julday(2024, 1, 15, 6)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        dret = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Venus",
            HELIACAL_RISING,
        )

        magnitude = dret[20]

        # Venus magnitude typically -4.5 to -3
        # But can vary, so just check it's a number
        assert isinstance(magnitude, float)

    def test_parallax_is_small(self):
        """Test that parallax is a small positive value."""
        jd = julday(2024, 1, 15, 6)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        dret = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Venus",
            HELIACAL_RISING,
        )

        parallax = dret[19]

        # Parallax for planets is typically < 1 degree
        assert 0 <= parallax < 1.0


class TestHeliacalPhenoEventTypes:
    """Test different event types for heliacal_pheno_ut."""

    def test_morning_first_event(self):
        """Test HELIACAL_RISING (morning first) event type."""
        jd = julday(2024, 3, 1, 5)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        dret = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Venus",
            HELIACAL_RISING,
        )

        assert len(dret) == 50

    def test_evening_last_event(self):
        """Test HELIACAL_SETTING (evening last) event type."""
        jd = julday(2024, 3, 1, 19)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        dret = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Venus",
            HELIACAL_SETTING,
        )

        assert len(dret) == 50

    def test_evening_first_event(self):
        """Test EVENING_FIRST event type."""
        jd = julday(2024, 3, 1, 19)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        dret = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Mercury",
            EVENING_FIRST,
        )

        assert len(dret) == 50

    def test_morning_last_event(self):
        """Test MORNING_LAST event type."""
        jd = julday(2024, 3, 1, 5)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        dret = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Mercury",
            MORNING_LAST,
        )

        assert len(dret) == 50


class TestHeliacalPhenoValidation:
    """Test input validation for heliacal_pheno_ut."""

    def test_invalid_body_raises_error(self):
        """Test that invalid body name raises ValueError."""
        jd = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)

        with pytest.raises(Error, match="not recognized"):
            heliacal_pheno_ut(
                jd,
                geopos,
                (1013.25, 15.0, 50.0, 0.0),
                (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                "InvalidBody",
                HELIACAL_RISING,
            )

    def test_invalid_event_type_raises_error(self):
        """Test that invalid event type raises ValueError."""
        jd = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)

        with pytest.raises(Error, match="Invalid event_type"):
            heliacal_pheno_ut(
                jd,
                geopos,
                (1013.25, 15.0, 50.0, 0.0),
                (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                "Venus",
                99,
            )


class TestHeliacalPhenoLocations:
    """Test heliacal_pheno_ut at various locations."""

    def test_northern_hemisphere(self):
        """Test at northern latitude."""
        jd = julday(2024, 6, 1, 3)  # Summer morning
        geopos = (24.9384, 60.1699, 0.0)  # Helsinki

        dret = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Jupiter",
            HELIACAL_RISING,
        )

        assert len(dret) == 50

    def test_southern_hemisphere(self):
        """Test at southern latitude."""
        jd = julday(2024, 1, 15, 5)  # Summer morning (S. hemisphere)
        geopos = (151.2093, -33.8688, 0.0)  # Sydney

        dret = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Saturn",
            HELIACAL_RISING,
        )

        assert len(dret) == 50

    def test_equatorial_location(self):
        """Test at equatorial location."""
        jd = julday(2024, 4, 1, 6)
        geopos = (32.5825, 0.3476, 0.0)  # Kampala, Uganda

        dret = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Mars",
            HELIACAL_RISING,
        )

        assert len(dret) == 50


class TestHeliacalPhenoAtmospheric:
    """Test heliacal_pheno_ut with different atmospheric conditions."""

    def test_with_altitude(self):
        """Test with observer at high altitude."""
        jd = julday(2024, 1, 15, 5)
        geopos_high = (12.4964, 41.9028, 2000.0)  # Rome, 2000 meters
        geopos_low = (12.4964, 41.9028, 0.0)  # Rome, sea level

        dret = heliacal_pheno_ut(
            jd,
            geopos_high,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Venus",
            HELIACAL_RISING,
        )

        assert len(dret) == 50
        # Extinction should be lower at high altitude
        k_high = dret[10]

        dret_low = heliacal_pheno_ut(
            jd,
            geopos_low,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Venus",
            HELIACAL_RISING,
        )
        k_low = dret_low[10]

        assert k_high < k_low  # Lower extinction at higher altitude

    def test_with_high_humidity(self):
        """Test with high humidity."""
        jd = julday(2024, 1, 15, 5)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        dret_low = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 15.0, 20.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Venus",
            HELIACAL_RISING,
        )
        dret_high = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 15.0, 90.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Venus",
            HELIACAL_RISING,
        )

        # Higher humidity should increase extinction
        k_low = dret_low[10]
        k_high = dret_high[10]

        assert k_high > k_low

    def test_pressure_does_not_affect_extinction(self):
        """Atmospheric pressure must not change the extinction coefficient.

        In Schaefer's VISLIMIT model (which the reference implements) the
        astronomical extinction coefficient depends on temperature,
        relative humidity, observer altitude and season, but NOT on the
        surface pressure (datm[0]); pressure only enters the refraction
        of the apparent altitude. The old library scaled Rayleigh
        extinction by pressure, which the reference does not do.
        """
        jd = julday(2024, 1, 15, 5)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        dret_normal = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Venus",
            HELIACAL_RISING,
        )
        dret_low = heliacal_pheno_ut(
            jd,
            geopos,
            (800.0, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Venus",
            HELIACAL_RISING,
        )

        # Extinction (dret[10]) is independent of surface pressure.
        assert dret_low[10] == pytest.approx(dret_normal[10], abs=1e-9)

    def test_with_higher_temperature(self):
        """Higher air temperature increases the extinction coefficient.

        The water-vapour extinction term scales as exp(T/15), so a warmer
        atmosphere (at fixed relative humidity) holds more water and is
        more extinctive, matching the reference.
        """
        jd = julday(2024, 1, 15, 5)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        dret_cold = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 0.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Venus",
            HELIACAL_RISING,
        )
        dret_warm = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 30.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Venus",
            HELIACAL_RISING,
        )

        assert dret_warm[10] > dret_cold[10]


class TestHeliacalPhenoAlias:
    """Test heliacal_pheno_ut alias."""

    def test_swe_alias_works(self):
        """Test that heliacal_pheno_ut and heliacal_pheno_ut return same data."""
        jd = julday(2024, 1, 15, 6)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 50.0, 0.0)
        dobs = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        result1 = heliacal_pheno_ut(jd, geopos, datm, dobs, "Venus", HELIACAL_RISING)

        result2 = heliacal_pheno_ut(jd, geopos, datm, dobs, "Venus", HELIACAL_RISING)

        # Both should return 50-element tuples with matching key fields
        assert len(result1) == 50
        assert len(result2) == 50
        # Altitudes and azimuths should match (same body, same time)
        assert abs(result1[0] - result2[0]) < 0.1  # AltO
        assert abs(result1[3] - result2[3]) < 0.1  # AziO
        assert abs(result1[20] - result2[20]) < 0.1  # Magnitude


class TestHeliacalPhenoMoon:
    """Test heliacal_pheno_ut specifically for Moon (crescent calculations)."""

    def test_moon_crescent_values(self):
        """Test that Moon-specific crescent values are calculated."""
        jd = julday(2024, 1, 12, 18)  # Near new moon
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        dret = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Moon",
            EVENING_FIRST,
        )

        w_moon = dret[16]  # Crescent width
        l_moon = dret[25]  # Crescent length
        q_yallop = dret[17]  # Yallop q-test
        illumination = dret[27]  # Illumination percentage

        # These values should be calculated for Moon
        assert isinstance(w_moon, float)
        assert isinstance(l_moon, float)
        assert isinstance(q_yallop, float)
        assert isinstance(illumination, float)


# =============================================================================
# PLANET-SPECIFIC HELIACAL TESTS
# =============================================================================


class TestHeliacalInnerOuterPlanets:
    """Test proper handling of inner vs outer planet geometry."""

    def test_is_inner_planet_mercury(self):
        """Test that Mercury is identified as an inner planet."""
        from libephemeris import is_inner_planet, INNER_PLANETS

        assert is_inner_planet(MERCURY)
        assert MERCURY in INNER_PLANETS

    def test_is_inner_planet_venus(self):
        """Test that Venus is identified as an inner planet."""
        from libephemeris import is_inner_planet, INNER_PLANETS

        assert is_inner_planet(VENUS)
        assert VENUS in INNER_PLANETS

    def test_is_inner_planet_mars(self):
        """Test that Mars is NOT an inner planet."""
        from libephemeris import is_inner_planet

        assert not is_inner_planet(MARS)

    def test_is_inner_planet_jupiter(self):
        """Test that Jupiter is NOT an inner planet."""
        from libephemeris import is_inner_planet

        assert not is_inner_planet(JUPITER)

    def test_is_inner_planet_saturn(self):
        """Test that Saturn is NOT an inner planet."""
        from libephemeris import is_inner_planet

        assert not is_inner_planet(SATURN)


class TestHeliacalOuterPlanetValidation:
    """Test that outer planets reject invalid event types."""

    def test_mars_evening_first_raises_error(self):
        """Test that Mars with EVENING_FIRST raises ValueError."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        with pytest.raises(Error, match="inner planets"):
            heliacal_ut(
                jd_start,
                geopos,
                (1013.25, 15.0, 50.0, 0.0),
                (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                "Mars",
                EVENING_FIRST,
            )

    def test_mars_morning_last_raises_error(self):
        """Test that Mars with MORNING_LAST raises ValueError."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        with pytest.raises(Error, match="inner planets"):
            heliacal_ut(
                jd_start,
                geopos,
                (1013.25, 15.0, 50.0, 0.0),
                (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                "Mars",
                MORNING_LAST,
            )

    def test_jupiter_evening_first_raises_error(self):
        """Test that Jupiter with EVENING_FIRST raises ValueError."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        with pytest.raises(Error, match="inner planets"):
            heliacal_ut(
                jd_start,
                geopos,
                (1013.25, 15.0, 50.0, 0.0),
                (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                "Jupiter",
                EVENING_FIRST,
            )

    def test_jupiter_morning_last_raises_error(self):
        """Test that Jupiter with MORNING_LAST raises ValueError."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        with pytest.raises(Error, match="inner planets"):
            heliacal_ut(
                jd_start,
                geopos,
                (1013.25, 15.0, 50.0, 0.0),
                (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                "Jupiter",
                MORNING_LAST,
            )

    def test_saturn_evening_first_raises_error(self):
        """Test that Saturn with EVENING_FIRST raises ValueError."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        with pytest.raises(Error, match="inner planets"):
            heliacal_ut(
                jd_start,
                geopos,
                (1013.25, 15.0, 50.0, 0.0),
                (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                "Saturn",
                EVENING_FIRST,
            )

    def test_saturn_morning_last_raises_error(self):
        """Test that Saturn with MORNING_LAST raises ValueError."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        with pytest.raises(Error, match="inner planets"):
            heliacal_ut(
                jd_start,
                geopos,
                (1013.25, 15.0, 50.0, 0.0),
                (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                "Saturn",
                MORNING_LAST,
            )


class TestHeliacalInnerPlanetEventTypes:
    """Test that inner planets accept all event types."""

    def test_mercury_heliacal_rising(self):
        """Test Mercury heliacal rising works."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Mercury",
            HELIACAL_RISING,
        )

        if result[0] > 0:
            assert result[0] > jd_start

    def test_mercury_heliacal_setting(self):
        """Test Mercury heliacal setting works."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Mercury",
            HELIACAL_SETTING,
        )

        if result[0] > 0:
            assert result[0] > jd_start

    def test_mercury_evening_first(self):
        """Test Mercury evening first visibility works."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Mercury",
            EVENING_FIRST,
        )

        # Should not raise an error - inner planets can use this event type
        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_mercury_morning_last(self):
        """Test Mercury morning last visibility works."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Mercury",
            MORNING_LAST,
        )

        # Should not raise an error - inner planets can use this event type
        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_venus_evening_first(self):
        """Test Venus evening first visibility works."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Venus",
            EVENING_FIRST,
        )

        # Should not raise an error - inner planets can use this event type
        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_venus_morning_last(self):
        """Test Venus morning last visibility works."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Venus",
            MORNING_LAST,
        )

        # Should not raise an error - inner planets can use this event type
        assert isinstance(result, tuple)
        assert len(result) == 3


class TestHeliacalOuterPlanetRisingSetting:
    """Test that outer planets work correctly with heliacal rising/setting."""

    def test_mars_heliacal_rising(self):
        """Test Mars heliacal rising works."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Mars",
            HELIACAL_RISING,
        )

        if result[0] > 0:
            assert result[0] > jd_start

    def test_mars_heliacal_setting(self):
        """Test Mars heliacal setting works."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Mars",
            HELIACAL_SETTING,
        )

        if result[0] > 0:
            assert result[0] > jd_start

    def test_jupiter_heliacal_rising(self):
        """Test Jupiter heliacal rising works."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Jupiter",
            HELIACAL_RISING,
        )

        if result[0] > 0:
            assert result[0] > jd_start

    def test_jupiter_heliacal_setting(self):
        """Test Jupiter heliacal setting works."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Jupiter",
            HELIACAL_SETTING,
        )

        if result[0] > 0:
            assert result[0] > jd_start

    def test_saturn_heliacal_rising(self):
        """Test Saturn heliacal rising works."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Saturn",
            HELIACAL_RISING,
        )

        if result[0] > 0:
            assert result[0] > jd_start

    def test_saturn_heliacal_setting(self):
        """Test Saturn heliacal setting works."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Saturn",
            HELIACAL_SETTING,
        )

        if result[0] > 0:
            assert result[0] > jd_start


class TestSweHeliacalUtOuterPlanetValidation:
    """Test that heliacal_ut also validates inner/outer planet constraints."""

    def test_mars_evening_first_via_swe_api_raises_error(self):
        """Test that Mars with EVENING_FIRST raises ValueError via swe API."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        with pytest.raises(Error, match="inner planets"):
            heliacal_ut(jd_start, geopos, datm, dobs, "Mars", EVENING_FIRST)

    def test_jupiter_morning_last_via_swe_api_raises_error(self):
        """Test that Jupiter with MORNING_LAST raises ValueError via swe API."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        with pytest.raises(Error, match="inner planets"):
            heliacal_ut(jd_start, geopos, datm, dobs, "Jupiter", MORNING_LAST)

    def test_saturn_evening_first_via_swe_api_raises_error(self):
        """Test that Saturn with EVENING_FIRST raises ValueError via swe API."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        with pytest.raises(Error, match="inner planets"):
            heliacal_ut(jd_start, geopos, datm, dobs, "Saturn", EVENING_FIRST)

    def test_mercury_evening_first_via_swe_api_works(self):
        """Test that Mercury with EVENING_FIRST works via swe API."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        # Should not raise an error
        result = heliacal_ut(jd_start, geopos, datm, dobs, "Mercury", EVENING_FIRST)

        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_venus_morning_last_via_swe_api_works(self):
        """Test that Venus with MORNING_LAST works via swe API."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        # Should not raise an error
        result = heliacal_ut(jd_start, geopos, datm, dobs, "Venus", MORNING_LAST)

        assert isinstance(result, tuple)
        assert len(result) == 3


class TestPhenoParallaxInAltitude:
    """dret[19] (ParO) is the parallax in altitude, not the horizontal parallax.

    Measured black-box: the reference's dret[19] equals dret[2] - dret[0]
    (GeoAlt - AltO) on every body/instant probed — i.e. ~HP*cos(alt), which
    only coincides with the horizontal parallax at the horizon. Frozen
    oracle value for the Moon at appreciable altitude (~26 deg): 0.920914.
    """

    def test_moon_parallax_in_altitude(self):
        jd = julday(2023, 3, 23, 16.5)
        dret = heliacal_pheno_ut(
            jd,
            (31.2, 29.9, 0.0),
            (1013.25, 15.0, 40.0, 0.0),
            (36.0, 1.0, 0.0, 0.0, 0.0, 0.0),
            "Moon",
            3,
        )
        assert abs(dret[19] - 0.920914) < 5e-4
        # Internal consistency: slot 19 == GeoAlt - AltO by construction.
        assert abs(dret[19] - (dret[2] - dret[0])) < 1e-9


class TestPhenoWindowObserver:
    """The visibility-window slots respect the caller's observer (round-W fix).

    Both pheno paths previously hardcoded the naked observer into the
    rise-window computation, so dret[12..14]/dret[24] never moved with the
    observer tuple. Structural assertion: an old presbyopic observer and a
    sharp young one get different windows for the same event (the absolute
    values carry the documented VISLIMIT model floor, so no oracle freeze).
    """

    def test_window_slots_vary_with_observer(self):
        jd = 2460539.6522778766
        args = ((12.5, 41.9, 0.0), (1013.25, 15.0, 50.0, 0.0))
        r_young = heliacal_pheno_ut(
            jd, args[0], args[1], (20.0, 1.5, 0.0, 0.0, 0.0, 0.0), "Sirius", 1
        )
        r_old = heliacal_pheno_ut(
            jd, args[0], args[1], (75.0, 0.7, 0.0, 0.0, 0.0, 0.0), "Sirius", 1
        )
        assert any(abs(r_young[i] - r_old[i]) > 1e-6 for i in (11, 12, 13, 14, 24))


@pytest.mark.unit
class TestPhenoWindowSpikeAtEdge:
    """The twilight visibility window is a narrow spike at the sunset edge
    of the 4-hour bracket; a bare golden-section over the bracket missed it
    (dret[12..14,24] came back as the no-window sentinel). The grid+refine
    search must find the real window (reference finds ~40 min for this
    case; exact instants carry the documented VISLIMIT model floor)."""

    def test_venus_evening_window_exists(self):
        import libephemeris as le

        jd = le.julday(2024, 9, 1, 4.0)
        d = le.heliacal_pheno_ut(
            jd,
            (12.5, 41.9, 0),
            (1013.25, 15, 40, 0),
            (36, 1, 0, 0, 0, 0),
            "Venus",
            3,
        )
        assert d[12] < 99999998.0, "TfirstVR missing"
        assert d[13] < 99999998.0, "TbVR missing"
        assert d[14] < 99999998.0, "TlastVR missing"
        # Window ~40 min (oracle 0.0274 d); allow the documented model floor.
        assert 0.01 < d[24] < 0.06
        assert d[12] < d[13] < d[14]
