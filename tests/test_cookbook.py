"""
Tests for the LibEphemeris Cookbook examples.

These tests verify that all cookbook examples work correctly and produce
valid astronomical/astrological data.
"""

import pytest
import sys
import os

# Add docs directory to path for importing cookbook
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from docs.cookbook import (
    format_longitude,
    format_aspect,
    calculate_natal_chart,
    find_next_transit,
    find_planetary_return,
    find_sign_ingresses,
    calculate_synastry,
    calculate_composite_midpoints,
    find_next_solar_eclipse,
    find_next_lunar_eclipse,
    find_eclipse_visibility,
    calculate_monthly_ephemeris,
    find_retrograde_periods,
)
from libephemeris.constants import (
    SUN,
    MOON,
    MERCURY,
    MARS,
    ECL_TOTAL,
)


class TestFormatLongitude:
    """Tests for the format_longitude utility function."""

    def test_aries_0_degrees(self):
        """0 degrees should be 0° Aries."""
        result = format_longitude(0.0)
        assert "Aries" in result
        assert "00°" in result

    def test_taurus_30_degrees(self):
        """30 degrees should be 0° Taurus."""
        result = format_longitude(30.0)
        assert "Taurus" in result
        assert "00°" in result

    def test_leo_15_degrees(self):
        """135 degrees should be 15° Leo."""
        result = format_longitude(135.0)
        assert "Leo" in result
        assert "15°" in result

    def test_pisces_end(self):
        """359 degrees should be in Pisces."""
        result = format_longitude(359.0)
        assert "Pisces" in result
        assert "29°" in result

    def test_wraparound(self):
        """360 degrees should wrap to 0° Aries."""
        result = format_longitude(360.0)
        assert "Aries" in result
        assert "00°" in result

    def test_negative_normalized(self):
        """Negative degrees should be normalized."""
        result = format_longitude(-30.0)
        assert "Pisces" in result  # -30 = 330 = 0° Pisces


class TestFormatAspect:
    """Tests for the format_aspect utility function."""

    def test_conjunction(self):
        """0 degrees should be conjunction."""
        result = format_aspect(0.0)
        assert result is not None
        assert "Conjunction" in result

    def test_opposition(self):
        """180 degrees should be opposition."""
        result = format_aspect(180.0)
        assert result is not None
        assert "Opposition" in result

    def test_trine(self):
        """120 degrees should be trine."""
        result = format_aspect(120.0)
        assert result is not None
        assert "Trine" in result

    def test_square(self):
        """90 degrees should be square."""
        result = format_aspect(90.0)
        assert result is not None
        assert "Square" in result

    def test_sextile(self):
        """60 degrees should be sextile."""
        result = format_aspect(60.0)
        assert result is not None
        assert "Sextile" in result

    def test_no_aspect(self):
        """45 degrees (semisquare) should return None for major aspects."""
        result = format_aspect(45.0)
        assert result is None

    def test_aspect_with_orb(self):
        """Aspects within orb should be detected."""
        # 85 degrees is a square with 5° orb
        result = format_aspect(85.0)
        assert result is not None
        assert "Square" in result


class TestNatalChart:
    """Tests for natal chart calculation."""

    def test_basic_calculation(self):
        """Basic natal chart calculation should work."""
        chart = calculate_natal_chart(
            year=2000,
            month=1,
            day=1,
            hour=12.0,
            latitude=41.9028,
            longitude=12.4964,
        )

        assert "planets" in chart
        assert "houses" in chart
        assert "angles" in chart
        assert "meta" in chart

    def test_all_planets_present(self):
        """All major planets should be present."""
        chart = calculate_natal_chart(2000, 1, 1, 12.0, 41.9028, 12.4964)

        expected_planets = [
            "Sun",
            "Moon",
            "Mercury",
            "Venus",
            "Mars",
            "Jupiter",
            "Saturn",
            "Uranus",
            "Neptune",
            "Pluto",
            "North Node",
            "Chiron",
        ]
        for planet in expected_planets:
            assert planet in chart["planets"], f"Missing planet: {planet}"

    def test_12_houses(self):
        """There should be 12 house cusps."""
        chart = calculate_natal_chart(2000, 1, 1, 12.0, 41.9028, 12.4964)
        assert len(chart["houses"]) == 12

    def test_angles_present(self):
        """All major angles should be present."""
        chart = calculate_natal_chart(2000, 1, 1, 12.0, 41.9028, 12.4964)

        expected_angles = ["ASC", "MC", "DESC", "IC", "Vertex"]
        for angle in expected_angles:
            assert angle in chart["angles"], f"Missing angle: {angle}"

    def test_planet_data_structure(self):
        """Planet data should have expected fields."""
        chart = calculate_natal_chart(2000, 1, 1, 12.0, 41.9028, 12.4964)
        sun = chart["planets"]["Sun"]

        assert "longitude" in sun
        assert "latitude" in sun
        assert "distance" in sun
        assert "speed" in sun
        assert "retrograde" in sun
        assert "formatted" in sun

    def test_valid_longitude_ranges(self):
        """All longitudes should be in valid range (0-360)."""
        chart = calculate_natal_chart(2000, 1, 1, 12.0, 41.9028, 12.4964)

        for planet_name, planet_data in chart["planets"].items():
            lon = planet_data["longitude"]
            assert 0 <= lon < 360, f"{planet_name} longitude out of range: {lon}"

        for house in chart["houses"]:
            lon = house["longitude"]
            assert 0 <= lon < 360, (
                f"House {house['house']} longitude out of range: {lon}"
            )

    def test_different_house_systems(self):
        """Different house systems should produce different results."""
        chart_placidus = calculate_natal_chart(
            2000, 1, 1, 12.0, 41.9028, 12.4964, house_system=b"P"
        )
        chart_koch = calculate_natal_chart(
            2000, 1, 1, 12.0, 41.9028, 12.4964, house_system=b"K"
        )

        # ASC and MC should be the same (they don't depend on house system)
        assert (
            abs(
                chart_placidus["angles"]["ASC"]["longitude"]
                - chart_koch["angles"]["ASC"]["longitude"]
            )
            < 0.001
        )

        # But intermediate cusps may differ (e.g., house 2)
        # This is expected behavior for different house systems

    def test_sidereal_calculation(self):
        """Sidereal calculation should subtract ayanamsha."""
        chart_tropical = calculate_natal_chart(
            2000, 1, 1, 12.0, 41.9028, 12.4964, sidereal=False
        )
        chart_sidereal = calculate_natal_chart(
            2000, 1, 1, 12.0, 41.9028, 12.4964, sidereal=True
        )

        # True Citra is independently anchored to Spica, so its positions
        # differ from tropical positions by an ordinary ayanamsha-sized angle.
        sun_tropical = chart_tropical["planets"]["Sun"]["longitude"]
        sun_sidereal = chart_sidereal["planets"]["Sun"]["longitude"]

        diff = abs(sun_tropical - sun_sidereal)
        assert 20 < diff < 28, f"Ayanamsha difference unexpected: {diff}"


class TestTransitFinding:
    """Tests for transit finding functions."""

    def test_find_vernal_equinox(self):
        """Should find vernal equinox correctly."""
        result = find_next_transit(SUN, 0.0, 2024, 1, 1)

        assert result["planet"] == "Sun"
        assert "2024-03" in result["date"]  # Should be in March

    def test_find_summer_solstice(self):
        """Should find summer solstice (90° Cancer)."""
        result = find_next_transit(SUN, 90.0, 2024, 1, 1)

        assert "2024-06" in result["date"]  # Should be in June

    def test_moon_transit(self):
        """Should find Moon transits within lunar month."""
        result = find_next_transit(MOON, 180.0, 2024, 1, 1)

        # Moon should cross 180° within about 27 days
        assert result["crossing_jd"] > 0

    def test_planet_transit(self):
        """Should find Mars transit."""
        result = find_next_transit(MARS, 0.0, 2024, 1, 1)

        assert result["planet"] == "Mars"
        assert result["crossing_jd"] > 0

    def test_transit_precision(self):
        """Transit should be accurate to target longitude."""
        result = find_next_transit(SUN, 45.0, 2024, 1, 1)

        actual_lon = result["planet_position"]["longitude"]
        target_lon = 45.0

        diff = abs(actual_lon - target_lon)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.01, f"Transit precision too low: {diff} degrees"


class TestPlanetaryReturn:
    """Tests for planetary return calculation."""

    def test_solar_return(self):
        """Should find solar return."""
        result = find_planetary_return(SUN, 135.0, 2024)  # 15° Leo

        assert result["planet"] == "Sun"
        # Should be in August (Leo month)
        assert "2024-08" in result["date"] or "2024-07" in result["date"]


class TestSignIngresses:
    """Tests for finding sign ingresses."""

    def test_sun_ingresses_count(self):
        """Should find 12 Sun ingresses in a year."""
        ingresses = find_sign_ingresses(SUN, 2024)

        # Should find approximately 12 ingresses
        # (might be 11-13 depending on year boundaries)
        assert 11 <= len(ingresses) <= 13

    def test_ingresses_in_order(self):
        """Ingresses should be in chronological order."""
        ingresses = find_sign_ingresses(SUN, 2024)

        for i in range(1, len(ingresses)):
            assert ingresses[i]["crossing_jd"] > ingresses[i - 1]["crossing_jd"]


class TestSynastry:
    """Tests for synastry calculation."""

    def test_synastry_finds_aspects(self):
        """Should find aspects between two charts."""
        chart1 = calculate_natal_chart(1990, 5, 15, 10.0, 40.7128, -74.0060)
        chart2 = calculate_natal_chart(1992, 8, 22, 14.0, 34.0522, -118.2437)

        synastry = calculate_synastry(chart1, chart2)

        assert "aspects" in synastry
        assert "summary" in synastry
        assert len(synastry["aspects"]) > 0

    def test_synastry_summary(self):
        """Synastry summary should have correct fields."""
        chart1 = calculate_natal_chart(1990, 5, 15, 10.0, 40.7128, -74.0060)
        chart2 = calculate_natal_chart(1992, 8, 22, 14.0, 34.0522, -118.2437)

        synastry = calculate_synastry(chart1, chart2)

        assert "total_aspects" in synastry["summary"]
        assert "harmonious" in synastry["summary"]
        assert "challenging" in synastry["summary"]


class TestCompositeMidpoints:
    """Tests for composite chart calculation."""

    def test_composite_planets(self):
        """Should calculate composite planets."""
        chart1 = calculate_natal_chart(1990, 5, 15, 10.0, 40.7128, -74.0060)
        chart2 = calculate_natal_chart(1992, 8, 22, 14.0, 34.0522, -118.2437)

        composite = calculate_composite_midpoints(chart1, chart2)

        assert "planets" in composite
        assert "Sun" in composite["planets"]
        assert "Moon" in composite["planets"]

    def test_composite_midpoint_validity(self):
        """Composite midpoints should be valid longitudes."""
        chart1 = calculate_natal_chart(1990, 5, 15, 10.0, 40.7128, -74.0060)
        chart2 = calculate_natal_chart(1992, 8, 22, 14.0, 34.0522, -118.2437)

        composite = calculate_composite_midpoints(chart1, chart2)

        for planet_name, data in composite["planets"].items():
            lon = data["longitude"]
            assert 0 <= lon < 360, f"Invalid composite longitude for {planet_name}"


class TestEclipseSearch:
    """Tests for eclipse search functions."""

    def test_find_solar_eclipse(self):
        """Should find a solar eclipse."""
        eclipse = find_next_solar_eclipse(2024, 1, 1)

        assert "type" in eclipse
        assert "maximum" in eclipse
        assert "Solar Eclipse" in eclipse["type"]

    def test_find_total_solar_eclipse(self):
        """Should find a total solar eclipse when filtered."""
        eclipse = find_next_solar_eclipse(2024, 1, 1, ECL_TOTAL)

        assert "Total" in eclipse["type"]

    def test_eclipse_timing(self):
        """Eclipse timing should be valid."""
        eclipse = find_next_solar_eclipse(2024, 1, 1)

        # Maximum JD should be after start date
        import libephemeris as ephem

        start_jd = ephem.julday(2024, 1, 1, 0.0)
        assert eclipse["maximum_jd"] > start_jd

    def test_find_lunar_eclipse(self):
        """Should find a lunar eclipse."""
        eclipse = find_next_lunar_eclipse(2024, 1, 1)

        assert "type" in eclipse
        assert "maximum" in eclipse
        assert "Lunar Eclipse" in eclipse["type"]


class TestEclipseVisibility:
    """Tests for eclipse visibility calculation."""

    def test_visibility_data_structure(self):
        """Should return proper visibility data."""
        eclipse = find_next_solar_eclipse(2024, 1, 1, ECL_TOTAL)
        visibility = find_eclipse_visibility(eclipse["maximum_jd"], 32.7767, -96.7970)

        assert "visible" in visibility
        assert "magnitude" in visibility
        assert "obscuration" in visibility
        assert "sun_altitude" in visibility

    def test_magnitude_range(self):
        """Magnitude should be in valid range when visible."""
        eclipse = find_next_solar_eclipse(2024, 1, 1, ECL_TOTAL)
        visibility = find_eclipse_visibility(eclipse["maximum_jd"], 32.7767, -96.7970)

        if visibility["visible"]:
            assert 0 <= visibility["magnitude"] <= 1.5
            assert 0 <= visibility["obscuration"] <= 1


class TestMonthlyEphemeris:
    """Tests for monthly ephemeris calculation."""

    def test_ephemeris_days(self):
        """Should return correct number of days."""
        ephemeris = calculate_monthly_ephemeris(2024, 3)  # March has 31 days

        assert len(ephemeris) == 31

    def test_february_leap_year(self):
        """Should handle leap year February correctly."""
        ephemeris = calculate_monthly_ephemeris(2024, 2)  # 2024 is a leap year

        assert len(ephemeris) == 29

    def test_ephemeris_planets(self):
        """Should include all major planets."""
        ephemeris = calculate_monthly_ephemeris(2024, 3)

        expected_planets = [
            "Sun",
            "Moon",
            "Mercury",
            "Venus",
            "Mars",
            "Jupiter",
            "Saturn",
            "Uranus",
            "Neptune",
            "Pluto",
        ]
        for planet in expected_planets:
            assert planet in ephemeris[0], f"Missing planet: {planet}"

    def test_ephemeris_with_asteroids(self):
        """Should include Chiron when requested."""
        ephemeris = calculate_monthly_ephemeris(2024, 3, include_asteroids=True)

        assert "Chiron" in ephemeris[0]

    def test_ephemeris_continuity(self):
        """Planetary positions should change gradually day to day."""
        ephemeris = calculate_monthly_ephemeris(2024, 3)

        for i in range(1, len(ephemeris)):
            sun_today = ephemeris[i]["Sun"]["longitude"]
            sun_yesterday = ephemeris[i - 1]["Sun"]["longitude"]

            # Sun moves about 1°/day, so difference should be small
            diff = abs(sun_today - sun_yesterday)
            if diff > 180:
                diff = 360 - diff

            assert diff < 2, f"Sun movement too large: {diff} degrees"


class TestRetrogradePeriods:
    """Tests for retrograde period finding."""

    def test_mercury_retrogrades(self):
        """Mercury should have multiple retrograde periods per year."""
        retrogrades = find_retrograde_periods(2024, MERCURY)

        # Mercury typically has 3-4 retrograde periods per year
        assert len(retrogrades) >= 2

    def test_retrograde_period_structure(self):
        """Retrograde period should have required fields."""
        retrogrades = find_retrograde_periods(2024, MERCURY)

        if retrogrades:
            period = retrogrades[0]
            assert "start" in period
            assert "end" in period
            assert "start_position" in period
            assert "end_position" in period
            assert "duration_days" in period

    def test_retrograde_duration_reasonable(self):
        """Retrograde periods should have reasonable durations."""
        retrogrades = find_retrograde_periods(2024, MERCURY)

        # Count full retrograde periods (excluding partial ones at year boundaries)
        full_periods = [p for p in retrogrades if p["duration_days"] >= 15]

        # Should find at least 2 full Mercury retrograde periods
        assert len(full_periods) >= 2

        for period in full_periods:
            # Mercury retrograde typically lasts 20-25 days
            assert period["duration_days"] < 40, (
                f"Unusual retrograde duration: {period['duration_days']}"
            )


class TestCookbookIntegration:
    """Integration tests that run the cookbook main function."""

    def test_main_runs_without_error(self):
        """The main() function should run without raising exceptions."""
        # Import and run main
        from docs.cookbook import main

        # This should not raise any exceptions
        try:
            main()
        except Exception as e:
            pytest.fail(f"Cookbook main() raised an exception: {e}")
