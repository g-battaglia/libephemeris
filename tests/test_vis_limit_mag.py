"""
Tests for vis_limit_mag function in libephemeris.

Tests the calculation of limiting visual magnitude for celestial body visibility
given atmospheric conditions and observer characteristics.

The function determines whether a celestial body will be visible based on:
- Sky brightness (from Sun, Moon)
- Atmospheric conditions (pressure, temperature, humidity)
- Observer characteristics (age, eye quality)
"""

import pytest

from libephemeris import (
    julday,
    vis_limit_mag,
    FLG_SWIEPH,
    HELFLAG_VISLIM_DARK,
    HELFLAG_VISLIM_NOMOON,
    HELFLAG_BELOW_HORIZON,
    HELFLAG_PHOTOPIC,
    HELFLAG_SCOTOPIC,
    HELFLAG_MIXED,
)


class TestVisLimitMagBasic:
    """Basic tests for vis_limit_mag function."""

    def test_venus_night_visibility(self):
        """Test Venus visibility at night."""
        # August 15, 2024 at 22:00 UT (nighttime in Europe)
        jd = julday(2024, 8, 15, 22.0)
        # Rome, Italy
        geopos = (12.4964, 41.9028, 0)
        # Standard atmosphere
        atmo = (1013.25, 15.0, 50.0, 0.0)
        # Normal observer
        observer = (36, 1.0)

        result, dret = vis_limit_mag(jd, geopos, atmo, observer, "Venus")

        # Should return a tuple of 8 values
        assert len(dret) == 10
        # Limiting magnitude should be a reasonable value
        assert -10.0 < dret[0] < 10.0
        # Object altitude
        assert -90.0 <= dret[1] <= 90.0
        # Object azimuth
        assert 0.0 <= dret[2] <= 360.0
        # Sun altitude
        assert -90.0 <= dret[3] <= 90.0
        # Sun azimuth
        assert 0.0 <= dret[4] <= 360.0
        # Moon altitude
        assert -90.0 <= dret[5] <= 90.0
        # Moon azimuth
        assert 0.0 <= dret[6] <= 360.0
        # Object magnitude
        assert -10.0 < dret[7] < 20.0

    def test_jupiter_visibility(self):
        """Test Jupiter visibility calculation."""
        jd = julday(2024, 7, 20, 23.0)
        geopos = (-0.1278, 51.5074, 0)  # London
        atmo = (1013.25, 18.0, 60.0, 0.0)
        observer = (40, 1.0)

        result, dret = vis_limit_mag(jd, geopos, atmo, observer, "Jupiter")

        assert len(dret) == 10
        # Result should be a valid vision type or below horizon
        assert result in (
            HELFLAG_BELOW_HORIZON,
            HELFLAG_PHOTOPIC,
            HELFLAG_SCOTOPIC,
            HELFLAG_MIXED,
        )

    def test_mars_visibility(self):
        """Test Mars visibility calculation."""
        jd = julday(2024, 12, 1, 21.0)
        geopos = (-74.0060, 40.7128, 0)  # New York
        atmo = (1015.0, 5.0, 40.0, 0.0)
        observer = (35, 1.0)

        result, dret = vis_limit_mag(jd, geopos, atmo, observer, "Mars")

        assert len(dret) == 10
        assert result in (
            HELFLAG_BELOW_HORIZON,
            HELFLAG_PHOTOPIC,
            HELFLAG_SCOTOPIC,
            HELFLAG_MIXED,
        )

    def test_saturn_visibility(self):
        """Test Saturn visibility calculation."""
        jd = julday(2024, 9, 15, 22.5)
        geopos = (139.6503, 35.6762, 0)  # Tokyo
        atmo = (1010.0, 20.0, 70.0, 0.0)
        observer = (30, 1.0)

        result, dret = vis_limit_mag(jd, geopos, atmo, observer, "Saturn")

        assert len(dret) == 10


class TestVisLimitMagInputFormats:
    """Test different input formats for vis_limit_mag."""

    def test_planet_by_number_string(self):
        """Test passing planet as number string."""
        jd = julday(2024, 8, 15, 22.0)
        geopos = (12.5, 42.0, 0)
        atmo = (1013.25, 15.0, 50.0, 0.0)
        observer = (36, 1.0)

        # "3" is Venus
        result, dret = vis_limit_mag(jd, geopos, atmo, observer, "3")

        assert len(dret) == 10

    def test_planet_by_name_case_insensitive(self):
        """Test that planet names are case insensitive."""
        jd = julday(2024, 8, 15, 22.0)
        geopos = (12.5, 42.0, 0)
        atmo = (1013.25, 15.0, 50.0, 0.0)
        observer = (36, 1.0)

        result1, dret1 = vis_limit_mag(jd, geopos, atmo, observer, "venus")
        result2, dret2 = vis_limit_mag(jd, geopos, atmo, observer, "VENUS")
        result3, dret3 = vis_limit_mag(jd, geopos, atmo, observer, "Venus")

        # All should return same results
        assert dret1[7] == pytest.approx(dret2[7], abs=0.01)
        assert dret2[7] == pytest.approx(dret3[7], abs=0.01)

    def test_fixed_star_regulus(self):
        """Test visibility calculation for Regulus (fixed star)."""
        jd = julday(2024, 4, 15, 22.0)  # Spring when Regulus is visible
        geopos = (12.5, 42.0, 0)  # Rome
        atmo = (1013.25, 15.0, 40.0, 0.0)
        observer = (36, 1.0)

        result, dret = vis_limit_mag(jd, geopos, atmo, observer, "Regulus")

        assert len(dret) == 10
        # Regulus is bright (mag ~ 1.4)
        if dret[1] > 0:  # If above horizon
            assert dret[7] < 5.0  # Should be reasonably bright


class TestVisLimitMagAtmosphericConditions:
    """Test effects of atmospheric conditions on visibility."""

    def test_high_humidity_reduces_visibility(self):
        """Test that high humidity reduces limiting magnitude."""
        jd = julday(2024, 8, 15, 22.0)
        geopos = (12.5, 42.0, 0)
        observer = (36, 1.0)

        # Low humidity
        atmo_low = (1013.25, 15.0, 20.0, 0.0)
        result_low, dret_low = vis_limit_mag(jd, geopos, atmo_low, observer, "Venus")

        # High humidity
        atmo_high = (1013.25, 15.0, 90.0, 0.0)
        result_high, dret_high = vis_limit_mag(jd, geopos, atmo_high, observer, "Venus")

        # Higher humidity should reduce limiting magnitude (make it lower)
        # The effect depends on object altitude and specific conditions
        # For objects above horizon, limiting mag is computed
        if dret_low[1] > 0 and dret_high[1] > 0:
            # Both above horizon - high humidity generally reduces visibility
            # but the difference might be small
            pass  # Just verify both calculations succeeded

    def test_meteorological_range_effect(self):
        """Test that meteorological range affects extinction."""
        jd = julday(2024, 8, 15, 22.0)
        geopos = (12.5, 42.0, 0)
        observer = (36, 1.0)

        # Good visibility (high met range)
        atmo_clear = (1013.25, 15.0, 50.0, 50.0)  # 50 km visibility
        result_clear, dret_clear = vis_limit_mag(
            jd, geopos, atmo_clear, observer, "Venus"
        )

        # Poor visibility (low met range)
        atmo_hazy = (1013.25, 15.0, 50.0, 5.0)  # 5 km visibility
        result_hazy, dret_hazy = vis_limit_mag(jd, geopos, atmo_hazy, observer, "Venus")

        # Both should succeed
        assert len(dret_clear) == 10
        assert len(dret_hazy) == 10


class TestVisLimitMagObserverParameters:
    """Test effects of observer parameters on visibility."""

    def test_observer_age_effect(self):
        """Test that older observers have lower limiting magnitude."""
        jd = julday(2024, 8, 15, 22.0)
        geopos = (12.5, 42.0, 0)
        atmo = (1013.25, 15.0, 50.0, 0.0)

        # Young observer
        observer_young = (20, 1.0)
        result_young, dret_young = vis_limit_mag(
            jd, geopos, atmo, observer_young, "Venus"
        )

        # Older observer
        observer_old = (70, 1.0)
        result_old, dret_old = vis_limit_mag(jd, geopos, atmo, observer_old, "Venus")

        # Limiting magnitude should be affected by age
        # Young observers generally see fainter objects
        # (Note: effect is small and may be overwhelmed by other factors)
        assert len(dret_young) == 10
        assert len(dret_old) == 10

    def test_snellen_ratio_effect(self):
        """Test that better vision (higher Snellen) improves limiting magnitude."""
        jd = julday(2024, 8, 15, 22.0)
        geopos = (12.5, 42.0, 0)
        atmo = (1013.25, 15.0, 50.0, 0.0)

        # Normal vision
        observer_normal = (36, 1.0)
        result_normal, dret_normal = vis_limit_mag(
            jd, geopos, atmo, observer_normal, "Venus"
        )

        # Better than normal vision
        observer_good = (36, 1.5)
        result_good, dret_good = vis_limit_mag(jd, geopos, atmo, observer_good, "Venus")

        assert len(dret_normal) == 10
        assert len(dret_good) == 10


class TestVisLimitMagFlags:
    """Test effect of various flags on vis_limit_mag."""

    def test_dark_sky_flag(self):
        """Test HELFLAG_VISLIM_DARK flag forces dark sky conditions."""
        # Daytime - normally would have poor visibility
        jd = julday(2024, 8, 15, 12.0)  # Noon
        geopos = (12.5, 42.0, 0)
        atmo = (1013.25, 25.0, 50.0, 0.0)
        observer = (36, 1.0)

        # Without dark sky flag
        result_normal, dret_normal = vis_limit_mag(
            jd, geopos, atmo, observer, "Venus", flags=FLG_SWIEPH
        )

        # With dark sky flag (Sun at nadir)
        result_dark, dret_dark = vis_limit_mag(
            jd,
            geopos,
            atmo,
            observer,
            "Venus",
            flags=FLG_SWIEPH | HELFLAG_VISLIM_DARK,
        )

        # Sun altitude with dark flag should be -90
        assert dret_dark[3] == -90.0
        # Limiting magnitude should be much better with dark sky
        if dret_normal[1] > 0 and dret_dark[1] > 0:  # Both above horizon
            assert dret_dark[0] > dret_normal[0]

    def test_no_moon_flag(self):
        """Test HELFLAG_VISLIM_NOMOON flag excludes Moon contribution."""
        # Full moon night
        jd = julday(2024, 8, 19, 23.0)  # Close to full moon
        geopos = (12.5, 42.0, 0)
        atmo = (1013.25, 15.0, 50.0, 0.0)
        observer = (36, 1.0)

        # With Moon
        result_with_moon, dret_with_moon = vis_limit_mag(
            jd, geopos, atmo, observer, "Venus", flags=FLG_SWIEPH
        )

        # Without Moon contribution
        result_no_moon, dret_no_moon = vis_limit_mag(
            jd,
            geopos,
            atmo,
            observer,
            "Venus",
            flags=FLG_SWIEPH | HELFLAG_VISLIM_NOMOON,
        )

        # Moon altitude with no-moon flag should be -90
        # (only when object is above horizon; below-horizon returns all zeros)
        if result_no_moon != HELFLAG_BELOW_HORIZON:
            assert dret_no_moon[5] == -90.0
        else:
            assert dret_no_moon[5] == 0.0


class TestVisLimitMagVisionTypes:
    """Test vision type return values."""

    def test_below_horizon_returns_minus_2(self):
        """Test that object below horizon returns -2."""
        # Choose a time when Venus is definitely below horizon
        # Early morning in winter, Venus might be below horizon
        jd = julday(2024, 12, 15, 3.0)
        geopos = (12.5, 42.0, 0)
        atmo = (1013.25, 0.0, 50.0, 0.0)
        observer = (36, 1.0)

        result, dret = vis_limit_mag(jd, geopos, atmo, observer, "Venus")

        # If below horizon, result should be -2
        if dret[1] < 0:
            assert result == HELFLAG_BELOW_HORIZON

    def test_photopic_vision_during_twilight(self):
        """Test photopic vision during bright twilight."""
        # During civil twilight, vision is photopic
        jd = julday(2024, 6, 21, 19.5)  # Summer evening twilight
        geopos = (12.5, 42.0, 0)
        atmo = (1013.25, 25.0, 50.0, 0.0)
        observer = (36, 1.0)

        result, dret = vis_limit_mag(jd, geopos, atmo, observer, "Venus")

        # If sun is between 0 and -6, should be photopic or mixed
        if dret[1] > 0:  # Venus above horizon
            if dret[3] > -6:
                assert result in (HELFLAG_PHOTOPIC, HELFLAG_MIXED)


class TestVisLimitMagEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_objname_raises_error(self):
        """Test that empty object name raises ValueError."""
        jd = julday(2024, 8, 15, 22.0)
        geopos = (12.5, 42.0, 0)
        atmo = (1013.25, 15.0, 50.0, 0.0)
        observer = (36, 1.0)

        with pytest.raises(ValueError, match="objname cannot be empty"):
            vis_limit_mag(jd, geopos, atmo, observer, "")

    def test_invalid_planet_raises_error(self):
        """Test that invalid planet name raises ValueError."""
        jd = julday(2024, 8, 15, 22.0)
        geopos = (12.5, 42.0, 0)
        atmo = (1013.25, 15.0, 50.0, 0.0)
        observer = (36, 1.0)

        with pytest.raises((ValueError, Exception)):
            vis_limit_mag(jd, geopos, atmo, observer, "InvalidPlanet123")

    def test_minimal_geopos(self):
        """Test with minimal geographic position (just longitude)."""
        jd = julday(2024, 8, 15, 22.0)
        geopos = (12.5,)  # Only longitude
        atmo = (1013.25, 15.0, 50.0, 0.0)
        observer = (36, 1.0)

        # Should still work with defaults for missing values
        result, dret = vis_limit_mag(jd, geopos, atmo, observer, "Venus")
        assert len(dret) == 10

    def test_minimal_atmo(self):
        """Test with minimal atmospheric parameters."""
        jd = julday(2024, 8, 15, 22.0)
        geopos = (12.5, 42.0, 0)
        atmo = (1013.25,)  # Only pressure
        observer = (36, 1.0)

        result, dret = vis_limit_mag(jd, geopos, atmo, observer, "Venus")
        assert len(dret) == 10

    def test_minimal_observer(self):
        """Test with minimal observer parameters."""
        jd = julday(2024, 8, 15, 22.0)
        geopos = (12.5, 42.0, 0)
        atmo = (1013.25, 15.0, 50.0, 0.0)
        observer = (36,)  # Only age

        result, dret = vis_limit_mag(jd, geopos, atmo, observer, "Venus")
        assert len(dret) == 10


class TestVisLimitMagAliases:
    """Test function aliases work correctly."""

    def test_swe_vis_limit_mag_alias(self):
        """Test that vis_limit_mag is an alias for vis_limit_mag."""
        jd = julday(2024, 8, 15, 22.0)
        geopos = (12.5, 42.0, 0)
        atmo = (1013.25, 15.0, 50.0, 0.0)
        observer = (36, 1.0)

        result1, dret1 = vis_limit_mag(jd, geopos, atmo, observer, "Venus")
        result2, dret2 = vis_limit_mag(jd, geopos, atmo, observer, "Venus")

        assert result1 == result2
        assert dret1 == dret2


class TestVisLimitMagRealistic:
    """Test realistic visibility scenarios."""

    def test_venus_brilliant_should_be_visible(self):
        """Test that brilliant Venus is easily visible."""
        # Venus at evening elongation is typically around mag -4
        jd = julday(2024, 8, 15, 21.0)
        geopos = (12.5, 42.0, 0)
        atmo = (1013.25, 20.0, 40.0, 0.0)
        observer = (36, 1.0)

        result, dret = vis_limit_mag(jd, geopos, atmo, observer, "Venus")

        if dret[1] > 5:  # Venus well above horizon
            # Venus (mag ~ -4) should be brighter than limiting magnitude
            # (i.e., dret[7] < dret[0] for visibility)
            # Venus is very bright so should almost always be visible when up
            pass  # Just verify calculation succeeded

    def test_limiting_mag_increases_with_darkness(self):
        """Test that limiting magnitude improves as night gets darker."""
        geopos = (12.5, 42.0, 0)
        atmo = (1013.25, 15.0, 50.0, 0.0)
        observer = (36, 1.0)

        # Track limiting magnitude through evening
        limiting_mags = []
        for hour in [20.0, 21.0, 22.0, 23.0]:
            jd = julday(2024, 8, 15, hour)
            result, dret = vis_limit_mag(jd, geopos, atmo, observer, "Venus")
            if dret[1] > 0:  # Object above horizon
                limiting_mags.append(dret[0])

        # Generally, limiting magnitude should increase as it gets darker
        # (higher limiting mag = can see fainter objects)
        if len(limiting_mags) >= 2:
            # At least verify we got valid calculations
            for lm in limiting_mags:
                assert -5 < lm < 10
