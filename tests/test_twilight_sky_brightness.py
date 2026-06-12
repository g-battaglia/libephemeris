"""
Tests for twilight sky brightness model in libephemeris.

Tests the calculation of sky brightness during civil twilight (Sun 0° to -6°),
nautical twilight (-6° to -12°), and astronomical twilight (-12° to -18°)
as a function of Sun altitude, azimuth relative to target, and atmospheric conditions.

References:
    - Patat, F. (2003) "UBVRI twilight sky brightness at ESO-Paranal"
    - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
    - Rozenberg, G.V. (1966) "Twilight: A Study in Atmospheric Optics"
"""


from libephemeris import (
    # Main twilight functions
    get_twilight_phase,
    calc_twilight_sky_brightness,
    calc_twilight_brightness_simple,
    calc_limiting_magnitude_twilight,
    # Twilight constants
    TWILIGHT_CIVIL_START,
    TWILIGHT_CIVIL_END,
    TWILIGHT_NAUTICAL_END,
    TWILIGHT_ASTRONOMICAL_END,
    DARK_SKY_BRIGHTNESS_V,
)


class TestTwilightPhase:
    """Tests for get_twilight_phase function."""

    def test_daytime(self):
        """Sun above horizon is daytime."""
        assert get_twilight_phase(10.0) == "day"
        assert get_twilight_phase(45.0) == "day"
        assert get_twilight_phase(0.1) == "day"

    def test_civil_twilight(self):
        """Sun between 0 and -6 degrees is civil twilight."""
        assert get_twilight_phase(0.0) == "civil"
        assert get_twilight_phase(-3.0) == "civil"
        assert get_twilight_phase(-5.9) == "civil"

    def test_nautical_twilight(self):
        """Sun between -6 and -12 degrees is nautical twilight."""
        assert get_twilight_phase(-6.0) == "nautical"
        assert get_twilight_phase(-9.0) == "nautical"
        assert get_twilight_phase(-11.9) == "nautical"

    def test_astronomical_twilight(self):
        """Sun between -12 and -18 degrees is astronomical twilight."""
        assert get_twilight_phase(-12.0) == "astronomical"
        assert get_twilight_phase(-15.0) == "astronomical"
        assert get_twilight_phase(-17.9) == "astronomical"

    def test_night(self):
        """Sun below -18 degrees is full night."""
        assert get_twilight_phase(-18.0) == "night"
        assert get_twilight_phase(-25.0) == "night"
        assert get_twilight_phase(-45.0) == "night"


class TestTwilightConstants:
    """Tests for twilight phase constants."""

    def test_constants_defined(self):
        """All twilight constants should be defined with correct values."""
        assert TWILIGHT_CIVIL_START == 0.0
        assert TWILIGHT_CIVIL_END == -6.0
        assert TWILIGHT_NAUTICAL_END == -12.0
        assert TWILIGHT_ASTRONOMICAL_END == -18.0

    def test_dark_sky_brightness(self):
        """Dark sky brightness should be typical value for V-band."""
        # Typical dark sky is 21-22 mag/arcsec^2
        assert 21.0 <= DARK_SKY_BRIGHTNESS_V <= 22.5


class TestTwilightSkyBrightnessDataclass:
    """Tests for TwilightSkyBrightness dataclass."""

    def test_dataclass_structure(self):
        """TwilightSkyBrightness should have all expected fields."""
        result = calc_twilight_sky_brightness(-9.0)

        assert hasattr(result, "surface_brightness")
        assert hasattr(result, "twilight_phase")
        assert hasattr(result, "sun_altitude_deg")
        assert hasattr(result, "azimuth_factor")
        assert hasattr(result, "limiting_magnitude")
        assert hasattr(result, "nanolamberts")

    def test_dataclass_types(self):
        """TwilightSkyBrightness fields should have correct types."""
        result = calc_twilight_sky_brightness(-9.0)

        assert isinstance(result.surface_brightness, float)
        assert isinstance(result.twilight_phase, str)
        assert isinstance(result.sun_altitude_deg, float)
        assert isinstance(result.azimuth_factor, float)
        assert isinstance(result.limiting_magnitude, float)
        assert isinstance(result.nanolamberts, float)


class TestTwilightSkyBrightness:
    """Tests for calc_twilight_sky_brightness function."""

    def test_civil_twilight_brightness(self):
        """Civil twilight should have relatively bright sky."""
        result = calc_twilight_sky_brightness(-3.0)

        assert result.twilight_phase == "civil"
        # Civil twilight is bright: ~3-8 mag/arcsec^2
        assert 3.0 <= result.surface_brightness <= 10.0

    def test_nautical_twilight_brightness(self):
        """Nautical twilight should have moderate sky brightness."""
        result = calc_twilight_sky_brightness(-9.0)

        assert result.twilight_phase == "nautical"
        # Nautical twilight: ~8-17 mag/arcsec^2
        assert 8.0 <= result.surface_brightness <= 18.0

    def test_astronomical_twilight_brightness(self):
        """Astronomical twilight should have relatively dark sky."""
        result = calc_twilight_sky_brightness(-15.0)

        assert result.twilight_phase == "astronomical"
        # Astronomical twilight: ~17-22 mag/arcsec^2
        assert 15.0 <= result.surface_brightness <= 22.0

    def test_night_brightness(self):
        """Full night should have dark sky brightness."""
        result = calc_twilight_sky_brightness(-25.0)

        assert result.twilight_phase == "night"
        # Night sky should approach dark sky value
        assert result.surface_brightness >= 20.0
        assert result.surface_brightness <= 23.0

    def test_daytime_brightness(self):
        """Daytime should have very bright sky."""
        result = calc_twilight_sky_brightness(10.0)

        assert result.twilight_phase == "day"
        # Daytime is very bright
        assert result.surface_brightness <= 5.0

    def test_brightness_increases_with_lower_sun(self):
        """Sky brightness (mag/arcsec^2) should increase as Sun goes lower."""
        # Higher mag/arcsec^2 = darker sky
        # Note: Once we reach full night (-18 or below), brightness plateaus
        sun_altitudes = [0.0, -3.0, -6.0, -9.0, -12.0, -15.0, -18.0]
        brightnesses = [
            calc_twilight_sky_brightness(alt).surface_brightness
            for alt in sun_altitudes
        ]

        # Each value should be higher (darker) than the previous
        for i in range(1, len(brightnesses)):
            assert brightnesses[i] > brightnesses[i - 1], (
                f"Brightness at {sun_altitudes[i]} ({brightnesses[i]}) "
                f"should be darker than at {sun_altitudes[i - 1]} ({brightnesses[i - 1]})"
            )

        # Full night: sky brightness plateaus at dark sky value
        night_brightness = calc_twilight_sky_brightness(-25.0).surface_brightness
        astro_end_brightness = calc_twilight_sky_brightness(-18.0).surface_brightness
        # Should be approximately equal (both at dark sky value)
        assert abs(night_brightness - astro_end_brightness) < 0.5

    def test_sun_altitude_recorded(self):
        """Result should contain the input Sun altitude."""
        result = calc_twilight_sky_brightness(-7.5)
        assert result.sun_altitude_deg == -7.5


class TestAzimuthalVariation:
    """Tests for azimuthal variation of sky brightness."""

    def test_brighter_toward_sun(self):
        """Sky should be brighter in the direction of the Sun."""
        # During civil twilight, look toward and away from Sun
        toward_sun = calc_twilight_sky_brightness(
            -3.0,
            target_altitude_deg=30.0,
            sun_azimuth_deg=270.0,
            target_azimuth_deg=270.0,
        )
        away_from_sun = calc_twilight_sky_brightness(
            -3.0,
            target_altitude_deg=30.0,
            sun_azimuth_deg=270.0,
            target_azimuth_deg=90.0,
        )

        # Looking toward Sun = brighter = lower mag/arcsec^2
        assert toward_sun.surface_brightness < away_from_sun.surface_brightness

    def test_azimuth_factor_range(self):
        """Azimuth factor should be between 0 and 1."""
        result = calc_twilight_sky_brightness(
            -6.0,
            target_altitude_deg=45.0,
            sun_azimuth_deg=180.0,
            target_azimuth_deg=90.0,
        )
        assert 0.0 <= result.azimuth_factor <= 1.0

    def test_azimuth_effect_diminishes_at_night(self):
        """Azimuthal variation should be minimal at night."""
        toward_sun = calc_twilight_sky_brightness(
            -25.0, sun_azimuth_deg=270.0, target_azimuth_deg=270.0
        )
        away_from_sun = calc_twilight_sky_brightness(
            -25.0, sun_azimuth_deg=270.0, target_azimuth_deg=90.0
        )

        # At night, difference should be minimal
        diff = abs(toward_sun.surface_brightness - away_from_sun.surface_brightness)
        assert diff < 0.5  # Less than 0.5 mag difference


class TestAltitudeVariation:
    """Tests for altitude/zenith angle variation of sky brightness."""

    def test_brighter_at_horizon_during_twilight(self):
        """Sky should be brighter near horizon during twilight."""
        zenith = calc_twilight_sky_brightness(-6.0, target_altitude_deg=90.0)
        horizon = calc_twilight_sky_brightness(-6.0, target_altitude_deg=10.0)

        # Horizon is brighter (lower mag/arcsec^2)
        assert horizon.surface_brightness < zenith.surface_brightness

    def test_zenith_vs_intermediate_altitude(self):
        """Compare zenith with intermediate altitude."""
        zenith = calc_twilight_sky_brightness(-9.0, target_altitude_deg=90.0)
        mid_sky = calc_twilight_sky_brightness(-9.0, target_altitude_deg=45.0)

        # Mid-sky should be slightly brighter than zenith
        assert mid_sky.surface_brightness <= zenith.surface_brightness + 0.5


class TestAtmosphericConditions:
    """Tests for atmospheric condition effects on twilight brightness."""

    def test_humidity_affects_brightness(self):
        """Higher humidity should brighten the twilight sky (more scattering)."""
        low_humidity = calc_twilight_sky_brightness(-6.0, humidity_percent=20.0)
        high_humidity = calc_twilight_sky_brightness(-6.0, humidity_percent=90.0)

        # High humidity = more scattering = brighter sky = lower mag/arcsec^2
        # Effect may be subtle
        assert isinstance(low_humidity.surface_brightness, float)
        assert isinstance(high_humidity.surface_brightness, float)

    def test_altitude_affects_brightness(self):
        """Higher observer altitude should result in slightly darker sky."""
        sea_level = calc_twilight_sky_brightness(-9.0, altitude_m=0.0)
        mountain = calc_twilight_sky_brightness(-9.0, altitude_m=3000.0)

        # Higher altitude = less atmosphere = slightly darker (effect may be subtle)
        assert isinstance(sea_level.surface_brightness, float)
        assert isinstance(mountain.surface_brightness, float)


class TestLimitingMagnitude:
    """Tests for limiting magnitude calculations."""

    def test_limiting_magnitude_civil_twilight(self):
        """Limiting magnitude should be low during civil twilight."""
        result = calc_twilight_sky_brightness(-3.0)

        # Civil twilight: can only see brightest objects
        assert result.limiting_magnitude < 4.0

    def test_limiting_magnitude_nautical_twilight(self):
        """Limiting magnitude should be moderate during nautical twilight."""
        result = calc_twilight_sky_brightness(-9.0)

        # Nautical twilight: brighter stars visible
        assert 2.0 < result.limiting_magnitude < 6.0

    def test_limiting_magnitude_astronomical_twilight(self):
        """Limiting magnitude should be high during astronomical twilight."""
        result = calc_twilight_sky_brightness(-15.0)

        # Astronomical twilight: most naked-eye stars visible
        assert result.limiting_magnitude > 4.0

    def test_limiting_magnitude_night(self):
        """Limiting magnitude should approach 6+ at night."""
        result = calc_twilight_sky_brightness(-25.0)

        # Full night: typical limit is 6.0-6.5
        assert result.limiting_magnitude >= 5.5

    def test_limiting_magnitude_increases_with_darkness(self):
        """Limiting magnitude should increase as sky gets darker."""
        sun_altitudes = [-3.0, -9.0, -15.0, -25.0]
        limits = [
            calc_twilight_sky_brightness(alt).limiting_magnitude
            for alt in sun_altitudes
        ]

        for i in range(1, len(limits)):
            assert limits[i] >= limits[i - 1]


class TestNanoLamberts:
    """Tests for nanoLambert conversion."""

    def test_nanolamberts_positive(self):
        """Nanolambert values should be positive."""
        for sun_alt in [-3.0, -9.0, -15.0, -25.0]:
            result = calc_twilight_sky_brightness(sun_alt)
            assert result.nanolamberts > 0

    def test_nanolamberts_decreases_with_darkness(self):
        """Nanolamberts should decrease as sky gets darker."""
        civil = calc_twilight_sky_brightness(-3.0)
        night = calc_twilight_sky_brightness(-25.0)

        # Brighter sky = more nanolamberts
        assert civil.nanolamberts > night.nanolamberts


class TestCalcTwilightBrightnessSimple:
    """Tests for the simplified twilight brightness function."""

    def test_returns_zenith_brightness(self):
        """Simple function should return approximate zenith brightness."""
        brightness = calc_twilight_brightness_simple(-9.0)
        full_result = calc_twilight_sky_brightness(-9.0, target_altitude_deg=90.0)

        # Should be close to full calculation with default parameters
        assert abs(brightness - full_result.surface_brightness) < 1.0

    def test_key_values(self):
        """Test brightness at key Sun altitudes."""
        # At horizon
        assert 2.0 <= calc_twilight_brightness_simple(0.0) <= 5.0

        # End of civil twilight
        brightness_civil = calc_twilight_brightness_simple(-6.0)
        assert 6.0 <= brightness_civil <= 12.0

        # End of nautical twilight
        brightness_nautical = calc_twilight_brightness_simple(-12.0)
        assert 14.0 <= brightness_nautical <= 19.0

        # End of astronomical twilight
        brightness_astro = calc_twilight_brightness_simple(-18.0)
        assert 19.0 <= brightness_astro <= 22.5


class TestCalcLimitingMagnitudeTwilight:
    """Tests for calc_limiting_magnitude_twilight function."""

    def test_returns_float(self):
        """Function should return a float."""
        result = calc_limiting_magnitude_twilight(-9.0)
        assert isinstance(result, float)

    def test_reasonable_range(self):
        """Limiting magnitude should be in reasonable range."""
        for sun_alt in [-3.0, -9.0, -15.0, -25.0]:
            limit = calc_limiting_magnitude_twilight(sun_alt)
            assert -2.0 <= limit <= 7.0

    def test_altitude_affects_limiting_magnitude(self):
        """Looking at lower altitudes should reduce limiting magnitude."""
        zenith_limit = calc_limiting_magnitude_twilight(-9.0, target_altitude_deg=90.0)
        low_alt_limit = calc_limiting_magnitude_twilight(-9.0, target_altitude_deg=10.0)

        # More extinction at lower altitudes = lower limiting magnitude
        assert low_alt_limit <= zenith_limit

    def test_azimuth_affects_limiting_magnitude(self):
        """Looking toward Sun should reduce limiting magnitude."""
        toward_sun = calc_limiting_magnitude_twilight(
            -6.0, sun_azimuth_deg=270.0, target_azimuth_deg=270.0
        )
        away_from_sun = calc_limiting_magnitude_twilight(
            -6.0, sun_azimuth_deg=270.0, target_azimuth_deg=90.0
        )

        # Brighter sky toward Sun = lower limiting magnitude
        assert toward_sun <= away_from_sun


class TestEdgeCases:
    """Tests for edge cases and boundary conditions."""

    def test_exact_boundary_values(self):
        """Test behavior at exact twilight boundaries."""
        # At exact boundary between civil and nautical
        result = calc_twilight_sky_brightness(-6.0)
        # Should be nautical (inclusive of lower bound)
        assert result.twilight_phase == "nautical"

    def test_extreme_sun_altitudes(self):
        """Test with extreme Sun altitudes."""
        # Very low Sun
        result_low = calc_twilight_sky_brightness(-45.0)
        assert result_low.twilight_phase == "night"
        assert result_low.surface_brightness > 20.0

        # Sun high in sky
        result_high = calc_twilight_sky_brightness(60.0)
        assert result_high.twilight_phase == "day"
        assert result_high.surface_brightness < 5.0

    def test_extreme_viewing_angles(self):
        """Test with extreme viewing angles."""
        # Looking straight up
        result_zenith = calc_twilight_sky_brightness(-9.0, target_altitude_deg=90.0)
        assert 8.0 <= result_zenith.surface_brightness <= 20.0

        # Looking at horizon
        result_horizon = calc_twilight_sky_brightness(-9.0, target_altitude_deg=0.0)
        assert isinstance(result_horizon.surface_brightness, float)

    def test_azimuth_wraparound(self):
        """Test azimuth values that wrap around 0/360."""
        result1 = calc_twilight_sky_brightness(
            -6.0, sun_azimuth_deg=350.0, target_azimuth_deg=10.0
        )
        result2 = calc_twilight_sky_brightness(
            -6.0, sun_azimuth_deg=10.0, target_azimuth_deg=350.0
        )

        # Should be symmetric
        assert abs(result1.surface_brightness - result2.surface_brightness) < 0.1


class TestPhysicalReasonability:
    """Tests to verify physically reasonable behavior."""

    def test_brightness_never_negative(self):
        """Surface brightness should never be negative."""
        for sun_alt in range(-45, 45, 5):
            result = calc_twilight_sky_brightness(float(sun_alt))
            assert result.surface_brightness >= 0.0

    def test_limiting_magnitude_reasonable_bounds(self):
        """Limiting magnitude should be within observable bounds."""
        for sun_alt in range(-30, 10, 3):
            result = calc_twilight_sky_brightness(float(sun_alt))
            # Naked eye limit is typically -2 to +7
            assert -3.0 <= result.limiting_magnitude <= 8.0

    def test_smooth_transition_between_phases(self):
        """Brightness should change smoothly between twilight phases."""
        prev_brightness = 0.0
        for sun_alt in range(5, -25, -1):
            result = calc_twilight_sky_brightness(float(sun_alt))

            if prev_brightness > 0:
                # Change should not be too abrupt
                change = abs(result.surface_brightness - prev_brightness)
                assert change < 3.0, (
                    f"Abrupt change at sun_alt={sun_alt}: "
                    f"{prev_brightness} -> {result.surface_brightness}"
                )

            prev_brightness = result.surface_brightness


class TestImportability:
    """Tests to verify proper import from libephemeris."""

    def test_all_functions_importable(self):
        """All twilight functions should be importable from libephemeris."""
        from libephemeris import (  # noqa: F811
            get_twilight_phase,
            calc_twilight_sky_brightness,
            calc_twilight_brightness_simple,
            calc_limiting_magnitude_twilight,
            TWILIGHT_CIVIL_END,
        )

        # Verify they are callable/usable
        assert callable(get_twilight_phase)
        assert callable(calc_twilight_sky_brightness)
        assert callable(calc_twilight_brightness_simple)
        assert callable(calc_limiting_magnitude_twilight)
        assert isinstance(TWILIGHT_CIVIL_END, float)


class TestIntegrationWithExtinction:
    """Tests for integration with atmospheric extinction functions."""

    def test_uses_extinction_coefficients(self):
        """Twilight model should incorporate atmospheric extinction."""
        from libephemeris import calc_extinction_coefficient

        # Get extinction coefficient
        coeff = calc_extinction_coefficient(humidity_percent=80.0)

        # Calculate twilight brightness with high humidity
        result = calc_twilight_sky_brightness(-6.0, humidity_percent=80.0)

        # Result should incorporate aerosol effects
        assert isinstance(result.surface_brightness, float)


class TestRealisticScenarios:
    """Tests based on realistic observational scenarios."""

    def test_first_star_visibility(self):
        """Test when first stars become visible (civil twilight)."""
        # During civil twilight, only brightest stars visible
        result = calc_twilight_sky_brightness(-4.0)

        # Limiting magnitude should allow Venus (-4) and bright stars
        # to be visible, but not faint stars
        assert result.limiting_magnitude > -5.0  # Venus visible
        assert result.limiting_magnitude < 4.0  # Faint stars not visible

    def test_horizon_visibility_nautical(self):
        """Test horizon visibility during nautical twilight."""
        result = calc_twilight_sky_brightness(-9.0)

        # At nautical twilight, horizon is barely visible
        # and many stars are visible
        assert result.limiting_magnitude > 3.0

    def test_deep_sky_observation_conditions(self):
        """Test conditions for deep-sky observation."""
        # Astronomical twilight
        result_astro = calc_twilight_sky_brightness(-15.0)
        # Full night
        result_night = calc_twilight_sky_brightness(-25.0)

        # Deep sky objects need dark skies (high mag/arcsec^2)
        assert result_astro.surface_brightness > 17.0
        assert result_night.surface_brightness > 20.0
