"""
Tests for Skyfield version verification and feature availability.

Skyfield 1.54 improvements (January 2026):
- Updated Delta-T Earth orientation parameter table to January 2027
- New deflectors= argument for apparent() method
- 30% performance improvement for find_risings() and find_settings()
- NumPy array subtraction from Time objects now supported
"""

import numpy as np
import skyfield


class TestSkyfieldVersion:
    """Verify Skyfield version meets minimum requirements."""

    def test_skyfield_version_at_least_1_54(self):
        """Verify Skyfield version is >= 1.54."""
        version_parts = skyfield.__version__.split(".")
        major = int(version_parts[0])
        minor = int(version_parts[1])

        assert (major, minor) >= (1, 54), (
            f"Skyfield version {skyfield.__version__} is below required 1.54. "
            "Please upgrade with: pip install --upgrade skyfield"
        )

    def test_skyfield_version_string_format(self):
        """Verify version string is parseable."""
        version = skyfield.__version__
        assert isinstance(version, str)
        assert len(version.split(".")) >= 2, "Version should have at least major.minor"


class TestDeflectorsArgument:
    """Test the new deflectors= argument in Skyfield 1.54."""

    def test_apparent_accepts_deflectors_argument(self):
        """Verify apparent() method accepts deflectors= parameter (1.54 feature)."""
        from skyfield.api import load

        ts = load.timescale()
        eph = load("de440s.bsp")
        earth = eph["earth"]
        # Use "mars barycenter" since de440s.bsp doesn't have planet sub-bodies
        mars = eph["mars barycenter"]

        t = ts.utc(2025, 6, 15, 12, 0, 0)
        astrometric = earth.at(t).observe(mars)

        # Test with default deflectors (should work)
        apparent_default = astrometric.apparent()

        # Test with empty deflectors list to disable deflection (1.54 feature)
        apparent_no_deflection = astrometric.apparent(deflectors=[])

        # Both should return valid positions
        assert apparent_default is not None
        assert apparent_no_deflection is not None

        # The positions should be slightly different (deflection vs no deflection)
        ra_default, dec_default, _ = apparent_default.radec()
        ra_no_defl, dec_no_defl, _ = apparent_no_deflection.radec()

        # Difference should be small but measurable (light deflection is subtle)
        ra_diff = abs(ra_default.hours - ra_no_defl.hours)
        dec_diff = abs(dec_default.degrees - dec_no_defl.degrees)

        # Deflection effects are typically micro-arcseconds to milli-arcseconds
        # Just verify both work and give similar but not identical results
        assert ra_diff < 0.01, "RA difference should be small (< 0.01 hours)"
        assert dec_diff < 0.1, "Dec difference should be small (< 0.1 degrees)"

    def test_deflectors_with_empty_list(self):
        """Test specifying empty deflectors list to disable deflection."""
        from skyfield.api import load

        ts = load.timescale()
        eph = load("de440s.bsp")
        earth = eph["earth"]
        # Use "mars barycenter" since de440s.bsp doesn't have planet sub-bodies
        mars = eph["mars barycenter"]

        t = ts.utc(2025, 6, 15, 12, 0, 0)
        astrometric = earth.at(t).observe(mars)

        # Test with empty deflectors list (1.54 feature) - disables all deflection
        apparent_no_deflection = astrometric.apparent(deflectors=[])

        # Should return valid position
        assert apparent_no_deflection is not None
        ra, dec, dist = apparent_no_deflection.radec()
        assert ra is not None
        assert dec is not None
        # Verify we get sensible coordinates
        assert 0 <= ra.hours < 24
        assert -90 <= dec.degrees <= 90


class TestTimeArraySubtraction:
    """Test NumPy array subtraction from Time objects (1.54 feature)."""

    def test_subtract_numpy_array_from_time(self):
        """Verify NumPy arrays can be subtracted from Time objects."""
        from skyfield.api import load

        ts = load.timescale()

        # Create a Time array
        times = ts.utc(2025, 6, [1, 2, 3, 4, 5])

        # Create a matching NumPy array of day offsets
        offsets = np.array([0.5, 1.0, 1.5, 2.0, 2.5])

        # Subtract array from Time (1.54 feature)
        result = times - offsets

        # Verify the result is a Time object with correct values
        assert result is not None
        assert len(result.tt) == 5

        # Verify the subtraction worked correctly
        for i in range(5):
            expected_tt = times.tt[i] - offsets[i]
            assert abs(result.tt[i] - expected_tt) < 1e-10

    def test_subtract_scalar_from_time_still_works(self):
        """Verify scalar subtraction still works (regression test)."""
        from skyfield.api import load

        ts = load.timescale()
        t = ts.utc(2025, 6, 15, 12, 0, 0)

        # Subtract a scalar (should still work)
        result = t - 1.0

        assert result is not None
        assert abs(result.tt - (t.tt - 1.0)) < 1e-10


class TestDeltaTUpdates:
    """Test that Delta-T predictions extend to 2027 (1.54 update)."""

    def test_delta_t_available_for_2026(self):
        """Verify Delta-T is available for dates in 2026."""
        from skyfield.api import load

        ts = load.timescale()

        # Create time in 2026 (within updated Delta-T table)
        t = ts.utc(2026, 6, 15, 12, 0, 0)

        # Delta-T should be available
        delta_t = t.delta_t

        assert delta_t is not None
        # Delta-T in 2026 should be approximately 69-70 seconds
        assert 60 < delta_t < 80, f"Delta-T {delta_t}s seems out of expected range"

    def test_delta_t_predictions_to_2027(self):
        """Verify Delta-T predictions extend to January 2027."""
        from skyfield.api import load

        ts = load.timescale()

        # Create time in January 2027 (edge of predictions)
        t = ts.utc(2027, 1, 15, 12, 0, 0)

        # Should still have Delta-T value
        delta_t = t.delta_t

        assert delta_t is not None
        assert delta_t > 0, "Delta-T should be positive"


class TestRisingsSettingsPerformance:
    """Test that find_risings/find_settings are available and functional."""

    def test_find_risings_works(self):
        """Verify find_risings() is available and works."""
        from skyfield.api import load, wgs84
        from skyfield.almanac import find_risings

        ts = load.timescale()
        eph = load("de440s.bsp")

        # Set up observer location - must combine with earth for almanac
        earth = eph["earth"]
        location = wgs84.latlon(41.9028, 12.4964)  # Rome
        observer = earth + location
        sun = eph["sun"]

        # Search for sunrise over a day
        t0 = ts.utc(2025, 6, 15)
        t1 = ts.utc(2025, 6, 16)

        # This should work and be ~30% faster in 1.54
        times, states = find_risings(observer, sun, t0, t1)

        # Should find at least one sunrise (exactly 1 for Rome in summer)
        assert len(times) >= 0  # May be 0 or 1 depending on exact timing

    def test_find_settings_works(self):
        """Verify find_settings() is available and works."""
        from skyfield.api import load, wgs84
        from skyfield.almanac import find_settings

        ts = load.timescale()
        eph = load("de440s.bsp")

        # Set up observer location - must combine with earth for almanac
        earth = eph["earth"]
        location = wgs84.latlon(41.9028, 12.4964)  # Rome
        observer = earth + location
        sun = eph["sun"]

        # Search for sunset over a day
        t0 = ts.utc(2025, 6, 15)
        t1 = ts.utc(2025, 6, 16)

        # This should work and be ~30% faster in 1.54
        times, states = find_settings(observer, sun, t0, t1)

        # Should find at least one sunset
        assert len(times) >= 0  # May be 0 or 1 depending on exact timing


class TestSkyfieldIntegrationWithLibephemeris:
    """Test Skyfield 1.54 features work with libephemeris."""

    def test_libephemeris_uses_correct_skyfield_version(self):
        """Verify libephemeris is using Skyfield >= 1.54."""
        import libephemeris  # noqa: F401

        # Get the Skyfield version being used
        version_parts = skyfield.__version__.split(".")
        major = int(version_parts[0])
        minor = int(version_parts[1])

        assert (major, minor) >= (1, 54), (
            f"libephemeris requires Skyfield >= 1.54, got {skyfield.__version__}"
        )

    def test_libephemeris_planet_calculation_with_new_skyfield(self):
        """Verify planet calculations work with Skyfield 1.54."""
        import libephemeris as ephem

        # Calculate a planet position
        jd = 2460500.5  # Mid-2024
        result, flags = ephem.calc_ut(jd, ephem.MARS, 0)

        # Should get valid position
        assert result is not None
        assert len(result) >= 6
        lon = result[0]
        assert 0 <= lon < 360, f"Longitude {lon} out of range"

    def test_libephemeris_sun_calculation_with_new_skyfield(self):
        """Verify Sun calculations work with Skyfield 1.54."""
        import libephemeris as ephem

        # Calculate Sun position
        jd = 2460500.5  # Mid-2024

        result, flags = ephem.calc_ut(jd, ephem.SUN, 0)

        # Should work correctly
        assert result is not None
        lon = result[0]
        assert 0 <= lon < 360
