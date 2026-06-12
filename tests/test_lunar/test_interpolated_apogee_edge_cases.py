"""
Tests for edge cases in the interpolated lunar apogee calculation.

Edge cases tested:
1. 360/0 degree wrapping: When the apogee crosses the 0/360 boundary during
   the interpolation window, the unwrapping logic must handle it correctly.

2. Ephemeris boundary dates: When the target date is within 28 days of the
   ephemeris boundary, the sampling window must be adjusted to stay within
   the valid ephemeris range.
"""

from libephemeris import lunar


class TestLongitudeUnwrapping:
    """Tests for the _unwrap_longitudes helper function."""

    def test_unwrap_no_discontinuity(self):
        """Test unwrapping when there is no discontinuity."""
        longitudes = [10.0, 15.0, 20.0, 25.0, 30.0]
        result = lunar._unwrap_longitudes(longitudes)

        assert len(result) == 5
        assert result[0] == 10.0
        assert result[1] == 15.0
        assert result[2] == 20.0
        assert result[3] == 25.0
        assert result[4] == 30.0

    def test_unwrap_crossing_zero_forward(self):
        """Test unwrapping when crossing 0 degrees moving forward."""
        # Crossing from ~350 to ~10 (forward motion crossing 0/360)
        longitudes = [340.0, 350.0, 0.0, 10.0, 20.0]
        result = lunar._unwrap_longitudes(longitudes)

        # Should unwrap to continuous values > 360
        assert len(result) == 5
        assert result[0] == 340.0
        assert result[1] == 350.0
        assert result[2] == 360.0  # 0 becomes 360
        assert result[3] == 370.0  # 10 becomes 370
        assert result[4] == 380.0  # 20 becomes 380

    def test_unwrap_crossing_zero_backward(self):
        """Test unwrapping when crossing 0 degrees moving backward."""
        # Crossing from ~10 to ~350 (backward motion crossing 0/360)
        longitudes = [20.0, 10.0, 0.0, 350.0, 340.0]
        result = lunar._unwrap_longitudes(longitudes)

        # Should unwrap to continuous values < 0
        assert len(result) == 5
        assert result[0] == 20.0
        assert result[1] == 10.0
        assert result[2] == 0.0
        assert result[3] == -10.0  # 350 becomes -10
        assert result[4] == -20.0  # 340 becomes -20

    def test_unwrap_empty_list(self):
        """Test unwrapping an empty list."""
        result = lunar._unwrap_longitudes([])
        assert result == []

    def test_unwrap_single_value(self):
        """Test unwrapping a single value."""
        result = lunar._unwrap_longitudes([180.0])
        assert result == [180.0]

    def test_unwrap_multiple_crossings(self):
        """Test unwrapping with multiple crossings (shouldn't happen normally)."""
        # This tests the algorithm handles multiple wraps correctly
        longitudes = [350.0, 10.0, 350.0, 10.0]  # Two back-and-forth crossings
        result = lunar._unwrap_longitudes(longitudes)

        # Should maintain continuity
        assert len(result) == 4
        # Each transition should be smooth (diff < 180)
        for i in range(1, len(result)):
            diff = result[i] - result[i - 1]
            assert abs(diff) < 180, f"Jump at index {i}: {diff}"


class TestEphemerisRange:
    """Tests for the _get_ephemeris_range helper function."""

    def test_ephemeris_range_returns_valid_values(self):
        """Test that ephemeris range returns sensible min/max values."""
        min_jd, max_jd = lunar._get_ephemeris_range()

        # min_jd should be before max_jd
        assert min_jd < max_jd

        # Range depends on the loaded ephemeris file:
        # - DE421: ~1899-2053 (min_jd ~2415020, max_jd ~2471184)
        # - DE440: ~1549-2650 (min_jd ~2287184, max_jd ~2688976)
        # - DE441: ~-13200-17191 (min_jd ~-3100255, max_jd ~7999032)
        assert min_jd > -4000000  # Reasonable lower bound (accommodates DE441)
        assert max_jd < 9000000  # Reasonable upper bound (accommodates DE441)


class TestInterpolatedApogeeAtBoundary:
    """Tests for interpolated apogee near ephemeris boundaries."""

    def test_apogee_near_start_of_ephemeris(self):
        """Test interpolated apogee calculation near start of ephemeris."""
        min_jd, _ = lunar._get_ephemeris_range()

        # Test at a date within 28 days of the start boundary
        # (where the full sampling window would extend before the start)
        jd_near_start = min_jd + 10.0  # 10 days after start

        # Should not raise an exception
        lon, lat, ecc = lunar.calc_interpolated_apogee(jd_near_start)

        # Results should be valid
        assert 0 <= lon < 360, f"Longitude {lon} out of range"
        assert -10 < lat < 10, f"Latitude {lat} seems unreasonable"
        assert 0.002 < ecc < 0.003, f"Distance {ecc} AU seems unreasonable"

    def test_apogee_near_end_of_ephemeris(self):
        """Test interpolated apogee calculation near end of ephemeris."""
        _, max_jd = lunar._get_ephemeris_range()

        # Test at a date within 28 days of the end boundary
        # (where the full sampling window would extend beyond the end)
        jd_near_end = max_jd - 10.0  # 10 days before end

        # Should not raise an exception
        lon, lat, ecc = lunar.calc_interpolated_apogee(jd_near_end)

        # Results should be valid
        assert 0 <= lon < 360, f"Longitude {lon} out of range"
        assert -10 < lat < 10, f"Latitude {lat} seems unreasonable"
        assert 0.002 < ecc < 0.003, f"Distance {ecc} AU seems unreasonable"

    def test_apogee_at_exact_start_boundary(self):
        """Test interpolated apogee at the exact start of ephemeris range."""
        min_jd, _ = lunar._get_ephemeris_range()

        # Test at exactly the start boundary (or just after)
        jd_at_start = min_jd + 1.0  # 1 day after start to ensure it's valid

        # Should not raise an exception
        lon, lat, ecc = lunar.calc_interpolated_apogee(jd_at_start)

        # Results should be valid
        assert 0 <= lon < 360, f"Longitude {lon} out of range"

    def test_apogee_at_exact_end_boundary(self):
        """Test interpolated apogee at the exact end of ephemeris range."""
        _, max_jd = lunar._get_ephemeris_range()

        # Test at exactly the end boundary (or just before)
        jd_at_end = max_jd - 1.0  # 1 day before end to ensure it's valid

        # Should not raise an exception
        lon, lat, ecc = lunar.calc_interpolated_apogee(jd_at_end)

        # Results should be valid
        assert 0 <= lon < 360, f"Longitude {lon} out of range"


class TestInterpolatedPerigeeAtBoundary:
    """Tests for interpolated perigee near ephemeris boundaries."""

    def test_perigee_near_start_of_ephemeris(self):
        """Test interpolated perigee calculation near start of ephemeris."""
        min_jd, _ = lunar._get_ephemeris_range()

        # Test at a date within 28 days of the start boundary
        jd_near_start = min_jd + 10.0  # 10 days after start

        # Should not raise an exception
        lon, lat, ecc = lunar.calc_interpolated_perigee(jd_near_start)

        # Results should be valid
        assert 0 <= lon < 360, f"Longitude {lon} out of range"
        assert -10 < lat < 10, f"Latitude {lat} seems unreasonable"
        assert 0.002 < ecc < 0.003, f"Distance {ecc} AU seems unreasonable"

    def test_perigee_near_end_of_ephemeris(self):
        """Test interpolated perigee calculation near end of ephemeris."""
        _, max_jd = lunar._get_ephemeris_range()

        # Test at a date within 28 days of the end boundary
        jd_near_end = max_jd - 10.0  # 10 days before end

        # Should not raise an exception
        lon, lat, ecc = lunar.calc_interpolated_perigee(jd_near_end)

        # Results should be valid
        assert 0 <= lon < 360, f"Longitude {lon} out of range"
        assert -10 < lat < 10, f"Latitude {lat} seems unreasonable"
        assert 0.002 < ecc < 0.003, f"Distance {ecc} AU seems unreasonable"


class TestLongitudeWrappingInInterpolation:
    """Tests for 360/0 degree wrapping during interpolation."""

    def test_apogee_crossing_zero_returns_valid_result(self):
        """
        Test that interpolated apogee handles crossing the 0/360 boundary.

        This test finds a date where the apogee is near 0/360 degrees and
        verifies the interpolation produces a valid result.
        """
        # J2000.0 as a reference point
        jd_tt = 2451545.0

        # Sample across many dates to find one where apogee is near 0/360
        for day_offset in range(0, 365, 1):  # Sample across a year
            jd = jd_tt + day_offset
            lon, lat, ecc = lunar.calc_interpolated_apogee(jd)

            # Check that the result is always valid
            assert 0 <= lon < 360, f"Invalid longitude {lon} at JD {jd}"
            assert -10 < lat < 10, f"Unreasonable latitude {lat} at JD {jd}"
            assert 0.002 < ecc < 0.003, f"Unreasonable distance {ecc} AU at JD {jd}"

            # If we're near 0 or 360, the algorithm handled the wrap correctly
            if lon < 5 or lon > 355:
                # Found a case near the boundary - this is the important test
                pass

    def test_perigee_crossing_zero_returns_valid_result(self):
        """
        Test that interpolated perigee handles crossing the 0/360 boundary.

        This test finds a date where the perigee is near 0/360 degrees and
        verifies the interpolation produces a valid result.
        """
        jd_tt = 2451545.0

        # Sample across many dates to find one where perigee is near 0/360
        for day_offset in range(0, 365, 1):
            jd = jd_tt + day_offset
            lon, lat, ecc = lunar.calc_interpolated_perigee(jd)

            # Check that the result is always valid
            assert 0 <= lon < 360, f"Invalid longitude {lon} at JD {jd}"
            assert -10 < lat < 10, f"Unreasonable latitude {lat} at JD {jd}"
            assert 0.002 < ecc < 0.003, f"Unreasonable distance {ecc} AU at JD {jd}"

    def test_continuity_across_zero_boundary(self):
        """
        Test that interpolated apogee moves continuously when crossing 0/360.

        The interpolation should produce smooth results even when the underlying
        osculating values cross the 0/360 boundary.
        """
        jd_tt = 2451545.0
        dt = 0.5  # Half-day steps for finer resolution

        prev_lon = None
        for i in range(100):
            jd = jd_tt + i * dt
            lon, _, _ = lunar.calc_interpolated_apogee(jd)

            if prev_lon is not None:
                diff = lon - prev_lon
                # Handle wrap-around for comparison
                if diff > 180:
                    diff -= 360
                elif diff < -180:
                    diff += 360

                # The change should be small (interpolated apogee moves slowly)
                # Allow up to 0.6 degrees per half-day (slightly relaxed for
                # the analytical ELP2000-82B perturbation series which can have
                # slightly larger variations due to multiple oscillating terms)
                assert abs(diff) < 0.6, (
                    f"Large jump of {diff} degrees at step {i} "
                    f"(from {prev_lon} to {lon})"
                )

            prev_lon = lon


class TestSamplingWithFallback:
    """Tests for the _sample_osculating_apogee_with_fallback function."""

    def test_fallback_returns_valid_samples(self):
        """Test that fallback sampling returns valid data."""
        jd_tt = 2451545.0  # J2000.0 - well within ephemeris range
        half_window = 28.0
        num_samples = 9

        times, lons, lats, eccs, target_idx = (
            lunar._sample_osculating_apogee_with_fallback(
                jd_tt, half_window, num_samples
            )
        )

        # Should return all 9 samples for a date well within range
        assert len(times) == 9
        assert len(lons) == 9
        assert len(lats) == 9
        assert len(eccs) == 9

        # Target index should be near the center
        assert target_idx == 4  # Central sample

        # All longitudes should be valid
        for lon in lons:
            assert 0 <= lon < 360

    def test_fallback_near_boundary_adjusts_window(self):
        """Test that sampling near boundary adjusts the window."""
        min_jd, _ = lunar._get_ephemeris_range()
        half_window = 28.0
        num_samples = 9

        # Test at a date where the window would extend before min_jd
        jd_near_start = min_jd + 10.0

        times, lons, lats, eccs, target_idx = (
            lunar._sample_osculating_apogee_with_fallback(
                jd_near_start, half_window, num_samples
            )
        )

        # Should return samples (may be fewer than 9 if window is constrained)
        assert len(times) >= 2  # At least 2 for regression
        assert len(times) == len(lons) == len(lats) == len(eccs)

        # All sample times should be within ephemeris range
        for t in times:
            assert t >= min_jd, f"Sample time {t} is before min_jd {min_jd}"


class TestInterpolatedApogeeConsistency:
    """Tests for consistency of interpolated apogee with boundary handling."""

    def test_boundary_result_reasonable_compared_to_osculating(self):
        """
        Test that boundary interpolation gives reasonable results.

        Near the boundary, the interpolated apogee should still be within
        a reasonable range of the osculating apogee.
        """
        min_jd, _ = lunar._get_ephemeris_range()
        jd_near_start = min_jd + 15.0

        # Get interpolated result
        interp_lon, interp_lat, interp_ecc = lunar.calc_interpolated_apogee(
            jd_near_start
        )

        # Get osculating result
        oscu_lon, oscu_lat, oscu_ecc = lunar.calc_true_lilith(jd_near_start)

        # Difference should be within reasonable range
        # The interpolation smooths out oscillations, so difference can be up to ~30°
        diff = abs(interp_lon - oscu_lon)
        if diff > 180:
            diff = 360 - diff

        assert diff < 40, (
            f"Interpolated ({interp_lon}) and osculating ({oscu_lon}) "
            f"differ by {diff} degrees"
        )

    def test_perigee_opposite_apogee_at_boundary(self):
        """
        Test that perigee is roughly opposite apogee near ephemeris boundary.

        Note: With the analytical ELP2000-82B implementation for apogee and the
        polynomial regression for perigee, there can be larger differences since
        they use different computational approaches. As established in lunar mechanics
        literature (Chapront-Touzé & Chapront 1988), apogee and perigee are NOT
        exactly opposite (can differ by up to ~28°).
        """
        min_jd, _ = lunar._get_ephemeris_range()
        jd_near_start = min_jd + 15.0

        apogee_lon, _, _ = lunar.calc_interpolated_apogee(jd_near_start)
        perigee_lon, _, _ = lunar.calc_interpolated_perigee(jd_near_start)

        # Calculate angular difference
        diff = abs(apogee_lon - perigee_lon)
        if diff > 180:
            diff = 360 - diff

        # Should be approximately 180 degrees (within 28° tolerance)
        # The tolerance is larger because:
        # 1. Apogee uses fitted ELP2000-82B perturbations with 2D-2M' dominant term
        # 2. Perigee uses independent fitted perturbations
        # 3. Up to 28° deviation is a known physical property of the Moon's orbit
        assert abs(diff - 180) < 28.0, (
            f"Apogee ({apogee_lon}) and perigee ({perigee_lon}) "
            f"differ by {diff} degrees (expected ~180)"
        )
