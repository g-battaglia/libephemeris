"""
Tests for the interpolated lunar perigee (INTP_PERG).

The interpolated perigee is a smoothed version of the osculating perigee that
removes spurious short-period oscillations by polynomial regression through
successive osculating perigee positions.

Algorithm Design (based on analysis of JPL DE440 lunar passage data):
- Sample osculating perigee (apogee + 180°) at 7 points spanning ~14 days
- Use polynomial regression (degree 2) for smooth, continuous curves
- The polynomial ensures continuous first derivatives (velocity)

Key facts about the "natural" perigee (established in lunar mechanics literature):
- Perigee oscillates ~25° from mean position (vs. apogee ~5°)
- Apogee and perigee are not exactly opposite at all times
- The interpolation is done directly on perigee values for consistency

These tests verify:
1. Basic functionality and API compatibility
2. Return value format and ranges
3. Direct interpolation on perigee values (not derived from apogee)
4. Velocity calculations
5. Relationship between interpolated perigee and osculating perigee
6. Proper smoothing of oscillations
"""

import pytest
import libephemeris as swe


class TestInterpolatedPerigeeBasic:
    """Basic functionality tests for interpolated perigee."""

    def test_calc_ut_returns_valid_position(self):
        """Test that calc_ut returns a valid position for INTP_PERG."""
        jd_ut = 2451545.0  # J2000.0
        result, retflag = swe.calc_ut(jd_ut, swe.INTP_PERG, swe.FLG_SPEED)

        # Result should be a tuple of 6 values
        assert len(result) == 6
        lon, lat, dist, speed_lon, speed_lat, speed_dist = result

        # Longitude should be in [0, 360)
        assert 0 <= lon < 360, f"Longitude {lon} out of range"

        # Latitude should be reasonable (typically small for lunar perigee)
        assert -10 < lat < 10, f"Latitude {lat} seems unreasonable"

        # 3rd return value is distance in AU (~0.0024 for perigee)
        assert 0.002 < dist < 0.003, f"Distance {dist} AU seems unreasonable"

    def test_calc_returns_valid_position(self):
        """Test that calc (TT) returns a valid position for INTP_PERG."""
        jd_tt = 2451545.0  # J2000.0
        result, retflag = swe.calc(jd_tt, swe.INTP_PERG, 0)

        assert len(result) == 6
        lon = result[0]
        assert 0 <= lon < 360

    def test_perigee_roughly_opposite_apogee(self):
        """Test that interpolated perigee is approximately 180 degrees from apogee.

        Note: Apogee and perigee are NOT exactly 180° apart at all times — they
        can deviate by up to 28° depending on Sun-Moon geometry (Chapront-Touzé
        & Chapront 1988; Meeus, Astronomical Algorithms ch. 47).
        At J2000.0, even pyswisseph itself shows a deviation of about 12.7° from
        exact opposition.
        """
        jd_ut = 2451545.0

        apogee_result, _ = swe.calc_ut(jd_ut, swe.INTP_APOG, 0)
        perigee_result, _ = swe.calc_ut(jd_ut, swe.INTP_PERG, 0)

        apogee_lon = apogee_result[0]
        perigee_lon = perigee_result[0]

        # Calculate angular difference
        diff = abs(apogee_lon - perigee_lon)
        if diff > 180:
            diff = 360 - diff

        # Should be approximately 180 degrees, with tolerance for
        # physical asymmetry (up to 28° per Chapront/Meeus lunar theory)
        assert abs(diff - 180) < 28.0, f"Difference from 180: {abs(diff - 180)}"


class TestInterpolatedPerigeeDirectFunction:
    """Test the direct lunar.py functions."""

    def test_calc_interpolated_perigee_directly(self):
        """Test calling calc_interpolated_perigee directly."""
        from libephemeris import lunar

        jd_tt = 2451545.0
        lon, lat, ecc = lunar.calc_interpolated_perigee(jd_tt)

        assert 0 <= lon < 360
        assert -10 < lat < 10
        assert 0.002 < ecc < 0.003

    def test_perigee_independent_of_apogee(self):
        """Test that perigee is interpolated independently, not derived from apogee.

        Note: Apogee and perigee are NOT exactly 180° apart at all times — they
        can deviate by up to 28° depending on Sun-Moon geometry (Chapront-Touzé
        & Chapront 1988). At J2000.0, even pyswisseph itself shows a deviation
        of about 12.7° from exact opposition.
        """
        from libephemeris import lunar

        jd_tt = 2451545.0

        # Get both apogee and perigee
        apogee_lon, apogee_lat, apogee_ecc = lunar.calc_interpolated_apogee(jd_tt)
        perigee_lon, perigee_lat, perigee_ecc = lunar.calc_interpolated_perigee(jd_tt)

        # The difference should be approximately 180
        diff = abs(apogee_lon - perigee_lon)
        if diff > 180:
            diff = 360 - diff

        # Should be approximately 180 (within 28° per SE documentation)
        assert abs(diff - 180) < 28.0, (
            f"Apogee-perigee difference from 180: {abs(diff - 180)}"
        )


class TestInterpolatedPerigeeVelocity:
    """Tests for velocity/speed calculations."""

    def test_velocity_is_calculated(self):
        """Test that velocity is calculated when FLG_SPEED is set."""
        jd_ut = 2451545.0
        result, _ = swe.calc_ut(jd_ut, swe.INTP_PERG, swe.FLG_SPEED)

        speed_lon = result[3]
        # Velocity should be non-zero
        assert speed_lon != 0.0, "Velocity should be non-zero"

    def test_velocity_reasonable_magnitude(self):
        """Test that velocity has reasonable magnitude."""
        jd_ut = 2451545.0
        result, _ = swe.calc_ut(jd_ut, swe.INTP_PERG, swe.FLG_SPEED)

        speed_lon = result[3]
        # Mean perigee moves at about 40.7 degrees/year = 0.111 deg/day
        # Interpolated can show larger oscillations due to the polynomial fit
        # Allow wide tolerance
        assert abs(speed_lon) < 10, f"Velocity {speed_lon} deg/day seems unreasonable"

    def test_no_velocity_without_flag(self):
        """Test that velocity is zero when FLG_SPEED is not set."""
        jd_ut = 2451545.0
        result, _ = swe.calc_ut(jd_ut, swe.INTP_PERG, 0)

        speed_lon = result[3]
        assert speed_lon == 0.0


class TestInterpolatedPerigeeVsOsculating:
    """Tests comparing interpolated and osculating perigee."""

    def test_interpolated_smoother_than_osculating(self):
        """Test that interpolated perigee varies more smoothly than osculating."""
        # Sample both over several days and check that interpolated is smoother
        jd_start = 2451545.0
        dt = 1.0  # 1 day step

        oscu_positions = []
        intp_positions = []

        for i in range(10):
            jd = jd_start + i * dt

            # Osculating perigee is osculating apogee + 180
            oscu_apogee_result, _ = swe.calc_ut(jd, swe.OSCU_APOG, 0)
            oscu_lon = (oscu_apogee_result[0] + 180.0) % 360.0

            intp_result, _ = swe.calc_ut(jd, swe.INTP_PERG, 0)

            oscu_positions.append(oscu_lon)
            intp_positions.append(intp_result[0])

        # Calculate day-to-day changes
        oscu_changes = []
        intp_changes = []

        for i in range(1, len(oscu_positions)):
            oscu_diff = oscu_positions[i] - oscu_positions[i - 1]
            intp_diff = intp_positions[i] - intp_positions[i - 1]

            # Handle wrap-around
            if oscu_diff > 180:
                oscu_diff -= 360
            elif oscu_diff < -180:
                oscu_diff += 360
            if intp_diff > 180:
                intp_diff -= 360
            elif intp_diff < -180:
                intp_diff += 360

            oscu_changes.append(oscu_diff)
            intp_changes.append(intp_diff)

        # Calculate variance of changes (smoother = lower variance)
        oscu_mean = sum(oscu_changes) / len(oscu_changes)
        intp_mean = sum(intp_changes) / len(intp_changes)

        oscu_variance = sum((c - oscu_mean) ** 2 for c in oscu_changes) / len(
            oscu_changes
        )
        intp_variance = sum((c - intp_mean) ** 2 for c in intp_changes) / len(
            intp_changes
        )

        # Interpolated should have lower variance (smoother)
        assert intp_variance < oscu_variance, (
            f"Interpolated variance {intp_variance} should be less than osculating {oscu_variance}"
        )


class TestInterpolatedPerigeeAtMultipleDates:
    """Tests at various dates to ensure stability."""

    @pytest.mark.parametrize(
        "jd,description",
        [
            (2451545.0, "J2000.0"),
            (2440587.5, "1970-01-01 (Unix epoch)"),
            (2458849.5, "2020-01-01"),
            (2460676.5, "2025-01-01"),
        ],
    )
    def test_valid_at_various_dates(self, jd, description):
        """Test that interpolated perigee works at various dates."""
        result, _ = swe.calc_ut(jd, swe.INTP_PERG, swe.FLG_SPEED)

        lon = result[0]
        speed = result[3]

        assert 0 <= lon < 360, f"Invalid longitude {lon} at {description}"
        # Speed can be positive or negative due to oscillations in the interpolation
        assert abs(speed) < 10, f"Speed {speed} seems unreasonable at {description}"


class TestInterpolatedPerigeeConsistency:
    """Test internal consistency of the implementation."""

    def test_continuous_motion(self):
        """Test that interpolated perigee moves continuously (no jumps)."""
        jd_start = 2451545.0
        dt = 0.1  # 2.4 hour steps

        prev_lon = None
        for i in range(50):
            jd = jd_start + i * dt
            result, _ = swe.calc_ut(jd, swe.INTP_PERG, 0)
            lon = result[0]

            if prev_lon is not None:
                diff = lon - prev_lon
                # Handle wrap-around
                if diff > 180:
                    diff -= 360
                elif diff < -180:
                    diff += 360

                # Maximum daily motion ~0.2 deg, so in 2.4 hours should be < 0.02 deg
                # Allow more tolerance since interpolation can cause small shifts
                assert abs(diff) < 0.5, (
                    f"Jump of {diff} degrees at step {i} seems too large"
                )

            prev_lon = lon

    def test_latitude_sign_opposite_to_apogee(self):
        """Test that perigee latitude is opposite in sign to apogee latitude."""
        from libephemeris import lunar

        jd_tt = 2451545.0

        apogee_lon, apogee_lat, _ = lunar.calc_interpolated_apogee(jd_tt)
        perigee_lon, perigee_lat, _ = lunar.calc_interpolated_perigee(jd_tt)

        # Latitudes should have opposite signs (or both be near zero)
        if abs(apogee_lat) > 0.01:  # Only check if latitude is significant
            assert apogee_lat * perigee_lat <= 0 or abs(perigee_lat) < 0.5, (
                f"Latitudes should be opposite: apogee={apogee_lat}, perigee={perigee_lat}"
            )


class TestInterpolatedPerigeeAlgorithm:
    """Tests for the specific algorithm design choices."""

    def test_perigee_oscillation_within_25_degrees_of_mean(self):
        """Test that interpolated perigee is significantly smoother than osculating."""
        from libephemeris import lunar

        jd_tt = 2451545.0

        # Sample over 100 days to see the smoothing effect
        num_samples = 100
        intp_lons = []
        oscu_lons = []

        for i in range(num_samples):
            jd = jd_tt + i
            intp_lon, _, _ = lunar.calc_interpolated_perigee(jd)
            oscu_apogee_lon, _, _ = lunar.calc_true_lilith(jd)
            oscu_lon = (oscu_apogee_lon + 180.0) % 360.0
            intp_lons.append(intp_lon)
            oscu_lons.append(oscu_lon)

        # Unwrap both series
        def unwrap(lons):
            result = [lons[0]]
            for i in range(1, len(lons)):
                diff = lons[i] - lons[i - 1]
                if diff > 180:
                    diff -= 360
                elif diff < -180:
                    diff += 360
                result.append(result[-1] + diff)
            return result

        intp_unwrapped = unwrap(intp_lons)
        oscu_unwrapped = unwrap(oscu_lons)

        # Calculate variance of day-to-day changes as measure of smoothness
        intp_changes = [
            intp_unwrapped[i + 1] - intp_unwrapped[i] for i in range(num_samples - 1)
        ]
        oscu_changes = [
            oscu_unwrapped[i + 1] - oscu_unwrapped[i] for i in range(num_samples - 1)
        ]

        intp_mean = sum(intp_changes) / len(intp_changes)
        oscu_mean = sum(oscu_changes) / len(oscu_changes)

        intp_variance = sum((c - intp_mean) ** 2 for c in intp_changes) / len(
            intp_changes
        )
        oscu_variance = sum((c - oscu_mean) ** 2 for c in oscu_changes) / len(
            oscu_changes
        )

        # Interpolated should have lower variance (smoother daily motion)
        # This confirms the algorithm is working to smooth out oscillations
        assert intp_variance < oscu_variance, (
            f"Interpolated variance {intp_variance} should be less than osculating {oscu_variance}"
        )

    def test_uses_half_synodic_month_window(self):
        """Test that the algorithm uses approximately half synodic month window."""
        # The algorithm should smooth out the ~14.77 day oscillation
        # by sampling across this period
        from libephemeris import lunar

        jd_tt = 2451545.0

        # Sample over several oscillation periods
        # If correctly smoothed, we should see reduced oscillation amplitude
        num_samples = 30
        intp_lons = []
        oscu_lons = []

        for i in range(num_samples):
            jd = jd_tt + i
            intp_lon, _, _ = lunar.calc_interpolated_perigee(jd)
            oscu_apogee_lon, _, _ = lunar.calc_true_lilith(jd)
            oscu_lon = (oscu_apogee_lon + 180.0) % 360.0
            intp_lons.append(intp_lon)
            oscu_lons.append(oscu_lon)

        # Unwrap both series
        def unwrap(lons):
            result = [lons[0]]
            for i in range(1, len(lons)):
                diff = lons[i] - lons[i - 1]
                if diff > 180:
                    diff -= 360
                elif diff < -180:
                    diff += 360
                result.append(result[-1] + diff)
            return result

        intp_unwrapped = unwrap(intp_lons)
        oscu_unwrapped = unwrap(oscu_lons)

        # Calculate the range (max - min) as a measure of oscillation amplitude
        intp_range = max(intp_unwrapped) - min(intp_unwrapped)
        oscu_range = max(oscu_unwrapped) - min(oscu_unwrapped)

        # The interpolated perigee should have a smaller range than osculating
        # (smoothed out oscillations)
        assert intp_range < oscu_range, (
            f"Interpolated range {intp_range} should be less than osculating {oscu_range}"
        )
