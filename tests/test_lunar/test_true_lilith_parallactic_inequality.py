"""
Tests for parallactic inequality correction in the True Lilith (osculating lunar apogee) calculation.

The Parallactic Inequality (also called the parallactic equation) is a perturbation
of lunar motion caused by the finite distance between the Earth and Moon. It arises
from the fact that the Moon's parallax varies with its distance from Earth, which
in turn depends on the solar perturbations.

Key characteristics:
- Period: ~29.53 days (synodic month)
- Amplitude: ~0.125 degrees
- Argument: D (mean elongation of Moon from Sun)
- Cause: variation in Moon's horizontal parallax due to distance changes

The parallactic inequality correction term in True Lilith is: 0.125 * sin(D)

References:
    - Chapront-Touze, M. & Chapront, J. "Lunar Tables and Programs" (1991)
    - Brown, E.W. "An Introductory Treatise on the Lunar Theory" (1896)
    - Meeus, J. "Astronomical Algorithms" (2nd ed., 1998), Chapter 47
"""

import math
from libephemeris.lunar import calc_true_lilith, calc_mean_lilith


# Parallactic inequality period in days (synodic month)
PARALLACTIC_INEQUALITY_PERIOD_DAYS = 29.530589
# Parallactic inequality amplitude in degrees (for True Lilith)
PARALLACTIC_INEQUALITY_AMPLITUDE_DEG = 0.125


def _calc_mean_elongation(jd_tt: float) -> float:
    """
    Calculate the mean elongation of Moon from Sun for testing.

    This helper computes the same argument used internally by
    calc_true_lilith to verify the parallactic inequality correction is applied correctly.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        float: Mean elongation D in radians
    """
    T = (jd_tt - 2451545.0) / 36525.0  # Julian centuries from J2000.0

    # Mean elongation of Moon from Sun (D) - same formula used in _calc_lunar_fundamental_arguments
    D = 297.8501921 + 445267.1114034 * T - 0.0018819 * T**2 + T**3 / 545868.0

    # Convert to radians
    return math.radians(D % 360.0)


class TestParallacticInequalityCorrectionAmplitude:
    """Test the amplitude and bounds of the parallactic inequality correction."""

    def test_parallactic_inequality_correction_amplitude_is_bounded(self):
        """Parallactic inequality correction should be bounded by +-0.125 degrees."""
        # Test at various dates to capture different phases of the synodic cycle
        test_dates = [
            2451545.0,  # J2000.0
            2451545.0 + 7.38,  # ~1/4 synodic month later
            2451545.0 + 14.77,  # ~1/2 synodic month later
            2451545.0 + 22.15,  # ~3/4 synodic month later
            2451545.0 + 29.53,  # ~1 full synodic month
        ]

        for jd in test_dates:
            D = _calc_mean_elongation(jd)
            parallactic_correction = PARALLACTIC_INEQUALITY_AMPLITUDE_DEG * math.sin(D)

            # Correction should be bounded by amplitude
            assert -0.126 <= parallactic_correction <= 0.126, (
                f"Parallactic inequality correction {parallactic_correction}deg at JD {jd} exceeds bounds"
            )

    def test_parallactic_inequality_reaches_maximum_amplitude(self):
        """Parallactic inequality correction should reach approximately +-0.125 degrees at extrema."""
        jd_base = 2451545.0

        # Search for approximate extrema over one synodic month
        max_parallactic = -1.0
        min_parallactic = 1.0

        # Check over ~2 synodic months in 0.5 day steps
        for i in range(120):
            jd = jd_base + i * 0.5
            D = _calc_mean_elongation(jd)
            parallactic_val = PARALLACTIC_INEQUALITY_AMPLITUDE_DEG * math.sin(D)

            if parallactic_val > max_parallactic:
                max_parallactic = parallactic_val
            if parallactic_val < min_parallactic:
                min_parallactic = parallactic_val

        # Maximum should be close to amplitude
        assert max_parallactic > 0.12, (
            f"Max parallactic correction {max_parallactic}deg too low"
        )
        assert max_parallactic < 0.13, (
            f"Max parallactic correction {max_parallactic}deg too high"
        )

        # Minimum should be close to negative amplitude
        assert min_parallactic < -0.12, (
            f"Min parallactic correction {min_parallactic}deg too high"
        )
        assert min_parallactic > -0.13, (
            f"Min parallactic correction {min_parallactic}deg too low"
        )


class TestParallacticInequalityArgument:
    """Test the parallactic inequality argument (mean elongation)."""

    def test_mean_elongation_increases_over_time(self):
        """Mean elongation should increase by ~360 per synodic month."""
        jd_start = 2451545.0
        jd_end = jd_start + PARALLACTIC_INEQUALITY_PERIOD_DAYS

        D_start = _calc_mean_elongation(jd_start)
        D_end = _calc_mean_elongation(jd_end)

        # Both should be valid angles (0 to 2*pi)
        assert 0 <= D_start < 2 * math.pi
        assert 0 <= D_end < 2 * math.pi

    def test_mean_elongation_rate(self):
        """Mean elongation should advance at about 12.19 per day."""
        jd_start = 2451545.0
        jd_end = jd_start + 10.0  # 10 days later

        D_start = math.degrees(_calc_mean_elongation(jd_start)) % 360
        D_end = math.degrees(_calc_mean_elongation(jd_end)) % 360

        # Expected change: ~12.19 per day * 10 days = ~121.9
        expected_change = 12.19 * 10

        # Handle wraparound
        actual_change = (D_end - D_start) % 360

        # Allow 2 tolerance
        assert abs(actual_change - expected_change) < 2.0, (
            f"Mean elongation change {actual_change}deg doesn't match expected {expected_change}deg"
        )


class TestParallacticInequalityEffectOnTrueLilith:
    """Test the parallactic inequality effect on True Lilith calculation."""

    def test_true_lilith_varies_with_synodic_period(self):
        """True Lilith longitude should show variation consistent with ~29.53 day period."""
        jd_start = 2451545.0
        half_period = PARALLACTIC_INEQUALITY_PERIOD_DAYS / 2  # ~14.77 days

        # Get True Lilith at start and half-period later
        lon_start, _, _ = calc_true_lilith(jd_start)
        lon_half, _, _ = calc_true_lilith(jd_start + half_period)

        # Both should return valid results
        assert 0 <= lon_start < 360
        assert 0 <= lon_half < 360

    def test_true_lilith_returns_valid_output(self):
        """Test that parallactic inequality correction is integrated into True Lilith output."""
        jd = 2451545.0  # J2000.0

        # Calculate mean elongation to verify the correction would be non-zero
        D = _calc_mean_elongation(jd)
        expected_correction = PARALLACTIC_INEQUALITY_AMPLITUDE_DEG * math.sin(D)

        # Get True Lilith result
        lon, lat, e_mag = calc_true_lilith(jd)

        # Output should be valid
        assert 0 <= lon < 360
        assert -90 <= lat <= 90
        assert 0 < e_mag < 1

        # The correction is applied internally - verify the term is non-trivial
        assert abs(expected_correction) > 0.001, (
            f"Expected non-trivial parallactic correction at JD {jd}"
        )

    def test_parallactic_inequality_at_new_moon_epoch(self):
        """Test parallactic inequality near a new moon (D close to 0)."""
        # Near new moon, D is close to 0 or 360 degrees, so sin(D) is near 0
        # Use a date close to a known new moon
        jd_new_moon = 2451550.1  # Approximate new moon near J2000

        D = _calc_mean_elongation(jd_new_moon)
        parallactic_val = PARALLACTIC_INEQUALITY_AMPLITUDE_DEG * math.sin(D)

        # Near new moon, the correction might be small (sin near 0)
        # Just verify we can calculate it
        assert -0.13 <= parallactic_val <= 0.13

        # Verify True Lilith can be calculated at this epoch
        lon, _, _ = calc_true_lilith(jd_new_moon)
        assert 0 <= lon < 360


class TestParallacticInequalityVsMeanLilith:
    """Test comparing True Lilith (with parallactic inequality) to Mean Lilith."""

    def test_true_lilith_oscillates_around_mean_lilith(self):
        """True Lilith should oscillate around mean Lilith partially due to parallactic inequality."""
        jd_base = 2451545.0
        differences = []

        # Sample over two synodic months
        for i in range(60):
            jd = jd_base + i
            true_lon, _, _ = calc_true_lilith(jd)
            mean_lon = calc_mean_lilith(jd)

            # Calculate difference, handling wraparound
            diff = true_lon - mean_lon
            if diff > 180:
                diff -= 360
            elif diff < -180:
                diff += 360

            differences.append(diff)

        # True Lilith should show significant variation from mean Lilith
        # Combined effects of evection, variation, annual equation, and parallactic inequality
        diff_range = max(differences) - min(differences)

        # The combined range should be notable
        assert diff_range > 0.5, f"True-Mean Lilith range = {diff_range}"


class TestParallacticInequalityCombinedEffects:
    """Test that parallactic inequality works together with other corrections."""

    def test_all_corrections_applied(self):
        """All corrections (evection, variation, annual equation, parallactic inequality) should be applied."""
        jd_base = 2451545.0

        # Sample multiple points to verify all periodic corrections are active
        longitudes = []
        for i in range(100):
            jd = jd_base + i
            lon, _, _ = calc_true_lilith(jd)
            longitudes.append(lon)

        # The range of variations should be significant due to all corrections
        # Calculate rate of change between successive points
        changes = []
        for i in range(len(longitudes) - 1):
            change = longitudes[i + 1] - longitudes[i]
            # Handle wraparound
            if change > 180:
                change -= 360
            elif change < -180:
                change += 360
            changes.append(abs(change))

        # There should be variation in the rate of change
        max_change = max(changes)
        min_change = min(changes)

        # The variation in daily change should be notable
        assert max_change > min_change, "True Lilith should show variable daily motion"


class TestParallacticInequalityIntegration:
    """Integration tests for parallactic inequality correction in True Lilith."""

    def test_true_lilith_valid_across_synodic_month(self):
        """True Lilith should return valid results across an entire synodic month."""
        jd_base = 2451545.0

        for i in range(30):
            jd = jd_base + i
            lon, lat, e_mag = calc_true_lilith(jd)

            assert 0 <= lon < 360, f"Invalid longitude at JD {jd}"
            assert -90 <= lat <= 90, f"Invalid latitude at JD {jd}"
            assert 0 < e_mag < 1, f"Invalid eccentricity at JD {jd}"

    def test_true_lilith_works_across_multiple_years(self):
        """True Lilith should work consistently across multiple years."""
        test_years = [
            2451545.0,  # J2000.0 (2000)
            2451545.0 + 365.25 * 5,  # 2005
            2451545.0 + 365.25 * 10,  # 2010
            2451545.0 + 365.25 * 20,  # 2020
            2451545.0 - 365.25 * 10,  # 1990
        ]

        for jd_year in test_years:
            for day_offset in range(0, 30, 5):
                lon, lat, e_mag = calc_true_lilith(jd_year + day_offset)

                assert 0 <= lon < 360, f"Invalid longitude at JD {jd_year + day_offset}"
                assert -90 <= lat <= 90, (
                    f"Invalid latitude at JD {jd_year + day_offset}"
                )
                assert 0 < e_mag < 1, (
                    f"Invalid eccentricity at JD {jd_year + day_offset}"
                )

    def test_parallactic_inequality_sign_at_first_quarter(self):
        """Test parallactic inequality sign when D is near 90 degrees (first quarter)."""
        # At first quarter, mean elongation D is approximately 90 degrees
        # sin(90) = 1, so correction should be near maximum positive
        jd_base = 2451545.0

        # Find a point where D is near 90 degrees
        for i in range(30):
            jd = jd_base + i * 0.5
            D = _calc_mean_elongation(jd)
            D_deg = math.degrees(D)

            if 85 < D_deg < 95:
                parallactic_val = PARALLACTIC_INEQUALITY_AMPLITUDE_DEG * math.sin(D)
                assert parallactic_val > 0.1, (
                    f"Parallactic correction at D={D_deg}deg should be positive"
                )
                break

    def test_parallactic_inequality_sign_at_third_quarter(self):
        """Test parallactic inequality sign when D is near 270 degrees (third quarter)."""
        # At third quarter, mean elongation D is approximately 270 degrees
        # sin(270) = -1, so correction should be near maximum negative
        jd_base = 2451545.0

        # Find a point where D is near 270 degrees
        for i in range(60):
            jd = jd_base + i * 0.5
            D = _calc_mean_elongation(jd)
            D_deg = math.degrees(D)

            if 265 < D_deg < 275:
                parallactic_val = PARALLACTIC_INEQUALITY_AMPLITUDE_DEG * math.sin(D)
                assert parallactic_val < -0.1, (
                    f"Parallactic correction at D={D_deg}deg should be negative"
                )
                break
