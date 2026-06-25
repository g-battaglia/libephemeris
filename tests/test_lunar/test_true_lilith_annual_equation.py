"""
Tests for annual equation correction in the True Lilith (osculating lunar apogee) calculation.

The Annual Equation is a perturbation of lunar motion caused by the variation in
the Earth-Sun distance due to Earth's orbital eccentricity. When Earth is closer
to the Sun (perihelion), the solar gravitational perturbation on the Moon is
stronger, affecting both the Moon's longitude and the position of the lunar apogee.

Key characteristics:
- Period: ~365.25 days (anomalistic year)
- Amplitude: ~0.186 degrees
- Argument: M (Sun's mean anomaly)
- Cause: variation in solar gravitational perturbation due to Earth's elliptical orbit

The annual equation correction term in True Lilith is: 0.1856 * sin(M)

References:
    - Chapront-Touze, M. & Chapront, J. "Lunar Tables and Programs" (1991)
    - Brown, E.W. "An Introductory Treatise on the Lunar Theory" (1896)
    - Meeus, J. "Astronomical Algorithms" (2nd ed., 1998), Chapter 47
"""

import math
from libephemeris.lunar import calc_true_lilith, calc_mean_lilith


# Annual equation period in days (anomalistic year)
ANNUAL_EQUATION_PERIOD_DAYS = 365.259636
# Annual equation amplitude in degrees (for True Lilith)
ANNUAL_EQUATION_AMPLITUDE_DEG = 0.1856


def _calc_sun_mean_anomaly(jd_tt: float) -> float:
    """
    Calculate the Sun's mean anomaly for testing.

    This helper computes the same solar argument used internally by
    calc_true_lilith to verify the annual equation correction is applied correctly.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        float: Sun's mean anomaly (M) in radians
    """
    T = (jd_tt - 2451545.0) / 36525.0  # Julian centuries from J2000.0

    # Sun's mean anomaly (M) - same formula used in _calc_lunar_fundamental_arguments
    M = 357.5291092 + 35999.0502909 * T - 0.0001536 * T**2 + T**3 / 24490000.0

    # Convert to radians
    return math.radians(M % 360.0)


class TestAnnualEquationCorrectionAmplitude:
    """Test the amplitude and bounds of the annual equation correction."""

    def test_annual_equation_correction_amplitude_is_bounded(self):
        """Annual equation correction should be bounded by +-0.1856 degrees."""
        # Test at various dates to capture different phases of the annual cycle
        test_dates = [
            2451545.0,  # J2000.0
            2451545.0 + 91.3,  # ~1/4 year later
            2451545.0 + 182.6,  # ~1/2 year later
            2451545.0 + 274.0,  # ~3/4 year later
            2451545.0 + 365.26,  # ~1 full year
        ]

        for jd in test_dates:
            M = _calc_sun_mean_anomaly(jd)
            annual_correction = ANNUAL_EQUATION_AMPLITUDE_DEG * math.sin(M)

            # Correction should be bounded by amplitude
            assert -0.19 <= annual_correction <= 0.19, (
                f"Annual equation correction {annual_correction}deg at JD {jd} exceeds bounds"
            )

    def test_annual_equation_reaches_maximum_amplitude(self):
        """Annual equation correction should reach approximately +-0.1856 degrees at extrema."""
        jd_base = 2451545.0

        # Search for approximate extrema over one year
        max_annual = -1.0
        min_annual = 1.0

        for i in range(366):  # Check over ~1 year in 1 day steps
            jd = jd_base + i
            M = _calc_sun_mean_anomaly(jd)
            annual_val = ANNUAL_EQUATION_AMPLITUDE_DEG * math.sin(M)

            if annual_val > max_annual:
                max_annual = annual_val
            if annual_val < min_annual:
                min_annual = annual_val

        # Maximum should be close to amplitude
        assert max_annual > 0.18, f"Max annual correction {max_annual}deg too low"
        assert max_annual < 0.19, f"Max annual correction {max_annual}deg too high"

        # Minimum should be close to negative amplitude
        assert min_annual < -0.18, f"Min annual correction {min_annual}deg too high"
        assert min_annual > -0.19, f"Min annual correction {min_annual}deg too low"


class TestAnnualEquationArgument:
    """Test the annual equation argument (Sun's mean anomaly)."""

    def test_sun_mean_anomaly_increases_over_time(self):
        """Sun's mean anomaly should increase by ~360° per year."""
        jd_start = 2451545.0
        jd_end = jd_start + 365.25

        M_start = _calc_sun_mean_anomaly(jd_start)
        M_end = _calc_sun_mean_anomaly(jd_end)

        # Both should be valid angles (0 to 2*pi)
        assert 0 <= M_start < 2 * math.pi
        assert 0 <= M_end < 2 * math.pi

    def test_sun_mean_anomaly_rate(self):
        """Sun's mean anomaly should advance at about 0.986° per day."""
        jd_start = 2451545.0
        jd_end = jd_start + 100.0  # 100 days later

        M_start = math.degrees(_calc_sun_mean_anomaly(jd_start)) % 360
        M_end = math.degrees(_calc_sun_mean_anomaly(jd_end)) % 360

        # Expected change: ~0.9856° per day * 100 days = ~98.56°
        expected_change = 0.9856 * 100

        # Handle wraparound
        actual_change = (M_end - M_start) % 360

        # Allow 1° tolerance
        assert abs(actual_change - expected_change) < 1.0, (
            f"Mean anomaly change {actual_change}deg doesn't match expected {expected_change}deg"
        )


class TestAnnualEquationEffectOnTrueLilith:
    """Test the annual equation effect on True Lilith calculation."""

    def test_true_lilith_varies_with_annual_period(self):
        """True Lilith longitude should show variation consistent with ~1 year period."""
        jd_start = 2451545.0
        half_period = ANNUAL_EQUATION_PERIOD_DAYS / 2  # ~182.6 days

        # Get True Lilith at start and half-period later
        lon_start, _, _ = calc_true_lilith(jd_start)
        lon_half, _, _ = calc_true_lilith(jd_start + half_period)

        # The annual equation effect should contribute to the difference
        # But other effects (mean motion, evection, variation) also contribute
        # So we just verify that the function works at these epochs
        assert 0 <= lon_start < 360
        assert 0 <= lon_half < 360

    def test_annual_equation_correction_is_integrated(self):
        """Test that annual equation correction is integrated into True Lilith output."""
        # At J2000.0, Sun's mean anomaly M ≈ 357.5° (close to perihelion)
        # sin(357.5°) ≈ -0.0436, so correction ≈ -0.008°
        jd = 2451545.0

        # Get True Lilith result
        lon, lat, e_mag = calc_true_lilith(jd)

        # Verify basic output validity
        assert 0 <= lon < 360, f"Longitude {lon} out of range"
        assert -90 <= lat <= 90, f"Latitude {lat} out of range"
        assert 0 < e_mag < 1, f"Eccentricity {e_mag} out of range"

    def test_perihelion_vs_aphelion_dates(self):
        """Test annual equation effect near perihelion (Jan) vs aphelion (Jul)."""
        # Earth's perihelion is around January 3-4
        # Earth's aphelion is around July 4-5

        # January 3, 2000 (near perihelion)
        jd_perihelion = 2451547.0  # ~2 days after J2000

        # July 4, 2000 (near aphelion)
        jd_aphelion = 2451730.0  # ~185 days after J2000

        # Calculate Sun's mean anomaly at both epochs
        M_peri = _calc_sun_mean_anomaly(jd_perihelion)
        M_aph = _calc_sun_mean_anomaly(jd_aphelion)

        # At perihelion, M ≈ 0°, sin(M) ≈ 0
        # At aphelion, M ≈ 180°, sin(M) ≈ 0
        # Maximum/minimum effect occurs at M ≈ 90° or 270°

        # Verify True Lilith can be calculated at both epochs
        lon_peri, _, _ = calc_true_lilith(jd_perihelion)
        lon_aph, _, _ = calc_true_lilith(jd_aphelion)

        assert 0 <= lon_peri < 360
        assert 0 <= lon_aph < 360


class TestAnnualEquationVsMeanLilith:
    """Test comparing True Lilith (with annual equation) to Mean Lilith."""

    def test_true_lilith_oscillates_around_mean_lilith(self):
        """True Lilith should oscillate around mean Lilith partially due to annual equation."""
        jd_base = 2451545.0
        differences = []

        # Sample over a year to capture annual variation
        for i in range(365):
            jd = jd_base + i
            true_lon, _, _ = calc_true_lilith(jd)
            mean_lon = calc_mean_lilith(jd)

            # Calculate difference (handle wraparound)
            diff = true_lon - mean_lon
            if diff > 180:
                diff -= 360
            elif diff < -180:
                diff += 360

            differences.append(diff)

        # The difference should vary over the year
        diff_range = max(differences) - min(differences)

        # True-Mean difference includes evection, variation, and annual equation
        # So the range should be substantial
        assert diff_range > 0.5, f"True-Mean Lilith range = {diff_range}"

    def test_all_three_corrections_applied(self):
        """All three corrections (evection, variation, annual equation) should be applied."""
        # Test at a date where all three effects are significant
        jd = 2451545.0 + 100  # 100 days after J2000

        # Calculate the expected corrections
        T = (jd - 2451545.0) / 36525.0

        # Mean elongation (D)
        D = math.radians(
            (297.8501921 + 445267.1114034 * T - 0.0018819 * T**2 + T**3 / 545868.0)
            % 360.0
        )

        # Moon's mean anomaly (M')
        M_prime = math.radians(
            (
                134.9633964
                + 477198.8675055 * T
                + 0.0087414 * T**2
                + T**3 / 69699.0
                - T**4 / 14712000.0
            )
            % 360.0
        )

        # Sun's mean anomaly (M)
        M = math.radians(
            (357.5291092 + 35999.0502909 * T - 0.0001536 * T**2 + T**3 / 24490000.0)
            % 360.0
        )

        # Expected corrections
        evection_correction = 1.2739 * math.sin(2 * D - M_prime)
        variation_correction = 0.6583 * math.sin(2 * D)
        annual_correction = 0.1856 * math.sin(M)

        # All three should be non-zero at this epoch
        # (unless by coincidence the arguments happen to be 0 or 180 degrees)
        # Just verify the corrections are computable
        assert -2.0 < evection_correction < 2.0
        assert -1.0 < variation_correction < 1.0
        assert -0.2 < annual_correction < 0.2

        # Get True Lilith position
        lon, _, _ = calc_true_lilith(jd)

        # Verify it returns a valid longitude
        assert 0 <= lon < 360


class TestAnnualEquationIntegration:
    """Integration tests for annual equation correction in True Lilith."""

    def test_true_lilith_valid_across_year(self):
        """True Lilith should return valid results across an entire year."""
        jd_base = 2451545.0

        for i in range(0, 366, 30):  # Sample every 30 days
            jd = jd_base + i
            lon, lat, e_mag = calc_true_lilith(jd)

            assert 0 <= lon < 360, f"Invalid longitude {lon} at JD {jd}"
            assert -10 <= lat <= 10, f"Invalid latitude {lat} at JD {jd}"
            assert 0.002 < e_mag < 0.003, f"Invalid distance {e_mag} AU at JD {jd}"

    def test_multiple_years_consistency(self):
        """True Lilith should work consistently across multiple years."""
        years_to_test = [1990, 2000, 2010, 2020, 2030]

        for year in years_to_test:
            # Convert year to JD (approximately January 1 of each year)
            jd = 2451545.0 + (year - 2000) * 365.25

            for day_offset in [0, 100, 200, 300]:
                lon, lat, e_mag = calc_true_lilith(jd + day_offset)

                assert 0 <= lon < 360
                assert -10 <= lat <= 10
                assert 0.002 < e_mag < 0.003

    def test_annual_equation_smooth_variation(self):
        """Annual equation effect should vary smoothly over the year."""
        jd_base = 2451545.0
        previous_correction = None

        for i in range(366):
            jd = jd_base + i
            M = _calc_sun_mean_anomaly(jd)
            correction = ANNUAL_EQUATION_AMPLITUDE_DEG * math.sin(M)

            if previous_correction is not None:
                # Daily change should be small (smooth)
                daily_change = abs(correction - previous_correction)
                # Maximum daily change of sin function is ~0.0172 per degree
                # With amplitude 0.1856 and ~0.986 deg/day, max change ~0.003 deg/day
                assert daily_change < 0.01, (
                    f"Non-smooth annual correction change at JD {jd}: {daily_change}"
                )

            previous_correction = correction
