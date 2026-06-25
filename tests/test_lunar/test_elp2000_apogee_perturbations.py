"""
Tests for the ELP2000-82B apogee perturbation series implementation.

This module tests the _calc_elp2000_apogee_perturbations() function which
implements the analytical perturbation series for the lunar apsidal line.

The perturbation series models the ~5° oscillations of the apogee around its
mean position, following the ELP2000-82B lunar theory by Chapront-Touze & Chapront.

Key tests:
1. Perturbation returns reasonable values (within expected range)
2. Perturbation varies smoothly over time (no discontinuities)
3. Dominant terms have expected periods and amplitudes
4. Integration with calc_interpolated_apogee() produces improved results
"""

import pytest
from libephemeris.lunar import (
    _calc_elp2000_apogee_perturbations,
    calc_interpolated_apogee,
    calc_mean_lilith,
)


class TestELP2000ApogeePerturbationsBasic:
    """Basic functionality tests for apogee perturbation series."""

    def test_perturbation_returns_float(self):
        """Test that the perturbation function returns a float value."""
        jd_tt = 2451545.0  # J2000.0
        result = _calc_elp2000_apogee_perturbations(jd_tt)
        assert isinstance(result, float)

    def test_perturbation_within_expected_range(self):
        """Test that perturbation values are within expected range (~5° oscillation)."""
        # Sample across a lunar month (where perturbations vary significantly)
        jd_start = 2451545.0
        for i in range(30):
            jd = jd_start + i
            perturbation = _calc_elp2000_apogee_perturbations(jd)
            # Perturbation should be within approximately ±5° for apogee
            # (larger than this would indicate incorrect coefficients)
            assert -10.0 < perturbation < 10.0, (
                f"Perturbation {perturbation}° at JD {jd} outside expected range"
            )

    def test_perturbation_varies_over_time(self):
        """Test that perturbation values change over time (not constant)."""
        jd1 = 2451545.0
        jd2 = 2451552.0  # 7 days later (half synodic oscillation)

        pert1 = _calc_elp2000_apogee_perturbations(jd1)
        pert2 = _calc_elp2000_apogee_perturbations(jd2)

        # Values should differ by at least 0.1° over 7 days
        assert abs(pert2 - pert1) > 0.1, (
            f"Perturbation did not change significantly: {pert1}° to {pert2}°"
        )

    def test_perturbation_continuous(self):
        """Test that perturbation changes smoothly (no discontinuities)."""
        jd_start = 2451545.0
        previous = _calc_elp2000_apogee_perturbations(jd_start)

        # Check 100 consecutive days
        for i in range(1, 100):
            jd = jd_start + i
            current = _calc_elp2000_apogee_perturbations(jd)

            # Day-to-day change should not exceed 1.5° (perturbations are smooth)
            # The 2D term has period ~14.77 days and amplitude ~2.1°, so max
            # daily change is about 2π × 2.1 / 14.77 ≈ 0.9°, but multiple terms
            # can combine to produce slightly larger changes.
            change = abs(current - previous)
            assert change < 1.5, (
                f"Perturbation jumped by {change}° from day {i - 1} to {i}"
            )
            previous = current


class TestELP2000ApogeePerturbationDominantTerms:
    """Tests for dominant perturbation term characteristics."""

    def test_fortnightly_oscillation_present(self):
        """Test that the dominant 2D (~14.77 day) oscillation is present."""
        jd_start = 2451545.0
        synodic_half_period = 14.77  # Half synodic month

        # Sample at quarter-period intervals to detect oscillation
        values = []
        for i in range(8):  # 8 samples over ~2 periods
            jd = jd_start + i * synodic_half_period / 2
            values.append(_calc_elp2000_apogee_perturbations(jd))

        # Calculate max-min range (should reflect 2D amplitude of ~2.1°)
        value_range = max(values) - min(values)
        assert value_range > 1.0, (
            f"Fortnightly oscillation range {value_range}° smaller than expected"
        )

    def test_evection_related_terms_present(self):
        """Test that evection-related perturbations affect the result."""
        # Evection has ~31.8 day period, test over two periods
        jd_start = 2451545.0
        evection_period = 31.8

        values = []
        for i in range(12):  # 12 samples over ~2 evection periods
            jd = jd_start + i * evection_period / 6
            values.append(_calc_elp2000_apogee_perturbations(jd))

        # There should be variation due to evection terms
        value_range = max(values) - min(values)
        assert value_range > 0.5, (
            f"Evection-period variation {value_range}° smaller than expected"
        )


class TestCalcInterpolatedApogeeWithELP2000:
    """Tests for the integrated calc_interpolated_apogee using ELP2000-82B."""

    def test_interpolated_apogee_returns_valid_tuple(self):
        """Test that calc_interpolated_apogee returns valid values."""
        jd_tt = 2451545.0
        lon, lat, ecc = calc_interpolated_apogee(jd_tt)

        # Longitude in [0, 360)
        assert 0.0 <= lon < 360.0, f"Longitude {lon}° outside valid range"

        # Latitude should be small (near ecliptic)
        assert -5.0 <= lat <= 5.0, f"Latitude {lat}° unexpectedly large"

        # 3rd return value is distance in AU (not eccentricity)
        assert 0.002 < ecc < 0.003, f"Distance {ecc} AU outside expected range"

    def test_interpolated_apogee_differs_from_mean(self):
        """Test that interpolated apogee differs from mean by perturbation + correction."""
        jd_tt = 2451545.0

        mean_lon = calc_mean_lilith(jd_tt)
        interp_lon, _, _ = calc_interpolated_apogee(jd_tt)
        perturbation = _calc_elp2000_apogee_perturbations(jd_tt)

        # The difference should approximately equal the perturbation
        # plus the apse correction table residual (up to ~0.15° at various dates)
        diff = interp_lon - mean_lon
        # Handle wrap-around
        if diff > 180:
            diff -= 360
        elif diff < -180:
            diff += 360

        assert abs(diff - perturbation) < 0.2, (
            f"Interpolated-mean difference {diff}° doesn't match perturbation {perturbation}°"
        )

    def test_interpolated_apogee_smooth_motion(self):
        """Test that interpolated apogee moves smoothly over time."""
        jd_start = 2451545.0
        positions = []

        for i in range(30):
            jd = jd_start + i
            lon, _, _ = calc_interpolated_apogee(jd)
            positions.append(lon)

        # Calculate daily changes
        changes = []
        for i in range(1, len(positions)):
            diff = positions[i] - positions[i - 1]
            # Handle wrap-around
            if diff > 180:
                diff -= 360
            elif diff < -180:
                diff += 360
            changes.append(diff)

        # Mean daily motion should be ~0.11° (40°/year / 365 days)
        mean_change = sum(changes) / len(changes)
        assert 0.0 < mean_change < 0.3, f"Mean daily motion {mean_change}° unexpected"

        # Variance should be low (smooth motion)
        variance = sum((c - mean_change) ** 2 for c in changes) / len(changes)
        assert variance < 1.0, f"Motion variance {variance} too high (not smooth)"

    def test_interpolated_apogee_multi_date_consistency(self):
        """Test interpolated apogee across multiple dates."""
        test_dates = [
            2451545.0,  # J2000.0
            2458849.5,  # 2020-01-01
            2459580.5,  # 2022-01-01
            2460310.5,  # 2024-01-01
        ]

        for jd in test_dates:
            lon, lat, ecc = calc_interpolated_apogee(jd)

            # All values should be valid
            assert 0.0 <= lon < 360.0
            assert -6.0 <= lat <= 6.0
            assert 0.002 < ecc < 0.003


class TestELP2000ApogeePerturbationEdgeCases:
    """Edge case tests for apogee perturbation calculations."""

    def test_perturbation_at_j2000(self):
        """Test perturbation value at J2000.0 epoch."""
        jd_j2000 = 2451545.0
        perturbation = _calc_elp2000_apogee_perturbations(jd_j2000)

        # Should return a valid value within expected range
        assert -10.0 < perturbation < 10.0
        assert isinstance(perturbation, float)

    def test_perturbation_distant_past(self):
        """Test perturbation calculation for distant past date."""
        jd_1900 = 2415020.5  # 1900-01-01
        perturbation = _calc_elp2000_apogee_perturbations(jd_1900)

        # Should still return a valid value
        assert -15.0 < perturbation < 15.0

    def test_perturbation_distant_future(self):
        """Test perturbation calculation for distant future date."""
        jd_2100 = 2488069.5  # 2100-01-01
        perturbation = _calc_elp2000_apogee_perturbations(jd_2100)

        # Should still return a valid value
        assert -15.0 < perturbation < 15.0

    def test_interpolated_apogee_no_exceptions(self):
        """Test that calc_interpolated_apogee doesn't raise exceptions."""
        test_dates = [
            2415020.5,  # 1900
            2440587.5,  # 1970
            2451545.0,  # 2000
            2459580.5,  # 2022
            2488069.5,  # 2100
        ]

        for jd in test_dates:
            # Should not raise any exceptions
            try:
                lon, lat, ecc = calc_interpolated_apogee(jd)
                # Basic validation
                assert isinstance(lon, float)
                assert isinstance(lat, float)
                assert isinstance(ecc, float)
            except Exception as e:
                pytest.fail(
                    f"calc_interpolated_apogee raised {type(e).__name__} at JD {jd}: {e}"
                )
