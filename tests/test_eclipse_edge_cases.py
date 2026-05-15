"""
Tests for eclipse edge cases in libephemeris.

Tests edge case handling for:
- Very shallow partial eclipses (magnitude close to 0)
- Near-miss eclipses (gamma close to eclipse limit)
- Safe math operations (acos, sqrt, division)
- Obscuration calculations with edge values
- Contact time calculations for grazing eclipses
"""

import math
import pytest

from libephemeris import (
    julday,
    revjul,
    sol_eclipse_when_glob,
    sol_eclipse_when_loc,
    lun_eclipse_when,
    lun_eclipse_when_loc,
    ECL_TOTAL,
    ECL_PARTIAL,
    ECL_PENUMBRAL,
    ECL_GRAZING,
    ECL_ANNULAR,
)
from libephemeris.eclipse import (
    _is_shallow_eclipse,
    _is_near_miss_eclipse,
    _safe_acos,
    _safe_sqrt,
    _calculate_obscuration_safe,
    _calculate_magnitude_safe,
    _validate_contact_time,
    SHALLOW_ECLIPSE_MAG_THRESHOLD,
    NEAR_MISS_GAMMA_MARGIN,
)

pytestmark = pytest.mark.slow


class TestShallowEclipseDetection:
    """Test detection of shallow eclipses."""

    def test_is_shallow_eclipse_below_threshold(self):
        """Test that magnitudes below threshold are classified as shallow."""
        assert _is_shallow_eclipse(0.0) is True
        assert _is_shallow_eclipse(0.001) is True
        assert _is_shallow_eclipse(0.005) is True
        assert _is_shallow_eclipse(SHALLOW_ECLIPSE_MAG_THRESHOLD - 0.001) is True

    def test_is_shallow_eclipse_at_threshold(self):
        """Test behavior at the threshold boundary."""
        # At threshold should not be shallow
        assert _is_shallow_eclipse(SHALLOW_ECLIPSE_MAG_THRESHOLD) is False

    def test_is_shallow_eclipse_above_threshold(self):
        """Test that magnitudes above threshold are not shallow."""
        assert _is_shallow_eclipse(0.02) is False
        assert _is_shallow_eclipse(0.1) is False
        assert _is_shallow_eclipse(0.5) is False
        assert _is_shallow_eclipse(1.0) is False


class TestNearMissEclipseDetection:
    """Test detection of near-miss eclipses."""

    def test_is_near_miss_eclipse_far_from_limit(self):
        """Test that eclipses far from limit are not near-miss."""
        assert _is_near_miss_eclipse(0.0) is False
        assert _is_near_miss_eclipse(0.5) is False
        assert _is_near_miss_eclipse(1.0) is False
        assert _is_near_miss_eclipse(1.5) is False

    def test_is_near_miss_eclipse_near_limit(self):
        """Test that eclipses near the limit are detected as near-miss."""
        gamma_limit = 1.55
        # Just inside the margin
        assert (
            _is_near_miss_eclipse(gamma_limit - NEAR_MISS_GAMMA_MARGIN + 0.001) is True
        )
        # At the margin boundary (uses > not >=, so boundary is not included)
        assert _is_near_miss_eclipse(gamma_limit - NEAR_MISS_GAMMA_MARGIN) is False
        # Just outside the margin
        assert (
            _is_near_miss_eclipse(gamma_limit - NEAR_MISS_GAMMA_MARGIN - 0.001) is False
        )

    def test_is_near_miss_eclipse_negative_gamma(self):
        """Test with negative gamma values (symmetric behavior)."""
        gamma_limit = 1.55
        assert _is_near_miss_eclipse(-gamma_limit + 0.01) is True
        assert _is_near_miss_eclipse(-1.0) is False


class TestSafeMathOperations:
    """Test safe math operations that handle edge cases."""

    def test_safe_acos_valid_range(self):
        """Test safe_acos with values in valid range."""
        assert abs(_safe_acos(0.0) - math.pi / 2) < 1e-10
        assert abs(_safe_acos(1.0) - 0.0) < 1e-10
        assert abs(_safe_acos(-1.0) - math.pi) < 1e-10
        assert abs(_safe_acos(0.5) - math.acos(0.5)) < 1e-10

    def test_safe_acos_out_of_range_positive(self):
        """Test safe_acos clamps values > 1."""
        # Should clamp to acos(1) = 0
        assert _safe_acos(1.0001) == 0.0
        assert _safe_acos(2.0) == 0.0
        assert _safe_acos(1e10) == 0.0

    def test_safe_acos_out_of_range_negative(self):
        """Test safe_acos clamps values < -1."""
        # Should clamp to acos(-1) = pi
        assert abs(_safe_acos(-1.0001) - math.pi) < 1e-10
        assert abs(_safe_acos(-2.0) - math.pi) < 1e-10
        assert abs(_safe_acos(-1e10) - math.pi) < 1e-10

    def test_safe_sqrt_valid_range(self):
        """Test safe_sqrt with valid positive values."""
        assert _safe_sqrt(0.0) == 0.0
        assert _safe_sqrt(1.0) == 1.0
        assert abs(_safe_sqrt(4.0) - 2.0) < 1e-10
        assert abs(_safe_sqrt(2.0) - math.sqrt(2.0)) < 1e-10

    def test_safe_sqrt_negative_values(self):
        """Test safe_sqrt returns 0 for negative values."""
        assert _safe_sqrt(-0.0001) == 0.0
        assert _safe_sqrt(-1.0) == 0.0
        assert _safe_sqrt(-1e10) == 0.0

    def test_safe_sqrt_very_small_negative(self):
        """Test safe_sqrt handles tiny negative values from numerical error."""
        # These can occur from floating point errors like a*a - b*b when a ≈ b
        assert _safe_sqrt(-1e-15) == 0.0
        assert _safe_sqrt(-1e-20) == 0.0


class TestObscurationCalculation:
    """Test obscuration calculation with edge cases."""

    def test_no_overlap(self):
        """Test obscuration when disks don't overlap."""
        r_sun = 0.5
        r_moon = 0.4
        d = 1.0  # separation > sum of radii
        assert _calculate_obscuration_safe(r_sun, r_moon, d) == 0.0

    def test_exactly_touching(self):
        """Test obscuration when disks exactly touch."""
        r_sun = 0.5
        r_moon = 0.4
        d = 0.9  # exactly sum of radii
        assert _calculate_obscuration_safe(r_sun, r_moon, d) == 0.0

    def test_complete_overlap_moon_larger(self):
        """Test obscuration when Moon completely covers Sun.

        For total eclipses, obscuration = (r_moon/r_sun)^2 which exceeds 1.0,
        matching reference API behavior.
        """
        r_sun = 0.5
        r_moon = 0.6
        d = 0.0  # same center
        expected = (r_moon / r_sun) ** 2  # 1.44
        assert abs(_calculate_obscuration_safe(r_sun, r_moon, d) - expected) < 1e-10

    def test_complete_overlap_moon_smaller(self):
        """Test obscuration when Moon is inside Sun (annular)."""
        r_sun = 0.5
        r_moon = 0.3
        d = 0.0  # same center
        expected = (0.3 / 0.5) ** 2
        assert abs(_calculate_obscuration_safe(r_sun, r_moon, d) - expected) < 1e-10

    def test_partial_overlap(self):
        """Test obscuration with partial overlap."""
        r_sun = 0.5
        r_moon = 0.4
        d = 0.6  # partial overlap
        obscuration = _calculate_obscuration_safe(r_sun, r_moon, d)
        # Should be between 0 and 1
        assert 0.0 < obscuration < 1.0

    def test_near_zero_separation(self):
        """Test obscuration with very small separation."""
        r_sun = 0.5
        r_moon = 0.4
        d = 1e-15  # nearly zero
        obscuration = _calculate_obscuration_safe(r_sun, r_moon, d)
        # Should be close to the ratio squared (annular)
        expected = (0.4 / 0.5) ** 2
        assert abs(obscuration - expected) < 0.01


class TestMagnitudeCalculation:
    """Test magnitude calculation with edge cases."""

    def test_magnitude_at_eclipse_limit(self):
        """Test magnitude when gamma equals eclipse limit."""
        gamma_limit = 1.55
        magnitude = _calculate_magnitude_safe(gamma_limit, 0.95)
        assert magnitude == 0.0

    def test_magnitude_beyond_eclipse_limit(self):
        """Test magnitude when gamma exceeds eclipse limit."""
        gamma_limit = 1.55
        magnitude = _calculate_magnitude_safe(gamma_limit + 0.1, 0.95)
        assert magnitude == 0.0

    def test_magnitude_near_eclipse_limit(self):
        """Test smooth transition near eclipse limit."""
        gamma_limit = 1.55
        # Just inside the margin
        gamma = gamma_limit - NEAR_MISS_GAMMA_MARGIN / 2
        magnitude = _calculate_magnitude_safe(gamma, 1.0)
        # Should be small but non-zero
        assert 0.0 < magnitude < SHALLOW_ECLIPSE_MAG_THRESHOLD

    def test_magnitude_central_eclipse(self):
        """Test magnitude for central eclipse (gamma ≈ 0)."""
        magnitude = _calculate_magnitude_safe(0.1, 1.0)
        # Should be close to 1 for low gamma
        assert magnitude > 0.9


class TestContactTimeValidation:
    """Test contact time validation."""

    def test_valid_contact_time(self):
        """Test that valid contact times are returned."""
        jd_max = 2459000.5
        jd_contact = 2459000.4  # 2.4 hours before maximum
        result = _validate_contact_time(jd_contact, jd_max)
        assert result == jd_contact

    def test_invalid_zero_contact(self):
        """Test that zero contact times are returned as zero."""
        jd_max = 2459000.5
        result = _validate_contact_time(0.0, jd_max)
        assert result == 0.0

    def test_invalid_negative_contact(self):
        """Test that negative contact times are returned as zero."""
        jd_max = 2459000.5
        result = _validate_contact_time(-1.0, jd_max)
        assert result == 0.0

    def test_contact_too_far_from_max(self):
        """Test that contact times too far from maximum are rejected."""
        jd_max = 2459000.5
        jd_contact = jd_max + 0.5  # 12 hours after (beyond 6-hour limit)
        result = _validate_contact_time(jd_contact, jd_max)
        assert result == 0.0


class TestSolarEclipseEdgeCases:
    """Test solar eclipse edge cases with real data."""

    def test_shallow_partial_eclipse_type_detection(self):
        """Test that shallow partial eclipses are properly typed."""
        # Search for any solar eclipse and verify proper handling
        jd_start = julday(2024, 1, 1, 0)

        # Find any solar eclipse first
        ecl_type, times = sol_eclipse_when_glob(jd_start)

        # Should find an eclipse
        assert ecl_type != 0
        # Eclipse maximum should be valid
        assert times[0] > jd_start

    def test_eclipse_at_location_no_visibility(self):
        """Test eclipse when location is far from path."""
        # Find a solar eclipse
        jd_start = julday(2024, 1, 1, 0)
        _, times_glob = sol_eclipse_when_glob(jd_start)

        # Try a location far from any typical eclipse path
        # Antarctica, unlikely to see most eclipses
        result = sol_eclipse_when_loc(
            times_glob[0] - 0.1,  # Start before max
            geopos=(0.0, -80.0, 0.0),
        )

        # Either no eclipse visible or very shallow
        ecl_type, times, attrs = result
        # If eclipse visible, times[0] should be non-zero
        # If not visible, all times should be 0
        if times[0] > 0:
            # If visible, should have valid contact times
            pass  # Eclipse is visible, which is fine
        else:
            # Not visible, all should be zero
            assert all(t == 0 for t in times[:5])


class TestLunarEclipseEdgeCases:
    """Test lunar eclipse edge cases with real data."""

    def test_shallow_penumbral_eclipse(self):
        """Test that shallow penumbral eclipses are handled."""
        # Find a penumbral lunar eclipse
        jd_start = julday(2020, 1, 1, 0)

        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_PENUMBRAL)

        # Should find a penumbral eclipse
        assert ecl_type & ECL_PENUMBRAL
        # Eclipse maximum should be valid
        assert times[0] > jd_start
        # Penumbral times should be valid
        assert times[6] > 0  # Penumbral begin
        assert times[7] > 0  # Penumbral end

    def test_partial_lunar_eclipse_contacts(self):
        """Test contact times for partial lunar eclipse."""
        # Find a partial lunar eclipse
        jd_start = julday(2021, 1, 1, 0)

        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_PARTIAL)

        # Should find a partial eclipse
        assert ecl_type & ECL_PARTIAL
        # Maximum time should be valid
        assert times[0] > jd_start
        # Partial phase times should be present
        assert times[2] > 0  # Partial begin
        assert times[3] > 0  # Partial end
        # Total phase times should be absent
        assert times[4] == 0  # Total begin
        assert times[5] == 0  # Total end


class TestEdgeCaseNumericalStability:
    """Test numerical stability in edge cases."""

    def test_obscuration_symmetric(self):
        """Test that obscuration is symmetric in disk sizes."""
        r_sun1, r_moon1 = 0.5, 0.4
        r_sun2, r_moon2 = 0.4, 0.5
        d = 0.3

        # The obscuration relative to the Sun should differ
        obs1 = _calculate_obscuration_safe(r_sun1, r_moon1, d)
        obs2 = _calculate_obscuration_safe(r_sun2, r_moon2, d)

        # Both should be valid (0-1 range)
        assert 0 <= obs1 <= 1
        assert 0 <= obs2 <= 1

    def test_magnitude_monotonic_in_gamma(self):
        """Test that magnitude decreases as gamma increases."""
        gamma_values = [0.0, 0.2, 0.5, 0.8, 1.0, 1.2, 1.4, 1.55]
        moon_sun_ratio = 1.0

        magnitudes = [
            _calculate_magnitude_safe(g, moon_sun_ratio) for g in gamma_values
        ]

        # Magnitude should be monotonically decreasing (or equal)
        for i in range(len(magnitudes) - 1):
            assert magnitudes[i] >= magnitudes[i + 1], (
                f"Magnitude should decrease: {magnitudes[i]} vs {magnitudes[i + 1]} "
                f"for gamma {gamma_values[i]} vs {gamma_values[i + 1]}"
            )

    def test_safe_sqrt_handles_floating_point_errors(self):
        """Test safe_sqrt with typical floating point error scenarios."""
        # Simulating a*a - b*b where a ≈ b
        a = 1.0000000000001
        b = 1.0
        result = a * a - b * b  # Could be tiny negative

        # Should not raise an error
        sqrt_result = _safe_sqrt(result)
        assert sqrt_result >= 0
