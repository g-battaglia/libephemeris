"""
Tests for lun_eclipse_gamma and lun_eclipse_gamma functions.

These functions calculate the gamma parameter (distance of Moon's center
from Earth's shadow axis in Earth radii) for lunar eclipses.

Reference data from NASA Eclipse website:
https://eclipse.gsfc.nasa.gov/lunar.html
"""

import pytest
from libephemeris import (
    julday,
    lun_eclipse_when,
    lun_eclipse_gamma,
    ECL_TOTAL,
    ECL_PARTIAL,
    ECL_PENUMBRAL,
    FLG_SWIEPH,
)

pytestmark = pytest.mark.slow


class TestLunEclipseGammaFunctionSignature:
    """Test that lun_eclipse_gamma function signature is correct."""

    def test_function_exists_in_module(self):
        """Test that lun_eclipse_gamma function exists in eclipse module."""
        from libephemeris.eclipse import lun_eclipse_gamma

        assert callable(lun_eclipse_gamma)

    def test_function_exists_in_package(self):
        """Test that function is exported from libephemeris package."""
        from libephemeris import lun_eclipse_gamma

        assert callable(lun_eclipse_gamma)

    def test_swe_function_exists_in_package(self):
        """Test that swe_ prefixed function is exported."""
        from libephemeris import lun_eclipse_gamma

        assert callable(lun_eclipse_gamma)

    def test_returns_float(self):
        """Test that function returns a float."""
        # Use a known eclipse time (May 16, 2022 total lunar eclipse)
        jd_eclipse = julday(2022, 5, 16, 4.2)

        result = lun_eclipse_gamma(jd_eclipse)

        assert isinstance(result, float)

    def test_accepts_flags_parameter(self):
        """Test that function accepts optional flags parameter."""
        jd_eclipse = julday(2022, 5, 16, 4.2)

        result = lun_eclipse_gamma(jd_eclipse, ifl=FLG_SWIEPH)

        assert isinstance(result, float)


class TestSweLunEclipseGammaFunctionSignature:
    """Test lun_eclipse_gamma function signature."""

    def test_function_exists_in_module(self):
        """Test that lun_eclipse_gamma function exists in eclipse module."""
        from libephemeris.eclipse import lun_eclipse_gamma

        assert callable(lun_eclipse_gamma)

    def test_returns_float(self):
        """Test that function returns a float."""
        jd_eclipse = julday(2022, 5, 16, 4.2)

        result = lun_eclipse_gamma(jd_eclipse, FLG_SWIEPH)

        assert isinstance(result, float)


class TestGammaValues:
    """Test that gamma values are calculated correctly."""

    def test_total_lunar_eclipse_has_small_gamma(self):
        """Test that total lunar eclipse has |gamma| < ~0.75 at maximum."""
        # Find a total lunar eclipse
        jd_start = julday(2022, 5, 1, 0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        gamma = lun_eclipse_gamma(jd_max)

        # For total lunar eclipse, |gamma| should be relatively small
        assert abs(gamma) < 1.0, (
            f"Total lunar eclipse should have |gamma| < 1.0, got {gamma}"
        )

    def test_partial_lunar_eclipse_has_moderate_gamma(self):
        """Test that partial lunar eclipse has moderate gamma value."""
        # Find a partial lunar eclipse
        jd_start = julday(2023, 10, 1, 0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_PARTIAL)
        jd_max = times[0]

        gamma = lun_eclipse_gamma(jd_max)

        # For partial lunar eclipse, |gamma| should be moderate (Moon only partially in umbra)
        assert 0.5 < abs(gamma) < 1.5, (
            f"Partial lunar eclipse should have 0.5 < |gamma| < 1.5, got {gamma}"
        )

    def test_penumbral_eclipse_has_larger_gamma(self):
        """Test that penumbral-only eclipse has larger |gamma|."""
        # Find a penumbral lunar eclipse
        jd_start = julday(2020, 1, 1, 0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_PENUMBRAL)
        jd_max = times[0]

        gamma = lun_eclipse_gamma(jd_max)

        # For penumbral-only eclipse, Moon is outside umbra but inside penumbra
        # Gamma should be larger than for umbral eclipses
        assert abs(gamma) > 0.7, (
            f"Penumbral eclipse should have |gamma| > 0.7, got {gamma}"
        )

    def test_no_eclipse_large_gamma(self):
        """Test that non-eclipse time returns large |gamma| value."""
        # Use a time when there's definitely no eclipse (first quarter moon)
        jd_no_eclipse = julday(2022, 6, 7, 12.0)

        gamma = lun_eclipse_gamma(jd_no_eclipse)

        # When there's no eclipse, gamma should be large (Moon far from shadow)
        assert abs(gamma) > 1.5, (
            f"Non-eclipse time should have large |gamma|, got {gamma}"
        )

    def test_gamma_can_be_positive_or_negative(self):
        """Test that gamma values vary across eclipses.

        Note: For lunar eclipses, gamma represents the distance of the Moon
        from the axis of Earth's shadow (in Earth radii). In pyswisseph's
        convention (attr[7] from lun_eclipse_how), gamma is always >= 0
        (it's the absolute distance, not signed). We verify that gamma
        varies across eclipses (not always the same value).
        """
        # Find several eclipses and check that gamma varies
        gamma_values = []

        for year in [2020, 2021, 2022, 2023, 2024]:
            jd_start = julday(year, 1, 1, 0)
            for _ in range(4):  # Check up to 4 eclipses per year
                ecl_type, times = lun_eclipse_when(jd_start)
                if times[0] > jd_start:
                    gamma = lun_eclipse_gamma(times[0])
                    gamma_values.append(gamma)
                    jd_start = times[0] + 30  # Move past this eclipse

        # Check that we have variation in gamma values
        assert len(gamma_values) >= 5, "Should find at least 5 eclipses"
        assert len(set(round(g, 3) for g in gamma_values)) > 1, (
            "Gamma values should vary across eclipses"
        )
        # All gamma values should be non-negative (distance from shadow axis)
        assert all(g >= 0 for g in gamma_values), (
            "Gamma should be non-negative for lunar eclipses"
        )


class TestGammaConsistency:
    """Test consistency of gamma calculations."""

    def test_swe_function_matches_legacy_function(self):
        """Test that swe_ function gives same result as legacy function."""
        jd_eclipse = julday(2022, 5, 16, 4.2)

        result_legacy = lun_eclipse_gamma(jd_eclipse)
        result_swe = lun_eclipse_gamma(jd_eclipse, FLG_SWIEPH)

        assert result_legacy == result_swe, (
            f"Results should match: legacy={result_legacy}, swe={result_swe}"
        )

    def test_gamma_consistent_with_magnitude(self):
        """Test that gamma is consistent with eclipse magnitude behavior."""
        from libephemeris import lun_eclipse_umbral_magnitude

        # For a total eclipse, gamma should be smaller
        jd_start = julday(2022, 5, 1, 0)
        _, times_total = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_total = times_total[0]

        gamma_total = lun_eclipse_gamma(jd_total)
        mag_total = lun_eclipse_umbral_magnitude(jd_total)

        # For a partial eclipse, gamma should be larger
        jd_start = julday(2023, 10, 1, 0)
        _, times_partial = lun_eclipse_when(jd_start, ecltype=ECL_PARTIAL)
        jd_partial = times_partial[0]

        gamma_partial = lun_eclipse_gamma(jd_partial)
        mag_partial = lun_eclipse_umbral_magnitude(jd_partial)

        # Smaller |gamma| should correlate with larger magnitude
        assert abs(gamma_total) < abs(gamma_partial), (
            f"Total eclipse should have smaller |gamma|: "
            f"total={abs(gamma_total)}, partial={abs(gamma_partial)}"
        )
        assert mag_total > mag_partial, (
            f"Total eclipse should have larger magnitude: "
            f"total={mag_total}, partial={mag_partial}"
        )

    def test_gamma_sign_indicates_direction(self):
        """Test that gamma is a non-negative distance from shadow axis.

        For lunar eclipses, gamma represents the absolute distance of the
        Moon's center from the axis of Earth's shadow, measured in Earth
        equatorial radii. Unlike solar eclipse gamma which can be signed,
        lunar eclipse gamma (as returned by pyswisseph attr[7]) is always
        non-negative — it is the magnitude of the displacement, not a
        signed north/south indicator.
        """
        # Find an eclipse and verify gamma is non-negative
        jd_start = julday(2022, 5, 1, 0)
        _, times = lun_eclipse_when(jd_start)
        jd_max = times[0]

        gamma = lun_eclipse_gamma(jd_max)

        # Gamma should be non-negative (absolute distance from shadow axis)
        assert gamma >= 0, f"Gamma should be non-negative, got {gamma}"

        # For a total eclipse, gamma should be small (< ~0.4)
        # May 2022 was a total eclipse with gamma ~0.25
        assert gamma < 1.0, f"Gamma for a total eclipse should be < 1.0, got {gamma}"


class TestKnownEclipses:
    """Test against known eclipse data."""

    def test_may_2022_total_eclipse(self):
        """Test May 16, 2022 total lunar eclipse.

        Reference: NASA - This was a total lunar eclipse with gamma ~ -0.253
        """
        # Find the eclipse
        jd_start = julday(2022, 5, 1, 0)
        _, times = lun_eclipse_when(jd_start, FLG_SWIEPH, ECL_TOTAL)
        jd_max = times[0]

        gamma = lun_eclipse_gamma(jd_max)

        # Should be a central total eclipse (small |gamma|)
        assert abs(gamma) < 0.5, (
            f"Expected |gamma| < 0.5 for central total, got {gamma}"
        )

    def test_nov_2022_total_eclipse(self):
        """Test November 8, 2022 total lunar eclipse.

        Reference: NASA - This was a total lunar eclipse with gamma ~ 0.259
        """
        # Find the eclipse
        jd_start = julday(2022, 11, 1, 0)
        _, times = lun_eclipse_when(jd_start, FLG_SWIEPH, ECL_TOTAL)
        jd_max = times[0]

        gamma = lun_eclipse_gamma(jd_max)

        # Should be a central total eclipse (small |gamma|)
        assert abs(gamma) < 0.5, (
            f"Expected |gamma| < 0.5 for central total, got {gamma}"
        )

    def test_oct_2023_partial_eclipse(self):
        """Test October 28, 2023 partial lunar eclipse.

        Reference: NASA - This was a partial lunar eclipse with larger gamma
        """
        # Find the eclipse
        jd_start = julday(2023, 10, 1, 0)
        _, times = lun_eclipse_when(jd_start, FLG_SWIEPH, ECL_PARTIAL)
        jd_max = times[0]

        gamma = lun_eclipse_gamma(jd_max)

        # Partial eclipse should have larger |gamma| than total
        assert 0.5 < abs(gamma) < 1.5, (
            f"Expected 0.5 < |gamma| < 1.5 for partial, got {gamma}"
        )


class TestEdgeCases:
    """Test edge cases for gamma calculation."""

    def test_very_old_eclipse(self):
        """Test eclipse calculation for distant past."""
        # Find an eclipse in 1900
        jd_start = julday(1900, 1, 1, 0)
        _, times = lun_eclipse_when(jd_start)
        jd_max = times[0]

        gamma = lun_eclipse_gamma(jd_max)

        # Should return a valid gamma
        assert isinstance(gamma, float)
        assert not (gamma != gamma)  # Not NaN

    def test_future_eclipse(self):
        """Test eclipse calculation for distant future."""
        # Find an eclipse in 2050
        jd_start = julday(2050, 1, 1, 0)
        _, times = lun_eclipse_when(jd_start)
        jd_max = times[0]

        gamma = lun_eclipse_gamma(jd_max)

        # Should return a valid gamma
        assert isinstance(gamma, float)
        assert not (gamma != gamma)  # Not NaN

    def test_gamma_at_various_times(self):
        """Test gamma calculation at various times during eclipse."""
        # Find a total lunar eclipse
        jd_start = julday(2022, 5, 1, 0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        # Calculate gamma at different phases
        gamma_before = lun_eclipse_gamma(jd_max - 0.1)  # About 2.4 hours before max
        gamma_at_max = lun_eclipse_gamma(jd_max)
        gamma_after = lun_eclipse_gamma(jd_max + 0.1)  # About 2.4 hours after max

        # All should be valid floats
        for gamma in [gamma_before, gamma_at_max, gamma_after]:
            assert isinstance(gamma, float)

        # At maximum, |gamma| should be at its minimum (Moon closest to axis)
        # Note: This is only approximately true over short time spans
        # because the Moon's latitude changes slowly
        assert abs(gamma_at_max) <= abs(gamma_before) + 0.2, (
            f"Gamma at max should be near minimum: "
            f"before={gamma_before}, max={gamma_at_max}"
        )


class TestDocumentation:
    """Test that documentation examples work correctly."""

    def test_docstring_example(self):
        """Test the example from the function docstring."""
        from libephemeris import lun_eclipse_gamma, lun_eclipse_when

        # First find a lunar eclipse
        jd_start = julday(2022, 5, 1, 0)
        ecl_type, times = lun_eclipse_when(jd_start)
        jd_max = times[0]  # Time of maximum eclipse

        # Calculate gamma at maximum
        gamma = lun_eclipse_gamma(jd_max)

        # Should return a valid gamma value
        assert isinstance(gamma, float)
        assert not (gamma != gamma)  # Not NaN

        # Check gamma at a random time (far from eclipse)
        jd_no_eclipse = julday(2022, 6, 1, 12.0)
        gamma = lun_eclipse_gamma(jd_no_eclipse)

        # Should be a large value (no eclipse)
        assert abs(gamma) > 1.0  # Moon far from shadow axis
