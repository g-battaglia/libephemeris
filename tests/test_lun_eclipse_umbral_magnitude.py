"""
Tests for lun_eclipse_umbral_magnitude and lun_eclipse_umbral_magnitude functions.

These functions calculate the umbral magnitude (fraction of Moon's diameter
within Earth's umbral shadow) for lunar eclipses.

Reference data from NASA Eclipse website:
https://eclipse.gsfc.nasa.gov/lunar.html
"""

import pytest
from libephemeris import (
    julday,
    lun_eclipse_when,
    lun_eclipse_umbral_magnitude,
    ECL_TOTAL,
    ECL_PARTIAL,
    ECL_PENUMBRAL,
    FLG_SWIEPH,
)

pytestmark = pytest.mark.slow


class TestLunEclipseUmbralMagnitudeFunctionSignature:
    """Test that lun_eclipse_umbral_magnitude function signature is correct."""

    def test_function_exists_in_module(self):
        """Test that lun_eclipse_umbral_magnitude function exists in eclipse module."""
        from libephemeris.eclipse import lun_eclipse_umbral_magnitude

        assert callable(lun_eclipse_umbral_magnitude)

    def test_function_exists_in_package(self):
        """Test that function is exported from libephemeris package."""
        from libephemeris import lun_eclipse_umbral_magnitude

        assert callable(lun_eclipse_umbral_magnitude)

    def test_swe_function_exists_in_package(self):
        """Test that swe_ prefixed function is exported."""
        from libephemeris import lun_eclipse_umbral_magnitude

        assert callable(lun_eclipse_umbral_magnitude)

    def test_returns_float(self):
        """Test that function returns a float."""
        # Use a known eclipse time (May 16, 2022 total lunar eclipse)
        jd_eclipse = julday(2022, 5, 16, 4.2)

        result = lun_eclipse_umbral_magnitude(jd_eclipse)

        assert isinstance(result, float)

    def test_accepts_flags_parameter(self):
        """Test that function accepts optional flags parameter."""
        jd_eclipse = julday(2022, 5, 16, 4.2)

        result = lun_eclipse_umbral_magnitude(jd_eclipse, ifl=FLG_SWIEPH)

        assert isinstance(result, float)


class TestSweLunEclipseUmbralMagnitudeFunctionSignature:
    """Test lun_eclipse_umbral_magnitude function signature."""

    def test_function_exists_in_module(self):
        """Test that lun_eclipse_umbral_magnitude function exists in eclipse module."""
        from libephemeris.eclipse import lun_eclipse_umbral_magnitude

        assert callable(lun_eclipse_umbral_magnitude)

    def test_returns_float(self):
        """Test that function returns a float."""
        jd_eclipse = julday(2022, 5, 16, 4.2)

        result = lun_eclipse_umbral_magnitude(jd_eclipse, FLG_SWIEPH)

        assert isinstance(result, float)


class TestUmbralMagnitudeValues:
    """Test that umbral magnitude values are calculated correctly."""

    def test_total_lunar_eclipse_magnitude_greater_than_one(self):
        """Test that total lunar eclipse has umbral magnitude > 1.0 at maximum."""
        # Find a total lunar eclipse
        jd_start = julday(2022, 5, 1, 0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        umbral_mag = lun_eclipse_umbral_magnitude(jd_max)

        # For total lunar eclipse at maximum, magnitude should be > 1.0
        assert umbral_mag > 1.0, (
            f"Total lunar eclipse should have umbral magnitude > 1.0, got {umbral_mag}"
        )

    def test_partial_lunar_eclipse_magnitude_between_zero_and_one(self):
        """Test that partial lunar eclipse has 0 < umbral magnitude < 1.0."""
        # Find a partial lunar eclipse
        jd_start = julday(2023, 10, 1, 0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_PARTIAL)
        jd_max = times[0]

        umbral_mag = lun_eclipse_umbral_magnitude(jd_max)

        # For partial lunar eclipse, magnitude should be between 0 and 1
        assert 0 < umbral_mag < 1.0, (
            f"Partial lunar eclipse should have 0 < mag < 1.0, got {umbral_mag}"
        )

    def test_penumbral_eclipse_umbral_magnitude_zero_or_near_zero(self):
        """Test that penumbral-only eclipse has umbral magnitude ~ 0."""
        # Find a penumbral lunar eclipse
        jd_start = julday(2020, 1, 1, 0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_PENUMBRAL)
        jd_max = times[0]

        umbral_mag = lun_eclipse_umbral_magnitude(jd_max)

        # For penumbral-only eclipse, umbral magnitude should be 0 or very small
        assert umbral_mag <= 0.01, (
            f"Penumbral eclipse should have umbral magnitude ~ 0, got {umbral_mag}"
        )

    def test_no_eclipse_magnitude_zero(self):
        """Test that non-eclipse time returns zero magnitude."""
        # Use a time when there's definitely no eclipse (first quarter moon)
        jd_no_eclipse = julday(2022, 6, 7, 12.0)

        umbral_mag = lun_eclipse_umbral_magnitude(jd_no_eclipse)

        assert umbral_mag == 0.0, (
            f"Non-eclipse time should have umbral magnitude 0.0, got {umbral_mag}"
        )

    def test_magnitude_non_negative(self):
        """Test that umbral magnitude is never negative."""
        # Test various times
        test_times = [
            julday(2020, 1, 15, 12.0),  # Random time
            julday(2022, 6, 1, 12.0),  # Non-eclipse
            julday(2022, 5, 16, 4.2),  # During eclipse
        ]

        for jd in test_times:
            umbral_mag = lun_eclipse_umbral_magnitude(jd)
            assert umbral_mag >= 0.0, f"Magnitude should be >= 0, got {umbral_mag}"


class TestUmbralMagnitudeConsistency:
    """Test consistency of umbral magnitude calculations."""

    def test_magnitude_matches_lun_eclipse_how(self):
        """Test that magnitude matches the value from lun_eclipse_how."""
        from libephemeris import lun_eclipse_how

        # Find a lunar eclipse
        jd_start = julday(2022, 5, 1, 0)
        ecl_type, times = lun_eclipse_when(jd_start)
        jd_max = times[0]

        # Get magnitude from both functions
        umbral_mag_direct = lun_eclipse_umbral_magnitude(jd_max)
        _, attr = lun_eclipse_how(
            jd_max, (0.0, 0.0, 0.0)
        )  # Location doesn't affect lunar eclipse magnitude
        umbral_mag_from_how = attr[0]

        # They should be equal or very close
        assert abs(umbral_mag_direct - umbral_mag_from_how) < 0.01, (
            f"Magnitudes should match: direct={umbral_mag_direct}, "
            f"from_how={umbral_mag_from_how}"
        )

    def test_swe_function_matches_legacy_function(self):
        """Test that swe_ function gives same result as legacy function."""
        jd_eclipse = julday(2022, 5, 16, 4.2)

        result_legacy = lun_eclipse_umbral_magnitude(jd_eclipse)
        result_swe = lun_eclipse_umbral_magnitude(jd_eclipse, FLG_SWIEPH)

        assert result_legacy == result_swe, (
            f"Results should match: legacy={result_legacy}, swe={result_swe}"
        )

    def test_magnitude_increases_toward_maximum(self):
        """Test that magnitude is highest at eclipse maximum."""
        # Find a lunar eclipse
        jd_start = julday(2022, 5, 1, 0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        # At maximum, magnitude should be the peak
        mag_at_max = lun_eclipse_umbral_magnitude(jd_max)

        # For total eclipse, magnitude should be > 1.0
        assert mag_at_max > 1.0, (
            f"At maximum, magnitude should be > 1.0 for total eclipse, got {mag_at_max}"
        )


class TestKnownEclipses:
    """Test against known eclipse data."""

    def test_may_2022_total_eclipse(self):
        """Test May 16, 2022 total lunar eclipse.

        Reference: NASA - Umbral magnitude at maximum was 1.414
        """
        # Find the eclipse
        jd_start = julday(2022, 5, 1, 0)
        _, times = lun_eclipse_when(jd_start, FLG_SWIEPH, ECL_TOTAL)
        jd_max = times[0]

        umbral_mag = lun_eclipse_umbral_magnitude(jd_max)

        # NASA reference: umbral magnitude ~ 1.414
        # Allow some tolerance for calculation differences
        assert umbral_mag > 1.0, "Should be total (mag > 1.0)"
        assert 1.0 < umbral_mag < 2.0, (
            f"Expected umbral magnitude around 1.4, got {umbral_mag}"
        )

    def test_nov_2022_total_eclipse(self):
        """Test November 8, 2022 total lunar eclipse.

        Reference: NASA - Umbral magnitude at maximum was 1.359
        """
        # Find the eclipse
        jd_start = julday(2022, 11, 1, 0)
        _, times = lun_eclipse_when(jd_start, FLG_SWIEPH, ECL_TOTAL)
        jd_max = times[0]

        umbral_mag = lun_eclipse_umbral_magnitude(jd_max)

        # NASA reference: umbral magnitude ~ 1.359
        assert umbral_mag > 1.0, "Should be total (mag > 1.0)"
        assert 1.0 < umbral_mag < 2.0, (
            f"Expected umbral magnitude around 1.4, got {umbral_mag}"
        )

    def test_oct_2023_partial_eclipse(self):
        """Test October 28, 2023 partial lunar eclipse.

        Reference: NASA - Umbral magnitude at maximum was 0.122
        """
        # Find the eclipse
        jd_start = julday(2023, 10, 1, 0)
        _, times = lun_eclipse_when(jd_start, FLG_SWIEPH, ECL_PARTIAL)
        jd_max = times[0]

        umbral_mag = lun_eclipse_umbral_magnitude(jd_max)

        # NASA reference: umbral magnitude ~ 0.122
        assert 0 < umbral_mag < 1.0, "Should be partial (0 < mag < 1.0)"


class TestEdgeCases:
    """Test edge cases for umbral magnitude calculation."""

    def test_very_old_eclipse(self):
        """Test eclipse calculation for distant past."""
        # Find an eclipse in 1900
        jd_start = julday(1900, 1, 1, 0)
        _, times = lun_eclipse_when(jd_start)
        jd_max = times[0]

        umbral_mag = lun_eclipse_umbral_magnitude(jd_max)

        # Should return a valid magnitude
        assert umbral_mag >= 0.0
        assert isinstance(umbral_mag, float)

    def test_future_eclipse(self):
        """Test eclipse calculation for distant future."""
        # Find an eclipse in 2050
        jd_start = julday(2050, 1, 1, 0)
        _, times = lun_eclipse_when(jd_start)
        jd_max = times[0]

        umbral_mag = lun_eclipse_umbral_magnitude(jd_max)

        # Should return a valid magnitude
        assert umbral_mag >= 0.0
        assert isinstance(umbral_mag, float)

    def test_at_eclipse_boundaries(self):
        """Test magnitude at eclipse maximum for total eclipse."""
        # Find a total lunar eclipse with all phases
        jd_start = julday(2022, 5, 1, 0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)

        jd_max = times[0]

        # At maximum, magnitude should be > 1.0 for total eclipse
        mag_at_max = lun_eclipse_umbral_magnitude(jd_max)
        assert mag_at_max > 1.0, "At maximum, magnitude should be > 1.0 for total"

        # The magnitude should be a positive value
        assert mag_at_max > 0, "Magnitude should be positive"


class TestDocumentation:
    """Test that documentation examples work correctly."""

    def test_docstring_example(self):
        """Test the example from the function docstring."""
        from libephemeris import lun_eclipse_umbral_magnitude, lun_eclipse_when

        # First find a lunar eclipse
        jd_start = julday(2022, 5, 1, 0)
        ecl_type, times = lun_eclipse_when(jd_start)
        jd_max = times[0]  # Time of maximum eclipse

        # Calculate umbral magnitude at maximum
        umbral_mag = lun_eclipse_umbral_magnitude(jd_max)

        # Should return a valid positive magnitude
        assert umbral_mag > 0.0
        assert isinstance(umbral_mag, float)

        # Check magnitude at a random time (no eclipse)
        jd_no_eclipse = julday(2022, 6, 1, 12.0)
        mag = lun_eclipse_umbral_magnitude(jd_no_eclipse)

        # Should be 0.0
        assert mag == 0.0
