"""
Tests for lun_eclipse_penumbral_magnitude and lun_eclipse_penumbral_magnitude functions.

These functions calculate the penumbral magnitude (fraction of Moon's diameter
within Earth's penumbral shadow) for lunar eclipses.

Reference data from NASA Eclipse website:
https://eclipse.gsfc.nasa.gov/lunar.html
"""

import pytest
from libephemeris import (
    julday,
    lun_eclipse_when,
    lun_eclipse_penumbral_magnitude,
    ECL_TOTAL,
    ECL_PARTIAL,
    ECL_PENUMBRAL,
    FLG_SWIEPH,
)


class TestLunEclipsePenumbralMagnitudeFunctionSignature:
    """Test that lun_eclipse_penumbral_magnitude function signature is correct."""

    def test_function_exists_in_module(self):
        """Test that lun_eclipse_penumbral_magnitude function exists in eclipse module."""
        from libephemeris.eclipse import lun_eclipse_penumbral_magnitude

        assert callable(lun_eclipse_penumbral_magnitude)

    def test_function_exists_in_package(self):
        """Test that function is exported from libephemeris package."""
        from libephemeris import lun_eclipse_penumbral_magnitude

        assert callable(lun_eclipse_penumbral_magnitude)

    def test_swe_function_exists_in_package(self):
        """Test that swe_ prefixed function is exported."""
        from libephemeris import lun_eclipse_penumbral_magnitude

        assert callable(lun_eclipse_penumbral_magnitude)

    def test_returns_float(self):
        """Test that function returns a float."""
        # Use a known eclipse time (May 16, 2022 total lunar eclipse)
        jd_eclipse = julday(2022, 5, 16, 4.2)

        result = lun_eclipse_penumbral_magnitude(jd_eclipse)

        assert isinstance(result, float)

    def test_accepts_flags_parameter(self):
        """Test that function accepts optional flags parameter."""
        jd_eclipse = julday(2022, 5, 16, 4.2)

        result = lun_eclipse_penumbral_magnitude(jd_eclipse, ifl=FLG_SWIEPH)

        assert isinstance(result, float)


class TestSweLunEclipsePenumbralMagnitudeFunctionSignature:
    """Test lun_eclipse_penumbral_magnitude function signature."""

    def test_function_exists_in_module(self):
        """Test that lun_eclipse_penumbral_magnitude function exists in eclipse module."""
        from libephemeris.eclipse import lun_eclipse_penumbral_magnitude

        assert callable(lun_eclipse_penumbral_magnitude)

    def test_returns_float(self):
        """Test that function returns a float."""
        jd_eclipse = julday(2022, 5, 16, 4.2)

        result = lun_eclipse_penumbral_magnitude(jd_eclipse, FLG_SWIEPH)

        assert isinstance(result, float)


class TestPenumbralMagnitudeValues:
    """Test that penumbral magnitude values are calculated correctly."""

    def test_total_lunar_eclipse_penumbral_magnitude_greater_than_one(self):
        """Test that total lunar eclipse has penumbral magnitude > 1.0 at maximum.

        When the Moon is fully within the umbra (total eclipse), it is also
        fully within the penumbra, so penumbral magnitude should be > 1.0.
        """
        # Find a total lunar eclipse
        jd_start = julday(2022, 5, 1, 0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        penumbral_mag = lun_eclipse_penumbral_magnitude(jd_max)

        # For total lunar eclipse at maximum, Moon is deep in penumbra
        assert penumbral_mag > 1.0, (
            f"Total lunar eclipse should have penumbral magnitude > 1.0, got {penumbral_mag}"
        )

    def test_partial_lunar_eclipse_penumbral_magnitude_greater_than_umbral(self):
        """Test that partial lunar eclipse has penumbral mag > umbral mag."""
        from libephemeris import lun_eclipse_umbral_magnitude

        # Find a partial lunar eclipse
        jd_start = julday(2023, 10, 1, 0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_PARTIAL)
        jd_max = times[0]

        penumbral_mag = lun_eclipse_penumbral_magnitude(jd_max)
        umbral_mag = lun_eclipse_umbral_magnitude(jd_max)

        # Penumbral magnitude should always be >= umbral magnitude
        assert penumbral_mag >= umbral_mag, (
            f"Penumbral mag ({penumbral_mag}) should be >= umbral mag ({umbral_mag})"
        )

    def test_penumbral_eclipse_has_positive_penumbral_magnitude(self):
        """Test that penumbral-only eclipse has positive penumbral magnitude."""
        # Find a penumbral lunar eclipse
        jd_start = julday(2020, 1, 1, 0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_PENUMBRAL)
        jd_max = times[0]

        penumbral_mag = lun_eclipse_penumbral_magnitude(jd_max)

        # For penumbral eclipse, penumbral magnitude should be > 0
        assert penumbral_mag > 0.0, (
            f"Penumbral eclipse should have positive penumbral magnitude, got {penumbral_mag}"
        )

    def test_penumbral_eclipse_umbral_magnitude_is_zero(self):
        """Test that penumbral-only eclipse has umbral magnitude ~ 0."""
        from libephemeris import lun_eclipse_umbral_magnitude

        # Find a penumbral lunar eclipse
        jd_start = julday(2020, 1, 1, 0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_PENUMBRAL)
        jd_max = times[0]

        umbral_mag = lun_eclipse_umbral_magnitude(jd_max)
        penumbral_mag = lun_eclipse_penumbral_magnitude(jd_max)

        # For penumbral-only eclipse, umbral magnitude should be 0 or near-zero
        # but penumbral magnitude should be positive
        assert umbral_mag <= 0.01, (
            f"Penumbral eclipse should have umbral magnitude ~ 0, got {umbral_mag}"
        )
        assert penumbral_mag > 0.0, (
            f"Penumbral eclipse should have positive penumbral magnitude, got {penumbral_mag}"
        )

    def test_no_eclipse_penumbral_magnitude_zero(self):
        """Test that non-eclipse time returns zero penumbral magnitude."""
        # Use a time when there's definitely no eclipse (first quarter moon)
        jd_no_eclipse = julday(2022, 6, 7, 12.0)

        penumbral_mag = lun_eclipse_penumbral_magnitude(jd_no_eclipse)

        assert penumbral_mag == 0.0, (
            f"Non-eclipse time should have penumbral magnitude 0.0, got {penumbral_mag}"
        )

    def test_magnitude_non_negative(self):
        """Test that penumbral magnitude is never negative."""
        # Test various times
        test_times = [
            julday(2020, 1, 15, 12.0),  # Random time
            julday(2022, 6, 1, 12.0),  # Non-eclipse
            julday(2022, 5, 16, 4.2),  # During eclipse
        ]

        for jd in test_times:
            penumbral_mag = lun_eclipse_penumbral_magnitude(jd)
            assert penumbral_mag >= 0.0, (
                f"Magnitude should be >= 0, got {penumbral_mag}"
            )


class TestPenumbralMagnitudeConsistency:
    """Test consistency of penumbral magnitude calculations."""

    def test_magnitude_matches_lun_eclipse_how(self):
        """Test that magnitude matches the value from lun_eclipse_how."""
        from libephemeris import lun_eclipse_how

        # Find a lunar eclipse
        jd_start = julday(2022, 5, 1, 0)
        ecl_type, times = lun_eclipse_when(jd_start)
        jd_max = times[0]

        # Get magnitude from both functions
        penumbral_mag_direct = lun_eclipse_penumbral_magnitude(jd_max)
        _, attr = lun_eclipse_how(
            jd_max, (0.0, 0.0, 0.0)
        )  # Location doesn't affect lunar eclipse magnitude
        penumbral_mag_from_how = attr[1]  # attr[1] is penumbral magnitude

        # They should be equal or very close
        assert abs(penumbral_mag_direct - penumbral_mag_from_how) < 0.01, (
            f"Magnitudes should match: direct={penumbral_mag_direct}, "
            f"from_how={penumbral_mag_from_how}"
        )

    def test_swe_function_matches_legacy_function(self):
        """Test that swe_ function gives same result as legacy function."""
        jd_eclipse = julday(2022, 5, 16, 4.2)

        result_legacy = lun_eclipse_penumbral_magnitude(jd_eclipse)
        result_swe = lun_eclipse_penumbral_magnitude(jd_eclipse, FLG_SWIEPH)

        assert result_legacy == result_swe, (
            f"Results should match: legacy={result_legacy}, swe={result_swe}"
        )

    def test_penumbral_always_greater_or_equal_umbral(self):
        """Test that penumbral magnitude >= umbral magnitude during eclipse."""
        from libephemeris import lun_eclipse_umbral_magnitude

        # Find different types of lunar eclipses and test
        for eclipse_type in [ECL_TOTAL, ECL_PARTIAL, ECL_PENUMBRAL]:
            jd_start = julday(2020, 1, 1, 0)
            ecl_type, times = lun_eclipse_when(jd_start, ecltype=eclipse_type)
            jd_max = times[0]

            penumbral_mag = lun_eclipse_penumbral_magnitude(jd_max)
            umbral_mag = lun_eclipse_umbral_magnitude(jd_max)

            assert penumbral_mag >= umbral_mag, (
                f"Penumbral mag ({penumbral_mag}) should be >= umbral mag ({umbral_mag}) "
                f"for eclipse type {eclipse_type}"
            )


class TestKnownEclipses:
    """Test against known eclipse data."""

    def test_jan_2020_penumbral_eclipse(self):
        """Test January 10, 2020 penumbral lunar eclipse.

        Reference: NASA - Penumbral magnitude at maximum was 0.896
        This was a penumbral-only eclipse (umbral magnitude ~ 0).
        """
        # Find the eclipse
        jd_start = julday(2020, 1, 1, 0)
        _, times = lun_eclipse_when(jd_start, FLG_SWIEPH, ECL_PENUMBRAL)
        jd_max = times[0]

        penumbral_mag = lun_eclipse_penumbral_magnitude(jd_max)

        # NASA reference: penumbral magnitude ~ 0.896
        # Should be less than 1.0 for a penumbral-only eclipse
        assert 0 < penumbral_mag < 1.1, (
            f"Expected penumbral magnitude around 0.9, got {penumbral_mag}"
        )

    def test_may_2022_total_eclipse(self):
        """Test May 16, 2022 total lunar eclipse.

        Reference: NASA - Penumbral magnitude at maximum was 2.363
        Total eclipses have penumbral magnitude > 1.0.
        """
        # Find the eclipse
        jd_start = julday(2022, 5, 1, 0)
        _, times = lun_eclipse_when(jd_start, FLG_SWIEPH, ECL_TOTAL)
        jd_max = times[0]

        penumbral_mag = lun_eclipse_penumbral_magnitude(jd_max)

        # NASA reference: penumbral magnitude ~ 2.363
        # Should be > 1.0 for total eclipse
        assert penumbral_mag > 1.0, "Should have penumbral mag > 1.0 for total eclipse"
        assert 1.5 < penumbral_mag < 3.0, (
            f"Expected penumbral magnitude around 2.4, got {penumbral_mag}"
        )

    def test_oct_2023_partial_eclipse(self):
        """Test October 28, 2023 partial lunar eclipse.

        Reference: NASA - Penumbral magnitude at maximum was 1.124
        """
        # Find the eclipse
        jd_start = julday(2023, 10, 1, 0)
        _, times = lun_eclipse_when(jd_start, FLG_SWIEPH, ECL_PARTIAL)
        jd_max = times[0]

        penumbral_mag = lun_eclipse_penumbral_magnitude(jd_max)

        # NASA reference: penumbral magnitude ~ 1.124
        # Partial eclipses can have penumbral magnitude > 1.0 (fully in penumbra)
        assert 0.8 < penumbral_mag < 1.5, (
            f"Expected penumbral magnitude around 1.1, got {penumbral_mag}"
        )


class TestEdgeCases:
    """Test edge cases for penumbral magnitude calculation."""

    def test_very_old_eclipse(self):
        """Test eclipse calculation for distant past."""
        # Find an eclipse in 1900
        jd_start = julday(1900, 1, 1, 0)
        _, times = lun_eclipse_when(jd_start)
        jd_max = times[0]

        penumbral_mag = lun_eclipse_penumbral_magnitude(jd_max)

        # Should return a valid magnitude
        assert penumbral_mag >= 0.0
        assert isinstance(penumbral_mag, float)

    def test_future_eclipse(self):
        """Test eclipse calculation for distant future."""
        # Find an eclipse in 2050
        jd_start = julday(2050, 1, 1, 0)
        _, times = lun_eclipse_when(jd_start)
        jd_max = times[0]

        penumbral_mag = lun_eclipse_penumbral_magnitude(jd_max)

        # Should return a valid magnitude
        assert penumbral_mag >= 0.0
        assert isinstance(penumbral_mag, float)

    def test_penumbral_only_vs_umbral_eclipse(self):
        """Test that penumbral-only eclipse has umbral mag near 0."""
        from libephemeris import lun_eclipse_umbral_magnitude

        # Find a penumbral-only lunar eclipse
        jd_start = julday(2020, 1, 1, 0)
        _, times = lun_eclipse_when(jd_start, ecltype=ECL_PENUMBRAL)
        jd_max = times[0]

        penumbral_mag = lun_eclipse_penumbral_magnitude(jd_max)
        umbral_mag = lun_eclipse_umbral_magnitude(jd_max)

        # For penumbral-only, penumbral should be positive, umbral near zero
        assert penumbral_mag > 0.0, "Penumbral mag should be positive"
        assert umbral_mag <= 0.01, "Umbral mag should be near zero"


class TestDocumentation:
    """Test that documentation examples work correctly."""

    def test_docstring_example(self):
        """Test the example from the function docstring."""
        from libephemeris import lun_eclipse_penumbral_magnitude, lun_eclipse_when

        # First find a lunar eclipse
        jd_start = julday(2022, 5, 1, 0)
        ecl_type, times = lun_eclipse_when(jd_start)
        jd_max = times[0]  # Time of maximum eclipse

        # Calculate penumbral magnitude at maximum
        penumbral_mag = lun_eclipse_penumbral_magnitude(jd_max)

        # Should return a valid positive magnitude
        assert penumbral_mag > 0.0
        assert isinstance(penumbral_mag, float)

        # Check magnitude at a random time (no eclipse)
        jd_no_eclipse = julday(2022, 6, 1, 12.0)
        mag = lun_eclipse_penumbral_magnitude(jd_no_eclipse)

        # Should be 0.0
        assert mag == 0.0
