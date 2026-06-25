"""
Tests for ERFA Nutation Integration Evaluation.

This module tests the pyerfa integration and evaluates precision improvements
from using erfa.nut00a(), erfa.nut06a(), and erfa.pnm06a().

Test categories:
1. Basic functionality tests - verify functions work correctly
2. Precision comparison tests - compare ERFA vs Skyfield implementations
3. Edge case tests - extreme dates, cached values
4. Integration tests - verify the module works with existing libephemeris code
"""

import math

import pytest

# Check if pyerfa is available
try:
    import erfa

    HAS_ERFA = True
except ImportError:
    HAS_ERFA = False

from libephemeris.erfa_nutation import (
    clear_erfa_nutation_cache,
    compare_nutation_models,
    get_erfa_nutation_cache_info,
    get_erfa_nutation_cached,
    get_erfa_nutation_nut00a,
    get_erfa_nutation_nut06a,
    get_erfa_obliquity_iau2006,
    get_erfa_pnm06a_matrix,
    has_erfa,
)

# Standard test dates
J2000_JD = 2451545.0  # J2000.0 epoch (Jan 1, 2000, 12:00 TT)


class TestERFAAvailability:
    """Test ERFA availability detection."""

    def test_has_erfa_returns_bool(self):
        """has_erfa() should return a boolean."""
        result = has_erfa()
        assert isinstance(result, bool)

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_has_erfa_true_when_installed(self):
        """has_erfa() should return True when pyerfa is installed."""
        assert has_erfa() is True


class TestERFANutation:
    """Test ERFA nutation functions."""

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_nut00a_returns_tuple(self):
        """get_erfa_nutation_nut00a should return a tuple of two floats."""
        result = get_erfa_nutation_nut00a(J2000_JD)
        assert result is not None
        assert isinstance(result, tuple)
        assert len(result) == 2
        dpsi, deps = result
        assert isinstance(dpsi, float)
        assert isinstance(deps, float)

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_nut06a_returns_tuple(self):
        """get_erfa_nutation_nut06a should return a tuple of two floats."""
        result = get_erfa_nutation_nut06a(J2000_JD)
        assert result is not None
        assert isinstance(result, tuple)
        assert len(result) == 2
        dpsi, deps = result
        assert isinstance(dpsi, float)
        assert isinstance(deps, float)

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_pnm06a_returns_matrix(self):
        """get_erfa_pnm06a_matrix should return a 3x3 matrix."""
        result = get_erfa_pnm06a_matrix(J2000_JD)
        assert result is not None
        # pyerfa returns numpy array
        assert hasattr(result, "shape")
        assert result.shape == (3, 3)

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_obliquity_iau2006_returns_float(self):
        """get_erfa_obliquity_iau2006 should return a float in radians."""
        result = get_erfa_obliquity_iau2006(J2000_JD)
        assert result is not None
        assert isinstance(result, float)
        # Obliquity at J2000 should be approximately 23.44 degrees
        obliquity_deg = math.degrees(result)
        assert 23.4 < obliquity_deg < 23.5

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_nutation_values_reasonable(self):
        """Nutation values should be in reasonable ranges."""
        result = get_erfa_nutation_nut06a(J2000_JD)
        dpsi, deps = result

        # Nutation in longitude is typically ±20 arcseconds
        dpsi_arcsec = math.degrees(dpsi) * 3600
        assert -25 < dpsi_arcsec < 25

        # Nutation in obliquity is typically ±10 arcseconds
        deps_arcsec = math.degrees(deps) * 3600
        assert -15 < deps_arcsec < 15


class TestCachedNutation:
    """Test cached nutation functionality."""

    def setup_method(self):
        """Clear cache before each test."""
        clear_erfa_nutation_cache()

    def test_cached_nutation_returns_tuple(self):
        """get_erfa_nutation_cached should always return a tuple."""
        # This should work even without ERFA (uses Skyfield fallback)
        result = get_erfa_nutation_cached(J2000_JD)
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_cache_works(self):
        """Caching should store results."""
        clear_erfa_nutation_cache()
        info_before = get_erfa_nutation_cache_info()

        # First call
        get_erfa_nutation_cached(J2000_JD)
        info_after_first = get_erfa_nutation_cache_info()

        # Second call with same value
        get_erfa_nutation_cached(J2000_JD)
        info_after_second = get_erfa_nutation_cache_info()

        # After first call, should have 1 miss
        assert info_after_first["misses"] == info_before["misses"] + 1

        # After second call, should have 1 hit
        assert info_after_second["hits"] == info_after_first["hits"] + 1

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_model_selection_nut00a(self):
        """Should be able to select nut00a model."""
        result = get_erfa_nutation_cached(J2000_JD, model="nut00a")
        assert isinstance(result, tuple)
        assert len(result) == 2

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_model_selection_nut06a(self):
        """Should be able to select nut06a model."""
        result = get_erfa_nutation_cached(J2000_JD, model="nut06a")
        assert isinstance(result, tuple)
        assert len(result) == 2


class TestPrecisionComparison:
    """
    Test precision comparison between Skyfield and ERFA implementations.

    These tests verify the precision differences between:
    - Skyfield's IAU 2000A (iau2000a_radians)
    - ERFA's nut00a (should be identical to Skyfield's IAU 2000A)
    - ERFA's nut06a (IAU 2000A with IAU 2006 adjustments - most accurate)
    """

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_skyfield_vs_erfa_nut00a_equivalence(self):
        """
        Skyfield's iau2000a and ERFA's nut00a should produce nearly identical results.

        Both implement the same IAU 2000A nutation model (MHB2000).
        Any differences should be due to floating point implementation details,
        not algorithmic differences.
        """
        from skyfield.nutationlib import iau2000a_radians

        from libephemeris.state import get_timescale

        ts = get_timescale()

        test_dates = [
            J2000_JD,  # J2000.0
            J2000_JD + 3652.5,  # ~10 years later
            J2000_JD - 3652.5,  # ~10 years earlier
            J2000_JD + 36525,  # 100 years later
            J2000_JD - 36525,  # 100 years earlier
        ]

        for jd_tt in test_dates:
            t = ts.tt_jd(jd_tt)
            dpsi_skyfield, deps_skyfield = iau2000a_radians(t)
            dpsi_erfa, deps_erfa = get_erfa_nutation_nut00a(jd_tt)

            # Difference in microarcseconds
            dpsi_diff_uas = abs(math.degrees(dpsi_skyfield - dpsi_erfa)) * 3600 * 1e6
            deps_diff_uas = abs(math.degrees(deps_skyfield - deps_erfa)) * 3600 * 1e6

            # Should agree to within 1 microarcsecond (numerical precision)
            assert dpsi_diff_uas < 1.0, (
                f"dpsi diff at JD {jd_tt}: {dpsi_diff_uas:.3f} µas"
            )
            assert deps_diff_uas < 1.0, (
                f"deps diff at JD {jd_tt}: {deps_diff_uas:.3f} µas"
            )

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_nut00a_vs_nut06a_differences(self):
        """
        ERFA's nut00a and nut06a should differ by small amounts.

        nut06a includes:
        - J2 secular variation correction
        - IAU 2006 obliquity frame correction

        These differences should be on the order of 0.01-0.1 milliarcseconds
        for dates within a century of J2000.
        """
        test_dates = [
            J2000_JD,  # J2000.0
            J2000_JD + 3652.5,  # ~10 years later
            J2000_JD + 36525,  # 100 years later
        ]

        for jd_tt in test_dates:
            dpsi_00a, deps_00a = get_erfa_nutation_nut00a(jd_tt)
            dpsi_06a, deps_06a = get_erfa_nutation_nut06a(jd_tt)

            # Difference in milliarcseconds
            dpsi_diff_mas = abs(math.degrees(dpsi_00a - dpsi_06a)) * 3600 * 1000
            deps_diff_mas = abs(math.degrees(deps_00a - deps_06a)) * 3600 * 1000

            # Should differ by less than 1 milliarcsecond for recent dates
            years_from_j2000 = (jd_tt - J2000_JD) / 365.25
            max_expected_mas = 0.1 + 0.005 * abs(years_from_j2000)  # Grows with time

            assert dpsi_diff_mas < max_expected_mas, (
                f"dpsi diff at JD {jd_tt} ({years_from_j2000:.0f} years from J2000): "
                f"{dpsi_diff_mas:.4f} mas (expected < {max_expected_mas:.4f} mas)"
            )
            assert deps_diff_mas < max_expected_mas, (
                f"deps diff at JD {jd_tt} ({years_from_j2000:.0f} years from J2000): "
                f"{deps_diff_mas:.4f} mas (expected < {max_expected_mas:.4f} mas)"
            )

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_iau2000b_vs_iau2000a_differences(self):
        """
        IAU 2000B and IAU 2000A should differ by ~1 milliarcsecond.

        IAU 2000B is a truncated version of IAU 2000A with only 77 terms
        (vs 1365 terms). The difference indicates the precision improvement
        from using the full model.
        """
        from skyfield.nutationlib import iau2000a_radians, iau2000b_radians

        from libephemeris.state import get_timescale

        ts = get_timescale()

        test_dates = [
            J2000_JD,
            J2000_JD + 365.25 * 10,  # 10 years
            J2000_JD + 365.25 * 50,  # 50 years
        ]

        for jd_tt in test_dates:
            t = ts.tt_jd(jd_tt)
            dpsi_2000a, deps_2000a = iau2000a_radians(t)
            dpsi_2000b, deps_2000b = iau2000b_radians(t)

            # Difference in milliarcseconds
            dpsi_diff_mas = abs(math.degrees(dpsi_2000a - dpsi_2000b)) * 3600 * 1000
            deps_diff_mas = abs(math.degrees(deps_2000a - deps_2000b)) * 3600 * 1000

            # IAU 2000B should be accurate to ~1 mas, so differences can be 0-1 mas
            assert dpsi_diff_mas < 2.0, (
                f"dpsi diff at JD {jd_tt}: {dpsi_diff_mas:.4f} mas"
            )
            assert deps_diff_mas < 2.0, (
                f"deps diff at JD {jd_tt}: {deps_diff_mas:.4f} mas"
            )


class TestCompareNutationModels:
    """Test the comparison function."""

    def test_compare_function_returns_dict(self):
        """compare_nutation_models should return a dictionary."""
        result = compare_nutation_models(J2000_JD)
        assert isinstance(result, dict)
        assert "jd_tt" in result
        assert result["jd_tt"] == J2000_JD

    def test_compare_function_has_skyfield_results(self):
        """Compare function should always include Skyfield results."""
        result = compare_nutation_models(J2000_JD)
        assert "skyfield_iau2000b" in result
        assert "skyfield_iau2000a" in result
        assert "dpsi_arcsec" in result["skyfield_iau2000b"]
        assert "deps_arcsec" in result["skyfield_iau2000b"]

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_compare_function_has_erfa_results(self):
        """Compare function should include ERFA results when available."""
        result = compare_nutation_models(J2000_JD)
        assert "erfa_nut00a" in result
        assert "erfa_nut06a" in result
        assert "differences_mas" in result


class TestPNM06aMatrix:
    """Test the precession-nutation matrix functionality."""

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_pnm06a_is_rotation_matrix(self):
        """pnm06a should return a proper rotation matrix."""
        import numpy as np

        rbpn = get_erfa_pnm06a_matrix(J2000_JD)

        # A rotation matrix should be orthogonal: R @ R.T = I
        identity = rbpn @ rbpn.T
        assert np.allclose(identity, np.eye(3), atol=1e-15)

        # Determinant should be +1
        det = np.linalg.det(rbpn)
        assert np.isclose(det, 1.0, atol=1e-15)

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_pnm06a_at_j2000_near_identity(self):
        """At J2000.0, the matrix should be close to identity (small corrections)."""

        rbpn = get_erfa_pnm06a_matrix(J2000_JD)

        # The off-diagonal elements represent small rotations
        # For J2000, they should be very small but non-zero (frame bias)
        off_diag = [
            rbpn[0, 1],
            rbpn[0, 2],
            rbpn[1, 0],
            rbpn[1, 2],
            rbpn[2, 0],
            rbpn[2, 1],
        ]

        # At J2000, the matrix includes:
        # - Frame bias (~17 mas in x, ~7 mas in y, ~80 mas in z)
        # - Nutation (~10-20 arcseconds, i.e., ~1e-4 radians)
        # Off-diagonal elements represent rotations, so they should be < 1e-4
        for val in off_diag:
            assert abs(val) < 1e-3, f"Off-diagonal element {val} too large for J2000"

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_pnm06a_transforms_icrs_to_true_equator(self):
        """
        pnm06a should correctly transform ICRS coordinates to true equator of date.

        The celestial pole (ICRS z-axis) should be transformed to the
        true pole position (with precession and nutation).
        """
        import numpy as np

        # 50 years from J2000 to see significant precession
        jd = J2000_JD + 50 * 365.25
        rbpn = get_erfa_pnm06a_matrix(jd)

        # ICRS north pole
        icrs_pole = np.array([0.0, 0.0, 1.0])

        # Transform to true equator of date
        true_pole = rbpn @ icrs_pole

        # The pole should have moved slightly due to precession
        # Precession is ~50.3" per year in longitude, so after 50 years
        # the pole has moved significantly
        assert not np.allclose(true_pole, icrs_pole, atol=1e-6)

        # But magnitude should still be 1
        assert np.isclose(np.linalg.norm(true_pole), 1.0, atol=1e-15)


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_distant_past_date(self):
        """Functions should work for dates in the distant past."""
        # 2000 years ago
        jd_tt = J2000_JD - 2000 * 365.25
        result = get_erfa_nutation_nut06a(jd_tt)
        assert result is not None
        dpsi, deps = result
        assert math.isfinite(dpsi)
        assert math.isfinite(deps)

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_distant_future_date(self):
        """Functions should work for dates in the distant future."""
        # 2000 years in future
        jd_tt = J2000_JD + 2000 * 365.25
        result = get_erfa_nutation_nut06a(jd_tt)
        assert result is not None
        dpsi, deps = result
        assert math.isfinite(dpsi)
        assert math.isfinite(deps)

    def test_invalid_model_raises_error(self):
        """Invalid model name should raise ValueError."""
        clear_erfa_nutation_cache()
        # This depends on pyerfa being available to reach the error
        if HAS_ERFA:
            with pytest.raises(ValueError, match="Unknown nutation model"):
                get_erfa_nutation_cached(J2000_JD, model="invalid_model")


class TestPrecisionImprovementEvaluation:
    """
    Evaluation tests to quantify precision improvements from ERFA integration.

    These tests document and verify the expected precision improvements.
    """

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_precision_improvement_summary(self):
        """
        Document the precision hierarchy of nutation models.

        Expected precision (in milliarcseconds):
        - IAU 2000B: ~1 mas
        - IAU 2000A: ~0.1 mas (10x improvement)
        - IAU 2006/2000A (nut06a): ~0.01-0.05 mas (2-10x further improvement)
        """
        from skyfield.nutationlib import iau2000a_radians, iau2000b_radians

        from libephemeris.state import get_timescale

        ts = get_timescale()
        jd_tt = J2000_JD + 10 * 365.25  # 10 years from J2000

        t = ts.tt_jd(jd_tt)

        # Get all nutation values
        dpsi_2000b, deps_2000b = iau2000b_radians(t)
        dpsi_2000a, deps_2000a = iau2000a_radians(t)
        dpsi_06a, deps_06a = get_erfa_nutation_nut06a(jd_tt)

        # Convert to milliarcseconds for easier interpretation
        def to_mas(rad):
            return math.degrees(rad) * 3600 * 1000

        # Compute differences
        diff_2000b_to_2000a_dpsi = abs(to_mas(dpsi_2000b) - to_mas(dpsi_2000a))
        diff_2000b_to_2000a_deps = abs(to_mas(deps_2000b) - to_mas(deps_2000a))

        diff_2000a_to_06a_dpsi = abs(to_mas(dpsi_2000a) - to_mas(dpsi_06a))
        diff_2000a_to_06a_deps = abs(to_mas(deps_2000a) - to_mas(deps_06a))

        # Document the precision hierarchy
        print(f"\n=== Nutation Precision Evaluation at JD {jd_tt} ===")
        print("\nIAU 2000B (77 terms, ~1 mas precision):")
        print(f"  dpsi = {to_mas(dpsi_2000b):.4f} mas")
        print(f"  deps = {to_mas(deps_2000b):.4f} mas")

        print("\nIAU 2000A (1365 terms, ~0.1 mas precision):")
        print(f"  dpsi = {to_mas(dpsi_2000a):.4f} mas")
        print(f"  deps = {to_mas(deps_2000a):.4f} mas")
        print(f"  Improvement over 2000B: {diff_2000b_to_2000a_dpsi:.4f} mas (dpsi)")
        print(f"                         {diff_2000b_to_2000a_deps:.4f} mas (deps)")

        print("\nIAU 2006/2000A (nut06a, ~0.01-0.05 mas precision):")
        print(f"  dpsi = {to_mas(dpsi_06a):.4f} mas")
        print(f"  deps = {to_mas(deps_06a):.4f} mas")
        print(f"  Improvement over 2000A: {diff_2000a_to_06a_dpsi:.4f} mas (dpsi)")
        print(f"                         {diff_2000a_to_06a_deps:.4f} mas (deps)")

        # Verify expected precision levels
        # 2000B vs 2000A should differ by up to ~1 mas
        assert diff_2000b_to_2000a_dpsi < 2.0
        assert diff_2000b_to_2000a_deps < 2.0

        # 2000A vs 06a should differ by up to ~0.1 mas for dates near J2000
        assert diff_2000a_to_06a_dpsi < 0.5
        assert diff_2000a_to_06a_deps < 0.5

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_precision_growth_over_time(self):
        """
        Evaluate how nutation model differences grow over time from J2000.

        This test documents the error growth pattern to help users understand
        when pyerfa's higher precision models become important.
        """
        from skyfield.nutationlib import iau2000a_radians, iau2000b_radians

        from libephemeris.state import get_timescale

        ts = get_timescale()

        # Test at different time offsets from J2000
        time_offsets_years = [0, 10, 50, 100, 200, 500]

        def to_mas(rad):
            return math.degrees(rad) * 3600 * 1000

        print("\n=== Nutation Precision Growth Analysis ===")
        print(
            f"{'Years':>6} | {'2000B-2000A dpsi':>16} | {'2000A-06a dpsi':>14} | "
            f"{'2000B-2000A deps':>16} | {'2000A-06a deps':>14}"
        )
        print("-" * 85)

        for years in time_offsets_years:
            jd_tt = J2000_JD + years * 365.25
            t = ts.tt_jd(jd_tt)

            dpsi_2000b, deps_2000b = iau2000b_radians(t)
            dpsi_2000a, deps_2000a = iau2000a_radians(t)
            dpsi_06a, deps_06a = get_erfa_nutation_nut06a(jd_tt)

            diff_b_a_dpsi = abs(to_mas(dpsi_2000b) - to_mas(dpsi_2000a))
            diff_a_06_dpsi = abs(to_mas(dpsi_2000a) - to_mas(dpsi_06a))
            diff_b_a_deps = abs(to_mas(deps_2000b) - to_mas(deps_2000a))
            diff_a_06_deps = abs(to_mas(deps_2000a) - to_mas(deps_06a))

            print(
                f"{years:>6} | {diff_b_a_dpsi:>16.4f} | {diff_a_06_dpsi:>14.4f} | "
                f"{diff_b_a_deps:>16.4f} | {diff_a_06_deps:>14.4f}"
            )

            # Verify differences stay within expected bounds
            # IAU 2000B vs 2000A difference grows with time from J2000
            # Near J2000: < 1 mas; At 100 years: < 2 mas; At 500 years: can reach 20+ mas
            # This reflects that IAU 2000B is a truncated model optimized for near-J2000
            max_b_a_diff = 2.5 + 0.03 * years  # ~2.5 mas at J2000, grows with time
            assert diff_b_a_dpsi < max_b_a_diff, (
                f"2000B-2000A dpsi too large at {years} years"
            )
            assert diff_b_a_deps < max_b_a_diff, (
                f"2000B-2000A deps too large at {years} years"
            )

            # IAU 2006 correction grows slowly with time from J2000
            # Expected: ~0.001-0.01 mas near J2000, up to ~1 mas at 500 years
            max_06_diff = 0.1 + 0.002 * years
            assert diff_a_06_dpsi < max_06_diff, (
                f"2000A-06a dpsi too large at {years} years"
            )
            assert diff_a_06_deps < max_06_diff, (
                f"2000A-06a deps too large at {years} years"
            )

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_pnm06a_vs_separate_rotations(self):
        """
        Compare pnm06a combined matrix vs applying bias, precession, nutation separately.

        The combined matrix avoids cross-term errors that accumulate when rotations
        are applied separately. This test quantifies that improvement.
        """
        import numpy as np

        # Test at a distant date where differences are more pronounced
        jd_tt = J2000_JD + 100 * 365.25  # 100 years from J2000

        # Get the combined pnm06a matrix
        rbpn_combined = get_erfa_pnm06a_matrix(jd_tt)

        # Build separate matrices using ERFA
        # Frame bias matrix
        rb = erfa.bi00()  # Returns dx, dy, dpsi (frame bias parameters)

        # For proper comparison, we need individual matrices
        # Get precession angles - erfa.p06e returns 16 values
        result = erfa.p06e(2451545.0, jd_tt - 2451545.0)
        # result contains: eps0, psia, oma, bpa, bqa, pia, bpia, epsa, chia, za, zetaa, thetaa, pa, gam, phi, psi
        # We need epsa (index 7) and the bias-precession matrix from erfa.bp06
        epsa = result[7]

        # Get the bias-precession matrix directly
        rb_matrix, rp, rbp = erfa.bp06(2451545.0, jd_tt - 2451545.0)

        # Get nutation matrix
        dpsi, deps = erfa.nut06a(2451545.0, jd_tt - 2451545.0)
        epsa_new = erfa.obl06(2451545.0, jd_tt - 2451545.0)
        rn = erfa.numat(epsa_new, dpsi, deps)

        # Combined separately: R = Rn @ Rp @ Rb
        # Note: rbp already contains bias-precession
        rbpn_separate = rn @ rbp

        # Compare the matrices
        diff = np.abs(rbpn_combined - rbpn_separate)
        max_diff = np.max(diff)

        print("\n=== pnm06a vs Separate Rotations Comparison ===")
        print(f"Date: J2000 + 100 years (JD {jd_tt})")
        print(f"Maximum matrix element difference: {max_diff:.2e}")

        # The difference should be very small (numerical precision)
        # but demonstrates the consistency of pnm06a
        assert max_diff < 1e-12, f"Matrix difference too large: {max_diff}"

        # Check orthogonality of both matrices
        for name, matrix in [
            ("pnm06a", rbpn_combined),
            ("separate", rbpn_separate),
        ]:
            identity = matrix @ matrix.T
            orth_error = np.max(np.abs(identity - np.eye(3)))
            assert orth_error < 1e-14, f"{name} orthogonality error: {orth_error}"

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_celestial_pole_position_precision(self):
        """
        Evaluate precision improvement by computing celestial pole position.

        The CIP (Celestial Intermediate Pole) position is a practical measure
        of nutation accuracy. Different models should give slightly different
        pole positions, demonstrating real-world precision differences.
        """
        import numpy as np

        from skyfield.nutationlib import iau2000a_radians

        from libephemeris.state import get_timescale

        ts = get_timescale()
        jd_tt = J2000_JD + 20 * 365.25  # 20 years from J2000

        # ICRS pole (z-axis)
        pole_icrs = np.array([0.0, 0.0, 1.0])

        # Using pnm06a (most accurate)
        rbpn_06a = get_erfa_pnm06a_matrix(jd_tt)
        pole_06a = rbpn_06a @ pole_icrs

        # Using Skyfield's IAU 2000A nutation
        # Build a simplified rotation using Skyfield nutation
        t = ts.tt_jd(jd_tt)
        dpsi_2000a, deps_2000a = iau2000a_radians(t)

        # Get obliquity for rotation
        eps = get_erfa_obliquity_iau2006(jd_tt)

        # Build nutation rotation (simplified, without precession for comparison)
        dpsi_06a, deps_06a = get_erfa_nutation_nut06a(jd_tt)

        # Convert nutation angles to pole offset (in arcseconds)
        # X = dpsi * sin(eps), Y = -deps
        x_2000a_mas = math.degrees(dpsi_2000a * math.sin(eps)) * 3600 * 1000
        y_2000a_mas = math.degrees(-deps_2000a) * 3600 * 1000

        x_06a_mas = math.degrees(dpsi_06a * math.sin(eps)) * 3600 * 1000
        y_06a_mas = math.degrees(-deps_06a) * 3600 * 1000

        # Difference in pole position
        dx_mas = abs(x_2000a_mas - x_06a_mas)
        dy_mas = abs(y_2000a_mas - y_06a_mas)
        pole_diff_mas = math.sqrt(dx_mas**2 + dy_mas**2)

        print("\n=== Celestial Pole Position Precision ===")
        print(f"Date: J2000 + 20 years (JD {jd_tt})")
        print(
            f"\nIAU 2000A pole offset (X, Y): ({x_2000a_mas:.4f}, {y_2000a_mas:.4f}) mas"
        )
        print(f"IAU 2006/2000A pole offset: ({x_06a_mas:.4f}, {y_06a_mas:.4f}) mas")
        print(f"\nPole position difference: {pole_diff_mas:.4f} mas")
        print(f"  dX: {dx_mas:.4f} mas")
        print(f"  dY: {dy_mas:.4f} mas")

        # The pole difference should be small (sub-milliarcsecond)
        assert pole_diff_mas < 0.1, f"Pole difference too large: {pole_diff_mas} mas"

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_astrological_precision_context(self):
        """
        Evaluate pyerfa precision in the context of astrological calculations.

        For astrology, position errors translate to timing errors for aspects
        and house cusps. This test quantifies the practical impact.

        Key insight: 1 arcsecond of error in Moon's position represents
        approximately 2 seconds of time, since the Moon moves ~0.5"/second.
        """
        # Typical astrological precision requirements:
        # - Natal chart: 1 arcminute (60") acceptable for most purposes
        # - Precise event timing: 1 arcsecond (1") desirable
        # - Research/high precision: 0.1 arcsecond (100 mas) desirable

        from skyfield.nutationlib import iau2000a_radians, iau2000b_radians

        from libephemeris.state import get_timescale

        ts = get_timescale()

        # Test current epoch (typical astrological use case)
        jd_tt = 2460000.5  # Early 2023

        t = ts.tt_jd(jd_tt)
        dpsi_2000b, deps_2000b = iau2000b_radians(t)
        dpsi_2000a, deps_2000a = iau2000a_radians(t)
        dpsi_06a, deps_06a = get_erfa_nutation_nut06a(jd_tt)

        def to_arcsec(rad):
            return math.degrees(rad) * 3600

        # Calculate impact on zodiacal longitude
        # Nutation in longitude directly affects zodiacal positions
        dpsi_2000b_as = to_arcsec(dpsi_2000b)
        dpsi_2000a_as = to_arcsec(dpsi_2000a)
        dpsi_06a_as = to_arcsec(dpsi_06a)

        diff_b_a = abs(dpsi_2000b_as - dpsi_2000a_as)
        diff_a_06 = abs(dpsi_2000a_as - dpsi_06a_as)

        print("\n=== Astrological Precision Context ===")
        print(f"Date: JD {jd_tt} (early 2023)")
        print("\nNutation in longitude (dpsi):")
        print(f'  IAU 2000B: {dpsi_2000b_as:.4f}"')
        print(f'  IAU 2000A: {dpsi_2000a_as:.4f}"')
        print(f'  IAU 2006/2000A: {dpsi_06a_as:.4f}"')
        print("\nZodiacal longitude differences:")
        print(f'  2000B vs 2000A: {diff_b_a * 1000:.4f} mas ({diff_b_a:.4f}")')
        print(f'  2000A vs 2006/2000A: {diff_a_06 * 1000:.4f} mas ({diff_a_06:.4f}")')

        # Moon timing impact (Moon moves ~0.5"/second)
        moon_timing_b_a = diff_b_a / 0.5  # seconds
        moon_timing_a_06 = diff_a_06 / 0.5  # seconds
        print('\nMoon timing impact (Moon speed ~0.5"/s):')
        print(f"  2000B vs 2000A: {moon_timing_b_a:.3f} seconds")
        print(f"  2000A vs 2006/2000A: {moon_timing_a_06:.3f} seconds")

        # Practical conclusion
        print("\n=== Practical Recommendations ===")
        if diff_a_06 < 0.001:  # < 1 mas = 0.001 arcseconds
            print("IAU 2006/2000A provides sub-milliarcsecond improvement over 2000A.")
            print("For typical astrological work, IAU 2000A is already sufficient.")
            print("pyerfa's nut06a is recommended for:")
            print("  - Research requiring highest precision")
            print("  - Long-term ephemeris calculations")
            print("  - Validation against official IERS data")

        # Verify precision levels
        assert diff_b_a < 1.0, "2000B should be accurate to ~1 arcsecond"
        assert diff_a_06 < 0.01, "nut06a should improve by < 10 mas over 2000A"

    @pytest.mark.skipif(not HAS_ERFA, reason="pyerfa not installed")
    def test_obliquity_model_comparison(self):
        """
        Compare obliquity values from different models.

        The mean obliquity of the ecliptic affects the transformation between
        ecliptic and equatorial coordinates. Different models give slightly
        different values.
        """
        # Test at various dates
        test_dates = [
            (J2000_JD, "J2000.0"),
            (J2000_JD + 36525, "J2100.0"),
            (J2000_JD - 36525, "J1900.0"),
        ]

        print("\n=== Obliquity Model Comparison ===")

        for jd_tt, date_name in test_dates:
            # IAU 2006 obliquity from ERFA
            eps_iau2006 = get_erfa_obliquity_iau2006(jd_tt)

            # Laskar 1986 formula (used as fallback in libephemeris)
            T = (jd_tt - J2000_JD) / 36525.0
            eps_laskar = math.radians(
                23.439291111
                - 0.013004166667 * T
                - 1.638888889e-7 * T**2
                + 5.036111111e-7 * T**3
            )

            # Difference in arcseconds
            diff_as = abs(math.degrees(eps_iau2006 - eps_laskar)) * 3600

            print(f"\n{date_name} (JD {jd_tt}):")
            print(f"  IAU 2006: {math.degrees(eps_iau2006):.8f}°")
            print(f"  Laskar 1986: {math.degrees(eps_laskar):.8f}°")
            print(f'  Difference: {diff_as * 1000:.4f} mas ({diff_as:.6f}")')

            # The difference should be reasonable for dates near J2000
            # Note: The Laskar 1986 vs IAU 2006 difference at J2000 is ~42 mas
            # which is a known offset between the two models
            # The difference grows roughly linearly with time: ~0.2 mas per year from J2000
            years_from_j2000 = abs((jd_tt - J2000_JD) / 365.25)
            max_expected_mas = 50 + 0.25 * years_from_j2000  # Baseline ~42 mas + growth
            assert diff_as * 1000 < max_expected_mas, (
                f"Obliquity difference too large at {date_name}"
            )
