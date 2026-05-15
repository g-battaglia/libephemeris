"""
Tests for the LEB reader module (leb_reader.py).

Tests Clenshaw algorithm, eval_body, eval_nutation, delta_t, and edge cases.
"""

from __future__ import annotations

import numpy as np
import pytest
from numpy.polynomial.chebyshev import chebder, chebval

from libephemeris.constants import EARTH, MARS, MEAN_NODE, MOON, SUN
from libephemeris.leb_reader import (
    LEBReader,
    _clenshaw,
    _clenshaw_derivative,
    _clenshaw_with_derivative,
    _deriv_coeffs,
)


# =============================================================================
# CLENSHAW ALGORITHM TESTS
# =============================================================================


class TestClenshaw:
    """Test the pure-Python Clenshaw evaluation."""

    @pytest.mark.unit
    def test_constant_polynomial(self):
        """Chebyshev [5.0] should evaluate to 5.0 everywhere."""
        for tau in [-1.0, -0.5, 0.0, 0.5, 1.0]:
            assert _clenshaw((5.0,), tau) == 5.0

    @pytest.mark.unit
    def test_linear_polynomial(self):
        """Chebyshev [a, b] = a + b*x."""
        # T_0(x) = 1, T_1(x) = x, so [3, 2] -> 3 + 2x
        assert abs(_clenshaw((3.0, 2.0), 0.5) - 4.0) < 1e-14
        assert abs(_clenshaw((3.0, 2.0), -1.0) - 1.0) < 1e-14
        assert abs(_clenshaw((3.0, 2.0), 1.0) - 5.0) < 1e-14

    @pytest.mark.unit
    def test_quadratic_polynomial(self):
        """Chebyshev [a, b, c] = a + b*T_1(x) + c*T_2(x), T_2(x)=2x^2-1."""
        coeffs = (1.0, 0.5, -0.3)
        for tau in [-1.0, -0.5, 0.0, 0.5, 1.0]:
            expected = float(chebval(tau, list(coeffs)))
            actual = _clenshaw(coeffs, tau)
            assert abs(actual - expected) < 1e-14, f"tau={tau}: {actual} != {expected}"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "coeffs",
        [
            (1.0, 0.5, -0.3, 0.1),
            (1.0, 0.5, -0.3, 0.1, -0.05, 0.02),
            tuple(np.random.default_rng(42).standard_normal(14) * 0.1),
        ],
        ids=["degree3", "degree5", "degree13"],
    )
    def test_matches_numpy(self, coeffs):
        """Clenshaw matches numpy.polynomial.chebyshev.chebval."""
        for tau in np.linspace(-1, 1, 21):
            expected = float(chebval(tau, list(coeffs)))
            actual = _clenshaw(coeffs, float(tau))
            assert abs(actual - expected) < 1e-12, f"tau={tau}"


class TestClenshawDerivative:
    """Test the derivative Chebyshev evaluation."""

    @pytest.mark.unit
    def test_constant_derivative_is_zero(self):
        """Derivative of a constant is 0."""
        assert _clenshaw_derivative((5.0,), 0.0) == 0.0
        assert _clenshaw_derivative((5.0,), 0.7) == 0.0

    @pytest.mark.unit
    def test_linear_derivative(self):
        """Derivative of [a, b] = b (coefficient of T_1)."""
        assert abs(_clenshaw_derivative((3.0, 2.0), 0.0) - 2.0) < 1e-14
        assert abs(_clenshaw_derivative((3.0, 2.0), 0.5) - 2.0) < 1e-14

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "coeffs",
        [
            (1.0, 0.5, -0.3),
            (1.0, 0.5, -0.3, 0.1),
            (1.0, 0.5, -0.3, 0.1, -0.05, 0.02),
        ],
        ids=["degree2", "degree3", "degree5"],
    )
    def test_derivative_matches_numpy(self, coeffs):
        """Derivative matches numpy chebder + chebval."""
        d_coeffs = chebder(list(coeffs))
        for tau in np.linspace(-1, 1, 11):
            expected = float(chebval(tau, d_coeffs))
            actual = _clenshaw_derivative(coeffs, float(tau))
            assert abs(actual - expected) < 1e-12, f"tau={tau}"

    @pytest.mark.unit
    def test_degree13_derivative(self):
        """Degree 13 (typical planet) derivative accuracy."""
        rng = np.random.default_rng(42)
        coeffs = tuple(rng.standard_normal(14) * 0.1)
        d_coeffs = chebder(list(coeffs))
        max_err = 0.0
        for tau in np.linspace(-1, 1, 21):
            expected = float(chebval(tau, d_coeffs))
            actual = _clenshaw_derivative(coeffs, float(tau))
            max_err = max(max_err, abs(actual - expected))
        assert max_err < 1e-10, f"max derivative error = {max_err}"


class TestClenshawWithDerivative:
    """Test simultaneous value + derivative evaluation."""

    @pytest.mark.unit
    def test_constant(self):
        """Constant: value=c, derivative=0."""
        val, deriv = _clenshaw_with_derivative((5.0,), 0.0)
        assert val == 5.0
        assert deriv == 0.0

    @pytest.mark.unit
    def test_linear(self):
        """Linear: [3, 2] -> val=4.0 at x=0.5, deriv=2.0."""
        val, deriv = _clenshaw_with_derivative((3.0, 2.0), 0.5)
        assert abs(val - 4.0) < 1e-14
        assert abs(deriv - 2.0) < 1e-14

    @pytest.mark.unit
    def test_both_match_separate_calls(self):
        """Combined function matches separate value and derivative calls."""
        coeffs = (1.0, 0.5, -0.3, 0.1, -0.05, 0.02)
        for tau in np.linspace(-1, 1, 11):
            val, deriv = _clenshaw_with_derivative(coeffs, float(tau))
            val_sep = _clenshaw(coeffs, float(tau))
            deriv_sep = _clenshaw_derivative(coeffs, float(tau))
            assert abs(val - val_sep) < 1e-15
            assert abs(deriv - deriv_sep) < 1e-13

    @pytest.mark.unit
    def test_known_chebyshev(self):
        """Verify against numpy for a 4-term series."""
        coeffs = (1.0, 0.5, -0.3, 0.1)
        d_coeffs = chebder(list(coeffs))
        for tau in [-1.0, -0.5, 0.0, 0.5, 1.0]:
            expected_val = float(chebval(tau, list(coeffs)))
            expected_deriv = float(chebval(tau, d_coeffs))
            val, deriv = _clenshaw_with_derivative(coeffs, tau)
            assert abs(val - expected_val) < 1e-13
            assert abs(deriv - expected_deriv) < 1e-13


class TestDerivCoeffs:
    """Test the derivative coefficient computation."""

    @pytest.mark.unit
    def test_constant_gives_zero(self):
        """Derivative of degree-0 is (0.0,)."""
        assert _deriv_coeffs((5.0,)) == (0.0,)

    @pytest.mark.unit
    def test_linear(self):
        """Derivative of [a, b] is (b,)."""
        d = _deriv_coeffs((3.0, 2.0))
        assert len(d) == 1
        assert abs(d[0] - 2.0) < 1e-14

    @pytest.mark.unit
    def test_matches_numpy_chebder(self):
        """Derivative coefficients match numpy chebder for degree 5."""
        coeffs = (1.0, 0.5, -0.3, 0.1, -0.05, 0.02)
        d = _deriv_coeffs(coeffs)
        np_d = chebder(list(coeffs))
        assert len(d) == len(np_d)
        for i in range(len(d)):
            assert abs(d[i] - np_d[i]) < 1e-12, f"index {i}: {d[i]} != {np_d[i]}"


# =============================================================================
# LEBREADER TESTS
# =============================================================================


class TestLEBReaderInit:
    """Test LEBReader construction and file parsing."""

    @pytest.mark.unit
    def test_file_not_found(self, tmp_path):
        """LEBReader raises FileNotFoundError for missing files."""
        with pytest.raises(FileNotFoundError):
            LEBReader(str(tmp_path / "nonexistent.leb"))

    @pytest.mark.unit
    def test_invalid_magic(self, tmp_path):
        """LEBReader raises ValueError for wrong magic bytes."""
        bad_file = tmp_path / "bad.leb"
        bad_file.write_bytes(b"XXXX" + b"\x00" * 200)
        with pytest.raises(ValueError, match="Invalid LEB magic"):
            LEBReader(str(bad_file))

    @pytest.mark.integration
    def test_opens_valid_file(self, leb_reader):
        """LEBReader can open a valid .leb file."""
        assert leb_reader is not None
        jd_start, jd_end = leb_reader.jd_range
        assert jd_start < jd_end

    @pytest.mark.integration
    def test_has_expected_bodies(self, leb_reader):
        """Test .leb file contains the bodies from conftest."""
        assert leb_reader.has_body(SUN)
        assert leb_reader.has_body(MOON)
        assert leb_reader.has_body(MARS)
        assert leb_reader.has_body(EARTH)
        assert leb_reader.has_body(MEAN_NODE)

    @pytest.mark.integration
    def test_does_not_have_unknown_body(self, leb_reader):
        """has_body returns False for bodies not in the file."""
        assert not leb_reader.has_body(99999)

    @pytest.mark.integration
    def test_context_manager(self, test_leb_file):
        """LEBReader works as a context manager."""
        with LEBReader(test_leb_file) as reader:
            assert reader.has_body(SUN)


class TestEvalBody:
    """Test eval_body for ICRS and ecliptic bodies."""

    @pytest.mark.integration
    def test_sun_returns_3_components(self, leb_reader):
        """Sun position returns 3 position and 3 velocity components."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0
        pos, vel = leb_reader.eval_body(SUN, jd_mid)
        assert len(pos) == 3
        assert len(vel) == 3

    @pytest.mark.integration
    def test_sun_position_reasonable(self, leb_reader):
        """Sun ICRS barycentric position should have reasonable magnitude."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0
        pos, vel = leb_reader.eval_body(SUN, jd_mid)
        # ICRS barycentric: Sun is near SSB, typically <0.01 AU from it
        import math

        dist = math.sqrt(pos[0] ** 2 + pos[1] ** 2 + pos[2] ** 2)
        assert dist < 0.02, f"Sun barycentric distance = {dist} AU (too large)"

    @pytest.mark.integration
    def test_mars_position_reasonable(self, leb_reader):
        """Mars ICRS barycentric position should have reasonable magnitude."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0
        pos, vel = leb_reader.eval_body(MARS, jd_mid)
        # ICRS barycentric: Mars orbits at ~1.38-1.67 AU from Sun
        import math

        dist = math.sqrt(pos[0] ** 2 + pos[1] ** 2 + pos[2] ** 2)
        assert 1.0 < dist < 2.0, f"Mars barycentric distance = {dist} AU"

    @pytest.mark.integration
    def test_mean_node_longitude_range(self, leb_reader):
        """Mean node longitude should be in [0, 360)."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0
        pos, vel = leb_reader.eval_body(MEAN_NODE, jd_mid)
        assert 0.0 <= pos[0] < 360.0, f"Mean node lon = {pos[0]}"

    @pytest.mark.integration
    def test_mean_node_velocity_negative(self, leb_reader):
        """Mean node longitude should be retrograde (negative speed)."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0
        pos, vel = leb_reader.eval_body(MEAN_NODE, jd_mid)
        # Mean node regresses at ~-0.053 deg/day
        assert vel[0] < 0, f"Mean node velocity = {vel[0]} (expected negative)"

    @pytest.mark.integration
    def test_sun_matches_skyfield(self, leb_reader):
        """Sun ICRS barycentric from .leb matches Skyfield within tolerance.

        Note: This test uses a small on-the-fly LEB file with default segment
        parameters.  The tolerance is generous to accommodate Chebyshev
        fitting error.  Production accuracy is validated by the compare/ suite.
        """
        from libephemeris.constants import FLG_SPEED
        from libephemeris.fast_calc import fast_calc_ut

        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0

        # Full pipeline: LEB → fast_calc → geocentric ecliptic of date
        fast_result, _ = fast_calc_ut(leb_reader, jd_mid, SUN, FLG_SPEED)

        # Skyfield reference via calc_ut (geocentric ecliptic of date)
        import libephemeris as ephem

        ref, _ = ephem.calc_ut(jd_mid, SUN, FLG_SPEED)

        # Longitude
        lon_err = abs(fast_result[0] - ref[0])
        if lon_err > 180.0:
            lon_err = 360.0 - lon_err
        lon_err_arcsec = lon_err * 3600.0
        assert lon_err_arcsec < 5.0, (
            f"Sun longitude error = {lon_err_arcsec:.4f} arcsec"
        )

        # Latitude
        lat_err_arcsec = abs(fast_result[1] - ref[1]) * 3600.0
        assert lat_err_arcsec < 5.0, f"Sun latitude error = {lat_err_arcsec:.4f} arcsec"

        # Distance
        dist_err = abs(fast_result[2] - ref[2])
        assert dist_err < 1e-4, f"Sun distance error = {dist_err:.2e} AU"

    @pytest.mark.integration
    def test_body_not_found(self, leb_reader):
        """eval_body raises KeyError for unknown body."""
        jd_start, _ = leb_reader.jd_range
        with pytest.raises(KeyError, match="Body 99999"):
            leb_reader.eval_body(99999, jd_start + 1)

    @pytest.mark.integration
    def test_jd_out_of_range(self, leb_reader):
        """eval_body raises ValueError for JD outside range."""
        jd_start, jd_end = leb_reader.jd_range
        with pytest.raises(ValueError, match="outside range"):
            leb_reader.eval_body(SUN, jd_start - 100.0)
        with pytest.raises(ValueError, match="outside range"):
            leb_reader.eval_body(SUN, jd_end + 100.0)

    @pytest.mark.integration
    def test_eval_at_segment_boundaries(self, leb_reader):
        """Evaluation at segment boundaries should not crash."""
        jd_start, jd_end = leb_reader.jd_range
        body = leb_reader._bodies[SUN]
        interval = body.interval_days

        # Evaluate at exact segment boundaries
        for i in range(min(5, body.segment_count)):
            jd = jd_start + i * interval
            pos, vel = leb_reader.eval_body(SUN, jd)
            assert len(pos) == 3

    @pytest.mark.integration
    def test_velocity_consistent_with_finite_diff(self, leb_reader):
        """Analytical velocity should be close to finite-difference velocity."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0
        dt = 0.01  # 0.01 day

        pos0, vel0 = leb_reader.eval_body(SUN, jd_mid)
        pos_prev, _ = leb_reader.eval_body(SUN, jd_mid - dt)
        pos_next, _ = leb_reader.eval_body(SUN, jd_mid + dt)

        for c in range(3):
            fd_vel = (pos_next[c] - pos_prev[c]) / (2.0 * dt)
            rel_err = abs(vel0[c] - fd_vel) / (abs(vel0[c]) + 1e-20)
            assert rel_err < 0.01, (
                f"Component {c}: analytical={vel0[c]:.8e}, "
                f"finite-diff={fd_vel:.8e}, rel_err={rel_err:.2e}"
            )


class TestEvalNutation:
    """Test nutation evaluation."""

    @pytest.mark.integration
    def test_nutation_returns_two_values(self, leb_reader):
        """eval_nutation returns (dpsi, deps) tuple."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0
        dpsi, deps = leb_reader.eval_nutation(jd_mid)
        assert isinstance(dpsi, float)
        assert isinstance(deps, float)

    @pytest.mark.integration
    def test_nutation_reasonable_values(self, leb_reader):
        """Nutation values should be small (< 0.001 rad ≈ 0.057 deg)."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0
        dpsi, deps = leb_reader.eval_nutation(jd_mid)
        # Nutation is typically < 20 arcsec ≈ 0.0001 rad
        assert abs(dpsi) < 0.001, f"dpsi = {dpsi} rad (too large)"
        assert abs(deps) < 0.001, f"deps = {deps} rad (too large)"

    @pytest.mark.integration
    def test_nutation_matches_erfa(self, leb_reader):
        """Nutation from .leb should match erfa within ~0.001 arcsec."""
        import erfa

        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0

        dpsi_leb, deps_leb = leb_reader.eval_nutation(jd_mid)
        dpsi_erfa, deps_erfa = erfa.nut06a(2451545.0, jd_mid - 2451545.0)

        # Tolerance: 0.01 arcsec = 4.8e-8 rad
        tol_rad = 0.01 / 206265.0
        assert abs(dpsi_leb - dpsi_erfa) < tol_rad, (
            f"dpsi error = {abs(dpsi_leb - dpsi_erfa) * 206265:.4f} arcsec"
        )
        assert abs(deps_leb - deps_erfa) < tol_rad, (
            f"deps error = {abs(deps_leb - deps_erfa) * 206265:.4f} arcsec"
        )


class TestDeltaT:
    """Test Delta-T interpolation."""

    @pytest.mark.integration
    def test_delta_t_returns_float(self, leb_reader):
        """delta_t returns a float."""
        jd_start, _ = leb_reader.jd_range
        dt = leb_reader.delta_t(jd_start + 1)
        assert isinstance(dt, float)

    @pytest.mark.integration
    def test_delta_t_reasonable(self, leb_reader):
        """Delta-T should be positive and reasonable (< 0.002 days ≈ 170 sec)."""
        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0
        dt = leb_reader.delta_t(jd_mid)
        # In 2020s, Delta-T ≈ 69 seconds ≈ 0.0008 days
        assert 0.0 < dt < 0.002, f"Delta-T = {dt} days ({dt * 86400:.1f} sec)"

    @pytest.mark.integration
    def test_delta_t_matches_swe(self, leb_reader):
        """Delta-T from .leb should match deltat within 0.1 seconds."""
        from libephemeris.time_utils import deltat

        jd_start, jd_end = leb_reader.jd_range
        jd_mid = (jd_start + jd_end) / 2.0
        dt_leb = leb_reader.delta_t(jd_mid)
        dt_swe = deltat(jd_mid)

        err_sec = abs(dt_leb - dt_swe) * 86400.0
        assert err_sec < 0.1, f"Delta-T error = {err_sec:.4f} sec"

    @pytest.mark.integration
    def test_delta_t_clamps_at_edges(self, leb_reader):
        """Delta-T should clamp gracefully beyond the table range."""
        jd_start, jd_end = leb_reader.jd_range
        # Should not crash, just return edge values
        dt_before = leb_reader.delta_t(jd_start - 10000)
        dt_after = leb_reader.delta_t(jd_end + 10000)
        assert isinstance(dt_before, float)
        assert isinstance(dt_after, float)


class TestGetStar:
    """Test star catalog lookup."""

    @pytest.mark.integration
    def test_star_not_found(self, leb_reader):
        """get_star raises KeyError for unknown star."""
        with pytest.raises(KeyError, match="Star 99999"):
            leb_reader.get_star(99999)

    @pytest.mark.integration
    def test_star_catalog_populated(self, leb_reader):
        """Star catalog should have entries if generated."""
        # The test fixture includes star catalog
        if leb_reader._stars:
            star_id = next(iter(leb_reader._stars))
            star = leb_reader.get_star(star_id)
            assert 0.0 <= star.ra_j2000 < 360.0
            assert -90.0 <= star.dec_j2000 <= 90.0


class TestLEBReaderClose:
    """Test resource cleanup."""

    @pytest.mark.integration
    def test_close_releases_resources(self, test_leb_file):
        """close() sets internal state to None."""
        reader = LEBReader(test_leb_file)
        reader.close()
        assert reader._mm is None
        assert reader._file is None

    @pytest.mark.integration
    def test_double_close_is_safe(self, test_leb_file):
        """Calling close() twice should not crash."""
        reader = LEBReader(test_leb_file)
        reader.close()
        reader.close()  # Should not raise
