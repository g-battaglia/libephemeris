"""
Tests for Besselian element time derivative functions in libephemeris.

Tests the hourly rates of change (dx/dt, dy/dt, dd/dt, dl1/dt, dl2/dt, dmu/dt)
of Besselian elements for solar eclipses.

These derivatives are essential for interpolating element values during an
eclipse and for calculating eclipse contact times.

Reference data sources:
- NASA Eclipse website (eclipse.gsfc.nasa.gov)
- Explanatory Supplement to the Astronomical Almanac
- Meeus "Astronomical Algorithms" Chapter 54
"""

import math
from libephemeris import (
    julday,
    calc_besselian_x,
    calc_besselian_y,
    calc_besselian_d,
    calc_besselian_l1,
    calc_besselian_l2,
    calc_besselian_mu,
    calc_besselian_dx_dt,
    calc_besselian_dy_dt,
    calc_besselian_dd_dt,
    calc_besselian_dl1_dt,
    calc_besselian_dl2_dt,
    calc_besselian_dmu_dt,
    FLG_SWIEPH,
)


class TestBesselianDerivativesBasicFunctionality:
    """Test basic functionality of Besselian element derivative functions."""

    def test_dx_dt_function_exists(self):
        """Test that calc_besselian_dx_dt function exists and is callable."""
        from libephemeris.eclipse import calc_besselian_dx_dt

        assert callable(calc_besselian_dx_dt)

    def test_dy_dt_function_exists(self):
        """Test that calc_besselian_dy_dt function exists and is callable."""
        from libephemeris.eclipse import calc_besselian_dy_dt

        assert callable(calc_besselian_dy_dt)

    def test_dd_dt_function_exists(self):
        """Test that calc_besselian_dd_dt function exists and is callable."""
        from libephemeris.eclipse import calc_besselian_dd_dt

        assert callable(calc_besselian_dd_dt)

    def test_dl1_dt_function_exists(self):
        """Test that calc_besselian_dl1_dt function exists and is callable."""
        from libephemeris.eclipse import calc_besselian_dl1_dt

        assert callable(calc_besselian_dl1_dt)

    def test_dl2_dt_function_exists(self):
        """Test that calc_besselian_dl2_dt function exists and is callable."""
        from libephemeris.eclipse import calc_besselian_dl2_dt

        assert callable(calc_besselian_dl2_dt)

    def test_dmu_dt_function_exists(self):
        """Test that calc_besselian_dmu_dt function exists and is callable."""
        from libephemeris.eclipse import calc_besselian_dmu_dt

        assert callable(calc_besselian_dmu_dt)

    def test_all_functions_exported_from_main_module(self):
        """Test that all derivative functions are exported from main module."""
        from libephemeris import (
            calc_besselian_dx_dt,
            calc_besselian_dy_dt,
            calc_besselian_dd_dt,
            calc_besselian_dl1_dt,
            calc_besselian_dl2_dt,
            calc_besselian_dmu_dt,
        )

        assert callable(calc_besselian_dx_dt)
        assert callable(calc_besselian_dy_dt)
        assert callable(calc_besselian_dd_dt)
        assert callable(calc_besselian_dl1_dt)
        assert callable(calc_besselian_dl2_dt)
        assert callable(calc_besselian_dmu_dt)

    def test_all_functions_return_float(self):
        """Test that all derivative functions return float values."""
        jd = julday(2024, 4, 8, 18.0)

        assert isinstance(calc_besselian_dx_dt(jd), float)
        assert isinstance(calc_besselian_dy_dt(jd), float)
        assert isinstance(calc_besselian_dd_dt(jd), float)
        assert isinstance(calc_besselian_dl1_dt(jd), float)
        assert isinstance(calc_besselian_dl2_dt(jd), float)
        assert isinstance(calc_besselian_dmu_dt(jd), float)

    def test_all_functions_accept_flags_parameter(self):
        """Test that all functions accept optional flags parameter."""
        jd = julday(2024, 4, 8, 18.0)

        # Should not raise any exceptions
        calc_besselian_dx_dt(jd, flags=FLG_SWIEPH)
        calc_besselian_dy_dt(jd, flags=FLG_SWIEPH)
        calc_besselian_dd_dt(jd, flags=FLG_SWIEPH)
        calc_besselian_dl1_dt(jd, flags=FLG_SWIEPH)
        calc_besselian_dl2_dt(jd, flags=FLG_SWIEPH)
        calc_besselian_dmu_dt(jd, flags=FLG_SWIEPH)


class TestDerivativesConsistencyWithBaseElements:
    """Test that derivatives are consistent with base element calculations."""

    def test_dx_dt_matches_numerical_derivative(self):
        """Test dx/dt matches a numerical derivative calculated independently."""
        jd = julday(2024, 4, 8, 18.0)
        h = 0.5 / 24.0  # 30 minutes in days

        # Independent numerical derivative
        x_minus = calc_besselian_x(jd - h)
        x_plus = calc_besselian_x(jd + h)
        expected_dx_dt = (x_plus - x_minus) / (2 * h) / 24.0  # per hour

        actual_dx_dt = calc_besselian_dx_dt(jd)

        # Allow 1% relative tolerance due to different step sizes
        assert abs(actual_dx_dt - expected_dx_dt) < abs(expected_dx_dt) * 0.01 + 1e-8

    def test_dy_dt_matches_numerical_derivative(self):
        """Test dy/dt matches a numerical derivative calculated independently."""
        jd = julday(2024, 4, 8, 18.0)
        h = 0.5 / 24.0  # 30 minutes in days

        y_minus = calc_besselian_y(jd - h)
        y_plus = calc_besselian_y(jd + h)
        expected_dy_dt = (y_plus - y_minus) / (2 * h) / 24.0

        actual_dy_dt = calc_besselian_dy_dt(jd)

        assert abs(actual_dy_dt - expected_dy_dt) < abs(expected_dy_dt) * 0.01 + 1e-8

    def test_dd_dt_matches_numerical_derivative(self):
        """Test dd/dt matches a numerical derivative calculated independently."""
        jd = julday(2024, 4, 8, 18.0)
        h = 0.5 / 24.0  # 30 minutes in days

        d_minus = calc_besselian_d(jd - h)
        d_plus = calc_besselian_d(jd + h)
        expected_dd_dt = (d_plus - d_minus) / (2 * h) / 24.0

        actual_dd_dt = calc_besselian_dd_dt(jd)

        assert abs(actual_dd_dt - expected_dd_dt) < abs(expected_dd_dt) * 0.01 + 1e-8

    def test_dl1_dt_matches_numerical_derivative(self):
        """Test dl1/dt matches a numerical derivative calculated independently."""
        jd = julday(2024, 4, 8, 18.0)
        h = 0.5 / 24.0  # 30 minutes in days

        l1_minus = calc_besselian_l1(jd - h)
        l1_plus = calc_besselian_l1(jd + h)
        expected_dl1_dt = (l1_plus - l1_minus) / (2 * h) / 24.0

        actual_dl1_dt = calc_besselian_dl1_dt(jd)

        # l1 changes very slowly, use absolute tolerance
        assert abs(actual_dl1_dt - expected_dl1_dt) < 1e-6

    def test_dl2_dt_matches_numerical_derivative(self):
        """Test dl2/dt matches a numerical derivative calculated independently."""
        jd = julday(2024, 4, 8, 18.0)
        h = 0.5 / 24.0  # 30 minutes in days

        l2_minus = calc_besselian_l2(jd - h)
        l2_plus = calc_besselian_l2(jd + h)
        expected_dl2_dt = (l2_plus - l2_minus) / (2 * h) / 24.0

        actual_dl2_dt = calc_besselian_dl2_dt(jd)

        # l2 changes very slowly, use absolute tolerance
        assert abs(actual_dl2_dt - expected_dl2_dt) < 1e-7

    def test_dmu_dt_matches_numerical_derivative(self):
        """Test dmu/dt matches a numerical derivative calculated independently."""
        jd = julday(2024, 4, 8, 18.0)
        h = 0.5 / 24.0  # 30 minutes in days

        mu_minus = calc_besselian_mu(jd - h)
        mu_plus = calc_besselian_mu(jd + h)

        # Handle wraparound
        delta_mu = mu_plus - mu_minus
        if delta_mu < -180.0:
            delta_mu += 360.0
        elif delta_mu > 180.0:
            delta_mu -= 360.0

        expected_dmu_dt = delta_mu / (2 * h) / 24.0

        actual_dmu_dt = calc_besselian_dmu_dt(jd)

        # Allow 0.1% relative tolerance
        assert abs(actual_dmu_dt - expected_dmu_dt) < abs(expected_dmu_dt) * 0.001


class TestDerivativesPhysicalReasonableness:
    """Test that derivatives are within physically reasonable ranges."""

    def test_dx_dt_in_reasonable_range(self):
        """Test that dx/dt is within reasonable range during eclipse."""
        # April 8, 2024 total solar eclipse
        jd = julday(2024, 4, 8, 18.0)
        dx_dt = calc_besselian_dx_dt(jd)

        # dx/dt should typically be less than 1 Earth radius per hour
        # The Moon moves about 0.5 deg/hour relative to Sun
        assert abs(dx_dt) < 1.0, f"dx/dt = {dx_dt} seems unreasonably large"

    def test_dy_dt_in_reasonable_range(self):
        """Test that dy/dt is within reasonable range during eclipse."""
        jd = julday(2024, 4, 8, 18.0)
        dy_dt = calc_besselian_dy_dt(jd)

        # dy/dt represents N-S motion, typically < 1 Earth radius per hour
        assert abs(dy_dt) < 1.0, f"dy/dt = {dy_dt} seems unreasonably large"

    def test_dd_dt_in_reasonable_range(self):
        """Test that dd/dt is within reasonable range."""
        jd = julday(2024, 4, 8, 18.0)
        dd_dt = calc_besselian_dd_dt(jd)

        # Shadow axis declination changes slowly (< 1 deg/hour typically)
        assert abs(dd_dt) < 1.0, (
            f"dd/dt = {dd_dt} degrees/hour seems unreasonably large"
        )

    def test_dl1_dt_in_reasonable_range(self):
        """Test that dl1/dt is within reasonable range."""
        jd = julday(2024, 4, 8, 18.0)
        dl1_dt = calc_besselian_dl1_dt(jd)

        # Penumbral radius changes very slowly (order of 1e-3 Earth radii/hour)
        assert abs(dl1_dt) < 0.1, f"dl1/dt = {dl1_dt} seems unreasonably large"

    def test_dl2_dt_in_reasonable_range(self):
        """Test that dl2/dt is within reasonable range."""
        jd = julday(2024, 4, 8, 18.0)
        dl2_dt = calc_besselian_dl2_dt(jd)

        # Umbral radius changes very slowly (order of 1e-4 Earth radii/hour)
        assert abs(dl2_dt) < 0.01, f"dl2/dt = {dl2_dt} seems unreasonably large"

    def test_dmu_dt_near_earth_rotation_rate(self):
        """Test that dmu/dt is approximately Earth's rotation rate."""
        jd = julday(2024, 4, 8, 18.0)
        dmu_dt = calc_besselian_dmu_dt(jd)

        # dmu/dt should be approximately 15 degrees/hour (Earth's rotation)
        # with small variations due to shadow axis motion (typically ±0.1 deg/hr)
        assert 14.5 < dmu_dt < 15.5, f"dmu/dt = {dmu_dt}, expected ~15 degrees/hour"


class TestDerivativesDuringEclipses:
    """Test derivatives during known solar eclipses."""

    def test_derivatives_april_2024_eclipse(self):
        """Test derivatives during April 8, 2024 total solar eclipse."""
        jd_max = julday(2024, 4, 8, 18.3)

        dx_dt = calc_besselian_dx_dt(jd_max)
        dy_dt = calc_besselian_dy_dt(jd_max)
        dd_dt = calc_besselian_dd_dt(jd_max)
        dl1_dt = calc_besselian_dl1_dt(jd_max)
        dl2_dt = calc_besselian_dl2_dt(jd_max)
        dmu_dt = calc_besselian_dmu_dt(jd_max)

        # All should be finite
        assert math.isfinite(dx_dt)
        assert math.isfinite(dy_dt)
        assert math.isfinite(dd_dt)
        assert math.isfinite(dl1_dt)
        assert math.isfinite(dl2_dt)
        assert math.isfinite(dmu_dt)

    def test_derivatives_october_2023_annular_eclipse(self):
        """Test derivatives during October 14, 2023 annular eclipse."""
        jd_max = julday(2023, 10, 14, 18.0)

        dx_dt = calc_besselian_dx_dt(jd_max)
        dy_dt = calc_besselian_dy_dt(jd_max)
        dd_dt = calc_besselian_dd_dt(jd_max)
        dl1_dt = calc_besselian_dl1_dt(jd_max)
        dl2_dt = calc_besselian_dl2_dt(jd_max)
        dmu_dt = calc_besselian_dmu_dt(jd_max)

        # All should be finite
        assert math.isfinite(dx_dt)
        assert math.isfinite(dy_dt)
        assert math.isfinite(dd_dt)
        assert math.isfinite(dl1_dt)
        assert math.isfinite(dl2_dt)
        assert math.isfinite(dmu_dt)

    def test_derivatives_december_2021_eclipse(self):
        """Test derivatives during December 4, 2021 total eclipse (Antarctica)."""
        jd_max = julday(2021, 12, 4, 7.5)

        dx_dt = calc_besselian_dx_dt(jd_max)
        dy_dt = calc_besselian_dy_dt(jd_max)
        dd_dt = calc_besselian_dd_dt(jd_max)
        dl1_dt = calc_besselian_dl1_dt(jd_max)
        dl2_dt = calc_besselian_dl2_dt(jd_max)
        dmu_dt = calc_besselian_dmu_dt(jd_max)

        # All should be finite
        assert math.isfinite(dx_dt)
        assert math.isfinite(dy_dt)
        assert math.isfinite(dd_dt)
        assert math.isfinite(dl1_dt)
        assert math.isfinite(dl2_dt)
        assert math.isfinite(dmu_dt)


class TestDerivativesNumericalStability:
    """Test numerical stability of derivative calculations."""

    def test_derivatives_deterministic(self):
        """Test that derivatives give consistent results for same input."""
        jd = julday(2024, 4, 8, 18.0)

        # Calculate twice
        dx_dt_1 = calc_besselian_dx_dt(jd)
        dx_dt_2 = calc_besselian_dx_dt(jd)

        dmu_dt_1 = calc_besselian_dmu_dt(jd)
        dmu_dt_2 = calc_besselian_dmu_dt(jd)

        assert dx_dt_1 == dx_dt_2, "dx/dt should be deterministic"
        assert dmu_dt_1 == dmu_dt_2, "dmu/dt should be deterministic"

    def test_derivatives_continuous(self):
        """Test that derivatives change smoothly over time."""
        jd_start = julday(2024, 4, 8, 17.0)

        # Sample at small intervals
        dx_values = []
        dmu_values = []
        for i in range(10):
            jd = jd_start + i * 0.01  # ~15 minute intervals
            dx_values.append(calc_besselian_dx_dt(jd))
            dmu_values.append(calc_besselian_dmu_dt(jd))

        # Check for smooth variation (no discontinuities)
        for i in range(1, len(dx_values)):
            delta_dx = abs(dx_values[i] - dx_values[i - 1])
            delta_dmu = abs(dmu_values[i] - dmu_values[i - 1])

            # Changes should be small over ~15 minutes
            assert delta_dx < 0.1, f"dx/dt discontinuity: {delta_dx}"
            assert delta_dmu < 0.1, f"dmu/dt discontinuity: {delta_dmu}"

    def test_derivatives_away_from_eclipse(self):
        """Test derivatives at arbitrary times (not during eclipse)."""
        jd_random = julday(2024, 7, 15, 12.0)

        dx_dt = calc_besselian_dx_dt(jd_random)
        dy_dt = calc_besselian_dy_dt(jd_random)
        dd_dt = calc_besselian_dd_dt(jd_random)
        dl1_dt = calc_besselian_dl1_dt(jd_random)
        dl2_dt = calc_besselian_dl2_dt(jd_random)
        dmu_dt = calc_besselian_dmu_dt(jd_random)

        # All should be finite even away from eclipse
        assert math.isfinite(dx_dt)
        assert math.isfinite(dy_dt)
        assert math.isfinite(dd_dt)
        assert math.isfinite(dl1_dt)
        assert math.isfinite(dl2_dt)
        assert math.isfinite(dmu_dt)


class TestDerivativesEdgeCases:
    """Test edge cases and special conditions."""

    def test_derivatives_at_various_dates(self):
        """Test derivatives at various dates including equinoxes and solstices."""
        test_dates = [
            julday(2000, 1, 1, 12.0),  # Y2K
            julday(2024, 3, 20, 12.0),  # Vernal equinox
            julday(2024, 6, 21, 12.0),  # Summer solstice
            julday(2024, 9, 22, 12.0),  # Autumn equinox
            julday(2024, 12, 21, 12.0),  # Winter solstice
        ]

        for jd in test_dates:
            dx_dt = calc_besselian_dx_dt(jd)
            dy_dt = calc_besselian_dy_dt(jd)
            dd_dt = calc_besselian_dd_dt(jd)
            dl1_dt = calc_besselian_dl1_dt(jd)
            dl2_dt = calc_besselian_dl2_dt(jd)
            dmu_dt = calc_besselian_dmu_dt(jd)

            # All should be finite
            assert math.isfinite(dx_dt), f"dx/dt not finite at JD {jd}"
            assert math.isfinite(dy_dt), f"dy/dt not finite at JD {jd}"
            assert math.isfinite(dd_dt), f"dd/dt not finite at JD {jd}"
            assert math.isfinite(dl1_dt), f"dl1/dt not finite at JD {jd}"
            assert math.isfinite(dl2_dt), f"dl2/dt not finite at JD {jd}"
            assert math.isfinite(dmu_dt), f"dmu/dt not finite at JD {jd}"

    def test_dmu_dt_wraparound_handling(self):
        """Test that dmu/dt handles the 0/360 degree boundary correctly."""
        # Find a time where mu is near 360 degrees
        # This requires testing across a range of times
        for hour in range(24):
            jd = julday(2024, 4, 8, hour)
            mu = calc_besselian_mu(jd)
            dmu_dt = calc_besselian_dmu_dt(jd)

            # dmu/dt should always be approximately 15 deg/hour (positive)
            # regardless of where mu is
            assert 14.0 < dmu_dt < 16.0, (
                f"dmu/dt = {dmu_dt} at mu = {mu:.1f} deg, expected ~15 deg/hr"
            )


class TestInterpolationUseCase:
    """Test the use case of interpolating Besselian elements."""

    def test_interpolation_accuracy(self):
        """Test that derivatives improve interpolation accuracy."""
        # Base time
        jd0 = julday(2024, 4, 8, 18.0)
        dt = 0.25  # 15 minutes in hours

        # Get base values and derivatives
        x0 = calc_besselian_x(jd0)
        dx_dt = calc_besselian_dx_dt(jd0)

        # Interpolate x at t0 + dt using Taylor series
        jd1 = jd0 + dt / 24.0  # Convert hours to days
        x_interpolated = x0 + dx_dt * dt

        # Actual value
        x_actual = calc_besselian_x(jd1)

        # Error should be small (second-order in dt)
        error = abs(x_interpolated - x_actual)
        assert error < 0.01, (
            f"Interpolation error = {error:.6f} Earth radii, "
            f"expected < 0.01 for 15-minute extrapolation"
        )

    def test_all_elements_interpolation(self):
        """Test interpolation for all Besselian elements."""
        jd0 = julday(2024, 4, 8, 18.0)
        dt = 0.25  # 15 minutes in hours
        jd1 = jd0 + dt / 24.0

        # Get base values and derivatives
        elements = {
            "x": (
                calc_besselian_x(jd0),
                calc_besselian_dx_dt(jd0),
                calc_besselian_x(jd1),
            ),
            "y": (
                calc_besselian_y(jd0),
                calc_besselian_dy_dt(jd0),
                calc_besselian_y(jd1),
            ),
            "d": (
                calc_besselian_d(jd0),
                calc_besselian_dd_dt(jd0),
                calc_besselian_d(jd1),
            ),
            "l1": (
                calc_besselian_l1(jd0),
                calc_besselian_dl1_dt(jd0),
                calc_besselian_l1(jd1),
            ),
            "l2": (
                calc_besselian_l2(jd0),
                calc_besselian_dl2_dt(jd0),
                calc_besselian_l2(jd1),
            ),
            "mu": (
                calc_besselian_mu(jd0),
                calc_besselian_dmu_dt(jd0),
                calc_besselian_mu(jd1),
            ),
        }

        for name, (val0, deriv, actual) in elements.items():
            interpolated = val0 + deriv * dt

            # Handle mu wraparound
            if name == "mu":
                interpolated = interpolated % 360.0
                error = abs(interpolated - actual)
                if error > 180:
                    error = 360 - error
            else:
                error = abs(interpolated - actual)

            # Allow larger tolerance for mu due to wraparound effects
            tolerance = 0.1 if name == "mu" else 0.01
            assert error < tolerance, (
                f"{name}: interpolation error = {error:.6f}, "
                f"interpolated = {interpolated:.4f}, actual = {actual:.4f}"
            )
