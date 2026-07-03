"""Regression tests for two physical-correctness fixes in extinction.py.

These validate the extinction/visibility model by *physical correctness* and
internal continuity (this module has no reference-API oracle):

Fix 1 - Non-monotonic limiting magnitude across the photopic/mesopic sky
        boundary (16 mag/arcsec^2). The mesopic branch used to restart near
        2.0 while the photopic branch reached ~3.4, so a *darker* sky yielded
        a *smaller* limiting magnitude (physically impossible). Limiting
        magnitude must be monotonically non-decreasing as the sky darkens, and
        `is_object_visible` must be monotonic (a darker sky never hides an
        object visible under a brighter sky).

Fix 2 - `calc_rayleigh_coefficient` double-counted the reduced air column by
        multiplying BOTH the pressure factor (P/1013.25) and the altitude
        factor (exp(-h/H)). Passing a physically-consistent site pressure
        together with altitude under-counted the coefficient by ~20%. The
        reduced column must be applied exactly once.
"""

from __future__ import annotations

import math

import pytest

from libephemeris.extinction import (
    DEFAULT_PRESSURE_MBAR,
    SCALE_HEIGHT_M,
    calc_limiting_magnitude_for_sky,
    calc_rayleigh_coefficient,
    is_object_visible,
)

# Fine grid spanning photopic -> mesopic -> dark adaptation states, covering
# both the photopic/mesopic (16) and mesopic/dark (20) transitions.
SKY_GRID = [10.0 + 0.05 * i for i in range(int((22.0 - 10.0) / 0.05) + 1)]


class TestLimitingMagnitudeMonotonic:
    """Fix 1: limiting magnitude vs sky brightness must be non-decreasing."""

    def test_specific_boundary_regression(self):
        """The exact values from the bug report must now be ordered."""
        m_159 = calc_limiting_magnitude_for_sky(15.9)
        m_160 = calc_limiting_magnitude_for_sky(16.0)
        m_165 = calc_limiting_magnitude_for_sky(16.5)
        # A darker sky (higher mag/arcsec^2) must not lower the limiting mag.
        assert m_160 >= m_159 - 1e-9
        assert m_165 >= m_159 - 1e-9

    def test_continuous_at_photopic_mesopic_boundary(self):
        """Mesopic branch is continuous with photopic branch at 16.0."""
        just_below = calc_limiting_magnitude_for_sky(15.999)
        at_boundary = calc_limiting_magnitude_for_sky(16.0)
        # Continuous to within the linear-branch slope over 0.001 mag.
        assert at_boundary == pytest.approx(just_below, abs=0.01)
        # And matches the mesopic comment's ~3.4-4.0 target region.
        assert 3.3 <= at_boundary <= 4.0

    def test_monotonic_non_decreasing_across_range(self):
        """Limiting magnitude is non-decreasing across [10, 22]."""
        lim = [calc_limiting_magnitude_for_sky(s) for s in SKY_GRID]
        for i in range(1, len(lim)):
            assert lim[i] >= lim[i - 1] - 1e-9, (
                f"limiting mag decreased as sky darkened: "
                f"sky {SKY_GRID[i - 1]:.3f}->{SKY_GRID[i]:.3f} gave "
                f"{lim[i - 1]:.4f}->{lim[i]:.4f}"
            )

    def test_monotonic_for_expert_observer(self):
        """Monotonicity holds independent of the (uniform) skill offset."""
        lim = [calc_limiting_magnitude_for_sky(s, observer_skill=4) for s in SKY_GRID]
        for i in range(1, len(lim)):
            assert lim[i] >= lim[i - 1] - 1e-9

    def test_is_object_visible_monotonic(self):
        """A darker sky never hides an object visible under a brighter sky."""
        for mag in (-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0):
            vis = [is_object_visible(mag, s) for s in SKY_GRID]
            for i in range(1, len(vis)):
                if vis[i - 1]:
                    assert vis[i], (
                        f"object mag {mag} visible at sky {SKY_GRID[i - 1]:.3f} "
                        f"but hidden at darker sky {SKY_GRID[i]:.3f}"
                    )

    def test_reported_visibility_case(self):
        """is_object_visible(3.0, .) no longer flips off as the sky darkens."""
        assert is_object_visible(3.0, 15.9) is True
        # Previously False under the buggy mesopic restart; must now stay True.
        assert is_object_visible(3.0, 16.5) is True


class TestRayleighSingleReduction:
    """Fix 2: the reduced air column is applied exactly once."""

    def test_site_pressure_with_altitude_single_reduction(self):
        """Consistent site pressure + altitude -> ~0.114, not the ~0.090 double-count."""
        k = calc_rayleigh_coefficient(pressure_mbar=795.0, altitude_m=2000.0)
        expected_single = 0.1451 * (795.0 / DEFAULT_PRESSURE_MBAR)
        assert k == pytest.approx(expected_single, abs=1e-4)
        assert k == pytest.approx(0.114, abs=0.001)
        # Explicitly reject the old double-counted value.
        double = expected_single * math.exp(-2000.0 / SCALE_HEIGHT_M)
        assert abs(k - double) > 0.02

    def test_sea_level_default_unchanged(self):
        """Sea-level default (P=1013.25, h=0) stays at the reference 0.1451."""
        k = calc_rayleigh_coefficient(
            pressure_mbar=DEFAULT_PRESSURE_MBAR, altitude_m=0.0
        )
        assert k == pytest.approx(0.1451, abs=1e-4)

    def test_altitude_only_reduction_when_pressure_default(self):
        """With pressure left at default, altitude alone drives the reduction."""
        k = calc_rayleigh_coefficient(altitude_m=2000.0)
        expected = 0.1451 * math.exp(-2000.0 / SCALE_HEIGHT_M)
        assert k == pytest.approx(expected, abs=1e-4)

    def test_no_double_count_vs_pressure_only(self):
        """Supplying altitude alongside a site pressure must not reduce further."""
        k_pressure_only = calc_rayleigh_coefficient(pressure_mbar=795.0)
        k_pressure_and_alt = calc_rayleigh_coefficient(
            pressure_mbar=795.0, altitude_m=2000.0
        )
        assert k_pressure_and_alt == pytest.approx(k_pressure_only, abs=1e-9)

    def test_pressure_monotonicity_preserved(self):
        """Higher pressure still means more scattering across the default point."""
        k_low = calc_rayleigh_coefficient(pressure_mbar=800.0)
        k_std = calc_rayleigh_coefficient(pressure_mbar=1013.25)
        k_high = calc_rayleigh_coefficient(pressure_mbar=1050.0)
        assert k_low < k_std < k_high

    def test_altitude_monotonicity_preserved(self):
        """Higher altitude (default pressure) still means less scattering."""
        k_sea = calc_rayleigh_coefficient(altitude_m=0.0)
        k_mid = calc_rayleigh_coefficient(altitude_m=2000.0)
        k_high = calc_rayleigh_coefficient(altitude_m=4000.0)
        assert k_sea > k_mid > k_high
