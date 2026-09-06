"""The visibility allowance of the solar local circumstances, contribution by
contribution, and the two obscuration regimes behind ``attr[2]``.

The allowance that decides whether any part of an eclipsed Sun can stand above
an observer's horizon is a sum of three published quantities, each with its own
source and its own dependence on the observer's height:

(a) the astronomical refraction of a ray arriving on the horizon (Bennett 1982);
(b) the geometric dip of the horizon on the IERS equatorial radius (Bomford
    1980; IERS Conventions (2010), Table 1.1);
(c) the terrestrial refraction that lightens that dip (Bomford 1980, through
    :func:`libephemeris.refraction.calc_dip`).

Each is pinned here on its own, at a configuration whose expected value is
published, so a reader can check one at a time rather than a single number
nobody can account for.
"""

from __future__ import annotations

import math

import pytest

from libephemeris import ECL_PARTIAL, ECL_VISIBLE, sol_eclipse_how, sol_eclipse_where
from libephemeris.eclipse import (
    _dip_refraction_deg,
    _geometric_horizon_dip_deg,
    _horizon_refraction_deg,
    _how_pressure_mbar,
    _overlap_area_fraction,
    _visibility_allowance_deg,
)
from libephemeris.refraction import calc_dip

#: Standard conditions of the horizon fits: 1010 mbar and 10 degrees Celsius.
_FIT_PRESSURE = 1010.0
_TEMPERATURE_C = 10.0
#: Arcminutes per degree, so the pins below read as their published figures.
_ARCMIN = 60.0


class TestHorizonRefraction:
    """Contribution (a): refraction of a ray arriving on the horizon."""

    def test_at_the_conditions_of_the_fit(self):
        """Bennett's form at zero apparent altitude is 34.48 arcminutes.

        ``1 / tan(7.31 / 4.4)`` in arcminutes, the published value of the
        refraction at the horizon for 1010 mbar and 10 C.
        """
        arcmin = _horizon_refraction_deg(_FIT_PRESSURE, _TEMPERATURE_C) * _ARCMIN
        assert arcmin == pytest.approx(34.4775, abs=1e-4)

    def test_scaled_to_the_sea_level_atmosphere(self):
        """At 1013.25 mbar and 10 C the same ray is lifted 34.59 arcminutes."""
        assert _how_pressure_mbar(0.0) == 1013.25
        deg = _horizon_refraction_deg(1013.25, _TEMPERATURE_C)
        assert deg * _ARCMIN == pytest.approx(34.588, abs=1e-3)
        assert deg == pytest.approx(0.5765, abs=1e-4)

    def test_weakens_with_the_observer_height(self):
        """An observer above part of the atmosphere is lifted less.

        The standard atmosphere leaves 88.7 percent of the sea-level pressure
        at 1000 m, and the refraction at the horizon follows it down: about
        four arcminutes less there, about 0.4 at 100 m.
        """
        sea = _horizon_refraction_deg(_how_pressure_mbar(0.0), _TEMPERATURE_C)
        hundred = _horizon_refraction_deg(_how_pressure_mbar(100.0), _TEMPERATURE_C)
        kilometre = _horizon_refraction_deg(_how_pressure_mbar(1000.0), _TEMPERATURE_C)
        assert _how_pressure_mbar(1000.0) / 1013.25 == pytest.approx(0.887, abs=1e-3)
        assert (sea - hundred) * _ARCMIN == pytest.approx(0.41, abs=0.02)
        assert (sea - kilometre) * _ARCMIN == pytest.approx(3.91, abs=0.05)


class TestGeometricHorizonDip:
    """Contribution (b): the dip of the horizon of a sphere."""

    @pytest.mark.parametrize(
        "height_m, arcmin",
        [(1.0, 1.925), (100.0, 19.250), (1000.0, 60.871), (4000.0, 121.719)],
    )
    def test_grows_as_the_root_of_the_height(self, height_m, arcmin):
        """1.925 arcminutes per square root of a metre, on the IERS radius."""
        assert _geometric_horizon_dip_deg(height_m) * _ARCMIN == pytest.approx(
            arcmin, abs=1e-3
        )
        coefficient = (
            _geometric_horizon_dip_deg(height_m) * _ARCMIN / math.sqrt(height_m)
        )
        assert coefficient == pytest.approx(1.925, abs=1e-3)

    @pytest.mark.parametrize("height_m", [0.0, -1.0, -400.0])
    def test_nothing_at_or_below_sea_level(self, height_m):
        """The whole recorded observer grid stands here, and sees no dip."""
        assert _geometric_horizon_dip_deg(height_m) == 0.0


class TestDipRefraction:
    """Contribution (c): the terrestrial refraction that lightens the dip."""

    @pytest.mark.parametrize("height_m", [1.0, 100.0, 1000.0, 4000.0])
    def test_takes_back_about_a_tenth_of_the_dip(self, height_m):
        """Negative, and between an eighth and a seventh of contribution (b)."""
        pressure = _how_pressure_mbar(height_m)
        geometric = _geometric_horizon_dip_deg(height_m)
        correction = _dip_refraction_deg(height_m, pressure, _TEMPERATURE_C)
        assert correction < 0.0
        assert -0.15 < correction / geometric < -0.07

    def test_is_the_correction_the_dip_unit_applies(self):
        """(b) + (c) is the refracted dip of :func:`calc_dip`, to the last bit."""
        for height_m in (1.0, 100.0, 1000.0, 4000.0):
            pressure = _how_pressure_mbar(height_m)
            total = _geometric_horizon_dip_deg(height_m) + _dip_refraction_deg(
                height_m, pressure, _TEMPERATURE_C
            )
            assert total == -calc_dip(height_m, 0.0065, pressure, _TEMPERATURE_C)

    def test_coefficient_of_terrestrial_refraction_is_about_a_quarter(self):
        """``sqrt(1 - k)`` scales the dip; the standard atmosphere gives k=0.258."""
        height_m = 1.0
        pressure = _how_pressure_mbar(height_m)
        geometric = _geometric_horizon_dip_deg(height_m)
        refracted = geometric + _dip_refraction_deg(height_m, pressure, _TEMPERATURE_C)
        k = 1.0 - (refracted / geometric) ** 2
        assert k == pytest.approx(0.258, abs=0.005)

    @pytest.mark.parametrize("height_m", [0.0, -1.0])
    def test_nothing_at_or_below_sea_level(self, height_m):
        assert _dip_refraction_deg(height_m, 1013.25, _TEMPERATURE_C) == 0.0


class TestVisibilityAllowance:
    """The sum of the three, and the gate it decides."""

    def test_is_exactly_the_sum_of_the_three_contributions(self):
        for height_m in (0.0, 1.0, 100.0, 1000.0, 4000.0):
            pressure = _how_pressure_mbar(height_m)
            parts = (
                _horizon_refraction_deg(pressure, _TEMPERATURE_C)
                + _geometric_horizon_dip_deg(height_m)
                + _dip_refraction_deg(height_m, pressure, _TEMPERATURE_C)
            )
            assert (
                _visibility_allowance_deg(height_m, pressure, _TEMPERATURE_C) == parts
            )

    def test_at_sea_level_only_the_refraction_survives(self):
        allowance = _visibility_allowance_deg(0.0, 1013.25, _TEMPERATURE_C)
        assert allowance == _horizon_refraction_deg(1013.25, _TEMPERATURE_C)
        assert allowance == pytest.approx(0.5765, abs=1e-4)

    def test_admits_the_sun_down_to_five_sixths_of_a_degree(self):
        """A limb still counts: the floor is the allowance plus a semidiameter.

        The Sun's apparent radius is about 16 arcminutes, so at sea level the
        centre may sit 0.84 degrees below the horizontal plane and part of the
        disc still be above the horizon.
        """
        allowance = _visibility_allowance_deg(0.0, 1013.25, _TEMPERATURE_C)
        sun_radius = 16.0 / _ARCMIN
        floor = -(allowance + sun_radius)
        assert floor == pytest.approx(-0.843, abs=0.005)
        assert -0.84 + sun_radius + allowance > 0.0
        assert -0.85 + sun_radius + allowance < 0.0

    def test_grows_by_more_than_a_degree_over_four_kilometres(self):
        sea = _visibility_allowance_deg(0.0, _how_pressure_mbar(0.0), _TEMPERATURE_C)
        top = _visibility_allowance_deg(
            4000.0, _how_pressure_mbar(4000.0), _TEMPERATURE_C
        )
        assert sea == pytest.approx(0.576, abs=1e-3)
        assert top == pytest.approx(2.213, abs=1e-3)


class TestVisibilityGateOnBothSides:
    """The gate against the two witnesses on each side of the horizon."""

    # (Julian Day UT, (longitude, latitude, height)) of four instants whose
    # answers straddle the horizon: the first two report a partial phase with
    # the Sun's centre already below the horizontal plane, the last two report
    # nothing because the refracted Sun has set.
    VISIBLE = (
        (2431824.010919301, (179.999, -66.56, 0.0), -0.561867517, 0.014397424),
        (2420188.5088523636, (0.0, -80.0, 0.0), -0.525758311, 0.044117009),
    )
    HIDDEN = (
        (2411890.17757277, (179.999, 45.0, 0.0), -0.580455195),
        (2477034.6796099963, (-90.0, -80.0, 0.0), -0.582517904),
    )

    @pytest.mark.parametrize("jd, geopos, true_alt, apparent_alt", VISIBLE)
    def test_a_limb_above_the_horizon_reports_and_is_visible(
        self, jd, geopos, true_alt, apparent_alt
    ):
        retflag, attr = sol_eclipse_how(jd, geopos)
        assert retflag & ECL_PARTIAL
        assert retflag & ECL_VISIBLE
        assert attr[5] == pytest.approx(true_alt, abs=1e-8)
        assert attr[6] == pytest.approx(apparent_alt, abs=1e-8)
        # The gate itself has room to spare: the limb stands a quarter of a
        # degree above the floor the three contributions put it at.
        allowance = _visibility_allowance_deg(0.0, 1013.25, _TEMPERATURE_C)
        assert attr[5] + 16.0 / _ARCMIN + allowance > 0.25

    @pytest.mark.parametrize("jd, geopos, true_alt", HIDDEN)
    def test_a_sun_below_the_refracted_horizon_reports_nothing(
        self, jd, geopos, true_alt
    ):
        retflag, attr = sol_eclipse_how(jd, geopos)
        assert retflag == 0
        assert not retflag & ECL_VISIBLE
        # The refraction unit floors at the dipped horizon, so the two
        # altitudes coincide bit for bit below it.
        assert attr[5] == pytest.approx(true_alt, abs=1e-8)
        assert attr[6] == attr[5]

    def test_the_deepest_reporting_case_clears_the_gate_by_a_quarter_degree(self):
        """The publishing gate binds long before the allowance does."""
        jd, geopos, true_alt, _apparent = self.VISIBLE[0]
        _retflag, attr = sol_eclipse_how(jd, geopos)
        allowance = _visibility_allowance_deg(0.0, 1013.25, _TEMPERATURE_C)
        sun_radius = 16.0 / _ARCMIN
        assert attr[5] + sun_radius + allowance == pytest.approx(0.28, abs=0.02)


class TestObscurationRegimes:
    """The covered area, while the limbs cross and once one disc is inside."""

    def test_zero_at_external_contact(self):
        assert _overlap_area_fraction(0.26, 0.25, 0.51) == pytest.approx(0.0, abs=1e-12)

    def test_whole_disc_at_internal_contact(self):
        """With the Moon the larger disc the lens closes on the whole Sun."""
        assert _overlap_area_fraction(0.25, 0.26, 0.01) == pytest.approx(1.0, abs=1e-12)

    def test_concentric_equal_discs(self):
        assert _overlap_area_fraction(0.25, 0.25, 0.0) == 1.0

    @pytest.mark.parametrize(
        "r_body, r_moon, separation",
        [(0.2666, 0.2500, 0.2700), (0.2666, 0.2760, 0.1500), (0.2500, 0.1000, 0.2000)],
    )
    def test_lens_area_matches_a_quadrature(self, r_body, r_moon, separation):
        """The closed form against the area integrated strip by strip."""
        steps = 400_000
        lo = max(-r_body, separation - r_moon)
        hi = min(r_body, separation + r_moon)
        width = (hi - lo) / steps
        area = 0.0
        for i in range(steps):
            x = lo + (i + 0.5) * width
            half_body = math.sqrt(max(r_body * r_body - x * x, 0.0))
            dx = separation - x
            half_moon = math.sqrt(max(r_moon * r_moon - dx * dx, 0.0))
            area += 2.0 * min(half_body, half_moon) * width
        expected = area / (math.pi * r_body * r_body)
        got = _overlap_area_fraction(r_body, r_moon, separation)
        assert got == pytest.approx(expected, rel=1e-6)

    @pytest.mark.parametrize("jd", [2400432.5190255954, 2448449.295865731])
    def test_internal_contact_reports_the_squared_ratio(self, jd):
        """One disc inside the other: the covered area is the smaller disc."""
        _retflag, _geopos, attr = sol_eclipse_where(jd)
        assert attr[2] == attr[1] * attr[1]
