"""Tests for nutation and obliquity via ECL_NUT body calculation."""

from __future__ import annotations

import math
import pytest
import libephemeris as swe
from libephemeris.constants import (
    ECL_NUT,
    FLG_SWIEPH,
)

JD_J2000 = 2451545.0


@pytest.mark.unit
class TestNutationBody:
    """Test ECL_NUT pseudo-body returns nutation/obliquity values."""

    def test_ecl_nut_returns_tuple(self):
        """ECL_NUT returns 6-element position tuple."""
        result, flags = swe.calc_ut(JD_J2000, ECL_NUT, FLG_SWIEPH)
        assert len(result) == 6

    def test_ecl_nut_true_obliquity(self):
        """Index 0 = true obliquity (should be ~23.4° at J2000)."""
        result, _ = swe.calc_ut(JD_J2000, ECL_NUT, FLG_SWIEPH)
        eps_true = result[0]
        assert 23.0 < eps_true < 24.0, f"True obliquity {eps_true} out of range"

    def test_ecl_nut_mean_obliquity(self):
        """Index 1 = mean obliquity (should be ~23.4° at J2000)."""
        result, _ = swe.calc_ut(JD_J2000, ECL_NUT, FLG_SWIEPH)
        eps_mean = result[1]
        assert 23.0 < eps_mean < 24.0, f"Mean obliquity {eps_mean} out of range"

    def test_ecl_nut_nutation_longitude(self):
        """Index 2 = nutation in longitude (small, < 0.01°)."""
        result, _ = swe.calc_ut(JD_J2000, ECL_NUT, FLG_SWIEPH)
        nut_lon = result[2]
        assert abs(nut_lon) < 0.01, f"Nutation lon {nut_lon} too large"

    def test_ecl_nut_nutation_obliquity(self):
        """Index 3 = nutation in obliquity (small, < 0.01°)."""
        result, _ = swe.calc_ut(JD_J2000, ECL_NUT, FLG_SWIEPH)
        nut_obl = result[3]
        assert abs(nut_obl) < 0.01, f"Nutation obl {nut_obl} too large"

    def test_true_minus_mean_equals_nut_obliquity(self):
        """True obliquity = mean obliquity + nutation in obliquity."""
        result, _ = swe.calc_ut(JD_J2000, ECL_NUT, FLG_SWIEPH)
        eps_true = result[0]
        eps_mean = result[1]
        nut_obl = result[3]
        assert eps_true == pytest.approx(eps_mean + nut_obl, abs=1e-10)

    @pytest.mark.parametrize(
        "jd",
        [
            2415020.0,  # 1900
            2440587.5,  # 1970
            JD_J2000,  # 2000
            2458849.5,  # 2020
            2460000.0,  # 2023
        ],
    )
    def test_ecl_nut_various_dates(self, jd):
        """Nutation values are reasonable at various dates."""
        result, _ = swe.calc_ut(jd, ECL_NUT, FLG_SWIEPH)
        eps_true = result[0]
        eps_mean = result[1]
        nut_lon = result[2]
        nut_obl = result[3]
        # Obliquity range
        assert 22.0 < eps_true < 24.5
        assert 22.0 < eps_mean < 24.5
        # Nutation small
        assert abs(nut_lon) < 0.01
        assert abs(nut_obl) < 0.01

    def test_mean_obliquity_decreasing(self):
        """Mean obliquity decreases over centuries (known secular trend)."""
        result_1900, _ = swe.calc_ut(2415020.0, ECL_NUT, FLG_SWIEPH)
        result_2000, _ = swe.calc_ut(JD_J2000, ECL_NUT, FLG_SWIEPH)
        result_2100, _ = swe.calc_ut(JD_J2000 + 36525.0, ECL_NUT, FLG_SWIEPH)
        # Secular decrease ~47" per century
        assert result_1900[1] > result_2000[1] > result_2100[1]

    def test_nutation_period_18_6_years(self):
        """Nutation in longitude has ~18.6 year period (lunar nodal cycle)."""
        # Sample at 0, ~9.3 yr, ~18.6 yr
        jd0 = JD_J2000
        jd_half = jd0 + 9.3 * 365.25
        jd_full = jd0 + 18.6 * 365.25
        r0, _ = swe.calc_ut(jd0, ECL_NUT, FLG_SWIEPH)
        r_half, _ = swe.calc_ut(jd_half, ECL_NUT, FLG_SWIEPH)
        r_full, _ = swe.calc_ut(jd_full, ECL_NUT, FLG_SWIEPH)
        # After ~18.6 years, nutation in lon should return near initial value
        # Difference at full cycle should be less than at half cycle
        diff_half = abs(r0[2] - r_half[2])
        diff_full = abs(r0[2] - r_full[2])
        assert diff_full < diff_half, "Nutation should return near initial after 18.6yr"

    def test_j2000_obliquity_known(self):
        """J2000 mean obliquity should be ~23.4392911° (IAU 2006)."""
        result, _ = swe.calc_ut(JD_J2000, ECL_NUT, FLG_SWIEPH)
        eps_mean = result[1]
        # Known value: 23° 26' 21.448" = 23.439291°
        assert eps_mean == pytest.approx(23.4393, abs=0.001)

    def test_all_values_finite(self):
        """All 6 return values are finite floats."""
        result, _ = swe.calc_ut(JD_J2000, ECL_NUT, FLG_SWIEPH)
        for i, v in enumerate(result):
            assert math.isfinite(v), f"Index {i} is not finite: {v}"
            assert isinstance(v, float), f"Index {i} is not float: {type(v)}"
