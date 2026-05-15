"""
LEB vs Skyfield Comparison: Elongation Helpers.

Validates elongation functions that call calc_ut internally.
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import (
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
)

from .conftest import TOLS, CompareHelper, year_to_jd, generate_test_dates, angular_diff


ELONGATION_BODIES = [
    (MERCURY, "Mercury"),
    (VENUS, "Venus"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
]


@pytest.fixture(scope="module")
def elongation_dates():
    """50 dates in 2020-2030 for elongation tests."""
    return generate_test_dates(50, year_to_jd(2020), year_to_jd(2030))


class TestElongationAngle:
    """Elongation angle precision."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ELONGATION_BODIES)
    def test_elongation_angle(
        self,
        compare: CompareHelper,
        elongation_dates: list[float],
        body_id: int,
        body_name: str,
    ):
        """Elongation angle matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in elongation_dates:
            ref = compare.skyfield(ephem.get_elongation_from_sun, jd, body_id)
            leb = compare.leb(ephem.get_elongation_from_sun, jd, body_id)

            err = angular_diff(ref[0], leb[0]) * 3600.0
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS.ELONGATION_ARCSEC, (
            f'{body_name}: max elongation error = {max_err:.4f}" at JD {worst_jd:.1f}'
        )


class TestElongationClassification:
    """Morning/evening star classification."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ELONGATION_BODIES)
    def test_morning_star(
        self,
        compare: CompareHelper,
        elongation_dates: list[float],
        body_id: int,
        body_name: str,
    ):
        """Morning star classification matches exactly."""
        mismatches = []

        for jd in elongation_dates:
            ref = compare.skyfield(ephem.is_morning_star, jd, body_id)
            leb = compare.leb(ephem.is_morning_star, jd, body_id)

            if ref != leb:
                mismatches.append((jd, ref, leb))

        assert len(mismatches) == 0, (
            f"{body_name}: {len(mismatches)} morning star mismatches"
        )

    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ELONGATION_BODIES)
    def test_evening_star(
        self,
        compare: CompareHelper,
        elongation_dates: list[float],
        body_id: int,
        body_name: str,
    ):
        """Evening star classification matches exactly."""
        mismatches = []

        for jd in elongation_dates:
            ref = compare.skyfield(ephem.is_evening_star, jd, body_id)
            leb = compare.leb(ephem.is_evening_star, jd, body_id)

            if ref != leb:
                mismatches.append((jd, ref, leb))

        assert len(mismatches) == 0, (
            f"{body_name}: {len(mismatches)} evening star mismatches"
        )

    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ELONGATION_BODIES)
    def test_elongation_type(
        self,
        compare: CompareHelper,
        elongation_dates: list[float],
        body_id: int,
        body_name: str,
    ):
        """Elongation type matches exactly."""
        mismatches = []

        for jd in elongation_dates:
            ref = compare.skyfield(ephem.get_elongation_type, jd, body_id)
            leb = compare.leb(ephem.get_elongation_type, jd, body_id)

            if ref != leb:
                mismatches.append((jd, ref, leb))

        assert len(mismatches) == 0, (
            f"{body_name}: {len(mismatches)} elongation type mismatches"
        )
