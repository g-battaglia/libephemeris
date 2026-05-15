"""
LEB vs Skyfield Comparison: Nutation Values.

Validates that LEB nutation (dpsi, deps) matches Skyfield reference
across the medium tier range. Uses calc_ut with ECL_NUT.
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import ECL_NUT

from .conftest import (
    TOLS,
    CompareHelper,
    generate_test_dates,
    year_to_jd,
)


@pytest.fixture(scope="module")
def nutation_dates() -> list[float]:
    """100 dates across the medium tier range for nutation tests."""
    return generate_test_dates(100, year_to_jd(1560), year_to_jd(2640))


class TestNutationValues:
    """Direct nutation comparison via ECL_NUT."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    def test_nutation_longitude(
        self,
        compare: CompareHelper,
        nutation_dates: list[float],
    ):
        """Nutation in longitude (dpsi) matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in nutation_dates:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, ECL_NUT, 0)
            leb, _ = compare.leb(ephem.calc_ut, jd, ECL_NUT, 0)

            # result[2] = nutation in longitude (degrees)
            err = abs(ref[2] - leb[2]) * 3600.0  # to arcsec
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS.NUTATION_ARCSEC, (
            f'Nutation dpsi: max error = {max_err:.6f}" at JD {worst_jd:.1f}'
        )

    @pytest.mark.leb_compare
    @pytest.mark.slow
    def test_nutation_obliquity(
        self,
        compare: CompareHelper,
        nutation_dates: list[float],
    ):
        """Nutation in obliquity (deps) matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in nutation_dates:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, ECL_NUT, 0)
            leb, _ = compare.leb(ephem.calc_ut, jd, ECL_NUT, 0)

            # result[3] = nutation in obliquity (degrees)
            err = abs(ref[3] - leb[3]) * 3600.0  # to arcsec
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS.NUTATION_ARCSEC, (
            f'Nutation deps: max error = {max_err:.6f}" at JD {worst_jd:.1f}'
        )

    @pytest.mark.leb_compare
    @pytest.mark.slow
    def test_true_obliquity(
        self,
        compare: CompareHelper,
        nutation_dates: list[float],
    ):
        """True obliquity matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in nutation_dates:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, ECL_NUT, 0)
            leb, _ = compare.leb(ephem.calc_ut, jd, ECL_NUT, 0)

            # result[0] = true obliquity (degrees)
            err = abs(ref[0] - leb[0]) * 3600.0  # to arcsec
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS.NUTATION_ARCSEC, (
            f'True obliquity: max error = {max_err:.6f}" at JD {worst_jd:.1f}'
        )

    @pytest.mark.leb_compare
    @pytest.mark.slow
    def test_mean_obliquity(
        self,
        compare: CompareHelper,
        nutation_dates: list[float],
    ):
        """Mean obliquity matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in nutation_dates:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, ECL_NUT, 0)
            leb, _ = compare.leb(ephem.calc_ut, jd, ECL_NUT, 0)

            # result[1] = mean obliquity (degrees)
            err = abs(ref[1] - leb[1]) * 3600.0  # to arcsec
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS.NUTATION_ARCSEC, (
            f'Mean obliquity: max error = {max_err:.6f}" at JD {worst_jd:.1f}'
        )


class TestNutationAllComponents:
    """Combined test of all 4 nutation/obliquity components."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    def test_all_nutation_components(
        self,
        compare: CompareHelper,
        nutation_dates: list[float],
    ):
        """All nutation/obliquity components match within tolerance."""
        labels = ["true_obliquity", "mean_obliquity", "nut_longitude", "nut_obliquity"]
        max_errs = [0.0] * 4

        for jd in nutation_dates:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, ECL_NUT, 0)
            leb, _ = compare.leb(ephem.calc_ut, jd, ECL_NUT, 0)

            for i in range(4):
                err = abs(ref[i] - leb[i]) * 3600.0
                if err > max_errs[i]:
                    max_errs[i] = err

        for i, label in enumerate(labels):
            assert max_errs[i] < TOLS.NUTATION_ARCSEC, (
                f'{label}: max error = {max_errs[i]:.6f}"'
            )
