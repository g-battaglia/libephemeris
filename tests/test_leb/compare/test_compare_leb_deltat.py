"""
LEB vs Skyfield Comparison: Delta-T Values.

Validates that LEB Delta-T (TT-UT difference) matches Skyfield reference
across the medium tier range.
"""

from __future__ import annotations

import pytest

import libephemeris as ephem

from .conftest import (
    TOLS,
    CompareHelper,
    generate_test_dates,
    year_to_jd,
)


@pytest.fixture(scope="module")
def deltat_dates() -> list[float]:
    """200 dates across the medium tier range for Delta-T tests."""
    return generate_test_dates(200, year_to_jd(1560), year_to_jd(2640))


@pytest.fixture(scope="module")
def deltat_edge_dates() -> list[float]:
    """Dates near the edges of the medium tier range."""
    start = year_to_jd(1560)
    end = year_to_jd(2640)
    # 10 dates near start, 10 near end
    near_start = generate_test_dates(10, start, year_to_jd(1570), margin=5.0)
    near_end = generate_test_dates(10, year_to_jd(2630), end, margin=5.0)
    return near_start + near_end


class TestDeltaTValues:
    """Direct Delta-T comparison."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    def test_deltat(
        self,
        compare: CompareHelper,
        deltat_dates: list[float],
    ):
        """Delta-T matches Skyfield within tolerance."""
        max_err_sec = 0.0
        worst_jd = 0.0

        for jd in deltat_dates:
            ref_dt = compare.skyfield(ephem.deltat, jd)
            leb_dt = compare.leb(ephem.deltat, jd)

            # Delta-T is in fractional days, convert to seconds
            err_sec = abs(ref_dt - leb_dt) * 86400.0
            if err_sec > max_err_sec:
                max_err_sec = err_sec
                worst_jd = jd

        assert max_err_sec < TOLS.DELTAT_SEC, (
            f"Delta-T: max error = {max_err_sec:.4f}s at JD {worst_jd:.1f}"
        )


class TestDeltaTEdges:
    """Delta-T precision at the edges of the LEB range."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    def test_deltat_edges(
        self,
        compare: CompareHelper,
        deltat_edge_dates: list[float],
    ):
        """Delta-T at range edges matches Skyfield within tolerance."""
        max_err_sec = 0.0
        worst_jd = 0.0

        for jd in deltat_edge_dates:
            ref_dt = compare.skyfield(ephem.deltat, jd)
            leb_dt = compare.leb(ephem.deltat, jd)

            err_sec = abs(ref_dt - leb_dt) * 86400.0
            if err_sec > max_err_sec:
                max_err_sec = err_sec
                worst_jd = jd

        assert max_err_sec < TOLS.DELTAT_SEC, (
            f"Delta-T edge: max error = {max_err_sec:.4f}s at JD {worst_jd:.1f}"
        )


class TestDeltaTStatistics:
    """Statistical analysis of Delta-T errors."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    def test_deltat_distribution(
        self,
        compare: CompareHelper,
        deltat_dates: list[float],
    ):
        """Delta-T error distribution stays within tolerance."""
        errors_sec: list[float] = []

        for jd in deltat_dates:
            ref_dt = compare.skyfield(ephem.deltat, jd)
            leb_dt = compare.leb(ephem.deltat, jd)

            err_sec = abs(ref_dt - leb_dt) * 86400.0
            errors_sec.append(err_sec)

        max_err = max(errors_sec)
        mean_err = sum(errors_sec) / len(errors_sec)

        # Sorted for percentiles
        sorted_errs = sorted(errors_sec)
        p95_idx = int(len(sorted_errs) * 0.95)
        p99_idx = int(len(sorted_errs) * 0.99)
        p95 = sorted_errs[p95_idx] if p95_idx < len(sorted_errs) else max_err
        p99 = sorted_errs[p99_idx] if p99_idx < len(sorted_errs) else max_err

        assert max_err < TOLS.DELTAT_SEC, (
            f"Delta-T: max={max_err:.4f}s, mean={mean_err:.4f}s, "
            f"p95={p95:.4f}s, p99={p99:.4f}s"
        )
