"""
LEB vs Skyfield Comparison: Stations and Retrogrades.

Validates station-finding and retrograde detection.

NOTE: These tests are skipped because the required functions are not yet
implemented in libephemeris. They will be enabled when the station/retrograde
API is added.
"""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.skip(reason="Station/retrograde functions not yet implemented")

import libephemeris as ephem
from libephemeris.constants import MARS, JUPITER, SATURN

from .conftest import TOLS, CompareHelper, year_to_jd


STATION_BODIES = [
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
]


class TestFindStation:
    """Station timing precision."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", STATION_BODIES)
    def test_find_station_ut(
        self, compare: CompareHelper, body_id: int, body_name: str
    ):
        """Station timing matches within tolerance."""
        jd_start = year_to_jd(2024)
        jd_end = year_to_jd(2027)

        # Find stations in range
        ref_result = compare.skyfield(
            ephem.find_station_ut, body_id, jd_start, 0, jd_end
        )
        leb_result = compare.leb(ephem.find_station_ut, body_id, jd_start, 0, jd_end)

        if ref_result[0] != 0 and leb_result[0] != 0:
            ref_jd = ref_result[1][0]
            leb_jd = leb_result[1][0]
            diff_sec = abs(ref_jd - leb_jd) * 86400.0

            assert diff_sec < TOLS.STATION_TIMING_SEC, (
                f"{body_name} station: timing diff = {diff_sec:.1f}s"
            )


class TestIsRetrograde:
    """Retrograde detection."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", STATION_BODIES)
    def test_is_retrograde(self, compare: CompareHelper, body_id: int, body_name: str):
        """Retrograde status matches exactly."""
        from .conftest import generate_test_dates

        dates = generate_test_dates(50, year_to_jd(2020), year_to_jd(2030))
        mismatches = []

        for jd in dates:
            ref = compare.skyfield(ephem.is_retrograde, jd, body_id)
            leb = compare.leb(ephem.is_retrograde, jd, body_id)

            if ref != leb:
                mismatches.append((jd, ref, leb))

        assert len(mismatches) == 0, (
            f"{body_name}: {len(mismatches)} retrograde mismatches"
        )


class TestStationType:
    """Station type classification."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", STATION_BODIES)
    def test_get_station_type(
        self, compare: CompareHelper, body_id: int, body_name: str
    ):
        """Station type matches exactly."""
        from .conftest import generate_test_dates

        dates = generate_test_dates(50, year_to_jd(2020), year_to_jd(2030))
        mismatches = []

        for jd in dates:
            ref = compare.skyfield(ephem.get_station_type, jd, body_id)
            leb = compare.leb(ephem.get_station_type, jd, body_id)

            if ref != leb:
                mismatches.append((jd, ref, leb))

        assert len(mismatches) == 0, (
            f"{body_name}: {len(mismatches)} station type mismatches"
        )


class TestNextRetrograde:
    """Retrograde period finding."""

    @pytest.mark.leb_compare
    @pytest.mark.parametrize("body_id,body_name", STATION_BODIES)
    def test_next_retrograde_ut(
        self, compare: CompareHelper, body_id: int, body_name: str
    ):
        """Retrograde period timing matches within tolerance."""
        jd_start = year_to_jd(2024)

        ref_result = compare.skyfield(ephem.next_retrograde_ut, body_id, jd_start, 0)
        leb_result = compare.leb(ephem.next_retrograde_ut, body_id, jd_start, 0)

        if ref_result[0] != 0 and leb_result[0] != 0:
            ref_jd = ref_result[1][0]
            leb_jd = leb_result[1][0]
            diff_sec = abs(ref_jd - leb_jd) * 86400.0

            assert diff_sec < TOLS.STATION_TIMING_SEC, (
                f"{body_name} retrograde: timing diff = {diff_sec:.1f}s"
            )
