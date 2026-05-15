"""
LEB vs Skyfield Comparison: Gauquelin Sectors.

Validates gauquelin_sector which calls calc_ut for planet positions.
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import SUN, MOON, MARS, JUPITER

from .conftest import TOLS, CompareHelper, year_to_jd


GAUQUELIN_BODIES = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
]

TEST_LOCATIONS = [
    ("Rome", 41.9, 12.5, 0),
    ("NewYork", 40.7, -74.0, 0),
    ("Sydney", -33.9, 151.2, 0),
    ("Equator", 0.0, 0.0, 0),
]

TEST_DATES = [year_to_jd(2000), year_to_jd(2024), year_to_jd(2050)]

GAUQUELIN_METHODS = [
    (0, "Placidus"),
    (1, "Koch"),
]


class TestGauquelinSector:
    """Gauquelin sector precision."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", GAUQUELIN_BODIES)
    @pytest.mark.parametrize("name,lat,lon,alt", TEST_LOCATIONS)
    @pytest.mark.parametrize("jd", TEST_DATES)
    @pytest.mark.parametrize("method,method_name", GAUQUELIN_METHODS)
    def test_sector(
        self,
        compare: CompareHelper,
        body_id: int,
        body_name: str,
        name: str,
        lat: float,
        lon: float,
        alt: float,
        jd: float,
        method: int,
        method_name: str,
    ):
        """Gauquelin sector matches within tolerance."""
        geopos = (lat, lon, alt)
        ref_sector = compare.skyfield(
            ephem.gauquelin_sector, jd, body_id, method, geopos
        )
        leb_sector = compare.leb(
            ephem.gauquelin_sector, jd, body_id, method, geopos
        )

        # Function returns float directly (sector number with decimal position)
        if ref_sector > 0 and leb_sector > 0:
            sector_diff = abs(ref_sector - leb_sector)
            assert sector_diff < TOLS.GAUQUELIN_SECTOR, (
                f"{body_name} at {name} ({method_name}): sector diff = {sector_diff:.4f}"
            )
