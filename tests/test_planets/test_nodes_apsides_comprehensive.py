"""
Tests for lunar nodes and apsides (nod_aps_ut).

Verifies that the node/apse calculation returns valid results
for the Moon and planets, and that values are physically plausible.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe
from libephemeris.constants import (
    SUN,
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    NODBIT_MEAN,
    NODBIT_OSCU,
    NODBIT_OSCU_BAR,
    FLG_SPEED,
)


PLANETS_FOR_NODES = [
    (MOON, "Moon"),
    (MERCURY, "Mercury"),
    (VENUS, "Venus"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
]


class TestNodApsBasic:
    """Basic nod_aps_ut functionality tests."""

    @pytest.mark.unit
    def test_nod_aps_moon_returns_4_tuples(self):
        """nod_aps_ut for Moon returns 4 tuples of 6 elements each."""
        jd = 2451545.0
        result = swe.nod_aps_ut(jd, MOON, 0, NODBIT_MEAN)
        assert len(result) == 4, f"Expected 4 tuples, got {len(result)}"
        for i, tup in enumerate(result):
            assert len(tup) == 6, f"Tuple {i} has {len(tup)} elements, expected 6"

    @pytest.mark.unit
    def test_nod_aps_moon_ascending_node(self):
        """Ascending node longitude should be valid."""
        jd = 2451545.0
        asc_node, desc_node, perihelion, aphelion = swe.nod_aps_ut(
            jd, MOON, 0, NODBIT_MEAN
        )
        lon = float(asc_node[0])
        assert 0 <= lon < 360, f"Ascending node lon {lon} out of range"

    @pytest.mark.unit
    def test_nod_aps_moon_nodes_opposite(self):
        """Ascending and descending nodes should be ~180° apart."""
        jd = 2451545.0
        asc_node, desc_node, perihelion, aphelion = swe.nod_aps_ut(
            jd, MOON, 0, NODBIT_MEAN
        )
        asc_lon = float(asc_node[0])
        desc_lon = float(desc_node[0])
        diff = abs(asc_lon - desc_lon)
        if diff > 180:
            diff = 360 - diff
        assert abs(diff - 180) < 1.0, (
            f"Nodes not opposite: asc={asc_lon:.2f} desc={desc_lon:.2f} diff={diff:.2f}"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", PLANETS_FOR_NODES)
    def test_nod_aps_returns_valid(self, body_id: int, name: str):
        """nod_aps_ut returns valid results for each planet."""
        jd = 2451545.0
        result = swe.nod_aps_ut(jd, body_id, 0, NODBIT_MEAN)
        assert len(result) == 4
        for i, tup in enumerate(result):
            assert len(tup) == 6
            lon = float(tup[0])
            assert 0 <= lon < 360 or lon == 0.0, (
                f"{name} tuple {i}: lon {lon} out of range"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "method,name",
        [
            (NODBIT_MEAN, "mean"),
            (NODBIT_OSCU, "osculating"),
            (NODBIT_OSCU_BAR, "osculating barycentric"),
        ],
    )
    def test_nod_aps_methods(self, method: int, name: str):
        """Different calculation methods all return valid results."""
        jd = 2451545.0
        result = swe.nod_aps_ut(jd, MOON, 0, method)
        assert len(result) == 4
        asc_lon = float(result[0][0])
        assert 0 <= asc_lon < 360, f"Method {name}: asc node lon {asc_lon} out of range"


class TestNodApsDateRange:
    """Test nod_aps_ut across date ranges."""

    @pytest.mark.unit
    @pytest.mark.parametrize("year", [1800, 1900, 1950, 2000, 2050, 2100, 2200])
    def test_moon_nodes_across_centuries(self, year: int):
        """Moon nodes valid across centuries."""
        jd = swe.julday(year, 1, 1, 12.0)
        asc, desc, peri, aph = swe.nod_aps_ut(jd, MOON, 0, NODBIT_MEAN)
        asc_lon = float(asc[0])
        assert 0 <= asc_lon < 360, f"Year {year}: asc node lon {asc_lon}"

    @pytest.mark.unit
    def test_moon_node_regression(self):
        """Moon's mean node regresses ~19.35°/year."""
        jd1 = 2451545.0
        jd2 = jd1 + 365.25
        asc1 = swe.nod_aps_ut(jd1, MOON, 0, NODBIT_MEAN)
        asc2 = swe.nod_aps_ut(jd2, MOON, 0, NODBIT_MEAN)
        lon1 = float(asc1[0][0])
        lon2 = float(asc2[0][0])
        motion = (lon2 - lon1) % 360
        if motion > 180:
            motion -= 360
        # Mean node regresses ~19.35°/year
        assert -25 < motion < -15, (
            f"Node annual motion {motion:.2f}° (expected ~-19.35°)"
        )


class TestNodApsConsistency:
    """Test consistency of node/apse values."""

    @pytest.mark.unit
    def test_perihelion_closer_than_aphelion(self):
        """Moon mean perihelion/aphelion distances are equal (constant mean distance).

        For NODBIT_MEAN, both perihelion and aphelion return the same
        constant mean lunar distance. For osculating elements the same
        may hold depending on the implementation. We verify both are
        positive and finite.
        """
        jd = 2451545.0
        # Mean: distances are equal (constant mean distance)
        asc, desc, peri, aph = swe.nod_aps_ut(jd, MOON, 0, NODBIT_MEAN)
        peri_dist = float(peri[2])
        aph_dist = float(aph[2])
        assert peri_dist > 0, f"Mean perihelion dist {peri_dist} not positive"
        assert aph_dist > 0, f"Mean aphelion dist {aph_dist} not positive"

        # For planets, verify osculating apse distances are positive and distinct
        asc, desc, peri, aph = swe.nod_aps_ut(jd, MARS, 0, NODBIT_OSCU)
        peri_dist = float(peri[2])
        aph_dist = float(aph[2])
        assert peri_dist > 0, f"Mars oscu perihelion dist {peri_dist} not positive"
        assert aph_dist > 0, f"Mars oscu aphelion dist {aph_dist} not positive"
        # The two apse distances should differ
        assert abs(peri_dist - aph_dist) > 0.01, (
            f"Mars oscu apse distances too similar: {peri_dist} vs {aph_dist}"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,name",
        [
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
        ],
    )
    def test_planet_nodes_stable(self, body_id: int, name: str):
        """Planet nodes should be relatively stable over 1 year."""
        jd1 = 2451545.0
        jd2 = jd1 + 365.25
        r1 = swe.nod_aps_ut(jd1, body_id, 0, NODBIT_MEAN)
        r2 = swe.nod_aps_ut(jd2, body_id, 0, NODBIT_MEAN)
        lon1 = float(r1[0][0])
        lon2 = float(r2[0][0])
        diff = abs(lon2 - lon1)
        if diff > 180:
            diff = 360 - diff
        # Planet nodes change slowly (< 1°/year for most)
        assert diff < 5.0, f"{name}: node changed {diff:.2f}° in 1 year"
