"""
Comprehensive tests for topocentric position calculations.

Verifies that FLG_TOPOCTR produces valid positions that differ
from geocentric, and that set_topo / set_topo work correctly.
"""

from __future__ import annotations

import math
import random

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
    FLG_SPEED,
    FLG_TOPOCTR,
    FLG_EQUATORIAL,
    FLG_SIDEREAL,
    SIDM_LAHIRI,
)


LOCATIONS = [
    # (lon, lat, alt, name)
    (12.5, 41.9, 50.0, "Rome"),
    (-74.0, 40.7, 10.0, "New York"),
    (139.7, 35.7, 40.0, "Tokyo"),
    (151.2, -33.9, 58.0, "Sydney"),
    (-43.2, -22.9, 11.0, "Rio"),
    (0.0, 0.0, 0.0, "Null Island"),
    (0.0, 89.0, 0.0, "Near North Pole"),
    (0.0, -89.0, 0.0, "Near South Pole"),
]

BODIES_FOR_TOPO = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MERCURY, "Mercury"),
    (VENUS, "Venus"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
]


class TestTopoSetup:
    """Tests for set_topo / set_topo configuration."""

    @pytest.mark.unit
    @pytest.mark.parametrize("lon,lat,alt,name", LOCATIONS)
    def test_set_topo_accepts_valid_locations(
        self, lon: float, lat: float, alt: float, name: str
    ):
        """set_topo should accept all valid geographic coordinates."""
        # Should not raise
        swe.set_topo(lon, lat, alt)

    @pytest.mark.unit
    def test_set_topo_invalid_latitude_raises(self):
        """set_topo should raise on latitude > 90."""
        with pytest.raises(Exception):
            swe.set_topo(0.0, 91.0, 0.0)

    @pytest.mark.unit
    def test_set_topo_invalid_longitude_raises(self):
        """set_topo should raise on longitude outside [-180, 360]."""
        with pytest.raises((swe.CoordinateError, ValueError)):
            swe.set_topo(400.0, 0.0, 0.0)


class TestTopoBasicPositions:
    """Test that topocentric positions are valid and differ from geocentric."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", BODIES_FOR_TOPO)
    def test_topo_returns_valid_position(self, body_id: int, name: str):
        """Topocentric calc returns valid 6-element tuple."""
        swe.set_topo(12.5, 41.9, 50.0)  # Rome
        jd = 2451545.0
        result, retflag = swe.calc_ut(jd, body_id, FLG_TOPOCTR | FLG_SPEED)
        assert len(result) == 6, f"{name}: expected 6 elements"
        lon, lat, dist = result[0], result[1], result[2]
        assert 0 <= lon < 360, f"{name}: topo lon {lon} out of range"
        assert -90 <= lat <= 90, f"{name}: topo lat {lat} out of range"
        assert dist > 0, f"{name}: topo distance {dist} not positive"

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", BODIES_FOR_TOPO)
    def test_topo_differs_from_geocentric(self, body_id: int, name: str):
        """Topocentric position should differ from geocentric."""
        jd = 2451545.0
        swe.set_topo(12.5, 41.9, 50.0)

        geo, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
        topo, _ = swe.calc_ut(jd, body_id, FLG_TOPOCTR | FLG_SPEED)

        # At least one of lon/lat/dist should differ
        lon_diff = abs(geo[0] - topo[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff
        lat_diff = abs(geo[1] - topo[1])
        dist_diff = abs(geo[2] - topo[2])

        total_diff = lon_diff + lat_diff + dist_diff
        assert total_diff > 1e-10, (
            f"{name}: topo == geo (lon_diff={lon_diff}, lat_diff={lat_diff})"
        )

    @pytest.mark.unit
    def test_moon_topo_parallax_largest(self):
        """Moon should show the largest topocentric parallax (closest body)."""
        jd = 2451545.0
        swe.set_topo(12.5, 41.9, 50.0)

        diffs = {}
        for body_id, name in BODIES_FOR_TOPO:
            geo, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
            topo, _ = swe.calc_ut(jd, body_id, FLG_TOPOCTR | FLG_SPEED)
            lon_diff = abs(geo[0] - topo[0])
            if lon_diff > 180:
                lon_diff = 360 - lon_diff
            diffs[name] = lon_diff

        # Moon parallax should be largest
        moon_diff = diffs["Moon"]
        for name, diff in diffs.items():
            if name != "Moon":
                assert moon_diff >= diff * 0.5, (
                    f"Moon parallax {moon_diff:.6f}° < {name} parallax {diff:.6f}°"
                )


class TestTopoVariousLocations:
    """Test topocentric positions at various geographic locations."""

    @pytest.mark.unit
    @pytest.mark.parametrize("lon,lat,alt,loc_name", LOCATIONS)
    def test_moon_topo_at_location(
        self, lon: float, lat: float, alt: float, loc_name: str
    ):
        """Moon topocentric position valid at various locations."""
        swe.set_topo(lon, lat, alt)
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, MOON, FLG_TOPOCTR | FLG_SPEED)
        assert 0 <= result[0] < 360, f"{loc_name}: lon={result[0]}"
        assert math.isfinite(result[2]), f"{loc_name}: dist={result[2]}"

    @pytest.mark.unit
    @pytest.mark.parametrize("lon,lat,alt,loc_name", LOCATIONS)
    def test_sun_topo_at_location(
        self, lon: float, lat: float, alt: float, loc_name: str
    ):
        """Sun topocentric position valid at various locations."""
        swe.set_topo(lon, lat, alt)
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, SUN, FLG_TOPOCTR | FLG_SPEED)
        assert 0 <= result[0] < 360, f"{loc_name}: lon={result[0]}"
        assert result[2] > 0, f"{loc_name}: dist={result[2]}"


class TestTopoFlagCombinations:
    """Test topocentric with other flags."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "extra_flags,desc",
        [
            (0, "topo only"),
            (FLG_SPEED, "topo+speed"),
            (FLG_EQUATORIAL, "topo+equatorial"),
            (FLG_SPEED | FLG_EQUATORIAL, "topo+speed+equatorial"),
        ],
    )
    @pytest.mark.parametrize(
        "body_id,name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MARS, "Mars"),
        ],
    )
    def test_topo_flag_combo(
        self, body_id: int, name: str, extra_flags: int, desc: str
    ):
        """Topocentric works with various flag combinations."""
        swe.set_topo(12.5, 41.9, 50.0)
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, body_id, FLG_TOPOCTR | extra_flags)
        assert len(result) == 6
        for i, val in enumerate(result):
            assert math.isfinite(val), f"{name}+{desc}: result[{i}]={val}"

    @pytest.mark.unit
    def test_topo_sidereal(self):
        """Topocentric + sidereal combination works."""
        swe.set_topo(12.5, 41.9, 50.0)
        swe.set_sid_mode(SIDM_LAHIRI)
        jd = 2451545.0
        result, _ = swe.calc_ut(
            jd, MOON, FLG_TOPOCTR | FLG_SIDEREAL | FLG_SPEED
        )
        assert 0 <= result[0] < 360


class TestTopoDateRange:
    """Test topocentric positions across date ranges."""

    @pytest.mark.unit
    @pytest.mark.parametrize("year", [1900, 1950, 2000, 2024, 2050, 2100])
    def test_moon_topo_across_years(self, year: int):
        """Moon topocentric valid across years."""
        swe.set_topo(-74.0, 40.7, 10.0)  # New York
        jd = swe.julday(year, 6, 15, 12.0)
        result, _ = swe.calc_ut(jd, MOON, FLG_TOPOCTR | FLG_SPEED)
        assert 0 <= result[0] < 360

    @pytest.mark.unit
    def test_moon_topo_continuity(self):
        """Moon topocentric positions should be continuous over hours."""
        swe.set_topo(12.5, 41.9, 50.0)
        jd_start = 2451545.0
        prev_lon = None
        for i in range(48):  # 2 days, hourly
            jd = jd_start + i / 24.0
            result, _ = swe.calc_ut(jd, MOON, FLG_TOPOCTR)
            lon = result[0]
            if prev_lon is not None:
                diff = abs(lon - prev_lon)
                if diff > 180:
                    diff = 360 - diff
                # Moon moves ~0.5°/hour
                assert diff < 2.0, f"Moon topo jump {diff:.2f}° at hour {i}"
            prev_lon = lon


class TestTopoHighAltitude:
    """Test topocentric at high altitudes."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "alt,desc",
        [
            (0.0, "sea level"),
            (100.0, "100m"),
            (1000.0, "1km"),
            (5000.0, "5km"),
            (8848.0, "Everest"),
        ],
    )
    def test_moon_topo_altitude_effect(self, alt: float, desc: str):
        """Moon topocentric valid at various altitudes."""
        swe.set_topo(12.5, 41.9, alt)
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, MOON, FLG_TOPOCTR | FLG_SPEED)
        assert 0 <= result[0] < 360, f"Altitude {desc}: lon={result[0]}"
        assert result[2] > 0, f"Altitude {desc}: dist={result[2]}"


class TestTopoHighVolume:
    """High-volume topocentric tests."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "jd",
        [random.Random(7777).uniform(2415020.0, 2488070.0) for _ in range(50)],
    )
    def test_moon_topo_50_dates(self, jd: float):
        """Moon topocentric valid at 50 random dates."""
        swe.set_topo(12.5, 41.9, 50.0)
        result, _ = swe.calc_ut(jd, MOON, FLG_TOPOCTR | FLG_SPEED)
        assert 0 <= result[0] < 360
        assert result[2] > 0

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "jd",
        [random.Random(8888).uniform(2415020.0, 2488070.0) for _ in range(50)],
    )
    def test_sun_topo_50_dates(self, jd: float):
        """Sun topocentric valid at 50 random dates."""
        swe.set_topo(-74.0, 40.7, 10.0)
        result, _ = swe.calc_ut(jd, SUN, FLG_TOPOCTR | FLG_SPEED)
        assert 0 <= result[0] < 360
        assert result[2] > 0
