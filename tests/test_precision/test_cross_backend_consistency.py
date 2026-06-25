"""
Cross-backend consistency tests.

Verifies that Skyfield and LEB backends produce consistent results
for all core bodies across multiple dates and flag combinations.
"""

from __future__ import annotations

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
    URANUS,
    NEPTUNE,
    PLUTO,
    MEAN_NODE,
    TRUE_NODE,
    MEAN_APOG,
    FLG_SPEED,
    FLG_EQUATORIAL,
    FLG_SIDEREAL,
    FLG_J2000,
    FLG_NOABERR,
    SIDM_LAHIRI,
)


def _random_jds(n: int, seed: int = 42) -> list[float]:
    """Generate n random JDs in 1900-2100."""
    rng = random.Random(seed)
    return [rng.uniform(2415020.0, 2488070.0) for _ in range(n)]


def _angle_diff(a: float, b: float) -> float:
    """Absolute angular difference handling wrap-around."""
    d = abs(a - b) % 360
    return min(d, 360 - d)


# Bodies supported by both Skyfield and LEB
LEB_BODIES = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MERCURY, "Mercury"),
    (VENUS, "Venus"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
    (URANUS, "Uranus"),
    (NEPTUNE, "Neptune"),
    (PLUTO, "Pluto"),
    (MEAN_NODE, "MeanNode"),
    (TRUE_NODE, "TrueNode"),
    (MEAN_APOG, "MeanApog"),
]


class TestSkyfieldVsLEB:
    """Compare Skyfield backend against LEB backend."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,body_name", LEB_BODIES)
    def test_skyfield_leb_agree_at_j2000(self, body_id: int, body_name: str):
        """Skyfield and LEB agree at J2000 within 1 arcsecond."""
        jd = 2451545.0

        swe.set_calc_mode("skyfield")
        r_sky, _ = swe.calc_ut(jd, body_id, FLG_SPEED)

        swe.set_calc_mode("leb")
        r_leb, _ = swe.calc_ut(jd, body_id, FLG_SPEED)

        # Longitude agreement
        lon_diff = _angle_diff(r_sky[0], r_leb[0])
        tol = 1 / 3600  # 1 arcsecond
        assert lon_diff < tol, (
            f"{body_name}: Skyfield lon {r_sky[0]:.6f} vs LEB {r_leb[0]:.6f} "
            f'diff {lon_diff * 3600:.3f}"'
        )

        # Latitude agreement
        lat_diff = abs(r_sky[1] - r_leb[1])
        assert lat_diff < tol, (
            f"{body_name}: Skyfield lat {r_sky[1]:.6f} vs LEB {r_leb[1]:.6f} "
            f'diff {lat_diff * 3600:.3f}"'
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
        ],
    )
    def test_skyfield_leb_agree_20_dates(self, body_id: int, body_name: str):
        """Skyfield and LEB agree at 20 random dates."""
        jds = _random_jds(20, seed=body_id * 41)
        tol = 2 / 3600  # 2 arcseconds

        for jd in jds:
            swe.set_calc_mode("skyfield")
            r_sky, _ = swe.calc_ut(jd, body_id, FLG_SPEED)

            swe.set_calc_mode("leb")
            r_leb, _ = swe.calc_ut(jd, body_id, FLG_SPEED)

            lon_diff = _angle_diff(r_sky[0], r_leb[0])
            assert lon_diff < tol, (
                f'{body_name} @ JD {jd:.1f}: lon diff {lon_diff * 3600:.3f}"'
            )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MARS, "Mars"),
        ],
    )
    def test_skyfield_leb_distance_agree(self, body_id: int, body_name: str):
        """Skyfield and LEB distances agree within 1e-7 AU."""
        jd = 2451545.0

        swe.set_calc_mode("skyfield")
        r_sky, _ = swe.calc_ut(jd, body_id, 0)

        swe.set_calc_mode("leb")
        r_leb, _ = swe.calc_ut(jd, body_id, 0)

        dist_diff = abs(r_sky[2] - r_leb[2])
        assert dist_diff < 1e-7, f"{body_name}: distance diff {dist_diff} AU"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MARS, "Mars"),
        ],
    )
    def test_skyfield_leb_speed_agree(self, body_id: int, body_name: str):
        """Skyfield and LEB speeds agree within 1e-4 °/day."""
        jd = 2451545.0

        swe.set_calc_mode("skyfield")
        r_sky, _ = swe.calc_ut(jd, body_id, FLG_SPEED)

        swe.set_calc_mode("leb")
        r_leb, _ = swe.calc_ut(jd, body_id, FLG_SPEED)

        speed_diff = abs(r_sky[3] - r_leb[3])
        # Moon speed can differ slightly more between backends (~1e-3)
        tol = 1e-3 if body_id == MOON else 1e-4
        assert speed_diff < tol, f"{body_name}: speed diff {speed_diff}°/day"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "flags,desc",
        [
            (FLG_EQUATORIAL, "equatorial"),
            (FLG_J2000, "J2000"),
            (FLG_NOABERR, "no aberration"),
        ],
    )
    def test_skyfield_leb_flag_variants(self, flags: int, desc: str):
        """Skyfield and LEB agree with different flag variants."""
        jd = 2451545.0
        tol = 2 / 3600  # 2"

        swe.set_calc_mode("skyfield")
        r_sky, _ = swe.calc_ut(jd, MARS, flags)

        swe.set_calc_mode("leb")
        r_leb, _ = swe.calc_ut(jd, MARS, flags)

        lon_diff = _angle_diff(r_sky[0], r_leb[0])
        assert lon_diff < tol, f'Mars {desc}: lon diff {lon_diff * 3600:.3f}"'

    @pytest.mark.unit
    def test_skyfield_leb_sidereal_agree(self):
        """Skyfield and LEB agree in sidereal mode."""
        swe.set_sid_mode(SIDM_LAHIRI)
        jd = 2451545.0
        tol = 2 / 3600

        swe.set_calc_mode("skyfield")
        r_sky, _ = swe.calc_ut(jd, SUN, FLG_SIDEREAL)

        swe.set_calc_mode("leb")
        r_leb, _ = swe.calc_ut(jd, SUN, FLG_SIDEREAL)

        lon_diff = _angle_diff(r_sky[0], r_leb[0])
        assert lon_diff < tol, f'Sun sidereal: lon diff {lon_diff * 3600:.3f}"'
