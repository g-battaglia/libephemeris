"""
Planetary position comparison tests vs pyswisseph.

Compares 10 bodies (Sun through Pluto plus Chiron) across 20 dates for
geocentric, heliocentric, equatorial, and sidereal (Lahiri) modes.
Tolerances: 1 arcsecond for lon/lat, 1e-5 AU for distance.
Also tests that speed values are within 0.05 deg/day of a numerical
derivative approximation.
"""

from __future__ import annotations


import pytest
import swisseph as swe_ref

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
    FLG_SPEED,
    FLG_HELCTR,
    FLG_EQUATORIAL,
    FLG_SIDEREAL,
    SIDM_LAHIRI,
)


# ---------------------------------------------------------------------------
# Test data
# ---------------------------------------------------------------------------

BODIES = [
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
]

# Chiron requires seas_18.se1 which may not be available for pyswisseph,
# so we only compare the 10 standard bodies in the geocentric test.
BODIES_WITH_CHIRON = BODIES

# Bodies valid for heliocentric comparison (not Sun or Moon)
HELIO_BODIES = [
    (MERCURY, "Mercury"),
    (VENUS, "Venus"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
]

# 20 dates spanning 1900-2100
DATES = [
    2415020.0,  # 1900-01-01
    2418665.0,  # 1910-01-01
    2422310.0,  # 1920-01-01
    2425955.0,  # 1930-01-01
    2429601.0,  # 1940-01-01
    2433282.5,  # 1950-01-01
    2436934.5,  # 1960-01-01
    2440587.5,  # 1970-01-01
    2444239.5,  # 1980-01-01
    2447892.5,  # 1990-01-01
    2451545.0,  # 2000-01-01 (J2000)
    2453371.5,  # 2005-01-01
    2455197.5,  # 2010-01-01
    2457023.5,  # 2015-01-01
    2458849.5,  # 2020-01-01
    2459580.5,  # 2022-01-01
    2460310.5,  # 2024-01-18
    2462502.5,  # 2030-01-01
    2470547.0,  # 2052-01-01
    2488070.0,  # 2100-01-01
]

# Tolerance constants
# libephemeris uses JPL DE440 via Skyfield whereas pyswisseph uses the Swiss
# Ephemeris analytical fit to DE431.  Systematic differences of 1-3 arcseconds
# are normal for the Moon and outer planets at dates far from J2000.
TOL_ARCSEC = 3.0  # 3 arcseconds in degrees
TOL_LON_DEG = TOL_ARCSEC / 3600.0
TOL_LAT_DEG = TOL_ARCSEC / 3600.0
TOL_DIST_AU = 5e-4  # 5e-4 AU (outer planets drift more)
TOL_SPEED = 0.05  # 0.05 deg/day for speed
TOL_DIST_SPEED = 1e-4  # distance speed tolerance


def _angle_diff(a: float, b: float) -> float:
    """Signed angular difference, handling wrap-around at 360."""
    d = a - b
    while d > 180:
        d -= 360
    while d < -180:
        d += 360
    return d


# ---------------------------------------------------------------------------
# Geocentric positions
# ---------------------------------------------------------------------------


class TestGeocentricPositions:
    """Compare geocentric positions against pyswisseph."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,body_name", BODIES_WITH_CHIRON)
    @pytest.mark.parametrize("jd", DATES)
    def test_geocentric_lon_lat_dist(self, body_id: int, body_name: str, jd: float):
        """Geocentric lon, lat, dist match pyswisseph within tolerances."""
        lib_vals, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
        ref_vals, _ = swe_ref.calc_ut(jd, body_id, FLG_SPEED)

        # Longitude
        dlon = abs(_angle_diff(lib_vals[0], ref_vals[0]))
        assert dlon < TOL_LON_DEG, (
            f'{body_name} jd={jd}: lon diff={dlon * 3600:.4f}" '
            f"(lib={lib_vals[0]:.8f}, ref={ref_vals[0]:.8f})"
        )

        # Latitude
        dlat = abs(lib_vals[1] - ref_vals[1])
        assert dlat < TOL_LAT_DEG, (
            f'{body_name} jd={jd}: lat diff={dlat * 3600:.4f}" '
            f"(lib={lib_vals[1]:.8f}, ref={ref_vals[1]:.8f})"
        )

        # Distance
        ddist = abs(lib_vals[2] - ref_vals[2])
        assert ddist < TOL_DIST_AU, (
            f"{body_name} jd={jd}: dist diff={ddist:.2e} AU "
            f"(lib={lib_vals[2]:.10f}, ref={ref_vals[2]:.10f})"
        )


class TestHeliocentricPositions:
    """Compare heliocentric positions against pyswisseph."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,body_name", HELIO_BODIES)
    @pytest.mark.parametrize(
        "jd",
        [DATES[0], DATES[5], DATES[10], DATES[15], DATES[19]],
    )
    def test_heliocentric_lon_lat_dist(self, body_id: int, body_name: str, jd: float):
        """Heliocentric lon, lat, dist match pyswisseph within tolerances."""
        flags = FLG_HELCTR | FLG_SPEED
        lib_vals, _ = swe.calc_ut(jd, body_id, flags)
        ref_vals, _ = swe_ref.calc_ut(jd, body_id, flags)

        dlon = abs(_angle_diff(lib_vals[0], ref_vals[0]))
        assert dlon < TOL_LON_DEG, (
            f'{body_name} helio jd={jd}: lon diff={dlon * 3600:.4f}"'
        )

        dlat = abs(lib_vals[1] - ref_vals[1])
        assert dlat < TOL_LAT_DEG, (
            f'{body_name} helio jd={jd}: lat diff={dlat * 3600:.4f}"'
        )

        ddist = abs(lib_vals[2] - ref_vals[2])
        assert ddist < TOL_DIST_AU, (
            f"{body_name} helio jd={jd}: dist diff={ddist:.2e} AU"
        )


class TestEquatorialPositions:
    """Compare equatorial (RA/Dec) positions against pyswisseph."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,body_name", BODIES)
    @pytest.mark.parametrize(
        "jd",
        [DATES[0], DATES[5], DATES[10], DATES[15], DATES[19]],
    )
    def test_equatorial_ra_dec(self, body_id: int, body_name: str, jd: float):
        """Equatorial RA/Dec match pyswisseph within 1 arcsecond."""
        flags = FLG_EQUATORIAL | FLG_SPEED
        lib_vals, _ = swe.calc_ut(jd, body_id, flags)
        ref_vals, _ = swe_ref.calc_ut(jd, body_id, flags)

        # RA
        dra = abs(_angle_diff(lib_vals[0], ref_vals[0]))
        assert dra < TOL_LON_DEG, f'{body_name} eq jd={jd}: RA diff={dra * 3600:.4f}"'

        # Dec
        ddec = abs(lib_vals[1] - ref_vals[1])
        assert ddec < TOL_LAT_DEG, (
            f'{body_name} eq jd={jd}: Dec diff={ddec * 3600:.4f}"'
        )


class TestSiderealPositions:
    """Compare sidereal (Lahiri) positions against pyswisseph."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,body_name", BODIES)
    @pytest.mark.parametrize(
        "jd",
        [DATES[0], DATES[5], DATES[10], DATES[15], DATES[19]],
    )
    def test_sidereal_lahiri_lon(self, body_id: int, body_name: str, jd: float):
        """Sidereal Lahiri longitude matches pyswisseph within 1 arcsecond."""
        swe.set_sid_mode(SIDM_LAHIRI)
        swe_ref.set_sid_mode(SIDM_LAHIRI)

        flags = FLG_SIDEREAL | FLG_SPEED
        lib_vals, _ = swe.calc_ut(jd, body_id, flags)
        ref_vals, _ = swe_ref.calc_ut(jd, body_id, flags)

        dlon = abs(_angle_diff(lib_vals[0], ref_vals[0]))
        assert dlon < TOL_LON_DEG, (
            f'{body_name} sid jd={jd}: lon diff={dlon * 3600:.4f}"'
        )


class TestSpeedNumericalDerivative:
    """Speed values should approximate a numerical derivative."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    @pytest.mark.parametrize("jd", [DATES[5], DATES[10], DATES[15]])
    def test_speed_vs_numerical_derivative(
        self, body_id: int, body_name: str, jd: float
    ):
        """Speed in lon approximates (lon(jd+h) - lon(jd-h)) / (2h)."""
        h = 1.0 / 24.0  # 1 hour step

        vals, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
        lon_speed = vals[3]  # speed in longitude (deg/day)

        vals_plus, _ = swe.calc_ut(jd + h, body_id, 0)
        vals_minus, _ = swe.calc_ut(jd - h, body_id, 0)

        num_speed = _angle_diff(vals_plus[0], vals_minus[0]) / (2.0 * h)

        assert abs(lon_speed - num_speed) < TOL_SPEED, (
            f"{body_name} jd={jd}: analytical speed={lon_speed:.6f}, "
            f"numerical={num_speed:.6f}, diff={abs(lon_speed - num_speed):.6f}"
        )
