"""
Comprehensive tests for orbital elements calculation.

Verifies get_orbital_elements and get_orbital_elements_ut
return physically plausible orbital parameters for planets and asteroids.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe
from libephemeris.constants import (
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    URANUS,
    NEPTUNE,
    PLUTO,
    FLG_HELCTR,
)


# Known approximate orbital elements (heliocentric, J2000)
PLANET_ELEMENTS = {
    # (body_id, name, semi_major_AU, eccentricity, inclination_deg, period_years)
    MERCURY: ("Mercury", 0.387, 0.206, 7.0, 0.241),
    VENUS: ("Venus", 0.723, 0.007, 3.4, 0.615),
    MARS: ("Mars", 1.524, 0.093, 1.85, 1.881),
    JUPITER: ("Jupiter", 5.203, 0.048, 1.3, 11.86),
    SATURN: ("Saturn", 9.537, 0.054, 2.49, 29.46),
    URANUS: ("Uranus", 19.19, 0.047, 0.77, 84.01),
    NEPTUNE: ("Neptune", 30.07, 0.009, 1.77, 164.8),
    PLUTO: ("Pluto", 39.48, 0.249, 17.2, 248.1),
}


class TestOrbitalElementsBasic:
    """Basic orbital elements tests."""

    @pytest.mark.unit
    def test_returns_50_elements(self):
        """get_orbital_elements returns 50 floats."""
        jd = 2451545.0
        result = swe.get_orbital_elements(jd, MARS, FLG_HELCTR)
        assert len(result) == 50, f"Expected 50, got {len(result)}"

    @pytest.mark.unit
    def test_returns_native_floats(self):
        """All elements should be native Python float."""
        jd = 2451545.0
        result = swe.get_orbital_elements(jd, MARS, FLG_HELCTR)
        for i, val in enumerate(result):
            assert type(val) is float, f"Element {i} is {type(val).__name__}"

    @pytest.mark.unit
    def test_all_finite(self):
        """All elements should be finite."""
        jd = 2451545.0
        result = swe.get_orbital_elements(jd, MARS, FLG_HELCTR)
        for i in range(18):  # First 18 are the real elements
            assert math.isfinite(result[i]), f"Element {i} = {result[i]}"

    @pytest.mark.unit
    def test_ut_variant_returns_50(self):
        """get_orbital_elements_ut also returns 50 elements."""
        jd = 2451545.0
        result = swe.get_orbital_elements_ut(jd, MARS, FLG_HELCTR)
        assert len(result) == 50


class TestOrbitalElementsSemiMajorAxis:
    """Test semi-major axis values."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id",
        [MERCURY, VENUS, MARS, JUPITER, SATURN],
    )
    def test_semi_major_axis_plausible(self, body_id: int):
        """Semi-major axis should be close to known values."""
        name, expected_a, _, _, _ = PLANET_ELEMENTS[body_id]
        jd = 2451545.0
        result = swe.get_orbital_elements(jd, body_id, FLG_HELCTR)
        a = result[0]
        # Within 5% of known value
        assert abs(a - expected_a) / expected_a < 0.05, (
            f"{name}: a = {a:.4f} AU, expected ~{expected_a} AU"
        )


class TestOrbitalElementsEccentricity:
    """Test eccentricity values."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id",
        [MERCURY, VENUS, MARS, JUPITER, SATURN],
    )
    def test_eccentricity_in_range(self, body_id: int):
        """Eccentricity should be 0 <= e < 1 for bound orbits."""
        name = PLANET_ELEMENTS[body_id][0]
        jd = 2451545.0
        result = swe.get_orbital_elements(jd, body_id, FLG_HELCTR)
        e = result[1]
        assert 0 <= e < 1, f"{name}: eccentricity {e} out of range"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id",
        [MERCURY, VENUS, MARS, JUPITER],
    )
    def test_eccentricity_approximate(self, body_id: int):
        """Eccentricity should be close to known value."""
        name, _, expected_e, _, _ = PLANET_ELEMENTS[body_id]
        jd = 2451545.0
        result = swe.get_orbital_elements(jd, body_id, FLG_HELCTR)
        e = result[1]
        assert abs(e - expected_e) < 0.02, (
            f"{name}: e = {e:.4f}, expected ~{expected_e}"
        )


class TestOrbitalElementsInclination:
    """Test inclination values."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id",
        [MERCURY, VENUS, MARS, JUPITER, SATURN],
    )
    def test_inclination_positive(self, body_id: int):
        """Inclination should be positive."""
        name = PLANET_ELEMENTS[body_id][0]
        jd = 2451545.0
        result = swe.get_orbital_elements(jd, body_id, FLG_HELCTR)
        i = result[2]
        assert 0 <= i < 180, f"{name}: inclination {i}°"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id",
        [MERCURY, VENUS, MARS],
    )
    def test_inclination_approximate(self, body_id: int):
        """Inclination should be close to known value."""
        name, _, _, expected_i, _ = PLANET_ELEMENTS[body_id]
        jd = 2451545.0
        result = swe.get_orbital_elements(jd, body_id, FLG_HELCTR)
        i = result[2]
        assert abs(i - expected_i) < 2.0, (
            f"{name}: i = {i:.2f}°, expected ~{expected_i}°"
        )


class TestOrbitalElementsPeriod:
    """Test orbital period values."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id",
        [MERCURY, VENUS, MARS, JUPITER, SATURN],
    )
    def test_period_positive(self, body_id: int):
        """Sidereal period should be positive."""
        name = PLANET_ELEMENTS[body_id][0]
        jd = 2451545.0
        result = swe.get_orbital_elements(jd, body_id, FLG_HELCTR)
        period = result[10]  # Sidereal period in tropical years
        assert period > 0, f"{name}: period {period} not positive"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id",
        [MERCURY, VENUS, MARS, JUPITER],
    )
    def test_period_approximate(self, body_id: int):
        """Period should be close to known value."""
        name, _, _, _, expected_p = PLANET_ELEMENTS[body_id]
        jd = 2451545.0
        result = swe.get_orbital_elements(jd, body_id, FLG_HELCTR)
        period = result[10]
        ratio = period / expected_p
        assert 0.9 < ratio < 1.1, (
            f"{name}: period = {period:.3f} yr, expected ~{expected_p} yr"
        )


class TestOrbitalElementsAngles:
    """Test angular orbital elements."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id",
        [MERCURY, MARS, JUPITER],
    )
    def test_node_in_range(self, body_id: int):
        """Longitude of ascending node should be 0-360°."""
        name = PLANET_ELEMENTS[body_id][0]
        jd = 2451545.0
        result = swe.get_orbital_elements(jd, body_id, FLG_HELCTR)
        node = result[3]
        assert 0 <= node < 360, f"{name}: node = {node}°"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id",
        [MERCURY, MARS, JUPITER],
    )
    def test_arg_perihelion_in_range(self, body_id: int):
        """Argument of perihelion should be 0-360°."""
        name = PLANET_ELEMENTS[body_id][0]
        jd = 2451545.0
        result = swe.get_orbital_elements(jd, body_id, FLG_HELCTR)
        omega = result[4]
        assert 0 <= omega < 360, f"{name}: omega = {omega}°"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id",
        [MERCURY, MARS, JUPITER],
    )
    def test_mean_anomaly_in_range(self, body_id: int):
        """Mean anomaly should be 0-360°."""
        name = PLANET_ELEMENTS[body_id][0]
        jd = 2451545.0
        result = swe.get_orbital_elements(jd, body_id, FLG_HELCTR)
        M = result[6]
        assert 0 <= M < 360, f"{name}: M = {M}°"


class TestOrbitalElementsConsistency:
    """Test consistency between orbital elements."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id",
        [MARS, JUPITER],
    )
    def test_perihelion_aphelion_distance(self, body_id: int):
        """q < a < Q for elliptical orbits."""
        name = PLANET_ELEMENTS[body_id][0]
        jd = 2451545.0
        result = swe.get_orbital_elements(jd, body_id, FLG_HELCTR)
        a = result[0]
        q = result[15]  # Perihelion distance
        Q = result[16]  # Aphelion distance
        assert q < a, f"{name}: q={q} >= a={a}"
        assert a < Q, f"{name}: a={a} >= Q={Q}"
        # q = a(1-e), Q = a(1+e)
        e = result[1]
        assert abs(q - a * (1 - e)) < 0.01, f"{name}: q != a(1-e)"
        assert abs(Q - a * (1 + e)) < 0.01, f"{name}: Q != a(1+e)"

    @pytest.mark.unit
    def test_mean_daily_motion_consistent(self):
        """Mean daily motion n should be ~360/P_days."""
        jd = 2451545.0
        result = swe.get_orbital_elements(jd, MARS, FLG_HELCTR)
        n = result[11]  # Mean daily motion (deg/day)
        P_years = result[10]  # Sidereal period (years)
        P_days = P_years * 365.25
        expected_n = 360.0 / P_days
        assert abs(n - expected_n) / expected_n < 0.01, (
            f"n={n:.6f}, expected {expected_n:.6f}"
        )


class TestOrbitalElementsDateRange:
    """Test orbital elements across dates."""

    @pytest.mark.unit
    @pytest.mark.parametrize("year", [1900, 2000, 2100])
    def test_mars_elements_across_centuries(self, year: int):
        """Mars orbital elements valid across centuries."""
        jd = swe.julday(year, 1, 1, 12.0)
        result = swe.get_orbital_elements(jd, MARS, FLG_HELCTR)
        a = result[0]
        e = result[1]
        # Mars semi-major axis doesn't change much
        assert 1.4 < a < 1.6, f"Year {year}: a={a}"
        assert 0.05 < e < 0.15, f"Year {year}: e={e}"


class TestOrbitalElementsMoon:
    """Test Moon orbital elements (geocentric)."""

    @pytest.mark.unit
    def test_moon_geocentric_elements(self):
        """Moon elements should have geocentric parameters."""
        jd = 2451545.0
        result = swe.get_orbital_elements(jd, MOON, 0)
        a = result[0]
        e = result[1]
        i = result[2]
        # Moon semi-major axis ~0.00257 AU (~384,400 km)
        assert 0.001 < a < 0.005, f"Moon a = {a} AU"
        # Moon eccentricity ~0.055
        assert 0.01 < e < 0.1, f"Moon e = {e}"
        # Moon inclination to ecliptic ~5.15°
        assert 3 < i < 8, f"Moon i = {i}°"
