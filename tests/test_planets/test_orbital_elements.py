"""
Tests for orbital elements calculation (get_orbital_elements).

This function calculates Keplerian orbital elements:
- Semi-major axis (a)
- Eccentricity (e)
- Inclination (i)
- Longitude of ascending node (Omega)
- Argument of perihelion (omega)
- Mean anomaly (M)
- And other derived parameters

The function returns a 50-element tuple matching pyswisseph's format.
"""

import pytest
import libephemeris as ephem
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
)


class TestOrbitalElementsBasic:
    """Test basic functionality of get_orbital_elements."""

    @pytest.mark.unit
    def test_returns_50_elements(self):
        """get_orbital_elements should return a tuple of 50 elements."""
        jd = 2451545.0  # J2000
        elements = ephem.get_orbital_elements(jd, MARS, 0)

        assert isinstance(elements, tuple)
        assert len(elements) == 50
        assert isinstance(elements[0], float)

    @pytest.mark.unit
    def test_returns_50_elements_ut(self):
        """get_orbital_elements_ut should return a tuple of 50 elements."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements_ut(jd, MARS, 0)

        assert isinstance(elements, tuple)
        assert len(elements) == 50

    @pytest.mark.unit
    def test_all_elements_are_floats(self):
        """All elements should be floats."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements(jd, JUPITER, 0)

        for i, val in enumerate(elements):
            assert isinstance(val, (int, float)), f"Element {i} should be numeric"

    @pytest.mark.unit
    def test_sun_returns_zeros(self):
        """Sun should return zero elements (no heliocentric orbit)."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements(jd, SUN, 0)

        # All elements should be zero for Sun
        assert elements[0] == 0.0  # Semi-major axis
        assert elements[1] == 0.0  # Eccentricity


class TestOrbitalElementsValues:
    """Test that orbital elements match known values."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name,expected_a,tolerance",
        [
            (MERCURY, "Mercury", 0.387, 0.01),
            (VENUS, "Venus", 0.723, 0.01),
            (MARS, "Mars", 1.524, 0.01),
            (JUPITER, "Jupiter", 5.203, 0.05),
            (SATURN, "Saturn", 9.537, 0.1),
            (URANUS, "Uranus", 19.19, 0.2),
            (NEPTUNE, "Neptune", 30.07, 0.3),
            (PLUTO, "Pluto", 39.48, 0.5),
        ],
    )
    def test_semi_major_axis(self, planet_id, planet_name, expected_a, tolerance):
        """Semi-major axis should match known values."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements(jd, planet_id, 0)

        a = elements[0]
        assert abs(a - expected_a) < tolerance, (
            f"{planet_name} semi-major axis {a:.4f} differs from expected {expected_a}"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name,expected_e,tolerance",
        [
            (MERCURY, "Mercury", 0.2056, 0.01),
            (VENUS, "Venus", 0.0068, 0.005),
            (MARS, "Mars", 0.0934, 0.01),
            (JUPITER, "Jupiter", 0.0484, 0.01),
            (SATURN, "Saturn", 0.0539, 0.01),
            (PLUTO, "Pluto", 0.2488, 0.02),
        ],
    )
    def test_eccentricity(self, planet_id, planet_name, expected_e, tolerance):
        """Eccentricity should match known values."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements(jd, planet_id, 0)

        e = elements[1]
        assert abs(e - expected_e) < tolerance, (
            f"{planet_name} eccentricity {e:.4f} differs from expected {expected_e}"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name,expected_i,tolerance",
        [
            (MERCURY, "Mercury", 7.0, 0.5),
            (VENUS, "Venus", 3.4, 0.5),
            (MARS, "Mars", 1.85, 0.2),
            (JUPITER, "Jupiter", 1.3, 0.2),
            (SATURN, "Saturn", 2.5, 0.3),
            (PLUTO, "Pluto", 17.1, 1.0),
        ],
    )
    def test_inclination(self, planet_id, planet_name, expected_i, tolerance):
        """Inclination should match known values."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements(jd, planet_id, 0)

        i = elements[2]
        assert abs(i - expected_i) < tolerance, (
            f"{planet_name} inclination {i:.4f}° differs from expected {expected_i}°"
        )


class TestOrbitalElementsRelationships:
    """Test mathematical relationships between elements."""

    @pytest.mark.unit
    def test_perihelion_from_a_and_e(self):
        """Perihelion distance should equal a*(1-e)."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements(jd, MARS, 0)

        a = elements[0]
        e = elements[1]
        q = elements[15]  # Perihelion distance

        expected_q = a * (1 - e)
        assert abs(q - expected_q) < 0.001, (
            f"Perihelion {q:.4f} should equal a*(1-e) = {expected_q:.4f}"
        )

    @pytest.mark.unit
    def test_aphelion_from_a_and_e(self):
        """Aphelion distance should equal a*(1+e)."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements(jd, MARS, 0)

        a = elements[0]
        e = elements[1]
        Q = elements[16]  # Aphelion distance

        expected_Q = a * (1 + e)
        assert abs(Q - expected_Q) < 0.001, (
            f"Aphelion {Q:.4f} should equal a*(1+e) = {expected_Q:.4f}"
        )

    @pytest.mark.unit
    def test_perihelion_less_than_aphelion(self):
        """Perihelion should always be less than aphelion."""
        jd = 2451545.0

        for planet_id in [MERCURY, MARS, JUPITER, SATURN]:
            elements = ephem.get_orbital_elements(jd, planet_id, 0)
            q = elements[15]
            Q = elements[16]
            assert q < Q, f"Perihelion {q} should be less than aphelion {Q}"

    @pytest.mark.unit
    def test_longitude_of_perihelion(self):
        """Longitude of perihelion should equal Omega + omega."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements(jd, MARS, 0)

        Omega = elements[3]  # Longitude of ascending node
        omega = elements[4]  # Argument of perihelion
        varpi = elements[5]  # Longitude of perihelion

        expected_varpi = (Omega + omega) % 360.0
        diff = abs(varpi - expected_varpi)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.01, (
            f"varpi {varpi:.4f} should equal Omega+omega = {expected_varpi:.4f}"
        )

    @pytest.mark.unit
    def test_mean_longitude(self):
        """Mean longitude should equal Omega + omega + M."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements(jd, JUPITER, 0)

        Omega = elements[3]
        omega = elements[4]
        M = elements[6]  # Mean anomaly
        L = elements[9]  # Mean longitude

        expected_L = (Omega + omega + M) % 360.0
        diff = abs(L - expected_L)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.01, f"L {L:.4f} should equal Omega+omega+M = {expected_L:.4f}"


class TestOrbitalElementsAnglesRange:
    """Test that angles are in valid ranges."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id",
        [MERCURY, VENUS, MARS, JUPITER, SATURN],
    )
    def test_angles_in_valid_range(self, planet_id):
        """All angular elements should be in 0-360 range."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements(jd, planet_id, 0)

        # Indices of angular elements
        angle_indices = [2, 3, 4, 5, 6, 7, 8, 9]  # i, Omega, omega, M, nu, E, L, varpi
        angle_names = ["i", "Omega", "omega", "M", "nu", "E", "L", "varpi"]

        for idx, name in zip(angle_indices, angle_names):
            val = elements[idx]
            assert 0 <= val < 360, f"{name} ({val}) should be in [0, 360)"


class TestOrbitalElementsOrbitalPeriod:
    """Test orbital period calculations."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name,expected_period_years,tolerance",
        [
            (MERCURY, "Mercury", 0.241, 0.01),
            (VENUS, "Venus", 0.615, 0.02),
            (MARS, "Mars", 1.881, 0.05),
            (JUPITER, "Jupiter", 11.86, 0.2),
            (SATURN, "Saturn", 29.46, 0.5),
        ],
    )
    def test_orbital_period(
        self, planet_id, planet_name, expected_period_years, tolerance
    ):
        """Orbital period should match known values."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements(jd, planet_id, 0)

        P = elements[10]  # Sidereal orbital period in tropical years
        assert abs(P - expected_period_years) < tolerance, (
            f"{planet_name} period {P:.3f} years differs from "
            f"expected {expected_period_years}"
        )

    @pytest.mark.unit
    def test_mean_motion_from_period(self):
        """Mean daily motion should be consistent with orbital period."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements(jd, MARS, 0)

        n = elements[11]  # Mean daily motion (degrees/day)
        P = elements[10]  # Sidereal orbital period (years)

        # n should equal 360 / (P * 365.24219)
        expected_n = 360.0 / (P * 365.24219)
        assert abs(n - expected_n) < 0.001, (
            f"Mean motion {n:.6f} deg/day differs from expected {expected_n:.6f}"
        )


class TestOrbitalElementsMoon:
    """Test orbital elements for the Moon (geocentric orbit)."""

    @pytest.mark.unit
    def test_moon_has_geocentric_elements(self):
        """Moon should have valid geocentric orbital elements."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements(jd, MOON, 0)

        a = elements[0]  # Semi-major axis
        e = elements[1]  # Eccentricity
        i = elements[2]  # Inclination

        # Moon's semi-major axis is about 384,400 km = 0.00257 AU
        assert 0.002 < a < 0.003, f"Moon semi-major axis {a} should be ~0.00257 AU"

        # Moon's eccentricity is about 0.0549
        assert 0.04 < e < 0.07, f"Moon eccentricity {e} should be ~0.0549"

        # Moon's inclination to ecliptic is about 5.145°
        assert 4.0 < i < 6.0, f"Moon inclination {i} should be ~5.145°"

    @pytest.mark.unit
    def test_moon_orbital_period(self):
        """Moon's orbital period should be about 27.3 days."""
        jd = 2451545.0
        elements = ephem.get_orbital_elements(jd, MOON, 0)

        P_years = elements[10]  # Sidereal orbital period
        P_days = P_years * 365.24219

        # Sidereal month is about 27.32 days
        assert 25 < P_days < 30, (
            f"Moon orbital period {P_days:.2f} days should be ~27.3"
        )


class TestOrbitalElementsCurrentDistance:
    """Test current heliocentric distance."""

    @pytest.mark.unit
    def test_current_distance_positive(self):
        """Current heliocentric distance should be positive."""
        jd = 2451545.0

        for planet_id in [MARS, JUPITER, SATURN]:
            elements = ephem.get_orbital_elements(jd, planet_id, 0)
            r = elements[16]  # Current distance
            assert r > 0, f"Current distance should be positive, got {r}"

    @pytest.mark.unit
    def test_current_distance_between_q_and_Q(self):
        """Current distance should be between perihelion and aphelion."""
        jd = 2451545.0

        for planet_id in [MARS, JUPITER, SATURN]:
            elements = ephem.get_orbital_elements(jd, planet_id, 0)
            q = elements[11]
            Q = elements[12]
            r = elements[16]

            assert q <= r <= Q, (
                f"Distance {r} should be between perihelion {q} and aphelion {Q}"
            )


class TestOrbitalElementsAliases:
    """Test that function aliases work correctly."""

    @pytest.mark.unit
    def test_swe_get_orbital_elements_alias(self):
        """get_orbital_elements should be available."""
        result = ephem.get_orbital_elements(2451545.0, MARS, 0)
        assert len(result) == 50
        assert isinstance(result[0], float)

    @pytest.mark.unit
    def test_swe_get_orbital_elements_ut_alias(self):
        """get_orbital_elements_ut should be available."""
        result = ephem.get_orbital_elements_ut(2451545.0, MARS, 0)
        assert len(result) == 50
        assert isinstance(result[0], float)

    @pytest.mark.unit
    def test_get_orbital_elements_alias(self):
        """get_orbital_elements should be available."""
        result = ephem.get_orbital_elements(2451545.0, MARS, 0)
        assert len(result) == 50
        assert isinstance(result[0], float)

    @pytest.mark.unit
    def test_get_orbital_elements_ut_alias(self):
        """get_orbital_elements_ut should be available."""
        result = ephem.get_orbital_elements_ut(2451545.0, MARS, 0)
        assert len(result) == 50
        assert isinstance(result[0], float)
