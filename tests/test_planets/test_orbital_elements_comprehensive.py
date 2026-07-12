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
    CHIRON,
    PHOLUS,
    CERES,
    PALLAS,
    JUNO,
    VESTA,
    AST_OFFSET,
    FLG_HELCTR,
    FLG_J2000,
    FLG_TRUEPOS,
    FLG_SPEED,
)
from libephemeris.exceptions import Error


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


# Known heliocentric osculating ranges for the curated minor bodies, valid
# across 1950-2050 and independent of the active backend (LEB/SPK/Skyfield).
# (body_id, name, a_lo, a_hi, e_lo, e_hi, i_lo, i_hi)
ASTEROID_RANGES = [
    (CERES, "Ceres", 2.75, 2.78, 0.07, 0.09, 10.4, 10.7),
    (CHIRON, "Chiron", 13.5, 13.8, 0.37, 0.39, 6.8, 7.1),
    (PHOLUS, "Pholus", 20.0, 20.5, 0.55, 0.59, 24.4, 25.0),
    (PALLAS, "Pallas", 2.76, 2.79, 0.22, 0.24, 34.6, 35.1),
    (JUNO, "Juno", 2.65, 2.68, 0.25, 0.27, 12.9, 13.1),
    (VESTA, "Vesta", 2.35, 2.37, 0.08, 0.10, 7.0, 7.2),
]


class TestOrbitalElementsAsteroids:
    """Osculating orbital elements for the curated minor bodies.

    Regression for finding I3-B: get_orbital_elements/orbit_max_min_true_distance
    used to return a silent all-zeros tuple for asteroids (reading as a real
    orbit at 0 AU). They now return the body's real osculating elements,
    computed from the same heliocentric state the position pipeline serves.
    """

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", [(b, n) for b, n, *_ in ASTEROID_RANGES])
    def test_asteroid_elements_nonzero(self, body_id: int, name: str):
        """Curated asteroids return non-zero elements (pre-fix: all zeros)."""
        jd = 2451545.0
        result = swe.get_orbital_elements(jd, body_id, FLG_HELCTR)
        assert result[0] != 0.0, f"{name}: semi-major axis silently zero"
        assert result[1] != 0.0, f"{name}: eccentricity silently zero"
        assert result[2] != 0.0, f"{name}: inclination silently zero"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,name,a_lo,a_hi,e_lo,e_hi,i_lo,i_hi", ASTEROID_RANGES
    )
    def test_asteroid_elements_physical_ranges(
        self,
        body_id: int,
        name: str,
        a_lo: float,
        a_hi: float,
        e_lo: float,
        e_hi: float,
        i_lo: float,
        i_hi: float,
    ):
        """Semi-major axis, eccentricity, inclination match known values."""
        jd = 2451545.0
        r = swe.get_orbital_elements(jd, body_id, FLG_HELCTR)
        assert a_lo < r[0] < a_hi, f"{name}: a = {r[0]:.4f} AU"
        assert e_lo < r[1] < e_hi, f"{name}: e = {r[1]:.4f}"
        assert i_lo < r[2] < i_hi, f"{name}: i = {r[2]:.4f} deg"

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", [(b, n) for b, n, *_ in ASTEROID_RANGES])
    def test_asteroid_elements_native_floats_finite(self, body_id: int, name: str):
        """All 50 slots are native Python floats and finite (no NaN/inf)."""
        jd = 2451545.0
        result = swe.get_orbital_elements(jd, body_id, FLG_HELCTR)
        assert len(result) == 50
        for idx, val in enumerate(result):
            assert type(val) is float, f"{name}[{idx}] is {type(val).__name__}"
            assert math.isfinite(val), f"{name}[{idx}] = {val}"

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", [(b, n) for b, n, *_ in ASTEROID_RANGES])
    def test_asteroid_perihelion_aphelion_consistency(self, body_id: int, name: str):
        """q = a(1-e) and Q = a(1+e), with q < a < Q."""
        jd = 2451545.0
        r = swe.get_orbital_elements(jd, body_id, FLG_HELCTR)
        a, e, q, Q = r[0], r[1], r[15], r[16]
        assert q < a < Q, f"{name}: not q < a < Q ({q}, {a}, {Q})"
        assert abs(q - a * (1 - e)) < 1e-6, f"{name}: q != a(1-e)"
        assert abs(Q - a * (1 + e)) < 1e-6, f"{name}: Q != a(1+e)"

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", [(b, n) for b, n, *_ in ASTEROID_RANGES])
    def test_asteroid_internal_coherence(self, body_id: int, name: str):
        """Reconstructing r from the reported elements at the epoch matches the
        heliocentric distance the position pipeline reports to < 1e-4 AU."""
        jd = 2451545.0
        r = swe.get_orbital_elements(jd, body_id, FLG_HELCTR)
        a, e, M = r[0], r[1], math.radians(r[6])
        # Solve Kepler's equation for the eccentric anomaly.
        E = M
        for _ in range(60):
            E -= (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
        r_recon = a * (1 - e * math.cos(E))
        pos, _ = swe.calc(jd, body_id, FLG_HELCTR | FLG_J2000 | FLG_TRUEPOS | FLG_SPEED)
        assert abs(r_recon - pos[2]) < 1e-4, (
            f"{name}: r_recon={r_recon:.8f} vs r_calc={pos[2]:.8f}"
        )

    @pytest.mark.unit
    def test_numbered_asteroid_elements(self):
        """A numbered asteroid (Hygiea, AST_OFFSET+10) yields real elements."""
        jd = 2451545.0
        try:
            r = swe.get_orbital_elements(jd, AST_OFFSET + 10, FLG_HELCTR)
        except Error as exc:  # pragma: no cover - only if no state source
            pytest.skip(f"Hygiea heliocentric state unavailable here: {exc}")
        # Hygiea: a ~ 3.14 AU, e ~ 0.10-0.13, i ~ 3.8 deg.
        assert 3.10 < r[0] < 3.20, f"Hygiea a = {r[0]:.4f} AU"
        assert 0.09 < r[1] < 0.14, f"Hygiea e = {r[1]:.4f}"
        assert 3.5 < r[2] < 4.1, f"Hygiea i = {r[2]:.4f} deg"
        assert all(type(v) is float and math.isfinite(v) for v in r)


class TestOrbitMaxMinTrueDistanceAsteroids:
    """orbit_max_min_true_distance is a true 3D two-ellipse extremum.

    The geocentric max/min are the maximum and minimum distance between the
    body's osculating ellipse and the Earth-Moon barycenter's, so the current
    geocentric distance must fall inside the [min, max] bracket, and both
    extremes must respect the aphelion/perihelion bounds of the two orbits.
    """

    # EMB osculating perihelion/aphelion (AU); stable to a few 1e-4 over the era.
    _EMB_Q = 1.0167
    _EMB_q = 0.9833

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", [(b, n) for b, n, *_ in ASTEROID_RANGES])
    def test_asteroid_max_min_brackets_true_distance(self, body_id: int, name: str):
        """Min/max form a valid 3D bracket around the true distance."""
        # Asteroid SPK safe range is ~1920-2080 CE; sample 1950/2000/2050.
        for jd in (2433282.5, 2451545.0, 2469807.5):
            max_d, min_d, true_d = swe.orbit_max_min_true_distance(jd, body_id, 0)
            assert all(
                type(v) is float and math.isfinite(v) for v in (max_d, min_d, true_d)
            ), f"{name}: non-float/non-finite result"
            assert min_d >= 0.0, f"{name}: min distance {min_d} < 0"
            assert min_d < max_d, f"{name}: min {min_d} >= max {max_d}"
            assert true_d > 0.0, f"{name}: true distance not positive"
            # The instantaneous geocentric distance lies inside the bracket
            # (small tolerance for the Earth-vs-EMB offset).
            assert min_d - 1e-3 <= true_d <= max_d + 1e-3, (
                f"{name}@{jd}: true {true_d} outside [{min_d}, {max_d}]"
            )

            # Physical bounds from the body's own osculating aphelion/perihelion.
            r = swe.get_orbital_elements(jd, body_id, FLG_HELCTR)
            a, e = r[0], r[1]
            q_body, big_q_body = a * (1 - e), a * (1 + e)
            # Max distance never exceeds aphelion-to-aphelion, never below the
            # body's aphelion minus Earth's perihelion.
            assert (
                big_q_body - self._EMB_q - 1e-6
                <= max_d
                <= big_q_body + self._EMB_Q + 1e-6
            ), f"{name}@{jd}: max {max_d} outside aphelion bounds"
            # Min distance never below perihelion-difference, never above the
            # body's perihelion plus Earth's aphelion.
            assert min_d <= q_body + self._EMB_Q + 1e-6, (
                f"{name}@{jd}: min {min_d} above perihelion+EMB aphelion"
            )
