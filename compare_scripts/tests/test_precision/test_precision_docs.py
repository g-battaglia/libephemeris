"""
Tests for precision documentation (PRECISION.md).

These tests verify that the precision claims documented in docs/PRECISION.md
are accurate and hold true for libephemeris calculations.
"""

import pytest
import random
import swisseph as swe
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
    MEAN_NODE,
    TRUE_NODE,
    MEAN_APOG,
    OSCU_APOG,
    CHIRON,
    FLG_SPEED,
    FLG_HELCTR,
    SIDM_LAHIRI,
    SIDM_FAGAN_BRADLEY,
    SIDM_TRUE_CITRA,
    SIDM_GALCENT_0SAG,
)
from libephemeris.crossing import (
    NR_TOLERANCE,
    NR_TOLERANCE_SUN,
    NR_TOLERANCE_MOON,
)


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================


def angle_diff(a1: float, a2: float) -> float:
    """Calculate smallest difference between two angles in degrees."""
    diff = abs(a1 - a2)
    if diff > 180:
        diff = 360 - diff
    return diff


def generate_random_jd(n: int = 100, seed: int = 42) -> list:
    """Generate n random Julian Days in DE421 range."""
    random.seed(seed)
    jds = []
    for _ in range(n):
        year = random.randint(1900, 2050)
        month = random.randint(1, 12)
        day = random.randint(1, 28)
        hour = random.uniform(0, 24)
        jd = ephem.julday(year, month, day, hour)
        jds.append(jd)
    return jds


# =============================================================================
# DOCUMENTED PRECISION CONSTANTS
# =============================================================================


# From PRECISION.md - Planetary longitude tolerances in arcseconds
PLANET_TOLERANCE_ARCSEC = {
    SUN: 1.0,
    MOON: 5.0,
    MERCURY: 1.0,
    VENUS: 1.0,
    MARS: 2.0,
    JUPITER: 5.0,
    SATURN: 5.0,
    URANUS: 5.0,
    NEPTUNE: 5.0,
    PLUTO: 5.0,
}

# From PRECISION.md - House cusp tolerance
HOUSE_CUSP_TOLERANCE_DEG = 0.001  # ~3.6 arcsec

# From PRECISION.md - Ayanamsha tolerances
AYANAMSHA_STANDARD_TOLERANCE_DEG = 0.01
AYANAMSHA_STAR_BASED_TOLERANCE_DEG = 0.06

# From PRECISION.md - Crossing tolerances
CROSSING_SUN_TOLERANCE_ARCSEC = 0.001
CROSSING_MOON_TOLERANCE_ARCSEC = 0.001
CROSSING_PLANET_TOLERANCE_ARCSEC = 0.001

# From PRECISION.md - Lunar points tolerances (vary by point type)
LUNAR_POINT_TOLERANCE_DEG = {
    MEAN_NODE: 0.01,  # High precision
    TRUE_NODE: 2.0,  # Different oscillation model
    MEAN_APOG: 0.2,  # Minor formula differences
    OSCU_APOG: 0.1,  # Eccentricity vector method (~235 arcsec max)
}

# From PRECISION.md - Julian Day tolerance
JULIAN_DAY_TOLERANCE = 1e-10


# =============================================================================
# TEST CLASSES
# =============================================================================


class TestDocumentedPlanetaryPrecision:
    """Verify documented planetary position precision claims."""

    @pytest.mark.precision
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_documented_planet_longitude_tolerance(self, planet_id, planet_name):
        """
        Verify planet longitude matches pyswisseph within documented tolerance.

        Tests 100 random dates and verifies max difference is within the
        tolerance documented in PRECISION.md.
        """
        tolerance_arcsec = PLANET_TOLERANCE_ARCSEC[planet_id]
        jds = generate_random_jd(100)

        max_diff_arcsec = 0.0
        for jd in jds:
            pos_lib, _ = ephem.calc_ut(jd, planet_id, 0)
            pos_swe, _ = swe.calc_ut(jd, planet_id, 0)

            diff = angle_diff(pos_lib[0], pos_swe[0])
            diff_arcsec = diff * 3600
            max_diff_arcsec = max(max_diff_arcsec, diff_arcsec)

        assert max_diff_arcsec < tolerance_arcsec, (
            f"{planet_name} max diff {max_diff_arcsec:.4f} arcsec "
            f">= documented tolerance {tolerance_arcsec} arcsec"
        )

    @pytest.mark.precision
    def test_documented_latitude_precision(self):
        """Verify latitude precision matches documented claims."""
        jd = 2451545.0  # J2000

        # Sun latitude should be < 1 arcsec
        sun_lib, _ = ephem.calc_ut(jd, SUN, 0)
        sun_swe, _ = swe.calc_ut(jd, SUN, 0)
        sun_lat_diff = abs(sun_lib[1] - sun_swe[1]) * 3600
        assert sun_lat_diff < 1.0, f"Sun latitude diff {sun_lat_diff:.4f} arcsec"

        # Moon latitude should be < 5 arcsec
        moon_lib, _ = ephem.calc_ut(jd, MOON, 0)
        moon_swe, _ = swe.calc_ut(jd, MOON, 0)
        moon_lat_diff = abs(moon_lib[1] - moon_swe[1]) * 3600
        assert moon_lat_diff < 5.0, f"Moon latitude diff {moon_lat_diff:.4f} arcsec"

    @pytest.mark.precision
    def test_documented_distance_precision(self):
        """Verify distance precision < 0.0001 AU as documented."""
        jd = 2451545.0

        for planet_id in [SUN, MERCURY, VENUS, MARS, JUPITER]:
            pos_lib, _ = ephem.calc_ut(jd, planet_id, 0)
            pos_swe, _ = swe.calc_ut(jd, planet_id, 0)

            dist_diff = abs(pos_lib[2] - pos_swe[2])
            assert dist_diff < 0.0001, (
                f"Planet {planet_id} distance diff {dist_diff:.6f} AU >= 0.0001 AU"
            )


class TestDocumentedHeliocentricPrecision:
    """Verify documented heliocentric precision claims."""

    @pytest.mark.precision
    def test_heliocentric_relaxed_tolerance(self):
        """Verify heliocentric tolerance is <= 0.03 degrees as documented."""
        jd = 2451545.0

        for planet_id in [MERCURY, VENUS, MARS, JUPITER]:
            pos_lib, _ = ephem.calc_ut(jd, planet_id, FLG_HELCTR)
            pos_swe, _ = swe.calc_ut(jd, planet_id, FLG_HELCTR)

            diff = angle_diff(pos_lib[0], pos_swe[0])
            assert diff < 0.03, (
                f"Planet {planet_id} heliocentric diff {diff:.6f} degrees >= 0.03"
            )


class TestDocumentedHousePrecision:
    """Verify documented house system precision claims."""

    @pytest.mark.precision
    @pytest.mark.parametrize(
        "hsys,name",
        [
            (ord("P"), "Placidus"),
            (ord("K"), "Koch"),
            (ord("O"), "Porphyry"),
            (ord("E"), "Equal"),
            (ord("W"), "Whole Sign"),
        ],
    )
    def test_house_cusp_tolerance(self, hsys, name):
        """Verify house cusp tolerance <= 0.001 degrees as documented."""
        jd = 2451545.0
        lat, lon = 41.9, 12.5  # Rome

        cusps_lib, ascmc_lib = ephem.houses(jd, lat, lon, hsys)
        cusps_swe, ascmc_swe = swe.houses(jd, lat, lon, bytes([hsys]))

        # Check all 12 cusps (both are 0-indexed)
        for i in range(12):
            diff = angle_diff(cusps_lib[i], cusps_swe[i])
            assert diff < HOUSE_CUSP_TOLERANCE_DEG, (
                f"{name} cusp {i + 1} diff {diff:.6f} degrees >= 0.001"
            )

        # Check ASC and MC
        asc_diff = angle_diff(ascmc_lib[0], ascmc_swe[0])
        mc_diff = angle_diff(ascmc_lib[1], ascmc_swe[1])

        assert asc_diff < HOUSE_CUSP_TOLERANCE_DEG, (
            f"{name} ASC diff {asc_diff:.6f} degrees >= 0.001"
        )
        assert mc_diff < HOUSE_CUSP_TOLERANCE_DEG, (
            f"{name} MC diff {mc_diff:.6f} degrees >= 0.001"
        )

    @pytest.mark.precision
    def test_polar_latitude_fallback(self):
        """
        Verify polar latitude behavior as documented.

        At latitudes > 66.5 degrees, Placidus should fall back to Porphyry
        or return valid values without error.
        """
        jd = 2451545.0
        polar_lat = 70.0  # Above Arctic Circle
        lon = 0.0

        # Should not raise an error - falls back gracefully
        try:
            cusps_lib, ascmc_lib = ephem.houses(jd, polar_lat, lon, ord("P"))
            # Values should be returned (fallback to Porphyry)
            assert len(cusps_lib) == 12
            assert len(ascmc_lib) >= 2
            # ASC should be valid (between 0 and 360)
            assert 0 <= ascmc_lib[0] < 360
        except Exception as e:
            # Some implementations may raise an error, which is acceptable
            assert "polar" in str(e).lower() or "latitude" in str(e).lower()


class TestDocumentedAyanamshaPrecision:
    """Verify documented ayanamsha precision claims."""

    @pytest.mark.precision
    @pytest.mark.parametrize(
        "sid_mode,name",
        [
            (SIDM_FAGAN_BRADLEY, "Fagan-Bradley"),
            (SIDM_LAHIRI, "Lahiri"),
        ],
    )
    def test_standard_ayanamsha_tolerance(self, sid_mode, name):
        """Verify standard ayanamsha tolerance <= 0.01 degrees."""
        jd = 2451545.0

        ephem.set_sid_mode(sid_mode)
        swe.set_sid_mode(sid_mode)

        ayan_lib = ephem.get_ayanamsa_ut(jd)
        ayan_swe = swe.get_ayanamsa_ut(jd)

        diff = abs(ayan_lib - ayan_swe)
        assert diff < AYANAMSHA_STANDARD_TOLERANCE_DEG, (
            f"{name} ayanamsha diff {diff:.6f} degrees >= 0.01"
        )

    @pytest.mark.precision
    @pytest.mark.parametrize(
        "sid_mode,name",
        [
            (SIDM_TRUE_CITRA, "True Citra"),
            (SIDM_GALCENT_0SAG, "Galactic Center 0 Sag"),
        ],
    )
    def test_star_based_ayanamsha_tolerance(self, sid_mode, name):
        """Verify star-based ayanamsha tolerance <= 0.06 degrees."""
        jd = 2451545.0

        ephem.set_sid_mode(sid_mode)
        swe.set_sid_mode(sid_mode)

        ayan_lib = ephem.get_ayanamsa_ut(jd)
        ayan_swe = swe.get_ayanamsa_ut(jd)

        diff = abs(ayan_lib - ayan_swe)
        assert diff < AYANAMSHA_STAR_BASED_TOLERANCE_DEG, (
            f"{name} ayanamsha diff {diff:.6f} degrees >= 0.06"
        )


class TestDocumentedCrossingPrecision:
    """Verify documented crossing precision claims."""

    @pytest.mark.precision
    def test_crossing_tolerance_constants(self):
        """Verify crossing tolerance constants match documentation."""
        # NR_TOLERANCE for generic planets: 0.1 arcsec
        expected_planet_tol = CROSSING_PLANET_TOLERANCE_ARCSEC / 3600.0
        assert abs(NR_TOLERANCE - expected_planet_tol) < 1e-10, (
            f"NR_TOLERANCE {NR_TOLERANCE} != documented {expected_planet_tol}"
        )

        # NR_TOLERANCE_SUN: 0.001 arcsec
        expected_sun_tol = CROSSING_SUN_TOLERANCE_ARCSEC / 3600.0
        assert abs(NR_TOLERANCE_SUN - expected_sun_tol) < 1e-12, (
            f"NR_TOLERANCE_SUN {NR_TOLERANCE_SUN} != documented {expected_sun_tol}"
        )

        # NR_TOLERANCE_MOON: 0.05 arcsec
        expected_moon_tol = CROSSING_MOON_TOLERANCE_ARCSEC / 3600.0
        assert abs(NR_TOLERANCE_MOON - expected_moon_tol) < 1e-10, (
            f"NR_TOLERANCE_MOON {NR_TOLERANCE_MOON} != documented {expected_moon_tol}"
        )

    @pytest.mark.precision
    def test_sun_crossing_sub_milliarcsecond(self):
        """Verify Sun crossing achieves 0.001 arcsec precision as documented."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        target = 0.0  # Aries ingress

        jd_cross = ephem.solcross_ut(target, jd_start, 0)
        pos, _ = ephem.calc_ut(jd_cross, SUN, 0)

        diff = angle_diff(pos[0], target)
        diff_arcsec = diff * 3600

        assert diff_arcsec < CROSSING_SUN_TOLERANCE_ARCSEC, (
            f"Sun crossing diff {diff_arcsec:.6f} arcsec >= 0.001"
        )

    @pytest.mark.precision
    def test_moon_crossing_sub_arcsecond(self):
        """Verify Moon crossing achieves 0.05 arcsec precision as documented."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        target = 0.0  # Aries ingress

        jd_cross = ephem.mooncross_ut(target, jd_start, 0)
        pos, _ = ephem.calc_ut(jd_cross, MOON, 0)

        diff = angle_diff(pos[0], target)
        diff_arcsec = diff * 3600

        assert diff_arcsec < CROSSING_MOON_TOLERANCE_ARCSEC, (
            f"Moon crossing diff {diff_arcsec:.6f} arcsec >= 0.05"
        )


class TestDocumentedLunarPointsPrecision:
    """Verify documented lunar points precision claims."""

    @pytest.mark.precision
    @pytest.mark.parametrize(
        "point_id,name",
        [
            (MEAN_NODE, "Mean Node"),
            (TRUE_NODE, "True Node"),
            (MEAN_APOG, "Mean Lilith"),
            (OSCU_APOG, "True Lilith"),
        ],
    )
    def test_lunar_point_tolerance(self, point_id, name):
        """Verify lunar points tolerance matches documented values."""
        jd = 2451545.0

        pos_lib, _ = ephem.calc_ut(jd, point_id, 0)
        pos_swe, _ = swe.calc_ut(jd, point_id, 0)

        diff = angle_diff(pos_lib[0], pos_swe[0])
        tolerance = LUNAR_POINT_TOLERANCE_DEG[point_id]
        assert diff < tolerance, f"{name} diff {diff:.6f} degrees >= {tolerance}"


class TestDocumentedTimePrecision:
    """Verify documented time calculation precision claims."""

    @pytest.mark.precision
    def test_julian_day_precision(self):
        """Verify Julian Day precision < 1e-10 days as documented."""
        test_dates = [
            (2000, 1, 1, 12.0),  # J2000
            (1900, 1, 1, 0.0),  # DE421 start
            (2050, 12, 31, 23.99),  # DE421 end
        ]

        for year, month, day, hour in test_dates:
            jd_lib = ephem.julday(year, month, day, hour)
            jd_swe = swe.julday(year, month, day, hour)

            diff = abs(jd_lib - jd_swe)
            assert diff < JULIAN_DAY_TOLERANCE, (
                f"Julian Day diff at {year}-{month}-{day} {hour}h: {diff} >= 1e-10 days"
            )


class TestDocumentedVelocityPrecision:
    """Verify documented velocity precision claims."""

    @pytest.mark.precision
    def test_angular_velocity_tolerance(self):
        """Verify angular velocity tolerance < 0.01 degrees/day as documented."""
        jd = 2451545.0

        for planet_id in [SUN, MOON, MERCURY, VENUS, MARS]:
            pos_lib, _ = ephem.calc_ut(jd, planet_id, FLG_SPEED)
            pos_swe, _ = swe.calc_ut(jd, planet_id, FLG_SPEED)

            # Longitude speed (degrees/day)
            lon_speed_diff = abs(pos_lib[3] - pos_swe[3])
            assert lon_speed_diff < 0.01, (
                f"Planet {planet_id} longitude speed diff {lon_speed_diff:.6f} "
                f"degrees/day >= 0.01"
            )


class TestDocumentedAsteroidPrecision:
    """Verify documented asteroid precision claims."""

    @pytest.mark.precision
    def test_chiron_keplerian_approximation(self):
        """
        Verify Chiron uses Keplerian approximation with documented tolerance.

        Note: This test documents that Chiron uses simplified Keplerian
        propagation, not full numerical integration. The tolerance is
        approximately 1.5 degrees due to perturbations from giant planets
        that are not modeled in the Keplerian approximation.
        """
        jd = 2451545.0

        pos_lib, _ = ephem.calc_ut(jd, CHIRON, 0)

        # Try to get pyswisseph result, but skip if asteroid files not available
        try:
            pos_swe, _ = swe.calc_ut(jd, CHIRON, 0)
        except Exception as e:
            if "not found" in str(e).lower() or "seas" in str(e).lower():
                pytest.skip(f"Swiss Ephemeris asteroid files not available: {e}")
            raise

        diff = angle_diff(pos_lib[0], pos_swe[0])

        # Documented tolerance for asteroids: ~1.5 degrees (Keplerian)
        # Chiron's orbit is strongly perturbed by Saturn and Uranus,
        # which can cause deviations of 1-2 degrees from Keplerian propagation
        assert diff < 2.0, (
            f"Chiron diff {diff:.4f} degrees >= 2.0 (Keplerian tolerance)"
        )


class TestDocumentedLimitations:
    """Verify documented limitations are correctly described."""

    @pytest.mark.precision
    def test_fixed_star_velocity_is_small(self):
        """Verify fixed star velocities are small but non-zero due to precession."""
        jd = 2451545.0

        # fixstar_ut returns (position, star_name, retflag) — raises Error on failure
        result = ephem.fixstar_ut("Regulus", jd, FLG_SPEED)
        pos = result[0]  # First element is the position tuple

        # Fixed stars have small non-zero velocities due to precession of equinoxes
        # (approximately 3.8e-05 degrees/day), consistent with pyswisseph behavior
        assert abs(pos[3]) < 0.001, "Fixed star lon velocity should be very small"
        assert abs(pos[4]) < 0.001, "Fixed star lat velocity should be very small"
        assert abs(pos[5]) < 0.001, "Fixed star dist velocity should be very small"

    @pytest.mark.precision
    def test_equal_whole_sign_work_at_polar_latitudes(self):
        """Verify Equal and Whole Sign work at polar latitudes as documented."""
        jd = 2451545.0
        polar_lat = 85.0  # Very high latitude
        lon = 0.0

        # Equal houses should work
        cusps_equal, ascmc_equal = ephem.houses(jd, polar_lat, lon, ord("E"))
        assert len(cusps_equal) == 12
        assert 0 <= ascmc_equal[0] < 360

        # Whole Sign houses should work
        cusps_ws, ascmc_ws = ephem.houses(jd, polar_lat, lon, ord("W"))
        assert len(cusps_ws) == 12
        assert 0 <= ascmc_ws[0] < 360


class TestDocumentedInterpolatedApogeeWarning:
    """Verify interpolated apogee precision limitation warning exists in docs."""

    @pytest.mark.precision
    def test_interpolated_apogee_precision_warning_exists(self):
        """
        Verify PRECISION.md contains documentation about interpolated apogee/perigee
        precision and the ELP2000-82B perturbation series implementation.
        """
        import pathlib

        precision_doc = (
            pathlib.Path(__file__).parent.parent.parent.parent / "docs" / "PRECISION.md"
        )

        assert precision_doc.exists(), "PRECISION.md not found"

        content = precision_doc.read_text()

        # Check that the documentation mentions the ELP2000-82B implementation
        assert "ELP2000-82B perturbation series" in content, (
            "Missing documentation about ELP2000-82B perturbation series"
        )

        # Check that interpolated apogee precision is documented
        assert "Interpolated Apogee" in content, (
            "Missing interpolated apogee documentation"
        )

        # Check that interpolated perigee precision is documented
        assert "Interpolated Perigee" in content, (
            "Missing interpolated perigee documentation"
        )
