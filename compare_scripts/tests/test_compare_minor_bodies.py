"""
Minor Bodies (Asteroids, TNOs) Comparison Tests.

Compares minor body calculations between pyswisseph and libephemeris.
Tests cover all minor bodies defined in constants.py.

NOTE: With SPK auto-download enabled (via conftest), libephemeris uses JPL SPK
ephemeris files for common minor bodies (Ceres, Pallas, Juno, Vesta, Chiron,
Pholus, Eris, Sedna), providing high-precision results comparable to pyswisseph.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    # Main belt asteroids
    CERES,
    PALLAS,
    JUNO,
    VESTA,
    # Centaurs (original)
    CHIRON,
    PHOLUS,
    # Additional Centaurs
    NESSUS,
    ASBOLUS,
    CHARIKLO,
    # Trans-Neptunian Objects (TNOs)
    ERIS,
    SEDNA,
    HAUMEA,
    MAKEMAKE,
    IXION,
    ORCUS,
    QUAOAR,
    VARUNA,
    GONGGONG,
    # Large main belt asteroids
    HYGIEA,
    INTERAMNIA,
    DAVIDA,
    EUROPA_AST,
    SYLVIA,
    PSYCHE,
    # Near-Earth asteroids
    APOPHIS,
    EROS,
    AMOR,
    ICARUS,
    TORO,
    TOUTATIS,
    ITOKAWA,
    BENNU,
    RYUGU,
    # Astrologically significant asteroids
    SAPPHO,
    PANDORA_AST,
    LILITH_AST,
    HIDALGO,
)


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# ============================================================================
# TOLERANCES
# ============================================================================

MAIN_ASTEROID_TOL = 0.01  # degrees - main belt asteroids with good ephemeris
CENTAUR_TOL = 0.01  # degrees - centaurs with good ephemeris
TNO_TOL = 0.1  # degrees - relaxed for distant TNOs
NEA_TOL = 0.05  # degrees - Near-Earth asteroids (variable orbits)
DISTANT_TNO_TOL = 0.5  # degrees - very distant/uncertain TNOs like Sedna


# ============================================================================
# TEST DATA
# ============================================================================

MAIN_ASTEROIDS = [
    (CERES, "Ceres"),
    (PALLAS, "Pallas"),
    (JUNO, "Juno"),
    (VESTA, "Vesta"),
]

LARGE_MAIN_BELT = [
    (HYGIEA, "Hygiea"),
    (INTERAMNIA, "Interamnia"),
    (DAVIDA, "Davida"),
    (EUROPA_AST, "Europa (asteroid)"),
    (SYLVIA, "Sylvia"),
    (PSYCHE, "Psyche"),
]

CENTAURS = [
    (CHIRON, "Chiron"),
    (PHOLUS, "Pholus"),
]

ADDITIONAL_CENTAURS = [
    (NESSUS, "Nessus"),
    (ASBOLUS, "Asbolus"),
    (CHARIKLO, "Chariklo"),
]

TNOS_STANDARD = [
    (ERIS, "Eris"),
    (HAUMEA, "Haumea"),
    (MAKEMAKE, "Makemake"),
    (ORCUS, "Orcus"),
    (QUAOAR, "Quaoar"),
    (VARUNA, "Varuna"),
    (IXION, "Ixion"),
    (GONGGONG, "Gonggong"),
]

TNOS_DISTANT = [
    (SEDNA, "Sedna"),  # Very distant, relaxed tolerance
]

NEAR_EARTH_ASTEROIDS = [
    (EROS, "Eros"),
    (AMOR, "Amor"),
    (ICARUS, "Icarus"),
    (TORO, "Toro"),
    (APOPHIS, "Apophis"),
    (TOUTATIS, "Toutatis"),
    (ITOKAWA, "Itokawa"),
    (BENNU, "Bennu"),
    (RYUGU, "Ryugu"),
]

ASTROLOGICAL_ASTEROIDS = [
    (SAPPHO, "Sappho"),
    (PANDORA_AST, "Pandora (asteroid)"),
    (LILITH_AST, "Lilith (asteroid)"),
    (HIDALGO, "Hidalgo"),
]

TEST_DATES = [
    (2000, 1, 1, 12.0, "J2000"),
    (2024, 6, 15, 0.0, "Current"),
    (1980, 5, 20, 14.5, "Past"),
]


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestMainAsteroids:
    """Compare main belt asteroid calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", MAIN_ASTEROIDS)
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    def test_main_asteroid_position(
        self, body_id, body_name, year, month, day, hour, desc
    ):
        """Test main asteroid positions match within tolerance."""
        jd = swe.julday(year, month, day, hour)

        try:
            pos_swe, _ = swe.calc_ut(jd, body_id, 0)
            pos_py, _ = ephem.calc_ut(jd, body_id, 0)
        except Exception as e:
            pytest.skip(f"{body_name} not available: {e}")
            return

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < MAIN_ASTEROID_TOL, (
            f"{body_name} at {desc}: longitude diff {diff:.6f}° exceeds tolerance"
        )


class TestLargeMainBeltAsteroids:
    """Compare large main belt asteroids (Hygiea, Interamnia, Davida, etc.)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", LARGE_MAIN_BELT)
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    def test_large_main_belt_position(
        self, body_id, body_name, year, month, day, hour, desc
    ):
        """Test large main belt asteroid positions match within tolerance."""
        jd = swe.julday(year, month, day, hour)

        try:
            pos_swe, _ = swe.calc_ut(jd, body_id, 0)
            pos_py, _ = ephem.calc_ut(jd, body_id, 0)
        except Exception as e:
            pytest.skip(f"{body_name} not available: {e}")
            return

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < MAIN_ASTEROID_TOL, (
            f"{body_name} at {desc}: longitude diff {diff:.6f}° exceeds tolerance"
        )


class TestCentaurs:
    """Compare centaur calculations (Chiron, Pholus)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", CENTAURS)
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    def test_centaur_position(self, body_id, body_name, year, month, day, hour, desc):
        """Test centaur positions match within tolerance."""
        jd = swe.julday(year, month, day, hour)

        try:
            pos_swe, _ = swe.calc_ut(jd, body_id, 0)
            pos_py, _ = ephem.calc_ut(jd, body_id, 0)
        except Exception as e:
            pytest.skip(f"{body_name} not available: {e}")
            return

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < CENTAUR_TOL, (
            f"{body_name} at {desc}: longitude diff {diff:.6f}° exceeds tolerance"
        )


class TestAdditionalCentaurs:
    """Compare additional centaur calculations (Nessus, Asbolus, Chariklo)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", ADDITIONAL_CENTAURS)
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    def test_additional_centaur_position(
        self, body_id, body_name, year, month, day, hour, desc
    ):
        """Test additional centaur positions match within tolerance."""
        jd = swe.julday(year, month, day, hour)

        try:
            pos_swe, _ = swe.calc_ut(jd, body_id, 0)
            pos_py, _ = ephem.calc_ut(jd, body_id, 0)
        except Exception as e:
            pytest.skip(f"{body_name} not available: {e}")
            return

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < CENTAUR_TOL, (
            f"{body_name} at {desc}: longitude diff {diff:.6f}° exceeds tolerance"
        )


class TestTNOsStandard:
    """Compare Trans-Neptunian Object calculations (Eris, Makemake, etc.)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", TNOS_STANDARD)
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    def test_tno_position(self, body_id, body_name, year, month, day, hour, desc):
        """Test TNO positions match within relaxed tolerance."""
        jd = swe.julday(year, month, day, hour)

        try:
            pos_swe, _ = swe.calc_ut(jd, body_id, 0)
            pos_py, _ = ephem.calc_ut(jd, body_id, 0)
        except Exception as e:
            pytest.skip(f"{body_name} not available: {e}")
            return

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < TNO_TOL, (
            f"{body_name} at {desc}: longitude diff {diff:.6f}° exceeds tolerance"
        )


class TestTNOsDistant:
    """Compare very distant TNO calculations (Sedna) with extra relaxed tolerance."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", TNOS_DISTANT)
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    def test_distant_tno_position(
        self, body_id, body_name, year, month, day, hour, desc
    ):
        """Test distant TNO positions match within extra relaxed tolerance."""
        jd = swe.julday(year, month, day, hour)

        try:
            pos_swe, _ = swe.calc_ut(jd, body_id, 0)
            pos_py, _ = ephem.calc_ut(jd, body_id, 0)
        except Exception as e:
            pytest.skip(f"{body_name} not available: {e}")
            return

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < DISTANT_TNO_TOL, (
            f"{body_name} at {desc}: longitude diff {diff:.6f}° exceeds tolerance"
        )


class TestNearEarthAsteroids:
    """Compare Near-Earth Asteroid calculations (Apophis, Eros, etc.)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", NEAR_EARTH_ASTEROIDS)
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    def test_nea_position(self, body_id, body_name, year, month, day, hour, desc):
        """Test Near-Earth asteroid positions match within tolerance."""
        jd = swe.julday(year, month, day, hour)

        try:
            pos_swe, _ = swe.calc_ut(jd, body_id, 0)
            pos_py, _ = ephem.calc_ut(jd, body_id, 0)
        except Exception as e:
            pytest.skip(f"{body_name} not available: {e}")
            return

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < NEA_TOL, (
            f"{body_name} at {desc}: longitude diff {diff:.6f}° exceeds tolerance"
        )


class TestAstrologicalAsteroids:
    """Compare astrologically significant asteroids (Sappho, Lilith, etc.)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", ASTROLOGICAL_ASTEROIDS)
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    def test_astrological_asteroid_position(
        self, body_id, body_name, year, month, day, hour, desc
    ):
        """Test astrological asteroid positions match within tolerance."""
        jd = swe.julday(year, month, day, hour)

        try:
            pos_swe, _ = swe.calc_ut(jd, body_id, 0)
            pos_py, _ = ephem.calc_ut(jd, body_id, 0)
        except Exception as e:
            pytest.skip(f"{body_name} not available: {e}")
            return

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < MAIN_ASTEROID_TOL, (
            f"{body_name} at {desc}: longitude diff {diff:.6f}° exceeds tolerance"
        )


class TestChironSpecific:
    """Specific tests for Chiron (most commonly used asteroid)."""

    @pytest.mark.comparison
    def test_chiron_at_j2000(self):
        """Test Chiron position at J2000 epoch."""
        jd = 2451545.0

        try:
            pos_swe, _ = swe.calc_ut(jd, swe.CHIRON, 0)
            pos_py, _ = ephem.calc_ut(jd, CHIRON, 0)
        except Exception as e:
            pytest.skip(f"Chiron not available: {e}")
            return

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < 0.001, f"Chiron at J2000 diff {diff:.6f}° exceeds tight tolerance"

    @pytest.mark.comparison
    def test_chiron_with_speed(self):
        """Test Chiron position with velocity."""
        jd = swe.julday(2024, 1, 1, 12.0)

        try:
            pos_swe, _ = swe.calc_ut(jd, swe.CHIRON, swe.FLG_SPEED)
            pos_py, _ = ephem.calc_ut(jd, CHIRON, 256)  # FLG_SPEED
        except Exception as e:
            pytest.skip(f"Chiron not available: {e}")
            return

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_speed = abs(pos_swe[3] - pos_py[3])

        assert diff_lon < CENTAUR_TOL, f"Chiron longitude diff {diff_lon:.6f}°"
        assert diff_speed < 0.01, f"Chiron speed diff {diff_speed:.6f}°/day"


class TestCeresSpecific:
    """Specific tests for Ceres (dwarf planet)."""

    @pytest.mark.comparison
    def test_ceres_at_j2000(self):
        """Test Ceres position at J2000 epoch."""
        jd = 2451545.0

        try:
            pos_swe, _ = swe.calc_ut(jd, swe.CERES, 0)
            pos_py, _ = ephem.calc_ut(jd, CERES, 0)
        except Exception as e:
            pytest.skip(f"Ceres not available: {e}")
            return

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < 0.001, f"Ceres at J2000 diff {diff:.6f}° exceeds tight tolerance"


class TestErisSpecific:
    """Specific tests for Eris (largest known dwarf planet)."""

    @pytest.mark.comparison
    def test_eris_at_j2000(self):
        """Test Eris position at J2000 epoch."""
        jd = 2451545.0

        try:
            pos_swe, _ = swe.calc_ut(jd, ERIS, 0)
            pos_py, _ = ephem.calc_ut(jd, ERIS, 0)
        except Exception as e:
            pytest.skip(f"Eris not available: {e}")
            return

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < TNO_TOL, f"Eris at J2000 diff {diff:.6f}° exceeds tolerance"


class TestGonggongSpecific:
    """Specific tests for Gonggong (formerly 2007 OR10)."""

    @pytest.mark.comparison
    def test_gonggong_at_current_epoch(self):
        """Test Gonggong position at current epoch."""
        jd = swe.julday(2024, 6, 15, 12.0)

        try:
            pos_swe, _ = swe.calc_ut(jd, GONGGONG, 0)
            pos_py, _ = ephem.calc_ut(jd, GONGGONG, 0)
        except Exception as e:
            pytest.skip(f"Gonggong not available: {e}")
            return

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < TNO_TOL, f"Gonggong diff {diff:.6f}° exceeds tolerance"


class TestApophisSpecific:
    """Specific tests for Apophis (potentially hazardous asteroid)."""

    @pytest.mark.comparison
    def test_apophis_at_current_epoch(self):
        """Test Apophis position at current epoch."""
        jd = swe.julday(2024, 6, 15, 12.0)

        try:
            pos_swe, _ = swe.calc_ut(jd, APOPHIS, 0)
            pos_py, _ = ephem.calc_ut(jd, APOPHIS, 0)
        except Exception as e:
            pytest.skip(f"Apophis not available: {e}")
            return

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < NEA_TOL, f"Apophis diff {diff:.6f}° exceeds tolerance"


class TestAsteroidVelocity:
    """Test asteroid velocity calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", MAIN_ASTEROIDS + CENTAURS)
    def test_asteroid_velocity(self, body_id, body_name):
        """Test asteroid velocity calculations match."""
        jd = swe.julday(2024, 6, 15, 12.0)

        try:
            pos_swe, _ = swe.calc_ut(jd, body_id, swe.FLG_SPEED)
            pos_py, _ = ephem.calc_ut(jd, body_id, 256)  # FLG_SPEED
        except Exception as e:
            pytest.skip(f"{body_name} not available: {e}")
            return

        diff_speed = abs(pos_swe[3] - pos_py[3])

        assert diff_speed < 0.01, (
            f"{body_name} velocity diff {diff_speed:.6f}°/day exceeds tolerance"
        )


class TestAllBodiesCoverage:
    """Verify all defined minor bodies have test coverage."""

    @pytest.mark.comparison
    def test_all_bodies_in_test_lists(self):
        """Verify all minor body constants are included in test lists."""
        all_tested_bodies = set()
        for body_list in [
            MAIN_ASTEROIDS,
            LARGE_MAIN_BELT,
            CENTAURS,
            ADDITIONAL_CENTAURS,
            TNOS_STANDARD,
            TNOS_DISTANT,
            NEAR_EARTH_ASTEROIDS,
            ASTROLOGICAL_ASTEROIDS,
        ]:
            for body_id, _ in body_list:
                all_tested_bodies.add(body_id)

        # These are all the SE_* minor body constants defined in constants.py
        expected_bodies = {
            CERES,
            PALLAS,
            JUNO,
            VESTA,
            CHIRON,
            PHOLUS,
            NESSUS,
            ASBOLUS,
            CHARIKLO,
            ERIS,
            SEDNA,
            HAUMEA,
            MAKEMAKE,
            IXION,
            ORCUS,
            QUAOAR,
            VARUNA,
            GONGGONG,
            HYGIEA,
            INTERAMNIA,
            DAVIDA,
            EUROPA_AST,
            SYLVIA,
            PSYCHE,
            APOPHIS,
            EROS,
            AMOR,
            ICARUS,
            TORO,
            TOUTATIS,
            ITOKAWA,
            BENNU,
            RYUGU,
            SAPPHO,
            PANDORA_AST,
            LILITH_AST,
            HIDALGO,
        }

        missing_bodies = expected_bodies - all_tested_bodies
        assert not missing_bodies, (
            f"Missing test coverage for body IDs: {missing_bodies}"
        )
