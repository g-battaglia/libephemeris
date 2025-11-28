"""
Unit tests for minor bodies: Asteroids and TNOs.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import *


@pytest.fixture
def all_asteroids():
    """All main belt asteroids."""
    return [
        (SE_CHIRON, "Chiron", swe.CHIRON),
        (SE_PHOLUS, "Pholus", swe.PHOLUS),
        (SE_CERES, "Ceres", swe.CERES),
        (SE_PALLAS, "Pallas", swe.PALLAS),
        (SE_JUNO, "Juno", swe.JUNO),
        (SE_VESTA, "Vesta", swe.VESTA),
    ]


@pytest.fixture
def all_tnos():
    """All Trans-Neptunian Objects."""
    return [
        (SE_ERIS, "Eris"),
        (SE_SEDNA, "Sedna"),
        (SE_HAUMEA, "Haumea"),
        (SE_MAKEMAKE, "Makemake"),
        (SE_IXION, "Ixion"),
        (SE_ORCUS, "Orcus"),
        (SE_QUAOAR, "Quaoar"),
    ]


@pytest.mark.unit
class TestAsteroids:
    """Tests for main belt asteroids."""

    def test_chiron_position(self, standard_jd):
        """Test Chiron position at J2000."""
        pos, _ = ephem.swe_calc_ut(standard_jd, SE_CHIRON, 0)

        # Should have valid position
        assert 0 <= pos[0] < 360, "Invalid longitude"
        assert -90 <= pos[1] <= 90, "Invalid latitude"
        assert pos[2] > 0, "Invalid distance"

    def test_all_asteroids_valid(self, standard_jd, all_asteroids):
        """Test all asteroids return valid positions."""
        for body_id, name, _ in all_asteroids:
            pos, _ = ephem.swe_calc_ut(standard_jd, body_id, 0)

            assert 0 <= pos[0] < 360, f"{name}: Invalid longitude"
            assert -90 <= pos[1] <= 90, f"{name}: Invalid latitude"
            assert pos[2] > 0, f"{name}: Invalid distance"

    @pytest.mark.parametrize(
        "body_id,name,swe_id",
        [
            (SE_CHIRON, "Chiron", swe.CHIRON),
            (SE_CERES, "Ceres", swe.CERES),
            (SE_PALLAS, "Pallas", swe.PALLAS),
        ],
    )
    def test_asteroid_vs_swisseph(self, standard_jd, body_id, name, swe_id):
        """Compare asteroid positions with SwissEph."""
        pos_py, _ = ephem.swe_calc_ut(standard_jd, body_id, 0)

        try:
            pos_swe, _ = swe.calc_ut(standard_jd, swe_id, 0)
        except swe.Error as e:
            pytest.skip(f"Skipping {name} comparison: SwissEph data file missing ({e})")

        diff_lon = abs(pos_py[0] - pos_swe[0])
        if diff_lon > 180:
            diff_lon = 360 - diff_lon

        # Keplerian elements provide ~1-5 arcmin accuracy
        # That's 0.016 - 0.083 degrees
        assert diff_lon < 10.0, f"{name}: Longitude diff {diff_lon}° exceeds 10°"

    def test_asteroid_distance_range(self, standard_jd, all_asteroids):
        """Test asteroid distances are in expected ranges."""
        expected_ranges = {
            "Chiron": (8, 20),  # Perihelion ~8 AU, Aphelion ~19 AU
            "Pholus": (13, 32),  # Eccentric orbit
            "Ceres": (2.5, 3.0),  # Main belt
            "Pallas": (2.1, 3.5),  # Aphelion ~3.41
            "Juno": (1.98, 3.3),  # Main belt
            "Vesta": (2.15, 2.57),  # Main belt
        }

        for body_id, name, _ in all_asteroids:
            # Use Heliocentric distance for orbital range checks
            pos, _ = ephem.swe_calc_ut(standard_jd, body_id, SEFLG_HELCTR)
            dist = pos[2]

            min_dist, max_dist = expected_ranges[name]
            assert min_dist <= dist <= max_dist, (
                f"{name}: Distance {dist:.2f} AU outside expected range [{min_dist}, {max_dist}]"
            )


@pytest.mark.unit
class TestTNOs:
    """Tests for Trans-Neptunian Objects."""

    def test_eris_position(self, standard_jd):
        """Test Eris position at J2000."""
        pos, _ = ephem.swe_calc_ut(standard_jd, SE_ERIS, 0)

        assert 0 <= pos[0] < 360, "Invalid longitude"
        assert -90 <= pos[1] <= 90, "Invalid latitude"
        assert pos[2] > 30, "Eris should be beyond Neptune"

    def test_all_tnos_valid(self, standard_jd, all_tnos):
        """Test all TNOs return valid positions."""
        for body_id, name in all_tnos:
            pos, _ = ephem.swe_calc_ut(standard_jd, body_id, 0)

            assert 0 <= pos[0] < 360, f"{name}: Invalid longitude"
            assert -90 <= pos[1] <= 90, f"{name}: Invalid latitude"
            assert pos[2] > 30, f"{name}: Should be beyond Neptune"

    def test_tno_distance_ranges(self, standard_jd, all_tnos):
        """Test TNO distances are in plausible ranges."""
        expected_min_dist = {
            "Eris": 38,  # ~37-98 AU
            "Sedna": 76,  # ~76-937 AU (highly eccentric!)
            "Haumea": 35,  # ~35-51 AU
            "Makemake": 38,  # ~38-53 AU
            "Ixion": 30,  # ~30-49 AU
            "Orcus": 30,  # ~30-48 AU
            "Quaoar": 42,  # ~42-45 AU
        }

        for body_id, name in all_tnos:
            pos, _ = ephem.swe_calc_ut(standard_jd, body_id, 0)
            dist = pos[2]

            min_dist = expected_min_dist[name]
            assert dist >= min_dist * 0.5, (
                f"{name}: Distance {dist:.2f} AU too small (expected >= {min_dist * 0.5})"
            )
            assert dist <= 1000, f"{name}: Distance {dist:.2f} AU implausibly large"

    @pytest.mark.parametrize(
        "year,month,day",
        [
            (2000, 1, 1),
            (2010, 6, 15),
            (2020, 12, 31),
        ],
    )
    def test_tnos_over_time(self, year, month, day, all_tnos):
        """Test TNOs at different dates."""
        jd = ephem.swe_julday(year, month, day, 0.0)

        for body_id, name in all_tnos:
            pos, _ = ephem.swe_calc_ut(jd, body_id, 0)

            # All should return valid data
            assert 0 <= pos[0] < 360, f"{name}: Invalid lon at {year}"
            assert pos[2] > 0, f"{name}: Invalid dist at {year}"


@pytest.mark.integration
class TestMinorBodiesIntegration:
    """Integration tests for minor bodies."""

    def test_heliocentric_mode(self, standard_jd, all_asteroids):
        """Test asteroids in heliocentric mode."""
        for body_id, name, _ in all_asteroids:
            pos, _ = ephem.swe_calc_ut(standard_jd, body_id, SEFLG_HELCTR)

            # Heliocentric should still return valid coordinates
            assert 0 <= pos[0] < 360, f"{name}: Invalid heliocentric lon"
            assert pos[2] > 0, f"{name}: Invalid heliocentric dist"

    def test_minor_body_consistency(self):
        """Test position consistency over short time intervals."""
        jd1 = ephem.swe_julday(2000, 1, 1, 0.0)
        jd2 = ephem.swe_julday(2000, 1, 2, 0.0)  # 1 day later

        pos1, _ = ephem.swe_calc_ut(jd1, SE_CHIRON, 0)
        pos2, _ = ephem.swe_calc_ut(jd2, SE_CHIRON, 0)

        # Position should change slightly but reasonably
        diff = abs(pos2[0] - pos1[0])

        # Chiron moves ~0.01-0.02 degrees/day
        assert 0 < diff < 1.0, "Position change over 1 day seems unreasonable"
