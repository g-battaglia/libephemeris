"""
Tests for house cusp velocity calculations in houses_ex2 and houses_armc_ex2.

Velocities are always computed (matching the reference ephemeris behavior), regardless of
whether FLG_SPEED is set. The FLG_SPEED flag is accepted but has no effect
on speed computation for house cusps.
"""

import pytest
import libephemeris as ephem
from libephemeris.constants import FLG_SPEED, FLG_SIDEREAL, SIDM_LAHIRI


class TestHousesEx2SpeedFlag:
    """Tests for velocities in houses_ex2."""

    @pytest.mark.unit
    def test_houses_ex2_always_computes_velocities(self):
        """Velocities are always computed, even without FLG_SPEED."""
        jd = 2451545.0  # J2000.0
        lat, lon = 41.9, 12.5

        cusps, ascmc, cusps_speed, ascmc_speed = ephem.houses_ex2(
            jd,
            lat,
            lon,
            ord("P"),
            0,  # No FLG_SPEED
        )

        # All cusp velocities should be non-zero (speeds always computed)
        for i, speed in enumerate(cusps_speed):
            assert speed != 0.0, f"Cusp {i + 1} speed should be non-zero"

        # ASC and MC velocities should be non-zero
        assert ascmc_speed[0] != 0.0, "ASC speed should be non-zero"
        assert ascmc_speed[1] != 0.0, "MC speed should be non-zero"

    @pytest.mark.unit
    def test_houses_ex2_with_speed_flag_returns_nonzero_velocities(self):
        """When FLG_SPEED is set, velocities should be calculated."""
        jd = 2451545.0
        lat, lon = 41.9, 12.5

        cusps, ascmc, cusps_speed, ascmc_speed = ephem.houses_ex2(
            jd, lat, lon, ord("P"), FLG_SPEED
        )

        # All cusp velocities should be non-zero
        for i, speed in enumerate(cusps_speed):
            assert speed != 0.0, f"Cusp {i + 1} speed should be non-zero"
            # Reasonable range: ~100-600 deg/day at mid-latitudes
            assert 100 < abs(speed) < 600, f"Cusp {i + 1} speed {speed} out of range"

        # ASC and MC velocities should be non-zero
        assert ascmc_speed[0] != 0.0, "ASC speed should be non-zero"
        assert ascmc_speed[1] != 0.0, "MC speed should be non-zero"
        # MC should be close to 360 deg/day
        assert 300 < abs(ascmc_speed[1]) < 400, (
            f"MC speed {ascmc_speed[1]} out of range"
        )

    @pytest.mark.unit
    def test_houses_ex2_speeds_identical_with_or_without_flag(self):
        """Speeds should be identical regardless of FLG_SPEED flag."""
        jd = 2451545.0
        lat, lon = 41.9, 12.5

        cusps_a, ascmc_a, cusps_speed_a, ascmc_speed_a = ephem.houses_ex2(
            jd, lat, lon, ord("P"), 0
        )
        cusps_b, ascmc_b, cusps_speed_b, ascmc_speed_b = ephem.houses_ex2(
            jd, lat, lon, ord("P"), FLG_SPEED
        )

        # Positions should be identical
        for i in range(12):
            assert cusps_a[i] == cusps_b[i], f"Cusp {i + 1} position differs"

        for i in range(8):
            assert ascmc_a[i] == ascmc_b[i], f"ASCMC {i} position differs"

        # Speeds should also be identical
        assert cusps_speed_a == cusps_speed_b, "Cusp speeds differ with/without flag"
        assert ascmc_speed_a == ascmc_speed_b, "ASCMC speeds differ with/without flag"

    @pytest.mark.unit
    def test_houses_ex2_with_sidereal_and_speed_flag(self):
        """FLG_SPEED should work together with FLG_SIDEREAL."""
        jd = 2451545.0
        lat, lon = 41.9, 12.5

        # Set Lahiri ayanamsa
        ephem.set_sid_mode(SIDM_LAHIRI)

        # With both flags
        cusps, ascmc, cusps_speed, ascmc_speed = ephem.houses_ex2(
            jd, lat, lon, ord("P"), FLG_SIDEREAL | FLG_SPEED
        )

        # Velocities should be calculated
        for i, speed in enumerate(cusps_speed):
            assert speed != 0.0, f"Cusp {i + 1} speed should be non-zero"

        # With FLG_SIDEREAL but not FLG_SPEED — speeds still computed
        cusps2, ascmc2, cusps_speed2, ascmc_speed2 = ephem.houses_ex2(
            jd, lat, lon, ord("P"), FLG_SIDEREAL
        )

        # Velocities should also be non-zero (always computed)
        for i, speed in enumerate(cusps_speed2):
            assert speed != 0.0, (
                f"Cusp {i + 1} speed should be non-zero even without FLG_SPEED"
            )

        # Speeds should be identical with or without FLG_SPEED
        assert cusps_speed == cusps_speed2

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "hsys,name",
        [
            (ord("P"), "Placidus"),
            (ord("K"), "Koch"),
            (ord("O"), "Porphyry"),
            (ord("R"), "Regiomontanus"),
            (ord("C"), "Campanus"),
            (ord("E"), "Equal"),
            (ord("B"), "Alcabitius"),
            (ord("M"), "Morinus"),
        ],
    )
    def test_houses_ex2_speed_flag_various_systems(self, hsys, name):
        """Velocities should be computed for all house systems."""
        jd = 2451545.0
        lat, lon = 45.0, 12.0

        # Without FLG_SPEED
        cusps1, ascmc1, cusps_speed1, ascmc_speed1 = ephem.houses_ex2(
            jd, lat, lon, hsys, 0
        )

        # With FLG_SPEED
        cusps2, ascmc2, cusps_speed2, ascmc_speed2 = ephem.houses_ex2(
            jd, lat, lon, hsys, FLG_SPEED
        )

        # Positions should match
        assert cusps1 == cusps2, f"{name}: positions differ"
        assert ascmc1 == ascmc2, f"{name}: ascmc differ"

        # Both should have non-zero velocities (speeds always computed)
        assert any(s != 0.0 for s in cusps_speed1), (
            f"{name}: speed should be non-zero without flag"
        )
        assert any(s != 0.0 for s in cusps_speed2), (
            f"{name}: speed should be non-zero with flag"
        )

        # Speeds should be identical
        assert cusps_speed1 == cusps_speed2, f"{name}: speeds differ with/without flag"

    @pytest.mark.unit
    def test_houses_ex2_speed_flag_whole_sign(self):
        """Whole Sign cusps are fixed at sign boundaries, so velocity is zero.

        Whole Sign houses place cusps at exact sign boundaries (0, 30, 60, etc.).
        Within a short time interval, these don't change, so cusp velocities are zero.
        However, ASC and MC velocities should still be calculated.
        """
        jd = 2451545.0
        lat, lon = 45.0, 12.0

        # Without FLG_SPEED
        cusps1, ascmc1, cusps_speed1, ascmc_speed1 = ephem.houses_ex2(
            jd, lat, lon, ord("W"), 0
        )

        # With FLG_SPEED
        cusps2, ascmc2, cusps_speed2, ascmc_speed2 = ephem.houses_ex2(
            jd, lat, lon, ord("W"), FLG_SPEED
        )

        # Positions should match
        assert cusps1 == cusps2, "Whole Sign: positions differ"
        assert ascmc1 == ascmc2, "Whole Sign: ascmc differ"

        # Speeds should be identical with or without flag
        assert cusps_speed1 == cusps_speed2, "Whole Sign: speeds differ"
        assert ascmc_speed1 == ascmc_speed2, "Whole Sign: ascmc speeds differ"

        # Most cusp velocities are zero (fixed at sign boundaries)
        # but cusps 1,7 (ASC/DESC) and 4,10 (IC/MC) get ASC/MC speeds
        # to match the reference ephemeris behaviour.
        for i in [1, 2, 4, 5, 7, 8]:  # 0-indexed: cusps 2,3,5,6,8,9
            assert cusps_speed1[i] == 0.0, (
                f"Whole Sign cusp {i + 1} should be 0, got {cusps_speed1[i]}"
            )
        # Cusps at ASC/DESC/MC/IC positions have non-zero speeds
        assert cusps_speed1[0] != 0.0, "Cusp 1 (ASC) speed should be non-zero"
        assert cusps_speed1[3] != 0.0, "Cusp 4 (IC) speed should be non-zero"
        assert cusps_speed1[6] != 0.0, "Cusp 7 (DESC) speed should be non-zero"
        assert cusps_speed1[9] != 0.0, "Cusp 10 (MC) speed should be non-zero"

    """Tests for velocities in houses_armc_ex2."""

    @pytest.mark.unit
    def test_houses_armc_ex2_always_computes_velocities(self):
        """Velocities are always computed in houses_armc_ex2."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393

        cusps, ascmc, cusps_speed, ascmc_speed = ephem.houses_armc_ex2(
            armc,
            lat,
            eps,
            ord("P"),
        )

        # All cusp velocities should be non-zero (speeds always computed)
        for i, speed in enumerate(cusps_speed):
            assert speed != 0.0, f"Cusp {i + 1} speed should be non-zero"

        # ASC and MC velocities should be non-zero
        assert ascmc_speed[0] != 0.0, "ASC speed should be non-zero"
        assert ascmc_speed[1] != 0.0, "MC speed should be non-zero"

    @pytest.mark.unit
    def test_houses_armc_ex2_velocity_ranges(self):
        """Cusp velocities should be in reasonable ranges."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393

        cusps, ascmc, cusps_speed, ascmc_speed = ephem.houses_armc_ex2(
            armc, lat, eps, ord("P")
        )

        # All cusp velocities should be non-zero
        for i, speed in enumerate(cusps_speed):
            assert speed != 0.0, f"Cusp {i + 1} speed should be non-zero"
            # Reasonable range: ~100-600 deg/day at mid-latitudes
            assert 100 < abs(speed) < 600, f"Cusp {i + 1} speed {speed} out of range"

        # ASC and MC velocities should be non-zero
        assert ascmc_speed[0] != 0.0, "ASC speed should be non-zero"
        assert ascmc_speed[1] != 0.0, "MC speed should be non-zero"
        # ARMC velocity should be the sidereal rotation rate (~360.986 deg/day)
        assert abs(ascmc_speed[2] - 360.9856) < 0.01, (
            f"ARMC speed {ascmc_speed[2]} should be ~360.986"
        )

    @pytest.mark.unit
    def test_houses_armc_ex2_positions_independent_of_ascmc9(self):
        """Cusp positions should be identical regardless of ascmc9 value."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393

        cusps_a, ascmc_a, speed_a, aspeed_a = ephem.houses_armc_ex2(
            armc, lat, eps, ord("P")
        )
        cusps_b, ascmc_b, speed_b, aspeed_b = ephem.houses_armc_ex2(
            armc, lat, eps, ord("P"), 0.0
        )

        # Positions should be identical
        for i in range(12):
            assert cusps_a[i] == cusps_b[i], f"Cusp {i + 1} position differs"

        for i in range(8):
            assert ascmc_a[i] == ascmc_b[i], f"ASCMC {i} position differs"

        # Speeds should be identical
        assert speed_a == speed_b, "Cusp speeds differ"
        assert aspeed_a == aspeed_b, "ASCMC speeds differ"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "hsys,name",
        [
            (ord("P"), "Placidus"),
            (ord("K"), "Koch"),
            (ord("O"), "Porphyry"),
            (ord("R"), "Regiomontanus"),
            (ord("C"), "Campanus"),
            (ord("E"), "Equal"),
            (ord("B"), "Alcabitius"),
            (ord("M"), "Morinus"),
        ],
    )
    def test_houses_armc_ex2_speed_flag_various_systems(self, hsys, name):
        """Velocities should be computed for all house systems."""
        armc = 150.0
        lat = 45.0
        eps = 23.44

        cusps, ascmc, cusps_speed, ascmc_speed = ephem.houses_armc_ex2(
            armc, lat, eps, hsys
        )

        # Should have non-zero velocities (speeds always computed)
        assert any(s != 0.0 for s in cusps_speed), f"{name}: speed should be non-zero"

    @pytest.mark.unit
    def test_houses_armc_ex2_speed_flag_whole_sign(self):
        """Whole Sign cusps are fixed at sign boundaries, so velocity is zero.

        Whole Sign houses place cusps at exact sign boundaries (0, 30, 60, etc.).
        Within a short time interval, these don't change, so cusp velocities are zero.
        However, ASC and MC velocities should still be calculated.
        """
        armc = 150.0
        lat = 45.0
        eps = 23.44

        cusps, ascmc, cusps_speed, ascmc_speed = ephem.houses_armc_ex2(
            armc, lat, eps, ord("W")
        )

        # Most cusp velocities are zero (fixed at sign boundaries)
        # but cusps 1,7 (ASC/DESC) and 4,10 (IC/MC) get ASC/MC speeds
        # to match the reference ephemeris behaviour.
        for i in [1, 2, 4, 5, 7, 8]:  # 0-indexed: cusps 2,3,5,6,8,9
            assert cusps_speed[i] == 0.0, (
                f"Whole Sign cusp {i + 1} should be 0, got {cusps_speed[i]}"
            )
        # Cusps at ASC/DESC/MC/IC positions have non-zero speeds
        assert cusps_speed[0] != 0.0, "Cusp 1 (ASC) speed should be non-zero"
        assert cusps_speed[3] != 0.0, "Cusp 4 (IC) speed should be non-zero"
        assert cusps_speed[6] != 0.0, "Cusp 7 (DESC) speed should be non-zero"
        assert cusps_speed[9] != 0.0, "Cusp 10 (MC) speed should be non-zero"


class TestHouseCuspVelocityEdgeCases:
    """Edge case tests for house cusp velocity calculations."""

    @pytest.mark.edge_case
    def test_houses_ex2_equator_with_speed(self):
        """Velocities should be computed at the equator."""
        jd = 2451545.0
        lat, lon = 0.0, 0.0

        cusps, ascmc, cusps_speed, ascmc_speed = ephem.houses_ex2(
            jd, lat, lon, ord("P"), FLG_SPEED
        )

        # Velocities should be calculated
        assert any(s != 0.0 for s in cusps_speed)
        assert ascmc_speed[1] != 0.0  # MC velocity

    @pytest.mark.edge_case
    def test_houses_ex2_high_latitude_with_speed(self):
        """Velocities should be computed at high latitudes."""
        jd = 2451545.0
        lat, lon = 60.0, 0.0

        cusps, ascmc, cusps_speed, ascmc_speed = ephem.houses_ex2(
            jd, lat, lon, ord("P"), FLG_SPEED
        )

        # Velocities should be calculated
        assert any(s != 0.0 for s in cusps_speed)

    @pytest.mark.edge_case
    def test_houses_ex2_southern_hemisphere_with_speed(self):
        """Velocities should be computed in the southern hemisphere."""
        jd = 2451545.0
        lat, lon = -33.9, 151.2  # Sydney

        cusps, ascmc, cusps_speed, ascmc_speed = ephem.houses_ex2(
            jd, lat, lon, ord("P"), FLG_SPEED
        )

        # Velocities should be calculated
        assert any(s != 0.0 for s in cusps_speed)

    @pytest.mark.edge_case
    def test_houses_armc_ex2_equator_with_speed(self):
        """Velocities should be computed at the equator for ARMC-based calculation."""
        armc = 180.0
        lat = 0.0
        eps = 23.44

        cusps, ascmc, cusps_speed, ascmc_speed = ephem.houses_armc_ex2(
            armc, lat, eps, ord("P")
        )

        # Velocities should be calculated
        assert any(s != 0.0 for s in cusps_speed)


class TestHouseCuspVelocityProgression:
    """Tests for progressed chart applications."""

    @pytest.mark.unit
    def test_velocity_for_progressed_mc_crossing(self):
        """Velocities should allow calculating when MC crosses a point."""
        jd = 2451545.0
        lat, lon = 45.0, 0.0

        cusps, ascmc, cusps_speed, ascmc_speed = ephem.houses_ex2(
            jd, lat, lon, ord("P"), FLG_SPEED
        )

        # MC velocity should be positive (moving forward)
        mc_speed = ascmc_speed[1]
        assert mc_speed > 0, f"MC speed {mc_speed} should be positive"

        # Calculate time for MC to move 1 degree
        # time = distance / velocity (in days)
        time_for_1_deg = 1.0 / mc_speed

        # Should be approximately 4 minutes (1/360 of a day)
        expected_time = 1.0 / 360.0  # ~0.00278 days
        assert abs(time_for_1_deg - expected_time) < 0.001, (
            f"Time for 1 deg: {time_for_1_deg}, expected ~{expected_time}"
        )

    @pytest.mark.unit
    def test_velocities_consistent_across_time(self):
        """Velocities should be consistent when calculated at different times."""
        lat, lon = 45.0, 0.0
        jd1 = 2451545.0
        jd2 = 2451545.5  # 12 hours later

        _, _, _, ascmc_speed1 = ephem.houses_ex2(
            jd1, lat, lon, ord("P"), FLG_SPEED
        )
        _, _, _, ascmc_speed2 = ephem.houses_ex2(
            jd2, lat, lon, ord("P"), FLG_SPEED
        )

        # MC velocity should be similar at different times
        # (variation possible due to Earth's orbital motion)
        mc_speed1 = ascmc_speed1[1]
        mc_speed2 = ascmc_speed2[1]

        # Should be within 1% of each other
        rel_diff = abs(mc_speed1 - mc_speed2) / mc_speed1
        assert rel_diff < 0.01, f"MC speed varies too much: {mc_speed1} vs {mc_speed2}"
