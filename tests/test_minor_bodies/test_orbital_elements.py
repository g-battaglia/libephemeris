"""
Tests for minor body orbital elements (epoch 2025.0).

Tests verify:
- Orbital elements are correctly defined and match JPL SBDB
- Elements produce valid heliocentric positions
- Mean motion values are consistent with semi-major axes
"""

import pytest
import math
import libephemeris as ephem
from libephemeris.constants import (
    CHIRON,
    PHOLUS,
    CERES,
    PALLAS,
    JUNO,
    VESTA,
    ERIS,
    SEDNA,
    HAUMEA,
    MAKEMAKE,
    IXION,
    ORCUS,
    QUAOAR,
    FLG_HELCTR,
)
from libephemeris.minor_bodies import MINOR_BODY_ELEMENTS


# Constants for orbital mechanics
GM_SUN = 0.00029591220828559  # GM of Sun in AU^3/day^2
DAYS_PER_YEAR = 365.25


@pytest.fixture
def epoch_2025():
    """Epoch JD 2461000.5 (2025-Sep-19)."""
    return 2461000.5


@pytest.fixture
def all_minor_body_ids():
    """All minor body IDs in the database."""
    return [
        CHIRON,
        PHOLUS,
        CERES,
        PALLAS,
        JUNO,
        VESTA,
        ERIS,
        SEDNA,
        HAUMEA,
        MAKEMAKE,
        IXION,
        ORCUS,
        QUAOAR,
    ]


class TestOrbitalElementsEpoch:
    """Verify orbital elements are defined at epoch 2025.0."""

    def test_epoch_is_2025(self, all_minor_body_ids):
        """All elements should be at epoch JD 2461000.5."""
        expected_epoch = 2461000.5
        for body_id in all_minor_body_ids:
            elements = MINOR_BODY_ELEMENTS[body_id]
            assert elements.epoch == expected_epoch, (
                f"{elements.name}: epoch should be {expected_epoch}, got {elements.epoch}"
            )

    def test_all_bodies_have_elements(self, all_minor_body_ids):
        """All expected minor bodies should have orbital elements."""
        for body_id in all_minor_body_ids:
            assert body_id in MINOR_BODY_ELEMENTS, (
                f"Missing elements for body {body_id}"
            )
            elements = MINOR_BODY_ELEMENTS[body_id]
            assert elements.name is not None
            assert elements.a > 0
            assert 0 <= elements.e < 1
            assert 0 <= elements.i <= 180
            assert 0 <= elements.omega < 360
            assert 0 <= elements.Omega < 360
            assert 0 <= elements.M0 < 360
            assert elements.n > 0


class TestMeanMotionConsistency:
    """Verify mean motion values are consistent with Kepler's 3rd law."""

    @pytest.mark.parametrize(
        "body_id,name",
        [
            (CHIRON, "Chiron"),
            (PHOLUS, "Pholus"),
            (CERES, "Ceres"),
            (PALLAS, "Pallas"),
            (JUNO, "Juno"),
            (VESTA, "Vesta"),
            (ERIS, "Eris"),
            (SEDNA, "Sedna"),
            (HAUMEA, "Haumea"),
            (MAKEMAKE, "Makemake"),
            (IXION, "Ixion"),
            (ORCUS, "Orcus"),
            (QUAOAR, "Quaoar"),
        ],
    )
    def test_mean_motion_matches_period(self, body_id, name):
        """Mean motion should match orbital period from Kepler's 3rd law."""
        elements = MINOR_BODY_ELEMENTS[body_id]

        # Calculate expected period from semi-major axis (Kepler's 3rd law)
        # T = 2*pi * sqrt(a^3 / GM_sun)
        period_days = 2 * math.pi * math.sqrt(elements.a**3 / GM_SUN)

        # Mean motion = 360 / period (in degrees per day)
        expected_n = 360.0 / period_days

        # Allow 5% tolerance for rounding
        rel_diff = abs(elements.n - expected_n) / expected_n
        assert rel_diff < 0.05, (
            f"{name}: mean motion mismatch - expected {expected_n:.6f}, "
            f"got {elements.n:.6f} (diff: {rel_diff * 100:.1f}%)"
        )


class TestOrbitalElementsValues:
    """Verify specific orbital element values match JPL SBDB."""

    def test_chiron_elements(self):
        """Chiron elements should match JPL values."""
        el = MINOR_BODY_ELEMENTS[CHIRON]
        assert el.name == "Chiron"
        assert abs(el.a - 13.6922) < 0.01  # Semi-major axis
        assert abs(el.e - 0.378979) < 0.001  # Eccentricity
        assert abs(el.i - 6.926) < 0.1  # Inclination

    def test_ceres_elements(self):
        """Ceres elements should match JPL values."""
        el = MINOR_BODY_ELEMENTS[CERES]
        assert el.name == "Ceres"
        assert abs(el.a - 2.7656) < 0.01
        assert abs(el.e - 0.079576) < 0.001
        assert abs(el.i - 10.588) < 0.1

    def test_eris_elements(self):
        """Eris elements should match JPL values."""
        el = MINOR_BODY_ELEMENTS[ERIS]
        assert el.name == "Eris"
        assert abs(el.a - 67.9964) < 0.1
        assert abs(el.e - 0.436965) < 0.001
        assert abs(el.i - 43.869) < 0.1

    def test_sedna_elements(self):
        """Sedna elements should match JPL values."""
        el = MINOR_BODY_ELEMENTS[SEDNA]
        assert el.name == "Sedna"
        assert abs(el.a - 549.541) < 1.0  # Large semi-major axis
        assert abs(el.e - 0.861297) < 0.001  # Highly eccentric
        assert abs(el.i - 11.926) < 0.1


class TestPositionCalculation:
    """Test position calculations using updated elements."""

    def test_position_at_epoch(self, epoch_2025, all_minor_body_ids):
        """All bodies should return valid positions at epoch."""
        for body_id in all_minor_body_ids:
            pos, _ = ephem.calc_ut(epoch_2025, body_id, FLG_HELCTR)

            name = MINOR_BODY_ELEMENTS[body_id].name
            assert 0 <= pos[0] < 360, f"{name}: invalid longitude {pos[0]}"
            assert -90 <= pos[1] <= 90, f"{name}: invalid latitude {pos[1]}"
            assert pos[2] > 0, f"{name}: invalid distance {pos[2]}"

    def test_position_near_epoch(self, epoch_2025, all_minor_body_ids):
        """Positions near epoch should be consistent."""
        for body_id in all_minor_body_ids:
            # Test at epoch and 1 day later
            pos1, _ = ephem.calc_ut(epoch_2025, body_id, FLG_HELCTR)
            pos2, _ = ephem.calc_ut(epoch_2025 + 1.0, body_id, FLG_HELCTR)

            name = MINOR_BODY_ELEMENTS[body_id].name

            # Position should change slightly but reasonably
            lon_diff = abs(pos2[0] - pos1[0])
            if lon_diff > 180:
                lon_diff = 360 - lon_diff

            # Maximum reasonable motion is ~1 degree/day for main belt asteroids
            assert lon_diff < 2.0, f"{name}: unreasonable motion {lon_diff}deg/day"

    @pytest.mark.parametrize(
        "body_id,min_dist,max_dist",
        [
            (CHIRON, 8.0, 19.0),  # Perihelion ~8.5 AU, Aphelion ~18.9 AU
            (CERES, 2.5, 3.0),  # Main belt
            (VESTA, 2.1, 2.6),  # Main belt
            (ERIS, 38.0, 98.0),  # Trans-Neptunian
            (SEDNA, 76.0, 1000.0),  # Extreme orbit
        ],
    )
    def test_distance_in_expected_range(self, epoch_2025, body_id, min_dist, max_dist):
        """Heliocentric distances should be within orbital constraints."""
        pos, _ = ephem.calc_ut(epoch_2025, body_id, FLG_HELCTR)
        dist = pos[2]

        name = MINOR_BODY_ELEMENTS[body_id].name
        assert min_dist <= dist <= max_dist, (
            f"{name}: distance {dist:.2f} AU outside range [{min_dist}, {max_dist}]"
        )


class TestPeriodsAndOrbits:
    """Test orbital period calculations."""

    @pytest.mark.parametrize(
        "body_id,name,expected_period_years,tolerance_years",
        [
            (CHIRON, "Chiron", 51, 2),
            (PHOLUS, "Pholus", 91, 3),
            (CERES, "Ceres", 4.6, 0.2),
            (VESTA, "Vesta", 3.6, 0.2),
            (ERIS, "Eris", 561, 10),
            (HAUMEA, "Haumea", 282, 5),
            (MAKEMAKE, "Makemake", 307, 5),
        ],
    )
    def test_orbital_period(
        self, body_id, name, expected_period_years, tolerance_years
    ):
        """Orbital periods should match expected values."""
        elements = MINOR_BODY_ELEMENTS[body_id]

        # Period = 360 / n (in days), then convert to years
        period_days = 360.0 / elements.n
        period_years = period_days / DAYS_PER_YEAR

        assert abs(period_years - expected_period_years) < tolerance_years, (
            f"{name}: period {period_years:.1f} years differs from "
            f"expected {expected_period_years} years"
        )
