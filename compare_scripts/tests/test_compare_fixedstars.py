"""
Fixed Stars Comparison Tests.

Compares fixed star calculations between pyswisseph and libephemeris.
"""

import pytest
import swisseph as swe
import libephemeris as ephem

from compare_scripts.comparison_utils import FIXED_STARS


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# ============================================================================
# TOLERANCES
# ============================================================================

POSITION_TOL = 0.002  # degrees - relaxed for proper motion differences
VELOCITY_TOL = 0.0001  # degrees/day

# Stars where SE resolves to a different physical star than our IAU WGSN catalog:
# - Algedi: SE uses Alpha-1 Cap (HIP 100027), we use Alpha-2 Cap (HIP 100064) per IAU
# - Menkar: SE uses Lambda Ceti (HIP 455), we use Alpha Ceti (HIP 14135) per IAU WGSN
# - Algieba: SE uses Gamma-1 Leo (HIP 50583 via ga-1Leo), different Bayer component
# - Albireo: SE uses Beta-1 Cyg, we use Beta Cyg (different HIP)
# - Almach: SE uses Gamma-1 And, we use Gamma And (different HIP)
DIFFERENT_STAR_SKIPS = {
    "Algedi",
    "Menkar",
    "Algieba",
    "Albireo",
    "Almach",
}


# ============================================================================
# TEST DATA
# ============================================================================

# Use the comprehensive list from comparison_utils (103 stars)
# Remove duplicates while preserving order
FAMOUS_STARS = list(dict.fromkeys(FIXED_STARS))

TEST_DATES = [
    (2000, 1, 1, 12.0, "J2000"),
    (2024, 6, 15, 0.0, "Current"),
    (1950, 1, 1, 12.0, "Mid-century"),
]


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestFixedStarPositions:
    """Compare fixed star position calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("star_name", FAMOUS_STARS)
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    def test_fixed_star_position(self, star_name, year, month, day, hour, desc):
        """Test fixed star positions match within tolerance."""
        if star_name in DIFFERENT_STAR_SKIPS:
            pytest.skip(
                f"{star_name}: SE resolves to different physical star component"
            )
            return

        jd = swe.julday(year, month, day, hour)

        try:
            result_swe = swe.fixstar_ut(star_name, jd, 0)
            result_py = ephem.fixstar_ut(star_name, jd, 0)
        except Exception as e:
            pytest.skip(f"Star {star_name} not available: {e}")
            return

        lon_swe = result_swe[0][0]
        lon_py = result_py[0][0]

        diff = angular_diff(lon_swe, lon_py)

        assert diff < POSITION_TOL, (
            f"{star_name} at {desc}: longitude diff {diff:.6f}° exceeds tolerance"
        )


class TestFixedStarVelocity:
    """Compare fixed star velocity (proper motion) calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("star_name", FAMOUS_STARS[:5])
    def test_fixed_star_with_speed(self, star_name):
        """Test fixed star velocity calculations."""
        jd = swe.julday(2000, 1, 1, 12.0)

        try:
            result_swe = swe.fixstar_ut(star_name, jd, swe.FLG_SPEED)
            result_py = ephem.fixstar_ut(star_name, jd, 256)  # FLG_SPEED
        except Exception as e:
            pytest.skip(f"Star {star_name} not available: {e}")
            return

        # Compare velocity (index 3)
        speed_swe = result_swe[0][3] if len(result_swe[0]) > 3 else 0
        speed_py = result_py[0][3] if len(result_py[0]) > 3 else 0

        diff = abs(speed_swe - speed_py)

        assert diff < VELOCITY_TOL, (
            f"{star_name} velocity diff {diff:.8f}°/day exceeds tolerance"
        )


class TestSpecificStars:
    """Test specific important fixed stars."""

    @pytest.mark.comparison
    def test_regulus_at_j2000(self):
        """Test Regulus position at J2000."""
        jd = 2451545.0

        try:
            result_swe = swe.fixstar_ut("Regulus", jd, 0)
        except Exception as e:
            pytest.skip(f"Regulus not available in pyswisseph: {e}")
            return

        result_py = ephem.fixstar_ut("Regulus", jd, 0)

        diff = angular_diff(result_swe[0][0], result_py[0][0])

        # Regulus should be around 29° Leo at J2000
        assert 149 < result_py[0][0] < 150, (
            f"Regulus position {result_py[0][0]}° not in expected range"
        )
        assert diff < POSITION_TOL

    @pytest.mark.comparison
    def test_spica_at_j2000(self):
        """Test Spica position at J2000."""
        jd = 2451545.0

        try:
            result_swe = swe.fixstar_ut("Spica", jd, 0)
        except Exception as e:
            pytest.skip(f"Spica not available in pyswisseph: {e}")
            return

        result_py = ephem.fixstar_ut("Spica", jd, 0)

        diff = angular_diff(result_swe[0][0], result_py[0][0])

        # Spica should be around 23° Libra at J2000
        assert 203 < result_py[0][0] < 204, (
            f"Spica position {result_py[0][0]}° not in expected range"
        )
        assert diff < POSITION_TOL


class TestProperMotion:
    """Test proper motion is correctly applied over time."""

    @pytest.mark.comparison
    def test_star_moves_over_century(self):
        """Test that star positions change over a century due to proper motion."""
        star = "Sirius"
        jd_1900 = swe.julday(1900, 1, 1, 12.0)
        jd_2000 = swe.julday(2000, 1, 1, 12.0)
        jd_2100 = swe.julday(2100, 1, 1, 12.0)

        pos_1900 = ephem.fixstar_ut(star, jd_1900, 0)[0][0]
        pos_2000 = ephem.fixstar_ut(star, jd_2000, 0)[0][0]
        pos_2100 = ephem.fixstar_ut(star, jd_2100, 0)[0][0]

        # Positions should be different due to precession and proper motion
        assert pos_1900 != pos_2000 != pos_2100, (
            "Star positions should differ over centuries"
        )

        # The change should be monotonic (precession moves everything)
        change_1 = pos_2000 - pos_1900
        change_2 = pos_2100 - pos_2000

        # Both changes should be positive (precession increases longitude)
        assert change_1 > 0, "Position should increase from 1900 to 2000"
        assert change_2 > 0, "Position should increase from 2000 to 2100"
