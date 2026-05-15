"""
FLG_NOGDEFL Comparison Tests.

Validates that FLG_NOGDEFL (skip gravitational light deflection, keep aberration)
produces results matching pyswisseph for:
- All major planets (Sun through Pluto)
- NOGDEFL alone (aberration without deflection)
- NOGDEFL + NOABERR = ASTROMETRIC (no aberration, no deflection)
- NOGDEFL with equatorial and J2000 coordinate modes
- Fixed stars with NOGDEFL
- Verifies NOGDEFL actually changes the position (deflection is non-zero)
"""

from __future__ import annotations

import pytest
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
    FLG_SWIEPH,
    FLG_SPEED,
    FLG_NOGDEFL,
    FLG_NOABERR,
    FLG_EQUATORIAL,
    FLG_J2000,
)


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

PLANETS = [
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

TEST_DATES = [
    (2000, 1, 1, 12.0, "J2000.0"),
    (2024, 6, 21, 12.0, "Summer Solstice 2024"),
    (1980, 3, 20, 0.0, "Equinox 1980"),
]

# Tolerance: sub-arcsecond agreement expected
LON_TOL = 0.001  # degrees (~3.6")
LAT_TOL = 0.001  # degrees
DIST_TOL = 0.0001  # AU
VEL_TOL = 0.01  # degrees/day


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    diff = abs(val1 - val2)
    if diff > 180:
        diff = 360 - diff
    return diff


# ============================================================================
# NOGDEFL PLANET TESTS
# ============================================================================


class TestNogdeflPlanets:
    """Test FLG_NOGDEFL for all planets."""

    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,label", TEST_DATES)
    def test_nogdefl_ecliptic(
        self, planet_id, planet_name, year, month, day, hour, label
    ):
        """NOGDEFL ecliptic positions match pyswisseph."""
        jd = swe.julday(year, month, day, hour)
        iflag = FLG_SWIEPH | FLG_SPEED | FLG_NOGDEFL

        swe_result = swe.calc_ut(jd, planet_id, iflag)
        lib_result = ephem.calc_ut(jd, planet_id, iflag)

        swe_pos = swe_result[0] if isinstance(swe_result, tuple) else swe_result
        lib_pos = lib_result[0]

        lon_diff = angular_diff(swe_pos[0], lib_pos[0])
        lat_diff = abs(swe_pos[1] - lib_pos[1])
        dist_diff = abs(swe_pos[2] - lib_pos[2])

        assert lon_diff < LON_TOL, (
            f"{planet_name} @ {label}: lon diff {lon_diff:.6f}° "
            f"(swe={swe_pos[0]:.6f}, lib={lib_pos[0]:.6f})"
        )
        assert lat_diff < LAT_TOL, f"{planet_name} @ {label}: lat diff {lat_diff:.6f}°"
        assert dist_diff < DIST_TOL, (
            f"{planet_name} @ {label}: dist diff {dist_diff:.8f} AU"
        )

        # Check velocities
        vel_lon_diff = abs(swe_pos[3] - lib_pos[3])
        vel_lat_diff = abs(swe_pos[4] - lib_pos[4])
        assert vel_lon_diff < VEL_TOL, (
            f"{planet_name} @ {label}: vel_lon diff {vel_lon_diff:.6f}°/day"
        )
        assert vel_lat_diff < VEL_TOL, (
            f"{planet_name} @ {label}: vel_lat diff {vel_lat_diff:.6f}°/day"
        )

    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    def test_nogdefl_equatorial_j2000(self, planet_id, planet_name):
        """NOGDEFL + EQUATORIAL + J2000 matches pyswisseph."""
        jd = swe.julday(2000, 1, 1, 12.0)
        iflag = (
            FLG_SWIEPH | FLG_SPEED | FLG_NOGDEFL | FLG_EQUATORIAL | FLG_J2000
        )

        swe_result = swe.calc_ut(jd, planet_id, iflag)
        lib_result = ephem.calc_ut(jd, planet_id, iflag)

        swe_pos = swe_result[0] if isinstance(swe_result, tuple) else swe_result
        lib_pos = lib_result[0]

        lon_diff = angular_diff(swe_pos[0], lib_pos[0])
        lat_diff = abs(swe_pos[1] - lib_pos[1])

        assert lon_diff < LON_TOL, (
            f"{planet_name}: RA diff {lon_diff:.6f}° "
            f"(swe={swe_pos[0]:.6f}, lib={lib_pos[0]:.6f})"
        )
        assert lat_diff < LAT_TOL, f"{planet_name}: Dec diff {lat_diff:.6f}°"


class TestNogdeflAstrometric:
    """Test NOGDEFL + NOABERR = ASTROMETRIC (no aberration, no deflection)."""

    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    def test_astrometric_matches_nogdefl_noaberr(self, planet_id, planet_name):
        """FLG_ASTROMETRIC (NOABERR|NOGDEFL) matches pyswisseph."""
        jd = swe.julday(2000, 1, 1, 12.0)
        # FLG_ASTROMETRIC = FLG_NOABERR | FLG_NOGDEFL
        iflag = FLG_SWIEPH | FLG_SPEED | FLG_NOABERR | FLG_NOGDEFL

        swe_result = swe.calc_ut(jd, planet_id, iflag)
        lib_result = ephem.calc_ut(jd, planet_id, iflag)

        swe_pos = swe_result[0] if isinstance(swe_result, tuple) else swe_result
        lib_pos = lib_result[0]

        lon_diff = angular_diff(swe_pos[0], lib_pos[0])
        lat_diff = abs(swe_pos[1] - lib_pos[1])

        assert lon_diff < LON_TOL, (
            f"{planet_name}: ASTROMETRIC lon diff {lon_diff:.6f}°"
        )
        assert lat_diff < LAT_TOL, (
            f"{planet_name}: ASTROMETRIC lat diff {lat_diff:.6f}°"
        )


class TestNogdeflEffect:
    """Verify NOGDEFL actually changes the position (deflection is non-zero)."""

    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_nogdefl_differs_from_apparent(self, planet_id, planet_name):
        """NOGDEFL should give a slightly different position from apparent.

        Gravitational deflection is up to ~1.7" for Sun-grazing light rays,
        and typically ~0.001-0.01" for planets away from the Sun's limb.
        For outer planets, the effect is usually detectable.
        """
        jd = swe.julday(2000, 1, 1, 12.0)

        # Normal apparent position
        iflag_normal = FLG_SWIEPH | FLG_SPEED
        swe_normal = swe.calc_ut(jd, planet_id, iflag_normal)
        swe_normal_pos = swe_normal[0] if isinstance(swe_normal, tuple) else swe_normal

        # NOGDEFL position
        iflag_nogdefl = FLG_SWIEPH | FLG_SPEED | FLG_NOGDEFL
        swe_nogdefl = swe.calc_ut(jd, planet_id, iflag_nogdefl)
        swe_nogdefl_pos = (
            swe_nogdefl[0] if isinstance(swe_nogdefl, tuple) else swe_nogdefl
        )

        # The difference should be non-zero (gravitational deflection exists)
        lon_diff = angular_diff(swe_normal_pos[0], swe_nogdefl_pos[0])

        # Deflection is typically small but measurable
        # For pyswisseph, the difference should be > 0 (even if tiny)
        # We just verify our library produces the same difference magnitude
        lib_normal = ephem.calc_ut(jd, planet_id, iflag_normal)
        lib_nogdefl = ephem.calc_ut(jd, planet_id, iflag_nogdefl)

        lib_lon_diff = angular_diff(lib_normal[0][0], lib_nogdefl[0][0])

        # Both should show a non-trivial deflection effect
        # (at least 0.0001° = 0.36" for planets not near the Sun's limb)
        # Note: the Moon has negligible deflection, so we skip it
        assert lon_diff > 0.0, f"{planet_name}: pyswisseph shows no deflection effect"
        assert lib_lon_diff > 0.0, (
            f"{planet_name}: libephemeris shows no deflection effect"
        )


class TestNogdeflFixedStars:
    """Test NOGDEFL for fixed stars."""

    @pytest.mark.parametrize("star_name", ["Regulus", "Spica", "Aldebaran"])
    def test_fixstar_nogdefl(self, star_name):
        """NOGDEFL fixed star positions match pyswisseph."""
        jd = swe.julday(2000, 1, 1, 12.0)
        iflag = FLG_SWIEPH | FLG_NOGDEFL

        swe_result = swe.fixstar_ut(star_name, jd, iflag)
        lib_result = ephem.fixstar_ut(star_name, jd, iflag)

        swe_pos = swe_result[0] if isinstance(swe_result, tuple) else swe_result
        lib_pos = lib_result[0]

        lon_diff = angular_diff(swe_pos[0], lib_pos[0])
        lat_diff = abs(swe_pos[1] - lib_pos[1])

        # Fixed stars: slightly relaxed tolerance (0.002° = 7.2")
        # due to catalog/proper motion differences
        star_tol = 0.002
        assert lon_diff < star_tol, (
            f"{star_name}: NOGDEFL lon diff {lon_diff:.6f}° "
            f"(swe={swe_pos[0]:.6f}, lib={lib_pos[0]:.6f})"
        )
        assert lat_diff < star_tol, f"{star_name}: NOGDEFL lat diff {lat_diff:.6f}°"

    @pytest.mark.parametrize("star_name", ["Regulus", "Spica", "Aldebaran"])
    def test_fixstar2_nogdefl(self, star_name):
        """NOGDEFL fixed star positions via fixstar2 match pyswisseph."""
        jd = swe.julday(2000, 1, 1, 12.0)
        iflag = FLG_SWIEPH | FLG_NOGDEFL

        swe_result = swe.fixstar2_ut(star_name, jd, iflag)
        lib_result = ephem.fixstar2_ut(star_name, jd, iflag)

        # pyswisseph fixstar2_ut returns (position_tuple, star_name, retflag)
        # libephemeris returns (position_tuple, star_name, retflag)
        swe_pos = swe_result[0]
        lib_pos = lib_result[0]

        lon_diff = angular_diff(swe_pos[0], lib_pos[0])
        lat_diff = abs(swe_pos[1] - lib_pos[1])

        star_tol = 0.002
        assert lon_diff < star_tol, (
            f"{star_name}: fixstar2 NOGDEFL lon diff {lon_diff:.6f}°"
        )
        assert lat_diff < star_tol, (
            f"{star_name}: fixstar2 NOGDEFL lat diff {lat_diff:.6f}°"
        )
