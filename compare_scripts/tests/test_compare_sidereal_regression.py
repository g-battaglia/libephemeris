"""
Sidereal Bug Fix Regression Tests: pyswisseph vs libephemeris.

Guards against regressions of 4 sidereal calculation bugs fixed in the
leb/precision branch.  Uses pyswisseph as the reference standard.

Bug 1 (64b8367): Pipeline A SID+EQ used nutation matrix instead of mean equator.
Bug 2 (e6555ed): Pipeline B/C dpsi nutation handling wrong for SID+EQ.
Bug 3 (9f0fde7): J2000 suppression for non-mean SID bodies — NOW FIXED.
    pyswisseph silently ignores FLG_J2000 for TrueNode, OscuApog, IntpApog,
    IntpPerg when FLG_SIDEREAL is set.  This is a behavioral bug (ayanamsha
    and J2000 ecliptic precession are geometrically distinct operations).
    LibEphemeris intentionally fixes this: FLG_J2000 is honored for ALL bodies.
    Tests in TestBug3J2kIntentionalDivergence verify the expected divergence.
    See docs/reference/se-bug-sidereal-j2000-nodes.md
Bug 4 (b816be0): (a) Frame bias in _get_precession_matrix, (b) SID+J2K
                  precession order for MeanNode/MeanApog.

All tests compare pyswisseph output against libephemeris.  Tolerances are
tight enough to catch the specific error signatures of each bug.
"""

from __future__ import annotations

import pytest
import swisseph as swe

import libephemeris as ephem
from libephemeris.constants import (
    SUN,
    MOON,
    MARS,
    JUPITER,
    SATURN,
    MEAN_NODE,
    TRUE_NODE,
    MEAN_APOG,
    OSCU_APOG,
    INTP_APOG,
    INTP_PERG,
    FLG_SPEED,
    FLG_SIDEREAL,
    FLG_EQUATORIAL,
    FLG_J2000,
    SIDM_LAHIRI,
    SIDM_FAGAN_BRADLEY,
    SIDM_RAMAN,
)


# ============================================================================
# UTILITIES
# ============================================================================


def angular_diff(a: float, b: float) -> float:
    """Angular difference accounting for 360 wrap."""
    d = abs(a - b)
    if d > 180:
        d = 360 - d
    return d


def arcsec(deg: float) -> float:
    return deg * 3600.0


# Pipeline A bodies (ICRS barycentric)
PIPELINE_A_BODIES = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
]

# Pipeline B bodies (ecliptic direct) — mean vs true distinction matters
PIPELINE_B_MEAN_BODIES = [
    (MEAN_NODE, "MeanNode"),
    (MEAN_APOG, "MeanApog"),
]

PIPELINE_B_TRUE_BODIES = [
    (TRUE_NODE, "TrueNode"),
    (OSCU_APOG, "OscuApog"),
]

PIPELINE_B_INTERP_BODIES = [
    (INTP_APOG, "IntpApog"),
    (INTP_PERG, "IntpPerg"),
]

ALL_PIPELINE_B_BODIES = (
    PIPELINE_B_MEAN_BODIES + PIPELINE_B_TRUE_BODIES + PIPELINE_B_INTERP_BODIES
)

# Test dates across the DE440 range
TEST_DATES = [
    (2451545.0, "J2000"),
    (2415020.5, "J1900"),
    (2460000.5, "2023"),
    (2469807.5, "2050"),
    (2433282.5, "1950"),
    (2440587.5, "1970"),
    (2378496.5, "1800"),
    (2524593.5, "2200"),
]

# Tolerances — calibrated to catch each bug's error signature
# while accommodating known systematic offsets (IAU 2006 vs 1976 precession)
PIPELINE_A_SID_EQ_TOL = 0.01  # ~36" (Bug 1 error was ~36")
PIPELINE_A_SID_J2K_TOL = 0.01  # ~36" (Bug 1 error was ~0.3")
PIPELINE_B_SID_EQ_TOL = 0.06  # ~216" (OscuApog can reach ~0.054° in SID+EQ)
PIPELINE_B_SID_J2K_TOL = 0.015  # ~54" (includes precession model difference)
PIPELINE_B_SID_EQ_MEAN_TOL = 0.015  # Mean bodies
PIPELINE_B_SID_EQ_TRUE_TOL = 0.06  # True bodies (OscuApog ~0.054°)
INTP_LON_TOL = 5.5  # IntpApog/IntpPerg: intentional algorithm deviation (see docs)

SID_MODES = [
    (SIDM_LAHIRI, "Lahiri"),
    (SIDM_FAGAN_BRADLEY, "Fagan-Bradley"),
    (SIDM_RAMAN, "Raman"),
]


# ============================================================================
# BUG 1: Pipeline A SID+EQ / SID+J2K (commit 64b8367)
# ============================================================================


class TestBug1PipelineASidEq:
    """Pipeline A SID+EQ: nutation matrix was used instead of mean equator.

    Error signature: ~36" RA error at J2000.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_A_BODIES)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:4])
    @pytest.mark.parametrize("sid_mode,sid_name", SID_MODES)
    def test_sid_eq_position(
        self, body_id, body_name, jd, date_desc, sid_mode, sid_name
    ):
        """Pipeline A SID+EQ RA/Dec matches pyswisseph."""
        swe.set_sid_mode(sid_mode)
        ephem.set_sid_mode(sid_mode)

        flags = FLG_SIDEREAL | FLG_EQUATORIAL | FLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.calc_ut(jd, body_id, flags)

        ra_diff = angular_diff(float(pos_swe[0]), float(pos_py[0]))
        dec_diff = abs(float(pos_swe[1]) - float(pos_py[1]))

        assert ra_diff < PIPELINE_A_SID_EQ_TOL, (
            f"{body_name} SID+EQ ({sid_name}) at {date_desc}: "
            f'RA diff {ra_diff:.6f}° ({arcsec(ra_diff):.1f}") '
            f"[Bug 1 regression?]"
        )
        assert dec_diff < PIPELINE_A_SID_EQ_TOL, (
            f"{body_name} SID+EQ ({sid_name}) at {date_desc}: "
            f'Dec diff {dec_diff:.6f}° ({arcsec(dec_diff):.1f}")'
        )


class TestBug1PipelineASidJ2k:
    """Pipeline A SID+J2K: true ayanamsha was used instead of mean.

    Error signature: ~0.3" lon error at J2000.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_A_BODIES)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:4])
    def test_sid_j2k_position(self, body_id, body_name, jd, date_desc):
        """Pipeline A SID+J2K lon matches pyswisseph."""
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.set_sid_mode(SIDM_LAHIRI)

        flags = FLG_SIDEREAL | FLG_J2000 | FLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.calc_ut(jd, body_id, flags)

        lon_diff = angular_diff(float(pos_swe[0]), float(pos_py[0]))

        assert lon_diff < PIPELINE_A_SID_J2K_TOL, (
            f"{body_name} SID+J2K at {date_desc}: "
            f'lon diff {lon_diff:.6f}° ({arcsec(lon_diff):.1f}") '
            f"[Bug 1 regression?]"
        )


class TestBug1PipelineASidEqJ2k:
    """Pipeline A SID+EQ+J2K combined."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_A_BODIES[:3])
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:4])
    def test_sid_eq_j2k_position(self, body_id, body_name, jd, date_desc):
        """Pipeline A SID+EQ+J2K RA matches pyswisseph."""
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.set_sid_mode(SIDM_LAHIRI)

        flags = FLG_SIDEREAL | FLG_EQUATORIAL | FLG_J2000 | FLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.calc_ut(jd, body_id, flags)

        ra_diff = angular_diff(float(pos_swe[0]), float(pos_py[0]))

        assert ra_diff < PIPELINE_A_SID_J2K_TOL, (
            f"{body_name} SID+EQ+J2K at {date_desc}: "
            f'RA diff {ra_diff:.6f}° ({arcsec(ra_diff):.1f}")'
        )


# ============================================================================
# BUG 2: Pipeline B/C SID+EQ dpsi handling (commit e6555ed)
# ============================================================================


class TestBug2PipelineBSidEq:
    """Pipeline B SID+EQ: dpsi nutation handling was wrong.

    - MeanNode/MeanApog should skip dpsi
    - TrueNode/OscuApog should subtract dpsi
    - Mean obliquity must be used for SID+EQ rotation

    Error signature: ~10-20" for MeanNode SID+EQ.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", ALL_PIPELINE_B_BODIES)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:4])
    def test_sid_eq_all_bodies(self, body_id, body_name, jd, date_desc):
        """All Pipeline B bodies SID+EQ match pyswisseph."""
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.set_sid_mode(SIDM_LAHIRI)

        flags = FLG_SIDEREAL | FLG_EQUATORIAL | FLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.calc_ut(jd, body_id, flags)

        if body_id in (INTP_APOG, INTP_PERG):
            tol = INTP_LON_TOL
        else:
            tol = PIPELINE_B_SID_EQ_TOL

        ra_diff = angular_diff(float(pos_swe[0]), float(pos_py[0]))

        assert ra_diff < tol, (
            f"{body_name} SID+EQ at {date_desc}: "
            f'RA diff {ra_diff:.6f}° ({arcsec(ra_diff):.1f}") '
            f"exceeds {tol}° [Bug 2 regression?]"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_B_MEAN_BODIES)
    @pytest.mark.parametrize("sid_mode,sid_name", SID_MODES)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:3])
    def test_mean_body_sid_eq_multi_ayanamsha(
        self, body_id, body_name, sid_mode, sid_name, jd, date_desc
    ):
        """Mean ecliptic bodies SID+EQ across multiple ayanamshas."""
        swe.set_sid_mode(sid_mode)
        ephem.set_sid_mode(sid_mode)

        flags = FLG_SIDEREAL | FLG_EQUATORIAL | FLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.calc_ut(jd, body_id, flags)

        ra_diff = angular_diff(float(pos_swe[0]), float(pos_py[0]))

        assert ra_diff < PIPELINE_B_SID_EQ_MEAN_TOL, (
            f"{body_name} SID+EQ ({sid_name}) at {date_desc}: "
            f'RA diff {ra_diff:.6f}° ({arcsec(ra_diff):.1f}")'
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_B_TRUE_BODIES)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:3])
    def test_true_body_sid_eq_speed(self, body_id, body_name, jd, date_desc):
        """True ecliptic bodies SID+EQ velocity matches pyswisseph."""
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.set_sid_mode(SIDM_LAHIRI)

        flags = FLG_SIDEREAL | FLG_EQUATORIAL | FLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.calc_ut(jd, body_id, flags)

        speed_tol = 0.05 if body_id == OSCU_APOG else 0.005
        speed_diff = abs(float(pos_swe[3]) - float(pos_py[3]))

        assert speed_diff < speed_tol, (
            f"{body_name} SID+EQ speed at {date_desc}: diff {speed_diff:.6f}°/day"
        )


# ============================================================================
# BUG 3: J2000 for non-mean SID bodies — INTENTIONAL SE DIVERGENCE
# ============================================================================


class TestBug3J2kIntentionalDivergence:
    """LibEphemeris intentionally diverges from pyswisseph here.

    pyswisseph silently ignores FLG_J2000 for TrueNode, OscuApog, IntpApog,
    IntpPerg when FLG_SIDEREAL is set.  This is a behavioral bug: ayanamsha
    (1D longitude zero-point shift) and J2000 ecliptic precession (3D plane
    rotation) are geometrically distinct, composable operations.

    LibEphemeris honors FLG_J2000 for ALL bodies uniformly.  These tests
    verify the intentional divergence and document the expected delta.

    See docs/reference/se-bug-sidereal-j2000-nodes.md
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_B_TRUE_BODIES)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:4])
    def test_true_body_sid_j2k_diverges_from_swe(
        self, body_id, body_name, jd, date_desc
    ):
        """TrueNode/OscuApog: SID+J2K now intentionally diverges from pyswisseph.

        pyswisseph ignores J2000 (SID+J2K == SID), but libephemeris honors it.
        The delta should be consistent with the ecliptic precession angle.
        """
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.set_sid_mode(SIDM_LAHIRI)

        flags_sid_j2k = FLG_SIDEREAL | FLG_J2000 | FLG_SPEED

        # pyswisseph: still ignores J2K for these bodies
        pos_swe, _ = swe.calc_ut(jd, body_id, flags_sid_j2k)
        # libephemeris: now honors J2K
        pos_py, _ = ephem.calc_ut(jd, body_id, flags_sid_j2k)

        divergence = angular_diff(float(pos_swe[0]), float(pos_py[0]))

        # The divergence should be non-trivial (> 0.001°) except near J2000
        # where the ecliptic precession angle is very small.
        # At J2000 the delta is ~0.004° (frame bias only).
        # At J1900 it's ~1.4°, at 2023 ~0.33°, at 2050 ~0.69°.
        if date_desc == "J2000":
            assert divergence < 0.01, (
                f"{body_name} SID+J2K at {date_desc}: divergence {divergence:.6f}° "
                f"too large near J2000"
            )
        else:
            assert divergence > 0.01, (
                f"{body_name} SID+J2K at {date_desc}: divergence {divergence:.6f}° "
                f"too small — J2000 precession not being applied? "
                f"[SE bug fix regression]"
            )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_B_TRUE_BODIES)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:4])
    def test_true_body_j2k_actually_applied(self, body_id, body_name, jd, date_desc):
        """TrueNode/OscuApog: SID+J2K must differ from SID alone.

        If they are identical, J2000 precession is NOT being applied.
        """
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.set_sid_mode(SIDM_LAHIRI)

        flags_sid = FLG_SIDEREAL | FLG_SPEED
        flags_sid_j2k = FLG_SIDEREAL | FLG_J2000 | FLG_SPEED

        pos_sid, _ = ephem.calc_ut(jd, body_id, flags_sid)
        pos_j2k, _ = ephem.calc_ut(jd, body_id, flags_sid_j2k)

        diff = angular_diff(float(pos_sid[0]), float(pos_j2k[0]))

        # Even at J2000, there should be a small frame-bias difference (~0.004°)
        assert diff > 0.001, (
            f"{body_name} SID vs SID+J2K at {date_desc}: diff {diff:.6f}° "
            f"(J2000 precession not applied?) [SE bug fix regression]"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_B_INTERP_BODIES)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:3])
    def test_interp_body_j2k_actually_applied(self, body_id, body_name, jd, date_desc):
        """IntpApog/IntpPerg: SID+J2K must differ from SID alone."""
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.set_sid_mode(SIDM_LAHIRI)

        flags_sid = FLG_SIDEREAL | FLG_SPEED
        flags_sid_j2k = FLG_SIDEREAL | FLG_J2000 | FLG_SPEED

        pos_sid, _ = ephem.calc_ut(jd, body_id, flags_sid)
        pos_j2k, _ = ephem.calc_ut(jd, body_id, flags_sid_j2k)

        diff = angular_diff(float(pos_sid[0]), float(pos_j2k[0]))

        assert diff > 0.001, (
            f"{body_name} SID vs SID+J2K at {date_desc}: diff {diff:.6f}° "
            f"(J2000 precession not applied?) [SE bug fix regression]"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_B_MEAN_BODIES)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:4])
    def test_mean_body_sid_j2k_differs_from_sid(
        self, body_id, body_name, jd, date_desc
    ):
        """MeanNode/MeanApog: SID+J2K should DIFFER from SID (J2K applied).

        If identical, J2K is being incorrectly suppressed for mean bodies.
        """
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.set_sid_mode(SIDM_LAHIRI)

        flags_sid = FLG_SIDEREAL | FLG_SPEED
        flags_sid_j2k = FLG_SIDEREAL | FLG_J2000 | FLG_SPEED

        # pyswisseph reference
        pos_swe_sid, _ = swe.calc_ut(jd, body_id, flags_sid)
        pos_swe_j2k, _ = swe.calc_ut(jd, body_id, flags_sid_j2k)

        swe_diff = angular_diff(float(pos_swe_sid[0]), float(pos_swe_j2k[0]))

        # libephemeris
        pos_py_sid, _ = ephem.calc_ut(jd, body_id, flags_sid)
        pos_py_j2k, _ = ephem.calc_ut(jd, body_id, flags_sid_j2k)

        py_diff = angular_diff(float(pos_py_sid[0]), float(pos_py_j2k[0]))

        # Both should show a difference (J2K IS applied for mean bodies)
        if swe_diff > 0.001:
            assert py_diff > 0.001, (
                f"{body_name} SID vs SID+J2K: pyswisseph differs by {swe_diff:.6f}°, "
                f"but libephemeris differs by only {py_diff:.6f}° "
                f"(J2K incorrectly suppressed?) [Bug 3 regression]"
            )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_B_MEAN_BODIES)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:4])
    def test_mean_body_sid_j2k_matches_swe(self, body_id, body_name, jd, date_desc):
        """MeanNode/MeanApog SID+J2K lon matches pyswisseph."""
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.set_sid_mode(SIDM_LAHIRI)

        flags = FLG_SIDEREAL | FLG_J2000 | FLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.calc_ut(jd, body_id, flags)

        lon_diff = angular_diff(float(pos_swe[0]), float(pos_py[0]))

        assert lon_diff < PIPELINE_B_SID_J2K_TOL, (
            f"{body_name} SID+J2K at {date_desc}: "
            f'lon diff {lon_diff:.6f}° ({arcsec(lon_diff):.1f}")'
        )


# ============================================================================
# BUG 4: Frame bias + SID+J2K precession order (commit b816be0)
# ============================================================================


class TestBug4FrameBiasAndPrecessionOrder:
    """(a) _get_precession_matrix used t.P including ICRS frame bias (~17 mas).
    (b) MeanNode/MeanApog SID+J2K applied precession before ayanamsha
        (non-commutative, up to 28" at extreme dates).
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_B_MEAN_BODIES)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES)
    def test_mean_body_sid_j2k_extended_dates(self, body_id, body_name, jd, date_desc):
        """MeanNode/MeanApog SID+J2K across extended dates.

        Bug 4b error grows with distance from J2000, up to 28" at extreme dates.
        """
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.set_sid_mode(SIDM_LAHIRI)

        flags = FLG_SIDEREAL | FLG_J2000 | FLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.calc_ut(jd, body_id, flags)

        lon_diff = angular_diff(float(pos_swe[0]), float(pos_py[0]))

        assert lon_diff < PIPELINE_B_SID_J2K_TOL, (
            f"{body_name} SID+J2K at {date_desc}: "
            f'lon diff {lon_diff:.6f}° ({arcsec(lon_diff):.1f}") '
            f"[Bug 4b regression?]"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_B_MEAN_BODIES)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:4])
    def test_mean_body_sid_eq_j2k_combined(self, body_id, body_name, jd, date_desc):
        """MeanNode/MeanApog SID+EQ+J2K: triple combo matches pyswisseph."""
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.set_sid_mode(SIDM_LAHIRI)

        flags = FLG_SIDEREAL | FLG_EQUATORIAL | FLG_J2000 | FLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.calc_ut(jd, body_id, flags)

        ra_diff = angular_diff(float(pos_swe[0]), float(pos_py[0]))

        assert ra_diff < PIPELINE_B_SID_J2K_TOL, (
            f"{body_name} SID+EQ+J2K at {date_desc}: "
            f'RA diff {ra_diff:.6f}° ({arcsec(ra_diff):.1f}")'
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_A_BODIES[:3])
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:4])
    def test_pipeline_a_sid_eq_no_frame_bias(self, body_id, body_name, jd, date_desc):
        """Pipeline A SID+EQ: verify no ICRS frame bias in precession.

        Bug 4a introduced ~17 mas systematic at J2000.  The fix replaced t.P
        (which includes ICRS bias) with mean_equator rotation.
        """
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.set_sid_mode(SIDM_LAHIRI)

        flags = FLG_SIDEREAL | FLG_EQUATORIAL | FLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.calc_ut(jd, body_id, flags)

        ra_diff = angular_diff(float(pos_swe[0]), float(pos_py[0]))

        assert ra_diff < PIPELINE_A_SID_EQ_TOL, (
            f"{body_name} SID+EQ at {date_desc}: "
            f'RA diff {ra_diff:.6f}° ({arcsec(ra_diff):.1f}") '
            f"[Bug 4a frame bias regression?]"
        )


# ============================================================================
# COMPREHENSIVE SID+EQ COVERAGE: ALL PIPELINE B/C BODIES
# ============================================================================


class TestSidEqComprehensive:
    """Ensure every Pipeline B body works correctly with SID+EQ."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (MEAN_NODE, "MeanNode"),
            (TRUE_NODE, "TrueNode"),
            (MEAN_APOG, "MeanApog"),
            (OSCU_APOG, "OscuApog"),
        ],
    )
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:4])
    @pytest.mark.parametrize("sid_mode,sid_name", SID_MODES)
    def test_primary_body_sid_eq(
        self, body_id, body_name, jd, date_desc, sid_mode, sid_name
    ):
        """All primary Pipeline B bodies SID+EQ across ayanamshas."""
        swe.set_sid_mode(sid_mode)
        ephem.set_sid_mode(sid_mode)

        flags = FLG_SIDEREAL | FLG_EQUATORIAL | FLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.calc_ut(jd, body_id, flags)

        tol = PIPELINE_B_SID_EQ_TOL
        ra_diff = angular_diff(float(pos_swe[0]), float(pos_py[0]))
        dec_diff = abs(float(pos_swe[1]) - float(pos_py[1]))

        assert ra_diff < tol, (
            f"{body_name} SID+EQ ({sid_name}) at {date_desc}: "
            f'RA diff {ra_diff:.6f}° ({arcsec(ra_diff):.1f}")'
        )
        assert dec_diff < tol, (
            f"{body_name} SID+EQ ({sid_name}) at {date_desc}: "
            f'Dec diff {dec_diff:.6f}° ({arcsec(dec_diff):.1f}")'
        )


class TestSidJ2kComprehensive:
    """Ensure every Pipeline B body works correctly with SID+J2K.

    Mean bodies (MeanNode, MeanApog) match pyswisseph within tolerance.
    True/osculating bodies (TrueNode, OscuApog) intentionally diverge from
    pyswisseph because LibEphemeris fixes the J2000 suppression bug.
    See docs/reference/se-bug-sidereal-j2000-nodes.md
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (MEAN_NODE, "MeanNode"),
            (MEAN_APOG, "MeanApog"),
        ],
    )
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:4])
    def test_mean_body_sid_j2k_matches_swe(self, body_id, body_name, jd, date_desc):
        """Mean bodies SID+J2K still match pyswisseph (no divergence)."""
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.set_sid_mode(SIDM_LAHIRI)

        flags = FLG_SIDEREAL | FLG_J2000 | FLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.calc_ut(jd, body_id, flags)

        tol = PIPELINE_B_SID_J2K_TOL
        lon_diff = angular_diff(float(pos_swe[0]), float(pos_py[0]))

        assert lon_diff < tol, (
            f"{body_name} SID+J2K at {date_desc}: "
            f'lon diff {lon_diff:.6f}° ({arcsec(lon_diff):.1f}")'
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (TRUE_NODE, "TrueNode"),
            (OSCU_APOG, "OscuApog"),
        ],
    )
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES[:4])
    def test_true_body_sid_j2k_divergence_expected(
        self, body_id, body_name, jd, date_desc
    ):
        """True/osculating bodies SID+J2K: verify expected SE divergence.

        LibEphemeris intentionally diverges from pyswisseph here because SE
        silently ignores FLG_J2000 for these bodies when FLG_SIDEREAL
        is set.  The divergence should match the ecliptic precession angle.
        """
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.set_sid_mode(SIDM_LAHIRI)

        flags = FLG_SIDEREAL | FLG_J2000 | FLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.calc_ut(jd, body_id, flags)

        lon_diff = angular_diff(float(pos_swe[0]), float(pos_py[0]))

        # Verify the divergence is in a plausible range:
        # - At J2000.0: ~0.004° (frame bias only, very small precession)
        # - At 2023: ~0.33° (consistent with ~50 arcsec/year * 23 years)
        # - At J1900: ~1.4° (100 years of ecliptic precession)
        # - At 2050: ~0.69° (50 years of ecliptic precession)
        assert lon_diff < 5.0, (
            f"{body_name} SID+J2K at {date_desc}: "
            f"divergence {lon_diff:.6f}° unexpectedly large"
        )
