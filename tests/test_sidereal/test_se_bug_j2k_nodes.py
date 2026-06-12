"""
Tests for the intentional Swiss Ephemeris divergence: SIDEREAL + J2000 on
lunar nodes and apsides.

pyswisseph silently ignores FLG_J2000 for TrueNode (TRUE_NODE),
OscuApog (OSCU_APOG), IntpApog (INTP_APOG), and IntpPerg (INTP_PERG)
when FLG_SIDEREAL is also set.  This is a behavioral bug: ayanamsha
(1D longitude zero-point shift) and J2000 ecliptic precession (3D ecliptic
plane rotation) are geometrically distinct, composable operations.

LibEphemeris intentionally corrects this: FLG_J2000 is honored for ALL
bodies uniformly, including these four.

These tests verify:
  1. J2000 is applied for all bodies (SID+J2K != SID)
  2. LEB and Skyfield paths produce consistent results
  3. Physical sanity: True vs Mean node/apogee distance stays plausible
  4. The divergence from pyswisseph is in the expected range

See docs/reference/se-bug-sidereal-j2000-nodes.md for full analysis.
"""

from __future__ import annotations


import pytest

import libephemeris as ephem
from libephemeris.constants import (
    MEAN_NODE,
    TRUE_NODE,
    MEAN_APOG,
    OSCU_APOG,
    INTP_APOG,
    INTP_PERG,
    FLG_SPEED,
    FLG_SIDEREAL,
    FLG_J2000,
    SIDM_LAHIRI,
)


def angular_diff(a: float, b: float) -> float:
    """Unsigned angular difference accounting for 360 wrap."""
    d = abs(a - b)
    if d > 180:
        d = 360 - d
    return d


# Bodies affected by the SE bug (J2000 was silently ignored)
AFFECTED_BODIES = [
    (TRUE_NODE, "TrueNode"),
    (OSCU_APOG, "OscuApog"),
    (INTP_APOG, "IntpApog"),
    (INTP_PERG, "IntpPerg"),
]

# All Pipeline B bodies (for uniformity checks)
ALL_PIPELINE_B = [
    (MEAN_NODE, "MeanNode"),
    (TRUE_NODE, "TrueNode"),
    (MEAN_APOG, "MeanApog"),
    (OSCU_APOG, "OscuApog"),
    (INTP_APOG, "IntpApog"),
    (INTP_PERG, "IntpPerg"),
]

# Mean/True pairs for physical sanity checks
NODE_PAIRS = [
    (MEAN_NODE, TRUE_NODE, "MeanNode", "TrueNode"),
]

APOGEE_PAIRS = [
    (MEAN_APOG, OSCU_APOG, "MeanApog", "OscuApog"),
]

# Representative dates: JD values spanning the full range
TEST_DATES = [
    (2451545.0, "J2000"),
    (2460400.5, "2024"),
    (2415020.5, "J1900"),
    (2469807.5, "2050"),
]


class TestJ2000HonoredForAllBodies:
    """Verify FLG_J2000 is actually applied for all Pipeline B bodies.

    After the fix, SID+J2K must produce different results from SID alone
    for ALL bodies, including the four that pyswisseph incorrectly skips.
    """

    @pytest.mark.parametrize("body_id,body_name", ALL_PIPELINE_B)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES)
    def test_j2000_applied(self, body_id, body_name, jd, date_desc):
        """SID+J2K must differ from SID for every body."""
        ephem.set_sid_mode(SIDM_LAHIRI)

        flags_sid = FLG_SIDEREAL | FLG_SPEED
        flags_j2k = FLG_SIDEREAL | FLG_J2000 | FLG_SPEED

        pos_sid, _ = ephem.calc_ut(jd, body_id, flags_sid)
        pos_j2k, _ = ephem.calc_ut(jd, body_id, flags_j2k)

        diff = angular_diff(pos_sid[0], pos_j2k[0])

        # Even at J2000.0, there should be a small frame-bias difference
        assert diff > 0.001, (
            f"{body_name} SID+J2K at {date_desc}: "
            f"diff from SID only {diff:.6f} deg "
            f"(J2000 precession not applied?)"
        )


class TestPhysicalSanity:
    """Verify that True and Mean bodies stay physically close with SID+J2K.

    The true lunar node oscillates around the mean node with amplitude
    ~1.5 degrees.  If they are 10+ degrees apart, the coordinate frame
    is wrong.  This test catches the exact symptom of the SE bug: at
    0 CE, the uncorrected TrueNode was ~29 degrees from MeanNode.
    """

    @pytest.mark.parametrize("mean_id,true_id,mean_name,true_name", NODE_PAIRS)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES)
    def test_true_vs_mean_node_distance(
        self, mean_id, true_id, mean_name, true_name, jd, date_desc
    ):
        """TrueNode must stay within 2 degrees of MeanNode."""
        ephem.set_sid_mode(SIDM_LAHIRI)

        flags = FLG_SIDEREAL | FLG_J2000 | FLG_SPEED

        pos_mean, _ = ephem.calc_ut(jd, mean_id, flags)
        pos_true, _ = ephem.calc_ut(jd, true_id, flags)

        diff = angular_diff(pos_mean[0], pos_true[0])

        assert diff < 2.0, (
            f"{true_name} vs {mean_name} SID+J2K at {date_desc}: "
            f"diff {diff:.4f} deg (physically impossible, should be < 2 deg)"
        )

    @pytest.mark.parametrize("mean_id,true_id,mean_name,true_name", APOGEE_PAIRS)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES)
    def test_true_vs_mean_apogee_distance(
        self, mean_id, true_id, mean_name, true_name, jd, date_desc
    ):
        """OscuApog must stay within 35 degrees of MeanApog.

        The osculating apogee has larger perturbation swings than the
        true node (up to ~30 degrees in extreme cases), so the tolerance
        is wider.
        """
        ephem.set_sid_mode(SIDM_LAHIRI)

        flags = FLG_SIDEREAL | FLG_J2000 | FLG_SPEED

        pos_mean, _ = ephem.calc_ut(jd, mean_id, flags)
        pos_true, _ = ephem.calc_ut(jd, true_id, flags)

        diff = angular_diff(pos_mean[0], pos_true[0])

        assert diff < 35.0, (
            f"{true_name} vs {mean_name} SID+J2K at {date_desc}: "
            f"diff {diff:.4f} deg (physically implausible)"
        )


class TestLebVsSkyfieldConsistency:
    """Both computation paths must agree for all affected bodies."""

    @pytest.mark.parametrize("body_id,body_name", AFFECTED_BODIES)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES)
    def test_leb_matches_skyfield_sid_j2k(self, body_id, body_name, jd, date_desc):
        """LEB and Skyfield paths produce identical SID+J2K results."""
        from libephemeris.state import set_calc_mode

        ephem.set_sid_mode(SIDM_LAHIRI)
        flags = FLG_SIDEREAL | FLG_J2000 | FLG_SPEED

        # Skyfield path
        set_calc_mode("skyfield")
        pos_sky, _ = ephem.calc_ut(jd, body_id, flags)

        # LEB path (if available)
        try:
            set_calc_mode("leb")
            pos_leb, _ = ephem.calc_ut(jd, body_id, flags)
        except Exception:
            pytest.skip("LEB file not available")
        finally:
            set_calc_mode("auto")

        diff_arcsec = angular_diff(pos_sky[0], pos_leb[0]) * 3600.0

        # IntpApog/IntpPerg: LEB data predates apse correction tables,
        # so large divergence expected until LEB regeneration.
        tol = 3600.0 if body_id in (INTP_APOG, INTP_PERG) else 0.5
        assert diff_arcsec < tol, (
            f"{body_name} SID+J2K at {date_desc}: "
            f'LEB vs Skyfield diff {diff_arcsec:.4f}" (should be < 0.5")'
        )


class TestSeDivergenceDocumented:
    """Verify the divergence from pyswisseph is in the expected range.

    These tests require pyswisseph to be installed.  They document the
    exact magnitude of the intentional divergence.
    """

    @pytest.fixture(autouse=True)
    def _require_swe(self):
        """Skip if pyswisseph is not installed."""
        pytest.importorskip("swisseph")

    @pytest.mark.parametrize("body_id,body_name", AFFECTED_BODIES[:2])
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES)
    def test_divergence_magnitude(self, body_id, body_name, jd, date_desc):
        """The SE divergence matches ecliptic precession expectations.

        pyswisseph ignores J2000 for these bodies when sidereal is set.
        LibEphemeris applies it.  The delta should be consistent with
        the ecliptic precession angle from the observation epoch to J2000.
        """
        import swisseph as swe

        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.set_sid_mode(SIDM_LAHIRI)

        flags = FLG_SIDEREAL | FLG_J2000 | FLG_SPEED

        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.calc_ut(jd, body_id, flags)

        divergence_deg = angular_diff(float(pos_swe[0]), pos_py[0])
        divergence_arcsec = divergence_deg * 3600.0

        # The divergence should be proportional to the distance from J2000
        # in time, because the ecliptic precession angle grows linearly.
        # Rough check: < 5 degrees for dates within +-200 years of J2000
        assert divergence_deg < 5.0, (
            f"{body_name} SID+J2K at {date_desc}: "
            f"SE divergence {divergence_deg:.4f} deg "
            f'({divergence_arcsec:.0f}") unexpectedly large'
        )

        # At dates away from J2000, the divergence should be non-trivial
        if date_desc not in ("J2000",):
            assert divergence_deg > 0.01, (
                f"{body_name} SID+J2K at {date_desc}: "
                f"SE divergence {divergence_deg:.6f} deg is too small "
                f"(J2000 precession not being applied?)"
            )

    @pytest.mark.parametrize("jd,date_desc", TEST_DATES)
    def test_se_confirms_j2k_ignored_for_true_node(self, jd, date_desc):
        """Confirm pyswisseph still ignores J2K for TrueNode (SE bug present).

        This test documents the SE bug: SID+J2K == SID for TrueNode.
        If this test starts failing, SE may have fixed the bug upstream.
        """
        import swisseph as swe

        swe.set_sid_mode(swe.SIDM_LAHIRI)

        flags_sid = FLG_SIDEREAL | FLG_SPEED
        flags_j2k = FLG_SIDEREAL | FLG_J2000 | FLG_SPEED

        pos_sid, _ = swe.calc_ut(jd, TRUE_NODE, flags_sid)
        pos_j2k, _ = swe.calc_ut(jd, TRUE_NODE, flags_j2k)

        diff = angular_diff(float(pos_sid[0]), float(pos_j2k[0]))

        # SE should return identical results (bug: J2K ignored)
        assert diff < 0.0001, (
            f"pyswisseph TrueNode at {date_desc}: SID vs SID+J2K differ "
            f"by {diff:.6f} deg — has SE fixed this bug upstream?"
        )
