"""
LEB vs Skyfield Comparison: Sidereal Mode (Extended Tier).

Validates sidereal positions for all body pipelines (A, B, C) across
all flag combinations: SID, SID+EQ, SID+J2K, SID+EQ+J2K.

This fills the critical gap: no sidereal tests existed for the extended tier.
It also serves as regression coverage for all 4 sidereal bug fixes.
"""

from __future__ import annotations

import pytest

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
)

from tests.test_leb.compare.conftest import (
    CompareHelper,
    ECLIPTIC_TOLERANCES,
    SIDEREAL_MODE_SWEEP,
    lon_error_arcsec,
)

from .conftest import TOLS_EXT, NUTATION_FLOOR_ARCSEC


# ============================================================================
# BODY LISTS
# ============================================================================

# Pipeline A (ICRS barycentric)
PIPELINE_A_BODIES = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
]

# Pipeline B (ecliptic direct) — mean bodies
PIPELINE_B_MEAN = [
    (MEAN_NODE, "MeanNode"),
    (MEAN_APOG, "MeanApog"),
]

# Pipeline B (ecliptic direct) — true/osculating bodies
PIPELINE_B_TRUE = [
    (TRUE_NODE, "TrueNode"),
    (OSCU_APOG, "OscuApog"),
]

# Pipeline B (ecliptic direct) — interpolated bodies
PIPELINE_B_INTERP = [
    (INTP_APOG, "IntpApog"),
    (INTP_PERG, "IntpPerg"),
]

ALL_ECLIPTIC_BODIES = PIPELINE_B_MEAN + PIPELINE_B_TRUE + PIPELINE_B_INTERP

# Sidereal tolerance: LEB vs Skyfield should be near-zero since both
# use the same selected ayanamsha definition. The position tolerance is the
# Chebyshev approximation error, NOT ayanamsha model difference.
# Extended tier uses 0.005" instead of 0.001" because nutation polynomial
# degradation beyond +-20 centuries from J2000 adds ~0.003" (known limitation #12).
# This still catches real sidereal bugs (10-36" errors) with >2000x margin.
SID_POSITION_TOL = 0.005  # arcsec (relaxed for nutation degradation at extreme dates)
SID_ECLIPTIC_TOL = TOLS_EXT.ECLIPTIC_ARCSEC  # 0.1" (Meeus polynomial limits)

# Cross-mode sidereal speed tolerance.  The sidereal path uses mean precession
# rate (analytical) while Skyfield uses true ayanamsha rate (finite-diff),
# giving a cross-mode delta of ~5e-5 deg/day.  On top of that, the Chebyshev
# speed approximation error for Moon can reach ~0.001 deg/day at extreme dates.
# Combined: up to ~0.002 deg/day for Moon.  Other bodies are well below 0.001.
SID_SPEED_TOL = 0.002  # deg/day (Moon-dominated, other bodies <0.001)

# A subset of public sidereal modes for combinatorial tests
SID_MODES_SUBSET = [
    0,
    1,
    2,
    3,
    5,
]  # Fagan-Bradley, Lahiri, De Luce, Raman, Krishnamurti


# ============================================================================
# SIDEREAL ECLIPTIC (FLG_SIDEREAL only)
# ============================================================================


class TestExtSiderealEcliptic:
    """Sidereal ecliptic longitude: LEB vs Skyfield."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_A_BODIES)
    @pytest.mark.parametrize("sid_mode", SID_MODES_SUBSET)
    def test_pipeline_a_sidereal(
        self,
        compare: CompareHelper,
        ext_dates_50: list[float],
        body_id: int,
        body_name: str,
        sid_mode: int,
    ):
        """Pipeline A sidereal lon matches Skyfield within tolerance."""
        flags = FLG_SPEED | FLG_SIDEREAL
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_50:
            ephem.set_sid_mode(sid_mode, 2451545.0, 0.0)

            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, flags)

            err = lon_error_arcsec(ref[0], leb[0])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < SID_POSITION_TOL, (
            f'{body_name} SID mode={sid_mode}: max error = {max_err:.4f}" '
            f"at JD {worst_jd:.1f}"
        )

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ALL_ECLIPTIC_BODIES)
    def test_pipeline_b_sidereal(
        self,
        compare: CompareHelper,
        ext_dates_50: list[float],
        body_id: int,
        body_name: str,
    ):
        """Pipeline B sidereal lon matches Skyfield within tolerance."""
        flags = FLG_SPEED | FLG_SIDEREAL
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_50:
            ephem.set_sid_mode(1, 2451545.0, 0.0)  # Lahiri

            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, flags)

            err = lon_error_arcsec(ref[0], leb[0])
            if err > max_err:
                max_err = err
                worst_jd = jd

        per_body = ECLIPTIC_TOLERANCES.get(body_id, {}).get("lon", SID_ECLIPTIC_TOL)
        tol = max(per_body, NUTATION_FLOOR_ARCSEC)
        assert max_err < tol, (
            f'{body_name} SID Lahiri: max error = {max_err:.4f}" at JD {worst_jd:.1f}'
        )


# ============================================================================
# SIDEREAL + EQUATORIAL (Bug 1 & Bug 2 regression)
# ============================================================================


class TestExtSiderealEquatorial:
    """Sidereal equatorial: LEB vs Skyfield."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_A_BODIES)
    def test_pipeline_a_sid_eq(
        self,
        compare: CompareHelper,
        ext_dates_50: list[float],
        body_id: int,
        body_name: str,
    ):
        """Pipeline A SID+EQ RA matches Skyfield (Bug 1 regression)."""
        flags = FLG_SPEED | FLG_SIDEREAL | FLG_EQUATORIAL
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_50:
            ephem.set_sid_mode(1, 2451545.0, 0.0)

            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, flags)

            err = lon_error_arcsec(ref[0], leb[0])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < SID_POSITION_TOL, (
            f'{body_name} SID+EQ: max RA error = {max_err:.4f}" at JD {worst_jd:.1f}'
        )

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_B_MEAN + PIPELINE_B_TRUE)
    def test_pipeline_b_sid_eq(
        self,
        compare: CompareHelper,
        ext_dates_50: list[float],
        body_id: int,
        body_name: str,
    ):
        """Pipeline B SID+EQ RA matches Skyfield (Bug 2 regression)."""
        flags = FLG_SPEED | FLG_SIDEREAL | FLG_EQUATORIAL
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_50:
            ephem.set_sid_mode(1, 2451545.0, 0.0)

            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, flags)

            err = lon_error_arcsec(ref[0], leb[0])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < SID_ECLIPTIC_TOL, (
            f'{body_name} SID+EQ: max RA error = {max_err:.4f}" at JD {worst_jd:.1f}'
        )


# ============================================================================
# SIDEREAL + J2000 (Bug 3 & Bug 4 regression)
# ============================================================================


class TestExtSiderealJ2000:
    """Sidereal J2000: LEB vs Skyfield."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_A_BODIES)
    def test_pipeline_a_sid_j2k(
        self,
        compare: CompareHelper,
        ext_dates_50: list[float],
        body_id: int,
        body_name: str,
    ):
        """Pipeline A SID+J2K lon matches Skyfield."""
        flags = FLG_SPEED | FLG_SIDEREAL | FLG_J2000
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_50:
            ephem.set_sid_mode(1, 2451545.0, 0.0)

            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, flags)

            err = lon_error_arcsec(ref[0], leb[0])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < SID_POSITION_TOL, (
            f'{body_name} SID+J2K: max error = {max_err:.4f}" at JD {worst_jd:.1f}'
        )

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_B_TRUE)
    def test_true_body_sid_j2k_leb_vs_skyfield(
        self,
        compare: CompareHelper,
        ext_dates_50: list[float],
        body_id: int,
        body_name: str,
    ):
        """True bodies SID+J2K: LEB vs Skyfield (SE bug fix applied).

        LibEphemeris intentionally honors FLG_J2000 for TrueNode/OscuApog
        (the reference ephemeris silently ignores it — this is a behavioral bug).
        Both LEB and Skyfield paths must produce consistent results.
        See docs/comparison/intentional-divergences.md
        """
        flags = FLG_SPEED | FLG_SIDEREAL | FLG_J2000
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_50[:10]:
            ephem.set_sid_mode(1, 2451545.0, 0.0)

            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, flags)

            err = lon_error_arcsec(ref[0], leb[0])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < SID_ECLIPTIC_TOL, (
            f'{body_name} SID+J2K: max LEB vs Skyfield error = {max_err:.4f}" '
            f"at JD {worst_jd:.1f}"
        )

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_B_MEAN)
    def test_mean_body_sid_j2k(
        self,
        compare: CompareHelper,
        ext_dates_50: list[float],
        body_id: int,
        body_name: str,
    ):
        """Mean bodies SID+J2K: LEB vs Skyfield (Bug 4 regression)."""
        flags = FLG_SPEED | FLG_SIDEREAL | FLG_J2000
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_50:
            ephem.set_sid_mode(1, 2451545.0, 0.0)

            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, flags)

            err = lon_error_arcsec(ref[0], leb[0])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < SID_ECLIPTIC_TOL, (
            f'{body_name} SID+J2K: max error = {max_err:.4f}" at JD {worst_jd:.1f}'
        )


# ============================================================================
# SIDEREAL + EQUATORIAL + J2000 (all bugs combined)
# ============================================================================


class TestExtSiderealEqJ2k:
    """Sidereal equatorial J2000: LEB vs Skyfield."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_A_BODIES[:3])
    def test_pipeline_a_sid_eq_j2k(
        self,
        compare: CompareHelper,
        ext_dates_50: list[float],
        body_id: int,
        body_name: str,
    ):
        """Pipeline A SID+EQ+J2K RA matches Skyfield."""
        flags = FLG_SPEED | FLG_SIDEREAL | FLG_EQUATORIAL | FLG_J2000
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_50:
            ephem.set_sid_mode(1, 2451545.0, 0.0)

            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, flags)

            err = lon_error_arcsec(ref[0], leb[0])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < SID_POSITION_TOL, (
            f'{body_name} SID+EQ+J2K: max RA error = {max_err:.4f}" '
            f"at JD {worst_jd:.1f}"
        )

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_B_MEAN)
    def test_mean_body_sid_eq_j2k(
        self,
        compare: CompareHelper,
        ext_dates_50: list[float],
        body_id: int,
        body_name: str,
    ):
        """Mean bodies SID+EQ+J2K: LEB vs Skyfield."""
        flags = FLG_SPEED | FLG_SIDEREAL | FLG_EQUATORIAL | FLG_J2000
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_50:
            ephem.set_sid_mode(1, 2451545.0, 0.0)

            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, flags)

            err = lon_error_arcsec(ref[0], leb[0])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < SID_ECLIPTIC_TOL, (
            f'{body_name} SID+EQ+J2K: max RA error = {max_err:.4f}" '
            f"at JD {worst_jd:.1f}"
        )


# ============================================================================
# SIDEREAL SPEED
# ============================================================================


class TestExtSiderealSpeed:
    """Sidereal speed precision across all flag combos."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", PIPELINE_A_BODIES[:3])
    @pytest.mark.parametrize(
        "extra_flags,flag_name",
        [
            (0, "SID"),
            (FLG_EQUATORIAL, "SID+EQ"),
            (FLG_J2000, "SID+J2K"),
            (FLG_EQUATORIAL | FLG_J2000, "SID+EQ+J2K"),
        ],
    )
    def test_sidereal_speed(
        self,
        compare: CompareHelper,
        ext_dates_50: list[float],
        body_id: int,
        body_name: str,
        extra_flags: int,
        flag_name: str,
    ):
        """Sidereal speed matches Skyfield within tolerance."""
        flags = FLG_SPEED | FLG_SIDEREAL | extra_flags
        max_err = 0.0

        for jd in ext_dates_50:
            ephem.set_sid_mode(1, 2451545.0, 0.0)

            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, flags)

            err = abs(ref[3] - leb[3])
            max_err = max(max_err, err)

        assert max_err < SID_SPEED_TOL, (
            f"{body_name} {flag_name}: max speed error = {max_err:.6f} deg/day"
        )


# ============================================================================
# MULTI-AYANAMSHA SWEEP
# ============================================================================


class TestExtSiderealMultiAyanamsha:
    """Broad public ayanamsha-mode sweep for sidereal positions."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", [(SUN, "Sun"), (MARS, "Mars")])
    @pytest.mark.parametrize("sid_mode", SIDEREAL_MODE_SWEEP)
    def test_all_sidereal_modes(
        self,
        compare: CompareHelper,
        ext_dates_50: list[float],
        body_id: int,
        body_name: str,
        sid_mode: int,
    ):
        """Sidereal longitude matches the direct backend for every swept mode."""
        flags = FLG_SPEED | FLG_SIDEREAL
        max_err = 0.0

        for jd in ext_dates_50[:10]:  # Reduced for combinatorial test
            ephem.set_sid_mode(sid_mode, 2451545.0, 0.0)

            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, flags)

            err = lon_error_arcsec(ref[0], leb[0])
            max_err = max(max_err, err)

        assert max_err < SID_POSITION_TOL, (
            f'{body_name} SID mode={sid_mode}: max error = {max_err:.4f}"'
        )
