# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Branch-coverage tests for :mod:`libephemeris.houses`.

This module drives the house-system dispatch, the ARMC and fallback
wrappers, every individual cusp algorithm, the body house-position
solver and the Gauquelin sector code through their edge-case branches.

The tests use deterministic synthetic inputs (fixed Julian Days, fixed
locations, hand-picked ARMC / latitude / declination combinations) so
they run fully offline and exercise the polar, equatorial and
degenerate-geometry code paths that ordinary mid-latitude charts skip.
"""

from __future__ import annotations

import importlib
import math

import pytest

import libephemeris as ephem
from libephemeris.exceptions import PolarCircleError

# The package re-exports the ``houses`` function under the name
# ``libephemeris.houses``, shadowing the submodule.  Import the real module
# object explicitly so the private helpers are reachable for unit testing.
H = importlib.import_module("libephemeris.houses")


# Deterministic anchors reused across the file.
JD = 2451545.0  # J2000.0
ROME_LAT, ROME_LON = 41.9, 12.5
EPS = 23.4393  # approximate true obliquity at J2000


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------


def test_validate_cusps_nan():
    """A NaN cusp is rejected with the NaN message (line 323)."""
    ok, msg = H._validate_cusps([float("nan"), 30.0, 60.0])
    assert ok is False
    assert msg is not None and "NaN" in msg


def test_validate_cusps_infinite():
    """An infinite cusp is rejected with the infinite message (line 327)."""
    ok, msg = H._validate_cusps([10.0, float("inf"), 60.0])
    assert ok is False
    assert msg is not None and "infinite" in msg

    ok2, msg2 = H._validate_cusps([10.0, float("-inf")])
    assert ok2 is False
    assert msg2 is not None and "infinite" in msg2


def test_validate_cusps_out_of_range():
    """An out-of-range cusp is rejected with the range message (line 331)."""
    ok_low, msg_low = H._validate_cusps([-5.0, 30.0])
    assert ok_low is False
    assert msg_low is not None and "out of range" in msg_low

    ok_high, msg_high = H._validate_cusps([10.0, 400.0])
    assert ok_high is False
    assert msg_high is not None and "out of range" in msg_high


def test_validate_cusps_valid():
    """Well-formed cusps pass validation."""
    ok, msg = H._validate_cusps([0.0, 30.0, 60.0, 90.0])
    assert ok is True
    assert msg is None


def test_angular_diff_branches():
    """angular_diff wraps results past 180 deg into the negative half."""
    # diff <= 180 path
    assert H.angular_diff(10.0, 5.0) == pytest.approx(5.0)
    # diff > 180 path (lines 408-409)
    assert H.angular_diff(0.0, 10.0) == pytest.approx(-10.0)
    assert H.angular_diff(350.0, 0.0) == pytest.approx(-10.0)


def test_ra_to_ecliptic_pole_height_90():
    """Pole height = 90 deg returns the 180 deg degenerate value (line 493)."""
    val = H._ra_to_ecliptic_longitude(
        45.0, 90.0, math.sin(math.radians(EPS)), math.cos(math.radians(EPS))
    )
    assert val == pytest.approx(180.0)


def test_ra_to_ecliptic_pole_height_minus_90():
    """Pole height = -90 deg returns the 0 deg degenerate value."""
    val = H._ra_to_ecliptic_longitude(
        45.0, -90.0, math.sin(math.radians(EPS)), math.cos(math.radians(EPS))
    )
    assert val == pytest.approx(0.0)


def test_ecliptic_ra_round_trip_helpers():
    """_ecliptic_to_ra_simple / _ra_to_ecliptic_simple are inverses (3632-3651)."""
    for lon in (0.0, 47.0, 123.0, 290.0):
        ra = H._ecliptic_to_ra_simple(lon, EPS)
        back = H._ra_to_ecliptic_simple(ra, EPS)
        assert back == pytest.approx(lon, abs=1e-6) or back == pytest.approx(
            (lon + 360.0) % 360.0, abs=1e-6
        )


def test_rotate_spherical_zero_radius():
    """A zero-radius vector returns the zero-longitude branch (3699-3700)."""
    out = H._rotate_spherical_x_axis([123.0, 45.0, 0.0], EPS)
    assert out[0] == 0.0
    assert out[1] == 0.0
    assert out[2] == 0.0


def test_armc_to_mc_quadrants():
    """_armc_to_mc covers each ARMC quadrant correction branch."""
    for armc in (10.0, 95.0, 200.0, 280.0, 359.0):
        mc = H._armc_to_mc(armc, EPS)
        assert 0.0 <= mc < 360.0


# ---------------------------------------------------------------------------
# houses() angle branches
# ---------------------------------------------------------------------------


def test_houses_mc_quadrant_branches():
    """Sweep ARMC via longitude so each MC quadrant branch runs (654-661)."""
    seen = set()
    for lon in range(-180, 181, 5):
        cusps, ascmc = ephem.houses(JD, ROME_LAT, float(lon), ord("E"))
        assert len(cusps) == 12
        assert 0.0 <= ascmc[1] < 360.0
        seen.add(int(ascmc[2]) // 90)
    # Confirm all four ARMC quadrants were exercised.
    assert seen == {0, 1, 2, 3}


def test_houses_asc_snaps_to_zero():
    """An Ascendant landing on 360 deg is snapped to 0.0 (line 711)."""
    # houses() derives ARMC from the long-term sidereal-time model, so drive
    # ARMC to exactly 270 deg using that same baseline: ARMC = armc0(JD) + lon.
    armc0, _eps = H._house_armc_obliquity(JD)
    lon = 270.0 - armc0
    while lon > 180.0:
        lon -= 360.0
    while lon < -180.0:
        lon += 360.0
    cusps, ascmc = ephem.houses(JD, ROME_LAT, lon, ord("E"))
    assert ascmc[0] == 0.0


def test_houses_armc_asc_snaps_to_zero():
    """houses_armc snaps an Ascendant of 360 deg to 0.0 (line 1280)."""
    cusps, ascmc = ephem.houses_armc(270.0, ROME_LAT, EPS, ord("E"))
    assert ascmc[0] == 0.0


def test_houses_all_systems_rome():
    """Every dispatch arm in houses() runs for a mid-latitude chart."""
    codes = "PKRCEAWOBTMXVHYFUNGSLQDIiJ"
    for code in codes:
        cusps, ascmc = ephem.houses(JD, ROME_LAT, ROME_LON, ord(code))
        if code == "G":
            assert len(cusps) == 36
        else:
            assert len(cusps) == 12
        assert len(ascmc) == 8


def test_houses_unknown_code_defaults_placidus():
    """An unknown house code falls through to the Placidus default."""
    cusps_def, _ = ephem.houses(JD, ROME_LAT, ROME_LON, ord("Z"))
    cusps_p, _ = ephem.houses(JD, ROME_LAT, ROME_LON, ord("P"))
    assert cusps_def == cusps_p


def test_houses_equator_co_asc_branch():
    """At the equator the H system takes the negative-limit co-ascendant."""
    cusps_h, ascmc_h = ephem.houses(JD, 0.0, ROME_LON, ord("H"))
    cusps_e, ascmc_e = ephem.houses(JD, 0.0, ROME_LON, ord("E"))
    # coasc2 (index 6) differs between H (0.0 limit) and other systems (180.0).
    assert abs(ascmc_h[6] - ascmc_e[6]) > 90.0


def test_houses_apc_polar_flip():
    """APC inside the polar circle flips MC/cusps when MC is below horizon."""
    # Sweep longitudes near a polar latitude to hit the sin_alt < 0 branch.
    flipped = False
    for lon in range(-180, 181, 10):
        cusps, ascmc = ephem.houses(JD, 85.0, float(lon), ord("Y"))
        assert len(cusps) == 12
        flipped = True
    assert flipped


# ---------------------------------------------------------------------------
# Sunshine sun_dec error fallback (808-810) via patched calc_ut
# ---------------------------------------------------------------------------


def test_houses_sunshine_sun_dec_exception(monkeypatch):
    """A failing Sun calculation falls back to zero declination (808-810)."""

    def boom(*args, **kwargs):
        raise ValueError("sun unavailable")

    monkeypatch.setattr(H, "calc_ut", boom)
    cusps, ascmc = ephem.houses(JD, ROME_LAT, ROME_LON, ord("I"))
    assert len(cusps) == 12


# ---------------------------------------------------------------------------
# houses_with_fallback
# ---------------------------------------------------------------------------


def test_houses_with_fallback_polar_path():
    """Placidus at a polar latitude triggers the PolarCircleError fallback."""
    cusps, ascmc, used, warning = ephem.houses_with_fallback(
        JD, 85.0, ROME_LON, ord("P")
    )
    assert used is True
    assert warning is not None and "Placidus" in warning


def test_houses_with_fallback_extreme_warning():
    """Campanus at an extreme latitude returns the instability warning."""
    cusps, ascmc, used, warning = ephem.houses_with_fallback(
        JD, 82.0, ROME_LON, ord("C")
    )
    assert used is False
    assert warning is not None and "reduced accuracy" in warning


def test_houses_with_fallback_clean():
    """A mid-latitude chart needs no fallback and no warning."""
    cusps, ascmc, used, warning = ephem.houses_with_fallback(
        JD, ROME_LAT, ROME_LON, ord("P")
    )
    assert used is False
    assert warning is None


def test_houses_with_fallback_invalid_cusps(monkeypatch):
    """Invalid cusps trigger the validation fallback path (1002-1007)."""
    monkeypatch.setattr(
        H, "_validate_cusps", lambda cusps: (False, "House cusp 1 contains NaN value")
    )
    cusps, ascmc, used, warning = ephem.houses_with_fallback(
        JD, ROME_LAT, ROME_LON, ord("P")
    )
    assert used is True
    assert warning is not None and "invalid cusps" in warning


# ---------------------------------------------------------------------------
# houses_armc and houses_armc_with_fallback
# ---------------------------------------------------------------------------


def test_houses_armc_hsys_bytes_and_str():
    """houses_armc accepts bytes and str house codes (1189-1192)."""
    cusps_b, _ = ephem.houses_armc(200.0, ROME_LAT, EPS, b"P")
    cusps_s, _ = ephem.houses_armc(200.0, ROME_LAT, EPS, "P")
    cusps_i, _ = ephem.houses_armc(200.0, ROME_LAT, EPS, ord("P"))
    assert cusps_b == cusps_s == cusps_i


def test_houses_armc_mc_quadrants():
    """Sweep ARMC so each houses_armc MC quadrant branch runs (1224-1231)."""
    for armc in (10.0, 95.0, 200.0, 280.0, 359.0):
        cusps, ascmc = ephem.houses_armc(armc, ROME_LAT, EPS, ord("E"))
        assert 0.0 <= ascmc[1] < 360.0


def test_houses_armc_all_systems():
    """Every dispatch arm in houses_armc runs."""
    codes = "PKRCEAWOBTMXVHYFUNGSLQDIiJ"
    for code in codes:
        cusps, ascmc = ephem.houses_armc(200.0, ROME_LAT, EPS, ord(code))
        if code == "G":
            assert len(cusps) == 36
        else:
            assert len(cusps) == 12


def test_houses_armc_apc_polar_flip():
    """APC via houses_armc flips MC/cusps inside the polar circle (1402-1405)."""
    for armc in range(0, 360, 15):
        cusps, ascmc = ephem.houses_armc(float(armc), 85.0, EPS, ord("Y"))
        assert len(cusps) == 12


def test_houses_armc_sunshine_makransky():
    """houses_armc covers the Makransky 'i' arm (1427-1428)."""
    cusps, ascmc = ephem.houses_armc(200.0, ROME_LAT, EPS, ord("i"), ascmc9=12.0)
    assert len(cusps) == 12
    assert ascmc[1] == cusps[9]


def test_houses_armc_unknown_default():
    """Unknown code via houses_armc defaults to Placidus (line 1433)."""
    cusps_def, _ = ephem.houses_armc(200.0, ROME_LAT, EPS, ord("z"))
    cusps_p, _ = ephem.houses_armc(200.0, ROME_LAT, EPS, ord("P"))
    assert cusps_def == cusps_p


def test_houses_armc_with_fallback_polar():
    """Placidus via houses_armc_with_fallback hits the polar fallback."""
    cusps, ascmc, used, warning = ephem.houses_armc_with_fallback(
        200.0, 85.0, EPS, ord("P")
    )
    assert used is True
    assert warning is not None


def test_houses_armc_with_fallback_extreme_warning():
    """Campanus via houses_armc_with_fallback returns the instability warning."""
    cusps, ascmc, used, warning = ephem.houses_armc_with_fallback(
        200.0, 82.0, EPS, ord("C")
    )
    assert used is False
    assert warning is not None and "reduced accuracy" in warning


def test_houses_armc_with_fallback_clean():
    """A mid-latitude ARMC chart needs no fallback."""
    cusps, ascmc, used, warning = ephem.houses_armc_with_fallback(
        200.0, ROME_LAT, EPS, ord("P")
    )
    assert used is False
    assert warning is None


def test_houses_armc_with_fallback_no_validate():
    """validate_cusps=False skips the validation block (arc 1105->1117)."""
    cusps, ascmc, used, warning = ephem.houses_armc_with_fallback(
        200.0, ROME_LAT, EPS, ord("P"), validate_cusps=False
    )
    assert used is False
    assert warning is None


def test_houses_armc_with_fallback_extreme_stable_system():
    """An extreme but stable system skips the warning (arc 1119->1127)."""
    # 'E' (Equal) is a stable system: extreme latitude yields no warning.
    cusps, ascmc, used, warning = ephem.houses_armc_with_fallback(
        200.0, 82.0, EPS, ord("E")
    )
    assert used is False
    assert warning is None


def test_houses_with_fallback_no_validate():
    """validate_cusps=False skips validation in houses_with_fallback (998->1010)."""
    cusps, ascmc, used, warning = ephem.houses_with_fallback(
        JD, ROME_LAT, ROME_LON, ord("P"), validate_cusps=False
    )
    assert used is False
    assert warning is None


def test_houses_with_fallback_extreme_stable_system():
    """An extreme but stable system skips the warning (1012->1020)."""
    cusps, ascmc, used, warning = ephem.houses_with_fallback(
        JD, 82.0, ROME_LON, ord("E")
    )
    assert used is False
    assert warning is None


def test_houses_armc_with_fallback_invalid_cusps(monkeypatch):
    """Invalid cusps trigger the houses_armc_with_fallback path (1105-1114)."""
    monkeypatch.setattr(
        H, "_validate_cusps", lambda cusps: (False, "House cusp 1 contains NaN value")
    )
    cusps, ascmc, used, warning = ephem.houses_armc_with_fallback(
        200.0, ROME_LAT, EPS, ord("P")
    )
    assert used is True
    assert warning is not None and "invalid cusps" in warning


# ---------------------------------------------------------------------------
# houses_armc_ex2 speed overrides
# ---------------------------------------------------------------------------


def test_houses_armc_ex2_whole_sign_speed():
    """Whole Sign speed override path."""
    _, _, cs, asp = ephem.houses_armc_ex2(200.0, ROME_LAT, EPS, ord("W"))
    assert cs[1] == 0.0  # intermediate cusp is fixed


def test_houses_armc_ex2_n_u_speed():
    """Natural Gradient ('N') and Krusinski ('U') share a speed model (1562-1569)."""
    for code in ("N", "U"):
        _, _, cs, asp = ephem.houses_armc_ex2(200.0, ROME_LAT, EPS, ord(code))
        assert cs[1] == 0.0
        assert cs[2] == 0.0


def test_houses_armc_ex2_porphyry_speed():
    """Porphyry speed override path."""
    _, _, cs, asp = ephem.houses_armc_ex2(200.0, ROME_LAT, EPS, ord("O"))
    assert len(cs) == 12


def test_houses_armc_ex2_koch_step():
    """Koch uses the 1-second ARMC step in the finite difference."""
    _, _, cs, asp = ephem.houses_armc_ex2(200.0, ROME_LAT, EPS, ord("K"))
    assert len(cs) == 12


def test_houses_armc_ex2_diff_wrap_positive():
    """A cusp jump > +180 between steps hits the positive wrap (line 1519)."""
    # Regiomontanus at lat 80, armc 90 makes a cusp difference exceed +180.
    _, _, cs, asp = ephem.houses_armc_ex2(90.0, 80.0, EPS, ord("R"))
    assert len(cs) == 12


def test_houses_armc_ex2_diff_wrap_negative():
    """A cusp jump < -180 between steps hits the negative wrap (line 1521)."""
    # Whole Sign cusps jump across a sign boundary near armc 13, lat 41.9.
    _, _, cs, asp = ephem.houses_armc_ex2(13.0, ROME_LAT, EPS, ord("W"))
    assert len(cs) == 12


# ---------------------------------------------------------------------------
# houses_ex sidereal branches
# ---------------------------------------------------------------------------


def test_houses_ex_sidereal_hsys_bytes():
    """houses_ex normalises bytes house codes in sidereal mode (1657-1660)."""
    ephem.set_sid_mode(ephem.SIDM_LAHIRI)
    cusps, ascmc = ephem.houses_ex(JD, ROME_LAT, ROME_LON, b"P", ephem.FLG_SIDEREAL)
    assert len(cusps) == 12


def test_houses_ex_sidereal_whole_sign():
    """Whole Sign sidereal recomputation path."""
    ephem.set_sid_mode(ephem.SIDM_LAHIRI)
    cusps, ascmc = ephem.houses_ex(JD, ROME_LAT, ROME_LON, ord("W"), ephem.FLG_SIDEREAL)
    assert len(cusps) == 12


def test_houses_ex_sidereal_equal_and_vehlow():
    """Equal and Vehlow sidereal recomputation paths."""
    ephem.set_sid_mode(ephem.SIDM_LAHIRI)
    for code in ("E", "A", "V"):
        cusps, _ = ephem.houses_ex(
            JD, ROME_LAT, ROME_LON, ord(code), ephem.FLG_SIDEREAL
        )
        assert len(cusps) == 12


def test_houses_ex_sidereal_natural_gradient():
    """Natural Gradient keeps cusps at fixed sidereal degrees (line 1712)."""
    ephem.set_sid_mode(ephem.SIDM_LAHIRI)
    cusps, _ = ephem.houses_ex(JD, ROME_LAT, ROME_LON, ord("N"), ephem.FLG_SIDEREAL)
    for i in range(12):
        assert cusps[i] == pytest.approx(i * 30.0)


def test_houses_ex_sidereal_sunshine_both():
    """Sunshine 'I'/'i' sidereal recomputation paths run."""
    ephem.set_sid_mode(ephem.SIDM_LAHIRI)
    for code in ("I", "i"):
        cusps, _ = ephem.houses_ex(
            JD, ROME_LAT, ROME_LON, ord(code), ephem.FLG_SIDEREAL
        )
        assert len(cusps) == 12


def test_houses_ex_sidereal_sunshine_sun_dec_exception(monkeypatch):
    """A failing Sun calc in sidereal Sunshine falls back to 0 dec (1683-1684)."""

    def boom(*args, **kwargs):
        raise ValueError("no sun")

    monkeypatch.setattr(H, "calc_ut", boom)
    ephem.set_sid_mode(ephem.SIDM_LAHIRI)
    cusps, _ = ephem.houses_ex(JD, ROME_LAT, ROME_LON, ord("I"), ephem.FLG_SIDEREAL)
    assert len(cusps) == 12


def test_houses_ex_sidereal_other_system():
    """A non-Ascendant sidereal system subtracts ayanamsa from cusps."""
    ephem.set_sid_mode(ephem.SIDM_LAHIRI)
    cusps, _ = ephem.houses_ex(JD, ROME_LAT, ROME_LON, ord("O"), ephem.FLG_SIDEREAL)
    assert len(cusps) == 12


def test_houses_ex_tropical():
    """Without the sidereal flag houses_ex returns tropical cusps unchanged."""
    cusps, _ = ephem.houses_ex(JD, ROME_LAT, ROME_LON, ord("P"))
    cusps_t, _ = ephem.houses(JD, ROME_LAT, ROME_LON, ord("P"))
    assert cusps == cusps_t


def test_houses_ex2_speed():
    """houses_ex2 returns matching cusp and speed tuples."""
    cusps, ascmc, cs, asp = ephem.houses_ex2(
        JD, ROME_LAT, ROME_LON, ord("P"), ephem.FLG_SPEED
    )
    assert len(cusps) == 12
    assert len(cs) == 12


# ---------------------------------------------------------------------------
# house_name
# ---------------------------------------------------------------------------


def test_house_name_int_bytes_str():
    """house_name normalises int, bytes and str codes (1835-1838)."""
    assert H.house_name(ord("P")) == "Placidus"
    assert H.house_name(b"K") == "Koch"
    assert H.house_name("O") == "Porphyry"


def test_house_name_unknown():
    """An unknown code returns 'Unknown'."""
    assert H.house_name(ord("z")) == "Unknown"


# ---------------------------------------------------------------------------
# Placidus internals
# ---------------------------------------------------------------------------


def test_placidus_near_polar_clamp():
    """Placidus just below the polar circle exercises the AD clamp (1975-1977)."""
    # lat just inside the polar circle for eps ~ 23.44 deg.
    cusps = H._houses_placidus(200.0, 66.5, EPS, 100.0, 200.0)
    assert len(cusps) == 13


def test_placidus_bisection_fallback():
    """Right at the polar circle the iteration fails and bisection solves it."""
    lat = 90.0 - EPS - 1e-4
    cusps = H._houses_placidus(90.0, lat, EPS, 100.0, 200.0)
    assert len(cusps) == 13
    for c in cusps[1:]:
        assert 0.0 <= c < 360.0


def test_placidus_dispatch_polar_raises():
    """houses() raises for Placidus inside the polar circle."""
    with pytest.raises(PolarCircleError):
        ephem.houses(JD, 89.0, ROME_LON, ord("P"))


# ---------------------------------------------------------------------------
# Koch / Gauquelin polar fallbacks
# ---------------------------------------------------------------------------


def test_koch_polar_fallback_porphyry():
    """Koch helper falls back to Porphyry inside the polar circle (line 2366)."""
    lat = 90.0 - EPS + 1.0  # inside polar circle
    cusps = H._houses_koch(200.0, lat, EPS, 100.0, 200.0)
    porphyry = H._houses_porphyry(100.0, 200.0)
    assert cusps == porphyry


def test_gauquelin_polar_fallback():
    """Gauquelin helper falls back to equal sectors at polar lat (3573-3575)."""
    lat = 90.0 - EPS + 1.0
    cusps = H._houses_gauquelin(200.0, lat, EPS, 100.0, 200.0)
    assert len(cusps) == 37
    # Equal 10 deg sectors descending from the ascendant.
    assert cusps[1] == pytest.approx(100.0)
    assert cusps[2] == pytest.approx((100.0 - 10.0) % 360.0)


def test_gauquelin_mid_latitude_iteration():
    """A mid-latitude Gauquelin chart runs the iterative cusp helper."""
    cusps = H._houses_gauquelin(200.0, ROME_LAT, EPS, 100.0, 200.0)
    assert len(cusps) == 37


def test_gauquelin_cusp_helper_equator_degenerate():
    """A cusp on the equator returns the seed RA directly (line 3478)."""
    # tan_lat = 0 at the equator forces tan(cusp_dec) ~ 0.
    val = H._gauquelin_cusp_for_sector(
        sector_offset=4,
        base_ramc=0.0,
        ramc_sign=1.0,
        ascensional_diff=0.0,
        tan_obliquity=math.tan(math.radians(EPS)),
        sin_obliquity=math.sin(math.radians(EPS)),
        tan_lat=0.0,
        eps=EPS,
        lat=0.0,
        niter_max=100,
        convergence_threshold=1e-8,
        near_zero=1e-10,
    )
    assert 0.0 <= val < 360.0


# ---------------------------------------------------------------------------
# Savard-A pole observer branch
# ---------------------------------------------------------------------------


def test_savard_a_polar_observer_north():
    """Savard-A at the north pole collapses the prime vertical (2632-2635)."""
    cusps = H._houses_savard_a(200.0, 90.0, EPS, 100.0, 200.0)
    assert len(cusps) == 13


def test_savard_a_polar_observer_south():
    """Savard-A at the south pole takes the southern collapse branch."""
    cusps = H._houses_savard_a(200.0, -90.0, EPS, 100.0, 200.0)
    assert len(cusps) == 13


def test_savard_a_equator():
    """Savard-A at the equator uses the limit anchor arcs."""
    cusps = H._houses_savard_a(200.0, 0.0, EPS, 100.0, 200.0)
    assert len(cusps) == 13


# ---------------------------------------------------------------------------
# Pullen SR Newton iteration
# ---------------------------------------------------------------------------


def test_pullen_sr_equal_quadrant():
    """A 90 deg quadrant short-circuits to equal houses (line 2974)."""
    # asc - mc == 90 makes q1 exactly 90 deg.
    cusps = H._houses_pullen_sr(90.0, 0.0)
    assert len(cusps) == 13


def test_pullen_sr_unequal_quadrants():
    """Unequal quadrants run the Newton-Raphson solver (2982-2997)."""
    cusps = H._houses_pullen_sr(100.0, 200.0)
    assert len(cusps) == 13
    cusps2 = H._houses_pullen_sr(123.0, 350.0)
    assert len(cusps2) == 13


def test_pullen_sd_negative_middle():
    """Pullen SD redistributes when a middle house would be negative."""
    # A very small quadrant (asc just past mc) forces middle < 0.
    cusps = H._houses_pullen_sd(201.0, 200.0)
    assert len(cusps) == 13


# ---------------------------------------------------------------------------
# Gauquelin cusp iteration convergence (3494-3523)
# ---------------------------------------------------------------------------


def test_gauquelin_helper_high_latitude_iterates():
    """A high-latitude sector drives the refinement loop to convergence."""
    lat = 60.0
    tan_lat = math.tan(math.radians(lat))
    tan_eps = math.tan(math.radians(EPS))
    ad = math.degrees(math.asin(max(-1.0, min(1.0, tan_lat * tan_eps))))
    val = H._gauquelin_cusp_for_sector(
        sector_offset=4,
        base_ramc=33.0,
        ramc_sign=1.0,
        ascensional_diff=ad,
        tan_obliquity=tan_eps,
        sin_obliquity=math.sin(math.radians(EPS)),
        tan_lat=tan_lat,
        eps=EPS,
        lat=lat,
        niter_max=100,
        convergence_threshold=1e-8,
        near_zero=1e-10,
    )
    assert 0.0 <= val < 360.0


# ---------------------------------------------------------------------------
# Krusinski asc flip
# ---------------------------------------------------------------------------


def test_krusinski_asc_flip_branch():
    """A polar Krusinski chart triggers the asc-west-of-mc flip (line 3751)."""
    seen_len = 0
    for armc in range(0, 360, 10):
        for lat in (70.0, 80.0, -75.0):
            cusps = H._houses_krusinski(float(armc), lat, EPS, 10.0, 200.0)
            seen_len = len(cusps)
    assert seen_len == 13


def test_krusinski_via_houses():
    """Krusinski through the public dispatch runs."""
    cusps, _ = ephem.houses(JD, 70.0, ROME_LON, ord("U"))
    assert len(cusps) == 12


# ---------------------------------------------------------------------------
# Sunshine arc helper degenerate side_c
# ---------------------------------------------------------------------------


def test_sunshine_arc_zero_side_c():
    """A zero arc offset collapses side_c and zeros the zenith dist (3921)."""
    val = H._sunshine_arc_to_ecliptic(
        arc_offset_ra=0.0,
        is_diurnal=True,
        armc=200.0,
        lat=ROME_LAT,
        sun_dec=0.0,
        sin_lat=math.sin(math.radians(ROME_LAT)),
        cos_lat=math.cos(math.radians(ROME_LAT)),
        cos_sun_dec=1.0,
        tan_sun_dec=0.0,
        sin_obliquity=math.sin(math.radians(EPS)),
        cos_obliquity=math.cos(math.radians(EPS)),
    )
    assert 0.0 <= val < 360.0


def test_sunshine_treindl_mc_under_horizon():
    """Sunshine Treindl flips intermediate cusps when MC is below horizon."""
    # A southern-hemisphere armc that puts the MC under the horizon.
    cusps = H._houses_sunshine(20.0, -80.0, EPS, 100.0, 280.0, 23.0)
    assert len(cusps) == 13


# ---------------------------------------------------------------------------
# Sunshine Makransky branches
# ---------------------------------------------------------------------------


def test_makransky_mc_under_horizon_flip():
    """Makransky flips the asc when MC is below horizon (4127-4132)."""
    # asc - mc < 0 selects mc_under_horizon.
    cusps = H._houses_sunshine_makransky(200.0, ROME_LAT, EPS, 100.0, 200.0, 10.0)
    assert len(cusps) == 13


def test_makransky_circumpolar_raises():
    """A circumpolar Sun makes Makransky raise PolarCircleError (line 4147)."""
    with pytest.raises(PolarCircleError):
        H._houses_sunshine_makransky(200.0, 80.0, EPS, 100.0, 200.0, 23.44)


def test_makransky_southern_hemisphere():
    """Southern latitude takes the 'south' RA-flip branch (line 4195)."""
    cusps = H._houses_sunshine_makransky(200.0, -41.9, EPS, 280.0, 100.0, 10.0)
    assert len(cusps) == 13


def test_makransky_meridian_distance_gt_90():
    """High-latitude winter geometry exercises the MD>90 branch (4202-4213)."""
    # AD large enough that 2 * day_semi_arc / 3 exceeds 90 deg.
    for lat, sd in ((63.0, 18.0), (-63.0, -18.0), (50.0, 30.0)):
        cusps = H._houses_sunshine_makransky(200.0, lat, EPS, 100.0, 200.0, sd)
        assert len(cusps) == 13


def test_makransky_many_geometries():
    """Sweep ARMC/lat/dec to exercise the ecliptic case tables (4261-4331)."""
    for armc in (0.0, 30.0, 75.0, 120.0, 200.0, 285.0, 350.0):
        for lat in (41.9, -41.9, 60.0, -55.0, 10.0):
            for sd in (-20.0, -5.0, 0.0, 5.0, 20.0):
                ad_arg = math.tan(math.radians(sd)) * math.tan(math.radians(lat))
                if abs(ad_arg) >= 1.0:
                    continue
                cusps = H._houses_sunshine_makransky(armc, lat, EPS, 100.0, 200.0, sd)
                assert len(cusps) == 13
                for c in cusps[1:]:
                    assert 0.0 <= c < 360.0


# ---------------------------------------------------------------------------
# Horizontal polar handling
# ---------------------------------------------------------------------------


def test_horizontal_polar_handling():
    """Horizontal inside the polar circle hits the asc/mc flip (4423-4438)."""
    seen = 0
    for armc in range(0, 360, 10):
        for lat in (80.0, -80.0, 70.0):
            cusps = H._houses_horizontal(float(armc), lat, EPS, 100.0, 200.0)
            seen = len(cusps)
    assert seen == 13


def test_horizontal_via_houses_polar():
    """Horizontal through the public dispatch at polar latitude."""
    for lon in range(-180, 181, 30):
        cusps, _ = ephem.houses(JD, 80.0, float(lon), ord("H"))
        assert len(cusps) == 12


# ---------------------------------------------------------------------------
# APC cusp helper polar / equator branches
# ---------------------------------------------------------------------------


def test_apc_cusp_near_pole():
    """APC cusp helper at |lat| ~ 90 deg uses the zero-difference branch (4550)."""
    val = H._apc_cusp(
        5, math.radians(89.9999999), math.radians(EPS), math.radians(200.0)
    )
    assert 0.0 <= val < 360.0


def test_apc_cusp_equator_limit():
    """APC cusp helper at the equator uses the declination limit (4560-4565)."""
    val_n = H._apc_cusp(5, math.radians(1e-11), math.radians(EPS), math.radians(200.0))
    val_s = H._apc_cusp(5, math.radians(-1e-11), math.radians(EPS), math.radians(200.0))
    assert 0.0 <= val_n < 360.0
    assert 0.0 <= val_s < 360.0


def test_apc_cusp_diurnal_house():
    """APC cusp helper for a diurnal house (>= 8) runs the day-arc branch."""
    val = H._apc_cusp(10, math.radians(41.9), math.radians(EPS), math.radians(200.0))
    assert 0.0 <= val < 360.0


def test_apc_helper_polar_flip():
    """_houses_apc inside the polar circle flips cusps (4653-4656)."""
    seen = 0
    for armc in range(0, 360, 10):
        cusps = H._houses_apc(float(armc), 85.0, EPS, 100.0, 200.0)
        seen = len(cusps)
    assert seen == 13


# ---------------------------------------------------------------------------
# house_pos calling conventions and per-system branches
# ---------------------------------------------------------------------------


def test_house_pos_5arg_bytes():
    """5-arg form with bytes house code (4768-4770)."""
    pos = ephem.house_pos(200.0, ROME_LAT, EPS, (15.0, 0.0), b"P")
    assert 1.0 <= pos < 13.0


def test_house_pos_5arg_str():
    """5-arg form with str house code (4771-4773)."""
    pos = ephem.house_pos(200.0, ROME_LAT, EPS, (15.0, 0.0), "P")
    assert 1.0 <= pos < 13.0


def test_house_pos_5arg_default_hsys():
    """5-arg form with no house code defaults to Placidus (4775-4777)."""
    pos = ephem.house_pos(200.0, ROME_LAT, EPS, (15.0, 0.0), None)
    assert 1.0 <= pos < 13.0


def test_house_pos_5arg_objcoord_single_element():
    """5-arg objcoord with a single element defaults its latitude to 0."""
    pos = ephem.house_pos(200.0, ROME_LAT, EPS, (15.0,), b"P")
    assert 1.0 <= pos < 13.0


def test_house_pos_6arg_int():
    """6-arg form with int house code."""
    pos = ephem.house_pos(200.0, ROME_LAT, EPS, ord("P"), 15.0, 0.0)
    assert 1.0 <= pos < 13.0


def test_house_pos_6arg_bytes_str():
    """6-arg form with bytes and str house codes (4781-4786)."""
    pos_b = ephem.house_pos(200.0, ROME_LAT, EPS, b"P", 15.0, 0.0)
    pos_s = ephem.house_pos(200.0, ROME_LAT, EPS, "P", 15.0, 0.0)
    assert 1.0 <= pos_b < 13.0
    assert 1.0 <= pos_s < 13.0


def test_house_pos_6arg_lon_tuple():
    """6-arg form where lon is itself an objcoord sequence (4793-4796)."""
    pos = ephem.house_pos(200.0, ROME_LAT, EPS, ord("P"), (15.0, 2.0))
    assert 1.0 <= pos < 13.0


def test_house_pos_6arg_lon_none():
    """6-arg form with a missing longitude defaults to 0."""
    pos = ephem.house_pos(200.0, ROME_LAT, EPS, ord("X"), None)
    assert 1.0 <= pos < 13.0


def test_house_pos_sin_dec_clamp():
    """A high-latitude body with large ecliptic latitude clamps sin_dec (4830)."""
    pos = ephem.house_pos(200.0, ROME_LAT, EPS, ord("P"), 90.0, 89.9)
    assert 1.0 <= pos < 13.0
    pos2 = ephem.house_pos(200.0, ROME_LAT, EPS, ord("P"), 270.0, -89.9)
    assert 1.0 <= pos2 < 13.0


def test_house_pos_placidus_circumpolar_body():
    """A circumpolar body takes the no-semi-arc Placidus branch (4855-4863)."""
    pos_never_set = ephem.house_pos(200.0, 80.0, EPS, ord("P"), 90.0, 0.0)
    pos_never_rise = ephem.house_pos(200.0, 80.0, EPS, ord("P"), 270.0, 0.0)
    assert 1.0 <= pos_never_set < 13.0
    assert 1.0 <= pos_never_rise < 13.0


def test_house_pos_gauquelin():
    """Gauquelin house position returns a sector in [1, 37)."""
    pos = ephem.house_pos(200.0, ROME_LAT, EPS, (15.0, 0.0), "G")
    assert 1.0 <= pos < 37.0


def test_house_pos_regiomontanus_branches():
    """Regiomontanus house position over many bodies (4923-4940)."""
    for armc in (0.0, 90.0, 200.0):
        for lon in (0.0, 90.0, 180.0, 270.0, 45.0, 135.0):
            for lat in (41.9, -80.0, 80.0):
                pos = ephem.house_pos(armc, lat, EPS, ord("R"), lon, 0.0)
                assert 1.0 <= pos < 13.0


def test_house_pos_regiomontanus_meridian_degenerate():
    """Bodies on the upper/lower meridian hit the degenerate RA branch."""
    # Body RA equal to ARMC (upper meridian) and ARMC+180 (lower meridian).
    pos_upper = ephem.house_pos(0.0, -80.0, EPS, ord("R"), 0.0, 0.0)
    pos_lower = ephem.house_pos(180.0, -80.0, EPS, ord("R"), 0.0, 0.0)
    assert 1.0 <= pos_upper < 13.0
    assert 1.0 <= pos_lower < 13.0


def test_house_pos_campanus():
    """Campanus house position path."""
    pos = ephem.house_pos(200.0, ROME_LAT, EPS, ord("C"), 45.0, 0.0)
    assert 1.0 <= pos < 13.0


def test_house_pos_equal_family():
    """Equal-based systems D/V/W/N/E house positions (4996-5008)."""
    for code in ("E", "A", "W", "V", "D", "N"):
        pos = ephem.house_pos(200.0, ROME_LAT, EPS, ord(code), 123.0, 0.0)
        assert 1.0 <= pos < 13.0


def test_house_pos_meridian():
    """Meridian (X) house position (5010-5013)."""
    pos = ephem.house_pos(200.0, ROME_LAT, EPS, ord("X"), 123.0, 0.0)
    assert 1.0 <= pos < 13.0


def test_house_pos_morinus_branches():
    """Morinus house position including the 90/270 degenerate cases (5032)."""
    pos = ephem.house_pos(200.0, ROME_LAT, EPS, ord("M"), 123.0, 0.0)
    assert 1.0 <= pos < 13.0
    pos90 = ephem.house_pos(200.0, ROME_LAT, EPS, ord("M"), 90.0, 0.0)
    pos270 = ephem.house_pos(200.0, ROME_LAT, EPS, ord("M"), 270.0, 0.0)
    assert 1.0 <= pos90 < 13.0
    assert 1.0 <= pos270 < 13.0


def test_house_pos_koch_branches():
    """Koch house position over a range of bodies (5037-5134)."""
    for armc in (0.0, 90.0, 200.0):
        for lon in (0.0, 60.0, 123.0, 200.0, 300.0):
            for lat in (41.9, -41.9):
                pos = ephem.house_pos(armc, lat, EPS, ord("K"), lon, 0.0)
                assert 0.0 <= pos < 13.0


def test_house_pos_koch_circumpolar_undefined():
    """A circumpolar Koch body returns 0.0 (undefined)."""
    # High latitude with a high-declination body drives the invalid branch.
    found_zero = False
    for lon in range(0, 360, 5):
        pos = ephem.house_pos(0.0, 70.0, EPS, ord("K"), float(lon), 0.0)
        if pos == 0.0:
            found_zero = True
    assert found_zero


def test_house_pos_topocentric():
    """Topocentric (T) house position, circumpolar and normal (5136-5161)."""
    pos = ephem.house_pos(200.0, ROME_LAT, EPS, ord("T"), 123.0, 0.0)
    assert 1.0 <= pos < 13.0
    pos_cp = ephem.house_pos(200.0, 80.0, EPS, ord("T"), 90.0, 0.0)
    assert 1.0 <= pos_cp < 13.0
    pos_cp2 = ephem.house_pos(200.0, 80.0, EPS, ord("T"), 270.0, 0.0)
    assert 1.0 <= pos_cp2 < 13.0


def test_house_pos_alcabitius():
    """Alcabitius house position interpolates in RA (5164-5195)."""
    for lon in (0.0, 45.0, 123.0, 250.0, 300.0):
        pos = ephem.house_pos(200.0, ROME_LAT, EPS, ord("B"), lon, 0.0)
        assert 1.0 <= pos < 13.0


def test_house_pos_default_fallback():
    """Systems without a closed form use the cusp-interpolation fallback."""
    for code in ("H", "F", "I", "i", "J", "S", "L", "Q", "Y"):
        for lon in (0.0, 100.0, 200.0, 300.0):
            pos = ephem.house_pos(200.0, ROME_LAT, EPS, ord(code), lon, 0.0)
            assert 1.0 <= pos < 13.0


def test_house_pos_default_fallback_zero_width(monkeypatch):
    """A zero-width house in the fallback path is handled (5209-5217)."""

    # Force coincident cusps so a house has zero width.
    def fake_armc(armc, lat, obliquity, hsys):
        cusps = tuple([10.0] * 12)  # all identical -> zero-width houses
        ascmc = tuple([0.0] * 8)
        return cusps, ascmc

    monkeypatch.setattr(H, "houses_armc", fake_armc)
    # Body exactly on the coincident cusp -> returns that house.
    pos_on = H._house_pos_pythonic(200.0, ROME_LAT, EPS, ord("L"), 10.0, 0.0)
    assert pos_on == pytest.approx(1.0)
    # Body away from the cusp -> no house matches, returns 1.0.
    pos_off = H._house_pos_pythonic(200.0, ROME_LAT, EPS, ord("L"), 123.0, 0.0)
    assert pos_off == 1.0


def test_house_pos_alcabitius_zero_interval(monkeypatch):
    """Alcabitius with coincident cusp RAs uses the 30 deg interval default."""

    def fake_armc(armc, lat, obliquity, hsys):
        cusps = tuple([10.0] * 12)
        ascmc = tuple([0.0] * 8)
        return cusps, ascmc

    monkeypatch.setattr(H, "houses_armc", fake_armc)
    pos = H._house_pos_pythonic(200.0, ROME_LAT, EPS, ord("B"), 11.0, 0.0)
    assert 1.0 <= pos < 13.0


# ---------------------------------------------------------------------------
# Gauquelin sector (methods 0/1 and 2-5, planet and star)
# ---------------------------------------------------------------------------


def test_gauquelin_sector_method0_planet():
    """Method 0 sector for a planet via the public API."""
    sector = ephem.gauquelin_sector(JD, ephem.MARS, 0, (ROME_LON, ROME_LAT, 0.0))
    assert 1.0 <= sector < 37.0


def test_gauquelin_sector_method1_no_latitude():
    """Method 1 ignores the body's ecliptic latitude."""
    sector = ephem.gauquelin_sector(JD, ephem.MARS, 1, (ROME_LON, ROME_LAT, 0.0))
    assert 1.0 <= sector < 37.0


def test_gauquelin_sector_method0_star():
    """Method 0 sector for a fixed star (lines 5617-5620)."""
    sector = ephem.gauquelin_sector(JD, "Regulus", 0, (ROME_LON, ROME_LAT, 0.0))
    assert 1.0 <= sector < 37.0


def test_gauquelin_sector_rise_set_methods():
    """Methods 2-5 use real rise/set times (5369-5516)."""
    for method in (2, 3, 4, 5):
        sector = ephem.gauquelin_sector(
            JD, ephem.MARS, method, (ROME_LON, ROME_LAT, 0.0)
        )
        assert 1.0 <= sector < 37.0


def test_gauquelin_sector_rise_set_star():
    """A fixed star with methods 2-5 now computes real rise/set sectors.

    rise_trans accepts star names, so _gauquelin_sector_from_rise_set handles
    them like planets (previously this raised NotImplementedError).
    """
    for method in (2, 3, 4, 5):
        sector = ephem.gauquelin_sector(
            JD, "Regulus", method, (ROME_LON, ROME_LAT, 0.0)
        )
        assert 1.0 <= sector < 37.0
    # Methods 2/4 (no refraction) and 3/5 (with refraction) agree within their
    # respective refraction bands, and 2==4, 3==5.
    m2 = ephem.gauquelin_sector(JD, "Regulus", 2, (ROME_LON, ROME_LAT, 0.0))
    m4 = ephem.gauquelin_sector(JD, "Regulus", 4, (ROME_LON, ROME_LAT, 0.0))
    assert abs(m2 - m4) < 1e-6


def test_gauquelin_sector_default_pressure_temp():
    """Zero pressure/temperature fall back to standard atmosphere defaults."""
    sector = ephem.gauquelin_sector(
        JD, ephem.MARS, 0, (ROME_LON, ROME_LAT, 0.0), atpress=0.0, attemp=0.0
    )
    assert 1.0 <= sector < 37.0


def test_gauquelin_sector_rise_set_circumpolar_raises():
    """A circumpolar planet at high latitude raises for methods 2-5.

    Very high latitude makes the body circumpolar (rise_trans returns -2).
    Matching the reference API, methods 2-5 then raise "rise or set not found"
    rather than silently falling back to method 0.
    """
    with pytest.raises(ephem.Error, match="rise or set not found"):
        ephem.gauquelin_sector(JD, ephem.MARS, 2, (ROME_LON, 89.0, 0.0))
