"""Regression tests for the Round-2 "core" speed-derivative fixes.

Every fix here enforces the project's speed philosophy: the reported speed is
the true derivative of the reported position, and the LEB and Skyfield backends
agree on it. Covered:

* L1 True Node distance rate (``ddist``) is the derivative of the reconstructed
  osculating-node radius, not the retired proxy channel.
* L2/L3/L4/L5 node/apse LEB speeds match the Skyfield backend (J2000 equinox
  drift, nutation-in-longitude rate, TrueNode/OscuApog smoothing window,
  MeanApog latitude rate).
* H1 hypothetical of-date speeds pick up the general-precession rate.
* H3 fictitious bodies (49-58) carry the nutation-in-longitude rate.
* H4 topocentric parallax is applied to the hypothetical bodies (and raises
  without an observer, like the planets).
* D1 the TT and UT fast-path entry points are coherent.
* M1 sidereal ``nod_aps`` returns native Python floats.
"""

from __future__ import annotations

import math
import os

import pytest

import libephemeris as le
from libephemeris import constants as C
from libephemeris import fast_calc as fc
from libephemeris.fast_calc import fast_calc_tt, fast_calc_ut
from libephemeris.leb_reader import open_leb

LEB_BASE_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
    "data",
    "leb",
    "ephemeris_base.leb",
)

SKIP_NO_LEB = pytest.mark.skipif(
    not os.path.exists(LEB_BASE_PATH), reason="LEB base file not found"
)

# Dates inside the base tier (1860-2140).
JD_2000 = 2451545.0
JD_1987 = 2446896.3
JD_1900 = 2415020.0
JD_2008 = 2454545.0


def _wrap180(d: float) -> float:
    while d > 180.0:
        d -= 360.0
    while d < -180.0:
        d += 360.0
    return d


@pytest.fixture
def reader():
    with open_leb(LEB_BASE_PATH) as r:
        fc._set_active_reader(r)
        yield r


@pytest.fixture
def skyfield_backend():
    le.set_calc_mode("skyfield")
    yield
    le.set_calc_mode(None)
    le.close()


def _leb(reader, jd, ipl, iflag):
    (pos, _rf) = fast_calc_tt(reader, jd, ipl, iflag)
    return pos


def _sky(jd, ipl, iflag):
    pos, _rf = le.calc(jd, ipl, iflag)
    return pos


def _deriv(fn, jd, idx, h):
    lo = fn(jd - h)
    hi = fn(jd + h)
    d = hi[idx] - lo[idx]
    if idx == 0:
        d = _wrap180(d)
    return d / (2.0 * h)


# ---------------------------------------------------------------------------
# L1 — True Node distance rate is the derivative of the reported distance.
# ---------------------------------------------------------------------------


@SKIP_NO_LEB
def test_true_node_ddist_is_derivative_of_reported_distance(reader):
    F = C.FLG_SPEED
    ddist_rep = _leb(reader, JD_1987, C.TRUE_NODE, F)[5]
    ddist_num = _deriv(lambda j: _leb(reader, j, C.TRUE_NODE, F), JD_1987, 2, 0.5)
    # Was ~1.7e-5 reported vs ~3e-6 stale; now the reported rate matches.
    assert abs(ddist_rep - ddist_num) < 5e-9
    # The distance itself is the osculating node radius (~0.0024-0.0026 AU),
    # not the retired ~0.0015 AU proxy.
    assert 0.0022 < _leb(reader, JD_1987, C.TRUE_NODE, F)[2] < 0.0028


# ---------------------------------------------------------------------------
# L2/L3/L4/L5 — node/apse LEB speeds match the Skyfield backend.
# ---------------------------------------------------------------------------


@SKIP_NO_LEB
@pytest.mark.parametrize("ipl", [C.MEAN_NODE, C.TRUE_NODE, C.MEAN_APOG, C.OSCU_APOG])
@pytest.mark.parametrize("extra", [0, C.FLG_J2000, C.FLG_J2000 | C.FLG_EQUATORIAL])
def test_node_apse_speed_leb_matches_skyfield(reader, skyfield_backend, ipl, extra):
    flags = C.FLG_SPEED | extra
    jd = JD_2000
    a = _leb(reader, jd, ipl, flags)
    s = _sky(jd, ipl, flags)
    # OscuApog is the fast-oscillating osculating apse; its equatorial output
    # carries a tiny residual between the LEB finite-difference rotation and
    # the Skyfield analytic cotrans (both self-consistent to <0.01"/day). Every
    # other case agrees to <0.01"/day.
    tol = 0.05 if (ipl == C.OSCU_APOG and (extra & C.FLG_EQUATORIAL)) else 0.01
    # Longitude/RA speed parity (the smoothed TrueNode/OscuApog windows are
    # matched per body, so the two backends agree).
    assert abs(_wrap180(a[3] - s[3]) * 3600.0) < tol
    # Latitude/dec speed parity.
    assert abs((a[4] - s[4]) * 3600.0) < tol
    # Distance rate parity.
    assert abs(a[5] - s[5]) < 5e-9


@SKIP_NO_LEB
def test_mean_node_j2000_speed_differs_from_ofdate_by_precession(reader):
    """L2: the J2000 node speed carries the equinox drift (~0.14"/day)."""
    of_date = _leb(reader, JD_2000, C.MEAN_NODE, C.FLG_SPEED)[3]
    j2000 = _leb(reader, JD_2000, C.MEAN_NODE, C.FLG_SPEED | C.FLG_J2000)[3]
    diff = abs(of_date - j2000) * 3600.0
    assert 0.10 < diff < 0.18


@SKIP_NO_LEB
def test_mean_node_default_speed_carries_nutation_rate(reader, skyfield_backend):
    """L3: the of-date node speed includes d(Δψ)/dt (matched to Skyfield)."""
    jd = JD_1900
    a = _leb(reader, jd, C.MEAN_NODE, C.FLG_SPEED)[3]
    s = _sky(jd, C.MEAN_NODE, C.FLG_SPEED)[3]
    assert abs(_wrap180(a - s) * 3600.0) < 0.01
    # NONUT strips the nutation rate on both backends -> still agree.
    an = _leb(reader, jd, C.MEAN_NODE, C.FLG_SPEED | C.FLG_NONUT)[3]
    sn = _sky(jd, C.MEAN_NODE, C.FLG_SPEED | C.FLG_NONUT)[3]
    assert abs(_wrap180(an - sn) * 3600.0) < 0.01


@SKIP_NO_LEB
def test_mean_apog_latitude_rate_matches_skyfield(reader, skyfield_backend):
    """L5: MeanApog dlat is the derivative of the 3-harmonic latitude model."""
    for jd in (JD_2008, 2453045.0):
        a = _leb(reader, jd, C.MEAN_APOG, C.FLG_SPEED)
        s = _sky(jd, C.MEAN_APOG, C.FLG_SPEED)
        assert abs((a[4] - s[4]) * 3600.0) < 0.01


# ---------------------------------------------------------------------------
# H1 — hypothetical of-date speeds pick up the general-precession rate.
# ---------------------------------------------------------------------------


@SKIP_NO_LEB
@pytest.mark.parametrize(
    "flags",
    [
        C.FLG_SPEED,
        C.FLG_SPEED | C.FLG_HELCTR,
        C.FLG_SPEED | C.FLG_EQUATORIAL,
        C.FLG_SPEED | C.FLG_J2000,
    ],
)
def test_uranian_speed_leb_matches_skyfield(reader, skyfield_backend, flags):
    jd = JD_1900
    a = _leb(reader, jd, C.CUPIDO, flags)
    s = _sky(jd, C.CUPIDO, flags)
    assert abs(_wrap180(a[3] - s[3]) * 3600.0) < 0.01


@SKIP_NO_LEB
def test_uranian_ofdate_speed_is_derivative_of_position(reader):
    F = C.FLG_SPEED
    rep = _leb(reader, JD_1900, C.CUPIDO, F)[3]
    num = _deriv(lambda j: _leb(reader, j, C.CUPIDO, F), JD_1900, 0, 0.5)
    assert abs(_wrap180(rep - num) * 3600.0) < 0.05


# ---------------------------------------------------------------------------
# H1b — Uranian / Transpluto (40-48) carry nutation in longitude (Δψ).
#
# Regression against the former class-split: bodies 40-48 were reported on the
# MEAN ecliptic of date (Δψ omitted, so lon(default) == lon(NONUT)), while every
# other of-date body — and the reference API, which applies Δψ uniformly — puts
# them on the TRUE ecliptic. Both backends must now add Δψ to the position and
# d(Δψ)/dt to the speed.
# ---------------------------------------------------------------------------


def _dpsi_deg(jd_tt: float) -> float:
    """Nutation in longitude Δψ at ``jd_tt``, in degrees."""
    from libephemeris.cache import get_cached_nutation

    dpsi_rad, _ = get_cached_nutation(jd_tt)
    return math.degrees(dpsi_rad)


@SKIP_NO_LEB
@pytest.mark.parametrize("ipl", [C.CUPIDO, C.ISIS])
def test_uranian_default_minus_nonut_is_nutation(reader, skyfield_backend, ipl):
    """lon(default) - lon(NONUT) == Δψ for CUPIDO and ISIS on both backends."""
    for jd in (JD_1900, JD_2000, JD_2008):
        expected = _dpsi_deg(jd)  # degrees; ~ ±0.005 deg = ±17"
        for name, get in (
            ("leb", lambda j, f: _leb(reader, j, ipl, f)),
            ("sky", lambda j, f: _sky(j, ipl, f)),
        ):
            lon_def = get(jd, C.FLG_SPEED)[0]
            lon_non = get(jd, C.FLG_SPEED | C.FLG_NONUT)[0]
            dpsi = _wrap180(lon_def - lon_non)
            # The geo/precession reduction is identical for the two flag sets,
            # so their difference is exactly the added Δψ (to the nutation-model
            # agreement between the backends, <0.001").
            assert abs((dpsi - expected) * 3600.0) < 0.01, (
                f"{name} ipl={ipl} jd={jd}: default-NONUT="
                f"{dpsi * 3600.0:.4f}\" vs Δψ={expected * 3600.0:.4f}\""
            )
    # And Δψ is a real, non-trivial term at JD_1900 (~17"), not the former ~0.
    assert abs(_dpsi_deg(JD_1900) * 3600.0) > 5.0


@SKIP_NO_LEB
@pytest.mark.parametrize("ipl", [C.CUPIDO, C.ISIS])
def test_uranian_default_speed_carries_nutation_rate(reader, skyfield_backend, ipl):
    """The of-date speed carries d(Δψ)/dt: default vs NONUT dlon differ by the
    nutation-in-longitude rate, and the two backends agree on it."""
    jd = JD_1900
    for get in (
        lambda j, f: _leb(reader, j, ipl, f),
        lambda j, f: _sky(j, ipl, f),
    ):
        dlon_def = get(jd, C.FLG_SPEED)[3]
        dlon_non = get(jd, C.FLG_SPEED | C.FLG_NONUT)[3]
        # d(Δψ)/dt is ~0.05-0.17"/day; the default speed carries it, NONUT does
        # not, so the difference is a real, non-zero nutation rate.
        rate = abs(_wrap180(dlon_def - dlon_non)) * 3600.0
        assert rate > 0.02
    # Cross-backend agreement on the default of-date speed (< 0.01"/day).
    a = _leb(reader, jd, ipl, C.FLG_SPEED)[3]
    s = _sky(jd, ipl, C.FLG_SPEED)[3]
    assert abs(_wrap180(a - s) * 3600.0) < 0.01


# ---------------------------------------------------------------------------
# H2/H3 — fictitious bodies (49-58) speeds are the derivative of position.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "ipl",
    [C.NIBIRU, C.VULCAN, C.WALDEMATH, C.WHITE_MOON, C.PROSERPINA],
)
def test_fictitious_speed_is_derivative_of_position(skyfield_backend, ipl):
    F = C.FLG_SPEED
    for jd in (JD_1900, JD_2000):
        rep = _sky(jd, ipl, F)[3]
        num = _deriv(lambda j: _sky(j, ipl, F), jd, 0, 0.005)
        assert abs(_wrap180(rep - num) * 3600.0) < 0.06, (ipl, jd)


def test_vulcan_helctr_speed_is_derivative_of_position(skyfield_backend):
    F = C.FLG_SPEED | C.FLG_HELCTR
    rep = _sky(JD_1900, C.VULCAN, F)[3]
    num = _deriv(lambda j: _sky(j, C.VULCAN, F), JD_1900, 0, 0.005)
    assert abs(_wrap180(rep - num) * 3600.0) < 0.05


# ---------------------------------------------------------------------------
# H4 — topocentric parallax on the hypothetical bodies.
# ---------------------------------------------------------------------------


def test_topocentric_parallax_applied_to_hypotheticals(skyfield_backend):
    le.set_topo(0.0, 0.0, 0.0)
    # Waldemath is close (~0.006 AU) so the parallax is large (~0.4 deg).
    geo = _sky(JD_2000, C.WALDEMATH, C.FLG_SPEED)
    topo = _sky(JD_2000, C.WALDEMATH, C.FLG_SPEED | C.FLG_TOPOCTR)
    sep = math.hypot(
        _wrap180(topo[0] - geo[0]) * math.cos(math.radians(geo[1])),
        topo[1] - geo[1],
    )
    # Real, of order the horizontal parallax; certainly not the silent no-op.
    assert sep * 3600.0 > 300.0
    # A distant Uranian has a tiny (but non-zero) parallax.
    geo_c = _sky(JD_2000, C.CUPIDO, C.FLG_SPEED)
    topo_c = _sky(JD_2000, C.CUPIDO, C.FLG_SPEED | C.FLG_TOPOCTR)
    assert 0.0 < abs(_wrap180(topo_c[0] - geo_c[0])) * 3600.0 < 1.0


def test_topocentric_hypothetical_requires_observer(skyfield_backend):
    from libephemeris.exceptions import Error

    le.close()  # clears any observer set by a prior test
    with pytest.raises(Error):
        _sky(JD_2000, C.WALDEMATH, C.FLG_SPEED | C.FLG_TOPOCTR)


# ---------------------------------------------------------------------------
# D1 — the TT and UT fast-path entry points are coherent.
# ---------------------------------------------------------------------------


@SKIP_NO_LEB
def test_fast_calc_tt_ut_coherent_at_historical_date(reader):
    """fast_calc_tt(jd_ut + ΔT) must reproduce fast_calc_ut(jd_ut)."""
    from libephemeris.time_utils import deltat

    le.set_topo(0.0, 51.5, 0.0)
    jd_ut = 2398000.0  # ~1854, where the LEB file ΔT table was stale
    jd_tt = jd_ut + deltat(jd_ut)
    a = fast_calc_ut(reader, jd_ut, C.MOON, C.FLG_SPEED | C.FLG_TOPOCTR)[0]
    b = fast_calc_tt(reader, jd_tt, C.MOON, C.FLG_SPEED | C.FLG_TOPOCTR)[0]
    assert abs(_wrap180(a[0] - b[0]) * 3600.0) < 0.01


# ---------------------------------------------------------------------------
# M1 — sidereal nod_aps returns native Python floats.
# ---------------------------------------------------------------------------


def test_sidereal_nod_aps_returns_native_floats(skyfield_backend):
    le.set_sid_mode(C.SIDM_LAHIRI)
    try:
        for ipl in (C.SUN, C.MARS):
            _n, _d, peri, aphe = le.nod_aps(
                JD_2000, ipl, C.NODBIT_MEAN, C.FLG_SIDEREAL | C.FLG_SPEED
            )
            assert type(peri[0]) is float
            assert type(aphe[0]) is float
    finally:
        le.set_sid_mode(C.SIDM_FAGAN_BRADLEY)


# ---------------------------------------------------------------------------
# RR1 — deflection rate through solar conjunction.
#
# The Pipeline-A velocity central-differences the apparent place on an
# analytically-extrapolated state; the gravitational-deflection term must be
# re-evaluated with the DEFLECTORS (and observer) at the sample epochs. With
# the deflectors frozen at jd_tt, slow targets near solar conjunction were
# biased by the Sun's apparent motion (~0.96 deg/day): Mars -0.34"/day at
# elongation 0.9 deg and +21.7"/day through the limb transit; Saturn
# +0.16"/day at conjunction. The reported speed must stay the derivative of
# the reported position (<= 0.01"/day) even there. The degenerate Sun-target
# case keeps its frozen delta (its own deflector geometry is singular).
# ---------------------------------------------------------------------------

# Solar-conjunction epochs (TT), elongations 0.004-0.91 deg, base tier.
JD_MARS_CONJ = 2460266.75  # elong 0.004 deg (limb transit)
JD_MARS_NEAR = 2460263.75  # elong 0.90 deg
JD_SATURN_CONJ = 2460369.40  # elong 0.006 deg


def _deriv5_lon(fn, jd, h=0.01):
    """5-point derivative of the reported longitude, wrap-aware."""
    l0 = fn(jd)[0]

    def rel(j):
        return _wrap180(fn(j)[0] - l0)

    return (
        -rel(jd + 2 * h) + 8 * rel(jd + h) - 8 * rel(jd - h) + rel(jd - 2 * h)
    ) / (12 * h)


@SKIP_NO_LEB
@pytest.mark.parametrize(
    "ipl, jd",
    [
        (C.MARS, JD_MARS_CONJ),
        (C.MARS, JD_MARS_NEAR),
        (C.SATURN, JD_SATURN_CONJ),
    ],
)
def test_conjunction_speed_is_derivative_leb(reader, ipl, jd):
    rep = _leb(reader, jd, ipl, C.FLG_SPEED)[3]
    num = _deriv5_lon(lambda j: _leb(reader, j, ipl, C.FLG_SPEED), jd)
    err_arcsec = abs(rep - num) * 3600.0
    assert err_arcsec < 0.01, (
        f"body {ipl} at JD {jd} (solar conjunction): reported dlon differs "
        f'from the derivative of the reported lon by {err_arcsec:.4f}"/day'
    )


@pytest.mark.parametrize(
    "ipl, jd",
    [(C.MARS, JD_MARS_CONJ), (C.SATURN, JD_SATURN_CONJ)],
)
def test_conjunction_speed_is_derivative_skyfield(skyfield_backend, ipl, jd):
    rep = _sky(jd, ipl, C.FLG_SPEED)[3]
    num = _deriv5_lon(lambda j: _sky(j, ipl, C.FLG_SPEED), jd)
    err_arcsec = abs(rep - num) * 3600.0
    assert err_arcsec < 0.01, (
        f"body {ipl} at JD {jd} (solar conjunction): reported dlon differs "
        f'from the derivative of the reported lon by {err_arcsec:.4f}"/day'
    )
