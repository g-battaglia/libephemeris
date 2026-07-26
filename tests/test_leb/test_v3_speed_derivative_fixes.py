"""Regression tests for the Round-2 "core" speed-derivative fixes.

Every fix here enforces the project's speed philosophy: the reported speed is
the true derivative of the reported position, and the LEB and Skyfield backends
agree on it. Covered:

* L1 True Node distance rate (``ddist``) is the derivative of the reconstructed
  osculating-node radius, not the retired proxy channel.
* L2/L3/L4/L5 node/apse LEB speeds match the Skyfield backend (J2000 equinox
  drift, nutation-in-longitude rate, TrueNode/OscuApog smoothing window,
  MeanApog latitude rate).
* H1 independently sourced hypothetical orbits carry finite position and
  derivative channels.
* H2 recognised fictitious identifiers without verified elements fail closed,
  while every reviewed primary-source model remains calculable.
* H3 topocentric parallax is applied to Harrington (and requires an observer,
  like the planets).
* D1 the TT and UT fast-path entry points are coherent.
* M1 sidereal ``nod_aps`` returns native Python floats.
* S1 ``SPEED3`` is the centered derivative of the reported positions.
"""

from __future__ import annotations

import math
import os

import pytest

import libephemeris as le
from libephemeris import constants as C
from libephemeris import fast_calc as fc
from libephemeris.exceptions import UnknownBodyError
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
# S1 — SPEED3 is a standard centered finite difference.
# ---------------------------------------------------------------------------


def test_speed3_centered_stencil_unwraps_spherical_longitude():
    """The rate uses P(t+h)-P(t-h), including across longitude zero."""
    current = (359.999, 1.5, 2.0, 0.0, 0.0, 0.0)
    previous = (359.998, 1.4998, 1.9999, 0.0, 0.0, 0.0)
    following = (0.0, 1.5002, 2.0001, 0.0, 0.0, 0.0)

    result = fc._apply_central_speed3_stencil(
        current, previous, following, C.FLG_SPEED3
    )

    assert result[:3] == current[:3]
    assert result[3:] == pytest.approx((10.0, 2.0, 1.0), abs=1e-9)


def test_speed3_centered_stencil_differentiates_cartesian_axes():
    """Cartesian SPEED3 differentiates each reported coordinate directly."""
    h = fc._SPEED3_STEP_DAYS
    current = (1.0, -2.0, 3.0, 0.0, 0.0, 0.0)
    expected = (0.25, -0.5, 2.0)
    previous = tuple(current[index] - h * expected[index] for index in range(3)) + (
        0.0,
        0.0,
        0.0,
    )
    following = tuple(current[index] + h * expected[index] for index in range(3)) + (
        0.0,
        0.0,
        0.0,
    )

    result = fc._apply_central_speed3_stencil(
        current, previous, following, C.FLG_SPEED3 | C.FLG_XYZ
    )

    assert result[:3] == current[:3]
    assert result[3:] == pytest.approx(expected, abs=1e-12)


@SKIP_NO_LEB
@pytest.mark.parametrize("ipl", [C.MEAN_NODE, C.MEAN_APOG])
def test_fast_speed3_equals_runtime_centered_position_derivative(reader, ipl):
    """The public fast path samples its own positions at exactly t±h."""
    h = fc._SPEED3_STEP_DAYS
    reported = fast_calc_tt(reader, JD_2000, ipl, C.FLG_SPEED3)[0]
    previous = fast_calc_tt(reader, JD_2000 - h, ipl, 0)[0]
    following = fast_calc_tt(reader, JD_2000 + h, ipl, 0)[0]
    expected_lon = _wrap180(following[0] - previous[0]) / (2.0 * h)
    expected = (
        expected_lon,
        (following[1] - previous[1]) / (2.0 * h),
        (following[2] - previous[2]) / (2.0 * h),
    )
    assert reported[3:] == pytest.approx(expected, abs=1e-10)


# ---------------------------------------------------------------------------
# L1 — True Node distance rate is the derivative of the reported distance.
# ---------------------------------------------------------------------------


@SKIP_NO_LEB
def test_true_node_ddist_is_derivative_of_reported_distance(reader):
    F = C.FLG_SPEED
    ddist_rep = _leb(reader, JD_1987, C.TRUE_NODE, F)[5]
    # Stencil matched to the TrueNode speed window (0.05 d): the finer window
    # resolves sub-daily osculating wiggles that a 0.5 d stencil averages out.
    ddist_num = _deriv(lambda j: _leb(reader, j, C.TRUE_NODE, F), JD_1987, 2, 0.05)
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
# H1 — Harrington of-date speeds pick up the general-precession rate.
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
def test_harrington_speed_leb_matches_skyfield(reader, skyfield_backend, flags):
    if not reader.has_body(C.HARRINGTON):
        pytest.skip("Harrington not in the LEB test file")
    jd = JD_1900
    a = _leb(reader, jd, C.HARRINGTON, flags)
    s = _sky(jd, C.HARRINGTON, flags)
    assert abs(_wrap180(a[3] - s[3]) * 3600.0) < 0.01


@SKIP_NO_LEB
def test_harrington_ofdate_speed_is_derivative_of_position(reader):
    if not reader.has_body(C.HARRINGTON):
        pytest.skip("Harrington not in the LEB test file")
    F = C.FLG_SPEED
    rep = _leb(reader, JD_1900, C.HARRINGTON, F)[3]
    num = _deriv(lambda j: _leb(reader, j, C.HARRINGTON, F), JD_1900, 0, 0.5)
    assert abs(_wrap180(rep - num) * 3600.0) < 0.05


# ---------------------------------------------------------------------------
# H1b — Harrington carries nutation in longitude (Δψ) in the of-date frame.
# The default result is on the true ecliptic; NONUT is on the mean ecliptic.
# ---------------------------------------------------------------------------


def _dpsi_deg(jd_tt: float) -> float:
    """Nutation in longitude Δψ at ``jd_tt``, in degrees."""
    from libephemeris.cache import get_cached_nutation

    dpsi_rad, _ = get_cached_nutation(jd_tt)
    return math.degrees(dpsi_rad)


@SKIP_NO_LEB
def test_harrington_default_minus_nonut_is_nutation(reader, skyfield_backend):
    """lon(default) - lon(NONUT) equals Δψ on both backends."""
    if not reader.has_body(C.HARRINGTON):
        pytest.skip("Harrington not in the LEB test file")
    ipl = C.HARRINGTON
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
                f'{dpsi * 3600.0:.4f}" vs Δψ={expected * 3600.0:.4f}"'
            )
    # Δψ is a real, non-trivial term at this epoch.
    assert abs(_dpsi_deg(JD_1900) * 3600.0) > 5.0


@SKIP_NO_LEB
def test_harrington_default_speed_carries_nutation_rate(reader, skyfield_backend):
    """The of-date speed carries d(Δψ)/dt: default vs NONUT dlon differ by the
    nutation-in-longitude rate, and the two backends agree on it."""
    if not reader.has_body(C.HARRINGTON):
        pytest.skip("Harrington not in the LEB test file")
    ipl = C.HARRINGTON
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
# H2 — reviewed source models calculate; unsupported IDs fail closed.
# ---------------------------------------------------------------------------


def test_harrington_speed_is_derivative_of_position(skyfield_backend):
    F = C.FLG_SPEED
    for jd in (JD_1900, JD_2000):
        rep = _sky(jd, C.HARRINGTON, F)[3]
        num = _deriv(lambda j: _sky(j, C.HARRINGTON, F), jd, 0, 0.005)
        assert abs(_wrap180(rep - num) * 3600.0) < 0.06, jd


def test_harrington_j2000_heliocentric_state_rates(skyfield_backend):
    """Both angular rates survive the J2000/of-date frame round trip."""
    flags = C.FLG_HELCTR | C.FLG_TRUEPOS | C.FLG_J2000 | C.FLG_SPEED
    for component in (0, 1):
        reported = _sky(JD_1900, C.HARRINGTON, flags)[component + 3]
        numerical = _deriv(
            lambda jd: _sky(jd, C.HARRINGTON, flags), JD_1900, component, 0.1
        )
        delta = (
            _wrap180(reported - numerical) if component == 0 else reported - numerical
        )
        assert abs(delta * 3600.0) < 0.01, component


REVIEWED_FICTITIOUS_IDS = (*range(C.CUPIDO, C.POSEIDON + 1), 48, *range(50, 59))
UNVERIFIED_FICTITIOUS_IDS = (49,)


@pytest.mark.parametrize("ipl", REVIEWED_FICTITIOUS_IDS)
def test_reviewed_fictitious_ids_return_finite_state(skyfield_backend, ipl):
    """Each primary-source-backed model returns all six native channels."""
    state = _sky(JD_2000, ipl, C.FLG_SPEED)
    assert len(state) == 6
    assert all(math.isfinite(value) for value in state)


@pytest.mark.parametrize("ipl", UNVERIFIED_FICTITIOUS_IDS)
def test_unverified_fictitious_ids_fail_closed(skyfield_backend, ipl):
    """Recognised IDs remain named but cannot use unverified element sets."""
    with pytest.raises(UnknownBodyError) as raised:
        _sky(JD_2000, ipl, C.FLG_SPEED)
    assert raised.value.body_id == ipl


# ---------------------------------------------------------------------------
# H3 — topocentric parallax on Harrington.
# ---------------------------------------------------------------------------


def test_topocentric_parallax_applied_to_harrington(skyfield_backend):
    le.set_topo(0.0, 0.0, 0.0)
    # A distant orbit has a tiny but non-zero topocentric displacement.
    geo = _sky(JD_2000, C.HARRINGTON, C.FLG_SPEED)
    topo = _sky(JD_2000, C.HARRINGTON, C.FLG_SPEED | C.FLG_TOPOCTR)
    assert 0.0 < abs(_wrap180(topo[0] - geo[0])) * 3600.0 < 1.0


def test_topocentric_harrington_requires_observer(skyfield_backend):
    from libephemeris.exceptions import Error

    le.close()  # clears any observer set by a prior test
    with pytest.raises(Error):
        _sky(JD_2000, C.HARRINGTON, C.FLG_SPEED | C.FLG_TOPOCTR)


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

    return (-rel(jd + 2 * h) + 8 * rel(jd + h) - 8 * rel(jd - h) + rel(jd - 2 * h)) / (
        12 * h
    )


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


# ---------------------------------------------------------------------------
# H1-F1 (zero-bug Round H) — of-date SIDEREAL speed of the ecliptic-direct
# bodies (Pipeline B: nodes/apogees) drops the nutation-in-longitude rate
# dΔψ/dt, mirroring the position correction (true ayanamsha = mean + Δψ is
# applied to ALL bodies). Pre-fix the LEB speed was off by up to ~0.23"/day
# (the nutation rate), diverging from both the reported-position derivative
# and the Skyfield backend.
# ---------------------------------------------------------------------------

# Dates near nutation-rate extrema (|dΔψ/dt| ~0.15-0.23"/day) so a regression
# cannot hide in a small-rate epoch.
JD_NUT_PEAK = 2452997.0
JD_2010 = 2455197.0

_SID_LAHIRI = {"sid_mode": C.SIDM_LAHIRI, "sid_t0": 0.0, "sid_ayan_t0": 0.0}

_ECLIPTIC_DIRECT = [
    C.MEAN_NODE,
    C.TRUE_NODE,
    C.MEAN_APOG,
    C.OSCU_APOG,
    C.INTP_APOG,
    C.INTP_PERG,
]


@SKIP_NO_LEB
@pytest.mark.parametrize("jd", [JD_NUT_PEAK, JD_2010])
@pytest.mark.parametrize("ipl", _ECLIPTIC_DIRECT)
def test_ecliptic_direct_sidereal_speed_is_derivative(reader, ipl, jd):
    if not reader.has_body(ipl):
        pytest.skip(f"body {ipl} not in the base LEB test file")
    F = C.FLG_SPEED | C.FLG_SIDEREAL
    rep = fast_calc_tt(reader, jd, ipl, F, **_SID_LAHIRI)[0][3]
    num = _deriv(
        lambda j: fast_calc_tt(reader, j, ipl, C.FLG_SIDEREAL, **_SID_LAHIRI)[0],
        jd,
        0,
        0.05,
    )
    sid_err = (rep - num) * 3600.0
    # The interpolated apsides carry a small pre-existing smoothing residual
    # in their speed channel that is identical in the tropical frame (and in
    # the Skyfield backend). Measure the tropical baseline and assert the
    # SIDEREAL path adds no extra drift on top of it — the nutation-rate bug
    # showed up precisely as that extra drift (up to ~0.23"/day).
    rep_t = fast_calc_tt(reader, jd, ipl, C.FLG_SPEED)[0][3]
    num_t = _deriv(lambda j: fast_calc_tt(reader, j, ipl, 0)[0], jd, 0, 0.05)
    trop_err = (rep_t - num_t) * 3600.0
    assert abs(sid_err - trop_err) < 0.02, (
        f'body {ipl} at JD {jd}: sidereal self-error {sid_err:+.4f}"/day vs '
        f'tropical baseline {trop_err:+.4f}"/day - extra sidereal drift'
    )


@SKIP_NO_LEB
@pytest.mark.parametrize("extra", [0, C.FLG_J2000])
@pytest.mark.parametrize("ipl", [C.MEAN_NODE, C.TRUE_NODE, C.MEAN_APOG, C.OSCU_APOG])
def test_ecliptic_direct_sidereal_speed_leb_matches_skyfield(
    reader, skyfield_backend, ipl, extra
):
    # extra=J2000 exercises the deferred rebuild: its forward velocity sample
    # must precess from its own epoch, or the LEB speed lands off by exactly
    # the general-precession rate (~0.1377"/day) versus Skyfield.
    F = C.FLG_SPEED | C.FLG_SIDEREAL | extra
    a = fast_calc_tt(reader, JD_NUT_PEAK, ipl, F, **_SID_LAHIRI)[0]
    le.set_sid_mode(C.SIDM_LAHIRI)
    try:
        s = le.calc(JD_NUT_PEAK, ipl, F)[0]
    finally:
        le.set_sid_mode(C.SIDM_FAGAN_BRADLEY)
    assert abs(_wrap180(a[3] - s[3]) * 3600.0) < 0.02


@SKIP_NO_LEB
@pytest.mark.parametrize("extra", [C.FLG_J2000, C.FLG_NONUT])
def test_mean_ayanamsha_sidereal_speed_has_no_nutation_term(reader, extra):
    """J2000/NONUT sidereal uses the MEAN ayanamsha (no Δψ), so the fix must
    NOT leak dΔψ/dt into those speeds: they stay the derivative of the
    reported position."""
    F = C.FLG_SPEED | C.FLG_SIDEREAL | extra
    rep = fast_calc_tt(reader, JD_NUT_PEAK, C.MEAN_NODE, F, **_SID_LAHIRI)[0][3]
    num = _deriv(
        lambda j: fast_calc_tt(
            reader, j, C.MEAN_NODE, C.FLG_SIDEREAL | extra, **_SID_LAHIRI
        )[0],
        JD_NUT_PEAK,
        0,
        0.05,
    )
    assert abs(rep - num) * 3600.0 < 0.02
