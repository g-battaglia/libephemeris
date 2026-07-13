"""Independent invariants for two v3 frame/context fixes.

Fix 1 — FLG_J2000 velocity is precessed consistently with the position for
analytic bodies (lunar nodes, Lilith, apsides, and Harrington). Before the
fix, ``calc_ut(jd, MEAN_NODE, FLG_SPEED | FLG_J2000)`` returned the *of-date*
longitude speed (the velocity was left in the of-date frame while only the
position was precessed), and ``FLG_J2000 | FLG_EQUATORIAL`` rotated an of-date
velocity through the J2000 obliquity (a frame-mixed RA/Dec speed). The J2000
speed must instead be the time derivative of the J2000 position, differing from
the of-date speed by the general-precession rate (~0.14"/day for the node).

Fix 2 — The Context Horizons dispatch honors this context's per-instance
user-defined ayanamsha (SIDM_USER t0 / ayan_t0). Previously it passed only the
sidereal *mode* and ran the call without swapping in the context's sidereal
state, so ``get_ayanamsa`` read the reference epoch and offset from the module
globals instead of the context, and the read was unlocked.

No output from the reference API is stored in this module.
"""

from __future__ import annotations

import pytest

import libephemeris as le
from libephemeris import state
from libephemeris.constants import (
    FLG_EQUATORIAL,
    FLG_J2000,
    FLG_SIDEREAL,
    FLG_SPEED,
    MEAN_APOG,
    MEAN_NODE,
    SIDM_USER,
)
from libephemeris.context import EphemerisContext

# --------------------------------------------------------------------------- #
# Shared constants
# --------------------------------------------------------------------------- #

JD = 2451545.0  # J2000.0 epoch (UT); a stable, well-conditioned test date.

ARCSEC_PER_DEG = 3600.0


@pytest.fixture()
def skyfield_backend():
    """Force the Skyfield backend (the code path fixed in planets.py).

    The LEB fast path lives in fast_calc.py and is out of scope for these
    fixes, so the tests must exercise the Skyfield backend explicitly. Skips
    when the DE440 ephemeris is not available locally.
    """
    original = le.get_calc_mode()
    le.set_calc_mode("skyfield")
    try:
        from libephemeris.constants import SUN

        le.calc_ut(JD, SUN, 0)  # probe: raises if DE440 is unavailable
    except Exception as exc:  # pragma: no cover - environment dependent
        le.set_calc_mode(original)
        pytest.skip(f"Skyfield backend/DE440 unavailable: {exc!r}")
    yield
    le.set_calc_mode(original)


def _lon_diff_deg(a: float, b: float) -> float:
    """Shortest-path a - b in degrees, wrapped to (-180, 180]."""
    d = (a - b) % 360.0
    if d > 180.0:
        d -= 360.0
    return d


# --------------------------------------------------------------------------- #
# Fix 1 — FLG_J2000 velocity precession
# --------------------------------------------------------------------------- #


class TestJ2000NodeSpeed:
    def test_j2000_speed_differs_from_ofdate_by_precession_rate(self, skyfield_backend):
        """J2000 node speed must NOT equal the of-date speed.

        The difference is the general-precession rate in longitude
        (~0.1377"/day), i.e. the equinox drift the of-date frame carries and
        the fixed J2000 frame does not.
        """
        ofdate = le.calc_ut(JD, MEAN_NODE, FLG_SPEED)[0][3]
        j2000 = le.calc_ut(JD, MEAN_NODE, FLG_SPEED | FLG_J2000)[0][3]

        diff_arcsec = abs(j2000 - ofdate) * ARCSEC_PER_DEG
        # Precession rate in longitude ~= 50.29"/yr ~= 0.1377"/day. The of-date
        # bug would give ~0.0 here.
        assert 0.10 < diff_arcsec < 0.16, (
            f'J2000 vs of-date node speed differs by {diff_arcsec:.4f}"/day; '
            f'expected ~0.1377"/day (general precession)'
        )

    @pytest.mark.parametrize("ipl", [MEAN_NODE, MEAN_APOG])
    @pytest.mark.parametrize(
        "flags",
        [FLG_SPEED | FLG_J2000, FLG_SPEED | FLG_J2000 | FLG_EQUATORIAL],
    )
    def test_j2000_speed_is_derivative_of_j2000_position(
        self, skyfield_backend, ipl, flags
    ):
        """Reference-free: the returned J2000 speed is d(J2000 position)/dt.

        This is the frame-consistency invariant the fix restores and does not
        depend on any external reference.
        """
        dt = 0.05
        pos_m = le.calc_ut(JD - dt, ipl, flags)[0][0]
        pos_p = le.calc_ut(JD + dt, ipl, flags)[0][0]
        numerical = _lon_diff_deg(pos_p, pos_m) / (2.0 * dt)
        returned = le.calc_ut(JD, ipl, flags)[0][3]
        err_arcsec = abs(returned - numerical) * ARCSEC_PER_DEG
        assert err_arcsec < 0.02, (
            f"body {ipl} flags {flags}: returned speed {returned:.10f} vs "
            f'numerical {numerical:.10f}: {err_arcsec:.4f}"/day'
        )


# --------------------------------------------------------------------------- #
# Fix 2 — Context + Horizons + SIDM_USER
# --------------------------------------------------------------------------- #


class TestContextHorizonsSidmUser:
    def test_swap_installs_context_sidereal_state(self, monkeypatch):
        """The context's SIDM_USER t0/ayan_t0 are swapped in for the call.

        The Horizons dispatch must run under _swapped_context_state so
        get_ayanamsa (which reads the module-level sidereal globals for
        SIDM_USER) sees THIS context's reference epoch and offset — not the
        global ones — and reads them under the swap lock.
        """
        import libephemeris.horizons_backend as hb

        # Global sidereal state deliberately different from the context's.
        monkeypatch.setattr(state, "_SIDEREAL_MODE", SIDM_USER, raising=False)
        monkeypatch.setattr(state, "_SIDEREAL_T0", 2451545.0, raising=False)
        monkeypatch.setattr(state, "_SIDEREAL_AYAN_T0", 99.0, raising=False)

        ctx = EphemerisContext()
        ctx.set_sid_mode(SIDM_USER, t0=2440000.0, ayan_t0=12.3456)

        # Force the Horizons path: a non-None client, no LEB reader anywhere.
        monkeypatch.setattr(state, "get_horizons_client", lambda: object())
        monkeypatch.setattr(state, "get_leb_reader", lambda: None)
        monkeypatch.setattr(ctx, "get_leb_reader", lambda: None)

        captured: dict[str, float] = {}

        def fake_horizons_calc_ut(client, jd, ipl, iflag, sid_mode=None):
            # Read the GLOBALS as get_ayanamsa would, from inside the call.
            captured["t0"] = state._SIDEREAL_T0
            captured["ayan_t0"] = state._SIDEREAL_AYAN_T0
            captured["mode"] = state._SIDEREAL_MODE
            captured["sid_mode_arg"] = sid_mode
            return (10.0, 0.0, 1.0, -0.05, 0.0, 0.0), iflag

        monkeypatch.setattr(hb, "horizons_calc_ut", fake_horizons_calc_ut)

        ctx.calc_ut(JD, MEAN_NODE, FLG_SIDEREAL | FLG_SPEED)

        assert captured["t0"] == 2440000.0
        assert captured["ayan_t0"] == 12.3456
        assert captured["mode"] == SIDM_USER
        assert captured["sid_mode_arg"] == SIDM_USER

        # State restored after the call.
        assert state._SIDEREAL_T0 == 2451545.0
        assert state._SIDEREAL_AYAN_T0 == 99.0

    def test_context_ayanamsa_shifts_sidereal_longitude(self, monkeypatch):
        """Behavioral: two contexts differing only in ayan_t0 differ in output.

        Mean Node is analytic, so horizons_calc_ut computes it with no network.
        At the reference epoch (t0 == JD) the ayanamsha equals ayan_t0, so a
        10 deg difference in ayan_t0 shifts the sidereal longitude by 10 deg.
        """
        monkeypatch.setattr(state, "get_horizons_client", lambda: object())
        monkeypatch.setattr(state, "get_leb_reader", lambda: None)

        def make_ctx(ayan_t0: float) -> EphemerisContext:
            ctx = EphemerisContext()
            ctx.set_sid_mode(SIDM_USER, t0=JD, ayan_t0=ayan_t0)
            monkeypatch.setattr(ctx, "get_leb_reader", lambda: None)
            return ctx

        lon_a = make_ctx(10.0).calc_ut(JD, MEAN_NODE, FLG_SIDEREAL)[0][0]
        lon_b = make_ctx(20.0).calc_ut(JD, MEAN_NODE, FLG_SIDEREAL)[0][0]

        # sidereal_lon = tropical_lon - ayanamsha; larger ayan_t0 -> smaller lon.
        assert abs(_lon_diff_deg(lon_a, lon_b) - 10.0) < 1e-6

    def test_valueerror_from_horizons_falls_back(self, monkeypatch, skyfield_backend):
        """A ValueError from the Horizons path falls back to Skyfield.

        The widened except must not let a ValueError/OSError escape the context
        calc (the fallback must produce a valid Skyfield result).
        """
        import libephemeris.horizons_backend as hb

        monkeypatch.setattr(state, "get_horizons_client", lambda: object())
        monkeypatch.setattr(state, "get_leb_reader", lambda: None)

        def boom(client, jd, ipl, iflag, sid_mode=None):
            raise ValueError("simulated HELCTR/velocity path failure")

        monkeypatch.setattr(hb, "horizons_calc_ut", boom)

        ctx = EphemerisContext()
        monkeypatch.setattr(ctx, "get_leb_reader", lambda: None)

        # Must not raise: falls through to the Skyfield backend.
        pos, _ = ctx.calc_ut(JD, MEAN_NODE, FLG_SPEED)
        assert 0.0 <= pos[0] < 360.0
