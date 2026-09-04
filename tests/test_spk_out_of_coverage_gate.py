"""Regression tests for the strict-precision gate with a REGISTERED SPK kernel.

Background
----------
Strict precision mode refuses the low-accuracy Keplerian fallback for bodies
that need SPK-grade data. The gate must decide "does a source better than
Keplerian exist for this epoch?" -- and that decision must NOT depend on an
irrelevant network attempt.

The bug these tests lock in: with a kernel already REGISTERED for a body but
the requested epoch OUTSIDE its coverage window, strict mode used to

  * raise SPKRequiredError when auto-download was OFF, yet
  * serve the very same epoch from ASSIST when auto-download was ON

purely because the (doomed) download attempt flipped an internal flag that
exempted the raise. It also raised at epochs sitting exactly on the coverage
border, and the raised message told the user to "register an SPK" even though
one was already registered.

After the fix:
  * ASSIST availability is a purely local probe -- auto-download ON and OFF
    behave identically.
  * A registered-but-out-of-coverage epoch is served by ASSIST when ASSIST
    data is present, regardless of the auto-download setting.
  * When no source at all can serve the epoch, the raise message names the
    registered kernel and its coverage instead of "register an SPK".

The first class is fully hermetic (all I/O mocked -- no network, no SPK cache,
no multi-GB ASSIST data). The second class is an opt-in integration check that
runs only when the real Chiron kernel and ASSIST data happen to be installed.
"""

from __future__ import annotations

import glob
import os

import pytest
from unittest.mock import patch

import libephemeris as eph
from libephemeris import planets, spk, state, tracing
from libephemeris.constants import CHIRON, FLG_SPEED
from libephemeris.exceptions import EphemerisRangeError, SPKRequiredError
from libephemeris.rebound_integration import PropagationResult

# A registered kernel whose coverage does NOT include the probed epoch.
FAKE_SPK = "/nonexistent/2060_fake_coverage.bsp"
COV_START = 2305447.5
COV_END = 2634166.5
OUT_OF_COV_JD = COV_END + 500.0  # comfortably past the end border


@pytest.fixture(autouse=True)
def _isolate(monkeypatch):
    """Reset all relevant global state and neutralise the LEB fast path."""
    from libephemeris.rebound_integration import reset_assist_data_cache

    state._SPK_BODY_MAP.clear()
    monkeypatch.delenv("LIBEPHEMERIS_LEB", raising=False)
    monkeypatch.delenv("LIBEPHEMERIS_MODE", raising=False)
    monkeypatch.setattr(state, "_LEB_FILE", None)
    monkeypatch.setattr(state, "_LEB_READER", None)
    monkeypatch.setattr(state, "_discover_leb_file", lambda: None)
    eph.set_calc_mode("auto")
    eph.set_strict_precision(None)
    eph.set_auto_spk_download(None)
    reset_assist_data_cache()
    yield
    state._SPK_BODY_MAP.clear()
    eph.set_strict_precision(None)
    eph.set_auto_spk_download(None)
    reset_assist_data_cache()


def _mock_assist_result(jd_tt: float) -> PropagationResult:
    """A deterministic heliocentric state for a body at ~15 AU."""
    return PropagationResult(
        x=-12.0, y=8.0, z=-2.0, vx=0.0006, vy=-0.0009, vz=0.00005, jd_tt=jd_tt
    )


class _RegisteredOutOfCoverage:
    """Mixin: register a kernel for Chiron that never covers the epoch.

    ``get_spk_type21_target`` returns None and ``calc_spk_body_position``
    raises EphemerisRangeError -- exactly what a real type-21 kernel does for
    an out-of-coverage epoch -- while ``get_spk_body_info`` reports the kernel
    as registered.
    """

    def _enter(self, stack):
        state._SPK_BODY_MAP[CHIRON] = (FAKE_SPK, 20002060)
        stack.enter_context(
            patch("libephemeris.spk.get_spk_type21_target", return_value=None)
        )
        stack.enter_context(
            patch(
                "libephemeris.spk.calc_spk_body_position",
                side_effect=EphemerisRangeError("outside coverage"),
            )
        )
        stack.enter_context(
            patch(
                "libephemeris.spk.get_spk_coverage",
                return_value=(COV_START, COV_END),
            )
        )


class TestStrictGateSourceCheck(_RegisteredOutOfCoverage):
    """Hermetic tests of the strict gate's 'is there a better source?' logic."""

    def _run_out_of_coverage(self, autodl: bool, assist: bool):
        """Call Chiron out-of-coverage; return (source, position) or raise."""
        from contextlib import ExitStack

        with ExitStack() as stack:
            self._enter(stack)
            stack.enter_context(
                patch(
                    "libephemeris.rebound_integration.check_assist_data_available",
                    return_value=assist,
                )
            )
            stack.enter_context(
                patch(
                    "libephemeris.rebound_integration.propagate_orbit_assist",
                    side_effect=lambda el, a, b, **k: _mock_assist_result(b),
                )
            )
            # Neutralise the network entirely: a real auto-download must never
            # be reached, and if it is (autodl=True) it resolves nothing.
            stack.enter_context(
                patch(
                    "libephemeris.planets._try_auto_spk_download",
                    return_value=None,
                )
            )
            # Hermetic: the gate's source logic is under test, so the LEB
            # asteroid channels must not serve the request either (see
            # _calc() above for the same reasoning).
            from libephemeris.state import get_calc_mode, set_calc_mode

            saved_mode = get_calc_mode()
            set_calc_mode("skyfield")
            stack.callback(set_calc_mode, saved_mode)
            eph.set_strict_precision(True)
            eph.set_auto_spk_download(autodl)
            tok = tracing.start_tracing()
            try:
                pos, _ = eph.calc_ut(OUT_OF_COV_JD, CHIRON, FLG_SPEED)
                return tracing.get_trace_results().get(CHIRON), pos
            finally:
                tracing._trace_data.reset(tok)

    def test_assist_served_regardless_of_auto_download(self):
        """With ASSIST present, auto-download ON and OFF are IDENTICAL."""
        src_off, pos_off = self._run_out_of_coverage(autodl=False, assist=True)
        src_on, pos_on = self._run_out_of_coverage(autodl=True, assist=True)

        assert src_off == "ASSIST"
        assert src_on == "ASSIST"
        # Byte-identical: the source no longer depends on the network attempt.
        assert pos_off == pos_on
        assert 0.0 <= pos_off[0] < 360.0

    def test_no_raise_at_out_of_coverage_when_assist_available(self):
        """Auto-download OFF + ASSIST present must NOT raise (the bug)."""
        src, pos = self._run_out_of_coverage(autodl=False, assist=True)
        assert src == "ASSIST"
        assert pos[2] > 0.0

    def test_out_of_coverage_message_is_honest(self):
        """No ASSIST -> raise, but the message must cite the registered kernel."""
        with pytest.raises(SPKRequiredError) as exc:
            self._run_out_of_coverage(autodl=False, assist=False)

        msg = str(exc.value)
        assert exc.value.body_id == CHIRON
        assert exc.value.body_name == "Chiron"
        # Truthful: it does NOT tell the user to register an SPK...
        assert "download_and_register_spk" not in msg
        # ...it says the kernel is registered but the epoch is out of coverage.
        assert "does not cover" in msg
        assert "already registered" in msg
        assert f"{COV_START:.1f}" in msg and f"{COV_END:.1f}" in msg
        assert os.path.basename(FAKE_SPK) in msg

    def test_no_kernel_message_still_says_register(self):
        """With no kernel registered and no ASSIST, the classic message stands."""
        from libephemeris.state import get_calc_mode, set_calc_mode

        saved_mode = get_calc_mode()
        set_calc_mode("skyfield")  # exclude the LEB asteroid channels too
        try:
            with patch(
                "libephemeris.rebound_integration.check_assist_data_available",
                return_value=False,
            ):
                eph.set_strict_precision(True)
                eph.set_auto_spk_download(False)
                with pytest.raises(SPKRequiredError) as exc:
                    eph.calc_ut(2451545.0, CHIRON, FLG_SPEED)
        finally:
            set_calc_mode(saved_mode)

        msg = str(exc.value)
        assert "download_and_register_spk" in msg
        assert "does not cover" not in msg


class TestStrictSourceHelper:
    """Unit tests for the ASSIST-availability helper used by the gate."""

    def test_true_for_curated_body_with_assist(self):
        with patch(
            "libephemeris.rebound_integration.check_assist_data_available",
            return_value=True,
        ):
            assert planets._strict_source_better_than_keplerian(CHIRON) is True

    def test_false_without_assist(self):
        with patch(
            "libephemeris.rebound_integration.check_assist_data_available",
            return_value=False,
        ):
            assert planets._strict_source_better_than_keplerian(CHIRON) is False

    def test_false_for_body_without_curated_elements(self):
        from libephemeris.constants import SUN

        # SUN has no ASSIST initial conditions -> never a better source here.
        with patch(
            "libephemeris.rebound_integration.check_assist_data_available",
            return_value=True,
        ):
            assert planets._strict_source_better_than_keplerian(SUN) is False


# ---------------------------------------------------------------------------
# Opt-in integration check against the real Chiron kernel + real ASSIST data.
# ---------------------------------------------------------------------------


def _find_chiron_kernel():
    pattern = os.path.expanduser("~/.libephemeris/spk/2060_*.bsp")
    files = sorted(
        f for f in glob.glob(pattern) if not os.path.basename(f).startswith("._")
    )
    return files[-1] if files else None


def _real_assist_available() -> bool:
    from libephemeris.rebound_integration import check_assist_data_available

    return check_assist_data_available()


_KERNEL = _find_chiron_kernel()
_HAVE_REAL = _KERNEL is not None and _real_assist_available()


@pytest.mark.skipif(
    not _HAVE_REAL,
    reason="requires a cached Chiron SPK kernel and installed ASSIST data",
)
class TestStrictGateRealKernel:
    """End-to-end: real registered kernel, real ASSIST, no network."""

    def _register(self):
        targets = spk._get_spk_targets(_KERNEL)
        naif = next((c for c in (2002060, 20002060) if c in targets), None)
        assert naif is not None
        spk.register_spk_body(CHIRON, _KERNEL, naif)
        return spk.get_spk_coverage(_KERNEL)

    def _calc(self, jd_ut, autodl, monkeypatch):
        # Block real network downloads for the auto-download=True comparison.
        monkeypatch.setattr(
            planets, "_try_auto_spk_download", lambda t, ipl, iflag: None
        )
        # This gate is about the SPK source specifically: pin the mode so an
        # installed LEB (whose asteroid channels can serve Chiron) does not
        # win the dispatch and make the trace report LEB instead of SPK.
        from libephemeris.state import set_calc_mode

        set_calc_mode("skyfield")
        eph.set_strict_precision(True)
        eph.set_auto_spk_download(autodl)
        tok = tracing.start_tracing()
        try:
            pos, _ = eph.calc_ut(jd_ut, CHIRON, FLG_SPEED)
            return tracing.get_trace_results().get(CHIRON), pos
        finally:
            tracing._trace_data.reset(tok)

    def test_in_coverage_uses_kernel(self, monkeypatch):
        start, end = self._register()
        src, pos = self._calc(
            (start + end) / 2.0, autodl=False, monkeypatch=monkeypatch
        )
        assert src == "SPK"
        assert 0.0 <= pos[0] < 360.0

    def test_out_of_coverage_autodl_on_equals_off(self, monkeypatch):
        start, end = self._register()
        src_off, pos_off = self._calc(
            end + 500.0, autodl=False, monkeypatch=monkeypatch
        )
        self._register()
        src_on, pos_on = self._calc(end + 500.0, autodl=True, monkeypatch=monkeypatch)
        assert src_off == "ASSIST"
        assert src_on == "ASSIST"
        # Identical value -- no dependency on the (blocked) download attempt.
        assert abs(pos_off[0] - pos_on[0]) < 1e-9

    def test_spk_assist_continuity_at_handoff(self, monkeypatch):
        """Value must be continuous (<0.05") where the kernel hands off to ASSIST."""
        from libephemeris.state import get_planets, get_timescale

        start, end = self._register()
        ts = get_timescale()

        def spk_lon(jd_tt):
            # Outside the usable coverage the direct path raises the typed
            # range error (the light-time band at the start is excluded from
            # the reported coverage, never silently degraded).
            try:
                r = spk.calc_spk_body_position(ts.tt_jd(jd_tt), CHIRON, FLG_SPEED)
            except EphemerisRangeError:
                return None
            return None if r is None else r[0]

        def assist_lon(jd_tt):
            st = planets._assist_state_at(jd_tt, CHIRON)
            lon, _, _ = planets._assist_position_from_state(
                st, jd_tt, FLG_SPEED, get_planets()
            )
            return lon

        # Find the first TT epoch at or after the reported start that the
        # kernel serves. get_spk_coverage() already reports the usable start
        # (light-time band excluded), so this is normally the start itself;
        # the scan only absorbs the TT/TDB offset at the boundary.
        handoff = None
        for k in range(0, 400):
            jd = start + k * 0.005
            if jd > start + 2.0:
                break
            if spk_lon(jd) is not None:
                handoff = jd
                break
        assert handoff is not None

        for jd in (handoff, end):
            s = spk_lon(jd)
            a = assist_lon(jd)
            assert s is not None
            diff_arcsec = abs(((s - a + 180.0) % 360.0) - 180.0) * 3600.0
            assert diff_arcsec < 0.05, f'jd={jd}: {diff_arcsec:.4f}" > 0.05"'
