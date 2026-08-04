"""Deterministic tests for the two global-state guards.

- planets._calc_body_race_safe: the bounded close()/reload race retry.
  Exercised here with a monkeypatched _calc_body, so the retry budget and
  the propagate-vs-retry predicate are pinned without real concurrency.
- eclipse._observer_scope: the save/restore of the process-global observer
  around event searches, including restore on mid-search exceptions.
"""

from __future__ import annotations

import pytest

import libephemeris as swe
from libephemeris import planets, state
from libephemeris.eclipse import _observer_scope
from libephemeris.planets import _MAX_CLOSE_RACE_RETRIES, _calc_body_race_safe


@pytest.fixture(autouse=True)
def _clean_state():
    yield
    swe.set_topo(0.0, 0.0, 0.0)
    swe.close()


class _Flaky:
    """Callable raising `exc` for the first `failures` calls, then returning."""

    def __init__(self, exc: BaseException, failures: int):
        self.exc = exc
        self.failures = failures
        self.calls = 0

    def __call__(self, t, planet, calc_iflag):
        self.calls += 1
        if self.calls <= self.failures:
            raise self.exc
        return ("ok", self.calls)


class TestCalcBodyRaceRetry:
    @pytest.mark.unit
    def test_transient_keyerror_retries_then_succeeds(self, monkeypatch):
        flaky = _Flaky(KeyError("kernel 'de440.bsp' is missing 'EARTH'"), 3)
        monkeypatch.setattr(planets, "_calc_body", flaky)
        assert _calc_body_race_safe(None, 4, 0) == ("ok", 4)
        assert flaky.calls == 4

    @pytest.mark.unit
    def test_budget_exhaustion_propagates_the_error(self, monkeypatch):
        flaky = _Flaky(KeyError(0), _MAX_CLOSE_RACE_RETRIES + 1)
        monkeypatch.setattr(planets, "_calc_body", flaky)
        with pytest.raises(KeyError):
            _calc_body_race_safe(None, 4, 0)
        # 1 initial call + _MAX_CLOSE_RACE_RETRIES retries, no more.
        assert flaky.calls == _MAX_CLOSE_RACE_RETRIES + 1

    @pytest.mark.unit
    def test_last_budget_slot_still_succeeds(self, monkeypatch):
        flaky = _Flaky(KeyError(0), _MAX_CLOSE_RACE_RETRIES)
        monkeypatch.setattr(planets, "_calc_body", flaky)
        assert _calc_body_race_safe(None, 4, 0) == (
            "ok",
            _MAX_CLOSE_RACE_RETRIES + 1,
        )

    @pytest.mark.unit
    def test_closed_file_valueerror_is_retried(self, monkeypatch):
        flaky = _Flaky(ValueError("seek of closed file"), 1)
        monkeypatch.setattr(planets, "_calc_body", flaky)
        assert _calc_body_race_safe(None, 4, 0) == ("ok", 2)

    @pytest.mark.unit
    def test_non_race_valueerror_propagates_immediately(self, monkeypatch):
        flaky = _Flaky(ValueError("Invalid Time"), 5)
        monkeypatch.setattr(planets, "_calc_body", flaky)
        with pytest.raises(ValueError, match="Invalid Time"):
            _calc_body_race_safe(None, 4, 0)
        assert flaky.calls == 1  # no retry for a non-race shape


class TestObserverScopeRestore:
    OUTER = (12.5, 41.9, 100.0)
    INNER = (-70.0, -33.5, 520.0)

    @staticmethod
    def _topo_tuple():
        topo = state.get_topo()
        if topo is None:
            return None
        return (
            topo.longitude.degrees,
            topo.latitude.degrees,
            topo.elevation.m,
        )

    @pytest.mark.unit
    def test_normal_exit_restores_previous_observer(self):
        swe.set_topo(*self.OUTER)
        with _observer_scope(*self.INNER):
            assert self._topo_tuple() == pytest.approx(self.INNER)
        assert self._topo_tuple() == pytest.approx(self.OUTER)

    @pytest.mark.unit
    def test_exception_mid_search_restores_previous_observer(self):
        swe.set_topo(*self.OUTER)
        with pytest.raises(RuntimeError, match="mid-search"):
            with _observer_scope(*self.INNER):
                assert self._topo_tuple() == pytest.approx(self.INNER)
                raise RuntimeError("mid-search failure")
        assert self._topo_tuple() == pytest.approx(self.OUTER)

    @pytest.mark.unit
    def test_unset_observer_round_trips(self):
        # With no observer set, the scope must restore the unset state,
        # not leave the search position behind.
        swe.close()
        assert state.get_topo() is None
        with _observer_scope(*self.INNER):
            assert self._topo_tuple() == pytest.approx(self.INNER)
        assert state.get_topo() is None
