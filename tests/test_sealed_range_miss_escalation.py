"""Unit tests for fast_calc._escalate_sealed_range_miss.

The sealed-mode escalation names the REQUESTED body when a curated minor
body misses its LEB window and the Keplerian fall-through cannot run
because the Sun/Earth support states are out of range too. Without it the
internal Sun lookup surfaces an EphemerisRangeError naming "Body 0"
instead of the id the caller asked for.

These are pure unit tests over a fake reader: no LEB file, no network,
mode is monkeypatched.
"""

from __future__ import annotations

import pytest

from libephemeris.constants import CHIRON, MARS, SUN
from libephemeris.exceptions import EphemerisRangeError
from libephemeris.fast_calc import _escalate_sealed_range_miss

JD_IN = 2451545.0
JD_OUT = 2816788.0  # far outside the fake windows below


class _FakeReader:
    """Reader stub exposing only per-body coverage windows."""

    def __init__(self, coverage: dict[int, tuple[float, float]]):
        self._coverage = coverage

    def body_coverage(self, body_id: int):
        return self._coverage.get(body_id)


class _NoCoverageReader:
    """Reader stub without the body_coverage capability."""


@pytest.fixture
def sealed_mode(monkeypatch):
    from libephemeris import state

    monkeypatch.setattr(state, "get_calc_mode", lambda: "leb")


WINDOWS = {
    SUN: (2400000.5, 2500000.5),
    CHIRON: (2420000.5, 2470000.5),
}


@pytest.mark.unit
def test_non_leb_mode_returns_none(monkeypatch):
    from libephemeris import state

    monkeypatch.setattr(state, "get_calc_mode", lambda: "auto")
    reader = _FakeReader(WINDOWS)
    err = ValueError("out of range")
    assert _escalate_sealed_range_miss(reader, CHIRON, JD_OUT, err) is None


@pytest.mark.unit
def test_core_body_is_not_escalated(sealed_mode):
    # Core ids keep the canonical sealed-mode contract downstream.
    reader = _FakeReader(WINDOWS)
    err = ValueError("out of range")
    assert _escalate_sealed_range_miss(reader, MARS, JD_OUT, err) is None


@pytest.mark.unit
def test_reader_without_coverage_returns_none(sealed_mode):
    err = ValueError("out of range")
    assert _escalate_sealed_range_miss(_NoCoverageReader(), CHIRON, JD_OUT, err) is None


@pytest.mark.unit
def test_covered_body_failure_is_not_relabelled(sealed_mode):
    # The body IS in range: the failure is something else (corrupt section,
    # unsupported flag) and must propagate unchanged.
    reader = _FakeReader(WINDOWS)
    err = KeyError("unsupported flag")
    assert _escalate_sealed_range_miss(reader, CHIRON, JD_IN, err) is None


@pytest.mark.unit
def test_body_miss_with_sun_support_allows_keplerian_fallthrough(sealed_mode):
    # Body outside its window but the Sun still covers the date: the
    # documented out-of-range -> declared-local-model fall-through runs.
    reader = _FakeReader({SUN: (2400000.5, 2500000.5), CHIRON: (2420000.5, 2430000.5)})
    err = ValueError("out of range")
    assert _escalate_sealed_range_miss(reader, CHIRON, 2480000.0, err) is None


@pytest.mark.unit
def test_body_and_sun_both_out_names_requested_body(sealed_mode):
    reader = _FakeReader(WINDOWS)
    err = ValueError("out of range")
    escalated = _escalate_sealed_range_miss(reader, CHIRON, JD_OUT, err)
    assert isinstance(escalated, EphemerisRangeError)
    assert escalated.body_id == CHIRON
    assert f"Body {CHIRON}" in str(escalated)
    assert "Body 0" not in str(escalated)
    # The window reported is the requested body's own stored window.
    assert escalated.start_jd == pytest.approx(WINDOWS[CHIRON][0])
    assert escalated.end_jd == pytest.approx(WINDOWS[CHIRON][1])


@pytest.mark.unit
def test_channelless_body_with_sun_out_names_requested_body(sealed_mode):
    # An analytically modelled body with no channel of its own (e.g.
    # Transpluto 48) reading Sun/Earth support from the LEB: when the Sun
    # is out of range the miss IS a range miss for the requested id.
    reader = _FakeReader({SUN: (2400000.5, 2500000.5)})
    err = KeyError(48)
    escalated = _escalate_sealed_range_miss(reader, 48, JD_OUT, err)
    assert isinstance(escalated, EphemerisRangeError)
    assert escalated.body_id == 48
    assert "support states" in str(escalated)


@pytest.mark.unit
def test_channelless_body_with_sun_in_range_returns_none(sealed_mode):
    reader = _FakeReader({SUN: (2400000.5, 2500000.5)})
    err = KeyError(48)
    assert _escalate_sealed_range_miss(reader, 48, JD_IN, err) is None
