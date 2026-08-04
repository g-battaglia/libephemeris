"""Regression tests for mean_lunar_apse._active_ephemeris_range.

The mean lunar node/apogee are analytic (IERS 2003 Delaunay arguments):
no LEB file's interval may narrow their validity — in auto mode the
router falls back to Skyfield/JPL for dates outside every LEB interval.
A provenance fix made state.get_current_file_data() report the concrete
LEB artifact (with its own interval) for body 0, and the range gate here
started enforcing the LEB window against requests the pipeline could
serve (caught by the extended-tier full suite before the v3.0.0 tag).
A JPL kernel's range, by contrast, keeps the public direct-evaluator
contract and is still enforced.
"""

from __future__ import annotations

import pytest

from libephemeris import mean_lunar_apse
from libephemeris.exceptions import EphemerisRangeError

DEEP_PAST_JD = 260098.74  # ~-13000, far outside the medium tier
MEDIUM_RANGE = (2287185.5, 2688952.5)


@pytest.mark.unit
def test_leb_interval_does_not_gate_the_analytic_series(monkeypatch):
    monkeypatch.setattr(
        mean_lunar_apse,
        "_active_ephemeris_range",
        lambda: ("/data/leb/medium_core.leb2", 0.0, 0.0),
    )
    lon, _lat, _dist = mean_lunar_apse.mean_lunar_node_position(DEEP_PAST_JD)
    assert 0.0 <= lon < 360.0


@pytest.mark.unit
def test_leb_path_with_range_falls_through_unenforced(monkeypatch):
    from libephemeris import state

    calls = {"n": 0}

    def fake_file_data(body_id):
        calls["n"] += 1
        return ("/data/leb/medium_core.leb2", *MEDIUM_RANGE, 440)

    monkeypatch.setattr(state, "get_current_file_data", fake_file_data)
    monkeypatch.setattr(state, "get_calc_mode", lambda: "auto")

    class _Reader:
        path = "/data/leb/medium_core.leb2"

    monkeypatch.setattr(state, "get_leb_reader", lambda: _Reader())
    path, start, end = mean_lunar_apse._active_ephemeris_range()
    assert (start, end) == (0.0, 0.0)  # no enforcement window
    assert calls["n"] >= 1
    lon, _lat, _dist = mean_lunar_apse.mean_lunar_node_position(DEEP_PAST_JD)
    assert 0.0 <= lon < 360.0


@pytest.mark.unit
def test_jpl_kernel_range_is_still_enforced(monkeypatch):
    from libephemeris import state

    monkeypatch.setattr(
        state,
        "get_current_file_data",
        lambda body_id: ("/data/de440.bsp", *MEDIUM_RANGE, 440),
    )
    path, start, end = mean_lunar_apse._active_ephemeris_range()
    assert (start, end) == MEDIUM_RANGE
    with pytest.raises(EphemerisRangeError):
        mean_lunar_apse.mean_lunar_node_position(DEEP_PAST_JD)
    # In-range requests still compute against the enforced kernel window.
    lon, _lat, _dist = mean_lunar_apse.mean_lunar_node_position(2451545.0)
    assert 0.0 <= lon < 360.0
