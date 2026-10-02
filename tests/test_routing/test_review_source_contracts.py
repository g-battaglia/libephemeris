# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Synthetic regressions for credential, cache and coverage review findings."""

from __future__ import annotations

import math
from types import SimpleNamespace

import pytest
from click.testing import CliRunner

import libephemeris as ephe
from libephemeris import _config_toml, cache, eclipse, fast_calc, state, time_utils
from libephemeris.cli import cli
from libephemeris.db.contract import Series, validate_series
from libephemeris.db.kernel import segment_coordinates
from libephemeris.db.reader import DBReader
from libephemeris.exceptions import DBDataError, EphemerisRangeError
from libephemeris.inventory import get_reader_body_coverage
from libephemeris.leb_vector import _eval_reader_body
from libephemeris.mean_lunar_apse import _active_ephemeris_range
from libephemeris.routing import RoutedReader, TierRoute


@pytest.mark.parametrize(
    "dsn",
    [
        "postgresql://synthetic-user:synthetic-secret@test.invalid/example",
        "postgresql://synthetic-user:synthetic%2Dsecret@test.invalid/example",
        "postgresql://test.invalid/example?password=synthetic-secret",
        "host=test.invalid user=synthetic-user password='synthetic-secret'",
    ],
)
def test_config_never_displays_connection_credentials(monkeypatch, dsn):
    """Hide complete DSNs without needing a driver, connection or URL parser."""
    monkeypatch.setattr(_config_toml, "get_config_path", lambda: "synthetic.toml")
    monkeypatch.setattr(
        _config_toml, "get_all", lambda: {"db_url": dsn, "precision": "extended"}
    )
    result = CliRunner().invoke(cli, ["config"])
    assert result.exit_code == 0, result.output
    assert dsn not in result.output
    assert "synthetic-secret" not in result.output
    assert "synthetic%2Dsecret" not in result.output
    assert "db_url" in result.output and "<redacted>" in result.output
    assert "'extended'" in result.output


@pytest.fixture
def synthetic_besselian_states(monkeypatch):
    """Vary synthetic states while keeping epoch and flags identical."""
    mode = ["skyfield"]
    shift = [0.0]
    calls = []

    def calculate(jd, body, flags):
        calls.append((jd, body))
        if body == ephe.SUN:
            return (0.0, 0.0, 1.0, 0.0, 0.0, 0.0), flags
        return (90.0 + shift[0], 3.0, 0.0025, 0.0, 0.0, 0.0), flags

    monkeypatch.setattr(eclipse, "_BESSELIAN_CACHE", {})
    monkeypatch.setattr(eclipse, "get_calc_mode", lambda: mode[0])
    monkeypatch.setattr(eclipse, "calc_ut", calculate)
    monkeypatch.setattr(time_utils, "sidtime", lambda jd: 12.0)
    return mode, shift, calls


@pytest.mark.parametrize("mode_name", ["db", "routed"])
def test_remote_besselian_results_never_enter_global_cache(
    synthetic_besselian_states, mode_name
):
    """A second operation reads current states rather than prior DB results."""
    mode, shift, calls = synthetic_besselian_states
    mode[0] = mode_name
    first = eclipse._besselian_core(2451545.0)
    assert not eclipse._BESSELIAN_CACHE
    shift[0] = 7.0
    second = eclipse._besselian_core(2451545.0)
    assert len(calls) == 4 and second != first
    assert not eclipse._BESSELIAN_CACHE


@pytest.mark.parametrize("mode_name", ["db", "routed"])
def test_remote_besselian_calls_ignore_existing_file_cache(
    synthetic_besselian_states, mode_name
):
    """Switching source cannot return a cached file-derived state."""
    mode, shift, calls = synthetic_besselian_states
    first = eclipse._besselian_core(2451545.0)
    assert eclipse._BESSELIAN_CACHE
    mode[0] = mode_name
    shift[0] = 7.0
    second = eclipse._besselian_core(2451545.0)
    assert len(calls) == 4 and second != first


def test_file_besselian_results_keep_the_existing_cache(synthetic_besselian_states):
    """The ownership fix does not remove ordinary file-mode reuse."""
    _mode, _shift, calls = synthetic_besselian_states
    first = eclipse._besselian_core(2451545.0)
    assert eclipse._besselian_core(2451545.0) == first
    assert len(calls) == 2


@pytest.mark.parametrize("mode_name", ["db", "routed"])
def test_cached_file_geometry_cannot_hide_a_remote_failure(
    synthetic_besselian_states, monkeypatch, mode_name
):
    """An old success cannot conceal an outage or trigger a partial eclipse."""
    mode, _shift, _calls = synthetic_besselian_states
    eclipse._besselian_core(2451545.0)
    mode[0] = mode_name

    def failed(*args):
        raise ephe.DBError("Synthetic unavailable source")

    monkeypatch.setattr(eclipse, "calc_ut", failed)
    with pytest.raises(ephe.DBError, match="unavailable"):
        eclipse._besselian_core(2451545.0)


def test_clear_caches_invalidates_besselian_results(monkeypatch):
    """The general invalidation contract includes eclipse results."""
    monkeypatch.setattr(eclipse, "_BESSELIAN_CACHE", {(2451545.0, 2): (1.0,)})
    cache.clear_caches()
    assert not eclipse._BESSELIAN_CACHE


def test_changing_mode_invalidates_file_besselian_results(monkeypatch):
    """Ordinary file-mode transitions also invalidate state-dependent results."""
    monkeypatch.setattr(state, "_CALC_MODE", "skyfield")
    monkeypatch.setattr(eclipse, "_BESSELIAN_CACHE", {(2451545.0, 2): (1.0,)})
    state.set_calc_mode("leb")
    assert not eclipse._BESSELIAN_CACHE


def test_changing_leb_file_invalidates_besselian_results(monkeypatch):
    """A file replacement also invalidates geometry derived from its states."""
    monkeypatch.setattr(state, "_LEB_FILE", None)
    monkeypatch.setattr(state, "_LEB_READER", None)
    monkeypatch.setattr(eclipse, "_BESSELIAN_CACHE", {(2451545.0, 2): (1.0,)})
    state.set_leb_file(None)
    assert not eclipse._BESSELIAN_CACHE


@pytest.mark.parametrize(
    "count,end,interval",
    [(1, 100.0, 1.0), (2, 4.1, 2.0), (5, 1.0, 1.0), (2, 1.0, 1e308)],
)
def test_series_rejects_inconsistent_segment_coverage(count, end, interval):
    """Neither truncated grids nor unreachable extra segments are valid."""
    with pytest.raises(DBDataError):
        validate_series(Series(ephe.MARS, 0, count, 0.0, end, interval, 0, 3))


def test_kernel_does_not_clamp_an_invalid_grid_to_the_last_segment():
    """Overdeclared coverage is corruption, not endpoint interpolation."""
    series = Series(ephe.MARS, 0, 1, 0.0, 100.0, 1.0, 0, 3)
    with pytest.raises(DBDataError):
        segment_coordinates(series, 50.0)


def test_reader_rejects_bad_coverage_before_fetching_coefficients():
    """Reader construction validates the same portable contract."""
    series = Series(ephe.MARS, 0, 1, 0.0, 100.0, 1.0, 0, 3)
    store = SimpleNamespace(metadata=lambda dataset: {ephe.MARS: series})
    with pytest.raises(DBDataError):
        DBReader(store, "synthetic-dataset")


@pytest.mark.parametrize("end", [3.5, 4.0, math.nextafter(4.0, math.inf)])
def test_series_preserves_partial_final_segments_and_roundoff(end):
    """Native headers may end inside the last full-width fit interval."""
    series = Series(ephe.MARS, 0, 2, 0.0, end, 2.0, 0, 3)
    validate_series(series)
    assert segment_coordinates(series, end)[0] == 1


class SyntheticWindow:
    """Closed native-style body coverage with no astronomical source inputs."""

    source = "LEB"
    _manifest_verified = True

    def __init__(self, tier, bounds):
        self.tier = tier
        self.jd_range = bounds
        self.path = f"{tier}_synthetic.leb2"
        self._bodies = {
            ephe.MARS: Series(ephe.MARS, 0, 1, *bounds, bounds[1] - bounds[0], 0, 3)
        }
        self.calls = []

    def has_body(self, body):
        return body in self._bodies

    def has_nutation(self):
        return False

    def body_reader(self, body):
        return self if self.has_body(body) else None

    def body_coverage(self, body):
        return self.jd_range if self.has_body(body) else None

    def eval_body(self, body, jd):
        self.calls.append(jd)
        if not self.jd_range[0] <= jd <= self.jd_range[1]:
            raise ValueError("Synthetic date outside native body range")
        return (float(jd), 0.0, 0.0), (1.0, 0.0, 0.0)


@pytest.fixture
def gapped_reader():
    """Declare separate intervals, with a genuine unsupported middle gap."""
    reader = RoutedReader({"base": TierRoute("leb"), "medium": TierRoute("leb")}, None)
    reader._tier_readers = {
        "base": SyntheticWindow("base", (0.0, 10.0)),
        "medium": SyntheticWindow("medium", (20.0, 30.0)),
    }
    yield reader
    reader.close()
    fast_calc._reset_active_reader()


@pytest.mark.parametrize("requested_jd", [None, 15.0])
def test_body_coverage_never_claims_a_gap_is_covered(gapped_reader, requested_jd):
    """Envelope bounds do not imply availability between disjoint intervals."""
    coverage = get_reader_body_coverage(gapped_reader, ephe.MARS, requested_jd)
    assert coverage is not None
    assert (coverage.jd_start, coverage.jd_end) == (0.0, 30.0)
    assert not coverage.contains(15.0)
    assert not coverage.contains(math.nan)
    assert all(coverage.contains(jd) for jd in (0.0, 10.0, 20.0, 30.0))
    assert coverage.to_dict()["intervals"] == ((0.0, 10.0), (20.0, 30.0))


def test_tiered_file_adapter_also_types_a_gap(gapped_reader):
    """The vector classifier respects exact selection for file-only tiers too."""
    from libephemeris.leb_composite import TieredLEBReader

    for source in gapped_reader._tier_readers.values():
        source._readers = [source]
    reader = TieredLEBReader(gapped_reader._tier_readers)
    coverage = get_reader_body_coverage(reader, ephe.MARS)
    assert coverage is not None and not coverage.contains(15.0)
    with pytest.raises(EphemerisRangeError):
        _eval_reader_body(reader, ephe.MARS, 15.0)


@pytest.mark.parametrize("jd", [15.0, -1.0, 31.0, math.nan, math.inf])
@pytest.mark.parametrize("adapter", ["reader", "vector", "fast"])
def test_routed_range_misses_are_typed_without_evaluating_an_edge(
    gapped_reader, monkeypatch, jd, adapter
):
    """Every dispatch path rejects a gap before calling a native evaluator."""
    monkeypatch.setattr(time_utils, "deltat", lambda jd: 0.0)
    with pytest.raises(EphemerisRangeError) as caught:
        if adapter == "reader":
            gapped_reader.eval_body(ephe.MARS, jd)
        elif adapter == "vector":
            _eval_reader_body(gapped_reader, ephe.MARS, jd)
        else:
            gapped_reader.prepare(jd, ephe.MARS, 0)
            fast_calc.fast_calc_tt(
                gapped_reader,
                jd,
                ephe.MARS,
                ephe.FLG_BARYCTR | ephe.FLG_TRUEPOS | ephe.FLG_J2000 | ephe.FLG_NONUT,
            )
    assert caught.value.body_id == ephe.MARS
    assert all(not source.calls for source in gapped_reader._tier_readers.values())


@pytest.mark.parametrize("mode", ["db", "routed"])
@pytest.mark.parametrize("ifno", [0, 1])
def test_current_file_data_never_reports_an_inactive_jpl_kernel(
    monkeypatch, mode, ifno
):
    """A loaded kernel from an earlier backend cannot masquerade as active."""
    monkeypatch.setattr(state, "_CALC_MODE", mode)
    kernel = SimpleNamespace(
        path="synthetic_previous_kernel.bsp",
        spk=SimpleNamespace(segments=[SimpleNamespace(start_jd=100.0, end_jd=200.0)]),
    )
    monkeypatch.setattr(state, "_PLANETS", kernel)
    assert state.get_current_file_data(ifno) == ("", 0.0, 0.0, 0)


@pytest.mark.parametrize("mode", ["db", "routed"])
def test_analytic_lunar_range_bypasses_even_a_stale_file_report(monkeypatch, mode):
    """The local model checks backend policy before consulting file metadata."""
    monkeypatch.setattr(state, "_CALC_MODE", mode)
    monkeypatch.setattr(
        state,
        "get_current_file_data",
        lambda ifno: ("synthetic_previous.bsp", 100.0, 200.0, 440),
    )
    monkeypatch.setattr(
        state, "get_planets", lambda: pytest.fail("Loaded JPL for analytic range")
    )
    monkeypatch.setattr(
        state, "get_leb_reader", lambda: pytest.fail("Opened a file for analytic range")
    )
    assert _active_ephemeris_range() == (None, 0.0, 0.0)
