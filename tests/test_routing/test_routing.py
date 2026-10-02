# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Native-artifact routing, lazy transport and scoped ownership contracts."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor

import pytest

import libephemeris as ephe
from libephemeris import routing, state, fast_calc
from libephemeris.db import backend
from libephemeris.exceptions import ConfigurationError, DBDataError, DBError
from tests.test_db.conftest import ArtifactStore, BUNDLED_FILE, DATASET_ID
from libephemeris.leb2_reader import LEB2Reader


class VersionStore(ArtifactStore):
    """Several immutable logical versions backed by native test coefficients."""

    def dataset_info(self, dataset_id):
        return "medium", {"artifacts": []}


class WindowReader:
    """Restrict test-local coverage without altering native arithmetic."""

    def __init__(self, reader, start=2451544.0, end=2451546.0):
        self.reader = reader
        self.jd_range = (start, end)
        self._bodies = reader._bodies
        self._readers = reader._readers

    def __getattr__(self, name):
        return getattr(self.reader, name)

    def body_coverage(self, body):
        return self.jd_range if self.reader.has_body(body) else None

    def eval_body(self, body, jd):
        if not self.jd_range[0] <= jd <= self.jd_range[1]:
            raise ValueError(f"JD {jd} outside range {self.jd_range} for body {body}")
        return self.reader.eval_body(body, jd)

    def body_reader(self, body):
        return self if self.reader.has_body(body) else None


@pytest.fixture
def runtime(monkeypatch):
    old_mode = ephe.get_calc_mode()
    old_tier = ephe.get_precision_tier()
    old_policy = ephe.get_configured_network_policy()
    ephe.close()
    ephe.set_precision_tier("extended")
    ephe.set_tier_routes(
        {
            "base": ephe.TierRoute("leb", (str(BUNDLED_FILE),)),
            "medium": ephe.TierRoute("db", dataset_id=DATASET_ID),
        },
        db_url="postgresql://test.invalid/ephemeris",
    )
    ephe.set_calc_mode("routed")
    ephe.set_network_policy("auto")
    with LEB2Reader(str(BUNDLED_FILE)) as reader:
        store = VersionStore(reader)
        monkeypatch.setattr(backend, "PostgresStore", lambda dsn: store)
        try:
            yield store
        finally:
            ephe.close()
            ephe.set_tier_routes(None)
            ephe.set_calc_mode(old_mode)
            ephe.set_precision_tier(old_tier)
            ephe.set_network_policy(old_policy)


def restrict_local():
    routes, _url = routing.get_tier_routes()
    local = routing._local_reader("base", routes["base"])
    routing._local_readers[("base", routes["base"].files)] = WindowReader(local)


@pytest.mark.parametrize("body", list(range(10)) + [ephe.EARTH, ephe.TRUE_NODE])
@pytest.mark.parametrize(
    "flags",
    [
        ephe.FLG_SPEED,
        ephe.FLG_SPEED | ephe.FLG_EQUATORIAL,
        ephe.FLG_J2000 | ephe.FLG_TRUEPOS,
    ],
)
def test_base_is_exact_and_never_contacts_remote(runtime, monkeypatch, body, flags):
    def forbidden(*args):
        pytest.fail("Local request contacted DB")

    monkeypatch.setattr(backend, "get_store", forbidden)
    with ephe.calculation_session():
        actual = ephe.calc_ut(2451545.0, body, flags)
    ephe.set_calc_mode("leb")
    ephe.set_leb_file(str(BUNDLED_FILE))
    expected = ephe.calc_ut(2451545.0, body, flags)
    assert actual == expected
    assert runtime.metadata_calls == runtime.segment_calls == 0
    assert getattr(backend._operation_state, "reader", None) is None


def test_date_routing_and_session_reuse(runtime):
    restrict_local()
    with ephe.calculation_session():
        owner = backend._operation_state.reader
        first = ephe.calc(2451500.0, ephe.MARS, ephe.FLG_SPEED)
        assert owner.source == "DB"
        count = runtime.segment_calls
        assert ephe.calc(2451500.0, ephe.MARS, ephe.FLG_SPEED) == first
        assert runtime.segment_calls == count
        with ephe.calculation_session():
            ephe.calc(2451545.0, ephe.MARS, ephe.FLG_SPEED)
            assert backend._operation_state.reader is owner
            assert owner.source == "LEB"
        assert runtime.metadata_calls == 1
        owned = owner._tier_readers["medium"]
    assert owned._closed and not owned._segments and not owned._series
    ephe.calc(2451500.0, ephe.MARS, ephe.FLG_SPEED)
    assert runtime.metadata_calls == 2


def test_nested_preparation_retains_outer_mixed_provenance(runtime):
    restrict_local()
    with ephe.calculation_session():
        owner = backend._operation_state.reader
        owner.prepare(2451500.0, ephe.MARS, ephe.FLG_SPEED)
        owner.eval_body(ephe.MARS, 2451500.0)
        assert owner.source == "DB"
        with ephe.calculation_session():
            owner.prepare(2451545.0, ephe.MARS, ephe.FLG_SPEED)
            owner.eval_body(ephe.MARS, 2451545.0)
            assert owner.source == "LEB"
        assert owner.source == "Mixed"


def test_fatal_caught_error_cannot_become_success(runtime, monkeypatch):
    restrict_local()

    def failed(*args):
        raise DBError("Disposable transport failure")

    monkeypatch.setattr(runtime, "fetch_segments", failed)
    with pytest.raises(DBError, match="Disposable"):
        with ephe.calculation_session():
            try:
                ephe.calc(2451500.0, ephe.MARS)
            except DBError:
                pass
            # A local subsequent result may not hide the previous fatal error.
            with pytest.raises(DBError):
                ephe.calc(2451545.0, ephe.SUN)
    assert backend._operation_state.reader is None


def test_tier_identity_is_verified(runtime, monkeypatch):
    restrict_local()
    monkeypatch.setattr(runtime, "dataset_info", lambda dataset: ("extended", {}))
    with pytest.raises(DBDataError, match="tier"):
        ephe.calc(2451500.0, ephe.MARS)
    assert backend._operation_state.reader is None


def test_missing_declared_file_is_fatal(runtime, tmp_path):
    ephe.set_tier_routes(
        {
            "base": ephe.TierRoute("leb", (str(tmp_path / "base_missing.leb2"),)),
            "medium": ephe.TierRoute("db", dataset_id=DATASET_ID),
        },
        db_url="postgresql://test.invalid/db",
    )
    with pytest.raises(ephe.RoutingDataError):
        ephe.calc(2451545.0, ephe.MARS)
    assert runtime.metadata_calls == 0


def test_local_readiness_does_not_probe_remote(runtime, monkeypatch):
    monkeypatch.setattr(
        backend, "get_store", lambda *args: pytest.fail("Remote readiness probe")
    )
    report = ephe.get_runtime_inventory("base")
    assert report["ready"] and len(report["sources"]) == 1
    assert ephe.get_leb_inventory()["ready"]


def test_remote_inventory_reports_no_fake_file_or_url(runtime):
    report = ephe.get_runtime_inventory("medium")
    assert report["ready"]
    assert report["sources"][1]["source"] == "DB"
    assert report["sources"][1]["dataset_id"] == DATASET_ID
    assert "test.invalid" not in str(report)
    assert all(item["data_file"] is None for item in report["sources"][1]["bodies"])


def test_date_less_mixed_coverage_never_claims_one_tier_or_dataset(runtime):
    coverage = ephe.get_body_coverage(ephe.MARS)
    assert coverage is not None
    assert coverage.source == "Mixed"
    assert coverage.tier is None and coverage.dataset_id is None
    assert coverage.data_file is None


def test_context_local_reader_cannot_override_routes(runtime):
    restrict_local()
    ctx = ephe.EphemerisContext()
    ctx.set_leb_file(str(BUNDLED_FILE))
    assert ctx.calc(2451500.0, ephe.MARS, ephe.FLG_SPEED) == ephe.calc(
        2451500.0, ephe.MARS, ephe.FLG_SPEED
    )
    assert runtime.metadata_calls == 2


def test_thread_owners_are_independent(runtime):
    restrict_local()

    def request(index):
        with ephe.calculation_session():
            owner = backend._operation_state.reader
            result = ephe.calc(2451500.0 + index, ephe.MARS)
        assert backend._operation_state.reader is None
        return result, owner

    with ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(request, range(12)))
    assert len({id(owner) for _result, owner in results}) == 12
    assert runtime.metadata_calls == 12


def test_speed_and_light_time_cross_boundary_without_extrapolation(runtime):
    restrict_local()
    with ephe.calculation_session():
        reader = backend._operation_state.reader
        value = ephe.calc(2451544.001, ephe.PLUTO, ephe.FLG_SPEED)
        assert reader.source == "Mixed"
        assert runtime.metadata_calls == 1
    ephe.set_calc_mode("leb")
    ephe.set_leb_file(str(BUNDLED_FILE))
    assert value == ephe.calc(2451544.001, ephe.PLUTO, ephe.FLG_SPEED)


def test_vector_construction_is_lazy_for_unrequested_remote_bodies(
    runtime, monkeypatch
):
    from libephemeris.leb_vector import get_leb_vector_ephemeris

    monkeypatch.setattr(
        backend, "get_store", lambda *args: pytest.fail("Vector eager DB access")
    )
    with ephe.calculation_session():
        reader = backend._operation_state.reader
        vector = get_leb_vector_ephemeris(reader)
        vector["earth"].at(state.get_timescale().tt_jd(2451545.0))
    assert runtime.metadata_calls == 0


def test_input_limit_is_aggregated_across_remote_tiers(runtime, monkeypatch):
    extended_id = "01234567-89ab-cdef-0123-456789abcdee"
    monkeypatch.setattr(
        runtime,
        "dataset_info",
        lambda dataset: ("medium" if dataset == DATASET_ID else "extended", {}),
    )
    ephe.set_tier_routes(
        {
            "medium": ephe.TierRoute("db", dataset_id=DATASET_ID),
            "extended": ephe.TierRoute("db", dataset_id=extended_id),
        },
        db_url="postgresql://test.invalid/db",
    )
    with ephe.calculation_session():
        owner = backend._operation_state.reader
        readers = [owner._reader(tier) for tier in ("medium", "extended")]
        for index in range(400):
            reader = readers[index % 2]
            reader.eval_body(ephe.SUN, 2451545.0 + index * 40)
            assert sum(len(r._segments) for r in readers) <= 256
    assert all(not r._segments and r._closed for r in readers)


def test_db_inventory_uses_actual_publication_tier(runtime):
    ephe.set_db_config("postgresql://test.invalid/db", DATASET_ID)
    ephe.set_calc_mode("db")
    try:
        inventory = ephe.get_runtime_inventory("base")
        assert inventory["ready"]
        assert inventory["sources"][0]["tier"] == "medium"
        assert not ephe.get_runtime_inventory("extended")["ready"]
    finally:
        ephe.set_db_config(None)
        ephe.set_calc_mode("routed")


def test_global_caches_do_not_retain_routed_owner(runtime):
    from libephemeris import leb_vector

    ephe.set_tier_routes(
        {"medium": ephe.TierRoute("db", dataset_id=DATASET_ID)},
        db_url="postgresql://test.invalid/db",
    )
    with ephe.calculation_session():
        reader = backend._operation_state.reader
        ephe.calc(2451500.0, ephe.MARS)
        leb_vector.get_leb_vector_ephemeris(reader)
        assert leb_vector._CACHED_READER is not reader
        assert all(key[0] != id(reader) for key in fast_calc._leb_frame_cache)
        owned_db = reader._tier_readers["medium"]
        assert owned_db._frame_cache
    assert not owned_db._frame_cache
    # File-only frame entries may persist; no entry contains a DB/routed owner.
    assert all(key[0] != id(reader) for key in fast_calc._leb_frame_cache)


def test_interrupted_scope_cleans_up_without_masking_interrupt(runtime):
    restrict_local()
    with pytest.raises(KeyboardInterrupt):
        with ephe.calculation_session():
            ephe.calc(2451500.0, ephe.MARS)
            raise KeyboardInterrupt()
    assert backend._operation_state.reader is None


def test_ordinary_caught_input_error_does_not_poison_session(runtime):
    with ephe.calculation_session():
        with pytest.raises(ConfigurationError):
            ephe.calc(2451545.0, ephe.MARS, ephe.FLG_TOPOCTR)
        assert ephe.calc(2451545.0, ephe.SUN)


def test_selector_closed_endpoints_and_missing_metadata():
    tiers = ("base", "medium", "extended")
    assert routing.select_tier(tiers, {}, 12).needs_metadata
    assert routing.select_tier(tiers, {"base": (10, 14)}, 14).tier == "base"
    assert routing.select_tier(tiers, {"base": (10, 14)}, 15) == routing.RouteDecision(
        "medium", True
    )
    assert routing.select_tier(tiers, dict.fromkeys(tiers), 15).tier is None


@pytest.mark.parametrize(
    "route",
    [
        {},
        {"invalid": ephe.TierRoute("db", dataset_id=DATASET_ID)},
        {"base": ephe.TierRoute("db", dataset_id="bad")},
        {"base": ephe.TierRoute("leb", ("medium_core.leb2",))},
        {"base": ephe.TierRoute("leb", ("a", "a"))},
        {"base": {"backend": "db", "dataset_id": DATASET_ID, "unknown": True}},
    ],
)
def test_invalid_config_fails_atomically(route):
    previous = routing._override
    with pytest.raises(ConfigurationError):
        ephe.set_tier_routes(route)
    assert routing._override is previous


@pytest.mark.parametrize("entrypoint", ["calc", "calc_ut"])
def test_typed_minor_body_range_miss_retains_approved_local_model(
    runtime, monkeypatch, entrypoint
):
    """Typing a stored miss must not disable the existing curated model."""
    from libephemeris.db.contract import Series

    metadata = runtime.metadata

    def with_narrow_ceres(dataset):
        records = metadata(dataset)
        records[ephe.CERES] = Series(ephe.CERES, 0, 1, 2451500.0, 2451510.0, 10.0, 0, 3)
        return records

    monkeypatch.setattr(runtime, "metadata", with_narrow_ceres)
    calculate = getattr(ephe, entrypoint)
    actual = calculate(2451545.0, ephe.CERES, ephe.FLG_SPEED)
    ephe.set_calc_mode("leb")
    ephe.set_leb_file(str(BUNDLED_FILE))
    assert actual == calculate(2451545.0, ephe.CERES, ephe.FLG_SPEED)
