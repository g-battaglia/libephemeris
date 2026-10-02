# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Synthetic regressions for review findings, without reference outputs."""

from __future__ import annotations

import threading
from contextlib import nullcontext
from types import SimpleNamespace
from uuid import uuid4

import pytest

import libephemeris as ephe
from libephemeris import context, inventory, lunar, routing, state
from libephemeris.db import __main__ as provisioner
from libephemeris.db import backend
from libephemeris.db.contract import Series
from libephemeris.exceptions import ConfigurationError, EphemerisRangeError
from tests.test_routing.test_routing import runtime as _runtime_fixture


@pytest.fixture(name="runtime")
def _coefficient_runtime(monkeypatch):
    yield from _runtime_fixture.__wrapped__(monkeypatch)


@pytest.mark.parametrize("mode", ["leb", "db", "routed"])
def test_context_jpl_denial_names_actual_mode(monkeypatch, mode):
    monkeypatch.setattr(state, "_CALC_MODE", mode)
    with pytest.raises(RuntimeError, match=repr(mode)):
        context.EphemerisContext().get_planets()


@pytest.mark.parametrize(
    "mode,label", [("leb", "LEB"), ("db", "DB"), ("routed", "routed")]
)
def test_core_range_denial_names_actual_source(monkeypatch, mode, label):
    from libephemeris.planets import _raise_leb_range_miss

    monkeypatch.setattr(state, "_CALC_MODE", mode)
    monkeypatch.setattr(inventory, "get_body_coverage", lambda *args: None)
    with pytest.raises(ephe.UnknownBodyError, match=label):
        _raise_leb_range_miss(ephe.MARS, 2451545.0)


def test_provisioner_reports_local_file_failure_and_restores_policy(
    monkeypatch, capsys
):
    monkeypatch.setattr(
        "sys.argv",
        [
            "provision",
            "--dsn",
            "postgresql://synthetic.invalid/example",
            "import",
            "--dataset",
            str(uuid4()),
            "--tier",
            "base",
            "missing.leb2",
        ],
    )
    monkeypatch.setattr(
        provisioner, "connect_provisioner", lambda _: nullcontext(object())
    )
    original = ephe.get_configured_network_policy()
    with pytest.raises(SystemExit) as error:
        provisioner.main()
    assert error.value.code == 1
    assert "artifact" in capsys.readouterr().err.lower()
    assert ephe.get_configured_network_policy() == original


def test_provisioner_does_not_print_driver_oserror_secrets(monkeypatch, capsys):
    monkeypatch.setattr(
        "sys.argv", ["provision", "--dsn", "synthetic-secret", "schema"]
    )

    def failed(_):
        raise OSError("driver included synthetic-secret")

    monkeypatch.setattr(provisioner, "connect_provisioner", failed)
    original = ephe.get_configured_network_policy()
    with pytest.raises(SystemExit):
        provisioner.main()
    assert "synthetic-secret" not in capsys.readouterr().err
    assert ephe.get_configured_network_policy() == original


class LunarWindow:
    """Native-style metadata only; no astronomical coefficients are evaluated."""

    source = "LEB"
    _manifest_verified = True

    def __init__(self, tier, bounds):
        self.tier = tier
        self.path = f"{tier}_synthetic.leb2"
        self.jd_range = bounds
        self._bodies = {
            body: Series(body, 0, 1, *bounds, bounds[1] - bounds[0], 0, 3)
            for body in (ephe.MOON, ephe.EARTH)
        }

    def has_body(self, body):
        return body in self._bodies

    def body_reader(self, body):
        return self if self.has_body(body) else None

    def body_coverage(self, body):
        return self.jd_range if self.has_body(body) else None


@pytest.fixture
def lunar_routes(monkeypatch):
    monkeypatch.setattr(state, "_CALC_MODE", "routed")
    monkeypatch.setattr(state, "_PRECISION_TIER", "extended")
    windows = {}
    owners = []
    original = routing.RoutedReader

    def create(routes, url):
        owner = original(routes, url)
        owner._tier_readers.update(windows)
        owners.append(owner)
        return owner

    monkeypatch.setattr(routing, "RoutedReader", create)
    monkeypatch.setattr(
        routing,
        "get_tier_routes",
        lambda: ({tier: ephe.TierRoute("leb") for tier in windows}, None),
    )
    return windows, owners


def test_lunar_sampling_never_spans_a_hole(lunar_routes, monkeypatch):
    windows, owners = lunar_routes
    windows.update(
        base=LunarWindow("base", (0.0, 10.0)),
        medium=LunarWindow("medium", (20.0, 30.0)),
    )
    samples = []

    def calculate(jd):
        assert any(
            reader.jd_range[0] <= jd <= reader.jd_range[1]
            for reader in windows.values()
        )
        samples.append(jd)
        return jd, 0.0, 1.0

    monkeypatch.setattr(lunar, "calc_true_lilith", calculate)
    result = lunar._sample_osculating_apogee_with_fallback(9.0, 12.0, 9)
    assert samples and all(0.0 <= jd <= 10.0 for jd in result[0])
    assert all(owner._closed and not owner._tier_readers for owner in owners)
    assert getattr(backend._operation_state, "reader", None) is None


def test_lunar_range_rejects_a_target_in_a_hole(lunar_routes):
    windows, owners = lunar_routes
    windows.update(
        base=LunarWindow("base", (0.0, 10.0)),
        medium=LunarWindow("medium", (20.0, 30.0)),
    )
    with pytest.raises(EphemerisRangeError):
        lunar._get_ephemeris_range((9.0, 21.0))
    assert all(owner._closed for owner in owners)


def test_lunar_range_accepts_continuous_cross_tier_coverage(lunar_routes):
    windows, owners = lunar_routes
    windows.update(
        base=LunarWindow("base", (0.0, 10.0)), medium=LunarWindow("medium", (9.0, 30.0))
    )
    low, high = lunar._get_ephemeris_range((5.0, 20.0))
    assert low <= 5.0 and 20.0 <= high
    assert all(owner._closed for owner in owners)


def test_local_lunar_range_does_not_inspect_optional_remote(monkeypatch, runtime):
    def forbidden(*args):
        pytest.fail("Local sampling inspected an optional DB source")

    monkeypatch.setattr(backend, "get_store", forbidden)
    low, high = lunar._get_ephemeris_range((2451544.0, 2451546.0))
    assert low <= 2451544.0 and high >= 2451546.0
    assert getattr(backend._operation_state, "reader", None) is None


@pytest.mark.parametrize(
    "action", ["routes", "db", "close", "close_db", "mode", "precision", "file"]
)
def test_reconfiguration_is_rejected_before_closing_active_inputs(runtime, action):
    routes, url = routing.get_tier_routes()
    with ephe.calculation_session():
        owner = backend._operation_state.reader
        ephe.calc(2451545.0, ephe.MARS, ephe.FLG_SPEED)
        local = owner._tier_readers["base"]
        actions = {
            "routes": lambda: ephe.set_tier_routes(routes, db_url=url),
            "db": lambda: ephe.set_db_config(None),
            "close": ephe.close,
            "close_db": backend.close_db,
            "mode": lambda: ephe.set_calc_mode("leb"),
            "precision": lambda: ephe.set_precision_tier("base"),
            "file": lambda: ephe.set_leb_file(None),
        }
        with pytest.raises(ConfigurationError, match="session"):
            actions[action]()
        assert owner._tier_readers["base"] is local
        assert ephe.get_calc_mode() == "routed"
        ephe.calc(2451545.0, ephe.MARS, ephe.FLG_SPEED)
    assert getattr(backend._operation_state, "reader", None) is None


def test_cross_thread_reconfiguration_does_not_close_borrowed_files(runtime):
    entered = threading.Event()
    release = threading.Event()
    errors = []

    def calculate():
        try:
            with ephe.calculation_session():
                ephe.calc(2451545.0, ephe.MARS, ephe.FLG_SPEED)
                entered.set()
                assert release.wait(5)
                ephe.calc(2451545.0, ephe.MARS, ephe.FLG_SPEED)
        except BaseException as error:
            errors.append(error)
            entered.set()

    worker = threading.Thread(target=calculate)
    worker.start()
    try:
        assert entered.wait(5)
        with pytest.raises(ConfigurationError, match="session"):
            ephe.set_tier_routes(None)
    finally:
        release.set()
        worker.join(5)
    assert not worker.is_alive() and not errors
    ephe.set_tier_routes(None)


def test_db_metadata_io_does_not_hold_the_lifecycle_lock(runtime, monkeypatch):
    from concurrent.futures import ThreadPoolExecutor
    from tests.test_db.conftest import DATASET_ID

    ephe.set_calc_mode("db")
    ephe.set_db_config("postgresql://synthetic.invalid/example", DATASET_ID)
    barrier = threading.Barrier(2)
    original = runtime.metadata

    def metadata(dataset):
        barrier.wait(5)
        return original(dataset)

    monkeypatch.setattr(runtime, "metadata", metadata)

    def calculate():
        with ephe.calculation_session():
            assert backend._operation_state.reader is not None

    with ThreadPoolExecutor(max_workers=2) as workers:
        list(workers.map(lambda _: calculate(), range(2)))
    ephe.set_db_config(None)


@pytest.mark.parametrize("error_class", [ConfigurationError, KeyboardInterrupt])
def test_failed_begin_releases_its_process_lease(runtime, monkeypatch, error_class):
    def invalid():
        raise error_class("synthetic failure")

    monkeypatch.setattr(routing, "get_tier_routes", invalid)
    with pytest.raises(error_class):
        with ephe.calculation_session():
            pass
    assert not getattr(backend._operation_state, "owner", None)
    ephe.set_tier_routes(None)


@pytest.mark.parametrize(
    "record",
    [
        (42, 1.0),
        (43, *([1.0] * 7)),
        (42, *([float("nan")] * 7)),
        (42, *(["not-numeric"] * 7)),
    ],
)
def test_malformed_star_rows_fail_with_typed_data_error(record):
    from libephemeris.db.reader import DBReader

    store = SimpleNamespace(
        metadata=lambda _: {0: Series(0, 0, 1, 0.0, 10.0, 10.0, 0, 3)},
        star=lambda *args: record,
    )
    reader = DBReader(store, "synthetic-dataset")
    try:
        with pytest.raises(ephe.DBDataError, match="star record"):
            reader.get_star(42)
    finally:
        reader.close()


def test_auxiliary_db_inputs_contribute_to_routed_provenance():
    reader = routing.RoutedReader(
        {"base": ephe.TierRoute("leb"), "medium": ephe.TierRoute("db")}, None
    )
    local = SimpleNamespace(jd_range=(0.0, 10.0), delta_t=lambda _: 1.0)
    remote = SimpleNamespace(
        source="DB",
        jd_range=(20.0, 30.0),
        delta_t=lambda _: 2.0,
        get_star=lambda _: "synthetic-star",
        close=lambda: None,
    )
    reader._tier_readers = {"base": local, "medium": remote}
    assert reader.delta_t(5.0) == 1.0
    assert reader.delta_t(25.0) == 2.0
    assert reader.source == "Mixed"
    del reader._tier_readers["base"]
    del reader.routes["base"]
    reader._sources.clear()
    assert reader.get_star(42) == "synthetic-star"
    assert reader.source == "DB"
    reader.close()
