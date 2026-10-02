# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Explicit routing configuration precedence and compatibility controls."""

from __future__ import annotations

from uuid import uuid4

import pytest

import libephemeris as ephe
from libephemeris import _config_toml as config, routing, state
from libephemeris.net import require_network
from libephemeris.db.policy import require_database
from libephemeris.exceptions import ConfigurationError, NetworkSealedError


@pytest.fixture
def isolated(monkeypatch):
    previous = routing._override
    old_mode = ephe.get_calc_mode()
    old_precision = ephe.get_precision_tier()
    old_policy = ephe.get_configured_network_policy()
    routing.set_tier_routes(None)
    for name in (
        "LIBEPHEMERIS_DB_URL",
        "LIBEPHEMERIS_DB_DATASET_BASE",
        "LIBEPHEMERIS_DB_DATASET_MEDIUM",
        "LIBEPHEMERIS_DB_DATASET_EXTENDED",
    ):
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setattr(config, "_CONFIG", {})
    monkeypatch.setattr(config, "_CONFIG_LOADED", True)
    ephe.set_precision_tier("extended")
    try:
        yield
    finally:
        routing.close_routing()
        routing._override = previous
        ephe.set_calc_mode(old_mode)
        ephe.set_precision_tier(old_precision)
        ephe.set_network_policy(old_policy)


def test_toml_environment_setter_precedence_and_reset(isolated, monkeypatch, tmp_path):
    first, second, third = (str(uuid4()) for _ in range(3))
    path = tmp_path / "routes.toml"
    path.write_text(f'''[libephemeris]
mode = "routed"
db_url = "postgresql://toml.invalid/db"
[libephemeris.tier_routes.medium]
backend = "db"
dataset_id = "{first}"
''')
    assert config.load_config(path)
    assert routing.get_tier_routes()[0]["medium"].dataset_id == first
    monkeypatch.setenv("LIBEPHEMERIS_DB_DATASET_MEDIUM", second)
    monkeypatch.setenv("LIBEPHEMERIS_DB_URL", "postgresql://env.invalid/db")
    routes, url = routing.get_tier_routes()
    assert routes["medium"].dataset_id == second and "env.invalid" in url
    ephe.set_tier_routes(
        {"medium": ephe.TierRoute("db", dataset_id=third)},
        db_url="postgresql://setter.invalid/db",
    )
    routes, url = routing.get_tier_routes()
    assert routes["medium"].dataset_id == third and "setter.invalid" in url
    ephe.set_tier_routes(None)
    assert routing.get_tier_routes()[0]["medium"].dataset_id == second


@pytest.mark.parametrize(
    "value",
    [
        '"bad"',
        "{ unknown = {} }",
        '{ base = { backend = "unknown" } }',
        '{ medium = { backend = "db", files = false, dataset_id = "11111111-1111-4111-8111-111111111111" } }',
    ],
)
def test_invalid_toml_route_is_not_silently_discarded(isolated, tmp_path, value):
    path = tmp_path / "bad.toml"
    path.write_text(f"[libephemeris]\ntier_routes = {value}\n")
    with pytest.raises(ConfigurationError):
        config.load_config(path)


def test_duplicate_dataset_is_rejected_before_mutating_configuration(isolated):
    dataset = str(uuid4())
    ephe.set_tier_routes(
        {"medium": ephe.TierRoute("db", dataset_id=dataset)},
        db_url="postgresql://test.invalid/db",
    )
    with pytest.raises(ConfigurationError, match="own dataset"):
        ephe.set_tier_routes(
            {
                tier: ephe.TierRoute("db", dataset_id=dataset)
                for tier in ("medium", "extended")
            }
        )
    assert list(routing.get_tier_routes()[0]) == ["medium"]


def test_invalid_routed_dsn_type_is_rejected_in_toml(isolated, tmp_path):
    path = tmp_path / "bad-url.toml"
    path.write_text(
        '[libephemeris]\ndb_url = false\n[libephemeris.tier_routes.base]\nbackend = "leb"\nfiles = ["base_custom.leb"]\n'
    )
    with pytest.raises(ConfigurationError, match="db_url"):
        config.load_config(path)


@pytest.mark.parametrize(
    "manifest",
    [None, {}, {"artifacts": {}}, {"artifacts": [{"name": [], "sha256": "bad"}]}],
)
def test_arbitrary_manifests_do_not_claim_reviewed_science(manifest):
    assert not routing.manifest_reviewed(manifest)


def test_mode_is_not_changed_by_setting_routes_and_auto_never_uses_db(
    isolated, tmp_path, monkeypatch
):
    ephe.set_calc_mode("auto")
    ephe.set_tier_routes(
        {"medium": ephe.TierRoute("db", dataset_id=str(uuid4()))},
        db_url="postgresql://invalid/db",
    )
    assert ephe.get_calc_mode() == "auto"
    monkeypatch.setattr(
        routing, "RoutedReader", lambda *args: pytest.fail("auto selected routing")
    )
    state._get_coefficient_reader()


def test_policy_separates_db_from_http_and_preserves_explicit_sealed(isolated):
    ephe.set_calc_mode("routed")
    ephe.set_network_policy("auto")
    require_database()
    assert ephe.get_network_policy() == "sealed"
    with pytest.raises(NetworkSealedError):
        require_network("Unexpected HTTP")
    ephe.set_network_policy("sealed")
    with pytest.raises(NetworkSealedError):
        require_database()


def test_precision_is_a_ceiling_not_a_per_request_mutation(isolated):
    ephe.set_tier_routes(
        {
            "base": ephe.TierRoute("leb", ("base_custom.leb",)),
            "extended": ephe.TierRoute("db", dataset_id=str(uuid4())),
        },
        db_url="postgresql://invalid/db",
    )
    ephe.set_precision_tier("base")
    routes, _url = routing.get_tier_routes()
    assert list(routes) == ["base"]
    assert ephe.get_precision_tier() == "base"
