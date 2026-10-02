# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Giacomo Battaglia
"""Focused tests for tier-source configuration and factory validation."""

from __future__ import annotations

from collections.abc import Sequence

import pytest

import libephemeris.state as state
from libephemeris import CoefficientSourceError, Error
from libephemeris.leb_format import BodyEntry, COORD_ICRS_BARY
from libephemeris.segment_source import SegmentSource


class _Source(SegmentSource):
    """Minimal in-memory source used to exercise reader construction."""

    def fetch_body_segment(self, body_id: int, idx: int) -> Sequence[float]:
        return (float(body_id), 0.0, 0.0)

    def fetch_nutation_segment(self, idx: int) -> Sequence[float]:
        raise CoefficientSourceError("no nutation")


def _source(tier: str, group: str, body_id: int) -> _Source:
    """Create one valid source for a canonical group."""
    return _Source(
        artifact_name=f"{tier}_{group}.leb2",
        locator="memory://test",
        jd_range=(1.0, 3.0),
        bodies={
            body_id: BodyEntry(
                body_id,
                COORD_ICRS_BARY,
                1,
                1.0,
                3.0,
                2.0,
                0,
                3,
                0,
            )
        },
        reviewed=group == "core",
    )


def _all_sources(tier: str = "medium") -> list[_Source]:
    """Return all four groups in deliberately non-canonical order."""
    return [
        _source(tier, "exotics", 3),
        _source(tier, "core", 0),
        _source(tier, "apogee", 2),
        _source(tier, "asteroids", 1),
    ]


def test_tier_source_setter_precedes_environment_and_toml(monkeypatch) -> None:
    """Setter values win, while None restores env/TOML resolution."""
    monkeypatch.setenv("LIBEPHEMERIS_TIER_SOURCE_MEDIUM", "env.module:factory")
    monkeypatch.setattr(state, "_TIER_SOURCES", {"medium": "setter.module:factory"})
    monkeypatch.setattr(
        "libephemeris._config_toml.get_str",
        lambda key: "toml.module:factory" if key == "tier_source_medium" else None,
    )
    assert state.get_tier_source("medium") == "setter.module:factory"
    state._TIER_SOURCES.pop("medium")
    assert state.get_tier_source("medium") == "env.module:factory"
    monkeypatch.delenv("LIBEPHEMERIS_TIER_SOURCE_MEDIUM")
    assert state.get_tier_source("medium") == "toml.module:factory"


def test_tier_source_factory_is_ordered_core_first(monkeypatch) -> None:
    """Factories provide every group and the composite receives canonical order."""
    sources = _all_sources()
    monkeypatch.setattr(state, "_discover_reviewed_leb_tier_cores", lambda: {})
    reader = state._build_tier_source_reader(
        {"base": None, "medium": lambda _: sources}
    )
    try:
        assert [item.artifact_name for item in reader._readers] == [
            "medium_core.leb2",
            "medium_apogee.leb2",
            "medium_asteroids.leb2",
            "medium_exotics.leb2",
        ]
        assert reader._manifest_verified is False
    finally:
        reader.close()


def test_tier_source_factory_errors_are_wrapped_and_negative_cached(
    monkeypatch,
) -> None:
    """A failing configured provider is retried only after the cooldown."""
    calls = 0

    def factory(_tier: str):
        nonlocal calls
        calls += 1
        raise RuntimeError("provider unavailable")

    monkeypatch.setattr(state, "_TIER_SOURCE_FAILURE", None)
    specs = {"medium": factory}
    with pytest.raises(CoefficientSourceError, match="Could not initialize"):
        state._build_tier_source_reader(specs)
    with pytest.raises(CoefficientSourceError, match="Could not initialize"):
        state._build_tier_source_reader(specs)
    assert calls == 1
    monkeypatch.setattr(state, "_TIER_SOURCE_FAILURE", None)


def test_invalid_factory_result_is_a_source_error(monkeypatch) -> None:
    """Missing canonical groups fail closed instead of selecting a fallback."""
    monkeypatch.setattr(state, "_TIER_SOURCE_FAILURE", None)
    with pytest.raises(CoefficientSourceError, match="every canonical group"):
        state._build_tier_source_reader(
            {"medium": lambda _: [_source("medium", "core", 0)]}
        )
    monkeypatch.setattr(state, "_TIER_SOURCE_FAILURE", None)


def test_public_source_error_is_plain_exception() -> None:
    """The source error must bypass broad calculation error handlers."""
    assert CoefficientSourceError.__bases__ == (Exception,)
    assert not isinstance(CoefficientSourceError("x"), Error)


def test_unknown_tier_is_rejected() -> None:
    """Configuration accepts only the three declared precision tiers."""
    with pytest.raises(ValueError):
        state.set_tier_source("unknown", None)
    with pytest.raises(ValueError):
        state.get_tier_source("unknown")


@pytest.mark.parametrize("operation", ["calc_ut", "solcross_ut", "eclipse"])
def test_configured_fetch_failure_propagates(monkeypatch, operation) -> None:
    """A configured provider failure is not converted into backend fallback."""
    import libephemeris as ephe
    from libephemeris.leb_composite import CompositeLEBReader

    class FailingSource(_Source):
        def fetch_body_segment(self, body_id: int, idx: int) -> Sequence[float]:
            raise CoefficientSourceError("configured fetch failed")

    source = FailingSource(
        artifact_name="medium_core.leb2",
        locator="memory://failure",
        jd_range=(2450000.0, 2460000.0),
        bodies={
            body_id: BodyEntry(
                body_id, COORD_ICRS_BARY, 1000, 2450000.0, 2460000.0, 10.0, 0, 3, 0
            )
            for body_id in (0, 14)
        },
    )
    reader = CompositeLEBReader([source])
    monkeypatch.setattr(state, "get_leb_reader", lambda: reader)
    monkeypatch.setattr(state, "get_calc_mode", lambda: "leb")
    monkeypatch.setattr(
        "libephemeris.planets.get_calc_mode", lambda: "leb", raising=False
    )
    calls = {
        "calc_ut": lambda: ephe.calc_ut(2451545.0, 0, 0),
        "solcross_ut": lambda: ephe.solcross_ut(90.0, 2451545.0),
        "eclipse": lambda: ephe.sol_eclipse_when_glob(2451545.0),
    }
    try:
        with pytest.raises(CoefficientSourceError, match="configured fetch failed"):
            calls[operation]()
    finally:
        reader.close()
        from libephemeris import fast_calc

        fast_calc._reset_active_reader()
