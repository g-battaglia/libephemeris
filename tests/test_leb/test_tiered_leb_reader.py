"""Best-by-date routing across reviewed LEB precision tiers."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from libephemeris import state
from libephemeris.inventory import get_body_coverage, get_runtime_data_requirements
from libephemeris.leb_composite import CompositeLEBReader, TieredLEBReader
from libephemeris.leb_format import COORD_ICRS_BARY


class _Reader:
    def __init__(
        self,
        path: str,
        ranges: dict[int, tuple[float, float]],
        marker: float,
    ) -> None:
        self.path = path
        self._ranges = ranges
        self._marker = marker
        self._bodies = {
            body_id: SimpleNamespace(
                coord_type=COORD_ICRS_BARY,
                jd_start=bounds[0],
                jd_end=bounds[1],
            )
            for body_id, bounds in ranges.items()
        }
        self._delta_t_jds = [1.0]
        self._stars = {1: SimpleNamespace(name="Test")}
        self.closed = False

    @property
    def jd_range(self) -> tuple[float, float]:
        return (
            min(bounds[0] for bounds in self._ranges.values()),
            max(bounds[1] for bounds in self._ranges.values()),
        )

    def has_body(self, body_id: int) -> bool:
        return body_id in self._ranges

    def body_coverage(self, body_id: int):
        return self._ranges.get(body_id)

    def eval_body(self, body_id: int, jd: float):
        start, end = self._ranges[body_id]
        if not start <= jd <= end:
            raise ValueError(f"JD {jd} outside [{start}, {end}]")
        return (self._marker, jd, float(body_id)), (0.0, 0.0, 0.0)

    def has_nutation(self) -> bool:
        return True

    def eval_nutation(self, jd: float):
        return self._marker, jd

    def delta_t(self, jd: float) -> float:
        if not self.jd_range[0] <= jd <= self.jd_range[1]:
            raise ValueError("outside")
        return self._marker

    def get_star(self, star_id: int):
        return self._stars[star_id]

    def warm(self, jd_start: float, jd_end: float) -> None:
        return None

    def cool(self) -> None:
        return None

    def close(self) -> None:
        self.closed = True


def _tiered_reader() -> TieredLEBReader:
    base = _Reader("/reviewed/base_core.leb2", {0: (100.0, 200.0)}, 1.0)
    medium = _Reader(
        "/reviewed/medium_core.leb2", {0: (50.0, 250.0), 15: (80.0, 220.0)}, 2.0
    )
    extended = _Reader(
        "/reviewed/extended_core.leb2",
        {0: (0.0, 300.0), 15: (60.0, 240.0)},
        3.0,
    )
    reader = TieredLEBReader(
        {
            "base": CompositeLEBReader([base]),
            "medium": CompositeLEBReader([medium]),
            "extended": CompositeLEBReader([extended]),
        }
    )
    reader._manifest_verified = True
    for child in reader._readers:
        child._manifest_verified = True
    return reader


@pytest.mark.parametrize(
    ("jd", "expected_tier", "expected_marker"),
    [
        (150.0, "base", 1.0),
        (225.0, "medium", 2.0),
        (275.0, "extended", 3.0),
    ],
)
def test_core_uses_most_precise_covering_tier(
    jd: float, expected_tier: str, expected_marker: float
) -> None:
    reader = _tiered_reader()

    position, _velocity = reader.eval_body(0, jd)

    assert position[0] == expected_marker
    assert reader.selected_tier(0, jd) == expected_tier


def test_selection_is_per_body_not_per_chart_date() -> None:
    reader = _tiered_reader()

    assert reader.selected_tier(0, 150.0) == "base"
    assert reader.selected_tier(15, 150.0) == "medium"
    assert reader.eval_body(15, 150.0)[0][0] == 2.0


def test_union_coverage_and_outside_range_preserve_range_miss() -> None:
    reader = _tiered_reader()

    assert reader.body_coverage(0) == (0.0, 300.0)
    assert reader.selected_tier(0, 301.0) is None
    with pytest.raises(ValueError, match="outside"):
        reader.eval_body(0, 301.0)


def test_date_aware_inventory_reports_the_selected_file(monkeypatch) -> None:
    reader = _tiered_reader()
    monkeypatch.setattr("libephemeris.state.get_leb_reader", lambda: reader)

    modern = get_body_coverage(0, 150.0)
    historical = get_body_coverage(0, 225.0)
    deep_time = get_body_coverage(0, 275.0)

    assert modern is not None and modern.data_file == "/reviewed/base_core.leb2"
    assert historical is not None
    assert historical.data_file == "/reviewed/medium_core.leb2"
    assert deep_time is not None
    assert deep_time.data_file == "/reviewed/extended_core.leb2"


def test_union_inventory_does_not_claim_one_tier_file(monkeypatch) -> None:
    reader = _tiered_reader()
    monkeypatch.setattr("libephemeris.state.get_leb_reader", lambda: reader)

    union = get_body_coverage(0)
    outside = get_body_coverage(0, 301.0)

    assert union is not None
    assert union.data_file is None
    assert union.group is None
    assert (union.jd_start, union.jd_end) == (0.0, 300.0)
    assert outside is not None
    assert outside.data_file is None
    assert outside.group is None


def test_inventory_distinguishes_persisted_model_classes(monkeypatch) -> None:
    direct = _Reader("/reviewed/medium_exotics.leb2", {16: (80.0, 220.0)}, 2.0)
    numerical = _Reader("/reviewed/extended_exotics.leb2", {16: (0.0, 300.0)}, 3.0)
    analytical = _Reader("/reviewed/extended_apogee.leb2", {10: (0.0, 300.0)}, 4.0)
    reader = TieredLEBReader(
        {
            "medium": CompositeLEBReader([direct]),
            "extended": CompositeLEBReader([numerical, analytical]),
        }
    )
    reader._manifest_verified = True
    for child in reader._readers:
        child._manifest_verified = True
    monkeypatch.setattr("libephemeris.state.get_leb_reader", lambda: reader)

    assert get_body_coverage(16, 150.0).precision_class == "ephemeris"
    assert get_body_coverage(16, 250.0).precision_class == "numerical-model"
    assert get_body_coverage(16).precision_class == "mixed"
    assert get_body_coverage(10, 150.0).precision_class == "analytical"


@pytest.mark.parametrize(
    ("tier", "expected_count"),
    [("base", 4), ("medium", 8), ("extended", 12)],
)
def test_sealed_runtime_requirements_include_every_more_precise_tier(
    tier: str, expected_count: int
) -> None:
    requirements = get_runtime_data_requirements(tier)

    assert len(requirements) == expected_count
    assert requirements[-1].name == f"{tier}_apogee.leb2"


def _patch_auto_discovery(
    monkeypatch,
    tier_cores: dict[str, str],
    fallback: str | None,
    reviewed: set[str],
) -> None:
    monkeypatch.delenv("LIBEPHEMERIS_LEB", raising=False)
    monkeypatch.setattr(state, "_LEB_FILE", None)
    monkeypatch.setattr(state, "_LEB_READER", None)
    monkeypatch.setattr("libephemeris._config_toml.get_str", lambda key: None)
    monkeypatch.setattr(state, "get_precision_tier", lambda: "extended")
    monkeypatch.setattr(
        state, "_discover_reviewed_leb_tier_cores", lambda: dict(tier_cores)
    )
    monkeypatch.setattr(state, "_discover_leb_file", lambda: fallback)
    monkeypatch.setattr(state, "_is_reviewed_core", lambda path: path in reviewed)
    monkeypatch.setattr(state, "_release_when_unused", lambda reader: None)
    monkeypatch.setattr(state, "_maybe_warm_reader", lambda reader: None)


def test_missing_active_core_still_routes_through_lower_reviewed_tiers(
    monkeypatch,
) -> None:
    base_path = "/reviewed/base_core.leb2"
    medium_path = "/reviewed/medium_core.leb2"
    bundled_base = "/bundled/base_core.leb2"
    _patch_auto_discovery(
        monkeypatch,
        {"base": base_path, "medium": medium_path},
        fallback=bundled_base,
        reviewed={base_path, medium_path, bundled_base},
    )
    readers = {
        base_path: _Reader(base_path, {0: (100.0, 200.0)}, 1.0),
        medium_path: _Reader(medium_path, {0: (50.0, 250.0)}, 2.0),
    }
    built: list[dict[str, str]] = []

    def _from_tier_cores(cores):
        built.append(dict(cores))
        return TieredLEBReader(
            {tier: CompositeLEBReader([readers[path]]) for tier, path in cores.items()}
        )

    monkeypatch.setattr(TieredLEBReader, "from_tier_cores", _from_tier_cores)

    reader = state._get_leb_reader_locked("leb")

    assert isinstance(reader, TieredLEBReader)
    assert built == [{"base": base_path, "medium": medium_path}]
    # JD 225 lies outside the base core (100-200); only medium covers it.
    assert reader.selected_tier(0, 225.0) == "medium"
    assert reader.eval_body(0, 225.0)[0][0] == 2.0
    assert reader._manifest_verified is True


def test_no_reviewed_cores_preserves_existing_fallback(monkeypatch) -> None:
    _patch_auto_discovery(monkeypatch, {}, fallback=None, reviewed=set())

    assert state._get_leb_reader_locked("auto") is None

    monkeypatch.setattr(state, "_LEB_READER", None)
    with pytest.raises(RuntimeError, match="no .leb file configured"):
        state._get_leb_reader_locked("leb")


def test_legacy_standard_path_keeps_single_reader_semantics(
    monkeypatch,
) -> None:
    base_path = "/reviewed/base_core.leb2"
    medium_path = "/reviewed/medium_core.leb2"
    legacy = "/home/user/.libephemeris/leb/extended_core.leb"
    _patch_auto_discovery(
        monkeypatch,
        {"base": base_path, "medium": medium_path},
        fallback=legacy,
        reviewed={base_path, medium_path},
    )
    sentinel = _Reader(legacy, {0: (0.0, 300.0)}, 9.0)
    opened: list[str] = []
    monkeypatch.setattr(
        CompositeLEBReader,
        "from_file_with_companions",
        lambda path, **kwargs: opened.append(path) or sentinel,
    )
    monkeypatch.setattr(
        TieredLEBReader,
        "from_tier_cores",
        lambda cores: pytest.fail("legacy discovery must not build tiered routing"),
    )

    assert state._get_leb_reader_locked("auto") is sentinel
    assert opened == [legacy]
