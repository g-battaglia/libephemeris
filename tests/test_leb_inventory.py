"""Coverage and inventory contract for active LEB readers."""

from __future__ import annotations

from pathlib import Path

import pytest

import libephemeris as eph


@pytest.fixture(autouse=True)
def _reviewed_bundled_base_core():
    core = Path(eph.__file__).parent / "data" / "leb2" / "base_core.leb2"
    eph.close()
    eph.set_calc_mode("leb")
    eph.set_network_policy("sealed")
    eph.set_leb_file(str(core))
    yield core
    eph.close()
    eph.set_leb_file(None)
    eph.set_calc_mode(None)
    eph.set_network_policy(None)


def test_body_coverage_is_per_body_and_machine_readable() -> None:
    sun = eph.get_body_coverage(eph.SUN)
    assert sun is not None
    assert sun.source == "LEB"
    assert sun.precision_class == "ephemeris"
    assert sun.reviewed is True
    assert sun.contains((sun.jd_start + sun.jd_end) / 2.0)
    assert not sun.contains(sun.jd_end + 1.0)
    assert sun.to_dict()["body_id"] == eph.SUN


def test_inventory_reports_files_not_only_boolean() -> None:
    inventory = eph.get_leb_inventory()
    assert inventory["ready"] is True
    assert inventory["network_policy_effective"] == "sealed"
    assert inventory["body_count"] >= 14
    assert inventory["files"]
    assert all("bodies" in item for item in inventory["files"])
    assert all("reviewed" in item for item in inventory["files"])


def test_runtime_requirements_follow_manifest_and_canonical_groups(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setenv("LIBEPHEMERIS_DATA_DIR", str(tmp_path))
    requirements = eph.get_runtime_data_requirements("extended")

    assert [item.group for item in requirements if item.kind == "leb2"] == [
        group
        for _tier in ("base", "medium", "extended")
        for group in ("core", "asteroids", "exotics", "apogee")
    ]
    assert requirements[-1].name == "extended_apogee.leb2"
    assert len(requirements) == 12
    assert all(item.kind == "leb2" for item in requirements)
    assert all(len(item.sha256) == 64 for item in requirements)
    assert requirements[0].path == str(tmp_path / "leb" / "base_core.leb2")


def test_leb_mode_fails_closed_outside_known_body_range() -> None:
    sun = eph.get_body_coverage(eph.SUN)
    assert sun is not None
    with pytest.raises(eph.EphemerisRangeError) as caught:
        eph.calc_ut(sun.jd_end + 10.0, eph.SUN, 0)
    assert caught.value.start_jd == sun.jd_start
    assert caught.value.end_jd == sun.jd_end
    assert "does not silently substitute" in str(caught.value)


def test_auto_mode_prefers_leb_before_other_sources(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A usable LEB state wins before any JPL loader in auto mode."""
    from libephemeris import state
    from libephemeris.tracing import get_trace_results, start_tracing

    def forbidden_jpl():
        raise AssertionError("auto mode consulted JPL before a usable LEB state")

    state.set_calc_mode("auto")
    monkeypatch.setattr(state, "get_planets", forbidden_jpl)
    token = start_tracing()
    try:
        position, _flags = eph.calc_ut(eph.julday(2024, 6, 15, 12.0), eph.SUN, 0)
        traces = get_trace_results()
    finally:
        token.var.reset(token)

    assert position[2] > 0.0
    assert traces[eph.SUN] == "LEB"
