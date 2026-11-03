"""Coefficient inventory reflects actual files and immutable scalar records."""
from data_paths import CORE_LEB, PUBLISHED_CORE_LEB, is_published_core

from pathlib import Path
import inspect
import libephemeris as eph
import pytest


@pytest.fixture(autouse=True)
def configured():
    eph.set_calc_mode("leb")
    eph.set_precision_tier("base")
    eph.set_leb_file(str(CORE_LEB))


def test_coverage_record_and_file_inventory():
    record = eph.coverage(0, 2451545.0)
    assert isinstance(record, eph.BodyCoverage)
    assert record.source == "LEB" and record.contains(2451545.0)
    assert record.contains(record.jd_start) and record.contains(record.jd_end)
    assert not record.contains(float("nan"))
    assert record.to_dict()["body_id"] == 0
    with pytest.raises(AttributeError):
        record.jd_start = 0.0
    assert eph.get_body_coverage(40) is None
    info = eph.inspect_leb_file(Path(record.data_file))
    assert info["body_count"] == len(info["bodies"]) > 0
    active = eph.get_leb_inventory()
    assert active["ready"] and active["body_count"] > 0
    expected_reviewed = is_published_core(CORE_LEB)
    assert record.reviewed is expected_reviewed
    selected = next(row for row in active["files"]
                    if Path(row["path"]).resolve() == CORE_LEB.resolve())
    assert selected["reviewed"] is expected_reviewed
    assert list(inspect.signature(eph.BodyCoverage).parameters) == ["body_id", "source", "precision_class", "jd_start", "jd_end", "data_file", "group", "reviewed"]


def test_published_manifest_identity_requires_the_actual_immutable_bytes(tmp_path):
    """Positive and negative controls are independent of the owner's local file."""
    path = tmp_path / 'base_core.leb2'
    path.write_bytes(PUBLISHED_CORE_LEB.read_bytes())
    assert is_published_core(path)
    eph.set_leb_file(str(path))
    assert eph.coverage(eph.SUN, 2451545.0).reviewed

    replacement = tmp_path / 'replacement.leb2'
    replacement.write_bytes(path.read_bytes() + b'\0')
    replacement.replace(path)
    # A live image retains its own digest; only explicit reload acquires the
    # replacement. A same-named but different file must not acquire trust.
    assert eph.coverage(eph.SUN, 2451545.0).reviewed
    eph.reload_configuration()
    assert not is_published_core(path)
    assert not eph.coverage(eph.SUN, 2451545.0).reviewed


def test_cumulative_requirements_and_plain_constructor():
    counts = [len(eph.get_runtime_data_requirements(tier)) for tier in ("base", "medium", "extended")]
    assert counts == [4, 8, 12]
    record = eph.RuntimeDataRequirement("x", "leb2", None, "/tmp/x", "ab")
    assert record.to_dict() == dict(name="x", kind="leb2", group=None, path="/tmp/x", sha256="ab")
    with pytest.raises(ValueError):
        eph.get_runtime_data_requirements("unknown")
