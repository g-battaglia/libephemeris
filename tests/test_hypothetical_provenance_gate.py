"""Focused tests for the G-06 hypothetical source-record verifier."""

from __future__ import annotations

import ast
import copy
import hashlib
import importlib.util
import json
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import pytest

from libephemeris import hypothetical as hyp


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/check_hypothetical_provenance.py"
RECORDS = ROOT / "docs/methodology/hypothetical-source-records.json"


def _gate_module():
    spec = importlib.util.spec_from_file_location(
        "check_hypothetical_provenance", SCRIPT
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def gate():
    return _gate_module()


def _data() -> dict[str, object]:
    return json.loads(RECORDS.read_text(encoding="utf-8"))


def _write_mutation(
    tmp_path: Path, data: dict[str, object], name: str = "records.json"
) -> Path:
    path = tmp_path / name
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    return path


def test_approved_record_file_and_canonical_digest(gate):
    table = gate.load_source_records()
    assert len(table.records) == 19
    assert sum(len(record.tokens) for record in table.records) == 232
    assert sum(len(record.derived_rules) for record in table.records) == 27
    assert [record.body_id for record in table.records] == list(range(40, 59))
    assert len(table.canonical_bytes) == 32715
    assert (
        hashlib.sha256(table.canonical_bytes).hexdigest()
        == gate.SOURCE_CANONICAL_SHA256
    )
    assert hashlib.sha256(table.raw_bytes).hexdigest() == gate.SOURCE_FILE_SHA256


def test_records_are_frozen_slot_values(gate):
    table = gate.load_source_records()
    record = table.records[0]
    assert record.__slots__
    with pytest.raises(AttributeError):
        record.name = "mutated"


def test_source_mutation_fails_before_numeric_work(gate, tmp_path):
    data = _data()
    data["records"][0]["tokens"][0]["token"] = "Mutated"  # type: ignore[index]
    path = _write_mutation(tmp_path, data)
    with pytest.raises(gate.GateError, match="canonical|file-byte"):
        gate.load_source_records(path)


def test_quantum_mutation_fails_closed(gate, tmp_path):
    data = _data()
    data["records"][0]["tokens"][1]["quantum_num"] = "0.2"  # type: ignore[index]
    with pytest.raises(gate.GateError):
        gate.load_source_records(_write_mutation(tmp_path, data))


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda d: d["records"].pop(), "IDs 40 through 58"),  # type: ignore[union-attr]
        (
            lambda d: d["records"].append(copy.deepcopy(d["records"][0])),
            "duplicate body ID",
        ),  # type: ignore[union-attr,index]
        (
            lambda d: d["records"][0]["tokens"][0].__setitem__("field", "missing"),
            "canonical",
        ),  # type: ignore[index]
    ],
)
def test_missing_duplicate_and_renamed_records_fail(gate, tmp_path, mutation, message):
    data = _data()
    mutation(data)
    with pytest.raises(gate.GateError, match=message):
        gate.load_source_records(_write_mutation(tmp_path, data))


def test_unknown_variant_and_operand_namespace_fail(gate, tmp_path):
    data = _data()
    data["records"][0]["derived_rules"][0]["variant"] = "not-a-variant"  # type: ignore[index]
    with pytest.raises(gate.GateError, match="variant|digest"):
        gate.load_source_records(_write_mutation(tmp_path, data))

    data = _data()
    data["records"][0]["derived_rules"][0]["operands"][0]["namespace"] = "runtime"  # type: ignore[index]
    with pytest.raises(gate.GateError, match="namespace|digest"):
        gate.load_source_records(_write_mutation(tmp_path, data))


def test_wrong_arity_unit_category_and_cycle_fail(gate, tmp_path):
    data = _data()
    data["records"][56 - 40]["derived_rules"][3]["operands"].pop()  # type: ignore[index]
    with pytest.raises(gate.GateError, match="arity|digest"):
        gate.load_source_records(_write_mutation(tmp_path, data))

    data = _data()
    data["records"][48 - 40]["derived_rules"][0]["output_unit"] = "metre"  # type: ignore[index]
    with pytest.raises(gate.GateError, match="unit|output|digest"):
        gate.load_source_records(_write_mutation(tmp_path, data))

    data = _data()
    # Adams's first exact rule now depends on the final difference, creating a cycle.
    data["records"][52 - 40]["derived_rules"][0]["operands"][0] = {  # type: ignore[index]
        "namespace": "rule",
        "body_id": 52,
        "key": "runtime_mean_anomaly_deg",
    }
    with pytest.raises(gate.GateError, match="cycle|namespace|digest"):
        gate.load_source_records(_write_mutation(tmp_path, data))


def test_exact_sexagesimal_one_conversion(gate):
    table = gate.load_source_records()
    results = gate.evaluate_rules(table)
    leverrier = results[(51, "runtime_mean_anomaly_deg")]
    adams = results[(52, "runtime_mean_anomaly_deg")]
    assert float(leverrier.exact).hex() == "0x1.1041fdb97530fp+5"
    assert float(adams.exact).hex() == "0x1.7d9999999999ap+4"
    assert float(adams.exact) != (323.0 + 2.0 / 60.0) - (299.0 + 11.0 / 60.0)


def test_interval_schedule_has_three_precisions_and_common_intersection(gate):
    table = gate.load_source_records()
    results = gate.evaluate_rules(table)
    result = results[(56, "radius_au")]
    assert len(result.intervals) == 3
    assert gate._interval_certificate(result, hyp.HYPOTHETICAL_BODIES[56].a)[0]


def test_runtime_field_mutation_is_rejected_without_changing_verifier(
    gate, monkeypatch
):
    table = gate.load_source_records()
    results = gate.evaluate_rules(table)
    problems: list[str] = []
    original = hyp.HYPOTHETICAL_BODIES[40]
    monkeypatch.setitem(
        hyp.HYPOTHETICAL_BODIES, 40, replace(original, a=original.a + 1.0)
    )
    gate._check_runtime_records(problems, table, results)
    assert any("body 40.a" in problem for problem in problems)


def test_neely_printed_rate_and_a_are_independent(gate):
    problems: list[str] = []
    gate._check_neely_independence(problems)
    assert not problems


def test_ast_rejects_generic_close_and_runtime_expected_tables(gate):
    tree = ast.parse(SCRIPT.read_text(encoding="utf-8"))
    names = {node.name for node in ast.walk(tree) if isinstance(node, ast.FunctionDef)}
    assert "_close" not in names
    assert "_EXPECTED_CSV_ROWS" not in SCRIPT.read_text(encoding="utf-8")
    assert "Weston" not in SCRIPT.read_text(encoding="utf-8")
    assert "18.58415" not in SCRIPT.read_text(encoding="utf-8")
    assert "_WALDEMATH_DISTANCE_AU" not in SCRIPT.read_text(encoding="utf-8")


def test_cli_normal_and_explain():
    normal = subprocess.run(
        [sys.executable, str(SCRIPT)],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    explain = subprocess.run(
        [sys.executable, str(SCRIPT), "--explain"],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert normal.returncode == 0, normal.stdout + normal.stderr
    assert explain.returncode == 0, explain.stdout + explain.stderr
    gate_output = explain.stdout
    assert gate_output
    assert "32715 bytes" in gate_output
    assert (
        "cd8177155f5fa2ec0755c31b46f9aba755fadafa628c89094ac3c8bf283cd46e"
        in gate_output
    )
    assert "Neely Cupido" in gate_output
    assert "Vulcan period" not in gate_output
