"""Exhaustive negative graph-contract tests for G-06 source records."""

from __future__ import annotations

import copy
import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/check_hypothetical_provenance.py"
RECORDS = ROOT / "docs/methodology/hypothetical-source-records.json"


def _module():
    spec = importlib.util.spec_from_file_location("g06_exhaustive", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def gate():
    return _module()


def _data() -> dict[str, object]:
    return json.loads(RECORDS.read_text(encoding="utf-8"))


def _record(data: dict[str, object], body_id: int) -> dict[str, object]:
    return next(record for record in data["records"] if record["body_id"] == body_id)  # type: ignore[union-attr]


def _rule(data: dict[str, object], body_id: int, key: str) -> dict[str, object]:
    return next(
        rule
        for rule in _record(data, body_id)["derived_rules"]
        if rule["rule_key"] == key
    )  # type: ignore[index]


def _mutations(data: dict[str, object]):
    """Yield every incompatible rule and resolved operand metadata mutation."""
    incompatible = {
        "category": "unsupported",
        "unit": "text",
        "exactness": "structural",
        "method": "structural_exact",
        "namespace": "runtime",
        "body_id": 999,
        "key": "not-a-reference",
        "quantum": "0.1",
    }
    genuine = 0
    for record in data["records"]:  # type: ignore[union-attr]
        for rule in record["derived_rules"]:
            body_id = record["body_id"]
            key = rule["rule_key"]
            variant = rule["variant"]
            for field, value in (
                ("output_category", incompatible["category"]),
                ("output_unit", incompatible["unit"]),
                ("output_exactness", incompatible["exactness"]),
                ("interval_method", incompatible["method"]),
                ("output_quantum", incompatible["quantum"]),
                (
                    "variant",
                    "rate_cancellation"
                    if variant != "rate_cancellation"
                    else "difference",
                ),
            ):
                if rule[field] == value:
                    continue
                mutated = copy.deepcopy(data)
                _rule(mutated, body_id, key)[field] = value
                yield f"rule:{body_id}:{key}:{field}", mutated
            for index, operand in enumerate(rule["operands"]):
                for field, value in (
                    ("namespace", incompatible["namespace"]),
                    ("body_id", incompatible["body_id"]),
                    ("key", incompatible["key"]),
                ):
                    mutated = copy.deepcopy(data)
                    _rule(mutated, body_id, key)["operands"][index][field] = value
                    yield f"operand:{body_id}:{key}:{index}:{field}", mutated
                if operand["namespace"] == "token":
                    token = next(
                        token
                        for token in _record(data, body_id)["tokens"]
                        if token["field"] == operand["key"]
                    )
                    for field, value in (
                        ("unit", "text"),
                        ("category", "unsupported"),
                        ("quantum_num", "0.1"),
                    ):
                        if token[field] == value:
                            continue
                        mutated = copy.deepcopy(data)
                        token2 = next(
                            token
                            for token in _record(mutated, body_id)["tokens"]
                            if token["field"] == operand["key"]
                        )
                        token2[field] = value
                        yield f"token:{body_id}:{key}:{index}:{field}", mutated
                else:
                    child = _rule(data, body_id, operand["key"])
                    for field, value in (
                        ("output_unit", "text"),
                        ("output_category", "unsupported"),
                        ("output_exactness", "structural"),
                        ("output_quantum", "0.1"),
                        ("interval_method", "structural_exact"),
                    ):
                        if child[field] == value:
                            continue
                        mutated = copy.deepcopy(data)
                        child2 = _rule(mutated, body_id, operand["key"])
                        child2[field] = value
                        yield f"child:{body_id}:{key}:{index}:{field}", mutated


def test_every_rule_and_operand_mutation_fails_before_eval(gate, monkeypatch):
    mutations = list(_mutations(_data()))
    assert len(mutations) == 522
    accepted_invalid: list[str] = []
    evaluated: list[object] = []

    def forbidden(*args, **kwargs):
        evaluated.append((args, kwargs))
        raise AssertionError("numeric evaluation reached mutated graph")

    monkeypatch.setattr(gate, "_eval_rule", forbidden)
    for label, data in mutations:
        path = ROOT / "tests" / ".tmp-g06-mutated.json"
        path.write_text(json.dumps(data) + "\n", encoding="utf-8")
        try:
            with pytest.raises(gate.GateError):
                table = gate.load_source_records(path)
                gate.evaluate_rules(table)
        except AssertionError:
            accepted_invalid.append(label)
        finally:
            path.unlink(missing_ok=True)
    assert not evaluated
    assert accepted_invalid == []
    assert len(mutations) == 522


def test_structural_child_interval_rejected_before_eval(gate, monkeypatch, tmp_path):
    data = _data()
    child = _rule(data, 52, "source_mean_longitude_deg")
    child["output_exactness"] = "interval"
    child["interval_method"] = "arb_outward"
    path = tmp_path / "child-interval.json"
    path.write_text(json.dumps(data) + "\n", encoding="utf-8")
    monkeypatch.setattr(gate, "_eval_rule", lambda *a, **k: pytest.fail("eval reached"))
    with pytest.raises(gate.GateError):
        gate.load_source_records(path)


def test_no_zip_truncation_in_operand_contract(gate, tmp_path):
    data = _data()
    _rule(data, 56, "rate_deg_per_day")["operands"].pop()
    path = tmp_path / "truncated.json"
    path.write_text(json.dumps(data) + "\n", encoding="utf-8")
    with pytest.raises(gate.GateError, match="arity|quantum"):
        gate.load_source_records(path)


def test_normative_metadata_digest_is_independent_and_pinned(gate, monkeypatch):
    expected = gate._metadata_canonical_bytes()
    assert len(expected) == 2797
    assert (
        gate.VARIANT_METADATA_SHA256
        == "6931031d63c06f16086214f2cc9ddff339b50a0bdcafe0f1114edd5ef4ed0dc1"
    )
    assert gate._metadata_canonical_bytes() == expected
    for variant in tuple(gate._VARIANT_METADATA):
        original = gate._VARIANT_METADATA[variant]
        mutated = (
            ((("token", "text", {"unsupported"}, {"rounded"}),), original[1])
            if not original[0]
            else ((), original[1])
        )
        monkeypatch.setitem(gate._VARIANT_METADATA, variant, mutated)
        with pytest.raises(gate.GateError, match="metadata (byte count|digest)"):
            gate.load_source_records()
        monkeypatch.setitem(gate._VARIANT_METADATA, variant, original)


def test_interval_certificate_requires_exact_precision_keys(gate):
    table = gate.load_source_records()
    result = gate.evaluate_rules(table)[(56, "radius_au")]
    for keys in (
        (160, 256),
        (160, 512),
        (256,),
        (),
        (160, 160, 512),
        (256, 160, 512),
        (True, 256, 512),
        ("160", 256, 512),
    ):
        broken = gate._RuleResult(result.rule, result.exact, result.intervals, keys)
        with pytest.raises(gate.GateError, match="exactly one pass"):
            gate._interval_certificate(broken)
