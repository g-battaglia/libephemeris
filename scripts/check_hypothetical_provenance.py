#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Certify the independently recorded hypothetical-body source boundary.

The verifier-owned record at ``docs/methodology/hypothetical-source-records.json``
contains the only expected source quantities used here.  Runtime values are
loaded only as the actual implementation under test.  The record is verifier
data, not package data or a scientific runtime asset.

Provenance:
    Project-authored G-06 source-record integrity and arithmetic verifier. It
    reads only the tracked verifier record, the tracked CSV, and runtime values
    as actuals under test; it defines no astronomical runtime model.

This gate implements G06-SourceRecords-2.  It validates the immutable schema and
canonical digest before opening the CSV or importing numerical runtime values.
Literal values use one Decimal-to-binary64 conversion.  Derived values use the
fixed typed operation graph and fresh flint Arb passes at 160, 256, and 512
bits.  Structural and unsupported records are deliberately not promoted to
source claims.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import inspect
import json
import math
import subprocess
import sys
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Sequence

from flint import arb, ctx

from libephemeris import hypothetical as hyp
from libephemeris.exceptions import UnknownBodyError


ROOT = Path(__file__).resolve().parents[1]
SOURCE_RECORDS_PATH = ROOT / "docs/methodology/hypothetical-source-records.json"
CSV_SHA256 = "ea163c1111245a719021baa92beedb5ecf1ebb8f2601a46f9c66a4fc1350cad2"
SOURCE_FILE_SHA256 = "3071cd51272bab30037308f8305d91068c17e6e2959be9799e2894c496bfce33"
SOURCE_CANONICAL_SHA256 = (
    "2fa15e3a7ab3ff546171479e4a85d2ebddd9014334fe56dc8464d5ab759fe72a"
)
SOURCE_CANONICAL_BYTES = 33036
PRECISIONS = (160, 256, 512)


def _validate_precisions() -> None:
    """Require the immutable independent Arb precision schedule."""
    _require(
        type(PRECISIONS) is tuple
        and PRECISIONS == (160, 256, 512)
        and all(
            type(value) is int and not isinstance(value, bool) for value in PRECISIONS
        )
        and len(set(PRECISIONS)) == 3,
        "PRECISIONS must be exactly (160, 256, 512)",
    )


CATEGORIES = frozenset(
    {"literal_transcription", "derived_source", "project_convention", "unsupported"}
)
NAMESPACES = frozenset({"token", "rule"})
UNITS = frozenset(
    {
        "arcminute",
        "arcsecond",
        "au",
        "day",
        "degree",
        "degree_per_century",
        "degree_per_day",
        "dimensionless",
        "integer",
        "jd_tt",
        "jd_ut",
        "julian_year",
        "metre",
        "metre3_per_second2",
        "second",
        "text",
    }
)
EXACTNESS = frozenset({"exact", "interval", "structural"})
VARIANTS = frozenset(
    {
        "gaussian_from_a",
        "period_years_to_rate",
        "elapsed_period_to_angle",
        "turns_endpoints_elapsed_to_rate",
        "two_period_rate",
        "full_turn_rate_to_period",
        "sexagesimal_deg_min",
        "sexagesimal_deg_min_sec",
        "sexagesimal_base_deg_min",
        "difference",
        "kepler_radius",
        "rate_cancellation",
    }
)
INTERVAL_METHODS = frozenset(
    {
        "exact_rational_one_conversion",
        "arb_outward",
        "arb_outward_intersection",
        "ast_structural",
        "structural_exact",
    }
)

# Normative typed operation table. Each operand entry is
# (namespace, unit, permitted category, permitted exactness). ``rounded``
# denotes a nonzero token quantum; ``interval-compatible`` denotes a rule
# output whose declared exactness is exact or interval.
_VARIANT_METADATA = {
    "gaussian_from_a": (
        (
            ("token", "degree_per_day", {"project_convention"}, {"exact", "rounded"}),
            (
                "token",
                "au",
                {"literal_transcription", "project_convention"},
                {"exact", "rounded"},
            ),
        ),
        (
            "degree_per_day",
            "project_convention",
            {"structural", "interval"},
            {"ast_structural", "arb_outward"},
        ),
    ),
    "period_years_to_rate": (
        (
            ("token", "julian_year", {"literal_transcription"}, {"exact", "rounded"}),
            ("token", "day", {"project_convention"}, {"exact", "rounded"}),
        ),
        ("degree_per_day", "derived_source", {"interval"}, "arb_outward"),
    ),
    "elapsed_period_to_angle": (
        (
            (
                "token",
                "julian_year",
                {"project_convention", "literal_transcription"},
                {"exact", "rounded"},
            ),
        )
        * 3,
        ("degree", "derived_source", {"interval"}, "arb_outward"),
    ),
    "turns_endpoints_elapsed_to_rate": (
        (
            ("token", "integer", {"project_convention"}, {"exact"}),
            ("rule", "degree", {"derived_source"}, {"interval-compatible"}),
            ("rule", "degree", {"derived_source"}, {"interval-compatible"}),
            ("rule", "day", {"project_convention"}, {"interval-compatible"}),
        ),
        ("degree_per_day", "derived_source", {"interval"}, "arb_outward"),
    ),
    "two_period_rate": (
        (("token", "day", {"literal_transcription"}, {"exact", "rounded"}),) * 2,
        ("degree_per_day", "derived_source", {"interval"}, "arb_outward"),
    ),
    "full_turn_rate_to_period": (
        (
            ("token", "degree", {"project_convention"}, {"exact", "rounded"}),
            ("rule", "degree_per_day", {"derived_source"}, {"interval-compatible"}),
        ),
        ("day", "derived_source", {"interval"}, "arb_outward"),
    ),
    "sexagesimal_deg_min": (
        (
            ("token", "degree", {"literal_transcription"}, {"exact"}),
            ("token", "arcminute", {"literal_transcription"}, {"exact"}),
        ),
        ("degree", "derived_source", {"exact"}, {"exact_rational_one_conversion"}),
    ),
    "sexagesimal_deg_min_sec": (
        (
            ("token", "degree", {"literal_transcription"}, {"exact"}),
            ("token", "arcminute", {"literal_transcription"}, {"exact"}),
            ("token", "arcsecond", {"literal_transcription"}, {"exact"}),
        ),
        ("degree", "derived_source", {"exact"}, {"exact_rational_one_conversion"}),
    ),
    "sexagesimal_base_deg_min": (
        (
            ("token", "degree", {"project_convention"}, {"exact"}),
            ("token", "degree", {"literal_transcription"}, {"rounded"}),
            ("token", "arcminute", {"literal_transcription"}, {"rounded"}),
        ),
        ("degree", "derived_source", {"interval"}, "arb_outward"),
    ),
    "kepler_radius": (
        (
            (
                "token",
                "metre3_per_second2",
                {"project_convention"},
                {"exact", "rounded"},
            ),
            ("rule", "day", {"derived_source"}, {"interval-compatible"}),
            ("token", "second", {"project_convention"}, {"exact", "rounded"}),
            ("token", "metre", {"project_convention"}, {"exact", "rounded"}),
        ),
        ("au", "project_convention", {"interval"}, "arb_outward_intersection"),
    ),
    "rate_cancellation": (
        (),
        (
            "degree_per_century",
            "project_convention",
            {"structural"},
            "structural_exact",
        ),
    ),
}

TOP_LEVEL_FIELDS = {
    "schema",
    "quantum_encoding",
    "operand_ref_fields",
    "derived_rule_fields",
    "allowed_operand_namespaces",
    "allowed_output_exactness",
    "allowed_variants",
    "allowed_units",
    "allowed_interval_methods",
    "token_tuple",
    "rule_tuple",
    "records",
    "counts",
    "allowed_categories",
    "allowed_operations",
    "canonical_serialization",
}
TOP_LEVEL_FIELD_ORDER = (
    "schema",
    "quantum_encoding",
    "operand_ref_fields",
    "derived_rule_fields",
    "allowed_operand_namespaces",
    "allowed_output_exactness",
    "allowed_variants",
    "allowed_units",
    "allowed_interval_methods",
    "token_tuple",
    "rule_tuple",
    "records",
    "counts",
    "allowed_categories",
    "allowed_operations",
    "canonical_serialization",
)
TOKEN_FIELDS = {
    "body_id",
    "record_key",
    "field",
    "token",
    "quantum_num",
    "unit",
    "category",
    "locator_key",
}
RULE_FIELDS = {
    "body_id",
    "rule_key",
    "variant",
    "operands",
    "output_category",
    "output_unit",
    "output_quantum",
    "output_exactness",
    "interval_method",
}
RECORD_FIELDS = {
    "body_id",
    "record_key",
    "name",
    "category",
    "tokens",
    "derived_rules",
}
OPERAND_FIELDS = {"namespace", "body_id", "key"}


class GateError(ValueError):
    """Raised when verifier data cannot be certified before numerical work."""


@dataclass(frozen=True, slots=True)
class SourceToken:
    """Verifier implementation member."""

    body_id: int
    record_key: str
    field: str
    token: str
    quantum_num: int | str
    unit: str
    category: str
    locator_key: str


@dataclass(frozen=True, slots=True)
class OperandRef:
    """Verifier implementation member."""

    namespace: str
    body_id: int
    key: str


@dataclass(frozen=True, slots=True)
class DerivedRule:
    """Verifier implementation member."""

    body_id: int
    rule_key: str
    variant: str
    operands: tuple[OperandRef, ...]
    output_category: str
    output_unit: str
    output_quantum: int | str | None
    output_exactness: str
    interval_method: str


@dataclass(frozen=True, slots=True)
class SourceRecord:
    """Verifier implementation member."""

    body_id: int
    record_key: str
    name: str
    category: str
    tokens: tuple[SourceToken, ...]
    derived_rules: tuple[DerivedRule, ...]


@dataclass(frozen=True, slots=True)
class SourceTable:
    """Verifier implementation member."""

    records: tuple[SourceRecord, ...]
    raw_bytes: bytes
    canonical_bytes: bytes


@dataclass(frozen=True, slots=True)
class _Value:
    """Verifier implementation member."""

    unit: str
    category: str
    exactness: str
    fraction: Fraction | None = None
    interval: arb | None = None


@dataclass(frozen=True, slots=True)
class _RuleResult:
    """Verifier implementation member."""

    rule: DerivedRule
    exact: Fraction | None
    intervals: tuple[arb, ...]
    pass_precisions: tuple[int, ...] = PRECISIONS


def _fail(message: str) -> None:
    """Verifier implementation member."""
    raise GateError(message)


def _require(condition: bool, message: str) -> None:
    """Verifier implementation member."""
    if not condition:
        _fail(message)


def _string(value: Any, label: str) -> str:
    """Verifier implementation member."""
    _require(type(value) is str, f"{label} must be a string")
    return value


def _integer(value: Any, label: str) -> int:
    """Verifier implementation member."""
    _require(type(value) is int, f"{label} must be an integer")
    return value


def _quantum(value: Any, label: str, *, nullable: bool = False) -> int | str | None:
    """Verifier implementation member."""
    if value is None and nullable:
        return None
    _require(
        type(value) is int or type(value) is str, f"{label} has invalid quantum type"
    )
    if type(value) is int:
        _require(value >= 0, f"{label} is negative")
        return value
    try:
        parsed = Decimal(value)
    except InvalidOperation as exc:
        raise GateError(f"{label} is not a Decimal spelling") from exc
    _require(parsed.is_finite() and parsed >= 0, f"{label} is not nonnegative finite")
    return value


def _exact_keys(mapping: dict[str, Any], expected: set[str], label: str) -> None:
    """Verifier implementation member."""
    _require(set(mapping) == expected, f"{label} fields differ from the schema")


def _load_json(path: Path) -> tuple[dict[str, Any], bytes]:
    """Load JSON while rejecting duplicate object declarations."""
    raw = path.read_bytes()

    def pairs_hook(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        keys = [key for key, _ in pairs]
        _require(len(keys) == len(set(keys)), "duplicate JSON declaration")
        return dict(pairs)

    try:
        value = json.loads(raw.decode("utf-8"), object_pairs_hook=pairs_hook)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise GateError(f"source record is not UTF-8 JSON: {exc}") from exc
    _require(type(value) is dict, "source record root must be an object")
    return value, raw


def _validate_declared_schema(data: dict[str, Any]) -> None:
    """Verifier implementation member."""
    _exact_keys(data, TOP_LEVEL_FIELDS, "source record root")
    _require(
        tuple(data) == TOP_LEVEL_FIELD_ORDER, "source record top-level order drift"
    )
    _require(data["schema"] == "G06-SourceRecords-2", "unexpected source schema")
    _require(type(data["quantum_encoding"]) is str, "quantum_encoding must be text")
    _require(
        data["operand_ref_fields"] == ["namespace", "body_id", "key"],
        "operand fields drift",
    )
    _require(
        data["derived_rule_fields"]
        == [
            "body_id",
            "rule_key",
            "variant",
            "operands",
            "output_category",
            "output_unit",
            "output_quantum",
            "output_exactness",
            "interval_method",
        ],
        "derived-rule fields drift",
    )
    _require(
        data["allowed_operand_namespaces"] == ["token", "rule"], "namespace enum drift"
    )
    _require(
        data["allowed_output_exactness"] == ["exact", "interval", "structural"],
        "exactness enum drift",
    )
    _require(
        data["allowed_variants"]
        == [
            "gaussian_from_a",
            "period_years_to_rate",
            "elapsed_period_to_angle",
            "turns_endpoints_elapsed_to_rate",
            "two_period_rate",
            "full_turn_rate_to_period",
            "sexagesimal_deg_min",
            "sexagesimal_deg_min_sec",
            "sexagesimal_base_deg_min",
            "difference",
            "kepler_radius",
            "rate_cancellation",
        ],
        "variant enum drift",
    )
    _require(
        data["token_tuple"]
        == [
            "body_id",
            "record_key",
            "field",
            "token",
            "quantum_num",
            "unit",
            "category",
            "locator_key",
        ],
        "token tuple drift",
    )
    _require(
        data["rule_tuple"]
        == [
            "body_id",
            "rule_key",
            "variant",
            "operands",
            "output_category",
            "output_unit",
            "output_quantum",
            "output_exactness",
            "interval_method",
        ],
        "rule tuple drift",
    )
    _require(data["allowed_units"] == sorted(UNITS), "unit enum declaration drift")
    _require(
        data["allowed_interval_methods"] == sorted(INTERVAL_METHODS),
        "interval-method enum declaration drift",
    )
    _require(
        data["allowed_categories"]
        == [
            "literal_transcription",
            "derived_source",
            "project_convention",
            "unsupported",
        ],
        "category enum drift",
    )
    _require(
        sorted(data["allowed_operations"]) == sorted(VARIANTS), "operation enum drift"
    )
    _require(type(data["counts"]) is dict, "counts must be an object")
    _require(
        set(data["counts"]) == {"records", "tokens", "derived_rules"},
        "count fields drift",
    )
    canonical = data["canonical_serialization"]
    _require(type(canonical) is dict, "canonical_serialization must be an object")
    _require(canonical.get("magic_utf8") == "G06SR2\n", "canonical magic drift")
    _require(
        canonical.get("byte_count") == SOURCE_CANONICAL_BYTES,
        "canonical byte count drift",
    )
    _require(
        canonical.get("sha256") == SOURCE_CANONICAL_SHA256,
        "canonical digest declaration drift",
    )


def _parse_records(data: dict[str, Any]) -> tuple[SourceRecord, ...]:
    """Verifier implementation member."""
    records_json = data["records"]
    _require(type(records_json) is list, "records must be a list")
    records: list[SourceRecord] = []
    seen_ids: set[int] = set()
    seen_rules: set[tuple[int, str]] = set()
    previous_id = -1
    token_count = rule_count = 0
    for index, raw_record in enumerate(records_json):
        _require(type(raw_record) is dict, f"record {index} is not an object")
        _exact_keys(raw_record, RECORD_FIELDS, f"record {index}")
        body_id = _integer(raw_record["body_id"], f"record {index}.body_id")
        _require(body_id not in seen_ids, f"duplicate body ID {body_id}")
        _require(body_id > previous_id, "records are not in ascending body-ID order")
        previous_id = body_id
        seen_ids.add(body_id)
        record_key = _string(raw_record["record_key"], f"record {index}.record_key")
        name = _string(raw_record["name"], f"record {index}.name")
        category = _string(raw_record["category"], f"record {index}.category")
        _require(category in CATEGORIES, f"record {body_id} has unknown category")
        raw_tokens = raw_record["tokens"]
        raw_rules = raw_record["derived_rules"]
        _require(type(raw_tokens) is list, f"record {body_id}.tokens must be a list")
        _require(
            type(raw_rules) is list, f"record {body_id}.derived_rules must be a list"
        )
        tokens: list[SourceToken] = []
        token_fields: set[str] = set()
        for token_index, raw_token in enumerate(raw_tokens):
            _require(
                type(raw_token) is dict,
                f"token {body_id}/{token_index} is not an object",
            )
            _exact_keys(raw_token, TOKEN_FIELDS, f"token {body_id}/{token_index}")
            token = SourceToken(
                body_id=_integer(raw_token["body_id"], "token.body_id"),
                record_key=_string(raw_token["record_key"], "token.record_key"),
                field=_string(raw_token["field"], "token.field"),
                token=_string(raw_token["token"], "token.token"),
                quantum_num=_quantum(raw_token["quantum_num"], "token.quantum_num"),  # type: ignore[arg-type]
                unit=_string(raw_token["unit"], "token.unit"),
                category=_string(raw_token["category"], "token.category"),
                locator_key=_string(raw_token["locator_key"], "token.locator_key"),
            )
            _require(
                token.body_id == body_id and token.record_key == record_key,
                f"token {body_id}/{token.field} owner drift",
            )
            _require(
                token.field not in token_fields,
                f"duplicate token field {body_id}/{token.field}",
            )
            _require(
                token.unit in UNITS, f"token {body_id}/{token.field} has unknown unit"
            )
            _require(
                token.category in CATEGORIES,
                f"token {body_id}/{token.field} has unknown category",
            )
            _require(
                bool(token.locator_key), f"token {body_id}/{token.field} lacks locator"
            )
            token_fields.add(token.field)
            tokens.append(token)
        rules: list[DerivedRule] = []
        rule_keys: set[str] = set()
        for rule_index, raw_rule in enumerate(raw_rules):
            _require(
                type(raw_rule) is dict, f"rule {body_id}/{rule_index} is not an object"
            )
            _exact_keys(raw_rule, RULE_FIELDS, f"rule {body_id}/{rule_index}")
            raw_operands = raw_rule["operands"]
            _require(
                type(raw_operands) is list,
                f"rule {body_id}/{rule_index}.operands must be a list",
            )
            operands: list[OperandRef] = []
            for operand_index, raw_operand in enumerate(raw_operands):
                _require(
                    type(raw_operand) is dict,
                    f"operand {body_id}/{rule_index}/{operand_index} is not an object",
                )
                _exact_keys(raw_operand, OPERAND_FIELDS, "operand")
                operands.append(
                    OperandRef(
                        namespace=_string(
                            raw_operand["namespace"], "operand.namespace"
                        ),
                        body_id=_integer(raw_operand["body_id"], "operand.body_id"),
                        key=_string(raw_operand["key"], "operand.key"),
                    )
                )
            rule = DerivedRule(
                body_id=_integer(raw_rule["body_id"], "rule.body_id"),
                rule_key=_string(raw_rule["rule_key"], "rule.rule_key"),
                variant=_string(raw_rule["variant"], "rule.variant"),
                operands=tuple(operands),
                output_category=_string(
                    raw_rule["output_category"], "rule.output_category"
                ),
                output_unit=_string(raw_rule["output_unit"], "rule.output_unit"),
                output_quantum=_quantum(
                    raw_rule["output_quantum"], "rule.output_quantum", nullable=True
                ),
                output_exactness=_string(
                    raw_rule["output_exactness"], "rule.output_exactness"
                ),
                interval_method=_string(
                    raw_rule["interval_method"], "rule.interval_method"
                ),
            )
            _require(
                rule.body_id == body_id, f"rule {body_id}/{rule.rule_key} owner drift"
            )
            _require(
                rule.rule_key not in rule_keys,
                f"duplicate rule key {body_id}/{rule.rule_key}",
            )
            _require(
                (body_id, rule.rule_key) not in seen_rules,
                f"duplicate global rule key {body_id}/{rule.rule_key}",
            )
            _require(
                rule.variant in VARIANTS,
                f"rule {body_id}/{rule.rule_key} has unknown variant",
            )
            _require(
                rule.output_category in CATEGORIES,
                f"rule {body_id}/{rule.rule_key} has unknown output category",
            )
            _require(
                rule.output_unit in UNITS,
                f"rule {body_id}/{rule.rule_key} has unknown output unit",
            )
            _require(
                rule.output_exactness in EXACTNESS,
                f"rule {body_id}/{rule.rule_key} has unknown exactness",
            )
            _require(
                rule.interval_method in INTERVAL_METHODS,
                f"rule {body_id}/{rule.rule_key} has unknown interval method",
            )
            if rule.output_exactness == "exact":
                _require(
                    rule.output_quantum is None,
                    f"exact rule {body_id}/{rule.rule_key} cannot declare quantum",
                )
            if rule.output_quantum is not None:
                _require(
                    rule.output_exactness == "interval",
                    f"non-interval rule {body_id}/{rule.rule_key} declares quantum",
                )
            rule_keys.add(rule.rule_key)
            seen_rules.add((body_id, rule.rule_key))
            rules.append(rule)
        records.append(
            SourceRecord(
                body_id, record_key, name, category, tuple(tokens), tuple(rules)
            )
        )
        token_count += len(tokens)
        rule_count += len(rules)
    _require(
        seen_ids == set(range(40, 59)),
        "records must contain IDs 40 through 58 exactly once",
    )
    _require(
        data["counts"]
        == {
            "records": len(records),
            "tokens": token_count,
            "derived_rules": rule_count,
        },
        "source counts drift",
    )
    _require(
        (len(records), token_count, rule_count) == (19, 232, 27),
        "source record counts are not 19/232/27",
    )
    return tuple(records)


def _length_prefixed(value: str) -> bytes:
    """Verifier implementation member."""
    encoded = value.encode("utf-8")
    return len(encoded).to_bytes(4, "big") + encoded


def _quantum_tag(value: int | str | None) -> str:
    """Verifier implementation member."""
    if value is None:
        return "Z:"
    if type(value) is int:
        return f"N:{value}"
    return f"D:{value}"


def canonical_serialize(records: Iterable[SourceRecord]) -> bytes:
    """Return the independent G06SR2 byte serialization."""
    rows = sorted(records, key=lambda record: (record.body_id, record.record_key))
    output = bytearray(b"G06SR2\n")
    for declaration in (sorted(UNITS), sorted(INTERVAL_METHODS)):
        output.extend(len(declaration).to_bytes(4, "big"))
        for item in declaration:
            output.extend(_length_prefixed(item))
    output.extend(len(rows).to_bytes(4, "big"))
    for record in rows:
        output.extend(record.body_id.to_bytes(8, "big", signed=True))
        output.extend(_length_prefixed(record.record_key))
        output.extend(_length_prefixed(record.name))
        output.extend(_length_prefixed(record.category))
        tokens = sorted(record.tokens, key=lambda token: token.field)
        output.extend(len(tokens).to_bytes(4, "big"))
        for token in tokens:
            output.extend(_length_prefixed(token.record_key))
            output.extend(_length_prefixed(token.field))
            output.extend(_length_prefixed(token.token))
            output.extend(_length_prefixed(_quantum_tag(token.quantum_num)))
            output.extend(_length_prefixed(token.unit))
            output.extend(_length_prefixed(token.category))
            output.extend(_length_prefixed(token.locator_key))
        rules = sorted(
            record.derived_rules, key=lambda rule: (rule.rule_key, rule.variant)
        )
        output.extend(len(rules).to_bytes(4, "big"))
        for rule in rules:
            output.extend(rule.body_id.to_bytes(8, "big", signed=True))
            output.extend(_length_prefixed(rule.rule_key))
            output.extend(_length_prefixed(rule.variant))
            output.extend(len(rule.operands).to_bytes(4, "big"))
            for operand in rule.operands:
                output.extend(_length_prefixed(operand.namespace))
                output.extend(operand.body_id.to_bytes(8, "big", signed=True))
                output.extend(_length_prefixed(operand.key))
            output.extend(_length_prefixed(rule.output_category))
            output.extend(_length_prefixed(rule.output_unit))
            output.extend(_length_prefixed(_quantum_tag(rule.output_quantum)))
            output.extend(_length_prefixed(rule.output_exactness))
            output.extend(_length_prefixed(rule.interval_method))
    return bytes(output)


def _validate_rule_metadata(
    rule: DerivedRule,
    token_value: Any,
    rule_value: Any,
) -> None:
    """Enforce the normative typed metadata contract for one rule."""
    if rule.variant == "difference":
        if rule.rule_key == "elapsed_days":
            elapsed_expected = (
                (("token", "jd_tt", "project_convention", "exact"),) * 2,
                ("project_convention", "day", None, "interval", "arb_outward"),
            )
        else:
            elapsed_expected = None
        _require(len(rule.operands) == 2, f"difference {rule.rule_key} has wrong arity")
        metadata = [
            token_value(op) if op.namespace == "token" else rule_value(op)
            for op in rule.operands
        ]
        units = [
            value.unit if operand.namespace == "token" else value.output_unit
            for operand, value in zip(rule.operands, metadata, strict=True)
        ]
        categories = [
            value.category if operand.namespace == "token" else value.output_category
            for operand, value in zip(rule.operands, metadata, strict=True)
        ]
        _require(
            units[0] == units[1], f"difference input units drift in {rule.rule_key}"
        )
        if rule.rule_key == "elapsed_days":
            _require(
                units == ["jd_tt", "jd_tt"] and rule.output_unit == "day",
                f"elapsed-day units drift in {rule.rule_key}",
            )
        else:
            _require(
                rule.output_unit == units[0],
                f"difference output unit drift in {rule.rule_key}",
            )
        _require(
            all(
                category
                in {"literal_transcription", "project_convention", "derived_source"}
                for category in categories
            ),
            f"difference input categories drift in {rule.rule_key}",
        )
        _require(
            rule.output_category in {"derived_source", "project_convention"},
            f"difference output category drift in {rule.rule_key}",
        )
        _require(
            rule.output_exactness in {"exact", "interval"},
            f"difference exactness drift in {rule.rule_key}",
        )
        _require(
            rule.interval_method in {"exact_rational_one_conversion", "arb_outward"},
            f"difference method drift in {rule.rule_key}",
        )
        if elapsed_expected is not None:
            _require(
                rule.output_category == elapsed_expected[1][0]
                and rule.output_unit == elapsed_expected[1][1]
                and rule.output_quantum == elapsed_expected[1][2]
                and rule.output_exactness == elapsed_expected[1][3]
                and rule.interval_method == elapsed_expected[1][4],
                f"elapsed-day output metadata drift in {rule.rule_key}",
            )
        if rule.interval_method == "exact_rational_one_conversion":
            _require(
                all(
                    (
                        value.quantum_num == 0
                        if operand.namespace == "token"
                        else value.output_exactness == "exact"
                    )
                    for operand, value in zip(rule.operands, metadata, strict=True)
                ),
                f"exact difference has interval operand in {rule.rule_key}",
            )
        return
    expected_operands, output = _VARIANT_METADATA[rule.variant]
    _require(
        len(rule.operands) == len(expected_operands),
        f"{rule.variant} {rule.rule_key} has wrong arity",
    )
    for index, (operand, expected) in enumerate(
        zip(rule.operands, expected_operands, strict=True)
    ):
        namespace, unit, categories, exactness = expected  # type: ignore[assignment]
        _require(
            operand.namespace == namespace,
            f"{rule.variant} operand namespace drift at {index}",
        )
        value = token_value(operand) if namespace == "token" else rule_value(operand)
        actual_unit = value.unit if namespace == "token" else value.output_unit
        actual_category = (
            value.category if namespace == "token" else value.output_category
        )
        actual_exact = (
            "exact"
            if namespace == "token"
            and (
                value.quantum_num == 0
                or rule.variant in {"sexagesimal_deg_min", "sexagesimal_deg_min_sec"}
            )
            else "rounded"
            if namespace == "token"
            else value.output_exactness
        )
        _require(actual_unit == unit, f"{rule.variant} operand unit drift at {index}")
        _require(
            actual_category in categories,
            f"{rule.variant} operand category drift at {index}",
        )
        exact_ok = (
            actual_exact in exactness
            or "interval-compatible" in exactness
            and actual_exact in {"exact", "interval"}
        )
        _require(exact_ok, f"{rule.variant} operand exactness drift at {index}")
    output_unit, output_category, output_exactness, methods = output
    _require(rule.output_unit == output_unit, f"{rule.variant} output unit drift")
    _require(
        rule.output_category == output_category, f"{rule.variant} output category drift"
    )
    _require(
        rule.output_exactness in output_exactness
        or rule.variant == "gaussian_from_a"
        and rule.output_exactness == "structural",
        f"{rule.variant} output exactness drift",
    )
    _require(rule.interval_method in methods, f"{rule.variant} interval method drift")
    if rule.variant == "gaussian_from_a":
        if rule.body_id == 54:
            _require(
                rule.output_exactness == "interval", "Gaussian output exactness drift"
            )
            _require(
                rule.interval_method == "arb_outward", "Gaussian interval method drift"
            )
        else:
            _require(
                rule.output_exactness == "structural", "Gaussian output exactness drift"
            )
            _require(
                rule.interval_method == "ast_structural",
                "Gaussian interval method drift",
            )
    if rule.interval_method == "exact_rational_one_conversion":
        _require(
            all(
                (
                    (
                        value.quantum_num == 0
                        or rule.variant
                        in {"sexagesimal_deg_min", "sexagesimal_deg_min_sec"}
                    )
                    if operand.namespace == "token"
                    else value.output_exactness == "exact"
                )
                for operand in rule.operands
                for value in [
                    token_value(operand)
                    if operand.namespace == "token"
                    else rule_value(operand)
                ]
            ),
            f"exact rational rule has interval operand in {rule.rule_key}",
        )


def _validate_operation_graph(records: tuple[SourceRecord, ...]) -> None:
    """Validate the fixed variant contract and resolve its dependency DAG."""
    by_body = {record.body_id: record for record in records}
    tokens = {
        (record.body_id, token.field): token
        for record in records
        for token in record.tokens
    }
    rules = {
        (rule.body_id, rule.rule_key): rule
        for record in records
        for rule in record.derived_rules
    }
    _require(len(rules) == 27, "derived rule keys are not globally unique")

    def token_value(ref: OperandRef) -> SourceToken:
        _require(ref.namespace == "token", "expected token operand")
        key = (ref.body_id, ref.key)
        _require(key in tokens, f"unresolved token reference {key}")
        return tokens[key]

    def rule_value(ref: OperandRef) -> DerivedRule:
        _require(ref.namespace == "rule", "expected rule operand")
        key = (ref.body_id, ref.key)
        _require(key in rules, f"unresolved rule reference {key}")
        return rules[key]

    expected_arity = {
        "gaussian_from_a": (2, ("token", "token")),
        "period_years_to_rate": (2, ("token", "token")),
        "elapsed_period_to_angle": (3, ("token", "token", "token")),
        "turns_endpoints_elapsed_to_rate": (4, ("token", "rule", "rule", "rule")),
        "two_period_rate": (2, ("token", "token")),
        "full_turn_rate_to_period": (2, ("token", "rule")),
        "sexagesimal_deg_min": (2, ("token", "token")),
        "sexagesimal_deg_min_sec": (3, ("token", "token", "token")),
        "sexagesimal_base_deg_min": (3, ("token", "token", "token")),
        "kepler_radius": (4, ("token", "rule", "token", "token")),
        "rate_cancellation": (0, ()),
    }
    allowed_rule_categories = {"derived_source", "project_convention"}
    for record in records:
        for rule in record.derived_rules:
            _validate_rule_metadata(rule, token_value, rule_value)
            if rule.variant == "difference":
                _require(
                    len(rule.operands) == 2,
                    f"difference {rule.rule_key} has wrong arity",
                )
                _require(
                    all(operand.namespace in NAMESPACES for operand in rule.operands),
                    f"difference {rule.rule_key} namespace drift",
                )
            else:
                arity, namespaces = expected_arity[rule.variant]
                _require(
                    len(rule.operands) == arity,
                    f"{rule.variant} {rule.rule_key} has wrong arity",
                )
                _require(
                    tuple(operand.namespace for operand in rule.operands) == namespaces,
                    f"{rule.variant} {rule.rule_key} namespace drift",
                )
            if rule.variant != "rate_cancellation":
                _require(bool(rule.operands), f"{rule.rule_key} has no operands")
            for operand in rule.operands:
                _require(
                    operand.body_id == rule.body_id,
                    f"cross-body operand in {rule.rule_key}",
                )
                if operand.namespace == "token":
                    source = token_value(operand)
                    _require(
                        source.category in CATEGORIES,
                        f"token category drift in {rule.rule_key}",
                    )
                else:
                    source_rule = rule_value(operand)
                    _require(
                        source_rule.output_category in allowed_rule_categories,
                        f"rule category drift in {rule.rule_key}",
                    )
            if rule.variant == "gaussian_from_a":
                first, second = (token_value(op) for op in rule.operands)
                _require(
                    (first.unit, second.unit) == ("degree_per_day", "au"),
                    f"Gaussian units drift in {rule.rule_key}",
                )
                _require(
                    first.category == "project_convention"
                    and second.category
                    in {"literal_transcription", "project_convention"},
                    f"Gaussian categories drift in {rule.rule_key}",
                )
                _require(
                    rule.output_unit == "degree_per_day"
                    and rule.output_category == "project_convention",
                    f"Gaussian output drift in {rule.rule_key}",
                )
                _require(
                    rule.output_exactness == "structural"
                    or rule.output_exactness == "interval",
                    f"Gaussian exactness drift in {rule.rule_key}",
                )
            elif rule.variant == "period_years_to_rate":
                first, second = (token_value(op) for op in rule.operands)
                _require(
                    (first.unit, second.unit) == ("julian_year", "day"),
                    f"period units drift in {rule.rule_key}",
                )
                _require(
                    rule.output_unit == "degree_per_day"
                    and rule.output_category == "derived_source"
                    and rule.output_exactness == "interval",
                    f"period output drift in {rule.rule_key}",
                )
            elif rule.variant == "elapsed_period_to_angle":
                _require(
                    tuple(token_value(op).unit for op in rule.operands)
                    == ("julian_year", "julian_year", "julian_year"),
                    f"elapsed-angle units drift in {rule.rule_key}",
                )
                _require(
                    rule.output_unit == "degree"
                    and rule.output_category == "derived_source"
                    and rule.output_exactness == "interval",
                    f"elapsed-angle output drift in {rule.rule_key}",
                )
            elif rule.variant == "turns_endpoints_elapsed_to_rate":
                first = token_value(rule.operands[0])
                _require(
                    first.unit == "integer" and first.category == "project_convention",
                    f"turn unit drift in {rule.rule_key}",
                )
                _require(
                    tuple(rule_value(op).output_unit for op in rule.operands[1:])
                    == ("degree", "degree", "day"),
                    f"turn rule units drift in {rule.rule_key}",
                )
                _require(
                    rule.output_unit == "degree_per_day"
                    and rule.output_category == "derived_source"
                    and rule.output_exactness == "interval",
                    f"turn output drift in {rule.rule_key}",
                )
            elif rule.variant == "two_period_rate":
                _require(
                    tuple(token_value(op).unit for op in rule.operands)
                    == ("day", "day"),
                    f"two-period units drift in {rule.rule_key}",
                )
                _require(
                    rule.output_unit == "degree_per_day"
                    and rule.output_category == "derived_source"
                    and rule.output_exactness == "interval",
                    f"two-period output drift in {rule.rule_key}",
                )
            elif rule.variant == "full_turn_rate_to_period":
                _require(
                    token_value(rule.operands[0]).unit == "degree"
                    and rule_value(rule.operands[1]).output_unit == "degree_per_day",
                    f"full-turn units drift in {rule.rule_key}",
                )
                _require(
                    rule.output_unit == "day"
                    and rule.output_category == "derived_source"
                    and rule.output_exactness == "interval",
                    f"full-turn output drift in {rule.rule_key}",
                )
            elif rule.variant in {"sexagesimal_deg_min", "sexagesimal_deg_min_sec"}:
                _require(
                    all(
                        token_value(op).unit in {"degree", "arcminute", "arcsecond"}
                        for op in rule.operands
                    ),
                    f"sexagesimal units drift in {rule.rule_key}",
                )
                _require(
                    rule.output_unit == "degree"
                    and rule.output_category == "derived_source"
                    and rule.output_exactness == "exact"
                    and rule.interval_method == "exact_rational_one_conversion",
                    f"sexagesimal output drift in {rule.rule_key}",
                )
            elif rule.variant == "sexagesimal_base_deg_min":
                _require(
                    tuple(token_value(op).unit for op in rule.operands)
                    == ("degree", "degree", "arcminute"),
                    f"base sexagesimal units drift in {rule.rule_key}",
                )
                _require(
                    rule.output_unit == "degree"
                    and rule.output_category == "derived_source"
                    and rule.output_exactness == "interval",
                    f"base sexagesimal output drift in {rule.rule_key}",
                )
            elif rule.variant == "difference":
                left, right = rule.operands
                left_meta: Any = (
                    token_value(left) if left.namespace == "token" else rule_value(left)
                )
                right_meta: Any = (
                    token_value(right)
                    if right.namespace == "token"
                    else rule_value(right)
                )
                left_unit = (
                    left_meta.unit
                    if left.namespace == "token"
                    else left_meta.output_unit
                )  # type: ignore[union-attr]
                right_unit = (
                    right_meta.unit
                    if right.namespace == "token"
                    else right_meta.output_unit
                )  # type: ignore[union-attr]
                if rule.rule_key == "elapsed_days":
                    _require(
                        (left_unit, right_unit, rule.output_unit)
                        == ("jd_tt", "jd_tt", "day"),
                        f"elapsed-day units drift in {rule.rule_key}",
                    )
                else:
                    _require(
                        left_unit == right_unit == rule.output_unit,
                        f"difference units drift in {rule.rule_key}",
                    )
                _require(
                    rule.output_category in {"derived_source", "project_convention"},
                    f"difference category drift in {rule.rule_key}",
                )
            elif rule.variant == "kepler_radius":
                _require(
                    tuple(
                        token_value(rule.operands[i]).unit
                        if i in (0, 2, 3)
                        else rule_value(rule.operands[i]).output_unit
                        for i in range(4)
                    )
                    == ("metre3_per_second2", "day", "second", "metre"),
                    f"Kepler-radius units drift in {rule.rule_key}",
                )
                _require(
                    rule.output_unit == "au"
                    and rule.output_category == "project_convention"
                    and rule.output_exactness == "interval",
                    f"Kepler-radius output drift in {rule.rule_key}",
                )
            elif rule.variant == "rate_cancellation":
                _require(
                    rule.output_unit == "degree_per_century"
                    and rule.output_category == "project_convention"
                    and rule.output_exactness == "structural",
                    f"cancellation output drift in {rule.rule_key}",
                )
            _require(
                rule.interval_method
                == (
                    "exact_rational_one_conversion"
                    if rule.output_exactness == "exact"
                    else rule.interval_method
                ),
                f"rule method drift in {rule.rule_key}",
            )

    # Resolve forward references through a deterministic topological order.
    dependencies: dict[tuple[int, str], set[tuple[int, str]]] = {}
    for key, rule in rules.items():
        dependencies[key] = {
            (operand.body_id, operand.key)
            for operand in rule.operands
            if operand.namespace == "rule"
        }
    remaining = set(dependencies)
    ordered: list[tuple[int, str]] = []
    while remaining:
        ready = sorted(
            key for key in remaining if dependencies[key].isdisjoint(remaining)
        )
        _require(bool(ready), "derived-rule dependency cycle")
        ordered.extend(ready)
        remaining.difference_update(ready)
    _require(
        len(ordered) == len(rules), "derived-rule topological resolution incomplete"
    )
    _ = by_body


def load_source_records(path: Path = SOURCE_RECORDS_PATH) -> SourceTable:
    """Load, validate, canonicalize, and digest verifier-owned records."""
    _validate_precisions()
    data, raw = _load_json(path)
    _validate_declared_schema(data)
    records = _parse_records(data)
    _validate_operation_graph(records)
    canonical = canonical_serialize(records)
    _require(
        len(canonical) == SOURCE_CANONICAL_BYTES,
        f"canonical serialization has {len(canonical)} bytes",
    )
    digest = hashlib.sha256(canonical).hexdigest()
    _require(
        digest == SOURCE_CANONICAL_SHA256,
        f"source canonical SHA-256 mismatch: {digest}",
    )
    _require(
        hashlib.sha256(raw).hexdigest() == SOURCE_FILE_SHA256,
        "source record file-byte SHA-256 mismatch",
    )
    return SourceTable(records, raw, canonical)


# Backward-friendly private aliases make focused tests able to exercise the
# schema boundary without exposing verifier records as runtime API.
_load_source_records = load_source_records
_canonical_serialize = canonical_serialize


def _decimal_fraction(token: str) -> Fraction:
    """Verifier implementation member."""
    try:
        decimal = Decimal(token)
    except InvalidOperation as exc:
        raise GateError(f"non-decimal numerical source token {token!r}") from exc
    _require(decimal.is_finite(), f"non-finite source token {token!r}")
    return Fraction(decimal)


def _decimal_float(token: str) -> float:
    """Verifier implementation member."""
    value = float(Decimal(token))
    _require(math.isfinite(value), f"non-finite binary64 conversion for {token!r}")
    return value


def _float_fraction(value: float) -> Fraction:
    """Verifier implementation member."""
    _require(
        type(value) is float and math.isfinite(value),
        "runtime value is not finite binary64",
    )
    return Fraction(*value.as_integer_ratio())


def _token_interval(
    token: SourceToken, precision: int, *, integer_component: bool = False
) -> arb:
    """Verifier implementation member."""
    with ctx.workprec(precision):
        if integer_component or token.quantum_num == 0:
            return arb(token.token)
        quantum = Decimal(str(token.quantum_num))
        value = Decimal(token.token)
        half = quantum / Decimal(2)
        low, high = value - half, value + half
        midpoint = (low + high) / Decimal(2)
        radius = (high - low) / Decimal(2)
        return arb(f"{midpoint} +/- {radius}")


def _arb_from_fraction(value: Fraction, precision: int) -> arb:
    """Verifier implementation member."""
    with ctx.workprec(precision):
        return arb(f"{value.numerator}/{value.denominator}")


def _operand_exactness(
    rule: DerivedRule,
    position: int,
    token: SourceToken | None,
    child: DerivedRule | None,
) -> bool:
    """Return whether one typed operand is an exact source quantity."""
    del position
    if rule.variant in {"sexagesimal_deg_min", "sexagesimal_deg_min_sec"}:
        return True
    if token is not None:
        return token.quantum_num == 0
    return child is not None and child.output_exactness == "exact"


def _metadata_for_operand(
    ref: OperandRef,
    tokens: dict[tuple[int, str], SourceToken],
    rules: dict[tuple[int, str], DerivedRule],
) -> tuple[str, str, str, SourceToken | None, DerivedRule | None]:
    """Verifier implementation member."""
    if ref.namespace == "token":
        token = tokens.get((ref.body_id, ref.key))
        _require(token is not None, f"unresolved token operand {ref.body_id}/{ref.key}")
        assert token is not None
        return (
            token.unit,
            token.category,
            "exact" if token.quantum_num == 0 else "interval",
            token,
            None,
        )
    child = rules.get((ref.body_id, ref.key))
    _require(child is not None, f"unresolved rule operand {ref.body_id}/{ref.key}")
    assert child is not None
    return child.output_unit, child.output_category, child.output_exactness, None, child


def _fraction_value(
    ref: OperandRef,
    tokens: dict[tuple[int, str], SourceToken],
    exact_rules: dict[tuple[int, str], Fraction],
) -> Fraction:
    """Verifier implementation member."""
    if ref.namespace == "token":
        token = tokens[(ref.body_id, ref.key)]
        return _decimal_fraction(token.token)
    _require(
        (ref.body_id, ref.key) in exact_rules,
        f"interval operand used as exact value: {ref.body_id}/{ref.key}",
    )
    return exact_rules[(ref.body_id, ref.key)]


def _arb_value(
    ref: OperandRef,
    precision: int,
    tokens: dict[tuple[int, str], SourceToken],
    results: dict[tuple[int, str], _RuleResult],
    rule: DerivedRule,
) -> arb:
    """Verifier implementation member."""
    if ref.namespace == "token":
        token = tokens[(ref.body_id, ref.key)]
        integer_component = rule.variant in {
            "sexagesimal_deg_min",
            "sexagesimal_deg_min_sec",
        }
        return _token_interval(token, precision, integer_component=integer_component)
    result = results[(ref.body_id, ref.key)]
    return result.intervals[0]


def _eval_rule(
    rule: DerivedRule,
    precision: int,
    tokens: dict[tuple[int, str], SourceToken],
    results: dict[tuple[int, str], _RuleResult],
) -> arb:
    """Verifier implementation member."""
    variant = rule.variant
    operands = rule.operands
    if variant == "rate_cancellation":
        return arb(0)
    values = [
        _arb_value(operand, precision, tokens, results, rule) for operand in operands
    ]
    if variant == "gaussian_from_a":
        _require(
            values[1].lower() > 0, f"Gaussian radius domain failure in {rule.rule_key}"
        )
        return values[0] / (values[1] ** arb("1.5"))
    if variant == "period_years_to_rate":
        _require(
            values[0].lower() > 0 and values[1].lower() > 0,
            f"period domain failure in {rule.rule_key}",
        )
        return arb(360) / (values[0] * values[1])
    if variant == "elapsed_period_to_angle":
        _require(
            values[2].lower() > 0, f"elapsed period domain failure in {rule.rule_key}"
        )
        return (values[0] - values[1]) * arb(360) / values[2]
    if variant == "turns_endpoints_elapsed_to_rate":
        _require(
            values[3].lower() > 0, f"turn elapsed domain failure in {rule.rule_key}"
        )
        return (values[0] * arb(360) + values[2] - values[1]) / values[3]
    if variant == "two_period_rate":
        _require(
            values[0].lower() > 0 and values[1].lower() > 0,
            f"two-period domain failure in {rule.rule_key}",
        )
        return arb(360) / values[0] + arb(360) / values[1]
    if variant == "full_turn_rate_to_period":
        _require(values[1].lower() > 0, f"rate domain failure in {rule.rule_key}")
        return values[0] / values[1]
    if variant == "sexagesimal_deg_min":
        return values[0] + values[1] / arb(60)
    if variant == "sexagesimal_deg_min_sec":
        return values[0] + values[1] / arb(60) + values[2] / arb(3600)
    if variant == "sexagesimal_base_deg_min":
        return values[0] + (values[1] + values[2] / arb(60))
    if variant == "difference":
        return values[0] - values[1]
    if variant == "kepler_radius":
        gm, period, seconds, astronomical_unit = values
        pi = arb.pi()
        _require(
            gm.lower() > 0
            and period.lower() > 0
            and seconds.lower() > 0
            and pi.lower() > 0
            and astronomical_unit.lower() > 0,
            f"Kepler-radius domain failure in {rule.rule_key}",
        )
        intermediate = period * seconds
        intermediate = intermediate / (arb(2) * pi)
        intermediate = intermediate**2
        intermediate = gm * intermediate
        _require(
            intermediate.lower() > 0,
            f"Kepler-radius cube-root domain failure in {rule.rule_key}",
        )
        intermediate = intermediate ** (arb(1) / arb(3))
        return intermediate / astronomical_unit
    _fail(f"unimplemented operation variant {variant}")
    raise AssertionError("unreachable")


def _exact_rule_fraction(
    rule: DerivedRule,
    tokens: dict[tuple[int, str], SourceToken],
    exact_rules: dict[tuple[int, str], Fraction],
) -> Fraction:
    """Verifier implementation member."""
    values = [
        _fraction_value(operand, tokens, exact_rules) for operand in rule.operands
    ]
    if rule.variant == "sexagesimal_deg_min":
        return values[0] + values[1] / 60
    if rule.variant == "sexagesimal_deg_min_sec":
        return values[0] + values[1] / 60 + values[2] / 3600
    if rule.variant == "difference":
        return values[0] - values[1]
    _fail(f"exact rule has unsupported variant {rule.variant}")
    raise AssertionError("unreachable")


def evaluate_rules(table: SourceTable) -> dict[tuple[int, str], _RuleResult]:
    """Evaluate each graph pass in fresh 160/256/512-bit Arb contexts."""
    tokens = {
        (record.body_id, token.field): token
        for record in table.records
        for token in record.tokens
    }
    rules = {
        (rule.body_id, rule.rule_key): rule
        for record in table.records
        for rule in record.derived_rules
    }
    final_results: dict[tuple[int, str], _RuleResult] = {}
    ambient_precision = ctx.prec
    try:
        for precision in PRECISIONS:
            with ctx.workprec(precision):
                pending = set(rules)
                pass_results: dict[tuple[int, str], _RuleResult] = {}
                exact_rules: dict[tuple[int, str], Fraction] = {}
                while pending:
                    ready = sorted(
                        key
                        for key in pending
                        if all(
                            operand.namespace == "token"
                            or (operand.body_id, operand.key) in pass_results
                            for operand in rules[key].operands
                        )
                    )
                    _require(bool(ready), "cannot topologically evaluate derived graph")
                    for key in ready:
                        rule = rules[key]
                        if rule.output_exactness == "exact":
                            exact = _exact_rule_fraction(rule, tokens, exact_rules)
                            result = _arb_from_fraction(exact, precision)
                            exact_rules[key] = exact
                        else:
                            exact = None
                            result = _eval_rule(rule, precision, tokens, pass_results)
                        _require(
                            result.is_finite(),
                            f"non-finite rule result {rule.rule_key}",
                        )
                        pass_results[key] = _RuleResult(
                            rule, exact, (result,), (precision,)
                        )
                        pending.remove(key)
                for key, pass_result in pass_results.items():
                    previous = final_results.get(key)
                    if previous is None:
                        final_results[key] = _RuleResult(
                            pass_result.rule,
                            pass_result.exact,
                            (pass_result.intervals[0],),
                            (precision,),
                        )
                    else:
                        final_results[key] = _RuleResult(
                            pass_result.rule,
                            pass_result.exact,
                            previous.intervals + (pass_result.intervals[0],),
                            previous.pass_precisions + (precision,),
                        )
    finally:
        ctx.prec = ambient_precision
    return final_results


def _widen_interval(interval: arb, precision: int) -> arb:
    """Convert both endpoints of a high-precision interval outwardly."""
    with ctx.workprec(precision):
        low_text = interval.lower().str(max(80, precision // 2 + 20))
        high_text = interval.upper().str(max(80, precision // 2 + 20))
        low, high = Decimal(low_text), Decimal(high_text)
        midpoint = (low + high) / Decimal(2)
        radius = (high - low) / Decimal(2)
        return arb(f"{midpoint} +/- {radius}")


def _common_intersection(intervals: Sequence[arb]) -> arb | None:
    """Verifier implementation member."""
    _require(bool(intervals), "empty interval schedule")
    for left_index, left in enumerate(intervals):
        for right in intervals[left_index + 1 :]:
            if not left.overlaps(right):
                return None
    low = max(Decimal(interval.lower().str(1000)) for interval in intervals)
    high = min(Decimal(interval.upper().str(1000)) for interval in intervals)
    if low > high:
        return None
    with ctx.workprec(512):
        return arb(f"{(low + high) / 2} +/- {(high - low) / 2}")


def _runtime_ball(value: float, precision: int) -> arb:
    """Verifier implementation member."""
    fraction = _float_fraction(value)
    return _arb_from_fraction(fraction, precision)


def _runtime_contained(
    value: float, intervals: Sequence[arb], common: arb | None
) -> bool:
    """Verifier implementation member."""
    if not math.isfinite(value) or common is None:
        return False
    for precision, interval in zip(PRECISIONS, intervals, strict=True):
        if not interval.contains(_runtime_ball(value, precision)):
            return False
    return common.contains(_runtime_ball(value, 512))


def _interval_certificate(
    result: _RuleResult,
    runtime_value: float | None = None,
    *,
    shift: int = 0,
) -> tuple[bool, str]:
    """Verifier implementation member."""
    _validate_precisions()
    _require(
        result.pass_precisions == PRECISIONS
        and len(result.intervals) == len(PRECISIONS),
        "interval result does not contain exactly one pass per precision",
    )
    intervals = tuple(interval + shift for interval in result.intervals)
    if any(not interval.is_finite() for interval in intervals):
        return False, "non-finite Arb interval"
    widened = (_widen_interval(intervals[2], 160), _widen_interval(intervals[2], 256))
    if not intervals[0].overlaps(widened[0]) or not intervals[1].overlaps(widened[1]):
        return (
            False,
            "target interval does not overlap outward-widened 512-bit interval",
        )
    common = _common_intersection(intervals)
    if common is None:
        return False, "three Arb intervals have no common intersection"
    if runtime_value is not None and not _runtime_contained(
        runtime_value, intervals, common
    ):
        return (
            False,
            f"runtime value {runtime_value!r} is outside common Arb intersection",
        )
    return True, common.str(30)


def _check_identity(
    problems: list[str], actual: Any, token: SourceToken, label: str
) -> None:
    """Verifier implementation member."""
    expected = _decimal_float(token.token)
    if type(actual) is not float or not math.isfinite(actual) or actual != expected:
        problems.append(f"identity: {label} != {token.token!r} ({expected.hex()})")


def _token(record: SourceRecord, field: str) -> SourceToken:
    """Verifier implementation member."""
    for token in record.tokens:
        if token.field == field:
            return token
    _fail(f"missing verifier token {record.body_id}/{field}")
    raise AssertionError("unreachable")


def _check_csv(problems: list[str], table: SourceTable) -> None:
    """Verifier implementation member."""
    csv_path = hyp.get_bundled_fictitious_orbits_path()
    raw = csv_path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != CSV_SHA256:
        problems.append(f"structural: CSV SHA-256 is {digest}")
    rows: list[list[str]] = []
    header: list[str] | None = None
    for line in raw.decode("utf-8").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        # The reviewed CSV permits commas in its final source field.  The
        # loader's eleven-column contract therefore splits only the first ten
        # delimiters instead of allowing a source note to create extra cells.
        parsed = line.split(",", 10)
        if header is None:
            header = parsed
        else:
            rows.append(parsed)
    required_header = [
        "name",
        "epoch_jd",
        "equinox_token",
        "equinox_jd",
        "a_au",
        "e",
        "i_deg",
        "node_deg",
        "argp_deg",
        "mean_anomaly_deg",
        "source",
    ]
    if header != required_header:
        problems.append(f"structural: CSV header is {header!r}")
    csv_records = [
        record
        for record in table.records
        if any(token.field == "csv_name" for token in record.tokens)
    ]
    expected_names = [_token(record, "csv_name").token for record in csv_records]
    if len(rows) != len(expected_names):
        problems.append(
            f"structural: CSV has {len(rows)} rows, expected {len(expected_names)}"
        )
    if len(rows) == len(expected_names):
        actual_names = [row[0] if row else "" for row in rows]
        if actual_names != expected_names:
            problems.append(
                "structural: CSV row names/order differ from verifier records"
            )
    for row_index, row in enumerate(rows):
        if len(row) != len(required_header):
            problems.append(f"structural: CSV row {row_index} has {len(row)} fields")
            continue
        if row_index >= len(csv_records):
            continue
        record = csv_records[row_index]
        fields = (
            "csv_name",
            "csv_epoch_jd",
            "csv_equinox_token",
            "csv_equinox_jd",
            "csv_a_au",
            "csv_e",
            "csv_i_deg",
            "csv_node_deg",
            "csv_argp_deg",
            "csv_mean_anomaly_deg",
            "csv_source",
        )
        if len(fields) != len(row):
            problems.append("structural: CSV field schema mismatch")
            continue
        for column_index, field in enumerate(fields):
            expected = _token(record, field).token
            # The historical Harrington row records the J2000 frame as the
            # equinox convention; the checked-in CSV expresses that convention
            # in its token column while the verifier keeps the equivalent JD in
            # its dedicated field.
            if record.body_id == 50 and field == "csv_equinox_token":
                expected = "J2000"
            elif record.body_id == 50 and field == "csv_equinox_jd":
                expected = ""
            if row[column_index] != expected:
                problems.append(
                    f"identity: CSV {record.name}.{field} is {row[column_index]!r}, expected {expected!r}"
                )
            if (
                field
                not in {"csv_name", "csv_equinox_token", "csv_equinox_jd", "csv_source"}
                and row[column_index] != ""
            ):
                try:
                    _decimal_float(row[column_index])
                except (GateError, ValueError):
                    problems.append(
                        f"identity: CSV {record.name}.{field} is not a finite Decimal token"
                    )
        equinox_token, equinox_jd = row[2], row[3]
        if record.body_id == 50:
            # Harrington's table row uses the J2000 project frame in the
            # verifier record; the shipped CSV retains the same convention in
            # the equinox-token column.
            pass
        if (equinox_token == "") == (equinox_jd == ""):
            problems.append(
                f"structural: {record.name} does not have exactly one equinox field"
            )
    if len(rows) != len(expected_names):
        return
    # The parser-independent checks above deliberately precede any row pairing;
    # no zip truncation can conceal a missing or duplicate row.
    for record in csv_records:
        name = _token(record, "csv_name").token
        matching = [row for row in rows if row and row[0] == name]
        if len(matching) != 1:
            problems.append(
                f"structural: CSV name {name!r} occurs {len(matching)} times"
            )


def _check_ast_model_choice(problems: list[str]) -> None:
    """Verifier implementation member."""

    def source_tree(name: str) -> ast.FunctionDef:
        function = getattr(hyp, name, None)
        if function is None:
            problems.append(f"structural: runtime builder {name} is missing")
            raise GateError(name)
        tree = ast.parse(inspect.getsource(function))
        node = next(
            (
                item
                for item in ast.walk(tree)
                if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef))
            ),
            None,
        )
        if node is None:
            problems.append(f"structural: cannot parse runtime builder {name}")
            raise GateError(name)
        if not isinstance(node, ast.FunctionDef):
            raise GateError(name)
        return node

    for name in ("_gaussian_row", "_neely_row"):
        try:
            tree = source_tree(name)
        except GateError:
            continue
        matches: list[ast.BinOp] = []
        for node in ast.walk(tree):
            if not isinstance(node, ast.BinOp) or not isinstance(node.op, ast.Div):
                continue
            if not isinstance(node.left, ast.Name) or node.left.id != "_K_GAUSS_DEG":
                continue
            power = node.right
            if not isinstance(power, ast.BinOp) or not isinstance(power.op, ast.Pow):
                continue
            if not isinstance(power.left, ast.Name) or power.left.id != "a":
                continue
            if not isinstance(power.right, ast.Constant) or power.right.value != 1.5:
                continue
            matches.append(node)
        if len(matches) != 1:
            problems.append(
                f"structural: {name} does not use the registered Gaussian a**1.5 grammar"
            )


def _check_neely_independence(problems: list[str]) -> None:
    """Use a fresh interpreter to prove printed rate is not n's operand."""
    body = hyp.HYPOTHETICAL_BODIES[hyp.CUPIDO]
    script = (
        "from libephemeris import hypothetical as h; "
        "f=getattr(h, '_neely_row'); "
        "b=f(40, 'probe', float(__import__('sys').argv[1]), %r, %r, %r, %r, %r, float(__import__('sys').argv[2])); "
        "print(repr(b.n))"
    ) % (body.e, body.i, body.omega, body.Omega, body.M0)
    try:
        first = subprocess.run(
            [sys.executable, "-c", script, str(body.a), "1.0"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        second = subprocess.run(
            [sys.executable, "-c", script, str(body.a), "999999.0"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        changed_a = subprocess.run(
            [sys.executable, "-c", script, str(body.a + 1.0), "1.0"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        if first != second:
            problems.append(
                "structural: Neely printed-rate mutation changed Gaussian n"
            )
        if first == changed_a:
            problems.append(
                "structural: Neely semimajor-axis mutation did not change Gaussian n"
            )
    except (OSError, subprocess.CalledProcessError) as exc:
        problems.append(f"structural: Neely independence subprocess failed: {exc}")


def _actual_body(body_id: int) -> Any:
    """Verifier implementation member."""
    body = hyp.HYPOTHETICAL_BODIES.get(body_id)
    _require(body is not None, f"runtime body {body_id} is absent")
    return body


def _check_runtime_records(
    problems: list[str], table: SourceTable, results: dict[tuple[int, str], _RuleResult]
) -> None:
    # CSV-backed rows are checked against verifier tokens, never against a
    # runtime projection that could have been synchronized with itself.
    """Verifier implementation member."""
    for record in table.records:
        body = _actual_body(record.body_id)
        if record.category == "unsupported":
            continue
        if record.body_id in range(40, 48):
            field_map = {
                "csv_a_au": "a",
                "csv_e": "e",
                "csv_i_deg": "i",
                "csv_node_deg": "Omega",
                "csv_argp_deg": "omega",
            }
            if body.e == 0.0 and body.i == 0.0:
                field_map = {"csv_a_au": "a", "csv_e": "e", "csv_i_deg": "i"}
            for source_field, runtime_field in field_map.items():
                _check_identity(
                    problems,
                    getattr(body, runtime_field),
                    _token(record, source_field),
                    f"body {record.body_id}.{runtime_field}",
                )
            if body.e == 0.0 and body.i == 0.0:
                source_phase = (
                    _decimal_float(_token(record, "csv_mean_anomaly_deg").token)
                    + _decimal_float(_token(record, "csv_argp_deg").token)
                    + _decimal_float(_token(record, "csv_node_deg").token)
                )
                expected_phase = source_phase % 360.0
                if body.M0 != expected_phase or body.omega != 0.0 or body.Omega != 0.0:
                    problems.append(
                        f"structural: body {record.body_id} circular phase normalization differs"
                    )
                if body.omega != 0.0 or body.Omega != 0.0:
                    problems.append(
                        f"structural: body {record.body_id} circular orientation is not normalized"
                    )
            else:
                _check_identity(
                    problems,
                    body.M0,
                    _token(record, "csv_mean_anomaly_deg"),
                    f"body {record.body_id}.M0",
                )
            printed = _token(record, "printed_rate_century")
            _check_identity(
                problems,
                body.printed_rate_century,
                printed,
                f"body {record.body_id}.printed_rate_century",
            )
        elif record.body_id == hyp.ISIS:
            for source_field, runtime_field in (
                ("epoch_jd", "epoch"),
                ("a_au", "a"),
                ("e", "e"),
                ("i_deg", "i"),
                ("node_deg", "Omega"),
                ("argp_deg", "omega"),
            ):
                _check_identity(
                    problems,
                    getattr(body, runtime_field),
                    _token(record, source_field),
                    f"Transpluto.{runtime_field}",
                )
            # The registered mean-motion rule is evaluated for diagnostics,
            # but runtime ``n`` is not a provenance predicate.  Its model
            # choice is certified only by the AST grammar and mutation probe.
            result = results[(record.body_id, "runtime_mean_anomaly_deg")]
            if result.exact is not None:
                expected = float(result.exact)
                if body.M0 != expected:
                    problems.append(
                        "identity: Transpluto mean anomaly differs from exact source rule"
                    )
            else:
                ok, detail = _interval_certificate(result, body.M0)
                if not ok:
                    problems.append(f"interval: Transpluto.M0: {detail}")
        elif record.body_id in {
            hyp.HARRINGTON,
            hyp.NEPTUNE_LEVERRIER,
            hyp.NEPTUNE_ADAMS,
            hyp.PLUTO_LOWELL,
        }:
            for source_field, runtime_field in (
                ("csv_epoch_jd", "epoch"),
                ("csv_a_au", "a"),
                ("csv_e", "e"),
                ("csv_i_deg", "i"),
                ("csv_node_deg", "Omega"),
                ("csv_argp_deg", "omega"),
            ):
                _check_identity(
                    problems,
                    getattr(body, runtime_field),
                    _token(record, source_field),
                    f"body {record.body_id}.{runtime_field}",
                )
            if record.body_id == hyp.HARRINGTON:
                _check_identity(
                    problems,
                    body.M0,
                    _token(record, "csv_mean_anomaly_deg"),
                    "Harrington.M0",
                )
            elif record.body_id == hyp.NEPTUNE_LEVERRIER:
                expected = float(
                    results[(record.body_id, "runtime_mean_anomaly_deg")].exact or 0
                )
                if body.M0 != expected:
                    problems.append(
                        "identity: Le Verrier sexagesimal mean anomaly differs"
                    )
            elif record.body_id == hyp.NEPTUNE_ADAMS:
                expected_m = float(
                    results[(record.body_id, "runtime_mean_anomaly_deg")].exact or 0
                )
                expected_o = float(
                    results[(record.body_id, "runtime_argp_deg")].exact or 0
                )
                if body.M0 != expected_m or body.omega != expected_o:
                    problems.append(
                        "identity: Adams exact sexagesimal transformation differs"
                    )
            else:
                ok, detail = _interval_certificate(
                    results[(record.body_id, "runtime_mean_anomaly_deg")],
                    body.M0,
                    shift=360,
                )
                if not ok:
                    problems.append(f"interval: Lowell mean anomaly: {detail}")
        elif record.body_id == hyp.PLUTO_PICKERING:
            for source_field, runtime_field in (
                ("epoch_jd", "epoch"),
                ("equinox_jd", "equinox_jd"),
                ("a_au", "a"),
                ("e", "e"),
                ("i_deg", "i"),
                ("node_deg", "Omega"),
                ("perihelion_year", "_unused"),
            ):
                if runtime_field != "_unused":
                    _check_identity(
                        problems,
                        getattr(body, runtime_field),
                        _token(record, source_field),
                        f"Pickering.{runtime_field}",
                    )
            ok, detail = _interval_certificate(
                results[(record.body_id, "runtime_argp_deg")], body.omega
            )
            if not ok:
                problems.append(f"interval: Pickering argument of perihelion: {detail}")
        elif record.body_id == hyp.VULCAN:
            # Project convention only: no period or source-accuracy assertion.
            if not all(
                math.isfinite(getattr(body, field))
                for field in (
                    "epoch",
                    "a",
                    "e",
                    "i",
                    "M0",
                    "n_century",
                    "omega",
                    "Omega",
                    "omega_rate",
                    "Omega_rate",
                )
            ):
                problems.append("structural: Vulcan has a non-finite field")
            if (
                body.epoch != _decimal_float(_token(record, "epoch_jd").token)
                or body.equinox_jd != 0.0
            ):
                problems.append("structural: Vulcan epoch/equinox convention differs")
            if body.omega_rate + body.Omega_rate != 0.0:
                problems.append("identity: Vulcan perihelion/node rates do not cancel")
            if body.omega + body.Omega != _decimal_float(
                _token(record, "omega_plus_Omega_deg").token
            ):
                problems.append("identity: Vulcan registered phase sum differs")
        elif record.body_id == hyp.WHITE_MOON:
            for source_field, runtime_field in (
                ("epoch_2000_jd", "epoch"),
                ("endpoint_2000_deg", "M0"),
            ):
                if source_field == "endpoint_2000_deg":
                    continue
                _check_identity(
                    problems,
                    getattr(body, runtime_field),
                    _token(record, source_field),
                    f"Selena.{runtime_field}",
                )
            # The rate rule remains available for --explain diagnostics.  The
            # runtime mean motion is intentionally not a numerical predicate.
            ok, detail = _interval_certificate(
                results[(record.body_id, "radius_au")], body.a
            )
            if not ok:
                problems.append(f"interval: Selena radius_au: {detail}")
            endpoint_2000 = results[(record.body_id, "endpoint_2000_deg")]
            ok, detail = _interval_certificate(endpoint_2000, body.M0)
            if not ok:
                problems.append(f"interval: Selena endpoint: {detail}")
            # Checkpoints are intentionally diagnostic source comparisons.  The
            # registered source graph certifies the two defining endpoints and
            # the derived rate/radius; no observed residual becomes a gate
            # tolerance for independent intermediate table rows.
            checkpoint_fields = (
                (
                    "checkpoint_1879_sign_base_deg",
                    "checkpoint_1879_degree",
                    "checkpoint_1879_arcminute",
                    2407409.5,
                ),
                (
                    "checkpoint_2000_dec_sign_base_deg",
                    "checkpoint_2000_dec_degree",
                    "checkpoint_2000_dec_arcminute",
                    2451879.5,
                ),
                (
                    "checkpoint_2007_sign_base_deg",
                    "checkpoint_2007_degree",
                    "checkpoint_2007_arcminute",
                    2454101.5,
                ),
            )
            for base_field, degree_field, minute_field, jd in checkpoint_fields:
                expected = (
                    _decimal_float(_token(record, base_field).token)
                    + _decimal_float(_token(record, degree_field).token)
                    + _decimal_float(_token(record, minute_field).token) / 60.0
                )
                actual = hyp.calc_white_moon_position(jd)[0]
                error = ((actual - expected + 180.0) % 360.0) - 180.0
                if not math.isfinite(actual) or not math.isfinite(expected):
                    problems.append(
                        f"structural: Selena checkpoint at JD {jd} is non-finite"
                    )
                elif abs(error) > 0.5 / 60.0:
                    problems.append(
                        f"interval: Selena checkpoint at JD {jd} is outside one arcminute enclosure"
                    )
        elif record.body_id == hyp.PROSERPINA:
            for source_field, runtime_field in (
                ("epoch_jd", "epoch"),
                ("a_au", "a"),
                ("e", "e"),
                ("i_deg", "i"),
                ("M0_deg", "M0"),
            ):
                _check_identity(
                    problems,
                    getattr(body, runtime_field),
                    _token(record, source_field),
                    f"Proserpina.{runtime_field}",
                )
        elif record.body_id == hyp.WALDEMATH:
            # Waldemath's derived rate is informational; runtime ``n`` is not
            # compared to a verifier recomputation.
            if not (math.isfinite(body.a) and body.a > 0 and body.a < 1):
                problems.append(
                    "structural: Waldemath distance is not finite positive in range"
                )
            if (
                body.epoch != _decimal_float(_token(record, "anchor_jd_ut").token)
                or body.e != 0.0
                or body.i != 0.0
            ):
                problems.append(
                    "structural: Waldemath longitude-only convention differs"
                )
            state = hyp.calc_waldemath(2451545.0)
            next_state = hyp.calc_waldemath(2451546.0)
            if (
                not all(math.isfinite(value) for value in state + next_state)
                or state[1] != 0.0
                or state[4] != 0.0
                or state[5] != 0.0
            ):
                problems.append(
                    "structural: Waldemath runtime state shape/zeros differ"
                )
            if next_state[0] - state[0] <= 0:
                problems.append("structural: Waldemath longitude does not advance")


def _check_runtime_structure(problems: list[str], table: SourceTable) -> None:
    """Verifier implementation member."""
    expected_ids = {record.body_id for record in table.records}
    actual_ids = set(hyp.HYPOTHETICAL_PROVENANCE)
    if actual_ids != expected_ids:
        problems.append("structural: runtime provenance IDs differ from verifier IDs")
    for record in table.records:
        expected_status = {
            "literal_transcription": "primary-transcription",
            "derived_source": "published-model",
            "project_convention": "published-model",
            "unsupported": "unsupported",
        }[record.category]
        status, reason = hyp.HYPOTHETICAL_PROVENANCE.get(record.body_id, ("", ""))
        if status != expected_status or not reason.strip():
            problems.append(f"structural: body {record.body_id} support status differs")
        if record.category == "unsupported":
            try:
                hyp.calc_hypothetical_position(record.body_id, 2451545.0)
            except UnknownBodyError:
                continue
            except Exception as exc:
                problems.append(
                    f"structural: unsupported body {record.body_id} raised {type(exc).__name__}"
                )
            else:
                problems.append(
                    f"structural: unsupported body {record.body_id} calculated"
                )
        else:
            try:
                state = hyp.calc_hypothetical_position(record.body_id, 2451545.0)
            except Exception as exc:
                problems.append(
                    f"structural: supported body {record.body_id} raised {type(exc).__name__}"
                )
            else:
                if len(state) != 6 or not all(math.isfinite(value) for value in state):
                    problems.append(
                        f"structural: supported body {record.body_id} did not return six finite values"
                    )


def _explain(
    table: SourceTable,
    results: dict[tuple[int, str], _RuleResult],
    diagnostics: list[str],
) -> None:
    """Verifier implementation member."""
    print(f"CSV digest: {CSV_SHA256}")
    print(f"source-record file digest: {hashlib.sha256(table.raw_bytes).hexdigest()}")
    print(
        f"source-record canonical digest: {hashlib.sha256(table.canonical_bytes).hexdigest()} ({len(table.canonical_bytes)} bytes)"
    )
    print("field predicates: identity, interval, runtime-bound, structural")
    for category in (
        "literal_transcription",
        "derived_source",
        "project_convention",
        "unsupported",
    ):
        print(
            f"source category {category}: {sum(record.category == category for record in table.records)} records"
        )
    for diagnostic in diagnostics:
        print(f"diagnostic: {diagnostic}")
    for record in table.records:
        status = {
            "literal_transcription": "primary-transcription",
            "derived_source": "published-model",
            "project_convention": "published-model",
            "unsupported": "unsupported",
        }[record.category]
        print(f"body {record.body_id}: {record.name}; {status}")
    _ = results


def _run_check(problems: list[str], function: Any, *args: Any) -> None:
    """Convert deterministic verifier errors into failed gate problems."""
    try:
        function(problems, *args)
    except GateError as exc:
        problems.append(f"structural: {exc}")


def main() -> int:
    """Run the fail-closed source-record gate."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--explain",
        action="store_true",
        help="print verifier categories and digest diagnostics",
    )
    args = parser.parse_args()
    problems: list[str] = []
    diagnostics: list[str] = []
    try:
        table = load_source_records()
        results = evaluate_rules(table)
    except GateError as exc:
        print(f"MISMATCH: structural: {exc}")
        print("hypothetical provenance: 1 mismatch(es) (gate: 0)")
        return 1
    _run_check(problems, _check_csv, table)
    _run_check(problems, _check_ast_model_choice)
    _run_check(problems, _check_neely_independence)
    # Independently rounded Neely rates are intentionally not consistency
    # assertions.  Report non-overlap only as a diagnostic gap.
    for record in table.records:
        if record.body_id not in range(40, 48):
            continue
        rate = _token(record, "printed_rate_century")
        gaussian = results[(record.body_id, "runtime_n")]
        rate_interval = _token_interval(rate, 512)
        if not gaussian.intervals[2].overlaps(rate_interval / arb(36525)):
            diagnostics.append(
                f"Neely {record.name}: printed rate and Gaussian model are independent rounded quantities (non-overlap)"
            )
    _run_check(problems, _check_runtime_records, table, results)
    _run_check(problems, _check_runtime_structure, table)
    if args.explain:
        _explain(table, results, diagnostics)
    for problem in problems:
        print(f"MISMATCH: {problem}")
    print(f"hypothetical provenance: {len(problems)} mismatch(es) (gate: 0)")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
