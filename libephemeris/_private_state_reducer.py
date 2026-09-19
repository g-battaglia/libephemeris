# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Disconnected, diagnostic-only Sun/Moon reducer over owned snapshots.

Provenance:
    Project-authored adapter of the reviewed private Sun/Moon reducer contract
    (SHA-256 ad2d7ae454f90ab7a2378e6c2d93aca4aecbf14566bae254f5fde4f66fb09299).
    Numerical state reduction calls the existing project ``_fast_calc_core``;
    this module supplies admitted source routing, uncached frame dispatch,
    source/control observations and call-local cleanup. The transitive code
    and binary inventory is not complete, so results are diagnostics and
    never state-source receipts or evidence for an ordinary public call.
"""

from __future__ import annotations

import math
import os
import platform
import struct
import sys
from copy import deepcopy
from dataclasses import dataclass
from types import MappingProxyType
from typing import Any, Callable, Mapping, cast

import erfa
import numpy as np
import skyfield
import zstandard as zstd

from . import fast_calc
from ._leb2_snapshot import (
    _LEB2SnapshotReader,
    _SnapshotClosedError,
    _SnapshotCorruptionError,
    _SnapshotIntegrityError,
    _open_production_snapshot,
)
from ._time_snapshot import (
    _DefaultTimeDeclaration,
    _DefaultTimeEvaluation,
    _OwnedDefaultTime,
    _TimeSnapshotClosedError,
    _TimeSnapshotConfigurationError,
    _TimeSnapshotConstructionError,
    _TimeSnapshotDependencyError,
    _TimeSnapshotIntegrityError,
    _TimeSnapshotParseError,
    _TimeSnapshotUnsupportedError,
    _check_dependencies as _check_time_dependencies,
    _open_production_time,
)
from .constants import (
    EARTH,
    FLG_EQUATORIAL,
    FLG_SWIEPH,
    FLG_XYZ,
    JUPITER,
    MOON,
    SATURN,
    SUN,
)
from .exceptions import LEBCorruptionError
from .leb_format import COORD_ICRS_BARY, COORD_ICRS_BARY_SYSTEM
from .leb_reader import _split_epoch_compare
from .planets import _echo_request_bits, _implied_retflag_bits
from .precession_vondrak import vondrak_pn_matrix

_ROLES = (
    "base_core.leb2",
    "medium_core.leb2",
    "iers.npz",
    "historic_deltat.npy",
    "delta_t.npz",
)
_RAW_REQUEST = FLG_EQUATORIAL | FLG_XYZ
_REQUESTS = frozenset((_RAW_REQUEST, _RAW_REQUEST | FLG_SWIEPH))
_BODY_ROLES = MappingProxyType(
    {
        SUN: COORD_ICRS_BARY,
        MOON: COORD_ICRS_BARY,
        EARTH: COORD_ICRS_BARY,
        JUPITER: COORD_ICRS_BARY_SYSTEM,
        SATURN: COORD_ICRS_BARY_SYSTEM,
    }
)


class _StateUnsupportedError(RuntimeError):
    """An input, coverage interval, or branch is outside the first domain."""


class _StateAssetUnavailableError(RuntimeError):
    """A required pinned asset could not be opened or read."""


class _StateAssetIntegrityError(RuntimeError):
    """A required asset's exact count or digest differs from its pin."""


class _StateAssetCorruptionError(RuntimeError):
    """An admitted child reports malformed or corrupt source data."""


class _StateDependencyError(RuntimeError):
    """The guarded numerical runtime differs from the reviewed environment."""


class _StateConfigurationError(RuntimeError):
    """The private default-time branch declaration failed."""


class _StateInvariantError(RuntimeError):
    """A private source, control, or numerical invariant failed."""


class _StateClosedError(RuntimeError):
    """A privately owned source was closed or invalidated during reduction."""


@dataclass(frozen=True, slots=True)
class _CandidateDecision:
    """One concrete child's role, channel, bounds, and selection decision."""

    tier: str
    tag: tuple[str, str, int, str]
    channel: int | None
    coverage_bits: tuple[bytes, bytes] | None
    covers: bool


@dataclass(frozen=True, slots=True)
class _SourceEvent:
    """One observed body/split/nutation read from a concrete owned child."""

    scope_body: int
    kind: str
    source_body: int | None
    jd_bits: bytes
    offset_bits: bytes
    decisions: tuple[_CandidateDecision, ...]
    selected_tag: tuple[str, str, int, str]
    value_bits: tuple[bytes, ...]


@dataclass(frozen=True, slots=True)
class _ControlEvent:
    """An immutable branch observation that cannot feed the computation."""

    kind: str
    fields: tuple[tuple[str, str, bytes | int | str | bool], ...]


@dataclass(frozen=True, slots=True)
class _PrivateBodyState:
    """One body's native result and independent source/control history."""

    body: int
    raw_request: int
    normalized_request: int
    ut_bits: bytes
    tt_bits: bytes
    time_seconds_bits: bytes
    time_days_bits: bytes
    time_tag: tuple[tuple[str, int, str], ...]
    words: tuple[float, float, float, float, float, float]
    word_bits: tuple[bytes, bytes, bytes, bytes, bytes, bytes]
    inner_flag: int
    returned_flag: int
    sources: tuple[_SourceEvent, ...]
    controls: tuple[_ControlEvent, ...]


@dataclass(frozen=True, slots=True)
class _PrivateStatePair:
    """Two separate private states with shared owned time evaluation."""

    moon: _PrivateBodyState
    sun: _PrivateBodyState
    time: _DefaultTimeEvaluation

    @property
    def source_status(self) -> str:
        """Keep this incomplete source graph explicitly diagnostic."""
        return "diagnostic-only"


def _bits(value: float) -> bytes:
    """Return the complete binary64 word, preserving a signed zero."""
    return struct.pack("!d", value)


def _native_finite_words(values: object, count: int) -> tuple[float, ...]:
    """Convert returned scalars exactly as the public LEB result does."""
    try:
        result = tuple(float(word) for word in values)  # type: ignore[attr-defined]
    except (TypeError, ValueError, OverflowError) as exc:
        raise _StateInvariantError("source words are not native finite floats") from exc
    if len(result) != count or not all(math.isfinite(word) for word in result):
        raise _StateInvariantError("source word count or finiteness differs")
    return result


def _control_word(value: object) -> tuple[str, bytes | int | str | bool]:
    """Keep a branch fact's type and exact float bits in immutable form."""
    if type(value) is bool:
        return "bool", value
    if type(value) is int:
        return "int", value
    if type(value) is float:
        return "float", _bits(value)
    if type(value) is str:
        return "str", value
    raise _StateInvariantError("unrecordable private control fact")


def _check_runtime() -> None:
    """Guard direct runtime versions; a full source manifest is still absent."""
    if (
        sys.version_info[:3] != (3, 12, 13)
        or sys.platform != "darwin"
        or platform.machine() != "arm64"
        or np.__version__ != "2.3.5"
        or skyfield.__version__ != "1.54"
        or erfa.__version__ != "2.0.1.5"
        or zstd.__version__ != "0.25.0"
        or zstd.ZSTD_VERSION != (1, 5, 7)
    ):
        raise _StateDependencyError("private reducer numerical runtime differs")
    try:
        _check_time_dependencies()
    except _TimeSnapshotDependencyError as exc:
        raise _StateDependencyError("private time numerical source differs") from exc


def _paths_once(paths: Mapping[str, str | os.PathLike[str]]) -> dict[str, str]:
    """Validate types and resolve each locator once before any asset opens."""
    if not isinstance(paths, Mapping):
        raise TypeError("paths must be a mapping")
    if any(type(key) is not str for key in paths):
        raise TypeError("path roles must be native strings")
    converted: dict[str, str] = {}
    for key, value in paths.items():
        if not (type(value) is str or isinstance(value, os.PathLike)):
            raise TypeError("asset locators must be string paths")
        resolved = os.fspath(value)
        if type(resolved) is not str:
            raise TypeError("asset locators must resolve to string paths")
        converted[key] = resolved
    if set(paths) != set(_ROLES):
        raise ValueError("the five fixed asset roles are required")
    return {role: converted[role] for role in _ROLES}


def _translate(exc: Exception) -> RuntimeError:
    """Map prerequisite failures without allowing deflector catch clauses to hide them."""
    if isinstance(
        exc,
        (
            _StateUnsupportedError,
            _StateAssetUnavailableError,
            _StateAssetIntegrityError,
            _StateAssetCorruptionError,
            _StateDependencyError,
            _StateConfigurationError,
            _StateInvariantError,
            _StateClosedError,
        ),
    ):
        return exc
    if isinstance(exc, (_SnapshotClosedError, _TimeSnapshotClosedError)):
        return _StateClosedError(str(exc))
    if isinstance(
        exc, (_SnapshotCorruptionError, _TimeSnapshotParseError, LEBCorruptionError)
    ):
        return _StateAssetCorruptionError(str(exc))
    if isinstance(exc, _TimeSnapshotConstructionError):
        return _StateInvariantError(str(exc))
    if isinstance(exc, (_SnapshotIntegrityError, _TimeSnapshotIntegrityError)):
        if isinstance(exc.__cause__, OSError):
            return _StateAssetUnavailableError(str(exc))
        return _StateAssetIntegrityError(str(exc))
    if isinstance(exc, _TimeSnapshotDependencyError):
        return _StateDependencyError(str(exc))
    if isinstance(exc, _TimeSnapshotConfigurationError):
        return _StateConfigurationError(str(exc))
    if isinstance(exc, _TimeSnapshotUnsupportedError):
        return _StateUnsupportedError(str(exc))
    if isinstance(exc, OSError):
        return _StateAssetUnavailableError(str(exc))
    return _StateInvariantError(f"private reducer failed: {exc}")


class _ObservedTierReader:
    """Route base before medium and bind every observed read to its child."""

    def __init__(self, base: _LEB2SnapshotReader, medium: _LEB2SnapshotReader):
        self.children = (("base", base), ("medium", medium))
        self.scope_body = MOON
        self.events: list[_SourceEvent] = []
        metadata: dict[int, object] = {}
        for body_id, expected_channel in _BODY_ROLES.items():
            for _tier, child in self.children:
                if child.production_asset_tag is None:
                    raise _StateAssetIntegrityError("test-only LEB child")
                entry = child._bodies.get(body_id)
                if entry is None:
                    continue
                if entry.coord_type != expected_channel:
                    raise _StateUnsupportedError("candidate LEB channel differs")
                metadata.setdefault(body_id, entry)
        self._bodies = MappingProxyType(metadata)

    def _select(
        self, body_id: int, jd: float, offset: float
    ) -> tuple[_LEB2SnapshotReader, tuple[_CandidateDecision, ...]]:
        """Retain every candidate decision, including rejected child coverage."""
        decisions: list[_CandidateDecision] = []
        selected: _LEB2SnapshotReader | None = None
        for tier, child in self.children:
            tag = child.production_asset_tag
            if tag is None or tag[0] != tier:
                raise _StateAssetIntegrityError("selected child has no fixed tag")
            entry = child._bodies.get(body_id)
            coverage = child.body_coverage(body_id)
            channel = entry.coord_type if entry is not None else None
            expected = _BODY_ROLES.get(body_id)
            if entry is not None and channel != expected:
                raise _StateUnsupportedError("candidate channel differs")
            covers = bool(
                coverage is not None
                and _split_epoch_compare(jd, offset, coverage[0]) >= 0
                and _split_epoch_compare(jd, offset, coverage[1]) <= 0
            )
            decisions.append(
                _CandidateDecision(
                    tier,
                    tag,
                    channel,
                    (_bits(coverage[0]), _bits(coverage[1])) if coverage else None,
                    covers,
                )
            )
            if selected is None and covers:
                selected = child
        if selected is None:
            raise _StateUnsupportedError("no admitted body coverage")
        return selected, tuple(decisions)

    def has_body(self, body_id: int) -> bool:
        """Report only inspected source roles to the shared reducer."""
        return body_id in self._bodies

    def _body(
        self, body_id: int, jd: float, offset: float
    ) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
        if not all(math.isfinite(x) for x in (jd, offset)):
            raise _StateUnsupportedError("nonfinite source epoch")
        try:
            selected, decisions = self._select(body_id, jd, offset)
            pos, vel = (
                selected.eval_body(body_id, jd)
                if offset == 0.0
                else selected._eval_body_split(body_id, jd, offset)
            )
            words = _native_finite_words((*pos, *vel), 6)
            tag = selected.production_asset_tag
            if tag is None:
                raise _StateAssetIntegrityError("selected child lost admission")
        except Exception as exc:
            translated = _translate(exc)
            if translated is exc:
                raise
            raise translated from exc
        self.events.append(
            _SourceEvent(
                self.scope_body,
                "body_split" if offset != 0.0 else "body",
                body_id,
                _bits(jd),
                _bits(offset),
                decisions,
                tag,
                tuple(_bits(word) for word in words),
            )
        )
        return (words[0], words[1], words[2]), (words[3], words[4], words[5])

    def eval_body(
        self, body_id: int, jd: float
    ) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
        """Observe a one-part body read."""
        return self._body(body_id, jd, 0.0)

    def _eval_body_split(
        self, body_id: int, jd: float, offset: float
    ) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
        """Observe an exact two-part retarded body read."""
        return self._body(body_id, jd, offset)

    def has_nutation(self) -> bool:
        """Report fixed children's available auxiliary channel."""
        return any(child.has_nutation() for _tier, child in self.children)

    def eval_nutation(self, jd_tt: float) -> tuple[float, float]:
        """Route nutation separately, retaining both tier coverage decisions."""
        decisions: list[_CandidateDecision] = []
        selected: _LEB2SnapshotReader | None = None
        for tier, child in self.children:
            tag = child.production_asset_tag
            if tag is None:
                raise _StateAssetIntegrityError("nutation child has no fixed tag")
            bounds = child.jd_range
            nut = child._nutation
            nut_bounds = (nut.jd_start, nut.jd_end) if nut is not None else None
            covers = bounds[0] <= jd_tt <= bounds[1] and child.has_nutation()
            decisions.append(
                _CandidateDecision(
                    tier,
                    tag,
                    None,
                    (_bits(nut_bounds[0]), _bits(nut_bounds[1]))
                    if nut_bounds is not None
                    else None,
                    covers,
                )
            )
            if selected is None and covers:
                selected = child
        if selected is None:
            raise _StateUnsupportedError("no admitted nutation coverage")
        selected_nut = selected._nutation
        if (
            selected_nut is None
            or not selected_nut.jd_start <= jd_tt <= selected_nut.jd_end
        ):
            raise _StateUnsupportedError("selected nutation interval does not cover TT")
        try:
            angles = _native_finite_words(selected.eval_nutation(jd_tt), 2)
        except Exception as exc:
            raise _translate(exc) from exc
        tag = selected.production_asset_tag
        if tag is None:
            raise _StateAssetIntegrityError("nutation child lost admission")
        self.events.append(
            _SourceEvent(
                self.scope_body,
                "nutation",
                None,
                _bits(jd_tt),
                _bits(0.0),
                tuple(decisions),
                tag,
                tuple(_bits(word) for word in angles),
            )
        )
        return angles[0], angles[1]


def _restore_thread_slot(slot: object, previous: dict[str, object]) -> None:
    """Restore a thread-local's precise pre-call dictionary, including nesting."""
    values = vars(slot)
    values.clear()
    values.update(previous)


def _run_body(
    reader: _ObservedTierReader,
    time: _DefaultTimeEvaluation,
    tt: float,
    body: int,
    raw_request: int,
    normalized: int,
    time_tag: tuple[tuple[str, int, str], ...],
    *,
    test_hook: Callable[[_ControlEvent], None] | None = None,
) -> _PrivateBodyState:
    """Run the shared native reducer with call-local source and frame evidence."""
    reader.scope_body = body
    reader.events.clear()
    controls: list[_ControlEvent] = []

    def record(kind: str, fields: dict[str, object]) -> None:
        event = _ControlEvent(
            kind,
            tuple(
                sorted((key, *_control_word(value)) for key, value in fields.items())
            ),
        )
        controls.append(event)
        if test_hook is not None:
            original = deepcopy(event)
            test_hook(event)
            if event != original:
                raise _StateInvariantError("private control event was changed")
        if kind in ("light_time_early_exit", "deflection_zero_geocentric_vector"):
            raise _StateUnsupportedError("private zero-distance branch is unresolved")

    def frame(
        jd: float,
    ) -> tuple[tuple[tuple[float, float, float], ...], float, float, float]:
        dpsi, deps = reader.eval_nutation(jd)
        pn, eps = vondrak_pn_matrix(jd, dpsi, deps, _uncached=True)
        if not all(math.isfinite(x) for row in pn for x in row) or not math.isfinite(
            eps
        ):
            raise _StateInvariantError("private frame is nonfinite")
        record("frame_selected", {"jd": jd, "dpsi": dpsi, "deps": deps})
        return pn, dpsi, deps, eps

    previous_active = vars(fast_calc._active_local).copy()
    previous_private = vars(fast_calc._private_state_local).copy()
    try:
        fast_calc._private_state_local.frame = frame
        fast_calc._private_state_local.control = record
        result, inner_flag = fast_calc._fast_calc_core(
            cast(Any, reader), tt, time.ut, body, normalized
        )
    finally:
        _restore_thread_slot(fast_calc._private_state_local, previous_private)
        _restore_thread_slot(fast_calc._active_local, previous_active)

    if type(inner_flag) is not int or inner_flag != normalized:
        raise _StateInvariantError("inner reducer flag differs from normalized request")
    words = _native_finite_words(result, 6)
    kinds = {event.kind for event in controls}
    required = {
        "pipeline_icrs",
        "light_time_iteration",
        "deflection_branch",
        "aberration_branch",
        "frame_request",
        "frame_selected",
    }
    if not required <= kinds or not any(e.kind == "nutation" for e in reader.events):
        raise _StateInvariantError("private branch or nutation evidence missing")
    if any(e.kind.endswith("read_skipped") for e in controls):
        raise _StateInvariantError("deflector source read was silently skipped")
    if body == MOON:
        if any(e.kind == "deflector_selected" for e in controls):
            raise _StateInvariantError("Moon deflection skip was not preserved")
    elif len([e for e in controls if e.kind == "deflector_selected"]) != 3:
        raise _StateInvariantError("Sun deflector selection is incomplete")
    returned_flag = _echo_request_bits(
        inner_flag | _implied_retflag_bits(normalized), raw_request
    )
    return _PrivateBodyState(
        body,
        raw_request,
        normalized,
        time.ut_bits,
        _bits(tt),
        time.seconds_bits,
        time.days_bits,
        time_tag,
        (words[0], words[1], words[2], words[3], words[4], words[5]),
        (
            _bits(words[0]),
            _bits(words[1]),
            _bits(words[2]),
            _bits(words[3]),
            _bits(words[4]),
            _bits(words[5]),
        ),
        inner_flag,
        returned_flag,
        tuple(reader.events),
        tuple(controls),
    )


def _facts(event: _ControlEvent) -> dict[str, tuple[str, bytes | int | str | bool]]:
    """Expose one immutable event's encoded facts for local consistency checks."""
    return {key: (kind, value) for key, kind, value in event.fields}


def _validate_body_state(
    state: _PrivateBodyState,
    body: int,
    raw_request: int,
    normalized: int,
    time: _DefaultTimeEvaluation,
    tt: float,
    time_tag: tuple[tuple[str, int, str], ...],
    child_tags: tuple[tuple[str, str, int, str], ...],
) -> None:
    """Reject inconsistent test-hook records before any pair is returned."""
    if (
        state.body != body
        or state.raw_request != raw_request
        or state.normalized_request != normalized
        or type(state.inner_flag) is not int
        or state.inner_flag != normalized
        or state.ut_bits != time.ut_bits
        or state.tt_bits != _bits(tt)
        or state.time_seconds_bits != time.seconds_bits
        or state.time_days_bits != time.days_bits
        or state.time_tag != time_tag
        or state.word_bits != tuple(_bits(word) for word in state.words)
        or state.returned_flag
        != _echo_request_bits(
            state.inner_flag | _implied_retflag_bits(normalized), raw_request
        )
        or not all(math.isfinite(word) for word in state.words)
    ):
        raise _StateInvariantError("private state words, flags, or time differ")
    if not state.sources or not state.controls:
        raise _StateInvariantError("private source/control history is empty")
    for event in state.sources:
        selected = [
            decision
            for decision in event.decisions
            if decision.tag == event.selected_tag and decision.covers
        ]
        if (
            event.scope_body != body
            or event.selected_tag not in child_tags
            or len(event.decisions) != 2
            or tuple(decision.tag for decision in event.decisions) != child_tags
            or len(selected) != 1
            or tuple(decision.tier for decision in event.decisions)
            != ("base", "medium")
            or len(event.value_bits) != (2 if event.kind == "nutation" else 6)
            or any(len(bits) != 8 for bits in event.value_bits)
        ):
            raise _StateInvariantError("private source event differs")
        if event.kind == "nutation":
            if event.source_body is not None or selected[0].channel is not None:
                raise _StateInvariantError("nutation source role differs")
        elif (
            event.kind not in ("body", "body_split")
            or event.source_body not in _BODY_ROLES
            or selected[0].channel != _BODY_ROLES[cast(int, event.source_body)]
        ):
            raise _StateInvariantError("body source role or channel differs")
    kinds = {event.kind for event in state.controls}
    if not {
        "pipeline_icrs",
        "light_time_iteration",
        "deflection_branch",
        "aberration_branch",
        "frame_request",
        "frame_selected",
    } <= kinds or any(kind.endswith("read_skipped") for kind in kinds):
        raise _StateInvariantError("private control event differs")
    if not any(event.kind == "nutation" for event in state.sources):
        raise _StateInvariantError("frame has no selected LEB nutation")
    for control in state.controls:
        facts = _facts(control)
        if control.kind == "pipeline_icrs" and facts != {
            "body": ("int", body),
            "system_bary": ("bool", False),
        }:
            raise _StateInvariantError("pipeline control differs")
        if control.kind == "frame_request" and facts != {
            "flags": ("int", normalized),
            "xyz": ("bool", True),
        }:
            raise _StateInvariantError("frame control differs")
        if control.kind == "deflection_branch" and facts.get("applied") != (
            "bool",
            body == SUN,
        ):
            raise _StateInvariantError("deflection branch differs")
        if control.kind == "aberration_branch" and facts.get("applied") != (
            "bool",
            True,
        ):
            raise _StateInvariantError("aberration branch differs")
    frames = [event for event in state.controls if event.kind == "frame_selected"]
    nutations = [event for event in state.sources if event.kind == "nutation"]
    if len(frames) != len(nutations):
        raise _StateInvariantError("frame/nutation count differs")
    for frame, nutation in zip(frames, nutations):
        facts = _facts(frame)
        if (
            facts.get("jd") != ("float", nutation.jd_bits)
            or facts.get("dpsi") != ("float", nutation.value_bits[0])
            or facts.get("deps") != ("float", nutation.value_bits[1])
        ):
            raise _StateInvariantError("frame/nutation words differ")
    iterations = [e for e in state.controls if e.kind == "light_time_iteration"]
    target_reads = [
        e
        for e in state.sources
        if e.kind in ("body", "body_split") and e.source_body == body
    ]
    retarded = target_reads[1 : 1 + len(iterations)]
    if len(retarded) != len(iterations):
        raise _StateInvariantError("retarded source count differs")
    for iteration, source in zip(iterations, retarded):
        facts = _facts(iteration)
        if facts.get("rounded_jd") != ("float", source.jd_bits) or facts.get(
            "residual"
        ) != ("float", source.offset_bits):
            raise _StateInvariantError("two-part source epoch differs")
    deflectors = [
        _facts(event) for event in state.controls if event.kind == "deflector_selected"
    ]
    if body == MOON:
        if deflectors:
            raise _StateInvariantError("Moon unexpectedly selected deflectors")
        return
    if [item.get("body") for item in deflectors] != [
        ("int", SUN),
        ("int", JUPITER),
        ("int", SATURN),
    ]:
        raise _StateInvariantError("Sun deflector roles differ")
    closest = [
        _facts(event)
        for event in state.controls
        if event.kind == "deflector_closest_approach"
    ]
    if len(closest) != 3:
        raise _StateInvariantError("Sun closest-approach controls differ")
    source_suffix = [
        event
        for event in state.sources
        if event.kind == "body" and event.source_body in (SUN, JUPITER, SATURN)
    ][-6:]
    if [event.source_body for event in source_suffix] != [
        SUN,
        SUN,
        JUPITER,
        JUPITER,
        SATURN,
        SATURN,
    ]:
        raise _StateInvariantError("Sun deflector reads differ")
    for index, facts in enumerate(closest):
        if facts.get("body") != deflectors[index]["body"] or facts.get("epoch") != (
            "float",
            source_suffix[2 * index + 1].jd_bits,
        ):
            raise _StateInvariantError("Sun closest-approach epoch differs")
    terminal_kinds = {
        "deflector_near_skip",
        "deflector_line_of_sight_skip",
        "deflector_limiter",
    }
    for deflector_facts in deflectors:
        terminal = [
            event
            for event in state.controls
            if event.kind in terminal_kinds
            and _facts(event).get("body") == deflector_facts["body"]
        ]
        if len(terminal) != 1:
            raise _StateInvariantError("Sun deflector outcome is incomplete")


def _validate_input(
    ut: float, request: int, paths: Mapping[str, str | os.PathLike[str]]
) -> dict[str, str]:
    """Enforce the contract's pre-open validation and error precedence."""
    if type(ut) is not float or type(request) is not int:
        raise TypeError("UT and state request require native float and int")
    locators = _paths_once(paths)
    if not math.isfinite(ut):
        raise ValueError("UT must be finite")
    if request not in _REQUESTS:
        raise _StateUnsupportedError("state request is outside the two admitted words")
    return locators


def _reduce(
    ut: float,
    request: int,
    paths: Mapping[str, str | os.PathLike[str]],
    order: tuple[int, ...],
    *,
    test_hook: Callable[[_ControlEvent], None] | None = None,
    result_hook: Callable[[_PrivateBodyState], _PrivateBodyState] | None = None,
) -> tuple[tuple[_PrivateBodyState, ...], _DefaultTimeEvaluation]:
    """Own and close every child while returning no partial result on failure."""
    locators = _validate_input(ut, request, paths)
    _check_runtime()
    owned: list[_LEB2SnapshotReader | _OwnedDefaultTime] = []
    primary: BaseException | None = None
    try:
        base = _open_production_snapshot("base", locators["base_core.leb2"])
        owned.append(base)
        medium = _open_production_snapshot("medium", locators["medium_core.leb2"])
        owned.append(medium)
        time_paths = {role: locators[role] for role in _ROLES[2:]}
        clock = _open_production_time(time_paths)
        owned.append(clock)
        time_tag = clock.production_asset_tag
        if time_tag is None:
            raise _StateConfigurationError("private time has no production admission")
        declaration = _DefaultTimeDeclaration()
        if declaration != _DefaultTimeDeclaration(
            None, False, "smh2016", False, True, False
        ):
            raise _StateConfigurationError("private default declaration differs")
        evaluation = clock.evaluate(ut, declaration)
        tt = ut + evaluation.days
        if not math.isfinite(tt):
            raise _StateUnsupportedError("derived TT is nonfinite")
        reader = _ObservedTierReader(base, medium)
        child_tags = (base.production_asset_tag, medium.production_asset_tag)
        if any(tag is None for tag in child_tags):
            raise _StateAssetIntegrityError("private LEB admission lost")
        normalized = request | FLG_SWIEPH
        results_list: list[_PrivateBodyState] = []
        for body in order:
            state = _run_body(
                reader,
                evaluation,
                tt,
                body,
                request,
                normalized,
                time_tag,
                test_hook=test_hook,
            )
            if result_hook is not None:
                original = deepcopy(state)
                state = result_hook(state)
                if state != original:
                    raise _StateInvariantError("private body record was changed")
            _validate_body_state(
                state,
                body,
                request,
                normalized,
                evaluation,
                tt,
                time_tag,
                cast(tuple[tuple[str, str, int, str], ...], child_tags),
            )
            results_list.append(state)
        if clock.production_asset_tag != time_tag:
            raise _StateClosedError("private time admission lost at finalization")
        for (tier, child), initial_tag in zip(reader.children, child_tags):
            tag = child.production_asset_tag
            if tag is None or tag[0] != tier or tag != initial_tag:
                raise _StateClosedError("private LEB admission lost at finalization")
        return tuple(results_list), evaluation
    except Exception as exc:
        translated = _translate(exc)
        primary = translated
        if translated is exc:
            raise
        raise translated from exc
    except BaseException as exc:
        primary = exc
        raise
    finally:
        close_failure: Exception | None = None
        for source in reversed(owned):
            try:
                source.close()
            except Exception as exc:
                if close_failure is None:
                    close_failure = exc
        if primary is None and close_failure is not None:
            raise _StateInvariantError(
                "private source cleanup failed"
            ) from close_failure


def _private_sun_moon_states(
    ut: float,
    request: int,
    paths: Mapping[str, str | os.PathLike[str]],
) -> _PrivateStatePair:
    """Return Moon then Sun diagnostics from five fixed owned assets."""
    states, shared = _reduce(ut, request, paths, (MOON, SUN))
    moon, sun = states
    return _PrivateStatePair(moon, sun, shared)


def _private_test_body_states(
    ut: float,
    request: int,
    paths: Mapping[str, str | os.PathLike[str]],
    order: tuple[int, int],
    *,
    test_hook: Callable[[_ControlEvent], None] | None = None,
    result_hook: Callable[[_PrivateBodyState], _PrivateBodyState] | None = None,
) -> tuple[_PrivateBodyState, _PrivateBodyState]:
    """Test-only body-order probe; never issues a state-source tag."""
    if type(order) is not tuple or len(order) != 2 or set(order) != {MOON, SUN}:
        raise TypeError("test order must contain Moon and Sun exactly once")
    states, _time = _reduce(
        ut, request, paths, order, test_hook=test_hook, result_hook=result_hook
    )
    first, second = states
    return first, second
