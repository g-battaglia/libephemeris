# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Focused, ephemeral regression checks for the private state diagnostic."""

from __future__ import annotations

import os
import shutil
import struct
from dataclasses import replace
from pathlib import Path
from threading import Event, Thread
from typing import Any, Generator, cast

import libephemeris as le
import pytest
import skyfield.data

from libephemeris import (
    _private_state_reducer as reducer,
    fast_calc,
    leb_compression,
    precession_vondrak,
    sidereal_longterm,
    state,
)
from libephemeris._leb2_snapshot import (
    _PRODUCTION_MANIFEST,
    _SnapshotCorruptionError,
    _SnapshotIntegrityError,
    _open_test_snapshot,
)
from libephemeris.constants import FLG_EQUATORIAL, FLG_SWIEPH, FLG_XYZ, MOON, SUN
from libephemeris.leb2_reader import LEB2Reader

DATES = (
    float.fromhex("0x1.2492b4020c49cp+21"),
    float.fromhex("0x1.26cd600000000p+21"),
    float.fromhex("0x1.2b42c80000000p+21"),
    float.fromhex("0x1.2d7d7c0000000p+21"),
)
RAW = FLG_EQUATORIAL | FLG_XYZ
REQUESTS = (RAW, RAW | FLG_SWIEPH)


def _paths() -> dict[str, Path]:
    """Locate only fixed project and installed Skyfield asset roles."""
    root = Path(__file__).resolve().parents[2]
    sf = Path(skyfield.data.__file__).parent
    return {
        "base_core.leb2": root / "data/leb2/base_core.leb2",
        "medium_core.leb2": root / "data/leb2/medium_core.leb2",
        "iers.npz": sf / "iers.npz",
        "historic_deltat.npy": sf / "historic_deltat.npy",
        "delta_t.npz": sf / "delta_t.npz",
    }


def _word_bits(values: tuple[float, ...]) -> tuple[bytes, ...]:
    """Compare native binary64 words without saving astronomical vectors."""
    return tuple(struct.pack("!d", value) for value in values)


@pytest.fixture
def controlled_ordinary(monkeypatch: pytest.MonkeyPatch) -> Generator[None, None, None]:
    """Select the stable ordinary sealed-LEB/default-time comparison branch."""
    state.close()
    monkeypatch.setattr(state, "_CALC_MODE", "leb")
    monkeypatch.setattr(state, "_PRECISION_TIER", "medium")
    monkeypatch.setattr(state, "_LEB_FILE", None)
    monkeypatch.setattr(state, "_EPHEMERIS_FILE_EXPLICIT", False)
    monkeypatch.setattr(state, "_IERS_DELTA_T_ENABLED", False)
    monkeypatch.setattr(state, "_DELTA_T_MODEL", "smh2016")
    monkeypatch.setattr(state, "_DELTA_T_USERDEF", None)
    monkeypatch.setattr(state, "_TIDAL_ACCELERATION", None)
    yield
    state.close()


@pytest.mark.parametrize("ut", DATES)
@pytest.mark.parametrize("state_request", REQUESTS)
def test_bitwise_pair_and_reverse_order(
    ut: float, state_request: int, controlled_ordinary: None
) -> None:
    """Both bodies, flags, time and selected nutation survive each order."""
    paths = _paths()
    pair = reducer._private_sun_moon_states(ut, state_request, paths)
    reverse = reducer._private_test_body_states(ut, state_request, paths, (SUN, MOON))
    assert pair.source_status == "diagnostic-only"
    assert (pair.moon.body, pair.sun.body) == (MOON, SUN)
    for private, reversed_state in ((pair.moon, reverse[1]), (pair.sun, reverse[0])):
        ordinary, flag = le.calc_ut(ut, private.body, state_request)
        assert all(type(word) is float for word in private.words)
        assert private.word_bits == _word_bits(private.words)
        assert private.word_bits == _word_bits(tuple(ordinary))
        assert private.returned_flag == flag
        assert private.word_bits == reversed_state.word_bits
        assert private.sources != () and reversed_state.sources != ()
        assert any(event.kind == "nutation" for event in private.sources)
        assert any(event.kind == "nutation" for event in reversed_state.sources)
        assert private.time_tag == pair.sun.time_tag == pair.moon.time_tag
        assert private.ut_bits == pair.time.ut_bits
        assert private.time_seconds_bits == pair.time.seconds_bits
        assert private.time_days_bits == pair.time.days_bits
        assert (
            len([e for e in private.controls if e.kind == "light_time_iteration"]) == 3
        )
        assert not any(e.kind.endswith("read_skipped") for e in private.controls)
    assert any(e.kind == "deflector_limiter" for e in pair.sun.controls)
    assert len([e for e in pair.sun.controls if e.kind == "deflector_selected"]) == 3
    assert not any(e.kind == "deflector_selected" for e in pair.moon.controls)
    if ut == DATES[0]:
        tiers = {e.selected_tag[0] for e in pair.sun.sources}
        assert tiers == {"base", "medium"}
        iterations = [e for e in pair.sun.controls if e.kind == "light_time_iteration"]
        target_reads = [e for e in pair.sun.sources if e.source_body == SUN]
        for control, source in zip(iterations, target_reads[1:]):
            facts = reducer._facts(control)
            assert facts["rounded_jd"] == ("float", source.jd_bits)
            assert facts["residual"] == ("float", source.offset_bits)
            assert len(source.decisions) == 2


def _cache_image() -> tuple[object, object, object]:
    """Capture global frame values and cache hit/miss counters."""
    return (
        fast_calc._leb_frame_cache.copy(),
        precession_vondrak._precession_matrix.cache_info(),
        sidereal_longterm.mean_obliquity_rad.cache_info(),
    )


@pytest.mark.parametrize("warm", [False, True])
def test_global_cache_and_thread_local_isolation(warm: bool) -> None:
    """Private calls neither consult nor update ordinary numerical caches."""
    fast_calc._leb_frame_cache.clear()
    precession_vondrak._precession_matrix.cache_clear()
    sidereal_longterm.mean_obliquity_rad.cache_clear()
    if warm:
        precession_vondrak._precession_matrix(DATES[2], True)
        sidereal_longterm.mean_obliquity_rad(DATES[2])
        fast_calc._leb_frame_cache[("sentinel",)] = ("value",)  # type: ignore[index,assignment]
    before = _cache_image()
    old_active = vars(fast_calc._active_local).copy()
    old_private = vars(fast_calc._private_state_local).copy()
    reducer._private_test_body_states(DATES[2], RAW, _paths(), (SUN, MOON))
    assert _cache_image() == before
    assert vars(fast_calc._active_local) == old_active
    assert vars(fast_calc._private_state_local) == old_private


def test_public_result_unchanged_across_private_call(
    controlled_ordinary: None,
) -> None:
    """A private diagnostic cannot change a controlled ordinary LEB result."""
    before_words, before_flag = le.calc_ut(DATES[2], SUN, RAW)
    reducer._private_sun_moon_states(DATES[2], RAW, _paths())
    after_words, after_flag = le.calc_ut(DATES[2], SUN, RAW)
    assert _word_bits(tuple(before_words)) == _word_bits(tuple(after_words))
    assert before_flag == after_flag


def test_state_close_during_private_frame_never_falls_back(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Ordinary generation invalidation cannot redirect the private frame."""
    ordinary_frame = fast_calc._get_skyfield_frame_data

    def forbidden(_jd: float) -> object:
        raise AssertionError("ordinary Skyfield frame reached")

    forbidden.cache_clear = ordinary_frame.cache_clear  # type: ignore[attr-defined]
    monkeypatch.setattr(fast_calc, "_get_skyfield_frame_data", forbidden)
    seen = False
    closed = Event()
    thread_errors: list[BaseException] = []

    def ordinary_close() -> None:
        try:
            state.close()
        except BaseException as exc:
            thread_errors.append(exc)
        finally:
            closed.set()

    def hook(event: reducer._ControlEvent) -> None:
        nonlocal seen
        if event.kind == "frame_request" and not seen:
            seen = True
            closer = Thread(target=ordinary_close)
            closer.start()
            assert closed.wait(5.0)
            closer.join()
            assert not thread_errors

    pair = reducer._private_test_body_states(
        DATES[2], RAW, _paths(), (MOON, SUN), test_hook=hook
    )
    assert seen
    assert all(any(e.kind == "nutation" for e in body.sources) for body in pair)
    assert fast_calc._get_skyfield_frame_data is forbidden
    assert ordinary_frame is not forbidden


def test_private_decoder_and_bitwise_coefficients(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Owned decoder retains shared transforms without the ordinary global decoder."""
    path = _paths()["base_core.leb2"]
    with LEB2Reader(str(path)) as ordinary:
        expected = ordinary._decompress_chunk(SUN, 0)
    original = leb_compression._DECOMPRESSOR

    class Forbidden:
        def decompress(self, *_args: object, **_kwargs: object) -> bytes:
            raise AssertionError("ordinary decoder reached")

    monkeypatch.setattr(leb_compression, "_DECOMPRESSOR", Forbidden())
    with reducer._open_production_snapshot("base", path) as private:
        assert private._decompress_chunk(SUN, 0) == expected
    reducer._private_sun_moon_states(DATES[2], RAW, _paths())
    assert leb_compression._DECOMPRESSOR is not original


@pytest.mark.parametrize(
    "value,state_request,error",
    [
        (True, RAW, TypeError),
        (1, RAW, TypeError),
        (float("nan"), RAW, ValueError),
        (float("inf"), RAW, ValueError),
        (DATES[2], True, TypeError),
        (DATES[2], RAW | 1, reducer._StateUnsupportedError),
    ],
)
def test_preopen_input_order(
    value: object,
    state_request: object,
    error: type[Exception],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Bad scalar/state_request inputs fail before any asset factory opens."""

    def forbidden(*_args: object) -> object:
        raise AssertionError("asset opened")

    monkeypatch.setattr(reducer, "_open_production_snapshot", forbidden)
    with pytest.raises(error):
        reducer._private_sun_moon_states(value, state_request, _paths())  # type: ignore[arg-type]


def test_bad_path_roles_and_type_precede_open(monkeypatch: pytest.MonkeyPatch) -> None:
    """Path map and locator validation is deterministic and pre-open."""
    monkeypatch.setattr(
        reducer, "_open_production_snapshot", lambda *_: pytest.fail("opened")
    )
    with pytest.raises(TypeError):
        reducer._private_sun_moon_states(DATES[2], RAW, {1: "file"})  # type: ignore[arg-type,dict-item]
    with pytest.raises(TypeError):
        reducer._private_sun_moon_states(DATES[2], RAW, {"x": 4})  # type: ignore[arg-type,dict-item]
    with pytest.raises(ValueError):
        reducer._private_sun_moon_states(DATES[2], RAW, {})


@pytest.mark.parametrize("role", tuple(_paths()))
@pytest.mark.parametrize("change", ["truncate", "append", "flip"])
def test_each_asset_is_hash_and_count_pinned(
    role: str,
    change: str,
    tmp_path: Path,
) -> None:
    """A changed fixed role cannot enter through a matching path label."""
    paths = _paths()
    copied = tmp_path / role
    shutil.copyfile(paths[role], copied)
    paths[role] = copied
    if change == "truncate":
        with copied.open("r+b") as stream:
            stream.truncate(copied.stat().st_size - 1)
    elif change == "append":
        with copied.open("ab") as stream:
            stream.write(b"x")
    else:
        with copied.open("r+b") as stream:
            byte = stream.read(1)
            stream.seek(0)
            stream.write(bytes((byte[0] ^ 1,)))
    with pytest.raises(reducer._StateAssetIntegrityError):
        reducer._private_sun_moon_states(DATES[2], RAW, paths)


def test_first_unavailable_role_blocks(tmp_path: Path) -> None:
    """The base role opens before other missing roles and reports unavailability."""
    paths = _paths()
    paths["base_core.leb2"] = tmp_path / "missing_base"
    paths["iers.npz"] = tmp_path / "missing_time"
    with pytest.raises(reducer._StateAssetUnavailableError) as caught:
        reducer._private_sun_moon_states(DATES[2], RAW, paths)
    assert "missing_base" in str(caught.value)


@pytest.mark.parametrize("change", ["replace", "in_place"])
def test_path_replacement_and_in_place_mutation_after_admission(
    change: str,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Later filesystem mutation cannot change the owned source graph or words."""
    paths = _paths()
    expected = reducer._private_sun_moon_states(DATES[2], RAW, paths)
    copied = tmp_path / "base_core.leb2"
    shutil.copyfile(paths["base_core.leb2"], copied)
    paths["base_core.leb2"] = copied
    real_run = reducer._run_body
    changed = False

    def mutate_then_run(*args: object, **kwargs: object) -> reducer._PrivateBodyState:
        nonlocal changed
        if not changed:
            changed = True
            if change == "replace":
                replacement = tmp_path / "replacement"
                replacement.write_bytes(b"not a LEB2 file")
                os.replace(replacement, copied)
            else:
                reader = cast(reducer._ObservedTierReader, args[0])
                tt = cast(float, args[2])
                selected, _decisions = reader._select(SUN, tt, 0.0)
                assert selected is reader.children[0][1]
                chunk = selected._chunk_index[SUN][selected._find_chunk(SUN, tt)]
                with copied.open("r+b") as stream:
                    stream.seek(chunk.blob_offset + 5)
                    byte = stream.read(1)
                    stream.seek(chunk.blob_offset + 5)
                    stream.write(bytes((byte[0] ^ 1,)))
        return real_run(*args, **kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr(reducer, "_run_body", mutate_then_run)
    actual = reducer._private_sun_moon_states(DATES[2], RAW, paths)
    assert changed
    assert (actual.moon.word_bits, actual.sun.word_bits) == (
        expected.moon.word_bits,
        expected.sun.word_bits,
    )
    assert (actual.moon.sources, actual.sun.sources) == (
        expected.moon.sources,
        expected.sun.sources,
    )
    with pytest.raises(reducer._StateAssetIntegrityError):
        reducer._private_sun_moon_states(DATES[2], RAW, paths)


def test_alternate_candidate_channel_blocks_before_reduction(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A second candidate cannot borrow the first child's channel metadata."""
    original = reducer._open_production_snapshot

    def changed(tier: str, path: Path) -> reducer._LEB2SnapshotReader:
        child = original(tier, path)
        if tier == "medium":
            child._bodies[SUN] = replace(child._bodies[SUN], coord_type=4)
        return child

    monkeypatch.setattr(reducer, "_open_production_snapshot", changed)
    with pytest.raises(reducer._StateUnsupportedError):
        reducer._private_sun_moon_states(DATES[2], RAW, _paths())


def test_test_only_child_cannot_gain_production_identity(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Real fixed bytes opened through the test factory remain unadmitted."""
    original = reducer._open_production_snapshot

    def test_first(tier: str, path: Path) -> reducer._LEB2SnapshotReader:
        if tier == "base":
            return _open_test_snapshot(path, _PRODUCTION_MANIFEST[tier])
        return original(tier, path)

    monkeypatch.setattr(reducer, "_open_production_snapshot", test_first)
    with pytest.raises(reducer._StateAssetIntegrityError):
        reducer._private_sun_moon_states(DATES[2], RAW, _paths())


@pytest.mark.parametrize("failure", ["missing", "corrupt"])
def test_missing_or_corrupt_nutation_blocks(
    failure: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Frame dispatch has no ordinary or third-party nutation fallback."""
    original = reducer._open_production_snapshot

    def changed(tier: str, path: Path) -> reducer._LEB2SnapshotReader:
        child = original(tier, path)
        if failure == "missing":
            child._nutation = None
        else:

            def corrupt(_jd: float) -> tuple[float, float]:
                raise _SnapshotCorruptionError("injected nutation corruption")

            monkeypatch.setattr(child, "eval_nutation", corrupt)
        return child

    monkeypatch.setattr(reducer, "_open_production_snapshot", changed)
    expected = (
        reducer._StateUnsupportedError
        if failure == "missing"
        else reducer._StateAssetCorruptionError
    )
    with pytest.raises(expected):
        reducer._private_sun_moon_states(DATES[2], RAW, _paths())


@pytest.mark.parametrize("read_number", [1, 2])
@pytest.mark.parametrize(
    "failure", [KeyError("missing"), _SnapshotCorruptionError("corrupt")]
)
def test_deflector_read_failure_is_terminal(
    read_number: int, failure: Exception, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Observation and closest-read errors cannot be swallowed as skips."""
    original = reducer._LEB2SnapshotReader.eval_body
    seen = 0

    def fail(self: reducer._LEB2SnapshotReader, body: int, jd: float):
        nonlocal seen
        if body == 5:
            seen += 1
            if seen == read_number:
                raise failure
        return original(self, body, jd)

    monkeypatch.setattr(reducer._LEB2SnapshotReader, "eval_body", fail)
    expected = (
        reducer._StateAssetCorruptionError
        if isinstance(failure, _SnapshotCorruptionError)
        else reducer._StateInvariantError
    )
    with pytest.raises(expected):
        reducer._private_sun_moon_states(DATES[2], RAW, _paths())
    assert seen >= read_number


@pytest.mark.parametrize("read_number", [1, 2])
@pytest.mark.parametrize(
    "failure,expected",
    [
        (
            _SnapshotCorruptionError("selection corruption"),
            reducer._StateAssetCorruptionError,
        ),
        (
            _SnapshotIntegrityError("selection integrity"),
            reducer._StateAssetIntegrityError,
        ),
    ],
)
def test_deflector_selection_failure_is_terminal(
    read_number: int,
    failure: Exception,
    expected: type[RuntimeError],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Selection failures before observation and closest reads cannot be skipped."""
    original = reducer._ObservedTierReader._select
    seen = 0

    def fail(self: reducer._ObservedTierReader, body: int, jd: float, offset: float):
        nonlocal seen
        if self.scope_body == SUN and body == 5:
            seen += 1
            if seen == read_number:
                raise failure
        return original(self, body, jd, offset)

    monkeypatch.setattr(reducer._ObservedTierReader, "_select", fail)
    with pytest.raises(expected):
        reducer._private_sun_moon_states(DATES[2], RAW, _paths())
    assert seen >= read_number


def test_deflector_postread_tag_failure_is_terminal(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A lost tag after a successful deflector read is a typed integrity block."""
    original = reducer._LEB2SnapshotReader.eval_body
    seen = False

    def lose_tag(self: reducer._LEB2SnapshotReader, body: int, jd: float):
        nonlocal seen
        result = original(self, body, jd)
        if body == 5 and not seen:
            seen = True
            self._admission = None
        return result

    monkeypatch.setattr(reducer._LEB2SnapshotReader, "eval_body", lose_tag)
    with pytest.raises(reducer._StateAssetIntegrityError):
        reducer._private_sun_moon_states(DATES[2], RAW, _paths())
    assert seen


def test_altered_inner_flag_blocks(monkeypatch: pytest.MonkeyPatch) -> None:
    """A shared reducer flag must equal the exact normalized request word."""
    original = fast_calc._fast_calc_core

    def changed(*args: object, **kwargs: object):
        words, flag = original(*args, **kwargs)  # type: ignore[arg-type]
        return words, flag ^ FLG_SWIEPH

    monkeypatch.setattr(fast_calc, "_fast_calc_core", changed)
    with pytest.raises(reducer._StateInvariantError, match="inner reducer flag"):
        reducer._private_sun_moon_states(DATES[2], RAW, _paths())


@pytest.mark.parametrize(
    "kind", ["light_time_early_exit", "deflection_zero_geocentric_vector"]
)
def test_private_zero_distance_branch_is_unsupported(
    kind: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Synthetic activation of either unresolved branch blocks the pair."""

    def synthetic(*_args: object):
        fast_calc._record_private_control(kind)
        return (0.0,) * 6, RAW | FLG_SWIEPH

    monkeypatch.setattr(fast_calc, "_fast_calc_core", synthetic)
    with pytest.raises(reducer._StateUnsupportedError):
        reducer._private_sun_moon_states(DATES[2], RAW, _paths())


def test_actual_synthetic_zero_distance_events() -> None:
    """Shared math emits both zero branches on synthetic finite vectors."""

    class StopAtBranch(Exception):
        pass

    class ZeroReader:
        def eval_body(self, _body: int, _jd: float):
            return (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)

    seen: list[str] = []
    previous_private = vars(fast_calc._private_state_local).copy()
    previous_active = vars(fast_calc._active_local).copy()

    def stop(kind: str, _fields: dict[str, object]) -> None:
        seen.append(kind)
        if kind in ("light_time_early_exit", "deflection_zero_geocentric_vector"):
            raise StopAtBranch

    try:
        fast_calc._private_state_local.control = stop
        with pytest.raises(StopAtBranch):
            fast_calc._pipeline_icrs(
                cast(Any, ZeroReader()), DATES[2], MOON, RAW, want_xyz=True
            )
        with pytest.raises(StopAtBranch):
            fast_calc._apply_gravitational_deflection(
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                DATES[2],
                0.1,
                cast(Any, ZeroReader()),
            )
    finally:
        reducer._restore_thread_slot(fast_calc._private_state_local, previous_private)
        reducer._restore_thread_slot(fast_calc._active_local, previous_active)
    assert "light_time_early_exit" in seen
    assert "deflection_zero_geocentric_vector" in seen


@pytest.mark.parametrize(
    "field",
    [
        "ut_bits",
        "tt_bits",
        "returned_flag",
        "word_bits",
        "time_tag",
        "sources",
        "controls",
        "raw_request",
    ],
)
def test_mutated_result_hook_blocks(field: str) -> None:
    """The pre-finalization test hook cannot alter any result dependency."""

    def change(body: reducer._PrivateBodyState) -> reducer._PrivateBodyState:
        if field == "word_bits":
            return replace(body, word_bits=(b"\x80" + b"\0" * 7,) + body.word_bits[1:])
        if field == "sources":
            first = replace(body.sources[0], selected_tag=("test", "x", 0, "x"))
            return replace(body, sources=(first,) + body.sources[1:])
        if field == "controls":
            control = replace(body.controls[0], kind="wrong")
            return replace(body, controls=(control,) + body.controls[1:])
        old = getattr(body, field)
        new = (
            b"wrong"
            if isinstance(old, bytes)
            else old + 1
            if isinstance(old, int)
            else ()
        )
        return replace(body, **cast(Any, {field: new}))

    with pytest.raises(reducer._StateInvariantError):
        reducer._private_test_body_states(
            DATES[2], RAW, _paths(), (MOON, SUN), result_hook=change
        )


def test_event_failure_after_moon_and_thread_cleanup() -> None:
    """An event hook exception during Sun returns no partial pair or context."""
    prior_active = vars(fast_calc._active_local).copy()
    prior_private = vars(fast_calc._private_state_local).copy()

    def fail(event: reducer._ControlEvent) -> None:
        if event.kind == "deflector_selected":
            raise ArithmeticError("injected")

    with pytest.raises(reducer._StateInvariantError):
        reducer._private_test_body_states(
            DATES[2], RAW, _paths(), (MOON, SUN), test_hook=fail
        )
    assert vars(fast_calc._active_local) == prior_active
    assert vars(fast_calc._private_state_local) == prior_private


def test_control_event_mutation_blocks() -> None:
    """A recorder hook cannot rewrite even an immutable event by force."""

    def change(event: reducer._ControlEvent) -> None:
        object.__setattr__(event, "kind", "wrong")

    with pytest.raises(reducer._StateInvariantError):
        reducer._private_test_body_states(
            DATES[2], RAW, _paths(), (MOON, SUN), test_hook=change
        )


def test_coherent_signed_zero_result_mutation_blocks() -> None:
    """Changing +0.0 to -0.0 is visible even when tuple equality says equal."""

    def change(body: reducer._PrivateBodyState) -> reducer._PrivateBodyState:
        assert body.word_bits[3] == struct.pack("!d", 0.0)
        words = body.words[:3] + (-0.0,) + body.words[4:]
        return replace(
            body,
            words=words,
            word_bits=cast(
                tuple[bytes, bytes, bytes, bytes, bytes, bytes], _word_bits(words)
            ),
        )

    with pytest.raises(reducer._StateInvariantError):
        reducer._private_test_body_states(
            DATES[2], RAW, _paths(), (MOON, SUN), result_hook=change
        )


def test_dependency_mismatch_precedes_asset_open(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A changed numerical runtime blocks before touching the five paths."""
    monkeypatch.setattr(reducer.zstd, "__version__", "changed")
    monkeypatch.setattr(
        reducer, "_open_production_snapshot", lambda *_: pytest.fail("opened")
    )
    with pytest.raises(reducer._StateDependencyError):
        reducer._private_sun_moon_states(DATES[2], RAW, _paths())


def test_pathlike_resolved_once_before_admission() -> None:
    """No later evaluation consults a caller-controlled locator again."""

    class Once(os.PathLike[str]):
        def __init__(self, path: Path) -> None:
            self.path = path
            self.calls = 0

        def __fspath__(self) -> str:
            self.calls += 1
            if self.calls != 1:
                raise AssertionError("path resolved more than once")
            return str(self.path)

    paths = {role: Once(path) for role, path in _paths().items()}
    reducer._private_sun_moon_states(DATES[2], RAW, paths)
    assert all(locator.calls == 1 for locator in paths.values())


def test_close_during_selected_read_is_typed_and_cleans_context(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Closing an owned child during Moon yields no partial state or leaked slot."""
    original = reducer._LEB2SnapshotReader.eval_body
    seen = False
    prior_active = vars(fast_calc._active_local).copy()
    prior_private = vars(fast_calc._private_state_local).copy()

    def close_then_read(self: reducer._LEB2SnapshotReader, body: int, jd: float):
        nonlocal seen
        if not seen:
            seen = True
            self.close()
        return original(self, body, jd)

    monkeypatch.setattr(reducer._LEB2SnapshotReader, "eval_body", close_then_read)
    with pytest.raises(reducer._StateClosedError):
        reducer._private_sun_moon_states(DATES[2], RAW, _paths())
    assert seen
    assert vars(fast_calc._active_local) == prior_active
    assert vars(fast_calc._private_state_local) == prior_private


def test_cleanup_failure_prevents_success(monkeypatch: pytest.MonkeyPatch) -> None:
    """A close failure cannot yield an otherwise successful pair."""
    original = reducer._OwnedDefaultTime.close

    def fail_after_close(self: reducer._OwnedDefaultTime) -> None:
        original(self)
        raise OSError("injected close failure")

    monkeypatch.setattr(reducer._OwnedDefaultTime, "close", fail_after_close)
    with pytest.raises(reducer._StateInvariantError, match="cleanup failed"):
        reducer._private_sun_moon_states(DATES[2], RAW, _paths())


@pytest.mark.parametrize(
    "deflector,expected",
    [
        ((0.0, 0.0, 0.0), "deflector_near_skip"),
        ((-1.0, 0.0, 0.0), "deflector_line_of_sight_skip"),
        ((0.5, 0.0001, 0.0), "deflector_limiter"),
        ((0.5, 0.5, 0.0), "deflector_limiter"),
    ],
)
def test_synthetic_deflector_math_branches(
    deflector: tuple[float, float, float],
    expected: str,
) -> None:
    """Synthetic branch probes observe skips/limiter without a production tag."""

    class Reader:
        def eval_body(self, _body: int, _jd: float):
            return deflector, (0.0, 0.0, 0.0)

    events: list[tuple[str, dict[str, object]]] = []
    previous = vars(fast_calc._private_state_local).copy()
    try:
        fast_calc._private_state_local.control = lambda kind, fields: events.append(
            (kind, fields)
        )
        fast_calc._apply_gravitational_deflection(
            (1.0, 0.0, 0.0), (0.0, 0.0, 0.0), DATES[2], 0.1, Reader()
        )
    finally:
        reducer._restore_thread_slot(fast_calc._private_state_local, previous)
    assert any(kind == expected for kind, _fields in events)
    if expected == "deflector_limiter":
        limited = deflector[1] == 0.0001
        assert any(
            kind == expected and fields["limited"] is limited for kind, fields in events
        )
