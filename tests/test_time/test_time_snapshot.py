# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Falsification and ordinary-byte parity for the private time snapshot."""

from __future__ import annotations

import builtins
import hashlib
import io
import math
import os
import struct
import warnings
import zipfile
from dataclasses import replace
from pathlib import Path
from typing import Any, cast

import numpy as np
import pytest
import skyfield
import skyfield.data

import libephemeris._time_snapshot as snapshot
from libephemeris._time_snapshot import (
    _AssetPin,
    _DefaultTimeDeclaration,
    _TimeShape,
    _TimeSnapshotClosedError,
    _TimeSnapshotConfigurationError,
    _TimeSnapshotConstructionError,
    _TimeSnapshotDependencyError,
    _TimeSnapshotIntegrityError,
    _TimeSnapshotParseError,
    _TimeSnapshotUnsupportedError,
    _open_production_time,
    _open_test_time,
)
from libephemeris.state import _build_enhanced_timescale

_ROLES = ("iers.npz", "historic_deltat.npy", "delta_t.npz")
_SHAPE = _TimeShape(3, 2, 3, 3)
_DEFAULT = _DefaultTimeDeclaration()


def _bits(word: float) -> bytes:
    assert type(word) is float
    return struct.pack("!d", word)


def _npy(value: np.ndarray) -> bytes:
    buffer = io.BytesIO()
    np.save(buffer, value)
    return buffer.getvalue()


def _npz(**arrays: np.ndarray) -> bytes:
    buffer = io.BytesIO()
    cast(Any, np.savez)(buffer, **arrays)
    return buffer.getvalue()


def _synthetic_bytes() -> dict[str, bytes]:
    iers = {
        "tt_jd_minus_arange": np.array([2450000.0] * 3, dtype=np.float64),
        "delta_t_1e7": np.array([600e6, 610e6, 620e6], dtype=np.float64),
        "leap_dates": np.array([2440000.0, 2441000.0], dtype=np.float64),
        "leap_offsets": np.array([10.0, 11.0], dtype=np.float64),
    }
    historic = np.array(
        [[2390000.0, 2400000.0, 2450000.5], [10.0, 20.0, 30.0]],
        dtype=np.float64,
    )
    s15 = np.array(
        [
            [1800.0, 1900.0, 2000.0],
            [1900.0, 2000.0, 2100.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [10.0, 20.0, 30.0],
        ],
        dtype=np.float64,
    )
    return {
        "iers.npz": _npz(**iers),
        "historic_deltat.npy": _npy(historic),
        "delta_t.npz": _npz(**{"Table-S15.2020.txt": s15}),
    }


def _installed_paths() -> dict[str, Path]:
    data_dir = Path(skyfield.data.__file__).parent
    return {role: data_dir / role for role in _ROLES}


def _write_assets(
    tmp_path: Path, data: dict[str, bytes]
) -> tuple[dict[str, Path], dict[str, _AssetPin]]:
    paths: dict[str, Path] = {}
    pins: dict[str, _AssetPin] = {}
    for role in _ROLES:
        path = tmp_path / role
        path.write_bytes(data[role])
        paths[role] = path
        pins[role] = _AssetPin(
            role, len(data[role]), hashlib.sha256(data[role]).hexdigest()
        )
    return paths, pins


def _synthetic(
    tmp_path: Path,
) -> tuple[dict[str, Path], dict[str, _AssetPin]]:
    return _write_assets(tmp_path, _synthetic_bytes())


def _alter_npz(data: dict[str, bytes], role: str, **edits: np.ndarray) -> None:
    with np.load(io.BytesIO(data[role]), allow_pickle=False) as archive:
        arrays = {key: archive[key].copy() for key in archive.files}
    arrays.update(edits)
    data[role] = _npz(**arrays)


def _duplicate_npz_member(data: bytes, name: str) -> bytes:
    """Keep a matching test pin while making an archive member ambiguous."""
    buffer = io.BytesIO()
    with (
        zipfile.ZipFile(io.BytesIO(data)) as source,
        zipfile.ZipFile(buffer, "w") as target,
    ):
        for member in source.namelist():
            target.writestr(member, source.read(member))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            target.writestr(f"{name}.npy", source.read(f"{name}.npy"))
    return buffer.getvalue()


def test_production_bitwise_seconds_days_leaps_and_cache_isolation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    paths = _installed_paths()
    ordinary = _build_enhanced_timescale()
    with _open_production_time(paths) as private:
        assert private.production_asset_tag == tuple(
            (
                role,
                snapshot._PRODUCTION_PINS[role].byte_count,
                snapshot._PRODUCTION_PINS[role].sha256,
            )
            for role in _ROLES
        )
        assert private._test_leap_words() == (
            ordinary.leap_dates.tobytes(),
            ordinary.leap_offsets.tobytes(),
        )
        for date in (
            (2016, 12, 31, 23, 59, 59.0),
            (2016, 12, 31, 23, 59, 60.0),
            (2017, 1, 1, 0, 0, 0.0),
        ):
            t = ordinary.utc(*date)
            assert private._test_utc_words(*date) == (
                _bits(float(t.tt)),
                _bits(float(t.tai)),
            )
        assert ordinary.leap_dates[-1] == 2457754.5

        with np.load(paths["iers.npz"], allow_pickle=False) as archive:
            iers_first = float(archive["tt_jd_minus_arange"][0])
        with np.load(paths["delta_t.npz"], allow_pickle=False) as archive:
            s15 = archive["Table-S15.2020.txt"]
        historic_first = float(ordinary.delta_t_function.table_tt[0])
        iers_last = float(ordinary.delta_t_function.table_tt[-1])
        s15_first = float(s15[0, 0])
        spline_cut = (historic_first - 1721045.0) / 365.25
        retained = s15[:, : np.searchsorted(s15[0], spline_cut)]
        spline_midpoints = tuple(
            1721045.0 + float((lower + upper) / 2.0) * 365.25
            for lower, upper in zip(retained[0], retained[1])
        )
        right_join = (iers_last - 1721045.0) / 365.25
        right_patch_end = (right_join + 800.0) // 100.0 * 100.0
        coverage_neighbors = tuple(
            float(ordinary.ut1(year, 1, 1).ut1)
            for year in (1549, 1551, 1848, 1850, 2149, 2151, 2649, 2651)
        )
        boundary_dates = tuple(
            1721045.0 + year * 365.25
            for year in (
                s15_first - 900.0,
                s15_first - 800.0,
                s15_first,
                spline_cut,
                right_patch_end,
                right_patch_end + 100.0,
            )
        )
        joins = (historic_first, iers_first, iers_last, *boundary_dates)
        exact_tt_boundaries = (
            *joins,
            *(
                1721045.0 + float(year) * 365.25
                for year in (*retained[0], *retained[1])
            ),
        )
        for boundary in exact_tt_boundaries:
            for tt in (
                math.nextafter(boundary, -math.inf),
                boundary,
                math.nextafter(boundary, math.inf),
            ):
                assert private._test_tt_seconds_bits(tt) == _bits(
                    float(ordinary.delta_t_function(tt))
                )
        dates = (
            -1000000.0,
            500000.0,
            2100000.0,
            2300000.0,
            2390000.0,
            2451545.0,
            2460000.0,
            2600000.0,
            3000000.0,
            5000000.0,
            *coverage_neighbors,
            *spline_midpoints,
            *(
                jd
                for join in joins
                for jd in (
                    math.nextafter(join, -math.inf),
                    join,
                    math.nextafter(join, math.inf),
                )
            ),
        )
        for ut in (*dates, 2451545.0):
            t = ordinary.ut1_jd(ut)
            seconds = float(t.delta_t)
            days = (seconds + 0.0) / 86400.0
            got = private.evaluate(float(ut), _DEFAULT)
            assert type(got.seconds) is type(got.days) is type(got.ut) is float
            assert got.ut_bits == _bits(float(ut))
            assert got.seconds_bits == _bits(seconds)
            assert got.days_bits == _bits(days)

        def no_global(*_args: object, **_kwargs: object) -> None:
            raise AssertionError("private evaluator touched ordinary time cache")

        monkeypatch.setattr("libephemeris.cache.get_cached_time_ut1", no_global)
        monkeypatch.setattr("libephemeris.state.get_timescale", no_global)
        assert private.evaluate(2451545.0, _DEFAULT).days_bits == _bits(
            (float(ordinary.ut1_jd(2451545.0).delta_t) + 0.0) / 86400.0
        )


def test_private_builder_never_calls_ambient_skyfield_loader(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    paths, pins = _synthetic(tmp_path)

    def forbidden(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("ambient Skyfield bundled-data loader called")

    monkeypatch.setattr(snapshot._sf_time, "load_bundled_npy", forbidden)
    monkeypatch.setattr(snapshot._sf_functions, "load_bundled_npy", forbidden)
    with _open_test_time(paths, pins, _SHAPE) as private:
        assert private.production_asset_tag is None
        assert type(private.evaluate(2450000.0, _DEFAULT).days) is float


def test_real_asset_pins_do_not_promote_test_factory() -> None:
    with _open_test_time(
        _installed_paths(), snapshot._PRODUCTION_PINS, snapshot._PRODUCTION_SHAPE
    ) as private:
        assert private.production_asset_tag is None
        assert type(private.evaluate(2451545.0, _DEFAULT).days) is float


@pytest.mark.parametrize(
    "fault", ["short", "extra", "changed", "digest", "count", "role", "path_swap"]
)
def test_exact_byte_and_role_pins_reject(tmp_path: Path, fault: str) -> None:
    paths, pins = _synthetic(tmp_path)
    if fault in ("short", "extra", "changed"):
        path = paths["iers.npz"]
        data = bytearray(path.read_bytes())
        if fault == "short":
            data = data[:-1]
        elif fault == "extra":
            data += b"X"
        else:
            data[0] ^= 1
        path.write_bytes(data)
    elif fault == "digest":
        pins["iers.npz"] = replace(pins["iers.npz"], sha256="0" * 64)
    elif fault == "count":
        pins["iers.npz"] = replace(pins["iers.npz"], byte_count=1)
    elif fault == "path_swap":
        paths["iers.npz"], paths["delta_t.npz"] = (
            paths["delta_t.npz"],
            paths["iers.npz"],
        )
    else:
        pins["iers.npz"] = replace(pins["iers.npz"], role="delta_t.npz")
    with pytest.raises(_TimeSnapshotIntegrityError):
        _open_test_time(paths, pins, _SHAPE)


@pytest.mark.parametrize(
    "fault",
    [
        "zip",
        "npy",
        "missing",
        "duplicate_iers",
        "duplicate_spline",
        "shape",
        "nonfinite",
        "unordered_tt",
        "unordered_historic",
        "bad_leaps",
        "bad_spline",
        "pickle_extra",
        "pickle_required",
    ],
)
def test_matching_synthetic_pin_reaches_typed_schema_rejection(
    tmp_path: Path, fault: str
) -> None:
    data = _synthetic_bytes()
    if fault == "zip":
        data["iers.npz"] = b"bad ZIP"
    elif fault == "npy":
        data["historic_deltat.npy"] = b"bad NPY"
    elif fault == "missing":
        with np.load(io.BytesIO(data["iers.npz"])) as archive:
            arrays = {key: archive[key] for key in archive.files if key != "leap_dates"}
        data["iers.npz"] = _npz(**arrays)
    elif fault == "duplicate_iers":
        data["iers.npz"] = _duplicate_npz_member(data["iers.npz"], "tt_jd_minus_arange")
    elif fault == "duplicate_spline":
        data["delta_t.npz"] = _duplicate_npz_member(
            data["delta_t.npz"], "Table-S15.2020.txt"
        )
    elif fault == "shape":
        _alter_npz(data, "iers.npz", tt_jd_minus_arange=np.array([2450000.0]))
    elif fault == "nonfinite":
        _alter_npz(data, "iers.npz", delta_t_1e7=np.array([1.0, np.nan, 2.0]))
    elif fault == "unordered_tt":
        _alter_npz(
            data,
            "iers.npz",
            tt_jd_minus_arange=np.array([2450000.0, 2449990.0, 2449980.0]),
        )
    elif fault == "unordered_historic":
        hist = np.load(io.BytesIO(data["historic_deltat.npy"]))
        hist[0, 1] = hist[0, 0]
        data["historic_deltat.npy"] = _npy(hist)
    elif fault == "bad_leaps":
        _alter_npz(data, "iers.npz", leap_offsets=np.array([10.0, 10.5]))
    elif fault == "bad_spline":
        with np.load(io.BytesIO(data["delta_t.npz"])) as archive:
            s15 = archive["Table-S15.2020.txt"].copy()
        s15[0, 1] += 1.0
        data["delta_t.npz"] = _npz(**{"Table-S15.2020.txt": s15})
    elif fault == "pickle_extra":
        _alter_npz(data, "iers.npz", evil=np.array([object()], dtype=object))
    else:
        _alter_npz(data, "iers.npz", delta_t_1e7=np.array([object()] * 3))
    paths, pins = _write_assets(tmp_path, data)
    with pytest.raises(_TimeSnapshotParseError):
        _open_test_time(paths, pins, _SHAPE)


def test_validated_arrays_can_still_fail_typed_construction(tmp_path: Path) -> None:
    data = _synthetic_bytes()
    hist = np.load(io.BytesIO(data["historic_deltat.npy"]))
    hist[0] = [2100000.0, 2150000.0, 2200000.0]
    data["historic_deltat.npy"] = _npy(hist)
    paths, pins = _write_assets(tmp_path, data)
    with pytest.raises(_TimeSnapshotConstructionError):
        _open_test_time(paths, pins, _SHAPE)


@pytest.mark.parametrize("ut", [True, 2451545, np.float64(2451545.0), "2451545"])
def test_wrong_ut_type_is_distinct(tmp_path: Path, ut: object) -> None:
    paths, pins = _synthetic(tmp_path)
    with _open_test_time(paths, pins, _SHAPE) as private:
        with pytest.raises(TypeError):
            private.evaluate(ut, _DEFAULT)  # type: ignore[arg-type]


@pytest.mark.parametrize("ut", [math.nan, math.inf, -math.inf])
def test_nonfinite_ut_is_invalid(tmp_path: Path, ut: float) -> None:
    paths, pins = _synthetic(tmp_path)
    with _open_test_time(paths, pins, _SHAPE) as private:
        with pytest.raises(ValueError, match="must be finite"):
            private.evaluate(ut, _DEFAULT)


def test_extreme_finite_ut_is_unsupported(tmp_path: Path) -> None:
    paths, pins = _synthetic(tmp_path)
    with _open_test_time(paths, pins, _SHAPE) as private:
        for ut in (1e308, -1e308, math.nextafter(1e12, math.inf)):
            with pytest.raises(_TimeSnapshotUnsupportedError):
                private.evaluate(ut, _DEFAULT)


def test_ut_support_guard_includes_its_finite_endpoints(tmp_path: Path) -> None:
    paths, pins = _synthetic(tmp_path)
    with _open_test_time(paths, pins, _SHAPE) as private:
        for ut in (-1e12, 1e12):
            result = private.evaluate(ut, _DEFAULT)
            assert result.ut_bits == _bits(ut)
            assert math.isfinite(result.seconds) and math.isfinite(result.days)


@pytest.mark.parametrize(
    "change",
    [
        {"user_delta_t": 0.0},
        {"iers_enabled": True},
        {"model": "espenak_meeus"},
        {"model": np.array(["smh2016", "smh2016"])},
        {"model_override": True},
        {"tidal_automatic": False},
        {"ordinary_fallback": True},
    ],
)
def test_nondefault_declarations_are_rejected(
    tmp_path: Path, change: dict[str, object]
) -> None:
    paths, pins = _synthetic(tmp_path)
    with _open_test_time(paths, pins, _SHAPE) as private:
        with pytest.raises(_TimeSnapshotConfigurationError):
            private.evaluate(2451545.0, replace(_DEFAULT, **cast(Any, change)))


@pytest.mark.parametrize("fault", ["skyfield", "numpy", "source"])
def test_dependency_mismatch_is_typed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, fault: str
) -> None:
    paths, pins = _synthetic(tmp_path)
    if fault == "skyfield":
        monkeypatch.setattr(skyfield, "__version__", "1.55")
    elif fault == "numpy":
        monkeypatch.setattr(np, "__version__", "2.4.0")
    else:
        monkeypatch.setattr(snapshot, "_source_digest", lambda _path: "0" * 64)
    with pytest.raises(_TimeSnapshotDependencyError):
        _open_test_time(paths, pins, _SHAPE)


def test_success_and_parse_failure_close_descriptors_before_parsing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    paths, pins = _synthetic(tmp_path)
    original_open = builtins.open
    opened: list[Any] = []

    def tracked(path: object, mode: str = "r", *args: object, **kwargs: object) -> Any:
        handle = cast(Any, original_open)(path, mode, *args, **kwargs)
        if path in paths.values():
            opened.append(handle)
        return handle

    monkeypatch.setattr(builtins, "open", tracked)
    original_parse = snapshot._parse_arrays

    def checked_parse(
        owned: tuple[bytes, ...], shape: _TimeShape
    ) -> tuple[np.ndarray, ...]:
        assert len(opened) == 3 and all(handle.closed for handle in opened)
        return original_parse(owned, shape)

    monkeypatch.setattr(snapshot, "_parse_arrays", checked_parse)
    with _open_test_time(paths, pins, _SHAPE):
        pass
    assert all(handle.closed for handle in opened)

    opened.clear()
    paths["iers.npz"].write_bytes(b"bad ZIP")
    bad = paths["iers.npz"].read_bytes()
    pins["iers.npz"] = _AssetPin("iers.npz", len(bad), hashlib.sha256(bad).hexdigest())
    with pytest.raises(_TimeSnapshotParseError):
        _open_test_time(paths, pins, _SHAPE)
    assert len(opened) == 3 and all(handle.closed for handle in opened)

    opened.clear()
    paths, pins = _synthetic(tmp_path)
    pins["iers.npz"] = replace(pins["iers.npz"], sha256="0" * 64)
    with pytest.raises(_TimeSnapshotIntegrityError):
        _open_test_time(paths, pins, _SHAPE)
    assert len(opened) == 1 and opened[0].closed


def test_read_error_and_midcopy_mutation_close_descriptor(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    paths, pins = _synthetic(tmp_path)
    original_open = builtins.open
    tracked_handles: list[Any] = []
    monkeypatch.setattr(snapshot, "_COPY_CHUNK_SIZE", 16)

    class ChangingReader:
        def __init__(self, inner: Any, *, fail: bool) -> None:
            self.inner = inner
            self.fail = fail
            self.reads = 0
            self.closed = False

        def __enter__(self) -> ChangingReader:
            return self

        def __exit__(self, *_args: object) -> None:
            self.inner.close()
            self.closed = True

        def read(self, count: int) -> bytes:
            self.reads += 1
            if self.reads == 2 and self.fail:
                raise OSError("injected read failure")
            part: bytes = self.inner.read(count)
            if self.reads == 1 and not self.fail:
                data = bytearray(paths["iers.npz"].read_bytes())
                data[20] ^= 1
                paths["iers.npz"].write_bytes(data)
            if self.reads == 2 and not self.fail:
                changed = bytearray(part)
                changed[4] ^= 1
                return bytes(changed)
            return part

    fail = True

    def intercept(
        path: object, mode: str = "r", *args: object, **kwargs: object
    ) -> Any:
        inner = cast(Any, original_open)(path, mode, *args, **kwargs)
        if path == paths["iers.npz"]:
            wrapper = ChangingReader(inner, fail=fail)
            tracked_handles.append(wrapper)
            return wrapper
        return inner

    monkeypatch.setattr(builtins, "open", intercept)
    with pytest.raises(_TimeSnapshotIntegrityError):
        _open_test_time(paths, pins, _SHAPE)
    assert tracked_handles[-1].closed

    fail = False
    with pytest.raises(_TimeSnapshotIntegrityError):
        _open_test_time(paths, pins, _SHAPE)
    assert tracked_handles[-1].closed


def test_live_bytes_ignore_replacement_and_inplace_change(tmp_path: Path) -> None:
    paths, pins = _synthetic(tmp_path)
    private = _open_test_time(paths, pins, _SHAPE)
    before = private.evaluate(2450000.0, _DEFAULT).days_bits
    replacement = tmp_path / "new-iers.npz"
    replacement.write_bytes(b"replacement")
    os.replace(replacement, paths["iers.npz"])
    assert private.evaluate(2450000.0, _DEFAULT).days_bits == before
    with pytest.raises(_TimeSnapshotIntegrityError):
        _open_test_time(paths, pins, _SHAPE)
    private.close()

    paths, pins = _synthetic(tmp_path)
    private = _open_test_time(paths, pins, _SHAPE)
    before = private.evaluate(2450000.0, _DEFAULT).seconds_bits
    changed = bytearray(paths["delta_t.npz"].read_bytes())
    changed[-1] ^= 1
    paths["delta_t.npz"].write_bytes(changed)
    assert private.evaluate(2450000.0, _DEFAULT).seconds_bits == before
    with pytest.raises(_TimeSnapshotIntegrityError):
        _open_test_time(paths, pins, _SHAPE)
    private.close()


def test_admission_is_immutable_and_close_releases_references(tmp_path: Path) -> None:
    paths, pins = _synthetic(tmp_path)
    with pytest.raises(TypeError):
        snapshot._OwnedDefaultTime()
    private = _open_test_time(paths, pins, _SHAPE)
    assert private.production_asset_tag is None
    for name, value in (
        ("_production", True),
        ("_tier", "production"),
        ("_role", "iers.npz"),
        ("_byte_count", 62966),
        ("_digest", snapshot._PRODUCTION_PINS["iers.npz"].sha256),
    ):
        with pytest.raises(AttributeError):
            setattr(private, name, value)
        assert private.production_asset_tag is None
    with pytest.raises(TypeError):
        snapshot._PRODUCTION_PINS["iers.npz"] = pins["iers.npz"]  # type: ignore[index]
    assert not hasattr(private, "__dict__")
    assert private._arrays is not None
    assert all(not array.flags.writeable for array in private._arrays)
    admission = snapshot._ADMISSIONS[private]
    assert admission.owned is private._owned
    direct = object.__new__(snapshot._OwnedDefaultTime)
    direct._owned = private._owned
    direct._timescale = private._timescale
    direct._arrays = private._arrays
    with pytest.raises(_TimeSnapshotIntegrityError):
        direct.evaluate(2450000.0, _DEFAULT)
    direct.close()
    timescale = private._timescale
    assert timescale is not None
    assert all(
        not array.flags.writeable
        for array in (
            timescale._leap_utc,
            timescale._leap_offsets,
            timescale._leap_tai,
            timescale.delta_t_function.long_term_function._n,
        )
    )
    private.close()
    private.close()
    assert private not in snapshot._ADMISSIONS
    assert private._owned is private._arrays is private._timescale is None
    with pytest.raises(_TimeSnapshotClosedError):
        private.production_asset_tag  # noqa: B018
    with pytest.raises(_TimeSnapshotClosedError):
        private.evaluate(2450000.0, _DEFAULT)
    with pytest.raises(_TimeSnapshotClosedError):
        private._test_leap_words()
    with pytest.raises(_TimeSnapshotClosedError):
        private._test_tt_seconds_bits(2450000.0)
    with pytest.raises(_TimeSnapshotClosedError):
        private._test_utc_words(2017, 1, 1)


def test_production_tag_detects_evaluator_swap() -> None:
    private = _open_production_time(_installed_paths())
    owned_timescale = private._timescale
    private._timescale = object()  # type: ignore[assignment]
    with pytest.raises(_TimeSnapshotIntegrityError):
        private.production_asset_tag  # noqa: B018
    private._timescale = owned_timescale
    assert private.production_asset_tag is not None
    private.close()
