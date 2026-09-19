# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Private three-asset, byte-owned default Delta-T construction.

Provenance:
    Project-authored adapter for the reviewed
    ``default-time-three-asset-snapshot-contract.md`` (SHA-256
    bc161b1f48c832afd3a6c464ab668b1aa338487f83f8adb22e590f62274846bd).
    Its numerical operation order follows the permitted MIT Skyfield 1.54
    ``build_delta_t`` path, using Skyfield's spline and DeltaT primitives.
    This disconnected component neither selects nor attests an ordinary call.
"""

from __future__ import annotations

import hashlib
import io
import math
import os
import struct
import weakref
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Any, Mapping, cast

import numpy as np
import skyfield
import skyfield.curvelib as _sf_curve
import skyfield.functions as _sf_functions
import skyfield.timelib as _sf_time
from skyfield.timelib import Timescale

_ROLES = ("iers.npz", "historic_deltat.npy", "delta_t.npz")
_COPY_CHUNK_SIZE = 65536


@dataclass(frozen=True, slots=True)
class _AssetPin:
    """An independently supplied byte identity for one logical asset."""

    role: str
    byte_count: int
    sha256: str


@dataclass(frozen=True, slots=True)
class _TimeShape:
    """Exact column counts for the three versioned asset schemas."""

    daily: int
    leaps: int
    historic: int
    splines: int


_PRODUCTION_PINS = MappingProxyType(
    {
        "iers.npz": _AssetPin(
            "iers.npz",
            62_966,
            "c7d7536d898dfa9f8cd43e8044ff51e108cc8289675a13fee9822010a1c4935c",
        ),
        "historic_deltat.npy": _AssetPin(
            "historic_deltat.npy",
            10_576,
            "f5346b780b36a0325b1847dc6c0083d66edc7e88b7f648b4c98a67bbd02b5d3f",
        ),
        "delta_t.npz": _AssetPin(
            "delta_t.npz",
            1_547,
            "2d12bd3e789543b78a1f53c8b76ed7fecffdf7e5149cfb6a0aed21a8b3db5ff6",
        ),
    }
)
_PRODUCTION_SHAPE = _TimeShape(19745, 27, 656, 58)
_SOURCE_DIGESTS = MappingProxyType(
    {
        "timelib": "a949b768909835ce3b504e89321241fbb5647065cfb17b0db225a285896109fd",
        "curvelib": "7074e301829fc41065471a1ca6dd5dfe2d05a6891a9358896fe7366fe5156964",
        "functions": "671c0316723482371437d69c39499cbbb0dbabc65e9ca73c6bea7510a694d9c1",
    }
)


class _TimeSnapshotIntegrityError(ValueError):
    """A role, exact count, or digest fails independent admission."""


class _TimeSnapshotDependencyError(RuntimeError):
    """The installed numerical dependency differs from the reviewed code."""


class _TimeSnapshotParseError(ValueError):
    """Accepted bytes do not satisfy the fixed array schemas."""


class _TimeSnapshotConstructionError(ValueError):
    """Validated arrays cannot form the reviewed numerical curve."""


class _TimeSnapshotConfigurationError(ValueError):
    """A declaration selects a branch outside default SMH2016."""


class _TimeSnapshotUnsupportedError(ValueError):
    """A finite UT date is outside the supported numerical domain."""


class _TimeSnapshotClosedError(ValueError):
    """An operation needs a live admitted private evaluator."""


@dataclass(frozen=True, slots=True)
class _DefaultTimeDeclaration:
    """Caller assertion of the default branch; not evidence of an ordinary call."""

    user_delta_t: float | None = None
    iers_enabled: bool = False
    model: str = "smh2016"
    model_override: bool = False
    tidal_automatic: bool = True
    ordinary_fallback: bool = False


@dataclass(frozen=True, slots=True)
class _DefaultTimeEvaluation:
    """Native results and the exact binary64 input/output words."""

    ut: float
    seconds: float
    days: float
    ut_bits: bytes
    seconds_bits: bytes
    days_bits: bytes


@dataclass(frozen=True, slots=True)
class _Admission:
    """Factory-owned connection from all accepted bytes to one Timescale."""

    pins: tuple[_AssetPin, ...]
    owned: tuple[bytes, ...]
    timescale: Timescale


@dataclass(frozen=True, slots=True)
class _ProductionAdmission(_Admission):
    """Capability issued solely by the fixed production factory."""


_ADMISSIONS: weakref.WeakKeyDictionary[_OwnedDefaultTime, _Admission] = (
    weakref.WeakKeyDictionary()
)


def _source_digest(path: Path) -> str:
    """Hash one installed, directly reviewed Skyfield source file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _check_dependencies() -> None:
    """Fail closed if numerical code or binary package versions have changed."""
    if skyfield.__version__ != "1.54" or np.__version__ != "2.3.5":
        raise _TimeSnapshotDependencyError("Skyfield/NumPy version differs")
    for name, module in (
        ("timelib", _sf_time),
        ("curvelib", _sf_curve),
        ("functions", _sf_functions),
    ):
        try:
            path = Path(module.__file__ or "")
            actual = _source_digest(path)
        except OSError as exc:
            raise _TimeSnapshotDependencyError(
                f"cannot verify Skyfield {name} source"
            ) from exc
        if actual != _SOURCE_DIGESTS[name]:
            raise _TimeSnapshotDependencyError(f"Skyfield {name} source differs")


def _validated_pins(pins: Mapping[str, _AssetPin]) -> tuple[_AssetPin, ...]:
    """Require one correctly associated independent record per fixed role."""
    if set(pins) != set(_ROLES):
        raise _TimeSnapshotIntegrityError("asset roles differ from fixed set")
    ordered = tuple(pins[role] for role in _ROLES)
    for role, pin in zip(_ROLES, ordered):
        if (
            type(pin) is not _AssetPin
            or type(pin.role) is not str
            or pin.role != role
            or type(pin.byte_count) is not int
            or pin.byte_count < 0
            or type(pin.sha256) is not str
            or len(pin.sha256) != 64
            or any(c not in "0123456789abcdef" for c in pin.sha256)
        ):
            raise _TimeSnapshotIntegrityError(f"invalid pin for {role}")
    return ordered


def _validated_shape(shape: _TimeShape) -> None:
    """Reject nonpositive or nonintegral synthetic shape declarations."""
    if type(shape) is not _TimeShape or any(
        type(value) is not int or value < 1
        for value in (shape.daily, shape.leaps, shape.historic, shape.splines)
    ):
        raise _TimeSnapshotParseError("invalid expected array shape")


def _copy_one(path: str | os.PathLike[str], pin: _AssetPin) -> bytes:
    """Close the sole descriptor before checking and returning copied bytes."""
    parts: list[bytes] = []
    remaining = pin.byte_count + 1
    try:
        with open(path, "rb") as source:
            while remaining:
                count = min(remaining, _COPY_CHUNK_SIZE)
                part = source.read(count)
                if type(part) is not bytes or len(part) > count:
                    raise _TimeSnapshotIntegrityError("source read contract failed")
                if not part:
                    break
                parts.append(part)
                remaining -= len(part)
    except OSError as exc:
        raise _TimeSnapshotIntegrityError(f"cannot copy {pin.role}") from exc
    owned = b"".join(parts)
    if len(owned) != pin.byte_count:
        raise _TimeSnapshotIntegrityError(f"{pin.role} byte count differs")
    if hashlib.sha256(owned).hexdigest() != pin.sha256:
        raise _TimeSnapshotIntegrityError(f"{pin.role} SHA-256 differs")
    return owned


def _copy_all(
    paths: Mapping[str, str | os.PathLike[str]], pins: tuple[_AssetPin, ...]
) -> tuple[bytes, ...]:
    """Verify all three full copies before opening any array archive."""
    if set(paths) != set(_ROLES):
        raise _TimeSnapshotIntegrityError("asset path roles differ")
    return tuple(_copy_one(paths[role], pin) for role, pin in zip(_ROLES, pins))


def _float64_array(value: np.ndarray, shape: tuple[int, ...], name: str) -> np.ndarray:
    """Require exact native float64 shape and finite members."""
    array = np.array(value, copy=True)
    if array.dtype != np.dtype("float64") or array.shape != shape:
        raise _TimeSnapshotParseError(f"{name} dtype or shape differs")
    if not np.all(np.isfinite(array)):
        raise _TimeSnapshotParseError(f"{name} has nonfinite values")
    return array


def _increasing(values: np.ndarray) -> bool:
    """Check strict order without requiring raw offsets or values to rise."""
    return bool(np.all(values[1:] > values[:-1]))


def _parse_arrays(
    owned: tuple[bytes, ...], shape: _TimeShape
) -> tuple[np.ndarray, ...]:
    """Extract only fixed NPZ keys and validate the reviewed array schema."""
    try:
        with np.load(io.BytesIO(owned[0]), allow_pickle=False) as archive:
            keys = {
                "tt_jd_minus_arange",
                "delta_t_1e7",
                "leap_dates",
                "leap_offsets",
            }
            if len(archive.files) != len(keys) or set(archive.files) != keys:
                raise _TimeSnapshotParseError("IERS archive members differ")
            tt_offset = _float64_array(
                archive["tt_jd_minus_arange"], (shape.daily,), "IERS TT offset"
            )
            iers_dt_raw = _float64_array(
                archive["delta_t_1e7"], (shape.daily,), "IERS Delta-T"
            )
            leap_dates = _float64_array(
                archive["leap_dates"], (shape.leaps,), "leap dates"
            )
            leap_offsets = _float64_array(
                archive["leap_offsets"], (shape.leaps,), "leap offsets"
            )
        historic = _float64_array(
            np.load(io.BytesIO(owned[1]), allow_pickle=False),
            (2, shape.historic),
            "historic observations",
        )
        with np.load(io.BytesIO(owned[2]), allow_pickle=False) as archive:
            if len(archive.files) != 1 or set(archive.files) != {"Table-S15.2020.txt"}:
                raise _TimeSnapshotParseError("spline archive members differ")
            s15 = _float64_array(
                archive["Table-S15.2020.txt"], (6, shape.splines), "Table-S15"
            )
    except _TimeSnapshotParseError:
        raise
    except Exception as exc:
        raise _TimeSnapshotParseError("asset array parsing failed") from exc

    tt = tt_offset + np.arange(shape.daily, dtype=np.float64)
    if (
        not np.all(np.isfinite(tt))
        or not _increasing(tt)
        or not _increasing(historic[0])
        or not _increasing(leap_dates)
        or not _increasing(leap_offsets)
        or not np.all(leap_offsets == np.floor(leap_offsets))
        or not np.all(s15[1] > s15[0])
        or not _increasing(s15[0])
        or not np.all(s15[0, 1:] == s15[1, :-1])
    ):
        raise _TimeSnapshotParseError("epoch, leap, or spline order differs")
    return tt_offset, iers_dt_raw, leap_dates, leap_offsets, historic, s15


def _finite_scalar(value: object, label: str) -> None:
    """Reject invalid intermediate scalar words before curve publication."""
    if not bool(np.all(np.isfinite(cast(Any, value)))):
        raise _TimeSnapshotConstructionError(f"nonfinite {label}")


def _build_timescale(
    arrays: tuple[np.ndarray, ...],
) -> tuple[Timescale, tuple[np.ndarray, ...]]:
    """Reproduce Skyfield 1.54's ordered spline graph from owned arrays."""
    tt_offset, dt_raw, leap_dates, leap_offsets, historic, s15 = arrays
    try:
        with np.errstate(over="raise", divide="raise", invalid="raise"):
            n = len(tt_offset)
            iers_tt = tt_offset + np.arange(n, dtype=np.float64)
            iers_dt = dt_raw / 1e7
            mask = historic[0] < iers_tt[0]
            merged_tt = np.concatenate([historic[0][mask], iers_tt])
            merged_dt = np.concatenate([historic[1][mask], iers_dt])
            if not len(merged_tt) or not _increasing(merged_tt):
                raise _TimeSnapshotConstructionError("merged epochs differ")

            p = _sf_curve.Splines([1825.0, 1925.0, 0.0, 32.5, 0.0, -320.0])
            pd = p.derivative
            s = _sf_curve.Splines(s15)
            sd = s.derivative
            width = p.upper[0] - p.lower[0]

            x1 = s.lower[0]
            x0 = x1 - 800.0
            left = _sf_curve.build_spline_given_ends(
                x0, p(x0), pd(x0), x1, s(x1), sd(x1)
            )
            x1 = x0
            x0 = x1 - width
            far_left = _sf_curve.build_spline_given_ends(
                x0, p(x0), pd(x0), x1, p(x1), pd(x1)
            )

            x = (merged_tt[0] - 1721045.0) / 365.25
            i = np.searchsorted(s15[0], x)
            s15 = s15[:, :i]
            if not s15.shape[1]:
                raise _TimeSnapshotConstructionError("empty truncated spline table")
            desired_y = merged_dt[0]
            current_y = s(x)
            sx0, sx1, a3, a2, a1, a0 = s15[:, -1]
            t = (x - sx0) / (sx1 - sx0)
            if t == 0.0:
                raise _TimeSnapshotConstructionError("zero spline join fraction")
            a1 = a1 + (desired_y - current_y) / t
            s15[:, -1] = sx0, sx1, a3, a2, a1, a0

            x0 = (merged_tt[-1] - 1721045.0) / 365.25
            x1 = (x0 + 800.0) // 100.0 * 100.0
            y0 = merged_dt[-1]
            lookback = min(366, len(merged_dt))
            slope = (merged_dt[-1] - merged_dt[-lookback]) * lookback / 365.0
            right = _sf_curve.build_spline_given_ends(x0, y0, slope, x1, p(x1), pd(x1))
            x0 = x1
            x1 = x0 + width
            far_right = _sf_curve.build_spline_given_ends(
                x0, p(x0), pd(x0), x1, p(x1), pd(x1)
            )
            curve = _sf_curve.Splines(
                np.concatenate(
                    (
                        np.array([far_left]).T,
                        np.array([left]).T,
                        s15,
                        np.array([right]).T,
                        np.array([far_right]).T,
                    ),
                    axis=1,
                )
            )
            for label, value in (
                ("left patch", left),
                ("far-left patch", far_left),
                ("spline join", a1),
                ("right patch", right),
                ("far-right patch", far_right),
                ("curve", curve.table),
                ("merged Delta-T", merged_dt),
            ):
                _finite_scalar(value, label)
            function = _sf_time.DeltaT(merged_tt, merged_dt, curve)
            timescale = Timescale(function, leap_dates, leap_offsets)
    except _TimeSnapshotConstructionError:
        raise
    except (
        FloatingPointError,
        OverflowError,
        ValueError,
        TypeError,
        IndexError,
    ) as exc:
        raise _TimeSnapshotConstructionError("cannot construct Delta-T curve") from exc

    held = (*arrays, iers_tt, iers_dt, merged_tt, merged_dt, curve.table)
    for array in held:
        array.setflags(write=False)
    for array in (
        curve.lower,
        curve.upper,
        curve.coefficients,
        curve._width,
        curve._n,
    ):
        array.setflags(write=False)
    for array in (
        timescale._leap_utc,
        timescale._leap_offsets,
        timescale._leap_tai,
    ):
        array.setflags(write=False)
    return timescale, held


class _OwnedDefaultTime:
    """Disconnected evaluator bound to three byte snapshots.

    Finite UT inputs with ``abs(ut) <= 1e12`` may be evaluated; larger
    magnitudes are unsupported. This numerical guard is not a statement
    about physical validity or about any ordinary public call.
    """

    __slots__ = ("_owned", "_timescale", "_arrays", "__weakref__")
    _owned: tuple[bytes, ...] | None
    _timescale: Timescale | None
    _arrays: tuple[np.ndarray, ...] | None

    def __init__(self) -> None:
        raise TypeError("use a verified private time factory")

    def _checked_admission(self) -> _Admission:
        """Revalidate exact owned bytes and the evaluator backing a claim."""
        if self._owned is None or self._timescale is None:
            raise _TimeSnapshotClosedError("private default time is closed")
        admission = _ADMISSIONS.get(self)
        if (
            admission is None
            or self._owned is not admission.owned
            or self._timescale is not admission.timescale
        ):
            raise _TimeSnapshotIntegrityError("time evaluator lacks its admission")
        for pin, owned in zip(admission.pins, admission.owned):
            if (
                type(owned) is not bytes
                or len(owned) != pin.byte_count
                or hashlib.sha256(owned).hexdigest() != pin.sha256
            ):
                raise _TimeSnapshotIntegrityError("owned time bytes differ")
        return admission

    @property
    def production_asset_tag(self) -> tuple[tuple[str, int, str], ...] | None:
        """Expose fixed production identities only after full live revalidation."""
        admission = self._checked_admission()
        if not isinstance(admission, _ProductionAdmission):
            return None
        if admission.pins != tuple(_PRODUCTION_PINS[role] for role in _ROLES):
            raise _TimeSnapshotIntegrityError("production pins differ")
        return tuple((pin.role, pin.byte_count, pin.sha256) for pin in admission.pins)

    def evaluate(
        self, ut: float, declaration: _DefaultTimeDeclaration
    ) -> _DefaultTimeEvaluation:
        """Return native default seconds/days and exact words inside the UT guard."""
        timescale = self._checked_admission().timescale
        if type(ut) is not float:
            raise TypeError("UT Julian date must be a native Python float")
        if not math.isfinite(ut):
            raise ValueError("UT Julian date must be finite")
        if (
            type(declaration) is not _DefaultTimeDeclaration
            or declaration.user_delta_t is not None
            or declaration.iers_enabled is not False
            or type(declaration.model) is not str
            or declaration.model != "smh2016"
            or declaration.model_override is not False
            or declaration.tidal_automatic is not True
            or declaration.ordinary_fallback is not False
        ):
            raise _TimeSnapshotConfigurationError("non-default time declaration")
        if abs(ut) > 1e12:
            raise _TimeSnapshotUnsupportedError("UT outside private numerical domain")
        try:
            with np.errstate(over="raise", divide="raise", invalid="raise"):
                time = timescale.ut1_jd(ut)
                tt = float(time.tt)
                seconds = float(time.delta_t)
                days = (seconds + 0.0) / 86400.0
        except (FloatingPointError, OverflowError, ValueError, IndexError) as exc:
            raise _TimeSnapshotUnsupportedError(
                "finite UT cannot be evaluated in private domain"
            ) from exc
        if not all(math.isfinite(word) for word in (tt, seconds, days)):
            raise _TimeSnapshotUnsupportedError("nonfinite private time intermediate")
        return _DefaultTimeEvaluation(
            ut,
            seconds,
            days,
            struct.pack("!d", ut),
            struct.pack("!d", seconds),
            struct.pack("!d", days),
        )

    def _test_leap_words(self) -> tuple[bytes, bytes]:
        """Return copied leap words for comparison with a fresh ordinary scale."""
        timescale = self._checked_admission().timescale
        return timescale.leap_dates.tobytes(), timescale.leap_offsets.tobytes()

    def _test_tt_seconds_bits(self, tt: float) -> bytes:
        """Probe an exact TT spline boundary without UT-to-TT iteration."""
        timescale = self._checked_admission().timescale
        return struct.pack("!d", float(timescale.delta_t_function(tt)))

    def _test_utc_words(
        self,
        year: int,
        month: int,
        day: int,
        hour: int = 0,
        minute: int = 0,
        second: float = 0.0,
    ) -> tuple[bytes, bytes]:
        """Exercise owned leap offsets through UTC-to-TT and UTC-to-TAI."""
        timescale = self._checked_admission().timescale
        time = timescale.utc(year, month, day, hour, minute, second)
        return struct.pack("!d", float(time.tt)), struct.pack("!d", float(time.tai))

    def close(self) -> None:
        """Idempotently release admission, bytes, arrays, and Timescale."""
        _ADMISSIONS.pop(self, None)
        self._owned = None
        self._arrays = None
        self._timescale = None

    def __enter__(self) -> _OwnedDefaultTime:
        self._checked_admission()
        return self

    def __exit__(self, *_exc: object) -> None:
        self.close()


def _open(
    paths: Mapping[str, str | os.PathLike[str]],
    pins: Mapping[str, _AssetPin],
    shape: _TimeShape,
    *,
    production: bool,
) -> _OwnedDefaultTime:
    """Construct only after all byte, code, schema, and curve checks pass."""
    _check_dependencies()
    ordered_pins = _validated_pins(pins)
    _validated_shape(shape)
    owned = _copy_all(paths, ordered_pins)
    arrays = _parse_arrays(owned, shape)
    timescale, held = _build_timescale(arrays)
    reader = object.__new__(_OwnedDefaultTime)
    reader._owned = owned
    reader._timescale = timescale
    reader._arrays = held
    admission_type = _ProductionAdmission if production else _Admission
    try:
        _ADMISSIONS[reader] = admission_type(ordered_pins, owned, timescale)
    except BaseException:
        reader.close()
        raise
    return reader


def _open_production_time(
    paths: Mapping[str, str | os.PathLike[str]],
) -> _OwnedDefaultTime:
    """Admit only the three fixed Skyfield 1.54 / NumPy 2.3.5 assets."""
    return _open(paths, _PRODUCTION_PINS, _PRODUCTION_SHAPE, production=True)


def _open_test_time(
    paths: Mapping[str, str | os.PathLike[str]],
    pins: Mapping[str, _AssetPin],
    shape: _TimeShape,
) -> _OwnedDefaultTime:
    """Exercise synthetic byte identities without issuing a production tag."""
    return _open(paths, pins, shape, production=False)
