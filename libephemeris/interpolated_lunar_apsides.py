# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Alternative asset-free smoothed lunar apsides derived from NASA JPL states.

``INTP_APOG`` and ``INTP_PERG`` are LibEphemeris-defined smooth points.  At
each requested TT instant this module forms the osculating geocentric lunar
eccentricity vector from the active JPL DE state, then applies a symmetric
unit-area boxcar over one local mean anomalistic month.  The integral is
evaluated with Gauss-Legendre quadrature in the fixed ICRS frame.  No sampled
output table, fitted coefficient set, or external compatibility artifact is
used.

This evaluator is retained for research and internal validation. The public
``INTP_APOG`` and ``INTP_PERG`` path uses separate non-antipodal curves through
the DE440 apsis passages; it does not call this one-month boxcar evaluator.

The smoothing kernel is an independent LibEphemeris definition, not a claim
that another ephemeris uses the same interpolation.  The anomalistic period is
derived from ERFA's BSD-licensed implementation of the IERS mean lunar anomaly.
JPL DE440/DE441 are documented by Park et al. (2021):
https://ssd.jpl.nasa.gov/doc/de440_de441.html
"""

from __future__ import annotations

import math
import threading
from collections import OrderedDict
from collections.abc import Callable
from pathlib import Path
from typing import Any

import erfa
import numpy as np

from .exceptions import CalculationError, EphemerisRangeError

INTP_APOG = 21
INTP_PERG = 22
SUPPORTED_BODIES = (INTP_APOG, INTP_PERG)

_J2000_JD = 2451545.0
_JULIAN_CENTURY_DAYS = 36525.0
_PERIOD_RATE_HALF_SPAN_DAYS = 0.5
_QUADRATURE_ORDER = 12
_VELOCITY_DIVISOR = 4096.0

# Generated from the Legendre polynomial itself at import time; these are
# mathematical quadrature values, not ephemeris samples.
_nodes_array, _weights_array = np.polynomial.legendre.leggauss(_QUADRATURE_ORDER)
_QUADRATURE_NODES = tuple(float(value) for value in _nodes_array)
_QUADRATURE_WEIGHTS = tuple(float(value) for value in _weights_array)
del _nodes_array, _weights_array

StateVector = tuple[float, float, float, float, float, float]
CartesianVector = tuple[float, float, float]
_KernelContext = tuple[str | None, float, float, int]
LunarStateProvider = Callable[[float], tuple[CartesianVector, CartesianVector]]


def _angle_delta(after: float, before: float) -> float:
    return math.atan2(math.sin(after - before), math.cos(after - before))


def _mean_anomalistic_period_days(jd: float) -> float:
    """Return the local IERS mean-anomaly period in days."""
    half = _PERIOD_RATE_HALF_SPAN_DAYS
    before = float(erfa.fal03((jd - half - _J2000_JD) / _JULIAN_CENTURY_DAYS))
    after = float(erfa.fal03((jd + half - _J2000_JD) / _JULIAN_CENTURY_DAYS))
    rate = abs(_angle_delta(after, before) / (2.0 * half))
    if not math.isfinite(rate) or rate <= 0.0:
        raise CalculationError("Cannot derive the local lunar anomalistic period")
    return float(2.0 * math.pi / rate)


def _active_jpl_context() -> _KernelContext:
    """Return path, exact JD coverage, and DE number for the active kernel."""
    from .state import get_current_file_data, get_planets

    get_planets()
    path, start_jd, end_jd, denum = get_current_file_data(0)
    return path, float(start_jd), float(end_jd), int(denum)


def _check_range(jd: float, body_id: int, context: _KernelContext) -> None:
    path, start_jd, end_jd, _ = context
    body_name = (
        "interpolated lunar apogee"
        if body_id == INTP_APOG
        else ("interpolated lunar perigee")
    )
    if not math.isfinite(jd) or (start_jd and end_jd and not start_jd <= jd <= end_jd):
        raise EphemerisRangeError(
            f"Julian day {jd} is outside the active JPL range {start_jd} to {end_jd}",
            requested_jd=jd,
            start_jd=start_jd or None,
            end_jd=end_jd or None,
            body_id=body_id,
            body_name=body_name,
            ephemeris_file=path,
        )


def _moon_icrs_state(jd: float) -> tuple[CartesianVector, CartesianVector]:
    """Return geocentric lunar position and velocity in ICRS, AU and AU/day."""
    from .cache import get_cached_time_tt
    from .state import get_planets

    planets = get_planets()
    state = (planets["moon"] - planets["earth"]).at(get_cached_time_tt(jd))
    position = state.position.au
    velocity = state.velocity.au_per_d
    return (
        (float(position[0]), float(position[1]), float(position[2])),
        (float(velocity[0]), float(velocity[1]), float(velocity[2])),
    )


def _osculating_eccentricity_and_p(
    jd: float, state_provider: LunarStateProvider = _moon_icrs_state
) -> tuple[CartesianVector, float]:
    """Return the JPL-derived eccentricity vector and semilatus rectum."""
    from .lunar import GM_EARTH_MOON_AU3_DAY2

    position, velocity = state_provider(jd)
    rx, ry, rz = position
    vx, vy, vz = velocity
    hx = ry * vz - rz * vy
    hy = rz * vx - rx * vz
    hz = rx * vy - ry * vx
    radius = math.sqrt(rx * rx + ry * ry + rz * rz)
    vxh_x = vy * hz - vz * hy
    vxh_y = vz * hx - vx * hz
    vxh_z = vx * hy - vy * hx
    mu = GM_EARTH_MOON_AU3_DAY2
    eccentricity_vector = (
        float(vxh_x / mu - rx / radius),
        float(vxh_y / mu - ry / radius),
        float(vxh_z / mu - rz / radius),
    )
    semilatus_rectum = (hx * hx + hy * hy + hz * hz) / mu
    return eccentricity_vector, float(semilatus_rectum)


def _smoothed_elements(
    jd: float,
    context: _KernelContext,
    state_provider: LunarStateProvider = _moon_icrs_state,
) -> tuple[CartesianVector, float]:
    """Apply the symmetric one-anomalistic-month boxcar in ICRS."""
    _, start_jd, end_jd, _ = context
    nominal_half_span = 0.5 * _mean_anomalistic_period_days(jd)
    if start_jd and end_jd:
        half_span = min(nominal_half_span, jd - start_jd, end_jd - jd)
    else:
        half_span = nominal_half_span
    half_span = max(0.0, float(half_span))

    mean_ex = mean_ey = mean_ez = mean_p = 0.0
    for node, weight in zip(_QUADRATURE_NODES, _QUADRATURE_WEIGHTS):
        eccentricity_vector, semilatus_rectum = _osculating_eccentricity_and_p(
            jd + half_span * node, state_provider
        )
        normalized_weight = 0.5 * weight
        mean_ex += normalized_weight * eccentricity_vector[0]
        mean_ey += normalized_weight * eccentricity_vector[1]
        mean_ez += normalized_weight * eccentricity_vector[2]
        mean_p += normalized_weight * semilatus_rectum
    return (float(mean_ex), float(mean_ey), float(mean_ez)), float(mean_p)


def _icrs_to_ecliptic_of_date(vector: CartesianVector, jd: float) -> CartesianVector:
    from skyfield.framelib import ecliptic_frame

    from .cache import get_cached_time_tt

    rotation = ecliptic_frame.rotation_at(get_cached_time_tt(jd))
    transformed = rotation @ np.asarray(vector, dtype=np.float64)
    return float(transformed[0]), float(transformed[1]), float(transformed[2])


def _cartesian_position(
    jd: float,
    body_id: int,
    context: _KernelContext,
    state_provider: LunarStateProvider = _moon_icrs_state,
) -> CartesianVector:
    eccentricity_vector, semilatus_rectum = _smoothed_elements(
        jd, context, state_provider
    )
    eccentricity = math.sqrt(sum(value * value for value in eccentricity_vector))
    if not math.isfinite(eccentricity) or eccentricity <= 1e-15:
        raise CalculationError("Smoothed lunar eccentricity vector is degenerate")
    direction_sign = -1.0 if body_id == INTP_APOG else 1.0
    denominator = 1.0 - eccentricity if body_id == INTP_APOG else (1.0 + eccentricity)
    if denominator <= 0.0:
        raise CalculationError("Smoothed lunar orbit is not elliptic")
    distance = semilatus_rectum / denominator
    icrs = tuple(
        float(direction_sign * distance * value / eccentricity)
        for value in eccentricity_vector
    )
    return _icrs_to_ecliptic_of_date(icrs, jd)  # type: ignore[arg-type]


def _cartesian_state(
    jd: float,
    body_id: int,
    context: _KernelContext,
    state_provider: LunarStateProvider = _moon_icrs_state,
) -> StateVector:
    """Evaluate position and a bounded centered numerical derivative."""
    position = _cartesian_position(jd, body_id, context, state_provider)
    _, start_jd, end_jd, _ = context
    nominal_step = _mean_anomalistic_period_days(jd) / _VELOCITY_DIVISOR
    if not start_jd or not end_jd:
        before = _cartesian_position(
            jd - nominal_step, body_id, context, state_provider
        )
        after = _cartesian_position(jd + nominal_step, body_id, context, state_provider)
        velocity = tuple(
            (after[index] - before[index]) / (2.0 * nominal_step) for index in range(3)
        )
    elif jd - nominal_step >= start_jd and jd + nominal_step <= end_jd:
        before = _cartesian_position(
            jd - nominal_step, body_id, context, state_provider
        )
        after = _cartesian_position(jd + nominal_step, body_id, context, state_provider)
        velocity = tuple(
            (after[index] - before[index]) / (2.0 * nominal_step) for index in range(3)
        )
    elif jd + 2.0 * nominal_step <= end_jd:
        first = _cartesian_position(jd + nominal_step, body_id, context, state_provider)
        second = _cartesian_position(
            jd + 2.0 * nominal_step, body_id, context, state_provider
        )
        velocity = tuple(
            (-3.0 * position[index] + 4.0 * first[index] - second[index])
            / (2.0 * nominal_step)
            for index in range(3)
        )
    elif jd - 2.0 * nominal_step >= start_jd:
        first = _cartesian_position(jd - nominal_step, body_id, context, state_provider)
        second = _cartesian_position(
            jd - 2.0 * nominal_step, body_id, context, state_provider
        )
        velocity = tuple(
            (3.0 * position[index] - 4.0 * first[index] + second[index])
            / (2.0 * nominal_step)
            for index in range(3)
        )
    else:
        raise EphemerisRangeError(
            "The active JPL interval is too short for a lunar-apsis derivative",
            requested_jd=jd,
            start_jd=start_jd,
            end_jd=end_jd,
            body_id=body_id,
        )
    return (
        float(position[0]),
        float(position[1]),
        float(position[2]),
        float(velocity[0]),
        float(velocity[1]),
        float(velocity[2]),
    )


def state_to_polar(state: StateVector) -> StateVector:
    """Convert a non-degenerate Cartesian state to polar coordinates."""
    x, y, z, vx, vy, vz = state
    rho_squared = x * x + y * y
    rho = math.sqrt(rho_squared)
    distance = math.sqrt(rho_squared + z * z)
    if distance == 0.0 or rho == 0.0:
        raise CalculationError("Cannot convert a degenerate lunar-apsis state")
    radial_speed = (x * vx + y * vy + z * vz) / distance
    longitude_speed = math.degrees((x * vy - y * vx) / rho_squared)
    latitude_speed = math.degrees((vz * distance - z * radial_speed) / (distance * rho))
    return (
        float(math.degrees(math.atan2(y, x)) % 360.0),
        float(math.degrees(math.atan2(z, rho))),
        float(distance),
        float(longitude_speed),
        float(latitude_speed),
        float(radial_speed),
    )


def evaluate_smoothed_lunar_apsis_state(
    jd_tt: float,
    body_id: int,
    state_provider: LunarStateProvider,
    *,
    jd_start: float,
    jd_end: float,
) -> StateVector:
    """Evaluate with an injected backend state provider.

    Backends can use this entry point to supply geocentric lunar ICRS states
    from local JPL, LEB, or Horizons without changing the smoothing definition.
    ``jd_start`` and ``jd_end`` are the provider's exact inclusive coverage.
    """
    if body_id not in SUPPORTED_BODIES:
        raise ValueError(f"Unsupported interpolated lunar apsis body: {body_id}")
    jd = float(jd_tt)
    context: _KernelContext = ("injected-state-provider", jd_start, jd_end, 0)
    _check_range(jd, int(body_id), context)
    return state_to_polar(_cartesian_state(jd, int(body_id), context, state_provider))


class InterpolatedLunarApsidesReader:
    """Bounded, thread-safe evaluator retaining the former reader API.

    No file is opened.  ``path`` and ``verify_payload`` remain accepted so
    callers written for the former artifact reader fail neither at import nor
    at construction; a non-``None`` path is rejected because external model
    artifacts are no longer part of this implementation.
    """

    def __init__(
        self,
        path: str | Path | None = None,
        *,
        cache_size: int = 8,
        verify_payload: bool = False,
    ) -> None:
        del verify_payload
        if path is not None:
            raise ValueError("external interpolated-apsides artifacts are unsupported")
        if (
            not isinstance(cache_size, int)
            or isinstance(cache_size, bool)
            or cache_size < 1
        ):
            raise ValueError("cache_size must be at least 1")
        self._cache_size = int(cache_size)
        self._cache: OrderedDict[tuple[float, int, bool], StateVector] = OrderedDict()
        self._lock = threading.RLock()
        self._closed = False
        self._context = _active_jpl_context()
        self.jd_start = self._context[1]
        self.jd_end = self._context[2]
        self.speed3_jd_start = self.jd_start
        self.speed3_jd_end = self.jd_end

    @property
    def closed(self) -> bool:
        return self._closed

    @property
    def supports_speed3(self) -> bool:
        return True

    def _refresh_context(self) -> None:
        context = _active_jpl_context()
        if context != self._context:
            self._context = context
            self.jd_start = context[1]
            self.jd_end = context[2]
            self.speed3_jd_start = self.jd_start
            self.speed3_jd_end = self.jd_end
            self._cache.clear()

    def _require_open(self) -> None:
        if self._closed:
            raise ValueError("Interpolated lunar apsides evaluator is closed")

    def state(self, jd_tt: float, body_id: int, *, speed3: bool = False) -> StateVector:
        """Evaluate a smoothed JPL-derived Cartesian state."""
        if body_id not in SUPPORTED_BODIES:
            raise ValueError(f"Unsupported interpolated lunar apsis body: {body_id}")
        jd = float(jd_tt)
        with self._lock:
            self._require_open()
            self._refresh_context()
            _check_range(jd, body_id, self._context)
            key = (jd, int(body_id), bool(speed3))
            cached = self._cache.get(key)
            if cached is not None:
                self._cache.move_to_end(key)
                return cached
            # The scientific model already returns a smooth instantaneous
            # derivative, so SPEED and SPEED3 share the same state.
            result = _cartesian_state(jd, int(body_id), self._context)
            self._cache[key] = result
            self._cache.move_to_end(key)
            while len(self._cache) > self._cache_size:
                self._cache.popitem(last=False)
            return result

    def polar(self, jd_tt: float, body_id: int, *, speed3: bool = False) -> StateVector:
        """Evaluate a smoothed JPL-derived polar state."""
        return state_to_polar(self.state(jd_tt, body_id, speed3=speed3))

    def close(self) -> None:
        """Clear the bounded cache and close this evaluator."""
        with self._lock:
            self._cache.clear()
            self._closed = True

    def __enter__(self) -> InterpolatedLunarApsidesReader:
        self._require_open()
        return self

    def __exit__(self, exc_type: Any, exc: Any, traceback: Any) -> None:
        self.close()


def open_interpolated_lunar_apsides(
    *, cache_size: int = 8, verify_payload: bool = False
) -> InterpolatedLunarApsidesReader:
    """Open an asset-free smoothed lunar-apsides evaluator."""
    return InterpolatedLunarApsidesReader(
        cache_size=cache_size, verify_payload=verify_payload
    )


__all__ = [
    "INTP_APOG",
    "INTP_PERG",
    "SUPPORTED_BODIES",
    "InterpolatedLunarApsidesReader",
    "LunarStateProvider",
    "StateVector",
    "evaluate_smoothed_lunar_apsis_state",
    "open_interpolated_lunar_apsides",
    "state_to_polar",
]
