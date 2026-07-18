# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Skyfield vector adapters backed exclusively by an active LEB reader.

Several mature calculation paths in :mod:`libephemeris.planets` consume the
small ``VectorFunction`` protocol rather than a concrete SPICE kernel.  This
module lets those paths keep their existing light-time, apparent-place, and
frame-reduction code while sourcing every persisted state from LEB.

Provenance:
    This is project-authored dispatch glue.  It does not add astronomical
    coefficients or alter stored LEB states.  The Earth-Moon barycentre uses
    the same IAU/DE-consistent mass ratio already used by the lunar mechanics
    module.
"""

from __future__ import annotations

import threading
from typing import TYPE_CHECKING, Any, Callable

import numpy as np
from skyfield.vectorlib import VectorFunction

from .constants import (
    EARTH,
    JUPITER,
    MARS,
    MERCURY,
    MOON,
    NEPTUNE,
    PLUTO,
    SATURN,
    SUN,
    URANUS,
    VENUS,
)
from .exceptions import EphemerisRangeError, LEBCorruptionError

if TYPE_CHECKING:
    from .leb2_reader import LEB2Reader
    from .leb_composite import CompositeLEBReader, TieredLEBReader
    from .leb_reader import LEBReader

    LEBReaderLike = LEBReader | LEB2Reader | CompositeLEBReader | TieredLEBReader
else:
    LEBReaderLike = Any


_EMRAT = 81.3005691


def _serving_path(
    reader: LEBReaderLike, body_id: int, jd_tt: float | None = None
) -> str | None:
    """Return the concrete file serving a body, when the reader exposes it."""
    body_reader = getattr(reader, "body_reader", None)
    if body_reader is None:
        serving = reader
    else:
        try:
            serving = body_reader(body_id, jd_tt)
        except TypeError:
            serving = body_reader(body_id)
    path = getattr(serving, "path", None)
    return str(path) if path is not None else None


def _eval_reader_body(
    reader: LEBReaderLike, body_id: int, jd_tt: float
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """Evaluate one LEB state and turn a stored-range miss into a public error."""
    try:
        return reader.eval_body(body_id, jd_tt)
    except LEBCorruptionError:
        # Corruption is a provisioning failure, never a coverage condition.
        raise
    except ValueError as exc:
        coverage_fn = getattr(reader, "body_coverage", None)
        coverage = coverage_fn(body_id) if coverage_fn is not None else None
        if coverage is not None:
            jd_start, jd_end = (float(coverage[0]), float(coverage[1]))
            if jd_tt < jd_start or jd_tt > jd_end:
                raise EphemerisRangeError(
                    message=(
                        f"Body {body_id} at JD {jd_tt:.6f} is outside active "
                        f"LEB coverage range [{jd_start:.6f}, {jd_end:.6f}]."
                    ),
                    requested_jd=jd_tt,
                    start_jd=jd_start,
                    end_jd=jd_end,
                    body_id=body_id,
                    ephemeris_file=_serving_path(reader, body_id, jd_tt),
                ) from exc
        raise


class _LEBVectorTarget(VectorFunction):
    """One SSB-relative Skyfield vector evaluated from LEB coefficients."""

    center = 0

    def __init__(
        self,
        ephemeris: "LEBVectorEphemeris",
        target: int,
        evaluator: Callable[
            [float], tuple[tuple[float, float, float], tuple[float, float, float]]
        ],
        name: str,
    ) -> None:
        self.ephemeris = ephemeris
        self.target = target
        self._evaluator = evaluator
        self._name = name

    @property
    def vector_name(self) -> str:
        """Return a stable diagnostic name without consulting JPL name tables."""
        return f"LEB {self._name}"

    def _at(self, t):
        jd = np.asarray(t.tt, dtype=float)
        if jd.ndim == 0:
            position, velocity = self._evaluator(float(jd))
            return (
                np.asarray(position, dtype=float),
                np.asarray(velocity, dtype=float),
                None,
                None,
            )

        flat_states = [self._evaluator(float(value)) for value in jd.ravel()]
        positions = np.asarray([state[0] for state in flat_states], dtype=float)
        velocities = np.asarray([state[1] for state in flat_states], dtype=float)
        output_shape = (3, *jd.shape)
        return (
            positions.T.reshape(output_shape),
            velocities.T.reshape(output_shape),
            None,
            None,
        )


class LEBVectorEphemeris:
    """Dictionary-like collection of Skyfield vectors sourced from one LEB."""

    source = "LEB"

    def __init__(self, reader: LEBReaderLike) -> None:
        self.reader = reader
        self._targets: dict[str | int, _LEBVectorTarget] = {}
        self._install_reader_target(SUN, 10, "sun", ("sun", 10))
        self._install_reader_target(MOON, 301, "moon", ("moon", 301))
        self._install_reader_target(
            MERCURY, 199, "mercury", ("mercury", "mercury barycenter", 1, 199)
        )
        self._install_reader_target(
            VENUS, 299, "venus", ("venus", "venus barycenter", 2, 299)
        )
        self._install_reader_target(
            MARS, 4, "mars barycenter", ("mars", "mars barycenter", 4)
        )
        self._install_reader_target(
            JUPITER, 5, "jupiter barycenter", ("jupiter", "jupiter barycenter", 5)
        )
        self._install_reader_target(
            SATURN, 6, "saturn barycenter", ("saturn", "saturn barycenter", 6)
        )
        self._install_reader_target(
            URANUS, 7, "uranus barycenter", ("uranus", "uranus barycenter", 7)
        )
        self._install_reader_target(
            NEPTUNE, 8, "neptune barycenter", ("neptune", "neptune barycenter", 8)
        )
        self._install_reader_target(
            PLUTO, 9, "pluto barycenter", ("pluto", "pluto barycenter", 9)
        )
        self._install_reader_target(EARTH, 399, "earth", ("earth", 399))

        if reader.has_body(EARTH) and reader.has_body(MOON):
            emb = _LEBVectorTarget(
                self, 3, self._eval_earth_barycenter, "earth barycenter"
            )
            self._install_aliases(
                emb,
                (
                    "earth barycenter",
                    "earth-moon barycenter",
                    "earth moon barycenter",
                    3,
                ),
            )

    def _install_reader_target(
        self,
        body_id: int,
        target_id: int,
        name: str,
        aliases: tuple[str | int, ...],
    ) -> None:
        if not self.reader.has_body(body_id):
            return

        def evaluate(jd_tt: float, body: int = body_id):
            return _eval_reader_body(self.reader, body, jd_tt)

        target = _LEBVectorTarget(self, target_id, evaluate, name)
        self._install_aliases(target, aliases)

    def _install_aliases(
        self, target: _LEBVectorTarget, aliases: tuple[str | int, ...]
    ) -> None:
        for alias in aliases:
            key = alias.lower() if isinstance(alias, str) else alias
            self._targets[key] = target

    def _eval_earth_barycenter(
        self, jd_tt: float
    ) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
        earth_pos, earth_vel = _eval_reader_body(self.reader, EARTH, jd_tt)
        moon_pos, moon_vel = _eval_reader_body(self.reader, MOON, jd_tt)
        weight = 1.0 / (_EMRAT + 1.0)
        return (
            tuple(
                earth_pos[index] + (moon_pos[index] - earth_pos[index]) * weight
                for index in range(3)
            ),
            tuple(
                earth_vel[index] + (moon_vel[index] - earth_vel[index]) * weight
                for index in range(3)
            ),
        )

    @staticmethod
    def _key(key: str | int) -> str | int:
        return key.lower().strip() if isinstance(key, str) else key

    def __getitem__(self, key: str | int) -> _LEBVectorTarget:
        return self._targets[self._key(key)]

    def __contains__(self, key: object) -> bool:
        if not isinstance(key, (str, int)):
            return False
        return self._key(key) in self._targets

    def __repr__(self) -> str:
        return f"<LEBVectorEphemeris reader={self.reader!r}>"


_CACHE_LOCK = threading.RLock()
_CACHED_READER: LEBReaderLike | None = None
_CACHED_EPHEMERIS: LEBVectorEphemeris | None = None


def get_leb_vector_ephemeris(reader: LEBReaderLike) -> LEBVectorEphemeris:
    """Return the process-cached vector adapter for an active reader."""
    global _CACHED_READER, _CACHED_EPHEMERIS
    if _CACHED_READER is reader and _CACHED_EPHEMERIS is not None:
        return _CACHED_EPHEMERIS
    with _CACHE_LOCK:
        if _CACHED_READER is not reader or _CACHED_EPHEMERIS is None:
            _CACHED_READER = reader
            _CACHED_EPHEMERIS = LEBVectorEphemeris(reader)
        return _CACHED_EPHEMERIS


def reset_leb_vector_ephemeris() -> None:
    """Release references held by the vector-adapter cache."""
    global _CACHED_READER, _CACHED_EPHEMERIS
    with _CACHE_LOCK:
        _CACHED_READER = None
        _CACHED_EPHEMERIS = None


__all__ = [
    "LEBVectorEphemeris",
    "get_leb_vector_ephemeris",
    "reset_leb_vector_ephemeris",
]
