# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Analytical theories for planetary moon positions.

This module provides analytical ephemerides for planetary moons, used to
calculate the offset between a planet's barycenter and its center of body (COB).

Theories implemented:
- Galilean moons (Io, Europa, Ganymede, Callisto): Lieske E5 / Meeus algorithm
- Saturn moons (Mimas through Hyperion): TASS 1.7 (Vienne & Duriez)
- Triton: Keplerian elements with J2 precession
- Charon: Two-body Keplerian solution

References:
- Lieske, J.H. (1998) "Galilean satellite ephemerides E5", A&AS 129, 205
- Vienne, A. & Duriez, L. (1995) "TASS1.6", A&A 297, 588-605
- Jacobson, R.A. (2009) "The Orbits of the Neptunian Satellites"
- Brozović, M. & Jacobson, R.A. (2024) "Pluto system orbits", AJ 167:256

License: AGPL-3.0-only, except the
         TASS 1.7 port from Stellarium (`tass17.py`, `tass17_data.py`),
         which is MIT. The Galilean theory is implemented from
         Lieske 1998 / Meeus ch. 44 (see
         docs/methodology/galilean-e5-spec.md).

Provenance:
    This package facade contains no additional theory. Project-owned modules
    state their publication/JPL solution, elements, frame adaptation, units,
    and approximation limits in their docstrings. The two TASS files preserve
    their identified MIT upstream and are inventoried separately in
    ``THIRD_PARTY_NOTICES.md``.
"""

from .constants import (
    # Planet GM values (km³/s²)
    GM_JUPITER,
    GM_SATURN,
    GM_URANUS,
    GM_NEPTUNE,
    GM_PLUTO,
    # Moon GM values (km³/s²)
    GM_IO,
    GM_EUROPA,
    GM_GANYMEDE,
    GM_CALLISTO,
    GM_TITAN,
    GM_TRITON,
    GM_CHARON,
    # Utility functions
    get_cob_offset,
)

__all__ = [
    "GM_JUPITER",
    "GM_SATURN",
    "GM_URANUS",
    "GM_NEPTUNE",
    "GM_PLUTO",
    "GM_IO",
    "GM_EUROPA",
    "GM_GANYMEDE",
    "GM_CALLISTO",
    "GM_TITAN",
    "GM_TRITON",
    "GM_CHARON",
    "get_cob_offset",
]
