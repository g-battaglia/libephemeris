# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Canonical registry of the exotic minor bodies precomputed into LEB.

Single source of truth consumed by:
  - ``leb_format.BODY_PARAMS``         (Chebyshev interval/degree)
  - ``leb_compression.BODY_TARGET_AU`` (LEB2 mantissa-truncation target)
  - ``scripts/generate_leb.py``        (_ASTEROID_NAIF / BODY_GROUPS / BODY_NAMES)
  - ``scripts/generate_leb2.py``       (LEB2_GROUPS / BODY_NAMES)

Each body is served from a JPL SPK kernel over its coverage window (max
1600-2500 from Horizons) and falls back to the existing Keplerian path
outside it. ``body_id`` uses the runtime ``AST_OFFSET + number`` convention
(e.g. Eris = 10000 + 136199 = 146199), except Pholus which has the dedicated
id 16. ``naif`` = asteroid number + 2000000 (Pholus = 2005145).

``cls`` groups bodies by orbital class. It drives the future extended-tier
(N-body) inclusion: ``nea`` bodies are chaotic and are excluded from the
extended tier; ``tno``/``centaur``/``mainbelt`` are integrable over millennia.

Bennu is intentionally absent: JPL Horizons blocks SPK generation for it
(see ``constants.SPK_AUTO_DOWNLOAD_BLOCKED``).
"""

from __future__ import annotations

from typing import Dict, List, NamedTuple


class ExoticBody(NamedTuple):
    """One precomputed exotic minor body."""

    body_id: int
    naif: int
    interval_days: float
    degree: int
    target_au: float
    name: str
    cls: str  # "tno" | "centaur" | "mainbelt" | "nea"


# Per-class Chebyshev defaults (the fit is on barycentric Cartesian, so the
# binding quantity is the peak heliocentric angular rate). Anchors: the main
# asteroids + Chiron use (8, 13); the Moon uses (4, 13) at ~13 deg/day.
_DEFAULT_TARGET = 5e-9  # 0.001" @ 1 AU — fine for d_geo >= ~1 AU

_REGISTRY: List[ExoticBody] = [
    # --- TNOs (distant, slow): (32, 13), default target ---
    ExoticBody(146199, 2136199, 32, 13, _DEFAULT_TARGET, "Eris", "tno"),
    ExoticBody(100377, 2090377, 32, 13, _DEFAULT_TARGET, "Sedna", "tno"),
    ExoticBody(146108, 2136108, 32, 13, _DEFAULT_TARGET, "Haumea", "tno"),
    ExoticBody(146472, 2136472, 32, 13, _DEFAULT_TARGET, "Makemake", "tno"),
    ExoticBody(38978, 2028978, 32, 13, _DEFAULT_TARGET, "Ixion", "tno"),
    ExoticBody(100482, 2090482, 32, 13, _DEFAULT_TARGET, "Orcus", "tno"),
    ExoticBody(60000, 2050000, 32, 13, _DEFAULT_TARGET, "Quaoar", "tno"),
    ExoticBody(30000, 2020000, 32, 13, _DEFAULT_TARGET, "Varuna", "tno"),
    ExoticBody(235088, 2225088, 32, 13, _DEFAULT_TARGET, "Gonggong", "tno"),
    # --- Centaurs: (8, 13), default target ---
    ExoticBody(16, 2005145, 8, 13, _DEFAULT_TARGET, "Pholus", "centaur"),
    ExoticBody(17066, 2007066, 8, 13, _DEFAULT_TARGET, "Nessus", "centaur"),
    ExoticBody(18405, 2008405, 8, 13, _DEFAULT_TARGET, "Asbolus", "centaur"),
    ExoticBody(20199, 2010199, 8, 13, _DEFAULT_TARGET, "Chariklo", "centaur"),
    ExoticBody(10944, 2000944, 8, 13, _DEFAULT_TARGET, "Hidalgo", "centaur"),
    # --- Main-belt exotics: (8, 13), default target ---
    ExoticBody(10010, 2000010, 8, 13, _DEFAULT_TARGET, "Hygiea", "mainbelt"),
    ExoticBody(10704, 2000704, 8, 13, _DEFAULT_TARGET, "Interamnia", "mainbelt"),
    ExoticBody(10511, 2000511, 8, 13, _DEFAULT_TARGET, "Davida", "mainbelt"),
    ExoticBody(10052, 2000052, 8, 13, _DEFAULT_TARGET, "Europa", "mainbelt"),
    ExoticBody(10087, 2000087, 8, 13, _DEFAULT_TARGET, "Sylvia", "mainbelt"),
    ExoticBody(10016, 2000016, 8, 13, _DEFAULT_TARGET, "Psyche", "mainbelt"),
    ExoticBody(10080, 2000080, 8, 13, _DEFAULT_TARGET, "Sappho", "mainbelt"),
    ExoticBody(10055, 2000055, 8, 13, _DEFAULT_TARGET, "Pandora", "mainbelt"),
    ExoticBody(11181, 2001181, 8, 13, _DEFAULT_TARGET, "Lilith-ast", "mainbelt"),
    # --- NEAs (close geocentric → tighter AU targets; best-effort accuracy) ---
    ExoticBody(10433, 2000433, 4, 13, 1e-10, "Eros", "nea"),
    ExoticBody(11221, 2001221, 4, 13, 1e-10, "Amor", "nea"),
    ExoticBody(109942, 2099942, 4, 13, 1e-12, "Apophis", "nea"),  # 2029 flyby
    ExoticBody(35143, 2025143, 4, 13, 1e-10, "Itokawa", "nea"),
    ExoticBody(172173, 2162173, 4, 13, 1e-10, "Ryugu", "nea"),
    ExoticBody(11685, 2001685, 2, 13, 1e-11, "Toro", "nea"),  # eccentric
    ExoticBody(14179, 2004179, 2, 13, 1e-11, "Toutatis", "nea"),  # eccentric
    ExoticBody(11566, 2001566, 1, 15, 1e-11, "Icarus", "nea"),  # extreme q=0.19 AU
]

# Indexed view: body_id -> ExoticBody
EXOTIC_BODIES: Dict[int, ExoticBody] = {b.body_id: b for b in _REGISTRY}

# Convenience exports for the consumers.
EXOTIC_IDS: List[int] = sorted(EXOTIC_BODIES)
# Bodies eligible for the future N-body extended tier (regular, non-chaotic).
EXOTIC_EXTENDED_IDS: List[int] = sorted(b.body_id for b in _REGISTRY if b.cls != "nea")


def naif_map() -> Dict[int, int]:
    """``body_id -> NAIF id`` for the SPK lookup table (``_ASTEROID_NAIF``)."""
    return {b.body_id: b.naif for b in _REGISTRY}


def body_params() -> Dict[int, tuple]:
    """``body_id -> (interval_days, degree)`` for ``BODY_PARAMS`` extension."""
    return {b.body_id: (b.interval_days, b.degree) for b in _REGISTRY}


def target_au_map() -> Dict[int, float]:
    """``body_id -> target_au`` for non-default ``BODY_TARGET_AU`` entries."""
    return {b.body_id: b.target_au for b in _REGISTRY if b.target_au != _DEFAULT_TARGET}


def name_map() -> Dict[int, str]:
    """``body_id -> display name`` for ``BODY_NAMES``."""
    return {b.body_id: b.name for b in _REGISTRY}
