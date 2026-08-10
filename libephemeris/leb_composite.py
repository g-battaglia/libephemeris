# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Composite reader that wraps multiple LEB readers (LEB1 and/or LEB2).

Dispatches eval_body() to the reader that contains the requested body.
Auxiliary data (nutation, delta-T, stars) is served from the first reader
that has it.

This enables the modular packaging strategy where bodies are split across
multiple files (for example, core.leb, asteroids.leb, and apogee.leb).

Provenance:
    Project-authored routing infrastructure over already validated LEB readers.
    First-file precedence, body-to-reader mapping, discovery order, and auxiliary
    channel selection are explicit container policy. This module evaluates no
    polynomial itself and contains no astronomical model or coefficient.
"""

from __future__ import annotations

import glob
import os
from contextlib import suppress
from typing import Any, Dict, List, Mapping, Optional, Tuple

from .leb_format import StarEntry

# Filename group suffixes recognized by the composite heuristics
# ({tier}_{group} naming). ``uranians`` is legacy-only: pre-3.1.0 installs
# left those companions on disk, and discovery must classify the filename —
# but never open or attach the file. Its fictitious channels are retired;
# attaching one would resurface source="LEB" coverage for bodies 40-47 and
# break the 3.1.0 invariant that fictitious provenance is "Analytical".
_GROUP_SUFFIXES = frozenset({"core", "asteroids", "apogee", "uranians", "exotics"})
_RETIRED_GROUP_SUFFIXES = frozenset({"uranians"})


class CompositeLEBReader:
    """Wraps multiple LEB readers and dispatches by body_id.

    Usage:
        reader = CompositeLEBReader.from_directory("/path/to/leb/")
        pos, vel = reader.eval_body(SUN, jd_tt)
    """

    def __init__(self, readers: List) -> None:
        """Create a composite reader from a list of LEBReader/LEB2Reader instances.

        Args:
            readers: List of reader instances. Must have at least one.
        """
        if not readers:
            raise ValueError("CompositeLEBReader requires at least one reader")

        self._readers = readers
        self._body_map: Dict[int, Any] = {}  # body_id -> reader

        # Build body -> reader dispatch map
        for reader in readers:
            for body_id in self._get_body_ids(reader):
                if body_id not in self._body_map:
                    self._body_map[body_id] = reader

        # Expose _bodies for fast_calc.py compatibility (accesses reader._bodies[ipl])
        self._bodies: Dict[int, Any] = {}
        for reader in readers:
            for body_id, entry in reader._bodies.items():
                if body_id not in self._bodies:
                    self._bodies[body_id] = entry

        # Find first reader with auxiliary data
        self._nutation_reader = None
        self._delta_t_reader = None
        self._star_reader = None
        for reader in readers:
            if (
                self._nutation_reader is None
                and hasattr(reader, "has_nutation")
                and reader.has_nutation()
            ):
                self._nutation_reader = reader
            if (
                self._delta_t_reader is None
                and hasattr(reader, "_delta_t_jds")
                and reader._delta_t_jds
            ):
                self._delta_t_reader = reader
            if (
                self._star_reader is None
                and hasattr(reader, "_stars")
                and reader._stars
            ):
                self._star_reader = reader

    @staticmethod
    def _get_body_ids(reader) -> list:
        """Extract body IDs from a reader."""
        return list(reader._bodies.keys())

    @classmethod
    def from_directory(cls, directory: str) -> "CompositeLEBReader":
        """Discover and open all .leb/.leb2 files in a directory.

        Args:
            directory: Path to directory containing .leb or .leb2 files.

        Returns:
            CompositeLEBReader wrapping all discovered readers.
        """
        from .leb_reader import open_leb

        leb_files = sorted(
            glob.glob(os.path.join(directory, "*.leb"))
            + glob.glob(os.path.join(directory, "*.leb2"))
        )
        if not leb_files:
            raise FileNotFoundError(f"No .leb/.leb2 files found in {directory}")

        # Tier detection. Filenames encode the tier either as {tier}_{group}
        # (group-file scheme), as {custom}_{tier}_{group}, or as a bare
        # known-tier token (e.g. ephemeris_base.leb); files with no recognizable
        # tier are not constrained by the guard below.
        _KNOWN_TIERS = {"base", "medium", "extended"}

        def _file_tier(path: str) -> Optional[str]:
            parts = os.path.basename(path).rsplit(".", 1)[0].split("_")
            if len(parts) >= 2 and parts[-1] in _GROUP_SUFFIXES:
                prefix = "_".join(parts[:-1])
                if prefix in _KNOWN_TIERS:
                    return prefix
                # A custom multi-token prefix (e.g. "myset_base") is not itself a
                # known tier; fall through to the token scan so an embedded
                # tier token ("base") is still detected rather than returning None.
            for token in parts:
                if token in _KNOWN_TIERS:
                    return token
            return None

        # Open the readers first, then apply the tier guard to the files that
        # actually opened. Computing tiers from filenames alone would let a
        # corrupt/0-byte different-tier stub (e.g. an interrupted download) abort
        # a composite that the skip-invalid loop would otherwise have built.
        readers = []
        opened_paths = []
        for path in leb_files:
            name_parts = os.path.basename(path).rsplit(".", 1)[0].split("_")
            if name_parts[-1] in _RETIRED_GROUP_SUFFIXES:
                continue  # legacy uranians companion: ignored, never opened
            try:
                readers.append(open_leb(path))
                opened_paths.append(path)
            except (ValueError, OSError):
                continue  # skip invalid files

        if not readers:
            raise ValueError(f"No valid .leb files found in {directory}")

        # Tier guard: refuse to silently merge files from different tiers
        # (e.g. base_core.leb2 + medium_core.leb2) — they share body ids but
        # cover different date ranges with different fits, so a first-wins merge
        # would corrupt positions. Mirrors the tier check in
        # from_file_with_companions.
        tiers = {t for t in (_file_tier(p) for p in opened_paths) if t is not None}
        if len(tiers) > 1:
            for reader in readers:
                with suppress(OSError, ValueError, KeyError, AttributeError):
                    reader.close()
            raise ValueError(
                "Refusing to build a composite from mixed tiers in "
                f"{directory}: {sorted(tiers)}. All files must share one tier "
                "(base/medium/extended)."
            )

        return cls(readers)

    @classmethod
    def from_file_with_companions(
        cls, path: str, pinned_only: bool = False
    ) -> "CompositeLEBReader":
        """Open a .leb/.leb2 file and discover companion files in the same directory.

        Companion files share a common tier prefix. For example, if the primary
        file is ``base_core.leb2``, companions would be ``base_asteroids.leb2``,
        ``base_apogee.leb2``, and ``base_exotics.leb2``.

        If no companions are found, returns a composite with a single reader.

        Args:
            path: Path to the primary .leb or .leb2 file.
            pinned_only: When True (the reviewed-core trust unit), every
                companion is checked against the manifest to attach; a
                candidate that does not resolve is skipped with a log
                message.

        Returns:
            CompositeLEBReader wrapping the primary and companion readers.
        """
        from .leb_reader import open_leb

        directory = os.path.dirname(os.path.abspath(path))
        basename = os.path.basename(path)

        # Try to extract tier prefix (e.g., "base" from "base_core.leb2")
        name_no_ext = basename.rsplit(".", 1)[0]
        parts = name_no_ext.split("_")

        readers = [open_leb(path)]

        # Companions exist only for the group-file naming scheme
        # ({tier}_{group}.leb2).  A merged file (e.g. "ephemeris_base.leb")
        # is complete on its own — a bare first-token prefix match would
        # pull in other tiers ("ephemeris_medium.leb", ...) and stale
        # partials, silently mixing tiers in one composite.
        if len(parts) >= 2 and parts[-1] in _GROUP_SUFFIXES:
            prefix = "_".join(parts[:-1])  # e.g., "base"
            companions = sorted(
                glob.glob(os.path.join(directory, f"{prefix}_*.leb"))
                + glob.glob(os.path.join(directory, f"{prefix}_*.leb2"))
            )
            for companion_path in companions:
                cname = os.path.basename(companion_path).rsplit(".", 1)[0]
                cparts = cname.split("_")
                if cparts[:-1] != parts[:-1] or cparts[-1] not in _GROUP_SUFFIXES:
                    continue
                if cparts[-1] in _RETIRED_GROUP_SUFFIXES:
                    continue  # legacy uranians companion: ignored, never opened
                if os.path.abspath(companion_path) == os.path.abspath(path):
                    continue
                if pinned_only:
                    # Lazy import: state also imports this module lazily, so
                    # there is no module-level cycle.
                    from .logging_config import get_logger
                    from .state import _matches_pinned_data_file

                    if not _matches_pinned_data_file(
                        companion_path, os.path.basename(companion_path)
                    ):
                        get_logger().warning(
                            "Skipping LEB companion %s: it is not a usable "
                            "manifest artifact",
                            companion_path,
                        )
                        continue
                try:
                    readers.append(open_leb(companion_path))
                except (ValueError, OSError):
                    continue

        return cls(readers)

    @property
    def path(self) -> str:
        """Return the path of the first reader."""
        return self._readers[0].path

    @property
    def jd_range(self) -> Tuple[float, float]:
        """Return the widest JD range across all readers."""
        jd_start = min(r.jd_range[0] for r in self._readers)
        jd_end = max(r.jd_range[1] for r in self._readers)
        return (jd_start, jd_end)

    def warm(self, jd_start: float, jd_end: float) -> None:
        """Pre-fault mmap pages across all constituent readers.

        Args:
            jd_start: Start of the Julian Day range to pre-fault.
            jd_end: End of the Julian Day range to pre-fault.
        """
        for reader in self._readers:
            reader.warm(jd_start, jd_end)

    def cool(self) -> None:
        """Advise the kernel that mmap pages can be reclaimed.

        Delegates to all constituent readers.  Idempotent and safe to
        call on closed readers.
        """
        for reader in self._readers:
            try:
                reader.cool()
            except (OSError, AttributeError, ValueError):
                pass

    def has_body(self, body_id: int) -> bool:
        return body_id in self._body_map

    def body_coverage(self, body_id: int) -> Optional[Tuple[float, float]]:
        """Return coverage from the reader that actually serves ``body_id``."""
        reader = self._body_map.get(body_id)
        if reader is None:
            return None
        coverage = getattr(reader, "body_coverage", None)
        if coverage is not None:
            return coverage(body_id)
        entry = getattr(reader, "_bodies", {}).get(body_id)
        if entry is None:
            return None
        return (float(entry.jd_start), float(entry.jd_end))

    def body_reader(self, body_id: int) -> Optional[Any]:
        """Return the constituent reader selected for a body."""
        return self._body_map.get(body_id)

    def eval_body(
        self, body_id: int, jd: float
    ) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
        if body_id not in self._body_map:
            raise KeyError(f"Body {body_id} not in any LEB file")
        return self._body_map[body_id].eval_body(body_id, jd)

    def has_nutation(self) -> bool:
        """Return True if any constituent LEB file contains nutation data."""
        return self._nutation_reader is not None

    def eval_nutation(self, jd_tt: float) -> Tuple[float, float]:
        if self._nutation_reader is None:
            raise ValueError("No nutation data in any LEB file")
        return self._nutation_reader.eval_nutation(jd_tt)

    def delta_t(self, jd: float) -> float:
        if self._delta_t_reader is None:
            raise ValueError("No Delta-T data in any LEB file")
        return self._delta_t_reader.delta_t(jd)

    def get_star(self, star_id: int) -> StarEntry:
        if self._star_reader is None:
            raise KeyError("No star catalog in any LEB file")
        return self._star_reader.get_star(star_id)

    def close(self) -> None:
        for reader in self._readers:
            try:
                reader.close()
            except (OSError, ValueError, KeyError):
                pass
        self._readers.clear()
        self._body_map.clear()
        # Drop aux-reader references too: they point at now-closed readers,
        # and a post-close delta_t()/eval_nutation() should fail with the
        # clean "no data" ValueError rather than a struct/TypeError.
        self._nutation_reader = None
        self._delta_t_reader = None
        self._star_reader = None

    def __enter__(self) -> "CompositeLEBReader":
        return self

    def __exit__(self, *args) -> None:
        self.close()

    def __repr__(self) -> str:
        n_bodies = len(self._body_map)
        n_files = len(self._readers)
        return f"CompositeLEBReader({n_files} files, {n_bodies} bodies)"


class TieredLEBReader:
    """Select the most precise reviewed LEB tier for each body and date.

    The three LEB tiers overlap deliberately.  ``base`` is generated from the
    short DE440 kernel and has the same source precision as DE440 inside its
    modern interval; ``medium`` extends that DE440 solution; ``extended`` is
    generated from DE441 for dates outside DE440.  A first-file-wins composite
    cannot express that policy, because it would use one tier for the entire
    process.  This reader keeps each per-tier composite intact and chooses at
    evaluation time:

    ``base (DE440s) -> medium (DE440) -> extended (DE441)``.

    Selection is per body because companion coverage can be narrower than the
    core.  When no stored interval contains the date, the highest-priority
    reader is allowed to raise its normal range error so the sealed dispatcher
    can apply the body's explicitly declared local model.  No JPL/BSP source is
    opened by this class.
    """

    TIER_PRIORITY = ("base", "medium", "extended")

    def __init__(self, tier_readers: Mapping[str, CompositeLEBReader]) -> None:
        unknown = set(tier_readers) - set(self.TIER_PRIORITY)
        if unknown:
            raise ValueError(f"Unknown LEB tier(s): {sorted(unknown)}")
        if not tier_readers:
            raise ValueError("TieredLEBReader requires at least one tier")

        self._tier_readers = {
            tier: tier_readers[tier]
            for tier in self.TIER_PRIORITY
            if tier in tier_readers
        }
        # Inventory and finalizer consumers already understand ``_readers`` as
        # the concrete file list.  Flatten the per-tier composites without
        # changing their calculation dispatch.
        self._readers = [
            child
            for tier_reader in self._tier_readers.values()
            for child in tier_reader._readers
        ]

        self._body_candidates: Dict[int, List[Tuple[str, Any]]] = {}
        for tier, tier_reader in self._tier_readers.items():
            for body_id in tier_reader._bodies:
                serving = tier_reader.body_reader(body_id)
                if serving is not None:
                    self._body_candidates.setdefault(body_id, []).append(
                        (tier, serving)
                    )

        # ``fast_calc`` uses the entry only for coordinate-channel metadata;
        # every tier stores the same channel for a given body.  Numerical
        # evaluation still goes through ``eval_body`` and is date-selected.
        self._bodies: Dict[int, Any] = {}
        self._body_map: Dict[int, Any] = {}
        for body_id, candidates in self._body_candidates.items():
            serving = candidates[0][1]
            self._body_map[body_id] = serving
            self._bodies[body_id] = serving._bodies[body_id]

    @classmethod
    def from_tier_cores(cls, tier_cores: Mapping[str, str]) -> "TieredLEBReader":
        """Open reviewed core paths and their pinned same-tier companions."""
        tier_readers = {
            tier: CompositeLEBReader.from_file_with_companions(path, pinned_only=True)
            for tier, path in tier_cores.items()
        }
        return cls(tier_readers)

    @property
    def path(self) -> str:
        """Return the first (highest-priority) tier core path."""
        return next(iter(self._tier_readers.values())).path

    @property
    def jd_range(self) -> Tuple[float, float]:
        """Return the union range across installed tier composites."""
        return (
            min(reader.jd_range[0] for reader in self._tier_readers.values()),
            max(reader.jd_range[1] for reader in self._tier_readers.values()),
        )

    @staticmethod
    def _bounds(reader: Any, body_id: int) -> Optional[Tuple[float, float]]:
        coverage = getattr(reader, "body_coverage", None)
        if coverage is None:
            return None
        bounds = coverage(body_id)
        if bounds is None:
            return None
        return (float(bounds[0]), float(bounds[1]))

    def selected_body_reader(self, body_id: int, jd: float) -> Optional[Any]:
        """Return the highest-priority concrete reader covering ``jd``."""
        for _tier, reader in self._body_candidates.get(body_id, ()):
            bounds = self._bounds(reader, body_id)
            if bounds is not None and bounds[0] <= float(jd) <= bounds[1]:
                return reader
        return None

    def selected_tier(self, body_id: int, jd: float) -> Optional[str]:
        """Return the tier selected for ``body_id`` at ``jd``."""
        selected = self.selected_body_reader(body_id, jd)
        if selected is None:
            return None
        for tier, reader in self._body_candidates.get(body_id, ()):
            if reader is selected:
                return tier
        return None

    def body_reader(self, body_id: int, jd: float | None = None) -> Optional[Any]:
        """Return the serving file, optionally selected for a target date."""
        if jd is not None:
            return self.selected_body_reader(body_id, jd)
        candidates = self._body_candidates.get(body_id)
        return candidates[0][1] if candidates else None

    def body_coverage(self, body_id: int) -> Optional[Tuple[float, float]]:
        """Return the union interval available across all installed tiers."""
        ranges = [
            bounds
            for _tier, reader in self._body_candidates.get(body_id, ())
            if (bounds := self._bounds(reader, body_id)) is not None
        ]
        if not ranges:
            return None
        return (
            min(bounds[0] for bounds in ranges),
            max(bounds[1] for bounds in ranges),
        )

    def has_body(self, body_id: int) -> bool:
        return body_id in self._body_candidates

    def eval_body(
        self, body_id: int, jd: float
    ) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
        candidates = self._body_candidates.get(body_id)
        if not candidates:
            raise KeyError(f"Body {body_id} not in any installed LEB tier")
        selected = self.selected_body_reader(body_id, jd)
        if selected is None:
            # Preserve the established range-miss exception type.  The sealed
            # dispatcher decides whether this body is core (fail closed) or a
            # curated companion with an allowed local model.
            selected = candidates[0][1]
        return selected.eval_body(body_id, jd)

    def _select_aux_reader(self, jd: float, capability: str) -> Optional[Any]:
        for tier_reader in self._tier_readers.values():
            start, end = tier_reader.jd_range
            if (
                start <= float(jd) <= end
                and getattr(tier_reader, capability, lambda: False)()
            ):
                return tier_reader
        return None

    def has_nutation(self) -> bool:
        return any(reader.has_nutation() for reader in self._tier_readers.values())

    def eval_nutation(self, jd_tt: float) -> Tuple[float, float]:
        reader = self._select_aux_reader(jd_tt, "has_nutation")
        if reader is None:
            raise ValueError(f"No LEB nutation data covers JD {jd_tt}")
        return reader.eval_nutation(jd_tt)

    def delta_t(self, jd: float) -> float:
        for tier_reader in self._tier_readers.values():
            start, end = tier_reader.jd_range
            if start <= float(jd) <= end:
                try:
                    return tier_reader.delta_t(jd)
                except ValueError:
                    continue
        raise ValueError(f"No LEB Delta-T data covers JD {jd}")

    def get_star(self, star_id: int) -> StarEntry:
        for tier_reader in self._tier_readers.values():
            try:
                return tier_reader.get_star(star_id)
            except KeyError:
                continue
        raise KeyError(f"Star {star_id} not in any installed LEB tier")

    def warm(self, jd_start: float, jd_end: float) -> None:
        for reader in self._tier_readers.values():
            overlap_start = max(float(jd_start), float(reader.jd_range[0]))
            overlap_end = min(float(jd_end), float(reader.jd_range[1]))
            if overlap_start <= overlap_end:
                reader.warm(overlap_start, overlap_end)

    def cool(self) -> None:
        for reader in self._tier_readers.values():
            reader.cool()

    def close(self) -> None:
        for reader in self._tier_readers.values():
            reader.close()
        self._tier_readers.clear()
        self._readers.clear()
        self._body_candidates.clear()
        self._body_map.clear()
        self._bodies.clear()

    def __enter__(self) -> "TieredLEBReader":
        return self

    def __exit__(self, *args) -> None:
        self.close()

    def __repr__(self) -> str:
        tiers = ",".join(self._tier_readers)
        return (
            f"TieredLEBReader(tiers={tiers}, files={len(self._readers)}, "
            f"bodies={len(self._bodies)})"
        )
