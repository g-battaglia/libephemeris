# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Composite reader that wraps multiple LEB readers (LEB1 and/or LEB2).

Dispatches eval_body() to the reader that contains the requested body.
Auxiliary data (nutation, delta-T, stars) is served from the first reader
that has it.

This enables the modular packaging strategy where bodies are split across
multiple files (for example, core.leb, asteroids.leb, and apogee.leb).
"""

from __future__ import annotations

import glob
import os
from contextlib import suppress
from typing import Any, Dict, List, Optional, Tuple

from .leb_format import StarEntry


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
        # ``uranians`` is recognized only as a legacy filename for tier-safety;
        # fictitious-body channels are rejected by the public calculation path.
        _GROUP_SUFFIXES = {"core", "asteroids", "apogee", "uranians", "exotics"}
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
    def from_file_with_companions(cls, path: str) -> "CompositeLEBReader":
        """Open a .leb/.leb2 file and discover companion files in the same directory.

        Companion files share a common tier prefix. For example, if the primary
        file is ``base_core.leb2``, companions would be ``base_asteroids.leb2``,
        ``base_apogee.leb2``, and ``base_exotics.leb2``.

        If no companions are found, returns a composite with a single reader.

        Args:
            path: Path to the primary .leb or .leb2 file.

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
        # Keep recognizing the retired suffix so a legacy file cannot evade the
        # mixed-tier guard. Public calculations reject its fictitious channels.
        _GROUP_SUFFIXES = {"core", "asteroids", "apogee", "uranians", "exotics"}

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
                if os.path.abspath(companion_path) != os.path.abspath(path):
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
