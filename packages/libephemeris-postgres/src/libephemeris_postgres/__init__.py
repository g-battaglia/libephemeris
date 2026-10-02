# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Compose native local and PostgreSQL-backed LEB2 readers."""

from __future__ import annotations

import os

from libephemeris import CoefficientSourceError, get_precision_tier
from libephemeris.download import DATA_FILES
from libephemeris.leb2_remote import RemoteLEB2Reader
from libephemeris.leb_composite import CompositeLEBReader, TieredLEBReader
from libephemeris.leb_groups import LEB2_GROUPS
from libephemeris.state import _discover_reviewed_leb_tier_cores

from .config import RuntimeConfig, runtime_config
from .pool import get_pool
from .source import PgByteSource

__version__ = "0.1.0"


def open_tier(
    tier: str, *, config: RuntimeConfig | None = None
) -> list[RemoteLEB2Reader]:
    """Open the four complete files identified by the installed manifest pins."""
    if tier not in {"base", "medium", "extended"}:
        raise CoefficientSourceError("Invalid PostgreSQL tier")
    readers: list[RemoteLEB2Reader] = []
    try:
        config = config or runtime_config()
        names = [f"{tier}_{group}.leb2" for group in LEB2_GROUPS]
        pins = [DATA_FILES[name]["sha256"] for name in names]
        with get_pool(config).connection() as conn:
            rows = conn.execute(
                "SELECT sha256,name,size FROM libephemeris.files "
                "WHERE sha256=ANY(%s) AND complete",
                (pins,),
            ).fetchall()
        files = {sha: (name, int(size)) for sha, name, size in rows}
        if set(files) != set(pins) or any(
            files[sha][0] != name for name, sha in zip(names, pins)
        ):
            raise CoefficientSourceError("Pinned PostgreSQL LEB2 files are incomplete")
        for name, sha in zip(names, pins):
            readers.append(
                RemoteLEB2Reader(
                    PgByteSource(name, sha, files[sha][1], config), reviewed=True
                )
            )
        return readers
    except Exception:
        for reader in readers:
            reader.close()
        raise CoefficientSourceError(
            "PostgreSQL coefficient source unavailable"
        ) from None


def open_reader() -> TieredLEBReader:
    """Compose manually selected remote tiers and reviewed local lower tiers."""
    order = ("base", "medium", "extended")
    remote = {
        value.strip()
        for value in os.environ.get("LIBEPHEMERIS_PG_TIERS", "medium,extended").split(
            ","
        )
    }
    if not remote or not remote <= set(order):
        raise CoefficientSourceError("Invalid LIBEPHEMERIS_PG_TIERS")
    tiers: dict[str, CompositeLEBReader] = {}
    try:
        cores = _discover_reviewed_leb_tier_cores()
        for tier in order[: order.index(get_precision_tier()) + 1]:
            if tier in remote:
                children = sorted(
                    open_tier(tier),
                    key=lambda r: (not r.path.endswith("_core.leb2"), r.path),
                )
                tiers[tier] = CompositeLEBReader(children)
            elif tier in cores:
                local = CompositeLEBReader.from_file_with_companions(
                    cores[tier], pinned_only=True
                )
                for child in local._readers:
                    child._manifest_verified = True
                tiers[tier] = local
        result = TieredLEBReader(tiers)
        result._manifest_verified = all(  # type: ignore[attr-defined]
            getattr(child, "_manifest_verified", False) for child in result._readers
        )
        return result
    except Exception:
        for opened in tiers.values():
            opened.close()
        raise CoefficientSourceError(
            "Could not initialize PostgreSQL LEB2 reader"
        ) from None
