# SPDX-License-Identifier: AGPL-3.0-only
"""Compose native local and PostgreSQL-backed LEB2 readers."""

from __future__ import annotations

import os

from libephemeris import get_precision_tier
from libephemeris.download import DATA_FILES
from libephemeris.leb_composite import CompositeLEBReader, TieredLEBReader
from libephemeris.leb_groups import LEB2_GROUPS
from libephemeris.state import _discover_reviewed_leb_tier_cores

from .db import CoefficientSourceError, ping as ping, query
from .source import PgByteSource, PgLEB2Reader

__version__ = "0.1.0"


def open_reader() -> TieredLEBReader:
    """Open installed local tiers and manually selected, complete PostgreSQL tiers."""
    order = ("base", "medium", "extended")
    remote = {
        value.strip().lower()
        for value in os.environ.get("LIBEPHEMERIS_PG_TIERS", "medium,extended").split(
            ","
        )
    }
    opened = []
    tiers: dict[str, CompositeLEBReader] = {}
    try:
        if not remote or not remote <= set(order):
            raise ValueError("Invalid remote tiers")
        cores = _discover_reviewed_leb_tier_cores()
        for tier in order[: order.index(get_precision_tier()) + 1]:
            if tier not in remote:
                if tier in cores:
                    tiers[tier] = CompositeLEBReader.from_file_with_companions(
                        cores[tier], pinned_only=True
                    )
                    opened.extend(tiers[tier]._readers)
                continue
            names = sorted(
                (f"{tier}_{group}.leb2" for group in LEB2_GROUPS),
                key=lambda name: (not name.endswith("_core.leb2"), name),
            )
            pins = [DATA_FILES[name]["sha256"] for name in names]
            files = {
                sha: (name, size)
                for sha, name, size in query(
                    "SELECT sha256,name,size FROM libephemeris.files "
                    "WHERE sha256=ANY(%s) AND complete",
                    (pins,),
                )
            }
            if set(files) != set(pins) or any(
                files[sha][0] != name for name, sha in zip(names, pins)
            ):
                raise ValueError("Incomplete pinned files")
            children = []
            for name, sha in zip(names, pins):
                source = PgByteSource(sha, files[sha][1])
                try:
                    reader = PgLEB2Reader(f"postgres://{sha}/{name}", data=source)
                    if not reader._chunked:
                        raise ValueError("Only LEB2 v2 is supported")
                except Exception:
                    source.close()
                    raise
                opened.append(reader)
                children.append(reader)
            tiers[tier] = CompositeLEBReader(children)
        for child in opened:
            child._manifest_verified = True
        return TieredLEBReader(tiers)
    except Exception:
        for child in opened:
            child.close()
        raise CoefficientSourceError(
            "Could not initialize PostgreSQL LEB2 reader"
        ) from None
