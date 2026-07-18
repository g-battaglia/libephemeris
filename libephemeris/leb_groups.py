# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Canonical LEB body-group names — single source of truth.

The LEB system partitions bodies two ways: a generation/merge partition
(LEB1) and a distribution partition (LEB2, one ``{tier}_{group}.leb2``
companion per group). Several modules need these group names — the
generators, the merge workflow, the runtime downloader, the setup wizard,
and status diagnostics — and historically each kept its own hard-coded copy.
That let a new group (e.g. ``exotics``) be added to the data structures yet
silently omitted from a workflow. Importing the names from here keeps every
consumer in sync: add a group in ONE place and it propagates.

The body-to-group mapping itself lives with the generators
(``scripts/generate_leb.py`` ``BODY_GROUPS``, ``scripts/generate_leb2.py``
``LEB2_GROUPS``); the keys of those dicts must match the tuples below.

Provenance:
    Project-authored build/distribution metadata with no astronomical model or
    coefficient. Group names affect file partitioning only; a body's source,
    coordinate channel, and approximation parameters remain in the registered
    generator/format records.
"""

from __future__ import annotations

# LEB1 merged-main partition (order = generation + merge order).
LEB1_GROUPS = ("planets", "asteroids", "exotics", "analytical")

# Standalone LEB1 companions that must be generated and verified but never
# merged into ``ephemeris_<tier>.leb``.  They are independent conversion
# sources for the matching LEB2 companion.
LEB1_COMPANION_GROUPS = ("uranians",)
LEB1_GENERATION_GROUPS = LEB1_GROUPS + LEB1_COMPANION_GROUPS

# LEB2 distribution partition → ``{tier}_{group}.leb2`` companion files.
# ``uranians`` carries the eight Hamburg-school bodies (IDs 40-47),
# regenerated from the independently sourced Neely (1980) transcription in
# ``libephemeris.hypothetical``. Runtime trust is per-file: a companion is
# attached only when it matches its pinned manifest SHA-256, so pre-existing
# artifacts with the same name remain unused.
LEB2_GROUPS = ("core", "asteroids", "exotics", "apogee", "uranians")
