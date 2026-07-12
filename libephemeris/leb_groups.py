# SPDX-License-Identifier: Apache-2.0
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
"""

from __future__ import annotations

# LEB1 generation/merge partition (order = generation + merge order).
LEB1_GROUPS = ("planets", "asteroids", "exotics", "analytical")

# LEB2 distribution partition → ``{tier}_{group}.leb2`` companion files.
LEB2_GROUPS = ("core", "asteroids", "exotics", "apogee", "uranians")
