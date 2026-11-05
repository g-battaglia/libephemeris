# SPDX-License-Identifier: AGPL-3.0
"""Named artifact inventories, using native exotic-body metadata."""

from libephemeris.exotic_bodies import EXOTIC_EXTENDED_IDS, EXOTIC_IDS

LEB2_GROUPS = {
    "core": frozenset([*range(13), 14]),
    "asteroids": frozenset((15, 17, 18, 19, 20)),
    "exotics": frozenset(EXOTIC_IDS),
    "apogee": frozenset((13, 21, 22)),
}


def select_group(group: str, tier: str | None, available: set[int]) -> set[int]:
    """Require the exact declared group; a missing source is never discarded."""
    if tier not in (None, "base", "medium", "extended"):
        raise ValueError(f"unknown tier: {tier!r}")
    expected = frozenset(EXOTIC_EXTENDED_IDS) if group == "exotics" and tier == "extended" else LEB2_GROUPS[group]
    actual = available & LEB2_GROUPS[group]
    if actual != expected:
        raise ValueError(f"source inventory does not match {tier or 'generic'} {group}: expected {sorted(expected)}, found {sorted(actual)}")
    return set(expected)
