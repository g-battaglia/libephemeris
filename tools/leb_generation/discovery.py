# SPDX-License-Identifier: AGPL-3.0
"""Select existing generation images from native identities and chain coverage.

This is offline input assembly. Rust reads descriptors, composes geometric
sources and converts their coverage to the TT fitting grid. No filename or
file-size heuristic decides which body or physical interval an image serves.
"""

import hashlib
import math
from dataclasses import dataclass
from pathlib import Path

from .source import continuous_overlap


@dataclass(frozen=True)
class LocalSpk:
    """One retained candidate's actual NAIF target, continuous domains and bytes."""

    path: Path
    naif_id: int
    windows_tt: tuple[tuple[float, float], ...]
    seed_at_j2000: bool
    sha256: str


def target_aliases(body_id):
    """Native body metadata admits the two documented asteroid NAIF offsets."""
    from libephemeris.constants import get_spk_body_info_from_map
    from libephemeris.spk import NAIF_ASTEROID_OFFSET, NAIF_ASTEROID_OFFSET_HORIZONS

    identity = get_spk_body_info_from_map(body_id)
    if identity is None:
        raise ValueError(f"body {body_id}: no native SPK identity for automatic generation")
    _, naif = identity
    offsets = sorted((NAIF_ASTEROID_OFFSET, NAIF_ASTEROID_OFFSET_HORIZONS))
    if not offsets[0] <= naif < offsets[1]:
        raise ValueError(f"body {body_id}: automatic selection requires an asteroid NAIF identity")
    number = naif - offsets[0]
    return frozenset(offset + number for offset in offsets)


def scan_local_spks(source, body_ids, directories):
    """Read explicit directories once; preserve complete native center chains.

    Files are deduplicated by canonical path. Only a DAF/SPK header identifies
    an SPK: force binaries can also carry a .bsp suffix. Other headers are
    reported as unused cache files. Unreadable files or malformed identified
    SPK metadata are errors. Every candidate is hashed before its
    metadata is used; attaching the final selection checks that same identity.
    """
    import libephemeris as engine

    aliases = {body: target_aliases(body) for body in body_ids}
    result = {body: [] for body in body_ids}
    paths = set()
    for directory in directories:
        directory = Path(directory).resolve(strict=True)
        if not directory.is_dir():
            raise ValueError(f"SPK search path is not a directory: {directory}")
        paths.update(path.resolve(strict=True) for path in directory.glob("*.bsp") if path.is_file())
    for path in sorted(paths):
        with path.open("rb") as stream:
            if stream.read(8) != b"DAF/SPK ":
                source.ignored_non_spk_files.add(path)
                continue
            stream.seek(0)
            sha256 = hashlib.file_digest(stream, "sha256").hexdigest()
        image = engine.get_loader()(str(path))
        descriptors = image.segment_descriptors_tdb()
        targets = {row[0] for row in descriptors}
        matching = {body: identifiers & targets for body, identifiers in aliases.items()}
        if not any(matching.values()):
            continue
        for body, identifiers in matching.items():
            for target in sorted(identifiers):
                geometric = image.geometric_source(target, [source.image])
                windows = tuple(source.clock.tdb_interval_to_tt(low, high)
                                for low, high in geometric.coverage_tdb())
                seed = any(identifier == target and center == 10 and frame == 1 and low <= 2451545.0 <= high
                           for identifier, center, frame, low, high in descriptors)
                result[body].append(LocalSpk(path, target, windows, seed, sha256))
    return result


def select_local_spk(body, candidates, requested, *, available_domains=False, nbody=False):
    """Select one dominant continuous coverage; ambiguity requires --spk.

    Coverage is an input-selection policy, never a claim that a wider image is
    more accurate. Actual source hashes are printed and used for comparisons.
    For integration seeds, only the direct heliocentric J2000 state matters;
    distinct seed images require an explicit selection rather than a width rank.
    """
    if not (all(math.isfinite(value) for value in requested) and requested[0] < requested[1]):
        raise ValueError("requested generation domain must be finite with positive width")
    eligible = []
    for candidate in candidates:
        if nbody:
            if candidate.seed_at_j2000:
                eligible.append((candidate, None))
            continue
        try:
            overlap = continuous_overlap(candidate.windows_tt, *requested)
        except ValueError:
            # Metadata says this image cannot serve one continuous body record.
            # No failed numerical read is converted into another source choice.
            continue
        if not available_domains and overlap != requested:
            continue
        window = next((low, high) for low, high in candidate.windows_tt
                      if low <= overlap[0] and overlap[1] <= high)
        eligible.append((candidate, window))
    # Identical source bytes and target do not become an ambiguity merely
    # because a cache contains multiple copies of the same file.
    unique = {}
    for candidate, window in eligible:
        unique.setdefault((candidate.sha256, candidate.naif_id), (candidate, window))
    eligible = list(unique.values())
    if not eligible:
        kind = "heliocentric J2000 seed" if nbody else "continuous requested coverage"
        raise ValueError(f"body {body}: no local SPK has {kind}; supply an explicit --spk source")
    if nbody:
        winners = eligible
    else:
        winners = [(candidate, window) for candidate, window in eligible
                   if all(window[0] <= other[0] and other[1] <= window[1] for _, other in eligible)]
    if len(winners) != 1:
        paths = ", ".join(str(candidate.path) for candidate, _ in eligible)
        raise ValueError(f"body {body}: ambiguous local SPK sources ({paths}); select one with --spk")
    return winners[0][0]


def automatic_nbody_ids(body_ids, tier):
    """Retain the declared extended-tier policy, including perturber exclusions."""
    if tier != "extended":
        return set()
    from libephemeris.exotic_bodies import (
        EXOTIC_ASSIST_PERTURBER_IDS,
        EXOTIC_EXTENDED_IDS,
    )

    return set(body_ids) & (set(EXOTIC_EXTENDED_IDS) - set(EXOTIC_ASSIST_PERTURBER_IDS))


def assemble_local_sources(source, args, body_ids, domains, *, available_domains=False):
    """Complete missing native inputs; explicit --spk selections retain priority."""
    import libephemeris as engine

    from .options import resolve_body_tokens

    nbody_ids = (set(resolve_body_tokens(args.nbody, engine)) if args.nbody is not None
                 else automatic_nbody_ids(body_ids, args.tier))
    if not nbody_ids <= set(body_ids):
        raise ValueError("--nbody IDs must be present in the requested output inventory")
    missing = set(body_ids) - source.targets.keys() - source.angular.keys() - source.spk_sources.keys()
    if missing:
        root = args.data_dir or source.kernel.parent
        directories = args.spk_dir or [root / "spk"]
        candidates = scan_local_spks(source, sorted(missing), directories)
        for body in sorted(missing):
            selection = select_local_spk(body, candidates[body], domains[body],
                                         available_domains=available_domains,
                                         nbody=body in nbody_ids)
            source.add_spk_source(body, selection.path, selection.naif_id)
            if source.spk_inputs[body]["images"][str(selection.path)] != selection.sha256:
                raise ValueError(f"body {body}: selected SPK changed during discovery")
    if nbody_ids:
        root = args.data_dir or source.kernel.parent
        name = "linux_m13000p17000.441" if args.tier == "extended" else "linux_p1550p2650.440"
        planets = args.assist_planets or root / "assist" / name
        asteroids = args.assist_asteroids or root / "assist" / "sb441-n16.bsp"
        source.add_nbody_sources(sorted(nbody_ids), planets, asteroids)
    elif args.assist_planets or args.assist_asteroids:
        raise ValueError("ASSIST force files require a requested N-body body")
