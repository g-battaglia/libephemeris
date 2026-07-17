# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Machine-readable LEB coverage and active-data inventory.

Availability is a property of a body *and a date*. A core-level readiness
boolean cannot represent modular files whose companions and per-body ranges
differ, so this module exposes the ranges stored in active reader headers.

Provenance:
    Project-authored introspection over already-validated LEB metadata. It
    reports file identities and ranges without changing astronomical values.
"""

from __future__ import annotations

import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True, slots=True)
class BodyCoverage:
    """Coverage contract for one body from the active calculation source."""

    body_id: int
    source: str
    precision_class: str
    jd_start: float
    jd_end: float
    data_file: str | None
    group: str | None
    reviewed: bool

    def contains(self, jd: float) -> bool:
        """Return whether ``jd`` is inside this body's closed interval."""
        return self.jd_start <= float(jd) <= self.jd_end

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable representation."""
        return asdict(self)


@dataclass(frozen=True, slots=True)
class RuntimeDataRequirement:
    """One immutable file required by a sealed LEB runtime tier."""

    name: str
    kind: str
    group: str | None
    path: str
    sha256: str

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable representation."""
        return asdict(self)


def get_runtime_data_requirements(
    tier: str | None = None,
) -> tuple[RuntimeDataRequirement, ...]:
    """Return the manifest-derived file contract for sealed LEB operation.

    The result follows the canonical LEB2 partition and the reviewed download
    manifest. A missing pin is a packaging error and raises immediately instead
    of producing a partial readiness contract.
    """
    from .download import DATA_FILES, get_data_dir
    from .leb_groups import LEB2_GROUPS
    from .state import get_precision_tier

    active_tier = tier or get_precision_tier()
    if active_tier not in ("base", "medium", "extended"):
        raise ValueError(f"Unsupported precision tier: {active_tier!r}")

    data_dir = get_data_dir()
    requirements: list[RuntimeDataRequirement] = []
    for group in LEB2_GROUPS:
        name = f"{active_tier}_{group}.leb2"
        info = DATA_FILES.get(name)
        sha256 = info.get("sha256") if info is not None else None
        if not isinstance(sha256, str) or not sha256:
            raise RuntimeError(f"Reviewed manifest pin missing for {name}")
        requirements.append(
            RuntimeDataRequirement(
                name=name,
                kind="leb2",
                group=group,
                path=str(data_dir / "leb" / name),
                sha256=sha256,
            )
        )

    centers_name = f"planet_centers_{active_tier}.bsp"
    centers_info = DATA_FILES.get(centers_name)
    centers_sha256 = centers_info.get("sha256") if centers_info is not None else None
    if not isinstance(centers_sha256, str) or not centers_sha256:
        raise RuntimeError(f"Reviewed manifest pin missing for {centers_name}")
    requirements.append(
        RuntimeDataRequirement(
            name=centers_name,
            kind="planet_centers",
            group=None,
            path=str(data_dir / centers_name),
            sha256=centers_sha256,
        )
    )
    return tuple(requirements)


def _group_from_path(path: str | None) -> str | None:
    """Infer a canonical LEB group name from a manifest-style file path."""
    if not path:
        return None
    stem = Path(path).stem
    if "_" not in stem:
        return None
    return stem.rsplit("_", 1)[-1]


def _serving_reader(reader: Any, body_id: int) -> Any | None:
    """Return the concrete reader serving ``body_id`` from a composite."""
    body_reader = getattr(reader, "body_reader", None)
    if body_reader is not None:
        return body_reader(body_id)
    return reader if getattr(reader, "has_body", lambda _body: False)(body_id) else None


def get_body_coverage(body_id: int) -> BodyCoverage | None:
    """Return active LEB coverage for ``body_id``, independent of target date.

    ``None`` means that the active LEB reader does not contain the body. It
    does not imply that an analytical or online fallback exists. Callers in
    sealed LEB mode should treat a miss for an optional minor body as an
    unavailable capability rather than probing a fallback.
    """
    from .state import get_leb_reader

    try:
        reader = get_leb_reader()
    except RuntimeError:
        return None
    if reader is None:
        return None

    serving = _serving_reader(reader, int(body_id))
    if serving is None:
        return None
    coverage_fn = getattr(serving, "body_coverage", None)
    if coverage_fn is None:
        return None
    bounds = coverage_fn(int(body_id))
    if bounds is None:
        return None

    path = getattr(serving, "path", None)
    reviewed = bool(
        getattr(serving, "_manifest_verified", False)
        or getattr(reader, "_manifest_verified", False)
    )
    path_str = str(path) if path is not None else None
    return BodyCoverage(
        body_id=int(body_id),
        source="LEB",
        precision_class="ephemeris" if reviewed else "unverified-local",
        jd_start=float(bounds[0]),
        jd_end=float(bounds[1]),
        data_file=path_str,
        group=_group_from_path(path_str),
        reviewed=reviewed,
    )


def coverage(body_id: int) -> BodyCoverage | None:
    """Concise alias for :func:`get_body_coverage`."""
    return get_body_coverage(body_id)


def _reader_file_info(reader: Any, *, inherited_verified: bool) -> dict[str, Any]:
    """Serialize one reader's file identity and body coverage metadata."""
    path = str(getattr(reader, "path", ""))
    entries = getattr(reader, "_bodies", {})
    bodies = []
    for body_id in sorted(entries):
        bounds = getattr(reader, "body_coverage", lambda _body: None)(body_id)
        if bounds is None:
            continue
        bodies.append(
            {
                "body_id": int(body_id),
                "jd_start": float(bounds[0]),
                "jd_end": float(bounds[1]),
            }
        )
    try:
        size_bytes = os.path.getsize(path)
    except OSError:
        size_bytes = None
    return {
        "name": os.path.basename(path),
        "path": path,
        "group": _group_from_path(path),
        "size_bytes": size_bytes,
        "reviewed": bool(
            getattr(reader, "_manifest_verified", False) or inherited_verified
        ),
        "body_count": len(bodies),
        "bodies": bodies,
    }


def get_leb_inventory() -> dict[str, Any]:
    """Return active LEB files, per-body ranges, mode and network policy."""
    from .net import get_configured_network_policy, get_network_policy
    from .state import get_calc_mode, get_leb_reader, get_precision_tier

    result: dict[str, Any] = {
        "mode": get_calc_mode(),
        "precision_tier": get_precision_tier(),
        "network_policy_configured": get_configured_network_policy(),
        "network_policy_effective": get_network_policy(),
        "ready": False,
        "reader_type": None,
        "files": [],
        "body_count": 0,
    }
    try:
        reader = get_leb_reader()
    except RuntimeError as exc:
        result["error"] = str(exc)
        return result
    if reader is None:
        return result

    inherited_verified = bool(getattr(reader, "_manifest_verified", False))
    readers = list(getattr(reader, "_readers", (reader,)))
    files = [
        _reader_file_info(item, inherited_verified=inherited_verified)
        for item in readers
    ]
    result.update(
        {
            "ready": True,
            "reader_type": type(reader).__name__,
            "files": files,
            "body_count": len(getattr(reader, "_bodies", {})),
        }
    )
    return result


def inspect_leb_file(path: str | os.PathLike[str]) -> dict[str, Any]:
    """Open one LEB file and return its body/range metadata."""
    from .leb_reader import open_leb

    reader = open_leb(os.fspath(path))
    try:
        return _reader_file_info(reader, inherited_verified=False)
    finally:
        reader.close()


__all__ = [
    "BodyCoverage",
    "RuntimeDataRequirement",
    "coverage",
    "get_body_coverage",
    "get_leb_inventory",
    "get_runtime_data_requirements",
    "inspect_leb_file",
]
