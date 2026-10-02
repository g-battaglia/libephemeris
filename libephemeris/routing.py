# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Explicit tier policy and lazy adapters over unchanged coefficient readers.

Provenance:
    Project-authored storage routing using the existing base/medium/extended
    precedence. No polynomial, fit, temporal model or astronomical constants
    are introduced. The selector uses explicit records and closed intervals.
"""

from __future__ import annotations

import os
import hashlib
import threading
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from uuid import UUID

from .exceptions import ConfigurationError, DBDataError, EphemerisRangeError
from .leb_composite import CompositeLEBReader, TieredLEBReader
from .operations import _idle_worker_action, source_guard

TIERS = TieredLEBReader.TIER_PRIORITY


@dataclass(frozen=True, slots=True)
class TierRoute:
    """Immutable storage declaration, with no I/O or scientific behavior.

    Attributes:
        backend: Either leb or db.
        files: Explicit ordered native artifacts for one tier.
        dataset_id: Immutable published UUID for a DB tier.
    """

    backend: str
    files: tuple[str, ...] = ()
    dataset_id: str | None = None


@dataclass(frozen=True, slots=True)
class RouteDecision:
    """Selector result: selected tier or the next metadata requirement."""

    tier: str | None
    needs_metadata: bool = False


def select_tier(
    tiers: tuple[str, ...],
    coverage: Mapping[str, tuple[float, float] | None],
    jd: float,
) -> RouteDecision:
    """Select a closed interval without I/O, ambient state or extrapolation.

    Args:
        tiers: Enabled tiers in existing priority order.
        coverage: Inspected intervals; None means verified absence.
        jd: Exact requested evaluation epoch.

    Returns:
        Selection, first missing metadata requirement, or no coverage.
    """
    for tier in tiers:
        if tier not in coverage:
            return RouteDecision(tier, True)
        bounds = coverage[tier]
        if bounds is not None and bounds[0] <= jd <= bounds[1]:
            return RouteDecision(tier)
    return RouteDecision(None)


_override: tuple[dict[str, TierRoute], str | None] | None = None
_local_readers: dict[tuple[str, tuple[str, ...]], CompositeLEBReader] = {}
_lock = threading.RLock()


def _validate_routes(routes: Mapping[str, Any]) -> dict[str, TierRoute]:
    """Validate and copy declarations without opening any resource.

    Args:
        routes: Tier-name mapping of records or TOML-shaped dictionaries.

    Returns:
        Canonical immutable records in priority order.

    Raises:
        ConfigurationError: A declaration is malformed or ambiguous.
    """
    if not isinstance(routes, Mapping) or not routes or set(routes) - set(TIERS):
        raise ConfigurationError("Routing requires known explicit coefficient tiers")
    result = {}
    for tier in TIERS:
        if tier not in routes:
            continue
        item = routes[tier]
        if isinstance(item, Mapping):
            if set(item) - {"backend", "files", "dataset_id"}:
                raise ConfigurationError("Unknown tier route fields")
            item = TierRoute(
                item.get("backend", ""), item.get("files", ()), item.get("dataset_id")
            )
        if not isinstance(item, TierRoute) or item.backend not in ("leb", "db"):
            raise ConfigurationError("Tier backend must be leb or db")
        if item.backend == "leb":
            if (
                not isinstance(item.files, (tuple, list))
                or not item.files
                or not all(isinstance(p, str) and p.strip() for p in item.files)
                or item.dataset_id is not None
            ):
                raise ConfigurationError("LEB tier requires only an explicit file list")
            files = tuple(os.path.abspath(os.path.expanduser(p)) for p in item.files)
            if len(set(files)) != len(files):
                raise ConfigurationError("Duplicate routed LEB artifact")
            for path in files:
                encoded = set(Path(path).stem.split("_")) & set(TIERS)
                if encoded and encoded != {tier}:
                    raise ConfigurationError("Routed artifact belongs to another tier")
            result[tier] = TierRoute("leb", files)
        else:
            if (
                not isinstance(item.files, (tuple, list))
                or item.files
                or not isinstance(item.dataset_id, str)
            ):
                raise ConfigurationError(
                    "DB tier requires only an explicit dataset UUID"
                )
            try:
                dataset = str(UUID(item.dataset_id))
            except (ValueError, AttributeError):
                raise ConfigurationError("DB dataset must be a valid UUID") from None
            result[tier] = TierRoute("db", dataset_id=dataset)
    datasets = [route.dataset_id for route in result.values() if route.backend == "db"]
    if len(set(datasets)) != len(datasets):
        raise ConfigurationError("Each DB tier requires its own dataset UUID")
    return result


@_idle_worker_action
def set_tier_routes(
    routes: Mapping[str, TierRoute] | None, *, db_url: str | None = None
) -> None:
    """Configure tier storage atomically; mode is never changed implicitly.

    Args:
        routes: Explicit declarations, or None to restore environment/TOML.
        db_url: Optional common PostgreSQL DSN override.

    Raises:
        ConfigurationError: Route records or DSN type are invalid.
    """
    global _override
    validated = None if routes is None else _validate_routes(routes)
    if db_url is not None and (not isinstance(db_url, str) or not db_url.strip()):
        raise ConfigurationError("Routing requires a nonempty PostgreSQL URL")
    if routes is None and db_url is not None:
        raise ConfigurationError("DB URL override requires tier routes")
    with _lock:
        close_routing()
        _override = None if validated is None else (validated, db_url)


def get_tier_routes() -> tuple[dict[str, TierRoute], str | None]:
    """Resolve an operation snapshot with setter → environment → TOML priority.

    Returns:
        Enabled route declarations and redaction-sensitive DSN.

    Raises:
        ConfigurationError: Explicit routed mode lacks usable declarations.
    """
    from ._config_toml import get_all
    from .state import get_precision_tier

    config = get_all()
    if _override is not None:
        routes, override_url = _override
        routes = dict(routes)
    else:
        raw = config.get("tier_routes", {})
        if not isinstance(raw, Mapping):
            raise ConfigurationError("tier_routes must be a table")
        raw = dict(raw)
        for tier in TIERS:
            dataset = os.environ.get(f"LIBEPHEMERIS_DB_DATASET_{tier.upper()}")
            if dataset is not None:
                raw[tier] = {"backend": "db", "dataset_id": dataset}
        routes = _validate_routes(raw)
        override_url = None
    url = override_url or os.environ.get("LIBEPHEMERIS_DB_URL") or config.get("db_url")
    eligible = TIERS[: TIERS.index(get_precision_tier()) + 1]
    routes = {tier: route for tier, route in routes.items() if tier in eligible}
    if not routes:
        raise ConfigurationError("No routed tiers enabled at the configured precision")
    if any(route.backend == "db" for route in routes.values()) and not url:
        raise ConfigurationError("DB tier routes require db_url")
    return routes, url


def manifest_reviewed(manifest: object) -> bool:
    """Compare publication identity with reviewed distribution pins.

    Args:
        manifest: Published source-artifact metadata, not evaluated states.

    Returns:
        Whether every named artifact has its exact reviewed SHA-256 pin.
    """
    from .download import DATA_FILES

    artifacts = manifest.get("artifacts", []) if isinstance(manifest, dict) else []
    return (
        isinstance(artifacts, list)
        and bool(artifacts)
        and all(
            isinstance(artifact, dict)
            and isinstance(artifact.get("name"), str)
            and isinstance(artifact.get("sha256"), str)
            and artifact["sha256"] == DATA_FILES.get(artifact["name"], {}).get("sha256")
            for artifact in artifacts
        )
    )


def manifest_groups(manifest: object, tier: str) -> tuple[str, ...]:
    """Identify only canonical published groups, without exposing raw filenames.

    Args:
        manifest: Source-artifact publication metadata.
        tier: Validated dataset tier.

    Returns:
        Known group identifiers, never connection strings or arbitrary names.
    """
    from .leb_groups import LEB2_GROUPS

    artifacts = manifest.get("artifacts", []) if isinstance(manifest, dict) else []
    if not isinstance(artifacts, list):
        return ()
    names = {
        artifact["name"]
        for artifact in artifacts
        if isinstance(artifact, dict) and isinstance(artifact.get("name"), str)
    }
    return tuple(group for group in LEB2_GROUPS if f"{tier}_{group}.leb2" in names)


@source_guard
def _local_reader(tier: str, route: TierRoute) -> CompositeLEBReader:
    """Open exactly the declared files, never discover companions.

    Args:
        tier: Declared tier identity.
        route: Validated local declaration.

    Returns:
        Process-owned file-only composite.
    """
    from .leb_reader import LEBCorruptionError, open_leb
    from .download import DATA_FILES

    key = (tier, route.files)
    with _lock:
        if key in _local_readers:
            return _local_readers[key]
        readers = []
        try:
            seen: set[int] = set()
            for path in route.files:
                reader = open_leb(path)
                readers.append(reader)
                if seen & set(reader._bodies):
                    raise LEBCorruptionError("Duplicate body in routed artifacts")
                seen.update(reader._bodies)
                pin = DATA_FILES.get(Path(path).name, {}).get("sha256")
                reviewed = False
                if pin:
                    with open(path, "rb") as artifact:
                        reviewed = (
                            hashlib.file_digest(artifact, "sha256").hexdigest() == pin
                        )
                setattr(reader, "_manifest_verified", reviewed)
                setattr(reader, "tier", tier)
            composite = CompositeLEBReader(readers)
        except (OSError, ValueError) as error:
            for reader in readers:
                reader.close()
            raise LEBCorruptionError(
                "Declared routed LEB artifact is unavailable or invalid"
            ) from error
        except Exception:
            for reader in readers:
                reader.close()
            raise
        _local_readers[key] = composite
        return composite


def get_local_leb_reader() -> TieredLEBReader | None:
    """Expose only declared local routes to the unchanged file inventory API.

    Returns:
        File-only tier adapter or None. Never opens a DB pool.
    """
    routes, _url = get_tier_routes()
    local = {
        tier: _local_reader(tier, route)
        for tier, route in routes.items()
        if route.backend == "leb"
    }
    return TieredLEBReader(local) if local else None


@_idle_worker_action
def close_routing() -> None:
    """Close persistent files only without active coefficient owners."""
    from .db.backend import close_db

    with _lock:
        for reader in _local_readers.values():
            reader.close()
        _local_readers.clear()
        # Global caches below are file-only, but their identity keys must not
        # survive closing a file and later reusing its Python object address.
        from . import fast_calc
        from .leb_vector import reset_leb_vector_ephemeris

        fast_calc._reset_active_reader()
        fast_calc._leb_frame_cache.clear()
        reset_leb_vector_ephemeris()
        close_db()


class _BodyEntries:
    """Lazy Python protocol adapter, never a portable numerical record."""

    def __init__(self, reader: RoutedReader) -> None:
        self.reader = reader

    def __getitem__(self, body_id: int) -> Any:
        serving = self.reader.body_reader(body_id, self.reader._jd)
        if serving is None:
            # Coordinate metadata still exists at an uncovered epoch. Only
            # numerical evaluation must fail; do not call a range miss absence.
            serving = self.reader.body_reader(body_id)
        if serving is None:
            raise KeyError(f"Body {body_id} not in any installed LEB tier")
        return serving._bodies[body_id]


class RoutedReader:
    """Operation-owned mixed-source adapter with lazy remote metadata."""

    cacheable = False
    deferred_targets = True
    _manifest_verified = False

    def __init__(self, routes: dict[str, TierRoute], db_url: str | None) -> None:
        """Snapshot validated policy without opening files or transport.

        Args:
            routes: Enabled declarations in tier-priority order.
            db_url: Shared DSN, never included in public diagnostics.
        """
        self.routes = dict(routes)
        self.db_url = db_url
        self._tier_readers: dict[str, Any] = {}
        self._bodies = _BodyEntries(self)
        self._jd: float | None = None
        self._sources: set[str] = set()
        self._closed = False

    @property
    def source(self) -> str:
        """Report sources evaluated in this call, including nested support calls."""
        return "Mixed" if len(self._sources) > 1 else next(iter(self._sources), "LEB")

    @source_guard
    def _reader(self, tier: str) -> Any:
        """Acquire one source lazily and enforce its publication identity.

        Args:
            tier: An enabled policy tier.

        Returns:
            File-only process reader or operation-owned DB reader.
        """
        if self._closed:
            raise DBDataError("Routed operation reader is closed")
        if tier in self._tier_readers:
            return self._tier_readers[tier]
        route = self.routes[tier]
        if route.backend == "leb":
            reader = _local_reader(tier, route)
        else:
            from .db.backend import get_store
            from .db.reader import DBReader

            assert self.db_url is not None and route.dataset_id is not None
            store = get_store(self.db_url)
            declared_tier, manifest = store.dataset_info(route.dataset_id)
            if declared_tier != tier:
                raise DBDataError("DB dataset tier does not match configured route")
            reader = DBReader(store, route.dataset_id)
            reader.tier = tier
            reader._manifest_verified = manifest_reviewed(manifest)
            reader.artifact_groups = manifest_groups(manifest, tier)
            reader._input_budget = self._reserve_inputs
        # Coordinate representations must agree when an iteration crosses a
        # tier boundary; do not compare values or fit a different trajectory.
        for previous in self._tier_readers.values():
            for body in set(previous._bodies) & set(reader._bodies):
                if previous._bodies[body].coord_type != reader._bodies[body].coord_type:
                    if route.backend == "db":
                        reader.close()
                    raise DBDataError(
                        "Routed coordinate channels disagree across tiers"
                    )
        self._tier_readers[tier] = reader
        return reader

    def _reserve_inputs(self, incoming: int) -> None:
        """Evict transient inputs before crossing the aggregate operation budget.

        Args:
            incoming: Number of new remote segments about to be fetched.
        """
        readers = [
            r
            for r in self._tier_readers.values()
            if getattr(r, "source", "LEB") == "DB"
        ]
        if sum(len(r._segments) for r in readers) + incoming > 256:
            for reader in readers:
                reader._segments.clear()

    @source_guard
    def selected_body_reader(self, body_id: int, jd: float) -> Any | None:
        """Resolve only metadata needed for one exact evaluation epoch.

        Args:
            body_id: Stored channel identifier.
            jd: Julian day TT, including retarded/speed support dates.

        Returns:
            Concrete covering reader, or None after verified lack of coverage.
        """
        coverage: dict[str, tuple[float, float] | None] = {}
        tiers = tuple(self.routes)
        while True:
            decision = select_tier(tiers, coverage, jd)
            if not decision.needs_metadata:
                if decision.tier is None:
                    return None
                reader = self._reader(decision.tier)
                body_reader = getattr(reader, "body_reader", None)
                return body_reader(body_id) if body_reader else reader
            assert decision.tier is not None
            reader = self._reader(decision.tier)
            coverage[decision.tier] = reader.body_coverage(body_id)

    def body_reader(self, body_id: int, jd: float | None = None) -> Any | None:
        """Return exact-date selection or the native priority presence fallback.

        Args:
            body_id: Stored channel identifier.
            jd: Optional Julian day TT; omit only for metadata/presence queries.

        Returns:
            Concrete reader or None for a missing channel.
        """
        if jd is not None:
            return self.selected_body_reader(body_id, jd)
        for tier in self.routes:
            reader = self._reader(tier)
            if reader.has_body(body_id):
                body_reader = getattr(reader, "body_reader", None)
                return body_reader(body_id) if body_reader else reader
        return None

    def has_body(self, body_id: int) -> bool:
        """Preserve native presence semantics while preferring the prepared date.

        Args:
            body_id: Stored channel identifier.

        Returns:
            Whether any enabled source declares the channel.
        """
        if self._jd is not None:
            selected = self.selected_body_reader(body_id, self._jd)
            if selected is not None:
                return True
        return self.body_reader(body_id) is not None

    def body_coverage(self, body_id: int) -> tuple[float, float] | None:
        """Inspect an all-tier envelope; unlike evaluation this can open DB tiers.

        Args:
            body_id: Stored channel identifier.

        Returns:
            Native union envelope, or None for verified absence.
        """
        ranges = [
            bounds
            for tier in self.routes
            if (bounds := self._reader(tier).body_coverage(body_id)) is not None
        ]
        return (
            (min(b[0] for b in ranges), max(b[1] for b in ranges)) if ranges else None
        )

    @property
    def jd_range(self) -> tuple[float, float]:
        """Inspect the complete configured envelope, never for lazy dispatch."""
        ranges = [self._reader(tier).jd_range for tier in self.routes]
        return min(b[0] for b in ranges), max(b[1] for b in ranges)

    def has_nutation(self) -> bool:
        """Return a capability hint without opening an unknown remote tier."""
        return any(reader.has_nutation() for reader in self._tier_readers.values())

    @source_guard
    def frame_reader(self, jd_tt: float) -> Any | None:
        """Resolve the exact auxiliary source, allowing file-only frame caches.

        Args:
            jd_tt: Frame epoch in Julian days TT.

        Returns:
            Concrete nutation reader, or None when the optional series is absent.
        """
        from .db.contract import NUTATION_ID

        for tier in self.routes:
            reader = self._reader(tier)
            start, end = reader.jd_range
            if not start <= jd_tt <= end:
                continue
            serving = getattr(reader, "_nutation_reader", reader)
            if serving is None or not serving.has_nutation():
                continue
            series = (
                serving._series.get(NUTATION_ID)
                if hasattr(serving, "_series")
                else serving._nutation
            )
            if series is not None and series.jd_start <= jd_tt <= series.jd_end:
                self._sources.add(getattr(serving, "source", "LEB"))
                return serving
        return None

    @source_guard
    def prepare(self, jd_tt: float, body_id: int, flags: int) -> None:
        """Prefetch one target dataset without pinning its support epochs.

        Args:
            jd_tt: Target Julian day TT.
            body_id: Public body identifier.
            flags: Existing calculation flags, passed unchanged to read planning.
        """
        self._jd = jd_tt
        from .constants import (
            MEAN_NODE,
            MEAN_APOG,
            INTP_APOG,
            INTP_PERG,
        )

        # These points use existing local models rather than stored channels.
        # Their support states are still selected at exact evaluation epochs.
        if (
            body_id in (MEAN_NODE, MEAN_APOG, INTP_APOG, INTP_PERG)
            or 40 <= body_id <= 58
        ):
            return
        # Only prefetch the selected target's dataset; exact support and
        # retarded epochs are resolved independently, just as in TieredLEBReader.
        reader = self.selected_body_reader(body_id, jd_tt)
        prepare = getattr(reader, "prepare", None)
        if prepare is not None:
            prepare(jd_tt, body_id, flags)

    @source_guard
    def eval_body(self, body_id: int, jd: float) -> tuple[tuple, tuple]:
        """Delegate unchanged native evaluation at the exact requested epoch.

        Args:
            body_id: Stored channel identifier.
            jd: Julian day TT.

        Returns:
            Native three-component position and velocity.
        """
        reader = self.selected_body_reader(body_id, jd)
        if reader is None:
            bounds = self.body_coverage(body_id)
            if bounds is None:
                raise KeyError(f"Body {body_id} not in any installed LEB tier")
            raise EphemerisRangeError(
                f"No routed coefficient interval covers body {body_id} at JD {jd}",
                requested_jd=jd,
                start_jd=float(bounds[0]),
                end_jd=float(bounds[1]),
                body_id=body_id,
            )
        result = reader.eval_body(body_id, jd)
        self._sources.add(getattr(reader, "source", "LEB"))
        return result

    @source_guard
    def eval_nutation(self, jd_tt: float) -> tuple[float, float]:
        """Evaluate the selected optional nutation channel without extrapolation.

        Args:
            jd_tt: Julian day TT.

        Returns:
            Native nutation angles in radians.
        """
        reader = self.frame_reader(jd_tt)
        if reader is None:
            raise ValueError(f"No LEB nutation data covers JD {jd_tt}")
        return reader.eval_nutation(jd_tt)

    @source_guard
    def delta_t(self, jd: float) -> float:
        """Expose stored auxiliary data without overriding temporal-model policy.

        Args:
            jd: Requested native interpolation epoch.

        Returns:
            Native Delta-T value in days.
        """
        for tier in self.routes:
            reader = self._reader(tier)
            start, end = reader.jd_range
            if start <= jd <= end:
                try:
                    value = reader.delta_t(jd)
                    self._sources.add(getattr(reader, "source", "LEB"))
                    return value
                except ValueError:
                    continue
        raise ValueError(f"No LEB Delta-T data covers JD {jd}")

    @source_guard
    def get_star(self, star_id: int) -> Any:
        """Look up stored auxiliary stars without changing catalog precedence.

        Args:
            star_id: Native catalog identifier.

        Returns:
            Detached native star record.
        """
        from .exceptions import StarNotFoundError

        for tier in self.routes:
            try:
                reader = self._reader(tier)
                star = reader.get_star(star_id)
                self._sources.add(getattr(reader, "source", "LEB"))
                return star
            except (KeyError, StarNotFoundError):
                continue
        raise KeyError(f"Star {star_id} not in any installed LEB tier")

    def close(self) -> None:
        """Discard remote scientific inputs while preserving local file handles."""
        for tier, reader in self._tier_readers.items():
            if self.routes[tier].backend == "db":
                reader.close()
        self._tier_readers.clear()
        self._sources.clear()
        self._closed = True
