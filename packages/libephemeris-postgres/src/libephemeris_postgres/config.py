# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Environment configuration for the PostgreSQL coefficient source."""

from __future__ import annotations

import os
from dataclasses import dataclass


class ProviderConfigError(Exception):
    """Raised when provider configuration is missing or invalid."""


@dataclass(frozen=True, slots=True)
class RuntimeConfig:
    """Runtime settings read from environment variables."""

    dsn: str
    pool_max: int = 2
    timeout_seconds: float = 5.0
    cache_pages: int = 1024


def _positive_int(name: str, default: int) -> int:
    value = os.environ.get(name, str(default))
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ProviderConfigError(f"{name} must be an integer") from exc
    if parsed <= 0:
        raise ProviderConfigError(f"{name} must be positive")
    return parsed


def _positive_float(name: str, default: float) -> float:
    value = os.environ.get(name, str(default))
    try:
        parsed = float(value)
    except ValueError as exc:
        raise ProviderConfigError(f"{name} must be a number") from exc
    if parsed <= 0:
        raise ProviderConfigError(f"{name} must be positive")
    return parsed


def runtime_config() -> RuntimeConfig:
    """Read and validate runtime settings without exposing the DSN."""

    dsn = os.environ.get("LIBEPHEMERIS_PG_URL", "").strip()
    if not dsn:
        raise ProviderConfigError("LIBEPHEMERIS_PG_URL is not configured")
    return RuntimeConfig(
        dsn=dsn,
        pool_max=_positive_int("LIBEPHEMERIS_PG_POOL_MAX", 2),
        timeout_seconds=_positive_float("LIBEPHEMERIS_PG_TIMEOUT_SECONDS", 5.0),
        cache_pages=_positive_int("LIBEPHEMERIS_PG_CACHE_PAGES", 1024),
    )


def admin_dsn(explicit: str | None = None) -> str:
    """Return the importer DSN, never including it in an error message."""

    dsn = (explicit or os.environ.get("LIBEPHEMERIS_PG_ADMIN_URL", "")).strip()
    if not dsn:
        raise ProviderConfigError("PostgreSQL admin URL is not configured")
    return dsn
