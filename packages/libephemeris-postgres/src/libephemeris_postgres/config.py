# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Credential-safe PostgreSQL runtime settings."""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field

from libephemeris import CoefficientSourceError


@dataclass(frozen=True, slots=True)
class RuntimeConfig:
    """Bound the process pool, transport timeout and per-file byte cache."""

    dsn: str = field(repr=False)
    pool_max: int = 2
    timeout_seconds: float = 5.0
    cache_blocks: int = 256

    def __post_init__(self) -> None:
        if not self.dsn or self.pool_max <= 0 or self.cache_blocks <= 0:
            raise CoefficientSourceError("Missing or invalid PostgreSQL configuration")
        if not math.isfinite(self.timeout_seconds) or self.timeout_seconds <= 0:
            raise CoefficientSourceError(
                "PostgreSQL timeout must be finite and positive"
            )


def runtime_config() -> RuntimeConfig:
    """Read explicit bounds, suppressing malformed values in diagnostics."""
    try:
        return RuntimeConfig(
            os.environ.get("LIBEPHEMERIS_PG_URL", "").strip(),
            int(os.environ.get("LIBEPHEMERIS_PG_POOL_MAX", "2")),
            float(os.environ.get("LIBEPHEMERIS_PG_TIMEOUT_SECONDS", "5")),
            int(os.environ.get("LIBEPHEMERIS_PG_CACHE_BLOCKS", "256")),
        )
    except ValueError:
        raise CoefficientSourceError("Invalid PostgreSQL configuration") from None


def admin_dsn(explicit: str | None = None) -> str:
    """Use administrator credentials only for provisioning."""
    dsn = (explicit or os.environ.get("LIBEPHEMERIS_PG_ADMIN_URL", "")).strip()
    if not dsn:
        raise CoefficientSourceError("PostgreSQL admin URL is not configured")
    return dsn
