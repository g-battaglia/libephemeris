# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Optional PostgreSQL coefficient backend and offline provisioning tools.

Provenance:
    Project-authored storage integration. All DB-specific runtime and import
    logic lives in this package, independently of the astronomical engine.
    This entry point defines no astronomical model or coefficient.
"""

from __future__ import annotations

from .backend import set_db_config

__all__ = ["set_db_config"]
