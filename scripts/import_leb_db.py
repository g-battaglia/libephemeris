#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Thin repository entry point for the dedicated DB provisioning package.

Provenance:
    Project-authored CLI delegation; no data or scientific logic. Schema and
    coefficient handling remain in the dedicated libephemeris.db package.
"""

from __future__ import annotations

from libephemeris.db.__main__ import main

if __name__ == "__main__":
    main()
