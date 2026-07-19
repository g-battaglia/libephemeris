# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Explicit Skyfield/DE440 backend for the lunar dev-pipeline test modules.

Every module in this directory verifies the lunar apse model pipeline against
the JPL DE440 kernel (direct ``get_planets()`` reads and JPL-backed readers).
The sealed ``leb`` calculation mode forbids JPL/SPICE access by contract, so
these dev-pipeline suites select their source backend explicitly instead of
inheriting the ambient mode resolution (which resolves to ``leb`` with the
bundled reviewed base core).
"""

from __future__ import annotations

import pytest


@pytest.fixture(autouse=True, scope="module")
def _lunar_skyfield_backend(skyfield_dev_backend: None) -> None:
    """Apply the explicit Skyfield backend to every module in this directory."""
