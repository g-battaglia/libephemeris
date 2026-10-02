# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Narrow database capability, independent of HTTP/download permissions.

Provenance:
    Project-authored transport policy, not an astronomical algorithm. It only
    authorizes a configured connection and contains no scientific coefficient.
"""

from __future__ import annotations

from ..exceptions import NetworkSealedError


def require_database() -> None:
    """Authorize the configured database without authorizing HTTP downloads.

    Automatic policy permits DB transport in DB mode, while LEB remains
    completely offline. Explicit provisioning can opt into the allow policy;
    an explicit sealed policy always denies PostgreSQL connections as well.

    Raises:
        NetworkSealedError: Transport is explicitly sealed or LEB-only.
    """
    from ..net import get_configured_network_policy
    from ..state import get_calc_mode

    configured = get_configured_network_policy()
    if configured == "sealed" or (configured == "auto" and get_calc_mode() == "leb"):
        raise NetworkSealedError(
            message="Network access is sealed; PostgreSQL connection blocked",
            purpose="PostgreSQL coefficient access",
        )
