# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Process-wide network policy and the single runtime URL opener.

``mode = "leb"`` is intended to be an offline, deterministic calculation
mode.  Historically only the live Horizons backend honoured that contract;
SPK generation, SBDB lookups, IERS refreshes and Skyfield's implicit kernel
download could still leave the process.  This module is the choke point that
makes the contract explicit and testable.

The default policy is ``"auto"``: LEB mode is sealed and every other mode is
allowed.  ``"sealed"`` and ``"allow"`` are explicit process overrides.  The
``libephemeris download ...`` CLI deliberately selects ``"allow"`` because a
download command is an explicit provisioning action, not calculation-time
egress.

Provenance:
    Project-authored network-policy infrastructure.  It selects whether an
    HTTPS operation may start and does not alter astronomical data or models.
"""

from __future__ import annotations

import os
import http.client
import urllib.error
import urllib.parse
import urllib.request
from typing import Any, Literal

from skyfield.api import Loader as _SkyfieldLoader

from .exceptions import NetworkSealedError

NetworkPolicy = Literal["auto", "allow", "sealed"]

# Re-export the small request/error/encoding surface used by network clients so
# no other package module needs to import a networking stdlib module directly.
# Keeping construction as well as opening behind this module makes the static
# sealed-mode invariant mechanically auditable.
HTTPException = http.client.HTTPException
HTTPError = urllib.error.HTTPError
URLError = urllib.error.URLError
Request = urllib.request.Request
quote = urllib.parse.quote
urlencode = urllib.parse.urlencode

NETWORK_POLICY_ENV = "LIBEPHEMERIS_NETWORK_POLICY"
_VALID_NETWORK_POLICIES = ("auto", "allow", "sealed")
_NETWORK_POLICY: NetworkPolicy | None = None


class PolicyAwareLoader(_SkyfieldLoader):
    """Skyfield loader whose implicit download path obeys this policy.

    ``Loader.__call__`` invokes ``download()`` only when its requested local
    asset is absent. Overriding that method seals both the module API and
    :class:`EphemerisContext`, even if a future caller omits a preflight check.
    """

    def _assure(
        self,
        url: str | None,
        filename: str,
        reload: bool,
        backup: bool,
    ) -> str:
        """Gate Skyfield's internal ``Loader.__call__`` download branch."""
        if reload or not os.path.exists(self.path_to(filename)):
            require_network(f"Skyfield data download: {filename}")
        return super()._assure(url, filename, reload, backup)

    def download(
        self,
        url: str,
        filename: str | None = None,
        backup: bool = False,
    ) -> str:
        """Download through Skyfield only when network policy permits it."""
        label = filename or url.rsplit("/", 1)[-1]
        require_network(f"Skyfield data download: {label}")
        return super().download(url, filename=filename, backup=backup)


def set_network_policy(policy: NetworkPolicy | None) -> None:
    """Set the process network policy, or reset to env/TOML with ``None``."""
    global _NETWORK_POLICY
    if policy is not None:
        normalized = policy.lower().strip()
        if normalized not in _VALID_NETWORK_POLICIES:
            raise ValueError(
                f"Invalid network policy: {policy!r}. Must be one of: "
                f"{list(_VALID_NETWORK_POLICIES)}"
            )
        policy = normalized  # type: ignore[assignment]
    _NETWORK_POLICY = policy


def get_configured_network_policy() -> NetworkPolicy:
    """Return the configured policy before resolving ``"auto"``."""
    if _NETWORK_POLICY is not None:
        return _NETWORK_POLICY

    env_value = os.environ.get(NETWORK_POLICY_ENV, "").lower().strip()
    if env_value in _VALID_NETWORK_POLICIES:
        return env_value  # type: ignore[return-value]

    from ._config_toml import get_str

    toml_value = (get_str("network_policy") or "").lower().strip()
    if toml_value in _VALID_NETWORK_POLICIES:
        return toml_value  # type: ignore[return-value]
    return "auto"


def get_network_policy() -> Literal["allow", "sealed"]:
    """Return the effective policy for the current calculation mode."""
    configured = get_configured_network_policy()
    if configured != "auto":
        return configured

    # Lazy import avoids a state -> net -> state cycle.  Reading the mode is
    # lock-free by design in state.py and is safe from network call sites.
    from .state import get_calc_mode

    return "sealed" if get_calc_mode() == "leb" else "allow"


def require_network(purpose: str) -> None:
    """Raise before any socket work when the effective policy is sealed."""
    if get_network_policy() == "sealed":
        raise NetworkSealedError(
            message=(
                f"Network access is sealed; blocked operation: {purpose}. "
                "Provision required files explicitly with "
                "'libephemeris download ...' before starting the service."
            ),
            purpose=purpose,
        )


def open_url(
    url_or_request: Any,
    *,
    purpose: str,
    timeout: float | None = None,
    context: Any = None,
) -> Any:
    """Policy-gated equivalent of :func:`urllib.request.urlopen`."""
    require_network(purpose)
    kwargs: dict[str, Any] = {}
    if timeout is not None:
        kwargs["timeout"] = timeout
    if context is not None:
        kwargs["context"] = context
    return urllib.request.urlopen(url_or_request, **kwargs)


__all__ = [
    "NETWORK_POLICY_ENV",
    "HTTPError",
    "HTTPException",
    "NetworkPolicy",
    "PolicyAwareLoader",
    "Request",
    "URLError",
    "get_configured_network_policy",
    "get_network_policy",
    "open_url",
    "quote",
    "require_network",
    "set_network_policy",
    "urlencode",
]
