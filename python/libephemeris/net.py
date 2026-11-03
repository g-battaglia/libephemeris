"""Native net public API."""

from ._core import (
    set_network_policy,
    get_configured_network_policy,
    get_network_policy,
    require_network,
    NetworkSealedError,
)

from libephemeris._core import (
    _module_net_NETWORK_POLICY_ENV as NETWORK_POLICY_ENV,
)
