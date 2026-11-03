from pathlib import Path
from typing import Any, Optional, Union, Dict, Literal
import os
from .exceptions import NetworkSealedError as NetworkSealedError

NetworkPolicy = Literal["auto", "allow", "sealed"]

def set_network_policy(policy: NetworkPolicy | None) -> None: ...
def get_configured_network_policy() -> NetworkPolicy: ...
def get_network_policy() -> Literal['allow', 'sealed']: ...
def require_network(purpose: str) -> None: ...

from libephemeris._core import (
    _module_net_NETWORK_POLICY_ENV as NETWORK_POLICY_ENV,
)
