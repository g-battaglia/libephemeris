from typing import Any

from ._core import (
    _module_exotic_bodies_EXOTIC_IDS as EXOTIC_IDS,
    _module_exotic_bodies_EXOTIC_EXTENDED_IDS as EXOTIC_EXTENDED_IDS,
    _module_exotic_bodies_EXOTIC_ASSIST_PERTURBER_IDS as EXOTIC_ASSIST_PERTURBER_IDS,
)

class ExoticBody:
    body_id: int
    naif: int
    interval_days: float
    degree: int
    target_au: float
    name: str
    cls: str
    def __init__(self, body_id: int, naif: int, interval_days: float, degree: int, target_au: float, name: str, cls: str) -> None: ...
EXOTIC_BODIES: dict[int, ExoticBody]
def naif_map() -> dict[int, int]: ...
def body_params() -> dict[int, tuple[int, int]]: ...
def target_au_map() -> dict[int, float]: ...
def name_map() -> dict[int, str]: ...
def exotic_display_name(body_id: int) -> str: ...
