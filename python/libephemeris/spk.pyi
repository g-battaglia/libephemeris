from pathlib import Path
from typing import Any, Optional, Union, Dict, Literal
import os

def download_spk(body: str, start: str, end: str, path: Optional[str]=None, center: str='500@0', overwrite: bool=False, timeout: int=120) -> str: ...
def register_spk_body(ipl: int, spk_file: str, naif_id: Union[int, str]) -> None: ...
def unregister_spk_body(ipl: int) -> None: ...
def get_spk_body_info(ipl: int) -> Optional[tuple[str, int]]: ...
def list_spk_bodies() -> dict[int, tuple[str, int]]: ...
def get_spk_coverage(spk_file: str) -> Optional[tuple[float, float]]: ...
def download_and_register_spk(body: str, ipl: int, start: str, end: str, naif_id: Optional[int]=None, path: Optional[str]=None, center: str='500@0', overwrite: bool=False, timeout: int=120) -> str: ...

from libephemeris._core import (
    SPKNotFoundError,
    format_file_size,
    get_logger,
    require_network,
)

def calc_spk_body_position(t: float, ipl: int, iflag: int) -> tuple[float,float,float,float,float,float] | None: ...

from libephemeris._core import (
    get_timescale,
)

from libephemeris._core import (
    _module_spk_NAIF_ASTEROID_OFFSET as NAIF_ASTEROID_OFFSET,
    _module_spk_NAIF_ASTEROID_OFFSET_HORIZONS as NAIF_ASTEROID_OFFSET_HORIZONS,
    _module_spk_HORIZONS_API_URL as HORIZONS_API_URL,
)
