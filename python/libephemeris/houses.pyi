from collections.abc import Sequence
from typing import overload

def get_polar_latitude_threshold(obliquity: float = 23.44) -> float: ...
def get_extreme_latitude_info(lat: float, obliquity: float = 23.44) -> dict: ...
def houses(tjdut: float, lat: float, lon: float, hsys: int | str | bytes = 80, iflag: int = 0) -> tuple[tuple[float, ...], tuple[float, ...]]: ...
def houses_ex(tjdut: float, lat: float, lon: float, hsys: int | str | bytes = 80, flags: int = 0) -> tuple[tuple[float, ...], tuple[float, ...]]: ...
def houses_ex2(tjdut: float, lat: float, lon: float, hsys: int | str | bytes = 80, flags: int = 0) -> tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...], tuple[float, ...]]: ...
def houses_armc(armc: float, lat: float, eps: float, hsys: int | str | bytes = 80, ascmc9: float = 0.0) -> tuple[tuple[float, ...], tuple[float, ...]]: ...
def houses_armc_ex2(armc: float, lat: float, eps: float, hsys: int | str | bytes = 80, ascmc9: float = 0.0) -> tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...], tuple[float, ...]]: ...
def houses_with_fallback(jd_ut: float, lat: float, lon: float, hsys: int | str | bytes, fallback_hsys: int | str | bytes = 79, validate_cusps: bool = True) -> tuple[tuple[float, ...], tuple[float, ...], bool, str | None]: ...
def houses_armc_with_fallback(armc: float, lat: float, eps: float, hsys: int | str | bytes, fallback_hsys: int | str | bytes = 79, validate_cusps: bool = True, ascmc9: float = 0.0) -> tuple[tuple[float, ...], tuple[float, ...], bool, str | None]: ...
def house_name(hsys: int | str | bytes) -> str: ...
@overload
def house_pos(armc: float, lat: float, obliquity: float, hsys_or_objcoord: Sequence[float], lon_or_hsys: int | str | bytes | None = None, lat_body: float = 0.0) -> float: ...
@overload
def house_pos(armc: float, lat: float, obliquity: float, hsys_or_objcoord: int | str | bytes, lon_or_hsys: float | Sequence[float] | None = None, lat_body: float = 0.0) -> float: ...
def gauquelin_sector(tjdut: float, body: int | str, method: int, geopos: Sequence[float], atpress: float = 0.0, attemp: float = 0.0, flags: int = 32770) -> float: ...

from libephemeris._core import (
    FLG_RADIANS,
    FLG_SIDEREAL,
    FLG_SPEED,
    FLG_EQUATORIAL,
    SUN,
    FLG_JPLEPH,
    FLG_SWIEPH,
    FLG_TOPOCTR,
    SIDBIT_ECL_T0,
    SIDBIT_SSY_PLANE,
    calc_ut,
    _nf_get_cached_nutation as get_cached_nutation,
    CalculationError,
    EphemerisRangeError,
    Error,
    PolarCircleError,
    validate_coordinates,
    cotrans,
    degnorm,
    difdeg2n,
)

def angular_diff(a: float, b: float) -> float: ...

from libephemeris._core import (
    _module_houses_EXTREME_LATITUDE_THRESHOLD as EXTREME_LATITUDE_THRESHOLD,
    _module_houses_CUSP_BOUNDARY_OFFSET as CUSP_BOUNDARY_OFFSET,
)
