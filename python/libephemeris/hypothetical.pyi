from typing import Dict, List, Optional, Tuple, Union
from pathlib import Path

class TPolynomial:
    constant: float
    linear: float
    def __init__(self, constant: float = 0.0, linear: float = 0.0) -> None: ...
    def evaluate(self, T: float) -> float: ...

class OrbitalElements:
    name: str
    epoch_jd: float
    equinox_jd: Optional[float]
    equinox_is_jdate: bool
    mean_anomaly: TPolynomial
    semi_axis: float
    eccentricity: TPolynomial
    arg_perihelion: TPolynomial
    asc_node: TPolynomial
    inclination: TPolynomial
    is_geocentric: bool
    line_number: int
    def __init__(self, name: str, epoch_jd: float, equinox_jd: Optional[float], equinox_is_jdate: bool, mean_anomaly: TPolynomial, semi_axis: float, eccentricity: TPolynomial, arg_perihelion: TPolynomial, asc_node: TPolynomial, inclination: TPolynomial, is_geocentric: bool = False, line_number: int = 0) -> None: ...
    def get_mean_motion(self) -> float: ...

class TransplutoKeplerianElements:
    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    M0: float
    n: float
    def __init__(self, name: str, epoch: float, a: float, e: float, i: float, omega: float, Omega: float, M0: float, n: float) -> None: ...

class UranianKeplerianElements:
    name: str
    epoch: float
    equinox_jd: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    M0: float
    n: float
    def __init__(self, name: str, epoch: float, equinox_jd: float, a: float, e: float, i: float, omega: float, Omega: float, M0: float, n: float) -> None: ...

class VulcanElements:
    name: str
    epoch: float
    a: float
    e: float
    i: float
    M0: float
    n_century: float
    omega0: float
    omega_rate: float
    Omega0: float
    Omega_rate: float
    def __init__(self, name: str, epoch: float, a: float, e: float, i: float, M0: float, n_century: float, omega0: float, omega_rate: float, Omega0: float, Omega_rate: float) -> None: ...

class WaldemathElements:
    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float
    M0: float
    M_rate: float
    omega_rate: float
    Omega_rate: float
    equinox: float
    def __init__(self, name: str, epoch: float, a: float, e: float, i: float, omega: float, Omega: float, L0: float, n: float, M0: float = 0.0, M_rate: float = 0.0, omega_rate: float = 0.0, Omega_rate: float = 0.0, equinox: float = 0.0) -> None: ...

class PickeringPlanetXElements:
    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    M0: float
    n: float
    def __init__(self, name: str, epoch: float, a: float, e: float, i: float, omega: float, Omega: float, M0: float, n: float) -> None: ...

SeorbelElements = OrbitalElements
def calc_cupido(jd_tt: float) -> Tuple[float, float, float, float, float, float]: ...
def calc_hades(jd_tt: float) -> Tuple[float, float, float, float, float, float]: ...
def calc_zeus(jd_tt: float) -> Tuple[float, float, float, float, float, float]: ...
def calc_kronos(jd_tt: float) -> Tuple[float, float, float, float, float, float]: ...
def calc_apollon(jd_tt: float) -> Tuple[float, float, float, float, float, float]: ...
def calc_admetos(jd_tt: float) -> Tuple[float, float, float, float, float, float]: ...
def calc_vulkanus(jd_tt: float) -> Tuple[float, float, float, float, float, float]: ...
def calc_poseidon(jd_tt: float) -> Tuple[float, float, float, float, float, float]: ...
def calc_transpluto(jd_tt: float) -> Tuple[float, float, float, float, float, float]: ...
def calc_vulcan(jd_tt: float) -> Tuple[float, float, float, float, float, float]: ...
def calc_waldemath(jd_tt: float) -> Tuple[float, float, float, float, float, float]: ...
def calc_proserpina(jd_tt: float) -> Tuple[float, float, float, float, float, float]: ...
def calc_planet_x_pickering(jd_tt: float) -> Tuple[float, float, float, float, float, float]: ...
def calc_white_moon_position(jd_tt: float, use_true_lilith: bool=False) -> Tuple[float, float, float, float, float, float]: ...
def calc_uranian_planet(body_id: int, jd_tt: float) -> Tuple[float, float, float, float, float, float]: ...
def calc_hypothetical_position(ipl: int, jd_tt: float) -> Tuple[float, float, float, float, float, float]: ...
def parse_orbital_elements(filepath: Union[str, Path]) -> List[OrbitalElements]: ...
def get_orbital_body_by_name(elements: List[OrbitalElements], name: str) -> Optional[OrbitalElements]: ...
def calc_orbital_position(elem: OrbitalElements, jd_tt: float) -> Tuple[float, float, float, float, float, float]: ...
def get_bundled_fictitious_orbits_path() -> Path: ...
def load_bundled_fictitious_orbits() -> List[OrbitalElements]: ...
def parse_seorbel(filepath: Union[str, Path]) -> List[OrbitalElements]: ...
def get_seorbel_body_by_name(elements: List[OrbitalElements], name: str) -> Optional[OrbitalElements]: ...
def calc_seorbel_position(elem: OrbitalElements, jd_tt: float) -> Tuple[float, float, float, float, float, float]: ...
def get_bundled_seorbel_path() -> Path: ...
def load_bundled_seorbel() -> List[OrbitalElements]: ...

URANIAN_KEPLERIAN_ELEMENTS: dict[int, UranianKeplerianElements]
TRANSPLUTO_KEPLERIAN_ELEMENTS: TransplutoKeplerianElements
VULCAN_ELEMENTS: VulcanElements
WALDEMATH_ELEMENTS: WaldemathElements
PICKERING_PLANET_X_ELEMENTS: PickeringPlanetXElements

from typing import Optional, Tuple, List, Dict, Sequence

class UranianElements:
    name: str
    L0: float
    n: float
    amplitude: float = 0.0
    phase: float = 0.0
    phase_rate: float = 0.0
    def __init__(self, name:  str,L0:  float,n:  float,amplitude:  float = 0.0,phase:  float = 0.0,phase_rate:  float = 0.0) -> None: ...

class HypotheticalElements:
    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    M0: float
    n: float
    def __init__(self, name:  str,epoch:  float,a:  float,e:  float,i:  float,omega:  float,Omega:  float,M0:  float,n:  float) -> None: ...

class CupidoKeplerianElements:
    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float
    def __init__(self, name:  str,epoch:  float,a:  float,e:  float,i:  float,omega:  float,Omega:  float,L0:  float,n:  float) -> None: ...

class HadesKeplerianElements:
    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    M0: float
    n: float
    def __init__(self, name:  str,epoch:  float,a:  float,e:  float,i:  float,omega:  float,Omega:  float,M0:  float,n:  float) -> None: ...

class ZeusKeplerianElements:
    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float
    def __init__(self, name:  str,epoch:  float,a:  float,e:  float,i:  float,omega:  float,Omega:  float,L0:  float,n:  float) -> None: ...

class KronosKeplerianElements:
    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float
    def __init__(self, name:  str,epoch:  float,a:  float,e:  float,i:  float,omega:  float,Omega:  float,L0:  float,n:  float) -> None: ...

class ApollonKeplerianElements:
    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float
    def __init__(self, name:  str,epoch:  float,a:  float,e:  float,i:  float,omega:  float,Omega:  float,L0:  float,n:  float) -> None: ...

class AdmetosKeplerianElements:
    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float
    def __init__(self, name:  str,epoch:  float,a:  float,e:  float,i:  float,omega:  float,Omega:  float,L0:  float,n:  float) -> None: ...

class VulkanusKeplerianElements:
    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float
    def __init__(self, name:  str,epoch:  float,a:  float,e:  float,i:  float,omega:  float,Omega:  float,L0:  float,n:  float) -> None: ...

class PoseidonKeplerianElements:
    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float
    def __init__(self, name:  str,epoch:  float,a:  float,e:  float,i:  float,omega:  float,Omega:  float,L0:  float,n:  float) -> None: ...

class LowellPlanetXElements:
    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    M0: float
    n: float
    def __init__(self, name:  str,epoch:  float,a:  float,e:  float,i:  float,omega:  float,Omega:  float,M0:  float,n:  float) -> None: ...
def is_hypothetical_body(ipl: int) -> bool:
    ...
def get_hypothetical_name(ipl: int) -> str:
    ...
def calc_fictitious_position(ipl: int, jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    ...
def calc_uranian_longitude(ipl: int, jd_tt: float) -> float:
    ...
def calc_uranian_position(ipl: int, jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    ...
def calc_transpluto_position(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    ...
def calc_waldemath_position(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    ...
def calc_planet_x_lowell(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    ...
def list_hypothetical_bodies() -> Dict[int, str]:
    ...
URANIAN_ELEMENTS: dict
HYPOTHETICAL_ELEMENTS: dict
HYPOTHETICAL_PROVENANCE: dict
HYPOTHETICAL_NAMES: dict
FICTITIOUS_ORBITAL_ELEMENTS: dict
CUPIDO_KEPLERIAN_ELEMENTS: dict
HADES_KEPLERIAN_ELEMENTS: dict
ZEUS_KEPLERIAN_ELEMENTS: dict
KRONOS_KEPLERIAN_ELEMENTS: dict
APOLLON_KEPLERIAN_ELEMENTS: dict
ADMETOS_KEPLERIAN_ELEMENTS: dict
VULKANUS_KEPLERIAN_ELEMENTS: dict
POSEIDON_KEPLERIAN_ELEMENTS: dict
LOWELL_PLANET_X_ELEMENTS: dict

from .constants import (
    ADMETOS, APOLLON, CUPIDO, FICT_OFFSET, HADES, HARRINGTON, ISIS, KRONOS, NEPTUNE_ADAMS, NEPTUNE_LEVERRIER, NIBIRU, PLANET_X_ADAMS, PLANET_X_LEVERRIER, PLANET_X_LOWELL, PLANET_X_PICKERING, PLUTO_LOWELL, PLUTO_PICKERING, POSEIDON, PROSERPINA, TRANSPLUTO, VULCAN, VULKANUS, WALDEMATH, WHITE_MOON, ZEUS,
)
