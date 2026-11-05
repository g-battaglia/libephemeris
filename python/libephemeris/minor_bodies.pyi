class OrbitalElements:
    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    M0: float
    n: float
    def __init__(self, name: str, epoch: float, a: float, e: float, i: float,
                 omega: float, Omega: float, M0: float, n: float) -> None: ...
    def validate(self) -> list[str]: ...

class ResonanceInfo:
    p: int
    q: int
    name: str
    a_resonant: float
    tolerance: float
    def __init__(self, p: int, q: int, name: str, a_resonant: float, tolerance: float) -> None: ...

class ResonanceResult:
    is_resonant: bool
    resonance: ResonanceInfo | None
    deviation: float
    warning_message: str | None
    def __init__(self, is_resonant: bool, resonance: ResonanceInfo | None, deviation: float, warning_message: str | None) -> None: ...

def get_asteroid_number(name: str, timeout: float = 30.0) -> int | None: ...
def clear_asteroid_name_cache() -> int: ...
def calc_asteroid_by_number(asteroid_number: int, jd_tt: float, use_spk: bool = True) -> tuple[float, float, float]: ...
def fetch_orbital_elements_from_sbdb(asteroid_number: int, timeout: float = 30.0) -> OrbitalElements | None: ...
def clear_asteroid_elements_cache() -> int: ...
def auto_download_asteroid_spk(body_id: int, jd_start: float | None = None, jd_end: float | None = None, force: bool = False, boundary_padding_days: float = 0.0) -> str | None: ...
def is_spk_available_for_body(body_id: int) -> bool: ...
def ensure_major_asteroid_spk(body_id: int, jd: float | None = None) -> bool: ...
def is_major_asteroid(body_id: int) -> bool: ...
def get_major_asteroid_info(body_id: int) -> tuple[int, str, int, str] | None: ...
def list_major_asteroids() -> list[tuple[int, str]]: ...
def detect_mean_motion_resonance(elements: OrbitalElements, tolerance: float = 0.02) -> ResonanceResult: ...

from typing import Optional, Tuple, List, Dict, Sequence

class LibrationParameters:
    amplitude: float
    period: float
    center: float
    phase_j2000: float
    resonance_p: int
    resonance_q: int
    def __init__(self, amplitude:  float,period:  float,center:  float,phase_j2000:  float,resonance_p:  int,resonance_q:  int) -> None: ...
def validate_elements(elements: OrbitalElements) -> None:
    ...
def is_body_resonant(body_id: int, tolerance: float=DEFAULT_RESONANCE_TOLERANCE) -> bool:
    ...
def get_resonance_info(body_id: int) -> Optional[ResonanceResult]:
    ...
def calc_neptune_mean_longitude(jd_tt: float) -> float:
    ...
def calc_resonant_argument_plutino(elements: OrbitalElements, jd_tt: float, omega_pert: float, M_pert: float) -> float:
    ...
def calc_libration_correction(body_id: int, jd_tt: float) -> float:
    ...
def has_libration_model(body_id: int) -> bool:
    ...
def get_libration_parameters(body_id: int) -> Optional[LibrationParameters]:
    ...
def calc_secular_perturbation_rates(elements: OrbitalElements) -> Tuple[float, float, float]:
    ...
def apply_secular_perturbations(elements: OrbitalElements, jd_tt: float, include_perturbations: bool=True) -> Tuple[float, float, float, float, float, float]:
    ...
def solve_kepler_equation_elliptic(M: float, e: float, tol: float=1e-08) -> float:
    ...
def solve_kepler_equation_hyperbolic(M: float, e: float, tol: float=1e-08) -> float:
    ...
def solve_barker_equation(M: float) -> float:
    ...
def solve_kepler_equation(M: float, e: float, tol: float=1e-08) -> float:
    ...
def calc_minor_body_position(elements: OrbitalElements, jd_tt: float, include_perturbations: bool=True, body_id: Optional[int]=None) -> Tuple[float, float, float]:
    ...
def calc_minor_body_heliocentric(body_id: int, jd_tt: float, use_spk: bool=True) -> Tuple[float, float, float]:
    ...
def is_spk_downloadable(body_id: int) -> bool:
    ...
def is_in_jpl_major_body_index(body_id: int) -> bool:
    ...
MINOR_BODY_ELEMENTS: dict
PLUTINO_LIBRATION_PARAMS: dict
YARKOVSKY_DA_DT: dict

from .constants import (
    AMOR, APOPHIS, ASBOLUS, BENNU, CERES, CHARIKLO, CHIRON, DAVIDA, ERIS, EROS, EUROPA_AST, GONGGONG, HAUMEA, HIDALGO, HYGIEA, ICARUS, INTERAMNIA, ITOKAWA, IXION, JUNO, LILITH_AST, MAKEMAKE, NAIF_CERES, NAIF_CHIRON, NAIF_JUNO, NAIF_PALLAS, NAIF_VESTA, NESSUS, ORCUS, PALLAS, PANDORA_AST, PHOLUS, PSYCHE, QUAOAR, RYUGU, SAPPHO, SEDNA, SYLVIA, TORO, TOUTATIS, VARUNA, VESTA,
)

from ._core import (
    _Minor_GM_SUN as GM_SUN,
    _Minor_MASS_RATIO_JUPITER as MASS_RATIO_JUPITER,
    _Minor_MASS_RATIO_SATURN as MASS_RATIO_SATURN,
    _Minor_MASS_RATIO_URANUS as MASS_RATIO_URANUS,
    _Minor_MASS_RATIO_NEPTUNE as MASS_RATIO_NEPTUNE,
    _Minor_JUPITER_A as JUPITER_A,
    _Minor_JUPITER_E as JUPITER_E,
    _Minor_JUPITER_I as JUPITER_I,
    _Minor_JUPITER_OMEGA as JUPITER_OMEGA,
    _Minor_JUPITER_NODE as JUPITER_NODE,
    _Minor_JUPITER_N as JUPITER_N,
    _Minor_SATURN_A as SATURN_A,
    _Minor_SATURN_E as SATURN_E,
    _Minor_SATURN_I as SATURN_I,
    _Minor_SATURN_OMEGA as SATURN_OMEGA,
    _Minor_SATURN_NODE as SATURN_NODE,
    _Minor_SATURN_N as SATURN_N,
    _Minor_URANUS_A as URANUS_A,
    _Minor_URANUS_E as URANUS_E,
    _Minor_URANUS_I as URANUS_I,
    _Minor_URANUS_OMEGA as URANUS_OMEGA,
    _Minor_URANUS_NODE as URANUS_NODE,
    _Minor_URANUS_N as URANUS_N,
    _Minor_NEPTUNE_A as NEPTUNE_A,
    _Minor_NEPTUNE_E as NEPTUNE_E,
    _Minor_NEPTUNE_I as NEPTUNE_I,
    _Minor_NEPTUNE_OMEGA as NEPTUNE_OMEGA,
    _Minor_NEPTUNE_NODE as NEPTUNE_NODE,
    _Minor_NEPTUNE_N as NEPTUNE_N,
    _Minor_JUPITER_DE as JUPITER_DE,
    _Minor_JUPITER_DI as JUPITER_DI,
    _Minor_JUPITER_DOMEGA as JUPITER_DOMEGA,
    _Minor_JUPITER_DNODE as JUPITER_DNODE,
    _Minor_SATURN_DE as SATURN_DE,
    _Minor_SATURN_DI as SATURN_DI,
    _Minor_SATURN_DOMEGA as SATURN_DOMEGA,
    _Minor_SATURN_DNODE as SATURN_DNODE,
    _Minor_URANUS_DE as URANUS_DE,
    _Minor_URANUS_DI as URANUS_DI,
    _Minor_URANUS_DOMEGA as URANUS_DOMEGA,
    _Minor_URANUS_DNODE as URANUS_DNODE,
    _Minor_NEPTUNE_DE as NEPTUNE_DE,
    _Minor_NEPTUNE_DI as NEPTUNE_DI,
    _Minor_NEPTUNE_DOMEGA as NEPTUNE_DOMEGA,
    _Minor_NEPTUNE_DNODE as NEPTUNE_DNODE,
    _Minor_J2000_EPOCH as J2000_EPOCH,
    _Minor_NEPTUNE_MEAN_MOTION as NEPTUNE_MEAN_MOTION,
    _Minor_EPOCH_MERGE_TOL_DAYS as EPOCH_MERGE_TOL_DAYS,
    _Minor_SBDB_API_URL as SBDB_API_URL,
    _Minor_HORIZONS_SPK_JD_MIN as HORIZONS_SPK_JD_MIN,
    _Minor_HORIZONS_SPK_JD_MAX as HORIZONS_SPK_JD_MAX,
)

from libephemeris._core import (
    get_logger,
    EphemerisRangeError,
)

from libephemeris._core import (
    _module_minor_bodies_DEFAULT_RESONANCE_TOLERANCE as DEFAULT_RESONANCE_TOLERANCE,
    _module_minor_bodies_NEPTUNE_MEAN_LONGITUDE_J2000 as NEPTUNE_MEAN_LONGITUDE_J2000,
    _module_minor_bodies_SPK_DOWNLOADABLE_ASTEROIDS as SPK_DOWNLOADABLE_ASTEROIDS,
    _module_minor_bodies_MAJOR_ASTEROID_SPK_INFO as MAJOR_ASTEROID_SPK_INFO,
    _module_minor_bodies_HORIZONS_SPK_DATE_MIN as HORIZONS_SPK_DATE_MIN,
    _module_minor_bodies_HORIZONS_SPK_DATE_MAX as HORIZONS_SPK_DATE_MAX,
)

NEPTUNE_RESONANCES: list[ResonanceInfo]
MINOR_BODY_ELEMENTS_MULTI: dict[int,list[OrbitalElements]]
