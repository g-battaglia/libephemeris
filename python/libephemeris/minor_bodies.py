"""Native minor-body elements, catalog queries and SPK provisioning."""

from ._core import _MinorOrbitalElements as OrbitalElements
from ._core import _ResonanceInfo as ResonanceInfo, _ResonanceResult as ResonanceResult
from ._core import (
    auto_download_asteroid_spk, calc_asteroid_by_number,
    clear_asteroid_elements_cache, clear_asteroid_name_cache,
    detect_mean_motion_resonance, ensure_major_asteroid_spk,
    fetch_orbital_elements_from_sbdb, get_asteroid_number,
    get_major_asteroid_info, is_major_asteroid, is_spk_available_for_body,
    list_major_asteroids,
)

from ._core import (
    LibrationParameters,
    validate_elements,
    is_body_resonant,
    get_resonance_info,
    calc_neptune_mean_longitude,
    calc_resonant_argument_plutino,
    calc_libration_correction,
    has_libration_model,
    get_libration_parameters,
    calc_secular_perturbation_rates,
    apply_secular_perturbations,
    solve_kepler_equation_elliptic,
    solve_kepler_equation_hyperbolic,
    solve_barker_equation,
    solve_kepler_equation,
    calc_minor_body_position,
    calc_minor_body_heliocentric,
    is_spk_downloadable,
    is_in_jpl_major_body_index,
    MINOR_BODY_ELEMENTS,
    PLUTINO_LIBRATION_PARAMS,
    YARKOVSKY_DA_DT,
)

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

from ._core import NEPTUNE_RESONANCES, MINOR_BODY_ELEMENTS_MULTI
