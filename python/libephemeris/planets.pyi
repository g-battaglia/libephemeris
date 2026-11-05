"""Typed native planetary and sidereal interface."""

PosTuple = tuple[float, float, float, float, float, float]

class HeliocentricNodApsWarning(UserWarning): ...
class NutationFallbackWarning(UserWarning): ...
def calc_ut(tjdut: float, planet: int, flags: int = 258) -> tuple[PosTuple, int]: ...
def calc(tjdet: float, planet: int, flags: int = 258) -> tuple[PosTuple, int]: ...
def calc_pctr(tjdet: float, planet: int, center: int, flags: int = 258) -> tuple[PosTuple, int]:
    """TT planet-centered position; EphemerisContext.calc_pctr instead takes UT1."""
    ...
def nod_aps(tjdet: float, planet: int, method: int = 1, flags: int = 258) -> tuple[PosTuple, PosTuple, PosTuple, PosTuple]: ...
def nod_aps_ut(tjdut: float, planet: int, method: int = 1, flags: int = 258) -> tuple[PosTuple, PosTuple, PosTuple, PosTuple]: ...
def get_orbital_elements(tjdet: float, planet: int, flags: int) -> tuple[float, ...]:
    """Fifty-slot TT element record; see native help for indexed mixed units."""
    ...
def get_orbital_elements_ut(tjd_ut: float, ipl: int, iflag: int) -> tuple[float, ...]: ...
def orbit_max_min_true_distance(tjdet: float, planet: int, flags: int) -> tuple[float, float, float]: ...
def pheno(tjdet: float, planet: int, flags: int = 2) -> tuple[float, ...]:
    """Twenty-slot TT phenomena record; see native help for angles and magnitude."""
    ...
def pheno_ut(tjdut: float, planet: int, flags: int = 2) -> tuple[float, ...]: ...
def get_elongation_from_sun(tjd_ut: float, ipl: int, iflag: int = 0) -> tuple[float, bool]: ...
def get_signed_elongation(tjd_ut: float, ipl: int, iflag: int = 0) -> float: ...
def is_morning_star(tjd_ut: float, ipl: int, iflag: int = 0) -> bool: ...
def is_evening_star(tjd_ut: float, ipl: int, iflag: int = 0) -> bool: ...
def get_elongation_type(tjd_ut: float, ipl: int, iflag: int = 0) -> str: ...
def get_nutation_model() -> dict[str, str | int | bool]: ...
def set_sid_mode(mode: int, t0: float = 0.0, ayan_t0: float = 0.0) -> None: ...
def get_ayanamsa_ut(tjdut: float) -> float: ...
def get_ayanamsa(tjdet: float) -> float: ...
def get_ayanamsa_ex(tjdet: float, flags: int = 0) -> tuple[int, float]: ...
def get_ayanamsa_ex_ut(tjdut: float, flags: int = 0) -> tuple[int, float]: ...
def get_ayanamsa_name(sidmode: int) -> str: ...
def get_planet_name(planet: int) -> str: ...

from libephemeris._core import (
    _nf_method_b_accumulated_precession as method_b_accumulated_precession,
    _nf_vondrak_mean_obliquity_deg as vondrak_mean_obliquity_deg,
    _nf_vondrak_mean_obliquity_rad as vondrak_mean_obliquity_rad,
    _nf_vondrak_pn_matrix as vondrak_pn_matrix,
    _nf_vondrak_precession_matrix as vondrak_precession_matrix,
    EphemerisRangeError,
    LEBCorruptionError,
    UnknownBodyError,
    SUN,
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    URANUS,
    NEPTUNE,
    PLUTO,
    EARTH,
    MEAN_NODE,
    TRUE_NODE,
    MEAN_APOG,
    OSCU_APOG,
    INTP_APOG,
    INTP_PERG,
    CUPIDO,
    HADES,
    ZEUS,
    KRONOS,
    APOLLON,
    ADMETOS,
    VULKANUS,
    POSEIDON,
    ISIS,
    AST_OFFSET,
    NIBIRU,
    HARRINGTON,
    NEPTUNE_LEVERRIER,
    NEPTUNE_ADAMS,
    PLUTO_LOWELL,
    PLUTO_PICKERING,
    CHIRON,
    PHOLUS,
    CERES,
    PALLAS,
    JUNO,
    VESTA,
    FIXSTAR_OFFSET,
    WHITE_MOON,
    VULCAN,
    PROSERPINA,
    WALDEMATH,
    FLG_SPEED,
    FLG_SWIEPH,
    FLG_HELCTR,
    FLG_TOPOCTR,
    FLG_CENTER_BODY,
    FLG_ORBEL_AA,
    FLG_SIDEREAL,
    ANGLE_OFFSET,
    ARABIC_OFFSET,
    FLG_BARYCTR,
    FLG_TRUEPOS,
    FLG_NOABERR,
    FLG_EQUATORIAL,
    FLG_J2000,
    FLG_NONUT,
    FLG_XYZ,
    FLG_RADIANS,
    FLG_ICRS,
    FLG_SPEED3,
    FLG_NOGDEFL,
    PARS_FORTUNAE,
    PARS_SPIRITUS,
    PARS_AMORIS,
    PARS_FIDEI,
    SIDBIT_ECL_DATE,
    SIDBIT_ECL_T0,
    SIDBIT_NO_PREC_OFFSET,
    SIDBIT_PREC_ORIG,
    SIDBIT_SSY_PLANE,
    SIDBIT_USER_UT,
    SIDM_FAGAN_BRADLEY,
    SIDM_LAHIRI,
    SIDM_DELUCE,
    SIDM_RAMAN,
    SIDM_USHASHASHI,
    SIDM_KRISHNAMURTI,
    SIDM_DJWHAL_KHUL,
    SIDM_YUKTESHWAR,
    SIDM_JN_BHASIN,
    SIDM_BABYL_KUGLER1,
    SIDM_BABYL_KUGLER2,
    SIDM_BABYL_KUGLER3,
    SIDM_BABYL_HUBER,
    SIDM_BABYL_ETPSC,
    SIDM_BABYL_BRITTON,
    SIDM_ALDEBARAN_15TAU,
    SIDM_TRUE_CITRA,
    SIDM_TRUE_REVATI,
    SIDM_TRUE_PUSHYA,
    SIDM_TRUE_MULA,
    SIDM_TRUE_SHEORAN,
    SIDM_HIPPARCHOS,
    SIDM_SASSANIAN,
    SIDM_J2000,
    SIDM_J1900,
    SIDM_B1950,
    SIDM_SURYASIDDHANTA,
    SIDM_SURYASIDDHANTA_MSUN,
    SIDM_ARYABHATA,
    SIDM_ARYABHATA_MSUN,
    SIDM_ARYABHATA_522,
    SIDM_SS_REVATI,
    SIDM_SS_CITRA,
    SIDM_GALCENT_0SAG,
    SIDM_GALCENT_RGILBRAND,
    SIDM_GALCENT_MULA_WILHELM,
    SIDM_GALCENT_COCHRANE,
    SIDM_GALEQU_IAU1958,
    SIDM_GALEQU_TRUE,
    SIDM_GALEQU_MULA,
    SIDM_GALEQU_FIORENZA,
    SIDM_GALALIGN_MARDYKS,
    SIDM_VALENS_MOON,
    SIDM_LAHIRI_1940,
    SIDM_LAHIRI_VP285,
    SIDM_KRISHNAMURTI_VP291,
    SIDM_LAHIRI_ICRC,
    SIDM_USER,
    NODBIT_MEAN,
    NODBIT_OSCU_BAR,
    NODBIT_FOPOINT,
)

from libephemeris._core import (
    _module_ayanamsha_definitions_ALDEBARAN_TARGET_LON as ALDEBARAN_TARGET_LON,
    _module_ayanamsha_definitions_AYANAMSHA_DEFINING as AYANAMSHA_DEFINING,
    _module_ayanamsha_definitions_DYNAMIC_AYANAMSHA_MODES as DYNAMIC_AYANAMSHA_MODES,
    _module_ayanamsha_definitions_GALCENT_TARGET_LON as GALCENT_TARGET_LON,
    _module_ayanamsha_definitions_MARDYKS_DEFINING as MARDYKS_DEFINING,
    _module_ayanamsha_definitions_MULA_WILHELM_TARGET_LON as MULA_WILHELM_TARGET_LON,
    _module_ayanamsha_definitions_SHEORAN_TARGET_LON as SHEORAN_TARGET_LON,
    _module_ayanamsha_definitions_VALENS_MOON_AYAN_T0_DEG as VALENS_MOON_AYAN_T0_DEG,
    _module_ayanamsha_definitions_VALENS_MOON_T0_UT as VALENS_MOON_T0_UT,
)

from libephemeris._core import (
    get_sid_mode,
    get_timescale,
    get_topo,
)

class StarData:
    ra_j2000: float
    dec_j2000: float
    pm_ra: float
    pm_dec: float
    parallax: float
    radial_velocity: float
    def __init__(self, ra_j2000: float, dec_j2000: float, pm_ra: float, pm_dec: float, parallax: float=0.0, radial_velocity: float=0.0) -> None: ...
STARS: dict[str,StarData]
def consume_non_echoed_flags(flags: int, planet: int) -> int: ...
