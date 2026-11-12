from collections.abc import Sequence

class StarData:
    ra_j2000: float
    dec_j2000: float
    pm_ra: float
    pm_dec: float
    parallax_mas: float
    radial_km_per_s: float
    def __init__(self, ra_j2000: float, dec_j2000: float, pm_ra: float, pm_dec: float, parallax_mas: float = 0.0, radial_km_per_s: float = 0.0) -> None: ...

class PhotometricBand:
    Unspecified: PhotometricBand
    JohnsonV: PhotometricBand
    HipparcosHp: PhotometricBand

class PhotometricCatalogueField:
    Unqualified: PhotometricCatalogueField
    HipparcosMainVmag: PhotometricCatalogueField

class PhotometricStatistic:
    Unspecified: PhotometricStatistic
    Median: PhotometricStatistic

class PhotometricVariability:
    Unknown: PhotometricVariability
    Constant: PhotometricVariability
    Periodic: PhotometricVariability
    Unsolved: PhotometricVariability

class StarPhotometry:
    value: float
    band: PhotometricBand
    catalogue_field: PhotometricCatalogueField
    statistic: PhotometricStatistic
    epoch: str | None
    variability: PhotometricVariability
    historical_unqualified_value: float | None

class StarCatalogEntry:
    id: int
    name: str
    nomenclature: str
    hip_number: int
    data: StarData
    magnitude: float
    photometry: StarPhotometry
    def __init__(self, id: int, name: str, nomenclature: str, hip_number: int, data: StarData, magnitude: float) -> None: ...

def star_photometry(star_id: int) -> StarPhotometry | None: ...
def list_fixed_stars() -> tuple[StarCatalogEntry, ...]: ...
def propagate_proper_motion(ra_epoch: float, dec_epoch: float, pm_ra_cosdec: float, pm_dec: float, from_jd: float, to_jd: float) -> tuple[float, float]: ...
def fixstar(star: str, tjdet: float, flags: int = 2) -> tuple[tuple[float, float, float, float, float, float], str, int]: ...
def fixstar_ut(star: str, tjdut: float, flags: int = 2) -> tuple[tuple[float, float, float, float, float, float], str, int]: ...
def fixstar2(star: str, tjdet: float, flags: int = 2) -> tuple[tuple[float, float, float, float, float, float], str, int]: ...
def fixstar2_ut(star: str, tjdut: float, flags: int = 2) -> tuple[tuple[float, float, float, float, float, float], str, int]: ...
def batch_fixstars_ut(stars: Sequence[str], tjdut: float, flags: int = 2, *, skip_errors: bool = False) -> tuple[tuple[tuple[float, float, float, float, float, float], str, int] | None, ...]: ...
def fixstar_mag(star: str) -> tuple[float, str]: ...
def fixstar2_mag(star: str) -> tuple[float, str]: ...

STAR_NAME_TO_HIP: dict[str, int]
STAR_ALIASES: dict[str, int]
STAR_CATALOG: tuple[StarCatalogEntry, ...]
FIXED_STARS: dict[int, StarData]
def resolve_star_name(name: str) -> int | None: ...
def get_canonical_star_name(star_id: int) -> str | None: ...
def get_hip_from_star_name(name: str) -> int | None: ...
def calc_fixed_star_position(star_id: int, jd_tt: float, noaberr: bool = False, nogdefl: bool = False, j2000_frame: bool = False, topo: tuple[float, float, float] | None = None, center: int = 0) -> tuple[float, float, float]: ...
def calc_fixed_star_velocity(star_id: int, jd_tt: float, noaberr: bool = False, nogdefl: bool = False, j2000_frame: bool = False, topo: tuple[float, float, float] | None = None, center: int = 0) -> tuple[float, float, float, float, float, float]: ...

from .constants import (
    ACHERNAR, ACRUX, ACUBENS, ADHARA, AIN, ALBIREO, ALCOR, ALCYONE, ALDEBARAN, ALGEDI, ALGENIB, ALGIEBA, ALGOL, ALGORAB, ALHENA, ALIOTH, ALKAID, ALMACH, ALNILAM, ALNITAK, ALPHARD, ALPHECCA, ALPHERATZ, ALRESCHA, ALTAIR, ANTARES, ARCTURUS, ASCELLA, ASELLUS_AUSTRALIS, ASELLUS_BOREALIS, ASTEROPE, ATLAS, BELLATRIX, BETELGEUSE, CANOPUS, CAPELLA, CASTOR, CELAENO, DABIH, DAYS_PER_JULIAN_YEAR, DELTA_CRUCIS, DENEB, DENEBOLA, DENEB_ALGEDI, DIPHDA, DSCHUBBA, DUBHE, ELECTRA, ELNATH, EPSILON_CENTAURI, ETAMIN, ETA_CENTAURI, ETA_PISCIUM, FLG_BARYCTR, FLG_EQUATORIAL, FLG_HELCTR, FLG_J2000, FLG_MOSEPH, FLG_NOABERR, FLG_NOGDEFL, FLG_NONUT, FLG_RADIANS, FLG_SIDEREAL, FLG_SPEED, FLG_SPEED3, FLG_SWIEPH, FLG_TOPOCTR, FLG_TRUEPOS, FLG_XYZ, FOMALHAUT, GACRUX, GRAFFIAS, HADAR, HAMAL, J2000, KAUS_AUSTRALIS, KAUS_BOREALIS, KAUS_MEDIA, KOCHAB, LESATH, MAIA, MARKAB, MEGREZ, MEISSA, MENKAR, MENKENT, MERAK, MEROPE, MESARTHIM, MIMOSA, MINTAKA, MIRA, MIRACH, MIZAR, MUHLIFAIN, NASHIRA, NUNKI, PHECDA, PLEIONE, POLARIS, POLLUX, PRIMA_HYADUM, PROCYON, PROPUS, RASALGETHI, RASALHAGUE, REGULUS, RIGEL, RIGIL_KENT, SADALMELIK, SADALSUUD, SAIPH, SARGAS, SCHEAT, SECUNDA_HYADUM, SHAULA, SHERATAN, SIRIUS, SKAT, SPICA_STAR, TARF, TAYGETA, TEJAT, THETA_TAURI, THUBAN, UNUKALHAI, VEGA, VINDEMIATRIX, WASAT, WEZEN, ZETA_CENTAURI, ZOSMA, ZUBENELGENUBI, ZUBENESCHAMALI,
)

from ._core import (
    _module_fixed_stars_GREEK_LETTER_ABBREV as GREEK_LETTER_ABBREV,
    _module_fixed_stars_CONSTELLATION_ABBREV as CONSTELLATION_ABBREV,
)

from libephemeris._core import (
    cotrans,
    _nf_get_true_obliquity as get_true_obliquity,
    _nf_get_mean_obliquity as get_mean_obliquity,
    Error,
    LEBCorruptionError,
    StarNotFoundError,
)
