"""Native fixed-star catalog and apparent-position functions."""

from ._core import (
    PhotometricBand, PhotometricCatalogueField, PhotometricStatistic,
    PhotometricVariability, StarPhotometry, star_photometry,
    StarCatalogEntry, batch_fixstars_ut, fixstar, fixstar2, fixstar2_mag,
    fixstar2_ut, fixstar_mag, fixstar_ut, list_fixed_stars, propagate_proper_motion,
)
from ._core import _StarData as StarData

from ._core import (
    calc_fixed_star_position, calc_fixed_star_velocity, resolve_star_name,
    get_canonical_star_name, get_hip_from_star_name, STAR_NAME_TO_HIP,
    STAR_ALIASES, STAR_CATALOG, FIXED_STARS,
)

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
