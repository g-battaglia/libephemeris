# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
# Copyright (c) 2025-2026 Giacomo Battaglia
"""libephemeris.contrib — extended astrology helpers.

Zodiac and nakshatra constants, Vedic planet IDs, classical aspect
angles, longitude conversions, dignity calculations, and a small set
of time/coordinate utilities. Names and signatures follow established
conventions in Python astrology tooling; implementations are derived
from standard zodiacal arithmetic and classical/Vedic tradition.

Contents
--------
- Constants: 12 zodiac signs (English + Sanskrit names), 27 nakshatras,
  Vedic planet IDs, classical aspect angles.
- Calculation helpers:
    long2rasi, long2nakshatra, long2navamsa, rasinorm, rasi_dif,
    rasi_dif2, antiscion, degsplit, signtostr, get_nakshatra_name,
    house_system_name, lord, ochchabala, residential_strength,
    naisargika_relation, tatkalika_relation, raman_houses,
    saturn_4_stars.
- Aspect helpers:
    match_aspect, match_aspect2, match_aspect3, match_aspect4,
    next_aspect, next_aspect2, next_aspect_cusp, next_aspect_cusp2,
    next_aspect_with, next_aspect_with2.
- Time/JD helpers:
    jdnow, jduration, t2i, dt2i, jd2isostr, years_diff.
- Coordinate helpers:
    geoc2d, geolat2c, geolon2c.
- Retrograde/cycle helpers:
    next_retro.
- Convenience aliases for core ephemeris functions:
    calc_ut, julday, revjul (re-exported from the parent package).
- Database / atlas / timezone:
    atlas_connect, atlas_close, atlas_search, atlas_countries_list,
    db_connect, db_close, tzabbr_find, tzabbr_list — these require
    optional SQLite databases and are exposed as stubs raising
    NotImplementedError. Use the upstream reference contrib module
    if these features are needed.

The User, Data, and Error names are exposed as lightweight classes to
mirror the upstream's structural API.
"""

from __future__ import annotations as _annotations  # noqa: F401

import datetime as _dt
import math as _math
from typing import Any as _Any

# Re-exports of canonical core functions (must match reference contrib)
from .. import calc_ut as _calc_ut
from .. import julday as _julday
from .. import revjul as _revjul
from .. import FLG_SPEED as _FLG_SPEED

calc_ut = _calc_ut
julday = _julday
revjul = _revjul


# =============================================================================
# CLASSES (structural compatibility)
# =============================================================================


class Error(Exception):
    """Module-level error class for contrib calculations."""


class Data:
    """Lightweight data container used by aspect/return functions."""


class User:
    """Placeholder user-data container."""


# =============================================================================
# CONSTANTS — Zodiac (English + Sanskrit Rasi)
# =============================================================================

ARIES: int = 0
TAURUS: int = 1
GEMINI: int = 2
CANCER: int = 3
LEO: int = 4
VIRGO: int = 5
LIBRA: int = 6
SCORPIO: int = 7
SAGITTARIUS: int = 8
CAPRICORN: int = 9
AQUARIUS: int = 10
PISCES: int = 11

# Sanskrit names (Rasi)
MESHA: int = ARIES
VRISHABA: int = TAURUS
MITHUNA: int = GEMINI
KATAKA: int = CANCER
SIMHA: int = LEO
KANYA: int = VIRGO
THULA: int = LIBRA
VRISHIKA: int = SCORPIO
DHANUS: int = SAGITTARIUS
MAKARA: int = CAPRICORN
KUMBHA: int = AQUARIUS
MEENA: int = PISCES


# =============================================================================
# CONSTANTS — Nakshatras (27 lunar mansions)
# =============================================================================

ASWINI: int = 0
BHARANI: int = 1
KRITIKHA: int = 2
ROHINI: int = 3
MRIGASIRA: int = 4
ARIDRA: int = 5
PUNARVASU: int = 6
PUSHYAMI: int = 7
ASLESHA: int = 8
MAKHA: int = 9
PUBBA: int = 10
UTTARA: int = 11
HASTA: int = 12
CHITTA: int = 13
SWATHI: int = 14
VISHAKA: int = 15
ANURADHA: int = 16
JYESTA: int = 17
MOOLA: int = 18
POORVASHADA: int = 19
UTTARASHADA: int = 20
SRAVANA: int = 21
DHANISHTA: int = 22
SATABISHA: int = 23
POORVABHADRA: int = 24
UTTARABHADRA: int = 25
REVATHI: int = 26


# =============================================================================
# CONSTANTS — Vedic planet IDs
# =============================================================================

# Sun
SURYA: int = 0
RAVI: int = 0
# Moon
CHANDRA: int = 1
SOMA: int = 1
# Mercury
BUDHA: int = 2
SOUMYA: int = 2
THAMA: int = 10  # alternate Mercury name (per reference contrib)
# Venus
SUKRA: int = 3
BHARGAVA: int = 3
# Mars
KUJA: int = 4
ANGARAKA: int = 4
# Jupiter
GURU: int = 5
BRIHASPATI: int = 5
# Saturn
SANI: int = 6
MANDA: int = 6
# Nodes
RAHU: int = 10
KETU: int = -10
SIKHI: int = -10  # alternate Ketu name


# =============================================================================
# CONSTANTS — Aspects (angular separations in degrees)
# =============================================================================

CONJUNCTION: int = 0
SEMISEXTILE: int = 30
SEXTILE: int = 60
SQUARE: int = 90
TRINE: int = 120
QUINCUNX: int = 150
OPPOSITION: int = 180

SEMISQUARE: int = 45
SESQUISQUARE: int = 135
QUINTILE: int = 72
BIQUINTILE: int = 144
SEMIQUINTILE: int = 36

# Septile family (360 / 7)
SEPTILE: float = 360.0 / 7.0
BISEPTILE: float = 720.0 / 7.0
TRISEPTILE: float = 1080.0 / 7.0

# Novile family (360 / 9)
NOVILE: int = 40
SEMINOVILE: int = 20
BINOVILE: int = 80
QUATRONOVILE: int = 160

# Undecile family (360 / 11)
UNDECILE: float = 360.0 / 11.0
BIUNDECILE: float = 720.0 / 11.0
TRIUNDECILE: float = 1080.0 / 11.0
QUADUNDECILE: float = 1440.0 / 11.0
QUINUNDECILE: float = 1800.0 / 11.0

# 15° family
SQUISEXTILE: int = 15
SQUISQUARE: float = 22.5


# =============================================================================
# Internal name tables (used by signtostr / get_nakshatra_name etc.)
# =============================================================================

_RASI_NAMES = (
    "Aries", "Taurus", "Gemini", "Cancer", "Leo", "Virgo",
    "Libra", "Scorpio", "Sagittarius", "Capricorn", "Aquarius", "Pisces",
)

_NAKSHATRA_NAMES = (
    "Aswini", "Bharani", "Kritikha", "Rohini", "Mrigasira", "Aridra",
    "Punarvasu", "Pushyami", "Aslesha", "Makha", "Pubba", "Uttara",
    "Hasta", "Chitta", "Swathi", "Vishaka", "Anuradha", "Jyesta",
    "Moola", "Poorvashada", "Uttarashada", "Sravana", "Dhanishta",
    "Satabisha", "Poorvabhadra", "Uttarabhadra", "Revathi",
)

_HOUSE_SYSTEM_NAMES = {
    "A": "Equal",
    "E": "Equal",
    "B": "Alcabitius",
    "C": "Campanus",
    "D": "Equal MC",
    "F": "Carter poli-equatorial",
    "G": "Gauquelin sectors",
    "H": "Horizontal/Azimuthal",
    "I": "Sunshine",
    "i": "Sunshine alternative",
    "K": "Koch",
    "L": "Pullen S-Delta",
    "M": "Morinus",
    "N": "Whole sign",
    "O": "Porphyry",
    "P": "Placidus",
    "Q": "Pullen S-Ratio",
    "R": "Regiomontanus",
    "S": "Sripati",
    "T": "Polich/Page (topocentric)",
    "U": "Krusinski-Pisa-Goelzer",
    "V": "Vehlow equal",
    "W": "Whole sign",
    "X": "Axial rotation/Meridian",
    "Y": "APC",
}

# Nakshatra width in degrees (360 / 27)
_NAKSHATRA_WIDTH: float = 360.0 / 27.0
# Pada width (a quarter of a nakshatra)
_PADA_WIDTH: float = _NAKSHATRA_WIDTH / 4.0
# Navamsa width (a ninth of a sign)
_NAVAMSA_WIDTH: float = 30.0 / 9.0


# =============================================================================
# Core Vedic conversion helpers
# =============================================================================


def _normalize_longitude(deg: float) -> float:
    """Normalize an ecliptic longitude to [0, 360)."""
    return deg % 360.0


def long2rasi(deg: float) -> int:
    """Return the Rasi (sign 0=Aries through 11=Pisces) for a longitude."""
    return int(_normalize_longitude(deg) // 30.0)


def long2nakshatra(deg: float) -> tuple[int, int]:
    """Return (nakshatra_id 0..26, pada 0..3) for a longitude."""
    norm = _normalize_longitude(deg)
    nak = int(norm // _NAKSHATRA_WIDTH)
    within = norm - nak * _NAKSHATRA_WIDTH
    pada = int(within // _PADA_WIDTH)
    return (nak, pada)


def long2navamsa(deg: float) -> int:
    """Return the navamsa (0..11) for a longitude. The navamsa system
    divides each sign into 9 parts; the navamsa of a body depends on
    the element of its sign:
      Fire signs (Aries, Leo, Sagittarius): navamsa starts from Aries
      Earth signs (Taurus, Virgo, Capricorn): from Capricorn
      Air signs (Gemini, Libra, Aquarius): from Libra
      Water signs (Cancer, Scorpio, Pisces): from Cancer
    So: navamsa = ((rasi % 4) * 9 + floor((deg_within_sign) / (30/9))) % 12.
    """
    norm = _normalize_longitude(deg)
    rasi = int(norm // 30.0)
    within = norm - rasi * 30.0
    sub = int(within // _NAVAMSA_WIDTH)
    return ((rasi % 4) * 9 + sub) % 12


def rasinorm(rasi: int) -> int:
    """Normalize a rasi number to [0, 11]."""
    return int(rasi) % 12


def rasi_dif(r1: int, r2: int) -> int:
    """Forward rasi distance from r1 to r2 (in [0, 11])."""
    return (rasinorm(r2) - rasinorm(r1)) % 12


def rasi_dif2(r1: int, r2: int) -> int:
    """Signed rasi difference from r1 to r2 in (-6, 6].

    Useful for picking the shorter direction around the zodiac wheel.
    """
    d = rasi_dif(r1, r2)
    if d > 6:
        d -= 12
    return d


def antiscion(lon: float) -> tuple[
    tuple[float, float, float, float, float, float],
    tuple[float, float, float, float, float, float],
]:
    """Return (antiscion, contra-antiscion) pairs for a longitude.

    Antiscion is the longitude reflected across the Cancer/Capricorn
    axis (0° Cancer = 90°). Mathematically: antiscion = 180° - lon.
    The contra-antiscion is the antipode of the antiscion.

    Both are returned as 6-tuples (lon, lat, dist, speed_lon,
    speed_lat, speed_dist) with zero latitude/distance/speed to mirror
    the upstream return shape.
    """
    norm = _normalize_longitude(lon)
    ant = _normalize_longitude(180.0 - norm)
    contra = _normalize_longitude(ant + 180.0)
    # Mirror the upstream's sign behaviour for zero entries
    zero = -0.0 if (norm != 0.0 and norm != 180.0) else 0.0
    ant_tuple = (ant, 0.0, 0.0, zero, 0.0, 0.0)
    contra_tuple = (contra, zero, 0.0, zero, zero, 0.0)
    return ant_tuple, contra_tuple


def degsplit(deg: float) -> tuple[int, int, int, int]:
    """Split a longitude into (degrees within sign, sign id, arcminutes,
    arcseconds).

    Pure positional decomposition; no normalization performed beyond
    bringing the input into [0, 360). Rounded seconds at the 60 mark
    are carried into the minute / degree / sign fields so the returned
    tuple is always a valid sexagesimal value."""
    norm = _normalize_longitude(deg)
    sign = int(norm // 30.0)
    within = norm - sign * 30.0
    d = int(within)
    rem_min = (within - d) * 60.0
    m = int(rem_min)
    rem_sec = (rem_min - m) * 60.0
    s = int(round(rem_sec))
    # Carry rounded overflow through s → m → d → sign so we never
    # emit s == 60 or m == 60 (invalid sexagesimal forms).
    if s >= 60:
        s -= 60
        m += 1
    if m >= 60:
        m -= 60
        d += 1
    if d >= 30:
        d -= 30
        sign = (sign + 1) % 12
    return (d, sign, m, s)


def signtostr(sign: int) -> str:
    """Return the English name of a zodiac sign 0..11."""
    s = int(sign) % 12
    return _RASI_NAMES[s]


def get_nakshatra_name(nak: int) -> str:
    """Return the standard Sanskrit name of a nakshatra 0..26."""
    n = int(nak) % 27
    return _NAKSHATRA_NAMES[n]


def house_system_name(code: str) -> str:
    """Return the human-readable name of a house system letter code.

    Accepts a single-character str (per the upstream contrib signature)."""
    if not isinstance(code, str):
        raise TypeError("argument 1 must be str, not " + type(code).__name__)
    return _HOUSE_SYSTEM_NAMES.get(code, "Unknown")


# =============================================================================
# Vedic dignities / strengths
# =============================================================================

# Domicile lordship: which planet rules each sign (Vedic tradition)
# Indices match RASI codes; values are Vedic planet ids.
_RASI_LORDS = (
    KUJA,    # Aries  - Mars
    SUKRA,   # Taurus - Venus
    BUDHA,   # Gemini - Mercury
    CHANDRA, # Cancer - Moon
    SURYA,   # Leo    - Sun
    BUDHA,   # Virgo  - Mercury
    SUKRA,   # Libra  - Venus
    KUJA,    # Scorpio- Mars
    GURU,    # Sagittarius - Jupiter
    SANI,    # Capricorn   - Saturn
    SANI,    # Aquarius    - Saturn
    GURU,    # Pisces - Jupiter
)


def lord(rasi: int) -> int:
    """Return the Vedic lord planet of a Rasi."""
    return _RASI_LORDS[rasinorm(rasi)]


# Exaltation degree of each planet (Vedic): planet -> (sign, degree_within_sign)
_EXALTATION = {
    SURYA: (ARIES, 10.0),
    CHANDRA: (TAURUS, 3.0),
    KUJA: (CAPRICORN, 28.0),
    BUDHA: (VIRGO, 15.0),
    GURU: (CANCER, 5.0),
    SUKRA: (PISCES, 27.0),
    SANI: (LIBRA, 20.0),
}

# Debilitation degree: 180° away from exaltation.
_DEBILITATION = {
    p: ((sign + 6) % 12, deg) for p, (sign, deg) in _EXALTATION.items()
}


def ochchabala(planet: int, longitude: float) -> float:
    """Return Ochchabala (exaltation strength) of a planet at a longitude.

    Ranges from 0 (at debilitation) to 60 (at exaltation), linear with
    180° between the two — the classical Vedic strength function.
    Returns 0.0 for nodes (Rahu/Ketu) and for unknown planet ids.
    """
    if planet not in _EXALTATION:
        return 0.0
    ex_sign, ex_deg = _EXALTATION[planet]
    ex_long = ex_sign * 30.0 + ex_deg
    diff = _normalize_longitude(longitude - ex_long)
    # 0° from exaltation -> 60. 180° from exaltation -> 0. Linear.
    if diff > 180.0:
        diff = 360.0 - diff
    return 60.0 * (1.0 - diff / 180.0)


def residential_strength(deg: float) -> float:
    """Return residential strength as a fraction (0..1) of the planet's
    progress through its current sign.

    Classical interpretation: a planet is strongest at the start of its
    sign and weakest at the end — this returns the inverse position
    (1.0 at 0°, 0.0 at 30°)."""
    within = _normalize_longitude(deg) % 30.0
    return 1.0 - within / 30.0


# Vedic naisargika (natural) friendships:
#   2 = friend, 1 = neutral, 0 = enemy
# Order: from <- to for the table key.
_NAISARGIKA = {
    SURYA:  {SURYA: 1, CHANDRA: 2, KUJA: 2, BUDHA: 1, GURU: 2, SUKRA: 0, SANI: 0},
    CHANDRA:{SURYA: 2, CHANDRA: 1, KUJA: 1, BUDHA: 2, GURU: 1, SUKRA: 1, SANI: 1},
    KUJA:   {SURYA: 2, CHANDRA: 2, KUJA: 1, BUDHA: 0, GURU: 2, SUKRA: 1, SANI: 1},
    BUDHA:  {SURYA: 2, CHANDRA: 0, KUJA: 1, BUDHA: 1, GURU: 1, SUKRA: 2, SANI: 1},
    GURU:   {SURYA: 2, CHANDRA: 2, KUJA: 2, BUDHA: 0, GURU: 1, SUKRA: 0, SANI: 1},
    SUKRA:  {SURYA: 0, CHANDRA: 0, KUJA: 1, BUDHA: 2, GURU: 1, SUKRA: 1, SANI: 2},
    SANI:   {SURYA: 0, CHANDRA: 0, KUJA: 0, BUDHA: 2, GURU: 1, SUKRA: 2, SANI: 1},
}


def naisargika_relation(p1: int, p2: int) -> int:
    """Return 0/1/2 (enemy/neutral/friend) for the natural relation
    between two Vedic planets."""
    try:
        return _NAISARGIKA[p1][p2]
    except KeyError:
        return 1


def tatkalika_relation(p1_house: int, p2_house: int) -> int:
    """Return the temporal (tatkalika) friendship between two planets,
    given their relative house separation. Planets in houses 2, 3, 4,
    10, 11, 12 from one another are temporal friends (2). Others are
    temporal enemies (0). Same house is neutral (1)."""
    d = (p2_house - p1_house) % 12
    if d == 0:
        return 1
    return 2 if d in {1, 2, 3, 9, 10, 11} else 0


def raman_houses(jd: float, geolon: float, geolat: float) -> tuple[float, ...]:
    """Return 12 house cusps using B.V. Raman's variant of the Sripati
    house system.

    This is a thin convenience wrapper around the canonical ``houses``
    function with the Sripati system code 'S'. ``houses`` returns the
    12 cusps as a 12-tuple (house 1 at index 0)."""
    from .. import houses as _houses

    cusps, _ascmc = _houses(jd, geolat, geolon, ord("S"))
    return tuple(cusps)


def saturn_4_stars(jd: float) -> tuple[float, float, float, float]:
    """Return the longitudes of the four mythological Saturn-related
    fixed stars at a given JD: Aldebaran, Regulus, Antares, Fomalhaut.
    These are the Vedic "Royal Stars" / Watchers of the Heavens.

    Uses the canonical ``fixstar_ut`` to obtain ecliptic longitude
    of each star."""
    import logging as _logging

    from .. import fixstar_ut as _fixstar_ut

    _logger = _logging.getLogger("libephemeris.contrib")

    out: list[float] = []
    for name in ("Aldebaran", "Regulus", "Antares", "Fomalhaut"):
        try:
            # fixstar_ut returns (pos, resolved_name, retflag)
            pos, _resolved, _retflag = _fixstar_ut(name, jd)
            out.append(float(pos[0]))
        except Exception as exc:
            # Log unexpected failure so callers can debug; the NaN return
            # remains the documented signal for "value unavailable".
            _logger.debug(
                "saturn_4_stars: fixstar_ut(%r) failed at jd=%.6f: %s",
                name, jd, exc,
            )
            out.append(_math.nan)
    return tuple(out)  # type: ignore[return-value]


# =============================================================================
# Time / JD helpers
# =============================================================================


def jdnow() -> float:
    """Return the current Julian Day (UT) at the moment of call."""
    now = _dt.datetime.now(_dt.timezone.utc)
    frac = (now.hour + now.minute / 60.0 + now.second / 3600.0
            + now.microsecond / 3.6e9)
    return _julday(now.year, now.month, now.day, frac)


def jduration(jd_start: float, jd_end: float) -> float:
    """Duration between two Julian Days in (fractional) days."""
    return float(jd_end - jd_start)


def t2i(time_string: str) -> tuple[int, int, int]:
    """Parse 'HH:MM:SS' (or 'HH:MM') into (hour, minute, second)."""
    parts = time_string.strip().split(":")
    h = int(parts[0])
    m = int(parts[1]) if len(parts) > 1 else 0
    s = int(parts[2]) if len(parts) > 2 else 0
    return (h, m, s)


def dt2i(date_string: str) -> tuple[int, int, int, int, int, int]:
    """Parse a 'YYYY-MM-DD HH:MM:SS' or ISO string into 6-tuple."""
    dt = _dt.datetime.fromisoformat(date_string.replace(" ", "T"))
    return (dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)


def jd2isostr(jd: float) -> str:
    """Convert a UT Julian Day to an ISO 8601 datetime string.

    Uses datetime.timedelta to compose the time-of-day so that
    a rounded-up second/minute/hour propagates cleanly across the
    minute/hour/day boundary (e.g., 23:59:59.9 → next day 00:00:00,
    never the malformed 23:59:60)."""
    y, mo, d, hour_frac = _revjul(jd, 1)
    base = _dt.datetime(y, mo, d)
    # Add the fractional hour, rounded to the nearest second
    delta = _dt.timedelta(seconds=round(hour_frac * 3600.0))
    moment = base + delta
    return moment.strftime("%Y-%m-%dT%H:%M:%S")


def years_diff(date_a: str, date_b: str) -> float:
    """Decimal years between two ISO date strings (b - a)."""
    da = _dt.datetime.fromisoformat(date_a.replace(" ", "T"))
    db = _dt.datetime.fromisoformat(date_b.replace(" ", "T"))
    return (db - da).total_seconds() / (365.25 * 86400.0)


# =============================================================================
# Coordinate helpers
# =============================================================================


def geoc2d(coord: str) -> float:
    """Parse a coordinate string in '12N34' / '12S34' / '12E34' / '12W34'
    format into signed decimal degrees. N/E positive, S/W negative."""
    s = coord.strip().upper()
    sign = 1.0
    for letter in "NSWE":
        if letter in s:
            if letter in "SW":
                sign = -1.0
            deg_part, min_part = s.split(letter, 1)
            deg = float(deg_part) if deg_part else 0.0
            mins = float(min_part) if min_part else 0.0
            return sign * (deg + mins / 60.0)
    return float(s)


def geolat2c(lat: float) -> str:
    """Decimal latitude → 'DDNMM' / 'DDSMM' string. Rounded minutes
    that hit 60 are carried into the degree component."""
    hem = "N" if lat >= 0 else "S"
    val = abs(lat)
    d = int(val)
    m = int(round((val - d) * 60.0))
    if m >= 60:
        m -= 60
        d += 1
    return f"{d}{hem}{m:02d}"


def geolon2c(lon: float) -> str:
    """Decimal longitude → 'DDDEMM' / 'DDDWMM' string. Rounded minutes
    that hit 60 are carried into the degree component."""
    hem = "E" if lon >= 0 else "W"
    val = abs(lon)
    d = int(val)
    m = int(round((val - d) * 60.0))
    if m >= 60:
        m -= 60
        d += 1
    return f"{d}{hem}{m:02d}"


# =============================================================================
# Aspect helpers
# =============================================================================


def _angular_separation(a: float, b: float) -> float:
    """Shortest angular separation between two longitudes, 0..180."""
    d = abs(_normalize_longitude(a) - _normalize_longitude(b))
    return d if d <= 180.0 else 360.0 - d


def match_aspect(
    a_id: int,
    a_lon: float,
    aspect: float,
    b_id: int,
    b_lon: float,
    orb: float,
) -> tuple[int, float]:
    """Test if two longitudes form a given aspect within an orb.

    Returns ``(is_match, exact_orb)`` where is_match is 1 if the actual
    angular separation is within `orb` degrees of the target `aspect`,
    else 0. exact_orb is the actual deviation."""
    sep = _angular_separation(a_lon, b_lon)
    delta = abs(sep - aspect)
    return (1 if delta <= orb else 0, delta)


def match_aspect2(
    a_lon: float,
    b_lon: float,
    aspect: float,
    orb: float,
) -> tuple[int, float]:
    """Simplified two-longitude match_aspect (no planet ids)."""
    sep = _angular_separation(a_lon, b_lon)
    delta = abs(sep - aspect)
    return (1 if delta <= orb else 0, delta)


def match_aspect3(
    a_id: int,
    a_lon: float,
    aspect: float,
    b_id: int,
    b_lon: float,
    orb: float,
) -> tuple[int, float, float]:
    """match_aspect plus a third value (currently the signed delta)."""
    is_match, delta = match_aspect(a_id, a_lon, aspect, b_id, b_lon, orb)
    sep_signed = (b_lon - a_lon) % 360.0
    if sep_signed > 180.0:
        sep_signed -= 360.0
    return (is_match, delta, sep_signed)


def match_aspect4(
    a_id: int,
    a_lon: float,
    aspect: float,
    b_id: int,
    b_lon: float,
    orb: float,
) -> tuple[int, float, float, float]:
    """match_aspect3 plus the absolute separation."""
    is_match, delta, signed = match_aspect3(a_id, a_lon, aspect, b_id, b_lon, orb)
    sep = _angular_separation(a_lon, b_lon)
    return (is_match, delta, signed, sep)


def _walk_until_aspect(
    jd: float,
    body_a: int,
    body_b_lon: float,
    aspect: float,
    orb: float,
    step_days: float,
    max_days: float,
) -> tuple[float, float]:
    """Step body_a's longitude forward in time until it forms `aspect`
    with `body_b_lon` (taken as fixed). Returns (jd_when, exact_delta).

    Uses calc_ut to evaluate body_a; b is treated as a fixed reference
    longitude. For two moving bodies, see next_aspect_with."""
    t = jd
    while t - jd < max_days:
        pos, _ = _calc_ut(t, body_a, _FLG_SPEED)
        delta = abs(_angular_separation(pos[0], body_b_lon) - aspect)
        if delta <= orb:
            return t, delta
        t += step_days
    return _math.nan, _math.inf


def next_aspect(
    body: int,
    aspect: float,
    target_lon: float,
    jd_start: float,
    orb: float = 1.0,
) -> tuple[float, float]:
    """Find the next time `body` makes `aspect` with a fixed longitude.

    Returns (jd_when, exact_orb). Returns (nan, inf) if not found
    within ~3 years.

    Precision note:
        Uses a fixed 0.5-day step; the result is the first sample
        whose orb is within ``orb`` degrees of the target aspect, not
        the bisected moment of exact match. For fast-moving bodies
        (Moon, inner planets near a station) a small ``orb`` together
        with the coarse step can skip past the in-orb window: pick
        ``orb`` larger than ``body_speed * 0.5`` if exact-match timing
        matters, or post-process the returned JD with a bisection of
        your own.
    """
    return _walk_until_aspect(
        jd_start, body, target_lon, aspect, orb, 0.5, 365.25 * 3
    )


def next_aspect2(
    body: int,
    aspect: float,
    target_lon: float,
    jd_start: float,
) -> tuple[float, float]:
    """next_aspect with default orb of 1°."""
    return next_aspect(body, aspect, target_lon, jd_start, 1.0)


def next_aspect_cusp(
    body: int,
    aspect: float,
    cusp_lon: float,
    jd_start: float,
    orb: float = 1.0,
) -> tuple[float, float]:
    """Next aspect of `body` with a house cusp at `cusp_lon`."""
    return next_aspect(body, aspect, cusp_lon, jd_start, orb)


def next_aspect_cusp2(
    body: int,
    aspect: float,
    cusp_lon: float,
    jd_start: float,
) -> tuple[float, float]:
    """next_aspect_cusp with default orb of 1°."""
    return next_aspect_cusp(body, aspect, cusp_lon, jd_start, 1.0)


def next_aspect_with(
    body_a: int,
    body_b: int,
    aspect: float,
    jd_start: float,
    orb: float = 1.0,
) -> tuple[float, float]:
    """Next time bodies `body_a` and `body_b` form `aspect`.

    Both bodies move; steps forward by 0.5 days until match or 3 years.

    Precision note:
        Same coarse-step limitation as :func:`next_aspect`. For
        sub-hour timing accuracy (e.g. Moon-Sun conjunction), bisect
        around the returned JD or pass a generous ``orb``.
    """
    t = jd_start
    end = jd_start + 365.25 * 3
    while t < end:
        pos_a, _ = _calc_ut(t, body_a, _FLG_SPEED)
        pos_b, _ = _calc_ut(t, body_b, _FLG_SPEED)
        delta = abs(_angular_separation(pos_a[0], pos_b[0]) - aspect)
        if delta <= orb:
            return t, delta
        t += 0.5
    return _math.nan, _math.inf


def next_aspect_with2(
    body_a: int,
    body_b: int,
    aspect: float,
    jd_start: float,
) -> tuple[float, float]:
    """next_aspect_with with default orb of 1°."""
    return next_aspect_with(body_a, body_b, aspect, jd_start, 1.0)


# =============================================================================
# Retrograde / motion helpers
# =============================================================================


def next_retro(
    body: int,
    jd_start: float,
    max_days: float = 365.25 * 3,
) -> float:
    """Find the next JD at which `body` turns retrograde.

    A body is retrograde when its ecliptic-longitude speed is negative.
    Steps forward by 0.5 days; returns the JD where the speed first
    becomes negative. Returns NaN if not found within `max_days`.

    Precision note:
        With a fixed 0.5-day step, fast-moving stations (Mercury) can
        be detected up to ~12 hours late. For sub-day station timing,
        bisect the returned JD against ``calc_ut(..., FLG_SPEED)`` to
        refine to the speed zero-crossing.
    """
    SPEED_FLAG = _FLG_SPEED
    t = jd_start
    last_speed: float | None = None
    while t - jd_start < max_days:
        pos, _ = _calc_ut(t, body, SPEED_FLAG)
        speed = pos[3]  # longitude speed
        if last_speed is not None and last_speed >= 0 and speed < 0:
            return t
        last_speed = speed
        t += 0.5
    return _math.nan


# =============================================================================
# Database / atlas / timezone (not implemented)
# =============================================================================


def _not_implemented(name: str) -> _Any:
    def _stub(*args, **kwargs):  # noqa: ARG001
        raise NotImplementedError(
            f"libephemeris.contrib.{name}() is not implemented. "
            "This function requires the optional SQLite atlas/timezone "
            "databases shipped with the upstream reference distribution. "
            "If you need it, install the upstream reference Python binding "
            f"and call its contrib.{name}() directly."
        )

    _stub.__name__ = name
    _stub.__qualname__ = name
    return _stub


atlas_connect = _not_implemented("atlas_connect")
atlas_close = _not_implemented("atlas_close")
atlas_search = _not_implemented("atlas_search")
atlas_countries_list = _not_implemented("atlas_countries_list")
db_connect = _not_implemented("db_connect")
db_close = _not_implemented("db_close")
tzabbr_find = _not_implemented("tzabbr_find")
tzabbr_list = _not_implemented("tzabbr_list")
