# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Reference API-compatible constants for libephemeris.

This module defines all constants used for planetary calculations, including:
- Planet/Body IDs: Numeric identifiers for celestial bodies
- Calculation Flags: Bitwise flags controlling observation parameters
- Sidereal Modes: Ayanamsha systems for sidereal astrology
- Calendar Systems: Julian vs Gregorian calendar selection
- Eclipse Types: Classification of solar and lunar eclipses

Constants are organized into logical groups for easy navigation.
Values match reference API v2.x for API compatibility.
"""

from __future__ import annotations

# __all__ is generated dynamically at the bottom of this module to export
# all public constants and functions while excluding private names and imports.

# =============================================================================
# PLANET AND BODY IDENTIFIERS
# =============================================================================

# Special values
ECL_NUT: int = -1  # Nutation and obliquity calculation

# Major planets (traditional + modern)
SUN: int = 0
MOON: int = 1
MERCURY: int = 2
VENUS: int = 3
MARS: int = 4
JUPITER: int = 5
SATURN: int = 6
URANUS: int = 7
NEPTUNE: int = 8
PLUTO: int = 9

# Planet aliases without SE_ prefix for reference API compatibility

MEAN_NODE: int = 10  # Mean lunar node (Dragon's Head)
TRUE_NODE: int = 11  # True (osculating) lunar node
MEAN_APOG: int = 12  # Mean lunar apogee (Black Moon Lilith)
OSCU_APOG: int = 13  # Osculating (true) lunar apogee

# Earth and centaurs
EARTH: int = 14
CHIRON: int = 15  # Centaur between Saturn and Uranus
PHOLUS: int = 16  # Centaur beyond Saturn

# Main belt asteroids
CERES: int = 17
PALLAS: int = 18
JUNO: int = 19
VESTA: int = 20

# Interpolated lunar apsides
INTP_APOG: int = 21  # Interpolated apogee
INTP_PERG: int = 22  # Interpolated perigee

# Count and offsets
NPLANETS: int = 23  # Total number of standard planet IDs
AST_OFFSET: int = 10000  # Offset for asteroid catalog numbers
VARUNA: int = AST_OFFSET + 20000  # TNO Varuna
FICT_OFFSET: int = 40  # Offset for fictitious bodies
NFICT_ELEM: int = 15  # Number of fictitious elements
COMET_OFFSET: int = 1000  # Offset for comet IDs

# Hamburg School Uranian planets (fictitious bodies)
CUPIDO: int = FICT_OFFSET + 0  # 40 - First Uranian planet (Cupido)
HADES: int = FICT_OFFSET + 1  # 41 - Second Uranian planet (Hades)
ZEUS: int = FICT_OFFSET + 2  # 42 - Third Uranian planet (Zeus)
KRONOS: int = FICT_OFFSET + 3  # 43 - Fourth Uranian planet (Kronos)
APOLLON: int = FICT_OFFSET + 4  # 44 - Fifth Uranian planet (Apollon)
ADMETOS: int = FICT_OFFSET + 5  # 45 - Sixth Uranian planet (Admetos)
VULKANUS: int = FICT_OFFSET + 6  # 46 - Seventh Uranian planet (Vulkanus)
POSEIDON: int = FICT_OFFSET + 7  # 47 - Eighth Uranian planet (Poseidon)

# Reference API-compatible aliases for Uranian planets

ISIS: int = FICT_OFFSET + 8  # 48 - Transpluto/Isis
TRANSPLUTO: int = ISIS  # Alias for ISIS

# Vulcan - hypothetical intramercurial planet
VULCAN: int = FICT_OFFSET + 15  # 55 - Intramercurial hypothetical planet

# Waldemath - Dr. Waldemath's hypothetical second moon of Earth
# Note: This is different from Mean Lilith and True Lilith which are lunar apogee points
WALDEMATH: int = FICT_OFFSET + 18  # 58 - Waldemath's hypothetical Dark Moon

# Planet X Leverrier - Leverrier's calculated "Planet X" that led to Neptune's discovery
# This is the orbital elements Leverrier computed in 1846 to predict the position of Neptune.
# Historical note: Leverrier called it "Planet X" before Neptune was discovered at that position.
PLANET_X_LEVERRIER: int = (
    FICT_OFFSET + 11
)  # 51 - Leverrier's Planet X (Neptune prediction)

# Planet X Adams - Adams' calculated "Planet X" (independently derived, similar to Leverrier's)
# John Couch Adams independently predicted Neptune's position around the same time as Leverrier.
# This uses Adams' orbital elements prediction (FICT_OFFSET + 12 = 52).
PLANET_X_ADAMS: int = (
    FICT_OFFSET + 12
)  # 52 - Adams' Planet X (Neptune prediction)

# Planet X Lowell - Percival Lowell's predicted "Planet X" that led to Pluto's discovery
# Lowell predicted a trans-Neptunian planet based on perceived perturbations in Uranus's orbit.
# The search based on his predictions eventually led to Clyde Tombaugh's discovery of Pluto in 1930,
# though Pluto was too small to be Lowell's predicted Planet X.
# This uses Lowell's orbital elements prediction (FICT_OFFSET + 13 = 53).
# Orbital elements (1915 prediction): a=43.0 AU, e=0.202, i=10°
PLANET_X_LOWELL: int = (
    FICT_OFFSET + 13
)  # 53 - Lowell's Planet X (Pluto prediction)

# Planet X Pickering - William H. Pickering's predicted "Planet O" (1919)
# Pickering proposed several trans-Neptunian planets (Planet O, P, Q, R, S, T, U).
# His most famous prediction, "Planet O" (1919), had the following orbital elements:
# - Semi-major axis: 51.9 AU
# - Eccentricity: 0.31
# - Inclination: 15°
# - Orbital period: ~373.5 years
# Like Lowell's Planet X, these predictions were based on supposed perturbations in outer
# planet orbits, which later proved to be observational errors.
# This uses Pickering's orbital elements prediction (FICT_OFFSET + 14 = 54).
PLANET_X_PICKERING: int = (
    FICT_OFFSET + 14
)  # 54 - Pickering's Planet O/X prediction

# White Moon (Selena) - Point opposite to Black Moon Lilith (lunar perigee = apogee + 180°)
# Calculated as Mean Lilith + 180° (i.e., the mean lunar perigee)
# Note: Some systems use True Lilith + 180° instead; libephemeris supports both via calc_white_moon_position()
WHITE_MOON: int = (
    FICT_OFFSET + 16
)  # 56 - White Moon Selena (opposite to Black Moon Lilith)

SELENA: int = WHITE_MOON  # Alias - Selena is another name for White Moon

# Proserpina - hypothetical trans-Plutonian planet used by some astrologers
# This is a hypothetical body not in the standard fictitious bodies set
# Orbital elements are based on traditional astrological sources
PROSERPINA: int = FICT_OFFSET + 17  # 57 - Hypothetical trans-Plutonian planet

_NALL_NAT_POINTS_INTERNAL: int = NPLANETS + NFICT_ELEM + AST_OFFSET + COMET_OFFSET

# Trans-Neptunian Objects (TNOs) - Catalog number + offset
ERIS: int = 136199 + AST_OFFSET  # Largest known dwarf planet
SEDNA: int = 90377 + AST_OFFSET  # Detached TNO
HAUMEA: int = 136108 + AST_OFFSET  # Fast-rotating dwarf planet
MAKEMAKE: int = 136472 + AST_OFFSET  # Classical Kuiper belt object
IXION: int = 28978 + AST_OFFSET  # Plutino
ORCUS: int = 90482 + AST_OFFSET  # Plutino, "anti-Pluto"
QUAOAR: int = 50000 + AST_OFFSET  # Classical KBO
NESSUS: int = 7066 + AST_OFFSET  # Centaur, astrologically important
ASBOLUS: int = 8405 + AST_OFFSET  # Centaur, astrologically significant
CHARIKLO: int = 10199 + AST_OFFSET  # Centaur, largest known, has ring system
GONGGONG: int = (
    225088 + AST_OFFSET
)  # TNO, dwarf planet candidate (formerly 2007 OR10)
APOPHIS: int = 99942 + AST_OFFSET  # Near-Earth asteroid, close approach 2029
HYGIEA: int = 10 + AST_OFFSET  # Fourth largest asteroid, dwarf planet candidate
INTERAMNIA: int = 704 + AST_OFFSET  # Fifth largest asteroid, main belt
DAVIDA: int = 511 + AST_OFFSET  # Seventh largest asteroid, main belt
EUROPA_AST: int = (
    52 + AST_OFFSET
)  # 52 Europa (main belt asteroid, not Jupiter's moon)
SYLVIA: int = (
    87 + AST_OFFSET
)  # 87 Sylvia (triple asteroid system with moons Romulus and Remus)
PSYCHE: int = (
    16 + AST_OFFSET
)  # 16 Psyche (metallic M-type asteroid, NASA Psyche mission target)
EROS: int = (
    433 + AST_OFFSET
)  # 433 Eros (near-Earth asteroid, NEAR Shoemaker mission target)
AMOR: int = (
    1221 + AST_OFFSET
)  # 1221 Amor (prototype of the Amor near-Earth asteroid class)
ICARUS: int = (
    1566 + AST_OFFSET
)  # 1566 Icarus (Apollo asteroid, highly eccentric orbit reaches inside Mercury)
TORO: int = 1685 + AST_OFFSET  # 1685 Toro (Apollo asteroid, near-Earth asteroid)
SAPPHO: int = (
    80 + AST_OFFSET
)  # 80 Sappho (main belt asteroid, artistic expression and same-sex love)
PANDORA_AST: int = (
    55 + AST_OFFSET
)  # 55 Pandora (main belt asteroid, distinct from Saturn moon Pandora)
LILITH_AST: int = (
    1181 + AST_OFFSET
)  # 1181 Lilith (main belt asteroid, not to be confused with lunar apogee Lilith)
HIDALGO: int = (
    944 + AST_OFFSET
)  # 944 Hidalgo (Centaur-class asteroid with comet-like orbit, astrological research)
TOUTATIS: int = (
    4179 + AST_OFFSET
)  # 4179 Toutatis (Apollo PHA, radar and spacecraft target, tumbling rotation)
ITOKAWA: int = (
    25143 + AST_OFFSET
)  # 25143 Itokawa (Apollo PHA, Hayabusa sample return mission target)
BENNU: int = (
    101955 + AST_OFFSET
)  # 101955 Bennu (Apollo PHA, OSIRIS-REx sample return mission target)
RYUGU: int = (
    162173 + AST_OFFSET
)  # 162173 Ryugu (Apollo PHA, Hayabusa2 sample return mission target)

# =============================================================================
# NAIF IDS FOR SPK KERNELS
# =============================================================================
# NAIF IDs used in JPL SPK kernels for minor bodies.
# Convention: For numbered asteroids, NAIF ID = asteroid_number + 2000000
# These constants simplify registration of SPK bodies.

NAIF_ASTEROID_OFFSET: int = 2000000  # Add asteroid number to get NAIF ID

# Common minor body NAIF IDs (asteroid_number + 2000000)
NAIF_CERES: int = 2000001  # 1 Ceres
NAIF_PALLAS: int = 2000002  # 2 Pallas
NAIF_JUNO: int = 2000003  # 3 Juno
NAIF_VESTA: int = 2000004  # 4 Vesta
NAIF_CHIRON: int = 2002060  # 2060 Chiron
NAIF_PHOLUS: int = 2005145  # 5145 Pholus
NAIF_ERIS: int = 2136199  # 136199 Eris
NAIF_SEDNA: int = 2090377  # 90377 Sedna
NAIF_HAUMEA: int = 2136108  # 136108 Haumea
NAIF_MAKEMAKE: int = 2136472  # 136472 Makemake
NAIF_IXION: int = 2028978  # 28978 Ixion
NAIF_ORCUS: int = 2090482  # 90482 Orcus
NAIF_QUAOAR: int = 2050000  # 50000 Quaoar
NAIF_NESSUS: int = 2007066  # 7066 Nessus (Centaur)
NAIF_ASBOLUS: int = 2008405  # 8405 Asbolus (Centaur)
NAIF_CHARIKLO: int = 2010199  # 10199 Chariklo (Centaur, largest, has rings)
NAIF_GONGGONG: int = 2225088  # 225088 Gonggong (TNO, dwarf planet candidate)
NAIF_APOPHIS: int = 2099942  # 99942 Apophis (Near-Earth asteroid)
NAIF_HYGIEA: int = (
    2000010  # 10 Hygiea (fourth largest asteroid, dwarf planet candidate)
)
NAIF_EROS: int = 2000433  # 433 Eros (near-Earth asteroid, NEAR Shoemaker mission)
NAIF_PSYCHE: int = 2000016  # 16 Psyche (metallic M-type asteroid)
NAIF_EUROPA_AST: int = 2000052  # 52 Europa (main belt asteroid)
NAIF_PANDORA_AST: int = 2000055  # 55 Pandora (main belt asteroid)
NAIF_SAPPHO: int = 2000080  # 80 Sappho (main belt asteroid)
NAIF_SYLVIA: int = 2000087  # 87 Sylvia (triple asteroid system)
NAIF_DAVIDA: int = 2000511  # 511 Davida (main belt asteroid)
NAIF_INTERAMNIA: int = 2000704  # 704 Interamnia (main belt asteroid)
NAIF_HIDALGO: int = 2000944  # 944 Hidalgo (centaur-class asteroid)
NAIF_LILITH_AST: int = 2001181  # 1181 Lilith (main belt asteroid)
NAIF_AMOR: int = 2001221  # 1221 Amor (near-Earth asteroid)
NAIF_ICARUS: int = 2001566  # 1566 Icarus (Apollo asteroid)
NAIF_TORO: int = 2001685  # 1685 Toro (Apollo asteroid)
NAIF_TOUTATIS: int = 2004179  # 4179 Toutatis (Apollo PHA)
NAIF_ITOKAWA: int = 2025143  # 25143 Itokawa (Apollo asteroid)
NAIF_BENNU: int = 2101955  # 101955 Bennu (Apollo asteroid)
NAIF_RYUGU: int = 2162173  # 162173 Ryugu (Apollo asteroid)

# =============================================================================
# SPK BODY NAME MAPPING
# =============================================================================
# Mapping from libephemeris body IDs (SE_*) to JPL Horizons target designations.
# This is used for automatic SPK downloads. The tuple contains:
#   (horizons_id: str, naif_id: int)
# where horizons_id is the identifier used in JPL Horizons queries (typically
# the asteroid catalog number), and naif_id is the NAIF SPICE ID.
#
# Note: Some bodies in this map cannot be auto-downloaded from JPL Horizons.
# See SPK_AUTO_DOWNLOAD_BLOCKED below for the list and reasons.

SPK_BODY_NAME_MAP: dict[int, tuple[str, int]] = {
    CHIRON: ("2060", NAIF_CHIRON),  # 2060 Chiron (centaur)
    PHOLUS: ("5145", NAIF_PHOLUS),  # 5145 Pholus (centaur)
    CERES: (
        "Ceres;",
        NAIF_CERES,
    ),  # 1 Ceres (dwarf planet) - use name; to bypass major body index
    PALLAS: (
        "Pallas;",
        NAIF_PALLAS,
    ),  # 2 Pallas (main belt asteroid) - use name; to bypass major body index
    JUNO: (
        "Juno;",
        NAIF_JUNO,
    ),  # 3 Juno (main belt asteroid) - use name; to bypass major body index
    VESTA: (
        "Vesta;",
        NAIF_VESTA,
    ),  # 4 Vesta (main belt asteroid) - use name; to bypass major body index
    ERIS: ("136199", NAIF_ERIS),  # 136199 Eris (dwarf planet)
    SEDNA: ("90377", NAIF_SEDNA),  # 90377 Sedna (detached TNO)
    HAUMEA: ("136108", NAIF_HAUMEA),  # 136108 Haumea (dwarf planet)
    MAKEMAKE: ("136472", NAIF_MAKEMAKE),  # 136472 Makemake (dwarf planet)
    IXION: ("28978", NAIF_IXION),  # 28978 Ixion (plutino)
    ORCUS: ("90482", NAIF_ORCUS),  # 90482 Orcus (plutino)
    QUAOAR: ("50000", NAIF_QUAOAR),  # 50000 Quaoar (classical KBO)
    VARUNA: ("20000", 2020000),  # 20000 Varuna (classical KBO)
    NESSUS: ("7066", NAIF_NESSUS),  # 7066 Nessus (centaur)
    ASBOLUS: ("8405", NAIF_ASBOLUS),  # 8405 Asbolus (centaur)
    CHARIKLO: ("10199", NAIF_CHARIKLO),  # 10199 Chariklo (centaur, largest, rings)
    GONGGONG: (
        "225088",
        NAIF_GONGGONG,
    ),  # 225088 Gonggong (TNO, dwarf planet candidate)
    APOPHIS: ("99942", NAIF_APOPHIS),  # 99942 Apophis (Near-Earth asteroid)
    HYGIEA: (
        "Hygiea;",
        NAIF_HYGIEA,
    ),  # 10 Hygiea (fourth largest asteroid) - use name; to bypass major body index
    EROS: ("433", NAIF_EROS),  # 433 Eros (near-Earth asteroid, NEAR Shoemaker)
    PSYCHE: ("16", NAIF_PSYCHE),  # 16 Psyche (metallic M-type asteroid)
    EUROPA_AST: ("52", NAIF_EUROPA_AST),  # 52 Europa (main belt asteroid)
    PANDORA_AST: ("55", NAIF_PANDORA_AST),  # 55 Pandora (main belt asteroid)
    SAPPHO: ("80", NAIF_SAPPHO),  # 80 Sappho (main belt asteroid)
    SYLVIA: ("87", NAIF_SYLVIA),  # 87 Sylvia (triple asteroid system)
    DAVIDA: (
        "Davida;",
        NAIF_DAVIDA,
    ),  # 511 Davida (main belt asteroid) - use name; to bypass major body index
    INTERAMNIA: (
        "Interamnia;",
        NAIF_INTERAMNIA,
    ),  # 704 Interamnia (main belt asteroid) - use name; to bypass major body index
    HIDALGO: ("944", NAIF_HIDALGO),  # 944 Hidalgo (centaur-class asteroid)
    LILITH_AST: ("1181", NAIF_LILITH_AST),  # 1181 Lilith (main belt asteroid)
    AMOR: ("1221", NAIF_AMOR),  # 1221 Amor (near-Earth asteroid)
    ICARUS: ("1566", NAIF_ICARUS),  # 1566 Icarus (Apollo asteroid)
    TORO: ("1685", NAIF_TORO),  # 1685 Toro (Apollo asteroid)
    TOUTATIS: ("4179", NAIF_TOUTATIS),  # 4179 Toutatis (Apollo PHA)
    ITOKAWA: ("25143", NAIF_ITOKAWA),  # 25143 Itokawa (Apollo asteroid)
    BENNU: (
        "Bennu;",
        NAIF_BENNU,
    ),  # 101955 Bennu (Apollo asteroid) - use name; to bypass major body index
    RYUGU: ("162173", NAIF_RYUGU),  # 162173 Ryugu (Apollo asteroid)
}

def get_horizons_id(ipl: int) -> str | None:
    """
    Get the JPL Horizons target identifier for a libephemeris body ID.

    Args:
        ipl: libephemeris body ID (e.g., CHIRON, ERIS)

    Returns:
        The Horizons target identifier as a string (e.g., "2060" for Chiron),
        or None if the body is not in the mapping.

    Example:
        >>> from libephemeris.constants import get_horizons_id, CHIRON
        >>> get_horizons_id(CHIRON)
        '2060'
    """
    if ipl in SPK_BODY_NAME_MAP:
        return SPK_BODY_NAME_MAP[ipl][0]
    return None

def get_naif_id_from_ipl(ipl: int) -> int | None:
    """
    Get the NAIF SPICE ID for a libephemeris body ID.

    Args:
        ipl: libephemeris body ID (e.g., CHIRON, ERIS)

    Returns:
        The NAIF ID (e.g., 2002060 for Chiron), or None if the body
        is not in the mapping.

    Example:
        >>> from libephemeris.constants import get_naif_id_from_ipl, CHIRON
        >>> get_naif_id_from_ipl(CHIRON)
        2002060
    """
    if ipl in SPK_BODY_NAME_MAP:
        return SPK_BODY_NAME_MAP[ipl][1]
    return None

def get_spk_body_info_from_map(ipl: int) -> tuple[str, int] | None:
    """
    Get both Horizons ID and NAIF ID for a libephemeris body.

    Args:
        ipl: libephemeris body ID (e.g., CHIRON, ERIS)

    Returns:
        A tuple of (horizons_id, naif_id), or None if the body
        is not in the mapping.

    Example:
        >>> from libephemeris.constants import get_spk_body_info_from_map, ERIS
        >>> get_spk_body_info_from_map(ERIS)
        ('136199', 2136199)
    """
    return SPK_BODY_NAME_MAP.get(ipl)

# =============================================================================
# SPK AUTO-DOWNLOAD BLOCKED BODIES
# =============================================================================
# Bodies that are in SPK_BODY_NAME_MAP (they have valid Horizons/NAIF IDs and
# can accept manually-registered mission-specific kernels) but for which JPL
# Horizons blocks automatic SPK generation.
#
# At runtime these bodies skip the auto-download attempt and the strict-
# precision check, allowing the fallback chain (ASSIST / Keplerian) to proceed.
#
# The dict maps SE_* body ID -> human-readable reason string.
# Adding a new body is a one-line change.

SPK_AUTO_DOWNLOAD_BLOCKED: dict[int, str] = {
    BENNU: (
        "JPL Horizons blocks SPK generation for Bennu. "
        "Use ASSIST n-body or register a mission-specific kernel "
        "(e.g. bennu_refdrmc_v1.bsp from OSIRIS-REx)."
    ),
}

def is_spk_auto_download_blocked(ipl: int) -> bool:
    """Return True if JPL Horizons blocks SPK generation for this body.

    Bodies in this set are still in SPK_BODY_NAME_MAP (they have valid
    Horizons/NAIF identifiers for manual registration) but cannot be
    auto-downloaded via the standard Horizons SPK generation path.

    Args:
        ipl: libephemeris body ID (e.g., BENNU)

    Returns:
        True if auto-download is blocked for this body.
    """
    return ipl in SPK_AUTO_DOWNLOAD_BLOCKED

# =============================================================================
# REQUIRED SPK BODIES FOR HIGH-PRECISION CALCULATIONS
# =============================================================================
# Set of body IDs that require SPK kernels for high-precision ephemeris
# calculations. These bodies benefit significantly from SPK data:
#
# - Main belt asteroids (Ceres, Pallas, Juno, Vesta): Keplerian elements give
#   ~10-30 arcsec precision, SPK gives sub-arcsecond precision
# - Centaurs (Chiron, Pholus, Nessus): Highly perturbed orbits need SPK for
#   accurate positions over multi-decade timescales
# - Dwarf planets (Eris): Distant TNOs require SPK for research-grade precision
#
# Usage:
#     for body in REQUIRED_SPK_BODIES:
#         eph.ensure_major_asteroid_spk(body)
#
# Note: All bodies listed here are in SPK_BODY_NAME_MAP and can be downloaded
# via auto_download_asteroid_spk() or download_and_register_spk().

REQUIRED_SPK_BODIES: frozenset[int] = frozenset(
    {
        CHIRON,  # 2060 Chiron - Centaur, frequently used in astrology
        CERES,  # 1 Ceres - Dwarf planet, largest main belt asteroid
        PALLAS,  # 2 Pallas - Second largest main belt asteroid
        JUNO,  # 3 Juno - Main belt asteroid
        VESTA,  # 4 Vesta - Second most massive main belt asteroid
        PHOLUS,  # 5145 Pholus - Centaur with highly eccentric orbit
        NESSUS,  # 7066 Nessus - Centaur, astrologically significant
        ERIS,  # 136199 Eris - Largest known dwarf planet
    }
)

# =============================================================================
# VIRTUAL POINTS AND CALCULATED POSITIONS
# =============================================================================

# Fixed Stars (high offset to avoid ID collisions)
FIXSTAR_OFFSET: int = 1000000
REGULUS: int = FIXSTAR_OFFSET + 1  # Alpha Leonis
SPICA_STAR: int = FIXSTAR_OFFSET + 2  # Alpha Virginis
ALGOL: int = FIXSTAR_OFFSET + 3  # Beta Persei
SIRIUS: int = FIXSTAR_OFFSET + 4  # Alpha Canis Majoris
ALDEBARAN: int = FIXSTAR_OFFSET + 5  # Alpha Tauri
ANTARES: int = FIXSTAR_OFFSET + 6  # Alpha Scorpii
VEGA: int = FIXSTAR_OFFSET + 7  # Alpha Lyrae
POLARIS: int = FIXSTAR_OFFSET + 8  # Alpha Ursae Minoris
FOMALHAUT: int = FIXSTAR_OFFSET + 9  # Alpha Piscis Austrini
BETELGEUSE: int = FIXSTAR_OFFSET + 10  # Alpha Orionis
RIGEL: int = FIXSTAR_OFFSET + 11  # Beta Orionis
PROCYON: int = FIXSTAR_OFFSET + 12  # Alpha Canis Minoris
CAPELLA: int = FIXSTAR_OFFSET + 13  # Alpha Aurigae
ARCTURUS: int = FIXSTAR_OFFSET + 14  # Alpha Bootis
DENEB: int = FIXSTAR_OFFSET + 15  # Alpha Cygni
POLLUX: int = FIXSTAR_OFFSET + 16  # Beta Geminorum
CASTOR: int = FIXSTAR_OFFSET + 17  # Alpha Geminorum
ALTAIR: int = FIXSTAR_OFFSET + 18  # Alpha Aquilae
ACHERNAR: int = FIXSTAR_OFFSET + 19  # Alpha Eridani
CANOPUS: int = FIXSTAR_OFFSET + 20  # Alpha Carinae
ACRUX: int = FIXSTAR_OFFSET + 21  # Alpha Crucis
MIMOSA: int = FIXSTAR_OFFSET + 22  # Beta Crucis
GACRUX: int = FIXSTAR_OFFSET + 23  # Gamma Crucis
HADAR: int = FIXSTAR_OFFSET + 24  # Beta Centauri
RIGIL_KENT: int = FIXSTAR_OFFSET + 25  # Alpha Centauri
SHAULA: int = FIXSTAR_OFFSET + 26  # Lambda Scorpii
BELLATRIX: int = FIXSTAR_OFFSET + 27  # Gamma Orionis
ELNATH: int = FIXSTAR_OFFSET + 28  # Beta Tauri
MIRA: int = FIXSTAR_OFFSET + 29  # Omicron Ceti
ALNILAM: int = FIXSTAR_OFFSET + 30  # Epsilon Orionis
ALNITAK: int = FIXSTAR_OFFSET + 31  # Zeta Orionis
MINTAKA: int = FIXSTAR_OFFSET + 32  # Delta Orionis
SAIPH: int = FIXSTAR_OFFSET + 33  # Kappa Orionis
DIPHDA: int = FIXSTAR_OFFSET + 34  # Beta Ceti
ALPHARD: int = FIXSTAR_OFFSET + 35  # Alpha Hydrae
RASALHAGUE: int = FIXSTAR_OFFSET + 36  # Alpha Ophiuchi
ETAMIN: int = FIXSTAR_OFFSET + 37  # Gamma Draconis
KOCHAB: int = FIXSTAR_OFFSET + 38  # Beta Ursae Minoris
ALKAID: int = FIXSTAR_OFFSET + 39  # Eta Ursae Majoris
DUBHE: int = FIXSTAR_OFFSET + 40  # Alpha Ursae Majoris
MERAK: int = FIXSTAR_OFFSET + 41  # Beta Ursae Majoris
ALIOTH: int = FIXSTAR_OFFSET + 42  # Epsilon Ursae Majoris
MIZAR: int = FIXSTAR_OFFSET + 43  # Zeta Ursae Majoris
ALCOR: int = FIXSTAR_OFFSET + 44  # 80 Ursae Majoris
VINDEMIATRIX: int = FIXSTAR_OFFSET + 45  # Epsilon Virginis
ZUBENELGENUBI: int = FIXSTAR_OFFSET + 46  # Alpha Librae
ZUBENESCHAMALI: int = FIXSTAR_OFFSET + 47  # Beta Librae
UNUKALHAI: int = FIXSTAR_OFFSET + 48  # Alpha Serpentis
ALGIEBA: int = FIXSTAR_OFFSET + 49  # Gamma Leonis
DENEBOLA: int = FIXSTAR_OFFSET + 50  # Beta Leonis
MARKAB: int = FIXSTAR_OFFSET + 51  # Alpha Pegasi
SCHEAT: int = FIXSTAR_OFFSET + 52  # Beta Pegasi
ALCYONE: int = FIXSTAR_OFFSET + 53  # Eta Tauri (Pleiades) - Behenian star
ALGORAB: int = FIXSTAR_OFFSET + 54  # Delta Corvi - Behenian star
ALPHECCA: int = FIXSTAR_OFFSET + 55  # Alpha Coronae Borealis - Behenian star
DENEB_ALGEDI: int = FIXSTAR_OFFSET + 56  # Delta Capricorni - Behenian star

# Pleiades cluster stars (visible members)
ASTEROPE: int = FIXSTAR_OFFSET + 57  # 21 Tauri - Pleiades member
CELAENO: int = FIXSTAR_OFFSET + 58  # 16 Tauri - Pleiades member
ELECTRA: int = FIXSTAR_OFFSET + 59  # 17 Tauri - Pleiades member
MAIA: int = FIXSTAR_OFFSET + 60  # 20 Tauri - Pleiades member
MEROPE: int = FIXSTAR_OFFSET + 61  # 23 Tauri - Pleiades member
TAYGETA: int = FIXSTAR_OFFSET + 62  # 19 Tauri - Pleiades member
ATLAS: int = FIXSTAR_OFFSET + 63  # 27 Tauri - Pleiades member
PLEIONE: int = FIXSTAR_OFFSET + 64  # 28 Tauri - Pleiades member

# Hyades cluster stars (visible members in Taurus)
PRIMA_HYADUM: int = FIXSTAR_OFFSET + 65  # Gamma Tauri - Hyades member
SECUNDA_HYADUM: int = FIXSTAR_OFFSET + 66  # Delta^1 Tauri - Hyades member
THETA_TAURI: int = FIXSTAR_OFFSET + 67  # Theta^2 Tauri - Hyades member
AIN: int = FIXSTAR_OFFSET + 68  # Epsilon Tauri - Hyades member

# Orion constellation - completing the major stars
MEISSA: int = FIXSTAR_OFFSET + 69  # Lambda Orionis - Orion's head

# Ursa Major (Big Dipper) - completing the asterism
PHECDA: int = FIXSTAR_OFFSET + 70  # Gamma Ursae Majoris - bowl star
MEGREZ: int = FIXSTAR_OFFSET + 71  # Delta Ursae Majoris - bowl-handle junction

# Crux (Southern Cross) - completing the constellation
DELTA_CRUCIS: int = (
    FIXSTAR_OFFSET + 72
)  # Delta Crucis - fourth star of Southern Cross

# Centaurus - completing the bright stars of the constellation
MENKENT: int = FIXSTAR_OFFSET + 73  # Theta Centauri
MUHLIFAIN: int = FIXSTAR_OFFSET + 74  # Gamma Centauri
EPSILON_CENTAURI: int = FIXSTAR_OFFSET + 75  # Epsilon Centauri
ETA_CENTAURI: int = FIXSTAR_OFFSET + 76  # Eta Centauri
ZETA_CENTAURI: int = FIXSTAR_OFFSET + 77  # Zeta Centauri

# Scorpius constellation - completing the bright stars
SARGAS: int = FIXSTAR_OFFSET + 78  # Theta Scorpii
DSCHUBBA: int = FIXSTAR_OFFSET + 79  # Delta Scorpii
GRAFFIAS: int = FIXSTAR_OFFSET + 80  # Beta Scorpii
LESATH: int = FIXSTAR_OFFSET + 81  # Upsilon Scorpii

# Leo constellation - completing the bright stars
ZOSMA: int = FIXSTAR_OFFSET + 82  # Delta Leonis

# ======== ZODIACAL CONSTELLATION BRIGHT STARS ========
# Stars from zodiacal constellations used in astrological interpretation

# Aries constellation - the Ram
HAMAL: int = FIXSTAR_OFFSET + 83  # Alpha Arietis - brightest in Aries
SHERATAN: int = FIXSTAR_OFFSET + 84  # Beta Arietis
MESARTHIM: int = FIXSTAR_OFFSET + 85  # Gamma Arietis

# Cancer constellation - the Crab
ACUBENS: int = FIXSTAR_OFFSET + 86  # Alpha Cancri
TARF: int = FIXSTAR_OFFSET + 87  # Beta Cancri - brightest in Cancer
ASELLUS_BOREALIS: int = FIXSTAR_OFFSET + 88  # Gamma Cancri - Northern Donkey
ASELLUS_AUSTRALIS: int = FIXSTAR_OFFSET + 89  # Delta Cancri - Southern Donkey

# Sagittarius constellation - the Archer
KAUS_AUSTRALIS: int = FIXSTAR_OFFSET + 90  # Epsilon Sagittarii - brightest
NUNKI: int = FIXSTAR_OFFSET + 91  # Sigma Sagittarii
KAUS_MEDIA: int = FIXSTAR_OFFSET + 92  # Delta Sagittarii
KAUS_BOREALIS: int = FIXSTAR_OFFSET + 93  # Lambda Sagittarii
ASCELLA: int = FIXSTAR_OFFSET + 94  # Zeta Sagittarii

# Capricornus constellation - the Sea Goat (complementing Deneb Algedi)
ALGEDI: int = FIXSTAR_OFFSET + 95  # Alpha Capricorni - the Goat
DABIH: int = FIXSTAR_OFFSET + 96  # Beta Capricorni
NASHIRA: int = FIXSTAR_OFFSET + 97  # Gamma Capricorni - the Fortunate One

# Aquarius constellation - the Water Bearer
SADALSUUD: int = FIXSTAR_OFFSET + 98  # Beta Aquarii - brightest in Aquarius
SADALMELIK: int = FIXSTAR_OFFSET + 99  # Alpha Aquarii
SKAT: int = FIXSTAR_OFFSET + 100  # Delta Aquarii

# Pisces constellation - the Fishes
ETA_PISCIUM: int = FIXSTAR_OFFSET + 101  # Eta Piscium - brightest in Pisces
ALRESCHA: int = FIXSTAR_OFFSET + 102  # Alpha Piscium - the Knot

# Andromeda constellation
ALPHERATZ: int = FIXSTAR_OFFSET + 103  # Alpha Andromedae

# Pegasus constellation
ALGENIB: int = FIXSTAR_OFFSET + 104  # Gamma Pegasi

# Gemini constellation - additional stars
PROPUS: int = FIXSTAR_OFFSET + 105  # Eta Geminorum
ALHENA: int = FIXSTAR_OFFSET + 106  # Gamma Geminorum
TEJAT: int = FIXSTAR_OFFSET + 107  # Mu Geminorum
WASAT: int = FIXSTAR_OFFSET + 108  # Delta Geminorum

# Canis Major constellation - additional stars
ADHARA: int = FIXSTAR_OFFSET + 109  # Epsilon Canis Majoris
WEZEN: int = FIXSTAR_OFFSET + 110  # Delta Canis Majoris

# Draco constellation
THUBAN: int = FIXSTAR_OFFSET + 111  # Alpha Draconis

# Hercules constellation
RASALGETHI: int = FIXSTAR_OFFSET + 112  # Alpha1 Herculis

# Cygnus constellation
ALBIREO: int = FIXSTAR_OFFSET + 113  # Beta1 Cygni

# Andromeda constellation (additional)
MIRACH: int = FIXSTAR_OFFSET + 114  # Beta Andromedae
ALMACH: int = FIXSTAR_OFFSET + 115  # Gamma1 Andromedae

# Cetus constellation
MENKAR: int = FIXSTAR_OFFSET + 116  # Alpha Ceti

# Astrological Angles (requires observer location)
ANGLE_OFFSET: int = 2000000
ASCENDANT: int = ANGLE_OFFSET + 1  # Rising sign/degree
_MC_ANGLE_ID: int = ANGLE_OFFSET + 2  # Medium Coeli (Midheaven)
DESCENDANT: int = ANGLE_OFFSET + 3  # Setting point (Asc + 180°)
IC: int = ANGLE_OFFSET + 4  # Imum Coeli (MC + 180°)
_VERTEX_ANGLE_ID: int = ANGLE_OFFSET + 5  # Western intersection of prime vertical
ANTIVERTEX: int = ANGLE_OFFSET + 6  # Eastern intersection (Vertex + 180°)

# Arabic Parts (Lots) - Require pre-calculated planetary positions
ARABIC_OFFSET: int = 3000000
PARS_FORTUNAE: int = ARABIC_OFFSET + 1  # Part of Fortune
PARS_SPIRITUS: int = ARABIC_OFFSET + 2  # Part of Spirit
PARS_AMORIS: int = ARABIC_OFFSET + 3  # Part of Eros/Love
PARS_FIDEI: int = ARABIC_OFFSET + 4  # Part of Faith

# =============================================================================
# CALCULATION FLAGS
# =============================================================================
# Ephemeris selection
FLG_JPLEPH: int = (
    1  # Use JPL ephemeris (default: DE440 via Skyfield, range 1550-2650 CE)
)
FLG_SWIEPH: int = 2  # Use reference ephemeris (same as JPLEPH in libephemeris)
FLG_MOSEPH: int = 4  # Semi-analytical ephemeris flag (accepted for API compatibility, always uses JPL)

# Observer location and reference frame
FLG_HELCTR: int = 8  # Heliocentric position
FLG_TRUEPOS: int = 16  # True geometric position (no light time)
FLG_J2000: int = 32  # J2000.0 reference frame
FLG_NONUT: int = 64  # No nutation
FLG_SPEED3: int = 128  # High precision speed (3 calls)
FLG_SPEED: int = 256  # Calculate velocity
FLG_NOGDEFL: int = 512  # No gravitational deflection
FLG_NOABERR: int = 1024  # No aberration
FLG_ASTROMETRIC: int = FLG_NOABERR | FLG_NOGDEFL  # Astrometric position
FLG_EQUATORIAL: int = 2048  # Equatorial coordinates (RA/Dec)
FLG_XYZ: int = 4096  # Cartesian coordinates
FLG_RADIANS: int = 8192  # Return angles in radians
FLG_BARYCTR: int = 16384  # Barycentric position
FLG_TOPOCTR: int = 32768  # Topocentric position (requires set_topo)
FLG_SIDEREAL: int = 65536  # Sidereal positions
FLG_ICRS: int = 131072  # ICRS reference frame

# =============================================================================
# REFERENCE API-COMPATIBLE FLAG ALIASES (FLG_* instead of SEFLG_*)
# =============================================================================
# FLG_* prefix aliases for SEFLG_* flags
# These aliases provide full API compatibility

# SIDEREAL (AYANAMSHA) MODES
# =============================================================================

# Western sidereal traditions
SIDM_FAGAN_BRADLEY: int = 0  # Fagan-Bradley (Synetic Vernal Point)
SIDM_LAHIRI: int = 1  # Lahiri (Chitrapaksha, Indian standard)
SIDM_DELUCE: int = 2  # De Luce
SIDM_RAMAN: int = 3  # B.V. Raman
SIDM_USHASHASHI: int = 4  # Ushashashi
SIDM_KRISHNAMURTI: int = 5  # K.S. Krishnamurti (KP)
SIDM_DJWHAL_KHUL: int = 6  # Djwhal Khul (Alice Bailey)
SIDM_YUKTESHWAR: int = 7  # Yukteshwar
SIDM_JN_BHASIN: int = 8  # J.N. Bhasin

# Babylonian traditions
SIDM_BABYL_KUGLER1: int = 9  # Kugler variant 1
SIDM_BABYL_KUGLER2: int = 10  # Kugler variant 2
SIDM_BABYL_KUGLER3: int = 11  # Kugler variant 3
SIDM_BABYL_HUBER: int = 12  # Huber
SIDM_BABYL_ETPSC: int = 13  # ETPSC
SIDM_BABYL_BRITTON: int = 38  # Britton

# Star-based ayanamshas
SIDM_ALDEBARAN_15TAU: int = 14  # Aldebaran at 15° Taurus
SIDM_TRUE_CITRA: int = 27  # True position of Spica (180° Citra)
SIDM_TRUE_REVATI: int = 28  # True position of Revati
SIDM_TRUE_PUSHYA: int = 29  # True position of Pushya (Cancri)
SIDM_TRUE_MULA: int = 35  # True position of Mula (λ Scorpii)
SIDM_TRUE_SHEORAN: int = 39  # True Sheoran

# Historical epochs
SIDM_HIPPARCHOS: int = 15  # Hipparchos (128 BC)
SIDM_SASSANIAN: int = 16  # Sassanian
SIDM_J2000: int = 18  # J2000.0 (no ayanamsha)
SIDM_J1900: int = 19  # J1900.0
SIDM_B1950: int = 20  # B1950.0

# Indian traditions
SIDM_SURYASIDDHANTA: int = 21  # Suryasiddhanta
SIDM_SURYASIDDHANTA_MSUN: int = 22  # Suryasiddhanta (mean Sun)
SIDM_ARYABHATA: int = 23  # Aryabhata
SIDM_ARYABHATA_MSUN: int = 24  # Aryabhata (mean Sun)
SIDM_ARYABHATA_522: int = 37  # Aryabhata 522
SIDM_SS_REVATI: int = 25  # Suryasiddhanta Revati
SIDM_SS_CITRA: int = 26  # Suryasiddhanta Citra

# Galactic alignment systems
SIDM_GALCENT_0SAG: int = 17  # Galactic Center at 0° Sagittarius
SIDM_GALCENT_RGILBRAND: int = 30  # Galactic Center (Gil Brand)
SIDM_GALCENT_MULA_WILHELM: int = 36  # Galactic Center at Mula (Wilhelm)
SIDM_GALCENT_COCHRANE: int = 40  # Galactic Center (Cochrane)
SIDM_GALEQU_IAU1958: int = 31  # Galactic Equator (IAU 1958)
SIDM_GALEQU_TRUE: int = 32  # Galactic Equator (True)
SIDM_GALEQU_MULA: int = 33  # Galactic Equator at Mula
SIDM_GALEQU_FIORENZA: int = 41  # Galactic Equator (Fiorenza)
SIDM_GALALIGN_MARDYKS: int = 34  # Galactic Alignment (Mardyks)

# Other systems
SIDM_VALENS_MOON: int = 42  # Vettius Valens (Moon-based)
SIDM_LAHIRI_1940: int = 43  # Lahiri (1940 Lahiri Commission value)
SIDM_LAHIRI_VP285: int = 44  # Lahiri (Vernal Point 285 CE)
SIDM_KRISHNAMURTI_VP291: int = 45  # Krishnamurti (Vernal Point 291 CE)
SIDM_LAHIRI_ICRC: int = 46  # Lahiri (ICRC, Indian Calendar Reform Committee)
SIDM_USER: int = 255  # User-defined ayanamsha

# =============================================================================
# REFERENCE API-COMPATIBLE SIDEREAL MODE ALIASES (SIDM_* instead of SIDM_*)
# =============================================================================
# SIDM_* prefix aliases for SIDM_* constants

# Western sidereal traditions

# CALENDAR SYSTEMS
# =============================================================================

JUL_CAL: int = 0  # Julian calendar
GREG_CAL: int = 1  # Gregorian calendar

# reference API-compatible aliases (without SE_ prefix)

# ASTRONOMICAL CONSTANTS
# =============================================================================

AUNIT: float = 149597870.7  # Astronomical Unit in km (IAU 2012 standard)

_MOON_MEAN_DIST_KM: float = 384400.0  # Mean Earth-Moon distance in km
_MOON_MEAN_ECC: float = 0.054900489  # Mean lunar orbital eccentricity
_MOON_MEAN_DIST_AU: float = _MOON_MEAN_DIST_KM / AUNIT  # ~0.002569555 AU
_MOON_MEAN_APOG_DIST_AU: float = (
    _MOON_MEAN_DIST_KM * (1.0 + _MOON_MEAN_ECC) / AUNIT
)  # ~0.002710625 AU

# =============================================================================
# ECLIPSE TYPES AND FLAGS
# =============================================================================
# Eclipse geometry types
ECL_CENTRAL: int = 1  # Central eclipse (Moon's shadow axis crosses Earth)
ECL_NONCENTRAL: int = 2  # Non-central eclipse
ECL_TOTAL: int = 4  # Total eclipse (Sun/Earth completely covered)
ECL_ANNULAR: int = 8  # Annular eclipse (ring of fire)
ECL_PARTIAL: int = 16  # Partial eclipse
ECL_ANNULAR_TOTAL: int = (
    32  # Hybrid eclipse (annular at some locations, total at others)
)
ECL_PENUMBRAL: int = 64  # Penumbral lunar eclipse
ECL_GRAZING: int = 65536  # Grazing occultation (target passes near lunar limb)

# Composite eclipse type masks
ECL_ALLTYPES_SOLAR: int = (
    ECL_CENTRAL
    | ECL_NONCENTRAL
    | ECL_TOTAL
    | ECL_ANNULAR
    | ECL_PARTIAL
    | ECL_ANNULAR_TOTAL
)
ECL_ALLTYPES_LUNAR: int = ECL_TOTAL | ECL_PARTIAL | ECL_PENUMBRAL

# Eclipse visibility and contact flags
ECL_VISIBLE: int = 128  # Eclipse visible at location
ECL_MAX_VISIBLE: int = 256  # Maximum phase visible
ECL_1ST_VISIBLE: int = 512  # First contact visible
ECL_2ND_VISIBLE: int = 1024  # Second contact visible
ECL_3RD_VISIBLE: int = 2048  # Third contact visible
ECL_4TH_VISIBLE: int = 4096  # Fourth contact visible
ECL_ONE_TRY: int = 32768  # Try only once (optimization flag)

# reference API-compatible aliases (without SE_ prefix)

# NODAL/APSIDAL CALCULATION METHOD FLAGS
# =============================================================================
# Method flags for nod_aps() and nod_aps_ut() functions

NODBIT_MEAN: int = 1  # Mean orbital elements (averaged over perturbations)
NODBIT_OSCU: int = 2  # Osculating elements (instantaneous/perturbed)
NODBIT_OSCU_BAR: int = 4  # Barycentric osculating elements
NODBIT_FOPOINT: int = 256  # Focal point (second focus of ellipse)

# reference API-compatible aliases (without SE_ prefix)

# RISE/SET/TRANSIT CALCULATION FLAGS
# =============================================================================
# Event type flags for rise_trans() rsmi parameter

CALC_RISE: int = 1  # Calculate rise time
CALC_SET: int = 2  # Calculate set time
CALC_MTRANSIT: int = 4  # Calculate meridian transit (upper culmination)
CALC_ITRANSIT: int = 8  # Calculate lower transit (anti-culmination)

# Bitmask flags for additional rise/set options
BIT_DISC_CENTER: int = 256  # Use center of disc instead of upper limb
BIT_DISC_BOTTOM: int = 8192  # Use lower limb of disc
BIT_NO_REFRACTION: int = 512  # No atmospheric refraction
BIT_CIVIL_TWILIGHT: int = 1024  # Civil twilight (Sun at -6 degrees)
BIT_NAUTIC_TWILIGHT: int = 2048  # Nautical twilight (Sun at -12 degrees)
BIT_ASTRO_TWILIGHT: int = 4096  # Astronomical twilight (Sun at -18 degrees)
BIT_FIXED_DISC_SIZE: int = 16384  # Use fixed disc size (ignore parallax)

# reference API-compatible aliases (without SE_ prefix)

# HELIACAL EVENT TYPES
# =============================================================================
# Event types for heliacal_ut() function

HELIACAL_RISING: int = 1  # Heliacal rising (morning first)
HELIACAL_SETTING: int = 2  # Heliacal setting (evening last)
MORNING_FIRST: int = HELIACAL_RISING  # Alias: first visibility at morning
EVENING_LAST: int = HELIACAL_SETTING  # Alias: last visibility at evening
EVENING_FIRST: int = 3  # First visibility at evening (after superior conjunction)
MORNING_LAST: int = 4  # Last visibility at morning (before superior conjunction)

# reference API-compatible aliases (without SE_ prefix)

# HELIACAL VISIBILITY FLAGS
# =============================================================================
# Flags for vis_limit_mag() and heliacal calculations

HELFLAG_LONG_SEARCH: int = 128  # Long search (extended search period)
HELFLAG_HIGH_PRECISION: int = 256  # High precision mode
HELFLAG_OPTICAL_PARAMS: int = 1 << 9  # 512 - Use optical instrument parameters
HELFLAG_NO_DETAILS: int = 1 << 10  # 1024 - Skip detailed calculations
HELFLAG_SEARCH_1_PERIOD: int = 1 << 11  # 2048 - Search one synodic period only
HELFLAG_VISLIM_DARK: int = 1 << 12  # 4096 - Assume Sun at nadir (dark sky)
HELFLAG_VISLIM_NOMOON: int = 1 << 13  # 8192 - Exclude Moon's contribution
HELFLAG_VISLIM_PHOTOPIC: int = 1 << 14  # 16384 - Force photopic vision mode
HELFLAG_AV: int = 1 << 16  # 65536 - Arcus visionis method (VR)
HELFLAG_AVKIND_VR: int = 1 << 16  # 65536 - Alias for HELFLAG_AV
HELFLAG_AVKIND_PTO: int = 1 << 17  # 131072 - Ptolemy method
HELFLAG_AVKIND_MIN7: int = 1 << 18  # 262144 - Babylonian method (min 7°)
HELFLAG_AVKIND_MIN9: int = 1 << 19  # 524288 - Babylonian method (min 9°)
HELFLAG_AVKIND: int = (  # 983040 - Mask for all AV kind bits
    HELFLAG_AVKIND_VR
    | HELFLAG_AVKIND_PTO
    | HELFLAG_AVKIND_MIN7
    | HELFLAG_AVKIND_MIN9
)

# reference API-compatible aliases (without SE_ prefix)

HELFLAG_BELOW_HORIZON: int = -2  # Object is below horizon
HELFLAG_PHOTOPIC: int = 0  # OK, photopic (daylight) vision
HELFLAG_SCOTOPIC: int = 1  # OK, scotopic (night) vision
HELFLAG_MIXED: int = 2  # OK, near limit photopic/scotopic vision

# SPLIT_DEG FLAGS
# =============================================================================
# Flags for split_deg() function to control output format and rounding

SPLIT_DEG_ROUND_SEC: int = 1  # Round to seconds
SPLIT_DEG_ROUND_MIN: int = 2  # Round to minutes
SPLIT_DEG_ROUND_DEG: int = 4  # Round to degrees
SPLIT_DEG_ZODIACAL: int = 8  # Return zodiac sign number (0-11)
SPLIT_DEG_NAKSHATRA: int = 1024  # Return nakshatra number (0-26)
SPLIT_DEG_KEEP_SIGN: int = 16  # Don't round to next zodiac sign/nakshatra
SPLIT_DEG_KEEP_DEG: int = 32  # Don't round to next degree

# =============================================================================
# TIDAL ACCELERATION CONSTANTS
# =============================================================================
# Tidal acceleration of the Moon in arcsec/century^2, for Delta T calculations.
# Different JPL ephemeris files use different tidal acceleration values.
# These affect the polynomial extrapolation of Delta T for historical dates.

TIDAL_DE200: float = -23.8946  # DE200 (older ephemeris)
TIDAL_DE403: float = -25.580  # DE403
TIDAL_DE404: float = -25.580  # DE404
TIDAL_DE405: float = -25.826  # DE405
TIDAL_DE406: float = -25.826  # DE406
TIDAL_DE421: float = -25.85  # DE421 (legacy default)
TIDAL_DE422: float = -25.85  # DE422
TIDAL_DE430: float = -25.82  # DE430
TIDAL_DE431: float = -25.80  # DE431
TIDAL_DE440: float = -25.936  # DE440 (current default)
TIDAL_DE441: float = -25.936  # DE441 (latest, same as DE440)
TIDAL_DEFAULT: float = -25.8  # Default value (matches pyswisseph)
TIDAL_AUTOMATIC: int = 999999  # Let library choose based on ephemeris file

# reference API-compatible aliases (without SE_ prefix)

# STANDARD ASTRONOMICAL EPOCHS (JULIAN DAY)
# =============================================================================
# Reference epochs used for astrometric calculations and proper motion.
# All values are in Julian Day TT (Terrestrial Time).

# J2000.0 epoch: Jan 1, 2000, 12:00 TT - Standard reference epoch for ICRS
J2000: float = 2451545.0

# J1991.25 epoch: Apr 2, 1991, 13:30 TT - Hipparcos catalog reference epoch
# This is the mean epoch of observations for the Hipparcos mission.
# Proper motion values in Hipparcos are defined relative to this epoch.
# Reference: Hipparcos and Tycho Catalogues, ESA SP-1200, Vol. 1, Section 1.2.2
J1991_25: float = 2448349.0625

# J1900.0 epoch: Jan 0.5, 1900 TT - Historical reference epoch
J1900: float = 2415020.0

# B1950.0 epoch: Jan 0.923, 1950 - Besselian epoch for FK4 catalog
B1950: float = 2433282.4235

# Days per Julian year (exactly 365.25 days)
DAYS_PER_JULIAN_YEAR: float = 365.25

# Days per Julian century (exactly 36525 days)
DAYS_PER_JULIAN_CENTURY: float = 36525.0

# =============================================================================
# PLANETARY MOON IDENTIFIERS
# =============================================================================
# Canonical body IDs for planetary satellites follow the reference API 2.10+
# convention: ipl = PLMOON_OFFSET + NAIF satellite id, e.g. Io = 9501,
# Titan = 9606, Charon = 9901.  Any NAIF satellite id (401-998) is accepted
# this way; 9n99 addresses the planet's center of body (e.g. 9599 Jupiter).
# Satellite SPK files (jup365.bsp, sat441.bsp, ...) must be registered with
# register_moon_spk() before calculation.
#
# The MOON_* constants below use a legacy private numbering
# (MOON_OFFSET + small index, e.g. Io = 9001) that predates the canonical
# scheme.  They are kept as deprecated aliases; new code should use
# PLMOON_OFFSET + NAIF id.

MOON_OFFSET: int = 9000  # Base offset for legacy moon IDs (deprecated)

# Jupiter's Galilean Moons (discovered by Galileo Galilei in 1610)
MOON_IO: int = MOON_OFFSET + 1  # Jupiter I - innermost Galilean moon
MOON_EUROPA: int = MOON_OFFSET + 2  # Jupiter II - potential subsurface ocean
MOON_GANYMEDE: int = MOON_OFFSET + 3  # Jupiter III - largest moon in solar system
MOON_CALLISTO: int = MOON_OFFSET + 4  # Jupiter IV - heavily cratered

# Saturn's Major Moons
MOON_MIMAS: int = MOON_OFFSET + 11  # Saturn I - "Death Star" appearance
MOON_ENCELADUS: int = MOON_OFFSET + 12  # Saturn II - active geysers
MOON_TETHYS: int = MOON_OFFSET + 13  # Saturn III - icy body
MOON_DIONE: int = MOON_OFFSET + 14  # Saturn IV - trailing hemisphere features
MOON_RHEA: int = MOON_OFFSET + 15  # Saturn V - second largest Saturn moon
MOON_TITAN: int = MOON_OFFSET + 16  # Saturn VI - thick atmosphere, liquid lakes
MOON_HYPERION: int = MOON_OFFSET + 17  # Saturn VII - chaotic rotation
MOON_IAPETUS: int = MOON_OFFSET + 18  # Saturn VIII - two-toned surface

# Uranus' Major Moons
MOON_MIRANDA: int = MOON_OFFSET + 21  # Uranus V - extreme geological features
MOON_ARIEL: int = MOON_OFFSET + 22  # Uranus I - brightest Uranian moon
MOON_UMBRIEL: int = MOON_OFFSET + 23  # Uranus II - darkest Uranian moon
MOON_TITANIA: int = MOON_OFFSET + 24  # Uranus III - largest Uranian moon
MOON_OBERON: int = MOON_OFFSET + 25  # Uranus IV - outermost major Uranian moon

# Neptune's Major Moon
MOON_TRITON: int = MOON_OFFSET + 31  # Neptune I - retrograde orbit, captured KBO

# Mars' Moons
MOON_PHOBOS: int = MOON_OFFSET + 41  # Mars I - larger, closer moon
MOON_DEIMOS: int = MOON_OFFSET + 42  # Mars II - smaller, farther moon

# Pluto's Moon
MOON_CHARON: int = MOON_OFFSET + 51  # Pluto I - largest moon (binary system)

# Aliases without SE_ prefix for reference API compatibility

# NAIF IDS FOR PLANETARY MOONS
# =============================================================================
# Standard NAIF SPICE IDs for planetary satellites

# Jupiter system (5xx)
NAIF_IO: int = 501
NAIF_EUROPA: int = 502
NAIF_GANYMEDE: int = 503
NAIF_CALLISTO: int = 504

# Saturn system (6xx)
NAIF_MIMAS: int = 601
NAIF_ENCELADUS: int = 602
NAIF_TETHYS: int = 603
NAIF_DIONE: int = 604
NAIF_RHEA: int = 605
NAIF_TITAN: int = 606
NAIF_HYPERION: int = 607
NAIF_IAPETUS: int = 608

# Uranus system (7xx)
NAIF_MIRANDA: int = 705
NAIF_ARIEL: int = 701
NAIF_UMBRIEL: int = 702
NAIF_TITANIA: int = 703
NAIF_OBERON: int = 704

# Neptune system (8xx)
NAIF_TRITON: int = 801

# Mars system (4xx)
NAIF_PHOBOS: int = 401
NAIF_DEIMOS: int = 402

# Pluto system (9xx)
NAIF_CHARON: int = 901

# =============================================================================
# BARE ALIASES FOR SE_* CONSTANTS (pyswisseph compatibility)
# =============================================================================
# These provide the same constants without the SE_ prefix, matching the names
# that pyswisseph exposes as module-level attributes (e.g. swe.ECL_NUT,
# swe.MEAN_NODE, swe.CHIRON, etc.).

# Special values

# ADDITIONAL PYSWISSEPH-COMPATIBLE CONSTANTS
# =============================================================================
# These constants match pyswisseph module-level attributes that were not
# previously exported. Added for full API compatibility.

# House cusps and special points
ASC: int = 0
MC: int = 1
ARMC: int = 2
VERTEX: int = 3
EQUASC: int = 4
COASC1: int = 5
COASC2: int = 6
POLASC: int = 7
NASCMC: int = 8

# Additional body IDs
FIXSTAR: int = -10
HARRINGTON: int = 50
NEPTUNE_ADAMS: int = 52
NEPTUNE_LEVERRIER: int = 51
NIBIRU: int = 49
PLUTO_LOWELL: int = 53
PLUTO_PICKERING: int = 54

# Counts and offsets
FICT_MAX: int = 999
FICT_OFFSET_1: int = 39
MAX_STNAME: int = 256
NALL_NAT_POINTS: int = 38
PLMOON_OFFSET: int = 9000

# Canonical planetary-moon ids (PLMOON_OFFSET + NAIF satellite id).
# These match the ipl numbers the reference API uses for its sepm* files.
PLMOON_IO: int = PLMOON_OFFSET + NAIF_IO  # 9501
PLMOON_EUROPA: int = PLMOON_OFFSET + NAIF_EUROPA  # 9502
PLMOON_GANYMEDE: int = PLMOON_OFFSET + NAIF_GANYMEDE  # 9503
PLMOON_CALLISTO: int = PLMOON_OFFSET + NAIF_CALLISTO  # 9504
PLMOON_MIMAS: int = PLMOON_OFFSET + NAIF_MIMAS  # 9601
PLMOON_ENCELADUS: int = PLMOON_OFFSET + NAIF_ENCELADUS  # 9602
PLMOON_TETHYS: int = PLMOON_OFFSET + NAIF_TETHYS  # 9603
PLMOON_DIONE: int = PLMOON_OFFSET + NAIF_DIONE  # 9604
PLMOON_RHEA: int = PLMOON_OFFSET + NAIF_RHEA  # 9605
PLMOON_TITAN: int = PLMOON_OFFSET + NAIF_TITAN  # 9606
PLMOON_HYPERION: int = PLMOON_OFFSET + NAIF_HYPERION  # 9607
PLMOON_IAPETUS: int = PLMOON_OFFSET + NAIF_IAPETUS  # 9608
PLMOON_MIRANDA: int = PLMOON_OFFSET + NAIF_MIRANDA  # 9705
PLMOON_ARIEL: int = PLMOON_OFFSET + NAIF_ARIEL  # 9701
PLMOON_UMBRIEL: int = PLMOON_OFFSET + NAIF_UMBRIEL  # 9702
PLMOON_TITANIA: int = PLMOON_OFFSET + NAIF_TITANIA  # 9703
PLMOON_OBERON: int = PLMOON_OFFSET + NAIF_OBERON  # 9704
PLMOON_TRITON: int = PLMOON_OFFSET + NAIF_TRITON  # 9801
PLMOON_PHOBOS: int = PLMOON_OFFSET + NAIF_PHOBOS  # 9401
PLMOON_DEIMOS: int = PLMOON_OFFSET + NAIF_DEIMOS  # 9402
PLMOON_CHARON: int = PLMOON_OFFSET + NAIF_CHARON  # 9901

# Additional calculation flags
FLG_CENTER_BODY: int = 1048576
FLG_DEFAULTEPH: int = 2
FLG_DPSIDEPS_1980: int = 262144
FLG_JPLHOR: int = 262144
FLG_JPLHOR_APPROX: int = 524288
FLG_ORBEL_AA: int = 32768
FLG_TEST_PLMOON: int = 2228280
FLG_TROPICAL: int = 0

# Ephemeris file names
FNAME_DE200: str = "de200.eph"
FNAME_DE403: str = "de403.eph"
FNAME_DE404: str = "de404.eph"
FNAME_DE405: str = "de405.eph"
FNAME_DE406: str = "de406.eph"
FNAME_DFT: str = "de431.eph"
FNAME_DFT2: str = "de406.eph"

# Additional eclipse constants
ECL_HYBRID: int = 32
ECL_OCC_BEG_DAYLIGHT: int = 8192
ECL_OCC_END_DAYLIGHT: int = 16384
ECL_PARTBEG_VISIBLE: int = 512
ECL_PARTEND_VISIBLE: int = 4096
ECL_PENUMBBEG_VISIBLE: int = 8192
ECL_PENUMBEND_VISIBLE: int = 16384
ECL_TOTBEG_VISIBLE: int = 1024
ECL_TOTEND_VISIBLE: int = 2048

# Sidereal mode constants
NSIDM_PREDEF: int = 47
SIDBITS: int = 256
SIDBIT_ECL_DATE: int = 2048
SIDBIT_ECL_T0: int = 256
SIDBIT_NO_PREC_OFFSET: int = 4096
SIDBIT_PREC_ORIG: int = 8192
SIDBIT_SSY_PLANE: int = 512
SIDBIT_USER_UT: int = 1024

# Precession, nutation, and delta-T model constants
MODEL_BIAS: int = 4
MODEL_DELTAT: int = 0
MODEL_JPLHORA_MODE: int = 6
MODEL_JPLHOR_MODE: int = 5
MODEL_NUT: int = 3
MODEL_PREC_LONGTERM: int = 1
MODEL_PREC_SHORTTERM: int = 2
MODEL_SIDT: int = 7
NSE_MODELS: int = 8

MOD_BIAS_DEFAULT: int = 3
MOD_BIAS_IAU2000: int = 2
MOD_BIAS_IAU2006: int = 3
MOD_BIAS_NONE: int = 1
MOD_NBIAS: int = 3

MOD_DELTAT_DEFAULT: int = 5
MOD_DELTAT_ESPENAK_MEEUS_2006: int = 4
MOD_DELTAT_STEPHENSON_1997: int = 2
MOD_DELTAT_STEPHENSON_ETC_2016: int = 5
MOD_DELTAT_STEPHENSON_MORRISON_1984: int = 1
MOD_DELTAT_STEPHENSON_MORRISON_2004: int = 3
MOD_NDELTAT: int = 5

MOD_JPLHORA_1: int = 1
MOD_JPLHORA_2: int = 2
MOD_JPLHORA_3: int = 3
MOD_JPLHORA_DEFAULT: int = 3
MOD_NJPLHORA: int = 3

MOD_JPLHOR_DEFAULT: int = 1
MOD_JPLHOR_LONG_AGREEMENT: int = 1
MOD_NJPLHOR: int = 2

MOD_NUT_DEFAULT: int = 4
MOD_NUT_IAU_1980: int = 1
MOD_NUT_IAU_2000A: int = 3
MOD_NUT_IAU_2000B: int = 4
MOD_NUT_IAU_CORR_1987: int = 2
MOD_NUT_WOOLARD: int = 5
MOD_NNUT: int = 5

MOD_PREC_BRETAGNON_2003: int = 7
MOD_PREC_DEFAULT: int = 9
MOD_PREC_DEFAULT_SHORT: int = 9
MOD_PREC_IAU_1976: int = 1
MOD_PREC_IAU_2000: int = 6
MOD_PREC_IAU_2006: int = 8
MOD_PREC_LASKAR_1986: int = 2
MOD_PREC_NEWCOMB: int = 11
MOD_PREC_OWEN_1990: int = 10
MOD_PREC_SIMON_1994: int = 5
MOD_PREC_VONDRAK_2011: int = 9
MOD_PREC_WILLIAMS_1994: int = 4
MOD_PREC_WILL_EPS_LASK: int = 3
MOD_NPREC: int = 11

# Heliacal/rise-set bit flags
BIT_FORCE_SLOW_METHOD: int = 32768
BIT_HINDU_RISING: int = 896

# Rise/set type constants
ACRONYCHAL_RISING: int = 5
ACRONYCHAL_SETTING: int = 6
COSMICAL_SETTING: int = 6

# Tidal acceleration constants
TIDAL_26: float = -26.0
TIDAL_JPLEPH: float = -25.8
TIDAL_MOSEPH: float = -25.58
TIDAL_STEPHENSON_2016: float = -25.85
TIDAL_SWIEPH: float = -25.8

# File names and paths
ASTNAMFILE: str = "seasnam.txt"
EPHE_PATH: str = ".:/users/ephe2/:/users/ephe/"
FICTFILE: str = "seorbel.txt"
FNAME_DE431: str = "de431.eph"
# IMPORTANT — DO NOT REMOVE during another cleanup pass.
# The upstream reference distribution exports this single value under
# the SE_-prefixed name, so we mirror it for 1:1 parity. This is the
# one intentional exception to the "no SE_/SEFLG_/swe_ prefix" rule
# and is explicitly listed in ``ALLOWED_PREFIXED_NAMES`` of
# tests/test_api_compat/test_api_surface.py.
SE_FNAME_DE431: str = FNAME_DE431
STARFILE: str = "sefstars.txt"
STARFILE_OLD: str = "fixstars.cat"

# Unit conversion constants
AUNIT_TO_KM: float = 149597870.7
AUNIT_TO_LIGHTYEAR: float = 1.5812507409819728e-05
AUNIT_TO_PARSEC: float = 4.848136811095274e-06

# Miscellaneous
DE_NUMBER: int = 431
DELTAT_AUTOMATIC: float = -1e-10
MIXEDOPIC_FLAG: int = 2
PHOTOPIC_FLAG: int = 0
SCOTOPIC_FLAG: int = 1
SIMULATE_VICTORVB: int = 1
TJD_INVALID: float = 99999999.0

# =============================================================================
# __all__ — explicit public API (excludes 'annotations' and private names)
# =============================================================================
__all__ = [
    _name
    for _name in sorted(dir())
    if not _name.startswith("_") and _name != "annotations"
]
