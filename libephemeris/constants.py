"""
Swiss Ephemeris API-compatible constants for libephemeris.

This module defines all constants used for planetary calculations, including:
- Planet/Body IDs: Numeric identifiers for celestial bodies
- Calculation Flags: Bitwise flags controlling observation parameters  
- Sidereal Modes: Ayanamsha systems for sidereal astrology
- Calendar Systems: Julian vs Gregorian calendar selection
- Eclipse Types: Classification of solar and lunar eclipses

Constants are organized into logical groups for easy navigation.
Values match Swiss Ephemeris v2.x for API compatibility.
"""

# =============================================================================
# PLANET AND BODY IDENTIFIERS
# =============================================================================

# Special values
SE_ECL_NUT: int = -1  # Nutation and obliquity calculation

# Major planets (traditional + modern)
SE_SUN: int = 0
SE_MOON: int = 1
SE_MERCURY: int = 2
SE_VENUS: int = 3
SE_MARS: int = 4
SE_JUPITER: int = 5
SE_SATURN: int = 6
SE_URANUS: int = 7
SE_NEPTUNE: int = 8
SE_PLUTO: int = 9

# Lunar nodes and apsides
SE_MEAN_NODE: int = 10  # Mean lunar node (Dragon's Head)
SE_TRUE_NODE: int = 11  # True (osculating) lunar node
SE_MEAN_APOG: int = 12  # Mean lunar apogee (Black Moon Lilith)
SE_OSCU_APOG: int = 13  # Osculating (true) lunar apogee

# Earth and centaurs
SE_EARTH: int = 14
SE_CHIRON: int = 15  # Centaur between Saturn and Uranus
SE_PHOLUS: int = 16  # Centaur beyond Saturn

# Main belt asteroids
SE_CERES: int = 17
SE_PALLAS: int = 18
SE_JUNO: int = 19
SE_VESTA: int = 20

# Interpolated lunar apsides
SE_INTP_APOG: int = 21  # Interpolated apogee
SE_INTP_PERG: int = 22  # Interpolated perigee

# Count and offsets
SE_NPLANETS: int = 23  # Total number of standard planet IDs
SE_AST_OFFSET: int = 10000  # Offset for asteroid catalog numbers
SE_VARUNA: int = SE_AST_OFFSET + 20000  # TNO Varuna
SE_FICT_OFFSET: int = 40  # Offset for fictitious bodies
SE_NFICT_ELEM: int = 15  # Number of fictitious elements
SE_COMET_OFFSET: int = 1000  # Offset for comet IDs
SE_NALL_NAT_POINTS: int = (
    SE_NPLANETS + SE_NFICT_ELEM + SE_AST_OFFSET + SE_COMET_OFFSET
)

# Trans-Neptunian Objects (TNOs) - Catalog number + offset
SE_ERIS: int = 136199 + SE_AST_OFFSET  # Largest known dwarf planet
SE_SEDNA: int = 90377 + SE_AST_OFFSET  # Detached TNO
SE_HAUMEA: int = 136108 + SE_AST_OFFSET  # Fast-rotating dwarf planet
SE_MAKEMAKE: int = 136472 + SE_AST_OFFSET  # Classical Kuiper belt object
SE_IXION: int = 28978 + SE_AST_OFFSET  # Plutino
SE_ORCUS: int = 90482 + SE_AST_OFFSET  # Plutino, "anti-Pluto"
SE_QUAOAR: int = 50000 + SE_AST_OFFSET  # Classical KBO

# =============================================================================
# VIRTUAL POINTS AND CALCULATED POSITIONS
# =============================================================================

# Fixed Stars (high offset to avoid ID collisions)
SE_FIXSTAR_OFFSET: int = 1000000
SE_REGULUS: int = SE_FIXSTAR_OFFSET + 1  # Alpha Leonis
SE_SPICA_STAR: int = SE_FIXSTAR_OFFSET + 2  # Alpha Virginis

# Astrological Angles (requires observer location)
SE_ANGLE_OFFSET: int = 2000000
SE_ASCENDANT: int = SE_ANGLE_OFFSET + 1  # Rising sign/degree
SE_MC: int = SE_ANGLE_OFFSET + 2  # Medium Coeli (Midheaven)
SE_DESCENDANT: int = SE_ANGLE_OFFSET + 3  # Setting point (Asc + 180°)
SE_IC: int = SE_ANGLE_OFFSET + 4  # Imum Coeli (MC + 180°)
SE_VERTEX: int = SE_ANGLE_OFFSET + 5  # Western intersection of prime vertical
SE_ANTIVERTEX: int = SE_ANGLE_OFFSET + 6  # Eastern intersection (Vertex + 180°)

# Arabic Parts (Lots) - Require pre-calculated planetary positions
SE_ARABIC_OFFSET: int = 3000000
SE_PARS_FORTUNAE: int = SE_ARABIC_OFFSET + 1  # Part of Fortune
SE_PARS_SPIRITUS: int = SE_ARABIC_OFFSET + 2  # Part of Spirit
SE_PARS_AMORIS: int = SE_ARABIC_OFFSET + 3  # Part of Eros/Love
SE_PARS_FIDEI: int = SE_ARABIC_OFFSET + 4  # Part of Faith

# =============================================================================
# CALCULATION FLAGS
# =============================================================================
# Ephemeris selection (currently only SWIEPH/JPL mode supported)
SEFLG_JPLEPH: int = 1  # Use JPL ephemeris
SEFLG_SWIEPH: int = 2  # Use Swiss Ephemeris (libephemeris uses Skyfield/JPL)
SEFLG_MOSEPH: int = 4  # Use Moshier ephemeris (not supported)

# Observer location and reference frame
SEFLG_HELCTR: int = 8  # Heliocentric position
SEFLG_TRUEPOS: int = 16  # True geometric position (no light time)
SEFLG_J2000: int = 32  # J2000.0 reference frame
SEFLG_NONUT: int = 64  # No nutation
SEFLG_SPEED3: int = 128  # High precision speed (3 calls)
SEFLG_SPEED: int = 256  # Calculate velocity
SEFLG_NOGDEFL: int = 512  # No gravitational deflection
SEFLG_NOABERR: int = 1024  # No aberration
SEFLG_ASTROMETRIC: int = SEFLG_NOABERR | SEFLG_NOGDEFL  # Astrometric position
SEFLG_EQUATORIAL: int = 2048  # Equatorial coordinates (RA/Dec)
SEFLG_XYZ: int = 4096  # Cartesian coordinates
SEFLG_RADIANS: int = 8192  # Return angles in radians
SEFLG_BARYCTR: int = 16384  # Barycentric position
SEFLG_TOPOCTR: int = 32768  # Topocentric position (requires swe_set_topo)
SEFLG_SIDEREAL: int = 65536  # Sidereal positions
SEFLG_ICRS: int = 131072  # ICRS reference frame

# =============================================================================
# PYSWISSEPH-COMPATIBLE FLAG ALIASES (FLG_* instead of SEFLG_*)
# =============================================================================
# pyswisseph uses FLG_* prefix while Swiss Ephemeris C library uses SEFLG_*
# These aliases provide full API compatibility with pyswisseph

FLG_JPLEPH: int = SEFLG_JPLEPH
FLG_SWIEPH: int = SEFLG_SWIEPH
FLG_MOSEPH: int = SEFLG_MOSEPH
FLG_HELCTR: int = SEFLG_HELCTR
FLG_TRUEPOS: int = SEFLG_TRUEPOS
FLG_J2000: int = SEFLG_J2000
FLG_NONUT: int = SEFLG_NONUT
FLG_SPEED3: int = SEFLG_SPEED3
FLG_SPEED: int = SEFLG_SPEED
FLG_NOGDEFL: int = SEFLG_NOGDEFL
FLG_NOABERR: int = SEFLG_NOABERR
FLG_ASTROMETRIC: int = SEFLG_ASTROMETRIC
FLG_EQUATORIAL: int = SEFLG_EQUATORIAL  # Equatorial coordinates (RA/Dec)
FLG_XYZ: int = SEFLG_XYZ
FLG_RADIANS: int = SEFLG_RADIANS
FLG_BARYCTR: int = SEFLG_BARYCTR
FLG_TOPOCTR: int = SEFLG_TOPOCTR
FLG_SIDEREAL: int = SEFLG_SIDEREAL
FLG_ICRS: int = SEFLG_ICRS

# Other aliases
AST_OFFSET: int = SE_AST_OFFSET

# =============================================================================
# SIDEREAL (AYANAMSHA) MODES
# =============================================================================

# Western sidereal traditions
SE_SIDM_FAGAN_BRADLEY: int = 0  # Fagan-Bradley (Synetic Vernal Point)
SE_SIDM_LAHIRI: int = 1  # Lahiri (Chitrapaksha, Indian standard)
SE_SIDM_DELUCE: int = 2  # De Luce
SE_SIDM_RAMAN: int = 3  # B.V. Raman
SE_SIDM_USHASHASHI: int = 4  # Ushashashi
SE_SIDM_KRISHNAMURTI: int = 5  # K.S. Krishnamurti (KP)
SE_SIDM_DJWHAL_KHUL: int = 6  # Djwhal Khul (Alice Bailey)
SE_SIDM_YUKTESHWAR: int = 7  # Yukteshwar
SE_SIDM_JN_BHASIN: int = 8  # J.N. Bhasin

# Babylonian traditions
SE_SIDM_BABYL_KUGLER1: int = 9  # Kugler variant 1
SE_SIDM_BABYL_KUGLER2: int = 10  # Kugler variant 2
SE_SIDM_BABYL_KUGLER3: int = 11  # Kugler variant 3
SE_SIDM_BABYL_HUBER: int = 12  # Huber
SE_SIDM_BABYL_ETPSC: int = 13  # ETPSC
SE_SIDM_BABYL_BRITTON: int = 38  # Britton

# Star-based ayanamshas
SE_SIDM_ALDEBARAN_15TAU: int = 14  # Aldebaran at 15° Taurus
SE_SIDM_TRUE_CITRA: int = 27  # True position of Spica (180° Citra)
SE_SIDM_TRUE_REVATI: int = 28  # True position of Revati
SE_SIDM_TRUE_PUSHYA: int = 29  # True position of Pushya (Cancri)
SE_SIDM_TRUE_MULA: int = 35  # True position of Mula (λ Scorpii)
SE_SIDM_TRUE_SHEORAN: int = 39  # True Sheoran

# Historical epochs
SE_SIDM_HIPPARCHOS: int = 15  # Hipparchos (128 BC)  
SE_SIDM_SASSANIAN: int = 16  # Sassanian
SE_SIDM_J2000: int = 18  # J2000.0 (no ayanamsha)
SE_SIDM_J1900: int = 19  # J1900.0
SE_SIDM_B1950: int = 20  # B1950.0

# Indian traditions
SE_SIDM_SURYASIDDHANTA: int = 21  # Suryasiddhanta
SE_SIDM_SURYASIDDHANTA_MSUN: int = 22  # Suryasiddhanta (mean Sun)
SE_SIDM_ARYABHATA: int = 23  # Aryabhata
SE_SIDM_ARYABHATA_MSUN: int = 24  # Aryabhata (mean Sun)
SE_SIDM_ARYABHATA_522: int = 37  # Aryabhata 522
SE_SIDM_SS_REVATI: int = 25  # Suryasiddhanta Revati
SE_SIDM_SS_CITRA: int = 26  # Suryasiddhanta Citra

# Galactic alignment systems
SE_SIDM_GALCENT_0SAG: int = 17  # Galactic Center at 0° Sagittarius
SE_SIDM_GALCENT_RGILBRAND: int = 30  # Galactic Center (Gil Brand)
SE_SIDM_GALCENT_MULA_WILHELM: int = 36  # Galactic Center at Mula (Wilhelm)
SE_SIDM_GALCENT_COCHRANE: int = 40  # Galactic Center (Cochrane)
SE_SIDM_GALEQU_IAU1958: int = 31  # Galactic Equator (IAU 1958)
SE_SIDM_GALEQU_TRUE: int = 32  # Galactic Equator (True)
SE_SIDM_GALEQU_MULA: int = 33  # Galactic Equator at Mula
SE_SIDM_GALEQU_FIORENZA: int = 41  # Galactic Equator (Fiorenza)
SE_SIDM_GALALIGN_MARDYKS: int = 34  # Galactic Alignment (Mardyks)

# Other systems
SE_SIDM_VALENS_MOON: int = 42  # Vettius Valens (Moon-based)
SE_SIDM_USER: int = 255  # User-defined ayanamsha

# =============================================================================
# PYSWISSEPH-COMPATIBLE SIDEREAL MODE ALIASES (SIDM_* instead of SE_SIDM_*)
# =============================================================================
# pyswisseph uses SIDM_* prefix while Swiss Ephemeris C library uses SE_SIDM_*

# Western sidereal traditions
SIDM_FAGAN_BRADLEY: int = SE_SIDM_FAGAN_BRADLEY
SIDM_LAHIRI: int = SE_SIDM_LAHIRI
SIDM_DELUCE: int = SE_SIDM_DELUCE
SIDM_RAMAN: int = SE_SIDM_RAMAN
SIDM_USHASHASHI: int = SE_SIDM_USHASHASHI
SIDM_KRISHNAMURTI: int = SE_SIDM_KRISHNAMURTI
SIDM_DJWHAL_KHUL: int = SE_SIDM_DJWHAL_KHUL
SIDM_YUKTESHWAR: int = SE_SIDM_YUKTESHWAR
SIDM_JN_BHASIN: int = SE_SIDM_JN_BHASIN

# Babylonian traditions
SIDM_BABYL_KUGLER1: int = SE_SIDM_BABYL_KUGLER1
SIDM_BABYL_KUGLER2: int = SE_SIDM_BABYL_KUGLER2
SIDM_BABYL_KUGLER3: int = SE_SIDM_BABYL_KUGLER3
SIDM_BABYL_HUBER: int = SE_SIDM_BABYL_HUBER
SIDM_BABYL_ETPSC: int = SE_SIDM_BABYL_ETPSC
SIDM_BABYL_BRITTON: int = SE_SIDM_BABYL_BRITTON

# Star-based ayanamshas
SIDM_ALDEBARAN_15TAU: int = SE_SIDM_ALDEBARAN_15TAU
SIDM_TRUE_CITRA: int = SE_SIDM_TRUE_CITRA
SIDM_TRUE_REVATI: int = SE_SIDM_TRUE_REVATI
SIDM_TRUE_PUSHYA: int = SE_SIDM_TRUE_PUSHYA
SIDM_TRUE_MULA: int = SE_SIDM_TRUE_MULA
SIDM_TRUE_SHEORAN: int = SE_SIDM_TRUE_SHEORAN

# Historical epochs
SIDM_HIPPARCHOS: int = SE_SIDM_HIPPARCHOS
SIDM_SASSANIAN: int = SE_SIDM_SASSANIAN
SIDM_J2000: int = SE_SIDM_J2000
SIDM_J1900: int = SE_SIDM_J1900
SIDM_B1950: int = SE_SIDM_B1950

# Indian traditions
SIDM_SURYASIDDHANTA: int = SE_SIDM_SURYASIDDHANTA
SIDM_SURYASIDDHANTA_MSUN: int = SE_SIDM_SURYASIDDHANTA_MSUN
SIDM_ARYABHATA: int = SE_SIDM_ARYABHATA
SIDM_ARYABHATA_MSUN: int = SE_SIDM_ARYABHATA_MSUN
SIDM_ARYABHATA_522: int = SE_SIDM_ARYABHATA_522
SIDM_SS_REVATI: int = SE_SIDM_SS_REVATI
SIDM_SS_CITRA: int = SE_SIDM_SS_CITRA

# Galactic alignment systems
SIDM_GALCENT_0SAG: int = SE_SIDM_GALCENT_0SAG
SIDM_GALCENT_RGILBRAND: int = SE_SIDM_GALCENT_RGILBRAND
SIDM_GALCENT_MULA_WILHELM: int = SE_SIDM_GALCENT_MULA_WILHELM
SIDM_GALCENT_COCHRANE: int = SE_SIDM_GALCENT_COCHRANE
SIDM_GALEQU_IAU1958: int = SE_SIDM_GALEQU_IAU1958
SIDM_GALEQU_TRUE: int = SE_SIDM_GALEQU_TRUE
SIDM_GALEQU_MULA: int = SE_SIDM_GALEQU_MULA
SIDM_GALEQU_FIORENZA: int = SE_SIDM_GALEQU_FIORENZA
SIDM_GALALIGN_MARDYKS: int = SE_SIDM_GALALIGN_MARDYKS

# Other systems
SIDM_VALENS_MOON: int = SE_SIDM_VALENS_MOON
SIDM_USER: int = SE_SIDM_USER


# =============================================================================
# CALENDAR SYSTEMS
# =============================================================================

SE_JUL_CAL: int = 0  # Julian calendar
SE_GREG_CAL: int = 1  # Gregorian calendar

# =============================================================================
# ECLIPSE TYPES AND FLAGS
# =============================================================================
# Eclipse geometry types
SE_ECL_CENTRAL: int = 1  # Central eclipse (Moon's shadow axis crosses Earth)
SE_ECL_NONCENTRAL: int = 2  # Non-central eclipse
SE_ECL_TOTAL: int = 4  # Total eclipse (Sun/Earth completely covered)
SE_ECL_ANNULAR: int = 8  # Annular eclipse (ring of fire)
SE_ECL_PARTIAL: int = 16  # Partial eclipse
SE_ECL_ANNULAR_TOTAL: int = 32  # Hybrid eclipse (annular at some locations, total at others)
SE_ECL_PENUMBRAL: int = 64  # Penumbral lunar eclipse

# Composite eclipse type masks
SE_ECL_ALLTYPES_SOLAR: int = (
    SE_ECL_CENTRAL
    | SE_ECL_NONCENTRAL
    | SE_ECL_TOTAL
    | SE_ECL_ANNULAR
    | SE_ECL_PARTIAL
    | SE_ECL_ANNULAR_TOTAL
)
SE_ECL_ALLTYPES_LUNAR: int = SE_ECL_TOTAL | SE_ECL_PARTIAL | SE_ECL_PENUMBRAL

# Eclipse visibility and contact flags
SE_ECL_VISIBLE: int = 128  # Eclipse visible at location
SE_ECL_MAX_VISIBLE: int = 256  # Maximum phase visible
SE_ECL_1ST_VISIBLE: int = 512  # First contact visible
SE_ECL_2ND_VISIBLE: int = 1024  # Second contact visible
SE_ECL_3RD_VISIBLE: int = 2048  # Third contact visible
SE_ECL_4TH_VISIBLE: int = 4096  # Fourth contact visible
SE_ECL_ONE_TRY: int = 32768  # Try only once (optimization flag)
