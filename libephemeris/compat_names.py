# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Display strings the public name lookups answer with.

Four tables hold every label the API returns as a *name* rather than as a
number:

* :data:`HOUSE_SYSTEM_NAMES` maps a house-system letter to the label
  ``house_name`` returns.
* :data:`AYANAMSHA_NAMES` maps a predefined sidereal mode (``SIDM_*``) to the
  label ``get_ayanamsa_name`` returns.
* :data:`BODY_DISPLAY_NAMES` maps the identifier of a planet, lunar point,
  fictitious body or built-in asteroid to the label ``get_planet_name``
  returns.
* :data:`CONTRIB_HOUSE_SYSTEM_NAMES` maps a house-system letter to the label
  ``contrib.house_system_name`` returns: a second house table with its own
  wording and its own default.

Each table maps a selector to a string and nothing else. The lookup rules
stay with the entry points that own them: ``house_name`` normalises the
selector to a character and folds ``'g'`` to ``'G'``; ``get_ayanamsa_name``
masks the mode with ``0xFF`` so the ``SIDBIT_*`` projection flags leave the
name unchanged; ``get_planet_name`` falls through to the asteroid,
exotic-body and planetary-moon resolvers for identifiers outside its table;
and an unknown selector receives the entry point's own default (the empty
string, or ``"Unknown"`` for the contrib table).

Provenance:
    Interface labels kept for compatibility with the public API. Each string
    is the value a consumer receives from ``house_name``,
    ``get_ayanamsa_name``, ``get_planet_name`` or
    ``contrib.house_system_name``; charts and printed output built on those
    calls store the strings as they are, so every label is preserved byte for
    byte, including its spelling, capitalisation and punctuation. The keys are
    the public identifiers of ``libephemeris.constants`` and the house-system
    letters. The tables hold no numerical model, coefficient or computation:
    the definitions behind the names are registered with the modules that
    compute them (``houses``, ``ayanamsha_definitions``, ``planets`` and
    ``hypothetical``).
"""

from __future__ import annotations

from .constants import (
    ADMETOS,
    APOLLON,
    CERES,
    CHIRON,
    CUPIDO,
    EARTH,
    HADES,
    HARRINGTON,
    INTP_APOG,
    INTP_PERG,
    ISIS,
    JUNO,
    JUPITER,
    KRONOS,
    MARS,
    MEAN_APOG,
    MEAN_NODE,
    MERCURY,
    MOON,
    NEPTUNE,
    NEPTUNE_ADAMS,
    NEPTUNE_LEVERRIER,
    NIBIRU,
    OSCU_APOG,
    PALLAS,
    PHOLUS,
    PLUTO,
    PLUTO_LOWELL,
    PLUTO_PICKERING,
    POSEIDON,
    PROSERPINA,
    SATURN,
    SIDM_ALDEBARAN_15TAU,
    SIDM_ARYABHATA,
    SIDM_ARYABHATA_522,
    SIDM_ARYABHATA_MSUN,
    SIDM_B1950,
    SIDM_BABYL_BRITTON,
    SIDM_BABYL_ETPSC,
    SIDM_BABYL_HUBER,
    SIDM_BABYL_KUGLER1,
    SIDM_BABYL_KUGLER2,
    SIDM_BABYL_KUGLER3,
    SIDM_DELUCE,
    SIDM_DJWHAL_KHUL,
    SIDM_FAGAN_BRADLEY,
    SIDM_GALALIGN_MARDYKS,
    SIDM_GALCENT_0SAG,
    SIDM_GALCENT_COCHRANE,
    SIDM_GALCENT_MULA_WILHELM,
    SIDM_GALCENT_RGILBRAND,
    SIDM_GALEQU_FIORENZA,
    SIDM_GALEQU_IAU1958,
    SIDM_GALEQU_MULA,
    SIDM_GALEQU_TRUE,
    SIDM_HIPPARCHOS,
    SIDM_J1900,
    SIDM_J2000,
    SIDM_JN_BHASIN,
    SIDM_KRISHNAMURTI,
    SIDM_KRISHNAMURTI_VP291,
    SIDM_LAHIRI,
    SIDM_LAHIRI_1940,
    SIDM_LAHIRI_ICRC,
    SIDM_LAHIRI_VP285,
    SIDM_RAMAN,
    SIDM_SASSANIAN,
    SIDM_SS_CITRA,
    SIDM_SS_REVATI,
    SIDM_SURYASIDDHANTA,
    SIDM_SURYASIDDHANTA_MSUN,
    SIDM_TRUE_CITRA,
    SIDM_TRUE_MULA,
    SIDM_TRUE_PUSHYA,
    SIDM_TRUE_REVATI,
    SIDM_TRUE_SHEORAN,
    SIDM_USHASHASHI,
    SIDM_VALENS_MOON,
    SIDM_YUKTESHWAR,
    SUN,
    TRUE_NODE,
    URANUS,
    VENUS,
    VESTA,
    VULCAN,
    VULKANUS,
    WALDEMATH,
    WHITE_MOON,
    ZEUS,
)

__all__ = [
    "AYANAMSHA_NAMES",
    "BODY_DISPLAY_NAMES",
    "CONTRIB_HOUSE_SYSTEM_NAMES",
    "HOUSE_SYSTEM_NAMES",
]

#: House-system letter to the label ``house_name`` returns. ``'A'`` and
#: ``'E'`` share one label; ``'I'`` and ``'i'`` are two distinct systems.
#: ``house_name`` upper-cases the other lowercase letters and folds ``'g'``
#: to ``'G'`` before looking a selector up here.
HOUSE_SYSTEM_NAMES: dict[str, str] = {
    "P": "Placidus",
    "K": "Koch",
    "O": "Porphyry",
    "R": "Regiomontanus",
    "C": "Campanus",
    "E": "equal",
    "A": "equal",
    "W": "equal/ whole sign",
    "M": "Morinus",
    "B": "Alcabitius",
    "T": "Polich/Page",
    "U": "Krusinski-Pisa-Goelzer",
    "G": "Gauquelin sectors",
    "V": "equal/Vehlow",
    "X": "axial rotation system/Meridian houses",
    "H": "horizon/azimut",
    "F": "Carter poli-equ.",
    "S": "Sripati",
    "L": "Pullen SD",
    "Q": "Pullen SR",
    "N": "equal/1=Aries",
    "Y": "APC houses",
    "D": "equal (MC)",
    "I": "Sunshine",
    "i": "Sunshine/alt.",
    "J": "Savard-A",
}

#: Predefined sidereal mode to the label ``get_ayanamsa_name`` returns. The
#: table covers the modes 0-46; the user-defined mode ``SIDM_USER`` (255) and
#: the unassigned ids above the last predefined mode have no entry.
AYANAMSHA_NAMES: dict[int, str] = {
    SIDM_FAGAN_BRADLEY: "Fagan/Bradley",
    SIDM_LAHIRI: "Lahiri",
    SIDM_DELUCE: "De Luce",
    SIDM_RAMAN: "Raman",
    SIDM_USHASHASHI: "Usha/Shashi",
    SIDM_KRISHNAMURTI: "Krishnamurti",
    SIDM_DJWHAL_KHUL: "Djwhal Khul",
    SIDM_YUKTESHWAR: "Yukteshwar",
    SIDM_JN_BHASIN: "J.N. Bhasin",
    SIDM_BABYL_KUGLER1: "Babylonian/Kugler 1",
    SIDM_BABYL_KUGLER2: "Babylonian/Kugler 2",
    SIDM_BABYL_KUGLER3: "Babylonian/Kugler 3",
    SIDM_BABYL_HUBER: "Babylonian/Huber",
    SIDM_BABYL_ETPSC: "Babylonian/Eta Piscium",
    SIDM_BABYL_BRITTON: "Babylonian/Britton",
    SIDM_ALDEBARAN_15TAU: "Babylonian/Aldebaran = 15 Tau",
    SIDM_TRUE_CITRA: "True Citra",
    SIDM_TRUE_REVATI: "True Revati",
    SIDM_TRUE_PUSHYA: "True Pushya (PVRN Rao)",
    SIDM_TRUE_MULA: "True Mula (Chandra Hari)",
    SIDM_TRUE_SHEORAN: '"Vedic"/Sheoran',
    SIDM_HIPPARCHOS: "Hipparchos",
    SIDM_SASSANIAN: "Sassanian",
    SIDM_J2000: "J2000",
    SIDM_J1900: "J1900",
    SIDM_B1950: "B1950",
    SIDM_SURYASIDDHANTA: "Suryasiddhanta",
    SIDM_SURYASIDDHANTA_MSUN: "Suryasiddhanta, mean Sun",
    SIDM_ARYABHATA: "Aryabhata",
    SIDM_ARYABHATA_MSUN: "Aryabhata, mean Sun",
    SIDM_ARYABHATA_522: "Aryabhata 522",
    SIDM_SS_REVATI: "SS Revati",
    SIDM_SS_CITRA: "SS Citra",
    SIDM_GALCENT_0SAG: "Galact. Center = 0 Sag",
    SIDM_GALCENT_RGILBRAND: "Galactic Center (Gil Brand)",
    SIDM_GALCENT_MULA_WILHELM: "Dhruva/Gal.Center/Mula (Wilhelm)",
    SIDM_GALCENT_COCHRANE: "Cochrane (Gal.Center = 0 Cap)",
    SIDM_GALEQU_IAU1958: "Galactic Equator (IAU1958)",
    SIDM_GALEQU_TRUE: "Galactic Equator",
    SIDM_GALEQU_MULA: "Galactic Equator mid-Mula",
    SIDM_GALEQU_FIORENZA: "Galactic Equator (Fiorenza)",
    SIDM_GALALIGN_MARDYKS: "Skydram (Mardyks)",
    SIDM_VALENS_MOON: "Vettius Valens",
    SIDM_LAHIRI_1940: "Lahiri 1940",
    SIDM_LAHIRI_VP285: "Lahiri VP285",
    SIDM_KRISHNAMURTI_VP291: "Krishnamurti-Senthilathiban",
    SIDM_LAHIRI_ICRC: "Lahiri ICRC",
}

#: Body identifier to the label ``get_planet_name`` returns for the planets,
#: the lunar points, the fictitious bodies and the built-in asteroids. The
#: same labels name a body in the range and coverage error messages.
BODY_DISPLAY_NAMES: dict[int, str] = {
    SUN: "Sun",
    MOON: "Moon",
    MERCURY: "Mercury",
    VENUS: "Venus",
    MARS: "Mars",
    JUPITER: "Jupiter",
    SATURN: "Saturn",
    URANUS: "Uranus",
    NEPTUNE: "Neptune",
    PLUTO: "Pluto",
    MEAN_NODE: "mean Node",
    TRUE_NODE: "true Node",
    MEAN_APOG: "mean Apogee",
    OSCU_APOG: "osc. Apogee",
    INTP_APOG: "intp. Apogee",
    INTP_PERG: "intp. Perigee",
    EARTH: "Earth",
    ISIS: "Isis-Transpluto",
    NIBIRU: "Nibiru",
    HARRINGTON: "Harrington",
    # Historical trans-Uranian/trans-Neptunian predictions carry a
    # parenthetical annotation naming the planet each prediction targeted
    # (Le Verrier 1846 and Adams 1845-46 for Neptune; Lowell 1915 and
    # Pickering 1919 for Pluto).
    NEPTUNE_LEVERRIER: "Leverrier (Neptune)",
    NEPTUNE_ADAMS: "Adams (Neptune)",
    PLUTO_LOWELL: "Lowell (Pluto)",
    PLUTO_PICKERING: "Pickering (Pluto)",
    CUPIDO: "Cupido",
    HADES: "Hades",
    ZEUS: "Zeus",
    KRONOS: "Kronos",
    APOLLON: "Apollon",
    ADMETOS: "Admetos",
    # "Vulcanus" is the canonical Witte-Sieggrün (Hamburg School) spelling of
    # the seventh Uranian planet; the constant keeps its VULKANUS spelling
    # while the label carries the Hamburg spelling.
    VULKANUS: "Vulcanus",
    POSEIDON: "Poseidon",
    # Geocentric/heliocentric symbolic points the runtime models analytically
    # (ids 55-58). They resolve to their published astrological names because
    # calc_ut returns positions for them, so the label must agree.
    VULCAN: "Vulcan",
    WHITE_MOON: "Selena/White Moon",
    PROSERPINA: "Proserpina",
    WALDEMATH: "Waldemath",
    CHIRON: "Chiron",
    PHOLUS: "Pholus",
    CERES: "Ceres",
    PALLAS: "Pallas",
    JUNO: "Juno",
    VESTA: "Vesta",
}

#: House-system letter to the label ``contrib.house_system_name`` returns.
#: The contrib helper takes the letter as given (no case folding) and answers
#: ``"Unknown"`` for a letter that is not listed here.
CONTRIB_HOUSE_SYSTEM_NAMES: dict[str, str] = {
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
