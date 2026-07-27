# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Generate ``libephemeris/ayanamsha_definitions.py`` from audited definitions.

Provenance:
    Sole authoritative generator for the shared definition module. Each mode's
    value, epoch, civil-calendar conversion, and star/galactic anchor is tied to
    the source record in ``docs/reference/ayanamsha.md``. ERFA and the Vondrak
    chain supply public frame geometry; generator parameters are stated below.
    No target is fitted to another implementation.

Every predefined sidereal mode is represented by the source/status record in
``docs/reference/ayanamsha.md``, expressed as the pair

    (ayanamsha value at the defining epoch [deg], defining epoch [TT JD])

and propagated at runtime with the Vondrák Method-B accumulated precession —
the same contract as ``SIDM_USER``. For rows labelled Published or Geometry,
this script converts the cited statement or public catalogue construction into
that pair form. Rows labelled Secondary or Convention remain visibly labelled
in the generated comments; generation does not upgrade their evidence status.
No input value derives from another ephemeris implementation's output; the
derivation conventions are documented inline and in
``docs/reference/ayanamsha.md``.

Conventions:
- Equinox/solstice instants use the mean-equinox polynomial expressions of
  Meeus, *Astronomical Algorithms* (2nd ed.), chapter 27, tables 27.A/27.B.
  The neglected periodic terms move the instant by at most ~half a day,
  which changes an ayanamsha anchor by < 0.1 arcsec.
- Historical civil dates use the Julian calendar (proleptic where needed)
  at 12:00, treated as TT; the Delta-T offset at those epochs shifts the
  anchor by < 0.01 arcsec.
- Ujjain local mean time is UT + 5h03m (longitude 75°47' E).
- Star anchors are geometric mean places from this repository's catalog
  (Hipparcos/Gaia values), proper motion via ERFA, precessed with the
  project's Vondrák chain. Aberration and parallax are excluded from the
  *defining* anchor by construction.

Run:  uv run python scripts/generate_ayanamsha_definitions.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import erfa  # noqa: E402

from libephemeris.precession_vondrak import (  # noqa: E402
    method_b_accumulated_precession,
    vondrak_precession_matrix,
)
from libephemeris.sidereal_longterm import mean_obliquity_rad  # noqa: E402

J2000 = 2451545.0
B1900_TT = 2415020.31352  # Besselian epoch 1900.0 (Lieske 1979)
B1950_TT = 2433282.42345905  # Besselian epoch 1950.0 (Lieske 1979)
J1900_TT = 2415020.0  # 1900 January 0.5
UJJAIN_LMT_OFFSET_DAYS = (5.0 + 3.0 / 60.0) / 24.0  # 75°47' E


def dms(d: float, m: float = 0.0, s: float = 0.0) -> float:
    """Sexagesimal degrees to decimal degrees (sign taken from ``d``)."""
    sign = -1.0 if d < 0 else 1.0
    return sign * (abs(d) + m / 60.0 + s / 3600.0)


def julday_jd(year: int, month: int, day: int, hour: float, julian: bool) -> float:
    """Julian Day (TT) for a civil date, Julian or Gregorian calendar."""
    if julian:
        y = year
        m = month
        if m <= 2:
            y -= 1
            m += 12
        jd = (
            math.floor(365.25 * (y + 4716))
            + math.floor(30.6001 * (m + 1))
            + day
            - 1524.5
            + hour / 24.0
        )
        return jd
    dd, tf = erfa.dtf2d("TT", year, month, day, int(hour), int((hour % 1) * 60), 0.0)
    return float(dd) + float(tf)


def mean_march_equinox_jd(year: float) -> float:
    """Mean March-equinox instant (TT JD), Meeus ch. 27 (27.A / 27.B)."""
    if year < 1000:
        y = year / 1000.0
        return (
            1721139.29189
            + 365242.13740 * y
            + 0.06134 * y**2
            + 0.00111 * y**3
            - 0.00071 * y**4
        )
    y = (year - 2000.0) / 1000.0
    return (
        2451623.80984
        + 365242.37404 * y
        + 0.05169 * y**2
        - 0.00411 * y**3
        - 0.00057 * y**4
    )


def sun_mean_longitude_deg(jd_tt: float) -> float:
    """Geometric mean longitude of the Sun, equinox of date (Meeus 25.2)."""
    t = (jd_tt - 2451545.0) / 36525.0
    return (280.46646 + 36000.76983 * t + 0.0003032 * t * t) % 360.0


def nutation_longitude_deg(jd_tt: float) -> float:
    """IAU 2000A nutation in longitude at ``jd_tt``, in degrees.

    Used to convert a published TRUE (nutation-bearing) ayanamsha anchor to
    the MEAN value stored in the defining table: mean = true - dpsi(t0).
    """
    dpsi, _ = erfa.nut06a(J2000, jd_tt - J2000)
    return math.degrees(float(dpsi))


def mean_september_equinox_jd(year: float) -> float:
    """Mean September-equinox instant (TT JD), Meeus ch. 27 (27.B)."""
    y = (year - 2000.0) / 1000.0
    return (
        2451810.21715
        + 365242.01767 * y
        - 0.11575 * y**2
        + 0.00337 * y**3
        + 0.00078 * y**4
    )


def mean_december_solstice_jd(year: float) -> float:
    """Mean December-solstice instant (TT JD), Meeus ch. 27 (27.A / 27.B)."""
    if year < 1000:
        y = year / 1000.0
        return (
            1721414.39987
            + 365242.88257 * y
            - 0.00769 * y**2
            - 0.00933 * y**3
            - 0.00006 * y**4
        )
    y = (year - 2000.0) / 1000.0
    return (
        2451900.05952
        + 365242.74049 * y
        - 0.06223 * y**2
        - 0.00823 * y**3
        + 0.00032 * y**4
    )


def _star_mean_equatorial_of_date(
    star, jd_tt: float
) -> tuple[float, float, float, float]:
    """Mean-equator-of-date direction of a catalog star (geometric).

    Proper motion is applied with ERFA ``pmsafe`` (parallax and radial
    velocity from the catalog; a pure angular defining anchor), and the
    direction is precessed with the project's Vondrák matrix.

    Returns ``(x, y, z, eps)`` — unit vector on the mean equator of date and
    the pole-angle mean obliquity in radians.
    """
    ra = math.radians(star.ra_j2000)
    dec = math.radians(star.dec_j2000)
    pmr = math.radians(star.pm_ra / 3600.0) / math.cos(dec)
    pmd = math.radians(star.pm_dec / 3600.0)
    px = star.parallax if star.parallax > 0 else 1e-7
    rv = star.radial_velocity
    ra2, dec2, *_ = erfa.pmsafe(ra, dec, pmr, pmd, px, rv, J2000, 0.0, jd_tt, 0.0)
    v = erfa.s2c(float(ra2), float(dec2))
    import numpy as np

    vq = np.asarray(vondrak_precession_matrix(jd_tt)) @ np.asarray(v)
    eps = mean_obliquity_rad(jd_tt)
    return float(vq[0]), float(vq[1]), float(vq[2]), eps


def star_mean_ecliptic_longitude(star, jd_tt: float) -> float:
    """Geometric mean-ecliptic-of-date longitude of a catalog star."""
    xq, yq, zq, eps = _star_mean_equatorial_of_date(star, jd_tt)
    ce, se = math.cos(eps), math.sin(eps)
    ye = yq * ce + zq * se
    return math.degrees(math.atan2(ye, xq)) % 360.0


def star_polar_longitude(star, jd_tt: float) -> float:
    """Surya-Siddhanta polar longitude (dhruvaka) of a catalog star.

    The polar longitude is the ecliptic longitude of the point where the
    star's hour circle (declination circle through the celestial poles)
    cuts the ecliptic: ``tan L = tan(alpha) / cos(eps)`` with the mean RA
    of date. See Burgess, *Translation of the Surya-Siddhanta*, ch. VIII.
    """
    xq, yq, zq, eps = _star_mean_equatorial_of_date(star, jd_tt)
    alpha = math.atan2(yq, xq)
    lon = math.atan2(math.sin(alpha), math.cos(alpha) * math.cos(eps))
    return math.degrees(lon) % 360.0


def _stars():
    """Catalog star registry shared with the runtime (single source)."""
    from libephemeris.planets import STARS

    return STARS


# --------------------------------------------------------------------------
# Registered defining statements and conventions
#
# Each entry below carries its own evidence status in the generated source text
# and docs/reference/ayanamsha.md. This heading must not be read as claiming
# that every historical mode has a recovered primary statement from its author.
# --------------------------------------------------------------------------
E_499 = mean_march_equinox_jd(499.0)  # Kali 3600 elapsed: 21 March 499 CE

# Kali-yuga epoch and Siddhantic year lengths (published constants):
# - The Kali era begins at midnight (Ardharatrika school) or sunrise
#   (Audayika school) of 18 February 3102 BCE proleptic-Julian at
#   Ujjayini, longitude 75\N{DEGREE SIGN}47' E (Burgess; standard in the
#   history-of-astronomy literature). Greenwich midnight of that date is
#   JD 588465.5, so Ujjayini local mean midnight is 75.7833\N{DEGREE SIGN}/360
#   of a day earlier and sunrise is a quarter day after local midnight.
# - Surya Siddhanta sidereal year: 1,577,917,828 civil days per 4,320,000
#   solar revolutions (Burgess, ch. I).
# - Aryabhatiya sidereal year: 1,577,917,500 days per 4,320,000 revolutions
#   (Clark, Gitikapada 3).
# The "mean Sun" modes place the sidereal zero at the mean Sun of the
# Kali + 3600-year instant: the Siddhantic mean Sun starts at the zero
# point at the Kali epoch and completes exactly 3600 revolutions by that
# moment. Day counts (ahargana) are civil days, so the instants are UT;
# they are converted to TT with the project's Delta-T model before the
# mean solar longitude (Meeus 25.2) is evaluated.
UJJAYINI_LON_DEG = 75.0 + 47.0 / 60.0
KALI_EPOCH_MIDNIGHT_UT = 588465.5 - UJJAYINI_LON_DEG / 360.0
KALI_EPOCH_SUNRISE_UT = KALI_EPOCH_MIDNIGHT_UT + 0.25
SS_YEAR_DAYS = 1_577_917_828 / 4_320_000
ARYABHATA_YEAR_DAYS = 1_577_917_500 / 4_320_000


def _kali_3600_msun_pair(kali_epoch_ut: float, year_days: float) -> tuple[float, float]:
    """Defining (value, epoch_tt) for a Siddhantic mean-Sun ayanamsha."""
    from libephemeris import deltat

    epoch_ut = kali_epoch_ut + 3600.0 * year_days
    epoch_tt = epoch_ut + deltat(epoch_ut)
    lon = sun_mean_longitude_deg(epoch_tt)
    value = (lon + 180.0) % 360.0 - 180.0
    return value, epoch_tt


V_SS_MSUN, E_SS_MSUN = _kali_3600_msun_pair(KALI_EPOCH_MIDNIGHT_UT, SS_YEAR_DAYS)
V_ARYABHATA_MSUN, E_ARYABHATA_MSUN = _kali_3600_msun_pair(
    KALI_EPOCH_SUNRISE_UT, ARYABHATA_YEAR_DAYS
)

# Raman's published rule is linear: (year - 397) x 50 1/3 arcsec. The pair
# form anchors the rule in the era of his own tables (1900.0, Julian-epoch
# year count) and propagates with Method-B like every other mode.
RAMAN_AT_1900 = (1900.0 - 397.0) * (50.0 + 1.0 / 3.0) / 3600.0

DEFINITIONS: list[tuple[int, str, float, float, str]] = [
    # (mode, name, value_deg, epoch_tt_jd, source)
    (
        0,
        "FAGAN_BRADLEY",
        360.0 - dms(335, 57, 28.64),
        B1950_TT,
        "Fagan & Firebrace, Primer of Sidereal Astrology: SVP 335\N{DEGREE SIGN}57'28.64\" at B1950.0",
    ),
    (
        1,
        "LAHIRI",
        dms(23, 51, 25.53),
        2451545.0,
        "Indian Astronomical Ephemeris (Positional Astronomy Centre, Kolkata):"
        " ayanamsa 23\N{DEGREE SIGN}51'25\".53 at 2000 Jan 1.5 TT (J2000.0)",
    ),
    (
        2,
        "DELUCE",
        dms(26, 24, 47),
        B1900_TT,
        "De Luce, Constellational Astrology (1963) p.5: 26\N{DEGREE SIGN}24'47\" in 1900",
    ),
    (
        3,
        "RAMAN",
        RAMAN_AT_1900,
        J1900_TT,
        'Raman, Hindu Predictive Astrology App. A: (year-397) x 50 1/3", anchored at 1900.0',
    ),
    (
        4,
        "USHASHASHI",
        0.0,
        julday_jd(559, 3, 21, 12.0, julian=True),
        "Usha & Shashi (1978); traditional zero epoch 559 CE (primary transcription pending)",
    ),
    (
        5,
        "KRISHNAMURTI",
        dms(22, 21, 50),
        J1900_TT,
        "Krishnamurti Paddhati: 22\N{DEGREE SIGN}21'50\" at 1900 January 0.5",
    ),
    (
        6,
        "DJWHAL_KHUL",
        30.0,
        julday_jd(2117, 1, 1, 12.0, julian=False),
        "Bailey/Ageless Wisdom: Aquarian age (ayanamsha 30\N{DEGREE SIGN}) begins 2117",
    ),
    (
        7,
        "YUKTESHWAR",
        dms(20, 54, 36),
        mean_march_equinox_jd(1894.0),
        "Sri Yukteswar, The Holy Science (1894): 20\N{DEGREE SIGN}54'36\" at the spring 1894 equinox",
    ),
    (
        8,
        "JN_BHASIN",
        0.0,
        julday_jd(364, 3, 21, 12.0, julian=True),
        "J.N. Bhasin; traditional zero epoch 364 CE (primary transcription pending)",
    ),
    (
        9,
        "BABYL_KUGLER1",
        0.0,
        julday_jd(309, 3, 21, 12.0, julian=True),
        "Kugler, Sternkunde II, solution 1; zero epoch 309 CE (primary transcription pending)",
    ),
    (
        10,
        "BABYL_KUGLER2",
        0.0,
        julday_jd(208, 3, 21, 12.0, julian=True),
        "Kugler, Sternkunde II, solution 2; zero epoch 208 CE (primary transcription pending)",
    ),
    (
        11,
        "BABYL_KUGLER3",
        0.0,
        julday_jd(146, 3, 21, 12.0, julian=True),
        "Kugler, Sternkunde II, solution 3; zero epoch 146 CE (primary transcription pending)",
    ),
    (
        12,
        "BABYL_HUBER",
        -dms(4, 28),
        julday_jd(-100, 3, 21, 12.0, julian=True),
        "Huber, Centaurus 5 (1958): vernal point at Babylonian Aries 4\N{DEGREE SIGN}28' in -100",
    ),
    (
        13,
        "BABYL_ETPSC",
        0.0,
        julday_jd(237, 3, 21, 12.0, julian=True),
        "eta Piscium Babylonian norm (Mercier); zero epoch 237 CE (primary transcription pending)",
    ),
    (
        15,
        "HIPPARCHOS",
        0.0,
        julday_jd(545, 3, 21, 12.0, julian=True),
        "Mercier, Studies in the Medieval Conception of Precession: Hipparchan norm zero 545 CE",
    ),
    (
        16,
        "SASSANIAN",
        0.0,
        julday_jd(564, 3, 18, 12.0, julian=True),
        "Mercier: Sassanian (Zij-i Shah) zero epoch 18 March 564 CE",
    ),
    (18, "J2000", 0.0, 2451545.0, "Mean ecliptic/equinox of J2000 (frame epoch)"),
    (
        19,
        "J1900",
        0.0,
        J1900_TT,
        "Mean ecliptic/equinox of 1900 January 0.5 (frame epoch)",
    ),
    (20, "B1950", 0.0, B1950_TT, "Mean ecliptic/equinox of B1950.0 (frame epoch)"),
    (
        21,
        "SURYASIDDHANTA",
        0.0,
        E_499,
        "Surya Siddhanta (Burgess): zero at the Kali-3600 equinox, 21 March 499 CE",
    ),
    (
        22,
        "SURYASIDDHANTA_MSUN",
        V_SS_MSUN,
        E_SS_MSUN,
        "Surya Siddhanta mean Sun (Burgess): mean Sun at the sidereal zero,"
        " Ardharatrika Kali epoch (Ujjayini midnight) + 3600 SS years",
    ),
    (
        23,
        "ARYABHATA",
        0.0,
        E_499,
        "Aryabhatiya (Clark): zero at the Kali-3600 equinox, 21 March 499 CE",
    ),
    (
        24,
        "ARYABHATA_MSUN",
        V_ARYABHATA_MSUN,
        E_ARYABHATA_MSUN,
        "Aryabhatiya mean Sun (Clark): mean Sun at the sidereal zero,"
        " Audayika Kali epoch (Ujjayini sunrise) + 3600 Aryabhata years",
    ),
    (
        25,
        "SS_REVATI",
        None,  # filled in main(): polar longitude of zeta Psc at E_499 - 359 deg 50'
        E_499,
        "Surya Siddhanta VIII: Revati (zeta Psc) at polar longitude 359\N{DEGREE SIGN}50' at the Kali-3600 epoch",
    ),
    (
        26,
        "SS_CITRA",
        None,  # filled in main(): polar longitude of Spica at E_499 - 180 deg
        E_499,
        "Surya Siddhanta VIII: Citra (Spica) at polar longitude 180\N{DEGREE SIGN} at the Kali-3600 epoch",
    ),
    (
        37,
        "ARYABHATA_522",
        0.0,
        julday_jd(522, 3, 21, 12.0, julian=True),
        "Aryabhata 522 CE tradition: zero at the 522 CE spring equinox",
    ),
    (
        38,
        "BABYL_BRITTON",
        0.0,
        julday_jd(230, 3, 21, 12.0, julian=True),
        "Britton, AHES 64 (2010) uniform zodiac; zero epoch 230 CE (primary transcription pending)",
    ),
    # Mode 39 (TRUE_SHEORAN) is NOT an epoch mode: Sheoran's rigorous
    # definition anchors the sidereal zero to the live star Asellus Australis
    # (delta Cancri) at a fixed sidereal longitude, making "the actual rate of
    # precession irrelevant" (The Science of Time, 2017). It is emitted below
    # as a dynamic star mode with SHEORAN_TARGET_LON; the book's -60 deg at
    # the 4174 BCE winter solstice and the 1 deg / 71.75 yr rule are the
    # author's own derived checkpoints/approximation of that anchor.
    (
        41,
        "GALEQU_FIORENZA",
        25.0,
        julday_jd(2000, 1, 1, 0.0, julian=False),
        "Fiorenza: galactic-equator ayanamsha 25\N{DEGREE SIGN} at 2000-01-01",
    ),
    (
        43,
        "LAHIRI_1940",
        0.0,
        julday_jd(286, 3, 21, 12.0, julian=True),
        "Lahiri (1940) tradition; zero epoch 286 CE (primary transcription pending)",
    ),
    (
        44,
        "LAHIRI_VP285",
        0.0,
        mean_march_equinox_jd(285.0),
        "Vernal-point convention: zero at the 285 CE spring equinox",
    ),
    (
        45,
        "KRISHNAMURTI_VP291",
        0.0,
        mean_march_equinox_jd(291.0),
        "Vernal-point convention: zero at the 291 CE spring equinox",
    ),
    (
        46,
        "LAHIRI_ICRC",
        # The committee's 23 deg 15'0" is a Citra-referenced value of the same
        # family as the Indian Astronomical Ephemeris series, which publishes
        # TRUE (nutation-bearing) ayanamshas (Senthilathiban, "A Study of KP
        # Ayanamsa with Modern Precession Theories", p. 112). The defining
        # pair stores MEAN values, so the anchor is converted at its own
        # epoch with the IAU 2000A nutation in longitude (public model, no
        # output fitting): mean = 23 deg 15' - dpsi(t0).
        dms(23, 15, 0) - nutation_longitude_deg(julday_jd(1956, 3, 21, 0.0, False)),
        julday_jd(1956, 3, 21, 0.0, julian=False),
        "Calendar Reform Committee (1955): true ayanamsha 23\N{DEGREE SIGN}15' "
        "at 1956-03-21 00:00 ET, stored as mean (- dpsi at the epoch)",
    ),
]

GOLDEN_RATIO = (1.0 + math.sqrt(5.0)) / 2.0
GALCENT_TARGETS = {
    17: (240.0, "Galactic Center at 0\N{DEGREE SIGN} Sagittarius"),
    30: (
        240.0 + 30.0 / GOLDEN_RATIO**4,
        "Gil Brand, Himmlische Matrix: GC at the double golden section of "
        "Sagittarius, 240\N{DEGREE SIGN} + 30\N{DEGREE SIGN}/phi^4",
    ),
    40: (270.0, "Cochrane: Galactic Center at 0\N{DEGREE SIGN} Capricorn"),
}
MULA_WILHELM_TARGET = 246.0 + 40.0 / 60.0  # middle of Mula: 246\N{DEGREE SIGN}40'
MARDYKS_EPOCH = mean_september_equinox_jd(1998.0)


HEADER = '''# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Shared defining data for predefined sidereal modes.

Provenance:
    This is a generated scientific artifact. Its sole authoritative generator
    is ``scripts/generate_ayanamsha_definitions.py``; every emitted entry maps
    to the exact source/status row in ``docs/reference/ayanamsha.md``.
    Published and geometric definitions remain distinct from rows labelled
    Secondary or Convention. No zero point is fitted to compatibility output.
    Regeneration preserves this evidence boundary.

AUTO-GENERATED by scripts/generate_ayanamsha_definitions.py - edit that
script, not this module. Every entry carries a concise source description; its
full evidence status and derivation conventions are documented in the
generator and in docs/reference/ayanamsha.md. No value is fitted to another
ephemeris implementation's output.

Runtime model for epoch modes (identical to the ``SIDM_USER`` contract):

    ayanamsha(t) = value + method_b_accumulated_precession(t, epoch)

with the Vondrak, Capitaine & Wallace (2011) Method-B accumulation on the
mean ecliptic of the defining epoch. Dynamic star and galactic modes are
evaluated from live catalog geometry in planets.py against the published
target longitudes below.
"""

from __future__ import annotations

'''


def main() -> None:
    """Resolve dynamic anchors and emit the authoritative ayanamsha module.

    Publication-defined constant modes are transcribed from ``DEFINITIONS``.
    Revati and Spica anchor values are computed at their defining epochs from
    the public generated star catalogue and the documented polar-longitude
    geometry. The output retains a source comment for every mode and delegates
    live star/galactic modes to runtime geometry rather than freezing a foreign
    implementation's values.
    """
    stars = _stars()
    for i, (mode, name, value, epoch, source) in enumerate(DEFINITIONS):
        if value is None:
            if mode == 25:
                value = star_polar_longitude(stars["REVATI"], epoch) - dms(359, 50)
            elif mode == 26:
                value = star_polar_longitude(stars["SPICA"], epoch) - 180.0
            DEFINITIONS[i] = (mode, name, value, epoch, source)

    lines = [HEADER]
    lines.append(
        "# mode: (ayanamsha at defining epoch [deg], defining epoch [TT JD])\n"
    )
    lines.append("AYANAMSHA_DEFINING: dict[int, tuple[float, float]] = {\n")
    for mode, name, value, epoch, source in DEFINITIONS:
        lines.append(f"    # {name}: {source}\n")
        lines.append(f"    {mode}: ({value!r}, {epoch!r}),\n")
    lines.append("}\n\n")

    lines.append(
        "# Modes evaluated from live catalog/frame geometry in planets.py.\n"
        "DYNAMIC_AYANAMSHA_MODES = frozenset(\n"
        "    {14, 17, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 39, 40, 42}\n"
        ")\n\n"
    )

    lines.append(
        "# Galactic-center modes: published sidereal target longitude of Sgr A*.\n"
        "GALCENT_TARGET_LON: dict[int, float] = {\n"
    )
    for mode, (target, source) in GALCENT_TARGETS.items():
        lines.append(f"    # {source}\n")
        lines.append(f"    {mode}: {target!r},\n")
    lines.append("}\n\n")

    lines.append(
        "# Wilhelm: Sgr A* hour-circle (dhruva) projection at the middle of Mula\n"
        f"# (246\N{DEGREE SIGN}40').\n"
        f"MULA_WILHELM_TARGET_LON = {MULA_WILHELM_TARGET!r}\n\n"
    )
    lines.append(
        "# Aldebaran at 15\N{DEGREE SIGN} Taurus (45\N{DEGREE SIGN} absolute): the "
        "star-fixed anchor of the\n"
        '# "Aldebaran = 15 Tau" sidereal convention, the Babylonian '
        "exaltation-degree\n"
        '# tradition revived by Cyril Fagan, "Zodiacs Old and New" (Llewellyn, '
        "1950).\n"
        "# It defines sidereal mode SIDM_ALDEBARAN_15TAU (14) and is distinct "
        "from the\n"
        "# Fagan-Bradley SVP ayanamsha, from which it differs by ~1'.\n"
        "ALDEBARAN_TARGET_LON = 45.0\n\n"
    )
    lines.append(
        "# Sheoran, The Science of Time (2017): the sidereal zero is anchored\n"
        "# to the live star Asellus Australis (delta Cancri) held at the fixed\n"
        "# sidereal longitude 103\N{DEGREE SIGN}29'32.9375\", chosen so that "
        "the Sun stands at\n"
        "# sidereal 330\N{DEGREE SIGN} (Pisces 0) at the 4174 BCE winter "
        "solstice; the star\n"
        '# anchor makes "the actual rate of precession irrelevant" (the\n'
        "# 1\N{DEGREE SIGN}/71.75 yr rule is the author's stated "
        "approximation only). It\n"
        "# defines sidereal mode SIDM_TRUE_SHEORAN (39).\n"
        f"SHEORAN_TARGET_LON = {dms(103, 29, 32.9375)!r}\n\n"
    )
    lines.append(
        "# Mardyks, Sacred Astronomy (1991): ayanamsha exactly 30\N{DEGREE SIGN} at\n"
        "# the September equinox 1998 (Meeus mean-equinox instant); fixed-epoch\n"
        "# frame mode propagated with the IAU 2006 general-precession angle.\n"
        f"MARDYKS_DEFINING: tuple[float, float] = (30.0, {MARDYKS_EPOCH!r})\n\n"
    )
    lines.append(
        "# Vettius Valens lunar-epoch mode: defining UT epoch and offset.\n"
        "# Offset consistent with the ~-3 deg Hellenistic sidereal norm reported\n"
        "# in the scholarship (e.g. GRBS 56 (2016) 573-575 discussion of Valens'\n"
        "# longitudes); primary transcription of Anthologiae I.18 pending.\n"
        "VALENS_MOON_T0_UT = 1775845.5\n"
        "VALENS_MOON_AYAN_T0_DEG = -2.9422\n"
    )

    out = (
        Path(__file__).resolve().parents[1]
        / "libephemeris"
        / "ayanamsha_definitions.py"
    )
    text = "".join(lines)
    out.write_text(text, encoding="utf-8")
    print(f"wrote {out} ({len(text)} bytes)")

    # Delta report vs the propagated value at J2000 (diagnostic only).
    print("\nmode  name                 value@epoch      J2000 equivalent")
    for mode, name, value, epoch, _src in DEFINITIONS:
        j2000_val = value + method_b_accumulated_precession(J2000, epoch)
        print(f"{mode:>4}  {name:<20} {value:>14.7f}  {j2000_val:>14.7f}")


if __name__ == "__main__":
    main()
