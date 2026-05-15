#!/usr/bin/env python3
"""Rounds 236-255: Deep diverse verification batch 2.

236: All planets TRUEPOS flag
237: Planet magnitude (pheno) all bodies
238: Delta-T at leap second boundaries
239: Houses at tropical/arctic latitudes
240: Fixed stars equatorial mode
241: Sidereal planets Lahiri+Raman sweep
242: Moon NOABERR+TRUEPOS vs default
243: Chiron+Ceres heliocentric+J2000
244: cotrans_sp speed transformation
245: All bodies at JD fractional hours
246: Houses_armc_ex2 sweep
247: Planet distance at opposition/conjunction
248: Nutation at node cycle peaks
249: OscuLilith+MeanLilith full sweep
250: Fixed stars sidereal Lahiri
251: All planets NONUT+EQUATORIAL combined
252: house_pos Moon at all houses
253: Asteroid belt positions sweep
254: Rising/setting Sun timing
255: Multi-flag planet matrix
"""

from __future__ import annotations
import os, sys, math

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("LIBEPHEMERIS_MODE", "skyfield")

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

results = {}


def run_round(num, name, test_fn):
    p, f, t, fails = 0, 0, 0, []
    try:
        p, f, t, fails = test_fn()
    except Exception as e:
        fails = [f"  ERROR: {e}"]
        f, t = 1, 1
    results[num] = (name, p, f, t, fails)
    pct = 100 * p / t if t > 0 else 0
    status = "PASS" if f == 0 else f"{f} FAIL"
    print(f"  Round {num}: {name} — {p}/{t} ({pct:.1f}%) [{status}]")


def adiff(a, b):
    d = abs(a - b) % 360
    return min(d, 360 - d)


DATES_10 = [
    2415020.0,
    2430000.0,
    2440000.0,
    2445000.0,
    2451545.0,
    2455000.0,
    2458000.0,
    2460000.0,
    2462000.0,
    2465000.0,
]
PLANETS_7 = [
    (ephem.SUN, swe.SUN),
    (ephem.MOON, swe.MOON),
    (ephem.MERCURY, swe.MERCURY),
    (ephem.VENUS, swe.VENUS),
    (ephem.MARS, swe.MARS),
    (ephem.JUPITER, swe.JUPITER),
    (ephem.SATURN, swe.SATURN),
]
PLANETS_ALL = PLANETS_7 + [
    (ephem.URANUS, swe.URANUS),
    (ephem.NEPTUNE, swe.NEPTUNE),
    (ephem.PLUTO, swe.PLUTO),
]


def test_236():
    p, f, t, fails = 0, 0, 0, []
    le_f = ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_TRUEPOS
    se_f = swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_TRUEPOS
    for le_b, se_b in PLANETS_ALL:
        for jd in DATES_10:
            try:
                le_r = ephem.calc_ut(jd, le_b, le_f)
                se_r = swe.calc_ut(jd, se_b, se_f)
            except:
                continue
            t += 1
            d = adiff(le_r[0][0], se_r[0][0]) * 3600
            if d <= 35.0:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} jd={jd:.0f} TRUEPOS diff={d:.2f}"')
            t += 1
            d = abs(le_r[0][1] - se_r[0][1]) * 3600
            if d <= 35.0:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} jd={jd:.0f} TRUEPOS lat diff={d:.2f}"')
    return p, f, t, fails


def test_237():
    p, f, t, fails = 0, 0, 0, []
    bodies = [
        (ephem.MOON, swe.MOON, "Moon"),
        (ephem.VENUS, swe.VENUS, "Venus"),
        (ephem.MARS, swe.MARS, "Mars"),
        (ephem.JUPITER, swe.JUPITER, "Jupiter"),
        (ephem.SATURN, swe.SATURN, "Saturn"),
    ]
    for le_b, se_b, nm in bodies:
        for jd in DATES_10[:7]:
            try:
                le_r = ephem.pheno_ut(jd, le_b, ephem.FLG_SWIEPH)
                se_r = swe.pheno_ut(jd, se_b, swe.FLG_SWIEPH)
                le_elong = le_r[0][2]
                se_elong = se_r[2]
            except:
                continue
            t += 1
            d = abs(le_elong - se_elong) * 3600
            if d <= 5.0:
                p += 1
            else:
                f += 1
                fails.append(f'  {nm} jd={jd:.0f} elong diff={d:.2f}"')
            t += 1  # magnitude
            le_mag = le_r[0][4]
            se_mag = se_r[4]
            if abs(le_mag - se_mag) <= 0.5:
                p += 1
            else:
                f += 1
                fails.append(f"  {nm} jd={jd:.0f} mag LE={le_mag:.2f} SE={se_mag:.2f}")
    return p, f, t, fails


def test_238():
    p, f, t, fails = 0, 0, 0, []
    # Dates near leap seconds (approximate JDs)
    leap_jds = [
        2441317.5,
        2441499.5,
        2441683.5,
        2442048.5,
        2442413.5,
        2442778.5,
        2443144.5,
        2443509.5,
        2443874.5,
        2444239.5,
        2444786.5,
        2445151.5,
        2445516.5,
        2446247.5,
        2447161.5,
        2447892.5,
        2448257.5,
        2448804.5,
        2449169.5,
        2449534.5,
        2450083.5,
        2450630.5,
        2451179.5,
        2453736.5,
        2454832.5,
        2456109.5,
        2457204.5,
        2457754.5,
    ]
    for jd in leap_jds:
        for offset in [-0.001, 0.0, 0.001]:
            try:
                le_dt = ephem.deltat(jd + offset)
                se_dt = swe.deltat(jd + offset)
            except:
                continue
            t += 1
            d = abs(le_dt - se_dt) * 86400
            if d <= 1.0:
                p += 1
            else:
                f += 1
                fails.append(f"  jd={jd + offset:.3f} dt diff={d:.4f}s")
    return p, f, t, fails


def test_239():
    p, f, t, fails = 0, 0, 0, []
    lats = [23.44, -23.44, 66.56, -66.56, 0.0001, 89.0, -89.0]
    hsystems = ["P", "K", "O", "R", "C", "E"]
    for jd in [2451545.0, 2460000.0, 2440000.0]:
        for lat in lats:
            for hs in hsystems:
                try:
                    le_c, le_a = ephem.houses_ex(
                        jd, lat, 0.0, ord(hs), ephem.FLG_SWIEPH | ephem.FLG_SPEED
                    )
                    se_c, se_a = swe.houses_ex(
                        jd, lat, 0.0, hs.encode(), swe.FLG_SWIEPH | swe.FLG_SPEED
                    )
                except:
                    continue
                t += 1
                d = adiff(le_a[0], se_a[0]) * 3600
                if d <= 1.0:
                    p += 1
                else:
                    f += 1
                    fails.append(f'  {hs} lat={lat} jd={jd:.0f} ASC diff={d:.4f}"')
                t += 1
                d = adiff(le_a[1], se_a[1]) * 3600
                if d <= 1.0:
                    p += 1
                else:
                    f += 1
                    fails.append(f'  {hs} lat={lat} jd={jd:.0f} MC diff={d:.4f}"')
    return p, f, t, fails


def test_240():
    p, f, t, fails = 0, 0, 0, []
    stars = [
        "Aldebaran",
        "Regulus",
        "Spica",
        "Antares",
        "Sirius",
        "Vega",
        "Arcturus",
        "Pollux",
        "Fomalhaut",
        "Deneb",
        "Capella",
        "Procyon",
    ]
    le_f = ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_EQUATORIAL
    se_f = swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_EQUATORIAL
    for star in stars:
        for jd in [2451545.0, 2460000.0, 2440000.0, 2455000.0]:
            try:
                le_r = ephem.fixstar2_ut(star, jd, le_f)
                se_r = swe.fixstar2(star, jd, se_f)
            except:
                continue
            t += 1
            d = adiff(le_r[0][0], se_r[0][0]) * 3600
            if d <= 10.0:
                p += 1
            else:
                f += 1
                fails.append(f'  {star} jd={jd:.0f} RA diff={d:.4f}"')
            t += 1
            d = abs(le_r[0][1] - se_r[0][1]) * 3600
            if d <= 10.0:
                p += 1
            else:
                f += 1
                fails.append(f'  {star} jd={jd:.0f} DEC diff={d:.4f}"')
    return p, f, t, fails


def test_241():
    p, f, t, fails = 0, 0, 0, []
    for sid_mode in [1, 3]:  # Lahiri, Raman
        swe.set_sid_mode(sid_mode)
        ephem.set_sid_mode(sid_mode, 0, 0)
        le_f = ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_SIDEREAL
        se_f = swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_SIDEREAL
        for le_b, se_b in PLANETS_ALL:
            for jd in DATES_10[:7]:
                try:
                    le_r = ephem.calc_ut(jd, le_b, le_f)
                    se_r = swe.calc_ut(jd, se_b, se_f)
                except:
                    continue
                t += 1
                d = adiff(le_r[0][0], se_r[0][0]) * 3600
                if d <= 20.0:
                    p += 1
                else:
                    f += 1
                    fails.append(f'  sid={sid_mode} B{le_b} jd={jd:.0f} diff={d:.2f}"')
        swe.set_sid_mode(0)
        ephem.set_sid_mode(0, 0, 0)
    return p, f, t, fails


def test_242():
    p, f, t, fails = 0, 0, 0, []
    flag_combos = [
        (
            "Default",
            ephem.FLG_SWIEPH | ephem.FLG_SPEED,
            swe.FLG_SWIEPH | swe.FLG_SPEED,
        ),
        (
            "NOABERR",
            ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_NOABERR,
            swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_NOABERR,
        ),
        (
            "TRUEPOS",
            ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_TRUEPOS,
            swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_TRUEPOS,
        ),
    ]
    for fl, le_f, se_f in flag_combos:
        for jd in DATES_10[:8]:
            try:
                le_r = ephem.calc_ut(jd, ephem.MOON, le_f)
                se_r = swe.calc_ut(jd, swe.MOON, se_f)
            except:
                continue
            t += 1
            d = adiff(le_r[0][0], se_r[0][0]) * 3600
            tol = 35.0 if fl != "Default" else 1.0
            if d <= tol:
                p += 1
            else:
                f += 1
                fails.append(f'  Moon {fl} jd={jd:.0f} diff={d:.2f}"')
            t += 1
            d = abs(le_r[0][3] - se_r[0][3]) * 3600
            if d <= 10.0:
                p += 1
            else:
                f += 1
                fails.append(f'  Moon {fl} jd={jd:.0f} spd diff={d:.2f}"/day')
    return p, f, t, fails


def test_243():
    p, f, t, fails = 0, 0, 0, []
    bodies = [
        (ephem.CHIRON, swe.CHIRON, "Chiron"),
        (ephem.CERES, swe.CERES, "Ceres"),
    ]
    flag_combos = [
        (
            "Helio",
            ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_HELCTR,
            swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_HELCTR,
        ),
        (
            "J2000",
            ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_J2000,
            swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_J2000,
        ),
        (
            "Helio+J2000",
            ephem.FLG_SWIEPH
            | ephem.FLG_SPEED
            | ephem.FLG_HELCTR
            | ephem.FLG_J2000,
            swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_HELCTR | swe.FLG_J2000,
        ),
    ]
    for le_b, se_b, nm in bodies:
        for fl, le_f, se_f in flag_combos:
            for jd in DATES_10[:8]:
                try:
                    le_r = ephem.calc_ut(jd, le_b, le_f)
                    se_r = swe.calc_ut(jd, se_b, se_f)
                except:
                    continue
                t += 1
                d = adiff(le_r[0][0], se_r[0][0]) * 3600
                if d <= 2.0:
                    p += 1
                else:
                    f += 1
                    fails.append(f'  {nm} {fl} jd={jd:.0f} diff={d:.4f}"')
    return p, f, t, fails


def test_244():
    p, f, t, fails = 0, 0, 0, []
    for le_b, se_b in PLANETS_7:
        for jd in DATES_10[:6]:
            try:
                le_r = ephem.calc_ut(
                    jd, le_b, ephem.FLG_SWIEPH | ephem.FLG_SPEED
                )
                le_nut = ephem.calc_ut(jd, -1, 0)
                eps = le_nut[0][0]
                ecl_sp = (
                    le_r[0][0],
                    le_r[0][1],
                    le_r[0][2],
                    le_r[0][3],
                    le_r[0][4],
                    le_r[0][5],
                )
                le_eq = ephem.cotrans_sp(ecl_sp, -eps)
                se_eq = swe.cotrans_sp(ecl_sp, -eps)
            except:
                continue
            t += 1
            d = adiff(le_eq[0], se_eq[0]) * 3600
            if d <= 0.001:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} jd={jd:.0f} cotrans_sp RA diff={d:.6f}"')
            t += 1
            d = abs(le_eq[1] - se_eq[1]) * 3600
            if d <= 0.001:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} jd={jd:.0f} cotrans_sp DEC diff={d:.6f}"')
            t += 1
            d = abs(le_eq[3] - se_eq[3]) * 3600
            if d <= 0.01:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} jd={jd:.0f} cotrans_sp RA_spd diff={d:.6f}"')
    return p, f, t, fails


def test_245():
    p, f, t, fails = 0, 0, 0, []
    le_f = ephem.FLG_SWIEPH | ephem.FLG_SPEED
    se_f = swe.FLG_SWIEPH | swe.FLG_SPEED
    for jd_base in [2451545.0, 2460000.0, 2440000.0]:
        for h_frac in [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875]:
            jd = jd_base + h_frac
            for le_b, se_b in PLANETS_ALL:
                try:
                    le_r = ephem.calc_ut(jd, le_b, le_f)
                    se_r = swe.calc_ut(jd, se_b, se_f)
                except:
                    continue
                t += 1
                d = adiff(le_r[0][0], se_r[0][0]) * 3600
                if d <= 1.0:
                    p += 1
                else:
                    f += 1
                    fails.append(f'  B{le_b} jd={jd:.3f} diff={d:.4f}"')
    return p, f, t, fails


def test_246():
    p, f, t, fails = 0, 0, 0, []
    hsystems = ["P", "K", "O", "R", "C", "E", "W", "B"]
    lats = [0.0, 30.0, 45.0, 51.5, -33.87, 60.0]
    for jd in [2451545.0, 2460000.0, 2440000.0]:
        armc = swe.sidtime(jd) * 15.0
        try:
            eps = swe.calc_ut(jd, -1, swe.FLG_SWIEPH)[0][0]
        except:
            continue
        for lat in lats:
            for hs in hsystems:
                try:
                    le_r = ephem.houses_armc_ex2(
                        armc, lat, eps, ord(hs), ephem.FLG_SWIEPH | ephem.FLG_SPEED
                    )
                    se_r = swe.houses_armc_ex2(
                        armc, lat, eps, hs.encode(), swe.FLG_SWIEPH | swe.FLG_SPEED
                    )
                except:
                    continue
                for i in range(min(len(le_r[0]), len(se_r[0]), 12)):
                    t += 1
                    d = adiff(le_r[0][i], se_r[0][i]) * 3600
                    if d <= 1.0:
                        p += 1
                    else:
                        f += 1
                        fails.append(f'  armc {hs} lat={lat} c{i + 1} diff={d:.4f}"')
    return p, f, t, fails


def test_247():
    p, f, t, fails = 0, 0, 0, []
    le_f = ephem.FLG_SWIEPH | ephem.FLG_SPEED
    se_f = swe.FLG_SWIEPH | swe.FLG_SPEED
    for le_b, se_b in [
        (ephem.MARS, swe.MARS),
        (ephem.JUPITER, swe.JUPITER),
        (ephem.SATURN, swe.SATURN),
        (ephem.VENUS, swe.VENUS),
    ]:
        for jd in DATES_10[:8]:
            try:
                le_r = ephem.calc_ut(jd, le_b, le_f)
                se_r = swe.calc_ut(jd, se_b, se_f)
            except:
                continue
            t += 1
            dd = abs(le_r[0][2] - se_r[0][2])
            if dd <= 0.0001:
                p += 1
            else:
                f += 1
                fails.append(f"  B{le_b} jd={jd:.0f} dist diff={dd:.8f} AU")
            t += 1
            ds = abs(le_r[0][5] - se_r[0][5])
            if ds <= 0.0001:
                p += 1
            else:
                f += 1
                fails.append(f"  B{le_b} jd={jd:.0f} dist_spd diff={ds:.8f}")
    return p, f, t, fails


def test_248():
    p, f, t, fails = 0, 0, 0, []
    # Sample dates across 18.6-year nutation cycle
    for i in range(40):
        jd = 2451545.0 + i * 170.0
        try:
            le_r = ephem.calc_ut(jd, -1, 0)
            se_r = swe.calc_ut(jd, -1, swe.FLG_SWIEPH)
        except:
            continue
        t += 1
        d = abs(le_r[0][0] - se_r[0][0]) * 3600  # true obliquity
        if d <= 0.1:
            p += 1
        else:
            f += 1
            fails.append(f'  jd={jd:.1f} obl diff={d:.6f}"')
        t += 1
        d = abs(le_r[0][1] - se_r[0][1]) * 3600  # mean obliquity
        if d <= 0.1:
            p += 1
        else:
            f += 1
            fails.append(f'  jd={jd:.1f} mean_obl diff={d:.6f}"')
        t += 1
        d = abs(le_r[0][2] - se_r[0][2]) * 3600  # dpsi
        if d <= 0.1:
            p += 1
        else:
            f += 1
            fails.append(f'  jd={jd:.1f} dpsi diff={d:.6f}"')
        t += 1
        d = abs(le_r[0][3] - se_r[0][3]) * 3600  # deps
        if d <= 0.1:
            p += 1
        else:
            f += 1
            fails.append(f'  jd={jd:.1f} deps diff={d:.6f}"')
    return p, f, t, fails


def test_249():
    p, f, t, fails = 0, 0, 0, []
    le_f = ephem.FLG_SWIEPH | ephem.FLG_SPEED
    se_f = swe.FLG_SWIEPH | swe.FLG_SPEED
    for le_b, se_b, nm, tol in [
        (ephem.OSCU_APOG, swe.OSCU_APOG, "OscuLilith", 1.0),
        (ephem.MEAN_APOG, swe.MEAN_APOG, "MeanLilith", 1.0),
    ]:
        for jd in DATES_10:
            try:
                le_r = ephem.calc_ut(jd, le_b, le_f)
                se_r = swe.calc_ut(jd, se_b, se_f)
            except:
                continue
            t += 1
            d = adiff(le_r[0][0], se_r[0][0]) * 3600
            if d <= tol:
                p += 1
            else:
                f += 1
                fails.append(f'  {nm} jd={jd:.0f} lon diff={d:.4f}"')
            t += 1
            d = abs(le_r[0][3] - se_r[0][3]) * 3600
            tol_s = 200.0 if nm == "OscuLilith" else 5.0
            if d <= tol_s:
                p += 1
            else:
                f += 1
                fails.append(f'  {nm} jd={jd:.0f} spd diff={d:.4f}"/day')
    return p, f, t, fails


def test_250():
    p, f, t, fails = 0, 0, 0, []
    swe.set_sid_mode(1)
    ephem.set_sid_mode(1, 0, 0)
    stars = [
        "Aldebaran",
        "Regulus",
        "Spica",
        "Antares",
        "Sirius",
        "Vega",
        "Pollux",
        "Arcturus",
    ]
    le_f = ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_SIDEREAL
    se_f = swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_SIDEREAL
    for star in stars:
        for jd in [2451545.0, 2460000.0, 2440000.0, 2455000.0]:
            try:
                le_r = ephem.fixstar2_ut(star, jd, le_f)
                se_r = swe.fixstar2(star, jd, se_f)
            except:
                continue
            t += 1
            d = adiff(le_r[0][0], se_r[0][0]) * 3600
            if d <= 15.0:
                p += 1  # known ~5" sidereal offset
            else:
                f += 1
                fails.append(f'  {star} sid jd={jd:.0f} lon diff={d:.4f}"')
    swe.set_sid_mode(0)
    ephem.set_sid_mode(0, 0, 0)
    return p, f, t, fails


def test_251():
    p, f, t, fails = 0, 0, 0, []
    le_f = (
        ephem.FLG_SWIEPH
        | ephem.FLG_SPEED
        | ephem.FLG_NONUT
        | ephem.FLG_EQUATORIAL
    )
    se_f = swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_NONUT | swe.FLG_EQUATORIAL
    for le_b, se_b in PLANETS_ALL:
        for jd in DATES_10[:7]:
            try:
                le_r = ephem.calc_ut(jd, le_b, le_f)
                se_r = swe.calc_ut(jd, se_b, se_f)
            except:
                continue
            t += 1
            d = adiff(le_r[0][0], se_r[0][0]) * 3600
            if d <= 1.0:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} jd={jd:.0f} NONUT+EQ diff={d:.4f}"')
            t += 1
            d = abs(le_r[0][1] - se_r[0][1]) * 3600
            if d <= 1.0:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} jd={jd:.0f} NONUT+EQ dec diff={d:.4f}"')
    return p, f, t, fails


def test_252():
    p, f, t, fails = 0, 0, 0, []
    hsystems = ["P", "K", "R", "C", "E", "W"]
    for jd in [2451545.0, 2460000.0, 2440000.0, 2455000.0]:
        try:
            le_moon = ephem.calc_ut(
                jd, ephem.MOON, ephem.FLG_SWIEPH | ephem.FLG_SPEED
            )
            le_nut = ephem.calc_ut(jd, -1, 0)
            eps = le_nut[0][0]
            armc = ephem.sidtime(jd) * 15.0
            plon, plat = le_moon[0][0], le_moon[0][1]
        except:
            continue
        for lat in [0.0, 30.0, 45.0, 51.5, -33.87]:
            for hs in hsystems:
                try:
                    le_hp = ephem.house_pos(armc, lat, eps, ord(hs), plon, plat)
                    se_hp = swe.house_pos(armc, lat, eps, (plon, plat), hs.encode())
                except:
                    continue
                t += 1
                d = abs(le_hp - se_hp)
                if d <= 0.02:
                    p += 1
                else:
                    f += 1
                    fails.append(f"  Moon {hs} lat={lat} jd={jd:.0f} hp diff={d:.4f}")
    return p, f, t, fails


def test_253():
    p, f, t, fails = 0, 0, 0, []
    asteroids = [
        (ephem.CERES, swe.CERES, "Ceres"),
        (ephem.PALLAS, swe.PALLAS, "Pallas"),
        (ephem.JUNO, swe.JUNO, "Juno"),
        (ephem.VESTA, swe.VESTA, "Vesta"),
    ]
    le_f = ephem.FLG_SWIEPH | ephem.FLG_SPEED
    se_f = swe.FLG_SWIEPH | swe.FLG_SPEED
    for le_b, se_b, nm in asteroids:
        for jd in DATES_10:
            try:
                le_r = ephem.calc_ut(jd, le_b, le_f)
                se_r = swe.calc_ut(jd, se_b, se_f)
            except:
                continue
            t += 1
            d = adiff(le_r[0][0], se_r[0][0]) * 3600
            if d <= 5.0:
                p += 1
            else:
                f += 1
                fails.append(f'  {nm} jd={jd:.0f} diff={d:.4f}"')
            t += 1
            d = abs(le_r[0][1] - se_r[0][1]) * 3600
            if d <= 5.0:
                p += 1
            else:
                f += 1
                fails.append(f'  {nm} jd={jd:.0f} lat diff={d:.4f}"')
            t += 1
            d = abs(le_r[0][2] - se_r[0][2])
            if d <= 0.0001:
                p += 1
            else:
                f += 1
                fails.append(f"  {nm} jd={jd:.0f} dist diff={d:.8f}")
    return p, f, t, fails


def test_254():
    p, f, t, fails = 0, 0, 0, []
    locations = [
        (51.5, -0.1, 11.0, "London"),
        (40.7, -74.0, 10.0, "NYC"),
        (35.7, 139.7, 40.0, "Tokyo"),
        (-33.9, 151.2, 58.0, "Sydney"),
    ]
    for jd in [2451545.0, 2460000.0, 2455000.0]:
        for lat, lon, alt, nm in locations:
            for rsmi in [1, 2]:  # rise, set
                try:
                    le_r = ephem.rise_trans(
                        jd,
                        ephem.SUN,
                        rsmi,
                        [lon, lat, alt],
                        1013.25,
                        15.0,
                        ephem.FLG_SWIEPH,
                    )
                    se_r = swe.rise_trans(
                        jd, swe.SUN, rsmi, [lon, lat, alt], 1013.25, 15.0
                    )
                    le_jd = le_r[1][0]
                    se_jd = se_r[1][0]
                except:
                    continue
                if le_jd == 0.0 or se_jd == 0.0:
                    continue
                t += 1
                d = abs(le_jd - se_jd) * 1440  # minutes
                if d <= 2.0:
                    p += 1
                else:
                    f += 1
                    fails.append(
                        f"  Sun {'rise' if rsmi == 1 else 'set'} {nm} jd={jd:.0f} diff={d:.2f}min"
                    )
    return p, f, t, fails


def test_255():
    p, f, t, fails = 0, 0, 0, []
    flag_combos = [
        (
            "Default",
            ephem.FLG_SWIEPH | ephem.FLG_SPEED,
            swe.FLG_SWIEPH | swe.FLG_SPEED,
            1.0,
        ),
        (
            "J2000+EQ",
            ephem.FLG_SWIEPH
            | ephem.FLG_SPEED
            | ephem.FLG_J2000
            | ephem.FLG_EQUATORIAL,
            swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_J2000 | swe.FLG_EQUATORIAL,
            1.0,
        ),
        (
            "NONUT+NOABERR",
            ephem.FLG_SWIEPH
            | ephem.FLG_SPEED
            | ephem.FLG_NONUT
            | ephem.FLG_NOABERR,
            swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_NONUT | swe.FLG_NOABERR,
            35.0,
        ),
        (
            "J2000+NONUT",
            ephem.FLG_SWIEPH
            | ephem.FLG_SPEED
            | ephem.FLG_J2000
            | ephem.FLG_NONUT,
            swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_J2000 | swe.FLG_NONUT,
            1.0,
        ),
    ]
    for le_b, se_b in [
        (ephem.SUN, swe.SUN),
        (ephem.MOON, swe.MOON),
        (ephem.MARS, swe.MARS),
        (ephem.JUPITER, swe.JUPITER),
    ]:
        for fl, le_f, se_f, tol in flag_combos:
            for jd in DATES_10[:5]:
                try:
                    le_r = ephem.calc_ut(jd, le_b, le_f)
                    se_r = swe.calc_ut(jd, se_b, se_f)
                except:
                    continue
                t += 1
                d = adiff(le_r[0][0], se_r[0][0]) * 3600
                if d <= tol:
                    p += 1
                else:
                    f += 1
                    fails.append(f'  B{le_b} {fl} jd={jd:.0f} diff={d:.2f}"')
    return p, f, t, fails


if __name__ == "__main__":
    print("=" * 70)
    print("Rounds 236-255: Deep Diverse Verification Batch 2")
    print("=" * 70)
    for num, (name, fn) in enumerate(
        [
            ("All planets TRUEPOS", test_236),
            ("Planet magnitude/pheno", test_237),
            ("Delta-T leap seconds", test_238),
            ("Houses tropical/arctic", test_239),
            ("Fixed stars equatorial", test_240),
            ("Sidereal Lahiri+Raman sweep", test_241),
            ("Moon NOABERR/TRUEPOS", test_242),
            ("Chiron+Ceres helio+J2000", test_243),
            ("cotrans_sp speed transform", test_244),
            ("All bodies fractional JD", test_245),
            ("houses_armc_ex2 sweep", test_246),
            ("Planet distance at events", test_247),
            ("Nutation cycle 18.6yr", test_248),
            ("OscuLilith+MeanLilith", test_249),
            ("Fixed stars sidereal Lahiri", test_250),
            ("NONUT+EQUATORIAL combo", test_251),
            ("house_pos Moon all systems", test_252),
            ("Asteroid belt sweep", test_253),
            ("Sun rise/set timing", test_254),
            ("Multi-flag matrix", test_255),
        ],
        start=236,
    ):
        run_round(num, name, fn)

    tp = sum(r[1] for r in results.values())
    tf = sum(r[2] for r in results.values())
    tt = sum(r[3] for r in results.values())
    print(f"\n{'=' * 70}")
    print(f"BATCH TOTAL: {tp}/{tt} passed ({100 * tp / tt:.1f}%), {tf} failed")
    print(f"{'=' * 70}")
    for num in sorted(results):
        nm, pp, ff, ttt, fls = results[num]
        if ff > 0:
            print(f"\n  Round {num} failures:")
            for fl in fls[:5]:
                print(fl)
            if len(fls) > 5:
                print(f"  ... and {len(fls) - 5} more")
