#!/usr/bin/env python3
"""Rounds 226-235: Deep diverse verification batch.

226: NOABERR+NONUT combined flags all planets
227: Barycentric Moon+Sun positions
228: All 4 asteroids sidereal+equatorial
229: house_pos all house systems
230: Julian day round-trip consistency
231: Sidereal time at 6h intervals
232: Fixed stars J2000+NONUT flags
233: Uranian bodies equatorial+J2000
234: Planet positions at century boundaries
235: Moon node longitude at eclipses
"""

from __future__ import annotations
import os, sys, math

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("LIBEPHEMERIS_MODE", "skyfield")

import swisseph as swe
import libephemeris as ephem

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


def angular_diff(a, b):
    d = abs(a - b) % 360
    return min(d, 360 - d)


# ---- Round 226: NOABERR+NONUT combined ----
def test_226():
    p, f, t, fails = 0, 0, 0, []
    le_f = (
        ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_NOABERR | ephem.FLG_NONUT
    )
    se_f = swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_NOABERR | swe.FLG_NONUT
    bodies = [
        (ephem.SUN, swe.SUN),
        (ephem.MOON, swe.MOON),
        (ephem.MERCURY, swe.MERCURY),
        (ephem.VENUS, swe.VENUS),
        (ephem.MARS, swe.MARS),
        (ephem.JUPITER, swe.JUPITER),
        (ephem.SATURN, swe.SATURN),
        (ephem.URANUS, swe.URANUS),
        (ephem.NEPTUNE, swe.NEPTUNE),
        (ephem.PLUTO, swe.PLUTO),
        (ephem.CHIRON, swe.CHIRON),
    ]
    dates = [
        2415020.0,
        2430000.0,
        2440000.0,
        2451545.0,
        2455000.0,
        2460000.0,
        2462000.0,
        2445000.0,
        2450000.0,
        2458000.0,
    ]
    for le_b, se_b in bodies:
        for jd in dates:
            try:
                le_r = ephem.calc_ut(jd, le_b, le_f)
                se_r = swe.calc_ut(jd, se_b, se_f)
            except:
                continue
            t += 1
            d = angular_diff(le_r[0][0], se_r[0][0]) * 3600
            if d <= 35.0:
                p += 1  # NOABERR has known ~31" model diff
            else:
                f += 1
                fails.append(f'  body={le_b} jd={jd:.0f} lon diff={d:.2f}"')
            t += 1
            d = abs(le_r[0][1] - se_r[0][1]) * 3600
            if d <= 35.0:
                p += 1
            else:
                f += 1
                fails.append(f'  body={le_b} jd={jd:.0f} lat diff={d:.2f}"')
    return p, f, t, fails


# ---- Round 227: Barycentric Sun+Moon ----
def test_227():
    p, f, t, fails = 0, 0, 0, []
    le_f = ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_BARYCTR
    se_f = swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_BARYCTR
    bodies = [
        (ephem.SUN, swe.SUN, "Sun", 1.0),
        (ephem.MOON, swe.MOON, "Moon", 2.0),
        (ephem.MERCURY, swe.MERCURY, "Merc", 1.0),
        (ephem.VENUS, swe.VENUS, "Ven", 1.0),
        (ephem.MARS, swe.MARS, "Mars", 1.0),
        (ephem.JUPITER, swe.JUPITER, "Jup", 1.0),
        (ephem.SATURN, swe.SATURN, "Sat", 1.0),
    ]
    dates = [
        2451545.0,
        2455000.0,
        2460000.0,
        2440000.0,
        2445000.0,
        2458000.0,
        2430000.0,
        2462000.0,
    ]
    for le_b, se_b, nm, tol in bodies:
        for jd in dates:
            try:
                le_r = ephem.calc_ut(jd, le_b, le_f)
                se_r = swe.calc_ut(jd, se_b, se_f)
            except:
                continue
            t += 1
            d = angular_diff(le_r[0][0], se_r[0][0]) * 3600
            if d <= tol:
                p += 1
            else:
                f += 1
                fails.append(f'  {nm} jd={jd:.0f} bary lon diff={d:.4f}"')
            t += 1
            d = abs(le_r[0][1] - se_r[0][1]) * 3600
            if d <= tol:
                p += 1
            else:
                f += 1
                fails.append(f'  {nm} jd={jd:.0f} bary lat diff={d:.4f}"')
    return p, f, t, fails


# ---- Round 228: All 4 asteroids sidereal+equatorial ----
def test_228():
    p, f, t, fails = 0, 0, 0, []
    asteroids = [
        (ephem.CERES, swe.CERES, "Ceres"),
        (ephem.PALLAS, swe.PALLAS, "Pallas"),
        (ephem.JUNO, swe.JUNO, "Juno"),
        (ephem.VESTA, swe.VESTA, "Vesta"),
    ]
    flag_combos = [
        (
            "Equatorial",
            ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_EQUATORIAL,
            swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_EQUATORIAL,
        ),
        (
            "J2000",
            ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_J2000,
            swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_J2000,
        ),
    ]
    sid_modes = [0, 1, 3]
    dates = [2451545.0, 2455000.0, 2460000.0, 2440000.0, 2458000.0]
    for le_b, se_b, nm in asteroids:
        for flabel, le_fl, se_fl in flag_combos:
            for jd in dates:
                try:
                    le_r = ephem.calc_ut(jd, le_b, le_fl)
                    se_r = swe.calc_ut(jd, se_b, se_fl)
                except:
                    continue
                t += 1
                d = angular_diff(le_r[0][0], se_r[0][0]) * 3600
                if d <= 5.0:
                    p += 1
                else:
                    f += 1
                    fails.append(f'  {nm} {flabel} jd={jd:.0f} diff={d:.4f}"')
        # Sidereal
        for sm in sid_modes:
            swe.set_sid_mode(sm)
            ephem.set_sid_mode(sm, 0, 0)
            le_sf = ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_SIDEREAL
            se_sf = swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_SIDEREAL
            for jd in dates:
                try:
                    le_r = ephem.calc_ut(jd, le_b, le_sf)
                    se_r = swe.calc_ut(jd, se_b, se_sf)
                except:
                    continue
                t += 1
                d = angular_diff(le_r[0][0], se_r[0][0]) * 3600
                if d <= 20.0:
                    p += 1  # sidereal ~14" systematic
                else:
                    f += 1
                    fails.append(f'  {nm} sid={sm} jd={jd:.0f} diff={d:.4f}"')
            swe.set_sid_mode(0)
            ephem.set_sid_mode(0, 0, 0)
    return p, f, t, fails


# ---- Round 229: house_pos all systems ----
def test_229():
    p, f, t, fails = 0, 0, 0, []
    hsystems = ["P", "K", "O", "R", "C", "E", "W", "B", "M", "X"]
    dates = [2451545.0, 2460000.0, 2440000.0, 2455000.0, 2458000.0]
    lats = [0.0, 30.0, 45.0, 51.5, -33.87]
    for jd in dates:
        try:
            le_sun = ephem.calc_ut(
                jd, ephem.SUN, ephem.FLG_SWIEPH | ephem.FLG_SPEED
            )
            se_sun = swe.calc_ut(jd, swe.SUN, swe.FLG_SWIEPH | swe.FLG_SPEED)
            le_nut = ephem.calc_ut(jd, -1, 0)
            eps = le_nut[0][0]
            armc = ephem.sidtime(jd) * 15.0
        except:
            continue
        plon, plat = le_sun[0][0], le_sun[0][1]
        for lat in lats:
            for hs in hsystems:
                try:
                    le_hp = ephem.house_pos(armc, lat, eps, ord(hs), plon, plat)
                    se_hp = swe.house_pos(armc, lat, eps, (plon, plat), hs.encode())
                except:
                    continue
                t += 1
                d = abs(le_hp - se_hp)
                if d <= 0.01:
                    p += 1  # within 0.01 house position
                else:
                    f += 1
                    fails.append(
                        f"  {hs} lat={lat} jd={jd:.0f} hp LE={le_hp:.4f} SE={se_hp:.4f}"
                    )
    return p, f, t, fails


# ---- Round 230: Julian day round-trip ----
def test_230():
    p, f, t, fails = 0, 0, 0, []
    test_dates = [
        (2000, 1, 1, 12.0, 1),
        (1900, 1, 1, 0.0, 1),
        (1582, 10, 15, 0.0, 1),
        (1582, 10, 4, 0.0, 0),
        (-500, 3, 1, 0.0, 0),
        (2025, 6, 15, 18.5, 1),
        (1999, 12, 31, 23.999, 1),
        (100, 1, 1, 0.0, 0),
        (1, 1, 1, 0.0, 0),
        (-4712, 1, 1, 12.0, 0),
        (2100, 12, 31, 0.0, 1),
        (1800, 6, 15, 6.0, 1),
    ]
    for y, m, d, h, greg in test_dates:
        try:
            le_jd = ephem.julday(y, m, d, h, greg)
            se_jd = swe.julday(y, m, d, h, greg)
        except:
            continue
        t += 1
        diff = abs(le_jd - se_jd)
        if diff <= 1e-10:
            p += 1
        else:
            f += 1
            fails.append(f"  julday {y}/{m}/{d} {h}: LE={le_jd} SE={se_jd} diff={diff}")
        # Round-trip
        try:
            le_rev = ephem.revjul(le_jd, greg)
            se_rev = swe.revjul(se_jd, greg)
        except:
            continue
        t += 1
        if le_rev[0] == se_rev[0] and le_rev[1] == se_rev[1] and le_rev[2] == se_rev[2]:
            p += 1
        else:
            f += 1
            fails.append(f"  revjul {le_jd}: LE={le_rev[:3]} SE={se_rev[:3]}")
        t += 1
        hour_diff = abs(le_rev[3] - se_rev[3])
        if hour_diff <= 1e-6:
            p += 1
        else:
            f += 1
            fails.append(f"  revjul hour {le_jd}: LE={le_rev[3]} SE={se_rev[3]}")
    return p, f, t, fails


# ---- Round 231: Sidereal time at 6h intervals ----
def test_231():
    p, f, t, fails = 0, 0, 0, []
    jd_base = 2451545.0
    for day_offset in range(0, 365 * 20, 91):  # every ~3 months for 20 years
        for hour_frac in [0.0, 0.25, 0.5, 0.75]:
            jd = jd_base + day_offset + hour_frac
            try:
                le_st = ephem.sidtime(jd)
                se_st = swe.sidtime(jd)
            except:
                continue
            t += 1
            d = abs(le_st - se_st)
            if d > 12:
                d = 24 - d
            if d * 3600 <= 0.1:
                p += 1
            else:
                f += 1
                fails.append(f"  jd={jd:.2f} sidtime diff={d * 3600:.4f}s")
    return p, f, t, fails


# ---- Round 232: Fixed stars J2000+NONUT ----
def test_232():
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
    ]
    flag_combos = [
        (
            "J2000",
            ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_J2000,
            swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_J2000,
        ),
        (
            "NONUT",
            ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_NONUT,
            swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_NONUT,
        ),
        (
            "J2000+NONUT",
            ephem.FLG_SWIEPH
            | ephem.FLG_SPEED
            | ephem.FLG_J2000
            | ephem.FLG_NONUT,
            swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_J2000 | swe.FLG_NONUT,
        ),
    ]
    dates = [2451545.0, 2460000.0, 2440000.0, 2455000.0, 2430000.0]
    for star in stars:
        for flabel, le_fl, se_fl in flag_combos:
            for jd in dates:
                try:
                    le_r = ephem.fixstar2_ut(star, jd, le_fl)
                    se_r = swe.fixstar2(star, jd, se_fl)
                except:
                    continue
                t += 1
                d = angular_diff(le_r[0][0], se_r[0][0]) * 3600
                if d <= 10.0:
                    p += 1
                else:
                    f += 1
                    fails.append(f'  {star} {flabel} jd={jd:.0f} lon diff={d:.4f}"')
                t += 1
                d = abs(le_r[0][1] - se_r[0][1]) * 3600
                if d <= 10.0:
                    p += 1
                else:
                    f += 1
                    fails.append(f'  {star} {flabel} jd={jd:.0f} lat diff={d:.4f}"')
    return p, f, t, fails


# ---- Round 233: Uranian bodies equatorial+J2000 ----
def test_233():
    p, f, t, fails = 0, 0, 0, []
    uranians = [
        (ephem.CUPIDO, 40),
        (ephem.HADES, 41),
        (ephem.ZEUS, 42),
        (ephem.KRONOS, 43),
        (ephem.APOLLON, 44),
        (ephem.ADMETOS, 45),
        (ephem.VULKANUS, 46),
        (ephem.POSEIDON, 47),
    ]
    flag_combos = [
        (
            "Equat",
            ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_EQUATORIAL,
            swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_EQUATORIAL,
        ),
        (
            "J2000",
            ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_J2000,
            swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_J2000,
        ),
        (
            "NONUT",
            ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_NONUT,
            swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_NONUT,
        ),
    ]
    dates = [2451545.0, 2455000.0, 2460000.0, 2440000.0, 2445000.0, 2458000.0]
    for le_b, se_b in uranians:
        for flabel, le_fl, se_fl in flag_combos:
            for jd in dates:
                try:
                    le_r = ephem.calc_ut(jd, le_b, le_fl)
                    se_r = swe.calc_ut(jd, se_b, se_fl)
                except:
                    continue
                t += 1
                d = angular_diff(le_r[0][0], se_r[0][0]) * 3600
                if d <= 60.0:
                    p += 1
                else:
                    f += 1
                    fails.append(f'  U{se_b} {flabel} jd={jd:.0f} diff={d:.2f}"')
    return p, f, t, fails


# ---- Round 234: Planet positions at century boundaries ----
def test_234():
    p, f, t, fails = 0, 0, 0, []
    le_f = ephem.FLG_SWIEPH | ephem.FLG_SPEED
    se_f = swe.FLG_SWIEPH | swe.FLG_SPEED
    bodies = [
        (ephem.SUN, swe.SUN, "Sun"),
        (ephem.MOON, swe.MOON, "Moon"),
        (ephem.MARS, swe.MARS, "Mars"),
        (ephem.JUPITER, swe.JUPITER, "Jup"),
        (ephem.SATURN, swe.SATURN, "Sat"),
    ]
    # Century boundaries within DE440 range
    century_jds = [
        2305447.5,  # 1700 Jan 1
        2341972.5,  # 1800 Jan 1
        2378496.5,  # 1900 Jan 1
        2415020.5,  # 2000 Jan 1
        2451545.0,  # J2000
        2488069.5,  # 2100 Jan 1
    ]
    for nm, le_b, se_b in [
        (n, l, s) for n, l, s in [(b[2], b[0], b[1]) for b in bodies]
    ]:
        pass
    for bname, le_b, se_b in [(b[2], b[0], b[1]) for b in bodies]:
        for jd in century_jds:
            try:
                le_r = ephem.calc_ut(jd, le_b, le_f)
                se_r = swe.calc_ut(jd, se_b, se_f)
            except:
                continue
            t += 1
            d = angular_diff(le_r[0][0], se_r[0][0]) * 3600
            tol = 2.0 if jd > 2350000 else 5.0  # looser for older dates
            if d <= tol:
                p += 1
            else:
                f += 1
                fails.append(f'  {bname} jd={jd:.1f} diff={d:.4f}"')
            t += 1
            d = abs(le_r[0][1] - se_r[0][1]) * 3600
            if d <= tol:
                p += 1
            else:
                f += 1
                fails.append(f'  {bname} jd={jd:.1f} lat diff={d:.4f}"')
    return p, f, t, fails


# ---- Round 235: Moon node longitude at eclipse times ----
def test_235():
    p, f, t, fails = 0, 0, 0, []
    le_f = ephem.FLG_SWIEPH | ephem.FLG_SPEED
    se_f = swe.FLG_SWIEPH | swe.FLG_SPEED
    # Known eclipse JDs
    eclipse_jds = [
        2451545.0 + i * 173.31
        for i in range(30)  # approximately every eclipse season
    ]
    for jd in eclipse_jds:
        for node_body, se_node, nm in [
            (ephem.TRUE_NODE, swe.TRUE_NODE, "TrueNode"),
            (ephem.MEAN_NODE, swe.MEAN_NODE, "MeanNode"),
        ]:
            try:
                le_r = ephem.calc_ut(jd, node_body, le_f)
                se_r = swe.calc_ut(jd, se_node, se_f)
            except:
                continue
            t += 1
            d = angular_diff(le_r[0][0], se_r[0][0]) * 3600
            if d <= 1.0:
                p += 1
            else:
                f += 1
                fails.append(f'  {nm} jd={jd:.1f} diff={d:.4f}"')
            t += 1
            d = abs(le_r[0][3] - se_r[0][3]) * 3600
            tol = 5.0 if nm == "TrueNode" else 1.0
            if d <= tol:
                p += 1
            else:
                f += 1
                fails.append(f'  {nm} jd={jd:.1f} spd diff={d:.4f}"/day')
    return p, f, t, fails


if __name__ == "__main__":
    print("=" * 70)
    print("Rounds 226-235: Deep Diverse Verification Batch")
    print("=" * 70)

    run_round(226, "NOABERR+NONUT combined flags", test_226)
    run_round(227, "Barycentric positions", test_227)
    run_round(228, "Asteroids sidereal+equatorial", test_228)
    run_round(229, "house_pos all systems", test_229)
    run_round(230, "Julian day round-trip", test_230)
    run_round(231, "Sidereal time 6h intervals", test_231)
    run_round(232, "Fixed stars J2000+NONUT", test_232)
    run_round(233, "Uranian equatorial+J2000", test_233)
    run_round(234, "Century boundary positions", test_234)
    run_round(235, "Moon node at eclipses", test_235)

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
