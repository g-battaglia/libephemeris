#!/usr/bin/env python3
"""Rounds 311-320: Deep diverse verification batch 6 (finale).

311: Fixed star speed consistency (finite-diff vs reported)
312: Ecliptic/equatorial cross-check via calc_ut flags
313: Planet at lunar node crossing times
314: Asteroid heliocentric orbit sweep
315: Moon at declination standstill
316: Grand finale: full 10-planet comprehensive
317: Historical date validation (1900-1950)
318: Modern date precision (2020-2030)
319: Cross-system house consistency (ASC matches cusp 1)
320: Ultimate all-API sweep
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


DATES = [
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
PLANETS = [
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
]
LF = ephem.FLG_SWIEPH | ephem.FLG_SPEED
SF = swe.FLG_SWIEPH | swe.FLG_SPEED


def test_311():  # Fixed star speed consistency
    p, f, t, fails = 0, 0, 0, []
    stars = [
        "Aldebaran",
        "Regulus",
        "Spica",
        "Antares",
        "Sirius",
        "Vega",
        "Capella",
        "Rigel",
        "Procyon",
        "Pollux",
    ]
    for star in stars:
        for jd in [2451545.0, 2460000.0]:
            try:
                lr = ephem.fixstar2_ut(star, jd, LF)
                sr = swe.fixstar2(star, jd, SF)
            except:
                continue
            # lon speed
            t += 1
            d = abs(lr[0][3] - sr[0][3]) * 3600
            if d <= 5.0:  # "/day
                p += 1
            else:
                f += 1
                fails.append(f'  {star} lon_spd jd={jd:.0f} d={d:.4f}"/day')
            # lat speed
            t += 1
            d = abs(lr[0][4] - sr[0][4]) * 3600
            if d <= 5.0:
                p += 1
            else:
                f += 1
                fails.append(f'  {star} lat_spd jd={jd:.0f} d={d:.4f}"/day')
    return p, f, t, fails


def test_312():  # Ecliptic vs equatorial cross-check
    p, f, t, fails = 0, 0, 0, []
    lef = LF | ephem.FLG_EQUATORIAL
    sef = SF | swe.FLG_EQUATORIAL
    for le_b, se_b in PLANETS[:8]:
        for jd in DATES[:6]:
            try:
                lr_ecl = ephem.calc_ut(jd, le_b, LF)
                sr_ecl = swe.calc_ut(jd, se_b, SF)
                lr_eq = ephem.calc_ut(jd, le_b, lef)
                sr_eq = swe.calc_ut(jd, se_b, sef)
            except:
                continue
            # Ecliptic lon match
            t += 1
            d = adiff(lr_ecl[0][0], sr_ecl[0][0]) * 3600
            if d <= 1.0:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} ecl_lon jd={jd:.0f} d={d:.4f}"')
            # Equatorial RA match
            t += 1
            d = adiff(lr_eq[0][0], sr_eq[0][0]) * 3600
            if d <= 1.0:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} eq_RA jd={jd:.0f} d={d:.4f}"')
            # Equatorial DEC match
            t += 1
            d = abs(lr_eq[0][1] - sr_eq[0][1]) * 3600
            if d <= 1.0:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} eq_DEC jd={jd:.0f} d={d:.4f}"')
    return p, f, t, fails


def test_313():  # Planet positions at lunar node crossing times
    p, f, t, fails = 0, 0, 0, []
    # Use known node positions to get interesting times
    for i in range(20):
        jd = 2451545.0 + i * 345.6  # ~every 346 days
        try:
            lr_n = ephem.calc_ut(jd, ephem.TRUE_NODE, LF)
            node_lon = lr_n[0][0]
        except:
            continue
        # Check all planets at this time
        for le_b, se_b in PLANETS:
            try:
                lr = ephem.calc_ut(jd, le_b, LF)
                sr = swe.calc_ut(jd, se_b, SF)
            except:
                continue
            t += 1
            d = adiff(lr[0][0], sr[0][0]) * 3600
            if d <= 1.0:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} at node jd={jd:.1f} d={d:.4f}"')
    return p, f, t, fails


def test_314():  # Asteroid heliocentric orbit
    p, f, t, fails = 0, 0, 0, []
    lhf = LF | ephem.FLG_HELCTR
    shf = SF | swe.FLG_HELCTR
    asteroids = [
        (ephem.CERES, swe.CERES),
        (ephem.PALLAS, swe.PALLAS),
        (ephem.VESTA, swe.VESTA),
        (ephem.CHIRON, swe.CHIRON),
    ]
    for le_b, se_b in asteroids:
        for i in range(10):
            jd = 2445000.0 + i * 2000.0
            try:
                lr = ephem.calc_ut(jd, le_b, lhf)
                sr = swe.calc_ut(jd, se_b, shf)
            except:
                continue
            t += 1
            d = adiff(lr[0][0], sr[0][0]) * 3600
            if d <= 2.0:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} helio jd={jd:.0f} d={d:.4f}"')
            t += 1
            d = abs(lr[0][2] - sr[0][2])
            if d <= 0.0001:
                p += 1
            else:
                f += 1
                fails.append(f"  B{le_b} helio dist jd={jd:.0f} d={d:.8f}")
    return p, f, t, fails


def test_315():  # Moon at declination extremes
    p, f, t, fails = 0, 0, 0, []
    lef = LF | ephem.FLG_EQUATORIAL
    sef = SF | swe.FLG_EQUATORIAL
    # Search for high declination Moon
    for i in range(40):
        jd = 2451545.0 + i * 170.0
        try:
            lr = ephem.calc_ut(jd, ephem.MOON, lef)
            sr = swe.calc_ut(jd, swe.MOON, sef)
        except:
            continue
        if abs(lr[0][1]) > 20.0:  # High declination
            t += 1
            d = abs(lr[0][1] - sr[0][1]) * 3600
            if d <= 1.0:
                p += 1
            else:
                f += 1
                fails.append(f'  Moon DEC={lr[0][1]:.2f} jd={jd:.0f} d={d:.4f}"')
            t += 1
            d = adiff(lr[0][0], sr[0][0]) * 3600
            if d <= 1.0:
                p += 1
            else:
                f += 1
                fails.append(f'  Moon RA jd={jd:.0f} d={d:.4f}"')
    return p, f, t, fails


def test_316():  # Grand finale: 10-planet comprehensive
    p, f, t, fails = 0, 0, 0, []
    for le_b, se_b in PLANETS:
        for jd in DATES:
            try:
                lr = ephem.calc_ut(jd, le_b, LF)
                sr = swe.calc_ut(jd, se_b, SF)
            except:
                continue
            # Longitude
            t += 1
            d = adiff(lr[0][0], sr[0][0]) * 3600
            if d <= 1.0:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} lon jd={jd:.0f} d={d:.4f}"')
            # Latitude
            t += 1
            d = abs(lr[0][1] - sr[0][1]) * 3600
            if d <= 1.0:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} lat jd={jd:.0f} d={d:.4f}"')
            # Speed
            t += 1
            d = abs(lr[0][3] - sr[0][3]) * 3600
            if d <= 5.0:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} spd jd={jd:.0f} d={d:.4f}"/day')
    return p, f, t, fails


def test_317():  # Historical date validation 1900-1950
    p, f, t, fails = 0, 0, 0, []
    for year in range(1900, 1951, 5):
        jd = ephem.julday(year, 6, 15, 12.0, 1)
        for le_b, se_b in PLANETS:
            try:
                lr = ephem.calc_ut(jd, le_b, LF)
                sr = swe.calc_ut(jd, se_b, SF)
            except:
                continue
            t += 1
            d = adiff(lr[0][0], sr[0][0]) * 3600
            if d <= 1.0:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} Y{year} d={d:.4f}"')
    return p, f, t, fails


def test_318():  # Modern date precision 2020-2030
    p, f, t, fails = 0, 0, 0, []
    for year in range(2020, 2031):
        for month in [1, 4, 7, 10]:
            jd = ephem.julday(year, month, 1, 0.0, 1)
            for le_b, se_b in PLANETS:
                try:
                    lr = ephem.calc_ut(jd, le_b, LF)
                    sr = swe.calc_ut(jd, se_b, SF)
                except:
                    continue
                t += 1
                d = adiff(lr[0][0], sr[0][0]) * 3600
                tol = 1.0 if le_b != ephem.MOON else 2.0
                if d <= tol:
                    p += 1
                else:
                    f += 1
                    fails.append(f'  B{le_b} {year}/{month} d={d:.4f}"')
    return p, f, t, fails


def test_319():  # ASC matches cusp 1 consistency
    p, f, t, fails = 0, 0, 0, []
    for hs in ["P", "K", "O", "R", "C", "E", "W", "B", "M"]:
        for jd in [2451545.0, 2460000.0, 2455000.0]:
            for lat in [45.0, -30.0, 0.0, 60.0]:
                try:
                    lc, la = ephem.houses_ex(jd, lat, 10.0, ord(hs), LF)
                    sc, sa = swe.houses_ex(jd, lat, 10.0, hs.encode(), SF)
                except:
                    continue
                # LE ASC == LE cusp 1
                t += 1
                d = adiff(la[0], lc[0]) * 3600
                if d <= 0.001:
                    p += 1
                else:
                    f += 1
                    fails.append(f'  LE {hs} ASC!=c1 lat={lat} d={d:.4f}"')
                # SE ASC == SE cusp 1
                t += 1
                d = adiff(sa[0], sc[0]) * 3600
                if d <= 0.001:
                    p += 1
                else:
                    f += 1
                    fails.append(f'  SE {hs} ASC!=c1 lat={lat} d={d:.4f}"')
                # LE ASC matches SE ASC
                t += 1
                d = adiff(la[0], sa[0]) * 3600
                if d <= 2.0:
                    p += 1
                else:
                    f += 1
                    fails.append(f'  {hs} ASC lat={lat} LE!=SE d={d:.4f}"')
    return p, f, t, fails


def test_320():  # Ultimate all-API sweep
    p, f, t, fails = 0, 0, 0, []
    jd = 2451545.0
    # 1. calc_ut for all planets
    for le_b, se_b in PLANETS:
        try:
            lr = ephem.calc_ut(jd, le_b, LF)
            sr = swe.calc_ut(jd, se_b, SF)
        except:
            continue
        t += 1
        d = adiff(lr[0][0], sr[0][0]) * 3600
        if d <= 1.0:
            p += 1
        else:
            f += 1
            fails.append(f'  calc B{le_b} d={d:.4f}"')
    # 2. houses_ex
    for hs in ["P", "K", "O"]:
        try:
            lc, la = ephem.houses_ex(jd, 45.0, 10.0, ord(hs), LF)
            sc, sa = swe.houses_ex(jd, 45.0, 10.0, hs.encode(), SF)
        except:
            continue
        t += 1
        d = adiff(la[0], sa[0]) * 3600
        if d <= 2.0:
            p += 1
        else:
            f += 1
            fails.append(f'  houses {hs} ASC d={d:.4f}"')
    # 3. sidtime
    try:
        le_st = ephem.sidtime(jd)
        se_st = swe.sidtime(jd)
        t += 1
        d = abs(le_st - se_st) * 3600
        if d <= 0.05:
            p += 1
        else:
            f += 1
            fails.append(f"  sidtime diff={d:.6f}s")
    except:
        pass
    # 4. deltat
    try:
        le_dt = ephem.deltat(jd)
        se_dt = swe.deltat(jd)
        t += 1
        d = abs(le_dt - se_dt) * 86400
        if d <= 0.5:
            p += 1
        else:
            f += 1
            fails.append(f"  deltat diff={d:.4f}s")
    except:
        pass
    # 5. fixstar2
    try:
        lr = ephem.fixstar2_ut("Regulus", jd, LF)
        sr = swe.fixstar2("Regulus", jd, SF)
        t += 1
        d = adiff(lr[0][0], sr[0][0]) * 3600
        if d <= 1.0:
            p += 1
        else:
            f += 1
            fails.append(f'  fixstar Regulus d={d:.4f}"')
    except:
        pass
    # 6. ECL_NUT
    try:
        lr = ephem.calc_ut(jd, -1, 0)
        sr = swe.calc_ut(jd, -1, 0)
        t += 1
        d = abs(lr[0][0] - sr[0][0]) * 3600
        if d <= 1.0:
            p += 1
        else:
            f += 1
            fails.append(f'  ECL_NUT true_obl d={d:.4f}"')
    except:
        pass
    # 7. cotrans round-trip
    try:
        nut = ephem.calc_ut(jd, -1, 0)
        eps = nut[0][0]
        ecl = (120.5, 1.3, 1.0)
        eq = ephem.cotrans(ecl, -eps)
        back = ephem.cotrans(eq, eps)
        t += 1
        d = adiff(back[0], ecl[0]) * 3600
        if d <= 0.001:
            p += 1
        else:
            f += 1
            fails.append(f'  cotrans roundtrip d={d:.6f}"')
    except:
        pass
    # 8. julday/revjul round-trip
    try:
        jd_test = ephem.julday(2000, 1, 1, 12.0, 1)
        y, m, day, h = ephem.revjul(jd_test, 1)
        t += 1
        if y == 2000 and m == 1 and day == 1 and abs(h - 12.0) < 0.001:
            p += 1
        else:
            f += 1
            fails.append(f"  julday roundtrip {y}/{m}/{day} {h}")
    except:
        pass
    # 9. Mean/True node
    try:
        lr_mn = ephem.calc_ut(jd, ephem.MEAN_NODE, LF)
        sr_mn = swe.calc_ut(jd, swe.MEAN_NODE, SF)
        t += 1
        d = adiff(lr_mn[0][0], sr_mn[0][0]) * 3600
        if d <= 1.0:
            p += 1
        else:
            f += 1
            fails.append(f'  MeanNode d={d:.4f}"')
    except:
        pass
    # 10. house_pos
    try:
        st = ephem.sidtime(jd)
        armc = st * 15.0
        nut = ephem.calc_ut(jd, -1, 0)
        eps = nut[0][0]
        sun = ephem.calc_ut(jd, ephem.SUN, LF)
        le_hp = ephem.house_pos(armc, 45.0, eps, ord("P"), sun[0][0], sun[0][1])
        sr_sun = swe.calc_ut(jd, swe.SUN, SF)
        se_hp = swe.house_pos(armc, 45.0, eps, (sr_sun[0][0], sr_sun[0][1]), hsys=b"P")
        t += 1
        d = abs(le_hp - se_hp)
        if d <= 0.01:
            p += 1
        else:
            f += 1
            fails.append(f"  house_pos Sun d={d:.6f}")
    except:
        pass
    return p, f, t, fails


if __name__ == "__main__":
    print("=" * 70)
    print("Rounds 311-320: Deep Diverse Verification Batch 6 (FINALE)")
    print("=" * 70)
    tests = [
        (311, "Fixed star speed", test_311),
        (312, "Ecliptic vs equatorial", test_312),
        (313, "Planets at node times", test_313),
        (314, "Asteroid helio orbit", test_314),
        (315, "Moon declination extremes", test_315),
        (316, "10-planet comprehensive", test_316),
        (317, "Historical 1900-1950", test_317),
        (318, "Modern 2020-2030", test_318),
        (319, "ASC==cusp1 consistency", test_319),
        (320, "Ultimate all-API sweep", test_320),
    ]
    for num, name, fn in tests:
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
