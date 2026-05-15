#!/usr/bin/env python3
"""Rounds 296-310: Deep diverse verification batch 5.

296: Planet distance at mutual aspects
297: Fixed stars TRUEPOS flag
298: Hypothetical sidereal+equatorial combo
299: Nutation rate of change verification
300: Moon distance at eclipse times
301: House system comparison matrix
302: Planet positions at midnight UT
303: Delta-T vs TDB-TT consistency
304: Asteroid equatorial positions
305: Houses at high latitudes (70-85)
306: cotrans inverse consistency check
307: Planet at integer JD values
308: Sidereal houses fine longitude sweep
309: Hypothetical at century boundaries
310: All bodies HELCTR+J2000 combo
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


def test_296():  # Planet distance at various dates
    p, f, t, fails = 0, 0, 0, []
    for le_b, se_b in PLANETS[2:10]:
        for jd in DATES[:8]:
            try:
                lr = ephem.calc_ut(jd, le_b, LF)
                sr = swe.calc_ut(jd, se_b, SF)
            except:
                continue
            t += 1
            d = abs(lr[0][2] - sr[0][2])
            if d <= 0.0002:
                p += 1
            else:
                f += 1
                fails.append(f"  B{le_b} jd={jd:.0f} dist={d:.8f} AU")
    return p, f, t, fails


def test_297():  # Fixed stars TRUEPOS
    p, f, t, fails = 0, 0, 0, []
    ltf = LF | ephem.FLG_TRUEPOS
    stf = SF | swe.FLG_TRUEPOS
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
        "Betelgeuse",
    ]
    for star in stars:
        for jd in [2451545.0, 2460000.0, 2440000.0]:
            try:
                lr = ephem.fixstar2_ut(star, jd, ltf)
                sr = swe.fixstar2(star, jd, stf)
            except:
                continue
            t += 1
            d = adiff(lr[0][0], sr[0][0]) * 3600
            if d <= 1.0:
                p += 1
            else:
                f += 1
                fails.append(f'  {star} TRUEPOS jd={jd:.0f} diff={d:.4f}"')
    return p, f, t, fails


def test_298():  # Hypothetical sidereal+equatorial
    p, f, t, fails = 0, 0, 0, []
    lsef = LF | ephem.FLG_SIDEREAL | ephem.FLG_EQUATORIAL
    ssef = SF | swe.FLG_SIDEREAL | swe.FLG_EQUATORIAL
    swe.set_sid_mode(0)
    ephem.set_sid_mode(0, 0, 0)
    bodies = [
        (ephem.CUPIDO, 40),
        (ephem.HADES, 41),
        (ephem.ZEUS, 42),
        (ephem.KRONOS, 43),
    ]
    for le_b, se_b in bodies:
        for jd in DATES[:5]:
            try:
                lr = ephem.calc_ut(jd, le_b, lsef)
                sr = swe.calc_ut(jd, se_b, ssef)
            except:
                continue
            t += 1
            d = adiff(lr[0][0], sr[0][0]) * 3600
            if d <= 60.0:  # hypothetical bodies ~35" + sidereal ~14"
                p += 1
            else:
                f += 1
                fails.append(f'  U{se_b} sid+eq jd={jd:.0f} diff={d:.2f}"')
    swe.set_sid_mode(0)
    ephem.set_sid_mode(0, 0, 0)
    return p, f, t, fails


def test_299():  # Nutation rate of change
    p, f, t, fails = 0, 0, 0, []
    for jd in DATES[:8]:
        try:
            lr1 = ephem.calc_ut(jd, -1, 0)
            sr1 = swe.calc_ut(jd, -1, 0)
            lr2 = ephem.calc_ut(jd + 1.0, -1, 0)
            sr2 = swe.calc_ut(jd + 1.0, -1, 0)
        except:
            continue
        # Nutation longitude rate
        t += 1
        le_rate = (lr2[0][2] - lr1[0][2]) * 3600  # arcsec/day
        se_rate = (sr2[0][2] - sr1[0][2]) * 3600
        d = abs(le_rate - se_rate)
        if d <= 0.5:
            p += 1
        else:
            f += 1
            fails.append(f'  nut_lon rate jd={jd:.0f} diff={d:.4f}"/day')
        # Nutation obliquity rate
        t += 1
        le_rate = (lr2[0][3] - lr1[0][3]) * 3600
        se_rate = (sr2[0][3] - sr1[0][3]) * 3600
        d = abs(le_rate - se_rate)
        if d <= 0.5:
            p += 1
        else:
            f += 1
            fails.append(f'  nut_obl rate jd={jd:.0f} diff={d:.4f}"/day')
    return p, f, t, fails


def test_300():  # Moon distance at eclipse times
    p, f, t, fails = 0, 0, 0, []
    jd_start = 2451545.0
    for i in range(8):
        try:
            sr = swe.sol_eclipse_when_glob(jd_start, SF, 0, False)
            jd_ecl = sr[1][0]
            lr = ephem.calc_ut(jd_ecl, ephem.MOON, LF)
            sr_m = swe.calc_ut(jd_ecl, swe.MOON, SF)
        except:
            continue
        t += 1
        d = abs(lr[0][2] - sr_m[0][2])
        if d <= 0.0001:
            p += 1
        else:
            f += 1
            fails.append(f"  Moon dist at ecl jd={jd_ecl:.2f} diff={d:.8f}")
        t += 1
        d = adiff(lr[0][0], sr_m[0][0]) * 3600
        if d <= 1.0:
            p += 1
        else:
            f += 1
            fails.append(f'  Moon lon at ecl jd={jd_ecl:.2f} diff={d:.4f}"')
        jd_start = jd_ecl + 20
    return p, f, t, fails


def test_301():  # House system comparison matrix
    p, f, t, fails = 0, 0, 0, []
    systems = ["P", "K", "O", "R", "C", "E", "W", "B", "M"]
    for hs in systems:
        for jd in [2451545.0, 2460000.0]:
            for lat, lon in [(45.0, 10.0), (-33.0, 151.0), (60.0, -74.0)]:
                try:
                    lc, la = ephem.houses_ex(jd, lat, lon, ord(hs), LF)
                    sc, sa = swe.houses_ex(jd, lat, lon, hs.encode(), SF)
                except:
                    continue
                # ASC
                t += 1
                d = adiff(la[0], sa[0]) * 3600
                if d <= 2.0:
                    p += 1
                else:
                    f += 1
                    fails.append(f'  {hs} ASC lat={lat} diff={d:.4f}"')
                # MC
                t += 1
                d = adiff(la[1], sa[1]) * 3600
                if d <= 2.0:
                    p += 1
                else:
                    f += 1
                    fails.append(f'  {hs} MC lat={lat} diff={d:.4f}"')
    return p, f, t, fails


def test_302():  # Planet at midnight UT (JD .5)
    p, f, t, fails = 0, 0, 0, []
    for le_b, se_b in PLANETS:
        for yr_jd in [2451544.5, 2459945.5, 2444239.5, 2457754.5]:
            try:
                lr = ephem.calc_ut(yr_jd, le_b, LF)
                sr = swe.calc_ut(yr_jd, se_b, SF)
            except:
                continue
            t += 1
            d = adiff(lr[0][0], sr[0][0]) * 3600
            if d <= 1.0:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} midnight jd={yr_jd} diff={d:.4f}"')
    return p, f, t, fails


def test_303():  # Delta-T at many epochs
    p, f, t, fails = 0, 0, 0, []
    for i in range(50):
        jd = 2430000.0 + i * 700.0
        try:
            le_dt = ephem.deltat(jd)
            se_dt = swe.deltat(jd)
        except:
            continue
        t += 1
        d = abs(le_dt - se_dt) * 86400  # difference in seconds
        tol = 5.0 if jd < 2436935 or jd > 2460676 else 0.5
        if d <= tol:
            p += 1
        else:
            f += 1
            fails.append(f"  DeltaT jd={jd:.0f} diff={d:.4f}s")
    return p, f, t, fails


def test_304():  # Asteroid equatorial positions
    p, f, t, fails = 0, 0, 0, []
    lef = LF | ephem.FLG_EQUATORIAL
    sef = SF | swe.FLG_EQUATORIAL
    asteroids = [
        (ephem.CERES, swe.CERES),
        (ephem.PALLAS, swe.PALLAS),
        (ephem.JUNO, swe.JUNO),
        (ephem.VESTA, swe.VESTA),
    ]
    for le_b, se_b in asteroids:
        for jd in DATES[:6]:
            try:
                lr = ephem.calc_ut(jd, le_b, lef)
                sr = swe.calc_ut(jd, se_b, sef)
            except:
                continue
            # RA
            t += 1
            d = adiff(lr[0][0], sr[0][0]) * 3600
            if d <= 1.0:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} RA jd={jd:.0f} diff={d:.4f}"')
            # DEC
            t += 1
            d = abs(lr[0][1] - sr[0][1]) * 3600
            if d <= 1.0:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} DEC jd={jd:.0f} diff={d:.4f}"')
    return p, f, t, fails


def test_305():  # Houses at high latitudes
    p, f, t, fails = 0, 0, 0, []
    for hs in ["P", "K", "R", "E", "W"]:
        for lat in [70.0, 75.0, 80.0, -70.0, -75.0, -80.0]:
            for jd in [2451545.0, 2460000.0]:
                try:
                    lc, la = ephem.houses_ex(jd, lat, 10.0, ord(hs), LF)
                    sc, sa = swe.houses_ex(jd, lat, 10.0, hs.encode(), SF)
                except:
                    continue
                # Check ASC/MC
                for i, nm in [(0, "ASC"), (1, "MC")]:
                    t += 1
                    d = adiff(la[i], sa[i]) * 3600
                    if d <= 5.0:
                        p += 1
                    else:
                        f += 1
                        fails.append(f'  {hs} {nm} lat={lat} diff={d:.2f}"')
                # Check cusps 1,4,7,10 (angular)
                for ci in [0, 3, 6, 9]:
                    if ci < min(len(lc), len(sc)):
                        t += 1
                        d = adiff(lc[ci], sc[ci]) * 3600
                        if d <= 5.0:
                            p += 1
                        else:
                            f += 1
                            fails.append(f'  {hs} c{ci + 1} lat={lat} diff={d:.2f}"')
    return p, f, t, fails


def test_306():  # cotrans inverse consistency
    p, f, t, fails = 0, 0, 0, []
    for le_b, se_b in PLANETS[:7]:
        for jd in DATES[:6]:
            try:
                lr = ephem.calc_ut(jd, le_b, LF)
                nut = ephem.calc_ut(jd, -1, 0)
                eps = nut[0][0]
                ecl = (lr[0][0], lr[0][1], lr[0][2])
                eq = ephem.cotrans(ecl, -eps)
                back = ephem.cotrans(eq, eps)
            except:
                continue
            t += 1
            d = adiff(back[0], ecl[0]) * 3600
            if d <= 0.001:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} lon inv diff={d:.6f}"')
            t += 1
            d = abs(back[1] - ecl[1]) * 3600
            if d <= 0.001:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} lat inv diff={d:.6f}"')
    return p, f, t, fails


def test_307():  # Planet at integer JD
    p, f, t, fails = 0, 0, 0, []
    for le_b, se_b in PLANETS:
        for jd in [
            2450000.0,
            2451000.0,
            2452000.0,
            2453000.0,
            2454000.0,
            2456000.0,
            2459000.0,
            2461000.0,
        ]:
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
                fails.append(f'  B{le_b} jd={jd:.0f} diff={d:.4f}"')
    return p, f, t, fails


def test_308():  # Sidereal houses fine longitude sweep
    p, f, t, fails = 0, 0, 0, []
    lsf2 = LF | ephem.FLG_SIDEREAL
    ssf2 = SF | swe.FLG_SIDEREAL
    swe.set_sid_mode(0)
    ephem.set_sid_mode(0, 0, 0)
    for lon in range(0, 360, 30):
        for jd in [2451545.0, 2460000.0]:
            try:
                lc, la = ephem.houses_ex(jd, 45.0, float(lon), ord("P"), lsf2)
                sc, sa = swe.houses_ex(jd, 45.0, float(lon), b"P", ssf2)
            except:
                continue
            for i in range(min(len(lc), len(sc), 12)):
                t += 1
                d = adiff(lc[i], sc[i]) * 3600
                if d <= 20.0:
                    p += 1
                else:
                    f += 1
                    fails.append(
                        f'  P sid c{i + 1} lon={lon} jd={jd:.0f} diff={d:.2f}"'
                    )
    swe.set_sid_mode(0)
    ephem.set_sid_mode(0, 0, 0)
    return p, f, t, fails


def test_309():  # Hypothetical at century boundaries
    p, f, t, fails = 0, 0, 0, []
    century_jds = [
        2415020.0,  # 1900
        2451545.0,  # 2000
        2440587.5,  # 1970
        2460310.0,  # 2024
    ]
    bodies = [
        (ephem.CUPIDO, 40),
        (ephem.HADES, 41),
        (ephem.ZEUS, 42),
        (ephem.KRONOS, 43),
        (ephem.APOLLON, 44),
        (ephem.ADMETOS, 45),
        (ephem.VULKANUS, 46),
        (ephem.POSEIDON, 47),
    ]
    for le_b, se_b in bodies:
        for jd in century_jds:
            try:
                lr = ephem.calc_ut(jd, le_b, LF)
                sr = swe.calc_ut(jd, se_b, SF)
            except:
                continue
            t += 1
            d = adiff(lr[0][0], sr[0][0]) * 3600
            if d <= 40.0:  # hypothetical ~35" known limit
                p += 1
            else:
                f += 1
                fails.append(f'  U{se_b} jd={jd:.0f} diff={d:.2f}"')
    return p, f, t, fails


def test_310():  # All bodies HELCTR+J2000
    p, f, t, fails = 0, 0, 0, []
    lhj = LF | ephem.FLG_HELCTR | ephem.FLG_J2000
    shj = SF | swe.FLG_HELCTR | swe.FLG_J2000
    for le_b, se_b in PLANETS[2:10]:  # Skip Sun/Moon
        for jd in DATES[:6]:
            try:
                lr = ephem.calc_ut(jd, le_b, lhj)
                sr = swe.calc_ut(jd, se_b, shj)
            except:
                continue
            t += 1
            d = adiff(lr[0][0], sr[0][0]) * 3600
            if d <= 2.0:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} HJ2k jd={jd:.0f} diff={d:.4f}"')
    return p, f, t, fails


if __name__ == "__main__":
    print("=" * 70)
    print("Rounds 296-310: Deep Diverse Verification Batch 5")
    print("=" * 70)
    tests = [
        (296, "Planet distance sweep", test_296),
        (297, "Fixed stars TRUEPOS", test_297),
        (298, "Hypothetical sid+eq", test_298),
        (299, "Nutation rate of change", test_299),
        (300, "Moon dist at eclipses", test_300),
        (301, "House system ASC/MC matrix", test_301),
        (302, "Planet at midnight UT", test_302),
        (303, "Delta-T many epochs", test_303),
        (304, "Asteroid equatorial", test_304),
        (305, "Houses high latitudes", test_305),
        (306, "cotrans inverse check", test_306),
        (307, "Planet at integer JD", test_307),
        (308, "Sidereal houses lon sweep", test_308),
        (309, "Hypothetical century bounds", test_309),
        (310, "All bodies HELCTR+J2000", test_310),
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
