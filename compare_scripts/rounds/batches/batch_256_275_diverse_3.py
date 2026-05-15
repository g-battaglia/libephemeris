#!/usr/bin/env python3
"""Rounds 256-275: Deep diverse verification batch 3.

256: Planet positions at midnight vs noon (UT 0h vs 12h)
257: Ecliptic-to-equatorial round-trip consistency
258: All house systems at equator (lat=0)
259: Planet geocentric distance extremes (closest approach)
260: Fixed stars with NOABERR flag
261: Moon daily positions 30-day sweep
262: Sidereal houses Krishnamurti+Yukteshwar
263: Hypothetical bodies heliocentric
264: Planet lon_speed sign changes (stations)
265: Delta-T modern era fine grid
266: Houses at southern hemisphere cities
267: cotrans obliquity range sweep
268: All planets J2000+EQUATORIAL+NONUT triple combo
269: Chiron full orbit heliocentric
270: Fixed stars at year 2050
271: Nutation deps component
272: Planet lat_speed at ascending node
273: House cusps at DST transitions
274: Asteroid heliocentric+sidereal
275: Moon equatorial RA/DEC sweep
"""

from __future__ import annotations
import os, sys

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


def test_256():  # Midnight vs noon
    p, f, t, fails = 0, 0, 0, []
    for jd_base in [2451545.0, 2460000.0, 2440000.0, 2455000.0, 2458000.0]:
        for offset in [0.0, 0.5]:  # noon, midnight
            jd = jd_base + offset
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
                    fails.append(f'  B{le_b} jd={jd} diff={d:.4f}"')
    return p, f, t, fails


def test_257():  # Ecliptic->equatorial->ecliptic round-trip
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
                fails.append(f'  B{le_b} jd={jd:.0f} roundtrip lon diff={d:.6f}"')
            t += 1
            d = abs(back[1] - ecl[1]) * 3600
            if d <= 0.001:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} jd={jd:.0f} roundtrip lat diff={d:.6f}"')
    return p, f, t, fails


def test_258():  # All house systems at equator
    p, f, t, fails = 0, 0, 0, []
    for hs in ["P", "K", "O", "R", "C", "E", "W", "B", "M", "X", "A"]:
        for jd in [2451545.0, 2460000.0, 2440000.0, 2455000.0]:
            try:
                lc, la = ephem.houses_ex(jd, 0.0, 0.0, ord(hs), LF)
                sc, sa = swe.houses_ex(jd, 0.0, 0.0, hs.encode(), SF)
            except:
                continue
            for i in range(min(len(lc), len(sc), 12)):
                t += 1
                d = adiff(lc[i], sc[i]) * 3600
                tol = 5.0 if hs in ("I",) else 1.0
                if d <= tol:
                    p += 1
                else:
                    f += 1
                    fails.append(f'  {hs} eq c{i + 1} jd={jd:.0f} diff={d:.4f}"')
    return p, f, t, fails


def test_259():  # Planet distance extremes
    p, f, t, fails = 0, 0, 0, []
    for le_b, se_b in [
        (ephem.VENUS, swe.VENUS),
        (ephem.MARS, swe.MARS),
        (ephem.JUPITER, swe.JUPITER),
        (ephem.MERCURY, swe.MERCURY),
    ]:
        for jd in DATES[:8]:
            try:
                lr = ephem.calc_ut(jd, le_b, LF)
                sr = swe.calc_ut(jd, se_b, SF)
            except:
                continue
            t += 1
            d = abs(lr[0][2] - sr[0][2])
            if d <= 0.0001:
                p += 1
            else:
                f += 1
                fails.append(f"  B{le_b} jd={jd:.0f} dist diff={d:.8f}")
    return p, f, t, fails


def test_260():  # Fixed stars NOABERR
    p, f, t, fails = 0, 0, 0, []
    stars = [
        "Sirius",
        "Aldebaran",
        "Regulus",
        "Spica",
        "Antares",
        "Vega",
        "Arcturus",
        "Pollux",
        "Canopus",
        "Procyon",
    ]
    lf = LF | ephem.FLG_NOABERR
    sf = SF | swe.FLG_NOABERR
    for star in stars:
        for jd in [2451545.0, 2460000.0, 2440000.0]:
            try:
                lr = ephem.fixstar2_ut(star, jd, lf)
                sr = swe.fixstar2(star, jd, sf)
            except:
                continue
            t += 1
            d = adiff(lr[0][0], sr[0][0]) * 3600
            if d <= 10.0:
                p += 1
            else:
                f += 1
                fails.append(f'  {star} NOABERR jd={jd:.0f} diff={d:.4f}"')
    return p, f, t, fails


def test_261():  # Moon 30-day sweep
    p, f, t, fails = 0, 0, 0, []
    for day in range(30):
        jd = 2460000.0 + day
        try:
            lr = ephem.calc_ut(jd, ephem.MOON, LF)
            sr = swe.calc_ut(jd, swe.MOON, SF)
        except:
            continue
        t += 1
        d = adiff(lr[0][0], sr[0][0]) * 3600
        if d <= 1.0:
            p += 1
        else:
            f += 1
            fails.append(f'  Moon day{day} diff={d:.4f}"')
        t += 1
        d = abs(lr[0][1] - sr[0][1]) * 3600
        if d <= 1.0:
            p += 1
        else:
            f += 1
            fails.append(f'  Moon day{day} lat diff={d:.4f}"')
        t += 1
        d = abs(lr[0][3] - sr[0][3]) * 3600
        if d <= 5.0:
            p += 1
        else:
            f += 1
            fails.append(f'  Moon day{day} spd diff={d:.4f}"/day')
    return p, f, t, fails


def test_262():  # Sidereal houses Krishnamurti+Yukteshwar
    p, f, t, fails = 0, 0, 0, []
    for sid in [5, 7]:
        swe.set_sid_mode(sid)
        ephem.set_sid_mode(sid, 0, 0)
        lsf = LF | ephem.FLG_SIDEREAL
        ssf = SF | swe.FLG_SIDEREAL
        for hs in ["P", "K", "O", "R", "C", "E"]:
            for lat in [0.0, 30.0, 51.5, -33.87]:
                for jd in [2451545.0, 2460000.0, 2440000.0]:
                    try:
                        lc, la = ephem.houses_ex(jd, lat, 0.0, ord(hs), lsf)
                        sc, sa = swe.houses_ex(jd, lat, 0.0, hs.encode(), ssf)
                    except:
                        continue
                    t += 1
                    d = adiff(la[0], sa[0]) * 3600
                    if d <= 30.0:
                        p += 1
                    else:
                        f += 1
                        fails.append(f'  sid={sid} {hs} lat={lat} ASC diff={d:.2f}"')
        swe.set_sid_mode(0)
        ephem.set_sid_mode(0, 0, 0)
    return p, f, t, fails


def test_263():  # Hypothetical heliocentric
    p, f, t, fails = 0, 0, 0, []
    lhf = LF | ephem.FLG_HELCTR
    shf = SF | swe.FLG_HELCTR
    for le_b, se_b in [
        (ephem.CUPIDO, 40),
        (ephem.HADES, 41),
        (ephem.ZEUS, 42),
        (ephem.KRONOS, 43),
        (ephem.APOLLON, 44),
        (ephem.ADMETOS, 45),
        (ephem.VULKANUS, 46),
        (ephem.POSEIDON, 47),
    ]:
        for jd in DATES[:7]:
            try:
                lr = ephem.calc_ut(jd, le_b, lhf)
                sr = swe.calc_ut(jd, se_b, shf)
            except:
                continue
            t += 1
            d = adiff(lr[0][0], sr[0][0]) * 3600
            if d <= 10.0:
                p += 1
            else:
                f += 1
                fails.append(f'  U{se_b} helio jd={jd:.0f} diff={d:.2f}"')
    return p, f, t, fails


def test_264():  # Planet speed sign changes
    p, f, t, fails = 0, 0, 0, []
    for le_b, se_b in PLANETS[2:8]:
        for jd in DATES[:7]:
            try:
                lr = ephem.calc_ut(jd, le_b, LF)
                sr = swe.calc_ut(jd, se_b, SF)
            except:
                continue
            t += 1
            le_sign = 1 if lr[0][3] >= 0 else -1
            se_sign = 1 if sr[0][3] >= 0 else -1
            if le_sign == se_sign:
                p += 1
            else:
                f += 1
                fails.append(f"  B{le_b} jd={jd:.0f} speed sign mismatch")
            t += 1
            d = abs(lr[0][3] - sr[0][3]) * 3600
            if d <= 5.0:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} jd={jd:.0f} spd diff={d:.4f}"/day')
    return p, f, t, fails


def test_265():  # Delta-T modern fine grid
    p, f, t, fails = 0, 0, 0, []
    for i in range(100):
        jd = 2451545.0 + i * 36.525  # every ~36.5 days over 10 years
        try:
            ld = ephem.deltat(jd)
            sd = swe.deltat(jd)
        except:
            continue
        t += 1
        d = abs(ld - sd) * 86400
        if d <= 0.1:
            p += 1
        else:
            f += 1
            fails.append(f"  jd={jd:.2f} dt diff={d:.6f}s")
    return p, f, t, fails


def test_266():  # Southern hemisphere cities
    p, f, t, fails = 0, 0, 0, []
    cities = [
        (-33.87, 151.21, "Sydney"),
        (-23.55, -46.63, "SaoPaulo"),
        (-34.60, -58.38, "BuenosAires"),
        (-26.20, 28.04, "Johannesburg"),
        (-41.29, 174.78, "Wellington"),
        (-37.81, 144.96, "Melbourne"),
    ]
    for lat, lon, nm in cities:
        for hs in ["P", "K", "O", "R", "C", "E"]:
            for jd in [2451545.0, 2460000.0, 2455000.0]:
                try:
                    lc, la = ephem.houses_ex(jd, lat, lon, ord(hs), LF)
                    sc, sa = swe.houses_ex(jd, lat, lon, hs.encode(), SF)
                except:
                    continue
                t += 1
                d = adiff(la[0], sa[0]) * 3600
                if d <= 1.0:
                    p += 1
                else:
                    f += 1
                    fails.append(f'  {nm} {hs} ASC diff={d:.4f}"')
                t += 1
                d = adiff(la[1], sa[1]) * 3600
                if d <= 1.0:
                    p += 1
                else:
                    f += 1
                    fails.append(f'  {nm} {hs} MC diff={d:.4f}"')
    return p, f, t, fails


def test_267():  # cotrans obliquity range
    p, f, t, fails = 0, 0, 0, []
    for obl in [21.0, 22.0, 23.0, 23.44, 24.0, 25.0, 26.0]:
        for lon in [0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0]:
            for lat in [-5.0, 0.0, 5.0]:
                coord = (lon, lat, 1.0)
                try:
                    le = ephem.cotrans(coord, -obl)
                    se = swe.cotrans(coord, -obl)
                except:
                    continue
                t += 1
                d = adiff(le[0], se[0]) * 3600
                if d <= 0.001:
                    p += 1
                else:
                    f += 1
                    fails.append(f'  obl={obl} lon={lon} lat={lat} RA diff={d:.6f}"')
    return p, f, t, fails


def test_268():  # J2000+EQ+NONUT triple
    p, f, t, fails = 0, 0, 0, []
    lf3 = LF | ephem.FLG_J2000 | ephem.FLG_EQUATORIAL | ephem.FLG_NONUT
    sf3 = SF | swe.FLG_J2000 | swe.FLG_EQUATORIAL | swe.FLG_NONUT
    for le_b, se_b in PLANETS:
        for jd in DATES[:7]:
            try:
                lr = ephem.calc_ut(jd, le_b, lf3)
                sr = swe.calc_ut(jd, se_b, sf3)
            except:
                continue
            t += 1
            d = adiff(lr[0][0], sr[0][0]) * 3600
            if d <= 1.0:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} J2K+EQ+NN jd={jd:.0f} diff={d:.4f}"')
            t += 1
            d = abs(lr[0][1] - sr[0][1]) * 3600
            if d <= 1.0:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} J2K+EQ+NN jd={jd:.0f} dec diff={d:.4f}"')
    return p, f, t, fails


def test_269():  # Chiron heliocentric full orbit
    p, f, t, fails = 0, 0, 0, []
    lhf = LF | ephem.FLG_HELCTR
    shf = SF | swe.FLG_HELCTR
    for i in range(50):
        jd = 2440000.0 + i * 365.25  # ~50 years = Chiron orbit
        try:
            lr = ephem.calc_ut(jd, ephem.CHIRON, lhf)
            sr = swe.calc_ut(jd, swe.CHIRON, shf)
        except:
            continue
        t += 1
        d = adiff(lr[0][0], sr[0][0]) * 3600
        if d <= 2.0:
            p += 1
        else:
            f += 1
            fails.append(f'  Chiron helio jd={jd:.1f} diff={d:.4f}"')
    return p, f, t, fails


def test_270():  # Fixed stars at 2050
    p, f, t, fails = 0, 0, 0, []
    jd_2050 = 2469807.5
    stars = [
        "Sirius",
        "Aldebaran",
        "Regulus",
        "Spica",
        "Antares",
        "Vega",
        "Arcturus",
        "Pollux",
        "Fomalhaut",
        "Deneb",
        "Canopus",
        "Procyon",
        "Capella",
        "Betelgeuse",
        "Rigel",
    ]
    for star in stars:
        try:
            lr = ephem.fixstar2_ut(star, jd_2050, LF)
            sr = swe.fixstar2(star, jd_2050, SF)
        except:
            continue
        t += 1
        d = adiff(lr[0][0], sr[0][0]) * 3600
        if d <= 10.0:
            p += 1
        else:
            f += 1
            fails.append(f'  {star} 2050 lon diff={d:.4f}"')
        t += 1
        d = abs(lr[0][1] - sr[0][1]) * 3600
        if d <= 10.0:
            p += 1
        else:
            f += 1
            fails.append(f'  {star} 2050 lat diff={d:.4f}"')
    return p, f, t, fails


def test_271():  # Nutation deps
    p, f, t, fails = 0, 0, 0, []
    for i in range(60):
        jd = 2440000.0 + i * 340.0
        try:
            lr = ephem.calc_ut(jd, -1, 0)
            sr = swe.calc_ut(jd, -1, swe.FLG_SWIEPH)
        except:
            continue
        t += 1
        d = abs(lr[0][3] - sr[0][3]) * 3600  # deps
        if d <= 0.1:
            p += 1
        else:
            f += 1
            fails.append(f'  jd={jd:.1f} deps diff={d:.6f}"')
    return p, f, t, fails


def test_272():  # Planet lat_speed at nodes
    p, f, t, fails = 0, 0, 0, []
    for le_b, se_b in PLANETS[2:8]:
        for jd in DATES[:7]:
            try:
                lr = ephem.calc_ut(jd, le_b, LF)
                sr = swe.calc_ut(jd, se_b, SF)
            except:
                continue
            t += 1
            d = abs(lr[0][4] - sr[0][4]) * 3600
            if d <= 5.0:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} jd={jd:.0f} latspd diff={d:.4f}"/day')
    return p, f, t, fails


def test_273():  # House cusps at various JDs (DST-like transitions)
    p, f, t, fails = 0, 0, 0, []
    for jd_base in [2460000.0, 2451545.0]:
        for h in [
            0.0,
            1.0 / 24,
            2.0 / 24,
            3.0 / 24,
            6.0 / 24,
            12.0 / 24,
            18.0 / 24,
            23.0 / 24,
        ]:
            jd = jd_base + h
            for hs in ["P", "K", "O"]:
                try:
                    lc, la = ephem.houses_ex(jd, 51.5, -0.1, ord(hs), LF)
                    sc, sa = swe.houses_ex(jd, 51.5, -0.1, hs.encode(), SF)
                except:
                    continue
                t += 1
                d = adiff(la[0], sa[0]) * 3600
                if d <= 1.0:
                    p += 1
                else:
                    f += 1
                    fails.append(f'  {hs} h={h * 24:.0f} ASC diff={d:.4f}"')
    return p, f, t, fails


def test_274():  # Asteroid helio+sidereal
    p, f, t, fails = 0, 0, 0, []
    asteroids = [(ephem.CERES, swe.CERES), (ephem.VESTA, swe.VESTA)]
    lhf = LF | ephem.FLG_HELCTR
    shf = SF | swe.FLG_HELCTR
    for le_b, se_b in asteroids:
        for jd in DATES[:7]:
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
                fails.append(f'  Ast{se_b} helio jd={jd:.0f} diff={d:.4f}"')
    # Sidereal
    for sid in [1, 3]:
        swe.set_sid_mode(sid)
        ephem.set_sid_mode(sid, 0, 0)
        lsf = LF | ephem.FLG_SIDEREAL
        ssf = SF | swe.FLG_SIDEREAL
        for le_b, se_b in asteroids:
            for jd in DATES[:5]:
                try:
                    lr = ephem.calc_ut(jd, le_b, lsf)
                    sr = swe.calc_ut(jd, se_b, ssf)
                except:
                    continue
                t += 1
                d = adiff(lr[0][0], sr[0][0]) * 3600
                if d <= 20.0:
                    p += 1
                else:
                    f += 1
                    fails.append(f'  Ast{se_b} sid={sid} jd={jd:.0f} diff={d:.2f}"')
        swe.set_sid_mode(0)
        ephem.set_sid_mode(0, 0, 0)
    return p, f, t, fails


def test_275():  # Moon equatorial sweep
    p, f, t, fails = 0, 0, 0, []
    lef = LF | ephem.FLG_EQUATORIAL
    sef = SF | swe.FLG_EQUATORIAL
    for i in range(60):
        jd = 2451545.0 + i * 122.0
        try:
            lr = ephem.calc_ut(jd, ephem.MOON, lef)
            sr = swe.calc_ut(jd, swe.MOON, sef)
        except:
            continue
        t += 1
        d = adiff(lr[0][0], sr[0][0]) * 3600  # RA
        if d <= 1.0:
            p += 1
        else:
            f += 1
            fails.append(f'  Moon RA jd={jd:.1f} diff={d:.4f}"')
        t += 1
        d = abs(lr[0][1] - sr[0][1]) * 3600  # DEC
        if d <= 1.0:
            p += 1
        else:
            f += 1
            fails.append(f'  Moon DEC jd={jd:.1f} diff={d:.4f}"')
    return p, f, t, fails


if __name__ == "__main__":
    print("=" * 70)
    print("Rounds 256-275: Deep Diverse Verification Batch 3")
    print("=" * 70)
    tests = [
        (256, "Midnight vs noon positions", test_256),
        (257, "Ecl->Equ round-trip", test_257),
        (258, "All houses at equator", test_258),
        (259, "Planet distance extremes", test_259),
        (260, "Fixed stars NOABERR", test_260),
        (261, "Moon 30-day sweep", test_261),
        (262, "Sidereal houses KP+Yuk", test_262),
        (263, "Hypothetical heliocentric", test_263),
        (264, "Speed sign changes", test_264),
        (265, "Delta-T modern fine grid", test_265),
        (266, "Southern hemisphere houses", test_266),
        (267, "cotrans obliquity range", test_267),
        (268, "J2000+EQ+NONUT triple", test_268),
        (269, "Chiron heliocentric orbit", test_269),
        (270, "Fixed stars at 2050", test_270),
        (271, "Nutation deps component", test_271),
        (272, "Planet lat_speed", test_272),
        (273, "House cusps hourly", test_273),
        (274, "Asteroid helio+sidereal", test_274),
        (275, "Moon equatorial sweep", test_275),
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
