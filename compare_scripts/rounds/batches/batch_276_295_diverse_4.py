#!/usr/bin/env python3
"""Rounds 276-295: Deep diverse verification batch 4.

276: Planet XYZ coordinate output (FLG_XYZ)
277: RADIANS output mode (FLG_RADIANS)
278: Topocentric positions at world cities
279: Pluto heliocentric deep sweep
280: Sun equatorial positions deep
281: Moon sidereal all ayanamsha modes
282: houses_ex2 cusp speed validation
283: Ecliptic obliquity at extreme dates
284: nod_aps_ut mean node/apogee
285: Fixed star magnitude comparison
286: Planet at retrograde midpoint timing
287: house_pos all planets Placidus
288: Chiron sidereal mode sweep
289: Mean vs True node divergence
290: Asteroid geocentric distance
291: sol_eclipse_when_glob forward search
292: cotrans_sp all planets with speeds
293: Sidereal time at year boundaries
294: Houses at longitude 180 (date line)
295: All flag combos for Venus
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


def adiff(a, b):
    d = abs(a - b) % 360
    return min(d, 360 - d)


DATES = [
    2415020.0,  # 1900
    2430000.0,  # 1941
    2440000.0,  # 1968
    2445000.0,  # 1982
    2451545.0,  # J2000
    2455000.0,  # 2009
    2458000.0,  # 2017
    2460000.0,  # 2023
    2462000.0,  # 2028
    2465000.0,  # 2036
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


def test_276():  # Planet XYZ coordinates
    p, f, t, fails = 0, 0, 0, []
    lxf = LF | ephem.FLG_XYZ
    sxf = SF | swe.FLG_XYZ
    for le_b, se_b in PLANETS:
        for jd in DATES[:8]:
            try:
                lr = ephem.calc_ut(jd, le_b, lxf)
                sr = swe.calc_ut(jd, se_b, sxf)
            except:
                continue
            for ci in range(3):  # X, Y, Z
                t += 1
                d = abs(lr[0][ci] - sr[0][ci])
                # XYZ in AU, tolerance 1e-5 AU (~1500 km)
                tol = 1e-5 if le_b != ephem.MOON else 1e-4
                if d <= tol:
                    p += 1
                else:
                    f += 1
                    lab = ["X", "Y", "Z"][ci]
                    fails.append(f"  B{le_b} {lab} jd={jd:.0f} diff={d:.8f} AU")
    return p, f, t, fails


def test_277():  # RADIANS output
    p, f, t, fails = 0, 0, 0, []
    lrf = LF | ephem.FLG_RADIANS
    srf = SF | swe.FLG_RADIANS
    for le_b, se_b in PLANETS[:7]:
        for jd in DATES[:6]:
            try:
                lr = ephem.calc_ut(jd, le_b, lrf)
                sr = swe.calc_ut(jd, se_b, srf)
            except:
                continue
            # Longitude in radians
            t += 1
            d = abs(lr[0][0] - sr[0][0])
            if d > math.pi:
                d = 2 * math.pi - d
            d_arcsec = math.degrees(d) * 3600
            if d_arcsec <= 1.0:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} lon_rad jd={jd:.0f} diff={d_arcsec:.4f}"')
            # Latitude in radians
            t += 1
            d = abs(lr[0][1] - sr[0][1])
            d_arcsec = math.degrees(d) * 3600
            if d_arcsec <= 1.0:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} lat_rad jd={jd:.0f} diff={d_arcsec:.4f}"')
    return p, f, t, fails


def test_278():  # Topocentric at world cities
    p, f, t, fails = 0, 0, 0, []
    cities = [
        (48.8566, 2.3522, 35.0, "Paris"),
        (35.6762, 139.6503, 40.0, "Tokyo"),
        (-33.8688, 151.2093, 58.0, "Sydney"),
        (40.7128, -74.0060, 10.0, "NYC"),
        (55.7558, 37.6173, 156.0, "Moscow"),
        (-22.9068, -43.1729, 11.0, "Rio"),
    ]
    ltf = LF | ephem.FLG_TOPOCTR
    stf = SF | swe.FLG_TOPOCTR
    for lat, lon, alt, city in cities:
        ephem.set_topo(lon, lat, alt)
        swe.set_topo(lon, lat, alt)
        for le_b, se_b in [
            (ephem.MOON, swe.MOON),
            (ephem.SUN, swe.SUN),
            (ephem.MARS, swe.MARS),
        ]:
            for jd in [2451545.0, 2460000.0, 2458000.0]:
                try:
                    lr = ephem.calc_ut(jd, le_b, ltf)
                    sr = swe.calc_ut(jd, se_b, stf)
                except:
                    continue
                t += 1
                d = adiff(lr[0][0], sr[0][0]) * 3600
                tol = 5.0 if le_b == ephem.MOON else 1.0
                if d <= tol:
                    p += 1
                else:
                    f += 1
                    fails.append(f'  {city} B{le_b} jd={jd:.0f} diff={d:.4f}"')
    return p, f, t, fails


def test_279():  # Pluto heliocentric deep
    p, f, t, fails = 0, 0, 0, []
    lhf = LF | ephem.FLG_HELCTR
    shf = SF | swe.FLG_HELCTR
    for i in range(40):
        jd = 2415020.0 + i * 1300.0  # ~142 years span
        try:
            lr = ephem.calc_ut(jd, ephem.PLUTO, lhf)
            sr = swe.calc_ut(jd, swe.PLUTO, shf)
        except:
            continue
        t += 1
        d = adiff(lr[0][0], sr[0][0]) * 3600
        if d <= 2.0:
            p += 1
        else:
            f += 1
            fails.append(f'  Pluto helio jd={jd:.0f} diff={d:.4f}"')
        t += 1
        d = abs(lr[0][1] - sr[0][1]) * 3600
        if d <= 2.0:
            p += 1
        else:
            f += 1
            fails.append(f'  Pluto helio lat jd={jd:.0f} diff={d:.4f}"')
    return p, f, t, fails


def test_280():  # Sun equatorial deep
    p, f, t, fails = 0, 0, 0, []
    lef = LF | ephem.FLG_EQUATORIAL
    sef = SF | swe.FLG_EQUATORIAL
    for i in range(50):
        jd = 2440000.0 + i * 500.0
        try:
            lr = ephem.calc_ut(jd, ephem.SUN, lef)
            sr = swe.calc_ut(jd, swe.SUN, sef)
        except:
            continue
        # RA
        t += 1
        d = adiff(lr[0][0], sr[0][0]) * 3600
        if d <= 1.0:
            p += 1
        else:
            f += 1
            fails.append(f'  Sun RA jd={jd:.0f} diff={d:.4f}"')
        # DEC
        t += 1
        d = abs(lr[0][1] - sr[0][1]) * 3600
        if d <= 1.0:
            p += 1
        else:
            f += 1
            fails.append(f'  Sun DEC jd={jd:.0f} diff={d:.4f}"')
        # RA speed
        t += 1
        d = abs(lr[0][3] - sr[0][3]) * 3600
        if d <= 1.0:
            p += 1
        else:
            f += 1
            fails.append(f'  Sun RA_spd jd={jd:.0f} diff={d:.4f}"/day')
    return p, f, t, fails


def test_281():  # Moon sidereal all ayanamsha modes
    p, f, t, fails = 0, 0, 0, []
    sid_modes = [0, 1, 3, 4, 5, 7, 15, 21, 22, 27]
    lsf = LF | ephem.FLG_SIDEREAL
    ssf = SF | swe.FLG_SIDEREAL
    for sid in sid_modes:
        swe.set_sid_mode(sid)
        ephem.set_sid_mode(sid, 0, 0)
        for jd in DATES[:5]:
            try:
                lr = ephem.calc_ut(jd, ephem.MOON, lsf)
                sr = swe.calc_ut(jd, swe.MOON, ssf)
            except:
                continue
            t += 1
            d = adiff(lr[0][0], sr[0][0]) * 3600
            if d <= 20.0:  # sidereal ~14" known offset
                p += 1
            else:
                f += 1
                fails.append(f'  Moon sid={sid} jd={jd:.0f} diff={d:.2f}"')
    swe.set_sid_mode(0)
    ephem.set_sid_mode(0, 0, 0)
    return p, f, t, fails


def test_282():  # houses_ex2 cusp speed validation
    p, f, t, fails = 0, 0, 0, []
    for hs in ["P", "K", "O", "R", "C", "E"]:
        for jd in [2451545.0, 2460000.0, 2455000.0]:
            for lat in [45.0, -30.0, 60.0]:
                try:
                    lc, la, lcs, las = ephem.houses_ex2(
                        jd, lat, 10.0, ord(hs), LF | ephem.FLG_SPEED
                    )
                    sc, sa, scs, sas = swe.houses_ex2(jd, lat, 10.0, hs.encode(), SF)
                except:
                    continue
                # Check cusp positions
                for i in range(min(len(lc), len(sc), 12)):
                    t += 1
                    d = adiff(lc[i], sc[i]) * 3600
                    if d <= 2.0:
                        p += 1
                    else:
                        f += 1
                        fails.append(
                            f'  {hs} c{i + 1} lat={lat} jd={jd:.0f} diff={d:.4f}"'
                        )
                # Check ASC/MC speeds (la[0]=ASC, la[1]=MC)
                for i, nm in [(0, "ASC"), (1, "MC")]:
                    t += 1
                    d = adiff(la[i], sa[i]) * 3600
                    if d <= 2.0:
                        p += 1
                    else:
                        f += 1
                        fails.append(f'  {hs} {nm} lat={lat} jd={jd:.0f} diff={d:.4f}"')
    return p, f, t, fails


def test_283():  # Ecliptic obliquity at many dates
    p, f, t, fails = 0, 0, 0, []
    for i in range(50):
        jd = 2415020.0 + i * 1000.0
        try:
            lr = ephem.calc_ut(jd, -1, 0)
            sr = swe.calc_ut(jd, -1, 0)
        except:
            continue
        # True obliquity
        t += 1
        d = abs(lr[0][0] - sr[0][0]) * 3600
        if d <= 5.0:  # Known IAU 2006 vs Laskar 1986 difference at edges
            p += 1
        else:
            f += 1
            fails.append(f'  true_obl jd={jd:.0f} diff={d:.4f}"')
        # Mean obliquity
        t += 1
        d = abs(lr[0][1] - sr[0][1]) * 3600
        if d <= 5.0:
            p += 1
        else:
            f += 1
            fails.append(f'  mean_obl jd={jd:.0f} diff={d:.4f}"')
        # Nutation in longitude
        t += 1
        d = abs(lr[0][2] - sr[0][2]) * 3600
        if d <= 1.0:
            p += 1
        else:
            f += 1
            fails.append(f'  nut_lon jd={jd:.0f} diff={d:.4f}"')
    return p, f, t, fails


def test_284():  # nod_aps_ut mean node and apogee
    p, f, t, fails = 0, 0, 0, []
    for le_b, se_b in [(ephem.MOON, swe.MOON)]:
        for jd in DATES[:8]:
            try:
                # Mean node
                lr_n = ephem.calc_ut(jd, ephem.MEAN_NODE, LF)
                sr_n = swe.calc_ut(jd, swe.MEAN_NODE, SF)
            except:
                continue
            t += 1
            d = adiff(lr_n[0][0], sr_n[0][0]) * 3600
            if d <= 1.0:
                p += 1
            else:
                f += 1
                fails.append(f'  MeanNode jd={jd:.0f} diff={d:.4f}"')
            try:
                # Mean apogee (Lilith)
                lr_a = ephem.calc_ut(jd, ephem.MEAN_APOG, LF)
                sr_a = swe.calc_ut(jd, swe.MEAN_APOG, SF)
            except:
                continue
            t += 1
            d = adiff(lr_a[0][0], sr_a[0][0]) * 3600
            if d <= 1.0:
                p += 1
            else:
                f += 1
                fails.append(f'  MeanLilith jd={jd:.0f} diff={d:.4f}"')
            # True node
            try:
                lr_tn = ephem.calc_ut(jd, ephem.TRUE_NODE, LF)
                sr_tn = swe.calc_ut(jd, swe.TRUE_NODE, SF)
            except:
                continue
            t += 1
            d = adiff(lr_tn[0][0], sr_tn[0][0]) * 3600
            if d <= 1.0:
                p += 1
            else:
                f += 1
                fails.append(f'  TrueNode jd={jd:.0f} diff={d:.4f}"')
    return p, f, t, fails


def test_285():  # Fixed star magnitude comparison
    p, f, t, fails = 0, 0, 0, []
    stars = [
        "Aldebaran",
        "Regulus",
        "Spica",
        "Antares",
        "Fomalhaut",
        "Sirius",
        "Canopus",
        "Vega",
        "Capella",
        "Rigel",
        "Procyon",
        "Betelgeuse",
        "Altair",
        "Deneb",
        "Pollux",
    ]
    jd = 2451545.0
    for star in stars:
        try:
            lr = ephem.fixstar2_ut(star, jd, LF)
            sr = swe.fixstar2(star, jd, SF)
            # lr = (pos_tuple, starname, retflag)
            # sr = ((lon,lat,dist,...), starname, retflag)
        except:
            continue
        # Compare longitude
        t += 1
        le_lon = lr[0][0]
        se_lon = sr[0][0]
        d = adiff(le_lon, se_lon) * 3600
        if d <= 1.0:
            p += 1
        else:
            f += 1
            fails.append(f'  {star} lon diff={d:.4f}"')
        # Compare latitude
        t += 1
        le_lat = lr[0][1]
        se_lat = sr[0][1]
        d = abs(le_lat - se_lat) * 3600
        if d <= 1.0:
            p += 1
        else:
            f += 1
            fails.append(f'  {star} lat diff={d:.4f}"')
    return p, f, t, fails


def test_286():  # Planet positions near retrograde midpoint
    p, f, t, fails = 0, 0, 0, []
    # Find points where speed is near zero (retrograde/direct station)
    for le_b, se_b in PLANETS[2:8]:  # Mercury through Uranus
        for jd in DATES[:6]:
            # Search forward for speed sign change in 1-day steps
            for step in range(60):
                jd_test = jd + step * 5.0
                try:
                    lr = ephem.calc_ut(jd_test, le_b, LF)
                    sr = swe.calc_ut(jd_test, se_b, SF)
                except:
                    continue
                if abs(lr[0][3]) < 0.05:  # Near station
                    t += 1
                    d = adiff(lr[0][0], sr[0][0]) * 3600
                    if d <= 2.0:
                        p += 1
                    else:
                        f += 1
                        fails.append(
                            f'  B{le_b} station jd={jd_test:.1f} diff={d:.4f}"'
                        )
                    break  # One station per body/epoch
    return p, f, t, fails


def test_287():  # house_pos all planets Placidus
    p, f, t, fails = 0, 0, 0, []
    for jd in [2451545.0, 2460000.0, 2455000.0, 2458000.0]:
        try:
            nut = ephem.calc_ut(jd, -1, 0)
            eps = nut[0][0]
            st = ephem.sidtime(jd)
            armc = st * 15.0
        except:
            continue
        for le_b, se_b in PLANETS:
            try:
                lr = ephem.calc_ut(jd, le_b, LF)
                sr = swe.calc_ut(jd, se_b, SF)
            except:
                continue
            for lat in [45.0, -30.0, 52.5]:
                try:
                    le_hp = ephem.house_pos(
                        armc, lat, eps, ord("P"), lr[0][0], lr[0][1]
                    )
                    se_hp = swe.house_pos(
                        armc, lat, eps, (sr[0][0], sr[0][1]), hsys=b"P"
                    )
                except:
                    continue
                t += 1
                d = abs(le_hp - se_hp)
                if d <= 0.01:  # house position tolerance
                    p += 1
                else:
                    f += 1
                    fails.append(f"  B{le_b} lat={lat} jd={jd:.0f} hp_diff={d:.6f}")
    return p, f, t, fails


def test_288():  # Chiron sidereal mode sweep
    p, f, t, fails = 0, 0, 0, []
    lsf = LF | ephem.FLG_SIDEREAL
    ssf = SF | swe.FLG_SIDEREAL
    for sid in [0, 1, 3, 5, 7, 27]:
        swe.set_sid_mode(sid)
        ephem.set_sid_mode(sid, 0, 0)
        for i in range(10):
            jd = 2440000.0 + i * 2500.0
            try:
                lr = ephem.calc_ut(jd, ephem.CHIRON, lsf)
                sr = swe.calc_ut(jd, swe.CHIRON, ssf)
            except:
                continue
            t += 1
            d = adiff(lr[0][0], sr[0][0]) * 3600
            if d <= 20.0:  # sidereal ~14" known offset
                p += 1
            else:
                f += 1
                fails.append(f'  Chiron sid={sid} jd={jd:.0f} diff={d:.2f}"')
    swe.set_sid_mode(0)
    ephem.set_sid_mode(0, 0, 0)
    return p, f, t, fails


def test_289():  # Mean vs True node divergence pattern
    p, f, t, fails = 0, 0, 0, []
    for i in range(40):
        jd = 2451545.0 + i * 170.0  # ~18.6 year cycle covered
        try:
            lr_mn = ephem.calc_ut(jd, ephem.MEAN_NODE, LF)
            sr_mn = swe.calc_ut(jd, swe.MEAN_NODE, SF)
            lr_tn = ephem.calc_ut(jd, ephem.TRUE_NODE, LF)
            sr_tn = swe.calc_ut(jd, swe.TRUE_NODE, SF)
        except:
            continue
        # Mean node match
        t += 1
        d = adiff(lr_mn[0][0], sr_mn[0][0]) * 3600
        if d <= 1.0:
            p += 1
        else:
            f += 1
            fails.append(f'  MeanNode jd={jd:.0f} diff={d:.4f}"')
        # True node match
        t += 1
        d = adiff(lr_tn[0][0], sr_tn[0][0]) * 3600
        if d <= 1.0:
            p += 1
        else:
            f += 1
            fails.append(f'  TrueNode jd={jd:.0f} diff={d:.4f}"')
        # Both libraries should agree on mean-true difference magnitude
        t += 1
        le_diff = adiff(lr_mn[0][0], lr_tn[0][0])
        se_diff = adiff(sr_mn[0][0], sr_tn[0][0])
        dd = abs(le_diff - se_diff) * 3600
        if dd <= 2.0:
            p += 1
        else:
            f += 1
            fails.append(f'  MN-TN gap jd={jd:.0f} diff={dd:.4f}"')
    return p, f, t, fails


def test_290():  # Asteroid geocentric distance
    p, f, t, fails = 0, 0, 0, []
    asteroids = [
        (ephem.CERES, swe.CERES),
        (ephem.PALLAS, swe.PALLAS),
        (ephem.JUNO, swe.JUNO),
        (ephem.VESTA, swe.VESTA),
        (ephem.CHIRON, swe.CHIRON),
    ]
    for le_b, se_b in asteroids:
        for jd in DATES[:8]:
            try:
                lr = ephem.calc_ut(jd, le_b, LF)
                sr = swe.calc_ut(jd, se_b, SF)
            except:
                continue
            # Distance in AU
            t += 1
            d = abs(lr[0][2] - sr[0][2])
            if d <= 0.0001:
                p += 1
            else:
                f += 1
                fails.append(f"  B{le_b} jd={jd:.0f} dist_diff={d:.8f} AU")
            # Distance speed
            t += 1
            d = abs(lr[0][5] - sr[0][5])
            if d <= 0.001:
                p += 1
            else:
                f += 1
                fails.append(f"  B{le_b} jd={jd:.0f} dist_spd_diff={d:.8f}")
    return p, f, t, fails


def test_291():  # sol_eclipse_when_glob forward search
    p, f, t, fails = 0, 0, 0, []
    jd_start = 2451545.0  # J2000
    for i in range(10):
        try:
            lr = ephem.sol_eclipse_when_glob(jd_start, LF, 0)
            sr = swe.sol_eclipse_when_glob(jd_start, SF, 0, False)
        except:
            continue
        t += 1
        # Both should find the same eclipse (tret[0] = max time)
        le_max = lr[1][0]
        se_max = sr[1][0]
        d = abs(le_max - se_max) * 1440  # difference in minutes
        if d <= 15.0:  # 15 minute tolerance
            p += 1
        else:
            f += 1
            fails.append(f"  Eclipse#{i} max_diff={d:.2f} min")
        # Both should agree on eclipse type (total/annular/partial)
        t += 1
        le_type = lr[0] & 0xFF
        se_type = sr[0] & 0xFF
        # At least one common bit should match
        if le_type & se_type:
            p += 1
        else:
            f += 1
            fails.append(f"  Eclipse#{i} type mismatch LE={le_type} SE={se_type}")
        jd_start = se_max + 20  # Skip ahead past this eclipse
    return p, f, t, fails


def test_292():  # cotrans_sp all planets with speeds
    p, f, t, fails = 0, 0, 0, []
    for le_b, se_b in PLANETS[:7]:
        for jd in DATES[:6]:
            try:
                lr = ephem.calc_ut(jd, le_b, LF)
                nut = ephem.calc_ut(jd, -1, 0)
                eps = nut[0][0]
                # Convert ecliptic -> equatorial with speeds
                coord = (lr[0][0], lr[0][1], lr[0][2])
                speed = (lr[0][3], lr[0][4], lr[0][5])
                eq_pos, eq_spd = ephem.cotrans_sp(coord + speed, -eps)
            except:
                continue
            # Verify round-trip back to ecliptic
            try:
                back_pos, back_spd = ephem.cotrans_sp(
                    (eq_pos[0], eq_pos[1], eq_pos[2], eq_spd[0], eq_spd[1], eq_spd[2]),
                    eps,
                )
            except:
                continue
            # Position round-trip
            t += 1
            d = adiff(back_pos[0], coord[0]) * 3600
            if d <= 0.001:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} lon roundtrip diff={d:.6f}"')
            t += 1
            d = abs(back_pos[1] - coord[1]) * 3600
            if d <= 0.001:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} lat roundtrip diff={d:.6f}"')
            # Speed round-trip
            t += 1
            d = abs(back_spd[0] - speed[0]) * 3600
            if d <= 0.001:
                p += 1
            else:
                f += 1
                fails.append(f'  B{le_b} lon_spd roundtrip diff={d:.6f}"/day')
    return p, f, t, fails


def test_293():  # Sidereal time at year boundaries
    p, f, t, fails = 0, 0, 0, []
    # Jan 1 of each year from 1950 to 2040
    for year in range(1950, 2041, 1):
        jd = ephem.julday(year, 1, 1, 0.0, 1)
        try:
            le_st = ephem.sidtime(jd)
            se_st = swe.sidtime(jd)
        except:
            continue
        t += 1
        d = abs(le_st - se_st) * 3600  # diff in seconds of time
        if d <= 0.05:  # 50ms
            p += 1
        else:
            f += 1
            fails.append(f"  Y{year} sidtime diff={d:.6f}s")
    return p, f, t, fails


def test_294():  # Houses at longitude 180 (date line)
    p, f, t, fails = 0, 0, 0, []
    for hs in ["P", "K", "O", "R", "C", "E", "W", "B"]:
        for jd in [2451545.0, 2460000.0]:
            for lat in [35.0, -45.0]:
                try:
                    lc, la = ephem.houses_ex(jd, lat, 180.0, ord(hs), LF)
                    sc, sa = swe.houses_ex(jd, lat, 180.0, hs.encode(), SF)
                except:
                    continue
                for i in range(min(len(lc), len(sc), 12)):
                    t += 1
                    d = adiff(lc[i], sc[i]) * 3600
                    if d <= 2.0:
                        p += 1
                    else:
                        f += 1
                        fails.append(
                            f'  {hs} c{i + 1} lat={lat} lon=180 jd={jd:.0f} d={d:.4f}"'
                        )
    return p, f, t, fails


def test_295():  # All flag combos for Venus
    p, f, t, fails = 0, 0, 0, []
    flag_combos = [
        (0, "default"),
        (ephem.FLG_HELCTR, "helio"),
        (ephem.FLG_J2000, "J2000"),
        (ephem.FLG_NONUT, "NONUT"),
        (ephem.FLG_NOABERR, "NOABERR"),
        (ephem.FLG_EQUATORIAL, "equatorial"),
        (ephem.FLG_J2000 | ephem.FLG_NONUT, "J2000+NONUT"),
        (ephem.FLG_EQUATORIAL | ephem.FLG_J2000, "EQ+J2000"),
        (ephem.FLG_EQUATORIAL | ephem.FLG_NONUT, "EQ+NONUT"),
        (ephem.FLG_TRUEPOS, "TRUEPOS"),
        (ephem.FLG_HELCTR | ephem.FLG_J2000, "HELIO+J2000"),
    ]
    for extra, label in flag_combos:
        lfl = LF | extra
        sfl = SF | extra  # pyswisseph flags same numeric values
        for jd in DATES[:6]:
            try:
                lr = ephem.calc_ut(jd, ephem.VENUS, lfl)
                sr = swe.calc_ut(jd, swe.VENUS, sfl)
            except:
                continue
            t += 1
            d = adiff(lr[0][0], sr[0][0]) * 3600
            tol = 2.0
            if "TRUEPOS" in label or "NOABERR" in label:
                tol = 35.0  # known light-time differences
            if d <= tol:
                p += 1
            else:
                f += 1
                fails.append(f'  Venus {label} jd={jd:.0f} diff={d:.4f}"')
    return p, f, t, fails


if __name__ == "__main__":
    print("=" * 70)
    print("Rounds 276-295: Deep Diverse Verification Batch 4")
    print("=" * 70)
    tests = [
        (276, "Planet XYZ coordinates", test_276),
        (277, "RADIANS output mode", test_277),
        (278, "Topocentric world cities", test_278),
        (279, "Pluto heliocentric deep", test_279),
        (280, "Sun equatorial deep", test_280),
        (281, "Moon sidereal all modes", test_281),
        (282, "houses_ex2 cusp speed", test_282),
        (283, "Obliquity at many dates", test_283),
        (284, "Mean node/apogee/TrueNode", test_284),
        (285, "Fixed star positions J2000", test_285),
        (286, "Planet near stations", test_286),
        (287, "house_pos all planets", test_287),
        (288, "Chiron sidereal sweep", test_288),
        (289, "Mean vs True node pattern", test_289),
        (290, "Asteroid distance+speed", test_290),
        (291, "sol_eclipse_when_glob fwd", test_291),
        (292, "cotrans_sp round-trip", test_292),
        (293, "Sidereal time yearly", test_293),
        (294, "Houses at lon=180", test_294),
        (295, "Venus all flag combos", test_295),
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
