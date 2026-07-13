# Mega Verification Plan — 21,666 Items (All Passing)

Last updated: 2026-03-27

> Clean-room note (2026-07): any `swe_ref` comparison named in this historical
> plan must be executed as an ephemeral public-API comparison. Numeric results
> must not be captured as fixtures, golden files, coefficients, tables, or
> generated artifacts. Current regression coverage should prefer independent
> JPL/ERFA/catalog truth, defining conditions, and metamorphic invariants.

## Status Summary

| Group | Script | Checks | Pass | Fail | Status |
|-------|--------|--------|------|------|--------|
| G01: Time functions | mega_g01_g02.py | 1,120 | 1,120 | 0 | DONE |
| G02: calc_ut positions | mega_g01_g02.py | 3,500 | 3,500 | 0 | DONE |
| G03: Flag combos exhaustive | mega_g03_g04.py | 5,480 | 5,480 | 0 | DONE |
| G04: House systems | mega_g03_g04.py | 3,600 | 3,600 | 0 | DONE |
| G05: Ayanamsha & sidereal | mega_g05_g06_g07.py | 1,042 | 1,042 | 0 | DONE |
| G06: Eclipses & occultations | mega_g05_g06_g07.py | 386 | 386 | 0 | DONE |
| G07: Rise/set/transit | mega_g05_g06_g07.py | 507 | 507 | 0 | DONE |
| G08: Fixed stars | mega_g08_to_g12.py | 400 | 400 | 0 | DONE |
| G09: Coord transforms | mega_g08_to_g12.py | 585 | 585 | 0 | DONE |
| G10: Phenomena & heliacal | mega_g08_to_g12.py | 630 | 630 | 0 | DONE |
| G11: Nodes, apsides, elements | mega_g08_to_g12.py | 600 | 600 | 0 | DONE |
| G12: Crossings & stations | mega_g08_to_g12.py | 300 | 300 | 0 | DONE |
| G13: Utility functions | mega_g13_to_g18.py | 1,084 | 1,084 | 0 | DONE |
| G14: State & context | mega_g13_to_g18.py | 239 | 239 | 0 | DONE |
| G15: Constants & aliases | mega_g13_to_g18.py | 1,006 | 1,006 | 0 | DONE |
| G16: Arabic parts | mega_g13_to_g18.py | 300 | 300 | 0 | DONE |
| G17: LEB backend | mega_g13_to_g18.py | 400 | 400 | 0 | DONE |
| G18: Edge cases & stress | mega_g13_to_g18.py | 487 | 487 | 0 | DONE |
| **TOTAL** | | **21,666** | **21,666** | **0** | **ALL PASS** |

Combined with Wave 1+2 scripts: **177,486 + 21,666 = 199,152 total checks**

## Scripts

| Script | Time | Groups |
|--------|------|--------|
| `mega_g01_g02.py` | 4.6s | G01 (Time), G02 (Positions) |
| `mega_g03_g04.py` | 3.9s | G03 (Flags), G04 (Houses) |
| `mega_g05_g06_g07.py` | 30s | G05 (Ayanamsha), G06 (Eclipses), G07 (Rise/Set) |
| `mega_g08_to_g12.py` | 2.7s | G08 (Stars), G09 (Coords), G10 (Pheno), G11 (Nodes), G12 (Crossings) |
| `mega_g13_to_g18.py` | 2.4s | G13 (Utils), G14 (State), G15 (Constants), G16 (Arabic), G17 (LEB), G18 (Edge) |
| **Total** | **~44s** | **All 18 groups** |

## Detailed Breakdown

### G01: Time Functions (1,120 checks)
- G01.01: julday — 20 dates × 5 checks (vs swe_ref, alias, round-trip, native type)
- G01.02: revjul — 20 JDs × 5 checks (year, month, day, hour vs swe_ref)
- G01.03: julday/revjul round-trip — 100 random dates
- G01.04: deltat — 50 JDs vs swe_ref, positive, < 0.01 day, smooth
- G01.05: deltat_ex — 50 JDs vs swe_ref, ephemeris flag
- G01.06: sidtime — 50 JDs vs swe_ref, range, sidereal gain
- G01.07: sidtime0 — 25 JDs vs swe_ref
- G01.08: utc_to_jd — 20 dates, jd_et + jd_ut match
- G01.09: jdet_to_utc / jdut1_to_utc — 20 JDs each, all components
- G01.10: time_equ — 25 JDs vs swe_ref
- G01.11: day_of_week — 50 dates known + vs swe_ref
- G01.12: TAI — UTC→TAI→UTC, TT→TAI→TT, leap seconds, JD conversions

### G02: calc_ut Positions (3,500 checks)
- G02.01: 10 bodies × 100 dates × 3 coords (lon < 3", lat < 1", dist < 5e-5 AU)
- G02.02: 5 node/apse bodies × 40 dates
- G02.03: 5 asteroids × 40 dates (< 6")
- G02.04: 10 bodies × 10 dates speed vs numerical derivative

### G03: Flag Combinations (5,480 checks)
- G03.01: 13 flags × 10 bodies × 10 dates (shape, finiteness, range)
- G03.02: 6 critical pairs × 5 bodies × 7 dates post-warmup (BUG-008)

### G04: House Systems (3,600 checks)
- G04.01: 24 systems × 5 locations × 10 dates (cusps + ASC + MC vs swe_ref)

### G05: Ayanamsha & Sidereal (1,042 checks)
- G05.01: 43 modes at J2000 + 2 other epochs (value, finite, range, name, ex variant, distinctness)
- G05.02: 10 modes × 10 bodies × 3 dates sidereal positions

### G06: Eclipses (386 checks)
- G06.01: Solar eclipse search 2000-2010 (time, type, geometry)
- G06.02: Lunar eclipse search 2000-2010 (time, type, magnitude)
- G06.03: 10 eclipse details (where, how, Besselian)

### G07: Rise/Set/Transit (507 checks)
- G07.01: Sun 5 locations × 20 dates
- G07.02: Moon 3 locations × 25 dates
- G07.03: 5 planets × 10 dates
- G07.04: Twilight (civil/nautical/astronomical)

### G08: Fixed Stars (400 checks)
- 20 stars × 5 dates × 3 coords vs swe_ref
- Magnitudes, TT variant, proper motion

### G09: Coordinate Transforms (585 checks)
- cotrans round-trip + vs swe_ref
- cotrans_sp with velocities
- azalt/azalt_rev round-trips
- refrac/refrac_extended round-trips

### G10: Phenomena (630 checks)
- pheno_ut 9 bodies × 10 dates (phase, elongation, magnitude)
- Elongation helpers (Mercury/Venus)
- vis_limit_mag 10 bodies × 5 conditions

### G11: Nodes/Apsides/Elements (600 checks)
- nod_aps_ut 5 bodies × 10 dates × 4 outputs
- Orbital elements 10 bodies × 10 dates
- orbit_max_min_true_distance
- Lunar node/apse properties

### G12: Crossings & Stations (300 checks)
- solcross_ut 4 targets × 25 years
- mooncross_ut 12 targets
- mooncross_node_ut
- find_station_ut (Mercury, Mars, Jupiter, Saturn)

### G13: Utility Functions (1,084 checks)
- degnorm 200 values, radnorm 100, difdeg* 200, split_deg 400
- get_planet_name, csnorm, d2l, deg_midp, rad_midp

### G14: State & Context (239 checks)
- set/get round-trips, close/reinit, EphemerisContext isolation, version

### G15: Constants & Aliases (1,006 checks)
- Body constants, flag constants, 100+ function alias pairs
- Calendar/eclipse constants, house system constants
- Return type consistency (native Python types)

### G16: Arabic Parts (300 checks)
- 10 dates × 5 locations × checks (Fortunae, Spiritus, range)

### G17: LEB Backend (400 checks)
- LEB vs Skyfield 14 bodies × 100 dates
- LEB flag fallback (TOPOCTR, XYZ, RADIANS, NONUT)
- LEB reader API

### G18: Edge Cases (487 checks)
- Date boundaries (-5000 to +3000)
- Polar latitudes (±60 to ±85)
- Invalid inputs (NaN, Inf, body 999)
- 360° wrap-around
- Special bodies (ECL_NUT, Earth geocentric, Sun heliocentric)
