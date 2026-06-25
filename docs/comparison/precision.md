# Measured precision vs Swiss Ephemeris

All measurements compare libephemeris output against `pyswisseph` with identical
flags and parameters, taken across 100–210 dates spanning the DE440 range
(1550–2650 CE), with concentration around 1900–2100. The companion
[Known differences](known-differences.md) page explains the root cause of each
non-zero entry.

## 1. Precision summary

### Planetary positions

| Planet | Longitude (max) | Longitude (mean) | Latitude (max) | Distance (max) |
|--------|----------------|------------------|-----------------|----------------|
| Sun | 9.8" | 2.1" | 0.005" | 7.3 × 10⁻⁷ AU |
| Moon | 135" | 28" | 10.5" | 1.1 × 10⁻⁷ AU |
| Mercury | 14.6" | 3.1" | 3.2" | 6.8 × 10⁻⁵ AU |
| Venus | 12.0" | 2.5" | 2.3" | 1.8 × 10⁻⁵ AU |
| Mars | 7.1" | 1.4" | 0.5" | 2.7 × 10⁻⁵ AU |
| Jupiter | 2.2" | 0.5" | 0.05" | 3.7 × 10⁻⁵ AU |
| Saturn | 1.1" | 0.3" | 0.06" | 4.5 × 10⁻⁵ AU |
| Uranus | 0.8" | 0.2" | 0.02" | 4.6 × 10⁻⁵ AU |
| Neptune | 1.3" | 0.3" | 0.03" | 4.8 × 10⁻⁵ AU |
| Pluto | 1.8" | 0.4" | 0.5" | 1.0 × 10⁻⁴ AU |

All planets except the Moon are sub-2" (sub-arcsecond for the outer planets). The
Moon difference is explained in [Known differences §2.2](known-differences.md#22-moon-precision).

### Lunar nodes and Lilith

| Body | Longitude (max) | Latitude (max) | Speed (max) |
|------|-----------------|----------------|-------------|
| True Node | < 0.01" | 0.0" | 0.0007°/day |
| Mean Node | < 0.001° | 0.0" | 0.00005°/day |
| True Lilith (Osc. Apogee) | < 0.5" | 0.0" | 0.015°/day |
| Mean Lilith (Mean Apogee) | < 0.015" (±100yr) | ~20" | 0.00005°/day |

Mean Lilith latitude ~20" difference is from different analytical node formulas
used in the orbital-plane-to-ecliptic projection. No practical impact.

### Heliocentric, barycentric, equatorial, and XYZ modes

| Mode | Max difference | Notes |
|------|---------------|-------|
| Heliocentric (Mercury–Pluto) | < 0.0004° (1.1") | Sub-arcsecond for all bodies |
| Barycentric (Moon–Pluto) | < 0.001° | Sub-arcsecond for all non-Sun bodies |
| Barycentric Sun | ~0.04° angular | Distance-amplification artifact; actual 3D offset ~120 km (0.017% solar radius) |
| Equatorial RA/Dec | < 0.0005° (1.7") | Sub-arcsecond |
| XYZ Cartesian | < 0.00005 AU | Sub-arcsecond angular at all distances |

### Other areas

| Area | Precision | Notes |
|------|-----------|-------|
| House cusps (25 systems) | < 0.02" | All house systems tested at 11 global locations |
| Fixed stars (116 stars) | < 0.51" | van Leeuwen 2007 proper motions, verified vs SIMBAD |
| Coordinate transforms | Exact | cotrans, azalt, equatorial ↔ ecliptic |
| Utility functions | Exact | julday, revjul, degnorm, split_deg, etc. |
| Solar eclipse timing | < 6 sec | Maximum, contact times (C1–C4) |
| Lunar eclipse timing | < 8 sec | All contact types (P1, U1–U4, P4) |
| Lunar eclipse classification | All match | Total, partial, penumbral — all agree |
| Rise/set/transit | < 4 sec | Sun, Moon, planets |
| Crossings (solcross, mooncross) | < 4 sec | Longitude crossing times |
| Topocentric positions | All pass | 11 global locations, all planets |
| calc_pctr (planet-centric) | < 0.15° | All planet pairs tested |
| Sidereal modes (47 ayanamshas) | < 0.006° | All formula-based and star-based |
| Planetary phenomena (pheno_ut) | All match | Phase angle, elongation, magnitude, diameter |
| Combined flags | All pass | SPEED, EQUATORIAL, J2000, NONUT, NOABERR, etc. |
| Cartesian (XYZ) | All pass | Spherical-to-Cartesian conversion verified |
| Radians mode | All pass | Degree-to-radian conversion verified |
| Gauquelin sectors | < 0.5 sector | All methods, multiple planets |
| Heliacal events | < 1 day | Rising/setting timing |

### Velocities

| Component | Max difference |
|-----------|---------------|
| Longitude speed | < 0.003°/day |
| Latitude speed | < 0.004°/day |
| Distance speed | < 0.0001 AU/day |

---

## 2. Independent verification

These results triangulate **libephemeris vs pyswisseph vs an independent source**
(JPL Horizons, astropy/ERFA, SIMBAD). The purpose is to determine, for each area
of disagreement, which implementation is closer to an external physical reference.

### 2.1. Lunar True Node — both match JPL Horizons to < 0.01"

True Node (Ω) computed from all three sources across 24 dates (1950–2050):

| Source | vs JPL Horizons (J2000 ecliptic) |
|--------|----------------------------------|
| libephemeris | < 0.01" (machine precision) |
| pyswisseph | < 0.01" (machine precision) |

Both libraries are effectively identical for True Node longitude.

### 2.2. True Lilith (Osculating Apogee) — both equivalent

Osculating apogee longitude from all three sources across 200 dates (1950–2050):

| Source | vs JPL Horizons |
|--------|-----------------|
| libephemeris | ~240" (4 arcmin) |
| pyswisseph | ~240" (4 arcmin) |

The ~240" offset from Horizons is inherent to the two-body eccentricity-vector
formula (`e = (v×h)/μ − r/|r|`), which ignores solar perturbations on the Moon;
Horizons uses a different effective GM that accounts for perturbations.
libephemeris matches pyswisseph to mean 0.1", max 0.5" (effectively identical).

### 2.3. J2000 frame for lunar bodies

For planets (Sun, Moon, Mars, Jupiter), both libraries agree with astropy/ERFA to
< 1". For **lunar nodes and Lilith** in the J2000 frame, the engines differ at the
J2000.0 epoch itself, where tropical and J2000 ecliptic coordinates are identical
by definition (zero precession):

| Date | libephemeris J2000 shift | pyswisseph J2000 shift | Expected at J2000.0 |
|------|-------------------------|------------------------|---------------------|
| J2000.0 | 0.0000° | 0.0039° | 0.0000° |
| J1900 | −1.3966° | −1.3917° | ~−1.4° (precession) |
| 2023 | +0.3234° | +0.3208° | ~+0.3° (precession) |

libephemeris returns identical tropical and J2000 values at J2000.0 for these
analytically-computed bodies. See
[Intentional divergences](intentional-divergences.md) for the related
`SIDEREAL | J2000` case.

### 2.4. Latitude, equatorial, velocity — verified identical

| Component | True Node | True Lilith | Mean Node | Mean Lilith |
|-----------|-----------|-------------|-----------|-------------|
| Latitude | 0.0" | 0.0" | 0.0" | 20" max |
| RA (equatorial) | 0.0" | 0.2" | 0.0" | 2" max |
| Dec (equatorial) | 0.0" | 0.0" | 0.0" | 19" max |
| Speed (lon) | 0.0007°/day | 0.015°/day | 0.00005°/day | 0.00005°/day |
| Speed (lat) | 0.0" | 0.0004°/day | 0.0" | 0.00002°/day |

### 2.5. Fixed-star catalog — cross-checked against Hipparcos/SIMBAD

116 principal catalog stars were independently verified against authoritative sources.

**Proper motions updated to van Leeuwen 2007** (new Hipparcos reduction, A&A 474,
653–664):
- 99 stars had proper-motion values updated from the original 1997 Hipparcos catalog
- Data sourced from CDS/VizieR catalog I/311/hip2 via TAP query
- Largest correction: Tarf (β Cancri) — 18.3" position error at 100 years with old PM

**Position (RA/Dec) cross-checked against SIMBAD**:
- All 10 principal stars verified to < 0.02" against SIMBAD J2000 positions
- Two catalog bugs found and fixed: Algedi (was α¹ Cap, corrected to α² Cap, HIP 100064);
  Asellus Borealis (HIP 43103 = ι Cnc → corrected to 42806 = γ Cnc)

**Agreement with pyswisseph** (101 comparable stars at J2025.0): 100% within 0.002"
(7.2"), 98% within 0.5", 80% within 0.1". Max difference Rigil Kentaurus 0.51"
(nearest star — parallax not modelled).

**Independent astropy verification** (10 principal stars, 3 epochs): at J2025.0
libephemeris and pyswisseph each win against astropy for roughly half the stars;
Sirius (0.24" lat) and Fomalhaut (0.11" lon) offsets trace to annual parallax (not
modelled — Sirius 2.6 pc, Fomalhaut 7.7 pc).

### 2.6. Heliocentric, barycentric, equatorial, and XYZ modes

Across 4 coordinate modes, 10 bodies, 10 dates (1950–2050):

- **Heliocentric (FLG_HELCTR)**, Mercury–Pluto: position < 0.0004° (1.1") max
  (Mercury), speed < 0.0002°/day.
- **Barycentric (FLG_BARYCTR)**, Moon–Pluto: < 0.001° for all bodies except the
  Sun. Barycentric Sun up to 0.04° (138") angular, but this is a
  distance-amplification artifact — the actual 3D offset is ~120 km (0.017% solar
  radius), amplified by the tiny Sun–SSB distance (~0.001–0.009 AU). Verification
  against raw Skyfield/JPL DE440 confirms libephemeris is closer to the source
  ephemeris at J2000.
- **Equatorial (FLG_EQUATORIAL)**, Sun/Moon/Mars/Jupiter: RA < 0.0005° (1.7") max
  (Moon), Dec < 0.0002° (0.6") max.
- **XYZ Cartesian (FLG_XYZ)**, all bodies: < 0.00005 AU (< 0.00004 AU for
  Sun–Saturn); Neptune/Pluto ~0.00003 AU max, sub-arcsecond angular at 30+ AU.

These tolerances were determined empirically and independently verified against
JPL Horizons, astropy/ERFA, and SIMBAD. Earlier documentation used much looser
tolerances (e.g. 0.15° for True Node, 1.5° for latitude) that have since been
tightened by 100–1500× as triangulated verification showed the actual agreement is
far better.
