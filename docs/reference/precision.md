# LibEphemeris -- Scientific Precision

LibEphemeris is built on modern IAU standards and NASA JPL ephemerides. Every component -- from nutation to planet center corrections -- uses the most accurate models available in the astronomical literature. This document describes the scientific foundations, the specific models chosen, and measured precision for every calculation.

## Table of Contents

- [1. Ephemeris Foundation](#1-ephemeris-foundation)
- [2. Nutation Model](#2-nutation-model)
- [3. Precession Model](#3-precession-model)
- [4. Aberration and Light Deflection](#4-aberration-and-light-deflection)
- [5. Planet Center Corrections](#5-planet-center-corrections)
- [6. Velocity Computation](#6-velocity-computation)
- [7. True Ecliptic Coordinates](#7-true-ecliptic-coordinates)
- [8. House Calculations](#8-house-calculations)
- [9. Lunar Points](#9-lunar-points)
- [10. Delta T (TT − UT1)](#10-delta-t-tt--ut1)
- [11. Fixed Stars](#11-fixed-stars)
- [12. Ayanamsha (Sidereal Modes)](#12-ayanamsha-sidereal-modes)
- [13. Eclipses and Occultations](#13-eclipses-and-occultations)
- [14. Rise, Set, and Transit](#14-rise-set-and-transit)
- [15. Heliacal Events](#15-heliacal-events)
- [16. Planetary Positions -- Measured Precision](#16-planetary-positions----measured-precision)
- [17. Photometric Models (Phenomena)](#17-photometric-models-phenomena)
- [18. Minor Bodies](#18-minor-bodies)
- [19. Thread-safe Context API](#19-thread-safe-context-api)
- [20. Comprehensive Precision Audit](#20-comprehensive-precision-audit)
- [References](#references)

---

## 1. Ephemeris Foundation

### NASA JPL Development Ephemeris

LibEphemeris uses NASA JPL DE440 and DE441, the most recent planetary ephemerides produced by the Jet Propulsion Laboratory (Park et al. 2021). These are the same ephemerides used for Mars rover navigation and the James Webb Space Telescope.

| Property | DE440s | DE440 | DE441 |
|----------|--------|-------|-------|
| Date range | 1849--2150 CE | 1550--2650 CE | -13200--+17191 CE |
| File size | ~31 MB | ~119 MB | ~3.4 GB |
| Reference frame | ICRF 3.0 | ICRF 3.0 | ICRF 3.0 |
| Lunar model | Numerically integrated | Numerically integrated | Numerically integrated |
| Lunar Laser Ranging fit | ~1 milliarcsecond | ~1 milliarcsecond | ~1 milliarcsecond |
| Planetary radar fit | ~10 m (inner planets) | ~10 m (inner planets) | ~10 m (inner planets) |

DE440 and DE441 have identical precision -- DE441 is simply the extended-range version. DE440s is a reduced-size subset of DE440 with no loss of precision within its range.

### Precision tiers

LibEphemeris organizes these files into three precision tiers:

| Tier | File | Use Case |
|------|------|----------|
| `base` | de440s.bsp | Lightweight, modern-era usage |
| `medium` | de440.bsp | General purpose **(DEFAULT)** |
| `extended` | de441.bsp | Historical/far-future research |

Select by tier or by file directly:

```python
import libephemeris as eph

# By tier
eph.set_precision_tier("extended")   # uses de441.bsp

# Or by file
eph.set_ephemeris_file("de441.bsp")
```

Environment variables: `LIBEPHEMERIS_PRECISION=extended` or `LIBEPHEMERIS_EPHEMERIS=de441.bsp`.

Resolution priority (highest to lowest):
1. `LIBEPHEMERIS_EPHEMERIS` environment variable
2. `set_ephemeris_file()` / `set_jpl_file()` programmatic call
3. Precision tier (`LIBEPHEMERIS_PRECISION` env var or `set_precision_tier()`)
4. Default: `de440.bsp` (medium tier)

### FLG_MOSEPH handling

The `FLG_MOSEPH` flag — which in some engines selects a Moshier semi-analytical
fallback (VSOP87 for planets, ELP2000-82B for the Moon, with errors of ~1" for
inner planets and ~10+" for outer planets at historical dates) — is accepted for
API compatibility but ignored. Every libephemeris calculation always uses the full
JPL numerical integration.

---

## 2. Nutation Model

### IAU 2006/2000A (1365 terms)

LibEphemeris uses the **IAU 2006/2000A** nutation model via the IAU ERFA library (`erfa.nut06a()`). This is the highest-precision nutation model currently adopted by the International Astronomical Union.

| Property | Value |
|----------|-------|
| Lunisolar terms | 678 |
| Planetary terms | 687 |
| Total terms | **1365** |
| Precision | ~0.01--0.05 milliarcseconds |
| Standard | IERS Conventions 2010, ch. 5 |

The model computes nutation in longitude (Δψ) and nutation in obliquity (Δε) from the five Delaunay fundamental arguments plus nine planetary fundamental arguments, with corrections from the IAU 2006 precession-rate adjustments (Capitaine et al. 2005).

### Full-precision override

Skyfield's internal nutation (used before pyerfa integration) defaults to IAU
2000B, a truncated model with only **77 terms** and ~1 milliarcsecond precision.
LibEphemeris overrides this with the full 1365-term model via pyerfa for all planet
and house calculations.

### Impact on ecliptic longitude

Nutation in longitude (Δψ) directly shifts ecliptic longitudes. The maximum amplitude of Δψ is ~17.2 arcseconds (the 18.6-year lunar nodal cycle). Using a 77-term model vs. a 1365-term model introduces up to ~1 mas error in this correction, which over decades accumulates to measurable differences in planetary longitudes.

---

## 3. Precession Model

### Vondrák 2011 long-term precession

LibEphemeris uses the **Vondrák 2011** long-term precession model (Vondrák, Capitaine & Wallace 2011) via pyerfa's reference routines `erfa.ltp()` / `erfa.ltpb()` (see `libephemeris/precession_vondrak.py`). Unlike the IAU 2006 polynomial (Capitaine et al. 2003) — valid only a few centuries from J2000 and diverging at remote epochs (~36" for the Sun's longitude at year −3000) — the Vondrák model is fitted to a numerical integration and stays accurate over ±200,000 years. Near J2000 it agrees with IAU 2006 to better than 1 milliarcsecond, so modern results are unchanged.

| Property | Value |
|----------|-------|
| Model | Vondrák, Capitaine & Wallace (2011) long-term precession |
| Implementation | pyerfa `erfa.ltp()` / `erfa.ltpb()` |
| Valid range | ±200,000 years around J2000.0 |
| Agreement near J2000 | < 1 mas vs IAU 2006 (Capitaine et al. 2003) |
| J2000 obliquity | 84381.406 arcseconds (23°26'21.406") |
| Reference | Vondrák et al. (2011), A&A 534, A22 |

### Frame bias (ICRS to J2000)

The frame bias from the International Celestial Reference System (ICRS) to the mean equator and equinox is folded directly into the Vondrák precession matrix: `erfa.ltpb()` maps ICRS → mean equator/equinox of date **with** the ICRS frame bias included, while `erfa.ltp()` gives the same rotation without the bias (used on the ICRS-frame paths). The bias incorporates the frame bias angles (dα₀ = −14.6 mas, ξ₀ = −16.617 mas, η₀ = −6.819 mas) defined in the IERS Conventions.

---

## 4. Aberration and Light Deflection

### Annual aberration

For planet calculations, LibEphemeris uses Skyfield's full relativistic treatment via the `.apparent()` method, which implements IERS 2003 conventions for annual aberration (Bradley formula + relativistic corrections to order v²/c²).

The maximum annual aberration is ~20.5 arcseconds. Relativistic corrections to the classical Bradley formula amount to ~0.003 arcseconds.

### Gravitational light deflection

Skyfield's `.apparent()` method also applies gravitational light deflection by the Sun (maximum ~1.75 arcseconds at the solar limb, falling as 1/sin(θ)). This correction is significant for planets near solar conjunction.

### Diurnal aberration

LibEphemeris omits diurnal aberration (~0.3 arcseconds maximum), which is below the
threshold of astrological significance.

---

## 5. Planet Center Corrections

### The barycenter problem

JPL DE440/DE441 provide positions for outer planets as **system barycenters** -- the center of mass of the planet plus all its moons -- not the planet body center. The maximum angular offset between barycenter and body center:

| Planet | Max offset (geocentric) | Primary contributor |
|--------|------------------------|-------------------|
| Jupiter | ~0.6 arcsec | Ganymede + Callisto |
| Saturn | ~0.2 arcsec | Titan (96% of moon mass) |
| Uranus | ~0.01 arcsec | Titania + Oberon |
| Neptune | ~0.05 arcsec | Triton |
| Pluto | ~0.3 arcsec | Charon (binary system) |

For inner planets (Mercury, Venus, Earth, Mars), the barycenter effectively equals the body center because their moons (if any) have negligible mass relative to the planet.

### Three-tier correction strategy

LibEphemeris corrects barycenters to body centers **automatically** using a three-tier fallback:

#### Tier 1: SPK-based planet centers (<0.001 arcsec)

An auto-downloaded planet-centers kernel (~25 MB, `planet_centers_{tier}.bsp` per precision tier, legacy name `planet_centers.bsp`; not shipped in the wheel) contains precise center-of-body segments for NAIF IDs 599 (Jupiter), 699 (Saturn), 799 (Uranus), 899 (Neptune), and 999 (Pluto). These segments are extracted from JPL satellite ephemerides:

| NAIF ID | Source SPK | Coverage |
|---------|-----------|----------|
| 599 | jup204.bsp | ~1925--2025 |
| 699 | sat319.bsp | ~1950--2050 |
| 799 | ura083.bsp | ~1950--2050 |
| 899 | nep050.bsp | ~1989--2049 |
| 999 | plu017.bsp | ~1990--2050 |

The `_SpkCenterTarget` class in `planets.py` adds the SPK offset to the DE440 barycenter position at the retarded time (when light left the planet), ensuring the correction is applied correctly for light-time-corrected observations.

**If the date falls outside the SPK coverage**, Tier 1 silently falls back to Tier 2.

#### Tier 2: Analytical satellite theories (sub-arcsecond)

The `_CobCorrectedTarget` class computes the barycenter-to-center offset using analytical theories for major satellites. The correction formula:

```
offset = -Σ(m_i × r_i) / M_total
```

where m_i is the satellite mass and r_i is the satellite position relative to the planet center.

**Jupiter -- Lieske E5 theory (Meeus implementation)**

Theory for Io, Europa, Ganymede, and Callisto based on Lieske (1998), as presented in Meeus "Astronomical Algorithms" ch. 44.

| Property | Value |
|----------|-------|
| Reference epoch | JDE 2443000.5 (1976-08-10) |
| Terms per satellite | 21--50 (longitude), 7--11 (latitude), 7--16 (radius) |
| Coordinate frame | Jupiter equatorial → J2000 ecliptic (4 rotations) |
| Precession correction | From B1950.0: P = 1.3966626·T₀ + 0.0003088·T₀² |
| Positional precision | ~100--200 km (~0.05 arcsec at opposition) |

GM values from JUP365:

| Satellite | GM (km³/s²) | Mass ratio (m/M_Jupiter) |
|-----------|------------|--------------------------|
| Io | 5959.915 | 4.70 × 10⁻⁵ |
| Europa | 3202.712 | 2.53 × 10⁻⁵ |
| Ganymede | 9887.833 | 7.80 × 10⁻⁵ |
| Callisto | 7179.283 | 5.67 × 10⁻⁵ |
| **Total** | **26229.744** | **2.07 × 10⁻⁴** |

Jupiter GM = 126,686,531.9 km³/s² (JUP365).

**Saturn -- TASS 1.7 theory (Vienne & Duriez 1995)**

Full analytical theory for 8 Saturnian satellites, ported from the Stellarium implementation (Johannes Gajdosik, MIT license).

| Property | Value |
|----------|-------|
| Reference epoch | JD 2444240.0 (1980-01-01 TT) |
| Satellites | Mimas, Enceladus, Tethys, Dione, Rhea, Titan, Iapetus, Hyperion |
| Total terms | ~1500 (from `tass17_data.py`, 2964 lines of coefficients) |
| Orbital elements | 6 per satellite (n, λ, e·cos ω̃, e·sin ω̃, sin(i/2)·cos Ω, sin(i/2)·sin Ω) |
| Kepler solver | Newton's method, 15 iterations, convergence at 10⁻¹⁴ |
| Coordinate frame | TASS17 → VSOP87 (J2000 ecliptic) via fixed rotation matrix |
| Positional precision | ~50--100 km (~0.05 arcsec at opposition) |

Titan dominates the correction (GM = 8978.137 km³/s², 96% of total Saturnian moon mass). Saturn GM = 37,931,206.23 km³/s² (SAT441).

**Neptune -- Triton Keplerian model**

Simple Keplerian orbit with secular J2 nodal precession, based on NEP097 elements (Jacobson 2009).

| Property | Value |
|----------|-------|
| Semi-major axis | 354,759 km |
| Period | 5.877 days (retrograde) |
| Eccentricity | 0.000016 |
| Inclination | 156.865° (retrograde orbit) |
| Node precession | −360°/(688 × 365.25 days) (~688-year period) |
| Kepler solver | 1 iteration (sufficient for e ≈ 0) |
| Positional precision | ~20--50 km (~0.003 arcsec at opposition) |

Triton GM = 1428.495 km³/s², Neptune GM = 6,835,099.97 km³/s² (NEP097). Mass ratio = 2.09 × 10⁻⁴.

**Pluto -- Charon two-body Keplerian**

Two-body Keplerian solution from PLU060 (Brozovic & Jacobson 2024).

| Property | Value |
|----------|-------|
| Semi-major axis | 19,591 km |
| Period | 6.387 days |
| Eccentricity | 0.0002 (tidally locked) |
| Inclination | 96.145° |
| Node precession | 0 (negligible) |
| Kepler solver | 5 iterations |

Pluto-Charon is a **binary system**: the barycenter lies ~2130 km outside Pluto's center. Charon GM = 106.1 km³/s², Pluto GM = 869.3 km³/s² (PLU060). Mass ratio Charon/(Pluto+Charon) = **0.109**.

**Uranus**: Currently returns zero offset. The total mass ratio of Uranian moons to Uranus is ~1.0 × 10⁻⁴, corresponding to a maximum geocentric offset of ~0.01 arcsec. This is below the precision of most astrological applications.

#### Tier 3: Raw barycenters (last resort)

If both SPK and analytical methods fail, the raw barycenter position from DE440/DE441 is used. This path is not reached under normal operation.

### Velocity correction

The velocity component of the COB offset is computed via central difference numerical differentiation with a 1-second step:

```
v_offset = [offset(t + 0.5s) - offset(t - 0.5s)] / 1s
```

---

## 6. Velocity Computation

### Central difference method

LibEphemeris computes all velocities via numerical differentiation using the central difference formula:

```
dp/dt = [p(t + dt) - p(t - dt)] / (2 · dt)
```

| Context | Step size (dt) | Precision |
|---------|---------------|-----------|
| Planetary longitude/latitude | 1 second | ~0.0001°/day |
| COB offset velocity | 1 second | ~10⁻⁸ AU/day |
| Lunar node/Lilith | 0.5 days | ~0.001°/day |
| Fixed star apparent motion | 1 second | ~10⁻⁶°/day |

The central difference method is O(h²) accurate, meaning halving the step size reduces the error by a factor of 4.

### Accuracy

The central-difference method is accurate to <0.001°/day for all planets. For the
Moon, the numerical sampling introduces up to ~0.01°/day relative to an analytic
derivative.

---

## 7. True Ecliptic Coordinates

### From ICRS to ecliptic of date

The full pipeline from JPL ephemeris to ecliptic longitude:

1. **JPL state vector** (ICRS Cartesian, barycentric or geocentric)
2. **Light-time correction** (iterative, typically 3 Newton iterations)
3. **COB correction** (barycenter → planet center, Tier 1/2/3)
4. **Apparent position** via Skyfield: annual aberration + gravitational deflection
5. **Ecliptic of date**: rotate from ICRS to the true ecliptic of date by applying the Vondrák 2011 precession matrix combined with IAU 2006/2000A nutation (`precession_vondrak.vondrak_pn_matrix()`, built from `erfa.ltpb()`/`erfa.ltp()` and `erfa.nut06a()` via `erfa.numat()`)

For J2000 ecliptic (`FLG_J2000`), step 5 uses a fixed rotation by the J2000 obliquity (23°26'21.406") without precession or nutation.

### Mean obliquity of the ecliptic

The of-date mean obliquity uses the **Vondrák 2011** long-term obliquity series (evaluated in `libephemeris/sidereal_longterm.py`, exposed via `precession_vondrak.vondrak_mean_obliquity_rad()`) — the same long-term model as the precession above, valid over ±200,000 years. Near J2000 it matches the IAU 2006 polynomial (Capitaine et al. 2003) to sub-milliarcsecond precision:

```
ε₀ = 84381.406" − 46.836769"·T − 0.0001831"·T² + 0.00200340"·T³
     − 0.000000576"·T⁴ − 0.0000000434"·T⁵
```

where T is Julian centuries from J2000.0.

---

## 8. House Calculations

### Sidereal time

Houses require Greenwich Apparent Sidereal Time (GAST). LibEphemeris uses **Skyfield's `t.gast`**, which implements the IAU SOFA algorithm:

1. Earth Rotation Angle (ERA) from UT1 (IERS 2003 conventions)
2. Greenwich Mean Sidereal Time (GMST) via ERA + equation of origins
3. Equation of equinoxes: Δψ · cos(ε) (nutation in RA from IAU 2006/2000A)
4. GAST = GMST + equation of equinoxes

Precision: ~0.001 seconds of time (~0.015 arcseconds in hour angle).

### True obliquity

Houses use the true obliquity (mean obliquity + nutation in obliquity): the Vondrák 2011 of-date mean obliquity (`sidereal_longterm.mean_obliquity_deg()`) plus IAU 2006/2000A nutation (`erfa.nut06a()`). This is the same shared obliquity realization and nutation model used for planet calculations, keeping a chart's bodies and angles in one self-consistent frame.

### ARMC

```
ARMC = (GAST × 15°/h + geographic_longitude) mod 360°
```

### Supported systems (25)

All major house systems are implemented with standard spherical trigonometry:

`P` Placidus, `K` Koch, `O` Porphyry, `R` Regiomontanus, `C` Campanus, `E`/`A` Equal (from Asc), `D` Equal (from MC), `W` Whole Sign, `M` Morinus, `B` Alcabitius, `T` Polich-Page/Topocentric, `U` Krusinski, `G` Gauquelin (36-sector), `V` Vehlow, `X` Meridian, `H` Horizontal, `F` Carter, `S` Sripati, `L` Pullen SD, `Q` Pullen SR, `N` Natural Gradient, `Y` APC, `I` Sunshine (Treindl), `i` Sunshine (Makransky), `J` Savard-A.

### Convergence tolerance

Iterative systems (Placidus, Koch) use a convergence threshold of 10⁻⁷ degrees (~0.00036 arcseconds).

### Polar latitude behavior

Above the polar circle (~66.56° = 90° − obliquity), Placidus and Koch are geometrically undefined because some ecliptic points never rise or set. LibEphemeris raises `PolarCircleError` with the option to fall back to Porphyry via `houses_with_fallback()`.

### Measured precision

| Component | Max difference |
|-----------|---------------|
| Ascendant | <0.001° |
| MC (Midheaven) | <0.001° |
| ARMC | <0.001° |
| House cusps (typical) | 0.001°--0.01° |
| Vertex | <0.01° |

---

## 9. Lunar Points

### True Node (osculating lunar node)

The True Node is where the Moon's **instantaneous orbital plane** intersects the ecliptic. LibEphemeris uses a two-stage approach.

**Stage 1: Geometric osculating node from JPL DE440**

1. Get Moon's geocentric position **r** and velocity **v** from JPL DE440 (~1 milliarcsecond source precision)
2. Compute angular momentum vector: **h** = **r** × **v**
3. This vector is perpendicular to the instantaneous orbital plane
4. Find intersection with the ecliptic: Ω = atan2(h_x, −h_y)
5. Apply Vondrák 2011 precession (J2000 → date)
6. Apply IAU 2006/2000A nutation (1365 terms)

This computes **exactly** what the True Node is by definition: the intersection of the orbital plane with the ecliptic. Deriving the node directly from the state vectors avoids the approximation errors inherent in analytical series — the geometric construction is exact by definition.

**ELP2000-82B perturbation series (available but not applied)**

The codebase contains a comprehensive ELP2000-82B perturbation series (~170 terms, ~900 lines) organized into the following categories:

| Category | Terms | Dominant amplitude |
|----------|-------|-------------------|
| Solar (main) | 9 | −1.5233° · sin(2D) |
| Solar (2nd order) | 12 | 0.003°--0.04° |
| Solar (3rd order) | 10 | 0.001°--0.006° |
| Inclination (F-related) | 9 | 0.001°--0.01° |
| Venus perturbation | 9 | 0.001°--0.005° |
| Mars perturbation | 9 | 0.001°--0.004° |
| Jupiter perturbation | 7 | 0.001°--0.003° |
| Saturn perturbation | 7 | 0.001°--0.003° |
| Evection | 6 | up to 0.047° |
| Variation | 10 | up to 0.052° |
| Annual equation | 6 | up to 0.186° (E · sin M) |
| Parallactic inequality | 3 | up to 0.035° |
| Second-order coupling | 8 | 0.0001°--0.003° |
| Secular terms | 3+ | T-dependent corrections to T⁵ |

The Earth eccentricity factor E = 1.0 − 0.002516·T − 0.0000074·T² is applied to Sun-dependent terms per Meeus convention.

> **Note:** These perturbation corrections are **not currently applied** to the True Node calculation. Investigation revealed that the ELP2000-82B series was designed for the *mean* lunar node, not the geometric node derived from state vectors via `h = r × v`. Applying the series to the geometric node produced errors of tens of degrees, confirming the two approaches are incompatible. The geometric method from Stage 1 is used alone. The residual against an analytical mean-node series (~8.9 arcsec mean, ~0.14° max) reflects the fundamental methodological difference. See [Precision History](../development/precision-history.md) for the full investigation record.

### Mean Node

Meeus polynomial for the mean ascending node longitude, extended with T⁴ and T⁵ corrections from Chapront et al. (2002) and Simon et al. (1994).

### Mean Lilith (Black Moon)

Mean longitude of the lunar apogee, computed from Meeus polynomials (ch. 50). Includes latitude computation. Velocity via central difference with 0.5-day step.

### True Lilith (osculating apogee)

Computed from JPL DE440 Moon state vectors via orbital mechanics: the eccentricity vector **e** is derived from position and velocity, and the apogee direction points from the focus along **e**.

### Interpolated Apogee and Perigee

Uses the Moshier analytical method: ELP2000-82B perturbation series (~50 harmonic terms) fitted to DE404 over a 7000-year span (Moshier 1992). Removes the spurious ~30° oscillations present in the osculating apsides.

### Meeus polynomial validity ranges

| Range | |T| (centuries) | Accuracy |
|-------|---------------|----------|
| Optimal | <2 (~±200 years) | <0.001° |
| Valid | <10 (~±1000 years) | <0.01° |
| Maximum | <20 (~±2000 years) | Significant degradation |

Fundamental arguments include T⁵ corrections from Chapront et al. (2002).

### Measured precision

| Point | Max difference | Mean difference | Independent verification |
|-------|---------------|-----------------|------------------------|
| Mean Node | < 0.001° | < 0.001° | — |
| True Node | < 0.01" | < 0.01" | **Verified vs JPL Horizons** to < 0.01" (machine precision) across 24 dates (1950–2050). Both libraries match Horizons identically. |
| Mean Lilith | < 0.015" (lon, ±100yr) | < 0.01" | Latitude ~20" systematic difference from different analytical node formulas. No practical impact. |
| True Lilith | < 0.5" | ~0.1" | **Verified vs JPL Horizons**: both libraries show ~240" offset from Horizons (inherent two-body approximation). Libraries match each other to < 0.5". |
| Interpolated Apogee | ~1300 arcsec (0.36°) | ~350 arcsec (0.10°) | Genuine algorithm difference (JPL DE440 vs ELP2000). |
| Interpolated Perigee | ~9400 arcsec (2.6°) | ~1650 arcsec (0.46°) | Intentional — JPL DE440 physical passages vs truncated ELP2000. |

**J2000 frame (FLG_J2000) for lunar bodies:** at the J2000.0 epoch, tropical and
J2000 ecliptic coordinates are identical by definition (zero precession), and
LibEphemeris correctly returns zero shift for analytically-computed bodies, verified
independently against astropy/ERFA.

---

## 10. Delta T (TT − UT1)

### Model

LibEphemeris uses the **Stephenson, Morrison & Hohenkerk (2016)** Delta T model via Skyfield:

| Date range | Method |
|------------|--------|
| 720 BC -- ~2016 AD | Cubic spline interpolation from Table S15 |
| Outside spline | Parabolic: ΔT = −320 + 32.5 · u² (u = (year − 1825)/100) |
| 1973 -- present | IERS observed values (optional, via `set_iers_delta_t_enabled(True)`) |

### IERS observed Delta T

When enabled, LibEphemeris uses observed ΔT values from the International Earth Rotation and Reference Systems Service for recent dates (1973--present), achieving ~0.1 second precision.

```python
from libephemeris import set_iers_delta_t_enabled, download_delta_t_data
download_delta_t_data()
set_iers_delta_t_enabled(True)
```

### Typical Delta T values

| Year | ΔT (seconds) |
|------|-------------|
| −3000 | ~72,000 |
| 1000 | ~1,500 |
| 1800 | ~14 |
| 1900 | ~−3 |
| 2000 | ~64 |
| 2020 | ~69 |

---

## 11. Fixed Stars

### Star catalog

1,447 stars from the **Hipparcos catalog** (ESA 1997), with proper motions updated to the **van Leeuwen 2007 new Hipparcos reduction** (A&A 474, 653-664). The catalog covers all bright and astrologically significant stars: the 4 Royal Stars, Behenian stars, Pleiades cluster, Hyades, and full zodiacal constellation coverage.

| Property | Value |
|----------|-------|
| Epoch | J2000.0 (ICRS) |
| Source | Hipparcos catalog (HIP numbers) |
| Proper motions | van Leeuwen 2007 (I/311/hip2) for 99 stars; original Hipparcos 1997 for remainder |
| Count | 1,447 stars |
| Data per star | RA, Dec, PM_RA (with cos δ), PM_Dec, visual magnitude |

**Independent verification:** All stars cross-checked against SIMBAD J2000 positions. Principal stars verified to < 0.02" against SIMBAD. Two catalog bugs found and fixed during audit (Algedi wrong component, Asellus Borealis wrong HIP number).

### Proper motion

Rigorous space motion propagation from Hipparcos Vol. 1, Section 1.5.5, with second-order Taylor term for celestial sphere curvature:

```
correction = −0.5 · |V|² · P · t²
```

Precision: <0.01 arcsec over ±100 years; <1 arcsec over ±500 years.

### Position pipeline

1. Proper motion propagation (J2000 → date)
2. Create Skyfield `Star` object with propagated RA/Dec
3. Apparent position via `observer.at(t).observe(star).apparent()` (aberration + gravitational deflection)
4. Ecliptic coordinates via the Vondrák 2011 precession matrix + IAU 2006/2000A nutation (`precession_vondrak.vondrak_pn_matrix()` with `erfa.nut06a()`)

### Supported flags

The catalog uses the updated van Leeuwen 2007 proper motions. All meaningful FLG
flags are supported for fixed-star calculations: `FLG_SIDEREAL`, `FLG_J2000`,
`FLG_NONUT`, `FLG_XYZ`, `FLG_RADIANS`, `FLG_TRUEPOS`, `FLG_MOSEPH`, `FLG_SPEED3`,
`FLG_TOPOCTR`.

---

## 12. Ayanamsha (Sidereal Modes)

47 predefined ayanamsha systems are supported, covering the full set of standard sidereal modes.

### Formula-based ayanamshas

| Category | Examples | Max difference (measured) |
|----------|----------|---------------------------|
| Standard | Fagan-Bradley, Lahiri, Raman | <0.0002° |
| Epoch-based | J2000, J1900, B1950 | <0.0002° |
| Historical | Babylonian variants (Kugler 1-3, Huber) | <0.001° |

### Star-based ayanamshas

Star-based ayanamshas anchor the sidereal zodiac to a specific fixed star. The slightly higher differences reflect different proper motion and precession models:

| Ayanamsha | Max difference |
|-----------|---------------|
| True Citra (Spica, HIP 65474) | <0.006° |
| True Revati | <0.006° |
| True Pushya | <0.006° |
| True Mula | <0.006° |
| Galactic Center variants | <0.001° |

The Vondrák 2011 long-term precession model (used by LibEphemeris) is more accurate than the older Lieske (1977) model for computing the precession of the equinox, which directly affects ayanamsha values.

---

## 13. Eclipses and Occultations

### Solar eclipses

Calculated using Besselian elements with the full JPL DE ephemeris for Sun and Moon positions. Contact times (C1--C4), path width, central line coordinates, magnitude, and obscuration are computed.

| Property | Precision |
|----------|-----------|
| Eclipse maximum timing | <10 seconds (measured) |
| Contact times | <10 seconds |
| Eclipse magnitude | ~0.01 |
| Path width | ~1 km |

### Lunar eclipses

Full implementation including umbral and penumbral contact times (P1, U1, U2, U3, U4, P4), umbral/penumbral magnitude, gamma parameter, and duration.

### Saros and Inex series

Saros and Inex series numbers are computed from eclipse-to-eclipse relationships using the Saros period (6585.32 days) and Inex period (10571.95 days).

---

## 14. Rise, Set, and Transit

Computed using the Bennett (1982) atmospheric refraction formula with configurable atmospheric conditions (pressure, temperature).

| Event | Precision (measured) |
|-------|----------------------|
| Sunrise/sunset | <30 seconds |
| Moonrise/moonset | <30 seconds |
| Meridian transit | <30 seconds |

---

## 15. Heliacal Events

Schaefer (1990) atmospheric visibility model with:

1. **Atmospheric extinction**: Rayleigh scattering + aerosol + ozone + water vapor
2. **Twilight sky brightness**: gradient model with moonlight contribution
3. **Ptolemaic visibility thresholds**: arcus visionis values
4. **Contrast threshold**: Schaefer model with configurable observer skill levels

Timing precision: <1 day for heliacal rising/setting events.

---

## 16. Planetary Positions -- Measured Precision

Validated across 100+ random dates within the DE440 range (1550--2650 CE):

### Longitude

| Planet | Max difference | Mean difference |
|--------|---------------|-----------------|
| Sun | 0.20 arcsec | 0.04 arcsec |
| Moon | 3.32 arcsec | 0.70 arcsec |
| Mercury | 0.32 arcsec | 0.05 arcsec |
| Venus | 0.33 arcsec | 0.08 arcsec |
| Mars | 0.58 arcsec | 0.06 arcsec |
| Jupiter | 0.44 arcsec | 0.12 arcsec |
| Saturn | 0.51 arcsec | 0.13 arcsec |
| Uranus | 0.50 arcsec | 0.23 arcsec |
| Neptune | 1.17 arcsec | 0.24 arcsec |
| Pluto | 0.75 arcsec | 0.26 arcsec |

All differences are **sub-arcsecond** except the Moon (~3 arcseconds), which reflects the nutation-model and COB-correction pipeline choices. All values are well within the precision of the underlying DE ephemeris.

### Latitude

| Planet | Max difference |
|--------|---------------|
| Sun | <0.1 arcsec |
| Moon | <1.3 arcsec |
| Other planets | <0.6 arcsec |

### Velocity

| Component | Max difference |
|-----------|---------------|
| Angular velocity | <0.0004°/day |
| Radial velocity | <0.001 AU/day |

### Source of measured differences

These measured differences are **not errors**. They arise from intentional methodological choices:

| Source | Typical contribution |
|--------|---------------------|
| Nutation model (IAU 2006/2000A vs internal) | ~0.01--0.05 mas |
| COB correction (analytical vs none) | ~0.01--0.6 arcsec (outer planets) |
| DE440 vs DE431 ephemeris generation | ~0.001 arcsec |
| Velocity method (numerical vs analytical) | ~0.0001°/day |
| Light-time iteration tolerance | ~0.001 arcsec |

---

## 17. Photometric Models (Phenomena)

LibEphemeris implements `pheno_ut()` / `pheno()` for computing observable planetary phenomena: phase angle, phase (illuminated fraction), elongation, apparent diameter, and visual magnitude. The photometric models use peer-reviewed formulas from the astronomical literature, validated against astropy and the Astronomical Almanac.

### Phase Angle

The phase angle (Sun-Body-Earth angle) is computed using 3D vector geometry (dot product of geocentric position vectors), which is numerically stable for all configurations including the extremely elongated Sun-Moon-Earth triangle. This avoids the numerical instability of the law-of-cosines approach for bodies at very different distances.

| Body | Max diff (measured) | Notes |
|------|---------------------|-------|
| Sun | 0" | Exact (always 0) |
| Moon | < 1" | Irreducible ephemeris difference |
| Inner planets | < 20" | Irreducible ephemeris difference |
| Outer planets | < 24" | Irreducible ephemeris difference |

The 4--24" differences for planets arise from the different underlying ephemeris generation and are not correctable.

### Visual Magnitude

LibEphemeris uses **Mallama & Hilton (2018)** formulas for all planets, published in *The Astronomical Journal* and adopted by the Astronomical Almanac. These are the current standard for planetary magnitude computation.

| Body | Formula | Reference | Max diff (measured) |
|------|---------|-----------|---------------------|
| Sun | V(1,0) = -26.86 at 1 AU | Mallama & Hilton 2018 | 0.0000 mag |
| Moon | V = -12.73 + 0.026\|alpha\| + 4e-9\|alpha\|^4 | Astronomical Almanac, Allen's Astrophysical Quantities | 0.03 mag (normal), 0.2 mag (thin crescent) |
| Mercury | 6th-order polynomial in alpha | Mallama & Hilton 2018 | 0.001 mag |
| Venus | Piecewise polynomial | Mallama & Hilton 2018 | 0.0002 mag |
| Mars | Piecewise polynomial | Mallama & Hilton 2018 | 0.0001 mag |
| Jupiter | Quadratic in alpha | Mallama & Hilton 2018 | 0.0001 mag |
| Saturn | Ring-corrected (Meeus geometry) | Mallama & Hilton 2018 | 0.002 mag |
| Uranus | Linear phase coefficient | Mallama & Hilton 2018 | 0.009 mag |
| Neptune | Secular V(1,0) variation | Lockwood & Thompson 1991, Sromovsky et al. 2003 | 0.0000 mag |
| Pluto | Linear phase coefficient | Mallama & Hilton 2018 | 0.04 mag |

#### Neptune secular brightness variation

Neptune's albedo has been increasing since the 1980s due to seasonal atmospheric changes over its 165-year orbital period. LibEphemeris models this with a secular V(1,0) that transitions linearly from -6.89 (pre-1980) to -7.00 (by J2000.0), matching observational data from Lockwood & Thompson (1991) and Sromovsky et al. (2003). This produces exact magnitude agreement across all epochs.

#### Moon thin crescent limitation

For phase angles > 165 degrees (very thin crescents near new moon), the Astronomical Almanac quartic formula diverges from observations. The mean error across all phases is ~0.03 mag; for alpha > 165 degrees it can reach 0.2 mag. This is an inherent limitation of the quartic phase model and affects a narrow range near new moon that is astronomically insignificant (the Moon is unobservable at these phases).

### Apparent Diameter

LibEphemeris uses **IAU 2015 equatorial radii** for all bodies, which is the standard adopted by the Astronomical Almanac for computing apparent angular diameters.

| Body | LibEphemeris radius (km) | Standard |
|------|------------------------|----------|
| Sun | 695,700.0 | IAU 2015 nominal |
| Moon | 1,737.4 | IAU mean |
| Mercury | 2,439.7 | IAU equatorial |
| Venus | 6,051.8 | IAU equatorial |
| Mars | 3,396.2 | IAU equatorial |
| Jupiter | 71,492.0 | IAU equatorial |
| Saturn | 60,268.0 | IAU equatorial |
| Uranus | 25,559.0 | IAU equatorial |
| Neptune | 24,764.0 | IAU equatorial |
| Pluto | 1,188.3 | IAU mean |

The equatorial radius is the correct choice for apparent diameter because it
represents the maximum cross-section as seen from an external observer. A mean
volumetric radius produces systematically smaller diameters (2--3.5% for giant
planets); the IAU equatorial values are the standard used in the Astronomical
Almanac and professional observatory software.

---

## 18. Minor Bodies

### Strict precision mode (default)

For major asteroids with strongly perturbed orbits (Chiron, Ceres, Pallas, Juno, Vesta), LibEphemeris requires SPK kernels by default. Simple Keplerian propagation can produce errors of **1--10 degrees** for these bodies, which is unacceptable for astrological use.

```python
eph.set_auto_spk_download(True)  # Recommended: automatic SPK from JPL Horizons
```

### SPK kernel precision

With SPK kernels (downloaded from JPL Horizons), all minor body calculations achieve **sub-arcsecond** precision matching JPL Horizons exactly.

### Trans-Neptunian Objects

TNOs use Keplerian elements with first-order secular perturbations from Jupiter, Saturn, Uranus, and Neptune (Laplace-Lagrange theory). Typical accuracy: <10° over 50-year spans. For research-grade TNO work, use SPK kernels.

---

## 19. Thread-safe Context API

Each `EphemerisContext` instance has its own isolated state while sharing the expensive ephemeris data across instances. The `_CONTEXT_SWAP_LOCK` (RLock) ensures thread safety during state operations.

Memory overhead per context: ~1 KB. Ephemeris files (DE440: ~119 MB) are loaded once and shared.

---

## 20. Comprehensive Precision Audit

A deep cross-validation audit was performed across all calculation modes and coordinate systems, using `pyswisseph` as a black-box validation oracle. This section summarizes the findings.

### Planetary positions at extreme dates

Tested all 10 major bodies (Sun through Pluto) at dates spanning 1551--2649 CE (the full DE440 range).

| Era | Years | Max dLon | Notes |
|-----|-------|----------|-------|
| Modern | 1900--2100 | <1" | All bodies sub-arcsecond |
| 19th century | 1800--1850 | <1" | All bodies excellent |
| 18th century | 1700 | ~1" | All bodies fine |
| 16th--17th century | 1551--1600 | ~11" | Moon only exceeds threshold |
| 22nd--27th century | 2149--2649 | up to 283" | Moon most sensitive; outer planets remain <3" |

The divergence at extreme dates is driven by ΔT-extrapolation differences. The Moon's rapid motion (~13°/day) amplifies any time-base discrepancy. Outer planets are insensitive due to slow motion and remain within tolerance at all dates.

### Fixed stars

Tested 116 representative stars at 3 epochs (J2000, J2025, J2100). Proper motions updated to van Leeuwen 2007 (new Hipparcos reduction). All stars independently cross-checked against SIMBAD.

| Metric | Value |
|--------|-------|
| Max longitude diff | 0.51" (Rigil Kentaurus — nearest star, parallax not modeled) |
| Mean longitude diff | < 0.1" |
| Stars within 0.5" | 98% of 101 comparable stars |
| Stars within 0.1" | 80% of 101 comparable stars |
| Catalog bugs found & fixed | 2 (Algedi wrong component, Asellus Borealis wrong HIP) |

Remaining sub-arcsecond differences are from precession/nutation pipeline choices. High proper motion stars (Sirius, Rigil Kentaurus) show largest differences due to annual parallax (not modeled).

### House cusps at extreme latitudes

Tested 7 house systems at latitudes 0°--89°, including the arctic circle (66.56°).

| Result | Detail |
|--------|--------|
| Max cusp difference | 0.000005° (0.018") |
| Placidus/Koch above arctic | Both libraries correctly fail (geometrically undefined) |
| Geometric systems (Regiomontanus, Campanus, Equal, Whole Sign, Porphyry) | Valid to 89° with sub-arcsecond agreement |

### Velocities and speed flags

Tested all bodies with FLG_SPEED, FLG_SPEED|FLG_EQUATORIAL, and FLG_SPEED|FLG_NONUT at 20 dates.

| Metric | Value |
|--------|-------|
| Max longitude speed diff | 0.000065 °/day (Moon) |
| Max latitude speed diff | 0.000025 °/day |
| Retrograde sign mismatches | 0 (Mercury station resolved to <10⁻⁶ °/day) |
| Threshold (0.001 °/day) exceeded | Never |

### Heliocentric and barycentric modes

| Mode | Bodies tested | Max dLon | Status |
|------|--------------|----------|--------|
| Heliocentric (FLG_HELCTR) | Mercury--Pluto | 0.00032° | All OK |
| Barycentric (FLG_BARYCTR) | Sun--Pluto | <1" at J2000 | All OK |
| Equatorial (FLG_EQUATORIAL) | Sun, Moon, Mars, Jupiter | 0.00047° (Moon) | All OK |
| XYZ Cartesian (FLG_XYZ) | All bodies | 0.000033 AU (Pluto) | All OK |

Note: Barycentric Sun shows large angular differences (~100"+) at dates when the Sun is near the Solar System Barycenter (distance ~0.001 AU). This is a geometric amplification effect, not a precision issue -- the Cartesian positions agree to <0.000001 AU.

### Sidereal mode and ayanamshas

Tested Lahiri, Raman, Krishnamurti, and Fagan/Bradley ayanamshas.

| Metric | Value |
|--------|-------|
| Max sidereal longitude diff | 0.07" (Moon) |
| Ayanamsha value differences | <0.000001° (floating-point noise) |

### Delta-T

| Era | Max diff (seconds) | Notes |
|-----|-------------------|-------|
| Modern (1900--2025) | <1 sec | Both use IERS observed data |
| Historical (1700--1900) | <1 sec | Model blending differences |
| Pre-telescope (<1700) | up to 187 sec | Different deep-past ΔT reconstructions |
| Future (>2025) | up to 297 sec (at 2500) | Different future parabolic extrapolations |

LibEphemeris uses the Stephenson, Morrison & Hohenkerk (2016) reconstruction. For the modern era where observed data exists, agreement is sub-second.

### Ecliptic obliquity and nutation

| Metric | Max diff |
|--------|----------|
| True obliquity | 0.00086" |
| Mean obliquity | 0.00016" |
| Nutation in longitude | 0.0013" |
| Nutation in obliquity | 0.00087" |

Sub-milliarcsecond agreement across all tested dates (1800--2200), implementing IAU 2006/2000A correctly.

### FLG_MOSEPH behavior

The `FLG_MOSEPH` flag (value 4) is accepted for API compatibility but always
resolves to JPL DE440 — there is no Moshier semi-analytical fallback. Against an
engine that does fall back to Moshier the resulting differences (typically <0.2"
for modern dates) reflect the lower precision of the Moshier ephemeris, not a bug.

---

## See Also

- **[Swiss Ephemeris Comparison](../comparison/index.md)** — the centralized
  head-to-head with pyswisseph: measured precision tables, known differences with
  detailed explanations, intentional divergences, bugs found and fixed, and API
  signature differences.

---

## References

1. Park, R.S. et al. (2021). "The JPL Planetary and Lunar Ephemerides DE440 and DE441." *Astronomical Journal*, 161(3), 105.
2. Capitaine, N. et al. (2003). "Expressions for IAU 2000 precession quantities." *Astronomy & Astrophysics*, 412, 567-586.
3. Mathews, P.M. et al. (2002). "Modeling of nutation-precession: New nutation series for nonrigid Earth." *Journal of Geophysical Research*, 107(B4).
4. Capitaine, N. et al. (2005). "IAU 2006 precession." *Highlights of Astronomy*, 14, 474-475.
5. Lieske, J.H. (1998). "Galilean Satellites of Jupiter. Theory E5." *Astronomy & Astrophysics Supplement*, 129, 205-217.
6. Vienne, A. & Duriez, L. (1995). "TASS1.6: Ephemerides of the major Saturnian satellites." *Astronomy & Astrophysics*, 297, 588-605.
7. Jacobson, R.A. (2009). "The orbits of the Neptunian satellites." *Astronomical Journal*, 137, 4322-4329.
8. Brozovic, M. & Jacobson, R.A. (2024). "The orbits of the Pluto system." *Planetary Science Journal*.
9. Stephenson, F.R., Morrison, L.V. & Hohenkerk, C.Y. (2016). "Measurement of the Earth's rotation: 720 BC to AD 2015." *Proceedings of the Royal Society A*, 472, 20160404.
10. Meeus, J. (1998). *Astronomical Algorithms*, 2nd edition. Willmann-Bell.
11. Moshier, S.L. (1992). "Comparison of a 7000-year lunar ephemeris with analytical theory." *Astronomy & Astrophysics*, 262, 613-616.
12. Schaefer, B.E. (1990). "Telescopic limiting magnitudes." *Publications of the Astronomical Society of the Pacific*, 102, 212-229.
13. IERS Conventions (2010). IERS Technical Note No. 36, ed. Petit, G. & Luzum, B.
14. ESA (1997). *The Hipparcos and Tycho Catalogues*. ESA SP-1200.
15. Mallama, A. & Hilton, J.L. (2018). "Computing Apparent Planetary Magnitudes for The Astronomical Almanac." *Astronomy and Computing*, 25, 10-24.
16. Lockwood, G.W. & Thompson, D.T. (1991). "Solar cycle chromospheric variations and photometric variability of Neptune." *Nature*, 349, 593-594.
17. Sromovsky, L.A. et al. (2003). "Episodic bright and dark spots on Uranus." *Icarus*, 163, 256-261.
18. IAU (2015). "IAU 2015 Resolution B3: Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties." *IAU General Assembly*, Honolulu.
