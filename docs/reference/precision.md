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

DE440s is a reduced-size subset of DE440 with no loss of precision within its
range. DE441 is the long-span integration in the same modern ephemeris family,
but it is not numerically identical to DE440 throughout their overlap. The
runtime therefore prefers DE440s/DE440 where they cover a body/date and uses
DE441 to extend the time span.

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

The `FLG_MOSEPH` flag — which in some engines selects a semi-analytical
fallback — is accepted for API compatibility but does not switch LibEphemeris
to that model. Planetary bodies continue through the selected JPL/LEB backend;
standards-derived and analytical bodies retain their documented local models.

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

### Runtime model

Skyfield 1.54 also defaults to the full IAU 2000A series. LibEphemeris uses
pyerfa's IAU 2006-adjusted `nut06a()` realization consistently for planet and
house calculations, keeping nutation aligned with the rest of the ERFA/Vondrák
reduction chain.

### Impact on ecliptic longitude

Nutation in longitude (Δψ) directly shifts ecliptic longitudes. The maximum
amplitude of Δψ is ~17.2 arcseconds (the 18.6-year lunar nodal cycle). A
truncated 77-term IAU 2000B calculation can differ from the full series at the
milliarcsecond scale; it is available only for comparison, not selected by the
runtime pipeline.

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

### JPL center segments and explicit fallback

LibEphemeris uses a physical planet-center segment when the active JPL data
covers the epoch. Otherwise it explicitly uses the planetary system
barycenter; no analytical satellite theory is substituted.

In sealed `leb` mode no planet-center BSP is opened. Jupiter through Pluto are
therefore exactly the system barycentres stored in the selected LEB core, with
`LEB` provenance. This is a source contract, not an implicit precision
fallback. `auto` and `skyfield` retain the JPL center-segment behavior below.

#### Tier 1: SPK-based planet centers (<0.001 arcsec)

An auto-downloaded planet-centers kernel (`planet_centers_{tier}.bsp` per precision tier, legacy name `planet_centers.bsp`; not shipped in the wheel) contains precise center-of-body segments for NAIF IDs 599 (Jupiter), 699 (Saturn), 799 (Uranus), 899 (Neptune), and 999 (Pluto). A regenerated medium build uses the widest published JPL center products:

| NAIF ID | Source SPK | Coverage |
|---------|-----------|----------|
| 599 | jup365.bsp | 1600--2200 |
| 699 | sat441xl (two parts) | covers 1550--2650 |
| 799 | ura111xl-799.bsp | covers 1550--2650 |
| 899 | nep097xl-899.bsp | covers 1550--2650 |
| 999 | plu060.bsp | 1800--2200 |

The `_SpkCenterTarget` class in `planets.py` adds the SPK offset to the DE440 barycenter position at the retarded time (when light left the planet), ensuring the correction is applied correctly for light-time-corrected observations.

If the date falls outside the center segment's coverage, the calculation uses
the DE440/DE441 system barycenter and exposes that choice through tracing. For
geocentric, heliocentric, and barycentric requests, light time is always formed
from the actual observer to the selected target state (physical center when
available, system barycenter otherwise).

### Velocity correction

The general apparent-place pipeline uses a one-second central half-step (six
seconds for the Moon). The LEB fast path uses its documented `0.0005`-day
(`43.2`-second) half-step so the Chebyshev/reduction derivative remains stable:

```
v = [position(t + 1s) - position(t - 1s)] / 2s
v_LEB = [position(t + 43.2s) - position(t - 43.2s)] / 86.4s
```

If that stencil reaches a planet-center coverage boundary, both samples use
the system barycenter. This prevents a finite center/barycenter offset from
being divided by the short stencil and reported as a velocity spike.

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
| Mean node/apogee | One-day centered derivative of ERFA arguments | Smooth standards rate |
| True node/osculating apogee | 0.05 days | Rapid osculating curve |
| Interpolated apogee/perigee | 0.5 days | Complete compatibility series |
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

The of-date mean obliquity is the **Vondrák 2011** pole angle: the angle between the long-term ecliptic pole (`erfa.ltpecl()`) and equator pole (`erfa.ltpequ()`), computed in `libephemeris/sidereal_longterm.py` and exposed via `precession_vondrak.vondrak_mean_obliquity_rad()`. These are the same two poles ERFA builds the precession matrix above from, so precession and obliquity form one frame, valid over ±200,000 years. Near J2000 it matches the IAU 2006 polynomial (Capitaine et al. 2003) to sub-milliarcsecond precision:

```
ε₀ = 84381.406" − 46.836769"·T − 0.0001831"·T² + 0.00200340"·T³
     − 0.000000576"·T⁴ − 0.0000000434"·T⁵
```

where T is Julian centuries from J2000.0.

---

## 8. House Calculations

### Sidereal time

Houses require Greenwich Apparent Sidereal Time (GAST). LibEphemeris takes it
from its own long-term model (`sidereal_longterm.apparent_sidereal_time_deg()`,
see [sidereal-time-longterm.md](../methodology/sidereal-time-longterm.md)):

1. Greenwich Mean Sidereal Time: the IAU 2006 expression — Earth Rotation Angle
   from UT1 plus the precession in right ascension — in the era that expression
   was fitted for, and the long-term geometric construction (Vondrák 2011
   precession, Simon 1994 mean longitude of the Earth) outside it, joined so
   that the two form one continuous curve
2. Equation of equinoxes: Δψ · cos(ε) with the library's own IAU 2006/2000A
   nutation, the same one the planetary positions use
3. GAST = GMST + equation of equinoxes

Precision: sub-milliarcsecond agreement with the IAU 2006 sidereal time in the
modern era, and no polynomial divergence at remote epochs, where an IAU-2006
sidereal time is wrong by degrees. Only the ΔT model choice remains (§ Delta T).

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

The True Node is where the Moon's **instantaneous orbital plane** intersects the ecliptic. LibEphemeris derives it geometrically from the active JPL ephemeris:

1. Get Moon's geocentric position **r** and velocity **v** from JPL DE440 (~1 milliarcsecond source precision)
2. Compute angular momentum vector: **h** = **r** × **v**
3. This vector is perpendicular to the instantaneous orbital plane
4. Find intersection with the ecliptic: Ω = atan2(h_x, −h_y)
5. Transform into the true ecliptic frame of date using IAU 2006 precession and IAU 2000A nutation

This computes **exactly** what the True Node is by definition: the intersection of the orbital plane with the ecliptic. Deriving the node directly from the state vectors avoids the approximation errors inherent in analytical series — the geometric construction is exact by definition.

No analytical perturbation series is layered over the state-vector result: the
JPL state already contains the physical perturbations. The residual against an
analytical mean-node series reflects the difference between instantaneous and
mean orbital planes. See [Precision History](../development/precision-history.md)
for the investigation record.

### Mean Node

Evaluated from the IERS 2003 Delaunay argument for the mean longitude of the
lunar ascending node, through ERFA's standards implementation.
The reference radius is the conventional mean lunar distance of 384,400 km,
converted with the exact IAU astronomical unit.

### Mean Lilith (Black Moon)

Evaluated from the IERS Delaunay identities for the mean lunar apse on a
conventional `5.145°` lunar orbital plane. Distance uses the conventional
mean-apogee radius. Rates are centered derivatives of the complete analytical
state.

### True Lilith (osculating apogee)

Computed from JPL DE440/DE441 Moon state vectors via orbital mechanics: the
eccentricity vector **e** is derived from position and velocity, and the
apogee direction is opposite **e**.

### Interpolated Apogee and Perigee

Evaluated as separate Delaunay-argument series for apogee and perigee fitted
at the actual DE440 apsis passages, plus spline-interpolated residual tables
that make the curves exact at every passage. The model file
(`lunar_apse_model.py`) is generated by the committed
`scripts/generate_lunar_apse_model.py` and accepted only at its reviewed
SHA-256; modified payloads fail the lunar provenance gate. Latitude follows a
two-harmonic DE440 fit (within `0.25°` of the instantaneous passage
latitude). The distance channels are the mean DE440 passage distances
(`0.0027100 AU` for apogee and `0.0024236 AU` for perigee), not instantaneous
lunar distances.

### Verification

The lunar model tests cover the ERFA/IERS fundamental arguments,
angular-momentum and eccentricity-vector mechanics for osculating points,
finite-difference agreement of reported velocities, exact artifact hashing,
and continuity and rate checks on the interpolated curves. Apogee and
perigee are intentionally evaluated by different series and are not required
to be antipodal.

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

All predefined base IDs 0--46 compute again. Formula/epoch modes share one
defining table across the direct and LEB backends; stellar and galactic modes
use the independent catalogue pipeline at runtime. The source audit is not
uniformly complete: [Sidereal modes](ayanamsha.md) distinguishes definitions
verified in primary publications from modes whose defining epochs are
documented but whose full primary-source chain is still being traced.

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

Plain `refrac()` uses the published Sæmundsson (true-to-apparent) and Bennett
(apparent-to-true) closed forms (Meeus 1998, ch. 16), including the documented
branch and clamp behavior of the public API. The physical ICAO ray tracer
remains available internally. Pressure and temperature are configurable in both
directions.

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
3. **Arcus visionis**: root-solved over the limiting-magnitude model (no
   Ptolemaic coefficient table)
4. **Contrast threshold**: Schaefer model with configurable observer skill levels

Timing precision: <1 day for heliacal rising/setting events.

---

## 16. Planetary Positions -- Measured Precision

Angular agreement with an independent reference is graded by the project
validation suite, which compares public outputs body-by-body against JPL
DE440/DE441 states, ERFA/SOFA reductions, and JPL Horizons — each within a
documented per-body envelope — and currently reports no result outside those
envelopes across the modern range. The classes below summarise that sampled
comparison; they are not presented as a continuous maximum, and no per-date
reference output is retained.

### Longitude and latitude

| Body group | Modern-era agreement class |
|------------|----------------------------|
| Sun, Mercury–Pluto (major planets) | sub-arcsecond, at the hundredths-of-an-arcsecond class (observed a few ×0.01 arcsec) |
| Moon | largest of the major bodies; sub-arcsecond in the modern range, set by the nutation-model and COB-correction pipeline choices |
| All bodies toward the JPL kernel endpoints | agreement widens away from J2000, most visibly for the Moon and at the extended-tier edges |

The angular class is set by intentional model choices (nutation realization,
body-center convention, ΔT model), not by solver error; the contributing
families are broken down below.

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
| Planet-center availability | JPL physical center or explicit system barycenter |
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

**Module-level event searches mutate the global observer.** `rise_trans`,
`rise_trans_true_hor`, and the location-dependent eclipse searches
(`sol_eclipse_when_loc`, `lun_eclipse_when_loc`, and related chains) set the
process-global topocentric observer for the duration of the search, so the
underlying `calc_ut(FLG_TOPOCTR)` chain sees the requested geographic
position; the previous observer is restored on exit, including on error. This
is deliberately not locked — a lock around the observer alone would only
feign safety while the rest of the module state remains shared. A thread that
calls `calc_ut` with `FLG_TOPOCTR` while another thread runs one of these
searches observes the search's observer. Threads that need concurrent
topocentric isolation should run their calculations through their own
`EphemerisContext`.

---

## 20. Comprehensive Precision Audit

A comprehensive audit covers all calculation modes, coordinate systems, flags,
and body families. Accuracy acceptance is based on independent sources:

- JPL DE440/DE441 and Horizons for planetary, lunar, and minor-body states;
- ERFA/SOFA and IERS for frames, Earth orientation, and time scales;
- Hipparcos/Gaia/SIMBAD and Astropy for fixed-star astrometry;
- published geometric definitions for houses, nodes, apsides, eclipses,
  occultations, refraction, and visibility;
- backend equivalence, defining conditions, round trips, step halving, and
  metamorphic relations for numerical correctness.

Behavioral compatibility with the reference API is assessed on the public call
surface; the compatibility overview lives in the comparison section.

The main bounded difference families are remote-epoch Earth rotation,
different JPL ephemeris solutions, abstract lunar-point definitions,
catalog/space-motion policies, Keplerian minor-body fallback, atmosphere and
observer models, grazing event thresholds, and historical ayanamsha
definitions.

`FLG_MOSEPH` is accepted for API compatibility but still resolves to the JPL
pipeline; LibEphemeris has no Moshier semi-analytical fallback.

---

## See Also

- **[Swiss Ephemeris Comparison](../comparison/index.md)** — compatibility
  methodology, aggregate behavior, known differences, intentional divergences,
  and API signature differences.

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
