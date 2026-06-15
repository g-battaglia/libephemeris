# Known Divergences from pyswisseph

This document catalogs all known, inherent divergences between libephemeris and
pyswisseph 2.10.03. These are **not bugs** — they arise from fundamental
differences in the underlying computation engines (Skyfield/JPL DE440 vs Swiss
Ephemeris) and different implementations of secondary algorithms.

All divergences documented here have been verified through systematic
hyper-validation (4400+ comparison rounds across 29 API sections).

## Summary

| Category | Typical Divergence | Maximum Divergence | Cause |
|----------|-------------------|-------------------|-------|
| Planetary positions | 0.001–0.5" | ~1" (Moon) | Different ephemeris engines |
| Planetary speeds | 0.01–5" | ~20" (Moon speed) | Numerical differentiation methods |
| House cusps | <0.01" | ~0.01" | Obliquity/nutation model differences |
| Fixed star positions | <0.01" | <0.01" | Proper motion catalog differences |
| Fixed star distances | <0.01% at J2000 | ~0.1% at ±50y | Radial velocity models |
| Ayanamsha values | <0.1" | ~40" (exotic modes) | Reference star position differences |
| Delta-T | <0.001s (modern) | ~43s (year 1900) | Different delta-T polynomial models |
| Refraction | <1" | ~15" | Different atmospheric models |
| Phase angles | <1" (inner) | ~18" (outer planets) | Position errors amplified |
| Orbital elements | <1" (inner) | ~1° (giant `varpi`) | SE convention; lib matches exact two-body to 0.000000° (see §9) |
| Sidereal Moon | ~3" | ~14" | Lunar theory differences in sidereal frame |

## 1. Ephemeris Engine Differences

### 1.1 Planetary Positions (`calc_ut`)

libephemeris uses **Skyfield + JPL DE440/DE441** while pyswisseph uses the
**Swiss Ephemeris** internal integration. Both are based on JPL ephemerides but
use different interpolation and integration methods.

**Typical divergence:**
- Sun, Mercury–Neptune: 0.001–0.01" (sub-arcsecond)
- Moon: 0.01–1.0" (lunar theory differences)
- Pluto: 0.001–0.01"

**Speed divergence:**
- Most planets: 0.01–2.0"
- Moon speed: up to ~5" (different numerical differentiation)
- With `FLG_SPEED` flag: central finite difference vs analytical

**Future dates (>2050):** Up to 2" position divergence due to delta-T model
extrapolation differences.

**Extreme / BCE dates (the modern-ephemeris fidelity advantage).** Validated at
the *same* TT (so Delta-T cancels) against an exact JPL **DE441** oracle
(jplephem) from ~3000 BCE to +5000 CE. libephemeris (extended LEB tier, fit to
DE441) tracks the exact ephemeris to **< 0.01"** across the whole range for the
Sun, Moon, and inner planets (< 0.09" for the giants, the residual being the
center-vs-barycenter choice). pyswisseph diverges increasingly toward the past —
its geocentric apparent positions vs the same DE441 truth grow to ~11" (Moon, 1000
CE) → ~157" (1000 BCE) → **~550" (~0.15°, Moon, 3000 BCE)**; planets stay smaller
(a few arcsec to ~1'). This is the consequence of the two engines using different
underlying ephemerides (libephemeris: JPL DE441, 2020; pyswisseph: DE431, 2013,
plus its own lunar secular model): DE441 is JPL's latest long-range integration,
and libephemeris reproduces it faithfully over the full -13000 … +17000 span. For
historical/archaeo-astronomy work libephemeris is the more accurate choice.

### 1.2 Interpolated Apogee/Perigee (bodies 21, 22)

IntpApog and IntpPerg use semi-analytical ELP2000-82B perturbation theory.
libephemeris implements this independently from published coefficients, producing
results that can differ by several arcseconds from the Swiss Ephemeris
implementation. These bodies are classified as KNOWN divergence in all tests.

### 1.3 Pholus (body 16) Historical Dates

Pholus SPK coverage in libephemeris starts around 1600 CE. Queries before this
date raise "Invalid Time to evaluate" while pyswisseph may use internal
Keplerian fallback. This affects dates before ~1850.

## 2. House System Differences

### 2.1 House Cusps

House cusps typically agree within 0.01" (sub-arcsecond). The small divergence
comes from slightly different obliquity and nutation models used in the
intermediate calculations.

### 2.2 Vertex at the Equator

At geographic latitude 0° (equator), the Vertex calculation has a 1/tan(lat)
singularity. libephemeris clamps latitude to a tiny epsilon so the formula
evaluates to the correct limiting value, matching pyswisseph. The only remaining
divergence is `ascmc[6]` (CoAsc Munkasey) for the Horizon (H) house system at
lat=0, where pyswisseph returns 0° due to a C `tan(90°)` floating-point
artifact, while the mathematical limit from any positive latitude is 180°.

### 2.3 House Position (`house_pos`)

For most house systems: <0.01" divergence.

For **Alcabitius (B)**, **Koch (K)**, and **Topocentric (T)**: up to ~46"
divergence due to different cusp interpolation algorithms between the engines.

## 3. Time and Delta-T

### 3.1 Delta-T Model

libephemeris uses Skyfield's delta-T model (based on historical observations and
IERS data) while pyswisseph uses the Espenak & Meeus (2006) polynomial model.

| Period | Typical Divergence | Maximum |
|--------|-------------------|---------|
| 1950–2020 | <0.001s | ~0.01s |
| 1900–1950 | ~1s | ~43s |
| 2020–2050 | <0.1s | ~1s |
| >2050 | ~1–10s | depends on extrapolation |

This affects all functions that internally convert between UT and ET/TT:
`deltat`, `deltat_ex`, `utc_to_jd`, `jdet_to_utc`,
`jdut1_to_utc`.

### 3.2 Sidereal Time (`sidtime`)

At modern dates: <0.001s divergence.
At future dates (>2050): up to ~0.05s due to delta-T propagation into the
GMST calculation.

### 3.3 Equation of Time (`time_equ`)

At modern dates: <0.15s (matches well after the GAST-RA rewrite).
At future dates (>2050): up to ~1.4s due to compounding delta-T and Sun RA
model differences.

## 4. Fixed Stars

### 4.1 Positions

Ecliptic longitude and latitude agree within 0.01" for all 116 catalog stars.
This is achieved by using the same FK5/Hipparcos proper motion values.

### 4.2 Distances

At J2000.0 epoch: distances match exactly (same Hipparcos parallax values).
At other dates: up to ~0.1% divergence due to different radial velocity models
affecting the distance variation over time.

### 4.3 Speed in Distance (`speed_dist`)

`speed_dist` (index 5 of the position tuple) diverges significantly (20–90%)
because libephemeris computes it via central finite difference on the full
apparent distance (which includes annual parallax oscillation), while pyswisseph
separates the radial velocity component differently.

### 4.4 Magnitudes

Star magnitudes agree within 0.5 mag for all catalog stars. Small differences
arise from different catalog versions.

### 4.5 Star-name resolution

For a star name absent from this library's catalog, `fixstar2`/`fixstar2_ut`
return an explicit "could not find star name" error rather than a fuzzy-matched
unrelated star (older builds could silently return e.g. `kaRet` for "Chort").

A handful of *traditional* names are deliberately resolved to a different star
than pyswisseph. In each case below libephemeris follows the IAU/WGSN-approved
assignment (or the dominant naming tradition) and pyswisseph uses an older or
non-standard star; mainstream named stars all match to <1".

| Name | libephemeris | pyswisseph | Basis |
|---|---|---|---|
| Menkar | α Cet | λ Cet | IAU name "Menkar" approved for α Ceti (2016) |
| Sadalbari | μ Peg | λ Peg | IAU name "Sadalbari" approved for μ Pegasi |
| Algedi | α² Cap | α¹ Cap | IAU name "Algedi" approved for α² Capricorni |
| Kurhah | ξ Cep | ζ Cep | IAU/WGSN name "Kurhah" approved for ξ Cephei (2016) |
| Mirak | β And (Mirach) | ε Boo | "Mirak" is a documented spelling variant of Mirach (β And); ε Boo is Izar |

Genuinely ambiguous traditional names (no single authority — libephemeris keeps
its existing assignment, pyswisseph differs):

- **Dheneb** — α Cyg (Deneb) here; pyswisseph uses ζ Aql (Deneb el Okab).
- **Ruc** — δ Cas (Ruchbah) here; pyswisseph uses δ Cyg. "Rucbah" is attested
  for δ Cas / α Aqr, not δ Cyg.
- **Girtab** — θ Sco (Sargas) here; "Girtab" is sometimes κ Sco.
- **Deneb Kaitos** — β Cet (Diphda) here; "Deneb Kaitos Shemali" is η Cet.
- **Ukdah** — ι Hya here (IAU-CSN); pyswisseph uses τ² Hya.

The clearly-wrong cases (Alaraph → β Vir, Gienah Corvi → γ Crv, Atri → δ UMa,
Nash → γ² Sgr, Deli → η Aqr) were corrected to match both the reference and the
naming tradition.

### 4.6 Flamsteed designations: `fixstar` vs `fixstar2`

For **Flamsteed-style designations** (e.g. `29Psc`, `2Cet`, `7Cet` — a number
plus a 3-letter constellation), the legacy `fixstar`/`fixstar_ut` and the modern
`fixstar2`/`fixstar2_ut` can resolve the *same name to a different star* (~289 of
the catalog's ~1450 entries, by tens of degrees). This is the inherently fragile
Flamsteed lookup: pyswisseph's own `swe_fixstar` and `swe_fixstar2` likewise
disagree, and resolve these names differently again from libephemeris. **`fixstar2`
is the reliable path** — it matches the requested Flamsteed star correctly (e.g.
`29Psc` → 29 Piscium), where the legacy `fixstar` and both pyswisseph variants
may pick an unrelated row. **Proper names and Bayer designations** (Sirius,
Aldebaran, Regulus, … and all mainstream named stars) resolve **identically** in
`fixstar` and `fixstar2` (0.0"). Recommendation: prefer `fixstar2` for
Flamsteed-designation lookups.

## 5. Refraction and Horizontal Coordinates

### 5.1 Atmospheric Refraction (`refrac`)

Up to 15" divergence, primarily at low altitudes near the horizon. Both
implementations use Bennett's formula but with slightly different coefficients
and boundary handling.

### 5.2 Azimuth/Altitude (`azalt`)

Above-horizon bodies: typically <1" divergence.
Below-horizon bodies: up to ~1654" (~27') divergence due to fundamentally
different atmospheric refraction models for negative apparent altitudes. This
is a known limitation — refraction below the horizon is physically meaningless
and implementations differ in how they extrapolate.

## 6. Eclipse and Occultation Functions

### 6.1 Solar Eclipses

Eclipse timing (`tret[0]` maximum): typically <10s divergence.
Eclipse geography (`geopos`): typically <1° divergence.
Eclipse type flags may occasionally differ for borderline cases.

**Obscuration of a total eclipse (deliberate divergence).** Obscuration is the
fraction of the Sun's *area* covered by the Moon, so it is physically bounded by
1.0 (a total eclipse covers 100% of the Sun). libephemeris returns exactly
**1.0** for total eclipses (`sol_eclipse_how` `attr[2]`,
`sol_eclipse_obscuration_at_loc`, and `sol_eclipse_how_details`
`max_obscuration`). pyswisseph instead reports the lunar/solar disc *area ratio*
`(R_moon/R_sun)² ≈ 1.05–1.12 > 1` for total eclipses, which is not a fraction.
This is an intentional choice for physical correctness; the "how much larger the
Moon appears than the Sun" information remains available in the eclipse
*magnitude* (`attr[0]`/`attr[8]`). Annular eclipses are identical to pyswisseph
(`(R_moon/R_sun)² < 1`, the ring-residual area fraction) and partial eclipses
use the standard two-disc lens overlap; both agree to ~1e-3.

At a no-eclipse instant the obscuration is **0.0** from every entry point —
`sol_eclipse_how` (which zeroes its attr when `retflag == 0`) and now
`sol_eclipse_where` / `lun_occult_where` (which return the closest-approach attr
directly). This is the physically correct fraction; earlier builds leaked a
`1.0` fallback through the `where` path.

### 6.2 Lunar Eclipses

Similar to solar eclipses. Timing agrees within ~10s for most events.

### 6.3 Lunar Occultations

Timing agrees within ~0.001 day (~86s) for most events.
**Note:** `lun_occult_when_glob` is computationally expensive (~6s/call).

## 7. Nodes and Apsides (`nod_aps_ut`)

Nodal and apsidal positions can diverge by 20–700" (up to ~1°) for some bodies.
This is because osculating orbital elements are derived differently from the
two ephemeris engines. The divergence is largest for:
- Outer planets (different integration methods)
- Bodies with high eccentricity (osculating elements more sensitive)

**Lunar True Node / True Lilith (verdict: libephemeris is exact).** Adjudicated
against an independent osculating computation built from the DE441 Moon−Earth
state vector (`Ω = atan2(h_x,−h_y)`, h = r×v; True Lilith = −eccentricity-vector
direction). Over 1975–2022 libephemeris reproduces this exact osculating truth to
**0.000"** (True Node) and **0.000"** (True Lilith), while pyswisseph differs by
~0.13–0.25" (True Lilith) — i.e. libephemeris's osculating lunar apparatus is
numerically exact and the residual is pyswisseph's convention, not a libephemeris
error. (The ~±30° month-scale libration of True Lilith is the real, physical
swing of the instantaneous osculating apse, not noise.)

**True Node *distance* in LEB mode (backend limitation).** The True Node's
*longitude and latitude* are correct in every mode (LEB agrees with the
Skyfield/default backend and pyswisseph to 0.000"). Its **distance** (and hence
`FLG_XYZ` output) is correct only on the **Skyfield/default** backend, which
derives it from the osculating conic of the Moon's DE state and matches pyswisseph
exactly. The **LEB backend** returns a stored proxy (~0.0015 AU vs the correct
~0.0024 AU) rather than recomputing the osculating node radius, because doing so
standalone would require reconstructing the geocentric-ecliptic Moon state inside
the LEB pipeline. Use the default/Skyfield mode if you need the True Node distance;
the angular position is reliable everywhere. (Mean Node and Oscu Apogee distances
agree across backends.)

## 8. Phenomena (`pheno_ut`)

### 8.1 Phase Angle

Inner planets (Sun–Mars): <1" divergence.
Outer planets (Jupiter–Pluto): up to ~18" divergence. The phase angle calculation
amplifies the small position differences between the two engines.

### 8.2 Other Phenomena Values

Phase, elongation, apparent diameter, and magnitude generally agree well (<1").

## 9. Orbital Elements (`get_orbital_elements`)

Orbital elements for **inner planets** (Mercury–Mars) agree within ~1".

Orbital elements for **outer planets** (Jupiter–Neptune) can diverge by
100–2000" in angular elements (argument of perihelion, longitude of ascending
node, mean anomaly). This is because osculating elements are derived from
instantaneous position and velocity vectors, and small position differences
between the engines lead to large differences in the derived elements,
especially for nearly circular orbits where the argument of perihelion is
poorly defined.

**Adjudication vs an exact independent oracle (verdict: libephemeris is correct).**
The longitude of periapsis (`varpi`, element [5]) for the giant planets was
adjudicated against an independent two-body extraction performed on the *same*
DE state vector (jplephem reading DE441, with the library's own obliquity
`23.4392911` and `GM = k²` conventions). Across 1970–2020:

| Quantity | Result |
|---|---|
| max \|lib − independent two-body truth\| | **0.000000°** |
| max \|SE − the same truth\| | up to **1.12°** |

The library reproduces the exact osculating element from the JPL barycentric
state to machine precision; the lib-vs-SE difference is entirely Swiss
Ephemeris's own pipeline convention (state source / light-time / constants),
**not** a libephemeris error. Both engines use the planet *barycenter* for the
giants (the standard osculating-element convention); using the planet *center*
instead would shift `varpi` by up to ~1.6° (Neptune), which is the inherent
center-vs-barycenter ambiguity of a giant-planet osculating element rather than
a defect in either engine.

**Fictitious bodies** (IntpApog=15, IntpPerg=16) have meaningless orbital
elements since they are not real orbiting bodies.

## 10. Sidereal Calculations

### 10.1 Sidereal Positions

Most planets in sidereal mode agree within 1" (the ayanamsha subtraction
is consistent).

**Sidereal Moon:** 3–14" divergence. The sidereal frame amplifies the
underlying lunar theory differences because the ayanamsha correction interacts
with nutation differently in the two engines.

### 10.2 Ayanamsha Values

Standard modes (Lahiri, Fagan-Bradley, Raman, Krishnamurti, and the other
fixed-epoch table systems): **exact** agreement (0.00" at J2000 and across
±100y — the J2000 zero-points and IAU 2006 precession rate are identical).

Exotic/experimental modes diverge, and the divergence **grows with distance
from J2000**: ~0" at J2000 for the star-anchored "True" modes (True Citra,
True Revati, True Pushya) rising to ~40" at ±100y, and up to ~145" for the
galactic/calculated modes (galactic-equator, Aryabhata-class) at the extremes.
These modes are computed dynamically from specific reference-star or
galactic-frame positions, so they inherit the small fixed-star proper-motion /
galactic-frame-definition differences between the two engines (see §4). This is
a definitional difference in those niche modes, not an error in the ayanamsha
machinery: every fixed-epoch mode is exact.

## 11. Crossing Functions

### 11.1 Solar/Lunar Crossings

`solcross_ut`, `mooncross_ut`: typically <1s timing divergence.

### 11.2 Moon Node Crossings

`mooncross_node_ut`: up to ~69s divergence due to different lunar
node calculation methods.

## 12. Asteroid Pipeline (`AST_OFFSET`)

When using `AST_OFFSET + N` to access asteroids by number, libephemeris
remaps to dedicated body IDs and uses its Skyfield/SPK pipeline. pyswisseph
uses `.se1` ephemeris files with a different integration.

Typical divergence: ~0.2" for the major asteroids (Ceres, Pallas, Juno, Vesta).

**Missing .se1 files:** Chiron (2060) and Pholus (5145) require dedicated
`.se1` files in pyswisseph's ephemeris path. If these files are not present,
pyswisseph raises an error. libephemeris always has these bodies available
through its SPK pipeline.

## 13. Heliacal Events

`heliacal_ut` is computationally expensive (>90s per call for some
configurations). Timing divergence: up to ~2–3 days for heliacal rising/setting
events (measured: Venus rising ≈2 d, setting ≈3 d vs pyswisseph), due to
different atmospheric extinction and visibility (arcus-visionis) models. There is
no independent oracle for these events at sub-day precision, so this is recorded
as a model divergence, not adjudicated.

## 13.1 Fictitious / Uranian bodies (Hamburg School, 40–48)

The Uranian planets (Cupido…Poseidon) and Transpluto are propagated from their
published Hamburg-School Keplerian elements. libephemeris vs pyswisseph differ by
up to **~33"** (Cupido; most < 25"), reflecting different element handling — there
is no independent astronomical reference for these fictitious bodies, so neither
is "truth". (The CLAUDE.md "~0.000000\"" figure for Uranians is the LEB-vs-Skyfield
*internal* agreement, not a lib-vs-pyswisseph or lib-vs-reference accuracy.)

## 14. Constants and API

### 14.1 Version String

`version` intentionally differs (libephemeris reports its own version).

### 14.2 `contrib` Attribute

pyswisseph exposes a `contrib` attribute (contributor credits). libephemeris
does not expose this attribute.

### 14.3 `d2l` with Negative Values

`d2l()` for negative input values produces different results due to
unsigned integer overflow behavior in the C implementation of pyswisseph
vs Python's native integer handling.

### 14.4 `FLG_MOSEPH`

`FLG_MOSEPH` is accepted for API compatibility but silently ignored. All
calculations always use JPL DE440/DE441 via Skyfield — there is no Moshier
ephemeris fallback.

## 15. Gauquelin Sectors

Gauquelin sector values typically agree within 0.01 sectors. Small divergences
(<0.1) arise from the underlying planetary position and house cusp differences.

## Hyper-Validation Results

The hyper-validation script (`scripts/hyper_validate.py`) runs 4400+ comparison
rounds across 29 sections. With all tolerance classifications applied:

| Metric | Count | Percentage |
|--------|-------|------------|
| PASS | 3947 | 89.7% |
| KNOWN | 441 | 10.0% |
| FAIL | 0 | 0.0% |
| SKIP | 12 | 0.3% |

**0 failures.** All divergences are documented and classified as inherent
differences between the two computation engines.

The 12 SKIP results are due to missing `.se1` asteroid ephemeris files in the
pyswisseph configuration (not a libephemeris limitation).

Run the hyper-validation yourself:

```bash
.venv/bin/python3 scripts/hyper_validate.py --json report.json
```

Exclude slow sections (occultations, heliacal) with:

```bash
.venv/bin/python3 scripts/hyper_validate.py --section A,B,C,D,E,F,G,H,I,K,L,M,N,O,P,Q,R,S,T,V,W,X,Y,Z,AA,AB,AC
```
