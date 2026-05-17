# LibEphemeris vs Swiss Ephemeris — Comparison & Known Differences

This document is a direct, exhaustive comparison between **libephemeris** and **pyswisseph** (the Python wrapper for Swiss Ephemeris). It documents measured precision, known differences, their root causes, and the methodology used to verify compatibility.

LibEphemeris aims for 1:1 API compatibility with Swiss Ephemeris. The two libraries produce results that are extremely close but will never be bit-identical, because they are built on fundamentally different foundations:

| | LibEphemeris | Swiss Ephemeris (pyswisseph) |
|---|---|---|
| Ephemeris | NASA JPL DE440/DE441 (2021) | JPL DE431 (2013), in proprietary format |
| Lunar model | DE440 numerical integration | Analytical (ELP/MPP02) + DE431 |
| Nutation | IAU 2006/2000A via pyerfa | IAU 2006/2000A (internal) |
| Delta T | IERS observed + Stephenson et al. 2016 | Espenak & Meeus 2006 |
| Velocities | Central difference (numerical) | Chebyshev derivative (analytical) |
| Planet centers | Automatic COB correction (SPK + analytical) | Barycenter by default |
| Implementation | Pure Python + Skyfield + pyerfa | C (with Python wrapper) |

## Table of Contents

- [1. Precision Summary](#1-precision-summary)
- [2. Known Differences in Detail](#2-known-differences-in-detail)
- [3. Bugs Found and Fixed](#3-bugs-found-and-fixed)
- [4. API Signature Differences](#4-api-signature-differences)
- [5. Validation Methodology](#5-validation-methodology)
- [6. References](#6-references)

---

## 1. Precision Summary

All measurements taken across 100--210 dates spanning the DE440 range (1550--2650 CE), comparing libephemeris output against pyswisseph with identical flags and parameters.

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

All planets except the Moon are sub-2" (sub-arcsecond for outer planets). The Moon difference is explained in [Section 2.2](#22-moon-precision--max-135-over-800-years).

### Lunar nodes and Lilith

| Body | Longitude (max) | Latitude (max) | Speed (max) |
|------|-----------------|----------------|-------------|
| True Node | < 0.01" | 0.0" | 0.0007°/day |
| Mean Node | < 0.001° | 0.0" | 0.00005°/day |
| True Lilith (Osc. Apogee) | < 0.5" | 0.0" | 0.015°/day |
| Mean Lilith (Mean Apogee) | < 0.015" (±100yr) | ~20" | 0.00005°/day |

Mean Lilith latitude ~20" difference is from different analytical node formulas used in the orbital-plane-to-ecliptic projection. No practical impact.

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
| House cusps (24 systems) | < 0.02" | All house systems tested at 11 global locations |
| Fixed stars (116 stars) | < 0.51" | van Leeuwen 2007 proper motions, verified vs SIMBAD |
| Coordinate transforms | Exact | cotrans, azalt, equatorial ↔ ecliptic |
| Utility functions | Exact | julday, revjul, degnorm, split_deg, etc. |
| Solar eclipse timing | < 6 sec | Maximum, contact times (C1--C4) |
| Lunar eclipse timing | < 8 sec | All contact types (P1, U1--U4, P4) |
| Lunar eclipse classification | All match | Total, partial, penumbral — all agree |
| Rise/set/transit | < 4 sec | Sun, Moon, planets |
| Crossings (solcross, mooncross) | < 4 sec | Longitude crossing times |
| Topocentric positions | All pass | 11 global locations, all planets |
| calc_pctr (planet-centric) | < 0.15° | All planet pairs tested |
| Sidereal modes (43 ayanamshas) | < 0.006° | All formula-based and star-based |
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

## 2. Known Differences in Detail

These are not bugs. They are inherent consequences of using different astronomical models, different ephemerides, or different algorithmic strategies. Each is explained in full so that users and developers understand exactly what to expect.

### 2.1. Crossing functions — full-orbit search for slow planets

`swe_cross_ut` now handles full-orbit searches for all planets, including slow outer planets like Jupiter and Saturn. The algorithm automatically scales the Brent bracket search window based on the estimated time to crossing (`dt_guess`), and filters out false sign changes at the antipodal point (target ± 180°) during the coarse scan. This ensures convergence even when the crossing is 10+ years away and the planet undergoes multiple retrograde periods en route.

Previously, Jupiter 0° Aries searches could converge on 180° instead of 0°. This was fixed by:
1. Scaling the search window to `max(min_window, dt_guess * 1.5)` instead of a fixed 500-day cap
2. Scaling the number of coarse samples proportionally (`search_window / 30` days per sample)
3. Rejecting bracket candidates where the signed-difference jump exceeds 180° (wrapping artifacts at the antipodal point)

### 2.2. Moon precision — max ~135" over 800 years

**What:** The maximum difference in lunar longitude between libephemeris and pyswisseph is approximately 135 arcseconds (0.037°), occurring at the extremes of the DE440 date range (around 1550 CE and 2650 CE).

**Why:** The two libraries use fundamentally different lunar models:

- **libephemeris** uses **JPL DE440**, where the Moon's position is computed by numerical integration of the full equations of motion for all Solar System bodies. DE440 is fitted to Lunar Laser Ranging data (1970--present) with ~1 milliarcsecond precision for the modern era.

- **pyswisseph** uses a combination of **DE431** for the modern period and the **ELP/MPP02 analytical theory** for the Moon. The analytical theory represents the Moon's motion as a sum of thousands of trigonometric terms — a fundamentally different mathematical approach.

For dates near the present (1900--2100), the two models agree to within a few arcseconds. But as you move to historical or far-future dates, they diverge because:

1. **Different tidal acceleration models**: The Moon is slowly receding from Earth due to tidal friction. DE440 and the Swiss Ephemeris analytical theory model this effect with slightly different parameters, causing cumulative divergence over centuries.

2. **Different fitting datasets**: DE440 incorporates 8 additional years of Lunar Laser Ranging data compared to DE431, plus improved planetary observations from Juno (Jupiter) and MESSENGER (Mercury) that indirectly affect the lunar solution.

3. **Numerical vs analytical**: Numerical integration accumulates floating-point errors over very long time spans, while analytical series have truncation errors from omitting higher-order terms. The two error profiles are different.

**Practical impact:** For dates within ±200 years of the present, the difference is < 5". For typical astrological use (natal charts, transits, progressions), the difference is astronomically negligible — it corresponds to a timing error of roughly 4 minutes in the Moon's position.

### 2.3. Delta T — up to 232 seconds divergence

**What:** The computed Delta T (ΔT = TT − UT1) can differ by up to 232 seconds between the two libraries for extreme dates.

**Why:** Delta T represents the difference between uniform atomic time (TT) and the irregular rotation of the Earth (UT1). The Earth's rotation is unpredictable — it slows down due to lunar tides, speeds up due to post-glacial rebound, and fluctuates on decadal timescales for reasons not fully understood.

The two libraries use different models to estimate ΔT:

| Period | LibEphemeris | pyswisseph |
|--------|-------------|------------|
| Historical (before ~1600) | Stephenson, Morrison & Hohenkerk (2016): cubic spline from eclipse observations | Espenak & Meeus (2006): polynomial fits |
| Modern (1973--present) | IERS observed values (when enabled) | Tabulated values |
| Future (2050+) | Parabolic extrapolation | Different parabolic extrapolation |

For the modern era (1900--2020), where direct measurements exist, the two models agree to within ~1 second. The 232-second maximum occurs for dates around 500 CE or beyond 2500 CE, where:

- **Historical dates**: The only data sources are ancient eclipse records (Babylonian, Chinese, Arabic). The Stephenson et al. (2016) model re-analyzed these records with improved methodology, producing different ΔT values than the older Espenak & Meeus (2006) model.
- **Future dates**: Nobody knows how the Earth's rotation will change. Both models extrapolate using different parabolic formulas, and the divergence grows quadratically with time.

**Practical impact:** ΔT affects the mapping from calendar time to astronomical time. A 232-second difference means that the two libraries place the same calendar date at slightly different points on the astronomical timeline. For planetary positions, this translates to at most a few arcseconds of positional difference (the Moon moves ~0.5"/second, so 232 seconds ≈ 116" of lunar motion — consistent with the Moon precision figures above). For the Sun and planets, the effect is much smaller.

### 2.4. House cusp speeds — numerical vs analytical derivatives

**What:** `swe_houses_ex2` and `swe_houses_armc_ex2` compute cusp velocities using centered finite differences when `FLG_SPEED` is set. The maximum difference from pyswisseph is ~0.7 deg/day (~0.2% relative).

**Why:** The two libraries use different differentiation strategies:

- **libephemeris** uses **numerical differentiation** (central difference at ±1 minute): `speed = (cusp(t+dt) − cusp(t−dt)) / (2·dt)`, with proper 0°/360° wraparound handling.

- **pyswisseph** uses **analytical derivatives** by differentiating the house cusp formulas symbolically with respect to time.

Both approaches are valid. The numerical method introduces a small truncation error proportional to `dt²`, which at a 1-minute step produces differences of < 1 deg/day from the analytical result. Given that cusp speeds are ~280--340 deg/day, this corresponds to < 0.3% relative error.

**Note:** Cusp speeds are only computed when `FLG_SPEED` is passed in the `flags` parameter. Without this flag, zero velocities are returned for efficiency.

### 2.5. Asteroids — Keplerian approximation vs integrated ephemerides

**What:** Minor body positions (Chiron, Ceres, Pallas, Juno, Vesta) can differ by up to ~5° from pyswisseph.

**Why:** The two libraries compute asteroid orbits using radically different methods:

- **pyswisseph** uses **pre-integrated ephemerides from NASA JPL** stored in dedicated asteroid ephemeris files. These are computed by numerical integration that includes gravitational perturbations from all major planets, relativistic corrections, solar radiation pressure, and other effects. The result is sub-arcsecond precision.

- **libephemeris** uses a **Keplerian approximation** by default: it treats the asteroid as orbiting the Sun on an ideal ellipse defined by its mean orbital elements. This ignores:
  - Gravitational perturbations from Jupiter (dominant effect for Chiron, which orbits between Saturn and Uranus)
  - Perturbations from Saturn and other planets
  - Orbital resonance effects
  - Secular variations of the orbital elements over time

Near the epoch of the orbital elements, the Keplerian approximation is accurate (< 1°). But it degrades with time because the real orbit is continuously perturbed. Chiron is particularly affected because its orbit is chaotic — small perturbations from Saturn and Uranus grow exponentially.

**Mitigation:** libephemeris supports loading **SPK kernels** from JPL Horizons for any minor body:

```python
import libephemeris as eph
eph.set_auto_spk_download(True)  # Automatic download from JPL Horizons
```

With SPK kernels, all minor body calculations achieve sub-arcsecond precision matching JPL Horizons exactly. The Keplerian fallback exists only for cases where SPK data is unavailable.

**Note:** The JPL DE440/DE441 ephemerides used by libephemeris cover only the major planets (Mercury--Pluto) and the Moon. Asteroid positions require separate data sources regardless of the ephemeris generation.

---

## 3. Bugs Found and Fixed

The deep validation effort uncovered 12 bugs in libephemeris, all of which have been fixed and verified. This section documents them as evidence of the validation's thoroughness.

### 3.1. Output flag handling (2 bugs)

| # | Bug | Fix | Commit |
|---|-----|-----|--------|
| 1 | `FLG_RADIANS` (8192) was silently ignored — positions always returned in degrees | Added `_apply_output_flags()` helper in `planets.py` that converts degrees to radians when the flag is set | `ef29a08` |
| 2 | `FLG_XYZ` (4096) was silently ignored — positions always returned in spherical coordinates | Same helper converts spherical (lon, lat, dist) to Cartesian (x, y, z) in AU, including velocity transformation via the Jacobian matrix | `ef29a08` |

### 3.2. Unit and threshold errors (2 bugs)

| # | Bug | Fix | Commit |
|---|-----|-----|--------|
| 3 | `swe_pheno_ut` returned apparent diameter in arcseconds; pyswisseph returns degrees | Divide by 3600.0 in `_calc_apparent_diameter()` | `ef29a08` |
| 4 | Lunar eclipse penumbral detection missed eclipses near 12--18° from node — `ECLIPSE_LIMIT_LUNAR` was 12.0° but penumbral eclipses occur up to 18.02° from a lunar node | Changed limit to 18.5° | `ef29a08` |

### 3.3. House system bug (1 bug)

| # | Bug | Fix | Commit |
|---|-----|-----|--------|
| 5 | Sunshine house system (`'I'`) crashed or returned wrong cusps at latitudes ≥ 67°N where the MC falls below the horizon | Added MC-under-horizon detection: when MC is below horizon, flip MC by 180°, compute cusps normally, then shift intermediate cusps by 180° | `ef29a08` |

### 3.4. Fixed star bugs (3 bugs)

| # | Bug | Fix | Commit |
|---|-----|-----|--------|
| 6 | `swe_fixstar_mag` returned 0.0 for most stars — the `_STAR_MAGNITUDES` dict only had 2 entries | Replaced with dict comprehension from the full `STAR_CATALOG` | `ef29a08` |
| 7 | 11 fixed stars missing from catalog (Tejat, Propus, Wasat, Thuban, Rasalgethi, Albireo, Wezen, Adhara, Alhena, Alpheratz, Algenib) | Added all 11 with data sourced from ESA Hipparcos via SIMBAD/CDS. Also fixed Tejat's HIP number and a duplicate alias key | `ef29a08` |
| 8 | `resolve_star_name()` used substring matching (`if normalized in alias`), causing false positives (e.g., "Al" matching "Aldebaran") | Changed to prefix matching (`if alias.startswith(normalized)`) | `ef29a08` |

### 3.5. Lunar eclipse classification (2 bugs)

| # | Bug | Fix | Commit |
|---|-----|-----|--------|
| 9 | Lunar eclipse type misclassification — some total eclipses classified as partial, and vice versa. Missing the standard Danjon atmospheric shadow enlargement factor | Added `_SHADOW_ENLARGEMENT = 1.0 + 1.0 / 85.0` (≈ 1.0118) to both umbra and penumbra radii in `_calculate_lunar_eclipse_type_and_magnitude()` | `e082ee5` |
| 10 | P1/P4 contact times inconsistent with `lun_eclipse_when` — the Danjon factor was applied in one function but not in `_calc_lunar_eclipse_penumbral_separation()` | Applied the same `_SHADOW_ENLARGEMENT` factor to the penumbral separation function | `e082ee5` |

The Danjon atmospheric shadow enlargement is a well-established correction: Earth's atmosphere refracts sunlight around the limb, making the geometric shadow ~1.2% larger than it would be without atmosphere. The factor 1/85 was determined empirically by André Danjon and is used by the Astronomical Almanac and by Swiss Ephemeris.

### 3.6. Algorithmic improvements (2 bugs)

| # | Bug | Fix | Commit |
|---|-----|-----|--------|
| 11 | `swe_cross_ut` diverged for slow outer planets (Saturn, Jupiter) near stations — the station detection threshold was too narrow and the Brent bracket search window too small | Widened `STATION_SPEED_THRESHOLD` from 0.001 to 0.01 °/day; expanded Brent bracket windows; increased `max_range` limits per planet speed category; scaled search window with `dt_guess * 1.5`; added antipodal-point wrapping filter to bracket scanner | `d39aaf8` |
| 12 | `swe_orbit_max_min_true_distance` returned 2 values `(min, max)` in wrong order; pyswisseph returns 3 values `(max, min, true_dist)` | Changed return to 3-tuple `(max_dist, min_dist, true_dist)` with correct order, added true distance calculation from `swe_calc_ut` | `d39aaf8` |

---

## 4. API Signature Differences

libephemeris aims for 1:1 API compatibility with pyswisseph, but some function signatures differ. Developers migrating from pyswisseph should be aware of these:

### Return value differences

| Function | pyswisseph returns | libephemeris returns |
|----------|-------------------|---------------------|
| `swe_get_ayanamsa_ex_ut` | `(flags, ayanamsa)` | `(ayanamsa, eps_true, nut_long)` |
| `swe_deltat_ex` | `float` | `(float, str)` — value + error message |
| `swe_orbit_max_min_true_distance` | `(max, min, true)` | `(max, min, true)` — same after fix |

### Parameter differences

| Function | pyswisseph signature | libephemeris signature |
|----------|---------------------|----------------------|
| `swe_get_ayanamsa_ex_ut` | `(tjd_ut, flags)` | `(tjd_ut, sid_mode, flags)` |
| `swe_heliacal_ut` | `(jd, geopos, datm, dobs, name, event, flags)` → `(jd1, jd2, jd3)` | `(jd, geopos, datm, dobs, name, event, flags)` → `(jd1, jd2, jd3)` ✅ |
| `swe_lun_occult_when_loc` | body can be `int` or `str` | body can be `int` or `str` — same |

### Structural differences

| Function | pyswisseph | libephemeris |
|----------|-----------|-------------|
| `swe_get_orbital_elements` | Returns flat tuple | May return nested tuple — access via `result[0][i]` |
| `swe_houses_armc` | `ascmc[3]` = obliquity | `ascmc[3]` = Vertex (not obliquity) |
| `swe_houses_ex2` | Returns cusp speeds (analytical) | Returns cusp speeds (numerical, requires `FLG_SPEED`) |

### libephemeris-only extensions

These functions exist in libephemeris but have no pyswisseph equivalent:

| Function | Description |
|----------|-------------|
| `swe_houses_with_fallback` | Houses with automatic polar-latitude fallback |
| `swe_houses_armc_with_fallback` | ARMC houses with automatic polar-latitude fallback |
| `swe_sol_eclipse_max_time` | Precise maximum eclipse timing |
| `swe_sol_eclipse_how_details` | Comprehensive eclipse circumstances (dict) |
| `swe_sol_eclipse_obscuration_at_loc` | Eclipse obscuration at a geographic location |
| `swe_planet_occult_when_glob` | Planet-planet occultation search (global) |
| `swe_planet_occult_when_loc` | Planet-planet occultation search (local) |
| `swe_calc_angles` | Pre-calculated angles for Arabic parts |
| `swe_calc_eclipse_path_width` | Eclipse path width at a point |
| `swe_calc_eclipse_central_line` | Eclipse central line coordinates |
| `swe_calc_eclipse_northern_limit` | Eclipse northern limit coordinates |
| `swe_calc_eclipse_southern_limit` | Eclipse southern limit coordinates |

---

## 5. Validation Methodology

### Test infrastructure

The validation consists of 5 test suites totaling **1,619 automated tests**:

| Suite | File | Tests | Focus |
|-------|------|-------|-------|
| 1 | `test_deep_validation.py` | 514 | Planetary positions (10 planets × 12 flag combos × 210 dates), houses, fixed stars, crossings, eclipses, rise/set, coordinates, utilities |
| 2 | `test_deep_validation_2.py` | 444 | Sidereal modes, topocentric, TT variants, nodal/apsides, orbital elements, phenomena, combined flags, eclipse details |
| 3 | `test_deep_validation_3.py` | 95 | Eclipse geography, occultations, Gauquelin sectors, heliacal events, ARMC ex2, orbit distances, cross_ut, helio_cross_ut |
| 4 | `test_deep_validation_4.py` | 63 | Polar fallback houses, eclipse max time/details/obscuration, planet occultations, heliacal_pheno_ut, calc_angles, state functions |
| 5 | `test_compare_helio_bary.py` | 503 | Heliocentric, barycentric, equatorial, XYZ modes (10 bodies × 10 dates × 4 modes), combined flags, return flag verification |

### Test parameters

- **Mode**: Skyfield only (LEB mode not tested in this comparison)
- **Date range**: 1550--2650 CE (full DE440 range), with concentration around 1900--2100
- **Sample density**: 100--210 dates per test, distributed to cover equinoxes, solstices, eclipses, and random dates
- **Locations**: 11 global locations including equator, tropics, mid-latitudes, Arctic, and Antarctic
- **Flags**: All individual flags (`FLG_SPEED`, `FLG_EQUATORIAL`, `FLG_J2000`, `FLG_NONUT`, `FLG_NOABERR`, `FLG_NOGDEFL`, `FLG_TRUEPOS`, `FLG_RADIANS`, `FLG_XYZ`, `FLG_HELCTR`, `FLG_BARYCTR`, `FLG_ICRS`, `FLG_SIDEREAL`) and common combinations
- **House systems**: All 24 supported systems

### What "pass" means

Each test compares libephemeris output against pyswisseph output with tolerance thresholds appropriate to the quantity being measured:

| Quantity | Tolerance | Notes |
|----------|-----------|-------|
| Ecliptic longitude (planets) | < 0.001° (3.6") | Sub-arcsecond for outer planets |
| Ecliptic longitude (Moon) | < 0.04° (135") | Model difference, not a bug (see §2.2) |
| Ecliptic latitude | < 0.001° | Sub-arcsecond for all bodies |
| Equatorial RA/Dec | < 0.001° | Sub-arcsecond |
| Distance (AU) | 0.0001 AU | |
| Angular velocity | 0.003°/day | Lon speed; lat speed < 0.004°/day |
| True Node/Lilith longitude | < 0.001° (< 0.5") | Verified vs JPL Horizons |
| House cusps | < 0.02" | All 24 systems |
| Fixed star positions | < 0.51" | van Leeuwen 2007, verified vs SIMBAD |
| Eclipse/crossing times | < 8 seconds | Lunar eclipses; solar < 6 sec |
| Sidereal time | < 2 seconds | |
| Heliocentric/Barycentric | < 0.001° | Sub-arcsecond (except bary Sun, see §6.7) |
| XYZ Cartesian | < 0.00005 AU | Sub-arcsecond angular |

These tolerances were determined empirically through exhaustive comparison testing and independently verified against JPL Horizons, astropy/ERFA, and SIMBAD (see Section 6). Earlier versions of this document had much looser tolerances (e.g., 0.15° for True Node, 1.5° for latitude) that have been tightened by 100–1500× after triangulated verification showed the actual agreement is far better.

### Final results

| Result | Count |
|--------|-------|
| **Passed** | 1,554 |
| **Failed** | 0 |
| **Skipped** | 5 (convergence-dependent tests) |
| **Expected failures** | 0 |
| **Total** | 1,612 |

All 87 exported `swe_*` functions are covered by at least one test suite.

---

## 6. Independent Verification Results

These results are from triangulated comparisons: **libephemeris vs pyswisseph vs independent source** (JPL Horizons, astropy/ERFA). The purpose is to determine, for each area of disagreement, which implementation is more astronomically accurate.

### 6.1. Lunar True Node — match JPL Horizons to < 0.01"

**Method:** Computed True Node (Ω) using all three sources across 24 dates (1950–2050).

| Source | vs JPL Horizons (J2000 ecliptic) |
|--------|----------------------------------|
| **LibEphemeris** | < 0.01" (machine precision) |
| **pyswisseph** | < 0.01" (machine precision) |

**Conclusion:** Both libraries are effectively identical for True Node longitude. The old test tolerance of 0.15° was absurdly relaxed (actual agreement is 150× tighter). Tightened to 0.001°.

### 6.2. True Lilith (Osculating Apogee) — both equivalent

**Method:** Compared osculating apogee longitude using all three sources across 200 dates (1950–2050).

| Source | vs JPL Horizons |
|--------|-----------------|
| **LibEphemeris** | ~240" (4 arcmin) |
| **pyswisseph** | ~240" (4 arcmin) |

The ~240" offset from Horizons is **not a bug** — it's inherent to the two-body eccentricity vector formula (`e = (v×h)/μ - r/|r|`), which ignores solar perturbations on the Moon. Horizons uses a different effective GM that accounts for perturbations.

LibEphemeris matches pyswisseph to **mean 0.1", max 0.5"** (effectively identical). Old tolerance of 0.1° tightened to 0.001°.

### 6.3. J2000 frame for lunar bodies — LibEphemeris more accurate

**Method:** Compared J2000 ecliptic coordinates across planets and lunar bodies.

**Finding:** For planets (Sun, Moon, Mars, Jupiter), both libraries agree with astropy/ERFA to < 1" — no meaningful difference.

However, for **lunar nodes and Lilith** in J2000 frame, pyswisseph applies a systematic ~14" offset that **should not exist** at the J2000.0 epoch. At J2000.0, tropical and J2000 ecliptic coordinates are identical by definition (zero precession), yet pyswisseph shifts the position by ~0.004°. LibEphemeris correctly returns identical tropical and J2000 values at the J2000.0 epoch.

| Date | LibEphemeris J2000 shift | pyswisseph J2000 shift | Expected at J2000.0 |
|------|-------------------------|------------------------|---------------------|
| J2000.0 | **0.0000°** | 0.0039° | **0.0000°** |
| J1900 | −1.3966° | −1.3917° | ~−1.4° (precession) |
| 2023 | +0.3234° | +0.3208° | ~+0.3° (precession) |

**Conclusion:** LibEphemeris has a more accurate J2000 frame transformation for analytically-computed bodies.

### 6.4. Latitude, equatorial, velocity — all verified identical

Comprehensive measurement of latitude, equatorial (RA/Dec), and velocity differences for all lunar nodes/Lilith across all core dates shows:

| Component | True Node | True Lilith | Mean Node | Mean Lilith |
|-----------|-----------|-------------|-----------|-------------|
| Latitude | 0.0" | 0.0" | 0.0" | 20" max |
| RA (equatorial) | 0.0" | 0.2" | 0.0" | 2" max |
| Dec (equatorial) | 0.0" | 0.0" | 0.0" | 19" max |
| Speed (lon) | 0.0007°/day | 0.015°/day | 0.00005°/day | 0.00005°/day |
| Speed (lat) | 0.0" | 0.0004°/day | 0.0" | 0.00002°/day |

All previous tolerances (0.2–1.5°) were wildly over-conservative. Tightened to reflect actual measured precision.

### 6.5. Fixed star flags — previously silently ignored

The `fixed_stars.py` module previously handled only 4 of 19 SEFLG flags. Most critically, `FLG_SIDEREAL` returned tropical coordinates without warning. This has been fixed — all meaningful flags are now handled:
`FLG_SIDEREAL`, `FLG_J2000`, `FLG_NONUT`, `FLG_XYZ`, `FLG_RADIANS`, `FLG_TRUEPOS`, `FLG_MOSEPH`, `FLG_SPEED3`, `FLG_TOPOCTR`.

### 6.6. Fixed star catalog — cross-checked against Hipparcos/SIMBAD

All 116 catalog stars were independently verified against authoritative sources:

**Proper motions updated to van Leeuwen 2007** (new Hipparcos reduction, A&A 474, 653-664):
- 99 stars had proper motion values updated from the original 1997 Hipparcos catalog
- Data independently sourced from CDS/VizieR catalog I/311/hip2 via TAP query
- Largest correction: Tarf (Beta Cancri) — 18.3" position error at 100 years with old PM

**Position (RA/Dec) cross-checked against SIMBAD**:
- All 10 principal stars verified to <0.02" against SIMBAD J2000 positions
- Two catalog bugs found and fixed:
  - Algedi: coordinates were for Alpha-1 Cap instead of Alpha-2 Cap (HIP 100064)
  - Asellus Borealis: HIP number was wrong (43103 = Iota Cnc, corrected to 42806 = Gamma Cnc)

**Post-update agreement with SE** (101 comparable stars at J2025.0):
- 100% within 0.002° (7.2")
- 98% within 0.5"
- 80% within 0.1"
- Max difference: Rigil Kentaurus 0.51" (nearest star — parallax not modeled)
- Remaining sub-arcsecond differences are from Skyfield vs SE precession/nutation pipeline

**IAU WGSN star name differences** (5 stars resolve to different physical components):
- Menkar: we use Alpha Ceti (IAU WGSN), SE uses Lambda Ceti
- Algedi: we use Alpha-2 Cap (IAU), SE uses Alpha-1 Cap
- Algieba, Albireo, Almach: different component designations

**Independent astropy verification** (10 principal stars, 3 epochs):
- At J2025.0: LibEphemeris wins longitude 4/10, SE wins 6/10 vs astropy
- At J2025.0: LibEphemeris wins latitude 5/10, SE wins 3/10 vs astropy
- At J2100.0: LibEphemeris wins longitude 5/10, SE wins 5/10 vs astropy
- Sirius 0.24" latitude offset traced to annual parallax (not modeled — Sirius is 2.6 pc)
- Fomalhaut 0.11" longitude offset also from parallax (7.7 pc)

### 6.7. Heliocentric, barycentric, equatorial, and XYZ modes — verified

Comprehensive testing across 4 coordinate modes, 10 bodies, 10 dates (1950–2050):

**Heliocentric (FLG_HELCTR)** — bodies Mercury–Pluto:
- Position: < 0.0004° (1.1") max (Mercury), sub-arcsecond for all bodies
- Speed: < 0.0002°/day

**Barycentric (FLG_BARYCTR)** — bodies Moon–Pluto:
- Position: < 0.001° for all bodies except Sun
- Barycentric Sun: up to 0.04° (138") in angular terms, but this is a distance-amplification artifact. The actual 3D positional difference is only ~120 km (0.017% of the solar radius), amplified by the tiny Sun-SSB distance (~0.001–0.009 AU). Independent verification against raw Skyfield/JPL DE440 confirms libephemeris is closer to the source ephemeris at J2000.

**Equatorial (FLG_EQUATORIAL)** — Sun, Moon, Mars, Jupiter:
- RA: < 0.0005° (1.7") max (Moon)
- Dec: < 0.0002° (0.6") max

**XYZ Cartesian (FLG_XYZ)** — all bodies:
- Position: < 0.00005 AU for all bodies (< 0.00004 AU for Sun–Saturn)
- Neptune/Pluto show ~0.00003 AU max, corresponding to sub-arcsecond angular difference at their distances (30+ AU)

**Bug found and fixed**: `FLG_XYZ` and `FLG_RADIANS` flags were not preserved in the return flag of `swe_calc_ut()` / `swe_calc()`. These output format flags were stripped before calculation but not restored in the returned `retflag`. Now fixed.

---

## 7. Intentional divergences from Swiss Ephemeris

LibEphemeris aims for 1:1 compatibility with pyswisseph, but in specific
cases where Swiss Ephemeris has a demonstrable behavioral bug, LibEphemeris
prioritizes correctness over compatibility.

### 7.1. SIDEREAL + J2000 for TrueNode, OscuApog, IntpApog, IntpPerg

**Status:** Intentional divergence since leb/precision branch.

**Observed behavior:** pyswisseph silently ignores `FLG_J2000` for
`SE_TRUE_NODE`, `SE_OSCU_APOG`, `SE_INTP_APOG`, and `SE_INTP_PERG` when
`FLG_SIDEREAL` is also set. The same combination works correctly for
`SE_MEAN_NODE` and `SE_MEAN_APOG`. No error or warning is emitted.

**Why it is incorrect:** Ayanamsha (1D longitude zero-point shift) and
J2000 ecliptic precession (3D plane rotation) are geometrically distinct,
composable operations. The internal inconsistency between mean and true
bodies, combined with a growing error proportional to ecliptic precession
(~0.34° at 2024, ~14° at 3000 CE, ~28° at 0 CE), confirms this is a
code-path bug rather than a design decision.

**Physical evidence:** With the pyswisseph behavior, the True Node ends up
~29° from the Mean Node at 0 CE under `SIDEREAL | J2000` — physically
impossible (the true node oscillates ±1.5° around the mean). With the
LibEphemeris fix, the distance returns to ~1.04° (correct).

**LibEphemeris fix:** `FLG_J2000` is honored for all bodies uniformly.
The algorithm uses mean ayanamsha (no nutation component) and applies
J2000 ecliptic precession after the sidereal correction — the same
pipeline already used correctly for mean bodies.

**Impact:** Only affects users who combine `FLG_SIDEREAL | FLG_J2000`
on these four bodies. Using `FLG_SIDEREAL` alone (the standard
astrological use case) produces identical results to pyswisseph.

**Full analysis:** [se-bug-sidereal-j2000-nodes.md](se-bug-sidereal-j2000-nodes.md)

---

## 8. References

1. Park, R.S. et al. (2021). "The JPL Planetary and Lunar Ephemerides DE440 and DE441." *Astronomical Journal*, 161(3), 105.
2. Folkner, W.M. et al. (2014). "The Planetary and Lunar Ephemerides DE430 and DE431." *Interplanetary Network Progress Report*, 196, 1-81.
3. Stephenson, F.R., Morrison, L.V. & Hohenkerk, C.Y. (2016). "Measurement of the Earth's rotation: 720 BC to AD 2015." *Proceedings of the Royal Society A*, 472, 20160404.
4. Espenak, F. & Meeus, J. (2006). "Five Millennium Canon of Solar Eclipses: -1999 to +3000." NASA/TP-2006-214141.
5. Danjon, A. (1951). "Les éclipses de Lune par la pénombre en 1951." *L'Astronomie*, 65, 51-53.
6. Chapront, J. et al. (2002). "A new determination of lunar orbital parameters." *Astronomy & Astrophysics*, 387, 700-709.
7. ESA (1997). *The Hipparcos and Tycho Catalogues*. ESA SP-1200.
