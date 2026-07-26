# Precision History

Record of precision improvements applied to LibEphemeris (February 2026),
including investigation results and architectural decisions that inform
future development.

> **Historical record.** This file records several superseded investigations.
> The current mean points use published lunar mean-element polynomials, while
> interpolated points use separate Delaunay perturbation series and a
> versioned, hash-pinned compatibility refinement. Current behavior is
> documented in
> [Lunar Nodes and Apsides](../methodology/lunar-apsides.md).

> Goal: make LibEphemeris scientifically precise from independently published
> JPL and IAU models, using **pyerfa** (the Python binding of ERFA).

---

## Table of Contents

1. [Completed Fixes](#completed-fixes)
   - [Formula-anchor ayanamshas (retired)](#formula-anchor-ayanamshas-retired)
   - [Nutation Unification to IAU 2006/2000A](#nutation-unification-to-iau-20062000a)
   - [Obliquity Unification to IAU 2006](#obliquity-unification-to-iau-2006)
   - [Precession Upgrade to IAU 2006 + Frame Bias](#precession-upgrade-to-iau-2006--frame-bias)
   - [Annual Aberration in True Citra](#annual-aberration-in-true-citra)
   - [Central Difference Velocity](#central-difference-velocity)
   - [Cleanup and Code Restoration](#cleanup-and-code-restoration)
2. [Investigations: ELP2000 Perturbations Not Applied](#investigations-elp2000-perturbations-not-applied)
   - [True Node (TRUE_NODE)](#true-node-true_node)
   - [True Lilith (OSCU_APOG)](#true-lilith-oscu_apog)
3. [Open Opportunities](#open-opportunities)
   - [Analytical Chebyshev Velocities](#analytical-chebyshev-velocities)
   - [Lunar geometry improvements](#lunar-geometry-improvements)
4. [Overall Precision Impact](#overall-precision-impact)
5. [Files Modified](#files-modified)
6. [pyerfa Dependency](#pyerfa-dependency)

---

## Completed Fixes

### Formula-anchor ayanamshas (retired)

All predefined IDs 0–46 retain their rc7 calculations without a generic J2000
fallback. Published anchors are cited in the current reference page; exact
legacy encodings still awaiting an independent primary derivation are marked
for owner review and covered by compatibility tests. `SIDM_USER` remains
available for caller-supplied sourced parameters, while IAU/ERFA precession
provides the shared frame evolution.

---

### Nutation Unification to IAU 2006/2000A

**Files:** `libephemeris/planets.py`, `libephemeris/utils.py`

**Problem:**
The main nutation path (`cache.py`) had already been updated to `erfa.nut06a()`,
but 4 secondary code paths still used Skyfield's `iau2000b_radians` (77 terms,
~1 mas), creating a ~1 mas inconsistency with GAST and obliquity in house
calculations. Internally, Skyfield uses `iau2000a_radians` (1365 terms, ~0.1 mas)
for `t.gast` and `ecliptic_frame`, so mixing IAU 2000B in other paths introduced
systematic discrepancies.

**Fix applied:**

| Location | Before | After |
|----------|--------|-------|
| `planets.py:_get_true_ayanamsa()` | `iau2000b_radians(t_obj)` | `erfa.nut06a(2451545.0, t_obj.tt - 2451545.0)` |
| `planets.py:_calc_ayanamsa_ex()` | `iau2000b_radians(t_obj)` | `erfa.nut06a(2451545.0, tjd_tt - 2451545.0)` |
| `utils.py:azalt()` | `iau2000b_radians(t)` | `erfa.nut06a(2451545.0, jd_tt - 2451545.0)` |
| `utils.py:azalt_rev()` | `iau2000b_radians(t)` | `erfa.nut06a(2451545.0, jd_tt - 2451545.0)` |

The import `from skyfield.nutationlib import iau2000b_radians` was removed from
`planets.py` (no longer used). `import erfa` was added to `utils.py`.

**Impact:**

- **Before:** ~1 mas inconsistency between code paths (IAU 2000B, 77 terms)
- **After:** ~0.01–0.05 mas uniform precision (IAU 2006/2000A via pyerfa)

---

### Obliquity Unification to IAU 2006

**Files:** `libephemeris/planets.py`, `libephemeris/utils.py`

**Problem:**
Three different formulas for mean obliquity were in use:

| Location | Formula | Constant |
|----------|---------|----------|
| `cache.py` (houses, ayanamsha) | IAU 2006 via `erfa.obl06()` | 84381.406″ (already fixed) |
| `_calc_ayanamsa_ex()` | Laskar 1986 (`23.43929111` = `84381.448″`) | 84381.448″ |
| `utils.py:azalt/azalt_rev` | Mixed (IAU 2006 constant but only 3 terms) | 84381.406″ (incomplete) |

The critical path (houses, ayanamsha) used Laskar 1986, with a constant offset of
**0.042″** from IAU 2006 and missing terms (T⁴, T⁵).

**Fix applied:**

**1. `_calc_ayanamsa_ex()`** — replaced Laskar 1986 formula with `erfa.obl06()`:

```python
# Before:
eps0 = 23.43929111 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T / 3600.0
dpsi_rad, deps_rad = iau2000b_radians(t_obj)

# After:
eps0 = math.degrees(erfa.obl06(2451545.0, tjd_tt - 2451545.0))
dpsi_rad, deps_rad = erfa.nut06a(2451545.0, tjd_tt - 2451545.0)
```

**2. `utils.py:azalt()` and `azalt_rev()`** — replaced manual formula with
`erfa.obl06()`:

```python
# Before (3 terms, missing T⁴ and T⁵):
eps0 = (84381.406 - 46.836769 * T - 0.0001831 * T*T + 0.00200340 * T*T*T) / 3600.0
dpsi_rad, deps_rad = iau2000b_radians(t)

# After:
eps0_rad = erfa.obl06(2451545.0, jd_tt - 2451545.0)
eps0 = math.degrees(eps0_rad)
dpsi_rad, deps_rad = erfa.nut06a(2451545.0, jd_tt - 2451545.0)
```

**Impact:**

- **Before:** 0.042″ constant offset (Laskar vs IAU 2006) + missing T⁴/T⁵ terms
- **After:** Consistent IAU 2006 obliquity everywhere via pyerfa

---

### Precession Upgrade to IAU 2006 + Frame Bias

**File:** `libephemeris/planets.py`, function `_get_star_position_ecliptic()`

**Problem:**
The precession coefficients zeta/z/theta used values from **Lieske 1977 (IAU 1976)**,
not IAU 2006. The code comment erroneously declared "IAU 2006 precession formulas".
Additionally, the GCRS→J2000 frame bias (~23 mas) was missing.

```python
# Original code (Lieske 1977, NOT IAU 2006):
zeta  = (2306.2181 * T + 0.30188 * T**2 + 0.017998 * T**3) / 3600.0
z     = (2306.2181 * T + 1.09468 * T**2 + 0.018203 * T**3) / 3600.0
theta = (2004.3109 * T - 0.42665 * T**2 - 0.041833 * T**3) / 3600.0
```

The correct IAU 2006 coefficients (Capitaine et al. 2003, A&A 412, Table 1)
include constant terms (±2.650545″) that encode the frame bias:

```python
zeta  = (2.650545 + 2306.083227 * T + 1.0967790 * T**2
         + 0.01860606 * T**3 - 0.000013 * T**4
         - 0.0000005 * T**5) / 3600.0
z     = (-2.650545 + 2306.077181 * T + 1.0927348 * T**2
         + 0.01826837 * T**3 - 0.000028 * T**4
         - 0.0000003 * T**5) / 3600.0
theta = (2004.191903 * T - 0.4294934 * T**2
         - 0.04182264 * T**3 - 0.000007089 * T**4
         - 0.0000001274 * T**5) / 3600.0
```

Dead code was also present (scalar A/B/C calculation that was never used,
overwritten by the subsequent rotation matrix).

**Intermediate fix:**
The entire manual precession block (Lieske 1977 + scalar rotations + dead code)
was replaced with `erfa.pmat06()`:

```python
# Precession-bias matrix IAU 2006 (includes frame bias)
rbp = erfa.pmat06(2451545.0, tjd_tt - 2451545.0)

# Application: P_date = rbp @ P_J2000
x3 = rbp[0][0] * x0 + rbp[0][1] * y0 + rbp[0][2] * z0
y3 = rbp[1][0] * x0 + rbp[1][1] * y0 + rbp[1][2] * z0
z3 = rbp[2][0] * x0 + rbp[2][1] * y0 + rbp[2][2] * z0
```

`erfa.pmat06()` automatically includes:
- Frame bias GCRS→J2000 (~23 mas, the constant terms ±2.650545″)
- IAU 2006 precession with all terms up to T⁵
- Exact coefficients from Capitaine et al. 2003

**Final fix:**
This intermediate fix was subsequently superseded by the complete rewrite of
`_get_star_position_ecliptic()` as part of the annual aberration fix (see below).
The `erfa.pmat06()` code was maintained as an intermediate solution but then
replaced by the Skyfield pipeline, which handles precession internally.

**Impact:**

- **Before:** Lieske 1977 (~0.1″/century) + missing frame bias (~23 mas)
- **After:** Exact IAU 2006 with frame bias included

---

### Annual Aberration in True Citra

**File:** `libephemeris/planets.py`, function `_get_star_position_ecliptic()`

**Problem:**
The function computed Spica's ecliptic position for True Citra without applying
annual aberration (Bradley's independently published effect).
For comparison, `fixed_stars.py` correctly applied aberration using Skyfield's
`astrometric.apparent()`.

The original implementation consisted of ~130 lines of manual code:
1. Proper motion with 3D vector (Hipparcos Vol. 1, Sec. 1.5.5)
2. Precession with scalar formulas (Lieske 1977) + rotation matrix
3. Manual ecliptic conversion with obliquity

All without aberration, and with Lieske 1977 instead of IAU 2006.

**Fix applied (Option B — recommended in TODO):**
The entire function (~130 lines) was rewritten to use the Skyfield pipeline:

```python
def _get_star_position_ecliptic(star, tjd_tt, eps_true):
    star_obj = Star(
        ra_hours=star.ra_j2000 / 15.0,
        dec_degrees=star.dec_j2000,
        ra_mas_per_year=star.pm_ra * 1000.0,
        dec_mas_per_year=star.pm_dec * 1000.0,
        parallax_mas=star.parallax * 1000.0 if star.parallax > 0 else 0.0,
        radial_km_per_s=star.radial_velocity,
    )

    planets = get_planets()
    ts = get_timescale()
    t = ts.tt_jd(tjd_tt)
    earth = planets["earth"]

    pos = earth.at(t).observe(star_obj).apparent()
    lat, lon, dist = pos.frame_latlon(ecliptic_frame)

    return lon.degrees
```

The Skyfield pipeline `observe().apparent().frame_latlon()` automatically includes
7 effects:
- **Proper motion** (rigorous propagation with radial velocity)
- **Light-time correction** (light propagation time)
- **Annual aberration** (~20.5″ — the primary fix)
- **Gravitational deflection** (gravitational light bending)
- **Precession IAU 2006** (with GCRS→J2000 frame bias)
- **Nutation IAU 2000A** (1365 terms, ~0.1 mas)
- **Ecliptic transformation** (true ecliptic of date)

**Dead code removed:**

- ~80 lines of manual proper motion propagation
- ~50 lines of manual precession (Lieske 1977 + rotations)
- Scalar dead code A/B/C (including a `# Wait, this is incomplete` comment)
- Manual ecliptic conversion with obliquity

This fix simultaneously resolved the aberration issue, the precession upgrade,
the frame bias, and the dead code cleanup.

**Impact:**

- **Before:** up to 20.5″ error (missing aberration) + ~0.1″/cy (Lieske) + ~23 mas (frame bias)
- **After:** sub-milliarcsecond precision (all handled by Skyfield)

---

### Central Difference Velocity

**Files:** `planets.py`, `hypothetical.py`, `spk.py`, `planetary_moons.py`

**Problem:**
All non-planetary velocities used **forward difference** O(h):
```
f'(x) ≈ (f(x+h) - f(x)) / h
```
instead of **central difference** O(h²):
```
f'(x) ≈ (f(x+h) - f(x-h)) / (2h)
```

Central difference provides ~100× better precision for the same timestep.

Additionally, the velocity correction for the ayanamsha rate used forward difference
even in code paths where the primary velocity already used central difference.

**Fix applied:**

#### 7a — Central difference for non-planetary bodies

| File | Body | Method before | Method after |
|------|------|---------------|--------------|
| `hypothetical.py` | Harrington | Forward difference | Central difference |
| `spk.py` | SPK Type 2/3 fallback | Forward 1s | Central 1s |
| `planetary_moons.py` | Planetary moons | Forward 1s | Central 1s |

An earlier version of the table also listed hypothetical bodies whose complete
numerical model is not available from a primary publication; those IDs are no
longer computed (see
[missing hypothetical models](../methodology/missing-hypothetical-models.md)).

Example transformation (schematic):

```python
# Before (forward difference):
pos_next = _calc_supported_body_raw(jd_tt + dt_step)
dlon = pos_next[0] - longitude

# After (central difference):
pos_prev = _calc_supported_body_raw(jd_tt - dt_step)
pos_next = _calc_supported_body_raw(jd_tt + dt_step)
dlon = (pos_next[0] - pos_prev[0]) / (2.0 * dt_step)
```

#### 7b — Ayanamsha rate correction with central difference

**7 occurrences** in `planets.py` (lines 1046, 1142, 1176, 1226, 1258, 1309, 1772)
where the velocity correction for the ayanamsha rate used forward difference:

```python
# Before (forward difference):
ayanamsa_next = _get_true_ayanamsa(t.ut1 + dt)
da = (ayanamsa_next - ayanamsa) / dt

# After (central difference):
ayanamsa_prev = _get_true_ayanamsa(t.ut1 - dt)
ayanamsa_next = _get_true_ayanamsa(t.ut1 + dt)
da = (ayanamsa_next - ayanamsa_prev) / (2.0 * dt)
```

**Impact:**

- **Secondary body velocities:** ~100× precision improvement for the same timestep
- **Ayanamsha rate:** consistent O(h²) everywhere, eliminated inconsistency with primary velocity

---

### Cleanup and Code Restoration

#### Dead code in `_get_star_position_ecliptic()` — Resolved

The scalar A/B/C calculation (which also contained a `# Wait, this is incomplete`
comment) and all manual propagation code were removed as part of the Skyfield
rewrite (annual aberration fix above).

#### Dead function `_calc_star_based_ayanamsha()` — Not touched

The function still exists but is never called. It was not removed to minimize
risk, since it could serve as a future reference.

#### Stale comment in `planets.py:34` — Already fixed

Had already been corrected prior to this work.

#### Obliquity in `_calc_ayanamsa()` — Already fixed

Had already been corrected prior to this work (uses `erfa.obl06()`).

#### Accidentally removed code — Restored

During editing of `planets.py`, three functions and one class were accidentally
removed (located between `_calc_body_special` and `_calc_body`):

- `NutationFallbackWarning` (class) — Warning for degraded precision
- `get_nutation_model()` — Check of the active nutation model
- `_calc_nutation_obliquity()` — Nutation/obliquity calculation for `ECL_NUT`
- `_maybe_equatorial_convert()` — Ecliptic→equatorial conversion

All four were **restored**, and `_calc_nutation_obliquity()` was updated to use
`erfa.obl06()` and `erfa.nut06a()` instead of Skyfield.

---

### Photometric Models (pheno_ut) — Phase Angle, Magnitude, Diameter

**File:** `libephemeris/planets.py`, functions `_calc_pheno()`, `_calc_planet_magnitude()`, `_calc_moon_magnitude()`

**Date:** March 2026

The `pheno_ut()` implementation was audited against Mallama & Hilton (2018),
published almanac formulae, vector identities, and IAU 2015 standards.

#### Phase Angle — 3D Vector Dot Product (Critical)

**Problem:**
The phase angle (Sun-Body-Earth angle) was computed using the law of cosines on
the triangle formed by heliocentric distance, geocentric distance, and Sun-Earth
distance. For the Moon this triangle is extremely elongated, making the
law-of-cosines form numerically ill-conditioned.

**Fix applied:**
Replaced with 3D vector dot product approach using geocentric position vectors:

```python
# Vector from body to Sun: sun_geo - body_geo
# Vector from body to Earth: -body_geo (Earth is at origin)
# Phase angle = angle between these vectors via dot product
cos_phase = dot(body_to_sun, body_to_earth) / (|body_to_sun| * |body_to_earth|)
phase_angle = arccos(clamp(cos_phase, -1, 1))
```

The same stable vector identity is applied to all planets and validated against
the independently computed JPL geometry.

#### Moon Magnitude — Astronomical Almanac Formula

**Problem:**
The previous empirical lunar formula did not follow the project's selected
published almanac convention.

**Fix applied:**
Replaced with the Astronomical Almanac formula (Allen's Astrophysical Quantities):

```python
V = -12.73 + 0.026 * |alpha| + 4e-9 * |alpha|^4
```

with distance correction `5 * log10(d / d_mean)`.

The published approximation has its own stated domain and becomes less reliable
for very thin crescents; that limitation is documented rather than corrected
from external comparison output.

#### Sun Magnitude — V(1,0) Correction

**Problem:** V(1,0) was -26.74, should be **-26.86** (Mallama & Hilton 2018).
Constant 0.12 mag offset.

**Fix:** Changed the constant to the published value.

#### Neptune Magnitude — Secular Brightness Variation

**Problem:** A time-invariant Neptune magnitude model does not represent the
published secular brightness changes described by Lockwood & Thompson (1991)
and Sromovsky et al. (2003).

Additionally, the `tjd` parameter was not being passed to `_calc_planet_magnitude()`
for non-Saturn planets (it was initialized to `0.0` instead of `t.tt`), so the
year calculation gave ~4712 BC regardless of actual date.

**Fix applied (2 changes):**

1. Implemented the independently published secular brightness behavior.
2. Changed `tjd = 0.0` to `tjd = t.tt` for all planets so the model receives the
   actual calculation epoch.

**References:** Lockwood & Thompson (1991), Sromovsky et al. (2003).

#### Uranus Magnitude — V(1,0) Correction

**Problem:** V(1,0) was -7.19, should be **-7.15** (Mallama & Hilton 2018).

**Fix:** Changed the constant and phase coefficient to the published
Mallama–Hilton values.

#### Apparent Diameter — IAU 2015 Equatorial Radii (Intentional Divergence)

**Investigation:** LibEphemeris uses **IAU 2015 equatorial radii** as the
convention for maximum apparent angular diameter. Mean volumetric radii answer
a different geometric question.

**Decision: NOT CHANGED.** The equatorial-radius convention is explicit and
sourced from IAU 2015 Resolution B3.

---

## Investigations: ELP2000 Perturbations Not Applied

### True Node (TRUE_NODE)

**File:** `libephemeris/lunar.py`, `calc_true_lunar_node()`

**Investigation:**
An experimental analytical perturbation series was conceived as a correction
to the geometric node. Investigation revealed that:

1. The ELP2000 perturbation series is designed to correct the **mean node**,
   not the geometric node calculated with `h = r × v`
2. The geometric node already captures perturbations through the JPL DE state
   vectors (which include all planetary perturbations)
3. Direct application of the series to the geometric node produced erroneous
   results (deviations of tens of degrees)

**Decision: NOT APPLIED.** The unused experimental series and its table-only
tests were removed during the 2026-07-13 provenance cleanup. The live geometric
implementation was unchanged.

### True Lilith (OSCU_APOG)

**File:** `libephemeris/lunar.py`, `calc_true_lilith()`

**Investigation:**
Analogous situation to True Node. The function `_calc_elp2000_apogee_perturbations()`
(~50 terms) is designed for the **interpolated apogee** (mean apogee, `INTP_APOG`),
not for the osculating apogee calculated with the eccentricity vector
`e = (v×h)/μ − r/|r|`.

The eccentricity-vector construction is an instantaneous two-body reduction of
the full JPL Earth–Moon state. It remains a conventional osculating element and
should not be conflated with the smoothed apogee definition.

**Decision: NOT APPLIED.** An explanatory comment was added:

```python
# Note: ELP2000/Moshier perturbation corrections (_calc_elp2000_apogee_perturbations)
# are available but not applied here. The perturbation series was designed for the
# interpolated (mean) apogee, not the osculating eccentricity vector.
```

## Open Opportunities

### Analytical Chebyshev Velocities

**Impact: Moon ~0.001 deg/day; planets smaller but systematic. Performance ~3×.**

**File:** `libephemeris/planets.py`

All planetary velocities are currently calculated with numerical finite-difference
(central difference O(h²)). Skyfield/jplephem provides analytical ICRS velocities
from differentiation of the DE440 Chebyshev polynomials, which are already used
for COB corrections but **not** for angular ecliptic velocities.

**Proposed approach** — obtain position and velocity from Skyfield and compute
angular velocity with exact derivatives:

```python
r_ecl, v_ecl = pos.frame_xyz_and_velocity(ecliptic_frame)
x, y, z = r_ecl.au
vx, vy, vz = v_ecl.au_per_d

xy_sq = x * x + y * y
r_sq = xy_sq + z * z
xy = math.sqrt(xy_sq)
r = math.sqrt(r_sq)

speed_lon = math.degrees((x * vy - y * vx) / xy_sq)         # dλ/dt
speed_lat = math.degrees((z * (x*vx + y*vy) / xy - xy * vz) / r_sq)  # dβ/dt
speed_dist = (x * vx + y * vy + z * vz) / r                  # dr/dt
```

This approach is already implemented for:
- SPK Type 21 (`spk.py:1006–1023`)
- Fixed stars (`fixed_stars.py:3491–3530`)

Extending it to all standard planets would eliminate the 2 recursive `_calc_body()`
calls currently used for finite-difference, yielding ~3× performance improvement.

---

### Lunar geometry improvements

Future changes to True Node or lunar apsides start from JPL states, IAU/ERFA
arguments, or citable primary literature, rather than from output-offset
calibration or externally fitted perturbation terms. A scientifically motivated
improvement can refine a published analytical definition or use a documented
multi-body dynamical construction validated against JPL identities while
preserving the public compatibility tolerances.

---

## Overall Precision Impact

The completed work unified nutation and obliquity on ERFA, adopted IAU/Vondrák
frame models, included the standard frame-bias transform, and replaced forward
finite differences with wrap-aware central differences where analytical state
derivatives are unavailable. These changes are validated against independent
standards and mathematical identities rather than retained comparison deltas.

The principal remaining opportunity is to propagate analytical JPL/LEB state
derivatives farther through the coordinate pipeline for performance and
numerical consistency.

### IAU Models Used

| Component | Model | Source | Precision |
|-----------|-------|--------|-----------|
| Nutation | IAU 2006/2000A | `erfa.nut06a()` | ~0.01–0.05 mas |
| Mean obliquity | IAU 2006 | `erfa.obl06()` | Sub-milliarcsecond |
| Stellar precession | IAU 2006 | Skyfield `ecliptic_frame` | Sub-milliarcsecond |
| Ayanamsha precession | IAU 2006 | Capitaine et al. 2003 (5 terms) | ~0.08 mas/cy |
| Aberration | Complete | Skyfield `.apparent()` | ~0.001″ |
| Ephemerides | JPL DE440/DE421 | Skyfield | ~1 mas |

---

## Files Modified

| File | Lines changed (approx.) | Issues resolved |
|------|------------------------|-----------------|
| `libephemeris/planets.py` | ~400 | Aberration, nutation, obliquity, PREC_RATE regression, precession, frame bias, central diff, cleanup |
| `libephemeris/utils.py` | ~40 | Nutation, obliquity |
| `libephemeris/lunar.py` | ~10 | True Node / True Lilith investigation notes |
| `libephemeris/hypothetical.py` | ~60 | Central difference velocity |
| `libephemeris/spk.py` | ~40 | Central difference velocity |
| `libephemeris/planetary_moons.py` | ~40 | Central difference velocity |

**Files not modified:** `astrometry.py` (already fixed), `fixed_stars.py` (already correct).

---

## pyerfa Dependency

**pyerfa** (`>=2.0.0`) is used as a required dependency.

| pyerfa function | Usage | Precision |
|-----------------|-------|-----------|
| `erfa.nut06a()` | Nutation in 4 code paths | ~0.01–0.05 mas |
| `erfa.obl06()` | Obliquity in 3 code paths | Sub-mas |
| `erfa.pmat06()` | Stellar precession (now via Skyfield) | Sub-mas |

Footprint is ~2 MB, pure C (IAU ERFA/SOFA wrapper), no significant transitive
dependencies.
