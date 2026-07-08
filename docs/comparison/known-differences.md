# Known differences vs Swiss Ephemeris

These are not bugs. They are the inherent consequences of using different
ephemerides, different astronomical models, or different algorithmic strategies.
Each is explained so that users and developers know exactly what to expect. Part A
gives the root-cause narrative for the headline differences; Part B is a granular,
per-API catalog verified through systematic hyper-validation (4400+ comparison
rounds across 29 API sections).

---

## Part A — Root causes

### 2.1 Crossing functions — full-orbit search for slow planets

`swe_cross_ut` handles full-orbit searches for all planets, including slow outer
planets like Jupiter and Saturn. The algorithm scales the Brent bracket search
window by the estimated time to crossing (`dt_guess`) and filters out false sign
changes at the antipodal point (target ± 180°) during the coarse scan, so it
converges even when the crossing is 10+ years away and the planet undergoes
multiple retrograde periods en route. (Earlier builds could converge on 180°
instead of 0° for Jupiter; fixed by scaling the window to
`max(min_window, dt_guess * 1.5)`, scaling coarse samples proportionally, and
rejecting bracket candidates whose signed-difference jump exceeds 180°.)

### 2.2 Moon precision

**What:** the maximum difference in lunar longitude is ~135" (0.037°), occurring at
the extremes of the DE440 range (around 1550 CE and 2650 CE).

**Why:** the two libraries use fundamentally different lunar models. libephemeris
uses **JPL DE440**, where the Moon's position comes from numerical integration of
the full equations of motion, fitted to Lunar Laser Ranging data (1970–present) at
~1 mas precision for the modern era. pyswisseph combines **DE431** for the modern
period with the **ELP/MPP02 analytical theory** (a sum of thousands of
trigonometric terms). Near the present (1900–2100) the two agree to a few
arcseconds; they diverge toward historical/far-future dates because of
(1) different tidal-acceleration parameters, (2) different fitting datasets (DE440
adds 8 years of LLR plus Juno/MESSENGER data), and (3) numerical-integration vs
analytical-truncation error profiles.

**Practical impact:** for dates within ±200 years the difference is < 5". For natal
charts, transits and progressions it is astronomically negligible (~4 minutes of
lunar motion).

### 2.3 Delta T

**What:** the computed ΔT (TT − UT1) can differ by up to ~232 s between the two
libraries for extreme dates.

**Why:** both libraries use the **same family** of ΔT model — Stephenson, Morrison
& Hohenkerk (2016) for the pre-1955 era and IERS observed values for the modern
era. (Swiss Ephemeris used Espenak & Meeus 2006 up to SE 1.77, then switched to
SMH-2016 at SE 2.06/2.08 in 2017–2018. Verified empirically: `swe.deltat` at 2100 =
93.2 s, matching libephemeris's 95.9 s and nowhere near the Espenak-Meeus value of
202.8 s.) The difference is therefore **not** a different model — it is the
unavoidable difference between two *extrapolations of an unknowable quantity*
(future Earth rotation, or deep-past rotation beyond the eclipse record). Where ΔT
is constrained by data the two agree to ≤ ~1 s.

| Era | libephemeris | pyswisseph 2.10 | difference |
|---|---|---|---|
| 1962 → present | IERS observed (ΔT(2000) = 63.83 s) | IERS observed (ΔT(2000) = 63.83 s) | identical |
| 1620 → 1955 (telescopic) | Skyfield ΔT table / SMH-2016 | tabulated values | within ~7 s (observational scatter) |
| −720 → 1620 | SMH-2016 | SMH-2016 | close (same family) |
| deep past / future | SMH-2016 long-term tidal tail | a different long-term extrapolation | the only meaningful divergence |

**The difference is purely ΔT, not the ephemeris.** Forcing the *same* ΔT into both
engines collapses the position difference to the ephemeris-model floor:

| date | Moon Δ (native ΔT) | Moon Δ (same ΔT) | worst planet (same ΔT) |
|---|---|---|---|
| 2100 | 1.6″ | 0.000″ | 0.04″ |
| 2200 | 29.3″ | 0.000″ | 0.08″ |
| 2300 | 87.0″ | 0.000″ | 0.06″ |

At a common ΔT, the engines agree to < 0.1″ (Moon to ~0.001″) even at year 2300 —
the far-epoch divergence is entirely the ΔT *extrapolation choice*. For any real
chart (natal/transit/progression, and even ancient astrology back to ~1000 BCE) the
ΔT difference is far below the arcminute working precision of astrology; it only
matters for far-future or deep-past research dates, where ΔT is genuinely uncertain
for everyone. libephemeris's ΔT model and selector are documented in
[../methodology/delta-t.md](../methodology/delta-t.md).

**Validation technique — forcing a matched ΔT.** The tables above were produced by
injecting the other engine's ΔT into libephemeris per date, so the comparison
isolates the pure ephemeris-model difference from the ΔT choice:

```python
import swisseph as swe
import libephemeris as L
L.set_delta_t_userdef(swe.deltat(jd))   # force pyswisseph's ΔT into libephemeris, per date
# ... compare L.calc_ut(jd, ...) vs swe.calc_ut(jd, ...) ...
L.set_delta_t_userdef(None)             # restore the native model
```

> ⚠️ Do **not** reset pyswisseph's own ΔT mid-loop with
> `swe.set_delta_t_userdef(-1.0)` / `AS_UNDEF`: it leaves pyswisseph in a state that
> injects a spurious constant ~3553″ (~0.987°) offset into subsequent
> `swe.sidtime` / `swe.deltat` calls. The clean method forces ΔT only into
> libephemeris and lets pyswisseph use its native `swe.deltat(jd)` (which equals
> the forced value). This injection works on libephemeris's LEB/fast/Horizons
> backends, which honour `deltat()`; the `"skyfield"` backend uses its own internal
> ΔT (see [../methodology/delta-t.md](../methodology/delta-t.md)).

### 2.4 House cusp speeds — numerical vs analytical derivatives

`houses_ex2` / `houses_armc_ex2` compute cusp velocities by centered finite
differences when `FLG_SPEED` is set; the maximum difference from pyswisseph is
~0.7°/day (~0.2% relative). pyswisseph differentiates the cusp formulas
analytically. Both approaches are valid; the numerical method has a truncation
error ∝ dt² that at a 1-minute step is < 1°/day against the analytical result, i.e.
< 0.3% of the ~280–340°/day cusp speeds. (Cusp speeds are returned only when
`FLG_SPEED` is passed.)

### 2.5 Asteroids — Keplerian approximation vs integrated ephemerides

Minor-body positions (Chiron, Ceres, Pallas, Juno, Vesta) can differ by up to ~5°.
pyswisseph uses pre-integrated NASA JPL asteroid ephemerides (sub-arcsecond).
libephemeris uses a **Keplerian approximation** by default — an ideal ellipse from
mean orbital elements — which ignores planetary perturbations, resonances and
secular element variation, so it degrades with time (Chiron is the most affected,
its orbit between Saturn and Uranus being chaotic).

**Mitigation:** libephemeris loads **SPK kernels** from JPL Horizons for any minor
body, giving sub-arcsecond accuracy matching Horizons:

```python
import libephemeris as eph
eph.set_auto_spk_download(True)  # automatic download from JPL Horizons
```

The Keplerian fallback exists only when SPK data is unavailable. (The DE440/DE441
ephemerides cover only the major planets and the Moon; asteroids require separate
data sources in any engine.)

### 2.6 House cusps and Ascendant at remote epochs

**What:** inside the modern range the house cusps (ASC/MC/ARMC) match pyswisseph to
~0.002″ (see [precision §1](precision.md), < 0.02″). Outside ≈ 1850–2050, *even
with the same ΔT forced on both sides*, the cusps diverge by a secular,
sign-changing amount: ASC residual ≈ +1.5″ (2050), +0.9″ (2100), −1.3″ (2200),
−5.6″ (2300), growing to ~0.1–1° at ±3000–5000 yr.

**Why:** over the range where real charts live the residual is essentially all in
the **apparent sidereal time (ARMC)**. In the modern era true/mean obliquity and
nutation (Δψ, Δε) agree with pyswisseph to < 0.002″, so the ASC residual in the
table below is driven by ARMC alone. (At deep-BCE epochs the of-date *mean
obliquity* intentionally differs — libephemeris uses the long-term Vondrák
pole-angle obliquity — but that deviation affects only ecliptic latitude, not
longitude or ARMC; see
[Intentional divergences §3](intentional-divergences.md#3-of-date-mean-obliquity-at-deep-bce-epochs).)
The ARMC difference is a *model* difference: libephemeris uses the IAU-2006 GMST
polynomial inside 1850–2050
(same branch as pyswisseph → ~0.002″) and a geometric Vondrák-2011 long-term
sidereal time outside it, whereas pyswisseph continues an IAU-2006-style
precession-in-RA realization that diverges at remote epochs. The long-term
construction is detailed in
[../methodology/sidereal-time-longterm.md](../methodology/sidereal-time-longterm.md).

Measured against pyswisseph **with the same ΔT forced on both sides**, the
Ascendant residual is:

| epoch | ASC residual | dominant component |
|---|---|---|
| 1850–2050 | ~0.002″ | identical IAU-2006 GMST branch |
| 2100 | +0.9″ | long-term sidereal-time model gap |
| 2200 | −1.3″ | (sign change near 2150) |
| 2300 | −5.6″ | grows secularly with \|epoch − 2050\| |

The MC residual tracks ARMC 1:1; the ASC is that residual times the latitude
Jacobian ∂Asc/∂ARMC, which *damps* the error at high northern latitude (≈ 0.56 at
65°N) and slightly amplifies it in the south (≈ 1.38 at 34°S) — there is no polar
blow-up.

**The ~1.9″ step seen near 2050 originates in Swiss Ephemeris's own model
boundary.** A sub-day kink test at JD 2469807.5 (2050-01-01 00:00) shows
libephemeris's ARMC incrementing smoothly across the boundary (constant
~1299.548″ per 0.001 d), while pyswisseph has a single anomalous increment short by
~1.908″ at that instant — an internal Swiss precession/sidereal-time model
boundary. libephemeris's own two branches join to 0.000000″ there (the continuity
offset pins the long-term branch to the IAU-2006 value at the boundary).

**Verdict:** a benign, expected model difference. Inside 1850–2050, where
essentially all real charts live, agreement is ~0.002″; the remote-epoch divergence
is the intentional Vondrák-2011 long-term design, which stays arcsecond-accurate
over ±200 millennia where a truncated IAU-2006 polynomial drifts by degrees. No
code change is warranted to "match" pyswisseph there.

---

## Part B — Per-API divergence catalog

All divergences below were verified through systematic hyper-validation against
pyswisseph 2.10.03 (4400+ comparison rounds across 29 API sections).

### Summary

| Category | Typical | Maximum | Cause |
|----------|---------|---------|-------|
| Planetary positions | 0.001–0.5" | ~1" (Moon) | Different ephemeris engines |
| Planetary speeds | 0.01–5" | ~20" (Moon speed) | Numerical differentiation methods |
| House cusps | <0.01" | ~0.01" | Obliquity/nutation model differences |
| Fixed star positions | <0.01" | <0.01" | Proper-motion catalog differences |
| Fixed star distances | <0.01% at J2000 | ~0.1% at ±50y | Radial velocity models |
| Ayanamsha values | <0.1" | ~40" (exotic modes) | Reference-star position differences |
| Delta-T | <0.001s (modern) | ~1s (1900-1950), ~232s (extremes) | Different ΔT extrapolation |
| Refraction | <1" | ~15" | Different atmospheric models |
| Phase angles | <1" (inner) | ~18" (outer planets) | Position errors amplified |
| Orbital elements | <1" (inner) | ~1° (giant `varpi`) | Convention; lib matches exact two-body to 0.000000° (see §9) |
| Sidereal Moon | ~3" | ~14" | Lunar theory differences in sidereal frame |

### 1. Ephemeris engine differences

**1.1 Planetary positions (`calc_ut`).** libephemeris uses Skyfield + JPL
DE440/DE441; pyswisseph uses its internal integration of DE431. Typical divergence:
Sun and Mercury–Neptune 0.001–0.01"; Moon 0.01–1.0"; Pluto 0.001–0.01". Speeds:
most planets 0.01–2.0", Moon speed up to ~5" (central finite difference vs
analytical). Future dates (>2050): up to 2" from ΔT extrapolation.

*Extreme / BCE dates.* Validated at the *same* TT (so ΔT cancels) against an exact
JPL **DE441** oracle (jplephem) from ~3000 BCE to +5000 CE: libephemeris (extended
LEB tier, fit to DE441) tracks the exact ephemeris to **< 0.01"** across the whole
range for Sun, Moon and inner planets (< 0.09" for the giants, the residual being
the center-vs-barycenter choice). pyswisseph's geocentric apparent positions vs the
same DE441 truth grow toward the past — ~11" (Moon, 1000 CE) → ~157" (1000 BCE) →
~550" (~0.15°, Moon, 3000 BCE); the planets stay smaller (a few arcsec to ~1′) — because it is built on DE431 (2013) plus its own
lunar secular model. For historical/archaeo-astronomy work this makes libephemeris
the closer match to JPL's latest long-range integration.

**1.2 Interpolated apogee/perigee (bodies 21, 22).** IntpApog/IntpPerg use
semi-analytical ELP2000-82B perturbation theory, implemented independently from
published coefficients; results can differ by several arcseconds. Classified as a
known divergence in all tests.

**1.3 Pholus (body 16) historical dates.** The SPK auto-download path requests
padding around the selected tier and verifies cached-kernel coverage before reuse,
so Pholus at 1900-01-01 is now served from a Horizons SPK kernel instead of a
fallback orbit. Against Horizons geocentric observed ecliptic coordinates, the
library is within ~0.25" in longitude and ~0.05" in latitude at JD 2415021.0
(`validation/verify/verify_pholus_horizons.py`, online verifier). The reference
API remains ~28" away at that date, apparently from a different historical orbit
solution; that residual is certified rather than copied.

### 2. House system differences

**2.1 House cusps.** Typically agree within 0.01"; the small divergence comes from
slightly different obliquity and nutation models in the intermediate calculations.

**2.2 Vertex at the equator.** At latitude 0° the Vertex has a 1/tan(lat)
singularity; libephemeris clamps latitude to a tiny epsilon to evaluate the
limiting value, matching pyswisseph. The only remaining divergence is `ascmc[6]`
(CoAsc Munkasey) for the Horizon (H) system at lat 0, where pyswisseph returns 0°
from a C `tan(90°)` artifact while the mathematical limit is 180°.

**2.3 House position (`house_pos`).** < 0.01" for most systems; up to ~46" for
**Alcabitius (B)**, **Koch (K)**, **Topocentric (T)** from different cusp
interpolation algorithms.

**2.4 MC "flip" for quadrant house systems inside the polar circle (parity
behavior, not a divergence).** At latitudes inside the polar circle where the
culminating point of the ecliptic is below the horizon, the quadrant house
systems **Regiomontanus (R), Campanus (C), Polich-Page/Topocentric (T), APC
(Y), Sunshine (I)** and **Savard-A (J)** return an MC that is 180° away from
the classical astronomical Midheaven, together with cusps whose consecutive
gaps run backwards (they sum to 3960° instead of 360°). The
ecliptic-proportional systems **Porphyry (O), Equal-from-MC (A), Meridian (X)**
do NOT flip and return the raw upper-meridian MC; **Placidus (P)** and **Koch
(K)** are undefined inside the polar circle and raise an error. All of this
exactly matches the reference API (0 mismatches over an 800-case
ARMC × latitude × system grid; every value in GitHub issue #46's reproduction
is reproduced by the reference API to the 6th decimal, including the onset
between −66.5° and −67°, the T-system gap sum of 1080° at the −66° boundary,
and the Sunshine `I` cusp collapse under a circumpolar Sun in *both*
hemispheres while its sibling `i` raises).

*Why:* the astronomical MC is the intersection of the ecliptic with the
*upper* meridian, i.e. the ecliptic point whose right ascension equals the
ARMC (`MC = atan2(tan(ARMC), cos(ε))`, upper-meridian branch). Inside the
polar circle that point can be circumpolar-below: its upper-transit altitude
`90° − |lat − δ_MC|` is negative (with `δ_MC = atan(sin(ARMC)·tan(ε))`). For
quadrant systems the reference API adopts the convention that the MC must be
the meridian∩ecliptic intersection *above the horizon*; when the upper-meridian
point is below it, the convention switches to the opposite branch
(`RA = ARMC + 180°`), which is geometrically the **IC** (lower-meridian,
anti-culmination point) but is the one above the horizon. Worked example at
lat = −85°, ARMC = 45°, ε = 23.44°:

| point               | ecliptic lon | RA   | δ       | meridian | altitude |
|---------------------|--------------|------|---------|----------|----------|
| raw MC (O/A/X)      | 47.464°      | 45°  | +17.04° | upper    | −12.04°  |
| flipped (C/R/T/Y/I) | 227.464°     | 225° | −17.04° | lower    | +12.04°  |

The two are exactly 180° apart. The flip touches the MC and cusp 10 only; the
Ascendant is computed independently from the unflipped ARMC and
hemisphere-corrected to the eastern horizon, so it is identical across all
house systems and unaffected by the flip — there is no associated ASC/DSC
inversion. For the flipped systems the reported MC equals cusp 10 (internal
consistency is preserved), but that cusp is the above-horizon lower-meridian
point rather than the geometric culmination.

*Geometric proof (not just parity):* Campanus and Regiomontanus are defined by
the family of great circles through the North and South points of the horizon;
cusp k is the ecliptic intersection of house circle k on the branch containing
the system's division point. Placing the reported (flipped) cusps back into
the true observer frame at lat −78° and −85° shows **all 12 cusps of both
systems lie exactly on their nominal house circles, on the nominal branches
(worst deviation 0.000000°)**: the flipped cusp 10 sits on the zenith branch
of the meridian (α = 90°), and the unflipped astronomical MC sits on the nadir
branch (α = 270°) — which is precisely where the system defines **cusp 4**,
and indeed the engine reports it there. The flipped output is therefore the
*unique solution of the systems' own defining equations*; reporting the
astronomical MC as cusp 10 in this regime would violate the definition of
these house systems. The backwards zodiacal ordering is likewise a theorem of
the geometry (the ecliptic crosses the circle family in reverse longitude
order there), not corruption. The one genuinely conventional residue is API
semantics: the `ascmc[1]` slot reports the quadrant cusp-10 anchor for
C/R/T/Y/I/J but the astronomical MC for O/A/X — inherited from the reference
API. (T and Y are not great-circle systems, so the proof covers C/R directly;
they anchor cusp 10 to the same above-horizon meridian point by
construction.)

*Verdict:* inherited convention, kept deliberately for 1:1 parity — and for
C/R proven to be the geometrically required branch under the systems' own
definitions. No default change is made. A future opt-in flag could expose the
purely astronomical upper-meridian MC in `ascmc[1]`, but only as opt-in, since
changing the default would break compatibility.

### 3. Time and Delta-T

**3.1 ΔT model.** See [Part A §2.3](#23-delta-t). Modern dates < 0.001 s;
1900–1950 ≤ ~1 s (measured max 0.98 s); >2050 depends on extrapolation (up to
~232 s at the extremes). Affects all UT↔ET/TT conversions (`deltat`,
`deltat_ex`, `utc_to_jd`, `jdet_to_utc`, `jdut1_to_utc`).

**3.2 Sidereal time (`sidtime`).** Modern < 0.001 s; >2050 up to ~0.05 s from ΔT
propagation into GMST.

**3.3 Equation of time (`time_equ`).** Modern < 0.15 s; >2050 up to ~1.4 s from
compounding ΔT and Sun-RA model differences.

### 4. Fixed stars

**4.1 Positions.** Ecliptic longitude/latitude agree within 0.01" for the 116 compared
stars (same FK5/Hipparcos proper motions).

**4.2 Distances.** Match exactly at J2000.0 for many catalog entries; away from
J2000 some entries differ because the reference API uses a different
radial-distance convention. The libephemeris distance channel is certified
against Astropy/SkyCoord; see
[Intentional divergences §7](intentional-divergences.md#7-fixed-star-distance-and-distance-speed-channels).

**4.3 Speed in distance (`speed_dist`).** Diverges 20–90% on some stars because
libephemeris uses a central finite difference on the full apparent geocentric
distance, while the reference API separates the radial-velocity component
differently. This is a certified semantic difference, not a copied convention.

**4.4 Magnitudes.** Agree within 0.5 mag (different catalog versions).

**4.5 Star-name resolution.** For a name absent from the catalog,
`fixstar2`/`fixstar2_ut` return an explicit "could not find star name" error rather
than a fuzzy-matched unrelated star. A handful of *traditional* names deliberately
resolve to a different star, following the IAU/WGSN-approved assignment:

| Name | libephemeris | pyswisseph | Basis |
|---|---|---|---|
| Menkar | α Cet | λ Cet | IAU name "Menkar" approved for α Ceti (2016) |
| Sadalbari | μ Peg | λ Peg | IAU name "Sadalbari" approved for μ Pegasi |
| Algedi | α² Cap | α¹ Cap | IAU name "Algedi" approved for α² Capricorni |
| Kurhah | ξ Cep | ζ Cep | IAU/WGSN name "Kurhah" approved for ξ Cephei (2016) |
| Mirak | β And (Mirach) | ε Boo | "Mirak" is a spelling variant of Mirach (β And); ε Boo is Izar |
| Algieba | γ Leo | γ¹ Leo | IAU/WGSN component (γ Leonis vs the γ¹ component) |
| Albireo | β Cyg | β¹ Cyg | IAU/WGSN component (β Cygni vs the β¹ component) |
| Almach | γ And | γ¹ And | IAU/WGSN component (γ Andromedae vs the γ¹ component) |

Genuinely ambiguous traditional names (no single authority; libephemeris keeps its
existing assignment): **Dheneb** α Cyg vs ζ Aql; **Ruc** δ Cas vs δ Cyg;
**Girtab** θ Sco vs κ Sco; **Deneb Kaitos** β Cet vs η Cet; **Ukdah** ι Hya vs
τ² Hya. Clearly-wrong cases (Alaraph → β Vir, Gienah Corvi → γ Crv, Atri → δ UMa,
Nash → γ² Sgr, Deli → η Aqr) were corrected to match both the naming tradition and
pyswisseph.

**4.6 Flamsteed designations: `fixstar` vs `fixstar2`.** For Flamsteed-style
designations (e.g. `29Psc`, `2Cet`) the legacy `fixstar` and the modern `fixstar2`
can resolve the same name to a different star (~289 of ~1450 entries). pyswisseph's
own `swe_fixstar` and `swe_fixstar2` likewise disagree. **`fixstar2` is the
reliable path** (e.g. `29Psc` → 29 Piscium). Proper names and Bayer designations
resolve identically in both (0.0"). Prefer `fixstar2` for Flamsteed lookups.

**4.7 Numeric star lookups (certified extension).** In `fixstar2`/`fixstar2_ut`/
`fixstar2_mag`, an all-digit search string (`"21421"`, `"HIP 21421"`, `",21421"`)
resolves by **Hipparcos number**; the reference API instead treats digits as a
sequential catalog index (the behaviour the legacy `fixstar` path keeps).
The HIP lookup is a deliberate extension: Hipparcos identifiers are stable
across catalog re-orderings while sequential indices are file-order-dependent.
Callers needing the reference's sequential-index semantics should use the
legacy `fixstar`. Comma-prefixed nomenclature (`",alTau"`) and every
non-numeric form resolve identically on both paths.

**4.8 Sidereal speed convention: legacy vs modern family (parity, per entry
point).** With `FLG_SIDEREAL | FLG_SPEED` the reference reports a *different*
first-coordinate speed on its two fixed-star families, and libephemeris now
mirrors each per entry point:

- **Modern** (`fixstar2` / `fixstar2_ut`): the true sidereal derivative — the
  ayanamsha drift rate d(ayanamsha)/dt is removed, so the reported longitude
  (or RA under `FLG_EQUATORIAL`) speed is the star's motion relative to the
  sidereal frame.
- **Legacy** (`fixstar` / `fixstar_ut`): the ayanamsha drift is **not** removed
  from the first-coordinate speed; the legacy speed carries the full
  mean-equinox motion and is therefore independent of the chosen ayanamsha.

Legacy speed = modern speed + d(ayanamsha)/dt, applied to the first coordinate
only (latitude/declination and distance speeds are identical on both families).
The offset is the *actual* ayanamsha rate at the epoch — ≈ the mean general
precession rate (~50.29″/yr at J2000) for precession-based modes such as
Fagan/Bradley and Lahiri, but the mode-specific rate (e.g. ~80.5″/yr at J2000)
for a fixed-target mode such as Galactic Centre 0° Sagittarius. Verified against
the reference oracle to < 0.001″/day across stars, ayanamshas and epochs.
Positions are identical on both families; only this speed channel differs.

### 5. Refraction and horizontal coordinates

**5.1 Atmospheric refraction (`refrac`).** Up to 15" near the horizon; both use
Bennett's formula with slightly different coefficients and boundary handling.

**5.2 Azimuth/altitude (`azalt`).** Above-horizon < 1"; below-horizon up to ~1654"
(~27') from fundamentally different refraction extrapolation for negative apparent
altitudes (physically meaningless region).

### 6. Eclipse and occultation functions

**6.1 Solar eclipses.** Timing (`tret[0]`) typically < 10 s; geography (`geopos`)
typically < 1°; type flags may differ for borderline cases. Total-eclipse
**obscuration** matches the reference (disc area ratio > 1 in totality) from the
reference-compatible entry points; a former clamp to 1.0 was removed — see
[Intentional divergences §2](intentional-divergences.md#2-total-eclipse-obscuration).

**6.2 Lunar eclipses.** Timing within ~10 s for most events.

**6.3 Lunar occultations.** Timing within ~0.001 day (~86 s) for most events.
`lun_occult_when_glob` is expensive (~6 s/call).

### 7. Nodes and apsides (`nod_aps_ut`)

Nodal/apsidal positions can diverge by 20–700" (up to ~1°) for some bodies, because
osculating elements are derived differently from the two engines (largest for outer
planets and high-eccentricity bodies).

**Lunar True Node / True Lilith.** Adjudicated against an independent osculating
computation from the DE441 Moon−Earth state vector (Ω = atan2(h_x, −h_y), h = r×v;
True Lilith = −eccentricity-vector direction). Over 1975–2022 libephemeris
reproduces this exact osculating truth to **0.000"** (both), while pyswisseph
differs by ~0.13–0.25" (True Lilith) — i.e. the residual is pyswisseph's
convention, not a libephemeris error. (The ~±30° month-scale libration of True
Lilith is the real physical swing of the instantaneous osculating apse.)

**True Node *distance* in LEB mode.** The True Node's longitude and latitude are
correct in every mode (0.000" across backends and vs the reference). Both
backends now agree on its **distance** as well: the LEB backend reconstructs the
geocentric lunar state from its Moon/Earth channels and evaluates the osculating
node radius at runtime (~0.0024 AU), matching the Skyfield/default backend to
~1e-11 AU. The reported distance rate (`ddist`) and `FLG_XYZ` output are the true
derivative of that reconstructed distance (a central difference over the same
±0.5 d window the reference smooths over), so speed and Cartesian output are
consistent across backends too.

**Minor bodies not yet supported (raises).** `nod_aps_ut`/`nod_aps` compute
nodes/apsides only for the Sun, Moon and the eight planets in this version. For
asteroids and centaurs (Chiron, Pholus, Ceres, Pallas, Juno, Vesta and any
numbered asteroid) the reference API computes real values, whereas this version
**raises `Error`** rather than returning a plausible-looking zero (which would
read as a node at 0° Aries). Planned for a future release.

**Sidereal `nod_aps` — uniform true-ayanamsha reduction (deliberate, ≤ ~17.5").**
For of-date sidereal output libephemeris subtracts the **true** ayanamsha
(mean + Δψ) uniformly, in `calc`/`calc_ut` AND in `nod_aps`/`nod_aps_ut` — the
two surfaces agree exactly (e.g. `nod_aps(MOON, MEAN)` node ≡
`calc(MEAN_NODE)` to 0.000" under `FLG_SIDEREAL`). The reference API is
internally inconsistent here: its `calc` subtracts the true ayanamsha
(measured shift 23.853203° at J2000, Lahiri) while its `nod_aps` subtracts
the **mean** ayanamsha (23.857073° — same date, same mode, all methods and
bodies). Sidereal `nod_aps` longitudes therefore differ from the reference by
exactly Δψ (nutation in longitude, ≤ ~17.5", −13.93" at J2000); tropical
`nod_aps` matches to < 0.01". libephemeris keeps the uniform reduction rather
than reproducing the reference's internal inconsistency.

### 8. Phenomena (`pheno_ut`)

Phase angle: inner planets < 1"; outer planets up to ~18" (position differences
amplified). Phase, elongation, apparent diameter and magnitude generally agree
< 1".

**By-number asteroids (limitation).** `pheno`/`pheno_ut` return real values for
the planets, the Moon, the Sun and the six curated minor bodies
(Chiron–Vesta). For other numbered asteroids (`AST_OFFSET + n`) and the exotic
LEB bodies this version returns zeros, whereas the reference API computes
phase/elongation (geometry) and magnitude (from its asteroid H/G data).
Universal H/G-based phenomena are planned for a future release (SBDB-sourced
photometric parameters).

### 9. Orbital elements (`get_orbital_elements`)

Inner planets (Mercury–Mars) agree within ~1". Outer planets (Jupiter–Neptune) can
diverge by 100–2000" in angular elements, because osculating elements derived from
instantaneous state vectors are highly sensitive for nearly circular orbits.

Adjudicated against an independent two-body extraction on the *same* DE441 state
(jplephem, library obliquity 23.4392911, GM = k²): across 1970–2020,
max |lib − independent two-body truth| = **0.000000°**, while max |SE − the same
truth| reaches **1.12°**. The lib-vs-SE difference is entirely Swiss Ephemeris's own
pipeline convention (state source / light-time / constants). Both engines use the
planet *barycenter* for the giants (the standard convention); using the planet
*center* instead would shift `varpi` by up to ~1.6° (Neptune) — the inherent
center-vs-barycenter ambiguity, not a defect in either engine. The lunar
nodes/apsides (10–14, 21–22) and the hypothetical bodies (40+) return zero
elements (no meaningful heliocentric orbit).

**Minor bodies.** The curated asteroids (Chiron, Pholus, Ceres, Pallas, Juno,
Vesta) and numbered asteroids (`AST_OFFSET + n`) return real osculating
elements derived from the same heliocentric state that serves their reported
positions (measured vs the reference API: Ceres/Vesta agree to ~1e-6 AU in
`a` and ~4e-5° in `M` at modern dates; centaurs within ~2e-4 AU, dominated by
the different orbit source). When no state source can serve the body, the
call raises a typed `Error` instead of returning silent zeros.
`orbit_max_min_true_distance` uses the same `a(1±e) ± a_Earth` approximation
for asteroids as for planets (internally uniform); the reference API computes
minor-body max/min slightly differently (~1%, e.g. 0.04 AU on Ceres), while
the "true" distance matches to ~1e-4 AU.

### 10. Sidereal calculations

**10.1 Sidereal positions.** Most planets agree within 1". **Sidereal Moon:**
3–14" — the sidereal frame amplifies the underlying lunar-theory difference because
the ayanamsha correction interacts with nutation differently in the two engines.

**10.2 Ayanamsha values.** Standard modes (Lahiri, Fagan-Bradley, Raman,
Krishnamurti, and the other fixed-epoch table systems) agree **exactly** (0.00" at
J2000, < 0.005" at ±100y). Exotic/experimental modes diverge and the divergence
grows with distance from J2000: ~0" at J2000 for star-anchored "True" modes rising
to ~40" at ±100y, up to ~145" for galactic/calculated modes at the extremes —
inherited from small fixed-star proper-motion / galactic-frame-definition
differences (§4). This is a definitional difference in those niche modes, not an
error: every fixed-epoch mode is exact.

A second, smaller family effect: libephemeris applies the same modern 5-term
IAU 2006 precession-rate polynomial to **all** formula-based ayanamsha modes,
while the reference API uses a different (older) accumulation model for a
group of classic/historical modes (De Luce, Suryasiddhanta variants, Aryabhata
variants, Kugler 1–3, Huber, Eta Piscium, Britton, Aldebaran-15Tau, J2000,
J1900, B1950, Lahiri-VP285, Krishnamurti-VP291). The residual is zero at
J2000, ~0.03–0.04" at ±100y, and reaches ~0.31" by 2500; its shape is
identical across the group (a pure precession-rate difference). Both
libephemeris backends agree with each other; the uniform modern model is
deliberate.

**10.3 SIDBIT projection flags (not yet supported).** The `SIDBIT_*` flags
(`ECL_T0`, `SSY_PLANE`, `USER_UT`, `ECL_DATE`, `NO_PREC_OFFSET`, `PREC_ORIG`; all
≥ 256) select alternative ecliptic/equinox projections that this version does not
implement. `set_sid_mode()` **strips them and emits a `UserWarning`**, keeping the
base ayanamsha mode — rather than the composite value silently falling back to
Lahiri. `SIDM_USER` (255) is unaffected. Planned for a future release.

**10.4 `SIDEREAL | EQUATORIAL` speed frame (parity).** Both engines report
the SID|EQ *position* on the mean equator of date and the accompanying
RA/Dec *speed* as the derivative of that reported position: measured against
the reference oracle, the Moon's SID|EQ speed agrees to 0.05"/day. An earlier
revision certified a divergence here based on a stale in-code note; direct
measurement disproved it. See `intentional-divergences.md` §9 (resolved).

### 11. Crossing functions

`solcross_ut`, `mooncross_ut`: typically < 1 s. `mooncross_node_ut`: up to ~69 s
from different lunar-node calculation methods.

### 12. Asteroid pipeline (`AST_OFFSET`)

`AST_OFFSET + N` remaps to dedicated body IDs through the Skyfield/SPK pipeline;
pyswisseph uses `.se1` files. Typical divergence ~0.2" for major asteroids.
**Missing .se1 files:** Chiron (2060) and Pholus (5145) require dedicated `.se1`
files in pyswisseph; if absent, pyswisseph raises an error, while libephemeris
always has these bodies via SPK.

### 13. Heliacal events

`heliacal_ut` is expensive for some configurations and depends on an empirical
atmospheric-extinction / arcus-visionis model.

The following **API-level behaviours now match the reference exactly** and are
verified against it:

- **Body/event acceptance matrix.** The Sun is rejected for every event; the
  Moon is accepted only for evening-first / morning-last (event types 3, 4) and
  rejected for rising/setting; Mercury and Venus accept all four events; the
  outer planets and fixed stars accept only rising/setting (1, 2). `vis_limit_mag`
  rejects the Sun and accepts every other body.
- **`heliacal_pheno_ut` illumination (element 27)** is 100.0 for fixed stars
  (point sources), not zero.
- **`vis_limit_mag` return flag** is the reference's binary encoding: −2 below
  the horizon, 0 photopic, 1 scotopic (never a separate "mixed" value).
- **Visibility-window ordering** (`jd1 ≤ jd2 ≤ jd3`) is guaranteed by
  construction; a failed crossing search clamps to the optimum or the search
  edge rather than an unrelated instant.
- **`refrac_extended` / `azalt`** reproduce the reference's below-horizon
  refraction curve and dip clamp (≤ ~10" over −3…+1° true altitude); see §refraction.

**Still divergent — the visibility magnitude model.** The sky-brightness,
limiting-magnitude, extinction (`kact`) and arcus-visionis (`minTAV`) model is a
library-specific empirical model, not the reference's Schaefer/Reijs pipeline.
Consequently the event *dates* (typically off by 1–27 days; e.g. Sirius heliacal
rising at Cairo 2024 returns 4 Aug vs the reference's 6 Aug), the visibility
*window widths* (often degenerate rather than the reference's few-minute spans),
the reported `kact`/`minTAV`, and `vis_limit_mag`'s limiting magnitude do **not**
match the reference. These remain recorded model differences for future
physical-model work (planned as a unified visibility engine reproducing the
Schaefer 1993 / Reijs extinction and threshold equations).

**13.1 Fictitious / Uranian bodies (Hamburg School, 40–48).** Cupido…Poseidon and
Transpluto are propagated from published Hamburg-School Keplerian elements. Both
engines place them on the **true ecliptic of date** — nutation in longitude (Δψ)
is applied uniformly, exactly as for every other body — so the remaining residual
is purely the difference in orbital-element sets, not a frame mismatch.
libephemeris vs the reference API now differ by up to ~25" in longitude
(most < 20"); Transpluto additionally carries up to ~25" in latitude from its
inclination/node elements. There is no independent astronomical reference for
these fictitious bodies, so neither is "truth".

Earlier libephemeris releases left bodies 40–48 on the *mean* ecliptic of date
(Δψ omitted), a ±17" split against the reference API's uniform nutation
treatment; that has been corrected in both backends (Skyfield path in
`planets.py`, LEB fast path in `fast_calc.py` Pipeline C), removing the whole
nutation term from this residual.

**13.2 Hypothetical bodies with `EQUATORIAL | J2000` (parity).** Both engines
rotate the J2000 ecliptic to the equator with the **J2000 obliquity** for
every body class including the hypotheticals (measured: the reference's
`EQ|J2000` matches its own eps_J2000 rotation to <0.01", vs 0.6-0.9" for an
of-date rotation at 2005). libephemeris briefly used a frame-mixed of-date
rotation for bodies 40-48 before v3 — that was a plain bug, fixed; an earlier
revision of this section mis-attributed the mixed frame to the reference.
See `intentional-divergences.md` §8 (resolved).

### 14. Constants and API

**14.1 Version string** intentionally differs (libephemeris reports its own
version). **14.2 `contrib` attribute** is not exposed by libephemeris. **14.3 `d2l`
with negative values** differs due to unsigned-integer overflow in the C
implementation vs Python integers. **14.4 `FLG_MOSEPH`** is accepted for API
compatibility but ignored — every calculation uses JPL DE440/DE441 via Skyfield,
with no Moshier fallback. **14.6 `cotrans_sp` speed for a non-physical `|lat| > 90°`** may differ from the
reference: the position (longitude/latitude) now matches for all inputs, and
the speed is exact for every physical `|lat| <= 90°`, but above 90° the two
engines reparametrise the (non-physical) point differently, so the returned
`lon`/`lat` speeds diverge. No real ecliptic/equatorial input reaches this
region. **14.5 `split_deg` with `SPLIT_DEG_NAKSHATRA |
SPLIT_DEG_ROUND_MIN`** returns the seconds field one arcsecond low (29 vs 30)
at a mathematically-exact nakshatra boundary (multiples of 13°20′: 40°, 160°,
200°…): the nakshatra width 360/27 is irrational in IEEE-754, so re-reducing
the rounded value reintroduces a sub-ULP deficit that truncates 29.9999″ to
29″. It is cosmetic (the minute/degree/nakshatra fields and the carry are all
correct) and unreachable with real ecliptic longitudes; the code keeps the
reduction that protects the far more common `90° + ROUND_MIN` case instead.

### 15. Gauquelin sectors

Typically agree within 0.01 sectors; small divergences (< 0.1) from the underlying
position and cusp differences.

### 16. pheno phase angle (attr[0])

The reference engine computes the Sun-body-Earth phase angle with an internal
vector recipe that black-box probing could not fully reproduce: no combination
of apparent / light-time-only / geometric positions matches it for every body
(apparent geocentric vectors match the Moon to <1" but leave 5-25" on planets
and asteroids; distance-based law-of-cosines behaves conversely). LibEphemeris
uses apparent geocentric vectors (Moon, planets) and the law of cosines on
distances (asteroids), leaving a residual of roughly 5-25 arcsec in attr[0]
(≲0.001 magnitudes through the H-G system, ≲1e-4 in illuminated fraction).
Elongation, diameter and magnitude agree to reference precision.

### Hyper-validation results

The hyper-validation script (`scripts/hyper_validate.py`) runs 4400+ comparison
rounds across 29 sections:

| Metric | Count | Percentage |
|--------|-------|------------|
| PASS | 3947 | 89.7% |
| KNOWN | 441 | 10.0% |
| FAIL | 0 | 0.0% |
| SKIP | 12 | 0.3% |

**0 failures.** All divergences are documented and classified as inherent
differences between the two computation engines. The 12 SKIP results are missing
`.se1` asteroid files in the pyswisseph configuration.

```bash
.venv/bin/python3 scripts/hyper_validate.py --json report.json
```
