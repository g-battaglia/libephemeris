# Keplerian Fallback — Possible Improvements

This document catalogs every known improvement opportunity for the Keplerian
fallback pipeline in `libephemeris/minor_bodies.py`. Items are organized by
estimated cost and expected precision gain.

## Table of Contents

- [Historical precision snapshot](#historical-precision-snapshot)
- [Current implementation summary](#current-implementation-summary)
- [Low-cost improvements](#low-cost-improvements)
- [Medium-cost improvements](#medium-cost-improvements)
- [High-cost improvements](#high-cost-improvements)
- [Structural issues (not precision improvements)](#structural-issues-not-precision-improvements)
- [Remaining investigation priorities](#remaining-investigation-priorities)
- [References](#references)

The Keplerian fallback is the last resort in the minor body calculation chain:

```
SPK kernel → Auto-download SPK → ASSIST n-body → Keplerian (this code)
```

It matters when SPK is unavailable AND ASSIST is not installed or its data
files are missing. The goal is to minimize error growth over time when this
fallback is the only option available.

---

## Historical precision snapshot

The table below is the original five-body benchmark that motivated the dense
multi-epoch work. It measured Ceres, Pallas, Juno, Vesta, and Chiron while only
six bodies had 50-year nodes. It is retained as dated engineering evidence,
not presented as the precision of the current ten-year, 36-body artifact.
Current measurements must be produced by
`tests/test_keplerian_precision_benchmark.py`, whose explicit body matrix now
contains all 37 configured bodies and skips only unavailable SPK truth inputs.

| Time from epoch | Max error | Dominant error source |
|-----------------|-----------|------------------------|
| At epoch        | 0.002"    | Floating point only |
| 1 month         | 7.5"     | Short-period perturbations |
| 6 months        | 49"      | Short-period perturbations |
| 1 year          | 1.9'     | Secular drift begins |
| 5 years         | 27'      | Secular + short-period |
| 10 years        | 45'      | Secular dominates |
| 25 years        | 2.6'     | Multi-epoch table helps |
| 50 years        | 3.6°     | Secular overwhelms |
| 100 years       | 3.5°     | Multi-epoch + secular |

---

## Current implementation summary

**37 bodies** have curated osculating elements: 13 main-belt objects, 6
centaurs, 9 TNOs, and 9 NEAs. Thirty-six use JD 2461000.5; Bennu uses the
separately documented OSIRIS-REx-derived heliocentric state at JD 2459000.5.

**Multi-epoch elements** cover **36 bodies** with 81 nodes each on the exact
ten-Julian-year grid 1650–2450. Bennu is the sole exception because no
compatible full-range heliocentric type-21 source kernel was available during
the reviewed generation. `scripts/generate_multi_epoch_elements.py` is the
authoritative generator; it now defaults to the exact shipped grid and rejects
non-solar centers or non-J2000 frames.

**Secular perturbations** from 4 planets (Jupiter, Saturn, Uranus, Neptune):
- omega (arg. perihelion): linear precession
- Omega (ascending node): linear regression
- e (eccentricity): oscillation via (h,k) vector formalism (Laplace-Lagrange)
- i (inclination): oscillation via (p,q) vector formalism (Laplace-Lagrange)
- n (mean motion): small, explicitly labelled project approximation coupled
  to the computed apsidal rate and eccentricity
- a (semi-major axis): constant, as in first-order secular theory

**Libration model** for 2 plutinos (Ixion, Orcus).

**Kepler equation solver**: Newton-Raphson with a Markley high-eccentricity
starter, tolerance 1e-8 and a 30-iteration elliptic cap; an independently
documented hyperbolic starter with a 50-iteration cap; and the Barker/Cardano
closed form for the parabolic limit. Both iterative solvers log non-convergence.

---

## Low-cost improvements

### L1. Denser multi-epoch table — completed at ten-year spacing

The initial 50-year, 17-node table was first proposed for reduction to twenty
years. The shipped artifact goes further: 81 nodes per covered body at ten-year
spacing from 1650 through 2450 inclusive. Away from the separately curated
present-day node, the nearest-node propagation distance is therefore at most
five Julian years.

The density is a LibEphemeris engineering choice, not a constant taken from a
publication. Every node is independently recomputed from a public JPL Horizons
type-21 position/velocity state by the equations in
`scripts/generate_multi_epoch_elements.py`; no compatibility output is fitted.

### L2. Multi-epoch coverage — completed for 36 of 37 bodies

The reviewed table expands the original six-body experiment to every
configured body for which compatible source coverage was available. Each of
the 36 lists has the same 81-node grid. Bennu deliberately remains outside the
dense table: its curated single node is derived from the identified
OSIRIS-REx/JPL kernel, while the generator fails closed rather than inventing
unavailable full-range states.

Dense osculating nodes reduce propagation distance but do not turn two-body
motion into a precision numerical integration. Close encounters and resonant
dynamics can still dominate, especially for NEAs and TNOs; the SPK and
ASSIST/REBOUND paths remain scientifically preferred when available.

### L3. Mean-motion drift approximation — implemented project choice

First-order Laplace-Lagrange secular theory keeps the semi-major axis constant,
so it does not itself supply a non-zero `d_a/dt`. The runtime nevertheless
applies a deliberately small project approximation

```
d_n = degrees(-1.5 * radians(d_omega) * e² / (1 - e²))
```

for `e < 0.99`. This couples the already computed apsidal rate to eccentricity
to limit phase drift. It is **not** transcribed from Murray & Dermott and must
not be represented as a literal second-order secular equation. The formula,
guard, and rationale are documented beside the implementation; changing it
requires a new SPK-backed validation, not an appeal to compatibility output.

### L4. Evolving planet elements — completed

`_calc_forced_elements()` accepts the target JD and evolves the four
perturbers' mean elements from J2000 with the published linear rates of Simon
et al. (1994). The static constants remain the defined J2000 intercepts. This
is a secular approximation with a documented validity limit, not a substitute
for evaluating the planets from JPL states.

### L5. Kepler solver convergence handling — completed

The elliptic and hyperbolic solvers now log the input and final residual when
their 30- or 50-iteration caps are reached. The elliptic branch also uses the
Markley (1995) high-eccentricity starter. No claim is made that a warning turns
a non-converged iterate into a valid result; callers needing stronger guarantees
should use the SPK or numerical-integration path.

### L6. Laplace coefficient integration density — completed

`_calc_laplace_coefficients()` evaluates the published integral with composite
trapezoidal quadrature: 200 panels normally and 500 when `alpha > 0.7`, where
the integrand becomes sharply peaked. The quadrature density is a project
accuracy choice. The defining integral and its `2/pi` normalization are cited
to Murray & Dermott section 6.4 beside the code.

---

## Medium-cost improvements

### M1. Short-period perturbation terms (analytical)

**Cost:** Medium-high — requires implementing perturbation series.
**Impact:** High at 1 month to 5 year timescales (currently 7-49" error).

This is the dominant error source at timescales < 5 years. Short-period
perturbations arise from conjunctions with Jupiter and Saturn and oscillate
with the synodic period (~400 days for Jupiter-Ceres).

**Theory:** Brouwer & Clemence (1961), Chapter 15. The first-order
short-period terms for the disturbing function give corrections to all
6 orbital elements as trigonometric series in the mean anomalies of the
asteroid and the perturbing planet.

**For the longitude of a main-belt asteroid, the dominant terms are:**

```
Δλ_short ≈ Σ_j A_j(a, e, i, a_J, e_J) × sin(k₁ M + k₂ M_J + k₃ ω + k₄ Ω)
```

where the amplitudes A_j depend on Laplace coefficients and eccentricity
functions. For Ceres, the dominant term has amplitude ~300" with period
~466 days (synodic period with Jupiter).

**What was tried and failed (Task 4.2):**
- Empirical Fourier fit over 800 years — secular drift dominated residuals
- Polynomial detrending + Fourier fit — correct amplitudes (~300" Ceres,
  ~780" Pallas) but WORSEN positions at short timescales because the
  amplitudes vary with time as elements drift. A static Fourier series
  captures the average amplitude which is wrong at any specific date.

**What would work (not yet attempted):**
- Analytical perturbation theory (Brouwer & Clemence): the amplitudes are
  functions of the current orbital elements, so they naturally evolve.
  Requires implementing the disturbing function expansion to second order
  in eccentricity and first order in inclination.
- Alternative: use the VSOP-like approach with precalculated coefficient
  tables per body from JPL data (Chapront-Touzé 1988 style).

**Estimated effort:** 2-4 days of implementation. ~200-400 lines of code.
Coefficient tables for Jupiter and Saturn perturbations on each body.

**Expected improvement:**
- 1 month: 7" → ~1"
- 6 months: 49" → ~5"
- 1 year: 1.9' → ~20"

### M2. Libration model for additional resonant TNOs

**Cost:** Medium — requires identifying resonances and fitting parameters.
**Impact:** High for specific bodies over decades.

Currently only Ixion (2:3) and Orcus (2:3) have libration corrections.
Other TNOs in mean-motion resonances with Neptune:

| Body | Resonance | Libration? | Currently modeled? |
|------|-----------|------------|-------------------|
| Ixion | 2:3 | Yes, ~78° amplitude | Yes |
| Orcus | 2:3 | Yes, ~68° amplitude | Yes |
| Haumea | 7:12 | Possibly | No |
| Makemake | Not resonant | N/A | N/A |
| Gonggong | 3:10 | Yes, ~30° amplitude | No |
| Eris | Not resonant | N/A | N/A |
| Varuna | Near 3:4? | Uncertain | No |
| Quaoar | Near 5:8? | Possible | No |
| Sedna | Not resonant | N/A | N/A |

**Changes required:**
- Identify which bodies are in confirmed resonances (from literature or
  numerical integration)
- Fit libration parameters (amplitude, period, phase) from SPK data
- Add entries to `PLUTINO_LIBRATION_PARAMS`
- The libration correction framework already exists (`calc_libration_correction()`)

**Expected improvement:** For resonant bodies, reduces error from degrees
to arcminutes over 100-year propagations.

### M3. Benchmark expansion — body matrix completed

`tests/test_keplerian_precision_benchmark.py` now enumerates all 37 configured
bodies, grouped as main-belt, centaur, TNO, and NEA objects. It compares only
when an independent SPK truth kernel is available and explicitly skips missing
truth inputs; therefore a green run is not evidence that every body was
measured unless the run report also confirms all kernels were present. The
historical five-body numbers at the top of this document remain labelled as
such until a controlled, fully provisioned 37-body matrix is archived.

### M4. Near-resonance rate scaling — heuristic implemented

The runtime identifies eight declared Jupiter/Neptune commensurabilities using
the dimensionless distance

```
delta = abs(p * n_body - q * n_planet) / n_body
```

and, only for `0.001 < delta < 0.05`, scales the non-resonant apsidal and nodal
rates by a capped project factor. Murray & Dermott chapter 8 supports the
physical warning that ordinary secular theory degrades near resonances; it
does **not** supply LibEphemeris's threshold, `1 / (20 delta)` scale, or caps.
Those numbers are transparent project heuristics and require empirical
validation against public SPK states. A true resonant Hamiltonian or numerical
integrator remains the scientifically preferred treatment.

---

## High-cost improvements

### H1. Second-order secular perturbation theory

**Cost:** High — substantial mathematical complexity.
**Impact:** Moderate improvement over first-order secular theory.

The current implementation uses first-order Laplace-Lagrange theory. The
second-order theory (Hori 1966, Yuasa 1973) includes:

- Cross-terms between different perturbing planets (Jupiter × Saturn
  coupling)
- Second-order eccentricity/inclination effects
- Long-period terms with periods of tens of thousands of years

These become important for:
- Bodies with high eccentricity (Hidalgo e=0.66, Icarus e=0.83)
- Bodies with high inclination (Pallas i=35°, Eris i=44°)
- Bodies crossing multiple planetary orbits (centaurs)

**Changes required:**
- Implement second-order averaging of the disturbing function
- Compute cross-coupling coefficients between planet pairs
- Add higher-order eccentricity functions d₁, d₂ to the Laplace
  coefficient calculation

**Estimated effort:** 3-5 days. ~500-800 lines of code.

**Expected improvement:** Factor of 2-3x at 100-year timescales for
high-e/high-i bodies. Negligible for low-e/low-i main belt asteroids.

### H2. Analytical short-period + long-period perturbation theory (full Brouwer)

**Cost:** High — complete Brouwer artificial satellite theory adapted for
heliocentric orbits.
**Impact:** High — would bring Keplerian to ~1" precision at 1-year scales.

Full implementation of Brouwer & Clemence (1961) perturbation theory for
heliocentric orbits:

1. **Short-period terms** (~synodic period): δa, δe, δi, δω, δΩ, δM
   as Fourier series in mean anomalies
2. **Long-period terms** (~secular, but with period ∝ 1/e): δω, δΩ
   corrections to secular rates
3. **Mixed secular-periodic terms**: amplitude modulation of short-period
   terms by secular evolution

This is what analytical orbit propagators (SGP4 for satellites, VSOP for
planets) implement. For asteroids, the key perturbers are Jupiter and Saturn.

**Changes required:**
- Disturbing function expansion in Legendre polynomials or Hansen
  coefficients
- Coefficient tables for each perturber-asteroid pair
- At evaluation time: compute all periodic terms and sum corrections

**Estimated effort:** 5-10 days. ~1000-2000 lines of code. Extensive
validation required.

**Expected improvement:**
- 1 month: 7" → < 0.5"
- 1 year: 1.9' → < 5"
- 10 years: 45' → < 1'

### H3. Semi-analytical propagation (Encke's method)

**Cost:** High — requires numerical integration at evaluation time.
**Impact:** Very high — ephemeris-quality for short arcs.

Instead of pure Keplerian + perturbation corrections, use Encke's method:
propagate the osculating orbit analytically (Kepler), then integrate the
perturbation residuals numerically with a low-order integrator (e.g., RK4
with large step size).

This combines the speed of Keplerian propagation (for the dominant 2-body
motion) with the accuracy of numerical integration (for the perturbative
residuals, which are small and smooth).

**Changes required:**
- Implement planetary position lookup (could use Skyfield or precomputed
  tables)
- RK4 or Störmer-Cowell integrator for the perturbation equations of motion
- Step size: ~10 days for main belt, ~30 days for TNOs
- Caching of intermediate results for repeated queries

**Estimated effort:** 3-5 days. But adds a runtime dependency on planetary
positions (Skyfield or LEB).

**Expected improvement:** Sub-arcsecond for propagation up to decades.
Essentially equivalent to ASSIST but without the ASSIST data files.

**Trade-off:** This blurs the line between Keplerian fallback and numerical
integration. If Skyfield is available (which it always is — it's a required
dependency), this approach would effectively be "lightweight ASSIST" without
the external data files.

### H4. Non-gravitational forces for NEAs

**Cost:** High — requires Yarkovsky/YORP modeling.
**Impact:** Critical for specific bodies (Apophis, Bennu, Ryugu) at
decade timescales.

Near-Earth asteroids experience measurable non-gravitational accelerations:

| Force | Magnitude (AU/day²) | Timescale | Bodies affected |
|-------|---------------------|-----------|-----------------|
| Yarkovsky | ~1e-15 | Decades | All NEAs |
| Solar radiation pressure | ~1e-13 | Days | Small NEAs (< 100m) |
| Outgassing | Variable | Perihelion | Cometary bodies |

For Apophis, the Yarkovsky effect shifts the orbit by ~200 km/year, which
translates to ~0.1" per year in geocentric longitude — small, but it
accumulates.

**Changes required:**
- Store Yarkovsky parameter (da/dt in AU/My) per body
- Apply as `a(t) = a₀ + (da/dt) × dt` in the Keplerian propagation
- Data source: JPL SBDB (https://ssd.jpl.nasa.gov/) provides measured
  da/dt for several NEAs

**Expected improvement:** ~1" per decade for NEAs with known Yarkovsky.
Irrelevant for main belt and TNOs.

### H5. Cubic Hermite element interpolation — tested and rejected

Cubic Hermite interpolation of osculating elements was prototyped, then
rejected because perturbation-driven element oscillations and angular wrapping
produced overshoot. The recorded tests degraded some TNO positions by up to
28× and main-belt positions by 3–4× at selected offsets. The production
consumer instead propagates the two nearest osculating nodes independently and
linearly cross-fades their **Cartesian positions** within a one-Julian-year
window centred on the epoch midpoint. This removes the step without claiming
that the elements themselves form a smooth interpolant.

---

## Structural issues (not precision improvements)

### S1. Solver convergence handling — remediated

Both iterative branches now emit a warning with the residual after exhausting
their iteration cap. A future bracketed fallback could strengthen this further,
but silent non-convergence is no longer the behavior.

### S2. Orbital-element sanity validation — remediated

`OrbitalElements.validate()` checks semi-major-axis sign, eccentricity/sign
consistency, inclination range, and bound-orbit mean motion. The runtime wrapper
logs every issue. These are diagnostic guards, not an orbit-determination
quality assessment.

### S3. Near-parabolic perturbation dispatch — remediated for bound orbits

Secular corrections are applied to bound and numerically near-parabolic inputs;
genuinely hyperbolic trajectories still bypass the bound-orbit perturbation
model by design. Applying Laplace-Lagrange terms to an unbound flyby would be
less defensible than exposing the limitation and preferring SPK/numerical
propagation.

### S4. Co-orbital detection too aggressive — remediated

The former implementation skipped a perturber whenever
`|a - a_planet| < 0.1 AU`. That fixed project heuristic could suppress terms
for Trojan or horseshoe objects and had no scale dependence. The current
implementation instead uses the public Hill-radius approximation
`r_H = a_planet * (mu_planet / 3)^(1/3)` and disables linear secular theory
inside that planet-specific domain. This is a validity guard: close/co-orbital
dynamics violate the small, non-resonant perturbation assumptions of
Laplace-Lagrange theory and require SPK or N-body propagation.

---

## Remaining investigation priorities

### If the goal is maximum evidence for minimum effort:

1. Provision independent SPK truth for every body in the existing 37-body
   benchmark and archive a dated result matrix.
2. Quantify the project `d_n` and near-resonance heuristics by ablation against
   those public states; remove either heuristic if it does not improve robustly.
3. Add a bracketed fallback for the rare non-converged Kepler solve.

### If the goal is pushing Keplerian to its theoretical limit:

1. **M1** — Analytical short-period perturbations (eliminates 7-49" errors)
2. **H2** — Full Brouwer theory (sub-arcsecond at 1 year)
3. **L4 follow-up** — validate evolved mean elements at century scale

### If the goal is making the Keplerian fallback irrelevant:

1. **H3** — Semi-analytical Encke propagation (sub-arcsecond without ASSIST)

This would effectively replace the Keplerian fallback with a lightweight
numerical integrator that uses Skyfield's planetary positions directly.
The advantage over ASSIST: no extra 1.3 GB download. The disadvantage:
slower than pure Keplerian (but still much faster than Skyfield per-body
evaluation).

---

## References

- Brouwer, D. & Clemence, G.M. (1961). Methods of Celestial Mechanics.
  Academic Press. — Chapters 11-16 (perturbation theory).
- Murray, C.D. & Dermott, S.F. (1999). Solar System Dynamics. Cambridge
  University Press. — Chapters 6-8 (secular and resonant dynamics).
- Hori, G. (1966). Theory of general perturbation with unspecified
  canonical variable. Publ. Astron. Soc. Japan, 18, 287.
- Markley, F.L. (1995). Kepler equation solver. Celestial Mechanics, 63, 101.
- Raposo-Pulido, V. & Pelaez, J. (2017). An efficient code to solve the
  Kepler equation. MNRAS, 467, 1702.
- Simon, J.L. et al. (1994). Numerical expressions for precession formulae
  and mean elements for the Moon and planets. A&A, 282, 663.
