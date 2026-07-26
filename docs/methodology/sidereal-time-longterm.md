# Long-term sidereal time, precession and obliquity

This note explains how libephemeris computes the **sidereal time (ARMC)** and the
**mean obliquity** that drive house cusps and the apparent-place reduction, and why
the chosen model stays correct over long time spans. The implementation is
`libephemeris/sidereal_longterm.py`.

## What these quantities are

A house chart is built from two geometric inputs:

* **ARMC** — the right ascension of the meridian, i.e. the *apparent sidereal
  time* at the observer's longitude. It fixes where the local meridian points on
  the celestial sphere.
* **the obliquity of the ecliptic ε** — the tilt between equator and ecliptic.

Both depend on **precession**, the ~25 800-year wobble of the Earth's axis that
moves the equinox. Sidereal time additionally depends on the Earth's rotation
(UT1) and therefore on **ΔT = TT − UT1**, the slowly varying, partly
unpredictable offset between uniform time and clock-from-rotation time.

## Why a long-term model is required

The classical IAU 1976 / IAU 2006 precession is a **polynomial in time** fitted
near J2000.0. A truncated polynomial is superb for a few centuries but **diverges
rapidly** outside its fit window: evaluated at ±8000 years it is wrong by
**degrees**, which corrupts every house cusp at historical or far-future epochs.
The same applies to a sidereal time read from an IAU-1976/2006 routine and an
obliquity read from the IAU-2006 polynomial — both are designed for the modern era
only.

This is not hypothetical for libephemeris: an earlier implementation derived the
house ARMC from a CIO-based IAU-2006 sidereal time, and the resulting house cusps
diverged from a long-term reference by up to **~3°** at ±8000 years. The current
model removes that error.

## The model we use

### Precession and obliquity — Vondrák 2011

We use the long-term precession of

> Vondrák, J., Capitaine, N. & Wallace, P. (2011), *New precession expressions,
> valid for long time intervals*, Astronomy & Astrophysics **534**, A22
> (DOI 10.1051/0004-6361/201117274; corrigendum A&A **541**, C1, 2012).

It is fitted to a numerical integration and stays accurate over **±200 000
years**, while agreeing with IAU 2006 to **sub-milliarcsecond** near J2000.0. So:

* modern results are **unchanged** (sub-mas agreement near J2000), and
* ancient/future results are **scientifically correct** (no polynomial blow-up).

The of-date **mean obliquity** is evaluated directly from the Vondrák obliquity
series. This single realization is shared by the house cusps **and** the
planetary/luminary position pipeline, so the bodies and the angles in one chart
sit in a single, self-consistent precession/obliquity frame — a chart can never
have its Sun and its Ascendant computed in mismatched frames.

### Sidereal time — a geometric construction, not a divergent polynomial

Sidereal time is the Earth-Rotation-Angle plus the accumulated precession in
right ascension. Evaluating the precession-in-RA as a polynomial diverges at
remote epochs for the same reason as above. We avoid that entirely with a
**geometric** construction that is stable everywhere:

1. take the mean longitude of the Earth from the published secular expression of
   Simon, J.L. et al. (1994), *Numerical expressions for precession formulae and
   mean elements for the Moon and the planets*, A&A **282**, 663;
2. correct it for Sun–Earth light-time;
3. form the corresponding unit direction on the ecliptic of J2000.0, rotate it to
   the J2000.0 equator, **precess it to the date with the Vondrák matrix**, and
   read off its ecliptic longitude of date;
4. add the equation of the equinoxes (nutation in longitude × cos ε) and the UT1
   hour angle.

Because step 3 transforms a *physical direction* through the long-term precession
matrix instead of summing a divergent RA polynomial, the result is stable over
the whole supported range. In the modern window 1850–2050 the well-established
IAU 2006 GMST polynomial is used directly (most precise there), with tiny
continuity offsets so the two branches join smoothly.

Time scales follow the standard convention: precession and obliquity are
evaluated at **TT**; the Earth-rotation hour angle uses **UT1**; the TT↔UT1
difference is the library's own ΔT, so houses and positions share one ΔT.

## Properties of this model

* **Validity range.** The model is correct to arcsecond level over ±200 millennia,
  so deep-historical charts (Babylonian, ancient-Egyptian) and far-future charts
  keep accurate house cusps, where a truncated IAU 1976/2006 sidereal time accrues
  degree-level errors.
* **Internal chart consistency.** House cusps and planetary positions use the
  *same* obliquity realization and the *same* ΔT. Using one model for positions and
  another for houses can place a body on the wrong side of a cusp at remote epochs
  even when each piece is individually "reasonable".
* **Verifiability.** Every coefficient comes from the cited peer-reviewed
  papers (Vondrák 2011 + corrigendum; Simon 1994) and the IAU 2006 GMST
  expression (Capitaine, Wallace & Chapront 2003, A&A 412, 567) — a direct,
  auditable implementation of published physics.

## The one honest caveat: ΔT at remote epochs

ΔT measures the Earth's rotation, which is **reconstructed** in the deep past
(from ancient eclipse records) and **unpredictable** in the future. libephemeris
uses the Stephenson–Morrison–Hohenkerk (2016) model plus observed IERS data,
applied consistently across the whole library (see
[delta-t.md](delta-t.md)).

At remote epochs the *sidereal time* is extremely sensitive to ΔT: the Sun's mean
longitude moves ~0.9856°/day ≈ 3548″/day, so a ΔT difference of even a fraction
of a day is amplified into arcminutes of ARMC. Consequently, comparing against any
engine that uses a different ΔT at ±8000 years can differ by arcminutes — but that
difference is **entirely the ΔT-model choice**, not a precession/obliquity error,
and it is:

* **zero** in the era of real use (1700–2300: sub-arcsecond), and
* negligible even for historical astrology (within a few arcseconds back to
  ~1000 BCE — far below the arcminute working precision of any astrologer).

When the same ΔT is used on both sides, the house cusps reproduce the
published-physics target to **< 0.05″ across the entire supported range** (and the
mean obliquity to < 0.001″) — i.e. the precession/sidereal-time model itself is an
exact, independent reproduction of the published physics; only the (physical,
unavoidable) ΔT choice remains.

## House cusp speeds (daily motion)

`houses_ex2` also returns the **speed** (daily motion) of each cusp and angle.
A cusp longitude λ is a function of time through the sidereal time, the obliquity
and (via the ecliptic frame) nutation, so its true speed is the **total time
derivative** dλ/dt. We compute it directly as a centered finite difference of the
full house solution in time,

    dλ/dt ≈ [ λ(jd + dt) − λ(jd − dt) ] / (2·dt),   dt = 2 seconds,

evaluated on `houses()` itself, so every time-dependent term is captured. The
2-second step is the measured optimum: the result is stable to ~10⁻³ °/day from
dt ≈ 30 s down to ≈ 1 s, while larger steps feel the cusp's curvature and much
smaller ones are dominated by floating-point noise.

This matters most for the **iteratively-solved** systems (Placidus, Koch). Their
intermediate cusps are defined by a transcendental (semidiurnal-arc) condition;
the genuine time derivative of that solution is what the finite difference above
returns, and integrating it reproduces the cusp's actual motion. An *analytic*
speed approximation of those systems can drift from the cusp's own motion — by a
fraction of a °/day at mid-latitude and by tens to hundreds of °/day near the
polar circle, where the iteration is most sensitive. The closed-form systems
(Regiomontanus, Campanus, Equal, Meridian) and the angles Asc/MC have no
iteration and are unaffected.

For **sign-locked** systems — Whole-Sign, Aries, Krusinski — the cusps sit at
fixed sign boundaries and jump discontinuously, so their instantaneous derivative
is meaningless. There we report the speed of the point that drives the wheel: the
Ascendant rate on cusps 1/7 and the MC rate on 4/10 (zero on the intermediates) —
the astrologically meaningful daily motion of the chart frame.

`houses_armc_ex2` is the one exception to the time-derivative method: it receives
only an ARMC (no Julian Day), so it differentiates with respect to ARMC and scales
by the sidereal rate — the most accurate speed obtainable from that input (the
missing obliquity-rate term is ~0.01 °/day). Callers that have the time should use
`houses_ex2`, which captures every term.

## Validation

Validation is source-based and reproducible in this repository:

- `tests/test_precession_vondrak.py` checks the published Vondrák/ERFA defining
  conditions and frame consistency.
- Sidereal and house tests verify periodicity, round trips, continuity,
  derivative consistency, finite behavior at high latitudes, and tier edges.
- Independent DE441 state vectors and ERFA frame routines provide a second
  implementation path without sharing LibEphemeris reduction code.
