# Long-term sidereal time, precession and obliquity

This note explains how libephemeris computes the **sidereal time (ARMC)** and the
**mean obliquity** that drive house cusps and the apparent-place reduction, why
the chosen model is the scientifically correct one over long time spans, and why
it is better than the approach most ephemeris/astrology engines use. The
implementation is `libephemeris/sidereal_longterm.py`.

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

## The problem with the usual choice

The classical IAU 1976 / IAU 2006 precession is a **polynomial in time** fitted
near J2000.0. A truncated polynomial is superb for a few centuries but **diverges
rapidly** outside its fit window. Evaluated at ±8000 years it is wrong by
**degrees**, which corrupts every house cusp at historical or far-future epochs.
Many astrology engines inherit exactly this limitation, because they read
sidereal time from an IAU-1976/2006 routine and obliquity from the IAU-2006
polynomial — both designed for the modern era only.

Concretely, before this work libephemeris derived the house ARMC from a CIO-based
IAU-2006 sidereal time; the resulting house cusps diverged from a long-term
reference by up to **~3°** at ±8000 years.

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

## Why this is better than the competitors

* **Validity range.** The model is correct to arcsecond level over ±200
  millennia. Engines built on IAU 1976/2006 sidereal time are correct only for a
  few centuries and accrue degree-level house errors at deep-historical dates
  (Babylonian, ancient-Egyptian charts) and far-future dates.
* **Internal chart consistency.** House cusps and planetary positions use the
  *same* obliquity realization and the *same* ΔT. An engine that uses one model
  for positions and another for houses can place a body on the wrong side of a
  cusp at remote epochs even when each piece is individually "reasonable".
* **Independence and verifiability.** Every coefficient comes from the cited
  peer-reviewed papers (Vondrák 2011 + corrigendum; Simon 1994) and the IAU 2006
  GMST expression (Capitaine, Wallace & Chapront 2003, A&A 412, 567). Nothing is
  fitted to, or copied from, any other software's output or source code — it is a
  direct, auditable implementation of published physics.

## The one honest caveat: ΔT at remote epochs

ΔT measures the Earth's rotation, which is **reconstructed** in the deep past
(from ancient eclipse records) and **unpredictable** in the future. libephemeris
uses the Stephenson–Morrison–Hohenkerk (2016) model plus observed IERS data,
applied consistently across the whole library.

At remote epochs the *sidereal time* is extremely sensitive to ΔT: the Sun's mean
longitude moves ~0.9856°/day ≈ 3548″/day, so a ΔT difference of even a fraction
of a day is amplified into arcminutes of ARMC. Consequently, comparisons against
another engine at ±8000 years can differ by arcminutes — but that difference is
**entirely the ΔT-model choice**, not a precession/obliquity error, and it is:

* **zero** in the era of real use (1700–2300: sub-arcsecond), and
* negligible even for historical astrology (within a few arcseconds back to
  ~1000 BCE — far below the arcminute working precision of any astrologer).

When the same ΔT is used on both sides, the house cusps reproduce the reference
model to **< 0.05″ across the entire supported range** (and the mean obliquity to
< 0.001″) — i.e. the precession/sidereal-time model itself is an exact,
independent reproduction of the published physics; only the (physical,
unavoidable) ΔT choice remains.

## Validation

The model is exhaustively validated in the separate `validation/` repository:

* `validation/compare_scripts/rounds/houses/lt*.py` (run with
  `sh run.sh houses`) — 25 scripts covering, across −13200…+17191: mean-obliquity
  model-purity (< 0.001″), sidereal-time model-purity at matched ΔT (< 0.05″),
  every house system at many latitudes (matched-ΔT cusps < 0.05″), vertex /
  auxiliary points, cusp speeds, high latitudes, a dense modern no-regression
  sweep, sidereal (ayanamsha) cusps, internal consistency, DE441 edges,
  round-trips, `house_pos`, and a massive systems×latitudes×epochs grid.
* `validation/compare_scripts/rounds/positions/lt*.py` (run with
  `sh run.sh positions`) — confirms the shared obliquity leaves modern positions
  unchanged and bounds the (smooth, expected) of-date obliquity shift at remote
  epochs.
