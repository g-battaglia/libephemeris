# Precision Summary

Compact precision summary for LibEphemeris vs pyswisseph 2.10.03.
For full details, models, and methodology see [reference/precision.md](reference/precision.md)
and [reference/swisseph-comparison.md](reference/swisseph-comparison.md).

## Planetary Positions (geocentric ecliptic, 1550–2650 CE)

| Body | Mean Diff | Max Diff | Notes |
|------|-----------|----------|-------|
| Sun | 0.04" | 0.20" | DE440 vs DE431 |
| Moon | 0.70" | 3.32" | Numerical vs analytical lunar theory |
| Mercury | 0.05" | 0.32" | |
| Venus | 0.08" | 0.33" | |
| Mars | 0.06" | 0.58" | |
| Jupiter | 0.12" | 0.44" | Includes COB correction |
| Saturn | 0.13" | 0.51" | Includes COB correction |
| Uranus | 0.23" | 0.50" | |
| Neptune | 0.24" | 1.17" | |
| Pluto | 0.26" | 0.75" | Includes COB correction |

All planets sub-arcsecond. Moon ~3" max reflects different lunar models (JPL DE440 numerical integration vs ELP/MPP02 + DE431).

## Velocities

| Component | Max Diff |
|-----------|----------|
| Longitude speed | < 0.003°/day |
| Latitude speed | < 0.004°/day |
| Distance speed | < 0.0001 AU/day |

### House cusp speeds (`houses_ex2`)

The cusp/angle speeds are the genuine time derivative dλ/dt of each cusp,
computed from the full house solution, so the reported speed integrates to the
cusp's real motion. Validated against an independent high-accuracy derivative
(no other engine in the loop):

| System group | Reported speed vs true dλ/dt |
|--------------|------------------------------|
| Angles (Asc/MC) + closed-form (Regiomontanus, Campanus, Equal, Meridian) | < 0.005°/day |
| Iterative (Placidus, Koch) | < 0.005°/day |

For the iterative systems this is **more accurate than an analytic speed
approximation**: such an approximation deviates from the cusp's own motion by a
fraction of a °/day at mid-latitude and by tens to hundreds of °/day near the
polar circle (e.g. Koch at 60° latitude). Sign-locked systems (Whole-Sign,
Aries, Krusinski) report the guiding-point (Asc/MC) rate on the angle cusps, the
astrologically meaningful daily motion of the chart frame.

## Lunar Points

| Point | Max Diff | Independent Verification |
|-------|----------|--------------------------|
| Mean Node | < 0.001° | — |
| True Node | < 0.01" | Verified vs JPL Horizons to machine precision |
| Mean Lilith | < 0.015" (lon) | Latitude ~20" systematic (different node formulas) |
| True Lilith | < 0.5" | Both libraries ~240" from Horizons (inherent two-body limit) |
| Interpolated Apogee | ~0.36° | Genuine algorithm difference (JPL DE440 vs ELP2000-82B perturbation series) |
| Interpolated Perigee | ~2.6° | JPL DE440 physical passages vs truncated ELP2000-82B perturbation series |

## House Cusps

< 0.02" for all 24 supported house systems, tested at 11 global locations.
Iterative systems (Placidus, Koch) use 10⁻⁷° convergence threshold.

**Long-term range (where LibEphemeris is more accurate).** House cusps are
driven by sidereal time (ARMC) and obliquity — both functions of precession.
LibEphemeris derives them from the long-term Vondrák 2011 model (valid ±200,000
years) through a geometric construction that stays stable everywhere, rather than
an IAU 1976/2006 sidereal-time polynomial that diverges by **degrees** outside a
few centuries. Verified across the full ephemeris range (−13200…+17191 CE):

| Quantity | Criterion | Result |
|----------|-----------|--------|
| Mean obliquity | vs reference, full range | < 0.001" |
| Sidereal time (ARMC) | vs reference, matched ΔT, full range | < 0.05" |
| Cusps, all systems | vs reference, matched ΔT, full range | < 0.05" |
| Cusps, modern era (1860–2040) | vs reference, own ΔT | < 0.05" |

At remote epochs a *whole-chart* comparison against another engine differs by the
ΔT-model choice (amplified through the Sun's mean longitude, ~3548"/day), which is
physical and shared by both house cusps and planetary positions — not a cusp
error. With the same ΔT the cusp model reproduces the reference exactly across the
whole range. See
[sidereal-time-longterm.md](methodology/sidereal-time-longterm.md).

## Fixed Stars

116 stars from Hipparcos catalog with van Leeuwen 2007 proper motions.
Max difference: 0.51" (Rigil Kentaurus — nearest star, parallax not modeled).
98% of 101 comparable stars within 0.5". Two catalog bugs found and fixed
(Algedi wrong component, Asellus Borealis wrong HIP number).

## Ayanamsha

- Standard modes (Lahiri, Fagan-Bradley, Raman): < 0.0002°
- Star-based modes (True Citra, True Revati): < 0.006°

## Eclipses

- Solar eclipse timing: < 6 seconds
- Lunar eclipse timing: < 8 seconds
- Rise/set/transit: < 30 seconds

## Delta T

- Modern (1900–2025): < 1 second
- Historical (< 1700): up to ~187 seconds (different models: SMH 2016 vs E&M 2006)
- Future (> 2050): grows with extrapolation divergence

## Minor Bodies

- With SPK kernels: sub-arcsecond (matching JPL Horizons)
- Keplerian fallback: ~10–30" for main belt asteroids near epoch, degrees over decades

## Hypothetical Planets

Uranian hypothetical planets: < 1" (Keplerian from published elements).

## Heliocentric / Barycentric / Equatorial / XYZ

| Mode | Max Diff |
|------|----------|
| Heliocentric | < 0.0004° (1.1") |
| Barycentric (non-Sun) | < 0.001° |
| Equatorial RA/Dec | < 0.0005° (1.7") |
| XYZ Cartesian | < 0.00005 AU |

## Hyper-Validation

4400+ comparison rounds across 29 API sections:
**3947 PASS, 441 KNOWN, 0 FAIL, 12 SKIP**.
All divergences documented in [divergences.md](divergences.md).
