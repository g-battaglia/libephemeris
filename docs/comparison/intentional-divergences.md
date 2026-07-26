# Intentional compatibility divergences

LibEphemeris aims for 1:1 public-API compatibility. A divergence is retained
only when reproducing the observed behavior would require an unsupported
numerical definition, would violate the meaning of a documented output
channel, or would conflict with an independently sourced scientific model.

This page records the semantic choices and the independent rationale behind
each retained divergence.

## `SIDEREAL | J2000` is honored uniformly for lunar points

`FLG_J2000` combined with `FLG_SIDEREAL` is applied to **every** lunar point
(`MEAN_NODE`, `TRUE_NODE`, `MEAN_APOG`, `OSCU_APOG`, `INTP_APOG`,
`INTP_PERG`) with one pipeline: the mean-ecliptic-of-date state is precessed
to the J2000 ecliptic, then the mean ayanamsha is subtracted in that frame.

Some external implementations silently discard `FLG_J2000` for the non-mean
lunar points under this combination while honoring it for the mean points.
That asymmetry is not reproduced, for physical reasons: ayanamsha and
ecliptic-plane precession are independent, composable rotations, and
discarding one of them for half of a body family makes the separation
between the true and mean node grow without bound away from J2000 — far
beyond the ±1.5° physical oscillation of the true node around the mean node.
LibEphemeris keeps |TrueNode − MeanNode| physically bounded at every epoch
under every flag combination.

## Of-date mean obliquity from the precession poles

The public of-date mean obliquity is the angle between the Vondrák (2011)
ecliptic-pole and equator-pole series — the same two poles that build the
long-term precession matrix. This keeps the precession and every
equator↔ecliptic-of-date rotation one self-consistent frame: a direction
lying in the mean ecliptic of date (the Sun, by definition) reduces to ~0
ecliptic latitude at every epoch.

An implementation that pairs pole-based precession with the separately
fitted direct `ε_A` series tilts its of-date ecliptic away from its own pole
by up to ~6.5″ at −3000, which surfaces as a spurious ecliptic latitude on
the Sun of the same size. LibEphemeris does not reproduce that
inconsistency. The two realizations agree to < 0.001″ across 1900–2100, so
modern-era output is unaffected; the direct series remains available as
`sidereal_longterm.mean_obliquity_series_rad`.

## `deltat_ex` is pure

`deltat_ex(jd, flag)` applies the flag-selected tidal acceleration
functionally, to the returned Delta-T value only. `get_tid_acc()` changes
exclusively through `set_tid_acc()`; no Delta-T query mutates library state.
An implementation whose Delta-T query rewrites the value observable through
its `get_tid_acc()` is exposing a side effect, not a documented result
channel, and that side effect is not reproduced. Only introspection through
`get_tid_acc()` after a `deltat_ex` call can observe the difference; every
returned Delta-T value matches.

## Total-eclipse obscuration slot carries the disc-area ratio

Obscuration — `attr[2]` of `sol_eclipse_how`, `sol_eclipse_when_loc`, and
`sol_eclipse_where` — reports the Moon/Sun disc *area ratio*
(`attr[1] ** 2`), which exceeds 1 during a total eclipse (measured
compatibility convention). Annular and partial eclipses are identical in
both conventions, and occultation slots stay at 1.0. Consumers who want the
physically bounded covered-area fraction in [0, 1] (the NASA/USNO glossary
definition) use the LibEphemeris-only helpers
(`sol_eclipse_obscuration_at_loc()`, `sol_eclipse_how_details()`), which
keep that contract.

## Sripati `house_pos` preserves the fraction in the house-12 wrap

`house_pos` with the Sripati system shifts the Porphyry placement forward by
half a house, which lands part of the house-12 region beyond 13.0 before
wrapping. LibEphemeris wraps only the true overflow back to [1.0, 1.5),
preserving the fractional position, so a body geometrically inside house 12
reports as 12.x — consistent with the engine's own Sripati cusps for the
same chart. A collapse of the whole wrap region onto exactly 1.0 loses the
fraction and contradicts the cusps returned for the same chart, so it is not
reproduced. Consumers that need the collapsed convention can apply
`1.0 if hpos > 12.0 else hpos`.

## Mean node/apsis speed semantics

For a MEAN node or apsis request, LibEphemeris keeps each speed slot as the
time derivative of the corresponding position slot. Derivatives are computed
from the independently defined mean geometry with a documented
finite-difference scheme when an analytical derivative is unavailable.

If an external implementation exposes a position-like quantity in a speed
slot, LibEphemeris does not reproduce that defect. No correction or special
value is derived from compatibility output. The same principle applies when
an external speed slot is not the derivative of that implementation's own
reported position: measured externally, MEAN `nod_aps*` speed slots can
deviate from the numerical derivative of the same implementation's position
series by degrees per day, while every LibEphemeris speed slot integrates
back to its position slot exactly.

## House-cusp speeds are the derivative of the reported cusps

The `houses_ex2` / `houses_armc_ex2` speed tuples are the total time
derivative of the reported cusp and angle functions — a centered
finite difference of the full house solution (ARMC rate, obliquity rate,
nutation, and the sidereal ayanamsha all included). Two consequences are
not reproduced from external implementations:

- For intermediate cusps of iteratively solved systems (Placidus, Koch),
  an external analytic speed approximation can deviate from the derivative
  of that implementation's own reported cusps by ~0.1–0.7°/day even when
  the cusps themselves agree to sub-arcsecond level. LibEphemeris reports
  the genuine derivative: integrating it reproduces the cusp motion.
- Whole Sign cusps are pinned to sign boundaries, so their time derivative
  is exactly 0.0 and LibEphemeris reports 0.0. An external convention that
  exposes the ascendant's speed in the cusp-1/7 slots (and the MC's in
  4/10) for pinned cusps is not reproduced; those rates remain available
  in the `ascmc_speed` tuple.

## Phase-angle geometry in `pheno`

`pheno` / `pheno_ut` compute the phase angle from the light-time-consistent
Sun–body–observer geometry, and the illuminated fraction follows from that
angle; the elongation channel agrees with external implementations at the
0.01″ level. The phase-angle channel differs from at least one external
implementation by a bounded 15–40″ for the inner planets and Mars (a
geometry-convention difference, not an ephemeris difference — no public
combination of apparent/geometric distances reproduces the external value
exactly). The difference in the illuminated fraction is below 1e-4.

The Sun itself reports 0.0 in all three phase channels (phase angle,
illuminated fraction, elongation): phase quantities are inapplicable to the
self-luminous disc.

## Rise/set events immediately after the search start

`rise_trans` returns the first event strictly following the start instant,
even when that event falls within minutes of it. At high latitudes near the
solstices (e.g. Reykjavík, where the midsummer sunset falls ~2 minutes
after local midnight) an external implementation can skip such an event and
report the following day's; LibEphemeris returns the real next event. When
comparing rise/set series externally, either start the search a few minutes
before the expected event or tolerate a one-day offset in this edge case.

## Fixed-star radial channels

For fixed stars, `dist` is the geocentric distance of the same propagated
star state used for the angular coordinates, and `dist_spd` is its time
derivative. This keeps position/speed channel semantics consistent across
stars and planets. The result is validated against the independently sourced
catalogue, Astropy space-motion geometry, and the JPL Earth state.

## Node positions without transient spikes

Planetary node/apsis points apply the aberration correction analytically, so
the reported ascending and descending nodes are exactly antipodal at all
times and their speeds stay smooth. Numerical per-point correction schemes
that occasionally produce isolated single-day jumps with an order-of-
magnitude speed spike are not reproduced. Validation comparisons of node
positions should either use `FLG_TRUEPOS` or tolerate such spikes in the
external data.

## Independent derivative and frame conventions

Several channels use one uniform, independently defined convention even when
another public implementation exposes a different edge-case result:

- `FLG_SPEED3` is the standard centered finite difference of the reported
  position, with angular unwrapping across the cyclic boundary.
- Every `nod_aps*` speed component is the physical time derivative of its
  corresponding reported node or apsis component.
- `FLG_HELCTR` and `FLG_BARYCTR` use the requested NASA JPL center and
  observer-to-target light time; they do not reuse a geocentric light-time
  convention.
- `fixstar*` and `fixstar2*` both differentiate the state in the frame
  actually returned to the caller, so position and speed cannot refer to
  different frames.
- `calc()` and `calc_ut()` apply the same default return-flag convention to
  the `ECL_NUT` pseudo-body.

These choices follow coordinate and derivative definitions, not retained or
fitted external observations.

## Numerically independent data models

Every predefined ayanamsha has an explicit defining pair or live geometric
anchor propagated with the project's Vondrák Method-B model. The evidence
record is deliberately per-mode: it distinguishes primary publications and
public catalogue/frame geometry from secondary historical attributions and
disclosed project conventions
([Sidereal modes](../reference/ayanamsha.md) gives the complete status table).
The interpolated lunar apsides are separately anchored to physical DE440 apsis
passages ([methodology](../methodology/interpolated-perigee.md)). Because no
value is tuned to external implementation output, small numeric offsets from
other engines are expected: sub-arcsecond for many modes (Lahiri on the Indian
Astronomical Ephemeris J2000 constant, the galactic-center family on the Reid &
Brunthaler 2004 Sgr A* position, and geometric frame modes), of order 10″–25″
for some historical-pair modes, up to a few arcminutes for Yukteshwar's primary
book value, and of order 0.01°–0.05° for `INTP_APOG` / `INTP_PERG` between
apsis passages. A numeric result can therefore be reproducible without being
misrepresented as a direct author transcription.
These offsets are properties of the independent derivation, not defects, and
they are stable and documented.

## LibEphemeris-only eclipse helpers

The LibEphemeris-only helpers `sol_eclipse_obscuration_at_loc()` and
`sol_eclipse_how_details()` expose a clearly named, physically bounded
covered-area fraction, matching the compatibility slots. Because these
functions are extensions, their contract is defined by LibEphemeris
documentation rather than by the reference API.

## Topocentric Moon reduction is internally consistent

The topocentric observer vector uses the full nutation-consistent frame
chain, and topocentric velocities are the exact time derivatives of the
reported topocentric positions. External implementations that build the
observer vector with a simplified rotation differ by up to ~0.16″ in the
Moon's topocentric position and by a few arcseconds/day in its topocentric
speed — and are there internally inconsistent with their own trajectory
(their reported speed does not differentiate their reported position).
LibEphemeris keeps the self-consistent reduction.

## Sidereal time is continuous across 2050

`sidtime`/`sidtime0` follow the ERFA GST06A realization at every date. At
least one external implementation switches its internal sidereal-time model
exactly at 2050-01-01, producing a ~0.13 s step (about 1.9″ of ARMC) that
propagates into its rise/set times after that date. LibEphemeris does not
reproduce the step: its sidereal time stays continuous and matches ERFA,
so a constant ~0.13 s offset against that implementation is expected for
dates from 2050 on.

## Independent astronomical models

Planetary states, Delta T, precession/nutation, atmosphere, photometry,
eclipses, and lunar points use the independent models listed in
[Known differences](known-differences.md). When a reference API uses a
different valid convention, LibEphemeris documents the bounded model
difference rather than fitting its output.
