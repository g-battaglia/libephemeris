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

`FLG_SPEED3` on the osculating lunar points (True Node, osculating apogee)
is the same case: the external coarse three-position value deviates from
the central difference of that implementation's own positions by ~10″/day
(True Node) up to several hundred ″/day (osculating apogee), and no two- or
three-point scheme at any step size reproduces it. LibEphemeris returns the
genuine derivative for `FLG_SPEED3` as well; for slow bodies (planets,
Moon, mean points) the two conventions agree.

## House-cusp speeds are the derivative of the reported cusps

For most systems the `houses_ex2` / `houses_armc_ex2` speed tuples are the
total time derivative of the reported cusp and angle functions — a centered
finite difference of the full house solution (ARMC rate, obliquity rate,
nutation, and the sidereal ayanamsha all included). One consequence is not
reproduced from external implementations:

- For intermediate cusps of iteratively solved systems (Placidus, Koch),
  an external analytic speed approximation can deviate from the derivative
  of that implementation's own reported cusps — by tens of degrees per day
  for Koch intermediate cusps (median ~27°/day measured) — even when the
  cusps themselves agree to sub-arcsecond level. LibEphemeris reports the
  genuine derivative: integrating it reproduces the cusp motion.
- Krusinski (`U`) is the extreme case: the external implementation reports
  0.0 for all eight intermediate cusp speeds while the true derivative
  reaches hundreds of degrees per day near high latitudes. LibEphemeris
  reports the genuine, self-consistent derivative here as well.
- Gauquelin (`G`, 36 sectors) belongs to the same class: sector-cusp
  positions agree to sub-arcsecond level, but the external analytic sector
  speeds deviate from the derivative of those same positions — by up to
  ~29°/day at latitude 60° (growing with latitude from a ~0.002°/day
  equatorial floor). LibEphemeris reports the genuine derivative.

Porphyry (`O`), Whole Sign (`W`) and Aries (`N`) are the exception: the
reference API reports an *analytic* cusp speed for these that is not the
derivative of its own cusp positions (its Porphyry second/fourth quadrants
are re-based on the ascendant speed, and its pinned Whole Sign / Aries cusps
expose the ascendant's speed in the 1/7 slots and the MC's in 4/10, with the
sidereal Whole Sign case exposing the ascendant speed in every slot).
LibEphemeris replicates that measured convention exactly for 1:1
compatibility.

## Phase-angle geometry in `pheno`

`pheno` / `pheno_ut` compute the phase angle from the **apparent geocentric
triangle**: both the body→observer leg (−P) and the illuminating body→Sun leg
(S − P) are formed from the apparent geocentric Sun and body position vectors
P, S that the same call already produces (both legs become geometric under
`FLG_TRUEPOS`). The illuminated fraction follows from that angle. The
elongation channel agrees with external implementations at the 0.01″ level,
and the Skyfield and LEB back ends agree with each other to below 0.01″.

The phase-angle channel differs from at least one external implementation at
the arcsecond level: across 1900–2100 the residual is bounded at roughly 19″
(Mercury, and the small-phase-angle outer planets Saturn and Pluto), about 7″
for Venus and 8″ for Mars, and below 1″ for the Moon. This is a
geometry-convention difference, not an ephemeris difference — no published
combination of apparent, astrometric, geometric, or single-/double-light-time
legs reproduces the external value for every body at once (the per-body best
fits contradict one another, and the external angle even responds to
aberration, shifting ~13–19″ when annual aberration is toggled off). The
apparent geocentric triangle is the documented convention chosen here; under
`FLG_TRUEPOS` the geometric triangle matches to <0.05″. The corresponding
difference in the illuminated fraction is below 1e-4.

The Sun itself reports 0.0 in all three phase channels (phase angle,
illuminated fraction, elongation): phase quantities are inapplicable to the
self-luminous disc.

## Heliocentric apparent-place light time

The apparent heliocentric place (`FLG_HELCTR` without `FLG_TRUEPOS`) retards the
target by the rigorous observer-to-target (Sun-to-target) light time: the
target state is re-evaluated at `t − lt` from the JPL ephemeris, the standard
IAU / NOVAS apparent-place reduction. The geometric heliocentric place
(`FLG_TRUEPOS`) and the reported heliocentric distance agree with an external
implementation to below 0.001″.

The apparent heliocentric *longitude* differs from that external
implementation by up to ~0.5″ for the fast inner planets (Mercury near
perihelion), ~0.06–0.3″ for Venus and Mars, and negligibly for the slow outer
planets. Its apparent displacement is a linear extrapolation along the
heliocentric velocity by an effective light time that departs from `r/c` by a
date-dependent amount — present even for a near-circular orbit, so it is not a
radial-velocity, Sun-retardation, or observer-aberration effect — which does
not reduce to a published closed form. The rigorous re-evaluated retardation is
the documented convention chosen here; no fitted per-date correction is carried.
Both the Skyfield and LEB back ends are identical.

## Outer-planet visual magnitude photometry

The `pheno` visual-magnitude channel differs from at least one external
implementation for three outer bodies. In every case the disk geometry
(diameter, phase, distance) is identical to the sub-arcsecond level and the
Skyfield and LEB back ends agree with each other exactly, so the difference is
purely a choice of photometric model, not an ephemeris difference.

- **Pluto** uses the modern Mallama & Hilton (2018), *Icarus* 306, 33
  photometry — V(1,0) = −1.024 mag with a linear phase coefficient
  β = 0.0362 mag/degree, i.e. V = V(1,0) + 5·log₁₀(r·Δ) + β·α. An external
  implementation uses the older flat photometry (V(1,0) = −1.00 mag, no phase
  term). Over the range of geocentric phase angles Pluto reaches the two
  disagree by up to ~47 mmag (worst near 1989); LibEphemeris keeps the modern
  phase-dependent model.

- **Uranus** uses the complete Mallama & Hilton (2018) photometric law — their
  Eq. 15 (*Astronomy and Computing* 25, 10; arXiv:1808.01973), the same work
  cited for Pluto and Saturn:

  V = 5·log₁₀(r·Δ) − 7.110 − 8.4×10⁻⁴·φ′ + 6.587×10⁻³·α + 1.045×10⁻⁴·α²

  Here φ′ is the photometric latitude — the average of the absolute
  planetographic sub-Earth and sub-solar latitudes. Those latitudes are the
  angles between Uranus's IAU (2009/2015) rotational pole (α₀ = 257.311°,
  δ₀ = −15.175°, J2000) and the Uranus→Earth / Uranus→Sun directions, each
  reduced to planetographic latitude with the flattening f = 0.0022927
  (Eq. 13). Because Uranus's poles are depleted in light-absorbing methane, the
  disk brightens as a pole turns toward the observer, which the φ′ term now
  models directly instead of the earlier single V(1,0). The sub-latitude term
  reproduces the ~42-year brightness half-cycle of the external model, but two
  structural differences remain against it: (a) the external model follows the
  *Astronomical Almanac* geocentric convention (Eq. 14) and drops the
  phase-angle polynomial, which for geocentric Uranus (α ≤ 3.1°) adds up to
  ~21 mmag and a ~+14 mmag mean bias; (b) the external sub-latitude realization
  is not reproducible without fitting. Over a full 1970–2059 orbit the residual
  is therefore bounded at ~46 mmag (rms ~25 mmag), larger than but of the same
  sub-latitude character as before. φ′ is evaluated by dotting the J2000 pole
  against the apparent ecliptic-of-date directions; the resulting frame mixing
  is ≤0.6 mmag, negligible beside the structural terms above.

- **Saturn** uses the full Mallama & Hilton (2018) formula with a Meeus ring-
  tilt reduction (the ring contribution dominates the disk near ring-plane
  crossings). The residual against the external model stays within ~2 mmag,
  concentrated near the low ring-opening epochs (e.g. ~2005), and reflects a
  ring-geometry convention rather than an ephemeris difference.

## Asteroid photometric data follow the current JPL SBDB

The H (absolute magnitude) and diameter values behind the asteroid
`pheno` photometry are transcribed from the JPL Small-Body Database
(retrieved 2026-07-13): for example Chiron H=5.54 and D=166 km. External
implementations that carry the older MPC-era values (Chiron H≈6.50)
differ by a constant magnitude offset per body — up to ~0.96 mag for
Chiron, ~0.14 mag Juno, ~0.05 mag Vesta — and a proportional apparent
diameter. The formula itself (IAU H-G system) is identical; only the
measured physical data vintage differs, and the current SBDB values are
kept.

## Rise/set events immediately after the search start

`rise_trans` returns the first event strictly following the start instant,
even when that event falls within minutes of it. At high latitudes near the
solstices (e.g. Reykjavík, where the midsummer sunset falls ~2 minutes
after local midnight) an external implementation can skip such an event and
report the following day's; LibEphemeris returns the real next event. When
comparing rise/set series externally, either start the search a few minutes
before the expected event or tolerate a one-day offset in this edge case.

The search is also complete at the far end: when the body stays below the
horizon for more than a day (the high-latitude Moon does so roughly once
per declination cycle), the next rise falls beyond the ~1-day forward
window some implementations scan, which then report a spurious
"circumpolar, no rise" result; LibEphemeris keeps searching and returns
the true event, however far ahead it lies.

## Meridian transits are solved to the geometric condition

`CALC_MTRANSIT`/`CALC_ITRANSIT` are refined until the body's apparent hour
angle is exactly 0°/180° to well under an arcsecond. At the instant
reported by at least one external implementation, the body's own hour
angle — evaluated with that implementation's functions — can still miss
the meridian: by 10–15 arcseconds (~0.6–1.0 s) for the fast-moving Moon's
lower transit, and by up to ~26 arcseconds (~1.8 s) for the upper transit
of a near-polar star (Polaris, growing as its declination approaches the
pole), where the flat culmination gives the residual nothing to push
against. LibEphemeris reports the geometrically exact instant rather than
reproducing the convergence residual.

## Lunar limb rise targets use the geocentric semidiameter

With `BIT_DISC_BOTTOM` (or the default upper limb) the Moon's limb
altitude target uses the semidiameter from the published IAU k = 0.2725076
and the geocentric distance; the event instant satisfies that target
exactly. An external implementation's limb target sits ~0.2–0.4" lower —
a composition of semidiameter and parallax-at-altitude that does not match
any canonical variant (geocentric, topocentric, or eclipse-k
semidiameter). The timing effect is below 0.15 s at ordinary latitudes and
reaches ~0.6 s only for grazing polar rises of the perigee Moon.

## Fixed-star speeds under `FLG_SPEED3`

For fixed stars, at least one external implementation returns zero speed
slots when only `FLG_SPEED3` is requested, while returning the real apparent
rates under `FLG_SPEED`. LibEphemeris reports the true derivative of its own
positions under both speed flags (verified against the numerical derivative
of both implementations' positions, which agree). The same
true-derivative-over-empty-slot choice documented for the lunar-apside
`SPEED3` case applies here.

## Robustness on degenerate inputs

Where the external implementation fails structurally on a legitimate or
degenerate input, LibEphemeris returns a defined result or a typed error
instead of reproducing the failure:

- `refrac_extended` with a below-sea-level observer (e.g. the Dead Sea)
  returns a zero horizon dip (the sea-horizon dip is not defined for a
  depressed observer) and applies the normal refraction; the external
  implementation emits NaN for the dip and silently skips the
  apparent-to-true correction.
- `rise_trans` with the latitude-zeroing bit (128) alone on a fixed star
  returns the stable zero-latitude projection event; the external result is
  non-deterministic (search-seed dependent, with round-hour sentinel
  values).
- A lunar occultation of the Moon itself raises the typed `Error`; the
  external implementation does not terminate.
- Event searches are stateless: a prior sequence of rise/transit calls at
  other locations (including polar cases that error) never changes a later
  result. In the external implementation such a sequence can corrupt a
  subsequent lunar meridian transit by ~31 s until the topocentric
  observer is explicitly reset.
- Geographic coordinates are validated against the geodetic domain before
  use: `set_topo` and the house functions raise the typed
  `CoordinateError` for finite out-of-range values (latitude outside
  [−90, 90], longitude outside [−180, 360]) as well as for NaN/infinite
  input. The external implementation accepts finite out-of-range values
  and returns the mathematically continued position (an observer latitude
  of −95° shifts the topocentric Sun by a few arcseconds); no geodetic
  datum defines an observer at such coordinates, so LibEphemeris rejects
  the input explicitly instead of computing from it.

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
- The `ECL_NUT` pseudo-body follows the ordinary `calc()`-versus-`calc_ut()`
  ephemeris-bit echo: the TT entry point echoes only the source bits the
  caller passed (a zero-flag request returns retflag 0), while the UT entry
  point injects the default `FLG_SWIEPH`. It is not treated as a special
  case.

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
