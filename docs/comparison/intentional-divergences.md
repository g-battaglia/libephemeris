# Intentional compatibility divergences

LibEphemeris aims for 1:1 public-API compatibility. A divergence is retained
only when reproducing the observed behavior would require an unsupported
numerical definition, would violate the meaning of a documented output
channel, or would conflict with an independently sourced scientific model.

External reference-API comparisons are ephemeral. This page records semantics
and independent rationale, never per-date output, fitted thresholds, or
inferred internals.

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

## Total-eclipse obscuration is a bounded fraction

Obscuration — `attr[2]` of `sol_eclipse_how`, `sol_eclipse_when_loc`,
`sol_eclipse_where`, and `lun_occult_where` — is the fraction of the Sun's
disc area covered by the Moon, which is bounded to [0, 1] by its published
definition (NASA/USNO eclipse glossaries). During a total eclipse
LibEphemeris reports exactly `1.0`. The > 1 Moon/Sun disc *area ratio* that
some implementations place in this slot remains derivable as
`attr[1] ** 2`. Annular and partial eclipses are identical in both
conventions. The LibEphemeris-only helpers
(`sol_eclipse_obscuration_at_loc()`, `sol_eclipse_how_details()`) report the
same bounded fraction.

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
value is derived from compatibility output.

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

Every predefined ayanamsha derives from its author's published defining
statement propagated with the project's Vondrák Method-B model
([Sidereal modes](../reference/ayanamsha.md) lists the per-mode sources),
and the interpolated lunar apsides are anchored to the actual DE440 apsis
passages ([methodology](../methodology/interpolated-perigee.md)). Because
these values are derived from the published sources rather than tuned to any
external implementation, small numeric offsets from other engines are
expected: sub-arcsecond for many modes (Fagan/Bradley, Krishnamurti, the
galactic-center family on the Reid & Brunthaler 2004 Sgr A* position, the
vernal-point and frame modes), of order 10″–20″ where the official defining
pair is propagated with the project's precession model (the Lahiri family),
up to a few arcminutes for modes anchored directly to their author's book
value (e.g. Yukteshwar's published 20°54′36″ at the 1894 equinox), and of
order 0.01°–0.05° for `INTP_APOG` / `INTP_PERG` between apsis passages.
These offsets are properties of the independent derivation, not defects, and
they are stable and documented.

## LibEphemeris-only eclipse helpers

The LibEphemeris-only helpers `sol_eclipse_obscuration_at_loc()` and
`sol_eclipse_how_details()` expose a clearly named, physically bounded
covered-area fraction, matching the compatibility slots. Because these
functions are extensions, their contract is defined by LibEphemeris
documentation rather than by the reference API.

## Independent astronomical models

Planetary states, Delta T, precession/nutation, atmosphere, photometry,
eclipses, and lunar points use the independent models listed in
[Known differences](known-differences.md). When a reference API uses a
different valid convention, LibEphemeris documents the bounded model
difference rather than fitting its output.
