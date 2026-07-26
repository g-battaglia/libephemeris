# Missing primary-source data for historical hypothetical bodies

This page is the fail-closed inventory for historical compatibility IDs whose
numerical model could **not** be reconstructed from an independently reviewed
primary source. It is intentionally explicit: a recognised name or numeric ID
is not evidence that LibEphemeris has the right data, frame convention, or time
law needed to calculate it. See [Hypothetical bodies](hypothetical-bodies.md)
for the models that *are* implemented and their sources.

No number from a legacy compatibility table, another ephemeris distribution,
or an unexplained secondary transcription may be used to fill these gaps. The
IDs remain public API constants, but `calc_hypothetical_position()` raises
`UnknownBodyError` for them. This is preferable to silently returning a
scientifically untraceable position.

## What constitutes a recoverable model

For the bundled CSV and two-body propagator, a source-complete heliocentric
row needs all of the following information, either printed directly or derived
by transparent arithmetic from quantities printed in the same source:

1. the propagation epoch and its calendar/time-scale convention;
2. the reference equinox or other definition of the angular frame;
3. mean anomaly at the epoch, or an equivalent phase such as mean longitude
   together with the longitude of perihelion;
4. semimajor axis in AU, or an unambiguous period/distance law from which it
   can be derived;
5. eccentricity;
6. argument of perihelion, or an equivalent longitude of perihelion when the
   source explicitly defines a planar orbit;
7. longitude of the ascending node;
8. inclination;
9. the center of motion (heliocentric or geocentric); and
10. every secular term, rate unit, time variable, and circular-orbit phase
    convention needed to reproduce the author's position at other dates.

A source may explicitly define a circular or coplanar model, in which case
some geometrically degenerate angles can be fixed by that source's stated
convention. Missing angles are **not** silently replaced by zero. A sky
position at one date is also insufficient: infinitely many orbits pass through
one direction, and it does not recover the historical model.

## Summary

| ID | Public name | Source-complete row | Runtime behavior | Blocking evidence |
|---:|---|---|---|---|
| 49 | Nibiru | Missing | `UnknownBodyError` | No credible primary publication defining a complete numerical orbit has been identified. |
| 58 | Waldemath | Missing | `UnknownBodyError` | Waltemath's 1898 announcements publish distance, period, diameter, and transit dates only; no printed source supplies the phase, orientation, and equinox of the orbital convention associated with this ID. |

IDs 48 (Transpluto/Isis), 54 (Pickering), 55 (Vulcan), and 57 (Proserpina)
were recovered from their primary publications and moved to
[Hypothetical bodies](hypothetical-bodies.md).

## ID 49 — Nibiru

No primary astronomical or author-defined technical publication was found that
specifies a complete, reproducible Keplerian orbit for the compatibility object
called Nibiru. Narrative descriptions, modern web calculators, and numerical
tables that do not identify an original derivation do not establish provenance.

All model fields are therefore missing: epoch/time scale, equinox/frame, phase,
semimajor axis or period law, eccentricity, perihelion orientation, node,
inclination, center, and secular terms. This ID cannot be admitted merely by
finding numbers that reproduce another implementation; the defining primary
publication and its stated methodology are required.

## ID 58 — Waldemath

Georg Waltemath's 1898 announcements (summarized in *Science* 8, No. 189,
1898, p. 185, and retold by Sepharial in *The Science of Foreknowledge*)
publish only a mean distance of about 1.03 million km, a sidereal period of
119 days (177-day synodic), an apparent diameter, and predicted transit
dates. The Delphine Jay ephemeris (AFA, 1983) associated with this point
realizes a different, uniform-plus-sinusoid longitude model whose published
fit (about 1.55 deg amplitude over a 354-day term) is incompatible with an
eccentric Keplerian orbit.

| Required field | Status |
|---|---|
| epoch and time scale | An epoch (JD 2414290.95827875) circulates in secondary tabulations without a printed source |
| reference equinox / geocentric plane | Missing |
| phase and its exact rate law | Missing |
| orbital radius / semimajor axis | Derivable from the printed 1.03 million km distance |
| eccentricity | A value (0.1587) circulates in secondary tabulations without a printed source |
| periapsis orientation | Missing |
| node and inclination | An inclination (2.5 deg) circulates in secondary tabulations without a printed source |
| center | Earth-centered per the historical claim |
| additional periodic terms | Missing |

Admission requires a page-level primary source for the complete model —
candidate leads are Waltemath's own Hamburg pamphlets (1898) and the
Delphine Jay volume (AFA, 1983). Modern retellings or a numerical
reconstruction tuned to reproduce legacy positions are explicitly
insufficient.

## Procedure when a missing source is found

A future contribution must include all of the following before an ID can move
from `unsupported` to `primary-transcription` or another accurately named,
source-backed model status:

1. stable bibliographic metadata and a public source URL or library record;
2. page numbers for every transcribed field and the SHA-256 of the exact scan
   reviewed, when a downloadable scan is legally available;
3. a literal source table in code or data, preserving sexagesimal units and
   the author's time variable before any conversion;
4. separately documented arithmetic for calendar-to-JD, mean-longitude-to-
   anomaly, period-to-semimajor-axis, or frame conversion;
5. tests that pin the transcription and verify scientifically meaningful
   invariants without persisting reference-engine outputs;
6. an update to `scripts/check_hypothetical_provenance.py`, this file,
   `THIRD_PARTY_NOTICES.md`, and
   `docs/methodology/hypothetical-bodies.md`; and
7. review confirming that the source and the resulting data may be distributed
   under the project's stated terms.

Until then, absence from `libephemeris/data/fictitious_orbits.csv` is
intentional and must be preserved.
