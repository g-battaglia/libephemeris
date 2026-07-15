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
| 48 | Transpluto / Isis | Missing | `UnknownBodyError` | The exact Strubell-attributed convention has not been recovered in a reviewable primary scan. |
| 49 | Nibiru | Missing | `UnknownBodyError` | No credible primary publication defining a complete numerical orbit has been identified. |
| 54 | Planet X Pickering | Missing | `UnknownBodyError` | Inspected Pickering publications give predictions/positions but not the complete element set needed by the runtime. |
| 55 | Vulcan | Missing | `UnknownBodyError` | The exact Weston-attributed numerical convention and its epoch/frame have not been recovered. |
| 57 | Proserpina | Missing | `UnknownBodyError` | No primary Abramov publication defining the complete numerical model has been recovered. |
| 58 | Waldemath | Missing | `UnknownBodyError` | Historical second-moon claims do not specify the complete computational reconstruction formerly associated with this ID. |

## ID 48 — Transpluto / Isis

### Bibliographic lead checked

Secondary descriptions attribute a Transpluto convention to an article by
Strubell in *Die Sterne*, issue 3 (1952), beginning around p. 70. During this
review, no openly inspectable primary scan or authoritative author-issued
transcription of the exact numerical element set was recovered. Because the
bibliographic details and the numbers could not be checked against the
publication itself, the secondary attribution is only a search lead, not a
source.

### Fields still missing or ambiguous

| Required field | Status |
|---|---|
| epoch and time scale | Missing |
| equinox / angular frame | Missing |
| mean anomaly or equivalent phase at epoch | Missing |
| semimajor axis or complete period law | Missing |
| eccentricity | Missing |
| argument/longitude of perihelion | Missing |
| ascending node | Missing |
| inclination | Missing |
| center of motion | Not independently confirmed from the primary publication |
| secular terms and units | Missing |

Admission requires a scan of the cited primary article, its complete
bibliographic metadata, and a field-by-field transcription showing the
source's epoch, frame, and time variable. A later secondary table is not a
substitute.

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

## ID 54 — Planet X Pickering

### Primary material inspected

- William H. Pickering, Harvard College Observatory Circular 215 (1919),
  [NASA ADS record](https://ui.adsabs.harvard.edu/abs/1919HarCi.215....1P/abstract),
  pp. 1–2.
- A public HathiTrust scan of Pickering material in the volume identified by
  HathiTrust ID `nyp.33433087548024`, including the relevant 1911 discussion.

The inspected material reports predicted sky information and discusses the
search for the trans-Neptunian object. Circular 215 points to a separate,
then-in-press memoir for the complete theory/elements; the inspected circular
does not itself supply a source-complete row. No separate memoir containing
all fields below was recovered during this review. Browser screenshots were
sufficient for this negative finding; no HathiTrust scan is bundled in the
repository.

### Fields still missing or ambiguous

| Required field | Status |
|---|---|
| propagation epoch and time scale | Missing as an element-set definition |
| equinox / angular frame | Missing as an element-set definition |
| mean anomaly or equivalent phase at epoch | Missing |
| semimajor axis or complete period law | Not recovered in a complete reviewed set |
| eccentricity | Not recovered in a complete reviewed set |
| argument/longitude of perihelion | Missing |
| ascending node | Missing |
| inclination | Missing |
| center of motion | Heliocentric context is apparent, but the exact model convention is not fully specified |
| secular terms and units | Missing |

Printed right ascension/declination or ecliptic position at a stated date must
not be inverted or fitted into substitute elements. Admission requires the
complete Pickering memoir (or another primary publication by Pickering that
prints the same complete model) and a documented literal transcription.

## ID 55 — Vulcan

Many nineteenth-century publications discuss proposed intra-Mercurial planets,
but that broad historical literature does not identify the exact numerical
**Weston** convention formerly associated with this compatibility ID. No
primary Weston publication containing its complete computation rule was
recovered and verified.

The missing items are the exact epoch/time scale, angular frame, phase,
semimajor axis or period law, eccentricity, perihelion orientation, node,
inclination, center convention, and any secular terms. Even if a historical
Vulcan orbit is found, it cannot be assumed to be the particular Weston model
without an explicit source chain connecting the author, edition, and fields.

## ID 57 — Proserpina

Secondary attributions connect this compatibility object with Abramov, but no
reviewable primary Abramov publication defining the precise model was
recovered. The author identity, publication title, edition/date, and page-level
location must be established before any numerical transcription can be
accepted.

All fields remain missing: epoch/time scale, frame/equinox, phase, semimajor
axis or period, eccentricity, perihelion orientation, node, inclination,
center, and secular terms. Numerically convenient circular defaults are not
evidence and must not be introduced.

## ID 58 — Waldemath

Historical literature reports Georg Waltemath/Waldemath's proposed second moon
of Earth, but a historical claim and a few observational predictions do not by
themselves recover the detailed computational convention formerly associated
with this API ID. No independently reviewed primary publication was found that
establishes every epoch, phase, rate, distance, and orientation needed for the
runtime reconstruction.

| Required field | Status |
|---|---|
| epoch and time scale | Missing |
| reference equinox / geocentric plane | Missing |
| phase and its exact rate law | Missing |
| orbital radius / semimajor axis | Missing as part of the exact reviewed convention |
| eccentricity | Missing |
| periapsis orientation | Missing |
| node and inclination | Missing |
| center | The historical second-moon claim implies Earth-centered motion, but the precise computational convention is unverified |
| additional periodic terms | Missing |

Admission requires a page-level primary source for the complete model. Modern
retellings or a numerical reconstruction tuned to reproduce legacy positions
are explicitly insufficient.

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
