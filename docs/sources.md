# Scientific sources and provenance

The following sources define models and data used by the native engine.
Citations identify origins, not universal accuracy. External components retain
[their licenses and notices](../THIRD_PARTY_NOTICES.md); scientific authorship
does not replace redistribution terms.

## Models and data

| Area | Defining sources | Implementation / retained input |
| --- | --- | --- |
| Planetary ephemerides | Park et al. (2021), *The JPL Planetary and Lunar Ephemerides DE440 and DE441*, [AJ 161, 105](https://ssd.jpl.nasa.gov/doc/Park.2021.AJ.DE440.pdf) | [SPK reader](../core/src/providers/spk/), selected external JPL files |
| Classic DE files | JPL [format](https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/ascii_format.txt) and [reading guide](https://ssd.jpl.nasa.gov/ftp/eph/planets/fortran/userguide.txt) | [Native classic-DE reader](../core/src/bodies/nbody/native/classic_de.rs) |
| Nutation and Earth orientation | IAU 2000A/2000B, IAU 2006; Simon et al. (1994), A&A 282, 663–683; [IERS Conventions 2010, Chapter 5](https://iers-conventions.obspm.fr/content/chapter5/icc5.pdf) | [ERFA module](../core/src/time/erfa.rs), licensed SOFA/ERFA-derived tables |
| Long-term precession | Vondrák, Capitaine & Wallace (2011), A&A 534 A22; corrigendum (2012), A&A 541 C1 | ERFA long-term pole series; [canonical sidereal model](../core/src/time/canonical_sidereal.rs) |
| Delta-T | Stephenson, Morrison & Hohenkerk (2016), historical rotation reconstruction; IERS observations and Skyfield numerical snapshots | [SMH2016 evaluator](../core/src/time/smh2016.rs), [table](../core/src/tables/smh2016_tables.rs) |
| Stars | ESA Hipparcos, XHIP, cross-identifications and IAU/WGSN names | [Catalog snapshot](../tools/generators/data/star_catalog.json), [compiled table](../core/src/tables/star_catalog_tables.rs) |
| Lunar mean quantities | Meeus, *Astronomical Algorithms*, 2nd ed. (1998), ch. 47; Chapront-Touzé & Chapront (1991); IAU 2012 Resolution B2 | [Lunar models](../core/src/bodies/lunar.rs); AU = 149,597,870.7 km exactly |
| Lunar apsides | DE440-fitted periodic/residual model identified in its table header | [Lunar apse tables](../core/src/tables/lunar_apse_tables.rs) |
| Saturnian satellites | TASS 1.7, later IMCCE revision of Vienne & Duriez (1995), A&A 297, 588–605 | [Snapshot](../tools/generators/data/tass17.json), Gajdosik/Stellarium MIT transcription; the paper describes TASS 1.6 |
| Refraction | Green (1985), eq. 4.25; ICAO atmosphere; Sæmundsson/Meeus altitude correction and Thom horizon dip | [Refraction](../core/src/phenomena/refraction.rs); [visibility convenience models](../core/src/phenomena/schaefer.rs) distinguish published equations from heuristics |
| Planetary photometry | Mallama & Hilton (2018), *Computing apparent planetary magnitudes for The Astronomical Almanac*; IAU 2015 physical constants | [Phenomena implementation](../core/src/phenomena/) |

Catalog astrometry, static magnitudes and variable-star brightness are different
quantities. The compiled catalog contains 1,447 stars; Spica uses Hipparcos H5
Johnson V = 0.98. Mean lunar points, osculating geometry and interpolated
passages likewise represent distinct constructions.

## Numerical snapshots and regeneration

Generators are in [`tools/generators/`](../tools/generators/). Retained snapshot
bytes and generated headers record exact input identities; regeneration must
not silently update sources or production ephemerides.

- `gen_nutation_tables.py` reads ERFA 2.0.1 sources from PyERFA 2.0.1.5:
  `nut00a.c`, `nut00b.c` and `s06.c` supply the 678 luni-solar/687 planetary
  terms, abbreviated 77-term series and CIO blocks. Headers retain SHA-256s.
- `gen_iau2000_tables.py` uses the separate
  [IAU snapshot](../tools/generators/data/iau2000_nutation.npz), SHA-256
  `db231254da3bcbef24faa994eadb80f378d39460aa407026c156e6627c30f11a`.
  Its Skyfield MIT notice remains applicable. This input carrier is not a
  Python runtime backend.
- Long-term pole inputs come from ERFA revision
  [`9915ba38c9365f8b0738269b8c2ac1fdd5f8dee3`](https://github.com/liberfa/erfa/tree/9915ba38c9365f8b0738269b8c2ac1fdd5f8dee3/src).
- `gen_smh2016_tables.py` retains the selected Skyfield 1.54 historical/IERS
  realization. Table dates are TT and Delta-T values are seconds; public
  `deltat` returns days. Observed rows and extrapolated predictions differ.
- `gen_star_catalog_tables.py` serializes the immutable catalog snapshot;
  `star_catalog_builder.py` is a separate acquisition tool, not a requirement
  to reproduce that snapshot. Catalog refreshes require separate review.
- `gen_tass17_tables.py` serializes 1,535 terms in 299 lists from the retained
  TASS snapshot, preserving its attribution.

The bundled LEB image's identity and stored units are documented in
[the format specification](data-format.md). External SPKs retain their own
segment centre, frame, timescale and coverage; inspect the actual file rather
than inferring these from its name.

## Native N-body engine

The engine is independently authored Rust; historical REBOUND/ASSIST API names
are compatibility labels, not runtime implementations. Its equations come from:

- Rein & Spiegel (2015), MNRAS 446, 1424–1437,
  [§2.1 equations 1–7 and §2.3](https://arxiv.org/html/1409.4779v2):
  Gauss–Radau collocation. `gen_native_radau.py` constructs the basis from
  Legendre polynomials, with independent Decimal100/140 rounding checks.
- Pham, Rein & Spiegel (2024), [equations 6–7](https://arxiv.org/html/2401.02849v1):
  A85 timescale from acceleration, jerk, snap and crackle norms. The project's
  sampling/anchoring policy is explicit, not an upstream implementation copy.
- Wisdom & Hernandez (2015), [universal-variable initial-value formulation](https://arxiv.org/html/1508.02699v1):
  Cartesian Kepler drift. Conic initial values use independent high-precision
  mathematical reference construction, not fitted engine outputs.
- [Published hybrid Hamiltonian partitions](https://academic.oup.com/mnras/article/533/3/3708/7735374):
  the exposed solar maps have the supported stationary-Sun/massless-target
  scope, not general collision or interacting-particle capabilities.
- Holman et al. (2023), [III.2–III.4](https://arxiv.org/html/2303.16246v1),
  and Park et al. (2021), III.1–III.2: PPN/EIH, nongravitational forces and
  zonal potentials. Constants come from selected files; source sets and pole
  policies are explicit. Production uses a fixed ICRF-north terrestrial axis
  and the IAU/WGCCRE solar axis, not a claim that Earth's physical axis is fixed.
- JPL [SB441 memorandum IOM392R-21-005](https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/SB441_IOM392R-21-005_perturbers.pdf):
  N16 trajectories are a subset of the mutually perturbed N373 solution,
  not a closed sixteen-body integration. Use consistent planetary/small-body
  generations and account for differing force models.

Current sources do not bundle REBOUND/ASSIST implementations. Earlier releases
retain the notices applicable to the code they included. The publications above
identify the mathematical basis of the independently authored Rust engine;
they do not imply that the cited software implementations were copied or
adapted. Copied or adapted third-party material remains subject to its own
terms. See [limitations](limitations.md#models-with-additional-scope) when
interpreting propagated trajectories.

## Independent implementation policy

Implement from identified publications, catalogs, project conventions and
appropriately licensed components. Separate published coefficients from
project interpolation, convergence and fallback choices. Source-incomplete
hypothetical models fail closed; names alone do not define an orbit.

Do not inspect, retrieve, translate, adapt or copy Swiss Ephemeris source,
comments, documentation, algorithms or data; introduce a wrapper; or fit model
constants/fixtures to another engine's outputs. Ephemeral public-API behavioral
comparisons must not persist Swiss outputs, residuals or fixtures. Explain
significant discrepancies with independent JPL observations or physical
invariants, and never weaken tolerances or substitute data to obtain agreement.
Compatibility terminology implies no affiliation or endorsement.

Preserve upstream notices with copied/adapted licensed material. Source
attribution is not a legal clearance or a physical-accuracy certificate.
[Contributing](../CONTRIBUTING.md) explains the local development workflow.
