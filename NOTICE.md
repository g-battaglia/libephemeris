# NOTICE — Provenance and Intellectual Property

## Independent Implementation

LibEphemeris is an **independent implementation** of an astronomical ephemeris
library for Python. Its computations are independently implemented and are not
derived from the Swiss Ephemeris (SE) source code by Astrodienst AG.

In June 2026 an internal provenance review found that a small cluster of
house-system routines in `houses.py` (Savard-A, Krusinski, APC, the
Sunshine/Makransky variant, and parts of `house_pos`) carried
implementation expression that followed the structure and identifiers of
SE's `swehouse.c` more closely than this project's independence standard
allows. Those routines were rewritten in place: identifiers, comments and
docstrings were re-derived from the published system definitions (Savard;
Krusinski 1995 / Pisa 1997 / Goelzer 1995; Knegt; Makransky 1988; Holden
1977; Koch), while the underlying mathematics — which is not protectable
subject matter — was preserved and verified bit-identical over a ~91,000
case grid against the pre-rewrite outputs. The full sweep record lives in
`docs/methodology/provenance-sweep-2026-06.md`.

All astronomical computations are based on:

- **JPL DE440/DE441** planetary and lunar ephemerides, accessed via the
  [Skyfield](https://rhodesmill.org/skyfield/) library (Brandon Rhodes, MIT license)
- **IAU 2006/2000A** precession-nutation model, via
  [pyerfa](https://github.com/liberfa/pyerfa) (BSD-3-Clause license)
- **Peer-reviewed academic sources**, including but not limited to:
  - Meeus, J. — *Astronomical Algorithms*, 2nd ed. (1998)
  - Simon, J.L. et al. — "Numerical expressions for precession formulae
    and mean elements for the Moon and the planets", A&A 282, 663 (1994)
  - Chapront, J. et al. — "A new determination of lunar orbital parameters,
    precession constant and tidal acceleration from LLR", A&A 387, 700 (2002)
  - Chapront-Touze, M. & Chapront, J. — ELP 2000-85 lunar theory
  - Park, R.S. et al. — "The JPL Planetary and Lunar Ephemerides
    DE440 and DE441", AJ 161, 105 (2021)
- **Primary historical sources** for hypothetical bodies:
  - Witte, A. & Sieggrun, F. — *Regelwerk fur Planetenbilder* (1928)
  - Neely, J. — refined orbital elements (1988)
  - Makransky, B. — *Primary Directions* (1988), for the Sunshine house system

## API Compatibility

LibEphemeris provides an API that is **signature-compatible** with
[pyswisseph](https://github.com/astrorigin/pyswisseph) (the Python binding
for Swiss Ephemeris). Function names follow the canonical bare-name form
of the upstream reference API (e.g., `calc_ut`, `houses`, `julday`), with
matching parameters and flag constants to allow drop-in migration.

API compatibility does not imply code derivation. The underlying algorithms,
data sources, and implementation are entirely independent. API signatures
and interface conventions are not copyrightable subject matter
(see *Google LLC v. Oracle America, Inc.*, 593 U.S. 1 (2021); EU Directive
2009/24/EC art. 1.2).

## Extended Astrology Helpers (`contrib` submodule)

The `libephemeris.contrib` submodule provides extended astrology helpers
(zodiac sign and nakshatra constants, aspect angles, Vedic dignity
calculations, longitude-to-rasi conversions, antiscion arithmetic, etc.).

These are **independently reimplemented in pure Python** from public
astrological tradition and standard zodiacal geometry. Specifically:

- Sign and nakshatra divisions follow the standard zodiacal partitions
  (12 × 30° and 27 × 13°20') documented in any classical astrology text.
- Vedic dignity tables (exaltation, debilitation, lordship, naisargika
  friendships) derive from traditional Vedic astrology (Parashara,
  Varahamihira) and are in the public domain.
- Aspect angles (conjunction, sextile, square, etc.) are geometric
  constants of pure mathematics.

No source code from the upstream pyswisseph contrib submodule (the C
``swephelp`` library) was consulted or copied during implementation.
API names and signatures match the upstream reference only at the
interface level, which is not copyrightable subject matter.

## Development History

During early development, some experimental branches temporarily included
data from Swiss Ephemeris sources (e.g., Moshier trigonometric tables,
`seorbel.txt` orbital element file). These were identified, removed, and
replaced with independently sourced alternatives before any stable release:

- Moshier analytical backend: removed entirely in favor of JPL DE440/DE441
  via Skyfield (no analytical approximations are used in production)
- `seorbel.txt`: replaced with `data/fictitious_orbits.csv`, compiled from
  primary published sources (Witte/Sieggrun 1928, Neely 1988)
- All algorithms reference peer-reviewed publications or JPL data products
  as their primary sources

The git history of this repository reflects this progression transparently.

## Vendored and Adapted Components

Three files shipped in the package are vendored or adapted third-party
code and keep their upstream licenses, all permissive (MIT):
`vendor/spktype21.py`, `moon_theories/tass17.py`, and the generated
periodic-term tables in `moon_theories/tass17_data.py`. The package
contains no copyleft code. Full inventory and license texts:
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

The Galilean satellite module (`moon_theories/galilean.py`) was, through
v2.1.0, adapted from PyMeeus and licensed LGPL-3.0. In June 2026 it was
rewritten clean-room from the published theory (Lieske 1998, A&AS 129,
205; Meeus 1998, *Astronomical Algorithms*, ch. 44) under strict
information barriers: a functional specification restating the published
mathematics and the project's frame adaptations was written from the
theory, and an independent implementation was produced from that
specification alone — without reference to PyMeeus or the prior module.
The rewrite reproduces the prior numeric output to floating-point
re-association level (sub-nanometre per moon component over 1800–2200) and
is now owned by the project under the dual license. Process and
independence record:
[docs/methodology/galilean-clean-room-2026-06.md](docs/methodology/galilean-clean-room-2026-06.md).

## Calibration Data Disclosure

Two generated data sets are calibrated against pyswisseph used strictly as
a **black-box oracle**: the residual tables in
`libephemeris/lunar_apse_corrections.py` (interpolated lunar apogee and
perigee, bodies INTP_APOG / INTP_PERG) and the trigonometric perturbation
coefficients in `lunar.py` produced by the `leph calibrate` workflow
(regenerable via `scripts/generate_lunar_apse_corrections.py`). These
tables contain numeric program *output* — computed positions of the lunar
apsides — not Swiss Ephemeris source expression; the fitting pipeline and
all runtime code are original. The INTP_* bodies are constructs defined by
the reference API, so 1:1 behavioral parity requires fitting to reference
output. This is disclosed for transparency and is part of the legal review
checklist for the commercial edition (see
[COMMERCIAL-LICENSE.md](COMMERCIAL-LICENSE.md)).

A third, smaller item: the **legacy Uranian display tables** in
`libephemeris/hypothetical.py` (`URANIAN_ELEMENTS` and the per-body
`*_KEPLERIAN_ELEMENTS` `L0` / Hades `M0` constants) are historical
pyswisseph-oracle-calibrated fits. They are **not consulted at runtime** —
the live Keplerian propagation uses the published Witte/Lefeldt 1928 +
Neely 1988 elements from `libephemeris/data/fictitious_orbits.csv` — and
are retained only for module-API stability (tests pin them). The
code↔CSV↔publication provenance is enforced by
`scripts/check_hypothetical_provenance.py` (`poe provenance:hypothetical`)
and documented in `docs/methodology/hypothetical-bodies.md`.

## AI-Assisted Development

Parts of this codebase were developed with AI assistance (Anthropic
Claude), always under the direction and review of the project author, with
behavior verified against published references and the black-box oracle
harness. Under Anthropic's commercial terms, rights in such output belong
to the customer; copyright in the work is held by Giacomo Battaglia.

## Trademarks

"Swiss Ephemeris" is a product name of Astrodienst AG. References to it in
this repository are nominative only: they describe interface compatibility
and black-box verification targets. LibEphemeris is not affiliated with,
or endorsed by, Astrodienst AG.

## License

LibEphemeris is **dual-licensed**: under the **GNU Affero General Public
License v3** (AGPL-3.0-only — see the [LICENSE](LICENSE) file) or,
alternatively, under a commercial license available from the copyright
holder (see [LICENSING.md](LICENSING.md) and
[COMMERCIAL-LICENSE.md](COMMERCIAL-LICENSE.md)). Owned source files carry
the SPDX expression `AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial`;
vendored/adapted modules keep their own identifiers.

This project has no license dependency on, and no license obligation toward,
the Swiss Ephemeris or Astrodienst AG.
