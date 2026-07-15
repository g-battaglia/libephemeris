# NOTICE

LibEphemeris is an independent implementation of an astronomical ephemeris
library for Python, licensed under the **Apache License, Version 2.0** (see
[LICENSE](LICENSE) and [LICENSING.md](LICENSING.md)).

## Data and Model Sources

All astronomical computations are based on:

- **JPL DE440/DE441** planetary and lunar ephemerides, accessed through the
  local [Skyfield](https://rhodesmill.org/skyfield/) path (Brandon Rhodes, MIT
  license), precomputed LEB coefficients, or NASA JPL Horizons state vectors
  ([bundled base-core build attestation](docs/leb/base-core-provenance.md))
- **IAU 2006/2000A** precession-nutation and **Vondrák 2011** long-term
  precession, via [pyerfa](https://github.com/liberfa/pyerfa) (BSD-3-Clause)
- **Peer-reviewed academic sources**, including:
  - Meeus, J. — *Astronomical Algorithms*, 2nd ed. (1998)
  - Simon, J.L. et al. — A&A 282, 663 (1994)
  - Chapront, J. et al. — A&A 387, 700 (2002), for independent lunar context
  - Park, R.S. et al. — "The JPL Planetary and Lunar Ephemerides DE440 and
    DE441", AJ 161, 105 (2021)
  - Lieske, J.H. — "Galilean satellite ephemerides E5", A&AS 129, 205 (1998)
- **Historical hypothetical-body conventions (IDs 40–58)**. Thirteen reviewed
  source-backed models are built in: the twelve orbital models represented by
  the Neely Hamburg-school points (IDs 40–47), Harrington (50), Le Verrier
  (51), Adams (52), and Lowell (53), plus the independently reconstructed
  uniform seven-year Selena convention (56), unwrapped from its published
  1800–2000 baseline. IDs 48, 49, 54, 55, 57, and 58
  remain named API constants but fail closed because a source-complete model
  was not recovered. Per-body sources and transformations are documented in the
  [source record](docs/methodology/hypothetical-bodies.md); exact missing fields
  are inventoried in the
  [missing-models inventory](docs/methodology/missing-hypothetical-models.md).
- **Sidereal zero points and catalogue geometry** documented per mode in the
  [ayanamsha source table](docs/reference/ayanamsha.md).

## API Compatibility

LibEphemeris provides an API signature-compatible with the reference
ephemeris API (function names, parameters and flag constants match, e.g.
`calc_ut`, `houses`, `julday`) to allow drop-in migration. API compatibility
does not imply code derivation: the algorithms and implementation are
independently written.

## Validation

The current validation policy permits only ephemeral calls to the reference
public API: raw outputs must not be recorded as fixtures, fitted into
coefficients, or encoded into artifacts. Current release artifacts do not ship
reference-distribution source, prose, algorithms, or data. Earlier development
history contained legacy/reference-derived material that the v3 remediation
removed; this notice does not rewrite the provenance or licensing of those
earlier revisions or distributions.

Current numerical models have model-specific source records and integrity
gates. Mean lunar points use ERFA's IERS 2003 Delaunay arguments; true and
osculating points use the active NASA JPL state; and supported historical
hypothetical models use the cited primary transcriptions in their source
record. A passing technical gate is not a legal opinion about historical
releases or third-party rights.

## Vendored Components

Three files are vendored or adapted third-party code and keep their upstream
licenses, all permissive (MIT): `vendor/spktype21.py`,
`moon_theories/tass17.py`, and `moon_theories/tass17_data.py`. The package
contains no vendored copyleft code in the current release artifact. Historical
versions and separately installed dependencies retain their own terms. Full
inventory:
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## AI-Assisted Development

Parts of this codebase were developed with AI assistance, including Anthropic
Claude and OpenAI Codex, under the direction and review of the project author.
Copyright in the project-authored work is held by Giacomo Battaglia.

## Trademarks

"Swiss Ephemeris" is a product name of Astrodienst AG. References to it in
this repository are nominative only (interface compatibility). LibEphemeris is
not affiliated with, or endorsed by, Astrodienst AG.
