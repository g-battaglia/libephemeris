# NOTICE

LibEphemeris is an astronomical ephemeris library for Python with a Swiss
Ephemeris-compatible API, licensed under the **GNU Affero General Public
License v3.0** (see [LICENSE](LICENSE) and [LICENSING.md](LICENSING.md)).

## Data and Model Sources

All astronomical computations are based on:

- **JPL DE440/DE441** planetary and lunar ephemerides, accessed through the
  local [Skyfield](https://rhodesmill.org/skyfield/) path (Brandon Rhodes, MIT
  license), precomputed LEB coefficients, or NASA JPL Horizons state vectors
  ([bundled base-core build attestation](docs/leb/base-core-provenance.md))
- **IAU 2006/2000A** precession-nutation and **Vondrak 2011** long-term
  precession, via [pyerfa](https://github.com/liberfa/pyerfa) (BSD-3-Clause)
- **Peer-reviewed academic sources**, including:
  - Meeus, J. -- *Astronomical Algorithms*, 2nd ed. (1998)
  - Simon, J.L. et al. -- A&A 282, 663 (1994)
  - Chapront, J. et al. -- A&A 387, 700 (2002), for independent lunar context
  - Park, R.S. et al. -- "The JPL Planetary and Lunar Ephemerides DE440 and
    DE441", AJ 161, 105 (2021)
  - Lieske, J.H. -- "Galilean satellite ephemerides E5", A&AS 129, 205 (1998)
- **Historical hypothetical-body conventions (IDs 40-58)**. Per-body source
  status is documented in the
  [per-body source record](docs/methodology/hypothetical-bodies.md).
- **Sidereal zero points and catalogue geometry** documented per mode in the
  [ayanamsha source table](docs/reference/ayanamsha.md).

## API Compatibility

LibEphemeris provides an API signature-compatible with the reference
ephemeris API (function names, parameters and flag constants match, e.g.
`calc_ut`, `houses`, `julday`) to allow drop-in migration. API compatibility
does not imply code derivation.

## Validation

Compatibility with the reference API is verified by calling its public
interface and comparing outputs. The development process aims to keep these
comparisons ephemeral. Reference-distribution source, prose, algorithms, and
data files should not be part of the shipped package.

All numerical data and models are reproducible from the sources
listed above. Mean lunar points use ERFA's IERS 2003 Delaunay arguments; true
and osculating points use the active NASA JPL state; and Harrington uses its
cited publication.

## Vendored Components

Three files are vendored or adapted third-party code and keep their upstream
licenses, all permissive (MIT): `vendor/spktype21.py`,
`moon_theories/tass17.py`, and `moon_theories/tass17_data.py`. The full
inventory: [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## AI-Assisted Development

Parts of this codebase were developed with AI assistance (Anthropic Claude)
under the direction and review of the project author. Copyright in the work
is held by Giacomo Battaglia.

## Trademarks

"Swiss Ephemeris" is a product name of Astrodienst AG. References to it in
this repository are nominative only (interface compatibility). LibEphemeris is
not affiliated with, or endorsed by, Astrodienst AG.
