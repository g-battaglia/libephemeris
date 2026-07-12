# NOTICE

LibEphemeris is an independent implementation of an astronomical ephemeris
library for Python, licensed under the **Apache License, Version 2.0** (see
[LICENSE](LICENSE) and [LICENSING.md](LICENSING.md)).

## Data and Model Sources

All astronomical computations are based on:

- **JPL DE440/DE441** planetary and lunar ephemerides, accessed through the
  local [Skyfield](https://rhodesmill.org/skyfield/) path (Brandon Rhodes, MIT
  license), precomputed LEB coefficients, or NASA JPL Horizons state vectors
- **IAU 2006/2000A** precession-nutation and **Vondrák 2011** long-term
  precession, via [pyerfa](https://github.com/liberfa/pyerfa) (BSD-3-Clause)
- **Peer-reviewed academic sources**, including:
  - Meeus, J. — *Astronomical Algorithms*, 2nd ed. (1998)
  - Simon, J.L. et al. — A&A 282, 663 (1994)
  - Chapront, J. et al. — A&A 387, 700 (2002); ELP 2000-85 lunar theory
  - Park, R.S. et al. — "The JPL Planetary and Lunar Ephemerides DE440 and
    DE441", AJ 161, 105 (2021)
  - Lieske, J.H. — "Galilean satellite ephemerides E5", A&AS 129, 205 (1998)
- **Primary historical sources** for hypothetical bodies: Witte, A. &
  Lefeldt, H. — *Regelwerk für Planetenbilder* (1928); Neely, J. (1988);
  Makransky, B. — *Primary Directions* (1988).

## API Compatibility

LibEphemeris provides an API signature-compatible with the reference
ephemeris API (function names, parameters and flag constants match, e.g.
`calc_ut`, `houses`, `julday`) to allow drop-in migration. API compatibility
does not imply code derivation: the algorithms and implementation are
independently written.

## Calibration Data Disclosure

Two generated data sets are fitted against the reference API used strictly
as a black-box oracle, and are disclosed for transparency: the interpolated
lunar-apse residual tables (`libephemeris/lunar_apse_corrections.py` and the
perturbation coefficients in `lunar.py` — the INTP_APOG/INTP_PERG bodies are
constructs defined by the reference API, so behavioral parity requires
fitting to its output), and four fictitious-body element rows in
`libephemeris/data/fictitious_orbits.csv` (Nibiru, Proserpina, the Selena
digits, the Waldemath reconstruction) that have no known publication and are
carried as interoperability values. These tables contain numeric program
*output*, not source expression; all runtime code is original.

## Vendored Components

Three files are vendored or adapted third-party code and keep their upstream
licenses, all permissive (MIT): `vendor/spktype21.py`,
`moon_theories/tass17.py`, and `moon_theories/tass17_data.py`. The package
contains no copyleft code. Full inventory:
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## AI-Assisted Development

Parts of this codebase were developed with AI assistance (Anthropic Claude)
under the direction and review of the project author. Copyright in the work
is held by Giacomo Battaglia.

## Trademarks

"Swiss Ephemeris" is a product name of Astrodienst AG. References to it in
this repository are nominative only (interface compatibility and black-box
verification targets). LibEphemeris is not affiliated with, or endorsed by,
Astrodienst AG.
