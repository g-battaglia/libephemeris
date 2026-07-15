# Third-Party Notices

This file inventories third-party software and data that LibEphemeris
**contains** (vendored or adapted code shipped inside the package) or
**depends on** (PyPI dependencies installed alongside it), together with
their licenses. It is maintained as part of the Apache-2.0 compliance
process — see [LICENSING.md](LICENSING.md) and [NOTICE.md](NOTICE.md).

## Code shipped inside the package (vendored / adapted)

| Component | Path | Upstream | License |
|---|---|---|---|
| spktype21 0.1.0 | `libephemeris/vendor/spktype21.py` | Shushi Uetsuki (whiskie14142) | MIT |
| TASS 1.7 port | `libephemeris/moon_theories/tass17.py`, `libephemeris/moon_theories/tass17_data.py` (generated periodic-term tables) | Johannes Gajdosik's Stellarium implementation of Vienne & Duriez TASS 1.7 | MIT |

Retained upstream copyright notices:

- spktype21 license: Copyright (c) 2018 The Python Packaging Authority
  ([upstream license](https://github.com/whiskie14142/spktype21/blob/master/LICENSE.txt)).
- Stellarium TASS 1.7 implementation: Copyright (c) 2005 Johannes Gajdosik
  ([upstream file header](https://github.com/Stellarium/stellarium/blob/master/src/core/planetsephems/tass17.c)).

Each of these files carries its own `SPDX-License-Identifier` header and
documents its upstream provenance in the module docstring. They are **not**
covered by the project's Apache-2.0 license, but they are all permissively
(MIT) licensed. The current release wheel contains no vendored copyleft code;
installed dependencies and optional extras remain governed by their own terms.

The Galilean satellite module (`libephemeris/moon_theories/galilean.py`)
was, through v2.1.0, adapted from PyMeeus and licensed LGPL-3.0. <!-- provenance-implementation-ok --> The
current module was rewritten in June 2026 from the published theory (Lieske,
J.H. 1998, "Galilean satellite ephemerides E5", A&AS 129, 205; Meeus 1998,
*Astronomical Algorithms*, ch. 44) and is project-authored Apache-2.0 code.
That statement describes the current file, not the licensing of earlier
distributions. Theory specification:
[docs/methodology/galilean-e5-spec.md](docs/methodology/galilean-e5-spec.md).

### MIT license text (spktype21, TASS 1.7 port)

> Permission is hereby granted, free of charge, to any person obtaining a
> copy of this software and associated documentation files (the
> "Software"), to deal in the Software without restriction, including
> without limitation the rights to use, copy, modify, merge, publish,
> distribute, sublicense, and/or sell copies of the Software, and to
> permit persons to whom the Software is furnished to do so, subject to
> the following conditions:
>
> The above copyright notice and this permission notice shall be included
> in all copies or substantial portions of the Software.
>
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
> OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
> MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
> IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
> CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
> TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
> SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Runtime dependencies (installed from PyPI, not redistributed here)

| Package | License | Notes |
|---|---|---|
| skyfield | MIT | Ephemeris access, time scales |
| skyfield-data | MIT | Bundled astronomy data files |
| jplephem | MIT | Transitive via skyfield; SPK kernel reader |
| sgp4 | MIT | Transitive via skyfield |
| numpy | BSD-3-Clause | Transitive |
| pyerfa | BSD-3-Clause | Bundles the ERFA C library, distributed under the ERFA license (permissive, SOFA-derived, with naming restrictions) |
| click | BSD-3-Clause | CLI framework |
| certifi | MPL-2.0 | Used unmodified; MPL file-level copyleft imposes no obligation on LibEphemeris code |
| zstandard | BSD-3-Clause | LEB2 compression |

Dependencies are installed from PyPI by the user or installer under their
own licenses; they are not part of LibEphemeris's license grant.

## Optional permissive tooling

| Package | License | Notes |
|---|---|---|
| astropy | BSD-3-Clause | Installed by the `stars` extra for star-catalog tooling; also a development dependency |
| astroquery | BSD-3-Clause | Installed by the `stars` extra and development environment for catalog/validation scripts; runtime SPK downloads use LibEphemeris's direct JPL Horizons HTTPS client |

## Optional `nbody` extra — GPL-3.0 (not bundled)

The opt-in `libephemeris[nbody]` extra (deliberately NOT part of
`libephemeris[all]`, which stays permissive-only) installs two
**GPL-3.0-or-later** packages used by the shipped `rebound_integration.py`
module for ephemeris-quality n-body propagation of minor bodies:

| Package | License | Notes |
|---|---|---|
| rebound | GPL-3.0-or-later | N-body integrator; imported lazily, only when the extra is installed |
| assist | GPL-3.0-or-later | Ephemeris-quality REBOUND extension (JPL small-body integrator) |

These packages are not part of the core install or bundled in the project
wheel; they are installed separately from PyPI and imported only on demand.
The effect of installing, combining, or redistributing them depends on the
exact versions, packaging and applicable law. Users and redistributors must
comply with those packages' GPL terms; this notice does not characterize every
runtime arrangement as a single combined work. The default/core required
dependency set contains no strong-copyleft dependency. `certifi` is MPL-2.0
and is installed unmodified as a separate dependency.

## Data sources

- **JPL DE440 / DE441** planetary and lunar ephemerides
  (Caltech/JPL under contract with NASA): SPK kernels are downloaded at
  runtime rather than included in the wheel. They must not be described
  categorically as U.S.-government public-domain works. The official
  [NAIF SPICE rules](https://naif.jpl.nasa.gov/naif/rules.html) permit anyone
  to download and use NAIF-hosted kernels and permit redistribution of
  unmodified NAIF kernels, subject to the stated attribution, modification,
  disclaimer and source-specific rules. The generic-kernel directory lists
  [DE440/DE441 artifacts](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/).
  For each downloaded or redistributed asset, retain its source URL, version,
  hash, embedded metadata and then-current provider terms.
- **ESA Hipparcos Catalogue** (ESA SP-1200, 1997) and **van Leeuwen 2007**
  reduction: star positions and proper motions (facts of nature) compiled
  into `libephemeris/fixed_stars.py`, with attribution in the module.
- **IAU Working Group on Star Names** (WGSN) designations.
- **IERS** Earth-orientation / Delta-T / leap-second tables: public data,
  fetched at runtime.
- `libephemeris/data/fictitious_orbits.csv`: 12 independently transcribed
  historical element rows. IDs 40–47 come from James Neely, “Orbital Elements
  for the Transneptunian Planets,” *Matrix Magazine* VII (1980), Table I,
  p. 10; ID 50 comes from Harrington (1988), *AJ* 96, p. 1478; ID 51 from
  Le Verrier (1846), *Comptes rendus* 23, p. 432 and his longer *Recherches*;
  ID 52 from Adams (1846), sections 47–48, cross-checked against Gould's 1850
  Smithsonian report; and ID 53 from Lowell's 1915 memoir, pp. 5, 9, and 105.
  The numerical facts were transcribed afresh and all project-authored
  conversions are documented in
  [the per-body source record](docs/methodology/hypothetical-bodies.md).
  Source publications retain whatever rights apply to them; the Apache-2.0
  grant covers LibEphemeris's project-authored code and presentation to the
  extent copyrightable.
- **White Moon / Selena (ID 56):** the uniform seven-year rule and five
  one-minute checkpoints are factual inputs reviewed at pp. 17, 18, 20, 29,
  and 45 of
  Felix Velichko and Max Larin, *Lilith, Selena, Proserpina: Articles and
  Ephemerides 1800–2050* (Mir Uranii, 2007), ISBN 5-900191-12-5. No scan or
  ephemeris table is bundled. LibEphemeris supplies its own explicit first-day
  TT convention, unwraps uniform motion from the publication's 200-year
  January baseline, and supplies an IAU-standards-derived display distance;
  full arithmetic is in the
  [per-body source record](docs/methodology/hypothetical-bodies.md). The source
  publication retains its rights; this notice does not purport to relicense it.
  IDs 48, 49, 54, 55, 57, and 58 have no bundled numerical model and fail
  closed; the exact unrecovered fields are listed in the
  [missing-models inventory](docs/methodology/missing-hypothetical-models.md).
- `libephemeris/data/leb2/base_core.leb2`: generated by this project from JPL
  DE440 state data and ERFA/IERS standards-derived abstract lunar channels.
  It is a project-generated representation, not an assertion that the
  underlying source data are relicensed under Apache-2.0; its source and build
  attestation are documented separately.

## Development-only tools (not shipped, never linked)

- **pyswisseph** (AGPL-3.0 / Swiss Ephemeris dual license): may be used
  for ephemeral compatibility comparisons in a separate validation environment.
  It is not declared among this package's dependencies (including the `dev`
  extra), is never imported by the shipped package, and is not a runtime
  dependency. See [NOTICE.md](NOTICE.md) for the independence statement.
- **Reference-distribution source and data artifacts are prohibited.** Neither
  this repository nor its separate validation worktree may contain or consume
  those files. Compatibility validation may call the installed Python API and
  compare outputs, but results must remain ephemeral. The provenance gate
  applies name checks to the physical project tree (including ignored and
  untracked paths) and reachable Git history, and scans project text outside
  its exact infrastructure-cache exclusions.
