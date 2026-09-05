# Third-Party Notices

This file inventories third-party software and data that LibEphemeris
**contains** (vendored or adapted code shipped inside the package) or
**depends on** (PyPI dependencies installed alongside it), together with
their licenses. It is maintained as part of the AGPL-3.0 compliance
process -- see [LICENSING.md](LICENSING.md) and [NOTICE.md](NOTICE.md).

## Code shipped inside the package (vendored / adapted)

| Component | Path | Upstream | License |
|---|---|---|---|
| spktype21 0.1.0 | `libephemeris/vendor/spktype21.py` | Shushi Uetsuki (whiskie14142) | MIT |
| TASS 1.7 port | `libephemeris/moon_theories/tass17.py`, `libephemeris/moon_theories/tass17_data.py` (generated periodic-term tables) | Johannes Gajdosik's Stellarium implementation of Vienne & Duriez TASS 1.7 | MIT |

Each of these files carries its own `SPDX-License-Identifier` header and
documents its upstream provenance in the module docstring. They are **not**
covered by the project's AGPL-3.0 license, but they are all permissively
(MIT) licensed. The two entries below record, for each component, the
retained upstream copyright notice, the license file that was consulted
with its SHA-256 (re-check a pin with `curl -sL <url> | shasum -a 256`),
the local modifications, and any third-party material the files carry
beyond the upstream author's own work.

### spktype21 0.1.0 (`libephemeris/vendor/spktype21.py`)

- **Upstream:** `spktype21` 0.1.0 on PyPI, by Shushi Uetsuki (GitHub
  `whiskie14142`), repository <https://github.com/whiskie14142/spktype21>.
  The PyPI metadata declares `License :: OSI Approved :: MIT License`.
- **License:** MIT. Retained copyright notice, as published upstream:
  Copyright (c) 2018 The Python Packaging Authority.
- **License file consulted (2026-09-05):**
  <https://raw.githubusercontent.com/whiskie14142/spktype21/9b1ef79aca2e06988bb7d5a2bb10c42066fc5186/LICENSE.txt>
  (commit `9b1ef79a` of 2018-10-15, the last change to that file; identical
  to `master` on the date consulted).
  SHA-256: `ec423cc5506eea1ffbfc9955c3ec8f86139996963d84ff306a5ee41eda8a4ff1`.
- **Local modifications:** two `.item()` calls added to the
  `daf.map_array()` results for numpy 2.x compatibility, and the vendoring
  record in the module header. The numerical implementation is otherwise the
  identified 0.1.0 upstream.
- **NASA/JPL/NAIF material carried by the file.** The upstream author states
  that the module "has been developed based on jplephem and FORTRAN source
  of the SPICE Toolkit of NASA/JPL/NAIF", and the file reproduces comments
  from that toolkit. Inventory (line numbers as of this revision):
  - Module docstring, lines 20, 26, 64 and 66: the statement that the input
    format is NASA/JPL/NAIF SPK type 21, the link to the NAIF "SPK Required
    Reading" (`naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/spk.html`),
    the sentence quoted above, and the SPICE Toolkit link
    (`naif.jpl.nasa.gov/naif/toolkit.html`).
  - Line 78: `MAXTRM = 25`, noted as included from `spk21.inc` of the
    FORTRAN source `spke21.f`.
  - Lines 256-266: the author's note that `spke21()` was translated from the
    SPICE Toolkit FORTRAN routine `spke21.f` and then modified, with the
    links at lines 259-260 (SPICE Toolkit for FORTRAN; SPK Required Reading,
    C edition).
  - Lines 269-421 (`#C` prefix): the SPICELIB header of routine `SPKE21`
    reproduced verbatim -- `$ Abstract`; `$ Disclaimer` (lines 274-297: the
    Caltech/JPL/NASA notice that the software was created under a U.S.
    Government contract, is publicly available under U.S. export laws, and
    is provided "as-is" without warranty, the recipient bearing all risk and
    indemnifying Caltech and NASA); `$ Required_Reading`; `$ Keywords`;
    `$ Declarations`; `$ Brief_I/O`; `$ Detailed_Input`; `$ Detailed_Output`;
    `$ Parameters`; `$ Exceptions`; `$ Files`; `$ Particulars`; `$ Examples`;
    `$ Restrictions`; `$ Literature_References` (line 401: NAIF Document
    168.0, "S- and P- Kernel (SPK) Specification and User's Guide");
    `$ Author_and_Institution` (N.J. Bachman, F.T. Krogh, W.L. Taber,
    I.M. Underwood -- JPL); `$ Version` (SPICELIB Version 1.0.0,
    03-FEB-2014); `$ Index_Entries`.
  - Lines 423-575: a further 67 `#C` comment lines inside the body of
    `spke21()`, the FORTRAN routine's own step-by-step comments carried
    along with the translation.

  The SPICE Toolkit is implemented and maintained by Caltech/JPL under
  contract to NASA and is distributed by NAIF under the rules published at
  <https://naif.jpl.nasa.gov/naif/rules.html>, which state that including
  Toolkit modules "as part of a package supporting a customer-built
  SPICE-based tool is entirely appropriate", ask that modifications be noted
  in the header, and designate the software "Technology and Software
  Publicly Available" for export purposes. The disclaimer is retained
  verbatim in the file, and the module header records both the translation
  and the local modifications.

### TASS 1.7 port (`libephemeris/moon_theories/tass17.py`, `libephemeris/moon_theories/tass17_data.py`)

- **Upstream:** the Stellarium implementation of the TASS 1.7 theory of the
  Saturnian satellites, file `src/core/planetsephems/tass17.c`, by Johannes
  Gajdosik. `tass17.py` is the Python port of the evaluation code;
  `tass17_data.py` is the periodic-term table generated from the same C
  file.
- **License:** MIT, granted in the header of `tass17.c` itself (the file
  carries its own permission notice; there is no separate license file for
  it). Retained copyright notice: Copyright (c) 2005 Johannes Gajdosik.
  The header states that this notice covers the author's compilation of the
  theory into software, not the underlying work of Alain Vienne and Luc
  Duriez.
- **License file consulted (2026-09-05):**
  <https://raw.githubusercontent.com/Stellarium/stellarium/eb189190cee37fdc14939c2407e4e3fd47e8b6a0/src/core/planetsephems/tass17.c>
  (commit `eb189190` of 2026-06-20, the last change to that file; identical
  to `master` on the date consulted). The MIT notice occupies lines 1-35.
  SHA-256 of the whole file:
  `bb001d2e0c2b2376c246264e3845bc59b4984b76dcde692a39a1c912bc70cbcc`.
- **Underlying theory:** Vienne, A. & Duriez, L. (1995), A&A 297, 588-605
  (TASS 1.6); TASS 1.7 is the later IMCCE revision distributed as Fortran
  code and data at `ftp://ftp.imcce.fr/pub/ephem/satel/tass17/`. The
  Stellarium header notes that its author "can neither allow nor forbid the
  usage of the TASS 1.7 theory"; the theory itself is the published
  scientific work of its authors.

The Galilean satellite module (`libephemeris/moon_theories/galilean.py`)
was, through v2.1.0, adapted from PyMeeus and licensed LGPL-3.0. <!-- provenance-implementation-ok --> It was
rewritten in June 2026 as an independent implementation of the published
theory (Lieske, J.H. 1998, "Galilean satellite ephemerides E5", A&AS 129,
205; Meeus 1998, *Astronomical Algorithms*, ch. 44) and is now owned by
the project under the AGPL-3.0 license. Theory specification:
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

## Optional `nbody` extra -- GPL-3.0 (not bundled)

The opt-in `libephemeris[nbody]` extra (deliberately NOT part of
`libephemeris[all]`) installs two **GPL-3.0-or-later** packages used by the
shipped `rebound_integration.py` module for ephemeris-quality n-body
propagation of minor bodies:

| Package | License | Notes |
|---|---|---|
| rebound | GPL-3.0-or-later | N-body integrator; imported lazily, only when the extra is installed |
| assist | GPL-3.0-or-later | Ephemeris-quality REBOUND extension (JPL small-body integrator) |

These packages are **not** part of the core install, are **never bundled or
redistributed** in any LibEphemeris artifact (they are installed separately
from PyPI under their own GPL license), and are imported only on demand.
GPL-3.0 and AGPL-3.0 are compatible per their respective Section 13
provisions; a user who installs `libephemeris[nbody]` forms a combined
work governed by the AGPL-3.0.

## Data sources

- **JPL DE440 / DE441** planetary and lunar ephemerides (NASA/JPL/Caltech):
  publicly distributed scientific data; downloaded at runtime, not bundled.
- **ESA Hipparcos Catalogue** (ESA SP-1200, 1997) and **van Leeuwen 2007**
  reduction: star positions and proper motions (facts of nature) compiled
  into `libephemeris/fixed_stars.py`, with attribution in the module.
- **IAU Working Group on Star Names** (WGSN) designations.
- **IERS** Earth-orientation / Delta-T / leap-second tables: public data,
  fetched at runtime.
- `libephemeris/data/fictitious_orbits.csv`: 12 transcribed historical
  element rows. IDs 40-47 come from James Neely, "Orbital Elements for the
  Transneptunian Planets," *Matrix Magazine* VII (1980), Table I, p. 8;
  ID 50 comes from Harrington (1988), *AJ* 96, p. 1478; ID 51 from
  Le Verrier (1846), *Comptes rendus* 23, p. 432 and his longer *Recherches*;
  ID 52 from Adams (1846), sections 47-48, cross-checked against Gould's 1850
  Smithsonian report; and ID 53 from Lowell's 1915 memoir, pp. 5, 9, and 105.
  Project-authored conversions are documented in
  [the per-body source record](docs/methodology/hypothetical-bodies.md).
  Source publications retain whatever rights apply to them.
- **White Moon / Selena (ID 56):** the uniform seven-year rule and its
  checkpoints are factual inputs reviewed at pp. 17, 18, 20, 29, and 45 of
  Felix Velichko and Max Larin, *Lilith, Selena, Proserpina: Articles and
  Ephemerides 1800-2050* (Mir Uranii, 2007), ISBN 5-900191-12-5. No scan or
  ephemeris table is bundled. LibEphemeris supplies its own explicit first-day
  TT convention, unwraps uniform motion from the publication's 200-year
  January baseline, and supplies an IAU-standards-derived display distance;
  full arithmetic is in the
  [per-body source record](docs/methodology/hypothetical-bodies.md).
  IDs 48, 49, 54, 55, 57, and 58 have no bundled numerical model and fail
  closed; the exact unrecovered fields are listed in the
  [missing-models inventory](docs/methodology/missing-hypothetical-models.md).
- `libephemeris/data/leb2/base_core.leb2`: generated by this project from JPL
  DE440 state data and ERFA/IERS standards-derived abstract lunar channels;
  its source and build attestation are documented in
  [docs/leb/base-core-provenance.md](docs/leb/base-core-provenance.md).

## Development-only tools (not shipped, never linked)

- **pyswisseph** (AGPL-3.0 / Swiss Ephemeris dual license): may be used
  for ephemeral compatibility comparisons in a separate validation environment.
  It is not declared among this package's dependencies (including the `dev`
  extra), is never imported by the shipped package, and is not a runtime
  dependency. See [NOTICE.md](NOTICE.md) for details.
