# Bundled base-tier LEB2 provenance

This is the build record for the tracked prebuilt LEB artifact in the
repository: the base core.

## Core artifact

- Path: `libephemeris/data/leb2/base_core.leb2`
- Build date: 2026-07-18
- Size: 10,232,283 bytes
- SHA-256: `5d708bdbe3e799e0802ba575984e57a3c5e44720dbfa1b4a01cf826640e0cb82`
- Body IDs: 0–12 and 14
- Format: LEB2 v2; populated section types 0, 2, 3, 4, and 6

The complete data-v3 matrix, including this artifact, was generated with:

```bash
./regenerate-leb.sh all -q
```

For the base core, that workflow generated and merged the LEB1 body groups,
verified 500 samples per body, converted the `core` group to LEB2, and verified
another 200 samples per body. The reviewed result was then installed at
`libephemeris/data/leb2/base_core.leb2` byte-for-byte.

Exact generation-input pins:

- DE440s SPK SHA-256:
  `c1c7feeab882263fc493a9d5a5b2ddd71b54826cdf65d8d17a76126b260a49f2`;
- `scripts/generate_leb.py` SHA-256:
  `85c8574c9acaea21b91d5dabb856800c29b717a36f1f592eab481da6c45458d5`;
- `scripts/generate_leb2.py` SHA-256:
  `8db3ff23e1f93d0b8b0af4c5937a2da8a8eb777638164af42b0ab9c210d3c4ea`;
- `libephemeris/mean_lunar_apse.py` SHA-256:
  `23285afe9817498b9bdd0f9daa1ea2403c4beb676068af2f20250620cdb55a9b`;
- `libephemeris/lunar.py` SHA-256:
  `9bb440497c9cf25f3eb43a3a425f84b940662bbf9bc4f14eb9a28c10529f83dc`.

These pins identify the bytes used for this artifact even if later commits
change a generator or lunar model. The source kernel remains an external JPL
input and is not relicensed by this project.

## Scientific inputs

- Planetary, Earth, and lunar state data: NASA JPL DE440s.
- Mean-node and mean-apogee inputs: ERFA's implementation of the IERS 2003
  Delaunay arguments plus conventional orbital-plane geometry, evaluated
  directly; no earlier LEB body channel is used as an input or target.
- Serialization and compression: the project-owned LEB2 generator and format
  implementation in this repository.

## Auxiliary sections

This record covers every payload in the file, not only the body channels:

- **Nutation:** `scripts.generate_leb.generate_nutation()` evaluates ERFA
  `nut06a`, the IAU 2006/2000A model, into 6,849 Chebyshev segments with a
  16-day interval, degree 16, and two components. The serialized section is
  1,862,968 bytes.
- **Delta T:** `scripts.generate_leb.generate_delta_t()` samples the library's
  default `smh2016`/IERS time model every 30 days. The section contains 3,654
  samples over JD 2396758.5–2506331.5 and is 58,472 bytes. A fresh generation
  with no user override and IERS live updates disabled reproduced every stored
  value bit for bit.
- **Fixed stars:** `scripts.generate_leb.generate_star_catalog()` serializes
  1,447 entries from `libephemeris.star_catalog_gen` and
  `libephemeris.fixed_stars`. Those modules build the catalog independently
  from the Hipparcos new reduction (VizieR I/311), ESA Hipparcos (I/239),
  Kostjuk's cross-index (IV/27A), XHIP (V/137D), and IAU WGSN names. The LEB
  records contain RA, declination, proper motion, and magnitude; parallax and
  radial-velocity fields are deliberately zero. The section is 92,608 bytes.
- **Indexes and compressed body data:** the remaining populated sections are
  generated entirely by the repository's LEB2 serializer from the body and
  auxiliary values described above.

## Hamburg-body companion (retired)

The formerly bundled `base_uranians.leb2` companion (Hamburg bodies 40–47)
was retired in 3.1.0; those bodies are always computed from the runtime
Neely (1980) propagation in `libephemeris.hypothetical`.

The core asset was regenerated from the inputs listed above rather than
converted from a previously published LEB file. The core artifact contains no
fictitious body IDs.

This page is an external build record. The command, hashes, and input
description are not embedded as metadata inside the LEB2 format.
