# Bundled base-tier LEB2 provenance

This is the build attestation for the two tracked prebuilt LEB artifacts in the
repository: the base core and the small Hamburg-body companion.

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
  Delaunay arguments plus conventional orbital-plane geometry. This build was
  performed after removal of the superseded compatibility baseline and does
  not use an earlier LEB body channel as an input or target.
- Serialization and compression: the project-owned LEB2 generator and format
  implementation in this repository.

## Auxiliary sections

The attestation covers every payload in the file, not only the body channels:

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

## Hamburg-body companion

- Path: `libephemeris/data/leb2/base_uranians.leb2`
- Build date: 2026-07-18
- Size: 46,148 bytes
- SHA-256: `2b052321672f995a2c2f6c6c3abe4dd623e2721eb17f7319873d0c8b0d65ee8c`
- Body IDs: 40–47
- Coverage: JD 2396758.5–2506331.5 (1850–2150)
- Format: LEB2 v2; populated section types 0 and 6

The standalone `ephemeris_base_uranians.leb` partial is generated from the
independently sourced Neely (1980) propagation in
`libephemeris.hypothetical`, verified over 500 samples per body, converted with
the `1e-12` native-component target, and verified over another 200 samples per
body. Shared nutation, Delta-T, and star-catalog tables are intentionally owned
by the core and are not duplicated in this or any other named companion.

The recorded build inputs contain no reference-distribution source,
documentation, algorithm, data file, generated artifact, or comparison output;
the build regenerated the asset rather than converting a previously published
LEB file. This attestation is scoped to the named build inputs and artifact. It
does not make a claim about every historical repository revision.
The hypothetical group, retired at the time of this build, was neither
generated into nor converted as part of this `base_core.leb2` artifact. The
Hamburg bodies (40–47) now ship separately as the independently sourced,
manifest-pinned `base_uranians.leb2` companion.

This page is an external build record. The command, hashes, input description,
and clean-room attestation are not embedded as metadata inside the LEB2 format.
