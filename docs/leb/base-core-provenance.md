# Bundled LEB2 base-core provenance

This is the build attestation for the only tracked prebuilt LEB artifact in the
repository.

## Artifact

- Path: `libephemeris/data/leb2/base_core.leb2`
- Build date: 2026-07-13
- Size: 10,841,639 bytes
- SHA-256: `a02b15344de946f8d0945c30c3ad47c3c1ce69f335af99e545c777fe9ec1bcfd`
- Body IDs: 0–12 and 14
- Format: LEB2 v2; populated section types 0, 2, 3, 4, and 6

The artifact was generated with:

```bash
uv run python scripts/generate_leb2.py generate --tier base --group core --workers 8 --output /tmp/libephemeris-base-core-clean.leb2
```

The reviewed result was then installed at
`libephemeris/data/leb2/base_core.leb2`.

## Scientific inputs

- Planetary, Earth, and lunar state data: NASA JPL DE440s.
- Mean-node and mean-apogee inputs: ERFA's implementation of the IERS 2003
  Delaunay arguments plus conventional orbital-plane geometry.
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

No Swiss Ephemeris source, documentation, algorithm, data file, or generated
artifact was used. The build made no call to the reference API and consumed no
external comparison output. It did not convert or reuse any previously published
LEB file.
The retired hypothetical group was neither generated nor converted; active
tooling excludes it from publication.

This page is an external build record. The command, hash, input description,
and clean-room attestation are not embedded as metadata inside the LEB2 format.
