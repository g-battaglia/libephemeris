# Bundled LEB2 base-core provenance

This is the build attestation for the only tracked prebuilt LEB artifact in the
repository.

## Artifact

- Path: `libephemeris/data/leb2/base_core.leb2`
- Build date: 2026-07-15
- Size: 10,841,639 bytes
- SHA-256: `e5a9730b09f4a21dd35c7adcc938767644eb242145807bd8a6bf7e6042f5b420`
- Body IDs: 0–12 and 14
- Format: LEB2 v2; populated section types 0, 2, 3, 4, and 6

The artifact was generated with:

```bash
uv run python scripts/generate_leb2.py generate --tier base --group core --workers 8 --output /tmp/libephemeris-base-core-clean.leb2
```

The reviewed result was then installed at
`libephemeris/data/leb2/base_core.leb2`.

Exact generation-input pins:

- DE440s SPK SHA-256:
  `c1c7feeab882263fc493a9d5a5b2ddd71b54826cdf65d8d17a76126b260a49f2`;
- `scripts/generate_leb.py` SHA-256:
  `a6b7376393f1431711d7fc63480bda580519aea94b4805fcc753c0484f05a322`;
- `scripts/generate_leb2.py` SHA-256:
  `8bcbcd70d3972b78bc4fce00980531c897c446152c8f8d38a108b24db7af518d`;
- `libephemeris/mean_lunar_apse.py` SHA-256:
  `ee2cfe81464e950da9760fed1942b4003178f7069fd75c058aa38a882d06b068`;
- `libephemeris/lunar.py` SHA-256:
  `83e48a4cf6c9eb18c088e91f4862bb781e68a53d3a942e005f2612919d5199af`.

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

The recorded build inputs contain no reference-distribution source,
documentation, algorithm, data file, generated artifact, or comparison output;
the build regenerated the asset rather than converting a previously published
LEB file. This attestation is scoped to the named build inputs and artifact. It
does not make a claim about every historical repository revision.
The hypothetical group, retired at the time of this build, was neither
generated into nor converted as part of this `base_core.leb2` artifact. The
Hamburg bodies (40–47) now ship separately as the independently sourced,
manifest-pinned `base_uranians.leb2` companion.

This page is an external build record. The command, hash, input description,
and clean-room attestation are not embedded as metadata inside the LEB2 format.
