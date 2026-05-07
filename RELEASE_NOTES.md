# Release Notes

## 1.2.0 — 2026-05-08

LibEphemeris 1.2.0 adds native fixed-star discovery support for downstream
consumers that need backend-specific catalogs without depending on Swiss
Ephemeris data files.

Highlights:

- `list_fixed_stars()` exposes libephemeris-owned fixed-star metadata.
- `swe_batch_fixstars_ut()` / `batch_fixstars_ut()` calculate multiple stars in
  one ordered batch and can preserve unresolved slots with `skip_errors=True`.
- Fixed-star calculations now reuse cached observer-at-time positions, reducing
  repeated Earth-position work in chart-level fixed-star scans.

This release is backward-compatible and additive.
