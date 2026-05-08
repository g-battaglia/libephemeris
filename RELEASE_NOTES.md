# Release Notes

## 1.2.0 - 2026-05-08

LibEphemeris 1.2.0 adds catalog listing and ordered batch calculation for
fixed-star discovery workflows.

Highlights:

- `list_fixed_stars()` returns the existing fixed-star catalog entries.
- `swe_batch_fixstars_ut()` and `batch_fixstars_ut()` calculate multiple stars
  in input order.
- `skip_errors=True` keeps unresolved stars as `None` in their original slots.
- Batch fixed-star calculations reuse time and observer positions locally within
  each call without adding persistent per-star caches.

This release is backward-compatible and additive.
