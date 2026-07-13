# TODO - Improvements and Technical Debt

## Deferred from the v3.0.0 release review (2026-06-25)

These are the **minor** findings from the v3.0.0 release-readiness audit. They
were deliberately deferred: the release fixes only the serious / release-gating
items (the `houses_ex2`/`houses_armc_ex2` cusp-speed `hsys`-normalization bug,
the clean-room provenance sweep, the public-API native-float leaks in
`heliacal.py`, the `__license__` SPDX expression, and the crashing manual
examples). Everything below is low-impact and can be batched later.

### Correctness / robustness (low impact, edge cases)

- `minor_bodies.py:~2733` — `asin(z / r)` is unclamped; a float round-off
  `|z/r| > 1` would raise `ValueError`. Clamp the argument to `[-1, 1]`.
- `time_utils.py:~760-789` — `jdut1_to_utc` round-trip yields `seconds = 59.99998`
  (~14 µs) for a clean input; tighten the rounding.
- `planetary_moons.py:~262-301` — `register_moon_spk` leaks the loaded SPK kernel
  when it finds no valid moons (error path only); close it before returning.
- `planets.py:~2406-2485` — `FLG_TOPOCTR` is silently ignored for
  Uranian/Transpluto bodies (returns geocentric); decide raise-vs-fallback and
  document it alongside the other per-body-class `FLG_TOPOCTR` policy.
- `state.py:~1913-1922` — `get_current_file_data()` reports `denum = 0` for the
  base tier (`de440s.bsp`); surface the real DE number.

### Performance (non-blocking)

- `spk.py:~1435-1466` — SPK type and coverage are re-read from disk on every
  position evaluation; memoize per (file, target).

### Hygiene / style

- Add `from __future__ import annotations` to the two `planetary_moons` data
  helpers that omit it.
- `astropy_integration.py:~222` — `parse_time_string` shadows the builtin
  `format`; rename the local.
- `astropy_integration.py:~364-380` — `compare_coordinate_transforms` /
  `compare_time_conversions` leak numpy floats in their returned dict (diagnostic
  helpers, not core API).
- Docstring corrections: `state.py` `get_tid_acc()` example claims `-25.936`
  (actual default `-25.8`); `eclipse.py:~2118` `_sol_eclipse_when_glob_pythonic`
  describes the wrong `tret` layout.

### Dead code (log in CLEAN.md — do NOT delete inline)

- `fast_calc.py:~654-690` — dead/partly-incorrect Fukushima-Williams matrix
  builder `_fw2m` / `_iau2006_precession_angles`.
- `lunar.py:384` — discarded computed z-component in
  `_calc_mean_apse_analytical`.
- `lunar.py` — unused helpers (osculating samplers, cubic spline, linear solver,
  unwrap) at ~2392-2552 / 2662-2921 / 3024-3070.

### Documentation (minor)

- `docs/leb/guide.md` — broken section anchor: link to
  `algorithms.md#gravitational-deflection` should be `#9-gravitational-deflection`
  (or drop the number from the header).
- `docs/comparison/intentional-divergences.md` §1 — the word "four" is reused for
  two different body sets; disambiguate.
- `AGENTS.md` / `CLAUDE.md` say base tier = `1849-2150`; README / getting-started
  say `1850-2150`. Harmonize (state de440s raw coverage 1849 vs usable LEB2
  base-core 1850 explicitly).
- `docs/manual/en/09-eclipses.md` — Italian variable names leak into a few
  English examples (`sol_eclipse_where` / `lun_eclipse_when_loc` /
  `lun_occult_when_loc`, ~lines 87/250/353+); translate for full en parity.

### Packaging decisions (confirm, then document)

- For manifest clarity, consider adding the `data/leb2/*.leb2` path to
  `MANIFEST.in` (it already ships correctly via setuptools package-data; this is
  cosmetic so the two manifests don't appear to disagree).

---

## Completed (pre-v3, as of March 2026)

- ~~#1~~ Thread safety — `_STATE_LOCK` RLock on critical setters/getters in state.py
- ~~#2~~ Bare `except Exception` — 101 -> 0 across entire codebase (all narrowed to specific types)
- ~~#3~~ MeeusRangeError — marked deprecated
- ~~#4~~ houses.py dedup — `_init_cardinal_cusps` (9x) + `_set_opposite_cusps` (6x) extracted
- ~~#5~~ constants.py `__all__` — 785 public names, `annotations` excluded
- ~~#6~~ PyPI Packaging — `base_core.leb` bundled in wheel, auto-discovered
- ~~#7~~ Download command — `download_leb2_for_tier()` in download.py, 12 files in DATA_FILES
- ~~#8~~ LEB1 regenerated — Pluto 64d/deg11, Uranians 256d/deg7
- ~~#9~~ GitHub Release data-v2 — 12 LEB2 files published
- ~~#10~~ Validation suite — 8/8 PASS, 61,547 checks, ALL PASS
- ~~#11a~~ Uranian geocentric — added geocentric path in `planets.py`
- ~~#11b~~ Sun heliocentric — returns (0,0,0) correctly now
- ~~#11c~~ True Node distance — full osculating orbit calc (2000x improvement)
- ~~#12~~ Skyfield TypeError — retry + radec removal, cache clear
- ~~#13~~ Developer CLI split — `leph` now has its own `status` command and a development bootstrap `download all` that fetches prerequisites and source inputs, while `libephemeris` stays focused on end-user runtime setup
