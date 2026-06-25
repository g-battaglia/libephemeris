# Provenance Fingerprint Sweep — June 2026

Record of the repository-wide sweep for Swiss Ephemeris C-source
fingerprints, performed as part of the June 2026 code review remediation
(WS1). Goal: results stay 1:1 with the reference API, but every line of
implementation *expression* (identifiers, comments, docstrings, code
structure) must be independently written.

> **Related:** the Galilean satellite module's separate clean-room rewrite
> (removing the last LGPL/PyMeeus-adapted code) is recorded in
> [galilean-clean-room-2026-06.md](galilean-clean-room-2026-06.md); its
> PyMeeus fingerprint class lives in the same `scripts/check_provenance.py`
> gate (class 3).

## Method

`grep -rniE` over `libephemeris/` (excluding `__pycache__`) for two
classes of fingerprints:

1. **Source-file references**: `swehouse`, `swecl.c`, `sweph.c`,
   `swemmoon`, `swemplan`, `swedate.c`
2. **swehouse.c/swecl.c implementation identifiers**: `swed`, `dgsect`,
   `xs1`, `xh1`, `fh1`, `mdd`, `mdn`, `adp`, `admc`, `samc`, `dfac`,
   `apc_sector`, `xeq0`, `xp0`, `VERY_SMALL`, `AltO`, `ArcusVisionis`

## Triage classes

- **A — allowed**: names of the public reference API or of its documented
  output slots (e.g. the heliacal `dret` slot names `AltO`/`AppAltO`
  documented in `heliacal.py`/`eclipse.py` docstrings, the test class
  `TestSwePheno20Values` naming the `pheno` API). Public interface names
  are not implementation expression.
- **B — coincidental**: generic identifiers that also happen to appear in
  the C source. Allowed outside the remediated cluster; renamed inside it
  for clarity.
- **C — implementation-derived**: identifiers/comments/structure traceable
  to the C implementation. Must be zero after remediation.

## Findings and remediation (2026-06-11)

Class-C hits were confined to `libephemeris/houses.py`, in five sites.
All five were rewritten in place — formulas preserved (verified
bit-identical), expression re-derived from published definitions:

| Site | Remediation |
|---|---|
| `_houses_savard_a` ('J') | Removed "Algorithm (from Swiss Ephemeris swehouse.c)" docstring; re-derived from the prime-vertical/declination-parallel geometry (Smart Ch. 3); renamed `xs1/xs2/xh1/xh2/fh1/fh2` → `pv_arc_*`/`ha_*`/`pole_*` |
| `_houses_krusinski` ('U') | Replaced the `(a−b+540)%360−180` idiom with the project's own `difdeg2n`; renamed spherical triples; sources already cited (Krusinski 1995; Pisa 1997; Goelzer 1995) |
| `_apc_sector` → `_apc_cusp` ('Y') | New descriptive signature (`house, lat_rad, eps_rad, armc_rad`); full derivation written from Knegt's published definition; renamed `a/k/ph/e/az` |
| `_house_pos_pythonic` P/G/K/R/C branches | Renamed `mdd/mdn/de/sinad/ad/sad/san/adp/admc/samc/dfac/xeq0/xp0/is_above_hor/is_invalid` → descriptive names; deleted the verbatim swehouse.c Koch formula comment block; wrote independent derivations (Holden 1977 semi-arc; position-circle closed form; GOH definition; horizon-test derivation); removed the vestigial write-only `is_circumpolar` flag; fixed the `cos_eps`/`sin_eps` shadowing in the Campanus branch (`cos_rot`/`sin_rot`) |
| `_houses_sunshine_makransky` ('i') | Renamed the one-letter locals (`a/b/c/f/zd/pole/q/w/m/z/r/cu/xh/rah/md/nsa/dsa`) to descriptive names; comments re-anchored to Makransky's published procedure (*Primary Directions*, 1988). Standard primary-directions symbols (MD, ZD, AD, pole height) retained — they are the literature's notation, not the C source's |

Additionally every local `VERY_SMALL` constant in `houses.py` was renamed
to `_NEAR_ZERO` (class-B hygiene).

## Verification

`scripts/verify_houses_clean_room.py` (26 house systems × ~91,000 grid
cases of `houses_armc` + `house_pos`, including polar latitudes, southern
hemisphere, cardinal ARMC, nonzero body latitudes):

- `--check` against the pre-remediation baseline: **0 of 91,240 cases
  drifted; worst delta 0.0** (bit-identical outputs).
- `compare` vs pyswisseph 2.10.03: unchanged from the pre-remediation
  summary (`tasks/results/houses_compare_pre_ws1.txt`); the known
  pre-existing deviations (Placidus within ~0.01° of the polar circle,
  Sunshine-alt polar raises, scattered polar `house_pos` cases) are
  tracked separately as parity work, not provenance work.
- `pytest tests/test_houses` — 1,290 passed.

## Final sweep state

After remediation, the sweep reports **zero class-C hits** across
`libephemeris/`. Remaining hits are class A (public API/dret slot names in
docstrings) only.

To re-run:

```bash
grep -rniE "swehouse|swecl\.c|sweph\.c|swemmoon|swemplan|swedate\.c" libephemeris/ --include="*.py"
grep -rnE "\bswed\b|\bdgsect\b|\bxs1\b|\bxh1\b|\bfh1\b|\bmdd\b|\bmdn\b|\badp\b|\badmc\b|\bsamc\b|\bdfac\b|apc_sector|\bxeq0\b|\bxp0\b" libephemeris/ --include="*.py"
.venv/bin/python scripts/verify_houses_clean_room.py --check tasks/results/houses_baseline.json
```
