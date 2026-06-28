# Galilean Module Clean-Room Rewrite — Evidence Record (June 2026)

## Summary

`libephemeris/moon_theories/galilean.py` was, through v2.1.0, an adaptation
of PyMeeus's `JupiterMoons.py` and carried `SPDX-License-Identifier:
LGPL-3.0-only`. To make every shipped line owned by the project (or
permissively licensed), the module was **rewritten clean-room** from the
published theory and the LGPL exception retired. (At the time of the rewrite
the project was dual-licensed AGPL-3.0/commercial; it was relicensed to
Apache-2.0 for v3.0.0.)

The rewrite reproduces the prior numeric output to **floating-point
re-association level** (worst per-component moon difference 2.85×10⁻⁹ km
≈ 3 picometres; worst Jupiter COB difference 1.9×10⁻²¹ AU) over a
9146-epoch grid spanning 1800–2200, so 1:1 pyswisseph parity is preserved.

## Theory sources (the only references for the rewrite)

- Lieske, J.H. (1998). "Galilean satellite ephemerides E5." *A&AS* 129,
  205–217.
- Meeus, J. (1998). *Astronomical Algorithms*, 2nd ed., Willmann-Bell,
  ch. 44, which prints the complete algorithm and coefficient tables.

The mathematics and the numeric coefficients are published facts; copyright
protects code *expression* (structure, identifiers, comments), not the
theory or its tabulated constants.

## Process — information barriers

The work followed the WS1 houses protocol
(`provenance-sweep-2026-06.md`) with strict role separation; each role ran
as a separate agent session so the transcripts are the audit trail.

| Phase | Role | Inputs allowed | Output |
|---|---|---|---|
| G0 | baseline | current module (numeric only) | `scripts/verify_galilean_clean_room.py`, `tasks/results/galilean_baseline.json` |
| GA | spec author (dirty room) | current module, Lieske/Meeus | `docs/methodology/galilean-e5-spec.md` |
| GB | implementer (clean room) | the spec + the checker's numeric output **only** | the new `galilean.py` |
| GC | auditor | new module + PyMeeus (fetched outside the tree) | this record, the provenance class-3 gate |

The **implementer never saw** the old `galilean.py`, its git history,
PyMeeus, or `THIRD_PARTY_NOTICES.md` — only the functional specification
(`galilean-e5-spec.md`), which restates the published theory and the
project's frame adaptations as mathematics, carrying no identifiers,
comments, or code structure from any prior implementation. The spec's
coefficient tables were mechanically verified against the recorded module
before the implementer began (per-series amplitude/argument fingerprint
multisets matched exactly for all 12 series; 220 periodic terms).

## Baseline & numeric acceptance

`scripts/verify_galilean_clean_room.py` records, per epoch, the four moon
vectors (km, Jupiter-centric J2000 ecliptic) **and** the Jupiter COB offset
from `moon_theories.constants` (AU, equatorial ICRF — this also pins the
GM weighting and the ecliptic→equatorial rotation downstream).

- **Grid:** 1800-01-01 .. ~2200-01 every 16 days (9132 epochs) + 10 golden
  epochs (E5 theory epoch, J2000, compare-suite dates) + 4 ±1 s pairs
  (velocity central-difference path). 9146 epochs total.
- **Gates:** moon |Δ| ≤ 1×10⁻³ km; COB |Δ| ≤ 1×10⁻¹² AU.
- **Result (new vs baseline):** 0 drifted; worst moon |Δ| = 2.852×10⁻⁹ km
  (Callisto.x @ JD 2392432.5); worst COB |Δ| = 1.906×10⁻²¹ AU.

The residual is broadband at the picometre level — pure floating-point
re-association from independently structured code — not a periodic
coefficient signature, confirming the coefficient tables are transcribed
correctly. The gate (1 m) is far tighter than pyswisseph parity needs:
through the ~2×10⁻⁴ COB dilution and Jupiter's geocentric distance, even
1 km of moon drift would move Jupiter's geocentric longitude by only
~7×10⁻¹¹ arcsec, vs the 3.6 arcsec compare-suite tolerance.

## Expression-independence audit (PyMeeus)

PyMeeus 0.5.12 was installed into `/tmp` (outside the repo) **after** the
implementation was complete. An AST identifier/function extraction over
`pymeeus/JupiterMoons.py` was intersected with the new module:

- **Function/class names — zero overlap.** PyMeeus: `JupiterMoons`,
  `rectangular_positions_jovian_equatorial`,
  `apparent_rectangular_coordinates`, `correct_rectangular_positions`,
  `jupiter_system_angles`, `check_coordinates`, … New module:
  `galilean_moon_positions`, `_E5Arguments`, `_fundamental_arguments`,
  `_sum_sine_series`, `_sum_cosine_series`, `_j2000_ecliptic_km`,
  `_to_radians`. A substring scan for every distinctive PyMeeus name over
  the new module returned no hits.
- **Identifier overlap — 4 names, all Class A (published facts/notation):**
  `callisto`, `europa`, `ganymede` (satellite names) and `psi` (the
  Lieske/Meeus symbol for the ascending node of Jupiter's equator on the
  ecliptic). No PyMeeus-distinctive (Class C) identifier appears.

The new module is **table-driven** (module-level coefficient tuples +
generic `_sum_sine_series` / `_sum_cosine_series` evaluators + an
`_E5Arguments` dataclass), structurally unlike PyMeeus's inline
per-coordinate expression style.

## Permanent gate

`scripts/check_provenance.py` gained **class 3** — a zero-hit sweep over
`libephemeris/**/*.py` for PyMeeus implementation identifiers (`pymeeus`,
`JupiterMoons`, `rectangular_positions_jovian_equatorial`,
`correct_rectangular_positions`, `apparent_rectangular_coordinates`,
`jupiter_system_angles`) — so the independence is enforced on every push,
alongside the existing Swiss Ephemeris classes.

## Verification commands

```
.venv/bin/python scripts/verify_galilean_clean_room.py --check tasks/results/galilean_baseline.json
pytest tests/test_planetary_moons.py tests/test_gas_giant_planet_center.py -q
poe license:check        # galilean.py now Apache-2.0; no SPDX exception
poe provenance:sweep     # 0 hits incl. the PyMeeus class
poe wheel:audit          # galilean.py still required-and-present
poe lint && poe typecheck
```

All green; the typecheck baseline is unchanged (galilean.py contributes no
errors).
