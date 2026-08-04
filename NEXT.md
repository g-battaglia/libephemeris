# NEXT — Roadmap beyond v3

This document records the work planned after the v3.0.0 release: the path to
the highest achievable precision and coverage across every surface, including
the niche ones. It is a living plan, not a commitment to dates.

Items that shipped code in v3 already point here from the source ("planned for
a future release; see NEXT.md"): notably the composite `SIDBIT` sidereal
projections (§Wave 3.1). Those guards currently raise or warn rather than return
a plausible-looking wrong value; the waves below turn them into full features.
(`nod_aps` for minor bodies — §Wave 1.3 — has since shipped: asteroids/centaurs
now compute osculating nodes/apsides from the shared position pipeline.)

---

## 1. Where we start

State at the v3.0.0 release:

- **Parity vs the reference engine: 89 GREEN + 22 CERTIFIED = 111/111, 0 RED.**
  The 22 CERTIFIED components are all deliberate, documented scientific
  divergences (see `docs/comparison/intentional-divergences.md`), not defects.
- **Unit suite: green** on both backends (Skyfield and the precomputed LEB
  Chebyshev backend).
- **Comparison suite: clean.**

### Binding philosophy

**We do not copy the reference: we stay independent, ideally more precise.**

- 1:1 parity applies to **behaviour, API and semantics**.
- Where we are more scientific / better aligned with modern astronomical
  standards (ΔT from IERS + a peer-reviewed long-term model, a continuous GMST
  join, true cusp derivatives, modern star catalogues, true rates in the speed
  slots) we **do not replicate the reference's data or its defects**.
- Precision is demonstrated against **independent arbiters** (JPL Horizons,
  IERS, Astropy/ERFA, peer-reviewed literature), never against the reference
  itself. Reference-API calls may check public behavior ephemerally, but their
  values may not be saved as fixtures, fitted into corrections, or used to
  derive shipped models.
- Published models (e.g. Schaefer 1993 for heliacal visibility, historical
  Ptolemaic/Babylonian criteria) are implemented **from the literature**.

---

## 2. Known scientific divergences (the state we defend, not "bugs to fix")

Documented in `docs/comparison/intentional-divergences.md` and
`docs/comparison/known-differences.md`. Summary:

| Topic | Bounded difference | Independent basis |
|---|---|---|
| Modern ΔT | Small interpolation residual | IERS observations |
| Historic/future ΔT | Growing extrapolation uncertainty | Stephenson–Morrison–Hohenkerk reconstruction |
| Remote-epoch GMST / houses | Smooth model-dependent drift | IAU 2006 plus continuous Vondrák long-term realization |
| House and apside speeds | Derivative-convention differences | Derivative of the reported position, checked by step halving |
| Planetary nodes | Occasional apparent-place convention differences | Independent orbital geometry and IAU aberration/deflection |
| Star distance/radial speed | Catalog and semantic differences | Hipparcos/Gaia/Astropy space motion |
| Minor-body historic orbits | Orbit-solution differences | JPL SPK/Horizons or optional ASSIST integration |
| Near-polar rise/set | Search-policy differences | Documented two-day physical search window |
| Heliacal visibility | Atmosphere/observer-model differences | Schaefer VISLIMIT, Reijs scale heights, and published optics relations |

---

## 3. The six waves

### Wave 1 — Minor bodies: total coverage
*Goal: any numbered asteroid is first-class, surpassing the reference (which
needs a per-body file).*

- **1.1 Cartesian-state helper** `calc_minor_body_state(ipl, jd_tt) -> (r, v)`
  heliocentric ecliptic-J2000, with the chain SPK type21 → SPK type2/3 →
  ASSIST n-body → Keplerian elements. Today the public surface returns only
  spherical coordinates.
- **1.2 Universal ASSIST**: generalise the curated-asteroid ASSIST path to any
  small body via the SBDB element fetch; add an integration cache and an on-disk
  SBDB cache with a TTL.
- **1.3 `nod_aps` for minor bodies** — DONE. The raise-guard for
  asteroids/centaurs is gone: the heliocentric geometric state is sourced from
  the shared position pipeline (`_minor_body_helio_icrs_state`), reduced with
  μ = GM_Sun through the same machinery as the planets (aberration, deflection,
  J2000/sidereal reduction), with every method falling back to osculating as
  Pluto does and NODBIT_OSCU_BAR re-referencing trans-jovian centaurs to the
  barycentre. Property/invariant tests cover six bodies × three epochs ×
  methods; the like-for-like reference agreement is pinned in the compare
  suite. A future refinement is universal ASSIST states (1.1/1.2) to tighten
  the osculating-realization residual on the apsidal longitudes.
- **1.4 Universal `pheno`**: H, G from SBDB for bodies outside the curated
  table; diameter from SBDB when available.
- **1.5 Edge contract**: with universal ASSIST the out-of-SPK-coverage region
  becomes accurate. Enforce with a `strict_precision` opt-in that raises when
  neither SPK nor ASSIST is available, with an explicit opt-out for the
  degraded Keplerian fallback. Publish a source→precision-by-epoch matrix.

### Wave 2 — Planetary moons out of the box
*The infrastructure already exists (NAIF id parity, the full light-time /
aberration / sidereal pipeline, SPK registration with auto-detect).*

- **2.1 NAIF kernel auto-download** keyed by parent planet, modelled on the
  existing data downloaders, hooked into the moon dispatch.
- **2.2 `pheno` for moons**: a radius table for the ~20 major moons plus a
  dedicated branch (exact phase/elongation/diameter; best-effort magnitude from
  published models, documented as such).
- **2.3 Missing names/constants** for unmapped ids and planet-center ids.
- **2.4 Validation** against **Horizons** (an independent arbiter — the
  reference needs satellite files nobody ships): six moons × three epochs,
  target < 1" geocentric.

### Wave 3 — SIDBIT, rise_trans, remaining conventions

- **3.1 `SIDBIT`** (resolves the composite-mode guard in `set_sid_mode`): today
  the projection bits are accepted with a warning and the base mode is used.
  Implement the `ECL_T0` / `ECL_DATE` / `SSY_PLANE` / `PREC_ORIG` projections
  from published geometric definitions and verify with independent frame
  invariants. Compatibility calls, if used, remain ephemeral.
- **3.2 `rise_trans`: full BIT_HINDU_RISING** — only the geocentric-no-ecliptic-
  latitude component is missing; the rest of the rsmi bits are already handled.
- **3.3 `FLG_JPLHOR` / `JPLHOR_APPROX` policy**: bits accepted and documented as
  no-ops (our positions are already DE440/ICRS-native); a true frame-fix
  emulation is optional and out of scope.
- **3.4 Close the MEAN planetary apside speed investigation**: search for an
  independent published convention; if none exists, keep the true derivative
  and document the bounded difference without deriving a correction from
  compatibility output.
- **3.5 Per-flag certification of minor surfaces**: micro-rounds for per-bit
  `rise_trans`, `refrac`/`refrac_extended` on a P/T/altitude grid, `azalt`
  round-trip, `gauquelin_sector` 0–5.

### Wave 4 — Heliacal: continued physical-model consolidation

The core now uses the published Schaefer VISLIMIT family, Reijs atmospheric
scale heights, a natural-sky baseline, and physical pupil/optics relations.
Future work can consolidate the remaining duplicated search and extinction
paths without changing that provenance boundary. Target architecture:

```
extinction.py (physics source) → visibility.py NEW (vectorised VisibilityEngine)
→ heliacal_search.py NEW (single search skeleton, LEB/Skyfield PositionProvider)
→ heliacal.py (API + dispatch) → schaefer.py (deprecated shim)
```

Each phase must land with independent defining-condition and metamorphic tests:
the unified search (twilight session as the unit, boundary bisection, coherent
`jd1`/`jd2`/`jd3`, conjunction tracking), backend equivalence, and validation
against the literature (Schaefer 1993, Krisciunas & Schaefer 1991, Yallop 1997,
and the Venus tablets of Ammisaduqa). No golden capture from a compatibility
oracle may enter the repository.

### Wave 5 — Certification: dual verdict

- **5.1 Dual-verdict validator**: alongside parity-vs-reference, an
  accuracy-vs-arbiter column (Horizons, astropy/ERFA, IERS, literature).
- **5.2 Divergence charter**: promote `intentional-divergences.md` to a
  versioned user document, with numeric arbiter evidence per entry.
- **5.3** — done in v3: the CERTIFIED verdict in the validation harness.

### Wave 6 — Standard-grade quality

- **6.1 Zero-setup**: `pip install` and it works to ±13000 years; unified
  downloaders; a `leph doctor` command that verifies/fetches everything and
  prints the precision-by-body-class matrix.
- **6.2 Continuous certification**: a CI validator profile plus a published
  per-release matrix. "Verified: N GREEN + M CERTIFIED / 111" badge.
- **6.3 Batch API**: `calc_ut_many(jds, planet, flags) -> ndarray (N, 6)`
  vectorising the LEB branch only, scalar fallback for Horizons/Skyfield.
  Target > 10× on the equivalent loop.
- **6.4 Hardening**: guard case for the Bennu-1980 boundary; cache hygiene
  (add the long-term sidereal LRU caches and the star-catalogue rebuild to
  `clear_caches()`).

---

## 4. Suggested order

Wave 1 → 2 → 3 → 5 → 4 (the long build) → 6, with 6.1 startable in parallel.

Each wave verifies with new targeted tests plus green `test:skyfield:fast` and
`test:leb:fast`, dedicated rounds/verify in the validation repo, a `./run.sh all`
at the end of the wave, and clean lint/format/typecheck on the touched files.
