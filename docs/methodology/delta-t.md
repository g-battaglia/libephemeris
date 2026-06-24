# Delta T (ΔT) in libephemeris

This note documents how libephemeris computes **ΔT = TT − UT1**, the slowly
varying, partly unpredictable offset between uniform atomic time (TT) and
clock-from-Earth-rotation time (UT1). ΔT is required by every UT-based
calculation: it maps a civil (UT) instant to the dynamical time at which the
ephemeris and precession/nutation are evaluated. A wrong ΔT shifts **all** body
positions (∝ each body's mean motion) and, through sidereal time, the **house
cusps / Ascendant**.

Implementation: `libephemeris/time_utils.py` (`deltat`, `_deltat_espenak_meeus`)
and the model selector in `libephemeris/state.py`
(`set_delta_t_model` / `get_delta_t_model`).

> **TL;DR.** libephemeris already uses a correct **multi-era (piecewise)** ΔT
> model — IERS-observed for modern dates, Stephenson-Morrison-Hohenkerk 2016 for
> the eclipse era, and a long-term tidal model at the extremes. It does **not**
> use a single diverging formula. A second, fully clean-room model (the classic
> NASA Espenak-Meeus 2006 polynomials) is selectable. **libephemeris never
> imports pyswisseph**; exact Swiss parity, when needed for validation, is
> injected from the *outside* via `set_delta_t_userdef`.

---

## 1. Why ΔT must be piecewise

ΔT is not given by one closed-form expression over all of history. Earth's
rotation is *measured* by atomic clocks since 1955, *reconstructed* from
telescopic observations back to ~1620, *inferred* from ancient eclipse/occultation
records back to ~720 BCE, and *unknowable* in the future (it depends on
unmodelled core-mantle coupling, post-glacial rebound, ice-mass redistribution,
etc.). A single polynomial applied outside its validated window diverges and
produces gross errors — for an Ascendant this can mean being off by a degree or
more. The standard solution everyone uses is a **block architecture**:

| Era | Method | libephemeris source of truth |
|---|---|---|
| 1955 → present (+ short prediction) | IERS / USNO observed (atomic clocks) | Skyfield ΔT table (IERS finals) |
| ~1620 → 1955 | Telescopic-era tabulated values | Skyfield ΔT table / SMH-2016 |
| −720 → ~1620 | Stephenson-Morrison-Hohenkerk 2016 spline (eclipse records) | SMH-2016 |
| before −720 / far future | Long-term tidal model (integrated length-of-day) | SMH-2016 long-term tail |

libephemeris obtains this piecewise model from **Skyfield**, whose `delta_t`
implements Stephenson, Morrison & Hohenkerk (2016) blended with IERS observations
and a long-term tail — the modern scientific standard.

## 2. Proof that libephemeris is already multi-era (not a single formula)

ΔT in seconds, libephemeris default model vs the *single* long-term parabola
`ΔT = −20 + 32·u²` (u = (year−1820)/100) and the observed/consensus value:

| year | **libe ΔT** | single parabola | observed / consensus | active block |
|---|---|---|---|---|
| 2000 | **63.8** | 83.7 | **63.83 (IERS)** | IERS observed |
| 2025 | **69.1** | 114.5 | ~69 (IERS) | IERS observed |
| 1900 | **−2.7** | 0.5 | **−2.8** | telescopic table |
| 0    | 10441 | 10580 | ~10580 (SMH-2016) | SMH-2016 spline |
| −500 | 16939 | 17204 | ~16800 (SMH-2016) | SMH-2016 spline |
| −6000 / +10000 | 198671 / 216871 | 195668 / 214100 | — | long-term tidal tail |

The decisive cell is **year 2000: libephemeris returns 63.8 s (the IERS observed
value), not the 83.7 s of the single parabola.** If libephemeris used "one
formula" it would read 83.7. It does not — it selects the correct block per
epoch, and the block boundaries are continuous (no steps in the 1600→2000 range).

## 3. The ΔT model selector

ΔT resolution order in `deltat()` is:

1. **user-defined** — `set_delta_t_userdef(value_in_days)` (highest priority);
2. **IERS observed** — if enabled via `set_iers_delta_t_enabled(True)`;
3. **the selected model** (this section).

The model is chosen with `set_delta_t_model(...)` / the
`LIBEPHEMERIS_DELTAT_MODEL` environment variable / the `deltat_model` TOML key:

| value | model | notes |
|---|---|---|
| `smh2016` *(default)* | Stephenson-Morrison-Hohenkerk 2016 via Skyfield | the modern scientific standard; piecewise as in §1 |
| `espenak_meeus` | classic NASA **Espenak & Meeus 2006** polynomial set | self-contained clean-room implementation; the multi-segment polynomials (−1999…+3000) used by NASA's eclipse predictions |

```python
import libephemeris as eph
eph.set_delta_t_model("espenak_meeus")   # or "smh2016"
eph.get_delta_t_model()                   # -> "espenak_meeus"
# or, by environment:  LIBEPHEMERIS_DELTAT_MODEL=espenak_meeus
```

Both models are **clean-room**: every coefficient comes from published sources
(the SMH-2016 / Skyfield implementation, and the NASA Espenak-Meeus polynomials).
The two agree closely in the well-observed range and differ mainly at the
extremes (deep past / far future), where ΔT is itself uncertain.

> **Scope — which calculations the selector (and `userdef` / IERS) affect.**
> All three ΔT priorities feed `deltat()` / `deltat_ex()`, and hence the TT used by
> the **LEB / fast / Horizons** position backends (the default is **LEB**), as well
> as eclipses, heliacal events and long-term sidereal time. The **`"skyfield"`**
> position backend is the exception: it derives TT from Skyfield's *own* internal
> ΔT model (SMH-2016 + IERS), so its planetary positions ignore `set_delta_t_model()`,
> `set_delta_t_userdef()` and the IERS toggle. Concretely, selecting `espenak_meeus`
> (or forcing a `userdef` value) shifts positions in the default LEB backend — e.g.
> a `userdef` step of 0.1 day moves the Sun by ~364″ — but leaves `"skyfield"`-mode
> positions unchanged. This is a pre-existing property of all three overrides, not
> specific to the model selector; the validation injection in §4 relies on it via
> the position backends that *do* honour `deltat()`.

> **Note — `espenak_meeus` is *not* "the Swiss model".** Modern Swiss Ephemeris
> (≥ 2.06/2.08, 2017-2018) **also** uses Stephenson-Morrison-Hohenkerk 2016 for
> the pre-1955 era, not Espenak-Meeus. So `espenak_meeus` reproduces the classic
> NASA model, which matches *neither* libephemeris's default *nor* current
> pyswisseph at remote epochs. It exists for users/tools that specifically expect
> the Espenak-Meeus polynomials, not as a Swiss-compatibility switch.

## 4. The clean-room / licensing constraint

**libephemeris must never import `swisseph` / pyswisseph** (which is AGPL / Swiss
dual-licensed). libephemeris is an independent, NASA-data-powered library; pulling
in Swiss code — even an optional `import` inside a ΔT branch — would contaminate
its license. There is therefore **no `swisseph` ΔT mode in the core.**

When exact Swiss-Ephemeris ΔT is needed *for validation*, it is injected from the
**outside** (in the `validation/` repo, which legitimately depends on pyswisseph):

```python
import swisseph as swe
import libephemeris as L
L.set_delta_t_userdef(swe.deltat(jd))   # force Swiss ΔT into libe, per date
# ... compare L.calc_ut(jd, ...) vs swe.calc_ut(jd, ...) ...
L.set_delta_t_userdef(None)             # restore the native model
```

> ⚠️ **Validation caveat.** Do **not** reset Swiss's own ΔT mid-loop with
> `swe.set_delta_t_userdef(-1.0)` / `AS_UNDEF`: it leaves pyswisseph in a polluted
> state that injects a spurious constant ~3553″ (~0.987°) offset into subsequent
> `swe.sidtime`/`swe.deltat` calls. The clean method forces ΔT only into
> libephemeris and lets Swiss use its native `swe.deltat(jd)` (which equals the
> forced value).

## 5. Comparison vs Swiss Ephemeris

Both engines use the **same SMH-2016 family** for the historical era, so they
agree closely where it matters:

| era | libephemeris vs pyswisseph 2.10 |
|---|---|
| 1962 → present | **identical** (both pinned to IERS observed; ΔT(2000)=63.83 s on both) |
| 1620 → 1955 (telescopic) | within the observational scatter (≤ ~7 s; both track the same records) |
| −720 → 1620 | close (both SMH-2016 family) |
| future > ~2027 / deep past | **diverge** — different extrapolations of an unknowable quantity |

The far-future divergence that originally looked alarming is **entirely ΔT
extrapolation**, not an ephemeris/position error. Forcing the *same* ΔT into both
engines (§4) collapses the difference:

| date | Moon Δ, default ΔT | Moon Δ, **same ΔT** | worst planet, same ΔT |
|---|---|---|---|
| 2100 | 1.6″ | **0.000″** | 0.04″ |
| 2200 | 29.3″ | **0.000″** | 0.08″ |
| 2300 | 87.0″ | **0.000″** | 0.06″ |

So at a common ΔT, libephemeris and Swiss positions match to **< 0.1″** (Moon to
~0.001″) even at year 2300 — the engines agree; only the ΔT *extrapolation
choice* differs, and for the future that choice is scientifically arbitrary (no
ground truth exists). libephemeris keeps the modern SMH-2016 model by default;
this is at least as good as Swiss everywhere and strictly more current for the
deep past. See also the Ascendant/sidereal-time discussion in
[sidereal-time-longterm.md](../methodology/sidereal-time-longterm.md) and
[swisseph-comparison.md](../reference/swisseph-comparison.md).

## 6. Validation

The ΔT model and its consequences are validated in the `validation/` repo:

- a dedicated **lib-vs-Swiss with Swiss ΔT** sweep (forces `swe.deltat` into
  libephemeris via `set_delta_t_userdef` and confirms positions match to < 0.1″
  and the Moon to ~0.001″ across the supported range, isolating the pure
  ephemeris-model difference from ΔT);
- the position/house long-term rounds (`run.sh positions` / `run.sh houses`)
  bound the (expected, physical) ΔT-driven divergence at remote epochs and prove
  the precession/sidereal-time model itself is an exact reproduction of the
  published physics at matched ΔT.

## 7. References

1. Stephenson, F.R., Morrison, L.V. & Hohenkerk, C.Y. (2016). *Measurement of the
   Earth's rotation: 720 BC to AD 2015.* Proc. Royal Society A, 472, 20160404.
2. Morrison, L.V., Stephenson, F.R., Hohenkerk, C.Y. & Zawilski, M. (2021).
   *Addendum 2020 to "Measurement of the Earth's rotation: 720 BC to AD 2015".*
   Proc. Royal Society A, 477, 20200776.
3. Espenak, F. & Meeus, J. (2006). *Five Millennium Canon of Solar Eclipses:
   −1999 to +3000.* NASA/TP-2006-214141 (polynomial ΔT expressions).
4. Morrison, L.V. & Stephenson, F.R. (2004). *Historical values of the Earth's
   clock error ΔT and the calculation of eclipses.* J. Hist. Astron. 35, 327.
5. IERS Earth Orientation Centre — observed ΔT / UT1 bulletins.
