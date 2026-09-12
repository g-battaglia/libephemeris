# Independent verifier D — Gauquelin oracle core

**Verdict: REJECTED (negative).** This report audits the standalone numerical core at the requested target state (the checkout contains no resolvable object for `0ab5c46bec31b6258bab9659e48a2c2406dfa270`; the audited working-tree files were inspected read-only). The approved specification is `bfcd0c0`, with review round 5 and acceptance addendum/reviews read. No production Gauquelin/root/chart helpers, restricted materials, reference outputs, or audit finding materials were inspected.

## Commands and test result

- `PYTHONPATH=validation uv run pytest validation/tests/test_golden_gauquelin_oracle.py -vv`: **40 passed, 10 skipped in 1.86 s**.
- `PYTHONPATH=validation uv run python -m py_compile validation/golden/gauquelin_oracle.py validation/tests/test_golden_gauquelin_oracle.py`: passed.
- `uv run python scripts/check_provenance.py`: not run because this standalone validation checkout has no such script at the requested path; no install or lock files were changed.

All ten skips come from the single parameterized ordinary-domain guard (`lat = ±66.56`, `eps = 60`), correctly excluded because `abs(lat)+eps > 90`; the tangent/polar tests are separate and do not claim ordinary-core coverage.

## Positive independent evidence

I decoded binary64 values with `struct`/`Fraction` and independently evaluated representative residuals with Arb, without calling oracle helpers for the residual check. For both latitude signs, ARMC values `0`, `17`, `90`, `179.999`, and `270`, and all 32 non-cardinal sectors, returned candidates had direct Arb residual magnitude below `1e-8` degree. Sectors 5 and 15 and equatorial (`lat=0`) complete-ring paths execute. Flint precision contexts are restored (`ctx.prec` returns to 53), and a two-thread probe showed per-thread precision isolation in the installed Flint build. Exact binary64 ingestion, signed zero, subnormal ARMC, exact `A=0`, nonlinear `lat=0`, endpoint signs, and synthetic unresolved midpoint/work-cap/disagreement tests are present and passed.

## Blocking defects

### 1. Ordinary valid ARMC turns fail, violating the explicit whole-turn contract

The core normalizes the alpha interval with a single `turn = floor(float(midpoint / 360))`, then `_select_longitude` only translates once when the resulting interval is wholly below zero or above 360. This uses native `float` for turn selection and does not robustly normalize arbitrary whole-turn ARMC values. Reproduction (all are mathematically equivalent modulo 360):

```text
certify_boundary(-943.3050469559873, 51.39275650212136,
                 3.5298003344289213, 29)
# OracleCertificationError: longitude rounding cell is unresolved

certify_boundary(720.0, -66.56, 23.4392911, 2)
# OracleCertificationError: longitude rounding cell is unresolved

certify_boundary(-1080.0, -66.56, 23.4392911, 2)
# OracleCertificationError: longitude rounding cell is unresolved
```

The same failure occurs for `derive_midheaven(1080.0, 23.4392911)`, while equivalent turns such as `0`, `360`, and `720` can inconsistently succeed/fail depending on the resulting enclosure. This is a direct failure of spec §3.81 and baseline-free check 17, and is not an acceptable fail-closed near-tangent case: the inputs are strict ordinary and the answer is unambiguously in a binary64 cell. The interval itself can be translated by +1440 to `[306.246..., 306.246...]` and `_candidate_from_interval` then returns `306.2464419855406`; the existing normalization simply never performs that required translation.

### 2. The advertised 160/256 passes are not independent same-proof certificates

`certify_boundary` invokes `_boundary_pass(..., 160)` and `_boundary_pass(..., 256)`, but `_boundary_pass` returns a `RootProof` whose interval is produced by `_root_for_sector`, and `_isolate_monotonic` stops at a precision-independent width threshold (`2**-max(96, start_precision//2)`). The root enclosures and discard histories differ materially between passes (for ARMC 0, lat -66.56, eps 23.4392911, sector 2: 104 discarded intervals at 160 bits versus 136 at 256; the `root` Fraction tuples are unequal). The final check compares only candidate and family, not root interval, endpoint signs, discarded certificates, derivative certificate, chart, seam, or longitude enclosure. Therefore the implementation does not satisfy the addendum/spec requirement that both complete proofs establish the same refined root interval and proof facts; it merely requires equal final float/family.

### 3. The generic isolator accepts a non-monotonic residual despite a caller-supplied derivative label

`certify_monotonic_root` accepts `derivative_lower=Fraction(1)` as metadata but never uses it to prove a derivative enclosure. Its generic `_isolate_monotonic` only evaluates the midpoint and endpoint signs. A synthetic residual with three roots,

```python
def residual(x, _precision):
    y = arb(x.numerator) / x.denominator
    return y * (y*y - arb(1)/4)
```

is accepted over `[-2, 2]` (and `[-1, 1]`) as an exact root at zero, despite multiple roots in the complete interval and no derivative proof. The production `_root_for_sector` does establish a global `|A|<1` lower bound, so this defect is primarily in the exposed certification utility and its test contract; nevertheless it demonstrates that the reported generic certificate is not rigorous and that `derivative_lower` is not evidence.

### 4. Root-discard records do not prove every discarded interval excludes zero

Each discard is recorded from the sign of the midpoint, e.g. `(mid, hi, +1)` or `(lo, mid, -1)`, but the implementation never evaluates an Arb residual enclosure over the discarded interval. A midpoint sign plus an asserted monotonicity label is not the required interval exclusion proof. The production path's global lower bound is checked only as `lower > 0`; `_isolate_monotonic` then records sign-only discarded intervals. This fails the explicit requirement that every discarded interval has a residual enclosure excluding zero and that derivative exclusion be tied to the complete interval.

### 5. Chart/rounding normalization uses native float decisions and has no nonlinear midpoint proof

`_seams_between`, chart-center selection, and longitude-turn selection use `float` conversions and `math.floor`; `_candidate_from_interval` obtains a native-float midpoint to enumerate cells. While exact cell containment is subsequently checked for ordinary cases, the implementation has no nonlinear exact-midpoint identity path. `_select_longitude` can only strict-cell certify or raise; it does not implement the spec's separately required exact nonlinear midpoint proof. The tests cover a synthetic unresolved midpoint and linear exact tie policy, but do not inject a nonlinear exact midpoint or verify rejection/identity at the chart-normalization boundary.

### 6. Cardinal anchor tests bypass the rising proof and hide direct anchor failures

`derive_anchors` special-cases `a % 180 == 90` and sets `asc` equal to normalized ARMC without running the rising-branch proof. Direct `derive_ascendant(90, 20, 23.4392911)` fails with `OracleCertificationError: unresolved exact midpoint sign`, while `derive_anchors(90, 20, 23.4392911)` returns `(90.0, 90.0)` solely through the bypass. Likewise ARMC 270 direct rising derivation fails for nonzero latitude, and `derive_anchors` still passes. The specification requires the rising branch/horizon-plane derivation and independent anchor proof; a test of the wrapper's hard-coded cardinal shortcut is not evidence for that proof. This also means the claimed “both independent anchor passes” do not certify the same mathematical rising derivation at cardinal ARMC.

### 7. Complete ring enforces only inequality, not half-turn mathematical agreement

`certify_gauquelin` assigns opposite cardinal values from rounded anchor floats and then checks only `values[index] != values[index+18]`. It does not independently certify the stated half-turn identity for non-cardinal pairs; the comment claims each opposite root is freshly certified, but no identity or independent relation is actually checked. Existing tests similarly check inequality and cardinal float arithmetic, not the mathematical half-turn relation against independently evaluated opposite roots.

## Test adequacy and skips

The focused tests are useful smoke tests but do not meaningfully cover all requested contracts. They do not independently verify all 32 sectors against a separately coded formula on a grid, do not test both signs across seam crossings with exact seam roots, do not test all whole-turn ARMC normalizations, do not inject a true 160/256 proof disagreement, do not exercise a nonlinear exact midpoint, and do not assert every discarded interval's interval enclosure. The 10 skips are valid strict-ordinary exclusions only; they do not excuse the missing required tangent/polar matrix, which is outside this standalone ordinary oracle and should be tested separately.

## Required disposition

**REJECTED.** Fix exact multi-turn/unwrapped normalization without native-float proof decisions; make both precision passes rerun and compare complete proof facts; certify discarded intervals by interval residual plus derivative exclusion; remove or rigorously constrain the generic non-monotone utility; prove cardinal rising anchors rather than bypassing them; implement or explicitly fail closed on nonlinear midpoint equality; and add independent grid/seam/multi-turn/precision-disagreement tests. No production or core/test files were edited. This report is the only file added by this audit.

Target audited files: `/Users/giacomo/dev/libephemeris/validation/golden/gauquelin_oracle.py`, `/Users/giacomo/dev/libephemeris/validation/tests/test_golden_gauquelin_oracle.py`.
