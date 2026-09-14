# Independent review — interval matrix primitives

Reviewer: independent bounded review  
Commits reviewed: `89898503` and `5b18669d`  
Scoped implementation: `libephemeris/intervals.py`  
Scoped tests: `tests/test_intervals.py`  
Specification basis: `validation/golden/specs/occultation_contact_residuals/ellipsoid.md`, sections 7–10 only

## Verdict: REJECTED

The Krawczyk validation changes are sound for the checked contracts, and the interval Cholesky recurrence is mathematically sound for finite symmetric interval inputs. One blocking input-validation hole remains: `interval_cholesky()` accepts non-finite Arb entries. The helper therefore does not uphold a finite certified SPD-input contract and its certification wording is too strong for all inputs it currently accepts.

## Findings

### 1. Blocking — Cholesky has no finite-input guard

`interval_cholesky()` validates non-empty square shape and compares every off-diagonal pair for exact interval symmetry, then requires each computed pivot to be strictly positive. It never checks that the matrix entries are finite. A tiny probe passed a matrix whose second diagonal entry was `arb("inf")`; the function returned a factor containing infinite intervals instead of rejecting the input. An infinite or non-finite matrix is not a finite positive-definite interval matrix, so this is not a valid Cholesky certificate.

Required disposition: reject any non-finite matrix entry before the recurrence (including NaN and either infinity), or otherwise make the finite-input precondition explicit and enforce it. Do not treat a positive infinite pivot as an SPD certificate.

### 2. Cholesky recurrence, symmetry, and pivots

For finite inputs, the scalar lower-triangular recurrence is correctly evaluated with Arb operations. The implementation uses the lower triangle only after requiring each corresponding upper/lower interval pair to compare equal; it does not silently symmetrize an asymmetric enclosure. Every pivot is tested with strict `> 0`, so zero-touching and indefinite pivots fail closed with `IntervalCertificationError`. Division is performed only by a previously certified positive diagonal pivot. The factor is lower triangular, and a probe confirmed that its product encloses the input on an uncertain finite symmetric example.

The symmetry test is conservative but sound: overlapping, non-identical lower and upper enclosures are not accepted as a proof of symmetry. The recurrence does not mutate the input matrix; the probe compared the matrix before and after execution and found it unchanged.

### 3. Krawczyk guards and operator

The Krawczyk helper has the correct operator:

```text
K = x - C F(x) + (I - C J(X)) (X - x).
```

Dimension checks require one non-empty square system. The center and preconditioner are required to be finite exact point quantities; function values, Jacobian entries, and domain entries are required to be finite intervals. The determinant check rejects a singular or determinant-uncertain preconditioner, and `domain.contains(center)` correctly permits a center on the closed-box boundary while rejecting a center outside the domain. The helper computes the image without mutating the center, matrices, or domain inputs; a probe confirmed all supplied objects remained unchanged. The exact linear contraction test and the invalid-input tests both pass.

No mathematical defect was found in the implemented formula or these Krawczyk guards.

### 4. Existence/uniqueness wording

The Krawczyk docstring is appropriately conditional: it says strict inclusion proves a root only “under the usual Krawczyk theorem,” and explicitly states that this helper computes only the enclosure while callers retain responsibility for the broader certificate. It does not silently claim that merely computing an image establishes existence or uniqueness without the theorem hypotheses. For maximum precision, “one root” could be expanded to distinguish existence from the standard uniqueness conclusion, but this is not the review blocker.

The Cholesky return wording (“a factor with `matrix = L * L.T`”) should be read as an enclosure statement for each finite SPD realization, not literal equality of interval objects. The missing finite guard makes that wording materially overbroad for the currently accepted non-finite input and is covered by blocker 1.

## Verification performed

- `uv run pytest tests/test_intervals.py -v` — **33 passed**.
- `uv run ruff check libephemeris/intervals.py` — **passed**.
- `uv run mypy libephemeris/intervals.py` — **passed** (one unrelated informational note from `libephemeris/state.py`).
- Tiny probes checked finite uncertain Cholesky enclosure and input immutability, Krawczyk input immutability, and non-finite Cholesky behavior. The non-finite probe reproduced the blocker described above.

No implementation files were edited. The verdict remains **REJECTED** until the Cholesky finite-entry guard is added and covered by a targeted test.

## Re-review of `fc7db5f5`

The finite-entry guard is now present before symmetry checks and the Cholesky recurrence. Parameterized tests cover `arb("nan")`, `arb("inf")`, and `arb("-inf")`; each is rejected with `ValueError` before recurrence. Verification passed: `uv run pytest tests/test_intervals.py -v` (**36 passed**), `uv run ruff check libephemeris/intervals.py tests/test_intervals.py`, and `uv run mypy libephemeris/intervals.py` (with the existing informational note in `libephemeris/state.py`). Verdict: **ACCEPTED**.
