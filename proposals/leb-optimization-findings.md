# LEB optimization: clean-room findings

**Status:** historical measurements sanitized; current guidance only.

The original document contained per-body residuals, generated-file sizes and
unsupported hypothetical-body parameters. Those details were removed because
the corresponding historical artifacts are retired and must not be used as
data or validation inputs.

## Findings that remain valid

- Chebyshev interval and degree choices depend on orbital frequency,
  eccentricity, perturbations, coordinate representation, and requested date
  range; one global setting is inappropriate.
- Optimization must be evaluated against independently sourced mathematical or
  JPL state truth, including worst-case geometry, not against an earlier float
  implementation or compatibility-oracle output.
- Compression and reader performance are separate concerns. Chunked lazy
  decompression improves cold access while the scalar hot path remains a small
  Clenshaw evaluation plus the apparent-place pipeline.
- Vectorization is useful for genuine multi-date workloads but adds overhead to
  isolated scalar calls.
- An optimization is accepted only if it preserves declared error bounds across
  coverage edges and all supported coordinate/flag paths.

## Provenance constraints

- Generate every non-bundled LEB file locally from NASA JPL kernels and cited
  independent models.
- Do not convert or benchmark retired published artifacts.
- Do not include unsupported hypothetical IDs in coefficient files. Harrington
  (ID 50) is the only built-in hypothetical orbit with accepted provenance.
- Never retain public compatibility-oracle values as sweep reports, fitted
  parameters, thresholds, or binary payloads.

## Recommended workflow

1. Select the tier and independently sourced body inventory.
2. Generate candidate coefficients locally.
3. Verify defining conditions, Cartesian state error and boundary continuity
   against the independent source path.
4. Run the focused LEB reader, compression and fast-path tests.
5. Run `uv run python scripts/check_provenance.py` before packaging.

Current commands and format details are maintained in `docs/leb/guide.md`.
