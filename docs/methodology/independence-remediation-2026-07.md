# Independence remediation — July 2026 (WS2)

This document records the July 2026 provenance remediation sweep ("WS2"),
the successor of the June 2026 sweep recorded in
[provenance-sweep-2026-06.md](provenance-sweep-2026-06.md) (WS1). It was
triggered by an external license audit of v3.0.0rc5 that identified code,
data and history items where the project's declared independence from the
Swiss Ephemeris was stronger than the evidence.

Method note: as in WS1, no Swiss Ephemeris source file was read at any
point during this remediation. The reference implementation is used
exclusively as a black-box oracle (its Python binding, driven in a
separate virtualenv, compared on inputs/outputs only).

## 1. `house_pos` placement branches (code)

**Finding.** The July 2026 parity work (81fca021) added placement branches
whose own comments declared them "1:1 ports" of the reference
implementation (topocentric 'T', horizon 'H', Carter 'F', Krusinski 'U',
Savard 'J', Sripati 'S', Sunshine/APC 'I'/'i'/'Y', a degree-trig macro
layer, the polar-ascendant correction and a coordinate-rotation helper).
A later pass (67d3dd0) only renamed identifiers, which is not a
provenance remedy.

**Remedy (WS1 method).** Every branch was re-derived as independent
expression from the published system definitions and restructured into
standalone documented functions; the coordinate rotation now delegates to
this project's public `utils.cotrans`; shared spherical-astronomy facts
were factored into helpers citing standard literature (Meeus). Public
sources used for the constructions: Polich & Page (Spica, 1964) for the
topocentric position-circle family; C.E.O. Carter for the poli-equatorial
definition; Krusinski/Pisa/Goelzer for the ascendant-zenith great circle;
J.F. Savard for the prime-vertical trisection; B.V. Raman for the Sripati
convention; Makransky (Sunshine) and Knegt (APC) for the semi-arc
placements.

**Verification.** Outputs bit-identical to the pre-rewrite code over a
9000-sample randomized grid (13 systems, latitudes ±75°, obliquity
22.5–24.5°, body latitude ±8°); 12/13 systems exact against the black-box
oracle on the same grid (the Sripati wrap-seam difference predates the
rewrite); full houses suite green.

## 2. Lunar-occultation search skeleton (code)

**Finding.** `lun_occult_when_glob` and the local-search twin shared an
inline search skeleton (type-filter cascade, mean-rate conjunction
stepping, latitude gates, conjunction hop) that tracked the reference
implementation's structure closely.

**Remedy.** The skeleton was re-derived as project expression: the filter
expansion became a documented standalone helper
(`_normalize_occultation_filter`), the conjunction seek a shared helper
(`_seek_moon_conjunction`) used by both searches, and every parameter is
now a named constant with its physical derivation (mean lunar sidereal
rate; sidereal-month spacing for the 20-day hop; Moon maximum ecliptic
latitude + parallax + semidiameter for the 7°/2° gates; IAU 1976 Earth
equatorial radius). Public error strings are unchanged (API contract).

**Verification.** Occultation suite (109) and eclipse suite (1034) green.

## 3. Fictitious-body orbital elements (data)

**Finding.** `data/fictitious_orbits.csv` carried four rows attributed to
private communications or unpublished formulations (Nibiru, Vulcan,
Proserpina, Waldemath) and one with a vague attribution (Selena),
matching the element sets that circulate with the reference distribution.

**Remedy.**

* Where a genuine public publication exists the row now cites it:
  Vulcan (L.H. Weston, "The Planet Vulcan", AFA monograph, published) and
  Waldemath (the 1898 published claims; the circulating fully-specified
  reconstruction is noted as such).
* Where no public source is known (Nibiru, Proserpina, the Selena
  digits), the rows are reclassified as **interoperability values** under
  the calibration-data disclosure of NOTICE.md: a black-box Keplerian
  fit against reference-API output positions recovers these element sets
  without access to any non-public source. Fit procedure: Gauss-Newton
  on (M0, a, e, ω, Ω, i) (static heliocentric rows), (M0, a) (circular
  Proserpina), (L0, rate, radius) (Selena), and the full 9-parameter set
  (M0, M_rate, a, e, ω, ω_rate, Ω, Ω_rate, i) for the geocentric
  Waldemath orbit; targets sampled from oracle output over multi-decade
  spans; recovered values agree with the CSV digits to the fit's model
  floor (RMS ≈ 2–3×10⁻³ deg, dominated by frame conventions in the quick
  forward model, orders below the bodies' own model differences).
* Numeric values were left unchanged, so runtime behaviour and parity
  are untouched.

## 4. Claims (docs)

NOTICE.md and LICENSING.md no longer assert a categorical absence of
license obligations; they state the working standard (independence,
black-box-only validation), the disclosures, and point here for the open
items.

## 5. Adversarial audit round (WS2-H)

After the remediation above, a 48-agent adversarial audit (6 finders with
distinct lenses + per-finding refutation verifiers) swept the tree for
residual taint. It confirmed 42 residual items — all remediated in the
follow-up commits: source-knowledge comment language replaced with
observed-output framing (houses `_apc_cusp`, eclipse shadow scales,
noon-transit flag note, `dcore` layout note, disc-radius table notes); the
`acmc`/`acmc_diff` identifiers renamed and the identifier gate extended
(`acmc`, `plaus_iflag`); the `_plaus_ephemeris_flags` helper renamed to
`_exclusive_ephemeris_bit` with black-box docstrings; magic factors named
and framed as fitted interop constants (`_LUN_UMBRA_SCALE`,
`_LUN_PENUMBRA_SCALE`, `_INNER_CONTACT_MOON_SCALE`); every remaining
`swe_`-prefixed name purged from docstrings; the Waldemath row's black-box
fit actually performed (9-parameter recovery, §3); all residual
categorical independence claims aligned (NOTICE, README,
THIRD_PARTY_NOTICES, comparison hub, release notes, v1.0.0 errata,
hypothetical-bodies.md); a working document containing a pointer into the
reference source distribution removed from the tree; the provenance gate
re-scoped (tracked-files filename gate, root legal docs allowlisted,
reviewed-line waivers) and re-verified at zero together with the SPDX
gate; `uv.lock` regenerated so the `all` extra resolves without GPL
packages.

### Round 2

A second, fresh 24-agent audit over the WS2-H tree confirmed 20 further
residuals (down from 42), all remediated: the remaining internal-mechanism
docstring language in the flag-normalization helpers replaced with
measured-output claims (including the tests that quoted the internal
routine name); the mean-ellipse apsis constants re-framed as textbook
values observable in output; the last product-name-and-version comment
and every ``swe.``-qualified name swept from shipped code; the
JPL-only calibration story in the lunar-apsides methodology docs, the
manuals (EN/IT), the methodology hub and known-differences §1.2 aligned
with the NOTICE oracle-calibration disclosure (the generated table's own
header states the black-box oracle); the unshipped v3.0.0 draft notes and
THIRD_PARTY_NOTICES `[all]` wording corrected; the provenance gate
further hardened (pyproject.toml and release-notes/ in scope with the
license-history classes expected in historical notes, waiver restricted
to the license-naming classes, foreign-data regex widened with data-file
targeting, `git ls-files` failure now fatal).

## 6. Open items

* **Repository history.** The repository history still contains: the
  upstream `seorbel.txt` data file added by ca32b36 (removed from the
  tree later, blob still reachable); pre-WS1 versions of `houses.py`;
  the pre-clean-room Galilean module. WS2-G prepared the remedy: a full
  pre-purge backup bundle, and a `git filter-repo` rewrite (executed on
  a local clone, verified: zero `seorbel.txt` blobs remain at either
  historical path, 1289 commits preserved with rewritten ids). Two
  further options are documented for the derived historical *versions*
  of `houses.py`/`eclipse.py`/the Galilean module: selective blob
  tombstoning, or a fresh-root squash. Executing any of these on the
  public repository rewrites published commit ids and requires a
  coordinated force-push (note: the hosting platform may keep
  previously-pushed objects reachable by hash until support removes
  them), so execution is gated on the maintainer's go-ahead.
* **Output-calibrated data sets.** The lunar apse residual tables and
  perturbation coefficients (NOTICE.md "Calibration Data Disclosure")
  remain output-calibrated: the INTP_* bodies are constructs defined by
  the reference API, so behavioural parity requires fitting to reference
  output. They contain computed positions (facts), not source
  expression. This is disclosed rather than claimed away; regeneration
  against a DE440-only ground truth is possible at a measured parity
  cost and remains an option if a stricter posture is ever required.
* **GPL-licensed optional extras.** `rebound`/`assist` (GPL-3.0-or-later)
  are optional N-body extras; WS2-F removes them from the `all` extra so
  the default install surface stays permissive-only.
* **Legal review.** Whether any residual similarity constitutes a
  derivative work is a legal question; this record documents the
  engineering facts for that review.
