# Hypothetical-Body Orbital Elements — Provenance Record

Provenance for the Hamburg-school Uranian planets and the trans-Plutonian
bodies computed by `libephemeris/hypothetical.py`. Orbital elements are
facts/data, not protectable expression; the project's policy is to cite
the published source per body and to enforce, by gate, that the values
used at runtime trace to those sources.

## Single source of truth

`libephemeris/data/fictitious_orbits.csv` holds the published elements. The
live runtime dictionary `URANIAN_KEPLERIAN_ELEMENTS` (and the Proserpina /
Transpluto entries in `HYPOTHETICAL_ELEMENTS`) must agree with the CSV;
this is checked on every run of `scripts/check_hypothetical_provenance.py`
(`poe provenance:hypothetical`, zero-mismatch gate).

`scripts/check_hypothetical_provenance.py --explain` prints the derivations
(it absorbed the role of the former `scripts/derive_hypothetical_elements.py`,
which carried stale hardcoded values and was removed).

## Per-body sources

### Hamburg-school Uranian planets (Cupido … Poseidon)

Witte, A. & Lefeldt, H. (1928). *Regelwerk für Planetenbilder.* Hamburg:
Ludwig Rudolph — original elements; refined by Neely, J. (1988), "The
Uranian Planets," *NCGR Research Journal* vol. 1.

| Body | a (AU) | e | i (°) | ω (°) | Ω (°) | M0 (°) @ J1900 |
|---|---|---|---|---|---|---|
| Cupido | 40.99837 | 0.00460 | 1.0833 | 171.4333 | 129.8325 | 163.7409 |
| Hades | 50.66744 | 0.00245 | 1.0500 | 148.1796 | 161.3339 | 27.6496 |
| Zeus | 59.21436 | 0.00120 | 0.0 | 299.0440 | 0.0 | 165.1232 |
| Kronos | 64.81690 | 0.00305 | 0.0 | 208.8801 | 0.0 | 169.0193 |
| Apollon | 70.29949 | 0.0 | 0.0 | 0.0 | 0.0 | 138.0533 |
| Admetos | 73.62765 | 0.0 | 0.0 | 0.0 | 0.0 | 351.3350 |
| Vulkanus | 77.25568 | 0.0 | 0.0 | 0.0 | 0.0 | 55.8983 |
| Poseidon | 83.66907 | 0.0 | 0.0 | 0.0 | 0.0 | 165.5163 |

Mean motion is the Gaussian value n = k/a^1.5 with k = 0.9856076686 °/day
(not an independent fitted parameter). Gate: exact equality live↔CSV per
field, plus the n formula.

### Transpluto (Isis)

Strubell, M. (1952). *Die Sterne* 3/1952, p. 70ff. CSV row preserves the
Strubell originals (epoch JD 2368547.66, equinox J1945 = JD 2431456.5,
a = 77.775 AU, e = 0.3, ω = 0.7°, M0 = 0.0). The runtime entry uses
elements referred to J2000, **derived** as:

- M0(J2000) = M0 + n·(JD_J2000 − 2368547.66), n = 360/(a^1.5·365.25)
  → 119.265936° (gate reproduces this to <1e-2°; observed Δ ≈ 0).
- ω(J2000) = 0.7° + general precession in longitude J1945→J2000
  = 0.7° + 5028.796195″/cy · 0.55 cy / 3600 → 1.46828° (gate reproduces
  the in-code 1.468282° to <1e-4°; observed Δ ≈ 3e-6°). The precession
  rate is the IAU 2006 value (Capitaine et al. 2003).

### Proserpina

Circular-orbit convention credited to V. Abramov (Tartu); no public
publication of the digits is known, so they are interoperability values
recovered by a black-box output fit against the reference API. CSV
row: epoch J1900, a = 79.225630 AU, circular, M0 = 170.73°. Gate: exact
equality live↔CSV.

### Other predicted/historical bodies

`fictitious_orbits.csv` carries full per-row provenance for Nibiru
(interoperability values; concept circulated via C. Wöltge), Harrington
(1988, AJ 96(4) 1476), the historical Neptune predictions (Leverrier/
Adams, via Hoyt 1980), the historical Pluto predictions (Lowell/
Pickering, via Hoyt 1980), Vulcan (L.H. Weston's published AFA monograph;
see also Le Verrier 1859), Selena/White Moon (Russian-school 7-year
circular convention; digits are interoperability values) and Waldemath
(1898 published claims; the fully-specified digits are interoperability
values). "Interoperability values" means: no public publication of the
digits is known and they are held as black-box fits against
reference-API output. These are propagated directly from the CSV.

## CSV ↔ publication cross-check

The gate verifies code ↔ CSV automatically. The CSV ↔ publication step is a
documented manual check: the Uranian and Transpluto/Proserpina values in
`fictitious_orbits.csv` were transcribed from the editions cited above and
spot-checked against the values in standard astrological-software element
sets that credit the same sources (Neely's published refinement is the
set the reference distribution also credits). For the rows classified as
interoperability values no publication exists to check against; their
provenance is the documented black-box output fit.

## Legacy display tables (disclosed, frozen)

`URANIAN_ELEMENTS` (a mean-longitude oscillation table) and the per-body
`*_KEPLERIAN_ELEMENTS` `L0` constants (and Hades `M0` = 26.850162) are
historical fits **calibrated against pyswisseph output**, not attested in
Hamburg-school literature. They are **no longer consulted at runtime**
(`calc_uranian_longitude` delegates to the published Keplerian set) but are
pinned by tests for module-API stability. The gate freezes them against a
snapshot so any future edit is a conscious reclassification, and they are
disclosed in NOTICE.md ("Calibration Data Disclosure"). The Hades figure
26.850162 is specifically the residue of a retired mean-longitude path
(old unified longitude 336.363662 − ω 148.1796 − Ω 161.3339); the published
mean anomaly is 27.6496, which the live dict uses.

## Verification

```
poe provenance:hypothetical
python scripts/check_hypothetical_provenance.py --explain
pytest tests/test_hypothetical.py -q
pytest tests/test_planets/test_uranian_bodies_comprehensive.py -q
```
