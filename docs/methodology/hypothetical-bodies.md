# Hypothetical bodies

LibEphemeris recognises all historical compatibility IDs 40–58. Thirteen IDs
have independently reconstructed numerical models: the eight Hamburg-school
points (40–47), Harrington (50), Le Verrier (51), Adams (52), Lowell (53), and
the symbolic White Moon / Selena convention (56). The remaining six IDs retain
their names and constants but deliberately
raise `UnknownBodyError`; their exact missing evidence is inventoried in
[the missing-models inventory](missing-hypothetical-models.md).

These are historical mathematical, predictive, or astrological conventions.
Their inclusion does not assert that the proposed objects physically exist.

## Per-body status

| ID(s) | Bodies | Built-in status and primary source |
|---:|---|---|
| 40–47 | Cupido, Hades, Zeus, Kronos, Apollon, Admetos, Vulkanus, Poseidon | Supported; James Neely, *Matrix Magazine* VII (1980), Table I, p. 10. |
| 48 | Transpluto / Isis | Unsupported; the exact Strubell-attributed element convention was not recovered. |
| 49 | Nibiru | Unsupported; no credible source-complete primary orbit was identified. |
| 50 | Harrington | Supported; Harrington (1988), *AJ* 96, p. 1478. |
| 51 | Le Verrier | Supported; Le Verrier (1846), *Comptes rendus* 23, p. 432, with the more precise values in his *Recherches*. |
| 52 | Adams | Supported; Adams (1846), sections 47–48, printed pp. 25–26. |
| 53 | Lowell | Supported; Lowell (1915), pp. 5, 9, and 105. |
| 54 | Pickering | Unsupported; the inspected circular refers elsewhere for complete elements and does not supply a complete row. |
| 55 | Vulcan | Unsupported; the exact Weston-attributed convention was not recovered. |
| 56 | White Moon / Selena | Supported; Velichko and Larin (2007), pp. 17, 18, 20, 29, and 45, plus explicit LibEphemeris time/radius conventions derived from IAU standards. |
| 57 | Proserpina | Unsupported; no complete primary Abramov publication was recovered. |
| 58 | Waldemath | Unsupported; no source-complete reconstruction of the historical second-moon convention was recovered. |

`HYPOTHETICAL_PROVENANCE` exposes the same boundary in machine-readable form.
Unsupported legacy element-container objects remain importable for API
compatibility, but every unavailable numeric field is `NaN`; they are not
registered with a calculation path.

## IDs 40–47: Neely's Hamburg-school elements

### Primary transcription

The source is James Neely, “Orbital Elements for the Transneptunian Planets,”
*Matrix Magazine*, issue VII (1980), pp. 6–12, Table I on p. 10. A public scan
is catalogued as [Alexandria iBase item 18](https://special-collections.alexandriaibase.org/items/show/18).
The reviewed 63-page scan has SHA-256
`f3d2e834a4a684736383ec387c9c2786064ca4668256724865da2f3d8fa03020`.

Neely defines `T` as Julian centuries from JD 2415020.0. Table I prints
semimajor axis `a`, eccentricity `e`, inclination `i`, argument of perihelion
`omega`, mean anomaly `M0 + rate*T`, and ascending node
`Omega0 + 1.3960*T`. The literal `T=0` transcription is:

| ID | Name | a (AU) | e | i (deg) | omega (deg) | M0 (deg) | dM/dT (deg/century) | Omega0 (deg) |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 40 | Cupido | 40.99837 | 0.00460 | 1.0833 | 171.4333 | 163.7409 | 137.13640 | 129.8325 |
| 41 | Hades | 50.66744 | 0.00245 | 1.0500 | 148.1796 | 27.6496 | 99.81803 | 161.3339 |
| 42 | Zeus | 59.21436 | 0.00120 | 0 | 299.0440 | 165.1232 | 79.00633 | 0 |
| 43 | Kronos | 64.81690 | 0.00305 | 0 | 208.8801 | 169.0193 | 68.98746 | 0 |
| 44 | Apollon | 70.29949 | 0 | 0 | 138.0533 | 0 | 61.0765 | 0 |
| 45 | Admetos | 73.62765 | 0 | 0 | 351.3350 | 0 | 56.98245 | 0 |
| 46 | Vulkanus | 77.25568 | 0 | 0 | 55.8983 | 0 | 53.015987 | 0 |
| 47 | Poseidon | 83.66907 | 0 | 0 | 165.5163 | 0 | 47.03868 | 0 |

### Frame and propagation decisions

The runtime keeps the printed `T=0` angles in a fixed J1900 ecliptic frame and
then applies LibEphemeris's independently documented Vondrák precession once.
Neely's common `+1.3960*T` node term represents the moving equinox of the
source convention; applying both that term and the project's frame precession
would count the same frame drift twice.

Mean motion is derived from the transcribed semimajor axis using the standard
Gaussian two-body relation `n = k / a^(3/2)`, with
`k = 0.9856076686 deg/day`. The resulting rates differ from Neely's rounded
printed rates by less than 0.0027 degree per Julian century. The printed rates
remain in the frozen source table and are checked by the provenance gate, so
this choice is visible rather than an undocumented replacement.

For Apollon through Poseidon, `e=i=0`. Perihelion and node are geometrically
undefined in that limit, and only `M + omega + Omega` affects position. The CSV
preserves Neely's literal convention (`M0=0`, phase in `omega`) rather than
moving the same phase into another angle.

## ID 50: Harrington's nominal Planet X

The source is R. S. Harrington,
[“The Location of Planet X,” *The Astronomical Journal* 96 (1988), 1476–1478](https://articles.adsabs.harvard.edu/pdf/1988AJ.....96.1476H),
with the orbital elements on p. 1478. The reviewed PDF has SHA-256
`c0f74fda659cb5fbd7c4db388b74c3169dc269d481338a1f6c4f163f3ab0d49c`.

Harrington prints perihelion date 1789-08-06, `a=101.2 AU`, `e=0.411`,
`omega=208.5 deg`, `Omega=275.4 deg`, and `i=32.4 deg`. The date converts by
the public Gregorian `julday()` rule to JD 2374696.5, and mean anomaly is zero
at perihelion. The paper does not name a modern machine frame; J2000 is an
explicit LibEphemeris propagation convention, not wording attributed to the
author.

## ID 51: Le Verrier's prediction

The primary short account is U.-J. Le Verrier, “Sur la planète qui produit les
anomalies observées dans le mouvement d'Uranus,” *Comptes rendus* 23
(31 August 1846), 428–438, especially
[p. 432 in the public facsimile](https://fr.wikisource.org/wiki/Page:Comptes_rendus_hebdomadaires_des_s%C3%A9ances_de_l%E2%80%99Acad%C3%A9mie_des_sciences,_tome_023,_1846.djvu/436).
It gives `a=36.154 AU`, `e=0.10761`, longitude of perihelion
`284 deg 45 arcmin`, and mean longitude `318 deg 47 arcmin` at 1847-01-01.
Le Verrier's longer *Recherches* table on pp. 234–236 provides the retained
precision `a=36.1539 AU` and `M=34 deg 1 arcmin 56 arcsec`.

The prediction is explicitly developed first in the ecliptic plane, so
`i=Omega=0` are source model assumptions, not guessed replacements. The
runtime epoch and equinox are 1847-01-01 (JD 2395662.5), and the sexagesimal
angles are converted by division by 60 and 3600.

## ID 52: Adams's second-distance prediction

The source is J. C. Adams, *An Explanation of the Observed Irregularities in
the Motion of Uranus* (read 13 November 1846; printed in the Appendix to the
*Nautical Almanac* for 1851), sections 47–48, printed pp. 25–26. A public scan
is available from the
[Internet Archive record](https://archive.org/details/explanationofobs00adam).
The reviewed PDF has SHA-256
`451b7668369911d6241128616a4dc6a108320e84a0b34bf9f4fbbfd6c6e873ec`.

Adams prints, for 1846-10-06, mean longitude `323 deg 02 arcmin`, longitude of
perihelion `299 deg 11 arcmin`, and `e=0.120615`. Direct subtraction gives
`M=23 deg 51 arcmin=23.85 deg`. His second hypothesis prints the distance
ratio `a(Uranus)/a(body)=0.515`. Benjamin A. Gould's 1850 Smithsonian
[*Report on the History of the Discovery of Neptune*](https://commons.wikimedia.org/wiki/File:Report_on_the_history_of_the_discovery_of_Neptune_(IA_reportonhistoryo00goulrich).pdf),
p. 30, explicitly reports the corresponding second mean distance as 37.25 AU
and cites Adams sections 47–48. This near-contemporary cross-check makes the AU
conversion auditable without relying on a later compatibility table. Adams's
section 60 explains that latitude observations did not determine a satisfactory
node or inclination, so the published prediction remains in the ecliptic plane.

The row uses epoch/equinox 1846-10-06 (JD 2395575.5). In particular,
`299 deg 11 arcmin = 299.183333... deg`; treating the arcminutes as decimal
hundredths would give 299.11 and is a scientific transcription error.

## ID 53: Lowell's preferred 1915 solution

The source is Percival Lowell, *Memoir on a Trans-Neptunian Planet* (1915),
pp. 5, 9, and 105, available through the
[HathiTrust catalogue record](https://catalog.hathitrust.org/Record/009026348).
Page 5 defines epsilon as mean longitude at the epoch and varpi as longitude
of perihelion. Page 9 identifies the calculation epoch as 1850.0. The preferred
solution on p. 105 prints `epsilon=22.1 deg`, `a=43.0 AU`, `e=0.202`, and
`varpi=203.8 deg`, so `M=(22.1-203.8) mod 360=178.3 deg`.

Lowell states that the latitude perturbations did not determine a trustworthy
orbital plane. His analogy-based expectation of about 10 degrees inclination
is not a determined element. The reproducible source model is therefore
planar (`i=Omega=0`) in the J1850 ecliptic, with calendar epoch 1850-01-01
(JD 2396757.5). Only page screenshots were needed for this transcription; no
HathiTrust scan is redistributed with LibEphemeris.

## ID 56: White Moon / Selena

### Published rule and page-level checkpoints

The reviewed source is Felix Velichko and Max Larin, *Lilith, Selena,
Proserpina: Articles and Ephemerides 1800–2050* (Mir Uranii, Moscow, 2007),
ISBN 5-900191-12-5, available through the
[public page viewer](https://djvu.online/file/lBJKQiKf9n4Gd). The source
publication retains its rights; LibEphemeris does not bundle its scan or copy
its ephemeris table.

The following page-level facts determine and independently test the convention:

1. p. 17 defines Selena as a symbolic zodiacal point with a seven-year cycle
   that may be treated as uniform motion;
2. p. 18 describes a distinct seven-years-minus-fifteen-days cycle and states
   that this alternative is not used in the printed ephemerides; and
3. p. 20 prints January 1800 as 6 degrees 08 arcminutes Taurus;
4. p. 29 prints March 1879 as 27 degrees 29 arcminutes Leo; and
5. p. 45 prints January 2000 as 2 degrees 18 arcminutes Sagittarius,
   December 2000 as 19 degrees 28 arcminutes Capricorn, and January 2007 as
   2 degrees 22 arcminutes Sagittarius.

Only those defining facts are transcribed. No compatibility-engine output and
no bulk extraction from the 1800–2050 table is stored or fitted.

### Continuous LibEphemeris realization

The book's month rows do not state a sub-day instant or dynamical time scale.
That ambiguity is resolved by an explicit project convention: each row is
interpreted at 00:00 TT on the first civil day of the month. Thus the January
1800 and January 2000 source epochs are JD 2378496.5 TT and JD 2451544.5 TT.
Their printed phases convert by ordinary zodiac arithmetic to

```text
lambda_1800 = 30 deg + 6 deg + 8/60 deg = 36.133333... deg
lambda_2000 = 240 deg + 2 deg + 18/60 deg = 242.3 deg.
```

Using the same nominal point within January at both endpoints makes their
73,048-day separation independent of the otherwise unknown within-month sample
instant. The source's seven-year statement selects the unwrap: 200 years
contain 28 complete turns plus the printed residual arc. Twenty-nine turns
would imply a roughly 6.76-year cycle and is inconsistent with that statement.
The entire rate derivation is therefore:

```text
delta_t = 2451544.5 - 2378496.5 = 73048 days
delta_lambda = 28*360 + (242 deg 18 arcmin - 36 deg 08 arcmin)
             = 10286 deg 10 arcmin
n = delta_lambda / delta_t = 0.14081380279633482 deg/day
P = 360 / n = 2556.5675583712755 days
lambda(JD) = (242.3 + n * (JD - 2451544.5)) mod 360.
```

January 1800 and January 2000 are the two defining endpoints, not fitted
samples from a compatibility engine. Three unused source rows then test the
result independently: March 1879 differs by about `0.0312` arcminute,
December 2000 by about `0.3574` arcminute, and January 2007 by about `0.3464`
arcminute. Every residual is within the half-arcminute rounding bound of a
table printed to the nearest arcminute. This long-baseline reconstruction is
also why LibEphemeris does not impose a modern `365.25`-day year that the
authors did not specify.

The source describes motion through the zodiac but supplies neither latitude
nor a physical radius. LibEphemeris consequently returns zero ecliptic
latitude as the explicit longitude-only convention. The finite positive
distance required by the six-component compatibility API is a display value,
not a physical claim: it is the circular Kepler radius for the
source-derived `P=2556.5675583712755 days`
under the exact nominal terrestrial mass parameter
`3.986004e14 m^3/s^2` from
[IAU 2015 Resolution B3](https://www.iau.org/static/resolutions/IAU2015_English.pdf),
converted with the exact `149597870700 m` astronomical unit from
[IAU 2012 Resolution B2](https://www.iau.org/static/resolutions/IAU2012_English.pdf).
The result, about `0.05279359825 AU`, exists solely to preserve the return
tuple's shape and finite-distance invariant.

## Shared independent orbital propagation

Every supported orbital row is propagated by project-authored, source-neutral
two-body mathematics:

1. advance mean anomaly with the row's source-derived mean motion;
2. solve `M = E - e sin(E)` by Newton iteration;
3. form the perifocal Cartesian vector;
4. rotate it by argument of perihelion, inclination, and ascending node;
5. precess the source ecliptic/equinox to J2000 with the project's published
   long-term precession implementation; and
6. compute speed as the centered numerical derivative of that same position
   model.

No correction term is fitted to a compatibility engine. No output from a
reference engine is stored as an element, fixture, or coefficient.

## Precomputed companion for IDs 40–47

The eight Hamburg-school bodies also ship as the `uranians` LEB2 companion
(`{tier}_uranians.leb2`). Its Chebyshev coefficients are fitted exclusively by
sampling the runtime propagation above (`generate_body_helio` in
`scripts/generate_leb.py` calls `_calc_uranian_planet_raw`); two-body orbits
are smooth enough that 10-year segments at degree 8 reproduce the model to
about `1e-12` degrees, and the compressed channels keep a `1e-12`
native-component target (the channels store degrees, so the default
AU-calibrated target would be too loose).

The runtime registry remains the source of truth:

- the partial `ephemeris_{tier}_uranians.leb` never merges into a main LEB1
  file, and conversion authenticates the exact {40..47} inventory;
- at runtime the companion is attached and consulted only when the file
  byte-matches its `DATA_FILES` SHA-256 pin, both at composite discovery and
  again when a calculation sources a fictitious body;
- dates outside the companion range, absent readers, or hash mismatches all
  fall back to the analytical model with identical public semantics — the
  Uranian branch keeps its own transform and centered-derivative chain and
  only swaps where raw body/Earth positions come from;
- `fast_calc` continues to reject IDs 40–58 from persisted channels, and IDs
  48–58 have no persisted representation at all.

## Custom orbits

The neutral orbital-elements parser and propagator remain available for
caller-supplied data. Callers are responsible for choosing data they are
entitled to use:

```python
import libephemeris as ephem

orbits = ephem.parse_orbital_elements("/path/to/orbits.csv")
body = ephem.get_orbital_body_by_name(orbits, "MyBody")
if body is not None:
    position = ephem.calc_orbital_position(body, 2451545.0)
```

## Verification

```bash
uv run python scripts/check_hypothetical_provenance.py --explain
uv run python scripts/check_provenance.py
uv run pytest tests/test_hypothetical.py tests/test_cov100_hypothetical.py \
  tests/test_orbital_elements_parser.py -q
```

The dedicated gate pins the 12-row CSV, literal source fields, reviewed-scan
hashes where a scan may be retained, the page-located Selena derivation,
source-rate consistency, registry coverage, finite behavior for all supported
IDs, `NaN` compatibility sentinels, and `UnknownBodyError` for all six
unsupported IDs. Reference-API calls may be used only as ephemeral
pass/fail checks; their outputs must never become rows, fixtures, fitted
coefficients, or generated artifacts.
