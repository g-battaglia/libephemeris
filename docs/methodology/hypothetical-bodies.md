# Hypothetical bodies

LibEphemeris calculates every historical compatibility ID 40--58. These are
mathematical or astrological conventions, not claims that the proposed objects
exist.

The bundled CSV carries the fixed elements for each body.
`scripts/check_hypothetical_provenance.py` pins the CSV by SHA-256 and
independently checks runtime values, coverage, and that all IDs return finite
states.

## Per-body sources

| IDs | Bodies | Source note |
|---:|---|---|
| 40--43 | Cupido, Hades, Zeus, Kronos | Alfred Witte's original articles; refined elements from the Neely table. |
| 44--47 | Apollon, Admetos, Vulkanus, Poseidon | Friedrich Sieggruen, “Nachrichten der Hamburger Schule. I. Nova” (1929), and *Ephemeride der transneptunischen Planeten Apollon, Admetos, Vulkanus und Poseidon* (1937). |
| 48 | Transpluto / Isis | M. Strubell, “Die Sterne” 3/1952, p. 70ff. The published 1772.76 orbit is propagated and precessed to J2000 by LibEphemeris. |
| 49 | Nibiru | Built-in compatibility convention. |
| 50 | Harrington | R. S. Harrington, “The Location of Planet X,” *Astronomical Journal* 96 (1988), p. 1478. |
| 51 | Le Verrier Neptune | Le Verrier's 1846 *Comptes rendus* papers; elements via the Hoyt (1980) compilation. |
| 52 | Adams Neptune | J. C. Adams, “Corrected Elements of Neptune,” *MNRAS* 7 (1847), pp. 244--245; elements via the Hoyt (1980) compilation. |
| 53 | Lowell Planet X | P. Lowell, *Memoir on a Trans-Neptunian Planet* (1915). |
| 54 | Pickering Planet O/X | W. H. Pickering's historical predictions. |
| 55 | Vulcan | Time-dependent L. H. Weston convention. |
| 56 | White Moon | Russian-school seven-year circular Selena convention. |
| 57 | Proserpina | Circular V. Abramov convention. |
| 58 | Waldemath | Waldemath's 1898 claim; built-in compatibility convention for the detailed elements. |

## Propagation

Fixed heliocentric rows use an independent two-body implementation: Gaussian
mean motion, Newton iteration for Kepler's equation, perifocal-to-ecliptic
rotation, and LibEphemeris' Vondrak precession. Vulcan, White Moon, and Waldemath
retain their time-dependent compatibility conventions. All reported rates are
derivatives of the corresponding local model, not stored observations.

## Primary-source bibliography

- Alfred Witte, [“Der erste Transneptunplanet Cupido?” (1923)](https://witte-verlag.com/pdf/source/1923-Der_erste_Transneptunplanet_Cupido.pdf).
- Alfred Witte, [“Der synodische Lauf des Planeten Cupido”](https://astrax.de/pdf/22_Synodischer%20Lauf%20des%20Planeten%20Cupido.pdf).
- Alfred Witte, [“Wahrscheinlicher Lauf des 2. Transneptun-Planeten Hades”](https://astrax.de/pdf/39_Wahrscheinlicher%20Lauf%20des%202.%20Transneptun-Planeten%20Hades.pdf).
- Alfred Witte, [“Der 4. Transneptun-Planet Kronos”](https://astrax.de/pdf/27_Der%204.%20Transneptun-Planet%20Kronos.pdf).
- Urbain Le Verrier, *Comptes rendus* 22 (1846), pp. 907--918, and 23 (1846), pp. 428--438; [primary scans indexed by Bibnum](https://www.bibnum.education.fr/physique/astronomie/la-decouverte-de-neptune-1846).
- John Couch Adams, [“Corrected Elements of Neptune,” *MNRAS* 7 (1847), 244--245](https://doi.org/10.1093/mnras/7.13.244).
- Percival Lowell, [*Memoir on a Trans-Neptunian Planet* (1915)](https://books.google.com/books?id=d9dLAAAAYAAJ).
- R. S. Harrington, [“The Location of Planet X,” *AJ* 96 (1988), 1476--1478](https://articles.adsabs.harvard.edu/pdf/1988AJ.....96.1476H).

## Verification

```bash
uv run python scripts/check_hypothetical_provenance.py
uv run python scripts/check_provenance.py
pytest tests/test_hypothetical.py -q
```

Reference-API calls may be used only as ephemeral pass/fail checks.
Their outputs must never be persisted as rows, fixtures, coefficients, or
generated artifacts.
