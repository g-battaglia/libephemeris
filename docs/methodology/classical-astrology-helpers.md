# Classical and Historical Astrology Helpers

This page is the source boundary for `libephemeris.arabic_parts` and
`libephemeris.contrib`. These helpers return conventional astrological
quantities, not observations or physical ephemeris states. That distinction
does not relax the provenance rule: a historical attribution needs a precise
text locator, an elementary identity must be shown as mathematics, and a value
without a verified historical locator must be labelled as a LibEphemeris
convention rather than described as ancient or standard.

The code uses no compatibility-output table as a source. Compatibility affects
names, signatures, tuple shapes, and error behavior only.

## Evidence labels

| Label | Meaning |
|---|---|
| **Published definition** | The implemented rule has a named, independently accessible text and section/verse locator. |
| **Mathematical identity** | The result follows completely from the equation printed in the code; no empirical or historical coefficient is involved. |
| **Project convention** | LibEphemeris preserves a documented API behavior, grouping, threshold, or normalization for which no sufficiently precise primary locator has been established. No claim of historical authority is made. |
| **Delegation** | The helper merely calls another source-mapped LibEphemeris calculation. |

## Lots in `arabic_parts.py`

All longitudes below are ecliptic degrees and are normalized modulo 360. In
the notation `ASC + A - B`, the arc from `B` to `A` is projected from the
Ascendant.

| Function | Exact implemented rule | Status and source |
|---|---|---|
| `calc_arabic_part_of_fortune` | day: `ASC + Moon - Sun`; night: `ASC + Sun - Moon` | **Published definition.** Al-Biruni, *The Book of Instruction in the Elements of the Art of Astrology* (1029), R. Ramsay Wright translation (1934), section/verse 476, “Fortune or Lunar horoscope.” |
| `calc_arabic_part_of_spirit` | day: `ASC + Sun - Moon`; night: `ASC + Moon - Sun` | **Published definition.** Same section, “Daemon and religion [Spirit].” |
| `calc_arabic_part_of_love` | `ASC + Venus - Sun`, unchanged by sect | **Project convention.** This is the library's explicitly retained convenience formula. It must not be attributed to al-Biruni: his section 476 “Friendship and love” instead uses Spirit and Fortune and reverses them by sect. |
| `calc_arabic_part_of_faith` | `ASC + Mercury - Moon`, unchanged by sect | **Project convention.** No page-level primary locator has been established for this exact name/formula pairing. It must not be described as al-Biruni's religion/Spirit lot. |
| `calc_all_arabic_parts` | evaluates the four rules above after validating five required inputs | **Project orchestration.** It introduces no additional lot formula. |

The public digital facsimile identifies Wright's edition and makes the source
independently inspectable:

- Al-Biruni, [*Book of Instruction* digital facsimile](https://www.skyscript.co.uk/albiruni_elements.html), source description and section navigation;
- [section 476 transcription and tabulation](https://skyscript.co.uk/alparts.html), where the Fortune and Spirit day/night expressions are shown explicitly.

These links identify the evidence; their surrounding website prose is not
copied into LibEphemeris.

### Sect determination

Sect is geometrically “Sun on or above the true horizon.” When Julian day and
observer coordinates are available at absolute latitudes above 60 degrees,
`is_day_chart` delegates to the source-mapped ecliptic-to-horizontal pipeline
and tests unrefracted altitude `>= 0`. The strict `60 degrees` switch is a
**project convention** chosen to retain the established low-latitude API path;
it is not an astronomical constant. Without complete observer inputs, the
fallback uses the open longitude arc from Ascendant to Descendant and therefore
cannot resolve every high-latitude horizon geometry. Sunrise and sunset
boundaries count as day by explicit project policy.

## Zodiac and Jyotisha helpers in `contrib`

### Published partitions and dignity rules

The following public transcription of *Brihat Parashara Hora Shastra* is used
only as a locator for the stated definitions, not as copied implementation
prose: [chapter 3, Sanskrit verses and English explanation](https://parashara.net/index.php?page=chapter3).

| Code group | Implemented definition | Locator |
|---|---|---|
| `long2rasi` and twelve sign IDs | `floor((longitude mod 360) / 30)` | Chapter 3, verse 5: the circle is divided equally into twelve rashis beginning with Mesha. |
| `long2nakshatra` and pada | 27 equal sectors of `360/27` degrees; each sector quartered | Chapter 3, verse 4 identifies the 27 nakshatras beginning with Asvini. Equal numerical subdivision and zero-based indexing are project representation choices. |
| `_EXALTATION` | Sun 10 Aries; Moon 3 Taurus; Mars 28 Capricorn; Mercury 15 Virgo; Jupiter 5 Cancer; Venus 27 Pisces; Saturn 20 Libra | Chapter 3, verses 49-50 state the signs and the degree sequence; debilitation at the same degree in the opposite sign is stated in verse 50. |
| `_NAISARGIKA` | expanded friend/neutral/enemy matrix | Chapter 3, verse 55 supplies the natural-friend rule. The checked-in matrix is its explicit finite expansion for the seven classical grahas. |
| `tatkalika_relation` | temporal friends at ordinal houses 2, 3, 4, 10, 11, and 12; otherwise enemies; same position returned as neutral by API policy | Chapter 3, verse 56 supplies the six friend houses. Zero-based offsets `{1,2,3,9,10,11}` are the exact program representation; same-house neutrality is a project API convention. |

`long2navamsa` implements the explicitly printed four-element starting-sign
mapping in its function docstring. Until a page/verse locator for that exact
mapping is added to this record, it is classified as a **project convention**,
not as a primary-source transcription.

The `_RASI_LORDS` tuple is the conventional seven-graha rulership table used by
the public API. Its exact values are visible beside every sign in the source.
Because this audit has not established one page-level primary locator for the
whole tuple, the registry treats the table as a **documented project
convention**. This is deliberately narrower than claiming that an unspecified
“Vedic tradition” is a source.

### Mathematical and project-defined helpers

| Helper | Provenance boundary |
|---|---|
| `rasinorm`, `rasi_dif`, `rasi_dif2` | **Mathematical identities:** integer modular arithmetic over twelve signs. |
| `antiscion` | **Mathematical identity:** reflection `180 degrees - longitude`; contra-antiscion is its antipode. Historical naming does not supply a numerical coefficient. |
| `degsplit`, `signtostr`, `get_nakshatra_name`, coordinate/text helpers | **Project representation:** normalization, sexagesimal carry, names, parsing, and formatting are fully stated in code. |
| aspect constants and `match_*`/`next_*` | **Mathematical/project convention:** named angles are literal degree separations or exact fractions of 360; orb, direction, scan, and error policy come from arguments and documented code. |
| `ochchabala` | Exaltation/debilitation points are published as above; the linear map from 0 to 60 is the function's explicit convention. |
| `residential_strength` | **Project convention:** `1 - degree_within_sign / 30`. LibEphemeris makes no claim that this normalized linear helper is a classical strength doctrine. |
| `raman_houses` | **Delegation/project naming:** calls the source-mapped Sripati code `S`; it contains no separate Raman coefficients or geometry. A primary locator for the alias name is not asserted. |
| `saturn_4_stars` | **Delegation/project grouping:** requests Aldebaran, Regulus, Antares, and Fomalhaut from the public-catalogue fixed-star pipeline. The function name and grouping are compatibility conventions; the code does not claim a physical Saturn relationship or a Vedic origin. |
| time/JD helpers | **Delegation or explicit civil arithmetic:** Julian-day conversion delegates to the registered time module; parsing/formatting is project code. |
| database/atlas/time-zone stubs | **Public-API boundary:** no database model is shipped; every stub raises `NotImplementedError` instead of inventing data. |

## Maintenance rule

A change to one of these helpers must update this table and its module
docstring in the same commit. A new historical attribution needs an edition,
section/page/verse, and exact mapping to the implemented fields. If that
evidence is unavailable, retain the behavior only as a clearly named project
convention or fail closed; never upgrade “commonly reported” into “published”
without evidence.
