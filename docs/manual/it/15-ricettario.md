# Capitolo 15 — Ricettario: calcoli pratici dalla A alla Z

## Cosa imparerai

Questo capitolo è una raccolta di **ricette complete**, pronte da copiare e incollare. Ogni ricetta risolve un problema pratico di uso comune: dal tema natale ai transiti, dalle retrogradazioni alle eclissi. Non c'è teoria da leggere — solo codice funzionante con il suo output.

---

## Ricetta 1 — Tema natale completo

**Problema**: data, ora e luogo di nascita, calcolare tutte le posizioni planetarie e le case.

**A cosa serve**: il tema natale è la base di ogni analisi astrologica. Questo codice produce la "fotografia del cielo" al momento della nascita: le posizioni dei pianeti nei segni, le dodici case e gli angoli principali (Ascendente e Medio Cielo). Da qui partono tutte le interpretazioni — dalla personalità (posizione del Sole) alle emozioni (Luna), al modo di comunicare (Mercurio), fino alla struttura delle aree di vita (case). Un'applicazione concreta: un software astrologico che genera automaticamente il tema natale partendo dall'anagrafica dell'utente.

```python
import libephemeris as ephem

# Roma, 8 aprile 2024, ore 14:30 UT
jd = ephem.julday(2024, 4, 8, 14.5)

segni = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

def formato_zodiacale(lon):
    segno_idx = int(lon / 30)
    gradi = lon % 30
    d, m, s, sf, si = ephem.split_deg(gradi, 0)
    return f'{d:2d}° {m:02d}\' {s:02d}" {segni[segno_idx]}'

corpi = [
    (ephem.SUN,       "Sole"),
    (ephem.MOON,      "Luna"),
    (ephem.MERCURY,   "Mercurio"),
    (ephem.VENUS,     "Venere"),
    (ephem.MARS,      "Marte"),
    (ephem.JUPITER,   "Giove"),
    (ephem.SATURN,    "Saturno"),
    (ephem.URANUS,    "Urano"),
    (ephem.NEPTUNE,   "Nettuno"),
    (ephem.PLUTO,     "Plutone"),
    (ephem.MEAN_NODE, "Nodo Nord"),
    (ephem.CHIRON,    "Chirone"),
]

print("--- Tema Natale: Roma, 8 aprile 2024, 14:30 ---")
print()
for body_id, nome in corpi:
    pos, _ = ephem.calc_ut(jd, body_id, ephem.FLG_SPEED)
    retro = " R" if pos[3] < 0 else ""
    print(f"{nome:11s}  {formato_zodiacale(pos[0])}{retro}")

# Case Placidus — Roma (lat 41.9, lon 12.5)
cusps, ascmc = ephem.houses(jd, 41.9, 12.5, ord("P"))
print()
print(f"Ascendente:  {formato_zodiacale(ascmc[0])}")
print(f"Medio Cielo: {formato_zodiacale(ascmc[1])}")
print()
for i in range(12):
    print(f"Casa {i+1:2d}:     {formato_zodiacale(cusps[i])}")
```

```
--- Tema Natale: Roma, 8 aprile 2024, 14:30 ---

Sole         19° 14' 34" Ari
Luna         16° 59' 41" Ari
Mercurio     24° 53' 58" Ari R
Venere        4° 14' 49" Ari
Marte        12° 55' 36" Psc
Giove        19° 00' 36" Tau
Saturno      14° 26' 16" Psc
Urano        21° 09' 46" Tau
Nettuno      28° 11' 03" Psc
Plutone       1° 57' 56" Aqr
Nodo Nord    15° 39' 19" Ari R
Chirone      19° 23' 45" Ari

Ascendente:  12° 15' 00" Vir
Medio Cielo:  9° 02' 28" Gem

Casa  1:     12° 15' 00" Vir
Casa  2:      6° 21' 02" Lib
Casa  3:      5° 32' 08" Sco
Casa  4:      9° 02' 28" Sgr
Casa  5:     13° 32' 25" Cap
Casa  6:     15° 02' 03" Aqr
Casa  7:     12° 15' 00" Psc
Casa  8:      6° 21' 02" Ari
Casa  9:      5° 32' 08" Tau
Casa 10:      9° 02' 28" Gem
Casa 11:     13° 32' 25" Cnc
Casa 12:     15° 02' 03" Leo
```

---

## Ricetta 2 — Transiti del giorno

**Problema**: data una carta natale, trovare quali pianeti di transito sono in aspetto con i pianeti natali oggi.

**A cosa serve**: i transiti sono il metodo più diffuso di previsione astrologica. Il principio è semplice: i pianeti "di oggi" formano angoli significativi (congiunzione, quadratura, trigono…) con le posizioni della carta natale, attivando i temi rappresentati da quei pianeti. Per esempio, Saturno in quadratura al Sole natale indica un periodo di sfide e ristrutturazione della propria identità — tipicamente dura qualche settimana. Un'applicazione concreta: un'app che ogni mattina mostra all'utente i transiti attivi sul suo tema, con l'orbe (quanto manca all'aspetto esatto).

```python
import libephemeris as ephem

# Posizioni natali: nato il 15 marzo 1990, 10:00 UT
jd_natale = ephem.julday(1990, 3, 15, 10.0)
# Data del transito: 8 aprile 2024
jd_transito = ephem.julday(2024, 4, 8, 12.0)

pianeti = [
    (ephem.SUN, "Sole"), (ephem.MOON, "Luna"),
    (ephem.MERCURY, "Mercurio"), (ephem.VENUS, "Venere"),
    (ephem.MARS, "Marte"), (ephem.JUPITER, "Giove"),
    (ephem.SATURN, "Saturno"),
]

# Calcola posizioni natali e di transito
natali = {}
for body_id, nome in pianeti:
    pos, _ = ephem.calc_ut(jd_natale, body_id, 0)
    natali[nome] = pos[0]

transiti = {}
for body_id, nome in pianeti:
    pos, _ = ephem.calc_ut(jd_transito, body_id, 0)
    transiti[nome] = pos[0]

# Definisci aspetti e orbi
aspetti = {0: "Congiunzione", 60: "Sestile", 90: "Quadratura",
           120: "Trigono", 180: "Opposizione"}
orbi = {0: 8, 60: 6, 90: 7, 120: 7, 180: 8}

print("--- Transiti dell'8 aprile 2024 su carta del 15/3/1990 ---")
print()
for t_nome, t_lon in transiti.items():
    for n_nome, n_lon in natali.items():
        diff = abs(ephem.difdeg2n(t_lon, n_lon))
        for asp_gradi, asp_nome in aspetti.items():
            orbe = abs(diff - asp_gradi)
            if orbe <= orbi[asp_gradi]:
                print(f"{t_nome:10s} {asp_nome:14s} {n_nome:10s}"
                      f"  (orbe: {orbe:.1f}°)")
                break
```

```
--- Transiti dell'8 aprile 2024 su carta del 15/3/1990 ---

Sole       Quadratura     Saturno     (orbe: 4.2°)
Mercurio   Quadratura     Saturno     (orbe: 1.6°)
Venere     Sestile        Venere      (orbe: 4.9°)
Venere     Sestile        Marte       (orbe: 1.3°)
Venere     Quadratura     Giove       (orbe: 2.7°)
Marte      Trigono        Luna        (orbe: 4.1°)
Giove      Sestile        Sole        (orbe: 5.6°)
Giove      Sestile        Mercurio    (orbe: 2.1°)
Giove      Trigono        Saturno     (orbe: 4.3°)
Saturno    Trigono        Luna        (orbe: 5.7°)
Saturno    Congiunzione   Mercurio    (orbe: 6.6°)
```

---

## Ricetta 3 — Ingressi del Sole nei segni

**Problema**: trovare la data e l'ora esatta in cui il Sole entra in ciascun segno zodiacale durante l'anno.

**A cosa serve**: gli ingressi del Sole nei segni segnano i cambiamenti di stagione (Ariete = equinozio di primavera, Cancro = solstizio d'estate, ecc.) e sono date fondamentali in astrologia mondiale. Il tema eretto per il momento dell'ingresso in Ariete (la cosiddetta "carta dell'anno") viene usato per fare previsioni sull'intero anno per una nazione. In astronomia, conoscere l'esatto momento dell'equinozio serve per calibrare orologi solari e almanacchi. Un'applicazione concreta: un calendario astrologico che mostra automaticamente le date di cambio segno per ogni anno.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 1, 1, 0.0)

segni = ["Ariete", "Toro", "Gemelli", "Cancro",
         "Leone", "Vergine", "Bilancia", "Scorpione",
         "Sagittario", "Capricorno", "Acquario", "Pesci"]

print("--- Ingressi del Sole nei segni (2024) ---")
print()
for i in range(12):
    target = i * 30.0
    jd_ingresso = ephem.solcross_ut(target, jd)
    anno, mese, giorno, ora = ephem.revjul(jd_ingresso)
    ore = int(ora)
    minuti = int((ora - ore) * 60)
    print(f"{segni[i]:12s}  {giorno:2d}/{mese:02d}/{anno}"
          f"  {ore:02d}:{minuti:02d} UT")
    jd = jd_ingresso + 1
```

```
--- Ingressi del Sole nei segni (2024) ---

Ariete        20/03/2024  03:06 UT
Toro          19/04/2024  13:59 UT
Gemelli       20/05/2024  12:59 UT
Cancro        20/06/2024  20:50 UT
Leone         22/07/2024  07:44 UT
Vergine       22/08/2024  14:55 UT
Bilancia      22/09/2024  12:43 UT
Scorpione     22/10/2024  22:14 UT
Sagittario    21/11/2024  19:56 UT
Capricorno    21/12/2024  09:20 UT
Acquario      19/01/2025  20:00 UT
Pesci         18/02/2025  10:06 UT
```

---

## Ricetta 4 — Retrogradazione di Mercurio

**Problema**: trovare tutte le stazioni retrograde e dirette di Mercurio nell'anno, e verificare se è retrogrado in una data specifica.

**A cosa serve**: la retrogradazione di Mercurio è uno degli eventi astrologici più popolari — e temuti. Quando Mercurio è retrogrado (cioè sembra muoversi all'indietro visto dalla Terra), la tradizione astrologica associa il periodo a malintesi nella comunicazione, problemi con contratti e tecnologia, ritardi nei viaggi. In pratica, molte persone evitano di firmare contratti o acquistare elettronica in quei periodi. Questo codice trova le date esatte di inizio (stazione retrograda) e fine (stazione diretta) di ogni retrogradazione, con il grado zodiacale preciso. Un'applicazione concreta: un'agenda digitale che evidenzia automaticamente i periodi di Mercurio retrogrado.

```python
import libephemeris as ephem
from libephemeris.crossing import find_station_ut, is_retrograde

jd = ephem.julday(2024, 1, 1, 0.0)
jd_fine = ephem.julday(2025, 1, 1, 0.0)

segni = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

print("--- Retrogradazioni di Mercurio nel 2024 ---")
print()
while jd < jd_fine:
    jd_station, tipo = find_station_ut(
        ephem.MERCURY, jd, "any"
    )
    if jd_station >= jd_fine:
        break
    anno, mese, giorno, ora = ephem.revjul(jd_station)
    ore = int(ora)
    minuti = int((ora - ore) * 60)
    pos, _ = ephem.calc_ut(jd_station, ephem.MERCURY, 0)
    segno = segni[int(pos[0] / 30)]
    gradi = pos[0] % 30
    label = ("Stazione Retrograda" if tipo == "SR"
             else "Stazione Diretta  ")
    print(f"{label}  {giorno:2d}/{mese:02d}/{anno}"
          f"  {ore:02d}:{minuti:02d} UT  a {gradi:5.1f}° {segno}")
    jd = jd_station + 1

# Verifica puntuale
print()
retro = is_retrograde(ephem.MERCURY,
                      ephem.julday(2024, 4, 8, 12.0))
print(f"Mercurio retrogrado l'8/4/2024? {retro}")
```

```
--- Retrogradazioni di Mercurio nel 2024 ---

Stazione Diretta     2/01/2024  03:07 UT  a  22.2° Sgr
Stazione Retrograda   1/04/2024  22:14 UT  a  27.2° Ari
Stazione Diretta    25/04/2024  12:54 UT  a  16.0° Ari
Stazione Retrograda   5/08/2024  04:56 UT  a   4.1° Vir
Stazione Diretta    28/08/2024  21:13 UT  a  21.4° Leo
Stazione Retrograda  26/11/2024  02:42 UT  a  22.7° Sgr
Stazione Diretta    15/12/2024  20:56 UT  a   6.4° Sgr

Mercurio retrogrado l'8/4/2024? True
```

> **Nota**: `find_station_ut` e `is_retrograde` si importano da `libephemeris.crossing` perché non sono esposti a livello di modulo principale.

---

## Ricetta 5 — Luna Nuova e Luna Piena

**Problema**: trovare tutte le Lune Nuove e Piene dell'anno con data, ora e segno zodiacale.

**A cosa serve**: le fasi lunari scandiscono il ritmo della vita quotidiana da millenni. La Luna Nuova (congiunzione Sole-Luna) è tradizionalmente il momento di seminare nuovi progetti, mentre la Luna Piena (opposizione) è il momento del raccolto e della massima energia — ma anche di tensione emotiva. In agricoltura biodinamica, le semine seguono ancora il calendario lunare. In astrologia, il segno in cui cade la Luna Nuova indica il "tema" del mese: una Luna Nuova in Ariete invita ad agire con coraggio, una in Pesci a lasciar andare. Un'applicazione concreta: un calendario lunare che mostra le fasi con i segni, utile per giardinieri, pescatori, o chi segue ritmi lunari.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 1, 1, 0.0)
jd_fine = ephem.julday(2025, 1, 1, 0.0)

segni = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

print("--- Lune Nuove e Piene del 2024 ---")
print()
while jd < jd_fine:
    # Longitudine del Sole adesso
    pos_sole, _ = ephem.calc_ut(jd, ephem.SUN, 0)
    lon_sole = pos_sole[0]

    # Luna Nuova: Luna congiunge il Sole
    jd_nuova = ephem.mooncross_ut(lon_sole, jd)
    if jd_nuova < jd_fine:
        anno, mese, giorno, ora = ephem.revjul(jd_nuova)
        ore = int(ora)
        minuti = int((ora - ore) * 60)
        pos_luna, _ = ephem.calc_ut(jd_nuova, ephem.MOON, 0)
        segno = segni[int(pos_luna[0] / 30)]
        gradi = pos_luna[0] % 30
        print(f"Luna Nuova   {giorno:2d}/{mese:02d}/{anno}"
              f"  {ore:02d}:{minuti:02d} UT"
              f"  in {gradi:5.1f}° {segno}")

    # Luna Piena: Luna opposta al Sole (180°)
    target_piena = (lon_sole + 180.0) % 360.0
    jd_piena = ephem.mooncross_ut(target_piena, jd)
    if jd_piena < jd_fine and jd_piena > jd:
        anno, mese, giorno, ora = ephem.revjul(jd_piena)
        ore = int(ora)
        minuti = int((ora - ore) * 60)
        pos_luna, _ = ephem.calc_ut(jd_piena, ephem.MOON, 0)
        segno = segni[int(pos_luna[0] / 30)]
        gradi = pos_luna[0] % 30
        print(f"Luna Piena   {giorno:2d}/{mese:02d}/{anno}"
              f"  {ore:02d}:{minuti:02d} UT"
              f"  in {gradi:5.1f}° {segno}")

    print()
    # Avanza di ~20 giorni per trovare la prossima coppia
    jd = min(jd_nuova, jd_piena) + 20
    if jd >= jd_fine:
        break
```

```
--- Lune Nuove e Piene del 2024 ---

Luna Nuova   10/01/2024  18:20 UT  in  10.0° Cap
Luna Piena   23/01/2024  16:58 UT  in  10.0° Cnc

Luna Nuova    9/02/2024  06:34 UT  in  10.3° Aqr
Luna Piena   22/02/2024  10:11 UT  in  10.3° Leo

Luna Nuova    9/03/2024  17:05 UT  in  10.2° Psc
Luna Piena   23/03/2024  04:18 UT  in  10.2° Vir

Luna Nuova    8/04/2024  02:31 UT  in   9.5° Ari
Luna Piena   21/04/2024  22:20 UT  in   9.5° Lib

Luna Nuova    7/05/2024  11:18 UT  in   8.3° Tau
Luna Piena   21/05/2024  14:50 UT  in   8.3° Sco

Luna Nuova    5/06/2024  19:59 UT  in   6.7° Gem
Luna Piena   20/06/2024  04:55 UT  in   6.7° Sgr

Luna Nuova    5/07/2024  05:23 UT  in   4.7° Cnc
Luna Piena   19/07/2024  16:32 UT  in   4.7° Cap

Luna Nuova    3/08/2024  16:25 UT  in   2.8° Leo
Luna Piena   18/08/2024  02:22 UT  in   2.8° Aqr

Luna Nuova    2/09/2024  05:50 UT  in   1.0° Vir
Luna Piena   16/09/2024  11:18 UT  in   1.0° Psc

Luna Nuova    1/10/2024  21:45 UT  in  29.7° Vir
Luna Piena   15/10/2024  20:07 UT  in  29.7° Psc

Luna Nuova   31/10/2024  15:26 UT  in  29.0° Lib
Luna Piena   14/11/2024  05:21 UT  in  29.0° Ari

Luna Nuova   30/11/2024  09:32 UT  in  28.8° Sco
Luna Piena   13/12/2024  15:23 UT  in  28.8° Tau

Luna Nuova   30/12/2024  02:45 UT  in  29.0° Sgr
```

---

## Ricetta 6 — Prossima eclissi solare dalla mia città

**Problema**: trovare le prossime eclissi solari visibili da una località specifica, con orari dei contatti e magnitudine.

**A cosa serve**: le eclissi solari sono eventi spettacolari ma rari per una data località — possono passare decenni tra un'eclissi totale e l'altra. Sapere in anticipo quando e come sarà visibile un'eclissi dalla propria città è essenziale per organizzare l'osservazione (filtri solari, viaggio, ferie). La magnitudine indica quanto disco solare sarà coperto: 1.0 = totale, 0.5 = metà, 0.01 = appena percettibile. Un'applicazione concreta: un sito web che, data la posizione dell'utente, mostra la prossima eclissi visibile con una mappa dei contatti.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 1, 1, 0.0)
lat, lon, alt = 41.9, 12.5, 0.0  # Roma

print("--- Prossime eclissi solari visibili da Roma ---")
print()
trovate = 0
for _ in range(20):  # cerca fino a 20 eclissi avanti
    ecl_type, times, attr = ephem.sol_eclipse_when_loc(
        jd, (lon, lat, alt)
    )
    if ecl_type > 0 and attr[0] > 0.01:
        jd_max = times[0]
        anno, mese, giorno, ora = ephem.revjul(jd_max)
        ore = int(ora)
        minuti = int((ora - ore) * 60)

        tipo = "Parziale"
        if ecl_type & 4: tipo = "Anulare"
        if ecl_type & 1: tipo = "Totale"

        print(f"{tipo:10s}  {giorno:2d}/{mese:02d}/{anno}"
              f"  {ore:02d}:{minuti:02d} UT"
              f"  magnitudine: {attr[0]:.3f}")

        # Orari di inizio e fine
        if times[1] > 0:
            _, _, _, o = ephem.revjul(times[1])
            print(f"  Inizio:  {int(o):02d}:"
                  f"{int((o - int(o)) * 60):02d} UT")
        if times[4] > 0:
            _, _, _, o = ephem.revjul(times[4])
            print(f"  Fine:    {int(o):02d}:"
                  f"{int((o - int(o)) * 60):02d} UT")
        print()
        trovate += 1
        if trovate >= 3:
            break

    jd = times[0] + 30 if times[0] > 0 else jd + 180
```

```
--- Prossime eclissi solari visibili da Roma ---

Parziale    29/03/2025  11:03 UT  magnitudine: 0.073
  Inizio:  10:35 UT
  Fine:    11:31 UT

Parziale     2/08/2027  09:12 UT  magnitudine: 0.786
  Inizio:  08:02 UT
  Fine:    10:26 UT

Parziale     1/06/2030  05:04 UT  magnitudine: 0.816
  Inizio:  04:02 UT
  Fine:    06:12 UT
```

---

## Ricetta 7 — Carta siderale (vedica)

**Problema**: calcolare le posizioni planetarie nello zodiaco siderale con ayanamsha Lahiri, includendo i Nakshatra.

**A cosa serve**: l'astrologia vedica (Jyotish) usa lo zodiaco siderale anziché quello tropicale — le posizioni sono sfasate di circa 24° rispetto a quelle occidentali. I Nakshatra sono 27 "stazioni lunari" di 13°20' ciascuna, ulteriormente divise in 4 Pada; la posizione della Luna nel Nakshatra è fondamentale per calcolare i Dasha (periodi planetari che governano la vita). Per esempio, chi ha la Luna nel Nakshatra di Ashwini vivrà il primo periodo sotto Ketu (7 anni), poi sotto Venere (20 anni), ecc. Un'applicazione concreta: un software di Jyotish che genera la carta siderale, i Nakshatra e il calcolo dei Dasha.

```python
import libephemeris as ephem
from libephemeris.constants import (SIDM_LAHIRI,
                                    FLG_SIDEREAL, FLG_SPEED)

jd = ephem.julday(2024, 4, 8, 14.5)
ephem.set_sid_mode(SIDM_LAHIRI)

nakshatras = [
    "Ashwini", "Bharani", "Krittika", "Rohini",
    "Mrigashira", "Ardra", "Punarvasu", "Pushya",
    "Ashlesha", "Magha", "Purva Phalguni",
    "Uttara Phalguni", "Hasta", "Chitra", "Swati",
    "Vishakha", "Anuradha", "Jyeshtha", "Mula",
    "Purva Ashadha", "Uttara Ashadha", "Shravana",
    "Dhanishta", "Shatabhisha", "Purva Bhadra",
    "Uttara Bhadra", "Revati",
]

segni = ["Mesh", "Vrish", "Mithun", "Kark", "Simha",
         "Kanya", "Tula", "Vrischik", "Dhanu", "Makar",
         "Kumbh", "Meen"]

corpi = [
    (ephem.SUN,       "Surya"),
    (ephem.MOON,      "Chandra"),
    (ephem.MARS,      "Mangal"),
    (ephem.MERCURY,   "Budha"),
    (ephem.JUPITER,   "Guru"),
    (ephem.VENUS,     "Shukra"),
    (ephem.SATURN,    "Shani"),
    (ephem.MEAN_NODE, "Rahu"),
]

ayan = ephem.get_ayanamsa_ut(jd)
print("--- Carta Siderale (Lahiri) ---")
print(f"Ayanamsha: {ayan:.4f}°")
print()

for body_id, nome in corpi:
    pos, _ = ephem.calc_ut(jd, body_id,
                           FLG_SPEED | FLG_SIDEREAL)
    lon = pos[0]
    segno = segni[int(lon / 30)]
    gradi = lon % 30
    nak_idx = int(lon / (360 / 27))
    nak = nakshatras[nak_idx]
    pada = int((lon % (360 / 27)) / (360 / 108)) + 1
    print(f"{nome:8s}  {gradi:5.1f}° {segno:10s}"
          f"  Nakshatra: {nak} (Pada {pada})")

# Ketu = opposto di Rahu
pos_rahu, _ = ephem.calc_ut(jd, ephem.MEAN_NODE,
                            FLG_SPEED | FLG_SIDEREAL)
ketu = (pos_rahu[0] + 180) % 360
segno_k = segni[int(ketu / 30)]
gradi_k = ketu % 30
nak_k = nakshatras[int(ketu / (360 / 27))]
pada_k = int((ketu % (360 / 27)) / (360 / 108)) + 1
print(f"Ketu      {gradi_k:5.1f}° {segno_k:10s}"
      f"  Nakshatra: {nak_k} (Pada {pada_k})")

ephem.set_sid_mode(0)  # ripristina
```

```
--- Carta Siderale (Lahiri) ---
Ayanamsha: 24.1961°

Surya      25.0° Meen        Nakshatra: Revati (Pada 3)
Chandra    22.8° Meen        Nakshatra: Revati (Pada 2)
Mangal     18.7° Kumbh       Nakshatra: Shatabhisha (Pada 4)
Budha       0.7° Mesh        Nakshatra: Ashwini (Pada 1)
Guru       24.8° Mesh        Nakshatra: Bharani (Pada 4)
Shukra     10.1° Meen        Nakshatra: Uttara Bhadra (Pada 3)
Shani      20.2° Kumbh       Nakshatra: Purva Bhadra (Pada 1)
Rahu       21.5° Meen        Nakshatra: Revati (Pada 2)
Ketu       21.5° Kanya       Nakshatra: Hasta (Pada 4)
```

---

## Ricetta 8 — Visibilità dei pianeti stasera

**Problema**: per ogni pianeta visibile a occhio nudo, calcolare altezza, direzione, magnitudine ed elongazione dal Sole un'ora dopo il tramonto.

**A cosa serve**: sapere quali pianeti sono visibili stasera è la domanda più pratica per chi osserva il cielo. L'altezza sull'orizzonte dice se il pianeta è alto abbastanza da vedere (sopra edifici e alberi); l'elongazione dal Sole dice se è immerso nel bagliore del crepuscolo; la magnitudine dice quanto è luminoso (sotto 6 = visibile a occhio nudo, sotto 0 = molto brillante). Per esempio, Giove a magnitudine −2.0 e altezza 15° a est è inconfondibile — è il "punto luminoso" che molti scambiano per un aereo. Un'applicazione concreta: un'app per astrofili che mostra "cosa c'è in cielo stasera" con direzione e luminosità.

```python
import libephemeris as ephem

# Roma, 8 aprile 2024
lat, lon = 41.9, 12.5
jd_noon = ephem.julday(2024, 4, 8, 12.0)

# Trova il tramonto esatto del Sole
_, tret = ephem.rise_trans(
    jd_noon, ephem.SUN, 2, (lon, lat, 0.0)
)
jd_tram = tret[0]
anno, mese, giorno, ora = ephem.revjul(jd_tram)
ore_t = int(ora)
min_t = int((ora - ore_t) * 60)
print(f"Tramonto del Sole: {ore_t:02d}:{min_t:02d} UT")
print()

# Un'ora dopo il tramonto
jd_sera = jd_tram + 1.0 / 24.0
ephem.set_topo(lon, lat, 0)

pos_sole, _ = ephem.calc_ut(jd_sera, ephem.SUN, 0)

pianeti = [
    (ephem.MERCURY, "Mercurio"),
    (ephem.VENUS,   "Venere"),
    (ephem.MARS,    "Marte"),
    (ephem.JUPITER, "Giove"),
    (ephem.SATURN,  "Saturno"),
]

dirs = ["N", "NE", "E", "SE", "S", "SO", "O", "NO"]

print("--- Visibilità pianeti (Roma, 8/4/2024) ---")
print()
for body_id, nome in pianeti:
    pos, _ = ephem.calc_ut(jd_sera, body_id, 0)

    # Altezza e azimut
    az, alt_v, alt_app = ephem.azalt(
        jd_sera, 0, (lon, lat, 0),
        1013.25, 15.0, (pos[0], pos[1], pos[2])
    )

    # Elongazione dal Sole
    elong = abs(ephem.difdeg2n(pos[0], pos_sole[0]))

    # Magnitudine apparente
    pheno = ephem.pheno_ut(jd_sera, body_id, 0)
    mag = pheno[4]

    # Direzione cardinale
    dir_idx = int((az + 22.5) / 45) % 8
    direzione = dirs[dir_idx]

    visibile = ("VISIBILE" if alt_app > 5 and elong > 15
                else "non visibile")
    print(f"{nome:10s}  alt: {alt_app:5.1f}°"
          f"  az: {az:5.1f}° ({direzione:2s})"
          f"  mag: {mag:+.1f}"
          f"  elong: {elong:5.1f}°"
          f"  {visibile}")
```

```
Tramonto del Sole: 17:43 UT

--- Visibilità pianeti (Roma, 8/4/2024) ---

Mercurio    alt:  -5.5°  az: 111.8° (E )  mag: +4.8  elong:   5.4°  non visibile
Venere      alt: -25.8°  az: 116.4° (SE)  mag: -3.9  elong:  15.0°  non visibile
Marte       alt: -44.6°  az: 129.0° (SE)  mag: +1.2  elong:  36.4°  non visibile
Giove       alt:  15.5°  az:  98.8° (E )  mag: -2.0  elong:  29.6°  VISIBILE
Saturno     alt: -43.6°  az: 127.4° (SE)  mag: +1.1  elong:  35.0°  non visibile
```

> In questa sera di aprile 2024, solo Giove è visibile — basso a est, molto luminoso (mag −2.0). Gli altri pianeti sono sotto l'orizzonte o troppo vicini al Sole.

---

## Ricetta 9 — Aspetti tra pianeti

**Problema**: calcolare la distanza angolare tra tutti i pianeti e individuare gli aspetti attivi (congiunzione, sestile, quadratura, trigono, opposizione).

**A cosa serve**: gli aspetti sono le relazioni angolari tra pianeti — il linguaggio fondamentale dell'astrologia. Una congiunzione (0°) fonde le energie dei due pianeti; un trigono (120°) le armonizza; una quadratura (90°) le mette in tensione. Per esempio, Sole congiunto Luna (come l'8 aprile 2024) indica Luna Nuova — un momento di reset emotivo e nuovi inizi. Marte congiunto Saturno indica una fase di frustrazione ma anche di disciplina costruttiva. Sapere se un aspetto è "applicante" (si sta formando) o "separante" (si sta sciogliendo) cambia l'interpretazione: l'applicante è più intenso. Un'applicazione concreta: un modulo di interpretazione automatica che legge gli aspetti e genera un testo descrittivo.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

pianeti = [
    (ephem.SUN,     "Sole"),
    (ephem.MOON,    "Luna"),
    (ephem.MERCURY, "Mercurio"),
    (ephem.VENUS,   "Venere"),
    (ephem.MARS,    "Marte"),
    (ephem.JUPITER, "Giove"),
    (ephem.SATURN,  "Saturno"),
]

# Calcola posizioni e velocità
posizioni = {}
for body_id, nome in pianeti:
    pos, _ = ephem.calc_ut(jd, body_id, ephem.FLG_SPEED)
    posizioni[nome] = (pos[0], pos[3])  # longitudine, velocità

aspetti = {0: "Congiunzione", 60: "Sestile",
           90: "Quadratura", 120: "Trigono",
           180: "Opposizione"}
orbi_max = {0: 8, 60: 5, 90: 6, 120: 6, 180: 8}

print("--- Aspetti tra pianeti (8 aprile 2024) ---")
print()
nomi = list(posizioni.keys())
for i in range(len(nomi)):
    for j in range(i + 1, len(nomi)):
        n1, n2 = nomi[i], nomi[j]
        lon1, vel1 = posizioni[n1]
        lon2, vel2 = posizioni[n2]
        diff = abs(ephem.difdeg2n(lon1, lon2))

        for asp_g, asp_nome in aspetti.items():
            orbe = abs(diff - asp_g)
            if orbe <= orbi_max[asp_g]:
                print(f"{n1:10s} {asp_nome:14s} {n2:10s}"
                      f"  orbe: {orbe:4.1f}°")
                break
```

```
--- Aspetti tra pianeti (8 aprile 2024) ---

Sole       Congiunzione   Luna        orbe:  3.7°
Sole       Congiunzione   Mercurio    orbe:  5.8°
Marte      Congiunzione   Saturno     orbe:  1.6°
Giove      Sestile        Saturno     orbe:  4.6°
```

---

## Ricetta 10 — Rivoluzione solare (Solar Return)

**Problema**: trovare il momento esatto in cui il Sole torna alla stessa longitudine della nascita, e calcolare la carta per quell'istante.

**A cosa serve**: la rivoluzione solare è una delle tecniche previsionali più usate in astrologia. Ogni anno, il Sole torna esattamente al grado in cui si trovava al momento della nascita — quel preciso istante genera una "carta dell'anno" che si sovrappone al tema natale. L'Ascendente della rivoluzione solare indica il "tono" dell'anno: per esempio, un Ascendente in Vergine suggerisce un anno dedicato al lavoro, all'organizzazione e alla salute. La posizione della Luna nella rivoluzione indica il clima emotivo dei prossimi dodici mesi. Un'applicazione concreta: un software che genera automaticamente la rivoluzione solare per ogni compleanno e la confronta con il tema natale.

```python
import libephemeris as ephem

# Nato il 15 marzo 1990, 10:00 UT, Roma
jd_natale = ephem.julday(1990, 3, 15, 10.0)
pos_natale, _ = ephem.calc_ut(jd_natale, ephem.SUN, 0)
lon_natale = pos_natale[0]

# Trova il ritorno solare per il 2024
jd_cerca = ephem.julday(2024, 3, 1, 0.0)
jd_ritorno = ephem.solcross_ut(lon_natale, jd_cerca)

anno, mese, giorno, ora = ephem.revjul(jd_ritorno)
ore = int(ora)
minuti = int((ora - ore) * 60)
secondi = int(((ora - ore) * 60 - minuti) * 60)

segni = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

print("--- Rivoluzione Solare 2024 ---")
print()
print(f"Sole natale: {lon_natale:.6f}°")
print(f"Ritorno:     {giorno:2d}/{mese:02d}/{anno}"
      f"  {ore:02d}:{minuti:02d}:{secondi:02d} UT")
print()

corpi = [
    (ephem.SUN,     "Sole"),
    (ephem.MOON,    "Luna"),
    (ephem.MERCURY, "Mercurio"),
    (ephem.VENUS,   "Venere"),
    (ephem.MARS,    "Marte"),
    (ephem.JUPITER, "Giove"),
    (ephem.SATURN,  "Saturno"),
]

for body_id, nome in corpi:
    pos, _ = ephem.calc_ut(jd_ritorno, body_id, 0)
    segno = segni[int(pos[0] / 30)]
    gradi = pos[0] % 30
    d, m, s, sf, si = ephem.split_deg(gradi, 0)
    print(f"{nome:10s}  {d:2d}° {m:02d}'  {segno}")

# Case Placidus per Roma
cusps, ascmc = ephem.houses(
    jd_ritorno, 41.9, 12.5, ord("P")
)
print()
print(f"Ascendente RS:  {ascmc[0]:.2f}°"
      f" ({segni[int(ascmc[0] / 30)]})")
print(f"Medio Cielo RS: {ascmc[1]:.2f}°"
      f" ({segni[int(ascmc[1] / 30)]})")
```

```
--- Rivoluzione Solare 2024 ---

Sole natale: 354.556761°
Ritorno:     14/03/2024  15:50:28 UT

Sole        24° 33'  Psc
Luna        23° 25'  Tau
Mercurio     8° 25'  Ari
Venere       3° 24'  Psc
Marte       23° 31'  Aqr
Giove       13° 48'  Tau
Saturno     11° 34'  Psc

Ascendente RS:  158.76° (Vir)
Medio Cielo RS: 64.82° (Gem)
```

---

## Ricetta 11 — Ore planetarie

**Problema**: calcolare le ore planetarie del giorno, dividendo il tempo tra alba e tramonto in 12 parti e assegnando a ciascuna un pianeta secondo la sequenza caldea.

**A cosa serve**: le ore planetarie sono uno dei sistemi più antichi di organizzazione del tempo, risalenti a Babilonia. L'idea è che ogni ora del giorno sia "governata" da un pianeta, seguendo la sequenza caldea (Saturno → Giove → Marte → Sole → Venere → Mercurio → Luna, e poi si ricomincia). La prima ora diurna è governata dal pianeta del giorno: Lunedì = Luna, Martedì = Marte, ecc. Un'ora planetaria non dura 60 minuti — dura 1/12 del periodo di luce (d'estate le ore diurne sono più lunghe). In astrologia oraria e magia cerimoniale, si sceglie l'ora del pianeta appropriato per iniziare un'attività: l'ora di Venere per questioni d'amore, l'ora di Mercurio per firmare contratti, l'ora di Giove per affari finanziari. Un'applicazione concreta: un'app che mostra in tempo reale quale ora planetaria è attiva adesso, con un countdown alla prossima.

```python
import libephemeris as ephem

lat, lon = 41.9, 12.5  # Roma
jd_noon = ephem.julday(2024, 4, 8, 12.0)

# Trova alba, tramonto e alba successiva
_, tret = ephem.rise_trans(
    jd_noon - 0.5, ephem.SUN, 1, (lon, lat, 0.0)
)
jd_alba = tret[0]
_, tret = ephem.rise_trans(
    jd_noon, ephem.SUN, 2, (lon, lat, 0.0)
)
jd_tram = tret[0]
_, tret = ephem.rise_trans(
    jd_tram, ephem.SUN, 1, (lon, lat, 0.0)
)
jd_alba_dom = tret[0]

_, _, _, o_alba = ephem.revjul(jd_alba)
_, _, _, o_tram = ephem.revjul(jd_tram)
print(f"Alba:     {int(o_alba):02d}:"
      f"{int((o_alba - int(o_alba)) * 60):02d} UT")
print(f"Tramonto: {int(o_tram):02d}:"
      f"{int((o_tram - int(o_tram)) * 60):02d} UT")

# Durata di un'ora planetaria
ora_diurna = (jd_tram - jd_alba) / 12
ora_notturna = (jd_alba_dom - jd_tram) / 12
print(f"Durata ora diurna:   {ora_diurna * 24 * 60:.1f} min")
print(f"Durata ora notturna: {ora_notturna * 24 * 60:.1f} min")
print()

# Sequenza caldea
caldea = ["Saturno", "Giove", "Marte", "Sole",
          "Venere", "Mercurio", "Luna"]
# 8 aprile 2024 = Lunedì → Luna (indice 6)
gov = 6

print("--- Ore Planetarie (Lunedì 8/4/2024, Roma) ---")
print()
print("Ore diurne:")
idx = gov
for i in range(12):
    ini = jd_alba + i * ora_diurna
    fin = jd_alba + (i + 1) * ora_diurna
    _, _, _, oi = ephem.revjul(ini)
    _, _, _, of_ = ephem.revjul(fin)
    print(f"  Ora {i+1:2d}  "
          f"{int(oi):02d}:{int((oi-int(oi))*60):02d}-"
          f"{int(of_):02d}:{int((of_-int(of_))*60):02d} UT"
          f"  {caldea[idx % 7]}")
    idx += 1

print()
print("Ore notturne:")
for i in range(12):
    ini = jd_tram + i * ora_notturna
    fin = jd_tram + (i + 1) * ora_notturna
    _, _, _, oi = ephem.revjul(ini)
    _, _, _, of_ = ephem.revjul(fin)
    print(f"  Ora {i+1:2d}  "
          f"{int(oi):02d}:{int((oi-int(oi))*60):02d}-"
          f"{int(of_):02d}:{int((of_-int(of_))*60):02d} UT"
          f"  {caldea[idx % 7]}")
    idx += 1
```

```
Alba:     04:40 UT
Tramonto: 17:43 UT
Durata ora diurna:   65.2 min
Durata ora notturna: 54.6 min

--- Ore Planetarie (Lunedì 8/4/2024, Roma) ---

Ore diurne:
  Ora  1  04:40-05:45 UT  Luna
  Ora  2  05:45-06:51 UT  Saturno
  Ora  3  06:51-07:56 UT  Giove
  Ora  4  07:56-09:01 UT  Marte
  Ora  5  09:01-10:06 UT  Sole
  Ora  6  10:06-11:12 UT  Venere
  Ora  7  11:12-12:17 UT  Mercurio
  Ora  8  12:17-13:22 UT  Luna
  Ora  9  13:22-14:27 UT  Saturno
  Ora 10  14:27-15:33 UT  Giove
  Ora 11  15:33-16:38 UT  Marte
  Ora 12  16:38-17:43 UT  Sole

Ore notturne:
  Ora  1  17:43-18:38 UT  Venere
  Ora  2  18:38-19:32 UT  Mercurio
  Ora  3  19:32-20:27 UT  Luna
  Ora  4  20:27-21:22 UT  Saturno
  Ora  5  21:22-22:16 UT  Giove
  Ora  6  22:16-23:11 UT  Marte
  Ora  7  23:11-00:05 UT  Sole
  Ora  8  00:05-01:00 UT  Venere
  Ora  9  01:00-01:55 UT  Mercurio
  Ora 10  01:55-02:49 UT  Luna
  Ora 11  02:49-03:44 UT  Saturno
  Ora 12  03:44-04:38 UT  Giove
```

> Nota come le ore diurne (65 min) sono più lunghe di quelle notturne (55 min) — ad aprile nell'emisfero nord il giorno è più lungo della notte.

---

## Ricetta 12 — Efemeride mensile stampabile

**Problema**: generare una tabella con la posizione di tutti i pianeti per ogni giorno del mese, in formato zodiacale.

**A cosa serve**: l'efemeride è lo strumento di lavoro dell'astrologo — una tabella che mostra "dove si trovano i pianeti" ogni giorno a mezzogiorno. Prima dei computer, gli astrologi consultavano voluminosi libri di efemeridi (come la celebre *Raphael's Ephemeris*) per calcolare manualmente le carte. Anche oggi, una tabella compatta è utile per avere una visione d'insieme del mese: si vede a colpo d'occhio quando un pianeta cambia segno, quando la Luna transita rapidamente attraverso lo zodiaco, quando Mercurio rallenta prima della retrogradazione. Un'applicazione concreta: un generatore di efemeridi PDF personalizzate per un anno qualsiasi.

```python
import libephemeris as ephem

segni = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

def fmt(lon):
    """Formatta longitudine in gradi + segno abbreviato."""
    s = segni[int(lon / 30)]
    g = lon % 30
    return f"{g:5.1f}{s}"

corpi = [
    (ephem.SUN,     "Sole"),
    (ephem.MOON,    "Luna"),
    (ephem.MERCURY, "Merc"),
    (ephem.VENUS,   "Vene"),
    (ephem.MARS,    "Mart"),
    (ephem.JUPITER, "Giov"),
    (ephem.SATURN,  "Satu"),
]

print("--- Efemeride Aprile 2024 (mezzogiorno UT) ---")
print()
header = "Data     "
for _, nome in corpi:
    header += f"  {nome:>8s}"
print(header)
print("-" * len(header))

for giorno in range(1, 31):
    jd = ephem.julday(2024, 4, giorno, 12.0)
    riga = f"{giorno:2d}/04/2024"
    for body_id, _ in corpi:
        pos, _ = ephem.calc_ut(jd, body_id, 0)
        riga += f"  {fmt(pos[0]):>8s}"
    print(riga)
```

```
--- Efemeride Aprile 2024 (mezzogiorno UT) ---

Data           Sole      Luna      Merc      Vene      Mart      Giov      Satu
-------------------------------------------------------------------------------
 1/04/2024   12.2Ari    4.4Cap   27.2Ari   25.5Psc    7.4Psc   17.5Tau   13.7Psc
 2/04/2024   13.2Ari   17.8Cap   27.2Ari   26.7Psc    8.2Psc   17.7Tau   13.8Psc
 3/04/2024   14.2Ari    1.7Aqr   27.1Ari   27.9Psc    9.0Psc   17.9Tau   13.9Psc
 4/04/2024   15.2Ari   15.9Aqr   26.8Ari   29.2Psc    9.7Psc   18.1Tau   14.0Psc
 5/04/2024   16.2Ari    0.5Psc   26.5Ari    0.4Ari   10.5Psc   18.3Tau   14.1Psc
 6/04/2024   17.2Ari   15.3Psc   26.1Ari    1.6Ari   11.3Psc   18.5Tau   14.2Psc
 7/04/2024   18.2Ari    0.4Ari   25.6Ari    2.9Ari   12.1Psc   18.8Tau   14.3Psc
 8/04/2024   19.1Ari   15.4Ari   25.0Ari    4.1Ari   12.8Psc   19.0Tau   14.4Psc
 9/04/2024   20.1Ari    0.4Tau   24.3Ari    5.4Ari   13.6Psc   19.2Tau   14.5Psc
10/04/2024   21.1Ari   15.1Tau   23.6Ari    6.6Ari   14.4Psc   19.4Tau   14.6Psc
11/04/2024   22.1Ari   29.4Tau   22.9Ari    7.8Ari   15.2Psc   19.7Tau   14.7Psc
12/04/2024   23.1Ari   13.4Gem   22.1Ari    9.1Ari   16.0Psc   19.9Tau   14.9Psc
13/04/2024   24.0Ari   26.8Gem   21.4Ari   10.3Ari   16.7Psc   20.1Tau   15.0Psc
14/04/2024   25.0Ari    9.9Cnc   20.6Ari   11.5Ari   17.5Psc   20.3Tau   15.1Psc
15/04/2024   26.0Ari   22.6Cnc   19.9Ari   12.8Ari   18.3Psc   20.6Tau   15.2Psc
16/04/2024   27.0Ari    4.9Leo   19.2Ari   14.0Ari   19.1Psc   20.8Tau   15.3Psc
17/04/2024   28.0Ari   17.0Leo   18.6Ari   15.2Ari   19.8Psc   21.0Tau   15.4Psc
18/04/2024   28.9Ari   28.9Leo   18.0Ari   16.5Ari   20.6Psc   21.2Tau   15.5Psc
19/04/2024   29.9Ari   10.8Vir   17.5Ari   17.7Ari   21.4Psc   21.5Tau   15.6Psc
20/04/2024    0.9Tau   22.5Vir   17.0Ari   18.9Ari   22.2Psc   21.7Tau   15.7Psc
21/04/2024    1.9Tau    4.4Lib   16.7Ari   20.2Ari   22.9Psc   21.9Tau   15.8Psc
22/04/2024    2.8Tau   16.3Lib   16.4Ari   21.4Ari   23.7Psc   22.2Tau   15.9Psc
23/04/2024    3.8Tau   28.3Lib   16.2Ari   22.6Ari   24.5Psc   22.4Tau   16.0Psc
24/04/2024    4.8Tau   10.5Sco   16.0Ari   23.9Ari   25.3Psc   22.6Tau   16.0Psc
25/04/2024    5.8Tau   22.9Sco   16.0Ari   25.1Ari   26.0Psc   22.8Tau   16.1Psc
26/04/2024    6.7Tau    5.5Sgr   16.0Ari   26.3Ari   26.8Psc   23.1Tau   16.2Psc
27/04/2024    7.7Tau   18.3Sgr   16.1Ari   27.6Ari   27.6Psc   23.3Tau   16.3Psc
28/04/2024    8.7Tau    1.3Cap   16.3Ari   28.8Ari   28.3Psc   23.5Tau   16.4Psc
29/04/2024    9.7Tau   14.6Cap   16.6Ari    0.0Tau   29.1Psc   23.8Tau   16.5Psc
30/04/2024   10.6Tau   28.1Cap   17.0Ari    1.3Tau   29.9Psc   24.0Tau   16.6Psc
```

> Dalla tabella si legge chiaramente: Mercurio retrocede da 27.2° Ari a 16.0° Ari (stazione diretta il 25/04), poi ricomincia ad avanzare. La Luna attraversa tutto lo zodiaco in un mese. Il Sole avanza di circa 1° al giorno attraverso Ariete e Toro.

---

## Riepilogo

Questo capitolo ha presentato 12 ricette complete per i calcoli astronomici e astrologici più comuni:

- **Ricetta 1** — Tema natale completo con pianeti, case e angoli
- **Ricetta 2** — Transiti del giorno su una carta natale
- **Ricetta 3** — Ingressi del Sole nei 12 segni (equinozi e solstizi)
- **Ricetta 4** — Retrogradazioni di Mercurio (stazioni e verifica)
- **Ricetta 5** — Lune Nuove e Piene dell'anno
- **Ricetta 6** — Prossime eclissi solari visibili da una città
- **Ricetta 7** — Carta siderale vedica con Nakshatra e Pada
- **Ricetta 8** — Visibilità dei pianeti stasera (altezza, magnitudine, direzione)
- **Ricetta 9** — Aspetti tra pianeti con orbe
- **Ricetta 10** — Rivoluzione solare (Solar Return)
- **Ricetta 11** — Ore planetarie (sequenza caldea)
- **Ricetta 12** — Efemeride mensile stampabile

Ogni ricetta è autocontenuta e pronta all'uso: basta copiare il codice, modificare data e coordinate, ed eseguirlo.

## Funzioni utilizzate nelle ricette

- `julday(y, m, d, h)` / `revjul(jd)` — conversione data ↔ Julian Day
- `calc_ut(jd, body, flag)` — posizione di un corpo celeste
- `houses(jd, lat, lon, hsys)` — cuspidi delle case e angoli
- `split_deg(deg, flag)` — scompone gradi decimali in °, ', "
- `difdeg2n(p1, p2)` — differenza angolare normalizzata [-180, 180]
- `solcross_ut(x, jd)` — attraversamento del Sole a una longitudine
- `mooncross_ut(x, jd)` — attraversamento della Luna a una longitudine
- `find_station_ut(body, jd, tipo)` — prossima stazione retrograda/diretta (da `libephemeris.crossing`)
- `is_retrograde(body, jd)` — verifica retrogradazione (da `libephemeris.crossing`)
- `sol_eclipse_when_loc(jd, geopos)` — eclissi solare locale
- `set_sid_mode(mode)` / `get_ayanamsa_ut(jd)` — zodiaco siderale
- `rise_trans(jd, body, rsmi, geopos)` — alba e tramonto
- `azalt(jd, flag, geopos, press, temp, xin)` — coordinate orizzontali
- `pheno_ut(jd, body, flag)` — fenomeni (magnitudine, fase, elongazione)
- `set_topo(lon, lat, alt)` — posizione dell'osservatore
