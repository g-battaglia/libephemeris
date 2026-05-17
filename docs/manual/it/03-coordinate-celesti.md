# Capitolo 3 — Dove si trova un pianeta: coordinate celesti

## Cosa imparerai

In questo capitolo imparerai i tre sistemi di coordinate usati in astronomia — eclittiche, equatoriali e orizzontali — e come convertire tra di essi. Capirai anche cosa sono i sistemi di riferimento (J2000, ICRS), e gli effetti di precessione, nutazione e aberrazione.

---

## 3.1 Coordinate eclittiche: longitudine e latitudine

Il sistema di coordinate più usato in astrologia. Il piano di riferimento è l'eclittica (il percorso del Sole), il punto zero è l'equinozio di primavera.

- **Longitudine eclittica**: 0°–360° lungo l'eclittica. È il "segno zodiacale".
- **Latitudine eclittica**: distanza dall'eclittica, da -90° a +90°.

Il Sole ha sempre latitudine ~0° (per definizione, è lui che traccia l'eclittica). La Luna può arrivare a ±5.1° di latitudine. I pianeti restano entro pochi gradi dall'eclittica, tranne Plutone che può raggiungere ~17°.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)
pos, flag = ephem.calc_ut(jd, ephem.MOON, 0)

lon = pos[0]   # longitudine eclittica (gradi)
lat = pos[1]   # latitudine eclittica (gradi)
dist = pos[2]  # distanza (UA)

segno_num = int(lon / 30)
segni = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]
print(f"Luna: {lon % 30:.2f}° {segni[segno_num]}, lat {lat:.2f}°")
# "Luna a 15° Cancro" = longitudine eclittica 105.x°
```

```
Luna: 15.43° Ari, lat -0.02°
```

---

## 3.2 Coordinate equatoriali: ascensione retta e declinazione

Il sistema usato dai telescopi. Il piano di riferimento è l'equatore celeste.

- **Ascensione Retta** (RA): 0°–360° lungo l'equatore celeste (a volte espressa in ore: 0h–24h).
- **Declinazione** (Dec): distanza dall'equatore celeste, da -90° a +90°.

I telescopi preferiscono questo sistema perché l'equatore celeste è allineato con la rotazione terrestre: basta ruotare attorno a un solo asse per inseguire un oggetto.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Flag FLG_EQUATORIAL: restituisce RA e Dec invece di lon e lat
pos, flag = ephem.calc_ut(jd, ephem.MOON, ephem.FLG_EQUATORIAL)

ra = pos[0]    # Ascensione Retta in gradi (0–360)
dec = pos[1]   # Declinazione in gradi (-90 a +90)

# Conversione RA da gradi a ore:minuti:secondi
ra_ore = ra / 15.0
ra_h = int(ra_ore)
ra_m = int((ra_ore - ra_h) * 60)
ra_s = ((ra_ore - ra_h) * 60 - ra_m) * 60

print(f"Luna: RA = {ra_h}h {ra_m}m {ra_s:.1f}s, Dec = {dec:+.2f}°")
```

```
Luna: RA = 0h 56m 52.2s, Dec = +6.06°
```

---

## 3.3 Coordinate orizzontali: altezza e azimut

Il sistema "pratico" — dove devo guardare per vedere un oggetto nel cielo? Questo sistema dipende dalla posizione dell'osservatore e dal momento.

- **Altezza**: angolo sopra l'orizzonte. 0° = orizzonte, 90° = zenit. Valori negativi = sotto l'orizzonte.
- **Azimut**: direzione sull'orizzonte. In LibEphemeris: 0° = Sud, 90° = Ovest, 180° = Nord, 270° = Est.

Servono la posizione dell'osservatore e, opzionalmente, pressione e temperatura per la rifrazione atmosferica.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 9, 15, 21.0)

# Giove in coordinate equatoriali
pos, _ = ephem.calc_ut(jd, ephem.JUPITER, ephem.FLG_EQUATORIAL)

# Roma (lon Est, lat Nord, altitudine in metri)
geopos = (12.4964, 41.9028, 50.0)

# Da equatoriali a orizzontali (SE_EQU2HOR = 1)
hor = ephem.azalt(jd, ephem.EQU2HOR, geopos, 1013.25, 15.0,
                  (pos[0], pos[1], pos[2]))

print(f"Azimut: {hor[0]:.1f}° (da Sud)")
print(f"Altezza vera: {hor[1]:.1f}°")
print(f"Altezza apparente: {hor[2]:.1f}° (con rifrazione)")

# Operazione inversa: da orizzontali a equatoriali
ra_dec = ephem.azalt_rev(jd, ephem.HOR2EQU, geopos, hor[0], hor[1])
print(f"RA recuperata: {ra_dec[0]:.4f}°, Dec: {ra_dec[1]:.4f}°")
```

```
Azimut: 235.7° (da Sud)
Altezza vera: -3.2°
Altezza apparente: -3.2° (con rifrazione)
RA recuperata: 79.6564°, Dec: 22.3890°
```

### 🌍 Vita reale

"Giove è visibile stasera?" — calcola l'altezza: se è > 0°, è sopra l'orizzonte. Se è > 10°, è abbastanza alto da non essere disturbato dall'atmosfera vicino all'orizzonte. L'azimut ti dice in che direzione guardare.

---

## 3.4 Sistemi di riferimento: ICRS, J2000, equinozio del giorno

Le coordinate di un oggetto celeste dipendono dal **sistema di riferimento** scelto. Il punto zero e gli assi possono essere definiti in modi diversi:

- **Equinozio del giorno** ("of date"): il default in LibEphemeris. Le coordinate tengono conto della precessione e della nutazione attuali. Il punto zero è l'equinozio di primavera *di quel momento*. È il sistema usato in astrologia.

- **J2000.0**: le coordinate sono fissate al 1° gennaio 2000, ore 12:00 TT. Il punto zero è l'equinozio di primavera come era nel 2000. È lo standard dei cataloghi astronomici.

- **ICRS**: il sistema moderno, basato sulle posizioni di quasar lontanissimi (praticamente fissi). È quasi identico a J2000 — la differenza ("frame bias") è di ~0.02".

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Coordinate eclittiche "of date" (default)
pos_date, _ = ephem.calc_ut(jd, ephem.MARS, 0)

# Coordinate eclittiche J2000
pos_j2000, _ = ephem.calc_ut(jd, ephem.MARS, ephem.FLG_J2000)

# Coordinate equatoriali ICRS
pos_icrs, _ = ephem.calc_ut(jd, ephem.MARS,
                            ephem.FLG_EQUATORIAL | ephem.FLG_ICRS)

print(f"Marte of-date: {pos_date[0]:.4f}°")
print(f"Marte J2000:   {pos_j2000[0]:.4f}°")
print(f"Differenza:    {pos_date[0] - pos_j2000[0]:.4f}° (≈ precessione)")
```

```
Marte of-date: 342.8458°
Marte J2000:   342.5082°
Differenza:    0.3376° (≈ precessione)
```

La differenza tra "of date" e J2000 è principalmente la precessione accumulata dal 2000 ad oggi — circa 0.34° nel 2024.

---

## 3.5 Precessione, nutazione e aberrazione

Tre effetti fisici influenzano le coordinate celesti. È importante sapere che esistono, anche se nella pratica la libreria li gestisce automaticamente.

### Precessione

L'asse di rotazione terrestre non punta sempre nella stessa direzione: disegna un cerchio nel cielo in circa 26000 anni, come una trottola che oscilla. L'effetto è che il punto vernale (0° Ariete) si sposta di circa 50" d'arco all'anno rispetto alle stelle.

In 2000 anni, la precessione ha spostato il punto vernale di circa 28° — quasi un intero segno zodiacale. È la ragione per cui zodiaco tropicale e costellazioni non coincidono più.

### Nutazione

Sovrapposta alla precessione, l'asse terrestre compie piccole oscillazioni con periodo principale di ~18.6 anni (legate al ciclo dei nodi lunari). L'ampiezza è di ~9" in longitudine e ~17" in obliquità.

Il flag `FLG_NONUT` calcola le coordinate "medie" (senza nutazione) invece di quelle "vere":

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Con nutazione (default = coordinate "vere")
pos_vera, _ = ephem.calc_ut(jd, ephem.SUN, 0)

# Senza nutazione (coordinate "medie")
pos_media, _ = ephem.calc_ut(jd, ephem.SUN, ephem.FLG_NONUT)

diff_arcsec = (pos_vera[0] - pos_media[0]) * 3600
print(f"Effetto nutazione: {diff_arcsec:.2f}\" d'arco")
```

```
Effetto nutazione: -5.30" d'arco
```

### Aberrazione

La luce impiega tempo ad arrivare, e nel frattempo la Terra si muove. L'effetto è uno spostamento apparente fino a ~20" nella direzione del moto terrestre. Il flag `FLG_NOABERR` lo disabilita; `FLG_TRUEPOS` restituisce la posizione geometrica pura (senza correzione per il tempo-luce né per l'aberrazione).

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Posizione apparente (default: con aberrazione)
pos_app, _ = ephem.calc_ut(jd, ephem.MARS, 0)

# Posizione geometrica vera (senza aberrazione né tempo-luce)
pos_true, _ = ephem.calc_ut(jd, ephem.MARS, ephem.FLG_TRUEPOS)

diff = (pos_app[0] - pos_true[0]) * 3600
print(f"Effetto aberrazione + tempo-luce: {diff:.1f}\" d'arco")
```

```
Effetto aberrazione + tempo-luce: -33.3" d'arco
```

---

## 3.6 Coordinate cartesiane (XYZ)

A volte servono le posizioni in coordinate cartesiane (X, Y, Z) invece che angolari. Utile per calcoli di distanza 3D, angoli di fase, o problemi geometrici.

Il flag `FLG_XYZ` cambia il significato della tupla restituita: invece di (lon, lat, dist, vel_lon, vel_lat, vel_dist) si ottiene (x, y, z, vx, vy, vz) in unità astronomiche e UA/giorno.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Posizione cartesiana equatoriale
pos, _ = ephem.calc_ut(jd, ephem.MARS,
                       ephem.FLG_XYZ | ephem.FLG_EQUATORIAL)

x, y, z = pos[0], pos[1], pos[2]
import math
dist = math.sqrt(x**2 + y**2 + z**2)
print(f"Marte: X={x:.4f}, Y={y:.4f}, Z={z:.4f} UA")
print(f"Distanza: {dist:.4f} UA = {dist * 149597870.7:.0f} km")
```

```
Marte: X=1.9696, Y=-0.5400, Z=-0.2829 UA
Distanza: 2.0618 UA = 308436754 km
```

---

## 3.7 Convertire tra sistemi

LibEphemeris offre funzioni dedicate per le conversioni:

| Funzione | Conversione |
|----------|------------|
| `cotrans(coord, -obliquità)` | Eclittica → Equatoriale |
| `cotrans(coord, +obliquità)` | Equatoriale → Eclittica |
| `cotrans_sp(coord, speed, obl)` | Come sopra, con velocità |
| `azalt(jd, SE_ECL2HOR, ...)` | Eclittica → Orizzontale |
| `azalt(jd, SE_EQU2HOR, ...)` | Equatoriale → Orizzontale |
| `azalt_rev(jd, SE_HOR2ECL, ...)` | Orizzontale → Eclittica |
| `azalt_rev(jd, SE_HOR2EQU, ...)` | Orizzontale → Equatoriale |

Il segno dell'obliquità in `cotrans` controlla la direzione: **negativo** = eclittica→equatoriale, **positivo** = equatoriale→eclittica.

### Esempio completo: dalla posizione eclittica all'altezza nel cielo

```python
import libephemeris as ephem

jd = ephem.julday(2024, 9, 15, 21.0)

# 1. Posizione eclittica di Saturno
pos, _ = ephem.calc_ut(jd, ephem.SATURN, 0)
lon, lat, dist = pos[0], pos[1], pos[2]

# 2. Otteniamo l'obliquità
nut, _ = ephem.calc_ut(jd, ephem.ECL_NUT, 0)
obliquita = nut[0]

# 3. Eclittica -> Equatoriale (obliquità negativa)
ra, dec, d = ephem.cotrans((lon, lat, dist), -obliquita)

# 4. Equatoriale -> Orizzontale (Roma)
geopos = (12.4964, 41.9028, 50.0)
hor = ephem.azalt(jd, ephem.EQU2HOR, geopos, 1013.25, 15.0,
                  (ra, dec, d))

print(f"Saturno: lon={lon:.2f}°, lat={lat:.2f}°")
print(f"         RA={ra:.2f}°, Dec={dec:+.2f}°")
print(f"         Altezza={hor[2]:.1f}°, Azimut={hor[0]:.1f}° (da S)")
```

```
Saturno: lon=345.44°, lat=-2.20°
         RA=347.46°, Dec=-7.76°
         Altezza=35.5°, Azimut=329.5° (da S)
```

> **Scorciatoia**: puoi anche passare direttamente le coordinate eclittiche ad `azalt` con `SE_ECL2HOR`, evitando la conversione manuale. La conversione esplicita via `cotrans` è utile quando ti servono i valori intermedi.

---

## Riepilogo

- **Eclittiche** (lon, lat): sistema dell'astrologia, basato sull'eclittica. Default di `calc_ut`.
- **Equatoriali** (RA, Dec): sistema dei telescopi, basato sull'equatore celeste. Flag `FLG_EQUATORIAL`.
- **Orizzontali** (azimut, altezza): sistema pratico, dipende dal luogo. Funzione `azalt`.
- **Of date** (default), **J2000** (`FLG_J2000`), **ICRS** (`FLG_ICRS`): tre sistemi di riferimento per lo zero delle coordinate.
- **Precessione**: spostamento lento del punto zero (~50"/anno). **Nutazione**: oscillazione dell'asse (~9"). **Aberrazione**: spostamento dovuto al moto terrestre (~20").
- `cotrans` converte tra eclittica ed equatoriale. `azalt` e `azalt_rev` convertono da/verso le coordinate orizzontali.

### Funzioni e costanti introdotte

| Funzione / Costante | Uso |
|---------------------|-----|
| `cotrans(coord, obliquità)` | Eclittica ↔ Equatoriale |
| `cotrans_sp(coord, speed, obliquità)` | Come sopra, con velocità |
| `azalt_rev(jd, flag, geopos, az, alt)` | Orizzontale → Ecl/Eq |
| `FLG_EQUATORIAL` | Coordinate equatoriali |
| `FLG_XYZ` | Coordinate cartesiane |
| `FLG_J2000`, `FLG_ICRS` | Sistemi di riferimento |
| `FLG_NONUT`, `FLG_NOABERR`, `FLG_TRUEPOS` | Disabilita correzioni |
| `SE_ECL2HOR`, `SE_EQU2HOR` | Direzione per `azalt` |
| `SE_HOR2ECL`, `SE_HOR2EQU` | Direzione per `azalt_rev` |
