# Capitolo 13 — Pianeti ipotetici e parti arabe

## Cosa imparerai

In questo capitolo scoprirai cosa sono i pianeti uraniani (punti matematici usati nella Scuola di Amburgo), altri corpi ipotetici come Transpluto e Vulcano, come definire orbite fittizie personalizzate, e come calcolare le parti arabe — formule antichissime che combinano posizioni planetarie e angoli.

---

## 13.1 I pianeti uraniani (Scuola di Amburgo)

Negli anni '20, l'astrologo tedesco Alfred Witte fondò la **Scuola di Amburgo** e postulò l'esistenza di otto "pianeti" trans-nettuniani. Questi corpi non sono mai stati osservati — non esistono fisicamente. Sono **punti matematici** con orbite definite da elementi kepleriani, usati esclusivamente nell'astrologia uraniana e nella tecnica dei "punti medi".

Gli otto pianeti uraniani sono:

- **Cupido** (`SE_CUPIDO`, ID 40) — associato alla famiglia, all'arte e ai gruppi
- **Hades** (`SE_HADES`, ID 41) — associato al passato, alla povertà, alla malattia
- **Zeus** (`SE_ZEUS`, ID 42) — associato al fuoco, alle macchine, alla forza dirigente
- **Kronos** (`SE_KRONOS`, ID 43) — associato all'autorità, al governo, all'eccellenza
- **Apollon** (`SE_APOLLON`, ID 44) — associato all'espansione, alla scienza, al commercio
- **Admetos** (`SE_ADMETOS`, ID 45) — associato alla profondità, alla concentrazione, ai blocchi
- **Vulkanus** (`SE_VULKANUS`, ID 46) — associato alla potenza, all'intensità, alla forza
- **Poseidon** (`SE_POSEIDON`, ID 47) — associato alla mente, all'illuminazione, alla verità

### Calcolare la posizione

Ogni pianeta uraniano ha la sua funzione dedicata:

```python
import libephemeris as ephem

# Convertiamo la data in TT (questi funzionano in TT, non UT)
jd_ut = ephem.julday(2024, 4, 8, 12.0)
delta_t = ephem.deltat(jd_ut)  # in giorni
jd_tt = jd_ut + delta_t

# Posizione di Kronos
pos = ephem.calc_kronos(jd_tt)
lon, lat, dist = pos[0], pos[1], pos[2]
vel = pos[3]  # velocità in gradi/giorno

segni = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]
segno = segni[int(lon / 30)]
gradi = lon % 30

print(f"Kronos: {gradi:.2f}° {segno}")
print(f"Velocità: {vel:.4f}°/giorno")
```

```
Kronos: 14.95° Cnc
Velocità: 0.0019°/giorno
```

Puoi anche usare la funzione generica `calc_uranian_planet` con l'ID del corpo:

```python
import libephemeris as ephem

jd_ut = ephem.julday(2024, 4, 8, 12.0)
jd_tt = jd_ut + ephem.deltat(jd_ut)

segni = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

uraniani = [
    (ephem.CUPIDO,   "Cupido"),
    (ephem.HADES,    "Hades"),
    (ephem.ZEUS,     "Zeus"),
    (ephem.KRONOS,   "Kronos"),
    (ephem.APOLLON,  "Apollon"),
    (ephem.ADMETOS,  "Admetos"),
    (ephem.VULKANUS, "Vulkanus"),
    (ephem.POSEIDON, "Poseidon"),
]

for body_id, nome in uraniani:
    pos = ephem.calc_uranian_planet(body_id, jd_tt)
    segno = segni[int(pos[0] / 30)]
    gradi = pos[0] % 30
    print(f"{nome:10s}  {gradi:5.1f}° {segno}")
```

```
Cupido        7.5° Cap
Hades        12.9° Cnc
Zeus         25.1° Lib
Kronos       14.9° Cnc
Apollon       6.5° Sco
Admetos       3.3° Gem
Vulkanus      3.8° Leo
Poseidon     16.3° Sco
```

I pianeti uraniani possono anche essere calcolati con `calc_ut` usando i loro ID, ma le funzioni dedicate accettano direttamente JD in TT.

---

## 13.2 Altri corpi ipotetici

Oltre agli uraniani, la libreria include altri corpi ipotetici che sono stati proposti nel corso della storia ma mai confermati osservativamente:

**Transpluto / Isis** (`SE_ISIS`, ID 48) — Un ipotetico pianeta oltre Plutone, postulato prima della scoperta di Eris. L'orbita usata è basata sulla proposta di Theodore Landscheidt.

```python
import libephemeris as ephem

jd_tt = ephem.julday(2024, 4, 8, 12.0) + ephem.deltat(ephem.julday(2024, 4, 8, 12.0))
pos = ephem.calc_transpluto(jd_tt)
print(f"Transpluto: {pos[0]:.2f}°")
```

```
Transpluto: 153.47°
```

**Vulcano** (`SE_VULCAN`, ID 55) — Un ipotetico pianeta tra Mercurio e il Sole, cercato per tutto il XIX secolo per spiegare le anomalie nell'orbita di Mercurio. La relatività generale di Einstein ha poi spiegato quelle anomalie senza bisogno di un pianeta aggiuntivo — ma il concetto resta nell'astrologia esoterica.

```python
pos = ephem.calc_vulcan(jd_tt)
print(f"Vulcano: {pos[0]:.2f}°")
```

```
Vulcano: 45.53°
```

**Luna Bianca / Selena** (`SE_WHITE_MOON`, ID 56) — Il punto diametralmente opposto alla Lilith Nera (l'apogeo lunare medio). Non è un corpo fisico, ma un punto simbolico usato in alcune scuole di astrologia come "complemento luminoso" di Lilith.

```python
pos = ephem.calc_white_moon_position(jd_tt)
print(f"Luna Bianca: {pos[0]:.2f}°")
```

```
Luna Bianca: 350.92°
```

**Luna di Waldemath** (ID 58) — Un ipotetico secondo satellite della Terra, "osservato" nel 1898 da Georg Waldemath. Mai confermato.

**Proserpina** (ID 57) — Un altro ipotetico trans-plutoniano.

**Pianeta X di Pickering** — La previsione del 1919 di William Pickering per un pianeta oltre Nettuno.

Tutti questi corpi si calcolano con la funzione generica `calc_hypothetical_position`:

```python
import libephemeris as ephem

jd_tt = ephem.julday(2024, 4, 8, 12.0) + ephem.deltat(ephem.julday(2024, 4, 8, 12.0))

# Qualsiasi corpo ipotetico dato il suo ID
pos = ephem.calc_hypothetical_position(ephem.WALDEMATH, jd_tt)
print(f"Luna di Waldemath: {pos[0]:.2f}°")
```

```
Luna di Waldemath: 61.36°
```

---

## 13.3 Orbite fittizie personalizzate

Se hai bisogno di un corpo ipotetico non incluso nella libreria, puoi definire la tua orbita. La libreria usa un formato di file con elementi orbitali (compatibile con il formato `seorbel.txt` della Swiss Ephemeris).

### Caricare orbite predefinite

La libreria include un file di orbite fittizie predefinite:

```python
import libephemeris as ephem

# Carica le orbite predefinite
orbite = ephem.load_bundled_fictitious_orbits()
print(f"Caricate {len(orbite)} orbite fittizie")

# Cerca un corpo per nome
corpo = ephem.get_orbital_body_by_name(orbite, "Cupido")
if corpo:
    print(f"Trovato: {corpo.name}")
    print(f"Semi-asse maggiore: {corpo.semi_axis:.2f} UA")
```

```
Caricate 24 orbite fittizie
Trovato: Cupido
Semi-asse maggiore: 41.00 UA
```

### Caricare un file personalizzato

```python
import libephemeris as ephem

# Carica orbite da un file custom
orbite = ephem.parse_orbital_elements("/percorso/al/mio/file.csv")

# Calcola la posizione di un corpo
jd_tt = ephem.julday(2024, 4, 8, 12.0) + ephem.deltat(ephem.julday(2024, 4, 8, 12.0))

corpo = ephem.get_orbital_body_by_name(orbite, "MioCorpo")
if corpo:
    pos = ephem.calc_orbital_position(corpo, jd_tt)
    print(f"MioCorpo: lon={pos[0]:.2f}°, lat={pos[1]:.2f}°")
```

---

## 13.4 Le parti arabe (lotti)

Le **parti arabe** (in greco *kleros*, "lotti") sono tra le tecniche più antiche dell'astrologia. Risalgono all'astrologia ellenistica (II–III secolo d.C.) e sono state ampiamente sviluppate dagli astrologi persiani e arabi nel Medioevo.

### Come funzionano

Ogni parte araba è un punto calcolato combinando tre posizioni: in genere l'Ascendente e due pianeti. La formula base è:

**Parte = Ascendente + Pianeta A − Pianeta B**

Il risultato (normalizzato a 0°–360°) è un punto sull'eclittica con un significato specifico.

### La Parte di Fortuna

La più famosa è la **Parte di Fortuna** (*Pars Fortunae*), associata alla prosperità, alla fortuna materiale e al benessere fisico:

- **Di giorno** (Sole sopra l'orizzonte): Fortuna = ASC + Luna − Sole
- **Di notte** (Sole sotto l'orizzonte): Fortuna = ASC + Sole − Luna

La formula si inverte tra giorno e notte perché la "luminare della setta" (il Sole di giorno, la Luna di notte) funge da punto di partenza.

### Calcolare tutte le parti

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 14.5)
lat, lon = 41.9028, 12.4964  # Roma

# Calcola le posizioni necessarie
cusps, ascmc = ephem.houses(jd, lat, lon, ord('P'))
sole, _ = ephem.calc_ut(jd, ephem.SUN, ephem.FLG_SPEED)
luna, _ = ephem.calc_ut(jd, ephem.MOON, ephem.FLG_SPEED)
mercurio, _ = ephem.calc_ut(jd, ephem.MERCURY, ephem.FLG_SPEED)
venere, _ = ephem.calc_ut(jd, ephem.VENUS, ephem.FLG_SPEED)

# Prepara le posizioni
positions = {
    "Asc": ascmc[0],
    "Sun": sole[0],
    "Moon": luna[0],
    "Mercury": mercurio[0],
    "Venus": venere[0],
}

# Calcola tutte le parti arabe
parti = ephem.calc_all_arabic_parts(
    positions,
    jd=jd,
    geo_lat=lat,
    geo_lon=lon,
)

segni = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

nomi_italiani = {
    "Pars_Fortunae": "Parte di Fortuna",
    "Pars_Spiritus": "Parte di Spirito",
    "Pars_Amoris":   "Parte dell'Amore",
    "Pars_Fidei":    "Parte della Fede",
}

for chiave, lon_parte in parti.items():
    nome = nomi_italiani.get(chiave, chiave)
    segno = segni[int(lon_parte / 30)]
    gradi = lon_parte % 30
    print(f"{nome:22s}  {gradi:5.1f}° {segno}")
```

```
Parte di Fortuna         14.5° Vir
Parte di Spirito         10.0° Vir
Parte dell'Amore         27.3° Leo
Parte della Fede         20.2° Vir
```

Le quattro parti calcolate sono:

- **Parte di Fortuna** (*Pars Fortunae*) — ASC + Luna − Sole (giorno) o ASC + Sole − Luna (notte). Prosperità e benessere materiale.

- **Parte di Spirito** (*Pars Spiritus*) — l'inverso della Parte di Fortuna. Rappresenta la volontà, lo spirito e la vocazione interiore.

- **Parte dell'Amore** (*Pars Amoris*) — ASC + Venere − Sole. Le relazioni affettive e l'attrazione.

- **Parte della Fede** (*Pars Fidei*) — ASC + Mercurio − Luna. La fede, la fiducia e le convinzioni.

---

## Riepilogo

In questo capitolo abbiamo esplorato i corpi celesti non fisici usati in diverse tradizioni astrologiche.

**Concetti chiave:**

- I **pianeti uraniani** sono otto punti matematici (Cupido, Hades, Zeus, Kronos, Apollon, Admetos, Vulkanus, Poseidon) con orbite ipotetiche, usati nella Scuola di Amburgo e nella tecnica dei punti medi
- Altri **corpi ipotetici** (Transpluto, Vulcano, Luna di Waldemath) sono stati proposti storicamente ma mai confermati; la libreria li calcola usando elementi orbitali definiti
- Le **orbite fittizie personalizzate** permettono di definire qualsiasi corpo ipotetico con i propri elementi orbitali
- Le **parti arabe** sono punti calcolati dalla formula ASC + Pianeta A − Pianeta B, con la Parte di Fortuna come la più importante

**Funzioni introdotte:**

- `calc_cupido(jd_tt)`, `calc_hades(jd_tt)`, ... `calc_poseidon(jd_tt)` — calcolano ciascun pianeta uraniano, restituendo una tupla di 6 valori (lon, lat, dist, vel_lon, vel_lat, vel_dist)
- `calc_uranian_planet(body_id, jd_tt)` — versione generica per qualsiasi pianeta uraniano
- `calc_transpluto(jd_tt)`, `calc_vulcan(jd_tt)`, `calc_waldemath(jd_tt)` — altri corpi ipotetici
- `calc_white_moon_position(jd_tt)` — la Luna Bianca (opposto della Lilith Nera media)
- `calc_hypothetical_position(body_id, jd_tt)` — funzione generica per qualsiasi corpo ipotetico
- `load_bundled_fictitious_orbits()` — carica le orbite fittizie predefinite
- `parse_orbital_elements(filepath)` — carica orbite fittizie da un file personalizzato
- `get_orbital_body_by_name(elements, nome)` — cerca un corpo per nome nella lista di orbite
- `calc_orbital_position(elem, jd_tt)` — calcola la posizione di un corpo dato i suoi elementi orbitali
- `calc_all_arabic_parts(positions, jd=..., geo_lat=..., geo_lon=...)` — calcola le quattro parti arabe principali
