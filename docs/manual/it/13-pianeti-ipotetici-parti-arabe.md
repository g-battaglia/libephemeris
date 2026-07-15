# Capitolo 13 — Pianeti ipotetici e parti arabe

## Cosa imparerai

In questo capitolo scoprirai cosa sono i pianeti uraniani (punti matematici
usati nella Scuola di Amburgo), come viene classificata la provenienza dei
modelli integrati, come definire orbite personalizzate e calcolare le
parti arabe — formule antichissime che combinano posizioni planetarie e angoli.

---

## 13.1 I pianeti uraniani (Scuola di Amburgo)

Negli anni '20, l'astrologo tedesco Alfred Witte fondò la **Scuola di Amburgo** e postulò l'esistenza di otto "pianeti" trans-nettuniani. Questi corpi non sono mai stati osservati — non esistono fisicamente. Sono **punti matematici** con orbite definite da elementi kepleriani, usati esclusivamente nell'astrologia uraniana e nella tecnica dei "punti medi".

Gli otto pianeti uraniani sono:

- **Cupido** (`CUPIDO`, ID 40) — associato alla famiglia, all'arte e ai gruppi
- **Hades** (`HADES`, ID 41) — associato al passato, alla povertà, alla malattia
- **Zeus** (`ZEUS`, ID 42) — associato al fuoco, alle macchine, alla forza dirigente
- **Kronos** (`KRONOS`, ID 43) — associato all'autorità, al governo, all'eccellenza
- **Apollon** (`APOLLON`, ID 44) — associato all'espansione, alla scienza, al commercio
- **Admetos** (`ADMETOS`, ID 45) — associato alla profondità, alla concentrazione, ai blocchi
- **Vulkanus** (`VULKANUS`, ID 46) — associato alla potenza, all'intensità, alla forza
- **Poseidon** (`POSEIDON`, ID 47) — associato alla mente, all'illuminazione, alla verità

### Disponibilità

Gli ID 40–47 vengono calcolati da una nuova trascrizione di “Orbital Elements
for the Transneptunian Planets” di James Neely, *Matrix Magazine* VII (1980),
Tabella I, p. 10. La fonte definisce il tempo in secoli giuliani da
JD 2415020.0; elementi letterali, scelta del frame e propagazione a due corpi
sono descritti nel
[registro di provenienza](../../methodology/hypothetical-bodies.md).

---

## 13.2 Corpi ipotetici integrati

Tredici ID vengono calcolati da ricostruzioni documentali verificate pagina per
pagina: **Cupido–Poseidon** (40–47), **Harrington** (50), **Le Verrier** (51),
**Adams** (52), **Lowell** (53) e **Luna Bianca / Selena** (56). Gli ID 48,
49, 54, 55, 57 e 58 restano nominati per compatibilità ma sollevano
`UnknownBodyError`; tutti i campi mancanti sono elencati
[nell'inventario dei modelli mancanti](../../methodology/missing-hypothetical-models.md).

```python
import libephemeris as ephem

jd_tt = ephem.julday(2024, 4, 8, 12.0) + ephem.deltat(ephem.julday(2024, 4, 8, 12.0))

cupido = ephem.calc_hypothetical_position(ephem.CUPIDO, jd_tt)
harrington = ephem.calc_hypothetical_position(ephem.HARRINGTON, jd_tt)
selena = ephem.calc_hypothetical_position(ephem.WHITE_MOON, jd_tt)
```

Nessun file di dati della distribuzione di riferimento viene incluso. Il
dataset contiene 12 righe orbitali trascritte indipendentemente; Selena è un
modello simbolico separato derivato dalla regola di moto uniforme e dai
checkpoint pubblicati. `HYPOTHETICAL_PROVENANCE` registra stato e confine
documentale di ogni ID.

---

## 13.3 Orbite fittizie personalizzate

Se hai bisogno di un corpo ipotetico non incluso nella libreria, puoi definire
la tua orbita usando il formato testuale documentato da LibEphemeris, composto da
nove campi di elementi orbitali. Usa valori provenienti da fonti che sei
autorizzato a utilizzare; i file di dati della distribuzione di riferimento non
sono richiesti né supportati come input inclusi.

### Caricare le orbite integrate verificate

La libreria include un file di orbite fittizie predefinite:

```python
import libephemeris as ephem

# Carica le orbite predefinite
orbite = ephem.load_bundled_fictitious_orbits()
print(f"Caricate {len(orbite)} orbite fittizie")

# Cerca uno dei corpi verificati per nome
corpo = ephem.get_orbital_body_by_name(orbite, "Harrington")
if corpo:
    print(f"Trovato: {corpo.name}")
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
- Gli ID 40–47, 50–53 e 56 hanno modelli ricostruiti indipendentemente e
  documentati; i sei ID non recuperati falliscono in modo esplicito
  conservando nomi e costanti
- Ogni modello storico supportato ha una derivazione pagina per pagina e un
  test di provenienza vincolante
- Le **orbite fittizie personalizzate** permettono di definire qualsiasi corpo ipotetico con i propri elementi orbitali
- Le **parti arabe** sono punti calcolati dalla formula ASC + Pianeta A − Pianeta B, con la Parte di Fortuna come la più importante

**Funzioni introdotte:**

- `calc_hypothetical_position(body_id, jd_tt)` — calcola un corpo ipotetico
  supportato (ID 40–47, 50–53 e 56)
- `load_bundled_fictitious_orbits()` — carica le orbite fittizie predefinite
- `parse_orbital_elements(filepath)` — carica orbite fittizie da un file personalizzato
- `get_orbital_body_by_name(elements, nome)` — cerca un corpo per nome nella lista di orbite
- `calc_orbital_position(elem, jd_tt)` — calcola la posizione di un corpo dato i suoi elementi orbitali
- `calc_all_arabic_parts(positions, jd=..., geo_lat=..., geo_lon=...)` — calcola le quattro parti arabe principali
