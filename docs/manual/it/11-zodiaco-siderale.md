# Capitolo 11 — Lo zodiaco siderale e le ayanamsha

## Cosa imparerai

In questo capitolo scoprirai perché esistono due zodiaci diversi (tropicale e siderale), cosa li ha fatti divergere nel corso dei millenni, cos'è l'ayanamsha, perché ci sono oltre 40 modi diversi di calcolarla, e come usare la libreria per lavorare con lo zodiaco siderale usato nell'astrologia vedica indiana.

---

## 11.1 Due zodiaci: tropicale e siderale

### Lo stesso cielo, due misure diverse

Nel Capitolo 1 abbiamo visto che lo zodiaco è una divisione dell'eclittica in 12 settori di 30° ciascuno. Ma da dove parte il conteggio? Qui nasce il problema — e una delle divisioni più profonde della storia dell'astrologia.

Lo **zodiaco tropicale** (usato in Occidente) fa partire il conteggio dall'**equinozio di primavera**: il punto dove il Sole attraversa l'equatore celeste salendo verso nord, attorno al 20 marzo. Questo punto è definito come 0° Ariete, indipendentemente da quali stelle si trovino in quella direzione.

Lo **zodiaco siderale** (usato in India e nell'astrologia vedica) fa partire il conteggio da un punto legato alle **stelle fisse**. Il punto zero è ancorato a una stella di riferimento (in genere Spica, fissata a 0° Bilancia = 180°).

Circa 2000 anni fa, i due zodiaci coincidevano quasi perfettamente. Ma da allora si sono spostati l'uno rispetto all'altro, e oggi la differenza è di circa **24 gradi**.

### Perché divergono: la precessione degli equinozi

La causa è la **precessione**: l'asse della Terra oscilla lentamente come una trottola, completando un giro in circa 26.000 anni. Questo fa sì che il punto dell'equinozio di primavera si sposti di circa 50 secondi d'arco all'anno *rispetto alle stelle*.

In pratica:
- Lo zodiaco tropicale "segue" l'equinozio, e quindi si sposta rispetto alle stelle
- Lo zodiaco siderale "segue" le stelle, e quindi si sposta rispetto all'equinozio

L'effetto pratico è notevole: se nel tuo tema natale tropicale (occidentale) il Sole è a 15° Ariete, nel tema siderale (vedico) sarà a circa 21° Pesci — quasi un segno intero di differenza.

```python
import libephemeris as ephem

# Confronto: posizione del Sole in tropicale e siderale
jd = ephem.julday(2024, 4, 8, 12.0)

# Posizione tropicale (default)
pos_trop, _ = ephem.calc_ut(jd, ephem.SUN, ephem.FLG_SPEED)

# Ayanamsha Lahiri
ephem.set_sid_mode(ephem.SIDM_LAHIRI)
ayan = ephem.get_ayanamsa_ut(jd)

# Posizione siderale = tropicale - ayanamsha
lon_sid = (pos_trop[0] - ayan) % 360

segni = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

s_trop = segni[int(pos_trop[0] / 30)]
g_trop = pos_trop[0] % 30

s_sid = segni[int(lon_sid / 30)]
g_sid = lon_sid % 30

print(f"Sole tropicale: {g_trop:.1f}° {s_trop}")
print(f"Sole siderale:  {g_sid:.1f}° {s_sid}")
print(f"Ayanamsha (Lahiri): {ayan:.4f}°")
```

```
Sole tropicale: 19.1° Ari
Sole siderale:  24.9° Psc
Ayanamsha (Lahiri): 24.1961°
```

---

## 11.2 Cos'è l'ayanamsha

L'**ayanamsha** (dal sanscrito *ayana* = solstizio, *amsha* = porzione) è semplicemente la differenza in gradi tra il punto zero tropicale e il punto zero siderale in un dato momento:

**longitudine siderale = longitudine tropicale − ayanamsha**

Oggi l'ayanamsha è circa 24°, e cresce di circa 50.3" all'anno (il tasso di precessione). Tra qualche millennio, i due zodiaci torneranno a coincidere — e poi divergeranno di nuovo nell'altra direzione.

### Calcolare l'ayanamsha

```python
import libephemeris as ephem

ephem.set_sid_mode(ephem.SIDM_LAHIRI)

# Ayanamsha in diverse epoche
for anno in [0, 500, 1000, 1500, 2000, 2024, 2100]:
    jd = ephem.julday(anno, 1, 1, 0.0)
    ayan = ephem.get_ayanamsa_ut(jd)
    print(f"Anno {anno:5d}: ayanamsha = {ayan:+7.2f}°")
```

```
Anno     0: ayanamsha = +356.04°
Anno   500: ayanamsha =   +2.97°
Anno  1000: ayanamsha =   +9.92°
Anno  1500: ayanamsha =  +16.88°
Anno  2000: ayanamsha =  +23.86°
Anno  2024: ayanamsha =  +24.19°
Anno  2100: ayanamsha =  +25.25°
```

Noterai che l'ayanamsha è negativa nell'antichità (lo zodiaco siderale era "avanti" rispetto al tropicale) e positiva oggi (il tropicale è "avanti").

### Usare il flag FLG_SIDEREAL

Invece di calcolare l'ayanamsha e sottrarla manualmente, puoi chiedere direttamente alla libreria di restituire le posizioni in coordinate siderali:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Imposta il sistema siderale
ephem.set_sid_mode(ephem.SIDM_LAHIRI)

# Calcola direttamente in siderale
pos_sid, _ = ephem.calc_ut(
    jd, ephem.SUN,
    ephem.FLG_SIDEREAL | ephem.FLG_SPEED
)

segni = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

segno = segni[int(pos_sid[0] / 30)]
gradi = pos_sid[0] % 30

print(f"Sole siderale (Lahiri): {gradi:.2f}° {segno}")
```

```
Sole siderale (Lahiri): 24.95° Psc
```

Il flag `FLG_SIDEREAL` (65536) può essere combinato con qualsiasi altro flag: `FLG_SPEED`, `FLG_EQUATORIAL`, ecc.

---

## 11.3 Le scuole di ayanamsha

### Il problema: nessuno è d'accordo

Se l'ayanamsha fosse un dato di fatto, non ci sarebbe discussione. Ma il problema è: **dove esattamente era il punto zero 2000 anni fa?** Non c'era un GPS celeste. Diversi astronomi e astrologi hanno proposto risposte leggermente diverse, basate su stelle diverse, testi antichi diversi, o calcoli diversi. La differenza tra le scuole più diffuse è di qualche grado — abbastanza da spostare un pianeta in un segno diverso.

La libreria supporta **47 sistemi di ayanamsha**. Ecco i più importanti:

### Ayanamsha indiane (le più usate)

**Lahiri** (`SE_SIDM_LAHIRI`, valore 1) è la più diffusa: è lo standard ufficiale del governo indiano, adottato nel 1955 dalla commissione del calendario indiano (N.C. Lahiri). Fissa Spica (Citra in sanscrito) a 180° di longitudine siderale. È usata dalla grande maggioranza degli astrologi vedici.

**Krishnamurti** (`SE_SIDM_KRISHNAMURTI`, valore 5) è usata nel sistema KP (Krishnamurti Paddhati), un metodo predittivo molto popolare in India del sud. Molto simile a Lahiri, con una differenza di pochi minuti d'arco.

**Raman** (`SE_SIDM_RAMAN`, valore 3) fu proposta da B.V. Raman, uno degli astrologi più influenti dell'India del XX secolo. Differisce da Lahiri di circa 1.5°.

**Varianti Lahiri**: ci sono anche `SE_SIDM_LAHIRI_1940` (43), `SE_SIDM_LAHIRI_VP285` (44), e `SE_SIDM_LAHIRI_ICRC` (46), che differiscono per pochi secondi d'arco e riflettono diverse interpretazioni dei dati originali.

### Ayanamsha occidentale siderale

**Fagan-Bradley** (`SE_SIDM_FAGAN_BRADLEY`, valore 0) fu sviluppata da Cyril Fagan e Donald Bradley negli anni '50 per l'astrologia siderale occidentale. Differisce da Lahiri di circa 1°. È poco usata al di fuori di un ristretto circolo di astrologi siderali occidentali.

### Ayanamsha basate su stelle vere

**True Citra** (`SE_SIDM_TRUE_CITRA`, valore 27) fissa la posizione *vera* (con moto proprio) di Spica a esattamente 180°. A differenza di Lahiri, che usa una formula polinomiale calcolata una volta per tutte, True Citra ricalcola la posizione reale di Spica ad ogni data, seguendone il moto proprio.

**True Revati** (`SE_SIDM_TRUE_REVATI`, valore 28) fissa la stella Revati (zeta Piscium) a 359°50'.

**True Pushya** (`SE_SIDM_TRUE_PUSHYA`, valore 29) fissa la stella Pushya (delta Cancri) a 106°.

### Ayanamsha galattiche

Per chi cerca un punto zero "cosmico", ci sono diversi sistemi basati sulla posizione del centro galattico o dell'equatore galattico:

**Galactic Center 0 Sag** (`SE_SIDM_GALCENT_0SAG`, valore 17) mette il centro galattico a 0° Sagittario.

### Ayanamsha babilonesi

Per la ricerca storica sull'astronomia babilonese: `SE_SIDM_BABYL_KUGLER1` (9), `SE_SIDM_BABYL_KUGLER2` (10), `SE_SIDM_BABYL_KUGLER3` (11), `SE_SIDM_BABYL_HUBER` (12), `SE_SIDM_BABYL_ETPSC` (13), `SE_SIDM_BABYL_BRITTON` (38).

### Confronto tra ayanamsha

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

sistemi = [
    (ephem.SIDM_LAHIRI,         "Lahiri"),
    (ephem.SIDM_FAGAN_BRADLEY,  "Fagan-Bradley"),
    (ephem.SIDM_RAMAN,          "Raman"),
    (ephem.SIDM_KRISHNAMURTI,   "Krishnamurti"),
    (ephem.SIDM_TRUE_CITRA,     "True Citra"),
    (ephem.SIDM_TRUE_REVATI,    "True Revati"),
    (ephem.SIDM_GALCENT_0SAG,   "Galactic 0 Sag"),
]

print(f"Ayanamsha al {8}/04/2024:\n")
for modo, nome in sistemi:
    ephem.set_sid_mode(modo)
    ayan = ephem.get_ayanamsa_ut(jd)
    print(f"  {nome:20s}  {ayan:.4f}°")

# Il nome di un sistema
nome = ephem.get_ayanamsa_name(ephem.SIDM_LAHIRI)
print(f"\nNome: {nome}")
```

```
Ayanamsha al 8/04/2024:

  Lahiri                24.1961°
  Fagan-Bradley         25.0793°
  Raman                 22.7498°
  Krishnamurti          24.0993°
  True Citra            24.1843°
  True Revati           20.3779°
  Galactic 0 Sag        27.1850°

Nome: Lahiri
```

---

## 11.4 Ayanamsha personalizzata

Se nessuno dei 47 sistemi predefiniti soddisfa le tue esigenze, puoi definire la tua ayanamsha con `SE_SIDM_USER`:

```python
import libephemeris as ephem

# Definisci un'ayanamsha custom:
# "Al 1° gennaio 2000 (J2000.0), l'ayanamsha vale 23.5°"
J2000 = 2451545.0

ephem.set_sid_mode(
    ephem.SIDM_USER,
    t0=J2000,       # epoca di riferimento
    ayan_t0=23.5     # valore dell'ayanamsha a quell'epoca
)

# Ora calcola l'ayanamsha per qualsiasi data
jd = ephem.julday(2024, 4, 8, 12.0)
ayan = ephem.get_ayanamsa_ut(jd)
print(f"Ayanamsha custom al 2024: {ayan:.4f}°")

# E usa FLG_SIDEREAL per le posizioni
pos, _ = ephem.calc_ut(jd, ephem.SUN,
    ephem.FLG_SIDEREAL | ephem.FLG_SPEED)
print(f"Sole siderale custom: {pos[0]:.4f}°")
```

```
Ayanamsha custom al 2024: 23.8390°
Sole siderale custom: 355.3029°
```

La libreria usa i polinomi di precessione IAU 2006 per calcolare come l'ayanamsha cambia nel tempo a partire dal valore che hai specificato.

---

## 11.5 Calcoli pratici in siderale

### Tema natale vedico

Ecco come calcolare un tema natale completo in siderale Lahiri — il formato usato nell'astrologia vedica (Jyotish):

```python
import libephemeris as ephem

# Data: 15 agosto 1947, 00:00 IST — Indipendenza dell'India
# IST = UT + 5:30, quindi UT = 18:30 del 14 agosto
jd = ephem.julday(1947, 8, 14, 18.5)
lat, lon = 28.6139, 77.2090  # Nuova Delhi

# Imposta siderale Lahiri
ephem.set_sid_mode(ephem.SIDM_LAHIRI)

segni = ["Mesha", "Vrishabha", "Mithuna", "Karka",
         "Simha", "Kanya", "Tula", "Vrischika",
         "Dhanus", "Makara", "Kumbha", "Meena"]

# Pianeti in siderale
pianeti = [
    (ephem.SUN,      "Surya  "),
    (ephem.MOON,     "Chandra"),
    (ephem.MARS,     "Mangal "),
    (ephem.MERCURY,  "Budha  "),
    (ephem.JUPITER,  "Guru   "),
    (ephem.VENUS,    "Shukra "),
    (ephem.SATURN,   "Shani  "),
    (ephem.MEAN_NODE,"Rahu   "),
]

print("Tema natale siderale (Lahiri):\n")
for body_id, nome in pianeti:
    pos, _ = ephem.calc_ut(jd, body_id,
        ephem.FLG_SIDEREAL | ephem.FLG_SPEED)
    segno = segni[int(pos[0] / 30)]
    gradi = pos[0] % 30
    print(f"  {nome}  {gradi:5.1f}°  {segno}")

# Ketu = opposto a Rahu
rahu_pos, _ = ephem.calc_ut(jd, ephem.MEAN_NODE,
    ephem.FLG_SIDEREAL | ephem.FLG_SPEED)
ketu_lon = (rahu_pos[0] + 180) % 360
segno_k = segni[int(ketu_lon / 30)]
gradi_k = ketu_lon % 30
print(f"  Ketu    {gradi_k:5.1f}°  {segno_k}")
```

```
Tema natale siderale (Lahiri):

  Surya    28.0°  Karka
  Chandra   4.0°  Karka
  Mangal    7.5°  Mithuna
  Budha    13.7°  Karka
  Guru     25.9°  Tula
  Shukra   22.6°  Karka
  Shani    20.5°  Karka
  Rahu      5.1°  Vrishabha
  Ketu      5.1°  Vrischika
```

### I Nakshatra

Nell'astrologia vedica, lo zodiaco è diviso anche in **27 Nakshatra** (asterismi lunari) di 13°20' ciascuno. Ogni Nakshatra ha un nome, un significato e un pianeta reggente. La Luna impiega circa un giorno per attraversare ciascuno.

```python
import libephemeris as ephem

nakshatra_nomi = [
    "Ashwini", "Bharani", "Krittika", "Rohini", "Mrigashira",
    "Ardra", "Punarvasu", "Pushya", "Ashlesha", "Magha",
    "Purva Phalguni", "Uttara Phalguni", "Hasta", "Chitra",
    "Swati", "Vishakha", "Anuradha", "Jyeshtha", "Mula",
    "Purva Ashadha", "Uttara Ashadha", "Shravana", "Dhanishtha",
    "Shatabhisha", "Purva Bhadrapada", "Uttara Bhadrapada", "Revati"
]

# Reggenti nella sequenza Vimshottari Dasha
reggenti = [
    "Ketu", "Venere", "Sole", "Luna", "Marte",
    "Rahu", "Giove", "Saturno", "Mercurio",
    "Ketu", "Venere", "Sole", "Luna", "Marte",
    "Rahu", "Giove", "Saturno", "Mercurio",
    "Ketu", "Venere", "Sole", "Luna", "Marte",
    "Rahu", "Giove", "Saturno", "Mercurio"
]

jd = ephem.julday(2024, 4, 8, 12.0)
ephem.set_sid_mode(ephem.SIDM_LAHIRI)

# Nakshatra della Luna
moon, _ = ephem.calc_ut(jd, ephem.MOON,
    ephem.FLG_SIDEREAL | ephem.FLG_SPEED)

nak_num = int(moon[0] / (360 / 27))  # 13°20' per Nakshatra
nak_pos = moon[0] % (360 / 27)       # posizione nel Nakshatra

# Pada (quarto) — ogni Nakshatra ha 4 pada di 3°20'
pada = int(nak_pos / (360 / 108)) + 1  # 1-4

print(f"Luna siderale: {moon[0]:.2f}°")
print(f"Nakshatra: {nakshatra_nomi[nak_num]} (n. {nak_num + 1})")
print(f"Pada: {pada}")
print(f"Reggente: {reggenti[nak_num]}")
```

```
Luna siderale: 351.24°
Nakshatra: Revati (n. 27)
Pada: 2
Reggente: Mercurio
```

Il Nakshatra della Luna alla nascita è fondamentale nell'astrologia vedica: determina il **Dasha** (il sistema di periodi planetari che governa la vita della persona).

---

## 11.6 La versione estesa: `get_ayanamsa_ex_ut`

Per calcoli avanzati dove hai bisogno di specificare flag aggiuntivi o vuoi il flag di ritorno, usa `get_ayanamsa_ex_ut`:

```python
import libephemeris as ephem

ephem.set_sid_mode(ephem.SIDM_LAHIRI)

jd = ephem.julday(2024, 4, 8, 12.0)

# Versione semplice
ayan_simple = ephem.get_ayanamsa_ut(jd)

# Versione estesa (restituisce anche il flag)
retflag, ayan_ex = ephem.get_ayanamsa_ex_ut(jd, ephem.FLG_SWIEPH)

print(f"Ayanamsha (semplice): {ayan_simple:.6f}°")
print(f"Ayanamsha (estesa):   {ayan_ex:.6f}°")
```

```
Ayanamsha (semplice): 24.196111°
Ayanamsha (estesa):   24.196111°
```

---

## Riepilogo

In questo capitolo abbiamo esplorato lo zodiaco siderale, fondamentale per l'astrologia vedica.

**Concetti chiave:**

- Esistono due zodiaci: il **tropicale** (punto zero = equinozio) e il **siderale** (punto zero legato alle stelle). Oggi differiscono di circa 24°
- La divergenza è causata dalla **precessione degli equinozi**: l'asse terrestre oscilla con un periodo di ~26.000 anni, spostando il punto vernale rispetto alle stelle
- L'**ayanamsha** è la differenza in gradi tra i due zodiaci: longitudine siderale = longitudine tropicale − ayanamsha
- Esistono oltre **47 sistemi di ayanamsha** perché non c'è accordo su dove fosse il punto zero 2000 anni fa
- **Lahiri** è la più usata (standard del governo indiano), **Fagan-Bradley** per l'astrologia siderale occidentale, **True Citra** usa la posizione reale di Spica
- Nell'astrologia vedica, i **27 Nakshatra** dividono lo zodiaco siderale in settori di 13°20', ognuno con un significato e un pianeta reggente

**Funzioni introdotte:**

- `set_sid_mode(sid_mode, t0=0.0, ayan_t0=0.0)` — imposta il sistema di ayanamsha. Usa `SE_SIDM_USER` con `t0` e `ayan_t0` per un'ayanamsha personalizzata
- `get_ayanamsa_ut(jd)` — restituisce l'ayanamsha in gradi per una data in UT
- `get_ayanamsa_ex_ut(jd, flags)` — versione estesa che restituisce anche il flag di ritorno
- `get_ayanamsa_name(sid_mode)` — restituisce il nome leggibile di un sistema di ayanamsha (es. `"Lahiri"`)
- `FLG_SIDEREAL` (65536) — flag da aggiungere a `calc_ut` per ottenere posizioni direttamente in coordinate siderali
