# Capitolo 2 — Misurare il tempo in astronomia

## Cosa imparerai

In questo capitolo scoprirai perché il tempo è molto più complicato di quanto sembri. Imparerai cos'è il Giorno Giuliano e perché gli astronomi lo usano, la differenza tra tempo universale (UT) e tempo terrestre (TT), cosa è il Delta-T e perché conta, e come funzionano il tempo siderale e il tempo locale.

---

## 2.1 Il Giorno Giuliano (JD)

### Il problema delle date

Quanti giorni sono passati tra il 15 marzo 44 a.C. (assassinio di Cesare) e il 20 luglio 1969 (allunaggio)? Per rispondere dovresti gestire il calendario giuliano, il calendario gregoriano, il fatto che nel 1582 sono "spariti" 10 giorni, gli anni bisestili con regole diverse nei due calendari, e il fatto che non esiste l'anno zero (dall'1 a.C. si passa direttamente all'1 d.C.).

Gli astronomi hanno risolto questo caos nel modo più semplice possibile: un contatore unico.

### La soluzione: un numero per ogni istante

Il **Giorno Giuliano** (Julian Day, JD) è un contatore continuo di giorni a partire dal mezzogiorno del 1° gennaio 4713 a.C. (calendario giuliano). Ogni istante della storia ha un unico numero JD:

| Evento | JD |
|--------|-----|
| Inizio del conteggio (1 gen 4713 a.C., mezzogiorno) | 0.0 |
| Nascita di Cristo (convenzionale) | ~1 721 424 |
| Riforma gregoriana (15 ott 1582) | 2 299 161 |
| J2000.0 (1 gen 2000, mezzogiorno TT) | 2 451 545.0 |
| Eclissi solare 8 apr 2024, mezzogiorno | 2 460 408.5 |

Un dettaglio importante: il giorno giuliano inizia a **mezzogiorno**, non a mezzanotte. Quindi JD 2460408**.0** è la mezzanotte dell'8 aprile 2024, e JD 2460408**.5** è il mezzogiorno. Questo è un residuo storico: gli astronomi lavoravano di notte, e far iniziare il giorno a mezzogiorno evitava di cambiare data nel mezzo di una sessione osservativa.

### 💻 Codice

```python
import libephemeris as ephem

# Data -> Giorno Giuliano
jd = ephem.julday(2024, 4, 8, 12.0)
print(f"JD = {jd}")

# Giorno Giuliano -> Data
anno, mese, giorno, ore = ephem.revjul(jd)
print(f"{giorno}/{mese}/{anno} ore {ore:.1f}")

# Calendario giuliano vs gregoriano
# Il 4 ottobre 1582 (giuliano) è il giorno prima del 15 ottobre 1582 (gregoriano)
jd_giuliano = ephem.julday(1582, 10, 4, 12.0, ephem.JUL_CAL)
jd_gregoriano = ephem.julday(1582, 10, 15, 12.0, ephem.GREG_CAL)
print(f"Differenza: {jd_gregoriano - jd_giuliano} giorni")

# Quanti giorni tra Cesare e l'allunaggio?
jd_cesare = ephem.julday(-43, 3, 15, 12.0, ephem.JUL_CAL)
jd_luna = ephem.julday(1969, 7, 20, 20.3)  # 20:17 UT
print(f"Giorni trascorsi: {jd_luna - jd_cesare:.0f}")
```

```
JD = 2460409.0
8/4/2024 ore 12.0
Differenza: 1.0 giorni
Giorni trascorsi: 734997
```

> **Nota sugli anni negativi**: LibEphemeris usa la convenzione astronomica dove l'anno 0 esiste. L'anno 1 a.C. = anno 0, il 2 a.C. = anno -1, e così via. Quindi il 44 a.C. = anno -43.

---

## 2.2 UT, TT e Delta-T

### Due tempi diversi per due scopi diversi

Il **Tempo Universale** (UT, più precisamente UT1) è basato sulla rotazione della Terra. Un giorno UT è il tempo che la Terra impiega per fare un giro completo rispetto al Sole medio. È il tempo dei nostri orologi, dei fusi orari, della vita quotidiana.

Il problema è che la Terra non è un orologio perfetto. La sua rotazione rallenta gradualmente (a causa delle maree lunari) e ha irregolarità imprevedibili. Un giorno UT non è sempre uguale.

Il **Tempo Terrestre** (TT) è un tempo perfettamente uniforme, basato su orologi atomici. Un secondo TT è sempre esattamente uguale a un altro. I pianeti si muovono secondo le leggi della fisica, che operano in tempo uniforme — per questo le efemeridi sono calcolate in TT.

### Delta-T: il ponte tra i due tempi

**Delta-T** (ΔT) è la differenza tra TT e UT:

    ΔT = TT − UT

Oggi Delta-T vale circa 69 secondi. Nel 1900 era circa 3 secondi. Nel 1800 circa 14 secondi. Andando indietro nei secoli, l'incertezza su Delta-T cresce rapidamente.

Perché conta? Immagina di calcolare la posizione della Luna per una data nel 1800. La Luna si muove di circa 0.5" al secondo. Se Delta-T ha un errore di 1 secondo, la posizione della Luna ha un errore di 0.5". Per date nel Medioevo, l'incertezza su Delta-T può essere di minuti, e l'errore sulla Luna di decine di arcminuti.

### Le funzioni `_ut` vs senza suffisso

LibEphemeris offre due versioni di molte funzioni:

- `calc_ut(jd_ut, ...)`: accetta un Giorno Giuliano in **UT** — converte internamente in TT aggiungendo Delta-T
- `calc(jd_tt, ...)`: accetta un Giorno Giuliano in **TT** — per chi vuole gestire Delta-T manualmente

Per la stragrande maggioranza degli usi, `calc_ut` è la scelta giusta: passi la data come la conosci (in UT, cioè il tempo "civile") e la libreria fa il resto.

### 💻 Codice

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Delta-T in giorni (la convenzione della libreria)
dt_giorni = ephem.deltat(jd)
dt_secondi = dt_giorni * 86400

print(f"Delta-T = {dt_secondi:.2f} secondi")

# Le due chiamate equivalenti:
pos_ut, _ = ephem.calc_ut(jd, ephem.MOON, 0)        # passa JD in UT
pos_tt, _ = ephem.calc(jd + dt_giorni, ephem.MOON, 0)  # passa JD in TT

print(f"Luna (via UT):  {pos_ut[0]:.6f}°")
print(f"Luna (via TT):  {pos_tt[0]:.6f}°")
```

```
Delta-T = 69.20 secondi
Luna (via UT):  15.429860°
Luna (via TT):  15.429860°
```

### 🌍 Vita reale

Per una carta natale del 1800, Delta-T era circa 14 secondi. La Luna si muove di circa 0.5" al secondo, quindi l'errore introdotto da un Delta-T sbagliato di 1 secondo è di ~0.5" sulla Luna — trascurabile per l'astrologia, ma non per il calcolo delle eclissi storiche. Per date antiche (prima del 1600), l'incertezza su Delta-T può essere di diversi minuti, rendendo impossibile calcolare posizioni lunari precise.

---

## 2.3 TAI, UTC e i secondi intercalari

Oltre a UT e TT, esistono altri due tempi importanti:

**UTC** (Coordinated Universal Time) è il tempo dei nostri orologi e dei computer. È basato su orologi atomici, ma viene periodicamente aggiustato con l'aggiunta di un **secondo intercalare** (leap second) per restare entro 0.9 secondi da UT1. Dal 1972 sono stati aggiunti 27 secondi intercalari.

**TAI** (International Atomic Time) è il tempo atomico puro, senza aggiustamenti. TAI corre continuamente senza mai fermarsi o saltare. La relazione è:

    TAI = UTC + N secondi    (dove N è il numero di leap second accumulati)
    TT  = TAI + 32.184 s     (costante fissa, per sempre)

Oggi (2024): TAI = UTC + 37 secondi, quindi TT = UTC + 69.184 secondi.

### 💻 Codice

```python
import libephemeris as ephem

# Convertire una data UTC precisa in JD TAI
# 8 aprile 2024, 18:23:45.123 UTC
jd_tai = ephem.utc_to_tai_jd(2024, 4, 8, 18, 23, 45.123)
print(f"JD (TAI) = {jd_tai:.8f}")

# Convertire da UTC a JD in TT e UT1
jd_tt, jd_ut = ephem.utc_to_jd(2024, 4, 8, 18, 23, 45.123)
print(f"JD (TT)  = {jd_tt:.8f}")
print(f"JD (UT1) = {jd_ut:.8f}")
print(f"TT - UT1 = {(jd_tt - jd_ut) * 86400:.2f} secondi")

# Riconvertire da JD TT a componenti UTC
anno, mese, giorno, ore, min, sec = ephem.jdet_to_utc(jd_tt)
print(f"UTC: {giorno}/{mese}/{anno} {ore}:{min}:{sec:.3f}")
```

```
JD (TAI) = 2460409.26692272
JD (TT)  = 2460409.26729522
JD (UT1) = 2460409.26649429
TT - UT1 = 69.20 secondi
UTC: 8/4/2024 18:23:45.123
```

### 🌍 Vita reale

Il tuo smartphone mostra l'ora in UTC (più il fuso orario locale). Ogni tanto — l'ultima volta il 31 dicembre 2016 — viene aggiunto un secondo intercalare: l'orologio segna 23:59:60 prima di passare a 00:00:00. Molti sistemi informatici non gestiscono bene questo caso, e nel corso degli anni i leap second hanno causato crash di server, rallentamenti di servizi cloud e bug nei sistemi di navigazione GPS.

---

## 2.4 Tempo siderale

### Il giorno siderale

Il **giorno siderale** è il tempo che la Terra impiega per fare un giro completo rispetto alle stelle (non rispetto al Sole). Dura **23 ore, 56 minuti e 4 secondi** — circa 4 minuti meno del giorno solare.

Perché la differenza? Mentre la Terra ruota su se stessa, si muove anche lungo la sua orbita attorno al Sole. Dopo un giro completo rispetto alle stelle, la Terra ha percorso un pezzetto della sua orbita, e deve ruotare ancora ~1° per "recuperare" e riallinearsi con il Sole. Quell'1° in più richiede circa 4 minuti.

Il **tempo siderale locale** misura quanta rotazione è avvenuta rispetto alle stelle. In pratica, indica quale parte della sfera celeste è al meridiano in quel momento. L'**ARMC** (Ascensione Retta del Medio Cielo) è l'orologio siderale locale, espresso in gradi (0°–360°).

Il tempo siderale è fondamentale per calcolare le case astrologiche: le case dipendono da come la sfera celeste è orientata rispetto all'orizzonte locale, e questo è esattamente ciò che il tempo siderale misura.

### 💻 Codice

```python
import libephemeris as ephem

# 8 aprile 2024, ore 21:00 UT
jd = ephem.julday(2024, 4, 8, 21.0)

# Tempo siderale di Greenwich (in ore)
st_greenwich = ephem.sidtime(jd)
print(f"Tempo siderale Greenwich: {st_greenwich:.4f} ore")

# Tempo siderale locale (Milano, lon 9.19° E)
st_milano = ephem.sidtime(jd, longitude=9.19)
print(f"Tempo siderale Milano:    {st_milano:.4f} ore")

# La differenza è proporzionale alla longitudine:
# 9.19° / 15 = 0.6127 ore ≈ 36 minuti e 46 secondi
```

```
Tempo siderale Greenwich: 10.1738 ore
Tempo siderale Milano:    10.7865 ore
```

### 🌍 Vita reale

Ecco perché le costellazioni "si alzano prima" ogni sera: ogni notte, alla stessa ora, la sfera celeste è ruotata di ~1° in più (i 4 minuti di differenza tra giorno solare e siderale). Dopo un mese, la differenza è di circa 2 ore. Dopo 6 mesi, il cielo notturno mostra costellazioni completamente diverse.

---

## 2.5 Tempo locale, fusi orari e LMT

### Il tempo solare locale

Prima dei fusi orari (introdotti nel 1884), ogni città usava il proprio **tempo solare locale**. Il mezzogiorno era quando il Sole attraversava il meridiano locale — e questo dipende dalla longitudine.

Il **Local Mean Time** (LMT) è il tempo solare medio per una data longitudine. La differenza rispetto a Greenwich è semplicemente la longitudine divisa per 15 (perché 360° / 24h = 15°/h):

    LMT = UT + longitudine / 15    (in ore)

A Milano (longitudine 9.19° E), il mezzogiorno solare medio arriva circa 37 minuti dopo quello di Greenwich.

Il **Local Apparent Time** (LAT) è il tempo della meridiana — il tempo solare *vero*, non medio. La differenza tra LAT e LMT è l'**equazione del tempo**: un'oscillazione di ±16 minuti durante l'anno, causata dall'eccentricità dell'orbita terrestre e dall'obliquità dell'eclittica.

### 💻 Codice

```python
import libephemeris as ephem

# 8 aprile 2024, ore 12:00 UT
jd = ephem.julday(2024, 4, 8, 12.0)

# Equazione del tempo (restituita in giorni)
eot_giorni = ephem.time_equ(jd)
eot_minuti = eot_giorni * 1440  # converti in minuti

print(f"Equazione del tempo: {eot_minuti:.2f} minuti")
# Se positivo: la meridiana è avanti rispetto all'orologio

# Conversione LMT -> LAT per Milano (lon 9.19° E)
jd_lmt = ephem.julday(2024, 4, 8, 12.0)  # mezzogiorno LMT
jd_lat = ephem.lmt_to_lat(jd_lmt, 9.19)

# La differenza è l'equazione del tempo
diff_minuti = (jd_lat - jd_lmt) * 1440
print(f"LAT - LMT = {diff_minuti:.2f} minuti")

# Conversione da UTC a ora locale
# Roma: UTC+2 in estate (CEST)
# utc_time_zone converte UTC -> ora locale aggiungendo l'offset
anno, mese, giorno, ore, min, sec = ephem.utc_time_zone(
    2024, 4, 8, 14, 30, 0.0, 2.0  # 14:30 UTC, offset +2 per CEST
)
print(f"Ora locale Roma: {ore}:{min}:{sec:.0f}")  # 16:30 CEST
```

```
Equazione del tempo: -0.42 minuti
LAT - LMT = -0.42 minuti
Ora locale Roma: 16:30:0
```

### 🌍 Vita reale

A Milano (longitudine 9.19° E, fuso CET = UTC+1), il mezzogiorno solare medio cade circa alle 12:23 CET in inverno (UTC+1) e alle 13:23 CEST in estate (UTC+2). Ma l'equazione del tempo aggiunge un'ulteriore oscillazione: a inizio novembre il mezzogiorno solare vero è circa 16 minuti prima del mezzogiorno solare medio, a metà febbraio circa 14 minuti dopo.

---

## 2.6 IERS e Delta-T osservato

### Il problema della previsione

Delta-T cambia nel tempo in modo non perfettamente prevedibile, perché dipende dalle irregolarità della rotazione terrestre. Per il passato recente (dal 1962 in poi), abbiamo **misurazioni precise** di Delta-T grazie all'IERS (International Earth Rotation and Reference Systems Service).

Per default, LibEphemeris usa il modello di Delta-T di Skyfield, che combina dati storici con previsioni. Ma per massima precisione su date recenti, puoi usare i dati IERS osservati:

### 💻 Codice

```python
import libephemeris as ephem

# Scarica i dati IERS (una tantum, ~1 MB)
ephem.download_delta_t_data()

# Abilita l'uso dei dati IERS osservati
ephem.set_iers_delta_t_enabled(True)

# Verifica che i dati siano disponibili per una data
jd = ephem.julday(2024, 4, 8, 12.0)
disponibile = ephem.is_iers_data_available(jd)
print(f"Dati IERS disponibili: {disponibile}")

# Confronto: Delta-T calcolato vs osservato
ephem.set_iers_delta_t_enabled(False)
dt_calcolato = ephem.deltat(jd) * 86400  # in secondi

ephem.set_iers_delta_t_enabled(True)
dt_osservato = ephem.deltat(jd) * 86400  # in secondi

print(f"Delta-T calcolato: {dt_calcolato:.3f} s")
print(f"Delta-T osservato: {dt_osservato:.3f} s")
print(f"Differenza:        {abs(dt_osservato - dt_calcolato):.3f} s")

# Per date recenti la differenza è piccola (< 0.1 s)
# Per date lontane nel futuro, non ci sono dati IERS
# e la libreria ricade sul modello calcolato
```

```
Dati IERS disponibili: False
Delta-T calcolato: 69.200 s
Delta-T osservato: 69.199 s
Differenza:        0.001 s
```

I dati IERS vengono aggiornati settimanalmente. La libreria può scaricarli automaticamente:

```python
# Abilita il download automatico dei dati IERS
ephem.set_iers_auto_download(True)

# Oppure scarica tutto manualmente
ephem.download_iers_finals()   # dati Earth Orientation
ephem.download_leap_seconds()  # tabella dei leap second
ephem.download_delta_t_data()  # serie storica Delta-T
```

---

## Riepilogo

- Il **Giorno Giuliano** (JD) è un contatore continuo di giorni, usato per evitare le complicazioni dei calendari. `julday()` converte una data in JD, `revjul()` fa il contrario.
- **UT** (Universal Time) è il tempo basato sulla rotazione terrestre — il tempo "civile". **TT** (Terrestrial Time) è il tempo uniforme degli orologi atomici.
- **Delta-T** = TT − UT, oggi ~69 secondi. Le funzioni `_ut` accettano JD in UT e convertono automaticamente.
- **UTC** è il tempo dei nostri orologi, con leap second. **TAI** è il tempo atomico puro. TT = TAI + 32.184 s.
- Il **tempo siderale** misura la rotazione rispetto alle stelle (~4 min più corto del giorno solare). È la base per il calcolo delle case astrologiche.
- Il **tempo locale** (LMT) dipende dalla longitudine. L'**equazione del tempo** è la differenza tra il tempo della meridiana e il tempo medio.
- I **dati IERS** forniscono Delta-T osservato per massima precisione su date recenti.

### Funzioni e costanti introdotte

| Funzione / Costante | Uso |
|---------------------|-----|
| `julday(anno, mese, giorno, ore, gregflag)` | Data → Giorno Giuliano |
| `revjul(jd, gregflag)` | Giorno Giuliano → data (anno, mese, giorno, ore) |
| `deltat(jd)` | Delta-T in giorni |
| `calc_ut(jd_ut, corpo, flag)` | Calcolo con JD in UT (converte internamente in TT) |
| `calc(jd_tt, corpo, flag)` | Calcolo con JD in TT |
| `utc_to_jd(anno, mese, giorno, ore, min, sec)` | UTC → JD in TT e UT1 |
| `utc_to_tai_jd(anno, mese, giorno, ore, min, sec)` | UTC → JD in TAI |
| `sidtime(jd, longitude)` | Tempo siderale locale in ore |
| `time_equ(jd)` | Equazione del tempo in giorni |
| `lmt_to_lat(jd_lmt, longitudine)` | Local Mean Time → Local Apparent Time |
| `utc_time_zone(anno, mese, giorno, ore, min, sec, offset)` | Fuso orario → UTC |
| `set_iers_delta_t_enabled(True/False)` | Abilita/disabilita dati IERS |
| `download_delta_t_data()` | Scarica i dati Delta-T osservati |
| `SE_GREG_CAL`, `SE_JUL_CAL` | Calendario gregoriano / giuliano |
