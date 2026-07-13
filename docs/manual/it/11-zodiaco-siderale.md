# Capitolo 11 — Lo zodiaco siderale e l'ayanamsha

## Cosa imparerai

Questo capitolo spiega la distinzione tra zodiaco tropicale e siderale, i modi
predefiniti implementati da LibEphemeris, lo stato dell'audit delle fonti e come
fornire una definizione propria documentata.

## 11.1 Due origini per la longitudine

Lo zodiaco tropicale usa come origine l'equinozio vernale mobile. Uno zodiaco
siderale usa invece un'origine legata a una stella o a un altro riferimento
celeste. La precessione modifica nel tempo l'angolo fra le due origini.

L'angolo sottratto alla longitudine tropicale è l'**ayanamsha**:

```text
longitudine siderale = longitudine tropicale - ayanamsha
```

Non esiste un'unica origine siderale universale. Le tradizioni storiche scelgono
stelle, epoche o costruzioni geometriche diverse. LibEphemeris richiede quindi
una definizione indipendente e riproducibile per ogni modo calcolato.

## 11.2 Classificazione dei modi

Tutte le costanti canoniche `SIDM_*` sono esportate e ogni ID base predefinito
0--46 viene calcolato senza warning di fallback. `SIDM_USER` resta disponibile
per epoca e offset definiti dall'utente.

I modi stellari propagano l'astrometria di cataloghi con provenienza indipendente
attraverso la pipeline Skyfield/ERFA. I modi galattici usano il frame IAU e la
costruzione geometrica pubblicata. I modi a epoca fissa sono proiezioni di frame,
non tabelle di longitudini prive di fonte.

L'audit numerico delle fonti è separato dal supporto runtime. La
[tabella autorevole](../../reference/ayanamsha.md) distingue le definizioni
verificate in pubblicazioni primarie dagli anchor di compatibilità ripristinati
che richiedono ancora revisione del proprietario.

## 11.3 Calcolare una posizione siderale

Scegli un modo e richiedi `FLG_SIDEREAL`:

```python
import libephemeris as ephem

jd_ut = ephem.julday(2024, 4, 8, 12.0)
ephem.set_sid_mode(ephem.SIDM_TRUE_CITRA)

ayanamsha = ephem.get_ayanamsa_ut(jd_ut)
posizione, retflag = ephem.calc_ut(
    jd_ut,
    ephem.SUN,
    ephem.FLG_SIDEREAL | ephem.FLG_SPEED,
)
```

`FLG_SIDEREAL` può essere combinato con i flag di coordinate e velocità
supportati. Il modo selezionato è globale nell'API classica; usa
`EphemerisContext` quando calcoli paralleli richiedono configurazioni diverse.

```python
ctx = ephem.EphemerisContext()
ctx.set_sid_mode(ephem.SIDM_J2000)
posizione, retflag = ctx.calc_ut(
    jd_ut,
    ephem.MOON,
    ephem.FLG_SIDEREAL | ephem.FLG_SPEED,
)
```

## 11.4 Ayanamsha definita dall'utente

Usa `SIDM_USER` quando disponi di un'epoca e di un offset provenienti da una
fonte indipendente:

```python
import libephemeris as ephem

ephem.set_sid_mode(
    ephem.SIDM_USER,
    t0=jd_riferimento_della_fonte,
    ayan_t0=offset_in_gradi_della_fonte,
)
```

LibEphemeris propaga quel punto zero con la propria implementazione indipendente
della precessione a lungo termine. L'utente è responsabile della provenienza e
della scala temporale dei valori. È preferibile passare esplicitamente l'epoca.

## 11.5 API per ayanamsha media e vera

`get_ayanamsa_ut(jd_ut)` e `get_ayanamsa(jd_tt)` restituiscono l'offset medio.
Le varianti estese restituiscono anche i flag e rispettano la convenzione di
nutazione richiesta:

```python
retflag, valore = ephem.get_ayanamsa_ex_ut(jd_ut, ephem.FLG_SWIEPH)
```

Le pipeline planetaria, stellare, Skyfield, Horizons e LEB usano lo stesso
dispatcher siderale. La scelta del backend non deve cambiare implicitamente
l'anchor.

## 11.6 Modi a epoca fissa

`SIDM_J2000`, `SIDM_J1900` e `SIDM_B1950` richiedono coordinate nel frame medio
dell'epoca indicata. L'output eclittico viene proiettato sull'eclittica e
l'equinozio medi di quell'epoca; l'output equatoriale usa il corrispondente frame
equatoriale medio. Non sono ayanamsha astrologiche tradizionali.

## Riepilogo

- La longitudine siderale sottrae un'ayanamsha dalla longitudine tropicale.
- Ogni modo base `SIDM_*` predefinito 0--46 viene calcolato senza fallback.
- Lo stato dell'audit delle fonti di ciascun anchor storico è documentato a parte.
- Usa `SIDM_USER` per una definizione tradizionale lecita e citabile.
- L'output di implementazioni esterne è usato solo per confronti effimeri,
  mai come sorgente di anchor.

Consulta [Modi siderali](../../reference/ayanamsha.md) per la lista autorevole e
le note di provenienza.
