# Piano di riordino documentazione

> **Stato: completato e superato.** Questo file conserva il piano storico
> precedente al riordino. L'inventario, i percorsi e i collegamenti proposti
> qui sotto non descrivono la struttura corrente e non sono un indice attivo;
> usare [`docs/README.md`](../docs/README.md).

## Stato attuale

La cartella `docs/` contiene 40+ file markdown distribuiti senza una logica chiara.
Ci sono file legacy in root (`AYANAMSHA.md`, `PRECISION.md`, etc.) che duplicano
contenuti gia' presenti nelle sottocartelle, e file nuovi (`horizons-backend.md`,
`flags.md`, `getting-started.md`) che non sono ancora integrati nell'indice.

### Inventario attuale

```
docs/
  README.md                     # Indice — INCOMPLETO (manca Horizons, flags, getting-started)
  getting-started.md            # NUOVO (v1.0.0a3)
  flags.md                      # NUOVO (v1.0.0a3)
  horizons-backend.md           # NUOVO (v1.0.0a3)
  divergences.md                # 311 righe — divergenze note
  PRECISION.md                  # 94 righe — DUPLICATO? c'e' anche reference/precision.md (913 righe)
  AYANAMSHA.md                  # 156 righe — DUPLICATO di reference/ayanamsha.md (997 righe)
  INTERPOLATED_APOGEE.md        # 96 righe — DUPLICATO di methodology/interpolated-apogee.md (399 righe)
  PYERFA_BENEFITS.md            # 81 righe — DUPLICATO di methodology/pyerfa-integration.md (156 righe)
  REBOUND_BENEFITS.md           # 101 righe — DUPLICATO di methodology/rebound-integration.md (180 righe)
  TRUE_LILITH_METHODS.md        # 68 righe — DUPLICATO di methodology/true-lilith.md (192 righe)
  guides/
    migration-guide.md          # 441 righe — OK
    precision-tuning.md         # 510 righe — OK
  reference/
    precision.md                # 913 righe — versione completa
    house-systems.md            # 710 righe — OK
    ayanamsha.md                # 997 righe — versione completa
    swisseph-comparison.md      # 527 righe — OK
    se-bug-sidereal-j2000-nodes.md  # 190 righe — OK
  methodology/
    overview.md                 # 164 righe — OK
    lunar-apsides.md            # 69 righe — OK
    interpolated-apogee.md      # 399 righe — versione completa
    interpolated-perigee.md     # 229 righe — OK
    true-lilith.md              # 192 righe — versione completa
    planet-centers-spk.md       # 247 righe — OK
    pyerfa-integration.md       # 156 righe — versione completa
    rebound-integration.md      # 180 righe — versione completa
  leb/
    guide.md                    # 2075 righe — OK (ha sezione LEB2)
    algorithms.md               # 1021 righe — OK
    testing.md                  # 475 righe — OK
  development/
    architecture-overview.md    # 382 righe — OK
    testing.md                  # 386 righe — OK
    roadmap.md                  # 197 righe — OK
    keplerian-improvements.md   # 562 righe — OK
    full-range-coverage.md      # 206 righe — OK
    precision-history.md        # 731 righe — OK
  manual/                       # IT manual (15 capitoli) — NON TOCCARE
  manual-eng/                   # EN manual (15 capitoli) — NON TOCCARE
```

## Problemi

1. **6 file DUPLICATI** in `docs/` root che sono versioni corte di file gia' presenti
   nelle sottocartelle (AYANAMSHA.md, PRECISION.md, INTERPOLATED_APOGEE.md,
   PYERFA_BENEFITS.md, REBOUND_BENEFITS.md, TRUE_LILITH_METHODS.md)

2. **Indice (`docs/README.md`) incompleto** — manca: getting-started, flags,
   horizons-backend, divergences

3. **File nuovi sparsi** in root (`getting-started.md`, `flags.md`,
   `horizons-backend.md`) invece di essere nelle sottocartelle appropriate

4. **Nessuna sezione "Architecture"** nell'indice per Horizons e i backend

## Piano di riordino

### Fase 1: Eliminare i duplicati

| File root (ELIMINARE) | File completo (MANTENERE) | Azione |
|------------------------|---------------------------|--------|
| `docs/PRECISION.md` | `docs/reference/precision.md` | Eliminare, redirect in README |
| `docs/AYANAMSHA.md` | `docs/reference/ayanamsha.md` | Eliminare |
| `docs/INTERPOLATED_APOGEE.md` | `docs/methodology/interpolated-apogee.md` | Eliminare |
| `docs/PYERFA_BENEFITS.md` | `docs/methodology/pyerfa-integration.md` | Eliminare |
| `docs/REBOUND_BENEFITS.md` | `docs/methodology/rebound-integration.md` | Eliminare |
| `docs/TRUE_LILITH_METHODS.md` | `docs/methodology/true-lilith.md` | Eliminare |

### Fase 2: Spostare file nelle sottocartelle corrette

| File attuale | Destinazione | Motivo |
|---|---|---|
| `docs/getting-started.md` | `docs/guides/getting-started.md` | Coerenza con guides/ |
| `docs/flags.md` | `docs/reference/flags.md` | E' un riferimento tecnico |
| `docs/horizons-backend.md` | `docs/architecture/horizons-backend.md` | Architettura backend |
| `docs/divergences.md` | `docs/reference/divergences.md` | Riferimento tecnico |

### Fase 3: Creare sottocartella `architecture/`

```
docs/architecture/
  horizons-backend.md           # (spostato da docs/)
  calculation-backends.md       # NUOVO — panoramica dei 3 backend (auto/leb/horizons/skyfield)
```

### Fase 4: Aggiornare `docs/README.md`

Struttura proposta:

```markdown
# LibEphemeris Documentation

## Getting Started
- [Getting Started](guides/getting-started.md)
- [Migration from PySwissEph](guides/migration-guide.md)
- [Precision Tuning](guides/precision-tuning.md)

## Architecture
- [Calculation Backends](architecture/calculation-backends.md)
- [Horizons API Backend](architecture/horizons-backend.md)
- [LEB Binary Ephemeris](leb/guide.md)
- [LEB Algorithms](leb/algorithms.md)
- [Architecture Overview](development/architecture-overview.md)

## Reference
- [Flag Reference](reference/flags.md)
- [Precision Report](reference/precision.md)
- [Known Divergences](reference/divergences.md)
- [House Systems](reference/house-systems.md)
- [Ayanamsha Modes](reference/ayanamsha.md)
- [SE Bug: Sidereal J2000 Nodes](reference/se-bug-sidereal-j2000-nodes.md)
- [Comparison Report](reference/swisseph-comparison.md)

## Methodology
- [Overview](methodology/overview.md)
- [Planet Centers](methodology/planet-centers-spk.md)
- [Lunar Apsides](methodology/lunar-apsides.md)
- [Interpolated Apogee](methodology/interpolated-apogee.md)
- [Interpolated Perigee](methodology/interpolated-perigee.md)
- [True Lilith](methodology/true-lilith.md)
- [pyerfa Integration](methodology/pyerfa-integration.md)
- [REBOUND Integration](methodology/rebound-integration.md)

## LEB Binary Ephemeris
- [Technical Guide](leb/guide.md)
- [Algorithms & Theory](leb/algorithms.md)
- [Comparison Testing](leb/testing.md)

## Development
- [Testing](development/testing.md)
- [Roadmap](development/roadmap.md)
- [Architecture Overview](development/architecture-overview.md)
- [Precision History](development/precision-history.md)
- [Keplerian Improvements](development/keplerian-improvements.md)
- [Full Range Coverage](development/full-range-coverage.md)

## Manuals
- [Manuale (IT)](manual/) — 15 capitoli
- [Manual (EN)](manual-eng/) — 15 chapters
```

### Fase 5: Aggiornare i link nel README.md principale

Aggiornare tutti i path dei link dopo gli spostamenti.

### Riepilogo impatto

| Azione | File |
|--------|------|
| Eliminare | 6 file duplicati |
| Spostare | 4 file |
| Creare | 1 nuovo (`architecture/calculation-backends.md`) |
| Aggiornare | `docs/README.md`, `README.md` (link) |
| NON TOCCARE | manual/, manual-eng/, tutti i file nelle sottocartelle |
