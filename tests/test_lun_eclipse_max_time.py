# SPDX-License-Identifier: AGPL-3.0-only
"""Prove sintetiche del raffinamento del massimo di eclissi lunare."""

from __future__ import annotations

import math
from collections.abc import Callable

import pytest

from libephemeris import eclipse
from libephemeris.constants import FLG_EQUATORIAL, FLG_MOSEPH, FLG_SPEED, FLG_XYZ
from libephemeris.constants import MOON, SUN
from libephemeris.exceptions import ConvergenceError, EphemerisRangeError


def _installa_stati_sintetici(
    monkeypatch: pytest.MonkeyPatch,
    seme: float,
    angolo: Callable[[float], float],
    distanza_sole: Callable[[float], float] = lambda _t: 1.0,
) -> list[tuple[float, int, int]]:
    """Installa stati apparenti coerenti e restituisce il registro chiamate."""
    chiamate: list[tuple[float, int, int]] = []
    distanza_luna = 0.00257

    def calc_ut_finto(t: float, corpo: int, flags: int):
        chiamate.append((t, corpo, flags))
        if corpo == MOON:
            posizione = (distanza_luna, 0.0, 0.0)
        elif corpo == SUN:
            alfa = angolo(t - seme)
            raggio = distanza_sole(t - seme)
            posizione = (
                distanza_luna - raggio * math.cos(alfa),
                raggio * math.sin(alfa),
                0.0,
            )
        else:  # pragma: no cover - rende esplicito il contratto del doppio stato
            raise AssertionError(f"corpo inatteso: {corpo}")
        return ((*posizione, 0.0, 0.0, 0.0), flags)

    monkeypatch.setattr(eclipse, "calc_ut", calc_ut_finto)
    return chiamate


def test_massimizza_la_sovrapposizione_apparente_nella_finestra_del_seme(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Il raggio solare variabile sposta il massimo dal minimo angolare."""
    seme = 2460000.0
    curvatura = 0.1
    pendenza_distanza = -0.5
    chiamate = _installa_stati_sintetici(
        monkeypatch,
        seme,
        lambda dt: curvatura * dt * dt,
        lambda dt: 1.0 + pendenza_distanza * dt,
    )

    ottenuto = eclipse._lun_eclipse_max_time(seme, FLG_MOSEPH | FLG_SPEED)

    # L'obiettivo pubblicato deve crescere fino al risultato e decrescere dopo;
    # ciò dimostra anche lo spostamento rispetto al minimo angolare posto al seme.
    assert isinstance(ottenuto, float)
    assert ottenuto > seme + 0.01
    valori = []
    for scarto in (-1.0e-5, 0.0, 1.0e-5):
        dt = ottenuto + scarto - seme
        distanza = 1.0 + pendenza_distanza * dt
        valori.append(math.asin(eclipse._ECL_RSUN_AU / distanza) - curvatura * dt * dt)
    assert valori[1] > valori[0]
    assert valori[1] > valori[2]
    assert chiamate
    assert all(seme - 0.3 <= t <= seme + 0.3 for t, _corpo, _flags in chiamate)
    assert {corpo for _t, corpo, _flags in chiamate} == {MOON, SUN}
    flags_stato = FLG_MOSEPH | FLG_EQUATORIAL | FLG_XYZ
    assert all(flags == flags_stato for _t, _corpo, flags in chiamate)


@pytest.mark.parametrize("scarto", [-0.3, 0.3])
def test_include_entrambi_gli_estremi_della_finestra_chiusa(
    monkeypatch: pytest.MonkeyPatch, scarto: float
) -> None:
    """Un massimo unico posto a un estremo viene restituito esattamente."""
    seme = 2460100.0
    estremo = seme + scarto
    _installa_stati_sintetici(
        monkeypatch,
        seme,
        lambda dt: 0.4 * abs((seme + dt) - estremo),
    )

    ottenuto = eclipse._lun_eclipse_max_time(seme)

    assert ottenuto == estremo


@pytest.mark.parametrize("obiettivo", [-0.295, 0.295])
@pytest.mark.parametrize("esponente", [1, 2])
def test_trova_un_massimo_nel_segmento_adiacente_a_un_estremo(
    monkeypatch: pytest.MonkeyPatch, obiettivo: float, esponente: int
) -> None:
    """Il primo e l'ultimo segmento non vengono scambiati per gli estremi."""
    seme = 2460000.0
    _installa_stati_sintetici(
        monkeypatch,
        seme,
        lambda dt: 0.4 * abs(dt - obiettivo) ** esponente,
    )

    ottenuto = eclipse._lun_eclipse_max_time(seme)

    assert ottenuto == pytest.approx(seme + obiettivo, abs=2.0e-8)


@pytest.mark.parametrize(
    "obiettivo",
    [
        segno * (0.275 + indice * 1.0e-9)
        for segno in (-1.0, 1.0)
        for indice in range(-20, 21)
    ],
)
def test_non_duplica_il_picco_unico_nella_transizione_del_bordo(
    monkeypatch: pytest.MonkeyPatch, obiettivo: float
) -> None:
    """Le regioni interna e di bordo descrivono una sola regione candidata."""
    seme = 2460000.0
    _installa_stati_sintetici(
        monkeypatch,
        seme,
        lambda dt: math.hypot(1.0e-5, 0.01 * (dt - obiettivo)),
    )

    ottenuto = eclipse._lun_eclipse_max_time(seme)

    assert ottenuto == pytest.approx(seme + obiettivo, abs=2.0e-9)


def test_rifiuta_due_massimi_separati(monkeypatch: pytest.MonkeyPatch) -> None:
    """Una valle campionata conserva due regioni di massimo indipendenti."""
    seme = 2460000.0
    _installa_stati_sintetici(
        monkeypatch,
        seme,
        lambda dt: 0.2 * min(abs(dt + 0.1), abs(dt - 0.1)),
    )

    with pytest.raises(ConvergenceError, match="non è unico"):
        eclipse._lun_eclipse_max_time(seme)


def test_propaga_la_mancanza_di_copertura_a_un_estremo(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """La copertura richiesta comprende l'intera finestra, estremi inclusi."""
    seme = 2460200.0
    limite = seme - 0.3
    errore = EphemerisRangeError("copertura sintetica assente", requested_jd=limite)

    def calc_ut_finto(t: float, corpo: int, flags: int):
        if t == limite and corpo == SUN:
            raise errore
        if corpo == MOON:
            posizione = (0.00257, 0.0, 0.0)
        else:
            posizione = (-0.99743, t - seme, 0.0)
        return ((*posizione, 0.0, 0.0, 0.0), flags)

    monkeypatch.setattr(eclipse, "calc_ut", calc_ut_finto)

    with pytest.raises(EphemerisRangeError) as catturato:
        eclipse._lun_eclipse_max_time(seme)

    assert catturato.value is errore


def test_rifiuta_massimo_non_unico_e_semi_non_finiti(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Un obiettivo costante non autorizza la scelta di un istante arbitrario."""
    seme = 2460300.0
    _installa_stati_sintetici(monkeypatch, seme, lambda _dt: 0.0)

    with pytest.raises(ConvergenceError, match="non è unico"):
        eclipse._lun_eclipse_max_time(seme)

    for seme_non_finito in (math.nan, math.inf, -math.inf):
        with pytest.raises(ValueError, match="finito"):
            eclipse._lun_eclipse_max_time(seme_non_finito)
