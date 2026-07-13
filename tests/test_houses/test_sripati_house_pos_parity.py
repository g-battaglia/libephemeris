"""Reference-free invariants for public Sripati ``house_pos`` handling."""

from __future__ import annotations

import pytest

import libephemeris as ephem


@pytest.mark.unit
@pytest.mark.parametrize(
    "armc,geolat,eps,lon,lat",
    [
        (0.0, 23.4, 23.4393, 50.0, 0.0),
        (123.0, 66.0, 23.4393, 169.0, 2.0),
        (271.0, 71.0, 23.4393, 335.0, -3.0),
        (34.0, -63.0, 23.4393, 68.0, 1.0),
    ],
)
def test_sripati_position_is_periodic_and_bounded(
    armc: float, geolat: float, eps: float, lon: float, lat: float
) -> None:
    position = ephem.house_pos(armc, geolat, eps, (lon, lat), b"S")
    wrapped = ephem.house_pos(armc, geolat, eps, (lon + 360.0, lat), b"S")
    assert isinstance(position, float)
    # A decimal house position in house 12 is in [12, 13), not capped at
    # the integer house number 12.
    assert 1.0 <= position < 13.0
    assert wrapped == pytest.approx(position, abs=1e-12)


@pytest.mark.unit
def test_sripati_calling_forms_agree() -> None:
    tuple_form = ephem.house_pos(0.0, 23.4, 23.4393, (75.0, 0.0), b"S")
    scalar_form = ephem.house_pos(0.0, 23.4, 23.4393, ord("S"), 75.0, 0.0)
    assert scalar_form == tuple_form


@pytest.mark.unit
def test_sripati_longitude_progression_is_locally_continuous() -> None:
    before = ephem.house_pos(123.0, 40.0, 23.4393, (169.0, 0.0), b"S")
    after = ephem.house_pos(123.0, 40.0, 23.4393, (169.001, 0.0), b"S")
    delta = (after - before + 6.0) % 12.0 - 6.0
    assert 0.0 <= delta < 0.01
