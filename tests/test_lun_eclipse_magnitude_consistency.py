"""Regression: lunar-eclipse magnitude convenience functions must agree with
the canonical lun_eclipse_how() shadow model and with NASA's published values.

Background (audit round v10): lun_eclipse_umbral_magnitude() and
lun_eclipse_penumbral_magnitude() used to route through a *second* shadow model
(_calculate_lunar_eclipse_type_and_magnitude, with a small-angle Danjon
approximation) that disagreed with the canonical exact-cone path. The canonical
path now implements the independently published Danjon equations used in
NASA/TP-2009-214173, without output-fitted shadow factors. The convenience
functions delegate to that same core, so:

  (1) lun_eclipse_umbral_magnitude(jd)   == lun_eclipse_how(jd, geo)[1][0]   (exact)
  (2) lun_eclipse_penumbral_magnitude(jd)== max(0, lun_eclipse_how(jd,geo)[1][1]) (exact)
  (3) the umbral magnitude reproduces NASA's independently published catalog
      to <0.003.

These are independent-reference checks (NASA values + cross-API identity), not
type/range guards. Reference: https://eclipse.gsfc.nasa.gov/lunar.html
"""

from __future__ import annotations

import pytest

from libephemeris import (
    ECL_PARTIAL,
    ECL_PENUMBRAL,
    ECL_TOTAL,
    FLG_SWIEPH,
    julday,
    lun_eclipse_how,
    lun_eclipse_penumbral_magnitude,
    lun_eclipse_umbral_magnitude,
    lun_eclipse_when,
)
from libephemeris.eclipse import _lun_eclipse_how_pythonic, _lun_how_core

pytestmark = pytest.mark.slow

_GEO = (0.0, 0.0, 0.0)  # lunar-eclipse magnitude is observer-independent

# NASA Five Millennium Canon of Lunar Eclipses: umbral magnitude at greatest
# eclipse. (year, month, umbral_magnitude, kind)
_NASA_UMBRAL = [
    (2018, 7, 1.6087, "total"),  # longest total eclipse of the 21st century
    (2019, 1, 1.1956, "total"),
    (2021, 5, 1.0095, "total"),
    (2022, 5, 1.4141, "total"),
    (2022, 11, 1.3589, "total"),
    (2021, 11, 0.9742, "partial"),  # deep partial
    (2023, 10, 0.1227, "partial"),  # shallow partial
]


def _find_max(year: int, month: int) -> float:
    """JD(UT) of the greatest lunar eclipse on or after the 1st of the month."""
    jd0 = julday(year, month, 1, 0.0)
    _rf, tret = lun_eclipse_when(jd0, FLG_SWIEPH, 0, False)
    return tret[0]


@pytest.mark.parametrize("year,month,_nasa,_kind", _NASA_UMBRAL)
def test_umbral_magnitude_matches_how_exactly(year, month, _nasa, _kind):
    """The convenience umbral magnitude is identical to lun_eclipse_how attr[0]."""
    jd_max = _find_max(year, month)
    conv = lun_eclipse_umbral_magnitude(jd_max, FLG_SWIEPH)
    _rc, attr = lun_eclipse_how(jd_max, _GEO, FLG_SWIEPH)
    assert conv == pytest.approx(attr[0], abs=1e-12), (
        f"{year}-{month}: umbral_magnitude {conv} != lun_eclipse_how[0] {attr[0]}"
    )


@pytest.mark.parametrize("year,month,_nasa,_kind", _NASA_UMBRAL)
def test_penumbral_magnitude_matches_how_exactly(year, month, _nasa, _kind):
    """The convenience penumbral magnitude is identical to lun_eclipse_how attr[1]."""
    jd_max = _find_max(year, month)
    conv = lun_eclipse_penumbral_magnitude(jd_max, FLG_SWIEPH)
    _rc, attr = lun_eclipse_how(jd_max, _GEO, FLG_SWIEPH)
    assert conv == pytest.approx(max(0.0, attr[1]), abs=1e-12), (
        f"{year}-{month}: penumbral_magnitude {conv} != "
        f"max(0, lun_eclipse_how[1] {attr[1]})"
    )


@pytest.mark.parametrize("year,month,nasa,kind", _NASA_UMBRAL)
def test_umbral_magnitude_matches_nasa(year, month, nasa, kind):
    """The umbral magnitude reproduces NASA's catalog value to <0.003."""
    jd_max = _find_max(year, month)
    mag = lun_eclipse_umbral_magnitude(jd_max, FLG_SWIEPH)
    assert mag == pytest.approx(nasa, abs=0.003), (
        f"{year}-{month} ({kind}): umbral magnitude {mag} vs NASA {nasa}"
    )
    # Classification must be on the correct side of the total/partial boundary.
    if kind == "total":
        assert mag >= 1.0
    else:
        assert 0.0 < mag < 1.0


# A spread of total, partial and penumbral-only eclipses for the pythonic check.
_PYTHONIC_ECLIPSES = [
    (2018, 7, "total"),
    (2022, 5, "total"),
    (2021, 11, "partial"),
    (2023, 10, "partial"),
    (2020, 1, "penumbral"),
    (2024, 3, "penumbral"),
]


@pytest.mark.parametrize("year,month,kind", _PYTHONIC_ECLIPSES)
def test_how_pythonic_draws_type_and_magnitude_from_core(year, month, kind):
    """_lun_eclipse_how_pythonic must take type+magnitudes from _lun_how_core.

    Regression for audit v10: the exported _lun_eclipse_how_pythonic() used to
    compute its magnitudes/type from the divergent
    _calculate_lunar_eclipse_type_and_magnitude() model (1/85 enlargement +
    small-angle), disagreeing with the canonical lun_eclipse_how() path by up
    to ~0.011. It now delegates to the same _lun_how_core shadow core. The
    magnitudes and distance-from-opposition are observer-independent, so they
    must equal the core's attr exactly; the eclipse type (after masking the
    VISIBLE flags pythonic adds for an above-horizon observer) must equal the
    core's geocentric classification.
    """
    jd_max = _find_max(year, month)
    retc_core, attr_core, _dc = _lun_how_core(jd_max, FLG_SWIEPH)
    rc_py, attr_py = _lun_eclipse_how_pythonic(jd_max, 41.9, 12.5)

    assert attr_py[0] == pytest.approx(attr_core[0], abs=1e-12)  # umbral mag
    assert attr_py[1] == pytest.approx(attr_core[1], abs=1e-12)  # penumbral mag
    assert attr_py[8] == pytest.approx(attr_core[8], abs=1e-12)  # == attr[0]
    assert attr_py[7] == pytest.approx(attr_core[7], abs=1e-12)  # dist. from opp.

    type_mask = ECL_TOTAL | ECL_PARTIAL | ECL_PENUMBRAL
    assert (rc_py & type_mask) == retc_core, (
        f"{year}-{month} ({kind}): pythonic type {rc_py & type_mask} "
        f"!= core {retc_core}"
    )
