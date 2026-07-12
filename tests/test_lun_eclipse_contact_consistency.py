"""Regression: lunar-eclipse contact functions agree with the canonical phase
times and satisfy their defining geometric conditions.

Background (audit round v10): the standalone contact helpers
calc_lunar_eclipse_{penumbral_first_contact_p1, penumbral_fourth_contact_p4,
umbral_first_contact_u1, umbral_second_contact_u2, umbral_third_contact_u3,
umbral_fourth_contact_u4} used to root a separate shadow model (1/85
atmospheric enlargement + small-angle separation) that disagreed with the
phase times tret[2..7] returned by lun_eclipse_when() (the canonical
_lun_how_core model) by tens to hundreds of seconds. They now delegate to the
same model, so:

  (1) each standalone contact == the matching lun_eclipse_when() tret entry
      (exact), and
  (2) each contact satisfies its DEFINING condition expressed through the
      magnitude functions (which share the same shadow core):
        P1, P4  -> penumbral magnitude == 0   (Moon limb tangent to penumbra)
        U1, U4  -> umbral magnitude    == 0   (Moon limb tangent to umbra)
        U2, U3  -> umbral magnitude    == 1   (Moon fully immersed / emerging)
  (3) the contacts are correctly ordered and roughly symmetric about maximum.

These are oracle-free internal-consistency + defining-condition guards (no
external reference required). The independent-truth validation (NASA's Five
Millennium Canon contact times, the erfa/JPL shadow-geometry oracle) lives in
the local validation sweeps.
"""

from __future__ import annotations

import pytest

from libephemeris import (
    FLG_SWIEPH,
    calc_lunar_eclipse_penumbral_first_contact_p1,
    calc_lunar_eclipse_penumbral_fourth_contact_p4,
    calc_lunar_eclipse_umbral_first_contact_u1,
    calc_lunar_eclipse_umbral_fourth_contact_u4,
    calc_lunar_eclipse_umbral_second_contact_u2,
    calc_lunar_eclipse_umbral_third_contact_u3,
    julday,
    lun_eclipse_penumbral_magnitude,
    lun_eclipse_umbral_magnitude,
    lun_eclipse_when,
)

pytestmark = pytest.mark.slow

# A spread of eclipses: total, partial and penumbral-only.
_ECLIPSES = [
    (2018, 7, "total"),
    (2019, 1, "total"),
    (2022, 5, "total"),
    (2022, 11, "total"),
    (2021, 11, "partial"),
    (2023, 10, "partial"),
    (2020, 1, "penumbral"),
    (2024, 3, "penumbral"),
]


def _eclipse(year, month):
    """(jd_max, tret) for the eclipse on or after the 1st of the month."""
    jd0 = julday(year, month, 1, 0.0)
    _rf, tret = lun_eclipse_when(jd0, FLG_SWIEPH, 0, False)
    return tret[0], tret


@pytest.mark.parametrize("year,month,kind", _ECLIPSES)
def test_contacts_match_when_phase_times_exactly(year, month, kind):
    """Each standalone contact equals the matching lun_eclipse_when tret entry."""
    jd_max, tret = _eclipse(year, month)
    # (standalone fn, tret index)
    pairs = [
        (calc_lunar_eclipse_penumbral_first_contact_p1, 6),
        (calc_lunar_eclipse_penumbral_fourth_contact_p4, 7),
        (calc_lunar_eclipse_umbral_first_contact_u1, 2),
        (calc_lunar_eclipse_umbral_second_contact_u2, 4),
        (calc_lunar_eclipse_umbral_third_contact_u3, 5),
        (calc_lunar_eclipse_umbral_fourth_contact_u4, 3),
    ]
    for fn, idx in pairs:
        got = fn(jd_max, FLG_SWIEPH)
        assert got == pytest.approx(tret[idx], abs=1e-9), (
            f"{year}-{month} {fn.__name__}: {got} != tret[{idx}] {tret[idx]}"
        )


@pytest.mark.parametrize("year,month,kind", _ECLIPSES)
def test_contacts_satisfy_defining_magnitude_conditions(year, month, kind):
    """At each contact the relevant magnitude hits its defining value."""
    jd_max, _tret = _eclipse(year, month)

    p1 = calc_lunar_eclipse_penumbral_first_contact_p1(jd_max, FLG_SWIEPH)
    p4 = calc_lunar_eclipse_penumbral_fourth_contact_p4(jd_max, FLG_SWIEPH)
    # P1/P4 exist for every lunar eclipse: penumbral magnitude tangent (== 0).
    assert p1 > 0.0 and p4 > 0.0
    assert lun_eclipse_penumbral_magnitude(p1, FLG_SWIEPH) == pytest.approx(
        0.0, abs=2e-3
    )
    assert lun_eclipse_penumbral_magnitude(p4, FLG_SWIEPH) == pytest.approx(
        0.0, abs=2e-3
    )

    if kind in ("total", "partial"):
        u1 = calc_lunar_eclipse_umbral_first_contact_u1(jd_max, FLG_SWIEPH)
        u4 = calc_lunar_eclipse_umbral_fourth_contact_u4(jd_max, FLG_SWIEPH)
        assert u1 > 0.0 and u4 > 0.0
        # U1/U4: umbral magnitude tangent (== 0).
        assert lun_eclipse_umbral_magnitude(u1, FLG_SWIEPH) == pytest.approx(
            0.0, abs=2e-3
        )
        assert lun_eclipse_umbral_magnitude(u4, FLG_SWIEPH) == pytest.approx(
            0.0, abs=2e-3
        )

    if kind == "total":
        u2 = calc_lunar_eclipse_umbral_second_contact_u2(jd_max, FLG_SWIEPH)
        u3 = calc_lunar_eclipse_umbral_third_contact_u3(jd_max, FLG_SWIEPH)
        assert u2 > 0.0 and u3 > 0.0
        # U2/U3: Moon fully immersed -> umbral magnitude == 1.
        assert lun_eclipse_umbral_magnitude(u2, FLG_SWIEPH) == pytest.approx(
            1.0, abs=2e-3
        )
        assert lun_eclipse_umbral_magnitude(u3, FLG_SWIEPH) == pytest.approx(
            1.0, abs=2e-3
        )


@pytest.mark.parametrize("year,month,kind", _ECLIPSES)
def test_contacts_ordered_and_symmetric(year, month, kind):
    """Contacts are time-ordered and roughly symmetric about maximum."""
    jd_max, _tret = _eclipse(year, month)
    p1 = calc_lunar_eclipse_penumbral_first_contact_p1(jd_max, FLG_SWIEPH)
    p4 = calc_lunar_eclipse_penumbral_fourth_contact_p4(jd_max, FLG_SWIEPH)

    assert p1 < jd_max < p4
    # Penumbral phase symmetric about maximum to within ~5 min.
    assert abs((jd_max - p1) - (p4 - jd_max)) * 86400.0 < 300.0

    if kind in ("total", "partial"):
        u1 = calc_lunar_eclipse_umbral_first_contact_u1(jd_max, FLG_SWIEPH)
        u4 = calc_lunar_eclipse_umbral_fourth_contact_u4(jd_max, FLG_SWIEPH)
        assert p1 < u1 < jd_max < u4 < p4
        assert abs((jd_max - u1) - (u4 - jd_max)) * 86400.0 < 300.0

    if kind == "total":
        u1 = calc_lunar_eclipse_umbral_first_contact_u1(jd_max, FLG_SWIEPH)
        u2 = calc_lunar_eclipse_umbral_second_contact_u2(jd_max, FLG_SWIEPH)
        u3 = calc_lunar_eclipse_umbral_third_contact_u3(jd_max, FLG_SWIEPH)
        u4 = calc_lunar_eclipse_umbral_fourth_contact_u4(jd_max, FLG_SWIEPH)
        assert u1 < u2 < jd_max < u3 < u4
        assert abs((jd_max - u2) - (u3 - jd_max)) * 86400.0 < 300.0
