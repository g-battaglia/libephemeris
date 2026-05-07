"""Tests for fixed-star catalog and batch APIs."""

from __future__ import annotations

import pytest

import libephemeris as swe


STARS = ("Regulus", "Spica", "Aldebaran", "Sirius", "Vega", "Polaris")
JD = 2460000.5


def _assert_batch_matches_loop(stars: tuple[str, ...], flags: int) -> None:
    batch = swe.swe_batch_fixstars_ut(stars, JD, flags)
    assert len(batch) == len(stars)

    for star, batch_result in zip(stars, batch):
        assert batch_result is not None
        loop_result = swe.swe_fixstar_ut(star, JD, flags)
        assert batch_result[1] == loop_result[1]
        assert batch_result[2] == loop_result[2]
        assert batch_result[0] == pytest.approx(loop_result[0], abs=1e-9)


def test_list_fixed_stars_returns_catalog_entries() -> None:
    catalog = swe.list_fixed_stars()
    names = {star.name for star in catalog}

    assert "Regulus" in names
    assert "Spica" in names
    assert all(star.name for star in catalog)
    assert all(star.nomenclature for star in catalog)
    assert all(isinstance(star.magnitude, float) for star in catalog)


@pytest.mark.parametrize(
    "flags",
    [
        0,
        swe.SEFLG_SPEED,
        swe.SEFLG_EQUATORIAL,
        swe.SEFLG_SPEED | swe.SEFLG_EQUATORIAL,
        swe.SEFLG_J2000,
        swe.SEFLG_RADIANS,
    ],
)
def test_batch_fixstars_ut_matches_individual_calls(flags: int) -> None:
    _assert_batch_matches_loop(STARS, flags)


def test_batch_fixstars_ut_matches_individual_sidereal_calls() -> None:
    swe.swe_set_sid_mode(swe.SE_SIDM_LAHIRI)
    try:
        _assert_batch_matches_loop(STARS, swe.SEFLG_SIDEREAL | swe.SEFLG_SPEED)
    finally:
        swe.reset_session()


def test_batch_fixstars_ut_preserves_input_slots_when_skipping_errors() -> None:
    batch = swe.swe_batch_fixstars_ut(
        ("Regulus", "Not A Real Star", "Spica"), JD, 0, skip_errors=True
    )

    assert batch[0] is not None
    assert batch[1] is None
    assert batch[2] is not None


def test_batch_fixstars_ut_raises_by_default_for_unknown_star() -> None:
    with pytest.raises(swe.Error):
        swe.swe_batch_fixstars_ut(("Not A Real Star",), JD, 0)


def test_batch_fixstars_ut_public_alias() -> None:
    assert swe.batch_fixstars_ut is swe.swe_batch_fixstars_ut
