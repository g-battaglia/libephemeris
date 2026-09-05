# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Regression tests for issue #62: sealed mode refuses on the date alone.

In sealed ``leb`` mode the interpolated apsides answered for every date
tested, down to year -13000 and up to +17000, while every other body raised
``EphemerisRangeError`` outside the loaded artefact's coverage.

The cause is not that these points skip a check but that they never reach
one. ``MEAN_NODE``, ``MEAN_APOG``, ``INTP_APOG`` and ``INTP_PERG`` are all
served from packaged runtime models rather than from a stored channel, so no
reader lookup fails for them and the sealed range policy — which is applied
where a reader exception is handled — was never invoked.

What made it look like only the apsides were affected is that the policy did
fire whenever the requested frame happened to need nutation, and the four
points need it under complementary conditions:

    MEAN_NODE / MEAN_APOG   refused of-date, answered under FLG_NONUT
    INTP_APOG / INTP_PERG   answered of-date, refused under FLG_NONUT

So all four leaked; the flag decided which. The refusal now depends on the
date alone.

The interpolated apsides carry a second bound. They live in the ``apogee``
group, so a configuration that installs only ``<tier>_core`` — the bundled
default, and the configuration the issue was reported in — declares no
coverage for them at all, and a check that only consults declared coverage
has nothing left to enforce. Their model was fitted over a bounded interval,
and that window is a property of the packaged residual grids, not of a file:
it is enforced whatever groups are installed and intersects with any declared
coverage rather than replacing it.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import libephemeris as eph
from libephemeris.exceptions import EphemerisRangeError
from libephemeris.inventory import get_body_coverage
from libephemeris.lunar import _interpolated_apse_model_window

# Well outside every canonical tier: year -3000, JD 0, and year 20000.
_OUT_OF_RANGE_JD = [625673.5, 0.0, 9000000.5]
_IN_RANGE_JD = 2451545.0

# The reviewed core shipped in the wheel, and nothing else: no apogee group.
_BUNDLED_CORE = Path(eph.__file__).parent / "data" / "leb2" / "base_core.leb2"
# Where `libephemeris download` installs the optional tiers and groups.
_HOME_LEB_DIR = Path("~/.libephemeris/leb").expanduser()

_INTERPOLATED_APSIDES = [
    ("INTP_APOG", eph.INTP_APOG),
    ("INTP_PERG", eph.INTP_PERG),
]

_MODEL_SERVED_POINTS = [
    ("MEAN_NODE", eph.MEAN_NODE),
    ("MEAN_APOG", eph.MEAN_APOG),
    ("INTP_APOG", eph.INTP_APOG),
    ("INTP_PERG", eph.INTP_PERG),
]

# The frames that decided the outcome before the guard became deterministic.
_FRAME_FLAGS = [
    ("of-date", 0),
    ("nonut", eph.FLG_NONUT),
    ("j2000", eph.FLG_J2000),
    ("speed", eph.FLG_SPEED),
    ("nonut-speed", eph.FLG_NONUT | eph.FLG_SPEED),
]


@pytest.fixture
def sealed_leb():
    """Run the body of the test in sealed leb mode, then restore the mode."""
    previous = eph.get_calc_mode()
    eph.set_calc_mode("leb")
    try:
        yield
    finally:
        eph.set_calc_mode(previous)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("name", "body"),
    _MODEL_SERVED_POINTS,
    ids=lambda v: v if isinstance(v, str) else "",
)
@pytest.mark.parametrize(
    ("frame", "flags"), _FRAME_FLAGS, ids=lambda v: v if isinstance(v, str) else ""
)
@pytest.mark.parametrize("jd", _OUT_OF_RANGE_JD)
def test_model_served_points_refuse_outside_coverage(
    sealed_leb, name, body, frame, flags, jd
):
    """Sealed mode promises a refusal; the frame must not decide it."""
    with pytest.raises(EphemerisRangeError):
        eph.calc_ut(jd, body, flags)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("name", "body"),
    _MODEL_SERVED_POINTS,
    ids=lambda v: v if isinstance(v, str) else "",
)
@pytest.mark.parametrize(
    ("frame", "flags"), _FRAME_FLAGS, ids=lambda v: v if isinstance(v, str) else ""
)
def test_model_served_points_still_answer_inside_coverage(
    sealed_leb, name, body, frame, flags
):
    """Narrowing the range must not disturb an in-coverage request."""
    position, _retflag = eph.calc_ut(_IN_RANGE_JD, body, flags)
    assert 0.0 <= position[0] < 360.0


@pytest.mark.unit
@pytest.mark.parametrize(
    ("name", "body"),
    _MODEL_SERVED_POINTS,
    ids=lambda v: v if isinstance(v, str) else "",
)
def test_the_tt_entry_point_agrees_with_the_ut_one(sealed_leb, name, body):
    """calc() carried the same gap as calc_ut() and gets the same guard."""
    with pytest.raises(EphemerisRangeError):
        eph.calc(_OUT_OF_RANGE_JD[0], body, 0)
    position, _retflag = eph.calc(_IN_RANGE_JD, body, 0)
    assert 0.0 <= position[0] < 360.0


@pytest.mark.unit
@pytest.mark.parametrize(
    ("name", "body"),
    _MODEL_SERVED_POINTS,
    ids=lambda v: v if isinstance(v, str) else "",
)
def test_the_error_describes_the_coverage_it_enforced(sealed_leb, name, body):
    """The typed error carries the bounds, as it does for every other body."""
    with pytest.raises(EphemerisRangeError) as excinfo:
        eph.calc_ut(_OUT_OF_RANGE_JD[0], body, 0)
    error = excinfo.value
    assert error.body_id == body
    assert error.start_jd is not None and error.end_jd is not None
    assert error.start_jd < error.end_jd


@pytest.mark.unit
def test_the_apsides_agree_with_the_other_bodies(sealed_leb):
    """The reported asymmetry, stated directly: same date, same outcome."""
    for body in (eph.SUN, eph.MOON, eph.MEAN_NODE, eph.INTP_APOG, eph.INTP_PERG):
        with pytest.raises(EphemerisRangeError):
            eph.calc_ut(_OUT_OF_RANGE_JD[0], body)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("name", "body"),
    _MODEL_SERVED_POINTS,
    ids=lambda v: v if isinstance(v, str) else "",
)
@pytest.mark.parametrize(
    ("frame", "flags"), _FRAME_FLAGS, ids=lambda v: v if isinstance(v, str) else ""
)
def test_the_context_api_enforces_the_same_range(sealed_leb, name, body, frame, flags):
    """EphemerisContext has its own LEB fast path and the same gap."""
    from libephemeris.context import EphemerisContext

    context = EphemerisContext()
    with pytest.raises(EphemerisRangeError):
        context.calc_ut(_OUT_OF_RANGE_JD[0], body, flags)
    position, _retflag = context.calc_ut(_IN_RANGE_JD, body, flags)
    assert 0.0 <= position[0] < 360.0


@pytest.mark.unit
@pytest.mark.parametrize("mode", ["auto", "skyfield"])
@pytest.mark.parametrize(
    ("name", "body"),
    [("INTP_APOG", eph.INTP_APOG), ("INTP_PERG", eph.INTP_PERG)],
    ids=lambda v: v if isinstance(v, str) else "",
)
def test_unsealed_modes_keep_the_documented_taper(mode, name, body):
    """Only sealed mode promises a refusal.

    Outside the fitted interval the model's residual table tapers to zero and
    the remaining trig-only longitude is documented with its own error bound.
    That degradation is deliberate and stays available in auto and skyfield.
    """
    previous = eph.get_calc_mode()
    eph.set_calc_mode(mode)
    try:
        position, _retflag = eph.calc_ut(_OUT_OF_RANGE_JD[0], body)
    finally:
        eph.set_calc_mode(previous)
    assert 0.0 <= position[0] < 360.0


# ---------------------------------------------------------------------------
# The model window: enforced with or without an apogee group
# ---------------------------------------------------------------------------


def _entry_point(entry: str):
    """Return the requested calculation entry point as ``f(jd, body, flags)``."""
    if entry == "calc_ut":
        return eph.calc_ut
    if entry == "calc":
        return eph.calc
    from libephemeris.context import EphemerisContext

    context = EphemerisContext()
    if entry == "context.calc_ut":
        return context.calc_ut
    if entry == "context.calc":
        return context.calc
    raise ValueError(entry)


_ENTRY_POINTS = ["calc_ut", "calc", "context.calc_ut", "context.calc"]


@pytest.fixture
def bundled_core_only():
    """Sealed leb mode on the bundled base core alone: the issue's setup.

    The bundled directory holds no ``base_apogee.leb2``, so the interpolated
    apsides have no declared coverage here — the inventory answers ``None``
    for them while the Sun and the mean node report the base interval.
    """
    previous_mode = eph.get_calc_mode()
    previous_leb = eph.state._LEB_FILE
    eph.set_leb_file(str(_BUNDLED_CORE))
    eph.set_calc_mode("leb")
    assert get_body_coverage(eph.INTP_APOG, _IN_RANGE_JD) is None
    assert get_body_coverage(eph.INTP_PERG, _IN_RANGE_JD) is None
    assert get_body_coverage(eph.SUN, _IN_RANGE_JD) is not None
    try:
        yield
    finally:
        eph.set_leb_file(previous_leb)
        eph.set_calc_mode(previous_mode)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("name", "body"),
    _INTERPOLATED_APSIDES,
    ids=lambda v: v if isinstance(v, str) else "",
)
def test_the_window_is_the_residual_grid_itself(name, body):
    """The bounds are read from the generated model, not written down twice."""
    from libephemeris import lunar_apse_model as model

    start, end = _interpolated_apse_model_window(body)
    if body == eph.INTP_APOG:
        grid = (
            model.APOGEE_CT_JD_START,
            model.APOGEE_CT_STEP_DAYS,
            model.APOGEE_CT_COUNT,
        )
    else:
        grid = (
            model.PERIGEE_CT_JD_START,
            model.PERIGEE_CT_STEP_DAYS,
            model.PERIGEE_CT_COUNT,
        )
    grid_start, step, count = grid
    assert start == grid_start
    assert end == grid_start + step * (count - 1)
    # The model documents ~1549--2651 CE; the last sample of a 2- or 10-day
    # grid falls in the closing days of 2650.
    assert eph.revjul(start)[:3] == (1549, 1, 1)
    assert eph.revjul(end)[:2] == (2650, 12)
    assert _interpolated_apse_model_window(eph.MEAN_NODE) is None
    assert _interpolated_apse_model_window(eph.SUN) is None


@pytest.mark.unit
@pytest.mark.parametrize(
    ("name", "body"),
    _INTERPOLATED_APSIDES,
    ids=lambda v: v if isinstance(v, str) else "",
)
@pytest.mark.parametrize("entry", _ENTRY_POINTS)
# Year -3000 (the value reported), JD 0, year 17000 (inside the extended
# tier's stored interval, outside the model), and year 20000.
@pytest.mark.parametrize("jd", [625673.5, 0.0, 7930183.0, 9000000.5])
def test_the_issue_configuration_refuses_without_an_apogee_group(
    bundled_core_only, name, body, entry, jd
):
    """Only base_core.leb2 installed: the apsides refuse like the Sun does.

    Before, ``calc_ut(625673.5, INTP_APOG)`` answered ~226.81 here while the
    Sun and the mean node raised. No group declares coverage for the apsides
    in this configuration, so the refusal can only come from the model
    window — and the typed error must say so.
    """
    calc = _entry_point(entry)
    with pytest.raises(EphemerisRangeError) as excinfo:
        calc(jd, body, 0)
    error = excinfo.value
    assert error.body_id == body
    assert error.requested_jd == jd
    assert (error.start_jd, error.end_jd) == _interpolated_apse_model_window(body)
    assert error.ephemeris_file is None
    # The asymmetry the issue reported, stated directly.
    for other in (eph.SUN, eph.MEAN_NODE):
        with pytest.raises(EphemerisRangeError):
            calc(jd, other, 0)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("name", "body"),
    _INTERPOLATED_APSIDES,
    ids=lambda v: v if isinstance(v, str) else "",
)
@pytest.mark.parametrize("entry", _ENTRY_POINTS)
@pytest.mark.parametrize(
    ("frame", "flags"), _FRAME_FLAGS, ids=lambda v: v if isinstance(v, str) else ""
)
def test_the_issue_configuration_still_answers_inside_the_window(
    bundled_core_only, name, body, entry, frame, flags
):
    """Inside the window the points keep answering, with or without a group.

    J2000 lies in the bundled core's interval; 1585 and 2625 lie outside it
    but inside the model window, and with no declared coverage to narrow the
    window the model still serves them, exactly as before.
    """
    calc = _entry_point(entry)
    for jd in (_IN_RANGE_JD, 2300000.5, 2680000.5):
        position, _retflag = calc(jd, body, flags)
        assert 0.0 <= position[0] < 360.0


@pytest.mark.unit
@pytest.mark.parametrize(
    ("name", "body"),
    _INTERPOLATED_APSIDES,
    ids=lambda v: v if isinstance(v, str) else "",
)
@pytest.mark.parametrize("entry", _ENTRY_POINTS)
def test_the_window_edges_are_inclusive_from_every_entry_point(
    bundled_core_only, name, body, entry
):
    """The exact bounds answer; one step beyond refuses — UT and TT alike.

    Each entry point compares the day number it was given, so the UT and TT
    entries place the edge delta-T apart (minutes against a grid sampled every
    few days). What must hold for both is the shape: inclusive at the bound,
    closed one step past it.
    """
    calc = _entry_point(entry)
    start, end = _interpolated_apse_model_window(body)
    step = 1e-6
    for edge in (start, end):
        position, _retflag = calc(edge, body, 0)
        assert 0.0 <= position[0] < 360.0
    for beyond in (start - step, end + step, start - 1.0, end + 1.0):
        with pytest.raises(EphemerisRangeError):
            calc(beyond, body, 0)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("name", "body"),
    _INTERPOLATED_APSIDES,
    ids=lambda v: v if isinstance(v, str) else "",
)
@pytest.mark.parametrize(
    ("frame", "flags"),
    [("of-date", 0), ("nonut", eph.FLG_NONUT)],
    ids=lambda v: v if isinstance(v, str) else "",
)
def test_the_guard_changes_nothing_inside_the_window(
    bundled_core_only, monkeypatch, name, body, frame, flags
):
    """Bit-identical results with the window guard present and absent.

    A grid of dates across 1850--2150 and 1560--2640, through the UT and the
    TT entry points: the guard only ever raises, so inside the window it must
    leave every returned float untouched.
    """
    from libephemeris import planets

    dates = [eph.julday(1850, 1, 1, 0.0) + i * (300 * 365.25 / 40) for i in range(41)]
    dates += [eph.julday(1560, 1, 1, 0.0) + i * (1080 * 365.25 / 12) for i in range(13)]

    def sample():
        return [
            (entry(jd, body, flags)[0], entry(jd, body, flags)[1])
            for entry in (eph.calc_ut, eph.calc)
            for jd in dates
        ]

    with_guard = sample()
    monkeypatch.setattr(
        planets, "_raise_leb_model_window_miss", lambda body_id, jd: None
    )
    without_guard = sample()
    assert with_guard == without_guard


@pytest.mark.unit
@pytest.mark.parametrize("mode", ["auto", "skyfield"])
@pytest.mark.parametrize(
    ("name", "body"),
    _INTERPOLATED_APSIDES,
    ids=lambda v: v if isinstance(v, str) else "",
)
def test_unsealed_modes_never_consult_the_window(monkeypatch, mode, name, body):
    """The window guard belongs to the sealed policy and is not asked elsewhere."""
    from libephemeris import planets

    def never(body_id, jd):
        raise AssertionError(f"window guard consulted in {mode} mode")

    monkeypatch.setattr(planets, "_raise_leb_model_window_miss", never)
    previous = eph.get_calc_mode()
    eph.set_calc_mode(mode)
    try:
        for jd in (_IN_RANGE_JD, _OUT_OF_RANGE_JD[0]):
            position, _retflag = eph.calc_ut(jd, body)
            assert 0.0 <= position[0] < 360.0
    finally:
        eph.set_calc_mode(previous)


@pytest.mark.parametrize("tier", ["base", "medium", "extended"])
@pytest.mark.parametrize(
    ("name", "body"),
    _INTERPOLATED_APSIDES,
    ids=lambda v: v if isinstance(v, str) else "",
)
def test_an_installed_apogee_group_intersects_with_the_window(tier, name, body):
    """With ``<tier>_apogee.leb2`` attached, both bounds apply.

    Outside the model window the refusal comes from the window, even when the
    group's stored interval is wider (the extended tier stores the channel
    down to -13200 and up to +17191 while the model was fitted on 1549--2650).
    Inside the window a group narrower than the window keeps refusing with
    its own file bounds, exactly as before.

    Needs the optional download in ``~/.libephemeris/leb``; skipped otherwise.
    """
    core = _HOME_LEB_DIR / f"{tier}_core.leb2"
    group = _HOME_LEB_DIR / f"{tier}_apogee.leb2"
    if not (core.is_file() and group.is_file()):
        pytest.skip(f"{core.name} + {group.name} not installed under {_HOME_LEB_DIR}")
    previous_mode = eph.get_calc_mode()
    previous_leb = eph.state._LEB_FILE
    eph.set_leb_file(str(core))
    eph.set_calc_mode("leb")
    try:
        coverage = get_body_coverage(body, _IN_RANGE_JD)
        if coverage is None or coverage.group != "apogee":
            pytest.skip(f"{group.name} is not attached to {core.name}")
        start, end = _interpolated_apse_model_window(body)

        # Outside the window: refused, whatever the file declares.
        for jd in (625673.5, 0.0, 7930183.0):
            with pytest.raises(EphemerisRangeError) as excinfo:
                eph.calc_ut(jd, body)
            assert (excinfo.value.start_jd, excinfo.value.end_jd) == (start, end)

        # Inside the window: the file's own interval still applies.
        position, _retflag = eph.calc_ut(_IN_RANGE_JD, body)
        assert 0.0 <= position[0] < 360.0
        if coverage.jd_start > start:
            inside_window_outside_file = (start + coverage.jd_start) / 2.0
            with pytest.raises(EphemerisRangeError) as excinfo:
                eph.calc_ut(inside_window_outside_file, body)
            assert excinfo.value.start_jd == coverage.jd_start
            assert excinfo.value.ephemeris_file == coverage.data_file
    finally:
        eph.set_leb_file(previous_leb)
        eph.set_calc_mode(previous_mode)
