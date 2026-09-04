# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Regression tests for issue #65: two helpers that answered instead of failing.

``calc_airmass()`` accepted any ``method`` and silently computed Kasten-Young
for it. An empty string, ``None``, ``0``, a misspelling and even a
differently-cased spelling of a real method all returned the same number, from
a model the caller did not ask for and with nothing in the result to say so.

``calc_eye_adaptation_state()`` and its documentation named different sets of
values. The prose presented "scotopic" as a returned state; the function
returns ``"dark"``, which is the token the whole module speaks. Code comparing
against the documented name never matched.

The same silent fallback existed one level down: an unrecognised
``eye_adaptation`` passed to ``calc_contrast_threshold`` landed on the
photopic branch.
"""

from __future__ import annotations

import pytest

import libephemeris as eph
from libephemeris import extinction as ext
from libephemeris.exceptions import Error, InputValidationError

# ---------------------------------------------------------------------------
# calc_airmass
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize("method", ext.AIRMASS_METHODS)
def test_every_implemented_method_is_accepted(method):
    """Narrowing the input must not cost a model the module implements."""
    airmass = eph.calc_airmass(30.0, method)
    assert 1.0 <= airmass <= 40.0


@pytest.mark.unit
def test_the_three_methods_are_actually_different_models():
    """If they all returned the same number the guard would be pointless."""
    values = {m: eph.calc_airmass(30.0, m) for m in ext.AIRMASS_METHODS}
    assert len(set(values.values())) == len(ext.AIRMASS_METHODS)


@pytest.mark.unit
@pytest.mark.parametrize(
    "method",
    [
        "",
        "nope",
        "kasten-young",  # right model, wrong separator
        "KASTEN_YOUNG",  # right model, wrong case
        "Rozenberg",
        None,
        0,
        1,
        [],
        {},
        set(),
    ],
)
def test_a_method_the_module_does_not_implement_is_refused(method):
    """The reported case: a misspelling used to return Kasten-Young."""
    with pytest.raises(InputValidationError):
        eph.calc_airmass(30.0, method)


@pytest.mark.unit
def test_the_message_names_the_methods_that_do_exist():
    with pytest.raises(InputValidationError) as excinfo:
        eph.calc_airmass(30.0, "kasten-young")
    message = str(excinfo.value)
    for method in ext.AIRMASS_METHODS:
        assert method in message


@pytest.mark.unit
def test_the_method_is_validated_before_the_horizon_shortcut():
    """An object below the horizon returns 40.0 — but not for a bad method."""
    assert eph.calc_airmass(-5.0, "kasten_young") == 40.0
    with pytest.raises(InputValidationError):
        eph.calc_airmass(-5.0, "nope")


@pytest.mark.unit
def test_the_default_is_unchanged():
    assert eph.calc_airmass(30.0) == eph.calc_airmass(30.0, "kasten_young")


# ---------------------------------------------------------------------------
# calc_eye_adaptation_state
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize(
    ("sky_brightness", "expected"),
    [
        (-5.0, "photopic"),
        (3.0, "photopic"),
        (15.999, "photopic"),
        (16.0, "mesopic"),
        (18.0, "mesopic"),
        (19.999, "mesopic"),
        (20.0, "dark"),
        (21.5, "dark"),
        (1e5, "dark"),
    ],
)
def test_all_three_states_are_reachable(sky_brightness, expected):
    """The issue sampled no value in [16, 20); the mesopic band is real."""
    assert eph.calc_eye_adaptation_state(sky_brightness) == expected


@pytest.mark.unit
def test_the_returned_set_is_exactly_what_the_module_declares():
    """Documentation and code have to name the same set."""
    produced = {
        eph.calc_eye_adaptation_state(v)
        for v in [-10, 0, 10, 15.9, 16, 17, 19.9, 20, 25, 1e6]
    }
    assert produced == set(ext.EYE_ADAPTATION_STATES)


@pytest.mark.unit
def test_scotopic_is_not_a_returned_value():
    """It is the regime's name; "dark" is the token, and the docs say so."""
    produced = {eph.calc_eye_adaptation_state(v) for v in range(-10, 40)}
    assert "scotopic" not in produced
    assert "dark" in produced


@pytest.mark.unit
def test_scotopic_is_accepted_where_a_state_is_taken_as_input():
    """Code written against the old prose keeps working."""
    assert ext.calc_contrast_threshold(21.5, "scotopic") == ext.calc_contrast_threshold(
        21.5, "dark"
    )


@pytest.mark.unit
def test_the_reported_state_is_the_canonical_token():
    result = ext.calc_visibility_threshold(3.0, 21.5, eye_adaptation="scotopic")
    assert result.eye_adaptation == "dark"


# ---------------------------------------------------------------------------
# The same silent fallback one level down
# ---------------------------------------------------------------------------


@pytest.mark.unit
# None is the documented auto-determine sentinel, not an unknown state.
@pytest.mark.parametrize("state", ["", "nope", "Dark", "DARK", 0, 2, [], {}, set()])
def test_an_unknown_eye_adaptation_is_refused(state):
    """It used to fall through to the photopic branch."""
    with pytest.raises(InputValidationError):
        ext.calc_contrast_threshold(21.5, state)


@pytest.mark.unit
@pytest.mark.parametrize("state", list(ext.EYE_ADAPTATION_STATES) + ["scotopic"])
def test_every_documented_state_is_accepted(state):
    assert ext.calc_contrast_threshold(21.5, state) > 0.0


@pytest.mark.unit
def test_the_three_states_give_three_different_thresholds():
    """Refusing an unknown state matters because the branches differ."""
    thresholds = {
        state: ext.calc_contrast_threshold(18.0, state)
        for state in ext.EYE_ADAPTATION_STATES
    }
    assert len(set(thresholds.values())) == 3


@pytest.mark.unit
def test_none_still_derives_the_state_from_the_sky():
    assert ext.calc_contrast_threshold(21.5, None) == ext.calc_contrast_threshold(
        21.5, "dark"
    )


@pytest.mark.unit
def test_the_limiting_magnitude_helper_takes_the_alias_and_refuses_the_rest():
    """It delegates to calc_visibility_threshold and inherits its contract."""
    assert ext.calc_limiting_magnitude_for_sky(
        21.5, eye_adaptation="scotopic"
    ) == ext.calc_limiting_magnitude_for_sky(21.5, eye_adaptation="dark")
    with pytest.raises(InputValidationError):
        ext.calc_limiting_magnitude_for_sky(21.5, eye_adaptation="nope")


# ---------------------------------------------------------------------------
# normalize_eye_adaptation, called directly
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize("state", ext.EYE_ADAPTATION_STATES)
def test_a_canonical_state_normalizes_to_itself(state):
    assert ext.normalize_eye_adaptation(state, "f") == state


@pytest.mark.unit
def test_scotopic_normalizes_to_dark():
    assert ext.normalize_eye_adaptation("scotopic", "f") == "dark"


@pytest.mark.unit
@pytest.mark.parametrize("state", [[], {}, set()])
def test_an_unhashable_state_is_refused_not_crashed_on(state):
    """A list, dict or set must raise the validation error, not TypeError."""
    with pytest.raises(InputValidationError):
        ext.normalize_eye_adaptation(state, "f")


@pytest.mark.unit
@pytest.mark.parametrize("state", [None, 0, 1, 2.0, b"dark", ("dark",), object()])
def test_a_non_string_state_is_refused(state):
    """None is the callers' sentinel; the normalizer itself takes strings only."""
    with pytest.raises(InputValidationError):
        ext.normalize_eye_adaptation(state, "f")


@pytest.mark.unit
@pytest.mark.parametrize(
    "state", ["", "Dark", "SCOTOPIC", "Scotopic", " dark", "dark "]
)
def test_matching_is_exact_and_case_sensitive(state):
    with pytest.raises(InputValidationError):
        ext.normalize_eye_adaptation(state, "f")


@pytest.mark.unit
def test_the_message_names_the_caller_and_the_accepted_states():
    with pytest.raises(InputValidationError) as excinfo:
        ext.normalize_eye_adaptation("nope", "calc_contrast_threshold")
    message = str(excinfo.value)
    assert message.startswith("calc_contrast_threshold:")
    for state in ext.EYE_ADAPTATION_STATES:
        assert state in message
    assert "scotopic" in message


@pytest.mark.unit
def test_the_errors_stay_inside_the_documented_hierarchy():
    assert issubclass(InputValidationError, Error)
    with pytest.raises(Error):
        eph.calc_airmass(30.0, "nope")
    with pytest.raises(Error):
        ext.calc_contrast_threshold(21.5, "nope")
