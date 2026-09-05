# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Property-based tests for ``split_deg`` and ``csroundsec``.

The two functions are pure arithmetic on the caller's argument: no date, no
frame, no ephemeris. That makes them a good fit for randomised testing, and
the invariants below are the ones the behavioural specification of the family
states — ranges and types, the direction and division reported in the fifth
element, recomposition of the parts, the rounding rules, how the ``KEEP_*``
flags hold the rounding back, and the monotonicity of both functions.

Angles are drawn over +/-1500 degrees (well past a turn in either direction),
including subnormals and the immediate neighbourhoods of every 30-degree sign
boundary and every 13d20' segment boundary, and each is split against all 128
combinations of the seven ``SPLIT_DEG_*`` bits. Centisecond values are drawn
over +/-2**32.
"""

from __future__ import annotations

import math
from collections.abc import Callable

import pytest
from hypothesis import assume, given, settings
from hypothesis import strategies as st

from libephemeris.utils import (
    SPLIT_DEG_KEEP_DEG,
    SPLIT_DEG_KEEP_SIGN,
    SPLIT_DEG_NAKSHATRA,
    SPLIT_DEG_ROUND_DEG,
    SPLIT_DEG_ROUND_MIN,
    SPLIT_DEG_ROUND_SEC,
    SPLIT_DEG_ZODIACAL,
    csroundsec,
    split_deg,
)

pytestmark = pytest.mark.unit

#: The seven flag bits, and the three that ask for a rounding.
SPLIT_BITS = (
    SPLIT_DEG_ROUND_SEC,
    SPLIT_DEG_ROUND_MIN,
    SPLIT_DEG_ROUND_DEG,
    SPLIT_DEG_ZODIACAL,
    SPLIT_DEG_KEEP_SIGN,
    SPLIT_DEG_KEEP_DEG,
    SPLIT_DEG_NAKSHATRA,
)
ROUNDING_MASK = SPLIT_DEG_ROUND_SEC | SPLIT_DEG_ROUND_MIN | SPLIT_DEG_ROUND_DEG
KEEP_MASK = SPLIT_DEG_KEEP_SIGN | SPLIT_DEG_KEEP_DEG
DIVISION_MASK = SPLIT_DEG_ZODIACAL | SPLIT_DEG_NAKSHATRA

#: Each rounding unit, in degrees, coarsest first.
ROUNDING_UNITS = (
    (SPLIT_DEG_ROUND_DEG, 1.0),
    (SPLIT_DEG_ROUND_MIN, 1.0 / 60.0),
    (SPLIT_DEG_ROUND_SEC, 1.0 / 3600.0),
)

#: One 30-degree sign and one 13d20' segment, in degrees.
SIGN_SPAN = 30.0
SEGMENT_SPAN = 360.0 / 27.0

#: Hundredths of an arcsecond in one arcsecond and in one 30-degree sign.
CS_PER_ARCSEC = 100
CS_PER_SIGN = 30 * 3600 * 100


def _flag_words() -> tuple[int, ...]:
    """All 128 subsets of the seven flag bits, as flag words."""
    words = []
    for mask in range(1 << len(SPLIT_BITS)):
        word = 0
        for index, bit in enumerate(SPLIT_BITS):
            if mask & (1 << index):
                word |= bit
        words.append(word)
    return tuple(sorted(words))


def _neighbourhood(value: float) -> tuple[float, ...]:
    """A value and its two immediate float neighbours."""
    return (
        math.nextafter(value, -math.inf),
        value,
        math.nextafter(value, math.inf),
    )


def _boundary_degrees() -> tuple[float, ...]:
    """Every sign and segment boundary within a few turns, plus its neighbours."""
    values: set[float] = set()
    for turn in range(-4, 5):
        for index in range(12):
            values.update(_neighbourhood(360.0 * turn + SIGN_SPAN * index))
        for index in range(27):
            values.update(_neighbourhood(360.0 * turn + SEGMENT_SPAN * index))
    values.update((1e-13, -1e-13, 5e-324, -5e-324, 1e-310, -1234.5678))
    return tuple(sorted(values)) + (-0.0,)


BOUNDARY_DEGREES = _boundary_degrees()

FLAG_WORDS = _flag_words()


def _degrees(low: float, high: float) -> st.SearchStrategy[float]:
    """Angles over a range: a dense draw plus the boundaries inside it."""
    edges = tuple(d for d in BOUNDARY_DEGREES if low <= d <= high)
    dense = st.floats(
        min_value=low, max_value=high, allow_nan=False, allow_infinity=False
    )
    return st.one_of(dense, st.sampled_from(edges)) if edges else dense


def _words(predicate: "Callable[[int], bool]") -> st.SearchStrategy[int]:
    """The flag words a property applies to, chosen rather than filtered."""
    return st.sampled_from(tuple(w for w in FLAG_WORDS if predicate(w)))


degrees = _degrees(-1500.0, 1500.0)
negative_degrees = _degrees(-1500.0, -5e-324)
non_negative_degrees = _degrees(0.0, 1500.0)
inside_one_turn = _degrees(-359.999999999999, 359.999999999999)
flag_words = _words(lambda word: True)
centiseconds = st.integers(min_value=-(2**32), max_value=2**32)


def _nakshatra_applies(degree: float, roundflag: int) -> bool:
    """Say whether this call takes the lunar-segment path."""
    return bool(roundflag & SPLIT_DEG_NAKSHATRA) and degree >= 0.0


def _coarsest_bit(roundflag: int) -> int:
    """The rounding bit in force: the coarsest one the word carries."""
    for bit, _ in ROUNDING_UNITS:
        if roundflag & bit:
            return bit
    return 0


def _requested_unit(roundflag: int) -> float:
    """The rounding unit in force, in degrees, or 0.0 when there is none."""
    for bit, unit in ROUNDING_UNITS:
        if roundflag & bit:
            return unit
    return 0.0


def _total_arcsec(parts: tuple[int, int, int, float, int]) -> float:
    """The first four fields of a result, added up as arcseconds."""
    ideg, imin, isec, secfr, _ = parts
    return ideg * 3600.0 + imin * 60.0 + isec + secfr


# ---------------------------------------------------------------------------
# Ranges and types
# ---------------------------------------------------------------------------
@settings(max_examples=400, deadline=None)
@given(degree=degrees, roundflag=flag_words)
def test_result_is_five_native_python_objects(degree: float, roundflag: int) -> None:
    """Property 1: four native ints and one native float, in that order."""
    result = split_deg(degree, roundflag)
    assert len(result) == 5
    for position in (0, 1, 2, 4):
        assert type(result[position]) is int
    assert type(result[3]) is float


@settings(max_examples=400, deadline=None)
@given(degree=degrees, roundflag=flag_words)
def test_minutes_and_seconds_are_sexagesimal(degree: float, roundflag: int) -> None:
    """Property 2: arcminutes and arcseconds stay in 0..59."""
    _, imin, isec, _, _ = split_deg(degree, roundflag)
    assert 0 <= imin <= 59
    assert 0 <= isec <= 59


@settings(max_examples=400, deadline=None)
@given(degree=degrees, roundflag=flag_words)
def test_degrees_are_bounded_by_the_division(degree: float, roundflag: int) -> None:
    """Property 3: degrees are never negative, and a division caps them."""
    ideg = split_deg(degree, roundflag)[0]
    assert ideg >= 0
    if _nakshatra_applies(degree, roundflag):
        assert ideg <= 13
    elif roundflag & SPLIT_DEG_ZODIACAL:
        assert ideg <= 29


@settings(max_examples=400, deadline=None)
@given(degree=degrees, roundflag=flag_words)
def test_sub_arcsecond_field_answers_the_rounding_bits(
    degree: float, roundflag: int
) -> None:
    """Property 4: a true fraction without a rounding bit, the arcseconds with one."""
    _, _, isec, secfr, _ = split_deg(degree, roundflag)
    if roundflag & ROUNDING_MASK:
        assert secfr == float(isec)
    else:
        assert 0.0 <= secfr < 1.0


# ---------------------------------------------------------------------------
# Direction and division
# ---------------------------------------------------------------------------
@settings(max_examples=400, deadline=None)
@given(degree=degrees, roundflag=_words(lambda word: not word & DIVISION_MASK))
def test_direction_is_reported_without_a_division_flag(
    degree: float, roundflag: int
) -> None:
    """Property 5: the fifth element is the direction, and nothing else."""
    sign = split_deg(degree, roundflag)[4]
    assert sign == (1 if degree >= 0.0 else -1)


@settings(max_examples=400, deadline=None)
@given(
    degree=degrees,
    roundflag=_words(
        lambda word: word & SPLIT_DEG_ZODIACAL and not word & SPLIT_DEG_NAKSHATRA
    ),
)
def test_zodiacal_reports_an_index_and_forgets_the_direction(
    degree: float, roundflag: int
) -> None:
    """Property 6: a zodiacal index is never negative and ignores the direction."""
    result = split_deg(degree, roundflag)
    assert result[4] >= 0
    assert result == split_deg(-degree, roundflag)


@settings(max_examples=400, deadline=None)
@given(
    degree=negative_degrees,
    roundflag=_words(
        lambda word: word & SPLIT_DEG_NAKSHATRA and not word & SPLIT_DEG_ZODIACAL
    ),
)
def test_negative_nakshatra_falls_back_to_the_signed_split(
    degree: float, roundflag: int
) -> None:
    """Property 7: a negative angle never enters the lunar-segment path."""
    assert split_deg(degree, roundflag) == split_deg(
        degree, roundflag & ~SPLIT_DEG_NAKSHATRA
    )


@settings(max_examples=400, deadline=None)
@given(
    degree=degrees,
    roundflag=_words(
        lambda word: word & SPLIT_DEG_ZODIACAL and word & SPLIT_DEG_NAKSHATRA
    ),
)
def test_both_division_flags_are_settled_by_the_direction(
    degree: float, roundflag: int
) -> None:
    """Property 8: the segments win for a non-negative angle, the signs otherwise."""
    dropped = SPLIT_DEG_ZODIACAL if degree >= 0.0 else SPLIT_DEG_NAKSHATRA
    assert split_deg(degree, roundflag) == split_deg(degree, roundflag & ~dropped)


# ---------------------------------------------------------------------------
# Recomposition
# ---------------------------------------------------------------------------
@settings(max_examples=400, deadline=None)
@given(
    degree=degrees,
    roundflag=_words(lambda word: not word & (ROUNDING_MASK | DIVISION_MASK)),
)
def test_plain_parts_recompose_the_magnitude(degree: float, roundflag: int) -> None:
    """Property 9: without rounding, the parts add back up to the magnitude."""
    ideg, imin, isec, secfr, _ = split_deg(degree, roundflag)
    recomposed = ideg + imin / 60.0 + (isec + secfr) / 3600.0
    assert abs(recomposed - abs(degree)) <= 1e-12


@settings(max_examples=400, deadline=None)
@given(
    degree=inside_one_turn,
    roundflag=_words(
        lambda word: word & SPLIT_DEG_ZODIACAL
        and not word & (SPLIT_DEG_NAKSHATRA | ROUNDING_MASK)
    ),
)
def test_zodiacal_parts_recompose_the_magnitude(degree: float, roundflag: int) -> None:
    """Property 10: the sign index and the position add back up to the magnitude."""
    ideg, imin, isec, secfr, index = split_deg(degree, roundflag)
    recomposed = SIGN_SPAN * index + ideg + imin / 60.0 + (isec + secfr) / 3600.0
    assert abs(recomposed - abs(degree)) <= 1e-12


@settings(max_examples=400, deadline=None)
@given(
    degree=_degrees(0.0, 359.999999999999),
    roundflag=_words(
        lambda word: word & SPLIT_DEG_NAKSHATRA and not word & ROUNDING_MASK
    ),
)
def test_nakshatra_parts_recompose_the_angle(degree: float, roundflag: int) -> None:
    """Property 11: the segment index and the position add back up to the angle.

    The comparison is made in arcseconds, the domain this path reduces in,
    where a segment is the exact integer 48000.
    """
    ideg, imin, isec, secfr, index = split_deg(degree, roundflag)
    recomposed = index * 48000.0 + ideg * 3600.0 + imin * 60.0 + isec + secfr
    assert abs(recomposed - degree * 3600.0) <= 3.6e-9


@settings(max_examples=600, deadline=None)
@given(
    degree=degrees,
    roundflag=_words(
        lambda word: word & ROUNDING_MASK and not word & (DIVISION_MASK | KEEP_MASK)
    ),
)
def test_rounding_matches_half_up_to_the_requested_unit(
    degree: float, roundflag: int
) -> None:
    """Property 12: the parts down to the requested unit are the half-up rounding.

    The identity is checked on the result rather than by recomposing, because
    everything below the requested unit is carry. Its right-hand side is a
    second, different, floating-point evaluation of the same real quantity:
    the two agree everywhere except within one unit in the last place of an
    exact tie, where they can name adjacent multiples of the unit.
    """
    unit = _requested_unit(roundflag)
    ideg, imin, isec, _, _ = split_deg(degree, roundflag)
    if unit == 1.0:
        reported = ideg
    elif unit == 1.0 / 60.0:
        reported = ideg * 60 + imin
    else:
        reported = ideg * 3600 + imin * 60 + isec
    scaled = abs(degree) / unit
    expected = math.floor(scaled + 0.5)
    assert abs(reported - expected) <= 1
    if abs(scaled - math.floor(scaled) - 0.5) > 4.0 * math.ulp(max(scaled, 1.0)):
        assert reported == expected


# ---------------------------------------------------------------------------
# Flags
# ---------------------------------------------------------------------------
@settings(max_examples=400, deadline=None)
@given(
    degree=degrees,
    roundflag=_words(
        lambda word: _coarsest_bit(word) and word & ROUNDING_MASK != _coarsest_bit(word)
    ),
)
def test_coarsest_rounding_bit_wins(degree: float, roundflag: int) -> None:
    """Property 13: the finer rounding bits have no effect at all."""
    coarsest = _coarsest_bit(roundflag)
    alone = (roundflag & ~ROUNDING_MASK) | coarsest
    assert split_deg(degree, roundflag) == split_deg(degree, alone)


@settings(max_examples=400, deadline=None)
@given(
    degree=degrees,
    roundflag=_words(lambda word: word & KEEP_MASK and not word & ROUNDING_MASK),
)
def test_keep_flags_do_nothing_without_a_rounding_bit(
    degree: float, roundflag: int
) -> None:
    """Property 14: with nothing to hold back, both KEEP_* flags are inert."""
    assert split_deg(degree, roundflag) == split_deg(degree, roundflag & ~KEEP_MASK)


@settings(max_examples=400, deadline=None)
@given(
    degree=degrees,
    roundflag=_words(
        lambda word: word & SPLIT_DEG_KEEP_SIGN and not word & SPLIT_DEG_KEEP_DEG
    ),
)
def test_keep_sign_never_moves_the_division(degree: float, roundflag: int) -> None:
    """Property 15: with KEEP_SIGN the fifth element is the unrounded one."""
    unrounded = split_deg(degree, roundflag & ~ROUNDING_MASK)
    assert split_deg(degree, roundflag)[4] == unrounded[4]


@settings(max_examples=400, deadline=None)
@given(degree=degrees, roundflag=flag_words)
def test_keep_deg_dominates_keep_sign(degree: float, roundflag: int) -> None:
    """Property 16: when both are set only the whole-degree test governs."""
    both = roundflag | SPLIT_DEG_KEEP_DEG | SPLIT_DEG_KEEP_SIGN
    alone = both & ~SPLIT_DEG_KEEP_SIGN
    assert split_deg(degree, both) == split_deg(degree, alone)


@settings(max_examples=400, deadline=None)
@given(degree=degrees, roundflag=_words(lambda word: word & SPLIT_DEG_KEEP_DEG))
def test_keep_deg_preserves_the_degrees_field(degree: float, roundflag: int) -> None:
    """Property 17: outside the lunar segments KEEP_DEG holds the degrees field.

    A segment boundary falls in the middle of a degree, so the lunar path can
    carry into the next segment without the whole-degree test ever firing;
    that case is pinned by example instead.
    """
    assume(not _nakshatra_applies(degree, roundflag))
    unrounded = split_deg(degree, roundflag & ~ROUNDING_MASK)
    assert split_deg(degree, roundflag)[0] == unrounded[0]


def test_keep_deg_can_still_change_the_degrees_field_in_a_segment() -> None:
    """The one exception to property 17, with the witness from the family."""
    flags = SPLIT_DEG_NAKSHATRA | SPLIT_DEG_ROUND_DEG | SPLIT_DEG_KEEP_DEG
    assert split_deg(106.344, SPLIT_DEG_NAKSHATRA)[0] == 13
    assert split_deg(106.344, flags) == (0, 10, 38, 38.0, 8)


# ---------------------------------------------------------------------------
# Monotonicity and idempotence
# ---------------------------------------------------------------------------
@settings(max_examples=400, deadline=None)
@given(first=non_negative_degrees, second=non_negative_degrees)
def test_split_is_monotone_in_the_magnitude(first: float, second: float) -> None:
    """Property 18: a larger angle never splits to a smaller total."""
    low, high = sorted((first, second))
    assert _total_arcsec(split_deg(low)) <= _total_arcsec(split_deg(high))


@settings(max_examples=400, deadline=None)
@given(cs=centiseconds)
def test_csroundsec_returns_whole_arcseconds_close_to_the_input(cs: int) -> None:
    """Property 19: whole arcseconds, and never a full arcsecond and a half away."""
    rounded = csroundsec(cs)
    assert type(rounded) is int
    assert rounded % CS_PER_ARCSEC == 0
    assert abs(rounded - cs) < 150


@settings(max_examples=400, deadline=None)
@given(first=centiseconds, second=centiseconds)
def test_csroundsec_is_monotone(first: int, second: int) -> None:
    """Property 20: rounding never reverses the order of two angles."""
    low, high = sorted((first, second))
    assert csroundsec(low) <= csroundsec(high)


@settings(max_examples=400, deadline=None)
@given(cs=st.integers(min_value=0, max_value=2**32))
def test_csroundsec_is_idempotent_on_whole_arcseconds(cs: int) -> None:
    """Property 21: a non-negative whole arcsecond is already rounded.

    A value already on a 30-degree boundary is left alone as well, which is
    what makes the hold below the boundary a one-sided rule.
    """
    whole = cs - cs % CS_PER_ARCSEC
    assert csroundsec(whole) == whole


def test_csroundsec_is_not_an_odd_function() -> None:
    """Property 21, second half: the negative side is not the mirror image."""
    assert csroundsec(12345) == 12300
    assert csroundsec(-12345) == -12200
    assert csroundsec(-12345) != -csroundsec(12345)


@settings(max_examples=400, deadline=None)
@given(cs=st.integers(min_value=0, max_value=2**32))
def test_csroundsec_keeps_a_reading_inside_its_sign(cs: int) -> None:
    """Property 22: rounding never promotes an angle into the next sign.

    Stated for the non-negative domain, which is the one a formatter reads:
    the negative half rounds toward positive infinity by contract, so an
    angle just below zero legitimately comes back as zero.
    """
    assert csroundsec(cs) // CS_PER_SIGN == cs // CS_PER_SIGN
