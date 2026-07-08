"""Regression tests: ``cross_ut`` must return the FIRST crossing for slow planets.

For a slow outer planet (Jupiter..Pluto) with a distant target, Newton-Raphson
seeded from a single linear estimate ``diff / speed`` can converge into a
crossing basin one or more retrograde loops PAST the first true crossing — the
net progress of the planet is retarded by its retrograde loops, so the seed
under-/over-shoots, and Newton then settles on a later crossing (or on the
2nd/3rd of a retrograde triple-crossing). The post-Newton guards only checked
">= tjdut" and "within max_range", never "earliest", so the contract "first
crossing at or after tjdut" was silently violated.

Three deterministic real cases (Saturn) that returned a later crossing:

  R1  target 318.248     start 2432957.5           returned the 3rd of a
      retrograde triple (+5413 d) instead of the 1st (+5145 d) — a 269 d skip.
  R2  target 222.706     start 2472807.6           a 263 d skip (same pattern).
  R3  target 84.35914... start 2446732.08028656    on the LEB backend, a whole
      Saturn revolution too far (+16652 d vs the true +5751 d) — a 10901 d skip.

The fix (``_guarded_first_crossing`` in crossing.py) re-validates a slow planet's
far Newton result against an earlier crossing via a provably skip-free adaptive
forward scan. These tests pin the first-crossing contract; they are backend
agnostic (they compare against an independent scan of the library's *own*
longitude curve, so they hold identically on the Skyfield and LEB backends).
"""

from __future__ import annotations

import pytest

import libephemeris as L
from libephemeris.constants import (
    JUPITER,
    NEPTUNE,
    SATURN,
    URANUS,
    FLG_SWIEPH,
)


def _lon(planet: int, jd: float) -> float:
    return L.calc_ut(jd, planet, FLG_SWIEPH)[0][0]


def _brute_first(
    planet: int, start: float, target: float, window: float, step: float = 2.0
):
    """First genuine sign-change crossing of the library's own longitude curve,
    via a fine forward scan + bisection. The independent oracle cross_ut must
    reproduce. Returns None if no crossing lies within ``window`` days."""
    prev = _lon(planet, start)
    jd = start
    while jd < start + window:
        jd += step
        cur = _lon(planet, jd)
        d1 = ((target - prev + 180.0) % 360.0) - 180.0
        d2 = ((target - cur + 180.0) % 360.0) - 180.0
        if (d1 < 0) != (d2 < 0) and abs(d2 - d1) < 180.0:
            lo, hi = jd - step, jd
            f_lo = ((target - _lon(planet, lo) + 180.0) % 360.0) - 180.0
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                f_mid = ((target - _lon(planet, mid) + 180.0) % 360.0) - 180.0
                if (f_lo < 0) == (f_mid < 0):
                    lo, f_lo = mid, f_mid
                else:
                    hi = mid
            return 0.5 * (lo + hi)
        prev = cur
    return None


# (name, planet, target, start, first-crossing-jd, old-buggy-jd)
_REPROS = [
    ("R1", SATURN, 318.248, 2432957.5, 2438102.3682, 2438370.9325),
    ("R2", SATURN, 222.706, 2472807.6, 2477859.2591, 2478122.7359),
    ("R3", SATURN, 84.3591459768107, 2446732.08028656, 2452482.8195, 2463383.7599),
]


@pytest.mark.parametrize("name,planet,target,start,first_jd,old_jd", _REPROS)
def test_repro_returns_first_crossing(name, planet, target, start, first_jd, old_jd):
    """Each reported case returns the FIRST crossing, not the later one Newton
    used to land on — on whichever backend the suite runs."""
    got = L.cross_ut(planet, target, start, FLG_SWIEPH)

    # Forward-only contract.
    assert got >= start - 1e-6, f"{name}: went backward by {start - got:.1f} d"

    # It is a genuine crossing (longitude equals the target).
    resid = abs(((_lon(planet, got) - target + 180.0) % 360.0) - 180.0)
    assert resid * 3600.0 < 1.0, f'{name}: not a crossing, {resid * 3600:.3f}"'

    # It is the FIRST crossing (matches an independent fine scan), NOT the later
    # basin Newton previously converged to (>= 260 d / a whole revolution later).
    assert abs(got - first_jd) < 1.0, (
        f"{name}: got +{got - start:.1f} d, expected first +{first_jd - start:.1f} d"
    )
    assert got < old_jd - 100.0, (
        f"{name}: returned the pre-fix late crossing (+{got - start:.1f} d)"
    )


@pytest.mark.parametrize("name,planet,target,start,first_jd,old_jd", _REPROS)
def test_repro_no_earlier_crossing(name, planet, target, start, first_jd, old_jd):
    """Independent confirmation: an oracle fine scan of the library's own curve
    agrees the returned crossing is the earliest one at or after the start."""
    got = L.cross_ut(planet, target, start, FLG_SWIEPH)
    truth = _brute_first(planet, start, target, window=(got - start) + 400.0)
    assert truth is not None
    assert abs(got - truth) < 1.0, (
        f"{name}: cross_ut +{got - start:.1f} d vs oracle +{truth - start:.1f} d"
    )


def test_forward_first_is_strictly_monotone():
    """Chaining forward from a result finds strictly later crossings each time
    (no oscillation back onto the same or an earlier basin)."""
    start = 2432957.5
    target = 318.248
    prev = L.cross_ut(SATURN, target, start, FLG_SWIEPH)
    seq = [prev]
    for _ in range(10):
        nxt = L.cross_ut(SATURN, target, prev + 1.0, FLG_SWIEPH)
        assert nxt > prev + 1e-6, f"not monotone: {nxt} <= {prev}"
        seq.append(nxt)
        prev = nxt
    # Each successive Saturn crossing of a fixed longitude is ~1 sidereal period
    # apart on average; strictly increasing is the contract we assert.
    assert seq == sorted(seq)


@pytest.mark.slow
@pytest.mark.parametrize(
    "planet,syn,lead_cap,seed",
    [
        (SATURN, 378.0, 12000.0, 1001),
        (JUPITER, 399.0, 9000.0, 1002),
        (NEPTUNE, 367.0, 12000.0, 1003),
        (URANUS, 370.0, 12000.0, 1004),
    ],
)
def test_slow_planet_far_targets_are_first_crossing(planet, syn, lead_cap, seed):
    """Breadth guard: for many distant targets (several synodic loops ahead —
    the geometry that trips Newton into a wrong loop), cross_ut reproduces the
    first crossing found by an independent fine scan. Targets are the planet's
    own longitude at a bounded future epoch, so they span 0-360 while keeping the
    oracle window tractable."""
    import random

    rng = random.Random(seed)
    starts = [2400000.5 + i * 3500.0 for i in range(40)]
    checked = 0
    for _ in range(24):
        start = rng.choice(starts) + rng.uniform(0.0, 3000.0)
        lead = rng.uniform(2.0 * syn, lead_cap)
        target = _lon(planet, start + lead) % 360.0
        got = L.cross_ut(planet, target, start, FLG_SWIEPH)
        assert got >= start - 1e-6
        window = min((got - start) + 500.0, lead + 800.0)
        truth = _brute_first(planet, start, target, window=window)
        if truth is None:
            continue
        checked += 1
        assert abs(got - truth) < 2.0, (
            f"planet {planet}: got +{got - start:.1f} d vs first +{truth - start:.1f} d"
        )
    assert checked > 0
