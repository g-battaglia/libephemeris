"""Regression tests for root-finder fixes in ``libephemeris.crossing``.

These pin three previously-broken behaviours:

1. ``helio_cross_ut`` / ``helio_cross`` forward search: an imminent forward
   crossing (planet a hair before, or exactly at, the target) must be returned
   as "now", not one full orbital period later. A forward dead-band used to push
   such a case a whole period into the future.

2. ``helio_cross_ut`` / ``helio_cross`` backward search: for a slow eccentric
   body (e.g. Pluto) the returned crossing must be the NEAREST one before the
   start time, not the oldest crossing inside the search window.

3. ``get_station_info``: the docstring promises the nearest station in either
   direction with a negative ``days_to_station`` when it is in the past. A
   just-passed station must therefore be reported with a negative value.

The reference used here is an independent brute-force scan of the library's own
heliocentric longitude / velocity (no external ephemeris dependency): the fixed
functions must agree with a direct root search of the same underlying model.
"""

from __future__ import annotations

import libephemeris as le
from libephemeris import crossing

FLG = le.FLG_SWIEPH
HELIO_FLAG = FLG | le.FLG_HELCTR | le.FLG_SPEED

# Approximate heliocentric orbital periods (days) for "is it a period later?"
# assertions. Only the order of magnitude matters here.
ORBITAL_PERIOD_DAYS = {
    le.MERCURY: 88.0,
    le.VENUS: 224.7,
    le.EARTH: 365.25,
    le.MARS: 687.0,
    le.JUPITER: 4332.6,
    le.NEPTUNE: 60190.0,
    le.PLUTO: 90560.0,
}


def _helio_lon(planet: int, jd: float) -> float:
    """Heliocentric ecliptic longitude (deg) from the library model."""
    pos, _ = le.calc_ut(jd, planet, HELIO_FLAG)
    return pos[0]


def _signed_diff(lon: float, target: float) -> float:
    """Signed angular distance lon - target wrapped to (-180, 180]."""
    d = (lon - target) % 360.0
    if d > 180.0:
        d -= 360.0
    return d


def _brute_nearest_past_crossing(
    planet: int,
    target: float,
    tjd: float,
    max_back_days: float,
    coarse: float = 300.0,
) -> float | None:
    """Nearest crossing of ``target`` strictly before ``tjd`` (brute force).

    Scans backward from ``tjd`` in ``coarse``-day steps for the first velocity-
    independent longitude sign change, then bisects to refine. Returns ``None``
    if none is found within ``max_back_days``.
    """
    prev = _signed_diff(_helio_lon(planet, tjd), target)
    steps = int(max_back_days / coarse)
    for i in range(1, steps + 1):
        cj = tjd - i * coarse
        cd = _signed_diff(_helio_lon(planet, cj), target)
        # Genuine crossing: sign change without an antipodal (~360°) wrap.
        if prev * cd <= 0 and abs(prev - cd) < 180.0:
            a, b = cj, cj + coarse
            fa = _signed_diff(_helio_lon(planet, a), target)
            for _ in range(70):
                m = 0.5 * (a + b)
                fm = _signed_diff(_helio_lon(planet, m), target)
                if fa * fm <= 0:
                    b = m
                else:
                    a, fa = m, fm
            return 0.5 * (a + b)
        prev = cd
    return None


# ---------------------------------------------------------------------------
# Bug 1: imminent forward heliocentric crossing must be "now", not a period on
# ---------------------------------------------------------------------------


class TestForwardImminentHelioCrossing:
    """A forward crossing a hair ahead must return ~now, not one period later."""

    def _check(self, func, planet: int, jd0: float) -> None:
        target = (_helio_lon(planet, jd0) + 20.0) % 360.0
        # The library's own next forward crossing of ``target``.
        jd_cross = func(planet, target, jd0, FLG)
        period = ORBITAL_PERIOD_DAYS[planet]

        # Start a hair (~0.086 s) before that crossing: the planet is now a
        # hair before the target, so an imminent forward crossing exists.
        start = jd_cross - 1e-6
        result = func(planet, target, start, FLG)

        # Must return the imminent crossing, not ~one orbital period later.
        assert abs(result - jd_cross) < 1.0, (
            f"planet {planet}: expected imminent crossing near {jd_cross}, "
            f"got {result} (dt={result - start:.3f} d)"
        )
        assert (result - start) < 0.5 * period, (
            f"planet {planet}: forward search jumped ~a period ahead "
            f"(dt={result - start:.1f} d, period={period} d)"
        )

    def test_helio_cross_ut_forward_imminent(self):
        jd0 = le.julday(2020, 6, 1, 0.0, le.GREG_CAL)
        # Earth/Mars/Venus exercise the Newton-Raphson path where the old
        # forward dead-band shifted the result a full period.
        for planet in (le.EARTH, le.MARS, le.VENUS, le.JUPITER):
            self._check(crossing.helio_cross_ut, planet, jd0)

    def test_helio_cross_tt_forward_imminent(self):
        jd0 = le.julday(2020, 6, 1, 0.0, le.GREG_CAL)
        for planet in (le.EARTH, le.MARS, le.VENUS, le.JUPITER):
            self._check(crossing.helio_cross, planet, jd0)


# ---------------------------------------------------------------------------
# Bug 2: backward heliocentric crossing must be the NEAREST past one
# ---------------------------------------------------------------------------


class TestBackwardNearestPastHelioCrossing:
    """Backward search on a slow eccentric body returns the nearest past crossing."""

    # (planet, year, longitude offset ahead of current position). Offsets ahead
    # of the body make the nearest past crossing lie almost a full orbit back,
    # which used to fall inside a window that returned the OLDER crossing.
    CASES = [
        (le.PLUTO, 2114, 180.0),
        (le.PLUTO, 2020, 90.0),
        (le.NEPTUNE, 2114, 0.5),
        (le.NEPTUNE, 1866, 5.0),
    ]

    def _check(self, func, planet: int, year: int, offset: float) -> None:
        jd0 = le.julday(year, 6, 1, 0.0, le.GREG_CAL)
        target = (_helio_lon(planet, jd0) + offset) % 360.0
        period = ORBITAL_PERIOD_DAYS[planet]

        result = func(planet, target, jd0, FLG, backwards=True)
        expected = _brute_nearest_past_crossing(planet, target, jd0, period * 2.5)
        assert expected is not None, "brute-force reference found no crossing"

        # Must be strictly in the past and match the nearest past crossing.
        assert result < jd0
        assert abs(result - expected) < 1.0, (
            f"planet {planet} {year} off={offset}: expected nearest past "
            f"crossing {expected} (dt={expected - jd0:.1f} d), got {result} "
            f"(dt={result - jd0:.1f} d) — likely the older in-window crossing"
        )
        # Guard against silently returning the previous-orbit crossing.
        assert abs(result - expected) < 0.5 * period

    def test_helio_cross_ut_backward_nearest(self):
        for planet, year, offset in self.CASES:
            self._check(crossing.helio_cross_ut, planet, year, offset)

    def test_helio_cross_tt_backward_nearest(self):
        for planet, year, offset in self.CASES:
            self._check(crossing.helio_cross, planet, year, offset)


# ---------------------------------------------------------------------------
# Bug 3: get_station_info reports a just-passed station with negative days
# ---------------------------------------------------------------------------


class TestStationInfoNearestDirection:
    """get_station_info considers the nearest station in either direction."""

    PLANETS = [le.MERCURY, le.MARS, le.JUPITER, le.PLUTO]

    def test_just_passed_station_is_negative(self):
        jd_start = le.julday(2020, 1, 1, 0.0, le.GREG_CAL)
        for planet in self.PLANETS:
            jd_station, _ = crossing.find_station_ut(planet, jd_start, "any", FLG)

            # Query 2 days AFTER the station: the nearest station is just-passed.
            info = crossing.get_station_info(planet, jd_station + 2.0, FLG)
            assert info["jd_station"] is not None
            assert info["days_to_station"] < 0.0, (
                f"planet {planet}: just-passed station reported with "
                f"days_to_station={info['days_to_station']}"
            )
            assert abs(info["days_to_station"] - (-2.0)) < 0.5
            assert abs(info["jd_station"] - jd_station) < 0.5

    def test_upcoming_station_still_positive(self):
        # Regression guard: a station a couple days ahead must remain positive.
        jd_start = le.julday(2020, 1, 1, 0.0, le.GREG_CAL)
        for planet in self.PLANETS:
            jd_station, _ = crossing.find_station_ut(planet, jd_start, "any", FLG)
            info = crossing.get_station_info(planet, jd_station - 2.0, FLG)
            assert info["days_to_station"] > 0.0
            assert abs(info["days_to_station"] - 2.0) < 0.5
            assert abs(info["jd_station"] - jd_station) < 0.5
