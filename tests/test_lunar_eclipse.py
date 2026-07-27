"""
Tests for lunar eclipse calculations in libephemeris.

Tests the lun_eclipse_when function which finds lunar eclipses
and calculates their phase times.

Reference data from NASA Eclipse website:
https://eclipse.gsfc.nasa.gov/lunar.html

Times tuple layout (10 elements, reference-API compatible):
    [0]: Time of maximum eclipse
    [1]: Reserved
    [2]: Time of partial eclipse beginning (Moon enters umbra)
    [3]: Time of partial eclipse ending (Moon leaves umbra)
    [4]: Time of totality beginning
    [5]: Time of totality ending
    [6]: Time of penumbral eclipse beginning
    [7]: Time of penumbral eclipse ending
    [8]: Reserved
    [9]: Reserved
"""

import pytest

pytestmark = pytest.mark.slow

from libephemeris import (
    julday,
    revjul,
    lun_eclipse_when,
    ECL_TOTAL,
    ECL_PARTIAL,
    ECL_PENUMBRAL,
)


class TestLunEclipseWhen:
    """Test suite for lun_eclipse_when function."""

    def test_finds_lunar_eclipse(self):
        """Test that function finds a lunar eclipse."""
        # Start from Jan 1, 2024
        jd_start = julday(2024, 1, 1, 0)

        ecl_type, times = lun_eclipse_when(jd_start)

        # Should find an eclipse
        assert ecl_type != 0
        assert times[0] > jd_start  # Maximum should be after start
        assert times[0] > 0  # Should have valid maximum time

    def test_returns_correct_tuple_structure(self):
        """Test that return values have correct structure."""
        jd_start = julday(2024, 1, 1, 0)

        ecl_type, times = lun_eclipse_when(jd_start)

        # Should return 10-element tuple (reference-API layout)
        assert len(times) == 10
        # All elements should be floats
        assert all(isinstance(t, float) for t in times)
        # Eclipse type should be int
        assert isinstance(ecl_type, int)

    def test_phase_times_ordering(self):
        """Test that phase times are in correct chronological order."""
        jd_start = julday(2024, 1, 1, 0)

        ecl_type, times = lun_eclipse_when(jd_start)

        # Penumbral begin should be first (if present)
        if times[6] > 0:
            assert times[6] < times[0]  # Penumbral begin < maximum

        # Penumbral end should be last (if present)
        if times[7] > 0:
            assert times[7] > times[0]  # Penumbral end > maximum

        # Partial phases should be between penumbral phases
        if times[2] > 0 and times[6] > 0:
            assert times[2] >= times[6]  # Partial begin >= penumbral begin
        if times[3] > 0 and times[7] > 0:
            assert times[3] <= times[7]  # Partial end <= penumbral end

        # Total phases should be between partial phases
        if times[4] > 0 and times[2] > 0:
            assert times[4] >= times[2]  # Total begin >= partial begin
        if times[5] > 0 and times[3] > 0:
            assert times[5] <= times[3]  # Total end <= partial end

    def test_filter_total_eclipse(self):
        """Test filtering for total lunar eclipses."""
        jd_start = julday(2020, 1, 1, 0)

        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)

        # Should find a total eclipse
        assert ecl_type & ECL_TOTAL
        # Total phase times should be present
        assert times[4] > 0  # Total begin
        assert times[5] > 0  # Total end

    def test_filter_partial_eclipse(self):
        """Test filtering for partial lunar eclipses."""
        jd_start = julday(2020, 1, 1, 0)

        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_PARTIAL)

        # Should find a partial eclipse
        assert ecl_type & ECL_PARTIAL
        # Partial phase times should be present
        assert times[2] > 0  # Partial begin
        assert times[3] > 0  # Partial end

    def test_filter_penumbral_eclipse(self):
        """Test filtering for penumbral lunar eclipses."""
        jd_start = julday(2020, 1, 1, 0)

        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_PENUMBRAL)

        # Should find a penumbral eclipse
        assert ecl_type & ECL_PENUMBRAL
        # Penumbral phase times should be present
        assert times[6] > 0  # Penumbral begin
        assert times[7] > 0  # Penumbral end

    def test_known_total_eclipse_may_2022(self):
        """Test finding the known total lunar eclipse of May 16, 2022.

        Reference: NASA - Total Lunar Eclipse of May 16, 2022
        Maximum eclipse: approximately 04:12 UTC
        """
        # Start searching from May 1, 2022
        jd_start = julday(2022, 5, 1, 0)

        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)

        # Should be a total eclipse
        assert ecl_type & ECL_TOTAL

        # Check maximum is on May 16, 2022
        year, month, day, hour = revjul(times[0])
        assert year == 2022
        assert month == 5
        assert day == 16
        # Maximum should be around 04:12 UTC (within 30 minutes tolerance)
        assert 3.5 < hour < 5.0

    def test_known_partial_eclipse_oct_2023(self):
        """Test finding the known partial lunar eclipse of October 28, 2023.

        Reference: NASA - Partial Lunar Eclipse of October 28, 2023
        Maximum eclipse: approximately 20:14 UTC
        """
        # Start searching from October 1, 2023
        jd_start = julday(2023, 10, 1, 0)

        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_PARTIAL)

        # Should be a partial eclipse
        assert ecl_type & ECL_PARTIAL

        # Check maximum is on October 28, 2023
        year, month, day, hour = revjul(times[0])
        assert year == 2023
        assert month == 10
        assert day == 28
        # Maximum should be around 20:14 UTC (within 30 minutes tolerance)
        assert 19.5 < hour < 21.0

    def test_swe_alias(self):
        """Test that lun_eclipse_when is an alias for lun_eclipse_when."""
        jd_start = julday(2024, 1, 1, 0)

        ecl_type1, times1 = lun_eclipse_when(jd_start)
        ecl_type2, times2 = lun_eclipse_when(jd_start)

        assert times1 == times2
        assert ecl_type1 == ecl_type2

    def test_sequential_eclipses(self):
        """Test finding multiple sequential lunar eclipses."""
        jd = julday(2024, 1, 1, 0)
        eclipses = []

        # Find 3 sequential eclipses
        for _ in range(3):
            ecl_type, times = lun_eclipse_when(jd)
            eclipses.append((times[0], ecl_type))
            jd = times[0] + 1  # Start after this eclipse

        # Each eclipse should be later than the previous
        for i in range(1, len(eclipses)):
            assert eclipses[i][0] > eclipses[i - 1][0]

        # Eclipses should be roughly 6 months apart (147-177 days)
        for i in range(1, len(eclipses)):
            diff = eclipses[i][0] - eclipses[i - 1][0]
            assert 140 < diff < 210  # Reasonable range for lunar eclipse intervals

    def test_penumbral_only_eclipse_has_no_umbral_phases(self):
        """Test that penumbral-only eclipses have zero umbral phase times."""
        jd_start = julday(2020, 1, 1, 0)

        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_PENUMBRAL)

        if ecl_type == ECL_PENUMBRAL:  # Pure penumbral, not partial
            # Umbral phase times should be zero
            assert times[2] == 0.0  # Partial begin
            assert times[3] == 0.0  # Partial end
            assert times[4] == 0.0  # Total begin
            assert times[5] == 0.0  # Total end
            # Penumbral times should be present
            assert times[6] > 0.0
            assert times[7] > 0.0


class TestLunEclipseEdgeCases:
    """Test edge cases for lunar eclipse calculations."""

    def test_early_20th_century(self):
        """Test finding eclipse in early 20th century."""
        # 1920 CE (within ephemeris range)
        jd_start = julday(1920, 1, 1, 0)

        ecl_type, times = lun_eclipse_when(jd_start)

        assert ecl_type != 0
        assert times[0] > jd_start

    def test_mid_21st_century(self):
        """Test finding eclipse in mid 21st century."""
        # 2050 CE (within ephemeris range)
        jd_start = julday(2050, 1, 1, 0)

        ecl_type, times = lun_eclipse_when(jd_start)

        assert ecl_type != 0
        assert times[0] > jd_start

    def test_eclipse_type_flags(self):
        """Test that eclipse type flags are set correctly."""
        jd_start = julday(2020, 1, 1, 0)

        # Total eclipse
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        assert ecl_type & ECL_TOTAL
        # Total eclipse should also have partial phase times
        assert times[2] > 0  # Partial begins
        assert times[3] > 0  # Partial ends


class TestLunEclipseWhenAdvanceMargin:
    """The `jd = tret[0]; when(jd)` idiom must advance, not stall.

    Measured reference behavior: a maximum within 1e-4 day (8.64 s) of the
    search epoch counts as already reached. A lunar eclipse's opposition (Full
    Moon) precedes its maximum by a few minutes, so the forward search also
    pre-probes the previous lunation to keep an in-progress eclipse whose Full
    Moon is already past but whose maximum is still ahead.
    """

    def test_start_on_maximum_advances_forward(self):
        _, t0 = lun_eclipse_when(julday(1968, 1, 1, 0))
        _, t1 = lun_eclipse_when(t0[0])
        assert t1[0] - t0[0] > 1.0

    def test_start_on_maximum_advances_backward(self):
        _, t0 = lun_eclipse_when(julday(1968, 1, 1, 0))
        _, tb = lun_eclipse_when(t0[0], backwards=True)
        assert t0[0] - tb[0] > 1.0

    def test_iteration_does_not_stall(self):
        jd = julday(1990, 1, 1, 0)
        prev = None
        for _ in range(8):
            _, t = lun_eclipse_when(jd)
            if prev is not None:
                assert t[0] - prev > 1.0
            prev = t[0]
            jd = t[0]

    def test_pre_probe_keeps_eclipse_between_opposition_and_maximum(self):
        # The 2003-05-16 eclipse's Full Moon precedes its maximum by ~4 min.
        # Starting a few minutes before the maximum (past the Full Moon) must
        # still return this eclipse, not skip forward to the next lunation.
        _, t0 = lun_eclipse_when(julday(2003, 5, 1, 0))
        jm = t0[0]
        _, cur = lun_eclipse_when(jm - 100.0 / 86400.0)  # 100 s before maximum
        assert abs(cur[0] - jm) < 0.5
