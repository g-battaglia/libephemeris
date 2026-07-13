"""Targeted regression tests for the v3.0.0 round-2 review fixes.

Each test pins the corrected behaviour of one review finding so a future
regression is caught. Grouped by subsystem; see release-notes for the
finding list.
"""

from __future__ import annotations

import os

import pytest

import libephemeris as le
from libephemeris import SUN
from libephemeris.extinction import calc_airmass
from libephemeris.fixed_stars import resolve_star_name
from libephemeris.state import _resolve_data_dir


class TestParityOneLiners:
    def test_pheno_sun_disc_is_fully_illuminated(self) -> None:
        """The emitting solar disc has zero phase angle and full illumination."""
        attr = le.pheno_ut(2451545.0, SUN, 0)
        assert attr[0] == 0.0  # phase angle
        assert attr[1] == 1.0

    def test_context_set_topo_alt_defaults_to_zero(self) -> None:
        """Context.set_topo parity with the module-level set_topo(alt=0.0)."""
        ctx = le.EphemerisContext()
        ctx.set_topo(12.5, 41.9)  # 2-arg form must not raise
        assert ctx.get_topo() is not None

    def test_regulus_flamsteed_alias_is_32_leo(self) -> None:
        """Regulus is 32 Leonis, not 87 Leonis."""
        assert resolve_star_name("32 Leonis") is not None
        assert resolve_star_name("87 Leonis") is None


class TestRobustness:
    def test_airmass_never_below_one_at_zenith(self) -> None:
        """Physical airmass floor is 1.0 at the zenith for every method."""
        assert calc_airmass(90.0, method="kasten_young") >= 1.0
        assert calc_airmass(90.0, method="rozenberg") >= 1.0
        assert calc_airmass(90.0, method="secant") >= 1.0

    def test_data_dir_expands_tilde_from_env(self) -> None:
        """A leading ~ in LIBEPHEMERIS_DATA_DIR is expanded, not literal."""
        prev = os.environ.get("LIBEPHEMERIS_DATA_DIR")
        os.environ["LIBEPHEMERIS_DATA_DIR"] = "~/customdata"
        try:
            resolved = _resolve_data_dir()
        finally:
            if prev is None:
                del os.environ["LIBEPHEMERIS_DATA_DIR"]
            else:
                os.environ["LIBEPHEMERIS_DATA_DIR"] = prev
        assert resolved == os.path.join(os.path.expanduser("~"), "customdata")
        assert "~" not in resolved


class TestStarResolution:
    def test_greek_single_letter_alias_resolves_in_all_resolvers(self) -> None:
        """'α Leo' must resolve to Regulus in resolve_star_name and _resolve_star2
        (the .upper() of a lowercase Greek letter differs from the key)."""
        from libephemeris.fixed_stars import _resolve_star2, resolve_star_name

        assert resolve_star_name("α Leo") == 1000001
        entry, err = _resolve_star2("α Leo")
        assert err is None
        assert entry is not None and entry.id == 1000001


class TestHypothetical:
    def test_undefined_fictitious_id_raises(self) -> None:
        """An id in the fictitious range but without an ephemeris must raise,
        not return a phantom (0,0,0) position."""
        from libephemeris.exceptions import UnknownBodyError
        from libephemeris.hypothetical import calc_hypothetical_position

        with pytest.raises(UnknownBodyError):
            calc_hypothetical_position(59, 2451545.0)


class TestTime:
    def test_date_conversion_reform_day_morning_not_shifted(self) -> None:
        """An instant in the 00:00-11:59 window of the Gregorian reform day
        (1582-10-15) must not be misclassified as Julian and shifted +10 days."""
        from libephemeris.time_utils import date_conversion

        year, month, day, hour = date_conversion(1582, 10, 15, 6.0, "g")
        assert (year, month, day) == (1582, 10, 15)

    def test_date_conversion_before_reform_is_julian(self) -> None:
        """A date before the reform is still treated as Julian input."""
        from libephemeris.time_utils import date_conversion

        year, month, day, hour = date_conversion(1582, 10, 4, 6.0, "g")
        # 1582-10-04 (Julian) -> 1582-10-14 (Gregorian)
        assert (year, month, day) == (1582, 10, 14)


class TestSplitDeg:
    def test_nakshatra_index_for_unnormalized_longitude(self) -> None:
        """Nakshatra indices use the longitude modulo one full circle."""
        nak = le.SPLIT_DEG_NAKSHATRA
        for angle in (0.0, 17.25, 359.9, 720.0):
            assert le.split_deg(angle, nak) == le.split_deg(angle + 360.0, nak)


class TestReaderResourceSafety:
    def test_zero_byte_leb_file_does_not_leak_fd(self, tmp_path) -> None:
        """A 0-byte stub raises without leaking the open file descriptor."""
        import gc
        import warnings

        from libephemeris.leb2_reader import LEB2Reader
        from libephemeris.leb_reader import LEBReader

        stub = tmp_path / "empty_base_core.leb"
        stub.write_bytes(b"")

        for reader_cls in (LEBReader, LEB2Reader):
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                with pytest.raises((ValueError, OSError)):
                    reader_cls(str(stub))
                gc.collect()
            assert not [w for w in caught if issubclass(w.category, ResourceWarning)]

    def test_process_global_leb_reader_closes_at_interpreter_exit(self) -> None:
        """Lazy global readers register an interpreter-shutdown finalizer."""
        import subprocess
        import sys
        from pathlib import Path

        leb_path = (
            Path(__file__).parents[1]
            / "libephemeris"
            / "data"
            / "leb2"
            / "base_core.leb2"
        )
        code = (
            "import libephemeris.state as s; "
            f"s.set_leb_file({str(leb_path)!r}); "
            "assert s.get_leb_reader() is not None"
        )
        result = subprocess.run(
            [sys.executable, "-W", "always::ResourceWarning", "-c", code],
            check=False,
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr
        assert "ResourceWarning" not in result.stderr


class TestCalcFlags:
    def test_ephemeris_bits_are_mutually_exclusive(self) -> None:
        """calc_ut never echoes two ephemeris bits at once.

        The public flag contract makes them mutually exclusive (JPLEPH >
        SWIEPH), so
        FLG_JPLEPH|FLG_SWIEPH must collapse to a single-bit retflag."""
        from libephemeris.constants import FLG_JPLEPH, FLG_MOSEPH, FLG_SWIEPH

        jd = 2451545.0
        for combo in (
            FLG_JPLEPH | FLG_SWIEPH,
            FLG_JPLEPH | FLG_SWIEPH | FLG_MOSEPH,
        ):
            _, retflag = le.calc_ut(jd, SUN, combo)
            assert not (retflag & FLG_JPLEPH and retflag & FLG_SWIEPH), (
                f"retflag {retflag} has both ephemeris bits set"
            )
            # JPLEPH has the documented public priority.
            assert retflag & FLG_JPLEPH and not (retflag & FLG_SWIEPH)

        # A lone ephemeris bit (or none) is preserved / defaulted as before.
        assert le.calc_ut(jd, SUN, FLG_JPLEPH)[1] & FLG_JPLEPH
        assert le.calc_ut(jd, SUN, 0)[1] & FLG_SWIEPH


class TestHeliacalMetRange:
    @pytest.mark.slow
    def test_heliacal_ut_honours_atmo_turbidity(self) -> None:
        """heliacal_ut must thread atmo[3] (meteorological range) into the
        event search: a hazy atmosphere delays first visibility relative to a
        clear one. Previously atmo[3] was dropped, so the date was identical."""
        jd0 = le.julday(2000, 1, 1, 0.0)
        geopos = (12.5, 41.9, 0.0)
        dobs = (36.0, 1.0, 1.0, 1.0, 1.0, 1.0)

        def rising(met_range: float) -> float:
            atmo = (1013.25, 15.0, 40.0, met_range)
            jd1, _, _ = le.heliacal_ut(
                jd0, geopos, atmo, dobs, "Venus", le.HELIACAL_RISING
            )
            return jd1

        clear = rising(40.0)  # 40 km meteorological range
        hazy = rising(5.0)  # 5 km -> much more extinction
        assert clear > 0 and hazy > 0
        # Hazier air pushes the heliacal rising later (and by more than a day,
        # so this cannot be numerical noise).
        assert hazy > clear + 1.0


class TestHeliacalSearchGates:
    """Regression tests for the conjunction-aware event-search gates."""

    GEOPOS = (31.0, 30.0, 0.0)  # Cairo
    ATMO = (1013.25, 15.0, 40.0, 0.0)
    DOBS = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    @pytest.mark.slow
    def test_evening_first_from_inside_conjunction_does_not_skip_apparition(
        self,
    ) -> None:
        """Starting the evening-first search while Venus is deep inside the
        2024 superior-conjunction gap (elongation ~5 deg, invisible) must
        return the imminent evening-first event, not the one a full synodic
        period (~584 days) later."""
        jd0 = le.julday(2024, 6, 21, 0.0)  # elong ~5.6 deg, invisible
        jd1, _, _ = le.heliacal_ut(
            jd0, self.GEOPOS, self.ATMO, self.DOBS, "Venus", le.EVENING_FIRST
        )
        assert jd1 > jd0
        assert jd1 - jd0 < 60.0  # the imminent apparition, not +584 d

    @pytest.mark.slow
    def test_evening_first_mid_apparition_returns_next_event(self) -> None:
        """Starting the search AFTER Venus is already visible as an evening
        star (ongoing apparition, low but climbing elongation) must return the
        NEXT evening-first — after the coming conjunction — not an ordinary
        visibility day of the current apparition."""
        # First resolve the current apparition's evening-first from inside
        # the conjunction gap, then start a few days after it.
        jd_gap = le.julday(2024, 6, 21, 0.0)
        ef, _, _ = le.heliacal_ut(
            jd_gap, self.GEOPOS, self.ATMO, self.DOBS, "Venus", le.EVENING_FIRST
        )
        assert ef > 0
        jd_mid = ef + 4.0  # visible evening star, elongation still < 10 deg
        nxt, _, _ = le.heliacal_ut(
            jd_mid, self.GEOPOS, self.ATMO, self.DOBS, "Venus", le.EVENING_FIRST
        )
        # The next true evening-first is a synodic period away (> 300 days),
        # never a re-report of the ongoing apparition (< 60 days).
        assert nxt - jd_mid > 300.0

    @pytest.mark.slow
    def test_rising_extended_lookback_blocks_wrong_apparition(self) -> None:
        """Mercury heliacal rising searched from 2024-01-01 Rome: the body is
        invisible at dawn with elongation ~18 deg. The extended lookback must
        keep the pre-window-gap clause from accepting a wrong apparition
        months out; the result must be the nearby one."""
        jd0 = le.julday(2024, 1, 1, 0.0)
        jd1, _, _ = le.heliacal_ut(
            jd0,
            (12.5, 41.9, 0.0),  # Rome
            self.ATMO,
            self.DOBS,
            "Mercury",
            le.HELIACAL_RISING,
        )
        assert jd1 > 0
        assert jd1 - jd0 < 90.0  # not the +243 d apparition


class TestAssistHeliocentricFrame:
    """Regression: the ASSIST minor-body heliocentric branch must return an
    of-date ecliptic longitude, not one frozen at the J2000 equinox.

    Beyond SPK and LEB coverage a minor body is propagated with ASSIST, whose
    state is heliocentric ecliptic J2000. Every other body path (and the
    geocentric branch of the same helper) precesses that to the ecliptic of
    date; the heliocentric branch used to return the raw J2000 longitude,
    drifting from every other path by the accumulated general precession
    (~7 deg by year 2500). Uses a public-API, reference-free invariant.
    """

    @staticmethod
    def _sph2cart(lon: float, lat: float, r: float) -> tuple[float, float, float]:
        import math

        lo, la = math.radians(lon), math.radians(lat)
        return (
            r * math.cos(la) * math.cos(lo),
            r * math.cos(la) * math.sin(lo),
            r * math.sin(la),
        )

    @pytest.mark.slow
    def test_chiron_heliocentric_is_of_date_beyond_coverage(self) -> None:
        """Year 2500 is past Chiron's SPK/LEB coverage, so it is propagated by
        ASSIST (or the Keplerian fallback). The of-date heliocentric longitude
        must differ from the J2000 one by the real precession since J2000
        (~7 deg), never be equal to it (equal == frozen-at-J2000 defect)."""
        from libephemeris import CHIRON, FLG_HELCTR, FLG_J2000

        jd = le.julday(2500, 1, 1, 0.0)
        of_date, _ = le.calc_ut(jd, CHIRON, FLG_HELCTR)
        j2000, _ = le.calc_ut(jd, CHIRON, FLG_HELCTR | FLG_J2000)
        precession = abs((of_date[0] - j2000[0] + 180.0) % 360.0 - 180.0)
        # Precession over 500 years is ~7 deg; the frozen bug gave ~0.
        assert precession > 5.0

    @pytest.mark.slow
    def test_chiron_heliocentric_reconstructs_geocentric(self) -> None:
        """Reference-free coherence: the geometric geocentric longitude must be
        reconstructable from the heliocentric position plus Earth's
        heliocentric position (Earth_helio = -Sun_geo), all of-date ecliptic.
        Before the fix this failed by ~7 deg in the ASSIST region."""
        import math

        from libephemeris import CHIRON, SUN, FLG_HELCTR, FLG_TRUEPOS

        jd = le.julday(2500, 1, 1, 0.0)
        hel, _ = le.calc_ut(jd, CHIRON, FLG_HELCTR)  # of-date ecliptic
        geo, _ = le.calc_ut(jd, CHIRON, FLG_TRUEPOS)  # geometric geocentric
        sun, _ = le.calc_ut(jd, SUN, FLG_TRUEPOS)  # Sun geocentric = -Earth_helio
        ah = self._sph2cart(hel[0], hel[1], hel[2])
        sg = self._sph2cart(sun[0], sun[1], sun[2])
        recon = [ah[i] + sg[i] for i in range(3)]  # geo = ast_helio + Sun_geo
        recon_lon = math.degrees(math.atan2(recon[1], recon[0])) % 360.0
        dlon = abs((recon_lon - geo[0] + 180.0) % 360.0 - 180.0)
        # Exact when both sources are geometric (ASSIST region); the frozen bug
        # produced ~7 deg. A generous bound stays robust across environments.
        assert dlon < 0.05


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
