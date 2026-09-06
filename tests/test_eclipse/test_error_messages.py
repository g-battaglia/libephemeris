"""The text of every error ``libephemeris.eclipse`` raises.

One test per ``raise`` statement in the module. Each pins the exception type
and the wording of the message: the type and the condition are the contract
callers rely on, and the wording is what a user reads in a traceback, so a
rewrite of either must show up here first.

The searches that walk to the ephemeris boundary would take minutes to fail
for real; they are made to fail on their first step by replacing the lunation
or conjunction stepper they call, which exercises the same ``raise`` without
the wait.
"""

from __future__ import annotations

import math

import pytest

import libephemeris
from libephemeris import eclipse
from libephemeris.eclipse import (
    _coerce_backwards,
    _eclipse_sampling_step_days,
    _normalize_occultation_filter,
    _OccSearchOffEphemeris,
    _raise_if_sealed_leb_miss,
    _seek_moon_conjunction,
    _sol_eclipse_when_glob_pythonic,
    _sol_eclipse_when_loc_pythonic,
    _sol_glob_reject_impossible,
    get_inex_number,
    get_saros_number,
    lun_eclipse_how,
    lun_eclipse_when,
    lun_eclipse_when_loc,
    lun_occult_when_glob,
    lun_occult_when_loc,
    planet_occult_when_glob,
    planet_occult_when_loc,
    rise_trans,
    rise_trans_true_hor,
    sol_eclipse_how,
    sol_eclipse_how_details,
    sol_eclipse_magnitude_at_loc,
    sol_eclipse_max_time,
    sol_eclipse_obscuration_at_loc,
    sol_eclipse_when_glob,
    sol_eclipse_when_loc,
)
from libephemeris.exceptions import (
    EphemerisRangeError,
    Error,
    IllegalBodyError,
    InputValidationError,
    UnknownBodyError,
)

JD = 2451545.0
ROME = (12.5, 41.9, 0.0)
GEOPOS_TEXT = (
    "geopos must be a sequence of at least three numbers: geographic "
    "longitude and latitude in degrees, then altitude in metres"
)


def _off_ephemeris(*_args, **_kwargs):
    """Stand-in for a lunation stepper that has walked off the ephemeris."""
    raise EphemerisRangeError("probe stepped outside the ephemeris", requested_jd=JD)


def _no_eclipse(*_args, **_kwargs):
    """Stand-in for a global search that found nothing."""
    raise Error("no event")


# ---------------------------------------------------------------------------
# sealed-mode range and body misses
# ---------------------------------------------------------------------------
class TestSealedLebMiss:
    @pytest.fixture
    def sealed(self):
        from libephemeris import state

        previous = state.get_calc_mode()
        libephemeris.set_calc_mode("leb")
        try:
            yield
        finally:
            libephemeris.set_calc_mode(previous)

    def test_range_miss_keeps_the_reader_text_and_names_the_policy(self, sealed):
        cause = ValueError(
            "JD 2458300.5 outside range [2459945.5, 2461771.5] for body 0"
        )
        with pytest.raises(
            EphemerisRangeError,
            match=r"for body 0\. LEB mode does not silently substitute a "
            r"lower-precision source\.",
        ):
            _raise_if_sealed_leb_miss(cause)

    def test_body_miss_names_the_body_and_the_policy(self, sealed):
        with pytest.raises(
            UnknownBodyError,
            match=r"Body 3 is not stored in the active LEB file\. LEB mode does "
            r"not silently substitute a non-LEB source for a missing body\.",
        ):
            _raise_if_sealed_leb_miss(KeyError("Body 3 not in LEB file"))


# ---------------------------------------------------------------------------
# the illegal-body contract of the lunar occultation entry points
# ---------------------------------------------------------------------------
class TestIllegalBodyContract:
    def test_unknown_body_is_re_raised_as_the_contract_type(self):
        with pytest.raises(IllegalBodyError) as excinfo:
            lun_occult_when_glob(JD, 25, 2, 0)
        cause = excinfo.value.__cause__
        assert isinstance(cause, UnknownBodyError)
        assert not isinstance(cause, ValueError)
        assert str(excinfo.value) == str(cause)
        assert excinfo.value.body_id == 25


# ---------------------------------------------------------------------------
# eclipse and occultation type filters
# ---------------------------------------------------------------------------
class TestSolarTypeFilter:
    CENTRAL_PARTIAL = (
        r"The eclipse type filter asks for a central partial eclipse, which "
        r"cannot occur: a partial solar eclipse has no central line\."
    )
    NONCENTRAL_HYBRID = (
        r"The eclipse type filter asks for a non-central hybrid \(annular-total\) "
        r"eclipse, which cannot occur: a hybrid eclipse always has a central "
        r"line\."
    )

    def test_central_partial(self):
        with pytest.raises(Error, match=self.CENTRAL_PARTIAL):
            _sol_glob_reject_impossible(
                libephemeris.ECL_CENTRAL | libephemeris.ECL_PARTIAL
            )

    def test_central_partial_from_the_public_search(self):
        with pytest.raises(Error, match=self.CENTRAL_PARTIAL):
            sol_eclipse_when_glob(
                JD, 2, libephemeris.ECL_CENTRAL | libephemeris.ECL_PARTIAL, False
            )

    def test_noncentral_hybrid(self):
        with pytest.raises(Error, match=self.NONCENTRAL_HYBRID):
            _sol_glob_reject_impossible(
                libephemeris.ECL_NONCENTRAL | libephemeris.ECL_ANNULAR_TOTAL
            )

    def test_noncentral_hybrid_from_the_public_search(self):
        with pytest.raises(Error, match=self.NONCENTRAL_HYBRID):
            sol_eclipse_when_glob(
                JD,
                2,
                libephemeris.ECL_NONCENTRAL | libephemeris.ECL_ANNULAR_TOTAL,
                False,
            )


class TestOccultationTypeFilter:
    CENTRAL_PARTIAL = (
        r"The occultation type filter asks for a central partial event, which "
        r"cannot occur: a partial occultation or eclipse has no central line\."
    )
    ANNULAR = (
        r"The occultation type filter asks for an annular event of body 2, "
        r"which cannot occur: only the Sun can be eclipsed annularly\."
    )

    def test_central_partial(self):
        with pytest.raises(Error, match=self.CENTRAL_PARTIAL):
            _normalize_occultation_filter(
                libephemeris.ECL_PARTIAL | libephemeris.ECL_CENTRAL, 2, False
            )

    def test_central_partial_from_the_public_search(self):
        with pytest.raises(Error, match=self.CENTRAL_PARTIAL):
            lun_occult_when_glob(
                JD, 2, 2, libephemeris.ECL_PARTIAL | libephemeris.ECL_CENTRAL
            )

    def test_annular_for_a_non_solar_body(self):
        with pytest.raises(Error, match=self.ANNULAR):
            _normalize_occultation_filter(libephemeris.ECL_ANNULAR, 2, False)

    def test_annular_from_the_public_search(self):
        with pytest.raises(Error, match=self.ANNULAR):
            lun_occult_when_glob(JD, 2, 2, libephemeris.ECL_ANNULAR)


class TestLunarTypeFilter:
    def test_annular_lunar_eclipse(self):
        with pytest.raises(
            Error,
            match=r"The eclipse type filter asks for an annular lunar eclipse, "
            r"which cannot occur: only solar eclipses can be annular\.",
        ):
            lun_eclipse_when(JD, 2, libephemeris.ECL_ANNULAR)


# ---------------------------------------------------------------------------
# the conjunction seeker
# ---------------------------------------------------------------------------
class TestConjunctionSeek:
    def test_ephemeris_boundary_on_the_first_step_is_signalled(self):
        with pytest.raises(_OccSearchOffEphemeris) as excinfo:
            _seek_moon_conjunction(JD, 1, _off_ephemeris, 0, None)
        assert isinstance(excinfo.value.__cause__, EphemerisRangeError)

    def test_star_beyond_the_moons_reach(self):
        with pytest.raises(
            Error,
            match=r"The Moon never occults the star Vega: its ecliptic latitude "
            r"of 61\.7 degrees is beyond the Moon's reach\.",
        ):
            _seek_moon_conjunction(JD, 1, lambda _t: (279.0, 61.7), 0, "Vega")

    def test_star_beyond_the_moons_reach_from_the_public_search(self):
        with pytest.raises(
            Error, match=r"The Moon never occults the star Vega: its ecliptic latitude"
        ):
            lun_occult_when_glob(JD, "Vega", 2, 0)


# ---------------------------------------------------------------------------
# search-direction coercion
# ---------------------------------------------------------------------------
class TestSearchDirection:
    UNKNOWN_WORD = (
        r"'bogus' is not a valid search direction: pass a bool, an int, or one "
        r"of 'backward', 'backwards', 'true', '1' \(search backward\) or "
        r"'forward', 'forwards', 'false', '0', '' \(search forward\)\."
    )

    def test_unknown_word(self):
        with pytest.raises(TypeError, match=self.UNKNOWN_WORD):
            _coerce_backwards("bogus")

    def test_unknown_word_from_the_public_search(self):
        with pytest.raises(TypeError, match=self.UNKNOWN_WORD):
            sol_eclipse_when_glob(JD, 2, 0, "bogus")

    def test_wrong_type(self):
        with pytest.raises(
            TypeError,
            match=r"The search direction must be a bool, an int or a str, not "
            r"float\.",
        ):
            _coerce_backwards(1.5)

    def test_wrong_type_from_the_public_search(self):
        with pytest.raises(
            TypeError,
            match=r"The search direction must be a bool, an int or a str, not "
            r"NoneType\.",
        ):
            lun_eclipse_when(JD, 2, 0, None)


# ---------------------------------------------------------------------------
# sol_eclipse_max_time
# ---------------------------------------------------------------------------
class TestMaxTime:
    def test_half_a_location(self):
        with pytest.raises(
            ValueError,
            match=r"Pass both lat and lon for the local maximum, or neither for "
            r"the global maximum\.",
        ):
            sol_eclipse_max_time(JD, 45.0, None)


# ---------------------------------------------------------------------------
# global solar search
# ---------------------------------------------------------------------------
class TestSolarGlobalSearch:
    def test_unknown_search_direction(self):
        with pytest.raises(
            ValueError,
            match=r"search_direction must be one of \('forward', 'backward', "
            r"'bidirectional'\), not 'sideways'\.",
        ):
            _sol_eclipse_when_glob_pythonic(JD, search_direction="sideways")

    def test_backward_search_reaches_the_boundary(self, monkeypatch):
        monkeypatch.setattr(eclipse, "_find_previous_new_moon", _off_ephemeris)
        with pytest.raises(
            Error,
            match=r"No matching solar eclipse was found searching backward from "
            r"JD 2451545\.0 to the ephemeris boundary\.",
        ):
            sol_eclipse_when_glob(JD, 2, 0, True)

    def test_forward_search_reaches_the_boundary(self, monkeypatch):
        monkeypatch.setattr(eclipse, "_find_next_new_moon", _off_ephemeris)
        with pytest.raises(
            Error,
            match=r"No matching solar eclipse was found searching forward from "
            r"JD 2451545\.0 to the ephemeris boundary\.",
        ):
            sol_eclipse_when_glob(JD, 2, 0, False)


# ---------------------------------------------------------------------------
# local solar search, legacy (lat, lon) signature
# ---------------------------------------------------------------------------
class TestSolarLocalSearchLegacy:
    NOT_VISIBLE = (
        r"No solar eclipse visible from latitude 41\.9, longitude 12\.5 was "
        r"found within 50 years after JD 2451545\.0\."
    )

    def test_no_global_eclipse_at_all(self, monkeypatch):
        monkeypatch.setattr(eclipse, "_sol_eclipse_when_glob_pythonic", _no_eclipse)
        with pytest.raises(Error, match=self.NOT_VISIBLE):
            _sol_eclipse_when_loc_pythonic(JD, 41.9, 12.5)

    def test_no_eclipse_visible_from_the_place(self, monkeypatch):
        monkeypatch.setattr(
            eclipse,
            "_sol_eclipse_when_glob_pythonic",
            lambda jd, *a, **k: (libephemeris.ECL_TOTAL, (jd + 10.0,) + (0.0,) * 9),
        )
        monkeypatch.setattr(
            eclipse, "_calculate_local_eclipse_phases", lambda *a, **k: (0.0,) * 8
        )
        with pytest.raises(Error, match=self.NOT_VISIBLE):
            _sol_eclipse_when_loc_pythonic(JD, 41.9, 12.5)


# ---------------------------------------------------------------------------
# local solar search, geopos signature
# ---------------------------------------------------------------------------
class TestSolarLocalSearch:
    def test_geopos_too_short(self):
        with pytest.raises(ValueError, match=GEOPOS_TEXT):
            sol_eclipse_when_loc(JD, (12.0,))

    def test_forward_search_finds_no_global_eclipse(self, monkeypatch):
        monkeypatch.setattr(eclipse, "_sol_eclipse_when_glob_pythonic", _no_eclipse)
        with pytest.raises(
            Error,
            match=r"No solar eclipse visible from longitude 12\.5, latitude 41\.9 "
            r"was found within 50 years after JD 2451545\.0\.",
        ):
            sol_eclipse_when_loc(JD, ROME, 2, False)

    def test_backward_search_exhausts_its_window(self, monkeypatch):
        monkeypatch.setattr(eclipse, "_sol_eclipse_when_glob_pythonic", _no_eclipse)
        with pytest.raises(
            Error,
            match=r"No solar eclipse visible from longitude 12\.5, latitude 41\.9 "
            r"was found within 50 years before JD 2451545\.0\.",
        ):
            sol_eclipse_when_loc(JD, ROME, 2, True)


# ---------------------------------------------------------------------------
# the observer-position argument of every local entry point
# ---------------------------------------------------------------------------
class TestGeopos:
    def test_sol_eclipse_how(self):
        with pytest.raises(ValueError, match=GEOPOS_TEXT):
            sol_eclipse_how(JD, (12.0,))

    def test_sol_eclipse_how_details(self):
        with pytest.raises(ValueError, match=GEOPOS_TEXT):
            sol_eclipse_how_details(JD, (12.0,))

    def test_lun_eclipse_when_loc(self):
        with pytest.raises(ValueError, match=GEOPOS_TEXT):
            lun_eclipse_when_loc(JD, (12.0,))

    def test_lun_eclipse_how(self):
        with pytest.raises(ValueError, match=GEOPOS_TEXT):
            lun_eclipse_how(JD, (12.0,))

    def test_lun_occult_when_loc(self):
        with pytest.raises(ValueError, match=GEOPOS_TEXT):
            lun_occult_when_loc(JD, 2, (12.0,))

    def test_rise_trans_too_short(self):
        with pytest.raises(ValueError, match=GEOPOS_TEXT):
            rise_trans(JD, 0, libephemeris.CALC_RISE, (12.0,))

    def test_rise_trans_not_numeric(self):
        with pytest.raises(ValueError, match=GEOPOS_TEXT) as excinfo:
            rise_trans(JD, 0, libephemeris.CALC_RISE, ("a", "b", "c"))
        assert isinstance(excinfo.value.__cause__, ValueError)

    def test_rise_trans_true_hor(self):
        with pytest.raises(ValueError, match=GEOPOS_TEXT):
            rise_trans_true_hor(JD, 0, libephemeris.CALC_RISE, (12.0,))

    def test_sol_eclipse_magnitude_at_loc(self):
        with pytest.raises(ValueError, match=GEOPOS_TEXT):
            sol_eclipse_magnitude_at_loc(JD, (12.0,))

    def test_sol_eclipse_obscuration_at_loc(self):
        with pytest.raises(ValueError, match=GEOPOS_TEXT):
            sol_eclipse_obscuration_at_loc(JD, (12.0,))


# ---------------------------------------------------------------------------
# lunar searches
# ---------------------------------------------------------------------------
class TestLunarSearch:
    def test_forward_search_reaches_the_boundary(self, monkeypatch):
        monkeypatch.setattr(eclipse, "_find_next_full_moon", _off_ephemeris)
        with pytest.raises(
            Error,
            match=r"No matching lunar eclipse was found searching forward from "
            r"JD 2451545\.0 to the ephemeris boundary\.",
        ):
            lun_eclipse_when(JD, 2, 0, False)

    def test_backward_search_reaches_the_boundary(self, monkeypatch):
        monkeypatch.setattr(eclipse, "_find_previous_full_moon", _off_ephemeris)
        with pytest.raises(
            Error,
            match=r"No matching lunar eclipse was found searching backward from "
            r"JD 2451545\.0 to the ephemeris boundary\.",
        ):
            lun_eclipse_when(JD, 2, 0, True)

    def test_local_search_exhausts_its_window(self, monkeypatch):
        monkeypatch.setattr(
            eclipse, "_lun_eclipse_when_pythonic", lambda *a, **k: (0, (0.0,) * 10)
        )
        with pytest.raises(
            Error,
            match=r"No lunar eclipse visible from longitude 12\.5, latitude 41\.9 "
            r"was found within the search range after JD 2451545\.0\.",
        ):
            lun_eclipse_when_loc(JD, ROME, 2, False)


# ---------------------------------------------------------------------------
# lunar occultations
# ---------------------------------------------------------------------------
class TestLunarOccultation:
    def test_the_moon_cannot_occult_itself(self):
        with pytest.raises(
            Error,
            match=r"A lunar occultation of the Moon itself \(body 1\) is undefined: "
            r"the Moon cannot occult itself\.",
        ):
            lun_occult_when_glob(JD, 1, 2, 0)

    def test_global_search_exhausts_its_window(self, monkeypatch):
        monkeypatch.setattr(eclipse, "_OCCULT_MAX_CONJUNCTIONS", 0)
        with pytest.raises(
            Error,
            match=r"No lunar occultation of 2 was found within the search window "
            r"after JD 2451545\.0\.",
        ):
            lun_occult_when_glob(JD, 2, 2, 0)

    def test_global_backward_search_exhausts_its_window(self, monkeypatch):
        monkeypatch.setattr(eclipse, "_OCCULT_MAX_CONJUNCTIONS", 0)
        with pytest.raises(
            Error,
            match=r"No lunar occultation of 2 was found within the search window "
            r"before JD 2451545\.0\.",
        ):
            lun_occult_when_glob(JD, 2, 2, 0, True)

    def test_local_search_exhausts_its_window(self, monkeypatch):
        monkeypatch.setattr(eclipse, "_OCCULT_MAX_CONJUNCTIONS", 0)
        with pytest.raises(
            Error,
            match=r"No lunar occultation of 2 observable from longitude 12\.5, "
            r"latitude 41\.9 was found after JD 2451545\.0\.",
        ):
            lun_occult_when_loc(JD, 2, ROME)


# ---------------------------------------------------------------------------
# rise, set and transit
# ---------------------------------------------------------------------------
class TestRiseTrans:
    def test_body_no_backend_can_place(self):
        with pytest.raises(
            IllegalBodyError,
            match=r"Body 9999 cannot be placed by the active ephemeris, so its "
            r"rise, set and transit are undefined\.",
        ) as excinfo:
            rise_trans(JD, 9999, libephemeris.CALC_RISE, ROME)
        assert isinstance(excinfo.value, Error)
        assert isinstance(excinfo.value, ValueError)
        assert excinfo.value.body_id == 9999

    def test_body_with_no_direction_in_the_sky(self):
        with pytest.raises(
            Error,
            match=r"Rise, set and transit are undefined for body -1: it has no "
            r"apparent direction in the sky\.",
        ):
            rise_trans(JD, libephemeris.ECL_NUT, libephemeris.CALC_RISE, ROME)


# ---------------------------------------------------------------------------
# Saros and Inex selectors
# ---------------------------------------------------------------------------
class TestSeriesSelectors:
    KIND = r"eclipse_type must be 'solar' or 'lunar', got "

    def test_saros(self):
        with pytest.raises(ValueError, match=self.KIND + "'x'"):
            get_saros_number(JD, "x")

    def test_inex_wrong_type(self):
        with pytest.raises(ValueError, match=self.KIND + "5"):
            get_inex_number(JD, 5)

    def test_inex_unknown_word(self):
        with pytest.raises(ValueError, match=self.KIND + "'x'"):
            get_inex_number(JD, "x")


# ---------------------------------------------------------------------------
# the sampling-step validator of the path and limit tracers
# ---------------------------------------------------------------------------
class TestSamplingStep:
    def test_non_finite_window(self):
        with pytest.raises(
            EphemerisRangeError,
            match=r"tracer: jd_start and jd_end must be finite Julian Days, got "
            r"jd_start=nan, jd_end=1\.0",
        ):
            _eclipse_sampling_step_days(math.nan, 1.0, 1.0, "tracer")

    def test_non_finite_step(self):
        with pytest.raises(
            InputValidationError,
            match=r"tracer: step_minutes must be a finite number of minutes, got nan",
        ):
            _eclipse_sampling_step_days(1.0, 2.0, math.nan, "tracer")

    def test_non_positive_step(self):
        with pytest.raises(
            InputValidationError,
            match=r"tracer: step_minutes must be greater than zero, got 0\.0",
        ):
            _eclipse_sampling_step_days(1.0, 2.0, 0.0, "tracer")

    def test_step_below_the_double_spacing(self):
        with pytest.raises(
            InputValidationError,
            match=r"tracer: step_minutes=1e-12 cannot advance the sampling instant "
            r"near Julian Day 2451546\.0, where consecutive double-precision "
            r"values are .* days apart; use a step larger than .* minutes",
        ):
            _eclipse_sampling_step_days(JD, JD + 1.0, 1e-12, "tracer")


# ---------------------------------------------------------------------------
# planetary occultations, global
# ---------------------------------------------------------------------------
class TestPlanetOccultGlobal:
    def test_no_occulted_body(self):
        with pytest.raises(
            ValueError,
            match=r"Specify the occulted body: pass a non-zero occulted_planet id "
            r"or a starname\.",
        ):
            planet_occult_when_glob(JD, 3, 0, "")

    def test_the_sun_as_occulting_body(self):
        with pytest.raises(
            ValueError,
            match=r"The Sun cannot be the occulting body here; use "
            r"sol_eclipse_when_glob for solar eclipses\.",
        ):
            planet_occult_when_glob(JD, libephemeris.SUN, 3)

    def test_the_moon_as_occulting_body(self):
        with pytest.raises(
            ValueError,
            match=r"The Moon cannot be the occulting body here; use "
            r"lun_occult_when_glob for lunar occultations\.",
        ):
            planet_occult_when_glob(JD, libephemeris.MOON, 3)

    def test_unknown_occulting_planet(self):
        with pytest.raises(
            ValueError, match=r"The occulting planet id 99 does not denote a planet\."
        ):
            planet_occult_when_glob(JD, 99, 3)

    def test_unknown_occulted_planet(self):
        with pytest.raises(
            ValueError, match=r"The occulted planet id 99 does not denote a planet\."
        ):
            planet_occult_when_glob(JD, 3, 99)

    def test_same_planet_twice(self):
        with pytest.raises(
            ValueError,
            match=r"The occulting and the occulted planet must be different bodies\.",
        ):
            planet_occult_when_glob(JD, 3, 3)

    def test_unknown_star(self):
        with pytest.raises(ValueError, match=r"no fixed star matches the search string 'NoSuchStar'"):
            planet_occult_when_glob(JD, 3, 0, "NoSuchStar")

    def test_search_exhausts_the_ephemeris(self, monkeypatch):
        def _outside(*_args, **_kwargs):
            raise ValueError("probe outside the ephemeris range")

        monkeypatch.setattr(eclipse, "calc_ut", _outside)
        with pytest.raises(
            Error,
            match=r"No planetary occultation of planet 4 by planet 3 was found "
            r"within 150 years after JD 2451545\.0\.",
        ):
            planet_occult_when_glob(JD, 3, 4)

    def test_backward_search_exhausts_the_ephemeris(self, monkeypatch):
        def _outside(*_args, **_kwargs):
            raise ValueError("probe outside the ephemeris range")

        monkeypatch.setattr(eclipse, "calc_ut", _outside)
        with pytest.raises(
            Error,
            match=r"No planetary occultation of Regulus by planet 3 was found "
            r"within 150 years before JD 2451545\.0\.",
        ):
            planet_occult_when_glob(JD, 3, 0, "Regulus", 2, -1)


# ---------------------------------------------------------------------------
# planetary occultations, local
# ---------------------------------------------------------------------------
class TestPlanetOccultLocal:
    def test_no_occulted_body(self):
        with pytest.raises(
            ValueError,
            match=r"Specify the occulted body: pass a non-zero occulted_planet id "
            r"or a star_name\.",
        ):
            planet_occult_when_loc(JD, 3, 0, "", 41.9, 12.5)

    def test_the_sun_as_occulting_body(self):
        with pytest.raises(
            ValueError,
            match=r"The Sun cannot be the occulting body here; use "
            r"sol_eclipse_when_loc for solar eclipses\.",
        ):
            planet_occult_when_loc(JD, libephemeris.SUN, 3, "", 41.9, 12.5)

    def test_the_moon_as_occulting_body(self):
        with pytest.raises(
            ValueError,
            match=r"The Moon cannot be the occulting body here; use "
            r"lun_occult_when_loc for lunar occultations\.",
        ):
            planet_occult_when_loc(JD, libephemeris.MOON, 3, "", 41.9, 12.5)

    def test_unknown_occulting_planet(self):
        with pytest.raises(
            ValueError, match=r"The occulting planet id 99 does not denote a planet\."
        ):
            planet_occult_when_loc(JD, 99, 3, "", 41.9, 12.5)

    def test_unknown_occulted_planet(self):
        with pytest.raises(
            ValueError, match=r"The occulted planet id 99 does not denote a planet\."
        ):
            planet_occult_when_loc(JD, 3, 99, "", 41.9, 12.5)

    def test_same_planet_twice(self):
        with pytest.raises(
            ValueError,
            match=r"The occulting and the occulted planet must be different bodies\.",
        ):
            planet_occult_when_loc(JD, 3, 3, "", 41.9, 12.5)

    def test_unknown_star_in_the_visibility_check(self, monkeypatch):
        fake_event = (1, (JD + 1.0, 0.0, JD + 0.9, JD + 1.1) + (0.0,) * 6)
        monkeypatch.setattr(
            eclipse, "planet_occult_when_glob", lambda *a, **k: fake_event
        )
        with pytest.raises(ValueError, match=r"no fixed star matches the search string 'NoSuchStar'"):
            eclipse._planet_occult_when_loc_impl(
                JD, 3, 0, "NoSuchStar", 41.9, 12.5, 0.0, 2, reader=None
            )

    def test_search_finds_no_global_event(self, monkeypatch):
        monkeypatch.setattr(eclipse, "planet_occult_when_glob", _no_eclipse)
        with pytest.raises(
            Error,
            match=r"No planetary occultation of planet 4 by planet 3 visible from "
            r"latitude 41\.9, longitude 12\.5 was found within 150 years after "
            r"JD 2451545\.0\.",
        ):
            planet_occult_when_loc(JD, 3, 4, "", 41.9, 12.5)
