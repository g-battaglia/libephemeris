"""Behavioural smoke tests for libephemeris.contrib helpers.

These complement the surface-only checks in ``tests/test_api_compat/`` by
exercising a handful of known-value inputs. They're intentionally small
and self-contained: their job is to catch obviously-wrong implementations
(off-by-one, swapped fields, broken carry, ...) rather than full
numerical verification.
"""

from __future__ import annotations


import pytest

import libephemeris.contrib as c


class TestSexagesimalConversions:
    """degsplit / geolat2c / geolon2c overflow handling."""

    @pytest.mark.parametrize(
        "deg,expected",
        [
            (0.0, (0, 0, 0, 0)),
            (15.5, (15, 0, 30, 0)),
            (30.0, (0, 1, 0, 0)),
            # Rounding overflow must carry through every field:
            (29.9999, (0, 1, 0, 0)),
            (359.99999, (0, 0, 0, 0)),  # wraps to start of zodiac
        ],
    )
    def test_degsplit(self, deg, expected):
        assert c.degsplit(deg) == expected

    @pytest.mark.parametrize(
        "lat,expected",
        [
            (0.0, "0N00"),
            (41.9, "41N54"),
            (-33.999, "34S00"),  # carry over the degree boundary
            (12.999, "13N00"),
        ],
    )
    def test_geolat2c(self, lat, expected):
        assert c.geolat2c(lat) == expected

    @pytest.mark.parametrize(
        "lon,expected",
        [
            (12.5, "12E30"),
            (-179.999, "180W00"),
            (0.0, "0E00"),
        ],
    )
    def test_geolon2c(self, lon, expected):
        assert c.geolon2c(lon) == expected


class TestZodiacalConversions:
    """long2rasi / long2nakshatra / long2navamsa / signtostr / round-trip."""

    @pytest.mark.parametrize(
        "deg,rasi",
        [
            (0.0, 0),  # 0° Aries
            (29.9, 0),  # still Aries
            (30.0, 1),  # exactly Taurus
            (35.0, 1),  # Taurus
            (180.0, 6),  # Libra
            (359.0, 11),  # Pisces
        ],
    )
    def test_long2rasi(self, deg, rasi):
        assert c.long2rasi(deg) == rasi

    @pytest.mark.parametrize(
        "deg,nak,pada",
        [
            (0.0, 0, 0),  # 0° → Aswini, pada 1 (0-indexed 0)
            (13.5, 1, 0),  # 13.5° → Bharani
            (360 / 27 - 0.01, 0, 3),  # last pada of Aswini
        ],
    )
    def test_long2nakshatra(self, deg, nak, pada):
        assert c.long2nakshatra(deg) == (nak, pada)

    @pytest.mark.parametrize(
        "deg,navamsa",
        [
            (0.0, 0),  # Aries 0° → navamsa Aries
            (15.0, 4),  # Aries mid → navamsa Leo
            (30.0, 9),  # Taurus 0° → navamsa Capricorn
        ],
    )
    def test_long2navamsa(self, deg, navamsa):
        assert c.long2navamsa(deg) == navamsa

    def test_signtostr_roundtrip(self):
        """signtostr(long2rasi(deg)) is the expected sign string."""
        assert c.signtostr(c.long2rasi(45.0)) == "Taurus"
        assert c.signtostr(c.long2rasi(15.0)) == "Aries"
        assert c.signtostr(c.long2rasi(180.0)) == "Libra"

    def test_nakshatra_name_lookup(self):
        assert c.get_nakshatra_name(0) == "Aswini"
        assert c.get_nakshatra_name(26) == "Revathi"
        # Modular: 27 → 0 again
        assert c.get_nakshatra_name(27) == "Aswini"


class TestVedicDignities:
    """lord / ochchabala / naisargika_relation / tatkalika_relation."""

    @pytest.mark.parametrize(
        "rasi,lord_id",
        [
            (c.ARIES, c.KUJA),  # Mars rules Aries
            (c.TAURUS, c.SUKRA),  # Venus rules Taurus
            (c.LEO, c.SURYA),  # Sun rules Leo
            (c.SCORPIO, c.KUJA),  # Mars also rules Scorpio
            (c.CAPRICORN, c.SANI),
            (c.PISCES, c.GURU),
        ],
    )
    def test_lord(self, rasi, lord_id):
        assert c.lord(rasi) == lord_id

    def test_ochchabala_extremes(self):
        """At exaltation degree, ochchabala = 60; at debilitation, ~0."""
        # Sun exalts at Aries 10° (longitude 10.0)
        ex = c.ochchabala(c.SURYA, 10.0)
        assert abs(ex - 60.0) < 1e-9
        # Sun debilitates 180° away (Libra 10° = 190°)
        deb = c.ochchabala(c.SURYA, 190.0)
        assert abs(deb) < 1e-9

    def test_naisargika_symmetric_anchors(self):
        # Sun and Moon are mutual friends in the natural relation table
        assert c.naisargika_relation(c.SURYA, c.CHANDRA) == 2
        # Sun and Venus are natural enemies
        assert c.naisargika_relation(c.SURYA, c.SUKRA) == 0

    def test_tatkalika_same_house(self):
        # Same house → neutral
        assert c.tatkalika_relation(5, 5) == 1


class TestAspectHelpers:
    """match_aspect family + antiscion."""

    def test_match_aspect_within_orb(self):
        # Sun at 0°, Mars at 90°: square aspect, exact
        is_match, delta = c.match_aspect(0, 0.0, c.SQUARE, 4, 90.0, 1.0)
        assert is_match == 1
        assert delta < 1e-9

    def test_match_aspect_outside_orb(self):
        # 91° separation vs 90° square with orb=0.5 → no match
        is_match, delta = c.match_aspect(0, 0.0, c.SQUARE, 4, 91.0, 0.5)
        assert is_match == 0
        assert abs(delta - 1.0) < 1e-9

    def test_match_aspect_wraps_through_180(self):
        # 359° vs 1° = 2° separation (not 358°)
        is_match, delta = c.match_aspect2(359.0, 1.0, c.CONJUNCTION, 5.0)
        assert is_match == 1

    def test_antiscion(self):
        # Antiscion is reflected across Cancer/Capricorn axis: 180° - lon
        ant, contra = c.antiscion(15.0)
        # ant[0] = antiscion longitude; contra[0] = contra-antiscion
        assert abs(ant[0] - 165.0) < 1e-9
        assert abs(contra[0] - 345.0) < 1e-9


class TestTimeJDHelpers:
    """jd2isostr / jduration / years_diff / t2i / dt2i."""

    def test_jd2isostr_j2000(self):
        # J2000.0 = 2000-01-01 12:00 UT
        assert c.jd2isostr(2451545.0) == "2000-01-01T12:00:00"

    def test_jd2isostr_handles_rollover(self):
        # 23:59:59.9 must NOT become :60 — it should roll to the next day
        jd = c.julday(2024, 1, 1, 23 + 59 / 60 + 59.9 / 3600)
        assert c.jd2isostr(jd) == "2024-01-02T00:00:00"

    def test_jduration(self):
        assert c.jduration(2451545.0, 2451550.5) == 5.5

    def test_t2i_parses_time(self):
        assert c.t2i("12:34:56") == (12, 34, 56)
        assert c.t2i("06:30") == (6, 30, 0)

    def test_dt2i_parses_iso(self):
        assert c.dt2i("2024-04-08 18:42:00") == (2024, 4, 8, 18, 42, 0)


class TestRasiDistance:
    """rasinorm / rasi_dif / rasi_dif2."""

    def test_rasinorm_mod12(self):
        assert c.rasinorm(0) == 0
        assert c.rasinorm(12) == 0
        assert c.rasinorm(15) == 3
        assert c.rasinorm(-1) == 11

    def test_rasi_dif_forward(self):
        # Forward distance from Aries to Libra is 6 rasis
        assert c.rasi_dif(c.ARIES, c.LIBRA) == 6

    def test_rasi_dif2_signed(self):
        # Signed distance in (-6, 6] — opposite signs go the short way
        assert c.rasi_dif2(c.ARIES, c.LIBRA) == 6  # ± exactly 6
        assert c.rasi_dif2(c.ARIES, c.GEMINI) == 2  # forward 2
        assert c.rasi_dif2(c.GEMINI, c.ARIES) == -2  # backward 2


class TestNotImplementedStubs:
    """Database/atlas/timezone helpers must raise NotImplementedError."""

    @pytest.mark.parametrize(
        "name",
        [
            "atlas_connect",
            "atlas_close",
            "atlas_search",
            "atlas_countries_list",
            "db_connect",
            "db_close",
            "tzabbr_find",
            "tzabbr_list",
        ],
    )
    def test_stub_raises(self, name):
        fn = getattr(c, name)
        with pytest.raises(NotImplementedError):
            fn()
