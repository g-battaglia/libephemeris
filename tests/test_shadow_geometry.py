"""
Tests for the ShadowGeometry record and its invariants.

The record carries the two shadow cones, the shadow axis and the shadowed
body from the eclipse module's geometry to its consumers. What is checked
here is the contract itself -- immutability, keyword construction, the
absence of any positional access, the derived half-angle sines, the
validating method -- and then the invariants on records the solar/occultation
producer actually builds: the penumbral section is a real width, the two cone
cosines are ordered, the axis offset is a distance, the surface pair travels
together, the miss flag agrees with the classification returned beside it, and
the sign of the inner cone's diameter agrees with the total/annular bit.
"""

from __future__ import annotations

import dataclasses
import math

import pytest

from libephemeris import (
    ECL_ANNULAR,
    ECL_CENTRAL,
    ECL_NONCENTRAL,
    ECL_PARTIAL,
    ECL_TOTAL,
    FLG_SWIEPH,
    julday,
)
from libephemeris.eclipse import _eclipse_where_core
from libephemeris.shadow_geometry import ShadowGeometry


# A well-formed record, of the shape the solar producer builds.
def _sample(**overrides) -> ShadowGeometry:
    fields = dict(
        axis_offset_km=5698.34,
        umbral_plane_diameter_km=-168.13,
        penumbral_plane_diameter_km=6801.11,
        cos_umbral_half_angle=0.9999895,
        cos_penumbral_half_angle=0.9999894,
        shadowed_radius_km=6378.14,
        shadow_misses_body=False,
        umbral_surface_diameter_km=-194.09,
        penumbral_surface_diameter_km=6827.06,
    )
    fields.update(overrides)
    offset = fields.pop("axis_offset_km")
    return ShadowGeometry(offset, **fields)


#: Instants of recorded solar events, one per class, plus one date with no
#: eclipse anywhere. Julian Days in UT.
_SOLAR_INSTANTS = (
    2435633.388979057,  # central, total
    2396800.770489158,  # central, annular
    2425386.0582784833,  # non-central, total
    2496466.4159777206,  # non-central, annular
    2400255.4144638283,  # non-central, partial
    2427807.7328914925,  # no eclipse anywhere on Earth
)


class TestRecordContract:
    """The shape of the type, independent of any eclipse."""

    def test_is_frozen(self):
        """A record cannot be edited after construction."""
        geometry = _sample()
        with pytest.raises(dataclasses.FrozenInstanceError):
            geometry.axis_offset_km = 1.0

    def test_equal_by_value(self):
        """Two records with the same fields are equal."""
        assert _sample() == _sample()
        assert _sample() != _sample(axis_offset_km=1.0)

    def test_only_the_first_field_is_positional(self):
        """Everything past the axis offset must be named."""
        with pytest.raises(TypeError):
            ShadowGeometry(5698.34, -168.13)  # type: ignore[call-arg]

    def test_is_not_a_sequence(self):
        """The record has no length, no index and no unpacking."""
        geometry = _sample()
        with pytest.raises(TypeError):
            len(geometry)  # type: ignore[arg-type]
        with pytest.raises(TypeError):
            geometry[0]  # type: ignore[index]
        with pytest.raises(TypeError):
            _first, _second = geometry  # type: ignore[misc]

    def test_sines_are_derived_from_the_cosines(self):
        """The sines are computed, never stored beside their cosine."""
        geometry = _sample()
        assert geometry.sin_umbral_half_angle == pytest.approx(
            math.sqrt(1.0 - geometry.cos_umbral_half_angle**2), abs=1e-15
        )
        assert geometry.sin_penumbral_half_angle == pytest.approx(
            math.sqrt(1.0 - geometry.cos_penumbral_half_angle**2), abs=1e-15
        )
        stored = {field.name for field in dataclasses.fields(geometry)}
        assert "sin_umbral_half_angle" not in stored
        assert "sin_penumbral_half_angle" not in stored

    def test_sines_are_non_negative_for_a_point_source(self):
        """A point-like source makes the two half angles equal, not negative."""
        geometry = _sample(
            cos_umbral_half_angle=0.99999,
            cos_penumbral_half_angle=0.99999,
        )
        assert geometry.sin_umbral_half_angle >= 0.0
        assert geometry.sin_umbral_half_angle == geometry.sin_penumbral_half_angle

    def test_the_surface_pair_may_be_left_unset(self):
        """A producer that resolves no surface point leaves both unset."""
        geometry = _sample(
            umbral_surface_diameter_km=None, penumbral_surface_diameter_km=None
        )
        assert geometry.umbral_surface_diameter_km is None
        assert geometry.penumbral_surface_diameter_km is None
        geometry.validate()


class TestValidate:
    """The invariants the validating method states."""

    def test_a_well_formed_record_passes(self):
        _sample().validate()

    def test_the_penumbral_section_must_be_a_width(self):
        with pytest.raises(ValueError, match="positive width"):
            _sample(penumbral_plane_diameter_km=-1.0).validate()

    def test_the_cosines_must_be_ordered(self):
        with pytest.raises(ValueError, match="opens at least as wide"):
            _sample(
                cos_umbral_half_angle=0.9,
                cos_penumbral_half_angle=0.99,
            ).validate()

    def test_the_axis_offset_is_a_distance(self):
        with pytest.raises(ValueError, match="is a distance"):
            _sample(axis_offset_km=-1.0).validate()

    def test_the_shadowed_body_has_a_positive_radius(self):
        with pytest.raises(ValueError, match="positive radius"):
            _sample(shadowed_radius_km=0.0).validate()

    def test_half_a_surface_pair_is_malformed(self):
        with pytest.raises(ValueError, match="together or not at all"):
            _sample(penumbral_surface_diameter_km=None).validate()


class TestProducedRecords:
    """Invariants on the records the solar/occultation producer builds."""

    @pytest.mark.parametrize("tjd_ut", _SOLAR_INSTANTS)
    def test_every_invariant_holds(self, tjd_ut):
        _retflag, _lon, _lat, geometry = _eclipse_where_core(tjd_ut, FLG_SWIEPH)
        geometry.validate()

    @pytest.mark.parametrize("tjd_ut", _SOLAR_INSTANTS)
    def test_fields_are_native_python_types(self, tjd_ut):
        """No numpy scalar and no integer may reach the record."""
        _retflag, _lon, _lat, geometry = _eclipse_where_core(tjd_ut, FLG_SWIEPH)
        for field in dataclasses.fields(geometry):
            value = getattr(geometry, field.name)
            if field.name == "shadow_misses_body":
                assert type(value) is bool
            else:
                assert type(value) is float
                assert math.isfinite(value)

    @pytest.mark.parametrize("tjd_ut", _SOLAR_INSTANTS)
    def test_the_miss_flag_agrees_with_the_classification(self, tjd_ut):
        """I7: the flag is set exactly when no shadow reaches the Earth."""
        retflag, _lon, _lat, geometry = _eclipse_where_core(tjd_ut, FLG_SWIEPH)
        assert geometry.shadow_misses_body == (retflag == 0)

    @pytest.mark.parametrize("tjd_ut", _SOLAR_INSTANTS)
    def test_the_sign_of_the_inner_cone_carries_the_class(self, tjd_ut):
        """I8: total implies a negative surface diameter, annular a positive."""
        retflag, _lon, _lat, geometry = _eclipse_where_core(tjd_ut, FLG_SWIEPH)
        diameter = geometry.umbral_surface_diameter_km
        assert diameter is not None
        if retflag & ECL_TOTAL:
            assert diameter < 0.0
        if retflag & ECL_ANNULAR:
            assert diameter > 0.0

    @pytest.mark.parametrize("tjd_ut", _SOLAR_INSTANTS)
    def test_the_shadowed_body_is_the_earth(self, tjd_ut):
        """The one radius the consumers read is the Earth's equatorial one."""
        _retflag, _lon, _lat, geometry = _eclipse_where_core(tjd_ut, FLG_SWIEPH)
        assert geometry.shadowed_radius_km == pytest.approx(6378.140, abs=1e-9)

    def test_a_point_source_opens_both_cones_alike(self):
        """A star has no disc: the two cones coincide (C-E01 spec, I4)."""
        tjd_ut = 2430031.9520767587  # a recorded occultation of Aldebaran
        retflag, _lon, _lat, geometry = _eclipse_where_core(
            tjd_ut, FLG_SWIEPH, "Aldebaran"
        )
        assert geometry.cos_umbral_half_angle == geometry.cos_penumbral_half_angle
        assert not retflag & (ECL_ANNULAR | ECL_PARTIAL)
        # The core cone is the Moon's own shadow, so its section at the ground
        # is minus the lunar diameter, 2 x 1738.15 km.
        assert geometry.umbral_surface_diameter_km == pytest.approx(-3476.30, abs=0.01)

    def test_the_class_is_one_of_the_documented_combinations(self):
        """Central and non-central never appear together (spec section 5)."""
        for tjd_ut in _SOLAR_INSTANTS:
            retflag, _lon, _lat, _geometry = _eclipse_where_core(tjd_ut, FLG_SWIEPH)
            central = bool(retflag & ECL_CENTRAL)
            noncentral = bool(retflag & ECL_NONCENTRAL)
            assert not (central and noncentral)
            assert (central or noncentral) == (retflag != 0)
            assert not (retflag & ECL_TOTAL and retflag & ECL_ANNULAR)
            if retflag & ECL_PARTIAL:
                assert noncentral
                assert not retflag & (ECL_TOTAL | ECL_ANNULAR)


class TestGroundPoint:
    """The coordinates the record travels beside."""

    @pytest.mark.parametrize("tjd_ut", _SOLAR_INSTANTS)
    def test_longitude_is_east_positive_in_the_half_open_turn(self, tjd_ut):
        _retflag, lon, lat, _geometry = _eclipse_where_core(tjd_ut, FLG_SWIEPH)
        assert -180.0 < lon <= 180.0
        assert -90.0 <= lat <= 90.0

    def test_a_central_eclipse_puts_the_point_on_the_axis(self):
        """The reported point of a central eclipse is where the axis lands.

        Checked through its observable signature: the Sun stands high above
        the horizon there, which is what a genuine sub-shadow point means
        (eclipse_where_core spec, section 6.2). The 1955 June 20 total
        eclipse is one of the recorded central events.
        """
        from libephemeris import sol_eclipse_where

        retflag, geopos, attr = sol_eclipse_where(2435633.388979057, FLG_SWIEPH)
        assert retflag & ECL_CENTRAL
        assert attr[5] > 1.0  # true altitude of the Sun at the point
        assert geopos[0] == pytest.approx(-140.7461, abs=0.01)
        assert geopos[1] == pytest.approx(-40.7513, abs=0.01)


def test_the_ground_point_matches_the_published_canon():
    """The 2017 August 21 greatest eclipse against the Five Millennium Canon.

    The ground point is the exact intersection of the shadow axis with the
    reference ellipsoid, reported as a geodetic latitude (Chauvenet;
    Explanatory Supplement, 3rd ed., ch. 11). NASA's Five Millennium Canon of
    Solar Eclipses puts greatest eclipse at 36 deg 58.0' N, 87 deg 40.3' W;
    the construction lands within a kilometre of it.
    """
    from libephemeris import sol_eclipse_when_glob, sol_eclipse_where

    _rf, tret = sol_eclipse_when_glob(julday(2017, 8, 21, 0.0), FLG_SWIEPH, 0, False)
    _retflag, geopos, _attr = sol_eclipse_where(tret[0], FLG_SWIEPH)
    assert geopos[1] == pytest.approx(36.0 + 58.0 / 60.0, abs=0.02)
    assert geopos[0] == pytest.approx(-(87.0 + 40.3 / 60.0), abs=0.02)
