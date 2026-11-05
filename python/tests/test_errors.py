"""Public exception classes use native constructors and plain scalar state."""
from data_paths import CORE_LEB

from pathlib import Path
import inspect

import pytest

import libephemeris as eph
from libephemeris import exceptions


@pytest.mark.parametrize("name", [name for name in exceptions.__all__ if isinstance(getattr(exceptions, name), type)])
def test_public_exception_identity(name):
    cls = getattr(exceptions, name)
    assert getattr(eph, name) is cls
    if name in ("LEBCorruptionError", "FictitiousRuntimeDispatch"):
        # These two internal v3 signals derive directly from built-ins.
        base = ValueError if name == "LEBCorruptionError" else KeyError
        assert cls.__bases__ == (base,)
    else:
        assert issubclass(cls, exceptions.Error)
    assert cls.__module__ == "libephemeris.exceptions"


def test_multiple_inheritance():
    illegal = exceptions.IllegalBodyError("unknown", body_id=99999)
    assert isinstance(illegal, exceptions.UnknownBodyError)
    assert isinstance(illegal, ValueError)
    assert illegal.args == ("unknown",)
    assert illegal.body_id == 99999
    sealed = exceptions.NetworkSealedError("disabled", purpose="SPK download")
    assert isinstance(sealed, exceptions.ConfigurationError)
    assert isinstance(sealed, RuntimeError)
    assert sealed.purpose == "SPK download"
    assert sealed.missing_config == "network_policy"


@pytest.mark.parametrize("name,attributes", [
    ("CoordinateError", dict(coordinate_name="latitude", value=91.0, min_value=-90.0, max_value=90.0)),
    ("InvalidBodyError", dict(body_id=0, body_name="Sun", operation="station")),
    ("UnknownBodyError", dict(body_id=99999)),
    ("StarNotFoundError", dict(star_id="missing", search_type="name")),
    ("SPKNotFoundError", dict(filepath="missing.bsp", body_name="Ceres", body_id="1;")),
    ("SPKRequiredError", dict(body_id=17, body_name="Ceres")),
    ("PolarCircleError", dict(latitude=89.0, threshold=66.5, obliquity=23.5, house_system="P")),
    ("EphemerisRangeError", dict(requested_jd=1.0, start_jd=2.0, end_jd=3.0, start_date="a", end_date="b", body_id=1, body_name="Moon", ephemeris_file="base_core.leb2")),
    ("ConvergenceError", dict(algorithm="brent", iterations=100, max_iterations=100, tolerance=1e-8, last_value=0.1)),
    ("SPKEvaluationError", dict(body_id=17, naif_id=2000001, spk_file="ceres.bsp", requested_jd=2451545.0)),
    ("ConfigurationError", dict(missing_config="observer", suggestion="set_topo")),
])
def test_typed_constructor_attributes(name, attributes):
    cls = getattr(exceptions, name)
    error = cls(message="test", **attributes)
    assert error.message == "test"
    assert str(error) == "test"
    assert error.args == ("test",)
    for field, expected in attributes.items():
        assert getattr(error, field) == expected
    defaults = cls("default")
    for field in attributes:
        assert getattr(defaults, field) is None


def test_spk_factory_methods():
    missing = exceptions.SPKNotFoundError.from_filepath("missing.bsp", "Ceres", "1;")
    assert missing.filepath == "missing.bsp"
    assert missing.body_id == "1;"
    required = exceptions.SPKRequiredError.for_body(17, "Ceres", "1;")
    assert required.body_id == 17
    assert required.body_name == "Ceres"


def test_native_failure_maps_public_class():
    eph.set_calc_mode("leb")
    eph.set_precision_tier("base")
    eph.set_leb_file(str(CORE_LEB))
    eph.set_iers_delta_t_enabled(False)
    eph.set_delta_t_userdef(64.0 / 86400.0)
    with pytest.raises(ValueError):
        eph.set_calc_mode("unknown-mode")
    with pytest.raises(exceptions.UnknownBodyError) as caught:
        eph.calc_ut(2451545.0, -999999, 258)
    assert caught.value.body_id == -999999


def test_typed_constructor_rejects_non_scalar_fields():
    with pytest.raises(TypeError):
        exceptions.UnknownBodyError("invalid", body_id=object())
    with pytest.raises(TypeError):
        exceptions.CoordinateError("invalid", value=object())


@pytest.mark.parametrize("name,optional_parameters", [
    ("CoordinateError", ["coordinate_name", "value", "min_value", "max_value"]),
    ("InvalidBodyError", ["body_id", "body_name", "operation"]),
    ("UnknownBodyError", ["body_id"]),
    ("IllegalBodyError", ["body_id"]),
    ("StarNotFoundError", ["star_id", "search_type"]),
    ("SPKNotFoundError", ["filepath", "body_name", "body_id"]),
    ("SPKRequiredError", ["body_id", "body_name"]),
    ("PolarCircleError", ["latitude", "threshold", "obliquity", "house_system"]),
    ("EphemerisRangeError", ["requested_jd", "start_jd", "end_jd", "start_date", "end_date", "body_id", "body_name", "ephemeris_file"]),
    ("ConvergenceError", ["algorithm", "iterations", "max_iterations", "tolerance", "last_value"]),
    ("SPKEvaluationError", ["body_id", "naif_id", "spk_file", "requested_jd"]),
    ("ConfigurationError", ["missing_config", "suggestion"]),
    ("NetworkSealedError", ["purpose"]),
])
def test_constructor_parameter_names_and_defaults(name, optional_parameters):
    signature = inspect.signature(getattr(exceptions, name))
    assert list(signature.parameters) == ["message", *optional_parameters]
    assert signature.parameters["message"].default is inspect.Parameter.empty
    for name in optional_parameters:
        assert signature.parameters[name].default is None


@pytest.mark.parametrize("name, value, field, low, high", [
    ("validate_latitude", 91.0, "latitude", -90.0, 90.0),
    ("validate_longitude", 361.0, "longitude", -180.0, 360.0),
    ("validate_latitude", float("nan"), "latitude", -90.0, 90.0),
    ("validate_longitude", float("inf"), "longitude", -180.0, 360.0),
])
def test_native_coordinate_validation_preserves_documented_attributes(name, value, field, low, high):
    with pytest.raises(exceptions.CoordinateError) as caught:
        getattr(exceptions, name)(value, func_name="observer")
    assert caught.value.coordinate_name == field
    assert caught.value.min_value == low
    assert caught.value.max_value == high


def test_coordinate_validation_accepts_both_longitude_conventions():
    assert exceptions.validate_latitude(-90.0) is None
    assert exceptions.validate_longitude(200.0) is None
    assert exceptions.validate_coordinates(90.0, 360.0) is None
    with pytest.raises(exceptions.CoordinateError) as caught:
        exceptions.validate_coordinates(91.0, 361.0)
    assert caught.value.coordinate_name == "latitude"
