"""Preserve the distinct frozen root/module date-conversion interfaces.

The module entry normalizes date fields; the root entry returns validity, JD and
normalized fields. Wildcard imports in the root stub must not shadow the root
signature with the module signature. These tests change typing, never runtime.
"""
import ast
import inspect
from pathlib import Path

import libephemeris as eph
from libephemeris import time_utils
import pytest


STUB_ROOT = Path(__file__).resolve().parents[1] / "libephemeris"


def stub_function(filename):
    functions = [node for node in ast.parse((STUB_ROOT / filename).read_text()).body
                 if isinstance(node, ast.FunctionDef) and node.name == "date_conversion"]
    assert len(functions) == 1
    return functions[0]


def test_installed_root_signature_and_nested_result():
    function = eph.date_conversion
    signature = inspect.signature(function)
    assert list(signature.parameters) == ["year", "month", "day", "hour", "cal"]
    assert signature.parameters["hour"].default == 12.0
    assert signature.parameters["cal"].default == b"g"
    for code in (b"g", "g"):
        valid, jd, normalized = function(2000, 1, 1, cal=code)
        assert type(valid) is bool and valid
        assert type(jd) is float and jd == 2451545.0
        assert normalized == (2000, 1, 1, 12.0)
        assert tuple(type(value) for value in normalized) == (int, int, int, float)
    assert function(2000, 1, 1) == function(2000, 1, 1, hour=12.0, cal=b"g")
    with pytest.raises(TypeError):
        function(2000, 1, 1, calendar="g")


def test_installed_module_signature_and_flat_date_result():
    function = time_utils.date_conversion
    signature = inspect.signature(function)
    assert list(signature.parameters) == ["year", "month", "day", "hour", "calendar"]
    assert all(parameter.default is inspect.Parameter.empty
               for parameter in signature.parameters.values())
    for code in (b"g", "g"):
        value = function(2000, 1, 1, hour=12.0, calendar=code)
        assert value == (2000, 1, 1, 12.0)
        assert tuple(type(word) for word in value) == (int, int, int, float)
    with pytest.raises(TypeError):
        function(2000, 1, 1)
    with pytest.raises(TypeError):
        function(2000, 1, 1, hour=12.0, cal="g")


def test_installed_stubs_match_the_reviewed_source():
    installed = Path(eph.__file__).resolve().parent
    for filename in ("__init__.pyi", "_core.pyi", "time_utils.pyi"):
        assert (installed / filename).read_bytes() == (STUB_ROOT / filename).read_bytes(), (
            f"installed {filename} is stale or missing the reviewed typing correction")


def test_stubs_preserve_both_distinct_contracts():
    root = stub_function("_core.pyi")
    assert [argument.arg for argument in root.args.args] == [
        "year", "month", "day", "hour", "cal"]
    assert [ast.literal_eval(value) for value in root.args.defaults] == [12.0, b"g"]
    assert ast.unparse(root.args.args[-1].annotation) == "str | bytes"
    assert ast.unparse(root.returns) == "tuple[bool, float, tuple[int, int, int, float]]"

    module = stub_function("time_utils.pyi")
    assert [argument.arg for argument in module.args.args] == [
        "year", "month", "day", "hour", "calendar"]
    assert module.args.defaults == []
    assert ast.unparse(module.args.args[-1].annotation) == "str | bytes"
    assert ast.unparse(module.returns) == "tuple[int, int, int, float]"

    # The final explicit import must restore root typing after all module stars.
    final = ast.parse((STUB_ROOT / "__init__.pyi").read_text()).body[-1]
    assert isinstance(final, ast.ImportFrom) and final.module == "_core"
    assert [(name.name, name.asname) for name in final.names] == [
        ("date_conversion", "date_conversion")]
