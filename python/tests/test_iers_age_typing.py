"""The public IERS cache-age constant has a concrete stub type."""

import ast
from pathlib import Path

from libephemeris import iers_data


STUB_ROOT = Path(__file__).resolve().parents[1] / "libephemeris"


def test_public_cache_age_is_a_native_integer():
    assert type(iers_data.DEFAULT_MAX_AGE_DAYS) is int
    assert iers_data.DEFAULT_MAX_AGE_DAYS == 30


def test_public_cache_age_has_a_concrete_stub_declaration():
    public_stub = ast.parse((STUB_ROOT / "iers_data.pyi").read_text())
    declarations = [
        node for node in public_stub.body
        if isinstance(node, ast.AnnAssign)
        and isinstance(node.target, ast.Name)
        and node.target.id == "DEFAULT_MAX_AGE_DAYS"
    ]
    assert len(declarations) == 1
    assert ast.unparse(declarations[0].annotation) == "int"
