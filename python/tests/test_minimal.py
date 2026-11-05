"""Keep the installed package's source boundary minimal."""

import importlib.util
from pathlib import Path


def test_minimal_package_and_bindings():
    root = Path(__file__).resolve().parents[2]
    spec = importlib.util.spec_from_file_location("check_minimal", root / "tools/release/check_minimal.py")
    checker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checker)
    assert checker.check(root, "P4") == []


def test_native_modules_cannot_hide_foreign_scalar_arithmetic(tmp_path):
    root = Path(__file__).resolve().parents[2]
    spec = importlib.util.spec_from_file_location("check_minimal", root / "tools/release/check_minimal.py")
    checker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checker)
    package = tmp_path / "python/libephemeris"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text("")
    core = tmp_path / "src/photometry.rs"
    core.parent.mkdir()
    core.write_text("pub type PositionalArithmetic = dyn Fn(f64) -> f64;")
    assert any("scalar operand adapter" in error for error in checker.check(tmp_path))
    core.write_text("pub fn cosine(value: f64) -> f64 { value.cos() }")
    assert checker.check(tmp_path) == []
