"""Regression tests for forbidden binding logic and phased core removal."""

import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location("check_minimal", Path(__file__).parents[1] / "release/check_minimal.py")
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


class MinimalPackageTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.package = self.root / "python/libephemeris"
        self.package.mkdir(parents=True)
        self.initializer = self.package / "__init__.py"
        self.initializer.write_text('"""Public imports."""\nfrom ._core import calc, __version__\n__all__ = ["calc"]\n')
        self.bindings = self.root / "python/src"
        self.bindings.mkdir()

    def test_imports_pass_but_validation_and_computed_exports_fail(self):
        self.assertEqual(checker.check(self.root), [])
        for code in ('def calc(x): return x + 1\n', '__all__ = list(globals())\n', 'if True: import math\n'):
            with self.subTest(code=code):
                self.initializer.write_text(code)
                self.assertTrue(checker.check(self.root))

    def test_embedded_execution_and_python_arithmetic_fail(self):
        for code in ('PyModule::from_code(py, source, file, name);', 'include_str!("runtime.py")',
                     'PyModule::import(py, "math")', 'obj.call_method1("__add__", (rhs,))',
                     'ffi::PyNumber_Multiply(a, b)'):
            with self.subTest(code=code):
                (self.bindings / "lib.rs").write_text(code)
                self.assertTrue(checker.check(self.root))

    def test_core_cleanup_is_required_at_p4(self):
        (self.root / "src/harness").mkdir(parents=True)
        self.assertEqual(checker.check(self.root, "P0"), [])
        self.assertTrue(checker.check(self.root, "P4"))

    def test_line_budget_is_enforced(self):
        (self.bindings / "lib.rs").write_text("// binding\n" * 20_001)
        self.assertTrue(checker.check(self.root))

    def test_remote_pipelines_are_forbidden_even_in_reset_phase(self):
        for name in (".github/workflows/test.yml", ".gitlab-ci.yml",
                     ".circleci/config.yml", "Jenkinsfile", "azure-pipelines.yaml"):
            with self.subTest(name=name):
                path = self.root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("# remote pipeline\n")
                for phase in ("P0", "P4"):
                    self.assertTrue(any("CI/CD is forbidden" in error
                                        for error in checker.check(self.root, phase)))
                path.unlink()

    def test_local_checks_and_workflow_documentation_are_allowed(self):
        (self.root / "tools").mkdir()
        (self.root / "tools/local-checks.sh").write_text("cargo test\n")
        docs = self.root / "docs"
        docs.mkdir()
        (docs / "workflows.md").write_text("Run the checks locally.\n")
        self.assertEqual(checker.check(self.root), [])

    def test_retired_nbody_paths_fail_without_reading_upstream_contents(self):
        for name in ("src/bodies/nbody/ffi.rs", "native/rebound", "native/assist"):
            with self.subTest(name=name):
                path = self.root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                if path.suffix == ".rs":
                    path.write_text("// Restored adapter.\n")
                else:
                    path.mkdir()
                    # A path check must not decode restored foreign contents.
                    (path / "unread-source.c").write_bytes(b"\xff")
                self.assertEqual(checker.check(self.root, "P0"), [])
                self.assertTrue(any("retired N-body implementation path" in error
                                    for error in checker.check(self.root)))
                if path.is_dir():
                    (path / "unread-source.c").unlink()
                    path.rmdir()
                else:
                    path.unlink()

    def test_legacy_c_build_and_link_routes_fail(self):
        build = self.root / "build.rs"
        for code in (
            'fn main() { cc::Build::new().include("native/rebound/src"); }',
            'fn main() { cc::Build::new().include(r#"native/assist/src"#); }',
            'fn main() { cc::Build::new().include(r"native\\rebound\\src"); }',
            'fn main() { cc::Build::new().compile("ephemeris_assist"); }',
            'fn main() { cc::Build::new().compile("rebound"); }',
            'fn main() { println!("cargo:rustc-link-lib=static=assist"); }',
        ):
            with self.subTest(code=code):
                build.write_text(code)
                self.assertTrue(any("retired N-body C compilation/link" in error
                                    for error in checker.check(self.root)))
        build.unlink()
        (self.root / "Cargo.toml").write_text('[package]\nbuild = "owned-build.rs"\n')
        (self.root / "owned-build.rs").write_text(
            'fn main() { cc::Build::new().compile("ephemeris_rebound"); }')
        self.assertTrue(any("retired N-body C compilation/link" in error
                            for error in checker.check(self.root)))

    def test_delivery_manifests_cannot_require_removed_sources_or_licenses(self):
        for name, text in (
            ("pyproject.toml", '[project]\nlicense-files = ["native/rebound/LICENSE"]\n'),
            ("pyproject.toml", '[tool.maturin]\ninclude = [{path="native/assist/**", format="sdist"}]\n'),
            ("Cargo.toml", '[package]\ninclude = ["native/rebound/**"]\n'),
            ("Cargo.toml", '[dependencies]\nhelper = {path="../native/assist"}\n'),
            ("python/Cargo.toml", '[package]\nbuild = "../native/rebound/abi.c"\n'),
        ):
            with self.subTest(name=name, text=text):
                path = self.root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(text)
                self.assertTrue(any("retired N-body source/license delivery" in error
                                    for error in checker.check(self.root)))
                path.unlink()

    def test_normal_runtime_build_and_tool_dependencies_cannot_restore_engines(self):
        for name, text in (
            ("pyproject.toml", '[project]\ndependencies = ["rebound>=4"]\n'),
            ("pyproject.toml", '[build-system]\nrequires = ["assist==1.2.3"]\n'),
            ("pyproject.toml", '[dependency-groups]\ntest = ["rebound"]\n'),
            ("tools/leb_generation/pyproject.toml", '[project]\ndependencies = ["assist"]\n'),
            ("Cargo.toml", '[dependencies]\nlegacy = {package="rebound-sys", version="1", optional=true}\n'),
            ("python/Cargo.toml", '[build-dependencies]\nassist-sys = "1"\n'),
        ):
            with self.subTest(name=name, text=text):
                path = self.root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(text)
                self.assertTrue(any("retired N-body delivery dependency" in error
                                    for error in checker.check(self.root)))
                path.unlink()

    def test_runtime_and_offline_python_helper_imports_fail(self):
        offline = self.root / "tools/leb_generation/source.py"
        offline.parent.mkdir(parents=True)
        runtime = self.package / "rebound_integration.py"
        for path in (runtime, offline):
            for code in (
                "import rebound\n", "import assist as legacy\n",
                "from rebound.simulation import Simulation\n",
                'import importlib\nimportlib.import_module("assist")\n',
                'from importlib import import_module as load\nload("rebound")\n',
                '__import__("rebound")\n',
                'from . import assist\n',
                'from .rebound import Simulation\n',
                "def source():\n    import assist\n",
            ):
                with self.subTest(path=path, code=code):
                    path.write_text(code)
                    self.assertTrue(any("retired N-body Python helper import" in error
                                        for error in checker.check(self.root)))
            path.unlink()

    def test_bindings_and_core_cannot_import_legacy_python_helpers(self):
        core = self.root / "src/nbody.rs"
        core.parent.mkdir()
        for path in (self.bindings / "nbody.rs", core):
            for code in ('PyModule::import(py, "rebound");',
                         'PyModule::import(py, "assist.ephem");',
                         'PyModule::import(py, intern!(py, "assist"));',
                         'py.import("assist");'):
                with self.subTest(path=path, code=code):
                    path.write_text(code)
                    self.assertTrue(any("retired N-body Python helper import" in error
                                        for error in checker.check(self.root)))
            path.write_text('pub fn get_rebound_version() {}\n'
                            'pub fn propagate_orbit_assist() {}\n'
                            'PyModule::import(py, "libephemeris")?.call_method1("mode", ("assist",));\n')
            self.assertEqual(checker.check(self.root), [])

    def test_permissive_c_compatibility_names_and_isolated_history_stay_allowed(self):
        (self.root / "build.rs").write_text(
            '// Retired include("native/rebound/src"); compile("ephemeris_assist");\n'
            '/* Historical "native/assist/src" path. */\n'
            'fn main() { cc::Build::new().file("native/zstd/decode.c").compile("zstd"); }\n')
        (self.root / "Cargo.toml").write_text(
            '[package]\nlicense = "AGPL-3.0"\nexclude = ["native/rebound/**", "native/assist/**"]\n'
            '[dependencies]\nzstd = "0.14"\nsha2 = "0.11"\n[build-dependencies]\ncc = "1"\n')
        (self.root / "pyproject.toml").write_text(
            '[project]\ndependencies = ["numpy>=1.26"]\nlicense-files = ["COPYING", "licenses/**/*"]\n'
            '[tool.maturin]\nexclude = ["native/rebound/**", "native/assist/**"]\n'
            '[dependency-groups]\nassist = ["numpy"]\ntest = [{include-group="assist"}]\n')
        (self.package / "rebound_integration.py").write_text(
            '"""Historical REBOUND/ASSIST public names."""\n'
            'from ._core import get_rebound_version, get_assist_version, ReboundIntegrator\n')
        offline = self.root / "tools/leb_generation/source.py"
        offline.parent.mkdir(parents=True)
        offline.write_text(
            'import importlib\nfrom libephemeris.rebound_integration import AssistEphemConfig\n'
            'def source():\n    return importlib.import_module("libephemeris.rebound_integration")\n')
        reference = self.root / "validation/isolated-reference"
        reference.mkdir(parents=True)
        (reference / "reference.py").write_text("import rebound\nimport assist\n")
        (reference / "pyproject.toml").write_text('[project]\ndependencies = ["rebound", "assist"]\n')
        history = self.root / "docs/archive"
        history.mkdir(parents=True)
        (history / "old.md").write_text('native/rebound/src; import assist; ffi.rs; build.rs\n')
        self.assertEqual(checker.check(self.root), [])


if __name__ == "__main__":
    unittest.main()
