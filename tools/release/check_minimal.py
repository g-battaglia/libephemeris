#!/usr/bin/env python3
"""Check the completed v4 package and core boundaries.

Package files may contain imports, docstrings, literal __all__ declarations and
__version__. Bindings convert values and call Rust; they must not execute Python
source or implement an arithmetic engine over Python objects. Core removal checks
also reject retired scalar-operand and foreign-exception adapters, even when
they were left in otherwise native modules. --phase P0 is for the reset only.
Normal delivery also excludes retired N-body implementations and helper imports.
Scientific comparison tooling is maintained in a separate repository.
"""

from __future__ import annotations

import argparse
import ast
from collections.abc import Iterator
from pathlib import Path
import re
import tomllib

ROOT = Path(__file__).resolve().parents[2]
CI_PATHS = (
    ".github/workflows/*.yml", ".github/workflows/*.yaml",
    ".gitlab-ci.yml", ".circleci/config.yml", ".circleci/config.yaml",
    ".buildkite/pipeline.yml", ".buildkite/pipeline.yaml", "Jenkinsfile*",
    "azure-pipelines.yml", "azure-pipelines.yaml", "bitbucket-pipelines.yml",
    ".travis.yml", "appveyor.yml", ".drone.yml", ".woodpecker.yml",
)
CORE_PATHS = (
    "src/harness", "src/foundation/python_exception*.rs",
    "src/positions/resource_*.rs", "src/positions/operand_spk.rs",
    "src/phenomena/operand_*.rs", "src/time/operand_*.rs",
    "src/providers/leb/source_*.rs", "bindings",
)
N_BODY_RUST_IMPORT = (
    r'\b(?:PyModule::import(?:_bound)?\s*\(\s*[^,;]*,\s*|py\.import\s*\(\s*)'
    r'(?:intern!\(\s*[^,;]*,\s*)?(?:r\#{0,255})?"'
    r'(?:rebound|assist|librebound|libassist)(?:\.[A-Za-z0-9_.]*)?"',
    "retired N-body Python helper import",
)
BINDING_PATTERNS = (
    (r"\bPyModule::from_code\b", "executed Python source"),
    (r"include_(?:str|bytes)!\s*\([^)]*\.py[\"']", "embedded Python file"),
    (r"\bPyModule::import\s*\([^;]*?[\"']math[\"']", "Python math import"),
    (r"[\"']__(?:add|sub|mul|truediv|floordiv|mod|pow|matmul|neg|pos|abs)__[\"']",
     "Python arithmetic protocol"),
    (r"\bPyNumber_(?:Add|Subtract|Multiply|TrueDivide|FloorDivide|Remainder|Power|MatrixMultiply)\b",
     "Python arithmetic API"),
    N_BODY_RUST_IMPORT,
)
CORE_PATTERNS = (
    (r"\b(?:PositionalPrimitive|PositionalArithmetic|OrbitPrimitive|OrbitPrimitives|StoredStarOperand)\b",
     "retired scalar operand adapter"),
    (r"\b(?:TraceExceptionCatch|RecordObject|capture_events_with_fallible_sink_and_catch|catches_exception)\b",
     "retired foreign object/exception transport"),
    (r"\b(?:StarPositionCompletion|StarSourceScope|StarSourceScopeKind|with_star_source_scope)\b",
     "retired language state replay"),
    (r"\b(?:skyfield_retarded_split|py_min|py_max|calc_mode_strip_char)\b",
     "retired interpreter-specific operation"),
    N_BODY_RUST_IMPORT,
)
N_BODY_RETIRED_PATHS = (
    "src/bodies/nbody/ffi.rs", "native/rebound", "native/assist",
)
N_BODY_HELPERS = {"rebound", "assist", "librebound", "libassist"}
RUST_STRINGS = re.compile(
    r'//[^\n]*|/\*.*?\*/|r(?P<hashes>\#{0,255})"(?P<raw>.*?)"(?P=hashes)'
    r'|"(?P<quoted>(?:\\.|[^"\\])*)"', re.DOTALL,
)


def legacy_dependency(name: str) -> bool:
    """Identify legacy engines, allowing compatibility API/module names."""
    normalized = re.sub(r"[-_.]+", "-", name).lower()
    return normalized in N_BODY_HELPERS or any(
        normalized == helper + suffix
        for helper in N_BODY_HELPERS for suffix in ("-sys", "-rs")
    )


def legacy_python_import_errors(path: Path) -> list[str]:
    """Reject helper imports in delivered Python, without importing the file."""
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    except (OSError, SyntaxError) as error:
        return [str(error)]
    errors = []
    loaders = {"__import__", "import_module"}
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module == "importlib" and node.level == 0:
            loaders.update(item.asname or item.name for item in node.names if item.name == "import_module")
    for node in ast.walk(tree):
        names = []
        if isinstance(node, ast.Import):
            names = [item.name for item in node.names]
        elif isinstance(node, ast.ImportFrom):
            names = [node.module] if node.module else [item.name for item in node.names]
        elif isinstance(node, ast.Call) and node.args:
            loader = (node.func.id if isinstance(node.func, ast.Name)
                      else node.func.attr if isinstance(node.func, ast.Attribute) else None)
            if loader in loaders:
                argument = node.args[0]
                if isinstance(argument, ast.Constant) and isinstance(argument.value, str):
                    names = [argument.value]
        if any(name.split(".", 1)[0] in N_BODY_HELPERS for name in names):
            errors.append(f"{path}:{node.lineno}: retired N-body Python helper import")
    return errors


def strings(value: object) -> Iterator[str]:
    """Yield strings from the selected delivery/dependency manifest fields."""
    if isinstance(value, str):
        yield value
    elif isinstance(value, list):
        for item in value:
            yield from strings(item)
    elif isinstance(value, dict):
        for item in value.values():
            yield from strings(item)


def legacy_delivery_path(value: str) -> bool:
    """Match removed source trees, not historical names or permissive C code."""
    normalized = re.sub(r"/+", "/", value.replace("\\", "/"))
    return bool(re.search(r"(?:^|/)native/(?:rebound|assist)(?:/|$)", normalized)
                or re.search(r"(?:^|/)src/bodies/nbody/ffi\.rs$", normalized))


def nbody_delivery_errors(root: Path) -> list[str]:
    """Check normal delivery paths; isolated references and reports stay outside."""
    errors = []
    # Inspect both the current stack and misplaced legacy paths: moving the
    # source tree must never make a restored forbidden implementation invisible.
    for engine_root in (root / "core", root / "rust", root):
        for relative in N_BODY_RETIRED_PATHS:
            path = engine_root / relative
            if path.exists() or path.is_symlink():
                errors.append(f"{path}: retired N-body implementation path")
    build_scripts = {root / "build.rs", root / "core/build.rs", root / "rust/build.rs", root / "python/build.rs"}
    for relative in ("Cargo.toml", "core/Cargo.toml", "core/data/Cargo.toml",
                     "rust/Cargo.toml", "rust/data/Cargo.toml",
                     "python/Cargo.toml", "pyproject.toml",
                     "tools/leb_generation/pyproject.toml"):
        path = root / relative
        if not path.is_file():
            continue
        try:
            manifest = tomllib.loads(path.read_text(encoding="utf-8"))
        except (OSError, tomllib.TOMLDecodeError) as error:
            errors.append(f"{path}: cannot check delivery manifest: {error}")
            continue
        requirements = []
        selected_paths = []
        if path.name == "Cargo.toml":
            package = manifest.get("package", {})
            selected_paths += [package.get("include", []), package.get("build", [])]
            if isinstance(package.get("build"), str) and not legacy_delivery_path(package["build"]):
                build_scripts.add(path.parent / package["build"])

            def cargo_dependencies(table):
                for key, value in table.items():
                    if key in {"dependencies", "build-dependencies", "dev-dependencies"}:
                        for name, dependency in value.items():
                            requirements.append(dependency.get("package", name)
                                                if isinstance(dependency, dict) else name)
                            if isinstance(dependency, dict):
                                selected_paths.append(dependency.get("path", []))
                    elif isinstance(value, dict):
                        cargo_dependencies(value)

            cargo_dependencies(manifest)
        else:
            project = manifest.get("project", {})
            requirements += list(strings(project.get("dependencies", [])))
            requirements += list(strings(manifest.get("build-system", {}).get("requires", [])))
            for group in manifest.get("dependency-groups", {}).values():
                # include-group names identify another group, not a package.
                requirements += [item for item in group if isinstance(item, str)]
            selected_paths += [project.get("license-files", []),
                               manifest.get("tool", {}).get("maturin", {}).get("include", [])]
        for requirement in requirements:
            name = re.match(r"\s*([A-Za-z0-9][A-Za-z0-9._-]*)", requirement)
            if name and legacy_dependency(name[1]):
                errors.append(f"{path}: retired N-body delivery dependency {name[1]}")
        if any(legacy_delivery_path(value) for value in strings(selected_paths)):
            errors.append(f"{path}: retired N-body source/license delivery requirement")
    for path in sorted(build_scripts):
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8")
        literals = [match.group("raw") if match.group("raw") is not None else match.group("quoted")
                    for match in RUST_STRINGS.finditer(text)
                    if match.group("raw") is not None or match.group("quoted") is not None]
        if any(legacy_delivery_path(value) or legacy_dependency(value) or re.search(
                r"(?:^|\b)(?:lib)?ephemeris_(?:rebound|assist)(?:\b|$)"
                r"|rustc-link-lib=(?:(?:static|dylib)=)?(?:rebound|assist)(?:\b|$)", value)
               for value in literals):
            errors.append(f"{path}: retired N-body C compilation/link requirement")
    for directory in (root / "python/libephemeris", root / "tools/leb_generation"):
        for path in sorted(directory.rglob("*.py")):
            errors.extend(legacy_python_import_errors(path))
    return errors


def package_errors(path: Path) -> list[str]:
    """Check an import-only Python package file, without importing it."""
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    except (OSError, SyntaxError) as error:
        return [str(error)]
    errors = []
    for node in tree.body:
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            continue
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
            continue
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            name = node.targets[0].id
            if name == "__all__" and isinstance(node.value, (ast.List, ast.Tuple)) and all(
                isinstance(item, ast.Constant) and isinstance(item.value, str) for item in node.value.elts
            ):
                continue
            if name == "__version__" and isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
                continue
        errors.append(f"{path}:{node.lineno}: package logic is forbidden ({type(node).__name__})")
    return errors


def check(root: Path, phase: str = "P4") -> list[str]:
    """Return violations in the package, binding and optionally completed core."""
    errors = []
    for pattern in CI_PATHS:
        errors.extend(f"{path}: CI/CD is forbidden; run checks locally"
                      for path in sorted(root.glob(pattern)))
    package = root / "python/libephemeris"
    bindings = root / "python/src"
    if not (package / "__init__.py").is_file():
        errors.append(f"{package}: missing package initializer")
    for path in sorted(package.rglob("*.py")):
        errors.extend(package_errors(path))
    line_count = 0
    for path in sorted(bindings.rglob("*.rs")):
        text = path.read_text(encoding="utf-8")
        line_count += len(text.splitlines())
        for pattern, label in BINDING_PATTERNS:
            for match in re.finditer(pattern, text, flags=re.DOTALL):
                line = text.count("\n", 0, match.start()) + 1
                errors.append(f"{path}:{line}: {label}")
    if line_count > 20_000:
        errors.append(f"{bindings}: {line_count} lines exceeds the 20,000-line limit")
    if phase == "P4":
        errors.extend(nbody_delivery_errors(root))
        for engine_root in (root / "core", root / "rust", root):
            for pattern in CORE_PATHS:
                errors.extend(f"{path}: retired emulation path"
                              for path in sorted(engine_root.glob(pattern)))
        sources = [path for directory in (root / "core/src", root / "core/data/src", root / "rust/src", root / "rust/data/src", root / "src")
                   for path in directory.rglob("*.rs")]
        for path in sorted(sources):
            text = path.read_text(encoding="utf-8")
            if re.search(r"Accelerate|\bBLAS\b|ILP64|einsum", text, flags=re.IGNORECASE):
                errors.append(f"{path}: replicated external arithmetic ordering")
            for pattern, label in CORE_PATTERNS:
                if re.search(pattern, text):
                    errors.append(f"{path}: {label}")
        for manifest in (root / "Cargo.toml", root / "core/Cargo.toml", root / "rust/Cargo.toml"):
            if manifest.exists() and re.search(r"^harness\s*=", manifest.read_text(), flags=re.MULTILINE):
                errors.append(f"{manifest}: retired harness feature")
    return errors


def main() -> int:
    """Print failures and return a nonzero status for any boundary violation."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--phase", choices=("P0", "P4"), default="P4")
    args = parser.parse_args()
    errors = check(args.root, args.phase)
    for error in errors:
        print(error)
    if not errors:
        print(f"Minimal package checks passed ({args.phase}).")
    return bool(errors)


if __name__ == "__main__":
    raise SystemExit(main())
