#!/usr/bin/env python3
"""Integrity gate for lunar-model architecture and integrity-pinned files.

Standards/JPL implementations are checked structurally.  Integrity-pinned files
are accepted only at their exact reviewed SHA-256; any modification changes the
hash and fails this gate.

Usage:
    uv run python scripts/check_lunar_model_provenance.py
"""

from __future__ import annotations

import ast
import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

RETIRED_PATHS = (
    "libephemeris/data/mean_lunar_apse.bin",
    "libephemeris/data/mean_lunar_apse.json",
    "libephemeris/data/interpolated_lunar_apsides.bin",
    "scripts/generate_mean_lunar_apse.py",
    "scripts/generate_interpolated_lunar_apsides.py",
)

PINNED_PATHS = {
    "libephemeris/lunar.py": (
        "b69348b6f3bf90f32aeaa1e34fbaaa1541d94c22147795bad38aff890c87cbff"
    ),
    "libephemeris/lunar_apse_model.py": (
        "c034bc9b7cde9aee7ee3d72c39a868371c5d275616f543d7fa837f5bd41ccf17"
    ),
}

LUNAR_FUNCTION_REQUIREMENTS = {
    "calc_mean_lunar_node": ("mean_lunar_node_position",),
    "calc_mean_lilith": ("mean_lunar_apogee_position",),
    "calc_mean_lilith_with_latitude": ("mean_lunar_apogee_position",),
    "calc_mean_lunar_node_state": ("mean_lunar_node_state",),
    "calc_mean_lilith_state": ("mean_lunar_apogee_state",),
    "calc_interpolated_apogee": (
        "_compat_legacy_mean_lilith",
        "_compat_legacy_mean_lunar_node",
    ),
    "calc_interpolated_perigee": (
        "_compat_legacy_mean_lilith",
        "_compat_legacy_mean_lunar_node",
    ),
}

SOURCE_REQUIREMENTS = {
    "libephemeris/mean_lunar_apse.py": (
        "erfa.fal03",
        "erfa.faf03",
        "erfa.faom03",
    ),
    "libephemeris/interpolated_lunar_apsides.py": (
        "_moon_icrs_state",
        "_osculating_eccentricity_and_p",
        "np.polynomial.legendre.leggauss",
    ),
}

RETIRED_RUNTIME_TOKENS = (
    "MODEL_FILE_SHA256",
    "MODEL_PAYLOAD_SHA256",
    "MODEL_SHA256",
    "mean_lunar_apse.bin",
    "interpolated_lunar_apsides.bin",
)


def main() -> int:
    problems: list[str] = []
    for relative in RETIRED_PATHS:
        if (ROOT / relative).exists():
            problems.append(f"retired lunar artifact still exists: {relative}")

    for relative, expected_sha256 in PINNED_PATHS.items():
        path = ROOT / relative
        if not path.is_file():
            problems.append(f"pinned compatibility artifact is missing: {relative}")
            continue
        actual_sha256 = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual_sha256 != expected_sha256:
            problems.append(
                f"{relative} SHA-256 is {actual_sha256}, expected {expected_sha256}"
            )

    for relative, requirements in SOURCE_REQUIREMENTS.items():
        text = (ROOT / relative).read_text(encoding="utf-8")
        for requirement in requirements:
            if requirement not in text:
                problems.append(f"{relative} is missing {requirement!r}")
        for token in RETIRED_RUNTIME_TOKENS:
            if token in text:
                problems.append(f"{relative} still references {token!r}")

    lunar_path = ROOT / "libephemeris/lunar.py"
    lunar_text = lunar_path.read_text(encoding="utf-8")
    lunar_tree = ast.parse(lunar_text)
    functions = {
        node.name: ast.get_source_segment(lunar_text, node) or ""
        for node in lunar_tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    for function_name, requirements in LUNAR_FUNCTION_REQUIREMENTS.items():
        function_source = functions.get(function_name)
        if function_source is None:
            problems.append(f"libephemeris/lunar.py is missing {function_name!r}")
            continue
        for requirement in requirements:
            if requirement not in function_source:
                problems.append(
                    "libephemeris/lunar.py function "
                    f"{function_name!r} is missing {requirement!r}"
                )

    for problem in problems:
        print(f"MISMATCH: {problem}")
    print(f"lunar model provenance: {len(problems)} mismatch(es) (gate: 0)")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
