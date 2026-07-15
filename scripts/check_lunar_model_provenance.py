#!/usr/bin/env python3
"""Integrity gate for lunar-model architecture and integrity-pinned files.

Standards/JPL implementations are checked structurally.  Integrity-pinned files
are accepted only at their exact reviewed SHA-256; any modification changes the
hash and fails this gate.

Usage:
    uv run python scripts/check_lunar_model_provenance.py

Provenance:
    Project-authored release gate over the documented IERS/ERFA/JPL lunar-model
    boundary. Pinned hashes attest that reviewed generated/source files did not
    change; they do not replace the publications, data provenance, or numerical
    validation. The checker parses project files only and creates no lunar
    coefficient.
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
    "libephemeris/lunar_apse_model.py": (
        "7ce63d2f390efceab53e6e27a8d4192ff23066af44f6af97360d70f82b4615ba"
    ),
}

LUNAR_FUNCTION_REQUIREMENTS = {
    "calc_mean_lunar_node": ("mean_lunar_node_position",),
    "calc_mean_lilith": ("mean_lunar_apogee_position",),
    "calc_mean_lilith_with_latitude": ("mean_lunar_apogee_position",),
    "calc_mean_lunar_node_state": ("mean_lunar_node_state",),
    "calc_mean_lilith_state": ("mean_lunar_apogee_state",),
    "calc_interpolated_apogee": (
        "_mean_lunar_apogee_position_unchecked",
        "_mean_lunar_node_position_unchecked",
    ),
    "calc_interpolated_perigee": (
        "_mean_lunar_apogee_position_unchecked",
        "_mean_lunar_node_position_unchecked",
    ),
    "_calc_de440_apogee_passage_terms": ("lunar_delaunay_arguments",),
    "_calc_de440_perigee_passage_terms": ("lunar_delaunay_arguments",),
}

SOURCE_REQUIREMENTS = {
    "libephemeris/mean_lunar_apse.py": (
        "erfa.fal03",
        "erfa.faf03",
        "erfa.faom03",
        "erfa.fad03",
        "erfa.falp03",
        "EARTH_ECCENTRICITY_J2000 = 0.016708634",
        "Simon et al. (1994)",
    ),
    "libephemeris/interpolated_lunar_apsides.py": (
        "_moon_icrs_state",
        "_osculating_eccentricity_and_p",
        "np.polynomial.legendre.leggauss",
    ),
    "scripts/generate_lunar_apse_model.py": (
        "lunar_delaunay_arguments",
        "mean_lunar_apogee_position",
        "mean_lunar_node_position",
        "denum != 440",
        "EARTH_ECCENTRICITY_J2000",
    ),
}

RETIRED_RUNTIME_TOKENS = (
    "MODEL_FILE_SHA256",
    "MODEL_PAYLOAD_SHA256",
    "MODEL_SHA256",
    "mean_lunar_apse.bin",
    "interpolated_lunar_apsides.bin",
)

RETIRED_LUNAR_FUNCTIONS = (
    "_calc_mean_apse_analytical",
    "_compat_legacy_mean_lilith",
    "_compat_legacy_mean_lunar_node",
)


def main() -> int:
    """Verify lunar generators, source markers, hashes, and retired boundaries."""
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
    for function_name in RETIRED_LUNAR_FUNCTIONS:
        if function_name in functions:
            problems.append(
                f"libephemeris/lunar.py still defines retired {function_name!r}"
            )
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
