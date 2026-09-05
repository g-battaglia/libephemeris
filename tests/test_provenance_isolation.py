# SPDX-License-Identifier: AGPL-3.0-only
"""Structural isolation guards for the shipped package.

The reference API may be used only for ephemeral behavioral comparison
from the private validation harness. The shipped package must therefore
never import the reference binding, reference the validation tree, or
read reference-distribution data files. These guards make that contract
executable.
"""

from __future__ import annotations

import pathlib
import re

PKG = pathlib.Path(__file__).resolve().parents[1] / "libephemeris"


def _sources():
    return sorted(PKG.rglob("*.py"))


def test_package_never_imports_the_reference_binding():
    pat = re.compile(r"^\s*(import\s+swisseph|from\s+swisseph)", re.M)
    hits = [p.name for p in _sources() if pat.search(p.read_text())]
    assert hits == [], f"reference binding imported by: {hits}"


def test_package_never_references_the_validation_tree():
    pat = re.compile(r"validation/(data|compare_scripts|\.venv)")
    hits = [p.name for p in _sources() if pat.search(p.read_text())]
    assert hits == [], f"validation tree referenced by: {hits}"


def test_package_never_opens_reference_data_files():
    # Reference-distribution compressed-ephemeris families must never be
    # read by the package. (The API-compat constant FICTFILE is name-level
    # API surface: the bundled path resolves to this project's own
    # fictitious_orbits.csv, checked below.)
    pat = re.compile(r"""['"][^'"]*(?:\.se1|sepl_|semo_|seas_)['"]""")
    hits = [p.name for p in _sources() if pat.search(p.read_text())]
    assert hits == [], f"reference data filename referenced by: {hits}"


def test_no_reference_data_file_is_shipped():
    data = PKG / "data"
    shipped = [q.name for q in data.rglob("*") if q.is_file()]
    for name in shipped:
        assert not name.endswith(".se1") and name != "seorbel.txt", name
    # The bundled path must resolve to this project's own data set.
    from libephemeris.hypothetical import get_bundled_fictitious_orbits_path

    assert get_bundled_fictitious_orbits_path().name == "fictitious_orbits.csv"
