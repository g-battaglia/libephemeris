"""Deterministic public-table fixtures exercise the real Rust epoch propagation."""

import copy
import json
import os
import struct
import sys
import types

import pytest

from tools.generators import star_catalog_builder as builder


def hip(ident, *, ra=10.0, dec=0.0, pmra=7200.0, pmdec=-3600.0, hp=5.2):
    return {"HIP": ident, "RArad": ra, "DErad": dec, "Plx": 2.5,
            "pmRA": pmra, "pmDE": pmdec, "Hpmag": hp}


def cached_tables(directory):
    rows = {
        "hip2": [hip(122), hip(999991), hip(999992, hp=9.0), hip(999993, hp=9.0), hip(999994, hp=9.0)],
        "hipmain": [{"HIP": ident, "HD": ident + 1000000, "Vmag": 4.25}
                    for ident in (122, 999991, 999992, 999993, 999994)],
        "cross": [{"HD": ident + 1000000, "Bayer": "alf", "Fl": None, "Cst": "And"}
                  for ident in (122, 999991)],
        "xhip": [{"HIP": ident, "RV": 0.0} for ident in (122, 999991, 999992, 999993, 999994)],
    }
    directory.mkdir()
    for key, values in rows.items():
        table = {"schema_version": 1, **copy.deepcopy(builder.QUERIES[key]),
                 "source_units": copy.deepcopy(builder.QUERIES[key]["units"]), "rows": values}
        (directory / f"{key}.json").write_text(json.dumps(table))
    # The full multiword name occupies the original fixed-width field.
    text = "# synthetic WGSN-format source rows\n"
    for name, ident in (("New Multi Word", 999992), ("Twin Name", 999993), ("Twin Name", 999994)):
        text += f"{name:<18} And {ident} 12345 10.0 0.0 extra columns\n"
    (directory / "IAU-CSN.txt").write_text(text)
    return rows


def curated_rows():
    return [[1000001, "Curated", "beAnd", 122, 1.0, 2.0, 0.1, 0.2, 3.0, -12.5, 1.25],
            [1000002, "Not In Query", "gaAnd", 145, 20.0, 30.0, 0.2, 0.3, 4.0, 2.0, 3.5]]


def test_local_join_uses_native_epochs_and_preserves_curated_identity_and_fallbacks(tmp_path):
    cache = tmp_path / "cache"
    cached_tables(cache)
    tables, names, sources = builder.acquire_inputs(cache)
    seed = curated_rows()
    rows, report = builder.build_rows(tables, names, seed)
    by_hip = {row[3]: row for row in rows}
    assert [row[3] for row in rows] == [122, 145, 999991, 999992, 999993, 999994]
    # Equatorial constant-rate invariant: 7.2/-3.6 arcsec per year for 8.75 yr.
    # The product code calls Rust, rather than repeating this formula in Python.
    assert by_hip[122][4:8] == [10.0175, -0.00875, 7.2, -3.6]
    assert by_hip[122][0:4] == seed[0][0:4]
    assert by_hip[122][9:] == [-12.5, 1.25]
    assert struct.pack("<7d", *by_hip[145][4:]) == struct.pack("<7d", *seed[1][4:])
    assert by_hip[999991][1:3] == ["alAnd", "alAnd"]
    assert by_hip[999992][1:3] == ["New Multi Word", "New Multi Word"]
    assert by_hip[999993][1:3] == ["Twin NameA", "Twin NameA"]
    assert by_hip[999994][1:3] == ["Twin NameB", "Twin NameB"]
    assert [row[0] for row in rows] == list(range(1000001, 1000007))
    assert report == {"total": 6, "iau_named": 3, "preserved": 2, "with_rv": 2}
    assert [record["key"] for record in sources] == ["hip2", "hipmain", "cross", "xhip", "wgsn"]
    assert all(len(record["sha256"]) == 64 for record in sources)


def test_curated_component_keeps_its_name_but_consumes_its_collision_letter(tmp_path):
    cache = tmp_path / "cache"
    cached_tables(cache)
    tables, names, _ = builder.acquire_inputs(cache)
    seed = curated_rows()
    seed[0][1] = "Twin Name"
    rows, _ = builder.build_rows(tables, names, seed)
    by_hip = {row[3]: row for row in rows}
    assert by_hip[122][1] == "Twin Name"
    assert by_hip[999993][1] == "Twin NameB"
    assert by_hip[999994][1] == "Twin NameC"


def test_final_valid_public_rows_and_field_specific_missing_values_keep_pinned_join_policy(tmp_path):
    cache = tmp_path / "cache"
    cached_tables(cache)
    tables, names, _ = builder.acquire_inputs(cache)
    tables["hipmain"].append({"HIP": 999991, "HD": None, "Vmag": None})
    tables["hip2"].append(hip(999991, ra=11.0))
    invalid = hip(999991, ra=None)
    tables["hip2"].append(invalid)
    rows, _ = builder.build_rows(tables, names, curated_rows())
    actual = next(row for row in rows if row[3] == 999991)
    assert actual[1:3] == ["alAnd", "alAnd"]
    assert actual[4] == 11.0175 and actual[10] == 4.25


@pytest.mark.parametrize("bayer,number,constellation,expected", [
    ("alf", None, "And", "alAnd"), ("the1", None, "Ori", "th01Ori"),
    ("mu.", None, "And", "muAnd"), ("omi", None, "And", "omiAnd"),
    ("ome", None, "And", "omeAnd"), (None, 32, "Leo", "32Leo"),
    ("invalid", 32, "Leo", None), ("alf", None, "", None),
])
def test_bayer_and_component_spelling(bayer, number, constellation, expected):
    assert builder.bayer_key(bayer, number, constellation) == expected


def test_native_polar_and_signed_zero_behavior_is_called_directly(tmp_path):
    from libephemeris.fixed_stars import propagate_proper_motion
    cache = tmp_path / "cache"
    cached_tables(cache)
    tables, names, _ = builder.acquire_inputs(cache)
    tables["hip2"][0] = hip(122, ra=-0.0, dec=90.0, pmra=1.0, pmdec=1000.0)
    rows, _ = builder.build_rows(tables, names, curated_rows())
    expected = propagate_proper_motion(-0.0, 90.0, 0.001, 1.0, builder.FROM_JD, builder.TO_JD)
    actual = next(row for row in rows if row[3] == 122)
    assert struct.pack("<2d", *actual[4:6]) == struct.pack("<2d", *(round(value, 9) for value in expected))


@pytest.mark.parametrize("invalid", ["catalog", "query", "units", "source-units", "columns", "nan", "coordinate-range"])
def test_local_cache_refuses_wrong_source_contract_or_nonfinite_values(tmp_path, invalid):
    cache = tmp_path / "cache"
    cached_tables(cache)
    path = cache / "hip2.json"
    table = json.loads(path.read_text())
    if invalid == "catalog":
        table["catalog"] = "different catalog"
    elif invalid == "query":
        table["constraints"] = {}
    elif invalid == "units":
        table["units"]["RArad"] = "rad"
    elif invalid == "source-units":
        table["source_units"].clear()
    elif invalid == "columns":
        del table["rows"][0]["Plx"]
    elif invalid == "nan":
        table["rows"][0]["Plx"] = float("nan")
    elif invalid == "coordinate-range":
        table["rows"][0]["DErad"] = 91.0
    path.write_text(json.dumps(table))
    with pytest.raises(ValueError):
        builder.acquire_inputs(cache)


def test_complete_candidate_validates_and_cache_reuse_never_fetches(tmp_path, monkeypatch):
    cache = tmp_path / "cache"
    cached_tables(cache)
    before = {path: path.read_bytes() for path in cache.iterdir()}
    def unexpected_fetch(key):
        raise AssertionError(f"unexpected remote query: {key}")
    monkeypatch.setattr(builder, "_fetch_table", unexpected_fetch)
    output = tmp_path / "candidate.json"
    assert builder.main(["--cache-dir", str(cache), "--out", str(output), "--fetch"]) == 0
    candidate, _ = builder.serializer.load_catalog(output)
    seed, _ = builder.serializer.load_catalog(builder.serializer.DEFAULT_INPUT)
    assert candidate["schema_version"] == 2
    assert len(candidate["rows"]) == len(seed["rows"]) + 4
    assert {row[0] for row in seed["rows"]} <= {row[0] for row in candidate["rows"]}
    for key in ("aliases", "name_to_hip", "hip_fix", "greek", "constellations", "decimal_ranges", "lineage", "photometry"):
        assert candidate[key] == seed[key]
    assert len(candidate["build"]["propagation"]["extension_sha256"]) == 64
    assert before == {path: path.read_bytes() for path in before}
    previous = output.read_bytes()
    assert builder.main(["--cache-dir", str(cache), "--out", str(output)]) == 2
    assert output.read_bytes() == previous
    assert builder.main(["--cache-dir", str(cache), "--out", str(output), "--overwrite"]) == 0
    assert output.read_bytes() == previous


def test_integral_json_measurements_are_emitted_as_binary64(tmp_path):
    cache = tmp_path / "cache"
    cached_tables(cache)
    for key, fields in (("hip2", ("Plx", "Hpmag")), ("hipmain", ("Vmag",)), ("xhip", ("RV",))):
        path = cache / f"{key}.json"
        table = json.loads(path.read_bytes())
        for row in table["rows"]:
            for field in fields:
                row[field] = 2
        path.write_text(json.dumps(table))
    output = tmp_path / "candidate.json"
    assert builder.main(["--cache-dir", str(cache), "--out", str(output)]) == 0
    candidate, _ = builder.serializer.load_catalog(output)
    assert all(type(value) is float for row in candidate["rows"] for value in row[4:])
    new = next(row for row in candidate["rows"] if row[3] == 999991)
    assert new[8:] == [2.0, 2.0, 2.0]


def test_output_aliases_missing_inputs_and_late_source_changes_preserve_previous_data(tmp_path, monkeypatch):
    cache = tmp_path / "cache"
    cached_tables(cache)
    source = cache / "hip2.json"
    original = source.read_bytes()
    assert builder.main(["--cache-dir", str(cache), "--out", str(source), "--overwrite"]) == 2
    alias = tmp_path / "alias.json"
    os.link(source, alias)
    assert builder.main(["--cache-dir", str(cache), "--out", str(alias), "--overwrite"]) == 2
    assert source.read_bytes() == original
    output = tmp_path / "previous.json"
    output.write_bytes(b"previous candidate")
    original_build = builder.build_rows
    def change_input(*args):
        result = original_build(*args)
        source.write_bytes(original + b"\n")
        return result
    monkeypatch.setattr(builder, "build_rows", change_input)
    assert builder.main(["--cache-dir", str(cache), "--out", str(output), "--overwrite"]) == 2
    assert output.read_bytes() == b"previous candidate"
    source.unlink()
    assert builder.main(["--cache-dir", str(cache), "--out", str(output), "--overwrite"]) == 2
    assert output.read_bytes() == b"previous candidate"


@pytest.mark.parametrize("filename", ["hip2.json", "IAU-CSN.txt"])
def test_hash_binds_exact_parsed_bytes_even_if_path_changes_immediately_after_read(tmp_path, monkeypatch, filename):
    cache = tmp_path / "cache"
    cached_tables(cache)
    changed = cache / filename
    output = tmp_path / "previous.json"
    output.write_bytes(b"previous candidate")
    original_read = builder.Path.read_bytes
    replaced = False
    def change_after_read(path):
        nonlocal replaced
        raw = original_read(path)
        if path == changed and not replaced:
            replaced = True
            path.write_bytes(raw + b"\n")
        return raw
    monkeypatch.setattr(builder.Path, "read_bytes", change_after_read)
    assert builder.main(["--cache-dir", str(cache), "--out", str(output), "--overwrite"]) == 2
    assert replaced
    assert output.read_bytes() == b"previous candidate"


def test_native_extension_identity_is_captured_before_propagation_and_rechecked(tmp_path, monkeypatch):
    import libephemeris._core as core

    cache = tmp_path / "cache"
    cached_tables(cache)
    output = tmp_path / "previous.json"
    output.write_bytes(b"previous candidate")
    extension = builder.Path(core.__file__)
    digest = builder._digest
    observed = []

    def changed_identity(path):
        actual = digest(path)
        if builder.Path(path) == extension:
            observed.append(actual)
            return actual if len(observed) == 1 else "0" * 64
        return actual

    original_build = builder.build_rows

    def checked_build(*args):
        assert len(observed) == 1
        return original_build(*args)

    monkeypatch.setattr(builder, "_digest", changed_identity)
    monkeypatch.setattr(builder, "build_rows", checked_build)
    assert builder.main(["--cache-dir", str(cache), "--out", str(output), "--overwrite"]) == 2
    assert len(observed) == 2
    assert output.read_bytes() == b"previous candidate"


@pytest.mark.parametrize("alias", [False, True])
def test_candidate_output_cannot_replace_its_native_extension(tmp_path, monkeypatch, alias):
    import libephemeris
    import libephemeris._core

    cache = tmp_path / "cache"
    cached_tables(cache)
    # A fake file makes the destructive regression safe to exercise. Native
    # propagation remains the real already-imported function, if reached.
    extension = tmp_path / "native.so"
    extension.write_bytes(b"native extension identity")
    fake = types.SimpleNamespace(__file__=str(extension))
    monkeypatch.setitem(sys.modules, "libephemeris._core", fake)
    monkeypatch.setattr(libephemeris, "_core", fake)
    output = extension
    if alias:
        output = tmp_path / "alias.so"
        os.link(extension, output)
    assert builder.main(["--cache-dir", str(cache), "--out", str(output), "--overwrite"]) == 2
    assert extension.read_bytes() == output.read_bytes() == b"native extension identity"


@pytest.mark.parametrize("target", ["snapshot", "rust"])
def test_explicit_alternative_seed_still_protects_both_shipped_catalog_files(tmp_path, monkeypatch, target):
    original_seed = builder.serializer.DEFAULT_INPUT
    snapshot = tmp_path / "shipped.json"
    rust = tmp_path / "core/src/tables/star_catalog_tables.rs"
    rust.parent.mkdir(parents=True)
    snapshot.write_bytes(b"shipped snapshot")
    rust.write_bytes(b"shipped Rust table")
    monkeypatch.setattr(builder.serializer, "DEFAULT_INPUT", snapshot)
    monkeypatch.setattr(builder.serializer, "REPO_ROOT", tmp_path)
    output = snapshot if target == "snapshot" else rust
    assert builder.main(["--cache-dir", str(tmp_path / "unused"), "--seed", str(original_seed),
                         "--out", str(output), "--overwrite"]) == 2
    assert snapshot.read_bytes() == b"shipped snapshot"
    assert rust.read_bytes() == b"shipped Rust table"


def test_fetch_normalizes_declared_astropy_units_instead_of_assuming_field_names(monkeypatch):
    calls = []
    class Quantity:
        def __init__(self, value, unit):
            self.value, self.unit = value, unit
        def to_value(self, target):
            calls.append((self.unit, target))
            return 30.0 if self.unit == "rad" else self.value
    class Unit:
        def __init__(self, name):
            self.name = name
        def __str__(self):
            return self.name
        def __rmul__(self, value):
            return Quantity(value, self.name)
    class Table:
        def __iter__(self):
            return iter([hip(122, ra=0.5, dec=0.5)])
        def __getitem__(self, key):
            name = "rad" if key in ("RArad", "DErad") else builder.QUERIES["hip2"]["units"][key]
            return types.SimpleNamespace(unit=Unit(name))
    class Vizier:
        def __init__(self, **kwargs):
            assert kwargs == {"columns": builder.QUERIES["hip2"]["columns"], "row_limit": -1}
        def query_constraints(self, **kwargs):
            assert kwargs == {"catalog": "I/311/hip2", "Hpmag": "<7.2", "cache": False}
            return [Table()]
    monkeypatch.setitem(sys.modules, "astroquery", types.ModuleType("astroquery"))
    monkeypatch.setitem(sys.modules, "astroquery.vizier", types.SimpleNamespace(Vizier=Vizier))
    normalized = builder._fetch_table("hip2")
    assert normalized["rows"][0]["RArad"] == normalized["rows"][0]["DErad"] == 30.0
    assert normalized["source_units"]["RArad"] == "rad"
    assert calls.count(("rad", "deg")) == 2


def test_real_astropy_columns_convert_radians_and_proper_motion_units(monkeypatch):
    import math
    from astropy import units
    from astropy.table import Table

    table = Table(rows=[hip(122, ra=0.5, dec=0.25, pmra=1000.0, pmdec=-2000.0)])
    for field, name in builder.QUERIES["hip2"]["units"].items():
        table[field].unit = units.Unit("rad" if field in ("RArad", "DErad") else name)

    class Vizier:
        def __init__(self, **kwargs):
            pass

        def query_constraints(self, **kwargs):
            return [table]

    monkeypatch.setitem(sys.modules, "astroquery.vizier", types.SimpleNamespace(Vizier=Vizier))
    normalized = builder._fetch_table("hip2")
    row = normalized["rows"][0]
    assert row["RArad"] == pytest.approx(0.5 * 180 / math.pi, abs=1e-13)
    assert row["DErad"] == pytest.approx(0.25 * 180 / math.pi, abs=1e-13)
    assert row["pmRA"] == 1000.0
    assert row["pmDE"] == -2000.0
    assert normalized["units"]["RArad"] == "deg"
    assert normalized["source_units"]["RArad"] == "rad"


def test_curated_h5_values_and_sparse_provenance_survive_a_public_rebuild(tmp_path):
    cache = tmp_path / "cache"
    cached_tables(cache)
    tables, names, _ = builder.acquire_inputs(cache)
    seed, _ = builder.serializer.load_catalog(builder.serializer.DEFAULT_INPUT)
    curated = [row for row in seed["rows"] if row[3] in (65474, 80763)]
    for hip_number, raw_vmag in ((65474, 0.9800000190734863), (80763, 1.059999942779541)):
        tables["hip2"].append(hip(hip_number))
        tables["hipmain"].append({"HIP": hip_number, "HD": None, "Vmag": raw_vmag})
    rows, _ = builder.build_rows(tables, names, curated)
    by_hip = {row[3]: row for row in rows}
    for annotation in seed["photometry"]:
        actual = by_hip[annotation["hip_number"]]
        assert actual[0] == annotation["star_id"]
        assert actual[10].hex() == annotation["value"].hex()
    assert by_hip[65474][10] == 0.98
    assert by_hip[80763][10] == 1.06
