"""Standalone serialization preserves the shipped catalog, not new astronomy."""

import ast
import copy
import hashlib
import json
import os
import re
import shutil
import struct
import subprocess
import sys
from pathlib import Path

import pytest

from tools.generators import gen_star_catalog_tables as generator


ROOT = Path(__file__).resolve().parents[2]
RUST_TABLE = ROOT / "core/src/tables/star_catalog_tables.rs"


def array(text, name):
    match = re.search(rf"pub\(crate\) const {name}:.*?= &(\[.*?\]);", text, re.S)
    assert match is not None, name
    return ast.literal_eval(match.group(1))


def row_values(text):
    return [ast.literal_eval("(" + match + ")") for match in
            re.findall(r"StarCatalogEntry::from_row\((.*?)\),", text, re.S)]


def test_canonical_measurements_and_ordered_maps_match_the_shipped_rust_snapshot():
    catalog, digest = generator.load_catalog(generator.DEFAULT_INPUT)
    actual = RUST_TABLE.read_text()
    rendered = generator.render_catalog(catalog, digest)
    original_rows, new_rows = row_values(actual), row_values(rendered)
    assert len(original_rows) == len(new_rows) == len(catalog["rows"]) == 1447
    for original, new, canonical in zip(original_rows, new_rows, catalog["rows"], strict=True):
        assert original[:4] == new[:4] == tuple(canonical[:4])
        assert struct.pack("<7d", *original[4:]) == struct.pack("<7d", *new[4:])
        assert struct.pack("<7d", *original[4:]) == struct.pack("<7d", *canonical[4:])
    for name in ("REF_SORTED_IDS", "NOMEN_SORTED_IDS", "NAME_TO_ID_UPPER", "ALIASES_UPPER",
                 "NOMEN_UPPER_KEYS", "STAR_ALIASES", "STAR_NAME_TO_HIP", "NAME_HIP_FIX",
                 "GREEK_LETTER_ABBREV", "CONSTELLATION_ABBREV", "UNICODE_ND_RANGES"):
        assert array(actual, name) == array(rendered, name), name
    assert array(actual, "STAR_ALIASES") == [tuple(pair) for pair in catalog["aliases"]]
    assert len(catalog["aliases"]) == 820
    assert ("ABSOLUTNO", -1) in array(rendered, "STAR_NAME_TO_HIP")
    assert digest == hashlib.sha256(generator.DEFAULT_INPUT.read_bytes()).hexdigest()


def test_ordered_duplicate_aliases_are_not_collapsed_in_the_canonical_input(tmp_path):
    catalog, _ = generator.load_catalog(generator.DEFAULT_INPUT)
    catalog["aliases"][1] = [catalog["aliases"][0][0], catalog["rows"][0][0]]
    source = tmp_path / "duplicates.json"
    source.write_text(json.dumps(catalog))
    decoded, digest = generator.load_catalog(source)
    rendered = generator.render_catalog(decoded, digest)
    assert array(rendered, "STAR_ALIASES")[:2] == [tuple(pair) for pair in catalog["aliases"][:2]]
    assert dict(array(rendered, "ALIASES_UPPER"))[catalog["aliases"][0][0].upper()] == catalog["rows"][0][0]


def test_default_generator_and_frozen_unicode_work_without_git_or_an_engine(tmp_path):
    output = tmp_path / "tables.rs"
    # Copy only the serializer and its data into a separate minimal tree.
    # Default discovery must work without the original checkout or its neighbors.
    directory = tmp_path / "standalone/tools/generators"
    (directory / "data").mkdir(parents=True)
    script = directory / "gen_star_catalog_tables.py"
    shutil.copyfile(ROOT / "tools/generators/gen_star_catalog_tables.py", script)
    shutil.copyfile(generator.DEFAULT_INPUT, directory / "data/star_catalog.json")
    environment = {**os.environ, "PATH": "", "PYTHONPATH": "", "PYTHONDONTWRITEBYTECODE": "1"}
    command = [sys.executable, "-I", "-B", str(script), "--out", str(output)]
    result = subprocess.run(command, cwd=tmp_path, env=environment, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "1447 rows, 820 aliases" in result.stdout
    assert "Unicode 16.0.0" in output.read_text()
    result = subprocess.run([*command, "--check"], cwd=tmp_path, env=environment,
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert len(array(output.read_text(), "UNICODE_ND_RANGES")) == 71


def test_explicit_input_check_and_drift_do_not_rewrite_files(tmp_path, capsys):
    source = tmp_path / "catalog.json"
    source.write_bytes(generator.DEFAULT_INPUT.read_bytes())
    output = tmp_path / "tables.rs"
    arguments = ["--input", str(source), "--out", str(output)]
    assert generator.main(arguments) == 0
    original = output.read_bytes()
    assert generator.main([*arguments, "--check"]) == 0
    assert output.read_bytes() == original
    output.write_bytes(original + b"// drift\n")
    assert generator.main([*arguments, "--check"]) == 1
    assert output.read_bytes() == original + b"// drift\n"
    assert "regeneration drift" in capsys.readouterr().err


@pytest.mark.parametrize("invalid", [
    "schema", "unknown-field", "row-count", "row-type", "boolean-id", "nan", "infinity",
    "duplicate-id", "duplicate-name", "aliases-count", "absent-alias-id", "duplicate-map-key",
    "control-string", "unicode-version", "unicode-count", "unicode-order", "unicode-width", "lineage",
    "lineage-hash", "lineage-concatenated-hash", "lineage-fields",
])
def test_invalid_inputs_preserve_existing_output(tmp_path, invalid):
    catalog, _ = generator.load_catalog(generator.DEFAULT_INPUT)
    catalog = copy.deepcopy(catalog)
    if invalid == "schema":
        catalog["schema_version"] = True
    elif invalid == "unknown-field":
        catalog["unrecognized"] = []
    elif invalid == "row-count":
        catalog["rows"].pop()
    elif invalid == "row-type":
        catalog["rows"][0][4] = 1
    elif invalid == "boolean-id":
        catalog["rows"][0][0] = True
    elif invalid == "nan":
        catalog["rows"][0][4] = float("nan")
    elif invalid == "infinity":
        catalog["rows"][0][4] = float("inf")
    elif invalid == "duplicate-id":
        catalog["rows"][1][0] = catalog["rows"][0][0]
    elif invalid == "duplicate-name":
        catalog["rows"][1][1] = catalog["rows"][0][1].upper()
    elif invalid == "aliases-count":
        catalog["aliases"].pop()
    elif invalid == "absent-alias-id":
        catalog["aliases"][0][1] = 42
    elif invalid == "duplicate-map-key":
        catalog["name_to_hip"][1][0] = catalog["name_to_hip"][0][0]
    elif invalid == "control-string":
        catalog["rows"][0][1] = "name\n"
    elif invalid == "unicode-version":
        catalog["unicode_version"] = "17.0.0"
    elif invalid == "unicode-count":
        catalog["decimal_ranges"].pop()
    elif invalid == "unicode-order":
        catalog["decimal_ranges"][1] = catalog["decimal_ranges"][0]
    elif invalid == "unicode-width":
        catalog["decimal_ranges"][0][1] -= 1
    elif invalid == "lineage":
        catalog["lineage"]["seed_commit"] = "not a commit"
    elif invalid == "lineage-hash":
        catalog["lineage"]["seed_blobs_sha256"]["libephemeris/constants.py"] = "z" * 64
    elif invalid == "lineage-concatenated-hash":
        catalog["lineage"]["seed_concatenated_sha256"] = 123
    elif invalid == "lineage-fields":
        del catalog["lineage"]["note"]
    source, output = tmp_path / "invalid.json", tmp_path / "preserved.rs"
    source.write_text(json.dumps(catalog))
    output.write_bytes(b"previous output")
    assert generator.main(["--input", str(source), "--out", str(output)]) == 2
    assert output.read_bytes() == b"previous output"


@pytest.mark.parametrize("payload", ["[]", "null", "true", "123", "{}", '{"rows":'])
def test_nonobject_or_malformed_json_is_a_clear_input_error(tmp_path, payload, capsys):
    source, output = tmp_path / "invalid.json", tmp_path / "preserved.rs"
    source.write_text(payload)
    output.write_bytes(b"previous output")
    assert generator.main(["--input", str(source), "--out", str(output)]) == 2
    assert "invalid catalog input/output" in capsys.readouterr().err
    assert output.read_bytes() == b"previous output"


def test_duplicate_json_fields_and_input_output_aliases_are_errors(tmp_path):
    source = tmp_path / "duplicate.json"
    source.write_text('{"schema_version": 1, "schema_version": 1}')
    with pytest.raises(ValueError, match="duplicate JSON field"):
        generator.load_catalog(source)
    source.write_bytes(generator.DEFAULT_INPUT.read_bytes())
    original = source.read_bytes()
    assert generator.main(["--input", str(source), "--out", str(source)]) == 2
    alias = tmp_path / "alias.json"
    os.link(source, alias)
    assert generator.main(["--input", str(source), "--out", str(alias)]) == 2
    assert source.read_bytes() == alias.read_bytes() == original


def public_candidate():
    """A synthetic acquisition record exercises the explicit candidate format."""
    catalog, _ = generator.load_catalog(generator.DEFAULT_INPUT)
    catalog["schema_version"] = 2
    catalog["rows"].append([1002000, "NewTestStar", "testNew", 1234567,
                            10.0, 20.0, 0.0, -0.0, 0.0, 0.0, 5.0])
    sources = [{"key": key, "catalog": name, "columns": ["HIP"],
                "constraints": {}, "row_limit": -1, "units": {}, "source_units": {},
                "input_path": f"{key}.json", "sha256": "a" * 64}
               for key, name in (("hip2", "I/311/hip2"), ("hipmain", "I/239/hip_main"),
                                 ("cross", "IV/27A/catalog"), ("xhip", "V/137D/XHIP"))]
    sources.append({"key": "wgsn", "url": "https://example.test/wgsn", "format": "test",
                    "input_path": "wgsn.txt", "sha256": "b" * 64})
    catalog["build"] = {
        "builder": "tools/generators/star_catalog_builder.py", "seed_sha256": "c" * 64,
        "sources": sources,
        "propagation": {"function": "libephemeris.fixed_stars.propagate_proper_motion",
                        "from_jd": 2448349.0625, "to_jd": 2451545.0,
                        "extension_sha256": "d" * 64},
        "selection": {"hp_limit": 5.2, "round_digits": [9, 9, 7, 7, 3, 2, 3]},
        "report": {"total": 1448, "iau_named": 1, "preserved": 1447, "with_rv": 0},
    }
    return catalog


def test_public_candidate_supports_new_rows_without_weakening_the_frozen_snapshot(tmp_path):
    catalog = public_candidate()
    source = tmp_path / "public.json"
    source.write_text(json.dumps(catalog))
    decoded, digest = generator.load_catalog(source)
    text = generator.render_catalog(decoded, digest)
    assert len(row_values(text)) == 1448
    assert "public-catalog candidate; review required" in text
    assert f"Public-catalog input SHA-256: `{digest}`" in text
    assert "from 1448 rows" in text
    assert array(text, "STAR_ALIASES") == [tuple(pair) for pair in catalog["aliases"]]
    catalog["schema_version"] = 1
    del catalog["build"]
    source.write_text(json.dumps(catalog))
    with pytest.raises(ValueError, match="rows requires"):
        generator.load_catalog(source)


def test_public_candidates_cannot_write_the_shipped_rust_table(tmp_path, monkeypatch):
    source = tmp_path / "public.json"
    source.write_text(json.dumps(public_candidate()))
    canonical = tmp_path / "core/src/tables/star_catalog_tables.rs"
    canonical.parent.mkdir(parents=True)
    canonical.write_bytes(b"existing shipped table")
    snapshot = tmp_path / "shipped.json"
    snapshot.write_bytes(b"existing shipped snapshot")
    monkeypatch.setattr(generator, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(generator, "DEFAULT_INPUT", snapshot)
    assert generator.main(["--input", str(source)]) == 2
    alias = tmp_path / "alias.rs"
    os.link(canonical, alias)
    assert generator.main(["--input", str(source), "--out", str(alias)]) == 2
    assert canonical.read_bytes() == alias.read_bytes() == b"existing shipped table"
    assert generator.main(["--input", str(source), "--out", str(snapshot)]) == 2
    snapshot_alias = tmp_path / "snapshot-alias.json"
    os.link(snapshot, snapshot_alias)
    assert generator.main(["--input", str(source), "--out", str(snapshot_alias)]) == 2
    assert snapshot.read_bytes() == snapshot_alias.read_bytes() == b"existing shipped snapshot"
    output = tmp_path / "review.rs"
    assert generator.main(["--input", str(source), "--out", str(output)]) == 0
    assert len(row_values(output.read_text())) == 1448


@pytest.mark.parametrize("invalid", [
    "missing-build", "source-count", "duplicate-source", "source-hash", "source-units",
    "source-constraints", "native-function", "native-epoch", "selection", "count", "boolean-count",
])
def test_public_candidates_require_complete_typed_provenance(tmp_path, invalid):
    catalog = public_candidate()
    build = catalog["build"]
    if invalid == "missing-build":
        del catalog["build"]
    elif invalid == "source-count":
        build["sources"].pop()
    elif invalid == "duplicate-source":
        build["sources"][1] = build["sources"][0]
    elif invalid == "source-hash":
        build["sources"][0]["sha256"] = "bad"
    elif invalid == "source-units":
        build["sources"][0]["units"] = {"RArad": 1}
    elif invalid == "source-constraints":
        build["sources"][0]["constraints"] = []
    elif invalid == "native-function":
        build["propagation"]["function"] = "python.math"
    elif invalid == "native-epoch":
        build["propagation"]["from_jd"] = 2451545.0
    elif invalid == "selection":
        build["selection"]["round_digits"] = [6] * 7
    elif invalid == "count":
        build["report"]["total"] -= 1
    elif invalid == "boolean-count":
        build["report"]["preserved"] = True
    source = tmp_path / "invalid.json"
    source.write_text(json.dumps(catalog))
    with pytest.raises(ValueError):
        generator.load_catalog(source)


def test_spica_h5_correction_preserves_all_other_canonical_words_and_lookup_data():
    catalog, _ = generator.load_catalog(generator.DEFAULT_INPUT)
    by_hip = {row[3]: row for row in catalog["rows"]}
    assert by_hip[65474][10].hex() == (0.98).hex()
    assert by_hip[80763][10].hex() == (1.06).hex()
    annotations = {item["hip_number"]: item for item in catalog["photometry"]}
    assert set(annotations) == {65474, 80763}
    for hip, source_value in ((65474, 0.9800000190734863), (80763, 1.059999942779541)):
        item = annotations[hip]
        assert item["star_id"] == by_hip[hip][0]
        assert item["value"].hex() == by_hip[hip][10].hex()
        assert item["band"] == "johnson_v"
        assert item["catalogue_field"] == "hipparcos_main_vmag"
        assert item["statistic"] == "unspecified" and item["epoch"] is None
        assert item["variability"] == "unknown"
        assert item["source"]["input_value"].hex() == source_value.hex()
        assert item["source"]["input_sha256"] == (
            "7fa20027153e894b10c45dc84d57b4039d009f587a2be6bea99e3e89e7c256ca")
        assert item["source"]["url"] == "https://cdsarc.cds.unistra.fr/ftp/cats/I/239/ReadMe"
        assert round(source_value, item["source"]["round_decimal_places"]).hex() == item["value"].hex()
    assert annotations[65474]["historical_unqualified_value"] == 1.04
    assert annotations[80763]["historical_unqualified_value"] is None
    # These digests bind the historical extraction's exact row words and
    # ordered lookup data. Undo only the explicitly documented Spica correction.
    historical = copy.deepcopy(catalog["rows"])
    next(row for row in historical if row[3] == 65474)[10] = 1.04
    assert hashlib.sha256(json.dumps(historical, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest() == (
        "b755a21221697982fd2e5137d3273fce56f6aa393cf2ec1e917ed552069593d7")
    keys = ("aliases", "name_to_hip", "hip_fix", "greek", "constellations", "decimal_ranges", "unicode_version")
    lookup = {key: catalog[key] for key in keys}
    assert hashlib.sha256(json.dumps(lookup, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest() == (
        "1ae58868663b815d03669fde26f3e4de3f837adf00e0b444b9fa5a2dadeee069")


@pytest.mark.parametrize("invalid", [
    "duplicate", "identity", "value", "band", "statistic", "epoch", "variability",
    "source-hash", "source-value", "historical",
])
def test_photometry_annotations_refuse_mismatched_or_inferred_qualifications(tmp_path, invalid):
    catalog, _ = generator.load_catalog(generator.DEFAULT_INPUT)
    item = catalog["photometry"][0]
    if invalid == "duplicate":
        catalog["photometry"].append(copy.deepcopy(item))
    elif invalid == "identity":
        item["hip_number"] += 1
    elif invalid == "value":
        item["value"] = 1.04
    elif invalid == "band":
        item["band"] = "hipparcos_hp"
    elif invalid == "statistic":
        item["statistic"] = "median"
    elif invalid == "epoch":
        item["epoch"] = "J1991.25"
    elif invalid == "variability":
        item["variability"] = "constant"
    elif invalid == "source-hash":
        item["source"]["input_sha256"] = "not a hash"
    elif invalid == "source-value":
        item["source"]["input_value"] = 1.04
    elif invalid == "historical":
        item["historical_unqualified_value"] = float("inf")
    source, output = tmp_path / "invalid.json", tmp_path / "preserved.rs"
    source.write_text(json.dumps(catalog))
    output.write_bytes(b"previous output")
    assert generator.main(["--input", str(source), "--out", str(output)]) == 2
    assert output.read_bytes() == b"previous output"


def test_unannotated_catalogs_remain_explicitly_unqualified(tmp_path):
    catalog, _ = generator.load_catalog(generator.DEFAULT_INPUT)
    del catalog["photometry"]
    source = tmp_path / "unqualified.json"
    source.write_text(json.dumps(catalog))
    decoded, digest = generator.load_catalog(source)
    text = generator.render_catalog(decoded, digest)
    assert "const STAR_PHOTOMETRY: &[(i32, StarPhotometry)] = &[\n];" in text
