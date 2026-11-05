"""The shipped Rust catalogue must match the authenticated public data-v4 release."""

import hashlib
import json
import math
from pathlib import Path
import struct

import libephemeris as eph
from libephemeris.download import DATA_FILES, LEB_RELEASES, PLANET_CENTER_RELEASES
import pytest

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "core/data/leb2/data-v4-manifest.json"
MANIFEST_SHA256 = "eadf1997ff9c3b373dc670d08812102602edbdfc7772e5c707910c010478e2bb"


def release_manifest():
    """Authenticate the unchanged published manifest, not a candidate-generated pin."""
    content = MANIFEST.read_bytes()
    assert hashlib.sha256(content).hexdigest() == MANIFEST_SHA256
    manifest = json.loads(content)
    assert manifest["release"] == "data-v4"
    assert len(manifest["files"]) == 12
    return manifest


def test_all_shipped_leb_downloads_match_the_public_data_v4_manifest():
    manifest = release_manifest()
    assert LEB_RELEASES.endswith("/data-v4")
    assert PLANET_CENTER_RELEASES.endswith("/data-v2rc1")
    assert {name for name in DATA_FILES if name.endswith(".leb2")} == set(manifest["files"])
    for name, expected in manifest["files"].items():
        actual = DATA_FILES[name]
        assert actual["sha256"] == expected["sha256"]
        assert actual["dest_subdir"] == "leb"
        assert actual["size_mb"] == round(expected["size"] / 1_000_000, 2)
        if name == "base_core.leb2":
            # The native binding omits URL for assets supplied by the bundle.
            assert "url" not in actual
            assert actual["bundled_resource"] == "native-embedded/base_core.leb2"
        else:
            assert actual["url"] == f"{LEB_RELEASES}/{name}"


def test_bundled_bytes_match_data_v4_identity_and_expanded_header_interval():
    expected = release_manifest()["files"]["base_core.leb2"]
    content = (ROOT / "core/data/leb2/base_core.leb2").read_bytes()
    assert len(content) == expected["size"]
    assert hashlib.sha256(content).hexdigest() == expected["sha256"]
    # The container interval is metadata, not a blanket per-body accuracy claim.
    start, end = struct.unpack_from("<dd", content, 16)
    assert (start, end) == tuple(expected["jd_range"]) == (2396757.5, 2506332.5)


def test_native_reader_serves_both_new_base_boundary_days(tmp_path, monkeypatch):
    """Use TT and TRUEPOS to avoid requiring a neighboring retarded epoch."""
    monkeypatch.setenv("LIBEPHEMERIS_DATA_DIR", str(tmp_path))
    monkeypatch.setenv("LIBEPHEMERIS_CONFIG", str(tmp_path / "absent.toml"))
    monkeypatch.setenv("LIBEPHEMERIS_ENV_FILE", str(tmp_path / "absent.env"))
    eph.reload_configuration()
    eph.set_calc_mode("leb")
    eph.set_network_policy("sealed")
    eph.set_iers_delta_t_enabled(False)
    context = eph.EphemerisContext()
    context.set_leb_file(str(ROOT / "core/data/leb2/base_core.leb2"))
    start, end = release_manifest()["files"]["base_core.leb2"]["jd_range"]
    flags = eph.FLG_TRUEPOS | eph.FLG_BARYCTR | eph.FLG_XYZ | eph.FLG_J2000 | eph.FLG_SPEED
    for epoch in (start, start + 0.5, end - 0.5, end):
        values, _ = context.calc(epoch, eph.SUN, flags)
        assert len(values) == 6 and all(math.isfinite(value) for value in values)
    for epoch in (math.nextafter(start, -math.inf), math.nextafter(end, math.inf)):
        with pytest.raises(eph.EphemerisRangeError):
            context.calc(epoch, eph.SUN, flags)
