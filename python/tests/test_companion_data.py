"""A selected group core must serve its actual companion coefficients."""
from data_paths import CORE_LEB, DATA_ROOT, LEB_HASHES
from pathlib import Path
import hashlib

import libephemeris as eph
import pytest

ROOT = DATA_ROOT / "leb"
HASHES = {f"base_{group}.leb2": digest for group, digest in LEB_HASHES['base'].items()}
# Actual pinned v3 results from the matched-data configurations grid.
REFERENCE = [
    (17, (184.4529148352954, 11.83754537505261, 2.256828151298012,
          0.2196661753259832, 0.07602292813047029, -0.01265917951265294)),
    (146199, (18.59127451963572, -15.603568099978439, 97.02360622901782,
              -0.0009858513294602744, 0.003259703746039122, 0.016558680769662715)),
]


@pytest.mark.parametrize("body,reference", REFERENCE)
def test_explicit_core_keeps_companion_precision(body, reference):
    for name, digest in HASHES.items():
        path = ROOT / name
        if not path.is_file():
            pytest.skip("optional local coefficient partitions are unavailable")
        with path.open("rb") as stream:
            assert hashlib.file_digest(stream, "sha256").hexdigest() == digest
    eph.set_calc_mode("leb")
    eph.set_precision_tier("base")
    eph.set_iers_delta_t_enabled(False)
    eph.set_delta_t_userdef(64. / 86400.)
    context = eph.EphemerisContext()
    context.set_leb_file(str(ROOT / "base_core.leb2"))
    position, flags = context.calc(2451545., body, 258)
    assert flags == 258
    for slot in (0, 1):
        assert abs(position[slot] - reference[slot]) * 3600. <= .001
    assert position[2] == pytest.approx(reference[2], rel=1e-8, abs=0.)
    for slot in (3, 4, 5):
        assert position[slot] == pytest.approx(reference[slot], rel=1e-6, abs=0.)


def test_curated_minor_without_coefficients_uses_the_local_model():
    path = ROOT / "base_core.leb2"
    if not path.is_file():
        pytest.skip("optional local coefficient partitions are unavailable")
    with path.open("rb") as stream:
        assert hashlib.file_digest(stream, "sha256").hexdigest() == HASHES[path.name]
    eph.set_calc_mode("leb")
    eph.set_precision_tier("base")
    eph.set_iers_delta_t_enabled(False)
    eph.set_delta_t_userdef(64. / 86400.)
    context = eph.EphemerisContext()
    context.set_leb_file(str(path))
    actual, flags = context.calc_ut(2451545., 111955, 258)
    reference = (215.98064696360825, 32.78973871336137, 0.17430990629822682,
                 -0.11827734097096254, -0.41330139024466916, 0.0013222793959766221)
    assert flags == 258
    for slot in (0, 1):
        assert abs(actual[slot] - reference[slot]) * 3600. <= .001
    assert actual[2] == pytest.approx(reference[2], rel=1e-8, abs=0.)
    for slot in (3, 4, 5):
        assert actual[slot] == pytest.approx(reference[slot], rel=1e-6, abs=0.)


def test_forced_leb_ignores_available_centers_while_auto_uses_them(tmp_path, monkeypatch):
    """Provider mode controls SPK offsets over identical stored barycenters."""
    coefficients = CORE_LEB
    centers = DATA_ROOT / "planet_centers_base.bsp"
    if not centers.is_file():
        pytest.skip("optional local planet-center SPK is unavailable")
    for path, digest in (
        (coefficients, LEB_HASHES["base"]["core"]),
        (centers, "a9ec744ff412b095129166587ea0814f81c850faebf92586a738cb5dc103c92a"),
    ):
        with path.open("rb") as stream:
            assert hashlib.file_digest(stream, "sha256").hexdigest() == digest
    (tmp_path / centers.name).symlink_to(centers)
    monkeypatch.setenv("LIBEPHEMERIS_DATA_DIR", str(tmp_path))
    monkeypatch.setenv("LIBEPHEMERIS_CONFIG", str(tmp_path / "absent.toml"))
    monkeypatch.setenv("LIBEPHEMERIS_ENV_FILE", str(tmp_path / "absent.env"))
    eph.set_precision_tier("base")
    eph.set_iers_delta_t_enabled(False)
    eph.set_delta_t_userdef(64.0 / 86400.0)
    eph.set_auto_spk_download(False)
    # Actual pinned Source calc(2433282.5, JUPITER, SPEED|TRUEPOS|SWIEPH).
    # The authenticated provider trace uses the selected inputs. The LEB row is
    # COORD_ICRS_BARY_SYSTEM=4; Source forbids its optional SPK offset in leb.
    references = {
        "leb": (306.5129995420849, -0.58403428789559, 5.935331520171296,
                0.225227884754986, -0.0007402404991196931, 0.006796374637296765),
        "auto": (306.5129896304244, -0.5840345045021424, 5.935330689412087,
                 0.22522530838386956, -0.0007404002099198337, 0.006797397190894117),
    }
    positions = {}
    try:
        for mode, reference in references.items():
            eph.set_calc_mode(mode)
            context = eph.EphemerisContext(ephe_path=str(tmp_path))
            context.set_leb_file(str(coefficients))
            actual, flags = context.calc(2433282.5, eph.JUPITER, 274)
            positions[mode] = actual
            assert flags == 1810
            for slot in (0, 1):
                assert abs(actual[slot] - reference[slot]) * 3600.0 <= 0.001
            assert actual[2] == pytest.approx(reference[2], rel=1e-8, abs=0.0)
            for slot in (3, 4, 5):
                assert actual[slot] == pytest.approx(reference[slot], rel=1e-6, abs=0.0)
        assert abs(positions["auto"][0] - positions["leb"][0]) * 3600.0 > 0.001
    finally:
        eph.set_calc_mode("leb")
