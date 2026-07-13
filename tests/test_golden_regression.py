"""Golden file regression tests.

Compares current libephemeris calculations against a stored reference
(golden file) to detect any unexpected changes in output. Any difference
triggers an explicit test failure requiring review.

Validation Plan v2, Section 7.1.

To regenerate the golden file after an intentional change:
    .venv/bin/python3 scripts/generate_golden.py
"""

from __future__ import annotations

import hashlib
import json
import os
import warnings
from pathlib import Path

import pytest

import libephemeris as swe
from libephemeris.download import DATA_FILES

warnings.filterwarnings("ignore")

GOLDEN_FILE = os.path.join(os.path.dirname(__file__), "golden", "golden_reference.json")


@pytest.fixture(scope="module")
def golden_data() -> dict:
    """Load the golden reference file."""
    if not os.path.exists(GOLDEN_FILE):
        pytest.skip(f"Golden file not found: {GOLDEN_FILE}")
    with open(GOLDEN_FILE) as f:
        return json.load(f)


@pytest.fixture(scope="module")
def golden_entries(golden_data: dict) -> list[dict]:
    """Return the list of golden entries."""
    return golden_data["entries"]


def entries_by_type(entries: list[dict], entry_type: str) -> list[dict]:
    """Filter entries by type."""
    return [e for e in entries if e["type"] == entry_type]


@pytest.fixture
def production_leb_source(monkeypatch):
    """Pin and isolate the exact source used by the golden generator."""
    from libephemeris import state

    reviewed_core = (
        Path(__file__).resolve().parents[1]
        / "libephemeris"
        / "data"
        / "leb2"
        / "base_core.leb2"
    )
    expected_hash = DATA_FILES["base_core.leb2"]["sha256"]
    actual_hash = hashlib.sha256(reviewed_core.read_bytes()).hexdigest()
    assert actual_hash == expected_hash

    previous_mode = state.get_calc_mode()
    previous_file = state._LEB_FILE
    monkeypatch.delenv("LIBEPHEMERIS_LEB", raising=False)
    monkeypatch.setattr(
        state,
        "get_planet_center_segment",
        lambda naif_id, jd=None: None,
    )
    state.set_leb_file(str(reviewed_core))
    state.set_calc_mode("leb")
    try:
        yield
    finally:
        state.set_leb_file(previous_file)
        state.set_calc_mode(previous_mode)


class TestGoldenCalcUt:
    """§7.1 Golden file regression for calc_ut results."""

    def test_calc_ut_count(self, golden_entries: list[dict]) -> None:
        """Verify expected number of calc_ut entries exist."""
        calc_entries = entries_by_type(golden_entries, "calc_ut")
        assert len(calc_entries) >= 96, (
            f"Expected >=96 calc_ut entries, got {len(calc_entries)}"
        )

    def test_calc_ut_positions_match(
        self, golden_entries: list[dict], production_leb_source
    ) -> None:
        """All calc_ut positions must match golden values within tolerance."""

        calc_entries = entries_by_type(golden_entries, "calc_ut")
        mismatches = []

        for entry in calc_entries:
            jd = entry["jd"]
            body = entry["body"]
            flags = entry["flags"]
            expected = entry["result"]

            pos, retflag = swe.calc_ut(jd, body, flags)

            for i in range(6):
                diff = abs(pos[i] - expected[i])
                # Tolerance: 1e-8 degrees for positions (~0.036 mas)
                # This is tight enough to catch real regressions but allows
                # for floating-point platform differences.
                tol = 1e-8
                if diff > tol:
                    mismatches.append(
                        f"body={body} jd={jd} flags={flags} elem={i}: "
                        f"got={pos[i]:.12f} expected={expected[i]:.12f} "
                        f"diff={diff:.2e}"
                    )

        assert not mismatches, (
            f"{len(mismatches)} regression(s) detected:\n" + "\n".join(mismatches[:20])
        )


class TestGoldenHouses:
    """§7.1 Golden file regression for house calculations."""

    def test_houses_match(self, golden_entries: list[dict]) -> None:
        """All house cusp positions must match golden values."""
        house_entries = entries_by_type(golden_entries, "houses")
        mismatches = []

        for entry in house_entries:
            jd = entry["jd"]
            lat = entry["lat"]
            lon = entry["lon"]
            hsys = ord(entry["hsys"])
            expected_cusps = entry["cusps"]
            expected_angles = entry["angles"]

            cusps, angles = swe.houses(jd, lat, lon, hsys)

            for i in range(min(len(cusps), len(expected_cusps))):
                diff = abs(cusps[i] - expected_cusps[i])
                if diff > 1e-8:
                    mismatches.append(
                        f"houses {entry['location']} {entry['hsys']} "
                        f"cusp {i}: diff={diff:.2e}"
                    )

            for i in range(min(len(angles), len(expected_angles))):
                diff = abs(angles[i] - expected_angles[i])
                if diff > 1e-8:
                    mismatches.append(
                        f"houses {entry['location']} {entry['hsys']} "
                        f"angle {i}: diff={diff:.2e}"
                    )

        assert not mismatches, f"{len(mismatches)} house regression(s):\n" + "\n".join(
            mismatches[:20]
        )


class TestGoldenSidereal:
    """§7.1 Golden file regression for sidereal calculations."""

    def test_sidereal_match(
        self, golden_entries: list[dict], production_leb_source
    ) -> None:
        """Sidereal positions must match golden values."""

        sid_entries = entries_by_type(golden_entries, "sidereal")
        mismatches = []

        for entry in sid_entries:
            jd = entry["jd"]
            body = entry["body"]
            mode = entry["sid_mode"]
            expected = entry["result"]

            swe.set_sid_mode(mode)
            pos, _ = swe.calc_ut(jd, body, swe.FLG_SIDEREAL | swe.FLG_SPEED)

            for i in range(6):
                diff = abs(pos[i] - expected[i])
                if diff > 1e-8:
                    mismatches.append(
                        f"sidereal body={body} mode={entry['sid_mode_name']} "
                        f"elem={i}: diff={diff:.2e}"
                    )

        # Reset
        swe.set_sid_mode(swe.SIDM_J2000)

        assert not mismatches, (
            f"{len(mismatches)} sidereal regression(s):\n" + "\n".join(mismatches[:20])
        )


class TestGoldenTime:
    """§7.1 Golden file regression for time conversions."""

    def test_julday_match(self, golden_entries: list[dict]) -> None:
        """Julian day conversions must match golden values."""
        jd_entries = entries_by_type(golden_entries, "julday")

        for entry in jd_entries:
            y, m, d, h = entry["input"]
            expected_jd = entry["jd"]

            jd = swe.julday(y, m, d, h)
            assert abs(jd - expected_jd) < 1e-10, (
                f"julday({y},{m},{d},{h}): got {jd}, expected {expected_jd}"
            )

    def test_sidtime_match(self, golden_entries: list[dict]) -> None:
        """Sidereal time must match golden values."""
        st_entries = entries_by_type(golden_entries, "sidtime")

        for entry in st_entries:
            jd = entry["jd"]
            expected = entry["result"]

            st = swe.sidtime(jd)
            diff = abs(st - expected)
            assert diff < 1e-10, (
                f"sidtime(jd={jd}): got {st}, expected {expected}, diff={diff:.2e}"
            )

    def test_deltat_match(self, golden_entries: list[dict]) -> None:
        """Delta T must match golden values."""
        dt_entries = entries_by_type(golden_entries, "deltat")

        for entry in dt_entries:
            jd = entry["jd"]
            expected = entry["result"]

            dt = swe.deltat(jd)
            diff = abs(dt - expected)
            assert diff < 1e-12, (
                f"deltat(jd={jd}): got {dt}, expected {expected}, diff={diff:.2e}"
            )


class TestGoldenEclipses:
    """§7.1 Golden file regression for eclipse calculations."""

    def test_solar_eclipse_match(self, golden_entries: list[dict]) -> None:
        """Solar eclipse times must match golden values."""
        ecl_entries = entries_by_type(golden_entries, "solar_eclipse")

        for entry in ecl_entries:
            search_jd = entry["search_jd"]
            expected_type = entry["ecl_type"]
            expected_times = entry["times"]

            ecl_type, times = swe.sol_eclipse_when_glob(
                search_jd, ecltype=swe.ECL_TOTAL
            )

            assert ecl_type == expected_type, (
                f"Eclipse type changed: got {ecl_type}, expected {expected_type}"
            )

            for i in range(len(expected_times)):
                if expected_times[i] > 0:
                    diff = abs(times[i] - expected_times[i]) * 86400
                    assert diff < 1.0, f"Eclipse time[{i}] diff {diff:.3f}s"

    def test_lunar_eclipse_match(self, golden_entries: list[dict]) -> None:
        """Lunar eclipse times must match golden values."""
        ecl_entries = entries_by_type(golden_entries, "lunar_eclipse")

        for entry in ecl_entries:
            search_jd = entry["search_jd"]
            expected_type = entry["ecl_type"]
            expected_times = entry["times"]

            ecl_type, times = swe.lun_eclipse_when(search_jd, ecltype=swe.ECL_TOTAL)

            assert ecl_type == expected_type, (
                f"Eclipse type changed: got {ecl_type}, expected {expected_type}"
            )

            for i in range(len(expected_times)):
                if expected_times[i] > 0:
                    diff = abs(times[i] - expected_times[i]) * 86400
                    assert diff < 1.0, f"Eclipse time[{i}] diff {diff:.3f}s"


class TestGoldenFileIntegrity:
    """Verify golden file metadata and structure."""

    def test_golden_file_exists(self) -> None:
        """Golden reference file must exist."""
        assert os.path.exists(GOLDEN_FILE), (
            f"Golden file missing: {GOLDEN_FILE}\n"
            f"Generate with: .venv/bin/python3 scripts/generate_golden.py"
        )

    def test_golden_file_version(self, golden_data: dict) -> None:
        """Golden file must have expected version."""
        assert golden_data["version"] == 1

    def test_golden_file_records_independent_provenance(
        self, golden_data: dict
    ) -> None:
        """The regression baseline must identify its independent generator."""
        assert golden_data["provenance"].endswith("no external-reference calls")
        assert golden_data["source_artifact"] == (
            "libephemeris/data/leb2/base_core.leb2"
        )
        assert golden_data["source_sha256"] == DATA_FILES["base_core.leb2"]["sha256"]
        assert golden_data["optional_planet_centers"] == (
            "disabled; base_core-only regression"
        )

    def test_golden_file_entry_count(self, golden_data: dict) -> None:
        """Golden file must have >= 100 entries."""
        assert golden_data["entry_count"] >= 100, (
            f"Expected >= 100 entries, got {golden_data['entry_count']}"
        )
