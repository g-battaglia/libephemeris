"""
Tests for the LEB generator script (generate_leb.py).

Tests that the generator produces valid .leb files and that round-trip
(generate → read → evaluate) produces correct results.
"""

from __future__ import annotations

import os
import sys
from argparse import Namespace
from types import SimpleNamespace

import numpy as np
import pytest

# Ensure scripts directory is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts"))

from libephemeris.constants import EARTH, MOON, SUN
from libephemeris.leb_format import (
    BODY_PARAMS,
    MAGIC,
)
from libephemeris.leb_reader import LEBReader
from libephemeris.time_utils import julday
from scripts import generate_leb


def _tier_args(**overrides: object) -> Namespace:
    values: dict[str, object] = {
        "tier": "extended",
        "start": None,
        "end": None,
        "start_jd": None,
        "end_jd": None,
        "output": "extended.leb",
    }
    values.update(overrides)
    return Namespace(**values)


def test_extended_tier_defaults_to_exact_de441_boundaries(monkeypatch) -> None:
    monkeypatch.setattr("libephemeris.set_jpl_file", lambda _path: None)

    jd_start, jd_end, _ = generate_leb._resolve_tier(_tier_args())

    assert jd_start == generate_leb.DE441_START_JD
    assert jd_end == generate_leb.DE441_END_JD


def test_explicit_extended_year_override_keeps_calendar_semantics(monkeypatch) -> None:
    monkeypatch.setattr("libephemeris.set_jpl_file", lambda _path: None)

    jd_start, jd_end, _ = generate_leb._resolve_tier(_tier_args(start=-5000, end=5000))

    assert jd_start == generate_leb._year_to_jd(-5000)
    assert jd_end == generate_leb._year_to_jd(5000)


def test_extended_nbody_coverage_is_clamped_to_assist_de441() -> None:
    requested_start = generate_leb._year_to_jd(-13200)
    requested_end = generate_leb._year_to_jd(17191)

    actual_start, actual_end = generate_leb._nbody_coverage_for_range(
        requested_start, requested_end
    )

    assert actual_start == generate_leb._year_to_jd(-13000)
    assert actual_end == generate_leb._year_to_jd(17000)


def test_nbody_coverage_rejects_non_overlapping_range() -> None:
    with pytest.raises(ValueError, match="does not overlap"):
        generate_leb._nbody_coverage_for_range(
            generate_leb._year_to_jd(17001),
            generate_leb._year_to_jd(17191),
        )


def test_required_fit_end_includes_the_complete_final_segment() -> None:
    assert generate_leb._required_fit_end(100.0, 121.0, 8.0) == 124.0
    assert generate_leb._required_fit_end(100.0, 124.0, 8.0) == 124.0


def test_partial_source_range_is_aligned_to_complete_segments() -> None:
    actual_start, actual_end = generate_leb._align_body_range_to_source(
        100.0,
        200.0,
        109.0,
        190.0,
        8.0,
        margin_days=1.0,
    )

    assert actual_start == 110.0
    assert actual_end == 182.0
    assert (actual_end - actual_start) % 8.0 == 0.0


def test_source_backed_grid_retains_exact_kernel_range() -> None:
    start = generate_leb.DE441_START_JD
    end = generate_leb.DE441_END_JD

    actual_start, actual_end, interval = generate_leb._source_backed_fit_grid(
        start,
        end,
        start,
        end,
        64.0,
    )

    assert actual_start == start
    assert actual_end == end
    assert interval <= 64.0
    segment_count = int(np.ceil((end - start) / interval))
    assert abs((start + segment_count * interval) - end) < 1e-6


def test_source_backed_grid_narrows_request_that_exceeds_source() -> None:
    actual_start, actual_end, interval = generate_leb._source_backed_fit_grid(
        100.0,
        201.0,
        99.0,
        200.0,
        8.0,
    )

    assert (actual_start, actual_end, interval) == (100.0, 196.0, 8.0)


def test_source_backed_grid_keeps_every_fit_node_inside_source() -> None:
    start = 100.0
    end = 201.0
    actual_start, actual_end, interval = generate_leb._source_backed_fit_grid(
        start,
        end,
        start,
        end,
        8.0,
    )

    all_jds, segment_count, _ = generate_leb._compute_all_segment_jds(
        actual_start,
        actual_end,
        interval,
        degree=12,
    )

    assert segment_count == 13
    assert float(np.min(all_jds)) >= start
    assert float(np.max(all_jds)) < end


def test_skyfield_vector_evaluator_accepts_source_backed_nodes() -> None:
    class Timescale:
        def __init__(self) -> None:
            self.received: np.ndarray | None = None

        def tt_jd(self, jds: np.ndarray) -> np.ndarray:
            self.received = np.asarray(jds)
            return self.received

    class Target:
        def at(self, jds: np.ndarray) -> SimpleNamespace:
            positions = np.vstack((jds, jds + 1.0, jds + 2.0))
            return SimpleNamespace(position=SimpleNamespace(au=positions))

    jds = np.array([100.25, 101.5, 199.75])
    timescale = Timescale()

    values = generate_leb._eval_target_vectorized(
        Target(),
        jds,
        timescale,
        100.0,
        200.0,
    )

    assert np.array_equal(timescale.received, jds)
    assert values.shape == (3, 3)
    assert np.array_equal(values[:, 0], jds)


@pytest.mark.parametrize(
    "jds",
    [
        np.array([99.999, 100.5]),
        np.array([199.5, 200.0]),
    ],
)
def test_skyfield_vector_evaluator_rejects_nodes_outside_source(
    jds: np.ndarray,
) -> None:
    with pytest.raises(ValueError, match="exceed the selected JPL source"):
        generate_leb._eval_target_vectorized(
            object(),
            jds,
            object(),
            100.0,
            200.0,
        )


def test_spk_coverage_gate_has_no_calendar_tolerance(monkeypatch) -> None:
    from libephemeris.vendor.spktype21 import SPKType21

    target = generate_leb._ASTEROID_NAIF[15]

    class Segment:
        center = 10

        def __init__(self, start_jd: float, end_jd: float) -> None:
            self.target = target
            self.start_second = (start_jd - 2451545.0) * 86400.0
            self.end_second = (end_jd - 2451545.0) * 86400.0

    class Kernel:
        def __init__(self, end_jd: float) -> None:
            self.segments = [Segment(100.0, end_jd)]

        def close(self) -> None:
            pass

    kernel = Kernel(199.0)
    monkeypatch.setattr(SPKType21, "open", lambda _path: kernel)

    assert generate_leb._spk_covers_range("source.bsp", 15, 100.0, 199.0)
    assert not generate_leb._spk_covers_range("source.bsp", 15, 100.0, 200.0)


class TestAssembleLeb:
    """Test the full LEB assembly pipeline."""

    @pytest.mark.integration
    def test_generates_valid_file(self, test_leb_file):
        """Generated .leb file should exist and be non-empty."""
        assert os.path.exists(test_leb_file)
        assert os.path.getsize(test_leb_file) > 0

    @pytest.mark.integration
    def test_file_has_correct_magic(self, test_leb_file):
        """Generated file should start with LEB1 magic bytes."""
        with open(test_leb_file, "rb") as f:
            magic = f.read(4)
        assert magic == MAGIC

    @pytest.mark.integration
    def test_reader_opens_generated_file(self, test_leb_file):
        """LEBReader should successfully parse the generated file."""
        reader = LEBReader(test_leb_file)
        assert reader is not None
        jd_start, jd_end = reader.jd_range
        assert jd_start < jd_end
        reader.close()

    @pytest.mark.integration
    def test_generated_file_has_expected_bodies(self, test_leb_file):
        """Generated file should contain the requested bodies."""
        with LEBReader(test_leb_file) as reader:
            # conftest requests: SUN(0), MOON(1), MARS(4),
            # EARTH(14), MEAN_NODE(10)
            assert reader.has_body(SUN)
            assert reader.has_body(MOON)
            assert reader.has_body(EARTH)

    @pytest.mark.integration
    def test_generated_file_has_nutation(self, test_leb_file):
        """Generated file should include nutation data."""
        with LEBReader(test_leb_file) as reader:
            assert reader._nutation is not None

    @pytest.mark.integration
    def test_generated_file_has_delta_t(self, test_leb_file):
        """Generated file should include Delta-T data."""
        with LEBReader(test_leb_file) as reader:
            assert len(reader._delta_t_jds) > 0

    @pytest.mark.integration
    def test_body_coord_types_correct(self, test_leb_file):
        """Body coordinate types should match BODY_PARAMS."""
        with LEBReader(test_leb_file) as reader:
            for body_id, body in reader._bodies.items():
                if body_id in BODY_PARAMS:
                    expected_coord = BODY_PARAMS[body_id][2]
                    assert body.coord_type == expected_coord, (
                        f"Body {body_id}: coord_type={body.coord_type}, "
                        f"expected={expected_coord}"
                    )


class TestGenerateMinimal:
    """Test minimal file generation (Sun + Earth only, 1 year)."""

    @pytest.mark.integration
    def test_minimal_file_works(self, test_leb_file_minimal):
        """Minimal .leb file should work correctly."""
        with LEBReader(test_leb_file_minimal) as reader:
            assert reader.has_body(SUN)
            assert reader.has_body(EARTH)
            # Should not have other bodies
            assert not reader.has_body(MOON)

    @pytest.mark.integration
    def test_minimal_roundtrip(self, test_leb_file_minimal):
        """Sun position should survive generate→read roundtrip."""
        with LEBReader(test_leb_file_minimal) as reader:
            jd_start, jd_end = reader.jd_range
            jd_mid = (jd_start + jd_end) / 2.0

            pos, vel = reader.eval_body(SUN, jd_mid)
            # ICRS barycentric: Sun is near SSB, typically <0.01 AU
            import math

            dist = math.sqrt(pos[0] ** 2 + pos[1] ** 2 + pos[2] ** 2)
            assert dist < 0.02, f"Sun barycentric distance = {dist} AU (too large)"


class TestGenerateSingleBody:
    """Test individual body generation."""

    @pytest.mark.integration
    @pytest.mark.slow
    def test_generate_sun(self, tmp_path):
        """Generate Sun-only LEB file and verify."""
        from scripts.generate_leb import assemble_leb

        path = str(tmp_path / "sun_only.leb")
        jd_start = julday(2024, 1, 1, 0.0)
        jd_end = julday(2024, 7, 1, 0.0)

        assemble_leb(
            output=path,
            jd_start=jd_start,
            jd_end=jd_end,
            bodies=[0],  # Sun only
            workers=1,
            verbose=False,
        )

        with LEBReader(path) as reader:
            assert reader.has_body(SUN)
            jd_mid = (jd_start + jd_end) / 2.0
            pos, vel = reader.eval_body(SUN, jd_mid)
            assert len(pos) == 3
            assert len(vel) == 3


class TestChebyshevFitting:
    """Test Chebyshev polynomial fitting quality."""

    @pytest.mark.integration
    def test_sun_fit_accuracy(self, test_leb_file):
        """Sun through full pipeline should match Skyfield within 5 arcsec.

        Note: Tolerance is generous (5") for on-the-fly test files with default
        segment parameters.  Production accuracy validated by compare/ tests.
        """
        import libephemeris as ephem
        from libephemeris.constants import FLG_SPEED
        from libephemeris.fast_calc import fast_calc_ut

        with LEBReader(test_leb_file) as reader:
            jd_start, jd_end = reader.jd_range
            jd_mid = (jd_start + jd_end) / 2.0

            # Full pipeline: LEB → fast_calc → geocentric ecliptic
            fast_result, _ = fast_calc_ut(reader, jd_mid, SUN, FLG_SPEED)

            # Skyfield reference (geocentric ecliptic of date)
            ref, _ = ephem.calc_ut(jd_mid, SUN, FLG_SPEED)

            lon_err = abs(fast_result[0] - ref[0])
            if lon_err > 180.0:
                lon_err = 360.0 - lon_err
            lon_err_arcsec = lon_err * 3600.0

            lat_err_arcsec = abs(fast_result[1] - ref[1]) * 3600.0

            max_err_arcsec = max(lon_err_arcsec, lat_err_arcsec)
            assert max_err_arcsec < 5.0, f"Sun fit error = {max_err_arcsec:.4f} arcsec"

    @pytest.mark.integration
    def test_moon_fit_accuracy(self, test_leb_file):
        """Moon through full pipeline should match Skyfield within 5 arcsec.

        Note: Tolerance is generous (5") for on-the-fly test files with default
        segment parameters.  Production accuracy validated by compare/ tests.
        """
        import libephemeris as ephem
        from libephemeris.constants import FLG_SPEED
        from libephemeris.fast_calc import fast_calc_ut

        with LEBReader(test_leb_file) as reader:
            jd_start, jd_end = reader.jd_range
            jd_mid = (jd_start + jd_end) / 2.0

            # Full pipeline: LEB → fast_calc → geocentric ecliptic
            fast_result, _ = fast_calc_ut(reader, jd_mid, MOON, FLG_SPEED)

            # Skyfield reference (geocentric ecliptic of date)
            ref, _ = ephem.calc_ut(jd_mid, MOON, FLG_SPEED)

            lon_err = abs(fast_result[0] - ref[0])
            if lon_err > 180.0:
                lon_err = 360.0 - lon_err
            lon_err_arcsec = lon_err * 3600.0

            lat_err_arcsec = abs(fast_result[1] - ref[1]) * 3600.0

            max_err_arcsec = max(lon_err_arcsec, lat_err_arcsec)
            assert max_err_arcsec < 5.0, f"Moon fit error = {max_err_arcsec:.4f} arcsec"


class TestGeneratorValidation:
    """Malformed generator inputs and verification results fail closed."""

    def test_fit_segment_rejects_non_finite_and_wrong_shape(self) -> None:
        with pytest.raises(ValueError, match="NaN or infinite"):
            generate_leb.fit_segment(
                lambda _jd: np.array([float("nan"), 0.0, 0.0]),
                0.0,
                1.0,
                2,
                3,
            )
        with pytest.raises(ValueError, match="shape"):
            generate_leb.fit_segment(lambda _jd: np.array([0.0, 0.0]), 0.0, 1.0, 2, 3)

    def test_segment_verifiers_reject_zero_samples_and_nan(self) -> None:
        coeffs = np.zeros((3, 2))

        def finite(_jd: float) -> np.ndarray:
            return np.zeros(3)

        def non_finite(_jd: float) -> np.ndarray:
            return np.array([0.0, float("nan"), 0.0])

        with pytest.raises(ValueError, match="n_test must be at least 1"):
            generate_leb.verify_segment(finite, coeffs, 0.0, 1.0, 3, n_test=0)
        with pytest.raises(ValueError, match="NaN or infinite"):
            generate_leb.verify_segment(non_finite, coeffs, 0.0, 1.0, 3, n_test=1)
        with pytest.raises(ValueError, match="n_test must be at least 1"):
            generate_leb._verify_segment_unwrapped(
                finite, coeffs, 0.0, 1.0, 3, n_test=0
            )

    def test_vectorized_verifier_rejects_zero_samples_and_nan(self) -> None:
        values = np.zeros((3, 3))
        with pytest.raises(ValueError, match="n_verify must be at least 1"):
            generate_leb._fit_and_verify_from_values(
                values, 0.0, 1.0, 1.0, 1, 3, 1, 3, n_verify=0
            )

        values[2, 0] = float("nan")
        with pytest.raises(ValueError, match="NaN or infinite"):
            generate_leb._fit_and_verify_from_values(
                values, 0.0, 1.0, 1.0, 1, 3, 1, 3, n_verify=1
            )

    def test_ecliptic_verifier_compares_distance_and_rejects_nan(self) -> None:
        angular, distance = generate_leb._verify_ecliptic_body(
            lambda _jd: np.array([359.99999999, 2.0, 1.25]),
            0.0,
            (0.00000001, 2.0, 1.5),
        )
        assert angular == pytest.approx(2e-8)
        assert distance == pytest.approx(0.25)

        with pytest.raises(ValueError, match="NaN or infinite"):
            generate_leb._verify_ecliptic_body(
                lambda _jd: np.array([0.0, 0.0, float("nan")]),
                0.0,
                (0.0, 0.0, 1.0),
            )

    def test_post_verifier_rejects_zero_samples(self) -> None:
        with pytest.raises(ValueError, match="n_samples must be at least 1"):
            generate_leb.verify_leb("unused.leb", n_samples=0, verbose=False)

    def test_post_verifier_opens_offline_jpl_source_boundary(self, monkeypatch) -> None:
        import libephemeris
        from libephemeris import leb_reader, state

        class Body:
            coord_type = generate_leb.COORD_ICRS_BARY
            jd_start = 0.0
            jd_end = 10.0

        class Reader:
            _bodies = {SUN: Body()}
            jd_range = (0.0, 10.0)

            def __init__(self, _path: str) -> None:
                pass

            def eval_body(self, _body_id: int, _jd: float):
                return (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)

            def close(self) -> None:
                pass

        def verify_source_boundary(_body_id: int, _jd: float, _position):
            assert state._JPL_SOURCE_ACCESS.get() is True
            return 0.0, 1.0

        monkeypatch.setattr(leb_reader, "LEBReader", Reader)
        monkeypatch.setattr(libephemeris, "set_auto_spk_download", lambda _on: None)
        monkeypatch.setattr(generate_leb, "_build_ecliptic_eval_funcs", lambda: {})
        monkeypatch.setattr(generate_leb, "_build_helio_eval_funcs", lambda: {})
        monkeypatch.setattr(generate_leb, "_verify_icrs_planet", verify_source_boundary)

        assert generate_leb.verify_leb("synthetic.leb", n_samples=1, verbose=False)
        assert state._JPL_SOURCE_ACCESS.get() is False

    def test_post_verifier_rejects_non_finite_reader_output(self, monkeypatch) -> None:
        import libephemeris
        from libephemeris import leb_reader

        class Body:
            coord_type = generate_leb.COORD_ICRS_BARY
            jd_start = 0.0
            jd_end = 10.0

        class Reader:
            _bodies = {SUN: Body()}
            jd_range = (0.0, 10.0)

            def __init__(self, _path: str) -> None:
                pass

            def eval_body(self, _body_id: int, _jd: float):
                return (float("nan"), 0.0, 0.0), (0.0, 0.0, 0.0)

            def close(self) -> None:
                pass

        monkeypatch.setattr(leb_reader, "LEBReader", Reader)
        monkeypatch.setattr(libephemeris, "set_auto_spk_download", lambda _on: None)
        monkeypatch.setattr(generate_leb, "_build_ecliptic_eval_funcs", lambda: {})
        monkeypatch.setattr(generate_leb, "_build_helio_eval_funcs", lambda: {})

        assert not generate_leb.verify_leb("synthetic.leb", n_samples=1, verbose=False)

    def test_post_verifier_always_samples_body_edges(self, monkeypatch) -> None:
        import libephemeris
        from libephemeris import leb_reader

        class Body:
            coord_type = generate_leb.COORD_ICRS_BARY
            jd_start = 0.0
            jd_end = 10.0
            interval_days = 2.0
            degree = 3

        class Reader:
            _bodies = {SUN: Body()}
            jd_range = (0.0, 10.0)

            def __init__(self, _path: str) -> None:
                pass

            def eval_body(self, _body_id: int, jd: float):
                edge_error = 0.1 if jd >= 9.0 else 0.0
                return (edge_error, 0.0, 0.0), (0.0, 0.0, 0.0)

            def close(self) -> None:
                pass

        def compare(_body_id: int, _jd: float, position):
            return abs(float(position[0])), 1.0

        monkeypatch.setattr(leb_reader, "LEBReader", Reader)
        monkeypatch.setattr(libephemeris, "set_auto_spk_download", lambda _on: None)
        monkeypatch.setattr(generate_leb, "_build_ecliptic_eval_funcs", lambda: {})
        monkeypatch.setattr(generate_leb, "_build_helio_eval_funcs", lambda: {})
        monkeypatch.setattr(generate_leb, "_verify_icrs_planet", compare)

        assert not generate_leb.verify_leb("synthetic.leb", n_samples=1, verbose=False)

    def test_post_verifier_rejects_ecliptic_distance_error(self, monkeypatch) -> None:
        import libephemeris
        from libephemeris import leb_reader

        body_id = 10

        class Body:
            coord_type = generate_leb.COORD_ECLIPTIC
            jd_start = 0.0
            jd_end = 10.0

        class Reader:
            _bodies = {body_id: Body()}
            jd_range = (0.0, 10.0)

            def __init__(self, _path: str) -> None:
                pass

            def eval_body(self, _body_id: int, _jd: float):
                return (0.0, 0.0, 1.0), (0.0, 0.0, 0.0)

            def close(self) -> None:
                pass

        monkeypatch.setattr(leb_reader, "LEBReader", Reader)
        monkeypatch.setattr(libephemeris, "set_auto_spk_download", lambda _on: None)
        monkeypatch.setattr(
            generate_leb,
            "_build_ecliptic_eval_funcs",
            lambda: {body_id: lambda _jd: np.zeros(3)},
        )
        monkeypatch.setattr(generate_leb, "_build_helio_eval_funcs", lambda: {})

        assert not generate_leb.verify_leb("synthetic.leb", n_samples=1, verbose=False)
