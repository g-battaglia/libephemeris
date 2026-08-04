"""Independent DE440 checks for the generated lunar-apse passage curves."""

from __future__ import annotations

import numpy as np
import pytest

from libephemeris.lunar import calc_interpolated_apogee, calc_interpolated_perigee
from libephemeris.state import get_planets
from scripts.generate_lunar_apse_model import _moon_states, find_passages

J2000 = 2451545.0


@pytest.mark.parametrize(
    ("passage_index", "evaluator", "max_error_arcsec"),
    [
        (0, calc_interpolated_apogee, 6.0),
        (1, calc_interpolated_perigee, 1.0),
    ],
)
def test_generated_longitudes_follow_de440_apsis_passages(
    passage_index: int,
    evaluator: object,
    max_error_arcsec: float,
) -> None:
    """Check a local passage set not used as a compatibility target."""
    passage_sets = find_passages(J2000 - 180.0, J2000 + 180.0)
    passages = passage_sets[passage_index]
    passage_longitudes, _, _ = _moon_states(passages)

    errors = []
    for jd, expected_longitude in zip(passages, passage_longitudes):
        actual_longitude = evaluator(float(jd))[0]  # type: ignore[operator]
        delta = (actual_longitude - float(expected_longitude) + 180.0) % 360.0 - 180.0
        errors.append(abs(delta) * 3600.0)

    assert len(errors) >= 10
    assert float(np.max(errors)) < max_error_arcsec


@pytest.fixture(scope="module")
def terminal_passage_errors() -> tuple[float, float]:
    """Return max longitude errors near both ends of DE440 coverage."""
    spk = get_planets()["moon"].ephemeris.spk
    moon_segments = [segment for segment in spk.segments if segment.target == 301]
    kernel_start = max(segment.start_jd for segment in moon_segments)
    kernel_end = min(segment.end_jd for segment in moon_segments)
    first = find_passages(kernel_start + 3.0, kernel_start + 200.0)
    last = find_passages(kernel_end - 200.0, kernel_end - 3.0)

    maxima = []
    for index, evaluator in enumerate(
        (calc_interpolated_apogee, calc_interpolated_perigee)
    ):
        passages = np.concatenate((first[index], last[index]))
        expected_longitudes, _, _ = _moon_states(passages)
        errors = []
        for jd, expected_longitude in zip(passages, expected_longitudes):
            actual_longitude = evaluator(float(jd))[0]
            delta = (
                actual_longitude - float(expected_longitude) + 180.0
            ) % 360.0 - 180.0
            errors.append(abs(delta) * 3600.0)
        maxima.append(float(np.max(errors)))
    return maxima[0], maxima[1]


def test_lunar_apse_model_matches_pinned_sha256() -> None:
    """The generated module is the exact SHA-256-pinned artifact.

    ``libephemeris/lunar_apse_model.py`` is an immutable generator output
    (see its header): any change must come from rerunning
    ``scripts/generate_lunar_apse_model.py`` against the reviewed DE440
    kernel, followed by re-pinning this hash in the same commit. A mismatch
    here means the artifact was hand-edited or the regeneration was not
    re-pinned.
    """
    import hashlib
    from pathlib import Path

    import libephemeris.lunar_apse_model as model

    digest = hashlib.sha256(Path(model.__file__).read_bytes()).hexdigest()
    assert digest == (
        "0132e75b0b39cd5df540f90acd5ebb740a3f39e593e292ab5c003106ccdc08d9"
    )


def test_apogee_edge_taper_preserves_terminal_passage_precision(
    terminal_passage_errors: tuple[float, float],
) -> None:
    assert terminal_passage_errors[0] < 20.0


def test_perigee_edge_taper_preserves_terminal_passage_precision(
    terminal_passage_errors: tuple[float, float],
) -> None:
    assert terminal_passage_errors[1] < 7.0
