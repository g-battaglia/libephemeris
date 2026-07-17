"""Focused tests for the LEB-backed Skyfield vector adapter."""

from __future__ import annotations

import numpy as np
import pytest

from libephemeris.constants import EARTH, MOON
from libephemeris.exceptions import EphemerisRangeError, LEBCorruptionError
from libephemeris.leb_vector import LEBVectorEphemeris
from libephemeris.state import get_timescale
from libephemeris.time_utils import julday


def test_scalar_and_array_vector_shapes(leb_reader):
    vectors = LEBVectorEphemeris(leb_reader)
    ts = get_timescale()
    jd = julday(2024, 6, 15, 12.0)

    scalar = vectors["earth"].at(ts.tt_jd(jd))
    array = vectors["earth"].at(ts.tt_jd([jd, jd + 1.0]))

    assert scalar.position.au.shape == (3,)
    assert scalar.velocity.au_per_d.shape == (3,)
    assert array.position.au.shape == (3, 2)
    assert array.velocity.au_per_d.shape == (3, 2)
    assert np.all(np.isfinite(array.position.au))


def test_aliases_expose_system_barycenters_not_planet_center_ids(leb_reader):
    vectors = LEBVectorEphemeris(leb_reader)

    assert vectors["jupiter"] is vectors["jupiter barycenter"]
    assert vectors[5] is vectors["jupiter"]
    assert 599 not in vectors
    assert vectors["earth barycenter"].target == 3


def test_earth_moon_barycenter_is_mass_weighted(leb_reader):
    vectors = LEBVectorEphemeris(leb_reader)
    jd = julday(2024, 6, 15, 12.0)
    ts = get_timescale()

    earth = np.asarray(leb_reader.eval_body(EARTH, jd)[0])
    moon = np.asarray(leb_reader.eval_body(MOON, jd)[0])
    expected = earth + (moon - earth) / (81.3005691 + 1.0)
    actual = vectors["earth barycenter"].at(ts.tt_jd(jd)).position.au

    assert np.allclose(actual, expected, rtol=0.0, atol=1e-15)


def test_adapter_range_miss_is_public_ephemeris_error(leb_reader):
    vectors = LEBVectorEphemeris(leb_reader)
    ts = get_timescale()
    outside = julday(2030, 1, 1, 0.0)

    with pytest.raises(EphemerisRangeError) as exc_info:
        vectors["earth"].at(ts.tt_jd(outside))

    assert exc_info.value.body_id == EARTH
    assert exc_info.value.start_jd is not None
    assert exc_info.value.end_jd is not None


def test_adapter_never_relabels_corruption_as_a_range_miss():
    class CorruptReader:
        def has_body(self, body_id):
            return body_id == EARTH

        def eval_body(self, body_id, jd_tt):
            raise LEBCorruptionError("corrupt coefficients outside valid bounds")

        def body_coverage(self, body_id):
            return (2450000.0, 2460000.0)

    vectors = LEBVectorEphemeris(CorruptReader())
    ts = get_timescale()

    with pytest.raises(LEBCorruptionError, match="corrupt coefficients"):
        vectors["earth"].at(ts.tt_jd(2440000.0))
