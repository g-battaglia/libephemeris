"""Useful module-level APIs execute native physics and publish actual provenance."""

from data_paths import CORE_LEB, DATA_ROOT
from pathlib import Path
import math

import libephemeris as eph
from libephemeris import fixed_stars, hypothetical, lunar, minor_bodies
import pytest


@pytest.fixture(autouse=True)
def configured_native_engine():
    eph.set_calc_mode("leb")
    eph.set_precision_tier("base")
    eph.set_leb_file(str(CORE_LEB))
    eph.set_iers_delta_t_enabled(False)
    eph.set_delta_t_userdef(64.0 / 86400.0)
    eph.set_auto_spk_download(False)
    eph.set_topo(12.5, 41.9, 150.0)
    yield
    eph.set_calc_mode("leb")


def test_stellar_trace_records_success_and_not_catalog_errors():
    token = eph.start_tracing()
    try:
        eph.fixstar_ut("Regulus", 2451545.0, eph.FLG_SPEED)
        assert eph.get_trace_results() == {eph.REGULUS: "LEB"}
        with pytest.raises(eph.StarNotFoundError):
            eph.fixstar_ut("unavailable-star-test", 2451545.0)
        assert eph.get_trace_results() == {eph.REGULUS: "LEB"}
    finally:
        token.reset()


def test_stellar_batch_trace_excludes_skipped_catalog_rows():
    token = eph.start_tracing()
    try:
        eph.batch_fixstars_ut(
            ["Regulus", "unavailable-star-test", "Spica"],
            2451545.0,
            eph.FLG_SPEED,
            skip_errors=True,
        )
        assert eph.get_trace_results() == {eph.REGULUS: "LEB", eph.SPICA_STAR: "LEB"}
    finally:
        token.reset()


def test_stellar_trace_names_the_actual_native_jpl_provider():
    main = DATA_ROOT / "de440.bsp"
    if not main.is_file():
        pytest.skip("the local DE440 image is required")
    eph.set_calc_mode("skyfield")
    eph.set_ephe_path(str(main.parent))
    eph.set_ephemeris_file(main.name)
    token = eph.start_tracing()
    try:
        eph.fixstar_ut("Regulus", 2451545.0, eph.FLG_SPEED)
        assert eph.get_trace_results() == {eph.REGULUS: "NativeJpl"}
        eph.batch_fixstars_ut(
            ["Spica", "unavailable-star-test"],
            2451545.0,
            eph.FLG_SPEED,
            skip_errors=True,
        )
        assert eph.get_trace_results() == {
            eph.REGULUS: "NativeJpl",
            eph.SPICA_STAR: "NativeJpl",
        }
    finally:
        token.reset()


def test_stellar_module_numeric_helpers_use_the_catalog_and_raw_place():
    identifier = fixed_stars.resolve_star_name("Regulus")
    assert identifier == eph.REGULUS
    assert fixed_stars.get_canonical_star_name(identifier) == "Regulus"
    assert (
        fixed_stars.get_hip_from_star_name("Regulus")
        == fixed_stars.STAR_NAME_TO_HIP["REGULUS"]
        == 49669
    )
    assert identifier in fixed_stars.FIXED_STARS
    position = fixed_stars.calc_fixed_star_position(
        identifier, 2451545.0, True, True, True
    )
    state = fixed_stars.calc_fixed_star_velocity(
        identifier, 2451545.0, True, True, True
    )
    assert state[:3] == position
    assert all(math.isfinite(value) for value in state) and state[2] > 1e6


def test_osculating_perigee_uses_opposite_eccentricity_vector_and_its_own_radius():
    apogee = lunar.calc_true_lilith(2451545.0)
    perigee = lunar.calc_osculating_perigee(2451545.0)
    assert abs(math.remainder(perigee[0] - apogee[0], 360.0)) == pytest.approx(
        180.0, abs=0.001 / 3600.0
    )
    assert perigee[1] == pytest.approx(-apogee[1], abs=0.001 / 3600.0)
    assert 0.0 < perigee[2] < apogee[2]
    comparison = lunar.compare_true_lilith_methods(2451545.0)
    assert comparison["eccentricity_vector"] == comparison["orbital_elements"]
    assert (
        comparison["lon_diff"] == comparison["lat_diff"] == comparison["e_diff"] == 0.0
    )
    assert lunar.calc_true_lilith_orbital_elements(2451545.0) == apogee
    for function in (lunar.calc_mean_lunar_node_state, lunar.calc_mean_lilith_state):
        assert len(function(2451545.0, speed3=True)) == 6
    assert len(lunar.calc_interpolated_apse_state(2451545.0, eph.INTP_APOG)) == 6
    assert lunar.MEEUS_OPTIMAL_CENTURIES == 2.0
    assert lunar.MEEUS_MAX_CENTURIES == 20.0


def test_hypothetical_module_helpers_and_constant_aliases():
    assert hypothetical.CUPIDO == eph.CUPIDO
    assert hypothetical.TRANSPLUTO == hypothetical.ISIS == eph.ISIS
    assert hypothetical.list_hypothetical_bodies()[eph.CUPIDO] == "Cupido"
    assert hypothetical.get_hypothetical_name(eph.CUPIDO) == "Cupido"
    assert hypothetical.is_hypothetical_body(eph.CUPIDO)
    assert hypothetical.calc_uranian_position(
        eph.CUPIDO, 2451545.0
    ) == hypothetical.calc_cupido(2451545.0)
    assert (
        hypothetical.calc_uranian_longitude(eph.CUPIDO, 2451545.0)
        == hypothetical.calc_cupido(2451545.0)[0]
    )
    assert hypothetical.calc_transpluto_position(
        2451545.0
    ) == hypothetical.calc_transpluto(2451545.0)
    assert hypothetical.calc_waldemath_position(
        2451545.0
    ) == hypothetical.calc_waldemath(2451545.0)
    assert hypothetical.calc_planet_x_lowell(
        2451545.0
    ) == hypothetical.calc_fictitious_position(eph.PLUTO_LOWELL, 2451545.0)
    assert (
        hypothetical.HYPOTHETICAL_ELEMENTS[eph.ISIS]
        is hypothetical.TRANSPLUTO_KEPLERIAN_ELEMENTS
    )


def test_minor_module_solver_residuals_and_typed_propagation():
    for mean, eccentricity in [(0.3, 0.2), (-1.2, 0.9)]:
        anomaly = minor_bodies.solve_kepler_equation_elliptic(mean, eccentricity)
        assert abs(anomaly - eccentricity * math.sin(anomaly) - mean) < 1e-8
    anomaly = minor_bodies.solve_kepler_equation_hyperbolic(0.3, 1.4)
    assert abs(1.4 * math.sinh(anomaly) - anomaly - 0.3) < 1e-8
    assert math.isfinite(minor_bodies.solve_barker_equation(0.3))
    assert minor_bodies.solve_kepler_equation(
        0.3, 0.2
    ) == minor_bodies.solve_kepler_equation_elliptic(0.3, 0.2)
    elements = minor_bodies.MINOR_BODY_ELEMENTS[eph.CHIRON]
    assert len(minor_bodies.calc_minor_body_position(elements, 2451545.0, False)) == 3
    assert len(minor_bodies.apply_secular_perturbations(elements, elements.epoch)) == 6
    assert len(minor_bodies.calc_secular_perturbation_rates(elements)) == 3
    assert minor_bodies.validate_elements(elements) is None
    assert (
        len(minor_bodies.calc_minor_body_heliocentric(eph.CHIRON, 2451545.0, False))
        == 3
    )
    assert minor_bodies.get_resonance_info(-1) is None
    assert not minor_bodies.is_in_jpl_major_body_index(eph.CERES)
    assert minor_bodies.is_spk_downloadable(eph.CERES)


def test_analytical_moon_module_paths_and_rectangular_units():
    from libephemeris.moon_theories import (
        charon,
        triton,
        galilean,
        uranian,
        tass17,
        constants,
    )

    jd = 2451545.0
    assert 0.0001 < math.dist(charon.charon_position(jd), (0, 0, 0)) < 0.0002
    assert 0.002 < math.dist(triton.triton_position(jd), (0, 0, 0)) < 0.003
    galileans = galilean.galilean_moon_positions(jd)
    uranians = uranian.uranian_moon_positions(jd)
    assert len(galileans) == 4 and len(uranians) == 5
    assert 400_000 < math.dist(galileans[0], (0, 0, 0)) < 450_000
    assert 180_000 < math.dist(uranians[0], (0, 0, 0)) < 200_000
    for body in range(8):
        position = tass17.saturn_moon_position(jd, body)
        state = tass17.saturn_moon_position_velocity(jd, body)
        assert state[:3] == position
        assert all(math.isfinite(value) for value in state)
    assert constants.get_cob_offset("earth", jd) == (0.0, 0.0, 0.0)
    assert len(tass17.SATELLITE_NAMES) == 8
    assert triton.AU_KM == 149597870.7


def test_exotic_module_native_registry_and_records():
    from libephemeris import exotic_bodies as exotic

    assert len(exotic.EXOTIC_BODIES) == 31
    assert len(exotic.EXOTIC_IDS) == 31
    assert exotic.naif_map()[146199] == 2136199
    assert exotic.body_params()[146199] == (32.0, 13)
    assert exotic.name_map()[11181] == "Lilith-ast"
    assert exotic.exotic_display_name(11181) == "Lilith"
    assert exotic.exotic_display_name(-1) == ""
    assert exotic.target_au_map()[109942] == 1e-12
    entry = exotic.EXOTIC_BODIES[146199]
    copy = exotic.ExoticBody(
        entry.body_id,
        entry.naif,
        entry.interval_days,
        entry.degree,
        entry.target_au,
        entry.name,
        entry.cls,
    )
    assert copy.body_id == 146199 and copy.name == "Eris"


@pytest.mark.parametrize(
    "jd,expected",
    [
        (2451545.0, (251.46102593759298, 10.855337650712395, 31.06450218711325)),
        (2460310.5, (299.36410251597465, -2.7674392832960355, 35.84805539250601)),
    ],
)
def test_direct_spk_retains_the_selected_images_own_center_chain(jd, expected):
    import hashlib
    from libephemeris import spk

    main = DATA_ROOT / "de440.bsp"
    satellite = DATA_ROOT / "tmp/planet_centers/plu017.bsp"
    if not main.is_file() or not satellite.is_file():
        pytest.skip("the same DE440 and Pluto satellite comparison images are required")
    for path, digest in [
        (main, "a4ce9bf9b3282becc9f4b2ac3cebe03a2ae7599981aabd7265fd8482fff7c4b5"),
        (satellite, "419eb90ad2ce723f8f61e2416c2ce9f01307de51d492ba268926a30cdf359bc0"),
    ]:
        with path.open("rb") as source:
            assert hashlib.file_digest(source, "sha256").hexdigest() == digest
    eph.set_calc_mode("skyfield")
    eph.set_ephe_path(str(main.parent))
    eph.set_ephemeris_file(main.name)
    try:
        eph.register_spk_body(42000, str(satellite), 901)
        actual = spk.calc_spk_body_position(jd, 42000, eph.FLG_TRUEPOS)
        # Values are the authenticated pinned-reference geometric places on
        # these exact two files; the tolerances are PLAN section 4 unchanged.
        assert abs(math.remainder(actual[0] - expected[0], 360.0)) * 3600 <= 0.001
        assert abs(actual[1] - expected[1]) * 3600 <= 0.001
        assert abs(actual[2] / expected[2] - 1) <= 1e-8
        eph.set_calc_mode("leb")
        with pytest.raises(RuntimeError):
            spk.calc_spk_body_position(jd, 42000, eph.FLG_TRUEPOS)
    finally:
        eph.unregister_spk_body(42000)
