"""Native N-body integration, physical invariants and reusable file ownership."""

from concurrent.futures import ThreadPoolExecutor
import cmath
import gc
import hashlib
import inspect
import json
import math
import os
from pathlib import Path
import struct
from threading import Barrier
import warnings

import libephemeris as eph
from libephemeris import rebound_integration as nbody
from libephemeris.minor_bodies import OrbitalElements
import pytest
from data_paths import DATA_ROOT


JD = 2451545.0
MU = 0.00029591220828559
ASSIST_DIR = DATA_ROOT / "assist"
ASSIST_HASHES = {
    "linux_p1550p2650.440": "29915576d0a6555766b99485ac3056ee415e86df4fce282611c31afb329ad062",
    "sb441-n16.bsp": "919d612ce3c72a78fc7158f9120156542d0f21e6b8b052e4c1339c759747fd90",
}


def elements():
    return OrbitalElements("comparison", JD, 2.7, 0.12, 7.0, 40.0, 65.0, 30.0, 0.222)


def verify_assist_files():
    """Bind comparisons to unchanged selected force files, before and after use."""
    for name, expected in ASSIST_HASHES.items():
        with (ASSIST_DIR / name).open("rb") as stream:
            assert hashlib.file_digest(stream, "sha256").hexdigest() == expected


@pytest.fixture(scope="module")
def assist_config():
    if any(not (ASSIST_DIR / name).is_file() for name in ASSIST_HASHES):
        pytest.skip("optional local ASSIST data are unavailable")
    verify_assist_files()
    yield nbody.AssistEphemConfig(
        str(ASSIST_DIR / "linux_p1550p2650.440"), str(ASSIST_DIR / "sb441-n16.bsp")
    )
    verify_assist_files()


@pytest.fixture(scope="module")
def alternate_assist_config(assist_config):
    """An actually different existing dataset, used only for isolation controls."""
    path = ASSIST_DIR / "linux_m13000p17000.441"
    expected = "476096486def4e41bfceb29aa27f50784da0bce318902bcf7b88caad058cd4da"
    if not path.is_file():
        pytest.skip("optional alternate classic-DE dataset is unavailable")

    def verify():
        with path.open("rb") as stream:
            assert hashlib.file_digest(stream, "sha256").hexdigest() == expected

    verify()
    yield nbody.AssistEphemConfig(str(path), assist_config.asteroids_file)
    verify()


def test_native_versions_and_named_integrators():
    assert nbody.check_rebound_available() is True
    assert nbody.check_assist_available() is True
    # Compatibility names identify the ordinary Rust implementations. These
    # established strings remain explicit; no environment selector can make an
    # installed legacy backend satisfy the current package contract.
    assert (nbody.get_rebound_version(), nbody.get_assist_version()) == (
        "libephemeris-native-radau15/4.0.0-experimental",
        "libephemeris-native-ppn-icrf-terrestrial-iau-solar/4.0.0-experimental",
    )
    for name in ("IAS15", "WHFAST", "MERCURIUS", "TRACE", "LEAPFROG"):
        value = getattr(nbody.ReboundIntegrator, name)
        assert value.name == name
        assert value.value == name.lower()
        assert nbody.ReboundIntegrator(value.value) == value
    with pytest.raises(ValueError):
        nbody.ReboundIntegrator("invalid")


def test_owned_elements_and_result_record_properties():
    orbit = elements()
    assert orbit.validate() == []
    orbit.a = -1.0
    assert orbit.validate()
    result = nbody.PropagationResult(3.0, 4.0, 0.0, 0.01, 0.02, 0.0, JD)
    assert result.distance == 5.0
    assert result.ecliptic_lat == 0.0
    assert abs(result.ecliptic_lon - math.degrees(math.atan2(4.0, 3.0))) < 1e-13
    assert result.to_tuple() == (3.0, 4.0, 0.0, 0.01, 0.02, 0.0)
    result.z = 5.0
    assert abs(result.distance - math.sqrt(50.0)) < 2e-15
    assert result.ecliptic_lat > 0.0
    zero = nbody.PropagationResult(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, JD)
    assert zero.distance == zero.ecliptic_lat == 0.0


def test_curated_assist_source_is_explicit_and_retained(assist_config):
    from libephemeris.minor_bodies import MINOR_BODY_ELEMENTS

    epoch = MINOR_BODY_ELEMENTS[eph.CERES].epoch
    with pytest.warns(
        UserWarning, match="all other configured asteroid forces remain active"
    ):
        selected = assist_config.propagate_curated(eph.CERES, epoch)
    with pytest.warns(UserWarning):
        repeated = assist_config.propagate_curated(eph.CERES, epoch)
    assert selected.to_tuple() == repeated.to_tuple()
    assert selected.jd_tt == epoch
    # Scalar callers share the retained seed/force owner, not mutable trajectory
    # state. Use different directions and duplicates to exercise that boundary.
    epochs = [epoch, epoch + 1.0, epoch - 1.0, epoch + 32.0] * 2
    with pytest.warns(UserWarning):
        expected = [
            assist_config.propagate_curated(eph.CERES, end).to_tuple() for end in epochs
        ]
    with ThreadPoolExecutor(max_workers=4) as pool:
        actual = list(
            pool.map(
                lambda end: assist_config.propagate_curated(eph.CERES, end).to_tuple(),
                epochs,
            )
        )
    assert actual == expected
    with pytest.warns(UserWarning):
        supplied = nbody.propagate_orbit_assist(
            MINOR_BODY_ELEMENTS[eph.CERES], epoch, epoch, assist_config
        )
    assert selected.to_tuple() != supplied.to_tuple()
    with pytest.raises(ValueError):
        assist_config.propagate_curated(eph.CERES, math.nan)


def test_catalogue_identity_survives_binding_without_scalar_identity_inference(
    assist_config,
):
    from libephemeris.minor_bodies import MINOR_BODY_ELEMENTS

    catalogue = MINOR_BODY_ELEMENTS[eph.CERES]
    fields = ("name", "epoch", "a", "e", "i", "omega", "Omega", "M0", "n")
    anonymous = OrbitalElements(*(getattr(catalogue, field) for field in fields))
    for record in (catalogue, anonymous):
        public_fields = {
            name
            for name in dir(record)
            if not name.startswith("_") and not callable(getattr(record, name))
        }
        assert public_fields == set(fields)

    with pytest.warns(UserWarning, match="Ceres"):
        identified = nbody.propagate_orbit_assist(
            catalogue, catalogue.epoch, catalogue.epoch, assist_config
        )
    with warnings.catch_warnings(record=True) as diagnostics:
        warnings.simplefilter("always")
        copied = nbody.propagate_orbit_assist(
            anonymous, anonymous.epoch, anonymous.epoch, assist_config
        )
    assert diagnostics == [], (
        "identical user scalars must not acquire a catalogue identity"
    )
    assert copied.to_tuple() == identified.to_tuple(), (
        "identity must not replace the supplied seed"
    )
    assert copied.jd_tt == identified.jd_tt == catalogue.epoch


def test_spk_seeded_source_orders_both_directions_and_retains_configuration(
    assist_config,
):
    import numpy as np

    main_path = DATA_ROOT / "de440.bsp"
    with main_path.open("rb") as stream:
        assert (
            hashlib.file_digest(stream, "sha256").hexdigest()
            == "a4ce9bf9b3282becc9f4b2ac3cebe03a2ae7599981aabd7265fd8482fff7c4b5"
        )
    context = eph.EphemerisContext(str(DATA_ROOT), main_path.name)
    main = context.get_loader()(str(main_path))
    primary = context.get_loader()(assist_config.asteroids_file)
    seed = primary.geometric_source(2000001, [main])
    config = nbody.AssistEphemConfig(
        assist_config.planets_file, assist_config.asteroids_file
    )
    source = config.nbody_source(seed)
    assert source.anchor_tdb == JD
    assert source.coverage_tdb() == (2287216.5, 2688944.5)
    assert source.states_tdb([]) == []
    epochs = [JD + 100, JD - 100, JD, JD + 100, JD - 200, JD + 200]
    with pytest.warns(UserWarning):
        rows = np.asarray(source.states_tdb(epochs))
    assert rows.shape == (len(epochs), 6)
    assert np.array_equal(rows[0], rows[3])
    for index, epoch in enumerate(epochs):
        with pytest.warns(UserWarning):
            scalar = np.asarray(source.states_tdb([epoch])[0])
        # Independent integration histories must satisfy the unchanged PLAN
        # position/rate budgets, rather than require solver-step bit identity.
        assert np.linalg.norm(rows[index, :3] - scalar[:3]) <= 1e-8 * np.linalg.norm(
            scalar[:3]
        )
        assert np.linalg.norm(rows[index, 3:] - scalar[3:]) <= 1e-6 * np.linalg.norm(
            scalar[3:]
        )
    for invalid in (math.nan, math.inf, source.coverage_tdb()[0] - 1):
        with pytest.raises(ValueError):
            source.states_tdb([JD, invalid])
    with pytest.raises((ValueError, eph.Error), match="center"):
        config.nbody_source(main.geometric_source(10))
    config.planets_file = "/libephemeris-nbody-missing-selected-force-file"
    context.close()
    with pytest.warns(UserWarning):
        assert len(source.states_tdb([JD + 1])) == 1
    with pytest.raises((OSError, eph.Error)):
        config.nbody_source(seed)


def test_offline_tt_source_applies_clock_chain_rule_and_recovers(assist_config):
    # Independently differentiate the published clock itself by complex step,
    # rather than call the native rate or duplicate its analytic derivative.
    # USNO Circular179, eq2.6 (seconds; T in TT centuries from J2000):
    # https://aa.usno.navy.mil/downloads/Circular_179.pdf
    terms = (
        (0.001657, 628.3076, 6.2401),
        (0.000022, 575.3385, 4.2970),
        (0.000014, 1256.6152, 6.1969),
        (0.000005, 606.9777, 4.0212),
        (0.000005, 52.9691, 0.4444),
        (0.000002, 21.3299, 5.5431),
    )

    def clock_seconds(tt):
        centuries = (tt - JD) / 36525.0
        return sum(
            amplitude * cmath.sin(frequency * centuries + phase)
            for amplitude, frequency, phase in terms
        ) + 0.000010 * centuries * cmath.sin(628.3076 * centuries + 4.2490)

    def rate(tt, step):
        return 1.0 + clock_seconds(complex(tt, step)).imag / (step * 86400.0)

    main_path = DATA_ROOT / "de440.bsp"
    expected_hash = "a4ce9bf9b3282becc9f4b2ac3cebe03a2ae7599981aabd7265fd8482fff7c4b5"

    def verify_main():
        with main_path.open("rb") as stream:
            assert hashlib.file_digest(stream, "sha256").hexdigest() == expected_hash

    def build():
        context = eph.EphemerisContext(str(DATA_ROOT), main_path.name)
        loader = context.get_loader()
        main = loader(str(main_path))
        primary = loader(assist_config.asteroids_file)
        seed = primary.geometric_source(2000001, [main])
        config = nbody.AssistEphemConfig(
            assist_config.planets_file, assist_config.asteroids_file
        )
        return config.nbody_source(seed)

    verify_main()
    source = build()
    gc.collect()  # No Python seed/loader/configuration object survives build().
    assert source.states_tt([]) == []
    epochs = [JD + 8.0, JD - 8.0, JD, JD + 8.0, JD + 182.25]
    tdb = [tt + clock_seconds(tt).real / 86400.0 for tt in epochs]
    rates = [rate(tt, 1e-4) for tt in epochs]
    assert rates == [rate(tt, 1e-5) for tt in epochs]
    assert min(rates) < 1.0 < max(rates), "clock-rate control must exercise both signs"
    with pytest.warns(UserWarning):
        expected = source.states_tdb(tdb)
    with pytest.warns(UserWarning):
        actual = source.states_tt(epochs)
    assert len(actual) == len(epochs)
    assert actual[0] == actual[3]
    for tt_row, tdb_row, factor in zip(actual, expected, rates, strict=True):
        # Identical physical epochs and integration order isolate the clock
        # unit conversion from orbit/model/solver error. Do not loosen PLAN.
        assert struct.pack(">3d", *tt_row[:3]) == struct.pack(">3d", *tdb_row[:3])
        for value, raw in zip(tt_row[3:], tdb_row[3:], strict=True):
            converted = raw * factor
            rounding = 4 * math.ulp(converted)
            assert abs(value - converted) <= rounding
            # Positive calibration: an omitted velocity scale must fail the
            # same arithmetic check, even though the factor is only ~1e-10.
            assert abs(raw - converted) > 100 * rounding
    for invalid in (math.nan, math.inf, source.coverage_tdb()[0] - 1.0):
        with pytest.raises(ValueError):
            source.states_tt([JD, invalid, JD + 1.0])
    with pytest.warns(UserWarning):
        assert source.states_tt(epochs) == actual
    verify_main()


@pytest.mark.parametrize("name", ["IAS15", "WHFAST", "MERCURIUS", "TRACE", "LEAPFROG"])
def test_circular_two_body_energy_angular_momentum_and_radius(name):
    orbit = OrbitalElements(
        "circular", JD, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, math.degrees(math.sqrt(MU))
    )
    result = nbody.propagate_orbit_rebound(
        orbit, JD, JD + 30.0, getattr(nbody.ReboundIntegrator, name), 0.05
    )
    radius = result.distance
    speed2 = result.vx**2 + result.vy**2 + result.vz**2
    energy = speed2 / 2.0 - MU / radius
    angular_momentum = math.sqrt(
        (result.y * result.vz - result.z * result.vy) ** 2
        + (result.z * result.vx - result.x * result.vz) ** 2
        + (result.x * result.vy - result.y * result.vx) ** 2
    )
    assert abs(radius - 1.0) < 1e-6
    assert abs(energy / (-MU / 2.0) - 1.0) < 1e-6
    assert abs(angular_momentum / math.sqrt(MU) - 1.0) < 1e-6
    assert result.jd_tt == JD + 30.0


def test_finite_step_and_zero_duration_owner_contract():
    orbit = elements()
    initial = nbody.propagate_orbit_rebound(orbit, JD, JD)
    zero_step = nbody.propagate_orbit_rebound(orbit, JD, JD, dt=0.0)
    assert initial.to_tuple() == zero_step.to_tuple()
    for step in (math.nan, math.inf, -math.inf):
        with pytest.raises(ValueError):
            nbody.propagate_orbit_rebound(orbit, JD, JD, dt=step)
    with pytest.raises(ValueError):
        nbody.propagate_orbit_rebound(orbit, JD, JD + 1.0, dt=0.0)
    with pytest.raises(ValueError):
        nbody.propagate_orbit_rebound(orbit, math.nan, JD)
    with pytest.raises(ValueError):
        nbody.PropagationResult(math.nan, 0.0, 0.0, 0.0, 0.0, 0.0, JD)
    with pytest.raises(TypeError):
        nbody.propagate_orbit_rebound(object(), JD, JD + 1.0)


def test_sample_epochs_forward_backward_and_empty():
    orbit = elements()
    for end in (JD + 6.0, JD - 6.0):
        results = nbody.propagate_trajectory(orbit, JD, end, 4, False)
        assert [item.jd_tt for item in results] == [
            JD + (end - JD) * k / 3 for k in range(4)
        ]
        assert all(
            math.isfinite(value) for result in results for value in result.to_tuple()
        )
    assert nbody.propagate_trajectory(orbit, JD, JD + 1.0, 0, False) == []
    assert nbody.propagate_trajectory(orbit, JD, JD + 1.0, -1, False) == []
    one = nbody.propagate_trajectory(orbit, JD, JD + 1.0, 1, False)
    assert len(one) == 1 and one[0].jd_tt == JD


def test_missing_assist_data_and_documented_fallback(tmp_path):
    config = nbody.AssistEphemConfig(
        str(tmp_path / "missing.440"), str(tmp_path / "missing.bsp")
    )
    with pytest.raises(FileNotFoundError):
        nbody.propagate_orbit_assist(elements(), JD, JD + 1.0, config)
    fallback = nbody.propagate_trajectory(elements(), JD, JD + 1.0, 3, True, config)
    rebound = nbody.propagate_trajectory(elements(), JD, JD + 1.0, 3, False)
    assert [value.to_tuple() for value in fallback] == [
        value.to_tuple() for value in rebound
    ]
    diagnostic = nbody.compare_with_keplerian(elements(), JD + 40.0, True, config)
    assert diagnostic["method"] == "assist"
    assert diagnostic["propagation_days"] == 40.0
    assert set(diagnostic["nbody"]) == {"x", "y", "z", "lon", "lat", "dist"}


def test_real_assist_reuse_non_gravitational_and_threads(assist_config):
    orbit = elements()
    expected = nbody.propagate_orbit_assist(orbit, JD, JD + 10.0, assist_config)
    assert all(math.isfinite(value) for value in expected.to_tuple())
    assert expected.jd_tt == JD + 10.0
    with ThreadPoolExecutor(max_workers=8) as pool:
        values = list(
            pool.map(
                lambda _: nbody.propagate_orbit_assist(
                    orbit, JD, JD + 10.0, assist_config
                ),
                range(8),
            )
        )
    assert all(value.to_tuple() == expected.to_tuple() for value in values)
    changed = nbody.propagate_orbit_assist(
        orbit, JD, JD + 10.0, assist_config, True, 1e-9, 2e-10, -1e-10
    )
    assert changed.to_tuple() != expected.to_tuple()
    trajectory = nbody.propagate_trajectory(
        orbit, JD, JD + 10.0, 3, True, assist_config
    )
    assert [point.jd_tt for point in trajectory] == [JD, JD + 5.0, JD + 10.0]
    assert all(
        abs(a - b) < 1e-12
        for a, b in zip(trajectory[-1].to_tuple(), expected.to_tuple())
    )


def test_eight_thread_rebound_matches_sequential():
    orbit = elements()
    ends = [JD + days for days in (1.0, 2.0, -1.0, 4.0, -3.0, 6.0, 7.0, 8.0)]
    expected = [
        nbody.propagate_orbit_rebound(orbit, JD, end).to_tuple() for end in ends
    ]
    with ThreadPoolExecutor(max_workers=8) as pool:
        results = list(
            pool.map(
                lambda end: nbody.propagate_orbit_rebound(orbit, JD, end).to_tuple(),
                ends,
            )
        )
    assert results == expected


def result_words(result):
    """Preserve all six PV words and the returned epoch, including signed zero."""
    values = (*result.to_tuple(), result.jd_tt)
    assert len(values) == 7 and all(math.isfinite(value) for value in values)
    return struct.pack(">7d", *values)


@pytest.mark.parametrize("curated", [False, True], ids=["elements", "explicit-seed"])
def test_all_catalogue_bodies_eight_threads_preserve_scalar_results(
    assist_config, curated
):
    from libephemeris.minor_bodies import MINOR_BODY_ELEMENTS

    assert len(MINOR_BODY_ELEMENTS) == 37
    # Nonmonotone directions and a repeated anchor expose shared trajectory
    # state accidentally leaking between otherwise independent scalar calls.
    requests = [
        (body, offset)
        for body in sorted(MINOR_BODY_ELEMENTS)
        for offset in (0.0, 1.0, -1.0, 32.0, -32.0, 0.0)
    ]

    def evaluate(request):
        body, offset = request
        orbit = MINOR_BODY_ELEMENTS[body]
        endpoint = orbit.epoch + offset
        if curated:
            result = assist_config.propagate_curated(body, endpoint)
        else:
            result = nbody.propagate_orbit_assist(
                orbit, orbit.epoch, endpoint, assist_config
            )
        return result_words(result)

    sequential = [evaluate(request) for request in requests]
    with ThreadPoolExecutor(max_workers=8) as pool:
        concurrent = list(pool.map(evaluate, requests))
    assert concurrent == sequential


def test_distinct_sessions_and_force_parameters_do_not_cross_contaminate(assist_config):
    sessions = [
        nbody.AssistEphemConfig(
            assist_config.planets_file, assist_config.asteroids_file
        )
        for _ in range(2)
    ]
    orbit = elements()
    requests = [
        (owner, acceleration, offset)
        for offset in (1.0, -1.0, 8.0, -8.0)
        for owner in (0, 1)
        for acceleration in (0.0, 1e-9)
    ]

    def evaluate(request):
        owner, acceleration, offset = request
        result = nbody.propagate_orbit_assist(
            orbit,
            JD,
            JD + offset,
            sessions[owner],
            True,
            acceleration,
            acceleration / 5,
            -acceleration / 10,
        )
        return result_words(result)

    sequential = [evaluate(request) for request in requests]
    # Calibration: differing nonzero force parameters really exercise distinct
    # numerical outcomes; merely sharing identical calls would miss leakage.
    assert sequential[0] != sequential[1]
    assert sequential[0] == sequential[2]
    with ThreadPoolExecutor(max_workers=8) as pool:
        concurrent = list(pool.map(evaluate, reversed(requests)))
    assert concurrent == list(reversed(sequential))
    assert [evaluate(request) for request in requests] == sequential


def test_concurrent_invalid_requests_leave_retained_session_usable(assist_config):
    orbit = elements()
    expected = result_words(
        nbody.propagate_orbit_assist(orbit, JD, JD + 2, assist_config)
    )

    def evaluate(valid):
        if valid:
            return result_words(
                nbody.propagate_orbit_assist(orbit, JD, JD + 2, assist_config)
            )
        with pytest.raises(ValueError):
            nbody.propagate_orbit_assist(orbit, JD, math.nan, assist_config)
        return None

    requests = [True, False, False, True] * 8
    with ThreadPoolExecutor(max_workers=8) as pool:
        actual = list(pool.map(evaluate, requests))
    assert actual == [expected if valid else None for valid in requests]
    assert evaluate(True) == expected


def test_offline_source_outlives_all_python_input_owners(assist_config):
    main_path = DATA_ROOT / "de440.bsp"
    expected_hash = "a4ce9bf9b3282becc9f4b2ac3cebe03a2ae7599981aabd7265fd8482fff7c4b5"
    with main_path.open("rb") as stream:
        assert hashlib.file_digest(stream, "sha256").hexdigest() == expected_hash
    epochs = [JD - 32, JD, JD + 32, JD - 32]

    def states(source):
        with pytest.warns(UserWarning):
            rows = source.states_tdb(epochs)
        assert len(rows) == len(epochs)
        assert all(
            len(row) == 6 and all(math.isfinite(value) for value in row) for row in rows
        )
        return [struct.pack(">6d", *row) for row in rows]

    def build():
        # Only the returned source and numerical bytes escape this scope. The
        # context, loader, kernels, seed and configuration must not be needed
        # as Python lifetime anchors for the native immutable owners.
        context = eph.EphemerisContext(str(DATA_ROOT), main_path.name)
        loader = context.get_loader()
        main = loader(str(main_path))
        primary = loader(assist_config.asteroids_file)
        seed = primary.geometric_source(2000001, [main])
        config = nbody.AssistEphemConfig(
            assist_config.planets_file, assist_config.asteroids_file
        )
        source = config.nbody_source(seed)
        return source, states(source)

    source, expected = build()
    gc.collect()
    assert states(source) == expected
    with pytest.raises(ValueError):
        source.states_tdb([JD, math.nan])
    assert states(source) == expected
    with main_path.open("rb") as stream:
        assert hashlib.file_digest(stream, "sha256").hexdigest() == expected_hash


def test_individual_nongravitational_channels_are_distinct_and_isolated(assist_config):
    orbit = elements()
    coefficients = [
        (0.0, 0.0, 0.0),
        (1e-9, 0.0, 0.0),
        (0.0, 1e-9, 0.0),
        (0.0, 0.0, 1e-9),
    ]

    def evaluate(index):
        return result_words(
            nbody.propagate_orbit_assist(
                orbit,
                JD,
                JD + 8,
                assist_config,
                True,
                *coefficients[index],
            )
        )

    expected = [evaluate(index) for index in range(4)]
    # This calibrates each public radial/transverse/normal parameter separately;
    # it does not replace the Rust force-law tests or claim physical accuracy.
    assert len(set(expected)) == 4
    requests = [3, 0, 2, 1, 1, 2, 0, 3] * 4
    with ThreadPoolExecutor(max_workers=8) as pool:
        actual = list(pool.map(evaluate, requests))
    assert actual == [expected[index] for index in requests]
    assert [evaluate(index) for index in range(4)] == expected


def test_different_force_datasets_and_fresh_owners_remain_isolated(
    assist_config, alternate_assist_config
):
    configs = [assist_config, alternate_assist_config]
    orbit = elements()
    endpoints = [JD + 1, JD - 1, JD + 32, JD - 32]
    requests = [(index, endpoint) for endpoint in endpoints for index in range(2)]

    def evaluate(config, endpoint):
        return result_words(nbody.propagate_orbit_assist(orbit, JD, endpoint, config))

    expected = [evaluate(configs[index], endpoint) for index, endpoint in requests]
    # These real DE440/DE441 endpoints can round to identical public words.
    # Therefore this checks fresh-owner concurrency and repeatability only;
    # the separate synthetic Rust GM test positively calibrates data isolation
    # through the complete force and integration route. Never infer distinct
    # physical outputs solely from distinct filenames.
    assert configs[0].planets_file != configs[1].planets_file
    rendezvous = Barrier(8)

    def fresh(request):
        index, endpoint = request
        selected = configs[index]
        owner = nbody.AssistEphemConfig(selected.planets_file, selected.asteroids_file)
        rendezvous.wait(timeout=30)
        return evaluate(owner, endpoint)

    # All eight fresh owners enter their first calls after the same barrier.
    # This exercises fresh-owner admission, not a claim of GIL/native overlap
    # or of cold filesystem/process caches.
    with ThreadPoolExecutor(max_workers=8) as pool:
        actual = list(pool.map(fresh, reversed(requests)))
    assert actual == list(reversed(expected))
    assert [
        evaluate(configs[index], endpoint) for index, endpoint in requests
    ] == expected


def test_repeated_force_owner_release_returns_file_descriptors(assist_config):
    if not Path("/dev/fd").is_dir():
        pytest.skip("descriptor census requires the local POSIX environment")
    orbit = elements()
    gc.collect()
    baseline = len(os.listdir("/dev/fd"))
    for _ in range(12):
        owner = nbody.AssistEphemConfig(
            assist_config.planets_file, assist_config.asteroids_file
        )
        result = nbody.propagate_orbit_assist(orbit, JD, JD + 0.25, owner)
        result_words(result)
        # An explicit invalid request must not prevent release of the owner.
        with pytest.raises(ValueError):
            nbody.propagate_orbit_assist(orbit, JD, math.nan, owner)
        del result, owner
        gc.collect()
        assert len(os.listdir("/dev/fd")) == baseline
