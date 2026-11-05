"""Native catalog snapshots, proper motion, star flags and batch ordering."""
from data_paths import CORE_LEB

from concurrent.futures import ThreadPoolExecutor
import inspect
from pathlib import Path

import libephemeris as eph
from libephemeris import fixed_stars as stars
import pytest


@pytest.fixture(autouse=True)
def configured_engine():
    eph.set_precision_tier("base")
    eph.set_calc_mode("leb")
    eph.set_leb_file(str(CORE_LEB))
    eph.set_iers_delta_t_enabled(False)
    eph.set_delta_t_userdef(64.0 / 86400.0)
    eph.set_topo(12.5, 41.9, 150.0)
    eph.set_sid_mode(0)


def test_catalog_value_records_and_named_constructor():
    catalog = stars.list_fixed_stars()
    assert isinstance(catalog, tuple) and len(catalog) == 1447
    assert all(isinstance(row, stars.StarCatalogEntry) for row in catalog)
    original = catalog[0]
    copy = stars.StarCatalogEntry(original.id, original.name, original.nomenclature,
                                  original.hip_number, original.data, original.magnitude)
    assert copy == original
    assert isinstance(copy.data, stars.StarData)
    assert hasattr(copy.data, "radial_km_per_s")
    with pytest.raises(AttributeError):
        copy.magnitude = 0.0
    assert inspect.signature(stars.fixstar_ut).parameters["flags"].default == 2
    assert inspect.signature(stars.batch_fixstars_ut).parameters["skip_errors"].kind is inspect.Parameter.KEYWORD_ONLY


def test_linear_proper_motion_has_angular_units_and_a_stationary_epoch():
    assert stars.propagate_proper_motion(120.0, 30.0, 0.1, -0.2, 2451545.0, 2451545.0) == (120.0, 30.0)
    ra, dec = stars.propagate_proper_motion(120.0, 0.0, 3600.0, -3600.0, 2451545.0, 2451910.25)
    assert (ra, dec) == (121.0, -1.0)
    ra, dec = stars.propagate_proper_motion(120.0, 90.0, 10.0, 10.0, 2451545.0, 2488070.0)
    assert (ra, dec) == (120.0, 90.0)


def test_scalar_families_and_magnitude_lookup():
    for function in (stars.fixstar, stars.fixstar_ut, stars.fixstar2, stars.fixstar2_ut):
        values, name, flags = function("Spica", 2451545.0, 258)
        assert isinstance(values, tuple) and len(values) == 6
        assert 0.0 <= values[0] < 360.0 and -90.0 <= values[1] <= 90.0
        assert values[2] > 1000.0
        assert "Spica" in name and isinstance(flags, int)
    assert stars.fixstar_mag("Spica") == stars.fixstar2_mag(",alVir")
    with pytest.raises(eph.StarNotFoundError):
        stars.fixstar_ut("not-a-catalog-star", 2451545.0)


@pytest.mark.parametrize('name', ['10CMa', '10Tau', '109Her', '109Vir'])
def test_digit_leading_names_select_their_own_catalog_astrometry(name):
    entry = next(row for row in stars.list_fixed_stars() if row.name == name)
    flags = eph.FLG_TRUEPOS | eph.FLG_EQUATORIAL | eph.FLG_J2000
    for function in (stars.fixstar, stars.fixstar_ut, stars.fixstar2, stars.fixstar2_ut):
        values, canonical, _ = function(name, 2451545., flags)
        # Numeric catalog IDs bypass name parsing. Compare the same public
        # clock/frame/observer route, including its actual frame-bias rotation.
        numeric = eph.calc_ut if function in (stars.fixstar_ut, stars.fixstar2_ut) else eph.calc
        expected, _ = numeric(2451545., entry.id, flags)
        assert canonical.split(',')[0] == name
        assert abs(values[0] - expected[0]) * 3600 < .001
        assert abs(values[1] - expected[1]) * 3600 < .001
    assert stars.fixstar_mag(name)[0] == stars.fixstar2_mag(name)[0] == entry.magnitude


def test_batch_preserves_order_and_explicit_missing_slots():
    names = ["Spica", "Regulus", "Polaris"]
    batch = stars.batch_fixstars_ut(names, 2451545.0, 258)
    assert isinstance(batch, tuple) and len(batch) == len(names)
    for output, name in zip(batch, names):
        assert name in output[1]
    skipped = stars.batch_fixstars_ut(["Spica", "missing", "Regulus"], 2451545.0, skip_errors=True)
    assert skipped[1] is None and "Spica" in skipped[0][1] and "Regulus" in skipped[2][1]
    assert stars.batch_fixstars_ut([], 2451545.0) == ()
    with pytest.raises(eph.StarNotFoundError):
        stars.batch_fixstars_ut(["Spica", "missing"], 2451545.0)


def test_fixed_star_calls_are_deterministic_across_eight_workers():
    names = ("Spica", "Regulus", "Sirius", "Polaris") * 4
    expected = [stars.fixstar_ut(name, 2451545.0, 258) for name in names]
    with ThreadPoolExecutor(max_workers=8) as pool:
        results = list(pool.map(lambda name: stars.fixstar_ut(name, 2451545.0, 258), names))
    assert results == expected


def test_native_photometry_exposes_only_two_reviewed_johnson_v_records():
    catalog = stars.list_fixed_stars()
    qualified = [row for row in catalog if row.photometry.catalogue_field !=
                 stars.PhotometricCatalogueField.Unqualified]
    assert {row.hip_number for row in qualified} == {65474, 80763}
    for row in qualified:
        value = row.photometry
        assert isinstance(value, stars.StarPhotometry)
        assert value == stars.star_photometry(row.id)
        assert value.value == row.magnitude == (0.98 if row.hip_number == 65474 else 1.06)
        assert value.band == stars.PhotometricBand.JohnsonV
        assert value.catalogue_field == stars.PhotometricCatalogueField.HipparcosMainVmag
        assert value.statistic == stars.PhotometricStatistic.Unspecified
        assert value.epoch is None
        assert value.variability == stars.PhotometricVariability.Unknown
        assert value.historical_unqualified_value == (1.04 if row.hip_number == 65474 else None)
        with pytest.raises(AttributeError):
            value.value = 0.0
    for row in catalog:
        if row.hip_number in (65474, 80763):
            continue
        value = row.photometry
        assert value.value == row.magnitude
        assert value.band == stars.PhotometricBand.Unspecified
        assert value.catalogue_field == stars.PhotometricCatalogueField.Unqualified
        assert value.statistic == stars.PhotometricStatistic.Unspecified
        assert value.epoch is None
        assert value.variability == stars.PhotometricVariability.Unknown
        assert value.historical_unqualified_value is None
    assert stars.PhotometricVariability.Unknown != stars.PhotometricVariability.Constant


def test_photometry_keeps_magnitude_aliases_missing_ids_and_constructor_signature():
    for name, designation, expected in (("Spica", ",alVir", 0.98), ("Antares", ",alSco", 1.06)):
        for function in (stars.fixstar_mag, stars.fixstar2_mag):
            assert function(name) == function(designation)
            assert function(name)[0] == expected
    assert stars.star_photometry(-1) is None
    assert stars.star_photometry(0) is None
    with pytest.raises(eph.StarNotFoundError):
        stars.fixstar_mag("missing-photometry")
    with pytest.raises(eph.StarNotFoundError):
        stars.fixstar2_mag("missing-photometry")
    original = next(row for row in stars.list_fixed_stars() if row.name == "Spica")
    copied = stars.StarCatalogEntry(original.id, original.name, original.nomenclature,
                                    original.hip_number, original.data, original.magnitude)
    assert copied == original
    assert copied.magnitude == original.magnitude
    # A user-supplied record has no authenticated catalogue provenance.
    assert copied.photometry.catalogue_field == stars.PhotometricCatalogueField.Unqualified
    assert copied.photometry.value == original.magnitude
    assert original.photometry.catalogue_field == stars.PhotometricCatalogueField.HipparcosMainVmag
    assert not hasattr(copied, "catalogue_origin")
    for row in (original, copied):
        with pytest.raises(AttributeError):
            row.photometry = copied.photometry
        with pytest.raises(AttributeError):
            row.magnitude = 7.0
