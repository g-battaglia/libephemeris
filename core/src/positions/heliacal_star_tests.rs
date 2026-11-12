//! Catalog-role, owned-clock and lazy-frame controls without an image.
use super::*;
use crate::constants::{SIRIUS, SUN};
use std::cell::Cell;

fn native(context: EphemerisContext) -> NativeHeliacalStar {
    NativeHeliacalStar::new(context, CatalogHeliacalStar::from_id(SIRIUS).unwrap()).unwrap()
}

#[test]
fn direct_star_keeps_only_constructor_fields_and_catalog_identity() {
    let catalog = CatalogHeliacalStar::from_id(SIRIUS).unwrap();
    let original = star_entry_by_id(SIRIUS).unwrap();
    let row = original.data();
    assert_eq!(catalog.id(), SIRIUS);
    assert_eq!(catalog.name(), original.name());
    assert_eq!(
        catalog.magnitude().to_bits(),
        original.magnitude().to_bits()
    );
    assert_eq!(
        catalog.data.ra_j2000_deg().to_bits(),
        row.ra_j2000_deg().to_bits()
    );
    assert_eq!(
        catalog.data.dec_j2000_deg().to_bits(),
        row.dec_j2000_deg().to_bits()
    );
    assert_eq!(
        catalog.data.pm_ra_arcsec_yr().to_bits(),
        row.pm_ra_arcsec_yr().to_bits()
    );
    assert_eq!(
        catalog.data.pm_dec_arcsec_yr().to_bits(),
        row.pm_dec_arcsec_yr().to_bits()
    );
    assert_eq!(catalog.data.parallax_mas().to_bits(), 0.0_f64.to_bits());
    assert_eq!(
        catalog.data.radial_velocity_km_s().to_bits(),
        0.0_f64.to_bits()
    );
    let irrelevant = StarData::from_row(
        row.ra_j2000_deg(),
        row.dec_j2000_deg(),
        row.pm_ra_arcsec_yr(),
        row.pm_dec_arcsec_yr(),
        f64::NAN,
        f64::INFINITY,
    );
    assert_eq!(direct_star_data(irrelevant), catalog.data);
}

#[test]
fn explicit_non_native_modes_do_not_acquire_a_fallback() {
    let catalog = CatalogHeliacalStar::from_id(SIRIUS).unwrap();
    for mode in [
        ProviderMode::Auto,
        ProviderMode::SealedLeb,
        ProviderMode::HorizonsRemote,
    ] {
        let result = NativeHeliacalStar::new(EphemerisContext::new().with_mode(mode), catalog);
        assert!(matches!(
            result,
            Err(Error::ConfigurationMissing {
                missing_config: "native_heliacal_direct_star",
                ..
            })
        ));
    }
}

#[test]
fn astrometry_retains_original_time_and_ignores_configured_model_clocks() {
    let context = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
    let mut other = context.clone();
    other.set_user_delta_t(Some(1000.0));
    other.set_tidal_acceleration(Some(-35.0));
    let first = native(context);
    let second = native(other);
    let time = TwoPartEpoch {
        hi: 2460463.125,
        lo: 0.0008123456789012345,
    };
    let observer = [0.375, -0.75, 0.03125];
    let a = first.astrometric(time, observer);
    let b = second.astrometric(time, observer);
    assert_eq!(a.ray.map(f64::to_bits), b.ray.map(f64::to_bits));
    assert_eq!(a.light_time_days.to_bits(), b.light_time_days.to_bits());
    let expected = astrometric_ray(
        &first.catalog.data,
        NativeStarObserver {
            time,
            position: observer,
            velocity: [123.0; 3],
            velocity_per_tt: [123.0; 3],
            surface: Some([456.0; 3]),
        },
    );
    assert_eq!(a.ray.map(f64::to_bits), expected.0.map(f64::to_bits));
    assert_eq!(a.light_time_days.to_bits(), expected.1.to_bits());
    assert_eq!(first.catalog().id(), SIRIUS);
}

#[test]
fn failed_apparent_support_is_raw_and_precedes_horizon_acquisition() {
    let owner = native(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
    let time = TwoPartEpoch {
        hi: 2460463.125,
        lo: 0.0008,
    };
    let horizon_read = Cell::new(false);
    let actual = owner.apparent_topocentric_place(
        time,
        [1.0, 0.0, 0.0],
        [0.0; 3],
        [0.00004, 0.0, 0.0],
        || {
            horizon_read.set(true);
            Err(Error::NonFiniteArgument { value: f64::NAN })
        },
    );
    assert!(matches!(
        actual,
        Err(Error::SpkUnregisteredBody { body: SUN })
    ));
    assert!(!horizon_read.get());
    assert!(matches!(
        owner.apparent_geocentric_ray(time, [1.0, 0.0, 0.0], [0.0; 3]),
        Err(Error::SpkUnregisteredBody { body: SUN })
    ));
}
