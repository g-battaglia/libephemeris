//! Public Session observer wiring, using pinned-source photometry values.
#[allow(
    dead_code,
    reason = "Shared test data helpers also serve native unit tests."
)]
#[path = "../src/foundation/test_data.rs"]
mod test_data;
use libephemeris_core::{
    constants::*, context::EphemerisContext, native as api, time::DeltaTModel,
};
use sha2::{Digest, Sha256};

#[test]
fn session_topocentric_photometry_uses_bound_observer() {
    // These unchanged reference words belong to data-v3. Use its authenticated
    // historical image, not the newer data-v4 image shipped with the beta.
    let path = test_data::historical_core_path();
    let bytes = std::fs::read(&path).expect("provision the pinned fixture");
    let sha: String = Sha256::digest(&bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect();
    assert_eq!(
        sha,
        "5d708bdbe3e799e0802ba575984e57a3c5e44720dbfa1b4a01cf826640e0cb82"
    );
    let ctx = EphemerisContext::new()
        .with_delta_t_model(Some(DeltaTModel::EspenakMeeus))
        .with_iers_delta_t(Some(false));
    let mut builder = api::SessionBuilder::from_context(ctx);
    builder.leb_file(&path).unwrap();
    let no_site = builder.build().unwrap();
    assert!(matches!(
        no_site.pheno_ut(2451545.0, MARS, FLG_TOPOCTR),
        Err(api::Error::ConfigurationMissing { .. })
    ));
    builder.topo(13.405, 52.52, 34.0);
    let site = builder.build().unwrap();
    let output = site
        .pheno_ut(2451545.0, MARS, FLG_SWIEPH | FLG_TOPOCTR)
        .unwrap();
    // Released pheno_ut, exact data-v3 SHA above, Espenak-Meeus/IERS off,
    // Berlin: also acquired by check_pheno_topocentric.py (UT/Mars/J2000).
    let source = [
        0x403f76e7f78a094e,
        0x3feda5bd85bb20c6,
        0x4047cd4f975272e6,
        0x3f56ff735f2465bc,
        0x3ff092c3bff3a5bd,
    ];
    for (x, word) in output.values.iter().zip(source) {
        assert!((x - f64::from_bits(word)).abs() <= 1e-9);
    }
    assert_eq!(&output.values[5..], &[0.0; 15]);
    assert_ne!(
        output.values,
        no_site
            .pheno_ut(2451545.0, MARS, FLG_SWIEPH)
            .unwrap()
            .values
    );
    assert_observer_and_time_contract(&site, &no_site);
}

#[test]
fn published_data_v4_preserves_observer_and_time_contract() {
    // Exercise the newly published data without fitting new expected numbers
    // from this implementation or reusing references captured with other bytes.
    let path = test_data::published_core_path();
    let context = EphemerisContext::new()
        .with_delta_t_model(Some(DeltaTModel::EspenakMeeus))
        .with_iers_delta_t(Some(false));
    let mut builder = api::SessionBuilder::from_context(context);
    builder.leb_file(&path).unwrap();
    let no_site = builder.build().unwrap();
    assert!(matches!(
        no_site.pheno_ut(2451545.0, MARS, FLG_TOPOCTR),
        Err(api::Error::ConfigurationMissing { .. })
    ));
    builder.topo(13.405, 52.52, 34.0);
    let site = builder.build().unwrap();
    let values = site
        .pheno_ut(2451545.0, MARS, FLG_SWIEPH | FLG_TOPOCTR)
        .unwrap()
        .values;
    assert!(values.iter().all(|value| value.is_finite()));
    assert_eq!(&values[5..], &[0.0; 15]);
    assert_observer_and_time_contract(&site, &no_site);
}

fn assert_observer_and_time_contract(site: &api::Session, no_site: &api::Session) {
    assert_ne!(
        site.pheno_ut(2451545.0, MARS, FLG_SWIEPH | FLG_TOPOCTR)
            .unwrap()
            .values,
        no_site
            .pheno_ut(2451545.0, MARS, FLG_SWIEPH)
            .unwrap()
            .values
    );
    // TT's public two-step model conversion is observable, not direct serving.
    let tt = 2451545.0;
    let ut = tt - site.deltat(tt).unwrap();
    assert_eq!(
        site.pheno(tt, MARS, FLG_TOPOCTR).unwrap().values,
        site.pheno_ut(ut, MARS, FLG_TOPOCTR).unwrap().values
    );
}
