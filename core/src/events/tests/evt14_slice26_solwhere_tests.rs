use super::*;
use crate::constants::{ECL_CENTRAL, ECL_TOTAL};

/// Global-solar fixture over the sealed asset (fixed-instant where
/// needs no site; Espenak-Meeus Delta-T). `SKIP`s (never fails)
/// when the asset is absent.
macro_rules! bind_solglob26 {
    ($provider:ident, $nutation_at:ident, $serving:ident, $deltat_at:ident, $bound:ident) => {
        let path = "/Users/giacomo/.libephemeris/leb/base_core.leb2";
        let bytes = match std::fs::read(path) {
            Ok(bytes) => bytes,
            Err(_) => {
                eprintln!("evt14_slice26_solwhere_tests: SKIP (missing {path})");
                return;
            }
        };
        let index = crate::leb::FileIndex::decode(&bytes).expect("core LEB decodes");
        let cfg = crate::time::DeltaTConfig {
            model: crate::time::DeltaTModel::EspenakMeeus,
            user_override_days: None,
            tidal_acceleration: crate::constants::TIDAL_DEFAULT,
            iers_enabled: false,
            iers_lookup_seconds: None,
            smh2016_lookup_seconds: None,
        };
        let $provider = crate::positions::SingleLeb::new(&index, &bytes);
        let $nutation_at = |jd: f64| index.eval_nutation(&bytes, jd);
        let $serving = crate::positions::SealedServing::new(&$provider, &$nutation_at);
        let $deltat_at = |jd_ut: f64| crate::time::deltat_with_config(jd_ut, &cfg);
        let $bound =
            SealedSolarGlobalServing::try_new(&$serving, &$deltat_at, crate::constants::FLG_SWIEPH)
                .expect("global-solar serving binds");
    };
}

#[test]
fn non_finite_instants_refuse() {
    bind_solglob26!(provider, nutation_at, serving, deltat_at, bound);
    let _ = (&provider, &nutation_at, &serving, &deltat_at);
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            sol_eclipse_where_sealed(&bound, bad),
            Err(Error::NonFiniteArgument { .. })
        ));
    }
}

#[test]
fn total_maximum_reports_central_shadow_center() {
    bind_solglob26!(provider, nutation_at, serving, deltat_at, bound);
    let _ = (&provider, &nutation_at, &serving, &deltat_at);
    // Banked 2024-04-08 total eclipse (Dallas local maximum
    // instant: the shadow axis is over North America): central +
    // total character, shadow-center lon/lat near the Dallas
    // totality site, totality magnitude, finite reserved zeros.
    let (retflag, geopos, attr) = sol_eclipse_where_sealed(&bound, 2460409.279538682).unwrap();
    assert_ne!(retflag & ECL_CENTRAL, 0, "retflag {retflag}");
    assert_ne!(retflag & ECL_TOTAL, 0, "retflag {retflag}");
    assert!(
        (geopos[0] - -96.7970).abs() < 5.0 && (geopos[1] - 32.7767).abs() < 5.0,
        "shadow center {:?}",
        &geopos[0..2]
    );
    assert!(attr[0] >= 1.0, "total magnitude {}", attr[0]);
    assert!(geopos[2..].iter().all(|w| *w == 0.0));
    assert!(geopos[0..2].iter().all(|w| w.is_finite()));
    assert!(attr.iter().all(|w| w.is_finite()));
}

#[test]
fn quiet_instant_reports_no_eclipse() {
    bind_solglob26!(provider, nutation_at, serving, deltat_at, bound);
    let _ = (&provider, &nutation_at, &serving, &deltat_at);
    // Quiet JD: retflag 0 is the normal no-eclipse result (never a
    // refusal); the carried closest-approach point stays finite.
    let (retflag, geopos, attr) = sol_eclipse_where_sealed(&bound, 2460140.0).unwrap();
    assert_eq!(retflag, 0);
    assert!(geopos[0..2].iter().all(|w| w.is_finite()));
    assert!(attr.iter().all(|w| w.is_finite()));
}
