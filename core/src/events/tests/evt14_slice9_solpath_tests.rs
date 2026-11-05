use super::*;

/// Global-solar fixture over the sealed asset (same shape as the
/// closure discriminant: sealed core, Espenak-Meeus Delta-T).
/// Returns the owned pieces the bound serving borrows; the caller
/// keeps them alive. `SKIP`s (never fails) when the asset is absent.
macro_rules! bind_solglob {
    ($bytes:ident, $index:ident, $cfg:ident, $provider:ident, $nutation_at:ident,
         $serving:ident, $deltat_at:ident, $bound:ident) => {
        let path = "/Users/giacomo/.libephemeris/leb/base_core.leb2";
        let $bytes = match std::fs::read(path) {
            Ok(bytes) => bytes,
            Err(_) => {
                eprintln!("evt14_slice9_solpath_tests: SKIP (missing {path})");
                return;
            }
        };
        let $index = crate::leb::FileIndex::decode(&$bytes).expect("core LEB decodes");
        let $cfg = crate::time::DeltaTConfig {
            model: crate::time::DeltaTModel::EspenakMeeus,
            user_override_days: None,
            tidal_acceleration: crate::constants::TIDAL_DEFAULT,
            iers_enabled: false,
            iers_lookup_seconds: None,
            smh2016_lookup_seconds: None,
        };
        let $provider = crate::positions::SingleLeb::new(&$index, &$bytes);
        let $nutation_at = |jd: f64| $index.eval_nutation(&$bytes, jd);
        let $serving = crate::positions::SealedServing::new(&$provider, &$nutation_at);
        let $deltat_at = |jd_ut: f64| crate::time::deltat_with_config(jd_ut, &$cfg);
        let $bound =
            SealedSolarGlobalServing::try_new(&$serving, &$deltat_at, crate::constants::FLG_SWIEPH)
                .expect("global-solar serving binds");
    };
}

#[test]
fn sealed_path_words_match_sealed_acquisitions() {
    bind_solglob!(
        bytes,
        index,
        cfg,
        provider,
        nutation_at,
        serving,
        deltat_at,
        bound
    );
    // Fixed probe word inside asset coverage (2024-04-08 greatest
    // eclipse neighbourhood, UT).
    let jd = 2460409.26;
    // `path_where` is the `sol_eclipse_where` triple: the where-core
    // character plus the shadow-center ground point.
    let (retflag, lon, lat) = bound.path_where(jd).expect("path where serves");
    let geom = bound
        .shadow_center_geometry_at_ut(jd)
        .expect("shadow center serves");
    assert_eq!(retflag, geom.retc());
    assert_eq!(lon.to_bits(), geom.lon_deg().to_bits());
    assert_eq!(lat.to_bits(), geom.lat_deg().to_bits());
    // The limit sampler reads the same Besselian core as the
    // central-line gamma gate.
    let (x, y) = bound.path_besselian_xy(jd).expect("besselian xy serves");
    let (lx, ly, d, mu, l2) = bound.limit_besselian(jd).expect("limit serves");
    assert_eq!(x.to_bits(), lx.to_bits());
    assert_eq!(y.to_bits(), ly.to_bits());
    assert!(d.is_finite() && mu.is_finite() && l2.is_finite());
    // 2024-04-08 is total: the umbral cone radius is negative while
    // the umbral vertex lies beyond the fundamental plane.
    assert!(l2 < 0.0);
    // Non-finite epochs are refused on every path entry, never walked.
    assert!(bound.path_where(f64::NAN).is_err());
    assert!(bound.limit_besselian(f64::NAN).is_err());
    assert!(sol_eclipse_central_line_sealed(&bound, f64::NAN, jd, 10.0).is_err());
    assert!(sol_eclipse_northern_limit_sealed(&bound, f64::NAN, jd, 10.0).is_err());
    assert!(sol_eclipse_southern_limit_sealed(&bound, f64::NAN, jd, 10.0).is_err());
    // The step gate fires before any provider call.
    assert!(sol_eclipse_central_line_sealed(&bound, jd, jd + 0.2, 0.0).is_err());
}

#[test]
fn sealed_path_covers_2024_total_eclipse() {
    bind_solglob!(
        bytes,
        index,
        cfg,
        provider,
        nutation_at,
        serving,
        deltat_at,
        bound
    );
    // ±2.8 h around the 2024-04-08 maximum (≈2460409.26) at
    // 10-minute steps.
    let (start, end, step) = (2460409.15, 2460409.38, 10.0);
    let (times, lats, lons) =
        sol_eclipse_central_line_sealed(&bound, start, end, step).expect("central serves");
    assert!(!times.is_empty(), "central line covers the 2024 eclipse");
    assert_eq!(times.len(), lats.len());
    assert_eq!(times.len(), lons.len());
    for w in times.windows(2) {
        assert!(w[0] < w[1], "central-line times ascend");
    }
    for ((t, la), lo) in times.iter().zip(&lats).zip(&lons) {
        assert!((*t >= start) && (*t <= end), "sample stays in window");
        assert!((*la >= -90.0) && (*la <= 90.0), "latitude in range");
        assert!((*lo > -180.0) && (*lo <= 180.0), "longitude in range");
    }
    let (_, north_lats, _) =
        sol_eclipse_northern_limit_sealed(&bound, start, end, step).expect("north serves");
    let (_, south_lats, _) =
        sol_eclipse_southern_limit_sealed(&bound, start, end, step).expect("south serves");
    assert!(!north_lats.is_empty() && !south_lats.is_empty());
    let mean = |xs: &[f64]| xs.iter().sum::<f64>() / xs.len() as f64;
    assert!(
        mean(&north_lats) > mean(&south_lats),
        "northern edge sits above the southern edge"
    );
}
