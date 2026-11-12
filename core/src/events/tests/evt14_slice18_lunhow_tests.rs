use super::*;

/// Global-lunar fixture over the sealed asset (same shape as the
/// slice-17 fixture: sealed core, Espenak-Meeus Delta-T).
macro_rules! bind_lunglob18 {
    ($bytes:ident, $index:ident, $cfg:ident, $provider:ident, $nutation_at:ident,
         $serving:ident, $deltat_at:ident, $bound:ident) => {
        let $bytes = super::tests::leb2_core_fixture()
            .expect("provision the hash-verified historical LEB2 core fixture");
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
            SealedLunarGlobalServing::try_new(&$serving, &$deltat_at, crate::constants::FLG_SWIEPH)
                .expect("global-lunar serving binds");
    };
}

/// Local-lunar fixture over the sealed asset at one observer triple
macro_rules! bind_lunloc18 {
    ($bytes:ident, $index:ident, $cfg:ident, $provider:ident, $nutation_at:ident,
         $serving:ident, $deltat_at:ident, $bound:ident, $lon:expr, $lat:expr, $alt:expr) => {
        let $bytes = super::tests::leb2_core_fixture()
            .expect("provision the hash-verified historical LEB2 core fixture");
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
        let $bound = SealedLunarLocalServing::try_new(
            &$serving,
            &$deltat_at,
            crate::constants::FLG_SWIEPH,
            $lon,
            $lat,
            $alt,
        )
        .expect("local-lunar serving binds");
    };
}

#[test]
fn non_finite_instants_refuse_on_all_slice18_arms() {
    bind_lunglob18!(
        bytes,
        index,
        cfg,
        provider,
        nutation_at,
        serving,
        deltat_at,
        glob
    );
    bind_lunloc18!(
        bytes2,
        index2,
        cfg2,
        provider2,
        nutation_at2,
        serving2,
        deltat_at2,
        loc,
        -99.1332,
        19.4326,
        2240.0
    );
    let _ = (
        &bytes2,
        &index2,
        &cfg2,
        &provider2,
        &nutation_at2,
        &serving2,
        &deltat_at2,
    );
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            lun_how_sealed(&loc, bad),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            lun_umbral_magnitude_sealed(&glob, bad),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            lun_penumbral_magnitude_sealed(&glob, bad),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            lun_gamma_sealed(&glob, bad),
            Err(Error::NonFiniteArgument { .. })
        ));
    }
}

#[test]
fn magnitude_and_gamma_shapes_follow_the_documented_contract() {
    bind_lunglob18!(
        bytes,
        index,
        cfg,
        provider,
        nutation_at,
        serving,
        deltat_at,
        bound
    );
    // Banked 2022-11-08 total maximum: umbral magnitude at or above
    // unity, penumbral above umbral.
    let jd_total = 2459891.957725105;
    let umb = lun_umbral_magnitude_sealed(&bound, jd_total).unwrap();
    let pen = lun_penumbral_magnitude_sealed(&bound, jd_total).unwrap();
    assert!(umb >= 1.0, "total umbral {umb}");
    assert!(pen > umb, "total penumbral {pen} umbral {umb}");
    // Banked 2024-03-25 penumbral-only maximum: no umbral phase,
    // positive penumbral magnitude.
    let jd_pen = 2460394.800551916;
    assert_eq!(lun_umbral_magnitude_sealed(&bound, jd_pen).unwrap(), 0.0);
    assert!(
        lun_penumbral_magnitude_sealed(&bound, jd_pen).unwrap() > 0.0,
        "penumbral-only magnitude"
    );
    // Off-eclipse JD: both magnitudes zero, gamma past the
    // no-eclipse threshold and ordered above the central value.
    let jd_off = 2459959.5;
    assert_eq!(lun_umbral_magnitude_sealed(&bound, jd_off).unwrap(), 0.0);
    assert_eq!(lun_penumbral_magnitude_sealed(&bound, jd_off).unwrap(), 0.0);
    let g_total = lun_gamma_sealed(&bound, jd_total).unwrap();
    let g_off = lun_gamma_sealed(&bound, jd_off).unwrap();
    assert!((0.0..1.0).contains(&g_total), "total gamma {g_total}");
    assert!(g_off > 1.5 && g_off > g_total, "off-eclipse gamma {g_off}");
}

#[test]
fn local_how_keeps_the_horizon_gate_and_shadow_words() {
    bind_lunloc18!(
        bytes,
        index,
        cfg,
        provider,
        nutation_at,
        serving,
        deltat_at,
        loc,
        -99.1332,
        19.4326,
        2240.0
    );
    bind_lunglob18!(
        bytes2,
        index2,
        cfg2,
        provider2,
        nutation_at2,
        serving2,
        deltat_at2,
        glob
    );
    let _ = (
        &bytes2,
        &index2,
        &cfg2,
        &provider2,
        &nutation_at2,
        &serving2,
        &deltat_at2,
    );
    // Banked 2022-11-08 total maximum over Mexico City: the local
    // triple stays finite, the horizon gate holds (retc zero at or
    // below the apparent horizon), and the shadow channels agree
    // bitwise with the global magnitude servings (same sealed
    // vectors through the same canonical core).
    let jd_total = 2459891.957725105;
    let (retc, attr) = lun_how_sealed(&loc, jd_total).unwrap();
    assert!(attr[4].is_finite() && attr[5].is_finite() && attr[6].is_finite());
    if attr[6] <= 0.0 {
        assert_eq!(retc, 0);
    }
    assert_eq!(
        attr[0],
        lun_umbral_magnitude_sealed(&glob, jd_total).unwrap()
    );
    assert_eq!(
        attr[1],
        lun_penumbral_magnitude_sealed(&glob, jd_total).unwrap()
    );
}
