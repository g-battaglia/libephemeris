use super::*;

/// Global-solar fixture over the sealed asset (same shape as the
/// slice-9 `bind_solglob!`: sealed core, Espenak-Meeus Delta-T).
/// `SKIP`s (never fails) when the asset is absent.
macro_rules! bind_solglob_contact {
    ($provider:ident, $nutation_at:ident, $serving:ident, $deltat_at:ident, $bound:ident) => {
        let path = "/Users/giacomo/.libephemeris/leb/base_core.leb2";
        let bytes = match std::fs::read(path) {
            Ok(bytes) => bytes,
            Err(_) => {
                eprintln!("evt14_slice20_solcontact_tests: SKIP (missing {path})");
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
fn sol_contact_non_finite_maxima_refuse() {
    bind_solglob_contact!(provider, nutation_at, serving, deltat_at, bound);
    for word in [
        sol_contact_c1_sealed(&bound, f64::INFINITY),
        sol_contact_c2_sealed(&bound, f64::NAN),
        sol_contact_c3_sealed(&bound, f64::NEG_INFINITY),
        sol_contact_c4_sealed(&bound, f64::NAN),
    ] {
        assert!(
            matches!(word, Err(Error::NonFiniteArgument { .. })),
            "non-finite jd_max refuses, got {word:?}"
        );
    }
}

#[test]
fn sol_contact_total_orders_c1_c2_max_c3_c4() {
    bind_solglob_contact!(provider, nutation_at, serving, deltat_at, bound);
    // Banked 2024-04-08 total global maximum (solglob reference values).
    let jd_max = 2460409.2619832945;
    let c1 = sol_contact_c1_sealed(&bound, jd_max).expect("C1 solves");
    let c2 = sol_contact_c2_sealed(&bound, jd_max).expect("C2 solves");
    let c3 = sol_contact_c3_sealed(&bound, jd_max).expect("C3 solves");
    let c4 = sol_contact_c4_sealed(&bound, jd_max).expect("C4 solves");
    assert!(
        c1 > 0.0 && c2 > 0.0 && c3 > 0.0 && c4 > 0.0,
        "all contacts solve"
    );
    assert!(
        c1 < c2 && c2 < jd_max && jd_max < c3 && c3 < c4,
        "C1 < C2 < max < C3 < C4, got {c1} {c2} {c3} {c4}"
    );
}

#[test]
fn sol_contact_partial_has_no_inner_contacts() {
    bind_solglob_contact!(provider, nutation_at, serving, deltat_at, bound);
    // Banked 2022-10-25 partial global maximum (retflag 18, non-central).
    let jd_max = 2459877.958411902;
    let c1 = sol_contact_c1_sealed(&bound, jd_max).expect("C1 solves");
    let c2 = sol_contact_c2_sealed(&bound, jd_max).expect("C2 degrades");
    let c3 = sol_contact_c3_sealed(&bound, jd_max).expect("C3 degrades");
    let c4 = sol_contact_c4_sealed(&bound, jd_max).expect("C4 solves");
    assert!(
        c1 > 0.0 && c4 > 0.0 && c1 < jd_max && jd_max < c4,
        "outer contacts solve"
    );
    assert_eq!(c2, 0.0, "partial C2 is 0.0");
    assert_eq!(c3, 0.0, "partial C3 is 0.0");
}

#[test]
fn sol_contact_quiet_maximum_is_all_zeros() {
    bind_solglob_contact!(provider, nutation_at, serving, deltat_at, bound);
    let jd_max = 2460140.0;
    assert_eq!(
        sol_contact_c1_sealed(&bound, jd_max).expect("C1 degrades"),
        0.0
    );
    assert_eq!(
        sol_contact_c2_sealed(&bound, jd_max).expect("C2 degrades"),
        0.0
    );
    assert_eq!(
        sol_contact_c3_sealed(&bound, jd_max).expect("C3 degrades"),
        0.0
    );
    assert_eq!(
        sol_contact_c4_sealed(&bound, jd_max).expect("C4 degrades"),
        0.0
    );
}
