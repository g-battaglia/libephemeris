use super::*;

/// Read the released core-data fixture, refusing unrelated coefficient files.
/// `LIBEPHEMERIS_TEST_LEB2` selects an explicit file; otherwise use
/// the owner-selected core under `~/.libephemeris/leb`. Absence is an optional-data skip,
/// but a present file with the wrong SHA-256 is a fixture error, not a numerical comparison.
fn released_core_fixture() -> Option<Vec<u8>> {
    use sha2::{Digest, Sha256};
    let path = std::env::var_os("LIBEPHEMERIS_TEST_LEB2")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(crate::test_data::core_leb_path);
    let bytes = match std::fs::read(&path) {
        Ok(bytes) => bytes,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
            eprintln!("Besselian fixture: SKIP (missing {})", path.display());
            return None;
        }
        Err(error) => panic!("cannot read fixture {}: {error}", path.display()),
    };
    assert!(
        crate::test_data::admitted_core_hash(
            &(Sha256::digest(&bytes)
                .iter()
                .map(|byte| format!("{byte:02x}"))
                .collect::<String>())
        ),
        "Besselian word fixture requires the released core coefficient file: {}",
        path.display()
    );
    Some(bytes)
}

/// Global-solar fixture over the sealed asset (no site: the
/// fundamental-plane words are geocentric; Espenak-Meeus Delta-T).
/// `SKIP`s (never fails) when the asset is absent.
macro_rules! bind_solglob22 {
    ($provider:ident, $nutation_at:ident, $serving:ident, $deltat_at:ident, $bound:ident) => {
        let Some(bytes) = released_core_fixture() else {
            return;
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
fn besselian_non_finite_instants_refuse() {
    bind_solglob22!(provider, nutation_at, serving, deltat_at, bound);
    for word in [
        calc_besselian_x_sealed(&bound, f64::NAN),
        calc_besselian_y_sealed(&bound, f64::INFINITY),
        calc_besselian_d_sealed(&bound, f64::NEG_INFINITY),
        calc_besselian_l1_sealed(&bound, f64::NAN),
        calc_besselian_l2_sealed(&bound, f64::NAN),
        calc_besselian_mu_sealed(&bound, f64::NAN),
        calc_besselian_dx_dt_sealed(&bound, f64::NAN),
        calc_besselian_dy_dt_sealed(&bound, f64::NAN),
        calc_besselian_dd_dt_sealed(&bound, f64::NAN),
        calc_besselian_dl1_dt_sealed(&bound, f64::NAN),
        calc_besselian_dl2_dt_sealed(&bound, f64::NAN),
        calc_besselian_dmu_dt_sealed(&bound, f64::NAN),
        besselian_elements_sealed(&bound, f64::NAN).map(|e| e.x()),
    ] {
        assert!(
            matches!(word, Err(Error::NonFiniteArgument { .. })),
            "non-finite instant refuses, got {word:?}"
        );
    }
}

#[test]
fn besselian_total_maximum_words_have_eclipse_shape() {
    bind_solglob22!(provider, nutation_at, serving, deltat_at, bound);
    // Banked 2024-04-08 total global maximum (solglob reference values).
    let jd_max = 2460409.2619832945;
    let x = calc_besselian_x_sealed(&bound, jd_max).expect("x solves");
    let y = calc_besselian_y_sealed(&bound, jd_max).expect("y solves");
    let d = calc_besselian_d_sealed(&bound, jd_max).expect("d solves");
    let mu = calc_besselian_mu_sealed(&bound, jd_max).expect("mu solves");
    let l1 = calc_besselian_l1_sealed(&bound, jd_max).expect("l1 solves");
    let l2 = calc_besselian_l2_sealed(&bound, jd_max).expect("l2 solves");
    assert!(
        x.hypot(y) < 1.0,
        "shadow axis pierces inside the disc at maximum, got x={x} y={y}"
    );
    assert!(
        (5.0..10.0).contains(&d),
        "axis declination tracks the April Sun, got {d}"
    );
    assert!(
        (0.0..360.0).contains(&mu),
        "hour angle normalizes to [0, 360), got {mu}"
    );
    assert!(
        l1 > 0.5 && l1 < 0.6,
        "penumbral radius near the April-2024 value, got {l1}"
    );
    assert!(l2 < 0.0, "total eclipse carries a negative l2, got {l2}");
    let dx = calc_besselian_dx_dt_sealed(&bound, jd_max).expect("dx/dt solves");
    let dmu = calc_besselian_dmu_dt_sealed(&bound, jd_max).expect("dmu/dt solves");
    assert!(
        dx > 0.4 && dx < 0.6,
        "shadow runs eastward near half a radius per hour, got {dx}"
    );
    assert!(
        (dmu - 15.04).abs() < 0.1,
        "mu rate tracks sidereal rotation, got {dmu}"
    );
}

#[test]
fn besselian_mu_uses_default_timescale_nutation_epoch() {
    bind_solglob22!(provider, nutation_at, serving, deltat_at, bound);
    let mu = calc_besselian_mu_sealed(&bound, 2460232.249616082).unwrap();
    // Pinned v3.2.1 value over the same identified coefficients and clock.
    let expected = f64::from_bits(0x40575736b35bcc34);
    let angular_error = (mu - expected + 180.0).rem_euclid(360.0) - 180.0;
    assert!(
        angular_error.abs() * 3600.0 <= 0.001,
        "Besselian hour angle {mu} vs pinned {expected}"
    );
}

#[test]
fn besselian_bundle_matches_the_twelve_scalars() {
    bind_solglob22!(provider, nutation_at, serving, deltat_at, bound);
    // Banked 2023-10-14 annular global maximum (solcontact reference values).
    let jd_max = 2460232.249616082;
    let bundle = besselian_elements_sealed(&bound, jd_max).expect("bundle solves");
    assert_eq!(bundle.t0(), jd_max);
    for (name, bundled, scalar) in [
        ("x", bundle.x(), calc_besselian_x_sealed(&bound, jd_max)),
        ("y", bundle.y(), calc_besselian_y_sealed(&bound, jd_max)),
        ("d", bundle.d(), calc_besselian_d_sealed(&bound, jd_max)),
        ("l1", bundle.l1(), calc_besselian_l1_sealed(&bound, jd_max)),
        ("l2", bundle.l2(), calc_besselian_l2_sealed(&bound, jd_max)),
        ("mu", bundle.mu(), calc_besselian_mu_sealed(&bound, jd_max)),
        (
            "dx_dt",
            bundle.dx_dt(),
            calc_besselian_dx_dt_sealed(&bound, jd_max),
        ),
        (
            "dy_dt",
            bundle.dy_dt(),
            calc_besselian_dy_dt_sealed(&bound, jd_max),
        ),
        (
            "dd_dt",
            bundle.dd_dt(),
            calc_besselian_dd_dt_sealed(&bound, jd_max),
        ),
        (
            "dl1_dt",
            bundle.dl1_dt(),
            calc_besselian_dl1_dt_sealed(&bound, jd_max),
        ),
        (
            "dl2_dt",
            bundle.dl2_dt(),
            calc_besselian_dl2_dt_sealed(&bound, jd_max),
        ),
        (
            "dmu_dt",
            bundle.dmu_dt(),
            calc_besselian_dmu_dt_sealed(&bound, jd_max),
        ),
    ] {
        assert_eq!(
            bundled,
            scalar.expect("scalar solves"),
            "{name} matches the scalar serving bitwise"
        );
    }
    assert!(
        bundle.l2() > 0.0,
        "annular eclipse carries a positive l2, got {}",
        bundle.l2()
    );
}
