use super::*;

/// Global-solar fixture over the sealed asset (same shape as the
/// slice-9 fixture: sealed core, Espenak-Meeus Delta-T).
macro_rules! bind_solglob10 {
    ($bytes:ident, $index:ident, $cfg:ident, $provider:ident, $nutation_at:ident,
         $serving:ident, $deltat_at:ident, $bound:ident) => {
        let path = "/Users/giacomo/.libephemeris/leb/base_core.leb2";
        let $bytes = match std::fs::read(path) {
            Ok(bytes) => bytes,
            Err(_) => {
                eprintln!("evt14_slice10_tests: SKIP (missing {path})");
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

/// Site-bound fixture over the same asset (Dallas by default: the
/// 2024-04-08 total-eclipse cover).
macro_rules! bind_solloc10 {
    ($bytes:ident, $index:ident, $cfg:ident, $provider:ident, $nutation_at:ident,
         $serving:ident, $deltat_at:ident, $bound:ident, $lon:expr, $lat:expr) => {
        let path = "/Users/giacomo/.libephemeris/leb/base_core.leb2";
        let $bytes = match std::fs::read(path) {
            Ok(bytes) => bytes,
            Err(_) => {
                eprintln!("evt14_slice10_tests: SKIP (missing {path})");
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
        let $bound = SealedSolarLocalServing::try_new(
            &$serving,
            &$deltat_at,
            crate::constants::FLG_SWIEPH,
            $lon,
            $lat,
            0.0,
        )
        .expect("local-solar serving binds");
    };
}

#[test]
fn sealed_path_width_covers_2024_total_eclipse() {
    bind_solglob10!(
        bytes,
        index,
        cfg,
        provider,
        nutation_at,
        serving,
        deltat_at,
        bound
    );
    // 2024-04-08 greatest-eclipse neighbourhood (UT).
    let jd = 2460409.26;
    let (retflag, clon, clat) = bound.path_width_where(jd).expect("where serves");
    assert_ne!(retflag & ECL_CENTRAL, 0, "greatest eclipse is central");
    // Central-line width is a positive sub-cap word.
    let central = sol_eclipse_path_width_sealed(&bound, jd, None, None).expect("width serves");
    assert!(central > 0.0 && central <= PATH_WIDTH_MAX_KM);
    // The explicit site at the shadow center reads the same words.
    let explicit = sol_eclipse_path_width_sealed(&bound, jd, Some(clat), Some(clon))
        .expect("explicit width serves");
    assert_eq!(explicit.to_bits(), central.to_bits());
    // A half-provided site falls back to the central line.
    let half = sol_eclipse_path_width_sealed(&bound, jd, Some(clat), None)
        .expect("half-site width serves");
    assert_eq!(half.to_bits(), central.to_bits());
    // Far from any eclipse there is no central shadow: 0.0.
    let quiet =
        sol_eclipse_path_width_sealed(&bound, 2460400.5, None, None).expect("quiet width serves");
    assert_eq!(quiet, 0.0);
    // Non-finite epochs are refused, never walked.
    assert!(sol_eclipse_path_width_sealed(&bound, f64::NAN, None, None).is_err());
    // A non-finite explicit site degrades to 0.0 through the
    // release non-range failure arm (FIX-C31(a)), never an error.
    let nan_site = sol_eclipse_path_width_sealed(&bound, jd, Some(f64::NAN), Some(clon))
        .expect("NaN site serves");
    assert_eq!(nan_site, 0.0);
}

#[test]
fn sealed_duration_covers_2024_total_eclipse() {
    bind_solglob10!(
        bytes,
        index,
        cfg,
        provider,
        nutation_at,
        serving,
        deltat_at,
        bound
    );
    // Maximum central duration of the 2024-04-08 total eclipse
    // .
    let minutes = sol_eclipse_duration_min_sealed(&bound, 2460409.26).expect("duration serves");
    assert!(
        minutes > 0.0 && minutes < 10.0,
        "2024 duration minutes = {minutes}"
    );
    // The duration contacts ride the local-phase core: C2 precedes
    // C3 at the shadow center.
    let (retflag, clon, clat) = bound.duration_where(2460409.26).expect("where serves");
    assert_ne!(retflag & ECL_CENTRAL, 0);
    let (c2, c3) = bound
        .duration_local_contacts(2460409.26, clat, clon)
        .expect("contacts serve");
    assert!(c2 > 0.0 && c3 > c2);
    assert!(((c3 - c2) * MINUTES_PER_DAY - minutes).abs() < 1e-9);
    // Far from any eclipse: 0.0. Non-finite epochs are refused.
    let quiet = sol_eclipse_duration_min_sealed(&bound, 2460400.5).expect("quiet serves");
    assert_eq!(quiet, 0.0);
    assert!(sol_eclipse_duration_min_sealed(&bound, f64::NAN).is_err());
}

#[test]
fn sealed_how_details_covers_dallas_2024() {
    bind_solloc10!(
        bytes,
        index,
        cfg,
        provider,
        nutation_at,
        serving,
        deltat_at,
        bound,
        -96.797,
        32.7767
    );
    // Dallas is in partial phase at 18:14 UT; the composition finds
    // the total local maximum (~18:29 UT).
    let out = sol_eclipse_how_details_sealed(&bound, 2460409.26).expect("details serve");
    assert!(out.is_visible, "Dallas sees the 2024 eclipse");
    assert!(out.is_total, "Dallas maximum is total");
    assert!(!out.is_partial && !out.is_annular);
    assert!(out.jd_c1 > 0.0);
    assert!(out.jd_c1 < out.jd_c2 && out.jd_c2 < out.jd_c3 && out.jd_c3 < out.jd_c4);
    assert!(out.duration_total_minutes > 0.0);
    assert!(out.duration_partial_minutes > out.duration_total_minutes);
    assert!(out.max_magnitude >= 1.0);
    assert!((out.max_obscuration - 1.0).abs() < 1e-6);
    // The provider words travel with the composition: instant slots
    // come from the initial how word, refined slots from the sealed
    // acquisitions.
    assert!(out.sun_angular_radius > 0.0 && out.moon_angular_radius > 0.0);
    assert!(out.sun_alt_max > 0.0);
    // Slice 12: the contact blocks compose their alt-az/position-angle
    // tails on the sealed path. Every solved contact carries a finite
    // daytime Sun position and a folded position angle.
    for (n, alt, az, pa) in [
        ("c1", out.sun_alt_c1, out.sun_az_c1, out.position_angle_c1),
        ("c2", out.sun_alt_c2, out.sun_az_c2, out.position_angle_c2),
        ("c3", out.sun_alt_c3, out.sun_az_c3, out.position_angle_c3),
        ("c4", out.sun_alt_c4, out.sun_az_c4, out.position_angle_c4),
    ] {
        assert!(alt.is_finite(), "sun_alt_{n} finite");
        assert!(az.is_finite(), "sun_az_{n} finite");
        assert!(alt > 0.0, "sun_alt_{n} above horizon = {alt}");
        assert!(
            (0.0..360.0).contains(&pa) || pa == 360.0,
            "position_angle_{n} folded = {pa}"
        );
    }
    // Far from any eclipse the composition takes the early return:
    // type 0 with zero contacts.
    let quiet = sol_eclipse_how_details_sealed(&bound, 2460400.5).expect("quiet serves");
    assert_eq!(quiet.eclipse_type, 0);
    assert_eq!(quiet.jd_c1, 0.0);
    assert_eq!(quiet.jd_c4, 0.0);
    // Non-finite epochs are refused.
    assert!(sol_eclipse_how_details_sealed(&bound, f64::NAN).is_err());
}

#[test]
fn sealed_how_details_partial_keeps_c1_c4_tails_only() {
    bind_solloc10!(
        bytes,
        index,
        cfg,
        provider,
        nutation_at,
        serving,
        deltat_at,
        bound,
        -74.006,
        40.7128
    );
    // New York sees the 2024-04-08 eclipse as a daytime partial:
    // outer contacts solve with their tails, inner contacts stay
    // `0.0`.
    let out = sol_eclipse_how_details_sealed(&bound, 2460409.30).expect("details serve");
    assert!(out.is_visible, "New York sees the 2024 eclipse");
    assert!(out.is_partial, "New York maximum is partial");
    assert!(!out.is_total && !out.is_annular);
    assert!(out.jd_c1 > 0.0 && out.jd_c4 > out.jd_c1);
    assert_eq!(out.jd_c2, 0.0);
    assert_eq!(out.jd_c3, 0.0);
    assert!(out.duration_partial_minutes > 0.0);
    assert_eq!(out.duration_total_minutes, 0.0);
    assert!(out.sun_alt_max > 0.0);
    // The solved outer blocks carry populated tails: a degraded
    // acquisition would leave the release `0.0` slots behind.
    for (n, alt, az, pa) in [
        ("c1", out.sun_alt_c1, out.sun_az_c1, out.position_angle_c1),
        ("c4", out.sun_alt_c4, out.sun_az_c4, out.position_angle_c4),
    ] {
        assert!(alt.is_finite(), "sun_alt_{n} finite");
        assert!(az.is_finite(), "sun_az_{n} finite");
        assert!(alt > 0.0, "sun_alt_{n} above horizon = {alt}");
        assert!(
            (0.0..360.0).contains(&pa),
            "position_angle_{n} folded = {pa}"
        );
    }
    // First and fourth contacts touch opposite limbs.
    assert!(
        (out.position_angle_c1 - out.position_angle_c4).abs() > 1.0,
        "C1/C4 position angles differ ({} vs {})",
        out.position_angle_c1,
        out.position_angle_c4
    );
    // Untouched inner tails keep the release `0.0` sentinels.
    assert_eq!(out.sun_alt_c2, 0.0);
    assert_eq!(out.sun_alt_c3, 0.0);
    assert_eq!(out.position_angle_c2, 0.0);
    assert_eq!(out.position_angle_c3, 0.0);
}
