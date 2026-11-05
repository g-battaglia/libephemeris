use super::*;
use crate::constants::{ECL_PARTIAL, ECL_PENUMBRAL, ECL_TOTAL, ECL_VISIBLE, FLG_SWIEPH};
use crate::events::{LunarGlobalEphemeris, LunarLocalEphemeris};

/// Global-lunar fixture over the sealed asset (same Espenak-Meeus
/// Delta-T shape as the solar slice-10 fixtures).
macro_rules! bind_lunglob_par {
    ($bytes:ident, $index:ident, $cfg:ident, $provider:ident, $nutation_at:ident,
         $serving:ident, $deltat_at:ident, $bound:ident) => {
        let path = "/Users/giacomo/.libephemeris/leb/base_core.leb2";
        let $bytes = match std::fs::read(path) {
            Ok(bytes) => bytes,
            Err(_) => {
                eprintln!("evt14_lunpar_tests: SKIP (missing {path})");
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
        let $bound = SealedLunarGlobalServing::try_new(&$serving, &$deltat_at, FLG_SWIEPH)
            .expect("global-lunar serving binds");
    };
}

/// Site-bound fixture over the same asset (Mexico City: the
/// 2024-lunar-eclipse cover site).
macro_rules! bind_lunloc_par {
    ($bytes:ident, $index:ident, $cfg:ident, $provider:ident, $nutation_at:ident,
         $serving:ident, $deltat_at:ident, $bound:ident) => {
        let path = "/Users/giacomo/.libephemeris/leb/base_core.leb2";
        let $bytes = match std::fs::read(path) {
            Ok(bytes) => bytes,
            Err(_) => {
                eprintln!("evt14_lunpar_tests: SKIP (missing {path})");
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
        let $bound = SealedLunarLocalServing::try_new(
            &$serving,
            &$deltat_at,
            FLG_SWIEPH,
            -99.1,
            19.4,
            2240.0,
        )
        .expect("local-lunar serving binds");
    };
}

#[test]
fn lunglob_sealed_finds_2024_event_with_ordered_phases() {
    bind_lunglob_par!(
        bytes,
        index,
        cfg,
        provider,
        nutation_at,
        serving,
        deltat_at,
        bound
    );
    // Forward walk from 2024-03-11 (UT) with mask 0 (any type).
    let (retflag, tret) = lun_eclipse_when_glob_sealed(&bound, 2460390.0, 0, false)
        .expect("2024 lunar eclipse serves");
    // Exactly one classified type bit.
    let type_bits = retflag & (ECL_TOTAL | ECL_PARTIAL | ECL_PENUMBRAL);
    assert!(
        type_bits == ECL_TOTAL || type_bits == ECL_PARTIAL || type_bits == ECL_PENUMBRAL,
        "single lunar type bit: {retflag:#x}"
    );
    // 2024 cover: the March-2024 penumbral eclipse maximum falls
    // inside this window (walk must not skip to a later lunation).
    assert!(
        tret[0] > 2460390.0 && tret[0] < 2460420.0,
        "2024 maximum: {}",
        tret[0]
    );
    // Penumbral contacts bracket the maximum; reserved slots stay
    // 0.0.
    assert!(tret[6] > 0.0 && tret[6] < tret[0] && tret[0] < tret[7]);
    assert_eq!((tret[1], tret[8], tret[9]), (0.0, 0.0, 0.0));
    // Phase-slot/type-bit consistency.
    if type_bits == ECL_TOTAL {
        assert!(tret[2] > 0.0 && tret[3] > tret[2]);
        assert!(tret[4] > 0.0 && tret[5] > tret[4]);
        assert!(tret[6] <= tret[4] && tret[5] <= tret[7]);
        assert!(tret[2] <= tret[4] && tret[5] <= tret[3]);
    } else if type_bits == ECL_PARTIAL {
        assert!(tret[2] > 0.0 && tret[3] > tret[2]);
        assert_eq!((tret[4], tret[5]), (0.0, 0.0));
    } else {
        assert_eq!((tret[2], tret[3], tret[4], tret[5]), (0.0, 0.0, 0.0, 0.0));
    }
    // The shadow core at the reported maximum classifies the same
    // type (walk/core agreement over identical sealed words).
    let (shadow_type, _) = bound.lun_shadow(tret[0]).expect("shadow serves");
    assert_eq!(shadow_type, type_bits);
}

#[test]
fn lunloc_sealed_how_attr_pins_contact_composition() {
    bind_lunloc_par!(
        bytes,
        index,
        cfg,
        provider,
        nutation_at,
        serving,
        deltat_at,
        bound
    );
    let (retflag, tret, attr) =
        lun_eclipse_when_loc_sealed(&bound, 2460390.0, false).expect("local serves");
    assert_ne!(
        retflag & ECL_VISIBLE,
        0,
        "2024 event visible from Mexico City"
    );
    assert!(tret[0] > 2460390.0 && tret[0] < 2460420.0);
    // The walk returns the `how_attr` words at its reported maximum
    // bitwise.
    let (retc, direct) = bound.how_attr(tret[0]).expect("how serves");
    assert_eq!(retc, retflag & (ECL_TOTAL | ECL_PARTIAL | ECL_PENUMBRAL));
    assert_eq!(direct, attr);
    // Contact composition: lunar saros pair filled at a real
    // eclipse, sealed Moon azalt triple finite in attr[4..6].
    assert_ne!(attr[9], 0.0, "lunar saros series serves");
    assert!(attr[4].is_finite() && attr[5].is_finite() && attr[6].is_finite());
    // Horizon gate: the found
    // maximum is up, and a below-horizon instant zeroes retc while
    // the attributes stay filled (lunar fill-preserving rule, unlike
    // the solar slot zeroing reviewed in slice 12). Note the gate is
    // one-way past the contacts: outside the shadow (past P4) retc
    // is 0 with the Moon still up, so the probe pins
    // `apparent <= 0 ⟹ retc == 0`, not the converse.
    assert!(attr[6] > 0.0);
    let mut down_found = false;
    let mut probe = tret[0];
    for _ in 0..30 {
        probe += 0.05;
        let (rc, probe_attr) = bound.how_attr(probe).expect("probe serves");
        assert!(
            probe_attr[4].is_finite() && probe_attr[5].is_finite() && probe_attr[6].is_finite(),
            "attributes stay filled at {probe}"
        );
        if probe_attr[6] <= 0.0 {
            assert_eq!(rc, 0, "horizon zeroes retc at {probe}");
            down_found = true;
            break;
        }
    }
    assert!(down_found, "a below-horizon probe inside 1.5 d");
}
