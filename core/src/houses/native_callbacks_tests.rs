//! Recorded provider controls for consumed sidereal house branches.
use super::*;
use crate::constants::*;
use std::cell::RefCell;

fn selection(mode: i32) -> crate::sidereal::SiderealSelection {
    let mut s = crate::sidereal::SiderealSelection::unset();
    s.set_sidereal_mode(mode, 2440587.5, 24.0);
    s
}

fn refusal() -> Error {
    Error::MathDomain { value: -1.0 }
}

#[test]
fn native_callbacks_tropical_does_not_acquire_sun_or_sidereal() {
    let fail_epoch = || Err(refusal());
    let sid = HousesSiderealCallbacks {
        ayanamsha_ex_at: &|_, _| Err(refusal()),
        mean_ut_plain_at: &|_| Err(refusal()),
        nutation_at: &|_| Err(refusal()),
        user_t0_tt_at: &fail_epoch,
        valens_epoch_tt_at: &fail_epoch,
    };
    let jd = HousesJdProviders {
        frame_at: &|_| Ok((100.0, 23.44)),
        sun_declination_at: &|_, _| Err(refusal()),
        deltat_days_at: &|_| Err(refusal()),
    };
    let actual = houses_ex_via_callbacks(
        2451545.125,
        52.52,
        13.405,
        'P',
        0,
        selection(SIDM_USER | SIDBIT_USER_UT | SIDBIT_ECL_T0),
        &jd,
        &sid,
    )
    .unwrap();
    let expected = houses_jd_via(2451545.125, 52.52, 13.405, 'P', 0, &jd).unwrap();
    assert_eq!(actual, expected);
    assert!(
        actual
            .0
            .iter()
            .chain(actual.1.iter())
            .all(|x| x.is_finite())
    );
}

#[test]
fn native_callbacks_plain_sidereal_masks_toggles_after_tropical_seed() {
    let calls = RefCell::new(Vec::new());
    let frame = |_| {
        calls.borrow_mut().push("frame");
        Ok((100.0, 23.44))
    };
    let aya = |_, flags| {
        calls.borrow_mut().push("aya");
        assert_eq!(flags, FLG_TRUEPOS | FLG_NOABERR | FLG_NOGDEFL);
        Ok((2, 24.0))
    };
    let fail_epoch = || Err(refusal());
    let sid = HousesSiderealCallbacks {
        ayanamsha_ex_at: &aya,
        mean_ut_plain_at: &|_| Err(refusal()),
        nutation_at: &|_| Err(refusal()),
        user_t0_tt_at: &fail_epoch,
        valens_epoch_tt_at: &fail_epoch,
    };
    let jd = HousesJdProviders {
        frame_at: &frame,
        sun_declination_at: &|_, _| Err(refusal()),
        deltat_days_at: &|_| Err(refusal()),
    };
    houses_ex_via_callbacks(
        2451545.125,
        40.0,
        10.0,
        'P',
        FLG_SIDEREAL | FLG_TRUEPOS | FLG_NOABERR | FLG_NOGDEFL | FLG_SPEED,
        selection(SIDM_LAHIRI | SIDBIT_USER_UT),
        &jd,
        &sid,
    )
    .unwrap();
    assert_eq!(*calls.borrow(), ["frame", "aya"]);
}

fn projected_order(mode: i32, expected: &[&str], mean_epoch: Option<f64>) {
    let calls = RefCell::new(Vec::new());
    let frame = |_| {
        calls.borrow_mut().push("frame");
        Ok((100.0, 23.44))
    };
    let delta = |_| {
        calls.borrow_mut().push("delta");
        Ok(0.0008)
    };
    let nut = |_| {
        calls.borrow_mut().push("nut");
        Ok((0.0, 0.0))
    };
    let user = || {
        calls.borrow_mut().push("user");
        Ok(2440587.5)
    };
    let valens = || {
        calls.borrow_mut().push("valens");
        Err(refusal())
    };
    let mean = |jd| {
        calls.borrow_mut().push("mean");
        assert_eq!(Some(jd), mean_epoch);
        Err(refusal())
    };
    let sid = HousesSiderealCallbacks {
        ayanamsha_ex_at: &|_, _| Err(refusal()),
        mean_ut_plain_at: &mean,
        nutation_at: &nut,
        user_t0_tt_at: &user,
        valens_epoch_tt_at: &valens,
    };
    let jd = HousesJdProviders {
        frame_at: &frame,
        sun_declination_at: &|_, _| Err(refusal()),
        deltat_days_at: &delta,
    };
    assert!(matches!(
        houses_ex_via_callbacks(
            2451545.125,
            40.0,
            10.0,
            'P',
            FLG_SIDEREAL,
            selection(mode),
            &jd,
            &sid
        ),
        Err(Error::MathDomain { .. })
    ));
    assert_eq!(*calls.borrow(), expected);
}

#[test]
fn native_callbacks_ecl_user_epoch_follows_frame_and_precedes_zero_point() {
    projected_order(
        SIDM_USER | SIDBIT_ECL_T0 | SIDBIT_USER_UT,
        &["frame", "delta", "nut", "user", "mean"],
        Some(2440587.5),
    );
}

#[test]
fn native_callbacks_ecl_valens_refusal_precedes_zero_point() {
    projected_order(
        SIDM_VALENS_MOON | SIDBIT_ECL_T0,
        &["frame", "delta", "nut", "valens"],
        None,
    );
}

#[test]
fn native_callbacks_ssy_never_reads_user_or_valens_plane_epochs() {
    projected_order(
        SIDM_USER | SIDBIT_SSY_PLANE | SIDBIT_USER_UT,
        &["frame", "delta", "nut", "mean"],
        Some(crate::sidereal::J2000_TT),
    );
}

#[test]
fn native_callbacks_classical_ecl_reads_table_without_model_epochs() {
    projected_order(
        SIDM_LAHIRI | SIDBIT_ECL_T0 | SIDBIT_USER_UT,
        &["frame", "delta", "nut", "mean"],
        Some(crate::sidereal::ecl_t0_epoch_jd(
            SIDM_LAHIRI,
            2440587.5,
            1775845.5,
        )),
    );
}

#[test]
fn native_callbacks_frame_refusal_precedes_sidereal_epochs() {
    let calls = RefCell::new(Vec::new());
    let frame = |_| {
        calls.borrow_mut().push("frame");
        Err(refusal())
    };
    let epoch = || {
        calls.borrow_mut().push("epoch");
        Err(refusal())
    };
    let sid = HousesSiderealCallbacks {
        ayanamsha_ex_at: &|_, _| Err(refusal()),
        mean_ut_plain_at: &|_| Err(refusal()),
        nutation_at: &|_| Err(refusal()),
        user_t0_tt_at: &epoch,
        valens_epoch_tt_at: &epoch,
    };
    let jd = HousesJdProviders {
        frame_at: &frame,
        sun_declination_at: &|_, _| Err(refusal()),
        deltat_days_at: &|_| Err(refusal()),
    };
    assert!(
        houses_ex_via_callbacks(
            2451545.125,
            40.0,
            10.0,
            'P',
            FLG_SIDEREAL,
            selection(SIDM_USER | SIDBIT_ECL_T0),
            &jd,
            &sid
        )
        .is_err()
    );
    assert_eq!(*calls.borrow(), ["frame"]);
}

#[test]
fn native_callbacks_sunshine_sun_is_consumed_before_plain_ayanamsha() {
    let calls = RefCell::new(Vec::new());
    let frame = |_| {
        calls.borrow_mut().push("frame");
        Ok((100.0, 23.44))
    };
    let sun = |_, flags| {
        calls.borrow_mut().push("sun");
        assert_eq!(flags, FLG_JPLEPH);
        Ok(11.25)
    };
    let aya = |_, _| {
        calls.borrow_mut().push("aya");
        Ok((2, 24.0))
    };
    let fail_epoch = || Err(refusal());
    let sid = HousesSiderealCallbacks {
        ayanamsha_ex_at: &aya,
        mean_ut_plain_at: &|_| Err(refusal()),
        nutation_at: &|_| Err(refusal()),
        user_t0_tt_at: &fail_epoch,
        valens_epoch_tt_at: &fail_epoch,
    };
    let jd = HousesJdProviders {
        frame_at: &frame,
        sun_declination_at: &sun,
        deltat_days_at: &|_| Err(refusal()),
    };
    houses_ex_via_callbacks(
        2451545.125,
        40.0,
        10.0,
        'I',
        FLG_SIDEREAL | FLG_JPLEPH,
        selection(SIDM_LAHIRI),
        &jd,
        &sid,
    )
    .unwrap();
    assert_eq!(*calls.borrow(), ["frame", "sun", "aya"]);
}

#[test]
fn native_callbacks_eager_api_keeps_identical_words() {
    let aya = |_, _| Ok((2, 24.0));
    let mean = |_| Ok(24.0);
    let nut = |_| Ok((0.0, 0.0));
    let user = || Ok(2440587.5);
    let valens = || Ok(1775845.5);
    let callbacks = HousesSiderealCallbacks {
        ayanamsha_ex_at: &aya,
        mean_ut_plain_at: &mean,
        nutation_at: &nut,
        user_t0_tt_at: &user,
        valens_epoch_tt_at: &valens,
    };
    let eager = HousesSiderealProviders {
        ayanamsha_ex_at: &aya,
        mean_ut_plain_at: &mean,
        nutation_at: &nut,
        user_t0_tt: 2440587.5,
        valens_epoch_tt: 1775845.5,
    };
    let jd = HousesJdProviders {
        frame_at: &|_| Ok((100.0, 23.44)),
        sun_declination_at: &|_, _| Ok(11.25),
        deltat_days_at: &|_| Ok(0.0008),
    };
    for mode in [
        SIDM_LAHIRI,
        SIDM_J2000,
        SIDM_USER | SIDBIT_ECL_T0,
        SIDM_USER | SIDBIT_SSY_PLANE,
    ] {
        let a = houses_ex_via(
            2451545.125,
            40.0,
            10.0,
            'P',
            FLG_SIDEREAL,
            selection(mode),
            &jd,
            &eager,
        )
        .unwrap();
        let b = houses_ex_via_callbacks(
            2451545.125,
            40.0,
            10.0,
            'P',
            FLG_SIDEREAL,
            selection(mode),
            &jd,
            &callbacks,
        )
        .unwrap();
        assert_eq!(
            a.0.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            b.0.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
        );
        assert_eq!(a.1.map(f64::to_bits), b.1.map(f64::to_bits));
    }
}
