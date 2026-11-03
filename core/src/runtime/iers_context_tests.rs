//! Model-clock consumption and lifecycle controls with owned explicit rows.
use super::*;
use crate::iers::{IersAcquisitionRequest, IersDataSnapshot, ResolvedIersAcquisition};
use crate::time::{DeltaTModel, deltat_with_config};
const ROWS: &str = "1973 2 1 43.4724\n1973 3 1 43.5648\n";
fn jd() -> f64 {
    crate::time::mjd_to_jd(crate::time::calendar_to_mjd(1973, 2, 1))
        .expect("finite monthly fixture epoch")
}
fn rows_context() -> EphemerisContext {
    let mut ctx = EphemerisContext::new();
    ctx.attach_iers_rows(IersDataSnapshot::from_texts(None, None, Some(ROWS)));
    ctx
}
fn must_not_resolve(_: &IersAcquisitionRequest) -> Result<ResolvedIersAcquisition> {
    panic!("unused IERS acquisition");
}

#[test]
fn user_override_and_disabled_iers_never_resolve_acquisition() {
    let mut ctx = EphemerisContext::new();
    ctx.set_iers_boundary_resolver(Some(must_not_resolve));
    let cfg = ctx.delta_t_config(DeltaTModel::EspenakMeeus, false);
    assert!(cfg.iers_lookup_seconds.is_none());
    assert!(deltat_with_config(jd(), &cfg).is_ok());
    ctx.set_user_delta_t(Some(0.25));
    assert_eq!(
        deltat_with_config(
            f64::NAN,
            &ctx.delta_t_config(DeltaTModel::EspenakMeeus, true)
        )
        .unwrap(),
        0.25
    );
}

#[test]
fn nonfinite_epoch_refuses_before_iers_callback() {
    let mut ctx = EphemerisContext::new();
    ctx.set_iers_boundary_resolver(Some(must_not_resolve));
    assert!(matches!(
        deltat_with_config(
            f64::INFINITY,
            &ctx.delta_t_config(DeltaTModel::EspenakMeeus, true)
        ),
        Err(Error::NonFiniteArgument { .. })
    ));
}

#[test]
fn observed_rows_are_consumed_before_selected_model() {
    let ctx = rows_context();
    let mut cfg = ctx.delta_t_config(DeltaTModel::Smh2016, true);
    cfg.smh2016_lookup_seconds = None;
    assert_eq!(
        deltat_with_config(jd(), &cfg).unwrap().to_bits(),
        (43.4724_f64 / 86400.0).to_bits()
    );
    assert!(matches!(
        deltat_with_config(jd() - 1.0, &cfg),
        Err(Error::DeltaTModelUnavailable { .. })
    ));
}

#[test]
fn reset_preserves_rows_close_releases_current_rows_old_clone_retains_them() {
    let mut ctx = rows_context();
    let serving = ctx.clone();
    ctx.reset_session();
    assert_eq!(
        ctx.delta_t_config(DeltaTModel::EspenakMeeus, true)
            .iers_lookup_seconds
            .unwrap()(jd()),
        Some(43.4724)
    );
    ctx.close();
    assert_eq!(
        ctx.delta_t_config(DeltaTModel::EspenakMeeus, true)
            .iers_lookup_seconds
            .unwrap()(jd()),
        None
    );
    assert_eq!(
        serving
            .delta_t_config(DeltaTModel::EspenakMeeus, true)
            .iers_lookup_seconds
            .unwrap()(jd()),
        Some(43.4724)
    );
}

#[test]
fn ordinary_resolution_error_falls_back_and_is_retried() {
    thread_local! { static CALLS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) }; }
    fn missing(_: &IersAcquisitionRequest) -> Result<ResolvedIersAcquisition> {
        CALLS.with(|x| x.set(x.get() + 1));
        Err(Error::ConfigurationMissing {
            missing_config: "test_iers_directory",
            suggestion: "Supply directory.",
        })
    }
    let mut ctx = EphemerisContext::new();
    ctx.set_iers_boundary_resolver(Some(missing));
    for _ in 0..2 {
        assert!(
            deltat_with_config(jd(), &ctx.delta_t_config(DeltaTModel::EspenakMeeus, true)).is_ok()
        );
    }
    CALLS.with(|x| assert_eq!(x.get(), 2));
}

#[test]
fn unrelated_panic_is_not_converted_to_model_fallback() {
    let mut ctx = EphemerisContext::new();
    ctx.set_iers_boundary_resolver(Some(must_not_resolve));
    let failed = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        deltat_with_config(jd(), &ctx.delta_t_config(DeltaTModel::EspenakMeeus, true))
    }));
    assert!(failed.is_err());
}
