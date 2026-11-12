//! Context-selected Delta-T policy through the high-level session calculations.
use libephemeris_core::{context::EphemerisContext, native as e, time::DeltaTModel};

#[test]
fn context_delta_t_model_survives_session_construction_and_rebuild() -> Result<(), e::Error> {
    for model in [DeltaTModel::Smh2016, DeltaTModel::EspenakMeeus] {
        let ctx = EphemerisContext::new().with_delta_t_model(Some(model));
        let direct = e::Session::from_context(ctx.clone());
        let built = e::SessionBuilder::from_context(ctx.clone()).build()?;
        let rebuilt = built.to_builder().build()?;
        for jd in [2378496.5, 2415020.5, 2451545.0] {
            let expected =
                libephemeris_core::time::deltat_with_config(jd, &ctx.delta_t_config(model, false))?;
            for session in [&direct, &built, &rebuilt] {
                assert_eq!(session.deltat(jd)?.to_bits(), expected.to_bits());
            }
        }
        let mut overridden = built.to_builder();
        overridden.delta_t_override(0.0012);
        assert_eq!(overridden.build()?.deltat(2451545.0)?, 0.0012);
        overridden.clear_delta_t_override();
        assert_eq!(
            overridden.build()?.deltat(2451545.0)?,
            built.deltat(2451545.0)?
        );
    }
    Ok(())
}

#[test]
fn ut_positions_use_the_selected_model_for_body_and_fixed_epoch_clocks() -> Result<(), e::Error> {
    let path = std::env::var_os("KERYKEION_LEB_FILE")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| {
            std::path::PathBuf::from(std::env::var_os("HOME").unwrap())
                .join(".libephemeris/leb/base_core.leb2")
        });
    let mut ctx = EphemerisContext::new()
        .with_mode(libephemeris_core::context::ProviderMode::SealedLeb)
        .with_delta_t_model(Some(DeltaTModel::EspenakMeeus));
    ctx.attach_leb_file(libephemeris_core::leb::LebFile::open(&path)?);
    let jd = 2451545.0;
    let expected_delta = libephemeris_core::time::deltat_with_config(
        jd,
        &ctx.delta_t_config(DeltaTModel::EspenakMeeus, false),
    )?;
    let mut builder = e::SessionBuilder::from_context(ctx);
    for sidereal in [false, true] {
        builder.sid_mode(e::SIDM_J2000, 0.0, 0.0);
        let session = builder.build()?;
        let flags = e::FLG_SWIEPH | e::FLG_SPEED | if sidereal { e::FLG_SIDEREAL } else { 0 };
        let ut = session.calc_ut(jd, e::MOON, flags)?;
        let tt = session.calc(jd + expected_delta, e::MOON, flags)?;
        assert_eq!(ut.values.map(f64::to_bits), tt.values.map(f64::to_bits));
    }
    Ok(())
}
