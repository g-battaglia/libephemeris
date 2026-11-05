//! Shared-image/cache lifetime and request configuration over explicit LEB data.
use libephemeris_core::{leb::LebFile, native as ephe};

#[test]
fn variants_share_one_image_and_cache_without_sharing_configuration() -> Result<(), ephe::Error> {
    fn send_sync<T: Send + Sync>() {}
    send_sync::<ephe::Session>();
    let path = std::env::var_os("KERYKEION_LEB_FILE")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| {
            std::path::PathBuf::from(std::env::var_os("HOME").unwrap())
                .join(".libephemeris/leb/base_core.leb2")
        });
    let file = LebFile::open(&path)?;
    let pointer = file.bytes()?.as_ptr();
    let mut builder = ephe::SessionBuilder::isolated();
    builder.leb_image(file.clone()).delta_t_override(0.0008);
    let first = builder.build()?;
    let mut variant = first.to_builder();
    variant
        .delta_t_override(0.0012)
        .sid_mode(ephe::SIDM_LAHIRI, 0.0, 0.0);
    let other = variant.build()?;
    assert_eq!(first.context().leb_files()[0].bytes()?.as_ptr(), pointer);
    assert_eq!(other.context().leb_files()[0].bytes()?.as_ptr(), pointer);
    assert_eq!(first.deltat(2451545.0)?, 0.0008);
    assert_eq!(other.deltat(2451545.0)?, 0.0012);
    let baseline = first
        .calc_ut(2451545.0, ephe::SUN, ephe::FLG_SWIEPH | ephe::FLG_SPEED)?
        .values
        .map(f64::to_bits);
    std::thread::scope(|scope| -> Result<(), ephe::Error> {
        let handles = (0..4)
            .map(|_| {
                scope.spawn(|| -> Result<(), ephe::Error> {
                    for _ in 0..8 {
                        assert_eq!(
                            baseline,
                            first
                                .calc_ut(2451545.0, ephe::SUN, ephe::FLG_SWIEPH | ephe::FLG_SPEED)?
                                .values
                                .map(f64::to_bits)
                        );
                        assert!(
                            other
                                .calc_ut(
                                    2451545.0,
                                    ephe::MOON,
                                    ephe::FLG_SWIEPH | ephe::FLG_SIDEREAL
                                )?
                                .values
                                .iter()
                                .all(|v| v.is_finite())
                        );
                    }
                    Ok(())
                })
            })
            .collect::<Vec<_>>();
        for handle in handles {
            handle.join().unwrap()?;
        }
        Ok(())
    })?;
    let (resident, budget) = file.chunk_cache_memory()?;
    assert!(resident > 0 && resident <= budget);
    assert!(file.chunk_cache_stats()?.0 > 0);
    file.cool();
    assert_eq!(other.context().leb_files()[0].chunk_cache_memory()?.0, 0);
    file.set_chunk_cache_byte_budget(0)?;
    assert_eq!(
        baseline,
        first
            .calc_ut(2451545.0, ephe::SUN, ephe::FLG_SWIEPH | ephe::FLG_SPEED)?
            .values
            .map(f64::to_bits)
    );
    assert_eq!(file.chunk_cache_memory()?.0, 0);
    Ok(())
}

#[test]
fn rebuilding_retains_sidereal_advisories_and_rejects_nonfinite_knobs() -> Result<(), ephe::Error> {
    let mut builder = ephe::SessionBuilder::isolated();
    builder.sid_mode(ephe::SIDM_LAHIRI | (1 << 20), 0.0, 0.0);
    let session = builder.build()?;
    assert!(!session.configuration_diagnostics().is_empty());
    assert_eq!(
        session.configuration_diagnostics(),
        session.to_builder().build()?.configuration_diagnostics()
    );
    assert!(
        ephe::SessionBuilder::isolated()
            .delta_t_override(f64::NAN)
            .build()
            .is_err()
    );
    assert!(
        ephe::SessionBuilder::isolated()
            .topo(0.0, 91.0, 0.0)
            .build()
            .is_err()
    );
    Ok(())
}
