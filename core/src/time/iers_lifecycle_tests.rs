//! Synthetic local-file/row lifecycle controls; no astronomical provider.
use super::*;
use std::sync::atomic::{AtomicU64, Ordering};

struct Local(PathBuf);
impl Local {
    fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "lephem-iers-lifecycle-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
    fn source(&self) -> IersContextSource {
        IersContextSource::default().with_request(IersAcquisitionRequest {
            cache_dir: Some(self.0.clone()),
            auto_download: Some(false),
            ..Default::default()
        })
    }
    fn monthly(&self, text: &str) {
        std::fs::write(self.0.join(DELTA_T_CACHE_FILENAME), text).unwrap();
    }
}
impl Drop for Local {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

#[test]
fn explicit_reload_reparses_without_mutating_retained_rows() {
    let local = Local::new();
    local.monthly("2000 1 1 20.0\n");
    let first = local.source().loaded(false).unwrap();
    let retained = first.clone();
    local.monthly("2000 1 1 30.0\n2000 2 1 31.0\n");
    let next = first.loaded(false).unwrap();
    assert_eq!(next.loaded_counts().delta_t, 2);
    assert_eq!(retained.loaded_counts().delta_t, 1);
    assert!(!Arc::ptr_eq(&first.cache, &next.cache));
}

#[test]
fn finals_file_view_survives_missing_reload_and_retained_clear() {
    let local = Local::new();
    let path = local.0.join(FINALS_CACHE_FILENAME);
    let line = |value| format!("230225 {:>8}{}I{:>10}\n", "60000.0", " ".repeat(42), value);
    std::fs::write(&path, line("nan")).unwrap();
    let first = local.source().loaded(false).unwrap();
    let retained = first.row_snapshot().unwrap();
    assert_eq!(first.loaded_counts().finals, 0);
    assert_eq!(first.loaded_source_counts().finals, 1);
    assert!(retained.source_finals_rows()[0].1.is_nan());
    std::fs::remove_file(&path).unwrap();
    let missing = first.loaded(false).unwrap();
    assert!(
        missing.row_snapshot().unwrap().source_finals_rows()[0]
            .1
            .is_nan()
    );
    std::fs::write(&path, line("inf")).unwrap();
    let next = missing.loaded(false).unwrap();
    assert_eq!(
        next.row_snapshot().unwrap().source_finals_rows()[0].1,
        f64::INFINITY
    );
    assert!(retained.source_finals_rows()[0].1.is_nan());
    assert_eq!(
        next.file_info_with_clock(|| Ok(0.0)).unwrap().data_points,
        0
    );
    assert_eq!(
        next.source_file_info_with_clock(|| Ok(0.0))
            .unwrap()
            .data_points,
        1
    );
    assert_eq!(next.loaded_ranges().0, None);
    assert_eq!(next.loaded_source_ranges().0, Some((2460000.5, 2460000.5)));
    assert_eq!(next.cleared().loaded_source_counts().finals, 0);
    assert_eq!(next.loaded_source_counts().finals, 1);
}

#[test]
fn leap_file_view_keeps_offsets_and_c04_sanity_across_reload_and_clear() {
    let local = Local::new();
    let path = local.0.join(LEAP_SECONDS_CACHE_FILENAME);
    std::fs::write(&path, "41317.0 1 1 1972 10\n60000.0 25 2 2023 nan\n").unwrap();
    local.monthly("2023 2 25 60000.0 0 0 0.1\n");
    let first = local.source().loaded(false).unwrap();
    let retained = first.row_snapshot().unwrap();
    assert!(retained.prepare_source_views(&[0, 1]));
    assert!(retained.source_tai_utc(60000.0).is_nan());
    assert_eq!(retained.leap_seconds().len(), 1);
    assert_eq!(retained.observed_rows().len(), 1);
    assert!(retained.source_observed_rows().is_empty());
    let info = first.source_file_info_with_clock(|| Ok(0.0)).unwrap();
    assert_eq!((info.leap_second_entries, info.delta_t_entries), (2, 0));
    std::fs::write(&path, "60000.0 25 2 2023 45\n").unwrap();
    let next = first.loaded(false).unwrap();
    let current = next.row_snapshot().unwrap();
    assert!(current.prepare_source_views(&[0]));
    assert_eq!(current.source_tai_utc(60000.0), 45.0);
    assert!(retained.source_tai_utc(60000.0).is_nan());
    assert_eq!(next.cleared().loaded_source_counts().leap_seconds, 0);
    assert_eq!(first.loaded_source_counts().leap_seconds, 2);
}

#[test]
fn offline_force_is_reload_not_implicit_network_permission() {
    let local = Local::new();
    local.monthly("2000 1 1 20.0\n");
    let next = local.source().loaded(true).unwrap();
    assert_eq!(next.loaded_counts().delta_t, 1);
    assert!(!local.0.join(FINALS_CACHE_FILENAME).exists());
}

#[test]
fn missing_finals_retained_but_missing_monthly_is_reset_on_reload() {
    let local = Local::new();
    let rows = IersDataSnapshot {
        source_leap_fields: leap_fallback_view(),
        source_observed_text: Some("2000 1 1 20.0\n".into()),
        source_views: OnceLock::new(),
        finals: vec![IersDataPoint::new(51_544.0, 2000, 1, 1, 0.0, None, false).unwrap()],
        source_finals: vec![(51_544.0, 0.0)],
        leaps: hardcoded_leap_seconds(),
        observed: parse_delta_t("2000 1 1 20.0\n", &hardcoded_leap_seconds()),
    };
    let next = local.source().with_rows(rows).loaded(false).unwrap();
    assert_eq!(next.loaded_counts().finals, 1);
    assert_eq!(next.loaded_counts().delta_t, 0);
    assert!(next.loaded_counts().leap_seconds > 0);
}

#[test]
fn cold_cache_info_counts_precede_lazy_range_load() {
    let local = Local::new();
    local.monthly("2000 1 1 20.0\n2000 2 1 21.0\n");
    let source = local.source();
    let cold = source.cache_info(0.0).unwrap();
    assert_eq!(cold.delta_t_entries, 0);
    assert_eq!(cold.leap_second_entries, 0);
    assert!(cold.delta_t_range.is_some());
    let warm = source.cache_info(0.0).unwrap();
    assert_eq!(warm.delta_t_entries, 2);
    assert!(warm.leap_second_entries > 0);
}

#[test]
fn memory_clear_does_not_delete_files_or_clear_retained_rows() {
    let local = Local::new();
    local.monthly("2000 1 1 20.0\n");
    let first = local.source().loaded(false).unwrap();
    let retained = first.clone();
    let cleared = first.cleared();
    assert_eq!(cleared.loaded_counts().delta_t, 0);
    assert_eq!(retained.loaded_counts().delta_t, 1);
    assert!(local.0.join(DELTA_T_CACHE_FILENAME).exists());
}

#[test]
fn disk_delete_clears_new_owner_and_preserves_old_snapshot() {
    let local = Local::new();
    local.monthly("2000 1 1 20.0\n");
    std::fs::write(local.0.join("unrelated.txt"), "kept").unwrap();
    let old = local.source().loaded(false).unwrap();
    let (next, deleted) = old.deleted_files().unwrap();
    assert_eq!(deleted, 1);
    assert_eq!(next.loaded_counts().delta_t, 0);
    assert_eq!(old.loaded_counts().delta_t, 1);
    assert!(local.0.join("unrelated.txt").exists());
}

#[test]
fn failed_reload_leaves_context_generation_and_policy_unchanged() {
    let local = Local::new();
    let blocked = local.0.join("file");
    std::fs::write(&blocked, "not a directory").unwrap();
    let mut context = crate::context::EphemerisContext::new();
    context.attach_iers_rows(IersDataSnapshot::from_texts(
        None,
        None,
        Some("2000 1 1 20.0\n"),
    ));
    context.set_iers_cache_dir(Some(blocked.join("child")));
    context.set_user_delta_t(Some(-0.0));
    assert!(context.load_iers_data(false).is_err());
    assert_eq!(context.iers_source().loaded_counts().delta_t, 1);
    assert_eq!(
        context.user_delta_t_days().unwrap().to_bits(),
        (-0.0_f64).to_bits()
    );
    assert!(!context.effective_iers_delta_t());
}

#[test]
fn immutable_session_clear_preserves_original_and_reset_stays_warm() {
    let mut context = crate::context::EphemerisContext::new();
    context.attach_iers_rows(IersDataSnapshot::from_texts(
        None,
        None,
        Some("2000 1 1 20.0\n"),
    ));
    context.reset_session();
    assert_eq!(context.iers_source().loaded_counts().delta_t, 1);
    let old = crate::engine_api::Session::from_context(context);
    let next = old.with_cleared_iers_cache();
    assert_eq!(old.iers_loaded_counts().delta_t, 1);
    assert_eq!(next.iers_loaded_counts().delta_t, 0);
}

#[test]
fn selected_policy_valueerror_falls_through_without_a_socket() {
    let local = Local::new();
    let policy = Err(Error::InvalidNetworkPolicy {
        value: "invalid".into(),
        source: "LIBEPHEMERIS_NETWORK_POLICY",
    });
    assert!(
        !IersContextSource::download_selected_in(
            &local.0,
            super::super::IersFileKind::Finals,
            policy
        )
        .unwrap()
    );
    let unknown = Err(Error::ConfigurationMissing {
        missing_config: "synthetic_unknown",
        suggestion: "not a caught download error",
    });
    assert!(
        IersContextSource::download_selected_in(
            &local.0,
            super::super::IersFileKind::Finals,
            unknown
        )
        .is_err()
    );
}

#[test]
fn independently_captured_paths_retained_after_settings_change() {
    let first = Local::new();
    let second = Local::new();
    first.monthly("2000 1 1 20.0\n");
    second.monthly("2000 1 1 30.0\n2000 2 1 31.0\n");
    let paths = super::super::iers_cache_filenames().map(|name| first.0.join(name));
    let now = second.source();
    let loaded = now.reloaded_paths(&paths);
    assert_eq!(loaded.loaded_counts().delta_t, 1);
    assert_eq!(loaded.request.cache_dir, Some(second.0.clone()));
}

#[test]
fn decimal_one_place_rounding_uses_exact_binary_input() {
    let round = super::super::round_age_days_tenth;
    assert_eq!(round(0.05).to_bits(), 0.1_f64.to_bits());
    assert_eq!(round(1.15).to_bits(), 1.1_f64.to_bits());
    assert_eq!(round(-0.05).to_bits(), (-0.1_f64).to_bits());
    assert_eq!(round(-0.0).to_bits(), (-0.0_f64).to_bits());
    assert_eq!(round(2.25), 2.2);
}

#[test]
fn valid_negative_unix_time_and_ordered_clock_samples() {
    use std::time::{Duration, UNIX_EPOCH};
    assert_eq!(
        super::super::iers_unix_seconds(UNIX_EPOCH - Duration::from_millis(1500)),
        -1.5
    );
    let local = Local::new();
    for name in super::super::iers_cache_filenames() {
        std::fs::write(local.0.join(name), "").unwrap();
    }
    let mut clocks = 0;
    let info = local
        .source()
        .cache_info_with_clock(|| {
            clocks += 1;
            Ok(clocks as f64)
        })
        .unwrap();
    assert_eq!(clocks, 3);
    assert!(info.finals.age_days.is_some());
    assert!(info.leap_seconds.age_days.is_some());
    assert!(info.delta_t.age_days.is_some());
}
