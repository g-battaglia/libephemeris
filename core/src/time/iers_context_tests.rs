//! Local-row acquisition, retry and owned-snapshot controls; no network/provider.
use super::*;
use std::sync::atomic::{AtomicU64, Ordering};

// Input rows copied from the admitted monthly data file, not native outcomes.
const ROWS: &str = "1973 2 1 43.4724\n1973 3 1 43.5648\n";
fn first_jd() -> f64 {
    crate::time::mjd_to_jd(crate::time::calendar_to_mjd(1973, 2, 1))
        .expect("finite monthly fixture epoch")
}
struct Temp(PathBuf);
impl Temp {
    fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "lephem-iers-owner-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
    fn source(&self) -> IersContextSource {
        IersContextSource::default().with_request(IersAcquisitionRequest {
            cache_dir: Some(self.0.clone()),
            ..Default::default()
        })
    }
}
impl Drop for Temp {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

#[test]
fn observed_snapshot_covers_inclusive_endpoints_without_extrapolation() {
    let rows = IersDataSnapshot::from_texts(None, None, Some(ROWS));
    assert_eq!(rows.counts().delta_t, 2);
    assert!(rows.counts().leap_seconds > 0);
    assert_eq!(rows.delta_t_seconds(first_jd()).unwrap(), Some(43.4724));
    let last = crate::time::mjd_to_jd(crate::time::calendar_to_mjd(1973, 3, 1))
        .expect("finite monthly fixture epoch");
    assert_eq!(rows.delta_t_seconds(last).unwrap(), Some(43.5648));
    assert_eq!(rows.delta_t_seconds(first_jd() - 1.0).unwrap(), None);
    assert_eq!(rows.delta_t_seconds(last + 1.0).unwrap(), None);
    assert!(matches!(
        rows.delta_t_seconds(f64::NAN),
        Err(Error::NonFiniteArgument { .. })
    ));
}

#[test]
fn observed_rows_win_over_covering_finals() {
    let rows = IersDataSnapshot {
        source_leap_fields: leap_fallback_view(),
        source_observed_text: Some(ROWS.into()),
        source_views: OnceLock::new(),
        source_finals: vec![(crate::time::calendar_to_mjd(1973, 2, 1), 0.5)],
        finals: vec![
            IersDataPoint::new(
                crate::time::calendar_to_mjd(1973, 2, 1),
                1973,
                2,
                1,
                0.5,
                None,
                false,
            )
            .unwrap(),
        ],
        leaps: hardcoded_leap_seconds(),
        observed: parse_delta_t(ROWS, &hardcoded_leap_seconds()),
    };
    assert_eq!(rows.delta_t_seconds(first_jd()).unwrap(), Some(43.4724));
}

#[test]
fn finals_file_view_retains_nonfinite_values_and_original_keys() {
    let line = |mjd, ut1, error| {
        format!(
            "230225 {:>8}{}I{:>10}{:>10}\n",
            mjd,
            " ".repeat(42),
            ut1,
            error
        )
    };
    let text = [
        line("60002.0", "0.5", "nan"),
        line("60000.0", "nan", ""),
        line("60001.0", "inf", "bad"),
        line("60002.0", "0.75", "bad"),
        line("nan", "1.0", ""),
        line("nan", "2.0", ""),
        line("-0.0", "0.2", ""),
        line("0.0", "0.3", ""),
    ]
    .concat();
    let snapshot = IersDataSnapshot::from_texts(Some(&text), None, None);
    let view = snapshot.source_finals_rows();
    assert_eq!(view.len(), 6);
    assert_eq!(view[0], (60002.0, 0.75));
    assert_eq!(view[1].0, 60000.0);
    assert!(view[1].1.is_nan());
    assert_eq!(view[2], (60001.0, f64::INFINITY));
    assert!(view[3].0.is_nan() && view[4].0.is_nan());
    assert_eq!((view[3].1, view[4].1), (1.0, 2.0));
    assert_eq!(view[5].0.to_bits(), (-0.0f64).to_bits());
    assert_eq!(view[5].1, 0.3);
    assert!(snapshot.finals_rows().iter().all(|row| {
        row.mjd().is_finite()
            && row.ut1_utc().is_finite()
            && row.ut1_utc_err().is_none_or(f64::is_finite)
    }));
    let first_nan = line("nan", "1.0", "") + &line("60000.0", "2.0", "");
    let view = IersDataSnapshot::from_texts(Some(&first_nan), None, None);
    let (lo, hi) = view.source_finals_range().unwrap();
    assert!(lo.is_nan() && hi.is_nan());
}

#[test]
fn latin1_comments_do_not_discard_monthly_rows() {
    let temp = Temp::new();
    let mut data = vec![b'#', 0xff, b'\n'];
    data.extend_from_slice(ROWS.as_bytes());
    std::fs::write(temp.0.join(DELTA_T_CACHE_FILENAME), data).unwrap();
    assert_eq!(temp.source().provider()(first_jd()), Some(43.4724));
}

#[test]
fn ordinary_missing_tables_retain_leap_fallback_until_explicit_clear() {
    let temp = Temp::new();
    let owner = temp.source();
    assert_eq!(owner.provider()(first_jd()), None);
    std::fs::write(temp.0.join(DELTA_T_CACHE_FILENAME), ROWS).unwrap();
    assert_eq!(owner.provider()(first_jd()), None);
    assert_eq!(owner.cleared().provider()(first_jd()), Some(43.4724));
}

#[test]
fn failed_directory_creation_is_retried_without_failed_snapshot() {
    let temp = Temp::new();
    let parent = temp.0.join("parent");
    std::fs::write(&parent, b"file blocks requested child").unwrap();
    let directory = parent.join("iers_cache");
    let owner = IersContextSource::default().with_request(IersAcquisitionRequest {
        cache_dir: Some(directory.clone()),
        ..Default::default()
    });
    assert_eq!(owner.provider()(first_jd()), None);
    std::fs::remove_file(&parent).unwrap();
    std::fs::create_dir_all(&directory).unwrap();
    std::fs::write(directory.join(DELTA_T_CACHE_FILENAME), ROWS).unwrap();
    assert_eq!(owner.provider()(first_jd()), Some(43.4724));
}

#[test]
fn cache_path_store_preserves_loaded_rows_then_clear_selects_new_path() {
    let a = Temp::new();
    let b = Temp::new();
    std::fs::write(a.0.join(DELTA_T_CACHE_FILENAME), ROWS).unwrap();
    std::fs::write(b.0.join(DELTA_T_CACHE_FILENAME), "1973 2 1 70.0\n").unwrap();
    let owner = a.source();
    assert_eq!(owner.provider()(first_jd()), Some(43.4724));
    let changed = owner.with_request(IersAcquisitionRequest {
        cache_dir: Some(b.0.clone()),
        ..Default::default()
    });
    assert_eq!(changed.provider()(first_jd()), Some(43.4724));
    assert_eq!(changed.cleared().provider()(first_jd()), Some(70.0));
    assert_eq!(owner.provider()(first_jd()), Some(43.4724));
}

#[test]
fn sealed_auto_missing_files_still_loads_existing_observed_rows() {
    let temp = Temp::new();
    std::fs::write(temp.0.join(DELTA_T_CACHE_FILENAME), ROWS).unwrap();
    let owner = IersContextSource::default().with_request(IersAcquisitionRequest {
        cache_dir: Some(temp.0.clone()),
        auto_download: Some(true),
        network_policy: Some(NetworkPolicy::Sealed),
        ..Default::default()
    });
    assert_eq!(owner.provider()(first_jd()), Some(43.4724));
    assert!(!temp.0.join(FINALS_CACHE_FILENAME).exists());
}

#[test]
fn successful_snapshot_owns_rows_after_files_are_removed() {
    let temp = Temp::new();
    std::fs::write(temp.0.join(DELTA_T_CACHE_FILENAME), ROWS).unwrap();
    let owner = temp.source();
    assert_eq!(owner.provider()(first_jd()), Some(43.4724));
    let clone = owner.clone();
    std::fs::remove_file(temp.0.join(DELTA_T_CACHE_FILENAME)).unwrap();
    assert_eq!(clone.provider()(first_jd()), Some(43.4724));
    assert_eq!(owner.cleared().provider()(first_jd()), None);
}

#[test]
fn sdk_without_explicit_rows_or_directory_does_not_resolve_process_home() {
    assert_eq!(IersContextSource::default().provider()(first_jd()), None);
    assert_eq!(IersContextSource::default().provider()(first_jd()), None);
}

#[test]
fn equality_distinguishes_owned_attached_row_values() {
    let a = IersContextSource::default().with_rows(IersDataSnapshot::from_texts(
        None,
        None,
        Some(ROWS),
    ));
    let b = IersContextSource::default().with_rows(IersDataSnapshot::from_texts(
        None,
        None,
        Some("1973 2 1 70.0\n"),
    ));
    assert_eq!(a, a.clone());
    assert_ne!(a, b);
}

#[test]
fn ordinary_policy_refusal_still_retains_local_rows_and_leap_fallback() {
    let temp = Temp::new();
    std::fs::write(temp.0.join(DELTA_T_CACHE_FILENAME), ROWS).unwrap();
    let rows = acquire(ResolvedIersAcquisition {
        directory: temp.0.clone(),
        auto_download: true,
        network: Err(Error::InvalidNetworkPolicy {
            value: "invalid".into(),
            source: "test",
        }),
    })
    .unwrap();
    assert_eq!(rows.delta_t_seconds(first_jd()).unwrap(), Some(43.4724));
    assert!(rows.counts().leap_seconds > 0);
    assert!(!temp.0.join(FINALS_CACHE_FILENAME).exists());
}
