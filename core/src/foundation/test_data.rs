//! Read-only owner-selected ephemeris assets for native tests.
use std::{path::PathBuf, sync::OnceLock};

/// Explicit test file, otherwise the current production core in the user's data root.
pub(crate) fn core_leb_path() -> PathBuf {
    std::env::var_os("LIBEPHEMERIS_TEST_LEB2")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            let root = std::env::var_os("LIBEPHEMERIS_TEST_DATA_DIR")
                .map(PathBuf::from)
                .unwrap_or_else(|| {
                    PathBuf::from(
                        std::env::var_os("HOME").expect("HOME or explicit test data root required"),
                    )
                    .join(".libephemeris")
                });
            root.join("leb/base_core.leb2")
        })
}

/// The exact published distribution used by manifest and pinned-value regressions.
///
/// This path is deliberately separate from `core_leb_path`: an owner-selected
/// file may contain different bytes and must never acquire the distribution's
/// trust pin merely by sharing its filename. Current-data qualification uses
/// the owner-selected path above and authenticates that input independently.
pub(crate) fn published_core_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data/leb2/base_core.leb2")
}

/// Exact historical data-v3 image used by unchanged, data-v3 numerical references.
///
/// Provision this input explicitly; never substitute the newer bundled image
/// merely because both files are named `base_core.leb2`.
#[allow(dead_code, reason = "Used by the historical integration regression.")]
pub(crate) fn historical_core_path() -> PathBuf {
    std::env::var_os("LIBEPHEMERIS_TEST_DATA_V3_LEB2")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            core_leb_path()
                .parent()
                .expect("core LEB path must have a parent")
                .join("historical/data-v3/base_core.leb2")
        })
}

/// One process-owned immutable image; static backing tests retain real static borrows.
pub(crate) fn core_bytes() -> &'static [u8] {
    static BYTES: OnceLock<Vec<u8>> = OnceLock::new();
    BYTES
        .get_or_init(|| std::fs::read(core_leb_path()).expect("owner-selected core LEB must exist"))
        .as_slice()
}

/// Admit the separately identified historical and owner-selected production distributions.
/// Numeric expectations and tolerances remain in each test unchanged.
pub(crate) fn admitted_core_hash(hash: &str) -> bool {
    matches!(
        hash,
        "5d708bdbe3e799e0802ba575984e57a3c5e44720dbfa1b4a01cf826640e0cb82"
            | "b891de02f82b68346046da67389fd3713aa81ab299287c632cb20cc396397964"
    )
}
