//! Native SPK discovery and tier acquisition summaries.
use crate::{errors::to_pyerr, positions::runtime};
use libephemeris_core::config::data::AcquisitionSummary;
use libephemeris_core::spk::EnsureVerdict;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::BTreeMap;
use std::path::PathBuf;

pub(crate) fn summary<'py>(
    py: Python<'py>,
    value: AcquisitionSummary,
) -> PyResult<Bound<'py, PyDict>> {
    let ephemeris_path = value.ephemeris_path;
    let value = value.summary;
    let result = PyDict::new(py);
    let ephemeris = PyDict::new(py);
    let spks = PyDict::new(py);
    let counts = PyDict::new(py);
    result.set_item("tier", &value.tier)?;
    ephemeris.set_item("file", &value.ephemeris_file)?;
    ephemeris.set_item("status", &value.ephemeris_status)?;
    if let Some(path) = ephemeris_path {
        ephemeris.set_item("path", path.to_string_lossy().as_ref())?;
    }
    for entry in &value.entries {
        spks.set_item(
            &entry.body_name,
            match entry.verdict {
                EnsureVerdict::Cached => "cached",
                EnsureVerdict::Downloaded => "downloaded",
                EnsureVerdict::Fallback => &entry.detail,
            },
        )?;
    }
    counts.set_item("total", value.entries.len() + 1)?;
    counts.set_item("cached", value.cached)?;
    counts.set_item("downloaded", value.downloaded)?;
    counts.set_item("fallback", value.fallback)?;
    counts.set_item("errors", value.errors)?;
    result.set_item("ephemeris", ephemeris)?;
    result.set_item("spks", spks)?;
    result.set_item("summary", counts)?;
    Ok(result)
}
/// Scan path for known minor-body .bsp targets and register local images.
///
/// path accepts str or os.PathLike and is a directory, not a kernel file.
/// Return a copied {body_name: status} dictionary, with statuses "registered",
/// "already_registered" or "error: ...". Wider file spans may replace an
/// existing mapping; this does not prove target-chain coverage at every epoch.
/// Non-directories return {}; unreadable/invalid candidates are skipped,
/// while directory-listing errors propagate. No download is attempted.
#[pyfunction]
fn discover_local_spks(py: Python<'_>, path: PathBuf) -> PyResult<BTreeMap<String, String>> {
    py.detach(|| runtime().discover_local_spks(&path))
        .map_err(|error| to_pyerr(py, error))
}
/// Provision the active tier's main JPL kernel and known minor-body SPKs.
///
/// Existing registered bodies are retained unless force_download=True;
/// forced provisioning requests replacement files. Downloads, validation and
/// registration happen eagerly under the current network policy. Completed
/// installations remain if a later operation fails. show_progress is accepted
/// for compatibility but does not enable a Python progress callback.
///
/// Return a fresh nested dictionary: tier (str); ephemeris with file/status
/// and path (str, only for a cached main file); spks mapping body names to
/// "cached", "downloaded" or a fallback reason; summary with integer total,
/// cached, downloaded, fallback and errors counts (including the main file).
/// Main acquisition errors can be reported in status, and ordinary minor-body
/// failures as fallbacks; a required sealed-network refusal raises
/// NetworkSealedError rather than returning a completed summary.
#[pyfunction(signature=(force_download=false,show_progress=true))]
fn ensure_all_ephemerides<'py>(
    py: Python<'py>,
    force_download: bool,
    show_progress: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let _ = show_progress;
    let result = py
        .detach(|| runtime().ensure_all_ephemerides(force_download))
        .map_err(|error| to_pyerr(py, error))?;
    summary(py, result)
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(discover_local_spks, module)?)?;
    module.add_function(wrap_pyfunction!(ensure_all_ephemerides, module)?)?;
    Ok(())
}
