//! Native coefficient metadata records and inventory bindings.
use crate::{errors::to_pyerr, positions::runtime};
use libephemeris_core::{config::data::BodyCoverageRecord, leb};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use std::path::PathBuf;

/// Immutable copied metadata for a library body's stored coefficient coverage.
///
/// Constructor parameters correspond to read-only properties: body_id (library
/// integer), source/precision_class (strings), jd_start/jd_end (TT Julian days),
/// data_file/group (strings or None), reviewed (bool). All parameters are
/// required, including those accepting None. Direct construction merely stores
/// values: it does not validate finite/ordered bounds, labels, paths or digests.
///
/// Query functions populate records from retained images. An unlabeled multi-file
/// envelope can bridge gaps; contains() compares bounds, not actual serving
/// availability. reviewed denotes manifest identity for query-produced records,
/// not scientific acceptance; a manually constructed record can assert any value.
/// The record owns its scalars/strings and remains usable after reader closure.
///
/// Example:
///
/// ```python
/// from libephemeris.inventory import BodyCoverage
/// row = BodyCoverage(0, 'custom', 'unverified-local', 2451545.,
///                    2451546., None, None, False)
/// assert row.contains(2451545.)
/// metadata = row.to_dict()  # An independent mutable dict.
/// ```
#[pyclass(module = "libephemeris._core", frozen, skip_from_py_object)]
#[derive(Clone)]
struct BodyCoverage {
    inner: BodyCoverageRecord,
}
#[pymethods]
impl BodyCoverage {
    /// Copy the supplied metadata without numerical/file validation; see class docs.
    #[new]
    #[pyo3(signature=(body_id,source,precision_class,jd_start,jd_end,data_file,group,reviewed))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        body_id: i32,
        source: String,
        precision_class: String,
        jd_start: f64,
        jd_end: f64,
        data_file: Option<String>,
        group: Option<String>,
        reviewed: bool,
    ) -> Self {
        Self {
            inner: BodyCoverageRecord {
                body_id,
                source,
                precision_class,
                jd_start,
                jd_end,
                data_file,
                group,
                reviewed,
            },
        }
    }
    /// Return jd_start <= jd <= jd_end, with no time conversion or I/O.
    ///
    /// jd is in the stored bounds' TT clock. NaN compares False; other IEEE
    /// values follow ordinary comparisons against even manually supplied bounds.
    /// True does not certify gap-free coverage or evaluable coefficients.
    fn contains(&self, py: Python<'_>, jd: f64) -> bool {
        py.detach(|| self.inner.contains(jd))
    }
    /// Return a new dict with all eight property names and their scalar values.
    ///
    /// Editing it does not change this record; no file is reopened or verified.
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let result = PyDict::new(py);
        result.set_item("body_id", self.inner.body_id)?;
        result.set_item("source", &self.inner.source)?;
        result.set_item("precision_class", &self.inner.precision_class)?;
        result.set_item("jd_start", self.inner.jd_start)?;
        result.set_item("jd_end", self.inner.jd_end)?;
        result.set_item("data_file", &self.inner.data_file)?;
        result.set_item("group", &self.inner.group)?;
        result.set_item("reviewed", self.inner.reviewed)?;
        Ok(result)
    }
    /// Library body identifier, not a NAIF target code.
    #[getter]
    fn body_id(&self) -> i32 {
        self.inner.body_id
    }
    /// Calculation-source/provenance label; direct construction does not validate it.
    #[getter]
    fn source(&self) -> String {
        self.inner.source.clone()
    }
    /// Provenance class string, not an independently measured error bound.
    ///
    /// Query labels include ephemeris, analytical, numerical-model, mixed and
    /// unverified-local; a direct constructor can supply another string.
    #[getter]
    fn precision_class(&self) -> String {
        self.inner.precision_class.clone()
    }
    /// Inclusive TT JD start; multi-file diagnostic envelopes can span gaps.
    #[getter]
    fn jd_start(&self) -> f64 {
        self.inner.jd_start
    }
    /// Inclusive TT JD end; no constructor ordering/finiteness check is implied.
    #[getter]
    fn jd_end(&self) -> f64 {
        self.inner.jd_end
    }
    /// Selected image's path label, or None for an unlabeled multi-file envelope.
    #[getter]
    fn data_file(&self) -> Option<String> {
        self.inner.data_file.clone()
    }
    /// Coefficient partition inferred from the path, or None when unknown.
    #[getter]
    fn group(&self) -> Option<String> {
        self.inner.group.clone()
    }
    /// Manifest-identity flag, not an accuracy certificate or a fresh file check.
    ///
    /// Query-produced True means its immutable image matches a manifest digest;
    /// directly constructed rows simply retain the caller's boolean.
    #[getter]
    fn reviewed(&self) -> bool {
        self.inner.reviewed
    }
}
impl From<BodyCoverageRecord> for BodyCoverage {
    fn from(inner: BodyCoverageRecord) -> Self {
        Self { inner }
    }
}

/// Immutable expected artifact metadata, not a verified installation result.
///
/// name/kind/path/sha256 are strings and group is a string or None; every
/// constructor parameter is required. Values are copied without validating the
/// digest syntax, probing the path, opening data or normalizing labels. Records
/// returned by get_runtime_data_requirements() describe manifest-pinned LEB2
/// artifacts in cumulative tier/group order. All properties are read-only.
#[pyclass(module = "libephemeris._core", frozen, get_all, skip_from_py_object)]
#[derive(Clone)]
struct RuntimeDataRequirement {
    /// Artifact basename; generated records combine tier/group with a .leb2 suffix.
    name: String,
    /// Artifact kind; generated requirements use 'leb2'.
    kind: String,
    /// Distribution partition string, or None in a caller-constructed record.
    group: Option<String>,
    /// Expected path label, not proof that a file exists or is readable.
    path: String,
    /// Expected SHA-256 hex string; not a computed digest of a local file.
    sha256: String,
}
#[pymethods]
impl RuntimeDataRequirement {
    /// Copy supplied expected-artifact strings without file/digest validation.
    #[new]
    #[pyo3(signature=(name,kind,group,path,sha256))]
    fn new(
        name: String,
        kind: String,
        group: Option<String>,
        path: String,
        sha256: String,
    ) -> Self {
        Self {
            name,
            kind,
            group,
            path,
            sha256,
        }
    }
    /// Return a new dict containing name, kind, group, path, sha256 and scalar copies.
    ///
    /// Editing the dict does not change this record or verify an installed file.
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let result = PyDict::new(py);
        result.set_item("name", &self.name)?;
        result.set_item("kind", &self.kind)?;
        result.set_item("group", &self.group)?;
        result.set_item("path", &self.path)?;
        result.set_item("sha256", &self.sha256)?;
        Ok(result)
    }
}
impl From<leb::RuntimeDataRequirement> for RuntimeDataRequirement {
    fn from(value: leb::RuntimeDataRequirement) -> Self {
        Self {
            name: value.name().to_owned(),
            kind: value.kind().to_owned(),
            group: Some(value.group().to_owned()),
            path: value.path().to_owned(),
            sha256: value.sha256().to_owned(),
        }
    }
}
fn file_info<'py>(py: Python<'py>, value: leb::FileInfo) -> PyResult<Bound<'py, PyDict>> {
    let result = PyDict::new(py);
    let bodies = PyList::empty(py);
    result.set_item("name", value.name())?;
    result.set_item("path", value.path())?;
    result.set_item("group", value.group())?;
    result.set_item("size_bytes", value.size_bytes())?;
    result.set_item("reviewed", value.reviewed())?;
    result.set_item("body_count", value.body_count())?;
    for body in value.bodies() {
        let item = PyDict::new(py);
        item.set_item("body_id", body.body_id())?;
        item.set_item("jd_start", body.jd_start())?;
        item.set_item("jd_end", body.jd_end())?;
        bodies.append(item)?;
    }
    result.set_item("bodies", bodies)?;
    Ok(result)
}
/// Return stored LEB coverage metadata for library body_id, or None if absent.
///
/// Optional jd is TT JD, used to select the actual serving image when covered.
/// No jd (or an uncovered jd) can produce a diagnostic multi-image envelope
/// without a file label; this is not a continuous serving-domain claim. Retired
/// fictitious channels return None. Only stored coverage is reported, not every
/// analytical/SPK route available to calc(). Results own their scalar metadata.
/// Reader acquisition/configuration, missing required data, malformed indices
/// and closed-resource failures propagate through the public exception mapper.
/// jd is compared as supplied rather than converted or eagerly validated.
#[pyfunction(signature=(body_id,jd=None))]
fn get_body_coverage(
    py: Python<'_>,
    body_id: i32,
    jd: Option<f64>,
) -> PyResult<Option<BodyCoverage>> {
    py.detach(|| runtime().body_coverage(body_id, jd))
        .map(|value| value.map(Into::into))
        .map_err(|error| to_pyerr(py, error))
}
/// Alias calculation for get_body_coverage(body_id, jd=None).
///
/// Accepts a library ID and optional TT JD; returns copied BodyCoverage or None,
/// with identical stored-reader selection, envelope limitations and errors.
#[pyfunction(signature=(body_id,jd=None))]
fn coverage(py: Python<'_>, body_id: i32, jd: Option<f64>) -> PyResult<Option<BodyCoverage>> {
    py.detach(|| runtime().body_coverage(body_id, jd))
        .map(|value| value.map(Into::into))
        .map_err(|error| to_pyerr(py, error))
}
/// Inspect body coverage in the supplied LebReader, without default-runtime selection.
///
/// body_id is a library ID; optional jd is a TT JD compared without conversion
/// or eager finiteness validation. Returns an owned BodyCoverage or None for an
/// absent/retired channel. Does not close/mutate the reader or select another
/// file. Closed/backing-resource and malformed-index errors propagate as native
/// ConfigurationError/DataNotFoundError subclasses. Coverage is metadata, not a
/// numerical evaluation or scientific-accuracy certificate.
#[pyfunction(signature=(reader,body_id,jd=None))]
fn get_reader_body_coverage(
    py: Python<'_>,
    reader: &crate::state::LebReader,
    body_id: i32,
    jd: Option<f64>,
) -> PyResult<Option<BodyCoverage>> {
    let image = reader.image();
    py.detach(|| libephemeris_core::config::data::stored_body_coverage(&[image], body_id, jd))
        .map(|value| value.map(Into::into))
        .map_err(|error| to_pyerr(py, error))
}
/// Return a new dict describing the default runtime's selected stored images.
///
/// Keys: mode, precision_tier, network_policy_configured/effective (strings),
/// ready (structural bool), reader_type (label or None), body_count (int), files
/// (list of new dicts), plus error when reported. Each file dict contains name,
/// path, group, size_bytes, reviewed, body_count, bodies; body rows have body_id,
/// jd_start/jd_end (TT JD). Size/group can be None. reviewed checks the retained
/// image's manifest identity, not scientific precision or current path contents.
/// May acquire configured readers and hash their immutable images; it does not
/// download or evaluate positions. Forced non-LEB modes report no stored files.
/// Configuration/reader admission failures propagate; ready is not proof that
/// every public calculation is supported. Returned metadata does not retain a
/// mutable view of the runtime.
#[pyfunction]
fn get_leb_inventory<'py>(py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
    let value = py
        .detach(|| runtime().leb_inventory())
        .map_err(|error| to_pyerr(py, error))?;
    let result = PyDict::new(py);
    let files = PyList::empty(py);
    result.set_item("mode", value.mode)?;
    result.set_item("precision_tier", value.precision_tier)?;
    result.set_item("network_policy_configured", value.network_policy_configured)?;
    result.set_item("network_policy_effective", value.network_policy_effective)?;
    result.set_item("ready", value.ready)?;
    result.set_item("reader_type", value.reader_type)?;
    result.set_item("body_count", value.body_count)?;
    for file in value.files {
        files.append(file_info(py, file)?)?;
    }
    result.set_item("files", files)?;
    if let Some(error) = value.error {
        result.set_item("error", error)?;
    }
    Ok(result)
}
/// Return a tuple of immutable RuntimeDataRequirement records for a tier.
///
/// tier is exactly 'base', 'medium' or 'extended'; None uses the resolved active
/// tier. Requirements include all lower eligible tiers and their groups, with
/// expected paths under the runtime data root and manifest SHA-256 pins. No
/// files are opened/downloaded and no installed digests are verified. Raises
/// ValueError for unknown tiers, DataNotFoundError for missing manifest pins,
/// and native configuration/path exceptions when resolving the data root.
///
/// Example:
///
/// ```python
/// from libephemeris.inventory import get_runtime_data_requirements
/// requirements = get_runtime_data_requirements('base')
/// expected_paths = [record.path for record in requirements]
/// ```
#[pyfunction(signature=(tier=None))]
fn get_runtime_data_requirements<'py>(
    py: Python<'py>,
    tier: Option<&str>,
) -> PyResult<Bound<'py, PyTuple>> {
    let values = py
        .detach(|| runtime().runtime_data_requirements(tier))
        .map_err(|error| to_pyerr(py, error))?;
    let records = values
        .into_iter()
        .map(|value| Py::new(py, RuntimeDataRequirement::from(value)))
        .collect::<PyResult<Vec<_>>>()?;
    PyTuple::new(py, records)
}
/// Read a regular LEB file's header/directory metadata into a new dict.
///
/// path accepts a filesystem path. Reads bounded index/body metadata, not
/// coefficient payloads; keep the file stable during inspection. Returns name,
/// path, group, size_bytes, reviewed=False, body_count and bodies (list of dicts
/// with library body_id and inclusive TT jd_start/jd_end). Does not acquire a
/// runtime reader, change configuration, hash the file or establish precision.
/// Missing files and corrupt/unsupported containers raise DataNotFoundError
/// subclasses through the native mapper, not necessarily built-in OSError.
#[pyfunction]
fn inspect_leb_file<'py>(py: Python<'py>, path: PathBuf) -> PyResult<Bound<'py, PyDict>> {
    let value = py
        .detach(|| leb::inspect_file(&path))
        .map_err(|error| to_pyerr(py, error))?;
    file_info(py, value)
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<BodyCoverage>()?;
    module.add_class::<RuntimeDataRequirement>()?;
    module.add_function(wrap_pyfunction!(get_body_coverage, module)?)?;
    module.add_function(wrap_pyfunction!(get_reader_body_coverage, module)?)?;
    module.add_function(wrap_pyfunction!(coverage, module)?)?;
    module.add_function(wrap_pyfunction!(get_leb_inventory, module)?)?;
    module.add_function(wrap_pyfunction!(get_runtime_data_requirements, module)?)?;
    module.add_function(wrap_pyfunction!(inspect_leb_file, module)?)?;
    Ok(())
}
