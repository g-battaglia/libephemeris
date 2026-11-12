// SPDX-License-Identifier: AGPL-3.0
//! Native exception classes, typed public constructors and core error conversion.
//!
//! Error conversion uses the core's primary failure and documented scalar
//! attributes. It does not copy exception chains or interpreter object state.

use libephemeris_core::errors::{Error as CoreError, ErrorCategory};
use pyo3::exceptions::{
    PyException, PyFileNotFoundError, PyImportError, PyOverflowError, PyRuntimeError, PyTypeError,
    PyValueError, PyZeroDivisionError,
};
use pyo3::prelude::*;
use pyo3::types::{PyCFunction, PyDict, PyTuple, PyType};

/// Validate finite latitude in degrees north on [-90, 90]; return None.
///
/// Raises CoordinateError with coordinate_name, value, min_value/max_value.
/// func_name is an accepted compatibility string and does not affect validation
/// or error formatting. No observer or runtime settings are changed.
#[pyfunction(signature = (lat, func_name=""))]
fn validate_latitude(py: Python<'_>, lat: f64, func_name: &str) -> PyResult<()> {
    let _ = func_name;
    py.detach(|| libephemeris_core::context::Observer::new(0.0, lat, 0.0))
        .map(|_| ())
        .map_err(|error| to_pyerr(py, error))
}

/// Validate finite longitude in degrees east on [-180, 360]; return None.
///
/// Both signed-west and east-positive conventions are accepted without wrapping.
/// Raises CoordinateError with value/bounds for invalid input. func_name is an
/// accepted but unused compatibility string; no runtime settings are changed.
#[pyfunction(signature = (lon, func_name=""))]
fn validate_longitude(py: Python<'_>, lon: f64, func_name: &str) -> PyResult<()> {
    let _ = func_name;
    py.detach(|| libephemeris_core::context::Observer::new(lon, 0.0, 0.0))
        .map(|_| ())
        .map_err(|error| to_pyerr(py, error))
}

/// Validate latitude then longitude, in degrees north/east; return None.
///
/// Requires finite lat in [-90, 90] and lon in [-180, 360]. The first invalid
/// component raises CoordinateError with its name/value/bounds. func_name is
/// accepted but unused. This does not select an observer or normalize coordinates.
///
/// Example:
///
/// ```python
/// from libephemeris.exceptions import CoordinateError, validate_coordinates
/// try:
///     validate_coordinates(91.0, 12.0)
/// except CoordinateError as error:
///     assert error.coordinate_name == 'latitude'
/// ```
#[pyfunction(signature = (lat, lon, func_name=""))]
fn validate_coordinates(py: Python<'_>, lat: f64, lon: f64, func_name: &str) -> PyResult<()> {
    let _ = func_name;
    py.detach(|| libephemeris_core::context::Observer::new(lon, lat, 0.0))
        .map(|_| ())
        .map_err(|error| to_pyerr(py, error))
}

/// Complete exception initialization without changing its constructed attributes.
///
/// Arguments are accepted and ignored for constructor compatibility. Create the
/// exception with its documented message/diagnostic parameters; calling __init__
/// again does not update them. Returns None.
#[pyfunction(signature = (*_args, **_kwargs))]
fn initialize(_args: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) {}

fn allocate<'py>(
    py: Python<'py>,
    cls: &Bound<'py, PyType>,
    message: &str,
) -> PyResult<Bound<'py, PyAny>> {
    let value = py
        .get_type::<PyException>()
        .getattr("__new__")?
        .call1((cls, message))?;
    value.setattr("message", message)?;
    Ok(value)
}

macro_rules! constructor {
    ($function:ident, $($field:ident: $kind:ty),* $(,)?) => {
        // The public constructor keeps its documented positional parameters.
        /// Create an exception with message and the class's optional diagnostics.
        ///
        /// message is a string; remaining parameters are typed scalar attributes
        /// defaulting to None. See the class documentation for names and units.
        /// Sets args to (message,) and retains a separate message attribute.
        /// Scalar conversion can raise TypeError/OverflowError; no physical
        /// validation, file lookup or calculation is performed.
        #[allow(clippy::too_many_arguments)]
        #[pyfunction(signature = (cls, message, $($field=None),*))]
        fn $function<'py>(
            py: Python<'py>,
            cls: &Bound<'py, PyType>,
            message: String,
            $($field: Option<$kind>),*
        ) -> PyResult<Bound<'py, PyAny>> {
            let value = allocate(py, cls, &message)?;
            $(value.setattr(stringify!($field), $field)?;)*
            Ok(value)
        }
    };
}

constructor!(coordinate_new, coordinate_name: String, value: f64, min_value: f64, max_value: f64);
constructor!(invalid_body_new, body_id: i64, body_name: String, operation: String);
constructor!(unknown_body_new, body_id: i64);
constructor!(star_not_found_new, star_id: String, search_type: String);
constructor!(spk_not_found_new, filepath: String, body_name: String, body_id: String);
constructor!(spk_required_new, body_id: i64, body_name: String);
constructor!(polar_circle_new, latitude: f64, threshold: f64, obliquity: f64, house_system: String);
constructor!(ephemeris_range_new, requested_jd: f64, start_jd: f64, end_jd: f64,
    start_date: String, end_date: String, body_id: i64, body_name: String, ephemeris_file: String);
constructor!(convergence_new, algorithm: String, iterations: i64, max_iterations: i64,
    tolerance: f64, last_value: f64);
constructor!(spk_evaluation_new, body_id: i64, naif_id: i64, spk_file: String, requested_jd: f64);
constructor!(configuration_new, missing_config: String, suggestion: String);

/// Construct a sealed-policy exception with an optional purpose string.
///
/// Stores message/args and fixed network-policy configuration attributes as
/// documented on NetworkSealedError. Does not attempt network access.
#[pyfunction(signature = (cls, message, purpose=None))]
fn network_sealed_new<'py>(
    py: Python<'py>,
    cls: &Bound<'py, PyType>,
    message: String,
    purpose: Option<String>,
) -> PyResult<Bound<'py, PyAny>> {
    let value = allocate(py, cls, &message)?;
    value.setattr("missing_config", "network_policy")?;
    value.setattr(
        "suggestion",
        "Provision data explicitly or select network_policy='allow'.",
    )?;
    value.setattr("purpose", purpose)?;
    Ok(value)
}

/// Create SPKNotFoundError with a download-hint message and supplied attributes.
///
/// filepath is a string; body_name/body_id are optional strings (not integer
/// library IDs). Returns a new exception without checking a file or downloading.
/// Invalid argument types raise TypeError through ordinary scalar conversion.
#[pyfunction(signature = (cls, filepath, body_name=None, body_id=None))]
fn from_filepath<'py>(
    cls: &Bound<'py, PyType>,
    filepath: String,
    body_name: Option<String>,
    body_id: Option<String>,
) -> PyResult<Bound<'py, PyAny>> {
    cls.call1((
        format!("SPK file not found: {filepath}. Download the kernel with download_spk or enable automatic SPK downloading."),
        filepath, body_name, body_id,
    ))
}

/// Create SPKRequiredError with library ID, name and a Horizons download hint.
///
/// body_id is an integer; body_name and horizons_id are strings. The Horizons
/// ID appears only in the generated message, not a separate attribute. Returns
/// a new exception without network/file I/O. TypeError/OverflowError can arise
/// during argument conversion; no target validation is performed.
#[pyfunction(signature = (cls, body_id, body_name, horizons_id))]
fn for_body<'py>(
    cls: &Bound<'py, PyType>,
    body_id: i64,
    body_name: String,
    horizons_id: String,
) -> PyResult<Bound<'py, PyAny>> {
    cls.call1((
        format!("SPK kernel required for {body_name} in strict precision mode. Download and register Horizons target {horizons_id}."),
        body_id, body_name,
    ))
}

fn define<'py>(
    module: &Bound<'py, PyModule>,
    name: &str,
    bases: &[Bound<'py, PyType>],
    constructor: Option<Bound<'py, PyCFunction>>,
    doc: &str,
) -> PyResult<Bound<'py, PyType>> {
    let py = module.py();
    let namespace = PyDict::new(py);
    namespace.set_item("__module__", "libephemeris.exceptions")?;
    let signature = if let Some(constructor) = &constructor {
        let signature = constructor
            .getattr("__text_signature__")?
            .extract::<String>()?;
        signature
            .strip_prefix("(cls, ")
            .map(|parameters| format!("({parameters}"))
    } else if let Some(base) = bases.first() {
        base.getattr("__text_signature__")?
            .extract::<Option<String>>()?
    } else {
        None
    };
    if let Some(signature) = signature {
        // CPython exposes the signature prefix of a native type's docstring
        // through __text_signature__. Use the actual typed callable metadata.
        namespace.set_item("__doc__", format!("{name}{signature}\n--\n\n{doc}"))?;
    } else {
        namespace.set_item("__doc__", doc)?;
    }
    if let Some(constructor) = constructor {
        let staticmethod = py.import("builtins")?.getattr("staticmethod")?;
        namespace.set_item("__new__", staticmethod.call1((constructor,))?)?;
        // A builtin native initializer leaves constructor reflection to the
        // precise class text signature, rather than a generic callable object.
        namespace.set_item("__init__", wrap_pyfunction!(initialize, module)?)?;
    }
    let bases = PyTuple::new(py, bases)?;
    let class = py
        .get_type::<PyType>()
        .call1((name, bases, namespace))?
        .cast_into::<PyType>()?;
    module.add(name, &class)?;
    Ok(class)
}

/// Check numeric jd against selected resource metadata; return None if admitted.
///
/// body_id is an optional library body ID; None uses the active image's range.
/// No timescale conversion occurs: use TT for stored LEB coverage and the main
/// SPK's TDB clock for its metadata. An absent LEB body is left to its other
/// provider/analytical route, not rejected solely by this check. Other modes may
/// acquire the selected main image; configuration, file/policy errors propagate.
/// Raises EphemerisRangeError for nonfinite or outside-range jd, with known
/// requested_jd/body/bounds attributes; a successful check is not proof that
/// every target/center chain can be evaluated. func_name is accepted but unused.
#[pyfunction(signature=(jd, body_id=None, func_name=""))]
fn validate_jd_range(
    py: Python<'_>,
    jd: f64,
    body_id: Option<i32>,
    func_name: &str,
) -> PyResult<()> {
    let _ = func_name;
    match py.detach(|| crate::positions::runtime().validate_jd_range(jd, body_id)) {
        Ok(()) => Ok(()),
        Err(error @ (CoreError::NonFiniteArgument { .. } | CoreError::EphemerisRange { .. })) => {
            let kwargs = PyDict::new(py);
            kwargs.set_item("requested_jd", jd)?;
            kwargs.set_item("body_id", body_id)?;
            kwargs.set_item(
                "body_name",
                body_id.map(libephemeris_core::python_support::get_planet_name),
            )?;
            if let CoreError::EphemerisRange {
                start_jd, end_jd, ..
            } = error
            {
                kwargs.set_item("start_jd", start_jd)?;
                kwargs.set_item("end_jd", end_jd)?;
            }
            let value = py
                .import("libephemeris._core")?
                .getattr("EphemerisRangeError")?
                .call((error.to_string(),), Some(&kwargs))?;
            Err(PyErr::from_value(value))
        }
        Err(error) => Err(to_pyerr(py, error)),
    }
}

/// Register the public exception hierarchy, including both multiple bases.
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(validate_latitude, module)?)?;
    module.add_function(wrap_pyfunction!(validate_longitude, module)?)?;
    module.add_function(wrap_pyfunction!(validate_coordinates, module)?)?;
    let py = module.py();
    let error = define(
        module,
        "Error",
        &[py.get_type::<PyException>()],
        None,
        concat!(
            "Base class for public ephemeris failures; inherits Exception.\n\n",
            "Catch Error for library-defined failures, but ordinary argument/native\n",
            "conversion failures can instead use built-in exceptions.\n",
            "Typed subclasses take a required message string plus optional scalar\n",
            "diagnostics, set args=(message,) and expose message. Diagnostics default\n",
            "to None and can remain None when a producer lacks that information.\n",
            "They are descriptive values, not physically validated inputs; scalar\n",
            "conversion can raise TypeError or OverflowError. Base Error and category\n",
            "classes use ordinary Exception arguments and do not promise .message.\n",
            "LEBCorruptionError and FictitiousRuntimeDispatch are separate built-in\n",
            "subclasses, not subclasses of Error. Native cause graphs are not replayed."
        ),
    )?;
    let input = define(
        module,
        "InputValidationError",
        std::slice::from_ref(&error),
        None,
        "Invalid calculation input; category base inheriting Error.\n\nPass ordinary Exception arguments. Specific subclasses provide typed optional\ndiagnostics; see CoordinateError and InvalidBodyError.",
    )?;
    define(
        module,
        "CoordinateError",
        std::slice::from_ref(&input),
        Some(wrap_pyfunction!(coordinate_new, module)?),
        concat!(
            "Invalid geographic coordinate; inherits InputValidationError.\n\n",
            "Parameters/attributes\n",
            "---------------------\n",
            "message : str\n    Human-readable failure; also the sole args item.\n",
            "coordinate_name : str or None\n    Component name, e.g. latitude/longitude/height_m.\n",
            "value, min_value, max_value : float or None\n",
            "    Reported value and bounds in degrees for latitude/longitude,\n",
            "    metres for height. Unavailable bounds remain None. All optional\n",
            "    attributes default to None; construction does not validate them.\n\n",
            "Observer admission uses latitude [-90, 90], longitude [-180, 360],\n",
            "and finite height. This exception also reports nonfinite coordinates."
        ),
    )?;
    define(
        module,
        "InvalidBodyError",
        &[input],
        Some(wrap_pyfunction!(invalid_body_new, module)?),
        concat!(
            "A body is unsupported by the requested operation; inherits InputValidationError.\n\n",
            "message is a required string. Optional attributes default to None:\n",
            "body_id (integer library ID, not NAIF), body_name (string), operation\n",
            "(string identifying the refused calculation). Constructor values are\n",
            "stored without checking whether a body/operation exists."
        ),
    )?;
    let data = define(
        module,
        "DataNotFoundError",
        std::slice::from_ref(&error),
        None,
        "Required ephemeris data is missing or unusable; category base inheriting Error.\n\nUses ordinary Exception arguments. This category can include malformed data,\nnot just nonexistent files. Specific subclasses supply typed diagnostics.",
    )?;
    let unknown = define(
        module,
        "UnknownBodyError",
        std::slice::from_ref(&data),
        Some(wrap_pyfunction!(unknown_body_new, module)?),
        concat!(
            "Unknown/unsupported body identifier; inherits DataNotFoundError.\n\n",
            "message is a required string. body_id is an optional integer library\n",
            "body ID (default None), not a NAIF code. No catalog lookup occurs in\n",
            "the constructor. IllegalBodyError inherits this constructor."
        ),
    )?;
    define(
        module,
        "IllegalBodyError",
        &[unknown, py.get_type::<PyValueError>()],
        None,
        "Target cannot be placed by an active backend.\n\nInherits both UnknownBodyError and ValueError. Uses the inherited required\nmessage string and optional integer library body_id attribute (default None).\nConstructing the exception does not calculate or resolve a target.",
    )?;
    define(
        module,
        "StarNotFoundError",
        std::slice::from_ref(&data),
        Some(wrap_pyfunction!(star_not_found_new, module)?),
        concat!(
            "Fixed-star catalog search did not resolve; inherits DataNotFoundError.\n\n",
            "message is a required string. star_id and search_type are optional\n",
            "strings defaulting to None: the requested star/search key and search\n",
            "method label. They are stored without performing a catalog search."
        ),
    )?;
    let missing = define(
        module,
        "SPKNotFoundError",
        std::slice::from_ref(&data),
        Some(wrap_pyfunction!(spk_not_found_new, module)?),
        concat!(
            "Requested SPK kernel is missing; inherits DataNotFoundError.\n\n",
            "message is a required string. filepath, body_name and body_id are\n",
            "optional strings defaulting to None; body_id is a descriptive/Horizons\n",
            "identifier here, not the integer library ID used by SPKRequiredError.\n",
            "from_filepath() supplies a download-hint message. Neither constructor\n",
            "checks filesystem contents or downloads anything."
        ),
    )?;
    let required = define(
        module,
        "SPKRequiredError",
        &[data],
        Some(wrap_pyfunction!(spk_required_new, module)?),
        concat!(
            "Selected precision requires a usable SPK; inherits DataNotFoundError.\n\n",
            "message is a required string. body_id (integer library ID) and\n",
            "body_name (string) are optional attributes defaulting to None. This\n",
            "can report a registered kernel outside usable coverage as well as a\n",
            "missing registration. for_body() builds a Horizons download hint\n",
            "without performing acquisition or storing a horizons_id attribute."
        ),
    )?;
    let classmethod = py.import("builtins")?.getattr("classmethod")?;
    missing.setattr(
        "from_filepath",
        classmethod.call1((wrap_pyfunction!(from_filepath, module)?,))?,
    )?;
    required.setattr(
        "for_body",
        classmethod.call1((wrap_pyfunction!(for_body, module)?,))?,
    )?;
    let calculation = define(
        module,
        "CalculationError",
        std::slice::from_ref(&error),
        None,
        "Numerical calculation/geometry failure; category base inheriting Error.\n\nUses ordinary Exception arguments. Range, polar geometry and convergence\nsubclasses provide typed optional diagnostics; this class has no fixed units.",
    )?;
    define(
        module,
        "PolarCircleError",
        std::slice::from_ref(&calculation),
        Some(wrap_pyfunction!(polar_circle_new, module)?),
        concat!(
            "Requested house geometry is undefined at this latitude; inherits CalculationError.\n\n",
            "message is a required string. Optional attributes default to None:\n",
            "latitude, threshold and obliquity (floats in degrees), house_system\n",
            "(string system code). Construction stores diagnostics without checking\n",
            "geometry or selecting a fallback house system."
        ),
    )?;
    define(
        module,
        "EphemerisRangeError",
        std::slice::from_ref(&calculation),
        Some(wrap_pyfunction!(ephemeris_range_new, module)?),
        concat!(
            "Requested epoch is outside ephemeris coverage; inherits CalculationError.\n\n",
            "message is a required string. Optional attributes default to None:\n",
            "requested_jd, start_jd, end_jd (floats in Julian days); start_date,\n",
            "end_date (string labels); body_id (integer library ID); body_name and\n",
            "ephemeris_file (strings). Clock scale belongs to the producing API:\n",
            "stored LEB bounds use TT, SPK bounds TDB; no universal UTC/TT scale\n",
            "or constructor conversion is implied. Some producers populate only\n",
            "the message. validate_jd_range() also uses this class for nonfinite jd."
        ),
    )?;
    define(
        module,
        "ConvergenceError",
        std::slice::from_ref(&calculation),
        Some(wrap_pyfunction!(convergence_new, module)?),
        concat!(
            "Numerical solver did not converge; inherits CalculationError.\n\n",
            "message is a required string. algorithm is an optional string;\n",
            "iterations/max_iterations optional integers; tolerance/last_value\n",
            "optional floats. All default to None. Numeric units depend on the\n",
            "reported algorithm (not uniformly days or degrees). These attributes\n",
            "describe a failed solve; constructing the exception performs no solve."
        ),
    )?;
    define(
        module,
        "SPKEvaluationError",
        &[calculation],
        Some(wrap_pyfunction!(spk_evaluation_new, module)?),
        concat!(
            "SPK data could not be evaluated; inherits CalculationError.\n\n",
            "Includes corrupt coefficients/directories, unsupported types/frames\n",
            "and ambiguous center chains, not just outside-coverage requests.\n",
            "message is a required string. Optional attributes default to None:\n",
            "body_id (integer library ID), naif_id (integer NAIF target), spk_file\n",
            "(string path), requested_jd (float, SPK TDB JD when supplied). Native\n",
            "producers may supply none of these details beyond the message."
        ),
    )?;
    let configuration = define(
        module,
        "ConfigurationError",
        &[error],
        Some(wrap_pyfunction!(configuration_new, module)?),
        concat!(
            "Required configuration/resource policy is unavailable; inherits Error.\n\n",
            "message is a required string. missing_config and suggestion are\n",
            "optional strings defaulting to None, identifying a setting and a\n",
            "possible remedy. Construction does not apply the suggested change."
        ),
    )?;
    define(
        module,
        "NetworkSealedError",
        &[configuration, py.get_type::<PyRuntimeError>()],
        Some(wrap_pyfunction!(network_sealed_new, module)?),
        concat!(
            "Network access refused before network I/O by sealed policy.\n\n",
            "Inherits both ConfigurationError and RuntimeError. message is a\n",
            "required string; purpose is an optional string defaulting to None.\n",
            "Sets message, args=(message,), missing_config='network_policy', and\n",
            "suggestion=\"Provision data explicitly or select network_policy='allow'.\"\n",
            "Creating the exception does not attempt a connection or change policy.\n",
            "Filesystem preparation may already have occurred before a network gate."
        ),
    )?;
    define(
        module,
        "LEBCorruptionError",
        &[py.get_type::<PyValueError>()],
        None,
        "Coefficient image is corrupt or truncated.\n\nDirect ValueError subclass, not an Error subclass. Uses ordinary built-in\nexception arguments; it has no typed diagnostic-attribute constructor.",
    )?;
    define(
        module,
        "FictitiousRuntimeDispatch",
        &[py.get_type::<pyo3::exceptions::PyKeyError>()],
        None,
        "Compatibility analytical-dispatch signal.\n\nDirect KeyError subclass, not an Error subclass. Uses ordinary KeyError\narguments and its quoted-string formatting; no typed body/epoch attributes\nare added. Native calculation selection uses typed routes, not this signal.",
    )?;
    module.add_function(wrap_pyfunction!(validate_jd_range, module)?)?;
    Ok(())
}

fn mapped_error(py: Python<'_>, error: &CoreError) -> PyResult<PyErr> {
    let error = error.primary_error();
    let message = error.to_string();
    let kwargs = PyDict::new(py);
    let class = match error {
        CoreError::InvalidObserver { field, value } => {
            let field = match *field {
                "lat_deg" => "latitude",
                "lon_deg" => "longitude",
                field => field,
            };
            kwargs.set_item("coordinate_name", field)?;
            kwargs.set_item("value", value)?;
            match field {
                "lat" | "latitude" => {
                    kwargs.set_item("min_value", -90.0)?;
                    kwargs.set_item("max_value", 90.0)?;
                }
                "lon" | "longitude" => {
                    kwargs.set_item("min_value", -180.0)?;
                    kwargs.set_item("max_value", 360.0)?;
                }
                _ => {}
            }
            "CoordinateError"
        }
        CoreError::UnknownBody { body_id } | CoreError::LebBodyNotFound { body_id } => {
            kwargs.set_item("body_id", body_id)?;
            "UnknownBodyError"
        }
        CoreError::InvalidBody { body_id, operation } => {
            kwargs.set_item("body_id", body_id)?;
            kwargs.set_item("operation", operation)?;
            "InvalidBodyError"
        }
        CoreError::StarNotFound {
            star_id,
            search_type,
        } => {
            kwargs.set_item("star_id", star_id)?;
            kwargs.set_item("search_type", search_type)?;
            "StarNotFoundError"
        }
        CoreError::SpkFileNotFound { path } => {
            kwargs.set_item("filepath", path)?;
            "SPKNotFoundError"
        }
        CoreError::SpkRequired { body } | CoreError::SpkRegisteredOutOfCoverage { body, .. } => {
            kwargs.set_item("body_id", body)?;
            "SPKRequiredError"
        }
        CoreError::EphemerisRange {
            requested_jd,
            start_jd,
            end_jd,
            body,
        } => {
            kwargs.set_item("requested_jd", requested_jd)?;
            kwargs.set_item("start_jd", start_jd)?;
            kwargs.set_item("end_jd", end_jd)?;
            kwargs.set_item("body_id", body)?;
            "EphemerisRangeError"
        }
        CoreError::LebOutOfRange {
            body_id,
            jd,
            start,
            end,
        } => {
            kwargs.set_item("requested_jd", jd)?;
            kwargs.set_item("start_jd", start)?;
            kwargs.set_item("end_jd", end)?;
            kwargs.set_item("body_id", body_id)?;
            "EphemerisRangeError"
        }
        CoreError::SpkOutOfRange { .. } => "EphemerisRangeError",
        CoreError::PolarCircle {
            latitude,
            threshold,
            obliquity,
            house_system,
        } => {
            kwargs.set_item("latitude", latitude)?;
            kwargs.set_item("threshold", threshold)?;
            kwargs.set_item("obliquity", obliquity)?;
            kwargs.set_item("house_system", house_system.to_string())?;
            "PolarCircleError"
        }
        CoreError::ConvergenceFailed {
            algorithm,
            iterations,
            max_iterations,
            tolerance,
            last_value,
        } => {
            kwargs.set_item("algorithm", algorithm)?;
            kwargs.set_item("iterations", iterations)?;
            kwargs.set_item("max_iterations", max_iterations)?;
            kwargs.set_item("tolerance", tolerance)?;
            kwargs.set_item("last_value", last_value)?;
            "ConvergenceError"
        }
        CoreError::SpkCorrupt { .. }
        | CoreError::SpkUnsupportedType { .. }
        | CoreError::SpkUnsupportedFrame { .. }
        | CoreError::SpkAmbiguousCenter { .. } => "SPKEvaluationError",
        CoreError::ConfigurationMissing {
            missing_config: "jpl_source_access",
            ..
        } => "Error",
        CoreError::ConfigurationMissing {
            missing_config,
            suggestion,
        } => {
            kwargs.set_item("missing_config", missing_config)?;
            kwargs.set_item("suggestion", suggestion)?;
            "ConfigurationError"
        }
        CoreError::NetworkSealed { purpose } => {
            kwargs.set_item("purpose", purpose)?;
            "NetworkSealedError"
        }
        CoreError::InvalidHouseSystem { .. } => return Ok(PyTypeError::new_err(message)),
        CoreError::InvalidHouseCuspIndex { index } => {
            return Ok(pyo3::exceptions::PyKeyError::new_err(*index));
        }
        CoreError::InvalidNutationModel { .. } => return Ok(PyValueError::new_err(message)),
        CoreError::EphemerisLoaderIo {
            kind: std::io::ErrorKind::NotFound,
            ..
        } => return Ok(PyFileNotFoundError::new_err(message)),
        CoreError::NBodyBackendUnavailable { .. } => return Ok(PyImportError::new_err(message)),
        CoreError::NBodyNative { .. } => return Ok(PyRuntimeError::new_err(message)),
        CoreError::IntegerOverflow | CoreError::NumericalOverflow { .. } => {
            return Ok(PyOverflowError::new_err(message));
        }
        CoreError::MathDomain { .. } => return Ok(PyValueError::new_err(message)),
        CoreError::RationalDivisionByZero => return Ok(PyZeroDivisionError::new_err(message)),
        CoreError::NonFiniteArgument { .. } => return Ok(PyValueError::new_err(message)),
        CoreError::InvalidEclipseSearch {
            entry: "geopos" | "planet_occult",
            ..
        } => return Ok(PyValueError::new_err(message)),
        CoreError::EclipseSearchExhausted {
            entry: "planet_occult_when_glob" | "planet_occult_when_loc",
            ..
        } => "Error",
        CoreError::InvalidCalcMode { .. }
        | CoreError::InvalidPadding { .. }
        | CoreError::InvalidOrbitElements { .. }
        | CoreError::InvalidEclipseType { .. }
        | CoreError::InvalidDeltaTModel { .. }
        | CoreError::InvalidTier { .. }
        | CoreError::InvalidNetworkPolicy { .. }
        | CoreError::InvalidCalendar { .. }
        | CoreError::InvalidCalendarCode { .. }
        | CoreError::NoLeapSecond { .. }
        | CoreError::BodyHasNoStations { .. }
        | CoreError::InvalidStarSearch { .. }
        | CoreError::UnknownIntegrator { .. } => return Ok(PyValueError::new_err(message)),
        CoreError::InvalidDate { .. }
        | CoreError::InvalidTime { .. }
        | CoreError::InvalidOrbitalObject { .. }
        | CoreError::InvalidEclipseSearch { .. }
        | CoreError::PublicCrossing { .. }
        | CoreError::UnsupportedFlagCombination { .. }
        | CoreError::StarMagnitudeUnavailable { .. } => "Error",
        _ => match error.category() {
            ErrorCategory::InputValidation => "InputValidationError",
            ErrorCategory::DataNotFound => "DataNotFoundError",
            ErrorCategory::Calculation => "CalculationError",
            ErrorCategory::Configuration => "ConfigurationError",
        },
    };
    let module = py.import("libephemeris._core")?;
    let value = module.getattr(class)?.call((message,), Some(&kwargs))?;
    Ok(PyErr::from_value(value))
}

/// Convert a native primary failure into a plain public Python exception.
///
/// Existing cause/context envelopes are intentionally not replayed. The
/// selected class and its documented scalar attributes remain available.
pub(crate) fn to_pyerr(py: Python<'_>, error: CoreError) -> PyErr {
    mapped_error(py, &error).unwrap_or_else(|error| error)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numeric_conversion_errors_use_standard_python_classes() {
        Python::initialize();
        Python::attach(|py| {
            let not_finite = to_pyerr(py, CoreError::NonFiniteArgument { value: f64::NAN });
            assert!(not_finite.is_instance_of::<PyValueError>(py));
            let overflow = to_pyerr(py, CoreError::IntegerOverflow);
            assert!(overflow.is_instance_of::<PyOverflowError>(py));
            let wrapped = to_pyerr(
                py,
                CoreError::IntegerOverflow.with_cause(CoreError::SingularMatrix),
            );
            assert!(wrapped.is_instance_of::<PyOverflowError>(py));
        });
    }

    #[test]
    fn planetary_exhaustion_does_not_flatten_other_native_failures() {
        Python::initialize();
        Python::attach(|py| {
            let modules = py.import("sys").unwrap().getattr("modules").unwrap();
            let previous = modules.get_item("libephemeris._core").ok();
            let previous_parent = modules.get_item("libephemeris").ok();
            // Own both namespace levels: the exception-only module must not
            // trigger unrelated re-exports from an installed or local wheel.
            let parent = PyModule::new(py, "libephemeris").unwrap();
            parent
                .add("__path__", pyo3::types::PyList::empty(py))
                .unwrap();
            let module = PyModule::new(py, "libephemeris._core").unwrap();
            register(&module).unwrap();
            parent.add("_core", &module).unwrap();
            modules.set_item("libephemeris", &parent).unwrap();
            modules.set_item("libephemeris._core", &module).unwrap();
            let checks = [
                (
                    CoreError::EclipseSearchExhausted {
                        entry: "planet_occult_when_glob",
                        detail: "no event".into(),
                    },
                    "Error",
                ),
                (
                    CoreError::EclipseSearchExhausted {
                        entry: "planet_occult_when_loc",
                        detail: "no event".into(),
                    },
                    "Error",
                ),
                (
                    CoreError::EclipseSearchExhausted {
                        entry: "other_event",
                        detail: "no event".into(),
                    },
                    "CalculationError",
                ),
                (
                    CoreError::ConvergenceFailed {
                        algorithm: "physical root",
                        iterations: 20,
                        max_iterations: 20,
                        tolerance: 1e-8,
                        last_value: 1.0,
                    },
                    "ConvergenceError",
                ),
                (
                    CoreError::SpkFileNotFound {
                        path: "absent.bsp".into(),
                    },
                    "SPKNotFoundError",
                ),
                (
                    CoreError::UnknownBody { body_id: -999999 },
                    "UnknownBodyError",
                ),
            ];
            let matches = checks.map(|(error, expected)| {
                let actual = to_pyerr(py, error).get_type(py);
                actual.is(module.getattr(expected).unwrap())
            });
            match previous {
                Some(previous) => modules.set_item("libephemeris._core", previous).unwrap(),
                None => modules.del_item("libephemeris._core").unwrap(),
            }
            match previous_parent {
                Some(previous) => modules.set_item("libephemeris", previous).unwrap(),
                None => modules.del_item("libephemeris").unwrap(),
            }
            assert!(matches.into_iter().all(|matched| matched));
        });
    }
}
