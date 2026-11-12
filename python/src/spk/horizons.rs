//! Plain vector records and the native Horizons transport boundary.

use crate::{errors::to_pyerr, positions::runtime};
use libephemeris_core::horizons;
use libephemeris_core::positions::horizons_serving::HorizonsOutcome;
use libephemeris_core::python_support::module_positions as native;
use pyo3::exceptions::PyKeyError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::sync::Arc;

type Position = (f64, f64, f64, f64, f64, f64);
type PositionResult = (Position, i32);

/// Mutable geometric Horizons state in ICRF axes, relative to the requested center.
///
/// x/y/z are AU and vx/vy/vz are AU/day. A vector fetched with center="@0"
/// is barycentric; vectors for other centers are not. No light-time,
/// aberration or apparent-frame corrections are included. The record stores
/// no epoch or center metadata: keep those alongside it. Manually constructed
/// records accept scalar values without physical/finite validation, and edits
/// to a returned record do not modify the client's cached native state.
#[derive(Clone, Debug)]
#[pyclass(
    module = "libephemeris.horizons_backend",
    get_all,
    set_all,
    skip_from_py_object
)]
struct StateVector {
    /// ICRF x position relative to the requested center, AU; writable.
    x: f64,
    /// ICRF y position relative to the requested center, AU; writable.
    y: f64,
    /// ICRF z position relative to the requested center, AU; writable.
    z: f64,
    /// ICRF x velocity relative to the requested center, AU/day; writable.
    vx: f64,
    /// ICRF y velocity relative to the requested center, AU/day; writable.
    vy: f64,
    /// ICRF z velocity relative to the requested center, AU/day; writable.
    vz: f64,
}

#[pymethods]
impl StateVector {
    /// Store x, y, z (AU) and vx, vy, vz (AU/day) without I/O or validation.
    /// All six arguments are required; no clock/frame conversion is performed.
    #[new]
    fn new(x: f64, y: f64, z: f64, vx: f64, vy: f64, vz: f64) -> Self {
        Self {
            x,
            y,
            z,
            vx,
            vy,
            vz,
        }
    }
    /// Position as a copied (x, y, z) tuple in AU, in the stored ICRF axes.
    #[getter]
    fn pos(&self) -> (f64, f64, f64) {
        (self.x, self.y, self.z)
    }
    /// Velocity as a copied (vx, vy, vz) tuple in AU/day.
    #[getter]
    fn vel(&self) -> (f64, f64, f64) {
        (self.vx, self.vy, self.vz)
    }
}

impl From<horizons::StateVector> for StateVector {
    fn from(value: horizons::StateVector) -> Self {
        Self::new(
            value.x(),
            value.y(),
            value.z(),
            value.vx(),
            value.vy(),
            value.vz(),
        )
    }
}

/// Synchronous native Horizons client with an owned bounded vector LRU cache.
///
/// Construction opens no connection. Fetches release the Python GIL and use
/// the process-default effective network policy; cache hits need no socket and
/// remain usable when sealed. Batch misses use bounded scoped native threads,
/// not Python workers/callbacks. clear_cache() and shutdown() only clear cached
/// values: neither permanently closes the client nor cancels/drains requests.
#[pyclass(module = "libephemeris.horizons_backend", frozen)]
struct HorizonsClient {
    native: Arc<horizons::HorizonsClient>,
}

#[pymethods]
impl HorizonsClient {
    /// Configure max_cache_size entries (4096), timeout seconds per attempt (30),
    /// and max_workers batch threads (8). All must be positive; zero capacities
    /// or nonpositive timeout raise InputValidationError. Negative/oversized
    /// unsigned capacities fail Python integer extraction. Retries can exceed
    /// timeout in total. No request or background worker starts here.
    #[new]
    #[pyo3(signature=(max_cache_size=4096, timeout=30, max_workers=8))]
    fn new(
        py: Python<'_>,
        max_cache_size: usize,
        timeout: i64,
        max_workers: usize,
    ) -> PyResult<Self> {
        py.detach(|| native::new_horizons_client(max_cache_size, timeout as f64, max_workers))
            .map(|client| Self {
                native: Arc::new(client),
            })
            .map_err(|error| to_pyerr(py, error))
    }

    /// Fetch one geometric ICRF state relative to center (default "@0").
    ///
    /// command is Horizons command text, not a library body ID. jd must be a
    /// finite Julian Day in time_type, exactly "TDB" (default) or "UT".
    /// Send the numeric JD unchanged: "TDB" asks the server to interpret TDB;
    /// no client TT-to-TDB conversion occurs. Per the JPL manual, "UT" means
    /// UT1 before 1962 and UTC from 1962 onward; future epochs use the closest
    /// known leap-second offset, not a prediction of future leap seconds.
    /// See <https://ssd.jpl.nasa.gov/horizons/manual.html>.
    /// Return an independent StateVector in AU/AU-day.
    ///
    /// Cache identity includes command, epoch, center and time_type. Cache
    /// misses may perform retried HTTP calls; no local model substitutes for
    /// a failed request. Nonfinite epochs raise ValueError; unsupported
    /// time_type raises InputValidationError. Sealed-policy, transport,
    /// API-refusal and malformed-response errors propagate.
    #[pyo3(signature=(command,jd,center="@0",time_type="TDB"))]
    fn fetch_state_vector(
        &self,
        py: Python<'_>,
        command: &str,
        jd: f64,
        center: &str,
        time_type: &str,
    ) -> PyResult<StateVector> {
        py.detach(|| {
            native::fetch_horizons_state(
                &self.native,
                command,
                jd,
                center,
                time_type,
                runtime().network_policy()?,
            )
        })
        .map(StateVector::from)
        .map_err(|error| to_pyerr(py, error))
    }

    /// Fetch (command, jd, center) triples and return successful states by triple.
    ///
    /// Every jd is a finite TDB-labelled Julian Day; there is no per-request
    /// time_type or TT-to-TDB conversion. Return a new dictionary mapping the
    /// original triples to independent StateVector records. Duplicate triples
    /// share dictionary keys; output order is not request order.
    ///
    /// Request conversion/validation and configuration errors abort the call.
    /// Individual fetch failures do not: successful requests remain in the
    /// result, and errors, if a dict, is updated with {triple: exception}.
    /// Existing unrelated entries are retained. With errors=None, failures are
    /// omitted silently from the result. Failures are not cached.
    ///
    /// The call releases the GIL, waits for its bounded native workers, then
    /// creates Python records/errors on the calling thread. No Python callback
    /// is accepted or run on workers. Cached successes require no network;
    /// uncached sealed requests become per-request failures.
    #[pyo3(signature=(requests,errors=None))]
    fn fetch_batch(
        &self,
        py: Python<'_>,
        requests: Vec<(String, f64, String)>,
        errors: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Py<PyDict>> {
        let (values, failures) = py
            .detach(|| {
                let requests = requests
                    .iter()
                    .map(|(command, jd, center)| {
                        horizons::VectorsRequest::new(command, *jd, center)
                    })
                    .collect::<libephemeris_core::errors::Result<Vec<_>>>()?;
                Ok::<_, libephemeris_core::errors::Error>(
                    self.native
                        .fetch_batch(&requests, runtime().network_policy()?),
                )
            })
            .map_err(|error| to_pyerr(py, error))?;
        let output = PyDict::new(py);
        for (key, value) in values {
            output.set_item(
                (key.command(), key.jd(), key.center()),
                StateVector::from(value),
            )?;
        }
        if let Some(errors) = errors {
            for (key, error) in failures {
                errors.set_item(
                    (key.command(), key.jd(), key.center()),
                    to_pyerr(py, error).value(py),
                )?;
            }
        }
        Ok(output.unbind())
    }

    /// Clear shared cached vectors; return None without changing configuration.
    /// Running requests are neither cancelled nor joined and may repopulate
    /// the cache after this call. The client remains usable.
    fn clear_cache(&self, py: Python<'_>) {
        py.detach(|| self.native.clear_cache());
    }
    /// Clear cached vectors, just like clear_cache(); return None.
    /// This is not a terminal close, cancellation request or worker drain.
    /// Subsequent fetches remain valid and in-flight work may refill the cache.
    fn shutdown(&self, py: Python<'_>) {
        py.detach(|| self.native.shutdown());
    }
}

/// Calculate a direct Horizons-backend place at UT1 Julian Day jd_ut.
///
/// client supplies the shared vector cache; body_id is a library ID, not a
/// NAIF command. iflag selects center/corrections/output. sid_mode=None uses
/// the default session's selection; an integer overrides it for this call.
/// Return ((a, b, distance, da, db, radial_rate), flags): normally ecliptic
/// degrees/AU and daily rates; EQUATORIAL, XYZ or RADIANS select the usual
/// representations. Rates require SPEED. Some analytical/hypothetical bodies
/// use local models rather than HTTP.
///
/// The native path forms scalar TT from configured Delta-T, then sends that
/// scalar with a TDB request label; it does not promise TT-to-TDB conversion.
/// Unsupported backend combinations (including TOPOCTR/ICRS and applicable
/// NONUT cases) raise KeyError instead of falling through to another provider.
/// Clock, frame, policy and remote errors propagate. Network/native work
/// releases the GIL; no Python callback is accepted by this helper.
#[pyfunction(signature=(client,jd_ut,body_id,iflag,sid_mode=None))]
fn horizons_calc_ut(
    py: Python<'_>,
    client: PyRef<'_, HorizonsClient>,
    jd_ut: f64,
    body_id: i32,
    iflag: i32,
    sid_mode: Option<i32>,
) -> PyResult<PositionResult> {
    let client = Arc::clone(&client.native);
    let outcome = py
        .detach(|| {
            runtime()
                .shared_session()?
                .horizons_backend_calc_ut(&client, jd_ut, body_id, iflag, sid_mode)
        })
        .map_err(|error| to_pyerr(py, error))?;
    match outcome {
        HorizonsOutcome::Served([a, b, c, d, e, f], flags) => Ok(((a, b, c, d, e, f), flags)),
        HorizonsOutcome::Fallback(message) => Err(PyKeyError::new_err(message)),
    }
}

/// Return the default runtime's shared HorizonsClient when its mode needs one.
///
/// Return None in forced LEB/native-JPL modes, or in auto when the selected
/// local main file/current main owner exists. Otherwise initialize lazily
/// without contacting Horizons. Configuration errors can raise. reset_session()
/// or close() retires the runtime's reference, not this returned client's
/// lifetime or ongoing requests; subsequent misses obey the current policy.
#[pyfunction]
fn get_horizons_client(py: Python<'_>) -> PyResult<Option<HorizonsClient>> {
    py.detach(|| runtime().get_horizons_client())
        .map(|client| client.map(|native| HorizonsClient { native }))
        .map_err(|error| to_pyerr(py, error))
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(get_horizons_client, module)?)?;
    module.add_class::<StateVector>()?;
    module.add_class::<HorizonsClient>()?;
    module.add_function(wrap_pyfunction!(horizons_calc_ut, module)?)?;
    module.add(
        "_HORIZONS_API_URL",
        libephemeris_core::spk::HORIZONS_API_URL,
    )?;
    Ok(())
}
