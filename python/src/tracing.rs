// SPDX-License-Identifier: AGPL-3.0
//! Native tracing sessions with explicit thread-owned reset tokens.

use libephemeris_core::trace;
use pyo3::prelude::*;
use std::collections::BTreeMap;

/// Thread-owned handle returned by start_tracing(); not directly constructible.
///
/// Call reset() explicitly on the calculating thread to close the innermost
/// trace. Keeping or discarding this object alone does not close its scope.
/// Nested traces retain independent maps and restore the parent on reset.
#[pyclass(module = "libephemeris.tracing", frozen, skip_from_py_object)]
pub(crate) struct TraceToken {
    native: trace::TraceToken,
}

#[pymethods]
impl TraceToken {
    /// Close this trace on its originating thread and return whether it closed.
    ///
    /// A successful reset restores the parent's results, or clears/deactivates
    /// the outermost trace. Repeated, stale, foreign-thread or out-of-order
    /// resets return False without changing the current trace.
    fn reset(&self, py: Python<'_>) -> bool {
        py.detach(|| self.native.reset())
    }
}

/// Start an empty computation-source trace on the calling native thread.
///
/// Return a token whose reset() closes this scope. Nested scopes do not merge
/// their results into their parent. Trace calculations on this same thread:
/// a token does not transfer tracing into a Python executor worker. Starting
/// a scope performs no calculation, file access or network acquisition.
/// For example: token = start_tracing(); ...; token.reset().
#[pyfunction]
fn start_tracing(py: Python<'_>) -> TraceToken {
    TraceToken {
        native: py.detach(trace::start_session),
    }
}

/// Return a copied {library_body_id: source_label} dictionary for this thread.
///
/// Each body's latest recorded successful computation replaces its previous
/// label. Keys are returned in ascending body-ID order. Inactive or empty
/// scopes return {}. Reading or mutating the copy neither clears the trace
/// nor changes provider selection; it performs no data acquisition.
#[pyfunction]
fn get_trace_results(py: Python<'_>) -> BTreeMap<i32, String> {
    py.detach(trace::get_trace_results)
}

/// Register the tracing functions and native reset token.
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("_TraceToken", module.py().get_type::<TraceToken>())?;
    module.add_function(wrap_pyfunction!(start_tracing, module)?)?;
    module.add_function(wrap_pyfunction!(get_trace_results, module)?)?;
    Ok(())
}
