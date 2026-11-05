//! Process-default network-policy settings over native configuration.
use crate::{errors::to_pyerr, positions::runtime};
use libephemeris_core::context::EffectiveNetworkPolicy;
use pyo3::prelude::*;

/// Set the process-default network override to "auto", "allow" or "sealed".
///
/// Pass None to restore captured environment/file configuration. "auto" seals
/// LEB mode and allows other calculation modes; "allow" does not itself start
/// a request. Return None; invalid names raise ValueError without changing
/// the selection.
#[pyfunction(signature=(policy))]
fn set_network_policy(py: Python<'_>, policy: Option<&str>) -> PyResult<()> {
    py.detach(|| runtime().set_network_policy(policy))
        .map_err(|error| to_pyerr(py, error))
}
/// Return the configured policy name: "auto", "allow" or "sealed".
///
/// Unlike get_network_policy(), this preserves "auto" rather than resolving
/// it against the calculation mode. Configuration errors propagate; no I/O
/// or connection is initiated by this query.
#[pyfunction]
fn get_configured_network_policy(py: Python<'_>) -> PyResult<String> {
    py.detach(|| runtime().configured_network_policy())
        .map(|policy| policy.as_policy_name().to_owned())
        .map_err(|error| to_pyerr(py, error))
}
/// Return the effective network gate, either "allow" or "sealed".
///
/// A configured "auto" resolves to "sealed" in LEB mode and "allow"
/// otherwise. This is permission, not a connectivity check. Invalid captured
/// configuration raises an error without making a network request.
#[pyfunction]
fn get_network_policy(py: Python<'_>) -> PyResult<&'static str> {
    py.detach(|| runtime().network_policy())
        .map(|policy| match policy {
            EffectiveNetworkPolicy::Allow => "allow",
            EffectiveNetworkPolicy::Sealed => "sealed",
        })
        .map_err(|error| to_pyerr(py, error))
}
/// Check network permission for the diagnostic purpose string; return None.
///
/// Raise NetworkSealedError if sealed, or a configuration error if resolution
/// fails. This does not contact a server, reserve permission for later calls,
/// or test whether any particular URL is reachable.
#[pyfunction(signature=(purpose))]
fn require_network(py: Python<'_>, purpose: &str) -> PyResult<()> {
    py.detach(|| runtime().effective_context()?.require_network(purpose))
        .map_err(|error| to_pyerr(py, error))
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(set_network_policy, module)?)?;
    module.add_function(wrap_pyfunction!(get_configured_network_policy, module)?)?;
    module.add_function(wrap_pyfunction!(get_network_policy, module)?)?;
    module.add_function(wrap_pyfunction!(require_network, module)?)?;
    Ok(())
}
