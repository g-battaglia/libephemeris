"""Native REBOUND and ASSIST propagation with owned typed records."""

from ._core import (
    AssistEphemConfig,
    NbodySpkSource,
    PropagationResult,
    ReboundIntegrator,
    check_assist_available,
    check_rebound_available,
    compare_with_keplerian,
    get_assist_version,
    get_rebound_version,
    propagate_orbit_assist,
    propagate_orbit_rebound,
    propagate_trajectory,
)

__all__ = [
    "AssistEphemConfig", "NbodySpkSource", "PropagationResult", "ReboundIntegrator",
    "check_assist_available", "check_rebound_available", "compare_with_keplerian",
    "get_assist_version", "get_rebound_version", "propagate_orbit_assist",
    "propagate_orbit_rebound", "propagate_trajectory",
]

from ._core import download_assist_data, check_assist_data_available, reset_assist_data_cache, elements_to_rebound_orbit
