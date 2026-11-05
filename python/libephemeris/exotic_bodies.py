"""Native physical model and fixed public parameters."""

from ._core import (
    _module_exotic_bodies_EXOTIC_IDS as EXOTIC_IDS,
    _module_exotic_bodies_EXOTIC_EXTENDED_IDS as EXOTIC_EXTENDED_IDS,
    _module_exotic_bodies_EXOTIC_ASSIST_PERTURBER_IDS as EXOTIC_ASSIST_PERTURBER_IDS,
)

from ._core import ExoticBody, EXOTIC_BODIES, naif_map, body_params, target_au_map, name_map, exotic_display_name
