# SPDX-License-Identifier: AGPL-3.0
"""Read the project LEB1 layout and convert it through the standalone writer."""

from pathlib import Path

from . import format as wire
from .groups import LEB2_GROUPS, select_group
from .read import read_leb1
from .writer import write_leb2


def convert_leb1(source: Path, output: Path, *, lossless: bool = True, target: float | None = None,
                 chunk_interval_days: float = wire.CHUNK_INTERVAL_DAYS, overwrite: bool = False,
                 body_ids: set[int] | None = None, group: str | None = None,
                 tier: str | None = None, validator=None) -> None:
    """Convert selected bodies and preserve auxiliary sections without v3.

    Coefficients and auxiliary bytes are streamed from the original read-only
    file. Metadata bounds are checked before coefficient reads. No replacement
    source, analytical model, file download or astronomical calculation occurs.
    """
    header, bodies, auxiliary = read_leb1(source)
    available = {body.entry.body_id for body in bodies}
    if group is not None:
        if body_ids is not None:
            raise ValueError("body IDs and a named group cannot both be selected")
        body_ids = select_group(group, tier, available)
    if body_ids is not None and not body_ids <= available:
        raise ValueError("requested body inventory is absent from the source")
    bodies = [body for body in bodies if body_ids is None or body.entry.body_id in body_ids]
    write_leb2(output, bodies, header, auxiliary if group in (None, "core") else {},
               lossless=lossless, target=target,
               chunk_interval_days=chunk_interval_days, overwrite=overwrite, validator=validator)


def convert_all_leb1(source: Path, output_directory: Path, tier: str, **options) -> list[Path]:
    """Split every present complete named group into its standard tier filename.

    The entire inventory is checked before the first output. An absent group
    is omitted; a partially present group is an error, never a silently smaller
    artifact. Each published file receives the same native verification and
    atomic replacement guarantees as one-file conversion.
    """
    _, bodies, _ = read_leb1(source)
    available = {body.entry.body_id for body in bodies}
    groups = [name for name, identifiers in LEB2_GROUPS.items() if available & identifiers]
    if not groups:
        raise ValueError("source has no recognized named body group")
    for group in groups:
        select_group(group, tier, available)
    outputs = [Path(output_directory) / f"{tier}_{group}.leb2" for group in groups]
    if not options.get("overwrite", False):
        for output in outputs:
            if output.exists():
                raise FileExistsError(output)
    for group, output in zip(groups, outputs, strict=True):
        convert_leb1(source, output, group=group, tier=tier, **options)
    return outputs
