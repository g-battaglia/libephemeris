# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""PostgreSQL coefficient source for LibEphemeris."""

from __future__ import annotations

import os
from collections.abc import Mapping
from typing import Any

import libephemeris
from libephemeris.leb_format import BodyEntry, NutationHeader, StarEntry
from libephemeris.leb_groups import LEB2_GROUPS
from libephemeris.segment_source import SEGMENT_SOURCE_API_VERSION

from .config import RuntimeConfig, runtime_config
from .pool import get_pool
from .source import PostgresSegmentSource

__version__ = "0.1.0"
if SEGMENT_SOURCE_API_VERSION != 1:
    raise ImportError("unsupported libephemeris segment-source API version")


def _error(message: str) -> Exception:
    return libephemeris.CoefficientSourceError(message)


def _pins() -> Mapping[str, Mapping[str, Any]]:
    from libephemeris.download import DATA_FILES

    return DATA_FILES


def _artifact_reviewed(name: str, sha256: str) -> bool:
    return sha256 == _pins().get(name, {}).get("sha256")


def _reviewed(tier: str, rows: list[tuple[Any, ...]], override: str | None) -> bool:
    if override:
        expected = {f"{tier}_{group}.leb2" for group in LEB2_GROUPS}
        return {str(row[1]) for row in rows} == expected and all(
            _artifact_reviewed(str(row[1]), str(row[2])) for row in rows
        )
    return all(_artifact_reviewed(str(row[1]), str(row[2])) for row in rows)


def open_tier(
    tier: str, *, config: RuntimeConfig | None = None
) -> list[PostgresSegmentSource]:
    """Open all artifacts for a complete, selected dataset tier."""

    if tier not in {"base", "medium", "extended"}:
        raise _error(f"invalid tier: {tier}")
    config = config or runtime_config()
    override = os.environ.get(f"LIBEPHEMERIS_PG_DATASET_{tier.upper()}")
    pool = get_pool(config)
    try:
        with pool.connection() as conn, conn.cursor() as cur:
            cur.execute("SELECT version FROM libephemeris.schema_version")
            schema_versions = [int(row[0]) for row in cur.fetchall()]
            if schema_versions != [1]:
                raise _error("unsupported PostgreSQL coefficient schema version")
            if override:
                cur.execute(
                    "SELECT dataset_id, tier, complete FROM libephemeris.datasets "
                    "WHERE dataset_id = %s AND tier = %s",
                    (override, tier),
                )
            else:
                cur.execute(
                    "SELECT dataset_id, tier, complete FROM libephemeris.datasets "
                    "WHERE tier = %s AND complete ORDER BY created_at DESC",
                    (tier,),
                )
            datasets = cur.fetchall()
            chosen = None
            artifact_rows: list[tuple[Any, ...]] = []
            for dataset_id, _dataset_tier, complete in datasets:
                cur.execute(
                    "SELECT artifact_no, name, sha256, jd_start, jd_end, delta_t_jd, delta_t_val "
                    "FROM libephemeris.artifacts WHERE dataset_id=%s ORDER BY artifact_no",
                    (dataset_id,),
                )
                rows = cur.fetchall()
                names = {str(row[1]) for row in rows}
                if complete and names == {
                    f"{tier}_{group}.leb2" for group in LEB2_GROUPS
                }:
                    if override or _reviewed(tier, rows, None):
                        chosen, artifact_rows = str(dataset_id), rows
                        break
            if chosen is None:
                raise _error(f"no complete reviewed dataset available for {tier}")
            result: list[PostgresSegmentSource] = []
            for row in artifact_rows:
                artifact_no, name, _sha, jd_start, jd_end, dt_jd, dt_values = row
                cur.execute(
                    "SELECT body_id, coord_type, segment_count, jd_start, jd_end, interval_days, "
                    "degree, components, page_size FROM libephemeris.series "
                    "WHERE dataset_id=%s AND artifact_no=%s ORDER BY body_id",
                    (chosen, artifact_no),
                )
                series = cur.fetchall()
                bodies: dict[int, BodyEntry] = {}
                nutation = None
                page_rows: dict[int, int] = {}
                for values in series:
                    (
                        body_id,
                        coord_type,
                        count,
                        start,
                        end,
                        interval,
                        degree,
                        components,
                        page_size,
                    ) = values
                    page_rows[int(body_id)] = int(page_size)
                    if int(body_id) == -1:
                        nutation = NutationHeader(
                            float(start),
                            float(end),
                            float(interval),
                            int(degree),
                            int(components),
                            int(count),
                            0,
                        )
                    else:
                        bodies[int(body_id)] = BodyEntry(
                            int(body_id),
                            int(coord_type),
                            int(count),
                            float(start),
                            float(end),
                            float(interval),
                            int(degree),
                            int(components),
                            0,
                        )
                cur.execute(
                    "SELECT star_id, ra_j2000, dec_j2000, pm_ra, pm_dec, parallax, rv, magnitude "
                    "FROM libephemeris.stars WHERE dataset_id=%s AND artifact_no=%s",
                    (chosen, artifact_no),
                )
                stars = {
                    int(row[0]): StarEntry(int(row[0]), *(float(v) for v in row[1:]))
                    for row in cur.fetchall()
                }
                result.append(
                    PostgresSegmentSource(
                        dataset_id=chosen,
                        artifact_no=int(artifact_no),
                        artifact_name=str(name),
                        locator=f"postgres://{chosen}",
                        jd_range=(float(jd_start), float(jd_end)),
                        bodies=bodies,
                        nutation=nutation,
                        delta_t=(list(dt_jd), list(dt_values)),
                        stars=stars,
                        reviewed=_artifact_reviewed(str(name), str(_sha)),
                        config=config,
                        page_rows=page_rows,
                    )
                )
            return result
    except libephemeris.CoefficientSourceError:
        raise
    except Exception:
        raise _error("PostgreSQL coefficient source unavailable") from None


__all__ = ["PostgresSegmentSource", "SEGMENT_SOURCE_API_VERSION", "open_tier"]
