"""Build configuration and the ordered dependency graph.

Provenance:
    Project workflow metadata derived from the canonical generator registries.
    Scientific parameters and date semantics remain in those generators.
"""

from __future__ import annotations

from argparse import Namespace
from dataclasses import asdict, dataclass
import os
from pathlib import Path
from typing import Any, Literal

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCHEMA_VERSION = 1

JobKind = Literal["generate", "verify1", "merge", "convert", "verify2"]


# ---------------------------------------------------------------------------
# Immutable configuration and jobs
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class BuildConfig:
    """Settings frozen for the lifetime of one output directory.

    All source directories are absolute. A resume uses the stored settings;
    changing tiers, sample counts or source locations requires a new build.
    """

    tiers: tuple[str, ...]
    data_dir: str
    spk_dir: str
    assist_dir: str
    verify_samples: int = 500
    leb2_verify_samples: int = 200

    def to_dict(self) -> dict[str, Any]:
        """Return JSON-compatible configuration."""
        return {**asdict(self), "tiers": list(self.tiers)}

    @classmethod
    def from_dict(cls, value: dict[str, Any]) -> BuildConfig:
        """Load and validate persisted configuration without coercing values."""
        from scripts.generate_leb import TIER_CONFIGS

        if not isinstance(value, dict) or set(value) != set(cls.__dataclass_fields__):
            raise ValueError("Invalid build configuration keys")
        tiers = value["tiers"]
        if (
            not isinstance(tiers, list)
            or not tiers
            or any(not isinstance(tier, str) for tier in tiers)
            or len(set(tiers)) != len(tiers)
            or any(tier not in TIER_CONFIGS for tier in tiers)
        ):
            raise ValueError("Invalid tiers in build configuration")
        for key in ("verify_samples", "leb2_verify_samples"):
            if type(value[key]) is not int or value[key] < 1:
                raise ValueError(f"{key} must be a positive integer")
        for key in ("data_dir", "spk_dir", "assist_dir"):
            if not isinstance(value[key], str) or not Path(value[key]).is_absolute():
                raise ValueError(f"{key} must be an absolute directory")
        return cls(**{**value, "tiers": tuple(tiers)})


@dataclass(frozen=True)
class Job:
    """One restartable phase, with paths relative to the build directory.

    Verification jobs reference an existing output and never write it. Export
    jobs are distinguished from per-body work files for the final checksum set.
    """

    id: str
    tier: str
    kind: JobKind
    output: str
    bodies: tuple[int, ...]
    dependencies: tuple[str, ...] = ()
    inputs: tuple[str, ...] = ()
    group: str | None = None
    aux_source: str | None = None
    export: bool = False

    @property
    def writes_output(self) -> bool:
        """Whether this job must write to a staging file."""
        return self.kind in ("generate", "merge", "convert")

    def to_dict(self) -> dict[str, Any]:
        """Return the exact job definition used to authenticate a resume."""
        value = asdict(self)
        for key in ("bodies", "dependencies", "inputs"):
            value[key] = list(value[key])
        return value


# ---------------------------------------------------------------------------
# Canonical inventories and date semantics
# ---------------------------------------------------------------------------


def tier_groups(tier: str) -> dict[str, tuple[int, ...]]:
    """Resolve LEB1 groups, including the declared extended-tier exclusions."""
    from libephemeris.exotic_bodies import EXOTIC_EXTENDED_IDS, EXOTIC_IDS
    from libephemeris.leb_groups import LEB1_GENERATION_GROUPS
    from scripts.generate_leb import BODY_GROUPS

    if tuple(BODY_GROUPS) != LEB1_GENERATION_GROUPS:
        raise ValueError("LEB1 registry mismatch; update the workflow explicitly")
    allowed = set(EXOTIC_EXTENDED_IDS if tier == "extended" else EXOTIC_IDS)
    return {
        group: tuple(b for b in bodies if group != "exotics" or b in allowed)
        for group, bodies in BODY_GROUPS.items()
    }


def tier_range(tier: str) -> tuple[float, float]:
    """Reuse the generator's margins and exact DE441 boundaries."""
    from scripts.generate_leb import _resolve_tier

    start, end, _ = _resolve_tier(
        Namespace(
            tier=tier,
            start=None,
            end=None,
            start_jd=None,
            end_jd=None,
            output="unused.leb",
        )
    )
    return start, end


def build_jobs(config: BuildConfig) -> list[Job]:
    """Return a topologically ordered plan for all selected canonical exports.

    Every body is generated and verified independently. Groups depend on those
    verification phases, and LEB2 conversion depends on the verified tier file.
    Auxiliary data is generated with Sun once, then copied to each LEB1 group.
    """
    from libephemeris.leb_groups import LEB1_GROUPS
    from scripts.generate_leb2 import LEB2_GROUPS

    jobs = []
    for tier in config.tiers:
        groups = tier_groups(tier)
        all_bodies = tuple(b for bodies in groups.values() for b in bodies)
        if len(set(all_bodies)) != len(all_bodies) or 0 not in groups["planets"]:
            raise ValueError("Invalid generation partition")
        aux = f"work/{tier}/body_0.leb"
        for body in all_bodies:
            output = f"work/{tier}/body_{body}.leb"
            generate_id = f"{tier}.body.{body}.generate"
            jobs.append(Job(generate_id, tier, "generate", output, (body,)))
            jobs.append(
                Job(
                    f"{tier}.body.{body}.verify",
                    tier,
                    "verify1",
                    output,
                    (body,),
                    (generate_id,),
                )
            )
        for group, bodies in groups.items():
            dependencies = tuple(
                dict.fromkeys(f"{tier}.body.{body}.verify" for body in (*bodies, 0))
            )
            jobs.append(
                Job(
                    f"{tier}.group.{group}.merge",
                    tier,
                    "merge",
                    f"leb/ephemeris_{tier}_{group}.leb",
                    tuple(sorted(bodies)),
                    dependencies,
                    tuple(f"work/{tier}/body_{b}.leb" for b in bodies),
                    group=group,
                    aux_source=aux,
                    export=True,
                )
            )
        merged = f"leb/ephemeris_{tier}.leb"
        merge_id = f"{tier}.merge"
        verify_id = f"{tier}.verify"
        jobs.append(
            Job(
                merge_id,
                tier,
                "merge",
                merged,
                tuple(sorted(all_bodies)),
                tuple(f"{tier}.group.{g}.merge" for g in LEB1_GROUPS),
                tuple(f"leb/ephemeris_{tier}_{g}.leb" for g in LEB1_GROUPS),
                export=True,
            )
        )
        jobs.append(
            Job(
                verify_id,
                tier,
                "verify1",
                merged,
                tuple(sorted(all_bodies)),
                (merge_id,),
            )
        )
        for group, registered in LEB2_GROUPS.items():
            bodies = tuple(sorted(set(registered) & set(all_bodies)))
            converted = f"leb2/{tier}_{group}.leb2"
            convert_id = f"{tier}.{group}.convert"
            jobs.append(
                Job(
                    convert_id,
                    tier,
                    "convert",
                    converted,
                    bodies,
                    (verify_id,),
                    (merged,),
                    group=group,
                    export=True,
                )
            )
            jobs.append(
                Job(
                    f"{tier}.{group}.verify2",
                    tier,
                    "verify2",
                    converted,
                    bodies,
                    (convert_id,),
                    (merged,),
                    group=group,
                )
            )
    return jobs


# ---------------------------------------------------------------------------
# Deterministic subprocess environment
# ---------------------------------------------------------------------------


def scientific_environment(config: BuildConfig, tier: str) -> dict[str, str]:
    """Disable user model overrides, automatic downloads and installed LEBs.

    The worker also applies public state setters because its imported package
    may otherwise retain process-global configuration. Other environment keys
    are inherited for normal Python/shared-library loading.
    """
    from scripts.generate_leb import TIER_CONFIGS

    env = {k: v for k, v in os.environ.items() if not k.startswith("LIBEPHEMERIS_")}
    env.update(
        {
            "PYTHONDONTWRITEBYTECODE": "1",
            "PYTHONUNBUFFERED": "1",
            "LIBEPHEMERIS_CONFIG": os.devnull,
            "LIBEPHEMERIS_ENV_FILE": os.devnull,
            "LIBEPHEMERIS_MODE": "skyfield",
            "LIBEPHEMERIS_NETWORK_POLICY": "sealed",
            "LIBEPHEMERIS_PRECISION": tier,
            "LIBEPHEMERIS_DATA_DIR": config.data_dir,
            "LIBEPHEMERIS_SPK_DIR": config.spk_dir,
            "LIBEPHEMERIS_EPHEMERIS": str(
                Path(config.data_dir) / TIER_CONFIGS[tier][0]
            ),
            "LIBEPHEMERIS_IERS_AUTO_DOWNLOAD": "0",
            "ASSIST_DIR": config.assist_dir,
        }
    )
    return env
