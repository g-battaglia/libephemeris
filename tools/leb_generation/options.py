# SPDX-License-Identifier: AGPL-3.0
"""Shared explicit source selections for the standalone LEB generators."""

from pathlib import Path

from .groups import LEB2_GROUPS
from .source import ANGULAR_IDS, PLANET_NAIF, NativeSource, configure_spk_sources

TIER_KERNELS = {"base": "de440s.bsp", "medium": "de440.bsp", "extended": "de441.bsp"}
TIER_YEARS = {"base": (1850, 2150), "medium": (1550, 2650)}
BODY_GROUPS = {**LEB2_GROUPS, "planets": frozenset(PLANET_NAIF), "analytical": ANGULAR_IDS}
BODY_GROUPS["all"] = frozenset().union(*LEB2_GROUPS.values())


def source_arguments(parser, *, required_tier=False, include_tier=True):
    """Require local, identified inputs rather than implicit downloads or fallbacks."""
    if include_tier:
        parser.add_argument("-t", "--tier", choices=TIER_KERNELS,
                            default=None if required_tier else "base", required=required_tier)
    parser.add_argument("--kernel", type=Path, help="Exact existing main JPL image; no fallback/download")
    parser.add_argument("--data-dir", type=Path, help="Directory containing the tier's existing JPL image")
    parser.add_argument("--auto-sources", action="store_true",
                        help="Assemble missing minor/extended N-body sources from existing native metadata; no downloads")
    parser.add_argument("--spk-dir", type=Path, action="append", default=[],
                        help="With --auto-sources: explicit BSP directory (repeatable); default DATA_DIR/spk")
    parser.add_argument("--spk", nargs=3, metavar=("BODY_ID", "PATH", "NAIF_ID"), action="append", default=[],
                        help="Exact existing minor-body source; repeat for each body")
    parser.add_argument("--nbody", nargs="+", help="Body IDs/names whose explicit SPK is used only as a J2000 N-body seed")
    parser.add_argument("--assist-planets", type=Path, help="Exact existing ASSIST planet force file")
    parser.add_argument("--assist-asteroids", type=Path, help="Exact existing ASSIST asteroid force file")
    parser.add_argument("--delta-t-model", choices=("smh2016", "espenak_meeus"), default="smh2016")
    parser.add_argument("--delta-t-userdef-days", type=float,
                        help="Explicit TT minus UT1 in days; overrides IERS and the selected model")
    parser.add_argument("--iers-delta-t", action="store_true",
                        help="Opt in to native local IERS Delta-T; uncovered dates retain the selected model")
    parser.add_argument("--iers-cache-dir", type=Path,
                        help="With --iers-delta-t: explicit existing local clock directory; no downloads")
    parser.add_argument("--verify-samples", type=int, default=64,
                        help="Native-reader samples per body; off-grid fit checks are always mandatory")


def require_native_engine(engine):
    """Refuse an unrelated same-named Python installation before any output."""
    if engine.__version__ != "4.0.0b1":
        raise ValueError("offline LEB tools require the LibEphemeris 4.0.0b1 native package")


def fitting_arguments(parser):
    parser.add_argument("--start", type=int, help="Gregorian start year (January 1)")
    parser.add_argument("--end", type=int, help="Gregorian end year (January 1)")
    parser.add_argument("--start-jd", type=float, help="First requested TT Julian day")
    parser.add_argument("--end-jd", type=float, help="Last requested TT Julian day")
    parser.add_argument("--group", choices=BODY_GROUPS, help="Exact artifact inventory; default core")
    parser.add_argument("--bodies", nargs="+", help="Explicit IDs or public constant names")
    parser.add_argument("--skip-aux", action="store_true")
    parser.add_argument("--available-domains", action="store_true",
                        help="Store each body's complete continuous source overlap; never omit a body or bridge a gap")
    parser.add_argument("--workers", type=int, default=1,
                        help="Concurrent isolated body processes; default 1 fits in the current process")
    parser.add_argument("--single", action="store_true",
                        help="Fit sequentially with a fresh process per body; requires --workers 1")


def resolve_domain(args, engine):
    if (args.start_jd is None) != (args.end_jd is None):
        raise ValueError("--start-jd and --end-jd must be supplied together")
    if args.start_jd is not None:
        if args.start is not None or args.end is not None:
            raise ValueError("JD and calendar-year ranges are mutually exclusive")
        return args.start_jd, args.end_jd
    if (args.start is None) != (args.end is None):
        raise ValueError("--start and --end must be supplied together")
    if args.start is None and args.tier == "extended":
        return -3100015.5, 8000016.5
    years = TIER_YEARS.get(args.tier) if args.start is None else (args.start, args.end)
    return tuple(engine.julday(year, 1, 1, 0.0, 1) for year in years)


def resolve_body_tokens(tokens, engine):
    """Resolve public identifiers without inventing a source mapping."""
    result = []
    for token in tokens:
        try:
            result.append(int(token))
        except ValueError:
            value = getattr(engine, token.upper(), None)
            if not isinstance(value, int):
                raise ValueError(f"unknown body: {token}") from None  # noqa: TRY004 -- invalid CLI name, not a Python operand type
            result.append(value)
    return result


def resolve_bodies(args, engine):
    if args.bodies:
        if args.group is not None:
            raise ValueError("--bodies and --group are mutually exclusive")
        return resolve_body_tokens(args.bodies, engine)
    group = args.group or "core"
    if group in ("all", "exotics") and args.tier == "extended":
        from libephemeris.exotic_bodies import EXOTIC_EXTENDED_IDS, EXOTIC_IDS
        return sorted((BODY_GROUPS[group] - set(EXOTIC_IDS)) | set(EXOTIC_EXTENDED_IDS))
    return sorted(BODY_GROUPS[group])


def validate_source_options(args):
    if args.iers_delta_t != (args.iers_cache_dir is not None):
        raise ValueError("--iers-delta-t and --iers-cache-dir must be supplied together")
    if args.spk_dir and not args.auto_sources:
        raise ValueError("--spk-dir requires --auto-sources")
    if args.nbody and not args.auto_sources and not (args.assist_planets and args.assist_asteroids):
        raise ValueError("--nbody requires both --assist-planets and --assist-asteroids")
    if not args.nbody and not args.auto_sources and (args.assist_planets or args.assist_asteroids):
        raise ValueError("ASSIST force files require --nbody")
    if args.verify_samples < 2:
        raise ValueError("verification needs at least two samples per body")


def create_source(args, engine, body_ids=None, domains=None, *, available_domains=False):
    """Retain exactly the requested main, minor and force images."""
    validate_source_options(args)
    if args.kernel is None and args.data_dir is None:
        raise ValueError("supply --kernel or --data-dir; generation uses existing explicit local data")
    kernel = args.kernel or args.data_dir / TIER_KERNELS[args.tier]
    source = NativeSource(kernel, args.tier, args.delta_t_model,
                          delta_t_userdef_days=args.delta_t_userdef_days,
                          iers_cache_dir=args.iers_cache_dir)
    configure_spk_sources(source, args.spk)
    if args.auto_sources:
        if body_ids is None or domains is None or set(body_ids) != set(domains):
            raise ValueError("automatic source assembly requires every actual requested body domain")
        from .discovery import assemble_local_sources
        assemble_local_sources(source, args, body_ids, domains, available_domains=available_domains)
    elif args.nbody:
        nbody_ids = resolve_body_tokens(args.nbody, engine)
        if body_ids is not None and not set(nbody_ids) <= set(body_ids):
            raise ValueError("--nbody IDs must be present in the requested output inventory")
        source.add_nbody_sources(nbody_ids, args.assist_planets, args.assist_asteroids)
    return source


def describe_source(source):
    print(f"Source: {source.kernel}; SHA-256 {source.sha256}")
    print(f"Delta-T: model={source.delta_t_model}; user override days={source.delta_t_userdef_days}; "
          f"local IERS directory={source.iers_cache_dir}")
    for path, digest in source.clock_hashes.items():
        print(f"Clock input: {path}; SHA-256 {digest if digest is not None else 'absent (native fallback)'}")
    for path in sorted(source.ignored_non_spk_files):
        print(f"Ignored non-SPK cache file: {path}")
    for body, data in source.spk_inputs.items():
        print(f"Minor-body source {body}: {data}")
    for body, data in source.nbody_inputs.items():
        print(f"N-body integration source {body}: {data}")


def describe_domains(path):
    """Show the actual published body indices, including their closed endpoints."""
    from .read import read_body_index

    for entry in read_body_index(path)[1]:
        print(f"Body {entry.body_id}: stored TT JD {entry.jd_start:.17g} to {entry.jd_end:.17g}")
