# SPDX-License-Identifier: AGPL-3.0
"""Coordinate complete offline LEB1 generation from explicit local sources."""

import argparse
from pathlib import Path

from .assemble import merge_leb1, replace_leb1
from .generate import fit_grid, generate_leb1, generation_header_start
from .options import (
    BODY_GROUPS,
    TIER_KERNELS,
    create_source,
    describe_domains,
    describe_source,
    fitting_arguments,
    require_native_engine,
    resolve_bodies,
    resolve_domain,
    source_arguments,
    validate_source_options,
)
from .read import read_body_index, read_leb1

PART_GROUPS = ("planets", "asteroids", "exotics", "analytical")


def _inventory(args, engine, mode):
    selection = argparse.Namespace(**vars(args))
    selection.group = "all" if mode in ("full", "single", "groups", "merge") else mode
    selection.bodies = None
    if mode == "body":
        selection.group = None
        selection.bodies = [token.strip() for token in args.body_list.split(",")]
        if any(not token for token in selection.bodies):
            raise ValueError("body lists must contain nonempty names or IDs")
    bodies = resolve_bodies(selection, engine)
    if not bodies or len(set(bodies)) != len(bodies):
        raise ValueError("generation requires a nonempty unique body inventory")
    if not set(bodies) <= BODY_GROUPS["all"]:
        raise ValueError("body list contains an unsupported LEB channel")
    return bodies


def _destinations(args):
    if args.update_main is not None:
        if args.tier == "all" or args.mode not in PART_GROUPS:
            raise ValueError("--update-main requires one tier and one named group")
        if not args.force:
            raise ValueError("--update-main requires --force to replace the existing main")
        read_leb1(args.update_main)  # Reject a missing or invalid main before fitting.
    if args.tier == "all":
        if args.mode is not None or args.kernel is not None or args.output is not None:
            raise ValueError("all takes no mode, --kernel or --output; use --data-dir and optionally --output-dir")
        if args.data_dir is None:
            raise ValueError("all requires --data-dir with each tier's local kernel")
        tiers = tuple(TIER_KERNELS)
        args.mode = "full"
    else:
        if args.mode is None:
            raise ValueError("supply a mode after the tier")
        tiers = (args.tier,)
    root = args.output_dir or (args.data_dir / "leb" if args.data_dir is not None else None)
    if root is None and args.output is None:
        raise ValueError("supply --output, --output-dir or --data-dir")
    plans = []
    for tier in tiers:
        suffix = f"_{args.mode}" if args.mode in BODY_GROUPS and args.mode != "all" else ""
        output = args.output or root / f"ephemeris_{tier}{suffix}.leb"
        parts = [output.parent / f"ephemeris_{tier}_{group}.leb" for group in PART_GROUPS]
        plans.append((tier, output, parts))
    # Check every destination before creating a source or publishing a part.
    outputs = [path for _, output, parts in plans
               for path in ([*parts, output] if args.mode == "groups" else [output])]
    inputs = [Path(item[1]) for item in args.spk]
    inputs.extend(path for path in (args.kernel, args.assist_planets, args.assist_asteroids) if path is not None)
    if args.update_main is not None:
        inputs.append(args.update_main)
    if args.data_dir is not None:
        inputs.extend(args.data_dir / TIER_KERNELS[tier] for tier in tiers)
    resolved = [path.resolve() for path in outputs]
    if len(set(resolved)) != len(resolved) or set(resolved) & {path.resolve() for path in inputs}:
        raise ValueError("output paths must be distinct and must not replace source or force files")
    if args.update_main is not None and any(path.exists() and path.samefile(args.update_main) for path in outputs):
        raise ValueError("group output must differ from the main file, including hard links")
    if args.mode == "merge" and any(output.resolve() in {part.resolve() for part in parts}
                                     for _, output, parts in plans):
        raise ValueError("merged output must differ from every partial input")
    for path in outputs:
        if path.is_dir() or (path.exists() and not args.force):
            raise FileExistsError(path)
    return plans


def _merge_parts(args, engine, parts, output):
    for group, path in zip(PART_GROUPS, parts, strict=True):
        _, entries = read_body_index(path)
        actual = [entry.body_id for entry in entries]
        expected = _inventory(args, engine, group)
        if len(actual) != len(expected) or set(actual) != set(expected):
            raise ValueError(f"{path}: expected complete {group} inventory {expected}, found {actual}")
    # The existing merger checks domains/auxiliary conflicts and verifies all
    # coefficients through the native reader before atomic publication.
    merge_leb1(parts, output, overwrite=args.force)


def _generate_tier(args, engine, output, parts):
    bodies = _inventory(args, engine, args.mode)
    start, end = resolve_domain(args, engine)
    fit_grid(start, end, 1.0)
    source = create_source(args, engine, bodies, dict.fromkeys(bodies, (start, end)),
                           available_domains=args.available_domains)
    destinations = [*parts, output] if args.mode == "groups" else [output]
    if args.update_main is not None:
        destinations.append(args.update_main)
    if any(path.resolve() in source.input_paths for path in destinations):
        raise ValueError("output must differ from every source and force file")
    # Admit the complete selected inventory before the first partial write;
    # fitting still performs its original per-body admission and checks.
    domains = {body: source.available_domain(body, start, end) if args.available_domains else (start, end)
               for body in bodies}
    for body, domain in domains.items():
        source.require_domain(body, *domain)
    header = {}
    if args.mode == "groups" or args.mode in PART_GROUPS:
        # Separately generated named parts must also merge at arbitrary starts.
        # Use the same analytic grid metadata without fitting extra channels.
        header_domains = {body: (start, end) for body in source.piecewise_grids}
        header_domains.update(domains)
        header = {"header_start": generation_header_start(source, header_domains, start)}
    if args.update_main is not None:
        main_header, _, _ = read_leb1(args.update_main)
        if (main_header.jd_start, main_header.jd_end) != (header["header_start"], end):
            raise ValueError("updated group and main must have the same global date range")
    jobs = [(group, path, _inventory(args, engine, group))
            for group, path in zip(PART_GROUPS, parts, strict=True)] if args.mode == "groups" else [
                (args.mode, output, bodies)]
    for group, path, selected in jobs:
        checks = generate_leb1(
            source, path, selected, start, end,
            skip_aux=args.skip_aux or (args.mode == "groups" and group != "planets"),
            overwrite=args.force, verify_samples=args.verify_samples,
            workers=args.workers, single=args.mode == "single",
            available_domains=args.available_domains,
            progress=None if args.quiet else lambda body: print(f"Fitted and checked body {body}", flush=True),
            **header,
        )
        if not args.quiet:
            print(f"Native-readable LEB1: {path}; off-grid position checks={checks.positions}, "
                  f"Cartesian rate checks={checks.cartesian_rates}")
            if args.available_domains:
                describe_domains(path)
    if args.mode == "groups":
        _merge_parts(args, engine, parts, output)
    if args.update_main is not None:
        replace_leb1(args.update_main, [output], args.update_main,
                     expected_bodies=set(bodies), overwrite=True)
        if not args.quiet:
            print(f"Native-readable updated LEB1: {args.update_main}; replaced bodies={bodies}; main auxiliaries retained")
    if not args.quiet:
        describe_source(source)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tier", type=str.lower, choices=(*TIER_KERNELS, "all"))
    parser.add_argument("mode", nargs="?", type=str.lower,
                        choices=("full", "single", "groups", "merge", "body", *BODY_GROUPS),
                        help="Complete generation, a named group, merge, or body NAME[,NAME...]")
    parser.add_argument("body_list", nargs="?", help="Comma-separated names/IDs for body mode")
    source_arguments(parser, include_tier=False)
    fitting_arguments(parser)
    outputs = parser.add_mutually_exclusive_group()
    outputs.add_argument("-o", "--output", type=Path)
    outputs.add_argument("--output-dir", type=Path)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--update-main", type=Path,
                        help="After fitting one named group, replace its bodies in this existing LEB1; requires --force and preserves main auxiliaries")
    parser.add_argument("-q", "--quiet", action="store_true")
    args = parser.parse_args(argv)
    import libephemeris as engine
    try:
        require_native_engine(engine)
        if (args.mode == "body") != (args.body_list is not None):
            raise ValueError("body mode requires one comma-separated body list; other modes take no body list")
        if args.group is not None or args.bodies is not None or args.single:
            raise ValueError("positional mode selects the inventory; --group, --bodies and --single cannot override it")
        validate_source_options(args)
        if args.workers < 1 or (args.mode == "single" and args.workers != 1):
            raise ValueError("workers must be positive; single mode requires exactly one worker")
        if args.mode == "merge" and (args.workers != 1 or args.available_domains or args.auto_sources or args.nbody):
            raise ValueError("merge does not accept fitting workers, available domains or source assembly")
        if args.mode == "merge" and (args.kernel or args.spk or args.spk_dir or args.assist_planets
                                    or args.assist_asteroids or args.skip_aux or args.start is not None
                                    or args.end is not None or args.start_jd is not None or args.end_jd is not None
                                    or args.delta_t_model != "smh2016" or args.delta_t_userdef_days is not None
                                    or args.iers_delta_t or args.iers_cache_dir or args.verify_samples != 64):
            raise ValueError("merge preserves input coefficients/auxiliaries and does not accept source or fitting options")
        plans = _destinations(args)
        for tier, output, parts in plans:
            selected = argparse.Namespace(**vars(args))
            selected.tier = tier
            if selected.mode == "merge":
                _merge_parts(selected, engine, parts, output)
            else:
                _generate_tier(selected, engine, output, parts)
            if not args.quiet and selected.mode in ("groups", "merge"):
                print(f"Native-readable merged LEB1: {output}; partial files retained")
    except (OSError, ValueError, engine.Error) as error:
        parser.exit(1, f"LEB generation failed: {error}\n")
    return 0
