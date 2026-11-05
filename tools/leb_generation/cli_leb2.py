#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0
"""Standalone v4 LEB2 generation, conversion and native-reader verification."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import format as wire
from .convert import convert_all_leb1, convert_leb1
from .generate import generate_leb2
from .groups import LEB2_GROUPS
from .options import (
    create_source,
    describe_domains,
    describe_source,
    fitting_arguments,
    require_native_engine,
    resolve_bodies,
    resolve_domain,
    source_arguments,
)
from .read import read_body_index
from .verify import verify_source_file
from .writer import verify_native_index


def compression_options(convert):
    """The same explicit storage choices apply to one file and a named split."""
    convert.add_argument("--precision", type=float, help="With --lossy: maximum coefficient-sum error; cannot relax a body's default")
    representation = convert.add_mutually_exclusive_group()
    representation.add_argument("--lossless", dest="lossless", action="store_true", help="Preserve every coefficient bit (default)")
    representation.add_argument("--lossy", dest="lossless", action="store_false", help="Truncate coefficient bits; component budgets alone do not guarantee astronomical rate precision")
    convert.set_defaults(lossless=True)
    convert.add_argument("--chunk-years", type=float, default=10.0)
    convert.add_argument("--force", action="store_true", help="Atomically replace an existing destination after verification")
    convert.add_argument("-q", "--quiet", action="store_true")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    convert = commands.add_parser("convert", help="Convert one LEB1 file into a verified LEB2 v2 file")
    convert.add_argument("input", type=Path)
    convert.add_argument("-o", "--output", type=Path, required=True)
    convert.add_argument("--group", choices=LEB2_GROUPS)
    convert.add_argument("--tier", choices=("base", "medium", "extended"))
    compression_options(convert)
    split = commands.add_parser("convert-all", help="Split all present complete groups into standard tier filenames")
    split.add_argument("input", type=Path)
    split.add_argument("-o", "--output", type=Path, required=True, help="Output directory")
    split.add_argument("--tier-name", choices=("base", "medium", "extended"), required=True)
    compression_options(split)
    generate = commands.add_parser("generate", help="Fit the selected native sources and publish a verified LEB2 file")
    generate.add_argument("-o", "--output", type=Path, required=True)
    source_arguments(generate)
    fitting_arguments(generate)
    compression_options(generate)
    verify = commands.add_parser("verify", help="Verify an existing LEB2 file against explicit native sources")
    verify.add_argument("input", type=Path)
    source_arguments(verify)
    verify.add_argument("--structural-only", action="store_true",
                        help="Check the selected file's native index and all chunks; makes no source-precision claim")
    verify.add_argument("-q", "--quiet", action="store_true")
    args = parser.parse_args(argv)
    import libephemeris as engine
    try:
        require_native_engine(engine)
        if args.command == "verify":
            header, entries = read_body_index(args.input)
            if header.magic != wire.LEB2_MAGIC:
                raise ValueError("LEB2 verification requires a LEB2 source file")
            if not entries:
                raise ValueError("verification requires a nonempty body inventory")
            if args.structural_only:
                if (args.kernel or args.data_dir or args.spk or args.nbody or args.auto_sources or args.spk_dir or
                        args.assist_planets or args.assist_asteroids or args.delta_t_userdef_days is not None or
                        args.iers_delta_t or args.iers_cache_dir):
                    raise ValueError("--structural-only does not accept astronomical source selections")
                source = None
            else:
                source = create_source(args, engine, [entry.body_id for entry in entries],
                                       {entry.body_id: (entry.jd_start, entry.jd_end) for entry in entries})
            verify_native_index(args.input, entries)
            counts = None if source is None else verify_source_file(args.input, source, args.verify_samples)
            if not args.quiet:
                if source is not None:
                    describe_source(source)
                    print(f"Native-reader source checks: {counts}")
                else:
                    print(f"Native structure and all payload chunks verified: {args.input}; source precision not checked")
            return 0
        else:
            options = {"target": args.precision, "lossless": args.lossless,
                       "chunk_interval_days": args.chunk_years * 365.25, "overwrite": args.force}
            if args.command == "generate":
                body_ids = resolve_bodies(args, engine)
                start, end = resolve_domain(args, engine)
                source = create_source(args, engine, body_ids, dict.fromkeys(body_ids, (start, end)),
                                       available_domains=args.available_domains)
                checks = generate_leb2(source, args.output, body_ids, start, end,
                                       skip_aux=args.skip_aux, verify_samples=args.verify_samples,
                                       workers=args.workers, single=args.single,
                                       available_domains=args.available_domains,
                                       progress=None if args.quiet else lambda body: print(f"Fitted and checked body {body}", flush=True),
                                       **options)
                if not args.quiet:
                    if args.available_domains:
                        describe_domains(args.output)
                    describe_source(source)
                    print(f"Off-grid position checks={checks.positions}; Cartesian rate checks={checks.cartesian_rates}")
                outputs = [args.output]
            elif args.command == "convert-all":
                outputs = convert_all_leb1(args.input, args.output, args.tier_name, **options)
            else:
                convert_leb1(args.input, args.output, group=args.group, tier=args.tier, **options)
                outputs = [args.output]
    except (OSError, ValueError, engine.Error) as error:
        parser.exit(1, f"LEB2 {args.command} failed: {error}\n")
    if not args.quiet:
        for output in outputs:
            print(f"Native-readable LEB2 ({'lossless' if args.lossless else 'lossy'}): {output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
