#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0
"""Generate and merge LEB1 using the installed v4 native engine and local JPL data."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .assemble import merge_leb1
from .generate import generate_leb1
from .options import (
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


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--output", type=Path)
    source_arguments(parser)
    fitting_arguments(parser)
    parser.add_argument("--merge", nargs="+", type=Path, help="Merge existing disjoint LEB1 parts without source astronomy")
    parser.add_argument("--verify-only", type=Path, help="Verify existing LEB1/LEB2 against the exact local native source")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("-q", "--quiet", action="store_true")
    args = parser.parse_args(argv)
    import libephemeris as engine
    try:
        require_native_engine(engine)
        validate_source_options(args)
        if (args.merge or args.verify_only) and (args.single or args.workers != 1):
            raise ValueError("--workers and --single apply only to generation")
        if (args.merge or args.verify_only) and args.available_domains:
            raise ValueError("--available-domains applies only to generation")
        if args.merge and args.auto_sources:
            raise ValueError("--auto-sources does not apply to merge")
        if args.merge and args.nbody:
            raise ValueError("merge does not accept N-body source options")
        if args.merge and (args.delta_t_userdef_days is not None or args.iers_delta_t or args.iers_cache_dir):
            raise ValueError("merge preserves auxiliary data and does not accept clock selections")
        if args.output is None and not args.verify_only:
            raise ValueError("generation and merge require --output")
        if args.merge:
            if args.verify_only:
                raise ValueError("--merge and --verify-only are mutually exclusive")
            merge_leb1(args.merge, args.output, overwrite=args.force)
            message = f"Native-readable merged LEB1: {args.output}"
        else:
            if args.verify_only:
                from .read import read_body_index
                from .verify import verify_source_file
                _, entries = read_body_index(args.verify_only)
                body_ids = [entry.body_id for entry in entries]
                domains = {entry.body_id: (entry.jd_start, entry.jd_end) for entry in entries}
                source = create_source(args, engine, body_ids, domains)
                counts = verify_source_file(args.verify_only, source, args.verify_samples)
                message = f"Native-reader source checks: {counts}"
            else:
                body_ids = resolve_bodies(args, engine)
                start, end = resolve_domain(args, engine)
                domains = dict.fromkeys(body_ids, (start, end))
                source = create_source(args, engine, body_ids, domains,
                                       available_domains=args.available_domains)
                checks = generate_leb1(source, args.output, body_ids, start, end,
                                       skip_aux=args.skip_aux, overwrite=args.force,
                                       workers=args.workers, single=args.single,
                                       available_domains=args.available_domains,
                                       verify_samples=args.verify_samples,
                                       progress=None if args.quiet else lambda body: print(f"Fitted and checked body {body}", flush=True))
                message = (f"Native-readable LEB1: {args.output}; off-grid position checks={checks.positions}, "
                           f"Cartesian rate checks={checks.cartesian_rates}")
                if args.available_domains and not args.quiet:
                    describe_domains(args.output)
            if not args.quiet:
                describe_source(source)
    except (OSError, ValueError, engine.Error) as error:
        parser.exit(1, f"LEB generation failed: {error}\n")
    if not args.quiet:
        print(message)
    return 0


if __name__ == "__main__":
    sys.exit(main())
