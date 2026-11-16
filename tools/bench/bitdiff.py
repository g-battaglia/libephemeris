"""Capture and compare installed Python APIs with exact binary64 encoding."""

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import gzip
import itertools
import json
import math
import os
from pathlib import Path
import random
import struct
import subprocess
import sys
import tempfile
import warnings

from common import (
    J2000, KERNELS, MODES, SCHEMA, SEED, SITE, TIERS,
    PrerequisiteError, configure, dump, environment, identity, prepare,
)

# Public request bits, including flags that may be rejected for particular bodies.
FLAGS = (
    0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096,
    8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2228280,
)
# Orthodox planets, nodes, apsides, Earth, major asteroids, fictitious,
# numbered asteroid, planetary moon, ecliptic/nutation and invalid ID.
BODIES = (*range(23), 40, 55, 10001, 100433, 9599, -1, -9999)
POSITION = ("calc", "calc_ut", "calc_pctr", "context.calc", "context.calc_ut", "context.calc_pctr")
SITE_ARGS = [SITE[0], SITE[1], SITE[2]]
SPEED = 258
FULL_PER_CONFIG = 350_000
SMOKE_COUNT = 10_000
AUTO_EDGE_COUNT = 30_000
PROFILES = {
    "smoke": (("leb",), ("base",), SMOKE_COUNT),
    "full": (MODES, TIERS, FULL_PER_CONFIG),
    "auto-edge": (("auto",), TIERS, AUTO_EDGE_COUNT),
    "auto-centers": (("auto",), TIERS, SMOKE_COUNT),
}


def encode(value):
    if isinstance(value, float):
        return {"f64": struct.pack(">d", value).hex()}
    if isinstance(value, (tuple, list)):
        return {type(value).__name__: [encode(item) for item in value]}
    if isinstance(value, dict):
        return {"dict": [[encode(key), encode(item)] for key, item in value.items()]}
    if value is None or isinstance(value, (bool, int, str)):
        return value
    raise TypeError(f"unsupported public result type: {type(value)}")


def position(api, jd, body, flags):
    return api, (jd, body, 14, flags) if "pctr" in api else (jd, body, flags)


def requests(mode, tier, inventory, count, seed):
    rng = random.Random(seed)
    core = inventory.get(f"leb/{tier}_core.leb2")
    if mode == "skyfield":
        segments = inventory[KERNELS[tier]]["coverage"]["segments_tdb"]
        lower = max(row[2] for row in segments)
        upper = min(row[3] for row in segments)
    else:
        lower, upper = core["coverage"]["range"]

    site = SITE[1], SITE[0]
    jd = J2000 + 123.5
    yield "calc_ut", (jd, 4, SPEED)
    yield "context.calc_ut", (jd, 4, SPEED)
    for api in ("houses", "houses_ex", "houses_ex2", "houses_with_fallback"):
        for system in ("P", "W", "G"):
            yield api, (jd, *site, system, "O") if api.endswith("fallback") else (jd, *site, system)
    for api in ("houses_armc", "houses_armc_ex2", "houses_armc_with_fallback"):
        args = (120., SITE[1], 23.4, "P")
        yield api, (*args, "O") if api.endswith("fallback") else args
    for api in ("sidtime", "get_ayanamsa", "get_ayanamsa_ut"):
        yield api, (jd,)
    yield "sidtime0", (jd, 23.4, 0.001)
    for api in ("get_ayanamsa_ex", "get_ayanamsa_ex_ut"):
        yield api, (jd, 2)
    for api in ("fixstar", "fixstar_ut", "fixstar2", "fixstar2_ut"):
        for star in ("Spica", "Sirius", "invalid star"):
            yield api, (star, jd, SPEED)
    for api in ("nod_aps", "nod_aps_ut"):
        for method in (1, 2, 4):
            yield api, (jd, 4, method, SPEED)
    for api in ("get_orbital_elements", "get_orbital_elements_ut", "pheno", "pheno_ut"):
        yield api, (jd, 4, 2)
    for event in (1, 2, 4, 8):
        yield "rise_trans", (jd, 0, event, SITE_ARGS)
        yield "rise_trans_true_hor", (jd, 0, event, SITE_ARGS)
    for api in ("solcross", "solcross_ut", "mooncross", "mooncross_ut", "helio_cross", "helio_cross_ut"):
        yield api, (4, 180., jd) if api.startswith("helio") else (180., jd)
    yield "cross_ut", (4, 180., jd)
    yield "mooncross_node", (jd,)
    yield "mooncross_node_ut", (jd,)
    for epoch in (2457987.25, 2460409.25):
        for api in ("sol_eclipse_where", "sol_eclipse_when_glob", "lun_eclipse_when"):
            yield api, (epoch,)
        for api in ("sol_eclipse_how", "sol_eclipse_when_loc", "lun_eclipse_how", "lun_eclipse_when_loc"):
            yield api, (epoch, SITE_ARGS)

    for api in POSITION:
        for body in BODIES:
            for flags in FLAGS:
                yield position(api, jd, body, flags)
        for bad_date in (math.nan, math.inf, -math.inf, lower - 1, upper + 1):
            yield position(api, bad_date, 4, SPEED)
    for api in ("set_topo", "set_sid_mode", "set_delta_t_userdef"):
        new, old = {
            "set_topo": ((-74., 40., 10.), SITE),
            "set_sid_mode": ((1,), (0,)),
            "set_delta_t_userdef": ((0.0007,), (None,)),
        }[api]
        yield api, new
        yield "calc_ut", (jd, 4, SPEED | 32768 | 65536)
        yield api, old
        yield "calc_ut", (jd, 4, SPEED | 32768 | 65536)

    if mode != "skyfield":
        for entry in core["coverage"]["bodies"]:
            boundaries = entry["boundaries"]
            if count == SMOKE_COUNT:
                boundaries = boundaries[::max(1, len(boundaries) // 5)]
            for boundary in boundaries:
                for date in (math.nextafter(boundary, -math.inf), boundary,
                             math.nextafter(boundary, math.inf)):
                    yield "calc", (date, entry["body"], SPEED | 16 | 16384)
    else:
        for target_id, _, left, right in segments:
            body = {10: 0, 301: 1, 199: 2, 299: 3, 4: 4, 5: 5, 6: 6,
                    7: 7, 8: 8, 9: 9, 399: 14}.get(target_id)
            if body is None:
                continue
            for boundary in (left, right):
                for date in (math.nextafter(boundary, -math.inf), boundary,
                             math.nextafter(boundary, math.inf)):
                    yield "calc", (date, body, SPEED | 16 | 16384)

    # Run each combination of the independent geometry/correction bits once.
    # The last block provides dense date coverage on common chart routes.
    bits = (8, 16, 32, 64, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072)
    for mask in range(1 << len(bits)):
        api = POSITION[mask % len(POSITION)]
        body = mask // len(POSITION) % 15
        flags = 2 | sum(bit for i, bit in enumerate(bits) if mask & (1 << i))
        yield position(api, rng.uniform(lower + 3, upper - 3), body, flags)
    for i in itertools.count():
        body = i % 10
        date = lower + 3 + (upper - lower - 6) * ((i // 10 % 40000) / 40000)
        flags = SPEED if i % 19 == 0 else 2
        yield position(("calc_ut", "calc")[i // 10 % 2], date, body, flags)


def auto_center_requests(count, seed):
    """Keep the scalar session active throughout the planet-center experiment."""
    rng = random.Random(seed)
    flags = (2, SPEED, SPEED | 32768, SPEED | 4096, SPEED | 8, SPEED | 16)
    for index in range(count):
        jd = J2000 if index < 2 else J2000 + rng.uniform(-36_525, 36_525)
        api = "calc_ut" if index % 2 == 0 else "calc"
        body = 5 + (index // 2) % 5
        yield position(api, jd, body, flags[(index // 10) % len(flags)])


def outcome(fn, args):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            raw = fn(*args)
        except Exception as exc:
            result = {"error": f"{type(exc).__module__}.{type(exc).__qualname__}"}
        else:
            result = {"value": encode(raw)}
    result["warnings"] = [[f"{w.category.__module__}.{w.category.__qualname__}", str(w.message)]
                          for w in caught]
    return result


def capture(options):
    import libephemeris as eph

    spec = json.loads(options.spec.read_text(encoding="utf-8"))
    build = identity(options.wheel, options.commit, options.build_options)
    build["leb_cache_mib"] = os.environ.get("LIBEPHEMERIS_LEB_CACHE_MB")
    build["positional_centers"] = (
        os.environ.get("LIBEPHEMERIS_EXPERIMENTAL_POSITIONAL_CENTERS") == "1"
    )
    if options.role == "baseline" and build["version"] != "4.0.0":
        raise PrerequisiteError("baseline interpreter must have the published 4.0.0 wheel")
    part = options.output.with_suffix(".partial")
    with gzip.open(part, "wt", encoding="utf-8", newline="\n") as out:
        print(dump({"schema": SCHEMA, "spec": spec["public"], "build": build}), file=out)
        total = 0
        coverage = Counter()
        for mode, tier in spec["matrix"]:
            print(f"{options.role} {mode}/{tier}: starting {spec['count']:,} calls",
                  file=sys.stderr, flush=True)
            configure(eph, mode, tier)
            ctx = eph.EphemerisContext(os.environ["LIBEPHEMERIS_DATA_DIR"], KERNELS[tier])
            ctx.set_topo(*SITE)
            context_apis = {api: getattr(ctx, api.removeprefix("context.")) for api in POSITION
                            if api.startswith("context.")}
            if spec["public"]["profile"] == "auto-centers":
                cases = auto_center_requests(spec["count"], spec["seed"])
            else:
                cases = requests(mode, tier, spec["public"]["data"], spec["count"], spec["seed"])
            for index, (api, args) in enumerate(itertools.islice(cases, spec["count"])):
                fn = context_apis[api] if api.startswith("context.") else getattr(eph, api)
                result = outcome(fn, args)
                if index < 2 and "error" in result:
                    fn(*args)  # Fail loudly if either required preflight call failed.
                if index == 0 and mode != "skyfield" and build["leb_cache_mib"] is not None:
                    usage = eph.get_leb_cache_usage()
                    if usage is None or usage[1] != int(build["leb_cache_mib"]) * 1024**2:
                        raise PrerequisiteError("candidate did not activate the requested LEB process cache")
                row = [mode, tier, index, api, encode(args), result]
                print(dump(row), file=out)
                coverage[f"{mode}/{tier}/{api}"] += 1
                total += 1
                if (index + 1) % 50_000 == 0 or index + 1 == spec["count"]:
                    print(f"{options.role} {mode}/{tier}: {index + 1:,}/{spec['count']:,} calls",
                          file=sys.stderr, flush=True)
        print(dump({"count": total, "coverage": dict(coverage)}), file=out)
    part.replace(options.output)


def compare(left, right):
    with gzip.open(left, "rt", encoding="utf-8") as a, gzip.open(right, "rt", encoding="utf-8") as b:
        header_a, header_b = json.loads(next(a)), json.loads(next(b))
        if header_a["schema"] != SCHEMA or header_b["schema"] != SCHEMA:
            raise ValueError("unsupported capture schema")
        spec = header_a["spec"]
        if spec != header_b["spec"]:
            raise ValueError("captures use different data or corpus configuration")
        compatibility = ("platform", "release", "machine", "libc", "python", "cpu")
        if any(header_a["build"][key] != header_b["build"][key] for key in compatibility):
            raise ValueError("captures must use the same platform and Python build")
        expected = spec["expected_count"]
        matrix = spec["matrix"]
        if not matrix or expected < len(matrix) or expected % len(matrix):
            raise ValueError("invalid capture matrix or count")
        per_config = expected // len(matrix)
        differences = 0
        first = []
        coverage = Counter()
        successes = Counter()
        error_types = Counter()
        for i in range(expected):
            row_a, row_b = json.loads(next(a)), json.loads(next(b))
            if not isinstance(row_a, list) or not isinstance(row_b, list):
                raise ValueError("capture ended before expected count")
            if row_a[:5] != row_b[:5]:
                raise ValueError(f"request identity mismatch: {row_a[:5]} versus {row_b[:5]}")
            if row_a[:3] != [*matrix[i // per_config], i % per_config]:
                raise ValueError("capture does not cover the declared matrix in order")
            mode, tier, _, api, _, result = row_a
            key = f"{mode}/{tier}/{api}"
            coverage[key] += 1
            if "value" in result:
                successes[key] += 1
            else:
                error_types[result["error"]] += 1
            if dump(result) != dump(row_b[-1]):
                differences += 1
                if len(first) < 10:
                    first.append({"request": row_a[:5], "baseline": result, "candidate": row_b[-1]})
        trailer = {"count": expected, "coverage": dict(coverage)}
        if json.loads(next(a)) != trailer or json.loads(next(b)) != trailer or a.read(1) or b.read(1):
            raise ValueError("invalid capture trailer or unexpected records")
        return {"schema": SCHEMA, "status": "different" if differences else "identical",
                "count": expected, "differences": differences, "first_differences": first,
                "coverage": dict(coverage), "baseline_successes": dict(successes),
                "baseline_error_types": dict(error_types),
                "baseline": header_a["build"], "candidate": header_b["build"]}


def run(options):
    options.output.mkdir(parents=True, exist_ok=True)
    report = options.output / "comparison.json"
    report.unlink(missing_ok=True)
    modes, tiers, count = PROFILES[options.profile]
    matrix = list(itertools.product(modes, tiers))
    with tempfile.TemporaryDirectory(prefix="ephemeris-bitdiff-") as tmp:
        data = Path(tmp) / "data"
        inventory = prepare(options.data, data, modes, tiers, companions=options.profile in ("full", "auto-edge"))
        spec_path = Path(tmp) / "spec.json"
        spec_path.write_text(dump({
            "matrix": matrix, "count": count, "seed": SEED,
            "public": {"data": inventory, "profile": options.profile, "matrix": matrix,
                       "seed": SEED, "expected_count": count * len(matrix)},
        }), encoding="utf-8")
        def acquire(role):
            args = [str(getattr(options, role).absolute()), "-E", "-s", str(Path(__file__).resolve()),
                    "capture", "--role", role, "--spec", str(spec_path),
                    "--output", str((options.output / f"{role}.jsonl.gz").absolute())]
            for name in ("wheel", "commit", "build_options"):
                value = getattr(options, f"{role}_{name}")
                if value is not None:
                    args.extend(["--" + name.replace("_", "-"), str(value)])
            env = environment(data)
            if role == "candidate" and options.candidate_cache_mib is not None:
                env["LIBEPHEMERIS_LEB_CACHE_MB"] = str(options.candidate_cache_mib)
            if role == "candidate" and options.candidate_positional_centers:
                env["LIBEPHEMERIS_EXPERIMENTAL_POSITIONAL_CENTERS"] = "1"
            with (options.output / f"{role}.log").open("w", encoding="utf-8", buffering=1) as log:
                with subprocess.Popen(args, cwd=tmp, env=env, stderr=subprocess.PIPE,
                                      encoding="utf-8", errors="replace") as worker:
                    for line in worker.stderr:
                        log.write(line)
                        print(line, end="", file=sys.stderr, flush=True)
                    status = worker.wait()
            if status:
                raise ValueError(f"{role} capture failed ({status}); see {log.name}")

        with ThreadPoolExecutor(max_workers=2) as pool:
            list(pool.map(acquire, ("baseline", "candidate")))
    print("Comparing the completed captures", file=sys.stderr, flush=True)
    result = compare(options.output / "baseline.jsonl.gz", options.output / "candidate.jsonl.gz")
    report.write_text(dump(result) + "\n", encoding="utf-8")
    print(f"{result['status']}: {result['count']} calls, {result['differences']} differences")
    return int(bool(result["differences"]))


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    sub = cli.add_subparsers(dest="action", required=True)
    runner = sub.add_parser("run")
    runner.add_argument("--profile", choices=PROFILES, default="smoke")
    runner.add_argument("--data", type=Path, required=True)
    runner.add_argument("--output", type=Path, required=True)
    runner.add_argument("--candidate-cache-mib", type=int)
    runner.add_argument("--candidate-positional-centers", action="store_true")
    for role in ("baseline", "candidate"):
        runner.add_argument(f"--{role}", type=Path, required=True)
        runner.add_argument(f"--{role}-wheel", type=lambda p: Path(p).resolve())
        runner.add_argument(f"--{role}-commit")
        runner.add_argument(f"--{role}-build-options")
    capturer = sub.add_parser("capture")
    capturer.add_argument("--role", choices=("baseline", "candidate"), required=True)
    capturer.add_argument("--spec", type=Path, required=True)
    capturer.add_argument("--output", type=Path, required=True)
    capturer.add_argument("--wheel", type=Path)
    capturer.add_argument("--commit")
    capturer.add_argument("--build-options")
    comparator = sub.add_parser("compare")
    comparator.add_argument("baseline", type=Path)
    comparator.add_argument("candidate", type=Path)
    options = cli.parse_args()
    try:
        if options.action == "run":
            if options.candidate_cache_mib is not None and options.candidate_cache_mib < 0:
                cli.error("the candidate cache size must be nonnegative")
            return run(options)
        if options.action == "capture":
            capture(options)
            return 0
        result = compare(options.baseline, options.candidate)
        print(dump(result))
        return int(bool(result["differences"]))
    except (PrerequisiteError, OSError, ValueError, StopIteration, subprocess.CalledProcessError) as exc:
        result = {"schema": SCHEMA, "status": "missing_prerequisite" if isinstance(exc, PrerequisiteError)
                  else "incomplete", "error": str(exc)}
        if options.action == "run":
            (options.output / "comparison.json").write_text(dump(result) + "\n", encoding="utf-8")
        print(dump(result), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
