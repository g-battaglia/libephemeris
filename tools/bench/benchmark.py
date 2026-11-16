"""Measure installed wheels and the native Rust example on the same requests."""

from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor
import ctypes
import itertools
import json
import mmap
import os
import platform
import math
from pathlib import Path
import random
import re
import statistics
from datetime import date
import subprocess
import sys
import tempfile
import time
import warnings

from common import (
    J2000, KERNELS, MODES, SCHEMA, SEED, SITE, TIERS, PrerequisiteError,
    command, configure, digest, dump, environment, identity, prepare,
)

CHART_BODIES = [*range(10), 11]
WORKLOADS = ("random_chart", "seq_chart", "same_date_chart", "chart_anatomy", "seq_position",
             "time_series", "houses", "sidtime", "ayanamsha", "fixstar", "topo", "threads",
             "cold_start", "memory")
POSITION_WORKLOADS = ("seq_position", "time_series", "topo")
SERIES = 100_000


def dates(inventory, tier, mode, workload, count, seed, window):
    ranges = []
    if mode != "skyfield":
        bodies = inventory[f"leb/{tier}_core.leb2"]["coverage"]["bodies"]
        ranges.extend((e["start"], e["stop"]) for e in bodies if e["body"] in CHART_BODIES)
    if mode != "leb":
        segments = inventory[KERNELS[tier]]["coverage"]["segments_tdb"]
        ranges.extend((start, stop) for _, _, start, stop in segments)
    lower = max(a for a, _ in ranges) + 3
    upper = min(b for _, b in ranges) - 3
    if window != "coverage":
        first, last = map(int, window.split("-"))
        def new_year_jd(year):
            return J2000 - 0.5 + (date(year, 1, 1) - date(2000, 1, 1)).days
        lower, upper = max(lower, new_year_jd(first)), min(upper, new_year_jd(last))
    if upper <= lower:
        raise PrerequisiteError("no common coverage for the selected dates")
    if workload in ("random_chart", "threads", "memory"):
        rng = random.Random(seed)
        return [rng.uniform(lower, upper) for _ in range(count)], [lower, upper]
    start = max(lower, min(J2000 + 0.5, upper - 0.37 * (count - 1)))
    if workload == "same_date_chart":
        return [start] * count, [lower, upper]
    if start + 0.37 * (count - 1) > upper:
        raise PrerequisiteError("consecutive series exceeds common coverage")
    return [start + i * 0.37 for i in range(count)], [lower, upper]


def operations(spec, eph):
    flags = spec["flags"]
    lat, lon = SITE[1], SITE[0]
    def house(jd):
        return eph.houses(jd, lat, lon, "P")

    def planets(jd, bodies):
        for body in bodies:
            eph.calc_ut(jd, body, flags)

    def chart(jd):
        planets(jd, CHART_BODIES)
        house(jd)

    return chart, {
        "seq_position": lambda jd: eph.calc_ut(jd, 4, flags),
        "time_series": lambda jd: eph.calc_ut(jd, 4, flags),
        "topo": lambda jd: eph.calc_ut(jd, 4, flags | 32768),
        "houses": house,
        "sidtime": eph.sidtime,
        "ayanamsha": eph.get_ayanamsa_ut,
        "fixstar": lambda jd: eph.fixstar_ut("Spica", jd, flags),
    }, planets


def elapsed(fn, *args):
    start = time.perf_counter_ns()
    fn(*args)
    return time.perf_counter_ns() - start


def emit(value):
    print(dump(value), flush=True)


def checkpoint(label):
    peak = None
    if sys.platform in ("darwin", "linux"):
        import resource
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        if sys.platform == "linux":
            peak *= 1024
    emit({"phase": label, "peak_rss_bytes": peak})
    if sys.stdin.readline() != "continue\n":
        raise RuntimeError("memory controller stopped")


def worker(path, *, probe=False):
    spec = json.loads(path.read_text(encoding="utf-8"))
    start = time.perf_counter_ns()
    import libephemeris as eph
    import_ns = time.perf_counter_ns() - start
    jds, workload = spec["dates"], spec["workload"]
    warnings.simplefilter("ignore")  # Timing excludes warning display.
    if workload == "memory" and not probe:
        checkpoint("after_import")
    configure(eph, spec["mode"], spec["tier"])
    chart, other, planets = operations(spec, eph)
    if probe:
        position_probe = workload in POSITION_WORKLOADS
        bodies = (4,) if position_probe else CHART_BODIES
        token = eph.start_tracing()
        try:
            if position_probe:
                other[workload](jds[0])
            else:
                planets(jds[0], bodies)
            routes = eph.get_trace_results()
        finally:
            token.reset()
        expected = "Skyfield" if spec["mode"] == "skyfield" else "LEB"
        if any(routes.get(body) != expected for body in bodies):
            raise PrerequisiteError(f"unexpected provider for {workload}: {routes}")
        emit(routes)
        return
    if workload in ("cold_start", "memory"):
        first = elapsed(eph.calc_ut, jds[0], 4, spec["flags"])
        if workload == "memory":
            checkpoint("after_first_call")
        first_chart = elapsed(chart, jds[0])
        if workload == "memory":
            checkpoint("after_first_chart")
            for jd in jds:
                chart(jd)
            checkpoint("steady_state")
        emit({"import_ns": import_ns, "first_call_ns": first, "first_chart_ns": first_chart})
        return
    op = other.get(workload, chart)
    op(jds[0] - 0.1)  # Warm initialization at a different epoch.
    if workload == "chart_anatomy":
        samples = []
        for jd in jds:
            samples.append([
                elapsed(planets, jd, CHART_BODIES[:1]),
                elapsed(planets, jd, CHART_BODIES[1:]),
                elapsed(other["houses"], jd),
            ])
        emit({"sample_ns": samples, "columns": ["first_body", "other_ten_bodies", "houses"],
              "sample_unit": "chart component"})
    elif workload == "threads":
        count = spec["parallelism"]
        def work(dates):
            for jd in dates:
                chart(jd)

        with ThreadPoolExecutor(count) as pool:
            start = time.perf_counter_ns()
            list(pool.map(work, (jds[i::count] for i in range(count))))
            duration = time.perf_counter_ns() - start
        emit({"sample_ns": [duration / len(jds)], "elapsed_ns": duration,
              "sample_unit": "batch average per chart"})
    else:
        emit({"sample_ns": [elapsed(op, jd) for jd in jds], "sample_unit": "operation"})


def kib_fields(path):
    fields = {}
    for line in path.read_text().splitlines():
        if match := re.fullmatch(r"(\w+):\s+(\d+) kB", line):
            fields[match[1]] = int(match[2]) * 1024
    return fields


def size_bytes(raw):
    if not (match := re.fullmatch(r"([\d.]+)([BKMGTP]?)", raw)):
        return None
    return round(float(match[1]) * 1024 ** ("BKMGTP".index(match[2] or "B")))


def snapshot(pid, raw_dir, label):
    values = dict(rss_bytes=None, anon_bytes=None, file_resident_bytes=None,
                  pss_bytes=None, shared_rss_bytes=None, footprint_bytes=None)
    if sys.platform == "linux":
        status = kib_fields(Path(f"/proc/{pid}/status"))
        totals = kib_fields(Path(f"/proc/{pid}/smaps_rollup"))
        values.update(rss_bytes=totals.get("Rss"), pss_bytes=totals.get("Pss"),
                      anon_bytes=status.get("RssAnon"), file_resident_bytes=status.get("RssFile"))
        if "Shared_Clean" in totals and "Shared_Dirty" in totals:
            values["shared_rss_bytes"] = totals["Shared_Clean"] + totals["Shared_Dirty"]
    elif sys.platform == "darwin":
        ps = subprocess.check_output(["ps", "-o", "rss=", "-p", str(pid)], text=True)
        values["rss_bytes"] = int(ps.strip()) * 1024
        output = subprocess.run(["vmmap", "-summary", str(pid)], text=True, capture_output=True)
        (raw_dir / f"{label}-{pid}.vmmap.txt").write_text(output.stdout + output.stderr, encoding="utf-8")
        for name, expression in (
            ("footprint_bytes", r"^Physical footprint:\s+(\S+)"),
            ("peak_footprint_bytes", r"^Physical footprint \(peak\):\s+(\S+)"),
            ("readonly_library_resident_bytes", r"^ReadOnly portion of Libraries:.*resident=(\S+?)\("),
        ):
            if match := re.search(expression, output.stdout, re.MULTILINE):
                values[name] = size_bytes(match[1])
        values["malloc_resident_bytes"] = sum(
            size_bytes(match) or 0
            for match in re.findall(r"^MALLOC[\w ()]+\s+[\d.]+[KMGB]\s+([\d.]+[KMGB])\s", output.stdout, re.MULTILINE)
        ) if output.returncode == 0 else None
    return values


def memory_run(args, env, cwd, workers, raw_dir):
    children = []
    try:
        for _ in range(workers):
            children.append(subprocess.Popen(args, cwd=cwd, env=env, text=True,
                                             stdin=subprocess.PIPE, stdout=subprocess.PIPE))
        phases = []
        for phase in ("after_import", "after_first_call", "after_first_chart", "steady_state"):
            readings = []
            for child in children:
                record = json.loads(child.stdout.readline())
                if record["phase"] != phase:
                    raise RuntimeError("memory worker did not reach expected checkpoint")
                record.update(snapshot(child.pid, raw_dir, phase))
                readings.append(record)

            def aggregate(key):
                values = [r[key] for r in readings]
                return sum(values) if all(v is not None for v in values) else None

            phases.append({"phase": phase, "processes": readings, "sum_rss_bytes": aggregate("rss_bytes"),
                           "sum_pss_bytes": aggregate("pss_bytes")})
            for child in children:
                child.stdin.write("continue\n")
                child.stdin.flush()
        for child in children:
            stdout, _ = child.communicate()
            if child.returncode != 0:
                raise subprocess.CalledProcessError(child.returncode, args)
            json.loads(stdout)
        return {"checkpoints": phases}
    finally:
        for child in children:
            if child.poll() is None:
                child.kill()
            child.communicate()


def summarize(samples):
    if not samples:
        return None
    if "sample_ns" in samples[0]:
        raw = list(itertools.chain.from_iterable(row["sample_ns"] for row in samples))
        names = samples[0].get("columns")
        columns = zip(*raw) if names else [raw]
        unit = samples[0]["sample_unit"]
    elif "first_call_ns" in samples[0]:
        names = ("import", "first_call", "first_chart")
        columns = ([sample[f"{name}_ns"] for sample in samples] for name in names)
        unit = "fresh process"
    else:
        return None
    result = []
    for values in columns:
        values = sorted(values)
        result.append({"median_ns": statistics.median(values),
                       "p95_ns": values[math.ceil(len(values) * 0.95) - 1]})
    summary = {"sample_unit": unit, "statistics": result}
    if names:
        summary["columns"] = names
    return summary


def warm_os_pages(data, names):
    """Prime the staged inputs immediately before a timed fresh process."""
    buffer = bytearray(1024 * 1024)
    for name in names:
        with (data / name).open("rb", buffering=0) as source:
            while source.readinto(buffer):
                pass


def linux_resident_bytes(fd, mincore):
    length = os.fstat(fd).st_size
    pages = (length + mmap.PAGESIZE - 1) // mmap.PAGESIZE
    residency = (ctypes.c_ubyte * pages)()
    with mmap.mmap(fd, length, access=mmap.ACCESS_COPY) as view:
        address = ctypes.addressof(ctypes.c_ubyte.from_buffer(view))
        if mincore(address, length, residency) != 0:
            error = ctypes.get_errno()
            raise OSError(error, os.strerror(error))
    return sum(value & 1 for value in residency) * mmap.PAGESIZE


def verify_cold_os_pages(data, names):
    """Evict and verify the private input copies immediately before startup."""
    resident = 0
    try:
        mincore = ctypes.CDLL(None, use_errno=True).mincore
        mincore.argtypes = (ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_ubyte))
        mincore.restype = ctypes.c_int
        for name in names:
            with (data / name).open("r+b") as source:
                os.fsync(source.fileno())
                os.posix_fadvise(source.fileno(), 0, 0, os.POSIX_FADV_DONTNEED)
                resident += linux_resident_bytes(source.fileno(), mincore)
    except (AttributeError, OSError, ValueError) as exc:
        raise PrerequisiteError(f"cannot check staged input page residency: {exc}") from exc
    if resident:
        raise PrerequisiteError(f"mincore reported {resident} resident bytes after eviction")
    return resident


def run(options):
    report = {"schema": SCHEMA, "status": "incomplete", "seed": options.seed,
              "runs": [], "site": SITE, "network_policy": "sealed",
              "delta_t_model": "smh2016", "iers_enabled": False}
    options.output.parent.mkdir(parents=True, exist_ok=True)
    try:
        with tempfile.TemporaryDirectory(prefix="ephemeris-bench-") as tmp:
            scratch = Path(tmp)
            data = scratch / "data"
            modes, tiers = options.mode or MODES, options.tier or TIERS
            if options.cold_condition == "page_cache_cold" and sys.platform != "linux":
                raise PrerequisiteError("verified input page eviction requires Linux")
            report["data"] = prepare(options.data, data, modes, tiers, companions=options.companions)
            if options.native:
                report["build"] = {"binary": str(options.native), "sha256": digest(options.native),
                                   "platform": platform.platform(), "rustc_on_host": command("rustc", "-Vv"),
                                   "commit": options.commit, "build_options": options.build_options}
            else:
                report["build"] = identity(options.wheel, options.commit, options.build_options)
            work = options.workload or ("random_chart", "seq_chart", "chart_anatomy")
            raw_dir = options.output.parent / (options.output.stem + "-vmmap")
            if "memory" in work and sys.platform == "darwin":
                raw_dir.mkdir(exist_ok=True)
            routes = {}
            for mode, tier, workload in itertools.product(modes, tiers, work):
                count = options.calls
                if workload == "memory":
                    count = options.memory_calls
                elif workload == "time_series":
                    count = options.series_calls
                elif workload == "cold_start":
                    count = 1
                jds, span = dates(report["data"], tier, mode, workload, count, options.seed, options.window)
                env = environment(data)
                position_probe = workload in POSITION_WORKLOADS
                route_key = (mode, tier, workload if position_probe else "chart")
                for parallel in ((options.workers or (1, 2, 4, 8)) if workload in ("threads", "memory") else (1,)):
                    spec = {"schema": SCHEMA, "mode": mode, "tier": tier, "workload": workload,
                            "flags": options.flags, "dates": jds, "parallelism": parallel}
                    job = scratch / "job.json"
                    job.write_text(dump(spec), encoding="utf-8")
                    args = ([str(options.native), str(job)] if options.native else
                            [sys.executable, "-E", "-s", str(Path(__file__).resolve()), "--worker", str(job)])
                    if options.native and workload in ("threads", "cold_start", "memory"):
                        raise PrerequisiteError(f"native example does not implement {workload}")
                    if not options.native and route_key not in routes:
                        routes[route_key] = json.loads(
                            subprocess.check_output([*args, "--probe"], env=env, cwd=tmp, text=True)
                        )
                    samples = []
                    if workload == "memory":
                        samples.append(memory_run(args, env, tmp, parallel, raw_dir))
                    else:
                        for _ in range(options.repeats):
                            resident_bytes = None
                            if workload == "cold_start" and options.cold_condition == "page_cache_warm":
                                warm_os_pages(data, report["data"])
                            if workload == "cold_start" and options.cold_condition == "page_cache_cold":
                                resident_bytes = verify_cold_os_pages(data, report["data"])
                            output = subprocess.check_output(args, env=env, cwd=tmp, text=True)
                            sample = json.loads(output)
                            if resident_bytes is not None:
                                sample["input_resident_bytes_before_spawn"] = resident_bytes
                            samples.append(sample)
                    report["runs"].append({
                        "mode": mode, "tier": tier, "workload": workload, "flags": options.flags,
                        "parallelism": parallel, "requests_per_run": count, "date_range": span,
                        "samples": samples, "summary": summarize(samples),
                        "chart_probe_routes": routes.get(route_key) if not position_probe else None,
                        "position_probe_routes": routes.get(route_key) if position_probe else None,
                        "cold_condition": options.cold_condition if workload == "cold_start" else None,
                        "startup_steps_per_process": {"calc_ut": 1, "charts": 1}
                        if workload in ("cold_start", "memory") else None,
                        "input_pages_primed": options.cold_condition == "page_cache_warm"
                        if workload == "cold_start" else None,
                    })
            report["status"] = "complete"
    except subprocess.CalledProcessError as exc:
        report.update(status="missing_prerequisite" if exc.returncode == 3 else "incomplete",
                      error=str(exc))
    except (PrerequisiteError, OSError, RuntimeError, ValueError) as exc:
        report.update(status="missing_prerequisite" if isinstance(exc, PrerequisiteError) else "incomplete",
                      error=str(exc))
    options.output.write_text(dump(report) + "\n", encoding="utf-8")
    if report["status"] != "complete":
        print(report.get("error"), file=sys.stderr)
        return 2
    print(f"{report['status']}: {len(report['runs'])} cases, {options.output}")
    return 0


def main():
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("--worker", type=Path)
    cli.add_argument("--probe", action="store_true")
    cli.add_argument("--data", type=Path)
    cli.add_argument("--output", type=Path)
    cli.add_argument("--wheel", type=Path)
    cli.add_argument("--commit")
    cli.add_argument("--build-options")
    cli.add_argument("--native", type=lambda p: Path(p).resolve())
    cli.add_argument("--mode", choices=MODES, action="append")
    cli.add_argument("--tier", choices=TIERS, action="append")
    cli.add_argument("--workload", choices=WORKLOADS, action="append")
    cli.add_argument("--window", choices=("1900-2100", "1550-2650", "coverage"), default="1900-2100")
    cli.add_argument("--cold-condition", choices=("first_use", "page_cache_warm", "page_cache_cold"),
                     default="page_cache_warm")
    cli.add_argument("--calls", type=int, default=1000)
    cli.add_argument("--memory-calls", type=int, default=SERIES)
    cli.add_argument("--series-calls", type=int, default=SERIES)
    cli.add_argument("--workers", type=int, choices=(1, 2, 4, 8), action="append")
    cli.add_argument("--companions", action="store_true", help="include the released companion groups")
    cli.add_argument("--repeats", type=int, default=5)
    cli.add_argument("--flags", type=int, default=258)
    cli.add_argument("--seed", type=int, default=SEED)
    options = cli.parse_args()
    if options.worker:
        try:
            worker(options.worker, probe=options.probe)
            return 0
        except PrerequisiteError as exc:
            print(str(exc), file=sys.stderr)
            return 3
    if options.data is None or options.output is None or min(options.calls, options.memory_calls, options.series_calls) < 1 or options.repeats < 5:
        cli.error("provide data, output, positive call counts and at least five repetitions")
    return run(options)


if __name__ == "__main__":
    raise SystemExit(main())
