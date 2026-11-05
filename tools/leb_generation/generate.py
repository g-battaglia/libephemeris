# SPDX-License-Identifier: AGPL-3.0
"""Bounded offline Chebyshev fitting using native v4 source astronomy."""

import math
import struct
import tempfile
import time
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np
from numpy.polynomial import chebyshev as cheb

from . import format as wire
from .assemble import write_leb1
from .checks import check_positions
from .source import ANGULAR_IDS
from .writer import BodyData


@dataclass
class FitChecks:
    """Counts of off-grid checks; angular rates need separate pipeline validation."""

    positions: int = 0
    cartesian_rates: int = 0


def fit_grid(start: float, end: float, preferred: float) -> tuple[int, float]:
    """Use complete fixed-width segments inside the explicit source domain.

    Shortening the nominal interval retains the complete requested domain and
    never increases the reviewed body's interval. One upward binary64 step
    prevents ceil(span/interval) from inventing an additional segment. Fitting
    nodes stay strictly interior; the final public endpoint is checked too.
    """
    if not all(math.isfinite(value) for value in (start, end, preferred)) or end <= start or preferred <= 0:
        raise ValueError("fit domain and interval must be finite with positive width")
    span = end - start
    if not math.isfinite(span):
        raise ValueError("fit domain span must be finite")
    count = math.ceil(span / preferred)
    interval = span / count
    if math.ceil(span / interval) > count:
        interval = math.nextafter(interval, math.inf)
    return count, interval


def piecewise_fit_grid(start: float, end: float, preferred: float, source_grid):
    """Align a uniform polynomial grid with every native residual/taper knot.

    The common spacing follows exact rational source metadata, never a fitted
    error or a chosen comparison epoch. A body may begin just before the user's
    requested start so an arbitrary start cannot place a knot inside a segment.
    The file header records that actual union; each body's domain stays explicit.
    """
    fit_grid(start, end, preferred)
    origin, table_end, step, taper = source_grid
    if not all(math.isfinite(value) for value in source_grid) or table_end <= origin or min(step, taper) <= 0:
        raise ValueError("invalid piecewise source-grid metadata")
    periods = [Fraction(value) for value in (preferred, step, taper, table_end - origin)]
    denominator = math.lcm(*(period.denominator for period in periods))
    numerator = math.gcd(*(period.numerator * (denominator // period.denominator) for period in periods))
    interval = numerator / denominator
    left = origin + math.floor((start - origin) / interval) * interval
    if left > start:
        left -= interval
    if left + interval == left:
        raise ValueError("source-grid spacing is below Julian-day resolution")
    count = math.ceil((end - left) / interval)
    if count > 2**32 - 1:
        raise ValueError("piecewise source requires more segments than the LEB format can represent")
    return left, count, interval


def generation_header_start(source, domains, requested_start):
    """Common global start for group files with independently aligned bodies.

    Residual/taper knots may move an angular body's first segment earlier than
    the requested start. All parts of one generation must advertise that same
    union in their global header; each body's actual closed domain is retained.
    This uses source-grid metadata only, without sampling or extrapolation.
    """
    starts = [requested_start]
    for body, (low, high) in domains.items():
        if grid := source.piecewise_grids.get(body):
            low = piecewise_fit_grid(low, high, wire.BODY_PARAMS[body][0], grid)[0]
        starts.append(low)
    return min(starts)


def fit_body(source, body_id: int, start: float, end: float, directory: Path,
             checks: FitChecks) -> BodyData:
    """Fit one body in bounded batches and check independent points per segment.

    No source value is clamped, extrapolated or replaced. Coefficients retain
    binary64. Cartesian rates come from analytic native SPK velocities. Angular
    channels are periodic in longitude. Every channel is sampled on the TT
    grid required by the LEB reader; Rust maps that clock into each source.
    """
    preferred, degree, coord_type, components = wire.BODY_PARAMS[body_id]
    count, interval = fit_grid(start, end, preferred)
    if grid := source.piecewise_grids.get(body_id):
        start, count, interval = piecewise_fit_grid(start, end, preferred, grid)
    source.require_domain(body_id, start, end)
    nodes = np.cos(np.pi * (np.arange(degree + 1) + 0.5) / (degree + 1))
    check_nodes = np.linspace(-1.0, 1.0, 9)
    angular = body_id in ANGULAR_IDS
    path = directory / f"body-{body_id}.f64"
    with path.open("wb") as destination:
        for first in range(0, count, 128):
            size = min(128, count - first)
            left = start + (first + np.arange(size)) * interval
            epochs = left[:, None] + (nodes[None, :] + 1.0) * interval / 2.0
            values = np.asarray(source.values(body_id, epochs.ravel()), dtype=np.float64)
            if values.shape != (size * (degree + 1), components) or not np.isfinite(values).all():
                raise ValueError(f"body {body_id}: invalid native source samples")
            values = values.reshape(size, degree + 1, components)
            if angular:
                values[:, :, 0] = np.rad2deg(np.unwrap(np.deg2rad(values[:, :, 0]), axis=1))
            # Fit at the actual representable source epochs, not ideal nodes
            # lost when added to a large JD. Otherwise that rounding becomes
            # polynomial noise and is amplified in differentiated channels.
            actual_nodes = 2.0 * (epochs - left[:, None]) / interval - 1.0
            fitted = np.linalg.solve(cheb.chebvander(actual_nodes, degree), values)
            coefficients = fitted.transpose(0, 2, 1)
            if not np.isfinite(coefficients).all():
                raise ValueError(f"body {body_id}: nonfinite fitted coefficients")
            test_epochs = left[:, None] + (check_nodes[None, :] + 1.0) * interval / 2.0
            # The stored closed end may differ by an interval-rounding ulp
            # from the nominal last segment end; verify the actual public end.
            if first + size == count:
                test_epochs[-1, -1] = end
            if angular:
                expected = source.values(body_id, test_epochs.ravel()).reshape(size, len(check_nodes), components)
            else:
                # One coherent trajectory batch supplies both observables;
                # restarting an integrator for each segment is unnecessary.
                expected_states = source.states(body_id, test_epochs.ravel()).reshape(size, len(check_nodes), 6)
                expected = expected_states[:, :, :3]
            for ordinal in range(size):
                actual_nodes = 2.0 * (test_epochs[ordinal] - left[ordinal]) / interval - 1.0
                actual = cheb.chebval(actual_nodes, coefficients[ordinal].T).T
                check_positions(expected[ordinal], actual, angular, body_id)
                checks.positions += len(check_nodes)
                if not angular:
                    expected_rates = expected_states[ordinal, :, 3:]
                    rates = cheb.chebval(actual_nodes, cheb.chebder(coefficients[ordinal].T, axis=0)).T * (2.0 / interval)
                    difference = np.linalg.norm(rates - expected_rates, axis=1)
                    limit = 1e-6 * np.linalg.norm(expected_rates, axis=1)
                    if np.any(difference > limit):
                        point = int(np.flatnonzero(difference > limit)[0])
                        raise ValueError(
                            f"body {body_id}: fitted Cartesian rate exceeds PLAN relative criterion; "
                            f"segment={first + ordinal}, start_tt={left[ordinal]:.17g}, "
                            f"interval_days={interval:.17g}, jd_tt={test_epochs[ordinal, point]:.17g}, "
                            f"error_au_per_tt_day={difference[point]:.17g}, "
                            f"limit_au_per_tt_day={limit[point]:.17g}, "
                            f"fitted={rates[point].tolist()}, source={expected_rates[point].tolist()}"
                        )
                    checks.cartesian_rates += len(check_nodes)
            destination.write(np.asarray(coefficients, dtype="<f8").tobytes())
    entry = wire.BodyEntry(body_id, coord_type, count, start, end, interval, degree, components, 0)
    return BodyData(entry, path=path)


def auxiliary_sections(start: float, end: float, directory: Path):
    """Native IAU 2006/2000A nutation, model Delta-T and catalogue snapshots.

    Nutation is radians at TT; Delta-T stores days at UT1 every 30 days, with
    both exact endpoints. Catalogue proper motions are deg/year. The LEB star
    snapshot retains the established zero-parallax/radial-velocity convention;
    the runtime's complete native catalogue remains authoritative for stars.
    """
    from libephemeris.fixed_stars import STAR_CATALOG
    from libephemeris.time_utils import deltat

    from .nutation import fit_nutation
    nutation = fit_nutation(start, end, directory)
    epochs = np.arange(start, end, 30.0).tolist() + [end]
    delta_t = bytearray(struct.pack(wire.DELTA_T_HEADER_FMT, len(epochs), 0))
    for epoch in epochs:
        value = deltat(epoch)
        if not math.isfinite(value):
            raise ValueError("nonfinite native Delta-T sample")
        delta_t.extend(struct.pack(wire.DELTA_T_ENTRY_FMT, epoch, value))
    stars = bytearray(len(STAR_CATALOG) * wire.STAR_ENTRY_SIZE)
    for index, star in enumerate(STAR_CATALOG):
        wire.write_star_entry(stars, index * wire.STAR_ENTRY_SIZE, wire.StarEntry(
            star.id, star.data.ra_j2000, star.data.dec_j2000,
            star.data.pm_ra / 3600.0, star.data.pm_dec / 3600.0,
            0.0, 0.0, star.magnitude,
        ))
    return {wire.SECTION_NUTATION: nutation,
            wire.SECTION_DELTA_T: bytes(delta_t), wire.SECTION_STARS: bytes(stars)}


def generate_leb1(source, output: Path, body_ids: list[int], start: float, end: float,
                  *, skip_aux=False, overwrite=False, verify_samples=64, progress=None,
                  workers=1, single=False, available_domains=False, header_start=None) -> FitChecks:
    """Generate from native sources, check fits, then verify and publish LEB1.

    Body coefficients and nutation are backed by operation-owned temporary
    files, with bounded fitting batches. A failed body leaves the destination
    unchanged. No v3 package, Skyfield or network access is required.
    Explicit ``available_domains`` retains each body's complete continuous
    source overlap in its index; the default requires the full requested range.
    ``header_start`` permits the common metadata/auxiliary start of a grouped
    generation; it never changes a body's fitting domain or source selection.
    """
    if not body_ids or len(set(body_ids)) != len(body_ids) or not set(body_ids) <= wire.BODY_PARAMS.keys():
        raise ValueError("generation requires a nonempty unique supported body inventory")
    if verify_samples < 2:
        raise ValueError("verification needs at least two samples per body")
    if workers < 1 or (single and workers != 1):
        raise ValueError("workers must be positive; --single requires exactly one worker")
    output = Path(output)
    if output.exists() and not overwrite:
        raise FileExistsError(output)
    if output.resolve() in source.input_paths:
        raise ValueError("output must differ from every source and force file")
    fit_grid(start, end, 1.0)  # Validate the requested file domain before sampling.
    if header_start is not None and (not math.isfinite(header_start) or header_start > start):
        raise ValueError("shared header start must be finite and no later than the requested start")
    domains = {body: source.available_domain(body, start, end) if available_domains else (start, end)
               for body in body_ids}
    for body, (low, high) in domains.items():
        source.require_domain(body, low, high)
    checks = FitChecks()
    with tempfile.TemporaryDirectory(prefix="libephemeris-fit-") as name:
        directory = Path(name)
        if single or workers > 1:
            from .processes import fit_bodies_isolated
            bodies, checks = fit_bodies_isolated(source, body_ids, domains, directory, workers, progress)
        else:
            bodies = []
            for body in body_ids:
                bodies.append(fit_body(source, body, *domains[body], directory, checks))
                if progress is not None:
                    progress(body)
        # The file retains the requested domain even if every body's declared
        # overlap is narrower. Knot alignment can extend an angular start.
        actual_start = min(start, *(body.entry.jd_start for body in bodies))
        if header_start is not None:
            if header_start > actual_start:
                raise ValueError("shared header starts after an aligned body domain")
            actual_start = header_start
        auxiliary = {} if skip_aux else auxiliary_sections(actual_start, end, directory)
        header = wire.FileHeader(wire.MAGIC, wire.VERSION, 0, len(bodies), actual_start, end,
                                 2440587.5 + time.time() / 86400.0, 0)
        from .verify import verify_source_file
        write_leb1(output, bodies, header, auxiliary, overwrite=overwrite,
                   validator=lambda path: verify_source_file(path, source, verify_samples))
    return checks


def generate_leb2(source, output: Path, body_ids: list[int], start: float, end: float,
                  *, skip_aux=False, overwrite=False, verify_samples=64, progress=None,
                  lossless=True, target=None, chunk_interval_days=wire.CHUNK_INTERVAL_DAYS,
                  workers=1, single=False, available_domains=False) -> FitChecks:
    """Fit native sources and publish LEB2 only after checking the final reader.

    The private LEB1 intermediate uses the same fitting and source checks as
    standalone LEB1 generation. Conversion is streamed, and its final compressed
    values must pass those source criteria before atomic publication. Lossy
    coefficient budgets never substitute for position/rate verification.
    """
    output = Path(output)
    if output.resolve() in source.input_paths:
        raise ValueError("output must differ from every source and force file")
    if output.exists() and not overwrite:
        raise FileExistsError(output)
    from .convert import convert_leb1
    from .verify import verify_source_file

    with tempfile.TemporaryDirectory(prefix="libephemeris-leb2-fit-") as name:
        intermediate = Path(name) / "coefficients.leb"
        checks = generate_leb1(source, intermediate, body_ids, start, end,
                               skip_aux=skip_aux, verify_samples=verify_samples, progress=progress,
                               workers=workers, single=single, available_domains=available_domains)
        convert_leb1(intermediate, output, lossless=lossless, target=target,
                     chunk_interval_days=chunk_interval_days, overwrite=overwrite,
                     validator=lambda path: verify_source_file(path, source, verify_samples))
    return checks
