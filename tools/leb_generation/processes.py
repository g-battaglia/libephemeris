# SPDX-License-Identifier: AGPL-3.0
"""Fit independent bodies in fresh processes with bounded native lifetimes."""

from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool
from multiprocessing import get_context

from .generate import FitChecks, fit_body
from .source import check_input_hashes


def _fit_body(selection, body_id, start, end, directory):
    """Open local native sources in the child, returning only file metadata."""
    import libephemeris as engine

    try:
        source = selection.open()
        checks = FitChecks()
        body = fit_body(source, body_id, start, end, directory, checks)
        selection.check_inputs()
        return body, checks
    except (OSError, ValueError, engine.Error) as error:
        # Extension exceptions need not be pickleable. Preserve their ordinary
        # diagnostic while returning a standard offline-tool error to the parent.
        raise ValueError(f"body {body_id}: isolated fitting failed: {error}") from None


def fit_bodies_isolated(source, body_ids, domains, directory, workers, progress=None):
    """One body per spawned process; preserve requested output order.

    Each process opens its own verified native handles and exits after its task.
    Only metadata and counts cross process boundaries; coefficient blocks stay
    in the operation's temporary directory. The caller still verifies the final
    assembled reader against its own source before atomic publication.
    """
    selections = {body: source.isolated_selection(body) for body in body_ids}
    for body in body_ids:
        source.require_domain(body, *domains[body])
    results = {}
    try:
        with ProcessPoolExecutor(max_workers=min(workers, len(body_ids)), mp_context=get_context("spawn"),
                                 max_tasks_per_child=1) as pool:
            futures = {pool.submit(_fit_body, selections[body], body, *domains[body], directory): body
                       for body in body_ids}
            try:
                for future in as_completed(futures):
                    body = futures[future]
                    results[body] = future.result()
                    if progress is not None:
                        progress(body)
            except Exception:
                for future in futures:
                    future.cancel()
                raise
    except BrokenProcessPool as error:
        raise ValueError("isolated fitting process exited before returning its result") from error
    check_input_hashes(source.input_hashes().items())
    checks = FitChecks(sum(result[1].positions for result in results.values()),
                       sum(result[1].cartesian_rates for result in results.values()))
    return [results[body][0] for body in body_ids], checks
