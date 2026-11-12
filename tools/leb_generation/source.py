# SPDX-License-Identifier: AGPL-3.0
"""Native source states and models for standalone offline coefficient fitting."""

import hashlib
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

PLANET_NAIF = {0: 10, 1: 301, 2: 199, 3: 299, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 14: 399}
ANGULAR_IDS = frozenset((10, 11, 12, 13, 21, 22))
IERS_FILES = ("finals2000A.data", "leap_seconds.dat", "deltat.data")


def _input_hash(path):
    try:
        with path.open("rb") as stream:
            return hashlib.file_digest(stream, "sha256").hexdigest()
    except FileNotFoundError:
        return None


def check_input_hashes(identities):
    """Reject changed bytes or changed absence of an optional local clock file."""
    for path, expected in identities:
        if _input_hash(path) != expected:
            raise ValueError(f"generation source changed: {path}")


@dataclass(frozen=True)
class NativeSourceSelection:
    """Serializable local selections; native handles never cross processes."""

    kernel: Path
    tier: str
    delta_t_model: str
    spk: tuple[tuple[int, Path, int], ...]
    nbody: tuple[tuple[int, Path, Path], ...]
    hashes: tuple[tuple[Path, str | None], ...]
    delta_t_userdef_days: float | None = None
    iers_cache_dir: Path | None = None

    def open(self):
        source = NativeSource(self.kernel, self.tier, self.delta_t_model,
                              delta_t_userdef_days=self.delta_t_userdef_days,
                              iers_cache_dir=self.iers_cache_dir)
        for body, path, target in self.spk:
            source.add_spk_source(body, path, target)
        for body, planets, asteroids in self.nbody:
            source.add_nbody_sources([body], planets, asteroids)
        if source.input_hashes() != dict(self.hashes):
            raise ValueError("isolated generation source identity changed")
        return source

    def check_inputs(self):
        check_input_hashes(self.hashes)


def merge_windows(windows):
    """Union overlapping closed intervals without bridging a coverage gap."""
    result = []
    for start, end in sorted(windows):
        if result and start <= result[-1][1]:
            result[-1] = result[-1][0], max(end, result[-1][1])
        else:
            result.append((start, end))
    return result


def continuous_overlap(windows, start, end):
    """Return the complete unique overlap, never bridge a gap or drop an island.

    A LEB body index describes one closed interval. Multiple disjoint source
    intervals therefore require separate files or an explicit narrower request;
    selecting just the largest would conceal available dates from the caller.
    """
    if not (math.isfinite(start) and math.isfinite(end) and start < end):
        raise ValueError("requested generation domain must be finite with positive width")
    overlaps = merge_windows((max(start, low), min(end, high)) for low, high in windows
                             if max(start, low) < min(end, high))
    if not overlaps:
        raise ValueError("requested range has no continuous source overlap")
    if len(overlaps) != 1:
        raise ValueError("requested range crosses a source gap; select one explicit continuous interval")
    return overlaps[0]


def chain_windows(descriptors, target, seen=()):
    """Intersect each target descriptor with its complete stored center chain.

    Metadata qualifies the requested domain before fitting. Actual native
    queries still enforce priority, frames, coefficient domains and finiteness.
    Neither a disjoint span nor a target-only span claims continuous support.
    """
    if target == 0:
        return [(-math.inf, math.inf)]
    if target in seen:
        raise ValueError("cyclic SPK center chain")
    result = []
    for identifier, center, frame, low, high in descriptors:
        if identifier != target:
            continue
        if frame != 1:
            raise ValueError(f"generation requires J2000 SPK frame, found {frame}")
        for start, end in chain_windows(descriptors, center, (*seen, target)):
            if max(start, low) <= min(end, high):
                result.append((max(start, low), min(end, high)))
    return merge_windows(result)


class NativeSource:
    """Own one exact local JPL image; evaluate all source astronomy in Rust.

    This class is for the standalone generation process. It explicitly selects
    its global native model configuration for angular/auxiliary APIs; it is
    not imported by the runtime package. Network and automatic downloads are
    disabled. Existing LEB coefficients cannot become circular fitting input.
    """

    def __init__(self, kernel: Path, tier: str, delta_t_model="smh2016", *,
                 delta_t_userdef_days=None, iers_cache_dir=None):
        import libephemeris as engine
        from libephemeris import iers_data, lunar, mean_lunar_apse

        self.kernel = Path(kernel).resolve(strict=True)
        self.tier = tier
        self.delta_t_model = delta_t_model
        if delta_t_userdef_days is not None and (not math.isfinite(delta_t_userdef_days)
                                                or delta_t_userdef_days == engine.DELTAT_AUTOMATIC):
            raise ValueError("--delta-t-userdef-days requires a finite explicit value, not DELTAT_AUTOMATIC")
        self.delta_t_userdef_days = delta_t_userdef_days
        self.iers_cache_dir = None if iers_cache_dir is None else Path(iers_cache_dir).resolve(strict=True)
        if self.iers_cache_dir is not None and not self.iers_cache_dir.is_dir():
            raise ValueError("--iers-cache-dir must be an existing local directory")
        # Bind absent inputs too: adding a monthly or leap file changes the
        # native selection just as replacing an existing file does.
        self.clock_hashes = {} if self.iers_cache_dir is None else {
            self.iers_cache_dir / name: _input_hash(self.iers_cache_dir / name) for name in IERS_FILES
        }
        with self.kernel.open("rb") as source:
            self.sha256 = hashlib.file_digest(source, "sha256").hexdigest()
        engine.close()
        engine.set_calc_mode("skyfield")
        engine.set_precision_tier(tier)
        engine.set_network_policy("sealed")
        engine.set_auto_spk_download(False)
        engine.set_iers_auto_download(False)
        engine.set_iers_delta_t_enabled(self.iers_cache_dir is not None)
        engine.set_delta_t_userdef(delta_t_userdef_days)
        engine.set_delta_t_model(delta_t_model)
        if self.iers_cache_dir is not None:
            # This native setter clears the old row generation. Reload alone
            # could retain preceding finals rows when the new file is absent.
            iers_data.set_iers_cache_dir(str(self.iers_cache_dir))
            iers_data.load_iers_data(force_download=False)
            check_input_hashes(self.clock_hashes.items())
        engine.set_leb_file(None)
        engine.set_ephe_path(str(self.kernel.parent))
        engine.set_ephemeris_file(str(self.kernel))
        self.image = engine.get_loader()(str(self.kernel))
        self.clock = engine.get_timescale()
        if Path(self.image.path).resolve() != self.kernel:
            raise ValueError("source loader selected a different kernel")
        self.descriptors = self.image.segment_descriptors_tdb()
        # The scalar lunar model APIs use the process's explicit provider map.
        # Attach this same image once so their successive snapshots share it;
        # otherwise their metadata factory would reopen a cold main per point.
        # NAIF aliases are local to this standalone generation process.
        for naif_id in sorted({row[0] for row in self.descriptors}):
            engine.register_spk_body(naif_id, str(self.kernel), naif_id)
        self.targets = dict(PLANET_NAIF)
        # Mars's center is used when this actual image contains it. Outer
        # planets deliberately store system barycenters (coordinate type 4),
        # leaving runtime center corrections to the existing native pipeline.
        if any(row[0] == 499 for row in self.descriptors):
            self.targets[4] = 499
        self.angular = {
            10: mean_lunar_apse.mean_lunar_node_position,
            11: lunar.calc_true_lunar_node,
            12: mean_lunar_apse.mean_lunar_apogee_position,
            13: lunar.calc_true_lilith,
            21: lunar.calc_interpolated_apogee,
            22: lunar.calc_interpolated_perigee,
        }
        self.piecewise_grids = {body: grid for body in ANGULAR_IDS
                               if (grid := lunar.interpolated_apse_source_grid(body)) is not None}
        self.spk_sources = {}
        self.spk_inputs = {}
        self.input_paths = {self.kernel, *(path.resolve() for path in self.clock_hashes)}
        self.nbody_sources = {}
        self.nbody_inputs = {}
        self.spk_selections = {}
        self.nbody_selections = {}
        self.ignored_non_spk_files = set()

    def input_hashes(self):
        """Main/target/force bytes and the selected local clock files (or absence)."""
        result = {self.kernel: self.sha256, **self.clock_hashes}
        groups = [data["images"] for data in self.spk_inputs.values()]
        groups.extend(data["force_images"] for data in self.nbody_inputs.values())
        for images in groups:
            for name, digest in images.items():
                path = Path(name)
                if path in result and result[path] != digest:
                    raise ValueError(f"inconsistent generation source identity: {path}")
                result[path] = digest
        return result

    def isolated_selection(self, body_id):
        """Retain only one body's explicit sources for a fresh worker process."""
        spk = ((body_id, *self.spk_selections[body_id]),) if body_id in self.spk_selections else ()
        nbody = ((body_id, *self.nbody_selections[body_id]),) if body_id in self.nbody_selections else ()
        paths = {self.kernel, *self.clock_hashes}
        if spk:
            paths.update(Path(name) for name in self.spk_inputs[body_id]["images"])
        if nbody:
            paths.update(Path(name) for name in self.nbody_inputs[body_id]["force_images"])
        identities = self.input_hashes()
        return NativeSourceSelection(self.kernel, self.tier, self.delta_t_model, spk, nbody,
                                     tuple((path, identities[path]) for path in sorted(paths)),
                                     self.delta_t_userdef_days, self.iers_cache_dir)

    def add_nbody_sources(self, bodies, planets: Path, asteroids: Path):
        """Use explicit SPK targets only as J2000 seeds for native integration.

        Both exact existing force files are mandatory. Rust owns the frame,
        anchor translation, guarded force coverage and directional integration.
        The native session is retained by every source; no global runtime
        fallback, catalogue initial condition or implicit data selection enters.
        """
        from libephemeris.rebound_integration import AssistEphemConfig

        bodies = list(bodies)
        if not bodies or len(set(bodies)) != len(bodies):
            raise ValueError("N-body sources require distinct explicit body IDs")
        if any(body not in self.spk_sources or body in self.nbody_sources for body in bodies):
            raise ValueError("every N-body source requires one unassigned explicit --spk seed")
        planets, asteroids = Path(planets).resolve(strict=True), Path(asteroids).resolve(strict=True)
        configuration = AssistEphemConfig(str(planets), str(asteroids))
        if (Path(configuration.planets_file).resolve(), Path(configuration.asteroids_file).resolve()) != (planets, asteroids):
            raise ValueError("ASSIST selected an unrequested force file")
        images = {}
        for path in (planets, asteroids):
            with path.open("rb") as stream:
                images[str(path)] = hashlib.file_digest(stream, "sha256").hexdigest()
        sources = {body: configuration.nbody_source(self.spk_sources[body]) for body in bodies}
        for body, source in sources.items():
            self.nbody_sources[body] = source
            self.nbody_inputs[body] = {"anchor_tdb": source.anchor_tdb,
                                      "coverage_tdb": source.coverage_tdb(),
                                      "force_images": images,
                                      "seed": self.spk_inputs[body]}
            self.nbody_selections[body] = (planets, asteroids)
        self.input_paths.update((planets, asteroids))

    def add_spk_source(self, body_id: int, path: Path, naif_id: int):
        """Attach an exact local minor-body SPK and the retained main support.

        Native code composes its geometric center chain and intersects actual
        stored coverage windows. Input images are hashed before fitting; no
        replacement image, TT conversion or model fallback enters this route.
        The same explicit body registration serves public validation calls.
        """
        import libephemeris as engine

        from .format import BODY_PARAMS

        if body_id not in BODY_PARAMS or body_id in self.targets or body_id in self.angular:
            raise ValueError("an explicit minor-body SPK requires a supported minor-body ID")
        if body_id in self.spk_sources:
            raise ValueError(f"duplicate minor-body source: {body_id}")
        if naif_id == 0:
            raise ValueError("minor-body generation requires an actual SPK target, not the SSB origin")
        path = Path(path).resolve(strict=True)
        primary = engine.get_loader()(str(path))
        if Path(primary.path).resolve() != path:
            raise ValueError("minor-body loader selected a different image")
        source = primary.geometric_source(naif_id, [self.image])
        inputs = {}
        for name in source.paths:
            image = Path(name).resolve(strict=True)
            if image not in (path, self.kernel):
                raise ValueError("geometric source selected an unrequested supporting image")
            with image.open("rb") as stream:
                inputs[str(image)] = hashlib.file_digest(stream, "sha256").hexdigest()
        engine.register_spk_body(body_id, str(path), naif_id)
        self.spk_sources[body_id] = source
        self.spk_inputs[body_id] = {"naif_id": naif_id, "images": inputs, "coverage_tdb": source.coverage_tdb()}
        self.spk_selections[body_id] = (path, naif_id)
        self.input_paths.update(Path(name) for name in inputs)

    def coverage(self, body_id: int):
        """Complete native source windows on the LEB format's TT grid.

        N-body windows already include the force-stage guard; SPK windows
        include every consumed center. An analytical angular model has no
        finite stored-image boundary, and is still sampled/checked by its API.
        """
        if body_id in self.nbody_sources:
            windows = [self.nbody_sources[body_id].coverage_tdb()]
        elif body_id in self.targets:
            windows = chain_windows(self.descriptors, self.targets[body_id])
        elif body_id in self.spk_sources:
            windows = self.spk_sources[body_id].coverage_tdb()
        elif body_id in self.angular:
            return [(-math.inf, math.inf)]
        else:
            raise ValueError(f"body {body_id}: supply its explicit minor-body --spk source")
        return [self.clock.tdb_interval_to_tt(start, end) for start, end in windows]

    def available_domain(self, body_id: int, start: float, end: float):
        """Explicitly requested per-body overlap; missing/gapped data is an error."""
        try:
            return continuous_overlap(self.coverage(body_id), start, end)
        except ValueError as error:
            raise ValueError(f"body {body_id}: {error}") from error

    def require_domain(self, body_id: int, start: float, end: float):
        if not (math.isfinite(start) and math.isfinite(end) and start < end):
            raise ValueError("requested generation domain must be finite with positive width")
        windows = self.coverage(body_id)
        if not any(low <= start < end <= high for low, high in windows):
            raise ValueError(f"body {body_id}: requested range lacks a complete continuous source chain")

    def values(self, body_id: int, epochs):
        """Positions on the TT fitting grid: ICRS AU, or angular deg/deg/AU."""
        if body_id in self.angular:
            return np.asarray([self.angular[body_id](float(epoch)) for epoch in epochs], dtype=np.float64)
        return self.states(body_id, epochs)[:, :3]

    def states(self, body_id: int, epochs):
        """Cartesian TT rows x/y/z/vx/vy/vz in AU and AU per TT day.

        Rust converts the clock before sampling SPK/integrated TDB states and
        applies the analytic velocity chain rule. Raw TDB APIs stay separate.
        """
        epochs = [float(epoch) for epoch in epochs]
        if body_id in self.nbody_sources:
            rows = self.nbody_sources[body_id].states_tt(epochs)
        elif body_id in self.spk_sources:
            rows = self.spk_sources[body_id].states_tt(epochs)
        elif body_id in self.targets:
            rows = self.image.geometric_states_tt(self.targets[body_id], epochs)
        else:
            raise ValueError(f"body {body_id}: supply its explicit minor-body --spk source")
        return np.asarray(rows, dtype=np.float64).reshape(-1, 6)


def configure_spk_sources(source: NativeSource, specifications):
    """Apply explicit CLI BODY_ID/PATH/NAIF_ID triples without inference."""
    for body, path, naif in specifications:
        source.add_spk_source(int(body), Path(path), int(naif))
