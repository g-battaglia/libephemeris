"""Safety and recovery contracts of the resumable LEB build command.

Scientific outputs are either synthetic toy polynomials or small registered JPL
calculations. No reference-distribution data or persisted comparison outputs.
"""

from __future__ import annotations

from dataclasses import replace
import errno
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import pytest

from scripts import regenerate_leb as cli
from scripts.leb_build import plan, runner, sources, storage, validation
from scripts.leb_build.plan import BuildConfig, Job
from scripts.generate_leb import DE441_START_JD, DE441_END_JD, merge_leb_files
from scripts.generate_leb2 import convert_leb1_to_leb2, verify_leb2
from tests.leb_build_helpers import write_synthetic_leb


# ---------------------------------------------------------------------------
# Small deterministic runner fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def config(tmp_path):
    return BuildConfig(
        ("base",),
        str(tmp_path / "data"),
        str(tmp_path / "spk"),
        str(tmp_path / "assist"),
    )


@pytest.fixture
def mini_build(tmp_path, config, monkeypatch):
    root = tmp_path / "build"
    root.mkdir()
    jobs = [
        Job("generate", "base", "generate", "work/body_0.leb", (0,)),
        Job("verify", "base", "verify1", "work/body_0.leb", (0,), ("generate",)),
        Job(
            "merge",
            "base",
            "merge",
            "leb/merged.leb",
            (0,),
            ("verify",),
            ("work/body_0.leb",),
            export=True,
        ),
    ]
    inputs = {"environment": {}, "files": {}}
    manifest = runner.new_manifest(config, jobs, inputs, {})
    monkeypatch.setattr(runner, "input_stamps", lambda _: {"fixed": (1,)})
    monkeypatch.setattr(runner, "attest_inputs", lambda _: inputs)
    monkeypatch.setattr(validation, "tier_range", lambda _: (100.0, 110.0))
    return root, jobs, manifest


def synthetic_executor(root, calls):
    """Execute toy generation and the real merge; verification is a fake phase."""

    def execute(job, target, log):
        calls.append(job.id)
        if job.kind == "generate":
            write_synthetic_leb(target)
        elif job.kind == "merge":
            merge_leb_files(
                [str(root / path) for path in job.inputs], str(target), False
            )
        return 0

    return execute


def run_mini(root, config, jobs, manifest, executor):
    with storage.build_lock(root) as descriptor:
        runner.BuildRunner(root, config, jobs, manifest, descriptor, executor).run()


# ---------------------------------------------------------------------------
# Inventory, settings and read-only command behavior
# ---------------------------------------------------------------------------


def test_canonical_plan_has_all_exports_and_per_body_checkpoints(config):
    jobs = plan.build_jobs(replace(config, tiers=("base", "medium", "extended")))
    exports = [job for job in jobs if job.export]
    assert len(exports) == 27
    assert sum(job.output.endswith(".leb") for job in exports) == 15
    assert sum(job.output.endswith(".leb2") for job in exports) == 12
    generated = [job for job in jobs if job.kind == "generate"]
    assert len(generated) == 151
    assert [
        sum(job.tier == tier for job in generated)
        for tier in config.tiers + ("medium", "extended")
    ] == [53, 53, 45]
    seen = set()
    for job in jobs:
        assert set(job.dependencies) <= seen
        seen.add(job.id)
        assert "uranians" not in job.output
    assert plan.tier_range("extended") == (DE441_START_JD, DE441_END_JD)


def test_scientific_environment_overrides_local_lebs_and_network(config, monkeypatch):
    monkeypatch.setenv("LIBEPHEMERIS_MODE", "leb")
    monkeypatch.setenv("LIBEPHEMERIS_LEB", "/old/model.leb")
    monkeypatch.setenv("LIBEPHEMERIS_NETWORK_POLICY", "allow")
    env = plan.scientific_environment(config, "base")
    assert env["LIBEPHEMERIS_MODE"] == "skyfield"
    assert env["LIBEPHEMERIS_NETWORK_POLICY"] == "sealed"
    assert env["LIBEPHEMERIS_CONFIG"] == os.devnull
    assert env["LIBEPHEMERIS_ENV_FILE"] == os.devnull
    assert env["LIBEPHEMERIS_PRECISION"] == "base"
    assert "LIBEPHEMERIS_LEB" not in env


def test_dry_run_does_not_create_output_or_access_sources(
    tmp_path, monkeypatch, capsys
):
    def forbidden(*args):
        pytest.fail("dry-run accessed source preflight")

    monkeypatch.setattr(cli, "preflight", forbidden)
    output = tmp_path / "preview"
    assert cli.main(["--output-dir", str(output), "--dry-run"]) == 0
    assert "27 canonical exports" in capsys.readouterr().out
    assert not output.exists()


def test_doctor_gate_failure_does_not_create_output(tmp_path, monkeypatch):
    def reject(*args):
        raise ValueError("provenance gate failed")

    monkeypatch.setattr(cli, "preflight", reject)
    output = tmp_path / "build"
    assert cli.main(["--output-dir", str(output), "--doctor"]) == 1
    assert not output.exists()


def test_resume_rejects_changed_options(config):
    args = cli.build_parser().parse_args(
        ["--output-dir", "/unused", "--resume", "--verify-samples", "1"]
    )
    with pytest.raises(ValueError, match="Resume options"):
        cli.resolve_config(args, config)


# ---------------------------------------------------------------------------
# Recovery: generation, verification, corruption and orphan files
# ---------------------------------------------------------------------------


def test_complete_resume_runs_no_scientific_phase(mini_build, config):
    root, jobs, manifest = mini_build
    calls = []
    execute = synthetic_executor(root, calls)
    run_mini(root, config, jobs, manifest, execute)
    assert manifest["status"] == "complete"
    assert not storage.invalidate(manifest, jobs, set())
    assert not runner.recover_checkpoints(root, manifest, jobs)
    run_mini(root, config, jobs, manifest, execute)
    assert calls == ["generate", "verify", "merge"]
    assert (root / "checksums.sha256").read_text().endswith("  leb/merged.leb\n")


def test_failed_verification_preserves_generated_checkpoint(mini_build, config):
    root, jobs, manifest = mini_build
    calls = []
    normal = synthetic_executor(root, calls)

    def fail_verify(job, target, log):
        return 7 if job.id == "verify" else normal(job, target, log)

    with pytest.raises(ValueError, match="exit 7"):
        run_mini(root, config, jobs, manifest, fail_verify)
    generated_hash = storage.sha256_file(root / jobs[0].output)
    assert manifest["jobs"]["generate"]["status"] == "done"
    assert runner.recover_checkpoints(root, manifest, jobs) == {"verify", "merge"}
    run_mini(root, config, jobs, manifest, normal)
    assert calls.count("generate") == 1
    assert storage.sha256_file(root / jobs[0].output) == generated_hash


@pytest.mark.parametrize("phase", ["generate", "verify", "merge"])
@pytest.mark.parametrize("signum", [signal.SIGINT, signal.SIGTERM])
def test_interruption_is_checkpointed_and_resumable(mini_build, config, phase, signum):
    root, jobs, manifest = mini_build
    calls = []
    normal = synthetic_executor(root, calls)

    def interrupt(job, target, log):
        if job.id == phase:
            if job.writes_output:
                target.write_bytes(b"unfinished")
            raise runner.BuildInterrupted(signum)
        return normal(job, target, log)

    with pytest.raises(runner.BuildInterrupted) as raised:
        run_mini(root, config, jobs, manifest, interrupt)
    assert raised.value.signum == signum
    assert manifest["status"] == "interrupted"
    assert manifest["jobs"][phase]["status"] == "interrupted"
    assert not list(root.rglob("*.part-*"))
    runner.recover_checkpoints(root, manifest, jobs)
    run_mini(root, config, jobs, manifest, normal)
    assert manifest["status"] == "complete"


@pytest.mark.parametrize("damage", ["corrupt", "missing"])
def test_damaged_checkpoint_invalidates_descendants(mini_build, config, damage):
    root, jobs, manifest = mini_build
    run_mini(root, config, jobs, manifest, synthetic_executor(root, []))
    path = root / jobs[0].output
    if damage == "missing":
        path.unlink()
    else:
        contents = bytearray(path.read_bytes())
        contents[250] ^= 1
        path.write_bytes(contents)
    assert runner.recover_checkpoints(root, manifest, jobs) == {
        "generate",
        "verify",
        "merge",
    }
    assert (root / "leb/merged.leb").exists()


def test_orphan_final_and_staging_are_not_adopted(mini_build, config):
    root, jobs, manifest = mini_build
    write_synthetic_leb(root / jobs[0].output)
    record = manifest["jobs"]["generate"]
    record.update(status="running", attempt=1)
    stage = storage.staging_path(root, jobs[0], 1)
    stage.write_bytes(b"stale")
    runner.recover_checkpoints(root, manifest, jobs)
    calls = []
    run_mini(root, config, jobs, manifest, synthetic_executor(root, calls))
    assert calls == ["generate", "verify", "merge"]
    assert not stage.exists()


def test_manifest_definition_cannot_redirect_output(mini_build):
    _, jobs, manifest = mini_build
    manifest["jobs"]["generate"]["definition"]["output"] = "../outside"
    with pytest.raises(ValueError, match="definition mismatch"):
        storage.validate_job_records(manifest, jobs)


def test_source_change_stops_before_next_phase(mini_build, config, monkeypatch):
    root, jobs, manifest = mini_build
    normal = synthetic_executor(root, [])

    def edit_source(job, target, log):
        code = normal(job, target, log)
        monkeypatch.setattr(runner, "input_stamps", lambda _: {"fixed": (2,)})
        return code

    with pytest.raises(ValueError, match="inputs changed"):
        run_mini(root, config, jobs, manifest, edit_source)
    assert manifest["jobs"]["verify"]["status"] == "pending"


def test_disk_full_keeps_previous_manifest_and_cleans_temp(tmp_path, monkeypatch):
    path = tmp_path / "manifest.json"
    path.write_text("old")

    def fail(*args):
        raise OSError(errno.ENOSPC, "No space left")

    monkeypatch.setattr(storage.os, "replace", fail)
    with pytest.raises(OSError):
        storage.atomic_write(path, "new")
    assert path.read_text() == "old"
    assert not list(tmp_path.glob(".manifest.json.*"))


# ---------------------------------------------------------------------------
# Source selection, conversion and structural acceptance
# ---------------------------------------------------------------------------


def test_anchor_only_spk_fails_preflight_coverage(config, monkeypatch):
    from scripts import generate_leb
    from libephemeris.minor_bodies import HORIZONS_SPK_JD_MIN, HORIZONS_SPK_JD_MAX

    monkeypatch.setattr(
        sources, "tier_groups", lambda _: {"asteroids": (15,), "exotics": ()}
    )
    monkeypatch.setattr(
        sources, "tier_range", lambda _: (HORIZONS_SPK_JD_MIN, HORIZONS_SPK_JD_MAX)
    )
    monkeypatch.setattr(
        sources, "minor_spk_candidates", lambda _: {15: [Path("/anchor.bsp")]}
    )
    monkeypatch.setattr(
        generate_leb, "_get_asteroid_spk_range", lambda *_: (2451000.0, 2452000.0)
    )
    with pytest.raises(ValueError, match="Missing usable local SPK"):
        sources.select_minor_spks(config)


def test_useful_partial_spk_coverage_is_allowed(config, monkeypatch):
    from scripts import generate_leb

    monkeypatch.setattr(
        sources, "tier_groups", lambda _: {"asteroids": (15,), "exotics": ()}
    )
    monkeypatch.setattr(sources, "tier_range", lambda _: (0.0, 100_000.0))
    monkeypatch.setattr(
        sources, "minor_spk_candidates", lambda _: {15: [Path("/partial.bsp")]}
    )
    monkeypatch.setattr(
        generate_leb, "_get_asteroid_spk_range", lambda *_: (1000.0, 10_000.0)
    )
    assert sources.select_minor_spks(config) == {"base": {"15": "/partial.bsp"}}


def test_extended_nbody_seed_requires_wide_verification_coverage(config, monkeypatch):
    from scripts import generate_leb

    body = min(
        set(generate_leb.EXOTIC_EXTENDED_IDS)
        - set(generate_leb.EXOTIC_ASSIST_PERTURBER_IDS)
    )
    monkeypatch.setattr(
        sources, "tier_groups", lambda _: {"asteroids": (), "exotics": (body,)}
    )
    monkeypatch.setattr(sources, "tier_range", lambda _: (DE441_START_JD, DE441_END_JD))
    monkeypatch.setattr(
        sources, "minor_spk_candidates", lambda _: {body: [Path("/seed.bsp")]}
    )
    monkeypatch.setattr(
        generate_leb, "_get_asteroid_spk_range", lambda *_: (2450000.0, 2460000.0)
    )
    with pytest.raises(ValueError, match="Missing usable local SPK"):
        sources.select_minor_spks(replace(config, tiers=("extended",)))


def test_minor_source_selection_is_per_tier(config, monkeypatch):
    from scripts import generate_leb

    paths = [Path("/base.bsp"), Path("/medium.bsp")]
    monkeypatch.setattr(
        sources, "tier_groups", lambda _: {"asteroids": (15,), "exotics": ()}
    )
    monkeypatch.setattr(
        sources,
        "tier_range",
        lambda tier: (0.0, 10_000.0) if tier == "base" else (20_000.0, 30_000.0),
    )
    monkeypatch.setattr(sources, "minor_spk_candidates", lambda _: {15: paths})
    monkeypatch.setattr(
        generate_leb,
        "_get_asteroid_spk_range",
        lambda path, _: (0.0, 10_000.0)
        if path == "/base.bsp"
        else (20_000.0, 30_000.0),
    )
    assert sources.select_minor_spks(replace(config, tiers=("base", "medium"))) == {
        "base": {"15": "/base.bsp"},
        "medium": {"15": "/medium.bsp"},
    }


def test_converter_and_verifier_accept_small_synthetic_group(tmp_path, monkeypatch):
    reference, output = tmp_path / "source.leb", tmp_path / "base_asteroids.leb2"
    bodies = (15, 17, 18, 19, 20)
    write_synthetic_leb(reference, bodies)
    convert_leb1_to_leb2(
        str(reference), str(output), "asteroids", "base", verbose=False
    )
    assert verify_leb2(
        str(output),
        str(reference),
        n_samples=10,
        verbose=False,
        expected_group="asteroids",
        expected_tier="base",
    )
    monkeypatch.setattr(validation, "tier_range", lambda _: (100.0, 110.0))
    job = Job("convert", "base", "convert", "unused", bodies, group="asteroids")
    assert validation.inspect_artifact(output, job)["format"] == "LEB2 v2"
    with pytest.raises(ValueError, match="inventory"):
        validation.inspect_artifact(output, replace(job, bodies=(15,)))


def test_worker_passes_reference_inventory_and_sample_count(
    tmp_path, config, monkeypatch
):
    from scripts import generate_leb2

    captured = {}
    monkeypatch.setattr(runner, "configure_runtime", lambda *args: None)

    def verify(path, **kwargs):
        captured.update(path=path, **kwargs)
        return False

    monkeypatch.setattr(generate_leb2, "verify_leb2", verify)
    job = Job(
        "verify2",
        "base",
        "verify2",
        "leb2/base_core.leb2",
        (0,),
        inputs=("leb/base.leb",),
        group="core",
    )
    with pytest.raises(ValueError, match="Scientific verification"):
        runner.execute_scientific_job(tmp_path, config, job, tmp_path / job.output, {})
    assert captured["reference_leb1"] == str(tmp_path / "leb/base.leb")
    assert captured["expected_group"] == "core"
    assert captured["expected_tier"] == "base"
    assert captured["n_samples"] == config.leb2_verify_samples


def test_small_local_jpl_pipeline_and_real_worker(tmp_path, monkeypatch):
    """Exercise fitting, streaming merge, compression and the worker entry point.

    Only Sun over 64 days is generated. A missing local kernel skips this test;
    the sealed policy prevents provisioning or network access during the test.
    """
    from libephemeris.state import _resolve_data_dir
    from scripts import generate_leb

    data_dir = Path(_resolve_data_dir()).resolve()
    if not (data_dir / "de440s.bsp").is_file():
        pytest.skip("Local DE440s not provisioned")
    config = BuildConfig(
        ("base",),
        str(data_dir),
        str(data_dir / "spk"),
        str(data_dir / "assist"),
        verify_samples=10,
        leb2_verify_samples=10,
    )
    for key, value in plan.scientific_environment(config, "base").items():
        if key.startswith("LIBEPHEMERIS_") or key == "ASSIST_DIR":
            monkeypatch.setenv(key, value)
    start = generate_leb.J2000
    monkeypatch.setattr(runner, "tier_range", lambda _: (start, start + 64))
    root = tmp_path / "small-science"
    root.mkdir()
    source = root / "work/base/body_0.leb"
    source.parent.mkdir(parents=True)
    generate_job = Job(
        "small.generate", "base", "generate", str(source.relative_to(root)), (0,)
    )
    runner.execute_scientific_job(root, config, generate_job, source, {})
    merged = root / "merged.leb"
    merge_job = Job(
        "small.merge",
        "base",
        "merge",
        "merged.leb",
        (0,),
        inputs=(generate_job.output,),
    )
    runner.execute_scientific_job(root, config, merge_job, merged, {})
    compressed = root / "small.leb2"
    convert_leb1_to_leb2(str(merged), str(compressed), verbose=False)
    assert verify_leb2(str(compressed), str(merged), n_samples=10, verbose=False)

    jobs = plan.build_jobs(config)
    manifest = runner.new_manifest(config, jobs, {"environment": {}, "files": {}}, {})
    record = manifest["jobs"]["base.body.0.verify"]
    record.update(status="running", attempt=1)
    manifest["jobs"]["base.body.0.generate"]["status"] = "done"
    storage.save_manifest(root, manifest)
    with storage.build_lock(root) as descriptor:
        env = plan.scientific_environment(config, "base")
        env["LEB_BUILD_LOCK_FD"] = str(descriptor)
        command = [
            sys.executable,
            "-B",
            str(plan.PROJECT_ROOT / "scripts/regenerate_leb.py"),
            "--output-dir",
            str(root),
            "--worker-job",
            "base.body.0.verify",
            "--worker-output",
            str(source),
        ]
        log = root / "worker.log"
        assert runner.run_process(command, env, log, descriptor) == 0, log.read_text()


# ---------------------------------------------------------------------------
# Real POSIX process lifetimes and path boundaries
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("relative", ["../escape", "/absolute"])
def test_managed_paths_reject_escapes(tmp_path, relative):
    with pytest.raises(ValueError):
        storage.managed_path(tmp_path, relative)


def test_output_cannot_overlap_sources_or_follow_symlinks(tmp_path, config):
    with pytest.raises(ValueError, match="overlaps"):
        storage.safe_root(Path(config.data_dir) / "build", config)
    with pytest.raises(ValueError, match="overlaps"):
        storage.safe_root(tmp_path, config)
    root = tmp_path / "build"
    root.mkdir()
    (root / "leb").symlink_to(tmp_path, target_is_directory=True)
    with pytest.raises(ValueError, match="Symlink"):
        storage.managed_path(root, "leb/file.leb")
    with pytest.raises(ValueError, match="path type"):
        storage.safe_root(root, config)


def test_child_retains_lock_when_parent_closes_descriptor(tmp_path):
    with storage.build_lock(tmp_path) as descriptor:
        child = subprocess.Popen(
            [sys.executable, "-B", "-c", "import sys; sys.stdin.read()"],
            stdin=subprocess.PIPE,
            pass_fds=(descriptor,),
            start_new_session=True,
        )
    try:
        with pytest.raises(ValueError, match="holds this build lock"):
            with storage.build_lock(tmp_path):
                pass
    finally:
        child.communicate(input=b"", timeout=5)
    with storage.build_lock(tmp_path):
        pass
    assert (tmp_path / ".lock").exists()


@pytest.mark.parametrize("signum", [signal.SIGINT, signal.SIGTERM])
def test_real_worker_is_stopped_before_parent_returns(tmp_path, signum):
    """Signal a runner controlling a real worker and prove its lock is released."""
    marker = tmp_path / "pid"
    child_code = "import os, pathlib, sys; pathlib.Path(sys.argv[1]).write_text(str(os.getpid())); sys.stdin.read()"
    # The worker must keep running rather than reading the runner's closed stdin.
    child_code = child_code.replace("sys.stdin.read()", "__import__('time').sleep(120)")
    parent_code = (
        "from pathlib import Path; import sys; "
        "from scripts.leb_build.storage import build_lock; "
        "from scripts.leb_build.runner import interruption_signals, run_process, BuildInterrupted\n"
        "root=Path(sys.argv[1])\n"
        "try:\n"
        " with build_lock(root) as fd, interruption_signals():\n"
        "  run_process([sys.executable,'-B','-c',sys.argv[2],str(root/'pid')], "
        "dict(__import__('os').environ),root/'log',fd)\n"
        "except BuildInterrupted as exc:\n"
        " sys.exit(128+exc.signum)\n"
    )
    parent = subprocess.Popen(
        [sys.executable, "-B", "-c", parent_code, str(tmp_path), child_code],
        cwd=plan.PROJECT_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )
    worker_pid = None
    try:
        deadline = time.monotonic() + 10
        while not marker.exists() and time.monotonic() < deadline:
            if parent.poll() is not None:
                pytest.fail(str(parent.communicate()))
            time.sleep(0.02)
        assert marker.exists()
        worker_pid = int(marker.read_text())
        parent.send_signal(signum)
        stdout, stderr = parent.communicate(timeout=10)
        assert parent.returncode == 128 + signum, (stdout, stderr)
        with pytest.raises(ProcessLookupError):
            os.kill(worker_pid, 0)
        with storage.build_lock(tmp_path):
            pass
    finally:
        if parent.poll() is None:
            parent.kill()
            parent.wait(timeout=5)
        if worker_pid:
            try:
                os.killpg(worker_pid, signal.SIGKILL)
            except ProcessLookupError:
                pass


def test_spk_target_is_read_from_file_with_repeated_segments(monkeypatch):
    """Horizons kernels repeat one 20000000+N target across many segments."""
    from libephemeris import spk

    monkeypatch.setattr(spk, "_get_spk_targets", lambda _: [20002060] * 43)
    assert sources.spk_target_id("chiron.bsp", 2060) == 20002060
    monkeypatch.setattr(spk, "_get_spk_targets", lambda _: [2060, 2060])
    assert sources.spk_target_id("chiron.bsp", 2060) == 2060
    monkeypatch.setattr(spk, "_get_spk_targets", lambda _: [])
    with pytest.raises(ValueError, match="single NAIF target"):
        sources.spk_target_id("empty.bsp", 2060)
