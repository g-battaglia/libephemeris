"""Offline orchestration tests; fitting and native source evaluation are replaced."""

from types import SimpleNamespace

import libephemeris as engine
import pytest

from tools.leb_generation import cli


@pytest.fixture
def workflow(tmp_path, monkeypatch):
    events = []
    inventories = {}

    def create_source(args, native, bodies, domains, **options):
        assert native is engine
        events.append(("source", args.tier, tuple(bodies), domains, options))
        return SimpleNamespace(
            input_paths={tmp_path / "source.bsp"}, tier=args.tier, piecewise_grids={},
            require_domain=lambda body, start, end: events.append(("admit", body, start, end)),
            available_domain=lambda body, start, end: (start, end),
        )

    def generate(source, output, bodies, start, end, **options):
        events.append(("fit", source.tier, output, tuple(bodies), start, end, options))
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_bytes(b"fitted")
        inventories[output] = list(bodies)
        return SimpleNamespace(positions=9, cartesian_rates=9)

    def merge(parts, output, **options):
        events.append(("merge", tuple(parts), output, options))
        output.write_bytes(b"merged")

    monkeypatch.setattr(cli, "create_source", create_source)
    monkeypatch.setattr(cli, "generate_leb1", generate)
    monkeypatch.setattr(cli, "merge_leb1", merge)
    monkeypatch.setattr(cli, "read_body_index", lambda path: (
        None, [SimpleNamespace(body_id=body) for body in inventories[path]]))
    monkeypatch.setattr(cli, "describe_source", lambda source: None)
    monkeypatch.setattr(cli, "describe_domains", lambda path: None)
    return events, inventories


def arguments(tmp_path, *selection):
    return [*selection, "--data-dir", str(tmp_path), "--start-jd", "2451545",
            "--end-jd", "2451546", "--quiet"]


@pytest.mark.parametrize("mode", ["full", "single"])
def test_full_and_single_keep_complete_inventory_and_original_fit_options(tmp_path, workflow, mode):
    events, _ = workflow
    assert cli.main(arguments(tmp_path, "base", mode) + ["--verify-samples", "17", "--available-domains"]) == 0
    fits = [event for event in events if event[0] == "fit"]
    assert len(fits) == 1
    _, tier, output, bodies, start, end, options = fits[0]
    assert tier == "base" and set(bodies) == cli.BODY_GROUPS["all"]
    assert output == tmp_path / "leb" / "ephemeris_base.leb"
    assert (start, end) == (2451545.0, 2451546.0)
    assert options == dict(skip_aux=False, overwrite=False, verify_samples=17, workers=1,
                           single=mode == "single", available_domains=True, progress=None)
    assert len([event for event in events if event[0] == "source"]) == 1
    assert all(event[0] != "merge" for event in events)


def test_groups_order_inventory_auxiliary_and_merge_preserve_partials(tmp_path, workflow, monkeypatch):
    events, _ = workflow
    output = tmp_path / "custom output.leb"

    def header_start(source, domains, requested):
        assert set(domains) == cli.BODY_GROUPS["all"]
        assert all(domain == (2451545.0, 2451546.0) for domain in domains.values())
        assert requested == 2451545.0
        return 2451544.5

    monkeypatch.setattr(cli, "generation_header_start", header_start)
    assert cli.main(arguments(tmp_path, "BASE", "GROUPS") + ["-o", str(output), "--workers", "2"]) == 0
    fits = [event for event in events if event[0] == "fit"]
    parts = [tmp_path / f"ephemeris_base_{group}.leb" for group in ("planets", "asteroids", "exotics", "analytical")]
    assert [event[2] for event in fits] == parts
    assert [set(event[3]) for event in fits] == [cli.BODY_GROUPS[group] for group in cli.PART_GROUPS]
    assert [event[6]["skip_aux"] for event in fits] == [False, True, True, True]
    assert all(event[6]["header_start"] == 2451544.5 for event in fits)
    assert all(event[4:6] == (2451545.0, 2451546.0) for event in fits)
    assert all(event[6]["workers"] == 2 and not event[6]["single"] for event in fits)
    assert events[-1] == ("merge", tuple(parts), output, {"overwrite": False})
    assert len([event for event in events if event[0] == "source"]) == 1
    assert set().union(*(set(event[3]) for event in fits)) == cli.BODY_GROUPS["all"]
    assert all(path.read_bytes() == b"fitted" for path in parts)
    assert output.read_bytes() == b"merged"


def test_merge_requires_each_named_inventory_without_acquiring_source(tmp_path, workflow):
    events, inventories = workflow
    for group in cli.PART_GROUPS:
        path = tmp_path / f"ephemeris_base_{group}.leb"
        path.write_bytes(b"existing part")
        inventories[path] = sorted(cli.BODY_GROUPS[group])
    command = ["base", "merge", "--output-dir", str(tmp_path), "--quiet"]
    assert cli.main(command) == 0
    assert [event[0] for event in events] == ["merge"]
    events.clear()
    inventories[tmp_path / "ephemeris_base_analytical.leb"].pop()
    with pytest.raises(SystemExit) as failure:
        cli.main(command + ["--force"])
    assert failure.value.code == 1 and events == []
    assert (tmp_path / "ephemeris_base.leb").read_bytes() == b"merged"


@pytest.mark.parametrize("mode", ["planets", "asteroids", "exotics", "analytical", "core", "apogee"])
def test_named_mode_selects_exact_inventory(tmp_path, workflow, mode):
    events, _ = workflow
    assert cli.main(arguments(tmp_path, "base", mode)) == 0
    fit, = [event for event in events if event[0] == "fit"]
    assert set(fit[3]) == cli.BODY_GROUPS[mode]
    assert fit[2] == tmp_path / "leb" / f"ephemeris_base_{mode}.leb"


def test_separate_named_parts_share_header_metadata_without_extra_fit_channels(tmp_path, workflow, monkeypatch):
    events, _ = workflow
    create = cli.create_source

    def source_with_analytic_grid(*args, **kwargs):
        source = create(*args, **kwargs)
        source.piecewise_grids = {21: (2451544.5, 2451550.5, 2.0, 365.25)}
        return source

    monkeypatch.setattr(cli, "create_source", source_with_analytic_grid)
    start = 2451545.123456789
    for mode in cli.PART_GROUPS:
        command = arguments(tmp_path, "base", mode)
        command[command.index("--start-jd") + 1] = str(start)
        assert cli.main(command) == 0
    fits = [event for event in events if event[0] == "fit"]
    assert len({event[6]["header_start"] for event in fits}) == 1
    assert fits[0][6]["header_start"] < start
    assert [set(event[3]) for event in fits] == [cli.BODY_GROUPS[group] for group in cli.PART_GROUPS]
    assert all(event[4:6] == (start, 2451546.0) for event in fits)


def test_case_insensitive_comma_bodies_and_explicit_output_without_data_dir(tmp_path, workflow):
    events, _ = workflow
    output = tmp_path / "two bodies.leb"
    assert cli.main(["base", "body", "mArS, 0", "--kernel", str(tmp_path / "source.bsp"),
                     "-o", str(output), "--start-jd", "2451545", "--end-jd", "2451546", "-q"]) == 0
    fit, = [event for event in events if event[0] == "fit"]
    assert fit[3] == (engine.MARS, 0) and fit[2] == output


def test_all_uses_sequential_tiers_and_stops_on_first_failure(tmp_path, workflow, monkeypatch):
    events, _ = workflow
    original = cli.generate_leb1

    def fail_medium(source, *args, **kwargs):
        if source.tier == "medium":
            raise ValueError("source verification refused")
        return original(source, *args, **kwargs)

    monkeypatch.setattr(cli, "generate_leb1", fail_medium)
    with pytest.raises(SystemExit) as failure:
        cli.main(arguments(tmp_path, "all"))
    assert failure.value.code == 1
    assert [event[1] for event in events if event[0] == "source"] == ["base", "medium"]
    assert (tmp_path / "leb" / "ephemeris_base.leb").exists()
    assert not (tmp_path / "leb" / "ephemeris_medium.leb").exists()
    monkeypatch.setattr(cli, "generate_leb1", original)
    events.clear()
    assert cli.main(arguments(tmp_path, "all") + ["--force", "--skip-aux"]) == 0
    fits = [event for event in events if event[0] == "fit"]
    assert [event[1] for event in fits] == ["base", "medium", "extended"]
    from libephemeris.exotic_bodies import EXOTIC_EXTENDED_IDS, EXOTIC_IDS
    assert set(fits[2][3]) == (cli.BODY_GROUPS["all"] - set(EXOTIC_IDS)) | set(EXOTIC_EXTENDED_IDS)
    assert all(event[6]["skip_aux"] and event[6]["overwrite"] for event in fits)


def test_preflight_refuses_late_existing_part_before_any_source_or_write(tmp_path, workflow):
    events, _ = workflow
    existing = tmp_path / "leb" / "ephemeris_base_analytical.leb"
    existing.parent.mkdir()
    existing.write_bytes(b"preserve")
    with pytest.raises(SystemExit) as failure:
        cli.main(arguments(tmp_path, "base", "groups"))
    assert failure.value.code == 1 and events == []
    assert list(existing.parent.iterdir()) == [existing]
    assert existing.read_bytes() == b"preserve"


def test_source_refusal_with_force_preserves_existing_outputs(tmp_path, workflow, monkeypatch):
    events, _ = workflow
    output = tmp_path / "ephemeris_base.leb"
    output.write_bytes(b"preserve")

    def refuse(*args, **kwargs):
        raise ValueError("missing required local source")

    monkeypatch.setattr(cli, "create_source", refuse)
    with pytest.raises(SystemExit) as failure:
        cli.main(arguments(tmp_path, "base", "full") + ["-o", str(output), "--force"])
    assert failure.value.code == 1 and events == []
    assert output.read_bytes() == b"preserve"


def test_wrong_engine_is_rejected_before_source_or_output(tmp_path, workflow, monkeypatch):
    events, _ = workflow
    monkeypatch.setattr(engine, "__version__", "3.2.1")
    with pytest.raises(SystemExit) as failure:
        cli.main(arguments(tmp_path, "base", "full"))
    assert failure.value.code == 1 and events == []
    assert not (tmp_path / "leb").exists()


@pytest.mark.parametrize("selection, extra", [
    (["base"], []),
    (["base", "full"], ["--group", "planets"]),
    (["base", "full"], ["--bodies", "MARS"]),
    (["base", "full"], ["--single"]),
    (["base", "full"], ["--tier", "medium"]),
    (["base", "single"], ["--workers", "2"]),
    (["base", "full"], ["--workers", "0"]),
    (["base", "body", "mars,4"], []),
    (["base", "body", "mars,"], []),
    (["base", "body", "unknown-body"], []),
    (["base", "body", "40"], []),
    (["base", "body"], []),
    (["all", "full"], []),
    (["all"], ["--kernel", "source.bsp"]),
    (["all"], ["--output", "one-tier.leb"]),
    (["base", "merge"], ["--workers", "2"]),
    (["base", "merge"], ["--kernel", "source.bsp"]),
    (["base", "full"], ["--output", "file.leb", "--output-dir", "directory"]),
])
def test_conflicting_selection_and_arity_never_reaches_source(tmp_path, workflow, selection, extra):
    events, _ = workflow
    with pytest.raises(SystemExit) as failure:
        cli.main(arguments(tmp_path, *selection) + extra)
    assert failure.value.code in (1, 2) and events == []


def test_force_cannot_replace_explicit_source(tmp_path, workflow):
    events, _ = workflow
    source = tmp_path / "source.bsp"
    source.write_bytes(b"immutable input")
    with pytest.raises(SystemExit) as failure:
        cli.main(arguments(tmp_path, "base", "full") + ["--kernel", str(source), "-o", str(source), "--force"])
    assert failure.value.code == 1 and events == []
    assert source.read_bytes() == b"immutable input"


def main_header(start=2451545.0, end=2451546.0):
    return SimpleNamespace(jd_start=start, jd_end=end), [], {}


def test_update_main_refits_one_complete_group_then_replaces_only_that_inventory(tmp_path, workflow, monkeypatch):
    events, _ = workflow
    main = tmp_path / "main.leb"
    main.write_bytes(b"original main")
    monkeypatch.setattr(cli, "read_leb1", lambda path: main_header())
    monkeypatch.setattr(cli, "replace_leb1", lambda main, parts, output, **options:
                        events.append(("replace", main, tuple(parts), output, options)))
    assert cli.main(arguments(tmp_path, "base", "planets") + ["--update-main", str(main), "--force"]) == 0
    fit, = [event for event in events if event[0] == "fit"]
    part = tmp_path / "leb/ephemeris_base_planets.leb"
    assert fit[2] == part and set(fit[3]) == cli.BODY_GROUPS["planets"]
    assert fit[6]["overwrite"] is True and fit[6]["header_start"] == 2451545.0
    assert events[-1] == ("replace", main, (part,), main,
                          {"expected_bodies": set(cli.BODY_GROUPS["planets"]), "overwrite": True})
    assert main.read_bytes() == b"original main"  # The actual publisher has separate native tests.


@pytest.mark.parametrize("selection, force", [
    (["base", "planets"], False),
    (["base", "full"], True),
    (["base", "single"], True),
    (["base", "groups"], True),
    (["base", "merge"], True),
    (["base", "body", "SUN"], True),
    (["base", "core"], True),
    (["all"], True),
])
def test_update_requires_explicit_overwrite_and_one_named_group(tmp_path, workflow, selection, force):
    events, _ = workflow
    command = arguments(tmp_path, *selection) + ["--update-main", str(tmp_path / "absent.leb")]
    with pytest.raises(SystemExit) as failure:
        cli.main(command + (["--force"] if force else []))
    assert failure.value.code == 1 and events == []


@pytest.mark.parametrize("exists", [False, True])
def test_update_rejects_missing_or_invalid_main_before_source_acquisition(tmp_path, workflow, exists):
    events, _ = workflow
    main = tmp_path / "invalid.leb"
    if exists:
        main.write_bytes(b"invalid main")
    with pytest.raises(SystemExit):
        cli.main(arguments(tmp_path, "base", "planets") + ["--update-main", str(main), "--force"])
    assert events == []
    assert not (tmp_path / "leb").exists()


@pytest.mark.parametrize("hardlink", [False, True])
def test_update_part_cannot_replace_main_before_overlay(tmp_path, workflow, monkeypatch, hardlink):
    events, _ = workflow
    main, part = tmp_path / "main.leb", tmp_path / "part.leb"
    main.write_bytes(b"original main")
    if hardlink:
        part.hardlink_to(main)
    else:
        part = main
    monkeypatch.setattr(cli, "read_leb1", lambda path: main_header())
    with pytest.raises(SystemExit):
        cli.main(arguments(tmp_path, "base", "planets") +
                 ["--update-main", str(main), "--force", "-o", str(part)])
    assert events == [] and main.read_bytes() == b"original main"


def test_update_refuses_global_domain_mismatch_before_any_fit(tmp_path, workflow, monkeypatch):
    events, _ = workflow
    main = tmp_path / "main.leb"
    main.write_bytes(b"original main")
    monkeypatch.setattr(cli, "read_leb1", lambda path: main_header(end=2451547.0))
    with pytest.raises(SystemExit):
        cli.main(arguments(tmp_path, "base", "planets") + ["--update-main", str(main), "--force"])
    assert all(event[0] != "fit" for event in events)
    assert main.read_bytes() == b"original main"


def test_update_failure_retains_verified_part_and_original_main(tmp_path, workflow, monkeypatch):
    events, _ = workflow
    main = tmp_path / "main.leb"
    main.write_bytes(b"original main")
    monkeypatch.setattr(cli, "read_leb1", lambda path: main_header())

    def refuse(*args, **options):
        raise ValueError("replacement input changed")

    monkeypatch.setattr(cli, "replace_leb1", refuse)
    with pytest.raises(SystemExit):
        cli.main(arguments(tmp_path, "base", "planets") + ["--update-main", str(main), "--force"])
    assert len([event for event in events if event[0] == "fit"]) == 1
    assert (tmp_path / "leb/ephemeris_base_planets.leb").read_bytes() == b"fitted"
    assert main.read_bytes() == b"original main"
