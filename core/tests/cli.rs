// SPDX-License-Identifier: AGPL-3.0
//! End-to-end checks for the shipped `libephemeris` binary (the CLI
//! inspection and provisioning commands).
//!
//! The hermetic cases run everywhere (version spelling, usage errors,
//! missing-file and corrupt-file refusals, all-missing status over an empty
//! directory, invalid download options, config
//! reports, a no-op ASSIST download that selects no files and touches no
//! socket, a defaults `init --non-interactive` writing to a scratch
//! path, and the interactive `init` wizard driven from piped stdin or a
//! closed stdin — the post-wizard report only checks files read-only and
//! prints the served download spellings, so no socket is opened). The `local_*`
//! cases additionally exercise a real installed distribution read-only when
//! one is present and skip otherwise, following the gated `local_leb_dir`
//! pattern in `src/leb.rs`.
//!
//! Real ASSIST fetching is never exercised here: the download pipeline is
//! covered by hermetic loopback tests in `src/download.rs`, and live remote
//! calls stay out of the suite. These cases pin the argument contract and
//! the offline `config` report only.

use std::path::PathBuf;
use std::process::Command;

fn binary() -> Command {
    let exe = env!("CARGO_BIN_EXE_libephemeris");
    let mut cmd = Command::new(exe);
    // Never let the developer's installed distribution leak into the
    // hermetic cases.
    cmd.env_remove("LIBEPHEMERIS_DATA_DIR");
    cmd
}

fn local_leb_dir() -> Option<PathBuf> {
    if let Ok(dir) = std::env::var("LIBEPHEMERIS_DATA_DIR") {
        let candidate = PathBuf::from(dir).join("leb");
        if candidate.is_dir() {
            return Some(candidate);
        }
    }
    if let Ok(home) = std::env::var("HOME") {
        let candidate = PathBuf::from(home).join(".libephemeris").join("leb");
        if candidate.is_dir() {
            return Some(candidate);
        }
    }
    None
}

#[test]
fn version_reports_the_crate_version() {
    let out = binary().arg("--version").output().expect("binary runs");
    assert!(out.status.success());
    let text = String::from_utf8_lossy(&out.stdout);
    assert!(text.starts_with("libephemeris, version "));
    assert!(text.contains(env!("CARGO_PKG_VERSION")));
}

#[test]
fn unknown_command_is_a_usage_error() {
    let out = binary().arg("bogus").output().expect("binary runs");
    assert_eq!(out.status.code(), Some(2));
}

#[test]
fn complete_bootstrap_requires_confirmation_before_acquiring_or_creating_data() {
    let directory = std::env::temp_dir().join(format!(
        "libephemeris-cli-bootstrap-declined-{}",
        std::process::id()
    ));
    assert!(!directory.exists());
    let result = binary()
        .args(["download", "all", "--data-dir"])
        .arg(&directory)
        .stdin(std::process::Stdio::null())
        .output()
        .unwrap();
    assert_eq!(result.status.code(), Some(1));
    let prompt = String::from_utf8_lossy(&result.stdout);
    assert!(prompt.contains("planet-center sources and ASSIST"));
    assert!(prompt.contains("[y/N]"));
    assert!(!directory.exists());
}

#[test]
fn command_help_never_acquires_data() {
    for arguments in [
        vec!["download", "--help"],
        vec!["download", "medium", "--help"],
        vec!["download", "all", "--help"],
        vec!["download", "planet-centers-sources", "--help"],
        vec!["status", "--help"],
        vec!["init", "--help"],
        vec!["generate", "planet-centers-base", "--help"],
        vec!["generate", "planet-centers-all", "--help"],
    ] {
        let output = binary().args(&arguments).output().unwrap();
        assert!(
            output.status.success(),
            "{arguments:?}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(String::from_utf8_lossy(&output.stdout).contains("download base|medium|extended"));
    }
    let output = binary()
        .args(["download", "unknown", "--help"])
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(2));
}

#[test]
fn detailed_status_uses_the_configured_tier_and_lists_files_without_installation() {
    let directory = std::env::temp_dir().join(format!(
        "libephemeris-cli-status-detail-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(directory.join("spk")).unwrap();
    std::fs::write(directory.join("spk/local.bsp"), b"presence-only metadata").unwrap();
    let output = binary()
        .env("LIBEPHEMERIS_PRECISION", "base")
        .args(["status", "-vv", "--data-dir"])
        .arg(&directory)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = String::from_utf8_lossy(&output.stdout);
    assert!(text.contains("tier: base"));
    assert!(text.contains(&directory.join("spk/local.bsp").display().to_string()));
    assert!(text.contains("de440s.bsp"));
    assert!(!directory.join("leb").exists());
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn download_needs_a_served_target() {
    for args in [
        vec!["download"],
        vec!["download", "leb2-huge"],
        vec!["download", "leb2-base", "--bogus"],
        vec!["download", "iers", "--bogus"],
        vec!["download", "iers", "--data-dir"],
        vec!["download", "assist", "--bogus"],
        vec!["download", "assist", "--target-dir"],
        vec!["download", "planet-centers-sources", "--tier"],
        vec!["download", "planet-centers-sources", "--tier", "unknown"],
        vec!["init", "--non-interactive", "--bogus"],
        vec!["init", "--output"],
    ] {
        let out = binary().args(&args).output().expect("binary runs");
        assert_eq!(out.status.code(), Some(2), "args: {args:?}");
    }
}

#[test]
fn planet_center_source_selection_errors_do_not_create_or_change_a_cache() {
    let directory = std::env::temp_dir().join(format!(
        "libephemeris-cli-center-source-usage-{}",
        std::process::id()
    ));
    std::fs::create_dir(&directory).unwrap();
    let sentinel = directory.join("source-preserved");
    std::fs::write(&sentinel, b"existing local source").unwrap();
    for arguments in [
        vec!["--tier", "bogus"],
        vec!["--tier", "base", "--tier", "medium"],
        vec!["--bogus"],
    ] {
        let result = binary()
            .args(["download", "planet-centers-sources", "--data-dir"])
            .arg(&directory)
            .args(arguments)
            .output()
            .unwrap();
        assert_eq!(result.status.code(), Some(2));
        assert_eq!(std::fs::read(&sentinel).unwrap(), b"existing local source");
        assert!(!directory.join("tmp").exists());
    }
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn provisioning_rejects_bad_options_before_any_network_work() {
    for target in ["medium", "base", "extended", "leb-base", "auto", "all"] {
        let args = ["download", target, "--bogus"];
        let out = binary().args(args).output().expect("binary runs");
        assert_eq!(out.status.code(), Some(2), "args: {args:?}");
        let text = String::from_utf8_lossy(&out.stderr);
        assert!(text.contains("--bogus"), "args: {args:?}, stderr: {text}");
    }
}

#[cfg(feature = "bundled-data")]
#[test]
fn legacy_core_download_installs_into_only_the_explicit_directory() {
    let dir =
        std::env::temp_dir().join(format!("libephemeris-cli-core-e2e-{}", std::process::id()));
    let out = binary()
        .args(["download:leb:base", "--no-progress", "-q", "--data-dir"])
        .arg(&dir)
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let installed = dir.join("leb/base_core.leb2");
    assert_eq!(
        libephemeris_core::download::file_sha256(&installed).unwrap(),
        libephemeris_core::download::asset_by_name("base_core.leb2")
            .unwrap()
            .sha256()
    );
    assert!(out.stdout.is_empty());
    assert!(!dir.join("spk").exists());
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn init_non_interactive_writes_defaults_and_respects_force() {
    let dir = std::env::temp_dir().join("libephemeris-cli-init");
    std::fs::create_dir_all(&dir).expect("scratch dir");
    let target = dir.join("libephemeris-config.toml");
    std::fs::remove_file(&target).ok();
    let out = binary()
        .arg("init")
        .arg("--non-interactive")
        .arg("-o")
        .arg(&target)
        .output()
        .expect("binary runs");
    assert!(out.status.success());
    let written = std::fs::read_to_string(&target).expect("config written");
    assert!(written.contains("precision = \"medium\""));
    assert!(written.contains("deltat_model = \"smh2016\""));
    // A second write without --force refuses; with --force it succeeds.
    let out = binary()
        .arg("init")
        .arg("--non-interactive")
        .arg("-o")
        .arg(&target)
        .output()
        .expect("binary runs");
    assert_eq!(out.status.code(), Some(1));
    let out = binary()
        .arg("init")
        .arg("--non-interactive")
        .arg("-o")
        .arg(&target)
        .arg("--force")
        .output()
        .expect("binary runs");
    assert!(out.status.success());
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn init_wizard_answers_defaults_from_piped_stdin() {
    use std::process::Stdio;
    let dir = std::env::temp_dir().join("libephemeris-cli-wizard");
    std::fs::create_dir_all(&dir).expect("scratch dir");
    let target = dir.join("libephemeris-config.toml");
    std::fs::remove_file(&target).ok();
    let input = dir.join("answers.txt");
    std::fs::write(&input, "yes\n").expect("answers file");
    let stdin = std::fs::File::open(&input).expect("answers readable");
    let out = binary()
        .arg("init")
        .arg("-o")
        .arg(&target)
        .stdin(Stdio::from(stdin))
        .output()
        .expect("binary runs");
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let written = std::fs::read_to_string(&target).expect("config written");
    assert!(written.contains("precision = \"medium\""));
    assert!(written.contains("deltat_model = \"smh2016\""));
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("Files needed for your configuration"));
    assert!(stdout.contains("download leb2-medium"));
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn init_wizard_aborts_cleanly_on_closed_stdin() {
    use std::process::Stdio;
    let dir = std::env::temp_dir().join("libephemeris-cli-wizard-eof");
    std::fs::create_dir_all(&dir).expect("scratch dir");
    let target = dir.join("libephemeris-config.toml");
    std::fs::remove_file(&target).ok();
    let out = binary()
        .arg("init")
        .arg("-o")
        .arg(&target)
        .stdin(Stdio::null())
        .output()
        .expect("binary runs");
    assert_eq!(out.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&out.stdout).contains("Aborted."));
    assert!(!target.exists(), "aborted wizard writes nothing");
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn download_assist_with_no_files_selected_is_a_hermetic_noop() {
    // Both files deselected: the engine loop requests nothing, so no
    // socket is opened and the run succeeds offline.
    let dir = std::env::temp_dir().join("libephemeris-cli-assist-noop");
    let out = binary()
        .arg("download")
        .arg("assist")
        .arg("--no-planets")
        .arg("--no-asteroids")
        .arg("--target-dir")
        .arg(&dir)
        .output()
        .expect("binary runs");
    assert!(out.status.success());
    let text = String::from_utf8_lossy(&out.stdout);
    assert!(text.contains("nothing requested"));
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn config_reports_data_dir_and_env_knobs() {
    let out = binary().arg("config").output().expect("binary runs");
    assert!(out.status.success());
    let text = String::from_utf8_lossy(&out.stdout);
    assert!(text.contains("data dir:"));
    assert!(text.contains("LIBEPHEMERIS_DATA_DIR"));
    assert!(text.contains("LIBEPHEMERIS_NETWORK_POLICY"));
    assert!(text.contains("tier: medium"));
}

#[test]
fn config_json_reports_nulls_for_unset_knobs() {
    let out = binary()
        .arg("config")
        .arg("--json")
        .output()
        .expect("binary runs");
    assert!(out.status.success());
    let text = String::from_utf8_lossy(&out.stdout);
    assert!(text.contains("\"data_dir\":"));
    assert!(text.contains("\"delta_t\": \"smh2016\""));
}

#[test]
fn leb_info_refuses_a_missing_file() {
    let out = binary()
        .args(["leb", "info", "/nonexistent-dir/missing.leb2"])
        .output()
        .expect("binary runs");
    assert_eq!(out.status.code(), Some(1));
    assert!(!String::from_utf8_lossy(&out.stderr).is_empty());
}

#[test]
fn leb_info_refuses_a_corrupt_file() {
    let dir = std::env::temp_dir().join("libephemeris-cli-test");
    std::fs::create_dir_all(&dir).expect("scratch dir");
    let corrupt = dir.join("corrupt.leb2");
    std::fs::write(&corrupt, b"not a leb file").expect("scratch file");
    let out = binary()
        .arg("leb")
        .arg("info")
        .arg(&corrupt)
        .output()
        .expect("binary runs");
    assert_eq!(out.status.code(), Some(1));
    std::fs::remove_file(&corrupt).ok();
}

#[test]
fn status_over_an_empty_dir_reports_every_asset_missing() {
    let dir = std::env::temp_dir().join("libephemeris-cli-empty");
    std::fs::create_dir_all(dir.join("leb")).expect("scratch dir");
    let out = binary()
        .arg("status")
        .arg("--tier")
        .arg("base")
        .arg("--data-dir")
        .arg(&dir)
        .arg("--json")
        .output()
        .expect("binary runs");
    // Status describes availability successfully; absent assets are metadata.
    assert!(out.status.success());
    let snapshot: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(snapshot["precision_tier"], "base");
    assert_eq!(snapshot["leb2_files"]["base_core.leb2"]["exists"], false);
    assert_eq!(snapshot["de_kernels"]["de440s.bsp"]["exists"], false);
    assert_eq!(snapshot["spk_cache"]["file_count"], 0);
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn leb_and_leb2_usage_lattice() {
    for args in [
        vec!["leb"],
        vec!["leb", "bogus"],
        vec!["leb", "generate"],
        vec!["leb", "generate", "bogus", "groups"],
        vec!["leb", "generate", "medium"],
        vec!["leb", "generate", "medium", "bogus"],
        vec!["leb", "generate", "medium", "body"],
        vec!["leb", "verify"],
        vec!["leb", "verify", "a", "b"],
        vec!["leb", "verify", "a", "--data-dir"],
        vec!["leb2"],
        vec!["leb2", "bogus"],
        vec!["leb2", "convert"],
        vec!["leb2", "convert", "medium", "extra"],
        vec!["leb2", "convert", "bogus"],
        vec!["leb2", "convert", "medium-bogus"],
        vec!["leb2", "verify"],
        vec!["leb2", "verify", "a", "--reference"],
        vec!["leb2", "verify", "a", "--samples"],
        vec!["leb2", "verify", "a", "--samples", "0"],
        vec!["leb2", "verify", "a", "--samples", "5000"],
        vec!["leb2", "verify", "a", "--samples", "bogus"],
    ] {
        let out = binary().args(&args).output().expect("binary runs");
        assert_eq!(out.status.code(), Some(2), "args: {args:?}");
    }
}

#[test]
fn leb_generation_requires_the_explicit_optional_application() {
    for args in [
        vec!["leb", "generate", "medium", "groups"],
        vec!["leb", "generate", "base", "full"],
        vec!["leb", "generate", "medium", "body", "moon"],
    ] {
        let out = binary()
            .env("PATH", "/no-such-libephemeris-tools")
            .args(&args)
            .output()
            .expect("binary runs");
        assert_eq!(out.status.code(), Some(1), "args: {args:?}");
        let text = String::from_utf8_lossy(&out.stderr);
        assert!(
            text.contains("cannot start offline LEB generator"),
            "{text}"
        );
        assert!(text.contains("libephemeris-leb-tools 4.0.0b1"), "{text}");
        assert!(out.stdout.is_empty());
    }
}

#[test]
fn multi_epoch_generation_reports_missing_inputs_without_writing_a_product() {
    let directory = std::env::temp_dir().join(format!(
        "libephemeris-elements-empty-{}",
        std::process::id()
    ));
    std::fs::create_dir(&directory).unwrap();
    let output = directory.join("elements.py");
    let out = binary()
        .args(["generate", "keplerian-elements"])
        .arg("--spk-dir")
        .arg(&directory)
        .args(["--body", "Ceres", "--output"])
        .arg(&output)
        .output()
        .expect("binary runs");
    assert_eq!(out.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&out.stderr).contains("no element rows generated"));
    assert!(!output.exists());
    let dry = binary()
        .args(["generate", "keplerian-elements", "--dry-run", "--quiet"])
        .arg("--spk-dir")
        .arg(&directory)
        .args(["--body", "Ceres"])
        .output()
        .unwrap();
    assert!(dry.status.success());
    let report: serde_json::Value = serde_json::from_slice(&dry.stdout).unwrap();
    assert_eq!(report["epoch_scale"], "TDB");
    assert_eq!(report["requested_tdb_jd"].as_array().unwrap().len(), 81);
    assert_eq!(report["omissions"][0]["body"], "Ceres");
    assert!(dry.stderr.is_empty());
    std::fs::remove_dir(&directory).unwrap();
}

#[cfg(unix)]
#[test]
fn native_generation_passes_literal_arguments_and_preserves_child_status() {
    use std::os::unix::fs::PermissionsExt;
    let directory =
        std::env::temp_dir().join(format!("libephemeris-cli-generator-{}", std::process::id()));
    std::fs::create_dir(&directory).unwrap();
    let helper = directory.join("libephemeris-generate");
    std::fs::write(&helper, "#!/bin/sh\nprintf '%s\\n' \"$@\"\nexit 23\n").unwrap();
    std::fs::set_permissions(&helper, std::fs::Permissions::from_mode(0o700)).unwrap();
    let out = binary()
        .env("PATH", &directory)
        .args([
            "leb",
            "generate",
            "medium",
            "body",
            "Sun,Moon,14",
            "--output",
            "literal ; $(false).leb",
        ])
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(23));
    assert_eq!(
        String::from_utf8(out.stdout).unwrap(),
        "medium\nbody\nSun,Moon,14\n--output\nliteral ; $(false).leb\n"
    );
    assert!(out.stderr.is_empty());
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn leb2_conversion_reports_missing_input_as_an_operational_error() {
    // Restored conversion is executable. Always provide an explicit absent
    // source, so this refusal test cannot acquire the developer's real data.
    for target in ["medium", "medium-core"] {
        let output = binary()
            .args([
                "leb2",
                "convert",
                target,
                "--input",
                "/nonexistent-dir/source.leb",
                "--output-dir",
                "/nonexistent-dir/products",
            ])
            .output()
            .expect("binary runs");
        assert_eq!(output.status.code(), Some(1));
        assert!(output.stdout.is_empty());
        assert!(String::from_utf8_lossy(&output.stderr).contains("source.leb"));
    }
}

#[test]
fn leb_verify_refuses_missing_and_corrupt_files() {
    let out = binary()
        .args(["leb", "verify", "/nonexistent-dir/missing.leb2"])
        .output()
        .expect("binary runs");
    assert_eq!(out.status.code(), Some(1));
    let dir = std::env::temp_dir().join("libephemeris-cli-verify");
    std::fs::create_dir_all(&dir).expect("scratch dir");
    let corrupt = dir.join("corrupt.leb2");
    std::fs::write(&corrupt, b"not a leb file").expect("scratch file");
    for args in [
        vec!["leb", "verify", corrupt.to_string_lossy().as_ref()],
        vec!["leb2", "verify", corrupt.to_string_lossy().as_ref()],
        vec![
            "leb2",
            "verify",
            corrupt.to_string_lossy().as_ref(),
            "--reference",
            corrupt.to_string_lossy().as_ref(),
        ],
    ] {
        let out = binary().args(&args).output().expect("binary runs");
        assert_eq!(out.status.code(), Some(1), "args: {args:?}");
        assert!(
            !String::from_utf8_lossy(&out.stderr).is_empty(),
            "args: {args:?}"
        );
    }
    let out = binary()
        .args(["leb2", "verify", "/nonexistent-dir/missing.leb2"])
        .output()
        .expect("binary runs");
    assert_eq!(out.status.code(), Some(1));
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn local_leb_verify_passes_over_the_installed_distribution() {
    let Some(dir) = local_leb_dir() else {
        eprintln!("SKIP: no local LEB distribution installed");
        return;
    };
    let entries: Vec<PathBuf> = std::fs::read_dir(&dir)
        .expect("leb dir is readable")
        .filter_map(|entry| entry.ok())
        .map(|entry| entry.path())
        .filter(|path| path.extension().is_some_and(|ext| ext == "leb2"))
        .collect();
    assert!(!entries.is_empty(), "local leb dir holds .leb2 files");
    let target = entries
        .iter()
        .find(|path| {
            path.file_name()
                .is_some_and(|name| name.to_string_lossy().ends_with("core.leb2"))
        })
        .unwrap_or(&entries[0]);
    for args in [
        vec!["leb", "verify", target.to_string_lossy().as_ref()],
        vec!["leb", "verify", target.to_string_lossy().as_ref(), "--json"],
        vec![
            "leb2",
            "verify",
            target.to_string_lossy().as_ref(),
            "--samples",
            "8",
        ],
        vec![
            "leb2",
            "verify",
            target.to_string_lossy().as_ref(),
            "--samples",
            "8",
            "--json",
        ],
    ] {
        let out = binary().args(&args).output().expect("binary runs");
        assert!(
            out.status.success(),
            "args: {args:?}, stderr: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let text = String::from_utf8_lossy(&out.stdout);
        assert!(
            text.contains("PASS") || text.contains("\"passed\": true"),
            "args: {args:?}, stdout: {text}"
        );
    }
}

#[test]
fn local_leb_info_and_status_agree_with_the_installed_distribution() {
    let Some(dir) = local_leb_dir() else {
        eprintln!("SKIP: no local LEB distribution installed");
        return;
    };
    let entries: Vec<PathBuf> = std::fs::read_dir(&dir)
        .expect("leb dir is readable")
        .filter_map(|entry| entry.ok())
        .map(|entry| entry.path())
        .filter(|path| path.extension().is_some_and(|ext| ext == "leb2"))
        .collect();
    assert!(!entries.is_empty(), "local leb dir holds .leb2 files");
    let target = entries.iter().find(|path| {
        path.file_name()
            .is_some_and(|name| name.to_string_lossy().ends_with("core.leb2"))
    });
    let target = target.unwrap_or(&entries[0]);
    let out = binary()
        .arg("leb")
        .arg("info")
        .arg(target)
        .output()
        .expect("binary runs");
    assert!(out.status.success());
    let text = String::from_utf8_lossy(&out.stdout);
    assert!(text.contains("bodies)"));
    assert!(text.contains("JD start"));
}
