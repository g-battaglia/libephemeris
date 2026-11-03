// SPDX-License-Identifier: AGPL-3.0
//! Shared command-line application for the native binary and Python console.
//!
//! Commands inspect/verify LEB images, report configured data, provision
//! explicitly requested data, run the setup wizard and inspect SPK element
//! grids. Astronomical decoding/evaluation delegates to the existing native
//! readers and kernels. Output goes to standard streams; returned exit codes
//! are 0 for success, 1 for operational failure and 2 for usage errors.
//!
//! The caller provides arguments without a program name. No command terminates
//! the embedding process. Downloads occur only through explicit provisioning
//! commands. Configuration uses the core environment/dotenv/TOML resolution.

use std::path::{Path, PathBuf};

use crate::leb::LebTier;

#[path = "cli/auxiliary_downloads.rs"]
mod auxiliary_downloads;
#[path = "cli/element_sources.rs"]
mod element_sources;
#[path = "cli/inspection.rs"]
mod inspection;
#[path = "cli/keplerian.rs"]
mod keplerian;
#[path = "cli/keplerian_availability.rs"]
mod keplerian_availability;
#[path = "cli/leb_conversion.rs"]
mod leb_conversion;
#[path = "cli/leb_generation.rs"]
mod leb_generation;
#[path = "cli/leb_verification.rs"]
mod leb_verification;
#[path = "cli/planet_center_sources.rs"]
mod planet_center_sources;
#[path = "cli/planet_centers.rs"]
mod planet_centers;
#[path = "cli/provision.rs"]
mod provision;
#[path = "cli/wizard.rs"]
mod wizard;

use keplerian_availability::{parse_dry_run_args, run_generate_keplerian_dry_run};

const EXIT_OK: i32 = 0;
const EXIT_FAILED: i32 = 1;
const EXIT_USAGE: i32 = 2;

const USAGE: &str = "\
usage:
  libephemeris [--version | -V | version]
  libephemeris leb info <path> [--json]
  libephemeris leb generate <tier> <groups|full|single|body NAME[,NAME...]|planets|asteroids|exotics|analytical|merge> [offline source/output options]
  libephemeris leb generate all --data-dir DIR [--output-dir DIR] [offline fitting options]
  libephemeris leb verify <tier|path> [--data-dir DIR] [--json]
  libephemeris leb2 convert <tier|tier-group> [--input PATH] [--output PATH | --output-dir DIR] [--data-dir DIR] [--max-chunk-bytes N] [--force] [--json]
  libephemeris leb2 verify <path> [--reference PATH] [--samples N] [--json]
  libephemeris status [--tier base|medium|extended] [--data-dir DIR] [--json] [-v|-vv]
  libephemeris download leb2-base|leb2-medium|leb2-extended [--force] [--quiet] [--data-dir DIR]
  libephemeris download base|medium|extended|leb-base|leb-medium|leb-extended|auto|all [--force] [--quiet] [--no-progress] [--data-dir DIR] [--yes]
  libephemeris download assist [--no-planets] [--no-asteroids] [--target-dir DIR] [--force] [--quiet] [--no-progress] [--data-dir DIR]
  libephemeris download iers [--force] [--quiet] [--data-dir DIR]
  libephemeris download planet-centers-sources [--tier base|medium|extended] [--force] [--quiet] [--data-dir DIR]
   libephemeris config [--json]
   libephemeris init [-o PATH] [--force] [--non-interactive]
   libephemeris generate planet-centers-base|planet-centers-medium|planet-centers-extended|planet-centers-all|planet-centers-spk [--cache-dir DIR] [--output-dir DIR] [--data-dir DIR] [--force] [--json] [--require-full-tier]
   libephemeris generate keplerian-elements [--spk-dir DIR]... [--kernel BODY=PATH]... [--body NAME]... [--start YEAR] [--end YEAR] [--spacing YEARS] [--output PATH] [--format python|json] [--require-complete] [--dry-run] [--quiet] [--force]
   libephemeris generate keplerian-dry-run [--spk-dir DIR]... [--data-dir DIR] [--start YEAR] [--end YEAR] [--spacing YEARS] [--body NAME]...";

/// Parsed command line (pure; unit-tested without touching the file system).
#[derive(Debug, Clone, PartialEq, Eq)]
enum Command {
    GenerateKeplerian(keplerian::Request),
    GeneratePlanetCenters(planet_centers::Request),
    Provision(provision::Request),
    Leb2Convert(leb_conversion::Request),
    LebGenerate(leb_generation::Request),
    Version,
    LebInfo {
        path: String,
        json: bool,
    },
    Status {
        tier: Option<LebTier>,
        data_dir: Option<String>,
        json: bool,
        verbose: i32,
    },
    DownloadAssist {
        target_dir: Option<String>,
        data_dir: Option<String>,
        planets: bool,
        asteroids: bool,
        force: bool,
        quiet: bool,
        no_progress: bool,
    },
    DownloadIers {
        data_dir: Option<String>,
        force: bool,
        quiet: bool,
    },
    Config {
        json: bool,
    },
    Init {
        output: String,
        force: bool,
        non_interactive: bool,
    },
    GenerateKeplerianDryRun {
        spk_dirs: Vec<String>,
        data_dir: Option<String>,
        start: i32,
        end: i32,
        spacing: i32,
        bodies: Vec<String>,
    },
    LebVerify {
        target: String,
        data_dir: Option<String>,
        json: bool,
    },
    Leb2Verify {
        path: String,
        reference: Option<String>,
        samples: usize,
        json: bool,
    },
}

/// Parse `argv` (without the program name) into a [`Command`].
///
/// Usage errors carry a human-readable message; the caller prints [`USAGE`]
/// beside it and exits with [`EXIT_USAGE`].
fn parse_args(argv: &[String]) -> Result<Command, String> {
    if let Some(first) = argv.first() {
        let alias = match first.as_str() {
            "download:base" => Some("base"),
            "download:medium" => Some("medium"),
            "download:extended" => Some("extended"),
            "download:leb:base" => Some("leb-base"),
            "download:leb:medium" => Some("leb-medium"),
            "download:leb:extended" => Some("leb-extended"),
            "download:assist" => Some("assist"),
            _ => None,
        };
        if let Some(target) = alias {
            let canonical = [vec!["download".into(), target.into()], argv[1..].to_vec()].concat();
            return parse_args(&canonical);
        }
    }
    let mut args = argv.iter().map(String::as_str).peekable();
    match args.next() {
        None => Err("no command given".to_string()),
        Some("--version" | "-V" | "version") => {
            if args.next().is_some() {
                return Err("version takes no arguments".to_string());
            }
            Ok(Command::Version)
        }
        Some("leb") => match args.next() {
            Some("info") => {
                let mut path: Option<String> = None;
                let mut json = false;
                for arg in args {
                    if arg == "--json" {
                        json = true;
                    } else if path.is_none() {
                        path = Some(arg.to_string());
                    } else {
                        return Err(format!("unexpected argument '{arg}'"));
                    }
                }
                match path {
                    Some(path) => Ok(Command::LebInfo { path, json }),
                    None => Err("leb info needs a file path".to_string()),
                }
            }
            Some("generate") => leb_generation::parse(args).map(Command::LebGenerate),
            Some("verify") => leb_verification::parse_leb_verify_args(&mut args),
            Some(other) => Err(format!("unknown leb subcommand '{other}'")),
            None => Err("leb needs a subcommand".to_string()),
        },
        Some("leb2") => parse_leb2_args(&mut args),
        Some("status") => {
            let mut tier = None;
            let mut data_dir: Option<String> = None;
            let mut json = false;
            let mut verbose = 0;
            while let Some(arg) = args.next() {
                if arg == "--json" {
                    json = true;
                } else if arg == "-v" || arg == "--verbose" {
                    verbose += 1;
                } else if arg == "-vv" {
                    verbose += 2;
                } else if arg == "--tier" {
                    match args.next() {
                        Some(token) => match LebTier::parse(token) {
                            Some(parsed) => tier = Some(parsed),
                            None => {
                                return Err(format!(
                                    "unknown tier '{token}' (base|medium|extended)"
                                ));
                            }
                        },
                        None => return Err("--tier needs a value".to_string()),
                    }
                } else if arg == "--data-dir" {
                    match args.next() {
                        Some(dir) => data_dir = Some(dir.to_string()),
                        None => return Err("--data-dir needs a value".to_string()),
                    }
                } else {
                    return Err(format!("unexpected argument '{arg}'"));
                }
            }
            Ok(Command::Status {
                tier,
                data_dir,
                json,
                verbose,
            })
        }
        Some("download") => parse_download_args(&mut args),
        Some("config") => {
            let mut json = false;
            for arg in args {
                if arg == "--json" {
                    json = true;
                } else {
                    return Err(format!("unexpected argument '{arg}'"));
                }
            }
            Ok(Command::Config { json })
        }
        Some("init") => wizard::parse_init_args(&mut args),
        Some("generate") => parse_generate_args(&mut args),
        Some(other) => Err(format!("unknown command '{other}'")),
    }
}

/// Parse the words after `generate` (pure; unit-tested).
///
/// Planet-center assembly and multi-epoch elements use local NASA/JPL sources.
fn parse_generate_args<'a>(
    args: &mut std::iter::Peekable<impl Iterator<Item = &'a str>>,
) -> Result<Command, String> {
    let subcommand = match args.next() {
        Some(subcommand) => subcommand,
        None => return Err("generate needs a subcommand".to_string()),
    };
    match subcommand {
        "planet-centers-base"
        | "planet-centers-medium"
        | "planet-centers-extended"
        | "planet-centers-all"
        | "planet-centers-spk" => {
            planet_centers::parse(subcommand, args).map(Command::GeneratePlanetCenters)
        }
        "keplerian-elements" => keplerian::parse(args).map(Command::GenerateKeplerian),
        "keplerian-dry-run" => parse_dry_run_args(args),
        other => Err(format!("unknown generate subcommand '{other}'")),
    }
}

/// Parse the `download` subcommand.
/// Supported targets include tiered assets, planet-center sources, ASSIST and IERS. Other spellings
/// return a usage error with the relevant provisioning command or limitation;
/// no unsupported target is reported as a successful download.
fn parse_download_args<'a>(
    args: &mut std::iter::Peekable<impl Iterator<Item = &'a str>>,
) -> Result<Command, String> {
    let target = match args.next() {
        Some(target) => target,
        None => {
            return Err(
                "download needs a target (leb2-base|leb2-medium|leb2-extended|assist|iers|planet-centers-sources)"
                    .to_string(),
            );
        }
    };
    if target == "assist" {
        return auxiliary_downloads::parse_download_assist_args(args);
    }
    if target == "iers" {
        return auxiliary_downloads::parse_download_iers_args(args);
    }
    if let Some(request) = provision::parse(target, args)? {
        return Ok(Command::Provision(request));
    }
    Err(format!("unknown download target '{target}'"))
}

/// Parse native lossless conversion and verification requests.
fn parse_leb2_args<'a>(
    args: &mut std::iter::Peekable<impl Iterator<Item = &'a str>>,
) -> Result<Command, String> {
    match args.next() {
        Some("convert") => leb_conversion::parse(args).map(Command::Leb2Convert),
        Some("verify") => leb_verification::parse_leb2_verify_args(args),
        Some(other) => Err(format!("unknown leb2 subcommand '{other}'")),
        None => Err("leb2 needs a subcommand".to_string()),
    }
}

/// Resolve an explicit CLI data path above the shared core configuration.
fn resolve_data_dir(
    flag: Option<&str>,
    sources: crate::config::ConfigurationSources,
) -> Result<PathBuf, String> {
    match flag {
        Some(path) => Ok(PathBuf::from(path)),
        None => crate::config::Runtime::with_sources(sources)
            .data_directory()
            .map_err(|error| error.to_string()),
    }
}

/// Minimal JSON string escaping for the machine-readable outputs.
fn json_escape(text: &str) -> String {
    let mut out = String::with_capacity(text.len() + 2);
    for ch in text.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            ch if (ch as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", ch as u32)),
            ch => out.push(ch),
        }
    }
    out
}

/// Format a Julian Day for JSON output with full round-trip precision.
fn json_f64(value: f64) -> String {
    format!("{value:?}")
}

/// Run the shared native command-line application.
///
/// `args` contains command/options without a program name. Commands use the
/// configured native data paths and print to standard output/error. Explicit
/// provisioning commands may access the network and filesystem; `init` may
/// read standard input. Returns 0 on success, 1 on operational failure or 2
/// for invalid arguments. Offline generation preserves its child application's
/// normal exit code. It never calls process::exit.
pub fn run<I, S>(args: I) -> i32
where
    I: IntoIterator<Item = S>,
    S: AsRef<str>,
{
    run_with_generator(args, None)
}

/// Run the CLI with an explicit standalone offline generator executable.
///
/// `generator` selects that exact path; `None` searches PATH for
/// `libephemeris-generate`. Arguments, standard streams and the working
/// directory are passed to the child without a shell. Only LEB generation
/// uses this process boundary. No dependency is installed and no checkout is
/// searched. Other commands and exit semantics are the same as [`run`].
/// An unavailable executable or abnormal child termination returns 1.
pub fn run_with_generator<I, S>(args: I, generator: Option<&Path>) -> i32
where
    I: IntoIterator<Item = S>,
    S: AsRef<str>,
{
    let argv: Vec<String> = args
        .into_iter()
        .map(|word| word.as_ref().to_owned())
        .collect();
    let help = argv
        .last()
        .is_some_and(|word| word == "--help" || word == "-h");
    let known_help_target = match &argv[..argv.len().saturating_sub(1)] {
        [] => true,
        [command] => matches!(
            command.as_str(),
            "download" | "status" | "config" | "init" | "leb" | "leb2" | "generate"
        ),
        [command, target] if command == "download" => {
            matches!(target.as_str(), "assist" | "iers")
                || provision::parse(target, &mut std::iter::empty())
                    .is_ok_and(|request| request.is_some())
        }
        [command, target] if command == "leb" => {
            matches!(target.as_str(), "info" | "verify" | "generate")
        }
        [command, target] if command == "leb2" => matches!(target.as_str(), "verify" | "convert"),
        [command, target] if command == "generate" => {
            planet_centers::parse(target, &mut std::iter::empty()).is_ok()
                || matches!(target.as_str(), "keplerian-elements" | "keplerian-dry-run")
        }
        _ => false,
    };
    if help && known_help_target {
        println!(
            "libephemeris {}\n{USAGE}",
            crate::identity::distribution_version()
        );
        return EXIT_OK;
    }
    match parse_args(&argv) {
        Ok(Command::GenerateKeplerian(request)) => keplerian::run(request),
        Ok(Command::GeneratePlanetCenters(request)) => planet_centers::run(request),
        Ok(Command::Provision(request)) => provision::run(request),
        Ok(Command::Leb2Convert(request)) => leb_conversion::run(request),
        Ok(Command::LebGenerate(request)) => leb_generation::run(request, generator),
        Ok(Command::Version) => {
            println!(
                "libephemeris, version {}",
                crate::identity::distribution_version()
            );
            EXIT_OK
        }
        Ok(Command::LebInfo { path, json }) => inspection::run_leb_info(&path, json),
        Ok(Command::Status {
            tier,
            data_dir,
            json,
            verbose,
        }) => inspection::run_status(tier, data_dir, json, verbose),
        Ok(Command::DownloadAssist {
            target_dir,
            data_dir,
            planets,
            asteroids,
            force,
            quiet,
            no_progress,
        }) => auxiliary_downloads::run_download_assist(
            target_dir,
            data_dir,
            planets,
            asteroids,
            force,
            quiet,
            no_progress,
        ),
        Ok(Command::DownloadIers {
            data_dir,
            force,
            quiet,
        }) => auxiliary_downloads::run_download_iers(data_dir, force, quiet),
        Ok(Command::Config { json }) => inspection::run_config(json),
        Ok(Command::Init {
            output,
            force,
            non_interactive,
        }) => {
            if non_interactive {
                wizard::run_init(&output, force)
            } else {
                wizard::run_init_wizard(&output, force)
            }
        }
        Ok(Command::GenerateKeplerianDryRun {
            spk_dirs,
            data_dir,
            start,
            end,
            spacing,
            bodies,
        }) => run_generate_keplerian_dry_run(&spk_dirs, data_dir, start, end, spacing, &bodies),
        Ok(Command::LebVerify {
            target,
            data_dir,
            json,
        }) => leb_verification::run_leb_verify(&target, data_dir, json),
        Ok(Command::Leb2Verify {
            path,
            reference,
            samples,
            json,
        }) => leb_verification::run_leb2_verify(&path, reference, samples, json),
        Err(message) => {
            eprintln!("libephemeris: {message}\n{USAGE}");
            EXIT_USAGE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn argv(words: &[&str]) -> Vec<String> {
        words.iter().map(ToString::to_string).collect()
    }

    #[test]
    fn version_spellings_parse() {
        for spelling in ["--version", "-V", "version"] {
            assert_eq!(parse_args(&argv(&[spelling])), Ok(Command::Version));
        }
    }

    #[test]
    fn leb_info_needs_a_path() {
        assert!(parse_args(&argv(&["leb", "info"])).is_err());
        assert_eq!(
            parse_args(&argv(&["leb", "info", "a.leb2"])),
            Ok(Command::LebInfo {
                path: "a.leb2".to_string(),
                json: false
            })
        );
        assert_eq!(
            parse_args(&argv(&["leb", "info", "a.leb2", "--json"])),
            Ok(Command::LebInfo {
                path: "a.leb2".to_string(),
                json: true
            })
        );
        assert!(parse_args(&argv(&["leb", "info", "a", "b"])).is_err());
        assert!(parse_args(&argv(&["leb", "bogus"])).is_err());
    }

    #[test]
    fn status_defaults_to_medium_and_parses_tiers() {
        assert_eq!(
            parse_args(&argv(&["status"])),
            Ok(Command::Status {
                tier: None,
                data_dir: None,
                json: false,
                verbose: 0,
            })
        );
        assert_eq!(
            parse_args(&argv(&["status", "--tier", "base", "--json"])),
            Ok(Command::Status {
                tier: Some(LebTier::Base),
                data_dir: None,
                json: true,
                verbose: 0,
            })
        );
        assert!(parse_args(&argv(&["status", "--tier", "huge"])).is_err());
        assert!(parse_args(&argv(&["status", "--tier"])).is_err());
        assert!(parse_args(&argv(&["bogus"])).is_err());
        assert!(parse_args(&[]).is_err());
    }

    #[test]
    fn data_directory_uses_the_shared_native_configuration() {
        let sources = || {
            crate::config::ConfigurationSources::load(
                std::collections::BTreeMap::from([(
                    "LIBEPHEMERIS_DATA_DIR".to_owned(),
                    "/env".to_owned(),
                )]),
                None,
                Some(PathBuf::from("/home")),
            )
        };
        assert_eq!(
            resolve_data_dir(Some("/flag"), sources()),
            Ok(PathBuf::from("/flag"))
        );
        assert_eq!(resolve_data_dir(None, sources()), Ok(PathBuf::from("/env")));
        let no_override = crate::config::ConfigurationSources::load(
            Default::default(),
            None,
            Some(PathBuf::from("/home")),
        );
        assert_eq!(
            resolve_data_dir(None, no_override),
            Ok(PathBuf::from("/home/.libephemeris"))
        );
    }

    #[test]
    fn embedding_entry_returns_exit_codes_without_exiting() {
        assert_eq!(run(["version"]), EXIT_OK);
        assert_eq!(run(["--help"]), EXIT_OK);
        assert_eq!(run(["unrecognized-command"]), EXIT_USAGE);
    }

    #[test]
    fn json_escape_covers_control_characters() {
        assert_eq!(
            json_escape("a\"b\\c\nd\te\x01f"),
            "a\\\"b\\\\c\\nd\\te\\u0001f"
        );
    }

    #[test]
    fn download_targets_parse_with_flags() {
        for target in ["leb2-base", "leb2-medium", "leb2-extended", "medium", "all"] {
            assert!(matches!(
                parse_args(&argv(&[
                    "download",
                    target,
                    "--force",
                    "-q",
                    "--no-progress",
                    "--data-dir",
                    "/data"
                ])),
                Ok(Command::Provision(_))
            ));
        }
        for args in [
            vec!["download"],
            vec!["download", "leb2-huge"],
            vec!["download", "leb2-base", "--data-dir"],
            vec!["download", "leb2-base", "--bogus"],
        ] {
            assert!(parse_args(&argv(&args)).is_err());
        }
    }

    #[test]
    fn assist_targets_parse_with_selection_flags() {
        assert_eq!(
            parse_args(&argv(&["download", "assist"])),
            Ok(Command::DownloadAssist {
                target_dir: None,
                data_dir: None,
                planets: true,
                asteroids: true,
                force: false,
                quiet: false,
                no_progress: false,
            })
        );
        assert_eq!(
            parse_args(&argv(&[
                "download",
                "assist",
                "--no-planets",
                "--no-asteroids",
                "--target-dir",
                "/assist",
                "--data-dir",
                "/data",
                "--force",
                "-q",
                "--no-progress",
            ])),
            Ok(Command::DownloadAssist {
                target_dir: Some("/assist".to_string()),
                data_dir: Some("/data".to_string()),
                planets: false,
                asteroids: false,
                force: true,
                quiet: true,
                no_progress: true,
            })
        );
        assert!(parse_args(&argv(&["download", "assist", "--target-dir"])).is_err());
        assert!(parse_args(&argv(&["download", "assist", "--data-dir"])).is_err());
        assert!(parse_args(&argv(&["download", "assist", "--bogus"])).is_err());
        assert!(parse_args(&argv(&["download", "assist", "--force", "extra"])).is_err());
        // Explicit provisioning targets are parsed without acquiring data.
        for target in [
            "base",
            "medium",
            "extended",
            "auto",
            "all",
            "leb-base",
            "leb-medium",
            "leb-extended",
        ] {
            assert!(
                parse_args(&argv(&["download", target])).is_ok(),
                "target: {target}"
            );
        }
    }

    #[test]
    fn iers_targets_parse_with_shared_flags() {
        assert_eq!(
            parse_args(&argv(&["download", "iers"])),
            Ok(Command::DownloadIers {
                data_dir: None,
                force: false,
                quiet: false,
            })
        );
        assert_eq!(
            parse_args(&argv(&[
                "download",
                "iers",
                "--force",
                "-q",
                "--data-dir",
                "/data"
            ])),
            Ok(Command::DownloadIers {
                data_dir: Some("/data".to_string()),
                force: true,
                quiet: true,
            })
        );
        assert!(parse_args(&argv(&["download", "iers", "--data-dir"])).is_err());
        assert!(parse_args(&argv(&["download", "iers", "--bogus"])).is_err());
        assert!(parse_args(&argv(&["download", "iers", "--force", "extra"])).is_err());
    }

    #[test]
    fn historical_download_aliases_share_the_canonical_parser() {
        for (alias, target) in [
            ("download:base", "base"),
            ("download:medium", "medium"),
            ("download:extended", "extended"),
            ("download:leb:base", "leb-base"),
            ("download:leb:medium", "leb-medium"),
            ("download:leb:extended", "leb-extended"),
            ("download:assist", "assist"),
        ] {
            assert_eq!(
                parse_args(&argv(&[alias, "--quiet"])),
                parse_args(&argv(&["download", target, "--quiet"]))
            );
        }
        assert!(parse_args(&argv(&["download", "leb2-base", "--no-progress"])).is_ok());
    }

    #[test]
    fn init_parses_file_flags_and_routes_bare_init_to_the_wizard() {
        assert_eq!(
            parse_args(&argv(&["init", "--non-interactive"])),
            Ok(Command::Init {
                output: "libephemeris-config.toml".to_string(),
                force: false,
                non_interactive: true,
            })
        );
        assert_eq!(
            parse_args(&argv(&[
                "init",
                "-o",
                "/tmp/cfg.toml",
                "--force",
                "--non-interactive"
            ])),
            Ok(Command::Init {
                output: "/tmp/cfg.toml".to_string(),
                force: true,
                non_interactive: true,
            })
        );
        // Bare `init` runs the interactive wizard; only unknown flags and
        // a missing `--output` value are usage errors.
        assert_eq!(
            parse_args(&argv(&["init"])),
            Ok(Command::Init {
                output: "libephemeris-config.toml".to_string(),
                force: false,
                non_interactive: false,
            })
        );
        assert!(parse_args(&argv(&["init", "--output"])).is_err());
        assert!(parse_args(&argv(&["init", "--non-interactive", "--bogus"])).is_err());
    }

    #[test]
    fn config_parses_an_optional_json_flag() {
        assert_eq!(
            parse_args(&argv(&["config"])),
            Ok(Command::Config { json: false })
        );
        assert_eq!(
            parse_args(&argv(&["config", "--json"])),
            Ok(Command::Config { json: true })
        );
        assert!(parse_args(&argv(&["config", "--tier", "base"])).is_err());
    }

    #[test]
    fn restored_planet_and_element_generators_are_explicit() {
        for spelling in [
            "planet-centers-base",
            "planet-centers-medium",
            "planet-centers-extended",
            "planet-centers-all",
            "planet-centers-spk",
        ] {
            assert!(matches!(
                parse_args(&argv(&["generate", spelling])),
                Ok(Command::GeneratePlanetCenters(_))
            ));
        }
        assert!(parse_args(&argv(&["generate"])).is_err());
        assert!(parse_args(&argv(&["generate", "bogus"])).is_err());
        assert!(parse_args(&argv(&["generate", "keplerian-elements", "extra"])).is_err());
        assert!(matches!(
            parse_args(&argv(&["generate", "keplerian-elements"])),
            Ok(Command::GenerateKeplerian(_))
        ));
    }

    #[test]
    fn generate_dry_run_parses_defaults_and_flags() {
        assert_eq!(
            parse_args(&argv(&["generate", "keplerian-dry-run"])),
            Ok(Command::GenerateKeplerianDryRun {
                spk_dirs: Vec::new(),
                data_dir: None,
                start: 1650,
                end: 2450,
                spacing: 10,
                bodies: Vec::new(),
            })
        );
        assert_eq!(
            parse_args(&argv(&[
                "generate",
                "keplerian-dry-run",
                "--spk-dir",
                "/kernels",
                "--data-dir",
                "/data",
                "--start",
                "1900",
                "--end",
                "2100",
                "--spacing",
                "50",
                "--body",
                "Ceres",
            ])),
            Ok(Command::GenerateKeplerianDryRun {
                spk_dirs: vec!["/kernels".to_string()],
                data_dir: Some("/data".to_string()),
                start: 1900,
                end: 2100,
                spacing: 50,
                bodies: vec!["Ceres".to_string()],
            })
        );
        assert!(parse_args(&argv(&["generate", "keplerian-dry-run", "--start"])).is_err());
        assert!(parse_args(&argv(&["generate", "keplerian-dry-run", "--start", "soon"])).is_err());
        assert!(
            parse_args(&argv(&[
                "generate",
                "keplerian-dry-run",
                "--start",
                "2100",
                "--end",
                "1900"
            ]))
            .is_err()
        );
        assert!(parse_args(&argv(&["generate", "keplerian-dry-run", "--spacing", "0"])).is_err());
        assert!(parse_args(&argv(&["generate", "keplerian-dry-run", "--body", "Nope"])).is_err());
        assert!(parse_args(&argv(&["generate", "keplerian-dry-run", "--bogus"])).is_err());
    }

    #[test]
    fn leb_generate_routes_every_supported_mode_to_the_offline_application() {
        for words in [
            vec!["leb", "generate", "medium", "groups"],
            vec!["leb", "generate", "base", "full"],
            vec!["leb", "generate", "extended", "single"],
            vec!["leb", "generate", "medium", "body", "moon,SUN,14"],
            vec!["leb", "generate", "medium", "planets"],
            vec!["leb", "generate", "medium", "asteroids"],
            vec!["leb", "generate", "medium", "exotics"],
            vec!["leb", "generate", "medium", "analytical"],
            vec!["leb", "generate", "medium", "merge"],
            vec!["leb", "generate", "all", "--data-dir", "/local/data"],
        ] {
            assert!(
                matches!(parse_args(&argv(&words)), Ok(Command::LebGenerate(_))),
                "{words:?}"
            );
        }
    }

    #[test]
    fn leb_generate_rejects_bad_tiers_modes_and_arities() {
        for words in [
            vec!["leb", "generate"],
            vec!["leb", "generate", "bogus", "groups"],
            vec!["leb", "generate", "medium"],
            vec!["leb", "generate", "medium", "bogus"],
            vec!["leb", "generate", "medium", "groups", "extra"],
            vec!["leb", "generate", "medium", "body"],
            vec!["leb", "generate", "medium", "body", "a", "b"],
            vec!["leb", "generate", "all", "extra"],
            vec!["leb", "generate", "medium", "full", "extra"],
        ] {
            assert!(parse_args(&argv(&words)).is_err(), "words: {words:?}");
        }
        assert!(
            parse_args(&argv(&["leb", "generate", "medium", "body"]))
                .expect_err("body needs a name")
                .contains("needs a body name")
        );
    }

    #[test]
    fn leb_verify_parses_tier_path_and_flags() {
        assert_eq!(
            parse_args(&argv(&["leb", "verify", "medium"])),
            Ok(Command::LebVerify {
                target: "medium".to_string(),
                data_dir: None,
                json: false,
            })
        );
        assert_eq!(
            parse_args(&argv(&[
                "leb",
                "verify",
                "base",
                "--data-dir",
                "/data",
                "--json"
            ])),
            Ok(Command::LebVerify {
                target: "base".to_string(),
                data_dir: Some("/data".to_string()),
                json: true,
            })
        );
        assert_eq!(
            parse_args(&argv(&["leb", "verify", "custom.leb"])),
            Ok(Command::LebVerify {
                target: "custom.leb".to_string(),
                data_dir: None,
                json: false,
            })
        );
        for words in [
            vec!["leb", "verify"],
            vec!["leb", "verify", "a", "b"],
            vec!["leb", "verify", "a", "--data-dir"],
        ] {
            assert!(parse_args(&argv(&words)).is_err(), "words: {words:?}");
        }
    }

    #[test]
    fn leb2_convert_spellings_dispatch_to_the_native_converter() {
        for target in [
            "base",
            "extended",
            "medium-core",
            "base-apogee",
            "extended-exotics",
        ] {
            assert!(
                matches!(
                    parse_args(&argv(&["leb2", "convert", target])),
                    Ok(Command::Leb2Convert(_))
                ),
                "target: {target}"
            );
        }
    }

    #[test]
    fn leb2_rejects_bad_subcommands_targets_and_arities() {
        for words in [
            vec!["leb2"],
            vec!["leb2", "bogus"],
            vec!["leb2", "convert"],
            vec!["leb2", "convert", "medium", "extra"],
            vec!["leb2", "convert", "bogus"],
            vec!["leb2", "convert", "medium-bogus"],
            vec!["leb2", "convert", "bogus-core"],
            vec!["leb2", "convert", "medium-core-extra"],
            vec!["leb2", "verify"],
            vec!["leb2", "verify", "a", "b"],
            vec!["leb2", "verify", "a", "--reference"],
            vec!["leb2", "verify", "a", "--samples"],
            vec!["leb", "bogus"],
            vec!["leb"],
        ] {
            assert!(parse_args(&argv(&words)).is_err(), "words: {words:?}");
        }
    }

    #[test]
    fn leb2_verify_parses_reference_samples_and_json() {
        assert_eq!(
            parse_args(&argv(&["leb2", "verify", "a.leb2"])),
            Ok(Command::Leb2Verify {
                path: "a.leb2".to_string(),
                reference: None,
                samples: leb_verification::LEB2_VERIFY_DEFAULT_SAMPLES,
                json: false,
            })
        );
        assert_eq!(
            parse_args(&argv(&[
                "leb2",
                "verify",
                "a.leb2",
                "--reference",
                "b.leb",
                "--samples",
                "128",
                "--json",
            ])),
            Ok(Command::Leb2Verify {
                path: "a.leb2".to_string(),
                reference: Some("b.leb".to_string()),
                samples: 128,
                json: true,
            })
        );
        for bad in ["0", "5000", "bogus", "-3"] {
            let message = parse_args(&argv(&["leb2", "verify", "a.leb2", "--samples", bad]))
                .expect_err("bad sample count");
            assert!(message.contains("--samples"), "samples={bad}: {message}");
        }
    }
}
