//! Explicit data provisioning shared by the native and Python console entries.

use crate::config::{ConfigurationSources, Runtime};
use crate::context::{EffectiveNetworkPolicy, ProviderMode};
use crate::errors::{Error, Result};
use crate::leb::LebTier;
use crate::net::TransportConfig;
use std::io::{self, Write};
use std::path::{Path, PathBuf};

/// Provisioning intent, independent of terminal spelling or calculation mode.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) enum Target {
    Tier(LebTier),
    Core(LebTier),
    Groups(LebTier),
    PlanetCenterSources(Option<LebTier>),
    Auto,
    All,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct Request {
    pub target: Target,
    pub directory: Option<PathBuf>,
    pub force: bool,
    pub quiet: bool,
    pub yes: bool,
}

/// Parse explicit provisioning requests without performing any filesystem I/O.
pub(super) fn parse<'a>(
    target: &str,
    arguments: &mut impl Iterator<Item = &'a str>,
) -> std::result::Result<Option<Request>, String> {
    let target = match target {
        "auto" => Target::Auto,
        "all" => Target::All,
        "planet-centers-sources" => Target::PlanetCenterSources(None),
        _ => {
            if let Some(tier) = LebTier::parse(target) {
                Target::Tier(tier)
            } else if let Some(tier) = target.strip_prefix("leb-").and_then(LebTier::parse) {
                Target::Core(tier)
            } else if let Some(tier) = target.strip_prefix("leb2-").and_then(LebTier::parse) {
                Target::Groups(tier)
            } else {
                return Ok(None);
            }
        }
    };
    let mut request = Request {
        target,
        directory: None,
        force: false,
        quiet: false,
        yes: false,
    };
    while let Some(argument) = arguments.next() {
        match argument {
            "--force" | "-f" => request.force = true,
            "--quiet" | "-q" => request.quiet = true,
            // Output reports completed files; there is no live terminal bar.
            "--no-progress" => {}
            "--yes" | "-y" if request.target == Target::All => request.yes = true,
            "--tier" if matches!(request.target, Target::PlanetCenterSources(_)) => {
                let Target::PlanetCenterSources(tier) = &mut request.target else {
                    unreachable!("guarded planet-center source request")
                };
                if tier.is_some() {
                    return Err("duplicate --tier".to_owned());
                }
                let value = arguments.next().ok_or("--tier needs a value")?;
                *tier = Some(
                    LebTier::parse(&value.to_ascii_lowercase())
                        .ok_or_else(|| format!("unknown tier '{value}'"))?,
                );
            }
            "--data-dir" => {
                request.directory = Some(PathBuf::from(
                    arguments.next().ok_or("--data-dir needs a value")?,
                ))
            }
            _ => return Err(format!("unexpected argument '{argument}'")),
        }
    }
    Ok(Some(request))
}

fn install_tier(runtime: &Runtime, tier: &str, force: bool, quiet: bool) -> Result<bool> {
    let report = runtime.download_for_tier(tier, force)?;
    if !quiet {
        println!(
            "{tier}: {} {}, {} cached, {} downloaded, {} fallback, {} errors",
            report.summary.ephemeris_file,
            report.summary.ephemeris_status,
            report.summary.cached,
            report.summary.downloaded,
            report.summary.fallback,
            report.summary.errors
        );
    }
    Ok(report.summary.errors == 0)
}

/// Require every published group; partial API results cannot imply CLI success.
fn install_groups(runtime: &Runtime, tier: &str, force: bool, quiet: bool) -> Result<()> {
    let paths = runtime.download_leb2_for_tier(tier, None, force, false)?;
    let directory = runtime.data_directory()?;
    let mut missing = Vec::new();
    for group in crate::leb::LEB2_GROUPS {
        if let Some(entry) = crate::download::leb2_asset(tier, group)? {
            if !paths.contains(&crate::download::dest_path_for(&directory, entry)) {
                missing.push(entry.name());
            }
        }
    }
    if !missing.is_empty() {
        return Err(Error::AssetRejected {
            filename: missing.join(", "),
            reason: "tier coefficient installation is incomplete",
        });
    }
    if !quiet {
        for path in paths {
            println!("  [OK] {}", path.display());
        }
    }
    Ok(())
}

/// Complete explicit bootstrap order. The coefficient/DE/SPK steps remain
/// interleaved per tier; the shared source cache and ASSIST follow IERS.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BootstrapPhase {
    Groups(LebTier),
    Tier(LebTier),
    Iers,
    PlanetCenterSources,
    Assist,
}

const ALL_PHASES: [BootstrapPhase; 9] = [
    BootstrapPhase::Groups(LebTier::Base),
    BootstrapPhase::Tier(LebTier::Base),
    BootstrapPhase::Groups(LebTier::Medium),
    BootstrapPhase::Tier(LebTier::Medium),
    BootstrapPhase::Groups(LebTier::Extended),
    BootstrapPhase::Tier(LebTier::Extended),
    BootstrapPhase::Iers,
    BootstrapPhase::PlanetCenterSources,
    BootstrapPhase::Assist,
];

/// Perform one bootstrap step using the captured runtime, directory and policy.
/// Typed errors propagate; no later phase runs after a failed step.
fn perform_phase(
    phase: BootstrapPhase,
    request: &Request,
    runtime: &Runtime,
    directory: &Path,
    policy: EffectiveNetworkPolicy,
    transport: &TransportConfig,
) -> Result<bool> {
    match phase {
        BootstrapPhase::Groups(tier) => {
            install_groups(runtime, tier.name(), request.force, request.quiet)?;
        }
        BootstrapPhase::Tier(tier) => {
            return install_tier(runtime, tier.name(), request.force, request.quiet);
        }
        BootstrapPhase::Iers => {
            runtime.set_leb_file(None);
            let cache = crate::download::iers_cache_dir(directory);
            if !request.quiet {
                println!("libephemeris download iers → {}", cache.display());
            }
            let report =
                crate::download::refresh_iers_cache(&cache, request.force, policy, transport)?;
            if !request.quiet {
                for file in &report.files {
                    println!("  [OK] {}: {}", file.name, file.path.display());
                }
                println!("download iers complete");
            }
        }
        BootstrapPhase::PlanetCenterSources => {
            super::planet_center_sources::install(runtime, None, request.force, request.quiet)?;
        }
        BootstrapPhase::Assist => {
            let destination = directory.join("assist");
            if !request.quiet {
                println!("libephemeris download assist → {}", destination.display());
            }
            let report = crate::download::download_assist_data(
                &destination,
                true,
                true,
                false,
                request.force,
                policy,
                transport,
                None,
            )?;
            if !request.quiet {
                println!(
                    "  [OK] assist: {} downloaded, {} kept",
                    report.downloaded, report.skipped
                );
                println!("download assist complete");
            }
        }
    }
    Ok(true)
}

fn perform(request: &Request, runtime: &Runtime) -> Result<bool> {
    match &request.target {
        Target::PlanetCenterSources(tier) => {
            super::planet_center_sources::install(runtime, *tier, request.force, request.quiet)?;
            Ok(true)
        }
        Target::Core(tier) => {
            let path = runtime.download_leb_for_tier(tier.name(), request.force, false)?;
            if !request.quiet {
                println!("  [OK] {}", path.display());
            }
            Ok(true)
        }
        Target::Tier(tier) => install_tier(runtime, tier.name(), request.force, request.quiet),
        Target::Groups(tier) => {
            install_groups(runtime, tier.name(), request.force, request.quiet)?;
            Ok(true)
        }
        Target::Auto => {
            let context = runtime.effective_context()?;
            let tier = context.tier().as_name();
            match context.mode() {
                ProviderMode::Auto | ProviderMode::SealedLeb => {
                    let tier = LebTier::parse(tier).expect("resolved core tier");
                    for eligible in crate::leb::eligible_tiers(tier) {
                        install_groups(runtime, eligible.name(), request.force, request.quiet)?;
                    }
                    runtime.set_leb_file(None);
                }
                _ => {}
            }
            match context.mode() {
                ProviderMode::Auto | ProviderMode::NativeJpl => {
                    install_tier(runtime, tier, request.force, request.quiet)
                }
                ProviderMode::HorizonsRemote => {
                    let name = crate::download::planet_centers_name(tier)?;
                    let entry = crate::download::asset_by_name(&name)
                        .ok_or_else(|| Error::DataManifestMissing { name: name.clone() })?;
                    crate::download::ensure_asset(
                        entry,
                        &crate::download::dest_path_for(&runtime.data_directory()?, entry),
                        request.force,
                        context.effective_network_policy(),
                        &crate::net::TransportConfig::default(),
                        crate::download::validator_for(&name),
                        None,
                    )?;
                    Ok(true)
                }
                ProviderMode::SealedLeb => Ok(true),
            }
        }
        Target::All => {
            let directory = runtime.data_directory()?;
            let policy = runtime.effective_context()?.effective_network_policy();
            let transport = TransportConfig::default();
            for phase in ALL_PHASES {
                if !perform_phase(phase, request, runtime, &directory, policy, &transport)? {
                    return Ok(false);
                }
            }
            Ok(true)
        }
    }
}

/// Capture configuration once and apply the CLI directory above configured roots.
/// No data are opened and the effective network policy is preserved.
pub(super) fn configured_runtime(directory: Option<&Path>) -> Result<Runtime> {
    let configured = Runtime::with_sources(ConfigurationSources::discover());
    let mut context = configured.effective_context()?;
    if let Some(directory) = directory {
        context.set_jpl_data_dir(Some(directory.to_path_buf()));
        context.set_jpl_dir(Some(directory.to_path_buf()));
        context.set_iers_cache_dir(Some(directory.join("iers_cache")));
        context.set_spk_cache_dir(Some(directory.join("spk")));
    }
    Ok(Runtime::from_context(context))
}

/// Execute one explicitly authorized provisioning request and return an exit code.
pub(super) fn run(request: Request) -> i32 {
    if request.target == Target::All && !request.yes {
        print!(
            "Download all tiers, IERS, planet-center sources and ASSIST data (several GB)? [y/N] "
        );
        if io::stdout().flush().is_err() {
            return super::EXIT_FAILED;
        }
        let mut answer = String::new();
        if io::stdin().read_line(&mut answer).is_err() {
            return super::EXIT_FAILED;
        }
        if !matches!(answer.trim().to_ascii_lowercase().as_str(), "y" | "yes") {
            return super::EXIT_FAILED;
        }
    }
    // Explicit provisioning is the network-authorizing CLI boundary.
    let outcome = configured_runtime(request.directory.as_deref()).and_then(|runtime| {
        runtime.set_network_policy(Some("allow"))?;
        perform(&request, &runtime)
    });
    match outcome {
        Ok(true) => super::EXIT_OK,
        Ok(false) => super::EXIT_FAILED,
        Err(error) => {
            eprintln!("libephemeris: download failed: {error}");
            super::EXIT_FAILED
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scratch(tag: &str) -> PathBuf {
        use std::sync::atomic::{AtomicU64, Ordering};
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "libephemeris-bootstrap-{}-{tag}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&path).unwrap();
        path
    }

    #[test]
    fn all_bootstrap_keeps_each_tier_then_iers_source_and_assist_prerequisites() {
        assert_eq!(
            ALL_PHASES,
            [
                BootstrapPhase::Groups(LebTier::Base),
                BootstrapPhase::Tier(LebTier::Base),
                BootstrapPhase::Groups(LebTier::Medium),
                BootstrapPhase::Tier(LebTier::Medium),
                BootstrapPhase::Groups(LebTier::Extended),
                BootstrapPhase::Tier(LebTier::Extended),
                BootstrapPhase::Iers,
                BootstrapPhase::PlanetCenterSources,
                BootstrapPhase::Assist,
            ]
        );
    }

    #[test]
    fn bootstrap_supplement_failures_keep_the_native_sealed_error() {
        let directory = scratch("sealed-supplements");
        let runtime = Runtime::from_context(
            crate::context::EphemerisContext::default().with_jpl_data_dir(Some(directory.clone())),
        );
        runtime.set_network_policy(Some("sealed")).unwrap();
        let request = Request {
            target: Target::All,
            directory: None,
            force: false,
            quiet: true,
            yes: true,
        };
        for phase in [
            BootstrapPhase::Iers,
            BootstrapPhase::PlanetCenterSources,
            BootstrapPhase::Assist,
        ] {
            let result = perform_phase(
                phase,
                &request,
                &runtime,
                &directory,
                EffectiveNetworkPolicy::Sealed,
                &TransportConfig::default(),
            );
            assert!(
                matches!(result, Err(Error::NetworkSealed { .. })),
                "{phase:?}: {result:?}"
            );
        }
        assert!(!directory.join("tmp/planet_centers").exists());
        assert_eq!(
            std::fs::read_dir(directory.join("assist")).unwrap().count(),
            0
        );
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[cfg(feature = "bundled-data")]
    #[test]
    fn failed_first_bootstrap_group_stops_before_later_phases() {
        let directory = scratch("stop-first");
        let sentinel = directory.join("preserved");
        std::fs::write(&sentinel, b"previous local data").unwrap();
        let runtime = Runtime::from_context(
            crate::context::EphemerisContext::default().with_jpl_data_dir(Some(directory.clone())),
        );
        runtime.set_network_policy(Some("sealed")).unwrap();
        let request = Request {
            target: Target::All,
            directory: None,
            force: false,
            quiet: true,
            yes: true,
        };
        assert!(matches!(
            perform(&request, &runtime),
            Err(Error::AssetRejected { .. })
        ));
        assert!(!directory.join("de440s.bsp").exists());
        assert!(!directory.join("iers_cache").exists());
        assert!(!directory.join("tmp/planet_centers").exists());
        assert!(!directory.join("assist").exists());
        assert_eq!(std::fs::read(sentinel).unwrap(), b"previous local data");
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn targets_options_and_invalid_requests_are_distinct() {
        for target in [
            "base",
            "medium",
            "extended",
            "leb-base",
            "leb-medium",
            "leb-extended",
            "auto",
            "all",
        ] {
            let request = parse(
                target,
                &mut ["-f", "-q", "--no-progress", "--data-dir", "/explicit"].into_iter(),
            )
            .unwrap()
            .unwrap();
            assert_eq!(request.directory.as_deref(), Some(Path::new("/explicit")));
            assert!(request.force && request.quiet);
        }
        assert!(parse("all", &mut ["-y"].into_iter()).unwrap().unwrap().yes);
        assert!(parse("medium", &mut ["--yes"].into_iter()).is_err());
        assert!(parse("auto", &mut ["--data-dir"].into_iter()).is_err());
        assert!(parse("leb-future", &mut [].into_iter()).unwrap().is_none());
    }

    #[test]
    fn planet_center_sources_defaults_to_all_tiers_and_rejects_invalid_selection() {
        let request = parse("planet-centers-sources", &mut [].into_iter())
            .unwrap()
            .unwrap();
        assert_eq!(request.target, Target::PlanetCenterSources(None));
        assert!(!request.force && !request.quiet);
        for (name, expected) in [
            ("BASE", LebTier::Base),
            ("medium", LebTier::Medium),
            ("Extended", LebTier::Extended),
        ] {
            let request = parse(
                "planet-centers-sources",
                &mut [
                    "--tier",
                    name,
                    "--force",
                    "--quiet",
                    "--data-dir",
                    "/explicit",
                ]
                .into_iter(),
            )
            .unwrap()
            .unwrap();
            assert_eq!(request.target, Target::PlanetCenterSources(Some(expected)));
            assert_eq!(request.directory.as_deref(), Some(Path::new("/explicit")));
            assert!(request.force && request.quiet);
        }
        for arguments in [
            vec!["--tier"],
            vec!["--tier", "unknown"],
            vec!["--tier", "base", "--tier", "medium"],
            vec!["--unknown"],
            vec!["--yes"],
        ] {
            assert!(parse("planet-centers-sources", &mut arguments.into_iter()).is_err());
        }
        assert!(parse("medium", &mut ["--tier", "base"].into_iter()).is_err());
    }

    #[cfg(feature = "bundled-data")]
    #[test]
    fn core_compatibility_target_installs_the_bundled_image_offline() {
        let directory =
            std::env::temp_dir().join(format!("libephemeris-cli-core-{}", std::process::id()));
        let runtime = Runtime::from_context(
            crate::context::EphemerisContext::default().with_jpl_data_dir(Some(directory.clone())),
        );
        runtime.set_network_policy(Some("sealed")).unwrap();
        let request = Request {
            target: Target::Core(LebTier::Base),
            directory: None,
            force: false,
            quiet: true,
            yes: false,
        };
        assert!(perform(&request, &runtime).unwrap());
        let path = directory.join("leb/base_core.leb2");
        assert_eq!(
            crate::download::file_sha256(&path).unwrap(),
            crate::download::asset_by_name("base_core.leb2")
                .unwrap()
                .sha256()
        );
        assert!(runtime.context().leb_file().is_none());
        assert!(perform(&request, &runtime).unwrap());
        // The Python API may return installed groups after an optional-group
        // failure. An all-group CLI request must expose that incompleteness.
        let groups = Request {
            target: Target::Groups(LebTier::Base),
            ..request
        };
        assert!(matches!(
            perform(&groups, &runtime),
            Err(Error::AssetRejected { .. })
        ));
        assert!(path.is_file());
        std::fs::remove_dir_all(directory).unwrap();
    }
}
