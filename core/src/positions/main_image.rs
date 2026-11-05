//! Owned acquisition of the selected main SPK image for a public native entry.
//!
//! This module performs no astronomical calculation. It selects paths from
//! caller-resolved configuration, applies the request's network policy, loads
//! an owned SPK image and registers its NAIF targets before a consuming clock.
//! An already populated registry is returned unchanged; other provider modes
//! never enter this helper. Filename fallback and flag fallback are separate.

use crate::context::{EffectiveNetworkPolicy, EphemerisContext};
use crate::errors::{Error, Result};
use std::path::{Path, PathBuf};

/// Download URL prefix for standard DE planetary SPK filenames.
const JPL_ROOT: &str = "https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/";
/// More-specific download prefix for filenames matching `jup*.bsp`.
const JUPITER_ROOT: &str = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/";

/// Entry-owned main image and its missing-filename fallback decision.
///
/// Registry images are owned by `context` and live through the synchronous
/// active calculation. No filesystem buffer or context borrow escapes.
/// `fallback_active` affects JPLEPH echo normalization, not physical formulas.
pub(super) struct AcquiredMainImage {
    pub(super) context: EphemerisContext,
    pub(super) fallback_active: bool,
    pub(super) main_owner: Option<super::main_owner::MainImageOwner>,
}

/// Find the first existing candidate: custom ephemeris directory, then the
/// explicitly resolved data directory. Absolute names override either prefix,
/// leaving an absolute pathname unprefixed. Presence includes directories: opening, not this
/// existence predicate, owns the later I/O/format refusal.
fn local_candidates<'a>(
    ctx: &'a EphemerisContext,
    filename: &'a str,
) -> impl Iterator<Item = PathBuf> + 'a {
    [
        ctx.paths().jpl_dir.as_deref(),
        ctx.jpl_data_dir().map(PathBuf::as_path),
    ]
    .into_iter()
    .flatten()
    .map(move |directory| directory.join(filename))
    .chain(std::iter::once(PathBuf::from(filename)).filter(|path| path.is_absolute()))
}

fn present_path(ctx: &EphemerisContext, filename: &str) -> Option<PathBuf> {
    local_candidates(ctx, filename).find(|candidate| candidate.exists())
}

/// Whether the literal selected main filename exists in its configured local
/// search roots. No open, fallback substitution, download or network work occurs.
pub(crate) fn selected_main_file_present(ctx: &EphemerisContext) -> bool {
    let requested = ctx
        .ephemeris_file_override()
        .unwrap_or_else(|| super::tier_de_kernel(ctx.tier()));
    local_candidates(ctx, requested).any(|path| path.is_file())
}

/// Filename/presence fallback decision without opening files or using network.
/// Explicit LEB mode seals JPL access and never reports a kernel fallback.
pub(crate) fn configured_jpl_fallback_active(ctx: &EphemerisContext) -> bool {
    if ctx.mode() == crate::context::ProviderMode::SealedLeb {
        return false;
    }
    let requested = ctx
        .ephemeris_file_override()
        .unwrap_or_else(|| super::tier_de_kernel(ctx.tier()));
    matches!(
        super::decide_jpl_file(
            super::is_recognized_de_kernel(requested),
            present_path(ctx, requested).is_some()
        ),
        super::JplFileDecision::TierDefault
    )
}

/// Resolve a filename into its literal download URL and opener pathname.
///
/// An explicit URL keeps its original bytes, with only its path's final
/// component used as the download filename. A pathname keeps its full bytes
/// in the URL rule: an absolute or nested name is not reduced to its basename.
/// The SPK extension is case-sensitive. Dot-only basenames do not name
/// supported files. Unsupported formats return `EphemerisLoaderValue` with
/// a native path diagnostic; no file or network access occurs here.
pub(crate) fn loader_target(filename: &str) -> Result<(String, String)> {
    let (url, name) = if let Some(scheme_end) = filename.find("://") {
        let rest = &filename[scheme_end + 3..];
        let path = rest.find('/').map_or("", |start| &rest[start..]);
        let path = path.split(['?', '#']).next().unwrap_or("");
        (
            filename.to_owned(),
            path.rsplit('/').next().unwrap_or("").to_owned(),
        )
    } else {
        let base = if filename.starts_with("jup") {
            JUPITER_ROOT
        } else {
            JPL_ROOT
        };
        (format!("{base}{filename}"), filename.to_owned())
    };
    let basename = name.rsplit('/').next().unwrap_or("");
    let bsp_extension = basename.rfind('.').is_some_and(|dot| {
        &basename[dot..] == ".bsp" && basename[..dot].bytes().any(|byte| byte != b'.')
    });
    if !bsp_extension {
        return Err(Error::EphemerisLoaderValue {
            message: format!("unsupported SPK file extension: {:?}", Path::new(&name)),
        });
    }
    Ok((url, name))
}

/// Load an optional, digest-verified center-offset image for the resolved tier.
/// Only existing local assets are considered; absent or rejected supplements
/// leave the DE system barycenters available. Exact opened bytes are checked.
pub(super) fn attach_automatic_centers(context: &mut EphemerisContext) -> Result<()> {
    let Some(directory) = context.jpl_data_dir().cloned() else {
        return Ok(());
    };
    let tier = match context.tier() {
        crate::context::Tier::Base => "base",
        crate::context::Tier::Medium => "medium",
        crate::context::Tier::Extended => "extended",
    };
    let mut names = vec![format!("planet_centers_{tier}.bsp")];
    if tier == "base" {
        names.push("planet_centers.bsp".to_owned());
    }
    for name in names {
        let path = directory.join(&name);
        if !path.is_file() {
            continue;
        }
        let Some(asset) = crate::download::asset_by_name(&name) else {
            continue;
        };
        let mut supplement = crate::spk::SpkRegistry::new();
        let key = match supplement.register_file(&path) {
            Ok(key) => key,
            Err(error) => {
                log::warn!("Ignoring planet-center image {}: {error}", path.display());
                continue;
            }
        };
        let digest = crate::download::hex_of(
            &supplement
                .kernel(&key)
                .expect("opened image")
                .image_sha256(),
        );
        if digest != asset.sha256() {
            log::warn!(
                "Ignoring planet-center image {}: digest does not match {name}",
                path.display()
            );
            continue;
        }
        if let Err(error) = supplement.register_automatic_planet_centers(&key) {
            log::warn!("Ignoring planet-center image {}: {error}", path.display());
            continue;
        }
        context
            .spk_registry_mut()
            .retain_automatic_centers_from(&supplement)?;
        break;
    }
    Ok(())
}

/// Acquire a selected main image, using the existing streaming atomic fetch.
///
/// Empty NativeJpl standard-planet or consumed stellar main-image registries
/// call this entry. Supplied nonempty registries remain authoritative.
/// A configured loader directory is initialized before filename selection,
/// including under sealed policy; local initialization errors precede policy
/// and opener failures. A present file is opened locally even under sealed
/// policy. Truly missing
/// selected files gate the exact `JPL {filename} ephemeris download`
/// purpose before HTTP/temporary-file work. Allowed requests perform the
/// actual exchange and then decode/register the owned result; no allowed
/// branch returns an empty registry as an apparent successful acquisition.
/// Downloads are streamed without a caller-provided expected digest. SPK structural decode,
/// I/O, loader-format, policy and transport errors propagate without retry.
pub(super) fn acquire_main_image(ctx: &EphemerisContext) -> Result<AcquiredMainImage> {
    acquire_main_image_with(ctx, |url, destination, purpose, policy| {
        crate::download::download_unpinned_file(
            url,
            destination,
            purpose,
            policy,
            &crate::net::TransportConfig::new()?,
        )
        .map(|_| ())
    })
}

/// Acquire an absent main image while retaining unrelated registered kernels.
/// A caller has already distinguished current main ownership from optional
/// minor/satellite registrations. Existing mappings remain owned; the normal
/// loader, policy, decode and target registration run on the configured main.
pub(super) fn acquire_missing_main_image(ctx: &EphemerisContext) -> Result<AcquiredMainImage> {
    acquire_missing_main_image_with_policy(ctx, &|| Ok(ctx.effective_network_policy()))
}

/// Acquire a main with transport policy resolved only for an absent file.
/// Local images consume no policy callback. A needed policy is admitted once
/// and retained through the loader/download stage; callback errors propagate.
pub(super) fn acquire_missing_main_image_with_policy(
    ctx: &EphemerisContext,
    policy: &dyn Fn() -> Result<EffectiveNetworkPolicy>,
) -> Result<AcquiredMainImage> {
    acquire_main_image_with_download_policy(
        ctx,
        false,
        |url, destination, purpose, policy| {
            crate::download::download_unpinned_file(
                url,
                destination,
                purpose,
                policy,
                &crate::net::TransportConfig::new()?,
            )
            .map(|_| ())
        },
        policy,
    )
}

/// Open a literal loader filename independently of current-main selection.
/// The acquired directory/filename select no tier fallback, body registration,
/// clock or main publication. A present file consumes no policy; a missing
/// file uses one live policy and atomic transport before SPK structural decode.
/// Format, directory/file I/O, policy, transport and decoder errors propagate.
pub(super) fn acquire_loader_spk_image_with_policy(
    directory: &Path,
    filename: &str,
    policy: &dyn Fn(&str) -> Result<EffectiveNetworkPolicy>,
) -> Result<(String, crate::spk::SpkRegistry)> {
    let (url, name) = loader_target(filename)?;
    let path = directory.join(&name);
    if !path.exists() {
        let purpose = format!("JPL data download: {name}");
        let selected = policy(&purpose)?;
        crate::context::require_network(selected, &purpose)?;
        crate::download::download_unpinned_file(
            &url,
            &path,
            &purpose,
            selected,
            &crate::net::TransportConfig::new()?,
        )?;
    }
    let mut registry = crate::spk::SpkRegistry::new();
    let key = registry.register_file(&path)?;
    Ok((key, registry))
}

/// Select and acquire the main image using a filesystem-producing fetch.
///
/// The callback receives a URL, destination path, exact selected
/// purpose and effective policy only after the network gate. It must complete
/// publication or return its actual failure. The production callback uses the
/// streaming atomic transport.
/// An explicitly configured loader root is created/validated before filename
/// selection. An SDK context without a data root
/// retains its absolute-local-file extension without reading process settings.
/// Local root I/O errors own their path, actual OS errno/kind and diagnostic.
/// All contexts/images are owned, with no global configuration or caches.
fn acquire_main_image_with(
    ctx: &EphemerisContext,
    fetch: impl FnMut(&str, &Path, &str, EffectiveNetworkPolicy) -> Result<()>,
) -> Result<AcquiredMainImage> {
    acquire_main_image_with_registry_policy(ctx, true, fetch)
}

/// Distinguish an already supplied main from unrelated optional registrations
/// without dropping their owners or altering policy. Ordinary consumers retain
/// their existing supplied-registry behavior through the preceding wrapper.
fn acquire_main_image_with_registry_policy(
    ctx: &EphemerisContext,
    retain_supplied: bool,
    fetch: impl FnMut(&str, &Path, &str, EffectiveNetworkPolicy) -> Result<()>,
) -> Result<AcquiredMainImage> {
    acquire_main_image_with_download_policy(ctx, retain_supplied, fetch, &|| {
        Ok(ctx.effective_network_policy())
    })
}

/// Shared local selection and loader, with one lazy transport-policy admission.
/// The policy is unused when the selected loader path already exists; its owned
/// value governs both absence checks and the actual streaming fetch.
fn acquire_main_image_with_download_policy(
    ctx: &EphemerisContext,
    retain_supplied: bool,
    mut fetch: impl FnMut(&str, &Path, &str, EffectiveNetworkPolicy) -> Result<()>,
    policy: &dyn Fn() -> Result<EffectiveNetworkPolicy>,
) -> Result<AcquiredMainImage> {
    if retain_supplied && !ctx.spk_registry().is_empty() {
        let mut context = ctx.clone();
        attach_automatic_centers(&mut context)?;
        return Ok(AcquiredMainImage {
            context,
            fallback_active: false,
            main_owner: None,
        });
    }
    if let Some(directory) = ctx.jpl_data_dir() {
        std::fs::create_dir_all(directory).map_err(|error| Error::EphemerisLoaderIo {
            path: directory.to_string_lossy().into_owned(),
            kind: error.kind(),
            errno: error.raw_os_error(),
            message: error.to_string(),
        })?;
    }
    let requested = ctx
        .ephemeris_file_override()
        .unwrap_or_else(|| super::tier_de_kernel(ctx.tier()));
    let requested_path = present_path(ctx, requested);
    let fallback_active = matches!(
        super::decide_jpl_file(
            super::is_recognized_de_kernel(requested),
            requested_path.is_some()
        ),
        super::JplFileDecision::TierDefault
    );
    let effective = if fallback_active {
        super::tier_de_kernel(ctx.tier())
    } else {
        requested
    };
    let selected = if fallback_active {
        present_path(ctx, effective)
    } else {
        requested_path
    };
    let purpose = format!("JPL {effective} ephemeris download");
    // A missing selected file requires network permission before an opener
    // is chosen; existing local images never consume network policy.
    let admitted_policy = if selected.is_none() {
        let value = policy()?;
        crate::context::require_network(value, &purpose)?;
        Some(value)
    } else {
        None
    };
    let call_name = selected.as_ref().map_or_else(
        || effective.to_owned(),
        |p| p.to_string_lossy().into_owned(),
    );
    let (url, loader_name) = loader_target(&call_name)?;
    let loader_path = if Path::new(&loader_name).is_absolute() {
        PathBuf::from(&loader_name)
    } else {
        ctx.jpl_data_dir()
            .ok_or(Error::ConfigurationMissing {
                missing_config: "native_main_data_dir",
                suggestion:
                    "set an explicit resolved JPL data directory before main-image acquisition",
            })?
            .join(&loader_name)
    };
    if !loader_path.exists() {
        // A relative custom-path call is resolved in the loader directory,
        // preserving its literal path/URL rather than guessing a basename.
        let value = match admitted_policy {
            Some(value) => value,
            None => policy()?,
        };
        crate::context::require_network(value, &purpose)?;
        fetch(&url, &loader_path, &purpose, value)?;
    }
    let mut context = ctx.clone();
    let key = context.spk_registry_mut().register_file(&loader_path)?;
    crate::spk::register_all_naif_targets(context.spk_registry_mut(), &key)?;
    // Retain caller registrations independently of automatic main aliases.
    for registration in ctx.spk_registry().registered_bodies() {
        context
            .spk_registry_mut()
            .retain_registration_from(ctx.spk_registry(), registration.body)?;
    }
    attach_automatic_centers(&mut context)?;
    let main_owner = super::main_owner::MainImageOwner::captured(
        context.clone(),
        key,
        requested.to_owned(),
        effective.to_owned(),
        fallback_active,
    )?;
    Ok(AcquiredMainImage {
        context,
        fallback_active,
        main_owner: Some(main_owner),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::context::{NetworkPolicy, ProviderMode, Tier};
    use std::sync::atomic::{AtomicU64, Ordering};

    /// Synthetic metadata-only image with complete Sun/Earth/Mars chains.
    /// Starting from the existing Rust fixture, change only each descriptor's
    /// target/center. No numerical directory or reference value is created.
    fn test_image() -> Vec<u8> {
        let mut image = crate::spk::moon_test_kernel_bytes();
        for (index, target) in [10i32, 399, 4].into_iter().enumerate() {
            let descriptor = 1048 + index * 40;
            image[descriptor + 16..descriptor + 20].copy_from_slice(&target.to_le_bytes());
            image[descriptor + 20..descriptor + 24].copy_from_slice(&0i32.to_le_bytes());
        }
        image
    }

    #[test]
    fn unreviewed_optional_center_image_does_not_replace_main_targets() {
        let data = Sandbox::new();
        data.image("de440.bsp");
        data.image("planet_centers_base.bsp");
        let acquired =
            acquire_main_image_with(&data.context(NetworkPolicy::Sealed), |_, _, _, _| {
                panic!("present images must not download")
            })
            .unwrap();
        let registry = acquired.context.spk_registry();
        assert!(registry.automatic_center_mappings().is_empty());
        assert_eq!(registry.kernel_count(), 1);
    }

    /// Owned unique sandbox; cleanup has no reference/consumer effects.
    struct Sandbox(PathBuf);
    impl Sandbox {
        fn new() -> Self {
            static ID: AtomicU64 = AtomicU64::new(0);
            let path = std::env::temp_dir().join(format!(
                "lephem-main-{}-{}",
                std::process::id(),
                ID.fetch_add(1, Ordering::Relaxed)
            ));
            std::fs::create_dir(&path).unwrap();
            Self(path)
        }
        fn context(&self, policy: NetworkPolicy) -> EphemerisContext {
            EphemerisContext::new()
                .with_mode(ProviderMode::NativeJpl)
                .with_network_policy(Some(policy))
                .with_jpl_data_dir(Some(self.0.clone()))
        }
        fn image(&self, name: &str) -> PathBuf {
            let path = self.0.join(name);
            if let Some(parent) = path.parent() {
                std::fs::create_dir_all(parent).unwrap();
            }
            std::fs::write(&path, test_image()).unwrap();
            path
        }
    }
    impl Drop for Sandbox {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    #[test]
    fn named_loader_does_not_recreate_retired_root_before_format_or_policy() {
        let data = Sandbox::new();
        std::fs::remove_dir(&data.0).unwrap();
        std::fs::write(&data.0, b"retired loader directory").unwrap();
        let result = acquire_loader_spk_image_with_policy(&data.0, "bad.txt", &|_| {
            panic!("unsupported opener must precede policy")
        });
        assert!(matches!(result, Err(Error::EphemerisLoaderValue { .. })));
        std::fs::remove_file(&data.0).unwrap();
        let result = acquire_loader_spk_image_with_policy(&data.0, "de441.bsp", &|purpose| {
            assert_eq!(purpose, "JPL data download: de441.bsp");
            Ok(EffectiveNetworkPolicy::Sealed)
        });
        assert!(matches!(result, Err(Error::NetworkSealed { .. })));
        assert!(!data.0.exists());
    }

    #[test]
    fn named_loader_keeps_literal_image_without_body_or_main_registration() {
        let data = Sandbox::new();
        let path = data.image("satellite.bsp");
        let (key, registry) =
            acquire_loader_spk_image_with_policy(&data.0, "satellite.bsp", &|_| {
                panic!("present loader image must not consume policy")
            })
            .unwrap();
        assert_eq!(key, path.to_string_lossy());
        assert!(registry.kernel(&key).is_some());
        assert!(registry.registered_bodies().is_empty());
    }

    #[test]
    fn named_loader_missing_de_does_not_substitute_present_tier_default() {
        let data = Sandbox::new();
        data.image("de440.bsp");
        let calls = std::cell::Cell::new(0);
        let result = acquire_loader_spk_image_with_policy(&data.0, "de441.bsp", &|_| {
            calls.set(calls.get() + 1);
            Ok(EffectiveNetworkPolicy::Sealed)
        });
        match result {
            Err(Error::NetworkSealed { purpose }) => {
                assert_eq!(purpose, "JPL data download: de441.bsp");
            }
            result => panic!("literal missing loader file must be refused: {result:?}"),
        }
        assert_eq!(calls.get(), 1);
        assert!(!data.0.join("de441.bsp").exists());
    }

    #[test]
    fn direct_main_local_image_never_consumes_download_policy() {
        let data = Sandbox::new();
        data.image("de440.bsp");
        let ctx = data.context(NetworkPolicy::Allow);
        let acquired = acquire_main_image_with_download_policy(
            &ctx,
            false,
            |_, _, _, _| panic!("local file must not fetch"),
            &|| panic!("local file must not resolve transport policy"),
        )
        .unwrap();
        assert!(acquired.main_owner.is_some());
    }

    #[test]
    fn direct_main_missing_image_admits_actual_lazy_policy_once() {
        let data = Sandbox::new();
        let ctx = data.context(NetworkPolicy::Allow);
        let calls = std::cell::Cell::new(0);
        let result = acquire_main_image_with_download_policy(
            &ctx,
            false,
            |_, _, _, _| panic!("actual sealed policy forbids fetch"),
            &|| {
                calls.set(calls.get() + 1);
                Ok(EffectiveNetworkPolicy::Sealed)
            },
        );
        assert!(matches!(result, Err(Error::NetworkSealed { .. })));
        assert_eq!(calls.get(), 1);
    }

    #[test]
    fn present_non_de_custom_precedes_data_and_never_fetches() {
        let data = Sandbox::new();
        let custom = Sandbox::new();
        let custom_path = custom.image("owned.bsp");
        std::fs::write(data.0.join("owned.bsp"), b"corrupt shadow").unwrap();
        let mut ctx = data
            .context(NetworkPolicy::Sealed)
            .with_ephemeris_file(Some("owned.bsp"));
        ctx.set_jpl_dir(Some(custom.0.clone()));
        let acquired = acquire_main_image_with(&ctx, |_, _, _, _| panic!("local fetch")).unwrap();
        assert!(!acquired.fallback_active);
        assert!(
            acquired
                .context
                .spk_registry()
                .kernel(custom_path.to_str().unwrap())
                .is_some()
        );
        assert_eq!(acquired.context.spk_registry().body_count(), 3);
        assert!(ctx.spk_registry().is_empty());
    }

    #[test]
    fn absolute_present_needs_no_data_directory() {
        let files = Sandbox::new();
        let path = files.image("neutral.bsp");
        let ctx = EphemerisContext::new()
            .with_mode(ProviderMode::NativeJpl)
            .with_network_policy(Some(NetworkPolicy::Sealed))
            .with_ephemeris_file(path.to_str());
        let acquired =
            acquire_main_image_with(&ctx, |_, _, _, _| panic!("absolute fetch")).unwrap();
        assert!(!acquired.fallback_active);
        assert_eq!(acquired.context.spk_registry().body_count(), 3);
    }

    #[test]
    fn absent_non_de_uses_local_tier_and_propagates_fallback_flag() {
        let files = Sandbox::new();
        files.image("de440.bsp");
        let ctx = files
            .context(NetworkPolicy::Sealed)
            .with_ephemeris_file(Some("missing.bsp"));
        let acquired =
            acquire_main_image_with(&ctx, |_, _, _, _| panic!("fallback fetch")).unwrap();
        assert!(acquired.fallback_active);
        assert_eq!(
            crate::positions::normalize_calc_flags(
                crate::constants::FLG_JPLEPH,
                acquired.fallback_active
            ) & 3,
            crate::constants::FLG_SWIEPH
        );
        let req = crate::positions::ServingRequest {
            ctx: &acquired.context,
            delta_t_model: crate::time::DeltaTModel::EspenakMeeus,
            iers_enabled: false,
            observer: None,
            tier_policy: crate::positions::TierPolicy::AttachOrder,
            skip_closed: false,
        };
        super::super::routing::with_local_jpl_entry_outputs_and_fallback(
            req,
            None,
            acquired.fallback_active,
            |engine, _| {
                assert!(engine.jpl_fallback_active);
                Ok(())
            },
        )
        .unwrap();
    }

    #[test]
    fn missing_sealed_is_exact_and_precedes_opener_and_fetch() {
        let files = Sandbox::new();
        for (name, selected) in [
            ("de441.bsp", "de441.bsp"),
            ("missing.bsp", "de440.bsp"),
            ("DE440.BSP", "DE440.BSP"),
        ] {
            let ctx = files
                .context(NetworkPolicy::Sealed)
                .with_ephemeris_file(Some(name));
            assert!(
                matches!(acquire_main_image_with(&ctx, |_,_,_,_| panic!("sealed fetch")),
                Err(Error::NetworkSealed { purpose }) if purpose == format!("JPL {selected} ephemeris download"))
            );
        }
    }

    #[test]
    fn permitted_acquisition_invokes_real_publication_then_owned_registration() {
        let files = Sandbox::new();
        let ctx = files.context(NetworkPolicy::Allow).with_tier(Tier::Base);
        let mut calls = 0;
        let acquired = acquire_main_image_with(&ctx, |url, destination, purpose, policy| {
            calls += 1;
            assert_eq!(url, format!("{JPL_ROOT}de440s.bsp"));
            assert_eq!(destination, files.0.join("de440s.bsp"));
            assert_eq!(purpose, "JPL de440s.bsp ephemeris download");
            assert_eq!(policy, EffectiveNetworkPolicy::Allow);
            std::fs::write(destination, test_image()).unwrap();
            Ok(())
        })
        .unwrap();
        assert_eq!(calls, 1);
        assert_eq!(acquired.context.spk_registry().body_count(), 3);
        assert!(ctx.spk_registry().is_empty());
    }

    #[test]
    fn populated_registry_never_resolves_bad_paths_or_fetches() {
        let files = Sandbox::new();
        let mut ctx = files
            .context(NetworkPolicy::Allow)
            .with_ephemeris_file(Some("bad.txt"));
        ctx.spk_registry_mut()
            .register_image("supplied", test_image())
            .unwrap();
        let before = ctx.clone();
        let acquired =
            acquire_main_image_with(&ctx, |_, _, _, _| panic!("supplied fetch")).unwrap();
        assert_eq!(acquired.context, before);
        assert!(!acquired.fallback_active);
    }

    #[test]
    fn corrupt_local_image_is_terminal_and_does_not_fetch() {
        let files = Sandbox::new();
        std::fs::write(files.0.join("de440.bsp"), b"not a kernel").unwrap();
        let ctx = files.context(NetworkPolicy::Allow);
        assert!(matches!(
            acquire_main_image_with(&ctx, |_, _, _, _| panic!("repair download")),
            Err(Error::SpkUnknownFormat { .. }) | Err(Error::SpkCorrupt { .. })
        ));
        assert!(ctx.spk_registry().is_empty());
    }

    #[test]
    fn allowed_fetch_failure_and_missing_publication_remain_errors() {
        let files = Sandbox::new();
        let ctx = files.context(NetworkPolicy::Allow);
        let failure = Error::HttpTransport {
            purpose: "probe".into(),
            detail: "injected".into(),
        };
        let got = acquire_main_image_with(&ctx, |_, _, _, _| Err(failure.clone()))
            .err()
            .unwrap();
        assert_eq!(got, failure);
        assert!(matches!(
            acquire_main_image_with(&ctx, |_, _, _, _| Ok(())),
            Err(Error::SpkFileNotFound { .. })
        ));
    }

    #[test]
    fn loader_literal_names_urls_and_format_refusals() {
        assert_eq!(
            loader_target("/missing/de440.bsp").unwrap(),
            (
                format!("{JPL_ROOT}/missing/de440.bsp"),
                "/missing/de440.bsp".into()
            )
        );
        assert_eq!(
            loader_target("nested/de440.bsp").unwrap().0,
            format!("{JPL_ROOT}nested/de440.bsp")
        );
        assert_eq!(
            loader_target("https://fixture.invalid/path/de440.bsp?q=1").unwrap(),
            (
                "https://fixture.invalid/path/de440.bsp?q=1".into(),
                "de440.bsp".into()
            )
        );
        assert_eq!(
            loader_target("jup365.bsp").unwrap().0,
            format!("{JUPITER_ROOT}jup365.bsp")
        );
        assert!(matches!(
            loader_target("bad.txt"),
            Err(Error::EphemerisLoaderValue { message }) if message.contains("bad.txt")
        ));
        let files = Sandbox::new();
        std::fs::write(files.0.join("DE440.BSP"), b"never opened").unwrap();
        assert!(matches!(
            acquire_main_image_with(
                &files
                    .context(NetworkPolicy::Sealed)
                    .with_ephemeris_file(Some("DE440.BSP")),
                |_, _, _, _| panic!("opener fetch")
            ),
            Err(Error::EphemerisLoaderValue { .. })
        ));
    }

    #[test]
    fn loader_opener_rejects_dot_only_basenames() {
        for name in [".bsp", "..bsp", "...bsp", "/present/.bsp", "/present/..bsp"] {
            assert!(matches!(
                loader_target(name),
                Err(Error::EphemerisLoaderValue { .. })
            ));
        }
        for name in ["ordinary.bsp", ".hidden.bsp", "/present/.hidden.bsp"] {
            assert_eq!(loader_target(name).unwrap().1, name);
        }
    }

    #[test]
    fn existing_noncanonical_spk_path_still_uses_its_actual_data() {
        let files = Sandbox::new();
        let name = "de٤٤٠.bsp";
        files.image(name);
        let result = acquire_main_image_with(
            &files
                .context(NetworkPolicy::Sealed)
                .with_ephemeris_file(Some(name)),
            |_, _, _, _| panic!("an existing local image must not download"),
        )
        .unwrap();
        assert!(!result.fallback_active);
        assert!(!result.context.spk_registry().is_empty());
    }

    #[test]
    fn configured_loader_root_failure_precedes_filename_policy_and_local_open() {
        let files = Sandbox::new();
        let blocked = files.0.join("not-a-directory");
        std::fs::write(&blocked, b"owned marker").unwrap();
        let actual = std::fs::create_dir_all(&blocked).unwrap_err();
        let expected = Error::EphemerisLoaderIo {
            path: blocked.to_string_lossy().into_owned(),
            kind: actual.kind(),
            errno: actual.raw_os_error(),
            message: actual.to_string(),
        };
        let local = files.image("local.bsp");
        for name in ["de441.bsp", "bad.txt", local.to_str().unwrap()] {
            let ctx = files
                .context(NetworkPolicy::Sealed)
                .with_jpl_data_dir(Some(blocked.clone()))
                .with_ephemeris_file(Some(name));
            assert_eq!(
                acquire_main_image_with(&ctx, |_, _, _, _| panic!("root error precedes fetch"))
                    .err()
                    .unwrap(),
                expected
            );
        }
    }

    #[test]
    fn configured_loader_root_is_created_before_sealed_missing_file() {
        let files = Sandbox::new();
        let root = files.0.join("new-loader-root");
        assert!(!root.exists());
        let ctx = files
            .context(NetworkPolicy::Sealed)
            .with_jpl_data_dir(Some(root.clone()));
        assert!(matches!(
            acquire_main_image_with(&ctx, |_, _, _, _| panic!("sealed fetch")),
            Err(Error::NetworkSealed { .. })
        ));
        assert!(root.is_dir());
    }
}
