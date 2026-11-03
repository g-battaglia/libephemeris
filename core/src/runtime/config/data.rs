//! Runtime-owned network, file acquisition, IERS and coverage operations.

use super::*;
use crate::context::{EffectiveNetworkPolicy, NetworkPolicy};
use crate::iers::{IersDataSnapshot, IersFileKind};
use std::collections::BTreeMap;
use std::path::Path;

/// Stored-body coverage metadata in TT Julian days.
///
/// An epoch-specific selection can identify a concrete image. An unlabeled
/// multi-file result is only the minimum/maximum envelope; it does not prove
/// that every intermediate epoch is served. `reviewed` identifies a manifest
/// digest match, not a scientific validation verdict for every calculation.
#[derive(Clone, Debug, PartialEq)]
pub struct BodyCoverageRecord {
    /// Public body identifier.
    pub body_id: i32,
    /// Calculation-source label.
    pub source: String,
    /// Precision class of the selected coefficients.
    pub precision_class: String,
    /// Inclusive stored start, TT Julian day.
    pub jd_start: f64,
    /// Inclusive stored end, TT Julian day.
    pub jd_end: f64,
    /// Selected file path, absent for a multi-file union without an epoch.
    pub data_file: Option<String>,
    /// Coefficient partition, if known from the path.
    pub group: Option<String>,
    /// True only when the immutable image matches its manifest digest.
    pub reviewed: bool,
}

impl BodyCoverageRecord {
    /// Whether `jd` lies in this closed TT-JD envelope.
    /// NaN compares false; no timescale conversion or actual provider query is
    /// performed. For a multi-file envelope, true does not rule out a gap.
    pub fn contains(&self, jd: f64) -> bool {
        self.jd_start <= jd && jd <= self.jd_end
    }
}
impl From<crate::leb::BodyCoverage> for BodyCoverageRecord {
    fn from(value: crate::leb::BodyCoverage) -> Self {
        Self {
            body_id: value.body_id(),
            source: value.source().to_owned(),
            precision_class: value.precision_class().as_str().to_owned(),
            jd_start: value.jd_start(),
            jd_end: value.jd_end(),
            data_file: value.data_file().map(str::to_owned),
            group: value.group().map(str::to_owned),
            reviewed: value.reviewed(),
        }
    }
}

/// Inspect retained coefficient images without opening or replacing resources.
///
/// `jd` is a coefficient TT epoch. A covered epoch selects its real image;
/// an uncovered epoch retains the stored bounds for diagnostics. Multiple
/// images without a selected serving interval report an unlabeled union.
/// Absent bodies and retired fictitious coefficient channels return `None`.
/// Closed images and malformed coefficient indices propagate their errors.
pub fn stored_body_coverage(
    files: &[crate::leb::LebFile],
    body: i32,
    jd: Option<f64>,
) -> Result<Option<BodyCoverageRecord>> {
    if (40..=58).contains(&body) {
        return Ok(None);
    }
    let mut bounds: Option<(f64, f64)> = None;
    for file in files {
        let Some((start, end)) = file.body_coverage(body)? else {
            continue;
        };
        bounds = Some(bounds.map_or((start, end), |(a, b)| (a.min(start), b.max(end))));
        if files.len() == 1 || jd.is_some_and(|jd| start <= jd && jd <= end) {
            return Ok(crate::leb::BodyCoverage::new(
                body,
                (start, end),
                Some(file.path()),
                Path::new(file.path())
                    .file_name()
                    .and_then(|name| name.to_str())
                    .and_then(crate::download::asset_by_name)
                    .is_some_and(|asset| asset.matches_image(file).unwrap_or(false)),
            )
            .map(Into::into));
        }
    }
    Ok(bounds
        .and_then(|bounds| crate::leb::BodyCoverage::new(body, bounds, None, false))
        .map(Into::into))
}

/// A tier acquisition report with the actual selected main-image path.
#[derive(Clone, Debug)]
pub struct AcquisitionSummary {
    /// Native per-body outcomes and counters.
    pub summary: crate::spk::EnsureSummary,
    /// Cached main image path, omitted for a newly acquired main image.
    pub ephemeris_path: Option<PathBuf>,
}

fn transport(timeout: f64) -> Result<crate::net::TransportConfig> {
    crate::net::TransportConfig::with_timeout_secs(timeout)
}

impl Runtime {
    /// Validate a scalar epoch against the active coefficient or main image.
    ///
    /// The epoch retains its caller-selected scale: this metadata check performs
    /// no UT/TT conversion. In coefficient mode, a body's stored interval is
    /// authoritative; an absent body is left to its analytical/provider route.
    /// With no body, the active file range is checked. Other modes acquire the
    /// actual main image before inspecting its range. Nonfinite epochs and
    /// outside-range epochs are errors; no extrapolation or network fallback
    /// can override an explicit coefficient-mode selection.
    pub fn validate_jd_range(&self, jd: f64, body: Option<i32>) -> Result<()> {
        if !jd.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd });
        }
        let context = self.effective_context()?;
        let bounds = if context.mode() == ProviderMode::SealedLeb {
            if let Some(body) = body {
                self.body_coverage(body, Some(jd))?
                    .map(|record| (record.jd_start, record.jd_end))
            } else {
                let session = self.shared_session()?;
                let mut bounds: Option<(f64, f64)> = None;
                for file in session.context().leb_files() {
                    let (start, end) = file.jd_range()?;
                    bounds = Some(bounds.map_or((start, end), |(a, b)| (a.min(start), b.max(end))));
                }
                bounds
            }
        } else {
            self.main_image()?;
            let (_, start, end, _) = self.current_file_data(0)?;
            Some((start, end))
        };
        if let Some((start, end)) = bounds {
            if (start != 0.0 || end != 0.0) && (jd < start || jd > end) {
                return Err(Error::EphemerisRange {
                    requested_jd: jd,
                    start_jd: start,
                    end_jd: end,
                    body,
                });
            }
        }
        Ok(())
    }

    /// Install requested ASSIST datasets atomically under the active policy.
    ///
    /// A supplied directory is used; otherwise the user data directory
    /// is `~/.libephemeris/assist`. Available valid files are retained unless
    /// `force` is true. Digests/format checks and transport are shared with the
    /// native installer. Policy, acquisition and filesystem failures propagate.
    pub fn download_assist_data(
        &self,
        target_dir: Option<&Path>,
        planets: bool,
        asteroids: bool,
        planets_de441: bool,
        force: bool,
    ) -> Result<PathBuf> {
        let target = match target_dir {
            Some(path) => path.to_path_buf(),
            None => self
                .sources()
                .home_directory()
                .map(|home| home.join(".libephemeris/assist"))
                .ok_or(Error::ConfigurationMissing {
                    missing_config: "assist_directory",
                    suggestion: "supply an explicit ASSIST target directory",
                })?,
        };
        let context = self.effective_context()?;
        crate::download::download_assist_data(
            &target,
            planets,
            asteroids,
            planets_de441,
            force,
            context.effective_network_policy(),
            &crate::net::TransportConfig::default(),
            None,
        )?;
        Ok(target)
    }

    /// Set this runtime's network override, not a process-global policy.
    ///
    /// Names are `auto`, `allow` and `sealed`, with surrounding whitespace and
    /// ASCII case ignored. `None` removes the explicit override and resumes
    /// captured environment/dotenv/TOML policy before the default. Effective
    /// `auto` seals LEB mode and permits the other provider modes.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidNetworkPolicy`] for an unknown name, preserving
    /// the selected policy. This call can still retire the cached session.
    /// Setting policy performs no socket work or cancellation of active calls.
    pub fn set_network_policy(&self, policy: Option<&str>) -> Result<()> {
        self.edit().context.set_network_policy_name(policy)
    }

    /// Configured auto/allow/sealed policy after explicit/environment resolution.
    pub fn configured_network_policy(&self) -> Result<NetworkPolicy> {
        Ok(self
            .effective_context()?
            .network_policy_override()
            .unwrap_or(NetworkPolicy::Auto))
    }

    /// Effective gate; auto seals LEB mode and allows the other provider modes.
    pub fn network_policy(&self) -> Result<EffectiveNetworkPolicy> {
        Ok(self.effective_context()?.effective_network_policy())
    }

    /// Merge dotenv values into the currently selected internal capture.
    ///
    /// `path` selects a file, or `None` uses the capture's normal dotenv
    /// selection. With `override_existing = false`, captured environment
    /// entries win; true allows replacement. See
    /// [`ConfigurationSources::load_dotenv`] for parsing and file semantics.
    /// No process environment variable is written.
    ///
    /// Returns false for an unreadable/missing file, true after reading one.
    /// On success, future resolution uses the updated capture and reader
    /// discovery is retired. Active snapshots survive. [`Self::sources`] still
    /// exposes the original construction-time capture, not this update.
    pub fn load_dotenv(&self, path: Option<&Path>, override_existing: bool) -> bool {
        let mut state = self.edit();
        let mut sources = state
            .reloaded_sources
            .clone()
            .unwrap_or_else(|| self.sources.clone());
        let loaded = sources.load_dotenv(path, override_existing);
        if loaded {
            state.reloaded_sources = Some(sources);
            state.center_images = None;
            state.invalidate_readers();
        }
        loaded
    }

    /// Library data root, independent of the explicit main-kernel search path.
    pub fn data_directory(&self) -> Result<PathBuf> {
        let context = self.effective_context()?;
        let path = context
            .jpl_data_dir()
            .cloned()
            .unwrap_or_else(|| PathBuf::from("."));
        std::path::absolute(&path).map_err(|error| Error::EphemerisLoaderIo {
            path: path.display().to_string(),
            kind: error.kind(),
            errno: error.raw_os_error(),
            message: error.to_string(),
        })
    }

    /// Explicit IERS cache override. None restores captured/default directory.
    pub fn set_iers_cache_dir(&self, path: Option<PathBuf>) {
        let mut state = self.edit();
        state.context.set_iers_cache_dir(path);
        state.context.clear_iers_cache();
    }

    /// Explicit IERS cache override, without resolving environment/default paths.
    pub fn iers_cache_dir_override(&self) -> Option<PathBuf> {
        self.lock().context.iers_cache_dir().cloned()
    }

    /// Override whether lazy IERS acquisition may download missing/forced files.
    pub fn set_iers_auto_download(&self, enabled: Option<bool>) {
        self.edit().context.set_iers_auto_download(enabled);
    }

    /// Resolved automatic IERS acquisition setting, independent of network gate.
    pub fn iers_auto_download(&self) -> Result<bool> {
        Ok(self.effective_context()?.effective_iers_auto_download())
    }

    /// Immutable parsed IERS rows, lazily loaded from the configured cache.
    /// The Arc remains alive after a reset, reload or directory change.
    pub fn iers_rows(&self) -> Result<Arc<IersDataSnapshot>> {
        self.effective_context()?.iers_source().row_snapshot()
    }

    /// Reparse IERS tables and publish a new immutable row generation.
    /// Force controls refresh only when automatic acquisition is enabled.
    pub fn load_iers_data(&self, force: bool) -> Result<bool> {
        let mut state = self.edit();
        let mut context = self.resolve(&state)?;
        let loaded = context.load_iers_data(force)?;
        state
            .context
            .attach_iers_source(context.iers_source().clone());
        Ok(loaded)
    }

    /// Release the default IERS generation without affecting active snapshots.
    pub fn clear_iers_cache(&self) {
        self.edit().context.clear_iers_cache();
    }

    /// Delete only the three selected IERS files and clear the default cache.
    /// Other contexts and in-flight row generations retain their owned data.
    pub fn delete_iers_cache_files(&self) -> Result<usize> {
        let mut state = self.edit();
        let mut context = self.resolve(&state)?;
        let deleted = context.delete_iers_cache_files()?;
        state
            .context
            .attach_iers_source(context.iers_source().clone());
        Ok(deleted)
    }

    /// Filesystem state, retained row counts and Julian-day ranges of IERS data.
    /// The native wall clock is used only for descriptive file ages in days.
    pub fn iers_cache_info(&self) -> Result<crate::iers::IersCacheInfo> {
        self.effective_context()?
            .iers_cache_info(crate::iers::iers_unix_seconds(std::time::SystemTime::now()))
    }

    /// Keep a fresh selected IERS file, or refresh primary then backup endpoints.
    /// A retained fresh file requires no network access. Download failures are
    /// returned as transport errors; failed data never replaces a valid file.
    pub fn download_iers_file(&self, kind: IersFileKind, force: bool) -> Result<PathBuf> {
        let context = self.effective_context()?;
        let directory = context.iers_source().prepare_directory()?;
        let path = directory.join(kind.filename());
        if !force
            && path.is_file()
            && crate::download::iers_file_age_days(&path)
                .is_some_and(crate::download::is_fresh_iers_age)
        {
            return Ok(path);
        }
        crate::download::refresh_iers_files(
            &directory,
            &[crate::download::IERS_FILES[kind.index()]],
            true,
            context.effective_network_policy(),
            &crate::net::TransportConfig::default(),
        )
        .map_err(|error| Error::HttpTransport {
            purpose: format!("download {}", kind.filename()),
            detail: error.to_string(),
        })?;
        Ok(path)
    }

    /// Resolve SPK input relative to cache root, data root, then current path.
    /// Missing paths remain explicit so opening reports their actual failure.
    pub fn resolve_spk_path(&self, path: &Path) -> Result<PathBuf> {
        if path.is_absolute() {
            return Ok(path.to_path_buf());
        }
        let context = self.effective_context()?;
        let data = self.data_directory()?;
        for directory in [
            context
                .spk_cache_dir()
                .cloned()
                .unwrap_or_else(|| data.join("spk")),
            data,
        ] {
            let candidate = directory.join(path);
            if candidate.exists() {
                return Ok(candidate);
            }
        }
        Ok(path.to_path_buf())
    }

    /// File-level TDB Julian-day envelope; absent/corrupt files return no span.
    pub fn spk_coverage(&self, path: &Path) -> Result<Option<(f64, f64)>> {
        let path = self.resolve_spk_path(path)?;
        let mut registry = crate::spk::SpkRegistry::new();
        Ok(registry
            .register_file(&path)
            .ok()
            .and_then(|key| registry.file_coverage_jd_tdb(&key)))
    }

    /// Register known SPK bodies discovered in a directory; returned labels
    /// report actual registration results. New immutable images are retained.
    pub fn discover_local_spks(&self, path: &Path) -> Result<BTreeMap<String, String>> {
        let mut state = self.edit();
        let entries = crate::spk::discover_local_spks(state.context.spk_registry_mut(), path)?;
        for record in state.context.spk_registry().registered_bodies() {
            state
                .spk_definitions
                .insert(record.body, (PathBuf::from(record.kernel), record.naif_id));
        }
        Ok(entries
            .into_iter()
            .map(|entry| (entry.body_name, entry.status))
            .collect())
    }

    /// Download a Horizons-generated SPK into an explicit or configured cache.
    /// Dates are Gregorian labels, center is a Horizons location, timeout is
    /// seconds. An existing valid file is kept unless overwrite is requested.
    #[allow(clippy::too_many_arguments)]
    pub fn download_spk(
        &self,
        body: &str,
        start: &str,
        end: &str,
        path: Option<&Path>,
        center: &str,
        overwrite: bool,
        timeout: f64,
    ) -> Result<PathBuf> {
        let context = self.effective_context()?;
        let cache = path
            .map(Path::to_path_buf)
            .or_else(|| context.spk_cache_dir().cloned())
            .unwrap_or(self.data_directory()?.join("spk"));
        let destination = crate::spk::spk_download_dest(&cache, body, start, end);
        crate::spk::download_horizons_spk_to(
            body,
            start,
            end,
            center,
            &destination,
            overwrite,
            context.effective_network_policy(),
            &transport(timeout)?,
        )?;
        Ok(destination)
    }

    /// Download and register a minor-body SPK, preserving valid cached coverage.
    /// Numeric target inference, coverage margins and transport stay in the
    /// native SPK controller; no registration is published after a failed load.
    #[allow(clippy::too_many_arguments)]
    pub fn download_and_register_spk(
        &self,
        body: &str,
        ipl: i32,
        start: &str,
        end: &str,
        naif: Option<i32>,
        path: Option<&Path>,
        center: &str,
        overwrite: bool,
        timeout: f64,
    ) -> Result<PathBuf> {
        let context = self.effective_context()?;
        let cache = path
            .map(Path::to_path_buf)
            .or_else(|| context.spk_cache_dir().cloned())
            .unwrap_or(self.data_directory()?.join("spk"));
        let mut state = self.edit();
        let result = crate::spk::download_and_register_spk_to(
            state.context.spk_registry_mut(),
            &cache,
            body,
            ipl,
            start,
            end,
            naif,
            center,
            overwrite,
            context.effective_network_policy(),
            &transport(timeout)?,
        )?;
        if let Some(record) = state.context.spk_registry().body_info(ipl) {
            state
                .spk_definitions
                .insert(record.body, (PathBuf::from(record.kernel), record.naif_id));
        }
        Ok(PathBuf::from(result))
    }

    /// Ensure the selected tier main image and registered minor-body SPKs.
    /// Files already present are used without replacement downloads. Transport
    /// respects the native policy; per-body fallbacks retain their real reason.
    pub fn ensure_all_ephemerides(&self, force: bool) -> Result<AcquisitionSummary> {
        let context = self.effective_context()?;
        let data = context
            .paths()
            .jpl_dir
            .clone()
            .unwrap_or(self.data_directory()?);
        let cache = context
            .spk_cache_dir()
            .cloned()
            .unwrap_or(self.data_directory()?.join("spk"));
        let tier = context.tier().as_name();
        let main = crate::context::tier_row(context.tier()).ephemeris_file;
        let mut state = self.edit();
        let report = crate::spk::ensure_all_ephemerides_to(
            tier,
            &data,
            &cache,
            state.context.spk_registry_mut(),
            force,
            context.effective_network_policy(),
            &crate::net::TransportConfig::default(),
        )?;
        for record in state.context.spk_registry().registered_bodies() {
            state
                .spk_definitions
                .insert(record.body, (PathBuf::from(record.kernel), record.naif_id));
        }
        let ephemeris_path = (report.ephemeris_status == "cached").then(|| data.join(main));
        Ok(AcquisitionSummary {
            summary: report,
            ephemeris_path,
        })
    }

    /// Install the tier's optional planet-center file and ensure its main/SPKs.
    /// This activates the tier; display/progress preferences never change data.
    pub fn download_for_tier(&self, tier: &str, force: bool) -> Result<AcquisitionSummary> {
        crate::download::check_tier(tier)?;
        self.set_precision_tier(tier)?;
        let context = self.effective_context()?;
        let data = self.data_directory()?;
        if let Some(entry) =
            crate::download::asset_by_name(&crate::download::planet_centers_name(tier)?)
        {
            if let Err(error) = crate::download::ensure_asset(
                entry,
                &crate::download::dest_path_for(&data, entry),
                force,
                context.effective_network_policy(),
                &crate::net::TransportConfig::default(),
                crate::download::validator_for(entry.name()),
                None::<&mut dyn FnMut(u64, Option<u64>)>,
            ) {
                log::warn!("planet-center provisioning for {tier}: {error}");
            }
        }
        self.ensure_all_ephemerides(force)
    }

    /// Ensure the reviewed tier-core image and optionally activate its path.
    /// The bundled base uses the same validated immutable image as calculations;
    /// other tiers use the existing pinned atomic installer. No substitute data.
    pub fn download_leb_for_tier(
        &self,
        tier: &str,
        force: bool,
        activate: bool,
    ) -> Result<PathBuf> {
        crate::download::check_tier(tier)?;
        let data = self.data_directory()?;
        let context = self.effective_context()?;
        #[cfg(feature = "bundled-data")]
        let installed = if tier == "base" {
            let file = crate::leb::bundled_base_core()?;
            let entry = crate::download::leb2_asset(tier, "core")?.ok_or_else(|| {
                Error::DataManifestMissing {
                    name: "base_core.leb2".into(),
                }
            })?;
            let path = crate::download::dest_path_for(&data, entry);
            if force
                || !path.is_file()
                || crate::download::file_sha256(&path).ok().as_deref() != Some(entry.sha256())
            {
                crate::download::install_bundled_bytes(
                    entry,
                    file.bytes()?,
                    &path,
                    crate::download::validator_for(entry.name()),
                )?;
            }
            path
        } else {
            crate::download::leb_core_path(&crate::download::ensure_leb_core_to(
                &data,
                tier,
                force,
                context.effective_network_policy(),
                &crate::net::TransportConfig::default(),
            )?)
            .to_path_buf()
        };
        #[cfg(not(feature = "bundled-data"))]
        let installed = crate::download::leb_core_path(&crate::download::ensure_leb_core_to(
            &data,
            tier,
            force,
            context.effective_network_policy(),
            &crate::net::TransportConfig::default(),
        )?)
        .to_path_buf();
        // Provisioning may add companions or atomically replace a selected
        // image even when activation leaves the path unchanged.
        self.lock().invalidate_readers();
        if activate {
            self.set_leb_file(Some(installed.clone()));
        }
        Ok(installed)
    }

    /// Coverage of one body in the active stored-reader hierarchy.
    /// An epoch selects its actual serving file; an undated multi-file query
    /// or an uncovered epoch returns diagnostic bounds without claiming a
    /// serving file. `None` means that no active image stores this body.
    pub fn body_coverage(&self, body: i32, jd: Option<f64>) -> Result<Option<BodyCoverageRecord>> {
        let session = self.shared_session()?;
        stored_body_coverage(session.context().leb_files(), body, jd)
    }

    /// File/body metadata of the actual configured stored images; no network.
    pub fn leb_inventory(&self) -> Result<crate::download::LebDirInventory> {
        let session = self.shared_session()?;
        Ok(crate::download::leb_inventory_for_context_with_verification(session.context(), true))
    }

    /// Manifest-derived cumulative tier requirements under the native data root.
    pub fn runtime_data_requirements(
        &self,
        tier: Option<&str>,
    ) -> Result<Vec<crate::leb::RuntimeDataRequirement>> {
        let selected = tier
            .map(str::to_owned)
            .unwrap_or(self.effective_context()?.tier().as_name().to_owned());
        let active = match selected.as_str() {
            "base" => crate::leb::LebTier::Base,
            "medium" => crate::leb::LebTier::Medium,
            "extended" => crate::leb::LebTier::Extended,
            _ => return Err(Error::InvalidTier { value: selected }),
        };
        let pins = crate::download::DATA_FILES
            .iter()
            .map(|entry| (entry.name(), entry.sha256()))
            .collect::<Vec<_>>();
        crate::leb::runtime_data_requirements(
            active,
            &crate::leb::LEB2_GROUPS,
            &pins,
            &self.data_directory()?.to_string_lossy(),
        )
    }
}
